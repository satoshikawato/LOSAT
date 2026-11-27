//! A Seq-table annotation as NCBI's object manager takes it: the checks of its columns when the
//! annotation is added (`CSeqTableInfo::x_Initialize`, `table_info`), the value of a column at a
//! row (`CSeqTable_column`, `CSeqTable_multi_data`, `CSeqTable_single_data`,
//! `CSeqTable_sparse_index`) and the index of its rows (`x_InitFeatTableKeys`, `index_table`).

use std::collections::HashMap;

use super::bitmagic::{self, BitVector, Fail};
use super::range::{Handle, HandleRangeMap, LiveHandles, Mode, Range};
use super::view::V;
use super::{Event, Indexer, RefElement, Refusal, Stop};
use crate::ncbi_exception::{NcbiException, Severity, ThrowSite};
// `ThrowUnassigned`, the getter of a member that is not set: the get check that
// `SERIAL_VERIFY_DATA_GET` turns off (the refusal of `super::index_entry`).
use crate::report::serial::object_frame::UNASSIGNED_MEMBER as UNASSIGNED;
use crate::utils::lower_bound::lower_bound_by;

/// A good feature table's `num-rows` above which NCBI's index depends on the machine's memory
/// (`x_InitFeatTable` makes an annotation object of about 135 bytes per row: 10,000,000 rows take
/// 1.37 GB; line RG D3-1 `om/SPEC.md` §8 row 26, §10 Q2).
const ROWS_LIMIT: i32 = 10_000_000;

// The places of NCBI's exceptions (line RG D3-1 `om/SPEC.md` §7.1, NCBI BLAST+ 2.17.0 runs).
// NCBI reference (598d8ae6): c++/src/objmgr/seq_table_info.cpp:351-357 (`SetColumn`)
const SET_COLUMN: ThrowSite = ThrowSite {
    file: "src/objmgr/seq_table_info.cpp",
    line: 355,
    kind: "CAnnotException::eBadLocation",
    function: Some("ncbi::objects::CSeqTableLocColumns::SetColumn"),
};

/// The throws of `CSeqTableLocColumns::ParseDefaults` (seq_table_info.cpp:459-521, cited at
/// `LocColumns::parse_defaults`) at `line`.
const fn parse_defaults_site(line: u32) -> ThrowSite {
    ThrowSite {
        file: "src/objmgr/seq_table_info.cpp",
        line,
        kind: "CAnnotException::eBadLocation",
        function: Some("ncbi::objects::CSeqTableLocColumns::ParseDefaults"),
    }
}

/// The throws of `CSeqTableInfo::x_Initialize` (seq_table_info.cpp:784-947) at `line`.
const fn initialize_site(line: u32) -> ThrowSite {
    ThrowSite {
        file: "src/objmgr/seq_table_info.cpp",
        line,
        kind: "CAnnotException::eOtherError",
        function: Some("ncbi::objects::CSeqTableInfo::x_Initialize"),
    }
}

// The generated `ThrowInvalidSelection` of the two CHOICEs (datatool's code in NCBI's build tree;
// NCBI BLAST+ 2.17.0 reports these places: probes `t30`, `t37`, `t17`, `t74`).
const MULTI_SELECTION: ThrowSite = ThrowSite {
    file: "src/objects/seq/../seqtable/SeqTable_multi_data_.cpp",
    line: 207,
    kind: "CInvalidChoiceSelection::eFail",
    function: None,
};
const SINGLE_SELECTION: ThrowSite = ThrowSite {
    file: "src/objects/seq/../seqtable/SeqTable_single_data_.cpp",
    line: 135,
    kind: "CInvalidChoiceSelection::eFail",
    function: None,
};

// NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_multi_data.cpp:435-444, 557-561
const MULTI_TOO_BIG: ThrowSite = ThrowSite {
    file: "src/objects/seq/../seqtable/SeqTable_multi_data.cpp",
    line: 441,
    kind: "CSeqTableException::eIncompatibleValueType",
    function: None,
};
const MULTI_CANNOT_CONVERT: ThrowSite = ThrowSite {
    line: 560,
    ..MULTI_TOO_BIG
};

// NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_single_data.cpp:54-70
const SINGLE_TOO_BIG: ThrowSite = ThrowSite {
    file: "src/objects/seq/../seqtable/SeqTable_single_data.cpp",
    line: 58,
    kind: "CSeqTableException::eIncompatibleValueType",
    function: Some("ncbi::objects::CSeqTable_single_data::ThrowOverflowError"),
};
const SINGLE_CANNOT_CONVERT: ThrowSite = ThrowSite {
    line: 66,
    function: None,
    ..SINGLE_TOO_BIG
};

/// The value type a getter asks for: `bool`, `int` (Int4) or `Int8`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Want {
    Bool,
    Int4,
    Int8,
}

impl Want {
    /// The name in `CSeqTable_multi_data::TryGet<name>()`.
    fn multi_name(self) -> &'static str {
        match self {
            Want::Bool => "Bool",
            Want::Int4 => "Int4",
            Want::Int8 => "Int8",
        }
    }

    /// The name in `CSeqTable_single_data::GetValue(<name>&)`.
    fn single_name(self) -> &'static str {
        match self {
            Want::Bool => "bool",
            Want::Int4 => "int",
            Want::Int8 => "Int8",
        }
    }

    /// `DstInt` of `value`: the value in the type (a bool is true when `value` is not 0).
    fn cast(self, value: i64) -> i64 {
        match self {
            Want::Bool => i64::from(value != 0),
            Want::Int4 => i64::from(value as i32),
            Want::Int8 => value,
        }
    }
}

/// The sums that NCBI's caches keep for a column's data, by the data's address: the prefix sums
/// of an `int-delta`'s elements (`CInt4DeltaSumCache`, `CInt8DeltaSumCache`), of a sparse index's
/// `indexes-delta` (`CIndexDeltaSumCache`) and the set bits before each byte of its `bit-set`
/// (`x_GetBitSetCache`). The values are those of NCBI's sums; the caches only keep a large table's
/// rows from being summed again for each row.
#[derive(Default)]
pub(super) struct Caches {
    int_delta: HashMap<usize, Vec<i64>>,
    index_delta: HashMap<usize, Vec<u64>>,
    bit_set: HashMap<usize, Vec<u64>>,
    /// `GetSize()` of each multi-data value asked (`Index::multi_size`).
    size: HashMap<usize, u64>,
    /// `GetIntSize()` of each `int-delta`'s data asked (`Index::int_size`).
    int_size: HashMap<usize, u64>,
    /// The vector of each `BVector-data` decoded (`Indexer::bit_vector`), or why it was not.
    bvector: HashMap<usize, Result<BitVector, Fail>>,
}

fn address(value: V<'_>) -> usize {
    value.value as *const _ as usize
}

/// A column of the table (`SeqTable-column`).
#[derive(Clone, Copy, Debug)]
pub(super) struct Column<'a>(V<'a>);

impl<'a> Column<'a> {
    fn header(self) -> V<'a> {
        self.0.field("header")
    }
    fn field_id(self) -> Option<i32> {
        self.header().get("field-id").map(|id| id.int() as i32)
    }
    fn field_name(self) -> Option<&'a [u8]> {
        self.header().get("field-name").map(V::str)
    }
    fn data(self) -> Option<V<'a>> {
        self.0.get("data")
    }
    fn sparse(self) -> Option<V<'a>> {
        self.0.get("sparse")
    }
    fn default(self) -> Option<V<'a>> {
        self.0.get("default")
    }
    fn sparse_other(self) -> Option<V<'a>> {
        self.0.get("sparse-other")
    }
    /// `CSeqTableColumnInfo::IsSingular`: a default and neither data nor a sparse index.
    fn is_singular(self) -> bool {
        self.default().is_some() && self.data().is_none() && self.sparse().is_none()
    }
}

/// `CSeqTableLocColumns`: the columns of the location (`loc`, base 0) or of the product
/// (`product`, base 10).
#[derive(Debug)]
pub(super) struct LocColumns<'a> {
    name: &'static str,
    base: i32,
    loc: Option<Column<'a>>,
    id: Option<Column<'a>>,
    gi: Option<Column<'a>>,
    from: Option<Column<'a>>,
    to: Option<Column<'a>>,
    strand: Option<Column<'a>>,
    extra: Vec<Column<'a>>,
    is_set: bool,
    real_loc: bool,
    simple: bool,
    probably_simple: bool,
    simple_point: bool,
    simple_interval: bool,
    default_id: Option<Handle>,
}

impl<'a> LocColumns<'a> {
    fn new(name: &'static str, base: i32) -> Self {
        LocColumns {
            name,
            base,
            loc: None,
            id: None,
            gi: None,
            from: None,
            to: None,
            strand: None,
            extra: Vec::new(),
            is_set: false,
            real_loc: false,
            simple: false,
            probably_simple: false,
            simple_point: false,
            simple_interval: false,
            default_id: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_table_info.cpp:351-372 (`SetColumn`, cited at
    // `SET_COLUMN`; `AddExtraColumn` sets `m_Is_set` too)
    fn set(&self, slot: Option<Column<'a>>, column: Column<'a>) -> Result<Column<'a>, Stop> {
        if slot.is_some() {
            return Err(Stop::exception(NcbiException::new(
                SET_COLUMN,
                format!("Duplicate {} column", self.name),
            )));
        }
        Ok(column)
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_table_info.cpp:371-454
    // `GetField_id()` minus `m_BaseValue` is `int` arithmetic: below 10 more than `INT_MIN` the
    // product's overflows (undefined behaviour), refused.
    fn add_column(&mut self, column: Column<'a>) -> Result<bool, Stop> {
        if let Some(id) = column.field_id() {
            let Some(field) = id.checked_sub(self.base) else {
                return Err(Stop::Refuse(Refusal::FieldIdOverflow));
            };
            if !(0..10).contains(&field) {
                return Ok(false);
            }
            match field {
                0 => self.loc = Some(self.set(self.loc, column)?),
                1 => self.id = Some(self.set(self.id, column)?),
                2 => self.gi = Some(self.set(self.gi, column)?),
                3 => self.from = Some(self.set(self.from, column)?),
                4 => self.to = Some(self.set(self.to, column)?),
                5 => self.strand = Some(self.set(self.strand, column)?),
                6 | 7 => self.extra.push(column),
                _ => return self.add_column_by_name(column),
            }
            self.is_set = true;
            return Ok(true);
        }
        self.add_column_by_name(column)
    }

    fn add_column_by_name(&mut self, column: Column<'a>) -> Result<bool, Stop> {
        let Some(field) = column.field_name() else {
            return Ok(false);
        };
        let name = self.name.as_bytes();
        if field == name {
            self.loc = Some(self.set(self.loc, column)?);
        } else if field.starts_with(name) && field.get(name.len()) == Some(&b'.') {
            let extra = &field[name.len() + 1..];
            if extra == b"id" || extra.ends_with(b".id") {
                self.id = Some(self.set(self.id, column)?);
            } else if extra == b"gi" || extra.ends_with(b".gi") {
                self.gi = Some(self.set(self.gi, column)?);
            } else if extra == b"pnt.point" || extra == b"int.from" {
                self.from = Some(self.set(self.from, column)?);
            } else if extra == b"int.to" {
                self.to = Some(self.set(self.to, column)?);
            } else if extra == b"strand" || extra.ends_with(b".strand") {
                self.strand = Some(self.set(self.strand, column)?);
            } else if extra == b"int.fuzz-from.lim"
                || extra == b"pnt.fuzz.lim"
                || extra == b"int.fuzz-to.lim"
            {
                self.extra.push(column);
            } else {
                return Ok(false);
            }
        } else {
            return Ok(false);
        }
        self.is_set = true;
        Ok(true)
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_table_info.cpp:457-521
    fn parse_defaults(&mut self, handles: &mut LiveHandles<'a>) -> Result<(), Stop> {
        if !self.is_set {
            return Ok(());
        }
        let fail = |line, message: String| {
            Err(Stop::exception(NcbiException::new(
                parse_defaults_site(line),
                message,
            )))
        };
        if self.loc.is_some() {
            self.real_loc = true;
            if self.id.is_some()
                || self.gi.is_some()
                || self.from.is_some()
                || self.to.is_some()
                || self.strand.is_some()
                || !self.extra.is_empty()
            {
                return fail(466, format!("Conflicting {} columns", self.name));
            }
            return Ok(());
        }
        if self.id.is_none() && self.gi.is_none() {
            return fail(473, format!("No {}.id column", self.name));
        }
        if self.id.is_some() && self.gi.is_some() {
            return fail(477, format!("Conflicting {} columns", self.name));
        }
        if let Some(default) = self.id.and_then(Column::default) {
            // The default handle lives as long as the table (the run).
            self.default_id = Some(Handle::of_run(single_get_id(default)?.seq_id(), handles));
        }
        if let Some(default) = self.gi.and_then(Column::default) {
            self.default_id = Some(Handle::gi(single_get(default, Want::Int8)?));
        }
        if self.to.is_some() {
            if self.from.is_none() {
                return fail(497, format!("column {0}.to without {0}.from", self.name));
            }
            self.simple_interval = true;
        } else if self.from.is_some() {
            self.simple_point = true;
        } else if self.strand.is_some() || !self.extra.is_empty() {
            return fail(510, format!("extra columns in whole {}", self.name));
        }
        if self.extra.is_empty() {
            self.simple = true;
        } else {
            self.probably_simple = true;
        }
        Ok(())
    }
}

/// `CSeqTableInfo`: what `x_Initialize` keeps of a table.
#[derive(Debug)]
pub(super) struct TableInfo<'a> {
    pub(super) is_feat: bool,
    num_rows: i32,
    table_loc: Option<V<'a>>,
    sorted: bool,
    loc: LocColumns<'a>,
    product: LocColumns<'a>,
    partial: Option<Column<'a>>,
    disabled: Option<Column<'a>>,
}

// NCBI reference (598d8ae6): c++/src/objmgr/seq_table_info.cpp:750-764
// `e_MaxChoice` is 23 (SeqFeatData_.hpp:174), `eSubtype_max` 106 (SeqFeatData.hpp:237).
fn is_good_feat_table(table: V<'_>) -> bool {
    let feat_type = table.field("feat-type").int();
    if !(1..23).contains(&feat_type) {
        return false;
    }
    table
        .get("feat-subtype")
        .is_none_or(|subtype| (1..106).contains(&subtype.int()))
}

impl<'p> Indexer<'p> {
    // NCBI reference (598d8ae6): c++/src/objmgr/seq_table_info.cpp:784-947
    // The maps of the columns by id and by name and the setters of the feature fields make no
    // message (the setters' constructors only keep the name: `CSeqTableSetExt` and the others
    // cut it with `CTempString::substr`, `CSeqTableSetAnyFeatField` throws only the
    // `CAnnotException`s caught here, seq_table_setters.cpp:180-188, 344-430).
    /// The making of a `CSeqTableInfo` from `table` (`x_InitFeatTable`,
    /// seq_annot_info.cpp:474-500): the checks of the columns, with their messages, then LOSAT's
    /// bound of the rows that NCBI makes an object for.
    pub(super) fn table_info(&mut self, table: V<'p>) -> Result<TableInfo<'p>, Stop> {
        let is_feat = is_good_feat_table(table);
        let mut info = TableInfo {
            is_feat,
            num_rows: table.field("num-rows").int() as i32,
            table_loc: None,
            sorted: false,
            loc: LocColumns::new("loc", 0),
            product: LocColumns::new("product", 10),
            partial: None,
            disabled: None,
        };
        let mut sorted_max: u32 = 0;
        for column in table.field("columns").items() {
            let column = Column(column.deref()?);
            if column.field_name() == Some(b"Seq-table location") {
                if info.table_loc.is_some() {
                    self.post(
                        Severity::Error,
                        "",
                        b"Duplicate 'Seq-table location' column",
                    );
                }
                match self.column_seq_loc(column, 0) {
                    Ok(loc) => info.table_loc = loc,
                    Err(Stop::Exception(e)) => {
                        self.ncbi_catch(&e, "Bad 'Seq-table location' column")
                    }
                    Err(Stop::Std(what)) => self.std_catch(what, "Bad 'Seq-table location' column"),
                    Err(other) => return Err(other),
                }
            }
            if !is_feat {
                continue;
            }
            if info.loc.add_column(column)? || info.product.add_column(column)? {
                continue;
            }
            let mut setter = false;
            if let Some(id) = column.field_id() {
                match id {
                    22 => {
                        if info.partial.is_some() {
                            return Err(Stop::exception(NcbiException::new(
                                initialize_site(835),
                                "Duplicate partial column",
                            )));
                        }
                        info.partial = Some(column);
                        continue;
                    }
                    23 | 30 | 31 | 40 => setter = true,
                    25..=27 => {
                        if column.field_name().is_none() {
                            return Err(Stop::exception(NcbiException::new(
                                UNASSIGNED,
                                "CSeqTable_column_info::GetField_name(): Attempt to get unassigned member NCBI-SeqTable::SeqTable-column-info.field-name",
                            )));
                        }
                        setter = true;
                    }
                    _ => {
                        if column.field_name().is_none() {
                            self.post(
                                Severity::Error,
                                "1213.8",
                                format!("SeqTable-column-info.field-id = {id}").as_bytes(),
                            );
                            continue;
                        }
                    }
                }
            } else if column.field_name().is_none() {
                self.post(
                    Severity::Error,
                    "1213.9",
                    b"SeqTable-column-info: neither field-id nor field-name is set",
                );
                continue;
            }
            if setter {
                continue;
            }
            let Some(name) = column.field_name() else {
                continue;
            };
            match name {
                b"" => {
                    self.post(
                        Severity::Error,
                        "1213.10",
                        b"SeqTable-column-info.field-name is empty",
                    );
                }
                [b'E' | b'D' | b'Q', ..] => {}
                b"partial" => {
                    if info.partial.is_some() {
                        return Err(Stop::exception(NcbiException::new(
                            initialize_site(891),
                            "Duplicate partial column",
                        )));
                    }
                    info.partial = Some(column);
                }
                b"disabled" => {
                    if info.disabled.is_some() {
                        return Err(Stop::exception(NcbiException::new(
                            initialize_site(899),
                            "Duplicate disabled column ",
                        )));
                    }
                    info.disabled = Some(column);
                }
                b"Sorted, max length" => {
                    if sorted_max != 0 {
                        self.post(
                            Severity::Error,
                            "",
                            b"Duplicate 'Sorted, max length' column",
                        );
                    }
                    match self.column_try_get(column, 0, Want::Int4) {
                        Ok(Some(value)) => sorted_max = value as i32 as u32,
                        Ok(None) => {}
                        Err(Stop::Exception(e)) => {
                            self.ncbi_catch(&e, "Bad 'Sorted, max length' column")
                        }
                        Err(Stop::Std(what)) => {
                            self.std_catch(what, "Bad 'Sorted, max length' column")
                        }
                        Err(other) => return Err(other),
                    }
                }
                _ => {}
            }
        }
        if is_feat {
            info.loc.parse_defaults(&mut self.handles)?;
            info.product.parse_defaults(&mut self.handles)?;
        }
        info.sorted = info.is_sorted(sorted_max);
        // NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:486-499 (`x_InitFeatTable`)
        if is_feat && !info.sorted && info.num_rows > ROWS_LIMIT {
            return Err(Stop::Refuse(Refusal::TableRows));
        }
        Ok(info)
    }

    /// `NCBI_CATCH` with `title` of `e` (`NcbiException::report_default`): one error post per line.
    fn ncbi_catch(&mut self, e: &NcbiException, title: &str) {
        let lines = e.report_default(self.program, title);
        for line in lines.split_inclusive(|&byte| byte == b'\n') {
            self.posts.post(crate::report::diag_log::ERROR, line);
        }
    }

    // NCBI reference (598d8ae6): c++/include/corelib/ncbiexpt.hpp:560-565, 580-584
    /// `NCBI_CATCH` with `title` of a `std::exception` whose `what()` is `what`: `Error: [<program>]
    /// [<title>] Exception: <what>` (line RG D3-1 `bm/` cases m12, m13: BitMagic's
    /// `BM::Invalid serialization format`).
    fn std_catch(&mut self, what: &str, title: &str) {
        let line = crate::ncbi_exception::error_post_line(
            None,
            Some(self.program),
            format!("[{title}] Exception: {what}").as_bytes(),
        );
        self.posts.post(crate::report::diag_log::ERROR, line);
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:1028-1170 (`x_InitFeatTableKeys`)
    /// `x_InitFeatTableKeys` of the table `info`.
    pub(super) fn index_table(&mut self, info: &TableInfo<'_>) -> Result<(), Stop> {
        if !info.is_feat {
            if let Some(loc) = info.table_loc {
                let map = HandleRangeMap::of_location(loc, Mode::None, &mut self.handles)?;
                self.map_location(&map, Event::TableLocation, loc)?;
                self.handles.end_object();
            }
            return Ok(());
        }
        if info.sorted {
            // `x_IsSorted` holds: the table location is an `int`.
            let (_, interval) = info.table_loc.expect("sorted").choice();
            let (id, range, _) = super::range::interval_parts(interval);
            // `key.m_Handle` holds the handle (mapped, unless the range is empty) until the table
            // has been indexed.
            let handle = Handle::of(id, &mut self.handles);
            let mapped = self.x_map(range, Some(handle.clone()));
            self.handles.set_key(Some(&handle));
            self.handles.end_object();
            self.handles.end_table();
            return mapped;
        }
        for row in 0..info.num_rows.max(0) as u64 {
            // NCBI reference (598d8ae6): c++/include/objmgr/impl/seq_table_info.hpp:92-96,241-246
            if let Some(disabled) = info.disabled {
                if self.column_try_get(disabled, row, Want::Bool)? == Some(1) {
                    continue;
                }
            }
            if let Some(partial) = info.partial {
                self.column_try_get(partial, row, Want::Bool)?;
            }
            for cols in [&info.loc, &info.product] {
                if !cols.is_set {
                    continue;
                }
                if cols.real_loc {
                    let column = cols.loc.expect("a real location");
                    // NCBI dereferences `GetLoc` of a row without a location: a null
                    // `CConstRef` (`CCoreException::eNullPtr`, a report of memory addresses).
                    let Some(loc) = self.column_seq_loc(column, row)? else {
                        return Err(self.null_pointer());
                    };
                    let map = HandleRangeMap::of_location(loc, Mode::None, &mut self.handles)?;
                    self.map_location(&map, Event::Row, loc)?;
                    // `key.m_Handle` takes each handle of a non-empty range (mapped).
                    if map.map.iter().any(|(_, hr)| !hr.overlapping().is_empty()) {
                        self.handles.set_key(None);
                    }
                } else {
                    let (id, range) = self.table_key(cols, row)?;
                    self.x_map(range, id.clone())?;
                    // `SetTableKeyAndIndex` puts the row's handle in `key.m_Handle`.
                    self.handles.set_key(id.as_ref());
                }
                // NCBI's row map is made for each location (c++/src/objmgr/seq_annot_info.cpp:1130).
                self.handles.end_object();
            }
        }
        self.handles.end_table();
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_table_info.cpp:701-746 (`SetTableKeyAndIndex`)
    // NCBI reference (598d8ae6): c++/src/objmgr/seq_table_info.cpp:538-603
    // The `_ASSERT`s are the debug build's: the release build reads a sparse index as any column.
    // The strand and the extra columns change no output, but their values are read (and can
    // throw).
    /// The key of a simple row (`SetTableKeyAndIndex`): its Seq-id handle (`None`: null) and range.
    fn table_key(
        &mut self,
        cols: &LocColumns<'_>,
        row: u64,
    ) -> Result<(Option<Handle>, Range), Stop> {
        let mut handle = None;
        if let Some(id) = cols.id {
            if id.data().is_some() {
                handle = self
                    .column_seq_id(id, row)?
                    .map(|id| Handle::of(id.seq_id(), &mut self.handles));
            }
        } else if let Some(gi) = cols.gi.filter(|gi| gi.data().is_some()) {
            handle = self.column_try_get(gi, row, Want::Int8)?.map(Handle::gi);
        }
        let handle = handle.or_else(|| cols.default_id.clone());
        let range = match cols.from {
            Some(from) => match self.column_try_get(from, row, Want::Int4)? {
                Some(from) => {
                    let mut to = from;
                    if let Some(column) = cols.to {
                        if let Some(value) = self.column_try_get(column, row, Want::Int4)? {
                            to = value;
                        }
                    }
                    Range::new(from as i32 as u32, to as i32 as u32)
                }
                None => Range::WHOLE,
            },
            None => Range::WHOLE,
        };
        if let Some(strand) = cols.strand {
            self.column_try_get(strand, row, Want::Int4)?;
        }
        if !cols.simple && cols.probably_simple {
            for &extra in &cols.extra {
                if self.column_is_set(extra, row)? {
                    break;
                }
            }
        }
        Ok((handle, range))
    }
}

impl TableInfo<'_> {
    // NCBI reference (598d8ae6): c++/src/objmgr/seq_table_info.cpp:972-1005
    // NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_interval.hpp:125-128
    fn is_sorted(&self, sorted_max: u32) -> bool {
        if self.product.is_set {
            return false;
        }
        let Some(id) = self
            .loc
            .id
            .filter(|_| self.loc.is_set && !self.loc.real_loc)
        else {
            return false;
        };
        if !id.is_singular() {
            return false;
        }
        if !self.loc.simple || !(self.loc.simple_point || self.loc.simple_interval) {
            return false;
        }
        let Some(("int", interval)) = self.table_loc.map(V::choice) else {
            return false;
        };
        let from = interval.field("from").int() as u32;
        let to = interval.field("to").int() as u32;
        let length = if to >= from {
            (to - from).wrapping_add(1)
        } else {
            0
        };
        sorted_max != 0 && sorted_max <= length / 16
    }
}

// The values of a column at a row.
impl Indexer<'_> {
    // NCBI reference (598d8ae6): c++/include/objects/seqtable/SeqTable_column.hpp:75-96
    /// `CSeqTable_column::TryGetValue` of the type `want` (a bool as 0 or 1).
    pub(super) fn column_try_get(
        &mut self,
        column: Column<'_>,
        row: u64,
        want: Want,
    ) -> Result<Option<i64>, Stop> {
        let mut index = row;
        if let Some(sparse) = column.sparse() {
            match self.sparse_index_at(sparse, row)? {
                Some(at) => index = at,
                None => {
                    return match column.sparse_other() {
                        Some(other) => single_get(other, want).map(Some),
                        None => Ok(None),
                    }
                }
            }
        }
        if let Some(data) = column.data() {
            if let Some(value) = self.multi_try_get(data, index, want)? {
                return Ok(Some(value));
            }
        }
        match column.default() {
            Some(default) => single_get(default, want).map(Some),
            None => Ok(None),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_column.cpp:62-75
    // NCBI reference (598d8ae6): c++/include/objects/seqtable/SeqTable_multi_data.hpp:71-73
    fn column_is_set(&mut self, column: Column<'_>, row: u64) -> Result<bool, Stop> {
        let mut index = row;
        if let Some(sparse) = column.sparse() {
            match self.sparse_index_at(sparse, row)? {
                Some(at) => index = at,
                None => return Ok(column.sparse_other().is_some()),
            }
        }
        if let Some(data) = column.data() {
            if index < self.multi_size(data) {
                return Ok(true);
            }
        }
        Ok(column.default().is_some())
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_column.cpp:126-174
    /// `CSeqTable_column::GetSeq_loc` of `row`.
    fn column_seq_loc<'a>(&mut self, column: Column<'a>, row: u64) -> Result<Option<V<'a>>, Stop> {
        self.column_object(column, row, "loc", "Loc")
    }

    /// `CSeqTable_column::GetSeq_id` of `row`.
    fn column_seq_id<'a>(&mut self, column: Column<'a>, row: u64) -> Result<Option<V<'a>>, Stop> {
        self.column_object(column, row, "id", "Id")
    }

    fn column_object<'a>(
        &mut self,
        column: Column<'a>,
        row: u64,
        variant: &str,
        getter: &str,
    ) -> Result<Option<V<'a>>, Stop> {
        let mut index = row;
        if let Some(sparse) = column.sparse() {
            match self.sparse_index_at(sparse, row)? {
                Some(at) => index = at,
                None => {
                    return match column.sparse_other() {
                        Some(other) => single_object(other, variant, getter).map(Some),
                        None => Ok(None),
                    }
                }
            }
        }
        if let Some(data) = column.data() {
            let (current, list) = data.choice();
            if current != variant {
                return Err(Stop::exception(NcbiException::new(
                    MULTI_SELECTION,
                    format!(
                        "CSeqTable_multi_data::Get{getter}(): Invalid choice selection: NCBI-SeqTable::SeqTable-multi-data.{current}"
                    ),
                )));
            }
            // An element stated as `NULL` is the row's empty reference, as a row without a value.
            if let Some(element) = usize::try_from(index)
                .ok()
                .and_then(|index| list.item(index))
            {
                return Ok((!element.is_null()).then_some(element));
            }
        }
        match column.default() {
            Some(default) => single_object(default, variant, getter).map(Some),
            None => Ok(None),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_multi_data.cpp:477-562
    // NCBI reference (598d8ae6): c++/src/objects/seqtable/Scaled_int_multi_data.cpp:77-84
    // NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_multi_data.cpp:474-484
    // A bit-bvector is BitMagic's serialized vector, which `GetBitVector()` decodes at its first
    // value (`Indexer::bit_vector`); its row is cast to BitMagic's 32-bit `size_type` first. An
    // `int-scaled` value that leaves `Int8` is a signed overflow: refused.
    /// `x_TryGetInt8` of the multi-data `data`, one level: its value, or the
    /// multi-data inside it to read first (`Get`).
    fn x_try_get_int8<'v>(
        &mut self,
        data: V<'v>,
        row: u64,
        type_name: &str,
    ) -> Result<Get<'v>, Stop> {
        let (variant, value) = data.choice();
        let element = |value: V<'_>| usize::try_from(row).ok().and_then(|row| value.int_at(row));
        Ok(Get::Value(match variant {
            "bit" => {
                let bytes = value.octets();
                let Some(&byte) = usize::try_from(row / 8).ok().and_then(|i| bytes.get(i)) else {
                    return Ok(Get::Value(None));
                };
                let shifted = u32::from(byte) << (row % 8);
                Some(i64::from(shifted & 0x80 != 0))
            }
            "bit-bvector" => {
                let vector = self.bit_vector(value)?;
                let row = row as u32;
                if row < vector.size() {
                    Some(i64::from(vector.get_bit(row).map_err(query_stop)?))
                } else {
                    None
                }
            }
            "int1" => usize::try_from(row)
                .ok()
                .and_then(|row| value.octets().get(row))
                .map(|&byte| i64::from(byte as i8)),
            "int2" | "int" | "int8" => element(value),
            "int-delta" => {
                if row >= self.multi_size(value) {
                    return Ok(Get::Value(None));
                }
                return Ok(Get::Delta(value));
            }
            "int-scaled" => {
                return Ok(Get::Scaled {
                    data: value.field("data"),
                    mul: value.field("mul").int() as i32,
                    add: value.field("add").int() as i32,
                })
            }
            _ => {
                return Err(Stop::exception(NcbiException::new(
                    MULTI_CANNOT_CONVERT,
                    format!(
                        "CSeqTable_multi_data::TryGet{type_name}(): value cannot be converted to requested type"
                    ),
                )))
            }
        }))
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_multi_data.cpp:435-516
    /// `CSeqTable_multi_data::TryGetValue` of the type `want`.
    fn multi_try_get(&mut self, data: V<'_>, row: u64, want: Want) -> Result<Option<i64>, Stop> {
        // The multi-data values waiting for one inside them, from a list on the heap (to any depth
        // the reader read), and the one to read next.
        let mut pending = vec![Pending::Cast(want)];
        let mut get = Some((data, row, want.multi_name()));
        let mut result = Ok(None);
        loop {
            if let Some((data, row, type_name)) = get.take() {
                result = match self.x_try_get_int8(data, row, type_name) {
                    Err(stop) => Err(stop),
                    Ok(Get::Value(value)) => Ok(value),
                    Ok(Get::Scaled { data, mul, add }) => {
                        pending.push(Pending::Scale { mul, add });
                        get = Some((data, row, "Int8"));
                        continue;
                    }
                    Ok(Get::Delta(deltas)) => match self.delta_sum(deltas, row) {
                        Err(stop) => Err(stop),
                        Ok(sum) => match self.delta_next(sum) {
                            DeltaNext::Read(sum) => {
                                get = Some((
                                    sum.deltas,
                                    sum.sums.len() as u64,
                                    sum.want.multi_name(),
                                ));
                                let want = sum.want;
                                pending.push(Pending::Delta(sum));
                                pending.push(Pending::Cast(want));
                                continue;
                            }
                            DeltaNext::Sum(sum) => Ok(Some(sum)),
                        },
                    },
                };
            }
            let Some(waiting) = pending.pop() else {
                return result;
            };
            result = match waiting {
                Pending::Scale { mul, add } => result.and_then(|value| {
                    value
                        .map(|v| {
                            v.checked_mul(i64::from(mul))
                                .and_then(|v| v.checked_add(i64::from(add)))
                                .ok_or(Stop::Refuse(Refusal::ScaledOverflow))
                        })
                        .transpose()
                }),
                Pending::Cast(want) => result.and_then(|value| {
                    let Some(value) = value else {
                        return Ok(None);
                    };
                    let cast = want.cast(value);
                    if cast != value {
                        return Err(Stop::exception(NcbiException::new(
                            MULTI_TOO_BIG,
                            format!(
                                "CSeqTable_multi_data::TryGet{}(): value is too big for requested type: {value}",
                                want.multi_name()
                            ),
                        )));
                    }
                    Ok(Some(cast))
                }),
                Pending::Delta(mut sum) => match result.and_then(|value| sum.add(value)) {
                    Err(stop) => {
                        self.caches.int_delta.insert(sum.key, sum.sums);
                        Err(stop)
                    }
                    Ok(()) => match self.delta_next(sum) {
                        DeltaNext::Read(sum) => {
                            get = Some((sum.deltas, sum.sums.len() as u64, sum.want.multi_name()));
                            let want = sum.want;
                            pending.push(Pending::Delta(sum));
                            pending.push(Pending::Cast(want));
                            // (`get` gives the next result.)
                            Ok(None)
                        }
                        DeltaNext::Sum(sum) => Ok(Some(sum)),
                    },
                },
            };
        }
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_multi_data.cpp:81-179,
    // 227-235, 305-316 (the sums of an `int-delta`, `USE_DELTA_CACHE` defined)
    // The value of row `r` reads every element of the blocks of 128 up to `r`'s, in order (as
    // `Int4` with `TryGetInt4` when the elements fit in 4 bytes, else as `Int8`), so the first
    // element that throws there throws for `r`; a sum that leaves the accumulator's type is a
    // signed overflow (refused). `CInt4DeltaSumCache`'s sum is the `Int8` that `GetDeltaSum8`
    // returns (`CIntDeltaSumCache::GetDeltaSum8` calls the virtual `GetDeltaSum4`).
    fn delta_sum<'v>(&mut self, deltas: V<'v>, row: u64) -> Result<DeltaSum<'v>, Stop> {
        const BLOCK: u64 = 128;
        let key = address(deltas);
        let size = self.multi_size(deltas);
        let int4 = self.int_size(deltas) <= 4;
        let want = if int4 { Want::Int4 } else { Want::Int8 };
        let end = ((row / BLOCK + 1) * BLOCK).min(size);
        let sums = self.caches.int_delta.remove(&key).unwrap_or_default();
        Ok(DeltaSum {
            deltas,
            key,
            row,
            want,
            int4,
            end,
            sum: sums.last().copied().unwrap_or(0),
            sums,
        })
    }

    /// `GetSize()` of the multi-data `data` (`multi_size_step`), found once for each value: NCBI
    /// walks an `int-delta`'s, `int-scaled`'s or `real-scaled`'s data to the bottom again at each
    /// of their levels where it asks (`x_TryGetInt8`, the sums' cache), the same size each time.
    fn multi_size(&mut self, data: V<'_>) -> u64 {
        Self::chain_size(&mut self.caches.size, data, multi_size_step)
    }

    /// `GetIntSize()` of the multi-data `data` (`int_size_step`), found once for each value (as
    /// `multi_size`).
    fn int_size(&mut self, data: V<'_>) -> u64 {
        Self::chain_size(&mut self.caches.int_size, data, int_size_step)
    }

    /// The size that `step` gives at the bottom of the values it passes to from `data` (its
    /// `Err`), from `cache` where a value of the way has one, kept for each value of the way.
    fn chain_size<'v>(
        cache: &mut HashMap<usize, u64>,
        data: V<'v>,
        step: fn(V<'v>) -> Result<u64, V<'v>>,
    ) -> u64 {
        let mut way = Vec::new();
        let mut at = data;
        let size = loop {
            let key = address(at);
            if let Some(&size) = cache.get(&key) {
                break size;
            }
            way.push(key);
            match step(at) {
                Ok(size) => break size,
                Err(inner) => at = inner,
            }
        };
        for key in way {
            cache.insert(key, size);
        }
        size
    }

    /// `delta_sum`'s next delta to read (`DeltaSum::add`), or the row's sum once the sums reach
    /// the end of its block (cached again).
    fn delta_next<'v>(&mut self, sum: DeltaSum<'v>) -> DeltaNext<'v> {
        if (sum.sums.len() as u64) < sum.end {
            return DeltaNext::Read(sum);
        }
        let value = sum.sums.get(sum.row as usize).copied().unwrap_or_default();
        self.caches.int_delta.insert(sum.key, sum.sums);
        DeltaNext::Sum(value)
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_sparse_index.cpp:533-580
    // NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_sparse_index.cpp:377-405
    // (`CIndexDeltaSumCache::FindDeltaSum`), 306-360 (`x_FindDeltaSum2`; `USE_DELTA_CACHE` is not
    // defined, line 232), 180-225 (`x_GetBitSetCache`); include/objects/seqtable/
    // SeqTable_sparse_index.hpp:96-97 (`kSkipped` and `kInvalidRow` are both the largest `size_t`).
    // `x_FindDeltaSum` of `row` is the first index whose sum of the deltas up to it reaches `row`,
    // when that sum is `row`, else `kInvalidRow`. The deltas are `unsigned` (`TIndexes_delta`), so
    // the sums never decrease, and the blocks that NCBI's cache fills give that index whatever rows
    // were asked before. `x_GetBitSetCache` of `b` is the number of set bits in the bytes before
    // byte `b`. The `indexes` are `unsigned` in the order of the data, which nothing checks: the
    // index found is the one NCBI's `lower_bound` finds (`utils::lower_bound`).
    // A `bit-set-bvector` is BitMagic data, decoded at its first index (`Indexer::bit_vector`).
    /// `CSeqTable_sparse_index::GetIndexAt` of `row`: the index in the column's data of the row's
    /// value, `None` for a row without one (`kSkipped`).
    fn sparse_index_at(&mut self, sparse: V<'_>, row: u64) -> Result<Option<u64>, Stop> {
        let (variant, value) = sparse.choice();
        Ok(match variant {
            "indexes" => {
                // The `unsigned` elements compared with the `size_t` row.
                let len = value.len();
                let at = lower_bound_by(0, len, |i| (value.int_in(i) as u64) < row);
                (at < len && value.int_in(at) as u64 == row).then_some(at as u64)
            }
            "indexes-delta" => {
                let sums = self
                    .caches
                    .index_delta
                    .entry(address(value))
                    .or_insert_with(|| {
                        let mut sum = 0u64;
                        value
                            .ints()
                            .map(|delta| {
                                sum = sum.wrapping_add(delta as u64);
                                sum
                            })
                            .collect()
                    });
                let at = lower_bound_by(0, sums.len(), |i| sums[i] < row);
                (at < sums.len() && sums[at] == row).then_some(at as u64)
            }
            "bit-set" => {
                let bytes = value.octets();
                let Some(byte_index) = usize::try_from(row / 8).ok().filter(|&i| i < bytes.len())
                else {
                    return Ok(None);
                };
                let byte = u32::from(bytes[byte_index]);
                // The most significant bit is row 0 of the byte.
                let bit = (row % 8) as u32;
                if (byte << bit) & 0x80 == 0 {
                    return Ok(None);
                }
                let before = self
                    .caches
                    .bit_set
                    .entry(address(value))
                    .or_insert_with(|| {
                        let mut count = 0u64;
                        bytes
                            .iter()
                            .map(|&b| {
                                let before = count;
                                count += u64::from(b.count_ones());
                                before
                            })
                            .collect()
                    });
                Some(before[byte_index] + u64::from((byte >> (8 - bit)).count_ones()))
            }
            "bit-set-bvector" => {
                let vector = self.bit_vector(value)?;
                if row >= u64::from(vector.size())
                    || !vector.get_bit(row as u32).map_err(query_stop)?
                {
                    return Ok(None);
                }
                Some(match row {
                    0 => 0,
                    _ => vector.count_to(row as u32 - 1).map_err(query_stop)?,
                })
            }
            _ => None,
        })
    }

    // NCBI reference (598d8ae6): c++/include/objects/seqtable/BVector_data.hpp:74-80
    // `x_CreateBitVector` (BVector_data.cpp:73-82, quoted in `bitmagic`) decodes the data into a
    // vector of the `size` member; a throw leaves no vector, and the next value decodes the data
    // again and throws again (`bm/SPEC.md` §0), as the cached failure here does.
    /// `GetBitVector()` of the `BVector-data` `value`: its vector (`bitmagic::decode`, once per
    /// `BVector-data`), NCBI's `std::logic_error` (`BM::Invalid serialization format`), or LOSAT's
    /// refusal of data that NCBI decodes outside BitMagic's checked form.
    fn bit_vector(&mut self, value: V<'_>) -> Result<&BitVector, Stop> {
        let decoded = self
            .caches
            .bvector
            .entry(address(value))
            .or_insert_with(|| {
                bitmagic::decode(
                    value.field("size").int() as u32,
                    value.field("data").octets(),
                )
            });
        match decoded {
            Ok(vector) => Ok(vector),
            Err(Fail::Invalid) => Err(Stop::std_exception(bitmagic::INVALID)),
            Err(Fail::Refuse(_)) => Err(Stop::Refuse(Refusal::BitVector)),
        }
    }
}

/// The refusal of a query (`BitVector::get_bit`, `count_to`) that reaches a block that NCBI keeps
/// as a GAP block that is not one, at the row that NCBI's index asks: NCBI reads such a block only
/// there, and its answer may depend on the build of BitMagic (`bitmagic::Rule::GapSearch`).
fn query_stop(fail: Fail) -> Stop {
    match fail {
        Fail::Invalid => Stop::std_exception(bitmagic::INVALID),
        Fail::Refuse(_) => Stop::Refuse(Refusal::BitVector),
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_multi_data.cpp:388-431
// NCBI reference (598d8ae6): datatool's include/objects/seqtable/BVector_data_.hpp:76, 97 (the generated
// getter of `BVector-data.size`, `unsigned` by seqtable.def:38-39)
// The size of a bit-bvector is its `size` member (a Uint4), read without decoding BitMagic's data
// (`bm/` probes p01, p03: NCBI is silent where only the size is asked); only a value decodes
// (`x_try_get_int8`, `sparse_index_at`, refused).
/// One level of `GetSize()`: the size, or the multi-data whose size it is (`Err`).
fn multi_size_step(data: V<'_>) -> Result<u64, V<'_>> {
    let (variant, value) = data.choice();
    Ok(match variant {
        "common-string" | "common-bytes" => value.field("indexes").len() as u64,
        "bit" => value.octets().len() as u64 * 8,
        "bit-bvector" => value.field("size").int() as u32 as u64,
        "int-delta" => return Err(value),
        "int-scaled" | "real-scaled" => return Err(value.field("data")),
        _ => value.len() as u64,
    })
}

// NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_multi_data.cpp:363-386
// NCBI reference (598d8ae6): c++/src/objects/seqtable/Scaled_int_multi_data.cpp:57-64
/// One level of `GetIntSize()`: the size, or the multi-data whose size it is (`Err`).
fn int_size_step(data: V<'_>) -> Result<u64, V<'_>> {
    let (variant, value) = data.choice();
    Ok(match variant {
        "int-delta" => return Err(value),
        "int" => 4,
        "int1" | "bit" | "bit-bvector" => 1,
        "int2" => 2,
        "int8" => 8,
        "int-scaled" => {
            let max = value.get("max").map(V::int);
            let min = value.get("min").map(V::int);
            if max.is_some_and(|max| max > i64::from(i32::MAX))
                || min.is_some_and(|min| min < i64::from(i32::MIN))
            {
                8
            } else {
                4
            }
        }
        _ => 0,
    })
}

// NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_single_data.cpp:85-180
// (`sx_DownCast` throws `ThrowOverflowError`, with the value and the type's name, when the value
// does not fit, quoted at `SINGLE_TOO_BIG`.)
/// `CSeqTable_single_data::GetValue` of the type `want` (a bool as 0 or 1).
fn single_get(single: V<'_>, want: Want) -> Result<i64, Stop> {
    let (variant, value) = single.choice();
    let number = match variant {
        "bit" => {
            return Ok(match value.value {
                crate::report::serial::asn_types::Value::Bool(bit) => i64::from(*bit),
                other => unreachable!("a BOOLEAN of the table, not {other:?}"),
            })
        }
        "int" | "int8" => value.int(),
        _ => {
            return Err(Stop::exception(NcbiException::new(
                SINGLE_CANNOT_CONVERT,
                format!(
                    "CSeqTable_single_data::GetValue({}&): {variant} value cannot be converted to requested type",
                    want.single_name()
                ),
            )))
        }
    };
    // An `int` is read as `int` without a check (the table's `int` is an Int4).
    let cast = want.cast(number);
    if cast != number {
        return Err(Stop::exception(NcbiException::new(
            SINGLE_TOO_BIG,
            format!(
                "CSeqTable_single_data::GetValue({}&): value is too big for requested type: {number}",
                want.single_name()
            ),
        )));
    }
    Ok(cast)
}

// NCBI reference (598d8ae6): src/objects/seqtable/SeqTable_single_data_.cpp (datatool's
// `GetLoc()` and `GetId()`: `CheckSelected` for `e_Loc` throws `ThrowInvalidSelection`; the place is
// `SINGLE_SELECTION`'s)
/// `CSeqTable_single_data::Get<getter>()`: the value of the variant `variant`.
fn single_object<'a>(single: V<'a>, variant: &str, getter: &str) -> Result<V<'a>, Stop> {
    let (current, value) = single.choice();
    if current != variant {
        return Err(Stop::exception(NcbiException::new(
            SINGLE_SELECTION,
            format!(
                "CSeqTable_single_data::Get{getter}(): Invalid choice selection: NCBI-SeqTable::SeqTable-single-data.{current}"
            ),
        )));
    }
    Ok(value)
}

/// `GetDefault().GetId()` of a location's id column (`ParseDefaults`).
fn single_get_id(default: V<'_>) -> Result<V<'_>, Stop> {
    single_object(default, "id", "Id")
}

#[cfg(test)]
#[path = "table_tests.rs"]
mod table_tests;

/// What `x_try_get_int8` gives for one multi-data value.
enum Get<'v> {
    Value(Option<i64>),
    /// `int-scaled`: its `data`, then the value times `GetMul()` plus `GetAdd()`.
    Scaled {
        data: V<'v>,
        mul: i32,
        add: i32,
    },
    /// `int-delta`: the sum of its deltas up to the row (`delta_sum`).
    Delta(V<'v>),
}

/// A multi-data value waiting for one inside it (`multi_try_get`).
enum Pending<'v> {
    /// `int-scaled`'s value times `GetMul()` plus `GetAdd()`; a value that leaves `Int8` is a
    /// signed overflow (refused).
    Scale { mul: i32, add: i32 },
    /// `TryGetValue`'s check that the value fits the type asked for.
    Cast(Want),
    /// `delta_sum`'s loop, waiting for the delta at `sums.len()`.
    Delta(DeltaSum<'v>),
}

/// The sums of an `int-delta`'s deltas being made (`delta_sum`).
struct DeltaSum<'v> {
    deltas: V<'v>,
    key: usize,
    row: u64,
    want: Want,
    int4: bool,
    end: u64,
    sum: i64,
    sums: Vec<i64>,
}

impl DeltaSum<'_> {
    /// The delta read (`None`: no value) added to the sums; a sum that leaves the accumulator's
    /// type is a signed overflow (refused).
    fn add(&mut self, value: Option<i64>) -> Result<(), Stop> {
        if let Some(value) = value {
            self.sum = match self.sum.checked_add(value) {
                Some(next) if !self.int4 || i32::try_from(next).is_ok() => next,
                _ => return Err(Stop::Refuse(Refusal::DeltaOverflow)),
            };
        }
        self.sums.push(self.sum);
        Ok(())
    }
}

enum DeltaNext<'v> {
    Read(DeltaSum<'v>),
    Sum(i64),
}
