//! A Seq-table annotation as NCBI's object manager takes it: the checks of its columns when the
//! annotation is added (`CSeqTableInfo::x_Initialize`, `table_info`), the value of a column at a
//! row (`CSeqTable_column`, `CSeqTable_multi_data`, `CSeqTable_single_data`,
//! `CSeqTable_sparse_index`) and the index of its rows (`x_InitFeatTableKeys`, `index_table`).

use std::collections::HashMap;

use super::range::{Handle, HandleRangeMap, LiveHandles, Mode, Range};
use super::view::V;
use super::{Event, Indexer, Refusal, Stop};
use crate::ncbi_exception::{NcbiException, Severity, ThrowSite};
use crate::utils::lower_bound::lower_bound_by;

/// A good feature table's `num-rows` above which NCBI's index depends on the machine's memory
/// (`x_InitFeatTable` makes an annotation object of about 135 bytes per row: 10,000,000 rows take
/// 1.37 GB; line RG D3-1 `om/SPEC.md` §8 row 26, §10 Q2).
const ROWS_LIMIT: i32 = 10_000_000;

// The places of NCBI's exceptions (line RG D3-1 `om/SPEC.md` §7.1, NCBI BLAST+ 2.17.0 runs).
// NCBI reference (598d8ae6): c++/src/objmgr/seq_table_info.cpp:351-357 (`SetColumn`)
// ```c++
//     if ( field ) {
//         NCBI_THROW_FMT(CAnnotException, eBadLocation,
//                        "Duplicate "<<m_FieldName<<" column");
//     }
// ```
const SET_COLUMN: ThrowSite = ThrowSite {
    file: "src/objmgr/seq_table_info.cpp",
    line: 355,
    kind: "CAnnotException::eBadLocation",
    function: Some("ncbi::objects::CSeqTableLocColumns::SetColumn"),
};

/// The throws of `CSeqTableLocColumns::ParseDefaults` (seq_table_info.cpp:459-521, quoted at
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

// NCBI reference (598d8ae6): c++/src/serial/serialobject.cpp:300-307 (`ThrowUnassigned`, the getter
// of a member that is not set; the report names no function)
const UNASSIGNED: ThrowSite = ThrowSite {
    file: "src/serial/serialobject.cpp",
    line: 306,
    kind: "CUnassignedMember::eGet",
    function: None,
};

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
// ```c++
//     v = DstInt(value);
//     if ( Int8(v) != value ) {
//         NCBI_THROW_FMT(CSeqTableException, eIncompatibleValueType,
//                        "CSeqTable_multi_data::TryGet"<<type_name<<"(): "
//                        "value is too big for requested type: "<<value);
//     }
// ...
//     default:
//         NCBI_THROW_FMT(CSeqTableException, eIncompatibleValueType,
//                        "CSeqTable_multi_data::TryGet"<<type_name<<"(): "
//                        "value cannot be converted to requested type");
// ```
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
// ```c++
// void CSeqTable_single_data::ThrowOverflowError(Int8 value,
//                                                const char* type_name)
// {
//     NCBI_THROW_FMT(CSeqTableException, eIncompatibleValueType,
//                    "CSeqTable_single_data::GetValue("<<type_name<<"&): "
//                    "value is too big for requested type: "<<value);
// }
//
//
// void CSeqTable_single_data::ThrowConversionError(const char* type_name) const
// {
//     NCBI_THROW_FMT(CSeqTableException, eIncompatibleValueType,
//                    "CSeqTable_single_data::GetValue("<<type_name<<"&): "
//                    <<SelectionName(Which())<<
//                    " value cannot be converted to requested type");
// }
// ```
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

    /// `DstInt(value)`: the value in the type (a bool is `value != 0`).
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

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_table_info.cpp:351-372 (`SetColumn`, quoted at
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

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_table_info.cpp:375-454
    // ```c++
    // bool CSeqTableLocColumns::AddColumn(const CSeqTable_column& column)
    // {
    //     const CSeqTable_column_info& type = column.GetHeader();
    //     if ( type.IsSetField_id() ) {
    //         int field = type.GetField_id() - m_BaseValue +
    //             CSeqTable_column_info::eField_id_location;
    //         if ( field < CSeqTable_column_info::eField_id_location ||
    //              field >= CSeqTable_column_info::eField_id_product ) {
    //             return false;
    //         }
    //         switch ( field ) {
    //         case CSeqTable_column_info::eField_id_location:
    //             SetColumn(m_Loc, column);
    //             return true;
    //         ...  (location-id m_Id, location-gi m_Gi, location-from m_From, location-to m_To,
    //               location-strand m_Strand, the two fuzz limits AddExtraColumn)
    //         default:
    //             break;
    //         }
    //     }
    //     if ( !type.IsSetField_name() ) {
    //         return false;
    //     }
    //
    //     CTempString field(type.GetField_name());
    //     if ( field == m_FieldName ) {
    //         SetColumn(m_Loc, column);
    //         return true;
    //     }
    //     else if ( NStr::StartsWith(field, m_FieldName) &&
    //               field[m_FieldName.size()] == '.' ) {
    //         CTempString extra = field.substr(m_FieldName.size()+1);
    //         if ( extra == "id" || NStr::EndsWith(extra, ".id") ) {
    //             SetColumn(m_Id, column);
    //             return true;
    //         }
    //         else if ( extra == "gi" || NStr::EndsWith(extra, ".gi") ) {
    //             SetColumn(m_Gi, column);
    //             return true;
    //         }
    //         else if ( extra == "pnt.point" || extra == "int.from" ) {
    //             SetColumn(m_From, column);
    //             return true;
    //         }
    //         else if ( extra == "int.to" ) {
    //             SetColumn(m_To, column);
    //             return true;
    //         }
    //         else if ( extra == "strand" ||
    //                   NStr::EndsWith(extra, ".strand") ) {
    //             SetColumn(m_Strand, column);
    //             return true;
    //         }
    //         else if ( extra == "int.fuzz-from.lim" ||
    //                   extra == "pnt.fuzz.lim" ) {
    //             AddExtraColumn(column, new CSeqTableSetLocFuzzFromLim());
    //             return true;
    //         }
    //         else if ( extra == "int.fuzz-to.lim" ) {
    //             AddExtraColumn(column, new CSeqTableSetLocFuzzToLim());
    //             return true;
    //         }
    //     }
    //     return false;
    // }
    // ```
    // `GetField_id() - m_BaseValue` is `int` arithmetic: below `INT_MIN + 10` the product's
    // overflows (undefined behaviour), refused.
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
    // ```c++
    // void CSeqTableLocColumns::ParseDefaults(void)
    // {
    //     if ( !m_Is_set ) {
    //         return;
    //     }
    //     if ( m_Loc ) {
    //         m_Is_real_loc = true;
    //         if ( m_Id || m_Gi || m_From || m_To || m_Strand ||
    //              !m_ExtraColumns.empty() ) {
    //             NCBI_THROW_FMT(CAnnotException, eBadLocation,
    //                            "Conflicting "<<m_FieldName<<" columns");
    //         }
    //         return;
    //     }
    //
    //     if ( !m_Id && !m_Gi ) {
    //         NCBI_THROW_FMT(CAnnotException, eBadLocation,
    //                        "No "<<m_FieldName<<".id column");
    //     }
    //     if ( m_Id && m_Gi ) {
    //         NCBI_THROW_FMT(CAnnotException, eBadLocation,
    //                        "Conflicting "<<m_FieldName<<" columns");
    //     }
    //     if ( m_Id ) {
    //         if ( m_Id->IsSetDefault() ) {
    //             m_DefaultIdHandle =
    //                 CSeq_id_Handle::GetHandle(m_Id->GetDefault().GetId());
    //         }
    //     }
    //     if ( m_Gi ) {
    //         if ( m_Gi->IsSetDefault() ) {
    //             TIntId gi;
    //             m_Gi->GetDefault().GetValue(gi);
    //             m_DefaultIdHandle = CSeq_id_Handle::GetGiHandle(GI_FROM(TIntId, gi));
    //         }
    //     }
    //
    //     if ( m_To ) {
    //         // interval
    //         if ( !m_From ) {
    //             NCBI_THROW_FMT(CAnnotException, eBadLocation,
    //                            "column "<<m_FieldName<<".to without "<<
    //                            m_FieldName<<".from");
    //         }
    //         m_Is_simple_interval = true;
    //     }
    //     else if ( m_From ) {
    //         // point
    //         m_Is_simple_point = true;
    //     }
    //     else {
    //         // whole
    //         if ( m_Strand || !m_ExtraColumns.empty() ) {
    //             NCBI_THROW_FMT(CAnnotException, eBadLocation,
    //                            "extra columns in whole "<<m_FieldName);
    //         }
    //         m_Is_simple_whole = true;
    //     }
    //     if ( m_ExtraColumns.empty() ) {
    //         m_Is_simple = true;
    //     }
    //     else {
    //         m_Is_probably_simple = true;
    //     }
    // }
    // ```
    fn parse_defaults(&mut self, handles: &mut LiveHandles) -> Result<(), Stop> {
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
// ```c++
// bool CSeqTableInfo::IsGoodFeatTable(const CSeq_table& table)
// {
//     if ( !table.IsSetFeat_type() ||
//          table.GetFeat_type() <= CSeqFeatData::e_not_set ||
//          table.GetFeat_type() >= CSeqFeatData::e_MaxChoice ) {
//         return false; // not a feature table
//     }
//     if ( table.IsSetFeat_subtype() &&
//          (table.GetFeat_subtype() <= CSeqFeatData::eSubtype_bad ||
//           table.GetFeat_subtype() >= CSeqFeatData::eSubtype_max) ) {
//         return false; // bad subtype
//     }
//     return true;
// }
// ```
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

impl Indexer<'_> {
    // NCBI reference (598d8ae6): c++/src/objmgr/seq_table_info.cpp:784-947
    // ```c++
    //     ITERATE ( CSeq_table::TColumns, it, feat_table.GetColumns() ) {
    //         const CSeqTable_column& col = **it;
    //         const CSeqTable_column_info& type = col.GetHeader();
    //         ...
    //         if ( type.IsSetField_name() ) {
    //             string name = type.GetField_name();
    //             ...
    //             if ( name == "Seq-table location" ) {
    //                 if ( m_TableLocation ) {
    //                     ERR_POST("Duplicate 'Seq-table location' column");
    //                 }
    //                 try {
    //                     CSeqTableColumnInfo info(col);
    //                     m_TableLocation = info.GetSeq_loc(0);
    //                 }
    //                 NCBI_CATCH("Bad 'Seq-table location' column");
    //             }
    //         }
    //         if ( !IsFeatTable() ) {
    //             continue;
    //         }
    //
    //         if ( m_Location.AddColumn(col) || m_Product.AddColumn(col) ) {
    //             continue;
    //         }
    //         CConstRef<CSeqTableSetFeatField> setter;
    //         if ( type.IsSetField_id() ) {
    //             int id = type.GetField_id();
    //             switch ( id ) {
    //             case CSeqTable_column_info::eField_id_partial:
    //                 if ( m_Partial ) {
    //                     NCBI_THROW_FMT(CAnnotException, eOtherError,
    //                                    "Duplicate partial column");
    //                 }
    //                 m_Partial = CSeqTableColumnInfo(col);
    //                 continue;
    //             case CSeqTable_column_info::eField_id_comment:
    //                 setter = new CSeqTableSetComment();
    //                 break;
    //             ...  (data-imp-key, data-region, ext-type: setters)
    //             case CSeqTable_column_info::eField_id_ext:
    //                 setter = new CSeqTableSetExt(type.GetField_name());
    //                 break;
    //             ...  (dbxref, qual: the same with the field name)
    //             default:
    //                 if ( !type.IsSetField_name() ) {
    //                     ERR_POST_X(8, "SeqTable-column-info.field-id = "<<id);
    //                     continue;
    //                 }
    //                 break;
    //             }
    //         }
    //         else if ( !type.IsSetField_name() ) {
    //             ERR_POST_X(9, "SeqTable-column-info: "
    //                        "neither field-id nor field-name is set");
    //             continue;
    //         }
    //         if ( !setter && type.IsSetField_name() ) {
    //             CTempString name(type.GetField_name());
    //             if ( name.empty() ) {
    //                 ERR_POST_X(10, "SeqTable-column-info.field-name is empty");
    //                 continue;
    //             }
    //             else if ( name[0] == 'E' ) {
    //                 setter = new CSeqTableSetExt(name);
    //             }
    //             ...  ('D' dbxref, 'Q' qual)
    //             else if ( name == "partial" ) {
    //                 if ( m_Partial ) {
    //                     NCBI_THROW_FMT(CAnnotException, eOtherError,
    //                                    "Duplicate partial column");
    //                 }
    //                 m_Partial = CSeqTableColumnInfo(col);
    //                 continue;
    //             }
    //             else if ( name == "disabled" ) {
    //                 if ( m_Disabled ) {
    //                     NCBI_THROW_FMT(CAnnotException, eOtherError,
    //                                    "Duplicate disabled column ");
    //                 }
    //                 m_Disabled = CSeqTableColumnInfo(col);
    //                 continue;
    //             }
    //             else if ( name == "Seq-table location" ) {
    //                 // already processed above
    //                 continue;
    //             }
    //             else if ( name == "Sorted, max length" ) {
    //                 if ( m_SortedMaxLength ) {
    //                     ERR_POST("Duplicate 'Sorted, max length' column");
    //                 }
    //                 try {
    //                     CSeqTableColumnInfo info(col);
    //                     int value;
    //                     if ( info.GetValue(0, value) ) {
    //                         m_SortedMaxLength = value;
    //                     }
    //                 }
    //                 NCBI_CATCH("Bad 'Sorted, max length' column");
    //                 continue;
    //             }
    //             if ( !setter ) {
    //                 try {
    //                     setter = new CSeqTableSetAnyFeatField(name);
    //                 }
    //                 catch ( CAnnotException& /*exc*/ ) {
    //                     // ignore invalid column names
    //                 }
    //             }
    //         }
    //         ...
    //     }
    //
    //     if ( IsFeatTable() ) {
    //         m_Location.ParseDefaults();
    //         m_Product.ParseDefaults();
    //     }
    //
    //     m_IsSorted = x_IsSorted();
    // ```
    // The maps of the columns by id and by name and the setters of the feature fields make no
    // message (the setters' constructors only keep the name: `CSeqTableSetExt` and the others
    // cut it with `CTempString::substr`, `CSeqTableSetAnyFeatField` throws only the
    // `CAnnotException`s caught here, seq_table_setters.cpp:180-188, 344-430).
    /// `new CSeqTableInfo(table)` (`x_InitFeatTable`, seq_annot_info.cpp:474-500): the checks of
    /// the columns, with their messages, then LOSAT's bound of the rows that NCBI makes an object
    /// for.
    pub(super) fn table_info<'a>(&mut self, table: V<'a>) -> Result<TableInfo<'a>, Stop> {
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
            let column = Column(column);
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
        // NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:486-498 (`x_InitFeatTable`)
        // ```c++
        //         if ( IsSortedTable() ) {
        //             ...
        //             m_ObjectIndex.AddInfo(CAnnotObject_Info(*this, kWholeAnnotIndex, type));
        //             return;
        //         }
        //
        //         // index each row separately
        //         TAnnotIndex rows = table.GetNum_rows();
        //         for ( TAnnotIndex index = 0; index < rows; ++index ) {
        //             m_ObjectIndex.AddInfo(CAnnotObject_Info(*this, index, type));
        //         }
        // ```
        if is_feat && !info.sorted && info.num_rows > ROWS_LIMIT {
            return Err(Stop::Refuse(Refusal::TableRows));
        }
        Ok(info)
    }

    /// `NCBI_CATCH(title)` of `e` (`NcbiException::report_default`): one error post per line.
    fn ncbi_catch(&mut self, e: &NcbiException, title: &str) {
        let lines = e.report_default(self.program, title);
        for line in lines.split_inclusive(|&byte| byte == b'\n') {
            self.posts.post(crate::report::diag_log::ERROR, line);
        }
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:1028-1170 (`x_InitFeatTableKeys`)
    // ```c++
    //     if ( !m_Table_Info->IsFeatTable() ) {
    //         // index whole Seq-table
    //         ...
    //         CConstRef<CSeq_loc> loc = m_Table_Info->GetTableLocation();
    //         if ( !loc ) {
    //             return;
    //         }
    //         hrmap.AddLocation(*loc);
    //         ...
    //         ITERATE ( CHandleRangeMap, hrit, hrmap ) {
    //             const CHandleRange& hr = hrit->second;
    //             key.m_Range = hr.GetOverlappingRange();
    //             if ( key.m_Range.Empty() ) {
    //                 ERR_POST_X(7, "Empty region in "<<GetDescription()<<
    //                            " Seq-table location "<<MSerial_AsnText<<*loc);
    //                 continue;
    //             }
    //             ...  (the circular location's two ranges, else the range, x_Map)
    //         }
    //         ...
    //         return;
    //     }
    //     ...
    //     for ( size_t row = 0; row < object_count; ++row, ++it ) {
    //         CAnnotObject_Info& info = *it;
    //         if ( info.GetAnnotIndex() == kWholeAnnotIndex ) {
    //             ...
    //             CConstRef<CSeq_loc> loc = m_Table_Info->GetTableLocation();
    //             _ASSERT(loc && loc->IsInt()); // checked in x_InitFeatTable()
    //             key.m_Handle = CSeq_id_Handle::GetHandle(loc->GetInt().GetId());
    //             key.m_Range = loc->GetTotalRange();
    //             x_Map(mapper, key, index);
    //             x_UpdateObjectKeys(info, keys_begin);
    //             continue;
    //         }
    //         ...
    //         if ( m_Table_Info->RowIsDisabled(row) ) {
    //             continue;
    //         }
    //         ...
    //         bool partial = m_Table_Info->IsPartial(row);
    //         for ( index.m_AnnotLocationIndex = 0;
    //               index.m_AnnotLocationIndex < 2;
    //               index.m_AnnotLocationIndex++ ) {
    //             // index by location, then by product
    //             const CSeqTableLocColumns& loc = ... GetLocation() : GetProduct();
    //             if ( loc.IsSet() ) { // if this Seq-loc field is set
    //                 if ( loc.IsRealLoc() ) { // full Seq-loc object
    //                     CHandleRangeMap hrmap;
    //                     hrmap.SetMasterSeq(master);
    //                     hrmap.AddLocation(*loc.GetLoc(row));
    //                     ...
    //                     ITERATE ( CHandleRangeMap, hrit, hrmap ) {
    //                         const CHandleRange& hr = hrit->second;
    //                         key.m_Range = hr.GetOverlappingRange();
    //                         if ( key.m_Range.Empty() ) {
    //                             ERR_POST_X(8, "Empty region in "<<GetDescription()<<
    //                                        " Seq-table row "<<row);
    //                             continue;
    //                         }
    //                         ...  (the circular location's two ranges, else the range, x_Map)
    //                     }
    //                 }
    //                 else { // simplified Seq-loc object
    //                     loc.SetTableKeyAndIndex(row, key, index);
    //                     index.SetPartial(partial);
    //                     x_Map(mapper, key, index);
    //                 }
    //             }
    //         }
    //         ...
    //     }
    // ```
    /// `x_InitFeatTableKeys` of the table `info`.
    pub(super) fn index_table(&mut self, info: &TableInfo<'_>) -> Result<(), Stop> {
        if !info.is_feat {
            if let Some(loc) = info.table_loc {
                let map = HandleRangeMap::of_location(loc, Mode::None, &mut self.handles);
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
            // ```c++
            //     bool GetBool(size_t row) const
            //         {
            //             bool val = false;
            //             return *this && m_Column->TryGetBool(row, val) && val;
            //         }
            // ...
            //     bool RowIsDisabled(size_t row) const {
            //         return m_Disabled.GetBool(row);
            //     }
            //     bool IsPartial(size_t row) const {
            //         return m_Partial.GetBool(row);
            //     }
            // ```
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
                    // `*loc.GetLoc(row)` of a row without a location dereferences a null
                    // `CConstRef` (`CCoreException::eNullPtr`, a report of memory addresses).
                    let Some(loc) = self.column_seq_loc(column, row)? else {
                        return Err(self.null_pointer());
                    };
                    let map = HandleRangeMap::of_location(loc, Mode::None, &mut self.handles);
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
                // NCBI's row map is made for each location (seq_annot_info.cpp:1130).
                self.handles.end_object();
            }
        }
        self.handles.end_table();
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_table_info.cpp:701-746 (`SetTableKeyAndIndex`)
    // ```c++
    //     key.m_Handle = GetIdHandle(row);
    //     key.m_Range = GetRange(row);
    //     ENa_strand strand = GetStrand(row);
    //     ...
    //     bool simple = m_Is_simple;
    //     if ( !simple && m_Is_probably_simple ) {
    //         simple = true;
    //         ITERATE ( TExtraColumns, it, m_ExtraColumns ) {
    //             if ( it->first.IsSet(row) ) {
    //                 simple = false;
    //                 break;
    //             }
    //         }
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objmgr/seq_table_info.cpp:538-603
    // ```c++
    // CSeq_id_Handle CSeqTableLocColumns::GetIdHandle(size_t row) const
    // {
    //     _ASSERT(!m_Loc);
    //     if ( m_Id ) {
    //         _ASSERT(!m_Id->IsSetSparse());
    //         if ( m_Id->IsSetData() ) {
    //             const CSeq_id* id = m_Id.GetSeq_id(row);
    //             if ( id ) {
    //                 return CSeq_id_Handle::GetHandle(*id);
    //             }
    //         }
    //     }
    //     else {
    //         _ASSERT(!m_Gi->IsSetSparse());
    //         if ( m_Gi->IsSetData() ) {
    //             TIntId gi;
    //             if ( m_Gi.GetValue(row, gi) ) {
    //                 return CSeq_id_Handle::GetGiHandle(GI_FROM(TIntId, gi));
    //             }
    //         }
    //     }
    //     return m_DefaultIdHandle;
    // }
    // ...
    // CRange<TSeqPos> CSeqTableLocColumns::GetRange(size_t row) const
    // {
    //     _ASSERT(!m_Loc);
    //     _ASSERT(m_From);
    //     int from;
    //     if ( !m_From || !m_From.GetValue(row, from) ) {
    //         return CRange<TSeqPos>::GetWhole();
    //     }
    //     int to = from;
    //     if ( m_To ) {
    //         m_To.GetValue(row, to);
    //     }
    //     return CRange<TSeqPos>(from, to);
    // }
    //
    //
    // ENa_strand CSeqTableLocColumns::GetStrand(size_t row) const
    // {
    //     _ASSERT(!m_Loc);
    //     int strand = eNa_strand_unknown;
    //     if ( m_Strand ) {
    //         m_Strand.GetValue(row, strand);
    //     }
    //     return ENa_strand(strand);
    // }
    // ```
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
    // NCBI reference (598d8ae6): c++/src/objmgr/seq_table_info.cpp:972-1004
    // ```c++
    // bool CSeqTableInfo::x_IsSorted(void) const
    // {
    //     ...
    //     if ( m_Product.IsSet() ) {
    //         return false;
    //     }
    //     if ( !m_Location.IsSet() || m_Location.IsRealLoc() || !m_Location.m_Id ) {
    //         return false;
    //     }
    //     if ( !m_Location.m_Id.IsSingular() ) {
    //         return false;
    //     }
    //     if ( !m_Location.m_Is_simple ||
    //          !(m_Location.m_Is_simple_point || m_Location.m_Is_simple_interval) ) {
    //         return false;
    //     }
    //     if ( !m_TableLocation || !m_TableLocation->IsInt() ) {
    //         return false;
    //     }
    //     if ( !m_SortedMaxLength ||
    //          m_SortedMaxLength > m_TableLocation->GetInt().GetLength()/16 ) {
    //         return false;
    //     }
    //     return true;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_interval.cpp:125-128
    // ```c++
    // TSeqPos CSeq_interval::GetLength(void) const
    // {
    //     return GetTo() >= GetFrom() ? ((GetTo() - GetFrom()) + 1) : 0;
    // }
    // ```
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
    // NCBI reference (598d8ae6): c++/include/objects/seqtable/SeqTable_column.hpp:76-96
    // ```c++
    //     template<class Value>
    //     bool TryGetValue(size_t row, Value& v) const {
    //         size_t index = row;
    //         if ( IsSetSparse() ) {
    //             index = GetSparse().GetIndexAt(row);
    //             if ( index == CSeqTable_sparse_index::kSkipped ) {
    //                 if ( !IsSetSparse_other() ) {
    //                     return false;
    //                 }
    //                 GetSparse_other().GetValue(v);
    //                 return true;
    //             }
    //         }
    //         if ( IsSetData() && GetData().TryGetValue(index, v) ) {
    //             return true;
    //         }
    //         if ( IsSetDefault() ) {
    //             GetDefault().GetValue(v);
    //             return true;
    //         }
    //         return false;
    //     }
    // ```
    /// `CSeqTable_column::TryGetValue(row, v)` of the type `want` (a bool as 0 or 1).
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
    // ```c++
    // bool CSeqTable_column::IsSet(size_t row) const
    // {
    //     size_t index = row;
    //     if ( IsSetSparse() ) {
    //         index = GetSparse().GetIndexAt(row);
    //         if ( index == CSeqTable_sparse_index::kSkipped ) {
    //             return IsSetSparse_other();
    //         }
    //     }
    //     if ( IsSetData() && GetData().IsSet(index) ) {
    //         return true;
    //     }
    //     return IsSetDefault();
    // }
    // ```
    // NCBI reference (598d8ae6): c++/include/objects/seqtable/SeqTable_multi_data.hpp:71-73
    // ```c++
    //     bool IsSet(size_t row_index) const {
    //         return row_index < GetSize();
    //     }
    // ```
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
    // ```c++
    // CConstRef<CSeq_id> CSeqTable_column::GetSeq_id(size_t row) const
    // {
    //     size_t index = row;
    //     if ( IsSetSparse() ) {
    //         index = GetSparse().GetIndexAt(row);
    //         if ( index == CSeqTable_sparse_index::kSkipped ) {
    //             if ( IsSetSparse_other() ) {
    //                 return ConstRef(&GetSparse_other().GetId());
    //             }
    //             return null;
    //         }
    //     }
    //     if ( IsSetData() ) {
    //         const CSeqTable_multi_data::TId& arr = GetData().GetId();
    //         if ( index < arr.size() ) {
    //             return arr[index];
    //         }
    //     }
    //     if ( IsSetDefault() ) {
    //         return ConstRef(&GetDefault().GetId());
    //     }
    //     return null;
    // }
    //
    //
    // CConstRef<CSeq_loc> CSeqTable_column::GetSeq_loc(size_t row) const
    // {
    //     ... (the same with GetSparse_other().GetLoc(), GetData().GetLoc(), GetDefault().GetLoc())
    // }
    // ```
    /// `CSeqTable_column::GetSeq_loc(row)`.
    fn column_seq_loc<'a>(&mut self, column: Column<'a>, row: u64) -> Result<Option<V<'a>>, Stop> {
        self.column_object(column, row, "loc", "Loc")
    }

    /// `CSeqTable_column::GetSeq_id(row)`.
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
            if let Some(element) = usize::try_from(index)
                .ok()
                .and_then(|index| list.item(index))
            {
                return Ok(Some(element));
            }
        }
        match column.default() {
            Some(default) => single_object(default, variant, getter).map(Some),
            None => Ok(None),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_multi_data.cpp:477-561
    // ```c++
    // bool CSeqTable_multi_data::x_TryGetInt8(size_t row, Int8& v,
    //                                         const char* type_name) const
    // {
    //     switch ( Which() ) {
    //     case e_Bit:
    //     {
    //         const TBit& bits = GetBit();
    //         size_t i = row/8;
    //         if ( i >= bits.size() ) {
    //             return false;
    //         }
    //         size_t j = row%8;
    //         Uint1 bb = bits[i];
    //         v = ((bb<<j)&0x80) != 0;
    //         return true;
    //     }
    //     case e_Bit_bvector:
    //         return sx_TryGet(GetBit_bvector().GetBitVector(),
    //                          static_cast<bm::bvector<>::size_type>(row), v);
    //     case e_Int1:
    //         return sx_TryGet(GetInt1(), row, v);
    //     case e_Int2:
    //         return sx_TryGet(GetInt2(), row, v);
    //     case e_Int:
    //         return sx_TryGet(GetInt(), row, v);
    //     case e_Int8:
    //         return sx_TryGet(GetInt8(), row, v);
    //     case e_Int_delta:
    //     {
    //         const CSeqTable_multi_data& deltas = GetInt_delta();
    //         if ( row >= deltas.GetSize() ) {
    //             return false;
    //         }
    //         CMutexGuard guard(sx_PrepareMutex_multi_data);
    //         v = x_GetIntDeltaCache().GetDeltaSum8(deltas, row);
    //         return true;
    //     }
    //     case e_Int_scaled:
    //         return GetInt_scaled().TryGetInt8(row, v);
    //     default:
    //         NCBI_THROW_FMT(CSeqTableException, eIncompatibleValueType,
    //                        "CSeqTable_multi_data::TryGet"<<type_name<<"(): "
    //                        "value cannot be converted to requested type");
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objects/seqtable/Scaled_int_multi_data.cpp:77-84
    // ```c++
    // bool CScaled_int_multi_data::TryGetInt8(size_t row, Int8& v) const
    // {
    //     if ( !GetData().TryGetInt8(row, v) ) {
    //         return false;
    //     }
    //     v = v*GetMul()+GetAdd();
    //     return true;
    // }
    // ```
    // A bit-bvector is BitMagic's serialized vector, which `GetBitVector()` decodes from the raw
    // bytes without a check (BVector_data.cpp:64-72): refused. An `int-scaled` value that leaves
    // `Int8` is a signed overflow: refused.
    /// `x_TryGetInt8(row, v, type_name)` of the multi-data `data`, one level: its value, or the
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
            "bit-bvector" => return Err(Stop::Refuse(Refusal::BitVector)),
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

    // NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_multi_data.cpp:435-475
    // ```c++
    // template<class DstInt>
    // static inline
    // bool sx_DownCastInt8(DstInt& v, const Int8& value, const char* type_name)
    // {
    //     v = DstInt(value);
    //     if ( Int8(v) != value ) {
    //         NCBI_THROW_FMT(... "value is too big for requested type: "<<value);
    //     }
    //     return true;
    // }
    // ...
    // bool CSeqTable_multi_data::TryGetBool(size_t row, bool& v) const
    // {
    //     Int8 value;
    //     return x_TryGetInt8(row, value, "Bool") &&
    //         sx_DownCastInt8(v, value, "Bool");
    // }
    // ...
    // bool CSeqTable_multi_data::TryGetInt4(size_t row, Int4& v) const
    // {
    //     Int8 value;
    //     return x_TryGetInt8(row, value, "Int4") &&
    //         sx_DownCastInt8(v, value, "Int4");
    // }
    // ```
    /// `CSeqTable_multi_data::TryGetValue(row, v)` of the type `want`.
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

    // NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_multi_data.cpp:81-94, 138-169,
    // 227-235, 305-316 (the sums of an `int-delta`, `USE_DELTA_CACHE` defined)
    // ```c++
    // CIntDeltaSumCache& CSeqTable_multi_data::x_GetIntDeltaCache(void) const
    // {
    //     CIntDeltaSumCache* info = m_Cache.GetNCPointerOrNull();
    //     if ( !info ) {
    //         if ( GetInt_delta().GetIntSize() <= sizeof(Int4) ) {
    //             info = new CInt4DeltaSumCache(GetInt_delta().GetSize());
    //         }
    //         else {
    //             info = new CInt8DeltaSumCache(GetInt_delta().GetSize());
    //         }
    //         m_Cache = info;
    //     }
    //     return *info;
    // }
    // ...
    //     if ( block_index != m_CacheBlockIndex ) {
    //         size_t size = deltas.GetSize();
    //         size_t block_pos = block_index*kBlockSize;
    //         ...
    //         size_t block_size = min(kBlockSize, size-block_pos);
    //         ...
    //         TValue sum = block_index == 0? 0: m_Blocks[block_index-1];
    //         for ( size_t i = 0; i < block_size; ++i ) {
    //             TValue v;
    //             if ( deltas.TryGetValue(block_pos+i, v) ) {
    //                 sum += v;
    //             }
    //             m_CacheBlockInfo[i] = sum;
    //         }
    //         ...
    //     }
    //     return m_CacheBlockInfo[block_offset];
    // ...
    // Int4 CInt4DeltaSumCache::GetDeltaSum4(const TDeltas& deltas,
    //                                       size_t index)
    // {
    //     _ASSERT(index < deltas.GetSize());
    //     size_t block_index  = index / kBlockSize;
    //     size_t block_offset = index % kBlockSize;
    //     while ( block_index >= m_BlocksFilled ) {
    //         x_GetDeltaSum2(deltas, m_BlocksFilled, 0);
    //     }
    //     return x_GetDeltaSum2(deltas, block_index, block_offset);
    // }
    // ```
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
    // ```c++
    // size_t CSeqTable_sparse_index::GetIndexAt(size_t row) const
    // {
    //     switch ( Which() ) {
    //     case e_Indexes:
    //     {
    //         const TIndexes& indexes = GetIndexes();
    //         TIndexes::const_iterator iter =
    //             lower_bound(indexes.begin(), indexes.end(), row);
    //         if ( iter != indexes.end() && *iter == row ) {
    //             return iter - indexes.begin();
    //         }
    //         else {
    //             return kSkipped;
    //         }
    //     }
    //     case e_Indexes_delta:
    //         return x_FindDeltaSum(row);
    //     case e_Bit_set:
    //     {
    //         const TBit_set& bytes = GetBit_set();
    //         size_t byte_index = row/8;
    //         if ( byte_index >= bytes.size() ) {
    //             return kSkipped;
    //         }
    //         Uint1 byte = bytes[byte_index];
    //         size_t bit_index = row%8; // most significant bit has index 0
    //         if ( !((byte<<bit_index)&0x80) ) {
    //             return kSkipped;
    //         }
    //         size_t count = sx_CalcByteBitCount(Uint1(byte>>(8-bit_index)));
    //         if ( byte_index ) {
    //             count += x_GetBitSetCache(byte_index);
    //         }
    //         return count;
    //     }
    //     case e_Bit_set_bvector:
    //     {
    //         const bm::bvector<>& bv = GetBit_set_bvector().GetBitVector();
    //         auto bv_row = static_cast<bm::bvector<>::size_type>(row);
    //         if ( row >= bv.size()  ||  !bv.get_bit(bv_row) ) {
    //             return kSkipped;
    //         }
    //         return row == 0 ? 0 : bv.count_range(0, bv_row - 1);
    //     }
    //     default:
    //         return kSkipped;
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_sparse_index.cpp:377-405
    // (`CIndexDeltaSumCache::FindDeltaSum`), 306-360 (`x_FindDeltaSum2`; `USE_DELTA_CACHE` is not
    // defined, line 232), 180-225 (`x_GetBitSetCache`); include/objects/seqtable/
    // SeqTable_sparse_index.hpp:96-97 (`kSkipped` and `kInvalidRow` are both `size_t(0)-1`).
    // `x_FindDeltaSum(row)` is the first index whose sum of the deltas up to it reaches `row`, when
    // that sum is `row`, else `kInvalidRow`. The deltas are `unsigned` (`TIndexes_delta`), so the
    // sums never decrease, and the blocks that NCBI's cache fills give that index whatever rows
    // were asked before. `x_GetBitSetCache(b)` is the number of set bits in the bytes before byte
    // `b`. The `indexes` are `unsigned` in the order of the data, which nothing checks: the index
    // found is the one NCBI's `lower_bound` finds (`utils::lower_bound`).
    /// `CSeqTable_sparse_index::GetIndexAt(row)`: the index in the column's data of the row's
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
            "bit-set-bvector" => return Err(Stop::Refuse(Refusal::BitVector)),
            _ => None,
        })
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_multi_data.cpp:388-429
// ```c++
// size_t CSeqTable_multi_data::GetSize(void) const
// {
//     switch ( Which() ) {
//     case e_Int:
//         return GetInt().size();
//     ...
//     case e_Common_string:
//         return GetCommon_string().GetIndexes().size();
//     case e_Common_bytes:
//         return GetCommon_bytes().GetIndexes().size();
//     case e_Bit:
//         return GetBit().size()*8;
//     ...
//     case e_Bit_bvector:
//         return GetBit_bvector().GetSize();
//     case e_Int_delta:
//         return GetInt_delta().GetSize();
//     case e_Int_scaled:
//         return GetInt_scaled().GetData().GetSize();
//     case e_Real_scaled:
//         return GetReal_scaled().GetData().GetSize();
//     default:
//         break;
//     }
//     return 0;
// }
// ```
// NCBI reference (598d8ae6): c++/include/objects/seqtable/BVector_data_.hpp:76, 97 (the generated
// getter of `BVector-data.size`, `unsigned` by seqtable.def:38-39)
// ```c++
//     typedef unsigned TSize;
//     ...
//     TSize GetSize(void) const;
// ```
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

// NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_multi_data.cpp:363-385
// ```c++
// size_t CSeqTable_multi_data::GetIntSize(void) const
// {
//     switch ( Which() ) {
//     case e_Int:
//         return sizeof(TInt::value_type);
//     case e_Int1:
//         return sizeof(TInt1::value_type);
//     case e_Int2:
//         return sizeof(TInt2::value_type);
//     case e_Int8:
//         return sizeof(TInt8::value_type);
//     case e_Bit:
//     case e_Bit_bvector:
//         return 1; // one char is enough to store bit value
//     case e_Int_delta:
//         return GetInt_delta().GetIntSize();
//     case e_Int_scaled:
//         return GetInt_scaled().GetIntSize();
//     default:
//         return 0;
//     }
// }
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqtable/Scaled_int_multi_data.cpp:57-64
// ```c++
// size_t CScaled_int_multi_data::GetIntSize(void) const
// {
//     if ( (IsSetMax() && GetMax() > kMax_Int) ||
//          (IsSetMin() && GetMin() < kMin_Int) ) {
//         return sizeof(Int8);
//     }
//     return sizeof(int);
// }
// ```
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
// ```c++
// void CSeqTable_single_data::GetValue(bool& v) const
// {
//     static const char* const type_name = "bool";
//     switch ( Which() ) {
//     case e_Bit:
//         v = GetBit();
//         break;
//     case e_Int:
//         sx_DownCast(v, GetInt(), type_name);
//         break;
//     case e_Int8:
//         sx_DownCast(v, GetInt8(), type_name);
//         break;
//     default:
//         ThrowConversionError(type_name);
//     }
// }
// ...
// void CSeqTable_single_data::GetValue(int& v) const
// {
//     static const char* const type_name = "int";
//     switch ( Which() ) {
//     case e_Bit:
//         v = GetBit();
//         break;
//     case e_Int:
//         v = GetInt();
//         break;
//     case e_Int8:
//         sx_DownCast(v, GetInt8(), type_name);
//         break;
//     default:
//         ThrowConversionError(type_name);
//     }
// }
//
//
// void CSeqTable_single_data::GetValue(Int8& v) const
// {
//     switch ( Which() ) {
//     case e_Bit:
//         v = GetBit();
//         break;
//     case e_Int:
//         v = GetInt();
//         break;
//     case e_Int8:
//         v = GetInt8();
//         break;
//     default:
//         ThrowConversionError("Int8");
//     }
// }
// ```
// (`sx_DownCast` throws `ThrowOverflowError(value, type_name)` when the value does not fit, quoted
// at `SINGLE_TOO_BIG`.)
/// `CSeqTable_single_data::GetValue(v)` of the type `want` (a bool as 0 or 1).
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

// NCBI reference (598d8ae6): c++/src/objects/seq/../seqtable/SeqTable_single_data_.cpp (datatool's
// `GetLoc()` and `GetId()`: `CheckSelected(e_Loc)` throws `ThrowInvalidSelection`; the place is
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
    /// `int-scaled`: its `data`, then `v*GetMul()+GetAdd()`.
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
    /// `int-scaled`'s `v*GetMul()+GetAdd()`; a value that leaves `Int8` is a signed overflow
    /// (refused).
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
