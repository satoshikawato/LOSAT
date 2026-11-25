//! The numeric id list files of the alias files' `GILIST`, `TILIST` and `TAXIDLIST` keys, read as
//! NCBI's `SeqDB_ReadMemoryGiList`, `SeqDB_ReadMemoryTiList` and `SeqDB_ReadMemoryTaxIdList` read
//! them: a binary list (its first byte 0xFF: a marker word, the big-endian count, then one
//! big-endian id per word) or a text list (decimal numbers, each ended by a space, a line feed or
//! a carriage return). The one owner of their layout, and of the lists the alias files attach to
//! the volumes they reach for NCBI's filtered header read (`AliasIdList`), and of NCBI's lookups of
//! a Seq-id or taxid in a list (`find_id`, `negative_find_id`, `defline_matches_taxids`).

use std::collections::BTreeSet;

use crate::ncbi_exception::NcbiException;
use crate::seqid::defline_set::BlastDefline;
use crate::seqid::model::ObjectId;
use crate::seqid::ParsedSeqId;

/// The four keys of an alias file that name an id list (`CSeqDB_AliasMask`'s `eGiList`,
/// `eTiList`, `eSiList`, `eTaxIdList`), in the order NCBI's `ComputeMasks` reads them.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum IdListKind {
    Gi,
    Ti,
    Si,
    TaxId,
}

impl IdListKind {
    /// Every kind, in NCBI's order.
    pub(crate) const ALL: [IdListKind; 4] = [
        IdListKind::Gi,
        IdListKind::Ti,
        IdListKind::Si,
        IdListKind::TaxId,
    ];

    /// The alias file's key.
    pub(crate) fn key(self) -> &'static str {
        match self {
            IdListKind::Gi => "GILIST",
            IdListKind::Ti => "TILIST",
            IdListKind::Si => "SEQIDLIST",
            IdListKind::TaxId => "TAXIDLIST",
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbalias.cpp:2016-2090
    /// The words of NCBI's error for a value with a space: `has multiple <words> (<value>).`.
    pub(crate) fn multiple_words(self) -> &'static str {
        match self {
            IdListKind::Gi => "GI lists",
            IdListKind::Ti => "TI lists",
            IdListKind::Si => "SEQID lists",
            IdListKind::TaxId => "Tax ids lists",
        }
    }
}

/// Compact storage for Seq-id byte strings (`AliasIdList::sis`, version 5 alias `SEQIDLIST`
/// filtering): one shared contiguous byte buffer plus `(start, end)` spans, avoiding one heap
/// allocation per id.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct PackedStrings {
    data: std::sync::Arc<[u8]>,
    spans: Vec<(usize, usize)>,
}

impl PackedStrings {
    pub(crate) fn from_parts(data: Vec<u8>, mut spans: Vec<(usize, usize)>) -> Self {
        spans.shrink_to_fit();
        Self {
            data: std::sync::Arc::from(data.into_boxed_slice()),
            spans,
        }
    }

    pub(crate) fn from_vecs(items: impl IntoIterator<Item = Vec<u8>>) -> Self {
        let mut data = Vec::new();
        let mut spans = Vec::new();
        for item in items {
            let start = data.len();
            data.extend_from_slice(&item);
            spans.push((start, data.len()));
        }
        Self::from_parts(data, spans)
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.spans.is_empty()
    }

    pub(crate) fn len(&self) -> usize {
        self.spans.len()
    }

    pub(crate) fn get(&self, index: usize) -> &[u8] {
        let (start, end) = self.spans[index];
        &self.data[start..end]
    }

    pub(crate) fn sort_unstable(&mut self) {
        let data = &self.data;
        self.spans
            .sort_unstable_by(|&(s1, e1), &(s2, e2)| data[s1..e1].cmp(&data[s2..e2]));
    }
}

/// Indexed view over Seq-id byte strings (`[Vec<u8>]`, `[&[u8]]`, or [`PackedStrings`]).
pub(crate) trait SeqIdSlice {
    fn len(&self) -> usize;
    fn get(&self, index: usize) -> &[u8];
}

impl SeqIdSlice for [Vec<u8>] {
    fn len(&self) -> usize {
        <[Vec<u8>]>::len(self)
    }

    fn get(&self, index: usize) -> &[u8] {
        self[index].as_slice()
    }
}

impl SeqIdSlice for [&[u8]] {
    fn len(&self) -> usize {
        <[&[u8]]>::len(self)
    }

    fn get(&self, index: usize) -> &[u8] {
        self[index]
    }
}

impl SeqIdSlice for PackedStrings {
    fn len(&self) -> usize {
        PackedStrings::len(self)
    }

    fn get(&self, index: usize) -> &[u8] {
        PackedStrings::get(self, index)
    }
}

/// An id list of the alias files as NCBI attaches it to a volume it reaches
/// (`CSeqDBVol::AttachVolumeGiList`): the ids its filtered header read finds a definition line's
/// Seq-ids in (`CSeqDBGiList::FindId`, [`find_id`]), each kind in NCBI's order, and the taxids of
/// a version 5 `TAXIDLIST` (`s_IncludeDefline_Taxid`, [`defline_matches_taxids`]).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct AliasIdList {
    /// The GIs of a `GILIST`.
    pub gis: Vec<i64>,
    /// The TIs of a `TILIST`, as NCBI's unsigned `TTi` (a number past `Int8` is none).
    pub tis: Vec<u64>,
    /// The Seq-id strings of a `SEQIDLIST`, as the last translation left them, in its order (out
    /// of order where a second simplification changed one after the sort; sorted for a version 5
    /// text list, file order for a version 5 binary list).
    pub sis: PackedStrings,
    /// The taxids of a version 5 `TAXIDLIST` (`GetTaxIdsList`).
    pub taxids: BTreeSet<i32>,
}

impl AliasIdList {
    /// `FindId`: whether the list holds `id` ([`find_id`]).
    pub(crate) fn find_id(&self, id: &ParsedSeqId) -> Result<bool, NcbiException> {
        find_id_in(&self.gis, &self.tis, &self.sis, id)
    }

    /// Whether the list holds taxids (positive `GetNumTaxIds` count).
    pub(crate) fn has_taxids(&self) -> bool {
        !self.taxids.is_empty()
    }

    /// Whether `defline`'s taxid or links hold one of this list's taxids (`s_IncludeDefline_Taxid`).
    pub(crate) fn keeps_taxid(&self, defline: &BlastDefline) -> bool {
        defline_matches_taxids(defline, &self.taxids)
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1910-1941
/// `s_IncludeDefline_Taxid`: whether `defline`'s `taxid` or `links` hold a taxid in `taxids`.
pub(crate) fn defline_matches_taxids(defline: &BlastDefline, taxids: &BTreeSet<i32>) -> bool {
    defline
        .taxid
        .into_iter()
        .chain(defline.links.iter().flatten().copied())
        .any(|taxid| taxids.contains(&taxid))
}

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbvol.hpp:957-968
/// `have_vol` of NCBI's id test of one Seq-id `id` of a definition line: one of `lists` holds it,
/// asked in their order; NCBI's exception of the first list whose lookup throws.
pub(crate) fn lists_hold(
    lists: &[std::sync::Arc<AliasIdList>],
    id: &ParsedSeqId,
) -> Result<bool, NcbiException> {
    for list in lists {
        if list.find_id(id)? {
            return Ok(true);
        }
    }
    Ok(false)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1633-1665
// A `ti` general id's tag is converted whatever the list holds: a string tag by
// `NStr::StringToInt8`, whose exception for a tag that is not a number ends NCBI's read (probe ti
// of line RG D1c-1: blastdbcmd's `Cannot convert string 'abc' to Int8`, rc 1; blastp's exception
// trace, rc 255; probe p3 of review 1 of line RG D1c-1 for a `-seqidlist`, which holds no ti).
/// `CSeqDBGiList::FindId`: whether a positive list of the GIs `gis`, the TIs `tis` (NCBI's
/// unsigned `TTi`) and the Seq-id strings `sis` holds `id`; NCBI's exception of a `ti` id whose
/// string tag is no number, and of the simplification of the Seq-id ([`find_si_forms`]). The
/// three lists are in NCBI's order (`sis` in the order its last translation left it).
pub(crate) fn find_id(
    gis: &[i64],
    tis: &[u64],
    sis: &[Vec<u8>],
    id: &ParsedSeqId,
) -> Result<bool, NcbiException> {
    find_id_in(gis, tis, sis, id)
}

fn find_id_in<S: SeqIdSlice + ?Sized>(
    gis: &[i64],
    tis: &[u64],
    sis: &S,
    id: &ParsedSeqId,
) -> Result<bool, NcbiException> {
    match id {
        ParsedSeqId::Gi(gi) => Ok(gis.binary_search(gi).is_ok()),
        ParsedSeqId::General(dbtag) if dbtag.db == b"ti" => {
            let ti = match &dbtag.tag {
                ObjectId::Id(value) => i64::from(*value),
                ObjectId::Str(text) => crate::ncbi_str::string_to_int8(text, 0)?,
            };
            Ok(tis.binary_search(&(ti as u64)).is_ok())
        }
        _ => find_si_forms_in(sis, id),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1564-1610
// The kind of the id decides which of the list's ids are asked (`match_type`): a gi only the
// GIs, a `ti` general id only the TIs (its tag is converted only then), any other id only the
// Seq-id strings. LOSAT's negative lists hold Seq-id strings alone.
/// `CSeqDBNegativeList::FindId` of a negative list of the Seq-id strings `sis` (sorted):
/// `(found, match_type)`; NCBI's exception of the simplification of the Seq-id
/// ([`find_si_forms`]).
pub(crate) fn negative_find_id(
    sis: &[Vec<u8>],
    id: &ParsedSeqId,
) -> Result<(bool, bool), NcbiException> {
    match id {
        ParsedSeqId::Gi(_) => Ok((false, false)),
        ParsedSeqId::General(dbtag) if dbtag.db == b"ti" => Ok((false, false)),
        _ => {
            let match_type = !sis.is_empty();
            Ok((match_type && find_si_forms(sis, id)?, match_type))
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1645-1662
// The string of the Seq-id with its version and without, then `str_id` of `SeqDB_SimplifySeqid`
// with no accession (`isam::seqdb_simplify_seqid`), with and without what follows its first
// '.': the simplification of a `BL_ORD_ID` or `PIG` general id with a string tag throws
// (`CObject_id::GetId`), as does that of a `ti` general id whose string tag is no number.
/// The lookups of a Seq-id that is not a gi or a ti in the Seq-id strings `sis` of a list
/// (`FindSi`, [`find_si`]); NCBI's exception of the simplification.
pub(crate) fn find_si_forms(sis: &[Vec<u8>], id: &ParsedSeqId) -> Result<bool, NcbiException> {
    find_si_forms_in(sis, id)
}

fn find_si_forms_in<S: SeqIdSlice + ?Sized>(
    sis: &S,
    id: &ParsedSeqId,
) -> Result<bool, NcbiException> {
    if find_si_in(sis, &blast_seq_id_string(id, true))
        || find_si_in(sis, &blast_seq_id_string(id, false))
    {
        return Ok(true);
    }
    let str_id = match super::isam::seqdb_simplify_seqid(id) {
        Ok(simplified) => simplified.str_id,
        Err(super::isam::SimplifyError::Exception(exception)) => return Err(exception),
        Err(super::isam::SimplifyError::SeqId(_)) => {
            unreachable!("a parsed Seq-id is simplified without parsing")
        }
    };
    if find_si_in(sis, &str_id) {
        return Ok(true);
    }
    Ok(match str_id.iter().position(|&c| c == b'.') {
        Some(pos) => find_si_in(sis, &str_id[..pos]),
        None => false,
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:745-767
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1542-1560
// `SiToOid` and `CSeqDBNegativeList::FindSi` halve the range of the ids by `std::string` order
// (bytes) until the middle id is `si`: on ids out of order (an alias file's list that a second
// simplification changed after its sort) the walk is NCBI's, not that of a search of sorted ids.
/// Whether NCBI's binary search of the Seq-id strings `sis`, in their order, finds `si`.
#[cfg(test)]
pub(crate) fn find_si(sis: &[Vec<u8>], si: &[u8]) -> bool {
    find_si_in(sis, si)
}

fn find_si_in<S: SeqIdSlice + ?Sized>(sis: &S, si: &[u8]) -> bool {
    let (mut b, mut e) = (0, sis.len());
    while b < e {
        let m = (b + e) / 2;
        match sis.get(m).cmp(si) {
            std::cmp::Ordering::Less => b = m + 1,
            std::cmp::Ordering::Greater => e = m,
            std::cmp::Ordering::Equal => return true,
        }
    }
    false
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2698-2705
/// `GetBlastSeqIdString` of `id`, with its version or without.
fn blast_seq_id_string(id: &ParsedSeqId, version: bool) -> Vec<u8> {
    let mut out = Vec::new();
    match id {
        ParsedSeqId::Pir(_) | ParsedSeqId::Prf(_) => crate::seqid::label::fasta(id, &mut out),
        _ => crate::seqid::label::content(id, version, &mut out),
    }
    out
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:854-887
// A first byte that is a digit or '#' makes a text list (where '#' is then an invalid byte), a
// first byte 0xFF in a file of 8 bytes or more a binary one; its fourth byte 0xFE or 0xFC names
// 8-byte ids.
/// Whether the bytes of a numeric list are a binary list (`Some`), and of 8-byte ids (its fourth
/// byte); NCBI's `CSeqDBException` text for an empty or unknown file.
fn binary_kind(bytes: &[u8]) -> Result<Option<bool>, String> {
    let Some(&first) = bytes.first() else {
        return Err("Specified file is empty.".to_string());
    };
    if first.is_ascii_digit() || first == b'#' {
        return Ok(None);
    }
    if bytes.len() >= 8 && first == 0xFF {
        return Ok(Some(matches!(bytes[3], 0xFE | 0xFC)));
    }
    Err("Specified file is not a valid GI/TI list.".to_string())
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1376-1403
/// `SeqDB_IsBinaryGiList` of a file's bytes (`s_ContainsBinaryNumericIdList` of `eGiList`):
/// whether they are a binary list, of GIs or TIs alike; NCBI's `CSeqDBException` text for an
/// empty or unknown file.
pub(crate) fn is_binary_gi_list(bytes: &[u8]) -> Result<bool, String> {
    binary_kind(bytes).map(|kind| kind.is_some())
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb.cpp:220-238
/// The first word of a binary GI list of 4-byte GIs, which the reader takes (`read_gi_list`) and
/// the writer writes (`writedb::binary_list`): -1.
pub(crate) const GI_LIST_MAGIC: i32 = -1;

/// The first word of a binary GI list of 8-byte GIs (-2), which the writer writes when a GI does
/// not fit 32 bits and NCBI's reader of GI lists does not take.
pub(crate) const LONG_GI_LIST_MAGIC: i32 = -2;

/// The big-endian 32-bit word `at` of `bytes`.
fn word(bytes: &[u8], at: usize) -> u32 {
    u32::from_be_bytes(bytes[4 * at..4 * at + 4].try_into().expect("4 bytes"))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:889-922
/// The digit of a text list's byte, `None` for a separator (a space, a line feed, a carriage
/// return); NCBI's `CSeqDBException` text for any other byte.
fn read_digit(byte: u8, list_type: &str) -> Result<Option<u8>, String> {
    match byte {
        b'0'..=b'9' => Ok(Some(byte - b'0')),
        b' ' | b'\n' | b'\r' => Ok(None),
        _ => Err(format!("Invalid byte in text{list_type} list [{byte}].")),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:983-1001
// The loop keeps a number at the separator after it and none after the last byte: a list whose
// last number has no separator after it loses that number. A 0 is no number.
/// The numbers of a text list in file order, `step` adding a digit to the number read so far.
fn read_text<T: Copy + Default + PartialEq>(
    bytes: &[u8],
    list_type: &str,
    step: impl Fn(T, u8) -> T,
) -> Result<Vec<T>, String> {
    let mut out = Vec::new();
    let mut elem = T::default();
    for &byte in bytes {
        match read_digit(byte, list_type)? {
            Some(digit) => elem = step(elem, digit),
            None => {
                if elem != T::default() {
                    out.push(elem);
                }
                elem = T::default();
            }
        }
    }
    Ok(out)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:924-1002
// The binary list's words are counted with the file's length rounded down to whole words; its
// first word is 0xFFFFFFFF and its second the count of the words after the two. A text list's
// number is an unsigned 32-bit value that wraps past 2^32 - 1, as NCBI's `Uint4` does.
/// The GIs of a GI list file (`SeqDB_ReadMemoryGiList`), in file order; NCBI's
/// `CSeqDBException` text for a file that is not one.
pub(crate) fn read_gi_list(bytes: &[u8]) -> Result<Vec<u32>, String> {
    if binary_kind(bytes)?.is_some() {
        let words = bytes.len() / 4;
        let num_gis = words - 2;
        if word(bytes, 0) != GI_LIST_MAGIC as u32 || word(bytes, 1) as usize != num_gis {
            return Err("Specified file is not a valid binary GI file.".to_string());
        }
        return Ok((2..words).map(|at| word(bytes, at)).collect());
    }
    read_text(bytes, "GI", |elem: u32, digit| {
        elem.wrapping_mul(10).wrapping_add(u32::from(digit))
    })
}

/// A number of a TI or taxid list that NCBI's signed accumulator cannot hold: NCBI's overflow is
/// undefined, and LOSAT keeps the number it reads past the limit (widened), which no TI or taxid
/// of a database is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum ListNumber {
    /// A number in the accumulator's range.
    Value(i64),
    /// A number past it.
    Past,
}

impl Default for ListNumber {
    fn default() -> Self {
        ListNumber::Value(0)
    }
}

/// `elem * 10 + digit` of NCBI's signed accumulator of `max` (`Int4` or `Int8`), widened.
fn widened_step(elem: ListNumber, digit: u8, max: i64) -> ListNumber {
    match elem {
        ListNumber::Value(value) => value
            .checked_mul(10)
            .and_then(|value| value.checked_add(i64::from(digit)))
            .filter(|&value| value <= max)
            .map_or(ListNumber::Past, ListNumber::Value),
        ListNumber::Past => ListNumber::Past,
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1149-1282
// A binary TI list's first word is -3 (4-byte TIs) or -4 (8-byte TIs, whose count is that of the
// 8-byte TIs and whose words must be even); its count is compared as an unsigned value. A text
// list's number is an `Int8`.
/// The TIs of a TI list file (`SeqDB_ReadMemoryTiList`), in file order; NCBI's
/// `CSeqDBException` text for a file that is not one.
pub(crate) fn read_ti_list(bytes: &[u8]) -> Result<Vec<ListNumber>, String> {
    if let Some(long_ids) = binary_kind(bytes)? {
        let words = bytes.len() / 4;
        let data = words - 2;
        let num_tis = if long_ids { data / 2 } else { data };
        let marker = word(bytes, 0) as i32;
        let bad = (marker != -3 && marker != -4)
            || word(bytes, 1) as usize != num_tis
            || (data % 2 == 1 && long_ids);
        if bad {
            return Err("Specified file is not a valid binary GI or TI file.".to_string());
        }
        if long_ids {
            return Ok((0..num_tis)
                .map(|at| {
                    let high = u64::from(word(bytes, 2 + 2 * at));
                    let low = u64::from(word(bytes, 3 + 2 * at));
                    ListNumber::Value(((high << 32) | low) as i64)
                })
                .collect());
        }
        return Ok((2..words)
            .map(|at| ListNumber::Value(i64::from(word(bytes, at) as i32)))
            .collect());
    }
    read_text(bytes, "TI", |elem, digit| {
        widened_step(elem, digit, i64::MAX)
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:1080-1121
// A binary list's ids are big-endian `Int4`s; a text list's number is an `Int4` (widened past
// `i32::MAX` to avoid C++ signed overflow UB, `ListNumber::Past`).
/// The taxids of a taxid list file (`SeqDB_ReadMemoryTaxIdList`), as the set NCBI keeps;
/// NCBI's `CSeqDBException` text for a file that is not one.
pub(crate) fn read_tax_id_list(bytes: &[u8]) -> Result<BTreeSet<ListNumber>, String> {
    if binary_kind(bytes)?.is_some() {
        let words = bytes.len() / 4;
        let num_taxids = words - 2;
        if word(bytes, 0) != 0xFFFF_FFFF || word(bytes, 1) as i32 != num_taxids as i32 {
            return Err("Specified file is not a valid binary Tax Id List file.".to_string());
        }
        return Ok((2..words)
            .map(|at| ListNumber::Value(i64::from(word(bytes, at) as i32)))
            .collect());
    }
    Ok(read_text(bytes, "TAXID", |elem, digit| {
        widened_step(elem, digit, i64::from(i32::MAX))
    })?
    .into_iter()
    .collect())
}

#[cfg(test)]
#[path = "id_list_tests.rs"]
mod tests;
