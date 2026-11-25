//! A database seen through a list of taxids (line DI-8a: `-taxids`, `-taxidlist`,
//! `-negative_taxids`, `-negative_taxidlist`) or of seqids (DI-8b: `-seqidlist`,
//! `-negative_seqidlist`): the OIDs that NCBI's `CSeqDB` visits when it is opened with a
//! user list or a negative list (`CSeqDBOIDList`), the totals it scans for them
//! (`CSeqDBImpl::x_ScanTotals`) and the definition lines its filtered header read keeps
//! (`CSeqDBVol::x_GetFilteredHeader`).
//!
//! A restricted [`SeqDb`] shows its included OIDs as the ordinal ids `0..k`, in OID order
//! (DI-8 DECISIONS 7): every table of a program that an ordinal id indexes stays dense, and
//! [`SeqDb::oid`] gives the OID that NCBI prints as `gnl|BL_ORD_ID|n`. The OIDs are held as a
//! bitmap with a rank table, as NCBI holds them (`CSeqDBBitVector`), not as a list.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::sync::atomic::Ordering;
use std::sync::Arc;

use super::file::DbFile;
use super::isam;
use super::lmdb_read::{Db, Env, ReadError};
use super::seqidlist::{self, SeqidListFile};
use super::{tax_files, SeqDb, SeqDbError, SeqType};
use crate::ncbi_exception::NcbiException;
use crate::seqid::accguide::SeqIdWarning;
use crate::seqid::{label, BlastDefline, DeflineSet, ParsedSeqId};

/// The taxids of a search's list after NCBI's expansion (`s_GetTaxIDList`): `-taxids` and
/// `-taxidlist` give a positive list, `-negative_taxids` and `-negative_taxidlist` a negative
/// one. An empty set is no list at all (`CSearchDatabase::x_InitializeDb` opens the database
/// without one); [`SeqDb::restrict_taxids`] is not called for it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TaxidList {
    /// Search only the OIDs of these taxids (`CSeqDBGiList` with taxids).
    Positive(BTreeSet<i32>),
    /// Search every OID except those whose taxids all are among these
    /// (`CSeqDBNegativeList` with taxids).
    Negative(BTreeSet<i32>),
}

/// The seqid list of a search as read from its file (`CSeqDBFileGiList` with `fname` and `eSiList`,
/// `seqidlist::read_si_list`): `-seqidlist` gives a positive list, `-negative_seqidlist`
/// a negative one. A file without an id is no list at all (`CSearchDatabase::x_InitializeDb`
/// opens the database without one); [`SeqDb::restrict`] is not called for it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeqidList {
    pub negative: bool,
    pub file: SeqidListFile,
}

/// The list a search restricts its database to: NCBI reads one of them
/// (`CBlastDatabaseArgs::ExtractAlgorithmOptions`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DbList {
    Taxids(TaxidList),
    Seqids(SeqidList),
}

/// A warning NCBI writes while it applies the list, to be written as it is.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ListWarning {
    /// `LOG_POST` of `text` as a `Warning` or a line written to `cerr`: the text alone.
    Log(Vec<u8>),
    /// `ERR_POST` of `text` as a `Warning`: with the diagnostic prefix of the program.
    Posted(String),
    /// `ERR_POST` of `text` at the default severity (an error that is reported before NCBI
    /// throws): with the diagnostic prefix of the program. The text is bytes (a seqid list's
    /// word keeps its bytes: probe R3 of line LIB-B).
    Error(Vec<u8>),
    /// `ERR_POST_X` of NCBI's Seq-id code, with `subcode` and `message` as a `Warning`: with its
    /// error code before the diagnostic prefix of the program.
    SeqId(SeqIdWarning),
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:206-226
const TEXT_LIST_WARNING: &str = "To obtain better run time performance, please run \
     blastdb_aliastool -seqid_file_in <INPUT_FILE_NAME> -seqid_file_out <OUT_FILE_NAME> and use \
     <OUT_FILE_NAME> as the argument to -seqidlist";
const BINARY_LIST_V4_ERROR: &str = "Seqidlist is not in BLAST db v4 format";
const DB_INFO_WARNING: &str = "Seqidlist file db info does not match input db";
const EMPTY_LIST_WARNING: &str = "Empty seqid list";

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:242-244
const V4_TAXID_ERROR: &str = "Taxonomy filtering is not supported in v4 BLAST dbs";

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:47-49
const VOLUME_COUNT_ERROR: &str = "Input db vol size does not match lmdb vol size";

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:65-68
const VOLUME_NAME_ERROR: &str = "Input db vol does not match lmdb vol";

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:356-385
// A search reads the OIDs of NCBI's OID list (the alias files' ranges, a user list) inside the
// iteration range `[m_RestrictBegin, m_RestrictEnd)` (`GetNextOIDChunk` caps every chunk at the
// range's end): those are the ordinal ids here, the rank of an OID counted from the range's
// start. The other readers of the list (`x_CheckOrFindOID`) see the OIDs outside the range too
// (`OidFilter::check_or_find`, `OidFilter::walk`).
/// NCBI's OID list of a database (`CSeqDBOIDList`: bit `oid % 64` of word `oid / 64`, with the
/// number of listed OIDs before every 512th OID) and its iteration range `[begin, end)`
/// (`m_RestrictBegin`, `m_RestrictEnd` once the list is set up). The included OIDs, which a
/// search reads and the ordinal ids number (`count`, `contains`, `rank`, `select`, `run_end`),
/// are the listed OIDs inside the range.
#[derive(Debug)]
pub(crate) struct OidFilter {
    words: Vec<u64>,
    /// `supers[s]`: the listed OIDs below OID `512 * s`.
    supers: Vec<u32>,
    num_oids: u32,
    /// The number of listed OIDs.
    total: u32,
    /// The iteration range, `begin <= end <= num_oids`.
    begin: u32,
    end: u32,
    /// The listed OIDs below `begin`.
    below: u32,
    /// The included OIDs: the listed OIDs of `begin..end`.
    count: u32,
}

/// Words of one entry of `OidFilter::supers`.
const SUPER_WORDS: usize = 8;

/// The OIDs of NCBI's loop from OID 0, one OID up at each step, while `CheckOrFindOID` holds for
/// the OID (`OidFilter::walk`): `first`, then the listed OIDs of `from..to`, then `last`,
/// ascending.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct OidWalk {
    pub first: Option<u32>,
    pub from: u32,
    pub to: u32,
    pub last: Option<u32>,
}

impl OidFilter {
    /// The list of the OIDs whose bits are set in `words` (bits at `num_oids` and above
    /// are cleared), with the iteration range of every OID.
    pub(super) fn new(mut words: Vec<u64>, num_oids: u32) -> Self {
        words.resize((num_oids as usize).div_ceil(64), 0);
        if !num_oids.is_multiple_of(64) {
            if let Some(last) = words.last_mut() {
                *last &= (1u64 << (num_oids % 64)) - 1;
            }
        }
        let mut supers = Vec::with_capacity(words.len().div_ceil(SUPER_WORDS));
        let mut count = 0u32;
        for (at, word) in words.iter().enumerate() {
            if at % SUPER_WORDS == 0 {
                supers.push(count);
            }
            count += word.count_ones();
        }
        Self {
            words,
            supers,
            num_oids,
            total: count,
            begin: 0,
            end: num_oids,
            below: 0,
            count,
        }
    }

    /// The same list with the iteration range `(begin, end)` (as `SetIterationRange` leaves it:
    /// `begin <= end <= num_oids`).
    pub(super) fn with_range(mut self, (begin, end): (u32, u32)) -> Self {
        debug_assert!(begin <= end && end <= self.num_oids, "{begin} {end}");
        self.begin = begin;
        self.end = end;
        self.below = self.listed_below(begin);
        self.count = self.listed_below(end) - self.below;
        self
    }

    /// The iteration range `(begin, end)` (of `SeqDb::included_length`).
    #[cfg(any(test, feature = "parallel"))]
    pub(crate) fn range(&self) -> (u32, u32) {
        (self.begin, self.end)
    }

    /// The number of included OIDs: the ordinal ids are `0..count`.
    pub(crate) fn count(&self) -> u32 {
        self.count
    }

    /// Whether `oid` is included.
    pub(crate) fn contains(&self, oid: u32) -> bool {
        (self.begin..self.end).contains(&oid) && self.listed(oid)
    }

    /// The number of included OIDs below `oid`: the ordinal id of `oid` when it is included.
    pub(crate) fn rank(&self, oid: u32) -> u32 {
        self.listed_below(oid.clamp(self.begin, self.end)) - self.below
    }

    /// The OID of ordinal id `ordinal` (below `count`).
    pub(crate) fn select(&self, ordinal: u32) -> u32 {
        let ordinal = self.below + ordinal;
        let block = self.supers.partition_point(|&rank| rank <= ordinal) - 1;
        let mut rest = ordinal - self.supers[block];
        let first = block * SUPER_WORDS;
        for (at, &word) in self.words[first..].iter().enumerate() {
            let ones = word.count_ones();
            if rest < ones {
                let mut word = word;
                for _ in 0..rest {
                    word &= word - 1;
                }
                return ((first + at) * 64) as u32 + word.trailing_zeros();
            }
            rest -= ones;
        }
        unreachable!("ordinal id {ordinal} of a filter of {} OIDs", self.count)
    }

    /// The first OID from `oid` on that is not included, `limit` when every OID of
    /// `oid..limit` is.
    pub(crate) fn run_end(&self, oid: u32, limit: u32) -> u32 {
        if oid >= limit {
            return limit;
        }
        if !(self.begin..self.end).contains(&oid) {
            return oid;
        }
        let limit = limit.min(self.end);
        let mut at = oid;
        while at < limit {
            let word = (at / 64) as usize;
            let Some(&bits) = self.words.get(word) else {
                return at;
            };
            let excluded = !bits >> (at % 64);
            if excluded != 0 {
                return (at + excluded.trailing_zeros()).min(limit);
            }
            at = (word as u32 + 1) * 64;
        }
        limit
    }

    /// Whether NCBI's OID list holds `oid`, inside the iteration range or not.
    pub(crate) fn listed(&self, oid: u32) -> bool {
        self.words
            .get((oid / 64) as usize)
            .is_some_and(|word| word >> (oid % 64) & 1 == 1)
    }

    /// The number of listed OIDs below `oid`.
    pub(crate) fn listed_below(&self, oid: u32) -> u32 {
        let word = (oid / 64) as usize;
        if word >= self.words.len() {
            return self.total;
        }
        let first = word - word % SUPER_WORDS;
        let mut rank = self.supers[first / SUPER_WORDS]
            + self.words[first..word]
                .iter()
                .map(|w| w.count_ones())
                .sum::<u32>();
        let bits = oid % 64;
        if bits > 0 {
            rank += (self.words[word] & ((1u64 << bits) - 1)).count_ones();
        }
        rank
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.cpp:90-140
    /// The first listed OID from `oid` on (`CSeqDBOIDList::CheckOrFindOID`).
    pub(crate) fn next_listed(&self, oid: u32) -> Option<u32> {
        let mut word = (oid / 64) as usize;
        let mut bits = *self.words.get(word)? & (u64::MAX << (oid % 64));
        while bits == 0 {
            word += 1;
            bits = *self.words.get(word)?;
        }
        Some(word as u32 * 64 + bits.trailing_zeros())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:272-297
    // Before the list is set up the range is the constructor's `[0, num_oids)`; `x_GetOidList`
    // sets `[begin, end)` (`alias::top_iteration_range`). The last check is `>`, not `>=`: the
    // OID at the range's end passes when the list holds it.
    /// `x_CheckOrFindOID` of `next`: whether the listed OID it moves `next` to (the first from
    /// `next` on, `next` first clamped to the range's start) passes NCBI's checks of the range,
    /// where `set_up` tells whether NCBI's list was set up before the call (it is after).
    pub(crate) fn check_or_find(&self, next: &mut u32, set_up: bool) -> bool {
        let (begin, end) = if set_up {
            (self.begin, self.end)
        } else {
            (0, self.num_oids)
        };
        if *next < begin {
            *next = begin;
        }
        if *next >= end {
            return false;
        }
        match self.next_listed(*next) {
            Some(oid) => {
                *next = oid;
                oid <= self.end
            }
            None => false,
        }
    }

    /// The OIDs of the loop from OID 0, one OID up at each step, while `x_CheckOrFindOID` holds for
    /// the OID (`check_or_find`), `set_up` whether NCBI's list was set up before the loop. Each call after
    /// the first sees the range `[begin, end)` and moves to the next listed OID, so the loop visits
    /// the listed OIDs of the range from where it is, and `end` itself when the list holds it and
    /// the OID before it was not visited. Before the first call the list is not set up: that call
    /// sees every OID and passes the first listed OID if it is not past `end`, inside the range or
    /// not.
    pub(crate) fn walk(&self, set_up: bool) -> OidWalk {
        let empty = OidWalk {
            first: None,
            from: 0,
            to: 0,
            last: None,
        };
        let (first, from) = if set_up {
            (None, self.begin)
        } else {
            let mut oid = 0;
            if !self.check_or_find(&mut oid, false) {
                return empty;
            }
            (Some(oid), (oid + 1).max(self.begin))
        };
        if from >= self.end {
            return OidWalk { first, ..empty };
        }
        let last = (self.listed(self.end) && !self.listed(self.end - 1)).then_some(self.end);
        OidWalk {
            first,
            from,
            to: self.end,
            last,
        }
    }
}

impl OidWalk {
    /// The OIDs of the walk of `filter`, ascending.
    pub(crate) fn oids(self, filter: &OidFilter) -> impl Iterator<Item = u32> + '_ {
        let mut next = self.from;
        let middle = std::iter::from_fn(move || {
            let oid = filter.next_listed(next).filter(|&oid| oid < self.to)?;
            next = oid + 1;
            Some(oid)
        });
        self.first.into_iter().chain(middle).chain(self.last)
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1910-1941
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1956-1968
// The sets are the taxids that the database holds of those requested: NCBI swaps the list's
// set for them (seqdblmdbset.cpp:315-333, 335-354).
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbvol.hpp:943-962
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbvol.hpp:997-1010
// A list of seqids holds no gi and no ti (`HasIdFilters` is its `GetNumSis()` above 0). The
// lists of the alias files attached to a volume are tested beside this rule ([`line_mask`]).
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbvol.hpp:918-922
// Every list the search opens the database with is a rule, also one left with no id and no taxid
// that NCBI holds (a version 5 seqid list of gis only, a negative taxid list of none the database
// holds): its volumes still run the id test (`x_HaveIdFilter` asks only for the list object), so
// a line without a Seq-id is erased under it (review 2 of line RG D1c, R2).
/// Which definition lines of an included OID a restricted database keeps.
#[derive(Debug)]
pub(crate) enum DeflineRule {
    /// A positive list: the lines whose `taxid` or `links` hold one of these taxids.
    Taxids(BTreeSet<i32>),
    /// A negative list: the lines with a taxid (`GetTaxIds`) outside these.
    NegativeTaxids(BTreeSet<i32>),
    /// A positive seqid list (sorted, unique, as looked up): the lines with a Seq-id in it
    /// (`CSeqDBGiList::FindId`).
    Seqids(Vec<Vec<u8>>),
    /// A negative seqid list (sorted, unique, as looked up): the lines with a Seq-id that is
    /// not a gi or a ti and not in it (`CSeqDBNegativeList::FindId`).
    NegativeSeqids(Vec<Vec<u8>>),
}

impl DeflineRule {
    /// Whether the rule reads the Seq-ids of a line (a seqid list), not only its taxids.
    pub(crate) fn reads_seqids(&self) -> bool {
        matches!(
            self,
            DeflineRule::Seqids(_) | DeflineRule::NegativeSeqids(_)
        )
    }

    // NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbvol.hpp:946-954
    /// The search list's test of one Seq-id `id` of a line (`x_ListIncludesId`): `None` for a
    /// list without ids (`HasIdFilters` false: every line passes), else whether the list includes
    /// `id` or NCBI's exception of its lookup (`super::id_list::find_id`,
    /// `super::id_list::negative_find_id`).
    fn includes_id(&self, id: &ParsedSeqId) -> Option<Result<bool, NcbiException>> {
        match self {
            DeflineRule::Seqids(list) if !list.is_empty() => {
                Some(super::id_list::find_id(&[], &[], list, id))
            }
            DeflineRule::NegativeSeqids(list) if !list.is_empty() => Some(
                super::id_list::negative_find_id(list, id)
                    .map(|(found, match_type)| !found && match_type),
            ),
            _ => None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2036-2039
    /// The positive user taxid test of a line whose ids passed (`s_IncludeDefline_Taxid`): true
    /// unless the search's list is a positive taxid list whose `GetNumTaxIds` count is positive.
    fn keeps_positive_taxid(&self, defline: &BlastDefline) -> bool {
        match self {
            DeflineRule::Taxids(found) if !found.is_empty() => {
                super::id_list::defline_matches_taxids(defline, found)
            }
            _ => true,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2052-2055
    /// The negative user taxid test of a line (`s_IncludeDefline_NegativeTaxid`), run after the
    /// alias files' `TAXIDLIST` rescue: true unless the search's list is a negative taxid list
    /// whose `GetNumTaxIds` count is positive.
    fn keeps_negative_taxid(&self, defline: &BlastDefline) -> bool {
        match self {
            DeflineRule::NegativeTaxids(found) if !found.is_empty() => {
                defline.tax_ids().iter().any(|taxid| !found.contains(taxid))
            }
            _ => true,
        }
    }

    #[cfg(test)]
    fn keeps(&self, defline: &BlastDefline) -> bool {
        self.keeps_positive_taxid(defline) && self.keeps_negative_taxid(defline)
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2016-2055
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/impl/seqdbvol.hpp:943-969
// The id test runs on each line while the volume has a list (`x_HaveIdFilter`: the search's, of
// any kind, or an alias file's): the line's Seq-ids in their order, each asked of the search's
// list until one passes it (`have_user`) and of the alias files' lists attached to the volume
// until one holds one (`have_vol`), the two tests of one id one after the other, until both
// passed; a line with no Seq-id passes neither. The exceptions of the lookups come in that order.
// The taxids of a positive user list are tested next, on the lines the ids keep. Then, for any
// line not yet kept, each alias file's `TAXIDLIST` attached to the volume can rescue the line by
// matching its `taxid` or `links` (`s_IncludeDefline_Taxid`, seqdbvol.cpp:2041-2049), even when
// the user list rejected it or the line holds no Seq-id. Finally, a negative user taxid list tests
// the kept lines (`s_IncludeDefline_NegativeTaxid`, seqdbvol.cpp:2052-2055).
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2059-2066
// NCBI reads no Seq-id and no title for an OID whose every line is erased (`GetBioseq`
// returns null, seqdbvol.cpp:1370-1376): its reports leave the sequence's hits out
// without a word, in every format (dropin-survey-20261029/DI8 probe R01-R07, R14, R21:
// a seqid list naming a PDB molecule without its chain, a local id in other letters, or
// a Swiss-Prot name on a version 4 database selects the OID through the lookup table
// but matches none of its lines). Under a taxid list it happens to an OID whose lines lack
// Seq-ids or whose taxids in the header are not those of the lookup files (makeblastdb writes
// both from the same lines and gives each line a taxid, 0 where none is known; review 2 of line
// RG D1c, R2): the hits leave the reports the same way (`SeqDb::drops_lineless_hits`).
/// Which of the definition lines of `set` NCBI's filtered header read keeps, in order, under
/// the search's list rule `rule` and the alias files' lists `lists` attached to the volume
/// (none: empty): `None` when a seqid rule passes no line and no alias taxid list rescues one;
/// NCBI's exception of a lookup.
pub(crate) fn line_mask(
    rule: Option<&DeflineRule>,
    lists: &[Arc<super::id_list::AliasIdList>],
    set: &DeflineSet,
) -> Result<Option<Vec<bool>>, NcbiException> {
    let mut user_passed_any = false;
    let has_vol_taxids = lists.iter().any(|list| list.has_taxids());
    let mut mask = Vec::with_capacity(set.deflines.len());
    for defline in &set.deflines {
        // Without a list of its kind a test passes; with one, a line passes it by a Seq-id.
        let (mut have_user, mut have_vol) = (rule.is_none(), lists.is_empty());
        for id in &defline.ids {
            if !have_user {
                have_user = match rule.and_then(|rule| rule.includes_id(id)) {
                    Some(includes) => includes?,
                    None => true,
                };
            }
            if !have_vol {
                have_vol = super::id_list::lists_hold(lists, id)?;
            }
            if have_user && have_vol {
                break;
            }
        }
        user_passed_any |= have_user;
        let mut keep = have_user && have_vol;
        if keep {
            keep = rule.is_none_or(|rule| rule.keeps_positive_taxid(defline));
        }
        if !keep && has_vol_taxids {
            keep = lists
                .iter()
                .any(|list| list.has_taxids() && list.keeps_taxid(defline));
        }
        if keep {
            keep = rule.is_none_or(|rule| rule.keeps_negative_taxid(defline));
        }
        mask.push(keep);
    }
    if rule.is_some_and(DeflineRule::reads_seqids) && !user_passed_any && !mask.contains(&true) {
        return Ok(None);
    }
    Ok(Some(mask))
}

/// LOSAT's refusal of an included OID whose every definition line a seqid list removes.
pub(crate) fn no_line_kept(oid: u32) -> String {
    format!(
        "the seqid list keeps sequence {oid} of the volume but none of its definition lines, \
         and NCBI leaves its hits out of the report without notice"
    )
}

/// The definition lines of one sequence in their stored order, as a reader decodes them
/// (`Volume::deflines`, `defline_taxids`, `defline_set`), shortened by a `DeflineRule`.
pub(crate) trait DeflineList {
    /// Keeps the lines `i` with `keep[i]`; a list of another length is left as it is.
    fn retain_lines(&mut self, keep: &[bool]);
}

impl<T> DeflineList for Vec<T> {
    fn retain_lines(&mut self, keep: &[bool]) {
        if keep.len() == self.len() {
            let mut at = 0;
            self.retain(|_| {
                at += 1;
                keep[at - 1]
            });
        }
    }
}

impl DeflineList for DeflineSet {
    fn retain_lines(&mut self, keep: &[bool]) {
        self.deflines.retain_lines(keep);
    }
}

/// The volumes of a database that share one LMDB file (`CSeqDBLMDBEntry`).
pub(super) struct LmdbEntry {
    /// The LMDB file (`.pdb`, `.ndb`).
    lmdb: PathBuf,
    /// The OID of the entry's first volume (`m_OIDStart`).
    oid_start: u32,
    /// The OIDs of the LMDB file's volumes, in its order: the end of each volume's OIDs
    /// (`max_oid`) and their count when the database does not hold the volume
    /// (`skipped_oids`, 0 for a held volume).
    volumes: Vec<(u32, u32)>,
    /// Whether the database holds only some of the LMDB file's volumes (`m_isPartial`).
    partial: bool,
    /// The number of OIDs of all the LMDB file's volumes (`m_NumOids`).
    num_oids: u32,
}

impl LmdbEntry {
    /// The file beside the LMDB file with the last two letters of its extension replaced
    /// (`GetFileNameFromExistingLMDBFile`, seqdb_lmdb.cpp:763-790): `tf`, `to`, `ot`.
    /// The name's bytes are the volume's header's (line LIB-B).
    fn lookup_file(&self, letters: &str) -> PathBuf {
        let mut name = crate::utils::os_bytes::path_bytes(&self.lmdb).into_owned();
        name.truncate(name.len().saturating_sub(2));
        name.extend_from_slice(letters.as_bytes());
        crate::utils::os_bytes::path_from_bytes(name)
    }

    /// The name NCBI's messages give the entry: the LMDB file's name without directory and
    /// extension (`GetDBFileNames`, `CSeqDB_Path::FindBaseName`), its bytes.
    fn base_name(path: &Path) -> Vec<u8> {
        path.file_stem()
            .map(|stem| crate::utils::os_bytes::os_str_bytes(stem).into_owned())
            .unwrap_or_default()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:115-143
    /// The database OIDs of the LMDB file's OIDs `oids`: the OIDs of volumes that the
    /// database does not hold are dropped.
    fn database_oids(&self, oids: Vec<u32>) -> Vec<u32> {
        if !self.partial {
            return oids.into_iter().map(|oid| oid + self.oid_start).collect();
        }
        let mut out = Vec::with_capacity(oids.len());
        for oid in oids {
            let mut skipped = 0;
            for &(max_oid, skipped_oids) in &self.volumes {
                if oid >= max_oid {
                    skipped += skipped_oids;
                } else {
                    if skipped_oids == 0 {
                        out.push(oid + self.oid_start - skipped);
                    }
                    break;
                }
            }
        }
        out
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:165-175
    /// The database OIDs of `tax_ids` in this entry (`negative`: the OIDs the list
    /// excludes) and the taxids the entry holds of `tax_ids`.
    fn tax_oids(
        &self,
        negative: bool,
        tax_ids: &BTreeSet<i32>,
    ) -> Result<tax_files::TaxOids, SeqDbError> {
        let ntf_path = self.lookup_file("tf");
        if !ntf_path.is_file() {
            return Err(missing_file(&ntf_path));
        }
        let ntf_file = DbFile::open(&ntf_path).map_err(lmdb_file_error)?;
        let ntf = Env::open(&ntf_file).map_err(|message| lmdb_open_error(&ntf_file, message))?;
        let nto_path = self.lookup_file("to");
        let nto = DbFile::open(&nto_path).map_err(lookup_file_error)?;
        let tax_error = |error: tax_files::TaxError, path: &Path| match error {
            tax_files::TaxError::Data(message) => SeqDbError::Format {
                path: path.to_path_buf(),
                message,
            },
            tax_files::TaxError::Unreadable(error) => error,
            other => SeqDbError::NotFound(other.message(&Self::base_name(&ntf_path))),
        };
        let mut found = if negative {
            let not_path = self.lookup_file("ot");
            let not = DbFile::open(&not_path).map_err(lookup_file_error)?;
            tax_files::negative_tax_ids_to_oids(&ntf, &nto, &not, self.num_oids, tax_ids)
                .map_err(|error| tax_error(error, &not_path))?
        } else {
            tax_files::get_oids_for_tax_ids(&ntf, &nto, self.num_oids, tax_ids)
                .map_err(|error| tax_error(error, &nto_path))?
        };
        found.oids = self.database_oids(std::mem::take(&mut found.oids));
        Ok(found)
    }

    /// Every taxid key of the entry's `.ntf` (`CSeqDBLMDBEntry::GetDBTaxIds`, the keys of
    /// `CSeqDBLMDB::GetDBTaxIds`, cited at `tax_files::get_db_tax_ids`), with the errors of
    /// `tax_oids`.
    fn db_tax_ids(&self) -> Result<Vec<i32>, SeqDbError> {
        let ntf_path = self.lookup_file("tf");
        if !ntf_path.is_file() {
            return Err(missing_file(&ntf_path));
        }
        let ntf_file = DbFile::open(&ntf_path).map_err(lmdb_file_error)?;
        let ntf = Env::open(&ntf_file).map_err(|message| lmdb_open_error(&ntf_file, message))?;
        tax_files::get_db_tax_ids(&ntf).map_err(|error| match error {
            tax_files::TaxError::Data(message) => SeqDbError::Format {
                path: ntf_path.clone(),
                message,
            },
            tax_files::TaxError::Unreadable(error) => error,
            other => SeqDbError::NotFound(other.message(&Self::base_name(&ntf_path))),
        })
    }

    /// The end of the database's OIDs that the entry holds (`m_OIDEnd`, seqdblmdbset.cpp:
    /// 55-73: the OIDs of the held volumes after `m_OIDStart`).
    fn oid_end(&self) -> u32 {
        let skipped: u32 = self.volumes.iter().map(|&(_, skipped)| skipped).sum();
        self.oid_start + (self.num_oids - skipped)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:182-201
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:737-746
    /// The taxids that the entry's `.pot`/`.not` file gives the OIDs `oids` (the entry's own
    /// numbering, ascending), added to `tax_ids`. An OID past the file's is LOSAT's refusal
    /// (NCBI reads past its index there).
    fn oid_tax_ids(&self, oids: &[u32], tax_ids: &mut BTreeSet<i32>) -> Result<(), SeqDbError> {
        let path = self.lookup_file("ot");
        let file = DbFile::open(&path).map_err(lookup_file_error)?;
        let format = |message: String| SeqDbError::Format {
            path: path.clone(),
            message,
        };
        let lookup = tax_files::OidTaxIds::open(&file).map_err(|error| error.into_seqdb(format))?;
        let (mut skipped, mut next) = (0u32, 0usize);
        for &oid in oids {
            let file_oid = if self.partial {
                while next < self.volumes.len()
                    && (self.volumes[next].1 > 0 || oid + skipped >= self.volumes[next].0)
                {
                    skipped += self.volumes[next].1;
                    next += 1;
                }
                oid + skipped
            } else {
                oid
            };
            tax_ids.extend(
                lookup
                    .taxids(file_oid)
                    .map_err(|error| error.into_seqdb(format))?,
            );
        }
        Ok(())
    }

    /// The LMDB file opened for a lookup, handed to `f`.
    fn with_env<T>(
        &self,
        f: impl FnOnce(&Env<DbFile>) -> Result<T, SeqDbError>,
    ) -> Result<T, SeqDbError> {
        let file = self.open_file()?;
        f(&self.open_env(&file)?)
    }

    /// The LMDB file, opened for lookups (`with_env`, `lookup::SeqidLookup`).
    pub(super) fn open_file(&self) -> Result<DbFile, SeqDbError> {
        DbFile::open(&self.lmdb).map_err(lmdb_file_error)
    }

    /// The LMDB environment over `file`, this entry's file (`open_file`).
    pub(super) fn open_env<'f>(&self, file: &'f DbFile) -> Result<Env<'f, DbFile>, SeqDbError> {
        Env::open(file).map_err(|message| lmdb_open_error(file, message))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:81-113
    /// The database OID of the LMDB file's OID `oid` of a seqid lookup: `None` when its
    /// volume is one the database does not hold.
    pub(super) fn database_oid(&self, oid: u32) -> Option<u32> {
        if !self.partial {
            return oid.checked_add(self.oid_start);
        }
        let mut skipped = 0;
        for &(max_oid, skipped_oids) in &self.volumes {
            if oid >= max_oid {
                skipped += skipped_oids;
            } else {
                return (skipped_oids == 0)
                    .then(|| oid.checked_add(self.oid_start)?.checked_sub(skipped))
                    .flatten();
            }
        }
        None
    }

    /// `CSeqDBLMDB::GetOids` of one id: the first OID of the key (`tax_files::acc_first_oid`),
    /// with NCBI's messages for a database without the table (`DB contains no accession
    /// info.`, the text of NCBI's LMDB manager for a `.pdb`/`.ndb` without `acc2oid`: a
    /// database made without `-parse_seqids`, probe S48) and for a damaged one.
    fn first_oid(&self, env: &Env<DbFile>, acc: &[u8]) -> Result<Option<u32>, SeqDbError> {
        tax_files::acc_first_oid(env, acc).map_err(|error| self.lookup_error(error))
    }

    /// The error of a lookup in the `acc2oid` table: NCBI's exception for a database without
    /// the table and for an empty key (`MDB_BAD_VALSIZE`, probe H01), LOSAT's for a damaged
    /// file.
    fn lookup_error(&self, error: ReadError) -> SeqDbError {
        let message = match error {
            ReadError::Format(message) => message,
            ReadError::Unreadable(error) => return error,
        };
        let name = Self::base_name(&self.lmdb);
        if message == tax_files::NO_ACCESSION_INFO {
            SeqDbError::NotFound(message.into_bytes())
        } else if message.starts_with("MDB_BAD_VALSIZE") {
            SeqDbError::NotFound(tax_files::accessions_lookup_error(&name))
        } else {
            SeqDbError::Format {
                path: self.lmdb.clone(),
                message: format!(
                    "{} ({message})",
                    String::from_utf8_lossy(&tax_files::accessions_lookup_error(&name))
                ),
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:146-151
    // `GetOid` fails on an accession of more than about 500 OIDs as `GetOids` does (probe
    // H107): LOSAT reads every duplicate (`PD-LOSAT-NCBI-DEFECTS` exception 6).
    /// All the OIDs of the key `acc` in the LMDB file, in LMDB's order.
    pub(super) fn oids(&self, acc: &[u8]) -> Result<Vec<u32>, SeqDbError> {
        self.with_env(|env| self.oids_in(env, acc))
    }

    /// `oids` in the entry's LMDB file opened as `env` (`open_env`).
    pub(super) fn oids_in(&self, env: &Env<DbFile>, acc: &[u8]) -> Result<Vec<u32>, SeqDbError> {
        tax_files::acc_oids(env, acc).map_err(|error| self.lookup_error(error))
    }

    /// The `acc2oid` database of the LMDB file opened as `env` (`tax_files::acc2oid`), for
    /// `oids_of` and `case_variant_oids_in`.
    pub(super) fn acc2oid_in(&self, env: &Env<DbFile>) -> Result<Db, SeqDbError> {
        tax_files::acc2oid(env).map_err(|error| self.lookup_error(error))
    }

    /// `oids_in` in the `acc2oid` database `db` of `env` (`acc2oid_in`).
    pub(super) fn oids_of(
        &self,
        env: &Env<DbFile>,
        db: &Db,
        acc: &[u8],
    ) -> Result<Vec<u32>, SeqDbError> {
        tax_files::acc_oids_of(env, db, acc).map_err(|error| self.lookup_error(error))
    }

    /// The OIDs of the other keys that differ from `acc` only in the case of letters at or after
    /// byte `from` (`tax_files::acc_case_variant_oids`, the prefixes found kept in `prefixes`), in
    /// the `acc2oid` database `db` of the LMDB file opened as `env`.
    pub(super) fn case_variant_oids_in(
        &self,
        env: &Env<DbFile>,
        db: &Db,
        acc: &[u8],
        from: usize,
        prefixes: &mut rustc_hash::FxHashMap<Vec<u8>, bool>,
    ) -> Result<Vec<u32>, SeqDbError> {
        tax_files::acc_case_variant_oids(env, db, acc, from, prefixes)
            .map_err(|error| self.lookup_error(error))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:153-157
    // `GetOids` (seqdb_lmdb.cpp:360-400) reads the first duplicate of each key with
    // `MDB_SET` and `MDB_GET_CURRENT`: an accession of more than about 500 OIDs, whose
    // duplicates LMDB keeps in a sub-database, makes NCBI's read fail (`Accessions to Oids
    // lookup error`); LOSAT reads its first duplicate like any other key (DI-8 DECISIONS 3,
    // PD-LOSAT-NCBI-DEFECTS).
    /// The database OID of each of `accs` in this entry (`None`: not found).
    fn accessions_to_oids(&self, accs: &[Vec<u8>]) -> Result<Vec<Option<u32>>, SeqDbError> {
        self.with_env(|env| {
            accs.iter()
                .map(|acc| {
                    Ok(self
                        .first_oid(env, acc)?
                        .and_then(|oid| self.database_oid(oid)))
                })
                .collect()
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:402-413
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:514-556
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:159-163
    // An excluded OID of a volume the database does not hold becomes `kSeqDBEntryNotFound`,
    // which NCBI then clears as a bit (`ClearBit` of -1, undefined): NCBI crashes (SIGSEGV,
    // exit 139, review of DI-8f: `-db mv/m.00`, or `-db "d1/m d2/m"` whose volumes
    // interleave). LOSAT refuses the list (PD-LOSAT-NCBI-DEFECTS, rule 3).
    /// The database OIDs that the negative list `ids` excludes in this entry: those whose
    /// every id of the `.nos`/`.pos` file is covered by the listed ids that have them as
    /// first OID.
    fn negative_seqids_to_oids(&self, ids: &[Vec<u8>]) -> Result<Vec<u32>, SeqDbError> {
        let oids = self.with_env(|env| {
            ids.iter()
                .map(|id| self.first_oid(env, id))
                .collect::<Result<Vec<_>, _>>()
        })?;
        let mut pairs: Vec<(u32, &[u8])> = oids
            .iter()
            .zip(ids)
            .filter_map(|(oid, id)| oid.map(|oid| (oid, id.as_slice())))
            .collect();
        if pairs.is_empty() {
            return Ok(Vec::new());
        }
        pairs.sort();
        let nos_path = self.lookup_file("os");
        let nos_file = DbFile::open(&nos_path).map_err(lookup_file_error)?;
        let lookup = tax_files::OidSeqIds::open(&nos_file).map_err(|error| {
            error.into_seqdb(|message| SeqDbError::Format {
                path: nos_path.clone(),
                message,
            })
        })?;
        let mut excluded = Vec::new();
        let mut i = 0;
        while i < pairs.len() {
            let current_oid = pairs[i].0;
            let file_ids = lookup.seqids(current_oid).map_err(|error| {
                error.into_seqdb(|message| SeqDbError::Format {
                    path: nos_path.clone(),
                    message,
                })
            })?;
            let mut input_ids = Vec::new();
            while i < pairs.len() && pairs[i].0 == current_oid {
                input_ids.push(pairs[i].1);
                i += 1;
            }
            if compare_id_list(&file_ids, &input_ids).map_err(|message| SeqDbError::Format {
                path: nos_path.clone(),
                message,
            })? {
                excluded.push(current_oid);
            }
        }
        excluded
            .into_iter()
            .map(|oid| {
                self.database_oid(oid)
                    .ok_or_else(|| SeqDbError::Unsupported {
                        path: self.lmdb.clone(),
                        message: format!(
                            "the negative seqid list excludes sequence {oid} of the LMDB file, in \
                         a volume that the database does not hold; NCBI BLAST+ then clears \
                         bit -1 of its OID mask and crashes"
                        ),
                    })
            })
            .collect()
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:461-513
// `find_first_of` of `tmp_pdb` being 0 asks whether the first character of the id is one of the
// characters of the molecule name, not whether the id starts with the name: kept as it is.
// A PDB id of the list with a molecule name that is not the file id's falls through to
// `break` (the `else` is not taken), as in NCBI.
/// Whether the sorted ids of the file (`file_ids`, an OID's `.nos` entries) are all covered
/// by the sorted listed ids `input_ids` that resolved to the OID. `Err`: an id that does not
/// parse (NCBI throws `CSeqIdException`).
fn compare_id_list(file_ids: &[Vec<u8>], input_ids: &[&[u8]]) -> Result<bool, String> {
    let parse = |text: &[u8]| {
        seqidlist::parse_list_id(text).map_err(|error| {
            format!(
                "{}: {}",
                String::from_utf8_lossy(text),
                String::from_utf8_lossy(&error.what()).trim_end()
            )
        })
    };
    let starts_with_one_of = |s: &[u8], set: &[u8]| s.first().is_some_and(|c| set.contains(c));
    let (mut f, mut i) = (0usize, 0usize);
    while f < file_ids.len() && i < input_ids.len() {
        if input_ids[i] == file_ids[f].as_slice() {
            i += 1;
            f += 1;
            continue;
        }
        let seq_id = parse(input_ids[i])?;
        if let ParsedSeqId::Pdb(pdb) = &seq_id {
            if pdb.chain_id.is_none() {
                let file_seqid = parse(&file_ids[f])?;
                if matches!(&file_seqid, ParsedSeqId::Pdb(file_pdb) if file_pdb.mol == input_ids[i])
                {
                    f += 1;
                    let tmp_pdb = input_ids[i];
                    while f < file_ids.len() && starts_with_one_of(&file_ids[f], tmp_pdb) {
                        f += 1;
                    }
                    while i < input_ids.len() && starts_with_one_of(input_ids[i], tmp_pdb) {
                        i += 1;
                    }
                    continue;
                }
            }
        } else {
            let file_seq_id = parse(&file_ids[f])?;
            if blast_seq_id_string_plain(&file_seq_id, false) == input_ids[i] {
                i += 1;
                if i < input_ids.len()
                    && blast_seq_id_string_plain(&file_seq_id, true) == input_ids[i]
                {
                    i += 1;
                }
                f += 1;
                continue;
            }
        }
        break;
    }
    Ok(f == file_ids.len())
}

/// `CSeq_id::GetSeqIdString` with `with_version`.
fn blast_seq_id_string_plain(id: &ParsedSeqId, with_version: bool) -> Vec<u8> {
    let mut out = Vec::new();
    label::content(id, with_version, &mut out);
    out
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:85-89
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:288-300
/// One OID for each of `accs`, the last LMDB file that has the key deciding. On an error
/// `oids` keeps the answers of the files before the one that failed.
pub(super) fn entries_accessions_to_oids(
    entries: &[LmdbEntry],
    accs: &[Vec<u8>],
    oids: &mut Vec<Option<u32>>,
) -> Result<(), SeqDbError> {
    oids.clear();
    oids.resize(accs.len(), None);
    for entry in entries {
        for (slot, oid) in oids.iter_mut().zip(entry.accessions_to_oids(accs)?) {
            if oid.is_some() {
                *slot = oid;
            }
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:91-96
// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:1461-1462
// (LMDB's notice and the OpenLDAP Public License 2.8: THIRD_PARTY_NOTICES, section 2 (LMDB).)
/// The error of an LMDB file (`.pdb`/`.ndb`, `.ptf`/`.ntf`) that opens but is not read: NCBI's
/// text for an empty file (probe SR4; an empty `.ptf` of a taxid list: line WQA D18 review M4,
/// NCBI BLAST+ 2.17.0 on 2026-11-21), a file of another format (SR3) and another LMDB version;
/// LOSAT's for the rest.
fn lmdb_open_error(file: &DbFile, error: ReadError) -> SeqDbError {
    let message = match error {
        ReadError::Format(message) => message,
        // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3683-3712
        // `mdb_env_open` reads the two meta pages with `pread` before it maps the file, so a read
        // that fails there is NCBI's `mdb_env_open` error with the system's text, as an open or a
        // first byte that fails (line WQA D18c, 2026-11-21).
        ReadError::Unreadable(SeqDbError::Read { path, source }) => {
            return match crate::utils::os_error_text::ncbi_strerror(&source) {
                Some(strerror) => mdb_env_open_error(&strerror),
                None => SeqDbError::Read { path, source },
            };
        }
        ReadError::Unreadable(error) => return lmdb_file_error(error),
    };
    let what = if file.is_empty() {
        "Bad file descriptor"
    } else if message.contains("MDB_VERSION_MISMATCH") {
        "MDB_VERSION_MISMATCH: Database environment version mismatch"
    } else if message.contains("MDB_INVALID") || message.contains("too small for a meta page") {
        "MDB_INVALID: File is not an LMDB file"
    } else {
        return SeqDbError::Format {
            path: file.path().to_path_buf(),
            message,
        };
    };
    mdb_env_open_error(what)
}

/// NCBI's error of an LMDB file that `mdb_env_open` refuses, with LMDB's or the system's text.
fn mdb_env_open_error(what: &str) -> SeqDbError {
    SeqDbError::NotFound(format!("LMDB runtime error: mdb_env_open: {what}").into())
}

/// The error of an LMDB file (`.pdb`/`.ndb`, `.ptf`/`.ntf`) that does not open or whose first
/// byte cannot be read (`DbFile::open`): `mdb_env_open`'s, with the system's text
/// (`utils::os_error_text::ncbi_strerror`, so LOSAT Web gives the native CLI's). NCBI BLAST+
/// 2.17.0 gives it for every code, EMFILE and ENFILE too (a file without read permission:
/// `Permission denied`, checked with NCBI on 2026-10-31; `open()`, `read()` or `mmap()` of the
/// file made to fail: line WQA D18 probes 4-7, 2026-11-21).
fn lmdb_file_error(error: SeqDbError) -> SeqDbError {
    if let SeqDbError::Open { source, .. } | SeqDbError::Map { source, .. } = &error {
        if let Some(strerror) = crate::utils::os_error_text::ncbi_strerror(source) {
            return mdb_env_open_error(&strerror);
        }
    }
    error
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:537-538
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:631-632
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:714-715
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:739-740
/// The error of a lookup file (`.pos`/`.nos`, `.pot`/`.not`, `.pto`/`.nto`) that does not open
/// or whose first byte cannot be read (`DbFile::open`): NCBI maps these with `CMemoryFile`
/// itself, outside the reader's atlas, so `CMemoryFile`'s exception reaches the program (`Error:
/// ` and the report, exit status 255; NCBI BLAST+ 2.17.0 with the file without read permission,
/// or its `open()` or `mmap()` made to fail, EMFILE included: line WQA D18 probes 5-7,
/// 2026-11-21).
fn lookup_file_error(error: SeqDbError) -> SeqDbError {
    match error.memory_file_exception() {
        Some(exception) => SeqDbError::Ncbi(exception),
        None => error,
    }
}

fn missing_file(path: &Path) -> SeqDbError {
    SeqDbError::NotFound(
        [
            b"File ".as_slice(),
            &crate::utils::os_bytes::path_bytes(path)[..],
            b" not found. If you renamed any BLAST database files, please use original file \
         names, and makeblastdb to rename the database. If you deleted any BLAST database files, \
         you need to recreate the database.",
        ]
        .concat(),
    )
}

/// A lookup file names an OID the database does not have.
fn bad_lookup_oid(lmdb: &Path, oid: u32, num_oids: u32) -> SeqDbError {
    SeqDbError::Format {
        path: lmdb.to_path_buf(),
        message: format!("the seqid lookup gives OID {oid} of a database of {num_oids} sequences"),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1300-1303
/// NCBI's text for a taxid lookup (`TaxIdsToOids`, `GetDBTaxIds`) on a version 4 database
/// opened without a list (the list's own text is `V4_TAXID_ERROR`).
const V4_TAXID_LIST_ERROR: &str = "Taxonomy list is not supported in v4 BLAST db";

/// The OIDs of `tax_ids` in the LMDB entries `entries` (`negative`: the OIDs the list
/// excludes) entry after entry, each entry's ascending, and the taxids of `tax_ids` the
/// entries hold (`CSeqDBLMDBSet::TaxIdsToOids`; the found taxids are their union, see
/// `restrict_taxids`).
pub(super) fn entries_tax_oids(
    entries: &[LmdbEntry],
    negative: bool,
    tax_ids: &BTreeSet<i32>,
) -> Result<(Vec<u32>, BTreeSet<i32>), SeqDbError> {
    let mut oids = Vec::new();
    let mut found = BTreeSet::new();
    for entry in entries {
        let part = entry.tax_oids(negative, tax_ids)?;
        oids.extend(part.oids);
        found.extend(part.tax_ids_found);
    }
    Ok((oids, found))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:369-397
/// `CSeqDBLMDBSet::GetTaxIdsForOids`: the taxids of the database's OIDs `oids` (ascending),
/// added to `tax_ids`. NCBI moves to the next entry once per OID, so an OID two entries on is
/// read in the numbering of the next entry, as here; past the last entry NCBI reads past its
/// vector, which LOSAT refuses.
fn entries_oid_tax_ids(
    entries: &[LmdbEntry],
    oids: &[u32],
    tax_ids: &mut BTreeSet<i32>,
) -> Result<(), SeqDbError> {
    let Some(first) = entries.first() else {
        return Ok(());
    };
    if entries.len() == 1 {
        return first.oid_tax_ids(oids, tax_ids);
    }
    let past_entries = || SeqDbError::Unsupported {
        path: entries[entries.len() - 1].lmdb.clone(),
        message: "NCBI BLAST+ reads the taxids of an OID past its LMDB entries".to_string(),
    };
    let (mut part, mut at) = (Vec::new(), 0usize);
    for &oid in oids {
        if oid >= entries[at].oid_end() {
            if !part.is_empty() {
                entries[at].oid_tax_ids(&part, tax_ids)?;
                part.clear();
            }
            at += 1;
        }
        let entry = entries.get(at).ok_or_else(past_entries)?;
        part.push(oid.checked_sub(entry.oid_start).ok_or_else(past_entries)?);
    }
    if !part.is_empty() {
        entries[at].oid_tax_ids(&part, tax_ids)?;
    }
    Ok(())
}

/// NCBI's text when a taxid list selects no OID of `entries` (`CSeqDBLMDBSet::TaxIdsToOids`,
/// cited at `SeqDb::restrict_taxids`).
pub(super) fn tax_ids_not_found(entries: &[LmdbEntry]) -> Vec<u8> {
    let names: Vec<Vec<u8>> = entries
        .iter()
        .map(|entry| LmdbEntry::base_name(&entry.lmdb))
        .collect();
    [
        b"Taxonomy ID(s) not found in the ".as_slice(),
        &names.join(b" ".as_slice())[..],
        b" database",
        if entries.len() > 1 {
            b"s.".as_slice()
        } else {
            b".".as_slice()
        },
    ]
    .concat()
}

/// The taxid lookup files of a version 5 database opened without a list, for the lookups
/// of `losat blastdbcmd` (`-taxids` without `-target_only`, `-tax_info`): NCBI's
/// `CSeqDBImpl::TaxIdsToOids` and `GetDBTaxIds` over the database's `CSeqDBLMDBSet`.
pub struct TaxIdLookup<'a> {
    db: &'a SeqDb,
    entries: Vec<LmdbEntry>,
}

impl SeqDb {
    /// The taxid lookup files of this database ([`TaxIdLookup`]); the LMDB files are read
    /// here. Errors: NCBI's for a version 4 database (`Taxonomy list is not supported in v4
    /// BLAST db`) and for the LMDB files ([`SeqDb::restrict_taxids`]).
    pub fn tax_id_lookup(&self) -> Result<TaxIdLookup<'_>, SeqDbError> {
        if !self.is_version5() {
            return Err(SeqDbError::NotFound(V4_TAXID_LIST_ERROR.to_string().into()));
        }
        Ok(TaxIdLookup {
            db: self,
            entries: self.lmdb_entries()?,
        })
    }
}

impl TaxIdLookup<'_> {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1285-1299
    /// `TaxIdsToOids` of `tax_ids`: the OIDs of the sequences that carry one of `tax_ids`, in the
    /// order of the lookup files, that pass NCBI's check (`SeqDb::check_oid`; the first one
    /// checked can set NCBI's OID list up). Errors: NCBI's when no OID carries one (`Taxonomy
    /// ID(s) not found in the <names> database[s].`), and for the lookup files
    /// ([`SeqDb::restrict_taxids`]).
    pub fn tax_ids_to_oids(&self, tax_ids: &BTreeSet<i32>) -> Result<Vec<u32>, SeqDbError> {
        let (oids, _) = entries_tax_oids(&self.entries, false, tax_ids)?;
        if oids.is_empty() {
            return Err(SeqDbError::NotFound(tax_ids_not_found(&self.entries)));
        }
        let mut kept = Vec::with_capacity(oids.len());
        for oid in oids {
            if self.db.check_oid(oid)? {
                kept.push(oid);
            }
        }
        Ok(kept)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1307-1325
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:356-366
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1317-1322
    // NCBI's OID list is that of the alias files' ranges here (`FIRST_OID`, `LAST_OID`: the
    // database's filter; a `losat blastdbcmd` database has no other list): the taxids of the
    // OIDs of the loop, which sets the list up first (`SeqDb::oid_walk`: the listed OIDs of the
    // iteration range, and its end when the list holds it). Without it they are every key of
    // every entry's `.ntf`, also those of the volumes an alias of some volumes leaves out.
    /// `GetDBTaxIds`: the taxids of the database's taxid lookup files (of the OIDs NCBI's loop
    /// visits where the alias files have ranges), ascending.
    pub fn db_tax_ids(&self) -> Result<BTreeSet<i32>, SeqDbError> {
        let mut tax_ids = BTreeSet::new();
        self.db.build_oid_list()?;
        if self.db.filter.is_some() {
            let oids: Vec<u32> = self.db.oid_walk()?.collect();
            entries_oid_tax_ids(&self.entries, &oids, &mut tax_ids)?;
            return Ok(tax_ids);
        }
        for entry in &self.entries {
            tax_ids.extend(entry.db_tax_ids()?);
        }
        Ok(tax_ids)
    }
}

/// The OIDs of the sequences `first..end` of a volume read per window of the index.
const SCAN_WINDOW: u32 = 1 << 16;

/// The OID of a database sequence: `n` of the `gnl|BL_ORD_ID|n` (`BL_ORD_ID:n` in SAM, the tag
/// of a Seq-align's general id) that NCBI prints for a sequence without parsed Seq-ids. A
/// database restricted by a list or by its alias files' range numbers the sequences that it
/// searches by ordinal ids, their rank in the list; only [`SeqDb::oid`] and [`OidMap::oid`] make
/// an OID of one, so that no report prints a rank as an OID.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Oid(u32);

impl Oid {
    /// The OID as a number.
    pub(crate) fn get(self) -> u32 {
        self.0
    }

    /// `gnl|BL_ORD_ID|n`, how NCBI's reports and the FASTA ids of LOSAT's database records name
    /// the sequence.
    pub(crate) fn fasta_id(self) -> String {
        format!("gnl|BL_ORD_ID|{}", self.0)
    }

    /// `BL_ORD_ID:n`, the content label of the general id `gnl|BL_ORD_ID|n`
    /// (`CSeq_id::GetLabel` with `eContent`, `GetSeqIdString` with true): the accession that the
    /// tabular formats print for it (`sacc`, `saccver`) and the name that SAM gives it.
    pub(crate) fn content_label(self) -> String {
        format!("BL_ORD_ID:{}", self.0)
    }

    /// The OID `oid` of a test's database that is not restricted.
    #[cfg(test)]
    pub(crate) fn of_unrestricted(oid: u32) -> Self {
        Oid(oid)
    }
}

/// The OID, by which LOSAT's refusals name a database sequence ("database sequence n", as
/// `blastdbcmd -outfmt %o` prints it), not by its rank in a list or an alias range
/// that restricts the search: each takes an `Oid`, a [`DbSequence`], or a closure that gives
/// one where it refuses.
impl std::fmt::Display for Oid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

/// The OID of each ordinal id of a database ([`SeqDb::oid`]), for the reports, which hold the
/// ordinal ids of their sequences but not the database.
#[derive(Clone, Copy, Debug)]
pub struct OidMap<'a>(Option<&'a OidFilter>);

impl OidMap<'_> {
    /// The OIDs of a test's database that is not restricted: each ordinal id's own.
    #[cfg(test)]
    pub(crate) fn unrestricted() -> Self {
        OidMap(None)
    }

    /// The OID of ordinal id `ordinal`.
    pub(crate) fn oid(self, ordinal: u32) -> Oid {
        Oid(match self.0 {
            Some(filter) => filter.select(ordinal),
            None => ordinal,
        })
    }

    /// The ordinal ids whose OIDs are of `from..to` (of the database's OIDs; past the database's
    /// sequences there are none, which the reader's blocks end).
    pub(crate) fn ordinals_in(self, from: u32, to: u32) -> std::ops::Range<u32> {
        match self.0 {
            Some(filter) => filter.rank(from)..filter.rank(to),
            None => from..to.max(from),
        }
    }

    /// Whether an ordinal id has an OID of `from..to` (of the database's OIDs).
    pub(crate) fn any_in(self, from: u32, to: u32) -> bool {
        match self.0 {
            Some(filter) => filter.rank(to) > filter.rank(from),
            None => from < to,
        }
    }
}

/// A database sequence that the report of another names (one that the lookup of the other's
/// Seq-align id gives, that also holds it, or that is a copy of its record): its ordinal id, by
/// which LOSAT reads it, and its OID, by which LOSAT's refusals name it (its `Display`). Only
/// [`SeqDb::sequence`] makes one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DbSequence {
    /// The ordinal id.
    ordinal: u32,
    /// The OID ([`SeqDb::oid`] of `ordinal`).
    oid: Oid,
}

impl DbSequence {
    /// The ordinal id, by which LOSAT reads the sequence.
    pub(crate) fn ordinal(self) -> u32 {
        self.ordinal
    }

    /// The sequence `ordinal` of a test's database that is not restricted.
    #[cfg(test)]
    pub(crate) fn of_unrestricted(ordinal: u32) -> Self {
        DbSequence {
            ordinal,
            oid: Oid(ordinal),
        }
    }
}

/// The OID, as LOSAT's refusals name the sequence.
impl std::fmt::Display for DbSequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.oid.fmt(f)
    }
}

impl SeqDb {
    /// The sequence of ordinal id `ordinal` with its OID ([`SeqDb::oid`]).
    pub(crate) fn sequence(&self, ordinal: u32) -> DbSequence {
        DbSequence {
            ordinal,
            oid: self.oid(ordinal),
        }
    }

    /// The OID of ordinal id `ordinal`: `ordinal` itself unless a list restricts the database
    /// ([`SeqDb::restrict`]) or its alias files' range starts past OID 0. NCBI prints it as `n`
    /// of `gnl|BL_ORD_ID|n`.
    pub(crate) fn oid(&self, ordinal: u32) -> Oid {
        self.oid_map().oid(ordinal)
    }

    /// The OIDs of the database's ordinal ids ([`SeqDb::oid`]).
    pub(crate) fn oid_map(&self) -> OidMap<'_> {
        OidMap(self.filter.as_deref())
    }

    /// The ordinal id of the OID `oid`, one the database's list includes (`oid` itself unless a
    /// list or an alias range restricts the database): the inverse of [`SeqDb::oid`].
    pub fn ordinal_of(&self, oid: u32) -> u32 {
        match &self.filter {
            Some(filter) => filter.rank(oid),
            None => oid,
        }
    }

    /// Whether the database is restricted by a list ([`SeqDb::restrict`]); the alias files'
    /// ranges alone do not make it so.
    pub fn is_restricted(&self) -> bool {
        self.list_applied
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:229-262
    // NCBI builds its OID list at its first use: `x_CheckOrFindOID`, and the readers that call
    // `x_GetOidList` first (`AccessionToOids`, `GetDBTaxIds`, `GetBioseq`, `GetHdr`,
    // `GetSeqIDs`, `OidToGi`, `GetNextOIDChunk`), or in the constructor that scans the totals
    // (`apply_filter`, `check_alias_values`). Only `x_CheckOrFindOID` sees the difference: before
    // the build its range is every OID. Where an alias file's OID list file fails and NCBI opened
    // the database without building the list (`read_alias_tree`), the build throws, and
    // `m_OidListSetup` stays false, so every later use throws again.
    /// NCBI's build of its OID list (`x_GetOidList`), the one gate of every reader that uses it:
    /// whether the list was built before this call, which builds it; the error of an alias
    /// file's OID list file where the build fails.
    pub fn build_oid_list(&self) -> Result<bool, SeqDbError> {
        if let Some(failure) = &self.oid_list_error {
            return Err(failure.to_error());
        }
        let set_up = self.oid_list_set_up.load(Ordering::Relaxed);
        if !set_up {
            self.oid_list_set_up.store(true, Ordering::Relaxed);
        }
        Ok(set_up)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:300-313
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:241-251 (at
    // `common::database::PreliminaryReadError`)
    // A search asks for its first chunk of OIDs (`GetNextOIDChunk`, which builds the list) just
    // before it reads their sequences, in its preliminary stage: a build that fails fails the
    // search where the first read of a sequence file fails, as one that cannot be opened does
    // (probes G2c, G8 of line RG D1b: `BLAST engine error: Database memory map file error`, exit
    // status 3; psiblast's full text with exit status 2). The other readers read the sequences
    // that NCBI reads without building the list (blastdbcmd `-ipg` with `%s`, `%h`, probe S3 of
    // review A of line RG D1b-1).
    /// Makes every read of a sequence fail with the failure of NCBI's build of its OID list
    /// ([`SeqDb::build_oid_list`]), as the reads of a search follow its first request of OIDs:
    /// the database's opening for a search (`common::database::open_search_database`).
    pub fn fail_search_reads_with_oid_list(&mut self) {
        if let Some(failure) = &self.oid_list_error {
            for volume in &mut self.volumes {
                volume.fail_sequence_reads(failure.clone());
            }
        }
    }

    /// `x_CheckOrFindOID` of `next` (`OidFilter::check_or_find`; without an OID list every OID of
    /// the volumes passes): whether `next`, moved to the first OID from it on that NCBI's list
    /// holds, passes NCBI's checks of the iteration range. Builds the list
    /// ([`SeqDb::build_oid_list`]).
    pub fn check_or_find_oid(&self, next: &mut u32) -> Result<bool, SeqDbError> {
        let set_up = self.build_oid_list()?;
        Ok(match &self.filter {
            Some(filter) => filter.check_or_find(next, set_up),
            None => *next < self.num_oids(),
        })
    }

    /// Whether the database includes `oid` ([`SeqDb::num_seqs`]: NCBI's list holds it, inside
    /// the iteration range): `x_CheckOrFindOID` of `oid2` together with the test that `oid2` is
    /// `oid` once NCBI's list is set up, as it is for the lookups after a search
    /// ([`SeqDb::check_oid`] at any time). Every NCBI caller of this answer (`SeqidToOids`,
    /// seqdbimpl.cpp:1388) sets the list up first; a lookup that NCBI makes before
    /// (`AccessionsToOids`, `TaxIdsToOids`) takes [`SeqDb::check_oid`], whose first check is the
    /// stale one.
    pub fn includes_oid(&self, oid: u32) -> bool {
        match &self.filter {
            Some(filter) => filter.contains(oid),
            None => oid < self.num_oids(),
        }
    }

    /// `x_CheckOrFindOID` of `oid2` together with the test that `oid2` is `oid`: whether NCBI's
    /// check keeps `oid`. Once the list is set up, the OIDs it keeps are the included ones
    /// ([`SeqDb::includes_oid`]); the first check of a list not set up keeps any listed OID up to
    /// the range's end.
    pub fn check_oid(&self, oid: u32) -> Result<bool, SeqDbError> {
        let mut next = oid;
        Ok(self.check_or_find_oid(&mut next)? && next == oid)
    }

    /// The OIDs of NCBI's loop from OID 0, one OID up at each step, while `CheckOrFindOID` holds
    /// for the OID, ascending (`OidFilter::walk`; without an OID list every OID of the volumes). The loop builds the
    /// list ([`SeqDb::build_oid_list`]).
    pub fn oid_walk(&self) -> Result<impl Iterator<Item = u32> + '_, SeqDbError> {
        let set_up = self.build_oid_list()?;
        let (filter, walk) = match &self.filter {
            Some(filter) => (Some(&**filter), Some(filter.walk(set_up))),
            None => (None, None),
        };
        let all = filter.is_none().then_some(0..self.num_oids());
        let listed = walk
            .zip(filter)
            .into_iter()
            .flat_map(|(walk, filter)| walk.oids(filter));
        Ok(all.into_iter().flatten().chain(listed))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:213-262
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:44-71
    // NCBI opens these entries with the database, LOSAT only for a taxid list. A volume of
    // the database that the LMDB file does not list, or lists with another number of OIDs,
    // leaves NCBI's OIDs off by the difference; LOSAT refuses such a database.
    /// The LMDB file of each run of volumes that name the same one, with the volumes of the
    /// run (`CSeqDBLMDBSet`'s grouping, cited above).
    fn lmdb_groups(&self) -> Vec<(PathBuf, std::ops::Range<usize>)> {
        let mut groups: Vec<(PathBuf, std::ops::Range<usize>)> = Vec::new();
        for (at, volume) in self.volumes.iter().enumerate() {
            let name = volume.index().lmdb_file.as_deref().unwrap_or_default();
            let lmdb = volume
                .base()
                .with_file_name(crate::utils::os_bytes::path_from_bytes(name));
            match groups.last_mut() {
                Some((path, range)) if *path == lmdb => range.end = at + 1,
                _ => groups.push((lmdb, at..at + 1)),
            }
        }
        groups
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:399-405
    /// `GetLMDBFileNames`: the LMDB file of each run of volumes of a version 5 database that
    /// name the same one, in the order of the volumes (a file twice when two runs name it);
    /// none for version 4. The files are not read.
    pub fn lmdb_file_names(&self) -> Vec<PathBuf> {
        if !self.is_version5() {
            return Vec::new();
        }
        self.lmdb_groups()
            .into_iter()
            .map(|(path, _)| path)
            .collect()
    }

    /// The LMDB files of a version 5 database with the volumes each one lists.
    pub(super) fn lmdb_entries(&self) -> Result<Vec<LmdbEntry>, SeqDbError> {
        let groups = self.lmdb_groups();
        let mut entries = Vec::with_capacity(groups.len());
        for (lmdb, range) in groups {
            let file = DbFile::open(&lmdb).map_err(lmdb_file_error)?;
            let env = Env::open(&file).map_err(|message| lmdb_open_error(&file, message))?;
            let listed = tax_files::volumes(&env).map_err(|error| {
                error.into_seqdb(|message| {
                    if tax_files::is_ncbi_volume_error(&message) {
                        SeqDbError::NotFound(message.into_bytes())
                    } else {
                        SeqDbError::Format {
                            path: lmdb.clone(),
                            message,
                        }
                    }
                })
            })?;
            let held = &self.volumes[range.clone()];
            if held.len() > listed.len() {
                return Err(SeqDbError::NotFound(VOLUME_COUNT_ERROR.into()));
            }
            let mut next = 0;
            let mut volumes = Vec::with_capacity(listed.len());
            let (mut total, mut included) = (0u64, 0u64);
            for (name, count) in &listed {
                total += u64::from(*count);
                let holds = held.get(next).is_some_and(|volume| {
                    volume
                        .base()
                        .file_name()
                        .is_some_and(|base| base.as_encoded_bytes() == name.as_slice())
                });
                if holds {
                    if held[next].num_oids() != *count {
                        return Err(SeqDbError::Unsupported {
                            path: lmdb.clone(),
                            message: format!(
                                "the LMDB file gives volume {} {count} sequences, its index {}",
                                String::from_utf8_lossy(name),
                                held[next].num_oids()
                            ),
                        });
                    }
                    included += u64::from(*count);
                    next += 1;
                }
                let max_oid = u32::try_from(total).map_err(|_| SeqDbError::Unsupported {
                    path: lmdb.clone(),
                    message: "the LMDB file lists more sequences than a signed 32-bit count"
                        .to_string(),
                })?;
                volumes.push((max_oid, if holds { 0 } else { *count }));
            }
            if included == 0 {
                return Err(SeqDbError::NotFound(VOLUME_NAME_ERROR.into()));
            }
            if next < held.len() {
                return Err(SeqDbError::Unsupported {
                    path: lmdb.clone(),
                    message: format!(
                        "the LMDB file does not list volume {}",
                        held[next].base().display()
                    ),
                });
            }
            entries.push(LmdbEntry {
                lmdb,
                oid_start: self.volume_starts[range.start],
                partial: included != total,
                num_oids: total as u32,
                volumes,
            });
        }
        Ok(entries)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:315-333
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdblmdbset.cpp:349-351
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:109-125
    // A union of the found taxids is what the loop keeps: they are taxids of `tax_ids`.
    /// Restricts the database to the OIDs of a taxid list, as NCBI's `CSeqDB` opened with the
    /// list: from here on the ordinal ids are those of the included OIDs, the number of
    /// sequences and the total length are NCBI's scanned totals of them (`num_seqs`,
    /// `total_length`, `max_length`; the reports print `printed_num_seqs` and
    /// `printed_total_length`, the whole database's), and the definition lines are those of
    /// NCBI's filtered header read. Returns the warning NCBI writes (a negative list of
    /// taxids the database does not hold, which excludes nothing), to be written as it is.
    ///
    /// Errors: NCBI's (`SeqDbError::NotFound`): a version 4 database, a positive list that
    /// selects no OID, a missing lookup file; the reader's for damaged lookup files.
    pub fn restrict_taxids(&mut self, list: &TaxidList) -> Result<Option<Vec<u8>>, SeqDbError> {
        let (negative, tax_ids) = match list {
            TaxidList::Positive(tax_ids) => (false, tax_ids),
            TaxidList::Negative(tax_ids) => (true, tax_ids),
        };
        if !self.is_version5() {
            return Err(SeqDbError::NotFound(V4_TAXID_ERROR.to_string().into()));
        }
        let entries = self.lmdb_entries()?;
        let (oids, found) = entries_tax_oids(&entries, negative, tax_ids)?;
        let mut warning = None;
        if oids.is_empty() {
            let text = tax_ids_not_found(&entries);
            if !negative {
                return Err(SeqDbError::NotFound(text));
            }
            warning = Some(text);
        }
        let num_oids = self.volume_starts[self.volumes.len()];
        let fill = if negative { u64::MAX } else { 0 };
        let mut words = vec![fill; (num_oids as usize).div_ceil(64)];
        for oid in oids {
            if oid >= num_oids {
                return Err(SeqDbError::Format {
                    path: entries[0].lmdb.clone(),
                    message: format!(
                        "the taxid lookup gives OID {oid} of a database of {num_oids} sequences"
                    ),
                });
            }
            let (word, bit) = ((oid / 64) as usize, oid % 64);
            if negative {
                words[word] &= !(1u64 << bit);
            } else {
                words[word] |= 1u64 << bit;
            }
        }
        // A positive list is the object that `tax_ids.swap` with `rv_tax_ids` narrows to the found
        // ids; a negative one is a copy (uniform_search.cpp:277-280), so the search's list keeps
        // all.
        self.list_taxids = if negative { tax_ids } else { &found }
            .iter()
            .copied()
            .collect();
        // A negative list of taxids the database does not hold excludes no OID and tests no
        // taxid, but it is still the volumes' list (`DeflineRule`).
        let rule = if negative {
            DeflineRule::NegativeTaxids(found)
        } else {
            DeflineRule::Taxids(found)
        };
        self.apply_filter(Some(words), num_oids, Some(rule))?;
        Ok(warning)
    }

    /// Restricts the database to `list` ([`SeqDb::restrict_taxids`], [`SeqDb::restrict_seqids`]);
    /// the warnings NCBI writes go to `warnings`, before an error too.
    pub fn restrict(
        &mut self,
        list: &DbList,
        warnings: &mut Vec<ListWarning>,
    ) -> Result<(), SeqDbError> {
        self.positive_list = match list {
            DbList::Taxids(list) => matches!(list, TaxidList::Positive(_)),
            DbList::Seqids(list) => !list.negative,
        };
        match list {
            DbList::Taxids(list) => {
                warnings.extend(self.restrict_taxids(list)?.map(ListWarning::Log));
                Ok(())
            }
            DbList::Seqids(list) => self.restrict_seqids(list, warnings),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:235-240
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:245-252
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:294-303
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:323-329
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:328-331
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:365-374
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.cpp:376-381
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.cpp:75-88
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbbitset.cpp:90-100
    // A text list that lost every id (all `gi|N`): NCBI warns `Empty seqid list`. A negative
    // one excludes nothing. A positive one is empty (`gis.Empty()`), and `x_ApplyUserGiList`
    // clears every OID of `m_AllBits` (the code above): no OID is searched and each query
    // gets `Filtering resulted in an empty database.` (two or more volumes, or one volume
    // whose every path in the alias files crosses an OID range that cuts it). When `m_AllBits`
    // is the special all-set (`alias::AliasOids::all_bits_special_set`), the clearing (`ClearBit` below
    // 24 OIDs, a `memset` from 24; a release build compiles the `_ASSERT` that `m_Special` is
    // `eNone` out) writes into the array its copy left behind, and the set still answers every OID:
    // NCBI searches the whole database
    // (probe S17, S18, R08, R09 of dropin-survey-20261029/DI8: 500 hits of 514 sequences, the
    // statistics of the whole database; line DBID gisl: 1, 2, 3 volumes and 14 alias shapes, 6
    // programs). LOSAT does what the runs show: the warnings, every OID kept or none, and, as
    // NCBI's `CSeqDBImpl` still scans the totals for a list (`x_ScanTotals`, approximate
    // nucleotide lengths, the shortest length), the scanned totals.
    // The ids of a binary list are used as written (the alias tool wrote them as the
    // processing writes them). The lists the volumes keep for the filtered header read are
    // the processed ids.
    /// Restricts the database to the OIDs of a seqid list, as NCBI's `CSeqDB` opened with the
    /// list (`CSeqDBGiListSet::x_ResolvePositiveList` / `x_ResolveNegativeList`, then
    /// `CSeqDBOIDList`): the ordinal ids, totals and definition lines as
    /// [`SeqDb::restrict_taxids`] leaves them. The warnings NCBI writes (a text list on a
    /// version 5 database, a list with no id, a binary list of another database) go to
    /// `warnings`, before an error too.
    ///
    /// Errors: NCBI's (`SeqDbError::NotFound`): a binary list on a version 4 database, an
    /// id that does not parse, a database without the accession table; the reader's for
    /// damaged lookup files.
    pub fn restrict_seqids(
        &mut self,
        list: &SeqidList,
        warnings: &mut Vec<ListWarning>,
    ) -> Result<(), SeqDbError> {
        if !self.is_version5() {
            return self.restrict_seqids_v4(list, warnings);
        }
        let ids: Vec<Vec<u8>> = if list.file.info.is_v4 {
            warnings.push(ListWarning::Posted(TEXT_LIST_WARNING.to_string()));
            let mut posted = Vec::new();
            let processed = seqidlist::process_for_v5(&list.file.ids, &mut |warning| {
                posted.push(ListWarning::SeqId(warning));
            });
            warnings.append(&mut posted);
            let ids = processed.map_err(|invalid| {
                warnings.push(ListWarning::Error(invalid.error.message.clone()));
                SeqDbError::NotFound(invalid.message())
            })?;
            if ids.is_empty() {
                // The emptied list is still the volumes' list (`DeflineRule`): no OID is
                // excluded by its ids, and its lines without a Seq-id are erased.
                warnings.push(ListWarning::Posted(EMPTY_LIST_WARNING.to_string()));
                let num_oids = self.volume_starts[self.volumes.len()];
                if list.negative {
                    let words = vec![u64::MAX; (num_oids as usize).div_ceil(64)];
                    return self.apply_filter(
                        Some(words),
                        num_oids,
                        Some(DeflineRule::NegativeSeqids(Vec::new())),
                    );
                }
                return self.apply_filter(None, num_oids, Some(DeflineRule::Seqids(Vec::new())));
            }
            ids
        } else {
            if list.file.info.db_vol_length > 0
                && list.file.info.db_vol_length != self.volume_set_length
            {
                warnings.push(ListWarning::Posted(DB_INFO_WARNING.to_string()));
            }
            let mut ids = list.file.ids.clone();
            ids.sort();
            ids
        };
        let num_oids = self.volume_starts[self.volumes.len()];
        let words = self.seqid_words_v5(&ids, list.negative, num_oids)?;
        let rule = if list.negative {
            DeflineRule::NegativeSeqids(ids)
        } else {
            DeflineRule::Seqids(ids)
        };
        self.apply_filter(Some(words), num_oids, Some(rule))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1580-1586
    // `GetSequenceLengthAndId` takes the Seq-ids of `CSeqDB::GetSeqIDs`, the filtered header
    // (seqdbvol.cpp:2106-2111): a sequence that a seqid list keeps while the header read
    // removes every line has none (a positive list's lookup found a pdb molecule without its
    // chain, a local id in other letters, a version 4 Swiss-Prot name or a `BL_ORD_ID` tag of a
    // later volume; a negative list keeps a line only with an id that is neither a gi nor a ti
    // nor listed), and its hits get no Seq-align: no report, prefetch, taxonomy lookup or PSSM
    // sees them, and the post is not printed. The sequence was searched and its length counts
    // in the statistics (`drops_lineless_hits`; under a negative list the search's fetches
    // exclude it earlier, `excludes_lineless_at_fetch`). (Line DBID items lineless and
    // lineless2; DI-8 `DECISIONS.md` 17 had refused such a list before the search.) The alias
    // files' `OID_MASK_TYPE` erases the lines of model sequences in the same read (line RG D1b-2):
    // a sequence that the OID list keeps with only such lines is searched, counts, and its hits
    // leave the reports the same way (probes A5, A6 of line RG D1b); it fills no id set, so the
    // fetches exclude nothing for it. A taxid list fills no id set either (`x_InitIdSet`): a
    // sequence whose lines it erases leaves at the Seq-align alone (review 2 of line RG D1c, R2).
    /// Whether the database's reads may keep a sequence without any of its lines (a list of any
    /// kind, or the alias files' erasure of model lines): such a sequence is
    /// searched and counts in the statistics, and its hits that reach NCBI's Seq-aligns leave the
    /// reports there (`keeps_line`, `algorithm::common::database::SeqAligned`).
    pub fn drops_lineless_hits(&self) -> bool {
        self.volumes.iter().any(|volume| volume.filters_lines())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:302-313
    // The fetches of the gapped traceback (blast_traceback.c:1567-1573), of the
    // composition-based redo (blast_kappa.c:1378, a match skipped at 3560-3575) and of the
    // ungapped re-evaluation of a nucleotide subject (blast_hits.c:2651-2665, the subject's
    // HSPs freed) ask for this check.
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1929-1943
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:108-113
    // `x_InitIdSet` (seqdbimpl.cpp:1909-1946) leaves the id set blank for a positive list of
    // Seq-ids (`-seqidlist` reads only strings, `CSeqDBFileGiList::eSiList`, and the branch
    // above it fills the set from a user list's gis or tis only), and fills it from the ids
    // of a negative list as read: the constructor makes the set before its first OID list
    // processes the list (`s_ProcessSeqIDsForV5`, which drops the gis on a version 5
    // database), so a list of gis only that the processing empties fills it too (review 1 of
    // line RG D1c-1b, F1), and every negative list does (a list file without an id is no list,
    // `DbFilterArgs::resolve`). The fetches with
    // `check_oid_exclusion` then exclude a sequence without a line: a negative list's such
    // sequence is searched, takes its place in the preliminary hit lists of a gapped search,
    // and leaves the search at the traceback's or the composition-based redo's fetch, before
    // the final hit lists are made (its hits take no place there); in an ungapped search of a
    // nucleotide database it leaves at the re-evaluation, before the preliminary hit lists. A
    // protein subject of an ungapped search is fetched by none of them and reaches the
    // Seq-align (`drops_lineless_hits`).
    /// Whether NCBI's fetches with `check_oid_exclusion` exclude the sequences that keep no
    /// line (`keeps_line`): the database's rule is a negative seqid list
    /// (`algorithm::common::database::FetchExcluded`).
    pub fn excludes_lineless_at_fetch(&self) -> bool {
        self.volumes
            .iter()
            .any(|volume| volume.has_negative_seqid_list())
    }

    /// Whether NCBI turns a hit on `ordinal` into a Seq-align: false where the database's list
    /// or its erasure of model lines (`OID_MASK_TYPE`) keeps the sequence but none of its
    /// definition lines, and where no line it keeps holds a Seq-id (an edited header; the lists
    /// erase such lines); the reader's error where a tag cannot be read as NCBI's header read
    /// adjusts it.
    pub fn keeps_line(&self, ordinal: u32) -> Result<bool, SeqDbError> {
        self.line_read(ordinal).map(|(kept, _)| kept)
    }

    /// [`SeqDb::keeps_line`], read as NCBI's `GetSeqIDs` reads the header, with whether the
    /// record holds a UTF8String where a VisibleString is expected (NCBI's read then posts its
    /// warning, `Utf8Reads`): false where the volume filters no line, whose walk of the record
    /// reads no string (`Volume::line_read`).
    pub fn line_read(&self, ordinal: u32) -> Result<(bool, bool), SeqDbError> {
        let at = self.find_volume(ordinal)?;
        self.volumes[at.volume].line_read(at.vol_oid)
    }

    /// The header file of the volume that holds `ordinal` (`None` outside the database, or where
    /// the file does not open): LOSAT's refusal of a header that NCBI's fetch cannot read names
    /// it (`algorithm::common::database::FetchExcluded::refusal`).
    #[cold]
    pub fn header_file(&self, ordinal: u32) -> Option<PathBuf> {
        let at = self.find_volume(ordinal).ok()?;
        let file = self.volumes[at.volume].hdr_file().ok()?;
        Some(file.path().to_path_buf())
    }

    /// [`SeqDb::keeps_line`] of the sequence of OID `oid`, included or not (the Seq-aligns of an
    /// RPS search read their profiles' headers by OID).
    pub fn keeps_line_of_oid(&self, oid: u32) -> Result<bool, SeqDbError> {
        let at = self.locate_oid(oid)?;
        self.volumes[at.volume]
            .line_read(at.vol_oid)
            .map(|(kept, _)| kept)
    }

    /// Whether the database's reads remove a definition line of the sequence of OID `oid`,
    /// included or not (`Volume::erases_a_line`: nothing posted or refused for its lines).
    pub fn erases_a_line_of_oid(&self, oid: u32) -> Result<bool, SeqDbError> {
        let at = self.locate_oid(oid)?;
        self.volumes[at.volume].erases_a_line(at.vol_oid)
    }

    /// The bitmap of the OIDs a seqid list keeps in a version 5 database (`ids`: sorted and
    /// unique, as the lookup wants them).
    fn seqid_words_v5(
        &self,
        ids: &[Vec<u8>],
        negative: bool,
        num_oids: u32,
    ) -> Result<Vec<u64>, SeqDbError> {
        let entries = self.lmdb_entries()?;
        let word_count = (num_oids as usize).div_ceil(64);
        if negative {
            let mut words = vec![u64::MAX; word_count];
            for entry in &entries {
                for oid in entry.negative_seqids_to_oids(ids)? {
                    if oid >= num_oids {
                        return Err(bad_lookup_oid(&entry.lmdb, oid, num_oids));
                    }
                    words[(oid / 64) as usize] &= !(1u64 << (oid % 64));
                }
            }
            return Ok(words);
        }
        let mut words = vec![0u64; word_count];
        let mut oids = Vec::new();
        entries_accessions_to_oids(&entries, ids, &mut oids)?;
        for oid in oids.into_iter().flatten() {
            if oid < num_oids {
                words[(oid / 64) as usize] |= 1u64 << (oid % 64);
            }
        }
        Ok(words)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:265-301
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:345-347
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbgilistset.cpp:369-375
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2540-2550
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:375-402
    // A text list is not rewritten for a version 4 database (`s_VerifySeqidlist` returns
    // false): every id, simplified and lower-cased, is a key of the string index; an id that
    // devolves to a number (a gi) becomes the empty key, which no line has (probe R12, R13).
    // An exception of `SeqDB_SimplifyAccession` (an id such as `sp||`) is not a database
    // error in NCBI (`Error: <what>`, exit status 255); LOSAT refuses it with the text.
    /// [`SeqDb::restrict_seqids`] on a version 4 database: the string ISAM lookups of each
    /// volume (`CSeqDBVol::IdsToOids`), after `PreprocessIdsForISAMSiLookup`.
    fn restrict_seqids_v4(
        &mut self,
        list: &SeqidList,
        _warnings: &mut [ListWarning],
    ) -> Result<(), SeqDbError> {
        if !list.file.info.is_v4 {
            return Err(SeqDbError::NotFound(
                BINARY_LIST_V4_ERROR.to_string().into(),
            ));
        }
        let protein = matches!(self.seqtype, SeqType::Protein);
        let mut ids: Vec<Vec<u8>> = Vec::with_capacity(list.file.ids.len());
        for id in &list.file.ids {
            ids.push(isam::simplify_accession_for_isam(id).map_err(|error| {
                SeqDbError::NcbiThrows {
                    path: self.volumes[0].base().to_path_buf(),
                    message: format!(
                        "the id {} of the seqid list throws in NCBI's SeqDB_SimplifyAccession: {}",
                        String::from_utf8_lossy(id),
                        error.to_string().trim_end()
                    ),
                }
            })?);
        }
        let string_isam = |at: usize| -> Result<isam::StringIsam, SeqDbError> {
            let base = self.volumes[at].base();
            if !isam::StringIsam::index_exists(base, protein) {
                return Err(SeqDbError::NotFound(
                    [
                        b"SI list specified but no ISAM file found for SI in ".as_slice(),
                        &crate::utils::os_bytes::path_bytes(base),
                    ]
                    .concat(),
                ));
            }
            isam::StringIsam::open_volume(base, protein)
        };
        let num_oids = self.volume_starts[self.volumes.len()];
        let word_count = (num_oids as usize).div_ceil(64);
        let words = if list.negative {
            let mut out = isam::NegativeOids::new();
            for at in 0..self.volumes.len() {
                string_isam(at)?.negative_scan(
                    self.volume_starts[at],
                    self.volume_starts[at + 1],
                    &mut ids,
                    &mut out,
                )?;
            }
            let mut words = vec![u64::MAX; word_count];
            let max = out.num_oids().min(num_oids);
            for oid in 0..num_oids {
                if oid >= max || !out.get_oid_status(oid) {
                    words[(oid / 64) as usize] &= !(1u64 << (oid % 64));
                }
            }
            words
        } else {
            let mut entries: Vec<isam::SiOid> =
                ids.iter().map(|id| isam::SiOid::new(id.clone())).collect();
            let mut sorted = false;
            let mut order: Vec<usize> = (0..self.volumes.len()).collect();
            order.sort_by_key(|&at| std::cmp::Reverse(self.volumes[at].num_oids()));
            for at in order {
                string_isam(at)?.ids_to_oids(self.volume_starts[at], &mut entries, &mut sorted)?;
            }
            let mut words = vec![0u64; word_count];
            for oid in entries.iter().filter_map(|entry| entry.oid) {
                if oid < num_oids {
                    words[(oid / 64) as usize] |= 1u64 << (oid % 64);
                }
            }
            words
        };
        ids.sort();
        let rule = if list.negative {
            DeflineRule::NegativeSeqids(ids)
        } else {
            DeflineRule::Seqids(ids)
        };
        self.apply_filter(Some(words), num_oids, Some(rule))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:100-113
    // The user list narrows the OID list of the alias files' ranges (`x_ApplyUserGiList`
    // intersects `m_AllBits`; a negative list clears bits of it): an OID is listed when both
    // include it.
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:144-151
    // The constructor scans the totals with the list (`x_ScanTotals`, cited at `scan_totals`):
    // its first `x_CheckOrFindOID` sets the list up, with the iteration range of the only alias
    // file of a `-db` (`OidFilter::walk` of a list not set up), and then the constructor's
    // `SetIterationRange` from 0 to 0 sets the range back to every OID, so the search reads every
    // listed OID (line RG D1a-2, survey A M4, M6).
    /// Makes `words` (bit `oid` set for an OID the list includes; `None` for a positive list
    /// left without ids, which clears every OID of NCBI's bit set: every OID stays where that set
    /// is the special all-set, `restrict_seqids`) and the alias files' ranges and lists under the
    /// search's list `rule` (`SeqDb::translate_alias_lists`) the database's OID list with
    /// the iteration range of every OID: the ordinal ids, NCBI's scanned totals and the
    /// definition lines `rule` keeps.
    fn apply_filter(
        &mut self,
        words: Option<Vec<u64>>,
        num_oids: u32,
        rule: Option<DeflineRule>,
    ) -> Result<(), SeqDbError> {
        if rule.is_some() {
            self.translate_alias_lists(rule.as_ref());
        }
        let mut words = words.unwrap_or_else(|| {
            let kept = if self.all_bits_special_set {
                u64::MAX
            } else {
                0
            };
            vec![kept; (num_oids as usize).div_ceil(64)]
        });
        // The constructor builds the OID list to scan: an alias file's OID list that fails is
        // its error, after the list's own (`CSeqDBGiListSet` resolves the list first).
        self.build_oid_list()?;
        let range = match &self.filter {
            Some(alias) => {
                for (word, kept) in words.iter_mut().zip(&alias.words) {
                    *word &= kept;
                }
                self.alias_range
            }
            None => (0, num_oids),
        };
        if let Some(rule) = rule.map(Arc::new) {
            for (volume, &start) in self.volumes.iter_mut().zip(&self.volume_starts) {
                volume.set_defline_rule(Arc::clone(&rule), start);
            }
        }
        self.list_applied = true;
        self.constructor_scan(OidFilter::new(words, num_oids).with_range(range))
            .map(|_| ())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:117-132
    // The alias files' own scan (`alias::needs_totals_scan`), as the constructor makes it with a
    // list (`apply_filter`): the reports, which open the database again without a list, scan
    // the same OIDs, so the printed counts are this scan (probes B1, B2, U1 of line RG D1b).
    /// Scans the totals of the alias files' OID list `list` (with its iteration range) as NCBI's
    /// constructor does without a list: the counts and the printed counts are the scan.
    pub(super) fn scan_alias(&mut self, list: OidFilter) -> Result<(), SeqDbError> {
        let scan = self.constructor_scan(list)?;
        self.alias_scanned = true;
        self.printed_num_seqs = scan.count as i32;
        self.printed_total_length = scan.total;
        Ok(())
    }

    /// The constructor's scan of the OID list `list` (`x_ScanTotals` over the OIDs of its first
    /// `x_CheckOrFindOID` loop, `OidFilter::walk`) and its `SetIterationRange` from 0 to 0: the
    /// list, with every OID in range, becomes the database's, and the counts are the scan.
    fn constructor_scan(&mut self, list: OidFilter) -> Result<ScanTotals, SeqDbError> {
        let num_oids = self.num_oids();
        let scan = self.scan_walk(&list, list.walk(false))?;
        let filter = list.with_range((0, num_oids));
        // `x_ScanTotals` counts the OIDs in an `int` (the volume set holds fewer than 2^31 OIDs:
        // `open_tree`) and keeps the shortest length without the `<= 0` fallback.
        self.num_seqs = filter.count();
        self.counts = super::Counts {
            num_seqs: scan.count as i32,
            total_length: scan.total,
            max_length: scan.longest,
            min_length: scan.shortest as i32,
        };
        self.oid_list_set_up.store(true, Ordering::Relaxed);
        // The scan is that of `included_length` when it visited every included OID.
        self.set_included_length((scan.count == filter.count()).then_some(scan.total));
        self.filter = Some(Box::new(filter));
        self.update_first_oid_length()?;
        Ok(scan)
    }

    /// The database restricted to the OIDs that `keep` takes, as a list that names them does:
    /// a test's list on a fixture without the ids or taxids of one (several volumes).
    #[cfg(test)]
    pub(crate) fn x_keep_oids(&mut self, keep: impl Fn(u32) -> bool) {
        let num_oids = self
            .volumes
            .iter()
            .map(|volume| volume.num_oids())
            .sum::<u32>();
        let mut words = vec![0u64; (num_oids as usize).div_ceil(64)];
        for oid in (0..num_oids).filter(|&oid| keep(oid)) {
            words[(oid / 64) as usize] |= 1 << (oid % 64);
        }
        self.apply_filter(Some(words), num_oids, None).unwrap();
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1569-1593
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:325-344
    // A protein's length is exact (`GetSeqLengthProt`: the bytes up to the next sequence
    // less the NUL between them); a nucleotide's is approximate, from its bytes and the
    // volume-relative OID, not from the count in its last byte.
    /// Keeps `length` as the length of the included sequences (`included_length`), or makes it
    /// at the first call (`None`). Only the parallel searches' schedules read it.
    pub(super) fn set_included_length(&mut self, length: Option<u64>) {
        #[cfg(any(test, feature = "parallel"))]
        {
            self.included_length =
                length.map_or_else(std::sync::OnceLock::new, std::sync::OnceLock::from);
        }
        #[cfg(not(any(test, feature = "parallel")))]
        let _ = length;
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:846-862
    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:545-555 (at
    // `blastdbcmd::write_info_head`)
    // A database whose constructor scanned the totals scans them again with exact lengths,
    // over the OIDs of its list set up with every OID in range (`x_CheckOrFindOID`), and keeps
    // the longest length of that scan, which blastdbcmd `-info` prints after the total.
    /// `GetExactTotalLength`, with the `GetMaxLength` after it: the total and the longest exact
    /// length of the included OIDs where the database scanned its totals (a list, an alias
    /// file's `OIDLIST`), else [`SeqDb::total_length`] and [`SeqDb::max_length`].
    pub fn exact_totals(&self) -> Result<(u64, u32), SeqDbError> {
        let filter = match &self.filter {
            Some(filter) if self.alias_scanned || self.list_applied => filter,
            _ => return Ok((self.counts.total_length, self.counts.max_length)),
        };
        let (mut total, mut longest) = (0u64, 0u32);
        for oid in filter.walk(true).oids(filter) {
            let at = self.locate_oid(oid)?;
            let length = match self.seqtype {
                SeqType::Protein => self.volumes[at.volume].protein_length(at.vol_oid)?,
                SeqType::Nucleotide => self.nucleotide_length_at(at)? as u32,
            };
            total += u64::from(length);
            longest = longest.max(length);
        }
        Ok((total, longest))
    }

    /// NCBI's totals of the OIDs of the walk `walk` of the list `list` (`OidFilter::walk`).
    fn scan_walk(&self, list: &OidFilter, walk: OidWalk) -> Result<ScanTotals, SeqDbError> {
        let mut totals = self.scan_totals(list, walk.from, walk.to)?;
        for oid in walk.first.into_iter().chain(walk.last) {
            totals.add(self.scan_totals(list, oid, oid + 1)?);
        }
        Ok(totals)
    }

    /// NCBI's totals of the listed OIDs of `from..to` of `list`: their number, the total
    /// length, the longest length and the shortest (`INT4_MAX` when there is none, as
    /// `min_count` starts).
    fn scan_totals(&self, list: &OidFilter, from: u32, to: u32) -> Result<ScanTotals, SeqDbError> {
        let mut totals = ScanTotals::default();
        for (at, volume) in self.volumes.iter().enumerate() {
            let start = self.volume_starts[at];
            let mut first = from.max(start) - start;
            let stop = to.min(start + volume.num_oids()).saturating_sub(start);
            while first < stop {
                let end = stop.min(first.saturating_add(SCAN_WINDOW));
                if list.listed_below(start + end) == list.listed_below(start + first) {
                    first = end;
                    continue;
                }
                let index = volume.index();
                let seq = index.sequence_offsets(first, end)?;
                let amb = index.ambiguity_offsets(first, end)?;
                for i in 0..(end - first) as usize {
                    let oid = first + i as u32;
                    if !list.listed(start + oid) {
                        continue;
                    }
                    totals.add_length(approx_length(self.seqtype, &seq, &amb, i, oid));
                }
                first = end;
            }
        }
        Ok(totals)
    }

    // No NCBI counterpart: the schedules of the searches; it changes no value NCBI computes.
    /// The length of the included sequences ([`SeqDb::num_seqs`] of them, those a search
    /// reads), which the searches' schedules read: the volumes' total length when every OID is
    /// included, else the total of NCBI's scan of the included OIDs (`scan_totals`, made at the
    /// first call). NCBI's `GetTotalLength` ([`SeqDb::total_length`]) is the walk of the alias
    /// files' `LENGTH` without a list, and the scan of `x_CheckOrFindOID`'s OIDs with one, which
    /// can differ from it (an OID range, a `LENGTH` key, the iteration range).
    #[cfg(any(test, feature = "parallel"))]
    pub fn included_length(&self) -> Result<u64, SeqDbError> {
        let Some(filter) = &self.filter else {
            return Ok(self.volume_set_length);
        };
        if let Some(&length) = self.included_length.get() {
            return Ok(length);
        }
        let (begin, end) = filter.range();
        let length = self.scan_totals(filter, begin, end)?.total;
        Ok(*self.included_length.get_or_init(|| length))
    }

    // No NCBI counterpart: NCBI's database reader hands its search threads chunks of sequences
    // whose exact lengths it reads from the mapped files (`CSeqDBImpl::x_FillSeqBuffer`); LOSAT
    // cuts its chunks before the search from these lengths (`scan_totals`' NCBI lengths), without
    // a read of the sequence file. Scheduling only.
    /// Calls `visit` with the length of every sequence in ordinal order as `scan_totals` reads it
    /// (NCBI's `GetSeqLengthApprox` for a nucleotide, from the index alone), `SCAN_WINDOW` index
    /// entries at a time: the chunks of TBLASTN's `-db` search (`XDbStream::x_approx_lengths`,
    /// `x_kept_cut`). An error of the index ends the visit; the sequences before it were visited.
    #[cfg(any(feature = "parallel", test))]
    pub fn x_approx_lengths(&self, mut visit: impl FnMut(u32)) -> Result<(), SeqDbError> {
        for (at, volume) in self.volumes.iter().enumerate() {
            let start = self.volume_starts[at];
            let mut first = 0;
            while first < volume.num_oids() {
                let end = volume.num_oids().min(first.saturating_add(SCAN_WINDOW));
                if let Some(filter) = &self.filter {
                    if filter.rank(start + end) == filter.rank(start + first) {
                        first = end;
                        continue;
                    }
                }
                let index = volume.index();
                let seq = index.sequence_offsets(first, end)?;
                let amb = index.ambiguity_offsets(first, end)?;
                for i in 0..(end - first) as usize {
                    let oid = first + i as u32;
                    if self
                        .filter
                        .as_ref()
                        .is_some_and(|filter| !filter.contains(start + oid))
                    {
                        continue;
                    }
                    visit(approx_length(self.seqtype, &seq, &amb, i, oid));
                }
                first = end;
            }
        }
        Ok(())
    }
}

/// What NCBI's `x_ScanTotals` adds up of the OIDs it visits.
#[derive(Debug, Clone, Copy)]
struct ScanTotals {
    count: u32,
    total: u64,
    longest: u32,
    shortest: u32,
}

impl Default for ScanTotals {
    fn default() -> Self {
        Self {
            count: 0,
            total: 0,
            longest: 0,
            shortest: i32::MAX as u32,
        }
    }
}

impl ScanTotals {
    fn add_length(&mut self, length: u32) {
        self.count += 1;
        self.total += u64::from(length);
        self.longest = self.longest.max(length);
        self.shortest = self.shortest.min(length);
    }

    fn add(&mut self, other: ScanTotals) {
        self.count += other.count;
        self.total += other.total;
        self.longest = self.longest.max(other.longest);
        self.shortest = self.shortest.min(other.shortest);
    }
}

/// The length of sequence `i` of a window of the index (`seq`, `amb`: its sequence and
/// ambiguity offsets; `oid`: its volume-relative OID) as NCBI's scan reads it: exact for a
/// protein, `GetSeqLengthApprox` for a nucleotide (`scan_totals`).
fn approx_length(seqtype: SeqType, seq: &[u32], amb: &[u32], i: usize, oid: u32) -> u32 {
    let length = match seqtype {
        SeqType::Protein => i64::from(seq[i + 1]) - i64::from(seq[i]) - 1,
        SeqType::Nucleotide => 4 * (i64::from(amb[i]) - i64::from(seq[i]) - 1) + i64::from(oid & 3),
    };
    u32::try_from(length.max(0)).unwrap_or(u32::MAX)
}

#[cfg(test)]
#[path = "filter_tests.rs"]
mod tests;
