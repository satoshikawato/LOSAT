//! A database sequence with parsed Seq-ids in the reports: what each format prints of its
//! definition lines (`seqid::DeflineSet`), chosen and printed by the functions of `seqid`.
//!
//! NCBI's formatters find a subject again by its Seq-align id: the scope gives the Bioseq of
//! the definition line that holds the id (`DeflineSet::bioseq_defline`, the "target") with
//! the whole set as a user object. The ids and title of the target are the Bioseq's (tabular
//! `sseqid`, `sacc`, `saccver`, the description titles of outfmt 0, the XML `Hit_accession`);
//! the formats that show every definition line read the set (the outfmt 0 alignment heading,
//! the XML `Hit_def`, the XML2 and JSON `HitDescr` list, `stitle`).
//!
//! A database without parsed Seq-ids keeps its own path (the title and `gnl|BL_ORD_ID|n`,
//! `algorithm::common::database::stored_title`); `SeqDb::read_header` tells the two apart.
//!
//! LOSAT writes the names of the sequence it searched. NCBI's formatters look the Seq-align id
//! up again in the database, and where the lookup finds another sequence first (an accession
//! without a version after the same accession with one, `gb|X.1|` and `emb|X.1|`, `lcl|X` after
//! `gb|X.1|`: PLAN §9.6 of line DBID) they write that sequence's names and letters, or stop.
//! LOSAT looks the id up as NCBI does (`seqdb::lookup::SeqidLookup`, `found_instead`) and refuses
//! such a subject where a format reads its Bioseq (`fetched`). The two-letter PDB chains of one
//! molecule (`A`, `AA`, `B`) are not such a case (probe of the ID2a review).

use crate::seqdb::DbSequence;
use crate::seqid::defline_set::{
    alignment_line_id, compare, description_id, seq_id_list, tabular_id, Sic, TabularId,
};
use crate::seqid::model::PatentNumber;
use crate::seqid::rank::{find_best, text_score};
use crate::seqid::{label, BlastDefline, DeflineSet, ParsedSeqId};

use super::blast_xml::HitTaxonomy;
use super::defline::{
    bioseq_facts, bioseq_title_from_source, generate_bioseq_defline, BioseqFacts, Defline,
    SourceTaxname, UnknownEncoding,
};
use super::html_ids::{self, HtmlHeading};

/// LOSAT's refusal of the organism report (outfmt 18) at a definition line without a Seq-id
/// (`ParsedSubject::tax_deflines`).
const LINELESS_TAX_REPORT: &str = "has a definition line without a Seq-id, at which NCBI BLAST+'s organism report (outfmt 18) stops for each HSP of the sequence and shows the lines before it once per HSP";

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:144-145
/// `kMaxDescrLength`: the length past which the outfmt 0 description's walk of the definition
/// lines stops (`ParsedSubject::description_meets_lineless_line`).
const MAX_DESCRIPTION_LENGTH: usize = 4096;

/// A reported database sequence with parsed Seq-ids: its definition lines, its Seq-align id
/// and the definition line of its Bioseq.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ParsedSubject {
    set: DeflineSet,
    alignment_id: ParsedSeqId,
    /// The definition line of the subject Bioseq: `DeflineSet::bioseq_defline` of the Seq-id that
    /// NCBI's data loader builds it for (`loader_id`; of `alignment_id` but where another line
    /// holds the Seq-id of its handle, `handle_lookup`).
    target: usize,
    /// Whether NCBI's formatters find no Bioseq by the Seq-align id whatever its lookup gives
    /// (`new`, `fetched`).
    unfetchable: bool,
    /// What NCBI's lookup of the Seq-align id gives (`resolved`, the one producer
    /// `common::database::HeaderReader::resolution`): the answer where it is not this sequence
    /// (`found`, refused by `fetched`), the sequence that it gives first where that is another
    /// one, the header with a UTF8String that it reads (`scope_utf8_read`), and in an API call
    /// after the first the holders listed from the database ([`ScopeResolution`]).
    resolution: ScopeResolution,
    /// Whether NCBI's data loader does not serve this sequence for its Seq-align id although
    /// its lookup gives it (`scope_finds_none`).
    scope_finds_none: bool,
    /// Whether the data loader's build of the Bioseq throws: another id of the sequence holds the
    /// Seq-id of the Seq-align id's handle, but no definition line does as the load reads the
    /// ids (`HandleLookup::Throws`, `handle_load_refusal`).
    load_throws: bool,
    /// The sequence of the database that this subject is (`in_database`; `None` for a `-subject`
    /// record), by which the run's loads name it ([`ScopeWalk`]).
    sequence: Option<DbSequence>,
    /// `-show_gis`: the outfmt 0 description line and heading show it (`description`, `heading`).
    show_gis: bool,
    /// The name of the organism of the subject Bioseq's first source, from which NCBI makes its
    /// title (`with_source`).
    source: SourceTaxname,
    /// Whether this is a `-subject` record (`of_record`): a Bioseq without definition lines.
    record: bool,
    /// Whether it is a record that the reports show as a record without `-parse_deflines`
    /// (`BelievedIds::subject_local`): only its tabular ids are its own (`shown_as_record`).
    local_record: bool,
    /// The strings of the Seq-aligns' `use_this_seqid` (`under_positive_list`); empty without a
    /// positive list.
    use_this_seqid: Vec<Vec<u8>>,
    /// `tabular_saccver`, kept for the reports of later query batches.
    saccver: std::sync::OnceLock<Option<std::sync::Arc<[u8]>>>,
    /// `query_scope_keys`, made at the first report that asks.
    query_scope_keys: std::sync::OnceLock<QueryScopeKeys>,
    /// `skipped_id`, made at the first line that asks.
    skipped_id: std::sync::OnceLock<Vec<u8>>,
}

/// The keys by which NCBI's scope finds a query for the Seq-align id of a reported database
/// sequence (`ParsedSubject::query_scope_keys`; `common::database::refuse_id_of_a_query`), and
/// for its display id, by which the description table looks its title up again
/// (`common::database::refuse_display_id_of_a_query`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QueryScopeKeys {
    /// Its own keys (`defline_set::scope_keys`).
    pub own: Vec<u64>,
    /// The keys of the ids that it matches weakly (`id_scope::weak_keys`).
    pub weak: Vec<u64>,
    /// The keys of the display id (`ParsedSubject::display_id`) and of the ids that it matches
    /// weakly, as `own` and `weak` of the Seq-align id: the same lookup of NCBI's scope.
    pub display: Vec<u64>,
}

/// The `use_this_seqid` strings of a database subject, `subject` when it has parsed Seq-ids: a
/// sequence without them has a `general` id of `BL_ORD_ID`, which is neither a text Seq-id nor a
/// GI (`ParsedSubject::under_positive_list`).
pub fn use_this_seqid(subject: Option<&ParsedSubject>) -> &[Vec<u8>] {
    subject.map_or(&[], ParsedSubject::use_this_seqid)
}

/// `AdjustSeqIdType` with `eGetId_ForceAcc` (`ParsedSubject::sam_id`) of a Bioseq with the ids
/// `ids` and the Seq-align id `alignment_id`.
pub fn forced_accession<'a>(
    alignment_id: &'a ParsedSeqId,
    ids: &'a [ParsedSeqId],
) -> &'a ParsedSeqId {
    let acc_ver =
        |id: &crate::seqid::model::TextseqId| id.accession.is_some() && id.version.is_some();
    if alignment_id.textseq().is_some_and(acc_ver) {
        return alignment_id;
    }
    let mut found = None;
    for id in ids {
        if let Some(text) = id.textseq() {
            found = Some(id);
            if acc_ver(text) {
                break;
            }
        }
    }
    found.unwrap_or(alignment_id)
}

/// What NCBI's lookup of a subject's Seq-align id in its database finds where that is not the
/// subject (`ParsedSubject::found_instead`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FoundInstead {
    /// Another sequence of the database (`DbSequence`) first: NCBI's formatters show it where
    /// its ids compare equal to the Seq-align id, else no sequence.
    Sequence(DbSequence),
    /// Another sequence of the database (`DbSequence`) with the same definition lines but for
    /// their taxids and the same letters (the same record in two databases made with other
    /// taxids), that the lookup gives first or that also holds the Seq-align id: NCBI's
    /// formatters can show it for the subject, which only the reports that write the subject's
    /// taxids show apart (`fetched`).
    OtherTaxids(DbSequence),
    /// Such a sequence, also with the same taxids in another volume, that every report that
    /// reads the Bioseq can show apart, for the reason given.
    CopyApart(DbSequence, CopyApart),
    /// No sequence.
    Nothing,
    /// The lookup fails on a file of the database (`LookupFailure`): NCBI's formatters, and its
    /// prefetch, stop with the `CSeqDBException`.
    Failed(LookupFailure),
    /// The lookup throws this exception on the id itself before it reads a file
    /// (`SeqidLookup::throws_on_id`: a `gnl|PIG|abc` or `gnl|BL_ORD_ID|zz` that makeblastdb stores
    /// with a text tag, `CObject_id::GetId`'s `CInvalidChoiceSelection`; a `gnl|ti|xyz`,
    /// `NStr::StringToInt8`'s `CStringException`): NCBI's prefetch catches it with its warning
    /// (`ParsedSubject::prefetch_fails`), which ends its loads, and its formatters, which then
    /// load the batch's subjects one by one, meet it where they load the subject
    /// (`ParsedSubject::lookup_throws`): each format writes it as where a lookup fails on a file
    /// with the prefetch off (`FormatThrow`, `HitLoad::Throws`; line DBID pfoff-7), and a run that
    /// it ends ends with `Error: ` and its `what()`, exit status 255 (`cli::NativeError::exception`).
    ThrowsOnId(crate::ncbi_exception::NcbiException),
}

// NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:832-939 (cited at `ScopeWalk::load`)
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:379-424
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/cached_sequence.cpp:98-110
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2295-2314
// In an API call after the first (`LoadsBefore::unrecorded`), LOSAT has no record of what the
// calls before loaded into the BLAST database data loader (`m_Ids`, kept for the run), so the
// sequences that hold an id by which NCBI's places can give one of them for a Seq-align id are
// listed from the database (`common::database::HeaderReader::unrecorded_holders`): the loader's,
// those with an id of the Seq-align id's handle on any definition line; the scope's, those whose
// definition line of their own Seq-align id holds an id that the scope takes for it or matches
// weakly (review 3 of line DBID pfoff-4, M1).
/// The other sequences of the database that hold an id by which NCBI's places may give one of
/// them for the Seq-align id of a subject, listed only in an API call after the first, whose
/// calls before LOSAT has no record of ([`ScopeWalk`] with `unrecorded`): any of them may have
/// been loaded there.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ScopeHolders {
    /// The holders of an id of the Seq-align id's own `CSeq_id_Handle`
    /// (`seqid::mask_handle::mask_handle_key`) in any of their definition lines, which the BLAST
    /// database data loader gives for it from the ids of the sequences that it loaded (`m_Ids`)
    /// for the rest of the run, in every format.
    pub loader: Holders,
    /// The holders of an id that NCBI's scope finds for it in its index, in the definition line by
    /// which their load made their Bioseq (an id that it takes for it,
    /// `defline_set::scope_keys`, or one that it matches weakly, `id_scope::weak_keys`), which the
    /// scope holds to the end of the run in outfmt 5 (`RunThrow::keeps_scope`).
    pub scope: Holders,
}

/// A list of [`ScopeHolders`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Holders {
    /// These sequences, in the order of their ordinal ids (none: no other sequence holds such an
    /// id).
    Listed(Vec<DbSequence>),
    /// LOSAT cannot make the list, for this reason (written only in LOSAT's refusal): it fails on
    /// a file, or it asks for the weak matches of a name, or of a `pir` or `prf` id, which a
    /// version 5 database cannot list.
    Unknown(String),
}

impl Default for Holders {
    fn default() -> Self {
        Holders::Listed(Vec::new())
    }
}

impl Holders {
    /// The first of these holders (`Some(None)` where LOSAT cannot list them).
    fn first(&self) -> Option<Option<DbSequence>> {
        match self {
            Holders::Listed(holders) => holders.first().copied().map(Some),
            Holders::Unknown(_) => Some(None),
        }
    }

    /// Why LOSAT cannot list these holders.
    fn unknown(&self) -> Option<&str> {
        match self {
            Holders::Unknown(why) => Some(why),
            Holders::Listed(_) => None,
        }
    }
}

/// What NCBI's lookup of a subject's Seq-align id gives (step 4 of [`ScopeWalk::load`], made
/// once per subject: `common::database::HeaderReader::resolution`, the one producer): the
/// answer where it is not the subject ([`FoundInstead`]), the sequence that it gives first where
/// that is another one (`gives`: also a copy of the record that the reports show as the
/// subject; its blob is what the load then holds), the header with a UTF8String that the lookup
/// reads (`utf8`, posted where the lookup is made), and in an API call after the first the
/// holders that LOSAT lists from the database ([`ScopeHolders`]).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ScopeResolution {
    pub found: Option<FoundInstead>,
    pub gives: Option<DbSequence>,
    /// The keys of the sequence that the lookup gives where its record differs
    /// (`FoundInstead::Sequence`: its blob registers its own lines, `ParsedSubject::shown_keys`).
    pub gives_keys: Option<std::sync::Arc<BlobKeys>>,
    pub utf8: Option<DbSequence>,
    pub unrecorded: ScopeHolders,
}

/// The key of a `CSeq_id_Handle` (`seqid::mask_handle::mask_handle_key`) in NCBI's places that
/// [`ScopeWalk`] keeps: 128 bits of the key's bytes (two 64-bit SipHash values), so a run keeps
/// small keys. Two handles share a key with a chance below 2^-100 per pair.
type HandleKey = u128;

/// [`HandleKey`] of the handle of `id`; `None` for an id that no handle holds again (a
/// Textseq-id with neither accession nor name) or whose handle depends on a tree's state (a
/// `giim` id, which NCBI's handles never take for another).
fn handle_key(id: &ParsedSeqId) -> Option<HandleKey> {
    use std::hash::{Hash, Hasher};
    let key = crate::seqid::mask_handle::mask_handle_key(id).ok()??;
    let mut low = std::hash::DefaultHasher::new();
    key.hash(&mut low);
    let mut high = std::hash::DefaultHasher::new();
    0x9e37_79b9_7f4a_7c15_u64.hash(&mut high);
    key.hash(&mut high);
    Some((u128::from(high.finish()) << 64) | u128::from(low.finish()))
}

/// The keys of a sequence with parsed Seq-ids that NCBI loads (`BlobKeys::of_deflines`): all exact
/// handles of all definition lines flattened into one slice, plus flat per-line scope and weak
/// keys for a holder loaded by `LoaderIds`, and the line of its own Seq-align id (`target`, the
/// line by which a load by that id makes its Bioseq). A sequence that NCBI's places give for a
/// subject that is not the subject or the sequence that the subject's lookup gives (a copy of its
/// record) is another sequence, also a copy of the record, which LOSAT refuses where a format
/// reads it (`ParsedSubject::fetched`).
#[derive(Debug, PartialEq, Eq)]
pub struct BlobKeys {
    sequence: Option<DbSequence>,
    target: usize,
    /// Exact handles (`mask_handle_key`) of every Seq-id of every definition line, flattened.
    handles: Box<[HandleKey]>,
    /// Line index for each entry in `handles`; empty when every line `0..handles.len()` holds
    /// exactly one handle (the common case, where handle `i` belongs to line `i`).
    handle_lines: Box<[u32]>,
    /// Flattened `defline_set::scope_keys` per line, with end offsets in `scope_ends`.
    scope: Box<[u64]>,
    scope_ends: Box<[u32]>,
    /// Flattened `id_scope::kind_keys` per line, with end offsets in `kinds_ends` (both empty when
    /// no line has weak kind keys, such as PDB and GI records).
    kinds: Box<[u64]>,
    kinds_ends: Box<[u32]>,
}

impl BlobKeys {
    /// Builds compact flat keys for `deflines` of `sequence`, targeting `target`.
    fn of_deflines(sequence: Option<DbSequence>, target: usize, deflines: &[BlastDefline]) -> Self {
        let mut handles = Vec::with_capacity(deflines.len());
        let mut handle_lines = Vec::new();
        let mut one_to_one = true;
        let mut scope = Vec::with_capacity(deflines.len());
        let mut scope_ends = Vec::with_capacity(deflines.len());
        let mut kinds = Vec::new();
        let mut kinds_ends = Vec::with_capacity(deflines.len());
        for (line_idx, defline) in deflines.iter().enumerate() {
            let start_h = handles.len();
            for id in &defline.ids {
                handles.extend(handle_key(id));
                crate::seqid::defline_set::scope_keys(id, &mut scope);
                crate::seqid::id_scope::kind_keys(id, &mut kinds);
            }
            let added_h = handles.len() - start_h;
            if one_to_one && added_h != 1 {
                one_to_one = false;
                handle_lines.reserve(handles.len());
                for prev in 0..line_idx {
                    handle_lines.push(prev as u32);
                }
            }
            if !one_to_one {
                for _ in 0..added_h {
                    handle_lines.push(line_idx as u32);
                }
            }
            scope_ends.push(scope.len() as u32);
            kinds_ends.push(kinds.len() as u32);
        }
        if kinds.is_empty() {
            kinds_ends.clear();
        }
        BlobKeys {
            sequence,
            target,
            handles: handles.into_boxed_slice(),
            handle_lines: handle_lines.into_boxed_slice(),
            scope: scope.into_boxed_slice(),
            scope_ends: scope_ends.into_boxed_slice(),
            kinds: kinds.into_boxed_slice(),
            kinds_ends: kinds_ends.into_boxed_slice(),
        }
    }

    /// All exact handles across every definition line, in order.
    fn all_handles(&self) -> &[HandleKey] {
        &self.handles
    }

    /// The first line of this sequence that holds an id of the handle `handle`: the line of the
    /// Bioseq that a load by that handle makes (`GetBioseq` with that target).
    fn line_of(&self, handle: HandleKey) -> Option<usize> {
        let pos = self.handles.iter().position(|&h| h == handle)?;
        Some(if self.handle_lines.is_empty() {
            pos
        } else {
            self.handle_lines[pos] as usize
        })
    }

    /// Exact handles of definition line `line`.
    fn line_handles(&self, line: usize) -> &[HandleKey] {
        if self.handle_lines.is_empty() {
            return self.handles.get(line..line + 1).unwrap_or(&[]);
        }
        let target = line as u32;
        let start = self.handle_lines.partition_point(|&l| l < target);
        let end = self.handle_lines.partition_point(|&l| l <= target);
        &self.handles[start..end]
    }

    /// `defline_set::scope_keys` of definition line `line`.
    fn line_scope(&self, line: usize) -> &[u64] {
        let Some(&end) = self.scope_ends.get(line) else {
            return &[];
        };
        let start = line
            .checked_sub(1)
            .map_or(0, |prev| self.scope_ends[prev] as usize);
        &self.scope[start..end as usize]
    }

    /// `id_scope::kind_keys` of definition line `line`.
    fn line_kinds(&self, line: usize) -> &[u64] {
        if self.kinds.is_empty() {
            return &[];
        }
        let Some(&end) = self.kinds_ends.get(line) else {
            return &[];
        };
        let start = line
            .checked_sub(1)
            .map_or(0, |prev| self.kinds_ends[prev] as usize);
        &self.kinds[start..end as usize]
    }
}

/// The keys of a subject's Seq-align id ([`ScopeWalk`]): its handle, the keys of the ids that
/// NCBI's scope takes for it (`defline_set::scope_keys`) and matches weakly
/// (`id_scope::weak_keys`), and its place in `CSeq_id_Handle`'s order (`order`,
/// `mask_handle_order` of `mask_handle_packed`; handles of the same place but different keys are
/// ordered by their info's address, which LOSAT does not know).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IdKeys {
    handle: Option<HandleKey>,
    scope: Vec<u64>,
    weak: Vec<u64>,
    order: u64,
}

impl IdKeys {
    /// The keys of the Seq-id `id`.
    fn of(id: &ParsedSeqId) -> Self {
        let mut scope = Vec::new();
        crate::seqid::defline_set::scope_keys(id, &mut scope);
        let mut weak = Vec::new();
        crate::seqid::id_scope::weak_keys(id, &mut weak);
        IdKeys {
            handle: handle_key(id),
            scope,
            weak,
            order: crate::seqid::mask_handle::mask_handle_order(
                crate::seqid::mask_handle::mask_handle_packed(id),
            ),
        }
    }
}

/// A sequence that registered a handle in the data loader's `m_Ids` last, and whether LOSAT knows
/// that it was the last (`sure`: not where several sequences of one record of loads registered
/// it, whose order the record does not keep, `ScopeRecord::note`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Registrant {
    sequence: DbSequence,
    sure: bool,
}

impl Registrant {
    /// `sequence`, surely the last.
    fn of(sequence: DbSequence) -> Self {
        Registrant {
            sequence,
            sure: true,
        }
    }
}

// NCBI reference (598d8ae6): c++/include/objtools/data_loaders/blastdb/bdbloader.hpp:155, 222
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:186, 379-423
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/cached_sequence.cpp:98-110
// NCBI reference (598d8ae6): c++/include/util/limited_size_map.hpp:199-206, 233-280
// The BLAST database data loader's `m_Ids` is a map of at most 1,000 handles with least recently
// used eviction: a find of a present handle marks it used, an insert adds a new one or marks a
// present one used without changing it, `RegisterIds`' assignment adds or marks used and assigns;
// each addition then drops the least recently used handles while the map holds more than 1,000
// (`x_GC`). So a handle is held exactly where fewer than 1,000 other handles were added or used
// after its last use (review 1 of line DBID holdrec, M1: `lru12x100`, 1,202 registrations after a
// holder's, is NCBI's own record; `lru9x100`, 902, and `lru12x1`, 12, still the holder's).
/// The loader's `m_Ids` as LOSAT knows it ([`ScopeWalk`], [`ScopeRecord`]): each handle's last
/// registrant ([`Registrant`]) and the times of its last use, a range where LOSAT does not know
/// when it was last used (`lo` the earliest, `hi` the latest; `lo == 0` where LOSAT does not know
/// that the handle was ever added: a load that may not have happened). A handle is held ([`Held`])
/// where fewer than 1,000 others may have been used after its earliest last use, and gone where at
/// least 1,000 others were used after its latest; in between LOSAT does not know. Timestamp counts
/// (`by_lo`, `by_hi`) store one entry per distinct timestamp rather than per handle, and `gc`
/// drops handles that are surely gone.
#[derive(Clone, Debug, Default)]
struct LoaderIds {
    entries: rustc_hash::FxHashMap<HandleKey, IdEntry>,
    /// Count of active entries with each earliest last-use time `lo > 0`.
    by_lo: std::collections::BTreeMap<u64, u32>,
    /// Count of active entries with each latest last-use time `hi > 0`.
    by_hi: std::collections::BTreeMap<u64, u32>,
    /// Count of active entries whose registrant is each sequence ordinal (so callers can drop
    /// `BlobKeys` of evicted sequences in O(1)).
    seq_counts: rustc_hash::FxHashMap<u32, u32>,
    /// For handles with `!registrant.sure`, all sequences that may have registered the handle last
    /// (for O(1) `ScopeRecord::registrants_of` without scanning earlier loads).
    ambiguous: rustc_hash::FxHashMap<HandleKey, Vec<DbSequence>>,
    /// The number of uses so far (each use takes the next time, starting at 1).
    clock: u64,
}

/// One handle of [`LoaderIds`].
#[derive(Clone, Copy, Debug)]
struct IdEntry {
    registrant: Registrant,
    /// `0` means `None` (maybe-only use); `>= 1` is `Some(lo)`.
    lo: u64,
    hi: u64,
}

/// Whether the loader's `m_Ids` holds a handle ([`LoaderIds::held`]): surely, with its
/// registrant; maybe (LOSAT does not know whether 1,000 other handles were used after it, or
/// whether it was added at all), with the registrant it holds if it does; surely not.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Held {
    Present(Registrant),
    Maybe(Registrant),
    Absent,
}

/// The size of the loader's `m_Ids` (bdbloader.cpp:186).
const LOADER_IDS: usize = 1000;

impl LoaderIds {
    /// Increments the count for timestamp `time` in `counts`.
    #[inline]
    fn inc_time(counts: &mut std::collections::BTreeMap<u64, u32>, time: u64) {
        *counts.entry(time).or_insert(0) += 1;
    }

    /// Decrements the count for timestamp `time` in `counts`, removing it at zero.
    #[inline]
    fn dec_time(counts: &mut std::collections::BTreeMap<u64, u32>, time: u64) {
        if let std::collections::btree_map::Entry::Occupied(mut e) = counts.entry(time) {
            let c = e.get_mut();
            *c -= 1;
            if *c == 0 {
                e.remove();
            }
        }
    }

    /// Increments the active handle count for sequence `ordinal`.
    #[inline]
    fn inc_seq(&mut self, ordinal: u32) {
        *self.seq_counts.entry(ordinal).or_insert(0) += 1;
    }

    /// Decrements the active handle count for sequence `ordinal`, removing it at zero.
    #[inline]
    fn dec_seq(&mut self, ordinal: u32) {
        if let std::collections::hash_map::Entry::Occupied(mut e) = self.seq_counts.entry(ordinal) {
            let c = e.get_mut();
            *c -= 1;
            if *c == 0 {
                e.remove();
            }
        }
    }

    /// Whether `m_Ids` still holds at least one entry registered to `ordinal`.
    #[inline]
    fn has_sequence(&self, ordinal: u32) -> bool {
        self.seq_counts.contains_key(&ordinal)
    }

    /// Records `sequence` as a possible last registrant of an unsure `handle`.
    fn push_ambiguous(&mut self, handle: HandleKey, sequence: DbSequence) {
        let list = self.ambiguous.entry(handle).or_default();
        if !list.contains(&sequence) {
            list.push(sequence);
        }
    }

    /// Whether `m_Ids` holds `handle` now.
    fn held(&self, handle: HandleKey) -> Held {
        let Some(entry) = self.entries.get(&handle) else {
            return Held::Absent;
        };
        if entry.lo != 0 {
            // At most `self.clock - entry.lo` handles were used after its earliest last use.
            let present = self.clock - entry.lo < LOADER_IDS as u64 || {
                let self_in_range = usize::from(entry.hi > entry.lo);
                let limit = LOADER_IDS + self_in_range;
                let mut total = 0usize;
                for (_, &cnt) in self.by_hi.range((entry.lo + 1)..) {
                    total += cnt as usize;
                    if total >= limit {
                        break;
                    }
                }
                total.saturating_sub(self_in_range) < LOADER_IDS
            };
            if present {
                return Held::Present(entry.registrant);
            }
        }
        // Invariant: `entry.lo <= entry.hi`, so `handle` itself is never in `by_lo` past `entry.hi`.
        let mut after_hi = 0usize;
        for (_, &cnt) in self.by_lo.range((entry.hi + 1)..) {
            after_hi += cnt as usize;
            if after_hi >= LOADER_IDS {
                return Held::Absent;
            }
        }
        Held::Maybe(entry.registrant)
    }

    /// Moves `handle`'s times to `lo` (`0` for `None`) and `hi`, with the registrant `registrant`.
    fn set(&mut self, handle: HandleKey, registrant: Registrant, lo: u64, hi: u64) {
        let new_ord = registrant.sequence.ordinal();
        if let Some(old) = self.entries.insert(handle, IdEntry { registrant, lo, hi }) {
            if old.lo != 0 {
                Self::dec_time(&mut self.by_lo, old.lo);
            }
            Self::dec_time(&mut self.by_hi, old.hi);
            let old_ord = old.registrant.sequence.ordinal();
            if old_ord != new_ord {
                self.dec_seq(old_ord);
                self.inc_seq(new_ord);
            }
        } else {
            self.inc_seq(new_ord);
        }
        if lo != 0 {
            Self::inc_time(&mut self.by_lo, lo);
        }
        Self::inc_time(&mut self.by_hi, hi);
        if registrant.sure {
            self.ambiguous.remove(&handle);
        }
    }

    /// A use of `handle` now (a find that marks it used, or an insert of a new one) whose
    /// registrant is `registrant`; then amortized eviction.
    fn use_now(&mut self, handle: HandleKey, registrant: Registrant) {
        self.clock += 1;
        self.set(handle, registrant, self.clock, self.clock);
        self.gc(false);
    }

    /// `RegisterIds` of a new blob of `sequence` whose lines hold the handles `handles`, in their
    /// order: each assigned (added, or marked used) now.
    fn register(&mut self, handles: impl IntoIterator<Item = HandleKey>, sequence: DbSequence) {
        for handle in handles {
            self.use_now(handle, Registrant::of(sequence));
        }
    }

    /// A use of `handle` with the registrant `registrant` that may not have happened (a blob that
    /// a failed prefetch may have loaded, one that NCBI may have dropped and loaded again): its
    /// latest last use is now, its registrant is `registrant`'s sequence, not sure where it held
    /// another.
    fn maybe_use(&mut self, handle: HandleKey, registrant: Registrant) {
        self.clock += 1;
        let (lo, old_seq, registrant) = match self.entries.get(&handle) {
            Some(entry) => (
                entry.lo,
                Some(entry.registrant.sequence),
                Registrant {
                    sure: entry.registrant.sure
                        && registrant.sure
                        && entry.registrant.sequence == registrant.sequence,
                    ..registrant
                },
            ),
            None => (0, None, registrant),
        };
        self.set(handle, registrant, lo, self.clock);
        if !registrant.sure {
            if let Some(old_seq) = old_seq {
                self.push_ambiguous(handle, old_seq);
            }
            self.push_ambiguous(handle, registrant.sequence);
        }
        self.gc(false);
    }

    /// Uses of pre-merged `merged` (handles with their registrants, plus optional extra sequences
    /// where a handle had multiple registrants) in an order that LOSAT does not know, advancing
    /// `clock` by `total_uses` (at least `merged.len()`).
    fn use_unordered_map(
        &mut self,
        merged: rustc_hash::FxHashMap<HandleKey, Registrant>,
        extra_ambiguous: rustc_hash::FxHashMap<HandleKey, Vec<DbSequence>>,
        total_uses: u64,
    ) {
        if merged.is_empty() {
            return;
        }
        let start = self.clock + 1;
        self.clock += total_uses.max(merged.len() as u64);
        let end = self.clock;
        for (handle, registrant) in merged {
            self.set(handle, registrant, start, end);
        }
        for (handle, seqs) in extra_ambiguous {
            if self
                .entries
                .get(&handle)
                .is_some_and(|e| !e.registrant.sure)
            {
                for seq in seqs {
                    self.push_ambiguous(handle, seq);
                }
            }
        }
        self.gc(true);
    }

    /// Convenience wrapper over [`Self::use_unordered_map`] for unit tests.
    #[cfg(test)]
    fn use_unordered(&mut self, registered: &[(HandleKey, Registrant)]) {
        let mut merged: rustc_hash::FxHashMap<HandleKey, Registrant> = Default::default();
        let mut extra_ambiguous: rustc_hash::FxHashMap<HandleKey, Vec<DbSequence>> =
            Default::default();
        for &(handle, registrant) in registered {
            merged
                .entry(handle)
                .and_modify(|existing| {
                    if !existing.sure
                        || !registrant.sure
                        || existing.sequence != registrant.sequence
                    {
                        existing.sure = false;
                        let list = extra_ambiguous.entry(handle).or_default();
                        if !list.contains(&existing.sequence) {
                            list.push(existing.sequence);
                        }
                        if !list.contains(&registrant.sequence) {
                            list.push(registrant.sequence);
                        }
                    }
                })
                .or_insert(registrant);
        }
        self.use_unordered_map(merged, extra_ambiguous, registered.len() as u64);
    }

    // NCBI reference (598d8ae6): c++/include/util/limited_size_map.hpp:243-250 (`x_GC`)
    /// Drops the handles that are surely gone: those after whose latest last use at least 1,000
    /// handles were surely used (the 1,000th latest earliest use is after it). When `force` is
    /// false, runs only once at least `max(LOADER_IDS, entries.len() / 4)` handles are evictable
    /// so every use pays O(1) amortized cost; when `force` is true (batch boundary), drops every
    /// surely-evicted handle immediately. Returns whether any entry was removed.
    fn gc(&mut self, force: bool) -> bool {
        if !force && self.entries.len() <= 2 * LOADER_IDS {
            return false;
        }
        let mut seen = 0usize;
        let mut threshold = None;
        for (&lo, &cnt) in self.by_lo.iter().rev() {
            seen += cnt as usize;
            if seen >= LOADER_IDS {
                threshold = Some(lo);
                break;
            }
        }
        let Some(threshold) = threshold else {
            return false;
        };
        let Some((&min_hi, _)) = self.by_hi.first_key_value() else {
            return false;
        };
        if min_hi >= threshold {
            return false;
        }
        let evictable: usize = self
            .by_hi
            .range(..threshold)
            .map(|(_, &cnt)| cnt as usize)
            .sum();
        if evictable == 0 {
            return false;
        }
        if !force && evictable < (self.entries.len() / 4).max(LOADER_IDS) {
            return false;
        }
        let by_lo = &mut self.by_lo;
        let by_hi = &mut self.by_hi;
        let seq_counts = &mut self.seq_counts;
        let ambiguous = &mut self.ambiguous;
        self.entries.retain(|handle, entry| {
            if entry.hi >= threshold {
                return true;
            }
            if entry.lo != 0 {
                Self::dec_time(by_lo, entry.lo);
            }
            Self::dec_time(by_hi, entry.hi);
            let ord = entry.registrant.sequence.ordinal();
            if let std::collections::hash_map::Entry::Occupied(mut e) = seq_counts.entry(ord) {
                let c = e.get_mut();
                *c -= 1;
                if *c == 0 {
                    e.remove();
                }
            }
            ambiguous.remove(handle);
            false
        });
        true
    }
}

/// What a load of a Seq-align id gave in [`ScopeWalk`], as the scope's cache (`m_Seq_idMap`)
/// keeps it: the sequence whose blob it holds, how it came to it and the definition line by which
/// that blob's Bioseq was loaded (its `Bioseq.id` and title, the B-line; [`NO_LINE`] where LOSAT
/// does not hold the blob's lines), or a sequence whose place LOSAT does not know (`Unknown`:
/// another sequence, `None` where LOSAT does not name it, may have come first), or a lookup that
/// throws (never kept).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Resolved {
    Blob(DbSequence, ShownBy, u32),
    Unknown(Option<DbSequence>, ShownBy),
    Throws,
    /// No Bioseq (the loader gave none: the lookup found nothing, or a sequence whose lines do
    /// not hold the id), which the scope keeps as unresolved: a later load of the id in the
    /// scope gets none either, without asking the loader again (`Nothing` is kept in the cache
    /// only; [`ScopeWalk::resolve`] answers it as `None`).
    Nothing,
}

impl Resolved {
    /// This place where LOSAT does not know that NCBI's places reach it (whether `m_Ids` still
    /// holds a handle, whether a failed prefetch loaded a blob): a holder or a weak match answers
    /// [`ShownBy::Order`], so that [`ShownBy::Holder`] and [`ShownBy::Weak`] with a holder name
    /// the sequence that NCBI gives (whose record the writers write, `report::shown_bioseq`; line
    /// DBID hr-3 takes the same rule for psiblast's PSSM).
    fn unsure(self) -> Self {
        match self {
            Self::Blob(sequence, ShownBy::Holder | ShownBy::Weak, line) => {
                Self::Blob(sequence, ShownBy::Order, line)
            }
            Self::Unknown(holder, ShownBy::Holder | ShownBy::Weak) => {
                Self::Unknown(holder, ShownBy::Order)
            }
            other => other,
        }
    }
}

/// How a load came to another sequence than the subject ([`ShownInstead`]): the words of LOSAT's refusal.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum ShownBy {
    /// The lookup itself (step 4), the loader's `m_Ids` by the exact handle (step 3), or the
    /// scope's index by the exact handle or an id that it takes for the Seq-align id (step 2),
    /// where LOSAT knows that NCBI's places reach it.
    Holder,
    /// The scope's index by a weak match (an accession with a version for one without, a name),
    /// where LOSAT knows that NCBI's places reach it.
    Weak,
    /// An order that LOSAT does not know: of a prefetch, two handles that no Seq-id orders (their
    /// infos' addresses); of a record of loads, two sequences of one batch that registered the
    /// same handle; whether `m_Ids` still holds a handle, or a failed prefetch loaded a blob
    /// (`Resolved::unsure`); or whether NCBI keeps its scope for the run (a call whose formats
    /// LOSAT refuses a throw for, `FormatThrow::Refused`).
    Order,
    /// The calls of the run before this one, of which LOSAT has no record
    /// ([`ScopeHolders`]).
    Calls,
    /// The description table's second lookup, by the display id ([`ScopeWalk::display`]).
    Display,
}

/// The line of a [`ShownInstead`] (and of a `Resolved::Blob`) where LOSAT does not know by which
/// definition line the holder's blob was loaded (its keys are not held).
pub const NO_LINE: u32 = u32::MAX;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1305-1470
// The holder's Bioseq (`CSeqDBVol::GetBioseq` with the blob's handle as the target) has the ids
// and title of the target line, every line in its user object, and the holder's length and
// letters (design of line DBID holdrec, `R/hr/logs/ncbi_model.md` 1.2).
/// Another sequence that NCBI gives for a subject's Seq-align id at a load ([`HitLoad::Shown`]):
/// that sequence (`None` where LOSAT does not name it), how it came to it, and the definition
/// line of it by which its Bioseq was loaded (`line`, the B-line: its `Bioseq.id` and title;
/// [`NO_LINE`] where LOSAT does not know it), whose record the writers write for the subject
/// ([`ShownInstead::written`], `report::shown_bioseq`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ShownInstead {
    pub holder: Option<DbSequence>,
    pub by: ShownBy,
    pub line: u32,
}

impl ShownInstead {
    /// A sequence (or none) given by `by`, of a line that LOSAT does not know.
    pub fn unlined(holder: Option<DbSequence>, by: ShownBy) -> Self {
        ShownInstead {
            holder,
            by,
            line: NO_LINE,
        }
    }

    /// The holder and its B-line where LOSAT knows the Bioseq that NCBI's place gives, so that the
    /// record writers write it (`report::shown_bioseq::ShownBioseq`): a holder that the scope's
    /// index, the loader's `m_Ids`, the prefetch's order or the lookup gives (`ShownBy::Holder`,
    /// `ShownBy::Weak`); `None` where LOSAT does not name the holder or its line, or does not know
    /// the order (`Order`, `Calls`), and for the description table's second lookup (`Display`,
    /// a pairwise title).
    pub fn written(self) -> Option<(DbSequence, usize)> {
        match self.by {
            ShownBy::Holder | ShownBy::Weak if self.line != NO_LINE => {
                Some((self.holder?, self.line as usize))
            }
            _ => None,
        }
    }
}

/// Inverted index built lazily for [`ScopePlaces`] only when `holders()` is actually queried
/// (prefetch-off, failed prefetch, or pairwise display lookup): stays unallocated on prefetch-on
/// tabular and XML runs.
#[derive(Debug, Default)]
struct ScopeIndex {
    exact: rustc_hash::FxHashMap<HandleKey, Vec<u32>>,
    taken: rustc_hash::FxHashMap<u64, Vec<u32>>,
    kinds: rustc_hash::FxHashMap<u64, Vec<u32>>,
}

impl ScopeIndex {
    /// Indexes load `at` with its line's `handles`, `scope`, and `kinds`.
    fn add_load(&mut self, at: u32, handles: &[HandleKey], scope: &[u64], kinds: &[u64]) {
        for &handle in handles {
            self.exact.entry(handle).or_default().push(at);
        }
        for &key in scope {
            self.taken.entry(key).or_default().push(at);
        }
        for &key in kinds {
            self.kinds.entry(key).or_default().push(at);
        }
    }

    /// Builds the inverted index from `places`'s flat arenas in one linear pass.
    fn build(places: &ScopePlaces) -> Self {
        let mut index = ScopeIndex::default();
        let mut h_start = 0usize;
        let mut s_start = 0usize;
        let mut k_start = 0usize;
        for at in 0..places.loads.len() {
            let h_end = places.handle_ends[at] as usize;
            let s_end = places.scope_ends[at] as usize;
            let k_end = places.kinds_ends[at] as usize;
            index.add_load(
                at as u32,
                &places.line_handles[h_start..h_end],
                &places.line_scope[s_start..s_end],
                &places.line_kinds[k_start..k_end],
            );
            h_start = h_end;
            s_start = s_end;
            k_start = k_end;
        }
        index
    }
}

/// One scope's places ([`ScopeRecord`], [`ScopeWalk`]): the ids resolved (`m_Seq_idMap`) and
/// the sequences loaded, in the order of their loads (`m_TSE_BySeqId`), with the keys of the one
/// definition line by which each was loaded stored in flat arenas and indexed lazily on demand.
#[derive(Debug, Default)]
struct ScopePlaces {
    cache: rustc_hash::FxHashMap<HandleKey, Resolved>,
    /// The sequences loaded, in order (8 bytes per load).
    loads: Vec<DbSequence>,
    /// The definition line of each load (its B-line, 4 bytes per load).
    lines: Vec<u32>,
    /// Flat arenas of the ONE loaded definition line per entry in `loads`.
    line_handles: Vec<HandleKey>,
    handle_ends: Vec<u32>,
    line_scope: Vec<u64>,
    scope_ends: Vec<u32>,
    line_kinds: Vec<u64>,
    kinds_ends: Vec<u32>,
    /// Inverted index over `loads`, built only on the first `holders()` query.
    index: std::sync::OnceLock<ScopeIndex>,
    /// The ranges of loads from a failed prefetch ([`ScopeWalk::failed_prefetch`]), whose
    /// occurrence and order among themselves LOSAT does not know.
    maybe: Vec<std::ops::Range<u32>>,
}

impl Clone for ScopePlaces {
    fn clone(&self) -> Self {
        ScopePlaces {
            cache: self.cache.clone(),
            loads: self.loads.clone(),
            lines: self.lines.clone(),
            line_handles: self.line_handles.clone(),
            handle_ends: self.handle_ends.clone(),
            line_scope: self.line_scope.clone(),
            scope_ends: self.scope_ends.clone(),
            line_kinds: self.line_kinds.clone(),
            kinds_ends: self.kinds_ends.clone(),
            index: std::sync::OnceLock::new(),
            maybe: self.maybe.clone(),
        }
    }
}

impl ScopePlaces {
    /// Adds a load of `sequence` by its definition line `line`, directly from that line's Seq-ids
    /// `ids` (zero temporary heap allocations).
    fn add_ids(&mut self, sequence: DbSequence, line: usize, ids: &[ParsedSeqId]) {
        let at = u32::try_from(self.loads.len()).unwrap_or(u32::MAX);
        let h_start = self.line_handles.len();
        let s_start = self.line_scope.len();
        let k_start = self.line_kinds.len();
        for id in ids {
            self.line_handles.extend(handle_key(id));
            crate::seqid::defline_set::scope_keys(id, &mut self.line_scope);
            crate::seqid::id_scope::kind_keys(id, &mut self.line_kinds);
        }
        self.handle_ends.push(self.line_handles.len() as u32);
        self.scope_ends.push(self.line_scope.len() as u32);
        self.kinds_ends.push(self.line_kinds.len() as u32);
        self.loads.push(sequence);
        self.lines.push(line_number(line));
        if let Some(index) = self.index.get_mut() {
            index.add_load(
                at,
                &self.line_handles[h_start..],
                &self.line_scope[s_start..],
                &self.line_kinds[k_start..],
            );
        }
    }

    /// Adds a load of `keys` by the line `line`.
    fn add_blob(&mut self, keys: &BlobKeys, line: usize) {
        let Some(sequence) = keys.sequence else {
            return;
        };
        let at = u32::try_from(self.loads.len()).unwrap_or(u32::MAX);
        let handles = keys.line_handles(line);
        let scope = keys.line_scope(line);
        let kinds = keys.line_kinds(line);
        self.line_handles.extend_from_slice(handles);
        self.line_scope.extend_from_slice(scope);
        self.line_kinds.extend_from_slice(kinds);
        self.handle_ends.push(self.line_handles.len() as u32);
        self.scope_ends.push(self.line_scope.len() as u32);
        self.kinds_ends.push(self.line_kinds.len() as u32);
        self.loads.push(sequence);
        self.lines.push(line_number(line));
        if let Some(index) = self.index.get_mut() {
            index.add_load(at, handles, scope, kinds);
        }
    }

    /// The loads whose line holds the id of `keys` exactly (`exact`), else by an id that the scope
    /// takes for it or matches weakly, in the order of the loads, each with how it holds it.
    fn holders(&self, keys: &IdKeys, exact: bool) -> Vec<(u32, ShownBy)> {
        if self.loads.is_empty() {
            return Vec::new();
        }
        let index = self.index.get_or_init(|| ScopeIndex::build(self));
        let mut found: Vec<(u32, ShownBy)> = Vec::new();
        if exact {
            if let Some(at) = keys.handle.and_then(|handle| index.exact.get(&handle)) {
                found.extend(at.iter().map(|&at| (at, ShownBy::Holder)));
            }
        } else {
            for key in &keys.scope {
                if let Some(at) = index.taken.get(key) {
                    found.extend(at.iter().map(|&at| (at, ShownBy::Holder)));
                }
            }
            for key in &keys.weak {
                if let Some(at) = index.kinds.get(key) {
                    found.extend(at.iter().map(|&at| (at, ShownBy::Weak)));
                }
            }
        }
        found.sort_unstable();
        found.dedup_by_key(|(at, _)| *at);
        found
    }

    /// Whether the load `at` comes from a failed prefetch (`maybe`) and may not have happened.
    fn is_maybe(&self, at: u32) -> bool {
        self.maybe.iter().any(|range| range.contains(&at))
    }

    /// The first of the loads `found` (of [`ScopePlaces::holders`]) as the scope's index gives it,
    /// and whether a load in `found` is definite (not only from a failed prefetch, `maybe`):
    /// `Unknown` where another sequence's load could come first (among `maybe` loads before the
    /// first definite load, or against that definite load).
    fn first_of(&self, found: &[(u32, ShownBy)]) -> Option<(Resolved, bool)> {
        let first_definite = found.iter().position(|&(at, _)| !self.is_maybe(at));
        let (candidates, definite) = match first_definite {
            Some(index) => (&found[..=index], true),
            None => (found, false),
        };
        let &(first, by) = candidates.first()?;
        let sequence = self.sequence(first);
        if candidates
            .iter()
            .any(|&(at, _)| self.sequence(at) != sequence)
        {
            return Some((Resolved::Unknown(Some(sequence), ShownBy::Order), definite));
        }
        let line = one_line(candidates.iter().map(|&(at, _)| self.line(at)));
        Some((Resolved::Blob(sequence, by, line), definite))
    }

    /// The sequence of the load `at`.
    fn sequence(&self, at: u32) -> DbSequence {
        self.loads[at as usize]
    }

    /// The definition line by which the load `at` was made.
    fn line(&self, at: u32) -> u32 {
        self.lines[at as usize]
    }
}

/// The B-line number of a definition line index ([`NO_LINE`] past `u32`).
fn line_number(line: usize) -> u32 {
    u32::try_from(line).unwrap_or(NO_LINE)
}

/// The line of the loads of one sequence whose lines are `lines` where they are one line, else
/// [`NO_LINE`]: the scope then holds that sequence's blob by several lines, and which Bioseq it
/// gives is LOSAT's open point (the sequence is the same: the subject itself stays loaded, and no
/// holder's record is written).
fn one_line(mut lines: impl Iterator<Item = u32>) -> u32 {
    let Some(first) = lines.next() else {
        return NO_LINE;
    };
    match lines.all(|line| line == first) {
        true => first,
        false => NO_LINE,
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:438-447 (the blob
// id: the OID and the requested handle)
// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:1999-2070 (cited at `ScopeWalk`)
// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:1515-1560
// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3280-3365 (cited at `ScopeWalk`)
// A blob that the data source holds (locked by a scope, or among the unlocked blobs it keeps) is
// given again without `x_LoadData`, so its `RegisterIds` does not run again: only a new blob id
// registers. The data source keeps up to 10 unlocked blobs (`OBJMGR_BLOB_CACHE`) and a scope up to
// 10 that no handle locks (`OBJMGR_SCOPE_AUTORELEASE_SIZE`); a blob dropped from both is loaded
// again, with `RegisterIds`, where a later read of the scope reaches it (`SUnloadedInfo::LockTSE`,
// `GetBlobById`). Which blobs the formatters' handles keep locked LOSAT does not follow, so in a
// run with more than 10 blob ids each read of a blob loaded before may register its lines again
// (`LoaderIds::maybe_use`); with 10 or fewer, nothing is ever dropped.
/// The number of blobs that NCBI's data source keeps unlocked (`OBJMGR_BLOB_CACHE`'s default):
/// with no more blob ids in the run, no blob is dropped and loaded again.
const BLOBS_KEPT: usize = 10;

/// One recorded scope (`ScopeRecord`'s latest batch, or the whole run when the run keeps one
/// scope): its loads (`places`), Seq-align id resolutions (`by_id`), and unresolved handles
/// (`unresolved`). Queries read only the entries from before a walk's snapshot `(notes, loads)`,
/// without scanning or filtering unrelated batches.
#[derive(Clone, Debug, Default)]
struct RecordedBatchScope {
    batch: usize,
    places: ScopePlaces,
    /// By Seq-align id handle: `(load_index, sequence)` of each load in this scope.
    by_id: rustc_hash::FxHashMap<HandleKey, Vec<(u32, DbSequence)>>,
    /// Handles that resolved to no Bioseq in this scope, with the notes that recorded them.
    unresolved: rustc_hash::FxHashMap<HandleKey, Vec<u32>>,
}

impl RecordedBatchScope {
    /// Step 1 (`m_Seq_idMap`) of `handle` among the first `loads` loads and `notes` notes of this
    /// recorded scope.
    fn cached(&self, handle: HandleKey, notes: u32, loads: u32) -> Option<Resolved> {
        let unresolved = self
            .unresolved
            .get(&handle)
            .is_some_and(|at| at.iter().any(|&note| note < notes));
        let found = self
            .by_id
            .get(&handle)
            .into_iter()
            .flatten()
            .filter(|&&(at, _)| at < loads);
        let Some(&(_, first)) = found.clone().next() else {
            return unresolved.then_some(Resolved::Nothing);
        };
        Some(
            if !unresolved && found.clone().all(|&(_, sequence)| sequence == first) {
                let line = one_line(found.map(|&(at, _)| self.places.line(at)));
                Resolved::Blob(first, ShownBy::Holder, line)
            } else {
                Resolved::Unknown(None, ShownBy::Order)
            },
        )
    }

    /// Step 2 (`m_TSE_BySeqId`) among the first `loads` loads of this recorded scope.
    fn scope_match(&self, keys: &IdKeys, exact: bool, loads: u32) -> Option<Resolved> {
        if loads == 0 {
            return None;
        }
        let mut found: Vec<(u32, ShownBy)> = self.places.holders(keys, exact);
        found.retain(|&(at, _)| at < loads);
        let &(first, by) = found.first()?;
        let sequence = self.places.sequence(first);
        Some(
            if found
                .iter()
                .all(|&(at, _)| self.places.sequence(at) == sequence)
            {
                let line = one_line(found.iter().map(|&(at, _)| self.places.line(at)));
                Resolved::Blob(sequence, by, line)
            } else {
                Resolved::Unknown(Some(sequence), ShownBy::Order)
            },
        )
    }
}

/// What a run loaded before a walk of its reports' loads ([`LoadsBefore`], kept by
/// `common::database::TaxdbRun::note_loaded`, the one producer):
/// - the data loader's `m_Ids` ([`LoaderIds`], copy-on-write via `Arc`);
/// - the compact `BlobKeys` of sequences that `m_Ids` may still hold as registrants (plus kept
///   scope sequences when the run keeps one scope, `keeps_scope != Some(false)`);
/// - the blob ids loaded (ordinal id, requested handle, with the first note of each: a blob loaded
///   again registers nothing while at most [`BLOBS_KEPT`] blob ids exist);
/// - the current batch's scope (`current_batch`: dropped immediately on every new batch when
///   `keeps_scope == Some(false)` so memory and lookup time stay bounded per batch), and the
///   run-wide scope (`all_scope`: populated only when `keeps_scope != Some(false)`, i.e. outfmt 5
///   or unknown, where lookups never filter out earlier batches).
#[derive(Debug, Default)]
pub struct ScopeRecord {
    /// Whether the run keeps one scope across query batches (`Some(false)` for every format that
    /// calls `ResetScopeHistory`, `Some(true)` for outfmt 5, `None` when unknown / tests).
    keeps_scope: Option<bool>,
    ids: std::sync::Arc<LoaderIds>,
    /// Compact `BlobKeys` of sequences that `ids` may still hold (and, when `keeps_scope !=
    /// Some(false)`, sequences of the kept scope for `relock`). Copy-on-write via `Arc` so
    /// `LoadsBefore` reads it in O(1) without locking or scanning.
    loader_blobs: std::sync::Arc<rustc_hash::FxHashMap<u32, std::sync::Arc<BlobKeys>>>,
    /// The blob ids loaded (ordinal id, requested handle; 0 for none), with the first note of
    /// each.
    blob_ids: rustc_hash::FxHashMap<(u32, HandleKey), u32>,
    /// The first note after which the run had loaded more than [`BLOBS_KEPT`] blob ids.
    many_blobs: Option<u32>,
    /// Latest batch's recorded scope (dropped on every new batch when `keeps_scope ==
    /// Some(false)`, and snapshot-shared via `Arc` with `LoadsBefore`).
    current_batch: Option<std::sync::Arc<RecordedBatchScope>>,
    /// Run-wide recorded scope, kept only when `keeps_scope != Some(false)` (outfmt 5 or unknown).
    all_scope: Option<RecordedBatchScope>,
    notes: u32,
    total_loads: u32,
}

impl ScopeRecord {
    /// Records whether the run keeps one scope across batches (`Some(false)` drops earlier
    /// batches' scope places on every new batch; `Some(true)` or `None` keeps run-wide scope).
    pub fn set_keeps_scope(&mut self, keeps: Option<bool>) {
        self.keeps_scope = keeps;
    }

    /// Records the loads `loads` of the query batch `batch` (`None`: the loader's alone, as the
    /// loads before a stop), each by the subject's own Seq-align id: its lookup in `m_Ids` (the
    /// sequence that the notes before registered for it, [`ScopeRecord::held_instead`]), else
    /// the blob that its lookup gives ([`ParsedSubject::loaded_blob_ref`]), whose lines register
    /// every handle where its blob id is new (where it is not, and the run has more than
    /// [`BLOBS_KEPT`] blob ids, they may register again), and its id resolves to it; an id whose
    /// lookup gives no blob (none, or a record whose lines do not hold it) resolves to no Bioseq
    /// where no blob of the run may have registered its handle (where one may have, the order of
    /// the note's loads decides, which LOSAT does not know: a later read asks the loader).
    pub fn note<'p>(
        &mut self,
        batch: Option<usize>,
        subjects: impl IntoIterator<Item = &'p ParsedSubject>,
    ) {
        let mut loads: Vec<(&ParsedSubject, LoadedBlobRef)> = Vec::new();
        let mut unresolved: Vec<HandleKey> = Vec::new();
        // The handles of loads that `m_Ids` may hold for another sequence than the lookup's,
        // with that sequence: LOSAT does not know which the load gives.
        let mut unsure: Vec<(HandleKey, DbSequence)> = Vec::new();
        for subject in subjects {
            match subject.loaded_blob_ref() {
                Some(blob) => match self.held_instead(subject, &blob) {
                    Ok(held) => loads.push((subject, held.unwrap_or(blob))),
                    Err(other) => {
                        unsure.push(other);
                        loads.push((subject, blob));
                    }
                },
                None if !matches!(
                    subject.resolution.found,
                    Some(FoundInstead::Failed(_) | FoundInstead::ThrowsOnId(_))
                ) =>
                {
                    unresolved.extend(subject.id_handle());
                }
                None => {}
            }
        }
        // Neither before the note nor after it may `m_Ids` hold the handle.
        unresolved.retain(|&handle| self.ids.held(handle) == Held::Absent);
        if loads.is_empty() && unresolved.is_empty() {
            return;
        }
        let note = self.notes;
        self.notes += 1;
        let many = self.many_blobs.is_some();
        let mut merged: rustc_hash::FxHashMap<HandleKey, Registrant> = Default::default();
        let mut extra_ambiguous: rustc_hash::FxHashMap<HandleKey, Vec<DbSequence>> =
            Default::default();
        let mut total_uses = 0u64;
        let mut maybe: Vec<(HandleKey, Registrant)> = Vec::new();
        let mut record_use = |handle: HandleKey, registrant: Registrant| {
            total_uses += 1;
            merged
                .entry(handle)
                .and_modify(|existing| {
                    if !existing.sure
                        || !registrant.sure
                        || existing.sequence != registrant.sequence
                    {
                        existing.sure = false;
                        let list = extra_ambiguous.entry(handle).or_default();
                        if !list.contains(&existing.sequence) {
                            list.push(existing.sequence);
                        }
                        if !list.contains(&registrant.sequence) {
                            list.push(registrant.sequence);
                        }
                    }
                })
                .or_insert(registrant);
        };
        for (subject, blob) in &loads {
            let sequence = blob.sequence();
            let requested = subject.id_handle();
            if let Some(handle) = requested {
                record_use(handle, Registrant::of(sequence));
            }
            let id = (sequence.ordinal(), requested.unwrap_or(0));
            let is_new = *self.blob_ids.entry(id).or_insert(note) == note;
            if is_new {
                blob.for_each_handle(subject, |handle| {
                    record_use(handle, Registrant::of(sequence));
                });
            } else if many {
                blob.for_each_handle(subject, |handle| {
                    maybe.push((handle, Registrant::of(sequence)));
                });
            }
        }
        let ids = std::sync::Arc::make_mut(&mut self.ids);
        ids.use_unordered_map(merged, extra_ambiguous, total_uses);
        for (handle, registrant) in maybe {
            ids.maybe_use(handle, registrant);
        }
        for (handle, sequence) in unsure {
            ids.maybe_use(
                handle,
                Registrant {
                    sequence,
                    sure: false,
                },
            );
        }
        ids.gc(true);
        let keeps_all = self.keeps_scope != Some(false);
        // Retain compact `BlobKeys` for sequences that `ids` may still hold (and, where the run
        // keeps one scope, for `relock`).
        let blobs = std::sync::Arc::make_mut(&mut self.loader_blobs);
        if !keeps_all {
            blobs.retain(|&ord, _| ids.has_sequence(ord));
        }
        for (subject, blob) in &loads {
            let ord = blob.sequence().ordinal();
            if keeps_all || ids.has_sequence(ord) {
                blobs
                    .entry(ord)
                    .or_insert_with(|| blob.to_arc_blob_keys(subject));
            }
        }
        if self.many_blobs.is_none() && self.blob_ids.len() > BLOBS_KEPT {
            self.many_blobs = Some(note);
        }
        let Some(batch) = batch else {
            return;
        };
        unresolved.retain(|&handle| self.ids.held(handle) == Held::Absent);
        if loads.is_empty() && unresolved.is_empty() {
            return;
        }
        self.total_loads = self.total_loads.saturating_add(loads.len() as u32);
        // Update `current_batch` (dropping the previous batch whenever `batch` changes), unless
        // `keeps_scope == Some(true)` (where only `all_scope` is ever read).
        if self.keeps_scope != Some(true) {
            if self
                .current_batch
                .as_ref()
                .is_none_or(|current| current.batch != batch)
            {
                self.current_batch = Some(std::sync::Arc::new(RecordedBatchScope {
                    batch,
                    ..RecordedBatchScope::default()
                }));
            }
            let current = std::sync::Arc::make_mut(
                self.current_batch.as_mut().expect("current batch created"),
            );
            for &handle in &unresolved {
                current.unresolved.entry(handle).or_default().push(note);
            }
            for (subject, blob) in &loads {
                let at = u32::try_from(current.places.loads.len()).unwrap_or(u32::MAX);
                blob.add_to_places(&mut current.places, subject);
                if let Some(id) = subject.id_handle() {
                    current
                        .by_id
                        .entry(id)
                        .or_default()
                        .push((at, blob.sequence()));
                }
            }
        }
        // Update `all_scope` only when the run may keep one scope across batches (outfmt 5 or
        // unknown).
        if keeps_all {
            let all = self.all_scope.get_or_insert_with(Default::default);
            for handle in unresolved {
                all.unresolved.entry(handle).or_default().push(note);
            }
            for (subject, blob) in &loads {
                let at = u32::try_from(all.places.loads.len()).unwrap_or(u32::MAX);
                blob.add_to_places(&mut all.places, subject);
                if let Some(id) = subject.id_handle() {
                    all.by_id.entry(id).or_default().push((at, blob.sequence()));
                }
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:379-423
    /// Where `m_Ids` before this note holds the Seq-align id of `subject` for another sequence
    /// than its lookup's blob `blob` (step 3 of [`ScopeWalk::load`], which the lookup follows):
    /// `Ok(Some)` that sequence's blob by its first line that holds the id, where `m_Ids` surely
    /// holds it and LOSAT holds that sequence's keys (a holder that a batch before loaded, which
    /// the reports of the later batches write too: line DBID hr2c2, rows `db.hr2.mb.slH.*`);
    /// `Err` the handle with that sequence where LOSAT does not know whether the load gives it;
    /// `Ok(None)` where the lookup answers.
    fn held_instead(
        &self,
        subject: &ParsedSubject,
        blob: &LoadedBlobRef,
    ) -> Result<Option<LoadedBlobRef>, (HandleKey, DbSequence)> {
        let Some(handle) = subject.id_handle() else {
            return Ok(None);
        };
        let held = self.ids.held(handle);
        let (Held::Present(registrant) | Held::Maybe(registrant)) = held else {
            return Ok(None);
        };
        let sequence = registrant.sequence;
        if sequence == blob.sequence() {
            return Ok(None);
        }
        if !registrant.sure || held != Held::Present(registrant) {
            return Err((handle, sequence));
        }
        if subject.mine(sequence) {
            let line = subject.line_of_handle(handle).unwrap_or(subject.target);
            return Ok(Some(LoadedBlobRef::Own { sequence, line }));
        }
        let keys = self.loader_blobs.get(&sequence.ordinal());
        match keys.and_then(|keys| Some((keys, keys.line_of(handle)?))) {
            Some((keys, line)) => Ok(Some(LoadedBlobRef::External {
                blob: keys.clone(),
                sequence,
                line,
            })),
            None => Err((handle, sequence)),
        }
    }

    /// Whether the blob id `id` was loaded before the note `notes`.
    fn loaded_before(&self, id: (u32, HandleKey), notes: u32) -> bool {
        self.blob_ids.get(&id).is_some_and(|&note| note < notes)
    }

    /// The sequences that `m_Ids` as `ids` holds (a walk's start, [`LoadsBefore`]) may hold for
    /// the handle of `subject`'s Seq-align id: the last registrant where the handle may still be
    /// held, and where LOSAT does not know that it registered it last, every sequence recorded in
    /// `ids.ambiguous` for that handle (O(1), never scans earlier loads).
    fn registrants_of(ids: &LoaderIds, subject: &ParsedSubject) -> Vec<DbSequence> {
        let Some(handle) = subject.id_handle() else {
            return Vec::new();
        };
        let (Held::Present(last) | Held::Maybe(last)) = ids.held(handle) else {
            return Vec::new();
        };
        if last.sure {
            return vec![last.sequence];
        }
        let mut found = ids.ambiguous.get(&handle).cloned().unwrap_or_default();
        found.push(last.sequence);
        found.sort_unstable_by_key(|sequence| sequence.ordinal());
        found.dedup();
        found
    }

    /// The sequences that `m_Ids` may hold now, by ordinal id, in order (tests).
    #[cfg(test)]
    pub fn registered_sequences(&self) -> Vec<u32> {
        let mut ordinals: Vec<u32> = self
            .ids
            .entries
            .keys()
            .filter_map(|&handle| match self.ids.held(handle) {
                Held::Present(registrant) | Held::Maybe(registrant) => {
                    Some(registrant.sequence.ordinal())
                }
                Held::Absent => None,
            })
            .collect();
        ordinals.sort_unstable();
        ordinals.dedup();
        ordinals
    }

    /// The number of notes and loads now: a walk's start ([`LoadsBefore`]).
    fn now(&self) -> (u32, u32) {
        (self.notes, self.total_loads)
    }
}

/// The run's [`ScopeRecord`], shared by the record (`common::database::TaxdbRun`) and the walks
/// that read it as it was at their start.
pub type SharedRecord = std::sync::Arc<std::sync::RwLock<ScopeRecord>>;

/// What a run loaded before the reports that [`HitLoads::of`] walks
/// (`common::database::TaxdbRun::loads_before`, the one producer): the record of the batches,
/// psiblast's rounds and their PSSM reads before ([`ScopeRecord`], as it was: its first `notes`
/// notes and `loads` loads), the query batch of the reports' first report (`batch`, whose scope
/// the record's last batch is where they are the same; `None`: the record's last batch), and
/// whether calls of the run before this one loaded sequences that no record here holds
/// (`unrecorded`: `api::run_batches`' calls after the first, each a query batch of its own).
#[derive(Clone, Debug, Default)]
pub struct LoadsBefore {
    record: SharedRecord,
    notes: u32,
    loads: u32,
    /// Number of loads in `ScopeRecord::all_scope` at this snapshot.
    all_loads: u32,
    batch: Option<usize>,
    unrecorded: bool,
    /// The record's `m_Ids` as it was ([`ScopeRecord::ids`]).
    ids: std::sync::Arc<LoaderIds>,
    /// The record's `loader_blobs` as it was ([`ScopeRecord::loader_blobs`]).
    loader_blobs: std::sync::Arc<rustc_hash::FxHashMap<u32, std::sync::Arc<BlobKeys>>>,
    /// Snapshot of the matching latest batch's scope (when `batch` is `None` or equals the
    /// record's latest batch), with the number of loads in that batch at snapshot time.
    batch_scope: Option<(std::sync::Arc<RecordedBatchScope>, u32)>,
    /// The number of blob ids that the record had loaded.
    blobs: usize,
}

impl LoadsBefore {
    /// Nothing loaded before: the reports are the run's first.
    pub fn none() -> Self {
        LoadsBefore::default()
    }

    /// The loads that `record` holds now, before reports of the query batch `batch`, after calls
    /// whose loads no record holds where `unrecorded`.
    pub fn new(record: &SharedRecord, batch: Option<usize>, unrecorded: bool) -> Self {
        let read = record.read().expect("the run's record");
        let (notes, loads) = read.now();
        let all_loads = read.all_scope.as_ref().map_or(0, |s| {
            u32::try_from(s.places.loads.len()).unwrap_or(u32::MAX)
        });
        let batch_scope = read.current_batch.as_ref().and_then(|current| {
            batch.is_none_or(|b| b == current.batch).then(|| {
                (
                    std::sync::Arc::clone(current),
                    u32::try_from(current.places.loads.len()).unwrap_or(u32::MAX),
                )
            })
        });
        LoadsBefore {
            record: record.clone(),
            notes,
            loads,
            all_loads,
            batch,
            unrecorded,
            ids: read.ids.clone(),
            loader_blobs: read.loader_blobs.clone(),
            batch_scope,
            blobs: read.blob_ids.len(),
        }
    }

    /// No record, after calls whose loads no record holds where `unrecorded`.
    pub fn after_calls(unrecorded: bool) -> Self {
        LoadsBefore {
            unrecorded,
            ..Self::default()
        }
    }

    /// These loads, also after calls whose loads no record holds where `unrecorded`.
    pub fn with_calls_before(self, unrecorded: bool) -> Self {
        LoadsBefore {
            unrecorded: self.unrecorded || unrecorded,
            ..self
        }
    }

    /// The query batch of the reports (`None`: not named).
    pub fn batch(&self) -> Option<usize> {
        self.batch
    }

    /// The sequences that the loader's record held for the handle of `subject`'s Seq-align id
    /// before these reports ([`ScopeRecord::registrants_of`]).
    pub fn registrants_of(&self, subject: &ParsedSubject) -> Vec<DbSequence> {
        ScopeRecord::registrants_of(&self.ids, subject)
    }

    /// Whether the run had loaded the blob id `id` before these reports.
    fn loaded_before(&self, id: (u32, HandleKey)) -> bool {
        self.notes > 0 && self.read(|record| record.loaded_before(id, self.notes))
    }

    /// Whether `other` is the same record as these loads, as it was (no note since).
    fn same_record(&self, other: &LoadsBefore) -> bool {
        // Two records with no note (`TaxdbRun::calls_before` makes one each time) are the same.
        let empty = |loads: &LoadsBefore| (loads.notes, loads.loads) == (0, 0);
        (std::sync::Arc::ptr_eq(&self.record, &other.record) || (empty(self) && empty(other)))
            && (self.notes, self.loads, self.unrecorded)
                == (other.notes, other.loads, other.unrecorded)
    }

    /// Reads the record as it was.
    fn read<T>(&self, read: impl FnOnce(&ScopeRecord) -> T) -> T {
        read(&self.record.read().expect("the run's record"))
    }
}

/// Reference to the blob that a load of a [`ParsedSubject`] loads (step 4 of
/// [`ScopeWalk::load`]), without allocating a `BlobKeys` when the blob is the subject's own
/// `DeflineSet`.
#[derive(Clone, Debug)]
enum LoadedBlobRef {
    Own {
        sequence: DbSequence,
        line: usize,
    },
    External {
        blob: std::sync::Arc<BlobKeys>,
        sequence: DbSequence,
        line: usize,
    },
}

impl LoadedBlobRef {
    #[inline]
    fn sequence(&self) -> DbSequence {
        match self {
            LoadedBlobRef::Own { sequence, .. } | LoadedBlobRef::External { sequence, .. } => {
                *sequence
            }
        }
    }

    #[inline]
    fn line(&self) -> usize {
        match self {
            LoadedBlobRef::Own { line, .. } | LoadedBlobRef::External { line, .. } => *line,
        }
    }

    fn for_each_handle(&self, subject: &ParsedSubject, mut f: impl FnMut(HandleKey)) {
        match self {
            LoadedBlobRef::Own { .. } => {
                for handle in subject.line_handles() {
                    f(handle);
                }
            }
            LoadedBlobRef::External { blob, .. } => {
                for &handle in blob.all_handles() {
                    f(handle);
                }
            }
        }
    }

    fn add_to_places(&self, places: &mut ScopePlaces, subject: &ParsedSubject) {
        match self {
            LoadedBlobRef::Own { sequence, line } => {
                places.add_ids(*sequence, *line, &subject.set.deflines[*line].ids);
            }
            LoadedBlobRef::External { blob, line, .. } => {
                places.add_blob(blob, *line);
            }
        }
    }

    fn to_arc_blob_keys(&self, subject: &ParsedSubject) -> std::sync::Arc<BlobKeys> {
        match self {
            LoadedBlobRef::Own { sequence, .. } => subject.blob_keys_as(*sequence),
            LoadedBlobRef::External { blob, .. } => blob.clone(),
        }
    }
}

/// How a single tie item in [`ScopeWalk::bulk_ties`] registers its blob's handles.
#[derive(Clone, Copy, PartialEq, Eq)]
enum TieRegister {
    None,
    Definite,
    Maybe,
}

/// What resolving one tie item alone against an immutable `&ScopeWalk` would load.
struct AloneBlobLoad {
    sequence: DbSequence,
    line: usize,
    new_blob_id: Option<(u32, HandleKey)>,
    register: TieRegister,
    /// `None` when the blob's handles and loaded line come directly from the item's `ParsedSubject`.
    external: Option<std::sync::Arc<BlobKeys>>,
    /// `Some(handles)` only for a `relock` on a foreign sequence whose handles came from
    /// `scope_handles` (no `BlobKeys` needed, `add_to_scope` is false).
    relock_handles: Option<Box<[HandleKey]>>,
    add_to_scope: bool,
}

/// Outcome of resolving one tie item alone against an immutable `&ScopeWalk`.
struct AloneResult {
    resolved: Option<Resolved>,
    /// `(handle, registrant, maybe)` for the item's own handle use before `load_blob`.
    requested_use: Option<(HandleKey, Registrant, bool)>,
    blob_load: Option<AloneBlobLoad>,
}

// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:1769-1810, 2195-2230
// NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:832-939
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:229-286, 379-462
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/cached_sequence.cpp:98-110
// NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3280-3365 (the prefetch's bulk load)
// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:1999-2070
// NCBI's formatters find a subject again by its Seq-align id at each load, in four places, the
// first that answers winning (design of line DBID holdrec, `R/hr/logs/ncbi_model.md` 2.2): (1) the
// scope's cache of the ids it resolved (`m_Seq_idMap`; a lookup that threw is not kept); (2) the
// scope's index of the Bioseqs it loaded (`m_TSE_BySeqId`), by the ids of the one definition line
// each was loaded by: the id's exact handle first, then the handles that match it, weak matches
// included (`GetMatchingHandles` with `eAllowWeakMatch`), of which the one loaded first
// (`x_IsBetter`, P3.a of that probe: not the lower OID); (3) the BLAST database data loader's
// `m_Ids`, the exact handle of every Seq-id of every definition line of each new blob it loaded
// (`RegisterIds`, which overwrites: the last registered wins, P3.b), kept for the run up to its
// 1,000 most recently used handles ([`LoaderIds`]; a blob id loaded before registers nothing
// again, [`BLOBS_KEPT`]); (4) the
// lookup (`SeqidToOid`). The scope (1, 2) is emptied after each query batch (`ResetScopeHistory`),
// but for outfmt 5, and psiblast's is the query's. A batch whose prefetch runs resolves all its
// Seq-align ids at once before its reports, in `CSeq_id_Handle`'s order (a gi's value, the digits
// of an accession that packs, then the infos' addresses), each by (1), then (3) with the
// registrations of the ids before it in that pass, then (4); the scope's index plays no part
// there (the data source's answer replaces it), and the reports then find the pass's answers in
// the cache (review 4 of line DBID pfoff-4, M1: 100 prefetch-on cases).
/// NCBI's places that answer the loads of a run's reports ([`HitLoads::of`]), in NCBI's order of
/// the loads: the one owner of "which sequence does a load by a Seq-align id give". It starts
/// from the run's record ([`LoadsBefore`]) and adds its own loads; each load is answered by
/// [`ScopeWalk::load`] (a report that reads one by one) or [`ScopeWalk::bulk`] (a batch's
/// prefetch). Where LOSAT does not know whether the run keeps one scope (`None`), it walks both
/// ways (`other`) and answers `ShownBy::Order` where they differ.
#[derive(Debug, Default)]
pub struct ScopeWalk {
    before: LoadsBefore,
    /// The loader's `m_Ids`: the record's as it was, then the walk's own uses.
    ids: LoaderIds,
    /// The blob ids (ordinal id, requested handle; 0 for none) that the walk loaded, and those
    /// that a failed prefetch may have loaded (`maybe_blobs`): a blob loaded before registers
    /// nothing again ([`ScopeWalk::load_blob`]).
    blob_ids: rustc_hash::FxHashSet<(u32, HandleKey)>,
    maybe_blobs: rustc_hash::FxHashSet<(u32, HandleKey)>,
    /// The walk's own places of the current scope.
    scope: ScopePlaces,
    /// Compact `BlobKeys` of sequences that `self.ids` may still hold (pruned when evicted from
    /// `self.ids`), so `held_blob` and `first_blob` can load a holder by any of its lines.
    blobs: rustc_hash::FxHashMap<u32, std::sync::Arc<BlobKeys>>,
    /// Flat handles of sequences loaded in the current scope (cleared on `enter` when
    /// `!keeps_scope`), used only if `relock` needs a foreign holder whose handles were already
    /// evicted from `self.ids`.
    scope_handles: rustc_hash::FxHashMap<u32, Box<[HandleKey]>>,
    /// The query batch of the current scope (`None`: none named yet).
    batch: Option<usize>,
    /// Whether the current scope began in the record (`before.batch`'s scope).
    in_record: bool,
    /// Whether the run keeps one scope (outfmt 5, `RunThrow::keeps_scope`).
    keeps_scope: bool,
    /// The walk that keeps one scope where this one does not, where LOSAT does not know which.
    other: Option<Box<ScopeWalk>>,
}

impl ScopeWalk {
    /// A walk after the loads `before`, in a run that keeps one scope where `keeps_scope` (`None`:
    /// LOSAT does not know).
    pub fn new(before: LoadsBefore, keeps_scope: Option<bool>) -> Self {
        let other = keeps_scope
            .is_none()
            .then(|| Box::new(ScopeWalk::new(before.clone(), Some(true))));
        let batch = before.batch;
        ScopeWalk {
            ids: (*before.ids).clone(),
            before,
            blob_ids: Default::default(),
            maybe_blobs: Default::default(),
            scope: ScopePlaces::default(),
            blobs: Default::default(),
            scope_handles: Default::default(),
            batch,
            in_record: true,
            keeps_scope: keeps_scope.unwrap_or(false),
            other,
        }
    }

    /// The scope of a report of the query batch `batch` (`None`: of the batch of the report
    /// before): where it is another batch, NCBI's formatter emptied the scope at the end of that
    /// one (`ResetScopeHistory`), but where the run keeps one scope.
    pub fn enter(&mut self, batch: Option<usize>) {
        if let Some(other) = &mut self.other {
            other.enter(batch);
        }
        let Some(batch) = batch else {
            return;
        };
        if self.batch.is_some_and(|current| current != batch) && !self.keeps_scope {
            self.scope = ScopePlaces::default();
            self.scope_handles.clear();
            self.in_record = false;
            self.prune_blobs(false);
        }
        self.batch = Some(batch);
    }

    /// Drops cached `BlobKeys` for sequences no longer held by `self.ids`. When `save_handles` is
    /// true (still inside a live scope), preserves their handle lists in `scope_handles` so
    /// [`Self::relock`] stays exact for sequences still in the scope.
    fn prune_blobs(&mut self, save_handles: bool) {
        let ids = &self.ids;
        let scope_handles = &mut self.scope_handles;
        self.blobs.retain(|&ord, keys| {
            if ids.has_sequence(ord) {
                return true;
            }
            if save_handles {
                scope_handles
                    .entry(ord)
                    .or_insert_with(|| keys.all_handles().to_vec().into_boxed_slice());
            }
            false
        });
    }

    /// Step 1: the cache of the current scope.
    fn cached(&self, handle: HandleKey) -> Option<Resolved> {
        if let Some(&resolved) = self.scope.cache.get(&handle) {
            return Some(resolved);
        }
        if self.keeps_scope {
            if self.before.all_loads == 0 && self.before.notes == 0 {
                return None;
            }
            let notes = self.before.notes;
            let all_loads = self.before.all_loads;
            return self.before.read(|record| {
                record
                    .all_scope
                    .as_ref()
                    .and_then(|scope| scope.cached(handle, notes, all_loads))
            });
        }
        if !self.in_record {
            return None;
        }
        let (batch_scope, loads) = self.before.batch_scope.as_ref()?;
        if self.batch.is_some_and(|b| b != batch_scope.batch) {
            return None;
        }
        batch_scope.cached(handle, self.before.notes, *loads)
    }

    /// Whether the run has loaded more than [`BLOBS_KEPT`] blob ids, so that NCBI may have dropped
    /// one and loads it again where a read reaches it.
    fn many_blobs(&self) -> bool {
        self.before.blobs + self.blob_ids.len() + self.maybe_blobs.len() > BLOBS_KEPT
    }

    /// Whether the run loaded the blob id `id` before (the record's, or the walk's).
    fn loaded_before(&self, id: (u32, HandleKey)) -> bool {
        self.blob_ids.contains(&id) || self.before.loaded_before(id)
    }

    /// Step 3 as a place, without a load: the registrant of `handle` where `m_Ids` may hold it,
    /// `Unknown` where LOSAT does not know that it registered it last.
    fn registered_place(&self, handle: HandleKey) -> Option<Resolved> {
        match self.ids.held(handle) {
            Held::Present(registrant) | Held::Maybe(registrant) => {
                Some(self.place(registrant, handle))
            }
            Held::Absent => None,
        }
    }

    /// The place that `m_Ids`' registrant `registrant` gives for `handle`: its blob by that
    /// handle, whose Bioseq's line is the first of its lines that holds the handle
    /// (`CSeqDBVol::GetBioseq`'s target; [`NO_LINE`] where LOSAT does not hold its keys).
    fn place(&self, registrant: Registrant, handle: HandleKey) -> Resolved {
        match registrant.sure {
            true => Resolved::Blob(
                registrant.sequence,
                ShownBy::Holder,
                self.registered_line(registrant.sequence, handle),
            ),
            false => Resolved::Unknown(Some(registrant.sequence), ShownBy::Order),
        }
    }

    /// The first line of `sequence` that holds `handle`, where LOSAT holds its keys.
    fn registered_line(&self, sequence: DbSequence, handle: HandleKey) -> u32 {
        self.blob_of(sequence)
            .and_then(|keys| keys.line_of(handle))
            .map_or(NO_LINE, line_number)
    }

    /// Step 2: the scope's index, for the Seq-align id of `keys`: the exact handle first, then
    /// the ids that it takes for it or matches weakly; of each, the sequence loaded first (the
    /// record's before the walk's; among the record's, which LOSAT holds as sets, one where they
    /// are the same sequence, else `Unknown`), and whether a matching load is definite (not only
    /// from a failed prefetch, where a miss can fall through to a weak match or the loader).
    fn scope_match(&self, keys: &IdKeys) -> Option<(Resolved, bool)> {
        let exact = self.scope_match_pass(keys, true);
        if let Some((resolved, true)) = exact {
            return Some((resolved, true));
        }
        let weak = self.scope_match_pass(keys, false);
        match (exact, weak) {
            (None, weak) => weak,
            (Some(exact), None) => Some(exact),
            (
                Some((Resolved::Blob(a, by, line), false)),
                Some((Resolved::Blob(b, _, other_line), definite)),
            ) if a == b => Some((
                Resolved::Blob(a, by, one_line([line, other_line].into_iter())),
                definite,
            )),
            (Some((exact, false)), Some((weak, _))) => {
                let holder = match weak {
                    Resolved::Blob(sequence, ..) => Some(sequence),
                    Resolved::Unknown(holder, _) => holder,
                    Resolved::Throws | Resolved::Nothing => None,
                };
                let holder = holder.or(match exact {
                    Resolved::Blob(sequence, ..) => Some(sequence),
                    Resolved::Unknown(holder, _) => holder,
                    Resolved::Throws | Resolved::Nothing => None,
                });
                // Which of the two NCBI's scope gives depends on what the failed prefetch loaded.
                Some((Resolved::Unknown(holder, ShownBy::Order), true))
            }
            (Some((_, true)), _) => unreachable!("definite exact match returned above"),
        }
    }

    /// One pass of [`ScopeWalk::scope_match`] (`exact` or weak/taken), returning whether a
    /// matching load is definite.
    fn scope_match_pass(&self, keys: &IdKeys, exact: bool) -> Option<(Resolved, bool)> {
        let recorded = if self.keeps_scope {
            let all_loads = self.before.all_loads;
            (all_loads > 0).then(|| {
                self.before.read(|record| {
                    record
                        .all_scope
                        .as_ref()
                        .and_then(|scope| scope.scope_match(keys, exact, all_loads))
                })
            })
        } else if self.in_record {
            self.before
                .batch_scope
                .as_ref()
                .filter(|(scope, _)| self.batch.is_none_or(|b| b == scope.batch))
                .map(|(scope, loads)| scope.scope_match(keys, exact, *loads))
        } else {
            None
        }
        .flatten();
        if let Some(found) = recorded {
            return Some((found, true));
        }
        self.scope.first_of(&self.scope.holders(keys, exact))
    }

    /// Loads `subject`'s own blob (as `sequence`) by `line` and `requested`: registers its handles
    /// directly from `subject.set.deflines` without allocating a `BlobKeys` unless `self.ids`
    /// retains `sequence`.
    fn load_subject_blob(
        &mut self,
        subject: &ParsedSubject,
        sequence: DbSequence,
        line: usize,
        requested: Option<HandleKey>,
    ) {
        let ord = sequence.ordinal();
        let id = (ord, requested.unwrap_or(0));
        if self.loaded_before(id) {
            if self.many_blobs() {
                for handle in subject.line_handles() {
                    self.ids.maybe_use(handle, Registrant::of(sequence));
                }
            }
        } else if self.maybe_blobs.remove(&id) {
            for handle in subject.line_handles() {
                self.ids.maybe_use(handle, Registrant::of(sequence));
            }
            self.blob_ids.insert(id);
        } else {
            self.ids.register(subject.line_handles(), sequence);
            self.blob_ids.insert(id);
        }
        self.scope
            .add_ids(sequence, line, &subject.set.deflines[line].ids);
        if self.ids.has_sequence(ord) {
            self.blobs
                .entry(ord)
                .or_insert_with(|| subject.blob_keys_as(sequence));
        } else if !self.blobs.contains_key(&ord) {
            self.scope_handles.entry(ord).or_insert_with(|| {
                let handles: Vec<HandleKey> = subject.line_handles().collect();
                handles.into_boxed_slice()
            });
        }
        if self.blobs.len() > self.ids.seq_counts.len() * 2 + 64 {
            self.prune_blobs(true);
        }
    }

    /// The blob of `keys` loaded by the requested handle `requested` (its blob id): where it is
    /// new, `RegisterIds` (every handle of every line, overwriting); where the run loaded it
    /// before, nothing, unless NCBI may have dropped it (more than [`BLOBS_KEPT`] blob ids) or a
    /// failed prefetch may have loaded it, where its lines may register again
    /// ([`LoaderIds::maybe_use`]); then its Bioseq enters the scope by its line `line`.
    fn load_blob(
        &mut self,
        keys: &std::sync::Arc<BlobKeys>,
        line: usize,
        requested: Option<HandleKey>,
    ) {
        let Some(sequence) = keys.sequence else {
            return;
        };
        let ord = sequence.ordinal();
        let id = (ord, requested.unwrap_or(0));
        if self.loaded_before(id) {
            if self.many_blobs() {
                self.maybe_register(keys);
            }
        } else if self.maybe_blobs.remove(&id) {
            self.maybe_register(keys);
            self.blob_ids.insert(id);
        } else {
            self.ids
                .register(keys.all_handles().iter().copied(), sequence);
            self.blob_ids.insert(id);
        }
        self.scope.add_blob(keys, line);
        if self.ids.has_sequence(ord) {
            self.blobs.entry(ord).or_insert_with(|| keys.clone());
        } else {
            self.scope_handles
                .entry(ord)
                .or_insert_with(|| keys.all_handles().to_vec().into_boxed_slice());
        }
        if self.blobs.len() > self.ids.seq_counts.len() * 2 + 64 {
            self.prune_blobs(true);
        }
    }

    /// A `RegisterIds` of the blob of `keys` that may run (a blob that NCBI may have dropped and
    /// loads again, or that a failed prefetch may have loaded).
    fn maybe_register(&mut self, keys: &BlobKeys) {
        let Some(sequence) = keys.sequence else {
            return;
        };
        for &handle in keys.all_handles() {
            self.ids.maybe_use(handle, Registrant::of(sequence));
        }
    }

    /// Where NCBI may have dropped the blob of `sequence` (more than [`BLOBS_KEPT`] blob ids in
    /// the run), a read that the scope answers with it may load it again, which registers its
    /// lines ([`LoaderIds::maybe_use`]). When `subject.mine(sequence)`, reads `subject` directly
    /// with zero lookup or allocation.
    fn relock(&mut self, subject: &ParsedSubject, resolved: Option<Resolved>) {
        let Some(Resolved::Blob(sequence, ..) | Resolved::Unknown(Some(sequence), _)) = resolved
        else {
            return;
        };
        if !self.many_blobs() {
            return;
        }
        if subject.mine(sequence) {
            for handle in subject.line_handles() {
                self.ids.maybe_use(handle, Registrant::of(sequence));
            }
            return;
        }
        if let Some(keys) = self.blob_of(sequence) {
            self.maybe_register(&keys);
        } else if let Some(handles) = self.scope_handles.get(&sequence.ordinal()).cloned() {
            for &handle in &handles {
                self.ids.maybe_use(handle, Registrant::of(sequence));
            }
        }
    }

    /// Steps 3 and 4 of a load of `subject` (also the prefetch's, which skips step 2): the
    /// sequence that the loader registered for the handle last (where that is another one, its
    /// new blob registers its lines again), else in a call after the first a holder that the calls
    /// before may have registered ([`ScopeHolders::loader`]), else the lookup
    /// ([`ScopeResolution`]): its blob, which registers its lines, a failure (`Throws`), or
    /// nothing (no blob: the lookup gives none, or another record, whose refusal is the
    /// lookup's, `FoundInstead`).
    fn loader(&mut self, subject: &ParsedSubject, handle: Option<HandleKey>) -> Option<Resolved> {
        if let Some(handle) = handle {
            let registrant = match self.ids.held(handle) {
                Held::Present(registrant) => Some(registrant),
                // Where LOSAT does not know whether `m_Ids` still holds the handle, the lookup
                // that would answer instead gives the same sequence, or LOSAT does not know.
                Held::Maybe(registrant)
                    if registrant.sure && self.looked_up(subject) == Some(registrant.sequence) =>
                {
                    Some(registrant)
                }
                Held::Maybe(registrant) => {
                    self.ids.maybe_use(handle, registrant);
                    return Some(match subject.mine(registrant.sequence) {
                        true => Resolved::Unknown(Some(registrant.sequence), ShownBy::Order),
                        false => self.place(registrant, handle).unsure(),
                    });
                }
                Held::Absent => None,
            };
            if let Some(registrant) = registrant {
                // The find marks the handle used.
                self.ids.use_now(handle, registrant);
                if !registrant.sure {
                    return Some(Resolved::Unknown(Some(registrant.sequence), ShownBy::Order));
                }
                if subject.mine(registrant.sequence) {
                    // The subject's own blob, loaded before in the run: its Bioseq enters this
                    // scope by the line of the handle, where a later load finds it (review 4 of
                    // line DBID pfoff-4, `pq.zr1.V.U1.b1.0`: a description's load, then a failing
                    // lookup's alignment that this Bioseq holds weakly); it registers its lines
                    // only where its blob id is new (`load_subject_blob`).
                    let line = subject.line_of_handle(handle).unwrap_or(subject.target);
                    self.load_subject_blob(subject, registrant.sequence, line, Some(handle));
                    return Some(Resolved::Blob(
                        registrant.sequence,
                        ShownBy::Holder,
                        line_number(line),
                    ));
                }
                let line = self.held_blob(registrant.sequence, handle);
                return Some(Resolved::Blob(registrant.sequence, ShownBy::Holder, line));
            }
        }
        if self.before.unrecorded {
            if let Some(holder) = subject.resolution.unrecorded.loader.first() {
                return Some(Resolved::Unknown(holder, ShownBy::Calls));
            }
        }
        if matches!(
            subject.resolution.found,
            Some(FoundInstead::Failed(_) | FoundInstead::ThrowsOnId(_))
        ) {
            return Some(Resolved::Throws);
        }
        let blob = subject.loaded_blob_ref()?;
        let (sequence, line) = (blob.sequence(), line_number(blob.line()));
        // The lookup adds the handle with the sequence it gives.
        if let Some(handle) = handle {
            self.ids.use_now(handle, Registrant::of(sequence));
        }
        match blob {
            LoadedBlobRef::Own { sequence, line } => {
                self.load_subject_blob(subject, sequence, line, handle);
            }
            LoadedBlobRef::External { blob, line, .. } => {
                self.load_blob(&blob, line, handle);
            }
        }
        Some(Resolved::Blob(sequence, ShownBy::Holder, line))
    }

    /// The sequence whose blob the lookup of `subject`'s Seq-align id gives ([`ScopeResolution`]),
    /// where LOSAT knows it: `None` for none, a lookup that throws, and a call after the first
    /// whose calls before may hold another ([`ScopeHolders::loader`]).
    fn looked_up(&self, subject: &ParsedSubject) -> Option<DbSequence> {
        if self.before.unrecorded && subject.resolution.unrecorded.loader.first().is_some() {
            return None;
        }
        if matches!(
            subject.resolution.found,
            Some(FoundInstead::Failed(_) | FoundInstead::ThrowsOnId(_))
        ) {
            return None;
        }
        subject.loaded_blob_ref().map(|blob| blob.sequence())
    }

    /// Where step 3 gave another sequence, `holder`, for the handle `handle`: its new blob (by
    /// that handle) registers its lines again and puts its Bioseq into the scope by the first
    /// line that holds the handle, where LOSAT holds its keys (the walk's or the record's): that
    /// line ([`NO_LINE`] where LOSAT does not hold them).
    fn held_blob(&mut self, holder: DbSequence, handle: HandleKey) -> u32 {
        let Some(keys) = self.blob_of(holder) else {
            return NO_LINE;
        };
        let Some(line) = keys.line_of(handle) else {
            return NO_LINE;
        };
        self.load_blob(&keys, line, Some(handle));
        line_number(line)
    }

    /// Resolves the Seq-align id of `subject` at a load by a report that reads it (steps 1 to 4,
    /// [`ScopeWalk`]), records what the load holds, and gives the answer as the reports take it:
    /// `HitLoad::Loaded` where it is the subject or a copy of its record that the lookup gives
    /// (or no blob, where the lookup's answer refuses it), `HitLoad::Throws` where the lookup
    /// throws, `HitLoad::Shown` where it is another sequence.
    pub fn load(&mut self, subject: &ParsedSubject) -> HitLoad {
        let resolved = self.resolve(subject, true);
        let answer = Self::answer(subject, resolved);
        match &mut self.other {
            Some(other) => Self::either(answer, other.load(subject)),
            None => answer,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:432-433, 497-498
    /// The description table's second lookup of a described subject whose load gave its own
    /// record (outfmt 0 to 4: the title of the Bioseq of its display id,
    /// [`ParsedSubject::display_keys`]), after its load: the cache, then the scope's index (the
    /// Bioseq loaded first whose line holds it: the subject's own holds it, so the index answers
    /// at the latest there), then the loader's `m_Ids`; kept in the cache. `Some` where that is
    /// another sequence (`ShownBy::Display`; probe `R/hr/bin/hr1open.py` D.ab: one batch whose
    /// prefetch is off finds the first record's title for the second, the prefetch the first
    /// loaded by the bulk order for both).
    pub fn display(&mut self, subject: &ParsedSubject) -> Option<ShownInstead> {
        let one = self.display_one(subject);
        let Some(other) = &mut self.other else {
            return one;
        };
        let other = other.display_one(subject);
        if one == other {
            return one;
        }
        Some(ShownInstead::unlined(
            one.or(other).and_then(|shown| shown.holder),
            ShownBy::Display,
        ))
    }

    /// [`ScopeWalk::display`] of this walk alone.
    fn display_one(&mut self, subject: &ParsedSubject) -> Option<ShownInstead> {
        let keys = subject.display_keys()?;
        let handle = keys.handle?;
        let resolved = self
            .cached(handle)
            .or_else(|| self.scope_match(&keys).map(|(resolved, _)| resolved))
            .or_else(|| self.registered_place(handle))?;
        self.scope.cache.insert(handle, resolved);
        match Self::answer(subject, Some(resolved)) {
            HitLoad::Shown(shown) => Some(ShownInstead {
                by: ShownBy::Display,
                ..shown
            }),
            _ => None,
        }
    }

    /// Whether the answers `one` and `other` give the same sequence the same way, whatever the
    /// line of a holder's Bioseq ([`ScopeWalk::bulk_ties`]).
    fn same_sequence(one: HitLoad, other: HitLoad) -> bool {
        match (one, other) {
            (HitLoad::Shown(one), HitLoad::Shown(other)) => {
                (one.holder, one.by) == (other.holder, other.by)
            }
            (one, other) => one == other,
        }
    }

    /// The answer of two walks that differ in whether the run keeps one scope.
    fn either(one: HitLoad, other: HitLoad) -> HitLoad {
        if one == other {
            return one;
        }
        let holder = [one, other].into_iter().find_map(|load| match load {
            HitLoad::Shown(shown) => shown.holder,
            _ => None,
        });
        HitLoad::Shown(ShownInstead::unlined(holder, ShownBy::Order))
    }

    /// [`ScopeWalk::load`]'s resolution; `scope`: whether the scope's index answers (a report's
    /// load, not the prefetch's).
    fn resolve(&mut self, subject: &ParsedSubject, scope: bool) -> Option<Resolved> {
        let handle = subject.id_handle();
        if let Some(handle) = handle {
            if let Some(resolved) = self.cached(handle) {
                if resolved == Resolved::Nothing {
                    return None;
                }
                self.relock(subject, Some(resolved));
                return Some(resolved);
            }
        }
        let resolved = if scope {
            let keys = subject.id_keys();
            match self.scope_match(&keys) {
                Some((resolved, true)) => {
                    self.relock(subject, Some(resolved));
                    Some(resolved)
                }
                Some((resolved, false)) => {
                    // Only a failed prefetch's load (`ScopePlaces::maybe`) matched in the scope:
                    // where that blob was loaded before the failure, the scope gives `resolved`;
                    // where it was not, the read misses the scope and reaches the loader.
                    match resolved {
                        Resolved::Blob(sequence, ..) if subject.mine(sequence) => {
                            match self.resolve_after_scope(subject, handle, true) {
                                Some(Resolved::Blob(other, ..)) if other == sequence => {
                                    Some(resolved)
                                }
                                other => {
                                    let (holder, by) = match other {
                                        Some(Resolved::Blob(holder, by, _)) => (Some(holder), by),
                                        Some(Resolved::Unknown(holder, by)) => (holder, by),
                                        _ => (Some(sequence), ShownBy::Order),
                                    };
                                    Some(Resolved::Unknown(holder, by).unsure())
                                }
                            }
                        }
                        other => {
                            let other = other.unsure();
                            self.relock(subject, Some(other));
                            Some(other)
                        }
                    }
                }
                None => self.resolve_after_scope(subject, handle, true),
            }
        } else {
            self.resolve_after_scope(subject, handle, false)
        };
        // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:3300-3330
        // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:2200-2218
        // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:1115-1125
        // NCBI reference (598d8ae6): c++/include/objmgr/impl/scope_info.hpp:569-572
        // A Seq-id that the loader gives no Bioseq for is kept unresolved in the scope's cache,
        // and only a new data source (none in a BLAST run) makes the scope resolve it again
        // (review 1 of line DBID holdrec, `ahA1`: A's lookup gives Y, whose lines do not hold
        // A's id, before Z registers that id; the PSSM's read of A then gets nothing, and NCBI
        // warns).
        if let Some(handle) = handle {
            if resolved != Some(Resolved::Throws) {
                self.scope
                    .cache
                    .insert(handle, resolved.unwrap_or(Resolved::Nothing));
            }
        }
        resolved
    }

    /// Steps 3 and 4 of [`ScopeWalk::resolve`] after the scope's index misses (or only had a
    /// failed prefetch's `maybe` load): unrecorded-call scope holders for outfmt 5, then the
    /// loader.
    fn resolve_after_scope(
        &mut self,
        subject: &ParsedSubject,
        handle: Option<HandleKey>,
        scope: bool,
    ) -> Option<Resolved> {
        // A scope kept for the run (outfmt 5) after calls whose loads LOSAT has no record
        // of may hold a holder that LOSAT lists from the database (`ScopeHolders::scope`).
        if scope && self.keeps_scope && self.before.unrecorded {
            if let Some(holder) = subject.resolution.unrecorded.scope.first() {
                return Some(Resolved::Unknown(holder, ShownBy::Calls));
            }
        }
        self.loader(subject, handle)
    }

    /// The reports' answer of `resolved` for `subject` ([`ScopeWalk::load`]).
    fn answer(subject: &ParsedSubject, resolved: Option<Resolved>) -> HitLoad {
        match resolved {
            None | Some(Resolved::Nothing) => HitLoad::Loaded,
            Some(Resolved::Throws) => HitLoad::Throws,
            Some(Resolved::Blob(sequence, ..)) if subject.mine(sequence) => HitLoad::Loaded,
            Some(Resolved::Blob(sequence, by, line)) => HitLoad::Shown(ShownInstead {
                holder: Some(sequence),
                by,
                line,
            }),
            Some(Resolved::Unknown(holder, by)) => HitLoad::Shown(ShownInstead::unlined(
                holder.filter(|&holder| !subject.mine(holder)),
                by,
            )),
        }
    }

    /// The prefetch of a batch whose reported subjects are `subjects` (key and parsed record, in
    /// any order): each Seq-align id resolved once, in `CSeq_id_Handle`'s order ([`IdKeys`]), by
    /// the cache, then the loader with the registrations of the ids before it, then the lookup
    /// ([`ScopeWalk`]); where two different handles have the same place in that order and the
    /// order decides what one of them gives (another's blob registers its handle),
    /// `ShownBy::Order`. The answers by key, as [`ScopeWalk::load`] gives them; the reports find
    /// them in the cache.
    pub fn bulk<'p>(
        &mut self,
        subjects: impl IntoIterator<Item = (usize, &'p ParsedSubject)>,
    ) -> std::collections::HashMap<usize, HitLoad> {
        let mut keyed: Vec<(u64, Option<HandleKey>, usize, &ParsedSubject)> = subjects
            .into_iter()
            .map(|(key, subject)| {
                let (handle, order) = subject.id_handle_order();
                (order, handle, key, subject)
            })
            .collect();
        keyed.sort_unstable_by_key(|&(order, _, key, _)| (order, key));
        self.bulk_sorted(&keyed)
    }

    /// Non-generic core of [`ScopeWalk::bulk`] over pre-keyed, sorted subjects.
    fn bulk_sorted(
        &mut self,
        keyed: &[(u64, Option<HandleKey>, usize, &ParsedSubject)],
    ) -> std::collections::HashMap<usize, HitLoad> {
        let mut answers = std::collections::HashMap::with_capacity(keyed.len());
        let mut at = 0;
        while at < keyed.len() {
            let order = keyed[at].0;
            let end = at
                + keyed[at..]
                    .iter()
                    .take_while(|&&(o, _, _, _)| o == order)
                    .count();
            let group = &keyed[at..end];
            let first = group[0].1;
            if first.is_some() && group.iter().all(|&(_, h, _, _)| h == first) {
                for &(_, _, key, subject) in group {
                    let resolved = self.resolve(subject, false);
                    answers.insert(key, Self::answer(subject, resolved));
                }
            } else {
                self.bulk_ties(group, &mut answers);
            }
            at = end;
        }
        if let Some(other) = &mut self.other {
            let others = other.bulk_sorted(keyed);
            for (key, answer) in &mut answers {
                if let Some(&other) = others.get(key) {
                    *answer = Self::either(*answer, other);
                }
            }
        }
        answers
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:756-788
    // NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:1999-2070
    // NCBI's prefetch loads the batch's blobs in `CSeq_id_Handle`'s order until the lookup that
    // throws on the id, whose exception ends the bulk pass (caught with its warning): the blobs
    // loaded before it have registered their lines; which ones those are LOSAT does not follow
    // (the pass's pieces of 100 and 200 ids, line DBID pfoff-7).
    /// The prefetch of a batch that fails at a subject whose lookup throws on the id
    /// (`ParsedSubject::prefetch_fails`): each other subject of its reports, `subjects`, may have
    /// been loaded before the failure, in an order that LOSAT does not follow, so each blob that
    /// its load would give (the loader's registrant of its handle, else its lookup's) may have
    /// registered its lines (a handle that such a blob and another sequence may hold last is not
    /// sure afterwards) and entered the scope by the line of that load, where their order is not
    /// known either (`ScopePlaces::maybe`). No id is resolved: the reports load each one by one.
    pub fn failed_prefetch<'p>(
        &mut self,
        subjects: impl IntoIterator<Item = (usize, &'p ParsedSubject)> + Clone,
    ) {
        if let Some(other) = &mut self.other {
            other.failed_prefetch(subjects.clone());
        }
        let start = u32::try_from(self.scope.loads.len()).unwrap_or(u32::MAX);
        let mut additions: Vec<(&ParsedSubject, LoadedBlobRef)> = Vec::new();
        for (_, subject) in subjects {
            let Some(blob) = self.first_blob(subject) else {
                continue;
            };
            let sequence = blob.sequence();
            let requested = subject.id_handle();
            if let Some(handle) = requested {
                self.ids.maybe_use(handle, Registrant::of(sequence));
            }
            let id = (sequence.ordinal(), requested.unwrap_or(0));
            if !self.loaded_before(id) {
                blob.for_each_handle(subject, |handle| {
                    self.ids.maybe_use(handle, Registrant::of(sequence));
                });
                self.maybe_blobs.insert(id);
            } else if self.many_blobs() {
                blob.for_each_handle(subject, |handle| {
                    self.ids.maybe_use(handle, Registrant::of(sequence));
                });
            }
            let ord = sequence.ordinal();
            if self.ids.has_sequence(ord) {
                self.blobs
                    .entry(ord)
                    .or_insert_with(|| blob.to_arc_blob_keys(subject));
            }
            additions.push((subject, blob));
        }
        for (subject, blob) in additions {
            blob.add_to_places(&mut self.scope, subject);
        }
        let end = u32::try_from(self.scope.loads.len()).unwrap_or(u32::MAX);
        if end > start {
            self.scope.maybe.push(start..end);
        }
    }

    /// What the prefetch's resolution of the Seq-align id of `subject` gives now without loading
    /// anything (steps 1 and 3, [`ScopeWalk::bulk`]'s order of places): the scope's cache, else the
    /// loader's registrant of its handle, else nothing that LOSAT knows (`HitLoad::Loaded`; the
    /// lookup's own answer). For a failed prefetch ([`HitLoads::of`]): where a blob loaded before
    /// the failure registered the id of a subject whose lookup throws, the pass gave that blob.
    pub fn peek(&self, subject: &ParsedSubject) -> HitLoad {
        let resolved = subject.id_handle().and_then(|handle| {
            self.cached(handle)
                .or_else(|| self.registered_place(handle))
        });
        let answer = Self::answer(subject, resolved);
        match &self.other {
            Some(other) => Self::either(answer, other.peek(subject)),
            None => answer,
        }
    }

    /// The blob that a load of `subject` by the loader would give without the scope's index
    /// ([`ScopeWalk::loader`] without its records): none where the scope's cache resolved its
    /// Seq-align id, the blob of the sequence that the loader registered for its handle last (by
    /// the first line that holds it), else its lookup's.
    fn first_blob(&self, subject: &ParsedSubject) -> Option<LoadedBlobRef> {
        if let Some(handle) = subject.id_handle() {
            if self.cached(handle).is_some() {
                return None;
            }
            if let Held::Present(registrant) | Held::Maybe(registrant) = self.ids.held(handle) {
                if subject.mine(registrant.sequence) {
                    let line = subject.line_of_handle(handle)?;
                    return Some(LoadedBlobRef::Own {
                        sequence: registrant.sequence,
                        line,
                    });
                }
                let blob = self.blob_of(registrant.sequence)?;
                let line = blob.line_of(handle)?;
                return Some(LoadedBlobRef::External {
                    sequence: registrant.sequence,
                    blob,
                    line,
                });
            }
        }
        match subject.resolution.found {
            Some(FoundInstead::Failed(_) | FoundInstead::ThrowsOnId(_)) => None,
            _ => subject.loaded_blob_ref(),
        }
    }

    /// The keys of the sequence `sequence` where LOSAT holds them (the walk's blobs, or the
    /// record's snapshot `loader_blobs`, in O(1) without locking or scanning).
    fn blob_of(&self, sequence: DbSequence) -> Option<std::sync::Arc<BlobKeys>> {
        let ord = sequence.ordinal();
        self.blobs
            .get(&ord)
            .cloned()
            .or_else(|| self.before.loader_blobs.get(&ord).cloned())
    }

    /// Plans what `load_blob` / `load_subject_blob` would do for a single tie item against the
    /// immutable walk state `&self`.
    fn plan_alone_blob(
        &self,
        sequence: DbSequence,
        line: usize,
        requested: Option<HandleKey>,
        external: Option<std::sync::Arc<BlobKeys>>,
    ) -> AloneBlobLoad {
        let id = (sequence.ordinal(), requested.unwrap_or(0));
        let (register, new_blob_id) = if self.loaded_before(id) {
            (
                if self.many_blobs() {
                    TieRegister::Maybe
                } else {
                    TieRegister::None
                },
                None,
            )
        } else if self.maybe_blobs.contains(&id) {
            (TieRegister::Maybe, Some(id))
        } else {
            (TieRegister::Definite, Some(id))
        };
        AloneBlobLoad {
            sequence,
            line,
            new_blob_id,
            register,
            external,
            relock_handles: None,
            add_to_scope: true,
        }
    }

    /// Resolves one tie item of [`ScopeWalk::bulk_ties`] as if it came first against `&self`,
    /// without cloning `self.ids`, `self.scope`, `self.blobs`, or `self.blob_ids`.
    fn resolve_alone_bulk(
        &self,
        subject: &ParsedSubject,
        handle: Option<HandleKey>,
    ) -> AloneResult {
        if let Some(handle) = handle {
            if let Some(resolved) = self.cached(handle) {
                if resolved == Resolved::Nothing {
                    return AloneResult {
                        resolved: None,
                        requested_use: None,
                        blob_load: None,
                    };
                }
                let blob_load = if self.many_blobs() {
                    match resolved {
                        Resolved::Blob(seq, ..) | Resolved::Unknown(Some(seq), _) => {
                            if subject.mine(seq) {
                                Some(AloneBlobLoad {
                                    sequence: seq,
                                    line: 0,
                                    new_blob_id: None,
                                    register: TieRegister::Maybe,
                                    external: None,
                                    relock_handles: None,
                                    add_to_scope: false,
                                })
                            } else if let Some(ext) = self.blob_of(seq) {
                                Some(AloneBlobLoad {
                                    sequence: seq,
                                    line: 0,
                                    new_blob_id: None,
                                    register: TieRegister::Maybe,
                                    external: Some(ext),
                                    relock_handles: None,
                                    add_to_scope: false,
                                })
                            } else {
                                self.scope_handles.get(&seq.ordinal()).map(|handles| {
                                    AloneBlobLoad {
                                        sequence: seq,
                                        line: 0,
                                        new_blob_id: None,
                                        register: TieRegister::Maybe,
                                        external: None,
                                        relock_handles: Some(handles.clone()),
                                        add_to_scope: false,
                                    }
                                })
                            }
                        }
                        _ => None,
                    }
                } else {
                    None
                };
                return AloneResult {
                    resolved: Some(resolved),
                    requested_use: None,
                    blob_load,
                };
            }
        }
        // Steps 3 and 4 (`scope = false`): check `self.ids`, then `unrecorded`, then lookup.
        if let Some(handle) = handle {
            let registrant = match self.ids.held(handle) {
                Held::Present(registrant) => Some(registrant),
                Held::Maybe(registrant)
                    if registrant.sure && self.looked_up(subject) == Some(registrant.sequence) =>
                {
                    Some(registrant)
                }
                Held::Maybe(registrant) => {
                    return AloneResult {
                        resolved: Some(match subject.mine(registrant.sequence) {
                            true => Resolved::Unknown(Some(registrant.sequence), ShownBy::Order),
                            false => self.place(registrant, handle).unsure(),
                        }),
                        requested_use: Some((handle, registrant, true)),
                        blob_load: None,
                    };
                }
                Held::Absent => None,
            };
            if let Some(registrant) = registrant {
                if !registrant.sure {
                    return AloneResult {
                        resolved: Some(Resolved::Unknown(
                            Some(registrant.sequence),
                            ShownBy::Order,
                        )),
                        requested_use: Some((handle, registrant, false)),
                        blob_load: None,
                    };
                }
                let blob_load = if subject.mine(registrant.sequence) {
                    let line = subject.line_of_handle(handle).unwrap_or(subject.target);
                    Some(self.plan_alone_blob(registrant.sequence, line, Some(handle), None))
                } else {
                    self.blob_of(registrant.sequence).and_then(|ext| {
                        let line = ext.line_of(handle)?;
                        Some(self.plan_alone_blob(
                            registrant.sequence,
                            line,
                            Some(handle),
                            Some(ext),
                        ))
                    })
                };
                let line = blob_load
                    .as_ref()
                    .map_or(NO_LINE, |load| line_number(load.line));
                return AloneResult {
                    resolved: Some(Resolved::Blob(registrant.sequence, ShownBy::Holder, line)),
                    requested_use: Some((handle, registrant, false)),
                    blob_load,
                };
            }
        }
        if self.before.unrecorded {
            if let Some(holder) = subject.resolution.unrecorded.loader.first() {
                return AloneResult {
                    resolved: Some(Resolved::Unknown(holder, ShownBy::Calls)),
                    requested_use: None,
                    blob_load: None,
                };
            }
        }
        if matches!(
            subject.resolution.found,
            Some(FoundInstead::Failed(_) | FoundInstead::ThrowsOnId(_))
        ) {
            return AloneResult {
                resolved: Some(Resolved::Throws),
                requested_use: None,
                blob_load: None,
            };
        }
        let Some(blob) = subject.loaded_blob_ref() else {
            return AloneResult {
                resolved: None,
                requested_use: None,
                blob_load: None,
            };
        };
        let sequence = blob.sequence();
        let line = blob.line();
        let external = match blob {
            LoadedBlobRef::Own { .. } => None,
            LoadedBlobRef::External { blob, .. } => Some(blob),
        };
        AloneResult {
            resolved: Some(Resolved::Blob(sequence, ShownBy::Holder, line_number(line))),
            requested_use: handle.map(|h| (h, Registrant::of(sequence), false)),
            blob_load: Some(self.plan_alone_blob(sequence, line, handle, external)),
        }
    }

    /// [`ScopeWalk::bulk`] of ids with one place in `CSeq_id_Handle`'s order but different handles
    /// (or none), whose order LOSAT does not know: each resolved as if it came first against the
    /// unmodified walk; where another's blob registers its handle with another answer, the order
    /// decides (`ShownBy::Order`); a handle that two of their blobs register is not sure
    /// afterwards. Runs in O(total handles in `group`) time and O(`group.len()`) space with no
    /// per-item clones of the walk state.
    fn bulk_ties(
        &mut self,
        group: &[(u64, Option<HandleKey>, usize, &ParsedSubject)],
        answers: &mut std::collections::HashMap<usize, HitLoad>,
    ) {
        // Target handles that are not yet cached in this scope: only uses touching one of these
        // can move another tie item's answer.
        let mut target_handles: rustc_hash::FxHashSet<HandleKey> =
            rustc_hash::FxHashSet::with_capacity_and_hasher(group.len(), Default::default());
        for &(_, handle, _, _) in group {
            if let Some(handle) = handle {
                if self.cached(handle).is_none() {
                    target_handles.insert(handle);
                }
            }
        }

        let mut results: Vec<AloneResult> = Vec::with_capacity(group.len());
        let mut merged_uses: rustc_hash::FxHashMap<HandleKey, Registrant> = Default::default();
        let mut extra_ambiguous: rustc_hash::FxHashMap<HandleKey, Vec<DbSequence>> =
            Default::default();
        let mut total_uses = 0u64;
        // For each target handle, at most two `(index, sequence)` entries per distinct sequence
        // (preserving ascending `index` order), so Pass 2 checks every `other != index` in O(1).
        let mut target_hits: rustc_hash::FxHashMap<HandleKey, Vec<(usize, DbSequence)>> =
            Default::default();

        for (index, &(_, handle, _, subject)) in group.iter().enumerate() {
            let alone = self.resolve_alone_bulk(subject, handle);
            let definite_req = alone
                .requested_use
                .and_then(|(h, _, maybe)| (!maybe).then_some(h));
            let mut record_use = |h: HandleKey, reg: Registrant, maybe: bool| {
                let effective = if maybe && definite_req != Some(h) {
                    match self.ids.entries.get(&h) {
                        Some(entry) => Registrant {
                            sure: entry.registrant.sure
                                && reg.sure
                                && entry.registrant.sequence == reg.sequence,
                            ..reg
                        },
                        None => reg,
                    }
                } else {
                    reg
                };
                total_uses += 1;
                merged_uses
                    .entry(h)
                    .and_modify(|existing| {
                        if !existing.sure
                            || !effective.sure
                            || existing.sequence != effective.sequence
                        {
                            existing.sure = false;
                            let list = extra_ambiguous.entry(h).or_default();
                            if !list.contains(&existing.sequence) {
                                list.push(existing.sequence);
                            }
                            if !list.contains(&effective.sequence) {
                                list.push(effective.sequence);
                            }
                        }
                    })
                    .or_insert(effective);
                if target_handles.contains(&h) {
                    let hits = target_hits.entry(h).or_default();
                    let seq = effective.sequence;
                    if !hits.iter().any(|&(i, s)| i == index && s == seq)
                        && hits.iter().filter(|&&(_, s)| s == seq).count() < 2
                    {
                        hits.push((index, seq));
                    }
                }
            };
            if let Some((h, reg, maybe)) = alone.requested_use {
                record_use(h, reg, maybe);
            }
            if let Some(load) = &alone.blob_load {
                if load.register != TieRegister::None {
                    let is_maybe = load.register == TieRegister::Maybe;
                    let reg = Registrant::of(load.sequence);
                    if let Some(ext) = &load.external {
                        for &h in ext.all_handles() {
                            record_use(h, reg, is_maybe);
                        }
                    } else if let Some(handles) = &load.relock_handles {
                        for &h in handles {
                            record_use(h, reg, is_maybe);
                        }
                    } else {
                        for h in subject.line_handles() {
                            record_use(h, reg, is_maybe);
                        }
                    }
                }
            }
            results.push(alone);
        }

        // Pass 2: compute answers and `moved` flags in O(group.len()).
        let mut moved = vec![false; group.len()];
        for (index, &(_, handle, key, subject)) in group.iter().enumerate() {
            let resolved = results[index].resolved;
            let mut answer = Self::answer(subject, resolved);
            if let Some(handle) = handle {
                if target_handles.contains(&handle) {
                    if let Some(hits) = target_hits.get(&handle) {
                        for &(other, seq) in hits {
                            if other == index {
                                continue;
                            }
                            // The line of that blob's Bioseq does not move the answer: the same
                            // holder by its first line that holds the handle.
                            let later = Self::answer(
                                subject,
                                Some(Resolved::Blob(seq, ShownBy::Holder, NO_LINE)),
                            );
                            if !Self::same_sequence(later, answer) {
                                answer = Self::either(answer, later);
                                moved[index] = true;
                            }
                        }
                    }
                }
            }
            answers.insert(key, answer);
        }

        // Pass 3: commit cache, scope loads, blob ids, `m_Ids`, and active `blobs`.
        self.ids
            .use_unordered_map(merged_uses, extra_ambiguous, total_uses);
        for (index, &(_, handle, _, subject)) in group.iter().enumerate() {
            let alone = &results[index];
            // A tie that resolved to nothing alone, where another may have registered its
            // handle before it, is not known to be unresolved: a later read asks the loader.
            let kept = match alone.resolved {
                Some(Resolved::Throws) => None,
                None if moved[index] => None,
                _ => Some(alone.resolved.unwrap_or(Resolved::Nothing)),
            };
            if let (Some(handle), Some(kept)) = (handle, kept) {
                self.scope.cache.entry(handle).or_insert(kept);
            }
            if let Some(load) = &alone.blob_load {
                if let Some(id) = load.new_blob_id {
                    self.maybe_blobs.remove(&id);
                    self.blob_ids.insert(id);
                }
                if load.add_to_scope {
                    if let Some(ext) = &load.external {
                        self.scope.add_blob(ext, load.line);
                    } else {
                        self.scope.add_ids(
                            load.sequence,
                            load.line,
                            &subject.set.deflines[load.line].ids,
                        );
                    }
                }
                let ord = load.sequence.ordinal();
                if self.ids.has_sequence(ord) {
                    self.blobs
                        .entry(ord)
                        .or_insert_with(|| match &load.external {
                            Some(ext) => ext.clone(),
                            None => subject.blob_keys_as(load.sequence),
                        });
                } else if load.add_to_scope && !self.blobs.contains_key(&ord) {
                    self.scope_handles
                        .entry(ord)
                        .or_insert_with(|| match &load.external {
                            Some(ext) => ext.all_handles().to_vec().into_boxed_slice(),
                            None => subject
                                .line_handles()
                                .collect::<Vec<_>>()
                                .into_boxed_slice(),
                        });
                }
            }
        }
        self.prune_blobs(true);
    }
}

/// What NCBI's prefetch of a query batch reads for a subject that no report shows, by its
/// Seq-align id (`ParsedSubject::prefetched`): psiblast's traceback hit lists beyond
/// `-max_target_seqs`, and the Seq-annots with `PRE_FETCH_SEQS_LIMIT=2147483647`
/// (`common::database::prefetch_reads`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Prefetched {
    /// The sequence itself, or the same record whose lookups of its taxids write the same texts
    /// (`ParsedSubject::fetched` without the taxonomy formats).
    Itself,
    /// No sequence: the prefetch reads nothing for it.
    Nothing,
    /// Another sequence (`DbSequence`) where LOSAT knows it, whose taxids the prefetch can look up
    /// instead: `FoundInstead::Sequence` (that sequence where one of its ids compares equal to the
    /// Seq-align id, else none), and the same record whose lookups write other texts (`None`);
    /// `reason` as `fetched` gives it. Another sequence that NCBI's places give for the id
    /// without the lookup is [`ScopeWalk::bulk`]'s answer, which the callers ask.
    Another {
        other: Option<DbSequence>,
        reason: String,
    },
    /// A lookup that fails (the reason): on the id itself NCBI's prefetch catches the exception
    /// and writes the warning `Error pre-fetching sequence data`
    /// (`ParsedSubject::prefetch_fails`, `common::database::batch_prefetch`); on a file it
    /// rethrows the `CSeqDBException` and stops.
    Failed(String),
}

/// A lookup of a Seq-align id that fails on a file of the database (`FoundInstead::Failed`): the
/// message, and whether it is NCBI's own (`SeqDbError::ncbi_message`: the text of the
/// `CSeqDBException` that NCBI writes after `BLAST Database error: `, as "DB contains no accession
/// info." of a list that puts version 5 databases with and without the accession table side by
/// side), else LOSAT's text of an error NCBI does not have.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LookupFailure {
    message: Vec<u8>,
    ncbi: NcbiFailure,
}

/// What NCBI's lookup does where LOSAT's fails on a file ([`LookupFailure`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum NcbiFailure {
    /// It throws a `CSeqDBException` with the message (NCBI's own text).
    Message,
    /// It throws, with a message that LOSAT does not reproduce.
    Other,
    /// It crashes on the file (`SeqDbError::ncbi_crashes`).
    Crashes,
}

impl LookupFailure {
    /// What NCBI does where it reads a subject whose lookup fails on a file: its prefetch rethrows
    /// the `CSeqDBException` (`BLAST Database error`, exit status 2); where `PRE_FETCH_SEQS_LIMIT`
    /// turns the prefetch off, its formatters catch it and write `Unknown`, then stop or not as
    /// the format decides (DECISIONS 53 (1) of line DBID).
    pub const NCBI_READS: &'static str = "NCBI BLAST+ 2.17.0's prefetch stops with a BLAST Database error, and where the prefetch is off its formatters fail to read it";

    /// What NCBI does where it reads a subject whose lookup throws NCBI's exception at a place
    /// LOSAT knows ([`LookupFailure::exception`]): LOSAT writes each format's failure as NCBI does
    /// (`FormatThrow::of`), and refuses it only for the reason that the run gives
    /// ([`ThrowRefusal::words`], added by [`LookupFailure::what_ncbi_does`]).
    pub const NCBI_THROWS: &'static str = "NCBI BLAST+ 2.17.0's prefetch stops with a BLAST Database error, and where the prefetch is off its formatters fail to read it and stop or go on as each format does";

    /// What NCBI does where its lookup reads a file that it crashes on (an ISAM whose header does
    /// not fit its files, as an empty `.psi` or `.psd`: SIGSEGV; review 4 of line DBID loadord):
    /// a run whose lookups never reach it (the sequence was loaded before by another id) goes on.
    pub const NCBI_CRASHES: &'static str =
        "NCBI BLAST+ 2.17.0 crashes where its lookup reads that file";

    /// NCBI's own message.
    pub fn ncbi(message: impl Into<Vec<u8>>) -> Self {
        Self {
            message: message.into(),
            ncbi: NcbiFailure::Message,
        }
    }

    /// LOSAT's text of an error that NCBI does not have.
    pub fn losat(message: impl Into<String>) -> Self {
        Self {
            message: message.into().into_bytes(),
            ncbi: NcbiFailure::Other,
        }
    }

    /// The failure of a lookup with the error `error`: NCBI's message where it is one of NCBI's
    /// (`SeqDbError::ncbi_message`), NCBI's crash where LOSAT refuses a file that NCBI crashes on
    /// (`SeqDbError::ncbi_crashes`) or where a read fails after the file's first byte
    /// (`SeqDbError::Read`: NCBI reads bytes of its mapping that the system no longer gives,
    /// SIGBUS; line WQA D18, 2026-11-21), else LOSAT's text. A first byte that cannot be read
    /// where the file is opened (`SeqDbError::Map`) is NCBI's failed mapping, an error that NCBI
    /// reports: without NCBI's text for it (a code outside `utils::os_error_text` on WASI,
    /// Windows) it is LOSAT's text, as a file that does not open (review r2 of D18, L1).
    pub fn of_error(error: &crate::seqdb::SeqDbError) -> Self {
        match error.ncbi_message() {
            Some(message) => Self::ncbi(message),
            None if error.ncbi_crashes()
                || matches!(error, crate::seqdb::SeqDbError::Read { .. }) =>
            {
                Self {
                    message: error.to_string().into_bytes(),
                    ncbi: NcbiFailure::Crashes,
                }
            }
            None => Self::losat(error.to_string()),
        }
    }

    /// The text NCBI writes after `BLAST Database error: `, where the message is NCBI's own.
    pub fn ncbi_message(&self) -> Option<&[u8]> {
        (self.ncbi == NcbiFailure::Message).then_some(self.message.as_slice())
    }

    /// The `CSeqDBException` that NCBI's lookup throws, where the message is NCBI's own
    /// (`ncbi_message`) and is one that NCBI's lookup of a Seq-align id was seen to throw, with
    /// its place (`seqdb::lookup_throw_site`: "DB contains no accession info." of
    /// `CBlastLMDBManager::CBlastEnv::GetDbi`, at a list that puts version 5 databases with and
    /// without the accession table side by side); `None` for every other text (also the other
    /// texts of NCBI's database reader whose places LOSAT knows, `seqdb::throw_site`, which no
    /// lookup was seen to throw), where only the prefetch's report of it (`BLAST Database error:
    /// <message>`) is known.
    pub fn exception(&self) -> Option<crate::ncbi_exception::NcbiException> {
        let message = self.ncbi_message()?;
        let site = crate::seqdb::lookup_throw_site(message)?;
        Some(crate::ncbi_exception::NcbiException::new(site, message))
    }

    /// What NCBI does where it reads the subject ([`LookupFailure::NCBI_READS`],
    /// [`LookupFailure::NCBI_THROWS`] where its exception is known,
    /// [`LookupFailure::NCBI_CRASHES`]), for the reason of LOSAT's refusal: where its exception is
    /// known and the run's formats refuse it (`refused`, [`FormatThrow::Refused`]), with the
    /// words of that situation ([`ThrowRefusal::words`]), the one text of each situation.
    pub fn what_ncbi_does(&self, refused: Option<ThrowRefusal>) -> std::borrow::Cow<'static, str> {
        match (self.ncbi, refused) {
            (NcbiFailure::Message, Some(refused)) if self.exception().is_some() => {
                format!("{}; {}", Self::NCBI_THROWS, refused.words()).into()
            }
            (NcbiFailure::Message, None) if self.exception().is_some() => Self::NCBI_THROWS.into(),
            (NcbiFailure::Message | NcbiFailure::Other, _) => Self::NCBI_READS.into(),
            (NcbiFailure::Crashes, _) => Self::NCBI_CRASHES.into(),
        }
    }
}

impl std::fmt::Display for LookupFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&String::from_utf8_lossy(&self.message))
    }
}

/// What NCBI does where it reads a subject whose lookup throws on the id itself
/// (`FoundInstead::ThrowsOnId`), for the reason of LOSAT's refusal: its prefetch writes its warning
/// and its formatters meet the exception where they load the subject, each format as where a
/// lookup fails on a file with the prefetch off ([`LookupFailure::what_ncbi_does`]'s situation);
/// where the run's formats refuse it (`refused`, [`FormatThrow::Refused`]), with the words of that
/// situation ([`ThrowRefusal::words`]).
fn throws_on_id_does(refused: Option<ThrowRefusal>) -> std::borrow::Cow<'static, str> {
    const DOES: &str = "NCBI BLAST+ 2.17.0's prefetch writes Error pre-fetching sequence data, and its formatters fail to read it and stop or go on as each format does";
    match refused {
        Some(refused) => format!("{DOES}; {}", refused.words()).into(),
        None => DOES.into(),
    }
}

/// Why the reports can show apart a copy of a subject with the same definition lines but for
/// their taxids (`FoundInstead::CopyApart`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CopyApart {
    /// A title can be made from the taxids (`ParsedSubject::title_reads_source`), which only
    /// the reports that show a generated title show apart
    /// (`ParsedSubject::generated_title_fetched`).
    Title,
    /// The taxonomy database writes a text at a lookup that the volume of the copy makes for
    /// its taxids.
    TaxonomyTexts,
}

/// One `HitDescr` of the XML2 and JSON formats (`s_SeqAlignToXMLHit`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HitDescription {
    pub id: Vec<u8>,
    /// Not written when empty.
    pub accession: Vec<u8>,
    /// Not written when the definition line has no title member.
    pub title: Option<Vec<u8>>,
    /// Not written when 0.
    pub taxid: i32,
    /// What the taxonomy database said about `taxid` (`Unnamed` where nothing was asked); a
    /// `Named` scientific name is written after the taxid.
    pub taxonomy: HitTaxonomy,
}

impl ParsedSubject {
    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:390-395
    // NCBI's formatters find the subject's sequence again from its Seq-align id, by the lookup that
    // `found_instead` ports, whatever kind of id it is: in version 5 `bbs|123`, `nat|X.1|` and
    // `lcl|777` by their LMDB keys, in version 4 the ISAM files that makeblastdb writes for neither
    // `bbs|` nor `nat|`, a `gnl|PIG|` id in the PIG files that it fills from no FASTA id, and a
    // `gnl|BL_ORD_ID|n` as the stored n plus the first OID of its volume
    // (`SeqDb::defline_set`), looked up from the start of each volume, so as its own sequence where
    // the stored n is its own number in its volume, also beside another id (probe of line DBID N2b,
    // `blastdb-ids-20261021/R/unfetch/p1`: 317 cases of NCBI 2.17.0, databases of NCBI's and
    // LOSAT's makeblastdb). The lookup does not decide one case: where the definition line of the
    // Bioseq is found by the content of a local id alone (`lcl|777` before the definition line of
    // the Seq-align id `gi|777`, `DeflineSet::bioseq_defline`), NCBI 2.17.0 finds no Bioseq, in
    // version 4 and 5 (outfmt 6 stops with a CCoreException, XML2 writes the Seq-align id alone;
    // probe of the ID2a review). LOSAT refuses that subject instead of repeating a defect of
    // NCBI's.
    /// The subject of the definition lines `set` of a database. `Err` names what LOSAT does not
    /// write as NCBI: a stored Seq-id with a byte outside ' ' ..= '~' (NCBI writes an error
    /// message each time it reads the definition line; such a title is the reports' refusal,
    /// `replaced_title_bytes`). A definition line without a Seq-id (makeblastdb writes none such)
    /// beside one with a Seq-id is a subject as any other: NCBI's scope loads it, and each format
    /// that walks its definition lines meets that line where its own code says
    /// (`has_lineless_line`). A set without any Seq-id gets no Seq-align
    /// (`algorithm::common::database::SeqAligned` drops its hits before the reports), so it never
    /// comes here: LOSAT's internal error.
    pub fn new(set: DeflineSet) -> Result<Self, String> {
        if set
            .deflines
            .iter()
            .any(|defline| defline.replaced_id_bytes > 0)
        {
            return Err("has a stored Seq-id with a byte outside ' ' to '~' (NCBI BLAST+ shows it as '#' and writes an error message each time it reads the definition line)".into());
        }
        let Some(alignment_id) = set.alignment_id().cloned() else {
            return Err(
                "has no Seq-id in any definition line, which gets no Seq-align (internal error)"
                    .into(),
            );
        };
        let Some(own_line) = set.bioseq_defline(&alignment_id) else {
            return Err("has no definition line that holds its Seq-align id".into());
        };
        // A giim Seq-align id is never found again: NCBI's handles never take it for one with an
        // equal id (`algorithm::psiblast::rounds::HandleIndex`), so its scope finds no Bioseq by it
        // (NCBI BLAST+ 2.17.0 writes `Unknown`; review 2 of line DBID pssmrows, H1).
        let unfetchable = matches!(alignment_id, ParsedSeqId::Giim(_))
            || set
                .display_defline(&alignment_id)
                .is_some_and(|shown| shown.index != own_line);
        let lookup = handle_lookup(&set, &alignment_id, own_line);
        let target = match lookup {
            Some(HandleLookup::Line(line)) => line,
            _ => own_line,
        };
        Ok(Self {
            set,
            alignment_id,
            target,
            unfetchable,
            resolution: ScopeResolution::default(),
            scope_finds_none: lookup == Some(HandleLookup::NoBioseq),
            load_throws: lookup == Some(HandleLookup::Throws),
            sequence: None,
            show_gis: false,
            source: SourceTaxname::Unnamed,
            record: false,
            local_record: false,
            use_this_seqid: Vec::new(),
            saccver: std::sync::OnceLock::new(),
            query_scope_keys: std::sync::OnceLock::new(),
            skipped_id: std::sync::OnceLock::new(),
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1305-1470
    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/cached_sequence.cpp:46-73
    /// The Bioseq of the sequence `holder` (its record as the reports read it: its source, GIs and
    /// sequence) that NCBI's places give for another subject's Seq-align id `alignment_id`, loaded
    /// by its line `line` (the B-line, `ShownInstead::written`): its ids and title are that
    /// line's (`tabular_ids`, `xml_hit`'s accession), every other writer reads all its lines, and
    /// the reports name the Seq-align id of the subject (`alignment_id`). Its lookup is the
    /// places' answer, so nothing of the holder's own lookup is refused for it. `Err` for a line
    /// that the holder does not have (LOSAT's internal error).
    pub fn shown_by(
        holder: &ParsedSubject,
        line: usize,
        alignment_id: &ParsedSeqId,
    ) -> Result<Self, String> {
        if line >= holder.set.deflines.len() {
            return Err(format!(
                "is shown by definition line {line} of a sequence of {} (internal error)",
                holder.set.deflines.len()
            ));
        }
        Ok(Self {
            set: holder.set.clone(),
            alignment_id: alignment_id.clone(),
            target: line,
            unfetchable: false,
            resolution: ScopeResolution::default(),
            scope_finds_none: false,
            load_throws: false,
            sequence: holder.sequence,
            show_gis: holder.show_gis,
            source: holder.source.clone(),
            record: false,
            local_record: false,
            use_this_seqid: holder.use_this_seqid.clone(),
            saccver: std::sync::OnceLock::new(),
            query_scope_keys: std::sync::OnceLock::new(),
            skipped_id: std::sync::OnceLock::new(),
        })
    }

    /// The sequence of the database whose Bioseq this record is (`in_database`; the holder of a
    /// record of [`ParsedSubject::shown_by`]), `None` for a `-subject` record: its length is the
    /// subject's `slen` (`report::shown_bioseq`).
    pub fn database_sequence(&self) -> Option<DbSequence> {
        self.sequence
    }

    /// Whether the header of this sequence was read with a UTF8String, whose load NCBI posts a
    /// warning at (`common::database::HeaderReader::note_utf8`).
    pub fn header_utf8(&self) -> bool {
        self.set.utf8
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:458-466
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:257-264
    /// The `-subject` record with the ids `ids` and the FASTA title `title` under
    /// `-parse_deflines`, as one definition line: its Bioseq has no Blast-def-line user object, so
    /// its titles are generated and `stitle` and `salltitles` are `N/A`, and NCBI's scope holds
    /// it, so it is always found.
    pub fn of_record(
        ids: &crate::blastinput::fasta_reader::believed::BelievedIds,
        title: &[u8],
    ) -> Self {
        Self {
            set: DeflineSet {
                deflines: vec![BlastDefline {
                    title: (!title.is_empty()).then(|| title.to_vec()),
                    ids: ids.ids.clone(),
                    ..BlastDefline::default()
                }],
                utf8: false,
            },
            alignment_id: ids.best.clone(),
            target: 0,
            unfetchable: false,
            resolution: ScopeResolution::default(),
            scope_finds_none: false,
            load_throws: false,
            sequence: None,
            show_gis: false,
            source: SourceTaxname::Unnamed,
            record: true,
            local_record: ids.subject_local(),
            use_this_seqid: Vec::new(),
            saccver: std::sync::OnceLock::new(),
            query_scope_keys: std::sync::OnceLock::new(),
            skipped_id: std::sync::OnceLock::new(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:482-498
    /// Whether the outfmt 0, XML and archive reports show this `-subject` record as one read
    /// without `-parse_deflines` (a local id `Subject_N`, `BelievedIds::subject_local`); its
    /// `sall*` are still its own ids (`s_ReplaceLocalId` keeps a local id with `parse_local`).
    pub fn shown_as_record(&self) -> bool {
        self.local_record
    }

    /// LOSAT's refusal of this `-subject` record (`of_record`) for `reason`, which names the
    /// record by its description id (`description_start`) and `program`.
    pub fn record_refusal(&self, reason: &str, program: &str) -> anyhow::Error {
        anyhow::anyhow!(
            "the subject record {} {reason}, which is not supported by LOSAT's {program}",
            String::from_utf8_lossy(&self.description_start()).trim_end()
        )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:393-416
    // NCBI's data loader serves the Bioseq of the first OID that `SeqidToOid` gives for the
    // Seq-align id (`seqdb::lookup::SeqidLookup::seqid_to_oid`), and only where one of that
    // sequence's ids compares equal to it: a version 5 database finds `gb|X|` under the key `X`,
    // which makeblastdb also writes for `gb|X.1|`, and `gb|X.1|` compares equal to `gb|X|` and
    // to `emb|X.1|`. NCBI 2.17.0 then writes the other sequence's names and counts the
    // alignment's identities against its letters (outfmt 6 `AAA00001.1 100.000 ... 0` for the
    // subject `gb|AAA00001|` of 99.561 % identity), or, for `lcl|X` after `gb|X.1|`, finds no
    // Bioseq and stops with a CCoreException in outfmt 6 (probe of line DBID N1,
    // `blastdb-ids-20261021/R/accoid/p1`; version 4 alike).
    /// The subject as NCBI's formatters find it again: `found` is what NCBI's lookup of the
    /// Seq-align id in the database finds where that is not this sequence (another one or
    /// none), which the formats that read the subject's Bioseq refuse (`fetched`).
    #[cfg(test)]
    pub fn found_instead(mut self, found: Option<FoundInstead>) -> Self {
        self.resolution.found = found;
        self
    }

    /// The subject with `resolution`, what NCBI's lookup of the Seq-align id gives
    /// ([`ScopeResolution`]; `common::database::HeaderReader::resolution`).
    pub fn resolved(mut self, resolution: ScopeResolution) -> Self {
        self.resolution = resolution;
        self
    }

    /// The subject whose lookup gives the sequence `other` first (another sequence, or a copy of
    /// its record that the reports show as the subject: `ScopeResolution::gives`).
    #[cfg(test)]
    pub fn lookup_gives(mut self, other: Option<DbSequence>) -> Self {
        self.resolution.gives = other;
        self
    }

    /// The subject with `holders`, the other sequences that NCBI's places may give for the
    /// Seq-align id after an API call before this one, listed from the database
    /// ([`ScopeHolders`]).
    #[cfg(test)]
    pub fn unrecorded_holders(mut self, holders: ScopeHolders) -> Self {
        self.resolution.unrecorded = holders;
        self
    }

    /// The subject as the sequence `sequence` of its database, by which the run's loads name it
    /// ([`ScopeWalk`]).
    pub fn in_database(mut self, sequence: DbSequence) -> Self {
        self.sequence = Some(sequence);
        self
    }

    /// Whether `sequence` is this subject as NCBI's load gives it: the subject itself or the
    /// sequence that its lookup gives first (a copy of its record, or another one whose refusal
    /// is the lookup's, `FoundInstead::Sequence`).
    fn mine(&self, sequence: DbSequence) -> bool {
        self.sequence == Some(sequence) || self.resolution.gives == Some(sequence)
    }

    /// Non-allocating reference to the blob that a load of this subject by its own Seq-align id
    /// loads where no place before the lookup answers (step 4 of [`ScopeWalk::load`]).
    fn loaded_blob_ref(&self) -> Option<LoadedBlobRef> {
        use FoundInstead as F;
        match &self.resolution.found {
            Some(F::Failed(_) | F::ThrowsOnId(_) | F::Nothing) => None,
            Some(F::Sequence(_)) => {
                let blob = self.resolution.gives_keys.clone()?;
                let sequence = blob.sequence?;
                let line = blob.target;
                Some(LoadedBlobRef::External {
                    blob,
                    sequence,
                    line,
                })
            }
            None | Some(F::OtherTaxids(_) | F::CopyApart(..)) => Some(LoadedBlobRef::Own {
                sequence: self.resolution.gives.or(self.sequence)?,
                line: self.target,
            }),
        }
    }

    /// Exact handles of every Seq-id of every definition line of this subject, without allocation.
    fn line_handles(&self) -> impl Iterator<Item = HandleKey> + '_ {
        self.set
            .deflines
            .iter()
            .flat_map(|defline| defline.ids.iter())
            .filter_map(handle_key)
    }

    /// First definition line of this subject that holds `handle`, without allocation.
    fn line_of_handle(&self, handle: HandleKey) -> Option<usize> {
        self.set
            .deflines
            .iter()
            .position(|defline| defline.ids.iter().any(|id| handle_key(id) == Some(handle)))
    }

    /// The exact handle key of this subject's Seq-align id, without allocating `IdKeys`.
    fn id_handle(&self) -> Option<HandleKey> {
        handle_key(&self.alignment_id)
    }

    /// The exact handle key and bulk sort order of this subject's Seq-align id, without
    /// allocating `IdKeys`.
    fn id_handle_order(&self) -> (Option<HandleKey>, u64) {
        (
            handle_key(&self.alignment_id),
            crate::seqid::mask_handle::mask_handle_order(
                crate::seqid::mask_handle::mask_handle_packed(&self.alignment_id),
            ),
        )
    }

    /// Compact flat [`BlobKeys`] for this subject as `sequence`.
    fn blob_keys_as(&self, sequence: DbSequence) -> std::sync::Arc<BlobKeys> {
        std::sync::Arc::new(BlobKeys::of_deflines(
            Some(sequence),
            self.target,
            &self.set.deflines,
        ))
    }

    /// The sequence whose blob a load of `other` by its own Seq-align id loads
    /// ([`ParsedSubject::loaded_blob_ref`]) where that blob registers the handle of this subject's
    /// Seq-align id (`RegisterIds`: a line of it holds an id of that handle; any blob where the id
    /// has no handle that a function can know, `id_handle_unknown`): the loads before the stop of
    /// a prefetch (`common::database::stop_lookup_unsure`).
    pub fn registered_by(&self, other: &ParsedSubject) -> Option<DbSequence> {
        let blob = other.loaded_blob_ref()?;
        let sequence = blob.sequence();
        match self.id_handle() {
            Some(handle) => match blob {
                LoadedBlobRef::Own { .. } => other.line_of_handle(handle).map(|_| sequence),
                LoadedBlobRef::External { blob, .. } => blob.line_of(handle).map(|_| sequence),
            },
            None => Some(sequence),
        }
    }

    // No NCBI counterpart: the keys of LOSAT's own check before a search whether any load can be
    // given another sequence (`common::database::reports_may_refuse`).
    /// The keys by which NCBI's places can find the blob that a load of this subject loads
    /// ([`ParsedSubject::loaded_blob_ref`]: each line's exact handles, and the scope's keys of the
    /// line of its Bioseq), each with that blob's ordinal id (`held`), and those by which a load of
    /// its Seq-align id looks for one: its handle in `m_Ids` (`sought`), and the keys of the
    /// scope's index (`scope_sought`, which only a load one by one reads, `ScopeWalk::bulk`): one
    /// key space, so that a key held by a blob that is neither the subject's nor its lookup's
    /// means that a load may give another sequence (a key of 64 bits: two keys that share one only
    /// make the check refuse more).
    pub fn place_keys(
        &self,
        held: &mut Vec<(u64, u32)>,
        (sought, scope_sought): (&mut Vec<u64>, &mut Vec<u64>),
    ) {
        fn tagged(tag: u64, key: u64) -> u64 {
            key.rotate_left(17 * tag as u32) ^ tag.wrapping_mul(0x9e37_79b9_7f4a_7c15)
        }
        if let Some(blob) = self.loaded_blob_ref() {
            let ordinal = blob.sequence().ordinal();
            let line = blob.line();
            match blob {
                LoadedBlobRef::Own { .. } => {
                    for (at, defline) in self.set.deflines.iter().enumerate() {
                        for id in &defline.ids {
                            if let Some(h) = handle_key(id) {
                                held.push((tagged(0, h as u64), ordinal));
                            }
                            if at == line {
                                let mut scope = Vec::new();
                                crate::seqid::defline_set::scope_keys(id, &mut scope);
                                held.extend(scope.into_iter().map(|k| (tagged(1, k), ordinal)));
                                let mut kinds = Vec::new();
                                crate::seqid::id_scope::kind_keys(id, &mut kinds);
                                held.extend(kinds.into_iter().map(|k| (tagged(2, k), ordinal)));
                            }
                        }
                    }
                }
                LoadedBlobRef::External { blob, .. } => {
                    held.extend(
                        blob.all_handles()
                            .iter()
                            .map(|&h| (tagged(0, h as u64), ordinal)),
                    );
                    held.extend(
                        blob.line_scope(line)
                            .iter()
                            .map(|&k| (tagged(1, k), ordinal)),
                    );
                    held.extend(
                        blob.line_kinds(line)
                            .iter()
                            .map(|&k| (tagged(2, k), ordinal)),
                    );
                }
            }
        }
        let keys = self.id_keys();
        sought.extend(keys.handle.map(|h| tagged(0, h as u64)));
        scope_sought.extend(keys.scope.iter().map(|&k| tagged(1, k)));
        scope_sought.extend(keys.weak.iter().map(|&k| tagged(2, k)));
    }

    /// The ordinal id of the sequence that the lookup of this subject's Seq-align id gives
    /// first where it is another one (`ScopeResolution::gives`).
    pub fn lookup_gives_ordinal(&self) -> Option<u32> {
        self.resolution.gives.map(DbSequence::ordinal)
    }

    /// Whether this subject's Seq-align id has no `CSeq_id_Handle` key that LOSAT knows
    /// (`seqid::mask_handle::mask_handle_key`: a giim id, or a Textseq-id with neither accession
    /// nor name), so that the stop of a prefetch takes any load for one that may register it.
    pub fn id_handle_unknown(&self) -> bool {
        self.id_handle().is_none()
    }

    /// The other sequences that NCBI's places may give for this subject's Seq-align id after the
    /// calls of the run before, of which LOSAT has no record ([`ScopeHolders`]; tests).
    #[cfg(test)]
    pub fn unrecorded(&self) -> &ScopeHolders {
        &self.resolution.unrecorded
    }

    /// The other sequences that, in an API call after the first, the data loader may give for
    /// this subject's Seq-align id from the calls before ([`ScopeHolders::loader`]), or why LOSAT
    /// cannot list them.
    pub fn unrecorded_loader(&self) -> Result<&[DbSequence], &str> {
        match &self.resolution.unrecorded.loader {
            Holders::Listed(holders) => Ok(holders),
            Holders::Unknown(why) => Err(why),
        }
    }

    /// The keys of the sequence `sequence`, whose definition lines are `set`, as a load by this
    /// subject's Seq-align id loads it ([`ScopeResolution::gives_keys`]): its Bioseq is made by
    /// the first line that holds an id of the Seq-align id's handle (the first line where none
    /// does).
    pub fn shown_keys(
        &self,
        set: &crate::seqid::DeflineSet,
        sequence: DbSequence,
    ) -> std::sync::Arc<BlobKeys> {
        let mut keys = BlobKeys::of_deflines(Some(sequence), 0, &set.deflines);
        keys.target = self
            .id_handle()
            .and_then(|handle| keys.line_of(handle))
            .unwrap_or(0);
        std::sync::Arc::new(keys)
    }

    /// The keys of the Seq-align id by which NCBI's places find a sequence for it ([`IdKeys`]).
    pub fn id_keys(&self) -> IdKeys {
        IdKeys::of(&self.alignment_id)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:432-433, 497-498
    /// The keys of the display id ([`ParsedSubject::display_id`]) where its handle is not the
    /// Seq-align id's: the description table (outfmt 0 to 4) looks the title up again by it in
    /// NCBI's places ([`ScopeWalk::display`]), which can give another sequence's Bioseq (probe
    /// `R/hr/bin/hr1open.py` D.ab: two records with other GIs and one accession). `None` where
    /// that lookup is the Seq-align id's own.
    pub fn display_keys(&self) -> Option<IdKeys> {
        let display = self.display_id()?;
        let display_handle = handle_key(display)?;
        if Some(display_handle) == self.id_handle() {
            return None;
        }
        Some(IdKeys::of(display))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:411-451 (cited at
    // `algorithm::psiblast::pssm::input_msa`)
    /// The reason of LOSAT's refusal of psiblast's PSSM read of this subject where NCBI's places
    /// give its Seq-align id another sequence (`shown`, [`ScopeWalk`]), whose letters the read
    /// takes for it (where they reach the end of an included Seq-align's range,
    /// `ListedTitles::check_pssm_holds` of BLASTP). Not the words of a report that shows the
    /// holder (`shown_reason`): the PSSM reads by the id, also beyond `-max_target_seqs` and in
    /// the formats that show no subject (review 5 of line DBID pfoff-4, L1).
    pub fn pssm_held_reason(&self, shown: ShownInstead) -> String {
        let id = String::from_utf8_lossy(&self.alignment_fasta()).into_owned();
        let failing = if self.lookup_throws().is_some() {
            "whose lookup in the database fails, and "
        } else {
            ""
        };
        match (shown.holder, shown.by) {
            (Some(other), ShownBy::Holder) => format!(
                "has the Seq-align id {id}, {failing}which database sequence {other} also holds: psiblast's PSSM reads the subject by that id (also beyond -max_target_seqs), and NCBI BLAST+ 2.17.0 takes that sequence's letters for it"
            ),
            (Some(other), ShownBy::Weak) => format!(
                "has the Seq-align id {id}, {failing}which can match an id of database sequence {other} weakly (an accession with a version for one without, a name with or without an accession, a pir or prf accession as a name): psiblast's PSSM reads the subject by that id (also beyond -max_target_seqs), and NCBI BLAST+ 2.17.0 takes that sequence's letters for it"
            ),
            _ => {
                let why = self.shown_unknown(shown);
                format!(
                    "has the Seq-align id {id}, {failing}for which LOSAT cannot list {why}: psiblast's PSSM reads the subject by that id (also beyond -max_target_seqs), and where such a sequence was loaded before, NCBI BLAST+ 2.17.0 takes its letters for it"
                )
            }
        }
    }

    /// What LOSAT cannot tell of a sequence that NCBI's places may give for the Seq-align id
    /// (`shown`, where [`ShownBy::Order`] or [`ShownBy::Calls`] gives no holder or an order LOSAT
    /// does not know): the words of its refusal.
    fn shown_unknown(&self, shown: ShownInstead) -> String {
        match shown.by {
            ShownBy::Calls => self
                .resolution
                .unrecorded
                .loader
                .unknown()
                .or(self.resolution.unrecorded.scope.unknown())
                .map_or_else(
                    || "the sequences that the calls of the run before loaded".to_owned(),
                    str::to_owned,
                ),
            ShownBy::Order | ShownBy::Holder | ShownBy::Weak | ShownBy::Display => match shown.holder {
                Some(other) => format!(
                    "whether database sequence {other}, which also holds it, is loaded before it (an order of NCBI BLAST+'s loads that their Seq-ids do not decide)"
                ),
                None => "the sequences that hold it, as NCBI BLAST+ loads them in an order that their Seq-ids do not decide".to_owned(),
            },
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:393-406
    /// The subject with `other`, the sequence whose header the lookup of the Seq-align id reads
    /// first where it holds a UTF8String where a VisibleString is expected (its version check of
    /// an accession, or the sequence it gives first where that is another one;
    /// `common::database::HeaderReader`): NCBI's data loader reads that header (`x_GetOid`'s
    /// `GetSeqIDs` of the first OID) wherever its scope resolves the id, and posts the warning
    /// there where the run has not (probe of line RG item T5d fix 3): the load of the id in LOSAT's
    /// model of the loads (`common::database::load_read_texts`: a batch's prefetch, a report that
    /// reads its subjects one by one, psiblast's PSSM).
    #[cfg(test)]
    pub fn scope_utf8(mut self, other: Option<DbSequence>) -> Self {
        self.resolution.utf8 = other;
        self
    }

    /// `scope_utf8`: the other sequence whose header NCBI reads with its warning where it
    /// resolves this subject's Seq-align id.
    pub fn scope_utf8_read(&self) -> Option<DbSequence> {
        self.resolution.utf8
    }

    /// Whether NCBI's formatters find no Bioseq by the Seq-align id whatever its lookup gives
    /// (see `new`), so that no lookup is needed. Read by tests only: the engine asks
    /// `common::database::scope_bioseq`, which looks it up first.
    #[cfg(test)]
    pub fn unfetchable(&self) -> bool {
        self.unfetchable
    }

    /// Whether NCBI's scope gets no Bioseq for the Seq-align id although the lookup gives this
    /// sequence (`handle_lookup`): the Bioseq handle that the formatters get is empty, and each
    /// format writes the `catch` of its read, or skips the subject, as NCBI does (outfmt 0
    /// `Unknown` and `no longer exists in database`, XML `Unknown`, the tabular rows of a failed
    /// lookup, `LookupFailure::empty_handle`); the prefetch and psiblast's PSSM read nothing.
    pub fn scope_finds_none(&self) -> bool {
        self.scope_finds_none
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1576-1589
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:886-892
    /// Whether the outfmt 0 to 4 description line of the subject prints a GI that nothing set:
    /// with an empty Bioseq handle (`scope_finds_none`) NCBI's description info keeps the GI it
    /// was allocated with, which `-show_gis` prints when it is positive (a different number in
    /// each run; probe `probe_gis` of line RG-KB item 12).
    pub fn description_gi_unset(&self) -> bool {
        self.scope_finds_none && self.show_gis
    }

    /// Whether the data loader's build of the Bioseq throws (`load_throws`): another id of the
    /// sequence holds the Seq-id of the Seq-align id's handle, while no definition line holds it
    /// as NCBI's load reads the ids.
    pub fn load_throws(&self) -> bool {
        self.load_throws
    }

    /// LOSAT's refusal of a read of the subject's Bioseq where NCBI's load of it throws
    /// (`load_throws`): the formats that read the Bioseq (`fetched`), the prefetch
    /// (`prefetched`) and psiblast's PSSM (`pssm_read`); the Seq-annots and the archive (outfmt 8,
    /// 9, 11, 12) read nothing of it.
    pub(crate) fn handle_load_refusal(&self) -> Option<String> {
        self.load_throws.then(|| {
            format!(
                "has the Seq-align id {}, whose Seq-id NCBI BLAST+'s data loader looks up with another chain that another of its ids matches, while no definition line holds it as NCBI BLAST+ 2.17.0's load reads the ids (its load throws: {} Error: oid headers do not contain target gi/seq_id.)",
                String::from_utf8_lossy(&self.alignment_fasta()),
                crate::wording::NCBI_DATABASE_ERROR
            )
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:389-411
    /// The Seq-id that NCBI's data loader looks the Seq-align id up as and compares with the ids
    /// of the sequence its lookup gives (`CSeq_id_Handle::GetSeqId` of the id's handle): the
    /// Seq-align id but for a PDB `chain` that the handle normalizes
    /// (`seqid::mask_handle::pdb_handle_seq_id`).
    pub fn loader_id(&self) -> std::borrow::Cow<'_, ParsedSeqId> {
        match handle_changes(&self.alignment_id) {
            Some(handle) => std::borrow::Cow::Owned(handle),
            None => std::borrow::Cow::Borrowed(&self.alignment_id),
        }
    }

    /// What NCBI's lookup of the Seq-align id finds where that is not this sequence
    /// (`found_instead`).
    pub fn found(&self) -> Option<&FoundInstead> {
        self.resolution.found.as_ref()
    }

    /// Whether unrecorded-call holders for this subject could not be listed
    /// (`Holders::Unknown`).
    pub fn holders_unknown(&self) -> bool {
        self.resolution.unrecorded.loader.unknown().is_some()
            || self.resolution.unrecorded.scope.unknown().is_some()
    }

    /// The subject as the outfmt 0 report of a run with or without `-show_gis` shows it.
    pub fn showing_gis(mut self, show_gis: bool) -> Self {
        self.show_gis = show_gis;
        self
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1146-1168
    // A title made from the source (`x_SetTitleFromBioSrc`, of the first source of the Bioseq,
    // whose volume puts the target's taxids first) names the organism only where the taxonomy
    // database gives a name for the taxid; else NCBI writes the title of a Bioseq without one
    // (review 1 of line DBID N1b L1, `blastdb-ids-20261021/R/titlesrc/review1/r2.tsv`).
    // A RefSeq chromosome, mRNA or ncRNA record gets its title from the source too: NCBI's own
    // functions for them make nothing of a database Bioseq (`report::defline::source_title`).
    // The end of a protein title ending in `]` is rewritten from the same name
    // (`x_AdjustProteinTitleSuffix`, `adjusted_title_shown`; probe of line RG T5c).
    /// The subject with the name of its Bioseq's first source (`names` of its taxid,
    /// `first_source_taxid`, asked only where `title_from_taxonomy` or, of a `protein` database,
    /// `adjusted_title_shown`), from which its generated title is made (`source`,
    /// `bioseq_title`). Other subjects ask nothing.
    pub fn with_source(mut self, protein: bool, names: impl FnOnce(i32) -> SourceTaxname) -> Self {
        if self.title_from_taxonomy() || self.adjusted_title_shown(protein) {
            if let Some(taxid) = self.first_source_taxid() {
                self.source = names(taxid);
            }
        }
        self
    }

    /// The taxid of the subject Bioseq's first source where NCBI fetches it by its Seq-align id
    /// (`DeflineSet::first_source_taxid`, the first of `volume_lookups`).
    pub fn first_source_taxid(&self) -> Option<i32> {
        self.set.first_source_taxid(Some(&*self.loader_id()))
    }

    /// The name of the organism of the subject Bioseq's first source (`with_source`).
    pub fn source(&self) -> &SourceTaxname {
        &self.source
    }

    /// The subject as a database with a positive list (`positive_list`,
    /// `SeqDb::has_positive_list`) reports it: its Seq-aligns name the ids of the definition
    /// lines the list keeps. `use_gi` is whether the Seq-align id is a GI.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqinfosrc_aux.cpp:113-150
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1604-1606
    pub fn under_positive_list(mut self, positive_list: bool) -> Self {
        if !positive_list {
            return self;
        }
        let use_gi = self.alignment_id.gi().is_some();
        for id in self.set.ids() {
            let prefix: &[u8] = if use_gi {
                if id.gi().is_none() {
                    continue;
                }
                b"gi:"
            } else {
                if !(id.textseq().is_some_and(|text| text.accession.is_some())
                    || matches!(id, ParsedSeqId::Pdb(_)))
                {
                    continue;
                }
                b"seqid:"
            };
            let mut string = prefix.to_vec();
            crate::seqid::label::content(id, true, &mut string);
            self.use_this_seqid.push(string);
        }
        self
    }

    /// The strings of the Seq-aligns' `use_this_seqid` (`under_positive_list`).
    pub fn use_this_seqid(&self) -> &[Vec<u8>] {
        &self.use_this_seqid
    }

    /// `Err` for a subject whose sequence NCBI's formatters cannot find by its Seq-align id
    /// (see `new` and `found_instead`): the formats that read the subject's Bioseq (outfmt 0, 5,
    /// 6, 7, 10, 13 to 16 and SAM, 17) refuse it; the Seq-annots and the archive (8, 9, 11, 12)
    /// write the Seq-align id alone. `taxonomy`: whether a report writes the subject's taxids
    /// or their names (outfmt 13 to 16, the tabular `staxid`, `staxids` and the names), which
    /// alone refuse a sequence found as the same record with other taxids
    /// (`FoundInstead::OtherTaxids`; also where a generated title can be made from them, which
    /// the reports that show it refuse, `generated_title_fetched`). `load`: what the reports'
    /// loads of the subject give (`HitLoads::load`, the one answer of [`ScopeWalk`]): a subject
    /// whose `Hit` throws (`HitLoad::Throws`: its lookup fails on a file with NCBI's exception in
    /// a batch whose prefetch is off, and every format of the run writes it as NCBI does,
    /// `FormatThrow::of_run`: BLAST XML (outfmt 5), the XML2 and JSON documents of outfmt 13 to
    /// 16, the tables of outfmt 6, 7, 10 and 20, the organism report of outfmt 18, the pairwise
    /// reports of outfmt 0 to 4) is not refused for what the lookup finds (it finds nothing: it
    /// throws before any definition line is read), and one that no report loads
    /// (`HitLoad::Unread`: shown only after such a `Hit`, or as a copy of the record in it) is not
    /// refused at all; one for which NCBI's places give another sequence (`HitLoad::Shown`: its
    /// scope's index, its data loader's `m_Ids`, the prefetch's order, whatever the lookup would
    /// give) is refused with the words of that sequence (`shown_reason`: LOSAT does not write
    /// another sequence's record for a subject yet); every other subject (`HitLoad::Loaded`,
    /// `HitLoad::Refused`) is refused as above, a lookup that fails on a file among them (outfmt
    /// 17 here, before any byte is written; the reason of a run whose formats refuse it,
    /// `HitLoad::Refused`).
    pub fn fetched(&self, taxonomy: bool, load: HitLoad) -> Result<(), String> {
        self.fetched_with_scope(taxonomy, load, false)
    }

    /// `fetched` with `seen_by_scope`: whether every reporting query's scope resolves this subject
    /// to a noted query Bioseq (`report::query_bioseq::QueryBioseqs::hides_ordinal`), in which
    /// case NCBI's scope returns the query Bioseq before calling the BLAST database data loader
    /// (`handle_load_refusal` does not apply).
    pub fn fetched_with_scope(
        &self,
        taxonomy: bool,
        load: HitLoad,
        seen_by_scope: bool,
    ) -> Result<(), String> {
        match load {
            HitLoad::Unread | HitLoad::Throws => Ok(()),
            HitLoad::Shown(shown) => Err(self.shown_reason(shown)),
            HitLoad::Loaded => self.found_fetched(taxonomy, None, seen_by_scope),
            HitLoad::Refused(refused) => self.found_fetched(taxonomy, Some(refused), seen_by_scope),
        }
    }

    /// The reason of LOSAT's refusal of a subject for which NCBI's places give another sequence
    /// at a load (`shown`, [`ScopeWalk`]), whatever its lookup gives: the words of a sequence that
    /// also holds the Seq-align id (its exact handle, or an id that the scope takes for it), of a
    /// weak match, or of what LOSAT cannot tell (an order its record does not keep, the calls of
    /// the run before).
    pub fn shown_reason(&self, shown: ShownInstead) -> String {
        let id = String::from_utf8_lossy(&self.alignment_fasta()).into_owned();
        match (shown.holder, shown.by) {
            (other, ShownBy::Display) => {
                let mut display = Vec::new();
                if let Some(id) = self.display_id() {
                    label::fasta(id, &mut display);
                }
                format!(
                    "has the display id {} (its description line's), by which the description table of NCBI BLAST+ 2.17.0's outfmt 0 to 4 looks the title up again and finds {}: it can write that sequence's title for it",
                    String::from_utf8_lossy(&display),
                    other.map_or_else(
                        || "another sequence".to_string(),
                        |other| format!("database sequence {other}")
                    )
                )
            }
            (Some(other), ShownBy::Holder | ShownBy::Calls) => format!(
                "has the Seq-align id {id}, which database sequence {other} also holds (where both are reported, NCBI BLAST+ 2.17.0 can write that sequence's names and letters for it)"
            ),
            (Some(other), ShownBy::Weak) => format!(
                "has the Seq-align id {id}, which can match an id of database sequence {other} weakly (an accession with a version for one without, a name with or without an accession, a pir or prf accession as a name): in a batch whose prefetch is off (PRE_FETCH_SEQS_LIMIT, or a format that NCBI BLAST+ 2.17.0 does not prefetch for), where NCBI BLAST+ 2.17.0 reads the reported sequences one by one and that sequence first, it can write that sequence's names and letters for it"
            ),
            _ => format!(
                "has the Seq-align id {id}, for which LOSAT cannot list {}: where such a sequence is also reported, NCBI BLAST+ 2.17.0 can write its names and letters for it",
                self.shown_unknown(shown)
            ),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:357-360, 412-417
    // Where the batch's prefetch is off, the `Hit` loads its subject by the Seq-align id: the
    // lookup throws on the file there, and again in the `catch` (`GetLength` loads by the same id),
    // outside any `catch` of the hit; the report function of outfmt 15 and 16 catches it, 13 and 14
    // end the program with it (probe of line DBID mixdb, `PU.off`; DECISIONS 62 (2)).
    /// The exception that NCBI's lookup of the Seq-align id throws where a formatter loads the
    /// subject one by one (a batch whose prefetch is off, or fails): a lookup that fails on a file
    /// with NCBI's `CSeqDBException` whose place LOSAT knows (`FoundInstead::Failed`,
    /// `LookupFailure::exception`; where the batch prefetches, its prefetch throws first and
    /// stops the run, `prefetch_stops`), or one that throws on the id itself
    /// (`FoundInstead::ThrowsOnId`, whose prefetch catches it, `prefetch_fails`). The one owner of
    /// "the load by this id throws": the walk of the reports' loads (`HitLoads`) and every
    /// writer's answer to it (`FormatThrow`) take it, whatever the cause.
    pub fn lookup_throws(&self) -> Option<crate::ncbi_exception::NcbiException> {
        match &self.resolution.found {
            Some(FoundInstead::Failed(failure)) => failure.exception(),
            Some(FoundInstead::ThrowsOnId(exception)) => Some(exception.clone()),
            _ => None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:440-450 (cited at
    // `algorithm::psiblast::pssm::input_msa`)
    /// LOSAT's refusal of what psiblast's PSSM reads by the Seq-align id of a subject that it
    /// includes (the `CSeqVector` of each Seq-align's range in NCBI's scope), where LOSAT does not
    /// read it as NCBI: what the scope gets is `common::database::scope_bioseq`'s, and this
    /// refuses a lookup that fails on a file with LOSAT's error (`FoundInstead::Failed`; with
    /// NCBI's exception, which the read catches, `prefetch_stops`: the caller reads no letters)
    /// and the same record whose read writes other taxonomy texts (`CopyApart::TaxonomyTexts`).
    /// Another sequence that NCBI's places give for the id is the walk's (`ListedTitles::
    /// check_pssm_holds` of BLASTP, [`ScopeWalk`]). The reasons are what the lookup finds, also
    /// for a subject whose definition line the formatters do not find (`unfetchable`), which the
    /// PSSM reads by the scope as any other (review 2 of line DBID pssmread, M2).
    pub fn pssm_read(&self) -> Result<(), String> {
        if let Some(reason) = self.handle_load_refusal() {
            return Err(reason);
        }
        match &self.resolution.found {
            Some(
                found @ (FoundInstead::Failed(_)
                | FoundInstead::CopyApart(_, CopyApart::TaxonomyTexts)),
            ) => Err(self.found_reason(found, false, None).unwrap_or_default()),
            _ => Ok(()),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1437-1439
    /// Whether the Bioseq that `GetBioseqNoData` builds of this sequence for its Seq-align id
    /// holds that id: the ids of the definition line its filter keeps (`target`,
    /// `DeflineSet::bioseq_defline`) hold one that compares equal to the Seq-id of the id's handle
    /// (`loader_id`; the scope's test on the Bioseq's ids, `common::database::scope_bioseq`). Not
    /// where that line is found by the content of a local id alone (`lcl|777` before `gi|777`),
    /// nor for a giim id.
    pub fn bioseq_holds_alignment_id(&self) -> bool {
        let loader_id = self.loader_id();
        !matches!(self.alignment_id, ParsedSeqId::Giim(_))
            && self
                .target()
                .ids
                .iter()
                .any(|id| compare(id, &loader_id) == Sic::Yes)
    }

    /// `fetched` of a subject that the reports load: what the lookup of the Seq-align id finds
    /// (`found_instead`), whatever the batch's prefetch.
    fn found_fetched(
        &self,
        taxonomy: bool,
        refused: Option<ThrowRefusal>,
        seen_by_scope: bool,
    ) -> Result<(), String> {
        if !seen_by_scope {
            if let Some(reason) = self.handle_load_refusal() {
                return Err(reason);
            }
        }
        if self.unfetchable {
            return Err(format!(
                "has the Seq-align id {}, by which NCBI BLAST+'s formatters do not find the sequence (NCBI BLAST+ 2.17.0 writes Unknown or another sequence's names, or stops)",
                String::from_utf8_lossy(&self.alignment_fasta())
            ));
        }
        match &self.resolution.found {
            None => Ok(()),
            Some(found) => self
                .found_reason(found, taxonomy, refused)
                .map_or(Ok(()), Err),
        }
    }

    /// The reason of a refusal of what the lookup of the Seq-align id finds, `found`, for a format
    /// that reads the Bioseq (`fetched`; `taxonomy` as there), or for psiblast's PSSM
    /// (`pssm_read`); `None` for a copy of the record that such a report shows as the subject.
    fn found_reason(
        &self,
        found: &FoundInstead,
        taxonomy: bool,
        refused: Option<ThrowRefusal>,
    ) -> Option<String> {
        let id = String::from_utf8_lossy(&self.alignment_fasta()).into_owned();
        Some(match found {
            FoundInstead::OtherTaxids(_) | FoundInstead::CopyApart(_, CopyApart::Title)
                if !taxonomy =>
            {
                return None
            }
            FoundInstead::OtherTaxids(other) | FoundInstead::CopyApart(other, CopyApart::Title) => format!(
                "has the Seq-align id {id}, which database sequence {other} also holds with the same definition lines and letters but other taxids (NCBI BLAST+ 2.17.0 can write that sequence's taxids and their names for it)"
            ),
            FoundInstead::CopyApart(other, CopyApart::TaxonomyTexts) => format!(
                "has the Seq-align id {id}, which database sequence {other} also holds with the same definition lines and letters: where NCBI BLAST+ 2.17.0 shows that sequence for it, the volume of that sequence looks its taxids up, and the taxonomy database writes a text at such a lookup"
            ),
            FoundInstead::Sequence(other) => format!(
                "has the Seq-align id {id}, by which NCBI BLAST+'s lookup gives database sequence {other} first (NCBI BLAST+ 2.17.0 writes that sequence's names and letters, or stops)"
            ),
            FoundInstead::Failed(failure) => format!(
                "has the Seq-align id {id}, whose lookup in the database fails ({failure}; {})",
                failure.what_ncbi_does(refused)
            ),
            FoundInstead::ThrowsOnId(exception) => format!(
                "has the Seq-align id {id}, whose lookup in the database throws on the id itself ({}; {})",
                String::from_utf8_lossy(exception.get_msg()),
                throws_on_id_does(refused)
            ),
            FoundInstead::Nothing => format!(
                "has the Seq-align id {id}, by which NCBI BLAST+'s formatters do not find the sequence (NCBI BLAST+ 2.17.0 writes Unknown or another sequence's names, or stops)"
            ),
        })
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:770-786
    // The prefetch loads each subject's Bioseq by its Seq-align id, as the formatters find it
    // (`fetched`), by the exact lookup only (the weak matches come where it is off): an id that
    // the lookup finds nowhere loads nothing, without a message (probe of DBID P1 review 1:
    // `gnl|PIG|124`, `gnl|BL_ORD_ID|40` of no sequence); a lookup that fails (a `gnl|BL_ORD_ID|0`
    // stored with a text tag) ends the prefetch with that warning; `lcl|555` that an earlier
    // sequence holds loads that sequence, whose taxids it looks up.
    // The loader reads by the Seq-align id what its lookup gives, also for a subject whose
    // definition line the formatters do not find (`unfetchable`): its own blob where the lookup
    // gives it, whose taxids it looks up (review 2 of DBID P1: `lcl|777` before `gi|777`,
    // `bad_magic`), another sequence's where it gives that one first (review 2 of line DBID
    // pssmread, H3: `gi|777|gb|YYY00001.1|` of taxid 99990 before `lcl|777 ^A gi|777`).
    /// What NCBI's prefetch reads for this subject where no report shows it (`Prefetched`).
    #[cfg(test)]
    pub fn prefetched(&self) -> Prefetched {
        self.prefetched_with_scope(false)
    }

    /// `prefetched` with `seen_by_scope`: where every reporting query's scope finds a noted query
    /// Bioseq for this subject, NCBI's prefetch finds the query in scope and loads nothing from
    /// the BLAST database data loader (`Prefetched::Nothing`).
    pub fn prefetched_with_scope(&self, seen_by_scope: bool) -> Prefetched {
        if seen_by_scope {
            return Prefetched::Nothing;
        }
        if let Some(reason) = self.handle_load_refusal() {
            return Prefetched::Failed(reason);
        }
        let reason = || {
            self.resolution
                .found
                .as_ref()
                .and_then(|found| self.found_reason(found, false, None))
                .unwrap_or_default()
        };
        match &self.resolution.found {
            // The loader does not serve the sequence that the lookup gives: nothing is loaded.
            None
            | Some(FoundInstead::OtherTaxids(_) | FoundInstead::CopyApart(_, CopyApart::Title))
                if self.scope_finds_none =>
            {
                Prefetched::Nothing
            }
            None
            | Some(FoundInstead::OtherTaxids(_) | FoundInstead::CopyApart(_, CopyApart::Title)) => {
                Prefetched::Itself
            }
            Some(FoundInstead::Nothing) => Prefetched::Nothing,
            // The prefetch catches the exception of a failed lookup and writes its warning, or
            // rethrows a database error (`CSeqDBException`) and stops; the formats that read the
            // Bioseq stop (`fetched`).
            Some(FoundInstead::ThrowsOnId(exception)) => Prefetched::Failed(format!(
                "has the Seq-align id {}, whose lookup in the database throws on the id itself ({}; NCBI BLAST+ 2.17.0's prefetch writes Error pre-fetching sequence data, and its later loads of the batch's subjects are made one by one)",
                String::from_utf8_lossy(&self.alignment_fasta()),
                String::from_utf8_lossy(exception.get_msg())
            )),
            Some(FoundInstead::Failed(failure)) => Prefetched::Failed(format!(
                "has the Seq-align id {}, whose lookup in the database fails ({failure}; {})",
                String::from_utf8_lossy(&self.alignment_fasta()),
                failure.what_ncbi_does(None)
            )),
            Some(FoundInstead::Sequence(other)) => Prefetched::Another {
                other: Some(*other),
                reason: reason(),
            },
            Some(FoundInstead::CopyApart(_, CopyApart::TaxonomyTexts)) => Prefetched::Another {
                other: None,
                reason: reason(),
            },
        }
    }

    /// Whether the lookup of this subject's Seq-align id throws on the id where NCBI's prefetch
    /// loads it (`FoundInstead::ThrowsOnId`): the batch's prefetch catches it and writes `Error
    /// pre-fetching sequence data` (`common::database::batch_prefetch`). A lookup that fails on a
    /// file stops NCBI (`FoundInstead::Failed`).
    pub fn prefetch_fails(&self) -> bool {
        matches!(self.resolution.found, Some(FoundInstead::ThrowsOnId(_)))
    }

    /// NCBI's message where the lookup of this subject's Seq-align id fails on a file of the
    /// database with NCBI's own `CSeqDBException` (`FoundInstead::Failed`,
    /// `LookupFailure::ncbi_message`; "DB contains no accession info." of a list that puts
    /// version 5 databases with and without the accession table side by side): NCBI's prefetch of
    /// a batch that reports the subject rethrows it, which stops the run with `BLAST Database
    /// error: <message>` (exit status 2) before the batch's reports
    /// (`common::database::batch_prefetch`).
    pub fn prefetch_stops(&self) -> Option<&[u8]> {
        match &self.resolution.found {
            Some(FoundInstead::Failed(failure)) => failure.ncbi_message(),
            _ => None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2250-2255
    // NCBI's formatters generate the title of the Bioseq they fetch (`CDeflineGenerator`, from
    // the sources of that sequence's volume) for the outfmt 0 description (showdefline.cpp:498)
    // and heading, for the first word that replaces an id in the tabular formats
    // (showdefline.cpp:227) and for the XML2 and JSON `title` (blastxml2_format.cpp:405); the
    // tabular `stitle` and `salltitles` (tabular.cpp:415-434) and the XML `Hit_def`
    // (`GetBioseqHandleDeflineAndId`, showdefline.cpp:283-299) write the titles of the
    // definition lines.
    /// `Err` for a subject found as the same record with other taxids from which NCBI's
    /// formatters can make the title they generate for it (`FoundInstead::CopyApart` with
    /// `CopyApart::Title`), for the reports that show a generated title: outfmt 0
    /// (`check_shown_parsed_subject`) and the tabular word of a replaced id
    /// (`check_printed_word`); outfmt 13 to 16 write the taxids (`fetched`). The callers ask
    /// it where the title they show can depend on the source (`pairwise_title_reads_source`, the
    /// word of a nucleotide or of a protein title ending in `]`).
    pub fn generated_title_fetched(&self) -> Result<(), String> {
        match &self.resolution.found {
            Some(FoundInstead::CopyApart(other, CopyApart::Title)) => Err(format!(
                "has the Seq-align id {}, which database sequence {other} also holds with the same definition lines and letters but other taxids, from which NCBI BLAST+ 2.17.0 can make that sequence's title for it (a title made from the taxids with the taxonomy database, or the end of a protein title rewritten)",
                String::from_utf8_lossy(&self.alignment_fasta())
            )),
            _ => Ok(()),
        }
    }

    /// The Seq-align id (outfmt 8, 9, 11, 12; `DeflineSet::alignment_id`).
    pub fn alignment_id(&self) -> &ParsedSeqId {
        &self.alignment_id
    }

    // NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:131-141
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:2743-2757 (`GetAccVer`, which
    // `GetId` with `eGetId_ForceAcc` calls: the id itself where it is an accession with a version,
    // else `x_GetAccVer` of the Bioseq's ids)
    // NCBI reference (598d8ae6): c++/src/objmgr/scope.cpp:425-445
    /// The Seq-id that NCBI's SAM formatter (outfmt 17) names the subject by
    /// (`AdjustSeqIdType` with `eGetId_ForceAcc`): the Seq-align id where it is an accession
    /// with a version, else the first such id of the subject Bioseq (its gis passed over), else
    /// its last Textseq-id, else the Seq-align id.
    pub fn sam_id(&self) -> &ParsedSeqId {
        forced_accession(&self.alignment_id, &self.target().ids)
    }

    /// The label of the subject's rows in a query-anchored report (outfmt 1 to 4): the Bioseq's
    /// Seq-ids as `x_FillSeqid` reads them, with `-show_gis` (`report::query_anchored::fill_seqid`).
    pub(crate) fn anchored_label(&self, out: &mut Vec<u8>) {
        crate::report::query_anchored::fill_seqid(&self.target().ids, self.show_gis, out);
    }

    /// The definition lines.
    pub fn deflines(&self) -> &[BlastDefline] {
        &self.set.deflines
    }

    /// The definition line of the subject Bioseq.
    fn target(&self) -> &BlastDefline {
        &self.set.deflines[self.target]
    }

    /// What `GenerateDefline` reads of the subject Bioseq: its Seq-ids, whether it has a source
    /// and the name of its first source (none where LOSAT refuses it, `source`).
    fn bioseq_facts(&self) -> BioseqFacts {
        let mut facts = bioseq_facts(&self.target().ids);
        facts.has_source = self.has_source();
        facts.taxname = self.source.name().ok().flatten();
        facts
    }

    /// `GenerateDefline` of the subject Bioseq: for the alignment heading
    /// (`leave_prefix_suffix` false) or the description table (true).
    pub fn bioseq_title(
        &self,
        protein: bool,
        leave_prefix_suffix: bool,
    ) -> Result<Defline, UnknownEncoding> {
        let title = self.target().title.as_deref().unwrap_or_default();
        generate_bioseq_defline(title, protein, leave_prefix_suffix, &self.bioseq_facts())
    }

    /// Whether the subject Bioseq has a source: some definition line has a positive taxid
    /// (`CSeqDBVol::x_GetTaxonomy` makes a source of each, seqdbvol.cpp:1136).
    fn has_source(&self) -> bool {
        self.first_source_taxid().is_some()
    }

    /// The lookups of the volume's taxonomy cache when NCBI fetches the subject Bioseq by its
    /// Seq-align id: the positive taxids of the definition lines in the order of
    /// `CSeqDBVol::x_GetTaxonomy` (`x_GetTaxDefline` puts the line that holds the id first,
    /// `DeflineSet::taxonomy_lookups`), a taxid of several lines once for each.
    pub fn volume_lookups(&self) -> Vec<i32> {
        self.set.taxonomy_lookups(Some(&*self.loader_id()))
    }

    /// The positive taxids that the volume looks up when it fetches the subject Bioseq, each once,
    /// in the order of their first lookup ([`ParsedSubject::volume_lookups`]).
    pub fn volume_taxids(&self) -> Vec<i32> {
        let mut taxids = Vec::new();
        for taxid in self.volume_lookups() {
            if !taxids.contains(&taxid) {
                taxids.push(taxid);
            }
        }
        taxids
    }

    /// Whether NCBI makes the subject Bioseq's title from its source (`bioseq_title_from_source`)
    /// and that source has a taxid: the title then depends on the taxonomy database.
    pub fn title_from_taxonomy(&self) -> bool {
        let title = self.target().title.as_deref().unwrap_or_default();
        bioseq_title_from_source(title, || self.bioseq_facts()) && self.has_source()
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3964-3972
    // `x_AdjustProteinTitleSuffix` (`report::defline::adjust_protein_title_suffix`) rewrites the
    // end of a protein title ending in `]` from the Bioseq's first source, also without a name
    // (a ` [ ]` end goes): `alpha, partial [Bacillus sp.]` is written `alpha [Bacillus sp.]`
    // (probe of line RG T5c, `T/t5c/out/ncbi/ppar.blastp.f0.tx`).
    /// Whether the subject Bioseq's title can lose its end to `x_AdjustProteinTitleSuffix`
    /// (`report::defline::protein_title_suffix_adjustable` of the target's title).
    fn title_ends_in_bracket(&self) -> bool {
        let title = self.target().title.as_deref().unwrap_or_default();
        crate::report::defline::protein_title_suffix_adjustable(title)
    }

    /// Whether some definition line has an id that `seq_id_list` replaces by the first word of
    /// the generated title (`replaces_an_id`; `all_tabular_ids`, `lookup_failure`).
    fn has_replaced_id(&self) -> bool {
        self.set
            .deflines
            .iter()
            .any(|defline| replaces_an_id(&defline.ids))
    }

    /// Whether a report shows the subject Bioseq's generated title (flags 0) where
    /// `x_AdjustProteinTitleSuffix` can rewrite its end from the source: a protein with a source
    /// whose title ends in `]`, with a line without a title (the outfmt 0 heading, outfmt 18) or
    /// an id that the title's first word replaces (the tabular ids). `with_source` names the
    /// source then.
    fn adjusted_title_shown(&self, protein: bool) -> bool {
        protein
            && self.title_ends_in_bracket()
            && self.has_source()
            && (self.heading_generates_title() || self.has_replaced_id())
    }

    /// Whether a title that the reports write of the subject can depend on its Bioseq's
    /// source, the taxids of its definition lines, whatever they are (`title_from_taxonomy`,
    /// `adjusted_title_shown`): in outfmt 0 (`pairwise_title_reads_source`), and the tabular word
    /// of an id that `seq_id_list` replaces, of a protein title ending in `]` (a title made from
    /// the source with the taxonomy database is also shown in outfmt 0; of a protein the word is
    /// `unnamed` either way).
    pub fn title_reads_source(&self, protein: bool, taxdb_found: bool) -> bool {
        self.pairwise_title_reads_source(protein, taxdb_found)
            || (protein && self.title_ends_in_bracket() && self.has_replaced_id())
    }

    /// Whether a title that outfmt 0 shows of the subject can depend on its Bioseq's source
    /// (`title_reads_source`): a title made from the source where `taxdb.bti` is found
    /// (`taxdb_found`; without it NCBI 2.17.0 writes the title of a Bioseq without a source,
    /// probe of line DBID N1b, `blastdb-ids-20261021/R/copytaxid/p1`) in the description and the
    /// heading, and a protein title ending in `]` in the heading of a set with a definition line
    /// without a title, which shows the generated title without `fLeavePrefixSuffix` (the
    /// description keeps the stored end: probe of line DBID N1b L1,
    /// `blastdb-ids-20261021/R/titlesrc/p1`).
    pub fn pairwise_title_reads_source(&self, protein: bool, taxdb_found: bool) -> bool {
        let title = self.target().title.as_deref().unwrap_or_default();
        (taxdb_found && bioseq_title_from_source(title, || self.bioseq_facts()))
            || (protein && self.title_ends_in_bracket() && self.heading_generates_title())
    }

    /// Whether the outfmt 0 alignment heading shows the subject Bioseq's generated title: some
    /// definition line has no title (`heading`).
    pub fn heading_generates_title(&self) -> bool {
        self.set
            .deflines
            .iter()
            .any(|defline| defline.title.as_deref().is_none_or(<[u8]>::is_empty))
    }

    /// The first word of the subject Bioseq's generated title, which `seq_id_list` puts in the
    /// place of a `gnl|BL_ORD_ID` or `lcl|Subject_` id of the Bioseq; empty when no such id is
    /// there (the title is then not made). `Err` where NCBI's `x_CleanAndCompress` reads past the
    /// end of the title (NCBI crashes) and the other refusals of `generated_title_word`.
    /// `Ok(None)` where `GenerateDefline` throws for that title (a title in no encoding that
    /// `CUtf8::GuessEncoding` recognises): `GetSeqIdList` makes no ids (`lookup_failure`).
    fn title_word(&self, protein: bool) -> Result<Option<Vec<u8>>, &'static str> {
        if !replaces_an_id(&self.target().ids) {
            return Ok(Some(Vec::new()));
        }
        self.generated_title_word(protein)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:222-230
    // The first word of a title that NCBI makes from the source with the taxonomy database is
    // the organism's for a nucleotide (`Homo` for an empty title with taxid 9606) and `unnamed`
    // for a protein, as without a source; of a protein title ending in `]`, the word before the
    // end that `x_AdjustProteinTitleSuffix` rewrites where the Bioseq has a source (probe of line
    // DBID N1b L1, `blastdb-ids-20261021/R/titlesrc/p1`). NCBI prints of the ids that
    // `GetSeqIdList` makes the GI and the best one (`GetSeqIdListString`, showdefline.cpp:155-187),
    // so a replaced id beside a better one is not printed.
    /// `Err` where a tabular report prints (`printed`) an id that the first word of a title
    /// NCBI can make from the source replaces: a title made from the taxids of the same record
    /// found for the subject (`generated_title_fetched`) for a nucleotide or a protein title
    /// ending in `]`. The subject's own title made from its source is made with its name
    /// (`source`); the other refusals of the word are `tabular_ids`' and `all_tabular_ids`'.
    pub fn check_printed_word(
        &self,
        protein: bool,
        printed: super::tabular::PrintedIds,
    ) -> Result<(), &'static str> {
        let problem = if self.generated_title_fetched().is_err()
            && (!protein || self.title_ends_in_bracket())
        {
            "a Seq-align id that another database sequence also holds with the same definition lines and letters but other taxids, from which NCBI BLAST+ 2.17.0 can make the title whose first word replaces a printed id (a title made from the taxids with the taxonomy database, or the end of a protein title rewritten)"
        } else {
            return Ok(());
        };
        // Whether the printed ids of `ids` change with the word.
        let shows_word = |ids: &[ParsedSeqId], kinds: &[TabularId]| {
            let (one, other) = (seq_id_list(ids, b"x"), seq_id_list(ids, b"y"));
            kinds.iter().any(|&kind| {
                let (mut a, mut b) = (Vec::new(), Vec::new());
                tabular_id(&one, kind, &mut a);
                tabular_id(&other, kind, &mut b);
                a != b
            })
        };
        let subject = printed.subject
            && shows_word(
                &self.target().ids,
                &[
                    TabularId::Full,
                    TabularId::Accession,
                    TabularId::AccessionVersion,
                ],
            );
        let all =
            printed.all
                && !self.record
                && self.set.deflines.iter().any(|defline| {
                    shows_word(&defline.ids, &[TabularId::Full, TabularId::Accession])
                });
        if subject || all {
            return Err(problem);
        }
        Ok(())
    }

    /// The first word of the subject Bioseq's generated title (`title_word` without the check
    /// of the target's ids): `Err` also where LOSAT refuses the source of a title made from it
    /// (`SourceTaxname::Refused`); `Ok(None)` where `GenerateDefline` throws.
    fn generated_title_word(&self, protein: bool) -> Result<Option<Vec<u8>>, &'static str> {
        self.source.name()?;
        let Ok(defline) = self.bioseq_title(protein, false) else {
            return Ok(None);
        };
        // A nucleotide title stops at the end of the string (approved exception 2 of
        // PD-LOSAT-NCBI-DEFECTS for BLASTN, TBLASTN and TBLASTX).
        if protein && defline.reads_past_end {
            return Err(concat!(
                "a title ",
                crate::report::defline::punctuation_title!()
            ));
        }
        let end = defline
            .text
            .iter()
            .position(|&byte| byte == b' ')
            .unwrap_or(defline.text.len());
        Ok(Some(defline.text[..end].to_vec()))
    }

    /// The tabular `sseqid`, `sacc` and `saccver` (`tabular_id` of the subject Bioseq's ids after
    /// `seq_id_list`); `Ok(None)` where `SetSubjectId`'s `GetSeqIdList` throws (`lookup_failure`).
    pub fn tabular_ids(&self, protein: bool) -> Result<Option<[Vec<u8>; 3]>, &'static str> {
        let Some(word) = self.title_word(protein)? else {
            return Ok(None);
        };
        let ids = seq_id_list(&self.target().ids, &word);
        let field = |kind| {
            let mut out = Vec::new();
            tabular_id(&ids, kind, &mut out);
            out
        };
        Ok(Some([
            field(TabularId::Full),
            field(TabularId::Accession),
            field(TabularId::AccessionVersion),
        ]))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:862-869
    /// Which of NCBI's lookups of the subject's tabular ids throw (`report::tabular::FailedLookup`):
    /// `SetSubjectId` where the subject Bioseq's ids hold one that `GetSeqIdList` replaces
    /// (`lcl|Subject_`, `gnl|BL_ORD_ID`) and `GenerateDefline` throws, `x_SetSubjectIds` where a
    /// definition line holds one (a `-subject` record's ids are not looked up there), with
    /// whether a definition line has no Seq-id (`LookupFailure::lineless`, at which `sallacc`
    /// stops); `None` where neither throws and every line has a Seq-id. `Err` as `tabular_ids`'.
    pub fn lookup_failure(
        &self,
        protein: bool,
    ) -> Result<Option<super::tabular::LookupFailure>, &'static str> {
        let subject_id = replaces_an_id(&self.target().ids);
        let all_ids = !self.record && self.has_replaced_id();
        let lineless = self.has_lineless_line();
        let failed = (subject_id || all_ids) && self.generated_title_word(protein)?.is_none();
        Ok(
            (failed || lineless).then_some(super::tabular::LookupFailure {
                subject_id: failed && subject_id,
                all_ids: failed && all_ids,
                bioseq: false,
                length: false,
                lineless,
            }),
        )
    }

    /// The Seq-align id as `sallseqid` and `sallacc` print it (`tabular_id` of the id alone,
    /// `eFullId` and `eAccession`): what the catch of `SetFields` pushes where a lookup throws
    /// (`report::tabular::FailedLookup::full`, `accession`).
    pub fn alignment_tabular_ids(&self) -> [Vec<u8>; 2] {
        let ids = std::slice::from_ref(&self.alignment_id);
        let field = |kind| {
            let mut out = Vec::new();
            tabular_id(ids, kind, &mut out);
            out
        };
        [field(TabularId::Full), field(TabularId::Accession)]
    }

    /// The Seq-align id as `sallgi` prints it (`FindGi` of the id alone): its GI, 0 for an id that
    /// is none.
    pub fn alignment_gi(&self) -> i64 {
        self.alignment_id.gi().unwrap_or(0)
    }

    /// The keys of the ids of queries that NCBI's scope would show for this sequence
    /// (`QueryScopeKeys`), made once.
    pub fn query_scope_keys(&self) -> &QueryScopeKeys {
        self.query_scope_keys.get_or_init(|| {
            let mut own = Vec::new();
            scope_keys(self, &mut own);
            let mut weak = Vec::new();
            crate::seqid::id_scope::weak_keys(&self.alignment_id, &mut weak);
            let mut display = Vec::new();
            if let Some(id) = self.display_id() {
                crate::seqid::defline_set::scope_keys(id, &mut display);
                crate::seqid::id_scope::weak_keys(id, &mut display);
            }
            QueryScopeKeys { own, weak, display }
        })
    }

    /// The `saccver` of `tabular_ids`, made once for a subject (`protein` is its database's);
    /// `None` where the lookup throws, as `tabular_ids` (the programs' id lists,
    /// `common::tabular::SubjectId`).
    pub fn tabular_saccver(
        &self,
        protein: bool,
    ) -> Result<Option<std::sync::Arc<[u8]>>, &'static str> {
        if let Some(saccver) = self.saccver.get() {
            return Ok(saccver.clone());
        }
        let saccver = self
            .tabular_ids(protein)?
            .map(|[_, _, saccver]| std::sync::Arc::from(saccver));
        Ok(self.saccver.get_or_init(|| saccver).clone())
    }

    /// `Err` as `tabular_ids`', which it does not make.
    pub fn check_tabular_ids(&self, protein: bool) -> Result<(), &'static str> {
        self.title_word(protein).map(drop)
    }

    /// The tabular `sgi`: the first GI of the subject Bioseq's ids (`FindGi`; `seq_id_list` only
    /// replaces ids that are no GI), 0 for none.
    pub fn tabular_gi(&self) -> i64 {
        self.target().gi()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:544-556
    /// The tabular ids of every definition line (`x_SetSubjectIds`: `sallseqid`, `sallgi`,
    /// `sallacc`), in the stored order: each line's ids after `seq_id_list` with the first word
    /// of the subject Bioseq's title, as `tabular_ids` makes them. Where `GenerateDefline` throws,
    /// the ids of the lines before the first line with an id to replace, which `x_SetSubjectIds`
    /// pushed before it threw (`AllTabularIds::failed`). `Err` as `tabular_ids`'.
    pub fn all_tabular_ids(&self, protein: bool) -> Result<AllTabularIds, &'static str> {
        // `GetSeqIdList` replaces the ids of each definition line, so the word is made when any
        // line has an id to replace, not only the target's.
        let replaced = !self.record && self.has_replaced_id();
        let word = if replaced {
            self.generated_title_word(protein)?
        } else {
            Some(Vec::new())
        };
        let Some(word) = word else {
            let before = self
                .set
                .deflines
                .iter()
                .take_while(|defline| !replaces_an_id(&defline.ids))
                .map(|defline| defline_tabular_ids(defline, b""))
                .collect();
            return Ok(AllTabularIds {
                ids: before,
                failed: true,
            });
        };
        Ok(AllTabularIds {
            ids: self
                .set
                .deflines
                .iter()
                .map(|defline| defline_tabular_ids(defline, &word))
                .collect(),
            failed: false,
        })
    }

    /// The title of the first definition line, empty for none (`stitle`; the writers show an
    /// empty one as `N/A`).
    pub fn first_title(&self) -> &[u8] {
        self.set
            .deflines
            .first()
            .and_then(|defline| defline.title.as_deref())
            .unwrap_or_default()
    }

    /// Whether a stored title of its definition lines held a byte outside ' ' ..= '~', read as
    /// '#' (`algorithm::common::database::check_reported_title`).
    pub fn replaced_title_bytes(&self) -> bool {
        self.set
            .deflines
            .iter()
            .any(|defline| defline.replaced_title_bytes > 0)
    }

    /// The start of the outfmt 0 description line (`description_id` of `GetDisplayIds`' id) and
    /// its title (`description_line`).
    pub fn description(&self, protein: bool) -> (Vec<u8>, Defline) {
        let (start, title, _) = self.description_line(protein);
        (start, title)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1585-1613
    /// The outfmt 0 description line: the start (`description_id` of `GetDisplayIds`' id) and the
    /// title, `GenerateDefline` of the subject Bioseq with `fLeavePrefixSuffix`, after two spaces
    /// where `gi|N|` was written; and whether it is the line of NCBI's `catch`
    /// (`x_GetDeflineInfo`): where the title throws (NCBI's `HtmlDecode` cannot decode it), the
    /// `catch` keeps the Seq-align id and `Unknown` (`thrown_description`; line RG T4a2 review O1:
    /// `>gi|777|ref|XP_000009.1| ...` with `-subject -parse_deflines` shows `gi|777 Unknown`).
    ///
    /// The `catch` also takes a subject whose walk of the definition lines meets one without a
    /// Seq-id (`description_meets_lineless_line`).
    pub fn description_line(&self, protein: bool) -> (Vec<u8>, Defline, bool) {
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:388-394
        // The scope's Bioseq handle is empty (`scope_finds_none`): `x_FillDeflineAndId` returns at
        // once, which leaves the Seq-align id and `Unknown` as the `catch` does (without the
        // `catch`'s links, `HtmlReport::description_links`; its GI is never set,
        // `description_gi_unset`).
        let mut title = if self.scope_finds_none {
            let (start, title) = self.thrown_description();
            return (start, title, true);
        } else {
            match self.bioseq_title(protein, true) {
                Ok(title) if !self.description_meets_lineless_line(&title.text) => title,
                _ => {
                    let (start, title) = self.thrown_description();
                    return (start, title, true);
                }
            }
        };
        let (start, wrote_gi) = self.description_id();
        if wrote_gi {
            title.text.splice(0..0, *b"  ");
        }
        (start, title, false)
    }

    /// Whether some definition line has no Seq-id (a hand-edited database; makeblastdb writes
    /// none such). NCBI's formatters take the best id of its empty list, a null reference, where
    /// they walk the lines: each format's `catch` decides what it writes (`description_line`,
    /// `heading`, `xml_hit`, `SubjectNames::of_parsed`, `lookup_failure`, `tax_deflines`).
    pub fn has_lineless_line(&self) -> bool {
        self.set
            .deflines
            .iter()
            .any(|defline| defline.ids.is_empty())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1585-1604
    /// Whether `description_line` can take NCBI's `catch` branch (whose start is the Seq-align
    /// id's, not the display id's): a definition line without a Seq-id, or a non-ASCII byte in the
    /// target's title, source name, or empty-title PDB/patent id fields (the only inputs for which
    /// `html_decode` can return `UnknownEncoding`).
    pub fn description_may_throw(&self) -> bool {
        if self.has_lineless_line() {
            return true;
        }
        let title = self.target().title.as_deref().unwrap_or_default();
        if !title.is_ascii() {
            return true;
        }
        if !title.iter().all(|byte| b".,;~ ".contains(byte)) {
            return false;
        }
        if self
            .source
            .name()
            .ok()
            .flatten()
            .is_some_and(|name| !name.is_ascii())
        {
            return true;
        }
        self.target().ids.iter().any(|id| match id {
            ParsedSeqId::Pdb(pdb) => pdb.chain_id.as_deref().is_some_and(|s| !s.is_ascii()),
            ParsedSeqId::Patent(patent) => {
                let number = match &patent.cit.id {
                    PatentNumber::Number(number) | PatentNumber::AppNumber(number) => number,
                };
                !patent.cit.country.is_ascii() || !number.is_ascii()
            }
            _ => false,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4748-4760
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:497-536
    // `GetDisplayIds` takes the label of each line's best id up to the line that holds the
    // Seq-align id, and the walk after the generated title takes it of every line before it looks
    // at the title: a line without a Seq-id throws in either (its best id is null). The walk adds
    // ` >`, `gi|N|` with `-show_gis` and a GI, the FASTA label of the line's `WorstRank`-best id,
    // a space and the title for each titled line (an empty title too) whose ids do not match the
    // display id, and only then stops once the text passes `kMaxDescrLength` (4096 bytes,
    // showdefline.cpp:145): a line after that stop is not reached.
    /// Whether the outfmt 0 description of the subject, whose description title (`GenerateDefline`
    /// with `fLeavePrefixSuffix`) is `title`, reaches a definition line without a Seq-id: one
    /// before the line of the display id (`display_defline`), or one that the walk of the lines
    /// meets before its text passes 4096 bytes. NCBI's `catch` then writes the line
    /// (`thrown_description`).
    fn description_meets_lineless_line(&self, title: &[u8]) -> bool {
        if !self.has_lineless_line() {
            return false;
        }
        let Some(shown) = self.set.display_defline(&self.alignment_id) else {
            return false;
        };
        if self.set.deflines[..shown.index]
            .iter()
            .any(|defline| defline.ids.is_empty())
        {
            return true;
        }
        let display = self.set.deflines[shown.index].display_id();
        let mut length = title.len();
        for defline in &self.set.deflines {
            if defline.ids.is_empty() {
                return true;
            }
            let Some(text) = &defline.title else {
                continue;
            };
            if display.is_some_and(|display| {
                defline
                    .ids
                    .iter()
                    .any(|id| compare(id, display) == Sic::Yes)
            }) {
                continue;
            }
            let gi = defline.gi();
            if self.show_gis && gi > 0 {
                length += b"gi|".len() + gi.to_string().len() + 1;
            }
            let mut label = Vec::new();
            if let Some(best) = defline.display_id() {
                label::fasta(best, &mut label);
            }
            length += b" >".len() + label.len() + 1 + text.len();
            if length > MAX_DESCRIPTION_LENGTH {
                return false;
            }
        }
        false
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1576-1604
    /// The outfmt 0 description line (`description`) of a subject whose lookup throws there
    /// (`HitLoad::Throws`): the `catch` keeps the Seq-align id and the title `Unknown`, and the
    /// GI of a GI id; the line's start is made of them as of any other (`description_id`: the bare
    /// id and a space, or with `-show_gis` the FASTA form of an id that is not a GI, which no space
    /// follows: `ref|NP_100001.1|Unknown`). A GI Seq-align id goes the same way (`description_id`,
    /// which the probes check with the loaded lines and the thrown accession, gnl, lcl and giim
    /// ones), though no database makes it throw: NCBI finds a GI by its ISAM index, never by the
    /// accession table whose absence throws (review 2 of pfoff-3, L3).
    pub fn thrown_description(&self) -> (Vec<u8>, Defline) {
        let mut start = Vec::new();
        let gi = self.alignment_id.gi().unwrap_or(0);
        let wrote_gi = description_id(&self.alignment_id, gi, self.show_gis, &mut start);
        let mut title = b"Unknown".to_vec();
        if wrote_gi {
            title.splice(0..0, *b"  ");
        }
        (
            start,
            Defline {
                text: title,
                reads_past_end: false,
            },
        )
    }

    /// Whether the Seq-align id is a local id, which `CAlignFormatUtil::GetDbType` does not look
    /// up (`AfterThrow::alignment_end`).
    pub fn alignment_id_is_local(&self) -> bool {
        matches!(self.alignment_id, ParsedSeqId::Local(_))
    }

    /// The start of the outfmt 0 description line (`description`) where its title does not make it
    /// the line of NCBI's `catch` (`description_line`).
    pub fn description_start(&self) -> Vec<u8> {
        self.description_id().0
    }

    /// `description`'s start and whether `gi|N|` was written.
    fn description_id(&self) -> (Vec<u8>, bool) {
        let mut start = Vec::new();
        let mut wrote_gi = false;
        if let Some((id, gi)) = self.description_ids() {
            wrote_gi = description_id(id, gi, self.show_gis, &mut start);
        }
        (start, wrote_gi)
    }

    /// The display id of the description line and its GI (`GetDisplayIds`' `wid` and `gi`: the
    /// GI of `display_defline` where it holds the Seq-align id, else 0).
    pub fn description_ids(&self) -> Option<(&ParsedSeqId, i64)> {
        let (defline, found) = self.display_defline()?;
        let id = defline.display_id()?;
        Some((id, if found { defline.gi() } else { 0 }))
    }

    /// The ids of the subject Bioseq (`handle.GetBioseqCore()->GetId()` of the Seq-align id: its
    /// definition line's, a `-subject` record's own), whose link `-html` writes
    /// (`report::html_ids::id_url`).
    pub fn bioseq_ids(&self) -> &[ParsedSeqId] {
        &self.target().ids
    }

    /// The taxid of the first definition line that has one (0 for none; a `-subject` record has
    /// no Blast-def-line), the `<@taxid@>` of a local id's `-html` link (showdefline.cpp:484-491).
    pub fn first_taxid(&self) -> i32 {
        self.set
            .deflines
            .iter()
            .find_map(|defline| defline.taxid)
            .unwrap_or(0)
    }

    /// The definition line whose id the description line shows (`DeflineSet::display_defline`
    /// of the Seq-align id) and whether it holds the Seq-align id.
    fn display_defline(&self) -> Option<(&BlastDefline, bool)> {
        let shown = self.set.display_defline(&self.alignment_id)?;
        Some((&self.set.deflines[shown.index], shown.found))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:432-433, 497-498
    /// The display id: the id of the description line (`GetDisplayIds`' `wid`, the `WorstRank`
    /// best id of `display_defline`), by which the description table (outfmt 0 to 4) looks the
    /// title up again in NCBI's scope, which finds a query of the batch with that id before the
    /// database sequence (`common::database::refuse_display_id_of_a_query`).
    pub fn display_id(&self) -> Option<&ParsedSeqId> {
        self.display_defline()?.0.display_id()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2380-2384,2451-2471
    // `out` holds the heading of this subject only; the title of a definition line is its own,
    // or the subject Bioseq's generated title when it is empty (`x_FillDispParams`).
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3166-3183
    /// The definition-line lines of the outfmt 0 alignment of the subject (each ending with a
    /// newline; `Length=` follows): the definition lines (`heading_shows`), the first after `>`
    /// and the others after a space, with its id (`alignment_line_id`) and its title wrapped at
    /// 60 bytes. `every_line`: the Seq-align has no `ext` when the report is made, so nothing
    /// limits the lines: a Seq-align of an ungapped search with two or more HSPs, which the
    /// formatter splits into new ones (`PrepareBlastUngappedSeqalign` above). The generated title
    /// is made once, for the first line without a title; where it throws, `Err` with the Seq-id
    /// text of NCBI's line in place of each alignment (`heading_skipped_id`).
    ///
    /// With `-html` (`html`): each line's id is its link, the id, `</a>` after a link and its
    /// anchor (`HtmlHeading::line_id`; a `-subject` record's anchor without the version, as
    /// `x_PrintDefLine` writes a Bioseq without definition lines), then always a space and the
    /// HTML-encoded title wrapped; a heading with more than 8 lines (NCBI's `maxNumBdl`: the
    /// definition lines, or the strings of the positive list) hides those after the third in the
    /// "more sequence titles" box (`html_ids::more_titles`, cited there). The id link precedes
    /// the title of each line (`x_FillAlnDispParams`); where the registry throws, `Err` with the
    /// Seq-id text too (`HtmlHeading::skips`). The outer error is the report's.
    pub fn heading(
        &self,
        protein: bool,
        every_line: bool,
        html: Option<HtmlHeading<'_>>,
    ) -> std::io::Result<Result<Vec<u8>, Vec<u8>>> {
        if self.scope_finds_none {
            return Ok(Err(self.skipped_id().to_vec()));
        }
        let mut out = Vec::new();
        let mut generated = None;
        let mut first = true;
        let max_lines = if every_line || self.use_this_seqid.is_empty() {
            self.set.deflines.len()
        } else {
            self.use_this_seqid.len()
        };
        let (mut shown, mut id_label, mut hides) = (0, Vec::new(), false);
        for defline in &self.set.deflines {
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2195-2198
            // `x_FillAlnDispParams` takes the label of each line's best id first, before its
            // link and title: a line without a Seq-id throws there (its best id is null), with
            // or without a limit of the lines (`heading_skipped_id`).
            if defline.ids.is_empty() {
                return Ok(Err(self.skipped_id().to_vec()));
            }
            if !every_line && !self.heading_shows(defline) {
                continue;
            }
            shown += 1;
            let start = out.len();
            out.push(if first { b'>' } else { b' ' });
            match html {
                None => alignment_line_id(defline, self.show_gis, &mut out),
                Some(html) => {
                    if max_lines > html_ids::MAX_SHOWN_LINES && shown == html_ids::SHOWN_LINES + 1 {
                        let hidden = max_lines - html_ids::SHOWN_LINES;
                        html_ids::more_titles(&id_label, hidden, &mut out);
                        hides = true;
                    }
                    if let Some(display) = defline.display_id() {
                        match html.line_id(
                            &defline.ids,
                            display,
                            defline.gi(),
                            defline.taxid.unwrap_or(0),
                            !self.record,
                            |out| alignment_line_id(defline, self.show_gis, out),
                            &mut out,
                        ) {
                            Ok(label) => id_label = label,
                            Err(error) => {
                                return HtmlHeading::skips(error, || self.skipped_id().to_vec())
                                    .map(Err)
                            }
                        }
                    }
                }
            }
            if !first || out.len() - start > 1 {
                out.push(b' ');
            }
            first = false;
            let title = match defline.title.as_deref() {
                Some(title) if !self.heading_generates(title) => title,
                _ => match &generated {
                    Some(text) => text,
                    None => match self.bioseq_title(protein, false) {
                        Ok(defline) => generated.insert(defline.text),
                        Err(UnknownEncoding) => return Ok(Err(self.skipped_id().to_vec())),
                    },
                },
            };
            match html {
                None => wrap_title(title, &mut out),
                Some(_) => wrap_title(&crate::report::html::html_encode(title), &mut out),
            }
            out.push(b'\n');
        }
        if hides {
            out.extend_from_slice(b"</div>");
        }
        Ok(Ok(out))
    }

    /// Whether the heading shows the generated title for a definition line with the title
    /// `title` (`x_FillDispParams`: an empty one; every line of a `-subject` record, whose
    /// title is the record's own, as `heading_skipped_id` has it).
    fn heading_generates(&self, title: &[u8]) -> bool {
        title.is_empty() || self.record
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1976-1981
    // `x_DisplayAlnvecInfo` fills each shown definition line (`x_FillDispParams`) before it
    // writes the heading, and a line without a title gets the Bioseq's generated title, whose
    // `NStr::HtmlDecode` throws where it cannot guess the encoding (a name of the taxonomy
    // database, line RG item T5; a `-subject` record's title, which its heading always generates:
    // probe `T/t5/sub`, `AB999999` for `gb|AB999999.1|`): no heading and no alignment are written,
    // the message instead, for each Seq-align of the subject (`previousId` is not set, so the
    // next one throws too).
    /// The Seq-id text (`CSeq_id::GetSeqIdString` of the Seq-align id) that outfmt 0 writes in
    /// place of the alignments of a subject whose heading (`heading`, with `every_line`) shows the
    /// generated title and NCBI cannot make it (`UnknownEncoding`), or that has a definition line
    /// without a Seq-id (`has_lineless_line`, which every heading meets); `None` otherwise.
    pub fn heading_skipped_id(&self, protein: bool, every_line: bool) -> Option<Vec<u8>> {
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1952-1981
        // `CAlnVec::GetBioseqHandle` throws on a subject whose Bioseq the scope does not get
        // (`scope_finds_none`), before any heading.
        if self.scope_finds_none || self.has_lineless_line() {
            return Some(self.skipped_id().to_vec());
        }
        let generated = self.set.deflines.iter().any(|defline| {
            (every_line || self.heading_shows(defline))
                && self.heading_generates(defline.title.as_deref().unwrap_or_default())
        });
        if !generated || self.bioseq_title(protein, false).is_ok() {
            return None;
        }
        Some(self.skipped_id().to_vec())
    }

    /// `CSeq_id::GetSeqIdString` of the Seq-align id: the id of NCBI's line in place of each
    /// alignment of a subject whose heading it cannot make (`heading_skipped_id`).
    pub fn skipped_id(&self) -> &[u8] {
        self.skipped_id.get_or_init(|| {
            let mut id = Vec::new();
            label::content(&self.alignment_id, false, &mut id);
            id
        })
    }

    /// Whether the `Score`-best id of the subject Bioseq is a local id (XML2 then writes no
    /// accession of a `-subject` record, `SubjectNames::of_believed`).
    pub fn has_local_best_id(&self) -> bool {
        matches!(
            find_best(&self.target().ids, text_score),
            Some(ParsedSeqId::Local(_))
        )
    }

    /// The XML `Hit_id` and `Hit_accession` of a `-subject` record (`of_record`; `Hit_def` is the
    /// generated title, `bioseq_title`).
    pub fn record_xml_ids(&self) -> (Vec<u8>, Vec<u8>) {
        let mut id = Vec::new();
        label::seqid_list_string(&self.target().ids, true, &mut id);
        let mut accession = Vec::new();
        if let Some(best) = find_best(&self.target().ids, text_score) {
            label::format_label(best, false, &mut accession);
        }
        (id, accession)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2196-2202
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4648-4664
    /// Whether the outfmt 0 alignment heading shows `defline`: under a positive list only a line
    /// that `use_this_seqid` names (a line of a local id alone is not shown, nor one without a
    /// GI when the strings are GIs; probe G05 of line DI-8).
    fn heading_shows(&self, defline: &BlastDefline) -> bool {
        if self.use_this_seqid.is_empty() {
            return true;
        }
        // `CAlignFormatUtil::GetLabel`: the tag of a `general` id, else `GetSeqIdString` with true.
        let mut label_of_id = Vec::new();
        match defline.display_id() {
            Some(id @ ParsedSeqId::General(_)) => {
                label::content_general_db_is_content(id, &mut label_of_id);
            }
            Some(id) => label::content(id, true, &mut label_of_id),
            None => {}
        }
        // `NStr::StringToNumeric` for a `TGi`, of `useThisSeq` (`NStr::StringToInt8`, flags 0), of
        // a GI that `under_positive_list` wrote in decimal.
        let gi = defline.gi();
        self.use_this_seqid
            .iter()
            .any(|string| match string.strip_prefix(b"gi:") {
                Some(value) => crate::ncbi_str::string_to_int8(value, 0).is_ok_and(|n| n == gi),
                None => string.strip_prefix(b"seqid:") == Some(label_of_id.as_slice()),
            })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:281-298
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:399-407
    // `GetBioseqHandleDeflineAndId` takes the id list of each line while no line had a title
    // (an empty list gives an empty string) and the label of the best id of each titled line
    // after the first titled one: a line without a Seq-id throws only there (its best id is
    // null), and the `Hit`'s `catch` writes the Seq-align id, `Unknown`, `Unknown` and 0.
    /// The outfmt 5 `Hit_id` and `Hit_def` (`DeflineSet::hit_id_and_def` with the GIs shown;
    /// an empty `Hit_def` is written `No definition line`) and `Hit_accession` (the label of
    /// the `Score`-best id of the subject Bioseq, without its version; a `general` id's tag).
    /// `None` where a definition line without a Seq-id has a title (an empty one too) and comes
    /// after the first line with a title: the `Hit` is NCBI's `catch`'s
    /// (`blast_xml::SubjectNames::of_parsed`).
    pub fn xml_hit(&self) -> Option<(Vec<u8>, Vec<u8>, Vec<u8>)> {
        let mut titled = false;
        for defline in &self.set.deflines {
            if defline.title.is_none() {
                continue;
            }
            if titled && defline.ids.is_empty() {
                return None;
            }
            titled = true;
        }
        let (id, def) = self.set.hit_id_and_def(true);
        let mut accession = Vec::new();
        if let Some(best) = find_best(&self.target().ids, text_score) {
            label::format_label(best, false, &mut accession);
        }
        Some((id, def, accession))
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:365-388
    /// The `HitDescr` of each definition line (XML2 and JSON), without the scientific names; those
    /// of the lines from the first one without a Seq-id on are not written
    /// (`blast_xml::SubjectNames::of_parsed`).
    pub fn hit_descriptions(&self) -> Vec<HitDescription> {
        self.set
            .deflines
            .iter()
            .map(|defline| {
                let mut id = Vec::new();
                label::seqid_list_string(&defline.ids, true, &mut id);
                let mut accession = Vec::new();
                if let Some(best) = defline.display_id() {
                    label::format_label(best, false, &mut accession);
                }
                HitDescription {
                    id,
                    accession,
                    title: defline.title.clone(),
                    taxid: defline.taxid.unwrap_or(0),
                    taxonomy: HitTaxonomy::Unnamed,
                }
            })
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:935-976
    /// The definition lines of the organism report (outfmt 18, `report::tax_report`): of a
    /// database sequence, each line that a positive list keeps (`heading_shows`) with the label of
    /// its `WorstRank`-best id without version, its title (the Bioseq's generated title where it
    /// has none) and its taxid; of a `-subject` record (no `Blast-def-line`), the Bioseq's
    /// `WorstRank`-best id, its generated title and taxid 0 (`TaxDefline::of_record`). A generated
    /// title that LOSAT does not write as NCBI is a `TaxTitle::Refused` (those of the outfmt 0
    /// heading, which shows the same title, `common::database::check_shown_parsed_subject`; the
    /// one made from a source that LOSAT refuses, `SourceTaxname::Refused`, only where the report
    /// shows it). A title made from the source names the organism of the Bioseq's first source
    /// (`with_source`). NCBI BLAST+ 2.17.0 keeps the positive list in an ungapped search too
    /// (probe `T/t2/fix1` `f1.d3.*.sl*`: tblastx and `blastn -ungapped` with `-seqidlist` show the
    /// listed line alone for a subject of 2 to 5 HSPs).
    pub fn tax_deflines(&self, protein: bool) -> Vec<crate::report::tax_report::TaxDefline> {
        use crate::report::tax_report::{generated_title, TaxDefline, TaxTitle};
        let label_of = |id: Option<&ParsedSeqId>| {
            let mut label = Vec::new();
            if let Some(id) = id {
                label::format_label(id, false, &mut label);
            }
            label
        };
        let refused = |reason: &str, shown_only: bool| TaxTitle::Refused {
            reason: reason.to_owned(),
            shown_only,
        };
        let generated = || -> TaxTitle {
            if let Err(reason) = self.source.name() {
                return refused(&format!("has {reason}"), true);
            }
            if let Err(reason) = self.generated_title_fetched() {
                return refused(&reason, false);
            }
            let title = self.target().title.as_deref().unwrap_or_default();
            generated_title(title, protein, &self.bioseq_facts())
        };
        if self.record {
            return vec![TaxDefline {
                label: label_of(find_best(&self.target().ids, text_score)),
                title: generated(),
                taxid: 0,
            }];
        }
        self.set
            .deflines
            .iter()
            .filter(|defline| self.heading_shows(defline))
            .map(|defline| TaxDefline {
                label: label_of(defline.display_id()),
                title: match defline.title.as_deref() {
                    // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:944-947
                    // `x_FillTaxDispParams` takes the label of the line's best id first, null for a
                    // line without a Seq-id: the `catch` of `x_InitTaxInfoMap` drops the rest of the
                    // HSP's subject without marking it shown, so its next HSP shows the lines before
                    // again (Owner question Q-D1c1cb-18 of line RG D1c-1c-b; LOSAT refuses it until
                    // it is answered).
                    _ if defline.ids.is_empty() => refused(LINELESS_TAX_REPORT, false),
                    Some(title) if !title.is_empty() => TaxTitle::Text(title.to_vec()),
                    _ => generated(),
                },
                taxid: defline.taxid.unwrap_or(0),
            })
            .collect()
    }

    /// The FASTA form of the Seq-align id (`CSeq_id::AsFastaString`): the id of the `HitDescr`
    /// that NCBI writes in place of the rest when the taxonomy database does not know a taxid
    /// (blastxml2_format.cpp:422-428), and the id that LOSAT's refusals name.
    pub fn alignment_fasta(&self) -> Vec<u8> {
        let mut out = Vec::new();
        label::fasta(&self.alignment_id, &mut out);
        out
    }
}

/// The Seq-id of the handle of the Seq-align id `id` where it is not `id` in a way that
/// `compare` sees (`ParsedSubject::loader_id`): a PDB id whose `chain` the handle normalizes
/// (`seqid::mask_handle::pdb_handle_seq_id`). A branch on the decoded id: other ids cost nothing.
fn handle_changes(id: &ParsedSeqId) -> Option<ParsedSeqId> {
    let ParsedSeqId::Pdb(pdb) = id else {
        return None;
    };
    let handle = crate::seqid::mask_handle::pdb_handle_seq_id(pdb);
    (handle.chain != pdb.chain).then_some(ParsedSeqId::Pdb(handle))
}

/// What NCBI's data loader makes of a sequence whose Seq-align id's handle has another chain
/// (`handle_lookup`, line RG-KB item 12).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum HandleLookup {
    /// No id of the sequence compares equal to the handle's Seq-id: the loader serves no Bioseq
    /// (`ParsedSubject::scope_finds_none`).
    NoBioseq,
    /// The Bioseq is made of this definition line, which holds the handle's Seq-id
    /// (`DeflineSet::bioseq_defline`), another than the Seq-align id's: its ids and title are the
    /// Bioseq's (`ParsedSubject::target`), which the description line, the titles that the
    /// heading and the organism report generate, and the tables read (probes `e15_twoline` and
    /// `e15nt` of the reviews of line RG-KB item 12: `>1AB1_A t1` for an untitled Seq-align id's
    /// line and a line `pdb|1AB1|A t1`).
    Line(usize),
    /// Another id of the sequence compares equal to it, but no definition line holds it as the
    /// load reads the ids (one definition line with the chains 321 and 65: its first pdb id is not
    /// the handle's): the load throws (`ParsedSubject::handle_load_refusal`).
    Throws,
}

// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:386-431
// NCBI reference (598d8ae6): c++/src/objects/seqloc/PDB_seq_id.cpp:53-67
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1328-1359
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/cached_sequence.cpp:58
/// What NCBI's data loader makes of the Seq-align id `id` of the sequence of the definition lines
/// `set` where its lookup gives that sequence, `None` where it serves the Bioseq of the Seq-align
/// id's line (`own_line`) as for any id: `x_GetOid` looks up the Seq-id of the id's handle
/// (`seqid::mask_handle::pdb_handle_seq_id`: of a PDB id, its normalized `chain`) and keeps the
/// sequence only where one of its ids compares equal to that Seq-id, which the stored id does not
/// where the handle changes its `chain`; the Bioseq is then built for that Seq-id, of the
/// definition line that holds it as the load reads the ids.
fn handle_lookup(set: &DeflineSet, id: &ParsedSeqId, own_line: usize) -> Option<HandleLookup> {
    let handle = handle_changes(id)?;
    if !set.ids().any(|held| compare(&handle, held) == Sic::Yes) {
        return Some(HandleLookup::NoBioseq);
    }
    match set.bioseq_defline(&handle) {
        Some(line) if line == own_line => None,
        Some(line) => Some(HandleLookup::Line(line)),
        None => Some(HandleLookup::Throws),
    }
}

/// The ids of every definition line of a subject (`ParsedSubject::all_tabular_ids`); where
/// `failed`, those of the lines before the one whose `GetSeqIdList` throws.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AllTabularIds {
    pub ids: Vec<super::tabular::DeflineIds>,
    pub failed: bool,
}

/// Whether `GetSeqIdList` replaces one of `ids` with the first word of the generated title
/// (`seq_id_list`: a `gnl|BL_ORD_ID` or `lcl|Subject_` id).
fn replaces_an_id(ids: &[ParsedSeqId]) -> bool {
    matches!(seq_id_list(ids, b"x"), std::borrow::Cow::Owned(_))
}

/// The tabular ids of one definition line after `seq_id_list` with `word`.
fn defline_tabular_ids(defline: &BlastDefline, word: &[u8]) -> super::tabular::DeflineIds {
    let ids = seq_id_list(&defline.ids, word);
    let field = |kind| {
        let mut out = Vec::new();
        tabular_id(&ids, kind, &mut out);
        out
    };
    super::tabular::DeflineIds {
        id: field(TabularId::Full),
        gi: defline.gi(),
        accession: field(TabularId::Accession),
    }
}

/// The subjects of a `-subject` run under `-parse_deflines`, by subject index: `records` gives the
/// ids (`None`: shown as a record without the option) and title of each. LOSAT's refusals of them
/// in outfmt 0 to 4 are `check_pairwise_records` (their titles) and
/// `algorithm::common::database::refuse_described_records` (their description lines).
pub fn parsed_records<'a>(
    records: impl Iterator<
        Item = (
            Option<&'a crate::blastinput::fasta_reader::believed::BelievedIds>,
            &'a [u8],
        ),
    >,
    show_gis: bool,
) -> Vec<Option<std::sync::Arc<ParsedSubject>>> {
    records
        .map(|(ids, title)| {
            ids.map(|ids| {
                std::sync::Arc::new(ParsedSubject::of_record(ids, title).showing_gis(show_gis))
            })
        })
        .collect()
}

/// LOSAT's refusals of outfmt 0 of the subjects of `parsed_records`, in their order: those that
/// outfmt 0 shows by their parsed ids (`ParsedSubject::shown_as_record`) go through
/// `check_pairwise_record`. A call passes the result to the refusal of each of its formats that
/// shows the titles (`FormatRefusal::check`), made right after `parsed_records`; `program` names
/// the program. What a record's description line shows is refused only where a report prints
/// that line (`algorithm::common::database::refuse_described_records`).
pub fn check_pairwise_records(
    subjects: &[Option<std::sync::Arc<ParsedSubject>>],
    protein: bool,
    program: &str,
) -> anyhow::Result<()> {
    subjects
        .iter()
        .flatten()
        .filter(|subject| !subject.shown_as_record())
        .try_for_each(|subject| check_pairwise_record(subject, protein, program))
}

/// LOSAT's refusal of a `-subject` record in outfmt 0, as of a database sequence
/// (`algorithm::common::database::check_shown_parsed_subject`): a protein title that NCBI's
/// `x_CleanAndCompress` reads past the end of.
fn check_pairwise_record(
    subject: &ParsedSubject,
    protein: bool,
    program: &str,
) -> anyhow::Result<()> {
    let (_, title) = subject.description(protein);
    if protein && title.reads_past_end {
        return Err(subject.record_refusal(
            "has a title of punctuation that NCBI BLAST+'s x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 crashes in outfmt 0)",
            program,
        ));
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:336-359
/// `s_WrapOutputLine`: `title`, with a newline after the first white space at or after each
/// multiple of 60 bytes.
pub fn wrap_title(title: &[u8], out: &mut Vec<u8>) {
    let mut do_wrap = false;
    for (i, &c) in title.iter().enumerate() {
        if i > 0 && i % 60 == 0 {
            do_wrap = true;
        }
        out.push(c);
        if do_wrap && (c.is_ascii_whitespace() || c == 11) {
            out.push(b'\n');
            do_wrap = false;
        }
    }
}

/// Whether NCBI's XML and outfmt 0 reports join the hits of `a` and `b` when they follow each
/// other in a hit list: their Seq-align ids match (`CSeq_id::Match`, which holds
/// when `Compare` gives `e_YES`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:478-487
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1750-1757
pub fn joined_in_reports(a: &ParsedSubject, b: &ParsedSubject) -> bool {
    compare(&a.alignment_id, &b.alignment_id) == Sic::Yes
}

// No NCBI counterpart: LOSAT joins the hits that NCBI joins only where they are copies of one
// record; it refuses the others (`common::database::refuse_joined_hits`).
/// Whether `later`, which follows `earlier` in a hit list, is shown in `earlier`'s hit by the
/// reports that join following hits (`joined_in_reports`) as LOSAT writes them: the two are
/// copies of one record (the same definition lines but for their taxids). NCBI shows that hit
/// with the names of the copy its lookup gives first, and each alignment with the letters of the
/// copy its scope holds for the Seq-align id; the reports that read the Bioseq refuse either
/// subject where another sequence holding its Seq-align id is not shown as it (`fetched`: other
/// letters, other taxids where they are written), so the copies' names and letters are the same.
pub fn joins_as_copy(earlier: &ParsedSubject, later: &ParsedSubject) -> bool {
    let without_taxid = |defline: &BlastDefline| BlastDefline {
        taxid: None,
        ..defline.clone()
    };
    joined_in_reports(earlier, later)
        && earlier.set.deflines.len() == later.set.deflines.len()
        && earlier
            .set
            .deflines
            .iter()
            .zip(&later.set.deflines)
            .all(|(a, b)| without_taxid(a) == without_taxid(b))
}

/// `joins_as_copy` of the subjects of ranks `earlier` and `later` of `parsed` (the subjects with
/// parsed Seq-ids by rank): false where either has none.
pub fn rank_joins_as_copy(
    parsed: &[Option<std::sync::Arc<ParsedSubject>>],
    earlier: usize,
    later: usize,
) -> bool {
    match (
        parsed.get(earlier).and_then(Option::as_deref),
        parsed.get(later).and_then(Option::as_deref),
    ) {
        (Some(earlier), Some(later)) => joins_as_copy(earlier, later),
        _ => false,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:476-506
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:357-360, 412-417
// (cited at `ParsedSubject::lookup_throws`: the `Hit` loads its subject by the Seq-align id of
// its first Seq-align, and a lookup that throws there throws again outside the hit's `catch`)
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:750-760
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:376-407 (cited at
// `report::blast_xml::writer::Xml5Hits`: outfmt 5's `catch` of the `Hit`'s lookup), 143-164 (cited
// at `report::blast_xml::writer::Xml5Hsps`: the `catch` of each HSP's lookup; neither throws again)
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:833-870 (cited at
// `report::tabular::TabularSpec::row_throw`: a table's row whose subject block catches the throw)
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:891-898 (cited at
// `report::tabular::TabularSpec::reads_subject`: `slen` read outside any `catch`)
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1586-1592 (cited at
// `ParsedSubject::thrown_description`: the description line's `catch`)
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4188-4197
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1881 (`x_InitAlignParams`,
// called by `DisplaySeqalign` before the alignments of outfmt 0 to 4, without a `catch`)
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1950-1981 (cited at
// `report::pairwise::write_subject_alignments_part`: outfmt 0's `catch` of each alignment, which
// writes the skipped line), 1998-2003 (cited at `report::x_pairwise_par::write_serial`: the merge
// of outfmt 1 to 4 loads each shown alignment's subject without a `catch`)
/// What a report does after a subject whose lookup throws where its batch's prefetch is off
/// (`ParsedSubject::lookup_throws`; DECISIONS 62 of line DBID): the one answer of each writer,
/// which [`HitLoads::of`] takes for the run ([`FormatThrow::of_run`]). The first four are ordered
/// by what the reports load after the throw (each loads at least what the one before it loads);
/// the pairwise reports' two (`is_pairwise`) load the subjects of their description lines and of
/// their alignments, an order of their own (sorted by `-sorthits`), which is in no such order with
/// the others ([`FormatThrow::and`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum AfterThrow {
    /// NCBI's application ends the run with the exception (`CATCH_ALL`'s `BLAST Database error:
    /// <message>`, exit status 2; `blast_xml::writer::run_ending_exception`): outfmt 13 and 14,
    /// and a table whose rows read the subject only after its subject block (`slen`, the letters
    /// and frames), at the subject's first row that reads it.
    EndsRun,
    /// The report function catches it and leaves the query's report out (outfmt 15 and 16,
    /// `BlastXML2_FormatReport`, `BlastJSON_FormatReport`); the next query's report goes on.
    EndsReport,
    /// BLAST XML, outfmt 5: the `catch` of the `Hit` writes it with the Seq-align id and
    /// `Unknown`, and the `catch` of each of its HSPs writes the HSP's Seq-align scores alone
    /// (`blast_xml::writer`, `XmlHit::lookup_throws`); the report goes on with the next `Hit`,
    /// which looks its own subject up. The copies of the record in the thrown `Hit` load nothing,
    /// as in outfmt 13 to 16, where a table's row of such a copy looks it up itself (`WritesRow`).
    WritesHit,
    /// The `catch` of the table's subject block ends the row's fields, the row is written with
    /// what the query's formatter kept (`report::tabular::TabularWriter::row`), and the report
    /// goes on: a table with a field of that block (outfmt 6, 7, 10 and 20); the organism report
    /// (outfmt 18, `CTaxFormat::x_InitTaxInfoMap`'s `catch`) leaves the HSP out and goes on.
    WritesRow,
    /// The pairwise report, outfmt 0: the `catch` of the description line writes it with the
    /// Seq-align id and `Unknown` (`ParsedSubject::thrown_description`) and the report goes on,
    /// but with `-html` (`html`), where the run ends there ([`AfterThrow::description_end`]);
    /// then `GetDbType` loads the subject of the first alignment shown, without a `catch` for an
    /// id that is not local: where it throws, the run ends after the query's description lines
    /// (`alignment_end`); otherwise the `catch` of each alignment writes NCBI's skipped line in
    /// place of each alignment of a subject that throws (`report::pairwise::DatabaseHeading::
    /// Thrown`), and the report goes on.
    SkipsAlignment { html: bool },
    /// The query-anchored reports, outfmt 1 to 4: the description lines as outfmt 0's
    /// (`SkipsAlignment`, with `html` alike: these formats write `-html`'s links since line RG
    /// item T4b, `blastinput::app::formats_link_ids`); then `GetDbType` and the merge load the
    /// subject of each alignment shown, in order, without a `catch`: the run ends after the
    /// query's description lines at the first that throws (`alignment_end`).
    EndsAtAlignment { html: bool },
}

impl AfterThrow {
    /// Whether the writer is the pairwise report's (outfmt 0 to 4: `SkipsAlignment`,
    /// `EndsAtAlignment`), whose reports read what their writers show ([`PairwiseReads`]).
    pub fn is_pairwise(self) -> bool {
        matches!(
            self,
            AfterThrow::SkipsAlignment { .. } | AfterThrow::EndsAtAlignment { .. }
        )
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1588-1611 (cited at
    // `report::pairwise::HtmlRows::end`)
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3904-3907
    // With `-html`, the `catch` of a description line's load makes the line's link, and
    // `GetIDUrl` loads the subject again by the same id outside any `catch` of the report: the
    // throw leaves the formatter and NCBI's application ends with it (`CATCH_ALL`, exit status 2,
    // `BLAST Database error: <message>`), with no alignment, and no report or later batch's
    // prefetch after it (line RG item T4a2; review 3 of line DBID pfoff-4, M2: 81 `-html` runs of
    // NCBI BLAST+ 2.17.0 end there; line RG item T4b review b F2: in formats 0 to 4, a later
    // batch's prefetch of an empty taxdb.btd is not reached).
    // The `catch` of an alignment makes no link (showalign.cpp:1950-1981), so the alignments do as
    // without `-html`, and so does a report without description lines (`-num_descriptions 0`).
    /// Where the description lines of a pairwise report end the run: the position, among the hits
    /// whose lines it shows, in their order (`described`), of the first whose load throws
    /// (`throws`, asked in that order up to it), with `-html` (`html` of `SkipsAlignment` and
    /// `EndsAtAlignment`); `None` without it, for the other writers, and where none throws. The one
    /// decision of the model of the loads ([`HitLoads::of`], with the lookups and holds) and of the
    /// writer (`report::pairwise::HtmlRows`, with the loads).
    pub fn description_end<K>(
        self,
        described: impl IntoIterator<Item = K>,
        throws: impl FnMut(K) -> bool,
    ) -> Option<usize> {
        match self {
            AfterThrow::SkipsAlignment { html: true }
            | AfterThrow::EndsAtAlignment { html: true } => described.into_iter().position(throws),
            _ => None,
        }
    }

    /// Where the alignments of a pairwise report (outfmt 0 to 4) end the run: the position, among
    /// the subjects of the alignments it shows, in order (`aligned`), of the one whose load ends it
    /// (`throws`: its load throws; `local`: its Seq-align id is a local id, which `GetDbType` does
    /// not look up). `SkipsAlignment`: the first, where it throws and is not local; `EndsAtAlignment`:
    /// the first that throws; `None` for the other writers and where none ends it. The one decision
    /// of both the model of the loads ([`HitLoads::of`], with the lookups) and the writer
    /// (`report::x_pairwise_par`, with the loads).
    pub fn alignment_end<K: Copy>(
        self,
        aligned: impl IntoIterator<Item = K>,
        throws: impl Fn(K) -> bool,
        local: impl Fn(K) -> bool,
    ) -> Option<usize> {
        let mut aligned = aligned.into_iter();
        match self {
            AfterThrow::SkipsAlignment { .. } => aligned
                .next()
                .filter(|&first| throws(first) && !local(first))
                .map(|_| 0),
            AfterThrow::EndsAtAlignment { .. } => aligned.position(throws),
            AfterThrow::EndsRun
            | AfterThrow::EndsReport
            | AfterThrow::WritesHit
            | AfterThrow::WritesRow => None,
        }
    }
}

/// The counts of the description lines and the alignments that a pairwise report (outfmt 0 to 4)
/// shows of each query (`-num_descriptions`, `-num_alignments`;
/// `report::x_pairwise_par::described_subjects`, `QueryGroup::aligned`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Shown {
    pub descriptions: usize,
    pub alignments: usize,
}

/// What one query's pairwise report (outfmt 0 to 4) reads where its batch's prefetch is off, as
/// its writer reads it (`report::x_pairwise_par::QueryGroup::read_order`, for [`HitLoads::of`] and
/// the texts of its reads): its hits in its order, each the subjects that it shows (its own
/// first, then those of the copies of its record shown in it, each once), the hits whose
/// description lines it shows (indices of `hits`, in the order of the lines:
/// `x_pairwise_par::described_subjects`), and the number of its first hits whose alignments it
/// shows (`QueryGroup::aligned`). Where its lines or alignments end the run is the run's answer
/// ([`AfterThrow::description_end`] with `-html`, [`AfterThrow::alignment_end`]).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PairwiseReads {
    pub hits: Vec<Vec<usize>>,
    pub described: Vec<usize>,
    pub aligned: usize,
}

impl PairwiseReads {
    /// The subjects in the order that the report reads them, each once: those of its description
    /// lines, then those of its alignments (the texts of the reads one by one,
    /// `common::database::subject_read_texts`). A subject of a hit that neither shows is never read.
    pub fn subjects(&self) -> Vec<usize> {
        let mut read = vec![false; self.hits.len()];
        let mut subjects = Vec::new();
        for hit in self.described.iter().copied().chain(0..self.aligned) {
            if !std::mem::replace(&mut read[hit], true) {
                subjects.extend(&self.hits[hit]);
            }
        }
        subjects
    }
}

/// The subjects of one query's report for [`HitLoads::of`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ReportKeys {
    /// The subject key of each of the report's entries (HSP lists or HSPs) in its order (`None`
    /// for an entry that shows nothing), of which the entries make `Hit`s as the XML writers make
    /// them (`blast_xml::hit_ranges`); of a pairwise report (outfmt 0 to 4), those of a batch that
    /// prefetches, which loads every subject, in any order.
    Entries(Vec<Option<usize>>),
    /// The reads of a pairwise report whose batch reads one by one ([`PairwiseReads`]).
    Pairwise(PairwiseReads),
}

impl ReportKeys {
    /// The subject keys of the report, each entry's (or each hit's subjects) in order.
    fn keys(&self) -> Vec<usize> {
        match self {
            ReportKeys::Entries(entries) => entries.iter().flatten().copied().collect(),
            ReportKeys::Pairwise(reads) => reads.hits.iter().flatten().copied().collect(),
        }
    }
}

impl From<Vec<Option<usize>>> for ReportKeys {
    fn from(entries: Vec<Option<usize>>) -> Self {
        ReportKeys::Entries(entries)
    }
}

/// What the reports of a run do after a subject whose lookup throws where its batch's prefetch is
/// off ([`FormatThrow::of_run`]): `after` ([`AfterThrow`]), at the rows that load the subject,
/// `rows`: every row, or for a table that reads the subject for `nident` alone the rows of HSPs
/// with identities (`report::tabular::SubjectFetch::WithIdentities`, `TabularSpec::throw_rows`);
/// a row of no identities there reads nothing of the subject, so its load neither throws nor ends
/// anything (it prints the 0 of `x_ResetFields`). `first_row_ends`: a table whose run NCBI ends
/// where a report's first row throws (`report::tabular::first_thrown_row_ends_run`: the null
/// pointer of `sacc` or `saccver`, exit 255), which LOSAT refuses
/// (`report::tabular::check_failed_lookups`); no report after it is reached, so no subject is
/// held in one (review 1 of line DBID pfoff-4, H2: the null pointer's words, not a later report's
/// holder's). `keeps_scope`: the run writes outfmt 5, whose formatter keeps NCBI's scope for the
/// whole run (`CBlastFormat::ResetScopeHistory` empties it after each query batch for every other
/// format), so that a sequence that NCBI's scope finds for a subject whose lookup fails holds it
/// in a later batch too ([`ScopeHolders::scope`]; review 2 of line DBID pfoff-4, M1). The formats
/// of one search of a call of the library, outfmt 5 and another, take the scope of the outfmt 5
/// run for all ([`FormatThrow::and`]): a subject that that run holds is refused whatever the other
/// format (LOSAT does not write the holder's record for it), and one that it does not hold is held
/// in no run of those formats, whose scope holds no more, so each run throws where the search
/// throws. The call's searches by hit list size each walk with their own formats' answer
/// ([`FormatThrow::of_search`] keeps it where the call does not refuse), each format's own NCBI
/// run (review 3 of line DBID pfoff-4, L1).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RunThrow {
    pub after: AfterThrow,
    pub rows: super::tabular::SubjectFetch,
    pub first_row_ends: bool,
    pub keeps_scope: bool,
}

impl RunThrow {
    /// What a format whose every row (or `Hit`) loads its subject does: `after`. BLAST XML's
    /// (`AfterThrow::WritesHit`, outfmt 5) keeps NCBI's scope for the run.
    pub fn every_row(after: AfterThrow) -> Self {
        RunThrow {
            after,
            rows: super::tabular::SubjectFetch::Every,
            first_row_ends: false,
            keeps_scope: after == AfterThrow::WritesHit,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:750-760 (cited at
    // `report::blast_xml2::framing::Xml2Kind::after_throw`: the `catch` of outfmt 15's and 16's
    // report posts the exception)
    /// Whether NCBI's runs of a format that does `self` and of one that does `other` end apart
    /// after a throw ([`FormatThrow::and`]): with another exit status, or other standard errors,
    /// which a call of the library cannot both write (it has one result and one stream of
    /// diagnostics, `api::run`). A pairwise report's run ends at an alignment
    /// (`AfterThrow::alignment_end`), an `EndsRun` format's at the first row that loads the
    /// subject (`rows`); outfmt 15 and 16 (`EndsReport`) go on but post NCBI's error on standard
    /// error for each report that they leave out (`Error: [blastp] DB contains no accession
    /// info.`); outfmt 5 and the formats that write the row (`WritesHit`, `WritesRow`) go on and
    /// write nothing there, as a format that reads nothing of the subject does
    /// ([`FormatThrow::Unread`], `None` here: review 2 of line DBID pfoff-4, M3). Two formats of
    /// one writer agree.
    fn stops_apart(self, other: Option<Self>) -> bool {
        let quiet = |run: Option<Self>| {
            run.is_none_or(|run| matches!(run.after, AfterThrow::WritesHit | AfterThrow::WritesRow))
        };
        let Some(other) = other else {
            return !quiet(Some(self));
        };
        if quiet(Some(self)) && quiet(Some(other)) {
            false
        } else {
            self.after != other.after
                || (self.after == AfterThrow::EndsRun && self.rows != other.rows)
        }
    }
}

/// Why LOSAT refuses a subject whose lookup throws where its batch's prefetch is off
/// ([`FormatThrow::Refused`]): the one owner of the words of each situation
/// ([`ThrowRefusal::words`], in the refusal's reason, `LookupFailure::what_ncbi_does`). Ordered:
/// a format that refuses it alone comes before a mix of formats (`FormatThrow::and`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum ThrowRefusal {
    /// A call of the library whose formats NCBI's runs of each end apart after the throw
    /// (`FormatThrow::and`, `RunThrow::stops_apart`): outfmt 0, outfmt 1 to 4, the formats that
    /// end the run there (outfmt 13 and 14, and the tables that read the subject after its subject
    /// block, at the rows that read it), outfmt 15 and 16 (which post the error on standard error)
    /// and the other formats (those that read nothing of the subject among them), in one search
    /// or in one search per hit list size (`FormatThrow::of_search`).
    Mixed,
    /// A format whose writer LOSAT does not make write it as NCBI does (outfmt 17).
    Format,
}

impl ThrowRefusal {
    /// What LOSAT does not write, after NCBI's part of the reason
    /// (`LookupFailure::NCBI_THROWS`): the command line meets `Format` alone (its one format),
    /// a call of the library either.
    pub fn words(self) -> &'static str {
        match self {
            ThrowRefusal::Format => "LOSAT does not write that for outfmt 17",
            ThrowRefusal::Mixed => "LOSAT does not write that in a call of its library whose formats NCBI BLAST+ stops at different places after it or answers with different standard errors (outfmt 0; outfmt 1 to 4; outfmt 13, 14 and the tables that read the subject after its subject block, each at the rows that read it; outfmt 15 and 16, which post the error on standard error; the other formats, those that read nothing of the subject among them), also where the call searches once for each hit list size: NCBI BLAST+ writes each format in a run of its own, and those runs end with different exit statuses or standard errors, where the call has one result and one standard error",
        }
    }
}

/// What one format of a run does with a subject whose lookup throws where its batch's prefetch is
/// off ([`RunThrow`]), as [`FormatThrow::of_run`] takes it for the run.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum FormatThrow {
    /// LOSAT refuses such a subject (`ParsedSubject::fetched`), for the reason given: a format
    /// whose writer LOSAT does not make write it as NCBI does (outfmt 17), or a call whose
    /// formats stop at different places after it (`and`, `of_search`).
    Refused(ThrowRefusal),
    /// The format reads nothing of the subject (a table of fields that read nothing of it; the
    /// Seq-annots and the archive), so NCBI never looks it up there.
    #[default]
    Unread,
    /// The format writes it as NCBI does, then does as the [`RunThrow`] says.
    Writes(RunThrow),
}

impl FormatThrow {
    /// What the format `format` does, with its field list `spec` (`None` for a format without
    /// one, and for the default outfmt 6 and 7, whose `saccver` is a field of the subject block),
    /// in `program`'s reports, with `-html` where `html` (the pairwise reports',
    /// `PairwiseView::after_throw`).
    pub fn of(
        format: crate::blastinput::app::ReportFormat,
        spec: Option<&super::tabular::TabularSpec>,
        program: crate::blastinput::app::Program,
        html: bool,
    ) -> Self {
        use crate::blastinput::app::ReportFormat;
        if let Some(kind) = super::blast_xml2::Xml2Kind::of(format) {
            return FormatThrow::Writes(RunThrow::every_row(kind.after_throw()));
        }
        if let Some(view) = super::pairwise::PairwiseView::of(format) {
            return FormatThrow::Writes(RunThrow::every_row(view.after_throw(html)));
        }
        match format {
            ReportFormat::Tabular
            | ReportFormat::TabularWithComments
            | ReportFormat::CommaSeparatedValues
            | ReportFormat::CommaSeparatedValuesWithHeader
            | ReportFormat::TaxFormat => spec.map_or(
                FormatThrow::Writes(RunThrow {
                    first_row_ends: super::tabular::first_thrown_row_ends_run(None),
                    ..RunThrow::every_row(AfterThrow::WritesRow)
                }),
                |spec| spec.format_throw(program),
            ),
            ReportFormat::Xml => FormatThrow::Writes(RunThrow::every_row(AfterThrow::WritesHit)),
            ReportFormat::AsnText
            | ReportFormat::AsnBinary
            | ReportFormat::ArchiveFormat
            | ReportFormat::JsonSeqalign => FormatThrow::Unread,
            _ => FormatThrow::Refused(ThrowRefusal::Format),
        }
    }

    /// What the `-outfmt` value `outfmt` does in `program`'s reports, as NCBI reads it
    /// (`app::parse_formatting_string`; a value that NCBI refuses is refused) and `app`'s
    /// formatter writes it (`TabularSpec::of_app`), with `-html` where `html`.
    pub fn of_outfmt(
        outfmt: impl AsRef<[u8]>,
        program: crate::blastinput::app::Program,
        app: crate::report::archive::options::RequestKind,
        html: bool,
    ) -> Self {
        match parsed_outfmt(outfmt, app) {
            Some((format, spec)) => Self::of(format, spec.as_ref(), program, html),
            None => FormatThrow::Refused(ThrowRefusal::Format),
        }
    }

    fn report_format(outfmt: impl AsRef<[u8]>) -> Option<crate::blastinput::app::ReportFormat> {
        use crate::blastinput::app::{parse_formatting_string, ReportFormat};
        ReportFormat::from_number(parse_formatting_string(outfmt).ok()?.number)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2295-2314
    /// Whether NCBI's formatter keeps its scope across the query batches of a run whose first
    /// format is `outfmt` (`CBlastFormat::ResetScopeHistory` returns without emptying the scope
    /// only for `eXml`, outfmt 5).
    pub fn outfmt_keeps_scope(outfmt: impl AsRef<[u8]>) -> bool {
        Self::report_format(outfmt) == Some(crate::blastinput::app::ReportFormat::Xml)
    }

    /// Whether `outfmt` is an XML2 or JSON format (outfmt 13 to 16), whose taxonomy lookups run
    /// through `read_reported_taxonomy`.
    pub fn is_xml2(outfmt: impl AsRef<[u8]>) -> bool {
        Self::report_format(outfmt)
            .and_then(super::blast_xml2::Xml2Kind::of)
            .is_some()
    }

    /// `(keeps_scope, separate_xml2)` for the first format `outfmt` of a run (`None` when a call
    /// has no output format): whether NCBI's formatter keeps its scope across query batches
    /// ([`Self::outfmt_keeps_scope`]), and whether any XML2/JSON format in the call is secondary
    /// (`!Self::is_xml2`).
    pub fn scope_flags(outfmt: Option<impl AsRef<[u8]>>) -> (bool, bool) {
        match outfmt {
            Some(outfmt) => (
                Self::outfmt_keeps_scope(outfmt.as_ref()),
                !Self::is_xml2(outfmt.as_ref()),
            ),
            None => (false, false),
        }
    }

    /// Whether the run's [`ScopeRecord`] keeps earlier batches' scopes across query batches:
    /// `Some(true)` where the sole format is outfmt 5, `Some(false)` where there is at most one
    /// known non-XML non-SAM format (so `ResetScopeHistory` empties the scope after each batch),
    /// and `None` where multiple formats or outfmt 17 / an unrecognised format may walk both
    /// scope lifetimes.
    pub fn record_keeps_scope<T: AsRef<[u8]>>(
        outfmts: impl IntoIterator<Item = T>,
    ) -> Option<bool> {
        let mut iter = outfmts.into_iter();
        let Some(first) = iter.next() else {
            return Some(false);
        };
        if iter.next().is_some() {
            return None;
        }
        match Self::report_format(first) {
            Some(crate::blastinput::app::ReportFormat::Xml) => Some(true),
            Some(crate::blastinput::app::ReportFormat::Sam) | None => None,
            Some(_) => Some(false),
        }
    }

    /// The answer of two formats of one run: a refusal by either refuses, a format that reads
    /// nothing leaves the other's answer, and of two that write, the one that loads more after the
    /// throw (a subject is loaded where some report loads it), at the rows that either loads at,
    /// with the scope of either that keeps it (`keeps_scope`: outfmt 5's, which holds what the
    /// other's holds). Two formats whose runs in NCBI BLAST+ end apart after the throw refuse
    /// (`ThrowRefusal::Mixed`, `RunThrow::stops_apart`; the library's calls of several formats,
    /// which have one result and one standard error): a pairwise report (outfmt 0 to 4,
    /// `AfterThrow::is_pairwise`) with a format of another writer (it loads in an order of its
    /// own, and its run ends at an alignment), a format whose run ends at the throw
    /// (`AfterThrow::EndsRun`) with one whose run goes on, two such formats that end at different
    /// rows (`rows`; review 1 of line DBID pfoff-4, M1: outfmt 5 writes its whole document with
    /// exit status 0 where `-outfmt "6 qseqid slen"` stops at the throw with exit status 2), and
    /// outfmt 15 or 16, which post the error on standard error, with one that posts nothing; a
    /// format that reads nothing goes on to its end with exit status 0 and posts nothing, so it
    /// is apart from each of these (review 2 of line DBID pfoff-4, M3: `["0", "8"]`). Formats
    /// whose runs all go on and post nothing (outfmt 5, the tables that write the row, the
    /// organism report, the formats that read nothing) agree, and so do 15 and 16; of two
    /// refusals, a format's own comes first.
    pub fn and(self, other: Self) -> Self {
        use super::tabular::SubjectFetch;
        match (self, other) {
            (FormatThrow::Refused(one), FormatThrow::Refused(other)) => {
                FormatThrow::Refused(one.max(other))
            }
            (FormatThrow::Refused(why), _) | (_, FormatThrow::Refused(why)) => {
                FormatThrow::Refused(why)
            }
            (FormatThrow::Unread, FormatThrow::Unread) => FormatThrow::Unread,
            (FormatThrow::Unread, FormatThrow::Writes(run))
            | (FormatThrow::Writes(run), FormatThrow::Unread) => {
                if run.stops_apart(None) {
                    FormatThrow::Refused(ThrowRefusal::Mixed)
                } else {
                    FormatThrow::Writes(run)
                }
            }
            (FormatThrow::Writes(one), FormatThrow::Writes(other))
                if one.stops_apart(Some(other)) =>
            {
                FormatThrow::Refused(ThrowRefusal::Mixed)
            }
            (FormatThrow::Writes(one), FormatThrow::Writes(other)) => {
                FormatThrow::Writes(RunThrow {
                    after: one.after.max(other.after),
                    rows: if one.rows == other.rows {
                        one.rows
                    } else {
                        SubjectFetch::Every
                    },
                    first_row_ends: one.first_row_ends || other.first_row_ends,
                    keeps_scope: one.keeps_scope || other.keeps_scope,
                })
            }
        }
    }

    /// The answer of the formats `formats` of a run ([`FormatThrow::and`] of them; `Unread`
    /// without one): the one answer of every program, which [`HitLoads::of`] takes. `Unread` is
    /// the answer of a format that reads nothing, not the answer of no format: a run of one
    /// format answers as that format does, so the first format starts the answer.
    pub fn of_formats(formats: impl IntoIterator<Item = Self>) -> Self {
        formats
            .into_iter()
            .reduce(FormatThrow::and)
            .unwrap_or(FormatThrow::Unread)
    }

    /// Whether NCBI keeps one scope for the run ([`ScopeWalk`], `RunThrow::keeps_scope`: outfmt
    /// 5); not where no format reads the subjects (outfmt 5 reads them), and `None` where a format
    /// refuses a throw (`Refused`: the formats' answer no longer tells).
    pub fn keeps_scope(self) -> Option<bool> {
        match self {
            FormatThrow::Writes(run) => Some(run.keeps_scope),
            FormatThrow::Unread => Some(false),
            FormatThrow::Refused(_) => None,
        }
    }

    /// What the reports do after a subject whose lookup throws: `None` where a format refuses it
    /// or no format reads it, so that every subject is loaded ([`HitLoads::of`]).
    pub fn run(self) -> Option<RunThrow> {
        match self {
            FormatThrow::Writes(run) => Some(run),
            FormatThrow::Refused(_) | FormatThrow::Unread => None,
        }
    }

    /// [`FormatThrow::run`] of the answer of the formats `formats` ([`FormatThrow::of_formats`]).
    pub fn of_run(formats: impl IntoIterator<Item = Self>) -> Option<RunThrow> {
        Self::of_formats(formats).run()
    }

    /// [`FormatThrow::of_formats`] of the `-outfmt` values `outfmts` of `program`'s run in `app`'s
    /// application, with `-html` where `html`.
    pub fn of_outfmts<O: AsRef<[u8]>>(
        outfmts: impl IntoIterator<Item = O>,
        program: crate::blastinput::app::Program,
        app: crate::report::archive::options::RequestKind,
        html: bool,
    ) -> Self {
        Self::of_formats(
            outfmts
                .into_iter()
                .map(|outfmt| Self::of_outfmt(outfmt, program, app, html)),
        )
    }

    /// The answer of one search of a call, whose own formats answer `self`: the call's refusal
    /// where the formats of the whole call refuse the throw together, `self` with those that the
    /// call's other searches write (`others`, their `-outfmt` values of `program` in `app`'s
    /// application, with the call's `-html` where `html`: `api::run`'s searches by hit list size,
    /// `ReportOutputs::other_outfmts`; none on the command line), else `self`. A call has one result: the rule is the call's, so that
    /// the split of its formats by `-num_descriptions` and `-num_alignments` does not change it
    /// (review 3 of line DBID pfoff-3, N2). A search whose own formats read nothing keeps
    /// `Unread` where the call does not refuse.
    pub fn of_search(
        self,
        others: &[&str],
        program: crate::blastinput::app::Program,
        app: crate::report::archive::options::RequestKind,
        html: bool,
    ) -> Self {
        let call = Self::of_formats(
            std::iter::once(self).chain(
                others
                    .iter()
                    .map(|outfmt| Self::of_outfmt(outfmt, program, app, html)),
            ),
        );
        match call {
            FormatThrow::Refused(why) => FormatThrow::Refused(why),
            FormatThrow::Unread | FormatThrow::Writes(_) => self,
        }
    }
}

/// A run's answer given as [`FormatThrow::run`] gives it (`None`: every subject is loaded).
impl From<Option<RunThrow>> for FormatThrow {
    fn from(run: Option<RunThrow>) -> Self {
        run.map_or(FormatThrow::Unread, FormatThrow::Writes)
    }
}

/// The report and field list of the `-outfmt` value `outfmt` as NCBI reads it
/// (`app::parse_formatting_string`) and `app`'s formatter writes it (`TabularSpec::of_app`);
/// `None` for a value that NCBI refuses.
fn parsed_outfmt(
    outfmt: impl AsRef<[u8]>,
    app: crate::report::archive::options::RequestKind,
) -> Option<(
    crate::blastinput::app::ReportFormat,
    Option<super::tabular::TabularSpec>,
)> {
    use crate::blastinput::app::{parse_formatting_string, ReportFormat};
    let choice = parse_formatting_string(outfmt).ok()?;
    let format = ReportFormat::from_number(choice.number)?;
    Some((
        format,
        super::tabular::TabularSpec::of_app(format, &choice, app),
    ))
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1146-1155 (cited at
// `common::database::taxdb_error`)
/// Whether `program`'s report `format` (with its field list `spec`; `None` for a format without
/// one and for the default outfmt 6 and 7, whose `saccver` loads the subject) writes what it
/// writes with a readable taxonomy database where a file of that database does not open
/// (`SeqDbError::TaxdbThrows`), so that every load of a database sequence with a positive taxid
/// and every taxid lookup throw: a table or organism report that LOSAT writes the same
/// (`tabular::TabularSpec::unchanged_by_load_throws`),
/// and the Seq-align and archive formats, which load no subject. Such a run's only change is the
/// batches' prefetch, which catches the throw with its warning
/// (`common::database::batch_prefetch`), and psiblast's PSSM, which reads no letters of such a
/// subject (`common::database::TaxdbRun::load_throws`). The other formats read such a subject as
/// their own `catch` decides (`Unknown` ids, `N/A` names, alignments without letters) or end the
/// run, which LOSAT refuses.
pub fn unchanged_by_load_throws(
    format: crate::blastinput::app::ReportFormat,
    spec: Option<&super::tabular::TabularSpec>,
    program: crate::blastinput::app::Program,
) -> bool {
    use crate::blastinput::app::ReportFormat;
    match format {
        ReportFormat::Tabular
        | ReportFormat::TabularWithComments
        | ReportFormat::CommaSeparatedValues
        | ReportFormat::CommaSeparatedValuesWithHeader
        | ReportFormat::TaxFormat => {
            spec.is_some_and(|spec| spec.unchanged_by_load_throws(program))
        }
        ReportFormat::AsnText
        | ReportFormat::AsnBinary
        | ReportFormat::ArchiveFormat
        | ReportFormat::JsonSeqalign => true,
        ReportFormat::Pairwise
        | ReportFormat::QueryAnchoredIdentities
        | ReportFormat::QueryAnchoredNoIdentities
        | ReportFormat::FlatQueryAnchoredIdentities
        | ReportFormat::FlatQueryAnchoredNoIdentities
        | ReportFormat::Xml
        | ReportFormat::Json
        | ReportFormat::Xml2
        | ReportFormat::JsonS
        | ReportFormat::Xml2S
        | ReportFormat::Sam => false,
    }
}

/// [`unchanged_by_load_throws`] of every `-outfmt` value of `outfmts` of `program`'s run in
/// `app`'s application (a value that NCBI refuses is not).
pub fn unchanged_by_load_throws_outfmts<O: AsRef<[u8]>>(
    outfmts: impl IntoIterator<Item = O>,
    program: crate::blastinput::app::Program,
    app: crate::report::archive::options::RequestKind,
) -> bool {
    outfmts.into_iter().all(|outfmt| {
        parsed_outfmt(outfmt, app)
            .is_some_and(|(format, spec)| unchanged_by_load_throws(format, spec.as_ref(), program))
    })
}

/// What NCBI's reports do at the `Hit` (outfmt 5, 13 to 16), the rows (outfmt 6, 7, 10 and 20) or
/// the organism report's HSPs (outfmt 18) of a subject: the one answer to "this subject's lookup
/// throws there" and "this subject is never loaded", which the refusals
/// (`ParsedSubject::fetched`, the tables' rows at which NCBI has no result,
/// `report::tabular::check_failed_lookups`), the names of the `Hit`s
/// (`blast_xml::subjects::HitNames`), the tables' rows (`report::tabular::SubjectLookup::Throws`)
/// and their names and taxonomy (`common::tabular::tabular_taxonomy`, the taxonomy database's
/// lookups), the organism report (`tax_report::SubjectTable`), the taxonomy reads
/// (`common::database::read_reported_taxonomy`, the texts of the subjects that a report reads one
/// by one) and the run's record of the loads (`TaxdbRun::note_loaded`, the stop's `before`) take
/// ([`HitLoads::of`]).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum HitLoad {
    /// A report loads the subject (or the batch's prefetch does), as every subject of a run that
    /// writes another format.
    #[default]
    Loaded,
    /// The subject's lookup throws where a report loads it (`lookup_throws`): NCBI loads nothing
    /// for it, its rows are written as NCBI's formatter leaves them, and the report goes on or
    /// ends at its first row that loads it as the run's [`RunThrow`] says.
    Throws,
    /// No report loads the subject: it is shown only after a subject that throws where the
    /// reports end there, or only in rows that read nothing of it (a table of `nident` alone, HSPs
    /// of no identities).
    Unread,
    /// A report loads the subject, as `Loaded`, in a run whose formats LOSAT refuses after a
    /// subject whose lookup throws (`FormatThrow::Refused`): such a lookup is refused there
    /// (`ParsedSubject::fetched`), for this reason.
    Refused(ThrowRefusal),
    /// A load of the subject gives another sequence ([`ScopeWalk`]: NCBI's scope or data loader
    /// holds it from a load before, or the prefetch's order gives it first; `ShownInstead::holder`, `None`
    /// where LOSAT cannot name it): nothing throws there, and NCBI writes that sequence's record
    /// for the subject, which LOSAT refuses (`ParsedSubject::fetched`). Where the subject throws in
    /// one report and is shown in a later one, it is shown.
    Shown(ShownInstead),
}

/// How one query's report reads its subjects in [`HitLoads::of`]: whether its batch's prefetch is
/// off or fails, so that the report loads each subject where it reaches it (`one_by_one`,
/// `common::database::BatchPrefetch::loads_one_by_one`), the query batch whose scope it reads in
/// (`batch`: a number that the program gives each query batch of the run, the same in its record,
/// `common::database::TaxdbRun::loads_before`; `None` for a report read in the scope of the report
/// before, as `From<bool>` makes it: one query batch), and whether the batch's prefetch ran and
/// failed at a subject whose lookup throws on the id (`prefetch_fails`,
/// `ParsedSubject::prefetch_fails`; `one_by_one` then holds too): the prefetch may have loaded
/// each other subject of the batch before the failure, in an order that LOSAT does not follow
/// (`common::database::batch_prefetch` refuses such a batch where a load writes a text), so the
/// walk takes every subject of the batch's reports whose load does not throw as loaded before
/// them, a superset of the loads that NCBI's scope may hold there (line DBID pfoff-7).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReportRead {
    pub one_by_one: bool,
    pub batch: Option<usize>,
    pub prefetch_fails: bool,
}

impl From<bool> for ReportRead {
    fn from(one_by_one: bool) -> Self {
        ReportRead {
            one_by_one,
            batch: None,
            prefetch_fails: false,
        }
    }
}

/// [`HitLoad`] of each subject of a call's reports, by the key that the caller gives its subjects
/// (a rank, or an ordinal id): those whose lookup throws where a report reaches them, with the
/// exception, those for which a load gives another sequence ([`ScopeWalk`]), and those that no
/// report loads; every other subject is loaded. Where a format of the run refuses a lookup that
/// throws, every subject that no other sequence stands for is loaded, with the refusal's reason.
/// With the reports that NCBI's run reaches ([`HitLoads::reaches`]).
#[derive(Debug, Default)]
pub struct HitLoads {
    throws: std::collections::HashMap<usize, crate::ncbi_exception::NcbiException>,
    /// The subjects for which a load gives another sequence (`HitLoad::Shown`).
    shown: std::collections::HashMap<usize, ShownInstead>,
    unread: std::collections::HashSet<usize>,
    /// The subjects that a batch's failed prefetch may have loaded (`ReportRead::prefetch_fails`)
    /// and that no report reaches: loaded, but read by no report (`HitLoads::reached`).
    prefetched_only: std::collections::HashSet<usize>,
    /// Why the run's formats refuse a lookup that throws (`FormatThrow::Refused`).
    refused: Option<ThrowRefusal>,
    /// The number of the reports given that NCBI's run reaches, where a report's load ended it
    /// (`None`: the run reaches every one).
    reached: Option<usize>,
    /// What each batch's prefetch gave for each subject of its prefetch set, by key
    /// ([`HitLoads::prefetched`]).
    prefetched: std::collections::HashMap<usize, HitLoad>,
    /// The subjects that a report of a batch without the prefetch reads ([`HitLoads::read_alone`]).
    read_alone: std::collections::HashSet<usize>,
    /// The reads after the reports ([`CallReads::after`]) for which a load gives another
    /// sequence, by report and key ([`HitLoads::read_after`], made by `HitLoads::load_after`).
    after_shown: std::collections::HashMap<(usize, usize), ShownInstead>,
    /// The sequence and its loaded line that the loads (each report's, or its batch's prefetch)
    /// give for a subject where that is another sequence ([`ShownInstead::written`]), by key:
    /// `None` where two loads give other ones, or one names no line ([`HitLoads::written_holder`]).
    /// `None` also where another load gives the subject itself (one subject can be itself in one
    /// batch's scope and another sequence in the next): the walk notes those loads only while it
    /// runs (`HitLoads::after`'s `own`), so that a call keeps nothing per subject that no load
    /// shows as another sequence.
    written: std::collections::HashMap<usize, Option<(DbSequence, usize)>>,
    /// Whether the run's formats write the record of the sequence that a load gives for a subject
    /// (`ShownInstead::written`, the record writers of line DBID hr-2a: the tables of outfmt 6,
    /// 7, 10 and 20, BLAST XML, outfmt 13 to 16 and 18); not the pairwise reports and SAM yet, nor a
    /// run whose formats refuse a throw (`HitLoads::unordered`, whose holders are only possible).
    writes_shown: bool,
}

/// What a call's walk ([`HitLoads::of_call`]) loads besides its reports' own subjects, by report
/// (the index of the report among those given): `prefetched`, the further subjects of the
/// report's batch's prefetch, which no report shows (psiblast's traceback hits beyond
/// `-max_target_seqs`, `report_prep`'s `psiblast_prefetched`; the batch's prefetch loads those
/// of all its reports at once); `after`, the subjects read by their Seq-align ids after the
/// report, in the order of the reads (psiblast's PSSM of the query, read after its report and
/// before the next one, `ListedTitles::check_pssm_holds`); and `stop`, the prefetch of the batch
/// after the reports whose lookup stops the run (`common::database::batch_prefetch`): its query
/// batch and subjects. Empty for a call that loads nothing else.
#[derive(Clone, Debug, Default)]
pub struct CallReads {
    pub prefetched: Vec<Vec<usize>>,
    pub after: Vec<Vec<usize>>,
    pub stop: Option<(Option<usize>, Vec<usize>)>,
}

impl CallReads {
    /// The further prefetched subjects of report `report`.
    fn prefetched(&self, report: usize) -> &[usize] {
        self.prefetched.get(report).map_or(&[], Vec::as_slice)
    }

    /// The reads after report `report`.
    fn after(&self, report: usize) -> &[usize] {
        self.after.get(report).map_or(&[], Vec::as_slice)
    }
}

impl HitLoads {
    /// The loads of the query reports `reports`, in the order of the reports: whether the
    /// query's batch has its prefetch off, so that its report loads each subject where it reaches
    /// it (`one_by_one`, `common::database::BatchPrefetch::loads_one_by_one`), and the subject key
    /// of each of its HSP lists (or HSPs) in the report's order (`None` for an entry that shows
    /// nothing). Each load asks NCBI's places ([`ScopeWalk`], after the loads `before` of the
    /// run: its record of the batches, rounds and calls before these reports; [`LoadsBefore::none`]
    /// for the run's first reports), which give the subject itself, another sequence
    /// (`HitLoad::Shown`) or a lookup that throws. The entries make `Hit`s as the XML writers make
    /// them (`blast_xml::hit_ranges`: a copy of one record that follows another joins its `Hit`,
    /// `joins_as_copy`; the entries of one subject that follow each other are one). In a report
    /// that reads one by one, a `Hit` loads its first subject at its first row that loads it: every
    /// row for the run's `rows` `Every`, the rows of HSPs with identities for `WithIdentities`
    /// (`identities(report, entry)`: the largest of the identities of the HSPs of entry `entry` of
    /// report `report`, the Seq-align's `num_ident`, asked only there); a `Hit` none of whose rows
    /// loads it is not loaded and ends nothing (`HitLoad::Unread`, unless another report loads it).
    /// A `Hit` whose load throws (`lookup_throws` of `parsed(key)`) loads nothing; what comes after
    /// it is the run's `after` ([`FormatThrow::of_run`]):
    /// - `WritesRow` (a table of a field of the subject block, the organism report): each row looks
    ///   its own Seq-align id up, so every subject whose lookup throws throws, every other one
    ///   loads;
    /// - `WritesHit` (outfmt 5): each `Hit` looks its own subject up, so every `Hit` whose
    ///   subject's lookup throws throws, every other one loads; the copies in a `Hit` that throws
    ///   load nothing;
    /// - `EndsReport` (outfmt 15 and 16): the `Hit`s after it in the report load nothing, the next
    ///   report goes on;
    /// - `EndsRun` (outfmt 13 and 14, a table that reads the subject after its subject block):
    ///   nothing after it in the run loads (`ended`, carried from the reports before: see
    ///   [`HitLoads::after`]);
    /// - `SkipsAlignment` and `EndsAtAlignment` (outfmt 0 to 4): a report that reads one by one
    ///   gives its reads as its writer makes them (`ReportKeys::Pairwise`,
    ///   `x_pairwise_par::QueryGroup::read_order`): each description line loads its subject or
    ///   throws (with `-html` the first that throws ends the run, `AfterThrow::description_end`),
    ///   then the alignments load theirs up to the one that ends the run, if any
    ///   (`pairwise_report`); one whose batch prefetches gives its entries in any order.
    ///
    /// A batch whose prefetch runs loads every subject of its reports at once before them
    /// ([`ScopeWalk::bulk`]; a lookup that fails on a file stops the prefetch first,
    /// `common::database::batch_prefetch`); one whose prefetch fails at a lookup that throws on the
    /// id (`ReportRead::prefetch_fails`) may have loaded each other subject of the batch before
    /// the failure ([`ScopeWalk::failed_prefetch`]), and its reports then read one by one. A
    /// subject is loaded where some report loads it; one that throws in a report throws in every
    /// report that reaches it. `run` (the run's answer, [`FormatThrow::of_formats`]; `None` for
    /// `Unread`): where a format of the run refuses the throw (`HitLoad::Refused`, with the
    /// reason), LOSAT does not know the order of the loads of each format's NCBI run, so a subject
    /// is shown wherever another sequence loaded in the run may stand for it
    /// ([`HitLoads::unordered`]); where no format reads the subjects, only the batches' prefetch
    /// loads them.
    ///
    /// NCBI's application ends at the report whose load ends the run (`CATCH_ALL`): no report,
    /// and no batch's prefetch, after it runs ([`HitLoads::reaches`]).
    pub fn of<'p, R: Into<ReportRead>, K: Into<ReportKeys>>(
        before: LoadsBefore,
        run: impl Into<FormatThrow>,
        reports: impl IntoIterator<Item = (R, K)>,
        parsed: impl Fn(usize) -> Option<&'p ParsedSubject>,
        identities: impl Fn(usize, usize) -> usize,
    ) -> Self {
        Self::of_call(
            before,
            run,
            reports,
            parsed,
            identities,
            &CallReads::default(),
        )
    }

    /// [`HitLoads::of`] in the one walk of a call that also loads the subjects of `call`
    /// ([`CallReads`]): each batch's prefetch loads the further subjects of its reports with
    /// theirs, the reads after a report come after its loads (none after the report whose load
    /// ends the run), and a stopping batch's prefetch after the reports.
    pub fn of_call<'p, R: Into<ReportRead>, K: Into<ReportKeys>>(
        before: LoadsBefore,
        run: impl Into<FormatThrow>,
        reports: impl IntoIterator<Item = (R, K)>,
        parsed: impl Fn(usize) -> Option<&'p ParsedSubject>,
        identities: impl Fn(usize, usize) -> usize,
        call: &CallReads,
    ) -> Self {
        let run = run.into();
        Self::after(
            &std::cell::Cell::new(false),
            &mut ScopeWalk::new(before, run.keeps_scope()),
            (run, call),
            reports,
            parsed,
            identities,
        )
    }

    /// [`HitLoads::of`] of the reports `reports` that follow those whose loads ended the run where
    /// `ended` holds true (an `EndsRun` subject that threw), which it holds after these reports
    /// too, with the walk `walk` of the loads before ([`RunLoads::loads`]).
    fn after<'p, R: Into<ReportRead>, K: Into<ReportKeys>>(
        ended: &std::cell::Cell<bool>,
        walk: &mut ScopeWalk,
        (run, call): (FormatThrow, &CallReads),
        reports: impl IntoIterator<Item = (R, K)>,
        parsed: impl Fn(usize) -> Option<&'p ParsedSubject>,
        identities: impl Fn(usize, usize) -> usize,
    ) -> Self {
        let mut loads = HitLoads::default();
        let reports: Vec<(ReportRead, ReportKeys)> = reports
            .into_iter()
            .map(|(read, keys)| (read.into(), keys.into()))
            .collect();
        let run = match run {
            FormatThrow::Writes(run) => Some(run),
            FormatThrow::Refused(why) => {
                loads.refused = Some(why);
                loads.unordered(walk, (&reports, call), &parsed);
                return loads;
            }
            FormatThrow::Unread => None,
        };
        loads.writes_shown = run.is_some_and(|run| !run.after.is_pairwise());
        // The answers of the prefetch of the batch of the reports read last (`ScopeWalk::bulk`).
        let mut batch_answers: std::collections::HashMap<usize, HitLoad> = Default::default();
        let mut loaded = std::collections::HashSet::new();
        // The subjects that some load gives as themselves (`HitLoads::written`).
        let mut own = std::collections::HashSet::new();
        let mut passed = std::collections::HashSet::new();
        // The subjects of the reports of each batch whose prefetch fails, which that prefetch may
        // have loaded before its failure (`ReportRead::prefetch_fails`), and those of them noted
        // (`failed_prefetched`): loaded for the holders and the run's record, but read only where
        // a report reaches them (review 1 of line DBID pf7, H1).
        let mut failed_prefetches: std::collections::HashMap<Option<usize>, Vec<usize>> =
            Default::default();
        let mut failed_prefetched = std::collections::HashSet::new();
        for (report, (read, keys)) in reports.iter().enumerate() {
            if read.prefetch_fails {
                let subjects = failed_prefetches.entry(read.batch).or_default();
                subjects.extend(keys.keys());
                subjects.extend_from_slice(call.prefetched(report));
            }
        }
        // The reports given, of which the run reaches those before the first that comes after
        // its end (`reached`).
        let mut given = 0;
        // The batch whose prefetch ran last (`ScopeWalk::bulk`), by the index of its first report.
        let mut prefetched = None;
        for (report, (read, keys)) in reports.iter().enumerate() {
            // The reads after the report before (psiblast's PSSM), in its batch's scope.
            if let Some(before) = report.checked_sub(1) {
                loads.load_after(before, call, (ended, &mut *walk), &parsed);
            }
            if ended.get() {
                loads.reached.get_or_insert(report);
            }
            given = report + 1;
            walk.enter(read.batch);
            // The batch's prefetch, which the run reaches, may have loaded every subject of the
            // batch whose load does not throw before it failed (`ScopeWalk::failed_prefetch`).
            if let Some(keys) = failed_prefetches
                .remove(&read.batch)
                .filter(|_| read.prefetch_fails && !ended.get())
            {
                failed_prefetched.extend(loads.failed_prefetch(walk, keys, &parsed));
            }
            if !read.one_by_one {
                // A batch's prefetch loads every subject of its reports before them, and the
                // further subjects of its prefetch set (`CallReads::prefetched`), unless the run
                // ended before.
                if ended.get() {
                    passed.extend(keys.keys());
                } else if prefetched.is_none_or(|first| !same_batch(&reports, first, report)) {
                    prefetched = Some(report);
                    let (shown_keys, further) = batch_prefetch_set(&reports, report, call);
                    let subjects: Vec<(usize, &ParsedSubject)> = shown_keys
                        .iter()
                        .chain(&further)
                        .filter_map(|&key| parsed(key).map(|subject| (key, subject)))
                        .collect();
                    let answers = walk.bulk(subjects.iter().copied());
                    for (&key, &answer) in &answers {
                        loads.prefetched.insert(key, answer);
                        if let HitLoad::Shown(shown) = answer {
                            if shown_keys.contains(&key) {
                                loads.shown.insert(key, shown);
                            }
                        }
                    }
                    batch_answers.clone_from(&answers);
                    // The description lines of a pairwise report look each title up again
                    // (`ScopeWalk::display`); the bulk resolved every Seq-align id, so the order
                    // of the lines does not change what each finds.
                    if run.is_some_and(|run| run.after.is_pairwise()) {
                        for &(key, subject) in &subjects {
                            if answers.get(&key) == Some(&HitLoad::Loaded) {
                                if let Some(shown) = walk.display(subject) {
                                    loads.shown.insert(key, shown);
                                }
                            }
                        }
                    }
                }
                if !ended.get() {
                    loaded.extend(keys.keys());
                    loads.note_batch_answers(keys, &batch_answers, &mut own);
                }
                continue;
            }
            // A report of a batch without the prefetch reads its subjects by its own lookups.
            loads.read_alone.extend(keys.keys());
            let Some(run) = run else {
                // No format reads the subjects: a report that reads one by one loads nothing.
                continue;
            };
            let subjects = match keys {
                ReportKeys::Pairwise(reads) => {
                    debug_assert!(run.after.is_pairwise());
                    loads.pairwise_report(
                        run,
                        reads,
                        &parsed,
                        (ended, &mut *walk),
                        (&mut loaded, &mut passed),
                    );
                    continue;
                }
                ReportKeys::Entries(subjects) => subjects,
            };
            // A table and the organism report look each row's Seq-align id up.
            let hits = if run.after == AfterThrow::WritesRow {
                (0..subjects.len()).map(|entry| entry..entry + 1).collect()
            } else {
                super::blast_xml::hit_ranges(
                    subjects,
                    |subject| *subject,
                    |earlier, later| match (parsed(earlier), parsed(later)) {
                        (Some(earlier), Some(later)) => joins_as_copy(earlier, later),
                        _ => false,
                    },
                )
            };
            // The report ended at a `Hit` that threw (`EndsReport`), or the run did (`EndsRun`).
            let mut thrown = ended.get();
            for range in hits {
                let Some(first) = subjects[range.start] else {
                    continue;
                };
                let keys = subjects[range.clone()].iter().flatten().copied();
                // The `Hit`'s first row that loads its subject (`rows`).
                let loads_it = run.rows == super::tabular::SubjectFetch::Every
                    || range.clone().any(|entry| {
                        subjects[entry].is_some() && run.rows.reads(identities(report, entry))
                    });
                if thrown {
                    passed.extend(keys);
                    continue;
                }
                let Some(subject) = parsed(first) else {
                    loaded.extend(keys);
                    continue;
                };
                if !loads_it {
                    // No row reads the subject: no load, so nothing is held or thrown; a subject
                    // whose lookup throws is passed, another counts as loaded (its rows read
                    // nothing of it).
                    if subject.lookup_throws().is_some() {
                        passed.extend(keys);
                    } else {
                        loaded.extend(keys);
                    }
                    continue;
                }
                let exception = match walk.load(subject) {
                    HitLoad::Shown(shown) => {
                        loads.shown.insert(first, shown);
                        loads.note_shown(first, shown);
                        loaded.extend(keys);
                        continue;
                    }
                    HitLoad::Throws => subject.lookup_throws(),
                    _ => None,
                };
                let Some(exception) = exception else {
                    own.insert(first);
                    loaded.extend(keys);
                    continue;
                };
                loads.throws.insert(first, exception);
                // Outfmt 13 to 16 load nothing after it in the report (or the run); BLAST XML's
                // next `Hit` and a table's next row look their own subjects up, but where the
                // report's first row throws, a table of `sacc` or `saccver` ends the run.
                let ends_run =
                    run.after == AfterThrow::EndsRun || (run.first_row_ends && range.start == 0);
                thrown = ends_run || run.after == AfterThrow::EndsReport;
                if ends_run {
                    ended.set(true);
                }
                passed.extend(keys.filter(|&key| key != first));
            }
        }
        if let Some(last) = given.checked_sub(1) {
            loads.load_after(last, call, (ended, &mut *walk), &parsed);
        }
        if ended.get() {
            loads.reached.get_or_insert(given);
        } else if let Some((batch, keys)) = &call.stop {
            // The stopping batch's prefetch, after the reports (`CallReads::stop`).
            walk.enter(*batch);
            let answers = walk.bulk(keys.iter().filter_map(|&key| Some((key, parsed(key)?))));
            loads.prefetched.extend(answers);
        }
        let (prefetched_only, unread) = passed
            .into_iter()
            .filter(|key| {
                !loaded.contains(key)
                    && !loads.throws.contains_key(key)
                    && !loads.shown.contains_key(key)
            })
            .partition(|key| failed_prefetched.contains(key));
        loads.prefetched_only = prefetched_only;
        loads.unread = unread;
        for (key, holder) in &mut loads.written {
            if own.contains(key) {
                *holder = None;
            }
        }
        loads
    }

    /// The answers of a batch's prefetch (`answers`, [`ScopeWalk::bulk`]) that a report of that
    /// batch finds in the scope's cache for its subjects `keys`: each another sequence
    /// (`note_shown`), or anything else (`own`: the subject itself, a throw, no load).
    fn note_batch_answers(
        &mut self,
        keys: &ReportKeys,
        answers: &std::collections::HashMap<usize, HitLoad>,
        own: &mut std::collections::HashSet<usize>,
    ) {
        for key in keys.keys() {
            match answers.get(&key) {
                Some(&HitLoad::Shown(shown)) => self.note_shown(key, shown),
                Some(_) => {
                    own.insert(key);
                }
                None => {}
            }
        }
    }

    /// A load's answer `shown` for the subject `key` ([`HitLoads::written`]).
    fn note_shown(&mut self, key: usize, shown: ShownInstead) {
        let holder = shown.written();
        self.written
            .entry(key)
            .and_modify(|written| {
                if *written != holder {
                    *written = None;
                }
            })
            .or_insert(holder);
    }

    /// The prefetch of a batch that fails at a subject whose lookup throws on the id
    /// (`ReportRead::prefetch_fails`), whose prefetch set is `prefetched` (the subjects of its
    /// reports and its further ones, `CallReads::prefetched`): each subject whose load does not
    /// throw may have been loaded before the failure ([`ScopeWalk::failed_prefetch`]); what the
    /// places give then for each subject's Seq-align id ([`ScopeWalk::peek`]: where an earlier
    /// blob of the pass registered the id of one that throws, the prefetch gave that blob for it)
    /// is the prefetch's answer ([`HitLoads::prefetched`]). The keys that may have been loaded.
    fn failed_prefetch<'p>(
        &mut self,
        walk: &mut ScopeWalk,
        prefetched: Vec<usize>,
        parsed: &impl Fn(usize) -> Option<&'p ParsedSubject>,
    ) -> Vec<usize> {
        let keys: Vec<usize> = prefetched
            .iter()
            .copied()
            .filter(|&key| parsed(key).is_none_or(|subject| subject.lookup_throws().is_none()))
            .collect();
        walk.failed_prefetch(keys.iter().filter_map(|&key| Some((key, parsed(key)?))));
        for key in prefetched {
            if let Some(subject) = parsed(key) {
                self.prefetched.insert(key, walk.peek(subject));
            }
        }
        keys
    }

    /// The reads after report `report` (`CallReads::after`, psiblast's PSSM), in the walk `walk`
    /// after its loads, unless the run ended (`ended`): each a load by the subject's Seq-align id
    /// ([`ScopeWalk::load`]), where its batch's prefetch ran found in the scope's cache.
    fn load_after<'p>(
        &mut self,
        report: usize,
        call: &CallReads,
        (ended, walk): (&std::cell::Cell<bool>, &mut ScopeWalk),
        parsed: &impl Fn(usize) -> Option<&'p ParsedSubject>,
    ) {
        if ended.get() {
            return;
        }
        for &key in call.after(report) {
            if let Some(subject) = parsed(key) {
                if let HitLoad::Shown(shown) = walk.load(subject) {
                    self.after_shown.insert((report, key), shown);
                }
            }
        }
    }

    /// The loads of a run whose formats refuse a lookup that throws (`FormatThrow::Refused`):
    /// every subject of `reports` is loaded, by the batches' prefetch ([`ScopeWalk::bulk`], whose
    /// pass does not depend on the formats: in the walk's order) or by its reports; LOSAT does not
    /// know the order of the loads of each format's run that reads one by one, so where another
    /// sequence loaded before (the record's, `walk`) or in these reports holds an id that NCBI's
    /// places take for such a subject's Seq-align id (on any line by its exact handle, on its
    /// loaded line by an id the scope takes for it or matches weakly), the subject is shown
    /// (`ShownBy::Order`).
    fn unordered<'p>(
        &mut self,
        walk: &mut ScopeWalk,
        (reports, call): (&[(ReportRead, ReportKeys)], &CallReads),
        parsed: &impl Fn(usize) -> Option<&'p ParsedSubject>,
    ) {
        let mut subjects: Vec<(usize, &ParsedSubject)> = Vec::new();
        let mut seen = std::collections::HashSet::new();
        let mut prefetched = None;
        for (report, (read, keys)) in reports.iter().enumerate() {
            walk.enter(read.batch);
            if read.one_by_one {
                // The reads after a report (`CallReads::after`) load one by one too.
                for key in keys
                    .keys()
                    .into_iter()
                    .chain(call.after(report).iter().copied())
                {
                    if seen.insert(key) {
                        subjects.extend(parsed(key).map(|subject| (key, subject)));
                    }
                }
            } else if prefetched.is_none_or(|first| !same_batch(reports, first, report)) {
                prefetched = Some(report);
                let (shown_keys, further) = batch_prefetch_set(reports, report, call);
                let batch: Vec<(usize, &ParsedSubject)> = shown_keys
                    .iter()
                    .chain(&further)
                    .filter_map(|&key| parsed(key).map(|subject| (key, subject)))
                    .collect();
                for (key, answer) in walk.bulk(batch) {
                    self.prefetched.insert(key, answer);
                    if let HitLoad::Shown(shown) = answer {
                        self.shown.insert(key, shown);
                    }
                }
            }
        }
        let mut places = ScopePlaces::default();
        let mut registry: rustc_hash::FxHashMap<HandleKey, Vec<DbSequence>> = Default::default();
        for (_, subject) in &subjects {
            if let Some(sequence) = subject.sequence {
                places.add_ids(
                    sequence,
                    subject.target,
                    &subject.set.deflines[subject.target].ids,
                );
                for handle in subject.line_handles() {
                    registry.entry(handle).or_default().push(sequence);
                }
            }
        }
        for (key, subject) in subjects {
            let keys = subject.id_keys();
            let mut others = places
                .holders(&keys, true)
                .into_iter()
                .chain(places.holders(&keys, false))
                .map(|(at, _)| places.sequence(at))
                .chain(
                    keys.handle
                        .and_then(|handle| registry.get(&handle))
                        .into_iter()
                        .flatten()
                        .copied(),
                )
                .filter(|&other| !subject.mine(other));
            let mut shown = others
                .next()
                .map(|holder| ShownInstead::unlined(Some(holder), ShownBy::Order));
            // A description line's second lookup by the display id (`ScopeWalk::display`).
            if shown.is_none() {
                shown = subject.display_keys().and_then(|keys| {
                    places
                        .holders(&keys, true)
                        .into_iter()
                        .map(|(at, _)| places.sequence(at))
                        .find(|&other| !subject.mine(other))
                        .map(|holder| ShownInstead::unlined(Some(holder), ShownBy::Display))
                });
            }
            if shown.is_none() {
                if let HitLoad::Shown(before) = walk.load(subject) {
                    shown = Some(before);
                }
            }
            if let Some(shown) = shown {
                self.shown.insert(key, shown);
            }
        }
        // A read after a report is shown where its subject is (in the reports or by the
        // prefetch, `HitLoads::prefetched`).
        for report in 0..reports.len() {
            for &key in call.after(report) {
                let shown = self
                    .shown
                    .get(&key)
                    .copied()
                    .or(match self.prefetched.get(&key) {
                        Some(&HitLoad::Shown(shown)) => Some(shown),
                        _ => None,
                    });
                if let Some(shown) = shown {
                    self.after_shown.insert((report, key), shown);
                }
            }
        }
        if let Some((batch, keys)) = &call.stop {
            walk.enter(*batch);
            let answers = walk.bulk(keys.iter().filter_map(|&key| Some((key, parsed(key)?))));
            self.prefetched.extend(answers);
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:686-700
    // The description lines load the subjects that they show, in their order: those of
    // `x_pairwise_par::described_subjects` (psiblast's later rounds, `two_part_lines`).
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1014-1040
    // (`PruneSeqalign` keeps the alignments of the first `m_NumAlignments` hits, cited at
    // `report::x_pairwise_par::QueryGroup::aligned`)
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1588-1611 (cited at
    // `AfterThrow::description_end`: with `-html` the run ends at the description line, before the
    // later lines, the alignments and the later batches with their prefetch)
    /// The loads of one pairwise report (outfmt 0 to 4, `AfterThrow::is_pairwise`, whose run does
    /// `run` after a throw) in a batch whose prefetch is off, whose reads are `reads` (its
    /// writer's own answers, `x_pairwise_par::QueryGroup::read_order`): unless the run ended
    /// before (`ended`), its description lines load the subjects of the hits that they show, up to
    /// the one whose load ends the run with `-html` (`AfterThrow::description_end`, the writer's
    /// answer too; then nothing after it loads), then its alignments those of its first
    /// `reads.aligned` hits up to the one whose load ends the run (`AfterThrow::alignment_end`),
    /// each in that order in the walk `walk` ([`ScopeWalk::load`]); `ended` records either end. A
    /// hit whose load throws (`lookup_throws`) loads nothing there, nor do the copies in it, and is
    /// loaded again at its next read; one for which a load gives another sequence is shown
    /// (`HitLoad::Shown`: an alignment of it ends nothing); every subject that neither loads is
    /// passed.
    fn pairwise_report<'p>(
        &mut self,
        run: RunThrow,
        reads: &PairwiseReads,
        parsed: &impl Fn(usize) -> Option<&'p ParsedSubject>,
        (ended, walk): (&std::cell::Cell<bool>, &mut ScopeWalk),
        (loaded, passed): (
            &mut std::collections::HashSet<usize>,
            &mut std::collections::HashSet<usize>,
        ),
    ) {
        /// What the report's reads of a hit gave so far.
        #[derive(Clone, Copy)]
        enum Read {
            Loads,
            Throws,
            Shown(ShownInstead),
        }
        let first = |hit: usize| reads.hits[hit].first().copied().and_then(parsed);
        let state = std::cell::RefCell::new(vec![None::<Read>; reads.hits.len()]);
        let walk = std::cell::RefCell::new(walk);
        // A read of the hit `hit` where the report reaches it, in the report's order: a hit that
        // threw at a read before is looked up again (a holder may have been loaded since).
        let read = |hit: usize| {
            let mut state = state.borrow_mut();
            if let Some(done @ (Read::Loads | Read::Shown(_))) = state[hit] {
                return done;
            }
            let now = match first(hit) {
                Some(subject) => match walk.borrow_mut().load(subject) {
                    HitLoad::Shown(shown) => Read::Shown(shown),
                    HitLoad::Throws if subject.lookup_throws().is_some() => Read::Throws,
                    _ => Read::Loads,
                },
                None => Read::Loads,
            };
            state[hit] = Some(now);
            now
        };
        // A description line's reads: its load, then where that gave the subject itself, the
        // second lookup by its display id (`ScopeWalk::display`), once.
        let displayed = std::cell::RefCell::new(vec![false; reads.hits.len()]);
        let describe = |hit: usize| {
            let now = read(hit);
            if !matches!(now, Read::Loads)
                || std::mem::replace(&mut displayed.borrow_mut()[hit], true)
            {
                return now;
            }
            match first(hit).and_then(|subject| walk.borrow_mut().display(subject)) {
                Some(shown) => {
                    state.borrow_mut()[hit] = Some(Read::Shown(shown));
                    Read::Shown(shown)
                }
                None => now,
            }
        };
        if !ended.get() {
            // With `-html` the description lines end the run at the first whose load throws
            // (`AfterThrow::description_end`, the writer's answer): nothing after it loads.
            let description_end = run
                .after
                .description_end(reads.described.iter().copied(), |hit| {
                    matches!(describe(hit), Read::Throws)
                });
            if description_end.is_some() {
                ended.set(true);
            } else {
                // The description lines that the end above did not read (all, without `-html`; where `-html`
                // read them all, no line threw, and the loop reads them again, which changes nothing:
                // `read` keeps `Loads` and `Shown`).
                for &hit in &reads.described {
                    describe(hit);
                }
                let end = run.after.alignment_end(
                    0..reads.aligned,
                    |hit| matches!(read(hit), Read::Throws),
                    |hit| first(hit).is_some_and(ParsedSubject::alignment_id_is_local),
                );
                // The alignments that the report reaches, each read once more where the end above
                // did not read it (`SkipsAlignment` reads the first alone there).
                for hit in 0..end.map_or(reads.aligned, |end| end + 1) {
                    read(hit);
                }
                ended.set(end.is_some());
            }
        }
        for (hit, keys) in reads.hits.iter().enumerate() {
            let now = state.borrow()[hit];
            match (now, keys.split_first()) {
                (Some(Read::Throws), Some((&first, copies))) => {
                    let exception = parsed(first)
                        .and_then(ParsedSubject::lookup_throws)
                        .expect("a subject that throws");
                    self.throws.insert(first, exception);
                    passed.extend(copies);
                }
                (Some(Read::Shown(shown)), Some((&first, _))) => {
                    self.shown.insert(first, shown);
                    loaded.extend(keys);
                }
                (Some(Read::Loads), _) => loaded.extend(keys),
                _ => passed.extend(keys),
            }
        }
    }

    /// These loads with each subject's key renamed by `key` (an ordinal id to its rank; a key
    /// that `key` does not name is dropped, which reads as loaded).
    pub fn renamed(self, key: impl Fn(usize) -> Option<usize>) -> Self {
        HitLoads {
            throws: self
                .throws
                .into_iter()
                .filter_map(|(old, exception)| Some((key(old)?, exception)))
                .collect(),
            shown: self
                .shown
                .into_iter()
                .filter_map(|(old, shown)| Some((key(old)?, shown)))
                .collect(),
            unread: self.unread.into_iter().filter_map(&key).collect(),
            prefetched_only: self.prefetched_only.into_iter().filter_map(&key).collect(),
            refused: self.refused,
            reached: self.reached,
            prefetched: self
                .prefetched
                .into_iter()
                .filter_map(|(old, load)| Some((key(old)?, load)))
                .collect(),
            read_alone: self.read_alone.into_iter().filter_map(&key).collect(),
            after_shown: self
                .after_shown
                .into_iter()
                .filter_map(|((report, old), shown)| Some(((report, key(old)?), shown)))
                .collect(),
            written: self
                .written
                .into_iter()
                .filter_map(|(old, holder)| Some((key(old)?, holder)))
                .collect(),
            writes_shown: self.writes_shown,
        }
    }

    /// What the reports do at the `Hit` of the subject `key`.
    pub fn load(&self, key: usize) -> HitLoad {
        if let Some(&shown) = self.shown.get(&key) {
            HitLoad::Shown(shown)
        } else if self.throws.contains_key(&key) {
            HitLoad::Throws
        } else if self.unread.contains(&key) {
            HitLoad::Unread
        } else if let Some(why) = self.refused {
            HitLoad::Refused(why)
        } else {
            HitLoad::Loaded
        }
    }

    /// What the prefetch of the batch of the subject `key` gave for its Seq-align id (its bulk
    /// load, [`ScopeWalk::bulk`]; for a prefetch that failed, what the places give after it,
    /// [`ScopeWalk::peek`]): `None` for a subject that no prefetch of the walk read (a batch whose
    /// prefetch is off, a subject without parsed Seq-ids). The one answer of the prefetch's checks
    /// (`common::database::check_prefetched_parsed`, `common::database::batch_prefetch`).
    pub fn prefetched(&self, key: usize) -> Option<HitLoad> {
        self.prefetched.get(&key).copied()
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:2935-2969
    // NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:1390-1437
    /// Whether a report of a batch without the prefetch (`ReportRead::one_by_one`: off, or
    /// failed) reads the subject `key`. A light lookup of its Seq-align id (`slen` without
    /// another field of the subject, `report::tabular::TabularSpec::light_length`) finds what the
    /// batch's prefetch put in the scope's cache ([`HitLoads::prefetched`], as the loads), but
    /// in such a batch the sequences that the data source loaded in the run, whose GC and order
    /// LOSAT does not follow, and then the subject's own record.
    pub fn read_alone(&self, key: usize) -> bool {
        self.read_alone.contains(&key)
    }

    /// The other sequence that a read of the subject `key` after report `report` gives
    /// ([`CallReads::after`]: psiblast's PSSM), where a load gives one.
    pub fn read_after(&self, report: usize, key: usize) -> Option<ShownInstead> {
        self.after_shown.get(&(report, key)).copied()
    }

    /// Whether a report (or a batch's prefetch) loads the subject `key`: by a read of its own, or
    /// as the other sequence whose record the writers write ([`HitLoads::writes_shown_records`]:
    /// its taxonomy is that sequence's); not another subject for which a load gives another
    /// sequence (`HitLoad::Shown`), which `ParsedSubject::fetched` refuses with the words of that
    /// sequence (the reads of the taxonomy texts, which BLASTN and TBLASTX make before they ask
    /// `fetched`, would refuse it with the words of its lookup); nor a subject that only a batch's
    /// failed prefetch may have loaded, which no report reaches ([`HitLoads::reached`]): no row,
    /// `Hit` or name of it is written, so none looks it up.
    pub fn loaded(&self, key: usize) -> bool {
        let loaded = match self.load(key) {
            HitLoad::Loaded | HitLoad::Refused(_) => true,
            HitLoad::Shown(_) => self.writes_shown_records(key),
            HitLoad::Throws | HitLoad::Unread => false,
        };
        loaded && !self.prefetched_only.contains(&key)
    }

    /// Whether NCBI's run reaches the report `report` of those given (counted from 0): false for
    /// every report after the one whose load ended the run, and for every report where an earlier
    /// call of the run ended it ([`RunLoads::loads`]). A batch whose first report the run does not
    /// reach is not prefetched either (NCBI's application ends before it, `CATCH_ALL`).
    pub fn reaches(&self, report: usize) -> bool {
        self.reached.is_none_or(|reached| report < reached)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp (`SeqidToOid`,
    // cited at `common::database::HeaderReader::resolution`): the lookup of a Seq-align id reads
    // the headers of its candidates (the version check) before it throws on a volume's file.
    /// Whether a report's load reaches the subject `key`, which makes the lookup of its Seq-align id
    /// (`HitLoad::Loaded` or `HitLoad::Throws`, not `HitLoad::Unread`): the lookup that throws has
    /// read the headers of the candidates before it, and posts NCBI's warning of a UTF8String
    /// header among them (`ParsedSubject::scope_utf8_read`, `common::database::load_read_texts`),
    /// though nothing of the subject is loaded (`loaded`: no taxonomy text, no definition line).
    /// A subject for which a load gives another sequence (`HitLoad::Shown`) is reached too:
    /// `ParsedSubject::fetched` refuses it, so no read is written, or the writers write that
    /// sequence's record, whose header LOSAT refuses where its load posts a text
    /// (`common::database::ShownHolders`).
    /// A subject that a batch's failed prefetch may have loaded is loaded (`load`), but reached
    /// only where a report reaches it: one shown only after a `Hit` that threw is read by no
    /// report, and its texts are the prefetch's or none (review 1 of line DBID pf7, H1).
    pub fn reached(&self, key: usize) -> bool {
        self.load(key) != HitLoad::Unread && !self.prefetched_only.contains(&key)
    }

    /// The exception that the `Hit` of the subject `key` throws (`HitLoad::Throws`).
    pub fn thrown(&self, key: usize) -> Option<&crate::ncbi_exception::NcbiException> {
        if self.shown.contains_key(&key) {
            return None;
        }
        self.throws.get(&key)
    }

    /// Whether the load of some subject throws (`HitLoad::Throws`).
    pub fn any_throws(&self) -> bool {
        !self.throws.is_empty()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:749-1096 (cited at
    // `report::shown_bioseq`: every record writer reads the Bioseq of the scope's handle)
    /// Whether the record writers write, for the subject `key`, the record of the other sequence
    /// that its loads give (`HitLoad::Shown`): the run's formats write such records
    /// (`writes_shown`: not the pairwise reports and SAM yet), every load of the subject gives
    /// that sequence and its line ([`HitLoads::written_holder`]: none gives the subject itself),
    /// and no report's load throws (`throws`: one report of a batch whose prefetch is off can
    /// fail the lookup where a later one finds the sequence that the first loaded). `false` where
    /// LOSAT refuses the subject for that sequence (`ParsedSubject::shown_reason`).
    pub fn writes_shown_records(&self, key: usize) -> bool {
        self.writes_shown
            && self.shown.contains_key(&key)
            && !self.throws.contains_key(&key)
            && matches!(self.written.get(&key), Some(Some(_)))
    }

    /// Whether the writers write another sequence's record for some subject
    /// ([`HitLoads::writes_shown_records`]).
    pub fn writes_any_shown(&self) -> bool {
        self.writes_shown
            && self
                .written
                .iter()
                .any(|(key, holder)| holder.is_some() && self.writes_shown_records(*key))
    }

    /// The other sequence, and its definition line that NCBI's Bioseq of it names, whose record
    /// the writers write for the subject `key` ([`HitLoads::writes_shown_records`]).
    pub fn written_holder(&self, key: usize) -> Option<(DbSequence, usize)> {
        if !self.writes_shown_records(key) {
            return None;
        }
        self.written.get(&key).copied().flatten()
    }
}

/// Whether the reports `one` and `other` (`other` after `one`) of `reports` are of one query
/// batch: `other` names no batch (`ReportRead::batch`: the batch of the report before) or the
/// batch of `one`, as do the reports between.
fn same_batch(reports: &[(ReportRead, ReportKeys)], one: usize, other: usize) -> bool {
    let batch = reports[one].0.batch;
    reports[one + 1..=other]
        .iter()
        .all(|(read, _)| read.batch.is_none() || read.batch == batch)
}

/// The prefetch set of the batch whose first report is `first` of `reports`: the subjects of its
/// reports (as a set), and its further subjects that no report of it shows (`call`,
/// [`CallReads::prefetched`]).
fn batch_prefetch_set(
    reports: &[(ReportRead, ReportKeys)],
    first: usize,
    call: &CallReads,
) -> (std::collections::HashSet<usize>, Vec<usize>) {
    let batch =
        (first..reports.len()).take_while(|&at| at == first || same_batch(reports, first, at));
    let mut shown = std::collections::HashSet::new();
    let mut further = Vec::new();
    for report in batch {
        shown.extend(reports[report].1.keys());
        further.extend_from_slice(call.prefetched(report));
    }
    further.retain(|key| !shown.contains(key));
    further.sort_unstable();
    further.dedup();
    (shown, further)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203 (the `CATCH_ALL` of
// every application, cited at `cli::NativeError::exception`): an exception that leaves a report
// ends the application; no later report, batch or load runs.
/// The run's record of what its reports do after a subject whose lookup throws where a batch's
/// prefetch is off (`run`, [`FormatThrow::of_run`]) and of whether a report's load ended the run
/// (an `EndsRun` subject that threw): the one record of the run's end where a program makes each
/// query's report's loads apart from the others (the texts of BLASTN's and TBLASTX's reports of a
/// batch whose prefetch is off, their `report_loads`), kept across the run's batches, so that no
/// report after the end loads anything and no batch after it is prefetched (review 2 of line DBID
/// pfoff-2, L1).
#[derive(Debug, Default)]
pub struct RunLoads {
    answer: FormatThrow,
    ended: std::cell::Cell<bool>,
    /// The walk of the reports' loads so far, after the record that it started from
    /// ([`RunLoads::loads`]).
    walk: std::cell::RefCell<Option<ScopeWalk>>,
}

impl RunLoads {
    /// The record of a run whose reports answer a throw with `answer` (the run's answer,
    /// [`FormatThrow::of_formats`]), before any report.
    pub fn new(answer: impl Into<FormatThrow>) -> Self {
        RunLoads {
            answer: answer.into(),
            ended: std::cell::Cell::new(false),
            walk: std::cell::RefCell::new(None),
        }
    }

    /// The run's answer to a subject whose lookup throws, as `new` took it (what its reports do
    /// after it, `FormatThrow::run`, or why it refuses it).
    pub fn answer(&self) -> FormatThrow {
        self.answer
    }

    /// Whether a report's load ended the run (NCBI's application ends there).
    pub fn ended(&self) -> bool {
        self.ended.get()
    }

    /// [`HitLoads::of`] of the reports `reports` after those of the run before (none after its
    /// end) and the loads `before` (the run's record before them: [`HitLoads::of`]), which records
    /// where these reports end the run. The walk goes on from the reports before where the record
    /// has not grown since (the reports of one batch, which the program records after them); a
    /// record that has grown holds the loads of the reports before (as a set), and the walk starts
    /// again from it.
    pub fn loads<'p, R: Into<ReportRead>, K: Into<ReportKeys>>(
        &self,
        before: LoadsBefore,
        reports: impl IntoIterator<Item = (R, K)>,
        parsed: impl Fn(usize) -> Option<&'p ParsedSubject>,
        identities: impl Fn(usize, usize) -> usize,
    ) -> HitLoads {
        let mut walk = self.walk_after(before);
        HitLoads::after(
            &self.ended,
            &mut walk,
            (self.answer, &CallReads::default()),
            reports,
            parsed,
            identities,
        )
    }

    /// The prefetch of a query batch whose prefetch set is `subjects` ([`ScopeWalk::bulk`]), in
    /// the run's walk after the loads `before` (as [`RunLoads::loads`], in the scope of the batch
    /// that `before` names): the answers by key, which the batch's checks read
    /// (`common::database::batch_prefetch`'s `PrefetchSubject::load`). None after the run's end,
    /// where no batch is prefetched.
    pub fn prefetched<'p>(
        &self,
        before: LoadsBefore,
        subjects: impl IntoIterator<Item = (usize, &'p ParsedSubject)>,
    ) -> std::collections::HashMap<usize, HitLoad> {
        if self.ended.get() {
            return Default::default();
        }
        let batch = before.batch();
        let mut walk = self.walk_after(before);
        walk.enter(batch);
        walk.bulk(subjects)
    }

    /// The run's walk after the loads `before`: the walk so far where the record has not grown
    /// since it started, else a new one from the record ([`RunLoads::loads`]).
    fn walk_after(&self, before: LoadsBefore) -> std::cell::RefMut<'_, ScopeWalk> {
        std::cell::RefMut::map(self.walk.borrow_mut(), |walk| {
            if walk
                .as_ref()
                .is_none_or(|walk| !walk.before.same_record(&before))
            {
                *walk = Some(ScopeWalk::new(before, self.answer.keeps_scope()));
            }
            walk.as_mut().expect("set above")
        })
    }
}

// No NCBI counterpart: LOSAT's check that no two sequences of a database are joined
// (`common::database::reports_may_refuse`).
/// Adds to `keys` the keys of the Seq-align id of `subject` (`defline_set::match_keys`): two
/// subjects that the reports join (`joined_in_reports`) share one.
pub fn join_keys(subject: &ParsedSubject, keys: &mut Vec<u64>) {
    crate::seqid::defline_set::match_keys(&subject.alignment_id, keys);
}

/// The keys of the Seq-align id of `subject` as NCBI's scope finds it
/// (`defline_set::scope_keys`).
pub fn scope_keys(subject: &ParsedSubject, keys: &mut Vec<u64>) {
    crate::seqid::defline_set::scope_keys(&subject.alignment_id, keys);
}

#[cfg(test)]
#[path = "db_subject_tests.rs"]
mod tests;
