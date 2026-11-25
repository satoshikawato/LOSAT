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
    bioseq_title_from_source, generate_bioseq_defline, BioseqFacts, Defline, SourceTaxname,
    UnknownEncoding,
};
use super::html_ids::{self, HtmlHeading};

/// A reported database sequence with parsed Seq-ids: its definition lines, its Seq-align id
/// and the definition line of its Bioseq.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ParsedSubject {
    set: DeflineSet,
    alignment_id: ParsedSeqId,
    /// `DeflineSet::bioseq_defline` of `alignment_id`.
    target: usize,
    /// Whether NCBI's formatters find no Bioseq by the Seq-align id whatever its lookup gives
    /// (`new`, `fetched`).
    unfetchable: bool,
    /// What NCBI's lookup of the Seq-align id in the database finds where that is not this
    /// sequence (`found_instead`, `fetched`).
    found_instead: Option<FoundInstead>,
    /// Another sequence that NCBI's scope can show for this one by a weak match where its
    /// formatters read the batch's sequences one by one (`weak_match`, `fetched`).
    weak_match: Option<WeakMatch>,
    /// Another sequence (`DbSequence`) that holds an id NCBI's scope takes for the Seq-align id,
    /// besides the sequence that the lookup gives first or where it gives none
    /// (`FoundInstead::Sequence`, `Nothing`): NCBI's scope can read it for this one where it has
    /// read it before (`also_held`, `pssm_read`).
    also_held: Option<DbSequence>,
    /// The sequence of the database that this subject is (`in_database`; `None` for a `-subject`
    /// record), by which the run's loads name it (`HitLoads`, `ScopeHolders`).
    sequence: Option<DbSequence>,
    /// Where the lookup of the Seq-align id fails on a file (`FoundInstead::Failed`): the other
    /// sequences that hold an id NCBI's scope takes for it (`scope_holders`, `held_before`).
    failed_holders: ScopeHolders,
    /// The sequence (`DbSequence`) that the lookup of the Seq-align id gives first where that is
    /// another one whose header record holds a UTF8String where a VisibleString is expected, read
    /// before the run posted NCBI's warning (`scope_utf8`): NCBI's scope reads that header
    /// wherever it resolves the id, and the warning comes there (`fetched`, `pssm_read`).
    scope_utf8: Option<DbSequence>,
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
    /// The sequence itself first (or the same record), and later another sequence of the
    /// database (`DbSequence`) that holds an id NCBI's scope takes for the Seq-align id: where
    /// both are reported, NCBI's formatters can show that one for it.
    AlsoHeld(DbSequence),
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
    /// LOSAT's own list of the other sequences that hold the Seq-align id fails on a file that
    /// NCBI's lookup does not read (`why`: what LOSAT cannot list, with the error), so whether
    /// another sequence holds an id that NCBI's scope takes for it is unknown: the list of every
    /// OID of the lookup (`SeqidLookup::seqid_to_oids` with `multi`) on an ISAM of a volume after
    /// the first one that holds the id (a version 4 string ISAM that does not open), where NCBI's
    /// lookup has returned (review 3 of line DBID loadord, H3), or the searches of the keys that
    /// differ from the Seq-align id's in the case of its letters (`SeqidLookup::case_folded_oids`
    /// and `case_variant_oids`: version 5). `instead`: the sequence that the lookup gives where it
    /// is another one (as `Sequence`), else the lookup gives the subject or the same record.
    HoldersUnknown {
        instead: Option<DbSequence>,
        why: String,
    },
    /// The lookup throws on the id itself before it reads a file (`SeqidLookup::throws_on_id`: a
    /// `gnl|BL_ORD_ID|0` that makeblastdb stores with a text tag), with this message: NCBI's
    /// formatters stop, its prefetch catches the exception with its warning
    /// (`ParsedSubject::prefetch_fails`).
    ThrowsOnId(String),
}

/// What NCBI's scope can show for a subject by a weak match where its formatters read the
/// reported sequences of a batch one by one (`PRE_FETCH_SEQS_LIMIT` turns the batch's prefetch
/// off; the prefetch loads each by its own id): `ParsedSubject::weak_match`, refused by
/// `fetched` with `one_by_one`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum WeakMatch {
    /// Another sequence of the database (`DbSequence`) that holds an id the Seq-align id
    /// matches weakly (`seqid::id_scope::weak_keys`): where that one is read first, NCBI's scope
    /// shows it for the subject.
    Held(DbSequence),
    /// A Seq-align id with weak keys that a version 5 database cannot look up: one with a name,
    /// or a `pir` or `prf` one (`common::database::HeaderReader::weakly_held`).
    Unlisted,
    /// The lookup of the sequences that could match weakly failed, with this message (LOSAT's
    /// lookup of the name alone; NCBI's scope makes none).
    Failed(String),
}

// NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:832-880 (cited at
// `common::database::HeaderReader::found_instead`)
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:379-424
// NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/cached_sequence.cpp:98-110
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2295-2314
// NCBI's scope answers a load by an id before its BLAST database data loader: from its index of
// the Bioseqs it holds, by the id's own handle, then by every handle that matches it, weak matches
// included (`GetMatchingHandles` with `eAllowWeakMatch`); a Bioseq's ids there are those of the
// definition line that its load targets, the line of its own Seq-align id
// (`common::database::Lines::Loaded`); the formatter empties that index after each query batch
// (`ResetScopeHistory`), but for outfmt 5. Then the loader answers from its `m_Ids`, the exact
// handle of every Seq-id of every definition line of each sequence it loaded (`RegisterIds`), kept
// for the run, before it looks the id up (`SeqidToOid`), which is what fails on a file (review 2
// of line DBID pfoff-4, H1 and M1; review 3, M1: `gi|7051|ref|XP_900001.1|` ^A
// `gi|7052|ref|NP_100001.1|` holds `ref|NP_100001.1|` by the loader and not `ref|NP_100001|` or
// `emb|NP_100001.1|` by the scope, which knows its first line; with `ref|XP_900002.1|` ^A
// `gi|7053|ref|NP_100001.1|`, whose Seq-align id is the gi of its second line, the scope holds
// both).
/// The other sequences of the database that hold an id by which NCBI's run may give one of them
/// for the Seq-align id of a subject whose lookup fails on a file (`FoundInstead::Failed`) without
/// that lookup, in the two places that answer before it, each with its own keys and lifetime, as
/// `common::database::HeaderReader::found_instead` lists them in the LMDB files that have the
/// accession table (a file without it holds no text id): where the run loaded one of them before a
/// load of the subject, within that place's lifetime, NCBI gives that one for it
/// (`ParsedSubject::held_before`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ScopeHolders {
    /// The holders of an id of the Seq-align id's own `CSeq_id_Handle`
    /// (`seqid::mask_handle::mask_handle_key`) in any of their definition lines, which the BLAST
    /// database data loader gives for it from the ids of the sequences that it loaded (`m_Ids`)
    /// for the rest of the run, in every format.
    pub loader: Holders,
    /// The holders of an id that NCBI's scope finds for it in its index, in the definition line by
    /// which their load made their Bioseq (the line of their own Seq-align id,
    /// `common::database::Lines::Loaded`; review 3 of line DBID pfoff-4, M1): an id that it takes
    /// for it (`defline_set::scope_keys`: the accession and version in any alternative and letter
    /// case, a name of an id without an accession; the loader's holders by that line among them)
    /// or one that it matches weakly (`id_scope::weak_keys`: the versions of an accession without
    /// one, a name with or without an accession, in any alternative), which the scope holds from
    /// their load to the end of the query batch (`ResetScopeHistory`), to the end of the run in
    /// outfmt 5 (`RunThrow::keeps_scope`).
    pub scope: Holders,
}

/// A list of [`ScopeHolders`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Holders {
    /// These sequences, in the order of their ordinal ids (none: no other sequence holds such an
    /// id).
    Listed(Vec<DbSequence>),
    /// LOSAT cannot make the list, for this reason (written only in LOSAT's refusal, as the `why`
    /// of `FoundInstead::HoldersUnknown`): it fails on a file, or it asks for the weak matches of
    /// a name, or of a `pir` or `prf` id, which a version 5 database cannot list
    /// (`WeakMatch::Unlisted`).
    Unknown(String),
}

impl Default for Holders {
    fn default() -> Self {
        Holders::Listed(Vec::new())
    }
}

impl Holders {
    /// The first of these holders that `loaded` says was loaded; where LOSAT cannot list them,
    /// `Some(None)` where any sequence with parsed Seq-ids may have been loaded (`any`).
    fn held(&self, loaded: impl Fn(DbSequence) -> bool, any: bool) -> Option<Option<DbSequence>> {
        match self {
            Holders::Listed(holders) => holders
                .iter()
                .copied()
                .find(|&holder| loaded(holder))
                .map(Some),
            Holders::Unknown(_) => any.then_some(None),
        }
    }
}

/// What a run loaded before the reports that [`HitLoads::of`] walks
/// (`common::database::TaxdbRun::loads_before`, the one producer): the sequences of the database
/// (ordinal ids, in their order) that the run's record of its loads holds
/// (`common::database::TaxdbRun::note_loaded`: the batches, psiblast's rounds and their PSSM reads
/// before; the record is kept where a lookup can fail on a file, `TaxdbRun::prefetch_can_stop`),
/// those of them loaded in the query batch `batch` (psiblast's rounds before of the batch's query,
/// which share its scope), and whether calls of the run before this one loaded sequences that no
/// record here holds (`unrecorded`: `api::run_batches`' calls after the first, each a query batch
/// of its own).
#[derive(Clone, Debug, Default)]
pub struct LoadsBefore {
    recorded: Vec<u32>,
    batch: Option<usize>,
    in_batch: Vec<u32>,
    unrecorded: bool,
}

impl LoadsBefore {
    /// Nothing loaded before: the reports are the run's first.
    pub const NONE: LoadsBefore = LoadsBefore {
        recorded: Vec::new(),
        batch: None,
        in_batch: Vec::new(),
        unrecorded: false,
    };

    /// The loads that the record holds: `recorded` (ordinal ids, sorted), of which `in_batch` in
    /// the query batch `batch`, after calls whose loads no record holds where `unrecorded`.
    pub fn new(
        recorded: Vec<u32>,
        (batch, in_batch): (Option<usize>, Vec<u32>),
        unrecorded: bool,
    ) -> Self {
        debug_assert!(recorded.is_sorted() && in_batch.is_sorted());
        LoadsBefore {
            recorded,
            batch,
            in_batch,
            unrecorded,
        }
    }

    /// These loads, then `run` (ordinal ids) of other query batches and `in_batch` of this
    /// record's batch: the loads of a walk before the next one (psiblast's PSSM reads after the
    /// round's reports, `ListedTitles::check_pssm_holds` of BLASTP).
    pub fn and_loads(mut self, run: Vec<u32>, in_batch: Vec<u32>) -> Self {
        self.recorded.extend(run.iter().chain(&in_batch));
        self.recorded.sort_unstable();
        self.recorded.dedup();
        self.in_batch.extend(in_batch);
        self.in_batch.sort_unstable();
        self.in_batch.dedup();
        self
    }

    /// No record, after calls whose loads no record holds where `unrecorded`.
    pub fn after_calls(unrecorded: bool) -> Self {
        LoadsBefore {
            unrecorded,
            ..Self::NONE
        }
    }

    /// These loads, also after calls whose loads no record holds where `unrecorded`.
    pub fn with_calls_before(self, unrecorded: bool) -> Self {
        LoadsBefore {
            unrecorded: self.unrecorded || unrecorded,
            ..self
        }
    }

    /// The query batch whose loads `in_batch` are (`None`: not named).
    pub fn batch(&self) -> Option<usize> {
        self.batch
    }
}

/// The sequences with parsed Seq-ids of one of the places of [`SeenLoads`].
#[derive(Debug, Default)]
struct Loads {
    /// The sequences (ordinal ids) loaded so far.
    ordinals: std::collections::HashSet<u32>,
    /// Whether a sequence with parsed Seq-ids was loaded so far (also one that names no database
    /// sequence, `ParsedSubject::sequence`).
    any: bool,
}

impl Loads {
    /// Adds the loads `ordinals`.
    fn extend(&mut self, ordinals: &[u32]) {
        self.ordinals.extend(ordinals);
        self.any |= !ordinals.is_empty();
    }

    /// Adds a load of the subject `subject`.
    fn note(&mut self, subject: &ParsedSubject) {
        self.any = true;
        if let Some(sequence) = subject.sequence {
            self.ordinals.insert(sequence.ordinal());
        }
    }
}

/// The sequences that a walk of [`HitLoads`] has seen loaded, in NCBI's order, in the two places
/// of [`ScopeHolders`]: from the run's start, the run's record before the walk ([`LoadsBefore`])
/// and the walk's own loads (the data loader's), and of those, the ones of the current query
/// batch (the scope's, emptied where a report of another batch comes; [`SeenLoads::enter`]).
/// [`ParsedSubject::held_before`] asks them.
#[derive(Debug, Default)]
struct SeenLoads {
    run: Loads,
    scope: Loads,
    /// The query batch of `scope` (`None`: none named yet).
    batch: Option<usize>,
    /// Whether calls of the run before loaded sequences that no record holds: any may be loaded.
    unrecorded: bool,
}

impl SeenLoads {
    /// The loads `before`, at the start of a walk.
    fn of(before: LoadsBefore) -> Self {
        let mut seen = SeenLoads::default();
        seen.add(before);
        seen
    }

    /// Adds the loads `before` (the run's record, which may hold loads seen already).
    fn add(&mut self, before: LoadsBefore) {
        self.run.extend(&before.recorded);
        self.unrecorded |= before.unrecorded;
        self.enter(before.batch);
        self.scope.extend(&before.in_batch);
    }

    /// The scope of a report of the query batch `batch` (`None`: of the batch of the report
    /// before): where it is another batch than the scope's, NCBI's formatter emptied the scope at
    /// the end of that one (`ResetScopeHistory`; the walk keeps the whole run's loads instead where
    /// the run's scope is kept, `RunThrow::keeps_scope`).
    fn enter(&mut self, batch: Option<usize>) {
        let Some(batch) = batch else {
            return;
        };
        if self.batch.is_some_and(|current| current != batch) {
            self.scope = Loads::default();
        }
        self.batch = Some(batch);
    }

    /// Records a load of the subject `subject` (nothing where it has no parsed Seq-ids: such a
    /// sequence holds no id that a lookup can fail on).
    fn note(&mut self, subject: Option<&ParsedSubject>) {
        if let Some(subject) = subject {
            self.run.note(subject);
            self.scope.note(subject);
        }
    }

    /// [`ParsedSubject::held_before`] of `subject` after these loads, in a run whose scope NCBI
    /// keeps for the whole run where `keeps_scope` (outfmt 5): the calls before loaded any
    /// sequence into the loader, and into the scope only where it is kept. A holder that either
    /// place gives is noted in the batch's scope from then on (the loader's answer loads its
    /// Bioseq there, of the definition line of the failing id, which LOSAT notes as a load of the
    /// holder: no output depends on it, `ParsedSubject::fetched` refuses every hold); one that
    /// LOSAT does not name counts as any sequence.
    fn held(&mut self, subject: &ParsedSubject, keeps_scope: bool) -> Option<Option<DbSequence>> {
        let held = self.held_now(subject, keeps_scope);
        if let Some(holder) = held {
            let ordinals = holder.map(DbSequence::ordinal);
            for place in [&mut self.run, &mut self.scope] {
                place.any = true;
                place.ordinals.extend(ordinals);
            }
        }
        held
    }

    /// The answer of [`SeenLoads::held`], which records nothing.
    fn held_now(&self, subject: &ParsedSubject, keeps_scope: bool) -> Option<Option<DbSequence>> {
        let (scope, scope_unrecorded) = if keeps_scope {
            (&self.run, self.unrecorded)
        } else {
            (&self.scope, false)
        };
        subject.held_before(
            (
                |holder: DbSequence| {
                    self.unrecorded || self.run.ordinals.contains(&holder.ordinal())
                },
                self.run.any || self.unrecorded,
            ),
            (
                |holder: DbSequence| scope_unrecorded || scope.ordinals.contains(&holder.ordinal()),
                scope.any || scope_unrecorded,
            ),
        )
    }
}

/// Records in `seen` and `loaded` the loads of the subjects `keys` (`parsed`: their parsed
/// records) by a walk of [`HitLoads`].
fn note_loads<'p>(
    seen: &mut SeenLoads,
    loaded: &mut std::collections::HashSet<usize>,
    keys: impl IntoIterator<Item = usize>,
    parsed: &impl Fn(usize) -> Option<&'p ParsedSubject>,
) {
    for key in keys {
        seen.note(parsed(key));
        loaded.insert(key);
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
    /// instead: `FoundInstead::Sequence` and `HoldersUnknown` with another sequence that the
    /// lookup gives (that sequence where one of its ids compares equal to the Seq-align id, else
    /// none), `AlsoHeld` (`or_itself`: the subject or that sequence, as the ids of the blobs read
    /// before in the run decide, `CBlastDbDataLoader::m_Ids`), and the same record whose lookups
    /// write other texts (`None`); `or_none`: or nothing, where the Bioseq of another sequence
    /// that holds the id, which the scope read before, can stand for it without a load (one that
    /// `also_held` names, or one that LOSAT cannot list, `HoldersUnknown`); `reason` as `fetched`
    /// gives it.
    Another {
        other: Option<DbSequence>,
        or_itself: bool,
        or_none: bool,
        reason: String,
    },
    /// The sequence itself, or nothing where the Bioseq of another sequence that holds the id,
    /// which the scope read before and which LOSAT cannot list, stands for it
    /// (`FoundInstead::HoldersUnknown` where the lookup gives the subject or the same record);
    /// `reason` as `fetched` gives it.
    ItselfOrNothing { reason: String },
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
    /// write as NCBI: a stored title or Seq-id with a byte outside ' ' ..= '~' (NCBI writes an
    /// error message each time it reads the definition line) and a definition line without a
    /// Seq-id (NCBI's formatters take the best of none; makeblastdb writes none such).
    pub fn new(set: DeflineSet) -> Result<Self, String> {
        if set
            .deflines
            .iter()
            .any(|defline| defline.replaced_title_bytes > 0 || defline.replaced_id_bytes > 0)
        {
            return Err("has a stored title or Seq-id with a byte outside ' ' to '~' (NCBI BLAST+ shows it as '#' and writes an error message each time it reads the definition line)".into());
        }
        if set.deflines.iter().any(|defline| defline.ids.is_empty()) {
            return Err("has a definition line without a Seq-id".into());
        }
        let Some(alignment_id) = set.alignment_id().cloned() else {
            return Err("has no Seq-id".into());
        };
        let Some(target) = set.bioseq_defline(&alignment_id) else {
            return Err("has no definition line that holds its Seq-align id".into());
        };
        // A giim Seq-align id is never found again: NCBI's handles never take it for one with an
        // equal id (`algorithm::psiblast::rounds::HandleIndex`), so its scope finds no Bioseq by it
        // (NCBI BLAST+ 2.17.0 writes `Unknown`; review 2 of line DBID pssmrows, H1).
        let unfetchable = matches!(alignment_id, ParsedSeqId::Giim(_))
            || set
                .display_defline(&alignment_id)
                .is_some_and(|shown| shown.index != target);
        Ok(Self {
            set,
            alignment_id,
            target,
            unfetchable,
            found_instead: None,
            weak_match: None,
            also_held: None,
            sequence: None,
            failed_holders: ScopeHolders::default(),
            scope_utf8: None,
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
            found_instead: None,
            weak_match: None,
            also_held: None,
            sequence: None,
            failed_holders: ScopeHolders::default(),
            scope_utf8: None,
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
    pub fn found_instead(mut self, found: Option<FoundInstead>) -> Self {
        self.found_instead = found;
        self
    }

    /// The subject with `weak`, what NCBI's scope can show for it by a weak match where its
    /// formatters read the batch's sequences one by one (`fetched`'s `one_by_one`).
    pub fn weak_match(mut self, weak: Option<WeakMatch>) -> Self {
        self.weak_match = weak;
        self
    }

    /// The subject with `holder`, another sequence that holds an id NCBI's scope takes for the
    /// Seq-align id where the lookup gives a third sequence first or none (`pssm_read`).
    pub fn also_held(mut self, holder: Option<DbSequence>) -> Self {
        self.also_held = holder;
        self
    }

    /// The subject as the sequence `sequence` of its database, by which the run's loads name it
    /// (`HitLoads`).
    pub fn in_database(mut self, sequence: DbSequence) -> Self {
        self.sequence = Some(sequence);
        self
    }

    /// The subject with `holders`, the other sequences that hold an id NCBI's scope takes for the
    /// Seq-align id where its lookup fails on a file (`held_before`).
    pub fn scope_holders(mut self, holders: ScopeHolders) -> Self {
        self.failed_holders = holders;
        self
    }

    /// The holders that `scope_holders` gave. Read by tests only: the walks ask `held_before`.
    #[cfg(test)]
    pub fn failed_holders(&self) -> &ScopeHolders {
        &self.failed_holders
    }

    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:379-396 (cited
    // at `common::database::stop_lookup_unsure`)
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_info.cpp:832-880 (cited at
    // `common::database::HeaderReader::found_instead`)
    // A load by the Seq-align id asks the scope first, which finds the Bioseq of a sequence that it
    // holds with an id of the Seq-align id's handle, or else with an id that it takes for that id
    // or that the id matches weakly (it holds what the query batch loaded, in outfmt 5 what the
    // run loaded), then the BLAST database data loader, which keeps the handle of every Seq-id of
    // each sequence that it loaded in the run (`RegisterIds`) and gives that sequence for the id
    // before it looks the id up: a lookup that would fail on a file is not made, nothing throws,
    // and NCBI 2.17.0 writes that sequence for the subject (review 1 of line DBID pfoff-4, H1 and
    // H2: `gi|5001|ref|NP_100001.1|`, loaded in an earlier query's report, for `ref|NP_100001.1|`
    // of a list with a database without the accession table, in outfmt 0, 5, 6, 7, 10, 13 to 16
    // and 18; review 2, H1: `ref|NP_100001|` and `gb||LOCUSN` by a weak match in the same batch;
    // M1: `emb|NP_100001.1|` of an earlier batch is not held, the same lookup before that load
    // throws).
    /// Where the lookup of the Seq-align id fails on a file with NCBI's exception
    /// (`lookup_throws`), what NCBI may give for this subject without the lookup at a load after
    /// the loads before it ([`ScopeHolders`]): `run` tells of a sequence whether the run may have
    /// loaded it before, and whether it may have loaded any with parsed Seq-ids (the data loader's
    /// loads, for its holders), `scope` the same of the loads that NCBI's scope holds there (for
    /// its holders): `Some(Some(holder))` for a holder that may be held (the first such, the
    /// loader's first), `Some(None)` where LOSAT cannot list the holders of a place and some
    /// sequence may have been loaded into it; `None` where NCBI makes the lookup (it throws), also
    /// for every other subject.
    pub fn held_before(
        &self,
        run: (impl Fn(DbSequence) -> bool, bool),
        scope: (impl Fn(DbSequence) -> bool, bool),
    ) -> Option<Option<DbSequence>> {
        self.lookup_throws()?;
        let holders = &self.failed_holders;
        holders
            .loader
            .held(run.0, run.1)
            .or_else(|| holders.scope.held(scope.0, scope.1))
    }

    /// The reason of LOSAT's refusal of a subject that NCBI may find without its lookup, which
    /// fails on a file (`HitLoad::Held` of `held_before`): the words of another sequence that
    /// holds the Seq-align id (`FoundInstead::AlsoHeld`), as for a lookup that does not fail, or of
    /// the holders that LOSAT cannot list (`FoundInstead::HoldersUnknown`, the reason of the first
    /// place whose list it cannot make).
    fn held_reason(&self, holder: Option<DbSequence>) -> String {
        let found = match (holder, self.unlisted_holders()) {
            (Some(other), _) => FoundInstead::AlsoHeld(other),
            (None, why) => FoundInstead::HoldersUnknown {
                instead: None,
                why: why.unwrap_or_else(|| "the sequences that hold it".to_owned()),
            },
        };
        self.found_reason(&found, false, None).unwrap_or_default()
    }

    /// Why LOSAT cannot list the sequences that NCBI may give for the Seq-align id without its
    /// failing lookup: the reason of the first place (the loader's handles, the scope's ids) whose
    /// list it cannot make; `None` where it lists both.
    fn unlisted_holders(&self) -> Option<String> {
        [&self.failed_holders.loader, &self.failed_holders.scope]
            .into_iter()
            .find_map(|holders| match holders {
                Holders::Unknown(why) => Some(why.clone()),
                Holders::Listed(_) => None,
            })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:411-451 (cited at
    // `algorithm::psiblast::pssm::input_msa`)
    /// The reason of LOSAT's refusal of psiblast's PSSM read of this subject, whose lookup fails on
    /// a file, where NCBI's scope gives its Seq-align id the sequence `holder`, loaded before
    /// (`None`: a holder that LOSAT cannot list, `held_before`), and the read takes that sequence's
    /// letters (its letters reach the end of an included Seq-align's range,
    /// `ListedTitles::check_pssm_holds` of BLASTP). Not the words of a report that shows the
    /// holder (`held_reason`): the PSSM reads by the id, also beyond `-max_target_seqs` and in
    /// the formats that show no subject (review 5 of line DBID pfoff-4, L1).
    pub fn pssm_held_reason(&self, holder: Option<DbSequence>) -> String {
        let id = String::from_utf8_lossy(&self.alignment_fasta()).into_owned();
        match holder {
            Some(other) => format!(
                "has the Seq-align id {id}, whose lookup in the database fails, and which database sequence {other}, loaded before, also holds: psiblast's PSSM reads the subject by that id (also beyond -max_target_seqs), and NCBI BLAST+ 2.17.0 takes that sequence's letters for it"
            ),
            None => {
                let why = self
                    .unlisted_holders()
                    .unwrap_or_else(|| "the sequences that hold it".to_owned());
                format!(
                    "has the Seq-align id {id}, whose lookup in the database fails, for which LOSAT cannot list {why}: psiblast's PSSM reads the subject by that id (also beyond -max_target_seqs), and where such a sequence was loaded before, NCBI BLAST+ 2.17.0 takes its letters for it"
                )
            }
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
    pub fn scope_utf8(mut self, other: Option<DbSequence>) -> Self {
        self.scope_utf8 = other;
        self
    }

    /// `scope_utf8`: the other sequence whose header NCBI reads with its warning where it
    /// resolves this subject's Seq-align id.
    pub fn scope_utf8_read(&self) -> Option<DbSequence> {
        self.scope_utf8
    }

    /// Whether NCBI's formatters find no Bioseq by the Seq-align id whatever its lookup gives
    /// (see `new`), so that no lookup is needed. Read by tests only: the engine asks
    /// `common::database::scope_bioseq`, which looks it up first.
    #[cfg(test)]
    pub fn unfetchable(&self) -> bool {
        self.unfetchable
    }

    /// What NCBI's lookup of the Seq-align id finds where that is not this sequence
    /// (`found_instead`).
    pub fn found(&self) -> Option<&FoundInstead> {
        self.found_instead.as_ref()
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
        self.set.first_source_taxid(Some(&self.alignment_id))
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
    /// the reports that show it refuse, `generated_title_fetched`). `one_by_one`: whether NCBI's
    /// formatters may read the subject one by one (a batch that reports it has its prefetch
    /// off, `common::database::ParsedShown::one_by_one`), which alone refuses a sequence that
    /// another can stand for by a weak match (`WeakMatch`), after what the lookup finds.
    /// `load`: what the reports do at the subject's `Hit` (`HitLoads::load`, the one answer): a
    /// subject whose `Hit` throws (`HitLoad::Throws`: its lookup fails on a file with NCBI's
    /// exception in a batch whose prefetch is off, and every format of the run writes it as NCBI
    /// does, `FormatThrow::of_run`: BLAST XML (outfmt 5), the XML2 and JSON documents of outfmt 13
    /// to 16, the tables of outfmt 6, 7, 10 and 20, the organism report of outfmt 18, the pairwise
    /// reports of outfmt 0 to 4) is not refused for what the lookup finds (it finds nothing: it
    /// throws before any definition line is read), and one that no report loads
    /// (`HitLoad::Unread`: shown only after such a `Hit`, or as a copy of the record in it) is not
    /// refused at all; one whose lookup fails but that NCBI's scope may find without it, by
    /// another sequence that holds an id it takes for the Seq-align id and that the run loaded
    /// before (`HitLoad::Held`, `held_before`), is refused as such a holder is where the lookup
    /// does not fail (`FoundInstead::AlsoHeld`); every other subject (`HitLoad::Loaded`,
    /// `HitLoad::Refused`) is refused as above, a lookup that fails on a file among them (outfmt
    /// 17 here, before any byte is written; the reason of a run whose formats refuse it,
    /// `HitLoad::Refused`).
    pub fn fetched(&self, taxonomy: bool, one_by_one: bool, load: HitLoad) -> Result<(), String> {
        match load {
            HitLoad::Unread => return Ok(()),
            HitLoad::Held(holder) => return Err(self.held_reason(holder)),
            HitLoad::Throws => {}
            HitLoad::Loaded => self.found_fetched(taxonomy, None)?,
            HitLoad::Refused(refused) => self.found_fetched(taxonomy, Some(refused))?,
        }
        self.weakly_fetched(one_by_one)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:357-360, 412-417
    // Where the batch's prefetch is off, the `Hit` loads its subject by the Seq-align id: the
    // lookup throws on the file there, and again in the `catch` (`GetLength` loads by the same id),
    // outside any `catch` of the hit; the report function of outfmt 15 and 16 catches it, 13 and 14
    // end the program with it (probe of line DBID mixdb, `PU.off`; DECISIONS 62 (2)).
    /// The exception that NCBI's lookup of the Seq-align id throws where a formatter loads the
    /// subject one by one (a batch whose prefetch is off): a lookup that fails on a file with
    /// NCBI's `CSeqDBException` whose place LOSAT knows (`FoundInstead::Failed`,
    /// `LookupFailure::exception`). Where the batch prefetches, its prefetch throws first and
    /// stops the run (`prefetch_stops`).
    pub fn lookup_throws(&self) -> Option<crate::ncbi_exception::NcbiException> {
        match &self.found_instead {
            Some(FoundInstead::Failed(failure)) => failure.exception(),
            _ => None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:440-450 (cited at
    // `algorithm::psiblast::pssm::input_msa`)
    /// LOSAT's refusal of what psiblast's PSSM reads by the Seq-align id of a subject that it
    /// includes (the `CSeqVector` of each Seq-align's range in NCBI's scope), where LOSAT does not
    /// read it as NCBI: what the scope gets is `common::database::scope_bioseq`'s, and this
    /// refuses another sequence that also holds an id the scope takes for the Seq-align id, which
    /// the scope reads for it where it has read it before (`FoundInstead::AlsoHeld` where the
    /// lookup gives the subject first, `also_held` where it gives another or none), a lookup that
    /// fails on a file with LOSAT's error (`FoundInstead::Failed`; with NCBI's exception, which the
    /// read catches, `prefetch_stops`: the caller reads no letters), the same record whose read writes other
    /// taxonomy texts (`CopyApart::TaxonomyTexts`), and with `one_by_one` a weak match
    /// (`fetched`). The reasons are what the lookup finds, also for a subject whose definition
    /// line the formatters do not find (`unfetchable`), which the PSSM reads by the scope as any
    /// other (review 2 of line DBID pssmread, M2).
    pub fn pssm_read(&self, one_by_one: bool) -> Result<(), String> {
        let id = || String::from_utf8_lossy(&self.alignment_fasta()).into_owned();
        match &self.found_instead {
            Some(
                found @ (FoundInstead::AlsoHeld(_)
                | FoundInstead::HoldersUnknown { .. }
                | FoundInstead::Failed(_)
                | FoundInstead::CopyApart(_, CopyApart::TaxonomyTexts)),
            ) => return Err(self.found_reason(found, false, None).unwrap_or_default()),
            Some(FoundInstead::Sequence(_) | FoundInstead::Nothing) => {
                if let Some(holder) = self.also_held {
                    return Err(format!(
                        "has the Seq-align id {}, which database sequence {holder} also holds (an id that NCBI BLAST+'s scope takes for it) besides what the lookup gives: NCBI BLAST+ 2.17.0 reads that sequence for it where it has read it before, as its reads and the blobs' layout decide",
                        id()
                    ));
                }
            }
            _ => {}
        }
        self.weakly_fetched(one_by_one)
    }

    /// Another sequence (`DbSequence`) that holds an id NCBI's scope takes for the Seq-align id
    /// besides what the lookup gives (`also_held`, `FoundInstead::AlsoHeld`), whose Bioseq the
    /// scope gets for this one where it holds it from an earlier read
    /// (`common::database::scope_bioseq`).
    pub fn other_holder(&self) -> Option<DbSequence> {
        match &self.found_instead {
            Some(FoundInstead::AlsoHeld(other)) => Some(*other),
            Some(FoundInstead::Sequence(_) | FoundInstead::Nothing) => self.also_held,
            _ => None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1437-1439
    /// Whether the Bioseq that `GetBioseqNoData` builds of this sequence for its Seq-align id
    /// holds that id: the ids of the definition line its filter keeps (`target`,
    /// `DeflineSet::bioseq_defline`) hold one that compares equal to it (the scope's test on the
    /// Bioseq's ids, `common::database::scope_bioseq`). Not where that line is found by the
    /// content of a local id alone (`lcl|777` before `gi|777`), nor for a giim id.
    pub fn bioseq_holds_alignment_id(&self) -> bool {
        !matches!(self.alignment_id, ParsedSeqId::Giim(_))
            && self
                .target()
                .ids
                .iter()
                .any(|id| compare(id, &self.alignment_id) == Sic::Yes)
    }

    /// `fetched`'s refusal of a weak match where NCBI's formatters may read the subject one by
    /// one (`one_by_one`).
    fn weakly_fetched(&self, one_by_one: bool) -> Result<(), String> {
        let id = || String::from_utf8_lossy(&self.alignment_fasta()).into_owned();
        match self.weak_match.as_ref().filter(|_| one_by_one) {
            None => Ok(()),
            Some(WeakMatch::Held(other)) => Err(format!(
                "has the Seq-align id {}, which can match an id of database sequence {other} weakly (an accession with a version for one without, a name with or without an accession, a pir or prf accession as a name): in a batch whose prefetch is off (PRE_FETCH_SEQS_LIMIT, or a format that NCBI BLAST+ 2.17.0 does not prefetch for), where NCBI BLAST+ 2.17.0 reads the reported sequences one by one and that sequence first, it can write that sequence's names and letters for it",
                id()
            )),
            Some(WeakMatch::Unlisted) => Err(format!(
                "has the Seq-align id {}, which ids of other sequences can match weakly in any letter case (a name with or without an accession, a pir or prf accession as a name) and which a version 5 database does not look up: in a batch whose prefetch is off (PRE_FETCH_SEQS_LIMIT, or a format that NCBI BLAST+ 2.17.0 does not prefetch for), where NCBI BLAST+ 2.17.0 reads the reported sequences one by one, it can write the names and letters of such a sequence read before",
                id()
            )),
            Some(WeakMatch::Failed(message)) => Err(format!(
                "has the Seq-align id {}, for which LOSAT cannot look up the sequences that can match it weakly ({message}): in a batch whose prefetch is off (PRE_FETCH_SEQS_LIMIT, or a format that NCBI BLAST+ 2.17.0 does not prefetch for), where NCBI BLAST+ 2.17.0 reads the reported sequences one by one, it can write the names and letters of such a sequence read before",
                id()
            )),
        }
    }

    /// `fetched` but for the weak matches: what the lookup of the Seq-align id finds
    /// (`found_instead`), whatever the batch's prefetch.
    fn found_fetched(&self, taxonomy: bool, refused: Option<ThrowRefusal>) -> Result<(), String> {
        if self.unfetchable {
            return Err(format!(
                "has the Seq-align id {}, by which NCBI BLAST+'s formatters do not find the sequence (NCBI BLAST+ 2.17.0 writes Unknown or another sequence's names, or stops)",
                String::from_utf8_lossy(&self.alignment_fasta())
            ));
        }
        match &self.found_instead {
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
            FoundInstead::AlsoHeld(other) => format!(
                "has the Seq-align id {id}, which database sequence {other} also holds (where both are reported, NCBI BLAST+ 2.17.0 can write that sequence's names and letters for it)"
            ),
            FoundInstead::HoldersUnknown { instead, why } => format!(
                "has the Seq-align id {id}, {}for which LOSAT cannot list {why}: where such a sequence is also reported, NCBI BLAST+ 2.17.0 can write its names and letters for it",
                instead.map_or_else(String::new, |other| format!(
                    "by which NCBI BLAST+'s lookup gives database sequence {other} first (NCBI BLAST+ 2.17.0 writes that sequence's names and letters, or stops), and "
                ))
            ),
            FoundInstead::Failed(failure) => format!(
                "has the Seq-align id {id}, whose lookup in the database fails ({failure}; {})",
                failure.what_ncbi_does(refused)
            ),
            FoundInstead::ThrowsOnId(message) => format!(
                "has the Seq-align id {id}, whose lookup in the database fails ({message}; NCBI BLAST+ 2.17.0 stops)"
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
    pub fn prefetched(&self) -> Prefetched {
        let reason = || {
            self.found_instead
                .as_ref()
                .and_then(|found| self.found_reason(found, false, None))
                .unwrap_or_default()
        };
        match &self.found_instead {
            None
            | Some(FoundInstead::OtherTaxids(_) | FoundInstead::CopyApart(_, CopyApart::Title)) => {
                Prefetched::Itself
            }
            Some(FoundInstead::Nothing) => Prefetched::Nothing,
            // The prefetch catches the exception of a failed lookup and writes its warning, or
            // rethrows a database error (`CSeqDBException`) and stops; the formats that read the
            // Bioseq stop (`fetched`).
            Some(FoundInstead::ThrowsOnId(message)) => Prefetched::Failed(format!(
                "has the Seq-align id {}, whose lookup in the database fails ({message}; NCBI BLAST+ 2.17.0's prefetch writes Error pre-fetching sequence data, and a report that reads it stops)",
                String::from_utf8_lossy(&self.alignment_fasta())
            )),
            Some(FoundInstead::Failed(failure)) => Prefetched::Failed(format!(
                "has the Seq-align id {}, whose lookup in the database fails ({failure}; {})",
                String::from_utf8_lossy(&self.alignment_fasta()),
                failure.what_ncbi_does(None)
            )),
            Some(FoundInstead::Sequence(other)) => Prefetched::Another {
                other: Some(*other),
                or_itself: false,
                or_none: self.also_held.is_some(),
                reason: reason(),
            },
            Some(FoundInstead::HoldersUnknown {
                instead: Some(other),
                ..
            }) => Prefetched::Another {
                other: Some(*other),
                or_itself: false,
                or_none: true,
                reason: reason(),
            },
            Some(FoundInstead::HoldersUnknown { instead: None, .. }) => {
                Prefetched::ItselfOrNothing { reason: reason() }
            }
            Some(FoundInstead::AlsoHeld(other)) => Prefetched::Another {
                other: Some(*other),
                or_itself: true,
                or_none: false,
                reason: reason(),
            },
            Some(FoundInstead::CopyApart(_, CopyApart::TaxonomyTexts)) => Prefetched::Another {
                other: None,
                or_itself: false,
                or_none: false,
                reason: reason(),
            },
        }
    }

    /// Whether the lookup of this subject's Seq-align id throws on the id where NCBI's prefetch
    /// loads it (`FoundInstead::ThrowsOnId`): the batch's prefetch catches it and writes `Error
    /// pre-fetching sequence data` (`common::database::batch_prefetch`). A lookup that fails on a
    /// file stops NCBI (`FoundInstead::Failed`).
    pub fn prefetch_fails(&self) -> bool {
        matches!(self.found_instead, Some(FoundInstead::ThrowsOnId(_)))
    }

    /// NCBI's message where the lookup of this subject's Seq-align id fails on a file of the
    /// database with NCBI's own `CSeqDBException` (`FoundInstead::Failed`,
    /// `LookupFailure::ncbi_message`; "DB contains no accession info." of a list that puts
    /// version 5 databases with and without the accession table side by side): NCBI's prefetch of
    /// a batch that reports the subject rethrows it, which stops the run with `BLAST Database
    /// error: <message>` (exit status 2) before the batch's reports
    /// (`common::database::batch_prefetch`).
    pub fn prefetch_stops(&self) -> Option<&[u8]> {
        match &self.found_instead {
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
        match &self.found_instead {
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
        self.set.taxonomy_lookups(Some(&self.alignment_id))
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
    /// definition line holds one (a `-subject` record's ids are not looked up there); `None`
    /// where neither throws. `Err` as `tabular_ids`'.
    pub fn lookup_failure(
        &self,
        protein: bool,
    ) -> Result<Option<super::tabular::LookupFailure>, &'static str> {
        let subject_id = replaces_an_id(&self.target().ids);
        let all_ids = !self.record && self.has_replaced_id();
        if !subject_id && !all_ids {
            return Ok(None);
        }
        Ok(self
            .generated_title_word(protein)?
            .is_none()
            .then_some(super::tabular::LookupFailure {
                subject_id,
                all_ids,
                bioseq: false,
            }))
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
    pub fn description_line(&self, protein: bool) -> (Vec<u8>, Defline, bool) {
        let Ok(mut title) = self.bioseq_title(protein, true) else {
            let (start, title) = self.thrown_description();
            return (start, title, true);
        };
        let (start, wrote_gi) = self.description_id();
        if wrote_gi {
            title.text.splice(0..0, *b"  ");
        }
        (start, title, false)
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

    /// The start of the outfmt 0 description line (`description`).
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
    /// generated title and NCBI cannot make it (`UnknownEncoding`); `None` otherwise.
    pub fn heading_skipped_id(&self, protein: bool, every_line: bool) -> Option<Vec<u8>> {
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

    /// The outfmt 5 `Hit_id` and `Hit_def` (`DeflineSet::hit_id_and_def` with the GIs shown;
    /// an empty `Hit_def` is written `No definition line`) and `Hit_accession` (the label of
    /// the `Score`-best id of the subject Bioseq, without its version; a `general` id's tag).
    pub fn xml_hit(&self) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
        let (id, def) = self.set.hit_id_and_def(true);
        let mut accession = Vec::new();
        if let Some(best) = find_best(&self.target().ids, text_score) {
            label::format_label(best, false, &mut accession);
        }
        (id, def, accession)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:365-388
    /// The `HitDescr` of each definition line (XML2 and JSON), without the scientific names.
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
/// in outfmt 0 are `check_pairwise_records`.
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
/// the program.
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

/// LOSAT's refusals of a `-subject` record in outfmt 0, as of a database sequence
/// (`algorithm::common::database::check_shown_parsed_subject`): a protein title that NCBI's
/// `x_CleanAndCompress` reads past the end of, and a description id longer than 65 bytes.
fn check_pairwise_record(
    subject: &ParsedSubject,
    protein: bool,
    program: &str,
) -> anyhow::Result<()> {
    let (start, title) = subject.description(protein);
    let name = || String::from_utf8_lossy(&start).trim_end().to_string();
    if protein && title.reads_past_end {
        anyhow::bail!(
            "the subject record {} has a title of punctuation that NCBI BLAST+'s \
             x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 crashes in outfmt 0), \
             which is not supported by LOSAT's {program}",
            name()
        );
    }
    if start.len() > 65 {
        anyhow::bail!(
            "the subject record {} has an outfmt 0 description id longer than 65 bytes (NCBI \
             BLAST+ 2.17.0 pads the description line without end), which is not supported by \
             LOSAT's {program}",
            name()
        );
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/accguide2.txt:241-266 (the rows of the RefSeq
// kinds that `x_SetFlags` tests; `3+6` is a prefix of three characters and six digits)
/// The RefSeq kind that `CSeq_id::IdentifyAccession` gives `accession` where `x_SetFlags` tests
/// it: a chromosome or ncRNA (`BioseqFacts::refseq_chromosome_or_ncrna`) or an mRNA
/// (`BioseqFacts::refseq_mrna`).
fn refseq_accession_kind(accession: &[u8]) -> Option<RefseqKind> {
    if accession.len() < 3 || accession[2] != b'_' {
        return None;
    }
    let digits = &accession[3..];
    if !digits.iter().all(u8::is_ascii_digit) {
        return None;
    }
    let prefix = [
        accession[0].to_ascii_uppercase(),
        accession[1].to_ascii_uppercase(),
    ];
    match &prefix {
        b"AC" | b"NC" | b"NR" if digits.len() == 6 => Some(RefseqKind::ChromosomeOrNcrna),
        b"NM" | b"XM" if digits.len() == 6 || digits.len() == 9 => Some(RefseqKind::Mrna),
        _ => None,
    }
}

/// A RefSeq kind whose empty title `GenerateDefline` makes with its own function
/// (`refseq_accession_kind`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RefseqKind {
    /// `m_IsNC` or `m_IsNR`.
    ChromosomeOrNcrna,
    /// `m_IsNM`.
    Mrna,
}

impl BioseqFacts {
    /// Sets the flag of the RefSeq kind of `accession` (`refseq_accession_kind`), if any.
    fn add_refseq_accession(&mut self, accession: Option<&[u8]>) {
        match accession.and_then(refseq_accession_kind) {
            Some(RefseqKind::ChromosomeOrNcrna) => self.refseq_chromosome_or_ncrna = true,
            Some(RefseqKind::Mrna) => self.refseq_mrna = true,
            None => {}
        }
    }
}

/// `BioseqFacts` of the Seq-ids `ids` of a Bioseq.
pub(crate) fn bioseq_facts(ids: &[ParsedSeqId]) -> BioseqFacts {
    let mut flags = BioseqFacts::default();
    for id in ids {
        match id {
            ParsedSeqId::Tpg(text) | ParsedSeqId::Tpe(text) | ParsedSeqId::Tpd(text) => {
                flags.third_party = true;
                flags.add_refseq_accession(text.accession.as_deref());
            }
            ParsedSeqId::Other(text)
            | ParsedSeqId::Genbank(text)
            | ParsedSeqId::Embl(text)
            | ParsedSeqId::Ddbj(text) => {
                flags.add_refseq_accession(text.accession.as_deref());
            }
            ParsedSeqId::Pdb(pdb) => {
                let (chain_id, chain) = flags.pdb.get_or_insert_with(|| (Vec::new(), 0));
                if let Some(id) = &pdb.chain_id {
                    *chain_id = id.clone();
                } else if let Some(value) = pdb.chain {
                    *chain = value;
                }
            }
            ParsedSeqId::Patent(patent) => {
                let number = match &patent.cit.id {
                    PatentNumber::Number(number) | PatentNumber::AppNumber(number) => {
                        number.clone()
                    }
                };
                flags.patent = Some((patent.seqid, patent.cit.country.clone(), number));
            }
            _ => {}
        }
    }
    flags
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
    pub fn keeps_scope(outfmt: impl AsRef<[u8]>) -> bool {
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
    /// ([`Self::keeps_scope`]), and whether any XML2/JSON format in the call is secondary
    /// (`!Self::is_xml2`).
    pub fn scope_flags(outfmt: Option<impl AsRef<[u8]>>) -> (bool, bool) {
        match outfmt {
            Some(outfmt) => (
                Self::keeps_scope(outfmt.as_ref()),
                !Self::is_xml2(outfmt.as_ref()),
            ),
            None => (false, false),
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
    /// of no identities) where its lookup would throw.
    Unread,
    /// A report loads the subject, as `Loaded`, in a run whose formats LOSAT refuses after a
    /// subject whose lookup throws (`FormatThrow::Refused`): such a lookup is refused there
    /// (`ParsedSubject::fetched`), for this reason.
    Refused(ThrowRefusal),
    /// A report reaches the subject, whose lookup fails on a file (`lookup_throws`), after the run
    /// loaded another sequence that NCBI's scope gives for it without the lookup
    /// (`ParsedSubject::held_before`: that sequence; `None` where LOSAT cannot list the holders):
    /// nothing throws there and NCBI writes that sequence's record for the subject, which LOSAT
    /// refuses (`ParsedSubject::fetched`). Where the subject throws in one report and is held in
    /// a later one, it is held.
    Held(Option<DbSequence>),
}

/// How one query's report reads its subjects in [`HitLoads::of`]: whether its batch's prefetch is
/// off, so that the report loads each subject where it reaches it (`one_by_one`,
/// `common::database::BatchPrefetch::reads_one_by_one`), and the query batch whose scope it reads
/// in (`batch`: a number that the program gives each query batch of the run, the same in its
/// record, `common::database::TaxdbRun::loads_before`; `None` for a report read in the scope of the
/// report before, as `From<bool>` makes it: one query batch).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReportRead {
    pub one_by_one: bool,
    pub batch: Option<usize>,
}

impl From<bool> for ReportRead {
    fn from(one_by_one: bool) -> Self {
        ReportRead {
            one_by_one,
            batch: None,
        }
    }
}

/// [`HitLoad`] of each subject of a call's reports, by the key that the caller gives its subjects
/// (a rank, or an ordinal id): those whose lookup throws where a report reaches them, with the
/// exception, and those that no report loads; every other subject is loaded. Empty (every subject
/// loaded) where a format of the run refuses a lookup that throws (with the refusal's reason), or
/// no report reads one by one. With the reports that NCBI's run reaches ([`HitLoads::reaches`]).
#[derive(Debug, Default)]
pub struct HitLoads {
    throws: std::collections::HashMap<usize, crate::ncbi_exception::NcbiException>,
    /// The subjects whose lookup fails where NCBI's scope gives another sequence for them
    /// (`HitLoad::Held`).
    held: std::collections::HashMap<usize, Option<DbSequence>>,
    unread: std::collections::HashSet<usize>,
    /// Why the run's formats refuse a lookup that throws (`FormatThrow::Refused`).
    refused: Option<ThrowRefusal>,
    /// The number of the reports given that NCBI's run reaches, where a report's load ended it
    /// (`None`: the run reaches every one).
    reached: Option<usize>,
}

impl HitLoads {
    /// The loads of the query reports `reports`, in the order of the reports: whether the
    /// query's batch has its prefetch off, so that its report loads each subject where it reaches
    /// it (`one_by_one`, `common::database::BatchPrefetch::reads_one_by_one`), and the subject key
    /// of each of its HSP lists (or HSPs) in the report's order (`None` for an entry that shows
    /// nothing). The entries make `Hit`s as the XML writers make them (`blast_xml::hit_ranges`: a
    /// copy of one record that follows another joins its `Hit`, `joins_as_copy`; the entries of
    /// one subject that follow each other are one). In a report that reads one by one, a `Hit`
    /// whose subject's lookup throws (`lookup_throws` of `parsed(key)`) throws (its subjects are
    /// not loaded) at its first row that loads the subject: every row for the run's `rows`
    /// `Every`, the rows of HSPs with identities for `WithIdentities` (`identities(report, entry)`:
    /// the largest of the identities of the HSPs of entry `entry` of report `report`, the
    /// Seq-align's `num_ident`, asked only there); a `Hit` none of whose rows loads it is not loaded
    /// and ends nothing (`HitLoad::Unread`, unless another report loads it or throws there). What
    /// comes after a `Hit` that throws is the run's `after` ([`FormatThrow::of_run`]):
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
    /// A report whose batch prefetches loads every subject (a lookup that throws stops the
    /// prefetch first, `common::database::batch_prefetch_stop`). A subject is loaded where some
    /// report loads it; one that throws in a report throws in every report that reaches it.
    /// `run` (the run's answer, [`FormatThrow::of_formats`]; `None` for `Unread`) that does not
    /// write: a format of the run refuses the throw (`HitLoad::Refused`, with the reason) or none
    /// reads the subjects, and every subject is loaded (`ParsedSubject::fetched` refuses such a
    /// subject where a format reads it).
    ///
    /// NCBI's application ends at the report whose load ends the run (`CATCH_ALL`): no report,
    /// and no batch's prefetch, after it runs ([`HitLoads::reaches`]).
    ///
    /// A subject whose lookup fails on a file is reached where the run has loaded before, in
    /// NCBI's order (the reports before it, the description lines and alignments before it in a
    /// pairwise report, the rows and `Hit`s before it in its report), another sequence that
    /// NCBI gives for it without the lookup ([`ParsedSubject::held_before`]): a holder of its
    /// handle loaded before in the run (the data loader's), or one that NCBI's scope finds for it
    /// loaded before in the report's query batch ([`ReportRead::batch`]; in the run where the
    /// run's scope is kept, [`RunThrow::keeps_scope`]): it is held there (`HitLoad::Held`), nothing
    /// throws, and the report goes on as after a load. The reports come after the loads `before`
    /// of the run (its record of the batches, rounds and calls before them; [`LoadsBefore::NONE`]
    /// for the run's first reports).
    pub fn of<'p, R: Into<ReportRead>, K: Into<ReportKeys>>(
        before: LoadsBefore,
        run: impl Into<FormatThrow>,
        reports: impl IntoIterator<Item = (R, K)>,
        parsed: impl Fn(usize) -> Option<&'p ParsedSubject>,
        identities: impl Fn(usize, usize) -> usize,
    ) -> Self {
        Self::after(
            &std::cell::Cell::new(false),
            &mut SeenLoads::of(before),
            run.into(),
            reports,
            parsed,
            identities,
        )
    }

    /// [`HitLoads::of`] of the reports `reports` that follow those whose loads ended the run where
    /// `ended` holds true (an `EndsRun` subject that threw), which it holds after these reports
    /// too, and after the loads `seen`, to which it adds theirs ([`RunLoads::loads`]).
    fn after<'p, R: Into<ReportRead>, K: Into<ReportKeys>>(
        ended: &std::cell::Cell<bool>,
        seen: &mut SeenLoads,
        run: FormatThrow,
        reports: impl IntoIterator<Item = (R, K)>,
        parsed: impl Fn(usize) -> Option<&'p ParsedSubject>,
        identities: impl Fn(usize, usize) -> usize,
    ) -> Self {
        let mut loads = HitLoads::default();
        let run = match run {
            FormatThrow::Writes(run) => run,
            FormatThrow::Refused(why) => {
                loads.refused = Some(why);
                return loads;
            }
            FormatThrow::Unread => return loads,
        };
        let mut loaded = std::collections::HashSet::new();
        let mut passed = std::collections::HashSet::new();
        // The reports given, of which the run reaches those before the first that comes after
        // its end (`reached`).
        let mut given = 0;
        for (report, (read, keys)) in reports.into_iter().enumerate() {
            if ended.get() {
                loads.reached.get_or_insert(report);
            }
            given = report + 1;
            let read: ReportRead = read.into();
            let one_by_one = read.one_by_one;
            seen.enter(read.batch);
            let subjects = match keys.into() {
                ReportKeys::Pairwise(reads) => {
                    debug_assert!(one_by_one && run.after.is_pairwise());
                    loads.pairwise_report(
                        run,
                        &reads,
                        &parsed,
                        (ended, &mut *seen),
                        (&mut loaded, &mut passed),
                    );
                    continue;
                }
                ReportKeys::Entries(subjects) => subjects,
            };
            if run.after.is_pairwise() {
                // A pairwise report whose batch prefetches loads every subject that it shows,
                // unless the run ended before.
                debug_assert!(!one_by_one);
                let keys = subjects.into_iter().flatten();
                if ended.get() {
                    passed.extend(keys);
                } else {
                    note_loads(seen, &mut loaded, keys, &parsed);
                }
                continue;
            }
            // A table and the organism report look each row's Seq-align id up.
            let hits = if run.after == AfterThrow::WritesRow {
                (0..subjects.len()).map(|entry| entry..entry + 1).collect()
            } else {
                super::blast_xml::hit_ranges(
                    &subjects,
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
                let exception = parsed(first)
                    .filter(|_| one_by_one && !thrown)
                    .and_then(ParsedSubject::lookup_throws);
                let Some(exception) = exception else {
                    if thrown {
                        passed.extend(keys);
                    } else {
                        note_loads(seen, &mut loaded, keys, &parsed);
                    }
                    continue;
                };
                // The `Hit`'s first row that loads its subject (`rows`).
                let loads_it = run.rows == super::tabular::SubjectFetch::Every
                    || range.clone().any(|entry| {
                        subjects[entry].is_some() && run.rows.reads(identities(report, entry))
                    });
                if !loads_it {
                    passed.extend(keys);
                    continue;
                }
                // NCBI gives a sequence loaded before for it: nothing throws.
                if let Some(holder) =
                    parsed(first).and_then(|subject| seen.held(subject, run.keeps_scope))
                {
                    loads.held.insert(first, holder);
                    loaded.extend(keys);
                    continue;
                }
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
        if ended.get() {
            loads.reached.get_or_insert(given);
        }
        loads.unread = passed
            .into_iter()
            .filter(|key| {
                !loaded.contains(key)
                    && !loads.throws.contains_key(key)
                    && !loads.held.contains_key(key)
            })
            .collect();
        loads
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
    /// each in that order after the loads `seen`, to which it adds its own; `ended` records either
    /// end. A subject whose lookup throws (`lookup_throws`)
    /// loads nothing there, nor do the copies in its hit, unless NCBI's scope gives a sequence
    /// loaded before for it at one of its reads (`held_before`: it is held, and an alignment of
    /// it ends nothing); every subject that neither loads is passed.
    fn pairwise_report<'p>(
        &mut self,
        run: RunThrow,
        reads: &PairwiseReads,
        parsed: &impl Fn(usize) -> Option<&'p ParsedSubject>,
        (ended, seen): (&std::cell::Cell<bool>, &mut SeenLoads),
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
            Held(Option<DbSequence>),
        }
        let first = |hit: usize| reads.hits[hit].first().copied().and_then(parsed);
        let state = std::cell::RefCell::new(vec![None::<Read>; reads.hits.len()]);
        let seen = std::cell::RefCell::new(seen);
        // A read of the hit `hit` where the report reaches it, in the report's order: a hit that
        // threw at a read before is looked up again (a holder may have been loaded since).
        let read = |hit: usize| {
            let mut state = state.borrow_mut();
            if let Some(done @ (Read::Loads | Read::Held(_))) = state[hit] {
                return done;
            }
            let mut seen = seen.borrow_mut();
            let now = match first(hit) {
                Some(subject) if subject.lookup_throws().is_some() => seen
                    .held(subject, run.keeps_scope)
                    .map_or(Read::Throws, Read::Held),
                _ => {
                    for &key in &reads.hits[hit] {
                        seen.note(parsed(key));
                    }
                    Read::Loads
                }
            };
            state[hit] = Some(now);
            now
        };
        if !ended.get() {
            // With `-html` the description lines end the run at the first whose load throws
            // (`AfterThrow::description_end`, the writer's answer): nothing after it loads.
            let description_end = run
                .after
                .description_end(reads.described.iter().copied(), |hit| {
                    matches!(read(hit), Read::Throws)
                });
            if description_end.is_some() {
                ended.set(true);
            } else {
                // The description lines that the end above did not read (all, without `-html`; where `-html`
                // read them all, no line threw, and the loop reads them again, which changes nothing:
                // `read` keeps `Loads` and `Held`).
                for &hit in &reads.described {
                    read(hit);
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
                (Some(Read::Held(holder)), Some((&first, _))) => {
                    self.held.insert(first, holder);
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
            held: self
                .held
                .into_iter()
                .filter_map(|(old, holder)| Some((key(old)?, holder)))
                .collect(),
            unread: self.unread.into_iter().filter_map(&key).collect(),
            refused: self.refused,
            reached: self.reached,
        }
    }

    /// What the reports do at the `Hit` of the subject `key`.
    pub fn load(&self, key: usize) -> HitLoad {
        if let Some(&holder) = self.held.get(&key) {
            HitLoad::Held(holder)
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

    /// Whether a report (or a batch's prefetch) loads the subject `key` by a read of its own: not
    /// a subject that NCBI's scope gives another sequence for without the lookup
    /// (`HitLoad::Held`), which reads nothing new and which `ParsedSubject::fetched` refuses with
    /// the words of that sequence (the reads of the taxonomy texts, which BLASTN and TBLASTX
    /// make before they ask `fetched`, would refuse it with the words of its failing lookup).
    pub fn loaded(&self, key: usize) -> bool {
        matches!(self.load(key), HitLoad::Loaded | HitLoad::Refused(_))
    }

    /// Whether NCBI's run reaches the report `report` of those given (counted from 0): false for
    /// every report after the one whose load ended the run, and for every report where an earlier
    /// call of the run ended it ([`RunLoads::loads`]). A batch whose first report the run does not
    /// reach is not prefetched either (NCBI's application ends before it, `CATCH_ALL`).
    pub fn reaches(&self, report: usize) -> bool {
        self.reached.is_none_or(|reached| report < reached)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp (`SeqidToOid`,
    // cited at `common::database::HeaderReader::found_instead`): the lookup of a Seq-align id reads
    // the headers of its candidates (the version check) before it throws on a volume's file.
    /// Whether a report's load reaches the subject `key`, which makes the lookup of its Seq-align id
    /// (`HitLoad::Loaded` or `HitLoad::Throws`, not `HitLoad::Unread`): the lookup that throws has
    /// read the headers of the candidates before it, and posts NCBI's warning of a UTF8String
    /// header among them (`ParsedSubject::scope_utf8_read`, `common::database::load_read_texts`),
    /// though nothing of the subject is loaded (`loaded`: no taxonomy text, no definition line).
    /// A subject that NCBI's scope gives a sequence loaded before for (`HitLoad::Held`, where no
    /// lookup is made) is reached too: `ParsedSubject::fetched` refuses it, so no read is written.
    pub fn reached(&self, key: usize) -> bool {
        self.load(key) != HitLoad::Unread
    }

    /// The exception that the `Hit` of the subject `key` throws (`HitLoad::Throws`).
    pub fn thrown(&self, key: usize) -> Option<&crate::ncbi_exception::NcbiException> {
        if self.held.contains_key(&key) {
            return None;
        }
        self.throws.get(&key)
    }

    /// Whether the load of some subject throws (`HitLoad::Throws`).
    pub fn any_throws(&self) -> bool {
        !self.throws.is_empty()
    }
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
    /// The loads of the run's reports so far, and of its record before them
    /// ([`HitLoads::of`]).
    seen: std::cell::RefCell<SeenLoads>,
}

impl RunLoads {
    /// The record of a run whose reports answer a throw with `answer` (the run's answer,
    /// [`FormatThrow::of_formats`]), before any report.
    pub fn new(answer: impl Into<FormatThrow>) -> Self {
        RunLoads {
            answer: answer.into(),
            ended: std::cell::Cell::new(false),
            seen: std::cell::RefCell::new(SeenLoads::of(LoadsBefore::NONE)),
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
    /// end) and the loads `before` (the run's record before them, [`HitLoads::of`]), which
    /// records where these reports end the run and what they load.
    pub fn loads<'p, R: Into<ReportRead>, K: Into<ReportKeys>>(
        &self,
        before: LoadsBefore,
        reports: impl IntoIterator<Item = (R, K)>,
        parsed: impl Fn(usize) -> Option<&'p ParsedSubject>,
        identities: impl Fn(usize, usize) -> usize,
    ) -> HitLoads {
        let mut seen = self.seen.borrow_mut();
        seen.add(before);
        HitLoads::after(
            &self.ended,
            &mut seen,
            self.answer,
            reports,
            parsed,
            identities,
        )
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
