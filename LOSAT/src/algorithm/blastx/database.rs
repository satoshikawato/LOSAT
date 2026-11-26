//! The database path of BLASTX (`-db NAME`): the search of an NCBI BLAST protein database,
//! as NCBI BLAST+ 2.17.0 `blastx -db NAME -num_threads 1` searches it
//! (PD-LOSAT-DATABASE-SEARCH, P1 and P1.1).
//!
//! What differs from the search of a FASTA subject set (`-subject`):
//!
//! - the subjects are the stored bytes of the database, read in blocks where they are searched
//!   (`x_db_stream.rs`, through `subjects::Subjects::Stream`) by ordinal id and checked where
//!   they are first read (`XStreamDb::x_check`, the checks of `BlastxDatabase::x_check_stores`);
//!   the set-level numbers are the index's (total, count, longest) and 10 for the shortest
//!   sequence (D1);
//! - a database that cannot be read fails the first query batch with a valid context, after
//!   the reports of the batches before it (NCBI does not search a batch without one), with the
//!   check's first error in ordinal order;
//! - after each batch's search, the titles of its reported sequences are read in NCBI's
//!   order and the batch's hit lists are renamed to ranks over records of those sequences
//!   (`ReportedSubjects::for_batch`), so the report writers keep indexing by subject;
//! - the reports show the database (`BlastxDatabase::pairwise_database`, `name`), the subject
//!   ids made from the stored titles, and `stitle`.
//!
//! The search, the composition redo and the report writers are those of the `-subject`
//! path. The opening of the database, the reader's errors and the rules for the stored
//! titles are shared with BLASTP, BLASTN and TBLASTN (`algorithm::common::database`).

use std::collections::BTreeMap;
use std::sync::Arc;

use anyhow::Result;

use super::input::FastaRecord;
use super::runtime::BatchResults;
use super::subjects::Subjects;
use crate::algorithm::common::database::{self as shared, PreliminaryReadError};
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::report::blast_xml::HitTaxonomy;
use crate::report::db_subject::{HitLoads, ParsedSubject};
use crate::report::pairwise::PairwiseDatabase;
use crate::report::tabular::SubjectFetch;
use crate::seqdb::{SeqDb, SeqType};

/// The program that LOSAT's refusals name.
/// The program of the texts of BLASTX's database when the run names none: the database's own
/// whole check from the stream (`x_db_stream`) and the tests.
const PROGRAM: &str = "BLASTX";

// The longest subject BLASTX searches in one piece (NCBI's `MAX_DBSEQ_LEN`, blast_gapalign.h:54;
// the one copy of the macro).
pub(super) use crate::core::subject_split::MAX_DBSEQ_LEN;

/// Sequences whose lengths `BlastxDatabase::x_open_stream` reads with one read of the index (64 KiB
/// of offsets). No NCBI counterpart; the size changes nothing but the number of reads.
const X_LENGTHS_PER_READ: u32 = 16_384;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:115-119
/// The letter of each NCBIstdaa byte.
const NCBISTDAA_TO_AMINOACID: &[u8; 28] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";

/// The database of a `-db` search: the opened database and its name as NCBI keeps it
/// (`BlastxArgs::database_name`: the `-db` value without white space at its ends).
pub(crate) struct BlastxDatabase {
    pub(crate) db: SeqDb,
    name: String,
    /// The taxonomy database of the run and the taxids its volumes have looked up (XML2 and
    /// JSON, and the batches' prefetch of outfmt 0, 5, 6 and 7), kept across the query batches
    /// of the run.
    taxdb: shared::TaxdbRun,
    /// Whether the database has parsed Seq-ids (`shared::check_seqids`).
    parsed_seqids: bool,
    /// The program of the run (`args::ProgramName::upper`: BLASTX, or RPSTBLASTN), which the
    /// database's messages name.
    pub(crate) program: &'static str,
}

impl BlastxDatabase {
    /// `-show_gis` for the subjects of the reports (`ParsedSubject::showing_gis`).
    pub(crate) fn showing_gis(mut self, show_gis: bool) -> Self {
        self.db.set_reports_show_gis(show_gis);
        self
    }

    /// The database of a call of the run after other calls whose loads its record does not hold
    /// where `unrecorded` (`shared::TaxdbRun::note_calls_before`: the call's first query is not the
    /// run's first, `first_query_ordinal` above 0).
    pub(crate) fn after_calls(mut self, unrecorded: bool) -> Self {
        self.taxdb.note_calls_before(unrecorded);
        self
    }

    /// Whether NCBI's formatter keeps its scope across query batches (`shared::TaxdbRun::
    /// set_keeps_scope`: BLAST XML, outfmt 5), and whether XML2/JSON taxonomy reads belong to a
    /// secondary format (`shared::TaxdbRun::set_separate_xml2`).
    pub(crate) fn keeping_scope(mut self, keeps: bool, separate_xml2: bool) -> Self {
        self.taxdb.set_keeps_scope(keeps);
        self.taxdb.set_separate_xml2(separate_xml2);
        self
    }

    /// Mirrors `CBlastFormat::ResetScopeHistory` at the end of a query batch
    /// (`shared::TaxdbRun::reset_scope`).
    pub(crate) fn reset_scope(&self) {
        self.taxdb.reset_scope();
    }

    /// Opens the database of `-db name` (`algorithm::common::database::open_database`)
    /// with the command line's errors.
    #[cfg(test)]
    pub(crate) fn open(name: &str, settings: &ApplicationSettings) -> Result<Self> {
        Self::open_restricted(name, settings, None, None, &mut std::io::sink(), PROGRAM)
    }

    /// `open` with the database restricted to the taxid or seqid list `list` of the command
    /// line, and with the id of the database mask `mask` checked (`-db_soft_mask`,
    /// `-db_hard_mask`; series NP, M4b) (`shared::restrict_database`, whose warnings go to
    /// `diagnostics`).
    pub(crate) fn open_restricted(
        name: &str,
        settings: &ApplicationSettings,
        list: Option<&crate::seqdb::DbList>,
        mask: Option<&crate::blastinput::db_mask::DbMaskOption>,
        diagnostics: &mut dyn std::io::Write,
        program: &'static str,
    ) -> Result<Self> {
        let mut db = shared::open_database(name, SeqType::Protein, settings, program)?;
        // A database with parsed Seq-ids is read (line DBID, ID2a).
        let parsed_seqids =
            shared::restrict_database(&mut db, name, list, mask, diagnostics, settings, program)?;
        Ok(Self {
            db,
            name: name.to_string(),
            taxdb: shared::TaxdbRun::default(),
            parsed_seqids,
            program,
        })
    }

    /// Whether a list kept no sequence of the database (`SeqDb::filtered_empty`).
    pub(crate) fn filtered_empty(&self) -> bool {
        self.db.filtered_empty()
    }

    /// The numbers the statistics read, with NCBI's types (`SeqDb::statistics`); the BLAST
    /// archive's database length reads them too.
    pub(crate) fn statistics(&self) -> crate::seqdb::DbStatistics {
        self.db.statistics()
    }

    /// The taxids of the search's list (`SeqDb::list_taxids`).
    pub(crate) fn list_taxids(&self) -> &[i32] {
        self.db.list_taxids()
    }

    /// The GIs of the search's list (`SeqDb::list_gis`).
    pub(crate) fn list_gis(&self) -> &[i64] {
        self.db.list_gis()
    }

    /// Whether the database has parsed Seq-ids: a batch can then refuse one of its subjects
    /// (`ReportedSubjects::for_batch`).
    pub(crate) fn parsed_seqids(&self) -> bool {
        self.parsed_seqids
    }

    /// Whether a batch's prefetch can stop at a lookup that fails on a file with NCBI's exception
    /// (`shared::TaxdbRun::prefetch_can_stop`, decided once for the run): a database with parsed
    /// Seq-ids whose lookups can fail so.
    pub(crate) fn prefetch_can_stop(&self) -> bool {
        self.parsed_seqids && self.taxdb.prefetch_can_stop(&self.db)
    }

    /// The `-db` value as NCBI keeps it.
    pub(crate) fn name(&self) -> &str {
        &self.name
    }

    /// Whether `taxdb.bti` resolves on NCBI's search path, so that the XML2 and JSON formats
    /// write the scientific name of a taxid (`SeqDb::taxdb_found`) and the batches' prefetch
    /// reads the taxonomy database (`ReportedSubjects::batch_prefetch`); asked once for the
    /// run (`shared::TaxdbRun::found`).
    pub(crate) fn taxdb_found(&self) -> bool {
        self.taxdb.found(&self.db)
    }

    /// NCBI's warning of a failed prefetch, which names the run's program
    /// (`shared::prefetch_warning`).
    pub(crate) fn prefetch_warning(&self) -> Vec<u8> {
        shared::prefetch_warning(self.program)
    }

    /// Whether a later query batch of an XML, XML2 or JSON run could still refuse a taxonomy
    /// database record (`algorithm::common::database::taxdb_may_refuse`: every sequence's taxid
    /// looked up once, before the first byte).
    pub(crate) fn taxdb_may_refuse(&self) -> bool {
        shared::taxdb_may_refuse(&self.db)
    }

    /// Whether the batches' fetches of the run may be refused for NCBI's taxonomy cache
    /// (`shared::TaxdbRun::cache_may_refuse`): the prefetch of outfmt 5 and of 13 to 16, or the
    /// reports' reads where it is off, where NCBI's order of the loads is not known, which an XML
    /// run holds its document for.
    pub(crate) fn cache_may_refuse(&self) -> bool {
        self.taxdb.cache_may_refuse(&self.db)
    }

    /// Whether a later query batch of a run that writes BLAST XML or an XML2 or JSON document
    /// could still be refused, which only the batch's hits show, so that the run holds what it
    /// writes until its epilog (the command line's document, `native::run_inner`; every format of
    /// a library call, `local::report_batches`) and the refusal follows no byte:
    /// - with `taxdb.bti` on the search path, a taxonomy database record (`taxdb_may_refuse`),
    ///   or a fetch where NCBI's taxonomy cache would repeat the offset error
    ///   (`cache_may_refuse`: the prefetch of outfmt 5 and of 13 to 16, or the reports' reads
    ///   where it is off);
    /// - with parsed Seq-ids, one of the batch's subjects (`ReportedSubjects::for_batch`).
    ///
    /// The taxids of every sequence are looked up here, once (DECISIONS 32); the usual run, in
    /// which none would be refused, writes batch by batch.
    pub(crate) fn later_batch_may_refuse(&self) -> bool {
        let taxdb_may_refuse =
            self.taxdb_found() && (self.taxdb_may_refuse() || self.cache_may_refuse());
        taxdb_may_refuse || self.parsed_seqids()
    }

    /// What the `Database:` lines of outfmt 0 show (title, date, real totals).
    pub(crate) fn pairwise_database(&self) -> PairwiseDatabase {
        shared::pairwise_database(&self.db)
    }

    /// `pairwise_database` of a search masked with `mask` (`None`: not masked): each name's
    /// algorithm, which the `Masked using:` lines show (`shared::pairwise_database_masked`).
    pub(crate) fn pairwise_database_masked(
        &self,
        mask: Option<&crate::blastinput::db_mask::ResolvedMask>,
    ) -> PairwiseDatabase {
        shared::pairwise_database_masked(&self.db, mask.map(|mask| mask.shown.as_slice()))
    }

    /// The formatter's step of the database masking `option` (`shared::resolve_database_masking`;
    /// series NP, M4b): its warning goes to `diagnostics` (as one message of the archive
    /// `archive`), and the program options to `archive` and `strategy`.
    pub(crate) fn resolve_mask_option(
        &self,
        option: &crate::blastinput::db_mask::DbMaskOption,
        archive: Option<&crate::report::archive::ArchiveHandle>,
        strategy: Option<&crate::report::archive::strategy::StrategyExport>,
        diagnostics: &mut dyn std::io::Write,
    ) -> Result<Option<crate::blastinput::db_mask::ResolvedMask>> {
        shared::resolve_database_masking(
            &self.db,
            &self.name,
            option,
            self.program,
            archive,
            strategy,
            &mut |text| match archive {
                Some(handle) => handle.write_warning(diagnostics, text),
                None => diagnostics.write_all(text),
            },
        )
    }

    /// LOSAT's refusals of the search masked with `mask` (`ResolvedMask::check_losat_limits`), for
    /// a search with the options `options` (its word and window).
    pub(crate) fn check_mask_limits(
        &self,
        mask: &crate::blastinput::db_mask::ResolvedMask,
        options: &super::args::ResolvedOptions,
    ) -> Result<()> {
        mask.check_losat_limits(
            &self.db,
            self.program,
            // LOSAT's translation gives no residue 0 (a codon with a base it does not know is
            // X, `utils::genetic_code::GeneticCode::get`), so the lookup table holds no word
            // with it.
            crate::blastinput::db_mask::HardScan::Protein {
                word: options.word_size,
                window: options.window_size,
                zero_words: false,
            },
        )
    }

    // A database without sequences (makeblastdb builds none; a volume whose index holds no
    // sequence): NCBI opens it, writes the outfmt 0 prolog, reads the first query batch and
    // fails in that batch's preliminary search, where its reader cannot map the sequence
    // file (BLASTP's `-db`, fixture `LOSAT/tests/seqdb/db/empty`).
    // No NCBI counterpart for the lengths (NCBI takes the longest from the volumes' index headers,
    // seqdbimpl.cpp:128-136, cited at `SeqDb::max_length`): the choices that read the longest
    // (`x_subject_flow::offsets_stay_below_2_30`, `search_core`'s limit, `x_bx_mem`) take the
    // value the check counted, also where a damaged header says less (with a `-window_size` near
    // 2^28 that header chose the kept tables of 2 GiB a thread).
    /// What the search needs of the database before the first query batch, without reading its
    /// sequences: the residues of the longest sequence, from the index's offsets
    /// (`SeqDb::x_protein_lengths`, the lengths of the check's reads; 4 bytes a sequence where the
    /// check read the residues). The empty database's error, and the check's (`x_check_stores`)
    /// where the longest is above `MAX_DBSEQ_LEN` or the index cannot be read, belong to the first
    /// query batch with a valid context, which the caller fails with it. The search checks each
    /// sequence where it first reads it (`XStreamDb::x_check`).
    pub(crate) fn x_open_stream(&self) -> Result<usize> {
        if self.db.first_oid_read_fails() {
            return Err(PreliminaryReadError::from_static("Database memory map file error").into());
        }
        let count = self.db.num_seqs();
        let mut longest = 0;
        let mut first = 0;
        while first < count {
            let end = count.min(first.saturating_add(X_LENGTHS_PER_READ));
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1425-1432
            // Proteins whose NCBI 32-bit length is <= 0 (`0` or `>= 2^31`) are skipped by the
            // search and do not contribute to the searchable longest sequence.
            let lengths = match self.db.x_protein_lengths_i32(first, end) {
                Ok(lengths) if !lengths.is_empty() => lengths,
                _ => return self.x_check_stores(),
            };
            longest = longest.max(
                lengths
                    .iter()
                    .copied()
                    .map(|len| len.max(0) as usize)
                    .max()
                    .unwrap_or(0),
            );
            first += lengths.len() as u32;
        }
        if longest > MAX_DBSEQ_LEN {
            return self.x_check_stores();
        }
        Ok(longest)
    }

    /// The whole check of the database (`x_check_whole`).
    pub(crate) fn x_check_stores(&self) -> Result<usize> {
        x_check_whole_named(&self.db, self.program)
    }

    /// The opened database.
    pub(crate) fn x_seqdb(&self) -> &SeqDb {
        &self.db
    }
}

/// Reads every sequence of the database once, in blocks (`x_check_protein_stores`), and keeps
/// nothing but the residues of the longest sequence. The error (a reader error of NCBI's
/// preliminary search, LOSAT's refusal of a sequence, or the empty database's) is that of the
/// first sequence in ordinal order: the error of a database whose longest is above
/// `MAX_DBSEQ_LEN` (`x_open_stream`), and the error the run reports for a read or a check
/// that fails during the search (`native::search_error`), whichever sequence the search's
/// threads read first.
pub(super) fn x_check_whole(db: &SeqDb) -> Result<usize> {
    x_check_whole_named(db, PROGRAM)
}

/// `x_check_whole` of a run of the program `program` (`BlastxDatabase::program`).
fn x_check_whole_named(db: &SeqDb, program: &str) -> Result<usize> {
    if db.first_oid_read_fails() {
        return Err(PreliminaryReadError::from_static("Database memory map file error").into());
    }
    let mut longest = 0;
    let mut first_too_long = None;
    shared::x_check_protein_stores(db, program, |ordinal, length| {
        longest = longest.max(length);
        if length > MAX_DBSEQ_LEN && first_too_long.is_none() {
            first_too_long = Some((ordinal, length));
        }
    })?;
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:53-54
    // NCBI searches a longer subject in overlapping parts, which BLASTX does not port (the
    // same limit stops `-subject` in `search::search_core`); the refusal belongs to the
    // first batch with a valid context, as the reader's errors.
    if let Some((ordinal, length)) = first_too_long {
        return Err(too_long(
            format_args!("database sequence {}", db.oid(ordinal)),
            length,
        ));
    }
    Ok(longest)
}

/// LOSAT's refusal of a subject longer than `MAX_DBSEQ_LEN`, not an error of NCBI's engine: it
/// gets no `BLAST engine error:` (`native::search_error_as_is`). The one refusal of a database
/// sequence and of a `-subject` record.
#[derive(Debug)]
pub(crate) struct LongSubjectRefusal(String);

impl std::fmt::Display for LongSubjectRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for LongSubjectRefusal {}

/// LOSAT's refusal of the subject `what` (`database sequence 0`, `subject sequence 1`) of
/// `length` residues, longer than `MAX_DBSEQ_LEN`.
pub(super) fn too_long(what: impl std::fmt::Display, length: usize) -> anyhow::Error {
    LongSubjectRefusal(format!(
        "{what} has {length} residues; a subject longer than {MAX_DBSEQ_LEN} residues, which NCBI BLAST+ searches in overlapping parts, is not supported by LOSAT's BLASTX"
    ))
    .into()
}

/// The reported sequences of one query batch, by rank (their order of ordinal ids): the
/// records the report writers index (`gnl|BL_ORD_ID|n`, the stored title, the letters of
/// the stored bytes), the stored titles, and the tabular ids when a tabular format is
/// written.
#[derive(Default)]
pub(crate) struct ReportedSubjects {
    records: Vec<FastaRecord>,
    titles: Vec<Arc<[u8]>>,
    /// `None` where the lookup fails (`common::tabular::SubjectId`).
    tabular_ids: Vec<Option<String>>,
    /// The ordinal id of each rank, ascending.
    ordinals: Vec<u32>,
    /// The OID of each rank, which names a sequence without parsed Seq-ids (the Seq-ids of the
    /// Seq-annots, `general { db "BL_ORD_ID", tag id <oid> }`).
    oids: Vec<crate::seqdb::Oid>,
    /// By rank: the Blast-def-line taxid and what the taxonomy database says of it, and the
    /// taxonomy database's messages of the batch (XML2 and JSON with `-db` only; empty
    /// otherwise, which reads as no taxid).
    taxids: Vec<i32>,
    taxonomy: Vec<HitTaxonomy>,
    /// What the taxonomy database says of each definition line's taxid of a sequence with
    /// parsed Seq-ids, by rank (`shared::read_reported_taxonomy`; empty as `taxonomy`).
    parsed_taxonomy: Vec<Vec<HitTaxonomy>>,
    taxdb_notes: shared::TaxdbNotes,
    /// The reported ranks whose taxonomy the XML2 and JSON reports read (`read_taxonomy`), those
    /// that the reports load; `None` for another format.
    xml_ranks: Option<std::collections::BTreeSet<usize>>,
    /// What the batch's reports do at the `Hit` of each subject, by rank (`HitLoads`, the one
    /// answer: a subject whose `Hit` throws, one that no report loads), which the refusals, the
    /// names, the taxonomy reads and the run's record of the loads take.
    loads: HitLoads,
    /// The sequences with parsed Seq-ids by rank (empty when the batch reports none): their
    /// ids and titles come from their definition lines (`report::db_subject`), and `titles`
    /// holds the title of their first definition line.
    parsed: Vec<Option<Arc<ParsedSubject>>>,
    /// The sources from which NCBI makes the empty stored titles of the others, by rank
    /// (`shared::StoredTitle::source`; empty when none has one).
    sources: Vec<crate::report::defline::SourceTaxname>,
    /// `sacc` and `saccver` of the sequences with parsed Seq-ids and of those with the general id
    /// (`shared::TabularSubjectId::General`) by rank, with `tabular` (empty without them; the
    /// others' are their `sseqid`).
    tabular_accessions: Vec<Option<[String; 2]>>,
    /// What a custom tabular field list reads beyond the default report's ids
    /// (`read_tabular`; empty and `None` otherwise).
    tabular: ReportedTabular,
    /// The taxonomy database's messages of each query's own fetches, by query ordinal, where
    /// `PRE_FETCH_SEQS_LIMIT` turned the batch's prefetch off (`fetch_by_query`): the writers
    /// write them after the query's warnings (`ShownSubjects::query_fetches`).
    query_fetches: std::collections::BTreeMap<usize, Vec<u8>>,
    /// The texts of the subjects that each query's report reads one by one (`fetch_by_query`),
    /// by query ordinal, where the format's writer places them at the reads
    /// (`report::query_warnings::ReadCursor`; `ShownSubjects::query_reads`).
    query_reads: std::collections::BTreeMap<usize, crate::report::query_warnings::QueryReads>,
    /// `-db_soft_mask` and `-db_hard_mask` with a pairwise report: the stored masks of the
    /// search's algorithm of each rank (`read_db_masks`); `None` otherwise.
    db_masks: Option<Vec<Vec<(u32, u32)>>>,
}

/// The values of a custom tabular report's fields that read more of the reported sequences than
/// the default report (`ReportedSubjects::read_tabular`), by rank.
#[derive(Default)]
struct ReportedTabular {
    /// The ids of each definition line of a sequence with parsed Seq-ids (`sallseqid`, `sallgi`,
    /// `sallacc`); empty for the others, whose `sall*` are their own ids. Where
    /// `x_SetSubjectIds` throws, the lines before (`AllTabularIds::failed`).
    all_ids: Vec<crate::report::db_subject::AllTabularIds>,
    /// `staxid` and `staxids` of the reported sequences.
    taxonomy: Vec<crate::report::tabular::TabularTaxonomy>,
    /// The taxonomy database's names that the rows look up (`None`: every lookup fails).
    names: Option<crate::report::tabular::TabularNames>,
}

/// The titles of a query batch's reported sequences (`ReportedSubjects::read_titles`), by ordinal
/// id: `None` for a listed sequence that no report shows; the sequences with parsed Seq-ids; the
/// sources from which NCBI makes the others' empty stored titles (only those not `Unnamed`).
pub(crate) struct BatchTitles {
    listed: BTreeMap<usize, Option<Arc<[u8]>>>,
    parsed: BTreeMap<usize, Arc<ParsedSubject>>,
    sources: BTreeMap<usize, crate::report::defline::SourceTaxname>,
    /// The first read of a record with a UTF8String where a VisibleString is expected
    /// (`shared::HeaderReader::first_utf8`): NCBI posts its warning there (`utf8_post`).
    pub(crate) utf8: Option<shared::Utf8Read>,
}

impl BatchTitles {
    /// The reported sequences with parsed Seq-ids, whose lookups decide whether the batch's
    /// prefetch fails (`shared::BatchPrefetch::prefetch_fails_among`, asked by
    /// `native::BatchPrefetch::report_read`).
    pub(crate) fn parsed_subjects(&self) -> impl Iterator<Item = &ParsedSubject> {
        self.parsed.values().map(Arc::as_ref)
    }
}

/// What one query's pairwise report (outfmt 0 to 4) of the hit lists `lists` reads, by the key
/// that `lists` names the subjects by, an ordinal id or a rank (`QueryGroup::read_order` of the
/// rows' hits as the writer reads them, `reading`: each list's HSPs that the reports show in the
/// reports' order, `HspList::evalue_order`, with the values of their rows, `report::hsp_scores`
/// of the query of `query_length` letters; `joins`: a copy joins the hit before).
fn pairwise_read_order(
    lists: &[super::results::HspList],
    query_length: usize,
    gapped: bool,
    joins: impl Fn(u32, u32) -> bool,
    reading: crate::report::pairwise::PairwiseReading,
) -> Result<crate::report::db_subject::PairwiseReads> {
    let mut hits = Vec::new();
    for list in lists {
        for hsp in list.evalue_order() {
            if super::report::has_seqalign(hsp, gapped) {
                hits.push((
                    list.oid as u32,
                    super::report::hsp_scores(hsp, query_length)?,
                ));
            }
        }
    }
    Ok(crate::report::x_pairwise_par::QueryGroup::read_order(
        hits,
        |scores| *scores,
        joins,
        reading,
        false,
        None,
    ))
}

impl ReportedSubjects {
    /// The hit lists of the query batch `batches` (subjects named by ordinal) without the
    /// subjects that NCBI makes no Seq-align of (`shared::SeqAligned`), where its traceback
    /// makes the batch's Seq-aligns: after the cut to the hit list size (`Runtime::batch_end`),
    /// before `read_titles`.
    pub(crate) fn drop_unaligned(
        database: &BlastxDatabase,
        batches: &mut [BatchResults],
        gapped: bool,
    ) -> Option<crate::seqdb::SeqLocation> {
        // A gapped search fetched every subject of the redo or the traceback with
        // `check_oid_exclusion` (`Runtime::batch_end`); an ungapped one fetches no subject.
        let mut aligned = shared::SeqAligned::of(&database.db, gapped)?;
        for batch in batches {
            for lists in &mut batch.queries {
                aligned.retain(lists, |list| u32::try_from(list.oid).ok());
            }
        }
        aligned.first_utf8()
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:264-275
    /// The titles of the reported sequences of the batches `batches` (one query batch), read
    /// query by query in hit-list order where NCBI's traceback makes the batch's Seq-aligns
    /// (their Seq-ids), before the batch's prefetch: NCBI's reader error stops the run at this
    /// batch (`ErrorAfterBatchReports`). A sequence with parsed Seq-ids is read with every Seq-id
    /// (`shared::read_header`). A title LOSAT does not show as NCBI is its refusal
    /// (`check_reported_title`, after the batch's reads, whose error comes first: NCBI's Seq-aligns
    /// end there, before any report; with `pairwise`, `StoredTitle::check_shown_source`). The rest
    /// of the batch's reported sequences is `for_batch`. `one_by_one`: the reports may read the
    /// reported sequences one by one (`shared::HeaderReader::new`).
    pub(crate) fn read_titles(
        database: &BlastxDatabase,
        batches: &[BatchResults],
        gapped: bool,
        (pairwise, title_reads): (bool, shared::TitleReads),
        (earlier, one_by_one): (shared::EarlierUtf8, bool),
    ) -> Result<BatchTitles> {
        // Every listed subject gets a record (the tabular writer looks its list's subject
        // up before it filters the HSPs); only the reported ones have their titles read.
        let mut listed: BTreeMap<usize, Option<Arc<[u8]>>> = BTreeMap::new();
        let mut parsed: BTreeMap<usize, Arc<ParsedSubject>> = BTreeMap::new();
        // The sources of the empty stored titles that NCBI makes from them.
        let mut sources: BTreeMap<usize, crate::report::defline::SourceTaxname> = BTreeMap::new();
        let headers = shared::HeaderReader::new(&database.db, one_by_one);
        headers.note_earlier(earlier);
        // The first title or shown-source refusal that the reports refuse (`check_reported_title`,
        // `StoredTitle::check_shown_source`), held until the batch's header reads finish so an
        // earlier refusal does not get overtaken by a later one or overtake a later header error.
        let mut refused = None;
        for batch in batches {
            for (qidx, lists) in batch.queries.iter().enumerate() {
                for list in lists {
                    let reported = super::report::reported(list, gapped);
                    if matches!(listed.get(&list.oid), Some(Some(_)))
                        || (!reported && listed.contains_key(&list.oid))
                    {
                        continue;
                    }
                    if !reported {
                        listed.insert(list.oid, None);
                        continue;
                    }
                    let ordinal = u32::try_from(list.oid)?;
                    match shared::read_header(
                        &headers,
                        ordinal,
                        batch.query_ordinal + qidx,
                        database.program,
                    )? {
                        shared::ReportedHeader::Stored(stored) => {
                            if let Err(refusal) = shared::check_reported_title(
                                &database.db,
                                ordinal,
                                stored.replaced_bytes,
                                title_reads,
                                database.program,
                            ) {
                                refused.get_or_insert(refusal);
                            } else if pairwise {
                                if let Err(refusal) = stored.check_shown_source(
                                    || database.db.oid(ordinal),
                                    database.program,
                                ) {
                                    refused.get_or_insert(refusal);
                                }
                            }
                            if stored.source != crate::report::defline::SourceTaxname::Unnamed {
                                sources.insert(list.oid, stored.source);
                            }
                            listed.insert(list.oid, Some(Arc::from(stored.title)));
                        }
                        shared::ReportedHeader::Parsed(subject) => {
                            if let Err(refusal) = shared::check_reported_title(
                                &database.db,
                                ordinal,
                                subject.replaced_title_bytes(),
                                title_reads,
                                database.program,
                            ) {
                                refused.get_or_insert(refusal);
                            }
                            listed.insert(list.oid, Some(Arc::from(subject.first_title())));
                            parsed.insert(list.oid, subject);
                        }
                    }
                }
            }
        }
        if let Some(refusal) = refused {
            return Err(refusal);
        }
        Ok(BatchTitles {
            listed,
            parsed,
            sources,
            utf8: headers.first_utf8(),
        })
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:286-291
    /// The reported sequences of the batches `batches` (one query batch: the results of
    /// one search), with their hit lists renamed from ordinal ids to ranks, from their titles
    /// (`titles`, `read_titles`), before any report of the batch: a title LOSAT does not show as
    /// NCBI is its refusal (with `tabular`, `tabular_subject_id`; with `pairwise`,
    /// `check_shown_protein_title`). A sequence with parsed Seq-ids has LOSAT's refusals of `ParsedSubject::fetched`, of
    /// `shared::check_shown_parsed_subject` (with `pairwise`), of `shared::refuse_joined_hits`
    /// (with `joined`: outfmt 0, 5 and 13 to 16) and of `ParsedSubject::tabular_ids` (with
    /// `tabular_reads`, and `ParsedSubject::check_printed_word` of the ids it prints, `words`).
    /// `tabular`: a tabular report, which takes every reported sequence's ids
    /// (`tabular_ids`); `tabular_reads`: one whose fields read the sequence found by the
    /// subject's Seq-align id (`report::tabular::reads_subject`), the only kind that refuses the
    /// ids. `structured`: outfmt 5 or 13 to 16, which read the Bioseq. `taxonomy`: whether the
    /// report writes the sequences' taxids or their names (13 to 16, a table of taxonomy fields;
    /// `ParsedSubject::fetched`). `prefetched`: whether the batch's prefetch reads the subjects
    /// (`shared::prefetches_subjects`, and the batch's Seq-aligns do not outnumber a positive
    /// `PRE_FETCH_SEQS_LIMIT`, `shared::prefetch_off`), which where no report reads their Bioseqs
    /// refuses what it cannot read as NCBI does (`shared::check_prefetched_parsed`).
    /// The batch's prefetch stops the run where it reads a subject whose lookup fails on a file
    /// with NCBI's exception (`prefetched`, `shared::batch_prefetch`).
    /// `read`: how the batch's reports read the subjects (`native::BatchPrefetch::report_read`, the
    /// answer of `shared::BatchPrefetch`): one by one where its prefetch is off, or fails (at an
    /// unranged row before it loads anything, or at a subject whose lookup throws on the id,
    /// `ReportRead::prefetch_fails`: after the subjects that it may have loaded,
    /// `shared::BatchPrefetch::loads_one_by_one`), in the scope of the batch that it numbers (its
    /// first query in the run). `cut`: the formatter's cut of each query's
    /// lists to the hit list size where RPSTBLASTN merges several volumes
    /// (`ResolvedOptions::shown_lists`; `usize::MAX` for every other search), which counts
    /// following lists whose Seq-align ids match as one subject (`BatchResults::shown`, which
    /// every report that NCBI cuts reads, `ResolvedOptions::shown_lists`: BLAST XML, XML2 and
    /// JSON too, whose `Hit`s `HitLoads` reads from the same lists, `BatchResults::cut_lists`).
    /// `after`: what the reports do after a subject whose lookup throws (the run's answer,
    /// `FormatThrow::of_formats`, or a refusal with its reason), for `HitLoads::of` (where `read`
    /// is one by one, which subjects throw and which no report loads; the refusals, the names, the
    /// tables' rows, the taxonomy reads and the run's record of the loads, `note_loaded`, take that
    /// one answer), with how a pairwise report reads the subjects (`pairwise_read_order`): the
    /// batch's reports come after the loads of the record (the batches before and the calls
    /// before, `TaxdbRun::loads_before`), in the batch's scope (`HitLoads::of`).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn for_batch(
        database: &BlastxDatabase,
        subjects: Subjects<'_>,
        batches: &mut [BatchResults],
        titles: BatchTitles,
        gapped: bool,
        tabular: bool,
        tabular_reads: bool,
        pairwise: bool,
        structured: bool,
        taxids: bool,
        joined: bool,
        taxonomy: bool,
        words: crate::report::tabular::PrintedIds,
        prefetched: bool,
        read: crate::report::db_subject::ReportRead,
        cut: usize,
        after: (
            crate::report::db_subject::FormatThrow,
            Option<crate::report::pairwise::PairwiseReading>,
        ),
    ) -> Result<Self> {
        let BatchTitles {
            listed,
            parsed,
            mut sources,
            utf8: _,
        } = titles;
        // A batch whose prefetch looks up a Seq-align id whose lookup fails on a file stops the
        // run there (`shared::batch_prefetch`), before any of its reports: none of their checks
        // and reads is made.
        let reports = !(prefetched
            && parsed
                .values()
                .any(|subject| subject.prefetch_stops().is_some()));
        if reports && joined && !parsed.is_empty() {
            for batch in batches.iter() {
                for (qidx, lists) in batch.queries.iter().enumerate() {
                    shared::refuse_joined_hits(
                        batch.query_ordinal + qidx,
                        lists
                            .iter()
                            .filter(|list| super::report::reported(list, gapped))
                            .map(|list| {
                                let ordinal = u32::try_from(list.oid).ok()?;
                                Some((ordinal, parsed.get(&list.oid)?.as_ref()))
                            }),
                        database.db.oid_map(),
                        database.program,
                    )?;
                }
            }
        }
        // RPSTBLASTN's merged lists of several volumes, longer than the hit list size: NCBI's
        // formatter keeps the Seq-aligns of its first `cut` subjects, counting following ones
        // whose Seq-ids match as one (`PruneSeqalign`, cited at `shared::joined_cut`): a copy of
        // one record that follows another is kept with it, as RPSBLAST's lists are
        // (`blastp::blast_engine::database::ListedTitles::joined_cut`); two other profiles that
        // the count reaches with matching Seq-align ids are LOSAT's refusal.
        if reports && cut != usize::MAX {
            for batch in batches.iter_mut() {
                let shown = batch
                    .queries
                    .iter()
                    .enumerate()
                    .map(|(qidx, lists)| {
                        if parsed.is_empty() || lists.len() <= cut {
                            return Ok(cut);
                        }
                        shared::joined_cut(
                            batch.query_ordinal + qidx,
                            cut,
                            lists.iter().map(|list| {
                                let ordinal = u32::try_from(list.oid).ok()?;
                                Some((ordinal, parsed.get(&list.oid)?.as_ref()))
                            }),
                            database.db.oid_map(),
                            database.program,
                        )
                    })
                    .collect::<Result<Vec<usize>>>()?;
                batch.shown = shown;
            }
        }
        // What the reports do at each subject's `Hit`, by ordinal id (renamed to ranks below),
        // with the `Hit`s that the writer makes (`blast_xml::write_blastx_xml`: the reported lists
        // of each query that the cut above keeps, `report::reported`).
        let report_lists: Vec<&[super::results::HspList]> = batches
            .iter()
            .flat_map(|batch| (0..batch.queries.len()).map(move |qidx| batch.cut_lists(qidx, cut)))
            .collect();
        let one_by_one = read.one_by_one;
        // A pairwise report's reads (`pairwise_read_order`), of the lists that the formatter's cut
        // keeps (`BatchResults::cut_lists`).
        let (after, pairwise_reading) = after;
        let pairwise_reading = pairwise_reading
            .filter(|_| one_by_one && after.run().is_some_and(|run| run.after.is_pairwise()));
        let mut pairwise_reads = Vec::new();
        if let Some(reading) = pairwise_reading {
            for batch in batches.iter() {
                for index in 0..batch.queries.len() {
                    pairwise_reads.push(pairwise_read_order(
                        batch.cut_lists(index, cut),
                        batch.prepared.original_lengths[index],
                        gapped,
                        |previous, oid| match (
                            parsed.get(&(previous as usize)),
                            parsed.get(&(oid as usize)),
                        ) {
                            (Some(earlier), Some(later)) => {
                                crate::report::db_subject::joins_as_copy(earlier, later)
                            }
                            _ => false,
                        },
                        reading,
                    )?);
                }
            }
        }
        let mut pairwise_reads = pairwise_reads.into_iter();
        // The reports of one query batch, in its scope.
        let loads = HitLoads::of(
            database.taxdb.loads_before(read.batch),
            after,
            report_lists.iter().map(|lists| {
                let keys = match pairwise_reads.next() {
                    Some(reads) => crate::report::db_subject::ReportKeys::Pairwise(reads),
                    None => crate::report::db_subject::ReportKeys::Entries(
                        lists
                            .iter()
                            .map(|list| super::report::reported(list, gapped).then_some(list.oid))
                            .collect(),
                    ),
                };
                (read, keys)
            }),
            |oid| parsed.get(&oid).map(Arc::as_ref),
            // A list's rows read its subject where one of the HSPs that they show has
            // identities (`report::reads_list`).
            |report, list| {
                report_lists[report][list]
                    .hsps
                    .iter()
                    .filter(|hsp| super::report::has_seqalign(hsp, gapped))
                    .map(|hsp| hsp.identity)
                    .max()
                    .unwrap_or(0)
            },
        );
        let taxdb_found = database.db.taxdb_found();
        let mut reported = ReportedSubjects::default();
        for (&oid, title) in &listed {
            let ordinal = u32::try_from(oid)?;
            let shown = title.is_some();
            let title = title.clone().unwrap_or_else(|| Arc::from(&b""[..]));
            let source = sources.remove(&oid).unwrap_or_default();
            if reports && pairwise && shown && !parsed.contains_key(&oid) {
                shared::check_shown_protein_title(
                    &title,
                    || database.db.oid(ordinal),
                    database.program,
                )?;
            }
            if let Some(subject) = parsed.get(&oid) {
                // The formats that read the subject's Bioseq; the Seq-annots and the archive
                // write its Seq-align id alone.
                let reads = tabular_reads || pairwise || structured || taxids || taxonomy;
                if reports && reads {
                    let sought = shared::weak_sought(&database.db, subject, ordinal, one_by_one);
                    sought
                        .as_ref()
                        .unwrap_or(subject)
                        .fetched(taxonomy, one_by_one, loads.load(oid))
                        .map_err(|reason| {
                            shared::parsed_refusal(
                                database.db.oid(ordinal),
                                &reason,
                                database.program,
                            )
                        })?;
                } else if reports && prefetched {
                    shared::check_prefetched_parsed(
                        &database.db,
                        subject,
                        ordinal,
                        database.program,
                    )?;
                }
                if reports && pairwise {
                    shared::check_pairwise_parsed(
                        subject,
                        loads.load(oid),
                        true,
                        taxdb_found,
                        || database.db.oid(ordinal),
                        database.program,
                    )?;
                }
                if tabular {
                    let ids = subject
                        .tabular_ids(true)
                        .and_then(|ids| subject.check_printed_word(true, words).map(|()| ids));
                    // A table of fields that read nothing of the subject prints none of its ids.
                    // No ids where the lookup throws (`ParsedSubject::tabular_ids`).
                    let found = match ids {
                        Ok(ids) => ids,
                        Err(_) if !tabular_reads || !reports => Some(Default::default()),
                        Err(problem) => {
                            return Err(shared::parsed_refusal(
                                database.db.oid(ordinal),
                                &format!("has {problem}"),
                                database.program,
                            ))
                        }
                    };
                    // The ids of a stored definition line are VisibleStrings: ASCII.
                    let text = |bytes: &[u8]| String::from_utf8_lossy(bytes).into_owned();
                    reported
                        .tabular_ids
                        .push(found.as_ref().map(|[id, _, _]| text(id)));
                    if reported.tabular_accessions.is_empty() {
                        reported.tabular_accessions = vec![None; listed.len()];
                    }
                    reported.tabular_accessions[reported.ordinals.len()] = found
                        .as_ref()
                        .map(|[_, accession, version]| [text(accession), text(version)]);
                }
                if reported.parsed.is_empty() {
                    reported.parsed = vec![None; listed.len()];
                }
                reported.parsed[reported.ordinals.len()] = Some(Arc::clone(subject));
            } else if tabular {
                // An unreported subject's id is never written, nor an id of a table of fields
                // that read nothing of the subject.
                let id = if reports && shown && tabular_reads {
                    let oid = database.db.oid(ordinal);
                    Some(shared::tabular_subject_id(
                        &title,
                        &source,
                        oid,
                        true,
                        database.program,
                    )?)
                } else {
                    None
                };
                // A stored title is a VisibleString (checked above), so its id is ASCII.
                let text = |bytes: &[u8]| String::from_utf8_lossy(bytes).into_owned();
                reported.tabular_ids.push(match &id {
                    Some(shared::TabularSubjectId::Unknown(_)) => None,
                    Some(id) => Some(text(&id.id())),
                    None => Some(String::new()),
                });
                if let Some(id @ shared::TabularSubjectId::General(_)) = &id {
                    if reported.tabular_accessions.is_empty() {
                        reported.tabular_accessions = vec![None; listed.len()];
                    }
                    let accession = text(&id.accession());
                    reported.tabular_accessions[reported.ordinals.len()] =
                        Some([accession.clone(), accession]);
                }
            }
            let encoded = subjects.read_encoded(oid)?;
            // A taxid list's ordinal ids are its included OIDs in order; NCBI prints the OID
            // (line DI-8).
            let db_oid = database.db.oid(ordinal);
            reported.records.push(FastaRecord {
                internal_id: db_oid.fasta_id(),
                title: title.to_vec(),
                sequence: encoded
                    .iter()
                    .map(|&residue| NCBISTDAA_TO_AMINOACID[usize::from(residue)])
                    .collect(),
                lowercase_masks: Vec::new(),
                source: None,
                around: None,
                skipped: None,
            });
            shared::keep_source(
                &mut reported.sources,
                reported.ordinals.len(),
                listed.len(),
                source,
            );
            reported.titles.push(title);
            reported.ordinals.push(ordinal);
            reported.oids.push(db_oid);
        }
        // Ranks keep the order of the ordinal ids, so every tie broken by subject index
        // is broken as before.
        let rank: BTreeMap<usize, usize> = listed
            .keys()
            .enumerate()
            .map(|(rank, &oid)| (oid, rank))
            .collect();
        for batch in batches.iter_mut() {
            for lists in &mut batch.queries {
                for list in lists.iter_mut() {
                    list.oid = rank[&list.oid];
                }
            }
        }
        reported.loads = loads.renamed(|oid| rank.get(&oid).copied());
        if reports && taxids {
            // The taxid of each reported sequence, read apart from the titles for the XML2 and
            // JSON formats only, after the batch's prefetch (`read_taxonomy`); the sequences with
            // a title are the reported ones. Of those, the ones that the reports load
            // (`HitLoads::loaded`): the reads of a subject whose `Hit` throws, and of one shown
            // only after it, are never made (review 1 of line DBID pfoff-1, H1).
            reported.xml_ranks = Some(
                listed
                    .iter()
                    .enumerate()
                    .filter(|(rank, (_, title))| title.is_some() && reported.loads.loaded(*rank))
                    .map(|(rank, _)| rank)
                    .collect(),
            );
        }
        Ok(reported)
    }

    /// NCBI's prefetch of the query batch whose hit lists are `batches` (renamed to ranks by
    /// `for_batch`) on standard error (`shared::batch_prefetch`): the messages about the
    /// taxonomy database of the sequences of the lists that the reports show
    /// (`report::reported`), the ones whose titles `for_batch` reads, with the definition lines
    /// of those with parsed Seq-ids, then the warning where the prefetch fails at a subject whose
    /// lookup fails (the protein subject has no strand). Without `taxdb.bti` on the search path
    /// no text, and nothing read. A taxonomy database that LOSAT cannot read as NCBI does is
    /// refused (`shared::TaxonomyRefusal`); a header whose taxid does not decode is the reader's
    /// error.
    ///
    /// A subject whose lookup fails on a file stops the run there (`shared::batch_prefetch`):
    /// BLASTX writes each batch's reports before it searches the next, so the error of the batch
    /// of `queries` names query 0 (`native::reported_batch`; `native::reported_error` takes the
    /// error alone), and the run's record holds what NCBI's data loader loaded for the batches before
    /// (`note_loaded`).
    pub(crate) fn batch_prefetch(
        &self,
        database: &BlastxDatabase,
        batches: &[BatchResults],
        gapped: bool,
        queries: &shared::PrefetchQueries<'_>,
    ) -> Result<crate::report::diag_log::Posts> {
        shared::batch_prefetch(
            &database.db,
            &database.taxdb,
            std::iter::empty(),
            self.prefetch_subjects(batches, gapped),
            database.program,
            // The model of each batch follows the reports of the batches before, so the texts
            // that their loads left are the model's (`ReadLoad::place`): no outfmt 5 reports
            // before it (`PrefetchQueries::after_xml_reports`).
            queries,
        )
    }

    /// The run's record of the sequences that NCBI's BLAST database data loader loaded for the
    /// query batch whose hit lists are `batches` (`shared::TaxdbRun::note_loaded`): BLASTX writes
    /// each batch's reports before it searches the next and keeps none of their lists, so each
    /// batch that does not stop is recorded for the stop of a later batch's prefetch (nothing for
    /// a database that cannot stop). `loads`: what the batch's loads read of the lists that the
    /// reports show (`shared::BatchLoads::Prefetched` where its prefetch read them, BLASTX's
    /// prefetched set being its reported set, `BatchLoads::reported_rows`; where it is off, or
    /// failed before it loaded anything, the rows of its reports; review 3 of line DBID loadord,
    /// H2), of each query whose report NCBI writes (`BatchLoads::report_spans`, the walk of
    /// `shared::LoadedBefore::subjects`: not one whose outfmt 5 report NCBI's XML stream cut
    /// before its hits, `cuts`, by query of the reader batch; review 6, H2). Not
    /// `LoadedBefore::subjects` itself: these loads are `BatchPrefetch::loads` (the rows where
    /// the prefetch fails at an empty subject row with the prefetch on, which the Seq-aligns'
    /// count that `subjects` asks cannot tell), and the lists nest in `BatchResults` by query
    /// of each search batch, not one slice by query of the call (review 7, L2). `batch`: the
    /// reader batch's number in the run (`for_batch`).
    pub(crate) fn note_loaded(
        &self,
        (database, batch): (&BlastxDatabase, usize),
        batches: &[BatchResults],
        gapped: bool,
        loads: shared::BatchLoads,
        cuts: shared::XmlCuts<'_>,
    ) {
        if database.parsed_seqids() {
            let rows = loads.reported_rows();
            database.taxdb.note_loaded(&database.db, batch, || {
                batches
                    .iter()
                    .flat_map(move |batch| {
                        let first = batch.query_ordinal;
                        loads
                            .report_spans(cuts, first..first + batch.queries.len())
                            .flat_map(move |queries| {
                                &batch.queries[queries.start - first..queries.end - first]
                            })
                    })
                    .flatten()
                    // Not a subject that no report loads (`HitLoads::loaded`: one whose `Hit`
                    // throws, or that is shown only after it).
                    .filter(move |list| {
                        super::report::reads_list(list, gapped, rows) && self.loads.loaded(list.oid)
                    })
                    .map(|list| {
                        shared::PrefetchSubject::unflipped(
                            self.ordinals[list.oid],
                            self.parsed_subject(list.oid),
                        )
                    })
            });
        }
    }

    /// The sequences that the prefetch of the query batch whose hit lists are `batches` reads: of
    /// the lists that the reports show (`report::reported`).
    fn prefetch_subjects<'s>(
        &'s self,
        batches: &'s [BatchResults],
        gapped: bool,
    ) -> impl Iterator<Item = shared::PrefetchSubject<'s>> + 's {
        batches
            .iter()
            .flat_map(|batch| batch.queries.iter().flatten())
            .filter(move |list| super::report::reported(list, gapped))
            .map(|list| {
                shared::PrefetchSubject::unflipped(
                    self.ordinals[list.oid],
                    self.parsed_subject(list.oid),
                )
            })
    }

    /// NCBI's messages about the taxonomy database where the reports of the query batch whose hit
    /// lists are `batches` load one by one (`shared::BatchPrefetch::loads_one_by_one`: its
    /// prefetch is off, `shared::prefetch_off`, or failed, after its posts: the reports read the
    /// subjects that they reach, `HitLoads::reached`):
    /// each query's report fetches the sequences it shows, where the first format's rows read
    /// their Bioseqs (`rows`, `TabularSpec::subject_fetch`, by the engine's identities of each
    /// shown HSP), after the query's own warnings. The messages are kept by query for the
    /// writers (`query_fetches`). Refusals as `batch_prefetch`; a batch that fetches nothing
    /// still asks for the run's refusal of NCBI's taxonomy cache, unless no batch of the run
    /// prefetches (`shared::prefetch_always_off`) and no row reads a subject. Rows that read a
    /// subject only for its identities ask in every batch: BLASTX writes each batch's reports
    /// before the next is searched, so a run whose later rows may look a tax id up refuses at
    /// its first batch, before any byte (also where none of them does, as before).
    /// `places`: the format's writer writes each query's texts where it reads each subject
    /// (`query_reads`; outfmt 0, 5, 13 to 16 and the tabular rows); `row_lookups`: the rows that
    /// look a taxid up (`shared::subject_read_texts`). The subjects are read in the order of the hit
    /// lists, or, for the pairwise report (`pairwise`: how it reads them), those that it shows in
    /// the order that it reads them (`pairwise_read_order`, `PairwiseReads::subjects`).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn fetch_by_query(
        &mut self,
        database: &BlastxDatabase,
        batches: &[BatchResults],
        gapped: bool,
        (rows, format): (SubjectFetch, shared::PrefetchFormat),
        places: bool,
        row_lookups: shared::RowLookups,
        pairwise: Option<crate::report::pairwise::PairwiseReading>,
    ) -> Result<()> {
        if rows == SubjectFetch::Never && shared::prefetch_always_off(format) {
            return Ok(());
        }
        for batch in batches {
            for (index, lists) in batch.queries.iter().enumerate() {
                let ranks: Vec<u32> = match pairwise {
                    Some(reading) => pairwise_read_order(
                        lists,
                        batch.prepared.original_lengths[index],
                        gapped,
                        |previous, rank| {
                            crate::report::db_subject::rank_joins_as_copy(
                                &self.parsed,
                                previous as usize,
                                rank as usize,
                            )
                        },
                        reading,
                    )?
                    .subjects()
                    .into_iter()
                    .map(|rank| rank as u32)
                    .collect(),
                    None => lists
                        .iter()
                        .filter(|list| super::report::reads_list(list, gapped, rows))
                        .map(|list| list.oid as u32)
                        .collect(),
                };
                // Not a subject that no report reaches after one whose load throws
                // (`HitLoads::reached`).
                let subjects = ranks
                    .into_iter()
                    .filter(|&rank| self.loads.reached(rank as usize))
                    .map(|rank| {
                        (
                            rank,
                            self.ordinals[rank as usize],
                            self.parsed_subject(rank as usize),
                        )
                    });
                let reads = shared::subject_read_texts(
                    &database.db,
                    &database.taxdb,
                    subjects,
                    row_lookups,
                    database.program,
                )?;
                if reads.is_empty() {
                    continue;
                }
                let ordinal = batch.query_ordinal + index;
                if places {
                    self.query_reads.insert(ordinal, reads);
                } else {
                    let messages = reads.into_iter().flat_map(|read| read.text).collect();
                    self.query_fetches.insert(ordinal, messages);
                }
            }
        }
        Ok(())
    }

    /// The messages of the query `ordinal`'s own fetches (`fetch_by_query`; empty where the
    /// batch prefetched).
    pub(crate) fn query_fetches(&self, ordinal: usize) -> &[u8] {
        self.query_fetches.get(&ordinal).map_or(&[], Vec::as_slice)
    }

    /// The texts of the subjects that the query `ordinal`'s report reads one by one, where its
    /// writer places them (`fetch_by_query` with `places`); empty otherwise.
    pub(crate) fn query_reads(&self, ordinal: usize) -> &crate::report::query_warnings::QueryReads {
        static NONE: crate::report::query_warnings::QueryReads = Vec::new();
        self.query_reads.get(&ordinal).unwrap_or(&NONE)
    }

    /// The messages of every query's own fetches, in query order, for a format that writes no
    /// query warning to standard error (outfmt 5 puts them in the Iteration's message).
    pub(crate) fn take_query_fetches(&mut self) -> Vec<u8> {
        std::mem::take(&mut self.query_fetches)
            .into_values()
            .flatten()
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:844-857
    /// Reads what the custom tabular report `spec` shows of the batch's reported sequences
    /// beyond the default report's ids, before the batch's first byte: the ids of every
    /// definition line of a sequence with parsed Seq-ids (`sall*`, with LOSAT's refusal of
    /// `ParsedSubject::all_tabular_ids`), the taxonomy (`staxid`, `staxids`; from the definition
    /// lines of a sequence with parsed Seq-ids, else from its stored ones, `SeqDb::defline_set`;
    /// a header that does not decode is the reader's error) and the names of the taxonomy
    /// database that the rows look up (`shared::tabular_names`, with its refusals). `batches`
    /// are the batch's hit lists, renamed to ranks by `for_batch`.
    pub(crate) fn read_tabular(
        &mut self,
        database: &BlastxDatabase,
        spec: &crate::report::tabular::TabularSpec,
        batches: &[BatchResults],
        gapped: bool,
    ) -> Result<()> {
        use crate::report::tabular::TabularTaxonomy;
        let needs = spec.needs();
        if !needs.all_deflines && !needs.taxonomy {
            return Ok(());
        }
        let ranks: std::collections::BTreeSet<usize> = batches
            .iter()
            .flat_map(|batch| batch.queries.iter().flatten())
            .filter(|list| super::report::reported(list, gapped))
            .map(|list| list.oid)
            // A subject that no report loads (`HitLoads::loaded`: one whose load throws, whose
            // rows keep the names of the rows before, and one that no report reaches): NCBI reads
            // none of its definition lines and looks none of its taxids up.
            .filter(|&rank| self.loads.loaded(rank))
            .collect();
        let count = self.ordinals.len();
        if needs.all_deflines {
            self.tabular.all_ids = vec![Default::default(); count];
        }
        if needs.taxonomy {
            self.tabular.taxonomy = vec![TabularTaxonomy::default(); count];
        }
        for rank in ranks {
            let ordinal = self.ordinals[rank];
            let parsed = self.parsed.get(rank).and_then(Option::as_deref);
            if let (true, Some(subject)) = (needs.all_deflines, parsed) {
                self.tabular.all_ids[rank] = subject.all_tabular_ids(true).map_err(|problem| {
                    shared::parsed_refusal(
                        self.oids[rank],
                        &format!("has {problem}"),
                        database.program,
                    )
                })?;
            }
            if needs.taxonomy {
                self.tabular.taxonomy[rank] = match parsed {
                    Some(subject) => TabularTaxonomy::of_deflines(subject.deflines()),
                    None => {
                        let set = database
                            .db
                            .defline_set(ordinal)
                            .map_err(|error| shared::database_error(error, database.program))?;
                        TabularTaxonomy::of_deflines(&set.deflines)
                    }
                };
            }
        }
        self.tabular.names = shared::tabular_names(
            &database.db,
            &database.taxdb,
            spec,
            &self.tabular.taxonomy,
            database.program,
        )?;
        Ok(())
    }

    /// The ids of each definition line of the sequence of rank `rank` with parsed Seq-ids
    /// (only after `read_tabular` with an `sall*` field; empty otherwise), and whether
    /// `x_SetSubjectIds` throws after them (`AllTabularIds::failed`).
    pub(crate) fn tabular_all_ids(
        &self,
        rank: usize,
    ) -> (&[crate::report::tabular::DeflineIds], bool) {
        self.tabular
            .all_ids
            .get(rank)
            .map_or((&[], false), |all| (all.ids.as_slice(), all.failed))
    }

    /// `staxid` and `staxids` of the sequence of rank `rank` (only after `read_tabular` with a
    /// taxonomy field; none otherwise).
    pub(crate) fn tabular_taxonomy(
        &self,
        rank: usize,
    ) -> Option<&crate::report::tabular::TabularTaxonomy> {
        self.tabular.taxonomy.get(rank)
    }

    /// The taxonomy database's names that the rows look up (`read_tabular`).
    pub(crate) fn tabular_names(&self) -> Option<&crate::report::tabular::TabularNames> {
        self.tabular.names.as_ref()
    }

    /// The records by rank.
    pub(crate) fn records(&self) -> &[FastaRecord] {
        &self.records
    }

    /// With `-db_soft_mask` or `-db_hard_mask` of the algorithm `algo_id` and a pairwise report:
    /// the stored masks of each rank, which the rows show (`shown_db_masks`). The records'
    /// letters stay those stored.
    #[cold]
    #[inline(never)]
    pub(crate) fn read_db_masks(&mut self, database: &BlastxDatabase, algo_id: i32) -> Result<()> {
        let mut reader = crate::algorithm::common::db_mask_search::SubjectMasks::with_window(1);
        let mut masks = Vec::with_capacity(self.ordinals.len());
        for &ordinal in &self.ordinals {
            masks.push(
                reader
                    .of(
                        &database.db,
                        ordinal,
                        algo_id,
                        shared::database_error_of,
                        database.program,
                    )?
                    .to_vec(),
            );
        }
        self.db_masks = Some(masks);
        Ok(())
    }

    /// The subject masks that the pairwise rows `hits` show in a search with `-db_soft_mask` or
    /// `-db_hard_mask` (`read_db_masks`), by (query, rank), inclusive, by NCBI's Seq-id
    /// (`db_mask_search::shown_db_subject_masks`); `None` for a search without them.
    #[cold]
    #[inline(never)]
    pub(crate) fn shown_db_masks(
        &self,
        hits: &[crate::report::pairwise::PairwiseHit],
    ) -> Option<std::collections::HashMap<(u32, u32), Vec<(i32, i32)>>> {
        use crate::algorithm::common::db_mask_search as db_mask;
        let masks = self.db_masks.as_ref()?;
        Some(db_mask::shown_db_subject_masks(
            hits,
            |rank| masks.get(rank as usize).map(Vec::as_slice),
            |rank| {
                db_mask::subject_id_keys(
                    self.parsed_subject(rank as usize),
                    self.oid(rank as usize),
                )
            },
        ))
    }

    /// The ordinal id of rank `rank`.
    pub(crate) fn ordinal(&self, rank: usize) -> u32 {
        self.ordinals[rank]
    }

    /// The OID of rank `rank` (`gnl|BL_ORD_ID|n`).
    pub(crate) fn oid(&self, rank: usize) -> crate::seqdb::Oid {
        self.oids[rank]
    }

    /// The OID of each rank (`oid`).
    pub(crate) fn oids(&self) -> &[crate::seqdb::Oid] {
        &self.oids
    }

    /// The Blast-def-line taxid of the sequence of rank `rank` (0 for none, and for a run that
    /// does not read taxids).
    pub(crate) fn taxid(&self, rank: usize) -> i32 {
        self.taxids.get(rank).copied().unwrap_or(0)
    }

    /// What the taxonomy database says of that taxid (`HitTaxonomy::Unnamed` when nothing).
    pub(crate) fn taxonomy(&self, rank: usize) -> HitTaxonomy {
        self.taxonomy.get(rank).cloned().unwrap_or_default()
    }

    /// What the taxonomy database says of the taxid of each definition line of the sequence of
    /// rank `rank` with parsed Seq-ids (empty where nothing was asked).
    pub(crate) fn parsed_taxonomy(&self, rank: usize) -> &[HitTaxonomy] {
        self.parsed_taxonomy.get(rank).map_or(&[], Vec::as_slice)
    }

    /// NCBI's messages about the taxonomy database that come before the `Hit`s' own, for the
    /// batch.
    pub(crate) fn taxdb_notes(&self) -> &shared::TaxdbNotes {
        &self.taxdb_notes
    }

    /// The XML2 and JSON formats: the taxid of each reported sequence that the reports load
    /// (`for_batch`), and what the taxonomy database says about it
    /// (`shared::read_reported_taxonomy`; a sequence with parsed Seq-ids has the taxid of each
    /// definition line), read after the batch's prefetch or its reports' reads one by one
    /// (`native::reported_batch`), whose loads write their texts, so that the texts of opening the
    /// database that a `Hit`'s lookup makes are the `Hit`'s (`TaxdbNotes`). A taxonomy database
    /// LOSAT cannot read as NCBI does is LOSAT's refusal of the XML2 document (a run that holds its
    /// output drops it). Nothing for another format.
    pub(crate) fn read_taxonomy(&mut self, database: &BlastxDatabase) -> Result<()> {
        let Some(ranks) = self.xml_ranks.take() else {
            return Ok(());
        };
        (
            self.taxids,
            self.taxonomy,
            self.parsed_taxonomy,
            self.taxdb_notes,
        ) = shared::read_reported_taxonomy(
            &database.db,
            &database.taxdb,
            &self.ordinals,
            &ranks,
            &self.parsed,
            database.program,
        )
        .map_err(super::blast_xml::XmlRefusal::of_taxonomy)?;
        Ok(())
    }

    /// What the batch's reports do at the `Hit` of each subject, by rank (`for_batch`).
    pub(crate) fn loads(&self) -> &HitLoads {
        &self.loads
    }

    /// The tabular `sseqid`, `sacc` and `saccver` of rank `rank` (only with `tabular` in
    /// `for_batch`): the same id for a sequence without parsed Seq-ids whose title gives a word;
    /// `None` where the lookup fails (`common::tabular::SubjectId`).
    pub(crate) fn tabular_ids(&self, rank: usize) -> Option<[&str; 3]> {
        let id = self.tabular_ids[rank].as_deref()?;
        Some(
            match self.tabular_accessions.get(rank).and_then(Option::as_ref) {
                Some([accession, version]) => [id, accession, version],
                None => [id, id, id],
            },
        )
    }

    /// With `-parse_deflines`, `common::database::refuse_id_of_a_query` of each reported sequence
    /// with parsed Seq-ids of the batches `batches` of one reader batch (their hit lists named by
    /// their ranks, `for_batch`), with the scope that `scopes` (numbered from the reader batch's
    /// first query) gives the query of its report (`BatchResults::query_ordinal`, which counts
    /// from the reader batch, and the query's index).
    pub(crate) fn refuse_ids_of_queries(
        &self,
        batches: &[BatchResults],
        gapped: bool,
        scopes: crate::blastinput::query_batch::QueryScopesFrom<'_>,
    ) -> Result<()> {
        if self.parsed.is_empty() {
            return Ok(());
        }
        for batch in batches {
            for (qidx, lists) in batch.queries.iter().enumerate() {
                for list in lists {
                    if !super::report::reported(list, gapped) {
                        continue;
                    }
                    if let Some(subject) = self.parsed_subject(list.oid) {
                        shared::refuse_id_of_a_query(
                            subject,
                            || self.oids[list.oid],
                            scopes.of_query(batch.query_ordinal + qidx),
                            PROGRAM,
                        )?;
                    }
                }
            }
        }
        Ok(())
    }

    /// The sequences with parsed Seq-ids by rank (empty when the batch reports none).
    pub(crate) fn parsed(&self) -> &[Option<Arc<ParsedSubject>>] {
        &self.parsed
    }

    /// The sequence of rank `rank` when it has parsed Seq-ids.
    pub(crate) fn parsed_subject(&self, rank: usize) -> Option<&ParsedSubject> {
        self.parsed.get(rank).and_then(Option::as_deref)
    }

    /// The stored titles by rank (the outfmt 0 alignment headings).
    pub(crate) fn stored_titles(&self) -> &[Arc<[u8]>] {
        &self.titles
    }

    /// The sources from which NCBI makes the empty stored titles, by rank (empty when none has
    /// one; `report::pairwise::stored_source`).
    pub(crate) fn sources(&self) -> &[crate::report::defline::SourceTaxname] {
        &self.sources
    }

    /// The description-table labels by rank: the title as NCBI's defline generator makes it
    /// for the description table (as BLASTP's, BLASTN's and TBLASTN's `-db` reports), or the
    /// Seq-align id and `Unknown` of one whose load throws (`HitLoads::thrown`).
    pub(crate) fn descriptions(&self) -> Vec<crate::report::pairwise::DescriptionLabel> {
        self.titles
            .iter()
            .enumerate()
            .map(|(rank, title)| {
                crate::report::pairwise::database_description_label(
                    title,
                    crate::report::pairwise::stored_source(&self.sources, rank as u32),
                    self.parsed_subject(rank),
                    true,
                    self.loads.thrown(rank).is_some(),
                )
            })
            .collect()
    }
}

#[cfg(test)]
#[path = "database_tests.rs"]
mod tests;
