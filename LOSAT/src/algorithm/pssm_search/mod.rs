//! The search of a query given as a position-specific score matrix (PSSM): psiblast `-in_pssm`
//! (series PS, items PS4a, `-comp_based_stats 0`, and PS4b, the position-based composition-based
//! statistics of modes 1 to 3).
//!
//! The kernels read the PSSM's rows by absolute query position where blastp's read the matrix
//! row of the query letter (`PssmRows`): the lookup table (`lookup.rs`), the ungapped extensions
//! (`ungapped.rs`), the gapped start, score-only and traceback alignments (`gapalign.rs`), all
//! scalar as NCBI's PSSM code is. The stages (`prelim.rs`, `traceback0.rs`, `redo1.rs`) run them for the one
//! query of a PSSM search; blastp's per-run driver calls them where a PSSM search differs
//! (`blastp::blast_engine`, the `psi_search` hooks) and does the rest (the subject sets and
//! databases, the hit lists, the reports) as for blastp.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2826
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:113-126

pub(crate) mod annot_index;
pub(crate) mod evalue;
pub(crate) mod gapalign;
pub(crate) mod in_pssm;
pub(crate) mod lookup;
pub(crate) mod prelim;
pub(crate) mod redo1;
pub(crate) mod rows;
pub(crate) mod traceback0;
pub(crate) mod ungapped;

use crate::blastinput::fasta_reader::FastaRecord;
use crate::stats::KarlinParams;

/// The statistical values of a PSSM's `finalData`, each `None` where NCBI's getter gives
/// `CPssm::kInvalidStat` (the member is absent, or 0.0).
///
/// NCBI reference (598d8ae6): c++/src/objects/scoremat/Pssm.cpp:101-155
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(crate) struct PssmStats {
    pub(crate) lambda_ungapped: Option<f64>,
    pub(crate) kappa_ungapped: Option<f64>,
    pub(crate) h_ungapped: Option<f64>,
    pub(crate) lambda: Option<f64>,
    pub(crate) kappa: Option<f64>,
    pub(crate) h: Option<f64>,
}

/// A query given as a PSSM (psiblast `-in_pssm`): its rows, the statistical values of its file,
/// the query record made from its Bioseq, and NCBI's errors of the PSSM before its search
/// (`stop_before_search`: the check of the PSSM, the PSSM engine's scores), which NCBI raises after
/// the outfmt 0 prolog.
#[derive(Debug, Clone)]
pub struct PsiSearch {
    pub(crate) pssm: rows::QueryPssm,
    /// The PSSM's frequency ratios, 28 per position (`psi_matrix->freq_ratios`), which the
    /// position-based composition-based statistics read; empty when the PSSM has none.
    pub(crate) freq_ratios: Vec<[f64; crate::core::blast_util::BLASTAA_SIZE]>,
    pub(crate) stats: PssmStats,
    /// The query record as the reports show it without `-parse_deflines` (`query_record`).
    pub(crate) query: FastaRecord,
    /// The letters of the query's Bioseq as NCBI's reports read them where they are not the
    /// record's (`query`): a PSSM query (`-in_pssm`) of iupacaa letters with one in lower case,
    /// which the search reads as its capital (the record's letter) and the reports show and
    /// compare as it is (`in_pssm::query_residues`, `report::row_letters::QueryLetters`). `None`
    /// otherwise. The query source is fixed for the run: a PSSM query's later rounds keep it
    /// (`blast_engine::psi_round`); a FASTA query's rounds have none.
    pub(crate) bioseq_letters: Option<std::sync::Arc<[u8]>>,
    /// The Seq-ids of the PSSM's query Bioseq of `-in_pssm` (at least one), the one owner of them:
    /// the first is the id of the query's Seq-loc (psiblast_app.cpp:540, `GetFirstId()`), which
    /// its Seq-aligns, the later rounds' PSSMs and the warnings take; the reports show them as
    /// `query_record` makes them. Empty for a FASTA query's later rounds (`from_engine`).
    pub(crate) query_ids: Vec<crate::seqid::ParsedSeqId>,
    /// `CPsiBlastValidate::Pssm`'s error (`in_pssm::validation_error`).
    pub(crate) validation_error: Option<String>,
    /// The PSSM engine's error for the scores of a PSSM of frequency ratios only
    /// (`PsiBlastComputePssmScores` in `CPsiBlastImpl`'s constructor; `in_pssm::psi_search_and_object`),
    /// which comes after the search strategy's request (`PssmQueryStop::Scores`).
    pub(crate) scores_error: Option<String>,
    /// The error of the query set-up of the PSSM's search (`SetupQueries_OMF`: the query's invalid
    /// residues, `in_pssm::query_residues`), which NCBI raises in `CPsiBlast::Run` after the search
    /// strategy's request and, in psiblast, after round 1's PSSM files (`PssmQueryStop::Setup`). A
    /// PSSM of frequency ratios only meets that set-up in `CPsiBlastImpl`'s constructor, where the
    /// PSSM engine reads the query: its error is then `scores_error`.
    pub(crate) setup_error: Option<String>,
    /// The warnings of `PsiBlastSetupScoreBlock` for the search's options (the frequency ratios
    /// all zero, the composition-based statistics reset to mode 1), which the query's warning
    /// line carries (`algorithm::psiblast::in_pssm::setup_messages`).
    pub(crate) setup_messages: Vec<String>,
    /// The title of the query's warning label when it is not the record's: a PSSM query's
    /// (`CBlastQuerySourceBioseqSet::GetTitle`, `algorithm::pssm_search::in_pssm`); `None` for a FASTA
    /// query's later rounds, whose label keeps the record's title.
    pub(crate) warning_title: Option<Vec<u8>>,
    /// What the search strategy's request holds of the PSSM (`PssmQueryRequest`): kept only for a
    /// run that writes one (`in_pssm::PssmReader::strategy`, `pssm_request_written`) of a PSSM that
    /// NCBI's check passes; `None` otherwise and for a FASTA query's later rounds.
    pub(crate) request: Option<PssmQueryRequest>,
    /// The PSSM that psiblast's archive (outfmt 11) of a round of this search holds in its request
    /// before the round's K/A write-back (`algorithm::blastp::archive::PsiArchive`): `-in_pssm`'s
    /// object after `CPsiBlastImpl`'s constructor (the scores of a PSSM of frequency ratios only
    /// computed, `in_pssm::InPssmFile::object`), or the engine's PSSM of a later round
    /// (`algorithm::psiblast::pssm::QueryPssm::object`). Kept only for a run that writes the
    /// archive; `None` otherwise and for tblastn.
    pub(crate) archive_pssm:
        Option<std::sync::Arc<crate::report::serial::pssm::PssmWithParameters>>,
    /// What NCBI's object manager makes of the annotations of the PSSM's query when the query
    /// joins the scope (`annot_index::index_entry`): its messages, the diagnostic stream it
    /// leaves failed and its exception, which the run writes where NCBI posts them
    /// (`annot_posts`, `annot_abort`); empty without annotations and for a FASTA query's rounds.
    pub(crate) annot: annot_index::AnnotOutcome,
    /// Whether NCBI's definition line of the query, which its XML and JSON reports make, reads
    /// the query's annotations (`in_pssm::defline_reads_annot`).
    pub(crate) defline_reads_annot: bool,
    /// The query that a later round of psiblast searches (`blast_engine::psi_round::RoundQuery`):
    /// the letters that the rounds before masked (X), filtered again without the lower-case mask
    /// that round 1 took; round 1's query for a query that NCBI splits into chunks
    /// (`split_setup`); `None` for a query's first round, whose query frame is its record's
    /// (`blast_engine::blastp_query_frames`).
    pub(crate) round_query: Option<crate::algorithm::tblastx::translation::QueryFrame>,
    /// A later round of a FASTA query that NCBI splits into chunks: the query of the round's
    /// set-up, whose `kbp_std` the round reads (`blast_engine::pssm_hooks::pssm_lookup`): the
    /// query that the round before searched, filtered again as `RoundQuery::next` filters it.
    /// The round's chunks search without the PSSM, and its traceback reads `round_query`, which
    /// NCBI's preliminary stage makes again from the record after the chunks
    /// (`blast_engine::psi_round::RoundQuery::split`). `None` otherwise: a query that NCBI does
    /// not split, and a PSSM query (`-in_pssm`), which NCBI never splits.
    pub(crate) split_setup: Option<crate::algorithm::tblastx::translation::QueryFrame>,
    /// The PSSM was made from psiblast's `-in_msa` alignment, not read from `-in_pssm`: the
    /// refusals that name the option name that one.
    pub(crate) from_msa: bool,
    /// The query's Seq-inst has the molecule `other` (psiblast only), with the coding number of
    /// its Seq-data (`CSeq_data::E_Choice`: iupacaa 2, ncbieaa 8, ncbistdaa 10). NCBI searches it
    /// as a protein; the tabular formats read the query's letters through the object manager for
    /// the fields `qseq`, `sseq`, `btop`, `positive`, `ppos`, `nident`, `mismatch` and `pident`
    /// (the default fields among them), which throws `CSeqVectorException` at the first row
    /// (`report::tabular::query_sequence_exception`).
    pub(crate) query_other_coding: Option<i32>,
}

/// The PSSM of `-in_pssm` as NCBI's application read it (before `CPsiBlastImpl` computes the
/// scores of a PSSM of frequency ratios only), which the search strategy's request holds
/// (`report::archive::requests::RunRequests::with_pssm_query`) and writes back whole
/// (`report::serial::pssm::write_value`).
#[derive(Debug, Clone)]
pub(crate) struct PssmQueryRequest {
    pub(crate) pssm: std::sync::Arc<crate::report::serial::pssm::PssmWithParameters>,
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:279-304
/// Whether NCBI writes the search strategy's request of a PSSM query (`s_InitializeExportStrategy`
/// of psiblast and tblastn): with `-export_search_strategy` (`strategy`) and a database
/// (`database`: `-db`); with `-subject` it stops instead ("Remote PSI-BL2SEQ is not supported").
pub(crate) fn pssm_request_written(strategy: bool, database: bool) -> bool {
    strategy && database
}

/// What stops a PSSM query's run before its search, after the outfmt 0 prolog, in NCBI's order
/// (`PsiSearch::stop_before_search`): the search strategy's request of a run with `-subject`, NCBI's
/// check of the PSSM, the PSSM engine's scores, then the search's set-up of the query.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum PssmQueryStop {
    /// "Remote PSI-BL2SEQ is not supported" (`s_InitializeExportStrategy`): the options error, no
    /// request.
    Bl2seq,
    /// `CPsiBlastValidate::Pssm`'s error, which the request makes too before it is written
    /// (`x_Process_Pssm`): the engine error, no request.
    Invalid(String),
    /// The PSSM engine's error for the scores of a PSSM of frequency ratios only, in
    /// `CPsiBlastImpl`'s constructor, after `SaveSearchStrategy` wrote the request: the engine
    /// error, after the request (`after_request`).
    Scores(String),
    /// The error of the query set-up (`PsiSearch::setup_error`), in `CPsiBlast::Run`: the engine
    /// error, after the request and after psiblast's round 1 wrote the input PSSM's files
    /// (`after_round_files`).
    Setup(String),
}

impl PssmQueryStop {
    /// NCBI's error (exit status 1 for the options error, 3 for the engine's).
    pub(crate) fn error(&self) -> anyhow::Error {
        match self {
            PssmQueryStop::Bl2seq => {
                crate::blastinput::app::options_error("Remote PSI-BL2SEQ is not supported")
            }
            PssmQueryStop::Invalid(message)
            | PssmQueryStop::Scores(message)
            | PssmQueryStop::Setup(message) => crate::blastinput::app::engine_error(message),
        }
    }

    /// Whether NCBI has written the search strategy's request when it stops (a run that writes
    /// one, `pssm_request_written`).
    pub(crate) fn after_request(&self) -> bool {
        matches!(self, PssmQueryStop::Scores(_) | PssmQueryStop::Setup(_))
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:246-262
    /// Whether psiblast's round 1 has written the input PSSM's files (`SavePssmToFile`, `-out_pssm`
    /// and `-out_ascii_pssm`) when NCBI stops: the query set-up fails in `CPsiBlast::Run`, after
    /// them; the other stops come before the rounds.
    pub(crate) fn after_round_files(&self) -> bool {
        matches!(self, PssmQueryStop::Setup(_))
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:537
// NCBI's scope holds the subjects that its reader read, then the PSSM's Bioseq; the command line
// reads its `-subject` records in the scope that the query joins (`IdScope`), while the records'
// entries (the library API) get records read in a scope of their own (`api::read_inputs`).
/// The scope of the ids of `subjects` (read under `-parse_deflines` in another scope) that the
/// query of a PSSM search joins in the records' entries, as the command line's query joins the
/// scope it read its subjects in (`PsiSearch::query_record_in`); `None` without
/// `-parse_deflines`. The subjects' marks are not kept: the PSSM's query has none, so a subject
/// with one of its ids is refused whatever its mark (`IdScope::add_record`), as on the command
/// line.
pub(crate) fn scope_of_subjects(
    parse_deflines: bool,
    subjects: &[FastaRecord],
    program: &str,
) -> anyhow::Result<Option<crate::seqid::id_scope::IdScope>> {
    if !parse_deflines {
        return Ok(None);
    }
    let scope = crate::seqid::id_scope::IdScope::default();
    for believed in subjects.iter().filter_map(FastaRecord::believed) {
        scope.add_record(&believed.ids, None, "subject", program)?;
    }
    Ok(Some(scope))
}

impl PsiSearch {
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:338-340 and
    // c++/src/app/blast/tblastn_app.cpp:307,324 (the request, then the search object)
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:533-534
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_impl.cpp:58-67
    /// What stops the run of this PSSM query before its search (psiblast's and tblastn's one rule):
    /// with `-export_search_strategy` (`strategy`) and no database (`database`), the request's
    /// "Remote PSI-BL2SEQ" before anything else; then NCBI's check of the PSSM (the request makes
    /// it, and `CPsiBlastImpl` without one); then the PSSM engine's scores of a PSSM of frequency
    /// ratios only, which `CPsiBlastImpl` computes after the request was written; then the query
    /// set-up that `CPsiBlast::Run` starts with (`setup_error`). `None` for a run that searches.
    pub(crate) fn stop_before_search(
        &self,
        strategy: bool,
        database: bool,
    ) -> Option<PssmQueryStop> {
        if strategy && !pssm_request_written(strategy, database) {
            return Some(PssmQueryStop::Bl2seq);
        }
        if let Some(message) = &self.validation_error {
            return Some(PssmQueryStop::Invalid(message.clone()));
        }
        if let Some(message) = &self.scores_error {
            return Some(PssmQueryStop::Scores(message.clone()));
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:149-151
        // The search's first step (`CBlastPrelimSearch`'s set-up in `CPsiBlast::Run`) sets the
        // query up (`SetupQueries_OMF`), before the score block, the lookup table and any read of
        // the database.
        self.setup_error.clone().map(PssmQueryStop::Setup)
    }

    /// The PSSM that the search strategy's request holds (`PssmQueryRequest`), `None` without one.
    pub(crate) fn request_pssm(
        &self,
    ) -> Option<std::sync::Arc<crate::report::serial::pssm::PssmWithParameters>> {
        self.request.as_ref().map(|request| request.pssm.clone())
    }

    /// `kbp_psi[0]`: the file's ungapped values, or those of `kbp_std` (the query's composition
    /// with the search's matrix, which NCBI also computes for `kbp_psi` before it overwrites them).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:89-107
    /// `kbp_psi[0]` holds `kbp_std[0]`'s values before (`Blast_ScoreBlkKbpUngappedCalc`,
    /// blast_stat.c:2799-2808: the same calculation), so either branch leaves those values.
    pub(crate) fn ungapped_karlin(&self, kbp_std: &KarlinParams) -> KarlinParams {
        KarlinParams {
            lambda: self.stats.lambda_ungapped.unwrap_or(kbp_std.lambda),
            k: self.stats.kappa_ungapped.unwrap_or(kbp_std.k),
            h: self.stats.h_ungapped.unwrap_or(kbp_std.h),
            ..*kbp_std
        }
    }

    /// Whether `kbp_psi[0]` (`ungapped_karlin`) is valid: a file's negative ungapped value makes
    /// it invalid, and the search's setup then stops with NCBI's engine error of an invalid query
    /// (`s_BlastFindValidKarlinBlk` over `sbp->kbp` set to `kbp_psi`, blast_stat.c:2826).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:43-50
    pub(crate) fn ungapped_karlin_valid(&self, kbp_std: &KarlinParams) -> bool {
        let kbp = self.ungapped_karlin(kbp_std);
        kbp.lambda > 0.0 && kbp.k > 0.0 && kbp.h > 0.0
    }

    /// `kbp_gap_psi[0]`: the file's gapped values, or those of `kbp_gap_std` (the table's values
    /// of the search's matrix and gap costs, which `kbp_gap_psi` is a copy of).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:109-127
    pub(crate) fn gapped_karlin(&self, kbp_gap_std: &KarlinParams) -> KarlinParams {
        KarlinParams {
            lambda: self.stats.lambda.unwrap_or(kbp_gap_std.lambda),
            k: self.stats.kappa.unwrap_or(kbp_gap_std.k),
            h: self.stats.h.unwrap_or(kbp_gap_std.h),
            ..*kbp_gap_std
        }
    }

    /// The query record made from the PSSM's Bioseq.
    pub fn query(&self) -> &FastaRecord {
        &self.query
    }

    /// The query record of the PSSM search of `-in_pssm` as the reports show it (psiblast's,
    /// tblastn's PSI-TBLASTN): the record of the PSSM's Bioseq (`in_pssm::psi_search_and_object`),
    /// whose one local id `-parse_deflines` shows as it shows a FASTA query's ids (NCBI's formatter
    /// takes `believe_query` of `-parse_deflines`, psiblast_app.cpp:482-484, tblastn_app.cpp; NCBI
    /// BLAST+ 2.17.0: `Query= Query_1 <title>`, `qseqid` `Query_1`).
    pub(crate) fn query_record(&self, parse_deflines: bool) -> FastaRecord {
        match (self.query_ids.as_slice(), parse_deflines) {
            ([crate::seqid::ParsedSeqId::Local(_)], true) => {
                self.query.clone().believing(self.query_ids.clone())
            }
            _ => self.query.clone(),
        }
    }

    /// `query_record`, whose ids, when the reports show them, join the scope `scope` of the
    /// `-subject` records that `-parse_deflines` reads, as a FASTA query's do (NCBI's scope holds
    /// the PSSM's Bioseq, psiblast_app.cpp:537: a subject with one of its ids shows as the query, or
    /// NCBI's tblastn fails; `IdScope::add_record` refuses it) in a run that searches (`searches`, of
    /// `PsiSearch::searches`):
    /// the subjects' Seq-ids matter to NCBI only in the reports, so a run that NCBI ends before the
    /// search (the object manager's stop, the PSSM's checks, the PSSM engine's error) meets no
    /// subject there and writes NCBI's stop (NCBI BLAST+ 2.17.0, line RG D3-1 review 4b S1 (c),
    /// review 4c N1).
    pub(crate) fn query_record_in(
        &self,
        parse_deflines: bool,
        scope: Option<&crate::seqid::id_scope::IdScope>,
        program: &str,
        searches: bool,
    ) -> anyhow::Result<FastaRecord> {
        let query = self.query_record(parse_deflines);
        if let (Some(scope), Some(believed), true) = (scope, query.believed(), searches) {
            scope.add_record(&believed.ids, None, "query", program)?;
        }
        Ok(query)
    }

    /// Whether the PSSM's query joins NCBI's scope: false where NCBI's object manager ends the run
    /// as it adds the query's Seq-entry (`annot.abort`, `scope->AddTopLevelSeqEntry`, and
    /// `GetSeq()` of a set, psiblast_app.cpp:537-540, tblastn_app.cpp:226-228). Such a run makes
    /// nothing that NCBI makes after it (tblastn's formatter and prolog, the search strategy, the
    /// search, its reports and PSSM files), so none of LOSAT's limits of those apply to it.
    pub(crate) fn joins_scope(&self) -> bool {
        self.annot.abort.is_none()
    }

    /// Whether the run searches the PSSM: its query joins the scope (`joins_scope`) and NCBI does
    /// not stop before the search (`stop_before_search`, `strategy` and `database` as there).
    pub(crate) fn searches(&self, strategy: bool, database: bool) -> bool {
        self.joins_scope() && self.stop_before_search(strategy, database).is_none()
    }

    /// The query that NCBI's scope holds for the PSSM of `-in_pssm`: the PSSM's Bioseq, which
    /// psiblast and tblastn add to the scope with its Seq-ids `query_ids` before the search
    /// (`query_batch::QueryScope::of_bioseq_ids`; psiblast_app.cpp:537, tblastn_app.cpp:226), for
    /// the stop of a batch's prefetch (`common::database::PrefetchQueries`). A later round's PSSM
    /// of a FASTA query (`from_engine`) has no ids of its own: the scope then holds the query's
    /// batch (`blastp::blast_engine::psi_round`).
    pub(crate) fn query_scope(&self) -> crate::blastinput::query_batch::QueryScope {
        crate::blastinput::query_batch::QueryScope::of_bioseq_ids(&self.query_ids)
    }

    /// The id of the query's Seq-loc of `-in_pssm` (`query_ids`), `None` for a FASTA query.
    pub(crate) fn query_seq_loc_id(&self) -> Option<&crate::seqid::ParsedSeqId> {
        self.query_ids.first()
    }

    /// The PSSM search of a later round of the FASTA query `query` (psiblast's rounds,
    /// `algorithm::psiblast::rounds`): the engine's PSSM of the round before, in full precision
    /// (NCBI's `CPssmWithParameters` in memory, never read back from text): its scores, its
    /// frequency ratios and the six statistical values of `x_PSIMatrix2Asn1`, each `None` where
    /// it is 0.0 (`CPssm::kInvalidStat`), as `psi_search_of` reads a file's; and the warnings of
    /// `PsiBlastSetupScoreBlock` for the round's options (`setup_messages`). The query stays the
    /// FASTA query: `CPsiBlastImpl::SetPssm` replaces only the PSSM.
    ///
    /// With `-query_loc` the query is round 1's range (`query`, cut from its record) while the
    /// PSSM is the whole record's (`s_GetQueryBioseq` gives the record's Bioseq to the PSSM
    /// engine, `blast_engine::psi_round::next_pssm`): the search reads the PSSM's rows at the
    /// query's offsets, 0 to the range's length, so the range's first letter meets the record's
    /// first row (NCBI 2.17.0 `-query_loc 20-300 -num_iterations 3`: round 2 aligns rows 19-190,
    /// round 1's alignment, as the record's letters 39-210 at 4.494 % identities; series Q, Q4a).
    /// No row past the range is read: the kernels, the lookup table of the query's offsets and
    /// the position-based composition-based statistics index the PSSM by the query's offset in
    /// its range (`pssm_search::rows`), and the Karlin values are the PSSM's own (NCBI's debug
    /// build asserts that the PSSM has the query's length in the redo, blast_kappa.c:3066; the
    /// release build, the reference, does not). Without `-query_loc` the query is the record and
    /// every row is read. `PsiBlastSetupScoreBlock`'s warnings look at every row
    /// (`setup_messages`, made by the caller from the engine's PSSM).
    ///
    /// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:303-310
    /// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:157-161 (`s_GetQueryBioseq`)
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_impl.cpp:173-181
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:648-660 (cited at
    /// `algorithm::psiblast::pssm::pssm_with_parameters`: the scores and values that the next
    /// round's `PsiBlastSetupScoreBlock` reads)
    pub(crate) fn from_engine(
        pssm: &crate::core::psi::PsiPssm,
        query: FastaRecord,
        setup_messages: Vec<String>,
    ) -> Self {
        let stat = |value: f64| (value != 0.0).then_some(value);
        let read = query.seq().len();
        Self {
            pssm: rows::QueryPssm {
                rows: pssm.scores[..read].to_vec(),
            },
            freq_ratios: pssm
                .diagnostics
                .frequency_ratios
                .as_deref()
                .map_or_else(Vec::new, |ratios| ratios[..read].to_vec()),
            stats: PssmStats {
                lambda_ungapped: stat(pssm.ungapped.lambda),
                kappa_ungapped: stat(pssm.ungapped.k),
                h_ungapped: stat(pssm.ungapped.h),
                lambda: stat(pssm.gapped.lambda),
                kappa: stat(pssm.gapped.k),
                h: stat(pssm.gapped.h),
            },
            query,
            bioseq_letters: None,
            query_ids: Vec::new(),
            validation_error: None,
            scores_error: None,
            setup_error: None,
            setup_messages,
            warning_title: None,
            request: None,
            archive_pssm: None,
            annot: annot_index::AnnotOutcome::default(),
            defline_reads_annot: false,
            round_query: None,
            split_setup: None,
            from_msa: false,
            query_other_coding: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:506, 537 and
    // c++/src/app/blast/tblastn_app.cpp:226, 270 (psiblast posts the object manager's messages after
    // its prolog, tblastn before it: `formatter.PrintProlog()` comes before
    // `scope->AddTopLevelSeqEntry` of the query (`pssm->SetQuery`) in psiblast, after it in tblastn)
    /// The messages of NCBI's object manager about the PSSM query's annotations, written to
    /// `diagnostics` (psiblast after its prolog, as the messages of its query batch,
    /// `blast_engine::psiblast_stop_at_pssm_query`; tblastn before its prolog), each saved as one
    /// message of its severity by the BLAST archive of the run (`archive`, outfmt 11: the error
    /// archive of a run that the object manager's exception ends, NCBI BLAST+ 2.17.0 tblastn
    /// `k02_failed_then_choice_exc`), after which a failed diagnostic stream silences the run's
    /// later messages (`ncbi_exception::set_main_diag_failbit`).
    pub(crate) fn write_annot_posts(
        &self,
        diagnostics: &mut dyn std::io::Write,
        archive: Option<&crate::report::archive::ArchiveHandle>,
    ) -> std::io::Result<()> {
        self.annot.posts.write(diagnostics, archive)?;
        self.fail_diag_stream();
        Ok(())
    }

    /// Marks the main thread's diagnostic stream failed when the object manager's messages
    /// failed it (`AnnotOutcome::muted`), once they are written.
    pub(crate) fn fail_diag_stream(&self) {
        if self.annot.muted {
            crate::ncbi_exception::set_main_diag_failbit();
        }
    }

    /// The exception with which NCBI's object manager ends the run (`CATCH_ALL`'s report).
    pub(crate) fn annot_abort(&self) -> Option<anyhow::Error> {
        self.annot.abort.as_ref().map(|abort| {
            crate::cli::NativeError {
                exit: abort.exit,
                message: abort.message.clone(),
            }
            .into()
        })
    }
}
