//! psiblast's PSSM after the last round (`-num_iterations 1 -save_pssm_after_last_round` with
//! `-out_pssm F` and/or `-out_ascii_pssm F`): the input of NCBI's PSSM engine from a query's
//! round-1 alignments (`CPsiBlastInputData`), the engine's score block (`CPssmEngine`), the
//! ASN.1 PSSM (`x_PSIMatrix2Asn1`, written by `report::serial::pssm`), the ASCII PSSM
//! (`CBlastFormatUtil::PrintAsciiPssm`) and their files (`CAutoOutputFileReset`). The engine
//! itself is `core::psi`. The blastp engine reaches this module only when a PSSM is requested
//! (`blast_engine.rs` `psiblast_pssms`).

use crate::wording::OPTIONS_ERROR;
use std::io::Write as _;

use anyhow::Result;

use super::rounds::SetKey;
use crate::algorithm::pssm_search::in_pssm::{score_block, NCBISTDAA_TO_AMINOACID};
use crate::common::GapEditOp;
use crate::config::ScoringMatrix;
use crate::core::psi::{
    self, KarlinBlk, PsiDiagnosticsRequest, PsiError, PsiMsa, PsiMsaCell, PsiPssm,
};
use crate::report::diag_log::Posts;
use crate::report::serial::encoder::ReplacedChar;
use crate::report::serial::pssm::{
    FormatRpsDbParameters, ObjectId, Pssm, PssmFinalData, PssmIntermediateData, PssmParameters,
    PssmQuery, PssmWithParameters, SeqData, SeqEntry, SeqInst, Seqdesc,
};
use crate::seqid::ParsedSeqId;
use crate::utils::matrix::{aa_char_to_ncbistdaa, ncbistdaa, BLASTAA_SIZE, BLAST_SCORE_MIN};

/// psiblast's PSSM options as given on the command line.
#[derive(Debug, Clone, Default)]
pub struct PsiblastPssmOptions {
    pub save_pssm_after_last_round: bool,
    pub save_each_pssm: bool,
    pub out_pssm: Option<std::path::PathBuf>,
    pub out_ascii_pssm: Option<std::path::PathBuf>,
    pub pseudocount: i32,
    /// `-in_pssm`: the query as a PSSM file (`algorithm::psiblast::in_pssm`).
    pub in_pssm: Option<std::path::PathBuf>,
    /// `-in_msa`: the query as an alignment, made into a PSSM with the options
    /// (`algorithm::psiblast::in_msa`).
    pub in_msa: Option<super::in_msa::MsaInput>,
    /// The run's `-save_each_pssm` numbers (`PsiblastArgs::pssm_versions`).
    pub versions: PssmVersions,
    /// The queries that the call's checks of its PSSM files cover (`PsiblastArgs::pssm_check`).
    pub check_scope: PssmCheckScope,
    /// The rounds that the round driver runs (`algorithm::psiblast::rounds`): `-num_iterations`
    /// when it is not 1 (0: until convergence), or 1 with `-in_pssm` and a PSSM file (the input
    /// PSSM is written before round 1); `None` for the one search of today's call.
    pub rounds: Option<u32>,
    /// Where the run writes its PSSM output files (`RunFiles::Disk` on CLI/library,
    /// `RunFiles::Receiver` in LOSAT Web). No NCBI counterpart: NCBI writes files.
    pub run_files: crate::report::run_files::RunFiles,
    /// `-phi_pattern`: the pattern file of a PHI-BLAST round 1 (`algorithm::phi_search`), read
    /// with the options (`BlastpArgs::check_options`).
    pub phi_pattern: Option<std::path::PathBuf>,
}

impl PartialEq for PsiblastPssmOptions {
    fn eq(&self, other: &Self) -> bool {
        self.save_pssm_after_last_round == other.save_pssm_after_last_round
            && self.save_each_pssm == other.save_each_pssm
            && self.out_pssm == other.out_pssm
            && self.out_ascii_pssm == other.out_ascii_pssm
            && self.pseudocount == other.pseudocount
            && self.in_pssm == other.in_pssm
            && self.in_msa == other.in_msa
            && self.versions == other.versions
            && self.check_scope == other.check_scope
            && self.rounds == other.rounds
            && self.phi_pattern == other.phi_pattern
    }
}

impl Eq for PsiblastPssmOptions {}

/// A run that writes the PSSM of each query whose round-1 results have alignments.
#[derive(Debug, Clone)]
pub struct PssmRequest {
    /// `-out_pssm`, as given: the ASN.1 text PSSM (`SaveCheckpoint`).
    pub asn_path: Option<std::path::PathBuf>,
    /// `-out_ascii_pssm`, as given (`SaveAsciiPssm`).
    pub ascii_path: Option<std::path::PathBuf>,
    /// `-save_each_pssm`: the files are `F.1`, `F.2`, ... (`CAutoOutputFileReset` with versions).
    pub save_each: bool,
    /// `-pseudocount` (`PSIBlastOptions::pseudo_count`).
    pub pseudocount: i32,
    /// The run's numbers of the files of `-save_each_pssm` (`PssmVersions`).
    pub versions: PssmVersions,
    /// Where the run writes its PSSM output files (`RunFiles::Disk` on CLI/library,
    /// `RunFiles::Receiver` in LOSAT Web). No NCBI counterpart: NCBI writes files.
    pub run_files: crate::report::run_files::RunFiles,
}

impl PartialEq for PssmRequest {
    fn eq(&self, other: &Self) -> bool {
        self.asn_path == other.asn_path
            && self.ascii_path == other.ascii_path
            && self.save_each == other.save_each
            && self.pseudocount == other.pseudocount
            && self.versions == other.versions
    }
}

impl Eq for PssmRequest {}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1304-1318
/// The next number of each of a run's two PSSM streams under `-save_each_pssm`
/// (`CAutoOutputFileReset::m_Version` of `-out_pssm`'s and of `-out_ascii_pssm`'s, from 1). NCBI's
/// argument object makes the two streams once per run, so their numbers go on over the run's
/// queries and rounds. The clones of one `PssmVersions` share its numbers: the library's
/// `run_batches`, whose query batches are calls of their own, gives every batch the run's
/// (`api::batches`); a command line's run has one.
#[derive(Debug, Clone)]
pub struct PssmVersions {
    asn: Version,
    ascii: Version,
}

/// The next number of one PSSM stream, shared by the clones of its `PssmVersions`.
type Version = std::sync::Arc<std::sync::atomic::AtomicU32>;

impl Default for PssmVersions {
    fn default() -> Self {
        let first = || Version::new(std::sync::atomic::AtomicU32::new(1));
        Self {
            asn: first(),
            ascii: first(),
        }
    }
}

impl PssmVersions {
    /// The next numbers of the two streams.
    fn next(&self) -> [u32; 2] {
        use std::sync::atomic::Ordering;
        [
            self.asn.load(Ordering::SeqCst),
            self.ascii.load(Ordering::SeqCst),
        ]
    }
}

/// Which queries a call's check of its PSSM files before its search covers
/// (`blastp::blast_engine::check_psiblast_pssm_files`: the names, `check_pssm_files`; a query title
/// with a byte that NCBI's stream replaces is written as NCBI writes it, `file_posts`). The command
/// line reads all its queries and checks them once, before its first byte; the library's
/// `run_batches`, whose batches are calls of their own, has its first batch make the run's check,
/// where the command line makes it, and the later batches none (`api::batches`).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub enum PssmCheckScope {
    /// The call's own records (the command line, the library's `run`, `run_collect`,
    /// `run_to_file`).
    #[default]
    Call,
    /// The run's: its number of query records.
    Run { queries: usize },
    /// Checked by the run's first batch.
    Checked,
}

/// Two `PssmVersions` are the same where their next numbers are.
impl PartialEq for PssmVersions {
    fn eq(&self, other: &Self) -> bool {
        self.next() == other.next()
    }
}

impl Eq for PssmVersions {}

impl PssmRequest {
    /// The run's streams of `-out_pssm` and of `-out_ascii_pssm` (`CAutoOutputFileReset`), each
    /// where it is requested: with `-save_each_pssm` they number their files with the run's
    /// numbers (`versions`).
    pub(crate) fn streams(&self) -> [Option<OutputFileReset>; 2] {
        let stream = |path: Option<&std::path::Path>, version: &Version| {
            path.map(|path| OutputFileReset::new(path, self.save_each.then(|| version.clone())))
        };
        [
            stream(self.asn_path.as_deref(), &self.versions.asn),
            stream(self.ascii_path.as_deref(), &self.versions.ascii),
        ]
    }

    /// The diagnostics the engine computes (`PSIDiagnosticsRequestNewEx` with `SaveAsciiPssm`).
    pub(crate) fn diagnostics(&self) -> PsiDiagnosticsRequest {
        PsiDiagnosticsRequest::new_ex(self.ascii_path.is_some())
    }
}

impl PsiblastPssmOptions {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1283-1303
    // `CPsiBlastArgs::ExtractAlgorithmOptions` runs after the composition-based statistics
    // handler's (psiblast_args.cpp:102-110, the order of `m_Args`) and before `Validate`
    // (blast_args.cpp:3826-3832).
    pub(crate) fn check(&self) -> Result<()> {
        let no_file = self.out_pssm.is_none() && self.out_ascii_pssm.is_none();
        if self.save_pssm_after_last_round && no_file {
            return Err(crate::blastinput::app::options_error(
                "save_pssm_after_last_round option requires out_pssm or out_ascii_pssm",
            ));
        }
        if self.save_each_pssm && no_file {
            return Err(crate::blastinput::app::options_error(
                "save_each_pssm option requires out_pssm or out_ascii_pssm",
            ));
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:304-321
    // With one round (`itr` is false after it), the PSSM is computed only with
    // `-save_pssm_after_last_round`, and `SavePssmToFile` writes the files of `-out_pssm` and
    // `-out_ascii_pssm` (`SaveCheckpoint`, `SaveAsciiPssm`).
    // With more rounds the driver writes the files (`rounds::RoundFiles`), not the search.
    pub(crate) fn request(&self) -> Option<PssmRequest> {
        // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:270-310
        // A PHI-BLAST round makes no PSSM (only a PSI round computes the next one), so one round
        // saves none.
        if self.rounds.is_some()
            || self.phi_pattern.is_some()
            || !self.save_pssm_after_last_round
            || (self.out_pssm.is_none() && self.out_ascii_pssm.is_none())
        {
            return None;
        }
        Some(PssmRequest {
            asn_path: self.out_pssm.clone(),
            ascii_path: self.out_ascii_pssm.clone(),
            save_each: self.save_each_pssm,
            pseudocount: self.pseudocount,
            versions: self.versions.clone(),
            run_files: self.run_files.clone(),
        })
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:246-248 (cited at
    // `in_pssm_files`: the loop writes the PSSM of each round at its start, the input PSSM of
    // `-in_pssm` in round 1)
    /// Whether the run can write a PSSM file: the PSSM after the last round (`request`), or one
    /// of the rounds' (`rounds`: at the start of each round after the first, and with `-in_pssm`
    /// before round 1).
    pub(crate) fn writes_files(&self) -> bool {
        self.request().is_some()
            || (self.rounds.is_some() && (self.out_pssm.is_some() || self.out_ascii_pssm.is_some()))
    }
}

/// The failure of the PSSM of the searched query `query` (an error that NCBI's psiblast raises
/// after that query's report, when it computes or saves its PSSM): the run reports the queries up
/// to it without the epilog, then fails with `error` (`blast_engine.rs` `search`).
#[derive(Debug)]
pub(crate) struct PssmStop {
    pub(crate) query: usize,
    pub(crate) error: anyhow::Error,
    /// The PSSM files of the queries before it, which NCBI has written: LOSAT writes them after
    /// the reports of the queries up to it (`PssmFiles::commit`).
    pub(crate) files: PssmFiles,
}

impl std::fmt::Display for PssmStop {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.error)
    }
}

impl std::error::Error for PssmStop {}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:216-220
pub(crate) fn blast_options_error(message: &str) -> anyhow::Error {
    crate::cli::NativeError {
        exit: 1,
        message: format!("{OPTIONS_ERROR} {message}\n").into(),
    }
    .into()
}

/// One Seq-align of a query's round-1 results, in the order of the Seq-align set: the subject
/// (`hsp.s_idx`), the e-value, the 0-based starts in the query and subject records and the
/// alignment's operations.
pub(crate) struct PssmAlignment<'a> {
    pub(crate) subject: u32,
    pub(crate) e_value: f64,
    pub(crate) q_start: usize,
    pub(crate) s_start: usize,
    pub(crate) q_length: usize,
    pub(crate) ops: Option<&'a [GapEditOp]>,
}

/// The E-value that the PSSM reads of an HSP (`GetLowestEvalue`: the lowest `e_value` or `sum_e`
/// score of the Seq-align, at most `BLAST_EXPECT_VALUE` 10, and 10 where it has none,
/// `seqalign_evalue_absent`); any other value of 1e-180 or less, a negative one of a database
/// with a negative size too, is 0 in it. The one owner of the inclusion test of a round
/// (`psi_round::included`) and of the PSSM's multiple alignment (`input_msa`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:512-545
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:158
pub(crate) fn seqalign_pssm_evalue(e_value: f64) -> f64 {
    const BLAST_EXPECT_VALUE: f64 = 10.0;
    if crate::report::serial::seqalign::seqalign_evalue_absent(e_value) {
        return BLAST_EXPECT_VALUE;
    }
    crate::report::serial::seqalign::seqalign_evalue(e_value).min(BLAST_EXPECT_VALUE)
}

impl<'a> PssmAlignment<'a> {
    /// The Seq-align of the BLASTP HSP `hsp` (1-based ends), whose query starts `q_shift` letters
    /// into the query's record (`-query_loc`).
    pub(crate) fn of_hsp(
        hsp: &'a crate::algorithm::blastp::hsp::BlastpHsp,
        q_shift: usize,
    ) -> Self {
        PssmAlignment {
            subject: hsp.s_idx,
            e_value: crate::report::serial::seqalign::seqalign_evalue(hsp.e_value),
            q_start: hsp.q_start - 1 + q_shift,
            s_start: hsp.s_start - 1,
            q_length: hsp.q_end + 1 - hsp.q_start,
            ops: hsp.gap_info.as_deref(),
        }
    }

    /// The fewest letters of a Bioseq from which NCBI's read of the Seq-align's subject range
    /// takes letters (`x_GetSubjectSequence`'s `CSeqVector` of that range, cited at `input_msa`):
    /// the end of the range; a shorter one throws, and the read gets none. The one rule of
    /// `input_msa` and of LOSAT's refusal of a read held by another sequence
    /// (`ListedTitles::pssm_holds` of BLASTP; review 5 of line DBID pfoff-4, H1).
    pub(crate) fn subject_end(&self) -> usize {
        self.s_start + subject_span(self)
    }
}

/// One query's Seq-aligns as the PSSM reads them ([`PssmAlignment`], in their order), kept apart
/// from the hit lists: psiblast makes a round's PSSMs after the round's reports, where NCBI reads
/// the subjects' letters for them (`ComputePssmForNextIteration`, cited at
/// `blast_engine::report_prep::report_values`), from the Seq-aligns before the formatter's cut,
/// which the reports' hit lists no longer hold by then. One flat list of the operations; the
/// call drops it once its PSSMs are made.
#[derive(Default)]
pub(crate) struct PssmSeqAligns {
    kept: Vec<KeptAlignment>,
    ops: Vec<GapEditOp>,
}

/// A [`PssmAlignment`] of [`PssmSeqAligns`]: its operations are `ops[start..end]` of the list
/// (`None` for an ungapped one).
struct KeptAlignment {
    subject: u32,
    e_value: f64,
    q_start: usize,
    s_start: usize,
    q_length: usize,
    ops: Option<(u32, u32)>,
}

impl PssmSeqAligns {
    /// The Seq-aligns `alignments`, in their order.
    pub(crate) fn of<'a>(alignments: impl IntoIterator<Item = PssmAlignment<'a>>) -> Self {
        let mut aligns = Self::default();
        for alignment in alignments {
            let ops = alignment.ops.map(|ops| {
                let start = aligns.ops.len() as u32;
                aligns.ops.extend_from_slice(ops);
                (start, aligns.ops.len() as u32)
            });
            aligns.kept.push(KeptAlignment {
                subject: alignment.subject,
                e_value: alignment.e_value,
                q_start: alignment.q_start,
                s_start: alignment.s_start,
                q_length: alignment.q_length,
                ops,
            });
        }
        aligns
    }

    /// The Seq-aligns, in their order.
    pub(crate) fn alignments(&self) -> Vec<PssmAlignment<'_>> {
        self.kept
            .iter()
            .map(|kept| PssmAlignment {
                subject: kept.subject,
                e_value: kept.e_value,
                q_start: kept.q_start,
                s_start: kept.s_start,
                q_length: kept.q_length,
                ops: kept
                    .ops
                    .map(|(start, end)| &self.ops[start as usize..end as usize]),
            })
            .collect()
    }

    /// Whether there is no Seq-align (a query without alignments makes no PSSM).
    pub(crate) fn is_empty(&self) -> bool {
        self.kept.is_empty()
    }
}

/// The Karlin-Altschul blocks of the round's ancillary data (`CBlastAncillaryData`, the blocks of
/// the query's context): `kbp_std`, `kbp_gap`, `kbp_psi` and `kbp_gap_psi`. In the first round
/// the PSI blocks are copies of the standard ones (`kbp_psi` and `kbp_gap_psi` of the search's
/// score block are its standard blocks, blast_setup.c:114-120); in a round that searched a PSSM
/// they are the PSSM's, and `kbp_gap` is `kbp_gap_psi` (blast_setup.c:113-126, cited at
/// `algorithm::pssm_search`).
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:99-111
#[derive(Debug, Clone, Copy)]
pub(crate) struct PssmAncillary {
    pub(crate) ungapped: KarlinBlk,
    pub(crate) gapped: KarlinBlk,
    pub(crate) psi_ungapped: KarlinBlk,
    pub(crate) psi_gapped: KarlinBlk,
}

/// The query of a PSSM: its letters (the whole record, without masks), its Seq-align id
/// (`fasta_query_id` of a FASTA query, or the first Seq-id of the PSSM query of `-in_pssm`) and
/// its title.
pub(crate) struct PssmQueryRecord<'a> {
    pub(crate) letters: &'a [u8],
    pub(crate) id: &'a ParsedSeqId,
    pub(crate) title: &'a [u8],
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:150-155
/// The Seq-id of the PSSM of FASTA query `record`, the `ordinal`-th of the run (1-based): its
/// Seq-align id, the local id `Query_<ordinal>`, or with `-parse_deflines` the id NCBI's Seq-aligns
/// take of those it read (`BelievedIds::best`: of `gi|4504347|ref|NP_000549.1|` the gi alone).
pub(crate) fn fasta_query_id(
    record: &crate::blastinput::fasta_reader::FastaRecord,
    ordinal: usize,
) -> ParsedSeqId {
    match record.believed() {
        Some(ids) => ids.best.clone(),
        None => ParsedSeqId::Local(ObjectId::Str(format!("Query_{ordinal}").into_bytes())),
    }
}

/// The scoring of the search that the PSSM records: the matrix as named (`-matrix`, any case)
/// and the gap costs.
#[derive(Debug, Clone, Copy)]
pub(crate) struct PssmScoring<'a> {
    pub(crate) matrix: ScoringMatrix,
    pub(crate) matrix_name: &'a str,
    pub(crate) gap_open: i32,
    pub(crate) gap_extend: i32,
}

/// The files of one query's PSSM: the ASN.1 text with `-out_pssm`, the ASCII PSSM with
/// `-out_ascii_pssm`.
#[derive(Debug)]
pub(crate) struct QueryPssm {
    /// The `PssmWithParameters` of the PSSM, made for its ASN.1 text (`-out_pssm`) and for
    /// psiblast's archive (outfmt 11), whose request of the round that searches the PSSM holds it
    /// (`PsiSearch::archive_pssm`); `None` when neither is written.
    pub(crate) object: Option<std::sync::Arc<PssmWithParameters>>,
    pub(crate) asn: Option<Vec<u8>>,
    /// NCBI's posts of the bytes that its stream replaces in the ASN.1 text (a query title
    /// outside the visible ASCII range; `report::serial::pssm::asn_text_posting`), which come
    /// with each write of the file (`file_posts`).
    pub(crate) asn_posts: Vec<ReplacedChar>,
    pub(crate) ascii: Option<Vec<u8>>,
    /// The engine's PSSM, which the next round searches (`algorithm::psiblast::rounds`).
    pub(crate) engine: PsiPssm,
}

/// The PSSM of one query (`ComputePssmForNextIteration`) and its files: the ASN.1 text
/// (`x_PSIMatrix2Asn1`, `PsiBlastAddAncillaryPssmData` and `MSerial_AsnText`) and the ASCII
/// PSSM (`PrintAsciiPssm`), as requested. The inputs are the query, its Seq-aligns, the letters
/// that NCBI's scope gives by their Seq-align ids (the subject of index `s_idx`: in NCBIstdaa, the
/// record's own or another sequence's, `None` for no sequence; `Err` for LOSAT's internal error of
/// a record without letters) and the keys of those ids (the key of index `s_idx`,
/// `rounds::SetKey`), the round's ancillary data, and the gapped block of the matrix and gap
/// costs (`gapped`, the score block's `kbp_gap_std[0]`). `unread` receives the subject of each
/// Seq-align whose letters cannot be read, in order (NCBI's warning `Failed to retrieve
/// sequence`, `input_msa`), also where the PSSM then fails. `archive`: the run's archive holds the
/// PSSM's object (`QueryPssm::object`).
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast.cpp:110-150
// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:380-404
// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:131-148
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1577-1600
// Every check of `PSIBlastOptionsValidate`, in its order (the options are never NULL here). The
// search before it does not check the inclusion threshold: nothing else calls this function, so
// a run with `-inclusion_ethresh 0` reports the query first. The ASCII PSSM prints the ancillary
// data's blocks; the ASN.1 PSSM keeps those of the scaling (`core::psi`).
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:894-903
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:98-99,166-170 and 166-171
// The MSA's query row and the PSSM's query Bioseq, which `PrintAsciiPssm` prints, are both the
// search's buffer (`GetBlastSequence` with `eBlastEncodingProtein`): an O is X in both. That is
// the only letter the buffer changes from the record's NCBIstdaa letters, and
// `encode_protein_sequence` makes the same change, so the letters are the buffer's (NCBI
// 2.17.0 and LOSAT give the same PSSMs for a query with lower-case letters, a lower-case o, and
// a `-` that the FASTA reader drops). The query has one title (the FASTA reader's), empty when
// the definition line has none, and `CPsiBlastInputData` writes it as the Bioseq's title also
// when it is empty.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn query_pssm(
    request: &PssmRequest,
    scoring: PssmScoring<'_>,
    gapped: KarlinBlk,
    inclusion_ethresh: f64,
    query: PssmQueryRecord<'_>,
    alignments: &[PssmAlignment<'_>],
    subject: &dyn Fn(u32) -> Result<Option<Vec<u8>>>,
    key: &dyn Fn(u32) -> SetKey,
    ancillary: PssmAncillary,
    unread: &mut Vec<u32>,
    archive: bool,
) -> Result<QueryPssm> {
    let mut sbp = score_block(scoring.matrix, gapped)?;
    sbp.set_ungapped_statistical_params(ancillary.psi_gapped);
    let encoded = crate::algorithm::blastp::encoding::encode_protein_sequence(query.letters);
    let letters = &encoded.aa_seq[1..=encoded.aa_len];
    let msa = input_msa(letters, alignments, subject, key, inclusion_ethresh, unread)?;
    if request.pseudocount < 0 {
        return Err(blast_options_error(
            "IPssmInputData returns invalid PSIBlastOptions: Pseudo count must be greater than or equal to 0",
        ));
    }
    if inclusion_ethresh <= 0.0 {
        return Err(blast_options_error(
            "IPssmInputData returns invalid PSIBlastOptions: Inclusion threshold must be greater than 0",
        ));
    }
    let pssm =
        psi::create_pssm_with_diagnostics(&msa, request.pseudocount, &sbp, request.diagnostics())
            .map_err(|error| crate::blastinput::app::engine_error(error_code_to_string(error)))?;
    let object = (request.asn_path.is_some() || archive).then(|| {
        std::sync::Arc::new(pssm_with_parameters(
            &pssm,
            letters,
            &query,
            scoring,
            request.pseudocount,
        ))
    });
    let (asn, asn_posts) = match (&request.asn_path, &object) {
        (Some(_), Some(value)) => {
            let (asn, posts) = crate::report::serial::pssm::asn_text_posting(value)
                .map_err(|error| anyhow::anyhow!("{error}"))?;
            (Some(asn), posts)
        }
        _ => (None, Vec::new()),
    };
    let ascii = request
        .ascii_path
        .as_ref()
        .map(|_| print_ascii_pssm(&pssm, letters, ancillary));
    Ok(QueryPssm {
        object,
        asn,
        asn_posts,
        ascii,
        engine: pssm,
    })
}

/// The `PssmWithParameters` of a PSSM (`x_PSIMatrix2Asn1` with its query, then
/// `PsiBlastAddAncillaryPssmData`).
// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:632-661
// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:684-778
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:426-433
// psiblast's options have no IMPALA scaling (`scalingFactor` is not set), and its diagnostics
// request never asks for the residue frequencies or the independent observations. The interval
// sizes and numbers of matching sequences are `Uint4` pushed to a `vector<int>`.
fn pssm_with_parameters(
    pssm: &PsiPssm,
    letters: &[u8],
    query: &PssmQueryRecord<'_>,
    scoring: PssmScoring<'_>,
    pseudocount: i32,
) -> PssmWithParameters {
    let entry = SeqEntry::Seq(Box::new(PssmQuery {
        ids: vec![query.id.clone()],
        descr: Some(vec![Seqdesc::Title(query.title.to_vec())]),
        inst: SeqInst::new(
            ("raw", 2),
            ("aa", 3),
            Some(letters.len() as u32),
            Some(SeqData::Ncbistdaa(letters.to_vec())),
        ),
        annot: None,
    }));
    PssmWithParameters {
        pssm: pssm_object(pssm, letters.len(), Some(entry)),
        params: Some(PssmParameters {
            pseudocount: Some(pseudocount),
            rpsdbparams: Some(FormatRpsDbParameters {
                matrix_name: scoring.matrix_name.to_ascii_uppercase().into_bytes(),
                gap_open: Some(scoring.gap_open),
                gap_extend: Some(scoring.gap_extend),
            }),
            ..PssmParameters::default()
        }),
    }
}

/// The `PssmWithParameters` that `CPssmEngine::Run` returns for an alignment (`-in_msa`): the
/// engine's `x_PSIMatrix2Asn1` with the options' diagnostics and the query Bioseq (`entry`, none
/// when the alignment holds no Bioseq of the query's letters), whose parameters are the
/// pseudocount of the default PSI-BLAST options (0) and the matrix's name in capitals; the gap
/// costs are `PsiBlastAddAncillaryPssmData`'s, which this path does not call.
pub(crate) fn msa_pssm_with_parameters(
    pssm: &PsiPssm,
    columns: usize,
    entry: Option<SeqEntry>,
    matrix_name: &str,
) -> PssmWithParameters {
    PssmWithParameters {
        pssm: pssm_object(pssm, columns, entry),
        params: Some(PssmParameters {
            pseudocount: Some(0),
            rpsdbparams: Some(FormatRpsDbParameters {
                matrix_name: matrix_name.to_ascii_uppercase().into_bytes(),
                gap_open: None,
                gap_extend: None,
            }),
            ..PssmParameters::default()
        }),
    }
}

/// `x_PSIMatrix2Asn1`'s `CPssm`: the scores, the Karlin-Altschul values and the diagnostics of
/// `pssm`, with the query `entry`.
fn pssm_object(pssm: &PsiPssm, columns: usize, entry: Option<SeqEntry>) -> Pssm {
    let by_column =
        |rows: &[[f64; BLASTAA_SIZE]]| -> Vec<f64> { rows.iter().flatten().copied().collect() };
    let diagnostics = &pssm.diagnostics;
    let intermediate_data = PssmIntermediateData {
        res_freqs_per_pos: None,
        weighted_res_freqs_per_pos: diagnostics.weighted_residue_freqs.as_deref().map(by_column),
        freq_ratios: diagnostics
            .frequency_ratios
            .as_deref()
            .map(by_column)
            .unwrap_or_default(),
        information_content: diagnostics.information_content.clone(),
        gapless_column_weights: diagnostics.gapless_column_weights.clone(),
        sigma: diagnostics.sigma.clone(),
        interval_sizes: diagnostics
            .interval_sizes
            .as_ref()
            .map(|sizes| sizes.iter().map(|&size| size as i32).collect()),
        num_matching_seqs: diagnostics
            .num_matching_seqs
            .as_ref()
            .map(|counts| counts.iter().map(|&count| count as i32).collect()),
        num_indept_obsr: None,
    };
    Pssm {
        is_protein: Some(true),
        identifier: None,
        num_rows: BLASTAA_SIZE as i32,
        num_columns: columns as i32,
        row_labels: None,
        by_row: Some(false),
        query: entry,
        intermediate_data: Some(intermediate_data),
        final_data: Some(PssmFinalData {
            scores: pssm.scores.iter().flatten().copied().collect(),
            lambda: pssm.gapped.lambda,
            kappa: pssm.gapped.k,
            h: pssm.gapped.h,
            scaling_factor: None,
            lambda_ungapped: Some(pssm.ungapped.lambda),
            kappa_ungapped: Some(pssm.ungapped.k),
            h_ungapped: Some(pssm.ungapped.h),
            word_score_threshold: None,
        }),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:218-282
pub(crate) fn error_code_to_string(error: PsiError) -> &'static str {
    match error {
        PsiError::BadSeqWeights => "Error computing sequence weights",
        PsiError::PositiveAvgScore => "PSSM has positive average score",
        PsiError::NoAlignedSeqs => {
            "No sequences left after purging biased sequences in multiple sequence alignment"
        }
        PsiError::GapInQuery => "Gap found in query sequence",
        PsiError::UnalignedColumn => "Found column with no sequences aligned in it",
        PsiError::ColumnOfGaps => "Found column with only GAP residues",
        PsiError::StartingGap => "Found flanking gap at start of alignment",
        PsiError::EndingGap => "Found flanking gap at end of alignment",
        // `x_ErrorCodeToString` of `PSIERR_UNKNOWN`, which is -255 (blast_psi_priv.h:392).
        PsiError::Unknown => "Unknown error code returned from PSSM engine: -255",
        PsiError::NoMatrixFrequencyData => {
            "No matrix frequency ratios were found for requested matrix"
        }
        // NCBI crashes here; `-in_msa` refuses the run before any text is written.
        PsiError::NoUngappedKarlinBlock => {
            "The ungapped Karlin-Altschul block of the query is not available"
        }
        // `x_ErrorCodeToString(1)`: `_IMPALAScaleMatrix`'s status (makeprofiledb).
        PsiError::ImpalaScoreRange => "Unknown error code returned from PSSM engine: 1",
        // makeprofiledb refuses the PSSM before this text is written.
        PsiError::ScalingUndefined => "The scaling of the PSSM leaves defined behaviour",
    }
}

/// The engine's input from the query and its Seq-aligns (`CPsiBlastInputData::Process`).
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:122-141
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:529-565
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:283-323
// The rows are the handles of the counted Seq-align ids (`key`: `rounds::SetKey`, which is the
// handle in the runs LOSAT accepts, `SetKey::check_handles`), and a Seq-align whose id `Match`es
// the previous Seq-align's id takes the previous row (`SetKey::matches`). So copies of one record
// that follow each other share a row (their later HSPs fill only cells left unaligned), and so do
// ids of two handles that `Match` (`gb|X|` and `gb|X.1|`), which leave the last row empty (line
// DBID weakbatch, review 2 H1 and DECISIONS 43 (4)).
// `msa_index` also moves on after a subject that is not counted, so an HSP below the threshold
// whose subject comes after one without such an HSP, before the first counted subject, can take
// `msa_index` past the rows NCBI allocates (one per counted handle), as can a copy of a record
// after another subject (one handle, two rows; NCBI BLAST+ 2.17.0 crashes, review 2 P3). NCBI then
// writes out of bounds (`_ASSERT` only), a failure without a checkable result: LOSAT refuses the
// PSSM there.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:343-354 and 440-450
// Each Seq-align reads its subject's range by the Seq-align id: the letters of the Bioseq that
// NCBI's scope gives for it (`subject`, the subject's own or another sequence's), and none where
// it gives no Bioseq or the range ends beyond that one's letters (probe of line DBID pssmread:
// one letter short fails, the same length and longer read that sequence's letters). Then the
// row becomes the query's letters, also over the cells that earlier Seq-aligns of the row filled,
// which the purge removes, and the subject goes to `unread` for the warning.
// The hit lists are sorted by their best e-value (`s_EvalueCompareHSPLists`), which puts the
// counted subjects first, but two e-values below 1.0e-180 compare as equal and those lists are
// then ordered by score, so with an inclusion threshold below 1.0e-180 the counted subjects are
// not shown to come first:
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1390-1395,3089-3101 and 3095-3102
fn input_msa(
    query: &[u8],
    alignments: &[PssmAlignment<'_>],
    subject: &dyn Fn(u32) -> Result<Option<Vec<u8>>>,
    key: &dyn Fn(u32) -> SetKey,
    inclusion_ethresh: f64,
    unread: &mut Vec<u32>,
) -> Result<PsiMsa> {
    let qualifying: std::collections::BTreeSet<SetKey> = alignments
        .iter()
        .filter(|alignment| seqalign_pssm_evalue(alignment.e_value) < inclusion_ethresh)
        .map(|alignment| key(alignment.subject))
        .collect();
    let mut msa = PsiMsa::new(query.len(), qualifying.len());
    for (cell, &letter) in msa.rows[0].iter_mut().zip(query) {
        *cell = PsiMsaCell {
            letter,
            is_aligned: true,
        };
    }
    let mut msa_index = 1;
    let mut last_sid: Option<SetKey> = None;
    for alignment in alignments {
        let current_sid = key(alignment.subject);
        if seqalign_pssm_evalue(alignment.e_value) < inclusion_ethresh {
            if last_sid
                .as_ref()
                .is_some_and(|last| !current_sid.matches(last))
            {
                msa_index += 1;
            }
            let Some(row) = msa.rows.get_mut(msa_index) else {
                anyhow::bail!(
                    "an alignment past the last row of the PSSM's multiple alignment (a subject included by -inclusion_ethresh after one that is not, or a copy of a record after another subject, where NCBI BLAST+'s psiblast writes out of bounds) is not supported by LOSAT's PSIBLAST"
                );
            };
            match subject(alignment.subject)?
                .filter(|letters| letters.len() >= alignment.subject_end())
            {
                Some(letters) => process_denseg(row, alignment, &letters),
                None => {
                    for (cell, &letter) in row.iter_mut().zip(query) {
                        *cell = PsiMsaCell {
                            letter,
                            is_aligned: true,
                        };
                    }
                    unread.push(alignment.subject);
                }
            }
        }
        last_sid = Some(current_sid);
    }
    Ok(msa)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:363-408
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:425-447
// The subject's letters from the first aligned one on: an HSP's operations walk them from
// `s_start` (an HSP starts with aligned letters). An HSP without operations is ungapped.
fn process_denseg(row: &mut [PsiMsaCell], alignment: &PssmAlignment<'_>, subject: &[u8]) {
    let whole = [GapEditOp::Sub(alignment.q_length as u32)];
    let ops = match alignment.ops {
        Some(ops) if !ops.is_empty() => ops,
        _ => &whole[..],
    };
    let mut q = alignment.q_start;
    let mut s = alignment.s_start;
    for op in ops {
        match *op {
            GapEditOp::Del(n) => s += n as usize,
            GapEditOp::Ins(n) => {
                for cell in &mut row[q..q + n as usize] {
                    if !cell.is_aligned {
                        *cell = PsiMsaCell {
                            letter: ncbistdaa::GAP,
                            is_aligned: true,
                        };
                    }
                }
                q += n as usize;
            }
            GapEditOp::Sub(n) => {
                for (cell, &letter) in row[q..q + n as usize]
                    .iter_mut()
                    .zip(&subject[s..s + n as usize])
                {
                    if !cell.is_aligned {
                        *cell = PsiMsaCell {
                            letter,
                            is_aligned: true,
                        };
                    }
                }
                q += n as usize;
                s += n as usize;
            }
        }
    }
}

/// The subject letters that an HSP aligns, from its first: the length of the range that
/// `x_GetSubjectSequence` reads (the segments without a gap in the subject, cited at
/// `process_denseg`).
fn subject_span(alignment: &PssmAlignment<'_>) -> usize {
    match alignment.ops {
        Some(ops) if !ops.is_empty() => ops
            .iter()
            .map(|op| match *op {
                GapEditOp::Sub(n) | GapEditOp::Del(n) => n as usize,
                GapEditOp::Ins(_) => 0,
            })
            .sum(),
        _ => alignment.q_length,
    }
}

/// A subject record's letters in NCBIstdaa as its Bioseq holds them (the `CSeqVector` of
/// `CPsiBlastInputData::x_GetSubjectSequence`, cited at `process_denseg`): upper case, and O
/// stays O. The query's letters are not these: the PSSM's query is the search's buffer, in which
/// O is X (`query_pssm`).
pub(crate) fn ncbistdaa_letters(letters: &[u8]) -> Vec<u8> {
    letters
        .iter()
        .map(|&letter| aa_char_to_ncbistdaa(letter.to_ascii_uppercase()))
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:145-166
pub(super) const RESIDUE_ORDER: [usize; 20] = [
    1, 16, 13, 4, 3, 15, 5, 7, 8, 9, 11, 10, 12, 6, 14, 17, 18, 20, 22, 19,
];

/// A double as the C++ stream writes it with `fixed` and `setprecision` of `precision` (the C
/// library's `%.*f`: correctly rounded, as Rust's formatting; `inf` and `nan` as glibc spells
/// them).
fn fixed(value: f64, precision: usize) -> String {
    if value.is_nan() {
        return if value.is_sign_negative() {
            "-nan"
        } else {
            "nan"
        }
        .to_string();
    }
    if value.is_infinite() {
        return if value < 0.0 { "-inf" } else { "inf" }.to_string();
    }
    format!("{value:.precision$}")
}

/// `CBlastFormatUtil::PrintAsciiPssm` of a computed PSSM.
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:191-327
// The scores matrix of `GetScores` has the 28 rows of NCBIstdaa (pssm_engine.cpp
// `x_PSIMatrix2Asn1` stores every row). The query's letters are those of the PSSM's query Bioseq
// (`GetQuerySequenceData`, blastfmtutil.cpp:253-255), which `CPsiBlastInputData` builds from the
// search's buffer, in which an O is X (cited at `query_pssm`): `query` is that buffer, so an O
// of the record prints as X and takes the X branch of the gapless weight. The diagnostics it
// reads are those that `-out_ascii_pssm` asks for (`PssmRequest::diagnostics`).
pub(crate) fn print_ascii_pssm(pssm: &PsiPssm, query: &[u8], ancillary: PssmAncillary) -> Vec<u8> {
    const REQUESTED: &str = "a diagnostic that -out_ascii_pssm requests";
    let diagnostics = &pssm.diagnostics;
    // `vector<int>` of `GetNumMatchingSeqs`: the engine's `Uint4` counts as `int`.
    let num_matching_seqs: Vec<i32> = diagnostics
        .num_matching_seqs
        .as_ref()
        .expect(REQUESTED)
        .iter()
        .map(|&count| count as i32)
        .collect();
    write_ascii_pssm(
        &pssm.scores,
        query,
        Some(AsciiPssmDiagnostics {
            weighted_residue_freqs: diagnostics
                .weighted_residue_freqs
                .as_ref()
                .expect(REQUESTED),
            information_content: diagnostics.information_content.as_ref().expect(REQUESTED),
            num_matching_seqs: &num_matching_seqs,
            gapless_column_weights: diagnostics
                .gapless_column_weights
                .as_ref()
                .expect(REQUESTED),
        }),
        &ancillary_blocks(&ancillary),
    )
}

/// The blocks that the footer of a computed PSSM prints: each of the round's ancillary blocks that
/// the results keep (`stats::tables::ancillary_kept`, for all four). The input PSSM's footer of
/// round 1 (`in_pssm_files`) prints its two PSI blocks as the file holds them: NCBI makes that
/// ancillary data from the PSSM's values without the check (`ExtractPssmAncillaryData`).
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:296-326
fn ancillary_blocks(ancillary: &PssmAncillary) -> Vec<(&'static str, KarlinBlk)> {
    [
        ("Standard Ungapped    ", ancillary.ungapped),
        ("Standard Gapped      ", ancillary.gapped),
        ("PSI Ungapped         ", ancillary.psi_ungapped),
        ("PSI Gapped           ", ancillary.psi_gapped),
    ]
    .into_iter()
    .filter(|(_, block)| crate::stats::tables::ancillary_kept(block.lambda))
    .collect()
}

/// The columns that the long form of `PrintAsciiPssm` prints (a PSSM with its information
/// content): `[position][residue]` weighted residue frequencies, the information content, the
/// numbers of matching sequences and the gapless column weights of each position.
pub(crate) struct AsciiPssmDiagnostics<'a> {
    pub(crate) weighted_residue_freqs: &'a [[f64; BLASTAA_SIZE]],
    pub(crate) information_content: &'a [f64],
    pub(crate) num_matching_seqs: &'a [i32],
    pub(crate) gapless_column_weights: &'a [f64],
}

/// `CBlastFormatUtil::PrintAsciiPssm` (cited at `print_ascii_pssm`) of the scores
/// `scores[position][residue]` of the query `query` (NCBIstdaa codes): the long form with
/// `long`, the short one (`Last position-specific scoring matrix computed` and the scores alone)
/// without, then the Karlin-Altschul blocks of the ancillary data that it has, in its order.
pub(crate) fn write_ascii_pssm(
    scores: &[[i32; BLASTAA_SIZE]],
    query: &[u8],
    long: Option<AsciiPssmDiagnostics<'_>>,
    blocks: &[(&str, KarlinBlk)],
) -> Vec<u8> {
    use std::fmt::Write as _;
    let mut out = String::new();
    if long.is_some() {
        out.push_str("\nLast position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts\n");
    } else {
        out.push_str("\nLast position-specific scoring matrix computed\n");
    }
    let mut max_score: i32 = 0;
    for &score in scores.iter().flatten() {
        if score <= BLAST_SCORE_MIN {
            continue;
        }
        if score > max_score {
            max_score = score;
        }
        if -score > max_score {
            max_score = -score;
        }
    }
    let mut num_digits = 0;
    while max_score > 0 {
        max_score /= 10;
        num_digits += 1;
    }
    let width = num_digits + 2;
    out.push_str("         ");
    for &c in &RESIDUE_ORDER {
        let _ = write!(out, "{:>width$}", NCBISTDAA_TO_AMINOACID[c] as char);
    }
    if long.is_some() {
        for &c in &RESIDUE_ORDER {
            let _ = write!(out, "   {}", NCBISTDAA_TO_AMINOACID[c] as char);
        }
    }
    for (i, &residue) in query.iter().enumerate() {
        let row = &scores[i];
        let _ = write!(
            out,
            "\n{:>5} {}  ",
            i + 1,
            NCBISTDAA_TO_AMINOACID[residue as usize] as char
        );
        for &c in &RESIDUE_ORDER {
            if row[c] == BLAST_SCORE_MIN {
                out.push_str("-I ");
            } else {
                let _ = write!(out, "{:>width$}", row[c]);
            }
        }
        out.push(' ');
        let Some(long) = &long else {
            continue;
        };
        for &c in &RESIDUE_ORDER {
            if row[c] != BLAST_SCORE_MIN {
                let mut value = 100.0;
                value *= long.weighted_residue_freqs[i][c];
                value = (value + if value > 0.0 { 0.5 } else { -0.5 }) as i32 as f64;
                let _ = write!(out, "{:>4}", value as i32);
            }
        }
        let _ = write!(out, "  {} ", fixed(long.information_content[i], 2));
        if long.num_matching_seqs[i] > 1 && residue != ncbistdaa::X {
            out.push_str(&fixed(long.gapless_column_weights[i], 2));
        } else {
            out.push_str("    0.00");
        }
    }
    out.push_str("\n\n");
    out.push_str("                      K         Lambda\n");
    for (label, block) in blocks {
        let _ = writeln!(
            out,
            "{label}{}     {}",
            fixed(block.k, 4),
            fixed(block.lambda, 4)
        );
    }
    out.into_bytes()
}

/// NCBI's posts of the bytes replaced in a PSSM's ASN.1 text (`QueryPssm::asn_posts`), which
/// psiblast makes while `SavePssmToFile` writes the file, one error post 803.7 of the application
/// `application` for each byte, after what the run wrote before (`cerr` is tied to `cout`), each
/// whole, which a BLAST archive saves as one message (`Posts::write`; NCBI 2.17.0 `psiblast
/// -outfmt 11 -out_pssm` saves them in the next archive's messages and the last ones in the error
/// archive at the end of the run; LOSAT's psiblast does not write outfmt 11). One title byte
/// outside 0x20..0x7E gives `Error: (803.7) [psiblast]
/// PssmWithParameters.pssm.query.seq.descr.E.title` / `Bad char [0xE9] in string at line 12` /
/// the whole title, `#` in its place in the file, status 0.
// NCBI reference (598d8ae6): c++/src/serial/objostrasn.cpp:413-444 (cited at
// `report::serial::seqalign::write_query_posting`)
// NCBI reference (598d8ae6): c++/src/serial/objistr.cpp:1855-1897 (cited at
// `report::serial::encoder::ReplacedChar`)
pub(crate) fn file_posts(application: &str, replaced: &[ReplacedChar]) -> Posts {
    let mut posts = Posts::new();
    for byte in replaced {
        posts.post(crate::report::diag_log::ERROR, byte.post(Some(application)));
    }
    posts
}

/// `CAutoOutputFileReset`: the file of one PSSM option, opened again at each write
/// (`GetStream`). Without `-save_each_pssm` an existing file of the name is deleted and the file
/// made again, so the last PSSM stays; with it the writes are `F.1`, `F.2`, ...
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:50-66
// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:54-55
#[derive(Debug)]
pub(crate) struct OutputFileReset {
    file_name: std::path::PathBuf,
    /// `m_Version`: `None` for version 0 (no `-save_each_pssm`), else the run's number of the
    /// stream (`PssmVersions`).
    version: Option<Version>,
}

impl OutputFileReset {
    fn new(file_name: &std::path::Path, version: Option<Version>) -> Self {
        Self {
            file_name: file_name.to_path_buf(),
            version,
        }
    }

    /// `GetStream`'s file: `F`, or `F.<version>` (the version counts the run's writes).
    pub(crate) fn next_file(&mut self) -> PssmFile {
        let Some(version) = self.version.as_ref() else {
            return PssmFile {
                path: self.file_name.clone(),
                replaced: true,
            };
        };
        let number = version.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        PssmFile {
            path: numbered_pssm_path(&self.file_name, number),
            replaced: false,
        }
    }
}

fn numbered_pssm_path(base: &std::path::Path, number: u32) -> std::path::PathBuf {
    let mut name = base.as_os_str().to_os_string();
    name.push(format!(".{number}"));
    std::path::PathBuf::from(name)
}

/// One file of `CAutoOutputFileReset::GetStream`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct PssmFile {
    pub(crate) path: std::path::PathBuf,
    /// Version 0 (no `-save_each_pssm`): an existing regular file of the name is removed before
    /// the file is opened, so the file is made again; a numbered file is only truncated, which
    /// keeps an existing file's inode and mode.
    pub(crate) replaced: bool,
}

impl PssmFile {
    /// `GetStream`: with version 0, `CFile::Exists` (a regular file, through symbolic links) and
    /// `Remove` (the name itself, the C library's `remove`; a failure is neither checked nor, by
    /// default, logged), then `ofstream` (made or truncated, through symbolic links).
    // NCBI reference (598d8ae6): c++/include/corelib/ncbifile.hpp:4039-4042
    // NCBI reference (598d8ae6): c++/src/corelib/ncbifile.cpp:2595-2610
    pub(crate) fn open(&self) -> std::io::Result<std::fs::File> {
        if self.replaced && self.path.is_file() {
            let _ = std::fs::remove_file(&self.path);
        }
        std::fs::File::create(&self.path)
    }
}

/// What a PSSM file's name is before the run, for NCBI's `GetStream` on it (`PssmFile::open`).
#[derive(Debug)]
enum TargetKind {
    /// No file: `ofstream` makes it in its directory.
    Absent,
    /// A regular file, also through a symbolic link.
    Regular,
    /// A character device (`/dev/null`, `/dev/full`): opened for writing as it is; `null`: it is
    /// `/dev/null`.
    #[cfg_attr(
        not(unix),
        expect(dead_code, reason = "constructed only on unix (special_target_kind)")
    )]
    Device { null: bool },
    /// A directory or a socket: `ofstream` cannot open it.
    Unopenable,
    /// A named pipe, a block device or a symbolic link to no file: not supported.
    Unsupported(&'static str),
    /// A name whose directories cannot be read: `ofstream` fails the same way.
    Unreachable(std::io::Error),
}

fn target_kind(path: &std::path::Path) -> TargetKind {
    match std::fs::metadata(path) {
        Ok(metadata) if metadata.is_file() => TargetKind::Regular,
        Ok(metadata) => special_target_kind(&metadata),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            if std::fs::symlink_metadata(path).is_ok() {
                TargetKind::Unsupported("a symbolic link to a missing file")
            } else {
                TargetKind::Absent
            }
        }
        Err(error) => TargetKind::Unreachable(error),
    }
}

#[cfg(unix)]
fn special_target_kind(metadata: &std::fs::Metadata) -> TargetKind {
    use std::os::unix::fs::FileTypeExt;
    let file_type = metadata.file_type();
    if file_type.is_fifo() {
        TargetKind::Unsupported("a named pipe")
    } else if file_type.is_block_device() {
        TargetKind::Unsupported("a block device")
    } else if file_type.is_char_device() {
        TargetKind::Device {
            null: device_number(metadata).is_some_and(|rdev| Some(rdev) == null_device()),
        }
    } else {
        TargetKind::Unopenable
    }
}

#[cfg(not(unix))]
fn special_target_kind(_metadata: &std::fs::Metadata) -> TargetKind {
    TargetKind::Unopenable
}

/// The device number of a character device.
#[cfg(unix)]
fn device_number(metadata: &std::fs::Metadata) -> Option<u64> {
    use std::os::unix::fs::{FileTypeExt, MetadataExt};
    metadata
        .file_type()
        .is_char_device()
        .then(|| metadata.rdev())
}

/// `/dev/null`'s device number.
#[cfg(unix)]
fn null_device() -> Option<u64> {
    std::fs::metadata("/dev/null")
        .ok()
        .as_ref()
        .and_then(device_number)
}

/// Whether the test open of an existing file failed with `EPERM`: an immutable or append-only file
/// (`open(2)`), whose name cannot be removed either (`unlink(2)`, `EPERM`).
#[cfg(unix)]
fn immutable(error: &std::io::Error) -> bool {
    rustix::io::Errno::from_io_error(error) == Some(rustix::io::Errno::PERM)
}

#[cfg(not(unix))]
fn immutable(_error: &std::io::Error) -> bool {
    false
}

/// Whether removing `path` succeeds, given that its directory takes a new file (the staging file
/// made there by the user `uid`): in a directory with the sticky bit only the owner of the file or
/// of the directory, or root, removes a name (`unlink(2)`, `EPERM`).
#[cfg(unix)]
fn removable(path: &std::path::Path, directory: &std::path::Path, uid: u32) -> bool {
    use std::os::unix::fs::MetadataExt;
    let (Ok(entry), Ok(directory)) = (
        std::fs::symlink_metadata(path),
        std::fs::metadata(directory),
    ) else {
        return false;
    };
    directory.mode() & 0o1000 == 0 || uid == 0 || entry.uid() == uid || directory.uid() == uid
}

/// Elsewhere a file that cannot be opened for writing is taken as one that cannot be removed
/// either (a refusal where NCBI might remove it).
#[cfg(not(unix))]
fn removable(_path: &std::path::Path, _directory: &std::path::Path, _uid: u32) -> bool {
    false
}

#[cfg(unix)]
fn owner(metadata: &std::fs::Metadata) -> u32 {
    std::os::unix::fs::MetadataExt::uid(metadata)
}

#[cfg(not(unix))]
fn owner(_metadata: &std::fs::Metadata) -> u32 {
    0
}

/// The PSSM files of a psiblast run (`-out_pssm`, `-out_ascii_pssm`): made with the PSSMs, before
/// the reports (`add`), and written after them (`commit`), each query's ASN.1 file then its ASCII
/// one, as psiblast's `SavePssmToFile` after the query's report. The ASCII PSSMs, whose writes
/// NCBI does not check, are kept in memory.
///
/// NCBI ends a run whose ASN.1 file cannot be opened or written after that query's report, so LOSAT
/// must know before its first report byte whether each file can be written; and a run that fails
/// before its reports are out (a report write error, a LOSAT refusal) must leave every file as it
/// was, as NCBI's, which stops before `SavePssmToFile`. So `add` writes the ASN.1 bytes to a
/// staging file beside the target (`.losat-<pid>-<n>`, a new file in the target's directory, whose
/// name is short so that any name the target can have is no reason for its failure), which shows
/// that the directory takes a file and the file system the bytes, and foresees NCBI's open of the
/// target without changing it; only `commit` writes the targets.
/// Dropped without `commit`, `PssmFiles` removes the staging files: no target has changed. A
/// failure in `add` is the refusal of a file that NCBI cannot write either (its stop after that
/// query's report, `PssmStop`) or that LOSAT cannot write as NCBI does; no partial file is left.
///
/// The rule for a target, from what `GetStream` does with it (`PssmFile::open`):
/// - no file of the name: made in its directory, which takes the staging file;
/// - a regular file (also through a symbolic link): removed and made again (version 0) or
///   truncated in place. LOSAT opens it for writing (without truncation); if that fails with
///   `EPERM` (an immutable or append-only file, whose name cannot be removed either), or fails and
///   the name cannot be removed (`removable`) or the file is numbered, NCBI's open fails there.
///   An existing file that opens for writing but beside which no staging file can be made (a
///   directory that takes no new file) is written in place, as NCBI's version 0, which cannot remove
///   its name, and a numbered file do: its bytes are kept in memory for `commit` (only the last
///   query's with version 0), which cannot show before the reports that they fit;
/// - a character device: opened for writing and written zero bytes, which a device that takes no
///   bytes refuses (`/dev/full`: ENOSPC), so the run stops at that query as NCBI's. `/dev/null`
///   gets nothing more; the PSSM of every query for another device is kept in memory (a few MB per
///   query of some thousand residues) and written by `commit`. A device that takes zero bytes but
///   not the PSSM's would fail only after the reports;
/// - a directory or a socket: NCBI's open fails, and so does LOSAT's;
/// - a named pipe, a block device, a symbolic link to no file: refused (also before the search,
///   `check_pssm_files`).
///
/// Not foreseen: a file system that fills between `add` and `commit`, or that has no room for the
/// bytes of a file written in place (LOSAT's peak use is the old file, at most two staged PSSMs of
/// version 0 and every numbered one); such a failure comes after the reports (`commit`).
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:385-407
#[derive(Debug)]
pub(crate) struct PssmFiles {
    asn: Option<OutputFileReset>,
    ascii: Option<OutputFileReset>,
    run_files: crate::report::run_files::RunFiles,
    pending: Vec<PendingPssm>,
    /// The entry whose staging file holds the last version-0 ASN.1 PSSM, which the next one
    /// replaces: of NCBI's writes of one name only the last one's bytes stay, made or truncated
    /// as each of them would make or truncate it.
    last_replaced: Option<usize>,
    staging_count: u64,
    /// NCBI's messages of each query's PSSM reads (`PssmReadMessages`), by searched query in the
    /// order of the queries: NCBI writes them after that query's report.
    read_warnings: Vec<(usize, PssmReadMessages)>,
    /// NCBI's posts of the bytes replaced in each query's ASN.1 PSSM (`file_posts`), by searched
    /// query in the order of the queries: NCBI writes them when it writes the file, after the
    /// query's read messages.
    file_posts: Vec<(usize, Posts)>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:440-450 (cited at
// `input_msa`)
/// What NCBI writes on standard error while one query's PSSM reads its subjects by their
/// Seq-align ids (`x_GetSubjectSequence`, after the query's report), in the order of the reads,
/// each line with the subject (by subject index) whose read writes it: the warning of each
/// Seq-align whose letters the read cannot get (`Failed to retrieve sequence`, from `unread` of
/// `query_pssm`), and the taxonomy database's texts where a read loads a sequence whose taxids
/// the volume looks up then (`add_load_texts`). The scope loads a subject's sequence at the
/// subject's first Seq-align (`CSeqVector`) before the letters are read, so the subject's texts
/// come before its warnings and after those of the subjects read before it.
#[derive(Debug, Default)]
pub(crate) struct PssmReadMessages {
    lines: Vec<(u32, Vec<u8>)>,
}

impl PssmReadMessages {
    /// One line (or text) that the read of subject `subject` writes, after those before it.
    pub(crate) fn push(&mut self, subject: u32, line: Vec<u8>) {
        self.lines.push((subject, line));
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.lines.is_empty()
    }

    /// The bytes of the messages, in their order.
    pub(crate) fn bytes(&self) -> Vec<u8> {
        self.lines
            .iter()
            .flat_map(|(_, line)| line.iter().copied())
            .collect()
    }

    /// The taxonomy texts of the PSSM's loads merged in: `loads` holds every subject that the
    /// PSSM includes, in the order of its reads (each once, at its first included Seq-align), with
    /// the texts of that read's load (empty for none). A subject's texts come before the warnings
    /// of its own Seq-aligns, which follow those of the subjects before it (`push` gives every
    /// warning in the order of the Seq-aligns, which keep a subject's together).
    pub(crate) fn add_load_texts(&mut self, loads: impl IntoIterator<Item = (u32, Vec<u8>)>) {
        let mut warnings = std::mem::take(&mut self.lines).into_iter().peekable();
        for (subject, text) in loads {
            if !text.is_empty() {
                self.lines.push((subject, text));
            }
            while let Some(line) = warnings.next_if(|(read, _)| *read == subject) {
                self.lines.push(line);
            }
        }
        debug_assert!(
            warnings.peek().is_none(),
            "every warning's subject is one of the PSSM's reads"
        );
        self.lines.extend(warnings);
    }
}

/// One query's PSSM files, to be written after the reports.
#[derive(Debug)]
struct PendingPssm {
    asn: Option<(PssmFile, Staged)>,
    ascii: Option<Vec<u8>>,
}

#[derive(Debug)]
enum Staged {
    /// The bytes in a staging file beside the target.
    File(std::path::PathBuf),
    /// The bytes of an existing file in a directory that takes no new file, written in place.
    InPlace(Vec<u8>),
    /// The bytes for a character device.
    Device(Vec<u8>),
    /// The host file name and bytes for `RunFiles::Receiver` (LOSAT Web).
    Receiver(Vec<u8>, Vec<u8>),
    /// Nothing: the target is `/dev/null`.
    Null,
    /// Written, or replaced by a later query's file of the same name.
    Done,
}

impl PssmFiles {
    pub(crate) fn new(request: &PssmRequest) -> Self {
        let [asn, ascii] = request.streams();
        Self {
            asn,
            ascii,
            run_files: request.run_files.clone(),
            pending: Vec::new(),
            last_replaced: None,
            staging_count: 0,
            read_warnings: Vec::new(),
            file_posts: Vec::new(),
        }
    }

    /// The warnings of the reads of the searched query `query`'s PSSM (`read_warnings`), made in
    /// the order of the queries.
    pub(crate) fn add_read_warnings(&mut self, query: usize, warnings: PssmReadMessages) {
        if !warnings.is_empty() {
            self.read_warnings.push((query, warnings));
        }
    }

    /// The messages of the reads of the searched query `query`'s PSSM, to add the taxonomy
    /// texts of its loads to (`PssmReadMessages::add_load_texts`); an empty entry is made in the
    /// order of the queries where it has no warning.
    pub(crate) fn read_messages(&mut self, query: usize) -> &mut PssmReadMessages {
        let at = match self
            .read_warnings
            .binary_search_by_key(&query, |&(of, _)| of)
        {
            Ok(at) => at,
            Err(at) => {
                self.read_warnings
                    .insert(at, (query, PssmReadMessages::default()));
                at
            }
        };
        &mut self.read_warnings[at].1
    }

    /// The posts of the write of the searched query `query`'s ASN.1 PSSM (`file_posts`), made in
    /// the order of the queries.
    pub(crate) fn add_file_posts(&mut self, query: usize, posts: Posts) {
        if !posts.is_empty() {
            self.file_posts.push((query, posts));
        }
    }

    /// What NCBI writes after each query's report as it makes and writes the query's PSSM, by
    /// searched query, in the order of the queries: the messages of the PSSM's reads (lines, each
    /// a warning), then the posts of its ASN.1 file's write.
    pub(crate) fn take_read_warnings(&mut self) -> Vec<(usize, Posts)> {
        let mut file_posts = std::mem::take(&mut self.file_posts).into_iter().peekable();
        let mut taken = Vec::new();
        for (query, messages) in std::mem::take(&mut self.read_warnings) {
            taken.extend(std::iter::from_fn(|| {
                file_posts.next_if(|&(of, _)| of < query)
            }));
            let mut posts = Posts::of_lines(messages.bytes());
            if let Some((_, written)) = file_posts.next_if(|&(of, _)| of == query) {
                posts.append(&written);
            }
            if !posts.is_empty() {
                taken.push((query, posts));
            }
        }
        taken.extend(file_posts);
        taken
    }

    /// One query's PSSM: its ASN.1 file staged and its ASCII PSSM kept. `Err` is LOSAT's refusal
    /// of an ASN.1 file that cannot be written as NCBI does (NCBI stops after that query's report).
    // NCBI's stream writes through a 4 KiB buffer (`COStreamBuffer`); the first write that fails
    // throws a `CIOException` (`eWrite`, "write fault"), the object stream reports it once as an
    // error of its own (`SetFailFlags`: `Error: (802.5) [psiblast] CObjectOStream: error at byte
    // N: stack is empty: write fault`, N being where the buffer first overflowed), and psiblast's
    // `CATCH_ALL` ends the run with status 0, without the epilog and the later queries. LOSAT
    // stops at the same place, with its own refusal instead of NCBI's message and status, and
    // without NCBI's partial file.
    // NCBI reference (598d8ae6): c++/src/util/strbuffer.cpp:871-924 and 921-924
    // NCBI reference (598d8ae6): c++/src/serial/objostr.cpp:373-384
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:242-247
    pub(crate) fn add(&mut self, asn: Option<&[u8]>, ascii: Option<Vec<u8>>) -> Result<()> {
        let file = asn.and(self.asn.as_mut()).map(OutputFileReset::next_file);
        let asn = match (file, asn) {
            (Some(file), Some(bytes)) => {
                let staged = self.stage(&file, bytes)?;
                Some((file, staged))
            }
            _ => None,
        };
        if let Some((file, Staged::File(_) | Staged::InPlace(_) | Staged::Receiver(_, _))) = &asn {
            if file.replaced {
                if let Some(previous) = self.last_replaced.replace(self.pending.len()) {
                    if let Some((_, staged)) = self.pending[previous].asn.as_mut() {
                        if let Staged::File(staging) = std::mem::replace(staged, Staged::Done) {
                            let _ = std::fs::remove_file(staging);
                        }
                    }
                }
            }
        }
        self.pending.push(PendingPssm { asn, ascii });
        Ok(())
    }

    /// The ASN.1 bytes of `file`, staged (`PssmFiles`).
    fn stage(&mut self, file: &PssmFile, asn: &[u8]) -> Result<Staged> {
        let path = &file.path;
        let cannot_write = |error: std::io::Error| {
            anyhow::anyhow!(
                "a -out_pssm file that cannot be written ({}: {error}), where NCBI BLAST+'s psiblast stops after the query's report with a serial-stream write error, is not supported by LOSAT's PSIBLAST",
                path.display()
            )
        };
        if matches!(
            self.run_files,
            crate::report::run_files::RunFiles::Receiver(_)
        ) {
            use crate::report::run_files::PssmTarget;
            return match self.run_files.pssm_target("-out_pssm", path.as_os_str()) {
                PssmTarget::File(host_name) => Ok(Staged::Receiver(host_name, asn.to_vec())),
                PssmTarget::Discard => Ok(Staged::Null),
                PssmTarget::Full => Err(cannot_write(std::io::Error::from(
                    std::io::ErrorKind::StorageFull,
                ))),
                PssmTarget::Unopened => Err(cannot_write(std::io::Error::new(
                    std::io::ErrorKind::NotFound,
                    "cannot open file",
                ))),
                PssmTarget::Refused(reason) => Err(anyhow::anyhow!("{reason}")),
            };
        }
        let existing = match target_kind(path) {
            TargetKind::Absent => false,
            TargetKind::Regular => true,
            TargetKind::Device { null } => {
                return test_device(path, null, asn).map_err(cannot_write)
            }
            TargetKind::Unopenable => return test_device(path, false, asn).map_err(cannot_write),
            TargetKind::Unsupported(kind) => return Err(unsupported_target(path, kind)),
            TargetKind::Unreachable(error) => return Err(cannot_write(error)),
        };
        if path.file_name().is_none() {
            return Err(cannot_write(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "not a file name",
            )));
        }
        let directory = match path.parent() {
            Some(parent) if !parent.as_os_str().is_empty() => parent,
            _ => std::path::Path::new("."),
        };
        let (staging, uid) = match self.write_staging(directory, asn) {
            Ok(staged) => staged,
            Err(_) if existing => {
                // Where the directory takes no new file, NCBI's version 0 cannot remove the file
                // either and opens it in place, as a numbered file.
                std::fs::OpenOptions::new()
                    .write(true)
                    .open(path)
                    .map_err(cannot_write)?;
                return Ok(Staged::InPlace(asn.to_vec()));
            }
            Err(error) => return Err(cannot_write(error)),
        };
        if existing {
            if let Err(error) = std::fs::OpenOptions::new().write(true).open(path) {
                if immutable(&error) || !(file.replaced && removable(path, directory, uid)) {
                    let _ = std::fs::remove_file(&staging);
                    return Err(cannot_write(error));
                }
            }
        }
        Ok(Staged::File(staging))
    }

    /// Writes `bytes` to a new staging file `.losat-<pid>-<n>` in the target's `directory`; its
    /// path and owner. The pid is `utils::process_id` (0 on the Wasm targets, which have none:
    /// `std::process::id()` panics there).
    fn write_staging(
        &mut self,
        directory: &std::path::Path,
        bytes: &[u8],
    ) -> std::io::Result<(std::path::PathBuf, u32)> {
        loop {
            self.staging_count += 1;
            let staging = directory.join(format!(
                ".losat-{}-{}",
                crate::utils::process_id::process_id(),
                self.staging_count
            ));
            let mut file = match std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&staging)
            {
                Ok(file) => file,
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => return Err(error),
            };
            let written = file
                .write_all(bytes)
                .and_then(|()| file.metadata())
                .map(|metadata| owner(&metadata));
            return match written {
                Ok(uid) => Ok((staging, uid)),
                Err(error) => {
                    drop(file);
                    let _ = std::fs::remove_file(&staging);
                    Err(error)
                }
            };
        }
    }

    /// After the reports: each query's ASN.1 file (`GetStream` on the target, then the staged
    /// bytes), then its ASCII file (not checked, as NCBI's `ofstream`). `Err`: an ASN.1 file that
    /// `add` found writable failed, after the reports (see `PssmFiles`).
    pub(crate) fn commit(mut self) -> Result<()> {
        for index in 0..self.pending.len() {
            if let Some((file, staged)) = self.pending[index].asn.as_mut() {
                let commit_err = |error: std::io::Error| {
                    anyhow::anyhow!(
                        "a -out_pssm file that could not be written after the reports ({}: {error}), although LOSAT found it writable before them, is not supported by LOSAT's PSIBLAST",
                        file.path.display()
                    )
                };
                match std::mem::replace(staged, Staged::Done) {
                    Staged::Receiver(host_name, bytes) => {
                        if let crate::report::run_files::RunFiles::Receiver(receiver) =
                            &self.run_files
                        {
                            receiver.receive(&host_name, &bytes).map_err(commit_err)?;
                        }
                    }
                    other => {
                        let bytes = match other {
                            Staged::File(staging) => {
                                let bytes = std::fs::read(&staging);
                                let _ = std::fs::remove_file(&staging);
                                Some(bytes)
                            }
                            Staged::InPlace(bytes) | Staged::Device(bytes) => Some(Ok(bytes)),
                            Staged::Receiver(_, _) | Staged::Null | Staged::Done => None,
                        };
                        if let Some(bytes) = bytes {
                            bytes
                                .and_then(|bytes| file.open()?.write_all(&bytes))
                                .map_err(commit_err)?;
                        }
                    }
                }
            }
            if let (Some(stream), Some(ascii)) =
                (self.ascii.as_mut(), self.pending[index].ascii.take())
            {
                let file = stream.next_file();
                match &self.run_files {
                    crate::report::run_files::RunFiles::Receiver(receiver) => {
                        if let crate::report::run_files::PssmTarget::File(host_name) = self
                            .run_files
                            .pssm_target("-out_ascii_pssm", file.path.as_os_str())
                        {
                            let _ = receiver.receive(&host_name, &ascii);
                        }
                    }
                    crate::report::run_files::RunFiles::Disk => {
                        if let Ok(mut out) = file.open() {
                            let _ = out.write_all(&ascii);
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

impl Drop for PssmFiles {
    // Out of line: it runs where the function that holds the files returns (`run_call`, or
    // `report_prep::report_values` on an error after the files are staged).
    #[cold]
    #[inline(never)]
    fn drop(&mut self) {
        for pending in &self.pending {
            if let Some((_, Staged::File(staging))) = &pending.asn {
                let _ = std::fs::remove_file(staging);
            }
        }
    }
}

/// A character device's test (`PssmFiles`): opened for writing and written zero bytes; what
/// `commit` writes to it.
fn test_device(path: &std::path::Path, null: bool, asn: &[u8]) -> std::io::Result<Staged> {
    let mut device = std::fs::OpenOptions::new().write(true).open(path)?;
    // A write of zero bytes is passed to the device (`write(2)`).
    let written = device.write(&[])?;
    debug_assert_eq!(written, 0);
    Ok(if null {
        Staged::Null
    } else {
        Staged::Device(asn.to_vec())
    })
}

fn unsupported_target(path: &std::path::Path, kind: &str) -> anyhow::Error {
    anyhow::anyhow!(
        "an -out_pssm file that is {kind} ({}) is not supported by LOSAT's PSIBLAST",
        path.display()
    )
}

/// Before the search: refuses PSSM files that LOSAT cannot write as NCBI does. NCBI writes each
/// query's files after its report through streams of their own, so a PSSM file that is also the
/// standard output, the standard error or the `-out` file gets its bytes among the reports' in an
/// order that their buffers make, which LOSAT does not reproduce (it writes the PSSM files after
/// all reports). Every name the run can write is checked: `F`, or with `-save_each_pssm` the next
/// `<number of queries>` names `F.N` of the run (`PssmVersions`: from `F.1` but in the library's
/// later query batches); a name that does not exist cannot be one of those streams. The same
/// file is the same device and inode (`metadata` follows `/dev/stdout` to the stream's file),
/// or for a character device the same device number (`/dev/null` excepted: no order of its bytes
/// shows); without inode numbers, the same path form as `-out` (`utils::same_file`: the canonical
/// path, or on WASI the path by name, which does not see a link to the file), or one of the names of the
/// standard streams (`ReportStreams`, taken once). Also refused, for `-out_pssm`: a named pipe, a
/// block device and a symbolic link to no file (`PssmFiles`).
pub(crate) fn check_pssm_files(
    request: &PssmRequest,
    queries: usize,
    out: Option<&std::path::Path>,
) -> Result<()> {
    let is_receiver = matches!(
        request.run_files,
        crate::report::run_files::RunFiles::Receiver(_)
    );
    let check = (!is_receiver).then(|| PssmFileCheck::new(out));
    for ((option, path), first) in [
        ("-out_pssm", request.asn_path.as_deref()),
        ("-out_ascii_pssm", request.ascii_path.as_deref()),
    ]
    .into_iter()
    .zip(request.versions.next())
    {
        let Some(path) = path else {
            continue;
        };
        let names: Vec<std::path::PathBuf> = if request.save_each {
            (first..first.saturating_add(u32::try_from(queries).unwrap_or(u32::MAX)))
                .map(|n| numbered_pssm_path(path, n))
                .collect()
        } else {
            vec![path.to_path_buf()]
        };
        for name in &names {
            if let Some(check) = &check {
                check.check(option, name)?;
            } else if let crate::report::run_files::PssmTarget::Refused(reason) =
                request.run_files.pssm_target(option, name.as_os_str())
            {
                anyhow::bail!("{reason}");
            }
        }
    }
    Ok(())
}

/// `check_pssm_files`' check of one name, also made by the round driver before each write
/// (`rounds::RoundFiles`), whose numbered names are not known before the search.
pub(crate) struct PssmFileCheck {
    streams: ReportStreams,
}

impl PssmFileCheck {
    pub(crate) fn new(out: Option<&std::path::Path>) -> Self {
        Self {
            streams: ReportStreams::new(out),
        }
    }

    /// The refusal of the file `name` of `option` (`-out_pssm` or `-out_ascii_pssm`) where it
    /// is a report stream, or for `-out_pssm` a named pipe, a block device or a symbolic link to
    /// no file.
    pub(crate) fn check(&self, option: &str, name: &std::path::Path) -> Result<()> {
        if let Some(stream) = self.streams.stream_of(name) {
            anyhow::bail!(
                "an {option} file that is also {stream} ({}) is not supported by LOSAT's PSIBLAST",
                name.display()
            );
        }
        if option == "-out_pssm" {
            if let TargetKind::Unsupported(kind) = target_kind(name) {
                return Err(unsupported_target(name, kind));
            }
        }
        Ok(())
    }
}

/// Whether `path` names a regular file (also through a symbolic link), which a failed write
/// leaves partly written (`rounds::RoundFiles`).
pub(crate) fn is_regular_file(path: &std::path::Path) -> bool {
    matches!(target_kind(path), TargetKind::Regular)
}

/// The report streams (`check_pssm_files`): the standard output, the standard error (`fstat` of a
/// copy of each descriptor) and the `-out` file, as files.
#[cfg(unix)]
struct ReportStreams {
    streams: Vec<(&'static str, std::fs::Metadata)>,
    null: Option<u64>,
}

#[cfg(unix)]
impl ReportStreams {
    fn new(out: Option<&std::path::Path>) -> Self {
        use std::os::fd::AsFd;
        let stream = |fd: std::os::fd::BorrowedFd<'_>| {
            fd.try_clone_to_owned()
                .ok()
                .and_then(|fd| std::fs::File::from(fd).metadata().ok())
        };
        let streams = [
            ("the standard output", stream(std::io::stdout().as_fd())),
            ("the standard error", stream(std::io::stderr().as_fd())),
            (
                "the -out file",
                out.and_then(|out| std::fs::metadata(out).ok()),
            ),
        ]
        .into_iter()
        .filter_map(|(name, metadata)| Some((name, metadata?)))
        .collect();
        Self {
            streams,
            null: null_device(),
        }
    }

    /// The report stream that the file `name` is; `/dev/null`, where no order of bytes shows, is
    /// none.
    fn stream_of(&self, name: &std::path::Path) -> Option<&'static str> {
        use std::os::unix::fs::MetadataExt;
        let target = std::fs::metadata(name).ok()?;
        let device = device_number(&target);
        if device.is_some() && device == self.null {
            return None;
        }
        self.streams.iter().find_map(|(stream, metadata)| {
            ((metadata.dev(), metadata.ino()) == (target.dev(), target.ino())
                || device.is_some_and(|rdev| device_number(metadata) == Some(rdev)))
            .then_some(*stream)
        })
    }
}

#[cfg(not(unix))]
struct ReportStreams {
    out: Option<std::path::PathBuf>,
}

#[cfg(not(unix))]
impl ReportStreams {
    fn new(out: Option<&std::path::Path>) -> Self {
        Self {
            out: out.and_then(crate::utils::same_file::path_identity),
        }
    }

    fn stream_of(&self, name: &std::path::Path) -> Option<&'static str> {
        match name.to_str() {
            Some("/dev/stdout" | "/dev/fd/1" | "/proc/self/fd/1") => {
                return Some("the standard output")
            }
            Some("/dev/stderr" | "/dev/fd/2" | "/proc/self/fd/2") => {
                return Some("the standard error")
            }
            _ => {}
        }
        let target = crate::utils::same_file::path_identity(name)?;
        (self.out.as_ref() == Some(&target)).then_some("the -out file")
    }
}

#[cfg(test)]
#[path = "pssm_tests.rs"]
mod tests;
