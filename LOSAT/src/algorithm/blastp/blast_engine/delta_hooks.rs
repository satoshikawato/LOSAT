//! deltablast's domain search and its PSSM of conserved domains (series PR, item PS13a): the part
//! of `CDeltaBlast::Run` before round 1. The domain search is the RPS engine on `-rpsdb`
//! (`rps_hooks`) run as a call of the preliminary stage only (`PrelimOnly::Domain`: for RPS-BLAST
//! the whole engine, traceback included, is that stage); `algorithm::rps_search::cdd_input` makes
//! the alignment of the domains of its hits, `core::psi::cd_stages` the PSSM, and the round driver
//! (`psi_round::run_queries`) searches it as a later round of psiblast. Cold and out of line:
//! blastp's default path only tests `ResolvedBlastpArgs::delta`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/deltablast.cpp:83-215

use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

use anyhow::{anyhow, bail, Result};

use super::{
    database, open_search_database, rps_hooks, run_resolved_in_pool, BlastpRanges, BlastpReports,
    BlastpSubjectSource, PrelimOnly,
};
use crate::algorithm::blastp::args::{BlastpArgs, DeltaOptions, ResolvedBlastpArgs};
use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::algorithm::psiblast::pssm::{cdd_query_pssm, PssmQueryRecord, PssmRequest, QueryPssm};
use crate::algorithm::pssm_search::in_pssm::{engine_score_block, query_ungapped_block};
use crate::algorithm::rps_search::cdd_input::{self, open_error, Hit};
use crate::api::local_blast::ReportOutputs;
use crate::blastinput::app::engine_error;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::blastinput::seq_range;
use crate::common::GapEditOp;
use crate::core::psi::create_pssm_from_cdd;
use crate::seqdb::rps::CdProfileFiles;
use crate::seqid::model::{Dbtag, ObjectId, ParsedSeqId};

/// deltablast's own part of a run: its options, and, once the first domain search needs them, the
/// domain database with the files of its profiles. Shared by the rounds of every query through
/// `ResolvedBlastpArgs::delta`.
pub(crate) struct DeltaRun {
    options: DeltaOptions,
    lcase_masking: bool,
    /// The application's settings (`search_cli`), which the domain database is opened with.
    settings: Mutex<Option<ApplicationSettings>>,
    side: Mutex<Option<Arc<DomainSide>>>,
}

impl std::fmt::Debug for DeltaRun {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DeltaRun")
            .field("options", &self.options)
            .field("lcase_masking", &self.lcase_masking)
            .finish_non_exhaustive()
    }
}

/// The domain database and the RPS search's options (NCBI's factory of the RPS-BLAST options
/// handle, with the domain inclusion threshold as its E-value and no filter).
struct DomainSide {
    args: ResolvedBlastpArgs,
    database: database::BlastpDatabase,
    /// `.wcounts` and `.obsr`, opened by the first query's PSSM (`x_FillHitsData`), after the
    /// domain search whether it found a domain or not.
    files: Mutex<Option<Arc<CdProfileFiles>>>,
}

/// One query's PSSM of conserved domains: what `CDeltaBlast::Run` makes for it before round 1.
pub(super) struct CddQuery<'a> {
    /// The query as the reader kept it (its lower case is `-lcase_masking`'s).
    pub(super) record: &'a FastaRecord,
    /// The query's Seq-align id (`Query_<n>`, or the id `-parse_deflines` read).
    pub(super) id: &'a ParsedSeqId,
    /// Whether the query is alone in its query batch.
    pub(super) alone: bool,
}

/// The PSSM of a query's conserved domains and the query of the PSSM's Bioseq, which every later
/// round reads: the record's letters with the lower case of `-lcase_masking` as X
/// (`CCddInputData`'s query, the domain search's query block). NCBI's later rounds start from that
/// Bioseq (deltablast_app.cpp:420-431): its later PSSMs are made of those letters
/// (`PsiBlastComputePssmFromAlignment`), and its round queries are made of them, without a
/// lower-case mask (`CPsiBlast` makes its query factory of the PSSM's Bioseq,
/// psiblast_impl.cpp:59-68 and 114-119).
pub(super) struct CddPssm {
    pub(super) pssm: QueryPssm,
    pub(super) query: FastaRecord,
}

impl DeltaRun {
    pub(crate) fn new(options: DeltaOptions, lcase_masking: bool) -> Self {
        Self {
            options,
            lcase_masking,
            settings: Mutex::new(None),
            side: Mutex::new(None),
        }
    }

    /// The application's settings, given by the command line before the first domain search.
    pub(super) fn set_settings(&self, settings: ApplicationSettings) {
        *self.settings.lock().expect("deltablast's settings") = Some(settings);
    }

    /// The domain search's side, opened by the first query.
    fn side(&self, diagnostics: &mut (dyn std::io::Write + Send)) -> Result<Arc<DomainSide>> {
        let mut side = self.side.lock().expect("deltablast's domain database");
        if let Some(side) = side.as_ref() {
            return Ok(Arc::clone(side));
        }
        let settings = self
            .settings
            .lock()
            .expect("deltablast's settings")
            .clone()
            .ok_or_else(|| {
                anyhow!("LOSAT internal error: deltablast has no application settings")
            })?;
        let opened = Arc::new(self.open_side(&settings, diagnostics)?);
        *side = Some(Arc::clone(&opened));
        Ok(opened)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/deltablast.cpp:195-215
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_rps_options.cpp:53-78
    fn open_side(
        &self,
        settings: &ApplicationSettings,
        diagnostics: &mut (dyn std::io::Write + Send),
    ) -> Result<DomainSide> {
        // deltablast made the domain search's arguments (`DeltaOptions::domain`). A threshold of 0
        // or less is NCBI's options error of the domain search, which comes after the domain
        // database is opened (below): the check of its options runs where the search does.
        let ethresh = self.options.domain_inclusion_ethresh;
        let mut args: BlastpArgs = (*self.options.domain).clone();
        args.matrix_search = Some(settings.matrix_search.clone());
        let mut resolved = args.check_options(diagnostics)?;
        let mut outputs = ReportOutputs {
            formats: Vec::new(),
            diagnostics,
            later_warnings: None,
            hits: None,
            observer: None,
            other_outfmts: &[],
            html_posts: None,
        };
        let database = open_search_database(&args, &mut resolved, None, settings, &mut outputs)?
            .ok_or_else(|| anyhow!("LOSAT internal error: deltablast opened no domain database"))?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1518-1523
        if ethresh <= 0.0 {
            return Err(crate::algorithm::psiblast::pssm::blast_options_error(
                "expect value or cutoff score must be greater than zero",
            ));
        }
        if rps_hooks::volume_count(&database) > 1 {
            bail!("a domain database of several volumes is not supported by LOSAT's DELTABLAST");
        }
        Ok(DomainSide {
            args: resolved,
            database,
            files: Mutex::new(None),
        })
    }

    /// The PSSM of the query's conserved domains (`CDeltaBlast::Run` before its round 1): the
    /// domain search, `CCddInputData` and the PSSM engine, with the PSSM's files as the options
    /// ask (`request`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/deltablast.cpp:93-150
    #[allow(clippy::too_many_arguments)]
    pub(super) fn cdd_pssm(
        &self,
        args: &ResolvedBlastpArgs,
        request: &PssmRequest,
        query: &CddQuery<'_>,
        pool: &crate::utils::threading::SearchPool<'_>,
        diagnostics: &mut (dyn std::io::Write + Send),
    ) -> Result<CddPssm> {
        let side = self.side(&mut *diagnostics)?;
        // The query's letters as the domain search's set-up leaves them in the query block that
        // the PSSM engine reads: with `-lcase_masking` its lower case is hard-masked.
        let letters = if self.lcase_masking {
            let frame = crate::algorithm::blastp::encoding::encode_protein_query_frame_lcase(
                query.record.seq(),
                None,
                false,
                false,
            );
            frame.aa_seq[1..=query.record.seq().len()].to_vec()
        } else {
            let encoded =
                crate::algorithm::blastp::encoding::encode_protein_sequence(query.record.seq());
            encoded.aa_seq[1..=encoded.aa_len].to_vec()
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/cdd_pssm_input.cpp:95-97
        // A query whose letters give no Karlin-Altschul parameters has no Seq-align set in the
        // domain search's results, which `CCddInputData` refuses (`NULL alignments`, exit 3; NCBI
        // 2.17.0 probes `x1`, `x1_it2`, `x2_it2`, `xs_it2` of PS13a's fix round 1) where the query is
        // alone in its batch. In a batch of several queries NCBI's deltablast crashes there
        // (exit 139, probes `x2`, `xs`, `sxs`).
        if query_ungapped_block(side.args.scoring.matrix, &letters).is_none() {
            if query.alone {
                return Err(engine_error("NULL alignments"));
            }
            bail!(
                "a query batch holding a query whose letters give no Karlin-Altschul parameters and other queries (with which NCBI BLAST+'s deltablast crashes) is not supported by LOSAT's DELTABLAST"
            );
        }
        let lists = side.search(query.record, pool, diagnostics)?;
        let hits = side.hits(&lists)?;
        let hits = cdd_input::included_hits(hits, self.options.domain_inclusion_ethresh);
        let hits = cdd_input::remove_multiple_cd_hits(hits)
            .map_err(|refusal| anyhow!(refusal.message()))?;
        // The files open after the domain search, also for a query without a hit.
        let files = side.profile_files()?;
        let msa = cdd_input::build_msa(&letters, &hits, &files)?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/pssm_engine.cpp:412-445
        let (matrix, sbp) = engine_score_block(
            &args.matrix_name,
            args.scoring.gap_open,
            args.scoring.gap_extend,
            &letters,
        )?
        .map_err(|message| engine_error(&message))?;
        let Some(sbp) = sbp else {
            bail!(
                "-matrix IDENTITY is not supported by LOSAT's DELTABLAST (the matrix has no frequency ratios)"
            );
        };
        // The search's matrix gives the query Karlin-Altschul parameters where the domain search's
        // did (above); a matrix where it does not is LOSAT's refusal, not reached with NCBI's
        // matrices.
        if query_ungapped_block(matrix, &letters).is_none() {
            bail!(
                "a query whose letters give Karlin-Altschul parameters with the domain search's matrix but not with -matrix is not supported by LOSAT's DELTABLAST"
            );
        }
        let pssm = create_pssm_from_cdd(&msa, 0, &sbp).map_err(|error| {
            engine_error(crate::algorithm::psiblast::pssm::error_code_to_string(
                error,
            ))
        })?;
        let defline = crate::report::archive::defline_of(query.record);
        let pssm = cdd_query_pssm(
            request,
            pssm,
            &letters,
            PssmQueryRecord {
                letters: query.record.seq(),
                id: query.id,
                title: query.record.title(),
            },
            &defline,
            &args.matrix_name,
        )?;
        let mut cdd_query = query.record.clone();
        cdd_query.warnings.clear();
        if self.lcase_masking {
            for letter in &mut cdd_query.sequence {
                if letter.is_ascii_lowercase() {
                    *letter = b'X';
                }
            }
        }
        Ok(CddPssm {
            pssm,
            query: cdd_query,
        })
    }
}

impl DomainSide {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/deltablast.cpp:195-215
    /// The domain search of one query: the hit lists of the RPS engine, final.
    fn search(
        &self,
        record: &FastaRecord,
        pool: &crate::utils::threading::SearchPool<'_>,
        diagnostics: &mut (dyn std::io::Write + Send),
    ) -> Result<Vec<Option<BlastpHitList>>> {
        let records = std::slice::from_ref(record);
        let query_input = seq_range::QueryInput::whole(records);
        let mut outputs = ReportOutputs {
            formats: Vec::new(),
            diagnostics,
            later_warnings: None,
            hits: None,
            observer: None,
            other_outfmts: &[],
            html_posts: None,
        };
        let mut lists = Vec::new();
        run_resolved_in_pool(
            self.args.clone(),
            records,
            BlastpSubjectSource::Database(&self.database),
            "",
            "",
            &mut outputs,
            pool,
            Some(PrelimOnly::Domain(&mut lists)),
            &BlastpRanges {
                input_records: records,
                query_input: &query_input,
                subject_placements: &seq_range::Placements::default(),
                x_ordinal_base: 0,
                searched_base: 0,
                split_search_spaces: None,
            },
            &BlastpReports::prelim_only(false),
            super::x_lookup_ahead::XAheadBlock::none(),
        )?;
        Ok(lists)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/cdd_pssm_input.cpp:151-165
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/cdd_pssm_input.cpp:335-376
    /// The Seq-aligns of the domain search (one for each HSP, in the lists' order) as `CHit`s: the
    /// domain's Seq-align id, the HSP's E-value as the Seq-align stores it, and its aligned runs
    /// (the gaps of either sequence are dropped).
    fn hits(&self, lists: &[Option<BlastpHitList>]) -> Result<Vec<Hit>> {
        let Some(list) = lists.first().and_then(Option::as_ref) else {
            return Ok(Vec::new());
        };
        let count = list.hsplist_count.min(list.hsplist_array.len());
        let mut ids: BTreeMap<u32, ParsedSeqId> = BTreeMap::new();
        let mut hits = Vec::new();
        for hsp_list in &list.hsplist_array[..count] {
            if !ids.contains_key(&hsp_list.oid) {
                ids.insert(hsp_list.oid, self.profile_id(hsp_list.oid)?);
            }
            let id = &ids[&hsp_list.oid];
            for hsp in &hsp_list.hsps {
                let mut aligned = Vec::new();
                match hsp.gap_info.as_deref() {
                    None => aligned.push((
                        hsp.q_start as i32 - 1,
                        hsp.s_start as i32 - 1,
                        (hsp.q_end + 1 - hsp.q_start) as i32,
                    )),
                    Some(ops) => {
                        let mut query = hsp.q_start as i32 - 1;
                        let mut subject = hsp.s_start as i32 - 1;
                        for op in ops {
                            match *op {
                                GapEditOp::Sub(n) => {
                                    aligned.push((query, subject, n as i32));
                                    query += n as i32;
                                    subject += n as i32;
                                }
                                GapEditOp::Del(n) => subject += n as i32,
                                GapEditOp::Ins(n) => query += n as i32,
                            }
                        }
                    }
                }
                hits.push(cdd_input::hit_of_segments(
                    id.clone(),
                    hsp_list.oid,
                    crate::report::serial::seqalign::seqalign_evalue(hsp.e_value),
                    aligned,
                ));
            }
        }
        Ok(hits)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqinfosrc_aux.cpp:51-69
    /// The Seq-align id of the profile `oid` (the id that `SeqidToOid` finds back, or the
    /// ordinal id of a database whose profiles have no parsed Seq-ids).
    fn profile_id(&self, oid: u32) -> Result<ParsedSeqId> {
        use crate::seqdb::SequenceHeader;
        let seqdb = self.database.seqdb();
        let header = seqdb
            .read_header_of_oid(oid)
            .map_err(|error| anyhow!("{error}"))?;
        let id = match header {
            SequenceHeader::Parsed(set) => set.alignment_id().cloned().ok_or_else(|| {
                anyhow!(
                    "a domain whose header has no Seq-id is not supported by LOSAT's DELTABLAST"
                )
            })?,
            SequenceHeader::Ordinal(_) => ParsedSeqId::General(Dbtag {
                db: b"BL_ORD_ID".to_vec(),
                tag: ObjectId::Id(oid as i32),
            }),
        };
        // `CCddInputData` reads the files by the OID that the id finds in the database.
        let found = seqdb
            .seqid_lookup()
            .map_err(|error| anyhow!("{error}"))?
            .seqid_to_oid(&id)
            .map_err(|error| anyhow!("{error}"))?;
        if found != Some(oid) {
            bail!(
                "a domain database in which the Seq-align id of a profile finds another profile or none is not supported by LOSAT's DELTABLAST"
            );
        }
        Ok(id)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/rps_aux.cpp:576-592
    /// `.wcounts` and `.obsr` of the domain database, opened once.
    fn profile_files(&self) -> Result<Arc<CdProfileFiles>> {
        let mut files = self.files.lock().expect("deltablast's profile files");
        if let Some(files) = files.as_ref() {
            return Ok(Arc::clone(files));
        }
        let base = self.database.seqdb().volumes()[0].base().to_path_buf();
        let opened = CdProfileFiles::open(&base, "DELTABLAST").map_err(open_error)?;
        let opened = Arc::new(opened);
        *files = Some(Arc::clone(&opened));
        Ok(opened)
    }
}
