use std::path::Path;

use crate::api::local_blast::{FormatOutput, FormatRefusal, OutputSink, ReportOutputs};
use crate::blastinput::fasta_reader::{read_all, FastaInputSource, ReaderConfig};
use crate::cli::NativeError;
use crate::common::GapEditOp::{self, Del, Ins, Sub};
use crate::core::composition_adjustment::redo_alignment::EMatrixAdjustRule;
use crate::stats::tables::KarlinParams;

use super::blast_xml::XmlValues;
use super::search_gapped::GappedHsp;
use super::stage_d_linking::{LinkedHsp, LinkedHspList};
use super::stage_d_pipeline::{LocalStageDScoring, SubjectSetNumbers};
use super::stage_d_results::{KappaHspPayload, KappaResultHitList, KappaResultList};
use super::stage_d_stats::{LocalContextCutoffs, LocalContextLength, LocalSubjectParameters};
use super::stage_e_report::{render, ReportRanges, ShownSubjects};
use super::x_subject_store::SubjectLetters;

/// The outcome of running a TBLASTN `-subject` report writer on one synthetic gapped HSP.
pub(crate) struct SubjectReportOutcome {
    pub(crate) stdout: String,
    pub(crate) stderr: String,
    pub(crate) error: Option<NativeError>,
}

/// Runs TBLASTN's `-subject` mapping (`put_in_report_order(true)`) and `stage_e_report::render`
/// (`"5"`, `"16"`, or `"7 qseqid sseqid qcovs"`) on a single gapped HSP with edit script `ops`.
pub(crate) fn run_subject_report(ops: Vec<GapEditOp>, outfmt: &str) -> SubjectReportOutcome {
    let [q_len, s_len] = ops.iter().fold([0i32, 0i32], |[q, s], op| match *op {
        Sub(n) => [q + n as i32, s + n as i32],
        Ins(n) => [q + n as i32, s],
        Del(n) => [q, s + n as i32],
    });
    let (align_length, gap_opens, gap_letters) =
        ops.iter()
            .fold((0usize, 0usize, 0usize), |(len, opens, letters), op| {
                let n = op.num() as usize;
                match op {
                    Sub(_) => (len + n, opens, letters),
                    Ins(_) | Del(_) => (len + n, opens + 1, letters + n),
                }
            });
    let q_fa = format!(">q1\n{}\n", "A".repeat((q_len as usize).max(20)));
    let s_fa = format!(">s1\n{}\n", "GCT".repeat((s_len as usize).max(20)));
    let queries = read_all(
        &mut FastaInputSource::from_bytes(
            q_fa.as_bytes(),
            ReaderConfig::query("TBLASTN", true, false),
        ),
        &mut |_| Ok(()),
    )
    .unwrap();
    let subjects = read_all(
        &mut FastaInputSource::from_bytes(
            s_fa.as_bytes(),
            ReaderConfig::subject("TBLASTN", false, false),
        ),
        &mut |_| Ok(()),
    )
    .unwrap();
    let subject_seqs: Vec<&[u8]> = subjects.iter().map(|r| r.seq()).collect();
    let mut hitlist = KappaResultHitList::new(500).unwrap();
    hitlist
        .update(KappaResultList {
            oid: 0,
            hsps: LinkedHspList {
                hsps: vec![LinkedHsp {
                    context: 0,
                    hsp: GappedHsp {
                        frame: 1,
                        score: 80,
                        q_start: 0,
                        q_end: q_len,
                        q_gapped_start: 0,
                        s_start: 0,
                        s_end: s_len,
                        s_gapped_start: 0,
                    },
                    num: 1,
                    evalue: 1e-6,
                    source_index: 0,
                }],
                best_evalue: 1e-6,
            },
            payloads: vec![KappaHspPayload {
                bit_score: 35.0,
                num_ident: 10,
                num_positives: 10,
                report_num_ident: 10,
                report_mismatches: 0,
                align_length,
                mismatches: 0,
                gap_opens,
                gap_letters,
                edit_script: ops,
                matrix_adjust_rule: EMatrixAdjustRule::DontAdjustMatrix,
                query_left_out: 0,
                nullb: None,
            }],
        })
        .unwrap();
    let mut hitlists = vec![hitlist];
    let parameters = LocalSubjectParameters {
        lengths: vec![LocalContextLength {
            length_adjustment: 0,
            eff_searchsp: 1000,
            subject_stat_length: (s_len as i64).max(20),
        }],
        cutoffs: vec![LocalContextCutoffs {
            hit_cutoff: 30,
            hit_cutoff_max: 45,
            word_cutoff: 12,
            word_xdrop: 16,
        }],
        link: None,
        hit_cutoff_min: 30,
        word_cutoff_min: 12,
        prelim_evalue: 10.0,
    };
    let ungapped_karlin = [KarlinParams {
        lambda: 0.3176,
        k: 0.134,
        h: 0.4012,
        alpha: 0.7916,
        beta: -3.2,
    }];
    let mut stdout = Vec::new();
    let mut stderr = Vec::new();
    let result = {
        let mut outputs = ReportOutputs {
            other_outfmts: &[],
            html_posts: None,
            formats: vec![FormatOutput {
                outfmt: outfmt.as_bytes(),
                refusal: FormatRefusal::FailsCall,
                sink: OutputSink::Writer(&mut stdout),
            }],
            diagnostics: &mut stderr,
            later_warnings: None,
            hits: None,
            observer: None,
        };
        render(
            &mut outputs,
            &queries,
            &subjects,
            SubjectLetters::Records(&subject_seqs),
            &ShownSubjects::Fasta {
                path: Path::new("subject.fna"),
                parsed: &[],
            },
            SubjectSetNumbers::default(),
            &mut hitlists,
            &parameters,
            &ungapped_karlin,
            &[true],
            &[false],
            LocalStageDScoring::default(),
            "BLOSUM62",
            Default::default(),
            Default::default(),
            1,
            None,
            None,
            false,
            &[],
            &[],
            &ReportRanges {
                queries: &Default::default(),
                subjects: &Default::default(),
                epilog: true,
                query_scopes: None,
            },
            XmlValues {
                expect: 10.0,
                soft_masking: false,
                gapped: true,
                cbs: 0,
                db_gencode: 1,
                identity_warned: 0,
                filtered_empty: false,
                pssm: false,
            },
            Default::default(),
            None,
            None,
            &Default::default(),
            None,
            false,
            &[],
        )
    };
    let error = match result {
        Ok(()) => None,
        Err(err) => Some(
            err.downcast::<NativeError>()
                .expect("expected NCBI NativeError"),
        ),
    };
    SubjectReportOutcome {
        stdout: String::from_utf8(stdout).unwrap(),
        stderr: String::from_utf8(stderr).unwrap(),
        error,
    }
}
