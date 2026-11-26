use std::path::Path;

use crate::api::local_blast::OutputSink;
use crate::cli::{Commands, NativeError, NcbiOstream};
use crate::common::GapEditOp::{self, Del, Ins, Sub};
use crate::report::blast_xml::BlastXmlWriter;
use crate::report::blast_xml2::{NoSideFiles, Xml2Kind};

use super::empty_rows::SubjectLoc;
use super::input::FastaRecord;
use super::parameters::ContextParameters;
use super::preliminary::PreliminaryHsp;
use super::query_setup::prepare_queries;
use super::report::ShownSubjects;
use super::results::{Hsp, HspList};
use super::runtime::BatchResults;

/// The outcome of running a BLASTX `-subject` report writer on one synthetic gapped HSP.
pub(crate) struct SubjectReportOutcome {
    pub(crate) stdout: String,
    pub(crate) stderr: String,
    pub(crate) error: Option<NativeError>,
}

fn options(outfmt: &str) -> super::args::ResolvedOptions {
    let cli: crate::cli::Cli = crate::cli::try_parse_from([
        "losat",
        "blastx",
        "-query",
        "unused.fna",
        "-subject",
        "unused.faa",
        "-seg",
        "no",
        "-comp_based_stats",
        "0",
        "-outfmt",
        outfmt,
    ])
    .unwrap();
    let Commands::Blastx(args) = cli.command else {
        unreachable!()
    };
    args.resolve().unwrap()
}

/// Runs BLASTX's `-subject` mapping (`put_in_report_order(SubjectLoc::Interval)`) and the
/// requested writer (`"5"`, `"16"`, or `"7 qseqid sseqid qcovs"`) on a single gapped HSP with edit
/// script `ops`.
pub(crate) fn run_subject_report(ops: Vec<GapEditOp>, outfmt: &str) -> SubjectReportOutcome {
    let [q_len, s_len] = ops.iter().fold([0i32, 0i32], |[q, s], op| match *op {
        Sub(n) => [q + n as i32, s + n as i32],
        Ins(n) => [q + n as i32, s],
        Del(n) => [q, s + n as i32],
    });
    let opts = options(outfmt);
    let query = FastaRecord {
        internal_id: "Query_1".into(),
        title: b"q1".to_vec(),
        sequence: b"GCT".repeat((q_len as usize).max(20)),
        lowercase_masks: Vec::new(),
        source: None,
        around: None,
        skipped: None,
    };
    let subject = FastaRecord {
        internal_id: "Subject_1".into(),
        title: b"s1".to_vec(),
        sequence: vec![b'A'; (s_len as usize).max(20)],
        lowercase_masks: Vec::new(),
        source: None,
        around: None,
        skipped: None,
    };
    let hsp = Hsp {
        hsp: PreliminaryHsp {
            context: 0,
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
        bit_score: 35.0,
        identity: 10,
        positive: 10,
        edit_script: ops,
        composition_method: 0,
        query_left_out: 0,
    };
    let params = vec![
        ContextParameters {
            valid: true,
            ungapped: crate::stats::tables::KarlinParams {
                lambda: 0.3176,
                k: 0.134,
                h: 0.4012,
                alpha: 0.7916,
                beta: -3.2,
            },
            search_space: 1000,
            length_adjustment: 0,
            word_cutoff: 12,
            word_xdrop_init: 16,
            word_xdrop: 38,
            hit_cutoff: 30,
            hit_cutoff_max: 45,
        };
        6
    ];
    let queries = [query];
    let subjects = [subject];
    let mut batch = BatchResults {
        query_ordinal: 0,
        in_report_order: false,
        queries: vec![vec![HspList {
            query_index: 0,
            oid: 0,
            hsps: vec![hsp],
            best_evalue: 1e-6,
        }]],
        prepared: prepare_queries(&queries, &opts).unwrap(),
        parameters: params,
        search_skipped: false,
        request_space: None,
        unlisted: Vec::new(),
        unlisted_subject_row_throws: false,
        shown: Vec::new(),
        gapped_k: None,
        fetched_utf8: None,
    };
    batch.put_in_report_order(SubjectLoc::Interval);

    let no_parsed = super::report::RecordsParsed::default();
    let shown = ShownSubjects::Fasta {
        path: Path::new("subject.faa"),
        parsed: &no_parsed,
    };
    let mut stdout = Vec::new();
    let mut stderr = Vec::new();
    let result = {
        let mut sink = OutputSink::Writer(&mut stdout);
        let mut writer = NcbiOstream::sink(&mut sink).unwrap();
        match outfmt {
            "5" => {
                let mut doc = BlastXmlWriter::new(0);
                super::blast_xml::write_blastx_xml(
                    &mut writer,
                    &mut stderr,
                    &mut doc,
                    &opts,
                    &queries,
                    &subjects,
                    &shown,
                    &[batch],
                    0,
                    false,
                    None,
                    None,
                )
                .and_then(|()| {
                    super::blast_xml::write_xml_epilog(&mut writer, &mut stderr, "blastx", 1, None)
                })
            }
            "16" => {
                let mut xml2 = super::blast_xml::xml2_writer(
                    Xml2Kind::Single16,
                    &opts,
                    &subjects,
                    &shown,
                    Box::new(NoSideFiles),
                );
                super::blast_xml::write_blastx_xml(
                    &mut writer,
                    &mut stderr,
                    &mut xml2,
                    &opts,
                    &queries,
                    &subjects,
                    &shown,
                    &[batch],
                    0,
                    false,
                    None,
                    None,
                )
                .and_then(|()| {
                    super::blast_xml::write_xml_epilog(
                        &mut writer,
                        &mut stderr,
                        "blastx",
                        1,
                        Some(&mut xml2),
                    )
                })
            }
            _ => {
                let spec = opts
                    .tabular_report()
                    .expect("7 qseqid sseqid qcovs has a TabularSpec");
                super::tabular::write_blastx_tabular(
                    &mut writer,
                    &mut stderr,
                    &opts,
                    spec,
                    &queries,
                    &subjects,
                    &shown,
                    &[batch],
                    None,
                    0,
                )
            }
        }
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
