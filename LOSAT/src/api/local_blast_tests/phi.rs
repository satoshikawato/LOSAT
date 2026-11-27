//! psiblast's PHI-BLAST search (`-phi_pattern`; line PS item PS8b) in process: the HSP numbers
//! that its formats give the observer.

use super::support::{assert_tabular_ranges, field, reader_records, run_formats};
use crate::api::local_blast::run_local_database_blastp;
use crate::api::QueryEnd;
use crate::blastinput::fasta_reader::ReaderConfig;
use crate::blastinput::ncbi_environment::ApplicationSettings;
use crate::cli::{try_parse_from, Cli, Commands};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3572-3632
// The search gives a result set per pattern occurrence, and the call's hit list holds them one
// after the other, the query batches' in order (`ps8_q3.faa`: three batches, the second with three
// occurrences). A batch's report and an occurrence's table number their HSPs from 0; the observer
// gets each HSP by its index in the hit list (`FormatProbe::from_hsp`).
#[test]
fn a_phi_search_reports_each_hsp_by_its_index_in_the_hit_list() {
    let manifest = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let query = manifest.join("tests/fasta/psiblast/ps8_q3.faa");
    let [query_name, database, pattern] = [
        &query,
        &manifest.join("tests/db_search/paja"),
        &manifest.join("tests/fixtures/phi/c2c.pat"),
    ]
    .map(|path| path.display().to_string());
    let records = reader_records(&query, ReaderConfig::query("PSIBLAST", true, true));
    let run = run_formats(&["6", "7", "0"], true, |outputs| {
        let argv = [
            "losat",
            "psiblast",
            "-query",
            &query_name,
            "-db",
            &database,
            "-phi_pattern",
            &pattern,
        ];
        let cli: Cli = try_parse_from(argv).expect("argv");
        let Commands::Psiblast(args) = cli.command else {
            unreachable!("psiblast's arguments")
        };
        run_local_database_blastp(
            args.into_blastp_args()?,
            &records,
            &QueryEnd::Input,
            ApplicationSettings::default(),
            outputs,
        )
    });
    let hits = &run.hits;
    assert_eq!(
        hits.len(),
        66,
        "the HSPs of the three queries' five occurrences"
    );
    // outfmt 6 and 7: the rows in hit-list order, with the coordinates of their HSP.
    let rows = assert_tabular_ranges(&run, 0, Some(1), "psiblast -phi_pattern");
    assert_eq!(rows, hits.len(), "outfmt 6 rows");
    for &(hsp, start, end) in &run.ranges[0] {
        let (hit, row) = (&hits[hsp].hit, &run.outputs[0][start..end]);
        assert_eq!(
            [6, 7, 8, 9].map(|index| field(row, index)),
            [hit.q_start, hit.q_end, hit.s_start, hit.s_end].map(|value| value.to_string()),
            "outfmt 6 HSP {hsp} coordinates"
        );
    }
    // outfmt 0: the section of each HSP once, with its raw score.
    let mut seen = vec![false; hits.len()];
    for &(hsp, start, end) in &run.ranges[2] {
        assert!(
            !std::mem::replace(&mut seen[hsp], true),
            "outfmt 0: HSP {hsp} twice"
        );
        let section = std::str::from_utf8(&run.outputs[2][start..end]).expect("utf-8 section");
        assert!(
            section.contains(&format!("bits ({}),", hits[hsp].hit.raw_score)),
            "outfmt 0: the raw score of HSP {hsp}"
        );
    }
    assert!(seen.iter().all(|&seen| seen), "outfmt 0: a section per HSP");
}
