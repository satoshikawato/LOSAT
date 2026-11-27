//! `losat deltablast -query_loc` past NCBI BLAST+'s alignment of conserved domains (series PR, item PS13b;
//! PD-LOSAT-NCBI-DEFECTS approved exception 35, `Q-DL-QLOC=B`): where some kept domain hit of a range NCBI writes past
//! the array of its alignment (NCBI 2.17.0 crashes), LOSAT completes the search on the range. Its domain section
//! (`-show_domain_hits`) is rpsblast `-query_loc`'s, in the record's coordinates, and its rounds are those of the
//! range's letters as a query of their own, at the range's coordinates: the frozen NCBI run of those letters
//! (`db.dl.ql.or.dd300.sh.s6` of `tests/db_search_manifest/deltablast.tsv`) with the record's coordinates in its
//! domain rows. The rule is each query's own: in a batch with such a range, a query whose cells stay inside NCBI's
//! array keeps NCBI's placement (`db.dl.ql.tail5.f6`, NCBI's run of that query alone with the range).

use std::path::Path;
use std::process::Command;

/// The fields of the rows: none reads a query letter.
const FIELDS: &str =
    "6 qseqid sseqid length gapopen gaps qstart qend sstart send evalue bitscore score";

/// The run of `losat deltablast` on `query` against the fixture's `dbs` and `deltatest` with `args`: its
/// standard output, after checking that it succeeds.
fn deltablast(root: &Path, query: &str, args: &[&str]) -> String {
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .current_dir(root)
        .env_remove("BLASTDB")
        .args(["deltablast", "-query", query])
        .args([
            "-db",
            "tests/db_search/delta/dbs",
            "-rpsdb",
            "tests/db_search/delta/deltatest",
        ])
        .args(args)
        .output()
        .expect("runs losat");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).expect("tabular rows are text")
}

/// The frozen standard output of the row `row` (`tests/fixtures/db_search`).
fn frozen(root: &Path, row: &str) -> String {
    std::fs::read_to_string(root.join(format!("tests/fixtures/db_search/{row}.out")))
        .expect("a frozen NCBI run")
}

#[test]
fn a_range_past_ncbis_alignment_writes_the_domain_rows_of_the_range_and_the_rows_of_its_letters() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let got = deltablast(
        root,
        "tests/db_search/delta/q_dd.faa",
        &[
            "-query_loc",
            "300-1294",
            "-show_domain_hits",
            "-outfmt",
            FIELDS,
        ],
    );
    let oracle = frozen(root, "db.dl.ql.or.dd300.sh.s6");
    // The domain rows (subjects `gnl|CDD|`) in the record's coordinates: the range starts at 300.
    let mut domain_rows = 0;
    let expected: String = oracle
        .lines()
        .map(|line| {
            let mut fields: Vec<String> = line.split('\t').map(str::to_string).collect();
            if fields[1].starts_with("gnl|CDD|") {
                domain_rows += 1;
                for field in &mut fields[5..7] {
                    *field = (field.parse::<u32>().expect("a query position") + 299).to_string();
                }
            }
            fields.join("\t") + "\n"
        })
        .collect();
    assert!(domain_rows > 0, "the oracle has domain rows");
    assert!(
        oracle.lines().count() > domain_rows,
        "the oracle has round rows"
    );
    assert_eq!(got, expected);
}

#[test]
fn a_batch_of_ranges_keeps_each_querys_placement() {
    // `q_three.faa` with the range 5-1094: qserp and YP_009220615.1, whose last kept domain hits end past the
    // range's length in their records (NCBI 2.17.0 crashes on the batch), then qtail, whose cells stay in their
    // rows (NCBI completes qtail alone).
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let got = deltablast(
        root,
        "tests/db_search/delta/q_three.faa",
        &["-query_loc", "5-1094", "-outfmt", FIELDS],
    );
    let is_qtail = |line: &&str| line.starts_with("qtail\t");
    // qserp and YP_009220615.1: NCBI's run of the ranges' letters as their own queries.
    let oracle = frozen(root, "db.dl.ql.or.three5.s6");
    let others: String = oracle
        .lines()
        .filter(|line| !is_qtail(line))
        .map(|line| format!("{line}\n"))
        .collect();
    // qtail: NCBI's run of the query alone with the range, whose fields hold these (`F6` of the manifest).
    let qtail: String = frozen(root, "db.dl.ql.tail5.f6")
        .lines()
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            [0, 1, 8, 9, 10, 4, 5, 6, 7, 11, 12, 13]
                .map(|index| fields[index])
                .join("\t")
                + "\n"
        })
        .collect();
    let own_query_qtail: String = oracle
        .lines()
        .filter(is_qtail)
        .map(|line| format!("{line}\n"))
        .collect();
    assert!(!qtail.is_empty(), "NCBI's run of qtail has rows");
    assert_ne!(
        qtail, own_query_qtail,
        "qtail's placement (cells shifted in their rows) is not that of its range's letters as a query"
    );
    assert_eq!(got, others + &qtail);
}
