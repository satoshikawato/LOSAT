//! `losat deltablast -query_loc` past NCBI BLAST+'s alignment of conserved domains (series PR, item PS13b;
//! PD-LOSAT-NCBI-DEFECTS approved exception 35, `Q-DL-QLOC=B`): where some kept domain hit of a range NCBI writes past
//! the array of its alignment (NCBI 2.17.0 crashes), LOSAT completes the search on the range. Its domain section
//! (`-show_domain_hits`) is rpsblast `-query_loc`'s, in the record's coordinates, and its rounds are those of the
//! range's letters as a query of their own, at the range's coordinates: the frozen NCBI run of those letters
//! (`db.dl.ql.or.dd300.sh.s6` of `tests/db_search_manifest/deltablast.tsv`) with the record's coordinates in its
//! domain rows.

use std::path::Path;
use std::process::Command;

/// The fields of the rows: none reads a query letter.
const FIELDS: &str =
    "6 qseqid sseqid length gapopen gaps qstart qend sstart send evalue bitscore score";

#[test]
fn a_range_past_ncbis_alignment_writes_the_domain_rows_of_the_range_and_the_rows_of_its_letters() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let output = Command::new(env!("CARGO_BIN_EXE_losat"))
        .current_dir(root)
        .env_remove("BLASTDB")
        .args([
            "deltablast",
            "-query",
            "tests/db_search/delta/q_dd.faa",
            "-query_loc",
            "300-1294",
            "-db",
            "tests/db_search/delta/dbs",
            "-rpsdb",
            "tests/db_search/delta/deltatest",
            "-show_domain_hits",
            "-outfmt",
            FIELDS,
        ])
        .output()
        .expect("runs losat");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let oracle =
        std::fs::read_to_string(root.join("tests/fixtures/db_search/db.dl.ql.or.dd300.sh.s6.out"))
            .expect("the frozen run of the range's letters");
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
    assert_eq!(String::from_utf8_lossy(&output.stdout), expected);
}
