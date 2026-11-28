//! psiblast `-phi_pattern` (PHI-BLAST) on subjects longer than NCBI's part (`MAX_DBSEQ_LEN`,
//! 5,000,000 residues), searched in NCBI's overlapping parts (overlap 100), and with
//! `-sum_stats true` (line PS, item PS8a-2).
//!
//! The inputs are generated deterministically (`write_inputs`): a 300-residue query whose pattern
//! `W-C-x(2)-H-W` occurs once (`phl_q.faa`) or twice, five letters apart (`phl_q2.faa`), and
//! subjects of letters without `W` (so the pattern occurs only where a copy of a query is put)
//! of 5,000,001 residues (a copy that starts the second part), 5,200,000 (a copy whose pattern is
//! in the overlap of the two parts, so that NCBI counts its hit in both), 5,100,000 (a pattern
//! across 5,000,000, in the second part only), 10,100,000 (three parts, patterns in both
//! overlaps) and 5,050,000 (the second query's two places in the overlap), with two short
//! subjects. They are searched with `-subject` and with a database built by
//! `losat makeblastdb`, at 1 and 4 threads. Expected rows in
//! `tests/fixtures/db_search_long/phi_long.*.6.out` were frozen from official NCBI BLAST+ 2.17.0
//! on the same inputs (`PHI_LONG_DIR=<dir>` keeps them there; README.md of that folder).
//!
//! The second test (`write_merge_inputs`, fix round 1 of PS8a-2) merges HSPs of two parts that
//! belong to different places of a pattern of variable length (`ZINC_FINGER`, eight places in
//! TF3A_YEAST): the merged HSP aligns the pattern region from one place's start
//! with another's length, so its Seq-align (the start plus the edit script) ends elsewhere than
//! its own ends: past the end of the query (`pz_q1`, NCBI's `qend` 442 of 429, shown with gaps
//! there) or inside it (`pz_q2`, the same query with 40 more residues that the subject does not
//! have). A third query (`pz_q3`) ends with a place of the pattern two residues shorter than
//! the place it merges with, so the merged HSP's pattern region runs two letters past its end,
//! into the zero bytes after the query in NCBI's buffer. Expected bytes in
//! `tests/fixtures/db_search_long/phi_merge.*.out` (outfmt 6 and 0, `-subject` and `-db`, with
//! and without `-sum_stats true`), frozen from official NCBI BLAST+ 2.17.0 (`PHI_MERGE_DIR=<dir>`
//! keeps the inputs).

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

/// The binary with no `LOSAT_*` variable that chooses a code path (`LOSAT_BLASTX_WORKER_LOG`
/// stays, as in `protein_db_hard_long.rs`).
fn losat() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command
}

/// The query: its pattern `WCAAHW` at 140.
const QUERY: &[u8] = b"IADGTVITSLEADGTNAQMEEYMITPLMMMQGYPANDGANYTRVRYSKPKIFKVRKESIS\
RVFTDVSFDVKQVIGTPDMRRQYPEAEGMLSYEKDTDAPTNVRKRLKTMKQQMAQAVPRA\
GKLKRGDKDQDSVQNKENYMWCAAHWKEGPGARLSDTDNGYMRSYLIYAGITIPQPYKNL\
DTVPPVGPVDYNENKPDTVQNRQLNVQASISEPQMMRNANIMYISDVDTTMFLVGYRNSL\
KKQPISFSIDDRPDTMNGIQKMPGSFMNQVPMNQKDVLVTKMMQRREPMVSLQLPNIPGL";

/// Letters without `W`.
const BACKGROUND: &[u8] = b"ACDEFGHIKLMNPQRSTVY";

/// `length` letters of `BACKGROUND` from a linear congruential generator seeded with `seed`
/// (bits 16 and up of each state).
fn background(seed: u32, length: usize) -> Vec<u8> {
    letters_of(BACKGROUND, seed, length)
}

/// `length` letters of `alphabet` from a linear congruential generator seeded with `seed` (bits
/// 16 and up of each state).
fn letters_of(alphabet: &[u8], mut seed: u32, length: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(length);
    for _ in 0..length {
        seed = seed.wrapping_mul(1_103_515_245).wrapping_add(12_345) & 0x7fff_ffff;
        out.push(alphabet[(seed >> 16) as usize % alphabet.len()]);
    }
    out
}

/// The query with its pattern twice, at 140 and 145 (`WCAAHWCAAHW`): 305 residues.
fn query2() -> Vec<u8> {
    let mut out = QUERY[..140].to_vec();
    out.extend_from_slice(b"WCAAHWCAAHW");
    out.extend_from_slice(&QUERY[146..]);
    out
}

fn record(seed: u32, length: usize, copies: &[(usize, &[u8])]) -> Vec<u8> {
    let mut seq = background(seed, length);
    for &(at, piece) in copies {
        seq[at..at + piece.len()].copy_from_slice(piece);
    }
    seq
}

fn write_fasta(path: &Path, records: &[(&str, &[u8])]) {
    let mut out = String::new();
    for (name, seq) in records {
        writeln!(out, ">{name}").unwrap();
        for line in seq.chunks(80) {
            out.push_str(std::str::from_utf8(line).unwrap());
            out.push('\n');
        }
    }
    std::fs::write(path, out).unwrap();
}

/// Writes `phl_q.faa`, `phl_q2.faa`, `phl_s.faa` and `phl.pat` into `dir`.
fn write_inputs(dir: &Path) {
    let q2 = query2();
    let subjects: Vec<(&str, Vec<u8>)> = vec![
        ("phl_a", record(21, 2_000, &[(1_000, QUERY)])),
        // The pattern at 4,999,900, where the second part starts.
        (
            "phl_5m1",
            record(
                13,
                5_000_001,
                &[(1_000_000, QUERY), (4_999_760, &QUERY[..241])],
            ),
        ),
        // The pattern at 4,999,940, in the overlap: one hit in each part.
        (
            "phl_ov",
            record(14, 5_200_000, &[(2_000_000, QUERY), (4_999_800, QUERY)]),
        ),
        // The pattern at 4,999,997 to 5,000,002: in the second part only.
        ("phl_c2", record(16, 5_100_000, &[(4_999_857, QUERY)])),
        // Three parts, patterns in both overlaps and in the middle part.
        (
            "phl_3",
            record(
                17,
                10_100_000,
                &[(4_999_800, QUERY), (7_000_000, QUERY), (9_999_700, QUERY)],
            ),
        ),
        // The second query's two places (4,999,950 and 4,999,955) in the overlap.
        ("phl_2", record(18, 5_050_000, &[(4_999_810, &q2)])),
        ("phl_b", record(22, 3_000, &[(500, &q2)])),
    ];
    let refs: Vec<(&str, &[u8])> = subjects
        .iter()
        .map(|(name, seq)| (*name, seq.as_slice()))
        .collect();
    write_fasta(&dir.join("phl_s.faa"), &refs);
    write_fasta(&dir.join("phl_q.faa"), &[("phl_q", QUERY)]);
    write_fasta(&dir.join("phl_q2.faa"), &[("phl_q2", &q2)]);
    std::fs::write(dir.join("phl.pat"), "PA  W-C-x(2)-H-W\n").unwrap();
}

fn expected(name: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/db_search_long")
        .join(format!("phi_long.{name}.6.out"));
    std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

const OUTFMT: &str =
    "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score nident";

struct RunDir {
    path: PathBuf,
    keep: bool,
}

impl Drop for RunDir {
    fn drop(&mut self) {
        if !self.keep {
            let _ = std::fs::remove_dir_all(&self.path);
        }
    }
}

#[test]
fn phi_pattern_on_subjects_over_5_000_000_residues_and_with_sum_statistics_match_ncbi() {
    let run_dir = match std::env::var_os("PHI_LONG_DIR") {
        Some(dir) => RunDir {
            path: PathBuf::from(dir),
            keep: true,
        },
        None => RunDir {
            path: std::env::temp_dir()
                .join(format!("losat-psiblast-phi-long-{}", std::process::id())),
            keep: false,
        },
    };
    let dir = run_dir.path.as_path();
    std::fs::create_dir_all(dir).unwrap();
    write_inputs(dir);
    let db = dir.join("phl").to_str().unwrap().to_string();
    let made = losat()
        .arg("makeblastdb")
        .arg("-in")
        .arg(dir.join("phl_s.faa"))
        .args(["-dbtype", "prot", "-out", &db])
        .output()
        .expect("run losat makeblastdb");
    assert!(
        made.status.success(),
        "makeblastdb: {}",
        String::from_utf8_lossy(&made.stderr)
    );
    let subject = dir.join("phl_s.faa").to_str().unwrap().to_string();
    let pattern = dir.join("phl.pat").to_str().unwrap().to_string();
    let q1 = dir.join("phl_q.faa").to_str().unwrap().to_string();
    let q2 = dir.join("phl_q2.faa").to_str().unwrap().to_string();
    // `-evalue 1e-138`: the linked e-values of `-sum_stats true` (the pattern table's Lambda,
    // `logK` 0) free the HSPs of scores 1,217 and 1,201 before the traceback, and HSPs whose sets
    // link HSPs of the other place of the second query; without it those HSPs stay.
    let cases: &[(&str, &str, [&str; 2], &[&str])] = &[
        ("sub", &q1, ["-subject", &subject], &[]),
        ("db", &q1, ["-db", &db], &[]),
        ("q2.sub", &q2, ["-subject", &subject], &[]),
        ("q2.db", &q2, ["-db", &db], &[]),
        (
            "e138.sub",
            &q1,
            ["-subject", &subject],
            &["-evalue", "1e-138"],
        ),
        (
            "ss.e138.sub",
            &q1,
            ["-subject", &subject],
            &["-evalue", "1e-138", "-sum_stats", "true"],
        ),
        (
            "ss.e138.q2.db",
            &q2,
            ["-db", &db],
            &["-evalue", "1e-138", "-sum_stats", "true"],
        ),
    ];
    let mut failed = Vec::new();
    for &(label, query, target, extra) in cases {
        let want = expected(label);
        for threads in [1, 4] {
            let output = losat()
                .args(["psiblast", "-query", query])
                .args(target)
                .args([
                    "-phi_pattern",
                    &pattern,
                    "-num_iterations",
                    "1",
                    "-outfmt",
                    OUTFMT,
                ])
                .args(extra)
                .args(["-num_threads", &threads.to_string()])
                .output()
                .expect("run losat psiblast");
            // Nothing on the standard error: LOSAT searches `-subject` with the threads asked for
            // and writes no line of NCBI's that it ignores them (docs/manual/src/differences.md).
            let stderr = String::from_utf8_lossy(&output.stderr);
            if !output.status.success() || !stderr.is_empty() {
                failed.push(format!(
                    "{label} t{threads}: status {:?}, stderr {stderr:?}",
                    output.status.code(),
                ));
            } else if output.stdout != want {
                failed.push(format!(
                    "{label} t{threads}: {} bytes vs expected {}:\n{}",
                    output.stdout.len(),
                    want.len(),
                    String::from_utf8_lossy(&output.stdout)
                ));
            }
        }
    }
    // Outfmt 8: NCBI's engine error (exit status 3, no output) once an alignment is found.
    let output = losat()
        .args(["psiblast", "-query", &q1, "-subject", &subject])
        .args([
            "-phi_pattern",
            &pattern,
            "-num_iterations",
            "1",
            "-outfmt",
            "8",
        ])
        .output()
        .expect("run losat psiblast");
    let stderr = String::from_utf8_lossy(&output.stderr);
    if output.status.code() != Some(3)
        || !output.stdout.is_empty()
        || stderr != "LOSAT engine error: Program type 'phiblastp' not supported\n"
    {
        failed.push(format!(
            "outfmt 8: status {:?}, {} bytes, stderr {stderr:?}",
            output.status.code(),
            output.stdout.len()
        ));
    }
    assert!(failed.is_empty(), "mismatches:\n{}", failed.join("\n"));
}

/// TF3A_YEAST (429 residues): eight places of the zinc-finger pattern `ZINC_FINGER`, of 21 to
/// 24 residues.
const TF3A: &[u8] = b"MGGEVLNNEGMPLAELKQETIPISRSESSESLNSLTSTRSSSSNRPKTYFCDYDGCDKAF\
TRPSILTEHQLSVHQGLRAFQCDKCAKSFVKKSHLERHLYTHSDTKPFQCSYCGKGVTTR\
QQLKRHEVTHTKSFICPEEGCNLRFYKHPQLRAHILSVHLHKLTCPHCNKSFQRPYRLRN\
HISKHHDPEVENPYQCTFAGCCKEFRIWSQLQSHIKNDHPKLKCPICSKPCVGENGLQMH\
MIIHDDSLVTKNWKCHICPDMSFSRKHDLLTHYGSIHTEEDIPLELKYKISDIQQLVQDH\
GVQLGNSKHSNEQDEEKISNRLRKRRKLTENNNVEFLQNEVDLEKRLESGENGLNLLLNT\
VGRKYRCFYNNCSRTFKTKEKYEKHIDKHKVHELKLKILQEKEENKTLVDQNHKEPFIIQ\
KETQSAGDK";

/// The residues after TF3A in the second query, without `C` or `H` (no new place).
const TF3A_TAIL: &[u8] = b"SNQLTPEVKAGDLSLRSMNEQAPKTIDSVLGNPEVKALAE";

const ZINC_FINGER: &str = "PA  C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H\n";

/// Letters without `C` and `H`, so that the pattern occurs only in the copy of a query.
const MERGE_BACKGROUND: &[u8] = b"ADEFGIKLMNPQRSTVWY";

/// Writes `pz_q1.faa` (TF3A), `pz_q2.faa` (TF3A and `TF3A_TAIL`), `pz_s.faa` (5,000,300 letters
/// of `MERGE_BACKGROUND`, TF3A put at 4,999,771, across the end of the first part, so that HSPs of
/// its places in the first part and in the second merge), `pz_q3.faa` (TF3A's first 81 residues,
/// with its place at 50 of 24 residues, then its place at 164 of 22, which ends the query),
/// `pz_r.faa` (as `pz_s`, with `pz_q3` put at 4,999,909) and `zf.pat` into `dir`. Some merged
/// HSPs' Seq-aligns run past their own ends: past the end of `pz_q1`, inside `pz_q2`, whose
/// letters after TF3A the subject does not have. `pz_q3`'s merged HSP aligns its pattern region
/// from its last place's start with the first place's length, two letters past its end.
fn write_merge_inputs(dir: &Path) {
    let q2 = [TF3A, TF3A_TAIL].concat();
    let q3 = [&TF3A[..81], &TF3A[164..186]].concat();
    write_fasta(&dir.join("pz_q1.faa"), &[("pz_q1", TF3A)]);
    write_fasta(&dir.join("pz_q2.faa"), &[("pz_q2", &q2[..])]);
    write_fasta(&dir.join("pz_q3.faa"), &[("pz_q3", &q3[..])]);
    for (name, at, copy) in [("pz_s", 4_999_771, TF3A), ("pz_r", 4_999_909, &q3[..])] {
        let mut subject = letters_of(MERGE_BACKGROUND, 31, 5_000_300);
        subject[at..at + copy.len()].copy_from_slice(copy);
        write_fasta(&dir.join(format!("{name}.faa")), &[(name, &subject[..])]);
    }
    std::fs::write(dir.join("zf.pat"), ZINC_FINGER).unwrap();
}

/// `output` with its line of the database's date (`Posted date:` of outfmt 0's footer, the date
/// the database was made) left out.
fn without_posted_date(output: &[u8]) -> Vec<u8> {
    output
        .split_inclusive(|&byte| byte == b'\n')
        .filter(|line| !line.trim_ascii_start().starts_with(b"Posted date:"))
        .flatten()
        .copied()
        .collect()
}

const MERGE_OUTFMT: &str = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart \
send evalue bitscore score nident positive gaps qcovhsp qcovs";

#[test]
fn phi_hsps_merged_across_places_of_the_pattern_match_ncbi() {
    let run_dir = match std::env::var_os("PHI_MERGE_DIR") {
        Some(dir) => RunDir {
            path: PathBuf::from(dir),
            keep: true,
        },
        None => RunDir {
            path: std::env::temp_dir()
                .join(format!("losat-psiblast-phi-merge-{}", std::process::id())),
            keep: false,
        },
    };
    let dir = run_dir.path.as_path();
    std::fs::create_dir_all(dir).unwrap();
    write_merge_inputs(dir);
    // The reports name the inputs as given (outfmt 0's database title): relative names in `dir`.
    for (input, db) in [("pz_s.faa", "pz"), ("pz_r.faa", "pr")] {
        let made = losat()
            .current_dir(dir)
            .args(["makeblastdb", "-in", input, "-dbtype", "prot", "-out", db])
            .output()
            .expect("run losat makeblastdb");
        assert!(
            made.status.success(),
            "makeblastdb {input}: {}",
            String::from_utf8_lossy(&made.stderr)
        );
    }
    const NONE: &[&str] = &[];
    const SUM_STATS: &[&str] = &["-sum_stats", "true"];
    // Outfmt 0 keeps the alignments of the best E-values (the frozen bytes stay under 64 KB),
    // those of the merged HSPs among them.
    const BEST: &[&str] = &["-evalue", "1e-180"];
    const BEST_SETS: &[&str] = &["-evalue", "1e-100", "-sum_stats", "true"];
    let mut cases: Vec<(String, String, [&str; 2], &str, &[&str])> = Vec::new();
    for n in [1, 2] {
        let (sub, db) = (["-subject", "pz_s.faa"], ["-db", "pz"]);
        for (label, target, outfmt, extra) in [
            ("sub.6", sub, MERGE_OUTFMT, NONE),
            ("db.6", db, MERGE_OUTFMT, NONE),
            ("ss.sub.6", sub, MERGE_OUTFMT, SUM_STATS),
            ("ss.db.6", db, MERGE_OUTFMT, SUM_STATS),
            ("sub.0", sub, "0", BEST),
            ("ss.db.0", db, "0", BEST_SETS),
        ] {
            cases.push((
                format!("q{n}.{label}"),
                format!("pz_q{n}.faa"),
                target,
                outfmt,
                extra,
            ));
        }
    }
    // The third query's merged HSP aligns its pattern region two letters past the query's end:
    // NCBI reads the zero bytes after the query's letters there.
    let (sub, db) = (["-subject", "pz_r.faa"], ["-db", "pr"]);
    for (label, target, outfmt, extra) in [
        ("sub.6", sub, MERGE_OUTFMT, NONE),
        ("db.6", db, MERGE_OUTFMT, NONE),
        ("ss.sub.6", sub, MERGE_OUTFMT, SUM_STATS),
        ("sub.0", sub, "0", NONE),
    ] {
        cases.push((
            format!("q3.{label}"),
            "pz_q3.faa".to_string(),
            target,
            outfmt,
            extra,
        ));
    }
    let mut failed = Vec::new();
    for (label, query, target, outfmt, extra) in cases {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/db_search_long")
            .join(format!("phi_merge.{label}.out"));
        let frozen =
            std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
        // NCBI's words in LOSAT's (outfmt 0's program and version lines).
        let want = without_posted_date(&losat::wording::losat_output_bytes(&frozen));
        for threads in [1, 4] {
            let output = losat()
                .current_dir(dir)
                .args(["psiblast", "-query", &query])
                .args(target)
                .args(["-phi_pattern", "zf.pat", "-num_iterations", "1"])
                .args(["-outfmt", outfmt])
                .args(extra)
                .args(["-num_threads", &threads.to_string()])
                .output()
                .expect("run losat psiblast");
            let stderr = String::from_utf8_lossy(&output.stderr);
            if !output.status.success() || !stderr.is_empty() {
                failed.push(format!(
                    "{label} t{threads}: status {:?}, stderr {stderr:?}",
                    output.status.code(),
                ));
            } else if without_posted_date(&output.stdout) != want {
                failed.push(format!(
                    "{label} t{threads}: {} bytes vs expected {}:\n{}",
                    output.stdout.len(),
                    want.len(),
                    String::from_utf8_lossy(&output.stdout)
                ));
            }
        }
    }
    assert!(failed.is_empty(), "mismatches:\n{}", failed.join("\n"));
}
