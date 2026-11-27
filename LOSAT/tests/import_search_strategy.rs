//! `-import_search_strategy FILE` of blastn, blastp, blastx, tblastn and tblastx (line RG, item
//! F2-2) against the bytes that NCBI BLAST+ 2.17.0 `-num_threads 1` writes: the standard output,
//! the standard error and the exit status of each row of
//! `tests/fixtures/strategy/import/import_cases.tsv` (its header says what the columns mean; the
//! strategy files are NCBI's own exports, `ss/`, and each output is stored once, `bytes/`). The rows are the RG-F2 probes of these programs
//! whose LOSAT result equals NCBI's: imports with `-db` and `-subject` strategies, the command
//! line's overrides (`-query`, `-db`, `-subject`, the formatting options), the warnings of the
//! options a strategy ignores, the strategies' query ranges and filters, the archive (outfmt 11)
//! of an imported run, the reader's errors, another program's strategy with options of the
//! command line (`C.*.from.*`), the parameters that the setters take (`S1.*`: best hits, database
//! masks, the repeat filter, the MegaBLAST index), the command line's index and database mask
//! over a strategy's (`B.*`), and the formatting options that a strategy's hit list size above
//! 500 leaves unread in the pairwise formats (`Y.*`), and the values of hand-edited strategies
//! outside the command line's ranges (`S3.*`: a negative window, the HSP filtering check, the
//! audit of the setters of F2's fix 3). The SAM documents (`-outfmt 17`) of a TBLASTX strategy that
//! blastn runs are compared without the CL tag of their `@PG` line, which holds the program's
//! own command line (NCBI's is removed from the frozen bytes, `F/f2/tools/freeze_b3.py`).
//!
//! A strategy's subjects are searched once whatever `-num_threads` says (Owner 2026-11-20 F2-Q1
//! "Q1=C"; NCBI prints every alignment once per thread or aborts): the runs at 2, 4 and 8 threads
//! write the frozen `-num_threads 1` bytes.

use std::path::PathBuf;
use std::process::{Command, Output, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_losat");

fn engine() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn fixtures() -> PathBuf {
    engine().join("tests/fixtures/strategy/import")
}

struct Case {
    name: String,
    program: String,
    args: Vec<String>,
    rc: i32,
    out: String,
    err: String,
}

fn cases() -> Vec<Case> {
    let table =
        std::fs::read_to_string(fixtures().join("import_cases.tsv")).expect("import_cases.tsv");
    table
        .lines()
        .filter(|line| !(line.starts_with('#') || line.starts_with("name\t") || line.is_empty()))
        .map(|line| {
            let cells: Vec<&str> = line.split('\t').collect();
            assert_eq!(cells.len(), 6, "import_cases.tsv row: {line}");
            Case {
                name: cells[0].to_string(),
                program: cells[1].to_string(),
                args: cells[2]
                    .split(' ')
                    .map(|arg| arg.replace("%20", " "))
                    .collect(),
                rc: cells[3].parse().expect("rc"),
                out: cells[4].to_string(),
                err: cells[5].to_string(),
            }
        })
        .collect()
}

/// NCBI's bytes of the file `name` of `bytes/` (`-`: none) in LOSAT's words.
fn frozen(name: &str) -> Vec<u8> {
    if name == "-" {
        return Vec::new();
    }
    crate::outfmt11_blastx::expected_losat_bytes(
        std::fs::read(fixtures().join("bytes").join(name)).expect("frozen bytes"),
    )
}

/// One run of the case from `LOSAT/` without LOSAT's switches and batch variables; standard
/// input is empty.
fn run(case: &Case, threads: Option<&str>) -> Output {
    let ss = fixtures().join("ss");
    let ss = ss.to_str().expect("UTF-8 path");
    let mut command = Command::new(BIN);
    command
        .current_dir(engine())
        .arg(&case.program)
        .args(case.args.iter().map(|arg| arg.replace("{ss}", ss)))
        .stdin(Stdio::null());
    if let Some(threads) = threads {
        command.args(["-num_threads", threads]);
    }
    for (name, _) in std::env::vars() {
        if name.starts_with("LOSAT_X_")
            || matches!(
                name.as_str(),
                "ARCHIVE_FORMAT" | "BATCH_SIZE" | "CHUNK_SIZE" | "BLASTDB"
            )
        {
            command.env_remove(name);
        }
    }
    command.output().expect("run losat")
}

/// The SAM document `document` without the CL tag of its `@PG` line: the bytes from the `\tCL:`
/// after `\n@PG\tID:0\tVN:` up to the line's `\tPN:` tag (the rule of `freeze_b3.py`
/// `without_sam_cl`).
fn without_sam_cl(document: &[u8]) -> Vec<u8> {
    let find = |needle: &[u8], from: usize| {
        document
            .get(from..)?
            .windows(needle.len())
            .position(|window| window == needle)
            .map(|at| at + from)
    };
    let Some(pg) = find(b"\n@PG\tID:0\tVN:", 0) else {
        return document.to_vec();
    };
    let line_end = find(b"\n", pg + 1).unwrap_or(document.len());
    let Some(cl) = find(b"\tCL:", pg + 1).filter(|&at| at < line_end) else {
        return document.to_vec();
    };
    let Some(pn) = find(b"\tPN:", cl).filter(|&at| at < line_end) else {
        return document.to_vec();
    };
    [&document[..cl], &document[pn..]].concat()
}

fn check(case: &Case, threads: Option<&str>) {
    let mut output = run(case, threads);
    let sam = case
        .args
        .iter()
        .position(|arg| arg == "-outfmt")
        .and_then(|at| case.args.get(at + 1))
        .is_some_and(|format| format.starts_with("17"));
    if sam {
        output.stdout = without_sam_cl(&output.stdout);
    }
    let label = format!("{} threads {threads:?}", case.name);
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert_eq!(
        output.status.code(),
        Some(case.rc),
        "{label}: exit status; stderr: {stderr}"
    );
    let want_out = frozen(&case.out);
    assert!(
        output.stdout == want_out,
        "{label}: stdout differs ({} bytes, NCBI wrote {}): {}",
        output.stdout.len(),
        want_out.len(),
        crate::outfmt11_blastx::first_difference(&output.stdout, &want_out)
    );
    let stderr = stderr.replace(fixtures().join("ss").to_str().unwrap(), "{SS}");
    assert_eq!(
        stderr,
        String::from_utf8_lossy(&frozen(&case.err)),
        "{label}: stderr"
    );
}

#[test]
fn imports_equal_ncbis() {
    let cases = cases();
    assert!(cases.len() >= 100, "import_cases.tsv rows");
    let failed: Vec<&str> = cases
        .iter()
        .filter(|case| std::panic::catch_unwind(|| check(case, None)).is_err())
        .map(|case| case.name.as_str())
        .collect();
    assert!(failed.is_empty(), "rows that differ: {failed:?}");
}

/// F2-Q1 (Owner 2026-11-20 "Q1=C"): with the strategy's subjects, `-num_threads` 2, 4 and 8
/// write NCBI's `-num_threads 1` bytes in each program (NCBI prints each alignment once per
/// thread, or aborts), also where blastn runs a TBLASTX strategy.
#[test]
fn a_strategys_subjects_are_searched_once_at_any_thread_count() {
    let cases = cases();
    let names = ["blastn", "blastp", "blastx", "tblastn", "tblastx"]
        .map(|program| format!("I.{program}.sub.alone"));
    for name in names
        .iter()
        .map(String::as_str)
        .chain(["C.blastn.from.tblastx.sub"])
    {
        let case = cases
            .iter()
            .find(|case| case.name == name)
            .unwrap_or_else(|| panic!("row {name}"));
        for threads in ["2", "4", "8"] {
            check(case, Some(threads));
        }
    }
}

/// F2 fix 3 (R3-B2): the values of a hand-edited strategy outside the command line's ranges that
/// LOSAT does not run where NCBI crashes (a genetic code without a table, composition-based
/// statistics 4, a greedy extension with dynamic programming's traceback) or runs a search that
/// LOSAT does not have (a word size of 1, unified P-values in TBLASTN, composition mode -1, the
/// refusals of the command line's values): LOSAT's refusal in its words, exit status 1, never a
/// panic (`docs/manual/src/refusals.md`). Refusals checked where the search starts (BLASTP
/// `-word_size 1`, BLASTP/TBLASTN infinite `-evalue`) follow NCBI's outfmt 0 prolog written before
/// the first query batch (T4b-f, `I.blastp.db.alone` / `I.tblastn.db.alone`). The strategies are
/// those of the audit, `ss/s3a.*`.
#[test]
fn a_strategy_value_that_losat_does_not_run_is_refused_in_losat_words() {
    let bp_prolog = Some("39de21324836dfef");
    let tn_prolog = Some("19f948e973ed7555");
    let refusals: &[(&str, &str, &str, Option<&str>)] = &[
        ("blastx", "blastx_QueryGeneticCode_0", "a search strategy with this query genetic code is not supported by LOSAT's BLASTX", None),
        ("blastx", "blastx_QueryGeneticCode_7", "a search strategy with this query genetic code is not supported by LOSAT's BLASTX", None),
        ("blastx", "blastx_QueryGeneticCode_256", "a search strategy with this query genetic code is not supported by LOSAT's BLASTX", None),
        ("tblastn", "tblastn_DbGeneticCode_0", "a search strategy with this database genetic code is not supported by LOSAT's TBLASTN", None),
        ("tblastn", "tblastn_DbGeneticCode_7", "a search strategy with this database genetic code is not supported by LOSAT's TBLASTN", None),
        ("tblastn", "tblastn_DbGeneticCode_-1", "a search strategy with this database genetic code is not supported by LOSAT's TBLASTN", None),
        ("tblastx", "tblastx_QueryGeneticCode_7", "a search strategy with this query genetic code is not supported by LOSAT's TBLASTX", None),
        ("tblastx", "tblastx_QueryGeneticCode_-1", "a search strategy with this query genetic code is not supported by LOSAT's TBLASTX", None),
        ("tblastx", "tblastx_DbGeneticCode_7", "a search strategy with this database genetic code is not supported by LOSAT's TBLASTX", None),
        ("blastp", "blastp_CompositionBasedStats_4", "a search strategy with composition-based statistics 4 (unified P 0) is not supported by LOSAT's BLASTP", None),
        ("blastp", "blastp_CompositionBasedStats_neg1", "a search strategy with composition-based statistics -1 (unified P 0) is not supported by LOSAT's BLASTP", None),
        ("tblastn", "tblastn_CompositionBasedStats_4", "a search strategy with composition-based statistics 4 (unified P 0) is not supported by LOSAT's TBLASTN", None),
        ("tblastn", "tblastn_UnifiedP_2", "a search strategy with composition-based statistics 2 (unified P 2) is not supported by LOSAT's TBLASTN", None),
        ("blastp", "blastp_WordSize_1", "-word_size 1 (the lookup table) is not supported by LOSAT's BLASTP", bp_prolog),
        ("blastx", "blastx_WordSize_1", "-word_size 1 (the lookup table) is not supported by LOSAT's BLASTX", None),
        ("tblastn", "tblastn_WordSize_1", "the matrix BLOSUM62 with word size 1 is not supported by LOSAT's TBLASTN", None),
        ("blastn", "blastn_GapExtnAlgorithm_0_1", "a search strategy with this pair of gapped extension algorithms is not supported by LOSAT's BLASTN", None),
        ("blastn", "blastn_MatchReward_0", "is not supported by LOSAT's BLASTN", None),
        ("blastp", "blastp_EvalueThreshold_inf", "is not supported by LOSAT's BLASTP", bp_prolog),
        ("tblastn", "tblastn_EvalueThreshold_inf", "is not supported by LOSAT's TBLASTN", tn_prolog),
        ("blastx", "blastx_EvalueThreshold_inf", "is not supported by LOSAT's BLASTX", None),
        ("blastx", "blastx_EvalueThreshold_nan", "is not supported by LOSAT's BLASTX", None),
    ];
    let mut failed = Vec::new();
    for (program, strategy, words, expected_prolog) in refusals {
        let case = Case {
            name: format!("refusal {program} {strategy}"),
            program: program.to_string(),
            args: vec![
                "-import_search_strategy".to_string(),
                format!("{{ss}}/s3a.{strategy}.ss"),
            ],
            rc: 1,
            out: "-".to_string(),
            err: "-".to_string(),
        };
        let output = run(&case, None);
        let stderr = String::from_utf8_lossy(&output.stderr);
        let expected_stdout = match expected_prolog {
            Some(name) => {
                // The prolog: the bytes of NCBI's run before its first query's report.
                let bytes = frozen(name);
                let end = bytes
                    .windows(8)
                    .position(|window| window == b"\n\nQuery=")
                    .expect("the first query's report");
                bytes[..end].to_vec()
            }
            None => Vec::new(),
        };
        if output.status.code() != Some(1)
            || !stderr.contains(words)
            || output.stdout != expected_stdout
        {
            failed.push(format!(
                "{}: exit {:?}, stdout_len {} (expected {}), stderr {stderr:?}",
                case.name,
                output.status.code(),
                output.stdout.len(),
                expected_stdout.len(),
            ));
        }
    }
    assert!(failed.is_empty(), "{failed:#?}");
}

struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "losat_import_search_strategy_{tag}_{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("create temp dir");
        Self(dir)
    }

    fn path(&self) -> &std::path::Path {
        &self.0
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// PR #1640 review 4 (5.B1): BLASTN with a negative `WindowSize` on a database containing subjects
/// shorter than the word size (`dbA`: s1, ten 5-base subjects, s3; `dbB`: s1, s3). Each non-empty
/// short subject must still exit the carried diagonal table (moving its offset by length plus
/// window), blocking s3 in `dbA` (1 row) while `dbB` reports both s3 and s1 (2 rows), matching
/// NCBI BLAST+ 2.17.0.
#[test]
fn blastn_negative_window_short_subjects_in_database_match_ncbi() {
    let dir = TempDir::new("bn_short_subj");
    let ss = fixtures().join("ss/s4.blastn.short_subj.ss");
    let want_a = frozen(
        cases()
            .iter()
            .find(|c| c.name == "S4.blastn.short_subj_A")
            .expect("S4.blastn.short_subj_A")
            .out
            .as_str(),
    );
    let want_b = frozen(
        cases()
            .iter()
            .find(|c| c.name == "S4.blastn.short_subj_B")
            .expect("S4.blastn.short_subj_B")
            .out
            .as_str(),
    );
    assert_ne!(want_a, want_b, "A must block s3 while B reports s3 and s1");
    for (label, fasta, want) in [
        ("dbA", "short_subj_A.fna", &want_a),
        ("dbB", "short_subj_B.fna", &want_b),
    ] {
        let db = dir.path().join(label);
        let mk = Command::new(BIN)
            .current_dir(engine())
            .args([
                "makeblastdb",
                "-in",
                &format!("tests/fixtures/strategy/import/in/{fasta}"),
                "-dbtype",
                "nucl",
                "-out",
                db.to_str().expect("utf8"),
            ])
            .stdin(Stdio::null())
            .output()
            .expect("makeblastdb");
        assert!(mk.status.success(), "{label} makeblastdb: {mk:?}");
        for threads in ["1", "4"] {
            let out = Command::new(BIN)
                .current_dir(engine())
                .args([
                    "blastn",
                    "-import_search_strategy",
                    ss.to_str().expect("utf8"),
                    "-db",
                    db.to_str().expect("utf8"),
                    "-outfmt",
                    "6",
                    "-num_threads",
                    threads,
                ])
                .stdin(Stdio::null())
                .output()
                .expect("blastn");
            assert_eq!(out.status.code(), Some(0), "{label} t={threads}: {out:?}");
            assert_eq!(
                out.stdout, *want,
                "{label} t={threads}: stdout differs from NCBI 2.17.0"
            );
        }
    }
}

/// PR #1640 review 4 (5.S1): BLASTP with a negative `WindowSize` (-40) over a subject longer than
/// `MAX_DBSEQ_LEN` (5,000,500 residues, with motifs in chunk 0, across the 4,999,900..5,000,000
/// overlap, and in chunk 1, followed by a second subject `s2`) matches NCBI BLAST+ 2.17.0 bytes
/// for both `-subject` and `-db` at 1 and 4 threads.
#[test]
fn blastp_negative_window_searches_subject_over_5m_residues_as_ncbi() {
    let dir = TempDir::new("bp_long5m");
    let motif = b"MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLMKYEFSSDIHPCNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHKYPKCMCVMVRPTDQSLRHWKIMHDANFEKLKIFHYMVGHGRTASEG";
    let mut seq = vec![b'A'; 5_000_500];
    for pos in [100usize, 4_999_880, 5_000_200] {
        seq[pos..pos + motif.len()].copy_from_slice(motif);
    }
    let fasta = dir.path().join("long5m.faa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta).expect("create long5m.faa");
        f.write_all(b">long1 5000500 aa\n").expect("write");
        f.write_all(&seq).expect("write");
        f.write_all(b"\n>s2 short after long\n").expect("write");
        f.write_all(motif).expect("write");
        f.write_all(b"\n").expect("write");
    }
    let db = dir.path().join("long5m");
    let mk = Command::new(BIN)
        .current_dir(engine())
        .args([
            "makeblastdb",
            "-in",
            fasta.to_str().expect("utf8"),
            "-dbtype",
            "prot",
            "-out",
            db.to_str().expect("utf8"),
        ])
        .stdin(Stdio::null())
        .output()
        .expect("makeblastdb");
    assert!(mk.status.success(), "makeblastdb: {mk:?}");

    // Frozen NCBI BLAST+ 2.17.0 `-num_threads 1` bytes (identical for `-subject` and `-db`).
    let want_out = concat!(
        "lcl|Query_1\tlong1\t100.000\t140\t0\t0\t1\t140\t101\t240\t3.92e-109\t338\n",
        "lcl|Query_1\tlong1\t100.000\t140\t0\t0\t1\t140\t4999881\t5000020\t3.92e-109\t338\n",
        "lcl|Query_1\tlong1\t100.000\t140\t0\t0\t1\t140\t5000201\t5000340\t3.92e-109\t338\n",
        "lcl|Query_1\ts2\t100.000\t140\t0\t0\t1\t140\t1\t140\t1.53e-104\t293\n",
    );
    let want_err = concat!(
        "Warning: [blastp] Overriding database/subject in saved strategy\n",
        "Warning: [blastp] mt_mode cannot be overridden when using a search strategy\n",
    );
    let ss = fixtures().join("ss/s3.blastp.window_neg40.ss");
    for (flag, target) in [
        ("-subject", fasta.to_str().expect("utf8")),
        ("-db", db.to_str().expect("utf8")),
    ] {
        for threads in ["1", "4"] {
            let out = Command::new(BIN)
                .current_dir(engine())
                .args([
                    "blastp",
                    "-import_search_strategy",
                    ss.to_str().expect("utf8"),
                    flag,
                    target,
                    "-outfmt",
                    "6",
                    "-num_threads",
                    threads,
                ])
                .stdin(Stdio::null())
                .output()
                .expect("blastp");
            assert_eq!(out.status.code(), Some(0), "{flag} t={threads}: {out:?}");
            assert_eq!(
                String::from_utf8_lossy(&out.stdout),
                want_out,
                "{flag} t={threads} stdout"
            );
            assert_eq!(
                String::from_utf8_lossy(&out.stderr),
                want_err,
                "{flag} t={threads} stderr"
            );
        }
    }
}

/// PR #1640 R2 blocker: TBLASTN `-db_hard_mask 11` over a database containing an unmasked
/// 15,600,003-bp nucleotide sequence (translated frame > 5,000,000 residues, with a 140-aa query
/// encoded across frame +1's 4,999,900..5,000,000 chunk overlap) and a dust-masked record must
/// merge the long frame's chunks as NCBI BLAST+ 2.17.0 at 1 and 4 threads rather than refusing it
/// as an overlapping-mask merge.
#[test]
fn tblastn_db_hard_mask_merges_frame_over_5m_residues_as_ncbi() {
    // MT19937 seeded with 1640 matching Python's choice over ACGT.
    struct Mt19937 {
        mt: [u32; 624],
        idx: usize,
    }
    impl Mt19937 {
        fn seed(key: u32) -> Self {
            let mut mt = [0u32; 624];
            mt[0] = 19_650_218;
            for mti in 1..624 {
                mt[mti] = 1_812_433_253u32
                    .wrapping_mul(mt[mti - 1] ^ (mt[mti - 1] >> 30))
                    .wrapping_add(mti as u32);
            }
            let (mut i, mut j) = (1usize, 0usize);
            for _ in 0..624 {
                mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1_664_525)))
                    .wrapping_add(key)
                    .wrapping_add(j as u32);
                i += 1;
                j += 1;
                if i >= 624 {
                    mt[0] = mt[623];
                    i = 1;
                }
                if j >= 1 {
                    j = 0;
                }
            }
            for _ in 0..623 {
                mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)).wrapping_mul(1_566_083_941)))
                    .wrapping_sub(i as u32);
                i += 1;
                if i >= 624 {
                    mt[0] = mt[623];
                    i = 1;
                }
            }
            mt[0] = 0x8000_0000;
            Self { mt, idx: 624 }
        }

        fn next_u32(&mut self) -> u32 {
            if self.idx >= 624 {
                for kk in 0..624 {
                    let y = (self.mt[kk] & 0x8000_0000) | (self.mt[(kk + 1) % 624] & 0x7fff_ffff);
                    let mag = if (y & 1) != 0 { 0x9908_b0df } else { 0 };
                    self.mt[kk] = self.mt[(kk + 397) % 624] ^ (y >> 1) ^ mag;
                }
                self.idx = 0;
            }
            let mut y = self.mt[self.idx];
            self.idx += 1;
            y ^= y >> 11;
            y ^= (y << 7) & 0x9d2c_5680;
            y ^= (y << 15) & 0xefc6_0000;
            y ^= y >> 18;
            y
        }

        fn next_base(&mut self) -> u8 {
            loop {
                let v = (self.next_u32() >> 29) as usize;
                if v < 4 {
                    return b"ACGT"[v];
                }
            }
        }
    }

    let motif = "MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLMKYEFSSDIHPCNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHKYPKCMCVMVRPTDQSLRHWKIMHDANFEKLKIFHYMVGHGRTASEG";
    let codon = |aa: u8| -> &'static [u8; 3] {
        match aa {
            b'A' => b"GCT",
            b'R' => b"CGT",
            b'N' => b"AAT",
            b'D' => b"GAT",
            b'C' => b"TGT",
            b'Q' => b"CAA",
            b'E' => b"GAA",
            b'G' => b"GGT",
            b'H' => b"CAT",
            b'I' => b"ATT",
            b'L' => b"CTG",
            b'K' => b"AAA",
            b'M' => b"ATG",
            b'F' => b"TTT",
            b'P' => b"CCT",
            b'S' => b"TCT",
            b'T' => b"ACT",
            b'W' => b"TGG",
            b'Y' => b"TAT",
            b'V' => b"GTT",
            _ => unreachable!("unexpected amino acid"),
        }
    };

    let dir = TempDir::new("tbn_hard_long5m");
    let mut rng = Mt19937::seed(1640);
    let mut seq: Vec<u8> = (0..15_600_003).map(|_| rng.next_base()).collect();
    let mut nt = Vec::with_capacity(motif.len() * 3);
    for &aa in motif.as_bytes() {
        nt.extend_from_slice(codon(aa));
    }
    let pos = 4_999_880 * 3;
    seq[pos..pos + nt.len()].copy_from_slice(&nt);

    let mut lc = Vec::with_capacity(800);
    lc.extend((0..300).map(|_| rng.next_base()));
    lc.extend(std::iter::repeat_n(b'A', 200));
    lc.extend((0..300).map(|_| rng.next_base()));

    let long_fna = dir.path().join("long.fna");
    let lc_fna = dir.path().join("lc.fna");
    let q_faa = dir.path().join("q.faa");
    let lc_asnb = dir.path().join("lc.asnb");
    let db_dir = dir.path().join("db");
    std::fs::create_dir_all(&db_dir).expect("mkdir db");
    let db = db_dir.join("long");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&long_fna).expect("create long.fna");
        f.write_all(b">long1\n").expect("write");
        f.write_all(&seq).expect("write");
        f.write_all(b"\n>lc1\n").expect("write");
        f.write_all(&lc).expect("write");
        f.write_all(b"\n").expect("write");
        let mut f = std::fs::File::create(&lc_fna).expect("create lc.fna");
        f.write_all(b">lc1\n").expect("write");
        f.write_all(&lc).expect("write");
        f.write_all(b"\n").expect("write");
        std::fs::write(&q_faa, format!(">q1\n{motif}\n")).expect("write q.faa");
    }

    let dust = Command::new(BIN)
        .current_dir(engine())
        .args([
            "dustmasker",
            "-in",
            lc_fna.to_str().expect("utf8"),
            "-infmt",
            "fasta",
            "-parse_seqids",
            "-outfmt",
            "maskinfo_asn1_bin",
            "-out",
            lc_asnb.to_str().expect("utf8"),
        ])
        .stdin(Stdio::null())
        .output()
        .expect("dustmasker");
    assert!(dust.status.success(), "dustmasker: {dust:?}");

    let mk = Command::new(BIN)
        .current_dir(engine())
        .args([
            "makeblastdb",
            "-in",
            long_fna.to_str().expect("utf8"),
            "-dbtype",
            "nucl",
            "-parse_seqids",
            "-mask_data",
            lc_asnb.to_str().expect("utf8"),
            "-out",
            db.to_str().expect("utf8"),
        ])
        .stdin(Stdio::null())
        .output()
        .expect("makeblastdb");
    assert!(mk.status.success(), "makeblastdb: {mk:?}");

    // Frozen NCBI BLAST+ 2.17.0 bytes (`tblastn -query q.faa -db db/long -db_hard_mask 11 -outfmt 6`).
    let want_out = concat!(
        "q1\tlong1\t100.000\t140\t0\t0\t1\t140\t14999641\t15000060\t3.78e-93\t293\n",
        "q1\tlong1\t24.194\t62\t43\t1\t12\t69\t11606142\t11606327\t1.3\t27.7\n",
        "q1\tlong1\t42.308\t26\t15\t0\t114\t139\t12996911\t12996834\t6.6\t25.8\n",
        "q1\tlc1\t33.333\t27\t18\t0\t22\t48\t51\t131\t4.6\t25.8\n",
    );
    for threads in ["1", "4"] {
        let out = Command::new(BIN)
            .current_dir(engine())
            .args([
                "tblastn",
                "-query",
                q_faa.to_str().expect("utf8"),
                "-db",
                db.to_str().expect("utf8"),
                "-db_hard_mask",
                "11",
                "-outfmt",
                "6",
                "-num_threads",
                threads,
            ])
            .stdin(Stdio::null())
            .output()
            .expect("tblastn");
        assert_eq!(out.status.code(), Some(0), "tblastn t={threads}: {out:?}");
        assert!(
            out.stderr.is_empty(),
            "tblastn t={threads} stderr: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        assert_eq!(
            String::from_utf8_lossy(&out.stdout),
            want_out,
            "tblastn t={threads} stdout differs from NCBI 2.17.0"
        );
    }
}
