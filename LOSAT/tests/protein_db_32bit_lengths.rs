//! Protein BLAST databases containing sequences of `2^31` or more residues or `0` residues
//! (KC-11f): NCBI BLAST+ 2.17.0 truncates protein lengths to a signed 32-bit `int`
//! (`CSeqDBVol::GetSeqLengthProt`, `x_GetSequence`), while nucleotide lengths stay widened.
//!
//! Each test builds sparse `.psq` files in a private temporary directory removed on drop so no
//! gigabytes of disk or RAM are allocated.

use std::fs::{self, File, OpenOptions};
use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_losat");
const TWO_POW_31: u64 = 1u64 << 31;
const AA_PATTERN: &str = "ACDEFGHIKLMNPQRSTVWY";
/// NCBIstdaa encoding of `AA_PATTERN` (`ACDEFGHIKLMNPQRSTVWY`).
const AA_PATTERN_STDAA: [u8; 20] = [
    1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
];
const COMBINED_DBS: [&str; 4] = ["dba", "dba0", "db_mix1", "db_mix2"];
const SEARCH_PROGS: [(&str, &str); 2] = [("blastp", "q.fa"), ("blastx", "qn.fa")];
const ERR_NO_SEQ: &str = "Error: [blastdbcmd] Error: could not get sequence or range.\n";
const ERR_BEYOND: &str = "Error: [blastdbcmd] Error: region beyond sequence range.\n";
const ERR_MKDB_NO_DATA: &str = "LOSAT database error: File error: could not get sequence data.\n";
const ERR_M_CREATE: &str = "Error: basic_string::_M_create\n";

struct TempDir(PathBuf);

impl TempDir {
    fn new(tag: &str) -> Self {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos());
        let dir =
            std::env::temp_dir().join(format!("losat-prot32-{tag}-{}-{nanos}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn losat(dir: &Path, args: &[&str]) -> Output {
    let mut command = Command::new(BIN);
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        // Keep `LOSAT_BLASTX_WORKER_LOG` so `--all-features` (`blastx-worker-probe`) runs stay valid.
        if (text.starts_with("LOSAT_") && text != "LOSAT_BLASTX_WORKER_LOG")
            || text.starts_with("NCBI_")
            || text.starts_with("BLAST")
        {
            command.env_remove(&name);
        }
    }
    let temp_worker_log = std::env::var_os("LOSAT_BLASTX_WORKER_LOG")
        .is_none()
        .then(|| {
            let path = dir.join(".blastx_worker.log");
            command.env("LOSAT_BLASTX_WORKER_LOG", &path);
            path
        });
    let out = command
        .env("NCBI_CONFIG_PATH", dir)
        .current_dir(dir)
        .args(args)
        .output()
        .unwrap();
    if let Some(path) = temp_worker_log {
        let _ = fs::remove_file(path);
    }
    out
}

fn out_str(bytes: &[u8]) -> std::borrow::Cow<'_, str> {
    String::from_utf8_lossy(bytes)
}

fn assert_ok_empty(out: &Output, label: impl std::fmt::Debug) {
    assert_eq!(
        (out.status.code(), out.stdout.len(), out.stderr.len()),
        (Some(0), 0, 0),
        "{label:?}"
    );
}

fn assert_rc_err(out: &Output, rc: i32, err: &str, label: &str) {
    assert_eq!(
        (out.status.code(), out_str(&out.stderr).as_ref()),
        (Some(rc), err),
        "{label}"
    );
}

/// Writes a sparse file of `total_len` bytes with `prefix` at offset 0 and `tail` at the end.
/// Returns `false` if Windows `fsutil sparse setflag` fails so callers skip rather than writing
/// multi-gigabyte real files.
fn make_sparse_file(path: &Path, prefix: &[u8], total_len: u64, tail: &[u8]) -> bool {
    {
        let mut f = File::create(path).unwrap();
        f.write_all(prefix).unwrap();
    }
    #[cfg(windows)]
    {
        let ok = Command::new("fsutil")
            .args(["sparse", "setflag"])
            .arg(path)
            .output()
            .is_ok_and(|out| out.status.success());
        if !ok {
            let _ = fs::remove_file(path);
            return false;
        }
    }
    let mut f = OpenOptions::new().write(true).open(path).unwrap();
    f.seek(SeekFrom::Start(total_len - tail.len() as u64))
        .unwrap();
    f.write_all(tail).unwrap();
    assert_eq!(fs::metadata(path).unwrap().len(), total_len);
    true
}

fn make_prot_db(dir: &Path, name: &str, fasta: &str) {
    let fa = dir.join(format!("{name}.fa"));
    fs::write(&fa, fasta).unwrap();
    let out = losat(
        dir,
        &[
            "makeblastdb",
            "-dbtype",
            "prot",
            "-in",
            fa.to_str().unwrap(),
            "-out",
            name,
            "-title",
            name,
        ],
    );
    assert!(out.status.success(), "{name}: {out:?}");
}

fn load_pin_and_psq(
    dir: &Path,
    name: &str,
    tail_bytes: usize,
) -> (PathBuf, PathBuf, Vec<u8>, usize, i32, i64) {
    let pin_path = dir.join(format!("{name}.pin"));
    let psq_path = dir.join(format!("{name}.psq"));
    let pin = fs::read(&pin_path).unwrap();
    let t = pin.len() - tail_bytes;
    let num = i32::from_be_bytes(pin[t - 16..t - 12].try_into().unwrap());
    let total = i64::from_le_bytes(pin[t - 12..t - 4].try_into().unwrap());
    (pin_path, psq_path, pin, t, num, total)
}

/// Patches a single-sequence 300-residue protein volume to `new_len` residues (`0` or `>= 2^31`)
/// using a sparse `.psq` file, placing `tail_residues` immediately before the trailing sentinel.
fn patch_single_prot_db(dir: &Path, name: &str, new_len: u64, tail_residues: &[u8]) -> bool {
    let (pin_path, psq_path, mut pin, t, num, total) = load_pin_and_psq(dir, name, 16);
    let s0 = u32::from_be_bytes(pin[t + 8..t + 12].try_into().unwrap());
    let s1 = u32::from_be_bytes(pin[t + 12..t + 16].try_into().unwrap());
    assert_eq!((num, total, s0, s1), (1, 300, 1, 302));

    let psq = fs::read(&psq_path).unwrap();
    if new_len == 0 {
        fs::write(&psq_path, [0u8, 0u8]).unwrap();
    } else {
        let mut tail = tail_residues.to_vec();
        tail.push(0u8);
        if !make_sparse_file(&psq_path, &psq[..301], new_len + 2, &tail) {
            return false;
        }
    }
    pin[t - 12..t - 4].copy_from_slice(&(new_len as i64).to_le_bytes());
    pin[t - 4..t].copy_from_slice(&(new_len as u32).to_be_bytes());
    pin[t + 12..t + 16].copy_from_slice(&((new_len + 2) as u32).to_be_bytes());
    fs::write(&pin_path, pin).unwrap();
    true
}

/// Patches a two-sequence protein volume (`[big(300), small(400)]` or `[small(400), big(300)]`)
/// so the 300-residue entry has `2^31` residues in a sparse `.psq`.
fn patch_mixed_prot_db(dir: &Path, name: &str, big_first: bool) -> bool {
    let (pin_path, psq_path, mut pin, t, num, total) = load_pin_and_psq(dir, name, 24);
    let s0 = u32::from_be_bytes(pin[t + 12..t + 16].try_into().unwrap());
    let s1 = u32::from_be_bytes(pin[t + 16..t + 20].try_into().unwrap());
    let s2 = u32::from_be_bytes(pin[t + 20..t + 24].try_into().unwrap());
    assert_eq!((num, total, s0), (2, 700, 1));

    let psq = fs::read(&psq_path).unwrap();
    let new_total = (total as u64 - 300 + TWO_POW_31) as i64;
    let new_s2 = if big_first {
        assert_eq!(s1 - s0 - 1, 300);
        let mut tail = vec![0u8];
        tail.extend_from_slice(&psq[s1 as usize..s2 as usize]);
        let new_s1 = s0 as u64 + TWO_POW_31 + 1;
        let end = new_s1 + (s2 - s1) as u64;
        if !make_sparse_file(&psq_path, &psq[..s1 as usize - 1], end, &tail) {
            return false;
        }
        pin[t + 16..t + 20].copy_from_slice(&(new_s1 as u32).to_be_bytes());
        end
    } else {
        assert_eq!(s2 - s1 - 1, 300);
        let end = s1 as u64 + TWO_POW_31 + 1;
        if !make_sparse_file(&psq_path, &psq[..s2 as usize - 1], end, &[0u8]) {
            return false;
        }
        end
    };
    pin[t + 20..t + 24].copy_from_slice(&(new_s2 as u32).to_be_bytes());
    pin[t - 12..t - 4].copy_from_slice(&new_total.to_le_bytes());
    pin[t - 4..t].copy_from_slice(&(TWO_POW_31 as u32).to_be_bytes());
    fs::write(&pin_path, pin).unwrap();
    true
}

fn write_oidlist(path: &Path, num_oids: u32) {
    let mut bytes = num_oids.to_be_bytes().to_vec();
    let words = (num_oids as usize).div_ceil(32);
    for w in 0..words {
        let mut word = 0u32;
        for bit in 0..32 {
            let oid = (w * 32 + bit) as u32;
            if oid < num_oids {
                word |= 1u32 << (31 - bit);
            }
        }
        bytes.extend_from_slice(&word.to_be_bytes());
    }
    fs::write(path, bytes).unwrap();
}

fn dna_for_protein(prot: &str) -> String {
    prot.bytes()
        .map(|b| match b {
            b'A' => "GCT",
            b'C' => "TGT",
            b'D' => "GAT",
            b'E' => "GAA",
            b'F' => "TTT",
            b'G' => "GGT",
            b'H' => "CAT",
            b'I' => "ATT",
            b'K' => "AAA",
            b'L' => "CTG",
            b'M' => "ATG",
            b'N' => "AAT",
            b'P' => "CCG",
            b'Q' => "CAG",
            b'R' => "CGT",
            b'S' => "TCT",
            b'T' => "ACC",
            b'V' => "GTT",
            b'W' => "TGG",
            b'Y' => "TAT",
            _ => unreachable!(),
        })
        .collect()
}

struct Fixture {
    dir: TempDir,
}

impl Fixture {
    fn build(tag: &str) -> Option<Self> {
        let dir = TempDir::new(tag);
        let big = AA_PATTERN.repeat(15);
        let small = format!(
            "{}{}{}",
            AA_PATTERN.repeat(15),
            "Q".repeat(40),
            AA_PATTERN.repeat(3)
        );
        let query = AA_PATTERN.repeat(3);

        let big_fa = format!(">big\n{big}\n");
        let small_fa = format!(">small\n{small}\n");
        fs::write(dir.0.join("q.fa"), format!(">q\n{query}\n")).unwrap();
        fs::write(
            dir.0.join("qn.fa"),
            format!(">qn\n{}\n", dna_for_protein(&query)),
        )
        .unwrap();

        make_prot_db(&dir.0, "dbh", &big_fa);
        if !patch_single_prot_db(&dir.0, "dbh", TWO_POW_31, &AA_PATTERN_STDAA) {
            return None;
        }

        make_prot_db(&dir.0, "dbg", &big_fa);
        if !patch_single_prot_db(&dir.0, "dbg", TWO_POW_31 + 10, &AA_PATTERN_STDAA) {
            return None;
        }

        make_prot_db(&dir.0, "db0", &big_fa);
        assert!(patch_single_prot_db(&dir.0, "db0", 0, &[]));

        make_prot_db(&dir.0, "dbs", &small_fa);

        make_prot_db(&dir.0, "dbz", &big_fa);
        if !patch_single_prot_db(&dir.0, "dbz", TWO_POW_31, &[]) {
            return None;
        }

        make_prot_db(&dir.0, "dbz0", &big_fa);
        assert!(patch_single_prot_db(&dir.0, "dbz0", 0, &[]));

        make_prot_db(
            &dir.0,
            "db_mix1",
            &format!(">big\n{big}\n>small\n{small}\n"),
        );
        if !patch_mixed_prot_db(&dir.0, "db_mix1", true) {
            return None;
        }

        make_prot_db(
            &dir.0,
            "db_mix2",
            &format!(">small\n{small}\n>big\n{big}\n"),
        );
        if !patch_mixed_prot_db(&dir.0, "db_mix2", false) {
            return None;
        }

        fs::write(dir.0.join("dba.pal"), "TITLE dba\nDBLIST dbh dbs\n").unwrap();
        fs::write(dir.0.join("dba0.pal"), "TITLE dba0\nDBLIST db0 dbs\n").unwrap();
        fs::write(dir.0.join("dba_sz.pal"), "TITLE dba_sz\nDBLIST dbs dbz\n").unwrap();
        fs::write(
            dir.0.join("dba_sz0.pal"),
            "TITLE dba_sz0\nDBLIST dbs dbz0\n",
        )
        .unwrap();
        fs::write(
            dir.0.join("dbh_len0.pal"),
            "TITLE dbh_len0\nDBLIST dbh\nLENGTH 0\n",
        )
        .unwrap();

        write_oidlist(&dir.0.join("ol1.msk"), 1);
        write_oidlist(&dir.0.join("ol2.msk"), 2);
        fs::write(
            dir.0.join("dbh_ol.pal"),
            "TITLE dbh_ol\nDBLIST dbh\nOIDLIST ol1.msk\n",
        )
        .unwrap();
        fs::write(
            dir.0.join("dba_ol.pal"),
            "TITLE dba_ol\nDBLIST dbh dbs\nOIDLIST ol2.msk\n",
        )
        .unwrap();

        Some(Self { dir })
    }

    fn run(&self, args: &[&str]) -> Output {
        losat(&self.dir.0, args)
    }
}

#[test]
fn blastdbcmd_and_segmasker_match_ncbi_32bit_protein_lengths() {
    let Some(fx) = Fixture::build("cmd-seg") else {
        eprintln!(
            "skipped: blastdbcmd_and_segmasker_match_ncbi_32bit_protein_lengths (sparse file creation unavailable)"
        );
        return;
    };

    // -info and -info -exact_length
    for (db, total, longest) in [
        ("dbh", "2,147,483,648", "0"),
        ("dbg", "2,147,483,658", "0"),
        ("dba", "2,147,484,048", "400"),
        ("db0", "0", "0"),
        ("dba0", "400", "400"),
        ("dbh_ol", "1,846,744,071,562,067,968", "0"),
        ("dba_ol", "1,846,744,071,562,068,368", "400"),
        ("db_mix1", "2,147,484,048", "0"),
        ("db_mix2", "2,147,484,048", "0"),
    ] {
        for extra in [Vec::new(), vec!["-exact_length"]] {
            let mut args = vec!["blastdbcmd", "-db", db, "-dbtype", "prot", "-info"];
            args.extend(extra);
            let out = fx.run(&args);
            assert_eq!(out.status.code(), Some(0), "{args:?}: {out:?}");
            let text = out_str(&out.stdout);
            assert!(
                text.contains(&format!("{total} total residues")),
                "{db}: {text}"
            );
            assert!(
                text.contains(&format!("Longest sequence: {longest} residues")),
                "{db}: {text}"
            );
        }
    }

    // %o %l
    for (db, expected) in [
        ("dbh", "0 -2147483648\n"),
        ("dbg", "0 -2147483638\n"),
        ("db0", "0 0\n"),
        ("dba", "0 -2147483648\n1 400\n"),
        ("dba0", "0 0\n1 400\n"),
        ("db_mix1", "0 -2147483648\n1 400\n"),
        ("db_mix2", "0 400\n1 -2147483648\n"),
    ] {
        let out = fx.run(&[
            "blastdbcmd",
            "-db",
            db,
            "-dbtype",
            "prot",
            "-entry",
            "all",
            "-outfmt",
            "%o %l",
        ]);
        assert_eq!(
            (
                out.status.code(),
                out_str(&out.stdout).as_ref(),
                out_str(&out.stderr).as_ref()
            ),
            (Some(0), expected, ""),
            "{db}"
        );
    }

    // Default FASTA and %b silently skip <= 0 proteins
    for db in ["dbh", "db0"] {
        for extra in [Vec::new(), vec!["-outfmt", "%b"]] {
            let mut args = vec!["blastdbcmd", "-db", db, "-dbtype", "prot", "-entry", "all"];
            args.extend(extra);
            let out = fx.run(&args);
            assert_ok_empty(&out, &args);
        }
    }
    for db in COMBINED_DBS {
        let out = fx.run(&["blastdbcmd", "-db", db, "-dbtype", "prot", "-entry", "all"]);
        assert_eq!(out.status.code(), Some(0), "fasta {db}");
        let fasta = out_str(&out.stdout);
        assert!(fasta.starts_with(">small\n"), "{db}: {fasta}");
        assert!(!fasta.contains(">big"), "{db}: {fasta}");
    }

    // %s, %h, and -range errors match CSeqDBVol::x_GetAmbigSeq
    for (db, fmt, range, expected_err) in [
        ("dbh", "%s", None, ERR_NO_SEQ),
        ("dbh", "%o %h", None, ERR_NO_SEQ),
        ("dbg", "%s", None, ERR_NO_SEQ),
        ("db0", "%s", None, ERR_NO_SEQ),
        ("db0", "%o %h", None, ERR_NO_SEQ),
        ("dba", "%s", None, ERR_NO_SEQ),
        ("dbh", "%s", Some("0-"), ERR_NO_SEQ),
        ("dbh", "%s", Some("1-10"), ERR_BEYOND),
        ("dbh", "%s", Some("1-1"), ERR_BEYOND),
        ("dbh", "%s", Some("1-"), ERR_NO_SEQ),
        ("dbg", "%s", Some("0-"), ERR_NO_SEQ),
        ("dbg", "%s", Some("1-10"), ERR_BEYOND),
        ("dbg", "%s", Some("2147483646-2147483647"), ERR_BEYOND),
        ("dbg", "%s", Some("1-"), ERR_NO_SEQ),
        ("dbg", "%s", Some("5-"), ERR_NO_SEQ),
        ("dbg", "%s", Some("11-"), ERR_NO_SEQ),
        (
            "dbg",
            "%s",
            Some("2147483649-2147483650"),
            "Error: [blastdbcmd] Cannot convert string '2147483649' to int, overflow\n",
        ),
        ("db0", "%s", Some("0-"), ERR_NO_SEQ),
        ("db0", "%s", Some("1-10"), ERR_NO_SEQ),
        ("db0", "%s", Some("1-"), ERR_NO_SEQ),
        ("db0", "%s", Some("10-20"), ERR_NO_SEQ),
    ] {
        let mut args = vec![
            "blastdbcmd",
            "-db",
            db,
            "-dbtype",
            "prot",
            "-entry",
            "all",
            "-outfmt",
            fmt,
        ];
        if let Some(r) = range {
            args.extend(["-range", r]);
        }
        let out = fx.run(&args);
        assert_eq!(
            (
                out.status.code(),
                out.stdout.len(),
                out_str(&out.stderr).as_ref()
            ),
            (Some(1), 0, expected_err),
            "{args:?}"
        );
    }

    // Open-ended `-range N-` on `2^31 + k` proteins with `N >= k + 2` wraps `end - begin`
    // positive in 32-bit `int` and prints the exact tail residues matching NCBI 2.17.0.
    for (db, range, expected_seq) in [
        ("dbg", "2147483647-", "KLMNPQRSTVWY\n"),
        ("dbg", "2147483640-", "CDEFGHIKLMNPQRSTVWY\n"),
        ("dbh", "2147483647-", "WY\n"),
        ("dbh", "2147483646-", "VWY\n"),
    ] {
        let out = fx.run(&[
            "blastdbcmd",
            "-db",
            db,
            "-dbtype",
            "prot",
            "-entry",
            "all",
            "-outfmt",
            "%s",
            "-range",
            range,
        ]);
        assert_eq!(
            (
                out.status.code(),
                out_str(&out.stdout).as_ref(),
                out_str(&out.stderr).as_ref()
            ),
            (Some(0), expected_seq, ""),
            "{db} -range {range}"
        );
    }

    // db_mix2 prints OID 0's sequence before failing on OID 1 (2^31)
    let mix2_seq = fx.run(&[
        "blastdbcmd",
        "-db",
        "db_mix2",
        "-dbtype",
        "prot",
        "-entry",
        "all",
        "-outfmt",
        "%s",
    ]);
    assert_eq!(mix2_seq.status.code(), Some(1));
    assert_eq!(mix2_seq.stdout.len(), 401);
    assert_eq!(out_str(&mix2_seq.stderr), ERR_NO_SEQ);

    // segmasker skips <= 0 proteins cleanly across formats
    for db in ["dbh", "db0"] {
        for outfmt in ["interval", "fasta"] {
            let out = fx.run(&[
                "segmasker",
                "-infmt",
                "blastdb",
                "-in",
                db,
                "-outfmt",
                outfmt,
            ]);
            assert_ok_empty(&out, (db, outfmt));
        }
    }
    for db in COMBINED_DBS {
        let out = fx.run(&["segmasker", "-infmt", "blastdb", "-in", db]);
        assert_eq!(
            (
                out.status.code(),
                out_str(&out.stdout).as_ref(),
                out_str(&out.stderr).as_ref()
            ),
            (Some(0), ">small\n300 - 339\n", ""),
            "seg {db}"
        );
    }
}

#[test]
fn blastp_blastx_psiblast_skip_non_positive_32bit_proteins() {
    let Some(fx) = Fixture::build("search") else {
        eprintln!(
            "skipped: blastp_blastx_psiblast_skip_non_positive_32bit_proteins (sparse file creation unavailable)"
        );
        return;
    };

    // Empty / 2^31-only DBs produce no hits with rc 0 for gapped searches
    for db in ["dbh", "db0", "dbh_ol"] {
        for args in [
            vec![
                "blastp",
                "-query",
                "q.fa",
                "-db",
                db,
                "-outfmt",
                "6 std slen",
            ],
            vec![
                "psiblast",
                "-query",
                "q.fa",
                "-db",
                db,
                "-num_iterations",
                "2",
                "-outfmt",
                "6 std slen",
            ],
        ] {
            let out = fx.run(&args);
            assert_ok_empty(&out, &args);
        }
    }
    for db in ["dbh", "db0"] {
        let out = fx.run(&[
            "blastx",
            "-query",
            "qn.fa",
            "-db",
            db,
            "-outfmt",
            "6 std slen",
        ]);
        assert_ok_empty(&out, ("blastx", db));
    }

    // Ungapped searches on dbh / dbh_ol succeed with 0 hits; db0 has 0 total length and fails setup
    for db in ["dbh", "dbh_ol"] {
        for (prog, q) in SEARCH_PROGS {
            let out = fx.run(&[
                prog,
                "-query",
                q,
                "-db",
                db,
                "-ungapped",
                "-comp_based_stats",
                "0",
                "-outfmt",
                "6 std slen",
            ]);
            assert_ok_empty(&out, (prog, db));
        }
    }
    for (prog, q) in SEARCH_PROGS {
        let out = fx.run(&[
            prog,
            "-query",
            q,
            "-db",
            "db0",
            "-ungapped",
            "-comp_based_stats",
            "0",
            "-outfmt",
            "6 std slen",
        ]);
        assert_eq!(
            (
                out.status.code(),
                out.stdout.len(),
                out_str(&out.stderr).as_ref()
            ),
            (
                Some(3),
                0,
                "LOSAT engine error: The average subject length is too short\n"
            ),
            "{prog} ug db0"
        );
    }

    // OID-0 statistics fallback when total_length <= 0 (dbh_len0 and dbh_ol)
    for db in ["dbh_len0", "dbh_ol"] {
        let out = fx.run(&["blastp", "-query", "q.fa", "-db", db]);
        assert_eq!(out.status.code(), Some(0), "len0 {db}");
        let report = out_str(&out.stdout);
        assert!(
            report.contains("Effective search space used: 60\n"),
            "{db}: {report}"
        );
    }

    // blastx gapped on dbh_ol / dba_ol matches NCBI's Unknown error code 1 (negative min_length)
    for db in ["dbh_ol", "dba_ol"] {
        let out = fx.run(&[
            "blastx",
            "-query",
            "qn.fa",
            "-db",
            db,
            "-outfmt",
            "6 std slen",
        ]);
        assert_rc_err(&out, 3, "LOSAT engine error: Unknown error code 1\n", db);
    }

    // Combined & single-volume mixed DBs: verified against official NCBI BLAST+ 2.17.0.
    // At 1 thread all 4 skip the non-positive OID and find hits on `small`. At 2 threads
    // (DB-split), NCBI's `CSeqDBImpl::x_FillSeqBuffer` (`seqdbimpl.cpp:699-700`) stops at the
    // first negative-length OID, so `dba` and `db_mix1` (OID 0 = -2147483648) return rc=0 with
    // empty stdout/stderr for both `blastp` and `blastx` (gapped and ungapped), while `db_mix2`
    // (OID 0 = 400) matches its 1-thread output byte for byte.
    // Open difference (tracked as a follow-up item, not parity): on `dba0` (whose volume 0 has
    // length 0), NCBI 2.17.0 aborts (rc 134, `CSeqDBException` at `seqdbimpl.cpp:660`) at 2+
    // threads, while LOSAT currently skips OID 0 and prints the 1-thread output.
    for db in COMBINED_DBS {
        let mt_stops_before_small = matches!(db, "dba" | "db_mix1");
        for base in [
            vec![
                "blastp",
                "-query",
                "q.fa",
                "-db",
                db,
                "-outfmt",
                "6 std slen",
            ],
            vec![
                "blastp",
                "-query",
                "q.fa",
                "-db",
                db,
                "-ungapped",
                "-comp_based_stats",
                "0",
                "-outfmt",
                "6 std slen",
            ],
            vec![
                "blastx",
                "-query",
                "qn.fa",
                "-db",
                db,
                "-outfmt",
                "6 std slen",
            ],
            vec![
                "blastx",
                "-query",
                "qn.fa",
                "-db",
                db,
                "-ungapped",
                "-comp_based_stats",
                "0",
                "-outfmt",
                "6 std slen",
            ],
        ] {
            let mut a1 = base.clone();
            a1.extend(["-num_threads", "1"]);
            let mut a2 = base.clone();
            a2.extend(["-num_threads", "2"]);
            let r1 = fx.run(&a1);
            let r2 = fx.run(&a2);
            assert_eq!(r1.status.code(), Some(0), "{a1:?}: {r1:?}");
            if mt_stops_before_small {
                assert_ok_empty(&r2, &a2);
            } else {
                assert_eq!(
                    (r1.status.code(), &r1.stdout, &r1.stderr),
                    (r2.status.code(), &r2.stdout, &r2.stderr),
                    "mt {base:?}"
                );
            }
            let hits = out_str(&r1.stdout);
            assert!(hits.contains("\tsmall\t100.000\t60\t"), "{base:?}: {hits}");
        }
        let psi = fx.run(&[
            "psiblast",
            "-query",
            "q.fa",
            "-db",
            db,
            "-num_iterations",
            "2",
            "-outfmt",
            "6 std slen",
        ]);
        assert_eq!(psi.status.code(), Some(0), "psi {db}");
        assert!(
            out_str(&psi.stdout).contains("\tsmall\t100.000\t60\t"),
            "psi {db}"
        );
    }
}

#[test]
fn makeblastdb_raw_copy_matches_ncbi_errors_and_partial_close() {
    let Some(fx) = Fixture::build("remk") else {
        eprintln!(
            "skipped: makeblastdb_raw_copy_matches_ncbi_errors_and_partial_close (sparse file creation unavailable)"
        );
        return;
    };

    // Failing at OID 0 creates no output DB files
    for (db, rc, err) in [
        ("dbh", 255, ERR_M_CREATE),
        ("dba", 255, ERR_M_CREATE),
        ("db_mix1", 255, ERR_M_CREATE),
        ("db0", 2, ERR_MKDB_NO_DATA),
        ("dba0", 2, ERR_MKDB_NO_DATA),
    ] {
        let outname = format!("out_{db}");
        let out = fx.run(&[
            "makeblastdb",
            "-in",
            db,
            "-input_type",
            "blastdb",
            "-dbtype",
            "prot",
            "-out",
            &outname,
        ]);
        assert_rc_err(&out, rc, err, db);
        assert!(
            !fx.dir.0.join(format!("{outname}.pin")).exists(),
            "{db} should not leave .pin"
        );
    }

    // Failing after a valid OID 0 publishes & closes the partial DB containing `small`
    for (db, rc, err) in [
        ("dba_sz", 255, ERR_M_CREATE),
        ("db_mix2", 255, ERR_M_CREATE),
        ("dba_sz0", 2, ERR_MKDB_NO_DATA),
    ] {
        let outname = format!("out_{db}");
        let out = fx.run(&[
            "makeblastdb",
            "-in",
            db,
            "-input_type",
            "blastdb",
            "-dbtype",
            "prot",
            "-out",
            &outname,
        ]);
        assert_rc_err(&out, rc, err, db);
        assert!(
            fx.dir.0.join(format!("{outname}.pin")).exists(),
            "{db} should close partial .pin"
        );
        let info = fx.run(&["blastdbcmd", "-db", &outname, "-dbtype", "prot", "-info"]);
        assert_eq!(info.status.code(), Some(0), "partial info {db}");
        let summary = out_str(&info.stdout);
        assert!(
            summary.contains("1 sequences; 400 total residues"),
            "partial info {db}: {summary}"
        );
    }
}
