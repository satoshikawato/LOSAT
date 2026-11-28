//! CLI integration tests for `-logfile` self-reference refusal (line WQA D23, Owner decision
//! `Q-D23=A` / WQA DECISIONS 49).
//!
//! Verifies that across search programs and tools, specifying `-logfile` pointing to an input
//! file read by the run (query, subject, `.ncbirc`, `-conffile`, `<program>.ini`, `-in`,
//! `-ustat`, `-in_pssm`, `-in_msa`, `-seqidlist`, `-gilist`, `-taxidlist`,
//! `-import_search_strategy`, `-mask_data`, `-entry_batch`, `-db` component, or via relative path
//! / symlink / hard link / leading space / trailing slash) is refused with exit code 1 before
//! creating or modifying the file, leaving existing inputs' bytes and mtime unchanged, while
//! non-self-referencing `-logfile` runs work as before.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::SystemTime;

const BIN: &str = env!("CARGO_BIN_EXE_losat");

const NUCL_FASTA: &[u8] = b">seq1\nATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG\n";
const PROT_FASTA: &[u8] = b">prot1\nMKTLLLTLVVVTIVCLDLGYTLKGSQATQVVGGVTGVTGVTGVTGVTGVTGVTGVTGVTG\n";

fn fresh_dir(name: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("losat-d23-self-log-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn run_in(dir: &Path, args: &[&str], extra_env: &[(&str, &str)]) -> Output {
    let mut cmd = Command::new(BIN);
    cmd.current_dir(dir)
        .args(args)
        .env_remove("NCBI_CONFIG_PATH")
        .env_remove("NCBI_DONT_USE_NCBIRC")
        .env_remove("NCBI_DONT_USE_LOCAL_CONFIG")
        .env_remove("BLAST_USAGE_REPORT")
        .env_remove("BLASTDB")
        .env("HOME", dir);
    for (k, v) in extra_env {
        cmd.env(k, v);
    }
    cmd.output().expect("failed to execute losat")
}

fn snapshot(path: &Path) -> (Vec<u8>, SystemTime) {
    let bytes = fs::read(path).unwrap();
    let mtime = fs::metadata(path).unwrap().modified().unwrap();
    (bytes, mtime)
}

fn assert_refused_and_unchanged(
    out: &Output,
    watched: &[(&Path, &(Vec<u8>, SystemTime))],
    label: &str,
) {
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(
        out.status.code(),
        Some(1),
        "{label}: expected rc 1, got {:?}, stderr={stderr}",
        out.status.code()
    );
    assert!(
        stderr.contains("Error: a -logfile file") && stderr.contains("is not supported by LOSAT's"),
        "{label}: unexpected stderr: {stderr}"
    );
    for (path, (old_bytes, old_mtime)) in watched {
        let (new_bytes, new_mtime) = snapshot(path);
        assert_eq!(
            &new_bytes,
            old_bytes,
            "{label}: bytes of {} were modified",
            path.display()
        );
        assert_eq!(
            &new_mtime,
            old_mtime,
            "{label}: mtime of {} changed",
            path.display()
        );
    }
}

#[test]
fn search_programs_refuse_self_referencing_logfile() {
    let dir = fresh_dir("searches");
    fs::create_dir_all(dir.join("sub")).unwrap();
    let nucl_q = dir.join("nq.fa");
    let nucl_s = dir.join("ns.fa");
    let prot_q = dir.join("pq.fa");
    let prot_s = dir.join("ps.fa");
    let ncbirc = dir.join(if cfg!(windows) { "ncbi.ini" } else { ".ncbirc" });
    let conf = dir.join("custom.ini");
    let list = dir.join("ids.list");
    let pssm = dir.join("check.pssm");
    let msa = dir.join("align.msa");
    let strat = dir.join("strat.asn");

    fs::write(&nucl_q, NUCL_FASTA).unwrap();
    fs::write(&nucl_s, NUCL_FASTA).unwrap();
    fs::write(&prot_q, PROT_FASTA).unwrap();
    fs::write(&prot_s, PROT_FASTA).unwrap();
    fs::write(&ncbirc, b"\xef").unwrap();
    fs::write(&conf, b"[BLAST]\nBLASTDB=.\n").unwrap();
    fs::write(&list, b"12345\n").unwrap();
    fs::write(&pssm, b"placeholder pssm\n").unwrap();
    fs::write(&msa, PROT_FASTA).unwrap();
    fs::write(&strat, b"Blast4-request ::= {}\n").unwrap();

    let nq_snap = snapshot(&nucl_q);
    let ns_snap = snapshot(&nucl_s);
    let pq_snap = snapshot(&prot_q);
    let ps_snap = snapshot(&prot_s);
    let rc_snap = snapshot(&ncbirc);
    let conf_snap = snapshot(&conf);
    let list_snap = snapshot(&list);
    let pssm_snap = snapshot(&pssm);
    let msa_snap = snapshot(&msa);
    let strat_snap = snapshot(&strat);

    // Every search program with -query == -logfile (and -subject == -logfile where supported).
    let query_subject_cases: &[(
        &str,
        &str,
        &str,
        &Path,
        &(Vec<u8>, SystemTime),
        &Path,
        &(Vec<u8>, SystemTime),
    )] = &[
        (
            "blastn", "nq.fa", "ns.fa", &nucl_q, &nq_snap, &nucl_s, &ns_snap,
        ),
        (
            "blastp", "pq.fa", "ps.fa", &prot_q, &pq_snap, &prot_s, &ps_snap,
        ),
        (
            "blastx", "nq.fa", "ps.fa", &nucl_q, &nq_snap, &prot_s, &ps_snap,
        ),
        (
            "tblastn", "pq.fa", "ns.fa", &prot_q, &pq_snap, &nucl_s, &ns_snap,
        ),
        (
            "tblastx", "nq.fa", "ns.fa", &nucl_q, &nq_snap, &nucl_s, &ns_snap,
        ),
        (
            "psiblast", "pq.fa", "ps.fa", &prot_q, &pq_snap, &prot_s, &ps_snap,
        ),
    ];

    for &(prog, q, s, q_path, q_snap, s_path, s_snap) in query_subject_cases {
        let out_q = run_in(
            &dir,
            &[prog, "-query", q, "-subject", s, "-logfile", q],
            &[],
        );
        assert_refused_and_unchanged(
            &out_q,
            &[(q_path, q_snap), (s_path, s_snap)],
            &format!("{prog} -logfile query"),
        );

        let out_s = run_in(
            &dir,
            &[prog, "-query", q, "-subject", s, "-logfile", s],
            &[],
        );
        assert_refused_and_unchanged(
            &out_s,
            &[(q_path, q_snap), (s_path, s_snap)],
            &format!("{prog} -logfile subject"),
        );

        let rc_name = if cfg!(windows) { "ncbi.ini" } else { ".ncbirc" };
        let out_rc = run_in(
            &dir,
            &[prog, "-query", q, "-subject", s, "-logfile", rc_name],
            &[],
        );
        assert_refused_and_unchanged(
            &out_rc,
            &[(&ncbirc, &rc_snap)],
            &format!("{prog} -logfile ncbirc"),
        );
    }

    // rpsblast and rpstblastn (-query == -logfile, -import_search_strategy == -logfile).
    for &(prog, q, q_path, q_snap) in &[
        ("rpsblast", "pq.fa", &prot_q, &pq_snap),
        ("rpstblastn", "nq.fa", &nucl_q, &nq_snap),
    ] {
        let out = run_in(
            &dir,
            &[prog, "-query", q, "-db", "no_db", "-logfile", q],
            &[],
        );
        assert_refused_and_unchanged(&out, &[(q_path, q_snap)], &format!("{prog} -logfile query"));

        let out_st = run_in(
            &dir,
            &[
                prog,
                "-import_search_strategy",
                "strat.asn",
                "-logfile",
                "strat.asn",
            ],
            &[],
        );
        assert_refused_and_unchanged(
            &out_st,
            &[(&strat, &strat_snap)],
            &format!("{prog} -logfile strat"),
        );
    }

    // Extra search inputs: -conffile, -gilist, -seqidlist, -taxidlist, -in_pssm, -in_msa,
    // -import_search_strategy, plus value-skipping (`-export_search_strategy -query -query nq.fa`).
    let out_conf = run_in(
        &dir,
        &[
            "blastn",
            "-query",
            "nq.fa",
            "-subject",
            "ns.fa",
            "-conffile",
            "custom.ini",
            "-logfile",
            "custom.ini",
        ],
        &[],
    );
    assert_refused_and_unchanged(&out_conf, &[(&conf, &conf_snap)], "blastn -conffile self");

    for flag in ["-gilist", "-seqidlist", "-taxidlist"] {
        let out = run_in(
            &dir,
            &[
                "blastn", "-query", "nq.fa", "-db", "no_db", flag, "ids.list", "-logfile",
                "ids.list",
            ],
            &[],
        );
        assert_refused_and_unchanged(&out, &[(&list, &list_snap)], &format!("blastn {flag} self"));
    }

    let out_pssm = run_in(
        &dir,
        &[
            "psiblast",
            "-in_pssm",
            "check.pssm",
            "-subject",
            "ps.fa",
            "-logfile",
            "check.pssm",
        ],
        &[],
    );
    assert_refused_and_unchanged(&out_pssm, &[(&pssm, &pssm_snap)], "psiblast -in_pssm self");

    let out_msa = run_in(
        &dir,
        &[
            "psiblast",
            "-in_msa",
            "align.msa",
            "-subject",
            "ps.fa",
            "-logfile",
            "align.msa",
        ],
        &[],
    );
    assert_refused_and_unchanged(&out_msa, &[(&msa, &msa_snap)], "psiblast -in_msa self");

    let out_strat = run_in(
        &dir,
        &[
            "blastn",
            "-import_search_strategy",
            "strat.asn",
            "-logfile",
            "strat.asn",
        ],
        &[],
    );
    assert_refused_and_unchanged(
        &out_strat,
        &[(&strat, &strat_snap)],
        "blastn -import_search_strategy self",
    );

    // Option value equal to "-query" must be skipped so real "-query nq.fa" is still detected.
    let out_skip = run_in(
        &dir,
        &[
            "tblastx",
            "-export_search_strategy",
            "-query",
            "-query",
            "nq.fa",
            "-subject",
            "ns.fa",
            "-logfile",
            "nq.fa",
        ],
        &[],
    );
    assert_refused_and_unchanged(
        &out_skip,
        &[(&nucl_q, &nq_snap)],
        "tblastx value=-query self",
    );

    // Path spellings: relative ./, sub/../, leading space, trailing slash, =form, absolute path.
    let abs_q = nucl_q.to_string_lossy().into_owned();
    let eq_q = "-logfile=nq.fa";
    for (label, log_args) in [
        ("dot-slash", vec!["-logfile", "./nq.fa"]),
        ("dot-dot", vec!["-logfile", "sub/../nq.fa"]),
        ("leading-space", vec!["-logfile", " nq.fa"]),
        ("trailing-slash", vec!["-logfile", "nq.fa/"]),
        ("equals-form", vec![eq_q]),
        ("absolute", vec!["-logfile", abs_q.as_str()]),
    ] {
        let mut argv = vec!["tblastx", "-query", "nq.fa", "-subject", "ns.fa"];
        argv.extend(log_args);
        let out = run_in(&dir, &argv, &[]);
        assert_refused_and_unchanged(&out, &[(&nucl_q, &nq_snap)], label);
    }

    // Missing input file equal to -logfile must be refused and NOT created.
    let missing = dir.join("not_here.fa");
    let out_missing = run_in(
        &dir,
        &[
            "blastn",
            "-query",
            "not_here.fa",
            "-subject",
            "ns.fa",
            "-logfile",
            "not_here.fa",
        ],
        &[],
    );
    assert_eq!(out_missing.status.code(), Some(1));
    assert!(
        !missing.exists(),
        "missing input file must not be created by -logfile"
    );

    // Non-self-referencing -logfile works normally (and still does not touch inputs).
    let ok_log = dir.join("normal.log");
    let out_ok = run_in(
        &dir,
        &[
            "blastn",
            "-query",
            "nq.fa",
            "-subject",
            "ns.fa",
            "-outfmt",
            "6",
            "-logfile",
            "normal.log",
        ],
        &[("NCBI_DONT_USE_NCBIRC", "1")],
    );
    assert_eq!(
        out_ok.status.code(),
        Some(0),
        "normal -logfile failed: {}",
        String::from_utf8_lossy(&out_ok.stderr)
    );
    assert!(ok_log.exists(), "normal -logfile should be created");

    // `-version` does not read command-line input files (`-query`), so `-logfile nq.fa` is not
    // refused by argv input collection (still succeeds with rc 0).
    let version_q = dir.join("ver_q.fa");
    fs::write(&version_q, NUCL_FASTA).unwrap();
    let out_ver = run_in(
        &dir,
        &[
            "blastn", "-version", "-query", "ver_q.fa", "-logfile", "ver_q.fa",
        ],
        &[("NCBI_DONT_USE_NCBIRC", "1")],
    );
    assert_eq!(
        out_ver.status.code(),
        Some(0),
        "blastn -version -query ver_q.fa -logfile ver_q.fa failed: {}",
        String::from_utf8_lossy(&out_ver.stderr)
    );

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn tools_refuse_self_referencing_logfile() {
    let dir = fresh_dir("tools");
    let fa = dir.join("in.fa");
    let counts = dir.join("u.counts");
    let batch = dir.join("batch.txt");
    let taxmap = dir.join("tax.map");
    let mask = dir.join("mask.asnb");
    let asn = dir.join("in.asn");

    fs::write(&fa, NUCL_FASTA).unwrap();
    fs::write(&counts, b"not counts\n").unwrap();
    fs::write(&batch, b"seq1\n").unwrap();
    fs::write(&taxmap, b"seq1 9606\n").unwrap();
    fs::write(&mask, b"mask\n").unwrap();
    fs::write(&asn, b"Blast4-archive ::= {}\n").unwrap();

    let fa_snap = snapshot(&fa);
    let counts_snap = snapshot(&counts);
    let batch_snap = snapshot(&batch);
    let taxmap_snap = snapshot(&taxmap);
    let mask_snap = snapshot(&mask);
    let asn_snap = snapshot(&asn);

    // makeblastdb: -in, -taxid_map, -mask_data, and leading-space raw vs normalized.
    for (label, argv, (path, snap)) in [
        (
            "makeblastdb -in",
            vec![
                "makeblastdb",
                "-in",
                "in.fa",
                "-dbtype",
                "nucl",
                "-out",
                "db",
                "-logfile",
                "in.fa",
            ],
            (&fa, &fa_snap),
        ),
        (
            "makeblastdb leading-space -logfile",
            vec![
                "makeblastdb",
                "-in",
                "in.fa",
                "-dbtype",
                "nucl",
                "-out",
                "db",
                "-logfile",
                " in.fa",
            ],
            (&fa, &fa_snap),
        ),
        (
            "makeblastdb -taxid_map",
            vec![
                "makeblastdb",
                "-in",
                "in.fa",
                "-dbtype",
                "nucl",
                "-parse_seqids",
                "-taxid_map",
                "tax.map",
                "-out",
                "db",
                "-logfile",
                "tax.map",
            ],
            (&taxmap, &taxmap_snap),
        ),
        (
            "makeblastdb -mask_data",
            vec![
                "makeblastdb",
                "-in",
                "in.fa",
                "-dbtype",
                "nucl",
                "-mask_data",
                "mask.asnb",
                "-out",
                "db",
                "-logfile",
                "mask.asnb",
            ],
            (&mask, &mask_snap),
        ),
    ] {
        let out = run_in(&dir, &argv, &[]);
        assert_refused_and_unchanged(&out, &[(path.as_path(), snap)], label);
    }

    // Build a real small DB so we can test blastdbcmd -db component and -entry_batch.
    let build = run_in(
        &dir,
        &[
            "makeblastdb",
            "-in",
            "in.fa",
            "-dbtype",
            "nucl",
            "-out",
            "mydb",
            "-logfile",
            "mk.log",
        ],
        &[],
    );
    assert_eq!(build.status.code(), Some(0));
    let nhr = dir.join("mydb.nhr");
    let nhr_snap = snapshot(&nhr);

    let out_db_comp = run_in(
        &dir,
        &[
            "blastn", "-query", "in.fa", "-db", "mydb", "-logfile", "mydb.nhr",
        ],
        &[],
    );
    assert_refused_and_unchanged(&out_db_comp, &[(&nhr, &nhr_snap)], "blastn -db component");

    // Masked nucleotide DB: `-db_soft_mask 11 -logfile masked_db.naa` must refuse before modifying
    // the mask column index `.naa`, while a distinct logfile succeeds.
    let dust_asnb = run_in(
        &dir,
        &[
            "dustmasker",
            "-in",
            "in.fa",
            "-outfmt",
            "maskinfo_asn1_bin",
            "-out",
            "dust.asnb",
        ],
        &[],
    );
    assert_eq!(dust_asnb.status.code(), Some(0));
    let build_masked = run_in(
        &dir,
        &[
            "makeblastdb",
            "-in",
            "in.fa",
            "-dbtype",
            "nucl",
            "-mask_data",
            "dust.asnb",
            "-out",
            "masked_db",
            "-logfile",
            "mk_masked.log",
        ],
        &[],
    );
    assert_eq!(build_masked.status.code(), Some(0));
    let naa = dir.join("masked_db.naa");
    let naa_snap = snapshot(&naa);
    let out_db_mask = run_in(
        &dir,
        &[
            "blastn",
            "-query",
            "in.fa",
            "-db",
            "masked_db",
            "-db_soft_mask",
            "11",
            "-logfile",
            "masked_db.naa",
        ],
        &[],
    );
    assert_refused_and_unchanged(
        &out_db_mask,
        &[(&naa, &naa_snap)],
        "blastn -db_soft_mask .naa",
    );
    let out_db_mask_ok = run_in(
        &dir,
        &[
            "blastn",
            "-query",
            "in.fa",
            "-db",
            "masked_db",
            "-db_soft_mask",
            "11",
            "-outfmt",
            "6",
            "-logfile",
            "masked_ok.log",
        ],
        &[],
    );
    assert_eq!(
        out_db_mask_ok.status.code(),
        Some(0),
        "blastn -db_soft_mask with distinct logfile failed: {}",
        String::from_utf8_lossy(&out_db_mask_ok.stderr)
    );

    // Alias pointing to `mydb` via DBLIST: both target volume (`mydb.nhr`) and alias file itself
    // (`myalias.nal`) must be refused, including when resolved through `BLASTDB`.
    let alias_nal = dir.join("myalias.nal");
    fs::write(&alias_nal, b"TITLE my alias\nDBLIST mydb\n").unwrap();
    let alias_snap = snapshot(&alias_nal);

    let out_alias_vol = run_in(
        &dir,
        &[
            "blastn", "-query", "in.fa", "-db", "myalias", "-logfile", "mydb.nhr",
        ],
        &[],
    );
    assert_refused_and_unchanged(
        &out_alias_vol,
        &[(&nhr, &nhr_snap), (&alias_nal, &alias_snap)],
        "blastn -db alias -> volume component",
    );

    let out_alias_self = run_in(
        &dir,
        &[
            "blastn",
            "-query",
            "in.fa",
            "-db",
            "myalias",
            "-logfile",
            "myalias.nal",
        ],
        &[],
    );
    assert_refused_and_unchanged(
        &out_alias_self,
        &[(&alias_nal, &alias_snap)],
        "blastn -db alias self",
    );

    let caller_dir = dir.join("caller");
    fs::create_dir_all(&caller_dir).unwrap();
    let out_blastdb_env = run_in(
        &caller_dir,
        &[
            "blastn",
            "-query",
            &fa.to_string_lossy(),
            "-db",
            "myalias",
            "-logfile",
            &nhr.to_string_lossy(),
        ],
        &[("BLASTDB", &dir.to_string_lossy())],
    );
    assert_refused_and_unchanged(
        &out_blastdb_env,
        &[(&nhr, &nhr_snap)],
        "blastn -db alias via BLASTDB",
    );

    let out_bdc_batch = run_in(
        &dir,
        &[
            "blastdbcmd",
            "-db",
            "mydb",
            "-entry_batch",
            "batch.txt",
            "-logfile",
            "batch.txt",
        ],
        &[],
    );
    assert_refused_and_unchanged(
        &out_bdc_batch,
        &[(&batch, &batch_snap)],
        "blastdbcmd -entry_batch",
    );

    let out_bdc_db = run_in(
        &dir,
        &["blastdbcmd", "-db", "mydb", "-info", "-logfile", "mydb.nhr"],
        &[],
    );
    assert_refused_and_unchanged(
        &out_bdc_db,
        &[(&nhr, &nhr_snap)],
        "blastdbcmd -db component",
    );

    // windowmasker: -in and -ustat.
    let out_wm_in = run_in(
        &dir,
        &[
            "windowmasker",
            "-mk_counts",
            "-in",
            "in.fa",
            "-out",
            "out.counts",
            "-logfile",
            "in.fa",
        ],
        &[],
    );
    assert_refused_and_unchanged(&out_wm_in, &[(&fa, &fa_snap)], "windowmasker -in");

    let out_wm_ustat = run_in(
        &dir,
        &[
            "windowmasker",
            "-ustat",
            "u.counts",
            "-in",
            "in.fa",
            "-logfile",
            "u.counts",
        ],
        &[],
    );
    assert_refused_and_unchanged(
        &out_wm_ustat,
        &[(&counts, &counts_snap)],
        "windowmasker -ustat",
    );

    // dustmasker, segmasker, convert2blastmask, makembindex, rpsbproc.
    for (label, argv, watched) in [
        (
            "dustmasker -in",
            vec!["dustmasker", "-in", "in.fa", "-logfile", "in.fa"],
            (&fa, &fa_snap),
        ),
        (
            "segmasker -in",
            vec!["segmasker", "-in", "in.fa", "-logfile", "in.fa"],
            (&fa, &fa_snap),
        ),
        (
            "convert2blastmask -in",
            vec!["convert2blastmask", "-in", "in.fa", "-logfile", "in.fa"],
            (&fa, &fa_snap),
        ),
        (
            "makembindex -input",
            vec![
                "makembindex",
                "-iformat",
                "fasta",
                "-input",
                "in.fa",
                "-logfile",
                "in.fa",
            ],
            (&fa, &fa_snap),
        ),
        (
            "rpsbproc -i",
            vec!["rpsbproc", "-i", "in.asn", "-logfile", "in.asn"],
            (&asn, &asn_snap),
        ),
    ] {
        let out = run_in(&dir, &argv, &[]);
        assert_refused_and_unchanged(&out, &[(watched.0.as_path(), watched.1)], label);
    }

    // Normal non-self-referencing -logfile for dustmasker works normally.
    let out_dust_ok = run_in(
        &dir,
        &["dustmasker", "-in", "in.fa", "-logfile", "dust.log"],
        &[],
    );
    assert_eq!(out_dust_ok.status.code(), Some(0));
    assert!(dir.join("dust.log").exists());

    let _ = fs::remove_dir_all(&dir);
}

// Line WQA D26 round 1: a name that the run cannot open is the run's log, as NCBI writes it. The
// mask column of a masked database is read with `-db_soft_mask` only, the id-lookup files with an
// id list or a taxonomy option, the RPS files by `rpsblast` and `rpstblastn`: a plain `blastn` or
// `blastp` search may name them as its `-logfile` (the refusals with the option are the tests
// above).
#[test]
fn database_files_that_the_run_cannot_open_are_the_log() {
    let dir = fresh_dir("d26-unopened");
    fs::write(dir.join("in.fa"), NUCL_FASTA).unwrap();
    fs::write(dir.join("prot.fa"), PROT_FASTA).unwrap();
    let ok = |args: &[&str], label: &str| {
        let out = run_in(&dir, args, &[]);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert_eq!(out.status.code(), Some(0), "{label}: {stderr}");
        assert!(!stderr.contains("-logfile file"), "{label}: {stderr}");
    };
    ok(
        &[
            "makeblastdb",
            "-in",
            "in.fa",
            "-dbtype",
            "nucl",
            "-out",
            "nd",
            "-logfile",
            "nd.mk",
        ],
        "makeblastdb nucl",
    );
    ok(
        &[
            "makeblastdb",
            "-in",
            "prot.fa",
            "-dbtype",
            "prot",
            "-out",
            "pd",
            "-logfile",
            "pd.mk",
        ],
        "makeblastdb prot",
    );
    // Not opened by this search: the mask column and the GI index of the database, the id-lookup
    // files, the metadata.
    for name in ["nd.naa", "nd.nog", "nd.nnd", "nd.nos", "nd.nto"] {
        ok(
            &[
                "blastn", "-query", "in.fa", "-db", "nd", "-outfmt", "6", "-logfile", name,
            ],
            name,
        );
    }
    // The RPS files belong to the RPS searches.
    for name in ["pd.rps", "pd.loo", "pd.aux", "pd.freq", "pd.pnd"] {
        ok(
            &[
                "blastp", "-query", "prot.fa", "-db", "pd", "-outfmt", "6", "-logfile", name,
            ],
            name,
        );
    }
    // The files that the same searches do open stay refused.
    for (args, name) in [
        (
            [
                "blastn", "-query", "in.fa", "-db", "nd", "-logfile", "nd.nsq",
            ],
            "nd.nsq",
        ),
        (
            [
                "blastp", "-query", "prot.fa", "-db", "pd", "-logfile", "pd.phr",
            ],
            "pd.phr",
        ),
    ] {
        let before = fs::read(dir.join(name)).unwrap();
        let out = run_in(&dir, &args, &[]);
        assert_eq!(out.status.code(), Some(1), "{name}");
        assert_eq!(fs::read(dir.join(name)).unwrap(), before, "{name}");
    }
    let _ = fs::remove_dir_all(&dir);
}

// Line WQA D35 (WQA DECISIONS 62): the taxonomy files that a run looks for along its search path
// (the working directory first, then `BLASTDB`) are files it reads, as the native open trace shows
// them: `taxdb.bti` up to the first found and the `taxdb.btd` beside it, for every search of a
// database whatever its format and for a `-subject` search whose report looks every taxid up;
// `taxonomy4blast.sqlite3` likewise for a taxonomy filter that is widened. A `-logfile` of such a
// name would be found by the run in place of the user's file (the names of outfmt 18 lost), so it
// is refused, before the file is made; a name that the run does not look at is the run's log.
#[test]
fn taxonomy_files_that_the_run_looks_for_are_refused() {
    let dir = fresh_dir("d35-taxonomy");
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let folder = dir.join("protdb");
    fs::create_dir_all(&folder).unwrap();
    for entry in fs::read_dir(tests.join("db_search/taxfilter")).unwrap() {
        let path = entry.unwrap().path();
        let name = path.file_name().unwrap().to_str().unwrap().to_string();
        if name.starts_with("tfp.") || name == "taxonomy4blast.sqlite3" {
            fs::copy(&path, folder.join(&name)).unwrap();
        }
    }
    for name in ["taxdb.bti", "taxdb.btd"] {
        fs::copy(tests.join("seqdb/taxdb/good").join(name), folder.join(name)).unwrap();
    }
    fs::copy(tests.join("db_search/taxfilter/q.faa"), dir.join("q.faa")).unwrap();
    let blastdb = [("BLASTDB", "protdb")];
    let search = |extra: &[&str], log: &str| -> Vec<String> {
        let mut args: Vec<String> = ["blastp", "-query", "q.faa", "-num_threads", "1"]
            .iter()
            .map(|word| word.to_string())
            .collect();
        args.extend(extra.iter().map(|word| word.to_string()));
        args.extend(["-logfile".to_string(), log.to_string()]);
        args
    };
    let refused = |args: &[String], made: &str, label: &str| {
        let words: Vec<&str> = args.iter().map(String::as_str).collect();
        let out = run_in(&dir, &words, &blastdb);
        assert_refused_and_unchanged(&out, &[], label);
        assert!(!dir.join(made).exists(), "{label}: {made} was made");
    };
    let runs = |args: &[String], made: &str, label: &str| -> Vec<u8> {
        let words: Vec<&str> = args.iter().map(String::as_str).collect();
        let out = run_in(&dir, &words, &blastdb);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert_eq!(out.status.code(), Some(0), "{label}: {stderr}");
        let _ = fs::remove_file(dir.join(made));
        out.stdout
    };
    let db = ["-db", "protdb/tfp"];
    // A search of a database: every format looks for `taxdb.bti` in the working directory first.
    for fmt in [
        "0",
        "6",
        "7",
        "6 qseqid evalue",
        "8",
        "11",
        "18",
        "6 std sscinames",
    ] {
        let args = search(&[&db[..], &["-outfmt", fmt]].concat(), "taxdb.bti");
        refused(&args, "taxdb.bti", &format!("-db outfmt {fmt}: taxdb.bti"));
    }
    // The data file beside the index found in `BLASTDB`, and the user's own files there.
    for name in ["protdb/taxdb.btd", "protdb/taxdb.bti"] {
        let before = snapshot(&dir.join(name));
        let args = search(&[&db[..], &["-outfmt", "6"]].concat(), name);
        let words: Vec<&str> = args.iter().map(String::as_str).collect();
        let out = run_in(&dir, &words, &blastdb);
        assert_refused_and_unchanged(&out, &[(&dir.join(name), &before)], name);
    }
    // `taxdb.btd` in the working directory is not looked at (the index is found in `BLASTDB`).
    runs(
        &search(&[&db[..], &["-outfmt", "6"]].concat(), "taxdb.btd"),
        "taxdb.btd",
        "-db taxdb.btd",
    );
    // A widened taxonomy filter looks for `taxonomy4blast.sqlite3` the same way.
    let taxids = [&db[..], &["-outfmt", "6", "-taxids", "9606"]].concat();
    refused(
        &search(&taxids, "taxonomy4blast.sqlite3"),
        "taxonomy4blast.sqlite3",
        "-taxids: taxonomy4blast.sqlite3",
    );
    let kept = [&taxids[..], &["-no_taxid_expansion"]].concat();
    runs(
        &search(&kept, "taxonomy4blast.sqlite3"),
        "taxonomy4blast.sqlite3",
        "-no_taxid_expansion: taxonomy4blast.sqlite3",
    );
    runs(
        &search(
            &[&db[..], &["-outfmt", "6"]].concat(),
            "taxonomy4blast.sqlite3",
        ),
        "taxonomy4blast.sqlite3",
        "no taxonomy filter: taxonomy4blast.sqlite3",
    );
    // A `-subject` search looks for `taxdb.bti` only for a report that looks every taxid up; the
    // others run with that log and write the bytes they write with another.
    let subject = ["-subject", "q.faa"];
    for fmt in ["18", "6 std sscinames"] {
        let args = search(&[&subject[..], &["-outfmt", fmt]].concat(), "taxdb.bti");
        refused(
            &args,
            "taxdb.bti",
            &format!("-subject outfmt {fmt}: taxdb.bti"),
        );
    }
    for fmt in ["6", "0", "6 std ssciname staxids"] {
        let format = [&subject[..], &["-outfmt", fmt]].concat();
        let with_bti = runs(&search(&format, "taxdb.bti"), "taxdb.bti", fmt);
        let with_log = runs(&search(&format, "run.log"), "run.log", fmt);
        assert_eq!(with_bti, with_log, "-subject outfmt {fmt}");
    }
    let _ = fs::remove_dir_all(&dir);
}

// Line WQA D26 item 7 (D19b review follow-up): the entries that `windowmasker -mk_counts
// -fa_list true` reads from its list file are inputs of the run, as the list file is.
#[test]
fn windowmasker_fa_list_entries_are_refused() {
    let dir = fresh_dir("wm-fa-list");
    let list = dir.join("list.txt");
    let entry = dir.join("m.fa");
    fs::write(&list, b"m.fa\n").unwrap();
    fn argv(log: Option<&str>) -> Vec<&str> {
        let mut argv = vec![
            "windowmasker",
            "-mk_counts",
            "-fa_list",
            "true",
            "-in",
            "list.txt",
            "-out",
            "c.out",
        ];
        if let Some(log) = log {
            argv.extend(["-logfile", log]);
        }
        argv
    }
    let list_snap = snapshot(&list);

    // A missing entry: refused before `-logfile` creates it.
    let out = run_in(&dir, &argv(Some("m.fa")), &[]);
    assert_refused_and_unchanged(
        &out,
        &[(&list, &list_snap)],
        "windowmasker -fa_list missing entry",
    );
    assert!(!entry.exists(), "m.fa created");
    assert!(!dir.join("c.out").exists(), "c.out created");

    // An existing entry: refused, its bytes and mtime unchanged.
    fs::write(&entry, NUCL_FASTA).unwrap();
    let entry_snap = snapshot(&entry);
    let out = run_in(&dir, &argv(Some("m.fa")), &[]);
    assert_refused_and_unchanged(
        &out,
        &[(&entry, &entry_snap), (&list, &list_snap)],
        "windowmasker -fa_list entry",
    );

    // The list file itself: refused.
    let out = run_in(&dir, &argv(Some("list.txt")), &[]);
    assert_refused_and_unchanged(
        &out,
        &[(&entry, &entry_snap), (&list, &list_snap)],
        "windowmasker -fa_list list file",
    );

    // Without `-logfile` the run counts the entry.
    let out = run_in(&dir, &argv(None), &[]);
    assert_eq!(
        out.status.code(),
        Some(0),
        "windowmasker -fa_list without -logfile: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(snapshot(&entry), entry_snap, "m.fa changed");
    assert!(dir.join("c.out").exists(), "c.out not written");

    let _ = fs::remove_dir_all(&dir);
}

#[cfg(unix)]
#[test]
fn symlinks_and_hardlinks_to_inputs_are_refused() {
    use std::os::unix::fs::symlink;

    let dir = fresh_dir("symlinks");
    fs::create_dir_all(dir.join("real")).unwrap();
    let q = dir.join("real").join("q.fa");
    let s = dir.join("s.fa");
    fs::write(&q, NUCL_FASTA).unwrap();
    fs::write(&s, NUCL_FASTA).unwrap();
    let q_snap = snapshot(&q);

    symlink(&q, dir.join("q_sym.fa")).unwrap();
    symlink(dir.join("real"), dir.join("dir_sym")).unwrap();
    fs::hard_link(&q, dir.join("q_hard.fa")).unwrap();

    for (label, log_name) in [
        ("symlink to query", "q_sym.fa"),
        ("symlinked dir to query", "dir_sym/q.fa"),
        ("hard link to query", "q_hard.fa"),
    ] {
        let out = run_in(
            &dir,
            &[
                "tblastx",
                "-query",
                "real/q.fa",
                "-subject",
                "s.fa",
                "-logfile",
                log_name,
            ],
            &[],
        );
        assert_refused_and_unchanged(&out, &[(&q, &q_snap)], label);
    }

    // Reversed: query is a symlink, logfile is the real file.
    let out_rev = run_in(
        &dir,
        &[
            "tblastx",
            "-query",
            "q_sym.fa",
            "-subject",
            "s.fa",
            "-logfile",
            "real/q.fa",
        ],
        &[],
    );
    assert_refused_and_unchanged(&out_rev, &[(&q, &q_snap)], "query symlink vs real logfile");

    // Existing directory as -logfile falls back to stderr (not refused).
    let out_dir = run_in(
        &dir,
        &[
            "blastn",
            "-query",
            "real/q.fa",
            "-subject",
            "s.fa",
            "-outfmt",
            "6",
            "-logfile",
            "real",
        ],
        &[],
    );
    assert_eq!(out_dir.status.code(), Some(0));

    let _ = fs::remove_dir_all(&dir);
}
