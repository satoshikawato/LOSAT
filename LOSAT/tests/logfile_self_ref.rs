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

/// Copies NCBI's strategy exports of `tests/fixtures/strategy/import/ss/` and the databases that
/// they name (`tests/db_search/<name>.*`) into `dir`, at the paths that the strategies name.
fn strategy_inputs(dir: &Path, strategies: &[&str], databases: &[&str]) {
    let engine = Path::new(env!("CARGO_MANIFEST_DIR"));
    let db_dir = dir.join("tests/db_search");
    fs::create_dir_all(&db_dir).unwrap();
    for entry in fs::read_dir(engine.join("tests/db_search")).unwrap() {
        let entry = entry.unwrap();
        let file = entry.file_name();
        let file_text = file.to_string_lossy();
        if databases
            .iter()
            .any(|name| file_text.starts_with(&format!("{name}.")))
            && entry.path().is_file()
        {
            fs::copy(entry.path(), db_dir.join(&file)).unwrap();
        }
    }
    for strategy in strategies {
        fs::copy(
            engine
                .join("tests/fixtures/strategy/import/ss")
                .join(strategy),
            dir.join(strategy),
        )
        .unwrap();
    }
}

/// Line RG-F item F2f: a `-logfile` naming a file that an imported strategy makes the run read
/// (its database, opened as the type of its program's subjects; the repeat filter's database; the
/// window masker's counts of its taxid; its MegaBLAST index, new-style and old-style) is refused
/// as the command line's inputs are, before the file is opened (NCBI's application opens it before
/// its run reads the strategy). A strategy whose database the command line overrides, a strategy
/// read from the standard input (read once, by the import) and a file that is not a strategy (the
/// import's own error, unchanged) add nothing.
#[test]
fn a_logfile_naming_what_an_imported_strategy_reads_is_refused() {
    let dir = fresh_dir("strategy");
    strategy_inputs(
        &dir,
        &[
            "blastp.db.ss",
            "tblastx.db.ss",
            "bn_rpt.ss",
            "blastn.useidx_force.ss",
            "s3a.blastn_WindowMaskerTaxId_0.ss",
        ],
        &[
            "arch_sprot",
            "arch_tbx_snucl",
            "arch_snucl",
            "rpt_d4",
            "wm_s60k",
        ],
    );
    // Hand-edited strategies: a forced index of another name than the database, and a window
    // masker taxid whose counts file exists in the run's folder (`[WINDOW_MASKER]` unset: `.`).
    let edit = |from: &str, to: &str, old: &str, new: &str| {
        let text = fs::read_to_string(dir.join(from)).unwrap();
        assert!(text.contains(old), "{from}: {old}");
        fs::write(dir.join(to), text.replace(old, new)).unwrap();
    };
    edit(
        "blastn.useidx_force.ss",
        "mbidx.ss",
        "value string \"tests/db_search/arch_snucl\"",
        "value string \"tests/db_search/rpt_d4\"",
    );
    edit(
        "s3a.blastn_WindowMaskerTaxId_0.ss",
        "wm9606.ss",
        "name \"WindowMaskerTaxId\",\n        value integer 0",
        "name \"WindowMaskerTaxId\",\n        value integer 9606",
    );
    fs::create_dir_all(dir.join("9606")).unwrap();
    fs::write(dir.join("9606/wmasker.obinary"), b"counts\n").unwrap();
    fs::write(dir.join("ps.fa"), PROT_FASTA).unwrap();
    fs::write(dir.join("junk.ss"), b"not a strategy\n").unwrap();

    let refused: &[(&str, &str, &str, &str)] = &[
        (
            "db",
            "blastp",
            "blastp.db.ss",
            "tests/db_search/arch_sprot.pin",
        ),
        (
            "db sequences",
            "blastp",
            "blastp.db.ss",
            "tests/db_search/arch_sprot.psq",
        ),
        // blastx running a TBLASTX strategy searches its nucleotide database.
        (
            "nucleotide db of blastx",
            "blastx",
            "tblastx.db.ss",
            "tests/db_search/arch_tbx_snucl.nin",
        ),
        (
            "nucleotide db of blastn",
            "blastn",
            "tblastx.db.ss",
            "tests/db_search/arch_tbx_snucl.nsq",
        ),
        (
            "repeat filter db",
            "blastn",
            "bn_rpt.ss",
            "tests/db_search/rpt_d4.nhr",
        ),
        (
            "new-style index",
            "blastn",
            "mbidx.ss",
            "tests/db_search/rpt_d4.nin",
        ),
        (
            "old-style index",
            "blastn",
            "mbidx.ss",
            "tests/db_search/rpt_d4.00.idx",
        ),
        (
            "window masker taxid",
            "blastn",
            "wm9606.ss",
            "9606/wmasker.obinary",
        ),
    ];
    for &(label, program, strategy, log_name) in refused {
        let watched = dir.join(log_name);
        let before = watched.exists().then(|| snapshot(&watched));
        let out = run_in(
            &dir,
            &[
                program,
                "-import_search_strategy",
                strategy,
                "-logfile",
                log_name,
            ],
            &[],
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            stderr.contains(&format!(
                "Error: a -logfile file '{log_name}' that is also an input file"
            )),
            "{label}: {stderr}"
        );
        match &before {
            Some(snap) => assert_refused_and_unchanged(&out, &[(&watched, snap)], label),
            None => {
                assert_refused_and_unchanged(&out, &[], label);
                assert!(!watched.exists(), "{label}: {log_name} was made");
            }
        }
    }

    // The command line's `-subject` overrides the strategy's database: not refused, and the
    // search runs.
    let out = run_in(
        &dir,
        &[
            "blastp",
            "-import_search_strategy",
            "blastp.db.ss",
            "-subject",
            "ps.fa",
            "-outfmt",
            "6",
            "-logfile",
            "tests/db_search/arch_sprot.pin",
        ],
        &[],
    );
    assert_eq!(
        out.status.code(),
        Some(0),
        "override: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    // A strategy read from the standard input is read once, by the import: the check does not
    // read it, and the run equals the one that names the file.
    let named = run_in(
        &dir,
        &[
            "blastp",
            "-import_search_strategy",
            "blastp.db.ss",
            "-outfmt",
            "6",
            "-logfile",
            "named.log",
        ],
        &[],
    );
    assert_eq!(
        named.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&named.stderr)
    );
    assert!(!named.stdout.is_empty());
    let mut cmd = Command::new(BIN);
    cmd.current_dir(&dir)
        .args([
            "blastp",
            "-import_search_strategy",
            "-",
            "-outfmt",
            "6",
            "-logfile",
            "stdin.log",
        ])
        .env_remove("BLASTDB")
        .env_remove("BLAST_USAGE_REPORT")
        .env("HOME", &dir)
        .stdin(fs::File::open(dir.join("blastp.db.ss")).unwrap());
    let piped = cmd.output().unwrap();
    assert_eq!(
        piped.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&piped.stderr)
    );
    assert_eq!(piped.stdout, named.stdout);
    assert_eq!(
        fs::read(dir.join("stdin.log")).unwrap(),
        fs::read(dir.join("named.log")).unwrap()
    );

    // A file that is not a strategy: the import's error, byte for byte as with another
    // `-logfile`, which gets it.
    let pin = dir.join("tests/db_search/arch_sprot.pin");
    let pin_before = fs::read(&pin).unwrap();
    let junk = |log_name: &str| {
        run_in(
            &dir,
            &[
                "blastp",
                "-import_search_strategy",
                "junk.ss",
                "-logfile",
                log_name,
            ],
            &[],
        )
    };
    let other = junk("junk.log");
    let at_db = junk("tests/db_search/arch_sprot.pin");
    assert_eq!(at_db.status.code(), other.status.code());
    assert_ne!(other.status.code(), Some(0));
    assert_eq!(at_db.stderr, other.stderr);
    assert_eq!(at_db.stdout, other.stdout);
    let appended = fs::read(&pin).unwrap()[pin_before.len()..].to_vec();
    assert_eq!(appended, fs::read(dir.join("junk.log")).unwrap());

    let _ = fs::remove_dir_all(&dir);
}

/// Line RG-F item F2f: the library API's inputs that its argv cannot name count as the command
/// line's do: the database of `Inputs::database` (the command line's `-db`) and what the strategy
/// of `Inputs::strategy` (`-import_search_strategy`) makes the run read (its import's names). A
/// `-logfile` naming a file of either database is refused before it is opened, and the file is
/// unchanged; a log of its own runs.
#[test]
fn the_library_apis_inputs_count_as_files_the_run_reads() {
    use losat::api::{read_inputs, run, validate, Error, Inputs, Program};
    let dir = fresh_dir("api_inputs");
    strategy_inputs(&dir, &["blastp.db.ss"], &["arch_sprot"]);
    let db = dir.join("tests/db_search/arch_sprot");
    let db_name = db.to_str().unwrap();
    // The strategy's database, named by the copy's whole path (the test's working directory is
    // the process's).
    let text = fs::read_to_string(dir.join("blastp.db.ss")).unwrap();
    assert!(text.contains("\"tests/db_search/arch_sprot\""));
    let text = text.replace("\"tests/db_search/arch_sprot\"", &format!("\"{db_name}\""));
    let run_with = |log: &Path, strategy: bool| -> Result<(), Error> {
        let search = validate(Program::Blastp, &["-logfile", log.to_str().unwrap()]).unwrap();
        let (mut report, mut diagnostics) = (Vec::new(), Vec::new());
        if strategy {
            let strategy = search.read_strategy(text.as_bytes())?;
            run(
                &search,
                Inputs::strategy(&strategy),
                &mut [("6", &mut report)],
                &mut diagnostics,
            )
        } else {
            let (queries, _) = read_inputs(&search, PROT_FASTA, b"").unwrap();
            run(
                &search,
                Inputs::database(&queries, db_name),
                &mut [("6", &mut report)],
                &mut diagnostics,
            )
        }
    };
    for file in ["arch_sprot.pin", "arch_sprot.psq"] {
        let log = dir.join("tests/db_search").join(file);
        let before = snapshot(&log);
        for strategy in [false, true] {
            let error = run_with(&log, strategy).unwrap_err();
            assert!(
                format!("{error:?}").contains("that is also an input file read by this run"),
                "{file} strategy {strategy}: {error:?}"
            );
            assert_eq!(snapshot(&log), before, "{file} strategy {strategy}");
        }
    }
    for strategy in [false, true] {
        let log = dir.join(format!("api_{strategy}.log"));
        run_with(&log, strategy).unwrap();
        assert!(log.is_file(), "strategy {strategy}");
    }
    let _ = fs::remove_dir_all(&dir);
}
