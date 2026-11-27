use std::ffi::{OsStr, OsString};
use std::path::PathBuf;

use super::*;
use crate::blastinput::ncbi_environment::ApplicationStart;

fn fresh_dir(tag: &str) -> PathBuf {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-log-inputs-{tag}-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn skips_option_values_that_look_like_input_flags() {
    let dir = fresh_dir("opt-value");
    let q = dir.join("q.fa");
    let s = dir.join("s.fa");
    std::fs::write(&q, b">q\nACGT\n").unwrap();
    std::fs::write(&s, b">s\nACGT\n").unwrap();

    let words: Vec<OsString> = vec![
        "-export_search_strategy".into(),
        "-query".into(),
        "-query".into(),
        q.as_os_str().to_os_string(),
        "-subject".into(),
        s.as_os_str().to_os_string(),
        "-logfile".into(),
        q.as_os_str().to_os_string(),
    ];
    let start = ApplicationStart::of_process("tblastx", &words);
    assert!(check_disk_logfile("tblastx", &words, &start, Some(q.as_os_str())).is_err());

    // Passing -logfile as a distinct file passes.
    let log = dir.join("run.log");
    assert!(check_disk_logfile("tblastx", &words, &start, Some(log.as_os_str())).is_ok());

    // Special targets always pass without error.
    for special in [
        None,
        Some(OsStr::new("-")),
        Some(OsStr::new("")),
        Some(OsStr::new("/dev/null")),
    ] {
        assert!(check_disk_logfile("tblastx", &words, &start, special).is_ok());
    }

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn collects_search_and_tool_inputs() {
    let dir = fresh_dir("collect");
    let q = dir.join("q.fa");
    let pssm = dir.join("p.smp");
    let msa = dir.join("m.fa");
    let gi = dir.join("gi.txt");
    let seqid = dir.join("seq.txt");
    let taxid = dir.join("tax.txt");
    let strat = dir.join("strat.asn");
    let search_files = [&q, &pssm, &msa, &gi, &seqid, &taxid, &strat];
    for f in search_files {
        std::fs::write(f, b"x").unwrap();
    }

    let words: Vec<OsString> = vec![
        "-query".into(),
        q.as_os_str().to_os_string(),
        "-in_pssm".into(),
        pssm.as_os_str().to_os_string(),
        "-in_msa".into(),
        msa.as_os_str().to_os_string(),
        "-gilist".into(),
        gi.as_os_str().to_os_string(),
        "-seqidlist".into(),
        seqid.as_os_str().to_os_string(),
        "-taxidlist".into(),
        taxid.as_os_str().to_os_string(),
        "-import_search_strategy".into(),
        strat.as_os_str().to_os_string(),
    ];
    let start = ApplicationStart::of_process("psiblast", &words);
    for f in search_files {
        assert!(
            check_disk_logfile("psiblast", &words, &start, Some(f.as_os_str())).is_err(),
            "expected refusal for {}",
            f.display()
        );
    }

    // windowmasker -ustat / -in / -ids / -exclude_ids
    let wm_in = dir.join("x.fa");
    let wm_ustat = dir.join("u.counts");
    let wm_ids = dir.join("ids.txt");
    let wm_ex = dir.join("ex.txt");
    let wm_files = [&wm_in, &wm_ustat, &wm_ids, &wm_ex];
    for f in wm_files {
        std::fs::write(f, b"x").unwrap();
    }
    let wm_words: Vec<OsString> = vec![
        "-in".into(),
        wm_in.as_os_str().to_os_string(),
        "-ustat".into(),
        wm_ustat.as_os_str().to_os_string(),
        "-ids".into(),
        wm_ids.as_os_str().to_os_string(),
        "-exclude_ids".into(),
        wm_ex.as_os_str().to_os_string(),
    ];
    let wm_start = ApplicationStart::of_process("windowmasker", &wm_words);
    for f in wm_files {
        assert!(
            check_disk_logfile("windowmasker", &wm_words, &wm_start, Some(f.as_os_str())).is_err(),
            "expected windowmasker refusal for {}",
            f.display()
        );
    }

    // makeblastdb comma-separated -mask_data and space-separated -in
    let mk_a = dir.join("a.fa");
    let mk_b = dir.join("b.fa");
    let mk_m1 = dir.join("m1.asnb");
    let mk_m2 = dir.join("m2.asnb");
    let mk_tm = dir.join("tax.map");
    let mk_files = [&mk_a, &mk_b, &mk_m1, &mk_m2, &mk_tm];
    for f in mk_files {
        std::fs::write(f, b"x").unwrap();
    }
    let mk_in = format!("\"{}\" \"{}\"", mk_a.display(), mk_b.display());
    let mk_mask = format!("{},{}", mk_m1.display(), mk_m2.display());
    let mk_words: Vec<OsString> = vec![
        "-in".into(),
        mk_in.into(),
        "-mask_data".into(),
        mk_mask.into(),
        "-taxid_map".into(),
        mk_tm.as_os_str().to_os_string(),
    ];
    let mk_start = ApplicationStart::of_process("makeblastdb", &mk_words);
    for f in mk_files {
        assert!(
            check_disk_logfile("makeblastdb", &mk_words, &mk_start, Some(f.as_os_str())).is_err(),
            "expected makeblastdb refusal for {}",
            f.display()
        );
    }

    let _ = std::fs::remove_dir_all(&dir);
}

fn blastn_start(words: &[OsString]) -> ApplicationStart {
    ApplicationStart::of_process("blastn", words)
}

#[test]
fn detects_database_component_files() {
    let dir = fresh_dir("db-comp");
    let nin = dir.join("mydb.nin");
    let nhr = dir.join("mydb.nhr");
    let naa = dir.join("mydb.naa");
    let shd = dir.join("mydb.shd");
    let idx00 = dir.join("mydb.00.idx");
    let bogus_nma = dir.join("mydb.nma");
    let other = dir.join("mydb.log");
    std::fs::write(&nin, b"index").unwrap();
    std::fs::write(&nhr, b"header").unwrap();
    std::fs::write(&naa, b"mask").unwrap();
    std::fs::write(&shd, b"shd").unwrap();
    std::fs::write(&idx00, b"idx").unwrap();
    std::fs::write(&bogus_nma, b"not-ncbi").unwrap();
    std::fs::write(&other, b"log").unwrap();
    let db = dir.join("mydb");
    let words: Vec<OsString> = vec!["-db".into(), db.as_os_str().to_os_string()];
    let start = blastn_start(&words);
    assert!(
        check_disk_logfile("blastn", &words, &start, Some(nhr.as_os_str())).is_err(),
        "expected blastn refusal for {}",
        nhr.display()
    );
    // The mask column is opened with a mask option only (line WQA D26 round 1): without one it
    // is the run's log; with one it is refused.
    let masked_words: Vec<OsString> = vec![
        "-db".into(),
        db.as_os_str().to_os_string(),
        "-db_soft_mask".into(),
        "0".into(),
    ];
    let masked_start = blastn_start(&masked_words);
    assert!(check_disk_logfile("blastn", &words, &start, Some(naa.as_os_str())).is_ok());
    assert!(check_disk_logfile(
        "blastn",
        &masked_words,
        &masked_start,
        Some(naa.as_os_str())
    )
    .is_err());
    // The MegaBLAST index files are read with `-use_index true` only (line WQA D26).
    let index_words: Vec<OsString> = vec![
        "-db".into(),
        db.as_os_str().to_os_string(),
        "-use_index".into(),
        "true".into(),
    ];
    let index_start = blastn_start(&index_words);
    for refused in [&shd, &idx00] {
        assert!(
            check_disk_logfile(
                "blastn",
                &index_words,
                &index_start,
                Some(refused.as_os_str())
            )
            .is_err(),
            "expected blastn -use_index true refusal for {}",
            refused.display()
        );
    }
    for allowed in [&other, &bogus_nma, &shd, &idx00] {
        assert!(
            check_disk_logfile("blastn", &words, &start, Some(allowed.as_os_str())).is_ok(),
            "expected blastn to allow {}",
            allowed.display()
        );
    }
    // Non-blastn nucleotide search (tblastx) does not read MegaBLAST index files.
    let tbx_start = ApplicationStart::of_process("tblastx", &words);
    assert!(check_disk_logfile("tblastx", &words, &tbx_start, Some(nhr.as_os_str())).is_err());
    assert!(check_disk_logfile("tblastx", &words, &tbx_start, Some(shd.as_os_str())).is_ok());
    assert!(check_disk_logfile("tblastx", &words, &tbx_start, Some(idx00.as_os_str())).is_ok());

    // Protein volume: PIG ISAM (.ppi/.ppd, with an id list), mask column (.paa, with a mask
    // option), RPS (.rps, for the RPS searches), and bogus .pma allowed.
    let pin = dir.join("protdb.pin");
    let ppi = dir.join("protdb.ppi");
    let paa = dir.join("protdb.paa");
    let rps = dir.join("protdb.rps");
    let bogus_pma = dir.join("protdb.pma");
    for f in [&pin, &ppi, &paa, &rps, &bogus_pma] {
        std::fs::write(f, b"x").unwrap();
    }
    let prot_words: Vec<OsString> =
        vec!["-db".into(), dir.join("protdb").as_os_str().to_os_string()];
    let bp_start = ApplicationStart::of_process("blastp", &prot_words);
    for (option, refused, program) in [
        ("-gilist", &ppi, "blastp"),
        ("-db_hard_mask", &paa, "blastp"),
        ("-db", &rps, "rpsblast"),
    ] {
        let mut words = prot_words.clone();
        if option != "-db" {
            words.extend([option.into(), "x".into()]);
        }
        let start = ApplicationStart::of_process(program, &words);
        assert!(
            check_disk_logfile(program, &words, &start, Some(refused.as_os_str())).is_err(),
            "expected {program} {option} refusal for {}",
            refused.display()
        );
        // The same name, not asked for, is the run's log.
        let start = ApplicationStart::of_process("blastp", &prot_words);
        assert!(
            check_disk_logfile("blastp", &prot_words, &start, Some(refused.as_os_str())).is_ok(),
            "expected plain blastp to allow {}",
            refused.display()
        );
    }
    assert!(check_disk_logfile(
        "blastp",
        &prot_words,
        &bp_start,
        Some(bogus_pma.as_os_str())
    )
    .is_ok());

    // Alias pointing via DBLIST to a differently named volume, with GILIST, OIDLIST, and MASKLIST:
    // alias file, target volume components, list files, GI-mask files (an index of one data
    // volume), and with `-use_index true` the old-style `<alias>.00.idx` (when `-index_name` is
    // absent) must all be refused.
    let real_nin = dir.join("real_db.nin");
    let real_nhr = dir.join("real_db.nhr");
    let gi_file = dir.join("alias_gis.txt");
    let oid_file = dir.join("alias_oids.msk");
    let mask_gni = dir.join("mymask.gni");
    let mask_gnd = dir.join("mymask.gnd");
    let alias_idx = dir.join("myalias.00.idx");
    let alias_nal = dir.join("myalias.nal");
    for f in [
        &real_nin, &real_nhr, &gi_file, &oid_file, &mask_gnd, &alias_idx,
    ] {
        std::fs::write(f, b"x").unwrap();
    }
    std::fs::write(&mask_gni, gi_mask_index_header(1)).unwrap();
    std::fs::write(
        &alias_nal,
        format!(
            "TITLE my alias\nDBLIST real_db\nGILIST alias_gis.txt\nOIDLIST alias_oids.msk\nMASKLIST {}\n",
            dir.join("mymask").display()
        ),
    )
    .unwrap();
    let alias_db = dir.join("myalias");
    let alias_words: Vec<OsString> = vec![
        "-db".into(),
        alias_db.as_os_str().to_os_string(),
        "-use_index".into(),
        "true".into(),
    ];
    let alias_start = blastn_start(&alias_words);
    // The GI-mask files are opened with a mask option only (line WQA D26 round 1).
    for allowed in [&mask_gni, &mask_gnd] {
        assert!(
            check_disk_logfile(
                "blastn",
                &alias_words,
                &alias_start,
                Some(allowed.as_os_str())
            )
            .is_ok(),
            "expected alias blastn to allow {} without a mask option",
            allowed.display()
        );
    }
    let mut mask_words = alias_words.clone();
    mask_words.extend(["-db_soft_mask".into(), "0".into()]);
    let mask_start = blastn_start(&mask_words);
    for refused in [&mask_gni, &mask_gnd] {
        assert!(
            check_disk_logfile(
                "blastn",
                &mask_words,
                &mask_start,
                Some(refused.as_os_str())
            )
            .is_err(),
            "expected alias blastn -db_soft_mask refusal for {}",
            refused.display()
        );
    }
    for refused in [&real_nhr, &alias_nal, &gi_file, &oid_file, &alias_idx] {
        assert!(
            check_disk_logfile(
                "blastn",
                &alias_words,
                &alias_start,
                Some(refused.as_os_str())
            )
            .is_err(),
            "expected alias blastn refusal for {}",
            refused.display()
        );
    }
    assert!(check_disk_logfile(
        "blastn",
        &alias_words,
        &alias_start,
        Some(other.as_os_str())
    )
    .is_ok());

    // Registry [BLAST] BLASTDB resolves bare database names and aliases.
    let conf = dir.join("blastdb.ini");
    std::fs::write(&conf, format!("[BLAST]\nBLASTDB={}\n", dir.display())).unwrap();
    let reg_words: Vec<OsString> = vec![
        "-conffile".into(),
        conf.as_os_str().to_os_string(),
        "-db".into(),
        "myalias".into(),
    ];
    let reg_start = blastn_start(&reg_words);
    assert!(
        check_disk_logfile("blastn", &reg_words, &reg_start, Some(real_nhr.as_os_str())).is_err()
    );
    assert!(check_disk_logfile(
        "blastn",
        &reg_words,
        &reg_start,
        Some(alias_nal.as_os_str())
    )
    .is_err());

    // A database that does not resolve: only the files its resolution checks are refused (line
    // WQA D26: `-logfile` would create them and the run then read them), not its volume files.
    let missing_nhr = dir.join("no_such_db.nhr");
    std::fs::write(&missing_nhr, b"stale").unwrap();
    let missing_words: Vec<OsString> = vec![
        "-db".into(),
        dir.join("no_such_db").as_os_str().to_os_string(),
    ];
    let missing_start = blastn_start(&missing_words);
    assert!(check_disk_logfile(
        "blastn",
        &missing_words,
        &missing_start,
        Some(missing_nhr.as_os_str())
    )
    .is_ok());
    for probed in ["no_such_db.nal", "no_such_db.nin", "index.alx"] {
        assert!(
            check_disk_logfile(
                "blastn",
                &missing_words,
                &missing_start,
                Some(dir.join(probed).as_os_str())
            )
            .is_err(),
            "expected refusal of the probed {probed}"
        );
    }

    let _ = std::fs::remove_dir_all(&dir);
}

/// The first 8 bytes of a version-1 GI-mask index naming `num_vols` data volumes.
fn gi_mask_index_header(num_vols: i32) -> Vec<u8> {
    [1i32.to_be_bytes(), num_vols.to_be_bytes()].concat()
}

fn refused(program: &str, words: &[OsString], log: &std::path::Path) -> bool {
    let start = ApplicationStart::of_process(program, words);
    check_disk_logfile(program, words, &start, Some(log.as_os_str())).is_err()
}

fn db_words(db: &std::path::Path, extra: &[&str]) -> Vec<OsString> {
    let mut words: Vec<OsString> = vec!["-db".into(), db.as_os_str().to_os_string()];
    words.extend(extra.iter().map(OsString::from));
    words
}

// Line WQA D26 item 1: the `OIDLIST` path beside the first volume is read only where the path
// beside the alias file is not a file (`compute_masks`), so it is refused only then.
#[test]
fn oidlist_fallback_is_refused_only_where_the_run_reads_it() {
    let dir = fresh_dir("d26-oidlist");
    let alias_dir = dir.join("al");
    let vol_dir = dir.join("vols");
    std::fs::create_dir_all(&alias_dir).unwrap();
    std::fs::create_dir_all(&vol_dir).unwrap();
    for ext in ["nin", "nhr", "nsq"] {
        std::fs::write(vol_dir.join(format!("v1.{ext}")), b"x").unwrap();
    }
    std::fs::write(
        alias_dir.join("m.nal"),
        b"TITLE m\nDBLIST v1\nOIDLIST o.msk\n",
    )
    .unwrap();
    let conf = dir.join("blastdb.ini");
    std::fs::write(&conf, format!("[BLAST]\nBLASTDB={}\n", vol_dir.display())).unwrap();
    let words = db_words(&alias_dir.join("m"), &["-conffile", conf.to_str().unwrap()]);
    let primary = alias_dir.join("o.msk");
    let fallback = vol_dir.join("o.msk");
    // The primary path is missing: the run reads the fallback, and `-logfile` naming the primary
    // would create it and make the run read it instead.
    assert!(refused("blastn", &words, &primary));
    assert!(refused("blastn", &words, &fallback));
    // The primary path is a file: the run never looks beside the volume.
    std::fs::write(&primary, b"x").unwrap();
    assert!(refused("blastn", &words, &primary));
    assert!(
        !refused("blastn", &words, &fallback),
        "the fallback beside the volume is not read when the primary exists"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

// Line WQA D26 item 2: the lookup candidates that the resolution checks before the one found
// (the group alias file, the alias file beside a volume found by its index file) are refused
// as LOSAT Web refuses them, since `-logfile` would create them and the run then read them.
#[test]
fn lookup_candidates_of_a_resolved_database_are_refused() {
    let dir = fresh_dir("d26-lookup");
    for ext in ["nin", "nhr", "nsq"] {
        std::fs::write(dir.join(format!("x.{ext}")), b"x").unwrap();
    }
    let words = db_words(&dir.join("x"), &[]);
    for name in ["index.alx", "x.nal", "x.nsq"] {
        assert!(refused("blastn", &words, &dir.join(name)), "{name}");
    }
    // The GI index is opened with a mask option only (round 1).
    assert!(!refused("blastn", &words, &dir.join("x.nog")));
    let masked = db_words(&dir.join("x"), &["-db_hard_mask", "0"]);
    assert!(refused("blastn", &masked, &dir.join("x.nog")));
    assert!(!refused("blastn", &words, &dir.join("x.log")));
    let _ = std::fs::remove_dir_all(&dir);
}

// Line WQA D26 item 3: `-index_name base,,first,last` reads the volumes of the whole range.
#[test]
fn old_style_index_range_is_refused_whole() {
    let dir = fresh_dir("d26-index-range");
    for ext in ["nin", "nhr", "nsq"] {
        std::fs::write(dir.join(format!("db.{ext}")), b"x").unwrap();
    }
    let db = dir.join("db");
    let range = format!("{},,0,150", db.display());
    let words = db_words(&db, &["-use_index", "true", "-index_name", &range]);
    for name in ["db.00.idx", "db.99.idx", "db.120.idx", "db.150.idx"] {
        assert!(refused("blastn", &words, &dir.join(name)), "{name}");
    }
    assert!(!refused("blastn", &words, &dir.join("db.151.idx")));
    // With `-index_name` the loader looks for no new-style index of the `-db` volumes.
    assert!(!refused("blastn", &words, &dir.join("db.shd")));
    // A first volume past the last reads no volume; a range NCBI's loop never finishes stops the
    // loader before it reads one.
    for (range, name) in [(",,5,2", "db.03.idx"), (",,0,-1", "db.00.idx")] {
        let index = format!("{}{range}", db.display());
        let words = db_words(&db, &["-use_index", "true", "-index_name", &index]);
        assert!(!refused("blastn", &words, &dir.join(name)), "{range}");
    }
    // Without `-use_index true`, with another task, or with a template, no index is read.
    for extra in [
        &["-index_name", &range][..],
        &["-use_index", "false", "-index_name", &range],
        &[
            "-use_index",
            "true",
            "-task",
            "blastn",
            "-index_name",
            &range,
        ],
        &[
            "-use_index",
            "true",
            "-index_name",
            &range,
            "-template_type",
            "coding",
            "-template_length",
            "16",
        ],
    ] {
        let words = db_words(&db, extra);
        assert!(
            !refused("blastn", &words, &dir.join("db.120.idx")),
            "{extra:?}"
        );
        assert!(
            !refused("blastn", &words, &dir.join("db.00.idx")),
            "{extra:?}"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

// Line WQA D26 item 4: the new-style index volumes of a database volume are read only where its
// super header loads, up to the first one not found; the old-style names of the `-db` value stay.
#[test]
fn new_style_index_volumes_are_refused_only_where_the_run_reads_them() {
    use crate::dbindex::format::SuperHeader;
    let dir = fresh_dir("d26-new-style");
    for volume in ["v1", "v2"] {
        for ext in ["nin", "nhr", "nsq"] {
            std::fs::write(dir.join(format!("{volume}.{ext}")), b"x").unwrap();
        }
    }
    std::fs::write(dir.join("al.nal"), b"TITLE al\nDBLIST v1 v2\n").unwrap();
    let shd = dir.join("v1.shd");
    assert!(!SuperHeader {
        num_seq: 1,
        num_vol: 3
    }
    .save(shd.as_os_str().as_encoded_bytes()));
    std::fs::write(dir.join("v1.00.idx"), b"x").unwrap();
    let words = db_words(&dir.join("al"), &["-use_index", "true"]);
    for name in [
        "v1.shd",
        "v1.00.idx",
        "v1.01.idx",
        "v2.shd",
        "al.00.idx",
        "al.99.idx",
    ] {
        assert!(refused("blastn", &words, &dir.join(name)), "{name}");
    }
    // The loader stops at the first index volume not found (v1.01.idx), and a volume without a
    // super header (v2) has no index volume that it reads.
    for name in ["v1.02.idx", "v1.03.idx", "v2.00.idx", "v2.99.idx"] {
        assert!(!refused("blastn", &words, &dir.join(name)), "{name}");
    }
    // `-filtering_db` is no index database.
    let filtering: Vec<OsString> = vec![
        "-subject".into(),
        dir.join("s.fa").into_os_string(),
        "-filtering_db".into(),
        dir.join("v1").into_os_string(),
        "-use_index".into(),
        "true".into(),
    ];
    assert!(refused("blastn", &filtering, &dir.join("v1.nsq")));
    assert!(!refused("blastn", &filtering, &dir.join("v1.shd")));
    let _ = std::fs::remove_dir_all(&dir);
}

// Line WQA D26 item 5: the GI-mask data files are read only after an index that reads, up to the
// first one not found; the GI-mask files are opened only with a mask option (round 1).
#[test]
fn gi_mask_data_files_are_refused_only_where_the_run_reads_them() {
    let dir = fresh_dir("d26-gi-mask");
    for ext in ["nin", "nhr", "nsq"] {
        std::fs::write(dir.join(format!("v.{ext}")), b"x").unwrap();
    }
    let mask = dir.join("mk");
    std::fs::write(
        dir.join("m.nal"),
        format!("TITLE m\nDBLIST v\nMASKLIST {}\n", mask.display()),
    )
    .unwrap();
    let plain = db_words(&dir.join("m"), &[]);
    for name in ["mk.gni", "mk.gno", "mk.gnd"] {
        assert!(!refused("blastn", &plain, &dir.join(name)), "{name}");
    }
    let words = db_words(&dir.join("m"), &["-db_soft_mask", "0"]);
    // No index: its name and the offset file's are refused, no data file.
    for name in ["mk.gni", "mk.gno"] {
        assert!(refused("blastn", &words, &dir.join(name)), "{name}");
    }
    assert!(!refused("blastn", &words, &dir.join("mk.gnd")));
    // An index of one data volume reads `mk.gnd`.
    std::fs::write(dir.join("mk.gni"), gi_mask_index_header(1)).unwrap();
    assert!(refused("blastn", &words, &dir.join("mk.gnd")));
    assert!(!refused("blastn", &words, &dir.join("mk.00.gnd")));
    // An index of three: `mk.00.gnd`, then the first missing one, `mk.01.gnd`.
    std::fs::write(dir.join("mk.gni"), gi_mask_index_header(3)).unwrap();
    std::fs::write(dir.join("mk.00.gnd"), b"x").unwrap();
    for name in ["mk.00.gnd", "mk.01.gnd"] {
        assert!(refused("blastn", &words, &dir.join(name)), "{name}");
    }
    for name in ["mk.02.gnd", "mk.gnd"] {
        assert!(!refused("blastn", &words, &dir.join(name)), "{name}");
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn version_request_does_not_collect_argv_inputs() {
    let dir = fresh_dir("version-run");
    let query = dir.join("ver_q.fa");
    std::fs::write(&query, b">vq\nACGT\n").unwrap();
    let words: Vec<OsString> = vec![
        "-version".into(),
        "-query".into(),
        query.as_os_str().to_os_string(),
        "-logfile".into(),
        query.as_os_str().to_os_string(),
    ];
    let start = blastn_start(&words);
    assert!(check_disk_logfile("blastn", &words, &start, Some(query.as_os_str())).is_ok());
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn missing_parent_directory_is_not_refused() {
    let dir = fresh_dir("missing-parent");
    let missing_dir_file = dir.join("no_such_dir").join("q.fa");
    let words: Vec<OsString> = vec!["-query".into(), missing_dir_file.as_os_str().to_os_string()];
    let start = blastn_start(&words);
    // open_log_file will not create a file in a missing parent dir; it falls back to stderr.
    assert!(
        check_disk_logfile("blastn", &words, &start, Some(missing_dir_file.as_os_str())).is_ok()
    );
    let _ = std::fs::remove_dir_all(&dir);
}
