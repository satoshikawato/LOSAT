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
    let nhr = dir.join("mydb.nhr");
    let other = dir.join("mydb.log");
    std::fs::write(&nhr, b"header").unwrap();
    std::fs::write(&other, b"log").unwrap();
    let db = dir.join("mydb");
    let words: Vec<OsString> = vec!["-db".into(), db.as_os_str().to_os_string()];
    let start = blastn_start(&words);
    assert!(check_disk_logfile("blastn", &words, &start, Some(nhr.as_os_str())).is_err());
    assert!(check_disk_logfile("blastn", &words, &start, Some(other.as_os_str())).is_ok());
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
