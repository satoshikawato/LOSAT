use std::ffi::OsString;
use std::fs;
use std::io;

use super::*;
use crate::cleanup_blastdb_volumes::{run_with, run_with_buffering};
use crate::utils::test_dir::TempDir;

fn invoke(words: &[&str], env: &Env) -> (i32, Vec<u8>, Vec<u8>) {
    let args: Vec<OsString> = words.iter().map(OsString::from).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let rc = run_with(&args, env, &mut out, &mut err);
    (rc, out, err)
}

struct FailingWriter {
    error_kind: io::ErrorKind,
    os_code: Option<i32>,
    writes: usize,
}

impl io::Write for FailingWriter {
    fn write(&mut self, _buf: &[u8]) -> io::Result<usize> {
        self.writes += 1;
        if let Some(code) = self.os_code {
            Err(io::Error::from_raw_os_error(code))
        } else {
            Err(io::Error::new(self.error_kind, "simulated write failure"))
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/cleanup-blastdb-volumes.py:61-65
#[test]
fn find_blastdb_direct_and_env_precedence() {
    let dir = TempDir::unique("cleanup-find");
    let direct = dir.join("direct");
    let blastdb = dir.join("blastdb");
    let cwd = dir.join("cwd");
    let home = dir.join("home");
    let ncbi = dir.join("ncbi");
    let ini_target = dir.join("ini_target");
    for d in [&direct, &blastdb, &cwd, &home, &ncbi, &ini_target] {
        fs::create_dir_all(d).unwrap();
    }

    let p = direct.join("only_pin");
    fs::write(direct.join("only_pin.pin"), b"").unwrap();
    let empty_env = Env::default();
    assert_eq!(
        find_blastdb(path_bytes(&p).as_ref(), true, &empty_env),
        Ok(Some(path_bytes(&p).into_owned()))
    );
    assert_eq!(
        find_blastdb(path_bytes(&p).as_ref(), false, &empty_env),
        Ok(None)
    );

    // BLASTDB env lookup.
    fs::write(blastdb.join("envdb.pal"), b"DBLIST envdb.00\n").unwrap();
    let env_blastdb = Env {
        blastdb: Some(path_bytes(&blastdb).into_owned()),
        ..Env::default()
    };
    assert_eq!(
        find_blastdb(b"envdb", true, &env_blastdb),
        Ok(Some(path_join(path_bytes(&blastdb).as_ref(), b"envdb")))
    );

    // Config fallback: cwd .ncbirc -> cwd ncbi.ini -> home .ncbirc -> ncbi .ncbirc.
    fs::write(ini_target.join("idb.nal"), b"DBLIST idb.00\n").unwrap();
    fs::write(cwd.join(".ncbirc"), b"[BLAST]\nBLASTDB=/nonexistent_dir\n").unwrap();
    fs::write(
        home.join(".ncbirc"),
        format!("[BLAST]\nBLASTDB={}\n", ini_target.display()),
    )
    .unwrap();
    let env_cfg = Env {
        cwd: Some(Ok(path_bytes(&cwd).into_owned())),
        home: Some(path_bytes(&home).into_owned()),
        ncbi: Some(path_bytes(&ncbi).into_owned()),
        ..Env::default()
    };
    assert_eq!(
        find_blastdb(b"idb", false, &env_cfg),
        Ok(Some(path_join(path_bytes(&ini_target).as_ref(), b"idb")))
    );

    // Removed cwd fails only when lookup falls through to config directories.
    let env_bad_cwd = Env {
        cwd: Some(Err(
            "cannot determine current working directory: No such file or directory".to_string(),
        )),
        blastdb: Some(path_bytes(&blastdb).into_owned()),
        ..Env::default()
    };
    assert_eq!(
        find_blastdb(b"envdb", true, &env_bad_cwd),
        Ok(Some(path_join(path_bytes(&blastdb).as_ref(), b"envdb")))
    );
    assert!(find_blastdb(b"missing_db", true, &env_bad_cwd).is_err());
}

#[test]
fn dry_run_verbose_and_deletion_preserve_referenced_and_unrelated_files() {
    let dir = TempDir::unique("cleanup-run");
    let db_path = dir.join("mydb");
    fs::write(dir.join("mydb.pal"), b"TITLE test\nDBLIST \"mydb.00\"\n").unwrap();
    for name in [
        "mydb.00.pin",
        "mydb.00.psq",
        "mydb.01.pin",
        "mydb.01.phr",
        "mydb.01.psq",
        "mydb.01.tar.gz.md5",
        "mydb.01.pa",
        "mydb.01.pinx",
        "mydb.01.nin",
        "other.01.pin",
    ] {
        fs::write(dir.join(name), b"").unwrap();
    }
    let db_str = db_path.to_str().unwrap();
    let env = Env::default();

    // Expected per-volume order matches `read_dir` order of `.p??`, followed by `.tar.gz.md5`.
    let mut expected_files: Vec<String> = fs::read_dir(&dir)
        .unwrap()
        .flatten()
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|name| matches!(name.as_str(), "mydb.01.phr" | "mydb.01.pin" | "mydb.01.psq"))
        .map(|name| dir.join(name).to_string_lossy().into_owned())
        .collect();
    expected_files.push(
        dir.join("mydb.01.tar.gz.md5")
            .to_string_lossy()
            .into_owned(),
    );

    // Dry-run without verbose.
    let (rc, out, err) = invoke(&["-db", db_str, "-dbtype", "prot", "-dry-run"], &env);
    assert_eq!(rc, 0);
    assert!(err.is_empty());
    assert_eq!(
        String::from_utf8(out).unwrap(),
        format!("Will remove extra volume {db_str}.01\n")
    );
    assert!(dir.join("mydb.01.pin").exists());

    // Dry-run with verbose: `read_dir` order of `.p??` then `.tar.gz.md5`.
    let (rc, out, err) = invoke(
        &["-db", db_str, "-dbtype", "prot", "-dry-run", "-verbose"],
        &env,
    );
    assert_eq!(rc, 0);
    assert!(err.is_empty());
    let mut want_dry = format!("Will remove extra volume {db_str}.01\n");
    for f in &expected_files {
        want_dry.push_str(&format!("Will remove {f}\n"));
    }
    assert_eq!(String::from_utf8(out).unwrap(), want_dry);

    // Actual deletion.
    let (rc, out, err) = invoke(&["-db", db_str, "-dbtype", "prot"], &env);
    assert_eq!(rc, 0);
    assert!(err.is_empty());
    let mut want_rm = String::new();
    for f in &expected_files {
        want_rm.push_str(&format!("Removed {f}\n"));
    }
    assert_eq!(String::from_utf8(out).unwrap(), want_rm);
    for removed in [
        "mydb.01.phr",
        "mydb.01.pin",
        "mydb.01.psq",
        "mydb.01.tar.gz.md5",
    ] {
        assert!(!dir.join(removed).exists(), "{removed} should be removed");
    }
    for kept in [
        "mydb.pal",
        "mydb.00.pin",
        "mydb.00.psq",
        "mydb.01.pa",
        "mydb.01.pinx",
        "mydb.01.nin",
        "other.01.pin",
    ] {
        assert!(dir.join(kept).exists(), "{kept} should be kept");
    }
}

#[test]
fn dblist_splits_on_python_whitespace_u001c_to_u001f() {
    let dir = TempDir::unique("cleanup-pyws");
    let db_path = dir.join("db");
    fs::write(
        dir.join("db.pal"),
        "DBLIST db.01\u{001c}db.02\u{001d}db.03\n",
    )
    .unwrap();
    for name in ["db.01.pin", "db.02.pin", "db.03.pin", "db.04.pin"] {
        fs::write(dir.join(name), b"").unwrap();
    }
    let db_str = db_path.to_str().unwrap();
    let (rc, out, err) = invoke(&["-db", db_str, "-dbtype", "prot"], &Env::default());
    assert_eq!(rc, 0);
    assert!(err.is_empty());
    assert_eq!(
        String::from_utf8(out).unwrap(),
        format!("Removed {db_str}.04.pin\n")
    );
    assert!(dir.join("db.01.pin").exists());
    assert!(dir.join("db.02.pin").exists());
    assert!(dir.join("db.03.pin").exists());
    assert!(!dir.join("db.04.pin").exists());
}

#[test]
fn alias_8kib_chunk_decoding_processes_earlier_chunks_before_invalid_utf8() {
    // Bad byte inside chunk 1 (byte 8191): nothing removed.
    let dir1 = TempDir::unique("cleanup-chunk-8191");
    let mut alias1 = b"DBLIST db.00\n".to_vec();
    alias1.resize(8191, b'#');
    alias1.push(0xff);
    fs::write(dir1.join("db.pal"), &alias1).unwrap();
    fs::write(dir1.join("db.00.pin"), b"").unwrap();
    fs::write(dir1.join("db.01.pin"), b"").unwrap();
    let db1 = dir1.join("db").to_str().unwrap().to_string();
    let (rc, out, err) = invoke(&["-db", &db1, "-dbtype", "prot"], &Env::default());
    assert_eq!(rc, 1);
    assert!(out.is_empty());
    assert!(String::from_utf8_lossy(&err).contains("invalid UTF-8"));
    assert!(dir1.join("db.01.pin").exists());

    // Bad byte at start of chunk 2 (byte 8192): chunk 1's DBLIST is processed first.
    let dir2 = TempDir::unique("cleanup-chunk-8192");
    let mut alias2 = b"DBLIST db.00\n".to_vec();
    alias2.resize(8191, b'#');
    alias2.push(b'\n');
    alias2.push(0xff);
    fs::write(dir2.join("db.pal"), &alias2).unwrap();
    fs::write(dir2.join("db.00.pin"), b"").unwrap();
    fs::write(dir2.join("db.01.pin"), b"").unwrap();
    let db2 = dir2.join("db").to_str().unwrap().to_string();
    let (rc, out, err) = invoke(&["-db", &db2, "-dbtype", "prot"], &Env::default());
    assert_eq!(rc, 1);
    assert_eq!(
        String::from_utf8(out).unwrap(),
        format!("Removed {db2}.01.pin\n")
    );
    assert!(String::from_utf8_lossy(&err).contains("invalid UTF-8"));
    assert!(!dir2.join("db.01.pin").exists());

    // Valid multibyte UTF-8 split across byte 8191..8194 succeeds.
    let dir3 = TempDir::unique("cleanup-chunk-split-utf8");
    let mut alias3 = b"DBLIST db.00\n".to_vec();
    alias3.resize(8191, b'#');
    alias3.extend_from_slice("\u{3042}\n".as_bytes());
    fs::write(dir3.join("db.pal"), &alias3).unwrap();
    fs::write(dir3.join("db.00.pin"), b"").unwrap();
    fs::write(dir3.join("db.01.pin"), b"").unwrap();
    let db3 = dir3.join("db").to_str().unwrap().to_string();
    let (rc, out, err) = invoke(&["-db", &db3, "-dbtype", "prot"], &Env::default());
    assert_eq!(rc, 0);
    assert!(err.is_empty());
    assert_eq!(
        String::from_utf8(out).unwrap(),
        format!("Removed {db3}.01.pin\n")
    );
    assert!(!dir3.join("db.01.pin").exists());

    // DBLIST ending in `\r` at byte 8191 peeks chunk 2 before yielding -> fails before deleting.
    let dir4 = TempDir::unique("cleanup-chunk-cr-8191");
    let dl = b"DBLIST db.00\r";
    let mut alias4 = vec![b'#'; 8192 - dl.len() - 1];
    alias4.push(b'\n');
    alias4.extend_from_slice(dl);
    assert_eq!(alias4.len(), 8192);
    alias4.push(0xff);
    fs::write(dir4.join("db.pal"), &alias4).unwrap();
    fs::write(dir4.join("db.00.pin"), b"").unwrap();
    fs::write(dir4.join("db.01.pin"), b"").unwrap();
    let db4 = dir4.join("db").to_str().unwrap().to_string();
    let (rc, out, err) = invoke(&["-db", &db4, "-dbtype", "prot"], &Env::default());
    assert_eq!(rc, 1);
    assert!(out.is_empty());
    assert!(String::from_utf8_lossy(&err).contains("invalid UTF-8"));
    assert!(dir4.join("db.01.pin").exists());
}

#[test]
fn stdout_buffering_matches_python_block_and_line_failure_modes() {
    // Block-buffered (< 8 KiB): removes all files, then fails at exit flush with exit 120.
    let dir = TempDir::unique("cleanup-buf-small");
    fs::write(dir.join("db.pal"), b"DBLIST db.01\n").unwrap();
    for name in [
        "db.01.pin",
        "db.02.pin",
        "db.02.psq",
        "db.03.pin",
        "db.03.psq",
    ] {
        fs::write(dir.join(name), b"").unwrap();
    }
    let db_str = dir.join("db").to_str().unwrap().to_string();
    let args: Vec<OsString> = ["-db", &db_str, "-dbtype", "prot"]
        .iter()
        .map(OsString::from)
        .collect();
    let mut out = FailingWriter {
        error_kind: io::ErrorKind::BrokenPipe,
        os_code: None,
        writes: 0,
    };
    let mut err = Vec::new();
    let rc = run_with_buffering(&args, &Env::default(), &mut out, false, &mut err);
    assert_eq!(rc, 120);
    assert_eq!(out.writes, 1);
    let err_text = String::from_utf8(err).unwrap();
    assert!(err_text.starts_with("Exception ignored on flushing sys.stdout:\nBrokenPipeError:"));
    for removed in ["db.02.pin", "db.02.psq", "db.03.pin", "db.03.psq"] {
        assert!(!dir.join(removed).exists(), "{removed} should be removed");
    }

    // Line-buffered (terminal): fails on the very first line (after unlinking 1 file), exit 1.
    let dir_tty = TempDir::unique("cleanup-buf-tty");
    fs::write(dir_tty.join("db.pal"), b"DBLIST db.01\n").unwrap();
    for name in ["db.01.pin", "db.02.pin", "db.03.pin"] {
        fs::write(dir_tty.join(name), b"").unwrap();
    }
    let db_tty = dir_tty.join("db").to_str().unwrap().to_string();
    let args_tty: Vec<OsString> = ["-db", &db_tty, "-dbtype", "prot"]
        .iter()
        .map(OsString::from)
        .collect();
    let mut out_tty = FailingWriter {
        error_kind: io::ErrorKind::BrokenPipe,
        os_code: None,
        writes: 0,
    };
    let mut err_tty = Vec::new();
    let rc_tty = run_with_buffering(&args_tty, &Env::default(), &mut out_tty, true, &mut err_tty);
    assert_eq!(rc_tty, 1);
    assert_eq!(out_tty.writes, 1);
    assert!(!dir_tty.join("db.02.pin").exists());
    assert!(dir_tty.join("db.03.pin").exists());
}

#[test]
fn single_volume_without_alias_exits_one_silently() {
    let dir = TempDir::unique("cleanup-single");
    let db_path = dir.join("sv");
    fs::write(dir.join("sv.pin"), b"").unwrap();
    fs::write(dir.join("sv.01.pin"), b"").unwrap();
    let (rc, out, err) = invoke(
        &["-db", db_path.to_str().unwrap(), "-dbtype", "prot"],
        &Env::default(),
    );
    assert_eq!(rc, 1);
    assert!(out.is_empty());
    assert!(err.is_empty());
    assert!(dir.join("sv.01.pin").exists());
}
