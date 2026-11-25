use std::ffi::OsString;
use std::fs;

use super::*;
use crate::cleanup_blastdb_volumes::run_with;
use crate::utils::test_dir::TempDir;

fn invoke(words: &[&str], env: &Env) -> (i32, Vec<u8>, Vec<u8>) {
    let args: Vec<OsString> = words.iter().map(OsString::from).collect();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let rc = run_with(&args, env, &mut out, &mut err);
    (rc, out, err)
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
        cwd: Some(path_bytes(&cwd).into_owned()),
        home: Some(path_bytes(&home).into_owned()),
        ncbi: Some(path_bytes(&ncbi).into_owned()),
        blastdb: None,
    };
    assert_eq!(
        find_blastdb(b"idb", false, &env_cfg),
        Ok(Some(path_join(path_bytes(&ini_target).as_ref(), b"idb")))
    );
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

    // Dry-run without verbose.
    let (rc, out, err) = invoke(&["-db", db_str, "-dbtype", "prot", "-dry-run"], &env);
    assert_eq!(rc, 0);
    assert!(err.is_empty());
    assert_eq!(
        String::from_utf8(out).unwrap(),
        format!("Will remove extra volume {db_str}.01\n")
    );
    assert!(dir.join("mydb.01.pin").exists());

    // Dry-run with verbose: sorted `.p??` then `.tar.gz.md5`.
    let (rc, out, err) = invoke(
        &["-db", db_str, "-dbtype", "prot", "-dry-run", "-verbose"],
        &env,
    );
    assert_eq!(rc, 0);
    assert!(err.is_empty());
    assert_eq!(
        String::from_utf8(out).unwrap(),
        format!(
            "Will remove extra volume {db_str}.01\n\
             Will remove {db_str}.01.phr\n\
             Will remove {db_str}.01.pin\n\
             Will remove {db_str}.01.psq\n\
             Will remove {db_str}.01.tar.gz.md5\n"
        )
    );

    // Actual deletion.
    let (rc, out, err) = invoke(&["-db", db_str, "-dbtype", "prot"], &env);
    assert_eq!(rc, 0);
    assert!(err.is_empty());
    assert_eq!(
        String::from_utf8(out).unwrap(),
        format!(
            "Removed {db_str}.01.phr\n\
             Removed {db_str}.01.pin\n\
             Removed {db_str}.01.psq\n\
             Removed {db_str}.01.tar.gz.md5\n"
        )
    );
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
