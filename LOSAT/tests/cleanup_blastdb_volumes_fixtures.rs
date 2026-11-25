//! `losat cleanup-blastdb-volumes.py` against NCBI BLAST+ 2.17.0 `cleanup-blastdb-volumes.py`
//! (run under Python 3.13): each row of `tests/fixtures/cleanup_blastdb_volumes/cases.tsv` runs in
//! a fresh copy of `in/` and checks stdout, stderr, exit status, and exact file deletions against
//! the frozen official NCBI oracle (`ncbi/<row>.{out,err,rc,files}`). NCBI is not run here.
//!
//! Within each extra volume, both NCBI's Python `glob.glob` (`os.scandir`) and LOSAT visit
//! `.p??`/`.n??` files in filesystem `read_dir` order before `.tar.gz.md5`. Because each test run
//! copies `in/` into a fresh directory whose `read_dir` order can differ from the frozen machine's,
//! the comparison checks both that each volume's lines equal NCBI's frozen lines as a multiset (with
//! `.tar.gz.md5` last) and that their order strictly matches the test copy's `fs::read_dir` order,
//! while volume order, non-verbose `-dry-run`, `-h`, `COLUMNS`, stderr, exit code, and removed/kept
//! files match byte for byte.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/cleanup_blastdb_volumes")
}

struct Case {
    name: String,
    env: Vec<(String, String)>,
    cwd: String,
    args: Vec<String>,
}

fn cases() -> Vec<Case> {
    let table = fs::read_to_string(fixtures().join("cases.tsv")).unwrap();
    table
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            assert_eq!(fields.len(), 4, "{line}");
            let env = if fields[1] == "-" {
                Vec::new()
            } else {
                fields[1]
                    .split(',')
                    .map(|pair| {
                        let (k, v) = pair.split_once('=').unwrap();
                        (k.to_string(), v.to_string())
                    })
                    .collect()
            };
            let args = if fields[3].is_empty() {
                Vec::new()
            } else {
                fields[3].split(' ').map(str::to_string).collect()
            };
            Case {
                name: fields[0].to_string(),
                env,
                cwd: fields[2].to_string(),
                args,
            }
        })
        .collect()
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, &byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

fn copy_tree(from: &Path, to: &Path) {
    fs::create_dir_all(to).unwrap();
    for entry in fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_tree(&entry.path(), &target);
        } else {
            fs::copy(entry.path(), target).unwrap();
        }
    }
}

fn rel_posix(path: &Path, root: &Path) -> String {
    path.strip_prefix(root)
        .unwrap()
        .components()
        .map(|part| part.as_os_str().to_string_lossy().into_owned())
        .collect::<Vec<_>>()
        .join("/")
}

fn listing(dir: &Path) -> BTreeMap<String, Vec<u8>> {
    let mut found = BTreeMap::new();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(at) = stack.pop() {
        for entry in fs::read_dir(&at).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                stack.push(path);
            } else {
                found.insert(rel_posix(&path, dir), fs::read(&path).unwrap());
            }
        }
    }
    found
}

/// Records `fs::read_dir` entry index per `(directory_rel, file_name)` for every directory in
/// `root`, so the test can verify that per-volume file lines follow the scratch copy's `read_dir`
/// order.
fn readdir_positions(root: &Path) -> BTreeMap<String, BTreeMap<String, usize>> {
    let mut out = BTreeMap::new();
    let mut stack = vec![root.to_path_buf()];
    while let Some(at) = stack.pop() {
        let rel_dir = rel_posix(&at, root);
        let mut order = BTreeMap::new();
        for (idx, entry) in fs::read_dir(&at).unwrap().enumerate() {
            let entry = entry.unwrap();
            let path = entry.path();
            if path.is_dir() {
                stack.push(path);
            } else {
                order.insert(entry.file_name().to_string_lossy().into_owned(), idx);
            }
        }
        out.insert(rel_dir, order);
    }
    out
}

fn normalize_run_dir(dir: &Path) -> String {
    let text = dir.to_string_lossy();
    #[cfg(windows)]
    let text = text.strip_prefix(r"\\?\").unwrap_or(&text);
    text.to_string()
}

fn replace_run_dir(bytes: &[u8], run_str: &str) -> Vec<u8> {
    let text = String::from_utf8_lossy(bytes);
    #[cfg(windows)]
    let replaced = {
        let verbatim = format!(r"\\?\{run_str}");
        let slash = run_str.replace('\\', "/");
        text.replace(&verbatim, "<RUN>")
            .replace(run_str, "<RUN>")
            .replace(&slash, "<RUN>")
            .replace('\\', "/")
    };
    #[cfg(not(windows))]
    let replaced = text.replace(run_str, "<RUN>");
    replaced.into_bytes()
}

fn scratch(name: &str) -> PathBuf {
    Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!(
        "cleanup-blastdb-volumes-{name}-{}",
        std::process::id()
    ))
}

fn base_command(cwd: &Path, empty_home: &Path) -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || text.starts_with("NCBI")
            || text.starts_with("BLAST")
            || text == "HOME"
            || text == "USERPROFILE"
            || text == "HOMEDRIVE"
            || text == "HOMEPATH"
            || text == "COLUMNS"
        {
            cmd.env_remove(&name);
        }
    }
    if let Some(log) = std::env::var_os("LOSAT_BLASTX_WORKER_LOG") {
        cmd.env("LOSAT_BLASTX_WORKER_LOG", log);
    }
    cmd.env("HOME", empty_home)
        .env("USERPROFILE", empty_home)
        .env("TZ", "UTC")
        .current_dir(cwd)
        .arg("cleanup-blastdb-volumes.py")
        .stdin(Stdio::null());
    cmd
}

/// Resolves a normalized stdout file path (`<RUN>/in/envdb/edb.01.pin` or `in/p5.09.pin`) into
/// `(rel_dir_under_run, file_name)` to look up its `fs::read_dir` position.
fn resolve_rel_dir_and_name(raw_path: &str, cwd: &str) -> (String, String) {
    let stripped = raw_path
        .strip_prefix("<RUN>/")
        .or_else(|| raw_path.strip_prefix("./"))
        .unwrap_or(raw_path);
    let combined = if raw_path.starts_with("<RUN>/") || cwd == "." {
        stripped.to_string()
    } else {
        format!("{cwd}/{stripped}")
    };
    let mut parts: Vec<&str> = combined
        .split('/')
        .filter(|p| !p.is_empty() && *p != ".")
        .collect();
    let file = parts.pop().unwrap_or_default().to_string();
    (parts.join("/"), file)
}

/// Groups consecutive per-file `Will remove <f>` / `Removed <f>` lines belonging to the same
/// volume stem, keeping volume headers (`Will remove extra volume ...`) and volume order intact.
/// Returns `(canonical_multiset_lines, readdir_order_ok)`. When `positions` is `Some`, checks that
/// `.p??`/`.n??` lines within each volume appear in the exact `fs::read_dir` order of the test
/// copy, with `.tar.gz.md5` always last.
fn canonicalize_stdout(
    bytes: &[u8],
    cwd: &str,
    positions: Option<&BTreeMap<String, BTreeMap<String, usize>>>,
) -> (Vec<String>, bool) {
    let raw = String::from_utf8_lossy(bytes);
    let lines: Vec<&str> = raw.lines().collect();
    let mut out = Vec::new();
    let mut order_ok = true;
    let mut i = 0;
    while i < lines.len() {
        let line = lines[i];
        let prefix = if line.starts_with("Removed ") {
            Some("Removed ")
        } else if line.starts_with("Will remove ") && !line.starts_with("Will remove extra volume ")
        {
            Some("Will remove ")
        } else {
            None
        };
        let Some(prefix) = prefix else {
            out.push(line.to_string());
            i += 1;
            continue;
        };
        fn vol_stem(l: &str, prefix_len: usize) -> &str {
            let path = &l[prefix_len..];
            if let Some(s) = path.strip_suffix(".tar.gz.md5") {
                s
            } else if let Some(dot) = path.rfind('.') {
                &path[..dot]
            } else {
                path
            }
        }
        let stem = vol_stem(line, prefix.len());
        let mut group = Vec::new();
        while i < lines.len()
            && lines[i].starts_with(prefix)
            && vol_stem(lines[i], prefix.len()) == stem
        {
            group.push(lines[i].to_string());
            i += 1;
        }
        let mut non_md5 = Vec::new();
        let mut md5 = Vec::new();
        let mut seen_md5 = false;
        let mut last_pos: Option<usize> = None;
        for item in &group {
            if item.ends_with(".tar.gz.md5") {
                seen_md5 = true;
                md5.push(item.clone());
            } else {
                if seen_md5 {
                    order_ok = false;
                }
                if let Some(pos_map) = positions {
                    let (rel_dir, file_name) = resolve_rel_dir_and_name(&item[prefix.len()..], cwd);
                    match pos_map
                        .get(&rel_dir)
                        .and_then(|dir_map| dir_map.get(&file_name))
                    {
                        Some(&pos) => {
                            if let Some(prev) = last_pos {
                                if pos <= prev {
                                    order_ok = false;
                                }
                            }
                            last_pos = Some(pos);
                        }
                        None => order_ok = false,
                    }
                }
                non_md5.push(item.clone());
            }
        }
        let mut expected_non_md5 = non_md5;
        expected_non_md5.sort();
        expected_non_md5.extend(md5);
        out.extend(expected_non_md5);
    }
    (out, order_ok)
}

#[test]
fn cleanup_blastdb_volumes_matches_frozen_ncbi_oracle() {
    let mut differ = Vec::new();
    for case in cases() {
        let dir = scratch(&case.name);
        let _ = fs::remove_dir_all(&dir);
        copy_tree(&fixtures().join("in"), &dir.join("in"));
        let dir = dir.canonicalize().unwrap();
        let run_str = normalize_run_dir(&dir);
        let dir = PathBuf::from(&run_str);
        // Expand `<RUN>` in config fixture files.
        for entry in fs::read_dir(dir.join("in")).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir()
                && path
                    .file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(|n| n.starts_with("cfg_"))
            {
                for cfg in fs::read_dir(&path).unwrap() {
                    let cfg_path = cfg.unwrap().path();
                    if cfg_path.is_file() {
                        let content = fs::read_to_string(&cfg_path).unwrap();
                        fs::write(&cfg_path, content.replace("<RUN>", &run_str)).unwrap();
                    }
                }
            }
        }
        let before = listing(&dir);
        let positions = readdir_positions(&dir);
        let empty_home = dir.join("_empty_home");
        let mut cmd = base_command(&dir.join(&case.cwd), &empty_home);
        for (k, v) in &case.env {
            let expanded = v.replace("<RUN>", &run_str);
            cmd.env(k, &expanded);
            #[cfg(windows)]
            if k == "HOME" {
                cmd.env("USERPROFILE", &expanded);
            }
        }
        cmd.args(&case.args);
        let output = cmd.output().unwrap();
        let after = listing(&dir);

        let mut files = String::new();
        for (path, bytes) in &after {
            if before.get(path) != Some(bytes) {
                files.push_str(&format!(
                    "{path}\t{}\t{:016x}\n",
                    bytes.len(),
                    fnv1a64(bytes)
                ));
            }
        }
        for path in before.keys().filter(|path| !after.contains_key(*path)) {
            files.push_str(&format!("{path}\tREMOVED\t-\n"));
        }

        let got_out = replace_run_dir(&output.stdout, &run_str);
        let got_err = replace_run_dir(&output.stderr, &run_str);
        let got_rc = output.status.code().unwrap();

        let ncbi_file = |ext: &str| fixtures().join("ncbi").join(format!("{}.{ext}", case.name));
        let want_out =
            losat::wording::losat_output_bytes(&fs::read(ncbi_file("out")).unwrap()).into_owned();
        let want_err =
            losat::wording::losat_bytes(&fs::read(ncbi_file("err")).unwrap()).into_owned();
        let want_rc: i32 = fs::read_to_string(ncbi_file("rc"))
            .unwrap()
            .trim()
            .parse()
            .unwrap();
        let want_files = fs::read_to_string(ncbi_file("files")).unwrap();

        let (got_canon, readdir_ok) = canonicalize_stdout(&got_out, &case.cwd, Some(&positions));
        let (want_canon, _) = canonicalize_stdout(&want_out, &case.cwd, None);

        let stdout_ok = readdir_ok && got_canon == want_canon;
        if !stdout_ok || got_err != want_err || got_rc != want_rc || files != want_files {
            eprintln!(
                "--- {}: rc got={} want={} readdir_ok={}\n\
                 got stdout:  {:?}\nwant stdout: {:?}\n\
                 got stderr:  {:?}\nwant stderr: {:?}\n\
                 got files:\n{}want files:\n{}",
                case.name,
                got_rc,
                want_rc,
                readdir_ok,
                String::from_utf8_lossy(&got_out),
                String::from_utf8_lossy(&want_out),
                String::from_utf8_lossy(&got_err),
                String::from_utf8_lossy(&want_err),
                files,
                want_files
            );
            differ.push(case.name);
        } else {
            let _ = fs::remove_dir_all(&dir);
        }
    }
    assert!(differ.is_empty(), "rows that differ from NCBI: {differ:?}");
}

#[test]
fn version_prints_losat_version() {
    let dir = scratch("version");
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let expected = format!(
        "cleanup-blastdb-volumes.py {}\n",
        losat::wording::LOSAT_VERSION
    );
    for flag in ["-version", "-vers"] {
        let out = base_command(&dir, &dir).arg(flag).output().unwrap();
        assert_eq!(out.status.code(), Some(0), "{flag}: {out:?}");
        assert_eq!(String::from_utf8(out.stdout).unwrap(), expected, "{flag}");
        assert!(out.stderr.is_empty(), "{flag}");
    }
    let _ = fs::remove_dir_all(&dir);
}

#[cfg(unix)]
#[test]
fn broken_symlinks_in_orphan_volume_are_removed() {
    let dir = scratch("symlinks");
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    fs::write(dir.join("db.pal"), b"DBLIST db.00\n").unwrap();
    fs::write(dir.join("db.00.pin"), b"x").unwrap();
    fs::write(dir.join("db.01.pin"), b"x").unwrap();
    std::os::unix::fs::symlink("/nonexistent_target", dir.join("db.01.psq")).unwrap();
    std::os::unix::fs::symlink("/nonexistent_target", dir.join("db.01.tar.gz.md5")).unwrap();

    let mut want = String::new();
    for entry in fs::read_dir(&dir).unwrap().flatten() {
        let name = entry.file_name().to_string_lossy().into_owned();
        if name == "db.01.pin" || name == "db.01.psq" {
            want.push_str(&format!("Removed {name}\n"));
        }
    }
    want.push_str("Removed db.01.tar.gz.md5\n");

    let out = base_command(&dir, &dir)
        .args(["-db", "db", "-dbtype", "prot"])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(0), "broken symlinks: {out:?}");
    assert_eq!(String::from_utf8(out.stdout).unwrap(), want);
    assert!(out.stderr.is_empty());
    assert!(fs::symlink_metadata(dir.join("db.01.psq")).is_err());
    assert!(fs::symlink_metadata(dir.join("db.01.tar.gz.md5")).is_err());
    assert!(dir.join("db.00.pin").exists());
    let _ = fs::remove_dir_all(&dir);
}

#[cfg(all(unix, not(target_os = "macos")))]
#[test]
fn non_utf8_directory_path_and_volume_sort_match_python() {
    use std::ffi::OsString;
    use std::os::unix::ffi::OsStringExt;

    let dir = scratch("nonutf8");
    let _ = fs::remove_dir_all(&dir);
    let sub = dir.join(OsString::from_vec(b"dir_\xff".to_vec()));
    fs::create_dir_all(&sub).unwrap();
    fs::write(sub.join("mydb.pal"), b"DBLIST mydb.00\n").unwrap();
    fs::write(sub.join("mydb.00.pin"), b"x").unwrap();
    fs::write(
        sub.join(OsString::from_vec(b"mydb.\x80.pin".to_vec())),
        b"x",
    )
    .unwrap();
    fs::write(sub.join("mydb.\u{00e9}.pin"), b"x").unwrap();
    fs::write(sub.join("mydb.\u{1f600}.pin"), b"x").unwrap();

    let out = base_command(&dir, &dir)
        .arg("-db")
        .arg(OsString::from_vec(b"dir_\xff/mydb".to_vec()))
        .args(["-dbtype", "prot", "-dry-run"])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(0), "non-utf8: {out:?}");
    let mut expected = Vec::new();
    expected.extend_from_slice(b"Will remove extra volume dir_\xff/mydb.\xc3\xa9\n");
    expected.extend_from_slice(b"Will remove extra volume dir_\xff/mydb.\x80\n");
    expected.extend_from_slice(b"Will remove extra volume dir_\xff/mydb.\xf0\x9f\x98\x80\n");
    assert_eq!(out.stdout, expected);
    assert!(out.stderr.is_empty());
    let _ = fs::remove_dir_all(&dir);
}

#[cfg(target_os = "linux")]
#[test]
fn closed_pipe_and_dev_full_match_python_buffering() {
    let dir_small = scratch("full-small");
    let _ = fs::remove_dir_all(&dir_small);
    fs::create_dir_all(&dir_small).unwrap();
    fs::write(dir_small.join("db.pal"), b"DBLIST db.01\n").unwrap();
    for name in [
        "db.01.pin",
        "db.02.pin",
        "db.02.phr",
        "db.02.psq",
        "db.03.pin",
        "db.03.phr",
        "db.03.psq",
    ] {
        fs::write(dir_small.join(name), b"").unwrap();
    }
    let dev_full = fs::OpenOptions::new()
        .write(true)
        .open("/dev/full")
        .unwrap();
    let out = base_command(&dir_small, &dir_small)
        .args(["-db", "db", "-dbtype", "prot"])
        .stdout(Stdio::from(dev_full))
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(120));
    assert_eq!(
        String::from_utf8(out.stderr).unwrap(),
        "Exception ignored on flushing sys.stdout:\nOSError: [Errno 28] No space left on device\n"
    );
    assert_eq!(fs::read_dir(&dir_small).unwrap().count(), 2);
    let _ = fs::remove_dir_all(&dir_small);

    // `-h` and `-version` to `/dev/full` fail on exit flush with exit 120.
    let dir_hv = scratch("full-hv");
    let _ = fs::remove_dir_all(&dir_hv);
    fs::create_dir_all(&dir_hv).unwrap();
    for flag in ["-h", "-version"] {
        let dev_full = fs::OpenOptions::new()
            .write(true)
            .open("/dev/full")
            .unwrap();
        let out = base_command(&dir_hv, &dir_hv)
            .arg(flag)
            .stdout(Stdio::from(dev_full))
            .output()
            .unwrap();
        assert_eq!(out.status.code(), Some(120), "{flag}");
        assert_eq!(
            String::from_utf8(out.stderr).unwrap(),
            "Exception ignored on flushing sys.stdout:\nOSError: [Errno 28] No space left on device\n",
            "{flag}"
        );
    }
    let _ = fs::remove_dir_all(&dir_hv);

    // Output exceeding 8 KiB stops mid-run at exact 8192-byte flush boundary (410th of 600 files).
    let dir_large = scratch("full-large");
    let _ = fs::remove_dir_all(&dir_large);
    fs::create_dir_all(&dir_large).unwrap();
    fs::write(dir_large.join("db.pal"), b"DBLIST db.0000\n").unwrap();
    fs::write(dir_large.join("db.0000.pin"), b"").unwrap();
    for v in 1..=200 {
        for ext in ["pin", "p01", "p02"] {
            fs::write(dir_large.join(format!("db.{v:04}.{ext}")), b"").unwrap();
        }
    }
    let dev_full = fs::OpenOptions::new()
        .write(true)
        .open("/dev/full")
        .unwrap();
    let out = base_command(&dir_large, &dir_large)
        .args(["-db", "db", "-dbtype", "prot"])
        .stdout(Stdio::from(dev_full))
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(1));
    let stderr = String::from_utf8(out.stderr).unwrap();
    assert!(stderr.contains("Error: failed to write output:"));
    assert!(!stderr.contains("Exception ignored"));
    // 2 kept (`db.pal`, `db.0000.pin`) + 600 - 410 = 192 remaining files (matches Python 3.13).
    assert_eq!(fs::read_dir(&dir_large).unwrap().count(), 192);
    let _ = fs::remove_dir_all(&dir_large);
}

#[cfg(unix)]
#[test]
fn removed_cwd_fails_when_config_fallback_runs() {
    let dir = scratch("removed-cwd");
    let _ = fs::remove_dir_all(&dir);
    let gone = dir.join("gone");
    fs::create_dir_all(&gone).unwrap();
    let bin = env!("CARGO_BIN_EXE_losat");
    let out = Command::new("/bin/sh")
        .arg("-c")
        .arg("cd \"$1\" && rmdir \"$1\" && exec \"$2\" cleanup-blastdb-volumes.py -db missing -dbtype prot")
        .arg("sh")
        .arg(&gone)
        .arg(bin)
        .env("HOME", &dir)
        .env_remove("BLASTDB")
        .env_remove("NCBI")
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(1));
    assert!(out.stdout.is_empty());
    assert!(
        String::from_utf8_lossy(&out.stderr).contains("cannot determine current working directory")
    );
    let _ = fs::remove_dir_all(&dir);
}
