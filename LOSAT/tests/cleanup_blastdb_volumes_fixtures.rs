//! `losat cleanup-blastdb-volumes.py` against NCBI BLAST+ 2.17.0 `cleanup-blastdb-volumes.py`:
//! each row of `tests/fixtures/cleanup_blastdb_volumes/cases.tsv` runs in a fresh copy of `in/`
//! and checks stdout, stderr, exit status, and exact file deletions against the frozen official
//! NCBI oracle (`ncbi/<row>.{out,err,rc,files}`). NCBI is not run here.
//!
//! Within each extra volume, NCBI's Python `glob.glob` visits `.p??`/`.n??` files in filesystem
//! directory order before `.tar.gz.md5`, while LOSAT sorts `.p??`/`.n??` lexicographically per
//! volume before `.tar.gz.md5`. The comparison below checks both that LOSAT's output is strictly
//! sorted per volume (with `.tar.gz.md5` last) and that every volume's lines equal NCBI's frozen
//! lines as a multiset, while volume order, non-verbose `-dry-run`, `-h`, stderr, exit code, and
//! removed/kept files match byte for byte.

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

fn listing(dir: &Path) -> BTreeMap<String, Vec<u8>> {
    let mut found = BTreeMap::new();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(at) = stack.pop() {
        for entry in fs::read_dir(&at).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                stack.push(path);
            } else {
                let relative = path.strip_prefix(dir).unwrap();
                let parts: Vec<String> = relative
                    .components()
                    .map(|part| part.as_os_str().to_string_lossy().into_owned())
                    .collect();
                found.insert(parts.join("/"), fs::read(&path).unwrap());
            }
        }
    }
    found
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

/// Groups consecutive per-file `Will remove <f>` / `Removed <f>` lines belonging to the same
/// volume stem, keeping volume headers (`Will remove extra volume ...`) and volume order intact.
/// Returns `(volume_groups, losat_is_sorted_with_md5_last)`.
fn canonicalize_stdout(bytes: &[u8]) -> (Vec<String>, bool) {
    let raw = String::from_utf8_lossy(bytes);
    let lines: Vec<&str> = raw.lines().collect();
    let mut out = Vec::new();
    let mut sorted_ok = true;
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
        // Check that non-md5 lines are sorted and md5 (if present) comes last.
        let mut non_md5 = Vec::new();
        let mut md5 = Vec::new();
        let mut seen_md5 = false;
        for item in &group {
            if item.ends_with(".tar.gz.md5") {
                seen_md5 = true;
                md5.push(item.clone());
            } else {
                if seen_md5 {
                    sorted_ok = false;
                }
                non_md5.push(item.clone());
            }
        }
        let mut expected_non_md5 = non_md5.clone();
        expected_non_md5.sort();
        if non_md5 != expected_non_md5 {
            sorted_ok = false;
        }
        expected_non_md5.extend(md5);
        out.extend(expected_non_md5);
    }
    (out, sorted_ok)
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

        let (got_canon, got_sorted) = canonicalize_stdout(&got_out);
        let (want_canon, _) = canonicalize_stdout(&want_out);

        let stdout_ok = got_sorted && got_canon == want_canon;
        if !stdout_ok || got_err != want_err || got_rc != want_rc || files != want_files {
            eprintln!(
                "--- {}: rc got={} want={} sorted_ok={}\n\
                 got stdout:  {:?}\nwant stdout: {:?}\n\
                 got stderr:  {:?}\nwant stderr: {:?}\n\
                 got files:\n{}want files:\n{}",
                case.name,
                got_rc,
                want_rc,
                got_sorted,
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

    let out = base_command(&dir, &dir)
        .args(["-db", "db", "-dbtype", "prot"])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(0), "broken symlinks: {out:?}");
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        "Removed db.01.pin\nRemoved db.01.psq\nRemoved db.01.tar.gz.md5\n"
    );
    assert!(out.stderr.is_empty());
    assert!(fs::symlink_metadata(dir.join("db.01.psq")).is_err());
    assert!(fs::symlink_metadata(dir.join("db.01.tar.gz.md5")).is_err());
    assert!(dir.join("db.00.pin").exists());
    let _ = fs::remove_dir_all(&dir);
}

#[cfg(all(unix, not(target_os = "macos")))]
#[test]
fn non_utf8_directory_path_works_with_raw_bytes() {
    use std::ffi::OsString;
    use std::os::unix::ffi::OsStringExt;

    let dir = scratch("nonutf8");
    let _ = fs::remove_dir_all(&dir);
    let sub = dir.join(OsString::from_vec(b"dir_\xff".to_vec()));
    fs::create_dir_all(&sub).unwrap();
    fs::write(sub.join("mydb.pal"), b"DBLIST mydb.00\n").unwrap();
    fs::write(sub.join("mydb.00.pin"), b"x").unwrap();
    fs::write(sub.join("mydb.01.pin"), b"x").unwrap();

    let out = base_command(&dir, &dir)
        .arg("-db")
        .arg(OsString::from_vec(b"dir_\xff/mydb".to_vec()))
        .args(["-dbtype", "prot", "-dry-run", "-verbose"])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(0), "non-utf8: {out:?}");
    assert_eq!(
        out.stdout,
        b"Will remove extra volume dir_\xff/mydb.01\nWill remove dir_\xff/mydb.01.pin\n"
    );
    assert!(out.stderr.is_empty());
    let _ = fs::remove_dir_all(&dir);
}
