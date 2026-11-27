//! Frozen official NCBI BLAST+ 2.17.0 comparison tests for
//! `losat makeblastdb -input_type asn1_txt` and `-input_type asn1_bin`.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use super::*;
use crate::utils::test_dir::TempDir;
use crate::writedb::tests::parse_makedb_words;

const FIXED_TIME: LocalTime = LocalTime {
    year: 2026,
    month: 11,
    day: 25,
    hour: 7,
    minute: 0,
    second: 0,
    unix_seconds: 1_795_503_600,
};

fn asn1_fixture_paths() -> (PathBuf, PathBuf, String, String) {
    let fix_a = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/asn1");
    let fix_m = fix_a.join("makeblastdb/in");
    let fix_a_str = fix_a.to_string_lossy().into_owned();
    let fix_m_str = fix_m.to_string_lossy().into_owned();
    (fix_a, fix_m, fix_a_str, fix_m_str)
}

fn dir_string(dir: &TempDir) -> String {
    dir.to_string_lossy().into_owned()
}

fn assert_refusal(label: &str, rc: i32, stderr: &[u8], expected_substr: &str) {
    let err_text = String::from_utf8_lossy(stderr);
    assert_eq!(rc, 1, "{label}: expected rc 1, err={err_text}");
    assert!(
        err_text.starts_with("losat makeblastdb:"),
        "{label}: {err_text}"
    );
    assert!(err_text.contains(expected_substr), "{label}: {err_text}");
}

fn sha256_hex(data: &[u8]) -> String {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    let bit_len = (data.len() as u64).wrapping_mul(8);
    let mut msg = data.to_vec();
    msg.push(0x80);
    while msg.len() % 64 != 56 {
        msg.push(0);
    }
    msg.extend_from_slice(&bit_len.to_be_bytes());
    for chunk in msg.chunks_exact(64) {
        let mut w = [0u32; 64];
        for i in 0..16 {
            w[i] = u32::from_be_bytes(chunk[i * 4..i * 4 + 4].try_into().unwrap());
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let [mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut hh] = h;
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ ((!e) & g);
            let t1 = hh
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[i])
                .wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let t2 = s0.wrapping_add(maj);
            hh = g;
            g = f;
            f = e;
            e = d.wrapping_add(t1);
            d = c;
            c = b;
            b = a;
            a = t1.wrapping_add(t2);
        }
        for (dst, add) in h.iter_mut().zip([a, b, c, d, e, f, g, hh]) {
            *dst = dst.wrapping_add(add);
        }
    }
    let mut out = String::with_capacity(64);
    for word in h {
        use std::fmt::Write;
        let _ = write!(out, "{word:08x}");
    }
    out
}

fn decode_hex(hex: &str) -> Vec<u8> {
    assert_eq!(hex.len() % 2, 0);
    (0..hex.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&hex[i..i + 2], 16).unwrap())
        .collect()
}

fn norm_index(b: &[u8]) -> Vec<u8> {
    if b.len() < 20 {
        return b.to_vec();
    }
    let ver = u32::from_be_bytes(b[0..4].try_into().unwrap());
    let mut off = 8 + if ver == 5 { 4 } else { 0 };
    if off + 4 > b.len() {
        return b.to_vec();
    }
    let tlen = u32::from_be_bytes(b[off..off + 4].try_into().unwrap()) as usize;
    off += 4 + tlen;
    if ver == 5 {
        if off + 4 > b.len() {
            return b.to_vec();
        }
        let lmdb_len = u32::from_be_bytes(b[off..off + 4].try_into().unwrap()) as usize;
        off += 4 + lmdb_len;
    }
    if off + 4 > b.len() {
        return b.to_vec();
    }
    let dlen = u32::from_be_bytes(b[off..off + 4].try_into().unwrap()) as usize;
    if off + 4 + dlen > b.len() {
        return b.to_vec();
    }
    let mut out = Vec::with_capacity(b.len() - dlen + 6);
    out.extend_from_slice(&b[..off]);
    out.extend_from_slice(b"<DATE>");
    out.extend_from_slice(&b[off + 4 + dlen..]);
    out
}

fn norm_col_index(b: &[u8]) -> Vec<u8> {
    if b.len() < 33 {
        return b.to_vec();
    }
    let meta_off = u32::from_be_bytes(b[24..28].try_into().unwrap()) as usize;
    let array_off = u32::from_be_bytes(b[28..32].try_into().unwrap()) as usize;
    let title_len = usize::from(b[32]);
    let title_end = 33 + title_len;
    if title_end > meta_off || meta_off > array_off || array_off > b.len() {
        return b.to_vec();
    }
    let mut meta_end = array_off;
    while meta_end > meta_off && b[meta_end - 1] == 0 {
        meta_end -= 1;
    }
    let mut out = Vec::with_capacity(b.len());
    out.extend_from_slice(&b[..24]);
    out.extend_from_slice(&b[32..title_end]);
    out.extend_from_slice(b"<DATE>");
    out.extend_from_slice(&b[meta_off..meta_end]);
    out.extend_from_slice(b"<PAD>");
    out.extend_from_slice(&b[array_off..]);
    out
}

fn norm_js(b: &[u8]) -> Vec<u8> {
    let needle = b"\"last-updated\": \"";
    let Some(start) = b.windows(needle.len()).position(|w| w == needle) else {
        return b.to_vec();
    };
    let val_start = start + needle.len();
    let Some(rel_end) = b[val_start..].iter().position(|&c| c == b'"') else {
        return b.to_vec();
    };
    let mut out = Vec::with_capacity(b.len() + 6);
    out.extend_from_slice(&b[..val_start]);
    out.extend_from_slice(b"<DATE>");
    out.extend_from_slice(&b[val_start + rel_end..]);
    out
}

fn norm_log(b: &[u8], dir_str: &str, fix_a_str: &str, fix_m_str: &str) -> Vec<u8> {
    let mut s = String::from_utf8_lossy(b).into_owned();
    s = s
        .replace(dir_str, "<DIR>")
        .replace(fix_m_str, "<M>")
        .replace(fix_a_str, "<A>")
        .replace("<DIR>\\", "<DIR>/")
        .replace("<M>\\", "<M>/")
        .replace("<A>\\", "<A>/");
    if let Some(pos) = s.find("Building a new DB, current time: ") {
        let start = pos + "Building a new DB, current time: ".len();
        if start + 19 <= s.len() {
            s.replace_range(start..start + 19, "<TIME>");
        }
    }
    let mut out = String::with_capacity(s.len());
    let mut rest = s.as_str();
    while let Some(idx) = rest.find(" in ") {
        let after = &rest[idx + 4..];
        if let Some(sec_pos) = after.find(" seconds.") {
            let num = &after[..sec_pos];
            if !num.is_empty()
                && num
                    .bytes()
                    .all(|c| c.is_ascii_digit() || matches!(c, b'.' | b'e' | b'E' | b'+' | b'-'))
            {
                out.push_str(&rest[..idx + 4]);
                out.push_str("<SEC> seconds.");
                rest = &after[sec_pos + " seconds.".len()..];
                continue;
            }
        }
        out.push_str(&rest[..idx + 4]);
        rest = after;
    }
    out.push_str(rest);
    out.into_bytes()
}

fn replace_bytes(haystack: &[u8], needle: &[u8], replacement: &[u8]) -> Vec<u8> {
    if needle.is_empty() {
        return haystack.to_vec();
    }
    let mut out = Vec::with_capacity(haystack.len());
    let mut i = 0;
    while i < haystack.len() {
        if haystack[i..].starts_with(needle) {
            out.extend_from_slice(replacement);
            i += needle.len();
        } else {
            out.push(haystack[i]);
            i += 1;
        }
    }
    out
}

fn snapshot_dir(
    dir: &Path,
    dir_str: &str,
    fix_a_str: &str,
    fix_m_str: &str,
) -> BTreeMap<String, (usize, String)> {
    let mut map = BTreeMap::new();
    for entry in std::fs::read_dir(dir).unwrap() {
        let entry = entry.unwrap();
        if !entry.file_type().unwrap().is_file() {
            continue;
        }
        let name = entry.file_name().to_string_lossy().into_owned();
        if name.ends_with("-lock") {
            continue;
        }
        let raw = std::fs::read(entry.path()).unwrap();
        let norm = if name.ends_with(".nin") || name.ends_with(".pin") {
            norm_index(&raw)
        } else if name.ends_with(".naa") || name.ends_with(".paa") {
            norm_col_index(&raw)
        } else if name.ends_with(".njs") || name.ends_with(".pjs") {
            norm_js(&raw)
        } else if name.ends_with(".log") {
            norm_log(&raw, dir_str, fix_a_str, fix_m_str)
        } else {
            raw
        };
        map.insert(name, (norm.len(), sha256_hex(&norm)));
    }
    map
}

fn run_case(args_words: &[String], env: MakedbEnv, stdin_bytes: &[u8]) -> (i32, Vec<u8>, Vec<u8>) {
    let args = parse_makedb_words(args_words);
    let mut stdin = std::io::Cursor::new(stdin_bytes);
    let mut stdout = Vec::new();
    let mut stderr = Vec::new();
    let rc = run_with_env(
        &args,
        env,
        &ApplicationSettings::default(),
        Some(&FIXED_TIME),
        Platform::LINUX_X86_64,
        &mut stdin,
        &mut stdout,
        &mut stderr,
    );
    (rc, stdout, stderr)
}

#[test]
fn all_frozen_ncbi_2_17_0_asn1_cases_match_byte_for_byte() {
    let (fix_a, _fix_m, fix_a_str, fix_m_str) = asn1_fixture_paths();
    let manifest_path = fix_a.join("makeblastdb/ncbi_2_17_0_manifest.tsv");
    let manifest = std::fs::read_to_string(&manifest_path).unwrap();

    let mut checked = 0usize;
    for line in manifest.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 11, "malformed manifest row: {line}");
        let [name, in_spec, itype, dbtype, extra_str, env_spec, stdin_file, rc_str, out_hex, err_hex, files_spec] =
            cols.as_slice()
        else {
            unreachable!();
        };
        let tmp = TempDir::unique(&format!("asn1-frozen-{name}"));
        let dir_str = dir_string(&tmp);
        let subst = |s: &str| {
            s.replace("@A", &fix_a_str)
                .replace("@M", &fix_m_str)
                .replace("@D", &dir_str)
        };
        let mut words = vec![
            "-in".into(),
            subst(in_spec),
            "-input_type".into(),
            (*itype).into(),
            "-dbtype".into(),
            (*dbtype).into(),
            "-title".into(),
            "T".into(),
            "-out".into(),
            format!("{dir_str}/db"),
        ];
        for tok in extra_str.split_whitespace() {
            words.push(subst(tok));
        }
        let env = MakedbEnv {
            skip_unverified: *env_spec == "SKIP_UNVER=1",
            scan_cfastareader: *env_spec != "NOSCAN=1",
        };
        let stdin_bytes = if stdin_file.is_empty() {
            Vec::new()
        } else {
            std::fs::read(subst(stdin_file)).unwrap()
        };

        let (rc, stdout, stderr) = run_case(&words, env, &stdin_bytes);
        let expected_rc: i32 = rc_str.parse().unwrap();
        assert_eq!(rc, expected_rc, "case {name}: rc mismatch");

        let actual_out = norm_log(&stdout, &dir_str, &fix_a_str, &fix_m_str);
        let expected_out = decode_hex(out_hex);
        assert_eq!(
            String::from_utf8_lossy(&actual_out),
            String::from_utf8_lossy(&expected_out),
            "case {name}: stdout mismatch"
        );

        let raw_expected_err = decode_hex(err_hex);
        let expected_err = crate::wording::losat_bytes(&raw_expected_err);
        let actual_err = replace_bytes(
            &replace_bytes(
                &replace_bytes(&stderr, dir_str.as_bytes(), b"<DIR>"),
                fix_m_str.as_bytes(),
                b"<M>",
            ),
            fix_a_str.as_bytes(),
            b"<A>",
        );
        assert_eq!(
            String::from_utf8_lossy(&actual_err),
            String::from_utf8_lossy(&expected_err),
            "case {name}: stderr mismatch"
        );

        let actual_files = snapshot_dir(&tmp, &dir_str, &fix_a_str, &fix_m_str);
        let mut expected_files = BTreeMap::new();
        if !files_spec.is_empty() {
            for item in files_spec.split(',') {
                let mut parts = item.splitn(3, ':');
                let fname = parts.next().unwrap().to_string();
                let flen: usize = parts.next().unwrap().parse().unwrap();
                let fsha = parts.next().unwrap().to_string();
                expected_files.insert(fname, (flen, fsha));
            }
        }
        assert_eq!(actual_files, expected_files, "case {name}: DB files differ");
        checked += 1;
    }
    assert!(
        checked >= 129,
        "expected >= 129 frozen cases, got {checked}"
    );
}

/// The objects that NCBI's makeblastdb reads with a report of memory addresses (an empty region's
/// post before the diagnostic stream failed, the `CCoreException::eNullPtr` report of an element
/// stated as NULL) or crashes on (a division by zero): LOSAT refuses each in its words, with exit
/// status 1, writes exactly the lines of NCBI's log that come before that point, and leaves no
/// file of the database. The table holds NCBI 2.17.0's exit status, log and standard error of
/// each, and the number of lines of the log written before the report or the crash
/// (`stdout_lines_before`, from a run with both streams on one pipe).
#[test]
fn object_manager_address_and_crash_cases_are_refused() {
    use crate::objmgr::annot_index::Refusal;
    let (fix_a, _fix_m, fix_a_str, fix_m_str) = asn1_fixture_paths();
    let table = std::fs::read_to_string(fix_a.join("makeblastdb/ncbi_2_17_0_refused.tsv")).unwrap();
    let mut checked = 0usize;
    for line in table.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        let [name, in_spec, itype, dbtype, refusal, ncbi_rc, ncbi_out_hex, ncbi_err_hex, lines_before] =
            cols.as_slice()
        else {
            panic!("malformed row: {line}");
        };
        let refusal = match *refusal {
            "EmptyRegion" => Refusal::EmptyRegion,
            "NullElement" => Refusal::NullElement,
            "DivideByZero" => Refusal::DivideByZero,
            other => panic!("case {name}: refusal {other}"),
        };
        // The row is one that NCBI writes with memory addresses or crashes on.
        let ncbi_err = String::from_utf8_lossy(&decode_hex(ncbi_err_hex)).into_owned();
        match refusal {
            Refusal::EmptyRegion => assert!(
                ncbi_err
                    .lines()
                    .find(|line| line.contains("Empty region in "))
                    .is_some_and(|line| line.contains("0x")),
                "case {name}"
            ),
            Refusal::NullElement => assert!(
                *ncbi_rc == "255"
                    && ncbi_err.contains("(CCoreException::eNullPtr)")
                    && ncbi_err.contains(" addr=0x"),
                "case {name}"
            ),
            _ => assert_eq!(*ncbi_rc, "136", "case {name}"),
        }
        let tmp = TempDir::unique(&format!("asn1-refused-{name}"));
        let dir_str = dir_string(&tmp);
        let words = vec![
            "-in".into(),
            in_spec.replace("@M", &fix_m_str).replace("@A", &fix_a_str),
            "-input_type".into(),
            (*itype).into(),
            "-dbtype".into(),
            (*dbtype).into(),
            "-title".into(),
            "T".into(),
            "-out".into(),
            format!("{dir_str}/db"),
        ];
        let env = MakedbEnv {
            skip_unverified: false,
            scan_cfastareader: true,
        };
        let (rc, stdout, stderr) = run_case(&words, env, &[]);
        assert_eq!(rc, 1, "case {name}: {}", String::from_utf8_lossy(&stderr));
        assert_eq!(
            String::from_utf8_lossy(&stderr),
            format!(
                "losat makeblastdb: {} is not supported by LOSAT's makeblastdb\n",
                refusal.text(crate::writedb::asn_source::ASN_OBJECT)
            ),
            "case {name}: stderr"
        );
        let actual_out = norm_log(&stdout, &dir_str, &fix_a_str, &fix_m_str);
        let ncbi_out = decode_hex(ncbi_out_hex);
        let lines_before: usize = lines_before.parse().unwrap();
        let expected_out: Vec<u8> = ncbi_out
            .split_inclusive(|&b| b == b'\n')
            .take(lines_before)
            .flatten()
            .copied()
            .collect();
        assert_eq!(
            String::from_utf8_lossy(&actual_out),
            String::from_utf8_lossy(&expected_out),
            "case {name}: the log before the refusal"
        );
        assert!(
            snapshot_dir(&tmp, &dir_str, &fix_a_str, &fix_m_str).is_empty(),
            "case {name}: files of the database were left"
        );
        checked += 1;
    }
    assert!(checked >= 45, "expected >= 45 refused cases, got {checked}");
}

/// An object refused after a sequence was written leaves no file of the build: neither those of
/// NCBI's `ListFiles` nor the others that a failed build of NCBI's leaves (the LMDB writers'
/// `.nos`, `.not`, `.ntf` and `.nto`, the GI-mask files and the alias file of a single volume
/// with GI masks; `WriteDb::every_file`). The build of the first object alone shows that each
/// of them is written.
#[test]
fn a_refused_object_after_a_written_sequence_leaves_no_file() {
    use crate::objmgr::annot_index::Refusal;
    let first =
        "Seq-entry ::= seq {\n  id { gi 5 },\n  descr { title \"t\" },\n  inst { repr raw, \
                 mol dna, length 60, seq-data iupacna \
                 \"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\" }\n}\n";
    let refused = "Seq-entry ::= seq {\n  id { local str \"s2\" },\n  descr { title \"t\" },\n  \
                   inst { repr raw, mol dna, length 60, seq-data iupacna \
                   \"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\" },\n  annot { { \
                   data ftable {\n    { data imp { key \"misc_feature\" }, location mix { empty \
                   local str \"X\", int { from 0, to 29, id local str \"s2\" } } }\n  } } }\n}\n";
    let masks = "Blast-db-mask-info ::= {\n  algo-id 0,\n  algo-program dust,\n  algo-options \
                 \"window=64; level=20; linker=1\",\n  masks {\n    masks {\n      int { from 10, \
                 to 40, id gi 5 }\n    },\n    more FALSE\n  }\n}\n";
    let inputs = TempDir::unique("asn1-refused-every-file-in");
    std::fs::write(inputs.join("one.asn"), first).unwrap();
    std::fs::write(inputs.join("two.asn"), [first, refused].concat()).unwrap();
    std::fs::write(inputs.join("masks.asn"), masks).unwrap();
    let path = |dir: &Path, name: &str| dir.join(name).to_string_lossy().into_owned();
    let build = |input: &str| {
        let out = TempDir::unique(&format!("asn1-refused-every-file-{input}"));
        let words: Vec<String> = vec![
            "-in".into(),
            path(&inputs, input),
            "-input_type".into(),
            "asn1_txt".into(),
            "-dbtype".into(),
            "nucl".into(),
            "-parse_seqids".into(),
            "-mask_data".into(),
            path(&inputs, "masks.asn"),
            "-gi_mask".into(),
            "-gi_mask_name".into(),
            path(&out, "gm"),
            "-title".into(),
            "T".into(),
            "-out".into(),
            path(&out, "db"),
        ];
        let env = MakedbEnv {
            skip_unverified: false,
            scan_cfastareader: true,
        };
        let (rc, _, stderr) = run_case(&words, env, &[]);
        let mut files: Vec<String> = std::fs::read_dir(&*out)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
            .filter(|name| !name.ends_with("-lock"))
            .collect();
        files.sort();
        (rc, String::from_utf8_lossy(&stderr).into_owned(), files)
    };
    let (rc, stderr, files) = build("one.asn");
    assert_eq!(rc, 0, "{stderr}");
    for name in [
        "db.nal", "db.nos", "db.not", "db.ntf", "db.nto", "gm.gmd", "gm.gmi", "gm.gmo", "gm.gnd",
        "gm.gni", "gm.gno",
    ] {
        assert!(files.iter().any(|file| file == name), "{name}: {files:?}");
    }
    let (rc, stderr, files) = build("two.asn");
    assert_eq!(rc, 1, "{stderr}");
    assert_eq!(
        stderr,
        format!(
            "losat makeblastdb: {} is not supported by LOSAT's makeblastdb\n",
            Refusal::EmptyRegion.text(crate::writedb::asn_source::ASN_OBJECT)
        )
    );
    assert!(files.is_empty(), "{files:?}");
}

#[test]
fn text_and_binary_asn1_produce_identical_databases() {
    let (fix_a, fix_m, fix_a_str, fix_m_str) = asn1_fixture_paths();
    let pairs = [
        ("e1", fix_a.join("e1.asn"), fix_a.join("e1.asnb"), "nucl"),
        (
            "np_n",
            fix_a.join("nucprot.asn"),
            fix_a.join("nucprot.asnb"),
            "nucl",
        ),
        (
            "np_p",
            fix_a.join("nucprot.asn"),
            fix_a.join("nucprot.asnb"),
            "prot",
        ),
        ("ids", fix_a.join("ids.asn"), fix_a.join("ids.asnb"), "nucl"),
        (
            "ptitles",
            fix_a.join("ptitles.asn"),
            fix_a.join("ptitles.asnb"),
            "prot",
        ),
        ("fr", fix_a.join("fr.asn"), fix_a.join("fr.asnb"), "nucl"),
        (
            "multi",
            fix_a.join("multi.asn"),
            fix_a.join("multi.asnb"),
            "nucl",
        ),
        (
            "big1",
            fix_m.join("big1.asn"),
            fix_m.join("big1.asnb"),
            "nucl",
        ),
        (
            "big4",
            fix_m.join("big4.asn"),
            fix_m.join("big4.asnb"),
            "nucl",
        ),
        (
            "gz",
            fix_m.join("e1.asn.gz"),
            fix_m.join("e1.asnb.gz"),
            "nucl",
        ),
        (
            "zst",
            fix_m.join("e1.asn.zst"),
            fix_m.join("e1.asnb.zst"),
            "nucl",
        ),
    ];
    let env = MakedbEnv {
        skip_unverified: false,
        scan_cfastareader: true,
    };
    for (label, txt_path, bin_path, dbtype) in pairs {
        let t_dir = TempDir::unique(&format!("asn1-twin-txt-{label}"));
        let b_dir = TempDir::unique(&format!("asn1-twin-bin-{label}"));
        let t_str = dir_string(&t_dir);
        let b_str = dir_string(&b_dir);
        let t_words = vec![
            "-in".into(),
            txt_path.to_string_lossy().into_owned(),
            "-input_type".into(),
            "asn1_txt".into(),
            "-dbtype".into(),
            dbtype.into(),
            "-title".into(),
            "T".into(),
            "-out".into(),
            format!("{t_str}/db"),
        ];
        let b_words = vec![
            "-in".into(),
            bin_path.to_string_lossy().into_owned(),
            "-input_type".into(),
            "asn1_bin".into(),
            "-dbtype".into(),
            dbtype.into(),
            "-title".into(),
            "T".into(),
            "-out".into(),
            format!("{b_str}/db"),
        ];
        let (trc, _, terr) = run_case(&t_words, env, &[]);
        let (brc, _, berr) = run_case(&b_words, env, &[]);
        assert_eq!(
            trc,
            0,
            "{label} txt failed: {}",
            String::from_utf8_lossy(&terr)
        );
        assert_eq!(
            brc,
            0,
            "{label} bin failed: {}",
            String::from_utf8_lossy(&berr)
        );
        assert_eq!(
            snapshot_dir(&t_dir, &t_str, &fix_a_str, &fix_m_str),
            snapshot_dir(&b_dir, &b_str, &fix_a_str, &fix_m_str),
            "text vs binary DB mismatch for {label}"
        );
    }
}

#[test]
fn d3c4_unsupported_title_and_ub_cases_are_refused_cleanly() {
    let (fix_a, fix_m, _, _) = asn1_fixture_paths();
    let env_default = MakedbEnv {
        skip_unverified: false,
        scan_cfastareader: true,
    };
    let file_cases = [
        (
            "molinfo",
            fix_a.join("titles.asn"),
            "nucl",
            env_default,
            "MolInfo",
        ),
        (
            "unver",
            fix_m.join("unver.asn"),
            "nucl",
            env_default,
            "Unverified",
        ),
        (
            "protfeat",
            fix_m.join("protfeat.asn"),
            "prot",
            env_default,
            "not supported by LOSAT's makeblastdb",
        ),
        (
            "llong",
            fix_m.join("llong.asn"),
            "nucl",
            env_default,
            "past the data",
        ),
    ];
    for (label, path, dbtype, env, expected_substr) in file_cases {
        let dir = TempDir::unique(&format!("asn1-refuse-{label}"));
        let dir_str = dir_string(&dir);
        let words = vec![
            "-in".into(),
            path.to_string_lossy().into_owned(),
            "-input_type".into(),
            "asn1_txt".into(),
            "-dbtype".into(),
            dbtype.into(),
            "-title".into(),
            "T".into(),
            "-out".into(),
            format!("{dir_str}/db"),
        ];
        let (rc, _, stderr) = run_case(&words, env, &[]);
        assert_refusal(label, rc, &stderr, expected_substr);
        for entry in std::fs::read_dir(&dir).unwrap() {
            let fname = entry.unwrap().file_name().to_string_lossy().into_owned();
            assert!(
                !fname.ends_with(".nin")
                    && !fname.ends_with(".pin")
                    && !fname.ends_with(".nhr")
                    && !fname.ends_with(".phr")
                    && !fname.ends_with(".nsq")
                    && !fname.ends_with(".psq")
                    && !fname.ends_with(".ndb")
                    && !fname.ends_with(".pdb"),
                "{label}: volume/LMDB file {fname} was left behind after refusal"
            );
        }
    }

    // Regression: Unverified on a genbank set (not skipped by SKIP_UNVER) whose explicit title
    // begins with "UNVERIFIED: " is still refused because NCBI strips the prefix and re-emits
    // UNVERIFIED_ORG: from the user object.
    let genbank_unver = br#"Seq-entry ::= set {
  class genbank,
  descr { user { type str "Unverified", data { { label str "Type", data str "Organism" } } } },
  seq-set {
    seq { id { local str "g1" }, descr { title "UNVERIFIED: ok" }, inst { repr raw, mol dna, length 4, seq-data iupacna "ACGT" } }
  }
}
"#;
    let empty_delta = br#"Seq-entry ::= seq {
  id { local str "ed" },
  descr { title "t" },
  inst { repr delta, mol dna, length 4, ext delta {} }
}
"#;
    let null_id = br#"Seq-entry ::= seq {
  id { NULL },
  descr { title "t" },
  inst { repr raw, mol dna, length 4, seq-data iupacna "ACGT" }
}
"#;
    for (label, bytes, env, substr) in [
        (
            "genbank_unver_prefix",
            genbank_unver.as_slice(),
            MakedbEnv {
                skip_unverified: true,
                scan_cfastareader: true,
            },
            "Unverified",
        ),
        (
            "empty_delta",
            empty_delta.as_slice(),
            env_default,
            "empty Delta-ext",
        ),
        (
            "null_id",
            null_id.as_slice(),
            env_default,
            "an ASN.1 input object with an element stated as NULL in a list that NCBI's makeblastdb reads",
        ),
    ] {
        let dir = TempDir::unique(&format!("asn1-refuse-inline-{label}"));
        let dir_str = dir_string(&dir);
        let words = vec![
            "-in".into(),
            "-".into(),
            "-input_type".into(),
            "asn1_txt".into(),
            "-dbtype".into(),
            "nucl".into(),
            "-title".into(),
            "T".into(),
            "-out".into(),
            format!("{dir_str}/db"),
        ];
        let (rc, _, stderr) = run_case(&words, env, bytes);
        assert_refusal(label, rc, &stderr, substr);
    }
}

#[test]
fn check_environment_allows_cfastareader_override_for_asn1_inputs() {
    const DONT_SCAN: &str = "NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS";
    let with_dont_scan = |name: &str| (name == DONT_SCAN).then(|| std::ffi::OsString::from("1"));
    for input_type in ["asn1_txt", "asn1_bin"] {
        let mut args = parse_makedb_words(&["-dbtype".into(), "nucl".into()]);
        args.input_type = input_type.to_string();
        assert!(
            check_environment(&args, &with_dont_scan).is_ok(),
            "{input_type} should allow {DONT_SCAN}"
        );
    }
}
