//! `makeblastdb -input_type blastdb` through `run_with`: a database written from FASTA and its
//! copy are the same files, and the errors of the command line and of the input's names are
//! makeblastdb 2.17.0's (probes of line RG item D3c, `D/d3/c/probe/out_090032`, Linux x86-64).

use std::path::Path;

use crate::writedb::makeblastdb::run_with;
use crate::writedb::tests::{parse_makedb_words, seqdb_fixtures, time_of, TempDir, LINUX};

/// The run of makeblastdb with `words`: exit status, standard output, standard error.
fn makeblastdb(words: &[&str]) -> (i32, String, String) {
    let words: Vec<String> = words.iter().map(|word| word.to_string()).collect();
    let args = parse_makedb_words(&words);
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let status = run_with(
        &args,
        None,
        Some(&time_of("Nov 21, 2026  9:00 AM")),
        LINUX,
        &mut std::io::empty(),
        &mut out,
        &mut err,
    );
    let text = |bytes: Vec<u8>| String::from_utf8(bytes).unwrap();
    (status, text(out), text(err))
}

/// The files of `dir`, sorted.
fn files_of(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .map(|entries| {
            entries
                .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
                .collect()
        })
        .unwrap_or_default();
    names.sort();
    names
}

/// The text mask file `file` of the fixtures with its Seq-ids `local id k` as the parsed ids of
/// the FASTA file's `k`th sequence, `local str "m0k"` (`mask_nucl.fna` names them `m01` to
/// `m08`).
fn masks_by_name(file: &str) -> String {
    let text = std::fs::read_to_string(seqdb_fixtures().join(file)).unwrap();
    let mut out = String::new();
    let mut rest = text.as_str();
    while let Some(at) = rest.find("local id ") {
        out.push_str(&rest[..at]);
        rest = &rest[at + "local id ".len()..];
        let digits = rest.bytes().take_while(u8::is_ascii_digit).count();
        out.push_str(&format!("local str \"m{:0>2}\"", &rest[..digits]));
        rest = &rest[digits..];
    }
    out.push_str(rest);
    out
}

/// segmasker's ranges for three proteins of `mask_prot.faa`, two on one of them.
const PROTEIN_MASKS: &str = r#"Blast-db-mask-info ::= {
  algo-id 0,
  algo-program seg,
  algo-options "window=12; locut=2.2; hicut=2.5",
  masks {
    masks {
      int { from 10, to 40, id local str "p01" },
      packed-int {
        { from 0, to 20, id local str "p03" },
        { from 100, to 130, id local str "p03" }
      },
      int { from 90, to 122, id local str "p05" }
    },
    more FALSE
  }
}
"#;

/// Writes the database of `fasta` with the mask files `masks` (name, text) and the options
/// `more` (with `-parse_seqids`: a copy has its Seq-ids parsed), copies it with `-input_type
/// blastdb` and the same options, and compares every file of the two.
fn assert_copy_is_the_original(
    label: &str,
    fasta: &str,
    dbtype: &str,
    masks: &[(&str, String)],
    more: &[&str],
) {
    let tmp = TempDir::new(&format!("raw-{label}"));
    let fasta = seqdb_fixtures().join(fasta).to_string_lossy().into_owned();
    let masks: Vec<String> = masks
        .iter()
        .map(|(name, text)| {
            std::fs::write(tmp.0.join(name), text).unwrap();
            tmp.base(name)
        })
        .collect();
    let (source, copy) = (tmp.base("src/db"), tmp.base("dst/db"));
    let mut words = vec![
        "-in",
        &fasta,
        "-dbtype",
        dbtype,
        "-out",
        &source,
        "-parse_seqids",
    ];
    let joined = masks.join(",");
    words.extend(["-mask_data", &joined]);
    words.extend(more);
    let (status, _, err) = makeblastdb(&words);
    assert_eq!((status, err.as_str()), (0, ""), "{label}: the FASTA build");

    let mut words = vec![
        "-in",
        &source,
        "-input_type",
        "blastdb",
        "-dbtype",
        dbtype,
        "-out",
        &copy,
    ];
    words.extend(more);
    let (status, log, err) = makeblastdb(&words);
    assert_eq!((status, err.as_str()), (0, ""), "{label}: the copy");
    assert!(log.contains("\nKeep Leaf Taxids: T\n"), "{label}: {log}");
    assert!(
        log.contains("\nAdding sequences from raw db source; added "),
        "{label}: {log}"
    );
    // The title is the input's, which is the name of the FASTA file.
    assert!(
        log.contains(&format!("\nNew DB title:  {fasta}\n")),
        "{label}: {log}"
    );

    let (source_dir, copy_dir) = (tmp.0.join("src"), tmp.0.join("dst"));
    let files = files_of(&source_dir);
    assert_eq!(files, files_of(&copy_dir), "{label}: files");
    let ext = if dbtype == "prot" { 'p' } else { 'n' };
    assert!(
        files
            .iter()
            .any(|file| file.ends_with(&format!(".{ext}aa"))),
        "{label}: a mask column in {files:?}"
    );
    for file in &files {
        let read = |dir: &Path| {
            crate::writedb::tests::mask_column_date(file, std::fs::read(dir.join(file)).unwrap())
        };
        let (original, copied) = (read(&source_dir), read(&copy_dir));
        assert!(
            original == copied,
            "{label}: {file}: {} bytes and {} bytes, first difference at byte {:?}",
            original.len(),
            copied.len(),
            original.iter().zip(&copied).position(|(a, b)| a != b)
        );
    }
}

/// A nucleotide database with ambiguities, two masking algorithms, a hash index and several
/// volumes (version 5 and 4), and a protein database with one algorithm, are copied file for
/// file: sequences, ambiguities, headers, masks, the column metadata, the indices, the volumes'
/// cut, the alias and metadata files and the LMDB files.
#[test]
#[cfg_attr(
    target_os = "wasi",
    ignore = "runs makeblastdb on files named from the working directory (tests/...): a WASI module has none, the runner preopens /"
)]
fn a_copied_database_is_the_database_written_from_its_fasta() {
    let nucl_masks = [
        ("dust.asnt", masks_by_name("mask_writer/dust.asnt")),
        ("wm.asnt", masks_by_name("mask_wm.asnt")),
    ];
    assert_copy_is_the_original(
        "n5",
        "mask_nucl.fna",
        "nucl",
        &nucl_masks,
        &["-hash_index", "-max_file_sz", "400"],
    );
    assert_copy_is_the_original(
        "n4",
        "mask_nucl.fna",
        "nucl",
        &nucl_masks,
        &[
            "-hash_index",
            "-max_file_sz",
            "400",
            "-blastdb_version",
            "4",
        ],
    );
    assert_copy_is_the_original(
        "p5",
        "mask_prot.faa",
        "prot",
        &[("seg.asnt", PROTEIN_MASKS.to_string())],
        &["-hash_index"],
    );
}

/// The probes' head of the log of `dir`'s database `db` titled `X`.
fn log_head(dir: &str) -> String {
    format!(
        "\n\nBuilding a new DB, current time: 11/21/2026 09:00:00\nNew DB name:   {dir}/db\n\
         New DB title:  X\nSequence type: Nucleotide\nKeep MBits: T\n\
         Maximum file size: 3000000000B\n"
    )
}

/// NCBI's error for the nucleotide database `name` that is not found in the search path, whose
/// directories start with the current one (LOSAT's text of it).
fn not_found(name: &str) -> String {
    let cwd = std::env::current_dir().unwrap();
    let ncbi = format!(
        "BLAST Database error: No alias or index file found for nucleotide database [{name}] in \
         search path [{}::]. Please verify the spelling of the BLAST database and its molecule \
         type.\n",
        cwd.display()
    );
    String::from_utf8(crate::wording::losat_bytes(ncbi.as_bytes()).into_owned()).unwrap()
}

/// The errors of `-input_type blastdb` before the copy: an output named as the input (`-out`
/// left out or the same), several names without `-title`, a database that cannot be opened for
/// its title (exit status 2, before the log), and one that cannot be opened with `-title` (each
/// name is tried, posted and kept, then the open of all fails; no file is written).
#[test]
#[cfg_attr(
    target_os = "wasi",
    ignore = "runs makeblastdb on files named from the working directory (tests/...): a WASI module has none, the runner preopens /"
)]
fn the_errors_before_the_copy_are_makeblastdb_s() {
    let tmp = TempDir::new("raw-errors");
    let dir = tmp.0.to_string_lossy().into_owned();
    let out = tmp.base("db");
    let same_name = "BLAST options error: Cannot create a BLAST database from an existing one \
                     without changing the output name, please provide a (different) database \
                     name using -out\n";
    let words = ["-input_type", "blastdb", "-dbtype", "nucl"];
    let run = |more: &[&str]| {
        let mut all = words.to_vec();
        all.extend(more);
        makeblastdb(&all)
    };
    // e_noout, e_sameout
    assert_eq!(
        run(&["-in", "x"]),
        (1, String::new(), same_name.to_string())
    );
    assert_eq!(
        run(&["-in", "x", "-out", "x", "-title", "T"]),
        (1, String::new(), same_name.to_string())
    );
    // e_two_notitle
    assert_eq!(
        run(&["-in", "a b", "-out", &out]),
        (
            1,
            String::new(),
            "BLAST options error: Please provide a title using -title\n".to_string()
        )
    );
    // e_missing
    assert_eq!(
        run(&["-in", "nonexist", "-out", &out]),
        (2, String::new(), not_found("nonexist"))
    );
    // e_missing_t, e_dash
    for name in ["nonexist", "-"] {
        assert_eq!(
            run(&["-in", name, "-out", &out, "-title", "X"]),
            (
                2,
                log_head(&dir),
                format!(
                    "Error: [makeblastdb] Unable to open input {name} as BLAST db\n{}",
                    not_found(name)
                )
            ),
            "{name}"
        );
    }
    assert_eq!(files_of(&tmp.0), Vec::<String>::new());
    // With no name at all NCBI reads the title of the first database of none (makeblastdb
    // 2.17.0: a segmentation fault, probe d3c2_empty/b).
    let (status, _, err) = run(&["-in", "", "-out", &out]);
    assert_eq!(status, 1);
    assert!(
        err.starts_with("losat makeblastdb: -in names no database"),
        "{err}"
    );
    // With a title, the empty list is NCBI's error after the head of the log (probe
    // d3c2_empty/a).
    let (status, log, err) = run(&["-in", "", "-out", &out, "-title", "X"]);
    assert_eq!(
        (status, log, err.as_str()),
        (
            1,
            log_head(&dir),
            "BLAST options error: No valid input FASTA file or BLAST db is found.\n"
        )
    );
    assert_eq!(files_of(&tmp.0), Vec::<String>::new());
}
