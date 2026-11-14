//! Release-facing CLI smoke tests.

use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

fn temp_path(stem: &str, extension: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system time before Unix epoch")
        .as_nanos();
    std::env::temp_dir().join(format!(
        "losat_cli_{stem}_{}_{}.{}",
        std::process::id(),
        nanos,
        extension
    ))
}

fn losat_binary() -> PathBuf {
    if let Some(path) = option_env!("CARGO_BIN_EXE_losat") {
        return PathBuf::from(path);
    }

    let mut path = std::env::current_exe().expect("current test executable");
    path.pop();
    if path.file_name().is_some_and(|name| name == "deps") {
        path.pop();
    }
    path.push(format!("losat{}", std::env::consts::EXE_SUFFIX));
    path
}

fn clean_losat_command() -> Command {
    let mut command = Command::new(losat_binary());
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") {
            command.env_remove(name);
        }
    }
    command
}

fn run_losat(args: &[&str]) -> Output {
    clean_losat_command()
        .args(args)
        .output()
        .expect("run LOSAT")
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blastn_args.cpp:48-53
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blastp_args.cpp:47-52
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/tblastx_args.cpp:47-52
// ```c
// m_ClientId = string(kProgram) + " " + CBlastVersion().Print();
// m_ClientId = kProgram + " " + CBlastVersion().Print();
// ```
#[test]
fn cli_help_and_version_are_available_without_inputs() {
    for args in [
        vec!["--help"],
        vec!["blastn", "--help"],
        vec!["tblastx", "--help"],
        vec!["blastp", "--help"],
    ] {
        let output = run_losat(&args);
        assert!(
            output.status.success(),
            "{args:?} failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            output.stderr.is_empty(),
            "{args:?} should not write help to stderr: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(
            String::from_utf8_lossy(&output.stdout).contains("USAGE\n"),
            "{args:?} help should open with NCBI's USAGE heading"
        );
    }

    let output = run_losat(&["--version"]);
    assert!(
        output.status.success(),
        "--version failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&output.stdout).trim(),
        format!("losat {}", env!("CARGO_PKG_VERSION"))
    );
    assert!(
        output.stderr.is_empty(),
        "--version should not write stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

// No NCBI counterpart: LOSAT_STARTUP_TRACE is a LOSAT diagnostic variable; it does not change any
// value NCBI computes.
// CI's `--all-features` run checks the feature side of this test (the trace line, no notice); a
// default `cargo test` checks the other (the notice line, no trace). `ci_fast_regressions.py`
// checks the notice on the default release binary.
#[test]
fn startup_trace_needs_the_diagnostics_feature() {
    let plain = clean_losat_command()
        .arg("--version")
        .output()
        .expect("run LOSAT");
    let traced = clean_losat_command()
        .arg("--version")
        .env("LOSAT_STARTUP_TRACE", "1")
        .output()
        .expect("run LOSAT");
    assert!(plain.status.success() && traced.status.success());
    assert_eq!(plain.stdout, traced.stdout);
    assert!(
        plain.stderr.is_empty(),
        "--version should not write stderr: {}",
        String::from_utf8_lossy(&plain.stderr)
    );
    let expected = if cfg!(feature = "diagnostics") {
        "[startup] enter main\n"
    } else {
        "LOSAT_STARTUP_TRACE is ignored: this losat was built without the diagnostics feature \
         (docs/agents/diagnostics.md, \"Diagnostic build\")\n"
    };
    assert_eq!(String::from_utf8_lossy(&traced.stderr), expected);
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3631-3639
// ```c
//     NON_CONST_ITERATE(TBlastCmdLineArgs, arg, m_Args) {
//         (*arg)->ExtractAlgorithmOptions(args, opts);
//     }
//
//     m_IsUngapped = !opts.GetGappedMode();
//     try { retval->Validate(); }
// ```
// NCBI's handlers read the subjects and open the query before the options are checked, so a
// missing file is reported first; LOSAT's rejection of an option that NCBI runs comes where
// NCBI checks the options.
#[test]
fn unsupported_blastp_option_fails_after_ncbis_file_checks() {
    let missing_query = temp_path("missing_query", "faa");
    let missing_subject = temp_path("missing_subject", "faa");
    let output = clean_losat_command()
        .args(["blastp", "-query"])
        .arg(&missing_query)
        .arg("-subject")
        .arg(&missing_subject)
        .args(["-word_size", "4"])
        .output()
        .expect("run unsupported blastp");
    assert!(!output.status.success(), "a missing subject should fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        !stderr.contains("not supported by LOSAT's BLASTP"),
        "the subject is read before the options are checked: {stderr}"
    );

    let query = temp_path("word4_query", "faa");
    let subject = temp_path("word4_subject", "faa");
    std::fs::write(&query, ">q\nMKTAYIAKQRQISFVKSHFSRQ\n").unwrap();
    std::fs::write(&subject, ">s\nMKTAYIAKQRQISFVKSHFSRQ\n").unwrap();
    let output = clean_losat_command()
        .args(["blastp", "-query"])
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .args(["-word_size", "4"])
        .output()
        .expect("run unsupported blastp");
    let _ = std::fs::remove_file(&query);
    let _ = std::fs::remove_file(&subject);
    assert!(!output.status.success(), "unsupported blastp should fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("-word_size 4 (the lookup table) is not supported by LOSAT's BLASTP"),
        "unsupported error should name program, option, and reason: {stderr}"
    );
    assert!(output.stdout.is_empty(), "nothing is written: {stderr}");
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_input_aux.cpp:242-246
// ```c
// CRef<CBlastFastaInputSource> fasta(new CBlastFastaInputSource(in, iconfig));
// CRef<CBlastInput> input(new CBlastInput(fasta));
// CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
// sequences = input->GetAllSeqs(*scope);
// ```
#[test]
fn missing_input_file_reports_ncbi_error() {
    let missing_query = temp_path("missing_query", "fa");
    let missing_subject = temp_path("missing_subject", "fa");
    let output = clean_losat_command()
        .args(["blastn", "-outfmt", "6", "-query"])
        .arg(&missing_query)
        .arg("-subject")
        .arg(&missing_subject)
        .output()
        .expect("run blastn with missing input");

    assert!(!output.status.success(), "missing input should fail");
    let stderr = String::from_utf8_lossy(&output.stderr);
    // NCBI opens the subject first (blastn_args.cpp:63-70), with NCBI's message
    // (ncbiargs.cpp:615-619).
    assert_eq!(output.status.code(), Some(1));
    assert_eq!(
        stderr,
        format!(
            "Command line argument error: Argument \"subject\". File is not accessible:  `{}'\n",
            missing_subject.display()
        ),
        "missing input error should be NCBI's: {stderr}"
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:751-758
// ```c++
//     if (bad >= good / 3  &&
//         (len_to_check > 3  ||  good == 0  ||  bad > good))
//     {
//         FASTA_ERROR( LineNumber(),
//             "CFastaReader: Near line " << LineNumber()
//             << ", there's a line that doesn't look like plausible data, "
//             "but it's not marked as defline or comment.",
//             CObjReaderParseException::eFormat);
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:181-184
// ```c++
//     catch (const CObjReaderParseException& e) {                             \
//         LOG_POST(Error << "BLAST query error: " << e.GetMsg());             \
//         exit_code = BLAST_INPUT_ERROR;                                      \
//     }                                                                       \
// ```
// A comment line before the first defline is skipped; a first data line that does not
// look like sequence data stops BLASTN with NCBI's reader error.
#[test]
fn malformed_fasta_reports_parse_error() {
    let query = temp_path("malformed_query", "fa");
    let subject = temp_path("valid_subject", "fa");
    fs::write(&query, "!a comment line\n@@@@ACGT\n").expect("write malformed FASTA");
    fs::write(&subject, ">s\nACGTACGTACGT\n").expect("write subject FASTA");

    let output = clean_losat_command()
        .args(["blastn", "-outfmt", "6", "-query"])
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .output()
        .expect("run blastn with malformed FASTA");

    let _ = fs::remove_file(&query);
    let _ = fs::remove_file(&subject);

    assert_eq!(output.status.code(), Some(1), "malformed FASTA should fail");
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "BLAST query error: CFastaReader: Near line 2, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.\n"
    );
}

/// Gzip members of `members`, each with bgzip's extra field, concatenated as `bgzip` and
/// `cat a.gz b.gz` write them.
fn gzip_members(members: &[&[u8]]) -> Vec<u8> {
    use std::io::Write;
    let mut bytes = Vec::new();
    for member in members {
        let mut encoder = flate2::GzBuilder::new()
            .extra(b"BC\x02\x00\x00\x00".to_vec())
            .write(&mut bytes, flate2::Compression::default());
        encoder.write_all(member).unwrap();
        encoder.finish().unwrap();
    }
    bytes
}

/// Bzip2 streams, concatenated as pbzip2 and `cat` write them.
fn bzip2_streams(streams: &[&[u8]]) -> Vec<u8> {
    use std::io::Write;
    let mut bytes = Vec::new();
    for stream in streams {
        let mut encoder = bzip2::write::BzEncoder::new(&mut bytes, bzip2::Compression::default());
        encoder.write_all(stream).unwrap();
        encoder.finish().unwrap();
    }
    bytes
}

/// Zstd frames of raw blocks (LOSAT only decodes zstd), each after a skippable frame as pzstd
/// writes one, concatenated.
fn zstd_frames(frames: &[&[u8]]) -> Vec<u8> {
    let mut bytes = Vec::new();
    for frame in frames {
        bytes.extend_from_slice(&0x184D_2A50u32.to_le_bytes());
        bytes.extend_from_slice(&4u32.to_le_bytes());
        bytes.extend_from_slice(b"skip");
        bytes.extend_from_slice(&0xFD2F_B528u32.to_le_bytes());
        bytes.push(0x80 | 0x20);
        bytes.extend_from_slice(&(frame.len() as u32).to_le_bytes());
        let mut blocks = frame.chunks(128 * 1024).peekable();
        if blocks.peek().is_none() {
            bytes.extend_from_slice(&[1, 0, 0]);
        }
        while let Some(block) = blocks.next() {
            let last = u32::from(blocks.peek().is_none());
            let header = last | ((block.len() as u32) << 3);
            bytes.extend_from_slice(&header.to_le_bytes()[..3]);
            bytes.extend_from_slice(block);
        }
    }
    bytes
}

/// The three compressions LOSAT reads by name: a suffix (upper or lower case, as the name rule
/// ignores case), the format's name, and an encoder of concatenated parts.
const COMPRESSIONS: [(&str, &str, &str, fn(&[&[u8]]) -> Vec<u8>); 3] = [
    (".gz", ".GZ", "gzip", gzip_members),
    (".bz2", ".Bz2", "bzip2", bzip2_streams),
    (".zst", ".ZST", "zstd", zstd_frames),
];

const GZ_QUERY: &str = ">q1 first\nACGTTGCAAGGCTTAACCGTATGCGATCCTAGGATCCAATGGTACCTTAGCAGTCGATCGACGTTGCAAGGCTTAACCGTATGCG\n>q2 second\nGGCATGCATGCATTTAGGCCAAGTACGATCGATCGACGTTGCAAGGCTTAACCGTATGCGATCCTAGGATCCAATGG\n";
const GZ_SUBJECT: &str = ">s1 the subject\nTTTTACGTTGCAAGGCTTAACCGTATGCGATCCTAGGATCCAATGGTACCTTAGCAGTCGATCGACGTTGCAAGGCTTAACCGTATGCGAAAAGGCATGCATGCATTTAGGCCAAGTACGATCGATCGACGTTGCAAGGCTTAACCGTATGCGATCCTAGGATCCAATGGCCCC\n";

fn blastn_output(query: &std::path::Path, subject: &std::path::Path, outfmt: &str) -> Output {
    clean_losat_command()
        .args(["blastn", "-outfmt", outfmt, "-query"])
        .arg(query)
        .arg("-subject")
        .arg(subject)
        .output()
        .expect("run blastn")
}

// No NCBI counterpart (line GZ, DECISIONS 1, 3 and 8; `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`
// exception 8): a `-query` or `-subject` named `.gz`, `.bz2` or `.zst` reads as the plain file,
// in any case of the suffix; the name the user gave is the one the reports show.
#[test]
fn compressed_query_and_subject_read_as_the_plain_files() {
    let query = temp_path("plain_query", "fa");
    let subject = temp_path("plain_subject", "fa");
    fs::write(&query, GZ_QUERY).unwrap();
    fs::write(&subject, GZ_SUBJECT).unwrap();
    let plain = blastn_output(&query, &subject, "6");
    let plain_7 = blastn_output(&query, &subject, "7");
    assert_eq!(
        plain.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&plain.stderr)
    );
    assert!(!plain.stdout.is_empty(), "the plain run has rows");
    let plain_7 = String::from_utf8_lossy(&plain_7.stdout).into_owned();
    let input_line = |text: &str, path: &std::path::Path| {
        let line = format!(
            "# Database: User specified sequence set (Input: {})",
            path.display()
        );
        assert!(text.contains(&line), "{text}");
        text.replace(
            &line,
            "# Database: User specified sequence set (Input: SUBJECT)",
        )
    };
    let plain_7 = input_line(&plain_7, &subject);
    for (suffix, upper, format, encode) in COMPRESSIONS {
        let query_c = temp_path("compressed_query", &format!("fa{suffix}"));
        let subject_c = temp_path("compressed_subject", &format!("fna{upper}"));
        fs::write(&query_c, encode(&[GZ_QUERY.as_bytes()])).unwrap();
        fs::write(&subject_c, encode(&[GZ_SUBJECT.as_bytes()])).unwrap();
        let inflated = blastn_output(&query_c, &subject_c, "6");
        let inflated_7 = blastn_output(&query_c, &subject_c, "7");
        let _ = fs::remove_file(&query_c);
        let _ = fs::remove_file(&subject_c);
        assert_eq!(
            inflated.status.code(),
            Some(0),
            "{format}: {}",
            String::from_utf8_lossy(&inflated.stderr)
        );
        assert_eq!(inflated.stdout, plain.stdout, "{format}");
        assert_eq!(inflated.stderr, plain.stderr, "{format}");
        let inflated_7 = String::from_utf8_lossy(&inflated_7.stdout).into_owned();
        assert_eq!(input_line(&inflated_7, &subject_c), plain_7, "{format}");
        // A gzip trailer is not checked, as makeblastdb does not check it
        // (`makeblastdb_reads_compressed_inputs_as_makeblastdb_does`): a wrong CRC-32 reads as
        // the intact file.
        if format == "gzip" {
            let mut bytes = encode(&[GZ_QUERY.as_bytes()]);
            let crc_at = bytes.len() - 8;
            bytes[crc_at] ^= 0xff;
            let badcrc = temp_path("compressed_query_badcrc", &format!("fa{suffix}"));
            fs::write(&badcrc, bytes).unwrap();
            let read = blastn_output(&badcrc, &subject, "6");
            let _ = fs::remove_file(&badcrc);
            assert_eq!(read.status.code(), Some(0), "{format}: bad CRC");
            assert_eq!(read.stdout, plain.stdout, "{format}: bad CRC");
            assert_eq!(read.stderr, plain.stderr, "{format}: bad CRC");
        }
    }
    let _ = fs::remove_file(&query);
    let _ = fs::remove_file(&subject);
}

// No NCBI counterpart (line GZ, DECISIONS 4 and 9): the members, streams and frames of a
// concatenated file (bgzip's blocks with their extra field and empty end block, pbzip2's
// streams, pzstd's skippable and data frames) are all read; a reader error inside the inflated
// file keeps the line number of the plain file.
#[test]
fn concatenated_parts_are_read_whole_with_the_plain_line_numbers() {
    let subject = temp_path("parts_subject", "fa");
    fs::write(&subject, GZ_SUBJECT).unwrap();
    let (first, second) = GZ_QUERY.split_at(GZ_QUERY.find(">q2").unwrap());
    for (suffix, _, format, encode) in COMPRESSIONS {
        let query_c = temp_path("parts_query", &format!("fa{suffix}"));
        fs::write(
            &query_c,
            encode(&[first.as_bytes(), second.as_bytes(), b""]),
        )
        .unwrap();
        let both = blastn_output(&query_c, &subject, "6");
        let bad_c = temp_path("parts_bad", &format!("fa{suffix}"));
        fs::write(&bad_c, encode(&[b"!a comment line\n", b"@@@@ACGT\n"])).unwrap();
        let bad = blastn_output(&bad_c, &subject, "6");
        let _ = fs::remove_file(&query_c);
        let _ = fs::remove_file(&bad_c);
        assert_eq!(
            both.status.code(),
            Some(0),
            "{format}: {}",
            String::from_utf8_lossy(&both.stderr)
        );
        let rows = String::from_utf8_lossy(&both.stdout).into_owned();
        assert!(
            rows.lines().any(|row| row.starts_with("q1\t")),
            "{format}: {rows}"
        );
        assert!(
            rows.lines().any(|row| row.starts_with("q2\t")),
            "{format}: {rows}"
        );
        assert_eq!(bad.status.code(), Some(1), "{format}");
        assert_eq!(
            String::from_utf8_lossy(&bad.stderr),
            "BLAST query error: CFastaReader: Near line 2, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.\n",
            "{format}"
        );
    }
    let _ = fs::remove_file(&subject);
}

// No NCBI counterpart (line GZ, DECISIONS 5): a compressed name whose bytes are not that format
// (for `.gz`: that makeblastdb's decompressor fails on) is refused before the search, with exit
// status 1, and `-` is read as it comes (never inflated).
#[test]
fn a_compressed_name_without_that_data_is_refused_before_the_search() {
    let subject = temp_path("text_subject", "fa");
    fs::write(&subject, GZ_SUBJECT).unwrap();
    for (suffix, _, format, _) in COMPRESSIONS {
        let query_c = temp_path("text_query", &format!("fa{suffix}"));
        fs::write(&query_c, GZ_QUERY).unwrap();
        let refused = blastn_output(&query_c, &subject, "6");
        let _ = fs::remove_file(&query_c);
        assert_eq!(refused.status.code(), Some(1), "{format}");
        assert_eq!(refused.stdout, b"", "{format}");
        let stderr = String::from_utf8_lossy(&refused.stderr).into_owned();
        assert!(stderr.starts_with("Error: -query: '"), "{format}: {stderr}");
        assert!(
            stderr.contains(&format!("' is named {suffix} but is not {format} data (")),
            "{format}: {stderr}"
        );
        assert!(
            stderr.trim_end().ends_with(&format!(
                "); LOSAT reads a file named {suffix} through {format}"
            )),
            "{format}: {stderr}"
        );
    }
    let _ = fs::remove_file(&subject);
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1037-1050 (an `-in` name that
// ends in `.gz`, `.bz2` or `.zst`, in any case, is read through a decompressing stream; the title
// defaults to the name as given). The database built from the inflated file is the one built
// from the plain file. Concatenated streams and frames are all read (makeblastdb reads only the
// first bzip2 stream or zstd frame: `PD-LOSAT-NCBI-DEFECTS`, line GZ DECISIONS 9); text under a
// compressed name stops the build by name (DECISIONS 5).
#[test]
fn makeblastdb_reads_compressed_inputs_as_makeblastdb_does() {
    let dir = temp_path("makedb_compressed", "dir");
    fs::create_dir_all(&dir).unwrap();
    let plain = dir.join("s.fna");
    let query = dir.join("q.fna");
    fs::write(&plain, GZ_SUBJECT).unwrap();
    fs::write(&query, GZ_QUERY).unwrap();
    let build = |input: &std::path::Path, out: &str| {
        clean_losat_command()
            .current_dir(&dir)
            .args(["makeblastdb", "-dbtype", "nucl", "-in"])
            .arg(input.file_name().unwrap())
            .args(["-out", out])
            .output()
            .expect("run makeblastdb")
    };
    let search = |db: &str| {
        clean_losat_command()
            .current_dir(&dir)
            .args(["blastn", "-outfmt", "6", "-query", "q.fna", "-db", db])
            .output()
            .expect("run blastn -db")
    };
    let from_plain = build(&plain, "db_plain");
    assert_eq!(
        from_plain.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&from_plain.stderr)
    );
    let plain_rows = search("db_plain");
    assert_eq!(
        plain_rows.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&plain_rows.stderr)
    );
    assert!(!plain_rows.stdout.is_empty());
    let (first, second) = GZ_SUBJECT.split_at(GZ_SUBJECT.find("AAAAGGCATG").unwrap() + 4);
    for (suffix, upper, format, encode) in COMPRESSIONS {
        let name = format!("s.fna{upper}");
        let compressed = dir.join(&name);
        fs::write(&compressed, encode(&[GZ_SUBJECT.as_bytes()])).unwrap();
        let out = format!("db{suffix}");
        let from_compressed = build(&compressed, &out);
        assert_eq!(
            from_compressed.status.code(),
            Some(0),
            "{format}: {}",
            String::from_utf8_lossy(&from_compressed.stderr)
        );
        assert_eq!(from_compressed.stderr, from_plain.stderr, "{format}");
        let log = String::from_utf8_lossy(&from_compressed.stdout).into_owned();
        assert!(
            log.contains(&format!("New DB title:  {name}\n")),
            "{format}: {log}"
        );
        assert!(log.contains("added 1 sequences in"), "{format}: {log}");
        for extension in ["nsq", "nhr"] {
            assert_eq!(
                fs::read(dir.join(format!("{out}.{extension}"))).unwrap(),
                fs::read(dir.join(format!("db_plain.{extension}"))).unwrap(),
                "{format}: {extension}"
            );
        }
        let rows = search(&out);
        assert_eq!(rows.stdout, plain_rows.stdout, "{format}");

        // Two parts (a record split over them) read as one input: "added 1 sequences".
        let parts = dir.join(format!("parts.fna{suffix}"));
        fs::write(&parts, encode(&[first.as_bytes(), second.as_bytes()])).unwrap();
        let from_parts = build(&parts, &format!("db_parts{suffix}"));
        assert_eq!(
            from_parts.status.code(),
            Some(0),
            "{format}: {}",
            String::from_utf8_lossy(&from_parts.stderr)
        );
        assert_eq!(
            fs::read(dir.join(format!("db_parts{suffix}.nsq"))).unwrap(),
            fs::read(dir.join("db_plain.nsq")).unwrap(),
            "{format}: the parts"
        );

        let text = dir.join(format!("text.fna{suffix}"));
        fs::write(&text, GZ_SUBJECT).unwrap();
        let refused = build(&text, &format!("db_text{suffix}"));
        assert_eq!(refused.status.code(), Some(1), "{format}");
        let stderr = String::from_utf8_lossy(&refused.stderr).into_owned();
        assert!(
            stderr.contains(&format!("losat makeblastdb: -in: 'text.fna{suffix}' is named {suffix} but is not {format} data (")),
            "{format}: {stderr}"
        );
        assert!(
            !dir.join(format!("db_text{suffix}.nsq")).exists(),
            "{format}: the failed build leaves no files"
        );

        // NCBI reference (598d8ae6): c++/src/util/compress/api/zlib.cpp:1513-1526
        // ```c++
        //             if ((errcode == Z_STREAM_END)  &&  m_IsGZ) {
        //                 // Skip .gz file footer (8 bytes)
        //                 if (STREAM->avail_in < 8) {
        //                     m_SkipInput = 8 - STREAM->avail_in;
        //                     STREAM->avail_in = 0;
        //                 } else {
        //                     STREAM->avail_in -= 8;
        //                 }
        // ```
        // makeblastdb skips a gzip member's trailer without checking it: a wrong CRC-32, or a
        // trailer cut short, builds the database of the intact file (makeblastdb 2.17.0,
        // 2026-11-08: exit 0, "added 1 sequences", no message).
        if format == "gzip" {
            let intact = encode(&[GZ_SUBJECT.as_bytes()]);
            let mut badcrc = intact.clone();
            let crc_at = intact.len() - 8;
            badcrc[crc_at] ^= 0xff;
            let cut = intact[..intact.len() - 4].to_vec();
            for (case, bytes) in [("badcrc", badcrc), ("cut4", cut)] {
                let input = dir.join(format!("{case}.fna.gz"));
                fs::write(&input, bytes).unwrap();
                let out = format!("db_{case}");
                let built = build(&input, &out);
                assert_eq!(
                    built.status.code(),
                    Some(0),
                    "{case}: {}",
                    String::from_utf8_lossy(&built.stderr)
                );
                assert_eq!(built.stderr, from_plain.stderr, "{case}");
                for extension in ["nsq", "nhr"] {
                    assert_eq!(
                        fs::read(dir.join(format!("{out}.{extension}"))).unwrap(),
                        fs::read(dir.join(format!("db_plain.{extension}"))).unwrap(),
                        "{case}: {extension}"
                    );
                }
            }
        }
    }
    let _ = fs::remove_dir_all(&dir);
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_input.cpp:140-152
// ```c++
//     while (size_read < GetBatchSize()) {
//
//         if (End())
//             break;
//
//         CRef<CBlastSearchQuery> q;
//         try { q.Reset(m_Source->GetNextSequence(scope)); }
//         catch (const CObjReaderParseException& e) {
//             if (e.GetErrCode() == CObjReaderParseException::eEOF) {
//                 break;
//             }
//             throw;
// ```
// NCBI reads the queries one batch at a time: a reader error in a later batch comes after
// the reports of the batches before, without the epilog (oracle BI b_q_b2_bad_7), and in
// the first batch before any report.
#[test]
fn a_reader_error_in_a_later_query_batch_follows_the_earlier_reports() {
    let query = temp_path("later_batch_query", "fa");
    let subject = temp_path("later_batch_subject", "fa");
    let sequence = "ACGTTGCAAGGCTTAACCGTATGCGATCCTAGGATCCAATGGTACCTTAGCAGTCGATCG".repeat(5);
    let mut lines = String::from(">q1\n");
    for chunk in sequence.as_bytes().chunks(60) {
        lines.push_str(std::str::from_utf8(chunk).unwrap());
        lines.push('\n');
    }
    lines.push_str(">q2\n@@@@\n");
    fs::write(&query, &lines).expect("write query FASTA");
    fs::write(&subject, format!(">s1\n{sequence}\n")).expect("write subject FASTA");
    let run = |batch_size: Option<&str>| {
        let mut command = clean_losat_command();
        command.args(["blastn", "-outfmt", "7", "-query"]);
        command.arg(&query).arg("-subject").arg(&subject);
        match batch_size {
            Some(size) => command.env("BATCH_SIZE", size),
            None => command.env_remove("BATCH_SIZE"),
        };
        command.output().expect("run blastn")
    };
    let later = run(Some("300"));
    let first = run(None);
    let _ = fs::remove_file(&query);
    let _ = fs::remove_file(&subject);
    let error = "BLAST query error: CFastaReader: Near line 8, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.\n";
    for output in [&later, &first] {
        assert_eq!(output.status.code(), Some(1));
        assert_eq!(String::from_utf8_lossy(&output.stderr), error);
    }
    let report = String::from_utf8_lossy(&later.stdout);
    assert!(
        report.starts_with("# BLASTN 2.17.0+\n# Query: q1\n")
            && !report.contains("# BLAST processed"),
        "{report}"
    );
    assert!(first.stdout.is_empty(), "{:?}", first.stdout);
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3265-3273
// ```c
// #if _BLAST_DEBUG
// arg_desc.AddFlag("verbose", "Produce verbose output (show BLAST options)",
//                  true);
// #endif /* _BLAST_DEBUG */
// ```
#[test]
fn normal_tblastx_run_keeps_stderr_clean_without_debug_env() {
    let query = temp_path("tblastx_query", "fa");
    let subject = temp_path("tblastx_subject", "fa");
    let out = temp_path("tblastx_out", "txt");
    fs::write(&query, ">q\nATGATGATGATGATGATGATGATGATGATG\n").expect("write query FASTA");
    fs::write(&subject, ">s\nATGATGATGATGATGATGATGATGATGATG\n").expect("write subject FASTA");

    let output = clean_losat_command()
        .args(["tblastx", "-query"])
        .arg(&query)
        .arg("-subject")
        .arg(&subject)
        .args(["-seg", "no", "-outfmt", "6", "-num_threads", "1", "-out"])
        .arg(&out)
        .output()
        .expect("run tblastx");

    let output_text = fs::read_to_string(&out).unwrap_or_default();
    let _ = fs::remove_file(&query);
    let _ = fs::remove_file(&subject);
    let _ = fs::remove_file(&out);

    assert!(
        output.status.success(),
        "tblastx smoke failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        output.stderr.is_empty(),
        "normal tblastx run should keep stderr clean: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        !output_text.is_empty(),
        "tblastx smoke input should produce at least one hit"
    );
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:972-977
// ```c++
//             } else if ( NStr::strcmp(argv[i], s_ArgVersion) == 0 ) {
//                 delete[] v;
//                 // Print VERSION
//                 cout << GetFullVersion().Print( appname, CVersionAPI::fVersionInfo | CVersionAPI::fPackageShort );
//                 diag_context.DiscardMessages();
//                 return 0;
// ```
// NCBI reads `-version` before any other argument: anywhere before `--`, in the value place
// of another option, and next to arguments that are invalid.
#[test]
fn ncbi_version_is_two_lines_for_every_program() {
    let package = format!(
        " Package: LOSAT {}, output compatible with blast 2.17.0\n",
        env!("CARGO_PKG_VERSION")
    );
    for (name, printed) in [
        ("blastn", "blastn"),
        ("blastp", "blastp"),
        ("blastx", "blastx"),
        ("tblastn", "tblastn"),
        ("tblastx", "tblastx"),
        ("psiblast", "psiblast"),
        ("makeblastdb", "makeblastdb"),
    ] {
        let expected = format!("{printed}: 2.17.0+\n{package}");
        for words in [
            vec![name, "-version"],
            vec![name, "-out", "-version"],
            vec![name, "-no_such_option", "-version"],
            vec![name, "-evalue", "x", "-version", "-bogus"],
        ] {
            let output = run_losat(&words);
            assert!(output.status.success(), "{words:?}: {:?}", output.status);
            assert_eq!(
                String::from_utf8_lossy(&output.stdout),
                expected,
                "{words:?}"
            );
            assert!(output.stderr.is_empty(), "{words:?}");
        }
        // After `--` the word is a positional argument, which no program has.
        let output = run_losat(&[name, "--", "-version"]);
        assert!(!output.status.success(), "{name} -- -version");
    }
    // The top level `--version` is clap's.
    let output = run_losat(&["--version"]);
    assert_eq!(
        String::from_utf8_lossy(&output.stdout).trim(),
        format!("losat {}", env!("CARGO_PKG_VERSION"))
    );
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:76-77
// ```c++
// static const char* s_AutoHelp     = "h";
// static const char* s_AutoHelpFull = "help";
// ```
#[test]
fn short_help_h_succeeds_for_every_program() {
    for name in [
        "blastn",
        "blastp",
        "blastx",
        "tblastn",
        "tblastx",
        "psiblast",
        "makeblastdb",
    ] {
        let words = vec![name, "-h"];
        let output = run_losat(&words);
        assert!(output.status.success(), "{words:?}: {:?}", output.status);
        assert!(output.stderr.is_empty(), "{words:?}");
        let text = String::from_utf8_lossy(&output.stdout);
        assert!(text.starts_with("USAGE\n"), "{words:?}: {text}");
        assert!(text.contains("\nDESCRIPTION\n"), "{words:?}: {text}");
        assert!(!text.contains("error:"), "{words:?}: {text}");
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:217-219
// ```c++
//     arg_desc->SetUsageContext(GetArguments().GetProgramBasename(),
//                   "Application to create BLAST databases, version "
//                   + CBlastVersion().Print());
// ```
// `losat makeblastdb` writes the database; its old name `makedb` is an unknown command.
#[test]
fn makeblastdb_writes_a_database_and_makedb_is_unknown() {
    let input = temp_path("mkdb_in", "faa");
    fs::write(
        &input,
        ">a first\nMKTAYIAKQRQISFVKSHFSRQ\n>b second\nLEERLGLIEVQAPILSRVGDGTQDN\n",
    )
    .unwrap();
    let base = temp_path("makeblastdb", "db");
    let (input_text, base_text) = (input.to_str().unwrap(), base.to_str().unwrap());
    let output = run_losat(&[
        "makeblastdb",
        "-in",
        input_text,
        "-dbtype",
        "prot",
        "-out",
        base_text,
    ]);
    assert!(output.status.success(), "{:?}", output);
    for extension in ["phr", "psq", "pin"] {
        let path = format!("{base_text}.{extension}");
        assert!(!fs::read(&path).unwrap().is_empty(), ".{extension}");
        fs::remove_file(path).unwrap();
    }
    // `makedb` is not a command, with or without a valid command line.
    for words in [
        vec![
            "makedb", "-in", input_text, "-dbtype", "prot", "-out", base_text,
        ],
        vec!["makedb", "-help"],
        vec!["makedb", "-version"],
    ] {
        let output = run_losat(&words);
        assert!(!output.status.success(), "{words:?}");
        assert!(output.stdout.is_empty(), "{words:?}");
        assert!(
            !std::path::Path::new(&format!("{base_text}.psq")).exists(),
            "{words:?}"
        );
    }
    fs::remove_file(&input).unwrap();
    // `-help` works, and the top level lists makeblastdb and not makedb.
    let output = run_losat(&["makeblastdb", "-help"]);
    assert!(output.status.success());
    assert!(String::from_utf8_lossy(&output.stdout).contains("-dbtype <molecule_type>"));
    let top = String::from_utf8_lossy(&run_losat(&["-help"]).stdout).into_owned();
    assert!(
        top.contains("makeblastdb") && !top.contains("makedb"),
        "{top}"
    );
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1203-1209
// ```c++
//     bool parse_seqids = x_ShouldParseSeqIds();
//     bool hash_index = args["hash_index"];
//     bool use_gi_mask = args["gi_mask"];
//
//     CWriteDB::TIndexType indexing = CWriteDB::eNoIndex;
//     indexing |= (hash_index ? CWriteDB::eAddHash : 0);
//     indexing |= (parse_seqids ? CWriteDB::eFullIndex : 0);
// ```
// `-parse_seqids`, `-hash_index` and `-taxid_map` on the command line: the id files of version 5
// (`.pos`, `.pog`, the hash index) and of version 4 (`.psi`/`.psd`), the argument rules (exit
// status 2 from LOSAT's parser), a version 5 input makeblastdb aborts on (exit status 1), and
// `[BLAST] MAX_SEQID_LENGTH` refused with `-parse_seqids` (makedb-v5 DECISIONS 12). The bytes are
// checked against makeblastdb in `writedb::tests::parse_seqids_databases_are_makeblastdb_s`.
#[test]
fn makeblastdb_parse_seqids_hash_index_and_taxid_map() {
    let input = temp_path("mkdb_parse_in", "faa");
    fs::write(
        &input,
        ">ref|NP_000001.1| first\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|second b\nLEERLGLIEVQAPILSRVGDGTQDN\n",
    )
    .unwrap();
    let map = temp_path("mkdb_parse_map", "txt");
    fs::write(&map, "NP_000001.1 9606\nsecond 10090\n").unwrap();
    let base = temp_path("mkdb_parse", "db");
    let (input_text, map_text, base_text) = (
        input.to_str().unwrap(),
        map.to_str().unwrap(),
        base.to_str().unwrap(),
    );
    let remove = |extensions: &[&str]| {
        for extension in extensions {
            let _ = fs::remove_file(format!("{base_text}.{extension}"));
        }
    };
    const V5: [&str; 12] = [
        "pdb", "phd", "phi", "phr", "pin", "pjs", "pog", "pos", "pot", "psq", "ptf", "pto",
    ];
    let output = run_losat(&[
        "makeblastdb",
        "-in",
        input_text,
        "-dbtype",
        "prot",
        "-out",
        base_text,
        "-parse_seqids",
        "-hash_index",
        "-taxid_map",
        map_text,
    ]);
    assert!(output.status.success(), "{output:?}");
    for extension in V5 {
        let path = format!("{base_text}.{extension}");
        assert!(!fs::read(&path).unwrap().is_empty(), ".{extension}");
    }
    let metadata = fs::read_to_string(format!("{base_text}.pjs")).unwrap();
    assert!(metadata.contains("\"number-of-taxids\": 2"), "{metadata}");
    remove(&V5);
    let output = run_losat(&[
        "makeblastdb",
        "-in",
        input_text,
        "-dbtype",
        "prot",
        "-out",
        base_text,
        "-parse_seqids",
        "-blastdb_version",
        "4",
    ]);
    assert!(output.status.success(), "{output:?}");
    for extension in ["psi", "psd", "pog"] {
        assert!(
            std::path::Path::new(&format!("{base_text}.{extension}")).exists(),
            ".{extension}"
        );
    }
    remove(&["phr", "pin", "pjs", "pog", "psd", "psi", "psq"]);
    // `-taxid_map` requires `-parse_seqids` and excludes `-taxid`: LOSAT's parser, exit 2.
    for words in [
        vec!["-taxid_map", map_text],
        vec!["-parse_seqids", "-taxid_map", map_text, "-taxid", "1"],
    ] {
        let mut argv = vec![
            "makeblastdb",
            "-in",
            input_text,
            "-dbtype",
            "prot",
            "-out",
            base_text,
        ];
        argv.extend(words.iter().copied());
        let output = run_losat(&argv);
        assert_eq!(output.status.code(), Some(2), "{words:?}: {output:?}");
        assert!(!std::path::Path::new(&format!("{base_text}.psq")).exists());
    }
    // A GI-only sequence between sequences with other ids: makeblastdb 2.17.0 aborts on version 5.
    fs::write(
        &input,
        ">ref|NP_1.1| a\nMKTAYIAKQR\n>gi|42 b\nMKTAYIAKQR\n>ref|NP_2.1| c\nMKTAYIAKQR\n",
    )
    .unwrap();
    let output = run_losat(&[
        "makeblastdb",
        "-in",
        input_text,
        "-dbtype",
        "prot",
        "-out",
        base_text,
        "-parse_seqids",
    ]);
    assert_eq!(output.status.code(), Some(1), "{output:?}");
    assert!(String::from_utf8_lossy(&output.stderr).starts_with(
        "BLAST Database creation error: Input id list not in ascending oid order: the sequence gi|42 (OID 1)"
    ));
    assert!(!std::path::Path::new(&format!("{base_text}.psq")).exists());
    remove(&V5);
    // `[BLAST] MAX_SEQID_LENGTH` is refused with `-parse_seqids` only.
    let refused = clean_losat_command()
        .env("NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH", "100")
        .args([
            "makeblastdb",
            "-in",
            input_text,
            "-dbtype",
            "prot",
            "-out",
            base_text,
            "-parse_seqids",
            "-blastdb_version",
            "4",
        ])
        .output()
        .unwrap();
    assert!(!refused.status.success(), "{refused:?}");
    assert!(
        String::from_utf8_lossy(&refused.stderr).contains("[BLAST] MAX_SEQID_LENGTH (100)"),
        "{refused:?}"
    );
    assert!(!std::path::Path::new(&format!("{base_text}.psq")).exists());
    let accepted = clean_losat_command()
        .env("NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH", "100")
        .args([
            "makeblastdb",
            "-in",
            input_text,
            "-dbtype",
            "prot",
            "-out",
            base_text,
        ])
        .output()
        .unwrap();
    assert!(accepted.status.success(), "{accepted:?}");
    remove(&V5);
    // A value makeblastdb ignores (StringToUInt cannot convert it, or 0) is accepted with
    // -parse_seqids (review B N4; makeblastdb 2.17.0: pr2/fix/n4probe).
    for ignored in ["abc", "0", ""] {
        let accepted = clean_losat_command()
            .env("NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH", ignored)
            .args([
                "makeblastdb",
                "-in",
                input_text,
                "-dbtype",
                "prot",
                "-out",
                base_text,
                "-parse_seqids",
                "-blastdb_version",
                "4",
            ])
            .output()
            .unwrap();
        assert!(accepted.status.success(), "{ignored:?}: {accepted:?}");
        remove(&V5);
    }
    // `-taxid_map -` reads the map from the standard input (CArg_InputFile maps `-` to cin;
    // review B S1).
    fs::write(
        &input,
        ">ref|NP_000001.1| first\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|second b\nLEERLGLIEVQAPILSRVGDGTQDN\n",
    )
    .unwrap();
    let mut child = clean_losat_command()
        .args([
            "makeblastdb",
            "-in",
            input_text,
            "-dbtype",
            "prot",
            "-out",
            base_text,
            "-parse_seqids",
            "-taxid_map",
            "-",
        ])
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    {
        use std::io::Write as _;
        let mut stdin = child.stdin.take().unwrap();
        stdin
            .write_all(b"NP_000001.1 9606\nsecond 10090\n")
            .unwrap();
    }
    let output = child.wait_with_output().unwrap();
    assert!(output.status.success(), "{output:?}");
    let metadata = fs::read_to_string(format!("{base_text}.pjs")).unwrap();
    assert!(metadata.contains("\"number-of-taxids\": 2"), "{metadata}");
    remove(&V5);
    // `-taxid_map ""`: "File is not accessible" without a name part (s_ArgExptMsg; review B
    // N1), exit status 1. A directory is read as an empty map (a read error ends the reading;
    // review B S2).
    let output = run_losat(&[
        "makeblastdb",
        "-in",
        input_text,
        "-dbtype",
        "prot",
        "-out",
        base_text,
        "-parse_seqids",
        "-taxid_map",
        "",
    ]);
    assert_eq!(output.status.code(), Some(1), "{output:?}");
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "Command line argument error: Argument \"taxid_map\". File is not accessible\n"
    );
    remove(&V5);
    #[cfg(target_os = "linux")]
    {
        let directory = std::env::temp_dir();
        let output = run_losat(&[
            "makeblastdb",
            "-in",
            input_text,
            "-dbtype",
            "prot",
            "-out",
            base_text,
            "-parse_seqids",
            "-taxid_map",
            directory.to_str().unwrap(),
        ]);
        assert!(output.status.success(), "{output:?}");
        assert!(output.stderr.is_empty(), "{output:?}");
        let metadata = fs::read_to_string(format!("{base_text}.pjs")).unwrap();
        assert!(!metadata.contains("number-of-taxids"), "{metadata}");
        remove(&V5);
    }
    fs::remove_file(&input).unwrap();
    fs::remove_file(&map).unwrap();
}

const MT_PROTEIN: &str = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV";
const MT_NUCLEOTIDE: &str = "ATGAAAACCGCCTATATTGCCAAACAGCGTCAGATTAGCTTTGTGAAAAGCCATTTTAGCCGTCAGCTGGAAGAACGTCTGGGCCTGATTGAAGTGCAGGCCCCGATTCTGAGCCGTGTGGGCGATGGCACCCAGGATAACCTGAGCGGCGCCGAAAAAGCCGTGCAGGTGAAAGTGAAAGCCCTGCCGGATGCCCAGTTTGAAGTGGTG";

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3176-3181
// ```c++
//     	arg_desc.SetConstraint(kArgMTMode,
//                                new CArgAllowValuesBetween(0, 2, true));
//     	arg_desc.SetDependency(kArgMTMode,
//                                CArgDescriptions::eRequires,
//                                kArgNumThreads);
// ```
// `-mt_mode` has no effect on LOSAT's output (CLI DECISIONS 5).
#[test]
fn mt_mode_does_not_change_the_output() {
    let protein = temp_path("mt_protein", "faa");
    let nucleotide = temp_path("mt_nucleotide", "fna");
    fs::write(&protein, format!(">p\n{MT_PROTEIN}\n")).unwrap();
    fs::write(&nucleotide, format!(">n\n{MT_NUCLEOTIDE}\n")).unwrap();
    let (protein, nucleotide) = (protein.to_str().unwrap(), nucleotide.to_str().unwrap());
    for (program, query, subject) in [
        ("blastn", nucleotide, nucleotide),
        ("blastp", protein, protein),
        ("blastx", nucleotide, protein),
        ("tblastn", protein, nucleotide),
    ] {
        let base = [
            program, "-query", query, "-subject", subject, "-outfmt", "6",
        ];
        let run = |words: &[&str]| {
            let mut command = clean_losat_command();
            // A build with the diagnostics feature needs this log for two or more threads.
            if cfg!(feature = "diagnostics") {
                command.env("LOSAT_BLASTX_WORKER_LOG", temp_path("mt_worker", "log"));
            }
            command.args(words).output().expect("run LOSAT")
        };
        let plain = run(&base);
        assert!(plain.status.success(), "{program}: {:?}", plain);
        assert!(!plain.stdout.is_empty(), "{program}");
        for mode in ["0", "1", "2"] {
            let mut words = base.to_vec();
            words.extend(["-mt_mode", mode, "-num_threads", "2"]);
            let output = run(&words);
            assert!(output.status.success(), "{program} {mode}: {:?}", output);
            assert_eq!(output.stdout, plain.stdout, "{program} -mt_mode {mode}");
        }
        // Out of range, and without `-num_threads`: argument errors.
        for extra in [
            &["-mt_mode", "3", "-num_threads", "2"][..],
            &["-mt_mode", "-1", "-num_threads", "2"],
            &["-mt_mode", "1"],
        ] {
            let mut words = base.to_vec();
            words.extend(extra);
            let output = run(&words);
            assert!(!output.status.success(), "{program} {extra:?}");
            assert!(output.stdout.is_empty(), "{program} {extra:?}");
        }
    }
    // NCBI's tblastx and psiblast have no -mt_mode.
    for program in ["tblastx", "psiblast"] {
        let output = run_losat(&[
            program,
            "-query",
            protein,
            "-subject",
            protein,
            "-mt_mode",
            "1",
            "-num_threads",
            "2",
        ]);
        assert!(!output.status.success(), "{program}");
    }
    fs::remove_file(protein).unwrap();
    fs::remove_file(nucleotide).unwrap();
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3381-3397
// ```c++
//     if( args.Exist( kArgUseIndex ) &&
//         !(args.Exist( kArgSubject ) && args[kArgSubject])) {
// ...
//         if( args.Exist( kTask ) && args[kTask] &&
//                 args[kTask].AsString() != "megablast" ) {
//             use_index = false;
//         }
// ```
// NCBI BLAST+ 2.17.0 (the official Linux x64 build, 2026-11-11, line LIB P-S3): `-use_index
// true` changes nothing with `-subject` or another task; a megablast search of a database
// without an index writes this warning and the output of the search without `-use_index`;
// with an index (makembindex) NCBI searches with it, which LOSAT refuses.
#[test]
fn use_index_true_runs_as_ncbi_and_an_indexed_database_is_refused() {
    let dir = temp_path("use_index", "dir");
    fs::create_dir_all(&dir).unwrap();
    fs::write(dir.join("q.fna"), format!(">n\n{MT_NUCLEOTIDE}\n")).unwrap();
    let run = |words: &[&str]| {
        clean_losat_command()
            .current_dir(&dir)
            .args(words)
            .output()
            .expect("run LOSAT")
    };
    let made = run(&[
        "makeblastdb",
        "-dbtype",
        "nucl",
        "-in",
        "q.fna",
        "-out",
        "db/plain",
    ]);
    assert!(made.status.success(), "{made:?}");
    let subject = [
        "blastn", "-query", "q.fna", "-subject", "q.fna", "-outfmt", "6",
    ];
    let plain = run(&subject);
    assert!(plain.status.success() && !plain.stdout.is_empty());
    for value in [
        "false", "FALSE", "f", "no", "0", "true", "t", "1", "yes", "on",
    ] {
        let output = run(&[&subject[..], &["-use_index", value]].concat());
        assert!(output.status.success(), "-use_index {value}: {output:?}");
        assert_eq!(output.stdout, plain.stdout, "-use_index {value}");
        assert!(output.stderr.is_empty(), "-use_index {value}");
    }
    let refused = run(&[&subject[..], &["-use_index", "maybe"]].concat());
    assert_eq!(refused.status.code(), Some(2));
    // A task other than megablast ignores it.
    let blastn = [
        "blastn", "-query", "q.fna", "-db", "db/plain", "-outfmt", "6", "-task", "blastn",
    ];
    let plain_blastn = run(&blastn);
    let output = run(&[&blastn[..], &["-use_index", "true"]].concat());
    assert!(
        output.status.success() && output.stderr.is_empty(),
        "{output:?}"
    );
    assert_eq!(output.stdout, plain_blastn.stdout);
    // Megablast without an index: NCBI's warning, then the search without the index.
    let root = "/home/coremake/release_build/build/PrepareRelease_Linux64-Alma_JSID_01_170042_130.14.18.128_9008__PrepareRelease_Linux64-Alma_1751374548/c++/compilers/unix/../..";
    let warning = format!(
        "Warning: (1.1) [blastn] NCBI C++ Exception:\n    T0 \"{root}/src/algo/blast/api/blast_dbindex.cpp\", line 795: Error: (CDbIndex_Exception::bad index creation option) BLAST::ncbi::blast::CIndexedDb_New::CIndexedDb_New() - no database volume has an index\n\nNCBI C++ Exception:\n    T0 \"{root}/src/algo/blast/api/blast_dbindex.cpp\", line 1026: Error: (CDbIndex_Exception::bad index creation option) BLAST::ncbi::blast::CIndexedDb_Old::CIndexedDb_Old() - no index file specified or index 'db/plain*' not found.\n Database index will not be used.\n"
    );
    let megablast = [
        "blastn", "-query", "q.fna", "-db", "db/plain", "-outfmt", "6",
    ];
    let plain_megablast = run(&megablast);
    assert!(plain_megablast.status.success() && !plain_megablast.stdout.is_empty());
    let output = run(&[&megablast[..], &["-use_index", "true"]].concat());
    assert!(output.status.success(), "{output:?}");
    assert_eq!(output.stdout, plain_megablast.stdout);
    assert_eq!(String::from_utf8_lossy(&output.stderr), warning);
    // The archive records the name and the warning (NCBI's `SetUseIndex` and diagnostics).
    let archive = run(&[
        "blastn",
        "-query",
        "q.fna",
        "-db",
        "db/plain",
        "-outfmt",
        "11",
        "-use_index",
        "true",
    ]);
    assert!(archive.status.success(), "{archive:?}");
    let text = String::from_utf8_lossy(&archive.stdout);
    assert!(
        text.contains("name \"ForceMbIndex\",\n          value boolean FALSE\n        },\n        {\n          name \"MbIndexName\",\n          value string \"db/plain\"\n        },\n        {\n          name \"HitlistSize\""),
        "{text}"
    );
    // One message (NCBI saves one diagnostic), its white space cleaned; the text encoding
    // breaks long lines.
    let message = format!(
        "message \"{}\"",
        warning
            .split_whitespace()
            .collect::<Vec<_>>()
            .join(" ")
            .replace('"', "\"\"")
    );
    assert!(text.replace('\n', "").contains(&message), "{text}");
    // A megablast search with a template looks for no index (setup_factory.cpp:440-446): the
    // one-line warning, also with an index on disk.
    fs::write(dir.join("db/plain.shd"), b"").unwrap();
    let template_args = [
        "-template_type",
        "coding",
        "-template_length",
        "16",
        "-word_size",
        "11",
    ];
    let output = run(&[&megablast[..], &["-use_index", "true"], &template_args].concat());
    assert!(output.status.success(), "{output:?}");
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "Warning: (1.1) [blastn] Database indexing is not available for discontiguous searches. Database index will not be used.\n"
    );
    let without = run(&[&megablast[..], &template_args].concat());
    assert_eq!(output.stdout, without.stdout);
    fs::remove_file(dir.join("db/plain.shd")).unwrap();
    // An index: a new-style superheader of a volume, or an old-style index volume.
    for file in ["db/plain.shd", "db/plain.00.idx", "db/plain.42.idx"] {
        fs::write(dir.join(file), b"").unwrap();
        let output = run(&[&megablast[..], &["-use_index", "true"]].concat());
        let error = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(1), "{file}: {error}");
        assert!(output.stdout.is_empty(), "{file}");
        assert!(
            error.contains("-use_index true with a database that has a MegaBLAST index")
                && error.contains("is not supported by LOSAT's BLASTN"),
            "{file}: {error}"
        );
        // NCBI does not look for it without `-use_index true`.
        let output = run(&megablast);
        assert_eq!(output.stdout, plain_megablast.stdout, "{file}");
        fs::remove_file(dir.join(file)).unwrap();
    }
    fs::remove_dir_all(&dir).unwrap();
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:87-88,106-107
// ```c++
//     arg.Reset(new CFilteringArgs( !kQueryIsProtein ));
//     m_Args.push_back(arg);
// ...
//     m_QueryOptsArgs.Reset(new CQueryOptionsArgs(kQueryIsProtein));
//     arg.Reset(m_QueryOptsArgs);
// ```
// NCBI BLAST+ 2.17.0 (line LIB P-S3): BLASTX's error of `-seg` comes from the filtering handler,
// before the query options handler (`-query_loc`) and the formatting handler (its warning of an
// ignored `-num_descriptions`).
#[test]
fn blastx_seg_error_comes_from_the_filtering_handler() {
    let dir = temp_path("blastx_seg", "dir");
    fs::create_dir_all(&dir).unwrap();
    fs::write(dir.join("q.fna"), format!(">n\n{MT_NUCLEOTIDE}\n")).unwrap();
    fs::write(
        dir.join("s.faa"),
        ">p\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV\n",
    )
    .unwrap();
    let error = "BLAST query/options error: Invalid number of arguments to filtering option\nPlease refer to the BLAST+ user manual.\n";
    for extra in [&["-num_descriptions", "5"][..], &["-query_loc", "5-1"], &[]] {
        let output = clean_losat_command()
            .current_dir(&dir)
            .args([
                "blastx", "-query", "q.fna", "-subject", "s.faa", "-outfmt", "6", "-seg", "12 2.2",
            ])
            .args(extra)
            .output()
            .expect("run LOSAT");
        assert_eq!(output.status.code(), Some(1), "{extra:?} {output:?}");
        assert!(output.stdout.is_empty(), "{extra:?}");
        assert_eq!(String::from_utf8_lossy(&output.stderr), error, "{extra:?}");
    }
    fs::remove_dir_all(&dir).unwrap();
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:460-469
// ```c++
// inline CArg_Double::CArg_Double(const string& name, const string& value)
//     : CArg_String(name, value)
// {
//     try {
//         m_Double = NStr::StringToDouble(value, NStr::fDecimalPosixOrLocal);
// ```
// NCBI BLAST+ 2.17.0 (swebok-improve-20261105/sw8/nfix/probe/cli): a real argument takes
// glibc's hexadecimal numbers and NCBI's mantissa with a final `e` (`blastn -evalue`, BLAST XML's
// `Parameters_expect`); `-threshold`'s constraint reads its value without the final `e` (an
// error), `-seg` takes hexadecimal cuts, and a PSSM with vertical tabs in a comment is text
// ASN.1 (psiblast searches it as without them: the frozen `ps.subj.inpssm.tk.6`).
#[test]
fn real_arguments_and_pssm_spaces_read_as_ncbis() {
    let dir = temp_path("ncbi_reals", "dir");
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("q.fa"),
        ">q1\nACGTACGTTAGCTAGCTAGGATCGATCGATGCTAGCTAGCTAGCATCGATCGATCGTAGCTAGCTAGCTAGC\n",
    )
    .unwrap();
    fs::write(
        dir.join("s.fa"),
        ">s1\nTTTTACGTACGTTAGCTAGCTAGGATCGATCGATGCTAGCTAGCTAGCATCGATCGATCGTAGCTAGCTAGCTAGCAAAA\n",
    )
    .unwrap();
    fs::write(
        dir.join("pq.faa"),
        ">p1 test protein\nMKVLAAGIVALLLAAQQQQQQQQQQQQQQQQPPPPPPPPGSGSGSGSGSGSAEDKLMNVWYRTSHCFIPKLEDGHKLMNRSTVWY\n",
    )
    .unwrap();
    fs::write(
        dir.join("ps.faa"),
        ">s1 subject protein\nMKVLAAGIVALLLAAQQQQQQQQQAQQQQQPPPPPPPPGSGSGSGSGAGSAEDKLMNVWYRTSHCFIPKLEDGHKLMNRSTVWYAAAA\n",
    )
    .unwrap();
    let run = |args: &[&str]| {
        clean_losat_command()
            .current_dir(&dir)
            .args(args)
            .output()
            .expect("run LOSAT")
    };
    let blastn = [
        "blastn", "-query", "q.fa", "-subject", "s.fa", "-outfmt", "5",
    ];
    for (value, expect) in [
        ("0x10", "16"),
        ("0x1p-3", "0.125"),
        ("1e", "1"),
        ("12345678901234567890e", "1.23456789012346e+19"),
        ("0.0000000000000001e", "1e-16"),
        ("-nan(1)", "NaN"),
    ] {
        let output = run(&[&blastn[..], &["-evalue", value]].concat());
        assert_eq!(output.status.code(), Some(0), "{value}: {output:?}");
        let xml = String::from_utf8(output.stdout).unwrap();
        assert!(
            xml.contains(&format!("<Parameters_expect>{expect}</Parameters_expect>")),
            "{value}: {xml}"
        );
    }
    for value in ["1x", "1e+", "2,5", "inf", "nan(5)"] {
        let output = run(&[&blastn[..], &["-evalue", value]].concat());
        assert_eq!(output.status.code(), Some(2), "{value}: {output:?}");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("expected a number"),
            "{value}: {output:?}"
        );
    }
    let blastp = [
        "blastp", "-query", "pq.faa", "-subject", "ps.faa", "-outfmt", "6",
    ];
    let threshold_16 = run(&[&blastp[..], &["-threshold", "16"]].concat());
    assert_eq!(threshold_16.status.code(), Some(0), "{threshold_16:?}");
    let threshold_hex = run(&[&blastp[..], &["-threshold", "0x10"]].concat());
    assert_eq!(threshold_hex.stdout, threshold_16.stdout);
    for value in ["1e", "12345678901234567890e"] {
        let output = run(&[&blastp[..], &["-threshold", value]].concat());
        assert_eq!(output.status.code(), Some(2), "{value}: {output:?}");
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("Illegal value"),
            "{value}: {output:?}"
        );
    }
    let seg_decimal = run(&[&blastp[..], &["-seg", "12 2.25 2.5"]].concat());
    let seg_hex = run(&[&blastp[..], &["-seg", "12 0x1.2p1 0x1.4p1"]].concat());
    assert_eq!(seg_decimal.status.code(), Some(0), "{seg_decimal:?}");
    assert_eq!(seg_hex.status.code(), Some(0), "{seg_hex:?}");
    assert!(!seg_hex.stdout.is_empty());
    assert_eq!(seg_hex.stdout, seg_decimal.stdout);
    // A comment line of 3,000 vertical tabs after the first line.
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut pssm = Vec::new();
    std::io::Read::read_to_end(
        &mut flate2::read::GzDecoder::new(
            fs::read(root.join("tests/fixtures/pssm/ps4_tk_paja.asn.gz"))
                .unwrap()
                .as_slice(),
        ),
        &mut pssm,
    )
    .unwrap();
    let first_line = pssm.iter().position(|&byte| byte == b'\n').unwrap() + 1;
    let mut with_tabs = pssm[..first_line].to_vec();
    with_tabs.extend_from_slice(b"-- ");
    with_tabs.extend(std::iter::repeat_n(0x0b, 3000));
    with_tabs.push(b'\n');
    with_tabs.extend_from_slice(&pssm[first_line..]);
    fs::write(dir.join("vt.asn"), with_tabs).unwrap();
    let subject = root.join("tests/fasta/PajaWSV.faa");
    let output = run(&[
        "psiblast",
        "-in_pssm",
        "vt.asn",
        "-subject",
        subject.to_str().unwrap(),
        "-comp_based_stats",
        "0",
        "-outfmt",
        "6",
    ]);
    assert_eq!(output.status.code(), Some(0), "{output:?}");
    assert_eq!(
        output.stdout,
        fs::read(root.join("tests/fixtures/outfmt0/ps.subj.inpssm.tk.6.out")).unwrap()
    );
    fs::remove_dir_all(&dir).unwrap();
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1147-1154
// ```c
//         sparamsp = SegParametersNewAa();
//         sparamsp->overlaps = TRUE;
//         if (seg_options->window > 0)
//             sparamsp->window = seg_options->window;
//         if (seg_options->locut > 0.0)
//             sparamsp->locut = seg_options->locut;
//         if (seg_options->hicut > 0.0)
//             sparamsp->hicut = seg_options->hicut;
// ```
// NCBI BLAST+ 2.17.0 (swebok-improve-20261105/sw8/nfix/probe/seg2): a `-seg` locut or hicut
// that is NaN keeps SEG's default (`-nan > 0.0` is false), an infinite locut makes every window
// low in complexity (no hit), and the filter string prints the values as glibc's `%1.1f`.
#[test]
fn seg_cuts_that_are_not_finite_search_as_ncbis() {
    let dir = temp_path("seg_not_finite", "dir");
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("pq.faa"),
        ">p1 test protein\nMKVLAAGIVALLLAAQQQQQQQQQQQQQQQQPPPPPPPPGSGSGSGSGSGSAEDKLMNVWYRTSHCFIPKLEDGHKLMNRSTVWY\n",
    )
    .unwrap();
    fs::write(
        dir.join("ps.faa"),
        ">s1 subject protein\nMKVLAAGIVALLLAAQQQQQQQQQAQQQQQPPPPPPPPGSGSGSGSGAGSAEDKLMNVWYRTSHCFIPKLEDGHKLMNRSTVWYAAAA\n",
    )
    .unwrap();
    let blastp = |seg: &str, outfmt: &str| {
        clean_losat_command()
            .current_dir(&dir)
            .args(["blastp", "-query", "pq.faa", "-subject", "ps.faa"])
            .args(["-seg", seg, "-outfmt", outfmt])
            .output()
            .expect("run LOSAT")
    };
    let default_rows = "p1\ts1\t100.000\t33\t0\t0\t53\t85\t52\t84\t1.70e-25\t77.0\n";
    for (seg, rows) in [
        ("12 2.2 2.5", default_rows),
        ("12 2.2 -nan", default_rows),
        ("12 +nan(1) 2.5", default_rows),
        ("12 -inf 2.5", default_rows),
        ("12 2.2 +inf", default_rows),
        ("12 +inf 2.5", ""),
        ("12 1e400 2.5", ""),
    ] {
        let output = blastp(seg, "6");
        assert_eq!(output.status.code(), Some(0), "{seg}: {output:?}");
        assert_eq!(String::from_utf8_lossy(&output.stdout), rows, "{seg}");
    }
    for (seg, filter) in [
        ("12 +inf 2.5", "S 12 inf 2.5;"),
        ("12 2.2 -nan", "S 12 2.2 -nan;"),
        ("12 +nan(1) 2.5", "S 12 nan 2.5;"),
        ("12 -inf 2.5", "S 12 -inf 2.5;"),
    ] {
        let output = blastp(seg, "5");
        assert_eq!(output.status.code(), Some(0), "{seg}: {output:?}");
        assert!(
            String::from_utf8_lossy(&output.stdout)
                .contains(&format!("<Parameters_filter>{filter}</Parameters_filter>")),
            "{seg}"
        );
    }
    // TBLASTN prints the options too, not the SEG parameters made of them (it printed
    // `S 12 inf inf;`, `S 12 3.0 3.0;` and `L;` for the others): NCBI BLAST+ 2.17.0 `tblastn
    // -outfmt 5`, also in the `filter` of outfmt 14 and 16.
    fs::write(
        dir.join("ns.fna"),
        ">ns1 back-translated subject\nATGAAAGTGCTGGCTGCTGGCATTGTGGCTCTGCTGCTGGCTGCTCAGCAGCAGCAGCAGCAGCAGCAGCAG\n",
    )
    .unwrap();
    for (seg, filter) in [
        ("12 +inf 2.5", "S 12 inf 2.5;"),
        ("12 -nan -nan", "S 12 -nan -nan;"),
        ("12 -inf 2.5", "S 12 -inf 2.5;"),
        ("12 3 2.5", "S 12 3.0 2.5;"),
        ("12 0 2.5", "S 12 0.0 2.5;"),
    ] {
        let output = clean_losat_command()
            .current_dir(&dir)
            .args(["tblastn", "-query", "pq.faa", "-subject", "ns.fna"])
            .args(["-seg", seg, "-outfmt", "5"])
            .output()
            .expect("run LOSAT");
        assert_eq!(output.status.code(), Some(0), "{seg}: {output:?}");
        assert!(
            String::from_utf8_lossy(&output.stdout)
                .contains(&format!("<Parameters_filter>{filter}</Parameters_filter>")),
            "tblastn {seg}"
        );
    }
    fs::remove_dir_all(&dir).unwrap();
}
