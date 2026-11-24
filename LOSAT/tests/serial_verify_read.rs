//! `SERIAL_VERIFY_DATA_READ` and `SERIAL_VERIFY_DATA_GET` outside rpsbproc (line RG D3c-3d). A
//! value that turns NCBI's check off (`NO`, `NEVER`, `DEFVALUE`, `DEFVALUE_ALWAYS`, in any case)
//! makes NCBI BLAST+ 2.17.0 post the member as missing and go on where an input fails that check;
//! LOSAT refuses the run there, after the output NCBI writes before its default error. Every other
//! value, the other variable, and every input that passes the check run as without the variable
//! (probes D3c-3d A, B and C: each input under each value, NCBI and LOSAT).

use std::ffi::OsStr;
use std::fmt::Debug;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

const READ: &str = "SERIAL_VERIFY_DATA_READ";
const GET: &str = "SERIAL_VERIFY_DATA_GET";

/// A directory of its own for one test.
fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "losat-serial-verify-read-{tag}-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// One run of `losat <args>` from `LOSAT/` without the variables that change a run, plus
/// `variable`; the registry is the empty directory `dir`.
fn losat<S: AsRef<OsStr>>(dir: &Path, args: &[S], variable: Option<(&str, &str)>) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if (text.starts_with("LOSAT_") && text != "LOSAT_BLASTX_WORKER_LOG")
            || text.starts_with("NCBI_")
            || text.starts_with("SERIAL_")
            || text == "BLASTDB"
        {
            command.env_remove(&name);
        }
    }
    command
        .env("NCBI_CONFIG_PATH", dir)
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .args(args)
        .stdin(Stdio::null());
    if let Some((name, value)) = variable {
        command.env(name, value);
    }
    command.output().expect("run LOSAT CLI")
}

/// The output without makeblastdb's clock (its "current time" and "seconds" lines).
fn clocked(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes)
        .lines()
        .filter(|line| !line.contains("current time") && !line.contains(" seconds"))
        .collect::<Vec<_>>()
        .join("\n")
}

/// `args` on an input that fails the check of `variable`: each value that keeps the checks, and
/// an off value of the other variable, give what the run without them gives; each value that
/// turns `variable`'s check off gives the same output and LOSAT's refusal, "the environment
/// variable <variable>=<value> with <refused>", with exit status 1.
fn refused<S: AsRef<OsStr> + Debug>(dir: &Path, args: &[S], variable: &str, refused: &str) {
    let unset = losat(dir, args, None);
    assert_ne!(unset.status.code(), Some(0), "{args:?}");
    let other = if variable == READ { GET } else { READ };
    for (name, value) in [
        (variable, "YES"),
        (variable, "ALWAYS"),
        (variable, "maybe"),
        (variable, ""),
        (other, "NO"),
    ] {
        let run = losat(dir, args, Some((name, value)));
        assert_eq!(
            (
                run.status.code(),
                clocked(&run.stdout),
                clocked(&run.stderr)
            ),
            (
                unset.status.code(),
                clocked(&unset.stdout),
                clocked(&unset.stderr)
            ),
            "{args:?} {name}={value}"
        );
    }
    for value in ["NO", "never", "DefValue_Always"] {
        let run = losat(dir, args, Some((variable, value)));
        let err = String::from_utf8_lossy(&run.stderr);
        assert_eq!(run.status.code(), Some(1), "{args:?} {value}: {err}");
        assert_eq!(
            clocked(&run.stdout),
            clocked(&unset.stdout),
            "{args:?} {value}"
        );
        let text = format!("the environment variable {variable}={value} with {refused}\n");
        assert!(err.ends_with(&text), "{args:?} {value}: {err}");
    }
}

/// `args` on an input that passes the checks: the values that turn them off change nothing.
fn unchanged<S: AsRef<OsStr> + Debug>(dir: &Path, args: &[S]) {
    let unset = losat(dir, args, None);
    assert_eq!(unset.status.code(), Some(0), "{args:?}");
    for (name, value) in [(READ, "NO"), (GET, "NEVER")] {
        let run = losat(dir, args, Some((name, value)));
        assert_eq!(
            (
                run.status.code(),
                clocked(&run.stdout),
                clocked(&run.stderr)
            ),
            (
                unset.status.code(),
                clocked(&unset.stdout),
                clocked(&unset.stderr)
            ),
            "{args:?} {name}={value}"
        );
    }
}

/// A PSSM whose read lacks a mandatory member (`c20`: a User-object without its type; `norows`:
/// the matrix without numRows, in ASN.1 text and in XML), and one whose query's feature table has
/// a column without the field-name that NCBI gets (`get`), for psiblast and tblastn `-in_pssm`.
#[test]
fn a_pssm_that_fails_a_check_is_refused_with_the_check_off() {
    let dir = temp_dir("pssm");
    let valid = std::fs::read_to_string("tests/fixtures/db_search/db.ps.asn.rnd.6.asn").unwrap();
    assert_eq!(valid.matches("    numRows 28,\n").count(), 1);
    let norows = dir.join("norows.asn");
    std::fs::write(&norows, valid.replacen("    numRows 28,\n", "", 1)).unwrap();
    let annotated = std::fs::read_to_string("tests/fixtures/pssm/d3a_o03_annot_seqtable.asn")
        .expect("the D3a feature table fixture");
    let columns =
        "columns { { header { field-id partial }, data bit-bvector { size 1, data '80'H } }, \
         { header { field-id partial }, data bit-bvector { size 1, data '80'H } } }";
    assert_eq!(annotated.matches(columns).count(), 1);
    let get = dir.join("get.asn");
    std::fs::write(
        &get,
        annotated.replacen(
            columns,
            "columns { { header { field-id ext }, data string { \"v\" } } }",
            1,
        ),
    )
    .unwrap();
    let one = std::fs::read_to_string("tests/fixtures/pssm/ps9_one.xml").unwrap();
    let xml_rows = "      <Pssm_numRows>28</Pssm_numRows>\n";
    assert_eq!(one.matches(xml_rows).count(), 1);
    let norows_xml = dir.join("norows.xml");
    std::fs::write(&norows_xml, one.replacen(xml_rows, "", 1)).unwrap();
    let (norows, get) = (norows.to_str().unwrap(), get.to_str().unwrap());
    let norows_xml = norows_xml.to_str().unwrap();
    let psiblast = |pssm: &str, outfmt: &str| {
        [
            "psiblast",
            "-db",
            "tests/db_search/paja",
            "-in_pssm",
            pssm,
            "-comp_based_stats",
            "0",
            "-outfmt",
            outfmt,
            "-num_threads",
            "1",
        ]
        .map(String::from)
    };
    let lacks = "a PSSM (-in_pssm) that lacks a mandatory member is not supported by LOSAT's";
    refused(
        &dir,
        &psiblast("tests/fixtures/pssm/d3a_c20_missing.asn", "6"),
        READ,
        &format!("{lacks} PSIBLAST"),
    );
    refused(
        &dir,
        &psiblast(norows, "6"),
        READ,
        &format!("{lacks} PSIBLAST"),
    );
    refused(
        &dir,
        &psiblast(norows_xml, "6"),
        READ,
        &format!("{lacks} PSIBLAST"),
    );
    let not_set =
        "a feature table of a PSSM's query (-in_pssm) of which NCBI BLAST+ gets a member \
         that is not set is not supported by LOSAT's";
    // NCBI writes the report's prolog before its exception (outfmt 0), and LOSAT the same bytes.
    refused(
        &dir,
        &psiblast(get, "0"),
        GET,
        &format!("{not_set} PSIBLAST"),
    );
    let tblastn = |pssm: &str| {
        [
            "tblastn",
            "-db",
            "tests/db_search/tslice",
            "-in_pssm",
            pssm,
            "-outfmt",
            "6",
            "-num_threads",
            "1",
        ]
        .map(String::from)
    };
    refused(
        &dir,
        &tblastn("tests/fixtures/pssm/d3a_c20_missing.asn"),
        READ,
        &format!("{lacks} TBLASTN"),
    );
    refused(
        &dir,
        &tblastn(norows_xml),
        READ,
        &format!("{lacks} TBLASTN"),
    );
    refused(&dir, &tblastn(get), GET, &format!("{not_set} TBLASTN"));
    unchanged(
        &dir,
        &psiblast("tests/fixtures/db_search/db.ps.asn.rnd.6.asn", "6"),
    );
    unchanged(&dir, &tblastn("tests/fixtures/pssm/ps9_one.xml"));
    let _ = std::fs::remove_dir_all(&dir);
}

/// A profile without numRows for makeprofiledb, and mask data without algo-options for
/// makeblastdb `-mask_data` (NCBI writes its first lines before the error; LOSAT the same).
#[test]
fn profiles_and_mask_data_that_lack_a_member_are_refused_with_the_check_off() {
    let dir = temp_dir("profiles");
    let profile = std::fs::read_to_string("tests/fixtures/makeprofiledb/in/pc.smp").unwrap();
    assert_eq!(profile.matches("    numRows 28,\n").count(), 1);
    let smp = dir.join("norows.smp");
    std::fs::write(&smp, profile.replacen("    numRows 28,\n", "", 1)).unwrap();
    let list = dir.join("norows.lst");
    std::fs::write(&list, format!("{}\n", smp.display())).unwrap();
    let out = dir.join("p").display().to_string();
    refused(
        &dir,
        &["makeprofiledb", "-in", list.to_str().unwrap(), "-out", &out],
        READ,
        "a profile (-in) that lacks a mandatory member is not supported by LOSAT's MAKEPROFILEDB",
    );
    let out = dir.join("m").display().to_string();
    let mask = |data: &str| {
        [
            "makeblastdb",
            "-in",
            "tests/seqdb/mask_nucl.fna",
            "-dbtype",
            "nucl",
            "-mask_data",
            data,
            "-out",
            &out,
        ]
        .map(String::from)
    };
    refused(
        &dir,
        &mask("tests/seqdb/mask_parse/h_member_missing_options.asnt"),
        READ,
        "mask data (-mask_data) that lacks a mandatory member is not supported by LOSAT's makeblastdb",
    );
    unchanged(&dir, &mask("tests/seqdb/mask_parse/dust.asnt"));
    let _ = std::fs::remove_dir_all(&dir);
}

const PROTEINS: &str = ">seqC protein C of test
YKNTARICGENSIPVEWIAHQKGPGDFYYRFFAAHHGGLMHVHGHPLANQGFKDMLYWAYMDLNLSMGSSGCKANPAVQNPWARCGYHEIRNTNTKREWNLCQDHMTYNFMKVDLMLGDF
>seqA protein A of test
LSGCDYVPCIYNKRQFCCSMHFWFQEGQNFCQLFRYGTRSMSKLSPFEPVGSMGDSKTVTNDNWCLNVKSKLMGWASVKMKRLTNNKNQNGRNMTFTGHNSLDQGYWTQLYVKAHGWRYG
";

const BASES: &str = ">seqC nucleotide C
AGCATCAATGATCGAGCCGTGGAAAAAACGTGACTCGCGGACCAGCCTTTAGGTCTTCTACTTAACTACAACTGTTCCGCGGCGGCATTGCC
CTTAACTAGCGTTACTAACTAGAGTTTTACTGACGGAAAGTGAGCAAAGGCTAACGTTATTCCGTGAGCACGGGACATCCATTCTTCGTGAG
>seqA nucleotide A
AACTTACATTTTACGCGCTAGCTTCGCTGGAACTAATATCCATGTCTCAGAACTAGCGGCCGAGAATGGGTTCCGAATCCTAAACTCCGACA
TGAGTTAAGGTTGCATACTAGGTCTGATACTAAAAGCGGGGTCAGGAGTCCGTCCAGAATATAATATTCAAAAATGAGATGGTGGAGTTTCC
";

/// A version 4 database of `fasta` (`dbtype` prot or nucl) at `dir/name`, built by LOSAT's
/// makeblastdb; with `cut`, the header record of OID 0 is rewritten at its length without its
/// Seq-ids (a longer title in their place), which NCBI's read fails with "member seqid expected".
fn database(dir: &Path, name: &str, fasta: &str, dbtype: &str, cut: bool) -> String {
    let input = dir.join(format!("{name}.fa"));
    std::fs::write(&input, fasta).unwrap();
    let out = dir.join(name).display().to_string();
    let built = losat(
        dir,
        &[
            "makeblastdb",
            "-in",
            input.to_str().unwrap(),
            "-dbtype",
            dbtype,
            "-blastdb_version",
            "4",
            "-out",
            &out,
        ],
        None,
    );
    assert_eq!(built.status.code(), Some(0), "{built:?}");
    if cut {
        let letter = &dbtype[..1];
        let index = std::fs::read(format!("{out}.{letter}in")).unwrap();
        let word = |at: usize| u32::from_be_bytes(index[at..at + 4].try_into().unwrap()) as usize;
        // Version, type, the title and the date (each a length and its bytes), the count, the
        // total length (8 bytes) and the longest (4), then the header offsets.
        let mut at = 8;
        at += 4 + word(at);
        at += 4 + word(at);
        let first_end = word(at + 4 + 12 + 4);
        let path = format!("{out}.{letter}hr");
        let mut records = std::fs::read(&path).unwrap();
        let title = first_end - 14;
        assert!(title < 128, "{title}");
        let mut record = vec![0x30, 0x80, 0x30, 0x80, 0xa0, 0x80, 0x1a, title as u8];
        record.extend(std::iter::repeat_n(b't', title));
        record.extend_from_slice(&[0; 6]);
        records[..first_end].copy_from_slice(&record);
        std::fs::write(&path, records).unwrap();
    }
    out
}

/// A database whose first header record lacks its seqid, read by each kind of reader of the
/// definition lines: a search (blastp, after NCBI's outfmt 0 prolog, and blastn), blastdbcmd,
/// makeblastdb `-input_type blastdb` and a masker.
#[test]
fn a_header_record_that_lacks_its_seqid_is_refused_with_the_check_off() {
    let dir = temp_dir("headers");
    let cut = database(&dir, "pc", PROTEINS, "prot", true);
    let whole = database(&dir, "pw", PROTEINS, "prot", false);
    let nucl = database(&dir, "nc", BASES, "nucl", true);
    let query = dir.join("q.faa");
    std::fs::write(
        &query,
        PROTEINS.split_inclusive('\n').take(2).collect::<String>(),
    )
    .unwrap();
    let nquery = dir.join("q.fna");
    std::fs::write(
        &nquery,
        BASES.split_inclusive('\n').take(3).collect::<String>(),
    )
    .unwrap();
    let (query, nquery) = (query.to_str().unwrap(), nquery.to_str().unwrap());
    let record =
        "a database header record that lacks a mandatory member is not supported by LOSAT's";
    for outfmt in ["0", "6"] {
        let blastp = [
            "blastp",
            "-db",
            &cut,
            "-query",
            query,
            "-outfmt",
            outfmt,
            "-num_threads",
            "1",
        ];
        refused(&dir, &blastp, READ, &format!("{record} BLASTP"));
        if outfmt == "0" {
            assert!(!losat(&dir, &blastp, None).stdout.is_empty());
        }
        let mut whole_args = blastp;
        whole_args[2] = &whole;
        unchanged(&dir, &whole_args);
    }
    refused(
        &dir,
        &[
            "blastn",
            "-db",
            &nucl,
            "-query",
            nquery,
            "-outfmt",
            "6",
            "-num_threads",
            "1",
        ],
        READ,
        &format!("{record} BLASTN"),
    );
    refused(
        &dir,
        &[
            "blastdbcmd",
            "-db",
            &cut,
            "-dbtype",
            "prot",
            "-entry",
            "all",
            "-outfmt",
            "%t",
        ],
        READ,
        &format!("{record} BLASTDBCMD"),
    );
    let copy = dir.join("copy").display().to_string();
    refused(
        &dir,
        &[
            "makeblastdb",
            "-input_type",
            "blastdb",
            "-in",
            &cut,
            "-dbtype",
            "prot",
            "-out",
            &copy,
        ],
        READ,
        &format!("{record} makeblastdb"),
    );
    refused(
        &dir,
        &["segmasker", "-in", &cut, "-infmt", "blastdb"],
        READ,
        &format!("{record} SEGMASKER"),
    );
    refused(
        &dir,
        &[
            "dustmasker",
            "-in",
            &nucl,
            "-infmt",
            "blastdb",
            "-outfmt",
            "interval",
        ],
        READ,
        &format!("{record} DUSTMASKER"),
    );
    let _ = std::fs::remove_dir_all(&dir);
}
