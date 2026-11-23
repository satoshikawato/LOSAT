use super::*;
use crate::report::serial::asn_types::Value;

/// The text of the stop `abort`: NCBI's report, or LOSAT's refusal.
fn abort_text(abort: &super::super::annot_index::Abort) -> String {
    match abort {
        super::super::annot_index::Abort::Native { message, .. } => {
            String::from_utf8(message.clone()).unwrap()
        }
        super::super::annot_index::Abort::Refused(text) => text.clone(),
    }
}

/// A stream that gives `self.0` and then fails, as a directory given as the standard input
/// does at once (`read`: EISDIR).
struct FailsAfter(&'static [u8]);

impl std::io::Read for FailsAfter {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        if self.0.is_empty() {
            return Err(std::io::Error::other("Is a directory"));
        }
        let count = self.0.len().min(buffer.len());
        buffer[..count].copy_from_slice(&self.0[..count]);
        self.0 = &self.0[count..];
        Ok(count)
    }
}

// O4 of the PS6 review: a standard input that cannot be read is an empty stream, not an
// inaccessible file (NCBI BLAST+ 2.17.0, `-in_pssm - < .`: "Unsupported format for PSSM";
// `-in_msa - < .`: the alignment reader's failure), and the bytes before a failure are kept.
#[test]
fn a_stream_that_fails_ends_where_it_fails() {
    assert_eq!(bytes_until_failure(FailsAfter(b"")), b"");
    assert_eq!(bytes_until_failure(FailsAfter(b"PssmWith")), b"PssmWith");
    assert!(matches!(guess_format(b""), PssmFileFormat::Unsupported));
}

// A file that cannot be opened is the argument's "File is not accessible"; a directory opens
// and gives no bytes (probe `pf.dir` of line PS PS7).
#[test]
fn only_a_file_that_does_not_open_is_inaccessible() {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat-in-pssm-input-{}",
        crate::utils::test_dir::process_id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let missing = dir.join("none.asn");
    let error = read_input_argument("in_pssm", &missing).unwrap_err();
    let expected = crate::cli::inaccessible("in_pssm", &missing);
    assert_eq!(error.to_string(), expected.to_string());
    #[cfg(unix)]
    assert_eq!(read_input_argument("in_msa", &dir).unwrap(), b"");
    let _ = std::fs::remove_dir_all(&dir);
}
use crate::report::serial::pssm::Seqdesc;

/// What `-in_pssm` makes of `bytes` up to the PSSM read: the guess's answer, or the reader's
/// (`read_in_pssm`): "read", NCBI's "Unrecognized format" of a `CSerialException`, the silent
/// "end" of a `CEofException`.
fn outcome(bytes: &[u8]) -> &'static str {
    use crate::blastinput::pssm_text::{read_pssm, read_pssm_binary};
    let read = match guess_format(bytes) {
        PssmFileFormat::TextAsn { start } => read_pssm(&bytes[start..]),
        PssmFileFormat::BinaryAsn { start } => read_pssm_binary(&bytes[start..]),
        PssmFileFormat::Unsupported => return "Unsupported format",
        PssmFileFormat::Exception(_) => return "exception",
        PssmFileFormat::Refused { what: XML_PSSM, .. } => return "XML refused",
        PssmFileFormat::Refused { .. } => return "refused",
    };
    match read.result {
        Ok(_) => "read",
        Err(PssmReadError::Ncbi(ReadError::Serial(_))) => "Unrecognized format",
        Err(PssmReadError::Ncbi(ReadError::Eof)) => "end",
        Err(_) => "other",
    }
}

// The expected answers are NCBI BLAST+ 2.17.0's (psiblast and tblastn `-in_pssm`, the probes of line RG D3-1
// `F/d3a/fg`, named after each case), and review 2's item S1 (`bom2`, `bom4`, `bom7`, `bom8`, `bom10`).
#[test]
fn the_format_guess_is_ncbis() {
    use PssmFileFormat as P;
    let binary = |format: &P| matches!(format, P::BinaryAsn { start: 0 });
    let xml = |format: &P| matches!(format, P::Refused { what: XML_PSSM, .. });
    assert!(matches!(guess_format(b""), P::Unsupported));
    assert!(matches!(guess_format(b"\x1f\x8b\x08\x00"), P::Unsupported));
    assert!(matches!(guess_format(b">q\nMKV\n"), P::Unsupported));
    assert!(binary(&guess_format(b"\x30\x80\xa0\x80")));
    assert!(xml(&guess_format(b"<?xml version=\"1.0\"?>\n")));
    // A vertical tab is C's white space (`isspace`): NCBI 2.17.0's psiblast reads a file that starts with
    // vertical tabs and form feeds before `<?xml` or `<!DOCTYPE` as XML ("Unrecognized format for PSSM ...
    // (must be PssmWithParameters)" for one that is not a PSSM), as without them
    // (swebok-improve-20261105/sw8/nfix/probe/b25; fg `x03_vt_xml`).
    assert!(xml(&guess_format(b"\x0b<?xml version=\"1.0\"?>\n")));
    assert!(xml(&guess_format(
        b"\x0b\x0c <!DOCTYPE PssmWithParameters>\n"
    )));
    // A UTF-8 byte order mark before the text: NCBI's binary reader, which stops at its first byte (review 1
    // `a04_bom`: "Unrecognized format for PSSM", rc 1), with a line end or a comment after it, and when the
    // split of the lines fails (more than 5 % of the bytes above 0x7F: `bom7`, `bom8`).
    let text = b"PssmWithParameters ::= {\n  pssm {\n    isProtein TRUE\n  }\n}\n";
    let bom = [&b"\xef\xbb\xbf"[..], text].concat();
    assert!(matches!(guess_format(text), P::TextAsn { start: 0 }));
    for case in [
        bom.clone(),
        [&b"\xef\xbb\xbf\n"[..], text].concat(),
        [&b"\xef\xbb\xbf-- c --\n"[..], text].concat(),
        b"\xef\xbb\xbfPssm".to_vec(),
        b"\xef\xbb\xbf PssmWithParameters ::= {}\n".to_vec(),
        [&b"\xef\xbb\xbf1"[..], text].concat(),
    ] {
        assert!(binary(&guess_format(&case)), "{case:?}");
        assert_eq!(outcome(&case), "Unrecognized format", "{case:?}");
    }
    // `\1` is not a binary byte (`bom10`); a NEXUS line makes it `eAlignment` (fg `n03_nexus_bom`).
    assert!(matches!(
        guess_format(&[&b"\x01"[..], text].concat()),
        P::Unsupported
    ));
    assert!(matches!(
        guess_format(&[&bom[..], b"-- #NEXUS\n"].concat()),
        P::Unsupported
    ));
    // A first byte `0` is the SEQUENCE tag: a binary PSSM.
    assert!(binary(&guess_format(&[&b"0\xef"[..], &text[1..]].concat())));
    // The text after a NEXUS line: TestFormatNewick reads the stream to its end, where the text reader starts
    // (fg `n01_nexus_text`: exit status 0, nothing written).
    let nexus = [&text[..], b"-- #NEXUS\n"].concat();
    assert!(matches!(guess_format(&nexus), P::TextAsn { start } if start == nexus.len()));
    // The binary reader's first tag (fg `b01`-`b08`): the end of the file in a long tag ends the run
    // without a word, a tag number of 2^24 or more and a tag other than SEQUENCE's are NCBI's
    // `CSerialException`.
    assert_eq!(outcome(b"\x3f\x80"), "end");
    assert_eq!(outcome(b"\x3f\x81\x80"), "end");
    assert_eq!(outcome(b"\x3f\xff\xff\xff\xff\x01"), "Unrecognized format");
    assert_eq!(outcome(b"\x3f\x11\x02"), "Unrecognized format");
    assert_eq!(outcome(b"\x20\x02"), "Unrecognized format");
    // SEQUENCE's tag in its long form begins the PSSM; its member's tag runs into the end of the
    // file inside the class (`HandleEOF`: a `CSerialException`).
    assert!(binary(&guess_format(b"\x3f\x10\x02")));
    assert_eq!(outcome(b"\x3f\x10\x02"), "Unrecognized format");
    assert_eq!(outcome(b"\x30\x80\xa0\x80"), "Unrecognized format");
    // A Glimmer3 frame beyond `int` (fg `f08`): NCBI's `CStringException`.
    assert!(matches!(
        guess_format(b">x\nA 1 2 99999999999 1.0\n"),
        P::Exception(e) if String::from_utf8(e.what()).unwrap().contains("Cannot convert string '99999999999' to int, overflow")
    ));
    // TestFormatFasta's `NStr::Find` past a buffer of four bytes (fg `f11`): refused.
    assert!(matches!(
        guess_format(b">ab\x02"),
        P::Refused {
            what: UNDECIDED_GUESS,
            ..
        }
    ));
}

fn fixture(name: &str) -> Vec<u8> {
    std::fs::read(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures")
            .join(name),
    )
    .unwrap()
}

fn gunzip(name: &str) -> Vec<u8> {
    use std::io::Read;
    let mut out = Vec::new();
    flate2::read::GzDecoder::new(fixture(name).as_slice())
        .read_to_end(&mut out)
        .unwrap();
    out
}

// psiblast's own PSSM (the frozen one of `db.ps.asn.rnd.6`, 30 letters): the query record,
// the 28 rows of each position and the six statistical values of its file.
#[test]
fn a_psiblast_pssm_becomes_a_search() {
    let pssm = crate::blastinput::pssm_text::read_pssm_with_parameters(&fixture(
        "db_search/db.ps.asn.rnd.6.asn",
    ))
    .unwrap();
    let scores = pssm.pssm.final_data.as_ref().unwrap().scores.clone();
    let search = psi_search_of(pssm).unwrap();
    assert_eq!(search.query.local_id, "Query_1");
    assert_eq!(search.query.title(), b"rnd random");
    assert_eq!(search.query.seq().len(), 30);
    assert_eq!(search.pssm.rows.len(), 30);
    // byRow FALSE: position 1's 28 scores follow position 0's.
    assert_eq!(&search.pssm.rows[1][..], &scores[28..56]);
    assert!(search.validation_error.is_none());
    let stats = search.stats;
    assert!(stats.lambda_ungapped.is_some() && stats.kappa.is_some() && stats.h.is_some());
    // The gapped Lambda of psiblast's PSSMs is the table's 0.267 as NCBI reads `{ 267, 10, -3 }`.
    assert_eq!(stats.lambda, Some(267.0 * 1e-3));
}

// Line RG D3-1 review 4 S1, S2: a PSSM whose query's annotations end the run as the query joins
// NCBI's scope (`d3om_t04_dup_partial_id`, two `partial` columns: CAnnotException) gets that stop
// whatever the values NCBI reads later; NCBI BLAST+ 2.17.0 writes plain t04's bytes for each
// variant below (`F/d3a/r3/s2`, review 4 `p3`), and neither a search strategy's request nor a
// PSSM file (`-out_pssm`, review 4 `p2`).
// Line RG D3-1 fix 4 (review 4b S1): the query's Seq-entry joins the scope as NCBI's object
// manager adds it (`annot_index::index_entry`): phase C, the attach of its ids, phase I; a set
// then stops at `GetSeq()`; Seq-data of another alphabet and the subjects' ids come after
// (NCBI BLAST+ 2.17.0: the fixture rows `db.ps.inpssm.d3.om.r4.*`).
#[test]
fn the_query_entry_joins_the_scope_as_ncbis_object_manager_adds_it() {
    // NCBI's report of an exception (the end of its place, its class and its message) in
    // LOSAT's words (`crate::wording::losat_text`): the message.
    let ncbi = |report: &str| {
        let text = format!("Error: NCBI C++ Exception:\n    T0 \"/src/{report}");
        crate::wording::losat_text(&text)
            .trim_start_matches("Error: ")
            .to_string()
    };
    let stop = |name: &str, database: bool| {
        let text = gunzip(&format!("pssm/d3om_r4_{name}.asn.gz"));
        let pssm = crate::blastinput::pssm_text::read_pssm_with_parameters(&text)
            .unwrap_or_else(|error| panic!("{name}: {error:?}"));
        let reader = PssmReader {
            database,
            ..PssmReader::PSIBLAST
        };
        psi_search_and_object(pssm, reader).map(|(search, file)| {
            assert!(!file.valid && search.request.is_none(), "{name}");
            assert!(!search.joins_scope(), "{name}");
            let abort = search.annot.abort.as_ref().unwrap();
            (
                String::from_utf8(search.annot.posts.text()).unwrap(),
                abort_text(abort),
            )
        })
    };
    let (posts, error) = stop("set_t04", true).unwrap();
    assert!(
        posts.is_empty() && error.contains("Duplicate partial column"),
        "{error}"
    );
    // The set's own annotation, then its member's, then the member's Std-seg (phase I).
    let (posts, error) = stop("set_order", true).unwrap();
    assert_eq!(
        posts,
        "Warning: (1207.4) [psiblast] Invalid 'lens' size in denseg\nWarning: (1207.6) [psiblast] Invalid 'starts' size in denseg\n"
    );
    assert!(error.contains("unsupported location type: feat"), "{error}");
    let (posts, error) = stop("set_getseq", true).unwrap();
    assert!(posts.contains("(1207.6)"), "{posts}");
    assert!(
        error.contains(&ncbi("Seq_entry_.cpp\", line 99: Error: (CInvalidChoiceSelection::eFail) CSeq_entry::GetSeq(): Invalid choice selection: NCBI-Seqset::Seq-entry.set")),
        "{error}"
    );
    // Seq-data that the search does not read: the object manager's stop comes first.
    let (_, error) = stop("ncbi8aa_t04", true).unwrap();
    assert!(error.contains("Duplicate partial column"), "{error}");
    // A Seq-id twice: phase C's exception, else the attach's before phase I's.
    let (_, error) = stop("dupid_t04", true).unwrap();
    assert!(error.contains("Duplicate partial column"), "{error}");
    let (_, error) = stop("dupid_a19", true).unwrap();
    assert!(
        error.contains(&ncbi("tse_info.cpp\", line 901: Error: (CObjMgrException::eAddDataError) ncbi::objects::CTSE_Info::x_SetBioseqIds() - duplicate Bioseq id lcl|Query_1 present in\n  seq1: lcl|Query_1 | lcl|Query_1\n  seq2: lcl|Query_1 | lcl|Query_1\n")),
        "{error}"
    );
    // With `-subject` a subject's Seq-id of the handle may have made its info: a local id's
    // text is its own whatever made it (NCBI BLAST+ 2.17.0 `-subject -parse_deflines` with
    // `lcl|query_1`, review 4c `out_x/` x01, x02), a GenBank accession's is not (an EMBL
    // subject's shows, `r5/out_ids5/ncbi/sub.ps.gb_twice_emb`).
    let (_, error) = stop("dupid_a19", false).unwrap();
    assert!(
        error.contains("duplicate Bioseq id lcl|Query_1 present in"),
        "{error}"
    );
    let plain = String::from_utf8(gunzip("pssm/d3om_r4_ncbi8aa_t04.asn.gz")).unwrap();
    let annot = ", annot { { data seq-table { feat-type 1, num-rows 1, columns { { header { field-id partial }, data bit '80'H }, { header { field-id partial }, data bit '80'H } } } } }";
    assert_eq!(plain.matches(annot).count(), 1);
    let go_on = plain.replace(annot, "");
    let pssm = crate::blastinput::pssm_text::read_pssm_with_parameters(go_on.as_bytes()).unwrap();
    let refusal = psi_search_of(pssm).unwrap_err().to_string();
    assert!(refusal.contains("whose Seq-data is ncbi8aa"), "{refusal}");
    // Two Seq-ids of one handle (GenBank and EMBL share an accession alone): NCBI's exception
    // with the spelling of the handle's first Seq-id (`r5/out_ids5/ncbi/ps.one_emb_gb`).
    let two = |ids: &str, database: bool| {
        let text = go_on
            .replace("seq-data ncbi8aa", "seq-data ncbistdaa")
            .replace("local str \"Query_1\"", ids);
        let pssm =
            crate::blastinput::pssm_text::read_pssm_with_parameters(text.as_bytes()).unwrap();
        let reader = PssmReader {
            database,
            ..PssmReader::PSIBLAST
        };
        psi_search_and_object(pssm, reader).map(|(search, _)| {
            search
                .annot
                .abort
                .as_ref()
                .map_or_else(String::new, abort_text)
        })
    };
    let ids =
        "embl { accession \"AAA00001\", version 1 }, genbank { accession \"AAA00001\", version 1 }";
    let error = two(ids, true).unwrap();
    assert!(
        error.contains("duplicate Bioseq id emb|AAA00001.1| present in\n  seq1: emb|AAA00001.1| | emb|AAA00001.1|\n  seq2: emb|AAA00001.1| | emb|AAA00001.1|\n"),
        "{error}"
    );
    // LOSAT's refusal stops the run at the same place (`annot_index::Abort::Refused`).
    let refusal = two(ids, false).unwrap();
    assert!(
        refusal.contains("may write as a subject record spelled it")
            && refusal.ends_with("is not supported by LOSAT's PSIBLAST"),
        "{refusal}"
    );
    // Two versions are two handles: no exception.
    let error = two(
        "genbank { accession \"AAA00001\", version 1 }, genbank { accession \"AAA00001\" }",
        true,
    );
    assert!(!error.unwrap_or_default().contains("duplicate Bioseq id"));
}

#[test]
fn the_object_managers_stop_comes_before_the_pssms_values() {
    let text = String::from_utf8(gunzip("pssm/d3om_t04_dup_partial_id.asn.gz")).unwrap();
    let variant = |from: &str, to: &str| {
        assert_eq!(text.matches(from).count(), 1, "{from}");
        text.replace(from, to)
    };
    let local = "local str \"Query_1\"";
    for (name, pssm) in [
        ("t04", text.clone()),
        ("dims", variant("numColumns 30,", "numColumns 31,")),
        (
            "range",
            variant(
                "scores {\n        -32768,\n        -1,",
                "scores {\n        -32768,\n        40000,",
            ),
        ),
        ("neglambda", variant("{ 267, 10, -3 }", "{ -267, 10, -3 }")),
        (
            "ids",
            variant(local, "local str \"Query_1\",\n        gi 5"),
        ),
        ("noid", variant(&format!("{local}\n"), "")),
        ("emptylocal", variant(local, "local str \"\"")),
        ("textnoacc", variant(local, "genbank { version 1 }")),
        ("molna", variant("mol aa,", "mol na,")),
        ("letters", variant("length 30,", "length 31,")),
    ] {
        let pssm = crate::blastinput::pssm_text::read_pssm_with_parameters(pssm.as_bytes())
            .unwrap_or_else(|error| panic!("{name}: {error:?}"));
        let reader = PssmReader {
            strategy: true,
            database: true,
            ..PssmReader::PSIBLAST
        };
        let (search, file) =
            psi_search_and_object(pssm, reader).unwrap_or_else(|error| panic!("{name}: {error}"));
        assert!(search.annot.abort.is_some(), "{name}");
        assert!(search.request.is_none(), "{name}");
        assert!(!file.valid, "{name}");
    }
}

// NCBI's unit-test PSSMs: one of frequency ratios only whose query has a gi and a swissprot
// Seq-id, read (line RG Q1b; NCBI BLAST+ 2.17.0 psiblast: `Query= P01013 Gene X protein
// (Ovalbumin-related)`, the warning label `129295`, the first id), and one without a query,
// which ends with NCBI's exception of its empty Seq-entry, exit status 255 (row
// db.ps.inpssm.err.noquery).
#[test]
fn ncbis_unit_test_pssms_are_read_or_refused() {
    let read = |name: &str| {
        crate::blastinput::pssm_text::read_pssm_with_parameters(&gunzip(name)).unwrap()
    };
    let search = psi_search_of(read("pssm/pssm_freq_ratios.asn.gz")).unwrap();
    for parse_deflines in [false, true] {
        let query = search.query_record(parse_deflines);
        assert_eq!(
            &*query.query_line(),
            b"P01013 Gene X protein (Ovalbumin-related)"
        );
        assert_eq!(query.warning_id(), b"129295");
    }
    let search = psi_search_of(read("pssm/pssm_zero_fratios.asn.gz")).unwrap();
    assert!(search.request.is_none() && !search.joins_scope());
    assert_eq!(
        search.annot.abort,
        Some(super::super::annot_index::Abort::Native {
            exit: 255,
            message: b"Error: CSeq_entry::GetSeq(): Invalid choice selection: \
          NCBI-Seqset::Seq-entry.not set\n\n"
                .to_vec()
        })
    );
}

// The search strategy's request (line RG F1b-2, review 4 R5-1, R5-4): the PSSM as NCBI's
// application read it (NCBI 2.17.0 writes a ratio-only PSSM without scores, `ps4_tk_ratios`,
// `ps5_ratio30_noparams`), and NCBI's order of the errors before the search
// (`PsiSearch::stop_before_search`).
#[test]
fn the_request_holds_the_pssm_as_read_and_the_engines_errors_come_in_order() {
    use crate::algorithm::pssm_search::PssmQueryStop;
    let read = |name: &str| {
        crate::blastinput::pssm_text::read_pssm_with_parameters(&gunzip(name)).unwrap()
    };
    let reader = |matrix_name, request| PssmReader {
        program: "TBLASTN",
        matrix_name,
        strategy: request,
        database: true,
        archive: None,
    };
    // A PSSM of frequency ratios only: the request holds it as read (no scores), the search
    // and the file the engine's scores (`with_computed_scores`).
    let as_read = read("pssm/ps5_ratio30_noparams.asn.gz");
    assert!(as_read.pssm.final_data.is_none());
    let (search, file) = psi_search_and_object(as_read.clone(), reader("BLOSUM62", true)).unwrap();
    assert_eq!(*search.request_pssm().unwrap(), as_read);
    assert!(file.object.pssm.final_data.is_some() && file.valid);
    assert_eq!(search.stop_before_search(true, true), None);
    // Without a request (no strategy, or `-subject`): no copy.
    let (search, _) = psi_search_and_object(as_read.clone(), reader("BLOSUM62", false)).unwrap();
    assert!(search.request.is_none());
    // A matrix without 11/1: the engine's error comes after the request (review 4 R5-1), the
    // strategy's "Remote PSI-BL2SEQ" before it, and without the strategy the error alone.
    let (search, file) = psi_search_and_object(as_read.clone(), reader("BLOSUM45", true)).unwrap();
    assert_eq!(*search.request_pssm().unwrap(), as_read);
    assert!(!file.valid);
    let scores = PssmQueryStop::Scores(
        "Gap existence and extension values of 11 and 1 not supported for BLOSUM45".to_string(),
    );
    assert_eq!(search.stop_before_search(true, true), Some(scores.clone()));
    assert_eq!(search.stop_before_search(false, true), Some(scores));
    assert_eq!(
        search.stop_before_search(true, false),
        Some(PssmQueryStop::Bl2seq)
    );
    // NCBI's check of the PSSM: no request.
    let scaled = read("pssm/ps4_rnd_scaled.asn.gz");
    let (search, _) = psi_search_and_object(scaled, reader("BLOSUM62", true)).unwrap();
    assert!(search.request.is_none());
    assert!(matches!(
        search.stop_before_search(true, true),
        Some(PssmQueryStop::Invalid(_))
    ));
    assert_eq!(
        search.stop_before_search(true, false),
        Some(PssmQueryStop::Bl2seq)
    );
    // A member that the search does not read (a comment descriptor): the request keeps the
    // object as read, which the PSSM writer writes back whole.
    let comment = read("pssm/ps5_rnd_comment.asn.gz");
    let (search, _) = psi_search_and_object(comment.clone(), reader("BLOSUM62", true)).unwrap();
    assert!(search.request_pssm().is_some_and(|pssm| *pssm == comment));
}

// A PSSM of frequency ratios only (the tk PSSM without its scores) gets NCBI's engine's
// scores, which are the tk PSSM's own (NCBI 2.17.0 searches both to the same bytes, fixture
// rows db.ps.inpssm.cbs.ratios.6 and db.ps.inpssm.cbs.tk.6), and the engine's gapped values;
// its ungapped values stay absent. By row it is refused.
#[test]
fn a_pssm_of_frequency_ratios_only_gets_the_engines_scores() {
    let read = |name: &str| {
        crate::blastinput::pssm_text::read_pssm_with_parameters(&gunzip(name)).unwrap()
    };
    let tk = psi_search_of(read("pssm/ps4_tk_paja.asn.gz")).unwrap();
    let ratios = psi_search_of(read("pssm/ps4_tk_ratios.asn.gz")).unwrap();
    assert_eq!(ratios.pssm, tk.pssm);
    assert_eq!(ratios.freq_ratios, tk.freq_ratios);
    assert_eq!(ratios.stats.lambda, tk.stats.lambda);
    assert_eq!(ratios.stats.h, tk.stats.h);
    let (kappa, file_kappa) = (ratios.stats.kappa.unwrap(), tk.stats.kappa.unwrap());
    assert!((kappa - file_kappa).abs() < 1e-15, "{kappa} {file_kappa}");
    assert_eq!(ratios.stats.lambda_ungapped, None);
    assert_eq!(ratios.stats.kappa_ungapped, None);
    let error = psi_search_of(read("pssm/ps4_rnd_ratios_byrow.asn.gz"))
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("a PSSM of frequency ratios only written by row"),
        "{error}"
    );
}

// `CPsiBlastValidate::Pssm`'s errors are kept for the search (after the prolog), in NCBI's
// order; a PSSM whose sizes do not match is LOSAT's refusal.
#[test]
fn ncbis_checks_and_the_size_refusal() {
    let read = || {
        crate::blastinput::pssm_text::read_pssm_with_parameters(&fixture(
            "db_search/db.ps.asn.rnd.6.asn",
        ))
        .unwrap()
    };
    let mut scaled = read();
    scaled.pssm.final_data.as_mut().unwrap().scaling_factor = Some(2);
    assert_eq!(
        psi_search_of(scaled).unwrap().validation_error.as_deref(),
        Some("PSSM has a scaling factor of 2. PSI-BLAST does not accept scaled PSSMs")
    );
    let mut nucleotide = read();
    nucleotide.pssm.is_protein = Some(false);
    assert_eq!(
        psi_search_of(nucleotide)
            .unwrap()
            .validation_error
            .as_deref(),
        Some("PSSM does not represent protein scoring matrix")
    );
    let mut empty = read();
    empty.pssm.final_data = None;
    empty.pssm.intermediate_data = None;
    assert_eq!(
        psi_search_of(empty).unwrap().validation_error.as_deref(),
        Some("PSSM data must contain either scores or frequency ratios")
    );
    let mut short = read();
    short.pssm.num_columns = 29;
    assert!(psi_search_of(short)
        .unwrap_err()
        .to_string()
        .contains("dimensions (28 rows, 29 columns)"));
}

fn rnd() -> PssmWithParameters {
    crate::blastinput::pssm_text::read_pssm_with_parameters(&fixture(
        "db_search/db.ps.asn.rnd.6.asn",
    ))
    .unwrap()
}

// The review of Q4c, S1: the columns' largest scores of a PSSM must add up to less than 2^30,
// NCBI's `MININT` margin (`Int4` sums past it can overflow); with every score `INT2_MAX`,
// 32,769 columns pass (2^30 - 1) and 32,770 are refused.
#[test]
fn a_pssm_whose_scores_can_add_up_to_2_30_is_refused() {
    use crate::report::serial::pssm::SeqData;
    let with_columns = |columns: usize| {
        let mut pssm = rnd();
        pssm.pssm.num_columns = columns as i32;
        pssm.pssm.intermediate_data = None;
        pssm.pssm.final_data.as_mut().unwrap().scores = vec![i32::from(i16::MAX); 28 * columns];
        let inst = &mut pssm
            .pssm
            .query
            .as_mut()
            .and_then(SeqEntry::seq_mut)
            .unwrap()
            .inst;
        inst.length = Some(columns as u32);
        match inst.seq_data.as_mut().unwrap() {
            SeqData::Iupacaa(letters) | SeqData::Ncbieaa(letters) | SeqData::Ncbistdaa(letters) => {
                *letters = vec![letters[0]; columns];
            }
            _ => unreachable!("the test's PSSM is NCBIstdaa"),
        }
        psi_search_of(pssm)
    };
    assert!(with_columns(32_769).is_ok());
    let error = with_columns(32_770).unwrap_err().to_string();
    assert!(
        error.contains("largest scores add up to 1073774590"),
        "{error}"
    );
}

// The review of PS4a, S1: titles that NCBI's Query= line and its tabular query id read
// differently (several titles; white space at a title's ends) are refused; one title is read.
#[test]
fn query_titles_that_ncbi_reads_two_ways_are_refused() {
    let with_titles = |titles: &[&[u8]]| {
        let mut pssm = rnd();
        pssm.pssm
            .query
            .as_mut()
            .and_then(SeqEntry::seq_mut)
            .unwrap()
            .descr = Some(
            titles
                .iter()
                .map(|title| Seqdesc::Title(title.to_vec()))
                .collect(),
        );
        psi_search_of(pssm).map(|search| search.query.title)
    };
    assert_eq!(with_titles(&[b"rnd random"]).unwrap(), b"rnd random");
    assert_eq!(with_titles(&[b""]).unwrap(), b"");
    let error = with_titles(&[b"rnd", b"second"]).unwrap_err().to_string();
    assert!(error.contains("a PSSM query with 2 titles"), "{error}");
    for title in [&b" rnd"[..], b"rnd ", b"  sp1  two  words  "] {
        let error = with_titles(&[title]).unwrap_err().to_string();
        assert!(error.contains("starts or ends with white space"), "{error}");
    }
}

// The review of PS5b, B1: the warning label's title is `CBlastQuerySourceBioseqSet::GetTitle`'s
// (bioseq_extract_data_priv.cpp:345-372, cited at `psi_search_and_object`): the title without
// its trailing `.` and spaces, none with a `molinfo` descriptor (NCBI 2.17.0 probes of
// `PS/tmp/rev5b`, `d.*` and `m.molinfo*`); the record keeps the title for the reports.
#[test]
fn the_warning_title_is_the_query_sources() {
    // A descriptor other than a title, of the variant `name` (its value does not matter here).
    let other = |name: &str| {
        let index = crate::report::serial::asn_types::member_index(
            crate::report::serial::asn_types::table::named::SEQDESC,
            name,
        );
        Seqdesc::Other(Value::Choice(index as u16, Box::new(Value::Null)))
    };
    let warning = |title: &[u8], others: Vec<&'static str>| {
        let mut pssm = rnd();
        let mut descr = vec![Seqdesc::Title(title.to_vec())];
        descr.extend(others.into_iter().map(|name| match name {
            "Seqdesc" => other("comment"),
            name => other(name),
        }));
        pssm.pssm
            .query
            .as_mut()
            .and_then(SeqEntry::seq_mut)
            .unwrap()
            .descr = Some(descr);
        let (search, _) = psi_search_and_object(pssm, PssmReader::PSIBLAST).unwrap();
        assert_eq!(search.query.title, title);
        search.warning_title.unwrap()
    };
    assert_eq!(warning(b"rnd random", vec![]), b"rnd random");
    assert_eq!(warning(b"rnd random.", vec![]), b"rnd random");
    assert_eq!(warning(b"rnd. random.. . .", vec![]), b"rnd. random");
    assert_eq!(warning(b".", vec![]), b"");
    assert_eq!(warning(b"rnd random", vec!["Seqdesc"]), b"rnd random");
    assert_eq!(warning(b"rnd random", vec!["Seqdesc", "molinfo"]), b"");
}

// Line RG D3-1 (`om/` cases i01, m03): annotations without objects give NCBI's object manager
// nothing to index; the search keeps no message (`annot_index` has the other cases).
#[test]
fn empty_annotations_give_no_messages() {
    use crate::report::serial::asn_types::{member_index, members, table::named};
    let data = member_index(named::SEQ_ANNOT, "data");
    let kinds = members(named::SEQ_ANNOT)[data].ty;
    let with_annot = |kinds_of_data: &[&str]| {
        let mut pssm = rnd();
        pssm.pssm
            .query
            .as_mut()
            .and_then(SeqEntry::seq_mut)
            .unwrap()
            .annot = Some(
            kinds_of_data
                .iter()
                .map(|kind| {
                    let variant = member_index(kinds, kind) as u16;
                    Value::Sequence(vec![(
                        data as u16,
                        Value::Choice(variant, Box::new(Value::List(Vec::new()))),
                    )])
                })
                .collect(),
        );
        psi_search_of(pssm)
    };
    for annot in [
        &[][..],
        &["ids"],
        &["ftable", "align", "graph", "locs", "ids"],
    ] {
        let search = with_annot(annot).unwrap();
        assert_eq!(search.annot, Default::default(), "{annot:?}");
        assert!(!search.defline_reads_annot);
    }
}

// The review of PS4a, S2 and N3: the record's local id is the Seq-id's `GetSeqIdString()`
// (NCBI 2.17.0: `local id 5` shows as `5`, `local id -3` as `-3`, in `qseqid` and the warnings).
// The record's Seq-id is the Bioseq's, with or without `-parse_deflines` (review of PS15b-1, B2:
// NCBI's Seq-aligns of `local id 5` write `local id 5`, not `local str "5"`).
#[test]
fn the_local_id_is_the_seq_id_string() {
    for (id, shown) in [
        (ObjectId::Str(b"myq".to_vec()), "myq"),
        (ObjectId::Id(5), "5"),
        (ObjectId::Id(-3), "-3"),
        (ObjectId::Id(0), "0"),
        (ObjectId::Id(i32::MAX), "2147483647"),
    ] {
        let mut pssm = rnd();
        pssm.pssm
            .query
            .as_mut()
            .and_then(SeqEntry::seq_mut)
            .unwrap()
            .ids = vec![ParsedSeqId::Local(id.clone())];
        let search = psi_search_of(pssm).unwrap();
        assert_eq!(search.query.local_id, shown);
        assert!(search.query.believed().is_none());
        assert_eq!(search.query.warning_id(), shown.as_bytes());
        assert!(search.query_record(true).believed().is_some());
        let parsed = ParsedSeqId::Local(id);
        for record in [search.query_record(false), search.query_record(true)] {
            // The Seq-id object that the Seq-aligns write (`SeqId::write`).
            let written = match record.seq_id() {
                crate::report::seqid::SeqId::Db(id) => id.clone(),
                crate::report::seqid::SeqId::Local(text) => {
                    ParsedSeqId::Local(ObjectId::Str(text.as_bytes().to_vec()))
                }
                other => panic!("{other:?}"),
            };
            assert_eq!(written, parsed);
        }
        assert_eq!(search.query_seq_loc_id(), Some(&parsed));
    }
}

// Line RG Q1b (`Q/SURVEY_Q1b.md`, NCBI BLAST+ 2.17.0): a PSSM query whose Seq-id is not local
// shows its ids with or without `-parse_deflines` (`Query= gi|4504347 ...`, `qseqid`
// `gi|4504347`, the warnings' label `4504347`), and its Seq-aligns take its first id
// (`m_gr`: of `gi 4504347, other NP_000549.1` the gi, while `Query=` shows `NP_000549.1`).
#[test]
fn a_query_with_other_seq_ids_shows_them() {
    use crate::seqid::model::TextseqId;
    let refseq = ParsedSeqId::Other(TextseqId {
        accession: Some(b"NP_000549".to_vec()),
        version: Some(1),
        ..TextseqId::default()
    });
    for (ids, line, shown, warning) in [
        (
            vec![ParsedSeqId::Gi(4504347)],
            "gi|4504347",
            "4504347",
            "4504347",
        ),
        (
            vec![ParsedSeqId::Gi(4504347), refseq.clone()],
            "NP_000549.1",
            "NP_000549.1",
            "4504347",
        ),
        (
            vec![refseq.clone(), ParsedSeqId::Gi(4504347)],
            "NP_000549.1",
            "NP_000549.1",
            "NP_000549",
        ),
    ] {
        let mut pssm = rnd();
        pssm.pssm
            .query
            .as_mut()
            .and_then(SeqEntry::seq_mut)
            .unwrap()
            .ids = ids.clone();
        let search = psi_search_of(pssm).unwrap();
        for record in [search.query_record(false), search.query_record(true)] {
            let believed = record.believed().expect("the ids are shown");
            assert_eq!(record.local_id, line);
            assert_eq!(record.shown_id(), shown.as_bytes());
            assert_eq!(record.warning_id(), warning.as_bytes());
            assert_eq!(believed.best, ids[0]);
            assert_eq!(believed.ids, ids);
        }
        assert_eq!(search.query_seq_loc_id(), Some(&ids[0]));
    }
}

// `Q/SURVEY_Q1b.md` of line RG: a query without a Seq-id, a Geninfo import id (NCBI: `Failed to
// resolve SeqId: gim|77`), a database ordinal, a local id among other ids (shown as no FASTA
// query is), an empty local id, a byte outside the visible range and a Textseq-id with neither
// an accession nor a name (NCBI: a database error) are refused.
#[test]
fn query_seq_ids_that_are_refused() {
    use crate::seqid::model::{Dbtag, GiimportId};
    for (ids, refusal) in [
        (vec![], "without a Seq-id"),
        (
            vec![ParsedSeqId::Giim(GiimportId {
                id: 77,
                db: None,
                release: None,
            })],
            "is a Geninfo import id",
        ),
        (
            vec![ParsedSeqId::General(Dbtag {
                db: b"BL_ORD_ID".to_vec(),
                tag: ObjectId::Id(5),
            })],
            "names a database ordinal",
        ),
        (
            vec![
                ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec())),
                ParsedSeqId::Gi(4504347),
            ],
            "a local Seq-id among other Seq-ids",
        ),
        (
            vec![ParsedSeqId::Local(ObjectId::Str(Vec::new()))],
            "an empty local id",
        ),
        (
            vec![ParsedSeqId::Other(crate::seqid::model::TextseqId {
                release: Some(b"x".to_vec()),
                ..Default::default()
            })],
            "neither an accession nor a name",
        ),
    ] {
        let mut pssm = rnd();
        pssm.pssm
            .query
            .as_mut()
            .and_then(SeqEntry::seq_mut)
            .unwrap()
            .ids = ids;
        let error = psi_search_of(pssm).unwrap_err().to_string();
        assert!(error.contains(refusal), "{error}");
    }
}

// The review of PS4a, N8: a score outside one letter pair's range (`Int2`) is refused;
// `BLAST_SCORE_MIN` and `BLAST_SCORE_MAX` are scores.
#[test]
fn scores_outside_int2_are_refused() {
    let with_score = |score: i32| {
        let mut pssm = rnd();
        pssm.pssm.final_data.as_mut().unwrap().scores[5] = score;
        psi_search_of(pssm)
    };
    for score in [BLAST_SCORE_MIN, BLAST_SCORE_MAX] {
        assert_eq!(with_score(score).unwrap().pssm.rows[0][5], score);
    }
    for score in [BLAST_SCORE_MIN - 1, BLAST_SCORE_MAX + 1, i32::MAX] {
        let error = with_score(score).unwrap_err().to_string();
        assert!(
            error.contains(&format!("a PSSM score of {score}, outside")),
            "{error}"
        );
    }
    // A PSSM that NCBI's check refuses keeps NCBI's error (after the prolog).
    let mut scaled = rnd();
    let data = scaled.pssm.final_data.as_mut().unwrap();
    data.scaling_factor = Some(2);
    data.scores[5] = i32::MAX;
    assert!(psi_search_of(scaled).unwrap().validation_error.is_some());
}

/// `rnd()` (30 letters of NCBIstdaa) with its query data replaced by `data(letters, codes)`.
fn rnd_with_query(data: impl Fn(Vec<u8>, Vec<u8>) -> SeqData) -> Result<PsiSearch> {
    let mut pssm = rnd();
    let inst = &mut pssm
        .pssm
        .query
        .as_mut()
        .and_then(SeqEntry::seq_mut)
        .unwrap()
        .inst;
    let Some(SeqData::Ncbistdaa(codes)) = inst.seq_data.clone() else {
        unreachable!("the rnd PSSM's query is NCBIstdaa");
    };
    let letters = codes
        .iter()
        .map(|&code| NCBISTDAA_TO_AMINOACID[usize::from(code)])
        .collect();
    inst.seq_data = Some(data(letters, codes));
    psi_search_of(pssm)
}

fn put(mut bytes: Vec<u8>, at: usize, byte: u8) -> Vec<u8> {
    bytes[at] = byte;
    bytes
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:894-917 and 632-651 (cited at `query_residues` and
// `query_warnings::invalid_residues_error`): a residue without an NCBIstdaa code stops the
// search at its query set-up, after the strategy's request and psiblast's round-1 files
// (NCBI BLAST+ 2.17.0, M4a-1c E probes `*_star_iupac.inv.*`, `old_ps_bad_*`).
#[test]
fn a_residue_without_a_code_stops_the_search_at_its_set_up() {
    let iupacaa: [u8; 4] = [b'*', b'#', 0xE9, b'-'];
    let ncbieaa: [u8; 3] = [b'#', 0xE9, b'a'];
    let ncbistdaa: [u8; 3] = [28, 30, 255];
    let cases = iupacaa
        .iter()
        .map(|&byte| rnd_with_query(|letters, _| SeqData::Iupacaa(put(letters, 26, byte))))
        .chain(
            ncbieaa
                .iter()
                .map(|&byte| rnd_with_query(|letters, _| SeqData::Ncbieaa(put(letters, 26, byte)))),
        )
        .chain(
            ncbistdaa
                .iter()
                .map(|&code| rnd_with_query(|_, codes| SeqData::Ncbistdaa(put(codes, 26, code)))),
        );
    for search in cases {
        let search = search.unwrap();
        assert_eq!(
            search.setup_error.as_deref(),
            Some("Warning: Invalid residues found at positions 26 ")
        );
        assert!(search.validation_error.is_none() && search.scores_error.is_none());
        let stop = search.stop_before_search(false, true).unwrap();
        assert!(matches!(
            stop,
            crate::algorithm::pssm_search::PssmQueryStop::Setup(_)
        ));
        assert!(stop.after_request() && stop.after_round_files());
        assert!(stop
            .error()
            .to_string()
            .starts_with("LOSAT engine error: Warning: Invalid residues found at positions 26 "));
    }
    let all = rnd_with_query(|letters, _| SeqData::Iupacaa(vec![b'#'; letters.len()])).unwrap();
    assert_eq!(
        all.setup_error.as_deref(),
        Some(
            "Warning: Invalid residues found at positions 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, \
                 11, 12, 13, 14, 15, 16, 17, 18, 19,... (only first 20 shown) "
        )
    );
}

// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1748-1781: `*` of ncbieaa is the stop residue (25), which NCBI
// searches (M4a-1c E probes `*_star_eaa`); the gap residue (`-`, the code 0 of ncbistdaa) is
// a residue too, which NCBI searches and its reports show as `-` and count as a residue
// (M4a-1c E and LG probes `*_dash_eaa`, `*_std0*`): the record keeps it as `-`.
#[test]
fn a_stop_letter_and_a_gap_residue_are_searched() {
    let star = rnd_with_query(|letters, _| SeqData::Ncbieaa(put(letters, 26, b'*'))).unwrap();
    assert!(star.setup_error.is_none() && star.stop_before_search(false, true).is_none());
    assert_eq!(star.query.seq()[26], b'*');
    for (gap, position) in [
        (
            rnd_with_query(|letters, _| SeqData::Ncbieaa(put(letters, 26, b'-'))),
            26,
        ),
        (
            rnd_with_query(|_, codes| SeqData::Ncbistdaa(put(codes, 26, 0))),
            26,
        ),
        (
            rnd_with_query(|_, codes| SeqData::Ncbistdaa(put(codes, 0, 0))),
            0,
        ),
    ] {
        let gap = gap.unwrap();
        assert!(gap.setup_error.is_none() && gap.stop_before_search(false, true).is_none());
        assert_eq!(gap.query.seq()[position], b'-');
        assert!(gap.bioseq_letters.is_none());
    }
    // With an invalid residue too, NCBI's stop comes first.
    let stopped =
        rnd_with_query(|letters, _| SeqData::Ncbieaa(put(put(letters, 26, b'-'), 3, b'#')))
            .unwrap();
    assert_eq!(
        stopped.setup_error.as_deref(),
        Some("Warning: Invalid residues found at positions 3 ")
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:466-473: a PSSM of frequency ratios only sets its query up as it
// computes the scores, so the stop comes there (before the request; no PSSM files).
#[test]
fn a_ratios_only_pssm_stops_at_its_scores_on_a_residue_without_a_code() {
    let pssm = crate::blastinput::pssm_text::read_pssm_with_parameters(&gunzip(
        "pssm/ps_ratios_std30.asn.gz",
    ))
    .unwrap();
    let search = psi_search_of(pssm).unwrap();
    assert_eq!(
        search.scores_error.as_deref(),
        Some("Warning: Invalid residues found at positions 5 ")
    );
    assert!(search.setup_error.is_none());
    let stop = search.stop_before_search(false, true).unwrap();
    assert!(stop.after_request() && !stop.after_round_files());
}

// NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:530-547,566-569 (cited at `query_residues`) and c++/src/objtools/align_format/tabular.cpp:996-1022
// (cited at `recount_shown_identities`): NCBI searches lower-case iupacaa as its capitals and
// its reports read the Bioseq's own letters (`PsiSearch::bioseq_letters`); a query that NCBI
// stops keeps NCBI's stop.
#[test]
fn lower_case_iupacaa_is_searched_as_capitals_and_shown_as_it_is() {
    let lower =
        rnd_with_query(|letters, _| SeqData::Iupacaa(letters.to_ascii_lowercase())).unwrap();
    assert!(lower.setup_error.is_none() && lower.stop_before_search(false, true).is_none());
    let bioseq = lower.bioseq_letters.as_deref().unwrap();
    assert_eq!(bioseq.to_ascii_uppercase(), lower.query.seq());
    assert!(bioseq.iter().all(|letter| !letter.is_ascii_uppercase()));
    let capitals = rnd_with_query(|letters, _| SeqData::Iupacaa(letters)).unwrap();
    assert!(capitals.bioseq_letters.is_none());
    assert_eq!(capitals.query.seq(), lower.query.seq());
    let stopped =
        rnd_with_query(|letters, _| SeqData::Iupacaa(put(letters.to_ascii_lowercase(), 3, b'*')))
            .unwrap();
    assert_eq!(
        stopped.setup_error.as_deref(),
        Some("Warning: Invalid residues found at positions 3 ")
    );
}

#[test]
fn scores_fill_28_rows_by_column_or_by_row() {
    // Two rows (A, B) of three positions, column by column.
    let rows = convert_2_matrix(&[1, 2, 3, 4, 5, 6], 2, 3, false, BLAST_SCORE_MIN);
    assert_eq!(rows.len(), 3);
    assert_eq!((rows[0][0], rows[0][1], rows[2][1]), (1, 2, 6));
    assert_eq!(rows[1][2], BLAST_SCORE_MIN);
    // Fewer values than cells: the cells past them keep the fill.
    let rows = convert_2_matrix(&[1, 2, 3], 2, 3, false, BLAST_SCORE_MIN);
    assert_eq!((rows[0][0], rows[0][1], rows[1][0]), (1, 2, 3));
    assert_eq!((rows[1][1], rows[2][0]), (BLAST_SCORE_MIN, BLAST_SCORE_MIN));
    assert_eq!(
        convert_2_matrix(&[1, 2], 0, 2, false, 0),
        vec![[0; BLASTAA_SIZE]; 2]
    );
    // The same values row by row.
    let rows = convert_2_matrix(&[1, 3, 5, 2, 4, 6], 2, 3, true, BLAST_SCORE_MIN);
    assert_eq!((rows[0][0], rows[0][1], rows[2][1]), (1, 2, 6));
    assert_eq!(rows[2][27], BLAST_SCORE_MIN);
}

// Line RG D3c-3a, finding 2: NCBI BLAST+ 2.17.0's psiblast ends by a signal on a binary PSSM whose
// query Seq-inst length (a Uint4) is not the count of its letters: `02 01 ff` (255 for 30 letters:
// SIGSEGV) and `02 05 00 ff ff ff ff` (2^32 - 1: SIGABRT); inputs kept in `D/d3/c3a/crash`. LOSAT
// refuses the query as it refuses that text PSSM.
#[test]
fn a_binary_query_whose_length_is_not_its_letters_is_refused() {
    let base = fixture("pssm/ps9_ber_base.asnb");
    let length = b"\xa2\x80\x02\x01\x1e\x00\x00";
    let at = base
        .windows(length.len())
        .position(|bytes| bytes == length)
        .unwrap();
    let read = crate::blastinput::pssm_text::read_pssm_binary(&base);
    assert!(psi_search_and_object(read.result.unwrap(), PssmReader::PSIBLAST).is_ok());
    for value in [&b"\x02\x01\xff"[..], b"\x02\x05\x00\xff\xff\xff\xff"] {
        let bytes = [&base[..at + 2], value, &base[at + 5..]].concat();
        let read = crate::blastinput::pssm_text::read_pssm_binary(&bytes);
        let Err(error) = psi_search_and_object(read.result.unwrap(), PssmReader::PSIBLAST) else {
            panic!("read");
        };
        assert_eq!(
            error.to_string(),
            "a PSSM query whose letters are missing or of another length than its Seq-inst is \
             not supported by LOSAT's PSIBLAST"
        );
    }
}

// Line RG D3c-3a, finding 3: a binary PSSM cut before the last byte of its root's end-of-contents
// ends inside the root class's frames, where NCBI's stack passes the `CEofException` on: psiblast
// ends without a word, exit status 0.
#[test]
fn a_binary_pssm_cut_in_its_last_end_of_contents_ends_silently() {
    let base = fixture("pssm/ps9_ber_base.asnb");
    let read = crate::blastinput::pssm_text::read_pssm_binary(&base[..base.len() - 1]);
    assert!(matches!(
        read.result,
        Err(PssmReadError::Ncbi(ReadError::Eof))
    ));
    assert_eq!(outcome(&base[..base.len() - 1]), "end");
}

// Line RG D3c-3a: a string whose length passes the end of a binary PSSM, as NCBI's reader reserves
// it before it reads the bytes (`CIStreamBuffer::GetChars`): 1 byte past the end, the end of the
// file's report ("Unrecognized format for PSSM", rc 1); 2^56 + 10 bytes, `std::bad_alloc` ("BLAST
// ran out of memory", rc 4); 1 GiB past the end, a reservation that succeeds or fails by the
// host's memory (NCBI's probe host: "Unrecognized format for PSSM"), which LOSAT refuses
// (refusals.md, PD-LOSAT-NCBI-DEFECTS). Inputs: the rows `db.ber.ps.past_end1`, `len88_huge` and
// `past_end1g` of `tests/db_search_manifest/psiblast.ber.tsv`.
#[test]
fn a_binary_string_past_the_end_is_reserved_as_ncbi_reserves_it() {
    use crate::blastinput::pssm_text::read_pssm_binary;
    assert_eq!(
        outcome(&fixture("pssm/ps9_ber_past_end1.asnb")),
        "Unrecognized format"
    );
    let read = read_pssm_binary(&fixture("pssm/ps9_ber_len88_huge.asnb"));
    assert!(matches!(
        read.result,
        Err(PssmReadError::Ncbi(ReadError::OutOfMemory))
    ));
    let read = read_pssm_binary(&fixture("pssm/ps9_ber_past_end1g.asnb"));
    let what = match read.result {
        Err(PssmReadError::Refused(what)) => what,
        other => panic!("{:?}", other.err()),
    };
    assert!(what.starts_with("a string whose length (1073741"), "{what}");
}

// Line RG D3c-3a, review B1: a byte block (the query's ncbistdaa OCTET STRING, or a BIT STRING)
// whose length passes the end of a binary PSSM, as NCBI's `vector<char>` reader reserves the
// block's length before its bytes: 2^62 and 2^63 - 1 bytes, `std::bad_alloc` ("BLAST ran out of
// memory", rc 4); 2^63 bytes and a BIT STRING of length 0 (2^64 - 1 bytes), `std::length_error`
// ("Error: vector::reserve", rc 255), each after the post (803.4) unended byte block; 1 GiB, a
// reservation that succeeds or fails by the host's memory, which LOSAT refuses; a BIT STRING of
// the letters is read as them. Inputs: the rows `db.ber.ps.block_*` of
// `tests/db_search_manifest/psiblast.ber.tsv`.
#[test]
fn a_binary_byte_block_past_the_end_is_reserved_as_ncbi_reserves_it() {
    use crate::blastinput::pssm_text::read_pssm_binary;
    use crate::report::serial::asn_types::Post;
    for (case, error) in [
        ("2p62", ReadError::OutOfMemory),
        ("2p63m1", ReadError::OutOfMemory),
        ("2p63", ReadError::Std("vector::reserve")),
        ("bit0", ReadError::Std("vector::reserve")),
    ] {
        let read = read_pssm_binary(&fixture(&format!("pssm/ps9_ber_block_{case}.asnb")));
        assert_eq!(read.posts, vec![Post::UnendedByteBlock], "{case}");
        assert_eq!(
            read.result.err(),
            Some(PssmReadError::Ncbi(error)),
            "{case}"
        );
    }
    let read = read_pssm_binary(&fixture("pssm/ps9_ber_block_1g.asnb"));
    assert!(read.posts.is_empty());
    let what = match read.result {
        Err(PssmReadError::Refused(what)) => what,
        other => panic!("{:?}", other.err()),
    };
    assert!(
        what.starts_with(
            "an OCTET STRING or BIT STRING whose length (1073741824 bytes) passes the end of the \
             input: NCBI's reader allocates 1073741824 bytes for it"
        ),
        "{what}"
    );
    let bit = read_pssm_binary(&fixture("pssm/ps9_ber_block_bit.asnb"));
    let base = read_pssm_binary(&fixture("pssm/ps9_ber_base.asnb"));
    assert!(bit.posts.is_empty());
    assert_eq!(bit.result, base.result);
}

// Line RG D3c-3a: Karlin-Altschul values that are not finite numbers, in ASN.1 text and in its
// binary twins (`ps9_ber_real_*`: `lambda` of PLUS-INFINITY `09 01 40`, MINUS-INFINITY `09 01 41`,
// NOT-A-NUMBER `09 01 42`, the NR text `INF`). NCBI BLAST+ 2.17.0 goes on with each (probes
// `D/d3/c3a/probe/realtext`, rows `db.ber.*.real_*`); LOSAT refuses where its search does not yet
// compute NCBI's hits, E-values and bit scores, and searches where it does (gapped `kappa` and `h`
// of PLUS-INFINITY; tblastn's `lambdaUngapped` of PLUS-INFINITY).
#[test]
fn a_karlin_altschul_value_that_is_not_finite_is_refused_where_losat_does_not_compute_it() {
    use crate::blastinput::pssm_text::{read_pssm, read_pssm_binary};
    let tblastn = PssmReader {
        program: "TBLASTN",
        matrix_name: "BLOSUM62",
        strategy: false,
        database: true,
        archive: None,
    };
    let refusal = |kind: &str, name: &str, what: &str, program: &str| {
        format!(
            "a PSSM whose {kind} Karlin-Altschul parameter {name} is {what} (LOSAT does not yet \
             compute the hits, E-values and bit scores that NCBI BLAST+ 2.17.0 computes from it) \
             is not supported by LOSAT's {program}"
        )
    };
    let search = |pssm: PssmWithParameters, reader: PssmReader<'static>| {
        psi_search_and_object(pssm, reader)
            .map(|_| ())
            .map_err(|error| error.to_string())
    };
    let base = String::from_utf8(fixture("db_search/db.ps.asn.rnd.6.asn")).unwrap();
    let text = |name: &str, value: &str| {
        let at = base.find(&format!("      {name} ")).unwrap();
        let end = at + base[at..].find('\n').unwrap();
        let comma = if base[at..end].ends_with(',') {
            ","
        } else {
            ""
        };
        let edited = format!("{}      {name} {value}{comma}{}", &base[..at], &base[end..]);
        read_pssm(edited.as_bytes()).result.unwrap()
    };
    for (name, value, what) in [
        ("lambda", "PLUS-INFINITY", "infinite"),
        ("lambda", "NOT-A-NUMBER", "not a number"),
        ("kappa", "NOT-A-NUMBER", "not a number"),
        ("h", "NOT-A-NUMBER", "not a number"),
    ] {
        for (reader, program) in [(PssmReader::PSIBLAST, "PSIBLAST"), (tblastn, "TBLASTN")] {
            assert_eq!(
                search(text(name, value), reader),
                Err(refusal("gapped", name, what, program)),
                "{name} {value}"
            );
        }
    }
    assert_eq!(
        search(
            text("lambdaUngapped", "PLUS-INFINITY"),
            PssmReader::PSIBLAST
        ),
        Err(refusal(
            "ungapped",
            "lambdaUngapped",
            "infinite",
            "PSIBLAST"
        ))
    );
    for (name, value) in [
        ("kappa", "PLUS-INFINITY"),
        ("h", "PLUS-INFINITY"),
        ("kappaUngapped", "PLUS-INFINITY"),
        ("hUngapped", "PLUS-INFINITY"),
    ] {
        assert_eq!(
            search(text(name, value), PssmReader::PSIBLAST),
            Ok(()),
            "{name}"
        );
    }
    assert_eq!(
        search(text("lambdaUngapped", "PLUS-INFINITY"), tblastn),
        Ok(())
    );
    // An ungapped value that is not a number makes NCBI's ungapped block invalid first: the
    // engine's error of a query without Karlin-Altschul parameters, not this refusal.
    let mut nan = text("lambda", "NOT-A-NUMBER");
    let ungapped = text("lambdaUngapped", "NOT-A-NUMBER");
    nan.pssm.final_data.as_mut().unwrap().lambda_ungapped =
        ungapped.pssm.final_data.unwrap().lambda_ungapped;
    assert_eq!(search(nan, PssmReader::PSIBLAST), Ok(()));
    // The binary twins: the same values, the same answers; MINUS-INFINITY is negative.
    let binary = |name: &str| {
        read_pssm_binary(&fixture(&format!("pssm/ps9_ber_{name}.asnb")))
            .result
            .unwrap()
    };
    for (name, what) in [
        ("real_pinf", "infinite"),
        ("real_inf_text", "infinite"),
        ("real_nan", "not a number"),
    ] {
        for (reader, program) in [(PssmReader::PSIBLAST, "PSIBLAST"), (tblastn, "TBLASTN")] {
            assert_eq!(
                search(binary(name), reader),
                Err(refusal("gapped", "lambda", what, program)),
                "{name}"
            );
        }
    }
    assert_eq!(
        search(binary("real_ninf"), PssmReader::PSIBLAST),
        Err(
            "a PSSM whose gapped Karlin-Altschul parameter lambda is negative (NCBI BLAST+ 2.17.0 \
             goes on with a Karlin-Altschul block that it holds invalid) is not supported by \
             LOSAT's PSIBLAST"
                .to_string()
        )
    );
    assert_eq!(search(binary("real_negzero"), PssmReader::PSIBLAST), Ok(()));
}
