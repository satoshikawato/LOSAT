//! Unit tests of the archive reader (`super`): every archive among LOSAT's fixtures read and
//! written back as read; the objects of the inputs that NCBI BLAST+ 2.17.0's blast_formatter was
//! given in the probes of line RG-F3 (`F/f3/ncbi/C`, SURVEY_F3_NCBI §2.3: P3.077-124), with the
//! objects delivered and how the stream ended; and what the typed model reads of the archives
//! that NCBI's programs wrote (`tests/fixtures/archive_read`), against the values of their text.

use std::ops::Range;
use std::path::Path;

use super::options::archive_request;
use super::request::{ParamSet, Queries, Request, Subject};
use super::results::{KarlinBlocks, MaskRange, Results, SearchStats};
use super::{seq_align_id, seq_loc_id, Archive, ArchiveReader, Message, Search, Severity, Stop};
use crate::report::serial::asn_text::AsnTextEncoder;
use crate::report::serial::asn_text_reader::real_of_parts;
use crate::report::serial::asn_types::stream::{guess_encoding, Encoding, ObjectStream};
use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::view::V;
use crate::report::serial::asn_types::write::write_document;
use crate::report::serial::asn_types::{read_document, ReadError, Value};
use crate::report::serial::encoder::Encoder;
use crate::seqid::ParsedSeqId;

/// A fixture of this test, by name (`LOSAT/tests/fixtures/archive_read`).
fn fixture(name: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/archive_read")
        .join(name);
    std::fs::read(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()))
}

/// The archive objects of a text until `EndOfData`, and the reader's error that ended them.
fn archives(text: &[u8]) -> (Vec<Archive>, Option<ReadError>) {
    let guessed = guess_encoding(text);
    assert_eq!(guessed.encoding.as_ref().ok(), Some(&Encoding::Text));
    let mut reader = ArchiveReader::new(guessed.into_text());
    let mut out = Vec::new();
    while let Some((archive, posts)) = reader.next_archive() {
        assert!(posts.is_empty(), "{posts:?}");
        match archive {
            Ok(archive) => out.push(archive),
            Err(error) => return (out, Some(error)),
        }
    }
    (out, None)
}

/// The archive objects of a fixture, read whole.
fn read_fixture(name: &str) -> Vec<Archive> {
    let (archives, error) = archives(&fixture(name));
    assert_eq!(error, None, "{name}");
    archives
}

/// The search of an archive object.
fn search_of(archive: &Archive) -> &Search {
    archive
        .search
        .as_ref()
        .expect("a search")
        .as_ref()
        .expect("not an error archive")
}

/// `CSeq_id::AsFastaString` of an id.
fn label(id: &ParsedSeqId) -> String {
    let mut out = Vec::new();
    crate::seqid::label::fasta(id, &mut out);
    String::from_utf8(out).expect("ASCII")
}

/// The value of a text of the type `root`, read whole without a message.
fn value(root: u16, text: &str) -> Value {
    let (value, posts) = read_document(text.as_bytes(), root);
    assert!(posts.is_empty(), "{posts:?}");
    value.unwrap_or_else(|error| panic!("{text}: {error:?}"))
}

/// The place, class and message of NCBI's exception that a model stops with.
fn exception(stop: Stop) -> (&'static str, u32, &'static str, String) {
    match stop {
        Stop::Exception(e) => (
            e.site().file,
            e.site().line,
            e.kind(),
            String::from_utf8_lossy(e.get_msg()).into_owned(),
        ),
        other => panic!("not an exception: {other:?}"),
    }
}

/// The archives of rpsbproc's tests that are copies of NCBI's edited by hand
/// (`fixtures/rpsbproc/README.md`: a member id in brackets, `length -0`, an `org` descriptor, no
/// `inst.length`, no query `seq-data`, a `delta` query): rpsbproc's own cases, not NCBI's writer's
/// layout (`length -0` is not even read).
const EDITED: [&str; 6] = [
    "fixtures/rpsbproc/root/bracket_member.asn",
    "fixtures/rpsbproc/root/neg0_length.asn",
    "fixtures/rpsbproc/root/org_descr.asn",
    "fixtures/rpsbproc/root/nolength.asn",
    "fixtures/rpsbproc/arch/q_nodata.asn",
    "fixtures/rpsbproc/arch/mfirst_delta.asn",
];

/// Every archive among LOSAT's fixtures (the files under `LOSAT/tests` that start with an archive
/// header: NCBI's `-outfmt 11` output and the files that LOSAT's tests freeze as NCBI's) is read
/// whole, object by object, and its values written back are its bytes. Left out: the malformed
/// archives of rpsbproc's tests (`fixtures/rpsbproc/arch/bad`), its copies edited by hand
/// (`EDITED`), and the `.merged` files of `archive_replaced` (a run's stdout and stderr in one
/// file, with NCBI's messages inside the archive's text). An archive whose paths the tests
/// replaced by a placeholder (`has_placeholder`) is compared without the line breaks inside its
/// strings: NCBI broke those lines where the real path put them.
/// Whether a fixture holds a placeholder that its test replaces by a path (`@CWD@`,
/// `{db_dir}`: an upper-case name between `@`, a lower-case one between braces).
fn has_placeholder(text: &[u8]) -> bool {
    let between = |open: u8, close: u8, letter: fn(&u8) -> bool| {
        text.iter().enumerate().any(|(at, &b)| {
            b == open
                && text[at + 1..]
                    .iter()
                    .position(|&b| b == close)
                    .is_some_and(|len| {
                        let name = &text[at + 1..at + 1 + len];
                        !name.is_empty() && name.iter().all(|b| *b == b'_' || letter(b))
                    })
        })
    };
    between(b'@', b'@', u8::is_ascii_uppercase) || between(b'{', b'}', u8::is_ascii_lowercase)
}

/// The text without the line breaks inside its strings (a doubled quote toggles twice).
fn unwrapped(text: &[u8]) -> Vec<u8> {
    let mut inside = false;
    text.iter()
        .copied()
        .filter(|&b| {
            if b == b'"' {
                inside = !inside;
            }
            !(inside && b == b'\n')
        })
        .collect()
}

#[test]
fn every_archive_of_the_fixtures_is_written_back_as_read() {
    let tests = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
    let bad = tests.join("fixtures/rpsbproc/arch/bad");
    let mut todo = vec![tests];
    let (mut files, mut objects) = (0, 0);
    let mut failures = Vec::new();
    while let Some(dir) = todo.pop() {
        for entry in std::fs::read_dir(&dir).expect("a directory") {
            let path = entry.expect("a directory entry").path();
            if path.is_dir() {
                if path != bad {
                    todo.push(path);
                }
                continue;
            }
            if path
                .extension()
                .is_some_and(|ext| ext == "rs" || ext == "merged")
            {
                continue;
            }
            if EDITED.iter().any(|name| path.ends_with(name)) {
                continue;
            }
            let text = std::fs::read(&path).expect("a file");
            if !text.starts_with(b"Blast4-archive ::=") {
                continue;
            }
            files += 1;
            let mut stream = ObjectStream::new(&text[..], b"");
            let mut out = Vec::new();
            let mut error = None;
            while !stream.end_of_data() {
                let (value, posts) = stream.read(named::BLAST4_ARCHIVE);
                if !posts.is_empty() {
                    error = Some(format!("messages {posts:?}"));
                }
                match value {
                    Ok(value) => {
                        let mut enc = AsnTextEncoder::new(&mut out);
                        write_document(&mut enc, named::BLAST4_ARCHIVE, &value);
                        assert_eq!(enc.error(), None);
                        objects += 1;
                    }
                    Err(e) => {
                        error = Some(format!("{e:?}"));
                        break;
                    }
                }
            }
            let placeholder = has_placeholder(&text);
            if error.is_none()
                && out != text
                && (!placeholder || unwrapped(&out) != unwrapped(&text))
            {
                let at = out
                    .iter()
                    .zip(&text)
                    .position(|(a, b)| a != b)
                    .unwrap_or(out.len().min(text.len()));
                let around = |bytes: &[u8]| {
                    String::from_utf8_lossy(
                        &bytes[at.saturating_sub(40)..(at + 40).min(bytes.len())],
                    )
                    .into_owned()
                };
                error = Some(format!(
                    "written otherwise from byte {at}: {:?} / {:?}",
                    around(&text),
                    around(&out)
                ));
            }
            if let Some(error) = error {
                failures.push(format!("{}: {error}", path.display()));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {files} archives: {:#?}",
        failures.len(),
        &failures[..failures.len().min(10)]
    );
    assert!(
        files >= 10 && objects >= files,
        "{files} archives, {objects} objects"
    );
}

/// How the objects of an archive text end.
#[derive(Debug, PartialEq)]
enum End {
    /// `EndOfData`: blast_formatter ends with 0.
    Data,
    /// The reader's `CEofException`, which blast_formatter's `CATCH_ALL` ends without a message
    /// and with 0.
    Eof,
    /// The reader's `CSerialException`: "Invalid input format for BLAST Archive.", 1.
    Serial,
}

/// The archive objects read of `text` and how they end; the encoding guessed for a text that is
/// not ASN.1 text.
fn outcome(text: &[u8]) -> Result<(usize, End), Encoding> {
    let guessed = guess_encoding(text);
    match guessed.encoding.as_ref() {
        Ok(Encoding::Text) => {}
        Ok(other) => return Err(*other),
        Err(stop) => panic!("{stop:?}"),
    }
    let mut reader = ArchiveReader::new(guessed.into_text());
    let mut count = 0;
    while let Some((archive, _)) = reader.next_archive() {
        match archive {
            Ok(_) => count += 1,
            Err(ReadError::Eof) => return Ok((count, End::Eof)),
            Err(ReadError::Serial(_)) => return Ok((count, End::Serial)),
            Err(other) => panic!("{other:?}"),
        }
    }
    Ok((count, End::Data))
}

/// The inputs of NCBI's probes, made from one archive (`blastn -outfmt 11` of three queries,
/// `blastn_multi.asn`), the same archive as binary ASN.1 and XML, and a batch of three
/// archives (`blastn_batch.asn`): the archives LOSAT reads before the end of the data or an
/// error, and the error, are the archives NCBI formatted and the code it ended with.
#[test]
fn the_objects_of_a_stream_end_as_ncbis_did() {
    let text = fixture("blastn_multi.asn");
    let batch = fixture("blastn_batch.asn");
    let binary = fixture("blastn_multi.bin");
    let xml = fixture("blastn_multi.xml");
    let with = |before: &[u8], after: &[u8]| [before, &text[..], after].concat();
    let crlf: Vec<u8> = text
        .iter()
        .flat_map(|&c| if c == b'\n' { vec![b'\r', c] } else { vec![c] })
        .collect();
    use End::{Data, Eof, Serial};
    let cases: Vec<(&str, Vec<u8>, Result<(usize, End), Encoding>)> = vec![
        // 0, the archives formatted (P3.087, 089-092, 096-097, 115-116).
        ("text", text.clone(), Ok((1, Data))),
        ("CRLF", crlf, Ok((1, Data))),
        ("leading white space", with(b"\n\n  ", b""), Ok((1, Data))),
        (
            "a comment before",
            with(b"-- comment\n", b""),
            Ok((1, Data)),
        ),
        (
            "trailing white space",
            with(b"", b"\n\n  \n"),
            Ok((1, Data)),
        ),
        ("twice", with(b"", &text[..]), Ok((2, Data))),
        ("no last newline", text[..12628].to_vec(), Ok((1, Data))),
        ("batch cut after 1", batch[..9366].to_vec(), Ok((1, Data))),
        ("batch cut after 2", batch[..11323].to_vec(), Ok((2, Data))),
        // 0 after the archive, no message: a line "junk" is a header cut short (P3.088).
        ("junk after", with(b"", b"junk\n"), Ok((1, Eof))),
        // 1 after the archives read before the error (P3.090, 117-118).
        ("NUL after", with(b"", b"\0"), Ok((1, Serial))),
        ("batch cut in 2", batch[..10344].to_vec(), Ok((1, Serial))),
        ("batch cut in 3", batch[..13187].to_vec(), Ok((2, Serial))),
        // 1, nothing formatted (P3.093, 098-104, 121).
        ("first 20 bytes", text[..20].to_vec(), Ok((0, Serial))),
        ("cut in the request", text[..1800].to_vec(), Ok((0, Serial))),
        (
            "cut before the results",
            text[..3601].to_vec(),
            Ok((0, Serial)),
        ),
        ("cut in the results", text[..8115].to_vec(), Ok((0, Serial))),
        ("no last 2 bytes", text[..12627].to_vec(), Ok((0, Serial))),
        ("no last 10 bytes", text[..12619].to_vec(), Ok((0, Serial))),
        ("Seq-entry", fixture("bad/other_type"), Ok((0, Serial))),
        ("search strategy", fixture("bad/strategy"), Ok((0, Serial))),
        // 1, nothing formatted: the guess is no text (P3.077-086, 099).
        ("first byte", text[..1].to_vec(), Err(Encoding::Other)),
        ("empty", fixture("bad/empty"), Err(Encoding::Other)),
        ("newline", fixture("bad/newline"), Err(Encoding::Other)),
        ("spaces", fixture("bad/spaces"), Err(Encoding::Other)),
        ("hello world", fixture("bad/garbage"), Err(Encoding::Other)),
        ("JSON", fixture("bad/json"), Err(Encoding::Other)),
        ("bytes 0-255", fixture("bad/binjunk"), Err(Encoding::Other)),
        ("gzip", fixture("bad/gzip"), Err(Encoding::Other)),
        (
            "UTF-8 BOM",
            with(b"\xef\xbb\xbf", b""),
            Err(Encoding::Binary),
        ),
        // NCBI's binary and XML readers (none in LOSAT yet): formatted (P3.001-015), 1 cut at
        // 50 % and 90 % or without the last 20 bytes, 0 without its last byte (binary: nothing
        // formatted; XML: formatted), 1 after the archive with an element after it.
        ("binary", binary.clone(), Err(Encoding::Binary)),
        (
            "binary cut at 50 %",
            binary[..2461].to_vec(),
            Err(Encoding::Binary),
        ),
        (
            "binary cut at 90 %",
            binary[..4429].to_vec(),
            Err(Encoding::Binary),
        ),
        (
            "binary no last byte",
            binary[..4921].to_vec(),
            Err(Encoding::Binary),
        ),
        ("XML", xml.clone(), Err(Encoding::Xml)),
        ("XML cut at 50 %", xml[..25305].to_vec(), Err(Encoding::Xml)),
        (
            "XML no last byte",
            xml[..50610].to_vec(),
            Err(Encoding::Xml),
        ),
        (
            "XML no last 20 bytes",
            xml[..50591].to_vec(),
            Err(Encoding::Xml),
        ),
        (
            "XML element after",
            [&xml[..], b"<junk/>"].concat(),
            Err(Encoding::Xml),
        ),
    ];
    for (name, bytes, expected) in cases {
        assert_eq!(outcome(&bytes), expected, "{name}");
    }
}

/// `GetSeqAlignSets`: the alignments of each query, cut where the id of the first row changes;
/// a query without alignments before the run of the next one gets none (here the queries
/// without hits, and NCBI's own form when a query between two others has none).
#[test]
fn the_alignments_are_grouped_by_query() {
    let cases: Vec<(&str, Vec<Vec<Range<usize>>>)> = vec![
        ("rpsblast_mixnh.asn", vec![vec![0..0, 0..5]]),
        ("rpstblastn_mixnh.asn", vec![vec![0..4, 4..14]]),
        ("blastn_subject_mixnh.asn", vec![vec![0..0, 0..7, 7..8]]),
        ("blastn_multi.asn", vec![vec![0..7, 7..7, 7..8]]),
        ("blastn_lcase.asn", vec![vec![0..7, 7..7, 7..8]]),
        ("blastx_strand_plus.asn", vec![vec![0..0]]),
        ("blastn_batch.asn", vec![vec![0..7], vec![0..0], vec![0..1]]),
        (
            "warnings_batch.asn",
            vec![vec![0..1], vec![0..0], vec![0..1]],
        ),
        ("psiblast_it3.asn", vec![vec![0..4], vec![0..4]]),
        ("psiblast_in_pssm.asn", vec![vec![0..4]]),
    ];
    for (name, expected) in cases {
        let archives = read_fixture(name);
        let sets: Vec<_> = archives
            .iter()
            .map(|archive| {
                let search = search_of(archive);
                search.results.alignment_sets(&search.request.queries)
            })
            .collect();
        assert_eq!(
            sets,
            expected.into_iter().map(Ok::<_, Stop>).collect::<Vec<_>>(),
            "{name}"
        );
    }
    // Queries as a list of locations (`blastn -parse_deflines`; LOSAT's fixture of NCBI's
    // output): the id of each location.
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/db_search/db.pd.blastn.nl.loc.11.out");
    let (archives, error) = archives(&std::fs::read(path).expect("the fixture"));
    assert_eq!(error, None);
    let request = &search_of(&archives[0]).request;
    assert!(matches!(request.queries, Queries::SeqLocList(_)));
    assert_eq!(request.queries.num_queries(), 4);
    let ids: Vec<_> = request
        .queries
        .query_ids()
        .expect("ids")
        .into_iter()
        .map(|id| id.map(label))
        .collect();
    assert_eq!(
        ids,
        [
            "gi|42",
            "ref|NC_000913.3|",
            "emb|X12345.1|N3",
            "tpg|BK000004.1|"
        ]
        .map(|id| Some(id.to_string()))
    );
    let sets: Vec<_> = archives
        .iter()
        .map(|archive| {
            let search = search_of(archive);
            search.results.alignment_sets(&search.request.queries)
        })
        .collect();
    assert_eq!(
        sets,
        vec![Ok::<_, Stop>(vec![0..2, 2..3, 3..3, 3..3]), Ok(vec![0..0])]
    );
}

/// The request as stored: the program and service, the queries (their Bioseqs, ids,
/// titles, deflines, molecule, length and data), the subject sequences or database, and the
/// parameters in the order of the file.
#[test]
fn the_request_is_read_as_stored() {
    let archives = read_fixture("blastn_subject_mixnh.asn");
    let request = &search_of(&archives[0]).request;
    assert_eq!(
        (request.program.as_slice(), request.service.as_slice()),
        (&b"blastn"[..], &b"megablast"[..])
    );
    assert!(!request.db_is_protein());
    assert_eq!(request.query_is_translated(), Ok(false));
    assert!(matches!(request.queries, Queries::BioseqSet(_)));
    assert_eq!(request.queries.num_queries(), 3);
    let queries = request.queries.bioseqs();
    assert_eq!(queries.len(), 3);
    let first = queries[0];
    assert_eq!(first.best_id().map(label).as_deref(), Some("lcl|Query_1"));
    assert_eq!(first.title(), Some(&b"nh1 random query without hits"[..]));
    assert_eq!(
        first.fasta_defline(),
        Some(&b">nh1 random query without hits"[..])
    );
    assert_eq!((first.molecule(), first.length()), (4, Some(300)));
    assert_eq!(
        first.seq_data().map(|data| data.choice().0),
        Some("ncbi2na")
    );
    assert_eq!(queries[1].length(), Some(600));
    assert_eq!(
        queries[1].fasta_defline(),
        Some(&b">mq1 query matching several subjects"[..])
    );
    let subjects = request.subject.bioseqs();
    assert_eq!(subjects.len(), 9);
    assert_eq!(
        subjects[0].best_id().map(label).as_deref(),
        Some("lcl|Subject_1")
    );
    assert_eq!(subjects[0].length(), Some(700));
    let params: Vec<_> = request
        .params()
        .iter()
        .map(|param| (param.set, String::from_utf8_lossy(&param.name).into_owned()))
        .collect();
    let algorithm = |name: &str| (ParamSet::Algorithm, name.to_string());
    assert_eq!(
        params,
        [
            algorithm("EvalueThreshold"),
            algorithm("MaskAtHash"),
            algorithm("UngappedMode"),
            algorithm("HitlistSize"),
            algorithm("EffectiveSearchSpace"),
            (ParamSet::Program, "LCaseMask".to_string()),
        ]
    );
    assert_eq!(hitlist(&request), Ok(500));

    // The database of a search, as stored, and the residue type `x_SetDatabase` gives it.
    for (name, program, protein, translated) in [
        ("rpsblast_mixnh.asn", "blastp", true, false),
        ("rpstblastn_mixnh.asn", "tblastn", true, true),
        ("blastx_strand_plus.asn", "blastx", true, true),
        ("blastn_multi.asn", "blastn", false, false),
    ] {
        let archives = read_fixture(name);
        let request = &search_of(&archives[0]).request;
        assert_eq!(request.program, program.as_bytes(), "{name}");
        assert_eq!(request.db_is_protein(), protein, "{name}");
        assert_eq!(request.query_is_translated(), Ok(translated), "{name}");
        let db: &[u8] = if name.starts_with("rps") {
            b"db/rps_mini"
        } else if name.starts_with("blastx") {
            b"db/bx"
        } else {
            b"db/bn"
        };
        assert!(
            matches!(&request.subject, Subject::Database(stored) if stored == db),
            "{name}"
        );
    }
}

/// psiblast's archives (`psiblast -num_iterations 3 -outfmt 11`: one object per iteration, the
/// first with the query's Bioseq-set, the next with the PSSM of the iteration before; `psiblast
/// -in_pssm`): the PSSM's query is the one query.
#[test]
fn the_psiblast_archives_are_read() {
    let iterations = read_fixture("psiblast_it3.asn");
    assert_eq!(iterations.len(), 2);
    for (archive, pssm) in iterations.iter().zip([false, true]) {
        let request = &search_of(archive).request;
        assert_eq!(
            (request.program.as_slice(), request.service.as_slice()),
            (&b"blastp"[..], &b"psi"[..])
        );
        assert_eq!(request.query_is_translated(), Ok(false));
        assert!(request.db_is_protein());
        assert!(matches!(&request.subject, Subject::Database(db) if db == b"pdb"));
        assert_eq!(request.queries.pssm().is_some(), pssm);
        assert_eq!(request.queries.num_queries(), 1);
        assert_eq!(
            query_ids(request),
            Ok(vec![Some("lcl|Query_1".to_string())])
        );
    }
    let archives = read_fixture("psiblast_in_pssm.asn");
    let search = search_of(&archives[0]);
    let queries = &search.request.queries;
    assert_eq!(
        queries
            .pssm()
            .map(|pssm| pssm.field("pssm").field("numColumns").int()),
        Some(30)
    );
    let bioseqs = queries.bioseqs();
    assert_eq!(bioseqs.len(), 1);
    assert_eq!(bioseqs[0].title(), Some(&b"q1 small query"[..]));
    assert_eq!(bioseqs[0].fasta_defline(), None);
    assert_eq!((bioseqs[0].molecule(), bioseqs[0].length()), (3, Some(30)));
    assert_eq!(
        bioseqs[0].seq_data().map(|data| data.choice().0),
        Some("ncbistdaa")
    );
    let masks = search.results.masks(queries.num_queries(), || {
        search.request.query_is_translated()
    });
    assert_eq!(masks.map(|masks| ranges(&masks)), Ok(vec![vec![]]));
    assert_eq!(
        search.results.karlin_blocks().gapped,
        (real(267, -3), real(44614779411758, -15), real(14, -2))
    );
    assert_eq!(
        search.results.search_stats(&search.request.service),
        SearchStats {
            effective_search_space: 4671,
            length_adjustment: 3,
            matrix: None,
        }
    );
}

/// A synthetic search request with the options `options` (written after the subject).
fn request(queries: &str, subject: &str, options: &str) -> Request {
    let text = format!(
        "Blast4-request ::= {{ ident \"x\", body queue-search {{ program \"blastn\", \
         service \"plain\", queries {queries}, subject {subject}{options} }} }}"
    );
    Request::from_value(value(named::BLAST4_REQUEST, &text))
        .expect("a search")
        .expect("not an error request")
}

/// A search request without `ident` stops NCBI's archive reader (`GetCreatedBy`); an error
/// archive's request needs none.
#[test]
fn a_search_request_has_an_ident() {
    let text = "Blast4-request ::= { body queue-search { program \"blastn\", service \"plain\", \
                queries seq-loc-list { whole local str \"q\" }, subject database \"db\" } }";
    assert_eq!(
        Request::from_value(value(named::BLAST4_REQUEST, text)).map(|request| request.is_some()),
        Err(Stop::Unassigned {
            type_name: "Blast4-request",
            member: "ident",
        })
    );
    let error = "Blast4-request ::= { body get-request-info { request-id \"Error\" } }";
    assert!(matches!(
        Request::from_value(value(named::BLAST4_REQUEST, error)),
        Ok(None)
    ));
}

/// The ids of the queries of a request (`x_ExtractQueryIds`), as FASTA labels.
fn query_ids(request: &Request) -> Result<Vec<Option<String>>, Stop> {
    Ok(request
        .queries
        .query_ids()?
        .into_iter()
        .map(|id| id.map(label))
        .collect())
}

/// The queries as each form gives them: a Bioseq-set's Bioseqs depth first (its number of
/// queries the entries at its top), a PSSM's query, the locations of a list; a PSSM without a
/// query stops NCBI (`GetQuery()` of an unset member). A subject list of locations is kept.
#[test]
fn the_queries_of_each_form_are_read() {
    let bioseq = |id: &str| {
        format!("seq {{ id {{ local str \"{id}\" }}, inst {{ repr raw, mol aa, length 1, seq-data ncbieaa \"A\" }} }}")
    };
    let nested = format!(
        "bioseq-set {{ seq-set {{ set {{ seq-set {{ {}, {} }} }}, {} }} }}",
        bioseq("a"),
        bioseq("b"),
        bioseq("c")
    );
    let set = request(&nested, "database \"db\"", "");
    assert_eq!(set.queries.num_queries(), 2);
    assert_eq!(
        query_ids(&set),
        Ok(["lcl|a", "lcl|b", "lcl|c"]
            .map(|id| Some(id.to_string()))
            .to_vec())
    );
    let pssm = request(
        &format!(
            "pssm {{ pssm {{ numRows 28, numColumns 1, query {} }} }}",
            bioseq("p")
        ),
        "database \"db\"",
        "",
    );
    assert_eq!(pssm.queries.num_queries(), 1);
    assert_eq!(query_ids(&pssm), Ok(vec![Some("lcl|p".to_string())]));
    assert!(pssm.queries.pssm().is_some());
    let no_query = request(
        "pssm { pssm { numRows 28, numColumns 1 } }",
        "database \"db\"",
        "",
    );
    let unassigned = Stop::Unassigned {
        type_name: "Pssm",
        member: "query",
    };
    assert_eq!(query_ids(&no_query), Err(unassigned.clone()));
    let results = Results::from_value(value(
        named::BLAST4_GET_SEARCH_RESULTS_REPLY,
        "Blast4-get-search-results-reply ::= { }",
    ));
    assert_eq!(results.alignment_sets(&no_query.queries), Err(unassigned));
    let locs = request(
        "seq-loc-list { int { from 0, to 9, id local str \"q\" }, null NULL }",
        "seq-loc-list { whole local str \"s\" }",
        "",
    );
    assert_eq!(locs.queries.num_queries(), 2);
    assert_eq!(query_ids(&locs), Ok(vec![Some("lcl|q".to_string()), None]));
    assert_eq!(locs.queries.seq_locs().count(), 2);
    assert!(matches!(&locs.subject, Subject::SeqLocList));
    assert!(locs.subject.bioseqs().is_empty());
    // A null query id is dereferenced only while a run of alignments is left.
    let align = "{ type partial, dim 2, segs denseg { numseg 1, ids { local str \"q\", \
                 local str \"s\" }, starts { 0, 0 }, lens { 5 } } }";
    let results = Results::from_value(value(
        named::BLAST4_GET_SEARCH_RESULTS_REPLY,
        &format!("Blast4-get-search-results-reply ::= {{ alignments {{ {align} }} }}"),
    ));
    assert_eq!(results.alignment_sets(&locs.queries), Ok(vec![0..1, 1..1]));
    let null_first = request(
        "seq-loc-list { null NULL, int { from 0, to 9, id local str \"q\" } }",
        "database \"db\"",
        "",
    );
    assert_eq!(
        results.alignment_sets(&null_first.queries),
        Err(Stop::NullPointer("GetSeqAlignSets: a null query id"))
    );
}

/// The hit list size of the options built for an archive's request.
fn hitlist(request: &Request) -> Result<i32, String> {
    archive_request(request)
        .map(|options| options.hitlist_size())
        .map_err(|e| String::from_utf8_lossy(e.get_msg()).into_owned())
}

/// `HitlistSize`: the last one of the algorithm and program options (the format options are
/// not read), else 500; `GetParamByName` the first of a list; another
/// type of value stops NCBI's option builder.
#[test]
fn the_hitlist_size_is_the_last_one_stored() {
    let locs = "seq-loc-list { int { from 0, to 9, id local str \"q\" } }";
    let db = "database \"db\"";
    assert_eq!(hitlist(&request(locs, db, "")), Ok(500));
    let options = ", algorithm-options { { name \"HitlistSize\", value integer 5 }, \
                   { name \"EvalueThreshold\", value cutoff e-value { 1, 10, 1 } }, \
                   { name \"HitlistSize\", value integer 6 } }, \
                   program-options { { name \"HitlistSize\", value integer 7 } }, \
                   format-options { { name \"HitlistSize\", value integer 9 } }";
    let stored = request(locs, db, options);
    assert_eq!(hitlist(&stored), Ok(7));
    // The format options reach the builder only where the service is psi, and only for the
    // program's adjustment: a PHI pattern there stops a PSI-BLAST request, and its hit list
    // size is still not read.
    let phi = ", program-options { { name \"HitlistSize\", value integer 7 } }, \
               format-options { { name \"HitlistSize\", value integer 9 }, \
               { name \"PHIPattern\", value string \"LLL\" } }";
    let blastp = |service: &str| {
        let text = format!(
            "Blast4-request ::= {{ ident \"x\", body queue-search {{ program \"blastp\", \
             service \"{service}\", queries {locs}, subject {db}{phi} }} }}"
        );
        Request::from_value(value(named::BLAST4_REQUEST, &text))
            .expect("a search")
            .expect("not an error request")
    };
    assert_eq!(hitlist(&blastp("plain")), Ok(7));
    assert_eq!(
        hitlist(&blastp("psi")),
        Err("Incorrect combination of option (PHIPattern) and program (blastp)".to_string())
    );
    let params: Vec<_> = stored
        .params()
        .iter()
        .map(|param| (param.set, param.value().choice().0))
        .collect();
    assert_eq!(
        params,
        [
            (ParamSet::Algorithm, "integer"),
            (ParamSet::Algorithm, "cutoff"),
            (ParamSet::Algorithm, "integer"),
            (ParamSet::Program, "integer"),
            (ParamSet::Format, "integer"),
        ]
    );
    let first = |set| {
        stored
            .param(set, b"HitlistSize")
            .map(|param| param.value().choice().1.int())
    };
    assert_eq!(
        (
            first(ParamSet::Algorithm),
            first(ParamSet::Program),
            first(ParamSet::Format)
        ),
        (Some(5), Some(7), Some(9))
    );
    let algorithm_only = ", algorithm-options { { name \"HitlistSize\", value integer 5 }, \
                          { name \"HitlistSize\", value integer 6 } }";
    assert_eq!(hitlist(&request(locs, db, algorithm_only)), Ok(6));
    let mistyped = ", algorithm-options { { name \"HitlistSize\", value string \"5\" } }";
    let error = archive_request(&request(locs, db, mistyped)).unwrap_err();
    assert_eq!(
        (error.site().file, error.site().line, error.kind()),
        (
            "src/algo/blast/api/blast_options_builder.cpp",
            541,
            "CRemoteBlastException::eServiceNotAvailable",
        )
    );
    assert_eq!(
        hitlist(&request(locs, db, mistyped)),
        Err("Internal: Error processing option [HitlistSize] type [12].".to_string())
    );
}

/// The masked ranges of each query as (id, from, to, frame).
fn ranges(masks: &[Vec<MaskRange<'_>>]) -> Vec<Vec<(String, i64, i64, i32)>> {
    masks
        .iter()
        .map(|ranges| {
            ranges
                .iter()
                .map(|range| (label(range.id), range.from, range.to, range.frame))
                .collect()
        })
        .collect()
}

/// The results of a synthetic reply whose members are `members`.
fn results(members: &str) -> Results {
    Results::from_value(value(
        named::BLAST4_GET_SEARCH_RESULTS_REPLY,
        &format!("Blast4-get-search-results-reply ::= {{ {members} }}"),
    ))
}

/// `GetMasks`: the masks of each query in order. A mask whose id is not the previous mask's
/// goes to the next query (so the masks of the first and third query of `blastn -lcase_masking`
/// go to the first and second, as NCBI gives them); the frames of a program that translates the
/// query are the masks' frames, else none.
#[test]
fn the_masks_are_taken_per_query() {
    let archives = read_fixture("blastn_lcase.asn");
    let search = search_of(&archives[0]);
    let masks = search
        .results
        .masks(search.request.queries.num_queries(), || {
            search.request.query_is_translated()
        })
        .expect("masks");
    let q = |id: &str, from: i64, to: i64, frame: i32| (format!("lcl|{id}"), from, to, frame);
    assert_eq!(
        ranges(&masks),
        vec![
            vec![q("Query_1", 10, 24, 0), q("Query_1", 588, 596, 0)],
            vec![
                q("Query_3", 10, 24, 0),
                q("Query_3", 92, 98, 0),
                q("Query_3", 288, 296, 0)
            ],
            vec![],
        ]
    );
    let archives = read_fixture("blastx_strand_plus.asn");
    let search = search_of(&archives[0]);
    let masks = search
        .results
        .masks(1, || search.request.query_is_translated())
        .expect("masks");
    let masks = ranges(&masks);
    assert_eq!(masks.len(), 1);
    let ends: Vec<_> = [0, 12, 13, 25, 26, 36]
        .iter()
        .map(|&at| masks[0][at].clone())
        .collect();
    assert_eq!(masks[0].len(), 37);
    assert_eq!(
        ends,
        [
            q("Query_1", 39, 69, 1),
            q("Query_1", 4647, 4692, 1),
            q("Query_1", 343, 388, 2),
            q("Query_1", 4903, 4945, 2),
            q("Query_1", 164, 197, 3),
            q("Query_1", 3734, 3806, 3),
        ]
    );
    let archives = read_fixture("rpsblast_mixnh.asn");
    let search = search_of(&archives[0]);
    let masks = search.results.masks(2, || panic!("asked without masks"));
    assert_eq!(masks.map(|masks| ranges(&masks)), Ok(vec![vec![], vec![]]));

    // Synthetic replies.
    let mask = |id: &str, frame: &str| {
        format!(
            "{{ locations {{ packed-int {{ {{ from 1, to 2, id local str \"{id}\" }} }} }}, \
             frame {frame} }}"
        )
    };
    let reply = |masks: &[String]| results(&format!("masks {{ {} }}", masks.join(", ")));
    let one = |id: &str, frame: i32| vec![(format!("lcl|{id}"), 1_i64, 2_i64, frame)];
    let abc = reply(&[mask("a", "plus1"), mask("b", "plus1"), mask("a", "plus1")]);
    assert_eq!(
        abc.masks(3, || Ok(false)).map(|masks| ranges(&masks)),
        Ok(vec![one("a", 0), one("b", 0), one("a", 0)])
    );
    assert_eq!(
        abc.masks(2, || Ok(false)).map(|masks| ranges(&masks)),
        Err(Stop::Undefined(
            "GetMasks: a mask of a query after the last (retval[++query_index])"
        ))
    );
    let minus = reply(&[
        mask("a", "minus1"),
        mask("a", "minus2"),
        mask("a", "minus3"),
    ]);
    assert_eq!(
        minus.masks(1, || Ok(true)).map(|masks| ranges(&masks)),
        Ok(vec![[one("a", -1), one("a", -2), one("a", -3)].concat()])
    );
    let notset = reply(&[mask("a", "notset")]);
    assert_eq!(
        notset.masks(1, || Ok(false)).map(|masks| ranges(&masks)),
        Ok(vec![one("a", 0)])
    );
    assert_eq!(
        notset.masks(1, || Ok(true)).map(|masks| ranges(&masks)),
        Err(Stop::Abort(
            "NetworkFrame2FrameNumber: a frame other than plus1-minus3"
        ))
    );
    assert_eq!(
        notset
            .masks(1, || Err(Stop::Abort("x")))
            .map(|masks| ranges(&masks)),
        Err(Stop::Abort("x"))
    );
    let int = results(
        "masks { { locations { int { from 1, to 2, id local str \"a\" } }, frame plus1 } }",
    );
    assert_eq!(
        int.masks(1, || Ok(false)).map(|masks| ranges(&masks)),
        Err(Stop::InvalidSelection {
            module: "NCBI-Seqloc",
            type_name: "Seq-loc",
            required: "packed-int",
            current: "int",
        })
    );
    let empty = results("masks { { locations { }, frame plus1 } }");
    assert_eq!(
        empty.masks(1, || Ok(false)).map(|masks| ranges(&masks)),
        Err(Stop::Undefined(
            "GetMasks: front() of a mask without locations"
        ))
    );
    let null = results("masks { { locations { null NULL }, frame plus1 } }");
    assert_eq!(
        null.masks(1, || Ok(false)).map(|masks| ranges(&masks)),
        Err(Stop::NullPointer(
            "GetMasks: the id of a mask's first location"
        ))
    );
    let none = results("masks { }");
    assert_eq!(
        none.masks(1, || panic!("asked without masks"))
            .map(|masks| ranges(&masks)),
        Ok(vec![vec![]])
    );
}

/// A REAL of ASN.1 text `{ mantissa, 10, exp }`.
fn real(mantissa: i64, exp: i32) -> f64 {
    real_of_parts(mantissa as f64, false, 10, exp).expect("base 10")
}

/// `GetResultSet`'s Karlin-Altschul parameters and search statistics.
#[test]
fn the_karlin_blocks_and_statistics_are_read() {
    let archives = read_fixture("rpsblast_mixnh.asn");
    let search = search_of(&archives[0]);
    assert_eq!(
        search.results.karlin_blocks(),
        KarlinBlocks {
            ungapped: (
                real(33214346218523, -14),
                real(14237250382963, -14),
                real(741982920340517, -15)
            ),
            gapped: (real(267, -3), real(4996776, -8), real(14, -2)),
        }
    );
    assert_eq!(
        search.results.search_stats(&search.request.service),
        SearchStats {
            effective_search_space: 14248,
            length_adjustment: 4,
            matrix: None,
        }
    );
    let archives = read_fixture("blastx_strand_plus.asn");
    let search = search_of(&archives[0]);
    assert_eq!(
        search.results.karlin_blocks(),
        KarlinBlocks {
            ungapped: (
                real(317605957635731, -15),
                real(133956144488482, -15),
                real(401214524497119, -15)
            ),
            gapped: (real(267, -3), real(41, -3), real(14, -2)),
        }
    );
    let stats = |name: &str| {
        let archives = read_fixture(name);
        let search = search_of(&archives[0]);
        let stats = search.results.search_stats(&search.request.service);
        (stats.effective_search_space, stats.length_adjustment)
    };
    assert_eq!(stats("blastx_strand_plus.asn"), (22502727, 67));
    assert_eq!(stats("blastn_lcase.asn"), (3102284, 14));
    assert_eq!(stats("rpstblastn_mixnh.asn"), (38532, 22));

    // The blocks in order until both kinds are found.
    let block = |n: i64, gapped: &str| {
        format!(
            "{{ lambda {{ {n}, 10, 0 }}, k {{ {}, 10, 0 }}, h {{ {}, 10, 0 }}, gapped {gapped} }}",
            n + 1,
            n + 2
        )
    };
    let blocks = results(&format!(
        "ka-blocks {{ {}, {}, {}, {} }}",
        block(1, "FALSE"),
        block(4, "FALSE"),
        block(7, "TRUE"),
        block(10, "FALSE")
    ));
    assert_eq!(
        blocks.karlin_blocks(),
        KarlinBlocks {
            ungapped: (4.0, 5.0, 6.0),
            gapped: (7.0, 8.0, 9.0),
        }
    );
    assert_eq!(results("").karlin_blocks(), KarlinBlocks::default());
    // Each text in turn, changed in place: the last one with a text gives its number (0 when the
    // rest is no number); the matrix of the rps services from the first one with "Matrix: ".
    let stats = results(
        "search-stats { \"Effective search space used: 7\", \
         \"Length adjustment: Effective search space used: 42\", \
         \"Matrix: Matrix: PAM30\", \"Matrix: BLOSUM62\" }",
    );
    let read = |service: &[u8]| stats.search_stats(service);
    let expected = |matrix: Option<&[u8]>| SearchStats {
        effective_search_space: 0,
        length_adjustment: 42,
        matrix: matrix.map(<[u8]>::to_vec),
    };
    assert_eq!(read(b"rpsblast"), expected(Some(&b"PAM30"[..])));
    assert_eq!(read(b"rpstblastn"), expected(Some(&b"PAM30"[..])));
    assert_eq!(read(b"plain"), expected(None));
    assert_eq!(
        results("").search_stats(b"rpsblast"),
        SearchStats::default()
    );
}

/// `LoadFromArchive`'s messages: each stored message with its code as the severity, "Unknown
/// Error Code: " before one whose code is no severity, none for an error without a message; an
/// error archive has no search, and a request of another kind stops `CImportStrategy`.
#[test]
fn the_stored_messages_are_read_by_their_code() {
    let archives = read_fixture("warnings_batch.asn");
    let messages: Vec<_> = archives.iter().map(|a| a.messages.clone()).collect();
    assert_eq!(
        messages,
        vec![
            vec![Message {
                severity: Severity::Warning,
                text: b"Warning: [blastn] Examining 5 or more matches is recommended".to_vec(),
            }],
            vec![],
            vec![],
        ]
    );
    let archives = read_fixture("error_nodb.asn");
    assert_eq!(archives.len(), 1);
    assert!(matches!(archives[0].search, Ok(None)));
    assert_eq!(
        archives[0].messages,
        vec![Message {
            severity: Severity::Error,
            text: concat!(
                "BLAST Database error: No alias or index file found for nucleotide database ",
                "[db/no_such_db] in search path [/dev/shm/rgff3/in::]. Please verify the ",
                "spelling of the BLAST database and its molecule type."
            )
            .as_bytes()
            .to_vec(),
        }]
    );

    let archive = |body: &str, messages: &str| {
        Archive::from_value(value(
            named::BLAST4_ARCHIVE,
            &format!("Blast4-archive ::= {{ request {{ body {body} }}, results {{ }}{messages} }}"),
        ))
    };
    let error = "get-request-info { request-id \"Error\" }";
    let coded = archive(
        error,
        ", messages { { code 0, message \"i\" }, { code 1, message \"w\" }, \
         { code 2, message \"e\" }, { code 3, message \"c\" }, { code 4, message \"f\" }, \
         { code 5, message \"t\" }, { code 6, message \"u\" }, { code -1, message \"n\" }, \
         { code 2 } }",
    );
    let message = |severity, text: &str| Message {
        severity,
        text: text.as_bytes().to_vec(),
    };
    assert_eq!(
        coded.messages,
        vec![
            message(Severity::Info, "i"),
            message(Severity::Warning, "w"),
            message(Severity::Error, "e"),
            message(Severity::Critical, "c"),
            message(Severity::Fatal, "f"),
            message(Severity::Trace, "t"),
            message(Severity::Error, "Unknown Error Code: u"),
            message(Severity::Error, "Unknown Error Code: n"),
        ]
    );
    assert!(matches!(coded.search, Ok(None)));
    let no_body = (
        "src/algo/blast/api/search_strategy.cpp",
        72,
        "CBlastException::eInvalidArgument",
        "No body in CBlast4_request".to_string(),
    );
    for body in ["get-request-info { request-id \"x\" }", "get-matrices NULL"] {
        let other = archive(body, "");
        assert!(other.messages.is_empty());
        assert_eq!(exception(other.search.unwrap_err()), no_body, "{body}");
    }
}

/// `CSeq_loc::GetId()` of a location as a FASTA label.
fn loc_id(text: &str) -> Option<String> {
    let loc = value(named::SEQ_LOC, &format!("Seq-loc ::= {text}"));
    seq_loc_id(V::new(named::SEQ_LOC, &loc)).map(label)
}

/// `CSeq_loc::GetId()`: the one id of every part (`null` has none and agrees with any), none
/// when two differ, none for a `feat`.
#[test]
fn a_location_has_one_id_or_none() {
    let a = Some("lcl|a".to_string());
    let cases = [
        ("null NULL", None),
        ("whole local str \"a\"", a.clone()),
        ("empty local str \"a\"", a.clone()),
        ("int { from 0, to 1, id local str \"a\" }", a.clone()),
        ("pnt { point 3, id local str \"a\" }", a.clone()),
        ("packed-pnt { id local str \"a\", points { 1, 2 } }", a.clone()),
        (
            "packed-int { { from 0, to 1, id local str \"a\" }, { from 5, to 6, id local str \"a\" } }",
            a.clone(),
        ),
        (
            "packed-int { { from 0, to 1, id local str \"a\" }, { from 5, to 6, id local str \"b\" } }",
            None,
        ),
        (
            "mix { null NULL, int { from 0, to 1, id local str \"a\" }, mix { whole local str \"a\" } }",
            a.clone(),
        ),
        ("mix { whole local str \"a\", mix { whole local str \"b\" } }", None),
        ("equiv { whole local str \"a\", whole local str \"a\" }", a.clone()),
        ("bond { a { point 1, id local str \"a\" } }", a.clone()),
        (
            "bond { a { point 1, id local str \"a\" }, b { point 2, id local str \"b\" } }",
            None,
        ),
        ("feat local id 1", None),
        ("mix { whole gi 5, whole gi 5 }", Some("gi|5".to_string())),
    ];
    for (text, expected) in cases {
        assert_eq!(loc_id(text), expected, "{text}");
    }
}

/// `CSeq_align::GetSeq_id` of a row of an alignment as a FASTA label, or the place, class and
/// message of NCBI's exception.
fn align_id(segs: &str, row: usize) -> Result<String, (&'static str, u32, &'static str, String)> {
    let align = value(
        named::SEQ_ALIGN,
        &format!("Seq-align ::= {{ type partial, dim 2, segs {segs} }}"),
    );
    seq_align_id(V::new(named::SEQ_ALIGN, &align), row)
        .map(label)
        .map_err(exception)
}

/// `CSeq_align::GetSeq_id` of a row of each kind of segments, and NCBI's exceptions.
#[test]
fn an_alignment_row_has_the_id_of_its_segments() {
    let invalid = |file: &'static str, line: u32, class: &str| {
        Err((
            file,
            line,
            "CSeqalignException::eInvalidRowNumber",
            format!("{class}::GetSeq_id(): can not get seq-id for the row requested."),
        ))
    };
    let align_cpp = "src/objects/seq/../seqalign/Seq_align.cpp";
    let ok = |id: &str| Ok(format!("lcl|{id}"));
    let denseg = "denseg { numseg 1, ids { local str \"q\", local str \"s\" }, \
                  starts { 0, 0 }, lens { 5 } }";
    assert_eq!(align_id(denseg, 0), ok("q"));
    assert_eq!(align_id(denseg, 1), ok("s"));
    assert_eq!(
        align_id(denseg, 2),
        invalid(
            "src/objects/seq/../seqalign/Dense_seg.cpp",
            159,
            "CDense_seg"
        )
    );
    let dendiag = "dendiag { { ids { local str \"q\" }, starts { 0 }, len 5 }, \
                   { ids { local str \"q\", local str \"s\" }, starts { 0, 0 }, len 5 } }";
    assert_eq!(align_id(dendiag, 1), ok("s"));
    assert_eq!(align_id(dendiag, 2), invalid(align_cpp, 396, "CSeq_align"));
    let std_seg = "std { { loc { null NULL, int { from 0, to 4, id local str \"s\" } } }, \
               { loc { int { from 0, to 4, id local str \"q\" }, null NULL } } }";
    assert_eq!(align_id(std_seg, 0), ok("q"));
    assert_eq!(align_id(std_seg, 1), ok("s"));
    let std_ids = "std { { ids { local str \"x\", local str \"y\" }, \
                   loc { whole local str \"q\", whole local str \"s\" } } }";
    assert_eq!(align_id(std_ids, 1), ok("y"));
    let disc = "disc { { type partial, segs denseg { dim 1, numseg 1, ids { local str \"a\" }, \
                starts { 0 }, lens { 5 } } }, { type partial, segs denseg { numseg 1, \
                ids { local str \"b\", local str \"c\" }, starts { 0, 0 }, lens { 5 } } } }";
    assert_eq!(align_id(disc, 0), ok("a"));
    assert_eq!(align_id(disc, 1), ok("c"));
    assert_eq!(align_id(disc, 2), invalid(align_cpp, 396, "CSeq_align"));
    let sparse = "sparse { rows { { first-id local str \"m\", second-id local str \"r\", \
                  numseg 1, first-starts { 0 }, second-starts { 0 }, lens { 5 } } } }";
    assert_eq!(align_id(sparse, 0), ok("m"));
    assert_eq!(align_id(sparse, 1), ok("r"));
    assert_eq!(
        align_id(sparse, 2),
        invalid(
            "src/objects/seq/../seqalign/Sparse_seg.cpp",
            87,
            "CSparse_seg"
        )
    );
    let spliced = "spliced { product-id local str \"p\", genomic-id local str \"g\", \
                   product-type transcript, exons { } }";
    assert_eq!(align_id(spliced, 0), ok("p"));
    assert_eq!(align_id(spliced, 1), ok("g"));
    assert_eq!(align_id(spliced, 2), invalid(align_cpp, 396, "CSeq_align"));
    let packed = "packed { numseg 1, ids { local str \"q\", local str \"s\" }, \
                  starts { 0, 0 }, present '03'H, lens { 5 } }";
    assert_eq!(
        align_id(packed, 0),
        Err((
            align_cpp,
            391,
            "CSeqalignException::eUnsupported",
            "CSeq_align::GetSeq_id() currently does not handle this type of alignment.".to_string()
        ))
    );
}

/// A strategy's request needs no `ident` (`CImportStrategy` never asks for it); an error
/// archive's request is no search there (`No body in CBlast4_request`).
#[test]
fn a_strategy_request_needs_no_ident() {
    let text = "Blast4-request ::= { body queue-search { program \"blastn\", service \"plain\", \
                queries seq-loc-list { whole local str \"q\" }, subject database \"db\" } }";
    let request = Request::from_strategy(value(named::BLAST4_REQUEST, text)).unwrap();
    assert_eq!(&request.program[..], &b"blastn"[..]);
    let error = "Blast4-request ::= { body get-request-info { request-id \"Error\" } }";
    match Request::from_strategy(value(named::BLAST4_REQUEST, error)) {
        Err(e) => assert_eq!(e.get_msg(), b"No body in CBlast4_request"),
        other => panic!("{other:?}"),
    }
}
