//! The Blast4-archive reader on NCBI rpsblast 2.17.0 archives (the fixtures of
//! `LOSAT/tests/fixtures/rpsbproc/arch/`) and on the edited archives that NCBI's rpsbproc read
//! as the end of its input (PS16a probe 3: `arch/bad/`).

use super::*;
use crate::report::serial::asn_text_reader::AsnTextErrorKind;

const ONE: &[u8] = include_bytes!("../../tests/fixtures/rpsbproc/arch/one.asn");

/// The objects of `text`, read as rpsbproc reads them: the outcome of each `next_archive`
/// until the end (`Ok(None)`) or a refusal.
fn read_all(text: &[u8]) -> Vec<std::result::Result<Archive, AsnTextError>> {
    let start = text
        .iter()
        .position(|&b| (b as i8) > b' ' as i8)
        .expect("a byte that is not white space");
    let mut stream = ArchiveStream::new(&text[start + 1..], &text[start..=start]);
    let mut out = Vec::new();
    loop {
        match stream.next_archive() {
            Ok(Some(archive)) => out.push(Ok(archive)),
            Ok(None) => return out,
            Err(error) => {
                out.push(Err(error));
                assert!(matches!(stream.next_archive(), Ok(None)));
                return out;
            }
        }
    }
}

#[test]
fn ncbi_archive_is_read() {
    let archives = read_all(ONE);
    assert_eq!(archives.len(), 1);
    let archive = archives[0].as_ref().unwrap();
    assert_eq!(archive.ident.as_deref(), Some(&b"2.17.0+"[..]));
    let RequestBody::QueueSearch(search) = &archive.body else {
        panic!("queue-search");
    };
    assert_eq!(search.program, b"blastp");
    assert!(matches!(&search.subject, Subject::Database(db) if db == b"Cdd"));
    assert!(search.has_algorithm_options);
    assert!(matches!(&search.matrix_name, Some(ParamValue::String(name)) if name == b"BLOSUM62"));
    assert!(matches!(
        search.evalue_threshold,
        Some(ParamValue::Cutoff(Cutoff::EValue(value))) if value == 0.01
    ));
    let Queries::BioseqSet(entries) = &search.queries else {
        panic!("bioseq-set");
    };
    assert_eq!(entries.len(), 1);
    let SeqEntry::Seq(query) = &entries[0] else {
        panic!("seq");
    };
    assert_eq!(
        query.ids,
        vec![ParsedSeqId::Local(ObjectId::Str(b"Query_1".to_vec()))]
    );
    assert_eq!(query.mol, 3);
    assert_eq!(query.length, Some(79));
    assert!(query
        .title
        .as_ref()
        .unwrap()
        .starts_with(b"NP_414783.1 protein YkfF"));
    assert!(query.unread_descriptor.is_none());
    let aligns = archive.results.alignments.as_ref().unwrap();
    assert_eq!(aligns.len(), 1);
    let Segs::Denseg(seg) = &aligns[0].segs else {
        panic!("denseg");
    };
    assert_eq!(seg.starts, vec![4, 0]);
    assert_eq!(seg.lens, vec![70]);
    assert!(matches!(&seg.ids[1], ParsedSeqId::General(tag) if tag.db == b"CDD"));
    assert_eq!(aligns[0].score, Some(337));
    assert_eq!(aligns[0].num_ident, Some(60));
    assert_eq!(aligns[0].e_value, Some(504742318164089.0 * 1e-57));
    assert!(aligns[0].score_of_other_type.is_none());
    assert!(archive.results.masks.as_ref().unwrap().is_empty());
}

// PS16a probe 3: NCBI's rpsbproc reads two objects without a line between them, and stops
// without an error at an object that its typed reader refuses.
#[test]
fn objects_follow_each_other_and_a_bad_object_ends_the_input() {
    let twice = include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/twice_nonl.asn");
    assert_eq!(read_all(twice).len(), 2);
    for bad in [
        &include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/dupmember.asn")[..],
        include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/unkmember.asn"),
        include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/trailcomma.asn"),
        include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/trunc.asn"),
        include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/noresults.asn"),
    ] {
        assert!(read_all(bad).is_empty());
    }
    let second_bad = include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/second_bad.asn");
    assert_eq!(read_all(second_bad).len(), 1);
    let trail = include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/trailgarb.asn");
    assert_eq!(read_all(trail).len(), 1);
    // Comments and CRLF line ends are white space.
    for good in [
        &include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/comment.asn")[..],
        include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/crlf.asn"),
    ] {
        assert_eq!(read_all(good).len(), 1);
    }
}

// A form that LOSAT does not read is refused once the object reads whole; when NCBI's reader
// fails later in the same object, the object ends the input as NCBI's does.
#[test]
fn refusals_wait_for_the_whole_object() {
    let text = String::from_utf8(ONE.to_vec()).unwrap();
    let long_real = text.replace(
        "{ 504742318164089, 10, -57 }",
        "{ 5047423181640891234, 10, -61 }",
    );
    assert_ne!(long_real, text);
    let read = read_all(long_real.as_bytes());
    assert_eq!(read.len(), 1);
    let error = read[0].as_ref().unwrap_err();
    assert_eq!(error.kind, AsnTextErrorKind::Unsupported);
    assert!(error.message.contains("18 significant digits"), "{error}");
    let then_bad = long_real.replace("gapped FALSE", "gapped MAYBE");
    assert_ne!(then_bad, long_real);
    assert!(read_all(then_bad.as_bytes()).is_empty());
    // A byte outside the visible range in a VisibleString, here the title (NCBI replaces it by
    // '#' and posts an error).
    let invisible = text.replacen(
        "title \"NP_414783.1 protein YkfF",
        "title \"NP_414783.1 protein Yk\u{e9}F",
        1,
    );
    assert_ne!(invisible, text);
    let read = read_all(invisible.as_bytes());
    assert!(read[0].as_ref().unwrap_err().message.contains("byte 0xc3"));
}

// A bad byte in a string that the shared readers read (here a Dbtag's `db` in an alignment's
// Seq-id, line 131: `db "C\xE9D"`) is refused as one in a title is: NCBI stores `#` and posts
// `Error: (803.7) [...] Blast4-archive...general.db` / `Bad char [0xE9] in string at line 131` /
// the chunk (review M4a-1b F1; the file is `one.asn` with that byte and without
// `search-stats`). The refusal names the line and the chunk.
#[test]
fn a_bad_byte_in_any_string_is_refused_where_it_is() {
    let text = include_bytes!("../../tests/fixtures/rpsbproc/arch/bad/bad_db_nostats.asn");
    let read = read_all(text);
    assert_eq!(read.len(), 1);
    let error = read[0].as_ref().unwrap_err();
    assert_eq!(error.kind, AsnTextErrorKind::Unsupported);
    assert_eq!(error.line, 131);
    assert!(
        error.message.contains("byte 0xe9 (\"C\\xe9D\""),
        "{}",
        error.message
    );
    // The same with `search-stats`; with another refusal of the object (a long REAL, line 102)
    // the first one in the text is named.
    let with_stats =
        String::from_utf8(ONE.to_vec())
            .unwrap()
            .replacen("db \"CDD\",", "db \"C\u{1}D\",", 1);
    let read = read_all(with_stats.as_bytes());
    assert!(read[0].as_ref().unwrap_err().message.contains("byte 0x1"));
    let then_real = with_stats.replace(
        "{ 504742318164089, 10, -57 }",
        "{ 5047423181640891234, 10, -61 }",
    );
    let read = read_all(then_real.as_bytes());
    let error = read[0].as_ref().unwrap_err();
    assert_eq!(error.line, 102, "{error}");
    assert!(error.message.contains("18 significant digits"), "{error}");
    // A bad byte in an object that NCBI's reader then fails to read: NCBI has posted the error,
    // which LOSAT's rpsbproc does not, so it is refused too.
    let then_bad = with_stats.replace("gapped FALSE", "gapped MAYBE");
    assert_ne!(then_bad, with_stats);
    let read = read_all(then_bad.as_bytes());
    assert_eq!(read.len(), 1);
    assert!(read[0].as_ref().unwrap_err().message.contains("byte 0x1"));
}

// The `CFastaReader` DefLine is a UTF8String, which NCBI reads with any byte (PS16a review
// S2-4, oracles O3-4 and O3-5): a TAB or UTF-8 there is read, and so are the other values of a
// User-field by their types.
#[test]
fn a_user_field_string_takes_any_byte() {
    let text = String::from_utf8(ONE.to_vec()).unwrap();
    let defline = text.replacen(
        "data str \">NP_414783.1 protein YkfF",
        "data str \">NP_414783.1 protein\t\u{3b1}-YkfF",
        1,
    );
    assert_ne!(defline, text);
    let read = read_all(defline.as_bytes());
    assert_eq!(read.len(), 1);
    assert!(read[0].is_ok());
    let fields = text.replacen(
        "label str \"DefLine\",",
        "label str \"Other\",\n data strs { \"a\tb\", \"\u{e9}\" } },\n {\n label id 5,\n num 2,\n data fields { { label str \"x\", data object { type id 1, data { { label id 2, data reals { { 5, 10, 0 } } } } } } } },\n {\n label str \"DefLine\",",
        1,
    );
    assert_ne!(fields, text);
    let read = read_all(fields.as_bytes());
    assert_eq!(read.len(), 1);
    assert!(read[0].is_ok(), "{:?}", read[0]);
    // A User-field value that does not read as its type is NCBI's error.
    let bad_int = fields.replacen("label id 5,", "label id 5, data int x }, { label id 6,", 1);
    assert!(read_all(bad_int.as_bytes()).is_empty());
}

/// A reader that gives one byte per `read`.
struct OneByte<'a>(&'a [u8]);

impl std::io::Read for OneByte<'_> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let Some((&first, rest)) = self.0.split_first() else {
            return Ok(0);
        };
        buf[0] = first;
        self.0 = rest;
        Ok(1)
    }
}

// The stream source refills at every byte, and a text of many blocks drops the bytes it has read
// (`StreamSource::release_before`); the objects read as from one slice.
#[test]
fn a_stream_refills_and_releases() {
    let three = include_bytes!("../../tests/fixtures/rpsbproc/arch/three.asn");
    let mut text = Vec::new();
    while text.len() <= 8 * 64 * 1024 {
        text.extend_from_slice(three);
    }
    let copies = text.len() / three.len();
    let whole = read_all(&text);
    assert_eq!(whole.len(), copies);
    let mut stream = ArchiveStream::new(OneByte(&text[1..]), &text[..1]);
    let mut largest = 0;
    for expected in &whole {
        let archive = stream.next_archive().unwrap().unwrap();
        assert_eq!(
            format!("{archive:?}"),
            format!("{:?}", expected.as_ref().unwrap())
        );
        largest = largest.max(stream.buffered());
    }
    assert!(matches!(stream.next_archive(), Ok(None)));
    // The bytes kept: at most the four blocks read before a release and what follows them.
    assert!(largest < 5 * 64 * 1024, "{largest} of {}", text.len());
}

// The type name of the root object: another type is NCBI's error, no name is read as the
// archive (`SkipFileHeader`, objistr.cpp:892-898).
#[test]
fn the_root_type_name_is_checked() {
    let text = String::from_utf8(ONE.to_vec()).unwrap();
    let other = text.replacen("Blast4-archive ::=", "Blast4-request ::=", 1);
    assert!(read_all(other.as_bytes()).is_empty());
    let unnamed = text.replacen("Blast4-archive ::=", "::=", 1);
    assert_eq!(read_all(unnamed.as_bytes()).len(), 1);
}

// `Seq-inst.ext`: the literals of a `delta` with their `seq-data` (review oracle R5:
// `mfirst_delta.asn`, two NCBIeaa literals); a `loc` item, a literal without data or with `gap`
// data and the other `Seq-ext` variants give none; a literal without its `length` ends the
// input as NCBI's typed reader does.
#[test]
fn delta_literals_are_read() {
    let delta: &[u8] = include_bytes!("../../tests/fixtures/rpsbproc/arch/mfirst_delta.asn");
    let literals = |text: &[u8]| {
        let archives = read_all(text);
        let archive = archives[0].as_ref().unwrap();
        let RequestBody::QueueSearch(search) = &archive.body else {
            panic!("queue-search");
        };
        let Queries::BioseqSet(entries) = &search.queries else {
            panic!("bioseq-set");
        };
        let SeqEntry::Seq(query) = &entries[0] else {
            panic!("seq");
        };
        assert!(query.seq_data.is_none());
        query
            .delta_literals
            .iter()
            .map(|(name, data)| (*name, data.len()))
            .collect::<Vec<_>>()
    };
    assert_eq!(literals(delta), [("ncbieaa", 100), ("ncbieaa", 85)]);
    let text = String::from_utf8(delta.to_vec()).unwrap();
    let first = "                literal {\n                  length 100,";
    assert!(text.contains(first));
    let edited = text.replacen(
        first,
        "                loc whole local str \"Query_1\",\n                literal {\n                  length 3\n                },\n                literal {\n                  length 5,\n                  seq-data gap {\n                    type unknown\n                  }\n                },\n                literal {\n                  length 100,\n                  fuzz lim unk,",
        1,
    );
    assert_eq!(
        literals(edited.as_bytes()),
        [("ncbieaa", 100), ("ncbieaa", 85)]
    );
    let start = text.find("ext delta {").unwrap();
    let end = text[start..].find("\n              }\n").unwrap() + start;
    let seg = format!(
        "{}ext seg {{\n                whole local str \"Query_1\"\n{}",
        &text[..start],
        &text[end + 1..]
    );
    assert!(literals(seg.as_bytes()).is_empty());
    let no_length = text.replacen("                  length 100,\n", "", 1);
    assert!(read_all(no_length.as_bytes()).is_empty());
}

// The `Std-seg`s of an rpstblastn archive (seqalign.asn:236-241, seqloc.asn:79-107): ids,
// `int` and `pnt` locations with their strands, `empty`, `null` (`NULL`) and `whole`; a mask
// with its frame; a nucleotide query's `ncbi2na` data.
#[test]
fn std_segs_of_a_translated_search() {
    let text = br#"Blast4-archive ::= {
  request {
    ident "2.17.0+",
    body queue-search {
      program "tblastn",
      service "rpsblast",
      queries bioseq-set {
        seq-set {
          seq {
            id { local str "Query_1" },
            descr { title "q" },
            inst { repr raw, mol na, length 8, seq-data ncbi2na '1BE4'H }
          }
        }
      },
      subject database "Cdd"
    }
  },
  results {
    alignments {
      {
        type partial,
        dim 2,
        score { { id str "e_value", value real { 1, 10, -5 } } },
        segs std {
          {
            dim 2,
            ids { local str "Query_1", general { db "CDD", tag id 5 } },
            loc {
              int { from 3, to 8, strand minus, id local str "Query_1" },
              int { from 0, to 1, strand unknown, id general { db "CDD", tag id 5 } }
            }
          },
          { loc { pnt { point 2, strand plus, id local str "Query_1" }, empty general { db "CDD", tag id 5 } } },
          { loc { null NULL, whole general { db "CDD", tag id 5 } } }
        }
      }
    },
    masks { { locations { packed-int { { from 0, to 2, id local str "Query_1" } } }, frame plus1 } }
  }
}
"#;
    let archives = read_all(text);
    assert_eq!(archives.len(), 1);
    let archive = archives[0].as_ref().unwrap();
    let RequestBody::QueueSearch(search) = &archive.body else {
        panic!("queue-search");
    };
    assert_eq!(search.program, b"tblastn");
    let Queries::BioseqSet(entries) = &search.queries else {
        panic!("bioseq-set");
    };
    let SeqEntry::Seq(query) = &entries[0] else {
        panic!("seq");
    };
    assert_eq!(query.mol, 4);
    assert_eq!(query.seq_data, Some(("ncbi2na", vec![0x1b, 0xe4])));
    let aligns = archive.results.alignments.as_ref().unwrap();
    let Segs::Std(std) = &aligns[0].segs else {
        panic!("std");
    };
    assert!(matches!(
        &std.query_id,
        QueryRowId::Id(ParsedSeqId::Local(ObjectId::Str(name))) if name == b"Query_1"
    ));
    assert!(matches!(&std.hit_id, Some(ParsedSeqId::General(tag)) if tag.db == b"CDD"));
    assert_eq!(
        std.segs,
        vec![
            StdPair {
                query: StdLoc::Int {
                    from: 3,
                    to: 8,
                    strand: Some(2)
                },
                hit: StdLoc::Int {
                    from: 0,
                    to: 1,
                    strand: Some(0)
                },
            },
            StdPair {
                query: StdLoc::Pnt {
                    point: 2,
                    strand: Some(1)
                },
                hit: StdLoc::Empty,
            },
            StdPair {
                query: StdLoc::Null,
                hit: StdLoc::Whole,
            },
        ]
    );
    let masks = archive.results.masks.as_ref().unwrap();
    assert!(matches!(masks[0].first, Some(MaskLocation::PackedInt(_))));
}

// `CSeq_align::GetSeq_id(0)` (Seq_align.cpp:339-358, quoted at `std_segs`) as read: segments
// without ids give the id of their first location, `null` going on to the next segment; a first
// location of another variant is named; no id at all is `Missing`.
#[test]
fn query_row_id_of_std_segs_without_ids() {
    let archive = |locs: &str| {
        format!(
            r#"Blast4-archive ::= {{
  request {{
    body queue-search {{
      program "tblastn",
      service "rpsblast",
      queries bioseq-set {{ seq-set {{ }} }},
      subject database "Cdd"
    }}
  }},
  results {{
    alignments {{
      {{
        type partial,
        segs std {{
          {locs}
        }}
      }}
    }}
  }}
}}
"#
        )
    };
    let read = |locs: &str| {
        let text = archive(locs);
        let mut archives = read_all(text.as_bytes());
        let archive = archives.remove(0).unwrap();
        let Segs::Std(std) = archive.results.alignments.unwrap().remove(0).segs else {
            panic!("std");
        };
        std
    };
    let std = read(
        r#"{ loc { null NULL, empty local str "x" } },
          { loc { int { from 1, to 3, id local str "Q2" }, empty local str "x" } },
          { ids { local str "Q3", local str "H" }, loc { null NULL, null NULL } }"#,
    );
    assert!(matches!(
        &std.query_id,
        QueryRowId::Id(ParsedSeqId::Local(ObjectId::Str(name))) if name == b"Q2"
    ));
    assert!(std.hit_id.is_none());
    let std = read(
        r#"{ loc { packed-int { { from 1, to 3, id local str "Q" } }, empty local str "x" } }"#,
    );
    assert!(matches!(std.query_id, QueryRowId::Location("packed-int")));
    assert_eq!(std.segs[0].query, StdLoc::Other);
    let std = read(r#"{ loc { null NULL } }"#);
    assert!(matches!(std.query_id, QueryRowId::Missing));
    assert_eq!(std.segs[0].hit, StdLoc::Missing);
}
