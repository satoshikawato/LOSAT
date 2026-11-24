//! The Blast4-archive reader on NCBI rpsblast 2.17.0 archives (the fixtures of
//! `LOSAT/tests/fixtures/rpsbproc/arch/`) and on the edited archives that NCBI's rpsbproc read
//! as the end of its input (PS16a probe 3: `arch/bad/`).

use super::*;
use crate::report::serial::asn_types::read::Post;

const ONE: &[u8] = include_bytes!("../../tests/fixtures/rpsbproc/arch/one.asn");

type Objects = Vec<std::result::Result<Archive, AsnTextError>>;

/// The objects of `text`, read as rpsbproc reads them: the outcome of each `next_archive`
/// until the end (`Ok(None)`) or a refusal.
fn read_all(text: &[u8]) -> Objects {
    read_with_posts(text).0
}

/// `read_all` with the lines that NCBI's reader posts while it reads the objects (`Post::line`,
/// each given to the stream's posts when it is posted), in order.
fn read_with_posts(text: &[u8]) -> (Objects, Vec<String>) {
    let start = text
        .iter()
        .position(|&b| (b as i8) > b' ' as i8)
        .expect("a byte that is not white space");
    let stream = ArchiveStream::new(&text[start + 1..], &text[start..=start], Vec::new());
    drain(stream)
}

/// The objects and posts of `stream` (`read_with_posts`).
fn drain<R: std::io::Read>(mut stream: ArchiveStream<R, Vec<Post>>) -> (Objects, Vec<String>) {
    let mut out = Vec::new();
    let lines = |stream: &mut ArchiveStream<R, Vec<Post>>| {
        let posts = std::mem::take(stream.posts());
        posts
            .iter()
            .map(|post| String::from_utf8_lossy(&post.line(None)).into_owned())
            .collect()
    };
    loop {
        match stream.next_archive() {
            Ok(Some(archive)) => out.push(Ok(archive)),
            Ok(None) => return (out, lines(&mut stream)),
            Err(error) => {
                out.push(Err(error));
                let posts = lines(&mut stream);
                assert!(matches!(stream.next_archive(), Ok(None)) && stream.posts().is_empty());
                return (out, posts);
            }
        }
    }
}

/// NCBI's post of the byte `byte` of a VisibleString at `line`, in the piece `piece`
/// (`ReplaceVisibleChar`), on the string at `trace`.
fn bad_char(trace: &str, byte: u8, line: usize, piece: &str) -> String {
    format!("Error: (803.7) {trace}\nBad char [0x{byte:X}] in string at line {line}\n{piece}\n")
}

/// The line of the first `needle` in `text` (from 1).
fn line_of(text: &str, needle: &str) -> usize {
    let at = text.find(needle).unwrap_or_else(|| panic!("{needle}"));
    text[..at].matches('\n').count() + 1
}

/// The trace of a query's descriptor (`Seqdesc`) of a `bioseq-set` query.
const DESCR: &str =
    "Blast4-archive.request.body.queue-search.queries.bioseq-set.seq-set.E.seq.descr.E";

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

// A byte outside 0x20..0x7E in a VisibleString is read as `#` with NCBI's post of it, its stack
// trace, line and piece (`ReplaceVisibleChar`: one post per byte, each with the piece as read, up
// to the end of the line); the posts of an object come also when NCBI's reader then fails in the
// object, which ends the input (item 9(a) of line RG-KB; probe agent1/p2 `cm_bad`).
#[test]
fn invisible_bytes_read_as_hash_with_posts() {
    let text = String::from_utf8(ONE.to_vec()).unwrap();
    let invisible = text.replacen(
        "title \"NP_414783.1 protein YkfF",
        "title \"NP_414783.1 protein Yk\u{e9}F",
        1,
    );
    assert_ne!(invisible, text);
    let (read, posts) = read_with_posts(invisible.as_bytes());
    assert_eq!(read.len(), 1);
    let archive = read[0].as_ref().unwrap();
    let RequestBody::QueueSearch(search) = &archive.body else {
        panic!("queue-search");
    };
    let Queries::BioseqSet(entries) = &search.queries else {
        panic!("bioseq-set");
    };
    let SeqEntry::Seq(query) = &entries[0] else {
        panic!("seq");
    };
    assert!(query
        .title
        .as_ref()
        .unwrap()
        .starts_with(b"NP_414783.1 protein Yk##F [Escherichia coli str. K-12 substr. MG1655]"));
    let line = line_of(&invisible, "title \"NP_414783.1");
    let piece = "NP_414783.1 protein Yk\u{e9}F [Escherichia coli str. K-12";
    let piece = String::from_utf8_lossy(piece.as_bytes()).into_owned();
    let trace = format!("{DESCR}.title");
    let expected: Vec<String> = [0xc3, 0xa9]
        .iter()
        .map(|&byte| bad_char(&trace, byte, line, &piece))
        .collect();
    assert_eq!(posts, expected);
    let then_bad = invisible.replace("gapped FALSE", "gapped MAYBE");
    assert_ne!(then_bad, invisible);
    let (read, posts) = read_with_posts(then_bad.as_bytes());
    assert!(read.is_empty());
    assert_eq!(posts, expected);
    // REALs that LOSAT once refused are read (line SW, SW8): a mantissa of more than 18 digits,
    // as NCBI reads it, and an exponent that overflows NCBI's `int` (undefined behaviour), in
    // LOSAT's widened arithmetic (0, which the word form's `fDecimalPosixFinite` makes DBL_MIN;
    // NCBI's wrapped exponent gives DBL_MAX, PD-LOSAT-NCBI-DEFECTS exception 11).
    for real in ["{ 5047423181640891234, 10, -61 }", "1.55e-2147483647"] {
        let read_real = text.replace("{ 504742318164089, 10, -57 }", real);
        assert_ne!(read_real, text);
        let read = read_all(read_real.as_bytes());
        assert_eq!(read.len(), 1);
        let archive = read[0].as_ref().unwrap_or_else(|e| panic!("{real}: {e}"));
        let e_value = archive.results.alignments.as_ref().unwrap()[0].e_value;
        if real == "1.55e-2147483647" {
            assert_eq!(e_value, Some(f64::MIN_POSITIVE));
        }
    }
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
    let mut stream = ArchiveStream::new(OneByte(&text[1..]), &text[..1], Vec::new());
    stream.fill_one_byte_at_a_time();
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

// `CSeq_align::GetSeq_id` for index 0 (Seq_align.cpp:339-358, cited at `std_segs`) as read: segments
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

/// `read_with_posts` of `text` with the input given one byte per read: every value is read
/// after the stream's buffer was filled again.
fn read_all_by_byte(text: &[u8]) -> (Objects, Vec<String>) {
    use crate::report::serial::asn_text_reader::tests::OneByte;
    let start = text
        .iter()
        .position(|&b| (b as i8) > b' ' as i8)
        .expect("a byte that is not white space");
    let mut stream = ArchiveStream::new(
        OneByte(&text[start + 1..]),
        &text[start..=start],
        Vec::new(),
    );
    stream.fill_one_byte_at_a_time();
    drain(stream)
}

/// The archive fixtures read with the input given one byte per read are read as with the input
/// given in blocks: the same objects, errors, refusals and posts, whole and cut at 24 lengths.
#[test]
fn archives_read_one_byte_per_read() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/rpsbproc/arch");
    let mut paths: Vec<_> = std::fs::read_dir(&dir)
        .unwrap()
        .chain(std::fs::read_dir(dir.join("bad")).unwrap())
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.is_file())
        .collect();
    paths.sort();
    let mut read = 0;
    for path in paths {
        let data = std::fs::read(&path).unwrap();
        let text = if path.extension().is_some_and(|e| e == "gz") {
            let mut out = Vec::new();
            std::io::Read::read_to_end(&mut flate2::read::MultiGzDecoder::new(&data[..]), &mut out)
                .unwrap();
            out
        } else {
            data
        };
        let Some(start) = text.iter().position(|&b| (b as i8) > b' ' as i8) else {
            continue;
        };
        let mut cuts: Vec<usize> = (1..24).map(|k| text.len() * k / 24).collect();
        cuts.push(text.len());
        for cut in cuts.into_iter().filter(|&cut| cut > start) {
            let blocks = format!("{:?}", read_with_posts(&text[..cut]));
            let bytes = format!("{:?}", read_all_by_byte(&text[..cut]));
            assert_eq!(blocks, bytes, "{} cut {cut}", path.display());
            read += 1;
        }
    }
    assert!(read > 1000, "{read}");
}

// The values that the grammar does not take are read by the typed reader below the grammar's
// trail, as NCBI's reader reads them: their strings post with their traces (NCBI rpsbproc 0.5.1,
// probe agent1/p2 of item 9(a): `cm_e9`, `nm_7f`, `an_e9`, `bv_e9`), and a value that does not
// read as its type ends the input.
#[test]
fn values_the_grammar_does_not_take_post_their_strings() {
    let text = String::from_utf8(ONE.to_vec()).unwrap();
    let title = "              title \"NP_414783.1";
    let comment = text.replacen(
        title,
        &format!("              comment \"c\u{7f}d\",\n              name \"\u{1}\",\n{title}"),
        1,
    );
    let (read, posts) = read_with_posts(comment.as_bytes());
    assert_eq!(read.len(), 1);
    assert!(read[0].is_ok());
    let line = line_of(&comment, "comment \"");
    assert_eq!(
        posts,
        [
            bad_char(&format!("{DESCR}.comment"), 0x7f, line, "c\u{7f}d"),
            bad_char(&format!("{DESCR}.name"), 0x1, line + 1, "\u{1}"),
        ]
    );
    // A Seq-annot of the query, with a descriptor of its own.
    let inst_end = "            }\n          }\n        }\n      },";
    assert_eq!(text.matches(inst_end).count(), 1);
    let annot = text.replacen(
        inst_end,
        "            },\n            annot {\n              {\n                name \"a\u{e9}\",\n                desc {\n                  comment \"\u{1}\"\n                },\n                data locs {\n                }\n              }\n            }\n          }\n        }\n      },",
        1,
    );
    let (read, posts) = read_with_posts(annot.as_bytes());
    assert_eq!(read.len(), 1);
    assert!(read[0].is_ok());
    let line = line_of(&annot, "name \"a");
    let seq = "Blast4-archive.request.body.queue-search.queries.bioseq-set.seq-set.E.seq.annot.E";
    let piece = String::from_utf8_lossy("a\u{e9}".as_bytes()).into_owned();
    assert_eq!(
        posts,
        [
            bad_char(&format!("{seq}.name"), 0xc3, line, &piece),
            bad_char(&format!("{seq}.name"), 0xa9, line, &piece),
            bad_char(&format!("{seq}.desc.E.comment"), 0x1, line + 2, "\u{1}"),
        ]
    );
    // A `Blast4-value` variant that rpsbproc does not read.
    let options = "      algorithm-options {\n";
    assert_eq!(text.matches(options).count(), 1);
    let value = text.replacen(
        options,
        &format!("{options}        {{\n          name \"Extra\",\n          value string-list {{ \"ok\", \"x\u{7f}\" }}\n        }},\n"),
        1,
    );
    let (read, posts) = read_with_posts(value.as_bytes());
    assert_eq!(read.len(), 1);
    assert!(read[0].is_ok());
    assert_eq!(
        posts,
        [bad_char(
            "Blast4-archive.request.body.queue-search.algorithm-options.E.value.string-list.E",
            0x7f,
            line_of(&value, "value string-list"),
            "x\u{7f}"
        )]
    );
    // A Seqdesc variant that does not read as its type: NCBI's reader fails, the input ends.
    let bad = text.replacen(title, &format!("              comment 5,\n{title}"), 1);
    let (read, posts) = read_with_posts(bad.as_bytes());
    assert!(read.is_empty());
    assert!(posts.is_empty());
}

// Each object is read by a new object stream of NCBI's, whose line count starts at 1 at the byte
// after the last one the object before it took (the line end after its `}` is the second
// object's line 1; NCBI rpsbproc 0.5.1, probe agent1/p2 `two_e9`), and its posts come after the
// object before it was read.
#[test]
fn the_lines_of_a_second_object_start_again() {
    let text = String::from_utf8(ONE.to_vec()).unwrap();
    let second = text.replacen(
        "title \"NP_414783.1 protein Yk",
        "title \"NP_414783.1 protein \u{7f}k",
        1,
    );
    assert_ne!(second, text);
    let line = line_of(&second, "\u{7f}k");
    let (read, posts) = read_with_posts(format!("{text}{second}").as_bytes());
    assert_eq!(read.len(), 2);
    assert_eq!(
        posts,
        [bad_char(
            &format!("{DESCR}.title"),
            0x7f,
            line + 1,
            "NP_414783.1 protein \u{7f}kfF [Escherichia coli str. K-12"
        )]
    );
}

// The names that the grammar gives its members, variants and enumerated values (its trail, the
// traces of its posts and the types of the typed reads) are the type table's.
#[test]
fn the_grammar_names_are_the_tables() {
    use crate::report::serial::asn_types::{self, table::named, Kind};
    let names =
        |ty: u16| -> Vec<&str> { asn_types::members(ty).iter().map(|m| m.name()).collect() };
    assert_eq!(names(named::BLAST4_REQUEST_BODY), REQUEST_BODY);
    assert_eq!(names(named::BLAST4_VALUE), BLAST4_VALUE);
    assert_eq!(names(named::SEQDESC), SEQDESC);
    assert_eq!(names(named::SEQ_DATA), SEQ_DATA);
    assert_eq!(names(named::SEQ_LOC), SEQ_LOC);
    assert_eq!(names(seq_align_segs()), SEQ_ALIGN_SEGS);
    let member = |ty: u16, name: &str| asn_types::members(ty)[asn_types::member_index(ty, name)].ty;
    let values = |ty: u16| -> Vec<(&str, i32)> {
        match asn_types::resolved(ty).kind {
            Kind::Enumerated { values, .. } => values
                .values()
                .iter()
                .map(|v| (v.name(), v.value))
                .collect(),
            _ => panic!("{ty} is not ENUMERATED"),
        }
    };
    assert_eq!(values(named::BLAST4_STRAND_TYPE), BLAST4_STRAND_TYPE);
    assert_eq!(values(named::BLAST4_FRAME_TYPE), BLAST4_FRAME_TYPE);
    assert_eq!(values(named::NA_STRAND), NA_STRAND);
    assert_eq!(values(member(named::BIOSEQ_SET, "class")), BIOSEQ_SET_CLASS);
    assert_eq!(
        values(member(named::SEQ_INST, "topology")),
        SEQ_INST_TOPOLOGY
    );
    assert_eq!(values(member(named::SEQ_INST, "strand")), SEQ_INST_STRAND);
    assert_eq!(values(member(named::SEQ_INST, "repr")), SEQ_INST_REPR);
    assert_eq!(values(member(named::SEQ_INST, "mol")), SEQ_INST_MOL);
    assert_eq!(values(member(named::SEQ_ALIGN, "type")), SEQ_ALIGN_TYPE);
}
