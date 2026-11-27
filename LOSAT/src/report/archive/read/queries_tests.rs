//! The query sets of probe `F/f2/tools/qprobe.py` (NCBI BLAST+ 2.17.0 on ws, 2026-11-20): each
//! case's Bioseq-set as NCBI's strategy held it (`<case>.asn`), and what NCBI did with it: the
//! FASTA text of its temporary file (`<case>.fasta`) or the exception it reported
//! (`<case>.exception`, the report's lines in the probe build's tree). `cases.tsv` says which
//! cases LOSAT refuses and why.

use super::*;
use crate::report::serial::asn_input::Names;
use crate::report::serial::asn_types::{read_document, table::named};

const FIXTURES: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/strategy_queries"
);

/// The tree that the probe build's exception reports name.
const PROBE_BUILD_ROOT: &str = "/opt/conda/conda-bld/blast_1787399212531/work/c++";

fn read(name: &str) -> Option<Vec<u8>> {
    std::fs::read(format!("{FIXTURES}/{name}")).ok()
}

/// The text of the query set `text` held whole, after the check that the text written to a file a
/// block at a time is the same, or the same error.
fn render(text: &[u8]) -> Result<Vec<u8>, RenderError> {
    let (value, posts) = read_document(text, named::BIOSEQ_SET, Names::STRATEGY);
    assert!(posts.is_empty(), "{posts:?}");
    let value = value.unwrap_or_else(|error| panic!("{error:?}"));
    let set = V::new(named::BIOSEQ_SET, &value);
    let held = render_queries(set, Text::held()).map(bytes);
    let file = crate::blastinput::input_files::temporary_output("query").expect("a file");
    let written = render_queries(set, Text::file(file)).map(bytes);
    assert_eq!(format!("{written:?}"), format!("{held:?}"));
    held
}

/// The bytes of a text written.
fn bytes(rendered: Rendered) -> Vec<u8> {
    match rendered {
        Rendered::Bytes(bytes) => bytes,
        Rendered::File(mut file) => {
            let mut bytes = Vec::new();
            std::io::Read::read_to_end(&mut file, &mut bytes).expect("the file read");
            bytes
        }
    }
}

/// The letters of `bioseq` as a row of the tabular formats over all of them on the plus strand
/// reads them (`seq_vector::RowVector`: one read from the first letter to the last), or its end.
fn bioseq_letters(bioseq: Bioseq<'_>) -> Result<Vec<u8>, RenderError> {
    let map = seq_vector::SeqMap::of(bioseq)?;
    let length = map.length()?;
    let mut out = Vec::new();
    if length > 0 {
        seq_vector::RowVector::new(&map, true).read(0, length - 1)?;
        map.letters(0, length, &mut out);
    }
    Ok(out)
}

/// [`bioseq_letters`]' end without its letters.
fn check_bioseq_letters(bioseq: Bioseq<'_>) -> Result<(), RenderError> {
    bioseq_letters(bioseq).map(|_| ())
}

#[test]
fn every_probed_query_set_is_written_as_ncbi_wrote_it() {
    let manifest = String::from_utf8(read("cases.tsv").expect("cases.tsv")).unwrap();
    let mut checked = 0;
    for line in manifest.lines().filter(|line| !line.starts_with('#')) {
        let mut fields = line.split('\t');
        let case = fields.next().unwrap();
        let losat = fields.next().unwrap();
        let text = read(&format!("{case}.asn")).unwrap_or_else(|| panic!("{case}.asn"));
        let got = render(&text);
        match (
            losat,
            read(&format!("{case}.fasta")),
            read(&format!("{case}.exception")),
        ) {
            ("refused", _, _) => match got {
                Err(RenderError::Unsupported(part)) => assert!(!part.is_empty(), "{case}"),
                other => panic!("{case}: {other:?}"),
            },
            ("same", Some(fasta), None) => match got {
                Ok(bytes) => assert!(
                    bytes == fasta,
                    "{case}:\n{}\n{}",
                    String::from_utf8_lossy(&bytes),
                    String::from_utf8_lossy(&fasta)
                ),
                Err(error) => panic!("{case}: {error:?}"),
            },
            ("same", None, Some(exception)) => match got {
                Err(RenderError::Exception(e)) => assert_eq!(
                    String::from_utf8_lossy(&e.what_in(PROBE_BUILD_ROOT)),
                    crate::wording::losat_text(&String::from_utf8(exception).unwrap()),
                    "{case}"
                ),
                other => panic!("{case}: {other:?}"),
            },
            _ => panic!("{case}: no expectation"),
        }
        checked += 1;
    }
    assert!(checked >= 100, "{checked} cases");
}

#[test]
fn fills_of_a_long_seq_data_end_at_1024_letter_steps() {
    // `reference out of range` names the end of the first fill past the data: 2000 letters of
    // which 1500 are held end their fills at 1024 and 2000; 3000 at 1024, 2048.
    for (length, end) in [(2000, 2000), (3000, 2048), (1501, 1501)] {
        let text = format!(
            "Bioseq-set ::= {{ seq-set {{ seq {{ id {{ local str \"q\" }}, inst {{ repr raw, mol na, length {length}, seq-data iupacna \"{}\" }} }} }} }}",
            "ACGT".repeat(375)
        );
        match render(text.as_bytes()) {
            Err(RenderError::Exception(e)) => assert_eq!(
                e.get_msg(),
                format!("reference out of range of Seq-inst data: {end}").as_bytes()
            ),
            other => panic!("{length}: {other:?}"),
        }
    }
}

#[test]
fn a_segment_that_passes_the_length_is_filled_up_to_it() {
    // A literal of 3000 letters that holds 1500, in a Bioseq of `length` letters: its fills end at
    // 1024 and at `length` (`dl_cut_short`); only a fill past the 1500 letters throws.
    let text = |length: u32| {
        format!(
            "Bioseq-set ::= {{ seq-set {{ seq {{ id {{ local str \"q\" }}, inst {{ repr delta, mol na, length {length}, ext delta {{ literal {{ length 3000, seq-data iupacna \"{}\" }} }} }} }} }} }}",
            "ACGT".repeat(375)
        )
    };
    for length in [1024, 1500] {
        let fasta = render(text(length).as_bytes()).unwrap();
        let letters = fasta.iter().filter(|b| b"ACGT".contains(b)).count();
        assert_eq!(letters, length as usize);
    }
    match render(text(1600).as_bytes()) {
        Err(RenderError::Exception(e)) => {
            assert_eq!(
                e.get_msg(),
                b"reference out of range of Seq-inst data: 1600"
            )
        }
        other => panic!("{other:?}"),
    }
}

/// Sequences of more letters than the blocks in which they are converted and written are the
/// lines of the letters written whole: a gap, letters copied and letters converted, held or
/// written to a file (`render` checks the two).
#[test]
fn long_sequences_are_written_in_blocks_as_whole() {
    let gap = 3 * LETTER_BLOCK + 5;
    let packed = FILE_BLOCK + 7;
    let copied = 2 * LETTER_BLOCK + 3;
    let set = format!(
        "Bioseq-set ::= {{ seq-set {{ \
         seq {{ id {{ local str \"Query_1\" }}, descr {{ title \"v\" }}, inst {{ repr virtual, mol na, length {gap} }} }}, \
         seq {{ id {{ local str \"Query_2\" }}, descr {{ title \"p\" }}, inst {{ repr raw, mol na, length {}, seq-data ncbi2na '{}'H }} }}, \
         seq {{ id {{ local str \"Query_3\" }}, descr {{ title \"c\" }}, inst {{ repr raw, mol aa, length {copied}, seq-data ncbieaa \"{}\" }} }} }} }}",
        4 * packed,
        "1B".repeat(packed),
        "M".repeat(copied),
    );
    let text = render(set.as_bytes()).unwrap_or_else(|error| panic!("{error:?}"));
    let mut expected = b">lcl|Query_1 v\n".to_vec();
    crate::sequence::write_fasta_lines(&vec![b'N'; gap], WIDTH, &mut expected);
    expected.extend_from_slice(b">lcl|Query_2 p\n");
    crate::sequence::write_fasta_lines(&b"ACGT".repeat(packed), WIDTH, &mut expected);
    expected.extend_from_slice(b">lcl|Query_3 c\n");
    crate::sequence::write_fasta_lines(&vec![b'M'; copied], WIDTH, &mut expected);
    assert!(expected.len() > 4 * FILE_BLOCK);
    assert!(text == expected, "the long sequences' lines");
}

/// A Bioseq that declares 2^32 - 2 letters and holds 8 gives the error of its data (the fill past
/// its 8 letters, `reference out of range`) after room for the letters that its Seq-data hold,
/// not for the length it declares (review R1 of PR 1683, finding 5: no reserve of 4 GiB, which
/// aborts under a memory limit and overflows the capacity on wasm32); one that declares 2^32 - 1
/// (`kInvalidSeqPos`) gives the sequence map's `Invalid sequence length`, as a delta literal of
/// that length does where the read reaches it; a delta's room is its gaps and the letters that
/// its Seq-data hold.
#[test]
fn a_declared_length_past_the_data_reserves_the_letters_held() {
    let read = |text: &str| {
        let (value, posts) = read_document(text.as_bytes(), named::BIOSEQ_SET, Names::STRATEGY);
        assert!(posts.is_empty(), "{posts:?}");
        value.unwrap_or_else(|error| panic!("{error:?}"))
    };
    let raw = |length: u32| {
        format!(
            "Bioseq-set ::= {{ seq-set {{ seq {{ id {{ local str \"q\" }}, inst {{ repr raw, mol na, \
             length {length}, seq-data iupacna \"ACGTACGT\" }} }} }} }}"
        )
    };
    // The message of the read's exception, which the read that only checks gives too.
    let message = |bioseq: Bioseq<'_>| {
        let checked = match check_bioseq_letters(bioseq) {
            Err(RenderError::Exception(e)) => e.get_msg().to_vec(),
            other => panic!("{other:?}"),
        };
        match bioseq_letters(bioseq) {
            Err(RenderError::Exception(e)) => {
                assert_eq!(e.get_msg(), checked);
                checked
            }
            other => panic!("{other:?}"),
        }
    };
    // 4294967294 letters declared, 8 held: the room of the 8, and NCBI's throw of the fill that
    // passes them (probe `of6_qlen_big` of agent E3, 2026-11-27).
    let value = read(&raw(u32::MAX - 1));
    let set = V::new(named::BIOSEQ_SET, &value);
    let entry = set.field("seq-set").item(0).expect("one entry");
    let bioseq = Bioseq::of(entry.choice().1);
    let pieces = sequence_map(bioseq.value().field("inst"), u32::MAX - 1).expect("a map");
    assert_eq!(written_letters(&pieces, u32::MAX - 1), 8);
    let out_of_range = b"reference out of range of Seq-inst data: 1024".to_vec();
    assert_eq!(message(bioseq), out_of_range);
    match render(raw(u32::MAX - 1).as_bytes()) {
        Err(RenderError::Exception(e)) => assert_eq!(e.get_msg(), out_of_range),
        other => panic!("{other:?}"),
    }
    // 4294967295, `kInvalidSeqPos`: no room, and the throw of the sequence map
    // (`of6_qlen_huge`).
    let value = read(&raw(u32::MAX));
    let set = V::new(named::BIOSEQ_SET, &value);
    let entry = set.field("seq-set").item(0).expect("one entry");
    let bioseq = Bioseq::of(entry.choice().1);
    let pieces = sequence_map(bioseq.value().field("inst"), u32::MAX).expect("a map");
    assert_eq!(written_letters(&pieces, u32::MAX), 0);
    assert_eq!(message(bioseq), b"Invalid sequence length");
    // A gap of 5 and a literal that holds 1500 of its 3000 letters: the room of 1505; a literal
    // of `kInvalidSeqPos` letters after the gap throws when the read reaches it.
    let delta = |second: u32, letters: usize| {
        format!(
            "Bioseq-set ::= {{ seq-set {{ seq {{ id {{ local str \"q\" }}, inst {{ repr delta, mol na, \
             length 1600, ext delta {{ literal {{ length 5 }}, literal {{ length {second}, seq-data \
             iupacna \"{}\" }} }} }} }} }} }}",
            "ACGT".repeat(letters / 4)
        )
    };
    let value = read(&delta(3000, 1500));
    let set = V::new(named::BIOSEQ_SET, &value);
    let entry = set.field("seq-set").item(0).expect("one entry");
    let bioseq = Bioseq::of(entry.choice().1);
    let pieces = sequence_map(bioseq.value().field("inst"), 1600).expect("a map");
    assert_eq!(written_letters(&pieces, 1600), 5 + 1500);
    // A delta that reads: its letters, and the read that only checks passes.
    let value = read(&delta(1596, 1596));
    let set = V::new(named::BIOSEQ_SET, &value);
    let entry = set.field("seq-set").item(0).expect("one entry");
    let bioseq = Bioseq::of(entry.choice().1);
    let letters = bioseq_letters(bioseq).unwrap_or_else(|error| panic!("{error:?}"));
    assert_eq!(
        letters,
        [b"NNNNN".as_slice(), &b"ACGT".repeat(399)[..1595]].concat()
    );
    assert!(check_bioseq_letters(bioseq).is_ok());
    let value = read(&delta(u32::MAX, 8));
    let set = V::new(named::BIOSEQ_SET, &value);
    let entry = set.field("seq-set").item(0).expect("one entry");
    let bioseq = Bioseq::of(entry.choice().1);
    let pieces = sequence_map(bioseq.value().field("inst"), 1600).expect("a map");
    assert_eq!(written_letters(&pieces, 1600), 5);
    assert_eq!(message(bioseq), b"Invalid sequence length");
}

/// The read of a `-subject` Bioseq by the formatter's setup (`check_setup_letters`): NCBI's
/// exceptions of the read from the first letter, and LOSAT's refusals of what the setup reads
/// otherwise (a residue past NCBIstdaa's 28 codes, a coding to convert); the row reads of the
/// tabular formats (`bioseq_letters`) refuse an error that depends on the row's range and keep
/// those of any read.
#[test]
fn the_setup_reads_a_subject_whole_and_the_rows_keep_only_range_free_errors() {
    let bioseq_of = |inst: &str| {
        let text = format!(
            "Bioseq-set ::= {{ seq-set {{ seq {{ id {{ local str \"s\" }}, inst {{ {inst} }} }} }} }}"
        );
        let (value, posts) = read_document(text.as_bytes(), named::BIOSEQ_SET, Names::STRATEGY);
        assert!(posts.is_empty(), "{posts:?}");
        value.unwrap_or_else(|error| panic!("{error:?}"))
    };
    let check = |inst: &str, protein: bool| {
        let value = bioseq_of(inst);
        let set = V::new(named::BIOSEQ_SET, &value);
        let entry = set.field("seq-set").item(0).expect("one entry");
        let bioseq = Bioseq::of(entry.choice().1);
        let setup = match check_setup_letters(bioseq, protein) {
            Ok(SubjectSetup::Read) => "ok".to_string(),
            Ok(SubjectSetup::Empty) => "empty".to_string(),
            Err(RenderError::Exception(e)) => String::from_utf8_lossy(e.get_msg()).into_owned(),
            Err(RenderError::Unsupported(part)) => format!("refused: {part}"),
            Err(other) => panic!("{other:?}"),
        };
        let rows = match check_bioseq_letters(bioseq) {
            Ok(()) => "ok".to_string(),
            Err(RenderError::Exception(e)) => String::from_utf8_lossy(e.get_msg()).into_owned(),
            Err(RenderError::Unsupported(part)) => format!("refused: {part}"),
            Err(other) => panic!("{other:?}"),
        };
        (setup, rows)
    };
    let out_of_range = "reference out of range of Seq-inst data: 1024";
    let range = "refused: letters that do not all read (NCBI's error depends on the row's range)";
    assert_eq!(
        check(
            "repr raw, mol na, length 8, seq-data iupacna \"ACGTACGT\"",
            false
        ),
        ("ok".into(), "ok".into())
    );
    assert_eq!(
        check(
            "repr raw, mol na, length 4294967295, seq-data iupacna \"ACGT\"",
            false
        ),
        (
            "Invalid sequence length".into(),
            "Invalid sequence length".into()
        )
    );
    assert_eq!(
        check(
            "repr raw, mol na, length 1500, seq-data iupacna \"ACGT\"",
            false
        ),
        (out_of_range.into(), range.into())
    );
    assert_eq!(
        check(
            "repr raw, mol na, length 8, seq-data ncbi2na '1B1B'H",
            false
        ),
        ("ok".into(), "ok".into())
    );
    // A protein Bioseq in a search of nucleotide subjects: the setup's NCBI4na is not probed;
    // the rows read it in ncbieaa.
    let (setup, rows) = check(
        "repr raw, mol aa, length 4, seq-data ncbieaa \"MKLV\"",
        false,
    );
    assert!(setup.starts_with("refused: a Seq-data coding"), "{setup}");
    assert_eq!(rows, "ok");
    assert_eq!(
        check(
            "repr raw, mol aa, length 4, seq-data ncbieaa \"MKLV\"",
            true
        ),
        ("ok".into(), "ok".into())
    );
    // NCBIstdaa codes 0 (`-`) and 27 (`J`) are residues; 28 is none.
    assert_eq!(
        check(
            "repr raw, mol aa, length 2, seq-data ncbistdaa '001B'H",
            true
        ),
        ("ok".into(), "ok".into())
    );
    let (setup, rows) = check(
        "repr raw, mol aa, length 2, seq-data ncbistdaa '011C'H",
        true,
    );
    assert!(setup.starts_with("refused: a residue"), "{setup}");
    assert_eq!(rows, "ok");
    // A delta whose second literal holds fewer letters than it declares: the setup's read from
    // the first letter throws at its fill; a row's range may not reach it.
    let (setup, rows) = check(
        "repr delta, mol na, length 9, ext delta { literal { length 4, seq-data iupacna \"ACGT\" }, \
         literal { length 5, seq-data iupacna \"AC\" } }",
        false,
    );
    assert_eq!(setup, "reference out of range of Seq-inst data: 5");
    assert_eq!(rows, range);
}
