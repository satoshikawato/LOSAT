//! Tests of the five scanners besides Clustal's, the defline modifiers and the title. The expected
//! answers are NCBI BLAST+ 2.17.0's (psiblast `-in_msa`, the probe corpus of line PS item PS6b,
//! `PS/ps6b`, cases named in the comments): the read of an alignment succeeds or fails as
//! `CAlnReader` does, whatever the message of its report.
use super::*;

const Q: &[u8] = b"MKVLA-GT";
const S1: &[u8] = b"MKVIAAGT";
const S2: &[u8] = b"-KVLAAGT";

fn rows(text: &str) -> Result<Vec<Vec<u8>>, AlnError> {
    read_alignment(text.as_bytes()).map(|alignment| alignment.rows)
}

fn three() -> Vec<Vec<u8>> {
    vec![Q.to_vec(), S1.to_vec(), S2.to_vec()]
}

fn fails(text: &str) -> bool {
    matches!(rows(text), Err(AlnError::NotAvailable))
}

const FASTA: &str = ">q1\nMKVLA-GT\n>s1\nMKVIAAGT\n>s2\n-KVLAAGT\n";
const PHYLIP: &str = " 3 8\nq1         MKVLA-GT\ns1         MKVIAAGT\ns2         -KVLAAGT\n";
const NEXUS: &str = "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=3 NCHAR=8;\nFORMAT DATATYPE=PROTEIN GAP=-;\nMATRIX\nq1 MKVLA-GT\ns1 MKVIAAGT\ns2 -KVLAAGT\n;\nEND;\n";
const SEQUIN: &str = "\n                   10         20         30         40         50\nq1         MKVLA-GT\ns1         MKVIAAGT\ns2         -KVLAAGT\n\n//\n";
// A MultAlin file's gap is the dot (`f6b.ma.dashgap` fails, `f6b.ma.dotgap` reads).
const MULTALIGN: &str = "//\n\n1 8\nq1         MKVLA.GT\ns1         MKVIAAGT\ns2         .KVLAAGT\nConsensus xxxxxxxx\n\n";

// Each format reads the same three rows (cases `f6b.fa.basic`, `f6b.ph.basic`, `f6b.nx.basic`,
// `f6b.sq.basic`, `f6b.ma.basic`).
#[test]
fn every_format_reads_the_same_rows() {
    for (name, text) in [
        ("fasta", FASTA),
        ("phylip", PHYLIP),
        ("nexus", NEXUS),
        ("sequin", SEQUIN),
    ] {
        let read = rows(text).unwrap_or_else(|error| panic!("{name}: {error:?}"));
        assert_eq!(read, three(), "{name}");
    }
    let dotted = rows(MULTALIGN).unwrap();
    assert_eq!(dotted[0], b"MKVLA.GT");
    assert_eq!(dotted[2], b".KVLAAGT");
}

#[test]
fn the_guess_of_each_format() {
    let guess = |text: &str| guess_format(&mut PeekAheadStream::new(text.as_bytes()));
    assert_eq!(guess(FASTA), AlignFormat::FastaGap);
    assert_eq!(guess(PHYLIP), AlignFormat::Phylip);
    assert_eq!(guess(NEXUS), AlignFormat::Nexus);
    assert_eq!(guess(SEQUIN), AlignFormat::Sequin);
    assert_eq!(guess(MULTALIGN), AlignFormat::MultAlign);
}

// FASTA-gap: a row is the lines up to the blank line or the next `>`, blanks inside a line are
// dropped (`f6b.fa.wrap7`, `f6b.fa.blankin`), a repeated id (also by case) fails (`f6b.fa.dup`,
// `f6b.fa.dup_case`), and so does a row that is not made of the alphabet (`f6b.fa.badchar_*`).
#[test]
fn fasta_rows_and_failures() {
    assert_eq!(
        rows(">q1\nMKVL\nA-GT\n>s1\nMKVI\nAAGT\n>s2\n-KVL\nAAGT\n").unwrap(),
        three()
    );
    assert_eq!(
        rows(">q1\nMKV LA-GT\n>s1\nMKVIAAGT\n>s2\n-KVLAAGT\n").unwrap(),
        three()
    );
    assert_eq!(
        rows(">q1\nMKVLA-GT\n\n>s1\nMKVIAAGT\n\n>s2\n-KVLAAGT\n").unwrap(),
        three()
    );
    assert!(fails(">q1\nMKVLA-GT\n>q1\nMKVIAAGT\n>s2\n-KVLAAGT\n"));
    assert!(fails(">q1\nMKVLA-GT\n>Q1\nMKVIAAGT\n>s2\n-KVLAAGT\n"));
    assert!(fails(">q1\nMKVLA-GT\n>s1\nMKV1AAGT\n>s2\n-KVLAAGT\n"));
    assert!(fails(">q1\nMKVLA-GT\n>s1\nMKV.AAGT\n>s2\n-KVLAAGT\n"));
    assert!(fails(">q1\nMKVLA-GT\n"));
    // Rows without letters are read (the query length is then 0, `f6b.fa.empty_rows`); a first row
    // alone with data and the others empty is not (`f6b.fa.empty_first`).
    assert!(rows(">q1\n>s1\n>s2\n").is_ok());
    // One row of another length passes only where it has the lines of the first row.
    assert!(fails(">q1\nMKVLA-GT\n>s1\nMKVIAA\n>s2\n-KVLAAGT\n"));
}

// PHYLIP: the header gives the count and the length (`f6b.ph.hdr*`), the id is the first ten
// characters of its line and the rest of the line is the row (`f6b.ph.id10`), a wrong count or
// length fails (`f6b.ph.hdr05`..), a count of 0 crashes NCBI (SIGFPE, `f6b.ph.hdr09`), and a
// number that does not fit an int is an uncaught `CStringException` (`f6b.ph.hdr12`).
#[test]
fn phylip_rows_and_failures() {
    assert!(fails(
        " 4 8\nq1         MKVLA-GT\ns1         MKVIAAGT\ns2         -KVLAAGT\n"
    ));
    assert!(fails(
        " 3 9\nq1         MKVLA-GT\ns1         MKVIAAGT\ns2         -KVLAAGT\n"
    ));
    assert!(fails(" 1 8\nq1         MKVLA-GT\n"));
    assert!(matches!(
        rows(" 0 8\nq1         MKVLA-GT\ns1         MKVIAAGT\n"),
        Err(AlnError::Undefined(_))
    ));
    assert!(matches!(
        rows(" 99999999999 8\nq1         MKVLA-GT\ns1         MKVIAAGT\n"),
        Err(AlnError::Exception(_))
    ));
}

// NEXUS: DIMENSIONS and FORMAT, the rows of MATRIX (`f6b.nx.*`).
#[test]
fn nexus_rows_and_failures() {
    // A GAP character of the FORMAT command is the gap of the alignment (`f6b.nx.dotgap`).
    let dotted = NEXUS
        .replace("GAP=-", "GAP=.")
        .replace("MKVLA-GT", "MKVLA.GT")
        .replace("-KVLAAGT", ".KVLAAGT");
    let read = rows(&dotted).unwrap();
    assert_eq!(read[0], b"MKVLA.GT");
    assert_eq!(read[2], b".KVLAAGT");
    assert!(fails(&NEXUS.replace("NTAX=3", "NTAX=4")));
    assert!(fails(&NEXUS.replace("NCHAR=8", "NCHAR=9")));
    assert!(fails(&NEXUS.replace("MATRIX\n", "")));
    // Comments in square brackets are cut out (`f6b.nx.comment_in_data`).
    assert_eq!(
        rows(&NEXUS.replace("MKVLA-GT", "MKV[hello]LA-GT")).unwrap(),
        three()
    );
    // The data type is not read (`f6b.nx.dna`): a DNA alignment reads, as protein letters.
    assert!(rows(&NEXUS.replace("PROTEIN", "DNA")).is_ok());
}

// Sequin and MultAlin (`f6b.sq.*`, `f6b.ma.*`).
#[test]
fn sequin_and_multalign_failures() {
    // The `//` that ends a Sequin file is optional (`f6b.sq.no_term`); a digit in a row fails
    // (`f6b.sq.digit`).
    assert!(rows(&SEQUIN.replace("\n//\n", "\n")).is_ok());
    assert!(fails(&SEQUIN.replace("MKVIAAGT", "MKV1AAGT")));
    assert!(fails(&MULTALIGN.replace(".KVLAAGT", "-KVLAAGT")));
}

// The defline of a row: the text outside the `[name=value]` modifiers is the title
// (`CTitleParser::Apply`, `f6b.md.*`).
fn title_of(defline: &str) -> Option<String> {
    let text = format!(">q1 {defline}\nMKVLA-GT\n>s1\nMKVIAAGT\n>s2\n-KVLAAGT\n");
    let alignment = read_alignment(text.as_bytes()).unwrap();
    alignment.bioseqs[0]
        .title
        .clone()
        .map(|title| String::from_utf8(title).unwrap())
}

#[test]
fn the_title_is_the_text_outside_the_modifiers() {
    assert_eq!(title_of("first title"), Some("first title".to_string()));
    assert_eq!(
        title_of("[organism=Homo sapiens] first title"),
        Some("first title".to_string())
    );
    assert_eq!(
        title_of("t1 [organism=Homo sapiens] t2"),
        Some("t1 t2".to_string())
    );
    assert_eq!(title_of("[organism=Homo sapiens]"), None);
    // Brackets without `=` are text, and the first one ends the search (`f6b.md.014`, `.040`).
    assert_eq!(title_of("[foo] t"), Some("[foo] t".to_string()));
    assert_eq!(
        title_of("t1 [a] t2 [c=d]"),
        Some("t1 [a] t2 [c=d]".to_string())
    );
    assert_eq!(title_of("t [foo"), Some("t [foo".to_string()));
    assert_eq!(title_of("   "), None);
}

#[test]
fn skipped_and_deprecated_modifiers_are_kept_in_the_title() {
    // Unrecognized (`f6b.md.006`..`.009`): the modifier goes back to the end of the title.
    assert_eq!(title_of("[foo=bar] t"), Some("t [foo=bar]".to_string()));
    assert_eq!(title_of("[foo=bar]"), Some("[foo=bar]".to_string()));
    assert_eq!(title_of("[ foo = bar ] t"), Some("t [foo=bar]".to_string()));
    // Deprecated (`f6b.md.051`): rejected modifiers come before the skipped ones.
    assert_eq!(
        title_of("[gene=a][dosage=b][foo=c] t"),
        Some("t [gene=a] [dosage=b] [foo=c]".to_string())
    );
    // A duplicate with the same value is dropped (`f6b.md.024`, `.029`).
    assert_eq!(title_of("[strain=a][Strain=A] t"), Some("t".to_string()));
}

#[test]
fn the_warnings_of_the_modifiers_come_in_ncbis_order() {
    let text = ">q1 [foo=1][gene=a][strain=a][strain=a][bar=2] t\nMKVLA-GT\n>s1\nMKVIAAGT\n>s2\n-KVLAAGT\n";
    let alignment = read_alignment(text.as_bytes()).unwrap();
    let warnings: Vec<String> = alignment
        .warnings
        .iter()
        .map(|w| String::from_utf8(w.clone()).unwrap())
        .collect();
    // The handler's own (deprecated, duplicate) first, then the adder's in the order of the
    // canonical names (bar before foo).
    assert_eq!(
        warnings,
        [
            "Use of the following modifier in a sequence file is discouraged and the information will be ignored: gene.",
            "Duplicated modifier value detected, ignoring duplicate, no action required: strain=a.",
            "Unrecognized modifier: bar.",
            "Unrecognized modifier: foo.",
        ]
    );
}

// An error of a modifier stops the read, after the warnings before it (`f6b.md.025`, `.064`,
// `.091`, `.m2.*`).
#[test]
fn modifier_errors_carry_their_message() {
    let error = |defline: &str| {
        let text = format!(">q1 {defline}\nMKVLA-GT\n>s1\nMKVIAAGT\n>s2\n-KVLAAGT\n");
        match read_alignment(text.as_bytes()) {
            Err(AlnError::Modifier(failure)) => (
                String::from_utf8(failure.message).unwrap(),
                failure.warnings.len(),
            ),
            other => panic!("{defline}: {other:?}"),
        }
    };
    assert_eq!(
        error("[organism=a][organism=b]").0,
        "Conflicting modifiers detected. Provide one modifier with one value for: organism."
    );
    assert_eq!(
        error("[pmid=abc] t").0,
        "Invalid value: pmid=abc. Expected integer value."
    );
    assert_eq!(
        error("[taxid=abc] t").0,
        "Invalid value: taxid=abc. Integer value expected."
    );
    assert_eq!(error("[location=bad] t").0, "Invalid value: location=bad.");
    assert_eq!(
        error("[environmental-sample=false] t").0,
        "Invalid value: environmental-sample=false."
    );
    assert_eq!(error("[focus=false] t").0, "Invalid value: focus=false.");
    // The warning of `foo` comes before the error of `location` (the names are in order).
    assert_eq!(
        error("[foo=1][location=bad] t"),
        ("Invalid value: location=bad.".to_string(), 1)
    );
    assert_eq!(error("[zzz=1][location=bad] t").1, 0);
}

// What the modifiers set besides the title is flagged (the run refuses where it would show),
// and a molecule of another kind than the row's is refused (`f6b.md.099`..`.109`).
#[test]
fn modifier_content_and_the_molecule() {
    let read = |defline: &str| {
        let text = format!(">q1 {defline}\nMKVLA-GT\n>s1\nMKVIAAGT\n>s2\n-KVLAAGT\n");
        read_alignment(text.as_bytes())
    };
    assert!(!read("[foo=bar] t").unwrap().bioseqs[0].content.has_any());
    assert!(read("[organism=a] t").unwrap().bioseqs[0].content.has_any());
    assert!(read("[protein=a] t").unwrap().bioseqs[0].content.has_any());
    assert!(read("[molecule=aa] t").is_ok());
    assert!(read("[molecule=dna] t").unwrap().bioseqs[0]
        .set_mol
        .is_some());
    assert!(read("[mol-type=genomic] t").unwrap().bioseqs[0]
        .set_mol
        .is_some());
    assert!(read("[molecule=other] t").unwrap().bioseqs[0]
        .set_mol
        .is_some());
    assert!(read("[molecule=aa] t").unwrap().bioseqs[0]
        .set_mol
        .is_none());
    // `molecule` comes after `mol-type` in the order of the names, so it decides (`f6b.md.111`).
    assert!(
        read("[mol-type=genomic][molecule=aa] t").unwrap().bioseqs[0]
            .set_mol
            .is_none()
    );
    // Only a protein row can change (`f6b.m2.dna_mol_aa`).
    let dna = format!(
        ">q1 [molecule=aa]\n{}\n>s1\n{}\n",
        "ACGT".repeat(15),
        "ACGA".repeat(15)
    );
    assert!(read_alignment(dna.as_bytes()).unwrap().bioseqs[0]
        .set_mol
        .is_none());
}
