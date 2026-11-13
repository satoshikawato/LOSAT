//! Tests of the organism report (outfmt 18) against NCBI BLAST+ 2.17.0's bytes (probe of line RG
//! T2, `remaining-gaps-20261112/T/t2/out`), with the taxonomy database `tests/seqdb/taxdb/of18`
//! (`of18.tsv`: the probe's `t2` database).

use super::*;
use crate::blastinput::app::{FormatChoice, ReportFormat};
use crate::report::tabular::{TabularSpec, TabularTaxonomy};
use crate::seqdb::taxdb::LazyTaxDb;
use crate::seqdb::SearchPath;
use std::path::{Path, PathBuf};

fn fixture(variant: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/seqdb/taxdb")
        .join(variant)
}

fn empty_cwd() -> PathBuf {
    static COUNTER: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "losat-tax-report-{}-{}",
        std::process::id(),
        COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// The lookups of `taxids` in the taxonomy database `variant` (`None`: no database).
fn names(variant: Option<&str>, taxids: &[i32]) -> Option<TabularNames> {
    let variant = variant?;
    let cwd = empty_cwd();
    let path = SearchPath::new(&cwd, Some(fixture(variant).as_os_str()), None).unwrap();
    let spec = tax_report_spec();
    let taxonomy = TabularTaxonomy {
        taxid: 0,
        taxids: taxids.to_vec(),
    };
    let names = TabularNames::look_up(&spec, &LazyTaxDb::new(path), [&taxonomy], "blastp").unwrap();
    std::fs::remove_dir_all(&cwd).unwrap();
    Some(names)
}

/// The spec of outfmt 18.
fn tax_report_spec() -> TabularSpec {
    let choice = FormatChoice {
        number: 18,
        spec: String::new(),
        delimiter: String::new(),
    };
    TabularSpec::of(ReportFormat::TaxFormat, &choice).unwrap()
}

struct Subjects(Vec<Vec<TaxDefline>>);

impl TaxSubjects for Subjects {
    fn deflines(&self, subject: usize) -> &[TaxDefline] {
        &self.0[subject]
    }
}

fn line(label: &str, title: &str, taxid: i32) -> TaxDefline {
    TaxDefline {
        label: label.as_bytes().to_vec(),
        title: title.as_bytes().to_vec(),
        taxid,
    }
}

fn hsp(subject: usize, evalue: f64, bit_score: f64) -> TaxHsp {
    TaxHsp {
        subject,
        evalue,
        bit_score,
    }
}

/// The organism report (without the caption and the `Query=` block).
fn report(
    subjects: &Subjects,
    hsps: &[TaxHsp],
    names: Option<&TabularNames>,
) -> (Vec<u8>, Vec<u8>) {
    let lookups = TaxLookups::new(names, "BLASTP").unwrap();
    let (mut out, mut err) = (Vec::new(), Vec::new());
    organism_report(&mut err, &mut out, hsps, subjects, &lookups)
        .unwrap()
        .write_to(&mut out)
        .unwrap();
    (out, err)
}

fn text(lines: &[&str]) -> Vec<u8> {
    lines.concat().into_bytes()
}

#[test]
fn add_spaces_centres_cuts_and_pads_as_ncbi() {
    // probe `blastp.d.ids.t2`: 42 + 16 + 42.
    let caption = add_spaces(b"Tax BLAST report", 100, Place::Center, true, true).flat();
    assert_eq!(
        caption,
        format!("\n{0}Tax BLAST report{0}\n", " ".repeat(42)).into_bytes()
    );
    assert_eq!(
        add_spaces(b"abc", 5, Place::End, false, false).flat(),
        b"abc   "
    );
    // Shorter than the value: cut to `width - 3`, then `...` and one space.
    assert_eq!(
        add_spaces(b"abcdefgh", 6, Place::End, false, false).flat(),
        b"abc... "
    );
    // Widths 0 to 2: `width - 3` wraps and `substr` keeps the whole value.
    for width in 0..3 {
        assert_eq!(
            add_spaces(b"Accession", width, Place::End, false, false).flat(),
            b"Accession... "
        );
    }
    assert_eq!(
        add_spaces(b"Description", 0, Place::Center, false, false).flat(),
        b" Description... "
    );
}

#[test]
fn a_query_without_hsps_has_every_width_zero() {
    // probe `blastp.d.ids.t2`, query pq2.
    let (out, err) = report(&Subjects(Vec::new()), &[], None);
    assert_eq!(
        out,
        text(&[
            "\n",
            &" ".repeat(43),
            "Organism Report",
            &" ".repeat(43),
            "\n",
            "Accession...  Description... Score... E-value... \n",
        ])
    );
    assert!(err.is_empty());
}

#[test]
fn rows_group_by_taxid_in_first_seen_order_with_ncbi_widths() {
    // probe `blastp.d.ids.t2` (excerpt): the template tokens of taxid 8888's names and of a
    // title are replaced by the later calls of the chain; a title longer than the description is
    // cut before the replacements, so `<@descr_text@>` survives where the cut keeps it.
    let subjects = Subjects(vec![
        vec![line("XP_000101", "subject one two HSPs", 9606)],
        vec![line("XP_000109", "subject nine", 8888)],
        vec![line("P99902", "subject two", 10090)],
        vec![line(
            "XP_000105",
            "subject five <@score@> s <@evalue@> e <@acc_hd@> a <@descr_text@> d <@taxid@> t",
            32630,
        )],
        vec![line("s7_not_in_taxdb_long_label", "subject seven", 777)],
        vec![line("1ABC_A", "subject six", 4932)],
        vec![line("BAA00110", "subject ten", 12908)],
    ]);
    let names = names(Some("of18"), &[9606, 8888, 10090, 32630, 777, 4932, 12908]);
    let hsps = [
        hsp(0, 5e-109, 313.2),
        hsp(0, 1e-3, 30.0),
        hsp(1, 7e-57, 168.3),
        hsp(2, 4e-53, 159.1),
        hsp(3, 9e-43, 131.7),
        hsp(4, 2e-40, 120.0),
        hsp(5, 1e-45, 139.4),
        hsp(6, 2e-28, 94.04),
    ];
    let (out, err) = report(&subjects, &hsps, names.as_ref());
    let want = text(&[
        "\n",
        &" ".repeat(43),
        "Organism Report",
        &" ".repeat(43),
        "\n",
        "Accession                                          Description                        Score E-value \n",
        "\n                                   Homo sapiens (human) [primates]                                   \n",
        " XP_000101                  subject one two HSPs                                       313   5e-109  \n",
        "\n Tok 8888 (b10090 c 1) sp. (b10090 c <@numhits@... \n",
        " XP_000109                  subject nine                                               168   7e-57   \n",
        "\n                                Mus musculus (house mouse) [rodents]                                \n",
        " P99902                     subject two                                                159   4e-53   \n",
        "\n                                     synthetic construct () [-]                                     \n",
        " XP_000105                  subject five 131    s 9e-43    e Accession                   a <@de... 131   9e-43   \n",
        "\n                               Saccharomyces cerevisiae  [ascomycetes]                               \n",
        " 1ABC_A                     subject six                                                139   1e-45   \n",
        "\n                              unclassified sequences () [unclassified]                              \n",
        " BAA00110                   subject ten                                                94.0  2e-28   \n",
    ]);
    assert_eq!(
        String::from_utf8_lossy(&out),
        String::from_utf8_lossy(&want)
    );
    assert!(err.is_empty());
}

#[test]
fn a_subject_whose_lookup_throws_feeds_every_hsp_and_keeps_earlier_rows() {
    // probe `s2.blastp.md2.good`: definition lines 9606, 777 (not in the database), 10090. The
    // first line's row is kept, the next HSP of the subject is read again (no previous id), and
    // the third line never widens the columns.
    let subjects = Subjects(vec![vec![
        line("XP_000206", "md2 ok first", 9606),
        line("XP_000207", "md2 not in taxdb", 777),
        line("md2_third_long_label_never_seen", "third", 10090),
    ]]);
    let names = names(Some("of18"), &[9606, 777, 10090]);
    let (out, _) = report(
        &subjects,
        &[hsp(0, 1e-86, 243.4), hsp(0, 9.6, 11.9)],
        names.as_ref(),
    );
    let want = text(&[
        "\n",
        &" ".repeat(43),
        "Organism Report",
        &" ".repeat(43),
        "\n",
        "Accession                                 Description                                Score E-value \n",
        "\n                                   Homo sapiens (human) [primates]                                   \n",
        " XP_000206 md2 ok first                                                                243   1e-86   \n",
        " XP_000206 md2 ok first                                                                11.9  9.6     \n",
    ]);
    assert_eq!(
        String::from_utf8_lossy(&out),
        String::from_utf8_lossy(&want)
    );
    // Without the database every lookup fails: no rows, every HSP widens the columns.
    let (out, _) = report(&subjects, &[hsp(0, 1e-86, 243.4), hsp(0, 9.6, 11.9)], None);
    assert!(
        out.ends_with(b"Score E-value \n"),
        "{}",
        String::from_utf8_lossy(&out)
    );
}

#[test]
fn consecutive_hsps_of_a_found_subject_give_one_row() {
    let subjects = Subjects(vec![
        vec![line("s1", "one", 9606)],
        vec![line("s2", "two", 9606)],
    ]);
    let names = names(Some("of18"), &[9606]);
    let hsps = [
        hsp(0, 1e-50, 200.0),
        hsp(0, 123456789.0, 10.0),
        hsp(1, 1e-5, 50.0),
    ];
    let (out, _) = report(&subjects, &hsps, names.as_ref());
    // The second HSP of s1 is skipped, so its 9-byte E-value widens nothing.
    assert!(
        out.windows(b"Score E-value \n".len())
            .any(|w| w == b"Score E-value \n"),
        "{}",
        String::from_utf8_lossy(&out)
    );
    assert_eq!(
        out.split(|&b| b == b'\n')
            .filter(|l| l.starts_with(b" s"))
            .count(),
        2
    );
}

#[test]
fn a_description_width_that_wraps_pads_with_ncbis_runs() {
    // probe `wrap.blastp.s.long_name85.none` (NCBI 2.17.0, 64 GiB): an 85-byte label makes the
    // description's `unsigned int` width 4294967295; the header is `Accession<SP*2147483719>
    // Description<SP*2147483642>Score E-value ` and a row of the `t2` database `<label><SP*1>
    // <title><SP*4294967262>180   7e-63   `.
    let label = "N".repeat(85);
    let subjects = Subjects(vec![vec![line(&label, "s", 0)]]);
    let (out, _) = report(&subjects, &[hsp(0, 7e-63, 180.3)], None);
    let header = add_spaces(b"Accession", 85, Place::End, false, false);
    assert_eq!(header.len(), 9 + 77);
    let mut widths = Widths::default();
    widths.take(label.as_bytes(), b"180", b"7e-63");
    assert_eq!(widths.descr, u32::MAX);
    let descr = add_spaces(
        b"Description",
        u64::from(widths.descr),
        Place::Center,
        false,
        false,
    );
    assert_eq!(
        descr.0,
        vec![
            Seg::Spaces(2147483642),
            Seg::Bytes(b"Description".to_vec()),
            Seg::Spaces(2147483642)
        ]
    );
    let title = add_spaces(
        b"subject title of 33 bytes is here",
        u64::from(widths.descr),
        Place::End,
        false,
        false,
    );
    assert_eq!(title.len(), 33 + 4294967295 - 33 + 1);
    // The report keeps the runs as counts (no 4 GiB in memory) and its bytes are those NCBI
    // writes: no rows without the database.
    let text = {
        let mut lines = Text::default();
        lines.push_bytes(b"\n");
        lines.push_text(&add_spaces(
            b"Organism Report",
            100,
            Place::Center,
            false,
            false,
        ));
        lines.push_bytes(b"\nAccession");
        lines.push_spaces(77 + 2147483642);
        lines.push_bytes(b"Description");
        lines.push_spaces(2147483642);
        lines.push_bytes(b"Score E-value \n");
        lines
    };
    let mut want = Vec::new();
    text.write_to(&mut want).unwrap();
    assert_eq!(out.len(), want.len());
    assert!(out == want);
}

#[test]
fn a_long_organism_line_is_cut_to_97_bytes() {
    // probe `blastp.d.p55555.t2`.
    let names = names(Some("of18"), &[55555]).unwrap();
    let line = organism_line(names.get(55555).unwrap()).flat();
    assert_eq!(line.len(), 1 + 1 + 97 + 3 + 1 + 1);
    assert!(line.starts_with(b"\n Longus nominus xxx"));
    assert!(line.ends_with(b"... \n"));
}

#[test]
fn bl_taxid_in_a_name_is_refused() {
    // NCBI writes `STaxInfo::blNameTaxid`, an uninitialised `int` (`-1784` in probe
    // `s2.blastp.bltaxid`).
    let tsv = "8889\tbl <@bl_taxid@> sp.\tc\tb\tEukaryota\n";
    let dir = empty_cwd();
    std::fs::write(dir.join("t.tsv"), tsv).unwrap();
    let names = {
        let index_data = crate_test_taxdb(&[(8889, b"bl <@bl_taxid@> sp.\tc\tb\tEukaryota")]);
        std::fs::write(dir.join("taxdb.bti"), index_data.0).unwrap();
        std::fs::write(dir.join("taxdb.btd"), index_data.1).unwrap();
        let cwd = empty_cwd();
        let path = SearchPath::new(&cwd, Some(dir.as_os_str()), None).unwrap();
        let taxonomy = TabularTaxonomy {
            taxid: 0,
            taxids: vec![8889],
        };
        TabularNames::look_up(
            &tax_report_spec(),
            &LazyTaxDb::new(path),
            [&taxonomy],
            "blastp",
        )
        .unwrap()
    };
    let error = TaxLookups::new(Some(&names), "BLASTP")
        .err()
        .unwrap()
        .to_string();
    assert!(error.contains("<@bl_taxid@>"), "{error}");
}

#[test]
fn the_tax_browser_url_reads_ncbirc_as_ncbi() {
    use crate::blastinput::ncbi_registry::{read_registry, Layer};
    assert_eq!(
        tax_browser_url(None),
        b"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi".to_vec()
    );
    let read = read_registry(b"[BLASTFMTUTIL]\nPROTOCOL=http:\n", Layer::Application).unwrap();
    assert_eq!(
        tax_browser_url(Some(&read)),
        b"http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi".to_vec()
    );
    // An empty value is no entry (`HasEntry`, `Get(...).empty()`).
    let read = read_registry(
        b"[BLASTFMTUTIL]\nPROTOCOL=\n[TAX_BROWSER]\nBLASTFMTUTIL=\n",
        Layer::Application,
    )
    .unwrap();
    assert_eq!(
        tax_browser_url(Some(&read)),
        b"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi".to_vec()
    );
    let read = read_registry(
        b"[TAX_BROWSER]\nBLASTFMTUTIL=<@protocol@>//x/<@rid@>\n",
        Layer::Application,
    )
    .unwrap();
    assert_eq!(tax_browser_url(Some(&read)), b"https://x/<@rid@>".to_vec());
}

/// A taxonomy database of `rows` (taxid, record), as `make_taxdb.py` writes it.
fn crate_test_taxdb(rows: &[(i32, &[u8])]) -> (Vec<u8>, Vec<u8>) {
    let mut index = Vec::new();
    index.extend(0x8739u32.to_be_bytes());
    index.extend((rows.len() as u32).to_be_bytes());
    index.extend([0u8; 16]);
    let mut data = Vec::new();
    for (taxid, record) in rows {
        index.extend((*taxid as u32).to_be_bytes());
        index.extend((data.len() as u32).to_be_bytes());
        data.extend(*record);
    }
    (index, data)
}
