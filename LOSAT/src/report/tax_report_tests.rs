//! Tests of the organism report (outfmt 18) against NCBI BLAST+ 2.17.0's bytes (probe of line RG
//! T2, `remaining-gaps-20261112/T/t2/out`), with the taxonomy database `tests/seqdb/taxdb/of18`
//! (`of18.tsv`: the probe's `t2` database).

use super::*;
use crate::seqdb::taxdb::LazyTaxDb;
use crate::seqdb::SearchPath;
use std::path::{Path, PathBuf};

fn fixture(variant: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/seqdb/taxdb")
        .join(variant)
}

/// The organism report's names of the taxonomy database in `dir` (as the search path's current
/// directory and `BLASTDB`), opened at the first lookup.
fn names_in(dir: &Path) -> TabularNames {
    let path = SearchPath::new(dir, Some(dir.as_os_str()), None).unwrap();
    TabularNames::for_report(
        TaxNamesDb::new(
            std::sync::Arc::new(LazyTaxDb::new(path)),
            None,
            |error: crate::seqdb::SeqDbError, _: &str| anyhow::anyhow!("{error}"),
        ),
        "blastp",
    )
}

/// The names of the fixture taxonomy database `variant` (`None`: no database).
fn names(variant: Option<&str>) -> Option<TabularNames> {
    Some(names_in(&fixture(variant?)))
}

/// The lookups of `names` without a `./.ncbirc` (the test process's directory is not read).
fn lookups(names: Option<&TabularNames>) -> TaxLookups {
    TaxLookups::with_registry(names, "BLASTP", None)
}

struct Subjects(Vec<Vec<TaxDefline>>);

impl TaxSubjects for Subjects {
    fn deflines(&self, subject: usize) -> &[TaxDefline] {
        &self.0[subject]
    }

    fn refused(&self, subject: usize, reason: &str) -> anyhow::Error {
        anyhow::anyhow!("subject {subject} {reason}")
    }
}

fn line(label: &str, title: &str, taxid: i32) -> TaxDefline {
    TaxDefline {
        label: label.as_bytes().to_vec(),
        title: TaxTitle::Text(title.as_bytes().to_vec()),
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

/// The organism report (without the caption and the `Query=` block) and the taxids it looks up.
fn report(
    subjects: &Subjects,
    hsps: &[TaxHsp],
    names: Option<&TabularNames>,
) -> (Vec<u8>, Vec<i32>) {
    let (taxids, text) = organism_report(hsps, subjects, &mut lookups(names)).unwrap();
    let mut out = Vec::new();
    text.write_to(&mut out).unwrap();
    (out, taxids)
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
    let (out, taxids) = report(&Subjects(Vec::new()), &[], None);
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
    assert!(taxids.is_empty());
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
    let names = names(Some("of18"));
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
    let (out, taxids) = report(&subjects, &hsps, names.as_ref());
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
    // Each taxid once, in the order of the definition lines that reach it.
    assert_eq!(taxids, [9606, 8888, 10090, 32630, 777, 4932, 12908]);
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
    let names = names(Some("of18"));
    let (out, taxids) = report(
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
    // Taxid 777 is looked up again at the next HSP; the third line's taxid never.
    assert_eq!(taxids, [9606, 777, 777]);
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
    let names = names(Some("of18"));
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
    // probe `wrap.blastp.s.long_name85.t2.as64` (NCBI 2.17.0, 64 GiB; `.sq` of its stdout): an
    // 85-byte label makes the description's `unsigned int` width 4294967295; the header is
    // `Accession<SP*2147483719>Description<SP*2147483642>Score E-value ` and the row
    // ` <label> <title><SP*4294967262>180   7e-63   `. The report keeps the runs as counts: the
    // segments are compared, no byte of the runs is made.
    let label = "N".repeat(85);
    let title = "name-only swissprot id of 85 bytes";
    let subjects = Subjects(vec![vec![line(&label, title, 0)]]);
    let names = names(Some("of18"));
    let (taxids, text) = organism_report(
        &[hsp(0, 7e-63, 180.3)],
        &subjects,
        &mut lookups(names.as_ref()),
    )
    .unwrap();
    assert_eq!(taxids, [0]);
    let row = format!(" {label} {title}<SP*4294967262>180   7e-63   \n");
    assert_eq!(
        String::from_utf8(squeezed(&text)).unwrap(),
        [
            &*format!("\n{0}Organism Report{0}\n", " ".repeat(43)),
            "Accession<SP*2147483719>Description<SP*2147483642>Score E-value \n",
            "\n                              Zero taxid sp. (zero common) [zeroblast]                              \n",
            row.as_str(),
        ]
        .concat()
    );
}

/// The text with every run of 64 spaces or more written `<SP*n>`, as the probe's `.sq` files
/// show NCBI's output (`T/t2/scripts/probe.py`).
fn squeezed(text: &Text) -> Vec<u8> {
    let mut out = Vec::new();
    let mut run = 0u64;
    let end_run = |out: &mut Vec<u8>, run: &mut u64| {
        if *run >= 64 {
            out.extend_from_slice(format!("<SP*{run}>").as_bytes());
        } else {
            out.resize(out.len() + *run as usize, b' ');
        }
        *run = 0;
    };
    for seg in &text.0 {
        match seg {
            Seg::Spaces(count) => run += count,
            Seg::Bytes(bytes) => {
                for &byte in bytes {
                    if byte == b' ' {
                        run += 1;
                    } else {
                        end_run(&mut out, &mut run);
                        out.push(byte);
                    }
                }
            }
        }
    }
    end_run(&mut out, &mut run);
    out
}

#[test]
fn a_long_organism_line_is_cut_to_97_bytes() {
    // probe `blastp.d.p55555.t2`.
    let names = names(Some("of18"));
    let mut lookups = lookups(names.as_ref());
    let line = organism_line(lookups.look_up(55555).unwrap().unwrap()).flat();
    assert_eq!(line.len(), 1 + 1 + 97 + 3 + 1 + 1);
    assert!(line.starts_with(b"\n Longus nominus xxx"));
    assert!(line.ends_with(b"... \n"));
}

#[test]
fn bl_taxid_in_a_name_is_refused() {
    // NCBI writes `STaxInfo::blNameTaxid`, an uninitialised `int` (`-1784` in probe
    // `s2.blastp.bltaxid`), where the report shows that taxid's group.
    static COUNTER: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
    let dir = std::env::temp_dir().join(format!(
        "losat-tax-report-{}-{}",
        std::process::id(),
        COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let (index, data) = crate_test_taxdb(&[
        (8889, b"bl <@bl_taxid@> sp.\tc\tb\tEukaryota"),
        (9606, b"Homo sapiens\thuman\tprimates\tEukaryota"),
    ]);
    std::fs::write(dir.join("taxdb.bti"), index).unwrap();
    std::fs::write(dir.join("taxdb.btd"), data).unwrap();
    let names = names_in(&dir);
    // A line of taxid 8889 that the report never reaches (its subject's earlier line throws)
    // is not refused.
    let unreached = Subjects(vec![vec![
        TaxDefline {
            title: TaxTitle::Throws,
            ..line("s0", "", 9606)
        },
        line("s1", "t", 8889),
    ]]);
    let (out, taxids) = report(&unreached, &[hsp(0, 1e-5, 50.0)], Some(&names));
    assert!(taxids.is_empty());
    assert!(out.ends_with(b"Accession...  Description... Score... E-value... \n"));
    let reached = Subjects(vec![vec![line("s1", "t", 8889)]]);
    let error = organism_report(&[hsp(0, 1e-5, 50.0)], &reached, &mut lookups(Some(&names)))
        .err()
        .unwrap()
        .to_string();
    assert!(error.contains("<@bl_taxid@>"), "{error}");
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn a_title_that_throws_drops_its_subject_before_it_widens() {
    // probe `f1.d1.blastp.s.pd.t2` (T/t2/fix1): a `-subject` title of byte 0x81, which NCBI's
    // HtmlDecode cannot read, under a 45-byte label: no row, the accession column stays 9 wide,
    // and every HSP of the subject is read again and dropped.
    let subjects = Subjects(vec![
        vec![TaxDefline {
            title: TaxTitle::Throws,
            ..line("a_long_local_label_of_the_bad_title_record", "", 0)
        }],
        vec![line("s2good", "two", 0)],
    ]);
    let names = names(Some("of18"));
    let (out, taxids) = report(
        &subjects,
        &[
            hsp(0, 0.0, 600.0),
            hsp(0, 1e-10, 60.0),
            hsp(1, 1e-53, 159.1),
        ],
        names.as_ref(),
    );
    let want = text(&[
        "\n",
        &" ".repeat(43),
        "Organism Report",
        &" ".repeat(43),
        "\n",
        "Accession                                 Description                                Score E-value \n",
        "\n                              Zero taxid sp. (zero common) [zeroblast]                              \n",
        " s2good    two                                                                         159   1e-53   \n",
    ]);
    assert_eq!(
        String::from_utf8_lossy(&out),
        String::from_utf8_lossy(&want)
    );
    assert_eq!(taxids, [0]);
}

#[test]
fn a_refused_title_is_refused_where_ncbi_makes_or_shows_it() {
    let refused = |shown_only| TaxDefline {
        title: TaxTitle::Refused {
            reason: "has a title".into(),
            shown_only,
        },
        ..line("s", "", 777)
    };
    let names = names(Some("of18"));
    let run = |subject: TaxDefline| {
        organism_report(
            &[hsp(0, 1e-5, 50.0)],
            &Subjects(vec![vec![subject]]),
            &mut lookups(names.as_ref()),
        )
    };
    // Taxid 777 is not in the database: NCBI drops the row and never shows the title.
    assert!(run(refused(true)).is_ok());
    assert!(run(refused(false)).is_err());
    let shown = TaxDefline {
        taxid: 9606,
        ..refused(true)
    };
    assert_eq!(
        run(shown).err().unwrap().to_string(),
        "subject 0 has a title"
    );
}

#[test]
fn the_tax_browser_url_reads_ncbirc_as_ncbi() {
    use crate::blastinput::ncbi_registry::{read_registry, Layer};
    assert_eq!(
        tax_browser_url(None),
        b"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi".to_vec()
    );
    let read = read_registry(b"[BLASTFMTUTIL]\nPROTOCOL=http:\n", Layer::SystemWide).unwrap();
    assert_eq!(
        tax_browser_url(Some(&read)),
        b"http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi".to_vec()
    );
    // An empty value is no entry (`HasEntry`, `Get(...).empty()`).
    let read = read_registry(
        b"[BLASTFMTUTIL]\nPROTOCOL=\n[TAX_BROWSER]\nBLASTFMTUTIL=\n",
        Layer::SystemWide,
    )
    .unwrap();
    assert_eq!(
        tax_browser_url(Some(&read)),
        b"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi".to_vec()
    );
    let read = read_registry(
        b"[TAX_BROWSER]\nBLASTFMTUTIL=<@protocol@>//x/<@rid@>\n",
        Layer::SystemWide,
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

#[test]
fn a_registry_read_posts_as_ncbi() {
    // NCBI BLAST+ 2.17.0 (probe `T/t2/probe_rc`, `rep` and `bom`): the posts of one read; each
    // query's report reads `./.ncbirc` twice.
    use crate::blastinput::ncbi_registry::{read_registry, Layer};
    let read = read_registry(b"[FOO]\na=1\na=2\n", Layer::SystemWide).unwrap();
    assert_eq!(
        registry_posts(&read, "blastp"),
        b"Warning: (110.8) [blastp] Found multiple [FOO] a settings; using the one from line 3\n"
            .to_vec()
    );
    let read = read_registry(b"\xef", Layer::SystemWide).unwrap();
    assert_eq!(
        registry_posts(&read, "blastp"),
        b"Error: (110.4) [blastp] Error reading the registry after line 1: \n".to_vec()
    );
}
