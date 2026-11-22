//! The format guess against NCBI's: each expected format is what NCBI BLAST+ 2.17.0's `CFormatGuess().Format(in)`
//! answered for the same bytes (a harness of `libxutil` 2.17.0 built from the release's source, comparison only;
//! line RG D3-1 `F/d3a/fg/unit_harness.tsv`, cases `F/d3a/fg/unit_cases.py`), and for the `-in_pssm` cases the
//! official psiblast and tblastn (`F/d3a/fg/out/ncbi`).

use super::*;

const TEXT: &[u8] = b"PssmWithParameters ::= {\n  pssm {\n    isProtein TRUE\n  }\n}\n";

fn guess(bytes: &[u8]) -> Result<Guessed, Stop> {
    FormatGuess::over(bytes).guess_format()
}

fn format_of(bytes: &[u8]) -> Format {
    guess(bytes).expect("a decided guess").format
}

#[test]
fn every_test_answers_as_ncbis_format_guess() {
    let cases: &[(&str, &[u8], Format)] = &[
        ("zip56", b"PK\x05\x06\x00\x00\x00\x00", Format::Zip),
        ("zip_bad", b"PK\x05\x07xx", Format::BinaryAsn),
        ("gzip", b"\x1f\x8b\x08\x00", Format::GZip),
        ("zstd", b"(\xb5/\xfd\x00", Format::Zstd),
        ("bzip2", b"BZh91AY", Format::BZip2),
        ("bzh_text", b"BZhello ::= {\n}\n", Format::TextAsn),
        ("lzo3", b"LZO", Format::Lzo),
        ("lzo_off", b"\x89LZO\x00\x0d", Format::Lzo),
        ("lzo_text", b"LZOx ::= {\n}\n", Format::UnknownOrJson),
        ("sra", b"NCBI.sra\x88\x19\x03\x05\x00\x00\x00\x00", Format::Sra),
        ("psl", b"1\t1\t1\t1\t1\t1\t1\t1\t+\tq\t1\t1\t1\tt\t1\t1\t1\t1\t1,\t1,\t1,\n", Format::Psl),
        ("psl_col", b"x\t1\t1\t1\t1\t1\t1\t1\t1\t-+\tq\t1\t1\t1\tt\t1\t1\t1\t1\t1,\t1,\t1,\n", Format::Psl),
        ("rmo_header", b"   SW perc query position matching\x02\nscore div. del. ins. sequence\nx\n", Format::Rmo),
        ("rmo_blank", b"  \n \t \n", Format::Rmo),
        ("rmo_line", b"1 2.5 0.0 1.0 q 10 20 (5) + r c 1 2 3\n", Format::Rmo),
        ("vcf", b"PssmWithParameters ::= {\n##fileformat=VCFv4.2\n  pssm {\n    isProtein TRUE\n  }\n}\n", Format::Vcf),
        ("gvf_version", b"##gvf-version 1\n\x02\n", Format::Gvf),
        ("gff3_version", b"##gff-version 3\n\xef\xbb\xbfPssmWithParameters ::= {\n  pssm {\n    isProtein TRUE\n  }\n}\n", Format::Gff3),
        ("gff3_line", b"c\ts\tf\t1\t2\t.\t+\t.\tID=x\n", Format::Gff3),
        ("gtf", b"c\ts\tf\t1\t2\t.\t+\t.\tgene_id x\x02\n", Format::Gtf),
        ("augustus", b"c\tA\tgene\t1\t2\t1\t+\t.\tg1\x02\n", Format::GffAugustus),
        ("gff2", b"c\ts\tf\t1\t2\t.\t+\t.\n", Format::Gff2),
        ("glimmer3", b">x\x02\nA 1 2 3 1.0\n", Format::Glimmer3),
        ("agp_blank", b"#a\x02\n  \n", Format::Agp),
        ("agp_line", b"c 1 10 1 W x 1 10 +\n", Format::Agp),
        ("xml", b"\n <?XmL version=\"1.0\"?>\n", Format::Xml),
        ("xml_blast4", b"<Blast4-request>", Format::Xml),
        ("wiggle", b"\xef\xbb\xbf\x02xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\nvariableStep\n", Format::Wiggle),
        ("wiggle_track", b"track type=wiggle_0\x02\n", Format::Wiggle),
        ("newick", b"(a,b)\x02;\n", Format::Newick),
        ("newick_unbalanced", b"(a,b)),c\x02\n", Format::BinaryAsn),
        ("nexus_trees", b"PssmWithParameters ::= {\n  pssm {\n  -- #NEXUS\n  -- begin trees;\n    isProtein TRUE\n  }\n}\n", Format::Newick),
        ("nexus_text", b"PssmWithParameters ::= {\n  pssm {\n    isProtein TRUE\n  }\n}\n-- #NEXUS\n", Format::TextAsn),
        ("fallback_fivecol", b"\x00\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff", Format::FiveColFeatureTable),
        ("fallback_bom_short", b"\xef\xbb\xbfPssm", Format::BinaryAsn),
        ("bed", b"chr1 1 2 x\x02\n", Format::Bed),
        ("bed_track", b"track\x02\n", Format::Bed),
        ("bed15", b"c\t1\t2\tn\t1\t+\t1\t2\t0\t1\t1\t0\t1\t1\tx\x02\n", Format::Bed15),
        ("hgvs", b"NC_1:g.1A>T\x02\n", Format::Hgvs),
        ("distmatrix", b"3\nA\x02\nB 1\n", Format::DistanceMatrix),
        ("flatseq", b"1 acgtacgtac acgtacgtac acgtacgtac acgtacgtac acgtacgtac \x02\n", Format::FlatFileSequence),
        ("uniprot", b"ID ::= x\nAC x\nDT x\nDT x\nDT x\nDE x\nOS x\nXX x\nXX x\nXX x\nXX x\nXX x\nXX x\nXX x\nXX x\nXX x\nXX x\n", Format::FlatFileUniProt),
        ("uniprot_short", b"ID ::= x\nAC x\nDT x\nDT x\nDT x\nDE x\nOS x\n", Format::TextAsn),
        ("ena", b"ID x\x02\nAC x\nDT x\nDT x\nDE x\nKW x\nOS x\nOC x\nXX\nXX\nXX\nXX\nXX\nXX\nXX\nXX\nXX\nXX\nXX\nXX\n", Format::FlatFileEna),
        ("genbank", b"LOCUS x\x02\nDEFINITION x\nACCESSION x\nVERSION x\nKEYWORDS x\nx\nx\nx\nx\nx\n", Format::FlatFileGenbank),
        ("genbank_cont", b"LOCUS x\x02\nDEFINITION x\n            more\nACCESSION x\nNID y\nVERSION x\nPROJECT p\nKEYWORDS x\nx\nx\nx\n", Format::FlatFileGenbank),
        ("fivecol", b">Feature x\x02\n", Format::FiveColFeatureTable),
        ("snp", b"\xef\xbb\xbf\nrs1 2 3\n", Format::Bed),
        ("snp_signs", b"rs+1-2\t+3\x02\n", Format::SnpMarkers),
        ("fasta", b">q\nMKV\n", Format::Fasta),
        ("fasta_bar", b">a|b\x02\n", Format::BinaryAsn),
        ("text", b"PssmWithParameters ::= {\n  pssm {\n    isProtein TRUE\n  }\n}\n", Format::TextAsn),
        ("clustal", b"CLUSTAL W\x02\n\nseq1 ACGT\nseq2 ACGT\n\n", Format::Alignment),
        ("alignment_nexus", b"\xef\xbb\xbfPssmWithParameters ::= {\n  pssm {\n    isProtein TRUE\n  }\n}\n-- #NEXUS\n", Format::Alignment),
        ("table", b"a b\x01\nc d\nd e\n", Format::UnknownOrJson),
        ("table_skip", b"\x02\x02\n\x02\nx y\nx y\nx y\nx y\nx y\n", Format::BinaryAsn),
        ("binary", b"\xef\xbb\xbfPssmWithParameters ::= {\n  pssm {\n    isProtein TRUE\n  }\n}\n", Format::BinaryAsn),
        ("one", b"\x01PssmWithParameters ::= {\n  pssm {\n    isProtein TRUE\n  }\n}\n", Format::UnknownOrJson),
        ("phrap", b"DNA x\x02\nACGTACGTACGT\n", Format::BinaryAsn),
        ("phrap_as", b"AS 1 2\x02\nacgtacgtac\n", Format::BinaryAsn),
        ("unknown", b"hello world\n", Format::UnknownOrJson),
        ("empty", b"", Format::UnknownOrJson),
        ("all_comment", b"-- only a comment\n", Format::UnknownOrJson),
    ];
    for &(name, bytes, want) in cases {
        let guessed = guess(bytes).expect(name);
        assert_eq!(guessed.format, want, "{name}");
        assert!(!guessed.fasta_undecided, "{name}");
        // A reader read piece by piece guesses as the bytes in memory do.
        let read = FormatGuess::new(Pieces(bytes)).guess_format().expect(name);
        assert_eq!(read, guessed, "{name}");
    }
}

/// A reader that gives at most 7 bytes a read.
struct Pieces<'a>(&'a [u8]);

impl std::io::Read for Pieces<'_> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.0.len().min(buf.len()).min(7);
        buf[..n].copy_from_slice(&self.0[..n]);
        self.0 = &self.0[n..];
        Ok(n)
    }
}

/// Bytes of every kind the tests look at, from a fixed seed.
fn mixed_bytes(seed: u64, len: usize) -> Vec<u8> {
    const ALPHABET: &[u8] =
        b"ACGTacgtNXZ>{}-#;!: \t\r\n\x0b\x00\x01\x7f\x80\xef\xbb\xbf0123456789|.,+";
    let mut state = seed;
    (0..len)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ALPHABET[(state >> 33) as usize % ALPHABET.len()]
        })
        .collect()
}

/// `EnsureStats` as NCBI writes it: line by line (`\r`, `\n`, either ends one), `\n` appended to each non-empty
/// line, each byte counted by its kind.
fn stats_line_by_line(data: &[u8]) -> Stats {
    let mut stats = Stats::default();
    for line in data
        .split(|&c| c == b'\r' || c == b'\n')
        .filter(|line| !line.is_empty())
    {
        let is_header = line[0] == b'>';
        for &c in line.iter().chain(std::iter::once(&b'\n')) {
            let type_ = SYMBOL_TYPE_TABLE[c as usize];
            if type_ & (F_ALPHA | F_DIGIT | F_SPACE) != 0 {
                stats.al_num_chars += 1;
            } else if c == b'{' || c == b'}' {
                stats.braces += 1;
            }
            if !is_header {
                stats.data += usize::from(type_ & F_SPACE == 0);
                stats.dna_chars += usize::from(type_ & F_DNA_MAIN_ALPHABET != 0);
                stats.aa_chars += usize::from(type_ & F_PROTEIN_ALPHABET != 0);
            }
        }
    }
    stats
}

#[test]
fn the_counts_and_splits_are_ncbis() {
    for seed in 0..200 {
        let data = mixed_bytes(seed, (seed as usize * 37) % 400);
        let (want, got) = (stats_line_by_line(&data), stats_of(&data));
        assert_eq!(
            (
                want.al_num_chars,
                want.braces,
                want.data,
                want.dna_chars,
                want.aa_chars
            ),
            (
                got.al_num_chars,
                got.braces,
                got.data,
                got.dna_chars,
                got.aa_chars
            ),
            "{data:?}"
        );
        // `NStr::Split(..., fSplit_Tokenize)`: the non-empty runs between delimiters.
        let runs: Vec<&[u8]> = data
            .split(|c| b" \t".contains(c))
            .filter(|run| !run.is_empty())
            .collect();
        let toks = tokens(&data, b" \t");
        assert_eq!(toks.len(), runs.len());
        for (i, run) in runs.iter().enumerate().take(KEPT_TOKENS) {
            assert_eq!(toks.at(i), *run);
        }
        // The lines of each split: the non-empty runs between line ends.
        for ends in [
            (|c: u8| c == b'\n') as fn(u8) -> bool,
            |c| c == b'\r',
            |c| c == b'\r' || c == b'\n',
        ] {
            let want = data
                .split(|&c| ends(c))
                .filter(|run| !run.is_empty())
                .count();
            assert_eq!(count_runs(&data, ends), want, "{data:?}");
        }
        // `NStr::Find`: the first occurrence.
        for pattern in [
            &b"#"[..],
            b"\r\n",
            b"AC",
            b"|.",
            b"\xef\xbb\xbf",
            b"ACGTACGT",
        ] {
            let want = data
                .windows(pattern.len())
                .position(|window| window == pattern);
            assert_eq!(find(&data, pattern), want, "{pattern:?} in {data:?}");
        }
    }
}

#[test]
fn the_test_buffer_doubles_while_it_is_all_comment() {
    let concat = |parts: &[&[u8]]| parts.concat();
    // 9,600 bytes of ASN.1 comments: the first 8,096 are all comment (the last line, cut, is dropped), so the
    // buffer doubles and holds the text (NCBI: text ASN.1; psiblast searches it, fg `t01_comment9k_pssm`).
    let comments = b"-- a comment line of ASN.1 text\n".repeat(300);
    assert_eq!(format_of(&concat(&[&comments, TEXT])), Format::TextAsn);
    // A line of 8,200 spaces: no line end in the first buffer, which `EnsureSplitLines` takes for one cut line
    // and `IsAllComment` for no line (fg `t04_spaces8k_pssm`).
    assert_eq!(
        format_of(&concat(&[&b" ".repeat(8200), b"\n", TEXT])),
        Format::TextAsn
    );
    // No line end at all: the doubled buffer holds the whole file as one line.
    assert_eq!(
        format_of(&concat(&[&b"a".repeat(9000), b" ::= {"])),
        Format::TextAsn
    );
    // `#` lines are comments to `IsAllComment` but not to `TestFormatTextAsn` (fg `t02_hash9k_pssm`).
    let hashes = b"# a hash comment line\n".repeat(450);
    assert_eq!(format_of(&concat(&[&hashes, TEXT])), Format::UnknownOrJson);
    // Comments that end before the buffer is full: no buffer, eUnknown (`TestFormatNewick` does not take it).
    assert_eq!(format_of(b"-- only a comment\n"), Format::UnknownOrJson);
    assert_eq!(format_of(&comments), Format::UnknownOrJson);
}

#[test]
fn the_nexus_scan_reads_the_stream_on() {
    let nexus = [TEXT, b"-- #NEXUS\n"].concat();
    // No `begin trees;`: TestFormatNewick reads to the end, eTextASN reads from there.
    assert_eq!(
        guess(&nexus).unwrap(),
        Guessed {
            format: Format::TextAsn,
            stream_pos: nexus.len(),
            fasta_undecided: false,
        }
    );
    // `BEGIN TREES;` in the second chunk of 16,384 bytes (case does not matter).
    let late = [&nexus[..], &b"-".repeat(20000), b"\nBEGIN TREES;\n"].concat();
    assert_eq!(format_of(&late), Format::Newick);
    // A NUL in the second chunk hides the rest of it from the search: eAlignment (the `#NEXUS` line).
    let hidden = [
        &nexus[..],
        &b"-".repeat(16370),
        b"\x00",
        &b"x".repeat(30),
        b"begin trees;\n",
    ]
    .concat();
    assert_eq!(format_of(&hidden), Format::Alignment);
    // A text without NEXUS is read from its start.
    assert_eq!(guess(TEXT).unwrap().stream_pos, 0);
}

#[test]
fn tests_that_read_past_their_data_or_throw() {
    // A GVF line of eight columns (`tokens[8]`) and an `AS` line of one token (`values[1]`): NCBI's harness crashed
    // on both (SIGSEGV).
    assert!(matches!(
        guess(b"c s snv 1 2 . + .\n"),
        Err(Stop::Undecided(_))
    ));
    assert!(matches!(guess(b"AS\n"), Err(Stop::Undecided(_))));
    assert!(matches!(guess(b"AS 1\n"), Err(Stop::Undecided(_))));
    // `AS x`: the second token is not a number, `&&` stops before the third.
    assert_eq!(format_of(b"AS x\n"), Format::UnknownOrJson);
    // A Glimmer3 frame beyond `int`: `NStr::StringToInt`'s exception, at its throw in ncbistr.cpp.
    let Err(Stop::Exception(exception)) = guess(b">x\nA 1 2 99999999999 1.0\n") else {
        panic!("an exception");
    };
    // NCBI's report in LOSAT's words (`crate::wording::losat_text`).
    assert_eq!(exception.site().line, 640);
    assert!(String::from_utf8(exception.what()).unwrap().contains(&*crate::wording::losat_text(
        "NCBI C++ Exception:\n    T0 \"/src/corelib/ncbistr.cpp\", line 640: Error: (CStringException::eConvert) ncbi::NStr::StringToInt() - Cannot convert string '99999999999' to int, overflow (m_Pos = 0)\n"
    )));
    // TestFormatFasta's search for `|` past four bytes: undecided, and the guess without eFasta goes on.
    let guessed = guess(b">ab\x02").unwrap();
    assert!(guessed.fasta_undecided);
    assert_eq!(guessed.format, Format::BinaryAsn);
    // Eleven bytes or more, or a NUL before the end, decide it.
    assert!(!guess(b">abcdefghij\x02").unwrap().fasta_undecided);
    assert!(!guess(b">ab\x00\x02").unwrap().fasta_undecided);
}

#[test]
fn split_in_two_keeps_ncbis_empty_parts() {
    assert_eq!(
        split_in_two(b"ID   x", b" ", true),
        (&b"ID"[..], &b"x"[..], true)
    );
    assert_eq!(
        split_in_two(b"ID   x", b" ", false),
        (&b"ID"[..], &b"  x"[..], true)
    );
    assert_eq!(
        split_in_two(b" ID", b" ", true),
        (&b""[..], &b"ID"[..], true)
    );
    assert_eq!(
        split_in_two(b"ID", b" ", true),
        (&b"ID"[..], &b""[..], false)
    );
    assert_eq!(
        split_in_two(b"", b" \t", false),
        (&b""[..], &b""[..], false)
    );
}
