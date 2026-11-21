use super::appears_fasta;
use std::io::{self, Read};

// Every expected value below is the answer of makeblastdb 2.17.0 (`-dbtype prot` and `nucl` agree): it printed
// "does NOT appear to be FASTA" for `false` and not for `true`. `name` is the case of
// db2/fg/cases/<name>.fa that holds the same bytes (made by scripts/db2_fg_gen.py).

fn check(name: &str, input: &[u8], want: bool) {
    assert_eq!(
        appears_fasta(input),
        want,
        "case {name} ({} bytes)",
        input.len()
    );
}

#[test]
fn small_inputs() {
    let cases: &[(&str, &[u8], bool)] = &[
        ("two_records", b">a\nACGT\n>b\nTTTT\n", true),
        ("protein", b">p\nMKVLAAGIVGLLLAQPAMA\nMKV\n", true),
        ("no_final_newline", b">s\nACGTACGTACGT", true),
        ("crlf", b">s\r\nACGTACGT\r\nACGTACGT\r\n", true),
        ("cr_only", b">s\rACGTACGT\rACGTACGT\r", true),
        ("lfcr", b">s\n\rACGTACGT\n\rACGTACGT\n\r", true),
        ("mixed_ends", b">s\rACGTACGT\r\nACGTACGT\nACGTACGT\n", true),
        ("seq_without_header", b"ACGTACGTACGT\nACGTACGT\n", false),
        ("plain_text", b"hello world, this is not a sequence file\nsecond line\n", false),
        ("fastq", b"@r1\nACGT\n+\nIIII\n", false),
        ("gff", b"##gff-version 3\nchr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1\n", false),
        ("text_asn", b"Seq-entry ::= seq {\n  id { local str \"a\" } ,\n  inst { repr raw, mol na, length 4 } }\n", false),
        ("leading_space_header_inline", b" >s\nACGTACGT\n", false),
        ("leading_space_header", b"  \n\t>s\nACGTACGT\n", false),
        ("hash_comment_then_fasta", b"#ACGTACGTACGTACGTACGT\n>s\nACGTACGT\nACGTACGT\n", true),
        ("semicolon_comment_then_fasta", b";ACGTACGTACGTACGTACGT\n>s\nACGTACGT\nACGTACGT\n", true),
        ("bang_comment_then_fasta", b"!ACGTACGTACGTACGTACGT\n>s\nACGTACGT\nACGTACGT\n", true),
        ("comment_two_then_fasta", b"#ACGTACGTACGT\n;ACGTACGTACGT\n   >s\nACGT\n", false),
        ("dash_comment_then_fasta", b"--ACGTACGTACGTACGTACGT\n>s\nACGTACGT\n", false),
        ("comment_nonl_only", b"#ACGTACGT", false),
        ("nul_in_seq", b">s\nACGT\x00ACGT\nACGTACGT\n", true),
        ("nul_first", b"\x00>s\nACGT\nACGT\n", false),
        ("bom_fasta", b"\xef\xbb\xbf>s\nACGT\nACGT\n", false),
        ("latin1_header_fasta", b">caf\xe9 s\nACGTACGTACGT\nACGTACGT\n", true),
        ("fasta_header_gt_in_seq", b">s\nAC>GT\nACGT\n", false),
        ("only_gt", b">", false),
        ("gt_bar", b">|", false),
        ("tiny_abc", b">abc", false),
        ("lonely_nl", b"\n", false),
        ("nul_only_file", &[0u8; 100], false),
        // defline-only files: the '|' must be at index <= 10 of the sample
        ("defline_bar_at_9", b">aaaaaaaa|bbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\n", true),
        ("defline_bar_at_10", b">aaaaaaaaa|bbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\n", true),
        ("defline_bar_at_11", b">aaaaaaaaaa|bbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\n", false),
        ("defline_bar_at_10_crlf", b">aaaaaaaaa|bbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\r\n", true),
        ("defline_bar_at_11_nonl", b">aaaaaaaaaa|bbbbbbbbbbbbbbbbbbbbbbbbbbbbbb", false),
        ("defline_no_bar", b">aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n", false),
        ("defline_ten_chars_bar", b">sp|P12345|FOO_HUMAN something\n", true),
        ("defline_bar_second_line_early", b">aaaaaaaaaaaaaaa\n>b|bb\n", false),
        ("nul_defline_only_before_bar", b">ab\x00cd|ef\n", false),
        ("empty_seq_spaces", b">s\n \n", false),
    ];
    assert!(cases.len() >= 30);
    for &(name, input, want) in cases {
        check(name, input, want);
    }
}

/// There is no test buffer for an empty stream, so `GuessFormat` cannot answer `eFasta`. (`makeblastdb` stops with
/// "is empty" before it asks, so there is no oracle answer for this.)
#[test]
fn empty_input_is_not_fasta() {
    assert!(!appears_fasta(&b""[..]));
}

// The generators below build the same bytes as scripts/db2_fg_gen.py.

fn repeat_cut(unit: &[u8], n: usize) -> Vec<u8> {
    let mut v = unit.repeat(n / unit.len() + 1);
    v.truncate(n);
    v
}

/// `n` bytes of `#AAAA...` comment lines (DNA-looking, so a later FASTA part still passes the 0.91 test).
fn comment_lines(n: usize, sep: &[u8]) -> Vec<u8> {
    let mut line = b"#".to_vec();
    line.extend_from_slice(&[b'A'; 79]);
    line.extend_from_slice(sep);
    repeat_cut(&line, n)
}

fn dna_lines(n: usize, sep: &[u8]) -> Vec<u8> {
    let mut line = b"ACGT".repeat(15);
    line.extend_from_slice(sep);
    repeat_cut(&line, n)
}

fn cat(parts: &[&[u8]]) -> Vec<u8> {
    parts.concat()
}

fn fasta_of_size(n: usize) -> Vec<u8> {
    cat(&[b">s\n", &dna_lines(n - 3, b"\n")])
}

/// `comment_lines(n)` with `k` evenly spread body bytes (not the '#', not the line end) replaced by `byte`.
fn comment_marked(n: usize, k: usize, byte: u8) -> Vec<u8> {
    let mut b = comment_lines(n, b"\n");
    let pos: Vec<usize> = (0..b.len())
        .filter(|&i| b[i] != b'\n' && b[i] != b'#')
        .collect();
    for j in 0..k {
        b[pos[j * pos.len() / k]] = byte;
    }
    b
}

/// FASTA of `n_total` bytes where `count` evenly spread sequence bytes are replaced by `junk`.
fn fasta_with(n_total: usize, junk: u8, count: usize) -> Vec<u8> {
    let mut body = cat(&[b">s\n", &dna_lines(n_total - 3, b"\n")]);
    let pos: Vec<usize> = (3..body.len()).filter(|&i| body[i] != b'\n').collect();
    for j in 0..count {
        body[pos[j * pos.len() / count]] = junk;
    }
    body
}

/// A header and 1000 data characters of which `good` are `letter` and the others the digit `1`, 50 to a line.
fn frac_case(good: usize, letter: u8) -> Vec<u8> {
    let mut chars = vec![letter; good];
    chars.extend_from_slice(&vec![b'1'; 1000 - good]);
    let lines: Vec<&[u8]> = chars.chunks(50).collect();
    cat(&[b">s\n", &lines.join(&b'\n'), b"\n"])
}

/// A header of `%` and `>` and `n - header` bytes of DNA: `alnum` alphanumeric-or-space characters in `n` bytes.
fn alnum_case(alnum: usize) -> Vec<u8> {
    let n = 1000;
    let mut header = b">".to_vec();
    header.extend_from_slice(&vec![b'%'; n - alnum - 1]);
    header.push(b'\n');
    cat(&[&header, &dna_lines(n - header.len(), b"\n")])
}

/// A defline-only file of 400 bytes with `alnum` alphanumeric-or-space characters (the '|' is at index 1).
fn defline_alnum_case(alnum: usize) -> Vec<u8> {
    let mut v = b">|".to_vec();
    v.extend_from_slice(&vec![b'%'; 400 - alnum - 2]);
    v.extend_from_slice(&vec![b'x'; alnum - 1]);
    v.push(b'\n');
    v
}

#[test]
fn sample_size_rule() {
    // the sample is the first 8096 bytes (more only for all-comment samples)
    for n in [8095, 8096, 8097, 8098, 16191, 16192, 16193] {
        check(&format!("size_{n}"), &fasta_of_size(n), true);
    }
    let junk = b"%%%%%%%%%%%%\n".repeat(5000);
    check(
        "head_fasta_tail_junk",
        &cat(&[&fasta_of_size(8096), &junk]),
        true,
    );
    check(
        "head_fasta_tail_junk_8000",
        &cat(&[&fasta_of_size(8000), &junk]),
        true,
    );
    let junk_head = b"%%%%%%%%%%%%\n".repeat(700);
    check(
        "head_junk_tail_fasta",
        &cat(&[&junk_head, &fasta_of_size(100000)]),
        false,
    );
    let binary: Vec<u8> = (0..=255u8).collect::<Vec<_>>().repeat(200);
    check(
        "head_fasta_tail_binary",
        &cat(&[&fasta_of_size(8096), &binary]),
        true,
    );
    let binary_short: Vec<u8> = (1..=255u8).collect::<Vec<_>>().repeat(3);
    check(
        "head_fasta_tail_binary_short",
        &cat(&[&fasta_of_size(4000), &binary_short]),
        false,
    );
}

#[test]
fn all_comment_sample_doubles() {
    // only comments: false at every size (the buffer doubles until the data ends before the buffer does)
    for n in [
        8095, 8096, 8097, 9000, 16191, 16192, 16193, 17000, 33000, 70000,
    ] {
        check(&format!("allcomment_{n}"), &comment_lines(n, b"\n"), false);
    }
    check(
        "allcomment_crlf_20000",
        &comment_lines(20000, b"\r\n"),
        false,
    );
    check("allcomment_cr_20000", &comment_lines(20000, b"\r"), false);
    let dashes = (b"--".iter().chain(&[b'A'; 77]).chain(b"\n"))
        .copied()
        .collect::<Vec<u8>>();
    check("allcomment_dash_20000", &dashes.repeat(250), false);
    // comments, then FASTA: found as soon as the doubled buffer reaches it
    let tail = cat(&[b"\n>s\n", &dna_lines(30000, b"\n")]);
    for (c, want) in [
        (8000, true),
        (8100, true),
        (9000, true),
        (16300, true),
        (17000, true),
        (40000, true),
        (4145152, true),
        (8000000, true),
        (8290300, true),
        (8290400, false),
        (8400000, false),
    ] {
        check(
            &format!("comment_{c}_then_fasta"),
            &cat(&[&comment_lines(c, b"\n"), &tail]),
            want,
        );
    }
    // the sample ends inside the comments: the cap of 1024 * 8096 bytes
    check("comment_4145152", &comment_lines(4145152, b"\n"), false);
    check("comment_8300000", &comment_lines(8300000, b"\n"), false);
    // CR and CR LF comment lines
    let cr_tail = cat(&[b">s\r", &dna_lines(30000, b"\r")]);
    check(
        "comment_9000_then_fasta_cr",
        &cat(&[&comment_lines(9000, b"\r"), &cr_tail]),
        false,
    );
    let crlf_tail = cat(&[b">s\r\n", &dna_lines(30000, b"\r\n")]);
    check(
        "comment_9000_then_fasta_crlf",
        &cat(&[&comment_lines(9000, b"\r\n"), &crlf_tail]),
        false,
    );
    check(
        "comment_cr_then_fasta",
        b"#ACGTACGTACGT\r>s\rACGT\rACGT\r",
        true,
    );
}

#[test]
fn printable_and_high_byte_limits() {
    let tail = cat(&[b"\n>s\n", &dna_lines(100000, b"\n")]);
    // at most 10% of the sample may be non-printable (the line ends count) for the comment test to apply
    for (c, want) in [
        (600, true),
        (709, true),
        (710, true),
        (711, false),
        (800, false),
    ] {
        check(
            &format!("ctrl_comment8096_{c}_then_fasta"),
            &cat(&[&comment_marked(8096, c, 1), &tail]),
            want,
        );
    }
    for (c, want) in [(1000, true), (1421, true), (1422, false), (1700, false)] {
        let second = comment_marked(8096, c, 1);
        check(
            &format!("ctrl_comment_2nd_{c}_then_fasta"),
            &cat(&[&comment_lines(8096, b"\n"), &second, &tail]),
            want,
        );
    }
    // more than 1 byte in 20 >= 0x80 makes EnsureSplitLines give up, and the sample counts as all comment
    for (h, want) in [
        (0, true),
        (200, true),
        (404, true),
        (405, true),
        (450, true),
        (451, false),
        (500, false),
    ] {
        check(
            &format!("hi_fasta_9000_{h}"),
            &fasta_with(9000, 0xE9, h),
            want,
        );
    }
    for (h, want) in [(350, true), (404, true), (405, false), (810, false)] {
        check(
            &format!("hi_fasta_8096_{h}"),
            &fasta_with(8096, 0xE9, h),
            want,
        );
    }
    for (h, want) in [(100, true), (404, false), (500, false)] {
        check(
            &format!("hi_fasta_5000_{h}"),
            &fasta_with(5000, 0xE9, h),
            want,
        );
    }
    check(
        "high_in_seq_few",
        &cat(&[
            b">s\n",
            &b"ACGT".repeat(100),
            b"\xff",
            &b"ACGT".repeat(100),
            b"\n",
        ]),
        true,
    );
}

#[test]
fn fraction_limits() {
    // alphanumeric-or-space characters / sample size must be >= 0.8
    for (a, want) in [
        (780, false),
        (798, false),
        (799, true),
        (800, true),
        (810, true),
    ] {
        check(&format!("alnum_{a}_of_1000"), &alnum_case(a), want);
    }
    // defline-only files: the same ratio must be >= 0.75
    for (a, want) in [(290, false), (299, false), (300, true), (301, true)] {
        check(
            &format!("defline_alnum_{a}_of_400"),
            &defline_alnum_case(a),
            want,
        );
    }
    // more than 91% of the data characters of DNA (ACGNTU-) or of protein letters
    for (g, want) in [
        (880, false),
        (909, false),
        (910, false),
        (911, true),
        (920, true),
    ] {
        check(&format!("dna_{g}_of_1000"), &frac_case(g, b'G'), want);
        check(&format!("aa_{g}_of_1000"), &frac_case(g, b'L'), want);
    }
    // J and O are not protein letters; '*' and '-' are not alphanumeric, so the 0.8 test fails first
    for (letter, name) in [(b'J', "junkletter"), (b'O', "o"), (b'*', "star")] {
        check(
            &format!("{name}_911_of_1000"),
            &frac_case(911, letter),
            false,
        );
    }
    check("n_911_of_1000", &frac_case(911, b'N'), true);
    check("x_911_of_1000", &frac_case(911, b'X'), true);
    check("dash_911_of_1000", &frac_case(911, b'-'), false);
}

/// A reader that gives its data in pieces of `chunk` bytes, then fails.
struct Chunked<'a> {
    data: &'a [u8],
    chunk: usize,
    fail_after_data: bool,
    interrupt_next: bool,
}

impl Read for Chunked<'_> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.interrupt_next {
            self.interrupt_next = false;
            return Err(io::ErrorKind::Interrupted.into());
        }
        self.interrupt_next = true;
        if self.data.is_empty() {
            return if self.fail_after_data {
                Err(io::Error::other("broken"))
            } else {
                Ok(0)
            };
        }
        let n = self.data.len().min(self.chunk).min(buf.len());
        buf[..n].copy_from_slice(&self.data[..n]);
        self.data = &self.data[n..];
        Ok(n)
    }
}

#[test]
fn reader_behaviour_does_not_matter() {
    let data = cat(&[
        &comment_lines(9000, b"\n"),
        b"\n>s\n",
        &dna_lines(3000, b"\n"),
    ]);
    for chunk in [1, 7, 4096, 100000] {
        for fail in [false, true] {
            let reader = Chunked {
                data: &data,
                chunk,
                fail_after_data: fail,
                interrupt_next: true,
            };
            assert!(appears_fasta(reader), "chunk {chunk} fail {fail}");
        }
    }
    // a read error is the end of the data: the same answer as for the shorter input
    let short = &data[..9004];
    let reader = Chunked {
        data: short,
        chunk: 100,
        fail_after_data: true,
        interrupt_next: false,
    };
    assert_eq!(appears_fasta(reader), appears_fasta(short));
}

/// Reads every file of the TSV named by `LOSAT_FORMAT_GUESS_CHECK` (`path<TAB>fasta|other`, `other` = makeblastdb
/// warned) and compares. The TSV is built from makeblastdb 2.17.0's answers by `scripts/db2_fg_check.py`.
#[test]
#[ignore]
fn agrees_with_makeblastdb_on_listed_files() {
    let list = std::env::var("LOSAT_FORMAT_GUESS_CHECK")
        .expect("LOSAT_FORMAT_GUESS_CHECK must name a TSV file (path<TAB>fasta|other)");
    let text = std::fs::read_to_string(&list).expect("read the TSV");
    let mut disagree: Vec<String> = Vec::new();
    let mut total = 0usize;
    for line in text.lines().filter(|l| !l.is_empty()) {
        let (path, answer) = line.split_once('\t').expect("path<TAB>fasta|other");
        let want = match answer {
            "fasta" => true,
            "other" => false,
            other => panic!("bad answer {other:?} in {line:?}"),
        };
        let file = std::fs::File::open(path).unwrap_or_else(|e| panic!("open {path}: {e}"));
        let got = appears_fasta(io::BufReader::new(file));
        total += 1;
        if got != want {
            disagree.push(format!(
                "{path}: makeblastdb {answer}, port {}",
                if got { "fasta" } else { "other" }
            ));
        }
    }
    println!("{total} files, {} disagreements", disagree.len());
    for d in &disagree {
        println!("DISAGREE {d}");
    }
    assert!(disagree.is_empty());
}
