//! A merged run (`2>&1`) of LOSAT that NCBI's frozen merged stream can judge at NCBI's own stream
//! offsets.
//!
//! NCBI's diagnostics go to `cerr`, which flushes `cout` first, so a post lands after the bytes
//! that NCBI's output buffer had written by then; which bytes those are depends on how many came
//! before (the buffer's size and the calls that fill it, `report::serial::out_buffer`). LOSAT
//! states its own version where NCBI states `2.17.0+` (`losat::wording`, AGENTS.md rule 8), which
//! adds bytes before the posts, so LOSAT's buffer reaches its output at other bytes of the content
//! than NCBI's did. The merged run of such a row therefore reads a query whose title is shorter by
//! as many bytes as the version adds before the posts, spread over the title's occurrences there
//! (an archive writes it twice: the sequence's title and the defline): every later byte is where
//! NCBI's was, and NCBI's frozen stream, in LOSAT's words and with the shorter title, is what
//! LOSAT writes. No position is moved by a rule of the test's own.

/// The merged run of a row with a shorter query title, and what LOSAT writes then.
pub struct BalancedTitle {
    /// The query file of the merged run: the row's with its first title shorter.
    pub query: String,
    /// NCBI's frozen merged stream in LOSAT's words (`losat::wording::losat_output_words`, no
    /// position moved) with the shorter title in each of its occurrences.
    pub expected: Vec<u8>,
}

/// The encoded forms of `text` in a stream: ASN.1 text (`"text"`), ASN.1 binary (a VisibleString,
/// tag 26 and a length byte) and XML (`>text<`, `&`, `<` and `>` escaped).
fn forms(text: &[u8]) -> Vec<Vec<u8>> {
    let mut xml = Vec::with_capacity(text.len());
    for &byte in text {
        match byte {
            b'&' => xml.extend_from_slice(b"&amp;"),
            b'<' => xml.extend_from_slice(b"&lt;"),
            b'>' => xml.extend_from_slice(b"&gt;"),
            _ => xml.push(byte),
        }
    }
    let mut forms = vec![
        [b"\"", text, b"\""].concat(),
        [b">", &xml[..], b"<"].concat(),
    ];
    if let Ok(length) = u8::try_from(text.len()) {
        if length < 0x80 {
            forms.push([&[0x1a, length][..], text].concat());
        }
    }
    forms
}

fn count(haystack: &[u8], needle: &[u8]) -> usize {
    haystack
        .windows(needle.len())
        .filter(|window| *window == needle)
        .count()
}

fn replace(haystack: &[u8], from: &[u8], to: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(haystack.len());
    let mut i = 0;
    while i < haystack.len() {
        if haystack[i..].starts_with(from) {
            out.extend_from_slice(to);
            i += from.len();
        } else {
            out.push(haystack[i]);
            i += 1;
        }
    }
    out
}

/// The merged run of a row whose frozen merged stream is `ncbi_merged`, whose standard output
/// before its last post is `stdout_before_posts` (NCBI's bytes), and whose query file is `query`
/// (FASTA; with `-parse_deflines` the first word of a definition line is the identifier and the
/// rest the title): the first query's title is the one made shorter. `None` when LOSAT's words add nothing before the posts (the frozen stream in
/// LOSAT's words is then the expectation), when the title is also in other text of the stream, or
/// when the title cannot take up the added bytes evenly.
pub fn balanced_title(
    ncbi_merged: &[u8],
    stdout_before_posts: &[u8],
    query: &str,
    parse_deflines: bool,
) -> Option<BalancedTitle> {
    let added = losat::wording::losat_output_words(stdout_before_posts)
        .len()
        .checked_sub(stdout_before_posts.len())
        .filter(|&added| added > 0)?;
    let (first, rest) = query.split_once('\n')?;
    let defline = first.strip_prefix('>')?;
    let title = if parse_deflines {
        defline.split_once(' ').map_or("", |(_, title)| title)
    } else {
        defline
    };
    // The title is the first query's alone (another query with the same title would keep it).
    if title.is_empty() || count(query.as_bytes(), title.as_bytes()) != 1 {
        return None;
    }
    let line = format!(">{defline}");
    let strings: Vec<&[u8]> = vec![line.as_bytes(), title.as_bytes()];
    // Every occurrence of the title in the stream is one of these forms (the definition line
    // holds the title), so the shorter title replaces them all; a title also in other text (a
    // message that quotes it) leaves the row to the plain mapping.
    let in_forms: usize = strings
        .iter()
        .flat_map(|text| forms(text))
        .map(|form| count(ncbi_merged, &form))
        .sum();
    if count(ncbi_merged, title.as_bytes()) != in_forms {
        return None;
    }
    let occurrences: usize = strings
        .iter()
        .flat_map(|text| forms(text))
        .map(|form| count(stdout_before_posts, &form))
        .sum();
    if occurrences == 0 || added % occurrences != 0 {
        return None;
    }
    let shorter = added / occurrences;
    if !title.is_ascii() || shorter >= title.len() {
        return None;
    }
    let mut expected = losat::wording::losat_output_words(ncbi_merged).into_owned();
    for text in &strings {
        let short = &text[..text.len() - shorter];
        for (from, to) in forms(text).iter().zip(forms(short)) {
            expected = replace(&expected, from, &to);
        }
    }
    Some(BalancedTitle {
        query: format!(">{}\n{rest}", &defline[..defline.len() - shorter]),
        expected,
    })
}
