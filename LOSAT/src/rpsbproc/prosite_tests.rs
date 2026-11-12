//! `ProSite` against NCBI's `CProSite`: the results of `Parse` (error and position),
//! `GetMinimalXMap` (length) and `Match` that RpsbProc/src/common/prosite.cpp gave, run
//! unmodified (g++ 13.3.0, `g++ -O2 -std=c++17`) on 3,103 patterns and texts
//! (`LOSAT/tests/unit/rpsbproc_prosite_vectors.tsv`: pattern, text, `tlen`, start, parse error
//! (0 none, 1 invalid char, 2 syntax error, 3 unexpected end), error position, minimal map
//! length and match end or `npos`; the generator is `prositegen.cpp` and `gen.py` of the task
//! folder of item PS16b: random mixes of the grammar's tokens, including malformed ones).
//! The patterns that make NCBI loop for ever (an anchor with a star count) are not in it.

use super::*;

#[test]
fn parse_and_match_give_ncbis_results() {
    let table = include_str!("../../tests/unit/rpsbproc_prosite_vectors.tsv");
    let (mut valid, mut matched, mut count) = (0, 0, 0);
    for line in table.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        let pattern = fields[0].as_bytes();
        let text = fields[1].as_bytes();
        let tlen: usize = fields[2].parse().unwrap();
        let start: usize = fields[3].parse().unwrap();
        let error: u32 = fields[4].parse().unwrap();
        let position: usize = fields[5].parse().unwrap();
        count += 1;
        match ProSite::parse(pattern) {
            Err((kind, at)) => {
                let code = match kind {
                    ParseError::InvalidChar => 1,
                    ParseError::SyntaxError => 2,
                    ParseError::UnexpectedEnd => 3,
                };
                assert_eq!((code, at), (error, position), "{line}");
            }
            Ok(site) => {
                assert_eq!(error, 0, "{line}");
                valid += 1;
                assert_eq!(site.minimal_map_len().to_string(), fields[6], "{line}");
                let result = site.match_at(text, tlen, start).unwrap();
                let expected = match fields[7] {
                    "npos" => None,
                    end => Some(end.parse::<usize>().unwrap()),
                };
                assert_eq!(result, expected, "{line}");
                matched += usize::from(expected.is_some());
            }
        }
    }
    assert_eq!((count, valid, matched), (3103, 2603, 596));
}

// The pattern errors that NCBI's rpsbproc printed for the motifs of probe 3 (PS16b, ws): `K#R`
// position 1, `]K` position 0, `x(2,a)` position 4, `[K]]` position 3; and the patterns it
// accepted: `K-`, `x(`, `x(2`, `{K}`, `<K`, `K>`, `[KR` (the end of the expression is not
// checked against the state).
#[test]
fn probe_patterns() {
    for (pattern, kind, at) in [
        ("K#R", ParseError::InvalidChar, 1),
        ("]K", ParseError::SyntaxError, 0),
        ("x(2,a)", ParseError::SyntaxError, 4),
        ("[K]]", ParseError::SyntaxError, 3),
    ] {
        assert_eq!(ProSite::parse(pattern.as_bytes()).unwrap_err(), (kind, at));
    }
    for pattern in ["K-", "x(", "x(2", "{K}", "<K", "K>", "[KR"] {
        assert!(ProSite::parse(pattern.as_bytes()).is_ok(), "{pattern}");
    }
    assert_eq!(
        ProSite::parse(b"").unwrap_err().0,
        ParseError::UnexpectedEnd
    );
}

// Every construct once, with the match NCBI's code gives by reading it.
#[test]
fn constructs() {
    let site = |pattern: &str| ProSite::parse(pattern.as_bytes()).unwrap();
    let end = |pattern: &str, text: &str| {
        site(pattern)
            .match_at(text.as_bytes(), text.len(), 0)
            .unwrap()
    };
    // Letters, either case; `-` and white space are ignored.
    assert_eq!(end("D-e-K", "DEK"), Some(3));
    assert_eq!(end("D E\tK", "dek"), Some(3));
    // Alternatives, with a repeated letter, and the wildcard in them.
    assert_eq!(end("[KRK]-D", "RD"), Some(2));
    assert_eq!(end("[KX]-D", "AD"), Some(2));
    // Negation.
    assert_eq!(end("{P}-D", "KD"), Some(2));
    assert_eq!(end("{P}-D", "PD"), None);
    // Repeats and ranges, the reversed range swapped, the star count.
    assert_eq!(end("K(2)-D", "KKD"), Some(3));
    assert_eq!(end("K(2,3)-D", "KKKD"), Some(4));
    assert_eq!(end("K(3,2)-D", "KKKD"), Some(4));
    assert_eq!(end("x(0,*)-D", "AAAD"), Some(4));
    assert_eq!(end("K*-D", "KKD"), Some(3));
    // Anchors.
    assert_eq!(end("<K", "KA"), Some(1));
    assert_eq!(end("<K", "AK"), None);
    assert_eq!(end("K>", "AK"), None); // `>` matches at the end only; the match starts at 0
    assert_eq!(end("x>", "A"), Some(1));
    // Minimal map: the sum of the minimal counts of the positions that have alternatives.
    assert_eq!(site("[KR]-x(2)-D(3,5)-{P}(0,2)").minimal_map_len(), 6);
    assert_eq!(site("<K>").minimal_map_len(), 3);
}

// `tlen` is the length of the query, longer than the text of a site (basedata.cpp:731): a
// wildcard matches beyond the text, the character after it is the terminating 0 (no match), and
// the next one is a read of heap bytes, which `match_at` reports.
#[test]
fn match_beyond_the_text() {
    let site = ProSite::parse(b"x(2)-K").unwrap();
    assert_eq!(site.match_at(b"AAK", 100, 0), Ok(Some(3)));
    assert_eq!(site.match_at(b"AA", 100, 0), Ok(None));
    let site = ProSite::parse(b"x(3)-K").unwrap();
    assert_eq!(site.match_at(b"AA", 100, 0), Err(ReadPast));
    assert_eq!(site.match_at(b"AA", 0, 0), Ok(None));
}
