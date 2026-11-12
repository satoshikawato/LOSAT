//! `ProSite` against NCBI's `CProSite`: the results of `Parse` (error and position),
//! `GetMinimalXMap` (length) and `Match` that RpsbProc/src/common/prosite.cpp gave, run
//! unmodified (g++ 13.3.0, `g++ -O2 -std=c++17`) on 3,103 patterns and texts
//! (`LOSAT/tests/unit/rpsbproc_prosite_vectors.tsv`: pattern, text, `tlen`, start, parse error
//! (0 none, 1 invalid char, 2 syntax error, 3 unexpected end), error position, minimal map
//! length and match end or `npos`; the generator is `prositegen.cpp` and `gen.py` of the task
//! folder of item PS16b: random mixes of the grammar's tokens, including malformed ones).
//! The patterns that make NCBI loop for ever (an anchor with a star count) are not in it.

use super::*;

impl ProSite {
    /// NCBI's `Match` where it reads no byte beyond the text's terminating 0: the end of the
    /// match (`Err(ReadPast)` otherwise).
    fn match_at(
        &self,
        text: &[u8],
        tlen: usize,
        start_pos: usize,
    ) -> Result<Option<usize>, MatchError> {
        match self.run(text, tlen, start_pos, Past::Fail) {
            (_, true) => Err(MatchError::ReadPast),
            (result, false) => result,
        }
    }

    fn minimal_map_len(&self) -> u64 {
        self.minimal_map().unwrap().len
    }
}

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
    assert_eq!(site.match_at(b"AA", 100, 0), Err(MatchError::ReadPast));
    assert_eq!(site.match_at(b"AA", 0, 0), Ok(None));
}

// NCBI's `Match` (prosite.cpp, unmodified, g++ 13.3.0 `-O2 -std=c++17`) on texts shorter than
// `tlen`, each run with the bytes past the text's terminating 0 filled by 20 fixed bytes and 12
// seeded random fills (`LOSAT/tests/unit/rpsbproc_prosite_past_vectors.tsv`: pattern, text,
// `tlen`, start, parse error, error position, minimal map length, and the result over the fills:
// the end every fill gives, `ends` (every fill matches at ends that differ), `npos` (no fill
// matches), `mixed`, or `exhausted` (NCBI ended under a 1 GiB address-space cap or a 3 s
// timeout); the generator is `gen2.cpp` and `gen2.py` of the task folder of item PS16b, fix
// round). `matches` answers only where the fills agree, and refuses where they do not.
#[test]
fn the_two_runs_give_ncbis_answer_for_every_byte_past_the_text() {
    let table = include_str!("../../tests/unit/rpsbproc_prosite_past_vectors.tsv");
    let (mut count, mut answered, mut mixed, mut refused, mut exhausted) = (0, 0, 0, 0, 0);
    for line in table.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        let text = fields[1].as_bytes();
        let tlen: usize = fields[2].parse().unwrap();
        let start: usize = fields[3].parse().unwrap();
        count += 1;
        let site = match ProSite::parse(fields[0].as_bytes()) {
            Ok(site) => site,
            Err((_, at)) => {
                assert_ne!(fields[4], "0", "{line}");
                assert_eq!(at.to_string(), fields[5], "{line}");
                continue;
            }
        };
        assert_eq!(fields[4], "0", "{line}");
        let result = fields[7];
        // NCBI's buffer holds `tlen + 64` bytes: no read leaves it.
        let answer = site.matches(text, tlen, start, tlen as u64 + 63);
        if result == "exhausted" {
            assert_eq!(answer, Err(MatchError::Exhausted), "{line}");
            exhausted += 1;
            continue;
        }
        assert_eq!(site.minimal_map_len().to_string(), fields[6], "{line}");
        match answer {
            Ok(true) => assert!(result != "npos" && result != "mixed", "{line}"),
            Ok(false) => assert_eq!(result, "npos", "{line}"),
            Err(MatchError::ReadPast) => {
                if result == "mixed" {
                    mixed += 1;
                } else {
                    refused += 1;
                }
            }
            Err(MatchError::Exhausted) => panic!("{line}"),
        }
        answered += usize::from(answer.is_ok());
        // Without a read past the text, the end is NCBI's for every fill.
        if let Ok(end) = site.match_at(text, tlen, start) {
            let end = end.map_or("npos".to_string(), |end| end.to_string());
            assert_eq!(end, result, "{line}");
        }
    }
    // 53 parse errors. Of the 1,022 answers, 199 come from runs that read past the text; the
    // other 31 refusals depend on the bytes too, where the 32 fills happened to agree (NCBI's
    // matcher with 200,000 fills drawn from each pattern's letters matched for some and not for
    // others, and gave the answer of each of the 199 for all of them; probe of the fix round).
    assert_eq!(
        (count, answered, mixed, refused, exhausted),
        (1415, 1022, 298, 31, 11)
    );
}

// Review oracles R1-R4 (`Hx(0,10)K`, `H{P}(0,10)K`, `Hx(0,10)Y`, `Hx(0,10)W` on the residues
// `HKY` of a query of 185): NCBI prints the site of the first three; the fourth depends on the
// 8 bytes after the text's terminating 0, in the string's 16-byte local buffer (capacity 15;
// NCBI's two builds dropped the site in 10 runs each).
#[test]
fn review_oracles_r1_to_r4() {
    let matches = |pattern: &str, buffer: u64| {
        ProSite::parse(pattern.as_bytes())
            .unwrap()
            .matches(b"HKY", 185, 0, buffer)
    };
    assert_eq!(matches("Hx(0,10)K", 15), Ok(true));
    assert_eq!(matches("H{P}(0,10)K", 15), Ok(true));
    assert_eq!(matches("Hx(0,10)Y", 15), Ok(true));
    assert_eq!(matches("Hx(0,10)W", 15), Err(MatchError::ReadPast));
    // The aggressive `x(0,10)` reaches 11 first: a buffer that ends before it is refused.
    assert_eq!(matches("Hx(0,10)K", 10), Err(MatchError::ReadPast));
    assert_eq!(matches("Hx(0,9)K", 10), Ok(true));
    // A repeat that may not advance and could reach 2^26 + 1 records is refused past the text.
    assert_eq!(
        matches("Hx(0,10)K-{<}(67108864)", 15),
        Err(MatchError::ReadPast)
    );
    assert_eq!(matches("Hx(0,10)K-{<}(67108863)", 15), Ok(true));
}

// `GetMinimalXMap`'s string (libstdc++): `reserve(ttl + ttl)`, then one `append` per position;
// `max_size()` is 2^62 - 1, an allocation is the capacity + 1 bytes, doubled when it grows by
// less than twice. Review oracles R20 (`Hx(99999999999999999999)K`, `atol` saturates:
// `length_error`) and R21 (`Hx(200000000000000)K`: `bad_alloc`).
#[test]
fn the_minimal_map_follows_ncbis_string() {
    let map = |pattern: &str| ProSite::parse(pattern.as_bytes()).unwrap().minimal_map();
    let string = |len, capacity| Ok(CxxString { len, capacity });
    // In the local buffer; `reserve(16)` of eight positions doubles 15 to 30; 16 doubles too.
    assert_eq!(map("HKY"), string(3, 15));
    assert_eq!(map("ACDEFGHI"), string(8, 30));
    assert_eq!(map("x(16)"), string(16, 30));
    assert_eq!(map("x(31)-x(2)"), string(33, 62));
    assert_eq!(map("<-x(0,5)"), string(1, 15));
    // The band of the host's memory starts at an allocation of 2^30 bytes.
    assert_eq!(map("x(1073741822)"), string(1073741822, 1073741822));
    // `H`, then `x(n)` in an allocation of `n + 2` bytes, then `K` doubles it.
    assert_eq!(map("Hx(536870910)K"), string(536870912, 1073741822));
    assert_eq!(
        map("Hx(536870911)K"),
        Err(MapEnd::HostDependent((1 << 30) + 1))
    );
    assert_eq!(map("x(1073741823)"), Err(MapEnd::HostDependent(1 << 30)));
    assert_eq!(
        map("x(536870912)-x(536870912)"),
        Err(MapEnd::HostDependent((1 << 30) + 1))
    );
    // `bad_alloc` from 2^47 bytes on.
    assert_eq!(
        map("x(140737488355326)"),
        Err(MapEnd::HostDependent((1 << 47) - 1))
    );
    assert_eq!(map("x(140737488355327)"), Err(MapEnd::Throw(BAD_ALLOC)));
    assert_eq!(map("Hx(200000000000000)K"), Err(MapEnd::Throw(BAD_ALLOC)));
    // `length_error` past `max_size()`, before any allocation.
    assert_eq!(map("x(4611686018427387903)"), Err(MapEnd::Throw(BAD_ALLOC)));
    assert_eq!(
        map("x(4611686018427387904)"),
        Err(MapEnd::Throw(LENGTH_ERROR))
    );
    assert_eq!(
        map("Hx(99999999999999999999)K"),
        Err(MapEnd::Throw(LENGTH_ERROR))
    );
}

// The records of a repeat that does not advance: `vector<size_t>` doubles from empty, and one
// allocation of 2^30 bytes or more (2^26 + 1 records) is refused; below, the match goes on.
#[test]
fn repeats_that_do_not_advance_refuse_one_gib_of_records() {
    let end = |pattern: &str| {
        ProSite::parse(pattern.as_bytes())
            .unwrap()
            .matches(b"HKY", 3, 0, 15)
    };
    assert_eq!(end("<(67108863)-H"), Ok(true));
    assert_eq!(end("<(67108864)-H"), Err(MatchError::Exhausted));
    assert_eq!(end("<*-H"), Err(MatchError::Exhausted));
    assert_eq!(end("H-{<}(67108863)-K"), Ok(true));
    assert_eq!(end("H-{<}(67108864)-K"), Err(MatchError::Exhausted));
    // A lazy `{X}*` at `tlen` is revisited for ever.
    assert_eq!(end("x(3)-{X}*-Q"), Err(MatchError::Exhausted));
    // The revisits of 2^26 - 1 records that did not advance end at once (NCBI pops them one by
    // one, and the rest of the pattern fails from the same position each time).
    assert_eq!(end("<(0,67108863)-Q"), Ok(false));
    assert_eq!(end("<(0,67108863)-H"), Ok(true));
    // Two of them (`<` only begins a pattern; `{<}` stays after its first position).
    assert_eq!(end("H-{<}(0,67108863)-{<}(0,67108863)-Y"), Ok(false));
    assert_eq!(end("H-{<}(0,67108863)-{<}(0,67108863)-K"), Ok(true));
    // The capacity that a copy back sets: growth from 3 is 3 * 2^k.
    let mut capacity = 3;
    assert_eq!(grow_records(&mut capacity, 100_663_296), Ok(()));
    assert_eq!(capacity, 100_663_296);
    assert_eq!(
        grow_records(&mut capacity, 100_663_297),
        Err(MatchError::Exhausted)
    );
    let mut capacity = 0;
    assert_eq!(grow_records(&mut capacity, 1 << 26), Ok(()));
    assert_eq!(
        grow_records(&mut capacity, (1 << 26) + 1),
        Err(MatchError::Exhausted)
    );
}
