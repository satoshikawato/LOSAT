//! `CProSite` of rpsbproc's library: the parser and matcher of the PROSITE-like patterns of
//! `cddannot.dat` / `cddannot_generic.dat` (`[KR]`, `x(2)`, `{P}`, `<` and `>` anchors, `-`
//! separators), as NCBI's `Parse`, `Match` and `GetMinimalXMap` read and run them, including the
//! parser's error returns and the quirks of its state machine (a pattern that ends inside a
//! bracket or a repeat is accepted; an escape stays set across `(` and `)`).
//!
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.hpp:7-66
//!
//! NCBI's `Match` reads `text[pos]` for every `pos < tlen`, and rpsbproc passes the length of the
//! whole query as `tlen` while the text holds only the mapped residues of a site (basedata.cpp:731).
//! A read at `text.size()` is the string's terminating 0 (defined); a read beyond it is
//! undefined in C++ (bytes of NCBI's memory). `ProSite::matches` gives NCBI's answer where it does
//! not depend on those bytes (DECISIONS 37, Owner-delegated 2026-11-12): the match is run once
//! with every such read failing and once with every such read passing; NCBI's search tries every
//! count of every position, so a match of the first run is a match for any bytes and a failure
//! of the second run a failure for any bytes. The two runs differing, or a read that may leave
//! NCBI's own string buffer (memory that may not be mapped), is refused (`MatchError::ReadPast`).
//!
//! NCBI's `Match` stores one record (`size_t`) per step of a position in a vector; a step that
//! does not advance (an anchor, a negation of anchors only) repeats at the same position up to the
//! maximal count, and a lazy position that does not advance is revisited for ever, until the
//! memory runs out. `Record` keeps the counts instead of the records. `ProSite::minimal_map`
//! follows `GetMinimalXMap`'s string by its length and capacity, with the exceptions that its
//! `append` throws (`MapEnd`). Where NCBI would hold 1 GiB or more at once for the motif check
//! (the string's old and new buffers while it grows; the vector of records, the copies of it on
//! the match's stack and the string that holds the text: `Live`), LOSAT refuses
//! (`MapEnd::HostDependent`, `MatchError::Exhausted`): NCBI's result then depends on the host's
//! memory.

/// NCBI's live memory of the motif check of this many bytes or more is refused: below the
/// allocation of 2^47 bytes, NCBI's allocations succeed or throw `std::bad_alloc` by the host's
/// memory and overcommit (DECISIONS 37 and 40, Owner-delegated 2026-11-12: NCBI's rpsbproc uses
/// about 130 MB, and a 64-bit host that processes CDD gives it 1 GiB; the bound counts what NCBI
/// holds at once, not one allocation: review 2 oracle O3, where 1.5 GiB of records in
/// allocations of 512 MB gave exit 0 under a 4 GiB address space and `std::bad_alloc` under
/// 1.5 GiB).
pub(crate) const HOST_DEPENDENT_BYTES: u64 = 1 << 30;

/// The user address space of x86-64 (4-level paging): `operator new` of this many bytes or more
/// fails on every host and throws `std::bad_alloc`.
pub(crate) const USER_SPACE_BYTES: u64 = 1 << 47;

// The growth of the C++ library's string and vector in NCBI's motif check (`CxxString`, `Live::grow`).
#[path = "prosite_growth.rs"]
mod growth;
#[cfg(test)]
pub(crate) use growth::{events as growth_events, BAD_ALLOC, LENGTH_ERROR};
pub(crate) use growth::{CxxString, MapEnd};

/// `CProSite::EParseError` without `eNoError`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ParseError {
    InvalidChar,
    SyntaxError,
    UnexpectedEnd,
}

const NEGATIVE_FILTER: u32 = 0x1;
const LAZY_MATCH: u32 = 0x2;
/// A `size_t` cast from -1: the biggest possible count. The counts are NCBI's `size_t` of x86-64,
/// `u64` on every target.
const INFINITE: u64 = u64::MAX;

/// `CProSite::TPatternPos`.
#[derive(Clone, Debug)]
struct PatternPos {
    flags: u32,
    alt_chars: Vec<u8>,
    min_count: u64,
    max_count: u64,
}

impl PatternPos {
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:11-13
    fn new() -> Self {
        PatternPos {
            flags: 0,
            alt_chars: Vec::new(),
            min_count: 1,
            max_count: 1,
        }
    }

    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:501-553
    fn normalize_alt_chars(&mut self) {
        if self.alt_chars.is_empty() {
            return;
        }
        let mut idx = 0;
        let mut total = self.alt_chars.len();
        while idx < total {
            let current = self.alt_chars[idx];
            if current == b'X' {
                self.alt_chars.clear();
                self.alt_chars.push(b'X');
                break;
            }
            let mut removed = 0;
            let mut idx2 = idx + 1;
            while idx2 < total {
                if self.alt_chars[idx2] == current {
                    removed += 1;
                } else if removed > 0 {
                    self.alt_chars[idx2 - removed] = self.alt_chars[idx2];
                }
                idx2 += 1;
            }
            total -= removed;
            self.alt_chars.truncate(total);
            idx += 1;
        }
    }

    fn negative(&self) -> bool {
        self.flags & NEGATIVE_FILTER != 0
    }
}

/// Why `ProSite::matches` or `ProSite::match_at` gives no answer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum MatchError {
    /// NCBI's answer depends on the bytes it reads beyond the text's terminating 0, or the
    /// match may read beyond NCBI's buffer of the text, or may hold `HOST_DEPENDENT_BYTES` for
    /// other bytes there (see the module comment).
    ReadPast,
    /// NCBI's live memory of the match (`Live`) reaches `HOST_DEPENDENT_BYTES`.
    Exhausted,
}

/// What a read beyond the text's terminating 0 gives in one of the two runs of `matches`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Past {
    /// The letter does not match; a negation fails.
    Fail,
    /// Any letter matches; a negation passes.
    Pass,
}

/// The text of a match with NCBI's `tlen` (the length that bounds the positions) and the read of
/// `text[pos]`.
struct Matcher<'a> {
    text: &'a [u8],
    tlen: usize,
    past: Past,
    read_past: std::cell::Cell<bool>,
}

impl Matcher<'_> {
    /// `text[pos]` of a `std::string`: the terminating 0 at `size()`; `None` beyond it.
    fn at(&self, pos: usize) -> Option<u8> {
        match pos.cmp(&self.text.len()) {
            std::cmp::Ordering::Less => Some(self.text[pos]),
            std::cmp::Ordering::Equal => Some(0),
            std::cmp::Ordering::Greater => {
                self.read_past.set(true);
                None
            }
        }
    }
}

/// The bytes of one match record, the size of `size_t`.
const RECORD_BYTES: u64 = 8;

/// What NCBI holds at once while its `Match` runs: the capacity of `curr.flex` (the one vector
/// that every position's `rMatchRec` uses), the records of the copies of it on `stkFlexStack`,
/// and the heap buffer of `minMap`, the string that holds the text (DECISIONS 40). Each
/// allocation for the records is checked with what lives during it; NCBI's live memory of
/// `HOST_DEPENDENT_BYTES` or more is `MatchError::Exhausted`.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:471-492
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:723-731
// (`x_TMatchRec` is the vector `flex` and the iterator, prosite.hpp:54-60. The stack's deque of
// 32-byte entries, at most one per position of the pattern, and the compiled pattern are not
// counted.)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Live {
    /// `curr.flex.capacity()`.
    capacity: u64,
    /// The records of the copies on the stack (a copy's capacity is its size).
    stacked: u64,
    /// `minMap`'s bytes on the heap.
    text: u64,
    /// The most records that the stack held in the revisits that LOSAT skipped (`retry`) since
    /// `ProSite::run` last took them: the copies of a revisited position and of the rest of the
    /// pattern above it, which an earlier position's revisits hold again.
    skipped: u64,
}

impl Live {
    fn new(text: u64) -> Self {
        Live {
            capacity: 0,
            stacked: 0,
            text,
            skipped: 0,
        }
    }

    /// NCBI holds `records` match records and the text.
    fn hold(&self, records: u64) -> Result<(), MatchError> {
        let bytes = u128::from(records) * u128::from(RECORD_BYTES) + u128::from(self.text);
        if bytes >= u128::from(HOST_DEPENDENT_BYTES) {
            Err(MatchError::Exhausted)
        } else {
            Ok(())
        }
    }

    /// A revisit of NCBI's that LOSAT does not make (`PatternPos::match_at`): the rest of the
    /// pattern runs again from the same position, as it has just run and failed, with a copy of
    /// `copy` records of the revisited position on the stack (0: none) and, at its most, `rest`
    /// records of its own copies above it. Its growths happened in the run before (the capacity
    /// never shrinks), so it holds the capacity, the stack and those copies at once. The stack
    /// it holds stays in `skipped`: a revisit of a position below runs these revisits again.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:480-492 (the
    // revisit that does not advance pushes its copy when counts are left, and the rest of the
    // pattern starts again from `curr.flex.back()`)
    fn retry(&mut self, copy: u64, rest: u64) -> Result<(), MatchError> {
        let stacked = self.stacked.saturating_add(copy).saturating_add(rest);
        self.skipped = self.skipped.max(stacked);
        self.hold(self.capacity.saturating_add(stacked))
    }
}

/// `rMatchRec` of one position, by its counts: the records are `start`, then one per step,
/// each the position after the step. The steps that advance (by one residue) come first; a
/// step that does not advance (an anchor, or a negation of anchors and `X` only) depends only
/// on the position, so every later step repeats it: those are counted in `still`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Record {
    /// `rMatchRec.size()`; 0 is empty.
    len: u64,
    start: usize,
    still: u64,
}

impl Record {
    /// `rMatchRec.back()`.
    fn last(&self) -> usize {
        // The steps that advanced are fewer than `tlen`.
        self.start + (self.len - 1 - self.still) as usize
    }

    /// `push_back` of `start_pos` on `rMatchRec`, at the first visit.
    fn begin(&mut self, start: usize) {
        *self = Record {
            len: 1,
            start,
            still: 0,
        };
    }

    /// The record of a step that advanced.
    fn advance(&mut self) {
        debug_assert_eq!(self.still, 0, "a step at a position that did not advance");
        self.len += 1;
    }

    /// The records of steps that do not advance, up to `len` records.
    fn stay_until(&mut self, len: u64) {
        self.still += len - self.len;
        self.len = len;
    }

    /// `rMatchRec.pop_back()`.
    fn pop(&mut self) {
        self.still = self.still.saturating_sub(1);
        self.len -= 1;
    }
}

impl PatternPos {
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:554-592
    fn neg_match(&self, m: &Matcher, pos: usize) -> Option<usize> {
        let mut adv = 0;
        for &ch in &self.alt_chars {
            match ch {
                b'<' => {
                    if pos == 0 {
                        return None;
                    }
                }
                b'>' => {
                    if pos == m.tlen {
                        return None;
                    }
                }
                b'X' => {
                    if pos < m.tlen {
                        return None;
                    }
                }
                _ => {
                    if pos < m.tlen {
                        adv = 1;
                        match m.at(pos) {
                            Some(mut c) => {
                                if c.is_ascii_lowercase() {
                                    c -= 0x20;
                                }
                                if ch == c {
                                    return None;
                                }
                            }
                            None => {
                                if m.past == Past::Fail {
                                    return None;
                                }
                            }
                        }
                    } else {
                        return None;
                    }
                }
            }
        }
        Some(adv)
    }

    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:594-627
    fn pos_match(&self, m: &Matcher, pos: usize) -> Option<usize> {
        for &ch in &self.alt_chars {
            match ch {
                b'<' => {
                    if pos == 0 {
                        return Some(0);
                    }
                }
                b'>' => {
                    if pos == m.tlen {
                        return Some(0);
                    }
                }
                b'X' => {
                    if pos < m.tlen {
                        return Some(1);
                    }
                }
                _ => {
                    if pos < m.tlen {
                        match m.at(pos) {
                            Some(mut c) => {
                                if c.is_ascii_lowercase() {
                                    c -= 0x20;
                                }
                                if ch == c {
                                    return Some(1);
                                }
                            }
                            None => {
                                if m.past == Past::Pass {
                                    return Some(1);
                                }
                            }
                        }
                    }
                }
            }
        }
        None
    }

    /// `m_pfnMatch`.
    fn step(&self, m: &Matcher, pos: usize) -> Option<usize> {
        if self.negative() {
            self.neg_match(m, pos)
        } else {
            self.pos_match(m, pos)
        }
    }

    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:630-698
    // (A lazy position whose `matched` is above `m_ulMaxCount` on a revisit leaves `retVal` 0
    // and the record as it is.) `live` follows NCBI's memory (`Live`): `rMatchRec` is
    // `curr.flex`.
    //
    // A step that does not advance gives the same at every later step (it depends only on the
    // position): NCBI's loops repeat it up to the count they stop at, which `Record::stay_until`
    // counts at once. A revisit starts where the rest of the pattern has just failed: one that
    // does not advance (a lazy step, or the pop of a record of a step that did not advance)
    // leaves the rest to fail again from the same position, so NCBI revisits until the count
    // runs out (a lazy `{X}*` for ever) and then fails; LOSAT goes there at once, and checks the
    // memory of the revisits it skips (`Live::retry`; `rest` is the most records that the rest
    // of the pattern's copies held at once in the run that failed).
    fn match_at(
        &self,
        m: &Matcher,
        record: &mut Record,
        live: &mut Live,
        mut start_pos: usize,
        rest: u64,
    ) -> Result<i32, MatchError> {
        let matched = record.len;
        let ret;
        if self.flags & LAZY_MATCH != 0 {
            if matched == 0 {
                record.begin(start_pos);
                live.grow(1)?;
                while record.len <= self.min_count {
                    let Some(adv) = self.step(m, start_pos) else {
                        return Ok(-1);
                    };
                    if adv == 0 {
                        let len = self.min_count.saturating_add(1);
                        live.grow(len)?;
                        record.stay_until(len);
                        break;
                    }
                    start_pos += adv;
                    record.advance();
                    live.grow(record.len)?;
                }
                ret = i32::from(self.max_count > self.min_count);
            } else if matched <= self.max_count {
                start_pos = record.last();
                match self.step(m, start_pos) {
                    None => ret = -1,
                    Some(0) => {
                        // NCBI's revisits add a record each up to `m_ulMaxCount` + 1; while
                        // counts are left they copy the records onto the stack, the most at
                        // `m_ulMaxCount` records, and the rest of the pattern fails again.
                        if matched < self.max_count {
                            live.grow(self.max_count)?;
                            live.retry(self.max_count, rest)?;
                        }
                        live.grow(self.max_count.saturating_add(1))?;
                        live.retry(0, rest)?;
                        ret = -1;
                    }
                    Some(adv) => {
                        debug_assert_eq!(adv, 1);
                        record.advance();
                        live.grow(record.len)?;
                        ret = i32::from(self.max_count > matched);
                    }
                }
            } else {
                ret = 0;
            }
        } else {
            let need = self.min_count.wrapping_add(1);
            if matched == 0 {
                record.begin(start_pos);
                live.grow(1)?;
                while record.len <= self.max_count {
                    let Some(adv) = self.step(m, start_pos) else {
                        break;
                    };
                    if adv == 0 {
                        let len = self.max_count.saturating_add(1);
                        live.grow(len)?;
                        record.stay_until(len);
                        break;
                    }
                    start_pos += adv;
                    record.advance();
                    live.grow(record.len)?;
                }
            } else {
                if record.still > 0 {
                    // The pops down to the advanced records keep `rMatchRec.back()`; they end
                    // the position (`retVal` 0 or -1 at `need`) or reach the advanced records.
                    // The first of them holds the most: a copy of the records left while they
                    // are more than `need` (`retVal` 1).
                    let first = record.len - 1;
                    live.retry(if first > need { first } else { 0 }, rest)?;
                    let advanced = record.len - record.still;
                    if need >= advanced {
                        return Ok(-1);
                    }
                    record.len = advanced;
                    record.still = 0;
                }
                record.pop();
            }
            let matched = record.len;
            ret = if matched < need {
                -1
            } else if matched > need {
                1
            } else {
                0
            };
        }
        Ok(ret)
    }

    /// Whether a step can succeed without advancing: a positive position with an anchor (`<` at
    /// 0, `>` at `tlen`), or a negation of anchors and `X` only.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:558-604,
    // 602-608 (`x_NegMatch`: `adv` becomes 1 only in `default`; `x_PosMatch`: an anchor
    // returns 0)
    fn may_stay(&self) -> bool {
        if self.negative() {
            !self.reads()
        } else {
            self.alt_chars.iter().any(|&c| c == b'<' || c == b'>')
        }
    }

    /// Whether a step reads `text` (an alternative that is not an anchor or `X`).
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:576-579,
    // 614-617 (only the `default` alternatives of `x_NegMatch` and `x_PosMatch` read `text`)
    fn reads(&self) -> bool {
        self.alt_chars
            .iter()
            .any(|&c| c != b'<' && c != b'>' && c != b'X')
    }

    /// Whether a step can advance: a positive position with an alternative that is not an
    /// anchor, or a negation that reads.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:575-612,
    // 610-620 (`x_NegMatch` returns `adv` = 1 only after a `default` alternative; `x_PosMatch`
    // returns 1 for `X` and the `default` letters)
    fn may_advance(&self) -> bool {
        if self.negative() {
            self.reads()
        } else {
            self.alt_chars.iter().any(|&c| c != b'<' && c != b'>')
        }
    }
}

/// `CProSite`: the compiled pattern.
#[derive(Clone, Debug, Default)]
pub(crate) struct ProSite {
    pattern: Vec<PatternPos>,
}

/// The parser's states (`__EParserStatus`).
#[derive(Clone, Copy, PartialEq, Eq)]
enum Status {
    Ready,
    PlainText,
    Alternates,
    AltRedund,
    Negations,
    NegRedund,
    Repeats,
    Range,
}

/// `atol` of the digits, cast to `size_t`: glibc's `strtol` saturates at `LONG_MAX`.
fn atol(digits: &[u8]) -> u64 {
    let mut value: i64 = 0;
    for &digit in digits {
        value = value
            .saturating_mul(10)
            .saturating_add(i64::from(digit - b'0'));
    }
    value as u64
}

impl ProSite {
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:49-466
    // (The cases that the quote shortens are ported from prosite.cpp:159-461, one branch for
    // one branch. The end of the expression is not checked against the state: `[KR` and `x(2`
    // parse without an error.) `Err((error, position))` carries `errorPos`.
    pub(crate) fn parse(expr: &[u8]) -> Result<ProSite, (ParseError, usize)> {
        use Status::*;
        if expr.is_empty() {
            // `errorPos` is not set by NCBI here; rpsbproc never parses an empty motif.
            return Err((ParseError::UnexpectedEnd, 0));
        }
        let mut status = Ready;
        let mut pattern: Vec<PatternPos> = Vec::new();
        // `iterCurr`: the position under construction (valid whenever the code uses it).
        let mut cur = 0usize;
        let mut lit: Vec<u8> = Vec::new();
        let mut escaped = false;
        let mut inf_lit = false;
        let mut at = 0usize;
        macro_rules! fail {
            ($error:expr) => {
                return Err(($error, at))
            };
        }
        // `emplace` of `__dummy` at the end, and the assignment to `iterCurr` that follows.
        macro_rules! new_pos {
            () => {{
                pattern.push(PatternPos::new());
                cur = pattern.len() - 1;
            }};
        }
        while at < expr.len() {
            let mut c = expr[at];
            if c.is_ascii_lowercase() {
                c -= 0x20;
            }
            // The cases `*` (escaped) and the letters share one branch in NCBI (fall through).
            let mut letter = matches!(
                c,
                b'A' | b'C'
                    | b'D'
                    | b'E'
                    | b'F'
                    | b'G'
                    | b'H'
                    | b'I'
                    | b'K'
                    | b'L'
                    | b'M'
                    | b'N'
                    | b'P'
                    | b'Q'
                    | b'R'
                    | b'S'
                    | b'T'
                    | b'U'
                    | b'V'
                    | b'W'
                    | b'Y'
            );
            if c == b'*' {
                if escaped {
                    letter = true;
                } else {
                    match status {
                        PlainText => {
                            pattern[cur].min_count = 0;
                            pattern[cur].max_count = INFINITE;
                            if pattern[cur].alt_chars == b"X" {
                                pattern[cur].flags |= LAZY_MATCH;
                            }
                        }
                        Alternates | Negations => pattern[cur].alt_chars.push(c),
                        Range | Repeats => {
                            if lit.is_empty() && !inf_lit {
                                inf_lit = true;
                            }
                        }
                        AltRedund | NegRedund => {}
                        Ready => fail!(ParseError::SyntaxError),
                    }
                    at += 1;
                    continue;
                }
            }
            if letter {
                match status {
                    Ready | PlainText => {
                        status = PlainText;
                        new_pos!();
                        pattern[cur].alt_chars.push(c);
                    }
                    Alternates | Negations => pattern[cur].alt_chars.push(c),
                    _ => fail!(ParseError::SyntaxError),
                }
                escaped = false;
                at += 1;
                continue;
            }
            match c {
                b'X' => {
                    match status {
                        Ready | PlainText => {
                            status = PlainText;
                            new_pos!();
                            pattern[cur].alt_chars.push(c);
                        }
                        Alternates | Negations => pattern[cur].alt_chars.push(c),
                        _ => fail!(ParseError::SyntaxError),
                    }
                    escaped = false;
                }
                b'0'..=b'9' => {
                    match status {
                        Repeats | Range if !inf_lit => lit.push(c),
                        AltRedund | NegRedund => {}
                        _ => fail!(ParseError::SyntaxError),
                    }
                    escaped = false;
                }
                b'(' => match status {
                    PlainText => {
                        status = Repeats;
                        lit.clear();
                    }
                    Alternates => status = AltRedund,
                    Negations => status = NegRedund,
                    _ => fail!(ParseError::SyntaxError),
                },
                b')' => match status {
                    Repeats => {
                        if !lit.is_empty() {
                            let n = atol(&lit);
                            pattern[cur].min_count = n;
                            pattern[cur].max_count = n;
                        } else if inf_lit {
                            pattern[cur].min_count = 0;
                            pattern[cur].max_count = INFINITE;
                            if pattern[cur].alt_chars == b"X" {
                                pattern[cur].flags |= LAZY_MATCH;
                            }
                            inf_lit = false;
                        } else {
                            fail!(ParseError::SyntaxError);
                        }
                        status = PlainText;
                    }
                    Range => {
                        if !lit.is_empty() {
                            pattern[cur].max_count = atol(&lit);
                            let pos = &mut pattern[cur];
                            if pos.max_count < pos.min_count {
                                std::mem::swap(&mut pos.max_count, &mut pos.min_count);
                            }
                        } else if inf_lit {
                            pattern[cur].max_count = INFINITE;
                            if pattern[cur].alt_chars == b"X" {
                                pattern[cur].flags |= LAZY_MATCH;
                            }
                            inf_lit = false;
                        } else {
                            fail!(ParseError::SyntaxError);
                        }
                        status = PlainText;
                    }
                    AltRedund => status = Alternates,
                    NegRedund => status = Negations,
                    _ => fail!(ParseError::SyntaxError),
                },
                b',' => {
                    if escaped {
                        fail!(ParseError::InvalidChar);
                    }
                    match status {
                        Repeats if !lit.is_empty() => {
                            pattern[cur].min_count = atol(&lit);
                            lit.clear();
                            status = Range;
                        }
                        AltRedund | NegRedund => {}
                        _ => fail!(ParseError::SyntaxError),
                    }
                }
                b'\\' => {
                    if escaped {
                        fail!(ParseError::InvalidChar);
                    }
                    escaped = true;
                }
                b'[' => {
                    if escaped {
                        fail!(ParseError::InvalidChar);
                    }
                    match status {
                        PlainText | Ready => {
                            new_pos!();
                            status = Alternates;
                        }
                        _ => fail!(ParseError::SyntaxError),
                    }
                }
                b']' => {
                    if escaped {
                        fail!(ParseError::InvalidChar);
                    }
                    match status {
                        Alternates if !pattern[cur].alt_chars.is_empty() => {
                            pattern[cur].normalize_alt_chars();
                            status = PlainText;
                        }
                        _ => fail!(ParseError::SyntaxError),
                    }
                }
                b'{' => {
                    if escaped {
                        fail!(ParseError::InvalidChar);
                    }
                    match status {
                        PlainText | Ready => {
                            new_pos!();
                            pattern[cur].flags |= NEGATIVE_FILTER;
                            status = Negations;
                        }
                        _ => fail!(ParseError::SyntaxError),
                    }
                }
                b'}' => {
                    if escaped {
                        fail!(ParseError::InvalidChar);
                    }
                    match status {
                        Negations if !pattern[cur].alt_chars.is_empty() => {
                            pattern[cur].normalize_alt_chars();
                            status = PlainText;
                        }
                        _ => fail!(ParseError::SyntaxError),
                    }
                }
                b'<' => {
                    if escaped {
                        fail!(ParseError::InvalidChar);
                    }
                    match status {
                        Alternates => {
                            // Not the first position of the pattern.
                            if cur != 0 {
                                fail!(ParseError::SyntaxError);
                            }
                            pattern[cur].alt_chars.push(c);
                        }
                        Negations => pattern[cur].alt_chars.push(c),
                        Ready => {
                            new_pos!();
                            pattern[cur].alt_chars.push(c);
                            status = PlainText;
                        }
                        _ => fail!(ParseError::SyntaxError),
                    }
                }
                b'>' => {
                    if escaped {
                        fail!(ParseError::InvalidChar);
                    }
                    match status {
                        PlainText | Ready => {
                            new_pos!();
                            pattern[cur].alt_chars.push(c);
                            status = PlainText;
                        }
                        Negations => pattern[cur].alt_chars.push(c),
                        _ => fail!(ParseError::SyntaxError),
                    }
                }
                b'-' | b' ' | b'\t' => {}
                _ => fail!(ParseError::InvalidChar),
            }
            at += 1;
        }
        Ok(ProSite { pattern })
    }

    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:468-500
    // `Ok(None)` is `string::npos`. A pattern without a position (only `-`) matches at
    // `start_pos`. The second value says whether the match read `text` beyond its terminating 0
    // (answered by `past`). `text_heap` is the heap bytes of NCBI's string of the text; `Live`
    // follows `curr.flex` and the copies of it on the stack.
    fn run(
        &self,
        text: &[u8],
        tlen: usize,
        mut start_pos: usize,
        past: Past,
        text_heap: u64,
    ) -> (Result<Option<usize>, MatchError>, bool) {
        if start_pos >= tlen {
            return (Ok(None), false);
        }
        #[cfg(test)]
        if past == Past::Fail {
            growth::events::match_starts();
        }
        let matcher = Matcher {
            text,
            tlen,
            past,
            read_past: std::cell::Cell::new(false),
        };
        // Each entry keeps the most records that the stack held while it was on it (in NCBI's
        // revisits that LOSAT skipped too).
        let mut stack: Vec<(usize, Record, u64)> = Vec::new();
        let mut pos = 0usize;
        let mut flex = Record::default();
        let mut live = Live::new(text_heap);
        let mut rest = 0u64;
        let result = loop {
            if pos == self.pattern.len() {
                break Ok(Some(start_pos));
            }
            let matched =
                match self.pattern[pos].match_at(&matcher, &mut flex, &mut live, start_pos, rest) {
                    Ok(matched) => matched,
                    Err(error) => break Err(error),
                };
            // The revisits that LOSAT skipped happen while the positions below stay on the
            // stack: their own revisits hold these copies again.
            let skipped = std::mem::take(&mut live.skipped);
            if let Some(below) = stack.last_mut() {
                below.2 = below.2.max(skipped);
            }
            if matched < 0 {
                let Some((top_pos, top_flex, top_most)) = stack.pop() else {
                    break Ok(None);
                };
                // The copies above it in the run of the rest of the pattern that failed.
                rest = top_most - live.stacked;
                if let Some(below) = stack.last_mut() {
                    below.2 = below.2.max(top_most);
                }
                if let Err(error) = live.pop(top_flex.len) {
                    break Err(error);
                }
                pos = top_pos;
                flex = top_flex;
            } else {
                rest = 0;
                if matched > 0 {
                    if let Err(error) = live.push(flex.len) {
                        break Err(error);
                    }
                    stack.push((pos, flex, live.stacked));
                }
                start_pos = flex.last();
                pos += 1;
                flex.len = 0;
            }
        };
        (result, matcher.read_past.get())
    }

    /// Whether NCBI's `Match` finds a match (`MotifCheck` asks no more), also where it reads
    /// bytes beyond the text's terminating 0 when its answer does not depend on them (the module
    /// comment; DECISIONS 37). `buffer` is NCBI's string that holds the text: its bytes
    /// `0..=capacity` are NCBI's own, a read beyond them may fault (refused even where both
    /// runs agree, DECISIONS 40), and its heap bytes count in NCBI's live memory.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:723-731
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:468-489,
    // 480-489 (every position is tried with every count before the match fails: a failure
    // pops the stack and revisits the last position that has counts left)
    pub(crate) fn matches(
        &self,
        text: &[u8],
        tlen: usize,
        start_pos: usize,
        buffer: CxxString,
    ) -> Result<bool, MatchError> {
        let text_heap = buffer.heap_bytes();
        let (result, read_past) = self.run(text, tlen, start_pos, Past::Fail, text_heap);
        if !read_past {
            return result.map(|end| end.is_some());
        }
        // The runs that NCBI makes for other bytes go through other positions and counts:
        // every read must stay in NCBI's buffer, and no run may hold the memory that
        // `Exhausted` refuses.
        if self.may_exhaust(tlen, text_heap) || self.max_read(tlen, start_pos) > buffer.capacity {
            return Err(MatchError::ReadPast);
        }
        if let Ok(Some(_)) = result {
            return Ok(true);
        }
        #[cfg(test)]
        growth::events::pause(true);
        let second = self.run(text, tlen, start_pos, Past::Pass, text_heap).0;
        #[cfg(test)]
        growth::events::pause(false);
        match second {
            Ok(None) => Ok(false),
            _ => Err(MatchError::ReadPast),
        }
    }

    /// The highest position that a step that reads can be at, for any bytes of the text: each
    /// step advances by one at most, each position up to its maximal count, a step that reads
    /// advances, and reads stop before `tlen`.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:575-612,
    // 610-620 (a step returns 0, 1 or `npos`; `text[pos]` is read only below `tlen`)
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:678-683 (an
    // aggressive position steps up to `m_ulMaxCount` times; a lazy one adds a step per revisit
    // while `matched` is at most `m_ulMaxCount`, :659)
    fn max_read(&self, tlen: usize, start_pos: usize) -> u64 {
        let advance = self
            .pattern
            .iter()
            .filter(|pos| pos.may_advance())
            .fold(0u64, |sum, pos| sum.saturating_add(pos.max_count));
        (start_pos as u64)
            .saturating_add(advance)
            .saturating_sub(1)
            .min((tlen as u64).saturating_sub(1))
    }

    /// Whether NCBI's live memory of the match (`Live`) can reach `HOST_DEPENDENT_BYTES` for some
    /// bytes of the text. A position holds at most `R` records: one that may not advance up to
    /// its maximal count + 1 (a lazy one for ever), one that advances up to `tlen` + 1. Its last
    /// growth holds the old buffer and the new one (3/2 of the power of two that holds `R`), with
    /// the stack's copies of the earlier positions (the stack holds positions in their order, at
    /// most one copy each, and only a position with counts left is copied: a maximal count above
    /// the minimal). A copy onto the stack adds to `curr.flex`'s capacity, the power of two that
    /// holds the largest `R`.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:656-694,
    // 679, 693-695 (`retVal` is 1, a copy onto the stack, only with counts left)
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:489-492
    fn may_exhaust(&self, tlen: usize, text_heap: u64) -> bool {
        let (mut stacked, mut largest, mut peak) = (0u128, 0u128, 0u128);
        for pos in &self.pattern {
            let records = if pos.may_stay() {
                if pos.flags & LAZY_MATCH != 0 {
                    return true;
                }
                u128::from(pos.max_count) + 1
            } else {
                u128::from(pos.max_count.min(tlen as u64)) + 1
            };
            let growth = if records <= 1 {
                records
            } else {
                3 * (records.next_power_of_two() / 2)
            };
            peak = peak.max(growth + stacked);
            largest = largest.max(records);
            if pos.max_count > pos.min_count {
                stacked += records;
            }
        }
        peak = peak.max(largest.next_power_of_two() + stacked);
        peak * u128::from(RECORD_BYTES) + u128::from(text_heap) >= u128::from(HOST_DEPENDENT_BYTES)
    }

    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:700-742
    // The map's letters are not used (basedata.cpp:696 takes its length; the residues are
    // pushed into the same string after `clear()`, basedata.cpp:723-728): this follows the
    // string's length and capacity through `reserve` and the `append`s, which end NCBI's
    // rpsbproc with an exception (`MapEnd`) for a length or an allocation that cannot be held.
    pub(crate) fn minimal_map(&self) -> Result<CxxString, MapEnd> {
        let mut map = CxxString::new();
        let total = self.pattern.len() as u64;
        if total > 0 {
            map.reserve(total + total)?;
            for pos in &self.pattern {
                if pos.min_count > 0 && !pos.alt_chars.is_empty() {
                    map.append(pos.min_count)?;
                }
            }
        }
        Ok(map)
    }
}

#[cfg(test)]
#[path = "prosite_tests.rs"]
mod tests;
