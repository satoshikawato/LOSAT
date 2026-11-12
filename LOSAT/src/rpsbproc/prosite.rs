//! `CProSite` of rpsbproc's library: the parser and matcher of the PROSITE-like patterns of
//! `cddannot.dat` / `cddannot_generic.dat` (`[KR]`, `x(2)`, `{P}`, `<` and `>` anchors, `-`
//! separators), as NCBI's `Parse`, `Match` and `GetMinimalXMap` read and run them, including the
//! parser's error returns and the quirks of its state machine (a pattern that ends inside a
//! bracket or a repeat is accepted; an escape stays set across `(` and `)`).
//!
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.hpp:7-76
//! ```c++
//! class CProSite
//! {
//! public:
//! 	enum EParseError
//! 	{
//! 		eNoError,
//! 		eInvalidChar,
//! 		eSyntaxError,
//! 		eUnexpectedEnd
//! 	};
//! 	...
//! 	EParseError Parse(const std::string &expr, size_t &errorPos);
//! 	size_t Match(const std::string &text, size_t tlen, size_t start_pos) const;
//! 	// -- return a vector of minimal match length, with 'X' as generic, 'A' as alternative and S as strict.
//! 	void GetMinimalXMap(std::string &minMap) const;
//! private:
//! 	struct TPatternPos
//! 	{
//! 		static const unsigned int NEGATIVE_FILTER = 0x1;
//! 		static const unsigned int LAZY_MATCH = 0x2;
//! 		unsigned int m_uiFlags;	//positive or negative select
//! 		std::string m_strAltChars;
//! 		size_t m_ulMinCount, m_ulMaxCount;
//! 		size_t (TPatternPos::*m_pfnMatch)(const std::string &text, size_t tlen, size_t pos) const;
//! ```
//!
//! NCBI's `Match` reads `text[pos]` for every `pos < tlen`, and rpsbproc passes the length of the
//! whole query as `tlen` while the text holds only the mapped residues of a site (basedata.cpp:731).
//! A read at `text.size()` is the string's terminating 0 (defined); a read beyond it is
//! undefined in C++ (heap bytes). `Matcher` reports such a read (`ReadPast`) instead of guessing.

/// `CProSite::EParseError` without `eNoError`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ParseError {
    InvalidChar,
    SyntaxError,
    UnexpectedEnd,
}

const NEGATIVE_FILTER: u32 = 0x1;
const LAZY_MATCH: u32 = 0x2;
/// `(size_t)-1`: the biggest possible count.
const INFINITE: usize = usize::MAX;

/// `CProSite::TPatternPos`.
#[derive(Clone, Debug)]
struct PatternPos {
    flags: u32,
    alt_chars: Vec<u8>,
    min_count: usize,
    max_count: usize,
}

impl PatternPos {
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:11-13
    // ```c++
    // CProSite::TPatternPos::TPatternPos(void):
    // 	m_uiFlags(0), m_strAltChars(kPlainEmptyString), m_ulMinCount(1), m_ulMaxCount(1), m_pfnMatch(&TPatternPos::x_PosMatch)
    // {};
    // ```
    fn new() -> Self {
        PatternPos {
            flags: 0,
            alt_chars: Vec::new(),
            min_count: 1,
            max_count: 1,
        }
    }

    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:501-553
    // ```c++
    // void CProSite::TPatternPos::x_NormalizeAltChars(void)
    // {
    // 	if (!m_strAltChars.empty())
    // 	{
    // 		size_t idx = 0, ttl = m_strAltChars.size();
    // 		vector<size_t> veciidx;
    // 		veciidx.reserve(ttl);
    // 		while (idx < ttl)
    // 		{
    // 			char currc = m_strAltChars[idx];
    // 			if ('X' == currc)	//no matter others, just being the one
    // 			{
    // 				m_strAltChars.clear();
    // 				m_strAltChars.push_back('X');
    // 				break;
    // 			}
    // 			veciidx.clear();
    //
    // 			size_t rmv = 0;
    // 			size_t idx2 = idx + 1;
    // 			while (idx2 < ttl)
    // 			{
    // 				if (m_strAltChars[idx2] == currc)
    // 					++rmv;
    // 				else if (rmv > 0)
    // 					m_strAltChars[idx2 - rmv] = m_strAltChars[idx2];
    //
    // 				++idx2;
    // 			}
    //
    // 			ttl -= rmv;
    // 			m_strAltChars.erase(ttl);	//erase m
    // 			++idx;
    // 		}
    // 	}
    // }
    // ```
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

/// A read of `text` beyond its terminating 0 (see the module comment).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ReadPast;

/// The text of a match with NCBI's `tlen` (the length that bounds the positions) and the read of
/// `text[pos]`.
struct Matcher<'a> {
    text: &'a [u8],
    tlen: usize,
    read_past: std::cell::Cell<bool>,
}

impl Matcher<'_> {
    /// `text[pos]` of a `std::string`: the terminating 0 at `size()`.
    fn at(&self, pos: usize) -> u8 {
        match pos.cmp(&self.text.len()) {
            std::cmp::Ordering::Less => self.text[pos],
            std::cmp::Ordering::Equal => 0,
            std::cmp::Ordering::Greater => {
                self.read_past.set(true);
                0
            }
        }
    }
}

impl PatternPos {
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:554-592
    // ```c++
    // // -- negative position match
    // size_t CProSite::TPatternPos::x_NegMatch(const string &text, size_t tlen, size_t pos) const
    // {
    // 	string::const_iterator iterChar = m_strAltChars.begin(), iterCharEnd = m_strAltChars.end();
    // 	size_t adv = 0;
    //
    // 	while (iterCharEnd != iterChar)
    // 	{
    // 		switch (*iterChar)
    // 		{
    // 		case '<':	//start anchor
    // 			if (0 == pos)	//matched, fail it
    // 				return string::npos;
    // 			break;
    // 		case '>':	//end anchor
    // 			if (pos == tlen)	//matched
    // 				return string::npos;
    // 			break;
    // 		case 'X':
    // 			if (pos < tlen)
    // 				return string::npos;
    // 			break;
    // 		default:
    // 			if (pos < tlen)
    // 			{
    // 				char c = text[pos];
    // 				adv = 1;
    // 				if ('a' <= c && 'z' >= c) c -= 0x20;
    // 				if (*iterChar == c)	//matched
    // 					return string::npos;
    // 			}
    // 			else
    // 				return string::npos;
    // 		}
    //
    // 		++iterChar;
    // 	}
    // 	return adv;	//nothing matched, successful
    // }
    // ```
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
                        let mut c = m.at(pos);
                        adv = 1;
                        if c.is_ascii_lowercase() {
                            c -= 0x20;
                        }
                        if ch == c {
                            return None;
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
    // ```c++
    // // -- positive position match
    // size_t CProSite::TPatternPos::x_PosMatch(const string &text, size_t tlen, size_t pos) const
    // {
    // 	string::const_iterator iterChar = m_strAltChars.begin(), iterCharEnd = m_strAltChars.end();
    // 	while (iterCharEnd != iterChar)
    // 	{
    // 		switch (*iterChar)
    // 		{
    // 		case '<':	//start anchor
    // 			if (0 == pos)	//matched
    // 				return 0;
    // 			break;
    // 		case '>':	//end anchor
    // 			if (pos == tlen)	//matched
    // 				return 0;
    // 			break;
    // 		case 'X':	//match any
    // 			if (pos < tlen)	//last,
    // 				return 1;
    // 			break;
    // 		default:
    // 			if (pos < tlen)
    // 			{
    // 				char c = text[pos];
    // 				if ('a' <= c && 'z' >= c) c -= 0x20;
    // 				if (*iterChar == c)	//matched
    // 					return 1;
    // 			}
    // 		}
    //
    // 		++iterChar;
    // 	}
    // 	return string::npos;
    // }
    // ```
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
                        let mut c = m.at(pos);
                        if c.is_ascii_lowercase() {
                            c -= 0x20;
                        }
                        if ch == c {
                            return Some(1);
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
    // ```c++
    // int CProSite::TPatternPos::Match(const string &text, size_t tlen, vector<size_t> &rMatchRec, size_t start_pos) const
    // {
    // 	size_t matched = rMatchRec.size();	//assume: at least one [0] as the start position
    //
    // 	int retVal = 0;	//return 0: match success but no more flex, -1: match failed. >1: # of flex left
    // 	if (m_uiFlags & LAZY_MATCH)	//do not clear rMatchRec
    // 	{
    // 		if (0 == matched)	//the first time
    // 		{
    // 			rMatchRec.push_back(start_pos);
    // 			// -- first, must match to minimal count
    // 			while (rMatchRec.size() <= m_ulMinCount)
    // 			{
    // 				size_t adv = (this->*m_pfnMatch)(text, tlen, start_pos);
    // 				if (string::npos == adv)	//fail
    // 				{
    // 					return -1;
    // 				}
    // 				start_pos += adv;
    // 				rMatchRec.push_back(start_pos);
    // 			}
    // 			if (m_ulMaxCount > m_ulMinCount) retVal = 1;
    // 			else retVal = 0;
    // 		}
    // 		else if (matched <= m_ulMaxCount)	//still flex	//revisit
    // 		{
    // 			start_pos = rMatchRec[matched - 1];
    // 			size_t adv = (this->*m_pfnMatch)(text, tlen, start_pos);
    // 			if (string::npos == adv)	//successful
    // 				retVal = -1;
    // 			else
    // 			{
    // 				start_pos += adv;
    // 				rMatchRec.push_back(start_pos);
    // 				if (m_ulMaxCount > matched) retVal = 1;
    // 				else retVal = 0;
    // 			}
    // 		}
    // 	}
    // 	else	//aggressive match
    // 	{
    // 		if (0 == matched)	//the first time
    // 		{
    // 			rMatchRec.push_back(start_pos);
    // 			while ((matched = rMatchRec.size()) <= m_ulMaxCount)
    // 			{
    // 				size_t adv = (this->*m_pfnMatch)(text, tlen, start_pos);
    // 				if (string::npos == adv)	//fail
    // 					break;
    // 				start_pos += adv;
    // 				rMatchRec.push_back(start_pos);
    // 			}
    // 		}
    // 		else
    // 		{
    // 			rMatchRec.pop_back();
    // 			matched = rMatchRec.size();
    // 		}
    // 		if (matched < m_ulMinCount + 1) retVal = -1;	// minimal match not reached
    // 		else if (matched > m_ulMinCount + 1) retVal = 1;	//with flex
    // 		else retVal = 0;
    // 	}
    // 	return retVal;
    // }
    // ```
    // (A lazy position whose `matched` is above `m_ulMaxCount` on a revisit leaves `retVal` 0
    // and the record as it is.)
    fn match_at(&self, m: &Matcher, record: &mut Vec<usize>, mut start_pos: usize) -> i32 {
        let mut matched = record.len();
        let ret;
        if self.flags & LAZY_MATCH != 0 {
            if matched == 0 {
                record.push(start_pos);
                while record.len() <= self.min_count {
                    let Some(adv) = self.step(m, start_pos) else {
                        return -1;
                    };
                    start_pos = start_pos.wrapping_add(adv);
                    record.push(start_pos);
                }
                ret = i32::from(self.max_count > self.min_count);
            } else if matched <= self.max_count {
                start_pos = record[matched - 1];
                match self.step(m, start_pos) {
                    None => ret = -1,
                    Some(adv) => {
                        start_pos = start_pos.wrapping_add(adv);
                        record.push(start_pos);
                        ret = i32::from(self.max_count > matched);
                    }
                }
            } else {
                ret = 0;
            }
        } else {
            if matched == 0 {
                record.push(start_pos);
                loop {
                    matched = record.len();
                    if matched > self.max_count {
                        break;
                    }
                    let Some(adv) = self.step(m, start_pos) else {
                        break;
                    };
                    start_pos = start_pos.wrapping_add(adv);
                    record.push(start_pos);
                }
            } else {
                record.pop();
                matched = record.len();
            }
            let need = self.min_count.wrapping_add(1);
            ret = if matched < need {
                -1
            } else if matched > need {
                1
            } else {
                0
            };
        }
        ret
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

/// `(size_t)atol(text)` of digits: glibc's `strtol` saturates at `LONG_MAX`.
fn atol(digits: &[u8]) -> usize {
    let mut value: i64 = 0;
    for &digit in digits {
        value = value
            .saturating_mul(10)
            .saturating_add(i64::from(digit - b'0'));
    }
    value as usize
}

impl ProSite {
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:49-466
    // ```c++
    // CProSite::EParseError CProSite::Parse(const string &expr, size_t &errorPos)
    // {
    // 	if (expr.empty()) return eUnexpectedEnd;
    // 	...
    // 	while (iterCharEnd != iterChar)
    // 	{
    // 		char c = *iterChar;
    // 		if ('a' <= c && 'z' >= c) c -= 0x20;	//turn to capital letters
    // 		switch (c)
    // 		{
    // 		case '*':	//stop codon
    // 			if (!bEscaped)	//not escaped, as equivalent to (0, inf)
    // 			{
    // 				switch (status)
    // 				{
    // 				case ePlainText:
    // 					iterCurr->m_ulMinCount = 0;
    // 					iterCurr->m_ulMaxCount = -1;	//the biggest possible
    // 					if (iterCurr->m_strAltChars == "X") iterCurr->m_uiFlags |= TPatternPos::LAZY_MATCH;
    // 					break;
    // 				case eAlternates:	//always treat as a literal * (stop codon)
    // 				case eNegations:
    // 					iterCurr->m_strAltChars.push_back(c);
    // 					break;
    // 				case eRange:	//accepting max
    // 				case eRepeats:
    // 					if (lit.empty() && !bInfLit)
    // 					{
    // 						bInfLit = true;
    // 						break;
    // 					}
    // 				case eAltRedund:
    // 				case eNegRedund:
    // 					break;
    // 				default:
    // 					eErrorType = eSyntaxError;
    // 					goto labelErrorRet;
    // 				}
    // 				break;
    // 			}
    //
    // 		case 'A': case 'C': ... case 'Y':	// the 20 amino acids and 'U'
    // 			switch (status)
    // 			{
    // 			case eParserReady:	//No current
    // 				status = ePlainText;
    // 			case ePlainText:	//current is valid
    // 				iterCurr = m_vecCompiledPattern.emplace(m_vecCompiledPattern.end(), __dummy);
    // 			case eAlternates:	//single char alternative, assume: iterCurrAltPattern is undefined
    // 			case eNegations:
    // 				iterCurr->m_strAltChars.push_back(c);
    // 				break;
    // 			default:
    // 				eErrorType = eSyntaxError;
    // 				goto labelErrorRet;
    // 			}
    // 			bEscaped = false;
    // 			break;
    // 		case 'X': ...
    // 		case '0': ... case '9': ...
    // 		case '(': ... case ')': ... case ',': ... case '\\': ... case '[': ... case ']': ...
    // 		case '{': ... case '}': ... case '<': ... case '>': ...
    // 		case '-':	//ignore
    // 		case ' ':
    // 		case '\t':
    // 			break;
    // 		default:
    // 			eErrorType = eInvalidChar;
    // 			goto labelErrorRet;
    // 		}
    // 		++iterChar;
    // 	}
    //
    // labelErrorRet:
    // 	errorPos = (iterChar - expr.begin());
    // 	return eErrorType;
    // }
    // ```
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
        // `emplace(end, __dummy)` and `iterCurr = ...`.
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
    // ```c++
    // size_t CProSite::Match(const string &text, size_t tlen, size_t start_pos) const
    // {
    // 	if (start_pos >= tlen) return string::npos;
    // 	stack<x_TMatchRec> stkFlexStack;
    //
    // 	x_TMatchRec curr(m_vecCompiledPattern.begin());
    //
    // 	vector<TPatternPos> :: const_iterator iterPosEnd = m_vecCompiledPattern.end();
    //
    // 	while (iterPosEnd != curr.iterPos)
    // 	{
    // 		int matched = curr.iterPos->Match(text, tlen, curr.flex, start_pos);
    // 		if (matched < 0)	//match failed
    // 		{
    // 			if (stkFlexStack.empty()) return string::npos;
    // 			curr = stkFlexStack.top();
    // 			stkFlexStack.pop();
    // 		}
    // 		else	//match success
    // 		{
    // 			if (matched > 0) stkFlexStack.push(curr);
    // 			start_pos = curr.flex.back();
    // 			++curr.iterPos;
    // 			curr.flex.clear();
    // 		}
    // 	}
    //
    // 	// -- successfully matched all PatternPosition
    // 	return start_pos;	//should always be valid.
    // }
    // ```
    // `Ok(None)` is `string::npos`. A pattern without a position (only `-`) matches at
    // `start_pos`. `Err(ReadPast)`: the match read `text` beyond its terminating 0.
    pub(crate) fn match_at(
        &self,
        text: &[u8],
        tlen: usize,
        mut start_pos: usize,
    ) -> Result<Option<usize>, ReadPast> {
        if start_pos >= tlen {
            return Ok(None);
        }
        let matcher = Matcher {
            text,
            tlen,
            read_past: std::cell::Cell::new(false),
        };
        let mut stack: Vec<(usize, Vec<usize>)> = Vec::new();
        let mut pos = 0usize;
        let mut flex: Vec<usize> = Vec::new();
        let result = loop {
            if pos == self.pattern.len() {
                break Some(start_pos);
            }
            let matched = self.pattern[pos].match_at(&matcher, &mut flex, start_pos);
            if matcher.read_past.get() {
                return Err(ReadPast);
            }
            if matched < 0 {
                let Some((top_pos, top_flex)) = stack.pop() else {
                    break None;
                };
                pos = top_pos;
                flex = top_flex;
            } else {
                if matched > 0 {
                    stack.push((pos, flex.clone()));
                }
                start_pos = *flex.last().expect("a matched position records its start");
                pos += 1;
                flex.clear();
            }
        };
        Ok(result)
    }

    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/prosite.cpp:700-742
    // ```c++
    // void CProSite::GetMinimalXMap(string &minMap) const
    // {
    // 	minMap.clear();
    // 	size_t ttl = m_vecCompiledPattern.size();
    // 	if (ttl > 0)
    // 	{
    // 		minMap.reserve(ttl + ttl);
    // 		for (size_t i = 0; i < ttl; ++i)
    // 		{
    // 			if (m_vecCompiledPattern[i].m_ulMinCount > 0)
    // 			{
    // 				switch (m_vecCompiledPattern[i].m_strAltChars.size())
    // 				{
    // 				case 0:
    // 					break;
    // 				case 1:
    // 					if ('X' == m_vecCompiledPattern[i].m_strAltChars[0]) minMap.append(m_vecCompiledPattern[i].m_ulMinCount, 'X');
    // 					else if (m_vecCompiledPattern[i].m_uiFlags & TPatternPos::NEGATIVE_FILTER) minMap.append(m_vecCompiledPattern[i].m_ulMinCount, 'A');
    // 					else minMap.append(m_vecCompiledPattern[i].m_ulMinCount, 'S');
    // 					break;
    // 				default:	//>1
    // 					minMap.append(m_vecCompiledPattern[i].m_ulMinCount, 'A');
    // 				}
    // 			}
    // 		}
    // 	}
    // }
    // ```
    // The only use of the map is its length (basedata.cpp:696), which is all this returns: the
    // sum of the minimal counts of the positions that have alternatives.
    pub(crate) fn minimal_map_len(&self) -> usize {
        self.pattern
            .iter()
            .filter(|pos| pos.min_count > 0 && !pos.alt_chars.is_empty())
            .fold(0usize, |sum, pos| sum.saturating_add(pos.min_count))
    }
}

#[cfg(test)]
#[path = "prosite_tests.rs"]
mod tests;
