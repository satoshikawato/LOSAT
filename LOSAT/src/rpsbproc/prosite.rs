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
//! memory runs out. `Record` keeps the counts instead of the records, and one allocation of
//! NCBI's of 1 GiB or more for them is refused (`MatchError::Exhausted`; NCBI's result depends on
//! the host's memory). `ProSite::minimal_map` follows `GetMinimalXMap`'s string by its length and
//! capacity, with the exceptions that its `append` throws (`MapEnd`).

/// One allocation of NCBI's of this many bytes or more is refused: below 2^47 bytes the
/// allocation succeeds or throws `std::bad_alloc` by the host's memory and overcommit
/// (DECISIONS 37, Owner-delegated 2026-11-12: NCBI's rpsbproc uses about 130 MB, and a 64-bit
/// host that processes CDD gives it 1 GiB).
pub(crate) const HOST_DEPENDENT_BYTES: u64 = 1 << 30;

/// The user address space of x86-64 (4-level paging): `operator new` of this many bytes or more
/// fails on every host and throws `std::bad_alloc`.
pub(crate) const USER_SPACE_BYTES: u64 = 1 << 47;

/// How NCBI's `GetMinimalXMap` ends instead of returning.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum MapEnd {
    /// An exception of libstdc++, by its `what()`: `LENGTH_ERROR` or `BAD_ALLOC`.
    Throw(&'static str),
    /// An allocation of this many bytes, from `HOST_DEPENDENT_BYTES` up to `USER_SPACE_BYTES`.
    HostDependent(u64),
}

/// `what()` of the `std::length_error` of `_M_check_length` in `append` (the length passes
/// `max_size()`).
pub(crate) const LENGTH_ERROR: &str = "basic_string::_M_replace_aux";

/// `what()` of `std::bad_alloc` (an allocation of `USER_SPACE_BYTES` or more).
pub(crate) const BAD_ALLOC: &str = "std::bad_alloc";

/// `operator new(bytes)` of NCBI's for a string.
fn allocate_string(bytes: u64) -> Result<(), MapEnd> {
    if bytes >= USER_SPACE_BYTES {
        Err(MapEnd::Throw(BAD_ALLOC))
    } else if bytes >= HOST_DEPENDENT_BYTES {
        Err(MapEnd::HostDependent(bytes))
    } else {
        Ok(())
    }
}

/// The length and capacity of a `std::string` of libstdc++ (the C++11 ABI) whose bytes LOSAT
/// does not need: `GetMinimalXMap`'s `minMap` and the text that `MotifCheck` builds in it.
/// `max_size()` is `(PTRDIFF_MAX - 1) / 2` = 2^62 - 1 (GCC 9 and later; NCBI's rpsbproc 0.5.1
/// throws `length_error` for an `append` of `LONG_MAX` characters, review oracle R20).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct CxxString {
    pub(crate) len: u64,
    /// `capacity()`: the bytes `0..=capacity` are NCBI's buffer (`_M_local_buf[16]` while it
    /// is 15).
    pub(crate) capacity: u64,
}

impl CxxString {
    const LOCAL_CAPACITY: u64 = 15;
    const MAX_SIZE: u64 = (1 << 62) - 1;

    /// An empty string (`string minMap(k_strEmptyString)`).
    pub(crate) fn new() -> Self {
        CxxString {
            len: 0,
            capacity: Self::LOCAL_CAPACITY,
        }
    }

    // libstdc++ (GCC 13): bits/basic_string.tcc:139-160
    // ```c++
    //   _M_create(size_type& __capacity, size_type __old_capacity)
    //   {
    //     // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //     // 83.  String::npos vs. string::max_size()
    //     if (__capacity > max_size())
    // 	std::__throw_length_error(__N("basic_string::_M_create"));
    //
    //     // The below implements an exponential growth policy, necessary to
    //     // meet amortized linear time requirements of the library: see
    //     // http://gcc.gnu.org/ml/libstdc++/2001-07/msg00085.html.
    //     if (__capacity > __old_capacity && __capacity < 2 * __old_capacity)
    // 	{
    // 	  __capacity = 2 * __old_capacity;
    // 	  // Never allocate a string bigger than max_size.
    // 	  if (__capacity > max_size())
    // 	    __capacity = max_size();
    // 	}
    //
    //     // NB: Need an array of char_type[__capacity], plus a terminating
    //     // null char_type() element.
    //     return _S_allocate(_M_get_allocator(), __capacity + 1);
    //   }
    // ```
    // (The callers never ask for more than `max_size()`: `append` checks the length first.)
    fn create(&mut self, mut capacity: u64) -> Result<(), MapEnd> {
        if capacity > self.capacity && capacity < 2 * self.capacity {
            capacity = (2 * self.capacity).min(Self::MAX_SIZE);
        }
        allocate_string(capacity + 1)?;
        self.capacity = capacity;
        Ok(())
    }

    // libstdc++ (GCC 13): bits/basic_string.tcc:305-320
    // ```c++
    //     reserve(size_type __res)
    //     {
    //       const size_type __capacity = capacity();
    //       // _GLIBCXX_RESOLVE_LIB_DEFECTS
    //       // 2968. Inconsistencies between basic_string reserve and
    //       // vector/unordered_map/unordered_set reserve functions
    //       // P0966 reserve should not shrink
    //       if (__res <= __capacity)
    // 	return;
    //
    //       pointer __tmp = _M_create(__res, __capacity);
    //       this->_S_copy(__tmp, _M_data(), length() + 1);
    //       _M_dispose();
    //       _M_data(__tmp);
    //       _M_capacity(__res);
    //     }
    // ```
    // (GCC before 11 also shrinks; `minMap` is empty and local when `GetMinimalXMap` reserves,
    // so both give this.)
    fn reserve(&mut self, res: u64) -> Result<(), MapEnd> {
        if res <= self.capacity {
            return Ok(());
        }
        self.create(res)
    }

    // libstdc++ (GCC 13): bits/basic_string.tcc:451-468
    // ```c++
    //     _M_replace_aux(size_type __pos1, size_type __n1, size_type __n2,
    // 		   _CharT __c)
    //     {
    //       _M_check_length(__n1, __n2, "basic_string::_M_replace_aux");
    //
    //       const size_type __old_size = this->size();
    //       const size_type __new_size = __old_size + __n2 - __n1;
    //
    //       if (__new_size <= this->capacity())
    // 	{
    // 	  ...
    // 	}
    //       else
    // 	this->_M_mutate(__pos1, __n1, 0, __n2);
    // ```
    // libstdc++ (GCC 13): bits/basic_string.h:397-401,1084-1085,1487-1488
    // ```c++
    //       _M_check_length(size_type __n1, size_type __n2, const char* __s) const
    //       {
    // 	if (this->max_size() - (this->size() - __n1) < __n2)
    // 	  __throw_length_error(__N(__s));
    //       }
    //       ...
    //       max_size() const _GLIBCXX_NOEXCEPT
    //       { return (_Alloc_traits::max_size(_M_get_allocator()) - 1) / 2; }
    //       ...
    //       append(size_type __n, _CharT __c)
    //       { return _M_replace_aux(this->size(), size_type(0), __n, __c); }
    // ```
    // libstdc++ (GCC 13): bits/basic_string.tcc:326-332
    // ```c++
    //     _M_mutate(size_type __pos, size_type __len1, const _CharT* __s,
    // 	      size_type __len2)
    //     {
    //       const size_type __how_much = length() - __pos - __len1;
    //
    //       size_type __new_capacity = length() + __len2 - __len1;
    //       pointer __r = _M_create(__new_capacity, capacity());
    // ```
    /// `append(n, c)`.
    pub(crate) fn append(&mut self, n: u64) -> Result<(), MapEnd> {
        if Self::MAX_SIZE - self.len < n {
            return Err(MapEnd::Throw(LENGTH_ERROR));
        }
        let len = self.len + n;
        if len > self.capacity {
            self.create(len)?;
        }
        self.len = len;
        Ok(())
    }

    // libstdc++ (GCC 13): bits/basic_string.h:1561-1568
    // ```c++
    //       push_back(_CharT __c)
    //       {
    // 	const size_type __size = this->size();
    // 	if (__size + 1 > this->capacity())
    // 	  this->_M_mutate(__size, size_type(0), 0, size_type(1));
    // 	traits_type::assign(this->_M_data()[__size], __c);
    // 	this->_M_set_length(__size + 1);
    //       }
    // ```
    pub(crate) fn push_back(&mut self) -> Result<(), MapEnd> {
        if self.len + 1 > self.capacity {
            self.create(self.len + 1)?;
        }
        self.len += 1;
        Ok(())
    }

    /// `clear()`: the capacity stays.
    pub(crate) fn clear(&mut self) {
        self.len = 0;
    }
}

/// `CProSite::EParseError` without `eNoError`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ParseError {
    InvalidChar,
    SyntaxError,
    UnexpectedEnd,
}

const NEGATIVE_FILTER: u32 = 0x1;
const LAZY_MATCH: u32 = 0x2;
/// `(size_t)-1`: the biggest possible count. The counts are NCBI's `size_t` of x86-64, `u64`
/// on every target.
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

/// Why `ProSite::matches` or `ProSite::match_at` gives no answer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum MatchError {
    /// NCBI's answer depends on the bytes it reads beyond the text's terminating 0, or the
    /// match may read beyond NCBI's buffer of the text (see the module comment).
    ReadPast,
    /// One allocation of NCBI's match records reaches `HOST_DEPENDENT_BYTES`.
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

/// The bytes of `sizeof(size_t)`, one match record.
const RECORD_BYTES: u64 = 8;

/// One allocation of `records` match records (`vector<size_t>`) of NCBI's.
fn allocate_records(records: u64) -> Result<(), MatchError> {
    if u128::from(records) * u128::from(RECORD_BYTES) >= u128::from(HOST_DEPENDENT_BYTES) {
        Err(MatchError::Exhausted)
    } else {
        Ok(())
    }
}

/// `push_back` of records into the vector of `capacity` until it holds `len`: libstdc++ makes
/// room for twice the records it holds when it is full (one when it is empty).
// libstdc++ (GCC 13): bits/stl_vector.h:1281-1293,1896-1903
// ```c++
//       push_back(const value_type& __x)
//       {
// 	if (this->_M_impl._M_finish != this->_M_impl._M_end_of_storage)
// 	  {
// 	    ...
// 	  }
// 	else
// 	  _M_realloc_insert(end(), __x);
//       }
//       ...
//       _M_check_len(size_type __n, const char* __s) const
//       {
// 	if (max_size() - size() < __n)
// 	  __throw_length_error(__N(__s));
//
// 	const size_type __len = size() + (std::max)(size(), __n);
// 	return (__len < size() || __len > max_size()) ? max_size() : __len;
//       }
// ```
// libstdc++ (GCC 13): bits/vector.tcc:446,453-455
// ```c++
//       _M_realloc_insert(iterator __position, _Args&&... __args)
//     ...
//     {
//       const size_type __len =
// 	_M_check_len(size_type(1), "vector::_M_realloc_insert");
// ```
// (`max_size()` is 2^60 - 1 records; the allocation of 1 GiB stops the growth long before.)
fn grow_records(capacity: &mut u64, len: u64) -> Result<(), MatchError> {
    while *capacity < len {
        *capacity = if *capacity == 0 { 1 } else { 2 * *capacity };
        allocate_records(*capacity)?;
    }
    Ok(())
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

    /// `rMatchRec.push_back(start_pos)` of the first visit.
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
    // and the record as it is.) `capacity` is the capacity of NCBI's `curr.flex`, the vector
    // that holds `rMatchRec`.
    //
    // A step that does not advance gives the same at every later step (it depends only on the
    // position): NCBI's loops repeat it up to the count they stop at, which `Record::stay_until`
    // counts at once. A revisit starts where the rest of the pattern has just failed: one that
    // does not advance (a lazy step, or the pop of a record of a step that did not advance)
    // leaves the rest to fail again from the same position, so NCBI revisits until the count
    // runs out (a lazy `{X}*` for ever) and then fails; LOSAT goes there at once.
    fn match_at(
        &self,
        m: &Matcher,
        record: &mut Record,
        capacity: &mut u64,
        mut start_pos: usize,
    ) -> Result<i32, MatchError> {
        let matched = record.len;
        let ret;
        if self.flags & LAZY_MATCH != 0 {
            if matched == 0 {
                record.begin(start_pos);
                grow_records(capacity, 1)?;
                while record.len <= self.min_count {
                    let Some(adv) = self.step(m, start_pos) else {
                        return Ok(-1);
                    };
                    if adv == 0 {
                        let len = self.min_count.saturating_add(1);
                        grow_records(capacity, len)?;
                        record.stay_until(len);
                        break;
                    }
                    start_pos += adv;
                    record.advance();
                    grow_records(capacity, record.len)?;
                }
                ret = i32::from(self.max_count > self.min_count);
            } else if matched <= self.max_count {
                start_pos = record.last();
                match self.step(m, start_pos) {
                    None => ret = -1,
                    Some(0) => {
                        grow_records(capacity, self.max_count.saturating_add(1))?;
                        ret = -1;
                    }
                    Some(adv) => {
                        debug_assert_eq!(adv, 1);
                        record.advance();
                        grow_records(capacity, record.len)?;
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
                grow_records(capacity, 1)?;
                while record.len <= self.max_count {
                    let Some(adv) = self.step(m, start_pos) else {
                        break;
                    };
                    if adv == 0 {
                        let len = self.max_count.saturating_add(1);
                        grow_records(capacity, len)?;
                        record.stay_until(len);
                        break;
                    }
                    start_pos += adv;
                    record.advance();
                    grow_records(capacity, record.len)?;
                }
            } else {
                if record.still > 0 {
                    // The pops down to the advanced records keep `rMatchRec.back()`; they end
                    // the position (`retVal` 0 or -1 at `need`) or reach the advanced records.
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
    fn may_stay(&self) -> bool {
        if self.negative() {
            !self.reads()
        } else {
            self.alt_chars.iter().any(|&c| c == b'<' || c == b'>')
        }
    }

    /// Whether a step reads `text` (an alternative that is not an anchor or `X`).
    fn reads(&self) -> bool {
        self.alt_chars
            .iter()
            .any(|&c| c != b'<' && c != b'>' && c != b'X')
    }

    /// Whether a step can advance: a positive position with an alternative that is not an
    /// anchor, or a negation that reads.
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

/// `(size_t)atol(text)` of digits: glibc's `strtol` saturates at `LONG_MAX`.
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
    // `start_pos`. The second value says whether the match read `text` beyond its terminating 0
    // (answered by `past`). `capacity` follows `curr.flex`: a copy onto the stack allocates the
    // records it holds, a copy back allocates when they do not fit (libstdc++'s `vector` copy
    // constructor and assignment, bits/stl_vector.h:601-609, bits/vector.tcc:235-246).
    fn run(
        &self,
        text: &[u8],
        tlen: usize,
        mut start_pos: usize,
        past: Past,
    ) -> (Result<Option<usize>, MatchError>, bool) {
        if start_pos >= tlen {
            return (Ok(None), false);
        }
        let matcher = Matcher {
            text,
            tlen,
            past,
            read_past: std::cell::Cell::new(false),
        };
        let mut stack: Vec<(usize, Record)> = Vec::new();
        let mut pos = 0usize;
        let mut flex = Record::default();
        let mut capacity = 0u64;
        let result = loop {
            if pos == self.pattern.len() {
                break Ok(Some(start_pos));
            }
            let matched =
                match self.pattern[pos].match_at(&matcher, &mut flex, &mut capacity, start_pos) {
                    Ok(matched) => matched,
                    Err(error) => break Err(error),
                };
            if matched < 0 {
                let Some((top_pos, top_flex)) = stack.pop() else {
                    break Ok(None);
                };
                if top_flex.len > capacity {
                    if let Err(error) = allocate_records(top_flex.len) {
                        break Err(error);
                    }
                    capacity = top_flex.len;
                }
                pos = top_pos;
                flex = top_flex;
            } else {
                if matched > 0 {
                    if let Err(error) = allocate_records(flex.len) {
                        break Err(error);
                    }
                    stack.push((pos, flex));
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
    /// comment; DECISIONS 37). `buffer` is the capacity of NCBI's string that holds the text: its
    /// bytes `0..=buffer` are NCBI's own, a read beyond them may fault.
    pub(crate) fn matches(
        &self,
        text: &[u8],
        tlen: usize,
        start_pos: usize,
        buffer: u64,
    ) -> Result<bool, MatchError> {
        let (result, read_past) = self.run(text, tlen, start_pos, Past::Fail);
        if !read_past {
            return result.map(|end| end.is_some());
        }
        // The runs that NCBI makes for other bytes go through other positions and counts:
        // every read must stay in NCBI's buffer, and no record may reach the allocation that
        // `Exhausted` refuses.
        if self.may_exhaust(tlen) || self.max_read(tlen, start_pos) > buffer {
            return Err(MatchError::ReadPast);
        }
        if let Ok(Some(_)) = result {
            return Ok(true);
        }
        match self.run(text, tlen, start_pos, Past::Pass).0 {
            Ok(None) => Ok(false),
            _ => Err(MatchError::ReadPast),
        }
    }

    /// The highest position that a step that reads can be at, for any bytes of the text: each
    /// step advances by one at most, each position up to its maximal count, a step that reads
    /// advances, and reads stop before `tlen`.
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

    /// Whether the records of a position can reach `2^26 + 1`, the fewest whose growth
    /// allocates 1 GiB (twice 2^26 records of 8 bytes): a position that may not advance
    /// repeats up to its maximal count (a lazy one for ever), one that advances up to `tlen`.
    fn may_exhaust(&self, tlen: usize) -> bool {
        const FEWEST: u64 = (1 << 26) + 1;
        self.pattern.iter().any(|pos| {
            let records = if pos.may_stay() {
                if pos.flags & LAZY_MATCH != 0 {
                    u64::MAX
                } else {
                    pos.max_count.saturating_add(1)
                }
            } else {
                pos.max_count.min(tlen as u64).saturating_add(1)
            };
            records >= FEWEST
        })
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
