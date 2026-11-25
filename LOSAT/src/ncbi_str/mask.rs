//! `NStr::MatchesMask` with `NStr::eCase`: a file name or a word against a shell-like mask, as
//! NCBI's directory listings (`CDir::GetEntries` with a `CMaskFileName`) and the accession guide's
//! wildcards match them.

/// What the match of a rest of the text against a rest of the mask found (`EMatchesMaskResult`):
/// a match, no match at this place (a `*` before it may try the next), or a mismatch that no
/// place can mend (the text ended, or the mask is malformed).
#[derive(Clone, Copy, PartialEq, Eq)]
enum Found {
    Match,
    NoMatch,
    Mismatch,
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:291-392
// `CTempString::operator[]` gives NUL past the end, so a NUL byte ends the text or the mask as
// its end does. `?` is any one byte; `*` any run (several in a row are one); `[...]` a set of
// bytes and ranges `a-b`, `[!...]` its complement (a `-` before `]` is a byte); `\` makes the
// next byte plain.
/// Whether `text` matches `mask`, byte for byte (case-sensitive).
pub(crate) fn matches_mask(text: &[u8], mask: &[u8]) -> bool {
    found(text, mask) == Found::Match
}

/// The byte at `at` of `bytes`, NUL past the end (`CTempString::operator[]`).
fn at(bytes: &[u8], at: usize) -> u8 {
    bytes.get(at).copied().unwrap_or(0)
}

fn found(text: &[u8], mask: &[u8]) -> Found {
    let (mut t, mut m) = (0, 0);
    loop {
        let mut mc = at(mask, m);
        if mc == 0 {
            break;
        }
        let s = at(text, t);
        if s == 0 && mc != b'*' {
            return Found::Mismatch;
        }
        match mc {
            b'?' => {}
            b'*' => {
                while at(mask, m) == b'*' {
                    m += 1;
                }
                if at(mask, m) == 0 {
                    return Found::Match;
                }
                let mut s = s;
                while s != 0 {
                    let rest = found(text.get(t..).unwrap_or(&[]), mask.get(m..).unwrap_or(&[]));
                    if rest != Found::NoMatch {
                        return rest;
                    }
                    s = at(text, t);
                    t += 1;
                }
                return Found::Mismatch;
            }
            b'[' => {
                m += 1;
                if at(mask, m) == 0 {
                    return Found::Mismatch;
                }
                let complement = at(mask, m) == b'!';
                if complement {
                    m += 1;
                }
                let mut seen = false;
                loop {
                    let a = at(mask, m);
                    m += 1;
                    if a == 0 {
                        return Found::Mismatch;
                    }
                    let b = if at(mask, m) == b'-' && at(mask, m + 1) != b']' {
                        m += 1;
                        let b = at(mask, m);
                        m += 1;
                        if b == 0 {
                            return Found::Mismatch;
                        }
                        b
                    } else {
                        a
                    };
                    // NCBI compares `char`s, signed on the oracle's platform.
                    if !seen && (a as i8) <= (s as i8) && (s as i8) <= (b as i8) {
                        seen = true;
                    }
                    if at(mask, m) == b']' {
                        break;
                    }
                }
                if complement == seen {
                    return Found::NoMatch;
                }
            }
            _ => {
                if mc == b'\\' {
                    m += 1;
                    mc = at(mask, m);
                    if mc == 0 {
                        return Found::Mismatch;
                    }
                }
                if s != mc {
                    return Found::NoMatch;
                }
            }
        }
        t += 1;
        m += 1;
    }
    if at(text, t) != 0 {
        return Found::NoMatch;
    }
    Found::Match
}

#[cfg(test)]
#[path = "mask_tests.rs"]
mod tests;
