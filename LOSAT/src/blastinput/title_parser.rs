//! `CTitleParser::x_FindBrackets` and `CTitleParser::HasMods`: the one place that finds the
//! modifiers (`[name=value]`) of a FASTA title, for the FASTA readers (`HasMods`: whether a title
//! has any) and the alignment reader (`CTitleParser::Apply`, `aln_reader::mods`).

use crate::utils::byte_search::find_first;

/// `CTitleParser::x_FindBrackets`: from `start`, the first `[` that is not inside brackets
/// already open (`start` is moved to it), and its closing `]` (`stop`), with the first `=` inside.
/// `None` where there is no closing bracket, a `]` comes before any `[`, or no `=` comes before
/// the closing bracket. Returns `(eq_pos, stop)`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/mod_reader.cpp:788-824
pub(crate) fn find_brackets(line: &[u8], start: &mut usize) -> Option<(usize, usize)> {
    let mut eq_pos = None;
    let mut open = 0usize;
    let mut i = *start;
    while i < line.len() {
        match line[i] {
            b'[' => {
                open += 1;
                if open == 1 {
                    *start = i;
                }
            }
            b'=' => {
                if open > 0 && eq_pos.is_none() {
                    eq_pos = Some(i);
                }
            }
            b']' => {
                if open == 1 {
                    return eq_pos.filter(|&eq| eq < i).map(|eq| (eq, i));
                } else if open == 0 {
                    return None;
                }
                open -= 1;
            }
            _ => {}
        }
        i += 1;
    }
    None
}

/// `CTitleParser::HasMods` of `title`: whether the title has a modifier, `[name=value]`. NCBI's
/// loop returns at the first bracket pair that `find_brackets` finds (true with an `=` before its
/// `]`, false without one) and false where it finds none, so it is true exactly where the first
/// call succeeds.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/mod_reader.cpp:793-850
/// Before the first `[` or `]` of the title `find_brackets` changes nothing (an `=` counts only
/// inside brackets), so it starts there, found with one search (LOSAT's speed).
pub(crate) fn has_mods(title: &[u8]) -> bool {
    find_first(title, |byte| (byte == b'[') | (byte == b']'))
        .is_some_and(|mut start| find_brackets(title, &mut start).is_some())
}

#[cfg(test)]
#[path = "title_parser_tests.rs"]
mod tests;
