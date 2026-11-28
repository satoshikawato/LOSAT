//! `BlastFilteringOptionsFromString`: the filtering options that a `FilterString` (the archive's
//! and the search strategy's `FilterString` parameter, `CBlastOptions::SetFilterString`) gives the
//! options handle, which its remote options then record (`options::RemoteOptions`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:40-294,435-605

use crate::ncbi_str::{scan_double, scan_long};

/// `BLASTOPTIONS_BUFFER_SIZE`: the bytes of one command that the parsers read.
const BUFFER_SIZE: usize = 128;

/// `kDefaultRepeatFilterDb` (blast_options.h), the database of `SRepeatFilterOptionsNew`.
const DEFAULT_REPEAT_FILTER_DB: &[u8] = b"repeat/repeat_9606";

/// `SDustOptionsNew`: level 20, window 64, linker 1.
const DEFAULT_DUST: [Option<i32>; 3] = [Some(20), Some(64), Some(1)];

/// `SSegOptionsNew`: window 12, locut 2.2, hicut 2.5.
const DEFAULT_SEG: Seg = Seg {
    window: Some(12),
    locut: Some(2.2),
    hicut: Some(2.5),
};

/// The SEG options of a `FilterString`; `None` is a value that `sscanf` did not store (NCBI then
/// reads a variable it did not set).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Seg {
    pub(crate) window: Option<i32>,
    pub(crate) locut: Option<f64>,
    pub(crate) hicut: Option<f64>,
}

/// The filtering options of a `FilterString` (`SBlastFilterOptions`), as far as an options
/// handle's remote options read them (the window masker's are not).
#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct FilterOptions {
    pub(crate) mask_at_hash: bool,
    /// The DUST level, window and linker; `None` is a value that `sscanf` did not store.
    pub(crate) dust: Option<[Option<i32>; 3]>,
    pub(crate) seg: Option<Seg>,
    /// The repeat filter's database.
    pub(crate) repeat_db: Option<Vec<u8>>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:435-605
/// `BlastFilteringOptionsFromString` of `instructions` for a program whose type is BLASTN's or the
/// mapper's (`blastn_or_mapping`: `L` and `T` are DUST, else SEG): `None` for NCBI's parse error
/// (the options are then those that `ClearFilterOptions` left).
pub(crate) fn filtering_options_from_string(
    instructions: &[u8],
    blastn_or_mapping: bool,
) -> Option<FilterOptions> {
    let mut options = FilterOptions::default();
    if instructions.eq_ignore_ascii_case(b"F") {
        return Some(options);
    }
    let mut rest = instructions;
    // Old-style filters, where `m` cannot be followed by the `;`.
    if let [b'm', b' ', after @ ..] = rest {
        options.mask_at_hash = true;
        rest = after;
    }
    while let [command, after @ ..] = rest {
        rest = after;
        match command {
            b'S' | b'D' | b'R' | b'W' => {
                let (buffer, next) = load_options_to_buffer(rest);
                rest = next;
                match command {
                    b'S' => {
                        options.seg = Some(if buffer.is_empty() {
                            DEFAULT_SEG
                        } else {
                            parse_seg_options(&buffer)?
                        })
                    }
                    b'D' => {
                        options.dust = Some(if buffer.is_empty() {
                            DEFAULT_DUST
                        } else {
                            parse_dust_options(&buffer)?
                        })
                    }
                    b'R' => {
                        options.repeat_db = Some(
                            parse_repeat_options(&buffer)
                                .unwrap_or_else(|| DEFAULT_REPEAT_FILTER_DB.to_vec()),
                        )
                    }
                    // The window masker's options, which no remote option records.
                    _ => {}
                }
            }
            b'L' | b'T' => {
                if blastn_or_mapping {
                    options.dust = Some(DEFAULT_DUST);
                } else {
                    options.seg = Some(DEFAULT_SEG);
                }
            }
            b'm' => options.mask_at_hash = true,
            _ => {}
        }
    }
    Some(options)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:54-97
/// `s_LoadOptionsToBuffer`: the bytes of one command, up to its `;` (read) or `BUFFER_SIZE`
/// bytes, without the blanks before and after; and the rest of `instructions`.
fn load_options_to_buffer(instructions: &[u8]) -> (Vec<u8>, &[u8]) {
    let mut buffer = Vec::new();
    let mut not_started = true;
    let mut read = 0;
    for _ in 0..BUFFER_SIZE {
        let Some(&byte) = instructions.get(read) else {
            break;
        };
        read += 1;
        if byte == b';' {
            break;
        }
        if not_started && byte == b' ' {
            continue;
        }
        not_started = false;
        buffer.push(byte);
    }
    while buffer.len() > 1 && buffer.last() == Some(&b' ') {
        buffer.pop();
    }
    (buffer, &instructions[read..])
}

/// The tokens that `s_ParseDustOptions` and `s_ParseSegOptions` read of `buffer` (a command
/// without blanks before and after): the bytes between spaces, each ended at a space or the end
/// within `BUFFER_SIZE` steps (one step a byte of a token or the end of one).
fn option_tokens(buffer: &[u8]) -> Vec<&[u8]> {
    let at = |pos: usize| buffer.get(pos).copied().unwrap_or(0);
    let mut tokens = Vec::new();
    let (mut start, mut pos) = (0, 0);
    for _ in 0..BUFFER_SIZE {
        if matches!(at(pos), b' ' | 0) {
            tokens.push(&buffer[start..pos]);
            while at(pos) == b' ' {
                pos += 1;
            }
            if at(pos) == 0 {
                break;
            }
            start = pos;
        } else {
            pos += 1;
        }
    }
    tokens
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:180-234
/// `s_ParseDustOptions`: the level, window and linker of three tokens (`sscanf` with `%ld`, cast
/// to `int`), -1 each without one; `None` for another count (NCBI's parse error).
fn parse_dust_options(buffer: &[u8]) -> Option<[Option<i32>; 3]> {
    match option_tokens(buffer)[..] {
        [] => Some([Some(-1); 3]),
        [level, window, linker] => Some([level, window, linker].map(|token| {
            // The `long` cast to an `int`: its low 32 bits (gcc).
            scan_long(token).map(|value| value as i32)
        })),
        _ => None,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:243-294
/// `s_ParseSegOptions`: the window (`%ld` stored in an `Int4`), locut and hicut (`%le`) of three
/// tokens, the caller's 0 each without one; `None` for another count (NCBI's parse error).
fn parse_seg_options(buffer: &[u8]) -> Option<Seg> {
    match option_tokens(buffer)[..] {
        [] => Some(Seg {
            window: Some(0),
            locut: Some(0.0),
            hicut: Some(0.0),
        }),
        [window, locut, hicut] => Some(Seg {
            // The `long` stored in an `Int4`: the low 32 bits (gcc).
            window: scan_long(window).map(|value| value as i32),
            locut: scan_double(locut),
            hicut: scan_double(hicut),
        }),
        _ => None,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:104-122
/// `s_ParseRepeatOptions`: the database after the first `-d` and its blanks and tabs, `None`
/// without a `-d`.
fn parse_repeat_options(buffer: &[u8]) -> Option<Vec<u8>> {
    let at = buffer.windows(2).position(|pair| pair == b"-d")?;
    let name = &buffer[at + 2..];
    let blanks = name
        .iter()
        .take_while(|&&c| matches!(c, b' ' | b'\t'))
        .count();
    Some(name[blanks..].to_vec())
}

#[cfg(test)]
#[path = "filter_string_tests.rs"]
mod tests;
