//! The tests of `CFormatGuess` that read the split lines (`m_TestLines`), in the order of `sm_CheckOrder`, with the
//! line tests they call.

use super::{
    c_isprint, c_string, contains, find, split_in_two, tokens, truncate_spaces, FormatGuess,
    Source, Stop, Tokens, F_DNA_MAIN_ALPHABET, SYMBOL_TYPE_TABLE,
};
use crate::ncbi_str::{is_space, string_to_int, string_to_non_negative_int, string_to_uint};

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:196-245
/// `s_IsTokenPosInt`.
fn is_token_pos_int(token: &[u8]) -> bool {
    match token {
        [] => false,
        [b'0'] => true,
        [first, rest @ ..] => (b'1'..=b'9').contains(first) && rest.iter().all(u8::is_ascii_digit),
    }
}

/// `s_IsTokenInteger`.
fn is_token_integer(token: &[u8]) -> bool {
    match token {
        [b'-' | b'+', rest @ ..] => is_token_pos_int(rest),
        _ => is_token_pos_int(token),
    }
}

/// `s_IsTokenDouble`: the first `.` (`ReplaceInPlace` of `.` by `1` in the token, at most once), then a first `-` and a first `0`
/// of a token of two bytes or more, count as `1`.
fn is_token_double(token: &[u8]) -> bool {
    let dot = token.iter().position(|&c| c == b'.');
    // The byte at `i` after the three replacements.
    let at = |i: usize| {
        let mut c = if dot == Some(i) { b'1' } else { token[i] };
        if i == 0 && token.len() > 1 && (c == b'-' || c == b'0') {
            c = b'1';
        }
        c
    };
    match token.len() {
        0 => false,
        1 if at(0) == b'0' => true,
        len => (b'1'..=b'9').contains(&at(0)) && (1..len).all(|i| at(i).is_ascii_digit()),
    }
}

/// Whether the token is one byte long and `find_first_of` finds it in `set` (a result other than `NPOS`).
fn one_of(token: &[u8], set: &[u8]) -> bool {
    token.len() == 1 && set.contains(&token[0])
}

impl<S: Source> FormatGuess<S> {
    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:2085-2122
    /// `TestFormatPsl`.
    pub(super) fn test_format_psl(&mut self) -> bool {
        if !self.ensure_lines() {
            return false;
        }
        let mut lines = self.lines().skip_while(|line| line.starts_with(b"#"));
        let Some(first) = lines.next() else {
            return false;
        };
        let ignore_first_column = if is_line_psl(first, false) {
            false
        } else if is_line_psl(first, true) {
            true
        } else {
            return false;
        };
        lines.all(|line| is_line_psl(line, ignore_first_column))
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:764-773
    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:2851-2943
    // NCBI reference (598d8ae6): c++/include/corelib/ncbistr.hpp:5503-5510
    /// `TestFormatRepeatMasker`. The header test trims the lines up to the first that is not blank in place: the
    /// later tests read them so.
    pub(super) fn test_format_repeat_masker(&mut self) -> bool {
        if !self.ensure_stats_and_lines() {
            return false;
        }
        self.repeat_masker_with_header() || self.repeat_masker_without_header()
    }

    fn repeat_masker_with_header(&mut self) -> bool {
        const LABELS_1ST_LINE: [&[u8]; 5] = [b"SW", b"perc", b"query", b"position", b"matching"];
        const LABELS_2ND_LINE: [&[u8]; 5] = [b"score", b"div.", b"del.", b"ins.", b"sequence"];
        // The lines up to the first that is not blank are trimmed (`extra`, added by a later test, is never here).
        debug_assert!(self.lines.extra.is_none());
        let count = self.lines.count;
        let at = self
            .lines()
            .position(|line| !truncate_spaces(line).is_empty())
            .unwrap_or(count);
        self.lines.trimmed = (at + 1).min(count);
        let labelled = |line: &[u8], labels: &[&[u8]]| {
            let mut current_offset = 0;
            labels.iter().all(|label| {
                match find(line.get(current_offset..).unwrap_or(&[]), label) {
                    Some(pos) => {
                        current_offset += pos;
                        true
                    }
                    None => false,
                }
            })
        };
        let mut lines = self.lines().skip(at);
        lines
            .next()
            .is_some_and(|line| labelled(line, &LABELS_1ST_LINE))
            && lines
                .next()
                .is_some_and(|line| labelled(line, &LABELS_2ND_LINE))
            && lines.next().is_some()
    }

    fn repeat_masker_without_header(&self) -> bool {
        self.lines().all(|line| {
            let str = truncate_spaces(line);
            str.is_empty() || is_line_rmo(str)
        })
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:2500-2518
    /// `TestFormatVcf`.
    pub(super) fn test_format_vcf(&mut self) -> bool {
        self.ensure_stats_and_lines()
            && self.lines_may_hold(b"##fileformat=VCFv")
            && self
                .lines()
                .any(|line| line.starts_with(b"##fileformat=VCFv"))
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:863-898
    /// `TestFormatGvf`.
    pub(super) fn test_format_gvf(&mut self) -> Result<bool, Stop> {
        if !self.ensure_lines() {
            return Ok(false);
        }
        let mut count = 0;
        for line in self.lines() {
            if line.is_empty() || line[0] == b'#' {
                if line.starts_with(b"##gvf-version") {
                    return Ok(true);
                }
                continue;
            }
            if count == 0 && (line.starts_with(b"browser ") || line.starts_with(b"track ")) {
                continue;
            }
            if !is_line_gvf(line)? {
                return Ok(false);
            }
            count += 1;
        }
        Ok(count != 0)
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:901-936
    /// `TestFormatGff3` (after `EnsureTestBuffer()` and then `EnsureSplitLines()`, both true).
    pub(super) fn test_format_gff3(&mut self) -> bool {
        if !self.ensure_lines() {
            return false;
        }
        let mut count = 0;
        for line in self.lines() {
            if count == 0 && line.starts_with(b"##gff-version") {
                return line.starts_with(b"##gff-version 3");
            }
            if line.is_empty() || line[0] == b'#' {
                continue;
            }
            if count == 0 && (line.starts_with(b"browser ") || line.starts_with(b"track ")) {
                continue;
            }
            if !is_line_gff3(line) {
                return false;
            }
            count += 1;
        }
        count != 0
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:829-860 (`TestFormatGtf`), :977-1009
    // (`TestFormatGff2`, the same loop with `IsLineGff2`)
    /// The loop of `TestFormatGtf` and `TestFormatGff2` with its line test.
    fn gff_lines(&mut self, is_line: fn(&[u8]) -> bool) -> bool {
        if !self.ensure_lines() {
            return false;
        }
        let mut count = 0;
        for line in self.lines() {
            if line.is_empty() || line[0] == b'#' {
                continue;
            }
            if count == 0 && (line.starts_with(b"browser ") || line.starts_with(b"track ")) {
                continue;
            }
            if !is_line(line) {
                return false;
            }
            count += 1;
        }
        count != 0
    }

    /// `TestFormatGtf`.
    pub(super) fn test_format_gtf(&mut self) -> bool {
        self.gff_lines(is_line_gtf)
    }

    /// `TestFormatGff2`.
    pub(super) fn test_format_gff2(&mut self) -> bool {
        self.gff_lines(is_line_gff2)
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:939-974
    /// `TestFormatAugustus` (after `EnsureTestBuffer()` and then `EnsureSplitLines()`, both true).
    pub(super) fn test_format_augustus(&mut self) -> bool {
        if !self.ensure_lines() {
            return false;
        }
        let mut count = 0;
        for line in self.lines() {
            if count == 0 && line.starts_with(b"##gff-version 3") {
                return false;
            }
            if line.is_empty() || line[0] == b'#' {
                continue;
            }
            if count == 0 && (line.starts_with(b"browser ") || line.starts_with(b"track ")) {
                return false;
            }
            if !is_line_augustus(line) {
                return false;
            }
            count += 1;
        }
        count != 0
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1012-1038
    /// `TestFormatGlimmer3`: `Err` when a frame is beyond `int` (`NStr::StringToInt` throws).
    pub(super) fn test_format_glimmer3(&mut self) -> Result<bool, Stop> {
        if !self.ensure_lines() {
            return Ok(false);
        }
        if self.line_at(0).and_then(<[u8]>::first) != Some(&b'>') || self.lines.len() < 2 {
            return Ok(false);
        }
        for line in self.lines().skip(1) {
            if !is_line_glimmer3(line)? {
                return Ok(false);
            }
        }
        Ok(true)
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1041-1058
    /// `TestFormatAgp` (`IsLineAgp` throws nothing).
    pub(super) fn test_format_agp(&mut self) -> bool {
        self.ensure_lines() && self.lines().all(is_line_agp)
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1861-1890
    /// `TestFormatWiggle`. `NStr::Find` as a condition is true unless it is 0 (`NPOS` is true too).
    pub(super) fn test_format_wiggle(&mut self) -> bool {
        if !self.ensure_stats_and_lines() {
            return false;
        }
        if ![&b"track"[..], b"fixedStep", b"variableStep"]
            .iter()
            .any(|start| self.lines_may_hold(start))
        {
            return false;
        }
        let found = |line: &[u8], pattern: &[u8]| find(line, pattern) != Some(0);
        self.lines().any(|line| {
            (line.starts_with(b"track")
                && (contains(line, b"type=wiggle_0") || contains(line, b"type=bedGraph")))
                || (line.starts_with(b"fixedStep")
                    && found(line, b"chrom=")
                    && found(line, b"start="))
                || line.starts_with(b"variableStep")
        })
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1747-1807
    /// `TestFormatBed`.
    pub(super) fn test_format_bed(&mut self) -> bool {
        if !self.ensure_stats_and_lines() {
            return false;
        }
        let mut track_line_found = false;
        let mut has_start_and_stop = false;
        let mut columncount = 0;
        for line in self.lines() {
            let mut str = std::borrow::Cow::Borrowed(truncate_spaces(line));
            if str.is_empty() {
                continue;
            }
            if str.starts_with(b"chr ") || str.starts_with(b"Chr ") || str.starts_with(b"CHR ") {
                str.to_mut().remove(3);
            }
            if str.starts_with(b"track") {
                track_line_found = true;
                continue;
            }
            if str.starts_with(b"browser") || str.starts_with(b"#") {
                continue;
            }
            let columns = tokens(&str, b" \t");
            if columns.len() < 3 || columns.len() > 12 {
                return false;
            }
            if columns.len() != columncount {
                if columncount == 0 {
                    columncount = columns.len();
                } else {
                    return false;
                }
            }
            if is_token_pos_int(columns.at(1)) && is_token_pos_int(columns.at(2)) {
                has_start_and_stop = true;
            }
        }
        has_start_and_stop || track_line_found
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1810-1858
    /// `TestFormatBed15`.
    pub(super) fn test_format_bed15(&mut self) -> bool {
        if !self.ensure_stats_and_lines() {
            return false;
        }
        let mut line_found = false;
        for line in self.lines() {
            if truncate_spaces(line).is_empty()
                || line.starts_with(b"track")
                || line.starts_with(b"browser")
                || line.starts_with(b"#")
            {
                continue;
            }
            let columns = tokens(line, b" \t");
            if columns.len() != 15
                || ![1, 2, 4, 6, 7]
                    .iter()
                    .all(|&i| is_token_pos_int(columns.at(i)))
                || !matches!(truncate_spaces(columns.at(5)), b"+" | b"-")
            {
                return false;
            }
            line_found = true;
        }
        line_found
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1893-1925
    /// `TestFormatHgvs` (its own buffer, as `TestFormatNewick`'s, is never read in `GuessFormat`, whose
    /// `TestFormatNewick` before it leaves lines).
    pub(super) fn test_format_hgvs(&mut self) -> bool {
        if !self.ensure_stats() || !self.ensure_split_lines() {
            const BUFFSIZE: usize = 1024;
            let read = self.stream.read(BUFFSIZE);
            self.stream.clear(); // in case we reached eof
            self.stream.stepback(&read);
            self.buffer = Some(read);
            self.terminated = true;
            let line = self.c_string_line();
            self.push_line(line);
        }
        let mut count = 0;
        for line in self.lines() {
            if line.is_empty() || line[0] == b'#' {
                continue;
            }
            if !is_line_hgvs(line) {
                return false;
            }
            count += 1;
        }
        count != 0
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1189-1234
    /// `TestFormatDistanceMatrix`. A last line without a token: libstdc++'s `++` of the end of an empty
    /// `std::list` is the end again (the list is a ring through its sentinel), so its value loop does not run.
    pub(super) fn test_format_distance_matrix(&mut self) -> bool {
        if !self.ensure_lines() {
            return false;
        }
        let first = tokens(self.line_at(0).unwrap_or_default(), b"\t ");
        if first.len() != 1 || !first.at(0).iter().all(u8::is_ascii_digit) {
            return false;
        }
        let rest = self.lines.len() - 1;
        for (i, line) in self.lines().skip(1).enumerate() {
            let mut toks = line
                .split(|&c| c == b'\t' || c == b' ')
                .filter(|token| !token.is_empty());
            let count = toks.clone().count();
            if count != i + 1 && i + 1 != rest {
                return false;
            }
            if !toks.by_ref().skip(1).all(is_token_double) {
                return false;
            }
        }
        true
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1237-1251
    /// `TestFormatFlatFileSequence`.
    pub(super) fn test_format_flat_file_sequence(&mut self) -> bool {
        self.ensure_lines() && self.lines().all(is_line_flat_file_sequence)
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:2400-2497
    /// `TestFormatFlatFileUniProt` (`SplitInTwo` throws nothing).
    pub(super) fn test_format_flat_file_uniprot(&mut self) -> bool {
        if !self.ensure_stats_and_lines() || self.lines.len() < 15 {
            return false;
        }
        let mut records = self.lines().map(|line| split_in_two(line, b" ", true).0);
        if records.next() != Some(b"ID".as_slice()) {
            return false;
        }
        let Some(mut line_code) = records.next() else {
            return false;
        };
        if line_code != b"AC" {
            return false;
        }
        while line_code == b"AC" {
            let Some(next) = records.next() else {
                return false;
            };
            line_code = next;
        }
        for _ in 0..3 {
            if line_code != b"DT" {
                return false;
            }
            let Some(next) = records.next() else {
                return false;
            };
            line_code = next;
        }
        if line_code != b"DE" {
            return false;
        }
        while line_code == b"DE" {
            let Some(next) = records.next() else {
                return true;
            };
            line_code = next;
        }
        line_code == b"GN" || line_code == b"OS"
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:2270-2397
    /// `TestFormatFlatFileEna`.
    pub(super) fn test_format_flat_file_ena(&mut self) -> bool {
        if !self.ensure_stats_and_lines() || self.lines.len() < 19 {
            return false;
        }
        if split_in_two(self.line_at(0).unwrap_or_default(), b" ", true).0 != b"ID" {
            return false;
        }
        let mut records = self
            .lines()
            .skip(1)
            .filter(|line| !line.starts_with(b"XX"))
            .map(|line| split_in_two(line, b" ", true).0);
        let Some(mut line_code) = records.next() else {
            return false;
        };
        if line_code != b"AC" {
            return false;
        }
        while line_code == b"AC" {
            let Some(next) = records.next() else {
                return false;
            };
            line_code = next;
        }
        while line_code == b"PR" {
            let Some(next) = records.next() else {
                return false;
            };
            line_code = next;
        }
        for _ in 0..2 {
            if line_code != b"DT" {
                return false;
            }
            let Some(next) = records.next() else {
                return false;
            };
            line_code = next;
        }
        for looking_for in [&b"DE"[..], b"KW", b"OS", b"OC"] {
            if line_code != looking_for {
                return false;
            }
            while line_code == looking_for {
                let Some(next) = records.next() else {
                    return true;
                };
                line_code = next;
            }
        }
        true
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:2125-2165
    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:2168-2267
    /// `TestFormatFlatFileGenbank`.
    pub(super) fn test_format_flat_file_genbank(&mut self) -> bool {
        if !self.ensure_stats_and_lines() || self.lines.len() < 9 {
            return false;
        }
        if split_in_two(self.line_at(0).unwrap_or_default(), b" ", true).0 != b"LOCUS" {
            return false;
        }
        let mut records = self.lines().skip(1);
        // `GenbankGetKeywordLine`: the keyword of the next line (after the continuation lines), the lines past it.
        let mut next_keyword = || -> Option<&[u8]> {
            let mut current = records.next()?;
            if current.len() > 79 {
                return None;
            }
            const VALID_INDENTS: [usize; 6] = [0, 2, 3, 5, 12, 21];
            let first_not = |line: &[u8], set: &[u8]| line.iter().position(|c| !set.contains(c));
            let mut first_not_blank = first_not(current, b" ");
            while first_not_blank != Some(0) {
                if !first_not_blank.is_some_and(|indent| VALID_INDENTS.contains(&indent))
                    && first_not(current, b" 1234567890") != Some(10)
                {
                    return None;
                }
                current = records.next()?;
                first_not_blank = first_not(current, b" ");
            }
            Some(split_in_two(current, b" ", true).0)
        };
        let Some(mut keyword) = next_keyword() else {
            return false;
        };
        for looking_for in [&b"DEFINITION"[..], b"ACCESSION"] {
            if keyword != looking_for {
                return false;
            }
            while keyword == looking_for {
                let Some(next) = next_keyword() else {
                    return false;
                };
                keyword = next;
            }
        }
        let nid_seen = keyword == b"NID";
        if nid_seen {
            let Some(next) = next_keyword() else {
                return false;
            };
            keyword = next;
        }
        if keyword != b"VERSION" {
            return false;
        }
        let Some(next) = next_keyword() else {
            return false;
        };
        keyword = next;
        if !nid_seen && keyword == b"NID" {
            let Some(next) = next_keyword() else {
                return false;
            };
            keyword = next;
        }
        for looking_for in [&b"PROJECT"[..], b"DBLINK"] {
            while keyword == looking_for {
                let Some(next) = next_keyword() else {
                    return false;
                };
                keyword = next;
            }
        }
        keyword == b"KEYWORDS"
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1254-1274
    /// `TestFormatFiveColFeatureTable`: the first line that is not empty starts with `>Feature ` or `>Features `
    /// (true when there is none: the empty line `TestFormatNewick` adds for a buffer that starts with NUL).
    pub(super) fn test_format_five_col_feature_table(&mut self) -> bool {
        self.ensure_lines()
            && self
                .lines()
                .find(|line| !line.is_empty())
                .is_none_or(|line| {
                    line.starts_with(b">Feature ") || line.starts_with(b">Features ")
                })
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1727-1744
    /// `TestFormatSnpMarkers`.
    pub(super) fn test_format_snp_markers(&mut self) -> bool {
        self.ensure_lines()
            && self.lines_may_hold(b"rs")
            && self
                .lines()
                .any(|line| line.starts_with(b"rs") && scans_snp_marker(c_string(line)))
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1316-1345
    /// `TestFormatAlignment`.
    pub(super) fn test_format_alignment(&mut self) -> bool {
        self.ensure_lines()
            && (self.test_format_clustal()
                || (self.lines_may_hold(b"#NEXUS")
                    && self.lines().any(|line| contains(line, b"#NEXUS"))))
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1434-1576
    /// `TestFormatCLUSTAL`. `NcbiGetline` with `\r\n` as the delimiters (cited at `stats_of`) takes `\r\n` as one line end; a last
    /// line without an end sets `eofbit` (with `failbit` too when it is empty), and the read after a last line end
    /// finds nothing (`eofbit` and `failbit`).
    fn test_format_clustal(&mut self) -> bool {
        if !self.ensure_test_buffer() {
            return false;
        }
        let data = self.data();
        let buffer_full = data.len() == self.buffer_size;
        let mut pos = 0;
        let mut eof = false;
        let mut in_block = false;
        let mut block_size = 0usize;
        let mut ids: Vec<&[u8]> = Vec::new();
        let mut has_valid_block = false;
        let mut seg_length_prev = 0usize;
        while !eof {
            let start = pos;
            while pos < data.len() && data[pos] != b'\r' && data[pos] != b'\n' {
                pos += 1;
            }
            let line = &data[start..pos];
            let fail;
            if pos == data.len() {
                eof = true;
                fail = line.is_empty();
            } else {
                fail = false;
                if data[pos] == b'\r' && data.get(pos + 1) == Some(&b'\n') {
                    pos += 1;
                }
                pos += 1;
            }
            if (buffer_full && eof) || fail {
                break;
            }
            if line.starts_with(b"CLUSTAL") {
                continue;
            }
            if line.iter().all(|&c| is_space(c)) {
                if in_block {
                    if block_size < 2 {
                        return false;
                    }
                    in_block = false;
                    block_size = 0;
                    ids.clear();
                }
                continue;
            }
            if looks_like_clustal_conserved_info(line) {
                if !in_block || block_size < 2 {
                    return false;
                }
                in_block = false;
                block_size = 0;
                ids.clear();
                continue;
            }
            let Some((seq_id, seg_length)) = try_process_clustal_seq_data(line) else {
                return false;
            };
            if seg_length > 60 {
                return false;
            }
            if in_block {
                if seg_length != seg_length_prev {
                    return false;
                }
                has_valid_block = true;
            }
            if ids.contains(&seq_id) {
                return false;
            }
            ids.push(seq_id);
            seg_length_prev = seg_length;
            in_block = true;
            block_size += 1;
        }
        has_valid_block
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1579-1621
    /// `TestFormatTable`.
    pub(super) fn test_format_table(&mut self) -> bool {
        if !self.ensure_lines() || !self.is_ascii_text() {
            return false;
        }
        [&b" "[..], b" \t", b"\t", b",", b"|"]
            .iter()
            .any(|delims| self.test_table_delimiter(delims))
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1521-1576 (`x_TestTableDelimiter`)
    /// `x_TestTableDelimiter` of `delims`: the count starts again at the line that gave the columns, so that line counts
    /// twice.
    fn test_table_delimiter(&self, delims: &[u8]) -> bool {
        let skip = (5..7).filter(|&i| self.lines.len() > i).count();
        let commented = |line: &[u8]| line.is_empty() || line[0] == b'#' || line[0] == b';';
        let Some((start, first)) = self
            .lines()
            .enumerate()
            .skip(skip)
            .find(|(_, line)| !commented(line))
        else {
            return false;
        };
        let ncols = tokens(first, delims).len();
        if ncols < 2 {
            return false;
        }
        let data_size = self.data().len();
        let mut nlines = 1;
        for (i, line) in self.lines().enumerate().skip(start) {
            if commented(line) {
                continue;
            }
            if tokens(line, delims).len() != ncols {
                if i + 1 != self.lines.len() || data_size < self.buffer_size {
                    return false;
                }
            } else {
                nlines += 1;
            }
            // Every token printable: every byte that is not a delimiter.
            if line.iter().any(|c| !delims.contains(c) && !c_isprint(*c)) {
                return false;
            }
        }
        nlines >= 3
    }

    // NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:801-826
    /// `TestFormatPhrapAce`.
    pub(super) fn test_format_phrap_ace(&mut self) -> Result<bool, Stop> {
        if !self.ensure_lines() || self.data().contains(&0) {
            return Ok(false);
        }
        let mut found_id = false;
        for line in self.lines() {
            if found_id {
                if looks_like_nuc_seq_data(line) {
                    return Ok(true);
                }
            } else if is_line_phrap_id(line)? {
                found_id = true;
            }
        }
        Ok(false)
    }
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3615-3677
/// `IsLinePsl` of `line` with `ignoreFirstLine`.
fn is_line_psl(line: &[u8], ignore_first_column: bool) -> bool {
    let tokens = tokens(line, b" \t");
    let first = usize::from(ignore_first_column);
    if tokens.len().wrapping_sub(first) != 21 {
        return false;
    }
    let token = tokens.at(first + 8);
    (first..first + 8).all(|column| is_token_pos_int(tokens.at(column)))
        && !token.is_empty()
        && token.len() <= 2
        && token.iter().all(|&c| c == b'-' || c == b'+')
        && (first + 10..first + 13).all(|column| is_token_pos_int(tokens.at(column)))
        && (first + 14..first + 18).all(|column| is_token_pos_int(tokens.at(column)))
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3542-3612
/// `IsLineRmo` of `line`.
fn is_line_rmo(line: &[u8]) -> bool {
    let values = tokens(line, b" \t");
    values.len() >= 14
        && is_token_pos_int(values.at(0))
        && (1..4).all(|i| is_token_double(values.at(i)))
        && is_token_pos_int(values.at(5))
        && is_token_pos_int(values.at(6))
        && matches!(values.at(8), b"+" | b"C")
}

/// The columns 3 to 7 of the GFF family (`tokens[3]` to `tokens[7]` of `IsLineGtf`, `IsLineGff3`, `IsLineGff2`):
/// two positions, a score, a strand of `strands`, a phase.
fn gff_columns(tokens: &Tokens<'_>, strands: &[u8]) -> bool {
    is_token_pos_int(tokens.at(3))
        && is_token_pos_int(tokens.at(4))
        && is_token_double(tokens.at(5))
        && one_of(tokens.at(6), strands)
        && one_of(tokens.at(7), b".0123")
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3242-3270
/// `IsLineGtf` of `line`.
fn is_line_gtf(line: &[u8]) -> bool {
    let tokens = tokens(line, b" \t");
    tokens.len() >= 9
        && gff_columns(&tokens, b".+-")
        && (contains(tokens.at(8), b"gene_id") || contains(tokens.at(8), b"transcript_id"))
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3273-3346
/// `IsLineGvf` of `line`: a line of eight columns that passes the column tests reads `tokens[8]` past the end of the
/// vector (`Stop::Undecided`).
fn is_line_gvf(line: &[u8]) -> Result<bool, Stop> {
    const TERMS: [&[u8]; 25] = [
        b"snv",
        b"cnv",
        b"copy_number_variation",
        b"gain",
        b"copy_number_gain",
        b"loss",
        b"copy_number_loss",
        b"loss_of_heterozygosity",
        b"complex",
        b"complex_substitution",
        b"complex_sequence_alteration",
        b"indel",
        b"insertion",
        b"inversion",
        b"substitution",
        b"deletion",
        b"duplication",
        b"translocation",
        b"upd",
        b"uniparental_disomy",
        b"maternal_uniparental_disomy",
        b"paternal_uniparental_disomy",
        b"tandom_duplication",
        b"structural_variation",
        b"sequence_alteration",
    ];
    let tokens = tokens(line, b" \t");
    if tokens.len() < 8
        || !is_token_pos_int(tokens.at(3))
        || !is_token_pos_int(tokens.at(4))
        || !TERMS
            .iter()
            .any(|term| term.eq_ignore_ascii_case(tokens.at(2)))
        || !is_token_double(tokens.at(5))
        || !one_of(tokens.at(6), b".+-")
        || !one_of(tokens.at(7), b".0123")
    {
        return Ok(false);
    }
    let Some(attrs) = tokens.get(8) else {
        return Err(Stop::Undecided(
            "TestFormatGvf reads the ninth column of a GVF line of eight",
        ));
    };
    Ok(contains(attrs, b"ID=") && contains(attrs, b"Variant_seq="))
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3349-3390
/// `IsLineGff3` of `line`.
fn is_line_gff3(line: &[u8]) -> bool {
    let tokens = tokens(line, b" \t");
    if tokens.len() < 9 || !gff_columns(&tokens, b".+-?") {
        return false;
    }
    let col9 = tokens.at(8);
    col9.len() <= 1
        || [
            &b"ID"[..],
            b"Parent",
            b"Target",
            b"Name",
            b"Alias",
            b"Note",
            b"Dbxref",
            b"Xref",
        ]
        .iter()
        .any(|key| contains(col9, key))
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3486-3511
/// `IsLineGff2` of `line`.
fn is_line_gff2(line: &[u8]) -> bool {
    let tokens = tokens(line, b" \t");
    tokens.len() >= 8 && gff_columns(&tokens, b".+-")
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3393-3483
/// `IsLineAugustus` of `line`: `SplitInTwo` without flags, so a run of delimiters gives empty columns.
fn is_line_augustus(line: &[u8]) -> bool {
    // The columns 0 to 7: an id, a method, the feature type, start, stop, score, strand and phase.
    let tests: [fn(&[u8]) -> bool; 8] = [
        |_| true,
        |_| true,
        |_| true,
        is_token_pos_int,
        is_token_pos_int,
        is_token_double,
        |head| one_of(head, b"+-.?"),
        |head| one_of(head, b".0123"),
    ];
    let mut remaining = line;
    let mut feature_type: &[u8] = &[];
    for (column, test) in tests.iter().enumerate() {
        let (head, tail, found) = split_in_two(remaining, b" \t", false);
        if !found || !test(head) {
            return false;
        }
        if column == 2 {
            feature_type = head;
        }
        remaining = tail;
    }
    if remaining.is_empty() {
        return false;
    }
    if feature_type == b"gene" || feature_type == b"transcript" {
        return !contains(remaining, b";") && !contains(remaining, b" ");
    }
    contains(remaining, b"transcript_id") && contains(remaining, b"gene_id")
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3201-3239
/// `IsLineGlimmer3` of `line`: `NStr::StringToInt` with its default flags throws for a frame beyond `int`
/// (`CStringException`, NCBI BLAST+ 2.17.0 psiblast `-in_pssm`: "Cannot convert string '99999999999' to int,
/// overflow", exit status 255).
fn is_line_glimmer3(line: &[u8]) -> Result<bool, Stop> {
    let toks = tokens(line, b"\t ");
    if toks.len() != 5
        || !is_token_integer(toks.at(1))
        || !is_token_integer(toks.at(2))
        || !is_token_integer(toks.at(3))
    {
        return Ok(false);
    }
    let frame = string_to_int(toks.at(3), 0).map_err(Stop::Exception)?;
    Ok((-3..=3).contains(&frame) && is_token_double(toks.at(4)))
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3129-3198
/// `IsLineAgp` of `strLine`.
fn is_line_agp(line: &[u8]) -> bool {
    let line = match find(line, b"#") {
        Some(comment) => &line[..comment],
        None => line,
    };
    let line = truncate_spaces(line);
    if line.is_empty() {
        return true;
    }
    let tokens = tokens(line, b" \t");
    if tokens.len() < 8 {
        return false;
    }
    // A first `-` counts as `1`: a number of 11 bytes or fewer is converted from a copy (`1` and ten more digits
    // are 10^10 or more, past `INT_MAX`, as any longer number).
    let non_negative_after_dash = |token: &[u8]| {
        if token.len() > 1 && token[0] == b'-' {
            if token.len() > 11 {
                return false;
            }
            let mut copy = [0u8; 11];
            copy[..token.len()].copy_from_slice(token);
            copy[0] = b'1';
            return string_to_non_negative_int(&copy[..token.len()]) != -1;
        }
        string_to_non_negative_int(token) != -1
    };
    if !(1..4).all(|i| non_negative_after_dash(tokens.at(i))) || !one_of(tokens.at(4), b"ADFGPNOW")
    {
        return false;
    }
    if tokens.at(4) == b"N" {
        string_to_non_negative_int(tokens.at(5)) != -1
    } else {
        string_to_non_negative_int(tokens.at(6)) != -1
            && string_to_non_negative_int(tokens.at(7)) != -1
            && tokens.len() == 9
            && one_of(tokens.at(8), b"+-")
    }
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3795-3843
/// `IsLineHgvs` of `line`.
fn is_line_hgvs(line: &[u8]) -> bool {
    let mut state = 0;
    let mut i = 0;
    while i < line.len() {
        let c = line[i];
        let next = line.get(i + 1).copied().unwrap_or(0);
        match state {
            0 if c.is_ascii_alphanumeric() => state = 1,
            1 if c == b':' => state = 2,
            2 => {
                if !matches!(c, b'g' | b'c' | b'r' | b'p' | b'n' | b'm') {
                    return false;
                }
                state = 3;
                if c == b'm' && next == b't' {
                    i += 1;
                }
            }
            3 => {
                if c != b'.' {
                    return false;
                }
                state = 4;
            }
            4 if c.is_ascii_alphanumeric() => state = 5,
            _ => {}
        }
        i += 1;
    }
    state == 5
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3074-3100
/// `IsLineFlatFileSequence` of `line`.
fn is_line_flat_file_sequence(line: &[u8]) -> bool {
    let Some(pos) = line.iter().position(|c| !b"0123456789 \t".contains(c)) else {
        return false;
    };
    if pos + 45 >= line.len() {
        return false;
    }
    line[pos..pos + 45].iter().enumerate().all(|(i, &c)| {
        if i % 11 == 10 {
            is_space(c)
        } else {
            c.is_ascii_alphabetic() || c == b'-' || c == b'*'
        }
    })
}

/// The glibc `sscanf` of `line` with the format `"rs%d\t%d\t%d"` that converts all three values, on
/// a C string: `rs`, then three conversions of `%d` (white space, a sign, one digit or more; the
/// value does not matter), each `\t` of the format any white space, none included.
fn scans_snp_marker(line: &[u8]) -> bool {
    let Some(mut rest) = line.strip_prefix(b"rs") else {
        return false;
    };
    for i in 0..3 {
        if i > 0 {
            rest = &rest[rest.iter().take_while(|&&c| is_space(c)).count()..];
        }
        rest = &rest[rest.iter().take_while(|&&c| is_space(c)).count()..];
        if let [b'+' | b'-', after @ ..] = rest {
            rest = after;
        }
        let digits = rest.iter().take_while(|c| c.is_ascii_digit()).count();
        if digits == 0 {
            return false;
        }
        rest = &rest[digits..];
    }
    true
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1348-1362
/// `x_LooksLikeCLUSTALConservedInfo` of `line`.
fn looks_like_clustal_conserved_info(line: &[u8]) -> bool {
    line.iter()
        .all(|&c| is_space(c) || matches!(c, b':' | b'*' | b'.'))
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:1366-1407
// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:569-590 (`S2N_CONVERT_ERROR`: with
// `fConvErr_NoThrow` the conversion returns 0)
/// `x_TryProcessCLUSTALSeqData` of `line`, which gives `id` and `seg_length`: the id and the length.
fn try_process_clustal_seq_data(line: &[u8]) -> Option<(&[u8], usize)> {
    let toks = tokens(line, b" \t");
    if toks.len() != 2 && toks.len() != 3 {
        return None;
    }
    let seqdata = toks.at(1);
    let mut cumulated_res = 0;
    if toks.len() == 3 {
        cumulated_res = string_to_uint(toks.at(2)).unwrap_or(0);
        if cumulated_res == 0 {
            return None;
        }
    }
    if !strict_sequence_type(seqdata) {
        return None;
    }
    if toks.len() == 3 {
        let num_gaps = seqdata.iter().filter(|&&c| c == b'-').count();
        if (seqdata.len() - num_gaps) as u64 > u64::from(cumulated_res) {
            return None;
        }
    }
    Some((toks.at(0), seqdata.len()))
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:307-374 (`CFormatGuess::SequenceType`, from `case
// eST_Strict`)
/// Whether `SequenceType` of `str` with its size and `eST_Strict` is other than `eUndefined`.
fn strict_sequence_type(str: &[u8]) -> bool {
    use super::{F_ALPHA, F_DIGIT, F_DNA_AMBIG_ALPHABET, F_PROTEIN_ALPHABET, F_SPACE};
    let (mut main_nuc, mut ambig, mut bad_nuc, mut amino_acid, mut exotic_aa, mut bad_aa) =
        (0u32, 0u32, 0u32, 0u32, 0u32, 0u32);
    for &c in str {
        let type_ = SYMBOL_TYPE_TABLE[c as usize];
        if type_ & F_DNA_MAIN_ALPHABET != 0 {
            main_nuc += 1;
        } else if type_ & F_DNA_AMBIG_ALPHABET != 0 {
            ambig += 1;
        } else if type_ & (F_SPACE | F_DIGIT) == 0 {
            bad_nuc += 1;
        }
        if type_ & F_PROTEIN_ALPHABET != 0 {
            amino_acid += 1;
        } else if type_ & F_ALPHA != 0 {
            exotic_aa += 1;
        } else if type_ & (F_SPACE | F_DIGIT) == 0 {
            bad_aa += 1;
        }
    }
    (bad_nuc == 0 && ambig <= main_nuc / 3) || (bad_aa == 0 && exotic_aa <= amino_acid / 9)
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:3514-3539
/// `IsLinePhrapId` of `line`: an `AS` line of fewer than three tokens reads past the end of the vector
/// (`Stop::Undecided`) where `&&` does not stop first.
fn is_line_phrap_id(line: &[u8]) -> Result<bool, Stop> {
    let values = tokens(line, b" \t");
    match values.get(0) {
        Some(b"DNA") => Ok(true),
        Some(b"AS") => {
            const PAST: &str = "TestFormatPhrapAce reads past the tokens of an AS line";
            let Some(second) = values.get(1) else {
                return Err(Stop::Undecided(PAST));
            };
            if string_to_non_negative_int(second) < 0 {
                return Ok(false);
            }
            let Some(third) = values.get(2) else {
                return Err(Stop::Undecided(PAST));
            };
            Ok(string_to_non_negative_int(third) >= 0)
        }
        _ => Ok(false),
    }
}

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:777-798
/// `s_LooksLikeNucSeqData` of `line`: `nucCount` divided by the line's size is an integer division,
/// so every byte is a letter of ACGNTU.
fn looks_like_nuc_seq_data(line: &[u8]) -> bool {
    if line.len() < 10 {
        return false;
    }
    let mut nuc_count = 0usize;
    for &c in line {
        if c.is_ascii_alphabetic() {
            if SYMBOL_TYPE_TABLE[c as usize] & F_DNA_MAIN_ALPHABET != 0 {
                nuc_count += 1;
            }
            continue;
        }
        if !is_space(c) {
            return false;
        }
    }
    (nuc_count / line.len()) as f64 > 0.9
}
