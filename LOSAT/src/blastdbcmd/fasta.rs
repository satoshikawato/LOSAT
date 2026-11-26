//! `-outfmt %f` (the default): NCBI's `CBlastDB_FastaFormatter`. The Bioseq of a sequence
//! (`CSeqDB::GetBioseq`: the ids of its first definition line, a title joined from every
//! definition line), the header of `CBlastDeflineUtil::ProcessFastaDeflines` (short ids) or
//! of `CFastaOstream` (`-long_seqids`), and the sequence lines of
//! `CFastaOstream::WriteSequence`.

use std::collections::HashSet;
use std::io::Write;

use crate::report::defline::{bioseq_title_from_source, generate_bioseq_defline};
use crate::report::ostream::CutBuffer;
use crate::seqid::accguide::parse_flags::PARTIAL_OK;
use crate::seqid::label::{bare, bioseq_fasta, fasta};
use crate::seqid::model::Dbtag;
use crate::seqid::parse::parse_ids;
use crate::seqid::rank::{find_best, text_score};
use crate::seqid::{DeflineSet, ParsedSeqId};

use super::bioseq::DbBioseq;
use super::letters::{apply_masks, reverse_complement};
use super::retrieve::{numeric_target, parse_target, Config, Output, Source, Stop};

/// `CBlastDB_FastaFormatter` with its one `CFastaOstream` (`m_fasta`), whose record of the
/// ids written keeps `-long_seqids` from writing one twice.
pub(super) struct FastaFormatter {
    width: u32,
    long_seqids: bool,
    /// `CFastaOstream::m_PreviousWholeIds`: the FASTA strings of the ids written without a
    /// range.
    previous_whole_ids: HashSet<Vec<u8>>,
    /// The record being written, kept from one sequence to the next (bytes and NCBI insertion cuts).
    units: CutBuffer,
}

impl FastaFormatter {
    pub fn new(width: u32, long_seqids: bool) -> Self {
        FastaFormatter {
            width,
            long_seqids,
            previous_whole_ids: HashSet::new(),
            units: CutBuffer::default(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:436-531
    /// Writes the FASTA of sequence `oid` to `out` (nothing for a sequence NCBI's `GetBioseq`
    /// gives no Bioseq for: no definition line, or no letters), after the texts of the taxonomy
    /// lookups of NCBI's two `GetBioseq` calls.
    pub fn write(
        &mut self,
        source: &Source,
        oid: u32,
        config: &Config,
        target: Option<&[u8]>,
        out: &mut Output<'_, '_>,
    ) -> Result<(), Stop> {
        DbBioseq::with(source, oid, false, target, target_id, out, |bioseq, out| {
            self.write_bioseq(source, oid, config, bioseq, out)
        })
    }

    /// The FASTA of the Bioseq `bioseq` of sequence `oid` (`write`).
    fn write_bioseq(
        &mut self,
        source: &Source,
        oid: u32,
        config: &Config,
        bioseq: &DbBioseq,
        out: &mut Output<'_, '_>,
    ) -> Result<(), Stop> {
        let (set, shown, length) = (&bioseq.set, bioseq.shown(), bioseq.length);
        let target_seqid = bioseq.target.as_ref();
        let range = config.range.clipped(length);
        // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:462-489
        // The masks of the algorithm are read first; a Seq-loc of the masks (where the sequence has
        // some) or of the range takes the Bioseq's best id, which a front line without Seq-ids
        // leaves null: NCBI's dereference throws before any byte of the sequence (line RG
        // D1c-1c-b; `-strand minus` alone reads no id).
        let masks = if config.filt_algo != -1 {
            Some(source.masks(oid, config.filt_algo)?)
        } else {
            None
        };
        let no_ids = shown
            .deflines
            .first()
            .is_some_and(|line| line.ids.is_empty());
        if no_ids && (range.is_some() || masks.as_ref().is_some_and(|masks| !masks.is_empty())) {
            return Err(super::retrieve::null_pointer());
        }
        let mut units = std::mem::take(&mut self.units);
        units.clear();
        // NCBI's insertion of the definition lines into m_Out (the long header
        // of `CFastaOstream::Write` is taken as one too, LOSAT's form of it).
        if self.long_seqids {
            // The Bioseq's first source, from every line the filter keeps (the target's first),
            // names the organism of a title made from it: the long header's alone.
            let source_taxid = set.first_source_taxid(target_seqid);
            self.long_header(source, shown, source_taxid, range, config, &mut units.bytes)?;
        } else {
            short_header(shown, range, config, &mut units.bytes)?;
        }
        let mut letters = source.fasta_letters(oid, range)?;
        // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:463-477
        // The approved differences 3, 5 and 6 of `PD-LOSAT-BLASTDBCMD-DEVIATIONS`, where the
        // letters are those of NCBI's own `%s` path: the mask goes on before the reverse
        // complement (NCBI lower-cases the output positions of the reversed letters), over
        // the whole length of a protein (NCBI's whole-sequence mapper keeps the first third
        // of it), and a sequence without masks of the algorithm stays in upper case (NCBI
        // keeps the mask of the batch line before when the Seq-id is the same).
        if let Some(masks) = &masks {
            apply_masks(&mut letters, masks, range);
        }
        if config.strand_minus {
            reverse_complement(&mut letters);
        }
        let width = self.width.max(1) as usize;
        let seq_bytes = letters.len() + letters.len().div_ceil(width);
        // When the whole record fits in NCBI's buffer without triggering a write (no flush occurs
        // inside a FASTA record), every intermediate line insertion and newline put is held whole;
        // skipping per-line cuts yields the identical buffer state.
        let written = if units.bytes.len() + seq_bytes < out.ncbi_room() {
            write_sequence_plain(&letters, width, &mut units.bytes);
            out.write_all(&units.bytes).map_err(Stop::write)
        } else {
            let changes = match &masks {
                Some(masks) => mask_changes(masks, range, letters.len(), config.strand_minus),
                None => Vec::new(),
            };
            units.end_insertion();
            write_sequence(&letters, self.width, &changes, &mut units);
            // `out` writes to a `files::OutputFile` (`retrieve::process`): a failed write does
            // nothing, as in NCBI's `ostream`; a closed pipe ends the retrieval.
            out.insertions(&mut units)
        };
        self.units = units;
        written
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:2927-2975
    // NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:2875-2925
    // NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:3000-3020
    /// The header of `-long_seqids`: the Bioseq's best id for FASTA (`FastaNARank`,
    /// `FastaAARank`) after its GI, the range as `:a-b` (`:cb-a` for the minus strand), and
    /// the title `CDeflineGenerator` makes of the joined title and the Bioseq's first source,
    /// whose taxid is `source_taxid` (line BDC DECISIONS 12: NCBI writes the forward letters
    /// under the `c` header of a range on the minus strand; LOSAT writes the reverse
    /// complement, as the short header's path does).
    fn long_header(
        &mut self,
        source: &Source,
        set: &DeflineSet,
        source_taxid: Option<i32>,
        range: Option<(u32, u32)>,
        config: &Config,
        out: &mut Vec<u8>,
    ) -> Result<(), Stop> {
        let ids = &set.deflines[0].ids;
        if range.is_none() {
            for id in ids {
                let mut key = Vec::new();
                fasta(id, &mut key);
                if !self.previous_whole_ids.insert(key.clone()) {
                    return Err(Stop::Ncbi(format!(
                        "Duplicate Seq-id {} in FASTA output",
                        String::from_utf8_lossy(&key)
                    )));
                }
            }
        }
        out.push(b'>');
        bioseq_fasta(ids, source.protein, out);
        if let Some((from, to)) = range {
            if config.strand_minus {
                out.extend_from_slice(format!(":c{}-{}", to + 1, from + 1).as_bytes());
            } else {
                out.extend_from_slice(format!(":{}-{}", from + 1, to + 1).as_bytes());
            }
        }
        let mut title = set.bioseq_title();
        if config.ctrl_a {
            title = rejoin_defline_title(&title, true)?;
        }
        let title = generated_title(source, set, source_taxid, &title)?;
        if !title.is_empty() {
            if title[0] != b' ' {
                out.push(b' ');
            }
            out.extend_from_slice(&title);
        }
        out.push(b'\n');
        Ok(())
    }
}

/// The id of a `-target_only` target as `GetBioseq` gets it: the GI of a number
/// (`GetBioseq` with `oid` and `gi`, which `x_GetTaxDefline` makes a `CSeq_id` of type `e_Gi` with
/// the number `gi`), else the Seq-id of the text.
fn target_id(target: &[u8]) -> Result<ParsedSeqId, Stop> {
    match numeric_target(target) {
        Some(number) => Ok(ParsedSeqId::Gi(number)),
        None => parse_target(target),
    }
}

/// The ids `CSeq_id::ParseIDs` gives for `text` with `fParse_PartialOK` (none when it throws).
/// `Err`: a parse that posts a warning NCBI prints (an id of several parts of which one
/// fails), whose line LOSAT does not reproduce; blastdbcmd filters out the warning of extra
/// parts (`CSeq_id` 1306.10) alone.
fn parse_ids_partial(text: &[u8]) -> Result<Vec<ParsedSeqId>, Stop> {
    let mut warned = None;
    let parsed = parse_ids(text, PARTIAL_OK, &mut |warning| {
        if warning.subcode != 10 {
            warned = Some(warning.message);
        }
    });
    if let Some(message) = warned {
        return Err(Stop::Refused(format!(
            "a Seq-id in a title that NCBI reads with a warning ({}) is not supported \
             by LOSAT's {}",
            // LOSAT's own refusal text.
            String::from_utf8_lossy(&message),
            super::PROGRAM
        )));
    }
    Ok(parsed
        .map(|ids| ids.into_iter().flatten().collect())
        .unwrap_or_default())
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:630-665
/// The title of the short header: each ` >` part whose first word parses as Seq-ids gets
/// the bare form of their best, after `\x01` with `-ctrl_a`.
fn configure_defline_title(title: &[u8], use_ctrl_a: bool) -> Result<Vec<u8>, Stop> {
    let separator: &[u8] = if use_ctrl_a { b"\x01" } else { b" >" };
    let mut retval = Vec::new();
    for (idx, token) in split_by_pattern(title, b" >").into_iter().enumerate() {
        if idx == 0 {
            retval.extend_from_slice(token);
            continue;
        }
        let pos = token.iter().position(|&c| c == b' ');
        let possible_id = &token[..pos.unwrap_or(token.len())];
        let seqids = parse_ids_partial(possible_id)?;
        match find_best(&seqids, text_score) {
            Some(id) => {
                retval.extend_from_slice(separator);
                bare(id, &mut retval);
                if let Some(pos) = pos {
                    retval.extend_from_slice(&token[pos..]);
                }
            }
            None => {
                retval.extend_from_slice(b" >");
                retval.extend_from_slice(token);
            }
        }
    }
    Ok(retval)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:361-398
/// The title of `-long_seqids -ctrl_a`: `\x01` before each ` >` part whose first word
/// parses as Seq-ids, the part itself unchanged.
fn rejoin_defline_title(title: &[u8], use_ctrl_a: bool) -> Result<Vec<u8>, Stop> {
    let separator: &[u8] = if use_ctrl_a { b"\x01" } else { b" >" };
    let mut retval = Vec::new();
    for (idx, token) in split_by_pattern(title, b" >").into_iter().enumerate() {
        if idx == 0 {
            retval.extend_from_slice(token);
            continue;
        }
        let pos = token.iter().position(|&c| c == b' ');
        let possible_id = &token[..pos.unwrap_or(token.len())];
        if parse_ids_partial(possible_id)?.is_empty() {
            retval.extend_from_slice(b" >");
        } else {
            retval.extend_from_slice(separator);
        }
        retval.extend_from_slice(token);
    }
    Ok(retval)
}

/// `NStr::Split` of `text` at `pattern` with the flag `fSplit_ByPattern`: the pieces between the
/// occurrences of the whole pattern, empty pieces kept.
fn split_by_pattern<'a>(text: &'a [u8], pattern: &[u8]) -> Vec<&'a [u8]> {
    let mut tokens = Vec::new();
    let mut start = 0;
    let mut at = 0;
    while at + pattern.len() <= text.len() {
        if &text[at..at + pattern.len()] == pattern {
            tokens.push(&text[start..at]);
            at += pattern.len();
            start = at;
        } else {
            at += 1;
        }
    }
    tokens.push(&text[start..]);
    tokens
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:999-1032
/// The header of the short ids (`ProcessFastaDeflines`): the title alone after
/// `gnl|BL_ORD_ID|n`, the local id's string, or the bare form of the Bioseq's best id
/// (`Score`) and the configured title; the range as `:a-b ` or `:cb-a ` in the place of the
/// space.
fn short_header(
    set: &DeflineSet,
    range: Option<(u32, u32)>,
    config: &Config,
    out: &mut Vec<u8>,
) -> Result<(), Stop> {
    let ids = &set.deflines[0].ids;
    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:959-985
    // A front line without Seq-ids gives no first id: NCBI writes no header line, the letters
    // follow (line RG D1c-1c-b).
    let Some(first) = ids.first() else {
        return Ok(());
    };
    let range_text = range.map(|(from, to)| {
        if config.strand_minus {
            format!(":c{}-{} ", to + 1, from + 1)
        } else {
            format!(":{}-{} ", from + 1, to + 1)
        }
    });
    let title = set.bioseq_title();
    out.push(b'>');
    match first {
        ParsedSeqId::General(Dbtag { db, .. }) if db == b"BL_ORD_ID" => {
            if let Some(text) = &range_text {
                out.extend_from_slice(text.as_bytes());
            }
            out.extend_from_slice(&title);
        }
        ParsedSeqId::Local(_) => {
            let mut text = Vec::new();
            fasta(first, &mut text);
            out.extend_from_slice(&text[4..]);
            out.extend_from_slice(range_text.as_deref().unwrap_or(" ").as_bytes());
            out.extend_from_slice(&title);
        }
        _ => {
            let best = find_best(ids, text_score).expect("the definition line has ids");
            bare(best, out);
            out.extend_from_slice(range_text.as_deref().unwrap_or(" ").as_bytes());
            out.extend_from_slice(&configure_defline_title(&title, config.ctrl_a)?);
        }
    }
    out.push(b'\n');
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:540
// NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:2983-2985
// NCBI reference (598d8ae6): c++/include/objmgr/util/create_defline.hpp:86-88
// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:621
// The Bioseq's first source (`CSeqDBVol::x_GetTaxonomy`: one for each definition line with a
// taxid above 0, in the order of `x_GetTaxDefline`, which puts the target's line first and keeps
// the others, also with `-target_only`) names its organism where the taxonomy database gives the
// taxid a name: an empty title is made from it (`report::defline::BioseqFacts`), without `unnamed
// protein product` (`m_LocalAnnotsOnly`): `[Homo sapiens]` for a protein, `Homo sapiens` for a
// nucleotide (probe of line RG T5, `T/t5/out/ncbi/c.*.long`; with `-target_only` of a line without
// a taxid, the other line's: review of RG T5, `T/t5/review`, cases `*.to.*`).
// The Bioseq has a source (`BioseqFacts::has_source`) wherever a definition line has a positive
// taxid, named or not: a protein title ending in `]` loses its end to `x_AdjustProteinTitleSuffix`
// then (`report::defline::adjust_protein_title_suffix`; a ` [ ]` also without a name: probe of line
// RG T5c, `T/t5c/out/bdc`, `zeta protein, partial [ ]` -> `zeta protein [Bacillus sp.]`, without
// a name `zeta protein`; `psi protein, partial []` stays, as an empty name matches nothing).
/// The title that `CDeflineGenerator::GenerateDefline` makes (with the flags `fFastaFormat` and
/// `fNoExpensiveOps`) of the Bioseq's title descriptor `title`
/// (`report::defline::generate_bioseq_defline`, with the suffix adjustment of a protein with a
/// source, and an empty title made from the name of the source), the Bioseq's ids those of `set`'s
/// first line and its first source's taxid `source_taxid` (`DeflineSet::first_source_taxid`), whose
/// names are looked up for a title made from the source and for every protein, as `GetBioseq` looks
/// them up.
fn generated_title(
    source: &Source,
    set: &DeflineSet,
    source_taxid: Option<i32>,
    title: &[u8],
) -> Result<Vec<u8>, Stop> {
    let ids = &set.deflines[0].ids;
    let mut bioseq_facts = crate::report::defline::bioseq_facts(ids);
    bioseq_facts.local_annots_only = true;
    bioseq_facts.has_source = source_taxid.is_some();
    let from_source = bioseq_title_from_source(title, || bioseq_facts.clone());
    if let Some(taxid) = source_taxid.filter(|_| from_source || source.protein) {
        bioseq_facts.taxname = source.tax_names(taxid)?.taxname().cloned();
    }
    let defline =
        generate_bioseq_defline(title, source.protein, false, &bioseq_facts).map_err(|_| {
            Stop::Refused(format!(
                "a title whose encoding NCBI cannot tell is not supported by LOSAT's {}",
                super::PROGRAM
            ))
        })?;
    if defline.reads_past_end {
        return Err(Stop::Refused(format!(
            "a title on which NCBI's defline generator reads past its end is not supported \
             by LOSAT's {}",
            super::PROGRAM
        )));
    }
    Ok(defline.text)
}

// NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:3074-3120
/// The places of the letters `len` of the output where NCBI's masking state changes
/// (`x_GetMaskingStates`: the masks `masks` mapped to the range `range`, overlapping and
/// abutting ones merged, each start and each end), as the letters show them: with `minus` the
/// letters were reverse-complemented after the masks went on (`apply_masks`), so are the places.
fn mask_changes(
    masks: &[(u32, u32)],
    range: Option<(u32, u32)>,
    len: usize,
    minus: bool,
) -> Vec<usize> {
    let from = range.map_or(0, |(from, _)| from) as usize;
    let mut spans: Vec<(usize, usize)> = masks
        .iter()
        .filter_map(|&(first, second)| {
            let start = (first as usize).max(from);
            let end = (second as usize).min(from + len);
            (start < end).then(|| (start - from, end - from))
        })
        .collect();
    spans.sort_unstable();
    let mut merged: Vec<(usize, usize)> = Vec::with_capacity(spans.len());
    for (start, end) in spans {
        match merged.last_mut() {
            Some(last) if start <= last.1 => last.1 = last.1.max(end),
            _ => merged.push((start, end)),
        }
    }
    if minus {
        merged.reverse();
        for span in &mut merged {
            *span = (len - span.1, len - span.0);
        }
    }
    merged
        .into_iter()
        .flat_map(|(start, end)| [start, end])
        .collect()
}

/// The letters that NCBI's sequence iterator holds at a time (`CSeqVector_CI`'s `kCacheSize`),
/// filled from the first letter of the output on.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/seq_vector_ci.cpp:56
const SEQ_VECTOR_CACHE: usize = 1024;

/// Plain FASTA line wrapping without per-line cuts, used only when the whole record is proven to
/// fit inside NCBI's stream buffer with zero writes (`ReportWrite::ncbi_room`).
fn write_sequence_plain(letters: &[u8], width: usize, out: &mut Vec<u8>) {
    out.reserve(letters.len() + letters.len().div_ceil(width));
    for line in letters.chunks(width) {
        out.extend_from_slice(line);
        out.push(b'\n');
    }
}

// NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:3164-3320
/// The sequence lines of `CFastaOstream::x_WriteSequence`: lines of `width` letters, the
/// last one shorter, each ended by a line feed, as NCBI's insertions into `out`. NCBI takes the
/// letters a piece at a time: up to the end of its iterator's cache (`SEQ_VECTOR_CACHE`) or to
/// the next change of the masking state (`changes`, increasing), whichever comes first. Each
/// line, or the part of it in a piece, is one insertion (`x_WriteBuffer`), and each line feed a
/// `put` (NCBI's `<<` of a `char`), the last one also after a shorter last line.
fn write_sequence(letters: &[u8], width: u32, changes: &[usize], out: &mut CutBuffer) {
    let width = width.max(1) as usize;
    let len = letters.len();
    out.bytes.reserve(len + len.div_ceil(width));
    out.cuts
        .reserve(2 * len.div_ceil(width) + len / SEQ_VECTOR_CACHE + changes.len() + 1);
    let mut changes = changes.iter().copied().peekable();
    let mut rem_line = width;
    let mut pos = 0;
    while pos < len {
        while changes.next_if(|&change| change <= pos).is_some() {}
        let mut end = ((pos / SEQ_VECTOR_CACHE + 1) * SEQ_VECTOR_CACHE).min(len);
        if let Some(&change) = changes.peek() {
            end = end.min(change);
        }
        let mut at = pos;
        while end - at >= rem_line {
            out.insertion(&letters[at..at + rem_line]);
            out.put(b'\n');
            at += rem_line;
            rem_line = width;
        }
        if end > at {
            out.insertion(&letters[at..end]);
            rem_line -= end - at;
        }
        pos = end;
    }
    if rem_line < width {
        out.put(b'\n');
    }
}

/// `ObjectId` of a `BL_ORD_ID` tag, for tests.
#[cfg(test)]
fn ordinal_id(n: i32) -> ParsedSeqId {
    use crate::seqid::model::ObjectId;
    ParsedSeqId::General(Dbtag {
        db: b"BL_ORD_ID".to_vec(),
        tag: ObjectId::Id(n),
    })
}

#[cfg(test)]
#[path = "fasta_tests.rs"]
mod tests;
