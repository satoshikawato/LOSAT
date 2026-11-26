//! `-outfmt %f` (the default): NCBI's `CBlastDB_FastaFormatter`. The Bioseq of a sequence
//! (`CSeqDB::GetBioseq`: the ids of its first definition line, a title joined from every
//! definition line), the header of `CBlastDeflineUtil::ProcessFastaDeflines` (short ids) or
//! of `CFastaOstream` (`-long_seqids`), and the sequence lines of
//! `CFastaOstream::WriteSequence`.

use std::collections::HashSet;

use crate::report::defline::{bioseq_title_from_source, generate_bioseq_defline};
use crate::report::ostream::StreamInserter;
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
    /// The formatted header line of the current record, reused across sequences.
    header: Vec<u8>,
}

impl FastaFormatter {
    pub fn new(width: u32, long_seqids: bool) -> Self {
        FastaFormatter {
            width,
            long_seqids,
            previous_whole_ids: HashSet::new(),
            header: Vec::new(),
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
        let mut header = std::mem::take(&mut self.header);
        header.clear();
        // NCBI's insertion of the definition lines into m_Out (the long header
        // of `CFastaOstream::Write` is taken as one too, LOSAT's form of it).
        if self.long_seqids {
            // The Bioseq's first source, from every line the filter keeps (the target's first),
            // names the organism of a title made from it: the long header's alone.
            let source_taxid = set.first_source_taxid(target_seqid);
            self.long_header(source, shown, source_taxid, range, config, &mut header)?;
        } else {
            short_header(shown, range, config, &mut header)?;
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
        let width = self.width;
        let written = out.with_inserter(|ins| {
            if !header.is_empty() {
                ins.insert(&header)?;
            }
            write_sequence_direct(
                &letters,
                width,
                || match &masks {
                    Some(masks) => mask_changes(masks, range, letters.len(), config.strand_minus),
                    None => Vec::new(),
                },
                ins,
            )
        });
        self.header = header;
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

// NCBI reference (598d8ae6): c++/src/objmgr/util/sequence.cpp:3164-3320
/// The sequence lines of `CFastaOstream::x_WriteSequence`: lines of `width` letters, the
/// last one shorter, each ended by a line feed, written straight into `out` (`StreamInserter`).
/// Whenever whole wrapped lines fit strictly inside `out.room()` (so no NCBI write can happen
/// across them), they are wrapped directly into NCBI's held buffer; at NCBI write boundaries
/// (and on terminals or cut-recording test buffers), each piece split at `SEQ_VECTOR_CACHE` or
/// a masking-state change (`changes`) is one insertion and each line feed is a `put`.
fn write_sequence_direct(
    letters: &[u8],
    width: u32,
    changes: impl FnOnce() -> Vec<usize>,
    out: &mut StreamInserter<'_>,
) -> std::io::Result<()> {
    let width_u32 = width;
    let width = width.max(1) as usize;
    let len = letters.len();
    if len == 0 {
        return Ok(());
    }
    if out.room() == usize::MAX {
        return write_sequence_unbuffered(letters, width, out);
    }
    let mut changes_opt: Option<Vec<usize>> = None;
    let mut make_changes = Some(changes);
    let mut changes_idx = 0;
    let mut rem_line = width;
    let mut pos = 0;
    while pos < len {
        let room = out.room();
        if rem_line == width && room > 1 {
            let avail = room - 1;
            let rem_letters = len - pos;
            let full_lines = rem_letters / width;
            let take_lines = (avail / (width + 1)).min(full_lines);
            let tail_letters = rem_letters - take_lines * width;
            let tail_bytes = if tail_letters > 0 {
                tail_letters + 1
            } else {
                0
            };
            let take_tail = take_lines == full_lines
                && tail_bytes > 0
                && take_lines * (width + 1) + tail_bytes <= avail;
            let batch_bytes = take_lines * (width + 1) + if take_tail { tail_bytes } else { 0 };
            if batch_bytes > 0 {
                let batch_letters = take_lines * width + if take_tail { tail_letters } else { 0 };
                let slice = &letters[pos..pos + batch_letters];
                if out.push_within_room(batch_bytes, |dst| {
                    crate::sequence::write_fasta_lines(slice, width_u32, dst);
                })? {
                    pos += batch_letters;
                    continue;
                }
            }
        }
        let changes_slice =
            changes_opt.get_or_insert_with(|| make_changes.take().map_or_else(Vec::new, |f| f()));
        while changes_idx < changes_slice.len() && changes_slice[changes_idx] <= pos {
            changes_idx += 1;
        }
        let mut end = ((pos / SEQ_VECTOR_CACHE + 1) * SEQ_VECTOR_CACHE).min(len);
        if let Some(&change) = changes_slice.get(changes_idx) {
            end = end.min(change);
        }
        let step_end = end.min(pos + rem_line);
        out.insert(&letters[pos..step_end])?;
        rem_line -= step_end - pos;
        pos = step_end;
        if rem_line == 0 {
            out.put(b'\n')?;
            rem_line = width;
        }
    }
    if rem_line < width {
        out.put(b'\n')?;
    }
    Ok(())
}

/// Batched line wrapping for library `CallerOutput` streams (`room == usize::MAX`) that have no
/// NCBI write-point boundaries to preserve.
fn write_sequence_unbuffered(
    letters: &[u8],
    width: usize,
    out: &mut StreamInserter<'_>,
) -> std::io::Result<()> {
    let mut chunk = [0u8; 4096];
    let mut used = 0;
    for line in letters.chunks(width) {
        if used + line.len() + 1 > chunk.len() && used > 0 {
            out.insert(&chunk[..used])?;
            used = 0;
        }
        if line.len() < chunk.len() {
            chunk[used..used + line.len()].copy_from_slice(line);
            chunk[used + line.len()] = b'\n';
            used += line.len() + 1;
        } else {
            out.insert(line)?;
            out.put(b'\n')?;
        }
    }
    if used > 0 {
        out.insert(&chunk[..used])?;
    }
    Ok(())
}

#[cfg(test)]
fn write_sequence(
    letters: &[u8],
    width: u32,
    changes: &[usize],
    out: &mut crate::report::ostream::CutBuffer,
) {
    crate::report::ostream::with_stream_inserter(out, &mut |ins| {
        write_sequence_direct(letters, width, || changes.to_vec(), ins)
    })
    .expect("writing to CutBuffer is infallible");
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
