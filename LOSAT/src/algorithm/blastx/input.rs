//! Local FASTA state before translation. Internal IDs and defline titles differ.
use anyhow::{bail, Result};
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
#[cfg(test)]
use std::path::Path;

use crate::blastinput::fasta_reader::RecordSource;
use crate::blastinput::seq_range::{
    record_interval, RecordInterval, RecordPlacement, SequenceRange,
};
use crate::ncbi_str::{is_space, truncate_spaces, truncate_spaces_begin, truncate_spaces_end};

#[derive(Debug, Clone)]
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:619-631
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:177-180
pub struct FastaRecord {
    pub internal_id: String,
    /// The title as read, bytes: NCBI's reader keeps a title in any encoding, and its formatters
    /// make what they show of it (`report::defline`).
    pub title: Vec<u8>,
    pub sequence: Vec<u8>,
    pub lowercase_masks: Vec<(i32, i32)>,
    pub warnings: Vec<InputWarning>,
    /// The defline line as read, for the BLAST archive (outfmt 11); `None` (8 bytes, no
    /// allocation) unless the reader was asked to keep it. A record without a defline line
    /// has none either (the archive then writes `>` and the title).
    pub source: Option<Box<crate::blastinput::fasta_reader::RecordSource>>,
    /// The record's letters on each side of the `-query_loc` interval that `cut_in_place` cut it
    /// to, where the reader keeps them (`keep_around`, `record_letter`): `None` (8 bytes, no
    /// allocation) for every other record.
    pub around: Option<Box<LettersAround>>,
    /// On the first record of a reader batch, the batch's records that the `-query_loc` range
    /// skipped (`SkippedRecord`), which NCBI's reader added to its scope, where the reader keeps
    /// them (`FastaStream::keep_skipped`); `None` (no allocation) for every other record.
    pub skipped: Option<Box<[SkippedRecord]>>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:146-155 (cited at
// `emit_fasta_batch`: the batch drops a record whose range starts past its end) after the
// reader added it to its scope (blast_fasta_input.cpp:395, cited at
// `query_batch::QueryScope::of_shown`)
/// A record of a reader batch that the `-query_loc` range skipped (`FastaRecord::skipped`): what
/// NCBI's scope holds of it at the batch's prefetch, its generated local id (`Query_<n>`, its
/// number of the run) and the Seq-ids that `-parse_deflines` read.
#[derive(Debug, Clone)]
pub struct SkippedRecord {
    pub internal_id: String,
    pub believed: Option<std::sync::Arc<crate::blastinput::fasta_reader::believed::BelievedIds>>,
    /// The letters of the whole record, by which NCBI's reader counts a later record with its id
    /// (`query_batch::BatchCount`).
    pub length: usize,
    /// The records of the batch that the reader kept before it (its place in the batch).
    pub kept_before: usize,
}

impl SkippedRecord {
    /// What NCBI's scope holds of `record`, which the range skips after `kept_before` records of
    /// its batch.
    pub fn of(record: &FastaRecord, kept_before: usize) -> Self {
        Self {
            internal_id: record.internal_id.clone(),
            believed: record
                .source
                .as_ref()
                .and_then(|source| source.believed.clone()),
            length: record.sequence.len(),
            kept_before,
        }
    }
}

/// The records of a reader batch `records` in the order NCBI's reader read them, the records
/// that the range skipped (`FastaRecord::skipped`) among them: the ids that `-parse_deflines`
/// read, the letters of the whole record, and the record's place in `records` (`None` for a
/// skipped record). NCBI's scope holds every record read, so the scopes of the queries
/// (`query_batch::QueryScopes::add_batch_ids`) take them all (line WQA D10 round 6, 2026-11-20:
/// NCBI BLAST+ 2.17.0 counts and sets up a `q1` after a `q1` that `-query_loc` skips by the
/// skipped record's letters, `r6/probe/skip2`, `xskip`).
pub fn scope_records(
    records: &[FastaRecord],
) -> Vec<(
    Option<&crate::blastinput::fasta_reader::believed::BelievedIds>,
    usize,
    Option<usize>,
)> {
    let skipped: &[SkippedRecord] = records
        .first()
        .and_then(|first| first.skipped.as_deref())
        .unwrap_or_default();
    let mut skipped = skipped.iter().peekable();
    let mut ordered = Vec::with_capacity(records.len());
    for place in 0..=records.len() {
        while let Some(record) = skipped.next_if(|record| record.kept_before <= place) {
            ordered.push((record.believed.as_deref(), record.length, None));
        }
        if let Some(record) = records.get(place) {
            ordered.push((record.believed(), record.placement().length, Some(place)));
        }
    }
    ordered
}

/// The local ids and believed Seq-ids (`QueryScope::of_shown`) of the query records of a reader
/// batch `records` as NCBI's reader added them to its scope: the records, then those that the
/// range skipped (`FastaRecord::skipped`).
pub fn shown_ids(
    records: &[FastaRecord],
) -> impl Iterator<
    Item = (
        &str,
        Option<&crate::blastinput::fasta_reader::believed::BelievedIds>,
    ),
> + '_ {
    records
        .iter()
        .map(|record| (record.internal_id.as_str(), record.believed()))
        .chain(records.iter().flat_map(|record| {
            record.skipped.iter().flat_map(|skipped| {
                skipped
                    .iter()
                    .map(|skipped| (skipped.internal_id.as_str(), skipped.believed.as_deref()))
            })
        }))
}

/// Attaches the records of a reader batch that the range skipped (`skipped`, taken) to the
/// batch's first record `records` (`FastaRecord::skipped`); a batch of no record (NCBI's query
/// factory refuses it, `Empty CBlastQueryVector`) drops them.
pub fn attach_skipped(records: &mut [FastaRecord], skipped: &mut Vec<SkippedRecord>) {
    if skipped.is_empty() {
        return;
    }
    let skipped = std::mem::take(skipped);
    if let Some(first) = records.first_mut() {
        first.skipped = Some(skipped.into_boxed_slice());
    }
}

/// The letters of a record on each side of the interval that a cut keeps (`FastaRecord::keep_around`):
/// `before` ends at the interval's first letter, `after` starts past its last.
#[derive(Debug, Clone)]
pub struct LettersAround {
    pub before: Vec<u8>,
    pub after: Vec<u8>,
}

impl FastaRecord {
    /// The Seq-ids `-parse_deflines` read, `None` without the option
    /// (`fasta_reader::FastaRecord::believed`).
    pub fn believed(&self) -> Option<&crate::blastinput::fasta_reader::believed::BelievedIds> {
        self.source.as_ref()?.believed.as_deref()
    }

    /// The record's names (`fasta_reader::RecordNames`, the one owner of their rules for every
    /// program's records): its local ID (`internal_id`), its title and the Seq-ids
    /// `-parse_deflines` read.
    pub fn names(&self) -> crate::blastinput::fasta_reader::RecordNames<'_> {
        crate::blastinput::fasta_reader::RecordNames {
            local_id: &self.internal_id,
            title: &self.title,
            believed: self.believed(),
        }
    }

    /// What `Query=` and `# Query:` show (`RecordNames::query_line`).
    pub fn query_line(&self) -> std::borrow::Cow<'_, [u8]> {
        self.names().query_line()
    }

    /// The id of the record in NCBI's warnings (`RecordNames::warning_id`).
    pub fn warning_id(&self) -> &[u8] {
        self.names().warning_id()
    }

    /// The Seq-id of the record in the Seq-aligns and the archive (`RecordNames::seq_id`).
    pub fn seq_id(&self) -> crate::report::seqid::SeqId<'_> {
        self.names().seq_id()
    }

    /// The ID that the tabular formats and the reports show (`RecordNames::shown_id`): the first
    /// word of the title, the local ID without a title, with `-parse_deflines` the `qaccver` of
    /// the record's ids.
    pub fn shown_id(&self) -> &[u8] {
        self.names().shown_id()
    }

    /// The tabular subject id (`sseqid`, `sacc`, `saccver`) of a `-subject` record read without
    /// `-parse_deflines` (`RecordNames::subject_tabular_id` of a protein record): `None` where
    /// `GenerateDefline` throws (`report::tabular::SubjectLookup::Failed`).
    pub fn subject_tabular_id(&self) -> Option<Vec<u8>> {
        // A record with parsed ids has its tabular ids in `report::RecordsParsed`, made (or refused)
        // when the subjects are read (`RecordsParsed::of`), and the reports read them there first:
        // they never ask this of one. Without parsed ids the shared rule has no refusal (`Err`).
        debug_assert!(
            self.believed().is_none(),
            "a record with parsed ids has its tabular ids in RecordsParsed"
        );
        self.names().subject_tabular_id(true).ok().flatten()
    }

    /// Where the letters lie in the record: the start of the `-query_loc` / `-subject_loc`
    /// interval that `cut_in_place` cut the record to, and the record's length (0 and the
    /// letters' length for a record searched whole).
    pub fn placement(&self) -> RecordPlacement {
        match self.source.as_deref() {
            Some(RecordSource {
                cut_from,
                uncut_length,
                ..
            }) if *uncut_length > 0 => RecordPlacement {
                offset: *cut_from,
                length: *uncut_length,
            },
            _ => RecordPlacement {
                offset: 0,
                length: self.sequence.len(),
            },
        }
    }

    /// The letters of the record as the archive's Bioseq holds them: the whole record when a range
    /// cut it keeping them (`cut_in_place`'s `keep_whole`), else the searched letters.
    pub fn archive_letters(&self) -> &[u8] {
        self.source
            .as_deref()
            .and_then(|source| source.uncut.as_deref())
            .unwrap_or(&self.sequence)
    }

    /// The record's lower-case runs `[left, right]` as the reader found them: in the whole record
    /// for a record that keeps its whole letters (those of `archive_letters`), else in the searched
    /// letters (`lowercase_masks`).
    pub fn archive_lowercase_runs(&self) -> Vec<(i32, i32)> {
        match self
            .source
            .as_deref()
            .and_then(|source| source.uncut.as_deref())
        {
            Some(uncut) => crate::algorithm::common::masks::lowercase_runs(uncut),
            None => self.lowercase_masks.clone(),
        }
    }

    /// Whether the record has lower-case runs (`archive_lowercase_runs`).
    pub fn has_lowercase_runs(&self) -> bool {
        match self
            .source
            .as_deref()
            .and_then(|source| source.uncut.as_deref())
        {
            Some(uncut) => uncut.iter().any(u8::is_ascii_lowercase),
            None => !self.lowercase_masks.is_empty(),
        }
    }

    /// The record as the Seq-loc of the archive's `seq-loc-list` (`-parse_deflines`, an id that is
    /// not local): the searched interval (`fasta_reader::FastaRecord::archive_seq_loc`).
    pub fn archive_seq_loc(&self) -> Option<crate::report::archive::bioseq::QuerySeqLoc<'_>> {
        let from = self.placement().offset as i32;
        self.believed()?
            .seq_loc(from, from + self.sequence.len() as i32 - 1)
    }

    /// Whether the letter before `from` lies in a lower-case run. NCBI restricts the runs of the
    /// whole record to each chunk of a split query (`split::x_chunk_masks`): of the letters
    /// outside an interval that starts at `from`, only that one can reach a chunk's masks (the
    /// first chunk's first letter), so a record cut there keeps its whole letters only then.
    pub fn lower_case_before(&self, from: usize) -> bool {
        let Some(letter) = from.checked_sub(1) else {
            return false;
        };
        self.lowercase_masks
            .iter()
            .any(|&(left, right)| left as usize <= letter && letter <= right as usize)
    }

    /// Keeps the record's letters within `around` of its interval `[from, to_exclusive)` on each
    /// side, which `cut_in_place` then drops (nothing for `around` 0): the letters outside a
    /// `-query_loc` interval that `-subject`'s Seq-loc mapper reads into a query row
    /// (`empty_rows::query_letters_around`), apart from the whole letters of `cut_in_place`'s
    /// `keep_whole` (which the archive, the lower case of the rows and a split query's masks read).
    pub fn keep_around(&mut self, from: usize, to_exclusive: usize, around: usize) {
        if around == 0 {
            return;
        }
        let before = &self.sequence[from.saturating_sub(around)..from];
        let after = &self.sequence[to_exclusive..(to_exclusive + around).min(self.sequence.len())];
        if before.is_empty() && after.is_empty() {
            return;
        }
        self.around = Some(Box::new(LettersAround {
            before: before.to_vec(),
            after: after.to_vec(),
        }));
    }

    /// The record's letter `at` (in the record), where the record holds it: a searched letter, or
    /// one that its cut kept around the interval (`keep_around`).
    pub fn record_letter(&self, at: usize) -> Option<u8> {
        let offset = self.placement().offset;
        let around = self.around.as_deref();
        match at.checked_sub(offset) {
            Some(inside) if inside < self.sequence.len() => Some(self.sequence[inside]),
            Some(inside) => around?.after.get(inside - self.sequence.len()).copied(),
            None => {
                let before = &around?.before;
                before.get(before.len().checked_sub(offset - at)?).copied()
            }
        }
    }

    /// Cuts the record to its letters `[from, to_exclusive)` in place, as NCBI searches a record
    /// with a range (`seq_range::record_cut`): the lower-case masks are clipped to the letters
    /// and moved with them, and the source keeps where they lie in the record (`placement`).
    /// `keep_whole`: the source keeps the record's whole letters too (`archive_letters`,
    /// `archive_lowercase_runs`), for the archive's Bioseq (outfmt 11), the lower case of the
    /// whole record that outfmt 0 shows on a `-subject_loc` subject's rows, and the run of the
    /// letter before the interval that the chunks of a split query restrict
    /// (`lower_case_before`, `split::x_chunk_masks`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1041-1046
    pub fn cut_in_place(&mut self, from: usize, to_exclusive: usize, keep_whole: bool) {
        let length = self.sequence.len();
        let source = self.source.get_or_insert_with(Default::default);
        source.uncut_length = length;
        source.cut_from = from;
        if keep_whole {
            // The searched letters are copied out and the record's letters become the whole
            // letters, so the record holds them once. As the shared reader keeps them: lower case
            // where the lower-case masks lie, which the archive's `LCaseMask` reads back
            // (`archive_lowercase_runs`).
            let mut uncut = std::mem::take(&mut self.sequence);
            self.sequence = uncut[from..to_exclusive].to_vec();
            for &(left, right) in &self.lowercase_masks {
                uncut[left as usize..=right as usize].make_ascii_lowercase();
            }
            source.uncut = Some(uncut);
        } else {
            self.sequence.truncate(to_exclusive);
            self.sequence.drain(..from);
        }
        let (from, last) = (from as i32, to_exclusive as i32 - 1);
        self.lowercase_masks.retain_mut(|(left, right)| {
            *left = (*left - from).max(0);
            *right = (*right).min(last) - from;
            *left <= *right
        });
    }
}
#[derive(Debug, Clone)]
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:974-995
pub struct InputWarning {
    pub line: usize,
    pub kind: &'static str,
    pub positions: Vec<usize>,
    /// The message of kind `defline_id` (the Seq-id warnings of `-parse_deflines`) and of kind
    /// `gap_line` (the warnings of a `>?` gap line, `fasta_reader::reader::gap_line_size`).
    pub text: Vec<u8>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:319-343
#[cfg(test)]
pub fn read_fasta(path: &Path, protein: bool, lowercase: bool) -> Result<Vec<FastaRecord>> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:291-301
    let file = anyhow::Context::with_context(std::fs::File::open(path), || {
        format!("cannot read BLASTX input {}", path.display())
    })?;
    parse_fasta_stream_with_warnings(
        &mut FastaStream::from_file(file),
        protein,
        lowercase,
        usize::MAX,
        &mut |_| Ok(()),
        &mut |_| Ok(()),
    )
}
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:350-375
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:378-389
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:941-966
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2038-2043
#[cfg(test)]
pub fn parse_fasta(data: &[u8], protein: bool, lowercase: bool) -> Result<Vec<FastaRecord>> {
    parse_fasta_batches(data, protein, lowercase, usize::MAX, &mut |_| Ok(()))
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:140-171
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:260-293
#[cfg(test)]
pub(crate) fn parse_fasta_batches(
    data: &[u8],
    protein: bool,
    lowercase: bool,
    batch_size: usize,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
) -> Result<Vec<FastaRecord>> {
    parse_fasta_batches_with_warnings(data, protein, lowercase, batch_size, deliver, &mut |_| {
        Ok(())
    })
}
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2213-2217
// Native caller posts each reader event immediately; retained warnings continue
// to serve the independent B input-state probe. Batch delivery is unchanged.
#[cfg(any(target_arch = "wasm32", test))]
pub(crate) fn parse_fasta_batches_with_warnings(
    data: &[u8],
    protein: bool,
    lowercase: bool,
    batch_size: usize,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
    warn: &mut dyn FnMut(&InputWarning) -> Result<()>,
) -> Result<Vec<FastaRecord>> {
    parse_fasta_lines(
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:302-314
        FastaMemoryLines { remaining: data },
        protein,
        lowercase,
        batch_size,
        true,
        false,
        None,
        None,
        0,
        false,
        None,
        deliver,
        warn,
    )
}
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:350-375
fn parse_fasta_lines<L: AsRef<[u8]>>(
    lines: impl Iterator<Item = L>,
    protein: bool,
    lowercase: bool,
    batch_size: usize,
    keep_records: bool,
    keep_defline: bool,
    parse_deflines: Option<BelievedRead<'_>>,
    range: Option<&SequenceRange>,
    around: usize,
    keep_skipped: bool,
    count: Option<crate::blastinput::query_batch::BatchCount>,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
    warn: &mut dyn FnMut(&InputWarning) -> Result<()>,
) -> Result<Vec<FastaRecord>> {
    let mut believed_keys = std::collections::HashSet::new();
    let mut believed_generated = 0;
    let mut believed_molecule = None;
    let mut batch = ReaderBatch {
        range,
        keep_whole: keep_defline,
        around,
        kept: keep_records.then(Vec::new),
        keep_skipped,
        count,
        ..ReaderBatch::default()
    };

    let mut records = Vec::new();
    let mut current: Option<FastaRecord> = None;
    let mut pending_title_warning = None;
    let mut raw_residues = 0;
    for (line_index, raw) in lines.enumerate() {
        let raw = raw.as_ref();
        let line = raw.strip_suffix(b"\r").unwrap_or(raw);
        // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:350-375
        let escaped_header = line
            .strip_prefix(b">?_")
            .map(|suffix| [b">".as_slice(), suffix].concat());
        let header = escaped_header.as_deref().unwrap_or(line);
        let is_header = line.starts_with(b">") && !header.starts_with(b">?");
        let line = if is_header {
            header
        } else {
            // `NStr::TruncateSpaces_Unsafe` of the next line of `GetLineReader()` (`eTrunc_Both`, fasta.cpp:376).
            truncate_spaces(line)
        };
        // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:399-408
        let rewritten_data;
        let line = if !is_header && line.starts_with(b">?_") {
            rewritten_data = [b">".as_slice(), &line[3..]].concat();
            rewritten_data.as_slice()
        } else {
            line
        };
        if line.is_empty() || matches!(line[0], b';' | b'#' | b'!') {
            continue;
        }
        // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:383-398
        // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:774-780
        // `-parse_deflines` reads without `fDLOptional`, so a first line that is not a defline,
        // a comment or blank stops the read before `ParseDataLine`: a `>?N` gap line too, and
        // whatever its size or modifiers.
        if !is_header && parse_deflines.is_some() && current.is_none() {
            return Err(NoDeflineError(line_index as u64 + 1).into());
        }
        if line.starts_with(b">?") {
            // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1098-1122
            // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:77-82
            // The gap line is read as the reader of the other programs reads it
            // (`fasta_reader::reader::gap_line_size`, the whole `ParseGapLine`): a size that is not
            // a number from 1 to 4294967295 is 1 after `Bad gap size`, and the modifiers only warn.
            // The warnings go out as the line is read, before those of the data lines after it.
            let mut messages = Vec::new();
            let count = crate::blastinput::fasta_reader::reader::gap_line_size(
                line,
                line_index as u64 + 1,
                &mut |message| {
                    messages.push(message.strip_suffix(b"\n").unwrap_or(message).to_vec());
                    Ok(())
                },
            )? as usize;
            let record = current
                .get_or_insert_with(|| new_record(batch.number(&records), protein, Vec::new()));
            for text in messages {
                let warning = InputWarning {
                    line: line_index + 1,
                    kind: "gap_line",
                    positions: Vec::new(),
                    text,
                };
                warn(&warning)?;
                record.warnings.push(warning);
            }
            // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1387-1394
            // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1411-1425
            if count > i32::MAX as usize - record.sequence.len() {
                bail!("BLASTX sequence exceeds Int4 length");
            }
            // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1082-1090
            // NCBI reference (598d8ae6): c++/include/objtools/readers/fasta.hpp:505-514
            if let Some(last) = record
                .lowercase_masks
                .last_mut()
                .filter(|r| r.1 + 1 == record.sequence.len() as i32)
            {
                last.1 += count as i32;
            }
            record.sequence.resize(
                record.sequence.len() + count,
                if protein { b'X' } else { b'N' },
            );
            continue;
        }
        if is_header {
            if let Some(record) = current.take() {
                finish_record(
                    record,
                    &mut records,
                    pending_title_warning.take(),
                    protein,
                    parse_deflines
                        .and_then(|read| read.id_scope.map(|scope| (scope, read.program.upper))),
                    warn,
                )?;
                emit_fasta_batch(&mut records, batch_size, &mut batch, deliver)?;
            }
            // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:353-359
            let believed = match parse_deflines {
                Some(read) => Some(read_believed_defline(
                    line,
                    protein,
                    read,
                    keep_defline,
                    line_index as u64 + 1,
                    &mut believed_keys,
                    &mut believed_generated,
                    warn,
                )?),
                None => None,
            };
            // The `isspace` skip before the title (`eTrunc_Begin`), then
            // `NStr::TruncateSpacesInPlace` of `processed_title` at its end (fasta.cpp:2038).
            let title_bytes = match &believed {
                Some(defline) => truncate_spaces_begin(&line[defline.title_search..]),
                None => truncate_spaces_begin(&line[1..]),
            };
            let end = title_bytes
                .iter()
                .skip(1)
                .position(|&b| b < b' ')
                .map_or(title_bytes.len(), |i| i + 1);
            let title = truncate_spaces_end(&title_bytes[..end]).to_vec();
            // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1623-1677
            let raw_title = &title_bytes[..end];
            let kind = if !protein
                && raw_title.len() > 20
                && raw_title[raw_title.len() - 20..]
                    .iter()
                    .all(|c| b"ACGTacgt".contains(c))
            {
                Some("title_nucleotide_residues")
            } else if protein
                && raw_title.len() > 50
                && raw_title[raw_title.len() - 50..]
                    .iter()
                    .all(u8::is_ascii_alphabetic)
            {
                Some("title_amino_acids")
            } else {
                None
            };
            pending_title_warning = kind.map(|kind| InputWarning {
                line: line_index + 1,
                kind,
                positions: Vec::new(),
                text: Vec::new(),
            });
            let mut record = new_record(batch.number(&records), protein, title);
            believed_molecule = None;
            if keep_defline || believed.is_some() {
                record.source = Some(Box::new(crate::blastinput::fasta_reader::RecordSource {
                    defline: if keep_defline {
                        header.to_vec()
                    } else {
                        Vec::new()
                    },
                    uncut: None,
                    believed: believed.as_ref().map(|defline| defline.ids.clone()),
                    cut_from: 0,
                    uncut_length: 0,
                }));
            }
            if let Some(defline) = believed {
                record.internal_id = defline.ids.bare_id();
                believed_molecule = defline.molecule;
            }
            current = Some(record);
            raw_residues = 0;
            continue;
        }
        let record =
            current.get_or_insert_with(|| new_record(batch.number(&records), protein, Vec::new()));
        if raw_residues == 0 {
            check_data_line(line, line_index + 1)?;
        }
        if believed_molecule.is_some_and(|named| named != protein)
            && !crate::blastinput::fasta_reader::believed::reads_the_same_in_both(line)
        {
            anyhow::bail!(
                "the Seq-id of the record on line {} names the other molecule under \
                 -parse_deflines (NCBI BLAST+ reads its letters as that molecule); this is not \
                 supported by LOSAT's BLASTX",
                line_index + 1
            );
        }
        let mut bad = Vec::new();
        let mut hyphens = false;
        for (pos, &c) in line.iter().enumerate() {
            if c == b';' {
                break;
            }
            // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:930-931
            if is_space(c) {
                continue;
            }
            if c == b'-' {
                hyphens = true;
                continue;
            }
            let upper = c.to_ascii_uppercase();
            let allowed = if protein {
                // The letters `ABCDEFGHIJKLMNOPQRSTUVWXYZ*`, tested without a table in memory:
                // `contains` on a static table first steps byte by byte to an aligned address,
                // so its cost per letter depended on where the linker placed the table.
                matches!(upper, b'A'..=b'Z' | b'*')
            } else {
                b"ABCDGHKMNRSTUVWY".contains(&upper)
            };
            if !allowed {
                bad.push(pos + 1);
                continue;
            }
            raw_residues += 1;
            let offset = record.sequence.len() as i32;
            record.sequence.push(if !protein && upper == b'U' {
                b'T'
            } else {
                upper
            });
            if lowercase && c.is_ascii_lowercase() {
                if let Some(last) = record
                    .lowercase_masks
                    .last_mut()
                    .filter(|r| r.1 + 1 == offset)
                {
                    last.1 = offset;
                } else {
                    record.lowercase_masks.push((offset, offset));
                }
            }
        }
        if hyphens {
            let warning = InputWarning {
                line: line_index + 1,
                kind: "hyphens",
                positions: Vec::new(),
                text: Vec::new(),
            };
            warn(&warning)?;
            record.warnings.push(warning);
        }
        if !bad.is_empty() {
            let warning = InputWarning {
                line: line_index + 1,
                kind: "invalid_residues",
                positions: bad,
                text: crate::blastinput::fasta_reader::believed::nuc_or_prot(believed_molecule)
                    .as_bytes()
                    .to_vec(),
            };
            warn(&warning)?;
            record.warnings.push(warning);
        }
    }
    if let Some(record) = current {
        finish_record(
            record,
            &mut records,
            pending_title_warning.take(),
            protein,
            parse_deflines.and_then(|read| read.id_scope.map(|scope| (scope, read.program.upper))),
            warn,
        )?;
        emit_fasta_batch(&mut records, batch_size, &mut batch, deliver)?;
    }
    // A last batch of records that the range skipped is delivered empty (NCBI's query factory
    // then refuses it).
    if batch.read > 0 {
        attach_skipped(&mut records, &mut batch.skipped_records);
        deliver(&records)?;
    }
    batch.release(&mut records);
    Ok(batch.kept.unwrap_or_default())
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:364-367
/// How `-parse_deflines` reads: the strand of a query's Seq-loc, and the search's ids of the other
/// role's records (`FastaStream::id_scope`).
#[derive(Clone, Copy)]
struct BelievedRead<'a> {
    strand: crate::report::archive::bioseq::SeqLocStrand,
    id_scope: Option<&'a crate::seqid::id_scope::IdScope>,
    /// The program of the reader's messages (`FastaStream::program`).
    program: super::args::ProgramName,
}

/// A defline under `-parse_deflines` (`fasta_reader::believed`).
struct BelievedLine {
    title_search: usize,
    ids: std::sync::Arc<crate::blastinput::fasta_reader::believed::BelievedIds>,
    molecule: Option<bool>,
}

/// `fasta_reader::reader::parse_believed_def_line` for BLASTX's reader.
fn read_believed_defline(
    line: &[u8],
    protein: bool,
    read: BelievedRead<'_>,
    keep_defline: bool,
    line_number: u64,
    keys: &mut std::collections::HashSet<u64>,
    generated: &mut usize,
    warn: &mut dyn FnMut(&InputWarning) -> Result<()>,
) -> Result<BelievedLine> {
    use crate::blastinput::fasta_reader::believed::{
        duplicate_id, generated_ids, named_molecule, parse_believed_defline,
        refuse_range_in_archive, refused_id, BelievedIds,
    };
    let mut written = Ok(());
    let parsed = parse_believed_defline(
        line,
        protein,
        read.program.lower,
        line_number,
        &mut |message| {
            if written.is_ok() {
                written = warn(&InputWarning {
                    line: line_number as usize,
                    kind: "defline_id",
                    positions: Vec::new(),
                    text: message.to_vec(),
                });
            }
        },
    );
    written?;
    let parsed = parsed.map_err(FastaParseError)?;
    let role = if protein { "subject" } else { "query" };
    if parsed.range && keep_defline {
        return Err(refuse_range_in_archive(role, read.program.upper));
    }
    let molecule = named_molecule(&parsed.ids);
    let ids = if parsed.ids.is_empty() || parsed.range {
        *generated += 1;
        generated_ids(&format!(
            "{}_{generated}",
            if protein { "Subject" } else { "Query" }
        ))
    } else {
        parsed.ids
    };
    // Two queries with one id are the scopes' of the query batches (`QueryScopes`, which the
    // search adds as the reader delivers the batches): NCBI's formatter empties the scope after
    // each batch against a database.
    let refusal = refused_id(&ids, role, "-parse_deflines", read.program.upper).or_else(|| {
        protein
            .then(|| duplicate_id(&ids, keys, role, read.program.upper))
            .flatten()
    });
    if let Some(refusal) = refusal {
        return Err(refusal);
    }
    Ok(BelievedLine {
        title_search: parsed.title_search,
        ids: std::sync::Arc::new(BelievedIds::new(ids, read.strand)),
        molecule,
    })
}

fn new_record(index: usize, protein: bool, title: Vec<u8>) -> FastaRecord {
    FastaRecord {
        internal_id: format!(
            "{}_{n}",
            if protein { "Subject" } else { "Query" },
            n = index + 1
        ),
        title,
        sequence: Vec::new(),
        lowercase_masks: Vec::new(),
        warnings: Vec::new(),
        source: None,
        around: None,
        skipped: None,
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:437-447
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:457-461
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:1361-1365
/// NCBI posts a record's title warning when it assembles the record, after the warnings of its
/// data lines; a record whose data stops the reading posts none. Under `-parse_deflines` the
/// record's ids then join the search's ids (`fasta_reader`'s `assemble_seq`, whose
/// `add_to_id_scope` follows `parse_title`), so LOSAT's refusal of a record that a record of the
/// other role names comes after the record's messages.
fn finish_record(
    mut record: FastaRecord,
    records: &mut Vec<FastaRecord>,
    title_warning: Option<InputWarning>,
    protein: bool,
    id_scope: Option<(&crate::seqid::id_scope::IdScope, &str)>,
    warn: &mut dyn FnMut(&InputWarning) -> Result<()>,
) -> Result<()> {
    if let Some(warning) = title_warning {
        warn(&warning)?;
        record.warnings.push(warning);
    }
    // A nucleotide query and a protein subject are never one record (`RecordMark`).
    if let (Some((scope, program)), Some(believed)) = (id_scope, record.believed()) {
        let role = if protein { "subject" } else { "query" };
        scope.add_record(&believed.ids, None, role, program)?;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:334-340
    // Empty records remain input objects; query/subject setup owns their
    if record.sequence.len() > i32::MAX as usize {
        bail!("BLASTX sequence exceeds Int4 length");
    }
    records.push(record);
    Ok(())
}
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:726-763
fn check_data_line(line: &[u8], number: usize) -> Result<()> {
    let len = line.len().min(70);
    let (mut good, mut bad) = (0, 0);
    for &c in &line[..len] {
        if c.is_ascii_alphabetic() || c == b'*' {
            good += 1;
        } else if c == b';' {
            break;
        } else if c != b'-' && !is_space(c) && !c.is_ascii_digit() {
            bad += 1;
        }
    }
    if bad >= good / 3 && (len > 3 || good == 0 || bad > good) {
        return Err(FastaParseError(format!("CFastaReader: Near line {number}, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.")).into());
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/include/objtools/readers/reader_exception.hpp:46-49
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:754-758
#[derive(Debug)]
pub(crate) struct FastaParseError(pub String);
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:181-184
impl std::fmt::Display for FastaParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
// NCBI reference (598d8ae6): c++/include/objtools/readers/reader_exception.hpp:46-49
impl std::error::Error for FastaParseError {}

/// The reader's `eNoDefline` at a line (1-based): under `-parse_deflines` the first line that
/// is neither blank, a comment nor a defline. Its text is the reader's (`CFastaReader: ...`, a
/// query's `BLAST query error`, `native::input_error`); the subjects end with `GetAllSeqs`'s
/// error instead (`native::subject_input_error`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:391-396
#[derive(Debug)]
pub(crate) struct NoDeflineError(pub u64);

impl std::fmt::Display for NoDeflineError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "CFastaReader: Input doesn't start with a defline or comment around line {}",
            self.0
        )
    }
}

impl std::error::Error for NoDeflineError {}

#[cfg(test)]
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:86-94
#[path = "input_tests.rs"]
mod tests;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:140-171
///
/// With a range (`-query_loc`) the batch counts the length of the whole record and cuts the record
/// to its interval; a record whose interval starts more than one letter past its end adds nothing
/// and is not kept (it keeps its number: `ReaderBatch::skipped`). Under `-parse_deflines` a
/// record counts as NCBI's reader counts it (`ReaderBatch::count`, `query_batch::BatchCount`: a
/// repeated id by its first record's letters).
fn emit_fasta_batch(
    records: &mut Vec<FastaRecord>,
    batch_size: usize,
    batch: &mut ReaderBatch<'_>,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
) -> Result<()> {
    let record = records.last_mut().expect("completed FASTA record");
    let length = record.sequence.len();
    batch.read += 1;
    let interval = batch.range.map(|range| record_interval(range, length));
    let skipped = matches!(interval, Some(RecordInterval::PastEnd));
    batch.length += match batch.count.as_mut() {
        Some(count) => count.count(record.believed(), length, skipped),
        None if skipped => 0,
        None => length,
    };
    match interval {
        None => {}
        Some(RecordInterval::Letters { from, to }) => {
            // The letter before the interval is the only one outside it that a split query's
            // chunks read (`split::x_chunk_masks`); a record keeps its whole letters for it.
            let keep_whole = batch.keep_whole || record.lower_case_before(from);
            record.keep_around(from, to + 1, batch.around);
            record.cut_in_place(from, to + 1, keep_whole);
        }
        Some(RecordInterval::Empty) => {
            record.cut_in_place(length, length, batch.keep_whole);
        }
        Some(RecordInterval::PastEnd) => {
            let skipped = records.pop().expect("completed FASTA record");
            batch.skipped += 1;
            if batch.keep_skipped {
                batch
                    .skipped_records
                    .push(SkippedRecord::of(&skipped, records.len()));
            }
        }
    }
    if batch.length >= batch_size {
        attach_skipped(records, &mut batch.skipped_records);
        deliver(records)?;
        // NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-261
        // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:931-938
        // NCBI's batch is local to one turn of the loop: the batch read for the next turn
        // replaces it. A read that keeps no record drops the batch once delivered likewise.
        batch.release(records);
        batch.length = 0;
        batch.read = 0;
        if let Some(count) = batch.count.as_mut() {
            count.end_batch();
        }
    }
    Ok(())
}

/// The reader batch that `emit_fasta_batch` fills.
#[derive(Default)]
struct ReaderBatch<'r> {
    /// The range that cuts each record (`-query_loc`), if any.
    range: Option<&'r SequenceRange>,
    /// Whether a cut record keeps its whole letters (the archive, `FastaRecord::cut_in_place`).
    keep_whole: bool,
    /// The letters on each side of its interval that a cut record keeps
    /// (`FastaRecord::keep_around`).
    around: usize,
    /// The records of the batches delivered, when the reader keeps them (and returns them all);
    /// `None` drops each batch once delivered, so that the reader holds one batch at a time.
    kept: Option<Vec<FastaRecord>>,
    /// The records of the batches delivered, kept or dropped.
    released: usize,
    /// The letters the batch has read (`size_read`).
    length: usize,
    /// The records the batch has read, kept or not.
    read: usize,
    /// The records the range skipped so far: NCBI's reader numbers them (`Query_<n>`).
    skipped: usize,
    /// The records of the batch that the range skipped, for its first record
    /// (`FastaRecord::skipped`), where they are kept (`keep_skipped`).
    skipped_records: Vec<SkippedRecord>,
    /// Whether the batch keeps the records that the range skipped (`FastaStream::keep_skipped`).
    keep_skipped: bool,
    /// Under `-parse_deflines`, NCBI's count of the records (`FastaStream::batch_count`).
    count: Option<crate::blastinput::query_batch::BatchCount>,
}

impl ReaderBatch<'_> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:364-367
    /// The index of the next record (its id is `Query_` and that index plus 1): the reader's one
    /// id generator counts every record read before it, kept, dropped once delivered, or skipped
    /// by the range.
    fn number(&self, records: &[FastaRecord]) -> usize {
        self.released + records.len() + self.skipped
    }

    /// Ends a delivered batch (`records`, the batch's records alone): they join the kept
    /// records, or are dropped.
    fn release(&mut self, records: &mut Vec<FastaRecord>) {
        self.released += records.len();
        match &mut self.kept {
            // The first batch (the only one of a `-subject` file) moves without a copy.
            Some(kept) if kept.is_empty() => std::mem::swap(kept, records),
            Some(kept) => kept.append(records),
            None => records.clear(),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:310-313
#[cfg(any(target_arch = "wasm32", test))]
struct FastaMemoryLines<'a> {
    remaining: &'a [u8],
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:353-381
#[cfg(any(target_arch = "wasm32", test))]
impl<'a> Iterator for FastaMemoryLines<'a> {
    type Item = &'a [u8];
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:353-381
    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining.is_empty() {
            return None;
        }
        let end = self
            .remaining
            .iter()
            .position(|&c| c == b'\r' || c == b'\n')
            .unwrap_or(self.remaining.len());
        let line = &self.remaining[..end];
        let next = if self.remaining.get(end..end + 2) == Some(b"\r\n") {
            end + 2
        } else if end < self.remaining.len() {
            end + 1
        } else {
            end
        };
        self.remaining = &self.remaining[next..];
        Some(line)
    }
}
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-873
pub(crate) fn stream_is_empty<R: Read + Seek>(stream: &mut FastaStream<R>) -> bool {
    if stream.no_position {
        return false;
    }
    let input = &mut stream.input;
    let position = match input.stream_position() {
        Ok(position) => position,
        Err(_) => return false,
    };
    let mut byte = [0];
    loop {
        match input.read(&mut byte) {
            Ok(0) | Err(_) => return true,
            Ok(_) if is_space(byte[0]) => {}
            Ok(_) => {
                // IsIStreamEmpty restores position only after successful extraction.
                // NCBI unconditionally clears seekg failure and restores orig_state.
                let _ = input.seek(SeekFrom::Start(position));
                return false;
            }
        }
    }
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:166-174
#[derive(Clone, Copy, PartialEq)]
enum EolStyle {
    Unknown,
    Cr,
    Lf,
    CrLf,
    Mixed,
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:81-89
// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:437-443
// NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:134-143
struct PushbackBuffer {
    bytes: Vec<u8>,
    position: usize,
    capacity: usize,
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:81-89
pub(crate) struct FastaStream<R: Read> {
    input: BufReader<R>,
    failed: bool,
    eol: EolStyle,
    pushback: Vec<PushbackBuffer>,
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
    backend_available: Option<fn(&mut R) -> std::io::Result<usize>>,
    /// Whether the input's buffered bytes count as readable (`in_avail` of a `filebuf`'s get
    /// area); false for `cin`, whose `stdio_sync_filebuf` has no get area
    /// (`fasta_reader::stream::InputBackend::of`).
    get_area: bool,
    /// Whether the input has no stream position, as a pipe has none (a negative `tellg()` in
    /// `IsIStreamEmpty`; `input_files::has_stream_position`): then it is never empty.
    no_position: bool,
    /// `-parse_deflines` (`fasta_reader::believed`).
    pub(crate) parse_deflines: bool,
    /// `-strand`, read with `-parse_deflines`.
    pub(crate) strand: crate::report::archive::bioseq::SeqLocStrand,
    /// The search whose records this stream reads under `-parse_deflines`.
    pub(crate) id_scope: Option<std::sync::Arc<crate::seqid::id_scope::IdScope>>,
    /// With NCBI's data loaders on (`ncbi_environment::ApplicationSettings::data_loaders`), the
    /// role ("query", "subject") named in LOSAT's refusal of a first line that NCBI tries as a
    /// Seq-id (`fasta_reader::seq_id::reject_seq_id_line`); `None` reads every first line as
    /// FASTA, as NCBI does with the data loaders off.
    pub(crate) seq_id_role: Option<&'static str>,
    /// The program of the reader's messages (`args::ProgramName`; BLASTX unless the run is
    /// RPSTBLASTN's, `QueryStreamSettings::program`).
    pub(crate) program: super::args::ProgramName,
    /// `CStreamLineReader::UngetLine`: the line that the next read gives again (the Seq-id step
    /// of `read_fasta_stream` read it and handed it back).
    unget_line: Option<Vec<u8>>,
    /// The letters on each side of its `-query_loc` interval that each record cut to it keeps
    /// (`FastaRecord::keep_around`, `empty_rows::query_letters_around`); 0 for none.
    pub(crate) around: usize,
    /// Whether the reader keeps the records that the `-query_loc` range skips
    /// (`FastaRecord::skipped`): NCBI's scope of a batch's prefetch holds them, which only the
    /// stop of a prefetch asks (`common::database::PrefetchQueries`), so the caller sets it where
    /// a prefetch of the run can stop (`BlastxDatabase::prefetch_can_stop`; review 5 of line DBID
    /// loadord, L2). Off by default.
    pub(crate) keep_skipped: bool,
    /// Under `-parse_deflines`, NCBI's count of the records into the reader batches
    /// (`query_batch::BatchCount`, `QueryStreamSettings::batch_count`); `None` counts each record
    /// by its own letters.
    pub(crate) batch_count: Option<crate::blastinput::query_batch::BatchCount>,
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:81-89
impl<R: Read> FastaStream<R> {
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:81-89
    pub(crate) fn new(input: R) -> Self {
        Self {
            // NCBI reference (598d8ae6): c++/include/corelib/ncbistre.hpp:437-440
            // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236
            // The pinned CNcbiIfstream exposes 8191 readable bytes per default
            // window (independent C++ file_window oracle; tested at EOF/pushback
            // boundaries). Preserve that window because reuse retains EOF state.
            input: BufReader::with_capacity(8191, input),
            failed: false,
            eol: EolStyle::Unknown,
            pushback: Vec::new(),
            // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:235-236
            backend_available: None,
            get_area: true,
            no_position: false,
            parse_deflines: false,
            id_scope: None,
            strand: crate::report::archive::bioseq::SeqLocStrand::Both,
            seq_id_role: None,
            program: super::args::ProgramName::BLASTX,
            unget_line: None,
            around: 0,
            keep_skipped: false,
            batch_count: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-873 (cited at the read
    // before the search, `blast_xml::prescan_queries`: NCBI reads ahead on the argument's one stream
    // and puts it back at `orig_p`)
    /// The input this stream reads, given back: the bytes it buffered are dropped, so the caller
    /// puts the input's position where its next reading starts (`blast_xml::prescan_queries`, which
    /// reads the query file before the search and hands the same open file to the search).
    pub(crate) fn into_input(self) -> R {
        self.input.into_inner()
    }
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:306-358
    fn raw_peek(&mut self) -> Option<u8> {
        while let Some(buffer) = self.pushback.last() {
            if let Some(&byte) = buffer.bytes.get(buffer.position) {
                return Some(byte);
            }
            if self.pushback.len() > 1 {
                let below = self.pushback.remove(self.pushback.len() - 2);
                if below.position < below.bytes.len() {
                    *self.pushback.last_mut().unwrap() = below;
                }
                continue;
            }
            // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:306-338
            let buffered = if self.get_area {
                self.input.buffer().len()
            } else {
                0
            };
            let available = if buffered != 0 {
                buffered
            } else if let Some(in_avail) = self.backend_available {
                in_avail(self.input.get_mut()).unwrap_or(0)
            } else {
                0
            };
            let capacity = self.pushback.last().unwrap().capacity.max(4096);
            let requested = capacity.min(available.max(1));
            let mut bytes = vec![0; requested];
            let mut length = 0;
            while length < requested {
                match self.input.read(&mut bytes[length..]) {
                    Ok(0) | Err(_) => break,
                    Ok(n) => length += n,
                }
            }
            if length == 0 {
                return None;
            }
            bytes.truncate(length);
            let buffer = self.pushback.last_mut().unwrap();
            buffer.capacity = capacity;
            buffer.bytes = bytes;
            buffer.position = 0;
        }
        match self.input.fill_buf() {
            Ok(bytes) => bytes.first().copied(),
            Err(_) => None,
        }
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:87-98
    fn raw_take(&mut self) -> Option<u8> {
        let byte = self.raw_peek()?;
        if let Some(buffer) = self.pushback.last_mut() {
            buffer.position += 1;
        } else {
            self.input.consume(1);
        }
        Some(byte)
    }
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:413-443
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:134-143
    fn push_back(&mut self, bytes: &[u8]) {
        let mut remaining = bytes.len();
        if let Some(buffer) = self.pushback.last_mut() {
            if bytes.len() <= (4096 >> 4) {
                let take = buffer.position.min(remaining);
                buffer.position -= take;
                remaining -= take;
                buffer.bytes[buffer.position..buffer.position + take]
                    .copy_from_slice(&bytes[remaining..]);
            }
        }
        if remaining != 0 {
            self.pushback.push(PushbackBuffer {
                bytes: bytes[..remaining].to_vec(),
                position: 0,
                capacity: remaining,
            });
            // Installing a new streambuf clears state; recycling the old one does not.
            self.failed = false;
        }
    }
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
    pub(crate) fn at_eof(&mut self) -> bool {
        if self.unget_line.is_some() {
            return false;
        }
        if self.failed {
            return true;
        }
        if self.raw_peek().is_none() {
            self.failed = true;
        }
        self.failed
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:151-166
    fn next_byte(&mut self) -> Option<u8> {
        if self.failed {
            return None;
        }
        let byte = self.raw_take();
        if byte.is_none() {
            self.failed = true;
        }
        byte
    }
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:87-105
    // NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:151-174
    fn getline(&mut self, delimiters: &[u8]) -> (Vec<u8>, Option<u8>) {
        let mut line = Vec::new();
        let mut last = None;
        loop {
            let byte = if delimiters.len() == 1 {
                self.next_byte()
            } else {
                self.raw_take()
            };
            let Some(byte) = byte else {
                self.failed = true;
                break;
            };
            last = Some(byte);
            if let Some(position) = delimiters.iter().position(|&delimiter| delimiter == byte) {
                if delimiters.len() > 1 {
                    if let Some(next) = self.raw_peek() {
                        if delimiters[position + 1..].contains(&next) {
                            last = self.raw_take();
                        }
                    }
                }
                break;
            }
            line.push(byte);
        }
        (line, last)
    }
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:245-267
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:437-443
    fn advance_simple(&mut self, eol: u8, alternate: u8) -> (Vec<u8>, EolStyle) {
        let (mut line, _) = self.getline(&[eol]);
        if let Some(position) = line.iter().position(|&byte| byte == alternate) {
            let position = position + 1;
            if eol != b'\n' || position != line.len() {
                self.push_back(&line[position..]);
                self.eol = EolStyle::Mixed;
            }
            line.truncate(position - 1);
            let style = if self.eol == EolStyle::Mixed {
                EolStyle::Mixed
            } else {
                EolStyle::CrLf
            };
            return (line, style);
        }
        if eol == b'\r' && !self.failed && self.raw_peek() == Some(alternate) {
            self.raw_take();
            return (line, EolStyle::CrLf);
        }
        (
            line,
            if eol == b'\r' {
                EolStyle::Cr
            } else {
                EolStyle::Lf
            },
        )
    }
}
// NCBI reference (598d8ae6): c++/include/corelib/ncbistre.hpp:437-440
impl FastaStream<crate::blastinput::input_files::InputFile> {
    // NCBI reference (598d8ae6): c++/src/corelib/stream_utils.cpp:224-236
    /// A stream over an opened input, whose refills read as NCBI's stream of it answers them
    /// (`fasta_reader::stream::InputBackend::of`: `cin` for `-`, a `CNcbiIfstream` for a path).
    pub(crate) fn from_file(input: impl Into<crate::blastinput::input_files::InputFile>) -> Self {
        let input = input.into();
        let no_position = !crate::blastinput::input_files::has_stream_position(&input);
        let backend = crate::blastinput::fasta_reader::stream::InputBackend::of(&input);
        let mut stream = Self::new(input);
        stream.no_position = no_position;
        stream.backend_available = Some(backend.available);
        stream.get_area = backend.get_area;
        stream
    }
}

/// How the search reads its query file: `-parse_deflines`, the search's ids and `-strand`
/// (`FastaStream`'s fields). The read of the file before the search (`blast_xml::prescan_queries`)
/// reads it the same way, so that it sees the search's records and titles and stops at the reader
/// error that stops the search (under `-parse_deflines`, a first line that is not a defline).
pub(crate) struct QueryStreamSettings {
    pub(crate) parse_deflines: bool,
    pub(crate) id_scope: Option<std::sync::Arc<crate::seqid::id_scope::IdScope>>,
    pub(crate) strand: crate::report::archive::bioseq::SeqLocStrand,
    /// NCBI's data loaders (`FastaStream::seq_id_role`).
    pub(crate) data_loaders: bool,
    /// The program of the reader's messages (`FastaStream::program`).
    pub(crate) program: super::args::ProgramName,
    /// NCBI's count of the records into the batches (`FastaStream::batch_count`), from the first
    /// record on.
    pub(crate) batch_count: Option<crate::blastinput::query_batch::BatchCount>,
}

impl QueryStreamSettings {
    /// The search's stream of the query file.
    pub(crate) fn stream(
        &self,
        file: impl Into<crate::blastinput::input_files::InputFile>,
    ) -> FastaStream<crate::blastinput::input_files::InputFile> {
        let mut stream = FastaStream::from_file(file);
        stream.parse_deflines = self.parse_deflines;
        stream.id_scope = self.id_scope.clone();
        stream.strand = self.strand;
        stream.seq_id_role = self.data_loaders.then_some("query");
        stream.program = self.program;
        stream.batch_count = self.batch_count.clone();
        // The records that `-query_loc` skips join the scopes of the queries
        // (`scope_records`).
        stream.keep_skipped = self.parse_deflines;
        stream
    }

    /// A stream of the query file for a read before the search: its ids go to a copy of the
    /// search's (`IdScope::copy`), so that it stops at the record the search refuses (a query
    /// with the id of a subject) and the search reads its records into its ids as before.
    pub(crate) fn stream_before_search(
        &self,
        file: impl Into<crate::blastinput::input_files::InputFile>,
    ) -> FastaStream<crate::blastinput::input_files::InputFile> {
        let mut stream = self.stream(file);
        stream.id_scope = self
            .id_scope
            .as_deref()
            .map(|scope| std::sync::Arc::new(scope.copy()));
        stream
    }
}
// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:154-174
impl<R: Read> Iterator for FastaStream<R> {
    type Item = Vec<u8>;
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:219-241
    // NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:271-281
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(line) = self.unget_line.take() {
            return Some(line);
        }
        if self.at_eof() {
            return None;
        }
        let line = match self.eol {
            EolStyle::Unknown => {
                let (line, last) = self.getline(b"\r\n");
                // Successful unget clears eofbit before get rereads the last byte.
                if self.failed && !line.is_empty() {
                    self.failed = false;
                }
                match last {
                    Some(b'\r') => self.eol = EolStyle::Cr,
                    Some(b'\n') => self.eol = EolStyle::CrLf,
                    _ => {}
                }
                line
            }
            EolStyle::Cr => self.advance_simple(b'\r', b'\n').0,
            EolStyle::Lf => self.advance_simple(b'\n', b'\r').0,
            EolStyle::CrLf => {
                let (line, style) = self.advance_simple(b'\n', b'\r');
                if style == EolStyle::Mixed {
                    self.eol = EolStyle::Cr;
                } else if style != EolStyle::CrLf {
                    self.eol = EolStyle::Lf;
                }
                line
            }
            EolStyle::Mixed => self.getline(b"\r\n").0,
        };
        Some(line)
    }
}
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:259-262
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:140-150
#[cfg(test)]
pub(crate) fn parse_fasta_stream_with_warnings<R: Read>(
    stream: &mut FastaStream<R>,
    protein: bool,
    lowercase: bool,
    batch_size: usize,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
    warn: &mut dyn FnMut(&InputWarning) -> Result<()>,
) -> Result<Vec<FastaRecord>> {
    parse_fasta_stream_keeping(
        stream, protein, lowercase, batch_size, false, None, deliver, warn,
    )
}
/// `parse_fasta_stream_with_warnings`, keeping each record's defline line (`keep_defline`: the
/// BLAST archive, outfmt 11, writes it) and cutting each record to `range` (`-query_loc`,
/// `emit_fasta_batch`). It keeps every record it reads and returns them all (a run's `-subject`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn parse_fasta_stream_keeping<R: Read>(
    stream: &mut FastaStream<R>,
    protein: bool,
    lowercase: bool,
    batch_size: usize,
    keep_defline: bool,
    range: Option<&SequenceRange>,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
    warn: &mut dyn FnMut(&InputWarning) -> Result<()>,
) -> Result<Vec<FastaRecord>> {
    read_fasta_stream(
        stream,
        protein,
        lowercase,
        batch_size,
        true,
        keep_defline,
        range,
        deliver,
        warn,
    )
}

/// `parse_fasta_stream_keeping` without keeping the records: the reader drops each batch once it
/// has delivered it, so that a read of a query file holds one reader batch, not the file (the
/// search's queries, as NCBI's loop holds one batch: `emit_fasta_batch`; and the read of them
/// before the search, `blast_xml::prescan_queries`, which has no NCBI counterpart: NCBI reads the
/// query file once, and this read is LOSAT's, for its refusals before the first byte). The
/// batches, their records and their numbers are the same.
pub(crate) fn read_fasta_stream_batches<R: Read>(
    stream: &mut FastaStream<R>,
    protein: bool,
    lowercase: bool,
    batch_size: usize,
    keep_defline: bool,
    range: Option<&SequenceRange>,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
    warn: &mut dyn FnMut(&InputWarning) -> Result<()>,
) -> Result<()> {
    read_fasta_stream(
        stream,
        protein,
        lowercase,
        batch_size,
        false,
        keep_defline,
        range,
        deliver,
        warn,
    )
    .map(drop)
}

/// The read of `parse_fasta_stream_keeping` (`keep_records`) and `read_fasta_stream_batches`.
fn read_fasta_stream<R: Read>(
    stream: &mut FastaStream<R>,
    protein: bool,
    lowercase: bool,
    batch_size: usize,
    keep_records: bool,
    keep_defline: bool,
    range: Option<&SequenceRange>,
    deliver: &mut dyn FnMut(&[FastaRecord]) -> Result<()>,
    warn: &mut dyn FnMut(&InputWarning) -> Result<()>,
) -> Result<Vec<FastaRecord>> {
    if stream.at_eof() {
        return Ok(Vec::new());
    }
    let id_scope = stream.id_scope.clone();
    let parse_deflines = stream.parse_deflines.then_some(BelievedRead {
        strand: stream.strand,
        id_scope: id_scope.as_deref(),
        program: stream.program,
    });
    // A file without records gives one empty batch; one whose records a range skipped gave it
    // already.
    let mut delivered = false;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:130-162
    // BLASTX reads its query and its subjects with the reader of the other programs
    // (`CBlastFastaInputSource`), whose every record first tries its first line as a Seq-id when
    // the data loaders are on (`seq_id_role`). Only the input's first line can be one: a record
    // ends before the next `>` line. LOSAT refuses such a line as the other programs' reader does
    // (`fasta_reader::reader::read_one_seq`), here, before `CFastaReader`'s read (and its
    // `-parse_deflines` stop at a first line that is not a defline), which gets the line back
    // (`UngetLine`).
    if let Some(role) = stream.seq_id_role {
        if let Some(line) = stream.next() {
            if let Some(refusal) = crate::blastinput::fasta_reader::seq_id::reject_seq_id_line(
                &line,
                stream.program.upper,
                role,
            ) {
                return Err(refusal);
            }
            stream.unget_line = Some(line);
        }
    }
    let (around, keep_skipped) = (stream.around, stream.keep_skipped);
    let count = stream.batch_count.take();
    let records = parse_fasta_lines(
        stream,
        protein,
        lowercase,
        batch_size,
        keep_records,
        keep_defline,
        parse_deflines,
        range,
        around,
        keep_skipped,
        count,
        &mut |records| {
            delivered = true;
            deliver(records)
        },
        warn,
    )?;
    if !delivered {
        deliver(&[])?;
    }
    Ok(records)
}

// NCBI reference (598d8ae6): c++/src/util/line_reader.cpp:100-104
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-873
#[cfg(test)]
#[path = "input_stream_tests.rs"]
mod stream_tests;

#[cfg(test)]
#[path = "input_range_tests.rs"]
mod range_tests;

#[cfg(test)]
#[path = "input_defline_tests.rs"]
mod defline_tests;

#[cfg(test)]
#[path = "input_query_settings_tests.rs"]
mod query_settings_tests;

#[cfg(test)]
#[path = "input_standard_input_tests.rs"]
mod standard_input_tests;
