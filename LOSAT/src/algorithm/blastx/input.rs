//! BLASTX's record of a FASTA input before translation: its local ID and title, its letters in
//! upper case with the lower-case masks, and where a `-query_loc` / `-subject_loc` cut left them.
//! The records are read by the readers of the other programs (`blastinput::fasta_reader`) and
//! given this form, and cut into reader batches, by `input_source`.
use crate::blastinput::fasta_reader::RecordSource;
use crate::blastinput::seq_range::RecordPlacement;

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
    /// them (`QueryStream::keep_skipped`); `None` (no allocation) for every other record.
    pub skipped: Option<Box<[SkippedRecord]>>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:146-155 (cited at
// `BatchCutter::add`: the batch drops a record whose range starts past its end) after the
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
            // BLASTX's records are FASTA records (no PSSM query).
            bioseq_ids: None,
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

// The reads that the tests of the BLASTX engine make of their inputs (`input_source`).
#[cfg(test)]
pub(crate) use super::input_source::{parse_fasta, read_fasta};

#[cfg(test)]
#[path = "input_tests.rs"]
mod tests;
