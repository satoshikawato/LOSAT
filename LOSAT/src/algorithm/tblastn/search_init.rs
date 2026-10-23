//! Internal NCBI TBLASTN two-hit WordFinder of the (frame, chunk) units of the preliminary search
//! (`XTwoHit`).

#[cfg(test)]
use super::search_seed::encode_tblastn_lookup_query;
use crate::algorithm::blastp::extension::extend_two_hit;
use crate::config::ScoringMatrix;

// NCBI reference: c++/include/algo/blast/core/blast_extend.h:142-163
// typedef struct BlastUngappedData { Int4 q_start, s_start, length, score; } ...;
// typedef struct BlastInitHSP { BlastOffsetPair offsets;
//                              BlastUngappedData* ungapped_data; } ...;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct InitHsp {
    pub frame: i8,
    // NCBI c++/src/algo/blast/core/blast_engine.c:241-250,572-586:
    // backup.offset is retained until HSP offsets are adjusted before merge.
    pub chunk_offset: u32,
    pub q_seed: u32,
    pub s_seed: u32,
    pub q_start: i32,
    pub s_start: i32,
    pub length: i32,
    pub score: i32,
}

// NCBI reference: c++/src/algo/blast/core/blast_extend.c:46-65,162-175
// while (diag_array_length < (qlen+window_size)) diag_array_length <<= 1;
// diag_table->offset = window_size;
// if (ewp->diag_table->offset >= INT4_MAX / 4) { ... clear ... }
// else ewp->diag_table->offset += subject_length + ewp->diag_table->window;
pub(super) struct Diagonals {
    entries: Vec<(i32, bool)>,
    mask: i32,
    offset: i32,
    window: i32,
}

impl Diagonals {
    fn new(query_length: usize, window: i32) -> Self {
        let size = (query_length + window as usize).next_power_of_two();
        Self {
            entries: vec![(0, false); size],
            mask: size as i32 - 1,
            offset: window,
            window,
        }
    }

    fn finish_frame(&mut self, subject_length: usize) {
        if self.offset >= i32::MAX / 4 {
            self.offset = self.window;
            self.entries.fill((0, false));
        } else {
            self.offset += subject_length as i32 + self.window;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-547
// ```c
//             diag_coord = (query_offset - subject_offset) & diag_mask;
// ...
//             if (diag_array[diag_coord].flag) {
//                 /* If we've already extended past this hit, skip it. */
//                 if ((Int4) (subject_offset + diag_offset) <
//                     diag_array[diag_coord].last_hit) {
//                     continue;
//                 }
// ...
//                 if (diff >= window) {
//                     /* We are beyond the window for this diagonal; start a
//                        new hit */
//                     diag_array[diag_coord].last_hit =
//                         subject_offset + diag_offset;
//                     continue;
//                 }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:162-173
// ```c
//     if (ewp->diag_table) {
// ...
//             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
// ```
// XTwoHit is the two-hit pass of NCBI's word finder with its own diagonal array per (frame, chunk)
// unit. Each unit starts from the state of a fresh table, which is what the offset update at the
// end of the previous WordFinder call gives (see the comment above), so units can run in any order
// or at the same time.
// The two-hit pass of one (frame, chunk) unit on its
// own diagonal array.  NCBI c++/src/algo/blast/core/blast_extend.c:162-175 ends a
// WordFinder call with `diag_table->offset += subject_length + window`, which
// makes every entry of the previous call older than a full window: the first hit
// of the next call on that diagonal only records `last_hit` (aa_ungapped.c:
// 516-547), exactly as on a cleared array.  A unit therefore depends on nothing
// that an earlier unit left behind.
pub(super) struct XTwoHit<'a> {
    diagonals: Diagonals,
    pub query_sequence: &'a [u8],
    pub context_offsets: &'a [i32],
    pub x_dropoffs: &'a [i32],
    pub cutoff_scores: &'a [i32],
    pub matrix: ScoringMatrix,
    pub window: i32,
    pub word_size: usize,
}

impl<'a> XTwoHit<'a> {
    #[cfg(test)]
    pub(super) fn new(
        query_sequence: &'a [u8],
        context_offsets: &'a [i32],
        x_dropoffs: &'a [i32],
        cutoff_scores: &'a [i32],
        matrix: ScoringMatrix,
        window: i32,
        word_size: usize,
    ) -> Self {
        Self::with_table(
            None,
            query_sequence,
            context_offsets,
            x_dropoffs,
            cutoff_scores,
            matrix,
            window,
            word_size,
        )
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1002-1003,1377-1379,1412-1413
    // ```c
    //     if ((status = BlastExtendWordNew(query->length, word_params,
    //                                     &aux_struct->ewp)) != 0)
    // ...
    //        s_BlastSetUpAuxStructures(seq_src, lookup_wrap, word_params,
    //           ext_options, hit_options, query, query_info, &aux_struct)) != 0)
    // ...
    //     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //            != BLAST_SEQSRC_EOF) {
    // ```
    // NCBI's search thread makes its diagonal table once, before its loop over the subjects, and
    // keeps it for every subject. A table of the same length and window is taken over as it is:
    // `exit` after each unit leaves every entry older than a full window (see `exit`).
    /// `new` on `table` (a table an earlier unit left) when it has the length and window this
    /// query needs, otherwise on a new table.
    pub(super) fn with_table(
        table: Option<Diagonals>,
        query_sequence: &'a [u8],
        context_offsets: &'a [i32],
        x_dropoffs: &'a [i32],
        cutoff_scores: &'a [i32],
        matrix: ScoringMatrix,
        window: i32,
        word_size: usize,
    ) -> Self {
        let query_length = query_sequence.len().saturating_sub(1);
        let size = (query_length + window as usize).next_power_of_two();
        let diagonals = match table {
            Some(table) if table.entries.len() == size && table.window == window => table,
            _ => Diagonals::new(query_length, window),
        };
        Self {
            diagonals,
            query_sequence,
            context_offsets,
            x_dropoffs,
            cutoff_scores,
            matrix,
            window,
            word_size,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:608-609
    // ```c
    //     /* increment the offset in the diagonal array */
    //     Blast_ExtendWordExit(ewp, subject->length);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:166-171
    // ```c
    //         if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //             ewp->diag_table->offset = ewp->diag_table->window;
    //             s_BlastDiagClear(ewp->diag_table);
    //         } else {
    //             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    //         }
    // ```
    // The end of NCBI's WordFinder call for one (frame, chunk). Every entry the unit wrote holds
    // at most `offset + subject_length - 1`, so after the step each one lies more than a window
    // behind any hit of the next unit: its first hit on a diagonal only records `last_hit` and
    // clears `flag` (the tests of `pairs`), as on a table made afresh by `Diagonals::new`. The
    // next unit sees this instead of a new table.
    /// End a unit of `subject_length` residues as NCBI does (the table is kept).
    pub(super) fn exit(&mut self, subject_length: usize) {
        self.diagonals.finish_frame(subject_length);
    }

    /// No NCBI counterpart: hands the diagonal table back for `with_table`; it does not change
    /// any value NCBI computes.
    pub(super) fn into_table(self) -> Diagonals {
        self.diagonals
    }

    /// The diagonal table, for tests.
    #[cfg(test)]
    pub(super) fn x_table(&self) -> &[(i32, bool)] {
        &self.diagonals.entries
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-606
    // ```c
    //             diag_coord = (query_offset - subject_offset) & diag_mask;
    // ...
    //                 score = s_BlastAaExtendTwoHit(matrix, subject, query,
    // ...
    //                 if (score >= cutoffs->cutoff_score)
    //                     BlastSaveInitHsp(ungapped_hsps, hsp_q, hsp_s,
    // ...
    //                 if (right_extend) {
    //                     diag_array[diag_coord].flag = 1;
    //                     diag_array[diag_coord].last_hit =
    //                         s_last_off - (wordsize - 1) + diag_offset;
    // ```
    // The body is the per-hit loop of the two-hit word finder (the same tests, extension and HSP
    // saving of aa_ungapped.c:516-606), run on one filled offset array.
    /// The loop body of the two-hit word finder for one filled offset array.
    #[inline]
    pub(super) fn pairs(
        &mut self,
        frame: i8,
        chunk_offset: u32,
        subject_sequence: &[u8],
        pairs: &[crate::algorithm::tblastx::blast_aascan::BlastOffsetPair],
        hits: &mut Vec<InitHsp>,
    ) {
        let window = self.window;
        let word_size = self.word_size;
        let word_size_i32 = word_size as i32;
        let diagonals = &mut self.diagonals;
        for pair in pairs {
            let q = pair.q_off as i32;
            let s = pair.s_off as i32;
            // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:516-547
            let index = ((q - s) & diagonals.mask) as usize;
            let (last_hit, flag) = &mut diagonals.entries[index];
            if *flag {
                if s + diagonals.offset < *last_hit {
                    continue;
                }
                *last_hit = s + diagonals.offset;
                *flag = false;
                continue;
            }
            let previous = *last_hit - diagonals.offset;
            let diff = s - previous;
            if diff >= window {
                *last_hit = s + diagonals.offset;
                continue;
            }
            if diff < word_size_i32 {
                continue;
            }
            // NCBI c++/src/algo/blast/core/aa_ungapped.c:562-579:
            // curr_context = BSearchContextInfo(query_offset, query_info);
            // if (query_offset - diff <
            //     query_info->contexts[curr_context].query_offset) continue;
            let context = self.context_offsets.partition_point(|&offset| offset <= q) - 1;
            if q - diff < self.context_offsets[context] {
                *last_hit = s + diagonals.offset;
                continue;
            }
            // NCBI c++/src/algo/blast/core/aa_ungapped.c:570-583,1089-1155:
            // s_BlastAaExtendTwoHit(matrix, ..., wordsize, ...);
            let result = extend_two_hit(
                self.matrix,
                self.query_sequence,
                subject_sequence,
                (previous + word_size_i32) as usize,
                s as usize,
                q as usize,
                self.x_dropoffs[context],
                word_size,
            );
            let Some(result) = result else {
                continue;
            };
            // NCBI c++/src/algo/blast/core/aa_ungapped.c:570-590:
            // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
            if result.ungapped_data.score >= self.cutoff_scores[context] {
                hits.push(InitHsp {
                    frame,
                    chunk_offset,
                    q_seed: pair.q_off,
                    s_seed: pair.s_off,
                    q_start: result.ungapped_data.q_start,
                    s_start: result.ungapped_data.s_start,
                    length: result.ungapped_data.length,
                    score: result.ungapped_data.score,
                });
            }
            // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:588-607
            if result.right_extend {
                *flag = true;
                *last_hit = result.s_last_off - (word_size_i32 - 1) + diagonals.offset;
            } else {
                *last_hit = s + diagonals.offset;
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-552 (the tests of
    // `pairs`; tblastn/x_tbn_hitbf.rs module doc)
    // EXPERIMENT (LOSAT_X_TBNHITBF): `x_hits_bf` takes the too-far and overlap outcomes; the
    // body of `pairs` runs on each hit it stops at (a flag-1 cell or a candidate), in order.
    /// `pairs` with the too-far and overlap outcomes taken by `x_tbn_hitbf::x_hits_bf`; `pairs`
    /// itself when `x_tbn_hitbf::fits` does not hold.
    #[inline(never)]
    pub(super) fn x_pairs_bf(
        &mut self,
        frame: i8,
        chunk_offset: u32,
        subject_sequence: &[u8],
        pairs: &[crate::algorithm::tblastx::blast_aascan::BlastOffsetPair],
        hits: &mut Vec<InitHsp>,
    ) {
        if !super::x_tbn_hitbf::fits(self.window, self.word_size) {
            self.pairs(frame, chunk_offset, subject_sequence, pairs, hits);
            return;
        }
        let word_size = self.word_size as i32;
        let mut done = 0;
        while done < pairs.len() {
            done += super::x_tbn_hitbf::x_hits_bf(
                &mut self.diagonals.entries,
                self.diagonals.offset,
                self.window,
                word_size,
                &pairs[done..],
            );
            if let Some(pair) = pairs.get(done) {
                self.pairs(
                    frame,
                    chunk_offset,
                    subject_sequence,
                    std::slice::from_ref(pair),
                    hits,
                );
                done += 1;
            }
        }
    }

    // No NCBI counterpart: runs `pairs` and `x_pairs_bf` on one offset array and compares them
    // (LOSAT_X_TBNHITBFSHADOW); the table and init HSPs kept are those of `x_pairs_bf`.
    /// EXPERIMENT (LOSAT_X_TBNHITBFSHADOW): `pairs` on a copy of the table, then `x_pairs_bf`;
    /// panics at the first difference of the init HSPs or of the tables.
    #[inline(never)]
    pub(super) fn x_pairs_bf_shadow(
        &mut self,
        frame: i8,
        chunk_offset: u32,
        subject_sequence: &[u8],
        pairs: &[crate::algorithm::tblastx::blast_aascan::BlastOffsetPair],
        hits: &mut Vec<InitHsp>,
    ) {
        let outside = !super::x_tbn_hitbf::fits(self.window, self.word_size);
        super::x_tbn_hitbf::x_shadow_count(pairs.len(), outside);
        if outside {
            self.pairs(frame, chunk_offset, subject_sequence, pairs, hits);
            return;
        }
        if pairs.is_empty() {
            return;
        }
        let before = self.diagonals.entries.clone();
        let mut reference_hits = Vec::new();
        self.pairs(
            frame,
            chunk_offset,
            subject_sequence,
            pairs,
            &mut reference_hits,
        );
        let reference_table = std::mem::replace(&mut self.diagonals.entries, before);
        let start = hits.len();
        self.x_pairs_bf(frame, chunk_offset, subject_sequence, pairs, hits);
        super::x_tbn_hitbf::x_shadow_compare(
            frame,
            chunk_offset,
            &hits[start..],
            &reference_hits,
            &self.diagonals.entries,
            &reference_table,
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:274-309
    // ```c
    //     if (0 == (result = BLAST_CMP(h2->ungapped_data->score,
    //                                  h1->ungapped_data->score)) &&
    //         0 == (result = BLAST_CMP(h1->ungapped_data->s_start,
    //                                  h2->ungapped_data->s_start)) &&
    //         0 == (result = BLAST_CMP(h2->ungapped_data->length,
    //                                  h1->ungapped_data->length)) &&
    //         0 == (result = BLAST_CMP(h1->ungapped_data->q_start,
    //                                  h2->ungapped_data->q_start))) {
    // ...
    //     qsort(init_hitlist->init_hsp_array, init_hitlist->total,
    //           sizeof(BlastInitHSP), score_compare_match);
    // ```
    // Same order as Blast_InitHitListSortByScore: score descending, then s_start, length
    // descending, q_start. The sort is stable, as glibc's qsort is under the pinned NCBI build.
    /// NCBI c++/src/algo/blast/core/blast_extend.c:273-313 (stable, as glibc's
    /// qsort under the pinned NCBI BLAST+).
    pub(super) fn sort(hits: &mut [InitHsp]) {
        hits.sort_by(|a, b| {
            b.score
                .cmp(&a.score)
                .then(a.s_start.cmp(&b.s_start))
                .then(b.length.cmp(&a.length))
                .then(a.q_start.cmp(&b.q_start))
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::tblastn::search_gapped::{preliminary_stages, word_finder_profile};
    use crate::sequence::packed_subject::SubjectNt;
    use std::fs;

    // NCBI c++/src/algo/blast/core/blast_engine.c:259-310,478-487:
    // subject->seq_ranges and subject->length are filled before WordFinder.
    // The comparison-only probe records those exact input rows at the call.
    fn ncbi_chunk_range_rows(path: &str) -> Vec<(i8, usize, Vec<(i32, i32)>)> {
        let trace = fs::read_to_string(path).unwrap();
        let mut calls: Vec<(i8, usize, Vec<(i32, i32)>)> = Vec::new();
        for line in trace.lines() {
            let fields: Vec<_> = line.split('\t').collect();
            match fields[0] {
                "RANGE_CALL" => {
                    assert_eq!(fields[1].parse::<usize>().unwrap(), calls.len());
                    calls.push((
                        fields[2].parse().unwrap(),
                        fields[3].parse().unwrap(),
                        Vec::new(),
                    ));
                }
                "RANGE" => {
                    let call = fields[1].parse::<usize>().unwrap();
                    assert_eq!(call, calls.len() - 1);
                    assert_eq!(fields[2].parse::<usize>().unwrap(), calls[call].2.len());
                    calls[call]
                        .2
                        .push((fields[3].parse().unwrap(), fields[4].parse().unwrap()));
                }
                "SET_RANGES" | "SET_RANGE" => {}
                _ => panic!("unexpected NCBI range trace row"),
            }
        }
        calls
    }

    fn read_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
        let text = fs::read_to_string(path).unwrap();
        let mut records: Vec<(String, Vec<u8>)> = Vec::new();
        for line in text.lines() {
            if let Some(id) = line.strip_prefix('>') {
                records.push((id.to_string(), Vec::new()));
            } else {
                records.last_mut().unwrap().1.extend(line.as_bytes());
            }
        }
        records
    }

    // NCBI c++/src/algo/blast/core/aa_ungapped.c:516-614:
    // the selected matrix, wordsize, x-dropoff, and per-context cutoff
    // determine two-hit extension before frame-local score sorting.
    #[test]
    fn alternate_matrix_word2_wordfinder_hsps_match_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/alternate_matrix_word2_20260924/"
        );
        let queries = read_fasta(&format!("{root}query.faa"));
        let refs: Vec<&[u8]> = queries.iter().map(|(_, seq)| seq.as_slice()).collect();
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let trace = fs::read_to_string(format!("{root}candidate.stderr")).unwrap();
        let observed_cutoffs: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("PARAM\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[2].parse::<i32>().unwrap(),
                    f[3].parse::<i32>().unwrap(),
                    f[4].parse::<i32>().unwrap(),
                )
            })
            .collect();
        assert_eq!(observed_cutoffs, [(21, 21, 0); 6]);
        let trace = fs::read_to_string(format!("{root}wordfinder.stderr")).unwrap();
        let frame_order = [1, 2, 3, -1, -2, -3];
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("INIT\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                InitHsp {
                    frame: frame_order[f[1].parse::<usize>().unwrap()],
                    chunk_offset: 0,
                    q_seed: f[3].parse().unwrap(),
                    s_seed: f[4].parse().unwrap(),
                    q_start: f[5].parse().unwrap(),
                    s_start: f[6].parse().unwrap(),
                    length: f[7].parse().unwrap(),
                    score: f[8].parse().unwrap(),
                }
            })
            .collect();
        let actual = preliminary_stages(
            &refs,
            subject,
            1,
            word_finder_profile(
                None,
                16,
                60,
                &[21, 21, 0],
                &[0, 0, i32::MAX],
                false,
                ScoringMatrix::Blosum45,
                2,
            ),
        )
        .unwrap()
        .initial;
        if actual != expected {
            let first = actual
                .iter()
                .zip(&expected)
                .position(|(left, right)| left != right)
                .unwrap_or(actual.len().min(expected.len()));
            panic!("first BLOSUM45/word-2 WordFinder difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
            actual.get(first), expected.get(first), actual.len(), expected.len());
        }
    }

    // NCBI c++/src/algo/blast/core/aa_ungapped.c:562-583:
    // BSearchContextInfo finds the query context before two-hit extension;
    // the left word must remain within that context's query_offset.
    // NCBI c++/src/algo/blast/core/aa_ungapped.c:200-234:
    // Blast_InitHitListSortByScore(init_hitlist) fixes frame-local HSP order.
    #[test]
    fn multi_query_wordfinder_hsps_match_ncbi_in_order() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/multi_query_20260924/"
        );
        let queries = read_fasta(&format!("{root}query.faa"));
        let refs: Vec<&[u8]> = queries
            .iter()
            .map(|(_, sequence)| sequence.as_slice())
            .collect();
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let trace = fs::read_to_string(format!("{root}wordfinder.stderr")).unwrap();
        let frame_order = [1, 2, 3, -1, -2, -3];
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("INIT\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                InitHsp {
                    frame: frame_order[f[1].parse::<usize>().unwrap()],
                    chunk_offset: 0,
                    q_seed: f[3].parse().unwrap(),
                    s_seed: f[4].parse().unwrap(),
                    q_start: f[5].parse().unwrap(),
                    s_start: f[6].parse().unwrap(),
                    length: f[7].parse().unwrap(),
                    score: f[8].parse().unwrap(),
                }
            })
            .collect();
        let actual = preliminary_stages(
            &refs,
            subject,
            1,
            word_finder_profile(
                None,
                13,
                40,
                &[16, 16, 0],
                &[0, 0, i32::MAX],
                false,
                ScoringMatrix::Blosum62,
                3,
            ),
        )
        .unwrap()
        .initial;
        if actual != expected {
            let first = actual
                .iter()
                .zip(&expected)
                .position(|(left, right)| left != right)
                .unwrap_or(actual.len().min(expected.len()));
            panic!(
                "first multi-query WordFinder difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
                actual.get(first),
                expected.get(first),
                actual.len(),
                expected.len()
            );
        }
    }

    // NCBI c++/src/algo/blast/core/aa_ungapped.c:570-590:
    // cutoffs = word_params->cutoffs + curr_context;
    // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
    // The saved natural multi-query fixture has score 20 in context 0 and
    // score 393 in context 1, so equality and just-above cutoffs are distinct.
    #[test]
    fn wordfinder_cutoff_is_selected_per_query_context_at_score_boundary() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/multi_query_20260924/"
        );
        let queries = read_fasta(&format!("{root}query.faa"));
        let refs: Vec<&[u8]> = queries.iter().map(|(_, seq)| seq.as_slice()).collect();
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let equal = preliminary_stages(
            &refs,
            subject,
            1,
            word_finder_profile(
                None,
                13,
                40,
                &[16, 16, 0],
                &[20, 393, i32::MAX],
                false,
                ScoringMatrix::Blosum62,
                3,
            ),
        )
        .unwrap()
        .initial;
        assert!(equal.iter().any(|h| h.score == 20 && h.q_seed < 121));
        assert!(equal
            .iter()
            .any(|h| h.score == 393 && (121..192).contains(&h.q_seed)));
        let above = preliminary_stages(
            &refs,
            subject,
            1,
            word_finder_profile(
                None,
                13,
                40,
                &[16, 16, 0],
                &[21, 394, i32::MAX],
                false,
                ScoringMatrix::Blosum62,
                3,
            ),
        )
        .unwrap()
        .initial;
        assert!(!above.iter().any(|h| h.score == 20 && h.q_seed < 121));
        assert!(!above
            .iter()
            .any(|h| h.score == 393 && (121..192).contains(&h.q_seed)));
    }

    // NCBI c++/src/algo/blast/core/aa_ungapped.c:562-590:
    // cutoffs = word_params->cutoffs + curr_context;
    // s_BlastAaExtendTwoHit(..., cutoffs->x_dropoff, ...);
    // The comparison-only real call changes context 1 x-drop from 16 to 1.
    #[test]
    fn real_wordfinder_distinct_context_xdrop_matches_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/"
        );
        let queries = read_fasta(&format!("{root}multi_query_20260924/query.faa"));
        let refs: Vec<&[u8]> = queries.iter().map(|(_, seq)| seq.as_slice()).collect();
        let subject = &read_fasta(&format!("{root}multi_query_20260924/subjects.fna"))[0].1;
        let trace =
            fs::read_to_string(format!("{root}word_xdrop_real_path_20260924/trace.tsv")).unwrap();
        let injections: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("WORD_XDROP_INJECT\t"))
            .collect();
        assert_eq!(injections.len(), 6);
        for (call, line) in injections.iter().enumerate() {
            assert_eq!(*line, format!("WORD_XDROP_INJECT\t{call}\t16\t1"));
        }
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("INIT\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                let call: usize = f[1].parse().unwrap();
                InitHsp {
                    frame: [1, 2, 3, -1, -2, -3][call],
                    chunk_offset: 0,
                    q_seed: f[3].parse().unwrap(),
                    s_seed: f[4].parse().unwrap(),
                    q_start: f[5].parse().unwrap(),
                    s_start: f[6].parse().unwrap(),
                    length: f[7].parse().unwrap(),
                    score: f[8].parse().unwrap(),
                }
            })
            .collect();
        assert_eq!(expected.len(), 32);
        let actual = preliminary_stages(
            &refs,
            subject,
            1,
            word_finder_profile(
                None,
                13,
                40,
                &[16, 1, 0],
                &[0, 0, i32::MAX],
                false,
                ScoringMatrix::Blosum62,
                3,
            ),
        )
        .unwrap()
        .initial;
        if actual != expected {
            let first = actual
                .iter()
                .zip(&expected)
                .position(|(left, right)| left != right)
                .unwrap_or(actual.len().min(expected.len()));
            panic!("first distinct-context x-drop HSP difference at {first}: Rust={:?}, NCBI={:?}; counts Rust={} NCBI={}",
                   actual.get(first), expected.get(first), actual.len(), expected.len());
        }
    }

    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:816-831
    // subject->seq_ranges[i].left = CONV_NUCL2PROT_COORDINATES(...);
    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:575-614
    // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
    #[test]
    fn lowercase_mask_init_hsps_match_ncbi_scores_and_internal_offsets() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/lowercase_20260923/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subjects = read_fasta(&format!("{root}subjects.fna"));
        let trace = fs::read_to_string(format!("{root}ncbi_init_trace.tsv")).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .skip(1)
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[1].to_string(),
                    InitHsp {
                        frame: f[2].parse().unwrap(),
                        chunk_offset: 0,
                        q_seed: f[4].parse().unwrap(),
                        s_seed: f[5].parse().unwrap(),
                        q_start: f[6].parse().unwrap(),
                        s_start: f[7].parse().unwrap(),
                        length: f[8].parse().unwrap(),
                        score: f[9].parse().unwrap(),
                    },
                )
            })
            .collect();
        let mut actual = Vec::new();
        for (id, subject) in &subjects {
            for hit in preliminary_stages(
                &[query],
                subject,
                1,
                word_finder_profile(None, 13, 40, &[16], &[0], true, ScoringMatrix::Blosum62, 3),
            )
            .unwrap()
            .initial
            {
                actual.push((id.clone(), hit));
            }
        }
        assert_eq!(actual, expected);
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:575-614
    // score = s_BlastAaExtendTwoHit(...);
    // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
    #[test]
    fn ambiguity_spectrum_init_hsp_matches_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/ambiguity_20260923/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{root}subjects.fna"))[0].1;
        let trace = fs::read_to_string(format!("{root}ncbi_init_trace.tsv")).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .skip(1)
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                InitHsp {
                    frame: f[2].parse().unwrap(),
                    chunk_offset: 0,
                    q_seed: f[4].parse().unwrap(),
                    s_seed: f[5].parse().unwrap(),
                    q_start: f[6].parse().unwrap(),
                    s_start: f[7].parse().unwrap(),
                    length: f[8].parse().unwrap(),
                    score: f[9].parse().unwrap(),
                }
            })
            .collect();
        let actual = preliminary_stages(
            &[query],
            subject,
            1,
            word_finder_profile(None, 13, 40, &[16], &[0], false, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .initial;
        assert_eq!(actual, expected);
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:200-234,575-614
    // status = s_BlastAaWordFinder_TwoHit(..., init_hitlist, ...);
    // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
    // Blast_InitHitListSortByScore(init_hitlist);
    // NCBI reference: c++/src/algo/blast/core/blast_engine.c:747-775
    // BLAST_GetAllTranslations(..., subject->gen_code_string, ...);
    #[test]
    fn code32_api_oracle_init_hsps_match_scores_internal_offsets_and_order() {
        let fixture = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_a/fixtures/"
        );
        let trace_root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/code32_20260923/"
        );
        let query = &read_fasta(&format!("{fixture}query.faa"))[0].1;
        let subject = &read_fasta(&format!("{fixture}subject_code32.fna"))[0].1;
        let trace = fs::read_to_string(format!("{trace_root}ncbi_init_trace.tsv")).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .skip(1)
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                InitHsp {
                    frame: f[0].parse().unwrap(),
                    chunk_offset: 0,
                    q_seed: f[2].parse().unwrap(),
                    s_seed: f[3].parse().unwrap(),
                    q_start: f[4].parse().unwrap(),
                    s_start: f[5].parse().unwrap(),
                    length: f[6].parse().unwrap(),
                    score: f[7].parse().unwrap(),
                }
            })
            .collect();
        let actual = preliminary_stages(
            &[query],
            subject,
            32,
            word_finder_profile(
                None,
                13,
                40,
                &[16],
                &[13],
                false,
                ScoringMatrix::Blosum62,
                3,
            ),
        )
        .unwrap()
        .initial;
        assert_eq!(actual, expected);
    }

    // NCBI reference: c++/src/algo/blast/core/aa_ungapped.c:200-234,575-614
    // score = s_BlastAaExtendTwoHit(..., cutoffs->x_dropoff, ...);
    // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
    // Blast_InitHitListSortByScore(init_hitlist);
    #[test]
    fn saved_init_hsps_match_ncbi_fixture_scores_and_internal_offsets() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/run_20260923/"
        );
        let query = &read_fasta(&format!("{root}query.faa"))[0].1;
        let subjects = read_fasta(&format!("{root}subjects.fna"));
        let trace = fs::read_to_string(format!("{root}ncbi_init_trace.tsv")).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .skip(1)
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    f[1].to_string(),
                    InitHsp {
                        frame: f[2].parse().unwrap(),
                        chunk_offset: 0,
                        q_seed: f[4].parse().unwrap(),
                        s_seed: f[5].parse().unwrap(),
                        q_start: f[6].parse().unwrap(),
                        s_start: f[7].parse().unwrap(),
                        length: f[8].parse().unwrap(),
                        score: f[9].parse().unwrap(),
                    },
                )
            })
            .collect();
        let params = fs::read_to_string(format!("{root}ncbi_word_params.tsv")).unwrap();
        let observed: Vec<_> = params
            .lines()
            .skip(1)
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (f[2].parse::<i32>().unwrap(), f[3].parse::<i32>().unwrap())
            })
            .collect();
        assert!(observed.iter().all(|pair| *pair == (16, 0)));
        let mut actual = Vec::new();
        for (id, subject) in &subjects {
            for hit in preliminary_stages(
                &[query],
                subject,
                1,
                word_finder_profile(None, 13, 40, &[16], &[0], false, ScoringMatrix::Blosum62, 3),
            )
            .unwrap()
            .initial
            {
                actual.push((id.clone(), hit));
            }
        }
        assert_eq!(actual, expected);
    }
    // NCBI c++/src/algo/blast/core/blast_engine.c:478-552:
    // s_GetNextSubjectChunk precedes each WordFinder call and the hit list
    // is reset before the next chunk. Subject offsets stay chunk-local.
    #[test]
    fn long_subject_chunk_wordfinder_hsps_match_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/"
        );
        let query = &read_fasta(&format!("{root}long_chunk_20260924/query.faa"))[0].1;
        let inserts = read_fasta(&format!("{root}run_20260923/subjects.fna"));
        let plus1 = &inserts.iter().find(|(id, _)| id == "plus1").unwrap().1;
        assert_eq!(plus1.len(), 362);
        let mut subject = Vec::with_capacity(15_000_962);
        for _ in 0..4_999_950 {
            subject.extend_from_slice(b"ATG");
        }
        subject.extend_from_slice(plus1);
        for _ in 0..250 {
            subject.extend_from_slice(b"ATG");
        }
        assert_eq!(subject.len(), 15_000_962);
        let trace =
            fs::read_to_string(format!("{root}long_chunk_20260924/frame_chunks.tsv")).unwrap();
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("INIT\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                let call: usize = f[1].parse().unwrap();
                InitHsp {
                    frame: [1, 2, 3, -1, -2, -3][call / 2],
                    chunk_offset: if call % 2 == 0 { 0 } else { 4_999_900 },
                    q_seed: f[3].parse().unwrap(),
                    s_seed: f[4].parse().unwrap(),
                    q_start: f[5].parse().unwrap(),
                    s_start: f[6].parse().unwrap(),
                    length: f[7].parse().unwrap(),
                    score: f[8].parse().unwrap(),
                }
            })
            .collect();
        // NCBI c++/src/algo/blast/core/aa_ungapped.c:570-590:
        // cutoffs = word_params->cutoffs + curr_context;
        // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
        // Retained long-fixture PARAM rows give cutoff 28.
        let actual = preliminary_stages(
            &[query],
            &subject,
            1,
            word_finder_profile(
                None,
                13,
                40,
                &[16],
                &[28],
                false,
                ScoringMatrix::Blosum62,
                3,
            ),
        )
        .unwrap()
        .initial;
        assert_eq!(actual, expected);
    }
    // NCBI c++/src/algo/blast/core/blast_engine.c:283-310:
    // a soft range ending at the second chunk's start is included, then
    // clipped to zero width before aa_ungapped.c:495-500 scans range zero.
    #[test]
    fn masked_chunk_boundary_wordfinder_hsp_matches_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/"
        );
        let query = &read_fasta(&format!("{root}masked_chunk_boundary_20260924/query.faa"))[0].1;
        let inserts = read_fasta(&format!("{root}run_20260923/subjects.fna"));
        let plus1 = &inserts.iter().find(|(id, _)| id == "plus1").unwrap().1;
        let mut subject = Vec::with_capacity(15_000_962);
        for _ in 0..4_999_950 {
            subject.extend_from_slice(b"ATG");
        }
        subject.extend_from_slice(plus1);
        for _ in 0..250 {
            subject.extend_from_slice(b"ATG");
        }
        subject[14_999_700..15_000_000].make_ascii_lowercase();
        // The (frame, chunk, scan ranges) units the search runs: the unit path's plan.
        let plan = super::super::search_seed::x_scan_plan(
            &[query],
            SubjectNt::Letters(&subject).into(),
            1,
            None,
            13,
            true,
            ScoringMatrix::Blosum62,
            3,
            &std::sync::OnceLock::new(),
            None,
        )
        .unwrap()
        .unwrap();
        let actual_ranges: Vec<_> = plan
            .subject
            .units
            .iter()
            .map(|(frame_index, chunk, ranges)| {
                (
                    plan.subject.frames[*frame_index].frame,
                    chunk.length,
                    ranges.clone(),
                )
            })
            .collect();
        assert_eq!(
            actual_ranges,
            ncbi_chunk_range_rows(&format!(
                "{root}masked_chunk_boundary_20260924/chunk_ranges.tsv"
            ))
        );
        drop(plan);
        let trace = fs::read_to_string(format!(
            "{root}masked_chunk_boundary_20260924/frame_chunks.tsv"
        ))
        .unwrap();
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("INIT\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                let call: usize = f[1].parse().unwrap();
                InitHsp {
                    frame: [1, 2, 3, -1, -2, -3][call / 2],
                    chunk_offset: if call % 2 == 0 { 0 } else { 4_999_900 },
                    q_seed: f[3].parse().unwrap(),
                    s_seed: f[4].parse().unwrap(),
                    q_start: f[5].parse().unwrap(),
                    s_start: f[6].parse().unwrap(),
                    length: f[7].parse().unwrap(),
                    score: f[8].parse().unwrap(),
                }
            })
            .collect();
        // NCBI c++/src/algo/blast/core/aa_ungapped.c:570-590:
        // cutoffs = word_params->cutoffs + curr_context;
        // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
        // Retained long-fixture PARAM rows give cutoff 28.
        let actual = preliminary_stages(
            &[query],
            &subject,
            1,
            word_finder_profile(None, 13, 40, &[16], &[28], true, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .initial;
        assert_eq!(actual, expected);
    }
    // NCBI c++/src/algo/blast/core/blast_engine.c:283-310,478-500:
    // SUBJECT_SPLIT_NO_RANGE skips the middle positive-frame chunk; the
    // later chunk still runs WordFinder with its chunk-local offsets.
    #[test]
    fn masked_no_range_middle_chunk_skips_and_later_hsp_matches_ncbi() {
        let root = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/"
        );
        let query = &read_fasta(&format!("{root}no_range_middle_20260924/query.faa"))[0].1;
        let inserts = read_fasta(&format!("{root}run_20260923/subjects.fna"));
        let plus1 = &inserts.iter().find(|(id, _)| id == "plus1").unwrap().1;
        let mut subject = Vec::with_capacity(30_001_262);
        for _ in 0..10_000_050 {
            subject.extend_from_slice(b"ATG");
        }
        subject.extend_from_slice(plus1);
        for _ in 0..250 {
            subject.extend_from_slice(b"ATG");
        }
        assert_eq!(subject.len(), 30_001_262);
        subject[14_997_000..30_000_000].make_ascii_lowercase();
        let trace =
            fs::read_to_string(format!("{root}no_range_middle_20260924/frame_chunks.tsv")).unwrap();
        let observed_calls: Vec<(i8, usize)> = trace
            .lines()
            .filter(|line| line.starts_with("FRAME_CHUNK\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (f[2].parse().unwrap(), f[3].parse().unwrap())
            })
            .collect();
        // The (frame, chunk, scan ranges) units the search runs: the unit path's plan.
        let plan = super::super::search_seed::x_scan_plan(
            &[query],
            SubjectNt::Letters(&subject).into(),
            1,
            None,
            13,
            true,
            ScoringMatrix::Blosum62,
            3,
            &std::sync::OnceLock::new(),
            None,
        )
        .unwrap()
        .unwrap();
        let actual_ranges: Vec<_> = plan
            .subject
            .units
            .iter()
            .map(|(frame_index, chunk, ranges)| {
                (
                    plan.subject.frames[*frame_index].frame,
                    chunk.length,
                    ranges.clone(),
                )
            })
            .collect();
        assert_eq!(
            actual_ranges
                .iter()
                .map(|(frame, length, _)| (*frame, *length))
                .collect::<Vec<_>>(),
            observed_calls
        );
        assert_eq!(
            actual_ranges,
            ncbi_chunk_range_rows(&format!("{root}no_range_middle_20260924/chunk_ranges.tsv"))
        );
        drop(plan);
        let expected: Vec<_> = trace
            .lines()
            .filter(|line| line.starts_with("INIT\t"))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                assert_eq!(f[1], "1");
                InitHsp {
                    frame: 1,
                    chunk_offset: 9_999_800,
                    q_seed: f[3].parse().unwrap(),
                    s_seed: f[4].parse().unwrap(),
                    q_start: f[5].parse().unwrap(),
                    s_start: f[6].parse().unwrap(),
                    length: f[7].parse().unwrap(),
                    score: f[8].parse().unwrap(),
                }
            })
            .collect();
        // NCBI c++/src/algo/blast/core/aa_ungapped.c:570-590:
        // cutoffs = word_params->cutoffs + curr_context;
        // if (score >= cutoffs->cutoff_score) BlastSaveInitHsp(...);
        // Retained long-fixture PARAM rows give cutoff 30.
        let actual = preliminary_stages(
            &[query],
            &subject,
            1,
            word_finder_profile(None, 13, 40, &[16], &[30], true, ScoringMatrix::Blosum62, 3),
        )
        .unwrap()
        .initial;
        assert_eq!(actual, expected);
    }

    /// The concatenated queries and their context offsets, as the unit path builds them.
    fn x_concatenated(queries: &[&[u8]]) -> (Vec<u8>, Vec<i32>) {
        let mut sequence = Vec::new();
        let mut offsets = Vec::new();
        for query in queries {
            offsets.push(sequence.len() as i32);
            sequence
                .extend_from_slice(&encode_tblastn_lookup_query(query, None, false).aa_seq[1..]);
        }
        (sequence, offsets)
    }

    /// Every exact word of three residues the subject shares with the queries, in subject order
    /// (the order of the subject scan).
    fn x_word_pairs(
        query: &[u8],
        subject: &[u8],
    ) -> Vec<crate::algorithm::tblastx::blast_aascan::BlastOffsetPair> {
        let mut pairs = Vec::new();
        for s in 0..subject.len().saturating_sub(2) {
            for q in 0..query.len().saturating_sub(2) {
                if query[q..q + 3] == subject[s..s + 3] && query[q..q + 3].iter().all(|&r| r != 0) {
                    pairs.push(crate::algorithm::tblastx::blast_aascan::BlastOffsetPair {
                        q_off: q as u32,
                        s_off: s as u32,
                    });
                }
            }
        }
        pairs
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:166-171
    // ```c
    //         if (ewp->diag_table->offset >= INT4_MAX / 4) {
    //             ewp->diag_table->offset = ewp->diag_table->window;
    //             s_BlastDiagClear(ewp->diag_table);
    //         } else {
    //             ewp->diag_table->offset += subject_length + ewp->diag_table->window;
    //         }
    // ```
    // The unit path's kept state: one table carried through many units with `exit` after each, as
    // NCBI's thread carries it through its subjects, gives every unit the hits of a cleared table
    // (a table made afresh): from a new table, across the clear near INT4_MAX / 4, and when
    // `with_table` takes it over for another query set of the same table length. The units are
    // fragments of the queries with substitutions, so that two-hit extensions (and the `flag` path
    // after a right extension) happen on many diagonals.
    #[test]
    fn x_two_hit_table_kept_across_units_gives_the_hits_of_a_cleared_table() {
        let mut seed = 0x2545_f491_u64;
        let mut next = move |bound: usize| {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((seed >> 33) as usize) % bound
        };
        const RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
        fn random(length: usize, next: &mut dyn FnMut(usize) -> usize) -> Vec<u8> {
            (0..length).map(|_| RESIDUES[next(20)]).collect()
        }
        let queries: Vec<Vec<u8>> = [150, 90, 210]
            .into_iter()
            .map(|length| random(length, &mut next))
            .collect();
        let mut units = Vec::new();
        for _ in 0..60 {
            let mut unit = random(next(40), &mut next);
            for _ in 0..1 + next(3) {
                let query = &queries[next(queries.len())];
                let start = next(query.len() - 30);
                let end = (start + 30 + next(80)).min(query.len());
                for &residue in &query[start..end] {
                    unit.push(if next(100) < 15 {
                        RESIDUES[next(20)]
                    } else {
                        residue
                    });
                }
                unit.extend(random(next(30), &mut next));
            }
            units.push(
                encode_tblastn_lookup_query(&unit, None, false).aa_seq[1..=unit.len()].to_vec(),
            );
        }
        let x_dropoffs = [16, 16, 16];
        let cutoffs = [24, 24, 24];
        let window = 40;
        /// The hits of every unit in turn, each on a table made afresh by `fresh` (a cleared table).
        fn cleared_hits<'q>(
            fresh: impl Fn() -> XTwoHit<'q>,
            query: &[u8],
            units: &[Vec<u8>],
        ) -> Vec<Vec<InitHsp>> {
            units
                .iter()
                .map(|unit| {
                    let mut hits = Vec::new();
                    fresh().pairs(1, 0, unit, &x_word_pairs(query, unit), &mut hits);
                    hits
                })
                .collect()
        }
        /// The hits of every unit in turn on one kept table (`exit` after each unit).
        fn kept_hits(
            two_hit: &mut XTwoHit<'_>,
            query: &[u8],
            units: &[Vec<u8>],
        ) -> Vec<Vec<InitHsp>> {
            units
                .iter()
                .map(|unit| {
                    let mut hits = Vec::new();
                    two_hit.pairs(1, 0, unit, &x_word_pairs(query, unit), &mut hits);
                    two_hit.exit(unit.len());
                    hits
                })
                .collect()
        }
        let set_a: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
        let (query_a, offsets_a) = x_concatenated(&set_a);
        let new_a = |table| {
            XTwoHit::with_table(
                table,
                &query_a,
                &offsets_a,
                &x_dropoffs,
                &cutoffs,
                ScoringMatrix::Blosum62,
                window,
                3,
            )
        };
        let cleared = cleared_hits(|| new_a(None), &query_a, &units);
        let found: usize = cleared.iter().map(Vec::len).sum();
        assert!(found >= 60, "the units must give hits ({found})");

        let mut kept = new_a(None);
        assert_eq!(
            kept_hits(&mut kept, &query_a, &units),
            cleared,
            "kept table"
        );

        let mut near_limit = new_a(Some(kept.into_table()));
        let start = i32::MAX / 4 - 5_000;
        near_limit.diagonals.offset = start;
        assert_eq!(
            kept_hits(&mut near_limit, &query_a, &units),
            cleared,
            "across the clear"
        );
        assert!(
            near_limit.diagonals.offset < start,
            "the offset must have passed INT4_MAX / 4 and been cleared"
        );

        // Another query set of the same table length (the same queries in another order) takes
        // the table over.
        let set_b: Vec<&[u8]> = queries.iter().rev().map(Vec::as_slice).collect();
        let (query_b, offsets_b) = x_concatenated(&set_b);
        let new_b = |table| {
            XTwoHit::with_table(
                table,
                &query_b,
                &offsets_b,
                &x_dropoffs,
                &cutoffs,
                ScoringMatrix::Blosum62,
                window,
                3,
            )
        };
        let cleared_b = cleared_hits(|| new_b(None), &query_b, &units);
        let table = near_limit.into_table();
        let entries = table.entries.as_ptr();
        let mut taken_over = new_b(Some(table));
        assert_eq!(
            taken_over.diagonals.entries.as_ptr(),
            entries,
            "the table must be taken over"
        );
        assert_eq!(
            kept_hits(&mut taken_over, &query_b, &units),
            cleared_b,
            "taken over"
        );

        // A table of another length is not taken over.
        let (short_query, short_offsets) = x_concatenated(&set_a[..1]);
        let other = XTwoHit::with_table(
            Some(taken_over.into_table()),
            &short_query,
            &short_offsets,
            &x_dropoffs[..1],
            &cutoffs[..1],
            ScoringMatrix::Blosum62,
            window,
            3,
        );
        assert_eq!(other.diagonals.offset, window);
        assert_eq!(
            other.diagonals.entries.len(),
            (short_query.len() - 1 + window as usize).next_power_of_two()
        );
    }

    // LOSAT_X_TBNHITBF: `x_pairs_bf` against `pairs` on units made as in the test above (fragments
    // of the queries with substitutions, so that extensions and flag-1 cells occur), on cleared
    // tables, on a kept table and across the clear near INT4_MAX / 4, for windows 40 and
    // wordsize + 1; each unit's offset pairs come in arrays of random lengths, as the scan fills
    // them. The init HSPs and every cell must be equal after every unit.
    #[test]
    fn x_pairs_bf_gives_the_hits_and_table_of_pairs() {
        let mut seed = 0x9e37_79b9_u64;
        let mut next = move |bound: usize| {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((seed >> 33) as usize) % bound
        };
        const RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
        fn random(length: usize, next: &mut dyn FnMut(usize) -> usize) -> Vec<u8> {
            (0..length).map(|_| RESIDUES[next(20)]).collect()
        }
        let queries: Vec<Vec<u8>> = [150, 90, 210]
            .into_iter()
            .map(|length| random(length, &mut next))
            .collect();
        let mut units = Vec::new();
        for _ in 0..60 {
            let mut unit = random(next(40), &mut next);
            for _ in 0..1 + next(3) {
                let query = &queries[next(queries.len())];
                let start = next(query.len() - 30);
                let end = (start + 30 + next(80)).min(query.len());
                for &residue in &query[start..end] {
                    unit.push(if next(100) < 15 {
                        RESIDUES[next(20)]
                    } else {
                        residue
                    });
                }
                unit.extend(random(next(30), &mut next));
            }
            units.push(
                encode_tblastn_lookup_query(&unit, None, false).aa_seq[1..=unit.len()].to_vec(),
            );
        }
        let x_dropoffs = [16, 16, 16];
        let cutoffs = [24, 24, 24];
        let set: Vec<&[u8]> = queries.iter().map(Vec::as_slice).collect();
        let (query, offsets) = x_concatenated(&set);
        let mut found = 0;
        // too far, overlap, flag 1, candidate (the last two go to `pairs` in `x_pairs_bf`).
        let mut outcomes = [0u64; 4];
        for window in [40, 4] {
            for (mode, start) in [
                ("cleared", None),
                ("kept", None),
                ("near the clear", Some(i32::MAX / 4 - 5_000)),
            ] {
                let new = || {
                    XTwoHit::new(
                        &query,
                        &offsets,
                        &x_dropoffs,
                        &cutoffs,
                        ScoringMatrix::Blosum62,
                        window,
                        3,
                    )
                };
                let mut reference = new();
                let mut branch_free = new();
                if let Some(start) = start {
                    reference.diagonals.offset = start;
                    branch_free.diagonals.offset = start;
                }
                for (u, unit) in units.iter().enumerate() {
                    if mode == "cleared" {
                        reference = new();
                        branch_free = new();
                    }
                    let pairs = x_word_pairs(&query, unit);
                    // The reference one hit at a time, its outcomes counted from the cell before.
                    let mut expected = Vec::new();
                    for pair in &pairs {
                        let mask = reference.x_table().len() - 1;
                        let (last_hit, flag) = reference.x_table()
                            [(pair.q_off.wrapping_sub(pair.s_off) as usize) & mask];
                        let diff = pair.s_off as i32 - (last_hit - reference.diagonals.offset);
                        let outcome = if flag {
                            2
                        } else if diff >= window {
                            0
                        } else if diff < 3 {
                            1
                        } else {
                            3
                        };
                        outcomes[outcome] += 1;
                        reference.pairs(1, 0, unit, std::slice::from_ref(pair), &mut expected);
                    }
                    let mut got = Vec::new();
                    let mut rest = &pairs[..];
                    while !rest.is_empty() {
                        let (array, tail) = rest.split_at(1 + next(rest.len().min(50)));
                        branch_free.x_pairs_bf(1, 0, unit, array, &mut got);
                        rest = tail;
                    }
                    assert_eq!(got, expected, "window {window} {mode} unit {u}: init HSPs");
                    assert_eq!(
                        branch_free.x_table(),
                        reference.x_table(),
                        "window {window} {mode} unit {u}: table"
                    );
                    found += expected.len();
                    if mode != "cleared" {
                        reference.exit(unit.len());
                        branch_free.exit(unit.len());
                    }
                }
            }
        }
        assert!(found >= 200, "the units must give hits ({found})");
        assert!(
            outcomes.iter().all(|&count| count > 0),
            "outcomes {outcomes:?}"
        );
    }

    // LOSAT_X_TBNHITBFSHADOW: the shadow keeps the branch-free run's init HSPs and table, which are
    // those of `pairs`.
    #[test]
    fn x_pairs_bf_shadow_keeps_the_hits_of_pairs() {
        let raw: &[u8] = b"MKVLAAGIVGLLLAQWERTYIPASDFGHKLCVNM";
        let (query, offsets) = x_concatenated(&[raw]);
        let mut encoded =
            encode_tblastn_lookup_query(raw, None, false).aa_seq[1..=raw.len()].to_vec();
        encoded.extend_from_within(..);
        let new = || {
            XTwoHit::new(
                &query,
                &offsets,
                &[16],
                &[10],
                ScoringMatrix::Blosum62,
                40,
                3,
            )
        };
        let pairs = x_word_pairs(&query, &encoded);
        let mut reference = new();
        let mut expected = Vec::new();
        reference.pairs(1, 0, &encoded, &pairs, &mut expected);
        let mut shadow = new();
        let mut got = Vec::new();
        shadow.x_pairs_bf_shadow(1, 0, &encoded, &pairs, &mut got);
        assert!(!expected.is_empty());
        assert_eq!(got, expected);
        assert_eq!(shadow.x_table(), reference.x_table());
    }
}
