use super::*;

/// The offsets of three profiles of 22, 18 and 20 columns (each with its sentinel row): `R` 63.
const OFFSETS: [i32; 4] = [0, 23, 42, 63];

/// The 63 rows, each row's scores its index.
fn rows() -> Vec<[i32; BLASTAA_SIZE]> {
    (0..63).map(|i| [i; BLASTAA_SIZE]).collect()
}

fn letters(n: usize) -> Vec<u8> {
    (0..n).map(|i| (i % 20) as u8 + 1).collect()
}

// blast_traceback.c:968-982 (cited at `rows_in_pointers`): the row pointers end at entry R, the
// row after the rows; the rescale reads L_pin of them from the profile's offset (read A).
#[test]
fn the_rows_from_the_offset_end_at_the_row_after_the_rows() {
    let rows = rows();
    let after = RowAfter::new([Some(7); BLASTAA_SIZE]);
    // c3 with 21 letters for 20 columns: its own rows and its sentinel row.
    let own = letters(21);
    let profile = TracebackProfile::new(2, 0, &OFFSETS, &own, true, "RPSBLAST");
    let read = profile.rows(&rows, &after, 28).unwrap();
    assert!(matches!(read, Cow::Borrowed(_)));
    assert_eq!(read.len(), 21);
    assert_eq!(read[20], [62; BLASTAA_SIZE]);
    // 22 letters: entry R, NCBI's row after the rows.
    let longer = letters(22);
    let profile = TracebackProfile::new(2, 0, &OFFSETS, &longer, true, "RPSBLAST");
    let read = profile.rows(&rows, &after, 28).unwrap();
    assert_eq!(read.len(), 22);
    assert_eq!(read[21], [7; BLASTAA_SIZE]);
    assert_eq!(profile.fault(), None);
    // Its letter 27 past the file's last page: a 28-column read faults, a 26-column one does not.
    let mut past = [Some(7); BLASTAA_SIZE];
    past[27] = None;
    let past = RowAfter::new(past);
    assert!(profile.rows(&rows, &past, 26).is_ok());
    let error = profile.rows(&rows, &past, 28).unwrap_err().to_string();
    assert_eq!(profile.fault(), Some(ProfileRead::Rows));
    assert!(error.contains("traceback reads past its rows"), "{error}");
    // 23 letters: one pointer past the R + 1 entries.
    let past_pointers = letters(23);
    let profile = TracebackProfile::new(2, 0, &OFFSETS, &past_pointers, true, "RPSBLAST");
    assert!(profile.rows(&rows, &after, 28).is_err());
    assert_eq!(profile.fault(), Some(ProfileRead::Rows));
    // Profile 0 with 24 letters reads profile 1's first row: in the rows.
    let into_next = letters(24);
    let profile = TracebackProfile::new(0, 0, &OFFSETS, &into_next, true, "RPSBLAST");
    let read = profile.rows(&rows, &after, 28).unwrap();
    assert_eq!(read[23], [23; BLASTAA_SIZE]);
}

// composition_adjustment.c:1162-1199 (cited at `composition`): the first letter of 28 or more
// (read B), the text of LOSAT's refusal, which names the profile's OID in the database (profile 1
// of a volume whose first OID is 3).
#[test]
fn a_letter_past_the_composition_tables_is_refused_with_the_profiles_numbers() {
    let packed = [143u8, 18, 159, 6, 210];
    let profile = TracebackProfile::new(1, 3, &OFFSETS, &packed, true, "RPSTBLASTN");
    assert_eq!(profile.length(), 5);
    let error = profile.composition().unwrap_err().to_string();
    assert_eq!(profile.fault(), Some(ProfileRead::Composition(143)));
    assert_eq!(
        error,
        "a hit of the profile of OID 4, whose consensus has 5 letters in .psq for 18 columns of scores in .rps (makeprofiledb writes it for a scored profile whose seq-data is nucleotide-coded or has more or fewer letters than its length), where NCBI BLAST+ 2.17.0's traceback reads letter 143 past its composition tables, is not supported by LOSAT's RPSTBLASTN"
    );
    // Letters below 28 pass (c1 as 22 A: 00 00 00 00 00 02).
    let below = [0u8, 0, 0, 0, 0, 2];
    let profile = TracebackProfile::new(0, 0, &OFFSETS, &below, true, "RPSBLAST");
    assert!(profile.composition().is_ok());
    assert_eq!(profile.fault(), None);
    // The volume's scan (`has_wide`, once per volume): a volume whose letters are all below 28
    // asks no profile's letters again.
    assert!(!has_wide(&[0, 27, 5, 0, 1]));
    assert!(has_wide(&[0, 27, 28, 0]));
    assert!(has_wide(&packed));
    let profile = TracebackProfile::new(0, 0, &OFFSETS, &below, false, "RPSBLAST");
    assert!(profile.composition().is_ok());
}

// The redo's reads C and D (`RedoRun::with_profile_reads` asks the profile's `trace_reads`) and
// the first read kept.
#[test]
fn the_redos_reads_name_the_copy_or_the_start_and_the_first_stays() {
    let own = letters(10);
    let profile = TracebackProfile::new(0, 0, &OFFSETS, &own, true, "RPSBLAST");
    let reads = profile.trace_reads(ProfileRead::Start);
    assert!(profile.check().is_ok());
    reads.copy_read_past();
    let error = profile.check().unwrap_err().to_string();
    assert!(
        error.contains("10 letters in .psq for 22 columns"),
        "{error}"
    );
    assert!(
        error.contains("traceback reads past its copy of the consensus"),
        "{error}"
    );
    // After the first read no alignment is made, and the first read stays.
    assert!(!reads.starts(0));
    assert!(!reads.starts(10));
    assert_eq!(profile.fault(), Some(ProfileRead::Copy));
    let profile = TracebackProfile::new(0, 0, &OFFSETS, &own, true, "RPSBLAST");
    let reads = profile.trace_reads(ProfileRead::Start);
    assert!(!reads.starts(10));
    reads.copy_read_past();
    let error = profile.check().unwrap_err().to_string();
    assert!(
        error.contains("traceback starts an alignment past the consensus"),
        "{error}"
    );
}

// blast_traceback.c:509-512 (cited at `ProfileTraceReads::starts`), read D without
// composition-based statistics: rpsblast's matrix of L_pin rows (a start at or past L_pin; the
// first read kept) and rpstblastn's mapped rows (the start at or past L_pin, its own arm).
#[test]
fn an_alignment_from_l_pin_on_is_refused() {
    let own = letters(6);
    let profile = TracebackProfile::new(0, 0, &OFFSETS, &own, true, "RPSBLAST");
    let reads = profile.trace_reads(ProfileRead::Start);
    assert!(reads.starts(5));
    assert_eq!(profile.fault(), None);
    assert!(!reads.starts(6));
    assert_eq!(profile.fault(), Some(ProfileRead::Start));
    assert!(!reads.starts(0));
    let error = profile.check().unwrap_err().to_string();
    assert!(error.contains("traceback starts an alignment past the consensus, is not supported by LOSAT's RPSBLAST"), "{error}");
    let profile = TracebackProfile::new(0, 0, &OFFSETS, &own, true, "RPSTBLASTN");
    let reads = profile.trace_reads(ProfileRead::StartInMapping);
    assert!(!reads.starts(6));
    assert_eq!(profile.fault(), Some(ProfileRead::StartInMapping));
    let error = profile.check().unwrap_err().to_string();
    assert!(error.contains("traceback starts an alignment past the consensus, is not supported by LOSAT's RPSTBLASTN"), "{error}");
}

// blast_traceback.c:1194-1198 (cited at `mapped_rows`): rpstblastn's rows from the offset run
// through the volume's rows, NCBI's row after them, then past the pointers (read A).
#[test]
fn the_mapped_rows_read_through_the_row_after_the_rows() {
    let rows = rows();
    let own = letters(20);
    let profile = TracebackProfile::new(2, 0, &OFFSETS, &own, true, "RPSTBLASTN");
    let after = RowAfter::new([Some(7); BLASTAA_SIZE]);
    let mapped = profile.mapped_rows(&rows, &after);
    assert_eq!(mapped.len(), 20);
    assert_eq!(mapped.row(0), &[42; BLASTAA_SIZE]);
    assert_eq!(mapped.row(20), &[62; BLASTAA_SIZE]);
    assert_eq!(mapped.row(21), &[7; BLASTAA_SIZE]);
    // A row after the rows inside the page is not read letter by letter.
    assert!(!mapped.reads_by_letter(mapped.row(21)));
    assert_eq!(profile.fault(), None);
    assert_eq!(mapped.row(22), &[0; BLASTAA_SIZE]);
    assert_eq!(profile.fault(), Some(ProfileRead::Rows));
}

// blast_gapalign.c:548-578 (cited at `MappedRows::reads_by_letter`): the row after the rows whose
// letter 3 lies past the `.rps` file's last page faults only where the alignment reads letter 3;
// its other letters are the file's.
#[test]
fn the_row_after_the_rows_past_the_page_is_read_letter_by_letter() {
    let rows = rows();
    let own = letters(20);
    let mut straddle = [Some(7); BLASTAA_SIZE];
    straddle[3] = None;
    let after = RowAfter::new(straddle);
    let profile = TracebackProfile::new(2, 0, &OFFSETS, &own, true, "RPSTBLASTN");
    let mapped = profile.mapped_rows(&rows, &after);
    let row = mapped.row(21);
    assert_eq!(row[2], 7);
    assert_eq!(row[3], 0);
    assert!(mapped.reads_by_letter(row));
    assert!(!mapped.reads_by_letter(mapped.row(20)));
    mapped.letters_read(row, [1u8, 2, 4, 25].into_iter());
    assert_eq!(profile.fault(), None);
    mapped.letters_read(row, [1u8, 3].into_iter());
    assert_eq!(profile.fault(), Some(ProfileRead::Rows));
    let error = profile.check().unwrap_err().to_string();
    assert!(error.contains("traceback reads past its rows"), "{error}");
}

// The gapped alignment through the row after the rows (profile 2 with 22 letters: 21 rows of the
// volume from its offset, then row R): letters inside the page give the alignment of the whole
// row; a band with letter 3 (past the page) is read A.
#[test]
fn an_alignment_through_a_straddling_row_after_the_rows_reads_its_band_letters() {
    use crate::algorithm::pssm_search::gapalign::{
        pssm_gapped_alignment_with_traceback, PssmGapScratch,
    };
    let rows = rows();
    let own = letters(22);
    let mut straddle = [Some(7); BLASTAA_SIZE];
    straddle[3] = None;
    let align = |after: &RowAfter, subject: &[u8]| {
        let profile = TracebackProfile::new(2, 0, &OFFSETS, &own, true, "RPSTBLASTN");
        let mapped = profile.mapped_rows(&rows, after);
        let mut scratch = PssmGapScratch::new();
        let traced = pssm_gapped_alignment_with_traceback(
            &mapped,
            subject,
            subject.len(),
            18,
            2,
            11,
            1,
            30,
            &mut scratch,
        );
        (traced, profile.fault())
    };
    let subject = [1u8; 12];
    let (whole, fault) = align(&RowAfter::new([Some(7); BLASTAA_SIZE]), &subject);
    assert_eq!(fault, None);
    // The alignment reaches the row after the rows (profile position 21).
    assert_eq!(whole.query_stop, 22);
    let (inside, fault) = align(&RowAfter::new(straddle), &subject);
    assert_eq!(fault, None);
    assert_eq!(inside, whole);
    let (_, fault) = align(&RowAfter::new(straddle), &[3u8; 12]);
    assert_eq!(fault, Some(ProfileRead::Rows));
}

// blast_hits.c:796-812 (cited at `ProfileTraceReads::traced`), read E: a profile letter of 28 or
// more in a substituted column; inserted columns (the profile's letters against a gap) and
// deleted ones read no matrix row.
#[test]
fn a_substituted_letter_past_the_score_matrix_is_refused() {
    let packed = [143u8, 18, 159, 6, 210];
    let profile = TracebackProfile::new(1, 0, &OFFSETS, &packed, true, "RPSTBLASTN");
    let reads = profile.trace_reads(ProfileRead::StartInMapping);
    assert!(reads.traced(1, &[GapEditOp::Sub(1)]));
    assert!(reads.traced(1, &[GapEditOp::Del(2), GapEditOp::Sub(1)]));
    assert!(reads.traced(
        1,
        &[GapEditOp::Sub(1), GapEditOp::Ins(1), GapEditOp::Sub(1)]
    ));
    assert_eq!(profile.fault(), None);
    assert!(!reads.traced(
        1,
        &[GapEditOp::Sub(1), GapEditOp::Del(2), GapEditOp::Sub(2)]
    ));
    assert_eq!(profile.fault(), Some(ProfileRead::Substitution(159)));
    let error = profile.check().unwrap_err().to_string();
    assert!(
        error.contains("traceback reads letter 159 past its score matrix"),
        "{error}"
    );
    // A profile of letters below 28 never walks its edit scripts.
    let own = letters(22);
    let profile = TracebackProfile::new(0, 0, &OFFSETS, &own, true, "RPSBLAST");
    let reads = profile.trace_reads(ProfileRead::Start);
    assert!(reads.traced(0, &[GapEditOp::Sub(22)]));
    assert_eq!(profile.fault(), None);
}
