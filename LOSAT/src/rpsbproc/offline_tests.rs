//! `ParseBlastOutput` and `_TOfflDomQuery` (offl_cd_align_proc.cpp:17-317, quoted at the
//! functions): the getters of OPTIONAL members, the E-value cutoff, the reading frames, the
//! commits and their memory, the segments and residues of `ParseAlignStrings`, libstdc++'s
//! `replace`, and what LOSAT refuses.

use super::*;

/// An HSP of `qseq` and `hseq` on the query from `from` (1-based) in frame `frame`.
fn hsp(from: i32, to: i32, frame: i32, qseq: &[u8], hseq: &[u8]) -> Hsp {
    Hsp {
        bit_score: 40.0,
        evalue: 1e-5,
        query_from: from,
        query_to: to,
        hit_from: 1,
        query_frame: Some(frame),
        identity: Some(1),
        align_len: Some(qseq.len() as i32),
        qseq: qseq.to_vec(),
        hseq: hseq.to_vec(),
    }
}

fn query(length: i32) -> OfflineQuery {
    let mut query = OfflineQuery::default();
    query.set_length(length);
    query
}

fn message(error: &ProcessError) -> String {
    match error {
        ProcessError::Exception(exception) => exception.get_msg().to_string(),
        ProcessError::Std(what) | ProcessError::Unsupported(what) => what.clone(),
    }
}

// The getters of `Hsp_align-len`, `Hsp_identity` and `Hsp_query-frame` throw in this order
// (`CUnassignedMember`, serialobject.cpp:266-306), after the E-value test; a NaN E-value
// passes the cutoff, as does every HSP under a NaN cutoff.
#[test]
fn getters_and_the_cutoff() {
    let mut missing = hsp(5, 7, 0, b"ABC", b"ABC");
    missing.align_len = None;
    missing.identity = None;
    let mut q = query(10);
    q.hsp(1, &missing, 0.01);
    assert_eq!(
        message(q.failure.as_ref().unwrap()),
        "CHsp::GetAlign_len(): Attempt to get unassigned member NCBI-BlastOutput::Hsp.align-len"
    );
    missing.align_len = Some(3);
    let mut q = query(10);
    q.hsp(1, &missing, 0.01);
    assert_eq!(
        message(q.failure.as_ref().unwrap()),
        "CHsp::GetIdentity(): Attempt to get unassigned member NCBI-BlastOutput::Hsp.identity"
    );
    missing.identity = Some(3);
    missing.query_frame = None;
    let mut q = query(10);
    q.hsp(1, &missing, 0.01);
    assert_eq!(
        message(q.failure.as_ref().unwrap()),
        "CHsp::GetQuery_frame(): Attempt to get unassigned member NCBI-BlastOutput::Hsp.query-frame"
    );
    // Above the cutoff, the getters are not reached.
    let mut q = query(10);
    q.hsp(1, &missing, 1e-6);
    assert!(q.failure.is_none() && q.aligns.is_empty());
    let mut nan = hsp(5, 7, 0, b"ABC", b"ABC");
    nan.evalue = f64::NAN;
    let mut q = query(10);
    q.hsp(1, &nan, 0.01);
    assert_eq!(q.aligns.len(), 1);
    let mut q = query(10);
    q.hsp(1, &hsp(5, 7, 0, b"ABC", b"ABC"), f64::NAN);
    assert_eq!(q.aligns.len(), 1);
}

// A protein HSP (`PRCommit`, frame 0): the segments where both strings have residues, the
// query's residues written at their positions (`-` elsewhere), a gap of the query or of the
// domain, and the leading gaps trimmed.
#[test]
fn protein_segments_and_residues() {
    // Query 3..10 (1-based), domain from 1: `--AB-CDE` over `QQAXYC-E`.
    let mut q = query(12);
    q.hsp(7, &hsp(3, 10, 0, b"--AB-CDE", b"QQAXYC-E"), 0.01);
    assert!(q.failure.is_none());
    let align = &q.aligns[0];
    assert!(!align.pr2na && align.reading_frame == 0);
    assert_eq!((align.from, align.to, align.pssm_id), (2, 9, 7));
    // Segments AB, C and E; the domain's start counts from `Hsp_hit-from` whatever the trimmed
    // columns hold, and steps back at a gap of the domain (`--hf`).
    assert_eq!(align.m_starts, vec![2, 4, 6]);
    assert_eq!(align.s_starts, vec![0, 3, 4]);
    assert_eq!(align.lens, vec![2, 1, 1]);
    assert_eq!(q.translated[0], b"--ABCDE-----".to_vec());
    assert!(q.seq_data_is_translated);
}

// `basic_string::replace` past the string's length throws libstdc++'s `out_of_range`; a replace
// that reaches past the end makes the string longer.
#[test]
fn libstdcxx_replace() {
    let mut text = b"----".to_vec();
    string_replace(&mut text, 2, 4, b"ABCDEF", 1, 4).unwrap();
    assert_eq!(text, b"--BCDE".to_vec());
    let error = string_replace(&mut text, -1, 1, b"A", 0, 1).unwrap_err();
    assert_eq!(
        message(&error),
        "basic_string::replace: __pos (which is 18446744073709551615) > this->size() (which is 6)"
    );
    let error = string_replace(&mut text, 7, 1, b"A", 0, 1).unwrap_err();
    assert_eq!(
        message(&error),
        "basic_string::replace: __pos (which is 7) > this->size() (which is 6)"
    );
    // A protein query of length 2 with a hit from residue 5.
    let mut q = query(2);
    q.hsp(1, &hsp(5, 7, 0, b"ABC", b"ABC"), 0.01);
    assert_eq!(
        message(q.failure.as_ref().unwrap()),
        "basic_string::replace: __pos (which is 4) > this->size() (which is 2)"
    );
}

// The frames of `ParseBlastOutput`: a plus frame's index, a minus frame's plus-strand index with
// its aligned length (`m_iTo / 3`) above it; `NACommit` makes the six strings of `-`.
#[test]
fn nucleotide_frames() {
    let mut q = query(20);
    // Frame +2: bases 2..10 (1-based), 3 residues from residue 0 of frame 1.
    q.hsp(5, &hsp(2, 10, 2, b"MKV", b"MKV"), 0.01);
    assert!(q.failure.is_none() && q.is_na);
    let lens: Vec<usize> = q.translated.iter().map(Vec::len).collect();
    assert_eq!(lens, vec![6, 6, 6, 6, 6, 6]);
    let plus = &q.aligns[0];
    assert!(plus.pr2na && plus.reading_frame == 1);
    assert_eq!(q.translated[1], b"MKV---".to_vec());
    // Frame -1: bases 4..12, m_iTo = 11: plus index (20 - 1 + 1) % 3 = 2, aligned length 3.
    q.hsp(6, &hsp(4, 12, -1, b"WY", b"WY"), 0.01);
    let minus = &q.aligns[1];
    assert_eq!(minus.reading_frame, 2 | (3 << 2));
    // Minus frame index (20 - 2) % 3 + 3 = 3; first residue (20 - 11 - 3 + 2) / 3 = 2.
    assert_eq!(q.translated[3], b"--WY--".to_vec());
    assert_eq!(minus.m_starts, vec![2 - (6 - 3)]);
    // A frame of 4 is frame index 3 (`rf - 1`), read as a plus frame of the string of frame -1;
    // from base 2 its first residue is `(1 - 3) / 3` in `unsigned int`, past the string.
    let mut q = query(20);
    q.hsp(5, &hsp(4, 12, 4, b"MKV", b"MKV"), 0.01);
    assert_eq!(q.aligns[0].reading_frame, 3);
    assert_eq!(q.translated[3], b"MKV---".to_vec());
    let mut q = query(20);
    q.hsp(5, &hsp(2, 10, 4, b"MKV", b"MKV"), 0.01);
    assert_eq!(
        message(q.failure.as_ref().unwrap()),
        "basic_string::replace: __pos (which is 1431655764) > this->size() (which is 6)"
    );
}

// The commits' strings: lengths in `unsigned int`; a nucleotide query of 0 or 1 base (more than
// a billion residues per frame) and a protein query of 2^29 residues are refused (NCBI holds
// 1 GiB or more); a query of 2 bases has empty frames, whose first replace throws.
#[test]
fn commits_and_their_memory() {
    for len in [0, 1] {
        let mut q = query(len);
        q.hsp(1, &hsp(1, 3, 1, b"M", b"M"), 0.01);
        assert!(
            matches!(q.failure, Some(ProcessError::Unsupported(_))),
            "{len}"
        );
    }
    let mut q = query(1 << 29);
    q.hsp(1, &hsp(1, 3, 0, b"M", b"M"), 0.01);
    assert!(matches!(q.failure, Some(ProcessError::Unsupported(_))));
    let mut q = query(2);
    q.hsp(1, &hsp(4, 6, 1, b"M", b"M"), 0.01);
    assert!(q.translated.iter().all(Vec::is_empty));
    assert_eq!(
        message(q.failure.as_ref().unwrap()),
        "basic_string::replace: __pos (which is 1) > this->size() (which is 0)"
    );
}

// "See if NA": a commit that leaves `m_dimTranslated[0]` empty (a protein query of length 0)
// decides the kind from the residues of the alignment against its range on the query.
#[test]
fn kind_from_the_alignment() {
    // 3 residues over 9 bases: nucleotide, then the frames of a query of length 0 (refused).
    let mut q = query(0);
    q.hsp(1, &hsp(1, 9, 0, b"MKV", b"MKV"), 0.01);
    assert!(matches!(q.failure, Some(ProcessError::Unsupported(_))));
    // 3 residues over 3: protein again, of length 0; the replace at its end inserts them.
    let mut q = query(0);
    q.hsp(1, &hsp(1, 3, 0, b"MKV", b"MKV"), 0.01);
    assert!(q.failure.is_none() && !q.is_na);
    assert_eq!(q.translated[0], b"MKV".to_vec());
}

// What NCBI leaves undefined: `hseq` read past its terminating 0, a residue written at or past
// its string's length.
#[test]
fn undefined_reads_and_writes() {
    // One letter short: the terminating 0 is a residue (defined).
    let mut q = query(10);
    q.hsp(1, &hsp(1, 3, 0, b"ABC", b"AB"), 0.01);
    assert!(q.failure.is_none());
    assert_eq!(q.aligns[0].lens, vec![3]);
    let mut q = query(10);
    q.hsp(1, &hsp(1, 4, 0, b"ABCD", b"AB"), 0.01);
    assert!(matches!(q.failure, Some(ProcessError::Unsupported(_))));
    // A gap of the domain after the query's end writes past the string.
    let mut q = query(3);
    q.hsp(1, &hsp(1, 4, 0, b"ABCD", b"ABC-"), 0.01);
    assert!(matches!(q.failure, Some(ProcessError::Unsupported(_))));
}

fn refused(q: &OfflineQuery) -> String {
    match &q.failure {
        Some(ProcessError::Unsupported(message)) => message.clone(),
        other => panic!("not a refusal: {other:?}"),
    }
}

// `GetQuery_from() - COORDSBASE`, `GetQuery_to() - COORDSBASE` and `GetHit_from() - COORDSBASE`
// of `-2147483648` overflow `int` (offl_cd_align_proc.cpp:166-167, 234, quoted at `parse_hsp`):
// refused for each of the three members.
#[test]
fn coordinates_of_int_min() {
    let message = "an HSP coordinate of -2147483648 (NCBI's `- COORDSBASE` overflows int, undefined) in the XML input is not supported by LOSAT's rpsbproc";
    for member in 0..3 {
        let mut hit = hsp(5, 7, 0, b"ABC", b"ABC");
        match member {
            0 => hit.query_from = i32::MIN,
            1 => hit.query_to = i32::MIN,
            _ => hit.hit_from = i32::MIN,
        }
        let mut q = query(10);
        q.hsp(1, &hit, 0.01);
        assert_eq!(refused(&q), message, "{member}");
    }
    // `-2147483647` is read: `m_iFrom` is then `INT_MIN`, and NCBI's `replace` at that `size_t`
    // throws.
    let mut q = query(10);
    q.hsp(1, &hsp(-i32::MAX, 7, 0, b"ABC", b"ABC"), 0.01);
    assert_eq!(
        message_of(&q),
        "basic_string::replace: __pos (which is 18446744071562067968) > this->size() (which is 10)"
    );
}

fn message_of(q: &OfflineQuery) -> String {
    message(q.failure.as_ref().unwrap())
}

// "See if NA" (offl_cd_align_proc.cpp:209-228): `m_iTo - m_iFrom + 1` is `int` arithmetic; a span
// of `INT_MAX + 1` overflows (refused). It is reached only while `m_dimTranslated[0]` is empty
// (a protein query of length 0).
#[test]
fn span_that_overflows_int() {
    let mut q = query(0);
    q.hsp(1, &hsp(0, i32::MAX, 0, b"ABC", b"ABC"), 0.01);
    assert_eq!(
        refused(&q),
        "an HSP whose Hsp_query-to - Hsp_query-from + 1 overflows NCBI's int (undefined) in the XML input is not supported by LOSAT's rpsbproc"
    );
}

// `++hf` and `--hf` of `ParseAlignStrings` at `int`'s limits (offl_cd_align_proc.cpp:51-122):
// a hit from `INT_MAX` with a query gap after it, and one whose `hf += segLen` wraps to `INT_MIN`
// before a domain gap. `qf` reaches no limit: `replace` checks it against the string first.
#[test]
fn domain_positions_that_overflow_int() {
    let message = "an alignment whose positions overflow NCBI's int (undefined) in the XML input is not supported by LOSAT's rpsbproc";
    // hf = INT_MAX - 1, then + 1 at the gap, then `++hf` in the run of query gaps.
    let mut up = hsp(1, 2, 0, b"A--B", b"AXYB");
    up.hit_from = i32::MAX;
    let mut q = query(10);
    q.hsp(1, &up, 0.01);
    assert_eq!(refused(&q), message);
    // hf = INT_MAX - 1, + 2 wraps to INT_MIN, then `--hf` at the domain gap.
    let mut down = hsp(1, 3, 0, b"ABC", b"AB-");
    down.hit_from = i32::MAX;
    let mut q = query(10);
    q.hsp(1, &down, 0.01);
    assert_eq!(refused(&q), message);
    // A `qf` past the string throws at the `replace` of its run.
    let mut q = query(10);
    q.hsp(1, &hsp(i32::MAX, i32::MAX, 0, b"A-B", b"AXB"), 0.01);
    assert_eq!(
        message_of(&q),
        "basic_string::replace: __pos (which is 2147483646) > this->size() (which is 10)"
    );
}

// `SeqLen_t alnLen = qseq.size()`: an `Hsp_qseq` of 2^31 letters or more is refused.
#[test]
fn alignment_lengths_beyond_int() {
    assert_eq!(
        alignment_length(i32::MAX as usize).ok(),
        Some(i32::MAX as u32)
    );
    let Err(ProcessError::Unsupported(message)) = alignment_length(1 << 31) else {
        panic!("refused");
    };
    assert_eq!(
        message,
        "an Hsp_qseq of 2^31 or more letters (beyond NCBI's int positions) in the XML input is not supported by LOSAT's rpsbproc"
    );
}

// `m_dimTranslated[rfidx][qf] = '\0'` at the string's length is the defined write of the
// terminator (C++17/20 [string.access]); another value there is undefined (refused).
#[test]
fn a_nul_written_at_the_end() {
    let mut text = b"--".to_vec();
    assert!(write_residue(&mut text, 2, 0).is_ok());
    assert_eq!(text, b"--".to_vec());
    assert!(matches!(
        write_residue(&mut text, 2, b'A'),
        Err(ProcessError::Unsupported(_))
    ));
    assert!(matches!(
        write_residue(&mut text, 3, 0),
        Err(ProcessError::Unsupported(_))
    ));
    // A NUL from a CDATA section of `Hsp_qseq`, written over the terminator of a query of 3.
    let mut q = query(3);
    q.hsp(1, &hsp(1, 4, 0, b"ABC\0", b"ABC-"), 0.01);
    assert!(q.failure.is_none());
    assert_eq!(q.translated[0], b"ABC".to_vec());
}

/// The data files of the fixtures (`tests/fixtures/rpsbproc/data`).
fn fixture_data() -> CddData {
    let file = |name: &str| {
        format!(
            "{}/tests/fixtures/rpsbproc/data/{name}",
            env!("CARGO_MANIFEST_DIR")
        )
    };
    CddData::load(&crate::rpsbproc::data::DataFiles {
        cddids: file("cddid.tbl"),
        spfeats: file("cddannot.dat"),
        genfeats: file("cddannot_generic.dat"),
        bscores: file("bitscore_specific.txt"),
        cdtrack: file("cdtrack.txt"),
        famlinks: file("family_superfamily_links"),
    })
    .unwrap()
}

// A kept HSP whose E-value or bit score is NaN (`Hsp_evalue` `nan`: NCBI's `HUGE_VAL/HUGE_VAL`,
// sign set) is classified as NCBI's `Calculate` does: its sort of one hit calls no comparison,
// and with a NaN bit score `m_dBitScore >= m_dMinBitScore` is false (PS16d review O1
// `s_doc2nan`: `-nan` on the DOMAINS line, exit 0).
#[test]
fn nan_scores_are_classified() {
    use crate::ncbi_str::NCBI_NAN;
    let data = fixture_data();
    let q_seq = b"VLLPPGPFTRRQAQAVTTTYSNITLEDDQGSHFRLVVRDTEGRMVWRAWNFEPDAGEGLNRYIRTSGIRT";
    for (e_value, bit_score) in [(NCBI_NAN, 133.904), (5e-43, NCBI_NAN)] {
        let mut hit = hsp(5, 74, 0, q_seq, q_seq);
        hit.evalue = e_value;
        hit.bit_score = bit_score;
        let mut q = query(79);
        q.hsp(488084, &hit, 0.01);
        let dom = q.finish(&data).unwrap();
        assert_eq!(dom.frames[0].concise, vec![0]);
        let align = &dom.aligns[0];
        assert!(align.rep);
        assert_eq!(align.spec_qualified, !bit_score.is_nan());
        assert!(align.e_value.to_bits() == e_value.to_bits());
    }
}
