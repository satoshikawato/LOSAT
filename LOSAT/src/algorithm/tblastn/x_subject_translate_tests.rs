use super::super::search_seed::x_scan_plan;
use super::*;
use crate::algorithm::tblastx::translation::generate_frames;
use crate::config::ScoringMatrix;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:95-103 and
// c++/src/algo/blast/api/blast_objmgr_tools.cpp:427-474,515-520
// else { int pick = random.GetRand() % bitcount; /* pick a set bit */ }
// The uncompressed ncbi4na sequence is retained by SetupSubjects_OMF for
// traceback reevaluation; these resolved bases are preliminary-search input.
fn resolve_local_subject_ncbi2na(subject: &[u8]) -> Result<Vec<u8>> {
    u32::try_from(subject.len())?;
    let ncbi4na = subject
        .iter()
        .map(|&base| {
            // TBLASTN reads the IUPAC letters without `U`.
            crate::core::blast_encoding::iupacna_to_ncbi4na(base)
                .filter(|_| !base.eq_ignore_ascii_case(&b'U'))
                .ok_or_else(|| {
                    anyhow::anyhow!("TBLASTN subject base '{}' is not implemented", base as char)
                })
        })
        .collect::<Result<Vec<u8>>>()?;
    Ok(
        crate::core::blast_encoding::resolve_ncbi4na_to_ncbi2na(&ncbi4na)
            .into_iter()
            .map(|code| b"ACGT"[code as usize])
            .collect(),
    )
}

// No NCBI counterpart: the three helpers below prepare the inputs of the tests; they compute
// no value NCBI computes.
/// The genetic codes of NCBI's `gc.prt` that LOSAT knows.
fn genetic_codes() -> Vec<u8> {
    let codes: Vec<u8> = (0..=u8::MAX)
        .filter(|&id| GeneticCode::try_from_id(id).is_ok())
        .collect();
    assert_eq!(codes.len(), 27);
    codes
}

/// The reference plan of one subject: the letter-wise translation of
/// `tblastx::translation::generate_frames` (after `resolve_local_subject_ncbi2na`), which the
/// tables must reproduce.
fn reference_plan(
    subject: &[u8],
    db_gencode: u8,
    mask_lowercase: bool,
) -> Result<Arc<XSubjectPlan>> {
    let code = GeneticCode::try_from_id(db_gencode).map_err(anyhow::Error::msg)?;
    let resolved_subject = resolve_local_subject_ncbi2na(subject)?;
    let frames = generate_frames(&resolved_subject, &code);
    let negative_first_length = frames
        .iter()
        .find(|frame| frame.frame == -1)
        .map(|frame| frame.aa_len)
        .unwrap_or(0);
    let mut units = Vec::new();
    for (index, frame) in frames.iter().enumerate() {
        for (chunk, ranges) in translated_chunk_scan_ranges(
            subject,
            frame.frame,
            frame.aa_len,
            negative_first_length,
            mask_lowercase,
        )? {
            units.push((index, chunk, ranges));
        }
    }
    Ok(Arc::new(XSubjectPlan { frames, units }))
}

// No NCBI counterpart: the comparison of two plans for the tests; it computes no value NCBI
// computes.
/// What differs between the plan from the tables and the plan of the letter-wise reference
/// (`reference_plan`) of one subject, `None` when nothing does. Two errors are the same when
/// their texts are.
fn x_plan_difference(
    tables: &Result<XSubjectPlan>,
    reference: &Result<Arc<XSubjectPlan>>,
) -> Option<String> {
    let (tables, reference) = match (tables, reference) {
        (Ok(tables), Ok(reference)) => (tables, reference),
        (Err(tables), Err(reference)) => {
            let (tables, reference) = (format!("{tables:#}"), format!("{reference:#}"));
            return (tables != reference)
                .then(|| format!("error \"{tables}\" for the reference's \"{reference}\""));
        }
        (Ok(_), Err(reference)) => {
            return Some(format!(
                "a plan for the reference's error \"{reference:#}\""
            ))
        }
        (Err(tables), Ok(_)) => {
            return Some(format!("error \"{tables:#}\" for a plan of the reference"))
        }
    };
    x_plan_fields_difference(tables, reference)
}

/// `length` letters of `alphabet`, drawn by a fixed linear congruential sequence.
fn letters(alphabet: &[u8], length: usize, seed: u64) -> Vec<u8> {
    let mut state = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            alphabet[(state >> 33) as usize % alphabet.len()]
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:441-450
// Every byte as a subject of one letter: the ncbi2na value is the letter the reference path
// resolves it to, and a byte that is no subject letter (`U` among them) gives the reference
// path's error text.
#[test]
fn x_every_byte_resolves_as_in_the_reference_path() {
    let mut letters = 0;
    for byte in 0..=u8::MAX {
        let subject = [byte];
        match (
            x_subject_ncbi2na(&subject),
            resolve_local_subject_ncbi2na(&subject),
        ) {
            (Ok(ncbi2na), Ok(resolved)) => {
                assert_eq!(ncbi2na.len(), 1);
                assert_eq!(b"ACGT"[ncbi2na[0] as usize], resolved[0], "byte {byte}");
                letters += 1;
            }
            (Err(ours), Err(theirs)) => {
                assert_eq!(ours.to_string(), theirs.to_string(), "byte {byte}");
            }
            (ours, theirs) => panic!("byte {byte}: {ours:?} for the reference's {theirs:?}"),
        }
    }
    // Fifteen IUPAC letters in either case and the gap.
    assert_eq!(letters, 31);
    assert!(x_subject_ncbi2na(b"ACGU").is_err());
    assert!(x_subject_ncbi2na(b"acgu").is_err());
    // The error names the first byte of the subject that is no letter, also after letters
    // that draw and before other such bytes.
    for (subject, named) in [
        (&b"A!A?"[..], '!'),
        (b"?A!", '?'),
        (b"NNRY!ACGT?", '!'),
        (b"ACGTNu*", 'u'),
        (b"-N-\xE9!", '\u{e9}'),
    ] {
        let expected = format!("TBLASTN subject base '{named}' is not implemented");
        let ours = x_subject_ncbi2na(subject).expect_err("no letter");
        let theirs = resolve_local_subject_ncbi2na(subject).expect_err("no letter");
        assert_eq!(ours.to_string(), expected);
        assert_eq!(theirs.to_string(), expected);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:999-1002,1022-1037
// The two tables of every genetic code against this loop over the genetic-code string
// (`GeneticCode::table`, in ncbieaa), and against frames 1 and -1 of the reference path's
// `generate_frames` on the codon as letters.
#[test]
fn x_codon_tables_are_ncbis_translation_tables_for_every_genetic_code() {
    let mapping = [2usize, 1, 3, 0];
    for db_gencode in genetic_codes() {
        let code = GeneticCode::try_from_id(db_gencode).unwrap();
        let tables = x_codon_tables(db_gencode).unwrap();
        let mut forward = [0u8; 64];
        let mut reverse = [0u8; 64];
        for index1 in 0..4usize {
            for index2 in 0..4usize {
                for index3 in 0..4usize {
                    let (bp1, bp2, bp3) = (3 - index1, 3 - index2, 3 - index3);
                    let codon = (mapping[bp1] << 4) + (mapping[bp2] << 2) + mapping[bp3];
                    reverse[(index3 << 4) + (index2 << 2) + index1] =
                        aa_char_to_ncbistdaa(code.table[codon]);
                    let codon = (mapping[index1] << 4) + (mapping[index2] << 2) + mapping[index3];
                    forward[(index1 << 4) + (index2 << 2) + index3] =
                        aa_char_to_ncbistdaa(code.table[codon]);
                }
            }
        }
        assert_eq!(tables.forward, forward, "genetic code {db_gencode}");
        assert_eq!(tables.reverse, reverse, "genetic code {db_gencode}");
        for index in 0..64usize {
            let codon = [
                b"ACGT"[index >> 4],
                b"ACGT"[(index >> 2) & 3],
                b"ACGT"[index & 3],
            ];
            let frames = generate_frames(&codon, &code);
            assert_eq!(frames.len(), 2);
            assert_eq!(
                (frames[0].frame, frames[0].aa_seq[1]),
                (1, tables.forward[index]),
                "genetic code {db_gencode}, codon {index}"
            );
            assert_eq!(
                (frames[1].frame, frames[1].aa_seq[1]),
                (-1, tables.reverse[index]),
                "genetic code {db_gencode}, codon {index}"
            );
        }
    }
    assert_eq!(
        x_codon_tables(0).err().map(|error| error.to_string()),
        GeneticCode::try_from_id(0).err()
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1080-1091,1100-1102
// The plan from the tables against the reference plan: every length from 0 to 40 bases (the
// subjects too short for a frame, and every remainder of the three frames on both strands)
// and longer ones, of the four bases, of the bases with lowercase runs, of every IUPAC
// letter with gaps (one draw each, in subject order), and of `N` only; with and without
// lowercase masking, for every genetic code on the long subjects.
#[test]
fn x_table_plans_are_the_reference_plans() {
    let alphabets: [&[u8]; 4] = [
        b"ACGT",
        b"ACGTACGTACGTacgtn",
        b"ACGTMRSVWYHKDBN-acgtmrsvwyhkdbn",
        b"N",
    ];
    let mut compared = 0usize;
    let mut residues = 0usize;
    let mut scan_plans = 0usize;
    // One lookup table for every call of `x_scan_plan` (as `XQuerySide` shares it).
    let query: &[u8] = b"MKVLAAGIVGLLLAQWSTRE";
    let lookup = OnceLock::new();
    for (which, alphabet) in alphabets.iter().enumerate() {
        let lengths = (0..=40).chain([99, 100, 101, 1_000, 4_999, 20_003]);
        for length in lengths {
            let subject = letters(alphabet, length, (which * 100_000 + length) as u64);
            let codes = if length >= 1_000 {
                genetic_codes()
            } else {
                vec![1, 4, 11]
            };
            for db_gencode in codes {
                for mask_lowercase in [false, true] {
                    let tables = x_table_plan(&subject, db_gencode, mask_lowercase);
                    let reference = reference_plan(&subject, db_gencode, mask_lowercase);
                    assert_eq!(
                        x_plan_difference(&tables, &reference),
                        None,
                        "alphabet {which}, {length} bases, genetic code {db_gencode}, lowercase masking {mask_lowercase}"
                    );
                    let plan = tables.unwrap();
                    assert_eq!(plan.frames.len(), 2 * length.saturating_sub(2).min(3));
                    residues += plan.frames.iter().map(|f| f.aa_len).sum::<usize>();
                    compared += 1;
                    // The plan `x_scan_plan` makes is the reference's too.
                    let scan = x_scan_plan(
                        &[query],
                        crate::sequence::packed_subject::SubjectNt::Letters(&subject).into(),
                        db_gencode,
                        None,
                        11,
                        mask_lowercase,
                        ScoringMatrix::Blosum62,
                        3,
                        &lookup,
                        &Default::default(),
                        None,
                        None,
                    )
                    .unwrap();
                    assert_eq!(scan.is_some(), length >= 11);
                    if let Some(scan) = scan {
                        let copy =
                            reference.map(|plan| Arc::try_unwrap(plan).ok().expect("one owner"));
                        assert_eq!(
                            x_plan_difference(&copy, &Ok(scan.subject)),
                            None,
                            "x_scan_plan: alphabet {which}, {length} bases, genetic code {db_gencode}, lowercase masking {mask_lowercase}"
                        );
                        scan_plans += 1;
                    }
                }
            }
        }
    }
    assert_eq!(compared, 4 * (44 * 3 + 3 * 27) * 2);
    assert_eq!(scan_plans, 4 * (33 * 3 + 3 * 27) * 2);
    assert!(residues > 1_000_000);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:431,461
// The generator is seeded with the length of the subject and drawn from once per gap or
// ambiguous letter, in subject order. A subject whose only such letter is its last one, and
// subjects that differ in length only, draw as the reference path does.
#[test]
fn x_draws_follow_the_seed_and_the_order_of_the_reference_path() {
    for length in [1usize, 2, 3, 33, 34, 257, 1_000] {
        for last in *b"NRYKMSWBDHV-" {
            let mut subject = vec![b'A'; length];
            subject[length - 1] = last;
            let ncbi2na = x_subject_ncbi2na(&subject).unwrap();
            let resolved = resolve_local_subject_ncbi2na(&subject).unwrap();
            let ours: Vec<u8> = ncbi2na.iter().map(|&code| b"ACGT"[code as usize]).collect();
            assert_eq!(ours, resolved, "{length} bases ending in {}", last as char);
        }
    }
    let subject = letters(b"NRYKMSWBDHV-", 50_000, 7);
    let ncbi2na = x_subject_ncbi2na(&subject).unwrap();
    let resolved = resolve_local_subject_ncbi2na(&subject).unwrap();
    let ours: Vec<u8> = ncbi2na.iter().map(|&code| b"ACGT"[code as usize]).collect();
    assert_eq!(ours, resolved);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775
// `x_subject_plan` gives the reference plan, and `x_scan_plan` gives it too. A subject with a
// byte that is no letter gives the reference path's error.
#[test]
fn x_subject_plans_agree_with_the_reference_and_the_scan_plan() {
    let subject = letters(b"ACGTACGTACGTacgtnRY", 3_001, 11);
    let lookup = OnceLock::new();
    for (db_gencode, mask_lowercase) in [(1, false), (1, true), (4, true), (11, false)] {
        let expected = reference_plan(&subject, db_gencode, mask_lowercase);
        let plan = x_subject_plan(&subject, db_gencode, mask_lowercase)
            .map(|plan| Arc::try_unwrap(plan).ok().expect("one owner"));
        assert_eq!(x_plan_difference(&plan, &expected), None);
        let query: &[u8] = b"MKVLAAGIVGLLLAQWSTRE";
        let scan = x_scan_plan(
            &[query],
            crate::sequence::packed_subject::SubjectNt::Letters(&subject).into(),
            db_gencode,
            None,
            11,
            mask_lowercase,
            ScoringMatrix::Blosum62,
            3,
            &lookup,
            &Default::default(),
            None,
            None,
        )
        .unwrap()
        .expect("a subject long enough to search");
        let tables = x_table_plan(&subject, db_gencode, mask_lowercase);
        assert_eq!(x_plan_difference(&tables, &Ok(scan.subject)), None);
    }
    let mut bad = subject.clone();
    bad[1_500] = b'u';
    let expected = reference_plan(&bad, 1, false)
        .err()
        .expect("no letter")
        .to_string();
    assert_eq!(expected, "TBLASTN subject base 'u' is not implemented");
    let error = x_subject_plan(&bad, 1, false).err().expect("no letter");
    assert_eq!(error.to_string(), expected);
}

// No NCBI counterpart: the comparison helper of these tests. It must see a changed residue, a
// changed field of a frame, a missing frame, changed units and a different error.
#[test]
fn x_plan_difference_sees_every_part_of_a_plan() {
    let subject = letters(b"ACGTacgt", 300, 3);
    let reference = reference_plan(&subject, 1, true);
    let table_plan = || x_table_plan(&subject, 1, true).unwrap();
    assert_eq!(x_plan_difference(&Ok(table_plan()), &reference), None);

    let mut plan = table_plan();
    plan.frames[4].aa_seq[17] ^= 1;
    let difference = x_plan_difference(&Ok(plan), &reference).expect("a residue");
    assert!(difference.contains("frame -2") && difference.contains("byte 17"));

    let mut plan = table_plan();
    plan.frames[2].aa_len += 1;
    assert!(x_plan_difference(&Ok(plan), &reference).is_some());
    let mut plan = table_plan();
    plan.frames[0].orig_len += 1;
    assert!(x_plan_difference(&Ok(plan), &reference).is_some());
    let mut plan = table_plan();
    plan.frames[5].frame = 3;
    assert!(x_plan_difference(&Ok(plan), &reference).is_some());
    let mut plan = table_plan();
    plan.frames[1].aa_seq_nomask = Some(Vec::new());
    assert!(x_plan_difference(&Ok(plan), &reference).is_some());
    let mut plan = table_plan();
    plan.frames[1].seg_masks.push((0, 1));
    assert!(x_plan_difference(&Ok(plan), &reference).is_some());
    let mut plan = table_plan();
    plan.frames.pop();
    assert!(x_plan_difference(&Ok(plan), &reference).is_some());
    let mut plan = table_plan();
    assert!(
        plan.units[3].2.len() > 1,
        "lowercase runs split the scan ranges"
    );
    plan.units[3].2[0].1 += 1;
    assert!(x_plan_difference(&Ok(plan), &reference).is_some());
    let mut plan = table_plan();
    plan.units.pop();
    assert!(x_plan_difference(&Ok(plan), &reference).is_some());

    let failed = || -> Result<XSubjectPlan> { bail!("one") };
    assert!(x_plan_difference(&failed(), &reference).is_some());
    assert!(x_plan_difference(&Ok(table_plan()), &Err(anyhow::anyhow!("one"))).is_some());
    assert!(x_plan_difference(&failed(), &Err(anyhow::anyhow!("two"))).is_some());
    assert_eq!(
        x_plan_difference(&failed(), &Err(anyhow::anyhow!("one"))),
        None
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1080-1090
// The plan of a packed subject against the reference plan of its letters without lowercase
// masking: every length from 0 to 40 bases, lengths around the blocks of `x_packed_frames`
// (every remainder of a block and of a frame), of the four bases, with lowercase runs, of
// every IUPAC letter with gaps, and of `N` only; for every genetic code on the long subjects.
// And as `x_scan_plan` makes it from the store.
#[test]
fn x_packed_plans_are_the_reference_plans() {
    use crate::sequence::packed_subject::{PackedNucleotide, SubjectNt};
    let alphabets: [&[u8]; 4] = [
        b"ACGT",
        b"ACGTACGTACGTacgtn",
        b"ACGTMRSVWYHKDBN-acgtmrsvwyhkdbn",
        b"N",
    ];
    let query: &[u8] = b"MKVLAAGIVGLLLAQWSTRE";
    let lookup = OnceLock::new();
    let mut compared = 0usize;
    for (which, alphabet) in alphabets.iter().enumerate() {
        let around_blocks = [1, 2, 3]
            .into_iter()
            .flat_map(|blocks| (0..7).map(move |extra| blocks * (3 << 14) + extra - 3));
        let lengths = (0..=40)
            .chain([99, 100, 101, 1_000, 4_999])
            .chain(around_blocks);
        for length in lengths {
            let subject = letters(alphabet, length, (which * 100_000 + length) as u64);
            let store = PackedNucleotide::from_subject_letters(&subject).expect("letters");
            let codes = if (1_000..10_000).contains(&length) {
                genetic_codes()
            } else {
                vec![1, 4, 11]
            };
            for db_gencode in codes {
                let reference = reference_plan(&subject, db_gencode, false);
                let packed = x_packed_plan(&store, db_gencode);
                assert_eq!(
                    x_plan_difference(&packed, &reference),
                    None,
                    "alphabet {which}, {length} bases, genetic code {db_gencode}"
                );
                compared += 1;
                if db_gencode != 11 {
                    continue;
                }
                let scan = x_scan_plan(
                    &[query],
                    SubjectNt::Packed(&store).into(),
                    db_gencode,
                    None,
                    11,
                    false,
                    ScoringMatrix::Blosum62,
                    3,
                    &lookup,
                    &Default::default(),
                    None,
                    None,
                )
                .unwrap();
                assert_eq!(scan.is_some(), length >= 11);
                if let Some(scan) = scan {
                    let reference =
                        reference.map(|plan| Arc::try_unwrap(plan).ok().expect("one owner"));
                    assert_eq!(x_plan_difference(&reference, &Ok(scan.subject)), None);
                }
            }
        }
    }
    assert_eq!(compared, 4 * (44 * 3 + 2 * 27 + 21 * 3));
    // a packed subject with lower-case letters is never searched with lowercase masking (the
    // letters are given); one without them (every `-db` sequence) gives the plan of the
    // search without the mask
    let scan = |store: &PackedNucleotide, mask_lowercase: bool| {
        x_scan_plan(
            &[query],
            SubjectNt::Packed(store).into(),
            1,
            None,
            11,
            mask_lowercase,
            ScoringMatrix::Blosum62,
            3,
            &lookup,
            &Default::default(),
            None,
            None,
        )
    };
    let mut lower = letters(b"ACGT", 300, 1);
    lower[100..160].make_ascii_lowercase();
    let store = PackedNucleotide::from_subject_letters(&lower).unwrap();
    assert!(scan(&store, true).is_err());
    let store = PackedNucleotide::from_subject_letters(&letters(b"ACGT", 300, 1)).unwrap();
    let unmasked = scan(&store, false).unwrap().unwrap();
    let masked = scan(&store, true).unwrap().unwrap();
    let masked = Arc::try_unwrap(masked.subject).ok().expect("one owner");
    assert_eq!(x_plan_difference(&Ok(masked), &Ok(unmasked.subject)), None);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:584-596 (the three
// bytes of four codons; see `x_phase_words`)
// The frames read from the packed bytes, by the word loop alone and with the vector loop
// first, are those of `x_generate_frames` on the store's codes: every length from 0 to 400
// bases (every phase of the last bytes, groups of 8 and 32 codons and the codons after
// them) and some long ones, for three genetic codes, and the vector loop takes the groups
// where the CPU has it.
#[test]
fn x_packed_frames_read_from_the_bytes_are_the_frames_of_the_codes() {
    use crate::sequence::packed_subject::PackedNucleotide;
    #[cfg(target_arch = "x86_64")]
    let mut vector_codons = 0usize;
    let lengths = (0..=400).chain([1_000, 4_999, 49_151, 49_153, 100_003]);
    for length in lengths {
        let subject = letters(b"ACGTACGTN", length, length as u64 + 7);
        let store = PackedNucleotide::from_subject_letters(&subject).expect("letters");
        let mut codes = vec![0u8; length];
        store.codes_into(0, &mut codes);
        for db_gencode in [1, 4, 11] {
            let tables = x_codon_tables(db_gencode).unwrap();
            let expected = x_generate_frames(&codes, tables);
            for simd in [false, true] {
                let frames = x_packed_frames_with(store.packed(), length, tables, simd);
                assert_eq!(frames.len(), expected.len(), "{length} bases, simd {simd}");
                for (frame, reference) in frames.iter().zip(&expected) {
                    assert_eq!(
                        (frame.frame, frame.aa_len, frame.orig_len, &frame.aa_seq),
                        (
                            reference.frame,
                            reference.aa_len,
                            reference.orig_len,
                            &reference.aa_seq
                        ),
                        "{length} bases, genetic code {db_gencode}, simd {simd}"
                    );
                }
            }
        }
        #[cfg(target_arch = "x86_64")]
        if x_avx2_available() && length >= 3 {
            let tables = x_codon_tables(1).unwrap();
            let codons = length / 3;
            let mut plus = vec![0u8; codons + 2];
            let mut minus = vec![0u8; codons + 2];
            // SAFETY: the CPU has AVX2.
            vector_codons += unsafe {
                x_phase_residues_avx2(store.packed(), 0, codons, tables, &mut plus, &mut minus)
            };
        }
    }
    #[cfg(target_arch = "x86_64")]
    if x_avx2_available() {
        assert!(vector_codons > 30_000, "{vector_codons}");
    }
}

// No NCBI counterpart: test of LOSAT_X_TBNXLATE128's loops (wasm32 simd128, and the x86-64
// SSSE3 / SSE4.1 build of the same code): for every phase of every length from 0 to 400 bases
// and five long ones, and three genetic codes, the residues the 128-bit loops write are those
// of the word loop (`simd` false): every codon below the count they return, and the `phase`
// codons after it when the count is not 0 (their groups start at codon `phase`), and no
// other byte.
#[test]
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn x_w128_loops_write_the_residues_of_the_word_loop() {
    use crate::sequence::packed_subject::PackedNucleotide;
    #[cfg(target_arch = "x86_64")]
    if !(std::arch::is_x86_feature_detected!("ssse3")
        && std::arch::is_x86_feature_detected!("sse4.1"))
    {
        return;
    }
    const UNTOUCHED: u8 = 0xA5;
    let mut vector_codons = 0usize;
    let lengths = (0..=400).chain([1_000, 4_999, 49_151, 49_153, 100_003]);
    for length in lengths {
        let subject = letters(b"ACGTACGTN", length, length as u64 + 11);
        let store = PackedNucleotide::from_subject_letters(&subject).expect("letters");
        let packed = store.packed();
        for db_gencode in [1, 4, 11] {
            let tables = x_codon_tables(db_gencode).unwrap();
            for phase in (0..3).filter(|&p| p + 3 <= length) {
                let codons = (length - phase) / 3;
                let mut plus = vec![0u8; codons + 2];
                let mut minus = vec![0u8; codons + 2];
                x_phase_residues(packed, phase, codons, tables, &mut plus, &mut minus, false);
                let mut w_plus = vec![UNTOUCHED; codons + 2];
                let mut w_minus = vec![UNTOUCHED; codons + 2];
                // SAFETY: on x86-64 the CPU has SSSE3 and SSE4.1 (checked above).
                let done = unsafe {
                    x_phase_residues_w128(packed, phase, codons, tables, &mut w_plus, &mut w_minus)
                };
                assert_eq!(done % 16, 0);
                vector_codons += done;
                let what = format!("{length} bases, phase {phase}, genetic code {db_gencode}");
                // the codons the groups from codon `phase` write
                let written = if done > 0 { done + phase } else { 0 };
                for j in 0..codons {
                    let (p, m) = (1 + j, codons - j);
                    if j < written {
                        assert_eq!(w_plus[p], plus[p], "plus codon {j}, {what}");
                        assert_eq!(w_minus[m], minus[m], "minus codon {j}, {what}");
                    } else {
                        assert_eq!(w_plus[p], UNTOUCHED, "plus codon {j} written, {what}");
                        assert_eq!(w_minus[m], UNTOUCHED, "minus codon {j} written, {what}");
                    }
                }
                for edge in [0, codons + 1] {
                    assert_eq!(
                        (w_plus[edge], w_minus[edge]),
                        (UNTOUCHED, UNTOUCHED),
                        "{what}"
                    );
                }
                for minus_strand in [false, true] {
                    let mut out = vec![0u8; codons + 2];
                    let mut w_out = vec![UNTOUCHED; codons + 2];
                    // SAFETY: on x86-64 the CPU has SSSE3 and SSE4.1 (checked above).
                    let one = if minus_strand {
                        x_phase_one::<true>(packed, phase, codons, tables, &mut out, false);
                        unsafe {
                            x_phase_one_w128::<true>(packed, phase, codons, tables, &mut w_out)
                        }
                    } else {
                        x_phase_one::<false>(packed, phase, codons, tables, &mut out, false);
                        unsafe {
                            x_phase_one_w128::<false>(packed, phase, codons, tables, &mut w_out)
                        }
                    };
                    assert_eq!(one, done, "{what}");
                    for j in 0..codons {
                        let at = if minus_strand { codons - j } else { 1 + j };
                        let expected = if j < written { out[at] } else { UNTOUCHED };
                        assert_eq!(
                            w_out[at], expected,
                            "codon {j}, minus {minus_strand}, {what}"
                        );
                    }
                    for edge in [0, codons + 1] {
                        assert_eq!(w_out[edge], UNTOUCHED, "{what}");
                    }
                }
            }
        }
    }
    assert!(vector_codons > 90_000, "{vector_codons}");
}

// No NCBI counterpart: test of LOSAT_X_TBNXLATE128's loops with hand-made codons: 1,024
// codons in which every codon index (0..64) stands at every one of the 16 places of a step,
// rotated by 0 to 40 codons and followed by 40 to 80 more (the end of the groups), after 0,
// 1 or 2 other bases (each phase within a byte), for every genetic code: codon `j` of the
// plus strand is `forward[index]` at `1 + j`, of the minus strand `reverse[index]` at
// `codons - j`, for every codon the loops count, and they count the 1,024.
#[test]
#[cfg(any(
    target_arch = "x86_64",
    all(target_arch = "wasm32", target_feature = "simd128")
))]
fn x_w128_loops_read_every_codon_at_every_place() {
    use crate::sequence::packed_subject::PackedNucleotide;
    #[cfg(target_arch = "x86_64")]
    if !(std::arch::is_x86_feature_detected!("ssse3")
        && std::arch::is_x86_feature_detected!("sse4.1"))
    {
        return;
    }
    // block `r` of 64 codons holds the indexes `(k + r) % 64`: place `k % 16` of a step
    // sees `L + 16 q + r` for every `q` and `r`
    let indexes: Vec<usize> = (0..16)
        .flat_map(|r| (0..64).map(move |k| (k + r) % 64))
        .collect();
    for db_gencode in genetic_codes() {
        let tables = x_codon_tables(db_gencode).unwrap();
        for phase in 0..3 {
            for more in [0, 1, 7, 15, 16, 17, 40] {
                let codon_of = |j: usize| indexes[(j + more) % indexes.len()];
                let codons = indexes.len() + 40 + more;
                let mut subject = b"TGA"[..phase].to_vec();
                for j in 0..codons {
                    let index = codon_of(j);
                    subject.extend([index >> 4, (index >> 2) & 3, index & 3].map(|b| b"ACGT"[b]));
                }
                let store = PackedNucleotide::from_subject_letters(&subject).expect("letters");
                let what = format!("phase {phase}, rotation {more}, genetic code {db_gencode}");
                let mut plus = vec![0u8; codons + 2];
                let mut minus = vec![0u8; codons + 2];
                let mut out = [vec![0u8; codons + 2], vec![0u8; codons + 2]];
                // SAFETY: on x86-64 the CPU has SSSE3 and SSE4.1 (checked above).
                let done = unsafe {
                    x_phase_residues_w128(
                        store.packed(),
                        phase,
                        codons,
                        tables,
                        &mut plus,
                        &mut minus,
                    )
                };
                // SAFETY: as above.
                let one = unsafe {
                    [
                        x_phase_one_w128::<false>(
                            store.packed(),
                            phase,
                            codons,
                            tables,
                            &mut out[0],
                        ),
                        x_phase_one_w128::<true>(
                            store.packed(),
                            phase,
                            codons,
                            tables,
                            &mut out[1],
                        ),
                    ]
                };
                assert!(done >= indexes.len(), "{done} of {codons}, {what}");
                assert_eq!(one, [done, done], "{what}");
                for j in 0..done {
                    let index = codon_of(j);
                    let (forward, reverse) = (tables.forward[index], tables.reverse[index]);
                    assert_eq!(
                        (plus[1 + j], out[0][1 + j]),
                        (forward, forward),
                        "codon {j}, {what}"
                    );
                    assert_eq!(
                        (minus[codons - j], out[1][codons - j]),
                        (reverse, reverse),
                        "codon {j}, {what}"
                    );
                }
            }
        }
    }
}

// No NCBI counterpart: test of the frame reuse; it does not change any value NCBI computes.
// The frames of a dropped plan make the next plan of the thread, shorter or longer, without
// being cleared, and that plan is still the reference plan of its letters.
#[test]
fn x_frames_of_a_dropped_plan_make_the_next_plan() {
    use crate::sequence::packed_subject::PackedNucleotide;
    let packed = |subject: &[u8], db_gencode| {
        x_packed_plan(
            &PackedNucleotide::from_subject_letters(subject).expect("letters"),
            db_gencode,
        )
    };
    let long = letters(b"ACGTACGTACGTacgtn", 5_000, 7);
    let first = packed(&long, 1).unwrap();
    let frames: Vec<*const u8> = first.frames.iter().map(|f| f.aa_seq.as_ptr()).collect();
    assert_eq!(frames.len(), 6);
    drop(first);
    let short = letters(b"ACGT", 1_001, 8);
    let next = packed(&short, 1);
    assert_eq!(
        x_plan_difference(&next, &reference_plan(&short, 1, false)),
        None
    );
    assert!(next
        .as_ref()
        .unwrap()
        .frames
        .iter()
        .all(|frame| frames.contains(&frame.aa_seq.as_ptr())));
    drop(next);
    let longer = letters(b"ACGTMRSVWYHKDBN-acgtmrsvwyhkdbn", 9_002, 9);
    assert_eq!(
        x_plan_difference(&packed(&longer, 4), &reference_plan(&longer, 4, false)),
        None
    );
    // A thread keeps one plan's frames at most.
    for _ in 0..3 {
        drop(packed(&short, 1));
    }
    assert!(X_FRAMES.with(|kept| kept.borrow().len()) <= X_REUSED_FRAMES);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:584-596 (see
// `x_phase_words`)
// One strand of one phase read by the word loop alone and with the vector loop first is the
// strand `x_phase_residues` writes: every length from 0 to 400 bases (every phase of the last
// bytes, groups of 8 and 32 codons and the codons after them) and some long ones, every phase,
// three genetic codes, written over a buffer of other bytes; the vector loop takes the groups
// where the CPU has it.
#[test]
fn x_one_strand_of_a_phase_is_that_strand_of_both() {
    use crate::sequence::packed_subject::PackedNucleotide;
    #[cfg(target_arch = "x86_64")]
    let mut vector_codons = 0usize;
    let lengths = (0..=400).chain([1_000, 4_999, 49_151, 49_153, 100_003]);
    for length in lengths {
        let subject = letters(b"ACGTACGTN", length, length as u64 + 11);
        let store = PackedNucleotide::from_subject_letters(&subject).expect("letters");
        for db_gencode in [1, 4, 11] {
            let tables = x_codon_tables(db_gencode).unwrap();
            for phase in (0..3).filter(|&phase| phase + 3 <= length) {
                let codons = (length - phase) / 3;
                let mut plus = vec![0u8; codons + 2];
                let mut minus = vec![0u8; codons + 2];
                x_phase_residues(
                    store.packed(),
                    phase,
                    codons,
                    tables,
                    &mut plus,
                    &mut minus,
                    false,
                );
                for simd in [false, true] {
                    let mut one = vec![0xAB; codons + 2];
                    x_phase_one::<false>(store.packed(), phase, codons, tables, &mut one, simd);
                    assert_eq!(one[1..=codons], plus[1..=codons], "{length} {phase} {simd}");
                    one.fill(0xCD);
                    x_phase_one::<true>(store.packed(), phase, codons, tables, &mut one, simd);
                    assert_eq!(
                        one[1..=codons],
                        minus[1..=codons],
                        "{length} {phase} {simd}"
                    );
                }
                #[cfg(target_arch = "x86_64")]
                if x_avx2_available() && db_gencode == 1 {
                    let mut one = vec![0u8; codons + 2];
                    // SAFETY: the CPU has AVX2.
                    vector_codons += unsafe {
                        x_phase_one_avx2::<true>(store.packed(), phase, codons, tables, &mut one)
                    };
                }
            }
        }
    }
    #[cfg(target_arch = "x86_64")]
    if x_avx2_available() {
        assert!(vector_codons > 90_000, "{vector_codons}");
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:805-812 (a unit reads
// `subject->sequence` of its frame from its chunk on)
// Positions `from..to` of every frame made by one strand are those of the frame
// `x_packed_frames` makes: every range of the short frames, and on the longer ones the ranges
// that start and end around the frame's ends, around a group of 32 codons and at a codon whose
// first base is a byte's last, with and without the vector loop, into one buffer whose bytes
// are those of the range before; every genetic code on one subject.
#[test]
fn x_frame_residues_of_a_range_are_the_frames() {
    use crate::sequence::packed_subject::PackedNucleotide;
    let mut buffer = Vec::new();
    let mut ranges = 0usize;
    for length in (0..=40).chain([99, 100, 101, 1_000, 4_999, 49_153]) {
        let subject = letters(b"ACGTACGTN", length, length as u64 + 13);
        let store = PackedNucleotide::from_subject_letters(&subject).expect("letters");
        let codes = if length == 4_999 {
            genetic_codes()
        } else {
            vec![1, 11]
        };
        for db_gencode in codes {
            let tables = x_codon_tables(db_gencode).unwrap();
            for frame in x_packed_frames_with(store.packed(), store.len(), tables, false) {
                let aa_len = frame.aa_len;
                let points: Vec<usize> = if aa_len <= 14 {
                    (1..=aa_len + 2).collect()
                } else {
                    let mut points: Vec<usize> = [1, 2, 3, 4, 5, 6, 31, 32, 33, 34, 35]
                        .into_iter()
                        .chain((aa_len - 12..=aa_len + 2).filter(|&p| p >= 1))
                        .chain((aa_len / 2..aa_len / 2 + 5).map(|p| p.max(1)))
                        .filter(|&p| p <= aa_len + 2)
                        .collect();
                    points.sort_unstable();
                    points.dedup();
                    points
                };
                for &from in &points {
                    for &to in points.iter().filter(|&&to| to >= from) {
                        for simd in [false, true] {
                            let (range, _, _) = x_frame_residues(
                                &store,
                                tables,
                                frame.frame,
                                aa_len,
                                from,
                                to,
                                &mut buffer,
                                simd,
                            );
                            assert_eq!(
                                buffer[range],
                                frame.aa_seq[from..to],
                                "{length} bases, code {db_gencode}, frame {}, {from}..{to}, simd {simd}",
                                frame.frame
                            );
                            ranges += 1;
                        }
                    }
                }
            }
        }
    }
    assert!(ranges > 100_000, "{ranges}");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:772-775,805-812 (the
// frames' lengths and the chunks of each frame)
// A plan without residues has the frames of `x_packed_plan` with an empty `aa_seq` and the
// same units, and the residues each unit makes are the frame's from its chunk to the margin or
// the frame's end sentinel: every length from 0 to 40 bases and some longer ones, two genetic
// codes; an unknown genetic code is refused as by `x_packed_plan`.
#[test]
fn x_lean_plans_are_the_plans_without_residues() {
    use crate::sequence::packed_subject::PackedNucleotide;
    let mut buffer = Vec::new();
    for length in (0..=40).chain([1_000, 4_999, 100_003]) {
        let subject = letters(b"ACGTMRSVWYHKDBN-acgtn", length, length as u64 + 17);
        let store = PackedNucleotide::from_subject_letters(&subject).expect("letters");
        for db_gencode in [1, 11] {
            let plan = x_packed_plan(&store, db_gencode).unwrap();
            let lean = x_packed_lean_plan(&store, db_gencode).unwrap();
            assert_eq!(lean.frames.len(), plan.frames.len(), "{length} bases");
            for (lean, frame) in lean.frames.iter().zip(&plan.frames) {
                assert!(lean.aa_seq.is_empty() && lean.aa_seq_nomask.is_none());
                assert_eq!(
                    (lean.frame, lean.aa_len, lean.orig_len, &lean.seg_masks),
                    (frame.frame, frame.aa_len, frame.orig_len, &frame.seg_masks),
                    "{length} bases"
                );
            }
            assert_eq!(lean.units, plan.units, "{length} bases");
            for (index, chunk, _) in &lean.units {
                let range = x_unit_residues(
                    &store,
                    db_gencode,
                    &lean.frames[*index],
                    *chunk,
                    &mut buffer,
                )
                .unwrap();
                let frame = &plan.frames[*index];
                let (from, to) = x_unit_range(frame.aa_len, *chunk);
                assert_eq!(to, frame.aa_len + 2, "one chunk ends at the sentinel");
                assert_eq!(buffer[range], frame.aa_seq[from..to], "{length} bases");
            }
        }
    }
    let store = PackedNucleotide::from_subject_letters(b"ACGTACGTACGTACG").unwrap();
    let refused = |plan: Result<Arc<XSubjectPlan>>| format!("{:#}", plan.err().unwrap());
    assert_eq!(
        refused(x_packed_lean_plan(&store, 7)),
        refused(x_packed_plan(&store, 7).map(Arc::new))
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:246-264 (chunks of
// MAX_DBSEQ_LEN residues overlapping by DBSEQ_CHUNK_OVERLAP)
// A subject whose frames have three chunks each (30 Mb of random bases): the plan without
// residues has the 18 units of `x_packed_plan`, and each unit's residues, made with and
// without the vector loop, are the frame's from the chunk's first residue to 4,096 past its
// last: the first chunk of each frame; the second (whose margin passes the short third) and
// the third end at the end sentinel. Its first frame is one the plan of a search leaves out
// (`x_units_translate`, with the switch).
#[test]
fn x_units_of_a_frame_of_three_chunks_read_the_frames_residues() {
    use crate::sequence::packed_subject::PackedNucleotide;
    let length: usize = 30_000_301;
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let packed: Vec<u8> = (0..length.div_ceil(4))
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 56) as u8
        })
        .collect();
    let store = PackedNucleotide::from_stored(&packed, length, &[]);
    let plan = x_packed_plan(&store, 1).unwrap();
    let lean = x_packed_lean_plan(&store, 1).unwrap();
    assert_eq!(lean.units, plan.units);
    assert_eq!(lean.units.len(), 18);
    assert!(lean.units.chunks(3).all(|units| units[0].0 == units[2].0));
    let tables = x_codon_tables(1).unwrap();
    let mut buffer = Vec::new();
    // units whose residues end 4,096 past the chunk, and at the end sentinel
    let (mut margin, mut sentinel) = (0, 0);
    for (index, chunk, _) in &lean.units {
        let frame = &plan.frames[*index];
        let (from, to) = x_unit_range(frame.aa_len, *chunk);
        if to == from + chunk.length + 4_096 {
            margin += 1;
        } else {
            assert_eq!(to, frame.aa_len + 2);
            assert!(from + chunk.length + 4_096 > to);
            sentinel += 1;
        }
        for simd in [true, false] {
            let (range, _, _) = x_frame_residues(
                &store,
                tables,
                frame.frame,
                frame.aa_len,
                from,
                to,
                &mut buffer,
                simd,
            );
            assert!(
                buffer[range] == frame.aa_seq[from..to],
                "frame {} chunk at {}, simd {simd}",
                frame.frame,
                chunk.offset
            );
        }
    }
    assert_eq!((margin, sentinel), (6, 12));
    assert_eq!(x_units_translate(length), x_units_on());
    assert!(!x_units_translate(3 * (X_MAPPED_FRAME - 2) - 1));
    assert_eq!(x_units_translate(3 * (X_MAPPED_FRAME - 2)), x_units_on());
}
