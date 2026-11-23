use super::*;

fn iupacna(letters: &[u8]) -> BinarySequence {
    seq_data_to_binary(
        &SeqData::Ncbi4na(iupacna_to_ncbi4na(letters)),
        letters.len() as u32,
        &mut Vec::new(),
    )
    .unwrap()
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

// Bytes that makeblastdb 2.17.0 wrote for these sequences (one-sequence
// and several-sequence databases of 2026-10-12; the generator is seeded
// with the sequence length, so each expectation fixes its output).
#[test]
fn ambiguous_bases_take_makeblastdb_s_random_bases() {
    // "ACGTNACGTN"
    let b = iupacna(b"ACGTNACGTN");
    assert_eq!(hex(&b.seq), "1b86c2");
    assert_eq!(hex(&b.ambig), "00000002f0000004f0000009");
    // "N" and "R", one base each
    let b = iupacna(b"N");
    assert_eq!(
        (hex(&b.seq), hex(&b.ambig)),
        ("81".into(), "00000001f0000000".into())
    );
    let b = iupacna(b"R");
    assert_eq!(
        (hex(&b.seq), hex(&b.ambig)),
        ("01".into(), "0000000150000000".into())
    );
    // "ACGTNN"
    let b = iupacna(b"ACGTNN");
    assert_eq!(
        (hex(&b.seq), hex(&b.ambig)),
        ("1b22".into(), "00000001f1000004".into())
    );
    // every IUPAC letter
    let b = iupacna(b"ACGTRYKMSWBDHVN");
    assert_eq!(hex(&b.seq), "1b994b0b");
    assert_eq!(
        hex(&b.ambig),
        "0000000b50000004a0000005c0000006300000076000000890000009e000000ad000000bb000000c7000000df000000e"
    );
}

#[test]
fn ambiguity_at_the_end_of_odd_and_even_lengths() {
    let cases: [(&[u8], &str, &str); 5] = [
        (b"ACGTN", "1b01", "00000001f0000004"),
        (b"ACGN", "1900", "00000001f0000003"),
        (b"ACN", "17", "00000001f0000002"),
        (b"AN", "22", "00000001f0000001"),
        (b"NA", "82", "00000001f0000000"),
    ];
    for (letters, seq, ambig) in cases {
        let b = iupacna(letters);
        assert_eq!(
            (hex(&b.seq), hex(&b.ambig)),
            (seq.into(), ambig.into()),
            "{letters:?}"
        );
    }
}

#[test]
fn a_text_sequence_without_ambiguity_gets_an_empty_ambiguity_block() {
    // makeblastdb keeps "G" and "GA" as text Seq-data: 4 bytes of count 0.
    let b = iupacna(b"G");
    assert_eq!(
        (hex(&b.seq), hex(&b.ambig)),
        ("81".into(), "00000000".into())
    );
    let b = iupacna(b"GA");
    assert_eq!(
        (hex(&b.seq), hex(&b.ambig)),
        ("82".into(), "00000000".into())
    );
}

#[test]
fn ncbi2na_has_no_ambiguity_data() {
    // "GAT", "GATC", "GATCA", "GATCGATCG" as ncbi2na
    let cases: [(&[u8], u32, &str); 4] = [
        (&[0x8c], 3, "8f"),
        (&[0x8d], 4, "8d00"),
        (&[0x8d, 0x00], 5, "8d01"),
        (&[0x8d, 0x8d, 0x80], 9, "8d8d81"),
    ];
    for (packed, length, seq) in cases {
        let b = seq_data_to_binary(&SeqData::Ncbi2na(packed.to_vec()), length, &mut Vec::new())
            .unwrap();
        assert_eq!(hex(&b.seq), seq);
        assert!(b.ambig.is_empty());
        assert_eq!(find_sequence_length(false, &b.seq), length);
    }
}

#[test]
fn runs_use_the_old_format_up_to_15_and_split_at_4095() {
    let run = |n: usize| {
        let mut letters = vec![b'A'];
        letters.extend(std::iter::repeat_n(b'N', n));
        letters.push(b'C');
        iupacna(&letters).ambig
    };
    // 15 N: old format, one word: value f, length - 1 = e, offset 1.
    assert_eq!(hex(&run(15)), "00000001fe000001");
    // 16 N: new format, two words.
    assert_eq!(hex(&run(16)), "80000002f00f000000000001");
    // 4095 N: one region; 4096 N: two regions (4095 + 1).
    assert_eq!(hex(&run(4095)), "80000002fffe000000000001");
    assert_eq!(hex(&run(4096)), "80000004fffe000000000001f000000000001000");
}

#[test]
fn protein_letters_become_ncbistdaa() {
    let b = seq_data_to_binary(&SeqData::Ncbieaa(b"MK*V-L*".to_vec()), 7, &mut Vec::new()).unwrap();
    assert_eq!(hex(&b.seq), "0c0a1913000b19");
    let b = seq_data_to_binary(
        &SeqData::Ncbieaa(b"ABCDEFGHIJKLMNOPQRSTUVWXYZ".to_vec()),
        26,
        &mut Vec::new(),
    )
    .unwrap();
    assert_eq!(
        hex(&b.seq),
        "0102030405060708091b0a0b0c0d1a0e0f101112181314151617"
    );
    assert!(b.ambig.is_empty());
}

#[test]
fn other_seq_data_is_not_processable() {
    assert_eq!(
        seq_data_to_binary(&SeqData::Other, 10, &mut Vec::new()),
        Err(UnprocessableSequence)
    );
}
/// `SeqDB_SequenceHash`: the keys makeblastdb 2.17.0 wrote into `.phd` / `.nhd` for
/// makedb-v5 survey B fixtures h01_prot and h02_nucl (`-hash_index`).
#[test]
fn sequence_hashes_are_makeblastdb_s() {
    for (letters, expected) in [
        (&b"ACDEFGHIKLMNPQRSTVWY"[..], 55_887_582u32),
        (b"XBZUOJ*ACDX", 329_294_776),
        (b"MKV-LLA--W", 1_592_154_558),
        (b"M", 12_357),
    ] {
        let data = SeqData::Ncbieaa(letters.to_vec());
        assert_eq!(sequence_hash(&data, letters.len() as u64), Some(expected));
        // The hash runs over the ncbistdaa codes of the sequence file.
        let stdaa = seq_data_to_binary(&data, letters.len() as u32, &mut Vec::new())
            .unwrap()
            .seq;
        let over_codes = stdaa.iter().fold(0u32, |h, &code| {
            h.wrapping_mul(1_103_515_245)
                .wrapping_add(u32::from(code).wrapping_add(12345))
        });
        assert_eq!(over_codes, expected);
    }
    // Nucleotides: NCBI4na per base, a gap read as N (0x0f, `s_FixBioseqDeltas`).
    for (letters, expected) in [
        (&b"ACGTACGTAC"[..], 4_241_009_515u32),
        (b"ACGTRYKMSWBDHVN", 3_170_882_851),
        (b"ACGTNACGTN", 2_994_364_842),
        (b"ACGTNNNNNNNNNNACGT", 3_309_723_638),
    ] {
        let na4 = SeqData::Ncbi4na(iupacna_to_ncbi4na(letters));
        assert_eq!(sequence_hash(&na4, letters.len() as u64), Some(expected));
    }
    let gaps = SeqData::Ncbi4na(vec![0x12, 0x48, 0xff, 0x12, 0x48, 0xf1]);
    assert_eq!(sequence_hash(&gaps, 12), Some(2_633_391_912));
    // ncbi2na: ACGTACGTAC as 00 01 10 11 ...
    let na2 = SeqData::Ncbi2na(vec![0x1b, 0x1b, 0x10]);
    assert_eq!(sequence_hash(&na2, 10), Some(4_241_009_515));
    assert_eq!(sequence_hash(&SeqData::Other, 3), None);
}

/// `x_ComputeHash` of a sequence in the form of the sequence file: a nucleotide of no base hashes
/// nothing before its ambiguity data are read (even malformed ones), as `SeqDB_UnpackAmbiguities`
/// returns there; one with bases and ambiguity data shorter than their count is the error that
/// the writer refuses.
#[test]
fn a_raw_nucleotide_of_no_base_hashes_before_its_ambiguities() {
    let malformed = 5u32.to_be_bytes();
    assert_eq!(raw_sequence_hash(false, &[0x00], &malformed), Ok(0));
    assert_eq!(raw_sequence_hash(false, &[0x00], &[]), Ok(0));
    assert_eq!(
        raw_sequence_hash(false, &[0x1b, 0x02], &malformed),
        Err("the ambiguity data are shorter than their count")
    );
    // ACGTAC: two bases in the last byte.
    let na2 = SeqData::Ncbi2na(vec![0x1b, 0x10]);
    assert_eq!(
        raw_sequence_hash(false, &[0x1b, 0x12], &[]).ok(),
        sequence_hash(&na2, 6)
    );
}

fn unhex(text: &str) -> Vec<u8> {
    (0..text.len())
        .step_by(2)
        .map(|at| u8::from_str_radix(&text[at..at + 2], 16).unwrap())
        .collect()
}

/// The sequence bytes and the standard output of `x_CookSequence` of `data` and `length`.
fn cooked(data: SeqData, length: u32) -> (String, String) {
    let mut cout = Vec::new();
    let binary = seq_data_to_binary(&data, length, &mut cout).unwrap();
    (hex(&binary.seq), String::from_utf8(cout).unwrap())
}

// NCBI makeprofiledb 2.17.0 on profiles whose query is in a nucleotide coding (row PS10d survey
// `PR/ps10d/probe`, 2026-11-23, two runs equal, valgrind clean): the sequence bytes that it
// wrote into the protein volume's `.psq` (the bytes between the sequence's offsets in the
// `.pin`, `probe/dbinfo.txt`) and its standard output, for the cases of the survey's §4 (`SA_*`, `VC_*`, `L*_*`). The query's bases
// are `ATGGCGTACCTGAAGTTCGACA` (22) and their first 21 or fewer.
#[test]
fn nucleotide_codings_give_ncbis_sequence_bytes_and_lines() {
    let iupacna = |letters: &str| SeqData::Iupacna(letters.as_bytes().to_vec());
    let na2 = |data: &str| SeqData::Ncbi2na(unhex(data));
    let na4 = |data: &str| SeqData::Ncbi4na(unhex(data));
    let c1 = "ATGGCGTACCTGAAGTTCGACA";
    let cases: Vec<(&str, SeqData, u32, &str, &str)> = vec![
        ("SA_na", iupacna(c1), 22, "3a6c5e0bd842", ""),
        ("SA_n2", na2("3A6C5E0BD840"), 22, "3a6c5e0bd842", ""),
        (
            "SA_n4",
            na4("1844248122841148824121"),
            22,
            "3a6c5e0bd842",
            "",
        ),
        // ncbi2na: bytes past the length's are not read; the last byte's unused bits stay.
        (
            "VC_n2_more2",
            na2("3A6C5E0BD8400000"),
            22,
            "3a6c5e0bd842",
            "",
        ),
        ("VC_n2_pad", na2("3A6C5E0BD84F"), 22, "3a6c5e0bd84e", ""),
        ("L21_n2_pad", na2("3A6C5E0BD87F"), 21, "3a6c5e0bd87d", ""),
        // ncbi4na: another byte count is a line; ambiguity codes and gaps are random bases (NCBI's
        // generator seeded with the length); a pad nibble that is a base beside a base leaks its
        // 2 bits.
        (
            "VC_n4_more1",
            na4("184424812284114882412111"),
            22,
            "3a6c5e0bd842",
            "ib=11,n4sz=12\n",
        ),
        (
            "VC_n4_more2",
            na4("18442481228411488241211111"),
            22,
            "3a6c5e0bd842",
            "ib=11,n4sz=13\n",
        ),
        (
            "VC_n4_N",
            na4("184F24812FF4114882412F"),
            22,
            "386c6e0bd842",
            "",
        ),
        (
            "VC_n4_RY",
            na4("184524A12284114552412A"),
            22,
            "38645e0a1852",
            "",
        ),
        (
            "VC_n4_gap",
            na4("0844240122841140024120"),
            22,
            "3a685e0b1842",
            "",
        ),
        (
            "VC_n4_allN",
            na4("FFFFFFFFFFFFFFFFFFFFFF"),
            22,
            "2c00c584c512",
            "",
        ),
        (
            "VC_n4_codes",
            na4("0123456789ABCDEF102030"),
            22,
            "04a4c5cd2442",
            "",
        ),
        (
            "L21_n4",
            na4("1844248122841148824120"),
            21,
            "3a6c5e0bd841",
            "",
        ),
        (
            "L21_n4_padA",
            na4("1844248122841148824121"),
            21,
            "3a6c5e0bd841",
            "",
        ),
        (
            "L21_n4_padN",
            na4("184424812284114882412F"),
            21,
            "3a6c5e0bd841",
            "",
        ),
        (
            "L21_n4_padT",
            na4("1844248122841148824128"),
            21,
            "3a6c5e0bd871",
            "",
        ),
        // iupacna: every letter converts (either case, U as T, any other byte N); the line
        // compares the NCBI4na bytes of the letters with the length's; one byte short reads the
        // string's terminating zero (N N).
        (
            "VC_na_lower",
            iupacna("atggcgtacctgaagttcgaca"),
            22,
            "3a6c5e0bd842",
            "",
        ),
        (
            "VC_na_U",
            iupacna("AUGGCGUACCUGAAGuuCGACA"),
            22,
            "3a6c5e0bd842",
            "",
        ),
        (
            "VC_na_amb",
            iupacna("ANRYKMSWBDHVNACGTACGTA"),
            22,
            "0384a246c6c2",
            "",
        ),
        (
            "VC_na_bad",
            iupacna("A*C0G9T-AXCZG.T AE#IQJ"),
            22,
            "06bc048c3162",
            "",
        ),
        (
            "VC_na_more1",
            iupacna("ATGGCGTACCTGAAGTTCGACAG"),
            22,
            "3a6c5e0bd842",
            "ib=11,n4sz=12\n",
        ),
        (
            "L21_na_more1T",
            iupacna("ATGGCGTACCTGAAGTTCGACT"),
            21,
            "3a6c5e0bd871",
            "",
        ),
        ("L21_na_more1A", iupacna(c1), 21, "3a6c5e0bd841", ""),
        ("VC_na_less1", iupacna(&c1[..21]), 22, "3a6c5e0bd842", ""),
        (
            "VC_na_less2",
            iupacna(&c1[..20]),
            22,
            "3a6c5e0bd822",
            "ib=11,n4sz=10\n",
        ),
        (
            "VC_na_less3",
            iupacna(&c1[..19]),
            22,
            "3a6c5e0bd8b2",
            "ib=11,n4sz=10\n",
        ),
        ("L1_naE", iupacna(""), 1, "81", "ib=1,n4sz=0\n"),
        ("L2_naE", iupacna(""), 2, "82", "ib=1,n4sz=0\n"),
        ("L1_na", iupacna("A"), 1, "01", ""),
        ("L2_na", iupacna("AT"), 2, "32", ""),
        ("L3_na", iupacna("ATG"), 3, "3b", ""),
        ("L1_n2", na2("00"), 1, "01", ""),
        ("L1_n4", na4("10"), 1, "01", ""),
        ("L4_n2", na2("3A"), 4, "3a00", ""),
        ("L5_n4", na4("184420"), 5, "3a41", ""),
    ];
    for (case, data, length, seq, line) in cases {
        assert_eq!(
            cooked(data, length),
            (seq.to_string(), line.to_string()),
            "{case}"
        );
    }
}

// NCBI makeprofiledb 2.17.0 reads past the data from these sizes on (row PS10d survey §6: valgrind,
// SIGSEGV for no bytes, bytes that change between runs); each boundary of the three codings.
#[test]
fn short_nucleotide_data_is_found_at_each_boundary() {
    let short = |data: SeqData, length: u32| {
        short_nucleotide_data(&data, length).map(|short| (short.coding, short.held, short.read))
    };
    let bytes = |n: usize| vec![0x1b; n];
    // ncbi2na: ceil(length / 4) bytes.
    for (length, need) in [(1, 1), (4, 1), (5, 2), (21, 6), (22, 6), (24, 6), (66, 17)] {
        assert_eq!(short(SeqData::Ncbi2na(bytes(need)), length), None);
        assert_eq!(short(SeqData::Ncbi2na(bytes(need + 1)), length), None);
        assert_eq!(
            short(SeqData::Ncbi2na(bytes(need - 1)), length),
            Some(("ncbi2na", need - 1, need)),
            "ncbi2na {length}"
        );
    }
    assert_eq!(
        short(SeqData::Ncbi2na(Vec::new()), 22),
        Some(("ncbi2na", 0, 6))
    );
    // ncbi4na: ceil(length / 2) bytes.
    for (length, need) in [(1, 1), (2, 1), (3, 2), (21, 11), (22, 11), (66, 33)] {
        assert_eq!(short(SeqData::Ncbi4na(bytes(need)), length), None);
        assert_eq!(
            short(SeqData::Ncbi4na(bytes(need - 1)), length),
            Some(("ncbi4na", need - 1, need)),
            "ncbi4na {length}"
        );
    }
    // iupacna: the NCBI4na bytes of its letters and the terminating zero.
    let letters = |n: usize| SeqData::Iupacna(vec![b'A'; n]);
    for (length, fewest) in [(1, 0), (2, 0), (3, 1), (4, 1), (21, 19), (22, 19), (66, 63)] {
        assert_eq!(short(letters(fewest), length), None, "iupacna {length}");
        if fewest > 0 {
            let read = length.div_ceil(2) as usize;
            assert_eq!(
                short(letters(fewest - 1), length),
                Some(("iupacna", fewest - 1, read)),
                "iupacna {length}"
            );
        }
    }
    assert_eq!(short(letters(0), 3), Some(("iupacna", 0, 2)));
    assert_eq!(short(letters(0), 22), Some(("iupacna", 0, 11)));
    // The protein codings and the others are not converted so.
    assert_eq!(short(SeqData::Ncbieaa(Vec::new()), 22), None);
    assert_eq!(short(SeqData::Other, 22), None);
}
