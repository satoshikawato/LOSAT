use super::*;

fn iupacna(letters: &[u8]) -> BinarySequence {
    seq_data_to_binary(
        &SeqData::Ncbi4na(iupacna_to_ncbi4na(letters)),
        letters.len() as u32,
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
        let b = seq_data_to_binary(&SeqData::Ncbi2na(packed.to_vec()), length).unwrap();
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
    let b = seq_data_to_binary(&SeqData::Ncbieaa(b"MK*V-L*".to_vec()), 7).unwrap();
    assert_eq!(hex(&b.seq), "0c0a1913000b19");
    let b = seq_data_to_binary(
        &SeqData::Ncbieaa(b"ABCDEFGHIJKLMNOPQRSTUVWXYZ".to_vec()),
        26,
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
        seq_data_to_binary(&SeqData::Other, 10),
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
        let stdaa = seq_data_to_binary(&data, letters.len() as u32).unwrap().seq;
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
