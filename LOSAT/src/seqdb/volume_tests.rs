use super::{overlay_ambiguity_runs, AmbiguityRun, NucleotideAmbiguity, SequenceBlock};
use crate::seqdb::SeqType;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:896-945
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:403-442
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2350-2385
/// `GetBioseq`'s data of three nucleotides of 5 bases whose last stored byte holds bits after
/// the last base (as NCBI's distributed databases have): without a word of ambiguity data (none,
/// or 3 bytes) `ncbi2na` keeps those bits and clears the count; with ambiguity data (a run of N
/// over the last base) `ncbi4na` clears them first, so the half byte after the last base is A
/// (1).
#[test]
fn seq_data_keeps_the_stored_bits_of_ncbi2na_and_clears_them_for_ncbi4na() {
    let block = SequenceBlock {
        seqtype: SeqType::Nucleotide,
        first_ordinal: 0,
        first_oid: 0,
        path: std::path::PathBuf::from("test.nsq"),
        data: vec![
            0x1B, 0xBD, // ACGT, G and the bits 1111, count 1
            0x1B, 0xBD, 0, 0, 0, 1, 0xF0, 0, 0, 4, // the same, one run: N at base 4
            0x1B, 0xBD, 0x80, 0, 0, // the same, 3 bytes of ambiguity data: no word
        ],
        data_start: 0,
        seq_offsets: vec![0, 2, 12, 17],
        amb_offsets: vec![2, 4, 14],
        stopped_at_negative: false,
        amb_error: None,
    };
    assert_eq!(
        block.nucleotide_seq_data(0).unwrap(),
        (true, vec![0x1B, 0xBC])
    );
    assert_eq!(
        block.nucleotide_seq_data(2).unwrap(),
        (true, vec![0x1B, 0xBC])
    );
    assert_eq!(
        block.nucleotide_seq_data(1).unwrap(),
        (false, vec![0x12, 0x48, 0xF1])
    );
    assert_eq!(block.nucleotide_iupac(1).unwrap(), b"ACGTN");
}

/// A small deterministic generator (xorshift), so the test needs no fixture and no crate.
struct Rng(u64);
impl Rng {
    fn next(&mut self, below: u64) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0 % below
    }
}

/// The letters of the whole sequence with every run applied in stored order (the decode
/// that `SequenceBlock` does), the reference for any range.
fn whole_decode(length: usize, runs: &[AmbiguityRun]) -> Vec<u8> {
    let mut letters = vec![b'A'; length];
    overlay_ambiguity_runs(runs, false, 0, length, &mut letters);
    letters
}

fn check(length: usize, runs: Vec<AmbiguityRun>, rng: &mut Rng, ranges: usize) -> bool {
    let whole = whole_decode(length, &runs);
    let ambiguity = NucleotideAmbiguity::new(length, 0, runs);
    for _ in 0..ranges {
        let from = rng.next(length as u64 + 1) as usize;
        let to = from + rng.next((length - from) as u64 + 1) as usize;
        let mut letters = vec![b'A'; to - from];
        ambiguity.overlay(from, to, &mut letters);
        assert_eq!(letters, whole[from..to], "range {from}..{to}");
    }
    ambiguity.ordered
}

/// The ambiguity loop of `nucleotide_blastna_in` over the whole sequence.
fn whole_runs(runs: &[AmbiguityRun], seq: &mut [u8]) {
    let region_end = seq.len() as i64;
    for run in runs {
        let mut position = i64::from(run.offset as i32);
        let row_len = i64::from(run.length);
        if position + row_len <= 0 {
            continue;
        }
        if position >= region_end {
            break;
        }
        for _ in 0..row_len {
            if position >= 0 && position < region_end {
                seq[position as usize] = run.ncbi4na;
            }
            position += 1;
        }
    }
}

// LOSAT_X_BNTBRANGE: the runs over ranges give the whole decode's bytes inside the ranges,
// for ordered runs, shuffled and overlapping ones, offsets of 2^31 and more (negative as
// NCBI's Int4), and runs past the end (which stop the whole decode's loop).
#[test]
fn runs_over_ranges_equal_the_whole_decode_inside_them() {
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);
    for case in 0..400 {
        let length = 1 + rng.next(3000) as usize;
        let mut runs = Vec::new();
        for _ in 0..rng.next(40) {
            let offset = match rng.next(20) {
                0 => 0x8000_0000u32 + rng.next(1 << 20) as u32,
                1 => length as u32 + rng.next(100) as u32,
                2 => u32::MAX - rng.next(50) as u32,
                _ => rng.next(length as u64) as u32,
            };
            let longest = if rng.next(10) == 0 { 5000 } else { 20 };
            runs.push(AmbiguityRun {
                offset,
                length: 1 + rng.next(longest) as u32,
                ncbi4na: rng.next(16) as u8,
            });
        }
        if case % 2 == 0 {
            runs.sort_by_key(|run| run.offset);
        }
        let mut ranges = Vec::new();
        let mut at = 0usize;
        while at < length {
            let from = at + rng.next(200) as usize;
            let to = (from + rng.next(300) as usize).min(length);
            if from < to {
                ranges.push((from, to));
            }
            at = to + 1 + rng.next(100) as usize;
        }
        let mut whole = vec![1u8; length];
        whole_runs(&runs, &mut whole);
        let mut ranged = vec![1u8; length];
        super::x_overlay_runs_in_ranges(&runs, &ranges, &mut ranged, &super::X_NCBI4NA);
        for &(from, to) in &ranges {
            assert_eq!(
                ranged[from..to],
                whole[from..to],
                "case {case}, range {from}..{to}"
            );
        }
    }
}

#[test]
fn many_runs_and_many_ranges_equal_the_whole_decode() {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    // 3000 ordered runs of 1..=20 bases with gaps of 0..=30 bases (NCBI's layout).
    let mut runs = Vec::new();
    let mut at = 5u32;
    for _ in 0..3000 {
        let length = 1 + rng.next(20) as u32;
        runs.push(AmbiguityRun {
            offset: at,
            length,
            ncbi4na: 2 + rng.next(14) as u8,
        });
        at += length + rng.next(31) as u32;
    }
    let length = at as usize + 10;
    assert!(check(length, runs.clone(), &mut rng, 5000));
    // Ranges that start inside a run, end inside one, or lie between two.
    for run in runs.iter().step_by(97) {
        let mut one = [b'A'; 2];
        let ambiguity = NucleotideAmbiguity::new(length, 0, runs.clone());
        ambiguity.overlay(run.offset as usize, run.offset as usize + 1, &mut one[..1]);
        assert_ne!(one[0], b'A');
    }
    // Stored runs that overlap or are out of order: every run applies in stored order, so
    // a later run still overwrites an earlier one.
    let mut messy = runs.clone();
    messy.swap(10, 2900);
    messy.push(AmbiguityRun {
        offset: 7,
        length: 5000,
        ncbi4na: 3,
    });
    messy.push(AmbiguityRun {
        offset: 100,
        length: 3,
        ncbi4na: 9,
    });
    assert!(!check(length, messy, &mut rng, 5000));
}

#[test]
fn a_range_outside_every_run_stays_unchanged() {
    let runs = vec![AmbiguityRun {
        offset: 10,
        length: 4,
        ncbi4na: 15,
    }];
    let ambiguity = NucleotideAmbiguity::new(40, 0, runs);
    for (from, to) in [(0, 10), (14, 40), (0, 0), (40, 40)] {
        let mut letters = vec![b'C'; to - from];
        ambiguity.overlay(from, to, &mut letters);
        assert!(letters.iter().all(|&letter| letter == b'C'), "{from}..{to}");
    }
    let mut letters = vec![b'C'; 6];
    ambiguity.overlay(8, 14, &mut letters);
    assert_eq!(letters, b"CCNNNN");
}

// LOSAT_X_BNDECODE: the by-byte decode gives the bytes of the decode `=0` runs, for the whole
// sequence and for ranges that start and end inside a stored byte, with every stored byte
// value, ordered and shuffled runs, runs of each code, offsets of 2^31 and more (negative as
// NCBI's Int4) and runs past the end.
#[test]
fn the_by_byte_decode_equals_the_by_base_decode() {
    let mut rng = Rng(0x0DDB_1A5E_5BAD_5EED);
    for case in 0..600 {
        let mut length = 1 + rng.next(if case % 3 == 0 { 9 } else { 4000 }) as usize;
        let mut packed: Vec<u8> = (0..length.div_ceil(4))
            .map(|_| rng.next(256) as u8)
            .collect();
        if case == 0 {
            // Every stored byte value once.
            packed = (0..=255u8).collect();
            length = packed.len() * 4;
        }
        let mut runs = Vec::new();
        for _ in 0..rng.next(30) {
            let offset = match rng.next(20) {
                0 => 0x8000_0000u32 + rng.next(1 << 20) as u32,
                1 => length as u32 + rng.next(100) as u32,
                2 => u32::MAX - rng.next(50) as u32,
                _ => rng.next(length as u64) as u32,
            };
            let longest = if rng.next(10) == 0 { 5000 } else { 20 };
            runs.push(AmbiguityRun {
                offset,
                length: 1 + rng.next(longest) as u32,
                ncbi4na: rng.next(16) as u8,
            });
        }
        if case % 2 == 0 {
            runs.sort_by_key(|run| run.offset);
        }
        let mut by_base = vec![0u8; length];
        super::x_blastna_by_base(&packed, &runs, &mut by_base);
        let mut by_byte = vec![0u8; length];
        super::x_blastna_by_byte(&packed, &runs, &[(0, length)], &mut by_byte);
        assert_eq!(by_byte, by_base, "case {case}, whole");
        let mut ranges = Vec::new();
        let mut at = rng.next(7) as usize;
        while at < length {
            let to = (at + 1 + rng.next(if case % 2 == 0 { 9 } else { 300 }) as usize).min(length);
            ranges.push((at, to));
            at = to + rng.next(12) as usize;
        }
        let mut by_base = vec![201u8; length];
        super::x_blastna_by_base_ranges(&packed, &runs, &ranges, &mut by_base);
        let mut by_byte = vec![201u8; length];
        super::x_blastna_by_byte(&packed, &runs, &ranges, &mut by_byte);
        assert_eq!(by_byte, by_base, "case {case}, ranges {ranges:?}");
    }
}

#[test]
fn ncbi_protein_length_wraps_signed_32_bit_boundaries() {
    use super::ncbi_protein_length;
    assert_eq!(ncbi_protein_length(0), 0);
    assert_eq!(ncbi_protein_length(1), 1);
    assert_eq!(ncbi_protein_length(i32::MAX as u32), i32::MAX);
    assert_eq!(ncbi_protein_length(1u32 << 31), i32::MIN);
    assert_eq!(ncbi_protein_length((1u32 << 31) + 10), -2_147_483_638);
    assert_eq!(ncbi_protein_length(u32::MAX), -1);
}

fn push_be_u32s(out: &mut Vec<u8>, vals: &[u32]) {
    for &v in vals {
        out.extend_from_slice(&v.to_be_bytes());
    }
}

fn assert_empty_sentinels(blk: &super::SequenceBlock, idx: usize) {
    assert_eq!(blk.protein_with_sentinels(idx), &[0, 0]);
}

/// Writes a minimal version-4 BLAST index file (`.pin` or `.nin`) with the given metadata and
/// offset arrays for unit testing length handling without allocating gigabyte buffers.
fn write_v4_index(
    path: &std::path::Path,
    seqtype: SeqType,
    volume_length: u64,
    max_length: u32,
    seq_offsets: &[u32],
    amb_offsets: Option<&[u32]>,
) {
    let num_oids = (seq_offsets.len() - 1) as u32;
    let mut out = Vec::new();
    out.extend_from_slice(&4u32.to_be_bytes());
    let code: u32 = match seqtype {
        SeqType::Protein => 1,
        SeqType::Nucleotide => 0,
    };
    out.extend_from_slice(&code.to_be_bytes());
    let title = b"unit";
    out.extend_from_slice(&(title.len() as u32).to_be_bytes());
    out.extend_from_slice(title);
    let date = b"Jan 1, 2026  12:00 PM";
    out.extend_from_slice(&(date.len() as u32).to_be_bytes());
    out.extend_from_slice(date);
    out.extend_from_slice(&num_oids.to_be_bytes());
    out.extend_from_slice(&volume_length.to_le_bytes());
    out.extend_from_slice(&max_length.to_be_bytes());
    for _ in 0..=num_oids {
        out.extend_from_slice(&0u32.to_be_bytes());
    }
    push_be_u32s(&mut out, seq_offsets);
    if let Some(amb) = amb_offsets {
        push_be_u32s(&mut out, amb);
    }
    std::fs::write(path, out).unwrap();
}

/// Writes a sparse `.psq` file of `total_len` bytes with sentinel `0` at offset 0 and any
/// `(offset, bytes)` segments placed inside it; the tail up to `total_len` is extended via
/// `set_len` so no multi-gigabyte buffer is allocated or written.
fn write_sparse_psq(path: &std::path::Path, total_len: u64, segments: &[(u64, &[u8])]) {
    use std::io::{Seek, SeekFrom, Write};
    let mut f = std::fs::File::create(path).unwrap();
    f.write_all(&[0u8]).unwrap();
    for &(offset, bytes) in segments {
        f.seek(SeekFrom::Start(offset)).unwrap();
        f.write_all(bytes).unwrap();
    }
    f.set_len(total_len).unwrap();
}

#[cfg(unix)]
fn write_oidlist_file(path: &std::path::Path, num_oids: u32) {
    let mut bytes = num_oids.to_be_bytes().to_vec();
    let words = (num_oids as usize).div_ceil(32);
    for w in 0..words {
        let mut word = 0u32;
        for bit in 0..32 {
            let oid = (w * 32 + bit) as u32;
            if oid < num_oids {
                word |= 1u32 << (31 - bit);
            }
        }
        bytes.extend_from_slice(&word.to_be_bytes());
    }
    std::fs::write(path, bytes).unwrap();
}

#[test]
fn sparse_protein_volumes_match_ncbi_32bit_lengths_and_keep_nucleotide_widened() {
    use crate::seqdb::{BlockRange, SearchPath, SeqDb};
    use crate::utils::test_dir::TempDir;

    const TWO_POW_31: u32 = 1u32 << 31;
    let dir = TempDir::unique("prot32-unit");
    let root: &std::path::Path = &dir;

    // 2. Single 0-residue protein volume (`db0`) and normal 40-residue volume (`dbs`).
    write_v4_index(&root.join("db0.pin"), SeqType::Protein, 0, 0, &[1, 2], None);
    write_sparse_psq(&root.join("db0.psq"), 2, &[]);
    std::fs::write(root.join("db0.phr"), b"").unwrap();

    let residues = [1u8; 40];
    write_v4_index(
        &root.join("dbs.pin"),
        SeqType::Protein,
        40,
        40,
        &[1, 42],
        None,
    );
    write_sparse_psq(&root.join("dbs.psq"), 42, &[(1, &residues)]);
    std::fs::write(root.join("dbs.phr"), b"").unwrap();

    // 3. Single volume with an empty protein between two normal proteins (`mix0`: [40, 0, 40]):
    // must use the single-read fast path and reuse the caller's `kept` buffer.
    write_v4_index(
        &root.join("mix0.pin"),
        SeqType::Protein,
        80,
        40,
        &[1, 42, 43, 84],
        None,
    );
    write_sparse_psq(
        &root.join("mix0.psq"),
        84,
        &[(1, &residues), (43, &residues)],
    );
    std::fs::write(root.join("mix0.phr"), b"").unwrap();

    // 4. Nucleotide volume whose `.nin` max_length is 2^31: must stay widened, NOT wrapped to 0.
    write_v4_index(
        &root.join("nucl.nin"),
        SeqType::Nucleotide,
        u64::from(TWO_POW_31),
        TWO_POW_31,
        &[0, 4],
        Some(&[4, 4]),
    );
    std::fs::write(root.join("nucl.nsq"), [0u8; 4]).unwrap();
    std::fs::write(root.join("nucl.nhr"), b"").unwrap();

    let path = SearchPath::new(root, None, None);

    // Check db0
    let db = SeqDb::open_search_in(b"db0", SeqType::Protein, &path).unwrap();
    assert_eq!(db.max_length(), 0);
    assert_eq!(db.protein_length_i32(0).unwrap(), 0);
    db.check_protein_range(0).unwrap();
    let blk = db.read_protein_search_sequence(0).unwrap();
    assert_empty_sentinels(&blk, 0);

    // Check mix0 ([40, 0, 40]) reuses `kept` on the single-read fast path
    let db = SeqDb::open_search_in(b"mix0", SeqType::Protein, &path).unwrap();
    assert_eq!(db.x_protein_lengths_i32(0, 3).unwrap(), vec![40, 0, 40]);
    let mut kept = Vec::with_capacity(256);
    let kept_ptr = kept.as_ptr();
    let blk = db
        .x_read_protein_search_block_into(
            BlockRange {
                volume: 0,
                first: 0,
                end: 3,
                first_ordinal: 0,
            },
            &mut kept,
        )
        .unwrap();
    assert_eq!(blk.protein_bytes().as_ptr(), kept_ptr);
    assert_eq!(blk.protein(0), &residues);
    assert_empty_sentinels(&blk, 1);
    assert_eq!(blk.protein(2), &residues);

    // Check nucleotide max_length >= 2^31 is NOT wrapped
    let ndb = SeqDb::open_search_in(b"nucl", SeqType::Nucleotide, &path).unwrap();
    assert_eq!(ndb.max_length(), TWO_POW_31);

    // 5. Sparse >= 2^31-residue protein volumes (Unix sparse holes; Windows CLI coverage lives in
    // `LOSAT/tests/protein_db_32bit_lengths.rs` via `fsutil sparse setflag`).
    #[cfg(unix)]
    {
        let dbh = root.join("dbh");
        write_v4_index(
            &root.join("dbh.pin"),
            SeqType::Protein,
            u64::from(TWO_POW_31),
            TWO_POW_31,
            &[1, TWO_POW_31 + 2],
            None,
        );
        write_sparse_psq(&root.join("dbh.psq"), u64::from(TWO_POW_31) + 2, &[]);
        std::fs::write(root.join("dbh.phr"), b"").unwrap();

        // Single volume mixing a 40-residue protein followed by a 2^31-residue protein (`mix2`).
        let m2_s1 = 42u32;
        let m2_s2 = m2_s1 + TWO_POW_31 + 1;
        write_v4_index(
            &root.join("mix2.pin"),
            SeqType::Protein,
            u64::from(TWO_POW_31) + 40,
            TWO_POW_31,
            &[1, m2_s1, m2_s2],
            None,
        );
        write_sparse_psq(&root.join("mix2.psq"), u64::from(m2_s2), &[(1, &residues)]);
        std::fs::write(root.join("mix2.phr"), b"").unwrap();

        // 2^31 + 10 volume (`dbg`) for `read_protein_slice_at` tail reads without allocating 2 GiB.
        let dbg_len = TWO_POW_31 + 10;
        let tail20: [u8; 20] = [
            1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22,
        ];
        write_v4_index(
            &root.join("dbg.pin"),
            SeqType::Protein,
            u64::from(dbg_len),
            dbg_len,
            &[1, dbg_len + 2],
            None,
        );
        write_sparse_psq(
            &root.join("dbg.psq"),
            u64::from(dbg_len) + 2,
            &[(u64::from(dbg_len) - 19, &tail20)],
        );
        std::fs::write(root.join("dbg.phr"), b"").unwrap();

        // Aliases: dba (dbh + dbs), dbh_len0 (LENGTH 0), dba_ol (OIDLIST triggering x_ScanTotals).
        std::fs::write(root.join("dba.pal"), "TITLE dba\nDBLIST dbh dbs\n").unwrap();
        std::fs::write(
            root.join("dbh_len0.pal"),
            "TITLE dbh_len0\nDBLIST dbh\nLENGTH 0\n",
        )
        .unwrap();
        write_oidlist_file(&root.join("ol2.msk"), 2);
        std::fs::write(
            root.join("dba_ol.pal"),
            "TITLE dba_ol\nDBLIST dbh dbs\nOIDLIST ol2.msk\n",
        )
        .unwrap();

        // Check dbh
        let db = SeqDb::open_search_in(b"dbh", SeqType::Protein, &path).unwrap();
        assert_eq!(db.max_length(), 0);
        assert_eq!(db.protein_length_i32(0).unwrap(), i32::MIN);
        db.check_protein_range(0).unwrap();
        let blk = db
            .read_protein_search_block(BlockRange {
                volume: 0,
                first: 0,
                end: 1,
                first_ordinal: 0,
            })
            .unwrap();
        assert_empty_sentinels(&blk, 0);
        assert!(dbh.with_extension("psq").exists());

        // Check dbg tail slice reads
        let db = SeqDb::open_search_in(b"dbg", SeqType::Protein, &path).unwrap();
        assert_eq!(db.protein_length_i32(0).unwrap(), -2_147_483_638);
        let at = db.locate_oid(0).unwrap();
        assert_eq!(
            db.read_protein_slice_at(at, u64::from(TWO_POW_31) - 2, 12)
                .unwrap(),
            &tail20[8..]
        );
        assert_eq!(
            db.read_protein_slice_at(at, u64::from(TWO_POW_31) - 9, 19)
                .unwrap(),
            &tail20[1..]
        );

        // Check single-volume mixed [40, 2^31] reads 40 residues for OID 0 and zero bytes for OID 1
        let db = SeqDb::open_search_in(b"mix2", SeqType::Protein, &path).unwrap();
        assert_eq!(db.max_length(), 0);
        assert_eq!(db.x_protein_lengths_i32(0, 2).unwrap(), vec![40, i32::MIN]);
        let blk = db
            .read_protein_search_block(BlockRange {
                volume: 0,
                first: 0,
                end: 2,
                first_ordinal: 0,
            })
            .unwrap();
        assert_eq!(blk.protein(0), &residues);
        assert_empty_sentinels(&blk, 1);

        // Check [2^31, 40] (sparse hole before the 40-residue sequence)
        let s1 = 1u32 + TWO_POW_31 + 1;
        let s2 = s1 + 41;
        write_v4_index(
            &root.join("mix1.pin"),
            SeqType::Protein,
            u64::from(TWO_POW_31) + 40,
            TWO_POW_31,
            &[1, s1, s2],
            None,
        );
        write_sparse_psq(
            &root.join("mix1.psq"),
            u64::from(s2),
            &[(u64::from(s1), &residues)],
        );
        std::fs::write(root.join("mix1.phr"), b"").unwrap();
        let db = SeqDb::open_search_in(b"mix1", SeqType::Protein, &path).unwrap();
        assert_eq!(db.max_length(), 0);
        assert_eq!(db.x_protein_lengths_i32(0, 2).unwrap(), vec![i32::MIN, 40]);
        let blk = db
            .read_protein_search_block(BlockRange {
                volume: 0,
                first: 0,
                end: 2,
                first_ordinal: 0,
            })
            .unwrap();
        assert_empty_sentinels(&blk, 0);
        assert_eq!(blk.protein(1), &residues);

        // Check alias dba max_length == 40
        let db = SeqDb::open_search_in(b"dba", SeqType::Protein, &path).unwrap();
        assert_eq!(db.max_length(), 40);

        // Check LENGTH 0 OID-0 fallback clamps negative protein length to -1
        let db = SeqDb::open_search_in(b"dbh_len0", SeqType::Protein, &path).unwrap();
        let stats = db.statistics();
        assert_eq!(stats.setup_total, -1);
        assert_eq!(stats.setup_count, 1);

        // Check OIDLIST x_ScanTotals wraps Uint8 base_count, keeps longest=40, min_length=i32::MIN,
        // exact_totals matches, and scheduling included_length stays 40 (not 18 exabytes).
        let db = SeqDb::open_search_in(b"dba_ol", SeqType::Protein, &path).unwrap();
        let expected_total = 0u64
            .wrapping_add(i64::from(i32::MIN) as u64)
            .wrapping_add(40);
        assert_eq!(db.total_length(), expected_total);
        assert_eq!(db.max_length(), 40);
        assert_eq!(db.min_length(), i32::MIN);
        assert_eq!(db.exact_totals().unwrap(), (expected_total, 40));
        assert_eq!(db.included_length().unwrap(), 40);
    }
}

#[test]
fn mt_negative_length_stop_and_irregular_nucleotide_ambiguity_offsets() {
    use crate::seqdb::{BlockRange, SearchPath, SeqDb, SeqDbError};
    use crate::utils::test_dir::TempDir;

    let dir = TempDir::unique("mt-neg-amb-unit");
    let root: &std::path::Path = &dir;
    let path = SearchPath::new(root, None, None);

    // 1. Protein volume with non-ascending sequence offset at OID 2 (len = -10).
    write_v4_index(
        &root.join("pneg.pin"),
        SeqType::Protein,
        21,
        13,
        &[1, 6, 20, 11, 16],
        None,
    );
    write_sparse_psq(
        &root.join("pneg.psq"),
        20,
        &[(1, &[1, 2, 3, 4]), (11, &[5, 6, 7, 8])],
    );
    std::fs::write(root.join("pneg.phr"), b"").unwrap();

    let mut db = SeqDb::open_search_in(b"pneg", SeqType::Protein, &path).unwrap();
    let full = BlockRange {
        volume: 0,
        first: 0,
        end: 4,
        first_ordinal: 0,
    };
    let blk_1t = db.read_protein_search_block(full).unwrap();
    assert!(!blk_1t.stopped_at_negative());
    assert_eq!(blk_1t.len(), 4);
    assert_empty_sentinels(&blk_1t, 2);
    assert_eq!(blk_1t.protein(3), &[5, 6, 7, 8]);

    db.set_mt_stop_at_negative(true);
    let blk_mt = db.read_protein_search_block(full).unwrap();
    assert!(blk_mt.stopped_at_negative());
    assert_eq!(blk_mt.len(), 2);
    assert_eq!(blk_mt.protein(0), &[1, 2, 3, 4]);

    // 1b. Filtered protein DB where the negative-length OID (OID 2) is excluded by a list:
    // NCBI's `CSeqDBImpl::x_FillSeqBuffer` still stops at OID 2 before applying the OID bitmap.
    // Block 0 (OIDs 0..2) precedes the gap; block 1 (OID 3) checks its preceding gap [2..3) and
    // returns an empty stopped block.
    let mut fdb = SeqDb::open_search_in(b"pneg", SeqType::Protein, &path).unwrap();
    fdb.x_keep_oids(|oid| oid != 2);
    assert_eq!(fdb.num_seqs(), 3);
    fdb.set_mt_stop_at_negative(true);
    let r0 = fdb.next_block(0, 1 << 20).unwrap();
    let b0 = fdb.read_protein_search_block(r0).unwrap();
    assert!(!b0.stopped_at_negative());
    assert_eq!(b0.len(), 2);
    let r1 = fdb.next_block(r0.end_ordinal(), 1 << 20).unwrap();
    let b1 = fdb.read_protein_search_block(r1).unwrap();
    assert!(b1.stopped_at_negative());
    assert_eq!(b1.len(), 0);
    // Excluded negative OID preceding all included OIDs stops block 0 with 0 scanned ordinals.
    let mut fdb_after = SeqDb::open_search_in(b"pneg", SeqType::Protein, &path).unwrap();
    fdb_after.x_keep_oids(|oid| oid == 3);
    fdb_after.set_mt_stop_at_negative(true);
    let r_after = fdb_after.next_block(0, 1 << 20).unwrap();
    let b_after = fdb_after.read_protein_search_block(r_after).unwrap();
    assert!(b_after.stopped_at_negative());
    assert_eq!(b_after.len(), 0);

    // 1c. `start == 0` protein cases:
    // - `start = 0, next = 1` -> length 0: skipped, NOT a stop in MT mode.
    // - `start = 0, next = 0` -> length -1: negative stop in MT mode.
    // - `start = 0, next = 5` -> length 4 > 0: refused (`NcbiCrashes`, reads before mapping).
    write_v4_index(
        &root.join("pstart0_len0.pin"),
        SeqType::Protein,
        4,
        4,
        &[0, 1, 6],
        None,
    );
    write_sparse_psq(&root.join("pstart0_len0.psq"), 6, &[(1, &[1, 2, 3, 4])]);
    std::fs::write(root.join("pstart0_len0.phr"), b"").unwrap();
    let mut s0_len0 = SeqDb::open_search_in(b"pstart0_len0", SeqType::Protein, &path).unwrap();
    assert_eq!(
        s0_len0
            .protein_sequence_len_at(s0_len0.locate_oid(0).unwrap())
            .unwrap(),
        0
    );
    s0_len0.set_mt_stop_at_negative(true);
    let blk_s0_len0 = s0_len0
        .read_protein_search_block(BlockRange {
            volume: 0,
            first: 0,
            end: 2,
            first_ordinal: 0,
        })
        .unwrap();
    assert!(!blk_s0_len0.stopped_at_negative());
    assert_eq!(blk_s0_len0.len(), 2);
    assert_empty_sentinels(&blk_s0_len0, 0);
    assert_eq!(blk_s0_len0.protein(1), &[1, 2, 3, 4]);

    write_v4_index(
        &root.join("pstart0_pos.pin"),
        SeqType::Protein,
        4,
        4,
        &[0, 5],
        None,
    );
    write_sparse_psq(&root.join("pstart0_pos.psq"), 6, &[(0, &[1, 2, 3, 4])]);
    std::fs::write(root.join("pstart0_pos.phr"), b"").unwrap();
    let s0_pos = SeqDb::open_search_in(b"pstart0_pos", SeqType::Protein, &path).unwrap();
    assert!(matches!(
        s0_pos.protein_sequence_len_at(s0_pos.locate_oid(0).unwrap()),
        Err(SeqDbError::NcbiCrashes { .. })
    ));
    assert!(matches!(
        s0_pos.read_protein_search_sequence(0),
        Err(SeqDbError::NcbiCrashes { .. })
    ));

    // 1d. Poked `start` past `.psq` EOF with non-positive length does not fail range checks before
    // skipping/stopping, matching NCBI pointer-arithmetic-only `x_GetSequence`.
    // OID 0: start=100_000, next=100_001 -> len=0 (skip, no stop);
    // OID 1: start=100_001, next=1 -> len=-100_001 (skip at 1t, stop at MT);
    // OID 2: start=1, next=6 -> len=4 (regular).
    write_v4_index(
        &root.join("ppast.pin"),
        SeqType::Protein,
        4,
        4,
        &[100_000, 100_001, 1, 6],
        None,
    );
    write_sparse_psq(&root.join("ppast.psq"), 6, &[(1, &[1, 2, 3, 4])]);
    std::fs::write(root.join("ppast.phr"), b"").unwrap();
    let mut ppast = SeqDb::open_search_in(b"ppast", SeqType::Protein, &path).unwrap();
    let ppast_range = BlockRange {
        volume: 0,
        first: 0,
        end: 3,
        first_ordinal: 0,
    };
    let ppast_1t = ppast.read_protein_search_block(ppast_range).unwrap();
    assert!(!ppast_1t.stopped_at_negative());
    assert_eq!(ppast_1t.len(), 3);
    assert!(ppast_1t.protein(0).is_empty());
    assert!(ppast_1t.protein(1).is_empty());
    assert_eq!(ppast_1t.protein(2), &[1, 2, 3, 4]);
    ppast.set_mt_stop_at_negative(true);
    let ppast_mt = ppast.read_protein_search_block(ppast_range).unwrap();
    assert!(ppast_mt.stopped_at_negative());
    assert_eq!(ppast_mt.len(), 1);

    // 2. Nucleotide volume with OID 0 regular (4 bases), OID 1 amb == start (len = -4),
    // OID 2 amb > next (reads 2 bytes past next, defers AMBIGUITY_DATA_ERROR), OID 3 regular.
    // seq_offsets = [0, 2, 4, 6, 8], amb = [2, 2, 8, 8, 8], nsq len = 8.
    write_v4_index(
        &root.join("nirr.nin"),
        SeqType::Nucleotide,
        16,
        4,
        &[0, 2, 4, 6, 8],
        Some(&[2, 2, 8, 8, 8]),
    );
    std::fs::write(
        root.join("nirr.nsq"),
        [0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
    )
    .unwrap();
    std::fs::write(root.join("nirr.nhr"), b"").unwrap();

    let mut ndb = SeqDb::open_search_in(b"nirr", SeqType::Nucleotide, &path).unwrap();
    assert_eq!(ndb.nucleotide_length_i32(0).unwrap(), 4);
    assert_eq!(ndb.nucleotide_length_i32(1).unwrap(), -4);
    assert_eq!(ndb.nucleotide_length_i32(2).unwrap(), 12);
    assert!(ndb
        .nucleotide_ambiguity(2)
        .unwrap_err()
        .is_ambiguity_data_error());

    let nrange = BlockRange {
        volume: 0,
        first: 0,
        end: 4,
        first_ordinal: 0,
    };
    let nblk_1t = ndb.read_nucleotide_search_block(nrange).unwrap();
    assert!(!nblk_1t.stopped_at_negative());
    assert_eq!(nblk_1t.len(), 4);
    assert_eq!(nblk_1t.nucleotide_len(1).unwrap(), 0);
    assert_eq!(nblk_1t.nucleotide_len(2).unwrap(), 12);
    assert!(nblk_1t
        .nucleotide_blastna(2)
        .unwrap_err()
        .is_ambiguity_data_error());

    ndb.set_mt_stop_at_negative(true);
    let nblk_mt = ndb.read_nucleotide_search_block(nrange).unwrap();
    assert!(nblk_mt.stopped_at_negative());
    assert_eq!(nblk_mt.len(), 1);

    // 2b. Filtered nucleotide DB where negative OID 1 is excluded still stops at block 1 (after
    // checking excluded gap [1..2)), and `exact_totals` sums signed `GetSeqLengthExact` lengths.
    let mut fndb = SeqDb::open_search_in(b"nirr", SeqType::Nucleotide, &path).unwrap();
    fndb.x_keep_oids(|oid| oid != 1);
    fndb.set_mt_stop_at_negative(true);
    let nr0 = fndb.next_block(0, 1 << 20).unwrap();
    let nb0 = fndb.read_nucleotide_search_block(nr0).unwrap();
    assert!(!nb0.stopped_at_negative());
    assert_eq!(nb0.len(), 1);
    let nr1 = fndb.next_block(nr0.end_ordinal(), 1 << 20).unwrap();
    let nb1 = fndb.read_nucleotide_search_block(nr1).unwrap();
    assert!(nb1.stopped_at_negative());
    assert_eq!(nb1.len(), 0);

    let mut fndb_all = SeqDb::open_search_in(b"nirr", SeqType::Nucleotide, &path).unwrap();
    fndb_all.x_keep_oids(|_| true);
    let (exact_total, exact_longest) = fndb_all.exact_totals().unwrap();
    let want_total = 0u64
        .wrapping_add(4)
        .wrapping_add((-4i64) as u64)
        .wrapping_add(12)
        .wrapping_add(4);
    assert_eq!((exact_total, exact_longest), (want_total, 12));

    // 2c. Poked nucleotide `start` past `.nsq` EOF with `whole_bytes < 0` and `last_pos` in-file:
    // NCBI reads only `file[start + whole_bytes]`, so `nucleotide_length_i32` and search block
    // reads cleanly stop/skip without checking `file[start]`.
    write_v4_index(
        &root.join("npast.nin"),
        SeqType::Nucleotide,
        4,
        4,
        &[100_000, 2, 4],
        Some(&[2, 4, 4]),
    );
    std::fs::write(root.join("npast.nsq"), [0u8; 4]).unwrap();
    std::fs::write(root.join("npast.nhr"), b"").unwrap();
    let mut npast = SeqDb::open_search_in(b"npast", SeqType::Nucleotide, &path).unwrap();
    // OID 0: whole_bytes = (2 - 100_000 - 1) = -99_999; last_pos = 100_000 + (-99_999) = 1.
    assert_eq!(npast.nucleotide_length_i32(0).unwrap(), -399_996);
    let npast_range = BlockRange {
        volume: 0,
        first: 0,
        end: 2,
        first_ordinal: 0,
    };
    let npast_1t = npast.read_nucleotide_search_block(npast_range).unwrap();
    assert!(!npast_1t.stopped_at_negative());
    assert_eq!(npast_1t.len(), 2);
    assert_eq!(npast_1t.nucleotide_len(0).unwrap(), 0);
    assert_eq!(npast_1t.nucleotide_len(1).unwrap(), 4);
    npast.set_mt_stop_at_negative(true);
    let npast_mt = npast.read_nucleotide_search_block(npast_range).unwrap();
    assert!(npast_mt.stopped_at_negative());
    assert_eq!(npast_mt.len(), 0);

    // 2d. Excluded nucleotide OID whose last-byte offset lies past the mapped 4 KiB page with
    // positive `whole_bytes` below 2^29 must fail the MT gap check because NCBI reads that tail
    // byte for every OID it scans (`seqdbvol.cpp:1835`).
    write_v4_index(
        &root.join("npast_excl.nin"),
        SeqType::Nucleotide,
        8,
        4,
        &[10_000, 0, 2],
        Some(&[10_002, 2, 2]),
    );
    std::fs::write(root.join("npast_excl.nsq"), [0u8; 4]).unwrap();
    std::fs::write(root.join("npast_excl.nhr"), b"").unwrap();
    let mut npast_excl = SeqDb::open_search_in(b"npast_excl", SeqType::Nucleotide, &path).unwrap();
    npast_excl.x_keep_oids(|oid| oid == 1);
    npast_excl.set_mt_stop_at_negative(true);
    let r_excl = npast_excl.next_block(0, 1 << 20).unwrap();
    let err_excl = npast_excl.read_nucleotide_search_block(r_excl).unwrap_err();
    assert!(
        matches!(
            &err_excl,
            SeqDbError::Format { message, .. }
                if message == "byte range 10001..10002 passes the end of the file (4 bytes)"
        ),
        "2d: {err_excl:?}"
    );

    // 3. Nucleotide volume with amb[0] == 0 -> reads before mmap -> NcbiCrashes refusal.
    write_v4_index(
        &root.join("namb0.nin"),
        SeqType::Nucleotide,
        4,
        4,
        &[0, 2],
        Some(&[0, 2]),
    );
    std::fs::write(root.join("namb0.nsq"), [0u8; 2]).unwrap();
    std::fs::write(root.join("namb0.nhr"), b"").unwrap();
    let zdb = SeqDb::open_search_in(b"namb0", SeqType::Nucleotide, &path).unwrap();
    assert!(matches!(
        zdb.nucleotide_length_i32(0),
        Err(SeqDbError::NcbiCrashes { .. })
    ));
    assert!(matches!(
        zdb.read_nucleotide_search_sequence(0),
        Err(SeqDbError::NcbiCrashes { .. })
    ));
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:324-373,684-710
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1805-1838
/// The excluded OIDs after the last included one: NCBI's multithreaded DB-split scan of a
/// list-filtered database reads them up to the iteration range's end, stopping at the first
/// negative length, and its read of a nucleotide's last byte past the mapping crashes. The block
/// that holds the last included OID returns that read's error (LOSAT's refusal); a list that
/// includes no OID gets it from `mt_check_scan_without_blocks`. Without the multithreaded stop
/// (1 thread, `-mt_mode 1`) nothing past the included OIDs is read.
#[test]
fn mt_scan_reads_the_excluded_oids_after_the_last_included_one() {
    use crate::seqdb::{SearchPath, SeqDb, SeqDbError};
    use crate::utils::test_dir::TempDir;

    let dir = TempDir::unique("mt-tail-excl-unit");
    let root: &std::path::Path = &dir;
    let path = SearchPath::new(root, None, None);
    let past_end = |err: &SeqDbError, range: &str, len: u64| {
        matches!(
            err,
            SeqDbError::Format { message, .. }
                if *message == format!("byte range {range} passes the end of the file ({len} bytes)")
        )
    };

    // 1. OID 0 regular (bytes 0..2); OID 1 (start 2, amb 10_002): `whole_bytes` 9_999, last byte
    // 10_001, past the 4-byte `.nsq` and its mapped page.
    write_v4_index(
        &root.join("ntail.nin"),
        SeqType::Nucleotide,
        8,
        4,
        &[0, 2, 4],
        Some(&[2, 10_002, 4]),
    );
    std::fs::write(root.join("ntail.nsq"), [0u8; 4]).unwrap();
    std::fs::write(root.join("ntail.nhr"), b"").unwrap();
    let open = || SeqDb::open_search_in(b"ntail", SeqType::Nucleotide, &path).unwrap();

    let mut db = open();
    db.x_keep_oids(|oid| oid == 0);
    let range = db.next_block(0, 1 << 20).unwrap();
    assert_eq!(range.end_ordinal(), db.num_seqs());
    let at_1t = db.read_nucleotide_search_block(range).unwrap();
    assert!(!at_1t.stopped_at_negative());
    assert_eq!(at_1t.len(), 1);
    db.mt_check_scan_without_blocks().unwrap();
    db.set_mt_stop_at_negative(true);
    let err = db.read_nucleotide_search_block(range).unwrap_err();
    assert!(past_end(&err, "10001..10002", 4), "1: {err:?}");
    let mut kept = Vec::new();
    let err = db
        .x_read_nucleotide_search_block_into(range, &mut kept)
        .unwrap_err();
    assert!(past_end(&err, "10001..10002", 4), "1 (into): {err:?}");
    // The list includes OID 0, so the block holds the check, not the scan without blocks.
    db.mt_check_scan_without_blocks().unwrap();

    // 2. A list that includes no OID: no block, but the scan reads OIDs 0 and 1.
    let mut none = open();
    none.x_keep_oids(|_| false);
    assert_eq!(none.num_seqs(), 0);
    assert!(none.next_block(0, 1 << 20).is_none());
    none.mt_check_scan_without_blocks().unwrap();
    none.set_mt_stop_at_negative(true);
    let err = none.mt_check_scan_without_blocks().unwrap_err();
    assert!(past_end(&err, "10001..10002", 4), "2: {err:?}");
    // Without a list there is no scan without blocks to check.
    let mut whole = open();
    whole.set_mt_stop_at_negative(true);
    whole.mt_check_scan_without_blocks().unwrap();

    // 3. A negative length before the bad OID ends NCBI's scan there: OID 1 (amb == start,
    // `whole_bytes` -1) stops it before OID 2 (last byte 10_003). The block is not stopped (the
    // negative length is past its OIDs).
    write_v4_index(
        &root.join("ntailneg.nin"),
        SeqType::Nucleotide,
        8,
        4,
        &[0, 2, 4, 6],
        Some(&[2, 2, 10_004, 6]),
    );
    std::fs::write(root.join("ntailneg.nsq"), [0u8; 6]).unwrap();
    std::fs::write(root.join("ntailneg.nhr"), b"").unwrap();
    let mut neg = SeqDb::open_search_in(b"ntailneg", SeqType::Nucleotide, &path).unwrap();
    neg.x_keep_oids(|oid| oid == 0);
    neg.set_mt_stop_at_negative(true);
    let range = neg.next_block(0, 1 << 20).unwrap();
    let block = neg.read_nucleotide_search_block(range).unwrap();
    assert!(!block.stopped_at_negative());
    assert_eq!(block.len(), 1);
    let mut neg_none = SeqDb::open_search_in(b"ntailneg", SeqType::Nucleotide, &path).unwrap();
    neg_none.x_keep_oids(|_| false);
    neg_none.set_mt_stop_at_negative(true);
    neg_none.mt_check_scan_without_blocks().unwrap();

    // 4. A last block that stops at a negative length inside it reads nothing after it: OID 1
    // stops the block of OIDs 0..3 before the excluded OID 3 (last byte 10_005).
    write_v4_index(
        &root.join("ntailstop.nin"),
        SeqType::Nucleotide,
        8,
        4,
        &[0, 2, 4, 6, 8],
        Some(&[2, 2, 6, 10_006, 8]),
    );
    std::fs::write(root.join("ntailstop.nsq"), [0u8; 8]).unwrap();
    std::fs::write(root.join("ntailstop.nhr"), b"").unwrap();
    let mut stop = SeqDb::open_search_in(b"ntailstop", SeqType::Nucleotide, &path).unwrap();
    stop.x_keep_oids(|oid| oid != 3);
    stop.set_mt_stop_at_negative(true);
    let range = stop.next_block(0, 1 << 20).unwrap();
    assert_eq!(range.end_ordinal(), 3);
    let block = stop.read_nucleotide_search_block(range).unwrap();
    assert!(block.stopped_at_negative());
    assert_eq!(block.len(), 1);
}
