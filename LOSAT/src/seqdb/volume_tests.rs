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
    assert_eq!(db.x_protein_mt_scan_count().unwrap(), 2);
    let blk_mt = db.read_protein_search_block(full).unwrap();
    assert!(blk_mt.stopped_at_negative());
    assert_eq!(blk_mt.len(), 2);
    assert_eq!(blk_mt.protein(0), &[1, 2, 3, 4]);

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
