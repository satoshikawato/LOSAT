//! Damaged copies of the query and subject inputs: four FASTA files, two of them also compressed
//! as gzip, bzip2 and zstd, and two zstd files written by the zstd CLI. Each input is cut at many
//! lengths, then damaged in rounds drawn from a fixed seed: bits flipped, a 4-byte word rewritten
//! (a length, a size or a checksum in the compressed forms), or some of its first 64 bytes
//! written (the headers of the compressed forms, the first defline of FASTA). Each case is read
//! as the programs read `-query` and `-subject`: a compressed input through
//! [`Compression::decoder`], up to `MAX_INFLATED` bytes, and the bytes through the FASTA reader
//! (`read_fasta`), with the Seq-id checks of the data loaders off and on. Any result but a panic
//! or a hang is accepted: an error or, where the damage leaves a readable input, records. A
//! panic is caught and named with its input and case, and the test fails listing them all. Each
//! input must also have seeded rounds whose outcome differs from the intact input's, and each
//! input that the tests' own encoders compress must have cases that inflate to other bytes than
//! the intact input (a cut at the end of its first member, stream or frame among them), so that
//! the damage reaches the FASTA reader through the decoder.
//!
//! The FASTA reader also has a seeded test of generated inputs
//! (`fasta_reader::tests::bytes_read_as_a_file_with_the_same_bytes`), which compares reading
//! from bytes, from a file and a byte at a time. The form here is after
//! `corrupt_files_give_errors_not_panics` (`seqdb/lmdb_read.rs`); the generator is the xorshift
//! of `seqdb/isam_tests.rs`.

use std::io::Read;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::{Path, PathBuf};

use super::compressed::test_data::{bzip2_streams, gzip_members, zstd_frames};
use super::compressed::Compression;
use super::fasta_reader::{
    read_queries, read_subjects, x_read_queries_packed, FastaInputSource, FastaRecord, QueryEnd,
    QueryRecords, ReaderConfig,
};
use super::x_packed_subjects::x_read_packed_subjects;

/// Seeded rounds of damage per input, after the cuts.
const ROUNDS: usize = 2000;
/// The most bytes inflated from a damaged compressed input.
const MAX_INFLATED: u64 = 1 << 20;

fn tests_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests")
}

/// Marsaglia's xorshift64, as in `seqdb/isam_tests.rs`, seeded from the bytes of `name`.
struct Rng(u64);

impl Rng {
    fn seeded(name: &str) -> Self {
        let seed = name.bytes().fold(0x5eed_0000_0000_0002u64, |seed, byte| {
            seed.rotate_left(7) ^ u64::from(byte)
        });
        Self(seed | 1)
    }

    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n.max(1) as u64) as usize
    }
}

/// How a seeded round damages an input.
#[derive(Clone, Copy, Debug)]
enum Damage {
    /// The input cut to a length drawn below its length.
    Cut,
    /// One to four bits flipped anywhere.
    Bits,
    /// A 4-byte word at a multiple of 4 set, in either byte order, to 0, 1, 2^31 - 1, 2^31,
    /// 2^32 - 1, the input length minus one, the input length, plus one or doubled, or the old
    /// value plus or minus one.
    Word,
    /// One to four of the first 64 bytes set to drawn values.
    Head,
}

fn damage(orig: &[u8], rng: &mut Rng, round: usize) -> (Damage, Vec<u8>) {
    let mut bytes = orig.to_vec();
    let kind = [Damage::Cut, Damage::Bits, Damage::Word, Damage::Head][round % 4];
    if bytes.is_empty() {
        return (kind, bytes);
    }
    match kind {
        Damage::Cut => bytes.truncate(rng.below(orig.len())),
        Damage::Bits => {
            for _ in 0..1 + rng.below(4) {
                let at = rng.below(bytes.len());
                bytes[at] ^= 1 << rng.below(8);
            }
        }
        Damage::Word => {
            if bytes.len() >= 4 {
                let at = rng.below(bytes.len() / 4) * 4;
                let big_endian = rng.next().is_multiple_of(2);
                let word: [u8; 4] = bytes[at..at + 4].try_into().unwrap();
                let old = if big_endian {
                    u32::from_be_bytes(word)
                } else {
                    u32::from_le_bytes(word)
                };
                let len = orig.len() as u32;
                let values = [
                    0,
                    1,
                    i32::MAX as u32,
                    1 << 31,
                    u32::MAX,
                    len - 1,
                    len,
                    len + 1,
                    len.wrapping_mul(2),
                    old.wrapping_add(1),
                    old.wrapping_sub(1),
                ];
                let value = values[rng.below(values.len())];
                let new = if big_endian {
                    value.to_be_bytes()
                } else {
                    value.to_le_bytes()
                };
                bytes[at..at + 4].copy_from_slice(&new);
            }
        }
        Damage::Head => {
            let head = bytes.len().min(64);
            for _ in 0..1 + rng.below(4) {
                let at = rng.below(head);
                bytes[at] = rng.next() as u8;
            }
        }
    }
    (kind, bytes)
}

/// What one case gave: the bytes inflated (`None` for plain FASTA or when the decoder failed),
/// the title and letter counts of the records each reading gave, the packed letter counts, and
/// the errors. A case whose summary differs from the intact input's read otherwise.
#[derive(Default, PartialEq, Eq, Debug)]
struct Summary {
    inflated: Option<usize>,
    records: Vec<(usize, usize)>,
    packed: Vec<usize>,
    errors: usize,
}

impl Summary {
    fn records(&mut self, records: &[FastaRecord]) {
        self.records.extend(
            records
                .iter()
                .map(|record| (record.title.len(), record.sequence.len())),
        );
    }

    fn queries(&mut self, queries: QueryRecords) {
        self.records(&queries.records);
        if let QueryEnd::Error { error, .. } = queries.end {
            let _ = error.to_string();
            self.errors += 1;
        }
    }

    fn error(&mut self, error: impl ToString) {
        let _ = error.to_string();
        self.errors += 1;
    }
}

/// Reads `bytes` as the programs read FASTA, with the Seq-id checks of the data loaders off and
/// on: `-query` through `read_queries`, and for nucleotide letters also through
/// `x_read_queries_packed` as BLASTN does; `-subject` through `read_subjects`, and for
/// nucleotide letters also through `x_read_packed_subjects` as TBLASTN and TBLASTX do.
fn read_fasta(bytes: &[u8], protein: bool, summary: &mut Summary) {
    let program = if protein { "BLASTP" } else { "BLASTN" };
    let mut warn = |_: &[u8]| -> std::io::Result<()> { Ok(()) };
    for data_loaders in [false, true] {
        let query = || {
            FastaInputSource::from_bytes(bytes, ReaderConfig::query(program, protein, data_loaders))
        };
        let subject = || {
            FastaInputSource::from_bytes(
                bytes,
                ReaderConfig::subject(program, protein, data_loaders),
            )
        };
        summary.queries(read_queries(&mut query()));
        match read_subjects(&mut subject(), None, &mut warn) {
            Ok(records) => summary.records(&records),
            Err(error) => summary.error(error),
        }
        if protein {
            continue;
        }
        let (queries, packed) = x_read_queries_packed(&mut query());
        summary.queries(queries);
        summary
            .packed
            .extend(packed.iter().map(|(_, letters)| letters.len()));
        // One thread: the packing runs on this thread, so a panic in it reaches `catch_unwind`
        // (with more, a panic of the packing thread comes back as an error).
        match x_read_packed_subjects(&mut subject(), None, &mut warn, 1, "TBLASTN") {
            Ok((records, packed, _)) => {
                summary.records(&records);
                summary
                    .packed
                    .extend(packed.iter().flatten().map(|letters| letters.len()));
            }
            Err(error) => summary.error(error),
        }
    }
}

/// Reads `bytes` as an input compressed by `compression` (or plain FASTA when `None`).
fn exercise(bytes: &[u8], compression: Option<Compression>, protein: bool) -> Summary {
    let mut summary = Summary::default();
    let Some(compression) = compression else {
        read_fasta(bytes, protein, &mut summary);
        return summary;
    };
    let mut inflated = Vec::new();
    match compression
        .decoder(bytes)
        .take(MAX_INFLATED)
        .read_to_end(&mut inflated)
    {
        Ok(_) => {
            summary.inflated = Some(inflated.len());
            read_fasta(&inflated, protein, &mut summary);
        }
        Err(error) => summary.error(error),
    }
    summary
}

/// One input of a sweep.
struct Input<'a> {
    name: String,
    bytes: &'a [u8],
    compression: Option<Compression>,
    protein: bool,
    /// The plain FASTA that the input inflates to.
    plain: &'a [u8],
    /// Lengths to cut the input at besides the spread of cuts (the end of a first member).
    cuts: Vec<usize>,
    /// Whether some case must inflate to other bytes than the intact input: false for the
    /// zstd CLI files, whose one frame's checksum refuses every other content.
    inflates_otherwise: bool,
}

/// Damages the input (cuts, then seeded rounds) and reads each case. Returns the problems: the
/// cases that panicked, and a note when the damage did not reach the reader.
fn sweep(input: &Input) -> Vec<String> {
    let Input { name, bytes, .. } = input;
    let intact = exercise(bytes, input.compression, input.protein);
    assert!(
        intact.errors == 0 && !intact.records.is_empty(),
        "{name}: the intact input must read: {intact:?}"
    );
    if input.compression.is_some() {
        assert_eq!(
            intact.inflated,
            Some(input.plain.len()),
            "{name}: inflates whole"
        );
    }
    let mut rng = Rng::seeded(name);
    let step = bytes.len().div_ceil(256).max(1);
    let mut cases: Vec<(Option<(usize, Damage)>, Vec<u8>)> = (0..bytes.len())
        .step_by(step)
        .chain(input.cuts.iter().copied())
        .map(|cut| (None, bytes[..cut].to_vec()))
        .collect();
    for round in 0..ROUNDS {
        let (kind, damaged) = damage(bytes, &mut rng, round);
        cases.push((Some((round, kind)), damaged));
    }
    let mut problems = Vec::new();
    let (mut seeded_changed, mut inflated_otherwise) = (0, 0);
    for (case, damaged) in cases {
        let outcome = catch_unwind(AssertUnwindSafe(|| {
            exercise(&damaged, input.compression, input.protein)
        }));
        match outcome {
            Ok(summary) => {
                if case.is_some() && summary != intact {
                    seeded_changed += 1;
                }
                if summary.inflated.is_some() && summary != intact {
                    inflated_otherwise += 1;
                }
            }
            Err(_) => {
                let what = match case {
                    None => format!("cut to {} bytes", damaged.len()),
                    Some((round, kind)) => format!("round {round} ({kind:?})"),
                };
                problems.push(format!("{name}: {what} panicked"));
            }
        }
    }
    if seeded_changed == 0 || (input.inflates_otherwise && inflated_otherwise == 0) {
        problems.push(format!(
            "{name}: the damage did not reach the reader (seeded rounds read otherwise \
             {seeded_changed} of {ROUNDS}, cases inflated otherwise {inflated_otherwise})"
        ));
    }
    problems
}

fn fixture(path: &str) -> Vec<u8> {
    std::fs::read(tests_dir().join(path)).unwrap()
}

fn assert_no_problems(problems: Vec<String>) {
    assert!(
        problems.is_empty(),
        "{} problems:\n{}",
        problems.len(),
        problems.join("\n")
    );
}

/// The FASTA files of the database fixtures, plain.
#[test]
fn damaged_fasta_inputs_never_panic() {
    let mut problems = Vec::new();
    for (path, protein) in [
        ("seqdb/nucl.fna", false),
        ("seqdb/mask_nucl.fna", false),
        ("seqdb/prot.faa", true),
        ("seqdb/mask_prot.faa", true),
    ] {
        let bytes = fixture(path);
        problems.extend(sweep(&Input {
            name: path.to_string(),
            bytes: &bytes,
            compression: None,
            protein,
            plain: &bytes,
            cuts: Vec::new(),
            inflates_otherwise: false,
        }));
    }
    assert_no_problems(problems);
}

/// Two of the FASTA files as two gzip members, two bzip2 streams and two zstd frames of raw
/// blocks (as `cat` joins them), and the zstd CLI's two files of one frame with a checksum each
/// (one with a content size, one without).
#[test]
fn damaged_compressed_inputs_never_panic() {
    let mut problems = Vec::new();
    for (path, protein) in [("seqdb/nucl.fna", false), ("seqdb/prot.faa", true)] {
        let plain = fixture(path);
        let half = plain.len() / 2;
        let parts = [&plain[..half], &plain[half..]];
        let forms: [(Compression, fn(&[&[u8]]) -> Vec<u8>); 3] = [
            (Compression::Gzip, gzip_members),
            (Compression::Bzip2, bzip2_streams),
            (Compression::Zstd, zstd_frames),
        ];
        for (compression, encode) in forms {
            let bytes = encode(&parts);
            problems.extend(sweep(&Input {
                name: format!("{path}{}", compression.suffix()),
                bytes: &bytes,
                compression: Some(compression),
                protein,
                plain: &plain,
                cuts: vec![encode(&parts[..1]).len()],
                inflates_otherwise: true,
            }));
        }
    }
    for (path, plain, protein) in [
        ("compressed/mask_nucl.fna.zst", "seqdb/mask_nucl.fna", false),
        ("compressed/prot.faa.zst", "seqdb/prot.faa", true),
    ] {
        let bytes = fixture(path);
        problems.extend(sweep(&Input {
            name: path.to_string(),
            bytes: &bytes,
            compression: Some(Compression::Zstd),
            protein,
            plain: &fixture(plain),
            cuts: Vec::new(),
            inflates_otherwise: false,
        }));
    }
    assert_no_problems(problems);
}

fn inflate(compression: Compression, bytes: &[u8]) -> std::io::Result<Vec<u8>> {
    let mut out = Vec::new();
    compression.decoder(bytes).read_to_end(&mut out)?;
    Ok(out)
}

/// Every single-bit flip of the zstd CLI's files is an error or inflates to the intact bytes:
/// the frame's checksum (and its content size, where the header has one) refuse every other
/// content. Before the checks of `ZstdFrames::check_frame`, flips of the compressed block
/// inflated to other letters without an error. Frames with and without a checksum or a size
/// still inflate one after the other.
#[test]
fn zstd_frames_whose_checksum_or_size_does_not_match_are_errors() {
    let cli = [
        ("compressed/mask_nucl.fna.zst", "seqdb/mask_nucl.fna"),
        ("compressed/prot.faa.zst", "seqdb/prot.faa"),
    ];
    for (path, plain) in cli {
        let bytes = fixture(path);
        let plain = fixture(plain);
        assert_eq!(inflate(Compression::Zstd, &bytes).unwrap(), plain, "{path}");
        let (mut refused, mut checksum) = (0, 0);
        for at in 0..bytes.len() {
            for bit in 0..8 {
                let mut damaged = bytes.clone();
                damaged[at] ^= 1 << bit;
                match inflate(Compression::Zstd, &damaged) {
                    Ok(out) => assert!(
                        out == plain,
                        "{path} byte {at} bit {bit}: inflated {} bytes other than the intact ones",
                        out.len()
                    ),
                    Err(error) => {
                        refused += 1;
                        checksum += usize::from(error.to_string().contains("checksum"));
                    }
                }
            }
        }
        assert!(
            refused > 0 && checksum > 0,
            "{path}: refused {refused}, checksum {checksum}"
        );
    }
    // A frame of raw blocks without a size or a checksum (descriptor 0, a Window_Descriptor of
    // 1 KiB), then the CLI's frames with a checksum, with and without a size, one after the
    // other: each frame is checked alone.
    let short: &[u8] = b">q1\nACGT\n";
    let mut bare = 0xFD2F_B528u32.to_le_bytes().to_vec();
    bare.extend_from_slice(&[0, 0]);
    bare.extend_from_slice(&(1 | ((short.len() as u32) << 3)).to_le_bytes()[..3]);
    bare.extend_from_slice(short);
    assert_eq!(inflate(Compression::Zstd, &bare).unwrap(), short);
    let (nucl, prot) = (fixture(cli[0].0), fixture(cli[1].0));
    let joined = [&bare[..], &nucl, &prot, &nucl, &bare].concat();
    let plain = [
        short,
        &fixture(cli[0].1),
        &fixture(cli[1].1),
        &fixture(cli[0].1),
        short,
    ]
    .concat();
    assert_eq!(inflate(Compression::Zstd, &joined).unwrap(), plain);
    // A frame of raw blocks has no checksum; its size in the header is compared.
    let frame = zstd_frames(&[b">q1\nACGT\n"]);
    // skippable frame (12 bytes), magic (4), descriptor (1), then the 4-byte content size
    let size_at = 12 + 4 + 1;
    assert_eq!(&frame[size_at..size_at + 4], &9u32.to_le_bytes());
    let mut damaged = frame.clone();
    damaged[size_at..size_at + 4].copy_from_slice(&10u32.to_le_bytes());
    let error = inflate(Compression::Zstd, &damaged)
        .unwrap_err()
        .to_string();
    assert_eq!(
        error,
        "a zstd frame holds 9 bytes where its header gives 10"
    );
    // A size below the block's is an error too (of ruzstd or of the check).
    damaged[size_at..size_at + 4].copy_from_slice(&8u32.to_le_bytes());
    assert!(inflate(Compression::Zstd, &damaged).is_err());
}
