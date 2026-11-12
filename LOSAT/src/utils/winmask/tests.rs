//! Tests of the WindowMasker port.
//!
//! The fixtures under `LOSAT/tests/winmask` were written by NCBI BLAST+ 2.17.0 `windowmasker`
//! (`README.md` there): counts files from `-mk_counts` and `-outfmt interval` masks. The
//! format and lookup tests build small counts files in the layout of
//! `CSeqMaskerOstatOptBin::write_out` (`seq_masker_ostat_opt_bin.cpp:76-108`) and
//! `CSeqMaskerOstatBin` (`seq_masker_ostat_bin.cpp`).

use std::io::Cursor;
use std::path::{Path, PathBuf};

use super::istat::{discover, hash_code, reverse_complement};
use super::*;

fn fixture_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/winmask")
}

/// FASTA records as the windowmasker application masks them: the title line (`CMaskWriter::
/// IdToString` prints the defline of a record read without `-parse_seqids`) and the letters in
/// upper case (`CBioseq_Handle::eCoding_Iupac`). Records of length 0 are skipped
/// (`win_mask_app.cpp:254-256`).
fn read_fasta(path: &Path) -> Vec<(String, Vec<u8>)> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines() {
        let line = line.trim_end_matches('\r');
        if let Some(title) = line.strip_prefix('>') {
            records.push((title.to_string(), Vec::new()));
        } else if let Some(last) = records.last_mut() {
            last.1.extend(
                line.bytes()
                    .filter(|b| !b.is_ascii_whitespace())
                    .map(|b| b.to_ascii_uppercase()),
            );
        }
    }
    records.retain(|r| !r.1.is_empty());
    records
}

/// The windowmasker application's flags that set `CSeqMasker` arguments
/// (`win_mask_config.cpp:331-335, 399-402`). `-window` goes into `Uint1 window_size`
/// (`win_mask_config.hpp:643`), so 256 is 0, the default.
fn params_from_flags(flags: &str) -> MaskerParams {
    let mut params = MaskerParams::default();
    if flags == "-" {
        return params;
    }
    let words: Vec<&str> = flags.split_whitespace().collect();
    for pair in words.chunks(2) {
        let value: u32 = pair[1].parse().expect("integer flag value");
        match pair[0] {
            "-window" => params.window_size = value as u8,
            "-t_extend" => params.textend = value,
            "-t_thres" => params.cutoff_score = value,
            "-t_high" => params.max_score = value,
            "-t_low" => params.min_score = value,
            "-set_t_high" => params.set_max_score = value,
            "-set_t_low" => params.set_min_score = value,
            other => panic!("unknown flag {other}"),
        }
    }
    params
}

/// NCBI reference (598d8ae6): c++/src/objtools/seqmasks_io/mask_writer_int.cpp:40-54
/// ```c
///     PrintId( bsh, parsed_id );
///     os << endl;
///     PrintMasks(os, mask);
/// ...
///     ITERATE(TMaskList, i, mask) {
///         os << i->first << " - " << i->second << "\n";
///     }
/// ```
fn interval_text(records: &[(String, Vec<u8>)], masker: &SeqMasker) -> String {
    let mut out = String::new();
    for (title, seq) in records {
        out.push('>');
        out.push_str(title);
        out.push('\n');
        for (a, b) in masker.mask(seq).unwrap_or_else(|e| panic!("{title}: {e}")) {
            out.push_str(&format!("{a} - {b}\n"));
        }
    }
    out
}

/// Runs every row of a manifest (`counts`, `fasta`, `flags`, `expected`, paths relative to the
/// manifest) and returns the rows whose text differs.
fn check_manifest(manifest: &Path) -> (usize, Vec<String>) {
    let dir = manifest.parent().unwrap();
    let text = std::fs::read_to_string(manifest).unwrap();
    let mut rows = 0;
    let mut failures = Vec::new();
    for line in text
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
    {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(cols.len(), 4, "{line}");
        let params = params_from_flags(cols[2]);
        let masker = SeqMasker::open(&dir.join(cols[0]), &params)
            .unwrap_or_else(|e| panic!("{}: {e}", cols[0]));
        let got = interval_text(&read_fasta(&dir.join(cols[1])), &masker);
        let want = std::fs::read_to_string(dir.join(cols[3])).unwrap();
        rows += 1;
        if got != want {
            let first = got
                .lines()
                .zip(want.lines())
                .position(|(a, b)| a != b)
                .unwrap_or(got.lines().count().min(want.lines().count()));
            failures.push(format!(
                "{} x {} [{}]: first difference at line {} (LOSAT {:?}, NCBI {:?})",
                cols[0],
                cols[1],
                cols[2],
                first + 1,
                got.lines().nth(first),
                want.lines().nth(first)
            ));
        }
    }
    (rows, failures)
}

#[test]
fn masks_equal_ncbi_windowmasker_intervals() {
    let (rows, failures) = check_manifest(&fixture_dir().join("manifest.tsv"));
    // 21 rows of the binary formats, and the 17 of mg1655_u9.obinary again with the ascii and
    // the optimized ascii counts (line NP M3d).
    assert_eq!(rows, 55);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// Compares counts files and queries outside the repository with NCBI's intervals.
///
/// `LOSAT_WINMASK_FULL_CHECK` names a manifest of the same columns as
/// `tests/winmask/manifest.tsv` (on the workstation: `~/losat-bench/npm/m1/full/manifest.tsv`,
/// counts of units of 12 to 16 letters).
#[test]
#[ignore = "needs LOSAT_WINMASK_FULL_CHECK (counts files outside the repository)"]
fn full_size_counts_match_ncbi() {
    let manifest = std::env::var("LOSAT_WINMASK_FULL_CHECK").expect("LOSAT_WINMASK_FULL_CHECK");
    let (rows, failures) = check_manifest(Path::new(&manifest));
    eprintln!("{rows} rows, {} differ", failures.len());
    assert!(rows > 0);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn t_low_below_the_file_warns_as_ncbi_does() {
    // `windowmasker -t_low 50` with a file whose t_low is 74: NCBI prints the warning through
    // its diagnostics ("Warning: " + text) and uses 74.
    let params = params_from_flags("-t_low 50");
    let masker = SeqMasker::open(&fixture_dir().join("counts/mg1655_u9.obinary"), &params).unwrap();
    let stderr =
        std::fs::read_to_string(fixture_dir().join("expected/mg1655_u9.p_low50.stderr")).unwrap();
    let warnings = masker.unit_counts().warnings();
    assert_eq!(warnings.len(), 1);
    assert_eq!(format!("Warning: {}\n", warnings[0]), stderr);
    assert_eq!(masker.unit_counts().min_count(), 74);
}

#[test]
fn window_below_the_unit_is_ncbis_validation_error() {
    let params = params_from_flags("-window 5");
    let err = match SeqMasker::open(&fixture_dir().join("counts/mg1655_u9.obinary"), &params) {
        Err(e) => e,
        Ok(_) => panic!("a window of 5 letters with units of 9 must fail"),
    };
    assert_eq!(
        err,
        WinMaskError::WindowSize {
            window: 5,
            unit: 9,
            warnings: Vec::new()
        }
    );
    let stderr =
        std::fs::read_to_string(fixture_dir().join("expected/mg1655_u9.p_window5.stderr")).unwrap();
    assert_eq!(
        stderr.lines().next().unwrap(),
        format!("Error: (CSeqMaskerException::validation error) {err}")
    );
}

#[test]
fn thresholds_come_from_the_file_unless_given() {
    let dir = fixture_dir();
    let ob = UnitCounts::open(
        &dir.join("counts/mg1655_u9.obinary"),
        &MaskerParams::default(),
    )
    .unwrap();
    assert_eq!(ob.unit_size(), 9);
    assert_eq!(
        (ob.min_count(), ob.textend(), ob.threshold(), ob.max_count()),
        (74, 147, 173, 210)
    );
    assert_eq!((ob.use_min_count(), ob.use_max_count()), (37, 210));
    assert!(ob.cache_bit_array().is_some());
    assert_eq!(
        ob.metadata()[0],
        b"##windowmasker-statistics-format-version:obinary 1.1.0"
    );
    let bin =
        UnitCounts::open(&dir.join("counts/dm2L_u9.binary"), &MaskerParams::default()).unwrap();
    assert_eq!(
        (
            bin.min_count(),
            bin.textend(),
            bin.threshold(),
            bin.max_count()
        ),
        (43, 106, 136, 184)
    );
    assert!(bin.cache_bit_array().is_none());
    let given = MaskerParams {
        textend: 100,
        cutoff_score: 120,
        max_score: 150,
        min_score: 80,
        set_max_score: 400,
        set_min_score: 10,
        window_size: 0,
    };
    let ob = UnitCounts::open(&dir.join("counts/mg1655_u9.obinary"), &given).unwrap();
    assert_eq!(
        (
            ob.min_count(),
            ob.textend(),
            ob.threshold(),
            ob.max_count(),
            ob.use_min_count(),
            ob.use_max_count()
        ),
        (80, 100, 120, 150, 10, 400)
    );
    assert!(ob.warnings().is_empty());
}

#[test]
fn a_unit_and_its_reverse_complement_share_a_count() {
    // `get_info` folds a unit to min(unit, reverse complement) before the lookup.
    let ob = UnitCounts::open(
        &fixture_dir().join("counts/mg1655_u9.obinary"),
        &MaskerParams::default(),
    )
    .unwrap();
    let mut nonzero = 0;
    for unit in (0..1u32 << 18).step_by(97) {
        let c = ob.true_count(unit).unwrap();
        assert_eq!(c, ob.true_count(reverse_complement(unit, 9)).unwrap());
        nonzero += (c != 0) as u32;
    }
    assert!(nonzero > 0);
}

#[test]
fn reverse_complement_of_two_bit_units() {
    // A=0 C=1 G=2 T=3, the first letter in the high bits.
    assert_eq!(reverse_complement(0b00_01_10_11, 4), 0b00_01_10_11); // ACGT
    assert_eq!(reverse_complement(0, 4), 0xFF); // AAAA -> TTTT
    assert_eq!(reverse_complement(0b00_00_01, 3), 0b10_11_11); // AAC -> GTT
    assert_eq!(reverse_complement(0, 16), 0xFFFF_FFFF);
}

#[test]
fn hash_code_splits_key_and_check_bits() {
    // key = bits [roff, roff + k); check = the bits above, shifted down to roff, plus the bits
    // below roff, as one byte.
    assert_eq!(hash_code(0b1011_0110, 4, 0), (0b0110, 0b1011));
    assert_eq!(hash_code(0b1011_0110, 4, 2), (0b1101, 0b10_10));
    // roff + k = 32 (units of 16): the shift by 32 is taken modulo 32, as NCBI's builds do.
    let unit = 0xDEAD_BEEFu32;
    let (key, check) = hash_code(unit, 26, 6);
    assert_eq!(key, (unit >> 6) & ((1 << 26) - 1));
    assert_eq!(check, ((unit << 6).wrapping_add(unit & 0x3F)) as u8);
}

/// A counts file in the optimized binary layout of `CSeqMaskerOstatOptBin::write_out`.
#[allow(clippy::too_many_arguments)]
fn obinary_file(
    metadata: Option<&[&str]>,
    fmt: u32,
    unit: u32,
    k: u32,
    roff: u32,
    bc: u32,
    t: [u32; 4],
    cba: Option<&[u32]>,
    ht: &[u32],
    vt: &[u16],
) -> Vec<u8> {
    let mut b = Vec::new();
    let w = |b: &mut Vec<u8>, x: u32| b.extend_from_slice(&x.to_ne_bytes());
    if let Some(md) = metadata {
        let mut block = Vec::new();
        for line in md {
            block.extend_from_slice(line.as_bytes());
            block.push(0);
        }
        w(&mut b, 3);
        w(&mut b, block.len() as u32);
        b.extend_from_slice(&block);
    }
    w(&mut b, fmt);
    w(&mut b, unit);
    w(&mut b, vt.len() as u32);
    w(&mut b, k);
    w(&mut b, roff);
    w(&mut b, bc);
    for x in t {
        w(&mut b, x);
    }
    if fmt == 2 {
        match cba {
            Some(cba) => {
                w(&mut b, 1);
                for &x in cba {
                    w(&mut b, x);
                }
            }
            None => w(&mut b, 0),
        }
    }
    for &x in ht {
        w(&mut b, x);
    }
    for &x in vt {
        b.extend_from_slice(&x.to_ne_bytes());
    }
    b
}

fn read_counts(bytes: Vec<u8>) -> Result<UnitCounts, WinMaskError> {
    UnitCounts::read(&mut Cursor::new(bytes), &MaskerParams::default())
}

#[test]
fn discover_stat_type_reads_the_first_word_after_the_metadata() {
    let kind = |bytes: &[u8]| discover(&mut Cursor::new(bytes.to_vec())).unwrap();
    let w = |x: u32| x.to_ne_bytes().to_vec();
    assert_eq!(kind(&w(0)).stat_type, StatType::Binary);
    assert_eq!(kind(&w(1)).stat_type, StatType::OBinary);
    assert_eq!(kind(&w(2)).stat_type, StatType::OBinary);
    assert_eq!(kind(b"").stat_type, StatType::Unknown);
    let with_md = [w(3), w(6), b"##a\0\0\0".to_vec(), w(2)].concat();
    let found = kind(&with_md);
    assert_eq!(found.stat_type, StatType::OBinary);
    assert_eq!(found.metadata, vec![b"##a".to_vec()]);
    assert_eq!(found.skip, 14);
    let text = kind(b"##windowmasker-statistics-format-version:oascii 1.1.0\n##x\nAAAA\n");
    assert_eq!(text.stat_type, StatType::OAscii);
    assert_eq!(text.metadata.len(), 2);
    assert_eq!(kind(b"##m\n0 12\n").stat_type, StatType::Ascii);
    assert_eq!(kind(b"##only\n").stat_type, StatType::Unknown);
    let err = discover_stat_type(Path::new("/nonexistent/wmasker.obinary")).unwrap_err();
    assert_eq!(
        err.to_string(),
        "could not open /nonexistent/wmasker.obinary"
    );
    match read_counts(b"##m\nAAAA\n".to_vec()) {
        Err(WinMaskError::Format {
            site: FormatSite::OAsciiFileTooShort,
            ..
        }) => {}
        Err(e) => panic!("{e}"),
        Ok(_) => panic!("a header of one line is too short"),
    }
}

#[test]
fn obinary_header_errors_are_ncbis() {
    let t = [1, 2, 3, 4];
    let msg = |bytes: Vec<u8>| match read_counts(bytes) {
        Err(e) => e.to_string(),
        Ok(_) => "ok".to_string(),
    };
    assert_eq!(
        msg(obinary_file(None, 1, 0, 1, 0, 1, t, None, &[0; 2], &[])),
        "unit size must be in [1,16]"
    );
    assert_eq!(
        msg(obinary_file(None, 1, 2, 4, 0, 1, t, None, &[0; 16], &[])),
        "hash key size must be in [1,2*unit_size - 1]"
    );
    assert_eq!(
        msg(obinary_file(None, 1, 16, 30, 3, 1, t, None, &[], &[])),
        "offset must by in [0,32 - hash_key_size]"
    );
    assert_eq!(
        msg(obinary_file(None, 1, 2, 3, 0, 30, t, None, &[0; 8], &[])),
        "shift must be in [1, 32 - hash_key_size]"
    );
    assert_eq!(
        msg(obinary_file(None, 1, 2, 3, 0, 2, t, None, &[0; 7], &[])),
        "not enough data to fill the hash table"
    );
    let mut short_vt = obinary_file(None, 1, 2, 3, 0, 2, t, None, &[0; 8], &[1, 2]);
    short_vt.pop();
    assert_eq!(msg(short_vt), "not enough data to fill the values table");
    // The unit word is past the end: NCBI's `readWord` returns its uninitialized word, read
    // here as 0.
    assert_eq!(
        msg(1u32.to_ne_bytes().to_vec()),
        "unit size must be in [1,16]"
    );
    assert_eq!(
        msg(obinary_file(None, 1, 2, 3, 0, 2, t, None, &[0; 8], &[])),
        "ok"
    );
}

#[test]
fn binary_header_errors_are_ncbis() {
    let w = |x: u32| x.to_ne_bytes().to_vec();
    let msg = |bytes: Vec<u8>| match read_counts(bytes) {
        Err(e) => e.to_string(),
        Ok(_) => "ok".to_string(),
    };
    assert_eq!(msg([w(0), w(2), w(0)].concat()), "wrong file size");
    let trailer = [w(1), w(2), w(3), w(4)].concat();
    assert_eq!(
        msg([w(0), w(0), trailer.clone()].concat()),
        "illegal unit size"
    );
    assert_eq!(
        msg([w(0), w(2), w(5), trailer.clone()].concat()),
        "wrong length"
    );
    let counts = read_counts([w(0), w(2), w(1), w(9), w(6), w(30), trailer].concat()).unwrap();
    // Units of 2 letters: AC = 0b0001 (its reverse complement GT = 0b1011 is larger).
    assert_eq!(counts.true_count(0b0001).unwrap(), 9);
    assert_eq!(counts.true_count(0b1011).unwrap(), 9);
    assert_eq!(counts.true_count(0b0110).unwrap(), 30);
    assert_eq!(counts.true_count(0b0000).unwrap(), 0);
    assert_eq!((counts.min_count(), counts.use_min_count()), (1, 1));
}

#[test]
fn obinary_lookup_follows_the_collision_field() {
    // Units of 2 letters, k = 3, roff = 1, bc = 2: key = bits 1-3, check = bit 0 (+ bits above
    // bit 3, none for canonical units below 16).
    let (k, roff, bc) = (3u32, 1u32, 2u32);
    let mut ht = vec![0u32; 8];
    let mut vt = Vec::new();
    // AC (0b0001): key 0, check 1, alone in its slot -> count in bits [bc, 24), check in 24-31.
    let (key, check) = hash_code(0b0001, k as u8, roff as u8);
    ht[key as usize] = 1 | (17 << bc) | ((check as u32) << 24);
    // AA (0b0000) and nothing else: key 0 is taken, so put AG (0b0010) and AT (0b0011) in a
    // collision list: key 1, checks 0 and 1, values (check << 9) | count.
    let (key_ag, check_ag) = hash_code(0b0010, k as u8, roff as u8);
    let (key_at, check_at) = hash_code(0b0011, k as u8, roff as u8);
    assert_eq!(key_ag, key_at);
    // Two collisions, list at index 0 of the values table.
    ht[key_ag as usize] = 2;
    vt.push(((check_ag as u16) << 9) | 5);
    vt.push(((check_at as u16) << 9) | 300);
    // A slot whose list runs past the values table.
    let (key_cc, _) = hash_code(0b0101, k as u8, roff as u8);
    ht[key_cc as usize] = 2 | (1 << bc);
    let counts = read_counts(obinary_file(
        None,
        1,
        2,
        k,
        roff,
        bc,
        [1, 2, 3, 400],
        None,
        &ht,
        &vt,
    ))
    .unwrap();
    assert_eq!(counts.true_count(0b0001).unwrap(), 17);
    assert_eq!(counts.true_count(0b1011).unwrap(), 17); // GT -> AC
    assert_eq!(counts.true_count(0b0000).unwrap(), 0); // same key as AC, other check
    assert_eq!(counts.true_count(0b0010).unwrap(), 5);
    assert_eq!(counts.true_count(0b0011).unwrap(), 300);
    assert_eq!(
        counts.true_count(0b0101).unwrap_err().to_string(),
        format!("bad index at key {key_cc} : {}", ht[key_cc as usize])
    );
}

#[test]
fn empty_and_short_sequences_have_no_masks() {
    let masker = SeqMasker::open(
        &fixture_dir().join("counts/mg1655_u9.obinary"),
        &MaskerParams::default(),
    )
    .unwrap();
    assert_eq!(masker.window_size(), 13);
    for seq in [&b""[..], b"A", b"ACGTACGTACGT", b"NNNNNNNNNNNNNNNNNNNN"] {
        assert_eq!(masker.mask(seq).unwrap(), Vec::new());
    }
}

#[test]
fn an_empty_cache_bit_array_fails_only_when_the_boost_reads_it() {
    // Units of 2 letters: NCBI writes a cache bit array of 4^2/32 = 0 words and keeps it in use
    // (`new Uint4[0]`). The manifest row `dm_u2` never scores below T_extend (500), so the boost
    // does not run and NCBI's intervals come out; here every unit is below T_low, so the first
    // window scores 1 < T_extend and the boost would read past the array.
    let file = obinary_file(
        None,
        2,
        2,
        3,
        0,
        2,
        [2, 1000, 2000, 3000],
        Some(&[]),
        &[0; 8],
        &[],
    );
    let counts = read_counts(file).unwrap();
    assert_eq!(counts.cache_bit_array(), Some(&[][..]));
    let masker = SeqMasker::new(Arc::new(counts), &MaskerParams::default()).unwrap();
    assert_eq!(
        masker.mask(b"ACGTACGTAC").unwrap_err(),
        WinMaskError::EmptyCacheBitArray { unit: 2 }
    );
    assert_eq!(masker.mask(b"ACGTA").unwrap(), Vec::new());
}

#[test]
fn a_short_metadata_block_is_an_unknown_format() {
    // CheckBin keeps the 3 of a short length word (skip 11) or reads the block as far as the
    // file goes; the skip of DiscoverStatType then fails: eUnknown.
    let kind = |bytes: Vec<u8>| discover(&mut Cursor::new(bytes)).unwrap().stat_type;
    let w = |x: u32| x.to_ne_bytes().to_vec();
    assert_eq!(kind(w(3)), StatType::Unknown);
    assert_eq!(
        kind([w(3), w(100), b"##a\0".to_vec()].concat()),
        StatType::Unknown
    );
    assert_eq!(
        kind([w(3), w(4), b"##a\0".to_vec()].concat()),
        StatType::Unknown
    );
    assert_eq!(
        kind([w(3), w(4), b"##a\0".to_vec(), w(0)].concat()),
        StatType::Binary
    );
    match read_counts(w(3)) {
        Err(WinMaskError::CreateFail { message }) => {
            assert_eq!(message, "unrecognized unit counts format")
        }
        Err(e) => panic!("{e}"),
        Ok(_) => panic!("must fail"),
    }
}

#[test]
fn a_file_ending_in_the_last_header_word_fails_at_the_hash_table() {
    let t = [1, 2, 3, 4];
    let header_only = |fmt: u32, cut: usize| {
        let mut b = obinary_file(None, fmt, 2, 3, 0, 2, t, None, &[], &[]);
        b.truncate(b.len() - cut);
        b
    };
    let err = |bytes: Vec<u8>| read_counts(bytes).err().expect("must fail");
    // Format 1 cut inside t_high, format 2 without its flag word: NCBI's next read is the
    // hash table.
    assert_eq!(
        err(header_only(1, 2)).to_string(),
        "not enough data to fill the hash table"
    );
    assert_eq!(
        err(header_only(2, 4)).to_string(),
        "not enough data to fill the hash table"
    );
    // Cut inside t_thres: the read of t_high finds the stream failed.
    assert_eq!(err(header_only(1, 6)).to_string(), "file too short");
    // A cache bit array cut short: NCBI's warning, then the hash table read fails.
    let mut short_cba = obinary_file(None, 2, 3, 5, 0, 2, t, Some(&[0, 0]), &[0; 32], &[]);
    short_cba.truncate(4 * 12);
    match err(short_cba) {
        WinMaskError::Format { site, warnings } => {
            assert_eq!(site, FormatSite::HashTable);
            assert_eq!(
                warnings,
                vec!["file read failed: bit array optimizations are not used.".to_string()]
            );
        }
        e => panic!("{e}"),
    }
}

#[test]
fn a_binary_data_part_of_minus_8_bytes_reads_the_trailer_after_the_header() {
    // A metadata block of 12 bytes, the two header words and two more words: datalen =
    // 36 - 8 - 16 - 20 = -8, which passes the unsigned `% 8`. t_low and t_extend are the two
    // words; the failed reads keep t_extend in `data` for t_threshold and t_high.
    let w = |x: u32| x.to_ne_bytes().to_vec();
    let md = b"##metadata!\0".to_vec();
    assert_eq!(md.len(), 12);
    let file = [w(3), w(12), md, w(0), w(2), w(7), w(40)].concat();
    let counts = read_counts(file).unwrap();
    assert_eq!(
        (
            counts.min_count(),
            counts.textend(),
            counts.threshold(),
            counts.max_count()
        ),
        (7, 40, 40, 40)
    );
    assert_eq!(counts.true_count(0b0001).unwrap(), 0);
}

#[test]
fn binary_lookup_probes_as_libstdcxx_lower_bound() {
    // Units out of order: AC (1), AT (3), CC (5) stored as 5, 1, 3. The halving loop finds
    // AT, misses AC and CC.
    let w = |x: u32| x.to_ne_bytes().to_vec();
    let file = [
        w(0),
        w(2),
        w(5),
        w(10),
        w(1),
        w(20),
        w(3),
        w(30),
        w(1),
        w(2),
        w(3),
        w(4),
    ]
    .concat();
    let counts = read_counts(file).unwrap();
    assert_eq!(counts.true_count(3).unwrap(), 30);
    assert_eq!(counts.true_count(1).unwrap(), 0);
    assert_eq!(counts.true_count(5).unwrap(), 0);
}

#[test]
fn a_bare_pct_line_is_the_factorys_substr_error() {
    let file = obinary_file(
        Some(&[
            "##windowmasker-statistics-format-version:obinary 1.1.0",
            "##pct:",
        ]),
        1,
        2,
        3,
        0,
        2,
        [1, 2, 3, 4],
        None,
        &[0; 8],
        &[],
    );
    assert_eq!(
        read_counts(file).err().expect("must fail").to_string(),
        "could not create a unit counts container: basic_string::substr: __pos (which is 7) > this->size() (which is 6)"
    );
    let file = obinary_file(
        Some(&["##pct: 0 0", "##pct:"]),
        1,
        2,
        3,
        0,
        2,
        [1, 2, 3, 4],
        None,
        &[0; 8],
        &[],
    );
    assert!(read_counts(file).is_ok());
}

// ---------------------------------------------------------------------------------------------
// The ascii and optimized ascii readers (line NP M3d). The texts are those of NCBI BLAST+ 2.17.0
// blastn `-window_masker_taxid` with such files (M/m3d/ws/run2.sh; the blastn fixtures
// `taxid_ascii_*` and `taxid_oascii_*` of tests/fixtures/winmask_blastn hold the same reports).

/// The report of an error as blastn prints it: the reader's exception, then the factory's rethrow.
fn reader_report(error: &WinMaskError) -> String {
    error
        .create_report()
        .unwrap_or_else(|| panic!("no report: {error}"))
}

fn rethrown_report(thrown: &str) -> String {
    let root = crate::ncbi_exception::NCBI_BUILD_ROOT;
    format!("NCBI C++ Exception:\n    T0 \"{root}/{thrown}\n    T0 \"{root}/src/algo/winmask/seq_masker_istat_factory.cpp\", line 358: Error: (Exception::creation failure) ncbi::CSeqMaskerIstatFactory::create() - could not create a unit counts container\n")
}

#[test]
fn the_text_formats_hold_the_counts_of_the_optimized_binary_file() {
    let open = |name: &str| {
        UnitCounts::open(
            &fixture_dir().join("counts").join(name),
            &MaskerParams::default(),
        )
        .unwrap()
    };
    let ob = open("mg1655_u9.obinary");
    let oascii_bytes = std::fs::read(fixture_dir().join("counts/mg1655_u9.oascii")).unwrap();
    let crlf = String::from_utf8(oascii_bytes)
        .unwrap()
        .replace('\n', "\r\n");
    let texts = [
        open("mg1655_u9.ascii"),
        open("mg1655_u9.oascii"),
        read_counts(crlf.into_bytes()).unwrap(),
    ];
    for text in &texts {
        assert_eq!(text.unit_size(), 9);
        assert_eq!(
            [
                text.min_count(),
                text.textend(),
                text.threshold(),
                text.max_count()
            ],
            [74, 147, 173, 210]
        );
        assert_eq!(
            [text.use_min_count(), text.use_max_count()],
            [ob.use_min_count(), ob.use_max_count()]
        );
        assert!(text.cache_bit_array().is_none() && text.warnings().is_empty());
        for unit in 0..1u32 << 18 {
            assert_eq!(
                text.true_count(unit).unwrap(),
                ob.true_count(unit).unwrap(),
                "{unit:x}"
            );
            assert_eq!(text.at(unit).unwrap(), ob.at(unit).unwrap(), "{unit:x}");
        }
    }
    assert!(ob.cache_bit_array().is_some());
}

#[test]
fn the_ascii_reader_reads_as_ncbis() {
    let md = "##windowmasker-statistics-format-version:ascii 1.1.0\n";
    let params = ">t_low 74\n>t_extend 147\n>t_threshold 173\n>t_high 210\n";
    let read = |text: String| read_counts(text.into_bytes());
    // Comments, blank lines, white space, upper case, parameters in bases 0, 8 and 16, a
    // parameter already set, one NCBI does not read, and a last line without its line feed.
    let counts = read(format!(
        "{md}#c\n\n+9\n  2\t76\n0x3 80\n#\n3FD00 100\n\n>t_low +74\n>t_extend 0223\n>t_threshold\t 0xad\n>t_threshold abc\n>t_x y\n>t_high 210"
    ))
    .unwrap();
    assert_eq!(counts.unit_size(), 9);
    assert_eq!(
        [
            counts.min_count(),
            counts.textend(),
            counts.threshold(),
            counts.max_count()
        ],
        [74, 147, 173, 210]
    );
    assert_eq!([counts.use_min_count(), counts.use_max_count()], [37, 210]);
    assert_eq!(
        [2, 3, 0x3fd00, 4].map(|unit| counts.true_count(unit).unwrap()),
        [76, 80, 100, 0]
    );
    // The parameters first: a unit below T_low is not kept, so its order is not checked.
    let counts = read(format!("{md}{params}9\n3 80\n2 5\n")).unwrap();
    assert_eq!(counts.true_count(2).unwrap(), 0);
    assert_eq!(counts.true_count(3).unwrap(), 80);
    // A second, higher t_low: NCBI's warning, and the higher value.
    let counts = read(format!("{md}9\n2 76\n{params}>t_low 80\n")).unwrap();
    assert_eq!(counts.min_count(), 80);
    assert_eq!(counts.warnings(), ["Requested value of t_low (74) is less than the one stored with the N-mer counts (80).The value 80 will be used."]);
    // No unit size line: the default 15; a unit size of 265 is 9 (`Uint1`).
    assert_eq!(read(format!("{md}{params}")).unwrap().unit_size(), 15);
    assert_eq!(read(format!("{md}265\n{params}")).unwrap().unit_size(), 9);
    // A threshold of 0 is read again from a later line ("0x" is 0).
    let counts = read(format!("{md}9\n>t_threshold 0x\n{params}")).unwrap();
    assert_eq!(counts.threshold(), 173);

    let ascii = "src/algo/winmask/seq_masker_istat_ascii.cpp\", line";
    let new = "ncbi::CSeqMaskerIstatAscii::CSeqMaskerIstatAscii()";
    let ncbistr = "src/corelib/ncbistr.cpp\", line";
    let cases = [
        (format!("{md}9\n2 76\n"), format!("{ascii} 169: Error: (Exception::bad parameter value) {new} - t_threhold t_extend t_high t_low ")),
        (format!("{md}9\n2 76\n>t_low 74\n>t_threshold 173\n"), format!("{ascii} 169: Error: (Exception::bad parameter value) {new} - t_extend t_high ")),
        (format!("{md}9\n>t_threshold\n"), format!("{ascii} 96: Error: (Exception::syntax error) {new} - at line 3")),
        (format!("{md}9\n>t_threshold \t \n"), format!("{ascii} 96: Error: (Exception::syntax error) {new} - at line 3")),
        (format!("{md}9\n2 76\n3\n"), format!("{ascii} 137: Error: (Exception::syntax error) {new} - at line 4")),
        (format!("{md}9\n2 76\n3 \n"), format!("{ascii} 137: Error: (Exception::syntax error) {new} - at line 4")),
        (format!("{md}9\n3 80\n3 76\n"), "src/algo/winmask/seq_masker_uset_simple.cpp\", line 62: Error: (Exception::bad unit order) ncbi::CSeqMaskerUsetSimple::add_info() - last unit: 3 ; adding 3".to_string()),
        (format!("{md}9\n1f 80\n1e 76\n"), "src/algo/winmask/seq_masker_uset_simple.cpp\", line 62: Error: (Exception::bad unit order) ncbi::CSeqMaskerUsetSimple::add_info() - last unit: 1f ; adding 1e".to_string()),
        (format!("{md}9x\n{params}"), format!("{ncbistr} 908: Error: (CStringException::eConvert) ncbi::NStr::StringToUInt8() - Cannot convert string '9x' to Uint8 (m_Pos = 1)")),
        (format!("{md}9\r\n{params}"), format!("{ncbistr} 908: Error: (CStringException::eConvert) ncbi::NStr::StringToUInt8() - Cannot convert string '9\\r' to Uint8 (m_Pos = 1)")),
        (format!("{md}9\n2 -5\n"), format!("{ncbistr} 908: Error: (CStringException::eConvert) ncbi::NStr::StringToUInt8() - Cannot convert string '-5' to Uint8 (m_Pos = 0)")),
        (format!("{md}9\n2 76 1\n"), format!("{ncbistr} 908: Error: (CStringException::eConvert) ncbi::NStr::StringToUInt8() - Cannot convert string '76 1' to Uint8 (m_Pos = 2)")),
        (format!("{md}9\n1g 76\n"), format!("{ncbistr} 980: Error: (CStringException::eConvert) ncbi::NStr::StringToUInt8() - Cannot convert string '1g' to Uint8 (m_Pos = 1)")),
        (format!("{md}9\n1ffffffff 76\n"), format!("{ncbistr} 652: Error: (CStringException::eConvert) ncbi::NStr::StringToUInt() - Cannot convert string '1ffffffff' to unsigned int, overflow (m_Pos = 0)")),
        (format!("{md}9\n>t_low abc\n"), format!("{ncbistr} 973: Error: (CStringException::eConvert) ncbi::NStr::StringToUInt8() - Cannot convert string 'abc' to Uint8 (m_Pos = 0)")),
        // The version of the statistics algorithm in the metadata (`ExtractStatAlgoVersion`).
        (format!("{md}##windowmasker-statistics-algorithm-version:99999999999.0.0\n9\n{params}"), format!("{ncbistr} 640: Error: (CStringException::eConvert) ncbi::NStr::StringToInt() - Cannot convert string '99999999999' to int, overflow (m_Pos = 0)")),
        (format!("{md}##windowmasker-statistics-algorithm-version:.0.0\n9\n{params}"), format!("{ncbistr} 861: Error: (CStringException::eConvert) ncbi::NStr::StringToInt8() - Cannot convert string '' to Int8 (m_Pos = 0)")),
    ];
    for (text, thrown) in cases {
        let error = read(text.clone())
            .err()
            .unwrap_or_else(|| panic!("{text:?} was read"));
        assert_eq!(reader_report(&error), rethrown_report(&thrown), "{text:?}");
    }
    // A version that is not all digits is skipped.
    assert!(read(format!(
        "{md}##windowmasker-statistics-algorithm-version:1.x.99999999999\n9\n{params}"
    ))
    .is_ok());
    // A warning posted before the error goes with it.
    let error = read(format!("{md}9\n>t_low 5\n>t_low 6\n2 x\n"))
        .err()
        .unwrap();
    assert_eq!(error.warnings().len(), 1);
    // Every `##pct:` line of an ascii file is read, the second one too.
    let error = read(format!("{md}##pct: 0 0\n##pct:\n9\n{params}"))
        .err()
        .unwrap();
    assert!(
        error
            .to_string()
            .ends_with("(which is 7) > this->size() (which is 6)"),
        "{error}"
    );
}

#[test]
fn the_optimized_ascii_reader_reads_as_ncbis() {
    let md = "##windowmasker-statistics-format-version:oascii 1.1.0\n";
    // Units of 2 letters, a hash key of 3 bits, no offset, a 1-bit collision count, 2 values.
    let table = "0\n0\n0\n0\n0\n0\n0\n0\n5\n6\n";
    let read = |text: String| read_counts(text.into_bytes());
    let counts = read(format!("{md}AAAA\n2\n2 3 0 1\n4\n5\n6\n7\n{table}")).unwrap();
    assert_eq!(counts.unit_size(), 2);
    assert_eq!(
        [
            counts.min_count(),
            counts.textend(),
            counts.threshold(),
            counts.max_count()
        ],
        [4, 5, 6, 7]
    );
    // atoi: white space, a sign, the digits up to another byte; `Uint1` and `Uint2` truncate.
    let counts = read(format!(
        "{md}AAAAxyz\n 258x\n2 3 0 1 9 x\n-1\n+5\n6abc\n99999999999\n{table}"
    ))
    .unwrap();
    assert_eq!(counts.unit_size(), 2);
    assert_eq!(
        [
            counts.min_count(),
            counts.textend(),
            counts.threshold(),
            counts.max_count()
        ],
        [u32::MAX, 5, 6, 0x4876_E7FF]
    );
    let oascii = "src/algo/winmask/seq_masker_istat_oascii.cpp\", line";
    let new = "ncbi::CSeqMaskerIstatOAscii::CSeqMaskerIstatOAscii()";
    let header = |unit: &str, hash: &str| format!("{md}AAAA\n{unit}\n{hash}\n4\n5\n6\n7\n{table}");
    let cases = [
        (header("0", "2 3 0 1"), format!("{oascii} 101: Error: (Exception::bad parameter) {new} - unit size must be in [1,16]")),
        (header("abc", "2 3 0 1"), format!("{oascii} 101: Error: (Exception::bad parameter) {new} - unit size must be in [1,16]")),
        (header("17", "2 3 0 1"), format!("{oascii} 101: Error: (Exception::bad parameter) {new} - unit size must be in [1,16]")),
        (header("2", "2 4 0 1"), format!("{oascii} 114: Error: (Exception::bad hash parameter) {new} - hash key size must be in [1,2*unit_size - 1]")),
        (header("2", "2 0 0 1"), format!("{oascii} 114: Error: (Exception::bad hash parameter) {new} - hash key size must be in [1,2*unit_size - 1]")),
        (header("2", "2 3 30 1"), format!("{oascii} 119: Error: (Exception::bad hash parameter) {new} - offset must by in [0,32 - hash_key_size]")),
        (header("2", "2 3 -1 1"), format!("{oascii} 119: Error: (Exception::bad hash parameter) {new} - offset must by in [0,32 - hash_key_size]")),
        (header("2", "2 3 0 0"), format!("{oascii} 124: Error: (Exception::bad hash parameter) {new} - shift must be in [1, 32 - hash_key_size]")),
        (header("2", "2 3 0 30"), format!("{oascii} 124: Error: (Exception::bad hash parameter) {new} - shift must be in [1, 32 - hash_key_size]")),
        (header("2", "2 3 0 4294967296"), format!("{oascii} 124: Error: (Exception::bad hash parameter) {new} - shift must be in [1, 32 - hash_key_size]")),
        // A value that `>>` fails to read is set (0, or UINT_MAX past it): NCBI's checks.
        (header("2", "2 x 0 1"), format!("{oascii} 114: Error: (Exception::bad hash parameter) {new} - hash key size must be in [1,2*unit_size - 1]")),
        (header("2", "2 3 99 x"), format!("{oascii} 119: Error: (Exception::bad hash parameter) {new} - offset must by in [0,32 - hash_key_size]")),
        (header("2", "2 3 0 x"), format!("{oascii} 124: Error: (Exception::bad hash parameter) {new} - shift must be in [1, 32 - hash_key_size]")),
        (header("2", "2 99"), format!("{oascii} 114: Error: (Exception::bad hash parameter) {new} - hash key size must be in [1,2*unit_size - 1]")),
        (format!("{md}AAAA\n2\n2 3 0 1\n4\n5\n6\n"), format!("{oascii} 145: Error: (Exception::format error) {new} - file too short")),
        (format!("{md}AAAA\n2\n2 3 0 1\n4\n5\n6\n7\n0\n0\n"), format!("{oascii} 176: Error: (Exception::format error) {new} - not enough lines to fill the hash table")),
        (format!("{md}AAAA\n2\n2 3 0 1\n4\n5\n6\n7\n0\n0\n0\n0\n0\n0\n0\n0\n5\n"), format!("{oascii} 192: Error: (Exception::format error) {new} - not enough lines to fill the values table")),
        // Sizes far beyond the file: no allocation of NCBI's size.
        (format!("{md}AAAA\n16\n4000000000 31 0 1\n4\n5\n6\n7\n1\n"), format!("{oascii} 176: Error: (Exception::format error) {new} - not enough lines to fill the hash table")),
        (format!("{md}AAAA\n2\n4000000000 3 0 1\n4\n5\n6\n7\n{table}"), format!("{oascii} 192: Error: (Exception::format error) {new} - not enough lines to fill the values table")),
    ];
    for (text, thrown) in cases {
        let error = read(text.clone())
            .err()
            .unwrap_or_else(|| panic!("{text:?} was read"));
        assert_eq!(reader_report(&error), rethrown_report(&thrown), "{text:?}");
    }
    // A value that `>>` does not reach and the checks read: NCBI's is uninitialized (its errors
    // change from one load to the next), refused.
    for (hash, field) in [
        ("2 3 0", "bc"),
        ("2 3 x 1", "bc"),
        ("2 3", "roff"),
        ("x 3 0 1", "k"),
        ("2", "k"),
        ("", "k"),
        (" ", "k"),
    ] {
        let error = read(header("2", hash)).err().unwrap();
        assert_eq!(
            error,
            WinMaskError::UnsetHashParameter { field },
            "{hash:?}"
        );
        assert!(error.create_report().is_none());
    }
}

#[test]
fn hash_parameters_are_read_as_libstdcxx_reads_unsigned_ints() {
    use super::istat::extract_uints;
    let read = |line: &[u8]| extract_uints::<4>(line);
    assert_eq!(read(b"3694 17 0 2"), ([3694, 17, 0, 2], 4));
    assert_eq!(read(b"\t3694\x0b17 0 2 99 x"), ([3694, 17, 0, 2], 4));
    assert_eq!(read(b"+1 -1 -4294967295 2x"), ([1, u32::MAX, 1, 2], 4));
    // A failure sets 0 (no digit) or UINT_MAX (overflow); nothing is read after it.
    assert_eq!(read(b"1 x 3 4"), ([1, 0, 0, 0], 2));
    assert_eq!(read(b"1 4294967296 3 4"), ([1, u32::MAX, 0, 0], 2));
    assert_eq!(read(b"1 - 3 4"), ([1, 0, 0, 0], 2));
    // Nothing is read after the end of the line.
    assert_eq!(read(b"1 2"), ([1, 2, 0, 0], 2));
    assert_eq!(read(b"1 2\r"), ([1, 2, 0, 0], 2));
    assert_eq!(read(b""), ([0; 4], 0));
    assert_eq!(read(b" "), ([0; 4], 0));
}

#[test]
fn the_masker_checks_the_unit_size_of_an_ascii_file() {
    let counts = |unit: &str| {
        Arc::new(
            read_counts(
                format!("{unit}\n>t_low 74\n>t_extend 147\n>t_threshold 173\n>t_high 210\n")
                    .into_bytes(),
            )
            .unwrap(),
        )
    };
    let params = MaskerParams::default();
    // 0 (and 256): a window of 4 letters, the units the last letter.
    assert_eq!(
        SeqMasker::new(counts("0"), &params).unwrap().window_size(),
        4
    );
    assert_eq!(
        SeqMasker::new(counts("16"), &params).unwrap().window_size(),
        20
    );
    // 252 to 255: the window (`Uint1`) wraps, NCBI's error (not rethrown by the factory).
    let error = SeqMasker::new(counts("252"), &params).err().unwrap();
    let root = crate::ncbi_exception::NCBI_BUILD_ROOT;
    assert_eq!(
        reader_report(&error),
        format!("NCBI C++ Exception:\n    T0 \"{root}/src/algo/winmask/seq_masker.cpp\", line 126: Error: (CSeqMaskerException::validation error) ncbi::CSeqMasker::CSeqMasker() - window size (0) must be greater or equal to unit size (252)\n")
    );
    assert_eq!(
        SeqMasker::new(counts("255"), &params).err().unwrap(),
        WinMaskError::WindowSize {
            window: 3,
            unit: 255,
            warnings: Vec::new()
        }
    );
    // The reader's warnings go with the error (review 1 F2).
    let warned = read_counts(
        b"255\n>t_low 74\n>t_extend 147\n>t_threshold 173\n>t_high 210\n>t_low 80\n".to_vec(),
    )
    .unwrap();
    let error = SeqMasker::new(Arc::new(warned), &params).err().unwrap();
    assert_eq!(error.warnings().len(), 1);
    // 17 to 251: refused (NCBI's reverse complement of the units is undefined).
    for unit in ["17", "251"] {
        let error = SeqMasker::new(counts(unit), &params).err().unwrap();
        assert!(matches!(error, WinMaskError::UnitSize { .. }), "{error}");
        assert!(error.create_report().is_none());
    }
}

// The metadata entries are NCBI's `std::string`s: bytes. `ExtractStatAlgoVersion` takes the name
// from the third byte on, whatever the first two are (review 1 F3).
#[test]
fn the_version_of_binary_metadata_is_read_from_its_bytes() {
    let data = std::fs::read(fixture_dir().join("counts/mg1655_u9.obinary")).unwrap();
    let md_len = u32::from_le_bytes(data[4..8].try_into().unwrap()) as usize;
    let body = &data[8 + md_len..];
    let with = |entry: &[u8]| {
        let block = [entry, b"\0"].concat();
        let len = (block.len() as u32).to_le_bytes();
        read_counts([&3u32.to_le_bytes()[..], &len, &block, body].concat())
    };
    let error = with(b"\xff\xffwindowmasker-statistics-algorithm-version:99999999999.0.0")
        .err()
        .unwrap();
    assert!(
        error
            .to_string()
            .ends_with("Cannot convert string '99999999999' to int, overflow"),
        "{error}"
    );
    assert!(with(b"\xff\xffwindowmasker-statistics-algorithm-version:1.0.0").is_ok());
}
