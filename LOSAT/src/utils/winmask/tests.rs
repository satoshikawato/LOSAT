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
    assert_eq!(rows, 21);
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
    assert_eq!(err, WinMaskError::WindowSize { window: 5, unit: 9 });
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
        "##windowmasker-statistics-format-version:obinary 1.1.0"
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
    let mut w = |b: &mut Vec<u8>, x: u32| b.extend_from_slice(&x.to_ne_bytes());
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
    assert_eq!(found.metadata, vec!["##a".to_string()]);
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
        Err(WinMaskError::UnsupportedFormat(StatType::OAscii)) => {}
        Err(e) => panic!("{e}"),
        Ok(_) => panic!("oascii must be refused"),
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
