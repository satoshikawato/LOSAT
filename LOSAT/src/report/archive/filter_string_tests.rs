//! `BlastFilteringOptionsFromString` of the forms that NCBI BLAST+ 2.17.0 gives `blast_formatter
//! -outfmt 11` (line RG-F3 review R6 M1: `remaining-gaps-20261112/F/f3/tmp/r6/loc/`).

use super::*;

fn seg(window: Option<i32>, locut: Option<f64>, hicut: Option<f64>) -> Option<Seg> {
    Some(Seg {
        window,
        locut,
        hicut,
    })
}

#[test]
fn the_commands_and_their_defaults() {
    let parse = |text: &str, blastn| filtering_options_from_string(text.as_bytes(), blastn);
    assert_eq!(parse("F", true), Some(FilterOptions::default()));
    assert_eq!(parse("f", false), Some(FilterOptions::default()));
    let lm = parse("L;m;", true).unwrap();
    assert_eq!(
        (lm.mask_at_hash, lm.dust, lm.seg),
        (true, Some(DEFAULT_DUST), None)
    );
    let lm = parse("m L", false).unwrap();
    assert_eq!(
        (lm.mask_at_hash, lm.dust, lm.seg),
        (true, None, Some(DEFAULT_SEG))
    );
    assert_eq!(
        parse("D 10 20 1;", true).unwrap().dust,
        Some([Some(10), Some(20), Some(1)])
    );
    assert_eq!(parse("D;", true).unwrap().dust, Some(DEFAULT_DUST));
    assert_eq!(parse("S;", false).unwrap().seg, Some(DEFAULT_SEG));
    let r = parse("R -d myrep;m;", true).unwrap();
    assert_eq!(
        (r.repeat_db, r.mask_at_hash),
        (Some(b"myrep".to_vec()), true)
    );
    assert_eq!(
        parse("R;", true).unwrap().repeat_db,
        Some(DEFAULT_REPEAT_FILTER_DB.to_vec())
    );
    assert_eq!(
        parse("R -x", true).unwrap().repeat_db,
        Some(DEFAULT_REPEAT_FILTER_DB.to_vec())
    );
    // The window masker's options leave the others as they are.
    let w = parse("L;W -t 9606;", true).unwrap();
    assert_eq!(
        (w.dust, w.seg, w.mask_at_hash),
        (Some(DEFAULT_DUST), None, false)
    );
    // Two values: NCBI's parse error.
    assert_eq!(parse("L;D 10 20;m;", true), None);
    assert_eq!(parse("S 10 1.8;", false), None);
    assert_eq!(parse("D 1 2 3 4", true), None);
}

#[test]
fn the_values_are_sscanf_s() {
    let parse = |text: &str, blastn| filtering_options_from_string(text.as_bytes(), blastn);
    // `%ld` beyond `long` is `LONG_MAX`, whose low 32 bits are -1; `LONG_MIN`'s are 0.
    assert_eq!(
        parse("D 99999999999999999999 64 1;", true).unwrap().dust,
        Some([Some(-1), Some(64), Some(1)])
    );
    assert_eq!(
        parse("S -99999999999999999999 1.8 2.1;", false)
            .unwrap()
            .seg,
        seg(Some(0), Some(1.8), Some(2.1))
    );
    assert_eq!(
        parse("S 10 0x1p1 2.5;", false).unwrap().seg,
        seg(Some(10), Some(2.0), Some(2.5))
    );
    let nan = parse("S 10 nan(1) 2.5;", false).unwrap().seg.unwrap();
    assert_eq!(nan.locut.map(f64::to_bits), Some(f64::NAN.to_bits()));
    assert_eq!(
        parse("S 0x10 inf 12x;", false).unwrap().seg,
        seg(Some(0), Some(f64::INFINITY), Some(12.0))
    );
    // A token that `sscanf` does not convert: no value.
    assert_eq!(
        parse("S 10 abc 2.5;", false).unwrap().seg,
        seg(Some(10), None, Some(2.5))
    );
    assert_eq!(
        parse("D abc 64 1;", true).unwrap().dust,
        Some([None, Some(64), Some(1)])
    );
    assert_eq!(
        parse("S x 0x 0xg;", false).unwrap().seg,
        seg(None, None, None)
    );
}

#[test]
fn the_buffers_end_where_ncbi_s_do() {
    // `s_LoadOptionsToBuffer` reads 128 bytes of a command; what follows is read as commands.
    let long = format!("R {}m", "x".repeat(127));
    let options = filtering_options_from_string(long.as_bytes(), true).unwrap();
    assert!(options.mask_at_hash, "{options:?}");
    assert_eq!(options.repeat_db, Some(DEFAULT_REPEAT_FILTER_DB.to_vec()));
    // One token of 128 bytes: the token's end is the 129th step, so no token is read (`arg` 0).
    let one = "9".repeat(128);
    assert_eq!(parse_dust_options(one.as_bytes()), Some([Some(-1); 3]));
    // Three tokens whose third ends past the 128th step: two are read (NCBI's parse error).
    let late = format!("1 2 {}", "3".repeat(124));
    assert_eq!(late.len(), 128);
    assert_eq!(parse_dust_options(late.as_bytes()), None);
    // ... and a fourth token that the steps do not reach the end of: three are read.
    let fourth = format!("1 2 3 {}", "4".repeat(122));
    assert_eq!(
        parse_dust_options(fourth.as_bytes()),
        Some([Some(1), Some(2), Some(3)])
    );
    // Blanks before and after a command are not its bytes.
    assert_eq!(
        load_options_to_buffer(b"   10 20 1   ;m"),
        (b"10 20 1".to_vec(), &b"m"[..])
    );
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:130-171 (cited at
// `parse_window_masker_options`). Not a port of an NCBI unit case.
#[test]
fn the_window_masker_command_keeps_its_database_or_taxid() {
    let w = |text: &str| {
        filtering_options_from_string(text.as_bytes(), true)
            .unwrap()
            .window_masker
    };
    let masker = |database: Option<&[u8]>, taxid| {
        Some(WindowMaskerFilter {
            database: database.map(<[u8]>::to_vec),
            taxid,
        })
    };
    assert_eq!(w("L;m;"), None);
    assert_eq!(w("W;"), masker(None, 0));
    assert_eq!(w("W -t 9606;m"), masker(None, 9606));
    assert_eq!(w("W -t\t 10090x"), masker(None, 10090));
    assert_eq!(w("W -t abc"), masker(None, 0));
    // The database is cut at a blank or a tab, and wins over a taxid.
    assert_eq!(w("W -d my wm -t 9606;"), masker(Some(b"my"), 0));
    assert_eq!(
        w("W -t 9606 -d\twm.obinary;"),
        masker(Some(b"wm.obinary"), 0)
    );
    // The last command's options.
    assert_eq!(w("W -d one;W -t 7;"), masker(None, 7));
    // Not BLASTN's: the window masker's options are read all the same.
    assert_eq!(
        filtering_options_from_string(b"W -t 9606", false)
            .unwrap()
            .window_masker,
        masker(None, 9606)
    );
}
