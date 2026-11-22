//! Tests of the Seq-table access (`table`): the values of a column at a row as NCBI's
//! `CSeqTable_column`, `CSeqTable_multi_data`, `CSeqTable_single_data` and
//! `CSeqTable_sparse_index` give them (line RG D3-1 `om/SPEC.md` §7).

use super::super::view::V;
use super::super::{Indexer, Stop};
use super::{Caches, Column, Want};
use crate::report::serial::asn_input::Names;
use crate::report::serial::asn_types::{read_document, table::named, Value};

fn column_value(text: &str) -> Value {
    let (value, posts) = read_document(
        format!("SeqTable-column ::= {{ header {{ field-id 3 }}, {text} }}").as_bytes(),
        named::SEQTABLE_COLUMN,
        Names::PSSM,
    );
    assert!(posts.is_empty());
    value.unwrap_or_else(|error| panic!("{text}: {error:?}"))
}

fn indexer() -> Indexer<'static> {
    Indexer {
        program: "psiblast",
        posts: crate::report::diag_log::Posts::new(),
        muted: false,
        caches: Caches::default(),
        handles: Default::default(),
    }
}

/// `TryGetValue(row, v)` of the column `text` at each of `rows`, as text (`-` for none, `E:` and
/// the message for an exception, `S:` and `what()` for a `std::exception`, `R` for a refusal).
fn values(text: &str, want: Want, rows: std::ops::Range<u64>) -> Vec<String> {
    let value = column_value(text);
    let column = Column(V::new(named::SEQTABLE_COLUMN, &value));
    let mut ix = indexer();
    rows.map(|row| match ix.column_try_get(column, row, want) {
        Ok(Some(v)) => v.to_string(),
        Ok(None) => "-".to_string(),
        Err(Stop::Exception(e)) => {
            format!("E:{}", String::from_utf8(e.get_msg().to_vec()).unwrap())
        }
        Err(Stop::Std(what)) => format!("S:{what}"),
        Err(Stop::Refuse(r)) => format!("R:{r:?}"),
        Err(other) => format!("{other:?}"),
    })
    .collect()
}

#[test]
fn data_then_default() {
    assert_eq!(
        values("data int { 4, -2 }", Want::Int4, 0..3),
        ["4", "-2", "-"]
    );
    assert_eq!(
        values("data int { 4 }, default int 7", Want::Int4, 0..2),
        ["4", "7"]
    );
    // `bit`: the most significant bit is row 0; rows past the bytes have no value.
    assert_eq!(
        values("data bit 'A0'H", Want::Bool, 0..9),
        ["1", "0", "1", "0", "0", "0", "0", "0", "-"]
    );
    // `int1` is signed bytes.
    assert_eq!(values("data int1 'FF7F'H", Want::Int4, 0..2), ["-1", "127"]);
    // A bool is 0 or 1; another value is too big (E14).
    assert_eq!(
        values("data int { 1, 2 }", Want::Bool, 0..2),
        [
            "1",
            "E:CSeqTable_multi_data::TryGetBool(): value is too big for requested type: 2"
        ]
    );
    // A type without integers (E13), a single value of another type (E15, E16).
    assert_eq!(
        values("data string { \"a\" }", Want::Int8, 0..1),
        ["E:CSeqTable_multi_data::TryGetInt8(): value cannot be converted to requested type"]
    );
    assert_eq!(
        values("default real { 1, 10, 0 }", Want::Int4, 0..1),
        ["E:CSeqTable_single_data::GetValue(int&): real value cannot be converted to requested type"]
    );
    assert_eq!(
        values("default int8 4294967296", Want::Int4, 0..1),
        ["E:CSeqTable_single_data::GetValue(int&): value is too big for requested type: 4294967296"]
    );
    assert_eq!(
        values("default int8 4294967296", Want::Int8, 0..1),
        ["4294967296"]
    );
}

// NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_sparse_index.cpp:533-580 (cited at
// `Indexer::sparse_index_at`)
#[test]
fn sparse_indexes() {
    // `indexes`: NCBI's `lower_bound` (`utils::lower_bound`), which on the unsorted { 5, 3, 9 }
    // finds 9 but neither 3 nor 5.
    assert_eq!(
        values(
            "data int { 10, 11, 12 }, sparse indexes { 1, 4, 6 }",
            Want::Int4,
            0..7
        ),
        ["-", "10", "-", "-", "11", "-", "12"]
    );
    assert_eq!(
        values(
            "data int { 10, 11, 12 }, sparse indexes { 5, 3, 9 }",
            Want::Int4,
            3..10
        ),
        ["-", "-", "-", "-", "-", "-", "12"]
    );
    // On { 0, 0, 1, 2, 0 } (observed with `observe_lower_bound one`): row 1 is index 2, and row
    // 2 has no value, where `slice::partition_point` would find index 3.
    assert_eq!(
        values(
            "data int { 10, 11, 12, 13, 14 }, sparse indexes { 0, 0, 1, 2, 0 }",
            Want::Int4,
            0..3
        ),
        ["10", "12", "-"]
    );
    // `indexes-delta`: the rows 2, 3 and 7.
    assert_eq!(
        values(
            "data int { 10, 11, 12 }, sparse indexes-delta { 2, 1, 4 }",
            Want::Int4,
            0..8
        ),
        ["-", "-", "10", "11", "-", "-", "-", "12"]
    );
    // `bit-set`: the set bits before the row count the data's index; a skipped row takes
    // `sparse-other`.
    assert_eq!(
        values(
            "data int { 10, 11, 12 }, sparse bit-set '4180'H, sparse-other int 0",
            Want::Int4,
            0..10
        ),
        ["0", "10", "0", "0", "0", "0", "0", "11", "12", "0"]
    );
    // `bit-set-bvector`: BitMagic data (`bitmagic`; the bits 1 and 5 of a vector of size 8, line RG
    // D3-1 `bm/` case s03), decoded at the first index: a row is in the set when it is below the
    // vector's size and its bit is set, and the set bits before it count its index. An unknown
    // token is BitMagic's exception (case m11); data that ends before its first block is refused
    // (a read past the data, cases m02, m09).
    let bits_1_5 = "'02018000000100020005080000004344720F0009'H";
    assert_eq!(
        values(
            &format!("data int {{ 10, 11 }}, sparse bit-set-bvector {{ size 8, data {bits_1_5} }}"),
            Want::Int4,
            0..9
        ),
        ["-", "10", "-", "-", "-", "11", "-", "-", "-"]
    );
    assert_eq!(
        values(
            "data int { 1 }, sparse bit-set-bvector { size 4, data '02018000000100020005040000000C'H }",
            Want::Int4,
            0..2
        ),
        [
            "S:BM::Invalid serialization format",
            "S:BM::Invalid serialization format"
        ]
    );
    assert_eq!(
        values(
            "data int { 1 }, sparse bit-set-bvector { size 8, data '00'H }",
            Want::Int4,
            0..1
        ),
        ["R:BitVector"]
    );
}

// NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_multi_data.cpp:474-485, 536-538
// (cited at `Indexer::x_try_get_int8`)
#[test]
fn bit_bvector_values() {
    // A row below the decoded vector's size has its bit; a row past it has no value. The size is
    // the larger of the `size` member and the stream's header (cases s01, s02).
    let bits_1_5 = "'02018000000100020005080000004344720F0009'H";
    assert_eq!(
        values(
            &format!("data bit-bvector {{ size 8, data {bits_1_5} }}"),
            Want::Bool,
            0..9
        ),
        ["0", "1", "0", "0", "0", "1", "0", "0", "-"]
    );
    assert_eq!(
        values(
            &format!("data bit-bvector {{ size 3, data {bits_1_5} }}"),
            Want::Int8,
            4..9
        ),
        ["0", "1", "0", "0", "-"]
    );
    assert_eq!(
        values(
            "data bit-bvector { size 3, data '02018000000100020005030000004344720F0009'H }",
            Want::Int4,
            0..4
        ),
        ["0", "1", "0", "-"]
    );
}

/// The rows of a `bm/` probe that NCBI's index reads.
enum Rows<'a> {
    Runs(&'a [(u64, u64)]),
    Every { first: u64, step: u64, count: u64 },
}

impl Rows<'_> {
    fn has(&self, row: u64) -> bool {
        match *self {
            Rows::Runs(runs) => runs.iter().any(|&(a, b)| (a..=b).contains(&row)),
            Rows::Every { first, step, count } => {
                row >= first && (row - first).is_multiple_of(step) && (row - first) / step < count
            }
        }
    }

    fn ends(&self) -> Vec<u64> {
        match *self {
            Rows::Runs(runs) => runs.iter().flat_map(|&(a, b)| [a, b]).collect(),
            Rows::Every { first, step, count } => vec![first, first + step * (count - 1)],
        }
    }
}

/// The rows of a table of `rows` rows to look at: all of a table of up to 200,000 rows; of a
/// larger one, the first 140,000, the ends of BitMagic's blocks, the rows around those that
/// NCBI's index reads, and the last two.
fn probe_rows(rows: u64, read: &Rows<'_>) -> std::collections::BTreeSet<u64> {
    if rows <= 200_000 {
        return (0..rows).collect();
    }
    let mut out: std::collections::BTreeSet<u64> = (0..140_000).collect();
    for block in 1..=rows >> 16 {
        out.extend([(block << 16) - 1, block << 16]);
    }
    for row in read.ends() {
        out.extend([row.saturating_sub(1), row, row + 1]);
    }
    out.extend([rows - 2, rows - 1]);
    out.retain(|&row| row < rows);
    out
}

/// Whether the rows `rows` of a feature table whose `disabled` column is `column` are read
/// (`RowIsDisabled(row)` false), against `read`.
fn check_disabled_rows(case: &str, column: &str, rows: u64, read: &Rows<'_>) {
    let value = column_value(column);
    let column = Column(V::new(named::SEQTABLE_COLUMN, &value));
    let mut ix = indexer();
    for row in probe_rows(rows, read) {
        let disabled = match ix.column_try_get(column, row, Want::Bool) {
            Ok(value) => value == Some(1),
            Err(stop) => panic!("{case}: row {row}: {stop:?}"),
        };
        assert_eq!(!disabled, read.has(row), "{case}: row {row}");
    }
}

// The BitMagic probes w01-w23, s01-s17 and p04 whose data is short (PSSMs whose query's Seq-tables
// hold bit vectors, run by NCBI BLAST+ 2.17.0 psiblast and tblastn; the rows are those of NCBI's
// `Empty region ... Seq-table row <row>` posts, which LOSAT refuses for the memory address, taken
// from NCBI's standard error case by case by a script outside this repository). s05 is a read error
// of the `size`; the cases of longer data are those of `bm_probes_of_a_directory`.
const PROBES: &[(&str, u64, &str, Rows<'static>)] = &[
    ("w01_1bit", 16, "data bit-bvector { size 16, data '020180000001000200051000000013030009'H }", Rows::Runs(&[(0, 2), (4, 15)])),
    ("w02_gap_v3s", 100, "data bit-bvector { size 100, data '02018000000100020005640000003E6608FE3F0900000009'H }", Rows::Runs(&[(1, 4), (6, 6), (8, 99)])),
    ("w03_gap_egamma_v3", 64, "data bit-bvector { size 64, data '02018000000100020005400000004332260E0009'H }", Rows::Runs(&[(0, 9), (21, 63)])),
    ("w04_gap_bienc_v3", 1300, "data bit-bvector { size 65536, data '02018000000100020005000001003D30D402E8FF25263209'H }", Rows::Runs(&[(0, 99), (201, 299), (401, 499), (601, 699), (801, 899), (1001, 1099), (1201, 1299)])),
    ("w08_arr_v3", 65536, "sparse bit-set-bvector { size 65536, data '02018000000100020005000001003FB10400C80000FF3F0000000009'H }, default bit FALSE, sparse-other bit TRUE", Rows::Every { first: 0, step: 4, count: 16384 }),
    ("w09_gap_v3_holes", 65536, "data bit-bvector { size 65536, data '02018000000100020005000001003D8012523501089199B09A3634E96D5EE995F5C4461F371E8C97783F1CD121C7F9ADDC2FAD62F02A0009'H }", Rows::Runs(&[(4914, 4914), (6328, 6328), (7602, 7602), (7747, 7747), (8108, 8108), (9156, 9156), (9494, 9494), (11265, 11265), (11889, 11889), (12337, 12337), (16226, 16226), (28140, 28140), (29260, 29260), (31544, 31544), (47931, 47931), (51750, 51750), (51993, 51993), (54810, 54810), (55642, 55642), (56838, 56838)])),
    ("w13_1one", 70000, "data bit-bvector { size 70000, data '020180000001000200057011010002430200E02D0200000009'H }", Rows::Runs(&[(65536, 65536)])),
    ("w14_8one", 327700, "data bit-bvector { size 327700, data '02018000000100020005140005000405439268020009'H }", Rows::Runs(&[(327680, 327680)])),
    ("w15_16one", 16777226, "data bit-bvector { size 16777226, data '020180000001000200050A0000010600014392C4000009'H }", Rows::Runs(&[(16777216, 16777216)])),
    ("w16_shortzero", 196628, "sparse bit-set-bvector { size 196628, data '02018000000100020005140003001305008213070009'H }, default bit FALSE, sparse-other bit TRUE", Rows::Runs(&[(5, 5), (196615, 196615)])),
    ("w17_1zero", 131092, "sparse bit-set-bvector { size 131092, data '02018000000100020005140002001305000113070009'H }, default bit FALSE, sparse-other bit TRUE", Rows::Runs(&[(5, 5), (131079, 131079)])),
    ("w18_8zero", 13107220, "sparse bit-set-bvector { size 13107220, data '020180000001000200051400C80013050003C713090009'H }, default bit FALSE, sparse-other bit TRUE", Rows::Runs(&[(5, 5), (13107209, 13107209)])),
    ("w19_16zero", 19660820, "sparse bit-set-bvector { size 19660820, data '0201800000010002000514002C01130500052B0113090009'H }, default bit FALSE, sparse-other bit TRUE", Rows::Runs(&[(5, 5), (19660809, 19660809)])),
    ("w20_sblock_v3", 1000010, "sparse bit-set-bvector { size 1000010, data '020180000001000200054A420F0013030013701113F0498C13404209'H }, default bit FALSE, sparse-other bit TRUE", Rows::Runs(&[(3, 3), (70000, 70000), (150000, 150000), (1000000, 1000000)])),
    ("w21_rank", 64, "data bit '55'H, sparse bit-set-bvector { size 64, data '02018000000100020005400000003E1409FCFF0300000009'H }, default bit FALSE, sparse-other bit TRUE", Rows::Runs(&[(1, 1), (5, 5), (9, 9), (13, 13)])),
    ("w22_size_gt_bits", 48, "data bit-bvector { size 40, data '0201800000010002000528000000433B00000009'H }", Rows::Runs(&[(8, 47)])),
    ("w23_noopt_bit", 64, "data bit-bvector { size 64, data '020180000001000200054000000043F301000009'H }", Rows::Runs(&[(63, 63)])),
    ("s01_size3_hdr3", 8, "data bit-bvector { size 3, data '02018000000100020005030000004344720F0009'H }", Rows::Runs(&[(0, 0), (2, 7)])),
    ("s02_size3_hdr8", 8, "data bit-bvector { size 3, data '02018000000100020005080000004344720F0009'H }", Rows::Runs(&[(0, 0), (2, 4), (6, 7)])),
    ("s03_size8_hdr3", 8, "data bit-bvector { size 8, data '02018000000100020005030000004344720F0009'H }", Rows::Runs(&[(0, 0), (2, 4), (6, 7)])),
    ("s04_default_hdr", 8, "data bit-bvector { size 3, data '010180000001000200054344720F0009'H }", Rows::Runs(&[(0, 0), (2, 7)])),
    ("s06_size0", 8, "data bit-bvector { size 0, data '010180000001000200054344720F0009'H }", Rows::Runs(&[(0, 7)])),
    ("s07_no_bo", 8, "data bit-bvector { size 8, data '0A8000000100020005080000004344720F0009'H }", Rows::Runs(&[(0, 0), (2, 4), (6, 7)])),
    ("s08_no_gapl", 8, "data bit-bvector { size 8, data '1201080000004344720F0009'H }", Rows::Runs(&[(0, 0), (2, 4), (6, 7)])),
    ("s09_no_bo_gapl", 8, "data bit-bvector { size 8, data '1A080000004344720F0009'H }", Rows::Runs(&[(0, 0), (2, 4), (6, 7)])),
    ("s10_be", 8, "data bit-bvector { size 8, data '020000800100020005000000000813000509'H }", Rows::Runs(&[(0, 4), (6, 7)])),
    ("s11_be_garbage_bo", 8, "data bit-bvector { size 8, data '023700800100020005000000000813000509'H }", Rows::Runs(&[(0, 4), (6, 7)])),
    ("s12_be_flag_le_data", 8, "data bit-bvector { size 8, data '020080000001000200050800000013050009'H }", Rows::Runs(&[(0, 7)])),
    ("s13_idlist", 24, "data bit-bvector { size 8, data '1E08000000020000000100000014000000'H }", Rows::Runs(&[(0, 0), (2, 19), (21, 23)])),
    ("s14_trailing", 8, "data bit-bvector { size 8, data '02018000000100020005080000004344720F0009FFFFFF'H }", Rows::Runs(&[(0, 0), (2, 4), (6, 7)])),
    ("s15_8one_zero", 8, "data bit-bvector { size 8, data '0201800000010002000508000000040009'H }", Rows::Runs(&[])),
    ("s16_sblock_gaps_skip", 16777230, "sparse bit-set-bvector { size 16777230, data '020180000001000200050E0000014513050009'H }, default bit FALSE, sparse-other bit TRUE", Rows::Runs(&[(16777221, 16777221)])),
    ("s17_end_token", 8, "data bit-bvector { size 8, data '02018000000100020005080000004344720F0000'H }", Rows::Runs(&[(0, 0), (2, 4), (6, 7)])),
    ("p04_size_u32max", 8, "data bit-bvector { size 4294967295, data '010180000001000200054344720F0009'H }", Rows::Runs(&[(0, 0), (2, 4), (6, 7)])),
];

#[test]
fn bitmagic_probes_read_ncbis_rows() {
    for (case, rows, column, read) in PROBES {
        check_disabled_rows(case, column, *rows, read);
    }
}

/// The probes of the directory `BM_PROBES`, kept outside this repository (`cases/index.tsv`, the
/// PSSMs `cases/<case>.asn` and NCBI psiblast's standard error `ncbi/psiblast/<case>.err`): each
/// case w*, s* (but s05) and p04 against the rows of NCBI's posts. Without `BM_PROBES` the test
/// does nothing.
#[test]
fn bm_probes_of_a_directory() {
    let Ok(dir) = std::env::var("BM_PROBES") else {
        return;
    };
    let index = std::fs::read_to_string(format!("{dir}/cases/index.tsv")).unwrap();
    for line in index.lines().skip(1) {
        let case = line.split('\t').next().unwrap();
        if !(case.starts_with('w') || case.starts_with('s') || case.starts_with("p04"))
            || case.starts_with("s05")
        {
            continue;
        }
        let asn = std::fs::read_to_string(format!("{dir}/cases/{case}.asn")).unwrap();
        let after = |text: &'static str| &asn[asn.find(text).unwrap() + text.len()..];
        let rows: u64 = after("num-rows ")
            .split(',')
            .next()
            .unwrap()
            .parse()
            .unwrap();
        let column = after("{ header { field-name \"disabled\" }, ");
        let column = &column[..column.find(" }, { header { field-id location }").unwrap()];
        let err = std::fs::read_to_string(format!("{dir}/ncbi/psiblast/{case}.err")).unwrap();
        let mut runs: Vec<(u64, u64)> = Vec::new();
        for row in err.split("Seq-table row ").skip(1) {
            let row: u64 = row.lines().next().unwrap().parse().unwrap();
            match runs.last_mut() {
                Some(last) if last.1 + 1 == row => last.1 = row,
                _ => runs.push((row, row)),
            }
        }
        check_disabled_rows(case, column, rows, &Rows::Runs(&runs));
        eprintln!("{case}: {rows} rows");
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqtable/SeqTable_multi_data.cpp:138-169 (cited at
// `Indexer::delta_sum`)
#[test]
fn int_delta_sums_read_whole_blocks() {
    assert_eq!(
        values("data int-delta int { 1, 2, 3 }", Want::Int4, 0..4),
        ["1", "3", "6", "-"]
    );
    // Elements without integers (E13, named after the accumulator's type).
    assert_eq!(
        values("data int-delta string { \"1\", \"2\" }", Want::Int4, 0..1),
        ["E:CSeqTable_multi_data::TryGetInt4(): value cannot be converted to requested type"]
    );
    // Row 0 reads every element of its block of 128: the second element, too big for the `Int4`
    // accumulator of an `int-scaled` without bounds, throws for row 0.
    assert_eq!(
        values(
            "data int-delta int-scaled { mul 1, add 0, data int8 { 1, 4294967296 } }",
            Want::Int4,
            0..1
        ),
        ["E:CSeqTable_multi_data::TryGetInt4(): value is too big for requested type: 4294967296"]
    );
    // The `Int8` sums of `int8` elements, and an `Int4` sum that overflows (refused).
    assert_eq!(
        values("data int-delta int8 { 4294967296, 1 }", Want::Int8, 0..2),
        ["4294967296", "4294967297"]
    );
    assert_eq!(
        values("data int-delta int { 2147483647, 1 }", Want::Int8, 0..1),
        ["R:DeltaOverflow"]
    );
    // An element of a later block is not read for a row of the first.
    assert_eq!(
        values(
            &format!(
                "data int-delta int {{ {} }}",
                (0..129)
                    .map(|i| if i == 128 { "2147483647" } else { "1" })
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
            Want::Int8,
            0..1
        ),
        ["1"]
    );
}

// NCBI reference (598d8ae6): c++/src/objects/seqtable/Scaled_int_multi_data.cpp:77-84 (cited at
// `Indexer::x_try_get_int8`)
#[test]
fn int_scaled_values() {
    assert_eq!(
        values(
            "data int-scaled { mul 3, add -1, data int { 1, 2 } }",
            Want::Int4,
            0..3
        ),
        ["2", "5", "-"]
    );
    assert_eq!(
        values(
            "data int-scaled { mul 2, add 0, data int8 { 9223372036854775807 } }",
            Want::Int8,
            0..1
        ),
        ["R:ScaledOverflow"]
    );
}

/// The multi-data inside multi-data are read from a list on the heap: 50,000 nested `int-scaled`
/// (each adding 1), 50,000 nested `int-delta` values and 25,000 `int-scaled` around `int-delta`
/// give their rows on a thread of 256 KiB, the sizes found once for each value (line RG D3-1
/// review 7 S1: NCBI walks to the bottom at each level).
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn nested_multi_data_are_read_without_recursion() {
    std::thread::Builder::new()
        .stack_size(256 << 10)
        .spawn(|| {
            let read = |text: String| {
                let (value, posts) =
                    read_document(text.as_bytes(), named::SEQTABLE_COLUMN, Names::PSSM);
                assert!(posts.is_empty());
                value.expect("read")
            };
            let wrap = |n: usize, open: &str, inner: &str, close: &str| {
                format!("{}{inner}{}", open.repeat(n), close.repeat(n))
            };
            let rows = |value: &Value| {
                let column = Column(V::new(named::SEQTABLE_COLUMN, value));
                let mut ix = indexer();
                (0..3)
                    .map(|row| match ix.column_try_get(column, row, Want::Int8) {
                        Ok(v) => format!("{v:?}"),
                        Err(other) => format!("{other:?}"),
                    })
                    .collect::<Vec<_>>()
            };
            let scaled = read(format!(
                "SeqTable-column ::= {{ header {{ field-id 3 }}, data {} }}",
                wrap(
                    50_000,
                    "int-scaled { mul 1, add 1, data ",
                    "int { 4, 9 }",
                    " }"
                )
            ));
            assert_eq!(rows(&scaled), ["Some(50004)", "Some(50009)", "None"]);
            // `int-delta` of `int-delta`s: the sums of the sums.
            let delta = read(format!(
                "SeqTable-column ::= {{ header {{ field-id 3 }}, data {} }}",
                wrap(50_000, "int-delta ", "int { 4, 9 }", "")
            ));
            assert_eq!(rows(&delta), ["Some(4)", "Some(200009)", "None"]);
            let both = read(format!(
                "SeqTable-column ::= {{ header {{ field-id 3 }}, data {} }}",
                wrap(
                    25_000,
                    "int-scaled { mul 1, add 0, data int-delta ",
                    "int { 4, 9 }",
                    " }"
                )
            ));
            assert_eq!(rows(&both), ["Some(4)", "Some(100009)", "None"]);
        })
        .expect("a thread")
        .join()
        .expect("no overflow");
}
