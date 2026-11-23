//! Tests of mask files whose Seq-locs nest to any depth (PD-LOSAT-NCBI-DEFECTS exception 28).
//!
//! NCBI's makeblastdb reads a `mix` or `equiv` by recursion and, with its default stack of 8 MiB,
//! crashes on nests of some thousands of levels; with a stack large enough it reads every depth.
//! The files here are those of `docs/evidence/mask_nest_depth/gen_mask_nests.py`: the builders
//! below write the same bytes for each of its kinds (`the_builders_write_the_fixtures_bytes`), and
//! the fixtures `n_*` of `tests/seqdb/mask_parse` are NCBI's output for them with no stack limit
//! (`expected.tsv`). The elements stated as `NULL` and the `:Seq-interval` element are the tests'
//! own kinds (no fixture: LOSAT refuses them, rule 3). The deep cases run on threads of 8 MiB and
//! 2 MiB, below what a recursion through them would take.

use super::*;
use crate::seqid::model::ObjectId;
use std::path::{Path, PathBuf};

/// The ways a Seq-loc nests in a mask file (`gen_mask_nests.py`'s nests).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Shape {
    /// `mix { mix { ... } }`
    Mix,
    /// `equiv { equiv { ... } }`
    Equiv,
    /// mix and equiv in turn, the outermost a mix
    Alt,
    /// mixes whose element is stated as `:Seq-loc`
    Ptr,
    /// mixes that each hold an interval of `local id 3` before the mix inside
    Wide,
    /// `Wide` whose innermost interval is of `local id 4`
    WideId,
}

impl Shape {
    fn is_equiv(self, level: usize) -> bool {
        self == Shape::Equiv || (self == Shape::Alt && level % 2 == 1)
    }

    /// The frame names that one level adds to a `CSerialException` (two levels for `Alt`).
    fn unit(self) -> &'static str {
        match self {
            Shape::Equiv => ".equiv.[].[]",
            Shape::Alt => ".mix.[].[].equiv.[].[]",
            Shape::Ptr => ".mix.[].[].Seq-loc",
            _ => ".mix.[].[]",
        }
    }
}

/// The innermost Seq-loc of a nest.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Inner {
    /// `int { from 5, to 14, id local id 3 }` (`local id 4` in a `WideId` nest)
    Interval,
    /// a variant that Seq-loc does not have
    BadVar,
    /// an interval whose strand is 7
    BadStrand,
    /// an interval whose `from` is too large
    BigNum,
    /// the interval, the file cut inside it
    Cut,
    /// a mix whose only element is stated as `NULL` (the tests' own)
    NullMix,
    /// an equiv whose only element is stated as `NULL` (the tests' own)
    NullEquiv,
}

/// The damage at the middle level `depth / 2`, after the levels inside it were read whole
/// (`gen_mask_nests.py`'s `mid_*`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Middle {
    /// text: an interval after the element, without a comma
    Comma,
    /// a second element of a variant that Seq-loc does not have
    BadVar,
    /// the file cut after the element
    Cut,
    /// the list opened wrongly (text `[`, binary a SET OF tag)
    Open,
    /// binary: the end of the list damaged
    ContainerEnd,
    /// binary: the end of the variant damaged
    VariantEnd,
    /// binary, a `Ptr` nest: the end of the `:Seq-loc` element damaged
    PointerEnd,
    /// text: the element stated as `:Seq-interval` (the tests' own: LOSAT refuses it, as at any
    /// depth)
    OtherType,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Nest {
    shape: Shape,
    inner: Inner,
    middle: Option<Middle>,
}

const fn nest(shape: Shape, inner: Inner) -> Nest {
    Nest {
        shape,
        inner,
        middle: None,
    }
}

const fn mid(shape: Shape, middle: Middle) -> Nest {
    Nest {
        shape,
        inner: Inner::Interval,
        middle: Some(middle),
    }
}

impl Nest {
    /// The kind's name in `gen_mask_nests.py`.
    fn name(self) -> String {
        let shape = match self.shape {
            Shape::Mix => "mix",
            Shape::Equiv => "equiv",
            Shape::Alt => "alt",
            Shape::Ptr => "ptr",
            Shape::Wide => "wide",
            Shape::WideId => "wideid",
        };
        if let Some(middle) = self.middle {
            let damage = match middle {
                Middle::Comma => "comma",
                Middle::BadVar => "badvar",
                Middle::Cut => "cut",
                Middle::Open => "open",
                Middle::ContainerEnd => "cend",
                Middle::VariantEnd => "vend",
                Middle::PointerEnd => return "mid_pend".to_string(),
                Middle::OtherType => panic!("the tests' own kind"),
            };
            return match self.shape {
                Shape::Mix => format!("mid_{damage}"),
                _ => format!("mid_{damage}_{shape}"),
            };
        }
        let error = match self.inner {
            Inner::Interval => return shape.to_string(),
            Inner::BadVar => "badvar",
            Inner::BadStrand => "badstrand",
            Inner::BigNum => "bignum",
            Inner::Cut => "cut",
            Inner::NullMix | Inner::NullEquiv => panic!("the tests' own kind"),
        };
        match self.shape {
            Shape::Mix => error.to_string(),
            _ => format!("{error}_{shape}"),
        }
    }

    /// Whether `gen_mask_nests.py` writes the kind in text (`TEXT_ONLY`, `BINARY_ONLY`).
    fn has_text(self) -> bool {
        !matches!(
            self.middle,
            Some(Middle::ContainerEnd | Middle::VariantEnd | Middle::PointerEnd)
        )
    }

    fn has_binary(self) -> bool {
        !matches!(self.middle, Some(Middle::Comma | Middle::OtherType))
    }
}

const TEXT_HEAD: &str = "Blast-db-mask-info ::= { algo-id 0, algo-program windowmasker, algo-options \"x\", masks { masks {\n";

fn interval_text(from: &str, strand: &str, id: i32) -> String {
    format!("int {{ from {from}, to 14, {strand}id local id {id} }}")
}

/// The ASN.1 text of a mask nested `depth` levels as `nest`, one level a line.
fn nest_text(nest: Nest, depth: usize) -> Vec<u8> {
    let id = if nest.shape == Shape::WideId { 4 } else { 3 };
    let inner = match nest.inner {
        Inner::BadStrand => interval_text("5", "strand 7, ", id),
        Inner::BigNum => interval_text("99999999999", "", id),
        Inner::BadVar => "mox { }".to_string(),
        Inner::NullMix => "mix { NULL }".to_string(),
        Inner::NullEquiv => "equiv { NULL }".to_string(),
        Inner::Interval | Inner::Cut => interval_text("5", "", id),
    };
    let m = depth / 2;
    let mut out = String::with_capacity(depth * 48 + 256);
    out.push_str(TEXT_HEAD);
    for level in 0..depth {
        let brace = if nest.middle == Some(Middle::Open) && level == m {
            "["
        } else {
            "{"
        };
        if nest.shape.is_equiv(level) {
            out.push_str(&format!("equiv {brace}\n"));
        } else if nest.middle == Some(Middle::OtherType) && level == m {
            out.push_str("mix { :Seq-interval\n");
        } else if nest.shape == Shape::Ptr {
            out.push_str(&format!("mix {brace} :Seq-loc\n"));
        } else if matches!(nest.shape, Shape::Wide | Shape::WideId) {
            out.push_str(&format!("mix {brace} {},\n", interval_text("5", "", 3)));
        } else {
            out.push_str(&format!("mix {brace}\n"));
        }
    }
    out.push_str(&inner);
    // the j-th `}` ends the list of level depth - 1 - j: that of level m + 1 is the last before
    // level m's next element
    let after = depth.wrapping_sub(2 + m);
    for j in 0..depth {
        out.push_str("\n}");
        if j == after {
            match nest.middle {
                Some(Middle::Comma) => out.push_str(&format!(" {}", interval_text("5", "", 3))),
                Some(Middle::BadVar) => out.push_str(", mox { }"),
                Some(Middle::Cut) => {
                    out.push('\n');
                    return out.into_bytes();
                }
                _ => {}
            }
        }
    }
    out.push_str("\n}, more FALSE } }\n");
    if nest.inner == Inner::Cut {
        out.truncate(out.rfind("to 14").unwrap());
    }
    out.into_bytes()
}

fn ber_mem(n: u8, content: &[u8]) -> Vec<u8> {
    [&[0xa0 | n, 0x80], content, &[0, 0]].concat()
}

fn ber_seq(tag: u8, members: &[Vec<u8>]) -> Vec<u8> {
    [&[tag, 0x80], members.concat().as_slice(), &[0, 0]].concat()
}

fn ber_int(n: u8) -> Vec<u8> {
    vec![0x02, 0x01, n]
}

fn interval_binary(from: Vec<u8>, id: u8, strand: Option<u8>) -> Vec<u8> {
    let mut members = vec![ber_mem(0, &from), ber_mem(1, &ber_int(14))];
    if let Some(strand) = strand {
        members.push(ber_mem(2, &[0x0a, 0x01, strand]));
    }
    members.push(ber_mem(3, &ber_mem(0, &ber_mem(0, &ber_int(id)))));
    ber_seq(0x30, &members)
}

/// `:Seq-loc`: the application long tag that spells the type name, with an indefinite length.
const OTHER_SEQ_LOC: &[u8] = b"\x7f\xd3\xe5\xf1\xad\xec\xefc\x80";

/// A Seq-loc of the variant [15], which Seq-loc does not have.
const BAD_VARIANT: &[u8] = &[0xaf, 0x80, 0x05, 0x00, 0x00, 0x00];

/// Where the Seq-loc of the one mask starts in a file of `nest_binary`: after the members
/// `algo-id`, `algo-program`, `algo-options` and the openings of `masks`, `Blast-mask-list.masks`
/// and its list.
const LOC_AT: usize = 31;

/// The ASN.1 binary (BER, indefinite lengths) of a mask nested `depth` levels as `nest`.
fn nest_binary(nest: Nest, depth: usize) -> Vec<u8> {
    let id = if nest.shape == Shape::WideId { 4 } else { 3 };
    let inner = match nest.inner {
        Inner::BadStrand => ber_mem(3, &interval_binary(ber_int(5), id, Some(7))),
        Inner::BigNum => ber_mem(
            3,
            &interval_binary([&[0x02, 0x09][..], &[1; 9]].concat(), id, None),
        ),
        Inner::BadVar => BAD_VARIANT.to_vec(),
        Inner::NullMix => ber_mem(7, &ber_seq(0x30, &[vec![0x05, 0x00]])),
        Inner::NullEquiv => ber_mem(8, &ber_seq(0x31, &[vec![0x05, 0x00]])),
        Inner::Interval | Inner::Cut => ber_mem(3, &interval_binary(ber_int(5), id, None)),
    };
    let wide = ber_mem(3, &interval_binary(ber_int(5), 3, None));
    let m = depth / 2;
    // the ends of a level: those of the `:Seq-loc` element (a `Ptr` nest), the list and the
    // variant
    let end = if nest.shape == Shape::Ptr { 6 } else { 4 };
    let mut loc = Vec::with_capacity(depth * (wide.len() + 4 + end) + inner.len());
    for level in 0..depth {
        let list = if nest.middle == Some(Middle::Open) && level == m {
            0x31
        } else {
            0x30
        };
        if nest.shape.is_equiv(level) {
            loc.extend_from_slice(&[0xa8, 0x80, 0x31, 0x80]);
            continue;
        }
        loc.extend_from_slice(&[0xa7, 0x80, list, 0x80]);
        match nest.shape {
            Shape::Ptr => loc.extend_from_slice(OTHER_SEQ_LOC),
            Shape::Wide | Shape::WideId => loc.extend_from_slice(&wide),
            _ => {}
        }
    }
    let head = loc.len();
    loc.extend_from_slice(&inner);
    loc.resize(loc.len() + depth * end, 0);
    // the ends of level m follow those of the levels inside it
    let level_m = head + inner.len() + (depth - 1 - m) * end;
    match (nest.middle, nest.shape) {
        (Some(Middle::ContainerEnd), Shape::Mix) => loc[level_m + 1] = 5,
        (Some(Middle::VariantEnd), Shape::Mix) => loc[level_m + 2] = 5,
        (Some(Middle::PointerEnd), Shape::Ptr) => loc[level_m] = 5,
        (Some(Middle::BadVar), _) => {
            let ends = loc.split_off(level_m);
            loc.extend_from_slice(BAD_VARIANT);
            loc.extend_from_slice(&ends);
        }
        (None | Some(Middle::Cut | Middle::Open), _) => {}
        other => panic!("no such nest: {other:?}"),
    }
    let masks = ber_seq(
        0x30,
        &[
            ber_mem(0, &ber_seq(0x30, &[loc])),
            ber_mem(1, &[0x01, 0x01, 0x00]),
        ],
    );
    let mut out = ber_seq(
        0x30,
        &[
            ber_mem(0, &ber_int(0)),
            ber_mem(1, &ber_int(30)),
            ber_mem(2, &[0x1a, 0x01, b'x']),
            ber_mem(3, &masks),
        ],
    );
    if nest.inner == Inner::Cut {
        out.truncate(head + 40);
    }
    if nest.middle == Some(Middle::Cut) {
        out.truncate(LOC_AT + level_m);
    }
    out
}

/// The file's first object, as `MaskFileReader` reads it, on a thread with `stack` bytes of stack.
fn read_on(stack: usize, data: Vec<u8>) -> Result<MaskInfo, MaskFileError> {
    std::thread::Builder::new()
        .stack_size(stack)
        .spawn(move || {
            let mut reader = MaskFileReader::from_bytes(Path::new("nest"), data)?;
            Ok(reader.next_info()?.expect("one object"))
        })
        .unwrap()
        .join()
        .expect("the reader does not overflow its stack")
}

const MIB: usize = 1 << 20;

/// The depth of the deep cases: past a million levels, and odd as 401 is, so that an `Alt` nest
/// ends in the same variant.
const DEEP: usize = 1_000_001;

fn local(n: i32) -> Option<crate::seqid::ParsedSeqId> {
    Some(crate::seqid::ParsedSeqId::Local(ObjectId::Id(n)))
}

const FIVE_TO_FOURTEEN: MaskPart = MaskPart::Interval(MaskInterval {
    from: 5,
    to: 14,
    strand: None,
});

/// The nests of the fixtures (`make_mask_parse_fixtures.py`'s `NESTS`), with their depths.
const FIXTURES: &[(Nest, usize)] = &[
    (nest(Shape::Mix, Inner::Interval), 401),
    (nest(Shape::Equiv, Inner::Interval), 401),
    (nest(Shape::Alt, Inner::Interval), 401),
    (nest(Shape::Ptr, Inner::Interval), 401),
    (nest(Shape::Wide, Inner::Interval), 401),
    (nest(Shape::WideId, Inner::Interval), 401),
    (nest(Shape::Mix, Inner::BadVar), 401),
    (nest(Shape::Mix, Inner::BadStrand), 401),
    (nest(Shape::Mix, Inner::Cut), 401),
    (nest(Shape::Mix, Inner::BigNum), 101),
    (nest(Shape::Mix, Inner::Interval), 20_000),
    (nest(Shape::Equiv, Inner::BadVar), 401),
    (nest(Shape::Alt, Inner::BadVar), 401),
    (nest(Shape::Ptr, Inner::BadVar), 401),
    (nest(Shape::Equiv, Inner::Cut), 401),
    (nest(Shape::Alt, Inner::Cut), 401),
    (nest(Shape::Ptr, Inner::Cut), 401),
    (nest(Shape::Alt, Inner::BigNum), 101),
    (nest(Shape::Ptr, Inner::BigNum), 101),
    (mid(Shape::Mix, Middle::Comma), 401),
    (mid(Shape::Mix, Middle::BadVar), 401),
    (mid(Shape::Mix, Middle::Cut), 401),
    (mid(Shape::Mix, Middle::Open), 401),
    (mid(Shape::Mix, Middle::ContainerEnd), 401),
    (mid(Shape::Mix, Middle::VariantEnd), 401),
    (mid(Shape::Ptr, Middle::PointerEnd), 401),
    (mid(Shape::Ptr, Middle::Comma), 401),
];

fn fixture_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/mask_parse")
}

/// The builders write the bytes of the fixtures (`gen_mask_nests.py`), so that the deep cases
/// below are the fixtures' nests made deeper; every fixture `n_*` is one of `FIXTURES`.
#[test]
fn the_builders_write_the_fixtures_bytes() {
    let mut names = Vec::new();
    for &(nest, depth) in FIXTURES {
        let name = format!("n_{}_{depth}", nest.name());
        if nest.has_text() {
            let file = fixture_dir().join(format!("{name}.asnt"));
            let committed = std::fs::read(file).unwrap_or_else(|e| panic!("{name}: {e}"));
            assert!(nest_text(nest, depth) == committed, "{name}.asnt");
            names.push(format!("{name}.asnt"));
        }
        if nest.has_binary() {
            let file = fixture_dir().join(format!("{name}.asnb"));
            let committed = std::fs::read(file).unwrap_or_else(|e| panic!("{name}: {e}"));
            assert!(nest_binary(nest, depth) == committed, "{name}.asnb");
            names.push(format!("{name}.asnb"));
        }
    }
    let mut files: Vec<String> = std::fs::read_dir(fixture_dir())
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|name| name.starts_with("n_"))
        .collect();
    files.sort();
    names.sort();
    assert_eq!(files, names);
}

#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn a_nest_of_a_million_levels_is_read_without_recursion() {
    let mix = nest(Shape::Mix, Inner::Interval);
    for stack in [8 * MIB, 2 * MIB] {
        for (format, data) in [
            ("text", nest_text(mix, 1_000_000)),
            ("binary", nest_binary(mix, 1_000_000)),
        ] {
            let info = read_on(stack, data).unwrap_or_else(|e| panic!("{format}: {e}"));
            assert_eq!(info.masks.masks.len(), 1, "{format}");
            let mask = &info.masks.masks[0];
            assert_eq!(mask.id, local(3), "{format}");
            assert_eq!(mask.parts, [FIVE_TO_FOURTEEN], "{format}");
            assert!(!info.masks.more && info.masks.null_element.is_none());
        }
    }
}

#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn every_nesting_path_is_read_at_depth() {
    const DEPTH: usize = 100_000;
    for shape in [
        Shape::Equiv,
        Shape::Alt,
        Shape::Ptr,
        Shape::Wide,
        Shape::WideId,
    ] {
        let kind = nest(shape, Inner::Interval);
        for (format, data) in [
            ("text", nest_text(kind, DEPTH)),
            ("binary", nest_binary(kind, DEPTH)),
        ] {
            let info = read_on(2 * MIB, data).unwrap_or_else(|e| panic!("{shape:?} {format}: {e}"));
            let mask = &info.masks.masks[0];
            match shape {
                Shape::Wide => {
                    assert_eq!(mask.id, local(3));
                    assert_eq!(mask.parts, vec![FIVE_TO_FOURTEEN; DEPTH + 1]);
                }
                // two Seq-ids: `GetId` has none ("Cannot get masked range Seq-id")
                Shape::WideId => {
                    assert_eq!(mask.id, None);
                    assert_eq!(mask.parts.len(), DEPTH + 1);
                }
                _ => {
                    assert_eq!(mask.id, local(3), "{shape:?} {format}");
                    assert_eq!(mask.parts, [FIVE_TO_FOURTEEN], "{shape:?} {format}");
                }
            }
        }
    }
}

/// The element stated as `NULL` deep inside a mix or an equiv is where `CheckId` dereferences it
/// (refused, rule 3, as at any depth).
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn a_null_element_deep_inside_is_found_where_ncbi_dereferences_it() {
    for (kind, list) in [
        (nest(Shape::Mix, Inner::NullMix), "a Seq-loc-mix (mix)"),
        (
            nest(Shape::Equiv, Inner::NullEquiv),
            "a Seq-loc-equiv (equiv)",
        ),
    ] {
        for (format, data) in [
            ("text", nest_text(kind, 1_000_000)),
            ("binary", nest_binary(kind, 1_000_000)),
        ] {
            let info = read_on(2 * MIB, data).unwrap();
            assert!(info.masks.masks.is_empty(), "{kind:?} {format}");
            let null = info.masks.null_element.expect("the NULL element");
            assert_eq!(null.list, list, "{kind:?} {format}");
        }
    }
}

/// What makeblastdb reports for the error of `nest`, `depth` levels deep.
fn report(nest: Nest, binary: bool, depth: usize) -> String {
    let data = if binary {
        nest_binary(nest, depth)
    } else {
        nest_text(nest, depth)
    };
    match read_on(2 * MIB, data) {
        Err(MaskFileError::Ncbi(e)) => e.to_string(),
        other => panic!("{nest:?} {depth}: {other:?}"),
    }
}

/// How many times `unit` repeats where it first occurs in `report`.
fn run_of(report: &str, unit: &str) -> usize {
    let Some(start) = report.find(unit) else {
        return 0;
    };
    report[start..]
        .as_bytes()
        .chunks(unit.len())
        .take_while(|chunk| *chunk == unit.as_bytes())
        .count()
}

/// `report` with the frame names of the levels of its nest (`count` times `unit`) as `<NEST>` and
/// the numbers of the stream's positions as `N` ([`positions`]).
fn normalized(report: &str, unit: &str, count: usize) -> String {
    positions(&report.replace(&unit.repeat(count), "<NEST>"))
}

/// `report` with the numbers of the stream's positions (`line N`, `byte N`) as `N`.
fn positions(report: &str) -> String {
    let mut out = String::with_capacity(report.len());
    let mut rest = report;
    while let Some(at) = ["line ", "byte "]
        .iter()
        .filter_map(|word| rest.find(word).map(|i| i + word.len()))
        .min()
    {
        out.push_str(&rest[..at]);
        rest = &rest[at..];
        let digits = rest.bytes().take_while(u8::is_ascii_digit).count();
        if digits > 0 {
            out.push('N');
        }
        rest = &rest[digits..];
    }
    out.push_str(rest);
    out
}

/// The levels that `nest` adds to the frame stack of its error between `small` and `deep` levels
/// (counted in `Shape::unit`s): every level for a damage innermost, those outside the middle
/// level for a damage there.
fn levels_added(nest: Nest, small: usize, deep: usize) -> usize {
    let (small, deep) = match nest.middle {
        Some(_) => (small / 2, deep / 2),
        None => (small, deep),
    };
    match nest.shape {
        Shape::Alt => deep / 2 - small / 2,
        _ => deep - small,
    }
}

/// The error of `nest` at `DEEP` levels is that of its fixture's depth (NCBI's report with no
/// stack limit, `expected.tsv`) with every level added: the frame stack of a `CSerialException`
/// (with the copies that `ReadWithSetFlagMember` freezes) holds each level's frames, built in time
/// linear in its length.
fn assert_deep_report(nest: Nest, small: usize) {
    for binary in [false, true] {
        if (binary && !nest.has_binary()) || (!binary && !nest.has_text()) {
            continue;
        }
        let unit = nest.shape.unit();
        let small_report = report(nest, binary, small);
        let deep_report = report(nest, binary, DEEP);
        let (a, b) = (run_of(&small_report, unit), run_of(&deep_report, unit));
        assert!(a > 0, "{nest:?} {binary}: {small_report}");
        assert_eq!(b - a, levels_added(nest, small, DEEP), "{nest:?} {binary}");
        assert_eq!(
            normalized(&deep_report, unit, b),
            normalized(&small_report, unit, a),
            "{nest:?} {binary}"
        );
    }
}

/// A damage innermost in a `mix` nest leaves every frame of it, innermost first.
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn an_error_deep_inside_leaves_every_frame_of_the_nest() {
    for inner in [Inner::BadVar, Inner::BadStrand, Inner::Cut] {
        assert_deep_report(nest(Shape::Mix, inner), 401);
    }
    // binary only: the text report takes a line per frame (below)
    let bignum = nest(Shape::Mix, Inner::BigNum);
    let (small, deep) = (report(bignum, true, 101), report(bignum, true, DEEP));
    let unit = Shape::Mix.unit();
    let (a, b) = (run_of(&small, unit), run_of(&deep, unit));
    assert_eq!(b - a, DEEP - 101);
    assert_eq!(normalized(&deep, unit, b), normalized(&small, unit, a));
}

/// A damage innermost in an `equiv`, `alt` or `ptr` nest: the frames of each variant and of each
/// `:Seq-loc` element.
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn an_error_deep_inside_other_nests_leaves_every_frame() {
    for shape in [Shape::Equiv, Shape::Alt, Shape::Ptr] {
        for inner in [Inner::BadVar, Inner::Cut] {
            assert_deep_report(nest(shape, inner), 401);
        }
    }
    for shape in [Shape::Alt, Shape::Ptr] {
        let bignum = nest(shape, Inner::BigNum);
        let (small, deep) = (report(bignum, true, 101), report(bignum, true, DEEP));
        let unit = shape.unit();
        let (a, b) = (run_of(&small, unit), run_of(&deep, unit));
        assert_eq!(b - a, levels_added(bignum, 101, DEEP), "{shape:?}");
        assert_eq!(
            normalized(&deep, unit, b),
            normalized(&small, unit, a),
            "{shape:?}"
        );
    }
}

/// A damage at the middle level, at each step that ends a level (the end of an element, of the
/// list, of the variant, of a `:Seq-loc` element) and at its start: the frames of the levels
/// outside it, those of the level itself as far as it got.
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn an_error_at_a_middle_level_leaves_the_frames_outside_it() {
    for (shape, middle) in [
        (Shape::Mix, Middle::Comma),
        (Shape::Mix, Middle::BadVar),
        (Shape::Mix, Middle::Cut),
        (Shape::Mix, Middle::Open),
        (Shape::Mix, Middle::ContainerEnd),
        (Shape::Mix, Middle::VariantEnd),
        (Shape::Ptr, Middle::PointerEnd),
        (Shape::Ptr, Middle::Comma),
    ] {
        assert_deep_report(mid(shape, middle), 401);
    }
}

/// A `CUtilException` (a number too large in text) takes one line for each frame it leaves: the
/// report of 20,101 levels is that of 101 (the fixtures) with each level's lines repeated.
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn a_number_too_large_takes_a_line_for_each_frame() {
    const MORE: usize = 20_000;
    for (shape, lines_per_level) in [(Shape::Mix, 4), (Shape::Alt, 4), (Shape::Ptr, 5)] {
        let bignum = nest(shape, Inner::BigNum);
        let small = positions(&report(bignum, false, 101));
        let deep = positions(&report(bignum, false, 101 + MORE));
        let small: Vec<&str> = small.lines().collect();
        let deep: Vec<&str> = deep.lines().collect();
        let extra = lines_per_level * MORE;
        assert_eq!(deep.len(), small.len() + extra, "{shape:?}");
        // the lines of the levels repeat (with the period of a level, two for `Alt`): taking
        // `extra` of them out of the middle gives the report of 101 levels
        let cut = small.len() / 2;
        let joined: Vec<&str> = deep[..cut]
            .iter()
            .chain(&deep[cut + extra..])
            .copied()
            .collect();
        assert_eq!(joined, small, "{shape:?}");
    }
}

/// A value of another type than Seq-loc stated as an element at the middle level is refused where
/// it stands, as at any depth (the frames around it change nothing).
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn a_value_of_another_type_at_a_middle_level_is_refused_where_it_stands() {
    for depth in [401, DEEP] {
        let data = nest_text(mid(Shape::Mix, Middle::OtherType), depth);
        match read_on(2 * MIB, data) {
            Err(MaskFileError::Unsupported(what)) => {
                let line = format!("(line {})", depth / 2 + 2);
                assert!(
                    what.starts_with("':Seq-interval'") && what.contains(&line),
                    "{depth}: {what}"
                );
            }
            other => panic!("{depth}: {other:?}"),
        }
    }
}

/// Strings that post in many elements of one object: each post gets its stack trace, and the
/// reader takes time linear in their number (the traces that wait for a frame are the last ones
/// made, `Rd::trace_frame`).
#[test]
fn many_strings_that_post_get_their_traces() {
    const ELEMENTS: usize = 100_000;
    let mut text = String::from(TEXT_HEAD);
    text.push_str("mix {\n");
    for i in 0..ELEMENTS {
        if i > 0 {
            text.push_str(",\n");
        }
        text.push_str("int { from 5, to 14, id local str \"a\tb\" }");
    }
    text.push_str("\n}\n}, more FALSE } }\n");
    let mut reader = MaskFileReader::from_bytes(Path::new("posts"), text.into_bytes()).unwrap();
    let info = reader.next_info().unwrap().expect("one object");
    assert_eq!(info.masks.masks.len(), 1);
    let posts = reader.take_posts();
    assert_eq!(posts.len(), ELEMENTS);
    for post in &posts {
        let MaskPost::BadChar(replaced) = post else {
            panic!("{post:?}")
        };
        assert_eq!(
            replaced.trace,
            "Blast-db-mask-info.masks.masks.E.mix.E.int.id.local.str"
        );
    }
}
