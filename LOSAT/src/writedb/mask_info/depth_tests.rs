//! Tests of mask files whose Seq-locs nest to any depth (PD-LOSAT-NCBI-DEFECTS exception 24).
//!
//! NCBI's makeblastdb reads a `mix` or `equiv` by recursion and, with its default stack of 8 MiB,
//! crashes on nests of some thousands of levels; with a stack large enough it reads every depth.
//! The files here are those of `docs/evidence/mask_nest_depth/gen_mask_nests.py` (the builders
//! below write the same bytes); the fixtures `n_*` of `tests/seqdb/mask_parse` are NCBI's output
//! for some of them with no stack limit (`expected.tsv`). The deep cases run on threads of 8 MiB
//! and 2 MiB, below what a recursion through them would take.

use super::*;
use crate::seqid::model::ObjectId;
use std::path::Path;

/// The ways a Seq-loc nests in a mask file (`gen_mask_nests.py`'s kinds).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    /// `mix { mix { ... int } }`
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
    /// around a variant that Seq-loc does not have
    BadVar,
    /// around an interval whose strand is 7
    BadStrand,
    /// around an interval whose `from` is too large
    BigNum,
    /// cut inside the innermost interval
    Cut,
    /// around a mix whose only element is stated as `NULL`
    Null,
}

const TEXT_HEAD: &str = "Blast-db-mask-info ::= { algo-id 0, algo-program windowmasker, algo-options \"x\", masks { masks {\n";

fn interval_text(from: &str, strand: &str, id: i32) -> String {
    format!("int {{ from {from}, to 14, {strand}id local id {id} }}")
}

/// The ASN.1 text of a mask nested `depth` levels as `kind`, one level a line.
fn nest_text(kind: Kind, depth: usize) -> Vec<u8> {
    let inner = match kind {
        Kind::BadStrand => interval_text("5", "strand 7, ", 3),
        Kind::BigNum => interval_text("99999999999", "", 3),
        Kind::WideId => interval_text("5", "", 4),
        Kind::BadVar => "mox { }".to_string(),
        Kind::Null => "mix { NULL }".to_string(),
        _ => interval_text("5", "", 3),
    };
    let mut out = String::with_capacity(depth * 48 + 256);
    out.push_str(TEXT_HEAD);
    for level in 0..depth {
        out.push_str(match kind {
            Kind::Equiv => "equiv {\n",
            Kind::Alt if level % 2 == 1 => "equiv {\n",
            Kind::Ptr => "mix { :Seq-loc\n",
            Kind::Wide | Kind::WideId => "mix { int { from 5, to 14, id local id 3 },\n",
            _ => "mix {\n",
        });
    }
    out.push_str(&inner);
    for _ in 0..depth {
        out.push_str("\n}");
    }
    out.push_str("\n}, more FALSE } }\n");
    if kind == Kind::Cut {
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

/// The ASN.1 binary (BER, indefinite lengths) of a mask nested `depth` levels as `kind`.
fn nest_binary(kind: Kind, depth: usize) -> Vec<u8> {
    let inner = match kind {
        Kind::BadStrand => ber_mem(3, &interval_binary(ber_int(5), 3, Some(7))),
        Kind::BigNum => ber_mem(
            3,
            &interval_binary([&[0x02, 0x09][..], &[1; 9]].concat(), 3, None),
        ),
        Kind::WideId => ber_mem(3, &interval_binary(ber_int(5), 4, None)),
        Kind::BadVar => ber_mem(15, &[0x05, 0x00]),
        Kind::Null => ber_mem(7, &ber_seq(0x30, &[vec![0x05, 0x00]])),
        _ => ber_mem(3, &interval_binary(ber_int(5), 3, None)),
    };
    let wide = [
        &[0xa7, 0x80, 0x30, 0x80][..],
        &ber_mem(3, &interval_binary(ber_int(5), 3, None)),
    ]
    .concat();
    let mut loc = Vec::with_capacity(depth * (wide.len() + 4) + inner.len());
    for level in 0..depth {
        match kind {
            Kind::Equiv => loc.extend_from_slice(&[0xa8, 0x80, 0x31, 0x80]),
            Kind::Alt if level % 2 == 1 => loc.extend_from_slice(&[0xa8, 0x80, 0x31, 0x80]),
            Kind::Ptr => {
                loc.extend_from_slice(&[0xa7, 0x80, 0x30, 0x80]);
                loc.extend_from_slice(OTHER_SEQ_LOC);
            }
            Kind::Wide | Kind::WideId => loc.extend_from_slice(&wide),
            _ => loc.extend_from_slice(&[0xa7, 0x80, 0x30, 0x80]),
        }
    }
    let head = loc.len();
    loc.extend_from_slice(&inner);
    let end = if kind == Kind::Ptr { 6 } else { 4 };
    loc.resize(loc.len() + depth * end, 0);
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
    if kind == Kind::Cut {
        out.truncate(head + 40);
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

fn local(n: i32) -> Option<crate::seqid::ParsedSeqId> {
    Some(crate::seqid::ParsedSeqId::Local(ObjectId::Id(n)))
}

const FIVE_TO_FOURTEEN: MaskPart = MaskPart::Interval(MaskInterval {
    from: 5,
    to: 14,
    strand: None,
});

#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn a_nest_of_a_million_levels_is_read_without_recursion() {
    for stack in [8 * MIB, 2 * MIB] {
        for (format, data) in [
            ("text", nest_text(Kind::Mix, 1_000_000)),
            ("binary", nest_binary(Kind::Mix, 1_000_000)),
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
    for kind in [Kind::Equiv, Kind::Alt, Kind::Ptr, Kind::Wide, Kind::WideId] {
        for (format, data) in [
            ("text", nest_text(kind, DEPTH)),
            ("binary", nest_binary(kind, DEPTH)),
        ] {
            let info = read_on(2 * MIB, data).unwrap_or_else(|e| panic!("{kind:?} {format}: {e}"));
            let mask = &info.masks.masks[0];
            match kind {
                Kind::Wide => {
                    assert_eq!(mask.id, local(3));
                    assert_eq!(mask.parts, vec![FIVE_TO_FOURTEEN; DEPTH + 1]);
                }
                // two Seq-ids: `GetId` has none ("Cannot get masked range Seq-id")
                Kind::WideId => {
                    assert_eq!(mask.id, None);
                    assert_eq!(mask.parts.len(), DEPTH + 1);
                }
                _ => {
                    assert_eq!(mask.id, local(3), "{kind:?} {format}");
                    assert_eq!(mask.parts, [FIVE_TO_FOURTEEN], "{kind:?} {format}");
                }
            }
        }
    }
}

/// The element stated as `NULL` deep inside is where `CheckId` dereferences it (refused, rule 3,
/// as at any depth).
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn a_null_element_deep_inside_is_found_where_ncbi_dereferences_it() {
    for (format, data) in [
        ("text", nest_text(Kind::Null, 1_000_000)),
        ("binary", nest_binary(Kind::Null, 1_000_000)),
    ] {
        let info = read_on(2 * MIB, data).unwrap();
        assert!(info.masks.masks.is_empty(), "{format}");
        let null = info.masks.null_element.expect("the NULL element");
        assert_eq!(null.list, "a Seq-loc-mix (mix)", "{format}");
    }
}

/// What makeblastdb reports for the error of a nest of `depth` levels of `kind`.
fn report(kind: Kind, binary: bool, depth: usize) -> String {
    let data = if binary {
        nest_binary(kind, depth)
    } else {
        nest_text(kind, depth)
    };
    match read_on(2 * MIB, data) {
        Err(MaskFileError::Ncbi(e)) => e.to_string(),
        other => panic!("{kind:?} {depth}: {other:?}"),
    }
}

/// `report` with the frame names of the `depth` levels of a mix nest (each level adds the names of
/// its variant, container and element frames to a `CSerialException`'s frame stack) as `<NEST>`
/// and the numbers of the stream's positions (`line N`, `byte N`) as `N`.
fn normalized(report: &str, depth: usize) -> String {
    let report = report.replace(&".mix.[].[]".repeat(depth), ".<NEST>");
    let mut out = String::with_capacity(report.len());
    let mut rest = report.as_str();
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

/// An error deep inside a nest leaves every frame of it, innermost first: the frame stack of a
/// `CSerialException` (with the copies that `ReadWithSetFlagMember` freezes) holds every level,
/// built in time linear in its length; another `CException` takes a line at each frame.
#[test]
#[cfg_attr(
    all(target_os = "wasi", not(target_feature = "atomics")),
    ignore = "needs threads: wasm32-wasip1 has none"
)]
fn an_error_deep_inside_leaves_every_frame_of_the_nest() {
    const DEPTH: usize = 1_000_000;
    for kind in [Kind::BadVar, Kind::BadStrand, Kind::Cut] {
        for binary in [false, true] {
            let small = report(kind, binary, 401);
            assert!(
                small.contains(&".mix.[].[]".repeat(401)),
                "{kind:?} {binary}: {small}"
            );
            let deep = report(kind, binary, DEPTH);
            assert!(
                deep.contains(&".mix.[].[]".repeat(DEPTH)),
                "{kind:?} {binary}"
            );
            assert_eq!(
                normalized(&deep, DEPTH),
                normalized(&small, 401),
                "{kind:?} {binary}"
            );
        }
    }
    // a CUtilException (text) takes one line for each frame it leaves
    let lines = |depth| report(Kind::BigNum, false, depth).lines().count();
    let (a, b) = (lines(401), lines(20_401));
    assert_eq!((b - a) % 20_000, 0, "{a} {b}");
    assert!(b - a >= 4 * 20_000, "{a} {b}");
}
