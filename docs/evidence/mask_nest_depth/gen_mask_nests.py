#!/usr/bin/env python3
"""A makeblastdb mask file (`-mask_data`, a Blast-db-mask-info) whose mask nests Seq-locs some levels deep.

usage: gen_mask_nests.py KIND FORMAT D > OUT

FORMAT is `text` (ASN.1 text, as `-outfmt maskinfo_asn1_text` writes it, one level a line) or `binary`
(ASN.1 binary, BER with indefinite lengths, as `-outfmt maskinfo_asn1_bin` writes it). The one mask of the
file is D Seq-locs of the KIND around an interval 5-14 of `local id 3` (the third sequence of
LOSAT/tests/seqdb/mask_nucl.fna, which makeblastdb without -parse_seqids calls so). The only type of a mask
file that holds values of its own type is the Seq-loc, through `mix` (SEQUENCE OF Seq-loc) and `equiv`
(SET OF Seq-loc); an element of either may also be stated with NCBI's pointer syntax (`:Seq-loc`, a frame
of its own in NCBI's reader).

The nests (level 0 is the outermost):

  mix        D mixes (`mix { mix { ... int { ... } } }`)
  equiv      D equivs
  alt        D Seq-locs, mix and equiv in turn (the outermost a mix)
  ptr        D mixes whose element is stated as `:Seq-loc`
  wide       D mixes, each holding an interval of `local id 3` before the mix inside it (D + 1 intervals)
  wideid     `wide` whose innermost interval is of `local id 4` (the mask has two Seq-ids)

Damaged files, the damage innermost (ERROR alone: in a `mix` nest; ERROR_NEST: in an `equiv`, `alt` or `ptr`
nest, for example `badvar_ptr`):

  badvar     a variant that Seq-loc does not have (text `mox`, binary tag [15])
  badstrand  an interval whose strand is 7, no Na-strand value
  bignum     an interval whose `from` is 99999999999 (text) or 9 bytes long (binary)
  cut        the file cut inside the innermost interval

Damaged files, the damage at the middle level M = D // 2, after the levels inside it were read whole (a `mix`
nest; `_ptr`: a `ptr` nest). These are the steps that end a level, which no damage innermost reaches:

  mid_comma      text: an interval after the element of level M, without a comma (`',' or '}' expected`)
  mid_comma_ptr  the same in a `ptr` nest (after the end of the `:Seq-loc` element)
  mid_badvar     a second element of level M of a variant that Seq-loc does not have
  mid_cut        the file cut after the element of level M
  mid_open       level M's list opened wrongly (text `[`, binary a SET OF tag where the SEQUENCE OF belongs)
  mid_cend       binary: the end of level M's list damaged (`00 05`)
  mid_vend       binary: the end of level M's variant damaged (`05 00`)
  mid_pend       binary, a `ptr` nest: the end of level M's `:Seq-loc` element damaged (`05 00`)

The inputs of docs/evidence/mask_nest_depth/ (its README), of the fixtures `n_*` of LOSAT/tests/seqdb/mask_parse
and of the tests' deep cases (`writedb/mask_info/depth_tests.rs`, whose builders write the same bytes) are made
with it.
"""
import sys

INT_TEXT = "int {{ from {from_}, to 14, {strand}id local id {id} }}"

SHAPES = ("mix", "equiv", "alt", "ptr", "wide", "wideid")
ERRORS = ("badvar", "badstrand", "bignum", "cut")
MIDDLE = ("comma", "badvar", "cut", "open", "cend", "vend", "pend")
TEXT_ONLY = ("mid_comma", "mid_comma_ptr")
BINARY_ONLY = ("mid_cend", "mid_vend", "mid_pend")


def parse(kind: str):
    """KIND as (nest shape, damage innermost or None, damage at the middle level or None)."""
    if kind in SHAPES:
        return kind, None, None
    if kind in ERRORS:
        return "mix", kind, None
    for error in ERRORS:
        for shape in ("equiv", "alt", "ptr"):
            if kind == f"{error}_{shape}":
                return shape, error, None
    if kind == "mid_pend":
        return "ptr", None, "pend"
    for damage in MIDDLE:
        if kind == f"mid_{damage}":
            return "mix", None, damage
        if kind == f"mid_{damage}_ptr" and damage in ("comma",):
            return "ptr", None, damage
    raise ValueError(kind)


def is_equiv(shape: str, level: int) -> bool:
    return shape == "equiv" or (shape == "alt" and level % 2 == 1)


def text(kind: str, d: int) -> bytes:
    if kind in BINARY_ONLY:
        raise ValueError(f"{kind} is binary only")
    shape, error, middle = parse(kind)
    from_, strand, seq = "5", "", 3
    if error == "badstrand":
        strand = "strand 7, "
    if error == "bignum":
        from_ = "99999999999"
    if shape == "wideid":
        seq = 4
    inner = INT_TEXT.format(from_=from_, strand=strand, id=seq)
    if error == "badvar":
        inner = "mox { }"
    m = d // 2
    opens, closes = [], []
    for level in range(d):
        brace = "[" if middle == "open" and level == m else "{"
        if is_equiv(shape, level):
            opens.append(f"equiv {brace}")
        elif shape == "ptr":
            opens.append(f"mix {brace} :Seq-loc")
        elif shape in ("wide", "wideid"):
            opens.append(f"mix {brace} " + INT_TEXT.format(from_="5", strand="", id=3) + ",")
        else:
            opens.append(f"mix {brace}")
        closes.append("}")
    # closes[j] closes level d - 1 - j: the list of level M + 1 is closes[d - 2 - M]
    after = d - 2 - m
    if middle == "comma":
        closes[after] += " " + INT_TEXT.format(from_="5", strand="", id=3)
    if middle == "badvar":
        closes[after] += ", mox { }"
    head = "Blast-db-mask-info ::= { algo-id 0, algo-program windowmasker, algo-options \"x\", masks { masks {\n"
    if middle == "cut":
        return (head + "\n".join(opens + [inner] + closes[: after + 1]) + "\n").encode("ascii")
    out = head + "\n".join(opens + [inner] + closes) + "\n}, more FALSE } }\n"
    if error == "cut":
        out = out[: out.rindex("to 14")]
    return out.encode("ascii")


def ber_mem(n: int, content: bytes) -> bytes:
    return bytes([0xA0 | n]) + b"\x80" + content + b"\x00\x00"


def ber_seq(*members: bytes, tag: int = 0x30) -> bytes:
    return bytes([tag]) + b"\x80" + b"".join(members) + b"\x00\x00"


def ber_int(n: int, size: int = 1) -> bytes:
    return b"\x02" + bytes([size]) + n.to_bytes(size, "big", signed=True)


def ber_interval(from_: bytes, seq: int, strand: bytes = b"") -> bytes:
    local = ber_mem(0, ber_mem(0, ber_int(seq)))
    return ber_seq(ber_mem(0, from_), ber_mem(1, ber_int(14)), *([ber_mem(2, strand)] if strand else []),
                   ber_mem(3, local))


# `:Seq-loc`: the application long tag whose bytes spell the type name (7 bits each, the last without
# the high bit), then the value with an indefinite length.
OTHER_SEQ_LOC = b"\x7f" + bytes(ord(c) | 0x80 for c in "Seq-lo") + b"c"

# A Seq-loc of the variant [15], which Seq-loc does not have.
BAD_VARIANT = ber_mem(15, b"\x05\x00")


def binary(kind: str, d: int) -> bytes:
    if kind in TEXT_ONLY:
        raise ValueError(f"{kind} is text only")
    shape, error, middle = parse(kind)
    from_, strand, seq = ber_int(5), b"", 3
    if error == "badstrand":
        strand = b"\x0a\x01\x07"
    if error == "bignum":
        from_ = b"\x02\x09" + b"\x01" * 9
    if shape == "wideid":
        seq = 4
    inner = ber_mem(3, ber_interval(from_, seq, strand))
    if error == "badvar":
        inner = BAD_VARIANT
    m = d // 2
    head, tail = [], []
    for level in range(d):
        container = b"\x31\x80" if middle == "open" and level == m else b"\x30\x80"
        if is_equiv(shape, level):
            head.append(b"\xa8\x80\x31\x80")
            tail.append(b"\x00\x00\x00\x00")
        elif shape == "ptr":
            head.append(b"\xa7\x80" + container + OTHER_SEQ_LOC + b"\x80")
            # the ends of the `:Seq-loc` element, the list and the variant
            ends = [b"\x00\x00", b"\x00\x00", b"\x00\x00"]
            if level == m and middle == "pend":
                ends[0] = b"\x05\x00"
            tail.append(b"".join(ends))
        elif shape in ("wide", "wideid"):
            head.append(b"\xa7\x80" + container + ber_mem(3, ber_interval(ber_int(5), 3)))
            tail.append(b"\x00\x00\x00\x00")
        else:
            head.append(b"\xa7\x80" + container)
            # the ends of the list and the variant
            ends = [b"\x00\x00", b"\x00\x00"]
            if level == m and middle == "cend":
                ends[0] = b"\x00\x05"
            if level == m and middle == "vend":
                ends[1] = b"\x05\x00"
            tail.append(b"".join(ends))
    # the levels' ends, innermost first; level M's follows the element of level M + 1
    ends = b"".join(reversed(tail[m + 1:]))
    if middle == "badvar":
        ends += BAD_VARIANT
    loc = b"".join(head) + inner + ends + b"".join(reversed(tail[: m + 1]))
    masks = ber_seq(ber_mem(0, ber_seq(loc)), ber_mem(1, b"\x01\x01\x00"))
    out = ber_seq(ber_mem(0, ber_int(0)), ber_mem(1, ber_int(30)), ber_mem(2, b"\x1a\x01x"), ber_mem(3, masks))
    if error == "cut":
        out = out[: len(b"".join(head)) + 40]
    if middle == "cut":
        at = out.index(loc)
        out = out[: at + len(b"".join(head)) + len(inner) + len(ends)]
    return out


KINDS = (SHAPES + ERRORS + tuple(f"{e}_{s}" for e in ERRORS for s in ("equiv", "alt", "ptr"))
         + tuple(f"mid_{m}" for m in MIDDLE) + ("mid_comma_ptr",))


def main() -> None:
    if len(sys.argv) != 4 or sys.argv[1] not in KINDS or sys.argv[2] not in ("text", "binary"):
        sys.exit(__doc__)
    kind, fmt, d = sys.argv[1], sys.argv[2], int(sys.argv[3])
    sys.stdout.buffer.write(text(kind, d) if fmt == "text" else binary(kind, d))


if __name__ == "__main__":
    main()
