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

  mix        D mixes (`mix { mix { ... int { ... } } }`)
  equiv      D equivs
  alt        D Seq-locs, mix and equiv in turn (the outermost a mix)
  ptr        D mixes whose element is stated as `:Seq-loc`
  wide       D mixes, each holding an interval of `local id 3` before the mix inside it (D + 1 intervals)
  wideid     `wide` whose innermost interval is of `local id 4` (the mask has two Seq-ids)
  badvar     D mixes around a variant that Seq-loc does not have (text `mox`, binary tag [15])
  badstrand  D mixes around an interval whose strand is 7, no Na-strand value
  bignum     D mixes around an interval whose `from` is 99999999999 (text) or 9 bytes long (binary)
  cut        D mixes, the file cut inside the innermost interval

The inputs of docs/evidence/mask_nest_depth/ (its README) and of the tests' deep cases are made with it.
"""
import sys

INT_TEXT = "int {{ from {from_}, to 14, {strand}id local id {id} }}"


def text(kind: str, d: int) -> bytes:
    from_, strand, seq = "5", "", 3
    if kind == "badstrand":
        strand = "strand 7, "
    if kind == "bignum":
        from_ = "99999999999"
    if kind == "wideid":
        seq = 4
    inner = INT_TEXT.format(from_=from_, strand=strand, id=seq)
    if kind == "badvar":
        inner = "mox { }"
    opens, closes = [], []
    for level in range(d):
        if kind == "equiv" or (kind == "alt" and level % 2 == 1):
            opens.append("equiv {")
        elif kind == "ptr":
            opens.append("mix { :Seq-loc")
        elif kind in ("wide", "wideid"):
            opens.append("mix { " + INT_TEXT.format(from_="5", strand="", id=3) + ",")
        else:
            opens.append("mix {")
        closes.append("}")
    body = "\n".join(opens + [inner] + closes)
    out = ("Blast-db-mask-info ::= { algo-id 0, algo-program windowmasker, algo-options \"x\", masks { masks {\n"
           + body + "\n}, more FALSE } }\n")
    if kind == "cut":
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


def binary(kind: str, d: int) -> bytes:
    from_, strand, seq = ber_int(5), b"", 3
    if kind == "badstrand":
        strand = b"\x0a\x01\x07"
    if kind == "bignum":
        from_ = b"\x02\x09" + b"\x01" * 9
    if kind == "wideid":
        seq = 4
    inner = ber_mem(3, ber_interval(from_, seq, strand))
    if kind == "badvar":
        inner = ber_mem(15, b"\x05\x00")
    head, tail = [], []
    for level in range(d):
        if kind == "equiv" or (kind == "alt" and level % 2 == 1):
            head.append(b"\xa8\x80\x31\x80")
            tail.append(b"\x00\x00\x00\x00")
        elif kind == "ptr":
            head.append(b"\xa7\x80\x30\x80" + OTHER_SEQ_LOC + b"\x80")
            tail.append(b"\x00\x00\x00\x00\x00\x00")
        elif kind in ("wide", "wideid"):
            head.append(b"\xa7\x80\x30\x80" + ber_mem(3, ber_interval(ber_int(5), 3)))
            tail.append(b"\x00\x00\x00\x00")
        else:
            head.append(b"\xa7\x80\x30\x80")
            tail.append(b"\x00\x00\x00\x00")
    loc = b"".join(head) + inner + b"".join(tail)
    masks = ber_seq(ber_mem(0, ber_seq(loc)), ber_mem(1, b"\x01\x01\x00"))
    out = ber_seq(ber_mem(0, ber_int(0)), ber_mem(1, ber_int(30)), ber_mem(2, b"\x1a\x01x"), ber_mem(3, masks))
    if kind == "cut":
        out = out[: len(b"".join(head)) + 40]
    return out


KINDS = ("mix", "equiv", "alt", "ptr", "wide", "wideid", "badvar", "badstrand", "bignum", "cut")


def main() -> None:
    if len(sys.argv) != 4 or sys.argv[1] not in KINDS or sys.argv[2] not in ("text", "binary"):
        sys.exit(__doc__)
    kind, fmt, d = sys.argv[1], sys.argv[2], int(sys.argv[3])
    sys.stdout.buffer.write(text(kind, d) if fmt == "text" else binary(kind, d))


if __name__ == "__main__":
    main()
