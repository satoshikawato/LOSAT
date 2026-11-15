#!/usr/bin/env python3
"""Writes the inputs of `tests/xml_post_order.rs` (item K4h, review 1 M3 and M2): where NCBI BLAST+ 2.17.0's
XML streams (outfmt 5, 14, 16) post a replaced control byte (803.7) or an exception after the bytes that their 4 KiB
output buffer had written, and where the taxonomy texts of the subjects' loads go when an outfmt 5 `Iteration` is
left out with the prefetch off.

Run from LOSAT/: python3 tests/db_search/make_xpost_inputs.py. The files are committed.

- `q1.faa`: one query; `s60*.faa`: 60 `-subject` records of the same 60 residues, so the report holds 60 hits, with
  record k's title `\\x01s<k> ab z` (`c01_at<k>`), `s40 a&#1;b z` (`ent1_at40`, NCBI decodes the reference to 0x01)
  or `s60 a\\xc3\\xa9b z` (`u2_at60`, the stream throws there).
- `q1_c01first_O.faa`: a query whose title starts with 0x01 and whose letters hold an `O`, so the message of the
  report holds the title (its post after 60 hits).
- `s2_c01_u2.faa`: a post, then a throw, in one report.
- `q3_u2mid.faa`: the second of three queries has a title of two-byte UTF-8 (its `Iteration` is left out).
- `m2_qa.faa`, `m2_qd.faa`: queries against `x16_taxid_99990` (its two sequences share taxid 99990, whose lookup in
  `tests/seqdb/taxdb/offset_last_past_end` writes the offset error, and `bad_magic` writes its open's text): the
  second query's `Iteration` is left out, and the third reads the other sequence first (`m2_qa`) or the same one
  (`m2_qd`).
"""
from pathlib import Path

OUT = Path(__file__).resolve().parent / "xpost"
SEQ = b"MDIVWSTFSARIQFAQYVPPKIIRFTWHNHPSLEGWNGRQFYDIGENTGMFRLGFPKFIG"
SEQ_O = b"MDIVWSTFSARIQFAQYVOPKIIRFTWHNHPSLEGWNGRQFYDIGENTGMFRLGFPKFIG"
ST1 = (b"MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLMKYEFSSDIHPCNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHK"
       b"YPKCMCVMVRPTDQSLRHWKIMHDANFEKLKIFHYM")
ST2 = (b"SVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLMKYEFSSDIHPCNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHKYPKCM"
       b"CVMVRPTDQSLRHWKIMHDANFEKLKIFHYMVGHGRTASEG")
NOHIT = b"WWWWWWWWWWHHHHHHHHHHWWWWWWWWWW"
U2 = b"a\xc3\xa9b z"


def rec(title, seq=SEQ):
    return b">" + title + b"\n" + seq + b"\n"


def main():
    OUT.mkdir(exist_ok=True)
    files = {
        "q1.faa": rec(b"q1 ab z"),
        "q1_c01first_O.faa": rec(b"\x01q1 ab z", SEQ_O),
        "q3_u2mid.faa": rec(b"q1 ab z") + rec(b"q2 " + U2) + rec(b"q3 ab z"),
        "s60.faa": b"".join(rec(b"s%d ab z" % i) for i in range(1, 61)),
        "s60_c01_at30.faa": b"".join(
            rec(b"\x01s%d ab z" % i if i == 30 else b"s%d ab z" % i) for i in range(1, 61)),
        "s60_ent1_at40.faa": b"".join(
            rec(b"s%d a&#1;b z" % i if i == 40 else b"s%d ab z" % i) for i in range(1, 61)),
        "s60_u2_at60.faa": b"".join(
            rec(b"s%d " % i + U2 if i == 60 else b"s%d ab z" % i) for i in range(1, 61)),
        "s2_c01_u2.faa": rec(b"\x01s1 ab z") + rec(b"s2 " + U2),
        "m2_qa.faa": rec(b"q0 none", NOHIT) + rec(b"q1 " + U2, ST1) + rec(b"q2 two", ST2),
        "m2_qd.faa": rec(b"q0 none", NOHIT) + rec(b"q1 " + U2, ST2) + rec(b"q2 two", ST2) + rec(b"q3 three", ST1),
    }
    for name, data in files.items():
        (OUT / name).write_bytes(data)
    print("wrote", len(files), "files to", OUT)


if __name__ == "__main__":
    main()
