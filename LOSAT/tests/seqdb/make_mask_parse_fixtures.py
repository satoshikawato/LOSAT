#!/usr/bin/env python3
r"""Build the mask-file fixtures of LOSAT/tests/seqdb/mask_parse (the files that `makeblastdb -mask_data` reads).

NCBI BLAST+ 2.17.0 is the comparison oracle only. dustmasker, segmasker, windowmasker and
convert2blastmask write mask files (`Blast-db-mask-info`) in both ASN.1 text
(`-outfmt maskinfo_asn1_text`, extension .asnt) and ASN.1 binary (`maskinfo_asn1_bin`, .asnb) from
the same input; hand-written text files and damaged binary files cover the rest of what
makeblastdb's reader (the format guess and NCBI's serial library) accepts or refuses. For
every file, makeblastdb is run on it and what it does is written to `expected.tsv`.
The Rust tests (`LOSAT/src/writedb/mask_info/tests.rs`) read the committed files and never run NCBI.

usage: make_mask_parse_fixtures.py --ncbi-bin DIR --work DIR [--out LOSAT/tests/seqdb]

Files written into `<out>/mask_parse/`:

  dust.asnt dust.asnb      dustmasker on mask_nucl.fna            (the sequences are `local id N`)
  seg.asnt seg.asnb        segmasker on mask_prot.faa
  wm.asnt wm.asnb          windowmasker (counts made with -mk_counts) on wm_nucl.fna
  c2m.asnt c2m.asnb        convert2blastmask on lower_nucl.fna    (lower-case runs)
  <tool>_ps.asnt/.asnb     the same with -parse_seqids on ids_nucl.fna / ids_prot.faa (or the file
                           above): the Seq-ids are `local str`, `general`, `other`, `genbank`, `pdb`, ...
  ids_nucl.fna ids_prot.faa lower_nucl.fna wm_nucl.fna   the inputs of the tools (ids_lower.fna, which
                           the script also makes, is the lower-case input of the c2m_ps files)
  h_*.asnt                 hand-written text cases
  b_*.asnb                 damaged binary files (made from dust.asnb)
  expected.tsv             what makeblastdb does with each file

expected.tsv: lines starting with '#' are comments. `RUN <file> <dbtype> <fasta> <rc> <stderr> <options>` is the
exit status and the standard error (all of it) of
`makeblastdb -in <fasta> -dbtype <dbtype> -out db -mask_data <file>` (no -parse_seqids; `<fasta>` is in
LOSAT/tests/seqdb or mask_parse), run in a directory that holds the file, so that its messages
name the file without a path. The last field of RUN is the options of the algorithm as `blastdbcmd -info`
prints them (the fixed-width table, "default options used" for none; empty when no database was
written or it has several algorithms). `NODB <file>` follows a RUN with status 0 that wrote no database (makeblastdb gave up
without an error status). When the database was written, one `MASK <file> <algorithm> <oid>
<ranges>` per algorithm that `blastdbcmd -info` lists and per sequence, from
`blastdbcmd -outfmt '%o\t%m<algorithm>'` (the format of expected/mask_*.masks.tsv). Fields use backslash
escapes (`\\`, `\t`, `\n`, `\r`, `\xNN` for a byte outside printable ASCII), `\s` for a space at either end
and `\e` for an empty field.
Rerunning writes the same files (the tools and the corpus are deterministic).
"""
import argparse
import random
import re
import shutil
import subprocess
from pathlib import Path


def write_fasta(path, records, width=70):
    with open(path, "w") as f:
        for title, seq in records:
            f.write(f">{title}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i:i + width] + "\n")


def run(cmd, cwd, check=True):
    result = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True)
    if check and result.returncode != 0:
        raise SystemExit(f"{' '.join(map(str, cmd))}\n{result.stdout!r}{result.stderr!r}")
    return result


# Defline ids of ids_nucl.fna and ids_prot.faa: what -parse_seqids turns into the Seq-ids of the
# masks (local str, general, other, genbank, embl, ddbj, pdb, patent, gi, local id, swissprot, pir,
# prf, tpg).
IDS = [
    "lcl|name1 first", "gnl|mydb|tag1 second", "ref|NC_000001.1| third", "gb|AB123456.1|ABC fourth",
    "emb|X12345.1| fifth", "dbj|AB000001.1| sixth", "pdb|1ABC|A seventh", "pat|US|RE33188|1 eighth",
    "gi|12345 ninth", "gnl|db|123 tenth", "lcl|123 eleventh", "sp|P12345|NAME_HUMAN twelfth",
    "pir||A12345 thirteenth", "prf||123456 fourteenth", "tpg|BK000001.1| fifteenth",
    "ref|NP_000001.2| sixteenth",
]


def corpus(work):
    """The tool inputs (deterministic)."""
    rng = random.Random(77)

    def rand(n, alphabet="ACGT"):
        return "".join(rng.choice(alphabet) for _ in range(n))

    aa = "ACDEFGHIKLMNPQRSTVWY"
    nucl, prot = [], []
    for i, title in enumerate(IDS):
        nucl.append((title, rand(120) + "AT" * 30 + rand(80) + ("CAG" * 20 if i % 3 == 0 else "") + rand(60)))
        prot.append((title, rand(60, aa) + "Q" * 18 + "P" * 4 + rand(40, aa)
                     + ("SSSSGSSSSGSSSSGSSSS" if i % 2 == 0 else "") + rand(30, aa)))
    write_fasta(work / "ids_nucl.fna", nucl)
    write_fasta(work / "ids_prot.faa", prot)
    lower = []
    for title in IDS[:8]:
        s = rand(150)
        lower.append((title, s[:30] + s[30:60].lower() + s[60:100] + s[100:120].lower() + s[120:]))
    write_fasta(work / "ids_lower.fna", lower)
    # lower_nucl.fna: the same runs under plain ids (`local id N` for the sequences of the database)
    write_fasta(work / "lower_nucl.fna", [(f"seq{i + 1}", s) for i, (_, s) in enumerate(lower)])
    # windowmasker learns from repeats: a 40-mer repeated inside random sequence.
    rng = random.Random(5)
    repeat = rand(40)
    wm = []
    for i in range(5):
        wm.append((f"lcl|w{i}", rand(600) + repeat * 16 + rand(300) + repeat * 10 + rand(150)))
    write_fasta(work / "wm_nucl.fna", wm)


def tool_outputs(bin_dir, work, committed):
    """Mask files written by the NCBI tools, text and binary, with and without -parse_seqids.
    Returns [(file name, fasta, dbtype)]."""
    cases = []
    tests = Path(committed)
    # (name, tool command without -outfmt/-out, input fasta, input also used with -parse_seqids)
    def both(name, cmd, fasta_plain, fasta_ps, dbtype):
        for suffix, fmt in ((".asnt", "maskinfo_asn1_text"), (".asnb", "maskinfo_asn1_bin")):
            run([*cmd(fasta_plain), "-outfmt", fmt, "-out", name + suffix], work)
            cases.append((name + suffix, fasta_plain, dbtype))
            run([*cmd(fasta_ps), "-parse_seqids", "-outfmt", fmt, "-out", name + "_ps" + suffix], work)
            cases.append((name + "_ps" + suffix, fasta_plain, dbtype))

    # The sequences of mask_nucl.fna / mask_prot.faa are `local id N`; the *_ps files use ids_*.
    shutil.copyfile(tests / "mask_nucl.fna", work / "mask_nucl.fna")
    shutil.copyfile(tests / "mask_prot.faa", work / "mask_prot.faa")
    both("dust", lambda f: [bin_dir / "dustmasker", "-in", f], "mask_nucl.fna", "ids_nucl.fna", "nucl")
    both("seg", lambda f: [bin_dir / "segmasker", "-in", f], "mask_prot.faa", "ids_prot.faa", "prot")
    run([bin_dir / "windowmasker", "-in", "wm_nucl.fna", "-mk_counts", "-out", "wm.counts"], work)
    both("wm", lambda f: [bin_dir / "windowmasker", "-in", f, "-ustat", "wm.counts"],
         "wm_nucl.fna", "wm_nucl.fna", "nucl")
    both("c2m", lambda f: [bin_dir / "convert2blastmask", "-in", f, "-masking_algorithm", "lowercase",
                           "-masking_options", "-x y"],
         "lower_nucl.fna", "ids_lower.fna", "nucl")
    return cases


def obj(masks, algo_id=0, program="windowmasker", options='"x"', more="FALSE", header="Blast-db-mask-info ::= "):
    """A Blast-db-mask-info in one line; `masks` is the text inside `masks { ... }`."""
    return (f"{header}{{ algo-id {algo_id}, algo-program {program}, algo-options {options}, "
            f"masks {{ masks {{ {masks} }}, more {more} }} }}\n")


def iv(a, b, i=3, extra=""):
    return f"{{ from {a}, to {b}, {extra} id local id {i} }}"


def hand_cases():
    """Hand-written text cases: name -> content (str, latin-1)."""
    ok = obj(f"packed-int {{ {iv(5, 14)} }}")
    c = {}
    c["ok"] = ok
    c["comments"] = (
        "-- a comment before the object\n"
        "Blast-db-mask-info ::= { -- inside\n"
        "  algo-id 0, -- two dashes end it -- algo-program windowmasker,\n"
        "  algo-options \"x\", masks { masks { packed-int { " + iv(5, 14) + " } }, more FALSE } } -- after\n"
        "-- comment at the end of the file without a line end")
    c["comment_triple_dash"] = ok.replace("algo-id 0,", "algo-id 0, --- three ---")
    c["doubled_quote"] = obj(f"packed-int {{ {iv(5, 14)} }}", options='"a""b"')
    c["string_two_lines"] = obj(f"packed-int {{ {iv(5, 14)} }}", options='"ab\n   cd"')
    c["string_nonprint"] = obj(f"packed-int {{ {iv(5, 14)} }}", options='"a\tb\xe9c"')
    c["two_objects"] = obj(f"packed-int {{ {iv(5, 14)} }}") + obj(f"packed-int {{ {iv(20, 29, 4)} }}", algo_id=1, program="dust", options='"y"')
    c["two_objects_headless"] = ok + "::= " + ok.split("::= ", 1)[1]
    c["more_true_continuation"] = (
        obj(f"packed-int {{ {iv(5, 14)} }}", more="TRUE")
        + "Blast-mask-list ::= { masks { packed-int { " + iv(20, 29, 4) + " } }, more FALSE }\n")
    c["more_true_two_continuations"] = (
        obj(f"packed-int {{ {iv(5, 14)} }}", more="TRUE")
        + "Blast-mask-list ::= { masks { packed-int { " + iv(20, 29, 4) + " } }, more TRUE }\n"
        + "Blast-mask-list ::= { masks { int " + iv(1, 3, 5) + " }, more FALSE }\n")
    c["more_true_no_continuation"] = obj(f"packed-int {{ {iv(5, 14)} }}", more="TRUE")
    c["more_true_continuation_header_only"] = obj(f"packed-int {{ {iv(5, 14)} }}", more="TRUE") + "Blast-mask-list ::="
    c["more_true_wrong_type"] = obj(f"packed-int {{ {iv(5, 14)} }}", more="TRUE") + ok
    c["empty_masks"] = obj("")
    c["packed_two_ids"] = obj(f"packed-int {{ {iv(5, 14, 3)}, {iv(5, 14, 4)} }}")
    c["packed_empty"] = obj("packed-int { }")
    c["same_other_id_spacing"] = obj(
        "packed-int { { from 5, to 14, id gi 5 }, {from 20,to 29,id   gi   5 -- same\n } }")
    c["different_other_ids"] = obj("packed-int { { from 5, to 14, id gi 5 }, { from 20, to 29, id gi 6 } }")
    c["local_and_other_id"] = obj("packed-int { " + iv(5, 14, 3) + ", { from 20, to 29, id gi 6 } }")
    c["strand_plus"] = obj("packed-int { " + iv(5, 14, 3, "strand plus,") + " }")
    c["strand_minus_two"] = obj("packed-int { " + iv(20, 29, 3, "strand minus,") + ", " + iv(5, 14, 3, "strand minus,") + " }")
    c["strand_number"] = obj("packed-int { " + iv(5, 14, 3, "strand 3,") + " }")
    c["strand_bad_number"] = obj("packed-int { " + iv(5, 14, 3, "strand 7,") + " }")
    c["fuzz_from"] = obj("packed-int { { from 5, to 14, id local id 3, fuzz-from p-m 2, fuzz-to lim gt } }")
    c["fuzz_wrong_order"] = obj("packed-int { { from 5, to 14, fuzz-from p-m 2, id local id 3 } }")
    c["algo_program_number"] = obj(f"packed-int {{ {iv(5, 14)} }}", program="30")
    c["algo_program_dust"] = obj(f"packed-int {{ {iv(5, 14)} }}", program="dust")
    c["algo_program_unknown_name"] = obj(f"packed-int {{ {iv(5, 14)} }}", program="foo")
    c["algo_program_not_set"] = obj(f"packed-int {{ {iv(5, 14)} }}", program="0")
    c["algo_program_bad_number"] = obj(f"packed-int {{ {iv(5, 14)} }}", program="31")
    c["local_str_1"] = obj("packed-int { { from 5, to 14, id local str \"1\" } }")
    c["int_variant"] = obj("int " + iv(5, 14))
    c["unsupported_whole"] = obj("whole local id 3")
    c["unsupported_pnt"] = obj("pnt { point 5, id local id 3 }")
    c["unsupported_mix"] = obj("mix { int " + iv(5, 14) + " }")
    c["unsupported_null"] = obj("null NULL")
    c["unsupported_feat"] = obj("feat gibb 5")
    c["cut_inside_object"] = ok[:len(ok) // 2]
    c["cut_after_brace"] = ok + "Blast-db-mask-info ::= {"
    c["cut_after_header"] = ok + "Blast-db-mask-info ::="
    c["cut_in_type_name"] = ok + "Blast-db-"
    c["garbage_after_object"] = ok + "garbage"
    c["garbage_after_object_2"] = ok + "garbage }\n"
    c["trailing_brace"] = ok + "}\n"
    c["unknown_member"] = ok.replace("algo-id 0,", "algo-id 0, foo 1,")
    c["member_order"] = ("Blast-db-mask-info ::= { algo-program windowmasker, algo-id 0, algo-options \"x\", "
                         "masks { masks { }, more FALSE } }\n")
    c["member_duplicate"] = ok.replace("algo-id 0,", "algo-id 0, algo-id 1,")
    c["member_missing_options"] = ("Blast-db-mask-info ::= { algo-id 0, algo-program windowmasker, "
                                   "masks { masks { }, more FALSE } }\n")
    c["member_missing_more"] = ("Blast-db-mask-info ::= { algo-id 0, algo-program windowmasker, algo-options \"x\", "
                                "masks { masks { } } }\n")
    c["negative_from"] = obj("packed-int { { from -1, to 14, id local id 3 } }")
    c["from_too_big"] = obj("packed-int { { from 4294967296, to 4294967297, id local id 3 } }")
    c["to_beyond_sequence"] = obj("packed-int { { from 5, to 5000, id local id 3 } }")
    c["from_after_to"] = obj("packed-int { " + iv(14, 5) + " }")
    c["id_too_big"] = obj("packed-int { " + iv(5, 14, 2147483648) + " }")
    c["id_negative"] = obj("packed-int { " + iv(5, 14, -3) + " }")
    c["algo_id_negative"] = obj(f"packed-int {{ {iv(5, 14)} }}", algo_id=-2147483648)
    c["algo_id_overflow"] = obj(f"packed-int {{ {iv(5, 14)} }}", algo_id=2147483648)
    c["more_lowercase"] = obj(f"packed-int {{ {iv(5, 14)} }}", more="false")
    c["wrong_type_name"] = ok.replace("Blast-db-mask-info", "Blast-mask-list", 1)
    c["tab_and_cr"] = ok.replace("\n", "\r").replace(" ", "\t")
    c["crlf"] = ok.replace("\n", "\r\n")
    c["bracket_member_names"] = ok.replace("algo-id", "[algo-id]").replace("masks {", "[3] {", 1)
    c["trailing_comma"] = ok.replace("more FALSE }", "more FALSE, }")
    c["vertical_tab_between_tokens"] = ok.replace("algo-id 0", "algo-id\x0b0")
    c["nul_in_comment.asnt"] = "-- a\x00b\n" + ok
    c["utf8_in_options"] = obj(f"packed-int {{ {iv(5, 14)} }}", options='"\xc3\xa9"')
    c["bom_before_object"] = "\xef\xbb\xbf" + ok
    c["long_comment_before_object"] = "-- " + "x" * 9000 + "\n" + ok
    c["comment_mostly_punctuation"] = "-- " + "#" * 400 + "\n" + ok
    # Files that makeblastdb takes for neither format ("Unknown encoding for mask data in ...").
    c["not_asn_fasta.fna"] = ">seq1 not a mask file\nACGTACGTACGT\nACGT\n"
    c["not_asn_text.txt"] = "hello\n"
    c["not_asn_empty.txt"] = ""
    c["not_asn_comment_only.txt"] = "-- only a comment\n"
    c["not_asn_comment_only_no_newline.txt"] = "-- only a comment"
    c["not_asn_control_one.bin"] = "\x01\x01"
    c["binary_control_two.bin"] = "\x01\x02"
    return c


# A small BER encoder for the hand-made binary cases ("automatic tags": every SEQUENCE member and
# every CHOICE variant is `[n]` constructed with an indefinite length).
def ber_len(n, extra_bytes=0):
    if n < 0x80 and not extra_bytes:
        return bytes([n])
    nb = max(1, (n.bit_length() + 7) // 8) + extra_bytes
    return bytes([0x80 + nb]) + n.to_bytes(nb, "big")


def ber_mem(n, content, indef=True, tag=None):
    tag = tag if tag is not None else bytes([0xA0 | n])
    return tag + b"\x80" + content + b"\x00\x00" if indef else tag + ber_len(len(content)) + content


def ber_seq(*members, indef=True, tag=0x30):
    body = b"".join(members)
    return bytes([tag]) + (b"\x80" + body + b"\x00\x00" if indef else ber_len(len(body)) + body)


def ber_int(n, size=1):
    return b"\x02" + bytes([size]) + n.to_bytes(size, "big", signed=True)


def ber_str(s, tag=0x1A):
    s = s.encode("latin-1")
    return bytes([tag]) + ber_len(len(s)) + s


def ber_local(n):
    return ber_mem(0, ber_mem(0, ber_int(n)))


def ber_interval(a, b, i, strand=None, extra=b""):
    members = [ber_mem(0, ber_int(a)), ber_mem(1, ber_int(b))]
    if strand is not None:
        members.append(ber_mem(2, b"\x0a\x01" + bytes([strand])))
    members.append(ber_mem(3, ber_local(i)))
    return ber_seq(*members, *([extra] if extra else []))


def ber_info(locs, more=False, options="x", program=30, algo_id=0):
    masks = ber_seq(ber_mem(0, ber_seq(*locs)), ber_mem(1, b"\x01\x01" + (b"\xff" if more else b"\x00")))
    return ber_seq(ber_mem(0, ber_int(algo_id)), ber_mem(1, ber_int(program)), ber_mem(2, ber_str(options)),
                   ber_mem(3, masks))


def damaged_binary(dust):
    """Hand-made and damaged binary files (bytes); `dust` is dust.asnb."""
    c = {}
    c["cut_10"] = dust[:10]
    c["cut_100"] = dust[:100]
    c["cut_before_last_eoc_byte"] = dust[:-1]
    c["indefinite_never_closed"] = dust[:-2]
    c["length_beyond_file"] = dust.replace(b"\x1a\x1d", b"\x1a\x84\x00\x10\x00\x00", 1)
    c["length_huge"] = dust.replace(b"\x1a\x1d", b"\x1a\x88\x7f\xff\xff\xff\xff\xff\xff\xff", 1)
    c["length_nine_bytes"] = dust.replace(b"\x1a\x1d", b"\x1a\x89\x01\x00\x00\x00\x00\x00\x00\x00\x00", 1)
    c["tag_ten_bytes"] = dust.replace(b"\xa0\x80\x02\x01\x02", b"\xbf\xff\xff\xff\xff\xff\xff\xff\xff\xff\x00\x80\x02\x01\x02", 1)
    c["tag_long_form"] = dust.replace(b"\xa0\x80\x02\x01\x02", b"\xbf\x00\x80\x02\x01\x02", 1)
    c["trailing_zero"] = dust + b"\x00"
    c["trailing_sequence_tag"] = dust + b"\x30"
    c["trailing_sequence_start"] = dust + b"\x30\x80"
    c["two_objects"] = dust + dust
    c["definite_top_long_form"] = b"\x30\x82" + (len(dust) - 4).to_bytes(2, "big") + dust[2:-2]
    c["utf8_options"] = dust.replace(b"\x1a\x1d", b"\x0c\x1d", 1)
    c["bad_options_tag"] = dust.replace(b"\x1a\x1d", b"\x13\x1d", 1)
    c["wrong_first_tag"] = b"\x31" + dust[1:]
    c["swapped_members"] = dust.replace(b"\xa0\x80\x02\x01\x02\x00\x00\xa1\x80\x02\x01\x0a\x00\x00",
                                        b"\xa1\x80\x02\x01\x0a\x00\x00\xa0\x80\x02\x01\x02\x00\x00", 1)
    one = ber_interval(5, 14, 3)
    c["small_packed"] = ber_info([ber_mem(4, ber_seq(one, ber_interval(20, 29, 3)))])
    c["small_int"] = ber_info([ber_mem(3, one)])
    c["small_strand"] = ber_info([ber_mem(3, ber_interval(5, 14, 3, strand=2))])
    c["small_strand_bad"] = ber_info([ber_mem(3, ber_interval(5, 14, 3, strand=7))])
    c["small_fuzz"] = ber_info([ber_mem(3, ber_interval(5, 14, 3, extra=ber_mem(4, ber_mem(0, ber_int(3)))))])
    c["small_two_ids"] = ber_info([ber_mem(4, ber_seq(one, ber_interval(5, 14, 4)))])
    c["small_empty"] = ber_info([])
    c["small_more_true"] = ber_info([ber_mem(3, one)], more=True)
    c["small_more_true_continuation"] = ber_info([ber_mem(3, one)], more=True) + ber_seq(
        ber_mem(0, ber_seq(ber_mem(3, ber_interval(20, 29, 4)))), ber_mem(1, b"\x01\x01\x00"))
    c["small_whole"] = ber_info([ber_mem(2, ber_local(3))])
    c["small_application_integer"] = ber_info([ber_mem(3, one)]).replace(ber_mem(0, ber_int(0)), ber_mem(0, b"\x42\x01\x00"))
    c["small_definite_members"] = ber_seq(
        ber_mem(0, ber_int(0), indef=False), ber_mem(1, ber_int(30), indef=False),
        ber_mem(2, ber_str("x"), indef=False),
        ber_mem(3, ber_seq(ber_mem(0, ber_seq(), indef=False), ber_mem(1, b"\x01\x01\x00", indef=False),
                           indef=False), indef=False), indef=False)
    c["small_definite_long_member"] = ber_info([]).replace(ber_mem(0, ber_int(0)), b"\xa0\x81\x03" + ber_int(0))
    c["small_unsigned_negative"] = ber_info([ber_mem(3, ber_seq(
        ber_mem(0, b"\x02\x01\xff"), ber_mem(1, ber_int(14)), ber_mem(3, ber_local(3))))])
    c["small_unsigned_5_bytes"] = ber_info([ber_mem(3, ber_seq(
        ber_mem(0, ber_int(5)), ber_mem(1, b"\x02\x05\x00\xff\xff\xff\xff"), ber_mem(3, ber_local(3))))])
    c["small_unsigned_4_bytes_high"] = ber_info([ber_mem(3, ber_seq(
        ber_mem(0, ber_int(5)), ber_mem(1, b"\x02\x04\x80\x00\x00\x00"), ber_mem(3, ber_local(3))))])
    c["small_gi"] = ber_info([ber_mem(3, ber_seq(
        ber_mem(0, ber_int(5)), ber_mem(1, ber_int(14)), ber_mem(3, ber_mem(11, ber_int(123456789012, 5)))))])
    c["small_missing_options"] = ber_seq(ber_mem(0, ber_int(0)), ber_mem(1, ber_int(30)),
                                         ber_mem(3, ber_seq(ber_mem(0, ber_seq()), ber_mem(1, b"\x01\x01\x00"))))
    return c


def escape(s):
    """Backslash escapes; a byte outside printable ASCII is `\\xNN`."""
    out = []
    for c in s.encode("latin-1"):
        ch = chr(c)
        if ch == "\\":
            out.append("\\\\")
        elif ch == "\t":
            out.append("\\t")
        elif ch == "\n":
            out.append("\\n")
        elif ch == "\r":
            out.append("\\r")
        elif 0x20 <= c < 0x7F:
            out.append(ch)
        else:
            out.append(f"\\x{c:02x}")
    text = "".join(out)
    if text == "":
        return "\\e"
    if text.startswith(" "):
        text = "\\s" + text[1:]
    if text.endswith(" "):
        text = text[:-1] + "\\s"
    return text


def algorithm_ids(info):
    """Algorithm ids listed by `blastdbcmd -info`."""
    m = re.search(r"Available filtering algorithms applied to database sequences:\n\n.*?\n(.*?)\nVolumes:", info, re.S)
    if not m:
        return []
    return [int(line.split()[0]) for line in m.group(1).splitlines() if line.strip()]


def algorithm_options(info):
    """The options column of the algorithm table of `blastdbcmd -info`: the third column of the fixed-width
    table (13 + 40 characters before it, padded on the right; an empty option string is printed as
    "default options used")."""
    m = re.search(r"Available filtering algorithms applied to database sequences:\n\n.*?\n(.*?)\nVolumes:", info, re.S)
    if not m:
        return []
    return [line[13 + 40:].rstrip(" ") for line in m.group(1).splitlines() if line.strip()]


def expected_rows(bin_dir, work, committed, files):
    """files: [(file name, fasta, dbtype)] with the mask files in `work/mask_parse`."""
    rows = ["# RUN\tfile\tdbtype\tfasta\trc\tstderr\toptions", "# NODB\tfile",
            "# MASK\tfile\talgorithm\toid\tranges"]
    runs = work / "runs"
    for name, fasta, dbtype in files:
        run_dir = runs / name
        run_dir.mkdir(parents=True)
        shutil.copyfile(work / "mask_parse" / name, run_dir / name)
        fasta_path = work / fasta if (work / fasta).exists() else Path(committed) / fasta
        shutil.copyfile(fasta_path, run_dir / fasta)
        result = run([bin_dir / "makeblastdb", "-in", fasta, "-dbtype", dbtype, "-out", "db",
                      "-mask_data", name], run_dir, check=False)
        stderr = result.stderr.decode("latin-1")
        run_row = ["RUN", name, dbtype, fasta, str(result.returncode), escape(stderr), escape("")]
        rows.append("\t".join(run_row))
        if result.returncode != 0:
            continue
        if not any(p.name.startswith("db.") for p in run_dir.iterdir()):
            rows.append("\t".join(["NODB", name]))
            continue
        info = run([bin_dir / "blastdbcmd", "-db", "db", "-dbtype", dbtype, "-info"], run_dir).stdout.decode("latin-1")
        options = algorithm_options(info)
        if len(options) == 1:
            run_row[6] = escape(options[0])
            rows[-1] = "\t".join(run_row)
        for algo in algorithm_ids(info):
            masks = run([bin_dir / "blastdbcmd", "-db", "db", "-dbtype", dbtype, "-entry", "all",
                         "-outfmt", f"%o\t%m{algo}"], run_dir).stdout.decode()
            for line in masks.splitlines():
                oid, ranges = line.split("\t", 1)
                rows.append("\t".join(["MASK", name, str(algo), oid, escape(ranges)]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncbi-bin", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent))
    args = ap.parse_args()
    bin_dir = Path(args.ncbi_bin)
    work = Path(args.work)
    out = Path(args.out)
    if work.exists():
        raise SystemExit(f"{work} exists; give a new work directory")
    (work / "mask_parse").mkdir(parents=True)
    corpus(work)
    files = tool_outputs(bin_dir, work, out)
    for name, _, _ in files:
        shutil.move(work / name, work / "mask_parse" / name)
    for name, content in hand_cases().items():
        file = f"h_{name}" if "." in name else f"h_{name}.asnt"
        (work / "mask_parse" / file).write_bytes(content.encode("latin-1"))
        files.append((file, "mask_nucl.fna", "nucl"))
    dust = (work / "mask_parse" / "dust.asnb").read_bytes()
    for name, content in damaged_binary(dust).items():
        file = f"b_{name}.asnb"
        (work / "mask_parse" / file).write_bytes(content)
        files.append((file, "mask_nucl.fna", "nucl"))
    rows = expected_rows(bin_dir, work, out, files)

    target = out / "mask_parse"
    target.mkdir(exist_ok=True)
    keep = {"ids_nucl.fna", "ids_prot.faa", "lower_nucl.fna", "wm_nucl.fna"}
    for path in sorted((work / "mask_parse").iterdir()):
        shutil.copyfile(path, target / path.name)
    for name in keep:
        shutil.copyfile(work / name, target / name)
    (target / "expected.tsv").write_text("\n".join(rows) + "\n")


if __name__ == "__main__":
    main()
