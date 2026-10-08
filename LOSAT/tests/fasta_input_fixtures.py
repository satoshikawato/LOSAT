#!/usr/bin/env python3
"""Frozen NCBI BLAST+ 2.17.0 outputs for the FASTA-input cases of the CFastaReader port (SF, E2h).

The rows cover how NCBI reads the FASTA files of -query and -subject (comparison oracle only):
deflines, sequence lines, line ends, text before the first defline, empty records, standard
input, query batches (BATCH_SIZE), -lcase_masking, for BLASTN (megablast, blastn, dc-megablast,
blastn-short), TBLASTX, TBLASTN (protein query, nucleotide subject) and BLASTP, with outfmt 0, 6
and 7 and, on a core subset, -num_threads 2 and 4. Rows of kind `losat-rejection` (first lines
that NCBI would parse as a Seq-id and fetch over the network) are never run on NCBI: LOSAT must
reject them, and the port fills their expectation (hash columns `pending`).

- generate: writes the inputs to LOSAT/tests/fixtures/fasta_input/ (deterministic; the files are
  committed, so later runs never regenerate them).
- freeze --ncbi-bin DIR --out DIR --jobs N: runs NCBI BLAST+ from LOSAT/ for every row of kind
  `ncbi` (twice: once with stdout and stderr separate, once with `2>&1`, which records the order
  of the warnings and the report), writes <row>.out/.err/.merged/.rc into --out and the hashes
  into LOSAT/tests/fasta_input_fixtures/<program>.tsv (row id, rc, stdout, stderr and combined
  hashes only, in a compact form; the row definitions stay in this script). `--rows` limits the rows by program (the
  other rows of the manifests stay as they are); `--verify` runs NCBI again and compares with the
  manifest. -num_threads is passed to LOSAT only (NCBI ignores it with -subject and prints a
  warning that LOSAT omits: PD-LOSAT-CLI-NONSEARCH-DIFFERENCES), as the other fixtures do.
- refresh: rewrites the manifests from their existing hashes and the current row definitions.
- check --losat BIN --out TSV --jobs N: runs `LOSAT <program> <argv>` on the same rows from
  LOSAT/ and compares stdout, stderr, the exit code and the combined `2>&1` stream with the
  manifest. One line per row: `same`, `differs` (with the first differing line, when the frozen
  files of --expected are present), `rejects` (LOSAT's explicit "not supported by LOSAT"
  rejection where NCBI reads the input: the port must turn it into `same`), `pending` (the
  expectation is not filled yet) or `error`. A last column `group` marks the rows whose option LOSAT rejects today
  (`-lcase_masking` of TBLASTX/BLASTP; port or explicit rejection is decided in SFb) and the summary counts them apart.

NCBI is run without ~/.ncbirc, with BLASTDB, BATCH_SIZE, CHUNK_SIZE, OVERLAP_CHUNK_SIZE and
BLASTINPUT_GEN_DELTA_SEQ unset (a row sets them in `env`), and no byte of its output is
normalised. The first line of every NCBI input is checked first: NCBI tries a first line that
starts with a letter or digit as a Seq-id and contacts the network, so only files that start
with `>`, a blank line, a non-alphanumeric byte or a letters-only line are run (and `unshare -rn`
cuts the network, where it works).

Usage:
  fasta_input_fixtures.py generate
  fasta_input_fixtures.py refresh
  fasta_input_fixtures.py freeze --ncbi-bin DIR --out DIR [--jobs N] [--rows blastn,tblastx] [--verify]
  fasta_input_fixtures.py check --losat BIN --out TSV [--jobs N] [--rows blastn] [--expected DIR]
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import os
import random
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass, fields
from pathlib import Path

ENGINE = Path(__file__).resolve().parents[1]
FIXTURES = ENGINE / "tests/fixtures/fasta_input"
REL = "tests/fixtures/fasta_input"
MANIFEST_DIR = ENGINE / "tests/fasta_input_fixtures"      # one <program>.tsv per program
EMPTY_SHA256 = hashlib.sha256(b"").hexdigest()
OPTION_NOTE = "option rejected today; port or keep as an explicit rejection is decided in SFb"
TIMEOUT = 180
# NCBI reads these and each one changes the batches or the report.
REPORT_ENV = ("BL2SEQ_LEGACY", "CTOOLKIT_COMPATIBLE", "OLD_FSC", "BATCH_SIZE", "CHUNK_SIZE", "ADAPTIVE_CBS",
              "OVERLAP_CHUNK_SIZE", "PRE_FETCH_SEQS_LIMIT", "BLASTDB", "BLASTINPUT_GEN_DELTA_SEQ")


@dataclass
class Row:
    row_id: str
    program: str
    task: str
    query: str          # path relative to LOSAT/, "-" for standard input, "" when -query is omitted
    subject: str
    args: str
    outfmt: str
    threads: str
    stdin: str          # "", "pipe:<path>", "file:<path>" or "empty" (an empty pipe)
    env: str
    kind: str           # ncbi | losat-rejection
    tier: str = ""      # core (all BLASTN tasks, outfmt 0/6/7), full (fewer), extra (stdin, batches, ...), seqid
    rc: str = "pending"
    stdout_sha256: str = "pending"
    stderr_sha256: str = "pending"
    combined_sha256: str = "pending"
    note: str = ""


FIELDS = [f.name for f in fields(Row)]
EXPECTED = ("rc", "stdout_sha256", "stderr_sha256", "combined_sha256")


# ---------------------------------------------------------------------------------------------
# Input material
# ---------------------------------------------------------------------------------------------
AA = "ACDEFGHIKLMNPQRSTVWY"
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS: dict[str, list[str]] = {}
for _i, _aa in enumerate(_AAS):
    CODONS.setdefault(_aa, []).append(_BASES[_i // 16] + _BASES[(_i // 4) % 4] + _BASES[_i % 4])


def make_base() -> dict:
    """Three protein fragments (60 aa) with a coding sequence each (180 nt). Queries are the
    fragments; subjects embed them between random flanks, so every program finds three hits."""
    rng = random.Random(20261008)
    prot = ["M" + "".join(rng.choice(AA) for _ in range(59)) for _ in range(3)]
    coding = ["".join(rng.choice(CODONS[a]) for a in p) for p in prot]

    def nuc(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    def aas(n):
        return "".join(rng.choice(AA) for _ in range(n))

    base = {("nuc", "q"): [c.encode() for c in coding], ("prot", "q"): [p.encode() for p in prot]}
    base[("nuc", "s")] = [(nuc(60) + c + nuc(60)).encode() for c in coding]
    base[("prot", "s")] = [(aas(30) + p + aas(30)).encode() for p in prot]
    # four more queries for the batch rows (the fragments again, other titles)
    base[("nuc", "bq")] = base[("nuc", "q")] + [c.encode() for c in coding[::-1]] + [coding[0].encode()]
    base[("prot", "bq")] = base[("prot", "q")] + [p.encode() for p in prot[::-1]] + [prot[0].encode()]
    return base


BASE = make_base()
WRAP = {"nuc": 60, "prot": 30}
EXT = {"nuc": "fa", "prot": "faa"}


def wrap_lines(seq: bytes, width: int) -> list[bytes]:
    return [seq[i:i + width] for i in range(0, len(seq), width)]


class Ctx:
    """The records of one clean file (role q or s, molecule nuc or prot) that a variant edits.
    recs is a list of [defline bytes without '>' (None: no defline), list of sequence lines]."""

    def __init__(self, kind: str, role: str):
        self.kind, self.role = kind, role
        self.p = "q" if role in ("q", "bq") else "s"
        word = {"q": ("first query", "second query", "third query"),
                "s": ("first subject", "second subject", "third subject")}[self.p]
        seqs = BASE[(kind, role)]
        titles = [f"{self.p}{i + 1} {word[i % 3]}" if i < 3 else f"{self.p}{i + 1} later query {i + 1}"
                  for i in range(len(seqs))]
        self.recs = [[t.encode(), wrap_lines(s, WRAP[kind])] for t, s in zip(titles, seqs)]

    def render(self, eol: bytes = b"\n", final: bool = True, recs=None) -> bytes:
        out = []
        for defline, lines in (self.recs if recs is None else recs):
            if defline is not None:
                out.append(b">" + defline)
            out.extend(lines)
        data = eol.join(out)
        return data + eol if final else data

    def seq(self, i: int) -> bytes:
        return b"".join(self.recs[i][1])


def ins(line: bytes, pos: int, what: bytes) -> bytes:
    return line[:pos] + what + line[pos:]


VARIANTS: dict[str, dict] = {}


def V(name: str, tier: str = "full", kinds=("nuc", "prot"), roles=("q", "s"), lcase: bool = False,
      threads: bool = False, note: str = ""):
    def deco(fn):
        assert name not in VARIANTS, name
        VARIANTS[name] = dict(fn=fn, tier=tier, kinds=kinds, roles=roles, lcase=lcase, threads=threads, note=note)
        return fn
    return deco


# --- deflines ----------------------------------------------------------------------------
@V("def_empty", "core", threads=True)
def _(c): c.recs[1][0] = b""; return c.render()
@V("def_ws_only", "core")
def _(c): c.recs[1][0] = b"   "; return c.render()
@V("def_lead_space", "core")
def _(c): c.recs[0][0] = b"   lead title"; return c.render()
@V("def_lead_tab")
def _(c): c.recs[0][0] = b"\t" + c.recs[0][0]; return c.render()
@V("def_tab_title", "core")
def _(c): c.recs[0][0] = c.recs[0][0].replace(b" ", b"\ttab ", 1); return c.render()
@V("def_ctl_title")
def _(c): c.recs[1][0] = c.recs[1][0].replace(b" ", b" ctl\x01tail ", 1); return c.render()
@V("def_ctl_first", "core")
def _(c): c.recs[0][0] = b"\x01" + c.recs[0][0]; return c.render()
@V("def_nul")
def _(c): c.recs[0][0] = c.recs[0][0].replace(b" ", b"\x00zz title ", 1); return c.render()
@V("def_utf8", "core", threads=True)
def _(c):
    c.recs[0][0] = f"{c.p}1é first café 日本語 title".encode()
    c.recs[1][0] = f"{c.p}2 second café 日本語".encode()
    return c.render()
@V("def_utf8_id_space")
def _(c): c.recs[0][0] = "éè first ü".encode(); return c.render()
@V("def_utf8_wrap3", "core")
def _(c): c.recs[0][0] = f"{c.p}1 ".encode() + "日".encode() * 50 + b" end"; return c.render()
@V("def_utf8_ends", "full")
def _(c): c.recs[0][0] = c.recs[0][0] + " café".encode(); c.recs[1][0] = c.recs[1][0] + " 日".encode(); return c.render()
@V("def_nbsp_u3000")
def _(c): c.recs[0][0] = f"{c.p}1 first　title".encode(); return c.render()
@V("def_nonutf8", "core")
def _(c):
    c.recs[0][0] = f"{c.p}1".encode() + b"\xff\xfe first \xe9"
    c.recs[1][0] = f"{c.p}2 second".encode() + b" \x93quoted\x94"
    return c.render()
@V("def_nonutf8_wrap", "core")
def _(c): c.recs[0][0] = f"{c.p}1 ".encode() + b"\xe9\xff\xfe" * 40 + b" end"; return c.render()
@V("def_gt_inside")
def _(c): c.recs[0][0] = c.recs[0][0] + b" a>b >c"; return c.render()
@V("def_mod_kv", "core")
def _(c):
    c.recs[0][0] = c.recs[0][0] + b" [organism=Escherichia coli] [strain=K-12 substr. MG1655]"
    c.recs[1][0] = c.recs[1][0] + b" [key=value]"
    return c.render()
@V("def_mod_odd")
def _(c):
    c.recs[0][0] = c.recs[0][0] + b" [organism=E. coli first"
    c.recs[1][0] = c.recs[1][0] + b" [] [=x] [organism] [a=b"
    c.recs[2][0] = b"[organism=Foo bar] " + c.recs[2][0]
    return c.render()


def title_of(length: int, tail: bytes = b"x") -> bytes:
    words, n = [], 0
    while n < length:
        words.append(b"word%03d" % (len(words) % 1000))
        n += 8
    return b" ".join(words)[:length - len(tail)] + tail


@V("def_title_1001", "core")
def _(c): c.recs[0][0] = f"{c.p}1 ".encode() + title_of(1001); return c.render()
@V("def_long_id")
def _(c): c.recs[0][0] = f"{c.p}".encode() + b"x" * 120 + b" long id title"; return c.render()
@V("def_seqid_like", "core")
def _(c):
    c.recs[0][0] = b"gb|AB123456.1| first"
    c.recs[1][0] = b"lcl|foo second"
    c.recs[2][0] = b"sp|P01308|INS_HUMAN third"
    return c.render()
@V("def_html")
def _(c): c.recs[0][0] = c.recs[0][0] + b" A &amp; B &lt;x&gt; &#65; &quot;q&quot; &unknown; &#x41;"; return c.render()
@V("def_trailing_punct")
def _(c): c.recs[0][0] = c.recs[0][0] + b" title, with trailing;~ ."; c.recs[1][0] = c.recs[1][0] + b" ,;~"; return c.render()
@V("def_prefix")
def _(c): c.recs[0][0] = b"TPA: " + c.recs[0][0]; c.recs[1][0] = c.recs[1][0] + b" MAG: x"; return c.render()
@V("def_spaces_compress")
def _(c): c.recs[0][0] = c.recs[0][0] + b"   many   spaces ( in ) title ,, here"; return c.render()
@V("def_trailing_ws")
def _(c): c.recs[0][0] = c.recs[0][0] + b"   "; c.recs[1][0] = c.recs[1][0] + b" \t "; return c.render()
@V("def_gt_only", "core", note="a lone > as the last line")
def _(c): c.recs.append([b"", []]); return c.render()
@V("def_gt_only_nofinal", "core")
def _(c): c.recs.append([b"", []]); return c.render()[:-1]


# 20 nucleotides / 50 amino acids at the end of the title (CFastaReader's warning at the end of the record)
@V("def_title_nuc20", "core", kinds=("nuc",), threads=True)
def _(c): c.recs[0][0] = c.recs[0][0] + b" ACGTACGTACGTACGTACGT"; c.recs[2][0] += b" acgtacgtacgtacgtacgtac"; return c.render()
@V("def_title_nuc20_ws", "core", kinds=("nuc",))
def _(c):
    c.recs[0][0] += b" ACGTACGTACGTACGTACGT  "
    c.recs[1][0] += b" ACGTACGTACGTACGTACGT\t"
    return c.render()
@V("def_title_nuc19", kinds=("nuc",))
def _(c): c.recs[0][0] += b" ACGTACGTACGTACGTACG"; return c.render()
@V("def_title_nuc_whole", kinds=("nuc",))
def _(c): c.recs[0][0] = b"ACGTACGTACGTACGTACGT"; c.recs[1][0] = b"ACGTACGTACGTACGTACGTA"; return c.render()
@V("def_title_nuc_n", kinds=("nuc",))
def _(c): c.recs[0][0] += b" ACGTACGTNCGTACGTACGTACGT"; c.recs[1][0] += b" UUUUUUUUUUUUUUUUUUUUUUU"; return c.render()
@V("def_title_aa50", "core", kinds=("prot",))
def _(c): c.recs[0][0] += b" " + b"ACDEFGHIKLMNPQRSTVWY" * 2 + b"ACDEFGHIKL"; return c.render()
@V("def_title_aa49_51", kinds=("prot",))
def _(c):
    c.recs[0][0] += b" " + (b"ACDEFGHIKLMNPQRSTVWY" * 3)[:49]
    c.recs[1][0] += b" " + (b"ACDEFGHIKLMNPQRSTVWY" * 3)[:51]
    return c.render()
@V("def_title_aa50_ws", kinds=("prot",))
def _(c): c.recs[0][0] += b" " + (b"ACDEFGHIKLMNPQRSTVWY" * 3)[:50] + b" "; return c.render()


# --- gap lines (">?") --------------------------------------------------------------------
def gap_at(c, line: bytes, rec: int = 1, pos: int = 1):
    c.recs[rec][1].insert(pos, line)
    return c.render()


@V("gap_100", "core")
def _(c): return gap_at(c, b">?100")
@V("gap_unk100", "core")
def _(c): return gap_at(c, b">?unk100")
@V("gap_bare", "core")
def _(c): return gap_at(c, b">?")
@V("gap_abc", "core")
def _(c): return gap_at(c, b">?abc")
@V("gap_underscore", "core")
def _(c): return gap_at(c, b">?_x pseudo")
@V("gap_space")
def _(c): return gap_at(c, b">? 100")
@V("gap_mods")
def _(c): return gap_at(c, b">?100 [gap-type=within scaffold] [linkage-evidence=paired-ends]")
@V("gap_mods_bad")
def _(c): return gap_at(c, b">?100 [gap-type=nonsense] [linkage-evidence=nope]")
@V("gap_lead_space")
def _(c): return gap_at(c, b" >?7")
@V("gap_start")
def _(c): c.recs.insert(0, [b"?100", []]); return c.render()
@V("gap_start_record")
def _(c): c.recs[0][1].insert(0, b">?30"); return c.render()
@V("gap_end")
def _(c): c.recs[1][1].append(b">?50"); return c.render()
@V("gap_two")
def _(c): c.recs[1][1].insert(1, b">?20"); c.recs[1][1].insert(3, b">?unk30"); return c.render()
@V("gap_lower", lcase=True)
def _(c):
    c.recs[1][1][0] = c.recs[1][1][0][:30] + c.recs[1][1][0][30:].lower()
    c.recs[1][1].insert(1, b">?40")
    c.recs[1][1][2] = c.recs[1][1][2][:20].lower() + c.recs[1][1][2][20:]
    return c.render()


# --- sequence lines ----------------------------------------------------------------------
def edit(c, rec: int, line: int, fn):
    c.recs[rec][1][line] = fn(c.recs[rec][1][line])
    return c.render()


@V("seq_digit", "core", threads=True)
def _(c): return edit(c, 1, 1, lambda b: ins(b, 10, b"1234"))
@V("seq_digit_first")
def _(c): return edit(c, 0, 0, lambda b: ins(b, 8, b"7"))
@V("seq_hyphen", "core")
def _(c): return edit(c, 1, 1, lambda b: ins(b, 12, b"-") if c.kind == "nuc" else ins(b, 12, b"--"))
@V("seq_hyphen_first")
def _(c): return edit(c, 0, 0, lambda b: ins(ins(b, 6, b"-"), 20, b"---"))
@V("seq_hyphen_lines", "core")
def _(c):
    for rec, line in ((0, 1), (1, 0), (1, 1)):
        c.recs[rec][1][line] = ins(c.recs[rec][1][line], 5, b"-")
    return c.render()
@V("seq_hyphen_first_only", note="first data line of the record is only hyphens (CheckDataLine)")
def _(c): c.recs[1][1].insert(0, b"----------"); return c.render()
@V("seq_semicolon", "core")
def _(c): return edit(c, 1, 1, lambda b: ins(b, 15, b";comment text"))
@V("seq_semicolon_short", note="AC;rest as the first data line (CheckDataLine)")
def _(c): c.recs[1][1].insert(0, b"AC;rest of line"); return c.render()
@V("seq_ws_inside", "core")
def _(c): return edit(c, 1, 1, lambda b: ins(ins(b, 10, b" "), 25, b"\t"))
@V("seq_ws_lead_trail")
def _(c): return edit(c, 1, 1, lambda b: b"  \t" + b + b" \t ")
@V("seq_ws_only_line")
def _(c): c.recs[1][1].insert(1, b"   \t  "); return c.render()
@V("seq_blank_lines", "core")
def _(c): c.recs[1][1].insert(1, b""); c.recs[1][1].insert(3, b""); c.recs[0][1].append(b""); return c.render()
@V("seq_nbsp", "core", threads=True)
def _(c): return edit(c, 1, 1, lambda b: ins(b, 10, b"\xc2\xa0"))
@V("seq_nbsp_first")
def _(c): return edit(c, 0, 0, lambda b: ins(b, 10, b"\xc2\xa0"))
@V("seq_nonascii", "core")
def _(c): return edit(c, 1, 1, lambda b: ins(ins(b, 10, "é".encode()), 30, b"\xff"))
@V("seq_star", "core")
def _(c): return edit(c, 1, 1, lambda b: ins(b, 10, b"*"))
@V("seq_star_end")
def _(c): return edit(c, 1, len(c.recs[1][1]) - 1, lambda b: b + b"*")
@V("seq_letters", "core", note="nucleotide: E F I J L O P Q Z X are invalid; protein: B J O U X Z are residues")
def _(c):
    extra = b"EFIJLOPQZX" if c.kind == "nuc" else b"BJOUXZ"
    return edit(c, 1, 1, lambda b: ins(b, 10, extra))
@V("seq_iupac", kinds=("nuc",))
def _(c): return edit(c, 1, 1, lambda b: ins(b, 10, b"RYKMSWBDHVN"))
@V("seq_u", "core", kinds=("nuc",), lcase=True, threads=True)
def _(c):
    c.recs[1][1] = [b.replace(b"T", b"U") for b in c.recs[1][1]]
    c.recs[0][1][0] = c.recs[0][1][0].replace(b"T", b"u")
    return c.render()
@V("seq_u_prot", kinds=("prot",))
def _(c): return edit(c, 1, 0, lambda b: ins(b, 10, b"UOU"))
@V("seq_lower_all", "core", lcase=True)
def _(c): c.recs[1][1] = [b.lower() for b in c.recs[1][1]]; return c.render()
@V("seq_lower_mixed", lcase=True)
def _(c): c.recs[0][1] = [b.lower() if i % 2 else b for i, b in enumerate(c.recs[0][1])]; return c.render()
@V("seq_lower_island", "core", lcase=True, threads=True)
def _(c): return edit(c, 1, 1, lambda b: b[:10] + b[10:25].lower() + b[25:])
@V("seq_lower_island_digit", "core", lcase=True)
def _(c): return edit(c, 1, 1, lambda b: b[:10] + b[10:16].lower() + b"7" + b[16:25].lower() + b[25:])
@V("seq_lower_island_strip", lcase=True)
def _(c): return edit(c, 1, 1, lambda b: b[:10] + b[10:16].lower() + b"-" + b[16:20].lower() + b" " + b[20:25].lower() + b";cmt")
@V("seq_lower_lines", lcase=True)
def _(c):
    c.recs[1][1][0] = c.recs[1][1][0][:40] + c.recs[1][1][0][40:].lower()
    c.recs[1][1][1] = c.recs[1][1][1][:15].lower() + c.recs[1][1][1][15:]
    return c.render()
@V("seq_lower_ends", lcase=True)
def _(c):
    c.recs[0][1][0] = c.recs[0][1][0][:12].lower() + c.recs[0][1][0][12:]
    c.recs[2][1][-1] = c.recs[2][1][-1][:-10] + c.recs[2][1][-1][-10:].lower()
    return c.render()
@V("seq_lower_ambig", kinds=("nuc",), lcase=True)
def _(c): return edit(c, 1, 1, lambda b: b[:10] + b"nnnnnrykm" + b[10:20].lower() + b[20:])
@V("seq_lower_star", kinds=("prot",), lcase=True)
def _(c): return edit(c, 1, 0, lambda b: b[:10] + b[10:20].lower() + b"*" + b[20:25].lower() + b[25:])
@V("seq_ambig_first", "core", note="the first data line is mostly N or X")
def _(c):
    amb = b"N" if c.kind == "nuc" else b"X"
    c.recs[0][1].insert(0, amb * 40 + b"ACGT" * 5)
    return c.render()
@V("seq_ambig_all_first_line")
def _(c):
    amb = b"N" if c.kind == "nuc" else b"X"
    c.recs[1][1].insert(0, amb * 60)
    return c.render()
# punctuation: CheckDataLine fails (words are only invalid residues in nucleotide text, and valid letters in protein text)
BAD_TEXT = {"nuc": b"$%&()+= <>?@ [] {} ^~", "prot": b"$%&()+= <>?@ [] {} ^~"}
WORDS = b"This is not a sequence line at all"


@V("check_text", "core", threads=True)
def _(c): c.recs[1][1].insert(0, BAD_TEXT[c.kind]); return c.render()
@V("check_text_first", "core")
def _(c): c.recs[0][1].insert(0, BAD_TEXT[c.kind]); return c.render()
@V("check_text_later_line", note="a bad line after a residue line is only a warning")
def _(c): c.recs[1][1].insert(1, BAD_TEXT[c.kind]); return c.render()
@V("check_words", note="words: invalid residues (warning) in nucleotide text, valid letters in protein text")
def _(c): c.recs[1][1].insert(0, WORDS); return c.render()
@V("check_digits_first")
def _(c): c.recs[1][1].insert(0, b"123456789012345"); return c.render()
@V("check_genbank_numbers")
def _(c): c.recs[1][1][0] = b"        1 " + b" ".join(wrap_lines(c.recs[1][1][0][:30], 10)); return c.render()
@V("check_third_record")
def _(c): c.recs[2][1].insert(0, BAD_TEXT[c.kind]); return c.render()
@V("seq_comment_lines", "core")
def _(c):
    c.recs[1][1].insert(1, b"; a comment line")
    c.recs[1][1].insert(3, b"# another comment")
    c.recs[0][1].insert(1, b"! third kind")
    c.recs[2][1].append(b";trailing comment")
    return c.render()
@V("seq_comment_between_records")
def _(c): c.recs[0][1].extend([b"#c1", b";c2", b"", b"!c3"]); return c.render()
@V("seq_gt_not_col0")
def _(c): c.recs[1][1].insert(1, b" >q9 not a defline"); return c.render()
@V("seq_gt_mid_line")
def _(c): return edit(c, 1, 1, lambda b: ins(b, 10, b">"))
@V("seq_nul")
def _(c): return edit(c, 1, 1, lambda b: ins(b, 10, b"\x00"))
@V("seq_ctrl")
def _(c): return edit(c, 1, 1, lambda b: ins(ins(b, 10, b"\x01"), 20, b"\x1a"))
@V("seq_long_line", kinds=("nuc",), roles=("q", "s"))
def _(c): c.recs[1][1] = [b"".join(c.recs[1][1])]; return c.render()
@V("seq_one_char_lines", kinds=("nuc",), roles=("q",))
def _(c): c.recs[1][1] = wrap_lines(c.seq(1), 1); return c.render()
@V("seq_lines_70_80")
def _(c):
    c.recs[0][1] = wrap_lines(c.seq(0), 70)
    c.recs[1][1] = wrap_lines(c.seq(1), 80)
    return c.render()


# --- no defline, text before the first defline --------------------------------------------
@V("nodef_seq_only", "core", note="a letters-only first line is read as FASTA data (no Seq-id)", threads=True)
def _(c): c.recs = [[None, c.recs[0][1]]]; return c.render()
@V("nodef_then_defline", "core")
def _(c): c.recs[0][0] = None; return c.render()
@V("nodef_two_chunks")
def _(c): c.recs[0][0] = None; c.recs[1][0] = None; return c.render()
@V("nodef_ws_lead")
def _(c): c.recs = [[None, [b"  " + c.recs[0][1][0]] + c.recs[0][1][1:]]]; return c.render()
@V("nodef_blank_first")
def _(c): c.recs = [[None, [b"", b""] + c.recs[0][1]]]; return c.render()
@V("nodef_comment_first")
def _(c): c.recs = [[None, [b"; comment"] + c.recs[0][1]]]; return c.render()
@V("nodef_digit_inside", kinds=("nuc", "prot"))
def _(c): c.recs = [[None, [c.recs[0][1][0], ins(c.recs[0][1][1], 5, b"1")] + c.recs[0][1][2:]]]; return c.render()


def pre(c, lines: list[bytes]) -> bytes:
    return b"\n".join(lines) + b"\n" + c.render()


@V("pre_blank", "core")
def _(c): return pre(c, [b"", b""])
@V("pre_ws_line")
def _(c): return pre(c, [b"   \t "])
@V("pre_comment_hash", "core")
def _(c): return pre(c, [b"# a comment"])
@V("pre_comment_semicolon", "core")
def _(c): return pre(c, [b"; a comment"])
@V("pre_comment_excl")
def _(c): return pre(c, [b"! a comment"])
@V("pre_comment_mixed", "core")
def _(c): return pre(c, [b"", b"# one", b"   ", b"; two", b"! three", b""])
@V("pre_bom", "core")
def _(c): return b"\xef\xbb\xbf" + c.render()
@V("pre_bom_nodef")
def _(c): c.recs = [[None, c.recs[0][1]]]; return b"\xef\xbb\xbf" + c.render()
@V("pre_bom_crlf")
def _(c): return b"\xef\xbb\xbf" + c.render(eol=b"\r\n")
@V("pre_lead_space_gt")
def _(c): return b" " + c.render()
@V("pre_text")
def _(c): return pre(c, [b"** free text line before the first record"])


# --- empty records and files -------------------------------------------------------------
@V("rec_empty_first", "core")
def _(c): c.recs[0][1] = []; return c.render()
@V("rec_empty_mid", "core", threads=True)
def _(c): c.recs[1][1] = []; return c.render()
@V("rec_empty_last", "core")
def _(c): c.recs[2][1] = []; return c.render()
@V("rec_empty_two", "core")
def _(c): c.recs[0][1] = []; c.recs[2][1] = []; return c.render()
@V("rec_empty_all", "core")
def _(c):
    for r in c.recs: r[1] = []
    return c.render()
@V("rec_empty_all_nofinal")
def _(c):
    for r in c.recs: r[1] = []
    return c.render(final=False)
@V("rec_empty_blank_lines")
def _(c): c.recs[1][1] = [b"", b"  ", b""]; return c.render()
@V("rec_empty_comment_only")
def _(c): c.recs[1][1] = [b"; only a comment", b"# and another"]; return c.render()
@V("rec_empty_all_bad")
def _(c): c.recs[1][1] = [b"1234567890", b"5555"]; return c.render()
@V("rec_empty_star_only", kinds=("prot",))
def _(c): c.recs[1][1] = [b"***"]; return c.render()
@V("rec_empty_untitled")
def _(c): c.recs[1] = [b"", []]; c.recs[2] = [b"", []]; return c.render()
@V("rec_empty_adjacent_gt")
def _(c): c.recs.insert(1, [b"e1 empty one", []]); c.recs.insert(2, [b"e2 empty two", []]); return c.render()
@V("rec_empty_long_title", note="the empty-record warning cuts a long ID and title")
def _(c): c.recs[1] = [f"{c.p}2 ".encode() + b"long title " * 5, []]; return c.render()
@V("rec_empty_n_only", kinds=("nuc",), note="a record of only N (all masked, not empty)")
def _(c): c.recs[1][1] = [b"N" * 60]; return c.render()


def solo(name, tier, data, note="", threads=False, roles=("q", "s")):
    @V(name, tier, kinds=("nuc", "prot"), roles=roles, threads=threads, note=note)
    def _(c, data=data):
        return data
    return _


solo("file_empty", "core", b"", threads=True)
solo("file_ws_only", "core", b"  \n\t\n \r\n\n")
solo("file_blank_lines", "full", b"\n\n\n")
solo("file_comment_only", "core", b"# one\n; two\n\n! three\n")
solo("file_gt_only", "core", b">\n")
solo("file_gt_only_nofinal", "full", b">")
solo("file_gt_title_only", "full", b">q1 title only\n")
solo("file_two_gt_only", "full", b">\n>\n")
solo("file_ff_only", "full", b"\x0c\n")
solo("file_nul_only", "full", b"\x00\n")
solo("file_nbsp_only", "full", b"\xc2\xa0\n")
solo("file_ctrlz_only", "full", b"\x1a\n")


# --- line ends ---------------------------------------------------------------------------
@V("eol_crlf", "core")
def _(c): return c.render(eol=b"\r\n")
@V("eol_cr", "core", threads=True)
def _(c): return c.render(eol=b"\r")
@V("eol_lf_nofinal", "core")
def _(c): return c.render(final=False)
@V("eol_crlf_nofinal", "core")
def _(c): return c.render(eol=b"\r\n", final=False)
@V("eol_cr_nofinal", "core")
def _(c): return c.render(eol=b"\r", final=False)


def mixed(c, eols: list[bytes], final: bytes = b"\n") -> bytes:
    lines = []
    for defline, ls in c.recs:
        if defline is not None:
            lines.append(b">" + defline)
        lines.extend(ls)
    out = b""
    for i, line in enumerate(lines):
        out += line + (eols[i % len(eols)] if i < len(lines) - 1 else final)
    return out


@V("eol_mixed_lf_crlf_cr", "core")
def _(c): return mixed(c, [b"\n", b"\r\n", b"\r"])
@V("eol_mixed_cr_lf", "core")
def _(c): return mixed(c, [b"\r", b"\n"])
@V("eol_mixed_crlf_lf")
def _(c): return mixed(c, [b"\r\n", b"\n"])
@V("eol_mixed_lf_cr")
def _(c): return mixed(c, [b"\n", b"\r"])
@V("eol_mixed_crlf_cr")
def _(c): return mixed(c, [b"\r\n", b"\r"])
@V("eol_lfcr")
def _(c): return c.render(eol=b"\n\r")
@V("eol_crcrlf")
def _(c): return c.render(eol=b"\r\r\n")
@V("eol_crlf_blank_tail")
def _(c): return c.render(eol=b"\r\n") + b"\r\n\r\n"
@V("eol_cr_blank_tail")
def _(c): return c.render(eol=b"\r") + b"\r\r"
@V("eol_lf_blank_tail")
def _(c): return c.render() + b"\n\n\n"
@V("eol_cr_comment_tail")
def _(c): return c.render(eol=b"\r") + b"#tail comment\r"
@V("eol_lf_comment_tail_nofinal")
def _(c): return c.render() + b"# tail comment"
@V("eol_cr_hyphen_lines", "core", note="line numbers of the hyphen warnings in a CR-only file")
def _(c):
    for rec, line in ((0, 1), (1, 0), (2, 1)):
        c.recs[rec][1][line] = ins(c.recs[rec][1][line], 5, b"-")
    return c.render(eol=b"\r")
@V("eol_crlf_hyphen_lines")
def _(c):
    for rec, line in ((0, 1), (1, 0), (2, 1)):
        c.recs[rec][1][line] = ins(c.recs[rec][1][line], 5, b"-")
    return c.render(eol=b"\r\n")
@V("eol_cr_title_warning", kinds=("nuc",))
def _(c): c.recs[1][0] += b" ACGTACGTACGTACGTACGTACGT"; return c.render(eol=b"\r")
@V("eol_cr_in_defline")
def _(c): c.recs[1][0] = c.recs[1][0].replace(b" ", b"\rcr ", 1); return c.render()
@V("eol_cr_in_data")
def _(c): return edit(c, 1, 1, lambda b: ins(b, 10, b"\r"))
@V("eol_lone_cr_end")
def _(c): return c.render() + b"\r"
@V("eol_lfcr_defline", note="LF file with a LF CR pair before the second defline")
def _(c): return c.render().replace(b"\n>" + c.recs[1][0], b"\n\r>" + c.recs[1][0], 1)
@V("eol_crlf_then_lf_in_data")
def _(c):
    out = c.render(eol=b"\r\n")
    return out.replace(c.recs[1][1][1] + b"\r\n", c.recs[1][1][1] + b"\n", 1)


# byte-exact files of the line-ends inventory (E2h LR): the line reader pushes back the rest of a
# line without its terminator when a second end-of-line style appears, so a break can be lost
LOST = {
    "eol_lost_a": b">f\rAC-GTACGTAGGCTAGCTAGGATCGATCG\rGG-ACACGTAGGCTAGCTAGGATCGATCG\nATTAC-CAGACGTAGGCTAGCTAGGATCGATCG",
    "eol_lost_b": b">f\rAC-GTACGTAGGCTAGCTAGGATCGATCG\rGG-ACACGTAGGCTAGCTAGGATCGATCG\nATTAC-CAGACGTAGGCTAGCTAGGATCGATCG\n",
    "eol_lost_c": b">f\rAC-GTACGTAGGCTAGCTAGGATCGATCG\rGG-ACACGTAGGCTAGCTAGGATCGATCG\nATTAC-CAGACGTAGGCTAGCTAGGATCGATCG\rTT-GGACGTAGGCTAGCTAGGATCGATCG\r",
    "eol_lost_d": b">f\nAC-GTACGTAGGCTAGCTAGGATCGATCG\nGG-ACACGTAGGCTAGCTAGGATCGATCG\rATTAC-CAGACGTAGGCTAGCTAGGATCGATCG",
    "eol_lost_g55": b">first\r\nTC-TTGGCTCAATCCTAGGTGGGCATGTTTCCTAATGCCC\rTT-TTTAACGTGAGGGTTCGCGTTTTTATCCCACCTAGC\r\r\nACGTACGT\n-ACGTACGTACGTACGTACGTAC",
    "eol_lost_g55_nl": b">first\r\nTC-TTGGCTCAATCCTAGGTGGGCATGTTTCCTAATGCCC\rTT-TTTAACGTGAGGGTTCGCGTTTTTATCCCACCTAGC\r\r\nACGTACGT\n-ACGTACGTACGTACGTACGTAC\n",
}
for _name, _data in LOST.items():
    VARIANTS[_name] = dict(fn=(lambda c, d=_data: d), tier="core" if _name in ("eol_lost_a", "eol_lost_g55") else "full",
                           kinds=("nuc",), roles=("q", "s"), lcase=False, threads=False,
                           note="byte-exact from the LR inventory (pushback quirk); nucleotide text")


# ---------------------------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------------------------
BLASTN_TASKS = ("megablast", "blastn", "dc-megablast", "blastn-short")
SEQID_LINES = (b"AB123456", b"gb|AB123456|", b"P01308", b"SRA:SRR000001")
# the clean partner file of a program and role
PARTNER = {("blastn", "q"): "base.s.nuc.fa", ("blastn", "s"): "base.q.nuc.fa",
           ("tblastx", "q"): "base.s.nuc.fa", ("tblastx", "s"): "base.q.nuc.fa",
           ("tblastn", "q"): "base.s.nuc.fa", ("tblastn", "s"): "base.q.prot.faa",
           ("blastp", "q"): "base.s.prot.faa", ("blastp", "s"): "base.q.prot.faa"}


def kind_of(program: str, role: str) -> str:
    """Molecule of the file in a role."""
    if program in ("blastn", "tblastx"):
        return "nuc"
    if program == "blastp":
        return "prot"
    return "prot" if role == "q" else "nuc"        # tblastn


def program_set(kind: str, role: str, tier: str) -> list[tuple[str, str, tuple]]:
    """(program, task, outfmts) for a variant of a molecule and role. Core variants run all four
    BLASTN tasks and outfmt 0, 6 and 7; the others run megablast and blastn and fewer outfmts."""
    out: list[tuple[str, str, tuple]] = []
    if tier == "core":
        if kind == "nuc":
            out += [("blastn", "megablast", (6, 0, 7)), ("blastn", "blastn", (6, 0)),
                    ("blastn", "dc-megablast", (0,)), ("blastn", "blastn-short", (0,)), ("tblastx", "", (6, 0))]
            if role == "s":
                out.append(("tblastn", "", (0,)))
        else:
            out.append(("blastp", "", (6, 0, 7)))
            if role == "q":
                out.append(("tblastn", "", (0,)))
        return out
    if kind == "nuc":
        out += [("blastn", "megablast", (6, 0)), ("blastn", "blastn", (0,)), ("tblastx", "", (0,) if role == "q" else (6,))]
        if role == "s":
            out.append(("tblastn", "", (0,)))
    else:
        out.append(("blastp", "", (6, 0) if role == "q" else (0,)))
        if role == "q":
            out.append(("tblastn", "", (0,)))
    return out


def make_row(row_id, program, task, query, subject, args="", outfmt=6, threads=1, stdin="", env="",
             kind="ncbi", note="", tier="extra") -> Row:
    return Row(row_id, program, task, query, subject, args, str(outfmt), str(threads), stdin, env, kind, tier, note=note)


def variant_file(name: str, role: str, kind: str) -> str:
    return f"@/{name}.{role}.{kind}.{EXT[kind]}"


def build_rows() -> list[Row]:
    rows: list[Row] = []
    for name, v in VARIANTS.items():
        argsets = [("", "")] + ([("-lcase_masking", ".lcm")] if v["lcase"] else [])
        for role in v["roles"]:
            for kind in v["kinds"]:
                for program, task, fmts in program_set(kind, role, v["tier"]):
                    path = variant_file(name, role, kind)
                    partner = f"@/{PARTNER[(program, role)]}"
                    query, subject = (path, partner) if role == "q" else (partner, path)
                    for args, tag in argsets:
                        if args and program in ("tblastx", "blastp") and v["tier"] != "core":
                            continue
                        for fmt in fmts:
                            rid = f"{program}{'.' + task if task else ''}.{name}{tag}.{role}.o{fmt}"
                            note = v["note"]
                            if args and program in ("tblastx", "blastp"):
                                note = (note + "; " if note else "") + "-lcase_masking: " + OPTION_NOTE
                            rows.append(make_row(rid, program, task, query, subject, args, fmt, note=note, tier=v["tier"]))
                            if v["threads"] and not args and (fmt == 0 or (fmt == 6 and task == "megablast")) \
                                    and task in ("", "megablast"):
                                for n in (2, 4):
                                    rows.append(make_row(f"{rid}.t{n}", program, task, query, subject, args, fmt, n,
                                                         note=v["note"], tier=v["tier"]))
    rows += extra_rows()
    seen = set()
    for r in rows:
        assert r.row_id not in seen, r.row_id
        seen.add(r.row_id)
    return rows


ALL_PROGRAMS = (("blastn", "megablast"), ("blastn", "blastn"), ("tblastx", ""), ("tblastn", ""), ("blastp", ""))


def base_pair(program: str) -> tuple[str, str]:
    """Clean query and subject of a program."""
    q = f"@/base.q.{'prot.faa' if program in ('blastp', 'tblastn') else 'nuc.fa'}"
    s = f"@/base.s.{'prot.faa' if program == 'blastp' else 'nuc.fa'}"
    return q, s


def extra_rows() -> list[Row]:
    rows: list[Row] = []

    def add(rid, program, task, query, subject, fmts=(6, 0, 7), **kw):
        for fmt in fmts:
            rows.append(make_row(f"{program}{'.' + task if task else ''}.{rid}.o{fmt}", program, task, query, subject,
                                 outfmt=fmt, **kw))

    # standard input
    for program, task in ALL_PROGRAMS:
        q, s = base_pair(program)
        kind = "prot" if program in ("blastp", "tblastn") else "nuc"
        add("stdin_q_pipe", program, task, "-", s, stdin=f"pipe:{q}")
        add("stdin_q_file", program, task, "-", s, stdin=f"file:{q}")
        add("stdin_q_omitted", program, task, "", s, stdin=f"pipe:{q}", fmts=(6, 0))
        add("stdin_s_pipe", program, task, q, "-", stdin=f"pipe:{s}")
        add("stdin_s_file", program, task, q, "-", stdin=f"file:{s}")
        add("stdin_qs_pipe", program, task, "-", "-", stdin=f"pipe:{s}", fmts=(6, 0, 7),
            note="the subject reads standard input to its end; the query then sees EOF (no query, no 'Query is Empty!')")
        add("stdin_qs_file", program, task, "-", "-", stdin=f"file:{s}", fmts=(6, 0))
        add("stdin_q_empty_pipe", program, task, "-", s, stdin="empty", note="E2g A-042: an empty pipe is not an empty query")
        add("stdin_q_devnull", program, task, "-", s, stdin="file:/dev/null", fmts=(6, 0))
        add("stdin_q_ws_pipe", program, task, "-", s, stdin=f"pipe:@/file_ws_only.q.{kind}.{EXT[kind]}")
        add("stdin_q_comment_pipe", program, task, "-", s,
            stdin=f"pipe:@/file_comment_only.q.{kind}.{EXT[kind]}", fmts=(6, 0))
        add("stdin_s_empty_pipe", program, task, q, "-", stdin="empty", fmts=(6, 0))
        add("stdin_s_ws_pipe", program, task, q, "-",
            stdin=f"pipe:@/file_ws_only.s.{'prot' if program == 'blastp' else 'nuc'}."
                  f"{EXT['prot' if program == 'blastp' else 'nuc']}", fmts=(6, 0))
        if kind == "nuc" and program != "tblastn":
            for mode in ("pipe", "file"):
                add(f"stdin_q_crlf_{mode}", program, task, "-", s,
                    stdin=f"{mode}:{variant_file('eol_crlf', 'q', 'nuc')}", fmts=(6, 0))
                add(f"stdin_q_lost_g55_{mode}", program, task, "-", s,
                    stdin=f"{mode}:{variant_file('eol_lost_g55', 'q', 'nuc')}", fmts=(6, 0),
                    note="the pushback quirk depends on how the bytes arrive (file window, pipe)")
                add(f"stdin_s_def_utf8_{mode}", program, task, q, "-",
                    stdin=f"{mode}:{variant_file('def_utf8', 's', 'nuc')}", fmts=(6, 0))
    # lost_g55 as plain files, to compare with standard input
    # files that are not files
    for program, task in ALL_PROGRAMS:
        q, s = base_pair(program)
        add("file_missing_q", program, task, f"@/does_not_exist.fa", s, fmts=(6,))
        add("file_missing_s", program, task, q, f"@/does_not_exist.fa", fmts=(6,))
        add("dir_q", program, task, "@", s, fmts=(6, 0))
        add("dir_s", program, task, q, "@", fmts=(6, 0))
        add("both_empty_files", program, task, f"@/file_empty.q.{'prot' if program in ('blastp', 'tblastn') else 'nuc'}."
            f"{'faa' if program in ('blastp', 'tblastn') else 'fa'}",
            f"@/file_empty.s.{'prot' if program == 'blastp' else 'nuc'}.{'faa' if program == 'blastp' else 'fa'}", fmts=(6,))
    # query batches (BATCH_SIZE: residues read before the next search; NCBI reads the query record by record)
    batch_sizes = {"nuc": (("each", "100"), ("two", "200")), "prot": (("each", "50"), ("two", "100"))}
    for program, task in ALL_PROGRAMS:
        kind = "prot" if program in ("blastp", "tblastn") else "nuc"
        _, s = base_pair(program)
        for name in ("batch_clean", "batch_late_bad", "batch_late_warn", "batch_late_empty", "batch_late_allempty",
                     "batch_late_title_warn", "batch_late_hyphen_cr"):
            path = f"@/{name}.{kind}.{EXT[kind]}"
            for size_name, size in batch_sizes[kind]:
                add(f"{name}.{size_name}", program, task, path, s, env=f"BATCH_SIZE={size}",
                    fmts=(6, 0, 7) if name in ("batch_late_bad", "batch_late_warn", "batch_late_empty",
                                               "batch_late_allempty") else (6, 0))
            if name != "batch_clean":
                add(f"{name}.default", program, task, path, s, fmts=(6, 0))
    # default batches of BLASTN: the first batch holds five queries of 1000 nucleotides or more
    for name in ("batchlong_late_bad", "batchlong_late_warn", "batchlong_late_empty", "batchlong_late_allempty", "batchlong_clean"):
        for task in ("megablast", "blastn"):
            add(name, "blastn", task, f"@/{name}.nuc.fa", f"@/batchlong.s.nuc.fa", fmts=(6, 0, 7))
    # a subject cannot be in a later batch, but its errors and warnings come before the query's
    for program, task in ALL_PROGRAMS:
        q, s = base_pair(program)
        kind_s = "prot" if program == "blastp" else "nuc"
        kind_q = "prot" if program in ("blastp", "tblastn") else "nuc"
        add("both_bad_subject_first", program, task, variant_file("check_text_first", "q", kind_q),
            variant_file("check_text_first", "s", kind_s), fmts=(6, 0))
        add("both_warn_order", program, task, variant_file("seq_hyphen_lines", "q", kind_q),
            variant_file("seq_hyphen_lines", "s", kind_s), fmts=(6, 0, 7))
        add("both_empty_records", program, task, variant_file("rec_empty_mid", "q", kind_q),
            variant_file("rec_empty_mid", "s", kind_s), fmts=(6, 0, 7))
        add("both_all_empty_records", program, task, variant_file("rec_empty_all", "q", kind_q),
            variant_file("rec_empty_all", "s", kind_s), fmts=(6, 0))
        add("both_empty_query_file_ws_subject", program, task, variant_file("file_ws_only", "q", kind_q),
            variant_file("rec_empty_all", "s", kind_s), fmts=(6, 0))
    # protein text in nucleotide programs and nucleotide text in BLASTP (AssignMolType never guesses)
    for task in ("megablast", "blastn"):
        add("prot_in_blastn", "blastn", task, f"@/base.q.prot.faa", f"@/base.s.nuc.fa", fmts=(6, 0))
        add("prot_in_blastn_subject", "blastn", task, f"@/base.q.nuc.fa", f"@/base.s.prot.faa", fmts=(6, 0))
    add("nuc_in_blastp", "blastp", "", f"@/base.q.nuc.fa", f"@/base.s.prot.faa", fmts=(6, 0))
    add("nuc_in_blastp_subject", "blastp", "", f"@/base.q.prot.faa", f"@/base.s.nuc.fa", fmts=(6, 0))
    add("prot_in_tblastn_subject", "tblastn", "", f"@/base.q.prot.faa", f"@/base.s.prot.faa", fmts=(6, 0))
    add("nuc_in_tblastn_query", "tblastn", "", f"@/base.q.nuc.fa", f"@/base.s.nuc.fa", fmts=(6, 0))
    add("prot_in_tblastx", "tblastx", "", f"@/base.q.prot.faa", f"@/base.s.nuc.fa", fmts=(6, 0))
    # Seq-id first lines: NCBI fetches these over the network, so they are never run on NCBI. LOSAT must
    # reject them with "not supported by LOSAT"; the port fills the expected output.
    for i, line in enumerate(SEQID_LINES):
        for program, task, role in (("blastn", "megablast", "q"), ("blastn", "blastn", "q"), ("blastn", "megablast", "s"),
                                    ("tblastx", "", "q"), ("tblastn", "", "q"), ("tblastn", "", "s"),
                                    ("blastp", "", "q"), ("blastp", "", "s")):
            if program in ("tblastx", "tblastn", "blastp") and i > 1 and not (program == "blastp" and i == 2):
                continue
            kind = kind_of(program, role)
            path = f"@/seqid{i}.{role}.{kind}.{EXT[kind]}"
            partner = f"@/{PARTNER[(program, role)]}"
            q, s = (path, partner) if role == "q" else (partner, path)
            rows.append(make_row(f"{program}{'.' + task if task else ''}.seqid{i}.{role}.o6", program, task, q, s, outfmt=6,
                                 kind="losat-rejection", note=f"first line {line.decode()}: NCBI fetches it (network)",
                                 tier="seqid"))
    return rows


# ---------------------------------------------------------------------------------------------
# generate
# ---------------------------------------------------------------------------------------------
def batch_files() -> dict[str, bytes]:
    out = {}
    for kind in ("nuc", "prot"):
        def bctx():
            return Ctx(kind, "bq")

        c = bctx()
        out[f"batch_clean.{kind}.{EXT[kind]}"] = c.render()
        c = bctx(); c.recs[5][1].insert(0, BAD_TEXT[kind])                       # record 6: CheckDataLine error
        out[f"batch_late_bad.{kind}.{EXT[kind]}"] = c.render()
        c = bctx(); c.recs[5][1][1] = ins(c.recs[5][1][1], 8, b"7")             # record 6: invalid residue warning
        c.recs[3][1][0] = ins(c.recs[3][1][0], 5, b"-")
        out[f"batch_late_warn.{kind}.{EXT[kind]}"] = c.render()
        c = bctx(); c.recs[5][1] = []                                            # record 6 empty
        out[f"batch_late_empty.{kind}.{EXT[kind]}"] = c.render()
        c = bctx(); c.recs[4][1] = []; c.recs[5][1] = []; c.recs[6][1] = []     # records 5, 6 and 7 empty: the last batch is empty
        out[f"batch_late_allempty.{kind}.{EXT[kind]}"] = c.render()
        c = bctx()
        if kind == "nuc":
            c.recs[4][0] += b" ACGTACGTACGTACGTACGTAC"
        else:
            c.recs[4][0] += b" " + (b"ACDEFGHIKLMNPQRSTVWY" * 3)[:55]
        out[f"batch_late_title_warn.{kind}.{EXT[kind]}"] = c.render()
        c = bctx(); c.recs[4][1][0] = ins(c.recs[4][1][0], 5, b"-")
        out[f"batch_late_hyphen_cr.{kind}.{EXT[kind]}"] = c.render(eol=b"\r")
    return out


def batchlong_files() -> dict[str, bytes]:
    """Eight queries of about 1100 nucleotides: blastn's first batch is five such queries."""
    rng = random.Random(20261009)
    subject_parts, queries = [], []
    for i in range(8):
        q = "".join(rng.choice("ACGT") for _ in range(1100))
        queries.append(q.encode())
        subject_parts.append(("".join(rng.choice("ACGT") for _ in range(100)) + q[100:400]).encode())
    out = {}

    def render(recs):
        return b"".join(b">" + d + b"\n" + b"\n".join(wrap_lines(s, 70)) + b"\n" if s is not None else b">" + d + b"\n"
                        for d, s in recs)

    recs = [(f"lq{i + 1} long query {i + 1}".encode(), q) for i, q in enumerate(queries)]
    out["batchlong_clean.nuc.fa"] = render(recs)
    out["batchlong_late_bad.nuc.fa"] = render(recs[:6]) + b">" + recs[6][0] + b"\n" + BAD_TEXT["nuc"] + b"\n" + \
        b"\n".join(wrap_lines(recs[6][1], 70)) + b"\n" + render(recs[7:])
    warn = list(recs); warn[6] = (warn[6][0], warn[6][1][:50] + b"7" + warn[6][1][50:]); warn[5] = (warn[5][0] + b" ACGTACGTACGTACGTACGTAC", warn[5][1])
    out["batchlong_late_warn.nuc.fa"] = render(warn)
    empty = list(recs); empty[6] = (empty[6][0], b""); empty[5] = (empty[5][0], b"")
    out["batchlong_late_empty.nuc.fa"] = render([(d, s if s else None) for d, s in empty])
    allempty = [(d, s if i < 5 else None) for i, (d, s) in enumerate(recs)]
    out["batchlong_late_allempty.nuc.fa"] = render(allempty)
    out["batchlong.s.nuc.fa"] = b"".join(b">ls%d long subject %d\n" % (i + 1, i + 1) + b"\n".join(wrap_lines(s, 70)) + b"\n"
                                         for i, s in enumerate(subject_parts))
    return out


def seqid_files() -> dict[str, bytes]:
    out = {}
    for i, line in enumerate(SEQID_LINES):
        for role in ("q", "s"):
            for kind in ("nuc", "prot"):
                c = Ctx(kind, role)
                out[f"seqid{i}.{role}.{kind}.{EXT[kind]}"] = line + b"\n" + c.render()
    return out


def command_generate(args) -> int:
    FIXTURES.mkdir(parents=True, exist_ok=True)
    files: dict[str, bytes] = {}
    for kind in ("nuc", "prot"):
        for role in ("q", "s"):
            files[f"base.{role}.{kind}.{EXT[kind]}"] = Ctx(kind, role).render()
    for name, v in VARIANTS.items():
        for role in v["roles"]:
            for kind in v["kinds"]:
                files[f"{name}.{role}.{kind}.{EXT[kind]}"] = v["fn"](Ctx(kind, role))
    files.update(batch_files())
    files.update(batchlong_files())
    files.update(seqid_files())
    for existing in FIXTURES.iterdir():
        if existing.is_file() and existing.name not in files:
            existing.unlink()
    total = 0
    for name, data in files.items():
        path = FIXTURES / name
        if not path.exists() or path.read_bytes() != data:
            path.write_bytes(data)
        total += len(data)
    print(f"{len(files)} input files, {total} bytes in {FIXTURES}")
    return 0


# ---------------------------------------------------------------------------------------------
# Running
# ---------------------------------------------------------------------------------------------
def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def clean_env(case_env: str = "") -> dict[str, str]:
    env = {k: v for k, v in os.environ.items()
           if k not in REPORT_ENV and not k.startswith(("LOSAT_", "RAYON_", "NCBI_CONFIG__"))}
    for item in shlex.split(case_env):
        key, _, value = item.partition("=")
        env[key] = value
    return env


def rp(value: str) -> str:
    """A path of a row: `@/name` is tests/fixtures/fasta_input/name, `@` the directory."""
    if value == "@":
        return REL
    return REL + value[1:] if value.startswith("@/") else value


def stdin_of(row: Row) -> tuple[str, str]:
    """(mode, path) of the standard input of a row; mode is "", "pipe", "file" or "empty"."""
    mode, _, path = row.stdin.partition(":")
    return mode, rp(path)


def row_argv(row: Row, losat: bool) -> list[str]:
    argv: list[str] = []
    if row.query != "":
        argv += ["-query", rp(row.query)]
    if row.subject != "":
        argv += ["-subject", rp(row.subject)]
    if row.task:
        argv += ["-task", row.task]
    argv += shlex.split(row.args)
    argv += ["-outfmt", row.outfmt]
    if losat and int(row.threads) > 1:
        argv += ["-num_threads", row.threads]
    return argv


def execute(command: list[str], row: Row, merged: bool) -> tuple[int, bytes, bytes]:
    kw: dict = dict(cwd=ENGINE, env=clean_env(row.env), timeout=TIMEOUT)
    handle = None
    try:
        mode, path = stdin_of(row)
        if mode == "pipe":
            kw["input"] = (ENGINE / path).read_bytes()
        elif mode == "empty":
            kw["input"] = b""
        elif mode == "file":
            handle = open(ENGINE / path, "rb")
            kw["stdin"] = handle
        else:
            kw["stdin"] = subprocess.DEVNULL
        if merged:
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, **kw)
            return result.returncode, result.stdout, b""
        result = subprocess.run(command, capture_output=True, **kw)
        return result.returncode, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return -999, b"", b"timeout"
    finally:
        if handle:
            handle.close()


def streams(command: list[str], row: Row) -> dict[str, object]:
    rc, out, err = execute(command, row, False)
    rc2, merged, _ = execute(command, row, True)
    return {"rc": rc, "rc_merged": rc2, "stdout": out, "stderr": err, "combined": merged}


def seqid_risk(data: bytes) -> bool:
    """True when NCBI would try the first line as a Seq-id (CBlastInputReader::ReadOneSeq): the line,
    trimmed of white space, starts with a letter or digit and is not letters only."""
    line = re.match(rb"[^\r\n]*", data).group().strip(b" \t\n\v\f\r")
    if not line or not line[:1].isalnum():
        return False
    return re.fullmatch(rb"[A-Za-z]+", line) is None


def network_safe(row: Row) -> str | None:
    """None when no input of the row starts with a Seq-id candidate, else a reason."""
    sources = []
    mode, spath = stdin_of(row)
    stdin_file = ENGINE / spath if mode in ("pipe", "file") and spath != "/dev/null" else None
    for role, value in (("query", row.query), ("subject", row.subject)):
        if value == "-" and stdin_file:
            sources.append((f"{role} (standard input)", stdin_file.read_bytes()))
        elif value and value != "-" and (ENGINE / rp(value)).is_file():
            sources.append((role, (ENGINE / rp(value)).read_bytes()))
    if row.query == "" and stdin_file:
        sources.append(("query (standard input)", stdin_file.read_bytes()))
    for role, data in sources:
        if seqid_risk(data):
            return f"the first line of the {role} could be a Seq-id"
    return None


def filter_rows(rows: list[Row], selector: str | None, match: str | None = None, tier: str | None = None) -> list[Row]:
    if tier:
        rows = [r for r in rows if r.tier in tier.split(",")]
    if selector:
        wanted = set(selector.split(","))
        rows = [r for r in rows if r.program in wanted]
    if match:
        rows = [r for r in rows if match in r.row_id]
    return rows


MANIFEST_HEAD = (
    "# Frozen outputs of {version} for FASTA-input cases (comparison oracle only; see fasta_input_fixtures.py).\n"
    "# Only what `freeze` measured is stored; the row definitions (inputs, args, outfmt, threads, stdin, env, kind, tier, note) are in the script's generator.\n"
    "# `@S <n> <sha256>` lines are a table of the distinct non-empty stderr hashes. Row lines: row_id, rc, stdout sha256,\n"
    "# stderr (`-` = empty stream, else a table index n), combined `2>&1` stream (`=` = same as stdout, `~` = same as stderr, else its sha256).\n"
    "# Rows of kind losat-rejection and rows not frozen yet have no line (their expectation is `pending`). Written by `freeze`/`refresh`; do not edit.\n")


def manifest_files() -> list[Path]:
    return sorted(MANIFEST_DIR.glob("*.tsv")) if MANIFEST_DIR.is_dir() else []


def manifest_version() -> str:
    for path in manifest_files():
        first = path.read_text().splitlines()[0]
        return first.removeprefix("# Frozen outputs of ").split(" for FASTA-input")[0]
    return "NCBI BLAST+ 2.17.0"


def load_frozen() -> dict[str, tuple[str, str, str, str]]:
    """row_id -> (rc, stdout_sha256, stderr_sha256, combined_sha256) from the per-program manifests."""
    frozen: dict[str, tuple[str, str, str, str]] = {}
    for path in manifest_files():
        table: dict[str, str] = {}
        for line in path.read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if cols[0] == "@S":
                table[cols[1]] = cols[2]
                continue
            row_id, rc, out, err, comb = cols
            err = EMPTY_SHA256 if err == "-" else table[err]
            comb = out if comb == "=" else err if comb == "~" else comb
            assert row_id not in frozen, row_id
            frozen[row_id] = (rc, out, err, comb)
    return frozen


def read_manifest() -> list[Row]:
    """The generator's rows with the frozen expectation filled in (`pending` for rows without one)."""
    frozen = load_frozen()
    rows = build_rows()
    for row in rows:
        if row.kind == "ncbi" and row.row_id in frozen:
            row.rc, row.stdout_sha256, row.stderr_sha256, row.combined_sha256 = frozen[row.row_id]
    return rows


def write_manifest(rows: list[Row], version: str) -> None:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    by_program: dict[str, list[Row]] = {}
    for row in rows:
        if row.kind == "ncbi" and row.rc != "pending":
            by_program.setdefault(row.program, []).append(row)
    for stale in manifest_files():
        if stale.stem not in by_program:
            stale.unlink()
    for program, group in by_program.items():
        table: dict[str, int] = {}
        for row in group:
            if row.stderr_sha256 != EMPTY_SHA256:
                table.setdefault(row.stderr_sha256, len(table) + 1)
        lines = [MANIFEST_HEAD.format(version=version)]
        lines += [f"@S\t{n}\t{h}\n" for h, n in table.items()]
        for row in group:
            err = "-" if row.stderr_sha256 == EMPTY_SHA256 else str(table[row.stderr_sha256])
            comb = "=" if row.combined_sha256 == row.stdout_sha256 else "~" if row.combined_sha256 == row.stderr_sha256 \
                else row.combined_sha256
            lines.append(f"{row.row_id}\t{row.rc}\t{row.stdout_sha256}\t{err}\t{comb}\n")
        (MANIFEST_DIR / f"{program}.tsv").write_text("".join(lines))


def ncbi_command(ncbi_bin: Path, row: Row, netns: bool) -> list[str]:
    command = [str(ncbi_bin / row.program), *row_argv(row, False)]
    return ["unshare", "-rn", *command] if netns else command


def freeze_one(ncbi_bin: Path, out: Path, row: Row, netns: bool) -> tuple[Row, str]:
    reason = network_safe(row)
    if reason:
        return row, f"refused: {reason}"
    result = streams(ncbi_command(ncbi_bin, row, netns), row)
    if result["rc"] == -999 or result["rc_merged"] == -999:
        return row, "error: timeout"
    (out / f"{row.row_id}.out").write_bytes(result["stdout"])
    (out / f"{row.row_id}.err").write_bytes(result["stderr"])
    (out / f"{row.row_id}.merged").write_bytes(result["combined"])
    (out / f"{row.row_id}.rc").write_text(f"{result['rc']}\n")
    note = "" if result["rc"] == result["rc_merged"] else f"rc differs between the two runs ({result['rc']} / {result['rc_merged']})"
    row.rc = str(result["rc"])
    row.stdout_sha256 = sha256(result["stdout"])
    row.stderr_sha256 = sha256(result["stderr"])
    row.combined_sha256 = sha256(result["combined"])
    return row, note


def command_freeze(args) -> int:
    if any(key in os.environ for key in REPORT_ENV) or (Path.home() / ".ncbirc").exists():
        raise SystemExit(f"unset {REPORT_ENV} and remove ~/.ncbirc before freezing")
    ncbi_bin = Path(args.ncbi_bin).resolve()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    version = subprocess.run([str(ncbi_bin / "blastn"), "-version"], capture_output=True, text=True).stdout.strip().replace("\n", "; ")
    netns = not args.no_netns and subprocess.run(["unshare", "-rn", "true"], capture_output=True).returncode == 0
    defined = build_rows()
    selected = filter_rows([r for r in defined if r.kind == "ncbi"], args.rows, args.match, args.tier)
    previous = {r.row_id: r for r in read_manifest()}
    if args.verify:
        selected = [r for r in selected if previous[r.row_id].rc != "pending"]
    print(f"{len(selected)} rows, netns={netns}, {version}", flush=True)
    problems = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(lambda r: freeze_one(ncbi_bin, out, r, netns), selected))
    done: dict[str, Row] = {}
    for row, note in results:
        if note.startswith(("refused", "error")):
            problems += 1
            print(f"{row.row_id}\t{note}", flush=True)
            continue
        if args.verify:
            old = previous[row.row_id]
            if EXPECTED_OF(old) != EXPECTED_OF(row):
                problems += 1
                print(f"{row.row_id}\tunstable: differs from the manifest", flush=True)
        elif note:
            print(f"{row.row_id}\t{note}", flush=True)
        done[row.row_id] = row
    if args.verify:
        print(f"verified {len(done)} rows, {problems} problems")
        return 1 if problems else 0
    merged = []
    for row in defined:
        if row.row_id in done:
            merged.append(done[row.row_id])
        else:
            old = previous.get(row.row_id)
            if old is not None and old.kind == row.kind:
                for f in EXPECTED:
                    setattr(row, f, getattr(old, f))
            merged.append(row)
    write_manifest(merged, version)
    print(f"{len(done)} frozen, {problems} problems, {sum(1 for r in merged if r.kind == 'losat-rejection')} losat-rejection rows")
    return 1 if problems else 0


def command_refresh(args) -> int:
    """Rewrite the manifests from the existing hashes (no NCBI run): normalises the files, drops stale rows."""
    rows = read_manifest()
    frozen = load_frozen()
    ids = {r.row_id for r in rows}
    stale = sorted(set(frozen) - ids)
    new = sum(1 for r in rows if r.kind == "ncbi" and r.rc == "pending")
    write_manifest(rows, manifest_version())
    print(f"{len(rows)} rows, {new} without a frozen expectation (run freeze for them), {len(stale)} stale frozen rows dropped")
    return 0


def EXPECTED_OF(row: Row) -> tuple:
    return tuple(getattr(row, f) for f in EXPECTED)


def first_difference(expected: bytes, actual: bytes) -> str:
    left, right = expected.split(b"\n"), actual.split(b"\n")
    for number, (a, b) in enumerate(zip(left, right), 1):
        if a != b:
            return f"line {number}: expected {a[:60]!r}, got {b[:60]!r}"
    return f"line {min(len(left), len(right))}: expected {len(left) - 1} lines, got {len(right) - 1}"


def check_one(losat: Path, expected_dir: Path | None, row: Row) -> tuple[Row, str, str]:
    result = streams([str(losat), row.program, *row_argv(row, True)], row)
    if result["rc"] == -999 or result["rc_merged"] == -999:
        return row, "error", "timeout"
    first_err = result["stderr"].split(b"\n")[0][:140].decode(errors="replace")
    if row.stdout_sha256 == "pending":
        return row, "pending", re.sub(r"[\t\r\n]+", " ", f"rc {result['rc']}: {first_err}")
    actual = {"rc": str(result["rc"]), "stdout_sha256": sha256(result["stdout"]),
              "stderr_sha256": sha256(result["stderr"]), "combined_sha256": sha256(result["combined"])}
    names = {"rc": "exit", "stdout_sha256": "stdout", "stderr_sha256": "stderr", "combined_sha256": "combined"}
    differing = [names[f] for f in EXPECTED if actual[f] != getattr(row, f)]
    if not differing:
        return row, "same", ""
    detail = []
    if "exit" in differing:
        detail.append(f"exit {actual['rc']}, expected {row.rc}")
    files = {"stdout": ("out", result["stdout"]), "stderr": ("err", result["stderr"]), "combined": ("merged", result["combined"])}
    for name in ("stdout", "stderr", "combined"):
        if name in differing and expected_dir is not None and (expected_dir / f"{row.row_id}.{files[name][0]}").exists():
            expected = (expected_dir / f"{row.row_id}.{files[name][0]}").read_bytes()
            detail.append(f"{name} {first_difference(expected, files[name][1])}")
            break
    else:
        detail.append("differs: " + ",".join(d for d in differing if d != "exit"))
    detail.append(f"losat stderr: {first_err}")
    detail = [re.sub(r"[\t\r\n]+", " ", d) for d in detail]
    rejected = b"not supported by LOSAT" in result["stderr"] and b"not supported by LOSAT" not in (
        (expected_dir / f"{row.row_id}.err").read_bytes() if expected_dir and (expected_dir / f"{row.row_id}.err").exists() else b"")
    return row, "rejects" if rejected else "differs", "; ".join(detail)


def command_check(args) -> int:
    losat = Path(args.losat).resolve()
    expected_dir = Path(args.expected) if args.expected else None
    if expected_dir is None:
        build = os.environ.get("BUILD_ROOT")
        default = Path(build) / "sf-e2h/fixtures/ncbi" if build else None
        expected_dir = default if default and default.is_dir() else None
    rows = filter_rows(read_manifest(), args.rows, args.match, args.tier)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(lambda r: check_one(losat, expected_dir, r), rows))
    def group(row: Row) -> str:
        return "option-rejected" if OPTION_NOTE in row.note else ""

    lines = ["row_id\tprogram\tkind\ttier\tresult\tdetail\tgroup"] + [
        f"{row.row_id}\t{row.program}\t{row.kind}\t{row.tier}\t{verdict}\t{detail}\t{group(row)}" for row, verdict, detail in results]
    if args.out:
        Path(args.out).write_text("\n".join(lines) + "\n")
    counts: dict[str, int] = {}
    for _, verdict, _ in results:
        counts[verdict] = counts.get(verdict, 0) + 1
    print(" ".join(f"{k}={v}" for k, v in sorted(counts.items())), f"rows={len(results)}")
    grouped: dict[str, int] = {}
    for row, verdict, _ in results:
        if group(row):
            grouped[verdict] = grouped.get(verdict, 0) + 1
    if grouped:
        print(f"separate group, {OPTION_NOTE}: " + " ".join(f"{k}={v}" for k, v in sorted(grouped.items())),
              f"rows={sum(grouped.values())}")
    return 0 if set(counts) <= {"same", "pending"} else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("generate")
    sub.add_parser("refresh", help="rewrite the manifest from the row definitions, keeping the frozen hashes")
    freeze = sub.add_parser("freeze")
    freeze.add_argument("--ncbi-bin", required=True)
    freeze.add_argument("--out", required=True)
    freeze.add_argument("--jobs", type=int, default=3)
    freeze.add_argument("--rows", help="comma-separated programs")
    freeze.add_argument("--match", help="only rows whose id contains this text")
    freeze.add_argument("--tier", help="comma-separated tiers: core, full, extra, seqid")
    freeze.add_argument("--verify", action="store_true", help="rerun NCBI and compare with the manifest; write nothing")
    freeze.add_argument("--no-netns", action="store_true", help="do not wrap NCBI in `unshare -rn`")
    check = sub.add_parser("check")
    check.add_argument("--losat", required=True)
    check.add_argument("--out")
    check.add_argument("--jobs", type=int, default=3)
    check.add_argument("--rows", help="comma-separated programs")
    check.add_argument("--match", help="only rows whose id contains this text")
    check.add_argument("--tier", help="comma-separated tiers: core, full, extra, seqid")
    check.add_argument("--expected", help="directory of the frozen files (default $BUILD_ROOT/sf-e2h/fixtures/ncbi)")
    args = parser.parse_args()
    return {"generate": command_generate, "refresh": command_refresh, "freeze": command_freeze,
            "check": command_check}[args.action](args)


if __name__ == "__main__":
    sys.exit(main())
