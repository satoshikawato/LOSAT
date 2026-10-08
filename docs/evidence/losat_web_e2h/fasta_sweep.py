#!/usr/bin/env python3
"""FASTA mutation sweep against NCBI BLAST+ for the CFastaReader port (E2h, session SFb; comparison only).

Shape of the E2g audit's fuzz_fasta.py and fuzz_defline.py (byte insertions and substitutions in
deflines and sequence lines, line-end kinds, record structure), made deterministic (every file
comes from a seeded generator and is written once), run for the whole matrix

    program   blastn (default task, megablast), blastn -task blastn, tblastx, tblastn, blastp
    role      query, subject  (the other role is a clean file; TBLASTN's query is protein, its subject
              nucleotide)
    outfmt    0 and 6

and compared with the frozen output of NCBI BLAST+ 2.17.0 (stdout, stderr, exit code and the
`2>&1` stream of a second run). Every case ends in one class:

    same                 NCBI succeeds; LOSAT's exit code, stdout, stderr and `2>&1` equal NCBI's
    same-error           NCBI fails; LOSAT has the same exit code and NCBI's stdout, stderr and `2>&1`
    explicit-rejection   LOSAT exits non-zero and stderr contains "not supported by LOSAT's <PROGRAM>".
                         Column `listed` says whether the message is one of the rejections that
                         docs/evidence/losat_web_e2h/AUTHORITY.md section J lists (LISTED_REJECTIONS);
                         an unlisted rejection is a failure of the port.
    approved-exception   NCBI dies of a signal and LOSAT exits 0 (approved exception 2 of
                         PD-LOSAT-NCBI-DEFECTS, a title of punctuation; the record is named in the
                         `approval` column). The content of the report is checked by E2g's title_sweep.py.
    differs              anything else, with the first difference
    timeout              a run exceeded the limit

NCBI runs inside `unshare -rn` (no network: a first line that NCBI reads as a Seq-id cannot reach a
data loader), without ~/.ncbirc and with the environment variables that change the report unset,
with the working directory = the sweep directory (so the paths in messages are relative). The few
files whose first line NCBI would try as a Seq-id are family `seqid`.

Sub-commands:
  generate     --dir DIR [--per-cell N]            write DIR/inputs/ and DIR/cases.tsv (deterministic)
  freeze-ncbi  --dir DIR --ncbi-bin BIN --jobs N [--out DIR/ncbi]    run NCBI, store rc/out/err/merged
  check        --dir DIR --losat BIN --jobs N --out TSV [--ncbi DIR/ncbi]   classify LOSAT against it
  summary      TSV [TSV ...]                       class counts per program x role x outfmt (markdown)
  compare-frozen A B                               compare two frozen directories (determinism)

Usage example (what README.md lists):
  flock "$BUILD_ROOT/oracle.lock" fasta_sweep.py freeze-ncbi --dir D --ncbi-bin "$NCBI_BIN" --jobs 3
"""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import os
import random
import re
import subprocess
import sys
import time
from pathlib import Path

TIMEOUT = 120
# NCBI reads these (and ~/.ncbirc); each one changes the batches or the report.
REPORT_ENV = ("BL2SEQ_LEGACY", "CTOOLKIT_COMPATIBLE", "OLD_FSC", "BATCH_SIZE", "CHUNK_SIZE", "ADAPTIVE_CBS",
              "OVERLAP_CHUNK_SIZE", "PRE_FETCH_SEQS_LIMIT", "BLASTDB", "BLASTINPUT_GEN_DELTA_SEQ",
              "BLASTINPUT_GEN_DELTA_SEQ", "NCBI", "BLAST_USAGE_REPORT")
SEED = 20261008
REJECT_RE = re.compile(rb"not supported by LOSAT's (BLASTN|TBLASTX|TBLASTN|BLASTP)")
# AUTHORITY.md section J: the explicit rejections LOSAT keeps (message fragments).
#   J-1 Seq-id lines, J-2 -parse_deflines, J-4 records over 2^31-1 letters, J-6 non-UTF-8 Subject_ titles
LISTED_REJECTIONS = (rb"may be a sequence identifier", rb"-parse_deflines", rb"longer than 2147483647 letters",
                     rb"non-UTF-8")
APPROVAL_ARG_ERROR = "PD-LOSAT-CLI-NONSEARCH-DIFFERENCES approved exception 1 (argument parser messages)"
APPROVAL_EXCEPTION_2 = "PD-LOSAT-NCBI-DEFECTS approved exception 2 (punctuation title; docs/evidence/losat_web_e2g/title_sweep.py)"

# program label -> (executable, -task value, (query molecule, subject molecule))
PROGRAMS = {
    "blastn": ("blastn", "", ("nuc", "nuc")),
    "blastn-task-blastn": ("blastn", "blastn", ("nuc", "nuc")),
    "tblastx": ("tblastx", "", ("nuc", "nuc")),
    "tblastn": ("tblastn", "", ("prot", "nuc")),
    "blastp": ("blastp", "", ("prot", "prot")),
}
OUTFMTS = ("0", "6")
ROLES = ("q", "s")


# ---------------------------------------------------------------------------------------------
# Base material: three protein fragments (60 aa) with a coding sequence each (180 nt)
# ---------------------------------------------------------------------------------------------
AA = "ACDEFGHIKLMNPQRSTVWY"
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS: dict[str, list[str]] = {}
for _i, _aa in enumerate(_AAS):
    CODONS.setdefault(_aa, []).append(_BASES[_i // 16] + _BASES[(_i // 4) % 4] + _BASES[_i % 4])
WRAP = {"nuc": 60, "prot": 30}


def make_base() -> dict[tuple[str, str], list[bytes]]:
    rng = random.Random(SEED)
    prot = ["M" + "".join(rng.choice(AA) for _ in range(59)) for _ in range(3)]
    coding = ["".join(rng.choice(CODONS[a]) for a in p) for p in prot]

    def nuc(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    def aas(n):
        return "".join(rng.choice(AA) for _ in range(n))

    return {("nuc", "q"): [c.encode() for c in coding], ("prot", "q"): [p.encode() for p in prot],
            ("nuc", "s"): [(nuc(60) + c + nuc(60)).encode() for c in coding],
            ("prot", "s"): [(aas(30) + p + aas(30)).encode() for p in prot]}


BASE = make_base()


def wrap_lines(seq: bytes, width: int) -> list[bytes]:
    return [seq[i:i + width] for i in range(0, len(seq), width)]


class Doc:
    """The records of a clean file that a mutation edits. recs: [defline without '>' (None: no
    defline), [sequence lines]]."""

    def __init__(self, kind: str, role: str):
        self.kind, self.role = kind, role
        prefix = "q" if role == "q" else "s"
        words = ("first", "second", "third")
        self.recs = [[f"{prefix}{i + 1} {words[i]} {'query' if role == 'q' else 'subject'}".encode(),
                      wrap_lines(seq, WRAP[kind])] for i, seq in enumerate(BASE[(kind, role)])]
        self.eol, self.final = b"\n", True
        self.raw: bytes | None = None          # set by a mutation that builds the bytes itself

    def copy(self) -> "Doc":
        other = Doc.__new__(Doc)
        other.kind, other.role = self.kind, self.role
        other.recs = [[d, list(lines)] for d, lines in self.recs]
        other.eol, other.final, other.raw = self.eol, self.final, self.raw
        return other

    def lines(self) -> list[bytes]:
        out: list[bytes] = []
        for defline, lines in self.recs:
            if defline is not None:
                out.append(b">" + defline)
            out.extend(lines)
        return out

    def render(self) -> bytes:
        if self.raw is not None:
            return self.raw
        data = self.eol.join(self.lines())
        return data + self.eol if self.final else data


# ---------------------------------------------------------------------------------------------
# Bytes to insert or substitute
# ---------------------------------------------------------------------------------------------
SPECIAL = {
    "ctl01": b"\x01", "bel": b"\x07", "vt": b"\x0b", "ff": b"\x0c", "sub": b"\x1a", "us": b"\x1f",
    "del": b"\x7f", "tab": b"\t", "cr": b"\r", "nul": b"\x00", "gt": b">", "qmark": b"?", "semi": b";",
    "hyphen": b"-", "star": b"*", "digit": b"7", "digits": b"1234", "space": b" ", "spaces": b"   ",
    "nbsp": b"\xc2\xa0", "latin1": b"\xe9", "ff-byte": b"\xff", "c1": b"\x80", "utf8-e": "é".encode(),
    "utf8-cjk": "日本".encode(), "bom": b"\xef\xbb\xbf", "lf": b"\n", "hash": b"#", "bang": b"!",
    "dot": b".", "bar": b"|", "bracket": b"[a=b]", "amp": b"&amp;", "lower": b"x",
}
# Bytes that are residues of one molecule only, or letters that are ambiguity codes.
RESIDUE = {"nuc": {"U": b"U", "N": b"N", "ambig": b"RYKMSWBDHV", "X": b"X", "lower-run": b"acgt", "O": b"O"},
           "prot": {"B": b"B", "Z": b"Z", "J": b"J", "X": b"X", "U": b"U", "O": b"O", "lower-run": b"acd"}}


def pick_special(rng: random.Random, kind: str) -> tuple[str, bytes]:
    names = list(SPECIAL) + [f"res-{k}" for k in RESIDUE[kind]]
    name = rng.choice(names)
    if name.startswith("res-"):
        return name, RESIDUE[kind][name[4:]]
    return name, SPECIAL[name]


# ---------------------------------------------------------------------------------------------
# Mutations. Each function takes (rng, doc) and returns (description, doc).
# ---------------------------------------------------------------------------------------------
def defline_positions(defline: bytes) -> dict[str, int]:
    space = defline.find(b" ")
    return {"start": 0, "in-id": max(1, (space if space > 0 else len(defline)) // 2),
            "before-space": space if space > 0 else len(defline), "after-space": space + 1 if space > 0 else len(defline),
            "mid-title": (len(defline) + (space if space > 0 else 0)) // 2, "end": len(defline)}


def m_defline_insert(rng, doc, kind):
    doc = doc.copy()
    index = rng.randrange(len(doc.recs))
    defline = doc.recs[index][0]
    where = rng.choice(sorted(defline_positions(defline)))
    pos = defline_positions(defline)[where]
    name, data = pick_special(rng, kind)
    doc.recs[index][0] = defline[:pos] + data + defline[pos:]
    return f"defline{index + 1} insert {name} at {where}", doc


def m_defline_replace(rng, doc, kind):
    doc = doc.copy()
    index = rng.randrange(len(doc.recs))
    defline = doc.recs[index][0]
    where = rng.choice(sorted(defline_positions(defline)))
    pos = min(defline_positions(defline)[where], len(defline) - 1)
    name, data = pick_special(rng, kind)
    doc.recs[index][0] = defline[:pos] + data + defline[pos + 1:]
    return f"defline{index + 1} replace byte {pos} by {name}", doc


def pick_line(rng, doc):
    candidates = [(r, i) for r, (_, lines) in enumerate(doc.recs) for i in range(len(lines))]
    return rng.choice(candidates)


def m_seq_insert(rng, doc, kind):
    doc = doc.copy()
    r, i = pick_line(rng, doc)
    line = doc.recs[r][1][i]
    where = rng.choice(["start", "early", "mid", "end"])
    pos = {"start": 0, "early": min(3, len(line)), "mid": len(line) // 2, "end": len(line)}[where]
    name, data = pick_special(rng, kind)
    doc.recs[r][1][i] = line[:pos] + data + line[pos:]
    return f"rec{r + 1} line{i + 1} insert {name} at {where}", doc


def m_seq_replace(rng, doc, kind):
    doc = doc.copy()
    r, i = pick_line(rng, doc)
    line = doc.recs[r][1][i]
    pos = rng.randrange(len(line))
    name, data = pick_special(rng, kind)
    doc.recs[r][1][i] = line[:pos] + data + line[pos + 1:]
    return f"rec{r + 1} line{i + 1} replace byte {pos} by {name}", doc


def m_seq_run(rng, doc, kind):
    doc = doc.copy()
    r, i = pick_line(rng, doc)
    line = doc.recs[r][1][i]
    name, data = rng.choice([("hyphen", b"-"), ("star", b"*"), ("digit", b"5"), ("space", b" "), ("tab", b"\t"),
                             ("semi", b";"), ("nul", b"\x00"), ("X", b"X"), ("N", b"N")])
    count = rng.choice([2, 5, 20, 70])
    pos = rng.randrange(len(line) + 1)
    doc.recs[r][1][i] = line[:pos] + data * count + line[pos:]
    return f"rec{r + 1} line{i + 1} insert {count} x {name} at {pos}", doc


def with_eols(doc: Doc, eols: list[bytes], final: bytes | None) -> Doc:
    """Join the lines with the given line ends in turn (the last item repeats), ending with `final`."""
    doc = doc.copy()
    lines = doc.lines()
    out = []
    for n, line in enumerate(lines):
        out.append(line)
        if n < len(lines) - 1:
            out.append(eols[min(n, len(eols) - 1)] if len(eols) > 1 else eols[0])
    if final is not None:
        out.append(final)
    doc.raw = b"".join(out)
    return doc


def eol_cases(doc: Doc, kind: str) -> list[tuple[str, Doc]]:
    n = len(doc.lines())

    def cyc(pattern):          # a pattern of line ends that repeats along the file
        return [pattern[k % len(pattern)] for k in range(n)]

    LF, CRLF, CR = b"\n", b"\r\n", b"\r"
    out = [
        ("eol LF, no final newline", with_eols(doc, [LF], None)),
        ("eol CRLF", with_eols(doc, [CRLF], CRLF)),
        ("eol CRLF, no final newline", with_eols(doc, [CRLF], None)),
        ("eol CRLF, final LF", with_eols(doc, [CRLF], LF)),
        ("eol CR only", with_eols(doc, [CR], CR)),
        ("eol CR only, no final newline", with_eols(doc, [CR], None)),
        ("eol mixed LF/CRLF", with_eols(doc, cyc([LF, CRLF]), LF)),
        ("eol mixed CRLF/LF", with_eols(doc, cyc([CRLF, LF]), CRLF)),
        ("eol mixed LF/CR", with_eols(doc, cyc([LF, CR]), LF)),
        ("eol mixed CR/LF", with_eols(doc, cyc([CR, LF]), CR)),
        ("eol mixed CRLF/CR", with_eols(doc, cyc([CRLF, CR]), CRLF)),
        ("eol mixed LF/CRLF/CR", with_eols(doc, cyc([LF, CRLF, CR]), LF)),
        ("eol LF but CRLF only after defline", with_eols(doc, [CRLF if line.startswith(b">") else LF for line in doc.lines()], LF)),
        ("eol LF+CR (LFCR)", with_eols(doc, [b"\n\r"], b"\n\r")),
        ("eol CRCRLF", with_eols(doc, [b"\r\r\n"], b"\r\r\n")),
        ("eol LF, three blank lines at the end", with_eols(doc, [LF], b"\n\n\n")),
        ("eol CRLF, two blank CRLF lines at the end", with_eols(doc, [CRLF], b"\r\n\r\n\r\n")),
        ("eol CR, two blank lines at the end", with_eols(doc, [CR], b"\r\r\r")),
        ("eol LF, final lone CR", with_eols(doc, [LF], b"\n\r")),
        ("eol CRLF, lone CR at the very end", with_eols(doc, [CRLF], b"\r\n\r")),
    ]
    # A lone CR or LF inside one line of a file of the other kind.
    lines = doc.lines()
    mid = lines[1] if len(lines) > 1 else lines[0]
    out.append(("lone CR inside a sequence line (LF file)", with_eols(replace_line(doc, 1, mid[:10] + b"\r" + mid[10:]), [LF], LF)))
    out.append(("lone LF inside a sequence line (CRLF file)", with_eols(replace_line(doc, 1, mid[:10] + b"\n" + mid[10:]), [CRLF], CRLF)))
    out.append(("lone CR inside the defline (LF file)", with_eols(replace_line(doc, 0, lines[0][:4] + b"\r" + lines[0][4:]), [LF], LF)))
    return out


def replace_line(doc: Doc, number: int, new: bytes) -> Doc:
    """A copy whose number-th line (counting deflines) is `new`."""
    doc = doc.copy()
    seen = 0
    for rec in doc.recs:
        if rec[0] is not None:
            if seen == number:
                rec[0] = new[1:] if new.startswith(b">") else new
                return doc
            seen += 1
        for i in range(len(rec[1])):
            if seen == number:
                rec[1][i] = new
                return doc
            seen += 1
    return doc


def structure_cases(doc: Doc, kind: str) -> list[tuple[str, Doc]]:
    out: list[tuple[str, Doc]] = []

    def variant(description, fn):
        d = doc.copy()
        fn(d)
        out.append((description, d))

    variant("empty first record", lambda d: d.recs.__setitem__(0, [d.recs[0][0], []]))
    variant("empty middle record", lambda d: d.recs.__setitem__(1, [d.recs[1][0], []]))
    variant("empty last record", lambda d: d.recs.__setitem__(2, [d.recs[2][0], []]))
    variant("empty first and last record", lambda d: (d.recs.__setitem__(0, [d.recs[0][0], []]), d.recs.__setitem__(2, [d.recs[2][0], []])))
    variant("all records empty", lambda d: [d.recs.__setitem__(i, [d.recs[i][0], []]) for i in range(3)])
    variant("empty records with empty deflines", lambda d: (d.recs.__setitem__(1, [b"", []]), d.recs.__setitem__(2, [b"", []])))
    variant("extra empty record at the end", lambda d: d.recs.append([b"empty tail", []]))
    variant("only blank sequence lines in a record", lambda d: d.recs[1].__setitem__(1, [b"", b"  ", b""]))
    variant("empty defline in record 2", lambda d: d.recs[1].__setitem__(0, b""))
    variant("blank defline in record 2", lambda d: d.recs[1].__setitem__(0, b"   "))
    variant("leading space in the title of record 1", lambda d: d.recs[0].__setitem__(0, b"   lead title"))
    variant("space after > in record 2", lambda d: d.recs[1].__setitem__(0, b" " + d.recs[1][0]))
    variant("defline only id, no title", lambda d: d.recs[1].__setitem__(0, d.recs[1][0].split(b" ")[0]))
    variant("defline of one byte", lambda d: d.recs[0].__setitem__(0, b"x"))
    variant("two deflines in a row", lambda d: d.recs.insert(1, [b"extra defline", []]))
    variant("headerless first record", lambda d: d.recs.__setitem__(0, [None, d.recs[0][1]]))
    variant("headerless only record", lambda d: setattr(d, "recs", [[None, d.recs[0][1]]]))
    variant("headerless first record, blank line first", lambda d: d.recs.__setitem__(0, [None, [b""] + d.recs[0][1]]))
    variant("headerless first record starts with spaces", lambda d: d.recs.__setitem__(0, [None, [b"  " + d.recs[0][1][0]] + d.recs[0][1][1:]]))
    for text, label in ((b"", "blank line"), (b"   \t ", "white-space line"), (b"# a comment", "# comment"),
                        (b"; a comment", "; comment"), (b"! a comment", "! comment"), (b"free text line", "free text line"),
                        (b"**", "star line"), (b"12345", "digits line"), (b"\t", "tab line")):
        variant(f"{label} before the first defline", lambda d, t=text: setattr(d, "raw", b"\n".join([t] + d.lines()) + b"\n"))
    variant("two blank lines then comments before the first defline",
            lambda d: setattr(d, "raw", b"\n".join([b"", b"", b"# one", b"; two", b"! three", b""] + d.lines()) + b"\n"))
    variant("BOM before the first defline", lambda d: setattr(d, "raw", b"\xef\xbb\xbf" + b"\n".join(d.lines()) + b"\n"))
    variant("BOM before a headerless first record", lambda d: setattr(d, "raw", b"\xef\xbb\xbf" + b"\n".join(d.recs[0][1] + [b">" + d.recs[1][0]] + d.recs[1][1]) + b"\n"))
    variant("BOM and CRLF", lambda d: setattr(d, "raw", b"\xef\xbb\xbf" + b"\r\n".join(d.lines()) + b"\r\n"))
    variant("space before the first >", lambda d: setattr(d, "raw", b" " + b"\n".join(d.lines()) + b"\n"))
    # '>?' lines: NCBI reads them as a gap in the sequence.
    for text in (b">?100", b">?", b">?unk100", b">?abc", b">? 100", b">?_x pseudo", b">?5 [gap-type=within scaffold] [linkage-evidence=paired-ends]"):
        variant(f"gap line {text.decode()} inside record 2", lambda d, t=text: d.recs[1][1].insert(1, t))
    variant("gap line >?30 at the start of record 1", lambda d: d.recs[0][1].insert(0, b">?30"))
    variant("gap line >?50 at the end of record 2", lambda d: d.recs[1][1].append(b">?50"))
    variant("first record is a gap line", lambda d: d.recs.insert(0, [b"?100", []]))
    variant("two gap lines in a row", lambda d: d.recs[1][1].insert(1, b">?20") or d.recs[1][1].insert(2, b">?30"))
    variant("comment lines in a record", lambda d: d.recs[1][1].insert(1, b"#c") or d.recs[1][1].insert(3, b";c") or d.recs[1][1].append(b"!c"))
    variant("> inside a sequence line", lambda d: d.recs[1][1].__setitem__(1, d.recs[1][1][1][:10] + b">" + d.recs[1][1][1][10:]))
    variant("space then > at line start", lambda d: d.recs[1][1].insert(1, b" >q9 not a defline"))
    variant("lone > record", lambda d: d.recs.insert(1, [b"", []]))
    # Very long titles.
    for length in (999, 1000, 1001, 5000, 40000):
        for style, make in (("words", lambda n: (b"w" * 7 + b" ") * (n // 8 + 1)), ("one word", lambda n: b"L" * n)):
            title = make(length)[:length]
            variant(f"title of {length} bytes ({style}) in record 2", lambda d, t=title: d.recs[1].__setitem__(0, d.recs[1][0].split(b" ")[0] + b" " + t))
    variant("long id (120 bytes)", lambda d: d.recs[0].__setitem__(0, b"i" * 120 + b" long id title"))
    variant("title of 25 nucleotide letters", lambda d: d.recs[0].__setitem__(0, d.recs[0][0] + b" ACGTACGTACGTACGTACGTACGTA"))
    variant("title of 20 amino acid letters", lambda d: d.recs[0].__setitem__(0, d.recs[0][0] + b" ACDEFGHIKLMNPQRSTVWY"))
    variant("title only letters (id of 20 letters)", lambda d: d.recs[0].__setitem__(0, b"ACGTACGTACGTACGTACGT"))
    variant("punctuation title ', ,' in record 2", lambda d: d.recs[1].__setitem__(0, b", ,"))
    variant("punctuation title '; ;' in record 2", lambda d: d.recs[1].__setitem__(0, b"; ;"))
    variant("title with HTML entities", lambda d: d.recs[0].__setitem__(0, d.recs[0][0] + b" A &amp; B &lt;x&gt; &#65; &#x41; &#0; &#xD800;"))
    variant("title ends with a byte >= 0x80", lambda d: d.recs[1].__setitem__(0, d.recs[1][0] + b" caf\xc3\xa9"))
    variant("title ends with a lone 0xE9", lambda d: d.recs[1].__setitem__(0, d.recs[1][0] + b" caf\xe9"))
    return out


# First lines that NCBI may try as a Seq-id (data loader). Run inside `unshare -rn` only.
SEQID_FIRST_LINES = [b"AB123456", b"gb|AB123456|", b"NC_000913.3", b"0123", b"XP_123456", b"ACGT1234", b"12.5", b"lcl|x",
                     b"ref|NC_000913|", b"sp|P12345|ABC_HUMAN", b"contig1", b"A12345", b"ACGTACGT", b"xx|abc"]


def seqid_cases(doc: Doc, kind: str, rng: random.Random, count: int) -> list[tuple[str, Doc]]:
    out = []
    for line in rng.sample(SEQID_FIRST_LINES, count):
        d = doc.copy()
        d.raw = b"\n".join([line] + d.lines()) + b"\n"
        out.append((f"first line {line.decode()!r} before the first defline", d))
        d = doc.copy()
        d.recs[0][1][0] = line
        d.recs[0][0] = None
        out.append((f"headerless first record starts with {line.decode()!r}", d))
    return out[:count]


def seqid_risk(data: bytes) -> bool:
    """NCBI tries the first line as a Seq-id when, trimmed of white space, it starts with a letter or
    digit and is not letters only (CBlastInputReader::ReadOneSeq)."""
    line = re.match(rb"[^\r\n]*", data).group().strip(b" \t\n\v\f\r")
    return bool(line) and line[:1].isalnum() and re.fullmatch(rb"[A-Za-z]+", line) is None


# ---------------------------------------------------------------------------------------------
# generate
# ---------------------------------------------------------------------------------------------
CORE_CASES = ("eol CRLF", "eol CR only", "eol mixed LF/CRLF", "eol LF, no final", "headerless first record",
              "gap line >?100 inside", "empty middle record", "BOM before the first defline", "blank line before",
              "title of 1001 bytes (words)", "punctuation title ', ,'", "title with HTML", "empty defline in record 2",
              "leading space in the title", "# comment before")


def generate_files(kind: str, role: str, per_cell: int) -> list[tuple[str, str, bytes]]:
    """(family, description, bytes) for one molecule and role, deterministic and without duplicates."""
    digest = hashlib.sha256(f"{SEED}:{kind}:{role}".encode()).digest()
    rng = random.Random(int.from_bytes(digest[:8], "big"))
    clean = Doc(kind, role)
    seen = {clean.render()}
    out: list[tuple[str, str, bytes]] = []

    def add(family, description, doc, allow_risk=False):
        data = doc.render()
        if data in seen or (seqid_risk(data) and not allow_risk):
            return False
        seen.add(data)
        out.append((family, description, data))
        return True

    # Deterministic families first, then the random ones fill the rest.
    # Half of the line-end and structure cases go to each role (CORE_CASES go to both), so that the
    # matrix stays near 2300 cases; every case still meets a query and a subject reader of its program.
    parity = 0 if role == "q" else 1
    for family, cases in (("eol", eol_cases(clean, kind)), ("structure", structure_cases(clean, kind))):
        for number, (description, doc) in enumerate(cases):
            if number % 2 == parity or any(description.startswith(core) for core in CORE_CASES):
                add(family, description, doc)
    for description, doc in seqid_cases(clean, kind, rng, 6):
        add("seqid", description, doc, allow_risk=True)
    # Mutations with a fixed share of the remaining slots.
    remaining = max(0, per_cell - len(out))
    shares = [("defline-insert", m_defline_insert, 0.30), ("defline-replace", m_defline_replace, 0.10),
              ("seq-insert", m_seq_insert, 0.32), ("seq-replace", m_seq_replace, 0.14), ("seq-run", m_seq_run, 0.06)]
    for family, fn, share in shares:
        want, attempts = round(remaining * share), 0
        done = 0
        while done < want and attempts < want * 50:
            attempts += 1
            description, doc = fn(rng, clean, kind)
            if add(family, description, doc):
                done += 1
    # Combinations: a defline mutation, a sequence mutation and a line-end kind together.
    for _ in range(max(2, remaining // 12)):
        d1, doc = m_defline_insert(rng, clean, kind)
        d2, doc = m_seq_insert(rng, doc, kind)
        description, eol_doc = rng.choice(eol_cases(doc, kind))
        add("combo", f"{d1}; {d2}; {description}", eol_doc)
    return out[:per_cell]


def program_cells() -> list[tuple[str, str, str, str]]:
    """(program label, role, kind of the mutated file, kind of the clean partner)."""
    rows = []
    for label, (_, _, (qk, sk)) in PROGRAMS.items():
        rows.append((label, "q", qk, sk))
        rows.append((label, "s", sk, qk))
    return rows


def command_generate(args) -> int:
    base = Path(args.dir).resolve()
    inputs = base / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    for kind in ("nuc", "prot"):
        for role in ROLES:
            doc = Doc(kind, role)
            (inputs / f"base.{kind}.{role}.fa").write_bytes(doc.render())
    files: dict[tuple[str, str], list[tuple[str, str, str]]] = {}
    for kind in ("nuc", "prot"):
        for role in ROLES:
            directory = inputs / f"{kind}.{role}"
            directory.mkdir(exist_ok=True)
            records = []
            for n, (family, description, data) in enumerate(generate_files(kind, role, args.per_cell), 1):
                name = f"{n:03d}.fa"
                (directory / name).write_bytes(data)
                records.append((name, family, description))
            files[(kind, role)] = records
    lines = ["case_id\tprogram\ttask\trole\toutfmt\tquery\tsubject\tfamily\tmutation"]
    for label, role, mutated_kind, partner_kind in program_cells():
        executable, task, (qk, sk) = PROGRAMS[label]
        for name, family, description in files[(mutated_kind, role)]:
            mutated = f"inputs/{mutated_kind}.{role}/{name}"
            partner = f"inputs/base.{partner_kind}.{'s' if role == 'q' else 'q'}.fa"
            query, subject = (mutated, partner) if role == "q" else (partner, mutated)
            for fmt in OUTFMTS:
                case_id = f"{label}.{mutated_kind}{role}.{name[:-3]}.o{fmt}"
                lines.append("\t".join([case_id, executable, task, role, fmt, query, subject, family, description]))
    (base / "cases.tsv").write_text("\n".join(lines) + "\n")
    print(f"{len(lines) - 1} cases from {sum(len(v) for v in files.values())} files in {base}")
    return 0


# ---------------------------------------------------------------------------------------------
# Running (shared with check_inputs.py)
# ---------------------------------------------------------------------------------------------
def clean_env(extra: dict[str, str] | None = None) -> dict[str, str]:
    env = {k: v for k, v in os.environ.items()
           if k not in REPORT_ENV and not k.startswith(("LOSAT_", "RAYON_", "NCBI_CONFIG__"))}
    env.update(extra or {})
    return env


def run_both(command: list[str], cwd: Path, env: dict[str, str], stdin_data: bytes = b"") -> dict:
    """Run a command twice: stdout and stderr apart, then `2>&1`."""
    result = {}
    try:
        p = subprocess.run(command, cwd=cwd, env=env, capture_output=True, input=stdin_data, timeout=TIMEOUT)
        m = subprocess.run(command, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           input=stdin_data, timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        return {"rc": -999, "rc_merged": -999, "out": b"", "err": b"<timeout>", "merged": b""}
    result.update(rc=p.returncode, rc_merged=m.returncode, out=p.stdout, err=p.stderr, merged=m.stdout)
    return result


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def wait_for_load(limit: float = 12.0) -> None:
    while os.getloadavg()[0] > limit:
        print(f"load average {os.getloadavg()[0]:.1f} > {limit}: waiting", file=sys.stderr, flush=True)
        time.sleep(30)


def first_difference(expected: bytes, actual: bytes) -> str:
    left, right = expected.split(b"\n"), actual.split(b"\n")
    for number, (a, b) in enumerate(zip(left, right), 1):
        if a != b:
            return f"line {number}: ncbi {a[:70]!r}, losat {b[:70]!r}"
    return f"line {min(len(left), len(right))}: ncbi {len(left) - 1} lines, losat {len(right) - 1} lines"


def classify(ncbi: dict, ours: dict, program_upper: str = "") -> tuple[str, str, str]:
    """(class, detail, approval). `listed` is encoded in the class suffix: explicit-rejection or
    explicit-rejection/unlisted."""
    if -999 in (ncbi["rc"], ours["rc"], ncbi["rc_merged"], ours["rc_merged"]):
        return "timeout", "", ""
    streams_equal = (ncbi["out"] == ours["out"] and ncbi["err"] == ours["err"] and ncbi["merged"] == ours["merged"]
                     and ncbi["rc"] == ours["rc"] and ncbi["rc_merged"] == ours["rc_merged"])
    if streams_equal:
        return ("same" if ncbi["rc"] == 0 else "same-error"), "", ""
    if ncbi["rc"] < 0 and ours["rc"] == 0:
        return "approved-exception", f"ncbi signal {-ncbi['rc']}, losat exit 0", APPROVAL_EXCEPTION_2
    # An argument that NCBI's and LOSAT's argument parsers both reject, each in its own words.
    if (ncbi["rc"] == 1 and b"CArgException" in ncbi["err"] and ours["rc"] == 2 and ours["err"].startswith(b"error: ")
            and ncbi["out"] == ours["out"]):
        return "approved-exception", "argument rejected by both parsers", APPROVAL_ARG_ERROR
    if ours["rc"] != 0 and REJECT_RE.search(ours["err"]):
        lines = ours["err"].decode(errors="replace").strip().splitlines()
        first = next((l for l in lines if "not supported by LOSAT" in l), lines[0])[:200]
        listed = any(re.search(p, ours["err"]) for p in LISTED_REJECTIONS)
        return ("explicit-rejection" if listed else "explicit-rejection/unlisted"), first, ""
    if ncbi["rc"] != ours["rc"]:
        first = ours["err"].decode(errors="replace").strip().splitlines()
        return "differs", f"exit ncbi {ncbi['rc']} losat {ours['rc']}; losat stderr: {first[0][:120] if first else ''}", ""
    for name, key in (("stdout", "out"), ("stderr", "err"), ("2>&1", "merged")):
        if ncbi[key] != ours[key]:
            return "differs", f"{name} {first_difference(ncbi[key], ours[key])}", ""
    return "differs", "exit code of the 2>&1 run", ""


def load_cases(base: Path) -> list[dict]:
    rows = [line.split("\t") for line in (base / "cases.tsv").read_text().splitlines()]
    return [dict(zip(rows[0], row)) for row in rows[1:]]


def argv_of(case: dict) -> list[str]:
    argv = ["-query", case["query"], "-subject", case["subject"]]
    if case["task"]:
        argv += ["-task", case["task"]]
    return argv + ["-outfmt", case["outfmt"]]


def freeze_dir(cases: list[tuple[str, list[str], str]], base: Path, ncbi_bin: Path, out: Path, jobs: int) -> dict:
    """cases: (id, argv, executable). Stores <id>.rc/.out/.err/.merged and manifest.tsv in `out`."""
    out.mkdir(parents=True, exist_ok=True)
    env = clean_env()

    def one(item):
        case_id, argv, executable = item
        result = run_both(["unshare", "-rn", str(ncbi_bin / executable), *argv], base, env)
        (out / f"{case_id}.out").write_bytes(result["out"])
        (out / f"{case_id}.err").write_bytes(result["err"])
        (out / f"{case_id}.merged").write_bytes(result["merged"])
        (out / f"{case_id}.rc").write_text(f"{result['rc']} {result['rc_merged']}\n")
        return case_id, result

    wait_for_load()
    started = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        results = dict(pool.map(one, cases))
    manifest = ["case_id\trc\trc_merged\tstdout_sha256\tstderr_sha256\tmerged_sha256"]
    for case_id, _, _ in cases:
        r = results[case_id]
        manifest.append("\t".join([case_id, str(r["rc"]), str(r["rc_merged"]), sha(r["out"]), sha(r["err"]), sha(r["merged"])]))
    (out / "manifest.tsv").write_text("\n".join(manifest) + "\n")
    counts: dict[str, int] = {}
    for r in results.values():
        counts[f"rc{r['rc']}"] = counts.get(f"rc{r['rc']}", 0) + 1
    print(f"froze {len(cases)} cases in {time.time() - started:.1f} s: {counts}")
    return results


def command_freeze(args) -> int:
    base = Path(args.dir).resolve()
    out = Path(args.out).resolve() if args.out else base / "ncbi"
    cases = [(c["case_id"], argv_of(c), c["program"]) for c in load_cases(base)]
    freeze_dir(cases, base, Path(args.ncbi_bin), out, args.jobs)
    return 0


def read_frozen(out: Path, case_id: str) -> dict:
    rc, rc_merged = (out / f"{case_id}.rc").read_text().split()
    return {"rc": int(rc), "rc_merged": int(rc_merged), "out": (out / f"{case_id}.out").read_bytes(),
            "err": (out / f"{case_id}.err").read_bytes(), "merged": (out / f"{case_id}.merged").read_bytes()}


def check_rows(rows: list[dict], base: Path, losat: Path, frozen: Path, jobs: int, env_of=None) -> list[dict]:
    """rows: dicts with case_id, program (executable), argv, plus columns to keep."""
    env = clean_env()

    def one(row):
        ours = run_both([str(losat), row["program"], *row["argv"]], base, env)
        ncbi = read_frozen(frozen, row["case_id"])
        klass, detail, approval = classify(ncbi, ours)
        return {**row, "ncbi_rc": ncbi["rc"], "losat_rc": ours["rc"], "class": klass, "detail": detail, "approval": approval}

    wait_for_load()
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        return list(pool.map(one, rows))


COLUMNS = ["case_id", "program", "task", "role", "outfmt", "family", "ncbi_rc", "losat_rc", "class", "approval", "detail", "mutation"]


def write_tsv(results: list[dict], path: Path, columns=COLUMNS) -> None:
    def clean(value):
        return re.sub(r"[\t\r\n]+", " ", str(value))
    path.write_text("\t".join(columns) + "\n" + "".join("\t".join(clean(r.get(c, "")) for c in columns) + "\n" for r in results))


def command_check(args) -> int:
    base = Path(args.dir).resolve()
    frozen = Path(args.ncbi).resolve() if args.ncbi else base / "ncbi"
    started = time.time()
    rows = [{**c, "argv": argv_of(c)} for c in load_cases(base)]
    results = check_rows(rows, base, Path(args.losat).resolve(), frozen, args.jobs)
    write_tsv(results, Path(args.out))
    counts: dict[str, int] = {}
    for r in results:
        counts[r["class"]] = counts.get(r["class"], 0) + 1
    print(" ".join(f"{k}={v}" for k, v in sorted(counts.items())), f"cases={len(results)} seconds={time.time() - started:.1f}")
    bad = counts.get("differs", 0) + counts.get("timeout", 0) + counts.get("explicit-rejection/unlisted", 0)
    return 1 if bad else 0


CLASS_ORDER = ["same", "same-error", "explicit-rejection", "explicit-rejection/unlisted", "approved-exception", "differs", "timeout"]


def summary_table(rows: list[dict], key_columns: list[str]) -> str:
    groups: dict[tuple, dict[str, int]] = {}
    for r in rows:
        counts = groups.setdefault(tuple(r.get(c, "") for c in key_columns), {})
        counts[r["class"]] = counts.get(r["class"], 0) + 1
    header = "| " + " | ".join(key_columns + CLASS_ORDER + ["total"]) + " |"
    lines = [header, "|" + "---|" * (len(key_columns) + len(CLASS_ORDER) + 1)]
    for key in sorted(groups):
        counts = groups[key]
        lines.append("| " + " | ".join([*key, *(str(counts.get(c, 0)) for c in CLASS_ORDER), str(sum(counts.values()))]) + " |")
    total: dict[str, int] = {}
    for counts in groups.values():
        for k, v in counts.items():
            total[k] = total.get(k, 0) + v
    lines.append("| " + " | ".join(["**all**"] + [""] * (len(key_columns) - 1) + [str(total.get(c, 0)) for c in CLASS_ORDER] + [str(sum(total.values()))]) + " |")
    return "\n".join(lines)


def read_tsv(path: Path) -> list[dict]:
    rows = [line.split("\t") for line in path.read_text().splitlines()]
    return [dict(zip(rows[0], row)) for row in rows[1:]]


def command_summary(args) -> int:
    rows = [r for p in args.tsv for r in read_tsv(Path(p))]
    for r in rows:
        r["program_task"] = r["program"] + (f" -task {r['task']}" if r.get("task") else "")
    print(summary_table(rows, ["program_task", "role", "outfmt"]))
    return 0


def command_compare(args) -> int:
    a, b = (Path(p) / "manifest.tsv" for p in (args.a, args.b))
    same = a.read_bytes() == b.read_bytes()
    print(f"manifest sha256 {sha(a.read_bytes())} / {sha(b.read_bytes())}: {'identical' if same else 'DIFFERENT'}")
    files_a = sorted(p.name for p in Path(args.a).iterdir())
    files_b = sorted(p.name for p in Path(args.b).iterdir())
    if files_a != files_b:
        print("file lists differ")
        same = False
    return 0 if same else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    g = sub.add_parser("generate")
    g.add_argument("--dir", required=True)
    g.add_argument("--per-cell", type=int, default=115)
    g.set_defaults(func=command_generate)
    f = sub.add_parser("freeze-ncbi")
    f.add_argument("--dir", required=True)
    f.add_argument("--ncbi-bin", required=True)
    f.add_argument("--jobs", type=int, required=True)
    f.add_argument("--out")
    f.set_defaults(func=command_freeze)
    c = sub.add_parser("check")
    c.add_argument("--dir", required=True)
    c.add_argument("--losat", required=True)
    c.add_argument("--jobs", type=int, required=True)
    c.add_argument("--out", required=True)
    c.add_argument("--ncbi")
    c.set_defaults(func=command_check)
    s = sub.add_parser("summary")
    s.add_argument("tsv", nargs="+")
    s.set_defaults(func=command_summary)
    k = sub.add_parser("compare-frozen")
    k.add_argument("a")
    k.add_argument("b")
    k.set_defaults(func=command_compare)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.dont_write_bytecode = True
    sys.exit(main())
