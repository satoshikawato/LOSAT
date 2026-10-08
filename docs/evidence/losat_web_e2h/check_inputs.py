#!/usr/bin/env python3
"""Compare the input handling of the four programs with NCBI (Session SFb, E2h; comparison only).

This is docs/evidence/losat_web_e2g/check_inputs.py (which builds on E2f's and E2c's) with the
expectations of E2h. What changed:

- The cases that E2g recorded as `losat-rejects` because LOSAT did not read the FASTA file as NCBI
  does (a tab, a leading space, `X`, `-` or digits in a file, an empty or blank defline, an empty
  record, text before the first defline, a byte order mark, a gap line, a no-break space, a title
  that ends in a space, non-ASCII titles, standard input given twice, an empty pipe) now expect what
  NCBI does: `same` when NCBI succeeds, `same-error` when it fails (the expectation `fixed`, decided
  by the frozen NCBI exit code). The cases that stay `explicit-rejection` are the ones about options
  (scoring limits, `-outfmt` delimiters, `-strand`, ...), which are not about reading the FASTA file.
- `arg-error` (NCBI's and LOSAT's argument parsers both reject, each in its own words: approved
  exception 1 of PD-LOSAT-CLI-NONSEARCH-DIFFERENCES) and `exception-2` (NCBI dies of a signal on a
  punctuation title with hits, LOSAT reports: approved exception 2 of PD-LOSAT-NCBI-DEFECTS) are the
  class `approved-exception`.
- The cases run for TBLASTX, TBLASTN (protein query, nucleotide subject) and BLASTP too: the same
  input conditions on the query and on the subject (e2h.* rows), with outfmt 0 and 6.
- NCBI's side is frozen once (`freeze-ncbi`, inside `unshare -rn`, stdout, stderr, exit code and the
  `2>&1` stream of a second run), so `check` runs only LOSAT. `check` compares stdout, stderr, the
  exit code and the `2>&1` stream (E2g compared stdout and stderr).

Classes: same, same-error, explicit-rejection, approved-exception, differs, timeout (see
fasta_sweep.py). A result that is not the expected class is listed as unexpected and the exit status
is 1.

Usage:
  check_inputs.py freeze-ncbi --ncbi-bin DIR --work DIR --out DIR --jobs N
  check_inputs.py check --losat BIN --work DIR --ncbi DIR --out TSV --jobs N
  check_inputs.py summary TSV [TSV ...]
"""
from __future__ import annotations

import argparse
import concurrent.futures
import importlib.util
import os
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fasta_sweep as fs  # noqa: E402

# E2g's script lives in docs/evidence/losat_web_e2g/; when this file is in docs/evidence/losat_web_e2h/
# the evidence directory is the parent, otherwise the engine worktree's.
EVIDENCE = HERE.parent if (HERE.parent / "losat_web_e2g").is_dir() else Path(os.environ["WT"]) / "docs/evidence"
SPEC = importlib.util.spec_from_file_location("e2g_check_inputs", EVIDENCE / "losat_web_e2g" / "check_inputs.py")
e2g = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(e2g)
e2c = e2g.e2c
ENGINE = e2c.ENGINE

# --------------------------------------------------------------------------------------------
# Expectations
# --------------------------------------------------------------------------------------------
# E2g `losat-rejects` cases that are about reading the FASTA file or the input stream: E2h expects
# NCBI's result (`fixed`).
FIXED = {
    "tab_after_id", "tab_in_title", "leading_space", "x_residue", "hyphen", "digits", "utf8_subject",
    "audit.empty_query.subject_x", "audit.empty_defline", "audit.blank_defline", "audit.empty_record",
    "audit.empty_records_only", "audit.empty_subject_record", "audit.piped_empty_query", "audit.leading_blank_line",
    "audit.byte_order_mark", "audit3.piped_white_space_query", "audit4.tab_subject",
    "audit4.leading_blank_subject.penalty_0", "audit4.stdin_both.file", "audit4.stdin_both.pipe",
    "audit4.fifo.empty_query", "audit5.empty_record_subject", "audit9.gap_line.subject", "audit9.gap_line.query",
    "audit9.gap_line.underscore", "audit9.title_warning.trailing_space", "audit10.nbsp.query",
    "audit10.ideographic_space.subject", "audit10.title_9.fmt0", "audit12.bom_query",
}
# The E2g `losat-rejects` cases that stay explicit rejections (options, not the FASTA reader):
# audit.invalid_first_batch.*, audit.invalid_strand.3_-1, audit.reward_*, audit.megablast_gap_max,
# audit2.reward_16bit_*, audit2.greedy_gap_limit, audit3.evalue_hex, audit4.outfmt.delimiter,
# audit5.outfmt.tabular_delim, audit5.hitlist.negative, audit5.threads.above_rayon, audit6.task.rmblastn,
# audit6.option.*, e2g.t7.batch_size_text.
E2H_EXPECT = {name: "fixed" for name in FIXED}
CLASS_OF_E2G = {"same": "same", "same-error": "same-error", "losat-rejects": "explicit-rejection",
                "arg-error": "approved-exception", "exception-2": "approved-exception", "fixed": "fixed", "auto": "auto"}


@dataclass
class Row:
    name: str
    program: str
    argv: list[str]
    expect: str                      # same | same-error | explicit-rejection | approved-exception | fixed | auto
    env: dict = field(default_factory=dict)
    stdin: object = None             # None (empty pipe), str (file piped in), ("file", path)
    fifos: list = field(default_factory=list)
    out_fifo: str | None = None
    role: str = "-"

    @property
    def task(self) -> str:
        return self.argv[self.argv.index("-task") + 1] if "-task" in self.argv else ""

    @property
    def outfmt(self) -> str:
        if "-outfmt" not in self.argv:
            return "0"
        words = self.argv[self.argv.index("-outfmt") + 1].split()
        return words[0] if words else ""


# --------------------------------------------------------------------------------------------
# New inputs and cases for TBLASTX, TBLASTN and BLASTP
# --------------------------------------------------------------------------------------------
def doc_of(kind: str, role: str) -> "fs.Doc":
    return fs.Doc(kind, role)


def ident(d) -> bytes:
    return d.recs[0][0].split(b" ")[0]


def _set(index, fn):
    def apply(d):
        d.recs[index][0] = fn(d.recs[index][0], d)
    return apply


def _line(rec, line, fn):
    def apply(d):
        d.recs[rec][1][line] = fn(d.recs[rec][1][line])
    return apply


def _raw(fn):
    def apply(d):
        d.raw = fn(d)
    return apply


CRLF_NO_FINAL = lambda d: b"\r\n".join(d.lines())  # noqa: E731
# (name, edit of a clean Doc, expectation by (role, outfmt) -> default "auto")
CONDITIONS = [
    ("tab_after_id", _set(0, lambda t, d: t.split(b" ")[0] + b"\tsecond field")),
    ("tab_in_title", _set(0, lambda t, d: t + b" a\tb c")),
    ("leading_space_title", _set(0, lambda t, d: b"   " + t)),
    ("x_residue", _line(1, 1, lambda b: b[:10] + b"XX" + b[10:])),
    ("hyphen", _line(1, 1, lambda b: b[:10] + b"--" + b[10:])),
    ("digits_and_spaces", _line(1, 0, lambda b: b"1 " + b[:30] + b" " + b"61 " + b[30:])),
    ("star_in_sequence", _line(1, 1, lambda b: b[:10] + b"*" + b[10:])),
    ("tab_in_sequence", _line(1, 1, lambda b: b[:10] + b"\t" + b[10:])),
    ("semicolon_comment_in_line", _line(1, 1, lambda b: b[:10] + b";comment")),
    ("invalid_text_line", lambda d: d.recs[1][1].insert(1, b"$%&() <> not a sequence")),
    ("nul_in_sequence", _line(1, 1, lambda b: b[:10] + b"\x00" + b[10:])),
    ("lowercase_sequence", lambda d: d.recs[1].__setitem__(1, [x.lower() for x in d.recs[1][1]])),
    ("utf8_title", _set(1, lambda t, d: t + " ümlaut 日本".encode())),
    ("non_utf8_title", _set(0, lambda t, d: t + b" \xe9\xff")),
    ("high_bit_title_end", _set(1, lambda t, d: t + b" caf\xc3\xa9")),
    ("control_in_title", _set(1, lambda t, d: t + b" ctl\x01tail")),
    ("nul_in_title", _set(1, lambda t, d: t.replace(b" ", b"\x00zz title ", 1))),
    ("empty_defline", _set(0, lambda t, d: b"")),
    ("blank_defline", _set(0, lambda t, d: b"   ")),
    ("empty_first_record", lambda d: d.recs.__setitem__(0, [d.recs[0][0], []])),
    ("empty_middle_record", lambda d: d.recs.__setitem__(1, [d.recs[1][0], []])),
    ("empty_last_record", lambda d: d.recs.__setitem__(2, [d.recs[2][0], []])),
    ("all_records_empty", lambda d: [d.recs.__setitem__(i, [d.recs[i][0], []]) for i in range(3)]),
    ("headerless_first_record", lambda d: d.recs.__setitem__(0, [None, d.recs[0][1]])),
    ("leading_blank_line", _raw(lambda d: b"\n" + d.render())),
    ("comment_before_first_defline", _raw(lambda d: b"# a comment\n; another\n" + d.render())),
    ("byte_order_mark", _raw(lambda d: b"\xef\xbb\xbf" + d.render())),
    ("gap_line", lambda d: d.recs[1][1].insert(1, b">?20")),
    ("gap_line_unknown", lambda d: d.recs[1][1].insert(1, b">?unk30")),
    ("gap_line_underscore", lambda d: d.recs[1][1].insert(1, b">?_x pseudo")),
    ("nbsp_in_sequence", _line(1, 1, lambda b: b[:10] + b"\xc2\xa0" + b[10:])),
    ("ideographic_space_line_at_end", _raw(lambda d: d.render() + "　\n".encode())),
    ("title_20_letters_trailing_space", _set(0, lambda t, d: t + b" ACGTACGTACGTACGTACGTA ")),
    ("title_20_letters", _set(0, lambda t, d: t + b" ACGTACGTACGTACGTACGT")),
    ("title_25_amino_acids", _set(0, lambda t, d: t + b" ACDEFGHIKLMNPQRSTVWYACDEF")),
    ("title_1001_bytes", _set(1, lambda t, d: t.split(b" ")[0] + b" " + (b"w" * 7 + b" ") * 126)),
    ("long_id_120", _set(0, lambda t, d: b"i" * 120 + b" long id")),
    ("punctuation_title", _set(1, lambda t, d: b", ,")),
    ("punctuation_title_semicolon", _set(1, lambda t, d: b"; ;")),
    ("html_entities_title", _set(0, lambda t, d: t + b" A &amp; B &lt;x&gt; &#65; &#x41; &#0;")),
    ("crlf", _raw(lambda d: b"\r\n".join(d.lines()) + b"\r\n")),
    ("crlf_no_final_newline", _raw(CRLF_NO_FINAL)),
    ("cr_only", _raw(lambda d: b"\r".join(d.lines()) + b"\r")),
    ("mixed_eol", _raw(lambda d: b"".join(line + (b"\r\n" if n % 3 == 0 else b"\n" if n % 3 == 1 else b"\r") for n, line in enumerate(d.lines())))),
    ("lone_cr_in_sequence", _raw(lambda d: b"\n".join(d.lines()[:2] + [d.lines()[2][:10] + b"\r" + d.lines()[2][10:]] + d.lines()[3:]) + b"\n")),
    ("no_final_newline", _raw(lambda d: b"\n".join(d.lines()))),
    ("blank_lines_between_records", lambda d: d.recs[0][1].append(b"") or d.recs[1][1].append(b"")),
    ("two_deflines_in_a_row", lambda d: d.recs.insert(1, [b"extra defline", []])),
]
# Subject titles of the first word with a non-UTF-8 byte in a tabular report: NCBI writes a partial
# row and dies (AUTHORITY.md section K-1, J-6): LOSAT keeps an explicit rejection of outfmt 6.
EXPECT_OVERRIDE = {("non_utf8_title", "s", "6"): "auto", ("non_utf8_title_id", "s", "6"): "explicit-rejection"}
NONUTF8_ID_TITLE = ("non_utf8_title_id", _set(0, lambda t, d: b"\x81bad" + b" " + t.split(b" ", 1)[1]))

NEW_PROGRAMS = (("tblastx", "nuc", "nuc"), ("tblastn", "prot", "nuc"), ("blastp", "prot", "prot"))


def new_input_path(work: Path, condition: str, kind: str, role: str) -> Path:
    return work / "e2h" / f"{condition}.{role}.{kind}.fa"


def make_inputs(work: Path) -> None:
    e2g.e2c.make_inputs(work)
    (work / "e2h").mkdir(exist_ok=True)
    (work / "outs").mkdir(exist_ok=True)
    for kind in ("nuc", "prot"):
        for role in ("q", "s"):
            (work / "e2h" / f"base.{role}.{kind}.fa").write_bytes(fs.Doc(kind, role).render())
            for name, edit in [*CONDITIONS, NONUTF8_ID_TITLE]:
                d = fs.Doc(kind, role)
                edit(d)
                new_input_path(work, name, kind, role).write_bytes(d.render())
    (work / "e2h" / "empty.fa").write_bytes(b"")
    (work / "e2h" / "white_space.fa").write_bytes(b"  \n\t\n")
    (work / "e2h" / "comments_only.fa").write_bytes(b"# only a comment\n; and another\n")
    (work / "e2h" / "defline_only.fa").write_bytes(b">only a defline\n")
    (work / "e2h" / "gap_only.fa").write_bytes(b">q\n>?100\n")


def new_rows(work: Path) -> list[Row]:
    rows: list[Row] = []
    w = f"{work}/e2h"
    for program, qkind, skind in NEW_PROGRAMS:
        base_q, base_s = f"{w}/base.q.{qkind}.fa", f"{w}/base.s.{skind}.fa"
        for name, _ in [*CONDITIONS, NONUTF8_ID_TITLE]:
            for role in ("q", "s"):
                kind = qkind if role == "q" else skind
                path = str(new_input_path(work, name, kind, role))
                inputs = ["-query", path, "-subject", base_s] if role == "q" else ["-query", base_q, "-subject", path]
                for fmt in ("0", "6"):
                    if name == "non_utf8_title_id" and role == "q":
                        continue
                    expect = EXPECT_OVERRIDE.get((name, role, fmt), "auto")
                    rows.append(Row(f"e2h.{program}.{name}.{role}.o{fmt}", program, [*inputs, "-outfmt", fmt], expect, role=role))
        # Files that are empty or have no residue, and the standard input in its forms.
        for role in ("q", "s"):
            for file in ("empty", "white_space", "comments_only", "defline_only", "gap_only"):
                inputs = ["-query", f"{w}/{file}.fa", "-subject", base_s] if role == "q" else ["-query", base_q, "-subject", f"{w}/{file}.fa"]
                for fmt in ("0", "6"):
                    rows.append(Row(f"e2h.{program}.file_{file}.{role}.o{fmt}", program, [*inputs, "-outfmt", fmt], "auto", role=role))
        for fmt in ("0", "6"):
            o = ["-outfmt", fmt]
            rows += [
                Row(f"e2h.{program}.stdin_query_pipe.o{fmt}", program, ["-query", "-", "-subject", base_s, *o], "auto", stdin=base_q, role="q"),
                Row(f"e2h.{program}.stdin_query_file.o{fmt}", program, ["-query", "-", "-subject", base_s, *o], "auto", stdin=("file", base_q), role="q"),
                Row(f"e2h.{program}.dev_stdin_query.o{fmt}", program, ["-query", "/dev/stdin", "-subject", base_s, *o], "auto", stdin=base_q, role="q"),
                Row(f"e2h.{program}.stdin_subject_pipe.o{fmt}", program, ["-query", base_q, "-subject", "-", *o], "auto", stdin=base_s, role="s"),
                Row(f"e2h.{program}.stdin_subject_file.o{fmt}", program, ["-query", base_q, "-subject", "-", *o], "auto", stdin=("file", base_s), role="s"),
                Row(f"e2h.{program}.stdin_empty_query.o{fmt}", program, ["-query", "-", "-subject", base_s, *o], "auto", role="q"),
                Row(f"e2h.{program}.stdin_empty_subject.o{fmt}", program, ["-query", base_q, "-subject", "-", *o], "auto", role="s"),
                Row(f"e2h.{program}.stdin_both_pipe.o{fmt}", program, ["-query", "-", "-subject", "-", *o], "auto", stdin=base_s, role="q"),
                Row(f"e2h.{program}.missing_query.o{fmt}", program, ["-query", f"{w}/missing.fa", "-subject", base_s, *o], "auto", role="q"),
                Row(f"e2h.{program}.missing_subject.o{fmt}", program, ["-query", base_q, "-subject", f"{w}/missing.fa", *o], "auto", role="s"),
                Row(f"e2h.{program}.directory_query.o{fmt}", program, ["-query", w, "-subject", base_s, *o], "auto", role="q"),
                Row(f"e2h.{program}.directory_subject.o{fmt}", program, ["-query", base_q, "-subject", w, *o], "auto", role="s"),
                Row(f"e2h.{program}.fifo_query.o{fmt}", program, ["-query", f"{work}/fifo_q", "-subject", base_s, *o], "auto",
                    fifos=[(f"{work}/fifo_q", base_q)], role="q"),
                Row(f"e2h.{program}.fifo_subject_then_query.o{fmt}", program, ["-query", f"{work}/fifo_q", "-subject", f"{work}/fifo_s", *o], "auto",
                    fifos=[(f"{work}/fifo_s", base_s), (f"{work}/fifo_q", base_q)], role="q"),
            ]
        # The pairs of the clean files, and -out.
        rows.append(Row(f"e2h.{program}.clean.o7", program, ["-query", base_q, "-subject", base_s, "-outfmt", "7"], "same"))
        rows.append(Row(f"e2h.{program}.clean_out_file", program, ["-query", base_q, "-subject", base_s, "-out", "{OUT}"], "same"))
    return rows


def all_rows(work: Path) -> list[Row]:
    rows = []
    for name, argv, expect, *extra in e2g.cases(work):
        expect = E2H_EXPECT.get(name, CLASS_OF_E2G[expect])
        env = extra[0] if extra else {}
        source = extra[1] if len(extra) > 1 else None
        fifos = extra[2] if len(extra) > 2 else []
        out_fifo = extra[3] if len(extra) > 3 else None
        rows.append(Row(name, "blastn", list(argv), expect, env, source, fifos, out_fifo))
    rows += new_rows(work)
    names = [r.name for r in rows]
    assert len(names) == len(set(names)), "duplicate case names"
    return rows


# --------------------------------------------------------------------------------------------
# Running
# --------------------------------------------------------------------------------------------
FIFO_LOCK = threading.Lock()


def execute(prefix: list[str], row: Row, work: Path, index: int, merged: bool) -> dict:
    """One run of a row. `{OUT}` is the -out file (per case, the same path for NCBI and LOSAT); its
    bytes join the stdout; a named pipe for -out is emptied by a reader."""
    out = work / "outs" / f"{index}.out"
    out.unlink(missing_ok=True)
    argv = [str(out) if word == "{OUT}" else word for word in row.argv]
    env = fs.clean_env(row.env)
    source = row.stdin
    stdin_file = source[1] if isinstance(source, tuple) else None
    stdin_data = (ENGINE / source).read_bytes() if isinstance(source, str) else b""
    writer = reader = None
    if row.fifos:
        for fifo, _ in row.fifos:
            Path(fifo).unlink(missing_ok=True)
            os.mkfifo(fifo)
        script = " && ".join(f"cat '{ENGINE / src}' > '{fifo}'" for fifo, src in row.fifos)
        writer = subprocess.Popen(["sh", "-c", script])
    if row.out_fifo:
        Path(row.out_fifo).unlink(missing_ok=True)
        Path(f"{row.out_fifo}.read").unlink(missing_ok=True)
        os.mkfifo(row.out_fifo)
        reader = subprocess.Popen(["sh", "-c", f"cat '{row.out_fifo}' > '{row.out_fifo}.read'"])
    handle = open(ENGINE / stdin_file, "rb") if stdin_file else None
    try:
        kw = dict(cwd=ENGINE, env=env, timeout=fs.TIMEOUT)
        kw.update(stdin=handle) if handle else kw.update(input=stdin_data)
        if merged:
            run = subprocess.run([*prefix, *argv], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, **kw)
            stdout, stderr = run.stdout, b""
        else:
            run = subprocess.run([*prefix, *argv], capture_output=True, **kw)
            stdout, stderr = run.stdout, run.stderr
        rc = run.returncode
    except subprocess.TimeoutExpired:
        rc, stdout, stderr = -999, b"", b"<timeout>"
    finally:
        if handle:
            handle.close()
        for process in (writer,):
            if process:
                process.kill()
                process.wait()
    if reader:
        try:
            reader.wait(timeout=10)
        except subprocess.TimeoutExpired:
            reader.kill()
            reader.wait()
        stdout += b"<fifo>" + Path(f"{row.out_fifo}.read").read_bytes()
    if "{OUT}" in row.argv:
        stdout += b"<out>" + (out.read_bytes() if out.exists() else b"<missing>")
    return {"rc": rc, "stdout": stdout, "stderr": stderr}


def run_row(prefix: list[str], row: Row, work: Path, index: int) -> dict:
    lock = FIFO_LOCK if (row.fifos or row.out_fifo) else None
    if lock:
        lock.acquire()
    try:
        apart = execute(prefix, row, work, index, False)
        joined = execute(prefix, row, work, index, True)
    finally:
        if lock:
            lock.release()
    return {"rc": apart["rc"], "rc_merged": joined["rc"], "out": apart["stdout"], "err": apart["stderr"], "merged": joined["stdout"]}


def hash_work(work: Path) -> dict[str, str]:
    """sha256 of every regular input file of the work directory (and of the engine files that the
    rows read), so that `check` can tell that the inputs are the frozen ones."""
    hashes = {}
    for path in sorted(work.rglob("*")):
        if (path.is_file() and "outs" not in path.relative_to(work).parts and not path.name.endswith(".read")
                and set(path.name) != {"o"}):          # (the 255-byte -out file that a case writes)
            hashes[str(path.relative_to(work))] = fs.sha(path.read_bytes())
    return hashes


def engine_inputs(rows: list[Row]) -> dict[str, str]:
    hashes = {}
    for row in rows:
        for word in row.argv:
            path = ENGINE / word
            if not word.startswith("/") and path.is_file():
                hashes[word] = fs.sha(path.read_bytes())
    return dict(sorted(hashes.items()))


def prepare(work: Path) -> list[Row]:
    work.mkdir(parents=True, exist_ok=True)
    make_inputs(work)
    return all_rows(work)


def command_freeze(args) -> int:
    work, out = Path(args.work).resolve(), Path(args.out).resolve()
    rows = prepare(work)
    inputs = {f"work/{k}": v for k, v in hash_work(work).items()} | {f"engine/{k}": v for k, v in engine_inputs(rows).items()}
    out.mkdir(parents=True, exist_ok=True)
    ncbi_bin = Path(args.ncbi_bin)
    started = time.time()
    fs.wait_for_load()

    def one(item):
        index, row = item
        result = run_row(["unshare", "-rn", str(ncbi_bin / row.program)], row, work, index)
        for suffix, key in ((".out", "out"), (".err", "err"), (".merged", "merged")):
            (out / f"{row.name}{suffix}").write_bytes(result[key])
        (out / f"{row.name}.rc").write_text(f"{result['rc']} {result['rc_merged']}\n")
        return row.name, result

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = dict(pool.map(one, enumerate(rows)))
    manifest = ["case_id\trc\trc_merged\tstdout_sha256\tstderr_sha256\tmerged_sha256"]
    for row in rows:
        r = results[row.name]
        manifest.append("\t".join([row.name, str(r["rc"]), str(r["rc_merged"]), fs.sha(r["out"]), fs.sha(r["err"]), fs.sha(r["merged"])]))
    (out / "manifest.tsv").write_text("\n".join(manifest) + "\n")
    (out / "inputs.sha256").write_text("".join(f"{v}  {k}\n" for k, v in inputs.items()) + f"# work {work}\n")
    print(f"froze {len(rows)} cases in {time.time() - started:.1f} s")
    return 0


def expected_class(row: Row, ncbi: dict) -> str:
    if row.expect in ("fixed", "auto"):
        if ncbi["rc"] < 0:
            return "approved-exception"
        return "same" if ncbi["rc"] == 0 else "same-error"
    return row.expect


def command_check(args) -> int:
    work, frozen = Path(args.work).resolve(), Path(args.ncbi).resolve()
    rows = prepare(work)
    recorded = {}
    for line in (frozen / "inputs.sha256").read_text().splitlines():
        if line.startswith("# work "):
            if line[7:] != str(work):
                print(f"WARNING: frozen with --work {line[7:]}, checking with {work}")
        else:
            digest, _, name = line.partition("  ")
            recorded[name] = digest
    now = {f"work/{k}": v for k, v in hash_work(work).items()} | {f"engine/{k}": v for k, v in engine_inputs(rows).items()}
    changed = sorted(k for k in set(recorded) | set(now) if recorded.get(k) != now.get(k))
    if changed:
        print(f"WARNING: {len(changed)} input files differ from the frozen ones, e.g. {changed[:3]}")
    losat = Path(args.losat).resolve()
    started = time.time()
    fs.wait_for_load()

    def one(item):
        index, row = item
        ours = run_row([str(losat), row.program], row, work, index)
        ncbi = fs.read_frozen(frozen, row.name)
        klass, detail, approval = fs.classify(ncbi, ours)
        wanted = expected_class(row, ncbi)
        # E2g's markers: an option that LOSAT does not support says so without naming a program.
        if (klass == "differs" and wanted == "explicit-rejection" and ours["rc"] != 0
                and any(m in ours["err"] for m in (b"not supported by LOSAT", b"which LOSAT does not reproduce"))):
            klass = "explicit-rejection"
        base = klass.split("/")[0]
        # An unlisted rejection is fine where the expectation is an explicit rejection (an option).
        ok = base == wanted
        return {"case_id": row.name, "program": row.program, "task": row.task, "role": row.role, "outfmt": row.outfmt,
                "expect": wanted, "ncbi_rc": ncbi["rc"], "losat_rc": ours["rc"], "class": base if wanted == "explicit-rejection" else klass,
                "approval": approval, "unexpected": "" if ok else "UNEXPECTED", "detail": detail}

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(one, enumerate(rows)))
    columns = ["case_id", "program", "task", "role", "outfmt", "expect", "ncbi_rc", "losat_rc", "class", "approval", "unexpected", "detail"]
    fs.write_tsv(results, Path(args.out), columns)
    counts: dict[str, int] = {}
    for r in results:
        counts[r["class"]] = counts.get(r["class"], 0) + 1
    unexpected = [r["case_id"] for r in results if r["unexpected"]]
    print(" ".join(f"{k}={v}" for k, v in sorted(counts.items())), f"cases={len(results)} unexpected={len(unexpected)} seconds={time.time() - started:.1f}")
    for name in unexpected[:20]:
        print("unexpected:", name)
    return 1 if unexpected or changed else 0


def command_summary(args) -> int:
    return fs.command_summary(args)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    f = sub.add_parser("freeze-ncbi")
    f.add_argument("--ncbi-bin", required=True)
    f.add_argument("--work", required=True)
    f.add_argument("--out", required=True)
    f.add_argument("--jobs", type=int, required=True)
    f.set_defaults(func=command_freeze)
    c = sub.add_parser("check")
    c.add_argument("--losat", required=True)
    c.add_argument("--work", required=True)
    c.add_argument("--ncbi", required=True)
    c.add_argument("--out", required=True)
    c.add_argument("--jobs", type=int, required=True)
    c.set_defaults(func=command_check)
    s = sub.add_parser("summary")
    s.add_argument("tsv", nargs="+")
    s.set_defaults(func=command_summary)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
