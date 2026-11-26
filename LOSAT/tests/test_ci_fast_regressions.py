#!/usr/bin/env python3
"""Tests for ci_fast_regressions.py and frozen_allowlist.py (no LOSAT binary is run)."""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import blastn_regression_fixtures as fixtures  # noqa: E402
import ci_fast_regressions as fast  # noqa: E402
import db_search_fixtures as db_fixtures  # noqa: E402
import fasta_input_fixtures as fasta_input  # noqa: E402
import range_regression_fixtures as range_fixtures  # noqa: E402
import tblastx_regression_fixtures as tblastx_fixtures  # noqa: E402
import xinclude_fixtures as xinclude  # noqa: E402
from frozen_allowlist import Allowlist  # noqa: E402

sys.path.insert(0, str(fast.FIXTURE_DIR))
import check_losat  # noqa: E402  (the -subject fixtures' checker)
import run_oracle  # noqa: E402  (the outfmt 0 manifest's reader and writer)

FROZEN, KNOWN, OTHER = "f" * 64, "a" * 64, "0" * 64


def allowlist() -> Allowlist:
    document = {"schema": "losat-frozen-mismatch-allowlist-v1", "entries": [
        {"program": "blastn", "case_id": "known", "frozen_sha256": FROZEN, "allowed_actual_sha256": KNOWN}]}
    handle = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
    json.dump(document, handle)
    handle.close()
    return Allowlist(Path(handle.name))


def row(case_id: str, output: str, expected: str = FROZEN, program: str = "blastn") -> dict[str, str]:
    return {"program": program, "case_id": case_id, "command": "[]", "exit": "0", "output_sha256": output,
            "stderr_sha256": "e", "expected_sha256": expected, "expected_source": "gate_a"}


class SamClRuleTests(unittest.TestCase):
    """The one rule of the SAM (outfmt 17) rows: the CL tag of the @PG line is not compared."""

    WITH = b"@HD\tVN:1.2\n@SQ\tSN:Query_1\tLN:3\n@PG\tID:0\tVN:2.17.0+\tCL:/x/blastn -a \tPN:blastn\nr\t0\n"
    WITHOUT = b"@HD\tVN:1.2\n@SQ\tSN:Query_1\tLN:3\n@PG\tID:0\tVN:2.17.0+\tPN:blastn\nr\t0\n"

    def test_a_row_that_differs_in_the_command_line_only_is_equal(self):
        other = self.WITH.replace(b"/x/blastn -a ", b"losat blastn -b -num_threads 4 ")
        self.assertEqual(db_fixtures.strip_sam_cl(self.WITH), self.WITHOUT)
        self.assertEqual(db_fixtures.strip_sam_cl(other), self.WITHOUT)
        self.assertEqual(db_fixtures.strip_sam_cl(self.WITHOUT), self.WITHOUT)

    def test_losats_version_form_is_stripped_the_same_way(self):
        # LOSAT states `VN:LOSAT 0.2.0` (AGENTS.md rule 8); NCBI's frozen bytes in LOSAT's words do too.
        losat = self.WITH.replace(b"VN:2.17.0+", b"VN:" + db_fixtures.losat_wording.LOSAT_VERSION)
        expected = self.WITHOUT.replace(b"VN:2.17.0+", b"VN:" + db_fixtures.losat_wording.LOSAT_VERSION)
        self.assertEqual(db_fixtures.strip_sam_cl(losat), expected)
        self.assertEqual(db_fixtures.strip_sam_cl(losat.replace(b"/x/blastn -a ", b"losat blastn -b ")), expected)

    def test_any_other_byte_is_compared(self):
        for old, new in ((b"r\t0", b"r\t1"), (b"LN:3", b"LN:4"), (b"PN:blastn", b"PN:blastx"),
                         (b"VN:2.17.0+", b"VN:2.17.1+"), (b"@HD", b"@XX")):
            changed = self.WITH.replace(old, new)
            self.assertNotEqual(db_fixtures.strip_sam_cl(changed), self.WITHOUT, (old, new))

    def test_a_freeze_keeps_a_frozen_sam_file_that_differs_in_the_command_line_only(self):
        # A freeze with another NCBI build (another argv[0] in CL) must leave the frozen file and its cells
        # unchanged, so the release's set E (freeze into a clean checkout, then `git status`) stays clean.
        other_build = self.WITH.replace(b"/x/blastn -a ", b"/opt/ncbi/bin/blastn -a ")
        with tempfile.TemporaryDirectory() as tmp, unittest.mock.patch.object(db_fixtures, "FROZEN", Path(tmp)):
            (Path(tmp) / "r.17.out").write_bytes(self.WITH)
            sam = {"row_id": "r.17", "outfmt": "17 SQ"}
            self.assertEqual(db_fixtures.frozen_stdout(sam, other_build), self.WITH)
            changed = other_build.replace(b"r\t0", b"r\t1")
            self.assertEqual(db_fixtures.frozen_stdout(sam, changed), changed)
            self.assertEqual(db_fixtures.frozen_stdout({"row_id": "new.17", "outfmt": "17"}, other_build),
                             other_build)
            (Path(tmp) / "t.6.out").write_bytes(self.WITH)
            self.assertEqual(db_fixtures.frozen_stdout({"row_id": "t.6", "outfmt": "6"}, other_build), other_build)

    def test_a_frozen_file_other_than_sam_is_replaced_by_ncbi_s_stdout(self):
        # Outside outfmt 17 the freeze writes NCBI's stdout as it is, also where it is the start of the frozen
        # file or the frozen file is its start (no prefix rule: a run that ends by a signal stops the freeze).
        prolog = b"BLASTX 2.17.0+\n\n\nReference: ...\n\nDatabase: x2g\n"
        with tempfile.TemporaryDirectory() as tmp, unittest.mock.patch.object(db_fixtures, "FROZEN", Path(tmp)):
            (Path(tmp) / "c.1.out").write_bytes(prolog)
            row = {"row_id": "c.1", "outfmt": "1"}
            for stdout in (prolog[:17], prolog + b"\nQuery= q1\n", prolog.replace(b"x2g", b"x3g"), b""):
                self.assertEqual(db_fixtures.frozen_stdout(row, stdout), stdout)

    def test_the_rule_applies_to_outfmt_17_only(self):
        for outfmt in ("17", "17 SQ", " 17 SR ", "17 SQ SR"):
            self.assertTrue(db_fixtures.is_sam_outfmt(outfmt), outfmt)
        for outfmt in ("0", "6", "5", "16", "7 qseqid", "170", ""):
            self.assertFalse(db_fixtures.is_sam_outfmt(outfmt), outfmt)


class AllowlistTests(unittest.TestCase):
    def test_classify(self):
        listed = allowlist()
        self.assertEqual(listed.classify("blastn", "known", FROZEN, KNOWN), "allowed")
        self.assertEqual(listed.classify("blastn", "known", FROZEN, FROZEN), "stale")
        self.assertEqual(listed.classify("blastn", "known", FROZEN, OTHER), "mismatch")
        self.assertEqual(listed.classify("blastn", "known", FROZEN, KNOWN, executed=False), "mismatch")
        self.assertEqual(listed.classify("blastn", "known", OTHER, OTHER), "match")
        self.assertEqual(listed.classify("blastn", "other", FROZEN, FROZEN), "match")
        self.assertEqual(listed.classify("blastn", "other", FROZEN, KNOWN), "mismatch")
        self.assertEqual(listed.unexecuted({("blastn", "other")}, {"blastn"}), [("blastn", "known")])
        self.assertEqual(listed.unexecuted(set(), {"blastp"}), [])

    def test_repository_allowlist_is_empty_and_gate_a_reads_the_v011_manifest(self):
        self.assertEqual(Allowlist().entries, {})
        # PD-NCBI-PLATFORM-VARIANCE 1.3: the approved Sakai output is the Gate A
        # expectation itself, read from the certification driver's manifest.
        expected = fast.capture_outputs.gate_a()
        self.assertEqual(len(expected), 43)
        self.assertEqual(
            expected[("blastn", "Sakai.MG1655.megablast")],
            "ac8177764f35b2fb92777f8b1d01ea646cc55ae21b05f1f486ee231348a34a45",
        )
        # The S02 baseline's outputs are the v0.1.1 expectations in all 43 cases.
        baseline = fast.read_tsv(fast.BASELINE)
        gate_a_rows = {key: row for key, row in baseline.items() if row["expected_source"] == "gate_a"}
        self.assertEqual(set(gate_a_rows), set(expected))
        self.assertEqual({key: row["output_sha256"] for key, row in gate_a_rows.items()}, expected)


class CheckRowsTests(unittest.TestCase):
    def test_baseline_and_frozen_checks(self):
        baseline = {("blastn", "known"): row("known", KNOWN), ("blastn", "same"): row("same", FROZEN)}
        failures, allowed = fast.check_rows([row("known", KNOWN), row("same", FROZEN)], baseline, allowlist(), {"blastn"}, True)
        self.assertEqual((failures, allowed), ([], ["blastn/known"]))
        failures, _ = fast.check_rows([row("known", FROZEN), row("same", FROZEN)], baseline, allowlist(), {"blastn"}, False)
        self.assertTrue(any("output_sha256 differs" in line for line in failures))
        self.assertTrue(any("remove it" in line for line in failures))
        failures, _ = fast.check_rows([row("same", FROZEN)], baseline, allowlist(), {"blastn"}, True)
        self.assertTrue(any("not executed" in line for line in failures))
        failures, _ = fast.check_rows([row("same", FROZEN)], baseline, allowlist(), {"blastn"}, False)
        self.assertEqual(failures, [])


class BlastnFixtureTests(unittest.TestCase):
    def test_manifest_matches_cases_and_files(self):
        for module in (fixtures, tblastx_fixtures):
            rows = module.read_manifest()
            self.assertEqual([(r["case_id"], r["argv"], r["losat_extra"], r.get("env") or "", r.get("oracle_env") or "")
                              for r in rows], module.CASES, module.__name__)
            for r in rows:
                data = (module.FIXTURES / f"{r['case_id']}.out").read_bytes()
                self.assertEqual(module.sha256(data), r["stdout_sha256"], r["case_id"])
                err = module.FIXTURES / f"{r['case_id']}.err"
                self.assertEqual(module.sha256(err.read_bytes()) if err.exists() else "", r["stderr_sha256"], r["case_id"])


class DbSearchFixtureTests(unittest.TestCase):
    def test_manifest_rows_match_frozen_files(self):
        rows = db_fixtures.read_manifest()[1]
        self.assertEqual({r["program"] for r in rows}, set(db_fixtures.PROGRAMS))
        # ci_fast runs the -db rows of every program the manifest has, at 1 and 4 threads.
        self.assertEqual(set(fast.DB_SEARCH_PROGRAMS), set(db_fixtures.PROGRAMS))
        self.assertTrue(set(fast.DB_SEARCH_PROGRAMS) <= set(fast.PROGRAMS))
        self.assertEqual(len({r["row_id"] for r in rows}), len(rows))
        by_id = {r["row_id"]: r for r in rows}
        for r in rows:
            if r["kind"] != "ncbi":
                self.assertFalse(db_fixtures.under_valgrind(r), r["row_id"])
                self.assertIsNone(db_fixtures.copies_cell(r), r["row_id"])
                # A rejection row whose stdout is a near input's NCBI run (`ncbi-near:`, the bytes LOSAT writes
                # before it refuses): the stdout cells of that `ncbi` row of its own run, or of its first N bytes
                # (`ncbi-near: ROW_ID first N bytes`), and no file of its own.
                near, cut = db_fixtures.near_row_id(r), db_fixtures.near_cut(r)
                if near:
                    source = by_id[near]
                    self.assertEqual((r["kind"], source["kind"], db_fixtures.near_row_id(source), source["outfmt"]),
                                     ("losat-rejection", "ncbi", None, r["outfmt"]), r["row_id"])
                    self.assertFalse((db_fixtures.FROZEN / f"{r['row_id']}.out").exists(), r["row_id"])
                    if cut is None:
                        self.assertEqual([source[k] for k in ("stdout_sha256", "stdout_bytes")],
                                         [r[k] for k in ("stdout_sha256", "stdout_bytes")], r["row_id"])
                    else:
                        self.assertLessEqual(cut, int(source["stdout_bytes"]), r["row_id"])
                        self.assertEqual(r["stdout_bytes"], str(cut), r["row_id"])
                # A rejection row that expects the first N bytes of NCBI's run of its own command
                # (`ncbi-first: N bytes`) has those N bytes as its own frozen stdout.
                own = db_fixtures.own_cut(r)
                if own is not None:
                    self.assertEqual((r["kind"], near, r["stdout_bytes"]), ("losat-rejection", None, str(own)),
                                     r["row_id"])
                if r["stdout_bytes"]:
                    self.assertTrue((db_fixtures.FROZEN / f"{db_fixtures.frozen_id(r)}.out").exists(), r["row_id"])
                    data = db_fixtures.expected_stdout(r)
                    self.assertEqual([db_fixtures.sha256(data), str(len(data))],
                                     [r["stdout_sha256"], r["stdout_bytes"]], r["row_id"])
                continue
            # A copy row's expectation is the cut of its whole NCBI output (`ncbi-copies:`).
            if db_fixtures.copies_cell(r) is not None:
                self.assertEqual(db_fixtures.copy_row_problems(r), [], r["row_id"])
                self.assertEqual(r["rc"], "0", r["row_id"])
                # Its stderr, where NCBI writes one, is a warning of the arguments that names no copy.
                err = db_fixtures.FROZEN / f"{r['row_id']}.err"
                self.assertFalse(err.exists() and db_fixtures.COPY_PREFIX.encode() in err.read_bytes(), r["row_id"])
            # An ncbi row takes all of its own run's bytes (`ncbi-first:` is a rejection row's).
            self.assertIsNone(db_fixtures.own_cut(r), r["row_id"])
            # A row that expects a near input's NCBI run (`ncbi-near:`) reads that row's files and has its columns.
            near = db_fixtures.near_row_id(r)
            # A row frozen from NCBI under valgrind (`ncbi-valgrind`) has frozen files of its own run.
            self.assertFalse(near and db_fixtures.under_valgrind(r), r["row_id"])
            if near:
                source = by_id[near]
                # An ncbi row takes all of its near row's bytes (rc and stderr too).
                self.assertIsNone(db_fixtures.near_cut(r), r["row_id"])
                self.assertEqual((source["kind"], db_fixtures.near_row_id(source), source["outfmt"]),
                                 ("ncbi", None, r["outfmt"]), r["row_id"])
                self.assertEqual([source[k] for k in ("rc", "stdout_sha256", "stdout_bytes", "stderr_sha256")],
                                 [r[k] for k in ("rc", "stdout_sha256", "stdout_bytes", "stderr_sha256")], r["row_id"])
                self.assertFalse((db_fixtures.FROZEN / f"{r['row_id']}.out").exists(), r["row_id"])
            frozen = db_fixtures.frozen_id(r)
            data = (db_fixtures.FROZEN / f"{frozen}.out").read_bytes()
            self.assertEqual(db_fixtures.sha256(data), r["stdout_sha256"], r["row_id"])
            self.assertEqual(str(len(data)), r["stdout_bytes"], r["row_id"])
            err = db_fixtures.FROZEN / f"{frozen}.err"
            self.assertEqual(db_fixtures.sha256(err.read_bytes()) if err.exists() else "", r["stderr_sha256"], r["row_id"])

    def test_near_marker_names_a_row_and_its_first_bytes(self):
        def near(covers):
            return db_fixtures.near_row({"row_id": "r", "covers": covers})
        self.assertIsNone(near("a crash | expect: frees memory"))
        self.assertEqual(near("a | ncbi-near: db.x.1"), ("db.x.1", None))
        self.assertEqual(near("a | ncbi-near: db.x.1 | expect: b"), ("db.x.1", None))
        self.assertEqual(near("a | ncbi-near: db.x.1 first 349 bytes | expect: b"), ("db.x.1", 349))
        self.assertEqual(near("ncbi-near: db.x.1 first 0 bytes"), ("db.x.1", 0))
        self.assertEqual((db_fixtures.near_row_id({"covers": "ncbi-near: db.x.1 first 9 bytes"}),
                          db_fixtures.near_cut({"covers": "ncbi-near: db.x.1 first 9 bytes"})), ("db.x.1", 9))
        self.assertEqual((db_fixtures.near_row_id({"covers": "x"}), db_fixtures.near_cut({"covers": "x"})),
                         (None, None))
        # Any other words in the marker's part stop the tool (a cut must not be read as prose, or prose as a cut).
        for bad in ("ncbi-near: db.x.1 first bytes", "ncbi-near: db.x.1 349", "ncbi-near: db.x.1 first 349 byte",
                    "ncbi-near: db.x.1 first -1 bytes", "ncbi-near: db.x.1 (the first 349 bytes)", "ncbi-near: | x"):
            with self.assertRaises(SystemExit, msg=bad):
                near(bad)

    def test_own_first_bytes_marker_is_a_part_of_covers(self):
        def own(covers):
            return db_fixtures.own_cut({"row_id": "r", "covers": covers})
        self.assertIsNone(own("staxid, expect: an uninitialised staxid"))
        self.assertEqual(own("staxid, expect: x | ncbi-first: 12 bytes"), 12)
        self.assertEqual(own("a | ncbi-first: 0 bytes | expect: b"), 0)
        # Only a part after a ` | `; any other words in that part stop the tool.
        self.assertIsNone(own("ncbi-first: 12 bytes"))
        for bad in ("a | ncbi-first: 12", "a | ncbi-first: twelve bytes", "a | ncbi-first: -1 bytes",
                    "a | ncbi-first: 12 byte"):
            with self.assertRaises(SystemExit, msg=bad):
                own(bad)

    def test_valgrind_marker_is_a_part_of_covers(self):
        self.assertTrue(db_fixtures.under_valgrind({"covers": "RG-J2A: windows past the buffers | ncbi-valgrind"}))
        self.assertTrue(db_fixtures.under_valgrind({"covers": "a | ncbi-valgrind | ncbi-near: db.x"}))
        self.assertFalse(db_fixtures.under_valgrind({"covers": "ncbi-valgrind"}))
        self.assertFalse(db_fixtures.under_valgrind({"covers": "frozen under ncbi-valgrind, see | the row"}))
        self.assertFalse(db_fixtures.under_valgrind({"covers": "ncbi-near: db.x | expect: valgrind"}))
        self.assertTrue(db_fixtures.valgrind_messages(b"Warning\n==1234== Invalid read of size 1\n"))
        self.assertFalse(db_fixtures.valgrind_messages(b"Warning: [blastx] Query_1: no hits == none\n"))

    def test_env_tokens_are_environment_not_arguments(self):
        row = {"extra_args": "BLASTDB=tests/seqdb/taxdb/good BATCH_SIZE=100 -evalue 1e-5 -word_size 3",
               "db": "d", "query": "q", "outfmt": "16"}
        self.assertEqual(db_fixtures.search_argv(row),
                         ["-db", "d", "-query", "q", "-evalue", "1e-5", "-word_size", "3", "-outfmt", "16"])
        env = db_fixtures.row_env(row, {"PATH": "/bin"})
        self.assertEqual((env["BLASTDB"], env["BATCH_SIZE"], env["PATH"]),
                         ("tests/seqdb/taxdb/good", "100", "/bin"))
        self.assertNotIn("CHUNK_SIZE", env)
        # a token of another name is an argument (and is not silently dropped)
        self.assertFalse(db_fixtures.is_env_token("-evalue"))
        self.assertFalse(db_fixtures.is_env_token("NCBI_CONFIG_PATH=/x"))
        self.assertFalse(db_fixtures.is_env_token("NCBI"))  # a name without `=` is a word of the argv
        # the variables of the `-html` formatter's configuration and of NCBI's registry files
        for token in ("FMTCFG=fmt.cfg", "NCBI=tests/db_search/html/ncbi", "GETURL_DEBUG=1", "NCBI_DONT_USE_NCBIRC=1"):
            self.assertTrue(db_fixtures.is_env_token(token), token)
        # the run directory: neither an environment variable nor an argument
        row = {"extra_args": "cwd=tests/x FMTCFG=fmt.cfg -html", "db": "d", "query": "q", "outfmt": "0"}
        self.assertFalse(db_fixtures.is_env_token("cwd=tests/x"))
        self.assertTrue(db_fixtures.is_cwd_token("cwd=tests/x"))
        self.assertNotIn("cwd=tests/x", db_fixtures.search_argv(row))
        self.assertNotIn("cwd", db_fixtures.row_env(row, {}))
        self.assertEqual(db_fixtures.row_env(row, {})["FMTCFG"], "fmt.cfg")  # a name in the run directory

    def test_a_cwd_row_names_its_paths_relative_to_its_run_directory(self):
        row = {"row_id": "db.html.x", "db": "tests/db_search/ids/prot tests/seqdb/db/swiss",
               "query": "tests/db_search/ids/ids_q.faa", "outfmt": "0",
               "extra_args": "cwd=tests/db_search/html/rc BLASTDB=tests/seqdb/taxdb/of18 NCBI=tests/db_search/html/n "
                             "FMTCFG=fmt.cfg -html"}
        self.assertEqual(db_fixtures.run_dir(row), "tests/db_search/html/rc")
        self.assertEqual(db_fixtures.row_cwd(row), db_fixtures.ENGINE / "tests/db_search/html/rc")
        self.assertEqual(db_fixtures.search_argv(row),
                         ["-db", "../../ids/prot ../../../seqdb/db/swiss", "-query", "../../ids/ids_q.faa", "-html",
                          "-outfmt", "0"])
        env = db_fixtures.row_env(row, {"PATH": "/bin"})
        self.assertEqual((env["BLASTDB"], env["NCBI"], env["FMTCFG"]), ("../../../seqdb/taxdb/of18", "../n", "fmt.cfg"))
        subject = {**row, "db": "subject:tests/fasta/outfmt0/s.faa"}
        self.assertEqual(db_fixtures.search_argv(subject)[:2], ["-subject", "../../../fasta/outfmt0/s.faa"])
        # the empty db cell stays `-db ""`; a row without the token keeps its paths and runs in LOSAT/
        self.assertEqual(db_fixtures.search_argv({**row, "db": ""})[:2], ["-db", ""])
        plain = {**row, "extra_args": "BLASTDB=tests/seqdb/taxdb/of18 NCBI=tests/n -html"}
        self.assertEqual(db_fixtures.search_argv(plain)[:4],
                         ["-db", "tests/db_search/ids/prot tests/seqdb/db/swiss", "-query", "tests/db_search/ids/ids_q.faa"])
        self.assertEqual(db_fixtures.row_env(plain, {})["NCBI"], "tests/n")
        self.assertEqual(db_fixtures.row_cwd(plain), db_fixtures.ENGINE)

    def test_a_cwd_row_takes_an_rps_gz_db_cell(self):
        row = {"row_id": "db.html.rps.x", "db": "rps-gz:rps_mini", "query": "tests/fasta/rps/xoj.faa", "outfmt": "0",
               "extra_args": "cwd=tests/db_search/html/rc/rep -html"}
        self.assertEqual(db_fixtures.run_dir(row), "tests/db_search/html/rc/rep")
        self.assertEqual(db_fixtures.search_argv(row), ["-db", "rps_mini", "-query", "../../../../fasta/rps/xoj.faa", "-html", "-outfmt", "0"])

    def test_a_cwd_row_takes_only_what_it_can_rebase(self):
        row = {"row_id": "db.html.x", "db": "tests/d", "query": "tests/q.faa", "outfmt": "0", "extra_args": "cwd=tests/r"}
        for change in ({"query": "pssm:tests/p.asn"}, {"query": "pssm-gz:tests/p.asn.gz"}, {"query": "stdin:empty"},
                       {"files": "xinclude"}, {"extra_args": "cwd=tests/r -seqidlist tests/l.txt"},
                       {"extra_args": "cwd=tests/r cwd=tests/s"}, {"extra_args": "cwd="}, {"extra_args": "cwd=/abs"}):
            with self.assertRaises(SystemExit, msg=str(change)):
                db_fixtures.search_argv({**row, **change})
        # every cwd= row of the manifest is one the runner takes, and its directory exists
        for r in db_fixtures.read_manifest()[1]:
            directory = db_fixtures.run_dir(r)
            if directory is not None:
                self.assertTrue((db_fixtures.ENGINE / directory).is_dir(), r["row_id"])

    def test_the_library_skips_the_rows_that_set_the_process_settings(self):
        # The formatter's settings of the process (the run directory's .ncbirc and FMTCFG file, $NCBI,
        # GETURL_DEBUG, NCBI_DONT_USE_NCBIRC) are not what a library call reads (LOSAT/src/api/mod.rs:78-86):
        # `check --library` skips such a row before it runs anything (the LOSAT path does not exist).
        row = {"row_id": "db.html.x", "program": "blastp", "kind": "ncbi", "db": "tests/d", "query": "tests/q.faa",
               "outfmt": "0", "covers": "", "rc": "0", "stdout_sha256": "", "stdout_bytes": "", "stderr_sha256": ""}
        missing = Path("/nonexistent/library_command_line")
        for extra in ("cwd=tests/db_search -html", "FMTCFG=fmt.cfg -html", "NCBI=tests/n -html", "GETURL_DEBUG=1 -html",
                      "NCBI_DONT_USE_NCBIRC=1 -html"):
            self.assertTrue(db_fixtures.process_settings({**row, "extra_args": extra}), extra)
            self.assertEqual(db_fixtures.check_row({**row, "extra_args": extra}, missing, 1, {}, library=True),
                             ("library-skip", ""), extra)
        for extra in ("-html", "BLASTDB=tests/seqdb/taxdb/of18 BATCH_SIZE=1 -html"):
            self.assertFalse(db_fixtures.process_settings({**row, "extra_args": extra}), extra)

    @unittest.skipUnless(os.name == "posix", "the stand-in programs are shell scripts")
    def test_freeze_and_check_run_a_cwd_row_in_its_directory(self):
        # NCBI's stand-in (freeze) and LOSAT's (check) print their directory, argv and variables: the same bytes,
        # so the row is `same`, and the bytes show the run directory and the paths relative to it.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp).resolve()  # the directory as the stand-in's `pwd` prints it
            engine, bin_dir, frozen = root / "LOSAT", root / "bin", root / "frozen"
            (engine / "tests/db_search/html/rc").mkdir(parents=True)
            bin_dir.mkdir()
            show = 'printf "%s|" "$(pwd)" "$@" "${FMTCFG-}" "${NCBI-}" "${BLASTDB-}"\n'
            (bin_dir / "blastp").write_text("#!/bin/sh\n" + show)
            (bin_dir / "losat").write_text("#!/bin/sh\nshift\n" + show)  # `losat blastp ...`
            for script in bin_dir.iterdir():
                script.chmod(0o755)
            row = {"row_id": "db.html.cwd", "program": "blastp", "kind": "ncbi", "db": "tests/db_search/ids/prot",
                   "query": "tests/db_search/ids/ids_q.faa", "outfmt": "0", "covers": "", "rc": "", "stdout_sha256": "",
                   "stdout_bytes": "", "stderr_sha256": "",
                   "extra_args": "cwd=tests/db_search/html/rc NCBI=tests/db_search/html FMTCFG=fmt.cfg -html"}
            with unittest.mock.patch.object(db_fixtures, "read_manifest", lambda: ({}, [row])), \
                    unittest.mock.patch.object(db_fixtures, "write_manifest", lambda *a: None), \
                    unittest.mock.patch.object(db_fixtures, "ENGINE", engine), \
                    unittest.mock.patch.object(db_fixtures, "FROZEN", frozen), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(db_fixtures.freeze(argparse.Namespace(ncbi_bin=str(bin_dir), program="", only="")), 0)
                self.assertEqual(db_fixtures.check_row(row, bin_dir / "losat", 1, {"PATH": os.environ["PATH"]}),
                                 ("same", ""))
                expected = (f"{engine / 'tests/db_search/html/rc'}|-db|../../ids/prot|-query|../../ids/ids_q.faa|-html|"
                            "-outfmt|0|-num_threads|1|fmt.cfg|..||")
                self.assertEqual((frozen / "db.html.cwd.out").read_text(), expected)
                # LOSAT/.ncbirc would change every row: both commands refuse to start
                (engine / ".ncbirc").write_text("[BLAST]\n")
                with self.assertRaises(SystemExit):
                    db_fixtures.freeze(argparse.Namespace(ncbi_bin=str(bin_dir), program="", only=""))
                with self.assertRaises(SystemExit):
                    db_fixtures.check(argparse.Namespace(losat=bin_dir / "losat"))

    def test_clean_env_drops_what_changes_a_result_and_seals_home(self):
        dirty = {"PATH": "/bin", "HOME": "/home/someone", "USERPROFILE": "C:/Users/someone",
                 "NCBI_CONFIG__BLAST__BLASTDB": "/x", "NCBI_CONFIG_PATH": "/y", "NCBI": "/z",
                 "BLASTDB": "/d", "BLAST_USE_NCBI_ENV": "yes", "LOSAT_X_SOMETHING": "1",
                 "BATCH_SIZE": "1", "CHUNK_SIZE": "2", "OLD_FSC": "1", "BL2SEQ_LEGACY": "1", "LANG": "C",
                 "ADAPTIVE_CBS": "1", "CTOOLKIT_COMPATIBLE": "", "FMTCFG": "fmt.cfg", "GETURL_DEBUG": "1",
                 "NCBI_DONT_USE_NCBIRC": "1"}
        clean = db_fixtures.clean_env(dirty, "/empty/home")
        self.assertEqual(clean, {"PATH": "/bin", "HOME": "/empty/home", "USERPROFILE": "/empty/home",
                                 "LANG": "C"})
        # nothing of the caller's environment is changed
        self.assertEqual(dirty["HOME"], "/home/someone")

    def test_ncbi_runs_have_the_usage_report_off(self):
        # NCBI BLAST+ 2.17.0 sends a run's parameters to NCBI at exit unless BLAST_USAGE_REPORT is false
        # (blast_usage_report.cpp:195-220): every NCBI run of the harness has it, and so has every LOSAT run
        # (`test_check_runs_losat_in_the_environment_of_the_freeze`).
        clean = db_fixtures.clean_env({"PATH": "/bin", "BLAST_USAGE_REPORT": "true"}, "/empty/home")
        self.assertNotIn("BLAST_USAGE_REPORT", clean)
        self.assertEqual(db_fixtures.ncbi_env(clean)["BLAST_USAGE_REPORT"], "false")
        self.assertEqual(db_fixtures.ncbi_env(clean), {**clean, "BLAST_USAGE_REPORT": "false"})
        # a row's tokens cannot turn it back on
        self.assertFalse(db_fixtures.is_env_token("BLAST_USAGE_REPORT=true"))

    @unittest.skipUnless(os.name == "posix", "the stand-in LOSAT is a shell script")
    def test_check_runs_losat_in_the_environment_of_the_freeze(self):
        # The variable decides what the programs read (without it, a second read of a malformed `.ncbirc` for
        # the usage report posts its Critical a second time: db.html.bp.di.f0.rc_bad.use, cfg.ncbibad.D.h/.p),
        # so `check` runs LOSAT, for every row, in the environment the row was frozen in (`ncbi_env`). The
        # stand-in prints the variable; the frozen file says `false`.
        rows = db_fixtures.read_manifest()[1]
        row = dict(next(r for r in rows if r["kind"] == "ncbi" and not xinclude.wants_files(r)
                        and not xinclude.wants_pssm(r) and db_fixtures.rps_cell(r["db"]) is None
                        and not db_fixtures.copies_cell(r) and r["outfmt"] == "0"))
        with tempfile.TemporaryDirectory() as tmp:
            frozen, losat = Path(tmp) / "frozen", Path(tmp) / "losat"
            frozen.mkdir()
            (frozen / f"{row['row_id']}.out").write_bytes(b"false")
            row.update(rc="0", stdout_sha256=db_fixtures.sha256(b"false"), stdout_bytes="5", stderr_sha256="")
            losat.write_text('#!/bin/sh\nprintf "%s" "${BLAST_USAGE_REPORT-unset}"\n')
            losat.chmod(0o755)
            out = io.StringIO()
            with unittest.mock.patch.object(db_fixtures, "read_manifest", lambda: ({}, [row])), \
                    unittest.mock.patch.object(db_fixtures, "FROZEN", frozen), \
                    unittest.mock.patch.dict(os.environ, {"BLAST_USAGE_REPORT": "true"}), \
                    contextlib.redirect_stdout(out):
                args = argparse.Namespace(losat=losat, threads=1, program="", only="", outfmts="", out=None,
                                          library=False)
                self.assertEqual(db_fixtures.check(args), 0, out.getvalue())
            self.assertIn(f"{row['row_id']}\tsame\t", out.getvalue())

    @unittest.skipUnless(os.name == "posix", "the stand-in NCBI programs are shell scripts")
    def test_freeze_runs_ncbi_with_the_usage_report_off(self):
        # `freeze` with stand-in NCBI programs that print the variable: a plain row, an RPS row
        # (`freeze_rps`) and a rejection row whose NCBI stdout is frozen; the manifest is not written.
        rows = db_fixtures.read_manifest()[1]
        def first(test):
            return dict(next(r for r in rows if test(r)))
        plain = first(lambda r: r["kind"] == "ncbi" and not xinclude.wants_files(r) and not xinclude.wants_pssm(r)
                      and db_fixtures.rps_cell(r["db"]) is None)
        rps = first(lambda r: r["kind"] == "ncbi" and db_fixtures.rps_cell(r["db"]) is not None)
        rejection = first(lambda r: r["kind"] == "losat-rejection" and r["stdout_bytes"]
                          and not xinclude.wants_pssm(r) and not db_fixtures.near_row_id(r))
        with tempfile.TemporaryDirectory() as tmp:
            bin_dir, frozen = Path(tmp) / "bin", Path(tmp) / "frozen"
            bin_dir.mkdir()
            for program in db_fixtures.PROGRAMS:
                script = bin_dir / program
                script.write_text('#!/bin/sh\nprintf "%s" "${BLAST_USAGE_REPORT-unset}"\n')
                script.chmod(0o755)
            written = []
            with unittest.mock.patch.object(db_fixtures, "read_manifest", lambda: ({}, [plain, rps, rejection])), \
                    unittest.mock.patch.object(db_fixtures, "write_manifest", lambda *a: written.append(a)), \
                    unittest.mock.patch.object(db_fixtures, "FROZEN", frozen), \
                    unittest.mock.patch.dict(os.environ, {"BLAST_USAGE_REPORT": "true"}), \
                    contextlib.redirect_stdout(io.StringIO()):
                args = argparse.Namespace(ncbi_bin=str(bin_dir), program="", only="")
                self.assertEqual(db_fixtures.freeze(args), 0)
            for row in (plain, rps, rejection):
                self.assertEqual((frozen / f"{row['row_id']}.out").read_bytes(), b"false", row["row_id"])
            self.assertEqual(len(written), 1)

    @unittest.skipUnless(os.name == "posix", "the stand-in NCBI program is a shell script")
    def test_freeze_gives_a_rejection_row_the_stdout_cells_of_its_near_row(self):
        # A rejection row with `ncbi-near:` takes the stdout cells of its near row and keeps its own rc; NCBI runs
        # only for the near row (the stand-in counts its runs), which may come after the rejection row in the
        # manifest, and no frozen file is written for the rejection row.
        rows = db_fixtures.read_manifest()[1]
        rejection = dict(next(r for r in rows if r["kind"] == "losat-rejection" and db_fixtures.near_row_id(r)
                              and db_fixtures.near_cut(r) is None))
        near = dict(next(r for r in rows if r["row_id"] == db_fixtures.near_row_id(rejection)))
        with tempfile.TemporaryDirectory() as tmp:
            bin_dir, frozen, runs = Path(tmp) / "bin", Path(tmp) / "frozen", Path(tmp) / "runs"
            bin_dir.mkdir()
            script = bin_dir / near["program"]
            script.write_text(f'#!/bin/sh\necho run >> "{runs}"\nprintf "the near run\\n"\nexit 3\n')
            script.chmod(0o755)
            written = []
            with unittest.mock.patch.object(db_fixtures, "read_manifest", lambda: ({}, [rejection, near])), \
                    unittest.mock.patch.object(db_fixtures, "write_manifest", lambda *a: written.append(a)), \
                    unittest.mock.patch.object(db_fixtures, "FROZEN", frozen), \
                    contextlib.redirect_stdout(io.StringIO()):
                args = argparse.Namespace(ncbi_bin=str(bin_dir), program="", only="", valgrind="")
                self.assertEqual(db_fixtures.freeze(args), 0)
            self.assertEqual(runs.read_text(), "run\n")
            self.assertFalse((frozen / f"{rejection['row_id']}.out").exists())
            self.assertEqual((frozen / f"{near['row_id']}.out").read_bytes(), b"the near run\n")
            self.assertEqual([rejection[k] for k in ("stdout_sha256", "stdout_bytes", "rc")],
                             [near["stdout_sha256"], "13", "1"])
            self.assertEqual(near["rc"], "3")
            self.assertEqual(len(written), 1)

    @staticmethod
    def stand_in(bin_dir: Path, program: str, body: str) -> None:
        bin_dir.mkdir(exist_ok=True)
        script = bin_dir / program
        script.write_text("#!/bin/sh\n" + body)
        script.chmod(0o755)

    def freeze_rows(self, rows: list[dict[str, str]], bin_dir: Path, frozen: Path) -> list[tuple]:
        """`freeze` of `rows` alone with the NCBI programs of `bin_dir`, the frozen files in `frozen`: the
        write_manifest calls (none when the freeze stops)."""
        written = []
        with unittest.mock.patch.object(db_fixtures, "read_manifest", lambda: ({}, rows)), \
                unittest.mock.patch.object(db_fixtures, "write_manifest", lambda *a: written.append(a)), \
                unittest.mock.patch.object(db_fixtures, "FROZEN", frozen), \
                contextlib.redirect_stdout(io.StringIO()):
            db_fixtures.freeze(argparse.Namespace(ncbi_bin=str(bin_dir), program="", only="", valgrind=""))
        return written

    @unittest.skipUnless(os.name == "posix", "the stand-in NCBI program is a shell script")
    def test_freeze_gives_a_cut_near_row_the_cells_of_the_first_bytes_of_its_near_run(self):
        # `ncbi-near: ROW_ID first N bytes` (the crash rows of db.x2g.bx.of, whose stdout is NCBI's prolog): NCBI
        # runs only for the near row; the rejection row's cells are the sha256 and size of the near run's first N
        # bytes, its rc stays LOSAT's, and a near run shorter than N stops the freeze.
        rows = db_fixtures.read_manifest()[1]
        rejection = dict(next(r for r in rows if r["kind"] == "losat-rejection" and db_fixtures.near_cut(r)))
        near = dict(next(r for r in rows if r["row_id"] == db_fixtures.near_row_id(rejection)))
        cut = db_fixtures.near_cut(rejection)
        with tempfile.TemporaryDirectory() as tmp:
            bin_dir, frozen, runs = Path(tmp) / "bin", Path(tmp) / "frozen", Path(tmp) / "runs"
            (Path(tmp) / "near.out").write_bytes(b"P" * cut + b"\n\nQuery= q\n")
            self.stand_in(bin_dir, near["program"], f'echo run >> "{runs}"\ncat "{tmp}/near.out"\nexit 255\n')
            self.assertEqual(len(self.freeze_rows([rejection, near], bin_dir, frozen)), 1)
            self.assertEqual(runs.read_text(), "run\n")
            self.assertFalse((frozen / f"{rejection['row_id']}.out").exists())
            self.assertEqual([rejection[k] for k in ("stdout_sha256", "stdout_bytes", "rc")],
                             [db_fixtures.sha256(b"P" * cut), str(cut), "1"])
            self.assertEqual((near["rc"], near["stdout_bytes"]), ("255", str(cut + 11)))
            (Path(tmp) / "near.out").write_bytes(b"P" * (cut - 1))
            with self.assertRaisesRegex(SystemExit, f"{re.escape(rejection['row_id'])}: ncbi-near .* fewer than"):
                self.freeze_rows([rejection, near], bin_dir, frozen)
            # An ncbi row takes all of its near row's bytes: a cut there stops the freeze.
            twin = dict(rejection, kind="ncbi")
            with self.assertRaisesRegex(SystemExit, "takes all of its near row's bytes"):
                self.freeze_rows([twin, near], bin_dir, frozen)

    @unittest.skipUnless(os.name == "posix", "the stand-in NCBI program is a shell script")
    def test_freeze_stops_where_ncbi_ends_a_run_by_a_signal(self):
        # A crash cuts NCBI's report where its heap and the environment decide: no row takes such a stdout. The
        # freeze stops before it writes the row's files or the manifest, for a rejection row with stdout cells and
        # no near row as for an ncbi row; a rejection row with a near row is not run at all.
        rows = db_fixtures.read_manifest()[1]
        rejection = dict(next(r for r in rows if r["kind"] == "losat-rejection" and r["stdout_bytes"]
                              and not xinclude.wants_pssm(r) and not db_fixtures.near_row_id(r)))
        plain = dict(next(r for r in rows if r["kind"] == "ncbi" and not xinclude.wants_files(r)
                          and not xinclude.wants_pssm(r) and db_fixtures.rps_cell(r["db"]) is None
                          and not db_fixtures.under_valgrind(r) and not db_fixtures.near_row_id(r)))
        with tempfile.TemporaryDirectory() as tmp:
            bin_dir, frozen = Path(tmp) / "bin", Path(tmp) / "frozen"
            for row in (rejection, plain):
                self.stand_in(bin_dir, row["program"], 'printf "BLASTP 2.17.0+\\n"\nkill -SEGV $$\n')
                with self.assertRaisesRegex(SystemExit, f"{re.escape(row['row_id'])}: NCBI ended its run by signal 11"):
                    self.freeze_rows([dict(row)], bin_dir, frozen)
                self.assertFalse((frozen / f"{row['row_id']}.out").exists(), row["row_id"])
                self.assertFalse((frozen / f"{row['row_id']}.err").exists(), row["row_id"])
            # The same rows with an exit status are frozen as NCBI wrote them.
            for row in (rejection, plain):
                self.stand_in(bin_dir, row["program"], 'printf "BLASTP 2.17.0+\\n"\nexit 3\n')
                row = dict(row)
                self.assertEqual(len(self.freeze_rows([row], bin_dir, frozen)), 1)
                self.assertEqual((frozen / f"{row['row_id']}.out").read_bytes(), b"BLASTP 2.17.0+\n")
                self.assertEqual(row["stdout_bytes"], "15")
            # A crash of the near row stops the freeze too; the rejection row that names it is never run.
            cut_row = dict(next(r for r in rows if r["kind"] == "losat-rejection" and db_fixtures.near_cut(r)))
            near = dict(next(r for r in rows if r["row_id"] == db_fixtures.near_row_id(cut_row)))
            self.stand_in(bin_dir, near["program"], 'kill -SEGV $$\n')
            with self.assertRaisesRegex(SystemExit, f"{re.escape(near['row_id'])}: NCBI ended its run by signal 11"):
                self.freeze_rows([cut_row, near], bin_dir, frozen)

    @unittest.skipUnless(os.name == "posix", "the stand-in LOSAT program is a shell script")
    def test_check_compares_a_cut_near_row_with_the_first_bytes_of_its_near_run(self):
        # `check` of a rejection row with `ncbi-near: ROW_ID first N bytes` expects exactly the first N bytes of
        # ROW_ID's frozen stdout (the prolog of the db.x2g.bx.of rows) in LOSAT's words (AGENTS.md rule 8), and
        # LOSAT's rc and expected stderr text.
        rows = db_fixtures.read_manifest()[1]
        row = dict(next(r for r in rows if r["kind"] == "losat-rejection" and db_fixtures.near_cut(r)))
        whole = (db_fixtures.FROZEN / f"{db_fixtures.near_row_id(row)}.out").read_bytes()
        cut = db_fixtures.near_cut(row)
        self.assertGreater(len(whole), cut)
        self.assertEqual(db_fixtures.expected_stdout(row), whole[:cut])
        with tempfile.TemporaryDirectory() as tmp:
            losat, out = Path(tmp) / "losat", Path(tmp) / "out"
            losat.write_text(f'#!/bin/sh\ncat "{out}"\nprintf "LOSAT engine error: %s\\n" '
                             f'{shlex.quote(db_fixtures.expected_text(row))} >&2\nexit {row["rc"]}\n')
            losat.chmod(0o755)
            env = db_fixtures.clean_env(dict(os.environ), tmp)
            words = db_fixtures.losat_wording.losat_output_bytes
            self.assertNotEqual(words(whole[:cut]), whole[:cut])  # the prolog states the version
            for stdout, status in ((words(whole[:cut]), "rejection"), (whole[:cut], "differs"),
                                   (words(whole), "differs"), (words(whole[:cut])[:-1], "differs"), (b"", "differs")):
                out.write_bytes(stdout)
                self.assertEqual(db_fixtures.check_row(row, losat, 1, env)[0], status, len(stdout))


def pairwise_report(query: str, hits: list[str]) -> str:
    """A query's report of NCBI's outfmt 0 (the shape of BLASTP 2.17.0's), with one alignment per hit."""
    table = "".join(f"{hit} a protein                                                   19.6    3.2  \n" for hit in hits)
    body = ("Sequences producing significant alignments:                          (Bits)  Value\n\n"
            f"{table}\n\n" + "".join(f">{hit} a protein\nLength=40\n\n Score = 19.6 bits (39),  Expect = 3.2\n\n"
                                       "Query  1  MKV  3\n           MKV\nSbjct  1  MKV  3\n\n\n" for hit in hits)
            if hits else "\n\n***** No hits found *****\n\n")
    return (f"Query= {query}\n\nLength=82\n" + body + "\n"
            "Lambda      K        H        a         alpha\n   0.318    0.136    0.374    0.792     4.96 \n\n"
            "Effective search space used: 1998016\n\n\n")


PAIRWISE_HEAD = ("BLASTP 2.17.0+\n\n\nReference: ...\n\n\n\nDatabase: PajaWSV proteins\n"
                 "           86 sequences; 66,738 total letters\n\n\n\n")
PAIRWISE_TAIL = ("  Database: PajaWSV proteins\n    Posted date:  Oct 13, 2026  1:44 PM\n"
                 "  Number of letters in database: 66,738\n  Number of sequences in database:  86\n\n\n\n"
                 "Matrix: BLOSUM62\nGap Penalties: Existence: 11, Extension: 1\n")


def tabular_block(query: str, rows: list[str]) -> str:
    """A query's block of NCBI's outfmt 7."""
    fields = "# Fields: query acc.ver, subject acc.ver, evalue\n" if rows else ""
    return (f"# BLASTP 2.17.0+\n# Query: {query}\n# Database: tests/db_search/paja\n{fields}# {len(rows)} hits found\n"
            + "".join(f"{query.split()[0]}\t{row}\n" for row in rows))


def xml_iteration(number: int, query: str, hits: list[str]) -> str:
    """A query's `<Iteration>` element of NCBI's outfmt 5."""
    return (f"<Iteration>\n  <Iteration_iter-num>{number}</Iteration_iter-num>\n"
            f"  <Iteration_query-ID>Query_{number}</Iteration_query-ID>\n  <Iteration_query-def>{query}</Iteration_query-def>\n"
            "  <Iteration_hits>\n" + "".join(f"<Hit>\n  <Hit_id>{hit}</Hit_id>\n</Hit>\n" for hit in hits)
            + "  </Iteration_hits>\n</Iteration>\n")


XML_HEAD = "<?xml version=\"1.0\"?>\n<BlastOutput>\n  <BlastOutput_query-def>A a protein</BlastOutput_query-def>\n<BlastOutput_iterations>\n"
XML_TAIL = "</BlastOutput_iterations>\n</BlastOutput>\n\n"


class CopyRowTests(unittest.TestCase):
    """Copy rows (`ncbi-copies:`): the cut of NCBI's output on a query file with copies, NCBI's query batches and
    the checks of the oracle file."""

    def test_query_batches_take_queries_until_the_batch_size_is_reached(self):
        self.assertEqual(db_fixtures.query_batches([250, 310, 200, 60], 260), [[0, 1], [2, 3]])
        self.assertEqual(db_fixtures.query_batches([250, 250, 200, 200], 300), [[0, 1], [2, 3]])
        self.assertEqual(db_fixtures.query_batches([100, 100], 10000), [[0, 1]])
        self.assertEqual(db_fixtures.query_batches([100, 100], 1), [[0], [1]])
        self.assertEqual(db_fixtures.query_batches([100, 100, 5], 200), [[0, 1], [2]])
        self.assertEqual(db_fixtures.query_batches([], 10), [])
        self.assertEqual(db_fixtures.batch_size({}), 10000)
        self.assertEqual(db_fixtures.batch_size({"BATCH_SIZE": "260"}), 260)

    def test_outfmt_6_drops_the_copies_rows(self):
        whole = b"A\ts1\t1e-5\nA\ts2\t3.2\nB\ts1\t0.1\nCOPY\ts1\t0.2\nCOPY2\ts3\t9\n"
        self.assertEqual(db_fixtures.cut_copies(whole, "6", 2), b"A\ts1\t1e-5\nA\ts2\t3.2\nB\ts1\t0.1\n")
        self.assertEqual(db_fixtures.cut_copies(b"", "6", 1), b"")
        self.assertEqual(db_fixtures.subjects_per_query(whole, "6"), [2, 1, 1, 1])

    def test_outfmt_7_drops_the_copies_blocks_and_counts_the_queries_kept(self):
        a, b = tabular_block("A a protein", ["s1\t1e-5", "s2\t3.2"]), tabular_block("B", [])
        copy, copy2 = tabular_block("COPY", ["s1\t0.2"]), tabular_block("COPY2", [])
        whole = (a + copy + b + copy2 + "# BLAST processed 4 queries\n").encode()
        self.assertEqual(db_fixtures.cut_copies(whole, "7", 2), (a + b + "# BLAST processed 2 queries\n").encode())
        trailing = (a + b + copy + "# BLAST processed 3 queries\n").encode()
        self.assertEqual(db_fixtures.cut_copies(trailing, "7", 1), (a + b + "# BLAST processed 2 queries\n").encode())
        self.assertEqual(db_fixtures.subjects_per_query(whole, "7"), [2, 1])
        for wrong in (whole.replace(b"4 queries", b"3 queries"), whole[:-len(b"# BLAST processed 4 queries\n")]):
            with self.assertRaises(ValueError):
                db_fixtures.cut_copies(wrong, "7", 2)
        with self.assertRaises(ValueError):
            db_fixtures.cut_copies(whole, "7", 1)

    def test_outfmt_0_drops_the_copies_reports_up_to_the_next_query_or_the_footer(self):
        a, b = pairwise_report("A a protein", ["s1", "s2"]), pairwise_report("B", [])
        copy, copy2 = pairwise_report("COPY", ["s1"]), pairwise_report("COPY2", [])
        whole = (PAIRWISE_HEAD + a + copy + b + copy2 + PAIRWISE_TAIL).encode()
        self.assertEqual(db_fixtures.cut_copies(whole, "0", 2), (PAIRWISE_HEAD + a + b + PAIRWISE_TAIL).encode())
        self.assertEqual(db_fixtures.subjects_per_query(whole, "0"), [2, 1, 0, 0])
        with self.assertRaises(ValueError):
            db_fixtures.cut_copies(whole, "0", 1)
        with self.assertRaises(ValueError):
            db_fixtures.cut_copies((PAIRWISE_HEAD + a + copy).encode(), "0", 1)

    def test_outfmt_5_drops_a_last_copy_only(self):
        a, b = xml_iteration(1, "A a protein", ["s1", "s2"]), xml_iteration(2, "B", [])
        whole = (XML_HEAD + a + b + xml_iteration(3, "COPY", ["s1"]) + XML_TAIL).encode()
        self.assertEqual(db_fixtures.cut_copies(whole, "5", 1), (XML_HEAD + a + b + XML_TAIL).encode())
        self.assertEqual(db_fixtures.subjects_per_query(whole, "5"), [2, 0, 1])
        middle = (XML_HEAD + a + xml_iteration(2, "COPY", []) + xml_iteration(3, "B", []) + XML_TAIL).encode()
        two = (XML_HEAD + a + xml_iteration(2, "COPY", []) + xml_iteration(3, "COPY2", []) + XML_TAIL).encode()
        for wrong, copies in ((middle, 1), (two, 2), (whole, 2)):
            with self.assertRaises(ValueError):
                db_fixtures.cut_copies(wrong, "5", copies)

    def test_other_formats_have_no_rule(self):
        for outfmt in ("11", "10", "6 qseqid sseqid", "17"):
            with self.assertRaises(ValueError):
                db_fixtures.cut_copies(b"A\ts1\n", outfmt, 1)

    def oracle_row(self, tmp: Path, query: list[str], oracle: list[str], **cells: str) -> dict[str, str]:
        (tmp / "q.faa").write_text("".join(query))
        (tmp / "o.faa").write_text("".join(oracle))
        words = cells.pop("words", "")
        row = {"row_id": "db.x.c.6", "program": "blastp", "kind": "ncbi", "query": str(tmp / "q.faa"),
               "extra_args": "-comp_based_stats 2u", "outfmt": "6",
               "covers": f"a copy row | ncbi-copies: {tmp / 'o.faa'}{words}"}
        return {**row, **cells}

    def test_the_oracle_file_is_checked_against_the_query_file_and_the_batches(self):
        a, b = ">A a protein\nMKVLA\nGG\n", ">B\nWWWW\n"
        copy_a, copy_b, copy2_b = ">COPY\nMKVLAGG\n", ">COPY\nWWWW\n", ">COPY2\nWWWW\n"
        with tempfile.TemporaryDirectory() as name:
            tmp = Path(name)
            self.assertEqual(db_fixtures.copy_oracle_problems(self.oracle_row(tmp, [a, b], [a, b, copy_b])), [])
            self.assertEqual(db_fixtures.copy_count(self.oracle_row(tmp, [a, b], [a, b, copy_b])), 1)
            for outfmt in ("0", "5", "7"):
                self.assertEqual(db_fixtures.copy_oracle_problems(
                    self.oracle_row(tmp, [a, b], [a, b, copy_b], outfmt=outfmt)), [], outfmt)
            # Two batches of LOSAT's run (BATCH_SIZE=1): each query is followed by its copy, and the oracle's
            # BATCH_SIZE (8: 7 + 7 and 4 + 4) makes each copy end a batch.
            two = self.oracle_row(tmp, [a, b], [a, copy_a, b, copy2_b], extra_args="BATCH_SIZE=1 -comp_based_stats 2u",
                                  words=" BATCH_SIZE=8")
            self.assertEqual(db_fixtures.copies_cell(two)[1], {"BATCH_SIZE": "8"})
            self.assertEqual(db_fixtures.copy_oracle_problems(two), [])
            # Without the oracle's own BATCH_SIZE each query is alone in its batch, and the copies too.
            self.assertTrue(db_fixtures.copy_oracle_problems({**two, "covers": f"ncbi-copies: {tmp / 'o.faa'}"}))
            # A hit-list limit is allowed: NCBI's redo ends no match of a query that another query follows.
            for limit, outfmt in (("-max_target_seqs 1", "6"), ("-num_descriptions 1 -num_alignments 1", "0")):
                self.assertEqual(db_fixtures.copy_oracle_problems(self.oracle_row(
                    tmp, [a, b], [a, b, copy_b], extra_args=f"{limit} -comp_based_stats 2u", outfmt=outfmt)), [], limit)
            # A batch-ending query without its copy, a copy of other residues or before its query, a copy row of
            # outfmt 5 whose copy is not last, another program, format or kind.
            wrong = [
                self.oracle_row(tmp, [a, b], [a, b]),
                self.oracle_row(tmp, [a, b], [a, b, copy_a]),
                self.oracle_row(tmp, [a, b], [copy_b, a, b]),
                self.oracle_row(tmp, [a, b], [a, copy_a, b, copy2_b]),
                self.oracle_row(tmp, [a, b], [a, copy_a, b, copy2_b], extra_args="BATCH_SIZE=1", words=" BATCH_SIZE=8",
                                outfmt="5"),
                self.oracle_row(tmp, [a, b], [a, b, copy_b], program="psiblast"),
                self.oracle_row(tmp, [a, b], [a, b, copy_b], outfmt="6 qseqid sseqid"),
                self.oracle_row(tmp, [a, b], [a, b, copy_b], outfmt="11"),
                self.oracle_row(tmp, [a, ">B other\nWWWW\n"], [a, b, copy_b]),
                self.oracle_row(tmp, [a, copy_b], [a, copy_b, ">COPY2\nWWWW\n"]),
            ]
            for row in wrong:
                self.assertTrue(db_fixtures.copy_oracle_problems(row), (row["query"], row["outfmt"], row["program"]))
            with self.assertRaises(SystemExit):
                db_fixtures.copies_cell({**two, "covers": f"ncbi-copies: {tmp / 'o.faa'} -evalue 1"})
            # Only the variables that place the queries in batches and chunks may differ in the oracle run.
            with self.assertRaises(SystemExit):
                db_fixtures.copies_cell({**two, "covers": f"ncbi-copies: {tmp / 'o.faa'} ADAPTIVE_CBS=1"})
            self.assertEqual(db_fixtures.copies_cell({**two, "covers": "ncbi-copies: o.faa OVERLAP_CHUNK_SIZE=50"}),
                             ("o.faa", {"OVERLAP_CHUNK_SIZE": "50"}))

    def test_query_chunks_follow_ncbis_split(self):
        self.assertEqual(db_fixtures.query_chunks(19799, {}), 1)
        self.assertEqual(db_fixtures.query_chunks(19800, {}), 2)
        self.assertEqual(db_fixtures.query_chunks(249, {"CHUNK_SIZE": "200"}), 2)
        self.assertEqual(db_fixtures.query_chunks(249, {"CHUNK_SIZE": "180"}), 3)
        self.assertEqual(db_fixtures.query_chunks(249, {"CHUNK_SIZE": "350"}), 1)
        self.assertEqual(db_fixtures.query_chunks(249, {"CHUNK_SIZE": "150", "OVERLAP_CHUNK_SIZE": "50"}), 2)
        # A chunk size not above the overlap, or a blank value (the default), splits nothing.
        self.assertEqual(db_fixtures.query_chunks(249, {"CHUNK_SIZE": "100"}), 1)
        self.assertEqual(db_fixtures.query_chunks(19800, {"CHUNK_SIZE": " "}), 2)

    def test_a_copy_stays_in_its_batch_and_each_batch_is_searched_in_one_chunk(self):
        a, b, copy_b = ">A a protein\nMKVLA\nGG\n", ">B\nWWWW\n", ">COPY\nWWWW\n"
        with tempfile.TemporaryDirectory() as name:
            tmp = Path(name)
            # [A, B] is one batch of LOSAT's run (11 residues); with BATCH_SIZE=11 for the oracle the batch ends at B
            # and the copy leaves it.
            self.assertEqual(db_fixtures.copy_oracle_problems(self.oracle_row(tmp, [a, b], [a, b, copy_b])), [])
            leaves = db_fixtures.copy_oracle_problems(self.oracle_row(tmp, [a, b], [a, b, copy_b], words=" BATCH_SIZE=11"))
            self.assertTrue(any("does not end with its one copy" in problem for problem in leaves), leaves)
            # CHUNK_SIZE=106 (6 residues a chunk beside the overlap of 100): LOSAT's batch of 11 residues is one chunk,
            # the oracle's of 15 two; the oracle's own CHUNK_SIZE=200 keeps it in one.
            row = self.oracle_row(tmp, [a, b], [a, b, copy_b], extra_args="CHUNK_SIZE=106 -comp_based_stats 2u")
            self.assertEqual(db_fixtures.copy_oracle_problems(row),
                             ["the oracle batch of queries [0, 1, 2] (15 residues) is searched in 2 chunks; a copy "
                              "row's batches are searched in one"])
            self.assertEqual(db_fixtures.copy_oracle_problems({**row, "covers": f"ncbi-copies: {tmp / 'o.faa'} "
                                                                                 "CHUNK_SIZE=200"}), [])
            # CHUNK_SIZE=105: LOSAT's own batch of 11 residues is split too.
            both = db_fixtures.copy_oracle_problems({**row, "extra_args": "CHUNK_SIZE=105 -comp_based_stats 2u"})
            self.assertEqual([problem.split(" (")[0] for problem in both],
                             ["the oracle batch of queries [0, 1, 2]", "the LOSAT batch of queries [0, 1]"])

    @unittest.skipUnless(os.name == "posix", "the stand-in NCBI program is a shell script")
    def test_freeze_keeps_the_whole_output_and_writes_its_cut(self):
        # A stand-in NCBI that prints the query file it was given as outfmt 6 rows: the oracle's, never the row's.
        a, b, copy_b = ">A\nMKV\n", ">B\nWW\n", ">COPY\nWW\n"
        with tempfile.TemporaryDirectory() as name:
            tmp = Path(name)
            row = {**self.oracle_row(tmp, [a, b], [a, b, copy_b]), "rc": "", "stdout_sha256": "",
                   "stdout_bytes": "", "stderr_sha256": "", "db": "d"}
            bin_dir, frozen = tmp / "bin", tmp / "frozen"
            bin_dir.mkdir()
            frozen.mkdir()
            script = bin_dir / "blastp"
            script.write_text('#!/bin/sh\nwhile [ "$1" != -query ]; do shift; done\n'
                              'grep ">" "$2" | sed "s/>//; s/$/\ts1\t${BATCH_SIZE-none}/"\n')
            script.chmod(0o755)
            written = []
            with unittest.mock.patch.object(db_fixtures, "read_manifest", lambda: ({}, [row])), \
                    unittest.mock.patch.object(db_fixtures, "write_manifest", lambda *a: written.append(a)), \
                    unittest.mock.patch.object(db_fixtures, "FROZEN", frozen), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(db_fixtures.freeze(argparse.Namespace(ncbi_bin=str(bin_dir), program="", only="")), 0)
                self.assertEqual(db_fixtures.copy_row_problems(row), [])
                (frozen / "db.x.c.6.out").write_bytes(b"A\ts1\tnone\n")
                self.assertEqual(db_fixtures.copy_row_problems(row), ["the cut of db.x.c.6.copies.out is not db.x.c.6.out"])
            self.assertEqual((frozen / "db.x.c.6.copies.out").read_bytes(), b"A\ts1\tnone\nB\ts1\tnone\nCOPY\ts1\tnone\n")
            self.assertEqual(row["stdout_bytes"], str(len(b"A\ts1\tnone\nB\ts1\tnone\n")))
            self.assertEqual((row["rc"], row["stderr_sha256"]), ("0", ""))
            # The copy's own report changes from run to run (NCBI reads past its array): a whole output whose cut is
            # the new run's is kept; one whose cut differs is replaced.
            (frozen / "db.x.c.6.copies.out").write_bytes(b"A\ts1\tnone\nB\ts1\tnone\nCOPY\ts9\tother\n")
            for whole, kept in ((b"A\ts1\tnone\nB\ts1\tnone\nCOPY\ts9\tother\n", True), (b"A\ts2\tnone\n", False)):
                (frozen / "db.x.c.6.copies.out").write_bytes(whole)
                with unittest.mock.patch.object(db_fixtures, "read_manifest", lambda: ({}, [row])), \
                        unittest.mock.patch.object(db_fixtures, "write_manifest", lambda *a: None), \
                        unittest.mock.patch.object(db_fixtures, "FROZEN", frozen), \
                        contextlib.redirect_stdout(io.StringIO()):
                    db_fixtures.freeze(argparse.Namespace(ncbi_bin=str(bin_dir), program="", only=""))
                self.assertEqual((frozen / "db.x.c.6.copies.out").read_bytes(),
                                 whole if kept else b"A\ts1\tnone\nB\ts1\tnone\nCOPY\ts1\tnone\n")
            # A stderr that names no copy (a warning of the arguments) is the row's own and is frozen; NCBI's oracle
            # run must exit 0 with a stderr that names no copy.
            script.write_text(script.read_text() + "echo 'Warning: [blastp] Examining 5 or more matches is "
                                                   "recommended' >&2\n")
            with unittest.mock.patch.object(db_fixtures, "read_manifest", lambda: ({}, [row])), \
                    unittest.mock.patch.object(db_fixtures, "write_manifest", lambda *a: None), \
                    unittest.mock.patch.object(db_fixtures, "FROZEN", frozen), \
                    contextlib.redirect_stdout(io.StringIO()):
                db_fixtures.freeze(argparse.Namespace(ncbi_bin=str(bin_dir), program="", only=""))
                self.assertEqual(db_fixtures.copy_row_problems(row), [])
            warning = b"Warning: [blastp] Examining 5 or more matches is recommended\n"
            self.assertEqual((frozen / "db.x.c.6.err").read_bytes(), warning)
            self.assertEqual(row["stderr_sha256"], db_fixtures.sha256(warning))
            for stand_in in ("#!/bin/sh\necho 'Warning: [blastp] Query_3 COPY: an invalid query' >&2\n",
                             "#!/bin/sh\nexit 3\n"):
                script.write_text(stand_in)
                with unittest.mock.patch.object(db_fixtures, "read_manifest", lambda: ({}, [row])), \
                        unittest.mock.patch.object(db_fixtures, "write_manifest", lambda *a: None), \
                        unittest.mock.patch.object(db_fixtures, "FROZEN", frozen), \
                        contextlib.redirect_stdout(io.StringIO()), self.assertRaises(SystemExit):
                    db_fixtures.freeze(argparse.Namespace(ncbi_bin=str(bin_dir), program="", only=""))


def fake_search(directory: Path, master: bytes, sides: list[bytes], stdout: bytes = b"", extra: str = "") -> list[str]:
    """A stand-in for a program that writes the files of outfmt 13/14: the command to run it."""
    script = directory / "fake.py"
    lines = [
        "import sys",
        "out = sys.argv[sys.argv.index('-out') + 1]",
        f"open(out, 'wb').write({master!r})",
        f"for n, side in enumerate({sides!r}, 1):",
        "    open(f'{out}_{n}.xml', 'wb').write(side)",
    ]
    if extra:
        lines.append(f"open(out + {extra!r}, 'wb').write(b'x')")
    lines.append(f"sys.stdout.buffer.write({stdout!r})")
    script.write_text("\n".join(lines) + "\n")
    return [sys.executable, str(script), "blastp", "-outfmt", "14"]


class XincludeFixtureTests(unittest.TestCase):
    """The `files` columns of the manifests (outfmt 13 and 14: a master file and side files)."""

    def manifest_lines(self, path: Path) -> tuple[list[str], list[list[str]]]:
        lines = [line for line in path.read_text().split("\n") if line and not line.startswith("#")]
        return lines[0].split("\t"), [line.split("\t") for line in lines[1:]]

    def test_rows_without_the_columns_keep_their_text(self):
        # the outfmt 0 manifest: a freeze writes every row as it was read (no trailing empty cells added)
        header, rows = self.manifest_lines(run_oracle.MANIFEST)
        self.assertEqual(header[run_oracle.BASE_COLUMNS:], xinclude.FILES_COLUMNS)
        for cells in rows:
            row = dict(zip(header, cells))
            self.assertEqual(run_oracle.row_cells(header, row), cells, cells[0])
        # the -db manifest likewise (each program's files: its file, then its row family files), and its
        # header is the old one only while no row has the columns
        parsed = db_fixtures.read_manifest()[1]
        listed = []
        for program in db_fixtures.PROGRAMS:
            for path in db_fixtures.program_files(program):
                header, rows = self.manifest_lines(path)
                self.assertEqual(header, db_fixtures.header_cells(parsed), path.name)
                listed += rows
        self.assertEqual(len(listed), len(parsed))
        for cells, row in zip(listed, parsed):
            self.assertEqual(db_fixtures.row_cells(row), cells, cells[0])
        old = {column: "" for column in db_fixtures.COLUMNS}
        self.assertEqual(db_fixtures.header_cells([old]), db_fixtures.COLUMNS)
        self.assertEqual(len(db_fixtures.row_cells(old)), len(db_fixtures.COLUMNS))
        self.assertEqual(db_fixtures.header_cells([{**old, "files": "xinclude"}]),
                         db_fixtures.COLUMNS + xinclude.FILES_COLUMNS)

    def test_files_rows_match_the_frozen_master_and_side_files(self):
        found = pssm_rows = 0
        for manifest_rows, frozen, key in (
                (run_oracle.read_manifest()[2], run_oracle.ENGINE / "tests/fixtures/outfmt0", "fixture_id"),
                (db_fixtures.read_manifest()[1], db_fixtures.FROZEN, "row_id")):
            for row in manifest_rows:
                if xinclude.wants_pssm(row):
                    # psiblast's PSSM rows: `side_sha256` names the PSSM files (empty: none); those of
                    # the `pssm` kind are all frozen, the others when 64 KB or less
                    pssm_rows += 1
                    self.assertTrue(xinclude.frozen_pssm_agrees(row, frozen), row[key])
                    if xinclude.files_kind(row) == xinclude.PSSM:
                        files, _ = xinclude.pssm_files(frozen, row[key])
                        self.assertEqual(xinclude.pssm_cell(files, row[key]), row["side_sha256"], row[key])
                    continue
                if not xinclude.wants_files(row):
                    self.assertEqual(row.get("side_sha256", ""), "", row[key])
                    continue
                found += 1
                row_id, outfmt = row[key], row["outfmt"]
                master = (frozen / f"{row_id}.out").read_bytes()
                self.assertEqual(xinclude.sha256(master), row["stdout_sha256"], row_id)
                self.assertEqual(str(len(master)), row["stdout_bytes"], row_id)
                hashes = row["side_sha256"].split(",") if row["side_sha256"] else []
                names = [xinclude.side_name(row_id, n, outfmt) for n in range(1, len(hashes) + 1)]
                self.assertEqual(sorted(path.name for path in frozen.glob(f"{row_id}.out_*")), sorted(names), row_id)
                self.assertEqual([xinclude.sha256((frozen / name).read_bytes()) for name in names], hashes, row_id)
                # a run that NCBI stops (rc other than 0) writes the side files of the batches before
                # and no master, whose list comes after the last query's file: a batch's prefetch
                # stops (line DBID mixdb2), an exception (BLASTX's empty rows, line RG X2g: rc 255)
                if row.get("rc", "0") != "0":
                    self.assertEqual(master, b"", row_id)
                    continue
                # the master lists the side files by those names (13: no end of line after the last)
                listed = re.findall(rb'href="([^"]*)"' if outfmt == "14" else rb'"File": "([^"]*)"', master)
                self.assertEqual([name.decode() for name in listed], names, row_id)
                self.assertEqual(master.endswith(b"\n"), outfmt == "14", row_id)
        self.assertGreaterEqual(found, 9)
        self.assertGreaterEqual(pssm_rows, 2)

    def test_pssm_files_are_listed_asn_first_and_frozen_by_kind(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            big = b"x" * (xinclude.FREEZE_LIMIT + 1)
            for name, data in (("r.pssm", b"ascii"), ("r.asn.2", b"two"), ("r.asn.1", big), ("r.asn", b"asn"),
                               ("other", b"")):
                (directory / name).write_bytes(data)
            files, stray = xinclude.pssm_files(directory, "r")
            self.assertEqual([name for name, _ in files], ["r.asn", "r.asn.1", "r.asn.2", "r.pssm"])
            self.assertEqual(stray, ["other"])
            result = xinclude.PssmResult(None, files)
            cell = xinclude.pssm_cell(files, "r", xinclude.PSSM_BOTH)
            self.assertTrue(cell.startswith(f".asn={xinclude.sha256(b'asn')}:3,"), cell)
            row = {"row_id": "r", "files": xinclude.PSSM_BOTH, "side_sha256": cell}
            frozen = directory / "frozen"
            frozen.mkdir()
            xinclude.freeze_pssm(result, frozen, "r", xinclude.PSSM_BOTH)
            self.assertEqual(sorted(path.name for path in frozen.iterdir()), ["r.asn", "r.asn.2", "r.pssm"])
            self.assertTrue(xinclude.frozen_pssm_agrees(row, frozen))
            # Every listed file of 64 KB or less must be frozen, and the cell must give the sizes.
            (frozen / "r.asn.2").unlink()
            self.assertFalse(xinclude.frozen_pssm_agrees(row, frozen))
            (frozen / "r.asn.2").write_bytes(b"two")
            self.assertFalse(xinclude.frozen_pssm_agrees({**row, "side_sha256": xinclude.pssm_cell(files, "r")}, frozen))
            first_difference = lambda want, got: "differs"  # noqa: E731
            self.assertEqual(xinclude.pssm_problems(result, row, frozen, first_difference), [])
            changed = xinclude.PssmResult(None, [(name, b"y" if name == "r.asn.1" else data) for name, data in files])
            problems = xinclude.pssm_problems(changed, row, frozen, first_difference)
            self.assertEqual(len(problems), 1)
            self.assertIn("r.asn.1 differs", problems[0])
            self.assertIn("not frozen", problems[0])
            # The `pssm` kind freezes every file, so the big one must be frozen too.
            self.assertFalse(xinclude.frozen_pssm_agrees({**row, "files": xinclude.PSSM}, frozen))

    def test_pssm_pre_files_exist_before_the_run(self):
        # `{pre:<suffix>}` writes `<tmp>/<row_id><suffix>` before the run and is not passed on.
        check = "import sys; from pathlib import Path; assert not any('{pre' in a for a in sys.argv); " \
                "assert Path(sys.argv[2]).with_suffix('.asn.1').read_bytes() == b'previous file\\n'"
        result = xinclude.run_with_pssm([sys.executable, "-c", check, "{pre:.asn.1}", "{pre:.pssm}"], "r",
                                        xinclude.PSSM_ASN)
        self.assertEqual(result.run.returncode, 0, result.run.stderr)
        self.assertEqual(result.files, [("r.asn.1", xinclude.PRE_BYTES), ("r.pssm", xinclude.PRE_BYTES)])
        self.assertEqual(result.stray, [])

    def test_side_files_are_collected_in_the_order_of_their_numbers(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            for n in (10, 2, 1):
                (directory / f"r.out_{n}.xml").write_bytes(str(n).encode())
            (directory / "r.out").write_bytes(b"master")
            master, sides, stray = xinclude.collect(directory, "r", "14")
            self.assertEqual((master, stray), (b"master", ["side file numbers [1, 2, 10] are not 1..3"]))
            self.assertEqual(sides, [b"1", b"2", b"10"])
            (directory / "r.out_10.xml").unlink()
            (directory / "other.txt").write_bytes(b"")
            (directory / "r.out_3.json").write_bytes(b"")
            _, sides, stray = xinclude.collect(directory, "r", "14")
            self.assertEqual((sides, stray), ([b"1", b"2"], ["other.txt", "r.out_3.json"]))

    def test_a_run_is_compared_with_the_frozen_master_side_files_and_standard_output(self):
        def first_difference(expected: bytes, actual: bytes) -> str:
            return f"expected {expected!r}, got {actual!r}"

        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            frozen = directory / "frozen"
            frozen.mkdir()
            row = {"row_id": "r", "outfmt": "14", "stdout_sha256": xinclude.sha256(b"M"), "stdout_bytes": "1",
                   "side_sha256": ",".join(xinclude.sha256(side) for side in (b"A", b"B"))}
            written = xinclude.run_with_files(fake_search(directory, b"M", [b"A", b"B"]), "r", "14")
            self.assertEqual(xinclude.observed(written)["side_sha256"], row["side_sha256"])
            xinclude.freeze(written, frozen, "r", "14")
            self.assertEqual(sorted(path.name for path in frozen.iterdir()), ["r.out", "r.out_1.xml", "r.out_2.xml"])
            self.assertEqual(xinclude.problems(written, row, frozen, first_difference), [])
            # a different master, side file, number of side files, a standard output, a stray file
            for master, sides, stdout, extra, expect in (
                    (b"N", [b"A", b"B"], b"", "", "master"),
                    (b"M", [b"A", b"C"], b"", "", "side file 2"),
                    (b"M", [b"A"], b"", "", "1 side files, expected 2"),
                    (b"M", [b"A", b"B"], b"x", "", "stdout is not empty"),
                    (b"M", [b"A", b"B"], b"", "_stray", "unexpected files")):
                other = xinclude.run_with_files(fake_search(directory, master, sides, stdout, extra), "r", "14")
                found_problems = xinclude.problems(other, row, frozen, first_difference)
                self.assertTrue(any(expect in line for line in found_problems), (expect, found_problems))
            # a freeze refuses standard output
            noisy = xinclude.run_with_files(fake_search(directory, b"M", [b"A"], b"x"), "r", "14")
            with self.assertRaises(SystemExit):
                xinclude.freeze(noisy, frozen, "r", "14")

    def test_a_comparison_counts_a_stray_file_as_a_difference_and_writes_nothing(self):
        # run_oracle.py in compare mode (no --freeze) with a stand-in for the NCBI program: a stray file
        # or standard output is a difference of that fixture (the comparison goes on), --out stays empty
        # and the manifest is not written.
        def compare(stray: bool, stdout: bool) -> tuple[int, str, list[str]]:
            with tempfile.TemporaryDirectory() as tmp:
                bin_dir, out = Path(tmp) / "bin", Path(tmp) / "out"
                bin_dir.mkdir()
                out.mkdir()
                program = bin_dir / "blastp"
                program.write_text("\n".join([
                    f"#!{sys.executable}", "import sys",
                    "if sys.argv[1] == '-version':", "    print('blastp: 2.17.0+'); sys.exit(0)",
                    "out = sys.argv[sys.argv.index('-out') + 1]",
                    "open(out, 'wb').write(b'M')", "open(out + '_1.xml', 'wb').write(b'A')",
                    "open(out + '_stray', 'wb').write(b'x')" if stray else "",
                    "sys.stdout.write('x')" if stdout else ""]) + "\n")
                program.chmod(0o755)
                good = {"fixture_id": "r", "program": "blastp", "outfmt": "14", "files": "xinclude",
                        "stdout_sha256": xinclude.sha256(b"M"), "stdout_bytes": "1", "stderr_sha256": "",
                        "side_sha256": xinclude.sha256(b"A")}
                rows = [good, {**good, "fixture_id": "s"}]
                saved = (sys.argv, run_oracle.read_manifest, run_oracle.search_argv)
                sys.argv = ["run_oracle.py", "--bin-dir", str(bin_dir), "--out", str(out)]
                run_oracle.read_manifest = lambda: ([], list(good), rows)
                run_oracle.search_argv = lambda row: ["-outfmt", "14"]
                captured = io.StringIO()
                try:
                    with contextlib.redirect_stdout(captured):
                        status = run_oracle.main()
                finally:
                    sys.argv, run_oracle.read_manifest, run_oracle.search_argv = saved
                return status, captured.getvalue(), [path.name for path in out.iterdir()]

        status, printed, written = compare(stray=False, stdout=False)
        self.assertEqual((status, written), (0, []), printed)
        for stray, stdout in ((True, False), (False, True)):
            status, printed, written = compare(stray, stdout)
            self.assertEqual((status, written), (1, []), printed)
            self.assertIn("differing=['r', 's']", printed)

    def test_the_files_column_and_the_format_are_validated(self):
        self.assertFalse(xinclude.wants_files({}))
        self.assertFalse(xinclude.wants_files({"files": ""}))
        self.assertTrue(xinclude.wants_files({"files": "xinclude"}))
        with self.assertRaises(SystemExit):
            xinclude.wants_files({"files": "yes"})
        self.assertEqual((xinclude.extension("13"), xinclude.extension("14")), ("json", "xml"))
        with self.assertRaises(SystemExit):
            xinclude.extension("16")
        self.assertEqual(xinclude.side_name("id", 3, "13"), "id.out_3.json")


class SelectionTests(unittest.TestCase):
    def test_paths_select_programs(self):
        everything = set(fast.PROGRAMS)
        self.assertEqual(fast.select_programs(["docs/README.md", "web/app/src/main.ts"]), set())
        self.assertEqual(fast.select_programs(["LOSAT/src/core/blast_util.rs"]), everything)
        self.assertEqual(fast.select_programs(["LOSAT/Cargo.lock"]), everything)
        self.assertEqual(fast.select_programs([".github/workflows/ci.yml"]), everything)
        # A program's file selects the programs whose code imports its directory, directly or through
        # another program, not every program: blastn imports no other program, so a tblastn file does
        # not select it.
        tblastn = fast.select_programs(["LOSAT/src/algorithm/tblastn/args.rs"])
        self.assertEqual(tblastn, fast.program_dependents()["tblastn"])
        self.assertIn("tblastn", tblastn)
        self.assertNotIn("blastn", tblastn)
        blastn = fast.select_programs(["LOSAT/src/algorithm/blastn/hsp.rs"])
        self.assertIn("blastn", blastn)
        # TBLASTX uses BLASTN's hit list (blastn/hsp.rs HitList) and blank-input check.
        self.assertIn("tblastx", blastn)
        for program in fast.PROGRAMS:
            self.assertIn(program, fast.program_dependents()[program])


class OutputFormatSelectionTests(unittest.TestCase):
    """`--outfmts` of the two fixture checkers, and the rows that check_wasm_threading.py runs with it: blastx's
    BLAST XML, XML2 and JSON rows through each WASI command (`blastx_prescan_checkers`)."""

    CELLS = ["5", "1", "16", "13", "6 qseqid sseqid", "15", "14", "0"]

    def test_a_row_is_kept_by_the_first_word_of_its_format_whole(self):
        rows = [{"row_id": f"r{index}", "program": "blastx", "outfmt": cell} for index, cell in enumerate(self.CELLS)]

        def kept(outfmts: str, program: str = "blastx") -> list[str]:
            args = argparse.Namespace(only="", program=program, outfmts=outfmts)
            return [row["outfmt"] for row in rows if db_fixtures.selected(row, args)]

        self.assertEqual(kept("5,13,14,15,16"), ["5", "16", "13", "15", "14"])
        self.assertEqual(kept("1"), ["1"])
        self.assertEqual(kept("6"), ["6 qseqid sseqid"])
        self.assertEqual(kept(""), self.CELLS)
        self.assertEqual(kept("5", "blastp"), [])
        # `freeze` has no --outfmts: every row.
        self.assertEqual(sum(db_fixtures.selected(row, argparse.Namespace(only="", program="")) for row in rows), 8)
        # check_losat.py: the same rule, an empty cell being outfmt 0.
        rows = [{"fixture_id": f"x{index}", "program": "blastx", "outfmt": cell}
                for index, cell in enumerate(["", *self.CELLS])]
        rows.append({"fixture_id": "y", "program": "blastp", "outfmt": "5"})

        def subject(outfmts: str, programs: str = "blastx", only: str = "") -> list[str]:
            return [row["fixture_id"] for row in check_losat.selected_rows(rows, programs, only, outfmts)]

        self.assertEqual(subject("0"), ["x0", "x8"])
        self.assertEqual(subject("5,16"), ["x1", "x3"])
        self.assertEqual(subject("5", "blastx,blastp"), ["x1", "y"])
        self.assertEqual(subject("", "blastx", "x1"), ["x1"])
        self.assertEqual(len(subject("")), 9)

    def test_the_wasi_section_runs_every_blastx_xml_and_json_row_of_both_sets(self):
        import check_wasm_threading as wasm

        formats = wasm.BLASTX_PRESCAN_FORMATS.split(",")
        self.assertEqual(formats, ["5", "13", "14", "15", "16"])
        subject = check_losat.selected_rows(run_oracle.read_manifest()[2], "blastx", "", wasm.BLASTX_PRESCAN_FORMATS)
        args = argparse.Namespace(only="", program="blastx", outfmts=wasm.BLASTX_PRESCAN_FORMATS)
        db = [row for row in db_fixtures.read_manifest()[1] if db_fixtures.selected(row, args)]
        # Every format, in the -subject and the -db rows; every blastx row of these formats and no other.
        for rows, every, number in [
                (subject, [row for row in run_oracle.read_manifest()[2] if row["program"] == "blastx"],
                 lambda row: run_oracle.outfmt(row).split()[0]),
                (db, [row for row in db_fixtures.read_manifest()[1] if row["program"] == "blastx"],
                 db_fixtures.format_number)]:
            self.assertEqual(sorted({number(row) for row in rows}, key=int), formats)
            self.assertEqual(rows, [row for row in every if number(row) in formats])
        # The section's command lines are the checkers' own: with a prefix that no row has, each selects nothing
        # and passes (no LOSAT runs).
        for name, checker, _ in wasm.blastx_prescan_checkers(Path("no-such-losat")):
            run = subprocess.run([sys.executable, *checker, "--only", "no-such-row"], capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, (name, run.stderr))
            self.assertRegex(run.stdout.splitlines()[-1], r"^# (fixtures=0 differing=\[\]|rows=0 differing=0 )", name)


class StartupTraceNoticeTests(unittest.TestCase):
    """`check_startup_trace_notice` against stub executables standing in for `losat`."""

    def check(self, body: str) -> list[str]:
        with tempfile.TemporaryDirectory() as tmp:
            program = Path(tmp) / "losat"
            program.write_text(f"#!{sys.executable}\nimport os, sys\n{body}\n")
            program.chmod(0o755)
            return fast.check_startup_trace_notice(program)

    def test_a_default_build_passes(self):
        self.assertEqual(self.check(
            "print('losat 1.0')\n"
            "if os.environ.get('LOSAT_STARTUP_TRACE'):\n"
            f"    sys.stderr.write({fast.STARTUP_TRACE_NOTICE!r})"), [])

    def test_a_build_that_reads_the_variable_fails(self):
        failures = self.check(
            "print('losat 1.0')\n"
            "if os.environ.get('LOSAT_STARTUP_TRACE'):\n"
            "    sys.stderr.write('[startup] enter main\\n')")
        self.assertEqual(len(failures), 1)
        self.assertIn("default build", failures[0])

    def test_stderr_without_the_variable_and_a_failing_exit_fail(self):
        failures = self.check("print('losat 1.0')\nsys.stderr.write('warning\\n')\nsys.exit(1)")
        self.assertTrue(any("exited" in failure for failure in failures))
        self.assertTrue(any("stderr without" in failure for failure in failures))


class RangeRegressionTaskTests(unittest.TestCase):
    """The frozen NCBI range errors (range_regression_fixtures.py) run whenever a classic program is selected."""

    def names(self, programs: set[str]) -> list[str]:
        with tempfile.TemporaryDirectory() as tmp:
            return [task.name for task in fast.fixture_tasks(Path("/x/losat"), programs, Path(tmp), None)]

    def test_each_classic_program_selects_the_range_fixtures_once(self):
        for program in fast.RANGE_PROGRAMS:
            self.assertEqual(self.names({program}).count("range regression"), 1, program)
        self.assertEqual(self.names(set(fast.RANGE_PROGRAMS)).count("range regression"), 1)

    def test_other_programs_do_not(self):
        self.assertNotIn("range regression", self.names({"psiblast", "rpsblast"}))

    def test_the_fixture_programs_are_the_manifest_programs(self):
        import range_regression_fixtures as ranges
        self.assertEqual({row["argv"].split(" ", 1)[0] for row in ranges.read_manifest()}, set(fast.RANGE_PROGRAMS))


class RunTasksTests(unittest.TestCase):
    """`run_tasks`: the fixture scripts run `jobs` at a time; their failure lines come in the order of the tasks."""

    @staticmethod
    def task(name: str, code: str, log: Path | None = None) -> fast.Task:
        return fast.Task(name, [sys.executable, "-c", code],
                         lambda result: [] if result.returncode == 0 else [f"{name} exit {result.returncode}"], log=log)

    def test_failures_in_task_order_logs_and_seconds(self):
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "b.tsv"
            # `a` ends last but its failure comes first.
            tasks = [self.task("a", "import time, sys; time.sleep(0.5); sys.exit(1)"),
                     self.task("b", "import sys; print('b-out'); sys.stderr.write('b-err')", log),
                     self.task("c", "import sys; sys.exit(2)")]
            failures, seconds = fast.run_tasks(tasks, 3)
            self.assertEqual(failures, ["a exit 1", "c exit 2"])
            self.assertEqual(log.read_text(), "b-out\nb-err")
            self.assertEqual(list(seconds), ["a", "b", "c"])

    def test_the_tasks_run_at_once(self):
        # Each task writes its flag, then waits (at most 60 s) for the other's: both pass only when they run
        # at the same time.
        with tempfile.TemporaryDirectory() as tmp:
            def waiter(mine: str, other: str) -> str:
                return (f"import pathlib, sys, time; d = pathlib.Path({tmp!r}); (d / {mine!r}).touch(); "
                        f"end = time.time() + 60\n"
                        f"while not (d / {other!r}).exists() and time.time() < end: time.sleep(0.05)\n"
                        f"sys.exit(0 if (d / {other!r}).exists() else 1)")
            failures, _ = fast.run_tasks([self.task("x", waiter("x", "y")), self.task("y", waiter("y", "x"))], 2)
            self.assertEqual(failures, [])

    def test_the_default_jobs_are_the_cpus_this_process_may_use(self):
        expected = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count() or 4
        self.assertEqual(fast.default_jobs(), expected)


STAND_IN_NCBI = '#!/bin/sh\nprintf "%s" "${BLAST_USAGE_REPORT-unset}"\n'


def stand_in_bin(directory: Path, programs) -> Path:
    """Stand-in NCBI programs that print the value of BLAST_USAGE_REPORT (`unset` without it)."""
    bin_dir = directory / "bin"
    bin_dir.mkdir()
    for program in programs:
        script = bin_dir / program
        script.write_text(STAND_IN_NCBI)
        script.chmod(0o755)
    return bin_dir


@unittest.skipUnless(os.name == "posix", "the stand-in NCBI programs are shell scripts")
class NcbiUsageReportOffTests(unittest.TestCase):
    """NCBI BLAST+ 2.17.0 sends a run's parameters to NCBI at exit unless BLAST_USAGE_REPORT is false
    (blast_usage_report.cpp:195-220): the harnesses' NCBI runs have it (db_search_fixtures.ncbi_env), the
    runs of LOSAT of these harnesses keep their environment (db_search_fixtures.py's `check` gives LOSAT the
    freeze's)."""

    def caller_env(self, home: str) -> dict[str, str]:
        # the freeze commands refuse to run with these set or with a ~/.ncbirc
        env = {k: v for k, v in os.environ.items() if k not in fasta_input.REPORT_ENV + range_fixtures.REPORT_ENV}
        return env | {"HOME": home, "BLAST_USAGE_REPORT": "true"}

    def test_fasta_input_freeze(self):
        one = next(r for r in fasta_input.build_rows() if r.kind == "ncbi" and fasta_input.network_safe(r) is None)
        with tempfile.TemporaryDirectory() as tmp:
            bin_dir = stand_in_bin(Path(tmp), {"blastn", one.program})
            written = []
            args = argparse.Namespace(ncbi_bin=str(bin_dir), out=str(Path(tmp) / "out"), no_netns=True, rows=None,
                                      match=None, tier=None, verify=False, jobs=1)
            with unittest.mock.patch.object(fasta_input, "build_rows", lambda: [one]), \
                    unittest.mock.patch.object(fasta_input, "read_manifest", lambda: []), \
                    unittest.mock.patch.object(fasta_input, "write_manifest", lambda *a: written.append(a)), \
                    unittest.mock.patch.dict(os.environ, self.caller_env(tmp), clear=True), \
                    contextlib.redirect_stdout(io.StringIO()):
                fasta_input.command_freeze(args)
                # LOSAT's run keeps the caller's environment
                self.assertEqual(fasta_input.execute([str(bin_dir / one.program)], one, False)[1], b"true")
            self.assertEqual((Path(tmp) / "out" / f"{one.row_id}.out").read_bytes(), b"false")
            self.assertEqual((Path(tmp) / "out" / f"{one.row_id}.merged").read_bytes(), b"false")
            self.assertEqual(written[0][1], "false")  # the version line: `blastn -version`'s output

    def test_range_regression_freeze(self):
        case = range_fixtures.CASES[0]
        program = case[1].split(" ", 1)[0]
        with tempfile.TemporaryDirectory() as tmp:
            bin_dir = stand_in_bin(Path(tmp), {"blastn", program})
            fixtures_dir = Path(tmp) / "fixtures"
            with unittest.mock.patch.object(range_fixtures, "CASES", [case]), \
                    unittest.mock.patch.object(range_fixtures, "FIXTURES", fixtures_dir), \
                    unittest.mock.patch.object(range_fixtures, "MANIFEST", fixtures_dir / "manifest.tsv"), \
                    unittest.mock.patch.dict(os.environ, self.caller_env(tmp), clear=True), \
                    contextlib.redirect_stdout(io.StringIO()):
                range_fixtures.command_freeze(argparse.Namespace(bin_dir=str(bin_dir)))
                # LOSAT's runs (`check_one`) go through run_case without `env_of`
                self.assertEqual(fixtures.run_case(str(bin_dir / program), [], "x", "", "").stdout, b"true")
            self.assertEqual((fixtures_dir / f"{case[0]}.out").read_bytes(), b"false")
            self.assertIn("of false (", (fixtures_dir / "manifest.tsv").read_text())



class WasmDbTotalRowsTests(unittest.TestCase):
    """The -db rows that check_wasm_threading.py runs through each WASI command for a database total of 2^32
    letters or more (`DB_TOTAL_ROWS`): NCBI's rows of one database whose alias file declares such a `LENGTH`
    (below 2^63: from there NCBI's `Int8` is negative, and the statistics read the length of OID 0)."""

    def test_the_rows_are_every_ncbi_row_of_one_database_of_2_to_the_32_letters_or_more(self):
        import check_wasm_threading as watcher

        def alias_length(row: dict[str, str]) -> int:
            names = row["db"].split()
            if len(names) != 1 or ":" in names[0]:
                return 0
            for suffix in (".nal", ".pal"):
                alias = db_fixtures.ENGINE / (names[0] + suffix)
                if alias.is_file():
                    for line in alias.read_text().splitlines():
                        key, _, value = line.partition(" ")
                        if key == "LENGTH":
                            return int(value)
            return 0

        rows = db_fixtures.read_manifest()[1]
        expected = [row["row_id"] for row in rows if row["kind"] == "ncbi" and 2**32 <= alias_length(row) < 2**63]
        self.assertEqual(sorted(watcher.DB_TOTAL_ROWS), sorted(expected))
        self.assertEqual(len(set(watcher.DB_TOTAL_ROWS)), 25)
        # The watcher selects each row by `--only` (a prefix of the row id): it names that row alone.
        for row_id in watcher.DB_TOTAL_ROWS:
            self.assertEqual([row["row_id"] for row in rows if row["row_id"].startswith(row_id)], [row_id])



class WasmLocalTimeTests(unittest.TestCase):
    """The local times that check_wasm_threading.py puts aside in the outputs of makeblastdb, blastdbcmd -info and
    a search report of its local-time-and-path block (`local_times`), and its test that each is the local time of
    its run (`is_local_time_of`)."""

    def test_each_form_is_read_and_put_aside(self):
        import calendar
        import check_wasm_threading as watcher

        log = (b"\n\nBuilding a new DB, current time: 11/18/2026 21:27:29\nNew DB name:   /w/db\n"
               b"Adding sequences from FASTA; added 3 sequences in 0.00400209 seconds.\n")
        self.assertEqual(watcher.local_times(log), (
            b"\n\nBuilding a new DB, current time: <time>\nNew DB name:   /w/db\n"
            b"Adding sequences from FASTA; added 3 sequences in <seconds> seconds.\n",
            [(calendar.timegm((2026, 11, 18, 21, 27, 29)), 1)]))
        # makeblastdb 2.17.0's form of a short build time.
        self.assertEqual(watcher.local_times(b"added 1 sequences in 5.5e-05 seconds.\n")[0],
                         b"added 1 sequences in <seconds> seconds.\n")
        info = b"Database: t\n\t3 sequences; 2,700 total bases\n\nDate: Nov 8, 2026  9:07 PM\tLongest sequence: 900 bases\n"
        self.assertEqual(watcher.local_times(info), (
            b"Database: t\n\t3 sequences; 2,700 total bases\n\nDate: <time>\tLongest sequence: 900 bases\n",
            [(calendar.timegm((2026, 11, 8, 21, 7, 0)), 60)]))
        report = b"  Database: t\n    Posted date:  Jan 1, 2027  12:05 AM\n  Number of letters in database: 2,700\n"
        self.assertEqual(watcher.local_times(report), (
            b"  Database: t\n    Posted date:  <time>\n  Number of letters in database: 2,700\n",
            [(calendar.timegm((2027, 1, 1, 0, 5, 0)), 60)]))
        self.assertEqual(watcher.local_times(b"Date: Dec 31, 2026  12:59 PM\t")[1],
                         [(calendar.timegm((2026, 12, 31, 12, 59, 0)), 60)])
        # Text that only looks like a time stays.
        other = b"Query= Date: Nov 8, 2026  9:07 PM\nbuilt 11/18/2026 21:27:29\n"
        self.assertEqual(watcher.local_times(other), (other, []))

    def test_a_time_is_the_local_time_of_its_run_in_the_zone(self):
        import check_wasm_threading as watcher

        started, ended = 1_795_000_000.4, 1_795_000_002.9
        jst = 32400
        # The log, to the second: any second of the run, one second of slack at each end.
        for printed in [1_795_000_000 + jst, 1_795_000_002 + jst, 1_794_999_999 + jst, 1_795_000_003 + jst]:
            self.assertTrue(watcher.is_local_time_of(printed, 1, jst, started, ended), printed)
        for printed in [1_794_999_998 + jst, 1_795_000_004 + jst, 1_795_000_001]:
            self.assertFalse(watcher.is_local_time_of(printed, 1, jst, started, ended), printed)
        # The date, to the minute: the minute of the run's start or end.
        minute = (1_795_000_000 + jst) // 60 * 60
        self.assertTrue(watcher.is_local_time_of(minute, 60, jst, started, ended))
        self.assertFalse(watcher.is_local_time_of(minute - 60, 60, jst, started, ended))
        self.assertFalse(watcher.is_local_time_of(minute - jst, 60, jst, started, ended))
        # UTC is the zone of offset 0.
        self.assertTrue(watcher.is_local_time_of(1_795_000_001, 1, 0, started, ended))


if __name__ == "__main__":
    unittest.main()
