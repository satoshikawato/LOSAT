#!/usr/bin/env python3
"""Tests for ci_fast_regressions.py and frozen_allowlist.py (no LOSAT binary is run)."""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import blastn_regression_fixtures as fixtures  # noqa: E402
import ci_fast_regressions as fast  # noqa: E402
import db_search_fixtures as db_fixtures  # noqa: E402
import tblastx_regression_fixtures as tblastx_fixtures  # noqa: E402
from frozen_allowlist import Allowlist  # noqa: E402

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
        for r in rows:
            if r["kind"] != "ncbi":
                continue
            data = (db_fixtures.FROZEN / f"{r['row_id']}.out").read_bytes()
            self.assertEqual(db_fixtures.sha256(data), r["stdout_sha256"], r["row_id"])
            self.assertEqual(str(len(data)), r["stdout_bytes"], r["row_id"])
            err = db_fixtures.FROZEN / f"{r['row_id']}.err"
            self.assertEqual(db_fixtures.sha256(err.read_bytes()) if err.exists() else "", r["stderr_sha256"], r["row_id"])


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
        self.assertFalse(db_fixtures.is_env_token("NCBI=1"))

    def test_clean_env_drops_what_changes_a_result_and_seals_home(self):
        dirty = {"PATH": "/bin", "HOME": "/home/someone", "USERPROFILE": "C:/Users/someone",
                 "NCBI_CONFIG__BLAST__BLASTDB": "/x", "NCBI_CONFIG_PATH": "/y", "NCBI": "/z",
                 "BLASTDB": "/d", "BLAST_USE_NCBI_ENV": "yes", "LOSAT_X_SOMETHING": "1",
                 "BATCH_SIZE": "1", "CHUNK_SIZE": "2", "OLD_FSC": "1", "BL2SEQ_LEGACY": "1", "LANG": "C"}
        clean = db_fixtures.clean_env(dirty, "/empty/home")
        self.assertEqual(clean, {"PATH": "/bin", "HOME": "/empty/home", "USERPROFILE": "/empty/home",
                                 "LANG": "C"})
        # nothing of the caller's environment is changed
        self.assertEqual(dirty["HOME"], "/home/someone")


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


if __name__ == "__main__":
    unittest.main()
