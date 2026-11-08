#!/usr/bin/env python3
"""Tests for frozen_paths.py (no LOSAT binary is run)."""
from __future__ import annotations

import json
import os
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import frozen_paths  # noqa: E402
from frozen_paths import FROZEN_PREFIX, RUN_ROOT, to_frozen, to_run  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
STAGE_G = REPO / "docs/evidence/tlosan_stage_g"
F, R = FROZEN_PREFIX.encode(), RUN_ROOT.encode()
SUBJECT = "docs/evidence/tlosan_stage_d/all_codes_20260925/fixtures/code1.fna"


class ToRunTests(unittest.TestCase):
    def test_a_path_under_the_prefix_starts_at_the_run_root(self):
        self.assertEqual(to_run(f"{FROZEN_PREFIX}/{SUBJECT}"), f"{RUN_ROOT}/{SUBJECT}")
        self.assertEqual(len(to_run(f"{FROZEN_PREFIX}/{SUBJECT}")), len(f"{FROZEN_PREFIX}/{SUBJECT}"))

    def test_the_run_root_is_the_working_directory(self):
        self.assertEqual(os.path.normpath(RUN_ROOT), ".")
        self.assertEqual(os.path.normpath(to_run(f"{FROZEN_PREFIX}/{SUBJECT}")), SUBJECT)

    def test_other_arguments_are_kept(self):
        for arg in ("-subject", "0", SUBJECT, FROZEN_PREFIX, f"{FROZEN_PREFIX}2/{SUBJECT}",
                    f"/x{FROZEN_PREFIX}/{SUBJECT}", f"{REPO}/{SUBJECT}"):
            self.assertEqual(to_run(arg), arg)

    def test_every_stage_g_input_is_a_file_of_this_checkout(self):
        rewritten = 0
        with open(STAGE_G / "matrix_162.jsonl") as handle:
            for line in handle:
                for arg in json.loads(line)["losat_command"][1:]:
                    if arg.startswith(FROZEN_PREFIX + "/"):
                        self.assertTrue((REPO / to_run(arg)).is_file(), arg)
                        rewritten += 1
        self.assertEqual(rewritten, 324)  # -query and -subject of 162 rows


class ToFrozenTests(unittest.TestCase):
    def test_bytes_without_the_run_root_are_kept(self):
        for data in (b"", b"# Database: User specified sequence set (Input: tests/fasta/a.fna)\n",
                     b"Query_1\tsubject\t100.000\n"):
            self.assertEqual(to_frozen(data), data)

    def test_every_run_root_is_replaced(self):
        data = ((b"# Database: User specified sequence set (Input: " + R + b"/a.fna)\n") * 3
                + b"Error: " + R + b"/b.faa and " + R + b"/c.faa\n")
        self.assertEqual(to_frozen(data), data.replace(R + b"/", F + b"/"))
        self.assertEqual(to_frozen(data).count(F + b"/"), 5)

    def test_look_alikes_are_kept(self):
        for data in (F + b"/a.fna\n", F + b"2/a.fna\n", b"/mnt/c/Users/genom/GitHub/LOSAT2/a.fna\n",
                     f"{REPO}/{SUBJECT}\n".encode(),
                     R + b"\n",  # the run root without a path under it
                     b"./" * 15 + b"a.fna\n",  # a shorter run
                     b"/x/" + R + b"/a.fna\n",  # inside another path
                     b"." + R + b"/a.fna\n"):  # inside ../././...
            self.assertEqual(to_frozen(data), data, data)

    def test_a_longer_run_is_replaced_once_at_its_start(self):
        # `./` x 17 is the run root's `./a.fna`; the run root inside it (after `/`) is kept.
        self.assertEqual(to_frozen(b"./" * 17 + b"a.fna"), F + b"/./a.fna")
        self.assertEqual(to_frozen(b"./" * 32 + b"a.fna"), F + b"/" + b"./" * 16 + b"a.fna")

    def test_a_wrapped_outfmt_0_block_gives_the_frozen_bytes(self):
        # How LOSAT writes the `Database:` lines for `-subject <RUN_ROOT>/<SUBJECT>`: wrapped at 68
        # columns (NStr::Wrap), where NCBI wrapped them for the frozen path, which has the same length.
        run = (b"Database: User specified sequence set (Input:\n"
               b"././././././././././././././././docs/evidence/tlosan_stage_d/all_cod\n"
               b"es_20260925/fixtures/code1.fna).\n")
        frozen = (STAGE_G / "batch_boundary/valid_then_skipped.fmt0.ncbi").read_bytes().splitlines(keepends=True)
        self.assertEqual(to_frozen(run), b"".join(frozen[10:13]))
        self.assertEqual(len(run.splitlines()[1]), 68)


class CaseTableTests(unittest.TestCase):
    def test_the_entry_is_to_frozen_as_data(self):
        entry = frozen_paths.case_table_entry()
        self.assertEqual(entry, {"run": RUN_ROOT + "/", "frozen": FROZEN_PREFIX + "/", "not_after": "/."})


if __name__ == "__main__":
    unittest.main()
