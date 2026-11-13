#!/usr/bin/env python3
"""Tests for frozen_paths.py (no LOSAT binary is run)."""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
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

    def test_every_stage_g_input_is_the_file_stage_g_recorded(self):
        rewritten = 0
        with open(STAGE_G / "matrix_162.jsonl") as handle:
            for line in handle:
                row = json.loads(line)
                recorded = {row["query"]: row["query_sha256"], row["subject"]: row["subject_sha256"]}
                for arg in row["losat_command"][1:]:
                    if arg.startswith(FROZEN_PREFIX + "/"):
                        local = REPO / to_run(arg)
                        self.assertEqual(hashlib.sha256(local.read_bytes()).hexdigest(), recorded[arg],
                                         f"{row['case']}: {local.relative_to(REPO)}")
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


# 2026-11-13: two fixtures named `nul.fa` made `actions/checkout` fail on Windows
# (`error: invalid path 'LOSAT/tests/convert2blastmask/in/nul.fa'`): Windows treats the base name
# of NUL, CON, ... as a device whatever the extension. Linux and macOS accept such paths, so only
# a check of the tracked names keeps them out.
WINDOWS_RESERVED = frozenset(["CON", "PRN", "AUX", "NUL"]
                             + [f"COM{n}" for n in range(1, 10)] + [f"LPT{n}" for n in range(1, 10)])


def is_windows_reserved(component: str) -> bool:
    """True when `component` (one path component) names a Windows device: the part before the
    first `.`, trailing spaces and dots ignored, any case."""
    return component.rstrip(" .").split(".", 1)[0].rstrip(" ").upper() in WINDOWS_RESERVED


def tracked_paths() -> list[str]:
    """Every tracked path of the repository (`git ls-files`); the directory walk without `.git`
    when git is not available."""
    try:
        out = subprocess.run(["git", "ls-files", "-z"], cwd=REPO, check=True,
                             stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout
        return [p for p in out.decode("utf-8", "surrogateescape").split("\0") if p]
    except (OSError, subprocess.CalledProcessError):
        found = []
        for root, dirs, files in os.walk(REPO):
            dirs[:] = [d for d in dirs if d != ".git"]
            found.extend(os.path.relpath(os.path.join(root, f), REPO).replace(os.sep, "/")
                         for f in dirs + files)
        return found


class WindowsReservedNameTests(unittest.TestCase):
    def test_reserved_names_are_found(self):
        for name in ("nul.fa", "NUL", "Com1.txt", "lpt9.tar.gz", "aux", "nul .fa", "nul."):
            self.assertTrue(is_windows_reserved(name), name)

    def test_other_names_are_not(self):
        for name in ("nul_bytes.fa", "null.fa", "console.txt", "com10", "lpt0", "auxiliary.fa", "a.nul", ""):
            self.assertFalse(is_windows_reserved(name), name)

    def test_no_tracked_path_has_a_windows_reserved_component(self):
        bad = [p for p in tracked_paths() if any(is_windows_reserved(c) for c in p.split("/"))]
        self.assertEqual(bad, [], "Windows cannot check these paths out; rename them: " + ", ".join(bad))


if __name__ == "__main__":
    unittest.main()
