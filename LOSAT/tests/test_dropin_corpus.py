#!/usr/bin/env python3
"""Tests for dropin_corpus.py (no LOSAT binary and no NCBI program is run)."""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))

import dropin_corpus as corpus  # noqa: E402


def field(text: bytes) -> bytes:
    return len(text).to_bytes(4, "big") + text


def index(version: int, date: bytes, rest: bytes = b"offsets") -> bytes:
    head = version.to_bytes(4, "big") + bytes(4)
    if version >= 5:
        return head + bytes(4) + field(b"title") + field(b"db.ndb") + field(date) + rest
    return head + field(b"title") + field(date) + rest


def row(check: str = "bytes", **cells: str) -> dict[str, str]:
    base = dict.fromkeys(corpus.COLUMNS, "")
    return base | {"row_id": "R1", "check": check, "rc": "0", "stderr_empty": "yes", "command": "blastn -query q"} | cells


class DatabaseBytes(unittest.TestCase):
    def test_the_index_drops_only_the_date_in_both_versions(self):
        for version in (4, 5):
            one, two = index(version, b"Nov 1, 2026  9:05 AM\0\0\0\0"), index(version, b"Nov 11, 2026  10:05 PM\0\0")
            self.assertEqual(corpus.database_bytes("db.nin", one), corpus.database_bytes("db.pin", two))
            self.assertNotEqual(corpus.database_bytes("db.nin", index(version, b"d", b"other")),
                                corpus.database_bytes("db.nin", index(version, b"d")))
            self.assertNotIn(b"2026", corpus.database_bytes("db.nin", one))

    def test_the_json_manifest_drops_the_time_and_the_sizes(self):
        one = {"dbname": "db", "last-updated": "a", "bytes-total": 1, "bytes-to-cache": 2}
        two = {"dbname": "db", "last-updated": "b", "bytes-total": 9, "bytes-to-cache": 8}
        self.assertEqual(corpus.database_bytes("db.njs", json.dumps(one).encode()),
                         corpus.database_bytes("db.pjs", json.dumps(two).encode()))
        self.assertNotEqual(corpus.database_bytes("db.njs", json.dumps(one).encode()),
                            corpus.database_bytes("db.njs", json.dumps(one | {"dbname": "x"}).encode()))

    def test_other_files_are_compared_whole(self):
        self.assertEqual(corpus.database_bytes("db.nsq", b"Nov 1"), b"Nov 1")

    def test_a_logfile_counts_by_name_in_a_database_row_only(self):
        run = corpus.Run(0, b"", b"", {"mk.log": b"time"})
        self.assertEqual(corpus.files_cell(row(command="makeblastdb -in a -dbtype nucl"), run), "mk.log:-")
        self.assertEqual(corpus.files_cell(row(), run), f"mk.log:{corpus.sha256(b'time')}")


class Compare(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        patcher = mock.patch.object(corpus, "FROZEN", Path(directory.name))
        patcher.start()
        self.addCleanup(patcher.stop)
        (corpus.FROZEN / "R1.out").write_bytes(b"one\ntwo\n")
        (corpus.FROZEN / "R1--a.tsv").write_bytes(b"row\n")
        self.bytes_row = row(stdout_sha256=corpus.sha256(b"one\ntwo\n"), files=f"a.tsv:{corpus.sha256(b'row')}")
        (corpus.FROZEN / "R1--a.tsv").write_bytes(b"row")

    def test_equal_bytes_are_the_same(self):
        self.assertEqual(corpus.compare(self.bytes_row, corpus.Run(0, b"one\ntwo\n", b"", {"a.tsv": b"row"})), "")

    def test_each_part_of_a_row_can_differ(self):
        same_file = {"a.tsv": b"row"}
        self.assertIn("exit status 1", corpus.compare(self.bytes_row, corpus.Run(1, b"one\ntwo\n", b"x", same_file)))
        self.assertIn("stderr", corpus.compare(self.bytes_row, corpus.Run(0, b"one\ntwo\n", b"x", same_file)))
        self.assertIn("line 2", corpus.compare(self.bytes_row, corpus.Run(0, b"one\nTWO\n", b"", same_file)))
        self.assertIn("a.tsv", corpus.compare(self.bytes_row, corpus.Run(0, b"one\ntwo\n", b"", {"a.tsv": b"other"})))
        self.assertIn("files", corpus.compare(self.bytes_row, corpus.Run(0, b"one\ntwo\n", b"", {})))

    def test_stderr_is_allowed_where_ncbi_writes_some(self):
        quiet = self.bytes_row | {"stderr_empty": "no"}
        self.assertEqual(corpus.compare(quiet, corpus.Run(0, b"one\ntwo\n", b"Warning", {"a.tsv": b"row"})), "")

    def test_a_frozen_file_that_left_the_manifest_is_reported(self):
        (corpus.FROZEN / "R1.out").write_bytes(b"changed")
        self.assertIn("does not match", corpus.compare(self.bytes_row, corpus.Run(0, b"changed", b"", {"a.tsv": b"row"})))

    def test_line1_wants_ncbis_first_line_and_one_more(self):
        one = row("line1", stdout_sha256=corpus.sha256(b"one\ntwo\n"))
        self.assertEqual(corpus.compare(one, corpus.Run(0, b"one\nnames LOSAT\n", b"", {})), "")
        self.assertTrue(corpus.compare(one, corpus.Run(0, b"other\nnames LOSAT\n", b"", {})))
        self.assertTrue(corpus.compare(one, corpus.Run(0, b"one\n", b"", {})))

    def test_contains_wants_every_text(self):
        one = row("contains:-num_threads |-mt_mode ")
        self.assertEqual(corpus.compare(one, corpus.Run(0, b"-num_threads 1 -mt_mode 0", b"", {})), "")
        self.assertIn("-mt_mode", corpus.compare(one, corpus.Run(0, b"-num_threads 1", b"", {})))


class Fixtures(unittest.TestCase):
    def test_the_manifest_and_the_frozen_files_agree(self):
        _, rows = corpus.read_manifest()
        expected = {".gitattributes", ".gitignore"}
        for one in rows:
            self.assertEqual(one["rc"], "0", one["row_id"])
            if one["check"].startswith("contains:"):
                continue
            expected.add(f"{one['row_id']}.out")
            stdout = (corpus.FROZEN / f"{one['row_id']}.out").read_bytes()
            self.assertEqual(corpus.sha256(stdout), one["stdout_sha256"], one["row_id"])
            if corpus.is_database_row(corpus.shlex.split(one["command"])):
                continue
            for cell in filter(None, one["files"].split(",")):
                name, digest = cell.split(":")
                expected.add(f"{one['row_id']}--{name}")
                frozen = corpus.FROZEN / f"{one['row_id']}--{name}"
                self.assertEqual(corpus.sha256(frozen.read_bytes()), digest, frozen.name)
        self.assertEqual({path.name for path in corpus.FROZEN.iterdir()}, expected)


if __name__ == "__main__":
    unittest.main()
