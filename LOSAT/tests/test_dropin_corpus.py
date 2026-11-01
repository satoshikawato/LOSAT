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


def index(version: int, date: bytes = b"Nov 1, 2026  9:05 AM", rest: bytes = b"offsets", aligned: bool = True) -> bytes:
    head = version.to_bytes(4, "big") + bytes(4)
    head += bytes(4) + field(b"title") + field(b"db.ndb") if version >= 5 else field(b"title")
    date += bytes(-(len(head) + 4 + len(date)) % 8 + (0 if aligned else 1))
    return head + field(date) + rest


def row(check: str = "bytes", **cells: str) -> dict[str, str]:
    base = dict.fromkeys(corpus.COLUMNS, "")
    return base | {"row_id": "R1", "check": check, "rc": "0", "stderr_empty": "yes", "command": "blastn -query q"} | cells


class DatabaseBytes(unittest.TestCase):
    def index_of(self, data: bytes) -> bytes:
        return corpus.database_bytes("db.nin", {"db.nin": data})

    def test_the_index_drops_only_the_date_in_both_versions(self):
        for version in (4, 5):
            self.assertEqual(self.index_of(index(version)), self.index_of(index(version, b"Nov 11, 2026  10:05 PM")))
            self.assertNotEqual(self.index_of(index(version)), self.index_of(index(version, rest=b"other")))
            self.assertNotIn(b"2026", self.index_of(index(version)))

    def test_another_date_form_or_a_lost_alignment_keeps_the_file_whole(self):
        for version in (4, 5):
            for bad in (index(version, b"2026-11-01T09:05:00"), index(version, b""), index(version, aligned=False)):
                self.assertEqual(self.index_of(bad), bad)

    def test_the_json_manifest_drops_the_time_and_compares_the_sizes(self):
        files = {"db.nin": b"12", "db.nsq": b"345", "db.not": b"6789"}

        def manifest(total=9, cache=5, time=b"2026-11-01T22:00:00", order=False):
            head = b'"dbname": "db",\n  "last-updated": "%s",\n' % time
            sizes = b'  "bytes-total": %d,\n  "bytes-to-cache": %d,\n' % (total, cache)
            listing = b'  "files": ["db.nin", "db.nsq", "db.not"]\n'
            return b"{\n  " + (sizes + b"  " + head if order else head + sizes) + listing + b"}\n"

        def view(data):
            return corpus.database_bytes("db.njs", files | {"db.njs": data})

        same = view(manifest())
        self.assertEqual(view(manifest(time=b"2027-01-02T03:04:05")), same)
        grown = files | {"db.nin": b"12" + bytes(8)}
        self.assertEqual(corpus.database_bytes("db.pjs", grown | {"db.pjs": manifest(total=17, cache=13)}), same)
        self.assertNotEqual(view(manifest(total=10)), same)
        self.assertNotEqual(view(manifest(cache=6)), same)
        self.assertNotEqual(view(manifest(order=True)), same)
        self.assertNotEqual(view(manifest(time=b"Nov 1, 2026")), same)
        self.assertEqual(view(b"not json"), b"not json")

    def test_other_files_are_compared_whole(self):
        self.assertEqual(corpus.database_bytes("db.nsq", {"db.nsq": b"Nov 1"}), b"Nov 1")

    def test_a_logfile_is_the_report_without_the_time_and_the_path(self):
        files = {"mk.log": b"\nBuilding a new DB, current time: x\nNew DB name:   /tmp/x\nNew DB title:  t\n"}
        self.assertEqual(corpus.database_bytes("mk.log", files), b"\nNew DB title:  t\n")
        database = row(command="makeblastdb -in a -dbtype nucl")
        filled = corpus.files_cell(database, corpus.Run(0, b"", b"", files))
        self.assertNotEqual(corpus.files_cell(database, corpus.Run(0, b"", b"", {"mk.log": b""})), filled)

    def test_a_search_row_compares_a_log_named_file_whole(self):
        self.assertEqual(corpus.files_cell(row(), corpus.Run(0, b"", b"", {"a.log": b"x"})), f"a.log:{corpus.sha256(b'x')}")


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

    def test_a_run_that_changes_an_input_differs(self):
        run = corpus.Run(0, b"one\ntwo\n", b"", {"a.tsv": b"row"}, ["ns.fna"])
        self.assertIn("ns.fna", corpus.compare(self.bytes_row, run))

    def test_stderr_is_allowed_where_ncbi_writes_some(self):
        quiet = self.bytes_row | {"stderr_empty": "no"}
        self.assertEqual(corpus.compare(quiet, corpus.Run(0, b"one\ntwo\n", b"Warning", {"a.tsv": b"row"})), "")

    def test_a_frozen_file_that_left_the_manifest_is_reported(self):
        (corpus.FROZEN / "R1.out").write_bytes(b"changed")
        self.assertIn("does not match", corpus.compare(self.bytes_row, corpus.Run(0, b"changed", b"", {"a.tsv": b"row"})))

    def test_line1_wants_ncbis_first_line_and_one_more(self):
        one = row("line1", stdout_sha256=corpus.sha256(b"one\ntwo\n"))
        self.assertEqual(corpus.compare(one, corpus.Run(0, b"one\nnames LOSAT\n", b"", {})), "")
        for stdout in (b"other\nnames LOSAT\n", b"one\n", b"one\nnames NCBI\n", b"one\nLOSAT\nmore\n"):
            self.assertTrue(corpus.compare(one, corpus.Run(0, stdout, b"", {})), stdout)

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
