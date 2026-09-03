"""Tests for the data-only benchmark renderer."""

from __future__ import annotations

import ast
import csv
import gzip
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RENDERER_PATH = ROOT / "scripts" / "render_benchmark_plots.py"
SPEC = importlib.util.spec_from_file_location("render_benchmark_plots", RENDERER_PATH)
assert SPEC and SPEC.loader
RENDERER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RENDERER)

# NCBI reference: ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:38-40,109-116
# defines the default tabular alignment fields and identifies `length` and
# `pident`; the fixture below exercises LOSAT's normalized representation only.
ALIGNMENT_HEADER = [
    "program",
    "case_id",
    "implementation",
    "classification",
    "primary_for_distribution",
    "source_row",
    "qseqid",
    "sseqid",
    "pident",
    "length",
]
TIMING_HEADER = [
    "provenance_id",
    "program",
    "mode",
    "case_id",
    "implementation",
    "thread_count",
    "seconds",
    "value_kind",
    "producing_sha",
    "ncbi_version",
    "recorded_date",
    "source_path",
    "source_sha256",
    "raw_elapsed",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class RendererTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.snapshot = Path(self.temporary_directory.name) / "snapshot"
        self.snapshot.mkdir()
        self._write_alignment_fixture()
        self._write_timing_fixture()
        metadata = {
            "schema_version": 1,
            "snapshot_id": "unit-test-snapshot",
            "datasets": {
                "alignment_results": {
                    "file": "alignment_results.tsv.gz",
                    "sha256": sha256(self.snapshot / "alignment_results.tsv.gz"),
                },
                "execution_times": {
                    "file": "execution_times.tsv",
                    "sha256": sha256(self.snapshot / "execution_times.tsv"),
                },
            },
        }
        (self.snapshot / "metadata.json").write_text(
            json.dumps(metadata), encoding="utf-8"
        )

    def _write_alignment_fixture(self) -> None:
        path = self.snapshot / "alignment_results.tsv.gz"
        with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(ALIGNMENT_HEADER)
            source_row = 0
            for program_index, program in enumerate(("blastn", "blastp", "tblastx")):
                for implementation_index, implementation in enumerate(
                    ("NCBI BLAST+", "LOSAT")
                ):
                    for point_index in range(3):
                        source_row += 1
                        writer.writerow(
                            [
                                program,
                                f"case-{program}",
                                implementation,
                                "EXACT_TEXT",
                                "true",
                                source_row,
                                "query",
                                "subject",
                                80 + program_index + implementation_index + point_index,
                                10 + program_index * 4 + point_index,
                            ]
                        )

    def _write_timing_fixture(self) -> None:
        path = self.snapshot / "execution_times.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(TIMING_HEADER)
            writer.writerow(
                [
                    "historical_plot_2026_06",
                    "tblastx",
                    "TBLASTX",
                    "historic-case",
                    "NCBI BLAST+",
                    "1",
                    "2.5",
                    "wall_time",
                    "",
                    "",
                    "",
                    "old.log",
                    "abc",
                    "0:02.50",
                ]
            )
            for case_id, baseline, candidate in (
                ("p11_avclpv_psclpv", "639.21", "636.70"),
                ("d06_ap027131_ap027133_db4", "60.26", "60.65"),
            ):
                for implementation, seconds in (
                    ("LOSAT baseline", baseline),
                    ("LOSAT candidate", candidate),
                ):
                    writer.writerow(
                        [
                            "pr93_performance_lineage",
                            "tblastx",
                            "TBLASTX",
                            case_id,
                            implementation,
                            "",
                            seconds,
                            "reported_value",
                            "candidate-sha",
                            "",
                            "",
                            "performance_lineage.json",
                            "def",
                            "",
                        ]
                    )

    def test_render_is_deterministic_within_one_environment(self) -> None:
        first = Path(self.temporary_directory.name) / "first"
        second = Path(self.temporary_directory.name) / "second"
        RENDERER.render(self.snapshot, first)
        RENDERER.render(self.snapshot, second)

        for filename in (
            "execution_time.png",
            "hit_distribution.png",
            "plot_data.json",
            "render_manifest.json",
        ):
            self.assertEqual(sha256(first / filename), sha256(second / filename))

    def test_dataset_checksum_mismatch_fails_closed(self) -> None:
        timing_path = self.snapshot / "execution_times.tsv"
        with timing_path.open("a", encoding="utf-8") as handle:
            handle.write("tampered\n")
        metadata, _ = RENDERER.load_metadata(self.snapshot)
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            RENDERER.verified_dataset_path(self.snapshot, metadata, "execution_times")

    def test_snapshot_file_cannot_escape_snapshot(self) -> None:
        with self.assertRaisesRegex(ValueError, "escapes snapshot directory"):
            RENDERER.snapshot_file(self.snapshot, "../outside.tsv")

    def test_renderer_has_no_process_or_network_imports(self) -> None:
        tree = ast.parse(RENDERER_PATH.read_text(encoding="utf-8"))
        imported_roots = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported_roots.update(alias.name.split(".")[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported_roots.add(node.module.split(".")[0])
        self.assertTrue(
            imported_roots.isdisjoint(
                {"http", "requests", "socket", "subprocess", "urllib"}
            )
        )


if __name__ == "__main__":
    unittest.main()
