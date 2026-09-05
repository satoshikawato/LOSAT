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
from unittest import mock


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
    "warmup_count",
    "sample_index",
    "wall_seconds",
    "output_sha256",
    "benchmark_losat_sha",
    "environment_id",
    "tool",
    "contract",
    "benchmark_timestamp",
    "material_environment_id",
    "collection_segment",
    "boot_id",
    "effective_thread_label",
    "effective_thread_evidence",
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
                    "provenance_groups": {
                        "unit_current": {
                            "benchmark_losat_sha": "candidate-sha",
                            "environment": {
                                "cpu_model": "unit-test CPU",
                                "os": "unit-test OS",
                            },
                            "ncbi_version": "2.17.0+",
                            "used_for_current_plot": True,
                        }
                    },
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
                + [""] * 14
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
                        + [""] * 14
                    )
            # NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:191-194
            # CLocalBlast lcl_blast(queries, opts_hndl, db_adapter);
            # lcl_blast.SetNumberOfThreads(m_CmdLineArgs->GetNumThreads());
            # results = lcl_blast.Run();
            for mode_index, (mode, implementation, threads, effective) in enumerate(
                (
                    ("ncbi_n1", "NCBI BLAST+", "1", "requested/configured 1"),
                    ("ncbi_n8", "NCBI BLAST+", "8", "requested/configured 8"),
                    ("losat_native_n1", "LOSAT native", "1", "requested/configured 1"),
                    ("losat_native_n8", "LOSAT native", "8", "requested/configured 8"),
                    ("losat_wasm_serial", "LOSAT serial Wasm", "1", "serial command-Wasm"),
                    (
                        "losat_wasm_threads_requested_n8",
                        "LOSAT threaded Wasm requested n8",
                        "8",
                        "effective serial",
                    ),
                ),
                start=1,
            ):
                for sample_index in range(1, 6):
                    seconds = mode_index + sample_index / 10
                    writer.writerow(
                        [
                            "unit_current",
                            "tblastx",
                            mode,
                            "current-case",
                            implementation,
                            threads,
                            str(seconds),
                            "wall_clock_sample",
                            "candidate-sha",
                            "2.17.0+",
                            "2026-09-03",
                            "",
                            "output-hash",
                            str(seconds),
                            "1",
                            str(sample_index),
                            str(seconds),
                            "output-hash",
                            "candidate-sha",
                            "unit-environment",
                            implementation,
                            "EXACT_TEXT",
                            f"2026-09-03T00:00:0{sample_index}+00:00",
                            "unit-material-environment",
                            "1",
                            "unit-boot",
                            effective,
                            "unit-probe",
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

    def test_current_timing_uses_all_five_samples_and_median(self) -> None:
        rows = RENDERER.load_timing_rows(self.snapshot / "execution_times.tsv")
        current = [row for row in rows if row["provenance_id"] == "unit_current"]
        summaries = RENDERER.summarize_current_timing(current)
        self.assertEqual(len(summaries), 6)
        self.assertEqual([summary["n"] for summary in summaries], [5] * 6)
        self.assertEqual(
            [summary["median"] for summary in summaries],
            [1.3, 2.3, 3.3, 4.3, 5.3, 6.3],
        )
        self.assertEqual(
            [summary["min"] for summary in summaries],
            [1.1, 2.1, 3.1, 4.1, 5.1, 6.1],
        )
        self.assertEqual(
            [summary["max"] for summary in summaries],
            [1.5, 2.5, 3.5, 4.5, 5.5, 6.5],
        )

    # NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blastn_args.cpp:48,55-61
    # static const char kProgram[] = "blastn";
    # static const char kDefaultTask[] = "megablast";
    # SetTask(kDefaultTask);
    # arg.Reset(new CTaskCmdLineArgs(tasks, kDefaultTask));
    # Validate presentation of the recorded task identities and samples only.
    def test_current_snapshot_has_independent_linear_facets_and_all_modes(self) -> None:
        snapshot = ROOT / "benchmarks" / "v0.1.0"
        metadata, _ = RENDERER.load_metadata(snapshot)
        expected = json.loads((snapshot / "plot_data.json").read_text())[
            "plots"
        ]["execution_time"]
        figures = []
        original_figure = RENDERER.plt.figure

        def capture_figure(*args, **kwargs):
            figure = original_figure(*args, **kwargs)
            figures.append(figure)
            return figure

        with mock.patch.object(
            RENDERER.plt, "figure", side_effect=capture_figure
        ):
            _, plot_data = RENDERER.render_execution_time(
                snapshot / "execution_times.tsv",
                Path(self.temporary_directory.name) / "execution-time.png",
                metadata["datasets"]["execution_times"],
            )

        self.assertEqual(len(figures), 1)
        figure = figures[0]
        expected_facets = {
            "TBLASTX": [
                "p03_mela_pemojnva",
                "d06_ap027131_ap027133_db4",
                "p11_avclpv_psclpv",
            ],
            "BLASTN": ["PesePMNV.MjPMNV.task_blastn"],
            "Megablast": ["Sakai.MG1655.megablast"],
            "BLASTP": ["pairwise_default_serial"],
        }
        self.assertEqual(
            [axis.get_title() for axis in figure.axes], list(expected_facets)
        )
        self.assertEqual(len(figure.legends), 1)
        self.assertEqual(
            [text.get_text() for text in figure.legends[0].get_texts()],
            [RENDERER.TIMING_MODE_LABELS[mode] for mode in RENDERER.TIMING_MODES],
        )
        self.assertEqual(len({axis.get_xlim() for axis in figure.axes}), 4)
        self.assertEqual(sum(len(axis.patches) for axis in figure.axes), 36)
        summaries = expected["current"]["summaries"]
        lookup = {(row["case_id"], row["mode"]): row for row in summaries}
        for axis in figure.axes:
            self.assertEqual(axis.get_xscale(), "linear")
            self.assertEqual(axis.get_yscale(), "linear")
            self.assertEqual(axis.get_xlabel(), "Execution time (s)")
            self.assertEqual(axis.get_xlim()[0], 0)
            self.assertEqual(list(axis.get_shared_x_axes().get_siblings(axis)), [axis])
            self.assertIsNone(axis.get_legend())
            self.assertEqual(axis.xaxis.get_offset_text().get_text(), "")
            cases = expected_facets[axis.get_title()]
            self.assertEqual(len(axis.get_yticklabels()), len(cases))
            self.assertEqual(len({bar.get_y() for bar in axis.patches}), len(cases) * 6)
            for mode in RENDERER.TIMING_MODES:
                bars = next(
                    container for container in axis.containers
                    if container.get_label() == RENDERER.TIMING_MODE_LABELS[mode]
                )
                self.assertEqual(len(bars), len(cases))
                segments = bars.errorbar.lines[2][0].get_segments()
                self.assertEqual(len(segments), len(cases))
                for case, bar, segment in zip(cases, bars, segments):
                    summary = lookup[(case, mode)]
                    self.assertEqual(bar.get_x(), 0)
                    self.assertAlmostEqual(bar.get_width(), summary["median"])
                    self.assertAlmostEqual(segment[0][0], summary["min"])
                    self.assertAlmostEqual(segment[1][0], summary["max"])
                    self.assertLess(summary["max"], axis.get_xlim()[1])
                    if axis.get_title() == "BLASTN":
                        self.assertGreater(bar.get_width() / axis.get_xlim()[1], 0.25)
            facet_max = max(
                lookup[(case, mode)]["max"]
                for case in cases for mode in RENDERER.TIMING_MODES
            )
            self.assertLess(axis.get_xlim()[1], facet_max * 1.25)
        self.assertEqual(len(plot_data["current"]["samples"]), 180)
        self.assertEqual(len(plot_data["current"]["summaries"]), 36)
        self.assertEqual([row["n"] for row in summaries], [5] * 36)
        self.assertEqual(plot_data, expected)

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
