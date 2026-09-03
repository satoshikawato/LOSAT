#!/usr/bin/env python3
"""Render committed LOSAT benchmark snapshots without collecting new data."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterator

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch


RENDERER_VERSION = 1
PROGRAMS = ("blastn", "blastp", "tblastx")
IMPLEMENTATIONS = ("NCBI BLAST+", "LOSAT")
PROGRAM_LABELS = {"blastn": "BLASTN", "blastp": "BLASTP", "tblastx": "TBLASTX"}
ALIGNMENT_FIELDS = {
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
}
TIMING_FIELDS = {
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
}

# NCBI reference: ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:38-40,109-116
# const char* kDfltArgTabularOutputFmt =
#     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
#     "evalue bitscore";
# SFormatSpec("length", "Alignment length", eAlignmentLength),
# SFormatSpec("pident", "Percentage of identical matches", ePercentIdentical),


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def snapshot_file(snapshot: Path, filename: str) -> Path:
    candidate = (snapshot / filename).resolve()
    snapshot_root = snapshot.resolve()
    if candidate != snapshot_root and snapshot_root not in candidate.parents:
        raise ValueError(f"snapshot file escapes snapshot directory: {filename}")
    if not candidate.is_file():
        raise FileNotFoundError(candidate)
    return candidate


def load_metadata(snapshot: Path) -> tuple[dict[str, object], Path]:
    metadata_path = snapshot_file(snapshot, "metadata.json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("schema_version") != 1:
        raise ValueError("metadata.json must use schema_version 1")
    if not isinstance(metadata.get("datasets"), dict):
        raise ValueError("metadata.json is missing datasets")
    return metadata, metadata_path


def verified_dataset_path(
    snapshot: Path, metadata: dict[str, object], dataset_name: str
) -> Path:
    datasets = metadata["datasets"]
    if not isinstance(datasets, dict) or dataset_name not in datasets:
        raise ValueError(f"metadata.json is missing dataset {dataset_name}")
    dataset = datasets[dataset_name]
    if not isinstance(dataset, dict):
        raise ValueError(f"dataset {dataset_name} metadata must be an object")
    path = snapshot_file(snapshot, str(dataset["file"]))
    actual = file_sha256(path)
    expected = str(dataset["sha256"])
    if actual != expected:
        raise ValueError(
            f"{dataset_name} checksum mismatch: metadata={expected}, actual={actual}"
        )
    return path


def alignment_rows(path: Path) -> Iterator[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if set(reader.fieldnames or ()) != ALIGNMENT_FIELDS:
            raise ValueError("alignment_results.tsv.gz has an unexpected schema")
        yield from reader


def checked_alignment_values(row: dict[str, str]) -> tuple[str, str, float, int]:
    program = row["program"]
    implementation = row["implementation"]
    if program not in PROGRAMS:
        raise ValueError(f"unknown program: {program}")
    if implementation not in IMPLEMENTATIONS:
        raise ValueError(f"unknown implementation: {implementation}")
    pident = float(row["pident"])
    length = int(row["length"])
    if not math.isfinite(pident) or not 0.0 <= pident <= 100.0:
        raise ValueError(f"invalid pident: {row['pident']}")
    if length <= 0:
        raise ValueError(f"invalid alignment length: {row['length']}")
    return program, implementation, pident, length


def aggregate_alignments(
    path: Path, max_scatter_points: int = 6000
) -> dict[str, object]:
    counts: dict[tuple[str, str], int] = defaultdict(int)
    length_bounds = {program: [math.inf, -math.inf] for program in PROGRAMS}

    for row in alignment_rows(path):
        if row["primary_for_distribution"] != "true":
            continue
        program, implementation, _, length = checked_alignment_values(row)
        counts[(program, implementation)] += 1
        length_bounds[program][0] = min(length_bounds[program][0], length)
        length_bounds[program][1] = max(length_bounds[program][1], length)

    for program in PROGRAMS:
        if not all(
            counts[(program, implementation)] for implementation in IMPLEMENTATIONS
        ):
            raise ValueError(
                f"alignment snapshot has no plottable {program} comparison rows"
            )

    length_edges: dict[str, np.ndarray] = {}
    for program, (lower, upper) in length_bounds.items():
        if lower == upper:
            lower = max(1.0, lower * 0.9)
            upper = upper * 1.1
        length_edges[program] = np.geomspace(lower, upper, 51)
    identity_edges = np.linspace(0.0, 100.0, 51)
    length_hist = {
        key: np.zeros(50, dtype=float)
        for key in ((program, impl) for program in PROGRAMS for impl in IMPLEMENTATIONS)
    }
    identity_hist = {key: np.zeros(50, dtype=float) for key in length_hist}
    scatter: dict[tuple[str, str], list[tuple[int, float]]] = {
        key: [] for key in length_hist
    }
    seen: dict[tuple[str, str], int] = defaultdict(int)
    strides = {
        key: max(1, math.ceil(value / max_scatter_points))
        for key, value in counts.items()
    }

    for row in alignment_rows(path):
        if row["primary_for_distribution"] != "true":
            continue
        program, implementation, pident, length = checked_alignment_values(row)
        key = (program, implementation)
        log_min = math.log(length_edges[program][0])
        log_span = math.log(length_edges[program][-1]) - log_min
        length_index = int((math.log(length) - log_min) / log_span * 50)
        length_index = min(49, max(0, length_index))
        identity_index = min(49, max(0, int(pident / 100.0 * 50)))
        length_hist[key][length_index] += length
        identity_hist[key][identity_index] += length
        ordinal = seen[key]
        if ordinal % strides[key] == 0 and len(scatter[key]) < max_scatter_points:
            scatter[key].append((length, pident))
        seen[key] += 1

    return {
        "counts": counts,
        "length_edges": length_edges,
        "identity_edges": identity_edges,
        "length_hist": length_hist,
        "identity_hist": identity_hist,
        "scatter": scatter,
    }


def chart_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.edgecolor": "#4b5563",
            "axes.labelcolor": "#1f2937",
            "axes.titlecolor": "#111827",
            "axes.facecolor": "#ffffff",
            "figure.facecolor": "#ffffff",
            "grid.color": "#e5e7eb",
            "grid.linewidth": 0.7,
            "text.color": "#111827",
            "xtick.color": "#4b5563",
            "ytick.color": "#4b5563",
        }
    )


def render_hit_distribution(
    path: Path, output: Path
) -> tuple[dict[str, int], dict[str, object]]:
    aggregated = aggregate_alignments(path)
    colors = {"NCBI BLAST+": "#4c72b0", "LOSAT": "#dd8452"}
    styles = {"NCBI BLAST+": "-", "LOSAT": "--"}
    markers = {"NCBI BLAST+": "o", "LOSAT": "x"}
    fig, axes = plt.subplots(3, 3, figsize=(16, 13), constrained_layout=False)

    for row_index, program in enumerate(PROGRAMS):
        length_ax, identity_ax, scatter_ax = axes[row_index]
        for implementation in IMPLEMENTATIONS:
            key = (program, implementation)
            length_ax.stairs(
                aggregated["length_hist"][key],
                aggregated["length_edges"][program],
                color=colors[implementation],
                linestyle=styles[implementation],
                linewidth=1.7,
                label=implementation,
            )
            identity_ax.stairs(
                aggregated["identity_hist"][key],
                aggregated["identity_edges"],
                color=colors[implementation],
                linestyle=styles[implementation],
                linewidth=1.7,
                label=implementation,
            )
            points = aggregated["scatter"][key]
            x_values = [point[0] for point in points]
            y_values = [point[1] for point in points]
            if implementation == "NCBI BLAST+":
                scatter_ax.scatter(
                    x_values,
                    y_values,
                    s=11,
                    marker=markers[implementation],
                    facecolors="none",
                    edgecolors=colors[implementation],
                    linewidths=0.65,
                    alpha=0.42,
                    label=implementation,
                )
            else:
                scatter_ax.scatter(
                    x_values,
                    y_values,
                    s=11,
                    marker=markers[implementation],
                    color=colors[implementation],
                    linewidths=0.65,
                    alpha=0.36,
                    label=implementation,
                )

        length_ax.set_xscale("log")
        length_ax.set_yscale("log")
        length_ax.set_title("Alignment-length distribution")
        length_ax.set_xlabel("Alignment length (bp or aa; log scale)")
        length_ax.set_ylabel("Accumulated aligned length (log scale)")
        identity_ax.set_title("Percent-identity distribution")
        identity_ax.set_xlabel("Percent identity")
        identity_ax.set_ylabel("Accumulated aligned length")
        scatter_ax.set_xscale("log")
        scatter_ax.set_title("Alignment length vs percent identity")
        scatter_ax.set_xlabel("Alignment length (bp or aa; log scale)")
        scatter_ax.set_ylabel("Percent identity")
        for axis in (length_ax, identity_ax, scatter_ax):
            axis.grid(True, axis="both", alpha=0.85)
        length_ax.text(
            -0.2,
            0.5,
            PROGRAM_LABELS[program],
            transform=length_ax.transAxes,
            rotation=90,
            va="center",
            ha="center",
            fontsize=13,
            fontweight="bold",
        )

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(0.5, 0.945),
    )
    fig.suptitle(
        "Certified v0.1.0 alignment-output distributions", fontsize=18, y=0.987
    )
    fig.text(
        0.5,
        0.957,
        "Recovered PR 5 native evidence · LOSAT 5845d22 · NCBI BLAST+ 2.17.0 · contract-output weighting",
        ha="center",
        fontsize=10,
        color="#4b5563",
    )
    fig.text(
        0.5,
        0.012,
        "Includes six approved local-subject non-default db-gencode deviation contracts and one source-undetermined BLASTN contract; this is not a parity score. "
        "Duplicate BLASTP thread-4 outputs are retained in the snapshot but excluded here.",
        ha="center",
        va="bottom",
        fontsize=8.2,
        color="#4b5563",
        wrap=True,
    )
    fig.subplots_adjust(
        left=0.08, right=0.985, top=0.91, bottom=0.07, hspace=0.42, wspace=0.3
    )
    fig.savefig(
        output,
        dpi=140,
        bbox_inches="tight",
        metadata={
            "Software": f"LOSAT data-only benchmark renderer v{RENDERER_VERSION}"
        },
    )
    plt.close(fig)
    counts = {
        f"{program}:{implementation}": int(
            aggregated["counts"][(program, implementation)]
        )
        for program in PROGRAMS
        for implementation in IMPLEMENTATIONS
    }
    plot_data = {
        "counts": counts,
        "identity_edges": aggregated["identity_edges"].tolist(),
        "programs": {
            program: {
                "length_edges": aggregated["length_edges"][program].tolist(),
                "series": {
                    implementation: {
                        "identity_hist": aggregated["identity_hist"][
                            (program, implementation)
                        ].tolist(),
                        "length_hist": aggregated["length_hist"][
                            (program, implementation)
                        ].tolist(),
                        "scatter_count": len(
                            aggregated["scatter"][(program, implementation)]
                        ),
                        "scatter_sha256": hashlib.sha256(
                            json.dumps(
                                aggregated["scatter"][(program, implementation)],
                                separators=(",", ":"),
                            ).encode("utf-8")
                        ).hexdigest(),
                    }
                    for implementation in IMPLEMENTATIONS
                },
            }
            for program in PROGRAMS
        },
    }
    return counts, plot_data


def load_timing_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if set(reader.fieldnames or ()) != TIMING_FIELDS:
            raise ValueError("execution_times.tsv has an unexpected schema")
        rows = list(reader)
    for row in rows:
        value = float(row["seconds"])
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"invalid execution time: {row['seconds']}")
    return rows


def history_series(row: dict[str, str]) -> str:
    return f"{row['implementation']} n{row['thread_count']}"


def render_historical_timing_panel(axis, rows: list[dict[str, str]], mode: str) -> None:
    mode_rows = [row for row in rows if row["mode"] == mode]
    tasks = list(dict.fromkeys(row["case_id"] for row in mode_rows))
    series_order = [
        "NCBI BLAST+ n1",
        "NCBI BLAST+ n8",
        "LOSAT native n1",
        "LOSAT native n8",
        "LOSAT wasm n1",
        "LOSAT wasm n8",
    ]
    series = [
        name
        for name in series_order
        if any(history_series(row) == name for row in mode_rows)
    ]
    colors = {
        "NCBI BLAST+ n1": "#4c72b0",
        "NCBI BLAST+ n8": "#4c72b0",
        "LOSAT native n1": "#dd8452",
        "LOSAT native n8": "#e6a700",
        "LOSAT wasm n1": "#8172b3",
        "LOSAT wasm n8": "#937860",
    }
    hatches = {
        "NCBI BLAST+ n1": "",
        "NCBI BLAST+ n8": "//",
        "LOSAT native n1": "",
        "LOSAT native n8": "//",
        "LOSAT wasm n1": "",
        "LOSAT wasm n8": "//",
    }
    value_lookup = {
        (row["case_id"], history_series(row)): float(row["seconds"])
        for row in mode_rows
    }
    y_positions = np.arange(len(tasks), dtype=float)
    bar_height = min(0.16, 0.78 / max(1, len(series)))
    center = (len(series) - 1) / 2
    for index, name in enumerate(series):
        values = [value_lookup.get((task, name), np.nan) for task in tasks]
        axis.barh(
            y_positions + (index - center) * bar_height,
            values,
            height=bar_height,
            color=colors[name],
            edgecolor="#374151",
            linewidth=0.35,
            hatch=hatches[name],
            label=name,
        )
    axis.set_yticks(y_positions, tasks, fontsize=7.2)
    axis.invert_yaxis()
    axis.set_xlim(left=0)
    axis.set_xlabel("Recorded wall time (seconds)")
    axis.set_title(mode)
    axis.grid(True, axis="x")
    axis.grid(False, axis="y")


def render_execution_time(
    path: Path, output: Path
) -> tuple[dict[str, int], dict[str, object]]:
    rows = load_timing_rows(path)
    history = [row for row in rows if row["provenance_id"] == "historical_plot_2026_06"]
    pr93 = [row for row in rows if row["provenance_id"] == "pr93_performance_lineage"]
    if not history or not pr93:
        raise ValueError(
            "execution snapshot must contain historical and PR 93 provenance groups"
        )

    fig = plt.figure(figsize=(18, 14.5))
    grid = fig.add_gridspec(
        3, 2, height_ratios=(1.15, 1.15, 0.62), hspace=0.42, wspace=0.38
    )
    modes = ("TBLASTX", "Megablast", "BLASTN", "BLASTP")
    history_axes = [fig.add_subplot(grid[index // 2, index % 2]) for index in range(4)]
    for axis, mode in zip(history_axes, modes):
        render_historical_timing_panel(axis, history, mode)

    legend_order = [
        ("NCBI BLAST+ n1", "#4c72b0", ""),
        ("NCBI BLAST+ n8", "#4c72b0", "//"),
        ("LOSAT native n1", "#dd8452", ""),
        ("LOSAT native n8", "#e6a700", "//"),
        ("LOSAT wasm n1", "#8172b3", ""),
        ("LOSAT wasm n8", "#937860", "//"),
    ]
    legend_handles = [
        Patch(facecolor=color, edgecolor="#374151", hatch=hatch, label=label)
        for label, color, hatch in legend_order
        if any(history_series(row) == label for row in history)
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        ncol=6,
        frameon=False,
        bbox_to_anchor=(0.5, 0.93),
        fontsize=8.5,
    )

    lineage_axis = fig.add_subplot(grid[2, :])
    cases = list(dict.fromkeys(row["case_id"] for row in pr93))
    labels = {
        "p11_avclpv_psclpv": "p11 AvCLPV vs PsCLPV",
        "d06_ap027131_ap027133_db4": "d06 AP027131 vs AP027133 db-gencode 4",
    }
    lineage_series = ("LOSAT baseline", "LOSAT candidate")
    lineage_colors = {"LOSAT baseline": "#6b7280", "LOSAT candidate": "#dd8452"}
    lineage_hatches = {"LOSAT baseline": "", "LOSAT candidate": "//"}
    lookup = {
        (row["case_id"], row["implementation"]): float(row["seconds"]) for row in pr93
    }
    y_positions = np.arange(len(cases), dtype=float)
    for index, series in enumerate(lineage_series):
        values = [lookup[(case, series)] for case in cases]
        bars = lineage_axis.barh(
            y_positions + (index - 0.5) * 0.32,
            values,
            height=0.3,
            color=lineage_colors[series],
            edgecolor="#374151",
            linewidth=0.5,
            hatch=lineage_hatches[series],
            label=series,
        )
        lineage_axis.bar_label(
            bars, labels=[f"{value:.2f} s" for value in values], padding=4, fontsize=8
        )
    lineage_axis.set_yticks(y_positions, [labels[case] for case in cases])
    lineage_axis.invert_yaxis()
    lineage_axis.set_xlim(left=0)
    lineage_axis.set_xlabel("Retained reported value (seconds)")
    lineage_axis.set_title(
        "PR 93 TBLASTX lineage — raw samples, statistic definition, thread count, and commands unavailable"
    )
    lineage_axis.grid(True, axis="x")
    lineage_axis.grid(False, axis="y")
    lineage_axis.legend(frameon=False, loc="lower right")

    fig.suptitle("Recovered execution-time records", fontsize=18, y=0.985)
    fig.text(
        0.5,
        0.952,
        "Historical panel: May–June 2026 plot inputs · producing LOSAT SHA and complete NCBI identity unknown · not final v0.1.0",
        ha="center",
        fontsize=10,
        color="#4b5563",
    )
    fig.text(
        0.5,
        0.012,
        "PR 93 values compare base 3ba0024 with candidate d929055. No certification wall-clock time or invented sample is included.",
        ha="center",
        va="bottom",
        fontsize=8.5,
        color="#4b5563",
    )
    fig.subplots_adjust(left=0.16, right=0.975, top=0.90, bottom=0.07)
    fig.savefig(
        output,
        dpi=140,
        bbox_inches="tight",
        metadata={
            "Software": f"LOSAT data-only benchmark renderer v{RENDERER_VERSION}"
        },
    )
    plt.close(fig)
    counts = {
        "historical_plot_2026_06": len(history),
        "pr93_performance_lineage": len(pr93),
    }
    plot_data = {
        "counts": counts,
        "historical": [
            {
                "case_id": row["case_id"],
                "mode": row["mode"],
                "seconds": float(row["seconds"]),
                "series": history_series(row),
            }
            for row in history
        ],
        "pr93": [
            {
                "case_id": row["case_id"],
                "seconds": float(row["seconds"]),
                "series": row["implementation"],
            }
            for row in pr93
        ],
    }
    return counts, plot_data


def render(snapshot: Path, output: Path) -> dict[str, object]:
    chart_style()
    metadata, metadata_path = load_metadata(snapshot)
    alignment_path = verified_dataset_path(snapshot, metadata, "alignment_results")
    timing_path = verified_dataset_path(snapshot, metadata, "execution_times")
    output.mkdir(parents=True, exist_ok=True)
    hit_path = output / "hit_distribution.png"
    time_path = output / "execution_time.png"
    hit_counts, alignment_plot_data = render_hit_distribution(alignment_path, hit_path)
    timing_counts, timing_plot_data = render_execution_time(timing_path, time_path)
    plot_data_path = output / "plot_data.json"
    plot_data = {
        "schema_version": 1,
        "renderer_version": RENDERER_VERSION,
        "renderer_source_sha256": file_sha256(Path(__file__).resolve()),
        "snapshot_id": metadata["snapshot_id"],
        "inputs": {
            "metadata.json": file_sha256(metadata_path),
            alignment_path.name: file_sha256(alignment_path),
            timing_path.name: file_sha256(timing_path),
        },
        "plots": {
            "execution_time": timing_plot_data,
            "hit_distribution": alignment_plot_data,
        },
    }
    plot_data_path.write_text(
        json.dumps(plot_data, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest = {
        "schema_version": 1,
        "renderer_version": RENDERER_VERSION,
        "snapshot_id": metadata["snapshot_id"],
        "inputs": {
            "metadata.json": file_sha256(metadata_path),
            alignment_path.name: file_sha256(alignment_path),
            timing_path.name: file_sha256(timing_path),
        },
        "plotted_rows": {"alignment": hit_counts, "execution_time": timing_counts},
        "outputs": {
            hit_path.name: file_sha256(hit_path),
            plot_data_path.name: file_sha256(plot_data_path),
            time_path.name: file_sha256(time_path),
        },
    }
    (output / "render_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    manifest = render(args.snapshot, args.output)
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
