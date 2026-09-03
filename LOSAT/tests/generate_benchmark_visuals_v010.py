#!/usr/bin/env python3
"""Generate fresh LOSAT-vs-NCBI diagnostic plots from the v0.1.0 manifests.

The plots produced here are diagnostics, not release certification authorities.
Every invocation runs the selected cases against the current LOSAT release binary
and the supplied official NCBI BLAST+ 2.17.0 binaries into a fresh output tree.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import statistics
import subprocess
import sys
import time
from typing import Any, Iterable

import audit_blastp_v010 as blastp
import audit_tblastx_v010 as tblastx
import compare_blastn_parity as blastn

NCBI_VERSION = "2.17.0"
PROGRAMS = ("blastn", "blastp", "tblastx")
TOOLS = ("NCBI BLAST+", "LOSAT")
REPRESENTATIVE_CASES = {
    "blastn": (
        "PesePMNV.MjPMNV.task_blastn",
        "Sakai.MG1655.megablast",
    ),
    "blastp": ("pairwise_default_serial",),
    "tblastx": (
        "p03_mela_pemojnva",
        "p14_ap027131_ap027133_query4",
        "d06_ap027131_ap027133_db4",
    ),
}
STANDARD_FIELDS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
)


@dataclass(frozen=True)
class CaseSpec:
    program: str
    case_id: str
    contract: str
    payload: Any


@dataclass(frozen=True)
class AlignmentRow:
    program: str
    case_id: str
    contract: str
    tool: str
    fields: tuple[str, ...]

    @property
    def pident(self) -> float:
        return float(self.fields[2])

    @property
    def length(self) -> int:
        return int(self.fields[3])


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        data_lines = (line for line in handle if not line.startswith("#"))
        return list(csv.DictReader(data_lines, delimiter="\t"))


def load_cases(repo_root: Path) -> list[CaseSpec]:
    crate_root = repo_root / "LOSAT"
    cases: list[CaseSpec] = []

    blastn_rows = read_tsv_manifest(crate_root / "tests" / "blastn_parity_manifest.tsv")
    for row in blastn_rows:
        case_id = row["case_id"].strip()
        contract = (
            "SOURCE_UNDETERMINED_ACCEPTED"
            if case_id == "Sakai.MG1655.megablast"
            else "parity"
        )
        cases.append(CaseSpec("blastn", case_id, contract, row))

    for case in blastp.load_manifest(
        crate_root / "tests" / "blastp_v010_parity_manifest.tsv", repo_root
    ):
        cases.append(CaseSpec("blastp", case.case_id, "parity", case))

    for case in tblastx.load_manifest(
        crate_root / "tests" / "tblastx_v010_parity_manifest.tsv", repo_root
    ):
        cases.append(CaseSpec("tblastx", case.case_id, case.contract, case))

    return cases


def select_cases(cases: Iterable[CaseSpec], scope: str) -> list[CaseSpec]:
    cases = list(cases)
    if scope == "full":
        return cases
    if scope != "representative":
        raise ValueError(f"unsupported scope: {scope}")

    by_key = {(case.program, case.case_id): case for case in cases}
    selected: list[CaseSpec] = []
    for program in PROGRAMS:
        for case_id in REPRESENTATIVE_CASES[program]:
            key = (program, case_id)
            if key not in by_key:
                raise ValueError(f"representative case missing from authority manifest: {key}")
            selected.append(by_key[key])
    return selected


def oracle_path(oracle_dir: Path, program: str) -> Path:
    path = oracle_dir / program
    if not path.is_file() or not os.access(path, os.X_OK):
        raise ValueError(f"NCBI oracle is missing or not executable: {path}")
    return path.resolve()


def validate_oracle(path: Path, program: str) -> dict[str, Any]:
    completed = subprocess.run(
        [str(path), "-version"], capture_output=True, text=True, check=True
    )
    version_text = completed.stdout.strip()
    expected = f"{program}: {NCBI_VERSION}+"
    if expected not in version_text:
        raise ValueError(f"unexpected NCBI oracle version for {program}: {version_text}")
    return {
        "path": str(path),
        "sha256": sha256_path(path),
        "version": version_text.splitlines(),
    }


def command_for(
    case: CaseSpec,
    tool: str,
    output_path: Path,
    repo_root: Path,
    losat_bin: Path,
    oracles: dict[str, Path],
) -> list[str]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    oracle = oracles[case.program]

    if case.program == "blastn":
        if tool == "NCBI BLAST+":
            return blastn.build_ncbi_command(case.payload, str(oracle), output_path)
        return blastn.build_losat_command(case.payload, losat_bin, output_path)

    if case.program == "blastp":
        ncbi_command, losat_command = blastp.build_commands(
            case.payload, oracle, losat_bin, output_path.parent
        )
        command = ncbi_command if tool == "NCBI BLAST+" else losat_command
        command = list(command)
        command[-1] = str(output_path)
        return command

    if case.program == "tblastx":
        scratch_ncbi = output_path.parent / f"{case.case_id}.ncbi.command.tsv"
        scratch_losat = output_path.parent / f"{case.case_id}.losat.command.tsv"
        ncbi_command, losat_command = tblastx.build_commands(
            case.payload, oracle, losat_bin, scratch_ncbi, scratch_losat
        )
        command = ncbi_command if tool == "NCBI BLAST+" else losat_command
        command = list(command)
        command[-1] = str(output_path)
        return command

    raise ValueError(f"unsupported program: {case.program}")


def command_environment(case: CaseSpec, tool: str) -> dict[str, str]:
    environment = os.environ.copy()
    if tool == "LOSAT" and case.program == "tblastx":
        environment["RAYON_NUM_THREADS"] = "1"
    return environment


def run_command(
    command: list[str],
    *,
    cwd: Path,
    environment: dict[str, str],
    stderr_path: Path,
) -> float:
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=cwd,
        env=environment,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    elapsed = time.perf_counter() - started
    stderr_path.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(command)}; "
            f"stderr={stderr_path}"
        )
    return elapsed


def parse_alignment_rows(
    path: Path, case: CaseSpec, tool: str
) -> list[AlignmentRow]:
    rows: list[AlignmentRow] = []
    with path.open(encoding="utf-8", errors="strict") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = tuple(line.split("\t"))
            if len(fields) < len(STANDARD_FIELDS):
                raise ValueError(
                    f"{path}:{line_number}: expected at least 12 tabular fields, "
                    f"got {len(fields)}"
                )
            rows.append(
                AlignmentRow(
                    program=case.program,
                    case_id=case.case_id,
                    contract=case.contract,
                    tool=tool,
                    fields=fields[: len(STANDARD_FIELDS)],
                )
            )
    return rows


def write_alignment_data(path: Path, rows: Iterable[AlignmentRow]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("program", "case_id", "contract", "tool", *STANDARD_FIELDS))
        for row in rows:
            writer.writerow((row.program, row.case_id, row.contract, row.tool, *row.fields))


def run_visual_snapshot(
    cases: list[CaseSpec],
    *,
    output_dir: Path,
    repo_root: Path,
    losat_bin: Path,
    oracles: dict[str, Path],
) -> tuple[list[AlignmentRow], list[dict[str, Any]]]:
    crate_root = repo_root / "LOSAT"
    raw_dir = output_dir / "raw" / "visual"
    all_rows: list[AlignmentRow] = []
    case_results: list[dict[str, Any]] = []

    for case in cases:
        outputs: dict[str, Path] = {}
        for tool in TOOLS:
            slug = "ncbi" if tool == "NCBI BLAST+" else "losat"
            output_path = raw_dir / case.program / case.case_id / f"{slug}.tsv"
            stderr_path = output_path.with_suffix(".stderr")
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.unlink(missing_ok=True)
            command = command_for(case, tool, output_path, repo_root, losat_bin, oracles)
            run_command(
                command,
                cwd=crate_root,
                environment=command_environment(case, tool),
                stderr_path=stderr_path,
            )
            outputs[tool] = output_path
            all_rows.extend(parse_alignment_rows(output_path, case, tool))

        ncbi_bytes = outputs["NCBI BLAST+"].read_bytes()
        losat_bytes = outputs["LOSAT"].read_bytes()
        case_results.append(
            {
                "program": case.program,
                "case_id": case.case_id,
                "contract": case.contract,
                "ncbi_sha256": hashlib.sha256(ncbi_bytes).hexdigest(),
                "losat_sha256": hashlib.sha256(losat_bytes).hexdigest(),
                "raw_byte_equal": ncbi_bytes == losat_bytes,
                "ncbi_rows": len(parse_alignment_rows(outputs["NCBI BLAST+"], case, "NCBI BLAST+")),
                "losat_rows": len(parse_alignment_rows(outputs["LOSAT"], case, "LOSAT")),
            }
        )

    return all_rows, case_results


def stable_sample(rows: list[AlignmentRow], limit: int, seed: str) -> list[AlignmentRow]:
    if len(rows) <= limit:
        return rows
    rng = random.Random(seed)
    return rng.sample(rows, limit)


def log_bins(values: list[int], count: int = 60) -> list[float]:
    positives = [value for value in values if value > 0]
    if not positives:
        return [1.0, 2.0]
    low = max(1.0, float(min(positives)))
    high = max(low * 1.01, float(max(positives)))
    low_log = math.log10(low)
    high_log = math.log10(high)
    step = (high_log - low_log) / count
    return [10 ** (low_log + index * step) for index in range(count + 1)]


def plot_hit_distribution(
    rows: list[AlignmentRow], path: Path, *, scope: str, candidate_sha: str
) -> None:
    import matplotlib.pyplot as plt

    labels = {"blastn": "BLASTN", "blastp": "BLASTP", "tblastx": "TBLASTX"}
    fig, axes = plt.subplots(len(PROGRAMS), 3, figsize=(18, 16))

    for row_index, program in enumerate(PROGRAMS):
        program_rows = [row for row in rows if row.program == program]
        for tool in TOOLS:
            tool_rows = [row for row in program_rows if row.tool == tool]
            if not tool_rows:
                continue
            lengths = [row.length for row in tool_rows]
            identities = [row.pident for row in tool_rows]
            axes[row_index][0].hist(
                lengths,
                bins=log_bins(lengths),
                weights=lengths,
                histtype="step",
                label=tool,
            )
            axes[row_index][1].hist(
                identities,
                bins=list(range(0, 102, 2)),
                weights=lengths,
                histtype="step",
                label=tool,
            )
            sampled = stable_sample(
                tool_rows, 20_000, f"{candidate_sha}:{program}:{tool}"
            )
            axes[row_index][2].scatter(
                [row.length for row in sampled],
                [row.pident for row in sampled],
                s=5,
                alpha=0.25,
                label=tool,
            )

        axes[row_index][0].set_xscale("log")
        axes[row_index][0].set_title(f"{labels[program]}: accumulated length vs alignment length")
        axes[row_index][0].set_xlabel("Alignment length (bp or aa)")
        axes[row_index][0].set_ylabel("Accumulated aligned length")
        axes[row_index][1].set_title(f"{labels[program]}: accumulated length vs identity")
        axes[row_index][1].set_xlabel("Identity (%)")
        axes[row_index][1].set_ylabel("Accumulated aligned length")
        axes[row_index][2].set_xscale("log")
        axes[row_index][2].set_title(f"{labels[program]}: length vs identity")
        axes[row_index][2].set_xlabel("Alignment length (bp or aa)")
        axes[row_index][2].set_ylabel("Identity (%)")
        for axis in axes[row_index]:
            axis.legend()
            axis.grid(alpha=0.2)

    fig.suptitle(
        f"Fresh LOSAT vs NCBI BLAST+ {NCBI_VERSION} hit distributions — {scope}\n"
        f"LOSAT candidate {candidate_sha[:12]}",
        fontsize=16,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def benchmark_case(
    case: CaseSpec,
    *,
    tool: str,
    output_dir: Path,
    repo_root: Path,
    losat_bin: Path,
    oracles: dict[str, Path],
    warmups: int,
    repetitions: int,
) -> dict[str, Any]:
    crate_root = repo_root / "LOSAT"
    samples: list[float] = []
    last_output: Path | None = None

    total_runs = warmups + repetitions
    for run_index in range(total_runs):
        phase = "warmup" if run_index < warmups else "sample"
        phase_index = run_index + 1 if phase == "warmup" else run_index - warmups + 1
        slug = "ncbi" if tool == "NCBI BLAST+" else "losat"
        output_path = (
            output_dir
            / "raw"
            / "performance"
            / case.program
            / case.case_id
            / f"{slug}.{phase}{phase_index}.tsv"
        )
        stderr_path = output_path.with_suffix(".stderr")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.unlink(missing_ok=True)
        command = command_for(case, tool, output_path, repo_root, losat_bin, oracles)
        elapsed = run_command(
            command,
            cwd=crate_root,
            environment=command_environment(case, tool),
            stderr_path=stderr_path,
        )
        last_output = output_path
        if phase == "sample":
            samples.append(elapsed)

    if last_output is None or not samples:
        raise RuntimeError(f"no timing samples produced for {case.program}/{case.case_id}/{tool}")

    return {
        "program": case.program,
        "case_id": case.case_id,
        "contract": case.contract,
        "tool": tool,
        "warmups": warmups,
        "repetitions": repetitions,
        "samples_seconds": samples,
        "median_seconds": statistics.median(samples),
        "last_output_sha256": sha256_path(last_output),
    }


def run_performance_snapshot(
    cases: list[CaseSpec],
    *,
    output_dir: Path,
    repo_root: Path,
    losat_bin: Path,
    oracles: dict[str, Path],
    warmups: int,
    repetitions: int,
) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    for case in cases:
        for tool in TOOLS:
            results.append(
                benchmark_case(
                    case,
                    tool=tool,
                    output_dir=output_dir,
                    repo_root=repo_root,
                    losat_bin=losat_bin,
                    oracles=oracles,
                    warmups=warmups,
                    repetitions=repetitions,
                )
            )
    return results


def write_performance_data(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "program",
                "case_id",
                "contract",
                "tool",
                "warmups",
                "repetitions",
                "median_seconds",
                "samples_seconds_json",
                "last_output_sha256",
            )
        )
        for row in rows:
            writer.writerow(
                (
                    row["program"],
                    row["case_id"],
                    row["contract"],
                    row["tool"],
                    row["warmups"],
                    row["repetitions"],
                    f"{row['median_seconds']:.9f}",
                    json.dumps(row["samples_seconds"], separators=(",", ":")),
                    row["last_output_sha256"],
                )
            )


def plot_execution_time(
    rows: list[dict[str, Any]], path: Path, *, scope: str, candidate_sha: str
) -> None:
    import matplotlib.pyplot as plt

    keys: list[tuple[str, str]] = []
    for row in rows:
        key = (row["program"], row["case_id"])
        if key not in keys:
            keys.append(key)
    labels = [f"{program}: {case_id}" for program, case_id in keys]
    positions = list(range(len(keys)))
    width = 0.36

    fig_height = max(6.0, 0.42 * len(keys) + 2.5)
    fig, axis = plt.subplots(figsize=(12, fig_height))
    values_by_tool = {
        tool: {
            (row["program"], row["case_id"]): row["median_seconds"]
            for row in rows
            if row["tool"] == tool
        }
        for tool in TOOLS
    }
    axis.barh(
        [position - width / 2 for position in positions],
        [values_by_tool["NCBI BLAST+"][key] for key in keys],
        height=width,
        label="NCBI BLAST+",
    )
    axis.barh(
        [position + width / 2 for position in positions],
        [values_by_tool["LOSAT"][key] for key in keys],
        height=width,
        label="LOSAT",
    )
    axis.set_yticks(positions, labels)
    axis.invert_yaxis()
    axis.set_xlabel("Median wall time (seconds)")
    axis.set_title(
        f"Diagnostic execution time — {scope} — LOSAT {candidate_sha[:12]} vs NCBI {NCBI_VERSION}"
    )
    axis.legend()
    axis.grid(axis="x", alpha=0.2)
    positive_values = [row["median_seconds"] for row in rows if row["median_seconds"] > 0]
    if positive_values and max(positive_values) / min(positive_values) >= 100:
        axis.set_xscale("log")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def write_case_results(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "program",
                "case_id",
                "contract",
                "raw_byte_equal",
                "ncbi_rows",
                "losat_rows",
                "ncbi_sha256",
                "losat_sha256",
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_summary(
    path: Path,
    *,
    scope: str,
    mode: str,
    candidate_sha: str,
    cases: list[CaseSpec],
    visual_case_results: list[dict[str, Any]],
    performance_rows: list[dict[str, Any]],
) -> None:
    lines = [
        "# Benchmark visualization snapshot",
        "",
        f"- Candidate SHA: `{candidate_sha}`",
        f"- NCBI BLAST+: `{NCBI_VERSION}`",
        f"- Scope: `{scope}`",
        f"- Mode: `{mode}`",
        f"- Selected cases: `{len(cases)}`",
        "- Authority: **diagnostic only; not a release certification gate**",
    ]
    if visual_case_results:
        exact = sum(1 for row in visual_case_results if row["raw_byte_equal"])
        lines.append(
            f"- Fresh raw-byte-equal visual cases: `{exact}/{len(visual_case_results)}` "
            "(known approved deviations remain visible rather than normalized)"
        )
    if performance_rows:
        lines.append(
            "- Timing: warmup + repeated same-runner wall-clock samples; medians plotted. "
            "GitHub-hosted timing remains diagnostic, not an official performance claim."
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv: list[str]) -> argparse.Namespace:
    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=repo_root)
    parser.add_argument(
        "--losat-bin",
        type=Path,
        default=repo_root / "LOSAT" / "target" / "release" / "LOSAT",
    )
    parser.add_argument("--oracle-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--scope", choices=("representative", "full"), default="representative"
    )
    parser.add_argument(
        "--mode", choices=("visual", "performance", "both"), default="visual"
    )
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--candidate-sha", default=os.environ.get("GITHUB_SHA", "unknown"))
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    repo_root = args.repo_root.resolve()
    crate_root = repo_root / "LOSAT"
    losat_bin = args.losat_bin.resolve()
    oracle_dir = args.oracle_dir.resolve()
    output_dir = args.output_dir.resolve()
    if not crate_root.is_dir():
        raise SystemExit(f"crate root missing: {crate_root}")
    if not losat_bin.is_file() or not os.access(losat_bin, os.X_OK):
        raise SystemExit(f"LOSAT release binary missing or not executable: {losat_bin}")
    if args.warmups < 0 or args.repetitions <= 0:
        raise SystemExit("warmups must be >= 0 and repetitions must be > 0")

    output_dir.mkdir(parents=True, exist_ok=True)
    all_cases = load_cases(repo_root)
    cases = select_cases(all_cases, args.scope)
    oracles = {program: oracle_path(oracle_dir, program) for program in PROGRAMS}
    oracle_identity = {
        program: validate_oracle(path, program) for program, path in oracles.items()
    }

    metadata: dict[str, Any] = {
        "schema": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "candidate_sha": args.candidate_sha,
        "scope": args.scope,
        "mode": args.mode,
        "ncbi_version": NCBI_VERSION,
        "oracle_identity": oracle_identity,
        "selected_cases": [
            {"program": case.program, "case_id": case.case_id, "contract": case.contract}
            for case in cases
        ],
        "runner": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "python": platform.python_version(),
            "github_runner_os": os.environ.get("RUNNER_OS"),
            "github_runner_arch": os.environ.get("RUNNER_ARCH"),
            "github_run_id": os.environ.get("GITHUB_RUN_ID"),
            "github_run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT"),
        },
        "performance_policy": {
            "warmups": args.warmups,
            "repetitions": args.repetitions,
            "summary_statistic": "median",
            "same_runner": True,
            "authority": "diagnostic_only",
        },
    }

    visual_case_results: list[dict[str, Any]] = []
    performance_rows: list[dict[str, Any]] = []

    if args.mode in {"visual", "both"}:
        alignment_rows, visual_case_results = run_visual_snapshot(
            cases,
            output_dir=output_dir,
            repo_root=repo_root,
            losat_bin=losat_bin,
            oracles=oracles,
        )
        write_alignment_data(output_dir / "hit_distribution_data.tsv", alignment_rows)
        write_case_results(output_dir / "visual_case_results.tsv", visual_case_results)
        plot_hit_distribution(
            alignment_rows,
            output_dir / "hit_distribution.png",
            scope=args.scope,
            candidate_sha=args.candidate_sha,
        )
        metadata["visual_case_results"] = visual_case_results

    if args.mode in {"performance", "both"}:
        performance_rows = run_performance_snapshot(
            cases,
            output_dir=output_dir,
            repo_root=repo_root,
            losat_bin=losat_bin,
            oracles=oracles,
            warmups=args.warmups,
            repetitions=args.repetitions,
        )
        write_performance_data(output_dir / "performance_data.tsv", performance_rows)
        plot_execution_time(
            performance_rows,
            output_dir / "execution_time.png",
            scope=args.scope,
            candidate_sha=args.candidate_sha,
        )
        metadata["performance_results"] = performance_rows

    (output_dir / "benchmark_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_summary(
        output_dir / "summary.md",
        scope=args.scope,
        mode=args.mode,
        candidate_sha=args.candidate_sha,
        cases=cases,
        visual_case_results=visual_case_results,
        performance_rows=performance_rows,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
