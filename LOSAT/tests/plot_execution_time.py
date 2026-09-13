#!/usr/bin/env python3
"""Render the simple Bash time logs, including serial and threaded Wasm."""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:75
# const string kArgNumThreads("num_threads");
# Thread suffixes describe requested counts, not inferred worker utilization.
from comparison_data import (
    CUSTOM_PALETTE, HUE_ORDER, LOSAT_THREADS, MODE_ORDER, PLOT_DIR,
    NATIVE_SINGLE, NATIVE_MULTI, WASM_SINGLE, WASM_MULTI,
    comparison_cases, completed_log, result_paths, successful_output, sha256,
)

OUTPUT_IMAGE = PLOT_DIR / "execution_time_comparison_all.png"


def parse_time(filepath):
    """Read Bash or GNU wall time; explicitly failed runs are never timings."""
    filepath = Path(filepath)
    if not filepath.is_file():
        return None
    content = filepath.read_text()
    if not completed_log(content):
        return None
    # The harness appends Bash time after command diagnostics; a command can
    # itself print a 'real' line, so use the final measurement.
    matches = list(re.finditer(r"^real\s+(?:(\d+)m)?(\d+(?:\.\d+)?)s?\s*$", content, re.M))
    if matches:
        match = matches[-1]
        return float(match[1] or 0) * 60 + float(match[2])
    match = re.search(r"(?:(\d+):)?(\d+):(\d+(?:\.\d+)?)elapsed", content)
    if match:
        return float(match[1] or 0) * 3600 + float(match[2]) * 60 + float(match[3])
    return None


def main():
    data = []
    # NCBI reference: c++/src/objtools/align_format/tabular.cpp:1100-1108
    # x_PrintField(*iter); ... m_Ostream << "\n";
    # Performance bars require raw equality, before any numeric plot parsing.
    for case in comparison_cases():
        reference = result_paths(case)["BLAST+"]
        if not successful_output(reference):
            continue
        expected_hash = sha256(reference)
        for tool, log in result_paths(case, "log").items():
            if not successful_output(log.with_suffix(".out")):
                continue
            if sha256(log.with_suffix(".out")) != expected_hash:
                print(f"[Raw mismatch; excluded from timing] {log}")
                continue
            seconds = parse_time(log)
            if seconds is None:
                print(f"[Missing/failed timing] {log}")
                continue
            data.append({
                "Task": case["name"], "Mode": case["mode"], "Tool": tool,
                "Time (s)": seconds, "Query": case["query"],
                "Subject": case["subject"], "Log": str(log),
            })
    if not data:
        print("No valid time data found.")
        return
    df = pd.DataFrame(data)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(PLOT_DIR / "execution_times.tsv", sep="\t", index=False)
    sns.set(style="whitegrid")
    modes = [mode for mode in MODE_ORDER if mode in df["Mode"].unique()]
    pair_count = df.groupby("Mode")["Task"].nunique().max()
    height = max(3.5, min(10, 1.5 + 0.45 * pair_count))
    g = sns.catplot(
        data=df, kind="bar", y="Task", x="Time (s)", hue="Tool", col="Mode",
        height=height, aspect=7 / height, sharex=False, sharey=False,
        palette=CUSTOM_PALETTE,
        hue_order=[tool for tool in HUE_ORDER if tool in df["Tool"].unique()],
        errorbar=None, col_wrap=min(2, len(modes)), col_order=modes,
    )
    for ax in g.axes.flat:
        ax.set_xlim(left=0)
    g.despine(left=True)
    g.set_axis_labels("Wall time (seconds)", "")
    g.fig.suptitle(
        "Execution Time: BLAST+ vs LOSAT Native/Wasm\n"
        "One run per condition; includes startup and Wasm compilation\n"
        f"BLAST+ requested threads: TBLASTX n{LOSAT_THREADS}, other tasks n1",
        y=1.12, fontsize=12,
    )
    g.fig.savefig(OUTPUT_IMAGE, bbox_inches="tight")
    plt.close(g.fig)
    print(f"Plot saved to {OUTPUT_IMAGE}")
    summary = df.pivot(index=["Mode", "Task"], columns="Tool", values="Time (s)")
    ratios = [(WASM_SINGLE, NATIVE_SINGLE), (WASM_MULTI, NATIVE_MULTI)]
    for numerator, denominator in ratios:
        if {numerator, denominator}.issubset(summary.columns):
            summary[f"Ratio ({numerator}/{denominator})"] = (
                summary[numerator] / summary[denominator].where(summary[denominator] > 0)
            ).round(2)
    print(summary.to_string())


if __name__ == "__main__":
    main()
