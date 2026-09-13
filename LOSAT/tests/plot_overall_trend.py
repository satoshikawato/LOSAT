#!/usr/bin/env python3
"""Aggregate the same pairs across every available native/Wasm series."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# NCBI reference: c++/src/objtools/align_format/format_flags.cpp:38-40
# "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
# "evalue bitscore";
# Compare these unchanged fields on a common cohort of query/subject pairs.
from comparison_data import CUSTOM_PALETTE, HUE_ORDER, PLOT_DIR, comparison_cases, load_case

OUTPUT_IMAGE = PLOT_DIR / "overall_trend_comparison.png"


def main():
    groups = {}
    for case in comparison_cases():
        frames = load_case(case)
        if frames:
            groups.setdefault(case["group"], []).append((case, frames))
    all_data = []
    for group, pairs in groups.items():
        tools = set().union(*(frames.keys() for _, frames in pairs))
        kept = 0
        for case, frames in pairs:
            missing = tools - frames.keys()
            if missing:
                print(f"[Unpaired] {case['name']} ({case['mode']}): missing {', '.join(sorted(missing))}")
                continue
            all_data.extend(frames.values())
            kept += 1
        print(f"{group}: {kept}/{len(pairs)} common pairs across {len(tools)} series")
    if not all_data:
        print("No common pairs to plot.")
        return
    df_all = pd.concat(all_data, ignore_index=True)
    if df_all.empty:
        print("All paired outputs have zero hits.")
        return
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    sns.set(style="whitegrid")

    mode_order = ["TBLASTX", "BLASTN (All Types)", "BLASTP"]
    modes = [mode for mode in mode_order if mode in df_all["Broad_Mode"].unique()]
    fig, axes = plt.subplots(len(modes), 3, figsize=(18, 6 * len(modes)))

    if len(modes) == 1:
        axes = [axes]

    for i, mode in enumerate(modes):
        data_subset = df_all[df_all["Broad_Mode"] == mode]

        axes[i][0].text(
            -0.2,
            0.5,
            mode,
            transform=axes[i][0].transAxes,
            fontsize=16,
            rotation=90,
            va="center",
            fontweight="bold",
        )

        sns.histplot(
            data=data_subset,
            x="length",
            hue="Tool",
            weights="length",
            bins=100,
            element="step",
            stat="count",
            common_norm=False,
            log_scale=True,
            ax=axes[i][0],
            palette=CUSTOM_PALETTE,
            hue_order=[tool for tool in HUE_ORDER if tool in data_subset["Tool"].unique()],
        )
        axes[i][0].set_title("Accumulated Length vs Alignment Length")
        axes[i][0].set_xlabel("Length (bp or aa)")
        axes[i][0].set_ylabel("Accumulated Length (bp/aa)")

        sns.histplot(
            data=data_subset,
            x="pident",
            hue="Tool",
            weights="length",
            bins=100,
            element="step",
            stat="count",
            common_norm=False,
            ax=axes[i][1],
            palette=CUSTOM_PALETTE,
            hue_order=[tool for tool in HUE_ORDER if tool in data_subset["Tool"].unique()],
        )
        axes[i][1].set_title("Accumulated Length vs Identity")
        axes[i][1].set_xlabel("Identity (%)")
        axes[i][1].set_ylabel("Accumulated Length (bp/aa)")

        if len(data_subset) > 100000:
            # Identical ordered outputs must select the same scatter points,
            # even when adding Wasm pushes the total above the display limit.
            per_tool = 100000 // data_subset["Tool"].nunique()
            plot_data = pd.concat([
                frame.sample(min(len(frame), per_tool), random_state=42)
                for _, frame in data_subset.groupby("Tool", sort=False)
            ])
            title_suffix = f"(At most {per_tool:,} points/series)"
        else:
            plot_data = data_subset
            title_suffix = ""

        sns.scatterplot(
            data=plot_data,
            x="length",
            y="pident",
            hue="Tool",
            style="Tool",
            alpha=0.5,
            ax=axes[i][2],
            palette=CUSTOM_PALETTE,
            hue_order=[tool for tool in HUE_ORDER if tool in plot_data["Tool"].unique()],
            style_order=[tool for tool in HUE_ORDER if tool in plot_data["Tool"].unique()],
        )
        axes[i][2].set_xscale("log")
        axes[i][2].set_title(f"Length vs Identity {title_suffix}")
        axes[i][2].set_xlabel("Length (bp or aa)")
        axes[i][2].set_ylabel("Identity (%)")

    plt.suptitle("Overall Distributions: BLAST+ vs LOSAT Native/Wasm (Common Pairs)", fontsize=20, y=1.02)
    plt.tight_layout()
    plt.savefig(OUTPUT_IMAGE, bbox_inches="tight")
    print(f"Overall trend plot saved to {OUTPUT_IMAGE}")
    plt.close(fig)


if __name__ == "__main__":
    main()
