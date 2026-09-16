#!/usr/bin/env python3
"""Plot each historical pair, with native and command-Wasm results."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# NCBI reference: c++/src/objtools/align_format/format_flags.cpp:38-40
# "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
# "evalue bitscore";
# Plot those fields without changing or filtering the alignments.
from comparison_data import CUSTOM_PALETTE, HUE_ORDER, PLOT_DIR, comparison_cases, load_case, require_plot_run


def generate_comparison_plot(config):
    name, mode_label = config["name"], config["mode"]
    frames = load_case(config)
    if not frames:
        return
    df_merged = pd.concat(frames.values(), ignore_index=True)
    if df_merged.empty:
        print(f"[Empty] {name} ({mode_label}): all available outputs have zero hits")
        return
    hue_order = [tool for tool in HUE_ORDER if tool in frames and not frames[tool].empty]
    print(f"{name} ({mode_label}): " + ", ".join(f"{tool}: {len(df)} hits" for tool, df in frames.items()))
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    sns.set(style="whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f"Comparison: {name} ({mode_label})", fontsize=16)

    sns.histplot(
        data=df_merged,
        x="length",
        hue="Tool",
        weights="length",
        bins=100,
        element="step",
        stat="count",
        common_norm=False,
        log_scale=True,
        ax=axes[0, 0],
        palette=CUSTOM_PALETTE,
        hue_order=hue_order,
    )
    axes[0, 0].set_title("Accumulated Length vs Alignment Length")
    axes[0, 0].set_xlabel("Alignment Length (bp/aa)")
    axes[0, 0].set_ylabel("Accumulated Length (bp/aa)")

    sns.histplot(
        data=df_merged,
        x="pident",
        hue="Tool",
        weights="length",
        bins=100,
        element="step",
        stat="count",
        common_norm=False,
        ax=axes[0, 1],
        palette=CUSTOM_PALETTE,
        hue_order=hue_order,
    )
    axes[0, 1].set_title("Accumulated Length vs Identity")
    axes[0, 1].set_xlabel("Identity (%)")
    axes[0, 1].set_ylabel("Accumulated Length (bp/aa)")

    sns.scatterplot(
        data=df_merged,
        x="length",
        y="pident",
        hue="Tool",
        alpha=0.5,
        style="Tool",
        ax=axes[1, 0],
        palette=CUSTOM_PALETTE,
        hue_order=hue_order,
        style_order=hue_order,
    )
    axes[1, 0].set_xscale("log")
    axes[1, 0].set_title("Alignment Length vs Identity")
    axes[1, 0].set_xlabel("Alignment Length (bp/aa)")
    axes[1, 0].set_ylabel("Identity (%)")

    sns.scatterplot(
        data=df_merged,
        x="length",
        y="bitscore",
        hue="Tool",
        style="Tool",
        alpha=0.5,
        ax=axes[1, 1],
        palette=CUSTOM_PALETTE,
        hue_order=hue_order,
        style_order=hue_order,
    )
    axes[1, 1].set_xscale("log")
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_title("Alignment Length vs Bit Score")
    axes[1, 1].set_xlabel("Alignment Length (bp/aa)")
    axes[1, 1].set_ylabel("Bit Score")

    output_filename = os.path.join(PLOT_DIR, f"compare_{name}_{mode_label}.png")
    plt.tight_layout()
    plt.savefig(output_filename)
    plt.close()
    print(f"  -> Saved to {output_filename}")


def main():
    for case in comparison_cases():
        generate_comparison_plot(case)


if __name__ == "__main__":
    # NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
    # CATCH_ALL(status) ... return status;
    require_plot_run()
    main()
