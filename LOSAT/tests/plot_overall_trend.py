#!/usr/bin/env python
# coding: utf-8

import os
import re

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# === Configuration ===
PLOT_DIR = "./plots"
OUTPUT_IMAGE = "./plots/overall_trend_comparison.png"
os.makedirs(PLOT_DIR, exist_ok=True)
# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:75
# ```c
# const string kArgNumThreads("num_threads");
# ```
LOSAT_THREADS = (
    os.environ.get("LOSAT_THREADS")
    or os.environ.get("LOSATP_THREADS")
    or os.environ.get("LOSAT_BLASTP_THREADS")
    or "8"
)
LOSAT_WASM_THREADS_VERIFIED = os.environ.get("LOSAT_WASM_THREADS_VERIFIED") == "1"
LOSAT_WASM_SINGLE_LABEL = (
    "LOSAT wasm n1" if LOSAT_WASM_THREADS_VERIFIED else "LOSAT wasm serial"
)
LOSAT_WASM_MULTI_LABEL = (
    f"LOSAT wasm n{LOSAT_THREADS}"
    if LOSAT_WASM_THREADS_VERIFIED
    else f"LOSAT wasm serial (requested n{LOSAT_THREADS})"
)

# === Color Settings (Seaborn Deep Palette) ===
CUSTOM_PALETTE = {
    "LOSAT": "#dd8452",
    "BLAST+": "#4c72b0",
    LOSAT_WASM_SINGLE_LABEL: "#8172b3",
    LOSAT_WASM_MULTI_LABEL: "#937860",
}
HUE_ORDER = ["BLAST+", "LOSAT", LOSAT_WASM_SINGLE_LABEL, LOSAT_WASM_MULTI_LABEL]

# === Column Definitions ===
COLUMNS = [
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
]

TBLASTX_COMPARISONS = [
    {"name": "NZ_CP006932_Self", "mode": "TBLASTX", "group": "TBLASTX", "ncbi": "./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n8.out", "losat": "./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n8.out"},
    {"name": "AP027132_vs_NZ", "mode": "TBLASTX", "group": "TBLASTX", "ncbi": "./blast_out/AP027132.NZ_CP006932.tblastx.n8.out", "losat": "./losat_out/AP027132.NZ_CP006932.tlosatx.n8.out"},
    {"name": "AP027078_vs_AP027131", "mode": "TBLASTX", "group": "TBLASTX", "ncbi": "./blast_out/AP027078.AP027131.tblastx.n8.out", "losat": "./losat_out/AP027078.AP027131.tlosatx.n8.out"},
    {"name": "AP027131_vs_AP027133", "mode": "TBLASTX", "group": "TBLASTX", "ncbi": "./blast_out/AP027131.AP027133.tblastx.n8.out", "losat": "./losat_out/AP027131.AP027133.tlosatx.n8.out"},
    {"name": "AP027133_vs_AP027132", "mode": "TBLASTX", "group": "TBLASTX", "ncbi": "./blast_out/AP027133.AP027132.tblastx.n8.out", "losat": "./losat_out/AP027133.AP027132.tlosatx.n8.out"},
    {"name": "AP027280_Self", "mode": "TBLASTX", "group": "TBLASTX", "ncbi": "./blast_out/AP027280.AP027280.tblastx.n8.out", "losat": "./losat_out/AP027280.AP027280.tlosatx.n8.out"},
    {"name": "MjeNMV_vs_MelaMJNV", "mode": "TBLASTX", "group": "TBLASTX", "ncbi": "./blast_out/MjeNMV.MelaMJNV.tblastx.n8.out", "losat": "./losat_out/MjeNMV.MelaMJNV.tlosatx.n8.out"},
    {"name": "MelaMJNV_vs_PemoMJNVA", "mode": "TBLASTX", "group": "TBLASTX", "ncbi": "./blast_out/MelaMJNV.PemoMJNVA.tblastx.n8.out", "losat": "./losat_out/MelaMJNV.PemoMJNVA.tlosatx.n8.out"},
    {"name": "PemoMJNVA_vs_PeseMJNV", "mode": "TBLASTX", "group": "TBLASTX", "ncbi": "./blast_out/PemoMJNVA.PeseMJNV.tblastx.n8.out", "losat": "./losat_out/PemoMJNVA.PeseMJNV.tlosatx.n8.out"},
    {"name": "PeseMJNV_vs_PemoMJNVB", "mode": "TBLASTX", "group": "TBLASTX", "ncbi": "./blast_out/PeseMJNV.PemoMJNVB.tblastx.n8.out", "losat": "./losat_out/PeseMJNV.PemoMJNVB.tlosatx.n8.out"},
    {"name": "PemoMJNVB_vs_LvMJNV", "mode": "TBLASTX", "group": "TBLASTX", "ncbi": "./blast_out/PemoMJNVB.LvMJNV.tblastx.n8.out", "losat": "./losat_out/PemoMJNVB.LvMJNV.tlosatx.n8.out"},
    {"name": "LvMJNV_vs_TrcuMJNV", "mode": "TBLASTX", "group": "TBLASTX", "ncbi": "./blast_out/LvMJNV.TrcuMJNV.tblastx.n8.out", "losat": "./losat_out/LvMJNV.TrcuMJNV.tlosatx.n8.out"},
    {"name": "TrcuMJNV_vs_MellatMJNV", "mode": "TBLASTX", "group": "TBLASTX", "ncbi": "./blast_out/TrcuMJNV.MellatMJNV.tblastx.n8.out", "losat": "./losat_out/TrcuMJNV.MellatMJNV.tlosatx.n8.out"},
    {"name": "MellatMJNV_vs_MeenMJNV", "mode": "TBLASTX", "group": "TBLASTX", "ncbi": "./blast_out/MellatMJNV.MeenMJNV.tblastx.n8.out", "losat": "./losat_out/MellatMJNV.MeenMJNV.tlosatx.n8.out"},
    {"name": "MeenMJNV_vs_MejoMJNV", "mode": "TBLASTX", "group": "TBLASTX", "ncbi": "./blast_out/MeenMJNV.MejoMJNV.tblastx.n8.out", "losat": "./losat_out/MeenMJNV.MejoMJNV.tlosatx.n8.out"},
    {"name": "AvCLPV_vs_PsCLPV", "mode": "TBLASTX", "group": "TBLASTX", "ncbi": "./blast_out/AvCLPV.PsCLPV.tblastx.n8.out", "losat": "./losat_out/AvCLPV.PsCLPV.tlosatx.n8.out"},
]

BLASTN_COMPARISONS = [
    {"name": "NZ_CP006932_Self(Mega)", "mode": "BLASTN (Megablast)", "group": "BLASTN (All Types)", "ncbi": "./blast_out/NZ_CP006932.NZ_CP006932.blastn.out", "losat": "./losat_out/NZ_CP006932.NZ_CP006932.losatn.megablast.out"},
    {"name": "EDL933_vs_Sakai", "mode": "BLASTN (Megablast)", "group": "BLASTN (All Types)", "ncbi": "./blast_out/EDL933.Sakai.blastn.megablast.out", "losat": "./losat_out/EDL933.Sakai.losatn.megablast.out"},
    {"name": "Sakai_vs_MG1655", "mode": "BLASTN (Megablast)", "group": "BLASTN (All Types)", "ncbi": "./blast_out/Sakai.MG1655.blastn.megablast.out", "losat": "./losat_out/Sakai.MG1655.losatn.megablast.out"},
    {"name": "NZ_CP006932_Self(Task)", "mode": "BLASTN (Task:blastn)", "group": "BLASTN (All Types)", "ncbi": "./blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out", "losat": "./losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.out"},
    {"name": "PesePMNV_vs_MjPMNV", "mode": "BLASTN (Task:blastn)", "group": "BLASTN (All Types)", "ncbi": "./blast_out/PesePMNV.MjPMNV.blastn.out", "losat": "./losat_out/PesePMNV.MjPMNV.losatn.blastn.out"},
    {"name": "MelaMJNV_vs_PemoMJNVA", "mode": "BLASTN (Task:blastn)", "group": "BLASTN (All Types)", "ncbi": "./blast_out/MelaMJNV.PemoMJNVA.blastn.out", "losat": "./losat_out/MelaMJNV.PemoMJNVA.losatn.blastn.out"},
    {"name": "SiNMV_vs_ChdeNMV", "mode": "BLASTN (Task:blastn)", "group": "BLASTN (All Types)", "ncbi": "./blast_out/SiNMV.ChdeNMV.blastn.out", "losat": "./losat_out/SiNMV.ChdeNMV.losatn.blastn.out"},
    {"name": "PmeNMV_vs_MjPMNV", "mode": "BLASTN (Task:blastn)", "group": "BLASTN (All Types)", "ncbi": "./blast_out/PmeNMV.MjPMNV.blastn.out", "losat": "./losat_out/PmeNMV.MjPMNV.losatn.blastn.out"},
    {"name": "PmeNMV_vs_PesePMNV", "mode": "BLASTN (Task:blastn)", "group": "BLASTN (All Types)", "ncbi": "./blast_out/PmeNMV.PesePMNV.blastn.out", "losat": "./losat_out/PmeNMV.PesePMNV.losatn.blastn.out"},
    {"name": "PeseMJNV_vs_PemoMJNVB", "mode": "BLASTN (Task:blastn)", "group": "BLASTN (All Types)", "ncbi": "./blast_out/PeseMJNV.PemoMJNVB.blastn.out", "losat": "./losat_out/PeseMJNV.PemoMJNVB.losatn.blastn.out"},
    {"name": "PemoMJNVA_vs_PeseMJNV", "mode": "BLASTN (Task:blastn)", "group": "BLASTN (All Types)", "ncbi": "./blast_out/PemoMJNVA.PeseMJNV.blastn.out", "losat": "./losat_out/PemoMJNVA.PeseMJNV.losatn.blastn.out"},
    {"name": "MjeNMV_vs_MelaMJNV", "mode": "BLASTN (Task:blastn)", "group": "BLASTN (All Types)", "ncbi": "./blast_out/MjeNMV.MelaMJNV.blastn.out", "losat": "./losat_out/MjeNMV.MelaMJNV.losatn.blastn.out"},
    {"name": "MjPMNV_vs_MlPMNV", "mode": "BLASTN (Task:blastn)", "group": "BLASTN (All Types)", "ncbi": "./blast_out/MjPMNV.MlPMNV.blastn.out", "losat": "./losat_out/MjPMNV.MlPMNV.losatn.blastn.out"},
]

BLASTP_COMPARISONS = [
    {"name": "WSSV.PajaWSV", "mode": "BLASTP", "group": "BLASTP", "ncbi": "./blast_out/WSSV.PajaWSV.BLASTP.out", "losat": "./losat_out/WSSV.PajaWSV.losatp.out"},
    {"name": "WSSV.SicyWSV", "mode": "BLASTP", "group": "BLASTP", "ncbi": "./blast_out/WSSV.SicyWSV.BLASTP.out", "losat": "./losat_out/WSSV.SicyWSV.losatp.out"},
    {"name": "WSSV.CoBV", "mode": "BLASTP", "group": "BLASTP", "ncbi": "./blast_out/WSSV.CoBV.BLASTP.out", "losat": "./losat_out/WSSV.CoBV.losatp.out"},
    {"name": "PajaWSV.SicyWSV", "mode": "BLASTP", "group": "BLASTP", "ncbi": "./blast_out/PajaWSV.SicyWSV.BLASTP.out", "losat": "./losat_out/PajaWSV.SicyWSV.losatp.out"},
    {"name": "PajaWSV.CoBV", "mode": "BLASTP", "group": "BLASTP", "ncbi": "./blast_out/PajaWSV.CoBV.BLASTP.out", "losat": "./losat_out/PajaWSV.CoBV.losatp.out"},
    {"name": "SicyWSV.CoBV", "mode": "BLASTP", "group": "BLASTP", "ncbi": "./blast_out/SicyWSV.CoBV.BLASTP.out", "losat": "./losat_out/SicyWSV.CoBV.losatp.out"},
    {"name": "AP027078.AP027131", "mode": "BLASTP", "group": "BLASTP", "ncbi": "./blast_out/AP027078.AP027131.BLASTP.out", "losat": "./losat_out/AP027078.AP027131.losatp.out"},
    {"name": "AP027078.AP027132", "mode": "BLASTP", "group": "BLASTP", "ncbi": "./blast_out/AP027078.AP027132.BLASTP.out", "losat": "./losat_out/AP027078.AP027132.losatp.out"},
    {"name": "AP027078.AP027133", "mode": "BLASTP", "group": "BLASTP", "ncbi": "./blast_out/AP027078.AP027133.BLASTP.out", "losat": "./losat_out/AP027078.AP027133.losatp.out"},
    {"name": "AP027078.NZ_CP006932", "mode": "BLASTP", "group": "BLASTP", "ncbi": "./blast_out/AP027078.NZ_CP006932.BLASTP.out", "losat": "./losat_out/AP027078.NZ_CP006932.losatp.out"},
    {"name": "AP027131.AP027132", "mode": "BLASTP", "group": "BLASTP", "ncbi": "./blast_out/AP027131.AP027132.BLASTP.out", "losat": "./losat_out/AP027131.AP027132.losatp.out"},
    {"name": "AP027131.AP027133", "mode": "BLASTP", "group": "BLASTP", "ncbi": "./blast_out/AP027131.AP027133.BLASTP.out", "losat": "./losat_out/AP027131.AP027133.losatp.out"},
    {"name": "AP027131.NZ_CP006932", "mode": "BLASTP", "group": "BLASTP", "ncbi": "./blast_out/AP027131.NZ_CP006932.BLASTP.out", "losat": "./losat_out/AP027131.NZ_CP006932.losatp.out"},
    {"name": "AP027132.AP027133", "mode": "BLASTP", "group": "BLASTP", "ncbi": "./blast_out/AP027132.AP027133.BLASTP.out", "losat": "./losat_out/AP027132.AP027133.losatp.out"},
    {"name": "AP027132.NZ_CP006932", "mode": "BLASTP", "group": "BLASTP", "ncbi": "./blast_out/AP027132.NZ_CP006932.BLASTP.out", "losat": "./losat_out/AP027132.NZ_CP006932.losatp.out"},
]

comparisons = TBLASTX_COMPARISONS + BLASTN_COMPARISONS + BLASTP_COMPARISONS


def strip_thread_suffix(filepath):
    return re.sub(r"\.n\d+\.out$", ".out", filepath)


def add_wasm_suffix(filepath, suffix=None):
    if not filepath.endswith(".out"):
        return f"{filepath}.wasm{'.' + suffix if suffix else ''}.out"
    stem = filepath[:-4]
    if suffix:
        return f"{stem}.wasm.{suffix}.out"
    return f"{stem}.wasm.out"


def losat_wasm_outputs(losat_path):
    threaded_suffix = f"n{LOSAT_THREADS}"
    base_path = strip_thread_suffix(losat_path)

    return [
        (LOSAT_WASM_SINGLE_LABEL, add_wasm_suffix(base_path)),
        (LOSAT_WASM_MULTI_LABEL, add_wasm_suffix(base_path, threaded_suffix)),
    ]


def wasm_log_path_from_output(filepath):
    return re.sub(r"\.out$", ".log", filepath)


def parse_effective_engine_threads(filepath):
    if not os.path.exists(filepath):
        return None
    try:
        with open(filepath, "r") as f:
            content = f.read()
    except Exception:
        return None
    match = re.search(r"effective_engine_threads=(\d+)", content)
    if match:
        return int(match.group(1))
    match = re.search(r"\[TIMING\] engine_threads: requested=\S+ effective=(\d+)", content)
    if match:
        return int(match.group(1))
    return None


def warn_if_serial_wasm_thread_label(tool_name, output_path):
    if re.search(r"\bwasm n\d+\b", tool_name, flags=re.IGNORECASE):
        log_path = wasm_log_path_from_output(output_path)
        effective_threads = parse_effective_engine_threads(log_path)
        if effective_threads == 1:
            print(
                f"[Warning] {tool_name} label was produced from serial Wasm log: {log_path}"
            )


def load_blast(filepath, tool_name, mode_name, task_name, group_name):
    """Load a BLAST file and add metadata columns."""
    try:
        df = pd.read_csv(filepath, sep="\t", comment="#", names=COLUMNS, header=None)
        numeric_cols = ["pident", "length", "bitscore", "evalue"]
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors="coerce")

        df["Tool"] = tool_name
        df["Mode"] = mode_name
        df["Task"] = task_name
        df["Broad_Mode"] = group_name
        return df
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None


def main():
    print("Loading all datasets...")
    all_data = []

    for item in comparisons:
        missing = []
        if not os.path.exists(item["ncbi"]):
            missing.append(f"BLAST+({item['ncbi']})")
        if not os.path.exists(item["losat"]):
            missing.append(f"LOSAT({item['losat']})")
        if missing:
            print(f"Skipping {item['name']} ({item['mode']}): missing {', '.join(missing)}")
            continue

        df_ncbi = load_blast(item["ncbi"], "BLAST+", item["mode"], item["name"], item["group"])
        if df_ncbi is not None:
            all_data.append(df_ncbi)

        df_losat = load_blast(item["losat"], "LOSAT", item["mode"], item["name"], item["group"])
        if df_losat is not None:
            all_data.append(df_losat)

        for tool_name, wasm_path in losat_wasm_outputs(item["losat"]):
            if not os.path.exists(wasm_path):
                print(f"Optional Wasm output missing for {item['name']}: {tool_name}({wasm_path})")
                continue
            warn_if_serial_wasm_thread_label(tool_name, wasm_path)
            df_wasm = load_blast(wasm_path, tool_name, item["mode"], item["name"], item["group"])
            if df_wasm is not None:
                all_data.append(df_wasm)

    if not all_data:
        print("No data found.")
        return

    df_all = pd.concat(all_data, ignore_index=True)
    print(f"Total alignment records loaded: {len(df_all)}")

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
            plot_data = data_subset.sample(100000, random_state=42)
            title_suffix = "(Subsampled 100k)"
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

    plt.suptitle("Overall Trends: LOSAT vs BLAST+ (Paired Tasks)", fontsize=20, y=1.02)
    plt.tight_layout()
    plt.savefig(OUTPUT_IMAGE, bbox_inches="tight")
    print(f"Overall trend plot saved to {OUTPUT_IMAGE}")


if __name__ == "__main__":
    main()
