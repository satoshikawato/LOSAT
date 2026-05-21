#!/usr/bin/env python
# coding: utf-8

import os
import re

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


# === Configuration: Output Directory ===
PLOT_DIR = "./plots"
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
LOSAT_WASM_EXECUTION_MODE = os.environ.get(
    "LOSAT_WASM_EXECUTION_MODE", "command-wasi-serial"
)
LOSAT_WASM_SIMD_LABEL = "LOSAT Wasm SIMD"
LOSAT_WASM_SCALAR_LABEL = "LOSAT Wasm scalar"


# NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1407-1427
# ```c
# db_length = BlastSeqSrcGetTotLen(seq_src);
# itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
# while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) {
#     if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
#         continue;
#     }
# }
# ```
def wasm_mode_labels():
    if LOSAT_WASM_EXECUTION_MODE == "command-wasi-serial":
        return (
            "LOSAT command-WASI serial",
            f"LOSAT command-WASI serial (requested n{LOSAT_THREADS})",
        )
    if LOSAT_WASM_EXECUTION_MODE == "browser-in-memory-serial":
        return (
            "LOSAT browser in-memory serial",
            f"LOSAT browser in-memory serial (requested n{LOSAT_THREADS})",
        )
    if LOSAT_WASM_EXECUTION_MODE == "browser-worker-parallel":
        if LOSAT_WASM_THREADS_VERIFIED:
            return ("LOSAT browser worker n1", f"LOSAT browser worker n{LOSAT_THREADS}")
        return (
            "LOSAT browser worker unverified",
            f"LOSAT browser worker unverified (requested n{LOSAT_THREADS})",
        )
    if LOSAT_WASM_EXECUTION_MODE == "future-wasi-threaded":
        if LOSAT_WASM_THREADS_VERIFIED:
            return ("LOSAT WASI threaded n1", f"LOSAT WASI threaded n{LOSAT_THREADS}")
        return (
            "LOSAT WASI threaded unverified",
            f"LOSAT WASI threaded unverified (requested n{LOSAT_THREADS})",
        )
    if LOSAT_WASM_THREADS_VERIFIED:
        return ("LOSAT wasm n1", f"LOSAT wasm n{LOSAT_THREADS}")
    return ("LOSAT wasm serial", f"LOSAT wasm serial (requested n{LOSAT_THREADS})")


LOSAT_WASM_SINGLE_LABEL, LOSAT_WASM_MULTI_LABEL = wasm_mode_labels()

# === Color Settings (Seaborn Deep Palette) ===
CUSTOM_PALETTE = {
    "LOSAT": "#dd8452",
    "BLAST+": "#4c72b0",
    LOSAT_WASM_SINGLE_LABEL: "#8172b3",
    LOSAT_WASM_MULTI_LABEL: "#937860",
    LOSAT_WASM_SIMD_LABEL: "#64b5cd",
    LOSAT_WASM_SCALAR_LABEL: "#da8bc3",
}
HUE_ORDER = [
    "BLAST+",
    "LOSAT",
    LOSAT_WASM_SINGLE_LABEL,
    LOSAT_WASM_MULTI_LABEL,
    LOSAT_WASM_SIMD_LABEL,
    LOSAT_WASM_SCALAR_LABEL,
]

# === Column Definitions (outfmt 7 / 6) ===
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


def load_blast(filepath, tool_name):
    """Load a BLAST format file."""
    if not os.path.exists(filepath):
        print(f"[Warning] File not found: {filepath}")
        return None
    try:
        df = pd.read_csv(filepath, sep="\t", comment="#", names=COLUMNS, header=None)
        numeric_cols = ["pident", "length", "bitscore", "evalue"]
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors="coerce")

        df["Tool"] = tool_name
        return df
    except Exception as e:
        print(f"[Error] Failed to load {filepath}: {e}")
        return None


def strip_thread_suffix(filepath):
    return re.sub(r"\.n\d+\.out$", ".out", filepath)


def add_wasm_suffix(filepath, suffix=None):
    if not filepath.endswith(".out"):
        return f"{filepath}.wasm{'.' + suffix if suffix else ''}.out"
    stem = filepath[:-4]
    if suffix:
        return f"{stem}.wasm.{suffix}.out"
    return f"{stem}.wasm.out"


# NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:846-921
# ```c
# status =
#     s_BlastAaExtendTwoHit(query, subject, word_params, ext_params,
#                           hit_params, init_hitlist, hsp_list);
# ```
def add_wasm_build_suffix(filepath, build_mode):
    if not filepath.endswith(".out"):
        return f"{filepath}.wasm.{build_mode}.out"
    return f"{filepath[:-4]}.wasm.{build_mode}.out"


def losat_wasm_outputs(losat_path, mode_label):
    threaded_suffix = f"n{LOSAT_THREADS}"
    base_path = strip_thread_suffix(losat_path)

    outputs = [
        (LOSAT_WASM_SINGLE_LABEL, add_wasm_suffix(base_path)),
        (LOSAT_WASM_MULTI_LABEL, add_wasm_suffix(base_path, threaded_suffix)),
    ]
    if mode_label == "TBLASTX":
        outputs.extend(
            [
                (LOSAT_WASM_SIMD_LABEL, add_wasm_build_suffix(base_path, "simd")),
                (LOSAT_WASM_SCALAR_LABEL, add_wasm_build_suffix(base_path, "scalar")),
            ]
        )
    return outputs


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


def generate_comparison_plot(config):
    """Generate and save a comparison plot for a single pair."""
    name = config["name"]
    ncbi_path = config["ncbi"]
    losat_path = config["losat"]
    mode_label = config["mode"]

    print(f"Processing: {name} ({mode_label})")

    df_ncbi = load_blast(ncbi_path, "BLAST+")
    df_losat = load_blast(losat_path, "LOSAT")
    optional_frames = []
    for tool_name, wasm_path in losat_wasm_outputs(losat_path, mode_label):
        if not os.path.exists(wasm_path):
            print(f"  [Optional] Missing {tool_name}: {wasm_path}")
            continue
        warn_if_serial_wasm_thread_label(tool_name, wasm_path)
        df_wasm = load_blast(wasm_path, tool_name)
        if df_wasm is not None:
            optional_frames.append(df_wasm)

    if df_ncbi is None or df_losat is None:
        print("  -> Skipped due to missing files.")
        return

    df_merged = pd.concat([df_ncbi, df_losat, *optional_frames])
    hue_order = [tool for tool in HUE_ORDER if tool in df_merged["Tool"].unique()]

    if df_merged.empty:
        print("  -> Skipped (No data rows found).")
        return

    sns.set(style="whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f"Comparison: {name} ({mode_label})", fontsize=16)

    try:
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
    except Exception as e:
        axes[0, 0].text(0.5, 0.5, f"Error plotting hist: {e}", ha="center")

    try:
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
    except Exception as e:
        axes[0, 1].text(0.5, 0.5, f"Error plotting hist: {e}", ha="center")

    try:
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
    except Exception:
        axes[1, 0].text(0.5, 0.5, "Data error", ha="center")

    try:
        sns.scatterplot(
            data=df_merged,
            x="length",
            y="bitscore",
            hue="Tool",
            alpha=0.5,
            ax=axes[1, 1],
            palette=CUSTOM_PALETTE,
            hue_order=hue_order,
        )
        axes[1, 1].set_xscale("log")
        axes[1, 1].set_yscale("log")
        axes[1, 1].set_title("Alignment Length vs Bit Score")
        axes[1, 1].set_xlabel("Alignment Length (bp/aa)")
        axes[1, 1].set_ylabel("Bit Score")
    except Exception:
        axes[1, 1].text(0.5, 0.5, "Data error", ha="center")

    output_filename = os.path.join(PLOT_DIR, f"compare_{name}_{mode_label}.png")
    plt.tight_layout()
    plt.savefig(output_filename)
    plt.close()
    print(f"  -> Saved to {output_filename}")


TBLASTX_COMPARISONS = [
    {"name": "AP027280_Self", "mode": "TBLASTX", "ncbi": "./blast_out/AP027280.AP027280.tblastx.n8.out", "losat": "./losat_out/AP027280.AP027280.tlosatx.n8.out"},
    {"name": "NZ_CP006932_Self", "mode": "TBLASTX", "ncbi": "./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n8.out", "losat": "./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n8.out"},
    {"name": "AP027132_vs_NZ_CP006932", "mode": "TBLASTX", "ncbi": "./blast_out/AP027132.NZ_CP006932.tblastx.n8.out", "losat": "./losat_out/AP027132.NZ_CP006932.tlosatx.n8.out"},
    {"name": "AP027078_vs_AP027131", "mode": "TBLASTX", "ncbi": "./blast_out/AP027078.AP027131.tblastx.n8.out", "losat": "./losat_out/AP027078.AP027131.tlosatx.n8.out"},
    {"name": "AP027131_vs_AP027133", "mode": "TBLASTX", "ncbi": "./blast_out/AP027131.AP027133.tblastx.n8.out", "losat": "./losat_out/AP027131.AP027133.tlosatx.n8.out"},
    {"name": "AP027133_vs_AP027132", "mode": "TBLASTX", "ncbi": "./blast_out/AP027133.AP027132.tblastx.n8.out", "losat": "./losat_out/AP027133.AP027132.tlosatx.n8.out"},
    {"name": "MjeNMV_vs_MelaMJNV", "mode": "TBLASTX", "ncbi": "./blast_out/MjeNMV.MelaMJNV.tblastx.n8.out", "losat": "./losat_out/MjeNMV.MelaMJNV.tlosatx.n8.out"},
    {"name": "MelaMJNV_vs_PemoMJNVA", "mode": "TBLASTX", "ncbi": "./blast_out/MelaMJNV.PemoMJNVA.tblastx.n8.out", "losat": "./losat_out/MelaMJNV.PemoMJNVA.tlosatx.n8.out"},
    {"name": "PemoMJNVA_vs_PeseMJNV", "mode": "TBLASTX", "ncbi": "./blast_out/PemoMJNVA.PeseMJNV.tblastx.n8.out", "losat": "./losat_out/PemoMJNVA.PeseMJNV.tlosatx.n8.out"},
    {"name": "PeseMJNV_vs_PemoMJNVB", "mode": "TBLASTX", "ncbi": "./blast_out/PeseMJNV.PemoMJNVB.tblastx.n8.out", "losat": "./losat_out/PeseMJNV.PemoMJNVB.tlosatx.n8.out"},
    {"name": "PemoMJNVB_vs_LvMJNV", "mode": "TBLASTX", "ncbi": "./blast_out/PemoMJNVB.LvMJNV.tblastx.n8.out", "losat": "./losat_out/PemoMJNVB.LvMJNV.tlosatx.n8.out"},
    {"name": "LvMJNV_vs_TrcuMJNV", "mode": "TBLASTX", "ncbi": "./blast_out/LvMJNV.TrcuMJNV.tblastx.n8.out", "losat": "./losat_out/LvMJNV.TrcuMJNV.tlosatx.n8.out"},
    {"name": "TrcuMJNV_vs_MellatMJNV", "mode": "TBLASTX", "ncbi": "./blast_out/TrcuMJNV.MellatMJNV.tblastx.n8.out", "losat": "./losat_out/TrcuMJNV.MellatMJNV.tlosatx.n8.out"},
    {"name": "MellatMJNV_vs_MeenMJNV", "mode": "TBLASTX", "ncbi": "./blast_out/MellatMJNV.MeenMJNV.tblastx.n8.out", "losat": "./losat_out/MellatMJNV.MeenMJNV.tlosatx.n8.out"},
    {"name": "MeenMJNV_vs_MejoMJNV", "mode": "TBLASTX", "ncbi": "./blast_out/MeenMJNV.MejoMJNV.tblastx.n8.out", "losat": "./losat_out/MeenMJNV.MejoMJNV.tlosatx.n8.out"},
    {"name": "AvCLPV_vs_PsCLPV", "mode": "TBLASTX", "ncbi": "./blast_out/AvCLPV.PsCLPV.tblastx.n8.out", "losat": "./losat_out/AvCLPV.PsCLPV.tlosatx.n8.out"},
]

BLASTN_COMPARISONS = [
    {"name": "NZ_CP006932_Self", "mode": "BLASTN_Task", "ncbi": "./blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out", "losat": "./losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.out"},
    {"name": "PesePMNV_vs_MjPMNV", "mode": "BLASTN", "ncbi": "./blast_out/PesePMNV.MjPMNV.blastn.out", "losat": "./losat_out/PesePMNV.MjPMNV.losatn.blastn.out"},
    {"name": "MelaMJNV_vs_PemoMJNVA", "mode": "BLASTN", "ncbi": "./blast_out/MelaMJNV.PemoMJNVA.blastn.out", "losat": "./losat_out/MelaMJNV.PemoMJNVA.losatn.blastn.out"},
    {"name": "SiNMV_vs_ChdeNMV", "mode": "BLASTN", "ncbi": "./blast_out/SiNMV.ChdeNMV.blastn.out", "losat": "./losat_out/SiNMV.ChdeNMV.losatn.blastn.out"},
    {"name": "PmeNMV_vs_MjPMNV", "mode": "BLASTN", "ncbi": "./blast_out/PmeNMV.MjPMNV.blastn.out", "losat": "./losat_out/PmeNMV.MjPMNV.losatn.blastn.out"},
    {"name": "PmeNMV_vs_PesePMNV", "mode": "BLASTN", "ncbi": "./blast_out/PmeNMV.PesePMNV.blastn.out", "losat": "./losat_out/PmeNMV.PesePMNV.losatn.blastn.out"},
    {"name": "PeseMJNV_vs_PemoMJNVB", "mode": "BLASTN", "ncbi": "./blast_out/PeseMJNV.PemoMJNVB.blastn.out", "losat": "./losat_out/PeseMJNV.PemoMJNVB.losatn.blastn.out"},
    {"name": "PemoMJNVA_vs_PeseMJNV", "mode": "BLASTN", "ncbi": "./blast_out/PemoMJNVA.PeseMJNV.blastn.out", "losat": "./losat_out/PemoMJNVA.PeseMJNV.losatn.blastn.out"},
    {"name": "MjeNMV_vs_MelaMJNV", "mode": "BLASTN", "ncbi": "./blast_out/MjeNMV.MelaMJNV.blastn.out", "losat": "./losat_out/MjeNMV.MelaMJNV.losatn.blastn.out"},
    {"name": "MjPMNV_vs_MlPMNV", "mode": "BLASTN", "ncbi": "./blast_out/MjPMNV.MlPMNV.blastn.out", "losat": "./losat_out/MjPMNV.MlPMNV.losatn.blastn.out"},
    {"name": "EDL933_vs_Sakai", "mode": "Megablast", "ncbi": "./blast_out/EDL933.Sakai.blastn.megablast.out", "losat": "./losat_out/EDL933.Sakai.losatn.megablast.out"},
    {"name": "Sakai_vs_MG1655", "mode": "Megablast", "ncbi": "./blast_out/Sakai.MG1655.blastn.megablast.out", "losat": "./losat_out/Sakai.MG1655.losatn.megablast.out"},
]

BLASTP_COMPARISONS = [
    {"name": "WSSV.PajaWSV", "mode": "BLASTP", "ncbi": "./blast_out/WSSV.PajaWSV.BLASTP.out", "losat": "./losat_out/WSSV.PajaWSV.losatp.out"},
    {"name": "WSSV.SicyWSV", "mode": "BLASTP", "ncbi": "./blast_out/WSSV.SicyWSV.BLASTP.out", "losat": "./losat_out/WSSV.SicyWSV.losatp.out"},
    {"name": "WSSV.CoBV", "mode": "BLASTP", "ncbi": "./blast_out/WSSV.CoBV.BLASTP.out", "losat": "./losat_out/WSSV.CoBV.losatp.out"},
    {"name": "PajaWSV.SicyWSV", "mode": "BLASTP", "ncbi": "./blast_out/PajaWSV.SicyWSV.BLASTP.out", "losat": "./losat_out/PajaWSV.SicyWSV.losatp.out"},
    {"name": "PajaWSV.CoBV", "mode": "BLASTP", "ncbi": "./blast_out/PajaWSV.CoBV.BLASTP.out", "losat": "./losat_out/PajaWSV.CoBV.losatp.out"},
    {"name": "SicyWSV.CoBV", "mode": "BLASTP", "ncbi": "./blast_out/SicyWSV.CoBV.BLASTP.out", "losat": "./losat_out/SicyWSV.CoBV.losatp.out"},
    {"name": "AP027078.AP027131", "mode": "BLASTP", "ncbi": "./blast_out/AP027078.AP027131.BLASTP.out", "losat": "./losat_out/AP027078.AP027131.losatp.out"},
    {"name": "AP027078.AP027132", "mode": "BLASTP", "ncbi": "./blast_out/AP027078.AP027132.BLASTP.out", "losat": "./losat_out/AP027078.AP027132.losatp.out"},
    {"name": "AP027078.AP027133", "mode": "BLASTP", "ncbi": "./blast_out/AP027078.AP027133.BLASTP.out", "losat": "./losat_out/AP027078.AP027133.losatp.out"},
    {"name": "AP027078.NZ_CP006932", "mode": "BLASTP", "ncbi": "./blast_out/AP027078.NZ_CP006932.BLASTP.out", "losat": "./losat_out/AP027078.NZ_CP006932.losatp.out"},
    {"name": "AP027131.AP027132", "mode": "BLASTP", "ncbi": "./blast_out/AP027131.AP027132.BLASTP.out", "losat": "./losat_out/AP027131.AP027132.losatp.out"},
    {"name": "AP027131.AP027133", "mode": "BLASTP", "ncbi": "./blast_out/AP027131.AP027133.BLASTP.out", "losat": "./losat_out/AP027131.AP027133.losatp.out"},
    {"name": "AP027131.NZ_CP006932", "mode": "BLASTP", "ncbi": "./blast_out/AP027131.NZ_CP006932.BLASTP.out", "losat": "./losat_out/AP027131.NZ_CP006932.losatp.out"},
    {"name": "AP027132.AP027133", "mode": "BLASTP", "ncbi": "./blast_out/AP027132.AP027133.BLASTP.out", "losat": "./losat_out/AP027132.AP027133.losatp.out"},
    {"name": "AP027132.NZ_CP006932", "mode": "BLASTP", "ncbi": "./blast_out/AP027132.NZ_CP006932.BLASTP.out", "losat": "./losat_out/AP027132.NZ_CP006932.losatp.out"},
]

comparisons = TBLASTX_COMPARISONS + BLASTN_COMPARISONS + BLASTP_COMPARISONS


def main():
    print("Starting batch plot generation...")
    print(f"Output directory: {PLOT_DIR}")

    for config in comparisons:
        generate_comparison_plot(config)

    print("All plots generated.")


if __name__ == "__main__":
    main()
