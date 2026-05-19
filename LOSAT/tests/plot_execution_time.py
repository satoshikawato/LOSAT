#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re

# === Configuration: Log File Directories ===
LOG_DIR_LOSAT = "./losat_out"
LOG_DIR_BLAST = "./blast_out"
OUTPUT_IMAGE = "./plots/execution_time_comparison_all.png"
LOSAT_THREADS = (
    os.environ.get("LOSAT_THREADS")
    or os.environ.get("LOSATP_THREADS")
    or os.environ.get("LOSAT_BLASTP_THREADS")
    or "8"
)
LOSAT_NATIVE_SINGLE_LABEL = "LOSAT native n1"
LOSAT_NATIVE_MULTI_LABEL = f"LOSAT native n{LOSAT_THREADS}"
LOSAT_WASM_SINGLE_LABEL = "LOSAT wasm n1"
LOSAT_WASM_MULTI_LABEL = f"LOSAT wasm n{LOSAT_THREADS}"

# === Color Settings (Seaborn Deep Palette) ===
# Explicitly define colors to match other plots
CUSTOM_PALETTE = {
    "BLAST+": "#4c72b0",
    "LOSAT": "#dd8452",
    LOSAT_NATIVE_SINGLE_LABEL: "#55a868",
    LOSAT_NATIVE_MULTI_LABEL: "#c44e52",
    LOSAT_WASM_SINGLE_LABEL: "#8172b3",
    LOSAT_WASM_MULTI_LABEL: "#937860",
}

# === Comparison List (Full Version) ===
comparisons = [
    # --- TBLASTX ---
    {
        "name": "NZ_CP006932_Self", "mode": "TBLASTX", 
        "losat_log": "NZ_CP006932.NZ_CP006932.tlosatx.n8.log", 
        "blast_log": "NZ_CP006932.NZ_CP006932.tblastx.n8.log"
    },
    {
        "name": "AP027132_vs_NZ", "mode": "TBLASTX", 
        "losat_log": "AP027132.NZ_CP006932.tlosatx.n8.log", 
        "blast_log": "AP027132.NZ_CP006932.tblastx.n8.log"
    },
    {
        "name": "AP027078_vs_AP027131", "mode": "TBLASTX", 
        "losat_log": "AP027078.AP027131.tlosatx.n8.log", 
        "blast_log": "AP027078.AP027131.tblastx.n8.log"
    },
    {
        "name": "AP027131_vs_AP027133", "mode": "TBLASTX", 
        "losat_log": "AP027131.AP027133.tlosatx.n8.log", 
        "blast_log": "AP027131.AP027133.tblastx.n8.log"
    },
    {
        "name": "AP027133_vs_AP027132", "mode": "TBLASTX", 
        "losat_log": "AP027133.AP027132.tlosatx.n8.log", 
        "blast_log": "AP027133.AP027132.tblastx.n8.log"
    },
    {
        "name": "AP027280_Self", "mode": "TBLASTX", 
        "losat_log": "AP027280.AP027280.tlosatx.n8.log", 
        "blast_log": "AP027280.AP027280.tblastx.n8.log"
    },
         {
        "name": "MjeNMV_vs_MelaMJNV", "mode": "TBLASTX", 
        "losat_log": "MjeNMV.MelaMJNV.tlosatx.n8.log", 
        "blast_log": "MjeNMV.MelaMJNV.tblastx.n8.log"
    },
    {
        "name": "MelaMJNV_vs_PemoMJNVA", "mode": "TBLASTX", 
        "losat_log": "MelaMJNV.PemoMJNVA.tlosatx.n8.log", 
        "blast_log": "MelaMJNV.PemoMJNVA.tblastx.n8.log"
    },
    {
        "name": "PemoMJNVA_vs_PeseMJNV", "mode": "TBLASTX", 
        "losat_log": "PemoMJNVA.PeseMJNV.tlosatx.n8.log", 
        "blast_log": "PemoMJNVA.PeseMJNV.tblastx.n8.log"
    },
    {
        "name": "PeseMJNV_vs_PemoMJNVB", "mode": "TBLASTX", 
        "losat_log": "PeseMJNV.PemoMJNVB.tlosatx.n8.log", 
        "blast_log": "PeseMJNV.PemoMJNVB.tblastx.n8.log"
    },
    {
        "name": "PemoMJNVB_vs_LvMJNV", "mode": "TBLASTX", 
        "losat_log": "PemoMJNVB.LvMJNV.tlosatx.n8.log", 
        "blast_log": "PemoMJNVB.LvMJNV.tblastx.n8.log"
    },
    {
        "name": "LvMJNV_vs_TrcuMJNV", "mode": "TBLASTX", 
        "losat_log": "LvMJNV.TrcuMJNV.tlosatx.n8.log", 
        "blast_log": "LvMJNV.TrcuMJNV.tblastx.n8.log"
    },
    {   "name": "TrcuMJNV_vs_MellatMJNV", "mode": "TBLASTX",
        "losat_log": "TrcuMJNV.MellatMJNV.tlosatx.n8.log", 
        "blast_log": "TrcuMJNV.MellatMJNV.tblastx.n8.log"
    },
    {
        "name": "MellatMJNV_vs_MeenMJNV", "mode": "TBLASTX",
        "losat_log": "MellatMJNV.MeenMJNV.tlosatx.n8.log", 
        "blast_log": "MellatMJNV.MeenMJNV.tblastx.n8.log"
    },
    {
        "name": "MeenMJNV_vs_MejoMJNV", "mode": "TBLASTX",
        "losat_log": "MeenMJNV.MejoMJNV.tlosatx.n8.log", 
        "blast_log": "MeenMJNV.MejoMJNV.tblastx.n8.log"
    },
    {
        "name": "AvCLPV_vs_PsCLPV", "mode": "TBLASTX",
        "losat_log": "AvCLPV.PsCLPV.tlosatx.n8.log", 
        "blast_log": "AvCLPV.PsCLPV.tblastx.n8.log"
    },  
    # --- BLASTN (Default / Megablast) ---
    {
        "name": "NZ_CP006932_Self", "mode": "Megablast", 
        "losat_log": "NZ_CP006932.NZ_CP006932.losatn.megablast.log", 
        "blast_log": "NZ_CP006932.NZ_CP006932.blastn.log"
    },
    {
        "name": "EDL933_vs_Sakai", "mode": "Megablast", 
        "losat_log": "EDL933.Sakai.losatn.megablast.log", 
        "blast_log": "EDL933.Sakai.blastn.megablast.log"
    },
    {
        "name": "Sakai_vs_MG1655", "mode": "Megablast", 
        "losat_log": "Sakai.MG1655.losatn.megablast.log", 
        "blast_log": "Sakai.MG1655.blastn.megablast.log"
    },

    # --- BLASTN (Task: blastn) ---
    {
        "name": "NZ_CP006932_Self", "mode": "BLASTN", 
        "losat_log": "NZ_CP006932.NZ_CP006932.losatn.blastn.log", 
        "blast_log": "NZ_CP006932.NZ_CP006932.task_blastn.log"
    },
    {
        "name": "PesePMNV_vs_MjPMNV", "mode": "BLASTN", 
        "losat_log": "PesePMNV.MjPMNV.losatn.blastn.log", 
        "blast_log": "PesePMNV.MjPMNV.blastn.log"
    },
    {
        "name": "MelaMJNV_vs_PemoMJNVA", "mode": "BLASTN", 
        "losat_log": "MelaMJNV.PemoMJNVA.losatn.blastn.log", 
        "blast_log": "MelaMJNV.PemoMJNVA.blastn.log"
    },
    {
        "name": "SiNMV_vs_ChdeNMV", "mode": "BLASTN", 
        "losat_log": "SiNMV.ChdeNMV.losatn.blastn.log", 
        "blast_log": "SiNMV.ChdeNMV.blastn.log"
    },
    {
        "name": "PmeNMV_vs_MjPMNV", "mode": "BLASTN", 
        "losat_log": "PmeNMV.MjPMNV.losatn.blastn.log", 
        "blast_log": "PmeNMV.MjPMNV.blastn.log"
    },
    {
        "name": "PmeNMV_vs_PesePMNV", "mode": "BLASTN", 
        "losat_log": "PmeNMV.PesePMNV.losatn.blastn.log", 
        "blast_log": "PmeNMV.PesePMNV.blastn.log"
    },
    {
        "name": "PeseMJNV_vs_PemoMJNVB", "mode": "BLASTN", 
        "losat_log": "PeseMJNV.PemoMJNVB.losatn.blastn.log", 
        "blast_log": "PeseMJNV.PemoMJNVB.blastn.log"
    },
    {
        "name": "PemoMJNVA_vs_PeseMJNV", "mode": "BLASTN", 
        "losat_log": "PemoMJNVA.PeseMJNV.losatn.blastn.log", 
        "blast_log": "PemoMJNVA.PeseMJNV.blastn.log"
    },
    {
        "name": "MjeNMV_vs_MelaMJNV", "mode": "BLASTN", 
        "losat_log": "MjeNMV.MelaMJNV.losatn.blastn.log", 
        "blast_log": "MjeNMV.MelaMJNV.blastn.log"
    },
    {
        "name": "MjPMNV_vs_MlPMNV", "mode": "BLASTN", 
        "losat_log": "MjPMNV.MlPMNV.losatn.blastn.log", 
        "blast_log": "MjPMNV.MlPMNV.blastn.log"
    },
    # --- BLASTP ---
    {
        "name": "WSSV.PajaWSV", "mode": "BLASTP",
        "losat_log": "WSSV.PajaWSV.losatp.log",
        "blast_log": "WSSV.PajaWSV.BLASTP.log"
    },
    {
        "name": "WSSV.SicyWSV", "mode": "BLASTP",
        "losat_log": "WSSV.SicyWSV.losatp.log",
        "blast_log": "WSSV.SicyWSV.BLASTP.log"
    },
    {
        "name": "WSSV.CoBV", "mode": "BLASTP",
        "losat_log": "WSSV.CoBV.losatp.log",
        "blast_log": "WSSV.CoBV.BLASTP.log"
    },
    {
        "name": "PajaWSV.SicyWSV", "mode": "BLASTP",
        "losat_log": "PajaWSV.SicyWSV.losatp.log",
        "blast_log": "PajaWSV.SicyWSV.BLASTP.log"
    },
    {
        "name": "PajaWSV.CoBV", "mode": "BLASTP",
        "losat_log": "PajaWSV.CoBV.losatp.log",
        "blast_log": "PajaWSV.CoBV.BLASTP.log"
    },
    {
        "name": "SicyWSV.CoBV", "mode": "BLASTP",
        "losat_log": "SicyWSV.CoBV.losatp.log",
        "blast_log": "SicyWSV.CoBV.BLASTP.log"
    },
    {
        "name": "AP027078.AP027131", "mode": "BLASTP",
        "losat_log": "AP027078.AP027131.losatp.log",
        "blast_log": "AP027078.AP027131.BLASTP.log"
    },
    {
        "name": "AP027078.AP027132", "mode": "BLASTP",
        "losat_log": "AP027078.AP027132.losatp.log",
        "blast_log": "AP027078.AP027132.BLASTP.log"
    },
    {
        "name": "AP027078.AP027133", "mode": "BLASTP",
        "losat_log": "AP027078.AP027133.losatp.log",
        "blast_log": "AP027078.AP027133.BLASTP.log"
    },
    {
        "name": "AP027078.NZ_CP006932", "mode": "BLASTP",
        "losat_log": "AP027078.NZ_CP006932.losatp.log",
        "blast_log": "AP027078.NZ_CP006932.BLASTP.log"
    },
    {
        "name": "AP027131.AP027132", "mode": "BLASTP",
        "losat_log": "AP027131.AP027132.losatp.log",
        "blast_log": "AP027131.AP027132.BLASTP.log"
    },
    {
        "name": "AP027131.AP027133", "mode": "BLASTP",
        "losat_log": "AP027131.AP027133.losatp.log",
        "blast_log": "AP027131.AP027133.BLASTP.log"
    },
    {
        "name": "AP027131.NZ_CP006932", "mode": "BLASTP",
        "losat_log": "AP027131.NZ_CP006932.losatp.log",
        "blast_log": "AP027131.NZ_CP006932.BLASTP.log"
    },
    {
        "name": "AP027132.AP027133", "mode": "BLASTP",
        "losat_log": "AP027132.AP027133.losatp.log",
        "blast_log": "AP027132.AP027133.BLASTP.log"
    },
    {
        "name": "AP027132.NZ_CP006932", "mode": "BLASTP",
        "losat_log": "AP027132.NZ_CP006932.losatp.log",
        "blast_log": "AP027132.NZ_CP006932.BLASTP.log"
    },
]

# NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94
# ```c
# const string kArgNumThreads("num_threads");
# ```
def add_thread_suffix(log_name, suffix):
    if not log_name.endswith(".log"):
        return f"{log_name}.{suffix}.log"
    return f"{log_name[:-4]}.{suffix}.log"


def strip_thread_suffix(log_name):
    return re.sub(r"\.n\d+\.log$", ".log", log_name)


def normalize_thread_suffix(log_name, suffix):
    if is_explicit_thread_log(log_name):
        return re.sub(r"\.n\d+\.log$", f".{suffix}.log", log_name)
    return add_thread_suffix(log_name, suffix)


def add_wasm_suffix(log_name, suffix=None):
    if not log_name.endswith(".log"):
        return f"{log_name}.wasm{'.' + suffix if suffix else ''}.log"
    stem = log_name[:-4]
    if suffix:
        return f"{stem}.wasm.{suffix}.log"
    return f"{stem}.wasm.log"

def is_explicit_thread_log(log_name):
    return re.search(r"\.n\d+\.log$", log_name) is not None

def parse_time(filepath):
    """Extract execution time (in seconds) from a log file."""
    if not os.path.exists(filepath):
        # Return None to handle missing files later
        return None
    
    try:
        with open(filepath, 'r') as f:
            content = f.read()
            
        # Match format: "real XmY.Ys" (bash time)
        match = re.search(r'real\s+(\d+)m([\d\.]+)s', content)
        if match:
            return float(match.group(1)) * 60 + float(match.group(2))
        
        # Match format: "X.XXuser ... X:XX.XXelapsed" (GNU time)
        match_elapsed = re.search(r'(\d+):([\d\.]+)elapsed', content)
        if match_elapsed:
             return float(match_elapsed.group(1)) * 60 + float(match_elapsed.group(2))

        return None
        
    except Exception as e:
        print(f"[Error] Reading {filepath}: {e}")
        return None

def main():
    data = []
    print(f"Checking {len(comparisons)} pairs...")

    for item in comparisons:
        missing = []
        threaded_suffix = f"n{LOSAT_THREADS}"
        candidate_logs = [
            ("BLAST+", LOG_DIR_BLAST, normalize_thread_suffix(item["blast_log"], threaded_suffix)
             if is_explicit_thread_log(item["blast_log"]) else item["blast_log"]),
        ]

        losat_log = item["losat_log"]
        if is_explicit_thread_log(losat_log):
            base_losat_log = strip_thread_suffix(losat_log)
            candidate_logs.append(
                (LOSAT_NATIVE_SINGLE_LABEL, LOG_DIR_LOSAT, normalize_thread_suffix(losat_log, "n1"))
            )
            candidate_logs.append(
                (LOSAT_NATIVE_MULTI_LABEL, LOG_DIR_LOSAT, normalize_thread_suffix(losat_log, threaded_suffix))
            )
            candidate_logs.append((LOSAT_WASM_SINGLE_LABEL, LOG_DIR_LOSAT, add_wasm_suffix(base_losat_log)))
            candidate_logs.append(
                (LOSAT_WASM_MULTI_LABEL, LOG_DIR_LOSAT, add_wasm_suffix(base_losat_log, threaded_suffix))
            )
        else:
            candidate_logs.append((LOSAT_NATIVE_SINGLE_LABEL, LOG_DIR_LOSAT, losat_log))
            candidate_logs.append(
                (LOSAT_NATIVE_MULTI_LABEL, LOG_DIR_LOSAT, add_thread_suffix(losat_log, threaded_suffix))
            )
            candidate_logs.append((LOSAT_WASM_SINGLE_LABEL, LOG_DIR_LOSAT, add_wasm_suffix(losat_log)))
            candidate_logs.append(
                (LOSAT_WASM_MULTI_LABEL, LOG_DIR_LOSAT, add_wasm_suffix(losat_log, threaded_suffix))
            )

        for tool, log_dir, log_name in candidate_logs:
            time_value = parse_time(os.path.join(log_dir, log_name))
            if time_value is None:
                missing.append(f"{tool}({log_name})")
                continue
            data.append({
                "Task": item["name"],
                "Mode": item["mode"],
                "Tool": tool,
                "Time (s)": time_value,
            })

        if missing:
            print(f"  [Partial] {item['name']}: Missing logs for {', '.join(missing)}")

    if not data:
        print("No valid time data found.")
        return

    df = pd.DataFrame(data)

    # === Create Plots ===
    sns.set(style="whitegrid")
    
    # Graph settings
    mode_order = ["TBLASTX", "Megablast", "BLASTN", "BLASTP"]
    tool_order = [
        "BLAST+",
        LOSAT_NATIVE_SINGLE_LABEL,
        LOSAT_NATIVE_MULTI_LABEL,
        LOSAT_WASM_SINGLE_LABEL,
        LOSAT_WASM_MULTI_LABEL,
        "LOSAT",
    ]
    g = sns.catplot(
        data=df, kind="bar",
        y="Task", x="Time (s)", hue="Tool", col="Mode",
        height=8, aspect=0.8,
        sharex=False, sharey=False,
        palette=CUSTOM_PALETTE, # Apply the custom palette
        hue_order=[tool for tool in tool_order if tool in df["Tool"].unique()],
        errorbar=None,
        col_wrap=4,
        col_order=[mode for mode in mode_order if mode in df["Mode"].unique()]
    )
    
    g.despine(left=True)
    g.set_axis_labels("Execution Time (seconds)", "")
    g.fig.suptitle("Execution Time Comparison: BLAST+ vs LOSAT Native/Wasm Threads", y=1.02)
    
    # Save
    os.makedirs(os.path.dirname(OUTPUT_IMAGE), exist_ok=True)
    plt.savefig(OUTPUT_IMAGE, bbox_inches='tight')
    print(f"\nPlot saved to {OUTPUT_IMAGE}")
    
    # Display Summary Table
    print("\n--- Summary Table ---")
    try:
        summary = df.pivot(index=['Mode', 'Task'], columns='Tool', values='Time (s)')
        if {LOSAT_NATIVE_SINGLE_LABEL, "BLAST+"}.issubset(summary.columns):
            summary['Ratio (native n1/BLAST)'] = (
                summary[LOSAT_NATIVE_SINGLE_LABEL] / summary['BLAST+']
            ).round(2)
        if {LOSAT_NATIVE_SINGLE_LABEL, LOSAT_NATIVE_MULTI_LABEL}.issubset(summary.columns):
            summary[f"Speedup (native n1/native n{LOSAT_THREADS})"] = (
                summary[LOSAT_NATIVE_SINGLE_LABEL] / summary[LOSAT_NATIVE_MULTI_LABEL]
            ).round(2)
        if {LOSAT_WASM_SINGLE_LABEL, LOSAT_WASM_MULTI_LABEL}.issubset(summary.columns):
            summary[f"Speedup (wasm n1/wasm n{LOSAT_THREADS})"] = (
                summary[LOSAT_WASM_SINGLE_LABEL] / summary[LOSAT_WASM_MULTI_LABEL]
            ).round(2)
        if {LOSAT_NATIVE_SINGLE_LABEL, LOSAT_WASM_SINGLE_LABEL}.issubset(summary.columns):
            summary["Ratio (wasm n1/native n1)"] = (
                summary[LOSAT_WASM_SINGLE_LABEL] / summary[LOSAT_NATIVE_SINGLE_LABEL]
            ).round(2)
        print(summary)
    except Exception as e:
        print("Could not generate summary table:", e)

if __name__ == "__main__":
    main()
