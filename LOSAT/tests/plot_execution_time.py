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

# === Color Settings (Seaborn Deep Palette) ===
# Explicitly define colors to match other plots
CUSTOM_PALETTE = {"LOSAT": "#dd8452", "BLAST+": "#4c72b0"}

# === Comparison List (Full Version) ===
comparisons = [
    # --- TBLASTX ---
    {
        "name": "NZ_CP006932_Self", "mode": "TBLASTX", 
        "losat_log": "NZ_CP006932.NZ_CP006932.tlosatx.n8.log", 
        "blast_log": "NZ_CP006932.NZ_CP006932.tblastx.n1.log"
    },
    {
        "name": "AP027132_vs_NZ", "mode": "TBLASTX", 
        "losat_log": "AP027132.NZ_CP006932.tlosatx.n8.log", 
        "blast_log": "AP027132.NZ_CP006932.tblastx.n1.log"
    },
    {
        "name": "AP027078_vs_AP027131", "mode": "TBLASTX", 
        "losat_log": "AP027078.AP027131.tlosatx.n8.log", 
        "blast_log": "AP027078.AP027131.tblastx.n1.log"
    },
    {
        "name": "AP027131_vs_AP027133", "mode": "TBLASTX", 
        "losat_log": "AP027131.AP027133.tlosatx.n8.log", 
        "blast_log": "AP027131.AP027133.tblastx.n1.log"
    },
    {
        "name": "AP027133_vs_AP027132", "mode": "TBLASTX", 
        "losat_log": "AP027133.AP027132.tlosatx.n8.log", 
        "blast_log": "AP027133.AP027132.tblastx.n1.log"
    },
    {
        "name": "AP027280_Self", "mode": "TBLASTX", 
        "losat_log": "AP027280.AP027280.tlosatx.n8.log", 
        "blast_log": "AP027280.AP027280.tblastx.n1.log"
    },
         {
        "name": "MjeNMV_vs_MelaMJNV", "mode": "TBLASTX", 
        "losat_log": "MjeNMV.MelaMJNV.tlosatx.n8.log", 
        "blast_log": "MjeNMV.MelaMJNV.tblastx.n1.log"
    },
    {
        "name": "MelaMJNV_vs_PemoMJNVA", "mode": "TBLASTX", 
        "losat_log": "MelaMJNV.PemoMJNVA.tlosatx.n8.log", 
        "blast_log": "MelaMJNV.PemoMJNVA.tblastx.n1.log"
    },
    {
        "name": "PemoMJNVA_vs_PeseMJNV", "mode": "TBLASTX", 
        "losat_log": "PemoMJNVA.PeseMJNV.tlosatx.n8.log", 
        "blast_log": "PemoMJNVA.PeseMJNV.tblastx.n1.log"
    },
    {
        "name": "PeseMJNV_vs_PemoMJNVB", "mode": "TBLASTX", 
        "losat_log": "PeseMJNV.PemoMJNVB.tlosatx.n8.log", 
        "blast_log": "PeseMJNV.PemoMJNVB.tblastx.n1.log"
    },
    {
        "name": "PemoMJNVB_vs_LvMJNV", "mode": "TBLASTX", 
        "losat_log": "PemoMJNVB.LvMJNV.tlosatx.n8.log", 
        "blast_log": "PemoMJNVB.LvMJNV.tblastx.n1.log"
    },
    {
        "name": "LvMJNV_vs_TrcuMJNV", "mode": "TBLASTX", 
        "losat_log": "LvMJNV.TrcuMJNV.tlosatx.n8.log", 
        "blast_log": "LvMJNV.TrcuMJNV.tblastx.n1.log"
    },
    {   "name": "TrcuMJNV_vs_MellatMJNV", "mode": "TBLASTX",
        "losat_log": "TrcuMJNV.MellatMJNV.tlosatx.n8.log", 
        "blast_log": "TrcuMJNV.MellatMJNV.tblastx.n1.log"
    },
    {
        "name": "MellatMJNV_vs_MeenMJNV", "mode": "TBLASTX",
        "losat_log": "MellatMJNV.MeenMJNV.tlosatx.n8.log", 
        "blast_log": "MellatMJNV.MeenMJNV.tblastx.n1.log"
    },
    {
        "name": "MeenMJNV_vs_MejoMJNV", "mode": "TBLASTX",
        "losat_log": "MeenMJNV.MejoMJNV.tlosatx.n8.log", 
        "blast_log": "MeenMJNV.MejoMJNV.tblastx.n1.log"
    },
    {
        "name": "AvCLPV_vs_PsCLPV", "mode": "TBLASTX",
        "losat_log": "AvCLPV.PsCLPV.tlosatx.n8.log", 
        "blast_log": "AvCLPV.PsCLPV.tblastx.n1.log"
    },  
    # --- BLASTN (Default / Megablast) ---
    {
        "name": "NZ_CP006932_Self", "mode": "Megablast", 
        "losat_log": "NZ_CP006932.NZ_CP006932.losatn.log", 
        "blast_log": "NZ_CP006932.NZ_CP006932.blastn.log"
    },
    {
        "name": "EDL933_vs_Sakai", "mode": "Megablast", 
        "losat_log": "EDL933.Sakai.blastn.megablast.log", 
        "blast_log": "EDL933.Sakai.blastn.megablast.log"
    },
    {
        "name": "Sakai_vs_MG1655", "mode": "Megablast", 
        "losat_log": "Sakai.MG1655.blastn.megablast.log", 
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
]

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
        # LOSAT
        losat_path = os.path.join(LOG_DIR_LOSAT, item['losat_log'])
        time_losat = parse_time(losat_path)
        
        # BLAST
        blast_path = os.path.join(LOG_DIR_BLAST, item['blast_log'])
        time_blast = parse_time(blast_path)

        if time_losat is not None and time_blast is not None:
            data.append({
                'Task': item['name'],
                'Mode': item['mode'],
                'Tool': 'LOSAT',
                'Time (s)': time_losat
            })
            data.append({
                'Task': item['name'],
                'Mode': item['mode'],
                'Tool': 'BLAST+',
                'Time (s)': time_blast
            })
        else:
            # If one of the logs is missing
            missing = []
            if time_losat is None: missing.append(f"LOSAT({item['losat_log']})")
            if time_blast is None: missing.append(f"BLAST({item['blast_log']})")
            print(f"  [Skip] {item['name']}: Missing logs for {', '.join(missing)}")

    if not data:
        print("No valid time data found.")
        return

    df = pd.DataFrame(data)

    # === Create Plots ===
    sns.set(style="whitegrid")
    
    # Define custom colors: LOSAT = Orange, BLAST+ = Blue
    # Using 'tab:blue' and 'tab:orange' to match standard matplotlib/seaborn defaults
    CUSTOM_PALETTE = {"LOSAT": "#dd8452", "BLAST+": "#4c72b0"}

    # Graph settings
    g = sns.catplot(
        data=df, kind="bar",
        y="Task", x="Time (s)", hue="Tool", col="Mode",
        height=8, aspect=0.8,
        sharex=False, sharey=False,
        palette=CUSTOM_PALETTE, # Apply the custom palette
        errorbar=None,
        col_wrap=3
    )
    
    g.despine(left=True)
    g.set_axis_labels("Execution Time (seconds)", "")
    g.fig.suptitle("Execution Time Comparison: LOSAT vs BLAST+ (All Tasks)", y=1.02)
    
    # Save
    os.makedirs(os.path.dirname(OUTPUT_IMAGE), exist_ok=True)
    plt.savefig(OUTPUT_IMAGE, bbox_inches='tight')
    print(f"\nPlot saved to {OUTPUT_IMAGE}")
    
    # Display Summary Table
    print("\n--- Summary Table (Ratio < 1.0 means LOSAT is faster) ---")
    try:
        summary = df.pivot(index=['Mode', 'Task'], columns='Tool', values='Time (s)')
        summary['Ratio (LOSAT/BLAST)'] = (summary['LOSAT'] / summary['BLAST+']).round(2)
        summary['Diff (s)'] = (summary['LOSAT'] - summary['BLAST+']).round(2)
        print(summary)
    except Exception as e:
        print("Could not generate summary table:", e)

if __name__ == "__main__":
    main()