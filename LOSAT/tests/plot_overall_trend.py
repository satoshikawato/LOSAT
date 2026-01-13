#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# === Configuration ===
PLOT_DIR = "./plots"
OUTPUT_IMAGE = "./plots/overall_trend_comparison.png"
os.makedirs(PLOT_DIR, exist_ok=True)

# === Color Settings (Seaborn Deep Palette) ===
CUSTOM_PALETTE = {"LOSAT": "#dd8452", "BLAST+": "#4c72b0"}

# === Column Definitions ===
COLUMNS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
           'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

# === Comparison List ===
comparisons = [
    # --- TBLASTX ---
    {"name": "NZ_CP006932_Self", "mode": "TBLASTX", "ncbi": "./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out", "losat": "./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n8.out"},
    {"name": "AP027132_vs_NZ", "mode": "TBLASTX", "ncbi": "./blast_out/AP027132.NZ_CP006932.tblastx.n1.out", "losat": "./losat_out/AP027132.NZ_CP006932.tlosatx.n8.out"},
    {"name": "AP027078_vs_AP027131", "mode": "TBLASTX", "ncbi": "./blast_out/AP027078.AP027131.tblastx.n1.out", "losat": "./losat_out/AP027078.AP027131.tlosatx.n8.out"},
    {"name": "AP027131_vs_AP027133", "mode": "TBLASTX", "ncbi": "./blast_out/AP027131.AP027133.tblastx.n1.out", "losat": "./losat_out/AP027131.AP027133.tlosatx.n8.out"},
    {"name": "AP027133_vs_AP027132", "mode": "TBLASTX", "ncbi": "./blast_out/AP027133.AP027132.tblastx.n1.out", "losat": "./losat_out/AP027133.AP027132.tlosatx.n8.out"},
    {"name": "AP027280_Self", "mode": "TBLASTX", "ncbi": "./blast_out/AP027280.AP027280.tblastx.n1.out", "losat": "./losat_out/AP027280.AP027280.tlosatx.n8.out"},
    {
        "name": "MjeNMV_vs_MelaMJNV", "mode": "TBLASTX", 
        "ncbi": "./blast_out/MjeNMV.MelaMJNV.tblastx.n1.out", 
        "losat": "./losat_out/MjeNMV.MelaMJNV.tlosatx.n8.out"
    },
    {
        "name": "MjeNMV_vs_MelaMJNV", "mode": "TBLASTX", 
        "ncbi": "./blast_out/AP027280.AP027280.tblastx.n1.out", 
        "losat": "./losat_out/AP027280.AP027280.tlosatx.n8.out"
    },
    {
        "name": "MelaMJNV_vs_PemoMJNVA", "mode": "TBLASTX", 
        "ncbi": "./blast_out/MelaMJNV.PemoMJNVA.tblastx.n1.out", 
        "losat": "./losat_out/MelaMJNV.PemoMJNVA.tlosatx.n8.out"
    },
    {
        "name": "PemoMJNVA_vs_PeseMJNV", "mode": "TBLASTX", 
        "ncbi": "./blast_out/PemoMJNVA.PeseMJNV.tblastx.n1.out", 
        "losat": "./losat_out/PemoMJNVA.PeseMJNV.tlosatx.n8.out"
    },
    {
        "name": "PeseMJNV_vs_PemoMJNVB", "mode": "TBLASTX", 
        "ncbi": "./blast_out/PeseMJNV.PemoMJNVB.tblastx.n1.out", 
        "losat": "./losat_out/PeseMJNV.PemoMJNVB.tlosatx.n8.out"
    },
    {
        "name": "PemoMJNVB_vs_LvMJNV", "mode": "TBLASTX", 
        "ncbi": "./blast_out/PemoMJNVB.LvMJNV.tblastx.n1.out", 
        "losat": "./losat_out/PemoMJNVB.LvMJNV.tlosatx.n8.out"
    },
    {
        "name": "LvMJNV_vs_TrcuMJNV", "mode": "TBLASTX", 
        "ncbi": "./blast_out/LvMJNV.TrcuMJNV.tblastx.n1.out", 
        "losat": "./losat_out/LvMJNV.TrcuMJNV.tlosatx.n8.out"
    },
    {   "name": "TrcuMJNV_vs_MellatMJNV", "mode": "TBLASTX",
        "ncbi": "./blast_out/TrcuMJNV.MellatMJNV.tblastx.n1.out", 
        "losat": "./losat_out/TrcuMJNV.MellatMJNV.tlosatx.n8.out"
    },
    {
        "name": "MellatMJNV_vs_MeenMJNV", "mode": "TBLASTX",
        "ncbi": "./blast_out/MellatMJNV.MeenMJNV.tblastx.n1.out", 
        "losat": "./losat_out/MellatMJNV.MeenMJNV.tlosatx.n8.out"
    },
    {
        "name": "MeenMJNV_vs_MejoMJNV", "mode": "TBLASTX",
        "ncbi": "./blast_out/MeenMJNV.MejoMJNV.tblastx.n1.out", 
        "losat": "./losat_out/MeenMJNV.MejoMJNV.tlosatx.n8.out"
    },
    {
        "name": "AvCLPV_vs_PsCLPV", "mode": "TBLASTX",
        "ncbi": "./blast_out/AvCLPV.PsCLPV.tblastx.n1.out", 
        "losat": "./losat_out/AvCLPV.PsCLPV.tlosatx.n8.out"
    },
    # --- BLASTN / Megablast ---
    {"name": "NZ_CP006932_Self(Mega)", "mode": "BLASTN (Megablast)", "ncbi": "./blast_out/NZ_CP006932.NZ_CP006932.blastn.out", "losat": "./losat_out/NZ_CP006932.NZ_CP006932.losatn.out"},
    {"name": "EDL933_vs_Sakai", "mode": "BLASTN (Megablast)", "ncbi": "./blast_out/EDL933.Sakai.blastn.megablast.out", "losat": "./losat_out/EDL933.Sakai.blastn.megablast.out"},
    {"name": "Sakai_vs_MG1655", "mode": "BLASTN (Megablast)", "ncbi": "./blast_out/Sakai.MG1655.blastn.megablast.out", "losat": "./losat_out/Sakai.MG1655.blastn.megablast.out"},

    # --- BLASTN (Task: blastn) ---
    {"name": "NZ_CP006932_Self(Task)", "mode": "BLASTN (Task:blastn)", "ncbi": "./blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out", "losat": "./losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.out"},
    {"name": "PesePMNV_vs_MjPMNV", "mode": "BLASTN (Task:blastn)", "ncbi": "./blast_out/PesePMNV.MjPMNV.blastn.out", "losat": "./losat_out/PesePMNV.MjPMNV.losatn.blastn.out"},
    {"name": "MelaMJNV_vs_PemoMJNVA", "mode": "BLASTN (Task:blastn)", "ncbi": "./blast_out/MelaMJNV.PemoMJNVA.blastn.out", "losat": "./losat_out/MelaMJNV.PemoMJNVA.losatn.blastn.out"},
    {"name": "SiNMV_vs_ChdeNMV", "mode": "BLASTN (Task:blastn)", "ncbi": "./blast_out/SiNMV.ChdeNMV.blastn.out", "losat": "./losat_out/SiNMV.ChdeNMV.losatn.blastn.out"},
    {"name": "PmeNMV_vs_MjPMNV", "mode": "BLASTN (Task:blastn)", "ncbi": "./blast_out/PmeNMV.MjPMNV.blastn.out", "losat": "./losat_out/PmeNMV.MjPMNV.losatn.blastn.out"},
    {"name": "PmeNMV_vs_PesePMNV", "mode": "BLASTN (Task:blastn)", "ncbi": "./blast_out/PmeNMV.PesePMNV.blastn.out", "losat": "./losat_out/PmeNMV.PesePMNV.losatn.blastn.out"},
    {"name": "PeseMJNV_vs_PemoMJNVB", "mode": "BLASTN (Task:blastn)", "ncbi": "./blast_out/PeseMJNV.PemoMJNVB.blastn.out", "losat": "./losat_out/PeseMJNV.PemoMJNVB.losatn.blastn.out"},
    {"name": "PemoMJNVA_vs_PeseMJNV", "mode": "BLASTN (Task:blastn)", "ncbi": "./blast_out/PemoMJNVA.PeseMJNV.blastn.out", "losat": "./losat_out/PemoMJNVA.PeseMJNV.losatn.blastn.out"},
    {"name": "MjeNMV_vs_MelaMJNV", "mode": "BLASTN (Task:blastn)", "ncbi": "./blast_out/MjeNMV.MelaMJNV.blastn.out", "losat": "./losat_out/MjeNMV.MelaMJNV.losatn.blastn.out"},
    {"name": "MjPMNV_vs_MlPMNV", "mode": "BLASTN (Task:blastn)", "ncbi": "./blast_out/MjPMNV.MlPMNV.blastn.out", "losat": "./losat_out/MjPMNV.MlPMNV.losatn.blastn.out"},
]

def load_blast(filepath, tool_name, mode_name, task_name):
    """Load a BLAST file and add metadata columns."""
    if not os.path.exists(filepath):
        return None
    try:
        df = pd.read_csv(filepath, sep='\t', comment='#', names=COLUMNS, header=None)
        numeric_cols = ['pident', 'length', 'bitscore', 'evalue']
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        df['Tool'] = tool_name
        df['Mode'] = mode_name 
        df['Task'] = task_name 
        return df
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None

def main():
    print("Loading all datasets...")
    all_data = []

    for item in comparisons:
        # Load NCBI BLAST+
        df_ncbi = load_blast(item['ncbi'], 'BLAST+', item['mode'], item['name'])
        if df_ncbi is not None: all_data.append(df_ncbi)
        
        # Load LOSAT
        df_losat = load_blast(item['losat'], 'LOSAT', item['mode'], item['name'])
        if df_losat is not None: all_data.append(df_losat)

    if not all_data:
        print("No data found.")
        return

    df_all = pd.concat(all_data, ignore_index=True)
    print(f"Total alignment records loaded: {len(df_all)}")

# === Plotting ===
    df_all['Broad_Mode'] = df_all['Mode'].apply(lambda x: 'TBLASTX' if 'TBLASTX' in x else 'BLASTN (All Types)')

    sns.set(style="whitegrid")
    
    modes = df_all['Broad_Mode'].unique()
    fig, axes = plt.subplots(len(modes), 3, figsize=(18, 6 * len(modes)))
    
    if len(modes) == 1:
        axes = [axes]

    for i, mode in enumerate(modes):
        data_subset = df_all[df_all['Broad_Mode'] == mode]
        
        # Row Title
        axes[i][0].text(-0.2, 0.5, mode, transform=axes[i][0].transAxes, 
                        fontsize=16, rotation=90, va='center', fontweight='bold')

        # 1. Accumulated Length vs Alignment Length
        # 修正: bins=50 を追加
        sns.histplot(data=data_subset, x='length', hue='Tool', 
                     weights='length', bins=100,
                     element="step", stat="count", common_norm=False, 
                     log_scale=True, ax=axes[i][0],
                     palette=CUSTOM_PALETTE)
        axes[i][0].set_title(f'Accumulated Length vs Alignment Length')
        axes[i][0].set_xlabel('Length (bp or aa)')
        axes[i][0].set_ylabel('Accumulated Length (bp/aa)')

        # 2. Accumulated Length vs Identity
        # 修正: bins=50 を追加
        sns.histplot(data=data_subset, x='pident', hue='Tool', 
                     weights='length', bins=100,
                     element="step", stat="count", common_norm=False, 
                     ax=axes[i][1],
                     palette=CUSTOM_PALETTE)
        axes[i][1].set_title(f'Accumulated Length vs Identity')
        axes[i][1].set_xlabel('Identity (%)')
        axes[i][1].set_ylabel('Accumulated Length (bp/aa)')

        # 3. Scatter: Length vs Identity
        if len(data_subset) > 100000:
            plot_data = data_subset.sample(100000, random_state=42)
            title_suffix = "(Subsampled 100k)"
        else:
            plot_data = data_subset
            title_suffix = ""
            
        sns.scatterplot(data=plot_data, x='length', y='pident', hue='Tool', 
                        style='Tool', alpha=0.5, ax=axes[i][2],
                        palette=CUSTOM_PALETTE)
        axes[i][2].set_xscale('log')
        axes[i][2].set_title(f'Length vs Identity {title_suffix}')
        axes[i][2].set_xlabel('Length (bp or aa)')
        axes[i][2].set_ylabel('Identity (%)')

    plt.suptitle("Overall Trends: LOSAT vs BLAST+ (All Tasks Combined)", fontsize=20, y=1.02)
    plt.tight_layout()
    plt.savefig(OUTPUT_IMAGE, bbox_inches='tight')
    print(f"Overall trend plot saved to {OUTPUT_IMAGE}")

if __name__ == "__main__":
    main()