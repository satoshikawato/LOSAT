#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os


# === Configuration: Output Directory ===
PLOT_DIR = "./plots"
os.makedirs(PLOT_DIR, exist_ok=True)

# === Color Settings (Seaborn Deep Palette) ===
CUSTOM_PALETTE = {"LOSAT": "#dd8452", "BLAST+": "#4c72b0"}

# === Column Definitions (outfmt 7 / 6) ===
COLUMNS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
           'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

def load_blast(filepath, tool_name):
    """Load a BLAST format file."""
    if not os.path.exists(filepath):
        print(f"[Warning] File not found: {filepath}")
        return None
    try:
        df = pd.read_csv(filepath, sep='\t', comment='#', names=COLUMNS, header=None)
        numeric_cols = ['pident', 'length', 'bitscore', 'evalue']
        for col in numeric_cols:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        
        df['Tool'] = tool_name
        return df
    except Exception as e:
        print(f"[Error] Failed to load {filepath}: {e}")
        return None

def generate_comparison_plot(config):
    """Generate and save a comparison plot for a single pair."""
    name = config['name']
    ncbi_path = config['ncbi']
    losat_path = config['losat']
    mode_label = config['mode']

    print(f"Processing: {name} ({mode_label})")

    df_ncbi = load_blast(ncbi_path, 'BLAST+')
    df_losat = load_blast(losat_path, 'LOSAT')

    if df_ncbi is None or df_losat is None:
        print(f"  -> Skipped due to missing files.")
        return

    df_merged = pd.concat([df_ncbi, df_losat])

    if df_merged.empty:
        print("  -> Skipped (No data rows found).")
        return

    # === Create Plots ===
    sns.set(style="whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f"Comparison: {name} ({mode_label})", fontsize=16)

    # 1. Accumulated Alignment Length by Length (Histogram)
    try:
        # 修正: bins=50 を追加
        sns.histplot(data=df_merged, x='length', hue='Tool', 
                     weights='length', bins=100,
                     element="step", stat="count", common_norm=False, 
                     log_scale=True, ax=axes[0,0],
                     palette=CUSTOM_PALETTE)
        axes[0,0].set_title('Accumulated Length vs Alignment Length')
        axes[0,0].set_xlabel('Alignment Length (bp/aa)')
        axes[0,0].set_ylabel('Accumulated Length (bp/aa)')
    except Exception as e:
        axes[0,0].text(0.5, 0.5, f"Error plotting hist: {e}", ha='center')

    # 2. Accumulated Alignment Length by Identity (Histogram)
    try:
        # 修正: bins=50 を追加
        sns.histplot(data=df_merged, x='pident', hue='Tool', 
                     weights='length', bins=100,
                     element="step", stat="count", common_norm=False, 
                     ax=axes[0,1],
                     palette=CUSTOM_PALETTE)
        axes[0,1].set_title('Accumulated Length vs Identity')
        axes[0,1].set_xlabel('Identity (%)')
        axes[0,1].set_ylabel('Accumulated Length (bp/aa)')
    except Exception as e:
        axes[0,1].text(0.5, 0.5, f"Error plotting hist: {e}", ha='center')

    # 3. Alignment Length vs Identity (Scatter)
    try:
        sns.scatterplot(data=df_merged, x='length', y='pident', hue='Tool', 
                        alpha=0.5, style='Tool', ax=axes[1,0],
                        palette=CUSTOM_PALETTE)
        axes[1,0].set_xscale('log')
        axes[1,0].set_title('Alignment Length vs Identity')
        axes[1,0].set_xlabel('Alignment Length (bp/aa)')
        axes[1,0].set_ylabel('Identity (%)')
    except Exception as e:
         axes[1,0].text(0.5, 0.5, "Data error", ha='center')

    # 4. Alignment Length vs Bit Score (Scatter)
    try:
        sns.scatterplot(data=df_merged, x='length', y='bitscore', hue='Tool', 
                        alpha=0.5, ax=axes[1,1],
                        palette=CUSTOM_PALETTE)
        axes[1,1].set_xscale('log')
        axes[1,1].set_yscale('log')
        axes[1,1].set_title('Alignment Length vs Bit Score')
        axes[1,1].set_xlabel('Alignment Length (bp/aa)')
        axes[1,1].set_ylabel('Bit Score')
    except Exception as e:
         axes[1,1].text(0.5, 0.5, "Data error", ha='center')

    # Save
    output_filename = os.path.join(PLOT_DIR, f"compare_{name}_{mode_label}.png")
    plt.tight_layout()
    plt.savefig(output_filename)
    plt.close()
    print(f"  -> Saved to {output_filename}")


# === Comparison List ===
comparisons = [
    # --- TBLASTX ---
    {
        "name": "AP027280_Self",
        "mode": "TBLASTX",
        "ncbi":  "./blast_out/AP027280.AP027280.tblastx.n8.out",
        "losat": "./losat_out/AP027280.AP027280.tlosatx.n8.out"
    },
    {
        "name": "NZ_CP006932_Self",
        "mode": "TBLASTX",
        "ncbi":  "./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n8.out",
        "losat": "./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n8.out"
    },
    {
        "name": "AP027132_vs_NZ_CP006932",
        "mode": "TBLASTX",
        "ncbi":  "./blast_out/AP027132.NZ_CP006932.tblastx.n8.out",
        "losat": "./losat_out/AP027132.NZ_CP006932.tlosatx.n8.out"
    },
    {
        "name": "AP027078_vs_AP027131",
        "mode": "TBLASTX",
        "ncbi":  "./blast_out/AP027078.AP027131.tblastx.n8.out",
        "losat": "./losat_out/AP027078.AP027131.tlosatx.n8.out"
    },
    {
        "name": "AP027131_vs_AP027133",
        "mode": "TBLASTX",
        "ncbi":  "./blast_out/AP027131.AP027133.tblastx.n8.out",
        "losat": "./losat_out/AP027131.AP027133.tlosatx.n8.out"
    },
    {
        "name": "AP027133_vs_AP027132",
        "mode": "TBLASTX",
        "ncbi":  "./blast_out/AP027133.AP027132.tblastx.n8.out",
        "losat": "./losat_out/AP027133.AP027132.tlosatx.n8.out"
    },
    {
        "name": "MjeNMV_vs_MelaMJNV", "mode": "TBLASTX", 
        "ncbi": "./blast_out/MjeNMV.MelaMJNV.tblastx.n8.out", 
        "losat": "./losat_out/MjeNMV.MelaMJNV.tlosatx.n8.out"
    },
    {
        "name": "MjeNMV_vs_MelaMJNV", "mode": "TBLASTX", 
        "ncbi": "./blast_out/AP027280.AP027280.tblastx.n8.out", 
        "losat": "./losat_out/AP027280.AP027280.tlosatx.n8.out"
    },
    {
        "name": "MelaMJNV_vs_PemoMJNVA", "mode": "TBLASTX", 
        "ncbi": "./blast_out/MelaMJNV.PemoMJNVA.tblastx.n8.out", 
        "losat": "./losat_out/MelaMJNV.PemoMJNVA.tlosatx.n8.out"
    },
    {
        "name": "PemoMJNVA_vs_PeseMJNV", "mode": "TBLASTX", 
        "ncbi": "./blast_out/PemoMJNVA.PeseMJNV.tblastx.n8.out", 
        "losat": "./losat_out/PemoMJNVA.PeseMJNV.tlosatx.n8.out"
    },
    {
        "name": "PeseMJNV_vs_PemoMJNVB", "mode": "TBLASTX", 
        "ncbi": "./blast_out/PeseMJNV.PemoMJNVB.tblastx.n8.out", 
        "losat": "./losat_out/PeseMJNV.PemoMJNVB.tlosatx.n8.out"
    },
    {
        "name": "PemoMJNVB_vs_LvMJNV", "mode": "TBLASTX", 
        "ncbi": "./blast_out/PemoMJNVB.LvMJNV.tblastx.n8.out", 
        "losat": "./losat_out/PemoMJNVB.LvMJNV.tlosatx.n8.out"
    },
    {
        "name": "LvMJNV_vs_TrcuMJNV", "mode": "TBLASTX", 
        "ncbi": "./blast_out/LvMJNV.TrcuMJNV.tblastx.n8.out", 
        "losat": "./losat_out/LvMJNV.TrcuMJNV.tlosatx.n8.out"
    },
    {   "name": "TrcuMJNV_vs_MellatMJNV", "mode": "TBLASTX",
        "ncbi": "./blast_out/TrcuMJNV.MellatMJNV.tblastx.n8.out", 
        "losat": "./losat_out/TrcuMJNV.MellatMJNV.tlosatx.n8.out"
    },
    {
        "name": "MellatMJNV_vs_MeenMJNV", "mode": "TBLASTX",
        "ncbi": "./blast_out/MellatMJNV.MeenMJNV.tblastx.n8.out", 
        "losat": "./losat_out/MellatMJNV.MeenMJNV.tlosatx.n8.out"
    },
    {
        "name": "MeenMJNV_vs_MejoMJNV", "mode": "TBLASTX",
        "ncbi": "./blast_out/MeenMJNV.MejoMJNV.tblastx.n8.out", 
        "losat": "./losat_out/MeenMJNV.MejoMJNV.tlosatx.n8.out"
    },
    {
        "name": "AvCLPV_vs_PsCLPV", "mode": "TBLASTX",
        "ncbi": "./blast_out/AvCLPV.PsCLPV.tblastx.n8.out", 
        "losat": "./losat_out/AvCLPV.PsCLPV.tlosatx.n8.out"
    },
    # --- BLASTN (Standard / Task:blastn) ---
    {
        "name": "NZ_CP006932_Self",
        "mode": "BLASTN_Task",
        "ncbi":  "./blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out",
        "losat": "./losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.out"
    },
    {
        "name": "PesePMNV_vs_MjPMNV",
        "mode": "BLASTN",
        "ncbi":  "./blast_out/PesePMNV.MjPMNV.blastn.out",
        "losat": "./losat_out/PesePMNV.MjPMNV.losatn.blastn.out"
    },
    {
        "name": "MelaMJNV_vs_PemoMJNVA",
        "mode": "BLASTN",
        "ncbi":  "./blast_out/MelaMJNV.PemoMJNVA.blastn.out",
        "losat": "./losat_out/MelaMJNV.PemoMJNVA.losatn.blastn.out"
    },
    {
        "name": "SiNMV_vs_ChdeNMV",
        "mode": "BLASTN",
        "ncbi":  "./blast_out/SiNMV.ChdeNMV.blastn.out",
        "losat": "./losat_out/SiNMV.ChdeNMV.losatn.blastn.out"
    },
    {
        "name": "PmeNMV_vs_MjPMNV",
        "mode": "BLASTN",
        "ncbi":  "./blast_out/PmeNMV.MjPMNV.blastn.out",
        "losat": "./losat_out/PmeNMV.MjPMNV.losatn.blastn.out"
    },
    {
        "name": "PmeNMV_vs_PesePMNV",
        "mode": "BLASTN",
        "ncbi":  "./blast_out/PmeNMV.PesePMNV.blastn.out",
        "losat": "./losat_out/PmeNMV.PesePMNV.losatn.blastn.out"
    },
    {
        "name": "PeseMJNV_vs_PemoMJNVB",
        "mode": "BLASTN",
        "ncbi":  "./blast_out/PeseMJNV.PemoMJNVB.blastn.out",
        "losat": "./losat_out/PeseMJNV.PemoMJNVB.losatn.blastn.out"
    },
    {
        "name": "PemoMJNVA_vs_PeseMJNV",
        "mode": "BLASTN",
        "ncbi":  "./blast_out/PemoMJNVA.PeseMJNV.blastn.out",
        "losat": "./losat_out/PemoMJNVA.PeseMJNV.losatn.blastn.out"
    },
    {
        "name": "MjeNMV_vs_MelaMJNV",
        "mode": "BLASTN",
        "ncbi":  "./blast_out/MjeNMV.MelaMJNV.blastn.out",
        "losat": "./losat_out/MjeNMV.MelaMJNV.losatn.blastn.out"
    },
    {
        "name": "MjPMNV_vs_MlPMNV",
        "mode": "BLASTN",
        "ncbi":  "./blast_out/MjPMNV.MlPMNV.blastn.out",
        "losat": "./losat_out/MjPMNV.MlPMNV.losatn.blastn.out"
    },

    # --- BLASTN (Megablast) ---
    {
        "name": "EDL933_vs_Sakai",
        "mode": "Megablast",
        "ncbi":  "./blast_out/EDL933.Sakai.blastn.megablast.out",
        "losat": "./losat_out/EDL933.Sakai.losatn.megablast.out"
    },
    {
        "name": "Sakai_vs_MG1655",
        "mode": "Megablast",
        "ncbi":  "./blast_out/Sakai.MG1655.blastn.megablast.out",
        "losat": "./losat_out/Sakai.MG1655.losatn.megablast.out"
    },
]

def main():
    print("Starting batch plot generation...")
    print(f"Output directory: {PLOT_DIR}")
    
    for config in comparisons:
        generate_comparison_plot(config)
        
    print("All plots generated.")

if __name__ == "__main__":
    main()