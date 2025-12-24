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

# === Column Definitions ===
COLUMNS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
           'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

# === Comparison List (Same as before) ===
comparisons = [
    # --- TBLASTX ---
    {"name": "NZ_CP006932_Self", "mode": "TBLASTX", "ncbi": "./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out", "losat": "./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.out"},
    {"name": "AP027132_vs_NZ", "mode": "TBLASTX", "ncbi": "./blast_out/AP027132.NZ_CP006932.tblastx.n1.out", "losat": "./losat_out/AP027132.NZ_CP006932.tlosatx.n1.out"},
    {"name": "AP027078_vs_AP027131", "mode": "TBLASTX", "ncbi": "./blast_out/AP027078.AP027131.tblastx.n1.out", "losat": "./losat_out/AP027078.AP027131.tlosatx.n1.out"},
    {"name": "AP027131_vs_AP027133", "mode": "TBLASTX", "ncbi": "./blast_out/AP027131.AP027133.tblastx.n1.out", "losat": "./losat_out/AP027131.AP027133.tlosatx.n1.out"},
    {"name": "AP027133_vs_AP027132", "mode": "TBLASTX", "ncbi": "./blast_out/AP027133.AP027132.tblastx.n1.out", "losat": "./losat_out/AP027133.AP027132.tlosatx.n1.out"},

    # --- BLASTN / Megablast ---
    {"name": "NZ_CP006932_Self(Mega)", "mode": "BLASTN (Megablast)", "ncbi": "./blast_out/NZ_CP006932.NZ_CP006932.blastn.out", "losat": "./losat_out/NZ_CP006932.NZ_CP006932.losatn.out"},
    {"name": "EDL933_vs_Sakai", "mode": "BLASTN (Megablast)", "ncbi": "./blast_out/EDL933.Sakai.blastn.megablast.out", "losat": "./losat_out/EDL933.Sakai.losatn.megablast.out"},
    {"name": "Sakai_vs_MG1655", "mode": "BLASTN (Megablast)", "ncbi": "./blast_out/Sakai.MG1655.blastn.megablast.out", "losat": "./losat_out/Sakai.MG1655.losatn.megablast.out"},

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
        df['Mode'] = mode_name # e.g., "TBLASTX"
        df['Task'] = task_name # e.g., "NZ_CP006932_Self"
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
    # We want to separate "TBLASTX" vs "BLASTN" because length/scores are different units/scales
    # We will treat "BLASTN (Megablast)" and "BLASTN (Task:blastn)" as basically the same category "Nucleotide"
    # for visualization scaling, or keep them distinct if distributions vary wildly.
    # Let's simplify into 2 main groups: Protein-level (TBLASTX) and Nucleotide-level (BLASTN variants)
    
    df_all['Broad_Mode'] = df_all['Mode'].apply(lambda x: 'TBLASTX' if 'TBLASTX' in x else 'BLASTN (All Types)')

    sns.set(style="whitegrid")
    
    # Create a FacetGrid-like figure manually with subplots
    # Rows: Broad Mode (TBLASTX vs BLASTN)
    # Cols: Length Dist, Identity Dist, Scatter (Len vs Ident)
    
    modes = df_all['Broad_Mode'].unique()
    fig, axes = plt.subplots(len(modes), 3, figsize=(18, 6 * len(modes)))
    
    # Handle case if only 1 mode exists (axes is 1D array)
    if len(modes) == 1:
        axes = [axes]

    for i, mode in enumerate(modes):
        data_subset = df_all[df_all['Broad_Mode'] == mode]
        
        # Row Title
        axes[i][0].text(-0.2, 0.5, mode, transform=axes[i][0].transAxes, 
                        fontsize=16, rotation=90, va='center', fontweight='bold')

        # 1. Alignment Length Distribution
        sns.histplot(data=data_subset, x='length', hue='Tool', 
                     element="step", stat="density", common_norm=False, 
                     log_scale=True, ax=axes[i][0])
        axes[i][0].set_title(f'Alignment Length Distribution')
        axes[i][0].set_xlabel('Length (bp or aa)')

        # 2. Identity Distribution
        sns.histplot(data=data_subset, x='pident', hue='Tool', 
                     element="step", stat="density", common_norm=False, 
                     ax=axes[i][1])
        axes[i][1].set_title(f'Percentage Identity Distribution')
        axes[i][1].set_xlabel('Identity (%)')

        # 3. Scatter: Length vs Identity
        # Sampling points if too many to prevent overplotting and slow rendering
        if len(data_subset) > 5000:
            plot_data = data_subset.sample(5000, random_state=42)
            title_suffix = "(Subsampled 5k)"
        else:
            plot_data = data_subset
            title_suffix = ""
            
        sns.scatterplot(data=plot_data, x='length', y='pident', hue='Tool', 
                        style='Tool', alpha=0.5, ax=axes[i][2])
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