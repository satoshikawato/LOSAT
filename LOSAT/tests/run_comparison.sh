#!/bin/bash

# --- Setup ---

# Create output directories
mkdir -p losat_out blast_out

# Path to LOSAT binary
LOSAT_BIN="../target/release/LOSAT"

# ==========================================
# Part 1: LOSAT Execution
# ==========================================
echo "Starting LOSAT commands..."


# --- TBLASTX Commands (Genetic Code: 4) ---

# NZ_CP006932 self
(time $LOSAT_BIN tblastx -q ./fasta/NZ_CP006932.fasta -s ./fasta/NZ_CP006932.fasta -o ./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1 )&>./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.log

# AP027132 vs NZ_CP006932
(time $LOSAT_BIN tblastx -q ./fasta/AP027132.fasta -s ./fasta/NZ_CP006932.fasta -o ./losat_out/AP027132.NZ_CP006932.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1 )&>./losat_out/AP027132.NZ_CP006932.tlosatx.n1.log

# AP027078 vs AP027131
(time $LOSAT_BIN tblastx -q ./fasta/AP027078.fasta -s ./fasta/AP027131.fasta -o ./losat_out/AP027078.AP027131.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1 )&>./losat_out/AP027078.AP027131.tlosatx.n1.log

# AP027131 vs AP027133
(time $LOSAT_BIN tblastx -q ./fasta/AP027131.fasta -s ./fasta/AP027133.fasta -o ./losat_out/AP027131.AP027133.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1 )&>./losat_out/AP027131.AP027133.tlosatx.n1.log

# AP027133 vs AP027132
(time $LOSAT_BIN tblastx -q ./fasta/AP027133.fasta -s ./fasta/AP027132.fasta -o ./losat_out/AP027133.AP027132.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1 )&>./losat_out/AP027133.AP027132.tlosatx.n1.log



# --- BLASTN Commands (Default / Megablast behavior) ---

# NZ_CP006932 self (Default/Megablast)
(time $LOSAT_BIN blastn -q ./fasta/NZ_CP006932.fasta -s ./fasta/NZ_CP006932.fasta -o ./losat_out/NZ_CP006932.NZ_CP006932.blastn.megablast.out -n 1 )&>./losat_out/NZ_CP006932.NZ_CP006932.blastn.megablast.log

# NZ_CP006932 self (Default, legacy name for compatibility)
(time $LOSAT_BIN blastn -q ./fasta/NZ_CP006932.fasta -s ./fasta/NZ_CP006932.fasta -o ./losat_out/NZ_CP006932.NZ_CP006932.losatn.out -n 1 )&>./losat_out/NZ_CP006932.NZ_CP006932.losatn.log

# EDL933 vs Sakai
(time $LOSAT_BIN blastn -q ./fasta/EDL933.fna -s ./fasta/Sakai.fna -o ./losat_out/EDL933.Sakai.blastn.megablast.out -n 1 )&>./losat_out/EDL933.Sakai.blastn.megablast.log

# Sakai vs MG1655
(time $LOSAT_BIN blastn -q ./fasta/Sakai.fna -s ./fasta/MG1655.fna -o ./losat_out/Sakai.MG1655.blastn.megablast.out -n 1 )&>./losat_out/Sakai.MG1655.blastn.megablast.log


# --- BLASTN Commands (Task: blastn) ---

# NZ_CP006932 self (Task: blastn)
(time $LOSAT_BIN blastn -q ./fasta/NZ_CP006932.fasta -s ./fasta/NZ_CP006932.fasta -o ./losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.out --task blastn -n 1 )&>./losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.log

# PesePMNV vs MjPMNV
(time $LOSAT_BIN blastn -q ./fasta/AP027152.fasta -s ./fasta/AP027202.fasta -o ./losat_out/PesePMNV.MjPMNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PesePMNV.MjPMNV.losatn.blastn.log

# MelaMJNV vs PemoMJNVA
(time $LOSAT_BIN blastn -q ./fasta/LC738874.fasta -s ./fasta/LC738870.fasta -o ./losat_out/MelaMJNV.PemoMJNVA.losatn.blastn.out --task blastn -n 1 )&>./losat_out/MelaMJNV.PemoMJNVA.losatn.blastn.log

# SiNMV vs ChdeNMV
(time $LOSAT_BIN blastn -q ./fasta/LC738884.fasta -s ./fasta/AP027155.fasta -o ./losat_out/SiNMV.ChdeNMV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/SiNMV.ChdeNMV.losatn.blastn.log

# PmeNMV vs MjPMNV
(time $LOSAT_BIN blastn -q ./fasta/LC738869.fasta -s ./fasta/AP027202.fasta -o ./losat_out/PmeNMV.MjPMNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PmeNMV.MjPMNV.losatn.blastn.log

# PmeNMV vs PesePMNV
(time $LOSAT_BIN blastn -q ./fasta/LC738869.fasta -s ./fasta/AP027152.fasta -o ./losat_out/PmeNMV.PesePMNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PmeNMV.PesePMNV.losatn.blastn.log

# PeseMJNV vs PemoMJNVB
(time $LOSAT_BIN blastn -q ./fasta/LC738873.fasta -s ./fasta/LC738871.fasta -o ./losat_out/PeseMJNV.PemoMJNVB.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PeseMJNV.PemoMJNVB.losatn.blastn.log

# PemoMJNVA vs PeseMJNV
(time $LOSAT_BIN blastn -q ./fasta/LC738870.fasta -s ./fasta/LC738873.fasta -o ./losat_out/PemoMJNVA.PeseMJNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/PemoMJNVA.PeseMJNV.losatn.blastn.log

# MjeNMV vs MelaMJNV
(time $LOSAT_BIN blastn -q ./fasta/LC738868.fasta -s ./fasta/LC738874.fasta -o ./losat_out/MjeNMV.MelaMJNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/MjeNMV.MelaMJNV.losatn.blastn.log

# MjPMNV vs MlPMNV
(time $LOSAT_BIN blastn -q ./fasta/AP027202.fasta -s ./fasta/LC738875.fasta -o ./losat_out/MjPMNV.MlPMNV.losatn.blastn.out --task blastn -n 1 )&>./losat_out/MjPMNV.MlPMNV.losatn.blastn.log



<<COMMENTOUT
# ==========================================
# Part 2: BLAST+ Execution
# ==========================================
echo "Starting BLAST+ commands..."

# --- TBLASTX Commands (gencode 4) ---

# NZ_CP006932 self
(time tblastx -query ./fasta/NZ_CP006932.fasta -subject ./fasta/NZ_CP006932.fasta -out ./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out -query_gencode 4 -db_gencode 4 -num_threads 32 -outfmt 7 )&>./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.log

# AP027132 vs NZ_CP006932
(time tblastx -query ./fasta/AP027132.fasta -subject ./fasta/NZ_CP006932.fasta -out ./blast_out/AP027132.NZ_CP006932.tblastx.n1.out -query_gencode 4 -db_gencode 4 -num_threads 32 -outfmt 7 )&>./blast_out/AP027132.NZ_CP006932.tblastx.n1.log

# AP027078 vs AP027131
(time tblastx -query ./fasta/AP027078.fasta -subject ./fasta/AP027131.fasta -out ./blast_out/AP027078.AP027131.tblastx.n1.out -query_gencode 4 -db_gencode 4 -num_threads 32 -outfmt 7 )&>./blast_out/AP027078.AP027131.tblastx.n1.log

# AP027131 vs AP027133
(time tblastx -query ./fasta/AP027131.fasta -subject ./fasta/AP027133.fasta -out ./blast_out/AP027131.AP027133.tblastx.n1.out -query_gencode 4 -db_gencode 4 -num_threads 32 -outfmt 7 )&>./blast_out/AP027131.AP027133.tblastx.n1.log

# AP027133 vs AP027132
(time tblastx -query ./fasta/AP027133.fasta -subject ./fasta/AP027132.fasta -out ./blast_out/AP027133.AP027132.tblastx.n1.out -query_gencode 4 -db_gencode 4 -num_threads 32 -outfmt 7 )&>./blast_out/AP027133.AP027132.tblastx.n1.log


# --- BLASTN Commands (Default / Megablast) ---

# NZ_CP006932 self (Default/Megablast)
(time blastn -query ./fasta/NZ_CP006932.fasta -subject ./fasta/NZ_CP006932.fasta -out ./blast_out/NZ_CP006932.NZ_CP006932.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/NZ_CP006932.NZ_CP006932.blastn.log

# EDL933 vs Sakai
(time blastn -query ./fasta/EDL933.fna -subject ./fasta/Sakai.fna -outfmt 7 -out ./blast_out/EDL933.Sakai.blastn.megablast.out )&>./blast_out/EDL933.Sakai.blastn.megablast.log

# Sakai vs MG1655
(time blastn -query ./fasta/Sakai.fna -subject ./fasta/MG1655.fna -outfmt 7 -out ./blast_out/Sakai.MG1655.blastn.megablast.out )&>./blast_out/Sakai.MG1655.blastn.megablast.log


# --- BLASTN Commands (Task: blastn) ---

# NZ_CP006932 self (Task: blastn)
(time blastn -task blastn -query ./fasta/NZ_CP006932.fasta -subject ./fasta/NZ_CP006932.fasta -out ./blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/NZ_CP006932.NZ_CP006932.task_blastn.log

# PesePMNV vs MjPMNV
(time blastn -task blastn -query ./fasta/AP027152.fasta -subject ./fasta/AP027202.fasta -out ./blast_out/PesePMNV.MjPMNV.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/PesePMNV.MjPMNV.blastn.log

# MelaMJNV vs PemoMJNVA
(time blastn -task blastn -query ./fasta/LC738874.fasta -subject ./fasta/LC738870.fasta -out ./blast_out/MelaMJNV.PemoMJNVA.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/MelaMJNV.PemoMJNVA.blastn.log

# SiNMV vs ChdeNMV
(time blastn -task blastn -query ./fasta/LC738884.fasta -subject ./fasta/AP027155.fasta -out ./blast_out/SiNMV.ChdeNMV.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/SiNMV.ChdeNMV.blastn.log

# PmeNMV vs MjPMNV
(time blastn -task blastn -query ./fasta/LC738869.fasta -subject ./fasta/AP027202.fasta -out ./blast_out/PmeNMV.MjPMNV.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/PmeNMV.MjPMNV.blastn.log

# PmeNMV vs PesePMNV
(time blastn -task blastn -query ./fasta/LC738869.fasta -subject ./fasta/AP027152.fasta -out ./blast_out/PmeNMV.PesePMNV.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/PmeNMV.PesePMNV.blastn.log

# PeseMJNV vs PemoMJNVB
(time blastn -task blastn -query ./fasta/LC738873.fasta -subject ./fasta/LC738871.fasta -out ./blast_out/PeseMJNV.PemoMJNVB.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/PeseMJNV.PemoMJNVB.blastn.log

# PemoMJNVA vs PeseMJNV
(time blastn -task blastn -query ./fasta/LC738870.fasta -subject ./fasta/LC738873.fasta -out ./blast_out/PemoMJNVA.PeseMJNV.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/PemoMJNVA.PeseMJNV.blastn.log

# MjeNMV vs MelaMJNV
(time blastn -task blastn -query ./fasta/LC738868.fasta -subject ./fasta/LC738874.fasta -out ./blast_out/MjeNMV.MelaMJNV.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/MjeNMV.MelaMJNV.blastn.log

# MjPMNV vs MlPMNV
(time blastn -task blastn -query ./fasta/AP027202.fasta -subject ./fasta/LC738875.fasta -out ./blast_out/MjPMNV.MlPMNV.blastn.out -num_threads 1 -outfmt 7 )&>./blast_out/MjPMNV.MlPMNV.blastn.log
COMMENTOUT
echo "All done."
