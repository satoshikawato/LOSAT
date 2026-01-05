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

# --- TLOSATX Commands (Genetic Code: 1) ---
echo "Starting AP027280 self..."
(time $LOSAT_BIN tblastx -q ./fasta/AP027280.fasta -s ./fasta/AP027280.fasta -o ./losat_out/AP027280.AP027280.tlosatx.n1.out --query-gencode 1 --db-gencode 1 -n 1)&>./losat_out/AP027280.AP027280.tlosatx.n1.log 
(time $LOSAT_BIN tblastx -q ./fasta/AP027280.fasta -s ./fasta/AP027280.fasta -o ./losat_out/AP027280.AP027280.tlosatx.n8.out --query-gencode 1 --db-gencode 1 -n 8 )&>./losat_out/AP027280.AP027280.tlosatx.n8.log 

# MjeNMV vs MelaMJNV
echo "Starting MjeNMV vs MelaMJNV (LOSAT)..."
(time $LOSAT_BIN tblastx -q ./fasta/MjeNMV.fasta -s ./fasta/MelaMJNV.fasta -o ./losat_out/MjeNMV.MelaMJNV.tlosatx.n1.out --query-gencode 1 --db-gencode 1 -n 1 ) &> ./losat_out/MjeNMV.MelaMJNV.tlosatx.n1.log 

# MelaMJNV vs PemoMJNVA
echo "Starting MelaMJNV vs PemoMJNVA (LOSAT)..."
(time $LOSAT_BIN tblastx -q ./fasta/MelaMJNV.fasta -s ./fasta/PemoMJNVA.fasta -o ./losat_out/MelaMJNV.PemoMJNVA.tlosatx.n1.out --query-gencode 1 --db-gencode 1 -n 1 ) &> ./losat_out/MelaMJNV.PemoMJNVA.tlosatx.n1.log 

# PemoMJNVA vs PeseMJNV
echo "Starting PemoMJNVA vs PeseMJNV (LOSAT)..."
(time $LOSAT_BIN tblastx -q ./fasta/PemoMJNVA.fasta -s ./fasta/PeseMJNV.fasta -o ./losat_out/PemoMJNVA.PeseMJNV.tlosatx.n1.out --query-gencode 1 --db-gencode 1 -n 1 ) &> ./losat_out/PemoMJNVA.PeseMJNV.tlosatx.n1.log 

# PeseMJNV vs PemoMJNVB
echo "Starting PeseMJNV vs PemoMJNVB (LOSAT)..."
(time $LOSAT_BIN tblastx -q ./fasta/PeseMJNV.fasta -s ./fasta/PemoMJNVB.fasta -o ./losat_out/PeseMJNV.PemoMJNVB.tlosatx.n1.out --query-gencode 1 --db-gencode 1 -n 1 ) &> ./losat_out/PeseMJNV.PemoMJNVB.tlosatx.n1.log 

# PemoMJNVB vs LvMJNV
echo "Starting PemoMJNVB vs LvMJNV (LOSAT)..."
(time $LOSAT_BIN tblastx -q ./fasta/PemoMJNVB.fasta -s ./fasta/LvMJNV.fasta -o ./losat_out/PemoMJNVB.LvMJNV.tlosatx.n1.out --query-gencode 1 --db-gencode 1 -n 1 ) &> ./losat_out/PemoMJNVB.LvMJNV.tlosatx.n1.log 

# LvMJNV vs TrcuMJNV
echo "Starting LvMJNV vs TrcuMJNV (LOSAT)..."
(time $LOSAT_BIN tblastx -q ./fasta/LvMJNV.fasta -s ./fasta/TrcuMJNV.fasta -o ./losat_out/LvMJNV.TrcuMJNV.tlosatx.n1.out --query-gencode 1 --db-gencode 1 -n 1 ) &> ./losat_out/LvMJNV.TrcuMJNV.tlosatx.n1.log 

# TrcuMJNV vs MellatMJNV
echo "Starting TrcuMJNV vs MellatMJNV (LOSAT)..."
(time $LOSAT_BIN tblastx -q ./fasta/TrcuMJNV.fasta -s ./fasta/MellatMJNV.fasta -o ./losat_out/TrcuMJNV.MellatMJNV.tlosatx.n1.out --query-gencode 1 --db-gencode 1 -n 1 ) &> ./losat_out/TrcuMJNV.MellatMJNV.tlosatx.n1.log 

# MellatMJNV vs MeenMJNV
echo "Starting MellatMJNV vs MeenMJNV (LOSAT)..."
(time $LOSAT_BIN tblastx -q ./fasta/MellatMJNV.fasta -s ./fasta/MeenMJNV.fasta -o ./losat_out/MellatMJNV.MeenMJNV.tlosatx.n1.out --query-gencode 1 --db-gencode 1 -n 1 ) &> ./losat_out/MellatMJNV.MeenMJNV.tlosatx.n1.log 

# MeenMJNV vs MejoMJNV
echo "Starting MeenMJNV vs MejoMJNV (LOSAT)..."
(time $LOSAT_BIN tblastx -q ./fasta/MeenMJNV.fasta -s ./fasta/MejoMJNV.fasta -o ./losat_out/MeenMJNV.MejoMJNV.tlosatx.n1.out --query-gencode 1 --db-gencode 1 -n 1 ) &> ./losat_out/MeenMJNV.MejoMJNV.tlosatx.n1.log 

# AvCLPV vs PsCLPV
echo "Starting AvCLPV vs PsCLPV (LOSAT)..."
(time $LOSAT_BIN tblastx -q ./fasta/AvCLPV.fasta -s ./fasta/PsCLPV.fasta -o ./losat_out/AvCLPV.PsCLPV.tlosatx.n1.out --query-gencode 1 --db-gencode 1 -n 1 )&>./losat_out/AvCLPV.PsCLPV.tlosatx.n1.log 

# --- TLOSATX Commands (Genetic Code: 4) ---
echo "Starting NZ_CP006932 self..."
# NZ_CP006932 self
(time $LOSAT_BIN tblastx -q ./fasta/NZ_CP006932.fasta -s ./fasta/NZ_CP006932.fasta -o ./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1 )&>./losat_out/NZ_CP006932.NZ_CP006932.tlosatx.n1.log 
echo "Starting AP027132 vs NZ_CP006932..."
# AP027132 vs NZ_CP006932
(time $LOSAT_BIN tblastx -q ./fasta/AP027132.fasta -s ./fasta/NZ_CP006932.fasta -o ./losat_out/AP027132.NZ_CP006932.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1 )&>./losat_out/AP027132.NZ_CP006932.tlosatx.n1.log 
echo "Starting AP027078 vs AP027131..."
# AP027078 vs AP027131
(time $LOSAT_BIN tblastx -q ./fasta/AP027078.fasta -s ./fasta/AP027131.fasta -o ./losat_out/AP027078.AP027131.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1 )&>./losat_out/AP027078.AP027131.tlosatx.n1.log 
echo "Starting AP027131 vs AP027133..."
# AP027131 vs AP027133
(time $LOSAT_BIN tblastx -q ./fasta/AP027131.fasta -s ./fasta/AP027133.fasta -o ./losat_out/AP027131.AP027133.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1 )&>./losat_out/AP027131.AP027133.tlosatx.n1.log 
echo "Starting AP027133 vs AP027132..."
# AP027133 vs AP027132
(time $LOSAT_BIN tblastx -q ./fasta/AP027133.fasta -s ./fasta/AP027132.fasta -o ./losat_out/AP027133.AP027132.tlosatx.n1.out --query-gencode 4 --db-gencode 4 -n 1 )&>./losat_out/AP027133.AP027132.tlosatx.n1.log 



<<COMMENTOUT
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
COMMENTOUT
echo "Finished LOSAT commands!"

<<COMMENTOUT
# ==========================================
# Part 2: BLAST+ Execution
# ==========================================
echo "Starting BLAST+ commands..."


# --- TBLASTX Commands (gencode 1) ---
# AP027280 self
echo "Starting AP027280_Self (BLAST)..."
(time tblastx -query ./fasta/AP027280.fasta -subject ./fasta/AP027280.fasta -out ./blast_out/AP027280.AP027280.tblastx.n1.out -query_gencode 1 -db_gencode 1 -num_threads 1 -outfmt 7 )&>./blast_out/AP027280.AP027280.tblastx.n1.log
(time tblastx -query ./fasta/AP027280.fasta -subject ./fasta/AP027280.fasta -out ./blast_out/AP027280.AP027280.tblastx.windowsize0.n1.out -query_gencode 1 -db_gencode 1 -num_threads 1 -outfmt 7 -window_size 0)&>./blast_out/AP027280.AP027280.tblastx.windowsize0.n1.log

# MjeNMV vs MelaMJNV
echo "Starting MjeNMV vs MelaMJNV (BLAST)..."
(time tblastx -query ./fasta/MjeNMV.fasta -subject ./fasta/MelaMJNV.fasta -out ./blast_out/MjeNMV.MelaMJNV.tblastx.n1.out -query_gencode 1 -db_gencode 1 -num_threads 1 -outfmt 7) &> ./blast_out/MjeNMV.MelaMJNV.tblastx.n1.log

# MelaMJNV vs PemoMJNVA
echo "Starting MelaMJNV vs PemoMJNVA (BLAST)..."
(time tblastx -query ./fasta/MelaMJNV.fasta -subject ./fasta/PemoMJNVA.fasta -out ./blast_out/MelaMJNV.PemoMJNVA.tblastx.n1.out -query_gencode 1 -db_gencode 1 -num_threads 1 -outfmt 7) &> ./blast_out/MelaMJNV.PemoMJNVA.tblastx.n1.log

# PemoMJNVA vs PeseMJNV
echo "Starting PemoMJNVA vs PeseMJNV (BLAST)..."
(time tblastx -query ./fasta/PemoMJNVA.fasta -subject ./fasta/PeseMJNV.fasta -out ./blast_out/PemoMJNVA.PeseMJNV.tblastx.n1.out -query_gencode 1 -db_gencode 1 -num_threads 1 -outfmt 7) &> ./blast_out/PemoMJNVA.PeseMJNV.tblastx.n1.log

# PeseMJNV vs PemoMJNVB
echo "Starting PeseMJNV vs PemoMJNVB (BLAST)..."
(time tblastx -query ./fasta/PeseMJNV.fasta -subject ./fasta/PemoMJNVB.fasta -out ./blast_out/PeseMJNV.PemoMJNVB.tblastx.n1.out -query_gencode 1 -db_gencode 1 -num_threads 1 -outfmt 7) &> ./blast_out/PeseMJNV.PemoMJNVB.tblastx.n1.log

# PemoMJNVB vs LvMJNV
echo "Starting PemoMJNVB vs LvMJNV (BLAST)..."
(time tblastx -query ./fasta/PemoMJNVB.fasta -subject ./fasta/LvMJNV.fasta -out ./blast_out/PemoMJNVB.LvMJNV.tblastx.n1.out -query_gencode 1 -db_gencode 1 -num_threads 1 -outfmt 7) &> ./blast_out/PemoMJNVB.LvMJNV.tblastx.n1.log

# LvMJNV vs TrcuMJNV
echo "Starting LvMJNV vs TrcuMJNV (BLAST)..."
(time tblastx -query ./fasta/LvMJNV.fasta -subject ./fasta/TrcuMJNV.fasta -out ./blast_out/LvMJNV.TrcuMJNV.tblastx.n1.out -query_gencode 1 -db_gencode 1 -num_threads 1 -outfmt 7) &> ./blast_out/LvMJNV.TrcuMJNV.tblastx.n1.log

# TrcuMJNV vs MellatMJNV
echo "Starting TrcuMJNV vs MellatMJNV (BLAST)..."
(time tblastx -query ./fasta/TrcuMJNV.fasta -subject ./fasta/MellatMJNV.fasta -out ./blast_out/TrcuMJNV.MellatMJNV.tblastx.n1.out -query_gencode 1 -db_gencode 1 -num_threads 1 -outfmt 7) &> ./blast_out/TrcuMJNV.MellatMJNV.tblastx.n1.log

# MellatMJNV vs MeenMJNV
echo "Starting MellatMJNV vs MeenMJNV (BLAST)..."
(time tblastx -query ./fasta/MellatMJNV.fasta -subject ./fasta/MeenMJNV.fasta -out ./blast_out/MellatMJNV.MeenMJNV.tblastx.n1.out -query_gencode 1 -db_gencode 1 -num_threads 1 -outfmt 7) &> ./blast_out/MellatMJNV.MeenMJNV.tblastx.n1.log

# MeenMJNV vs MejoMJNV
echo "Starting MeenMJNV vs MejoMJNV (BLAST)..."
(time tblastx -query ./fasta/MeenMJNV.fasta -subject ./fasta/MejoMJNV.fasta -out ./blast_out/MeenMJNV.MejoMJNV.tblastx.n1.out -query_gencode 1 -db_gencode 1 -num_threads 1 -outfmt 7) &> ./blast_out/MeenMJNV.MejoMJNV.tblastx.n1.log

# AvCLPV vs PsCLPV
echo "Starting AvCLPV vs PsCLPV (BLAST)..."
(time tblastx -query ./fasta/AvCLPV.fasta -subject ./fasta/PsCLPV.fasta -out ./blast_out/AvCLPV.PsCLPV.tblastx.n1.out -query_gencode 1 -db_gencode 1 -num_threads 1 -outfmt 7) &> ./blast_out/AvCLPV.PsCLPV.tblastx.n1.log



# --- TBLASTX Commands (gencode 4) ---

# NZ_CP006932 self
(time tblastx -query ./fasta/NZ_CP006932.fasta -subject ./fasta/NZ_CP006932.fasta -out ./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.out -query_gencode 4 -db_gencode 4 -num_threads 1 -outfmt 7 )&>./blast_out/NZ_CP006932.NZ_CP006932.tblastx.n1.log

# AP027132 vs NZ_CP006932
(time tblastx -query ./fasta/AP027132.fasta -subject ./fasta/NZ_CP006932.fasta -out ./blast_out/AP027132.NZ_CP006932.tblastx.n1.out -query_gencode 4 -db_gencode 4 -num_threads 1 -outfmt 7 )&>./blast_out/AP027132.NZ_CP006932.tblastx.n1.log

# AP027078 vs AP027131
(time tblastx -query ./fasta/AP027078.fasta -subject ./fasta/AP027131.fasta -out ./blast_out/AP027078.AP027131.tblastx.n1.out -query_gencode 4 -db_gencode 4 -num_threads 1 -outfmt 7 )&>./blast_out/AP027078.AP027131.tblastx.n1.log

# AP027131 vs AP027133
(time tblastx -query ./fasta/AP027131.fasta -subject ./fasta/AP027133.fasta -out ./blast_out/AP027131.AP027133.tblastx.n1.out -query_gencode 4 -db_gencode 4 -num_threads 1 -outfmt 7 )&>./blast_out/AP027131.AP027133.tblastx.n1.log

# AP027133 vs AP027132
(time tblastx -query ./fasta/AP027133.fasta -subject ./fasta/AP027132.fasta -out ./blast_out/AP027133.AP027132.tblastx.n1.out -query_gencode 4 -db_gencode 4 -num_threads 1 -outfmt 7 )&>./blast_out/AP027133.AP027132.tblastx.n1.log


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

wait
echo "All done."
