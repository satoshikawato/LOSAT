# Angle (b): fidelity of the transpile items
Work dir: /home/kawato/.cache/losat-web-gui-target/e2g-audit/b/ ; results FINDINGS.md and CASES.tsv (case_id, item, command, env, ncbi_exit, losat_exit, stdout_same, stderr_same, first_difference).
For each item T1–T14, R1–R3 and V1 (list and commits from build_inventory.py RESULTS and `git log --oneline 9c810a3d7..HEAD`), read the commit diff (`git -C <repo> show <sha>`) next to the cited NCBI source and judge whether the Rust is a faithful transpile (C integer widths and wrap, comparison operators, order of operations, sort stability, timing of messages, exact texts). Then exercise each item's branch against NCBI with inputs you build yourself (deterministic, seeded), in addition to the frozen ones under LOSAT/tests/fixtures/blastn_regression/:
- reuse the generators of docs/evidence/losat_web_e2f/audit_round2/{gen,cases,driver}.py (copy them to your work dir and change the hard-coded binary paths there; never edit the repository copy);
- include other organisms: the viral genomes in LOSAT/tests/fasta (e.g. MeenMJNV.fasta, MejoMJNV.fasta, PemoMJNVA.fasta, LvMJNV.fasta, AvCLPV.fasta, LC738874.fasta) and MG1655.fna, not only EDL933/Sakai;
- include IUPAC-heavy sequences (replace 5–20% of bases by R, Y, K, M, S, W, B, D, H, V, N);
- include BATCH_SIZE (100, 1000, 5000), CHUNK_SIZE (2000, 40000) and OVERLAP_CHUNK_SIZE (0, 50) runs on multi-query inputs, -task blastn with small -word_size (4–7), and multi-query files with palindromic queries (X + revcomp(X));
- for T8/T9 build multi-batch query files with titles that end in 20+ nucleotides and invalid (all-N) queries, with a scoring that has no Karlin-Altschul table (e.g. -reward 1 -penalty -6) and with a default one;
- for T14: -out /dev/full with -outfmt 0;
- for R1–R3 and the T7/T12 rejections: confirm LOSAT rejects explicitly and that NCBI indeed behaves as the rejection message says.
Aim for at least 400 compared cases in total. Verdict per item and overall.
