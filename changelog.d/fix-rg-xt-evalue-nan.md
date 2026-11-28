- Fixed: an `-evalue` of NaN (`-nan`, `+nan`, which NCBI BLAST+ accepts) gives NCBI's output in
  every program. BLASTP, PSIBLAST, RPSBLAST, TBLASTN and BLASTX without `-searchsp` now take NCBI's
  score cutoff of 0 for it (LOSAT computed the cutoff of the smallest double and wrote no hit), and
  TBLASTN, RPSBLAST and PHI-BLAST keep every HSP as NCBI does (LOSAT kept none); with the
  composition-based redo no match is kept, as in NCBI, in TBLASTN and BLASTX too (both kept every
  match). BLASTX and RPSTBLASTN no longer refuse it (line RG, X-evalue-nan).
- Changed: an infinite `-evalue` (`+inf`, `1e999`) runs without the composition-based redo
  (`-comp_based_stats 0`, `-matrix IDENTITY`, `-ungapped`) in BLASTP, PSIBLAST, RPSBLAST, TBLASTN,
  BLASTX and RPSTBLASTN and gives NCBI's output there; with the redo it stays refused, as DBL_MAX
  is, since NCBI crashes on some inputs (line RG, X-evalue-nan).
