- Fixed: an `-evalue` of NaN (`-nan`, `+nan`, which NCBI BLAST+ accepts) gives NCBI's output in
  every program. BLASTP, PSIBLAST, RPSBLAST, TBLASTN and BLASTX without `-searchsp` now take NCBI's
  score cutoff of 0 for it (LOSAT computed the cutoff of the smallest double and wrote no hit), and
  TBLASTN, RPSBLAST and PHI-BLAST keep every HSP as NCBI does (LOSAT kept none); with the
  composition-based redo no match is kept, as in NCBI, in TBLASTN and BLASTX too (both kept every
  match), while RPSBLAST and RPSTBLASTN keep every match as NCBI's redo hands it back, with the
  scores scaled by 32 and bit scores of 0 (LOSAT kept none). BLASTX and RPSTBLASTN no longer
  refuse it. A subnormal `-evalue` (such as `1e-310`, which NCBI accepts) takes NCBI's score
  cutoff too (LOSAT raised it to the smallest normal double first) (line RG, X-evalue-nan).
- Changed: an infinite `-evalue` (`+inf`, `1e999`) runs without the composition-based redo
  (`-comp_based_stats 0`, `-matrix IDENTITY`, `-ungapped`) in BLASTP, PSIBLAST, RPSBLAST, TBLASTN,
  BLASTX and RPSTBLASTN and gives NCBI's output there; with the redo it stays refused, as DBL_MAX
  is, since NCBI crashes on some inputs, except in RPSBLAST and RPSTBLASTN, which now run both
  values with the redo as NCBI does. BLASTX now refuses DBL_MAX with the redo too (it ran it;
  NCBI's blastx crashes on some inputs there) (line RG, X-evalue-nan).
- Fixed: LOSAT Web's ABI v1 reads the `-evalue` of BLASTP and TBLASTX as the command line and
  NCBI BLAST+ do: `+inf`, `-nan`, `+nan(1)` and `1e999` are read, a bare `inf` or `nan` is an
  argument error (line RG, X-evalue-nan).
- Fixed: `-outfmt 15` (BLAST JSON) and `-outfmt 13` give NCBI's output where a number of the
  report is NaN or infinite (an `-evalue` of NaN or infinite; the statistics of an `-in_pssm`
  PSSM whose values are NaN), in every program that writes them: each report ends after the
  number's key with NCBI's two error posts and the run goes on to the document's end (exit status
  0); `-outfmt 13` ends its first side file there and the run with NCBI's error (exit status 255).
  LOSAT stopped with an error of its own (line RG, X-evalue-nan).
