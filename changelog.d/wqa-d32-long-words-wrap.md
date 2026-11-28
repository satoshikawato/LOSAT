- LOSAT Web: a long name or ID no longer makes the page wider than the window. This fixes the line
  under Run LOSAT (an MSA, a PSSM or a database name; at 320 px even `hemoglobin_alpha_family_seed_alignment.aln`),
  the usage text that an unknown option in Additional options shows (its line of 72 `=`), the
  Alignments summary, the HSP list's heading and the Dot plot's caption for a sequence ID such as
  `NC_000913.3_Escherichia_coli_K12_MG1655`, Format Archive's chosen file and WindowMasker's output
  tab and download button. They wrap inside their box, as other long names do.

- LOSAT Web: a "?" help that is open while the window gets narrower stays inside the window. It kept the place it
  had when it opened, and the help of an option at the right edge (Mask lower case letters, Max HSPs per subject,
  Max target sequences) made the page scroll sideways by 3 to 9 px.
