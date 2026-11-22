- `losat windowmasker` with an `-ids` or `-exclude_ids` file matched as text: a line that gives an empty
  id (one that starts with a space or a tab, or a line `>`, `|`, `>|`, `>>` or `>>|`) gets NCBI's error
  `bad id: ...: ignoring` and is skipped, and the run goes on with NCBI BLAST+ 2.17.0's output for the
  file without that line. NCBI writes the same error and then crashes; this is a
  documented difference (the manual's "Differences from NCBI BLAST+" page).
