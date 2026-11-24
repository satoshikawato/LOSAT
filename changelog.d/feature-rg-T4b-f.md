- BLASTP, PSIBLAST, RPSBLAST and TBLASTN write the start of their outfmt 0 to 4 report (the
  program line, the references and the `Database:` lines) before they read the first query batch,
  where BLAST+ writes it, as BLASTN, BLASTX, TBLASTX and RPSTBLASTN do. With `-html` in a folder
  whose `./.ncbirc` is malformed, BLASTP and TBLASTN with `-export_search_strategy -` now end as
  BLAST+ does and as the run without the option: the first 159 (BLASTP) or 160 (TBLASTN) bytes of
  the HTML report, the registry's error and exit status 255. Before, they printed LOSAT's refusal of
  `-` with exit status 1 and nothing on standard output. A refusal that comes after that point (for
  example a database title or a query ID that BLAST+ shows wrongly) now follows the report's start on
  standard output, as BLAST+ writes it before the search. With `-out` on a full disk they now stop
  at once, without running the search first.
