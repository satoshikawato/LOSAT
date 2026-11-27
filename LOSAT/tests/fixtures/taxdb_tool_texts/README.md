# Texts of dustmasker, segmasker and blastdbcmd on a damaged taxonomy database (line RG item T5b)

`<case>.txt` is the merged standard output and standard error (`2>&1`) of NCBI BLAST+ 2.17.0's `dustmasker`, `segmasker` or
`blastdbcmd` for the case, and `<case>.log` its `-logfile` (see `tests/taxdb_tool_texts.rs`). The FASTA letters in the `.txt`
files are those of the databases `tests/db_search/t5b/` that these programs read; `tests/db_search/t5b/make_t5b.py` writes
those databases from the letters of `tests/db_search/ids/ndup.fna` and `ids/pdup.faa`. Added in commit 937985708
(2026-11-14).

`pipe.*`, `split.*` and `closed.*` (line RG item T8a, added 2026-11-25) are NCBI BLAST+ 2.17.0's bytes of the cases of
`PIPE_CASES` in `tests/taxdb_tool_texts.rs`, each run once on Linux with `BLAST_USAGE_REPORT=false`, `HOME` and
`NCBI_CONFIG_PATH` an empty directory and `BLASTDB` naming `tests/seqdb/taxdb/<variant>`: `pipe.<case>.txt` is what was
read from one pipe that was standard output and standard error, with `-out /dev/stdout`; `split.<case>.out` and `.err`
standard output and standard error as two pipes, without `-out`; `closed.<case>.err` standard error where standard output
was a pipe whose reader had closed before the start, with `-out /dev/stdout` (NCBI ended by SIGPIPE, exit status 141).
The inputs are those of `tests/db_search/t5b/` and of `tests/db_search/t8/` (`make_t8.py` there, letters of its own
generator and of `tests/db_search/ids/txr_q.fna`). The test maps them to LOSAT's words (`losat::wording`).
