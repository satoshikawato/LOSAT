# Texts of dustmasker, segmasker and blastdbcmd on a damaged taxonomy database (line RG item T5b)

`<case>.txt` is the merged standard output and standard error (`2>&1`) of NCBI BLAST+ 2.17.0's `dustmasker`, `segmasker` or
`blastdbcmd` for the case, and `<case>.log` its `-logfile` (see `tests/taxdb_tool_texts.rs`). The FASTA letters in the `.txt`
files are those of the databases `tests/db_search/t5b/` that these programs read; `tests/db_search/t5b/make_t5b.py` writes
those databases from the letters of `tests/db_search/ids/ndup.fna` and `ids/pdup.faa`. Added in commit 937985708
(2026-11-14).
