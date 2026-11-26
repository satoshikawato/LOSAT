- A database sequence whose Seq-id NCBI BLAST+ cannot look up by the id itself (a text tag of `gnl|PIG|` or
  `gnl|BL_ORD_ID|`, such as the `gnl|BL_ORD_ID|0` that `makeblastdb` stores for a text id, or a `gnl|ti|` tag that is not
  a number) no longer makes LOSAT refuse the run in every format. BLASTP, psiblast, BLASTN, BLASTX, TBLASTN and TBLASTX now
  write NCBI BLAST+'s results, with the batch's preload on or off: `Warning: [<program>] Error pre-fetching sequence
  data` where the batch's preload fails there, `Unknown`, length 0 and `N/A` for that sequence in BLAST XML, the archive,
  XML2, JSON and the tables that do not print its accession, and NCBI BLAST+'s `Error:` and exit status 255 where a format
  stops at it (the pairwise and query-anchored reports, outfmt 13 and 14, and the table field lists that NCBI BLAST+
  stops at). psiblast's PSSM reads none of its letters (`Failed to retrieve sequence`), as NCBI BLAST+'s does. LOSAT still
  refuses a table row of such a sequence that prints `sacc` or `saccver` (the default fields of outfmt 6, 7, 10 and 20),
  where NCBI BLAST+ crashes, and BLASTN's outfmt 17; the batch's other sequences are refused only as where its preload is
  off (an id that another sequence's id can match weakly, where NCBI BLAST+ can write that sequence's names for it).
- With `-outfmt 15` or `16` and a taxonomy database whose reads write a message (for example `Offset error at end of
  taxdb file`), LOSAT now writes that message where NCBI BLAST+ does, right after the document's opening lines and before
  the first query's report, so the output of `2>&1` matches NCBI BLAST+'s; it came after the first query's report before.
- With `-html` and outfmt 0 to 4, BLASTX now writes the messages of a batch's preload (`Error pre-fetching sequence data`,
  or a taxonomy database's message such as `Offset error at end of taxdb file`) after the report's reference and
  database lines, so the output of `2>&1` matches NCBI BLAST+'s; they came right after the program's name line before.
