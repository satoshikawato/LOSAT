- blastx and rpstblastn read a query from standard input (`-query -`, or a run without `-query`)
  as NCBI BLAST+ reads it, from a pipe or a redirected file: a FASTA file whose line ends mix CR
  and CRLF and whose last line has no line end now keeps the last letters that NCBI keeps (107
  letters where LOSAT read 85). The same file named by its path reads as before, as NCBI reads a
  named file. blastn, blastp, tblastn, tblastx, psiblast and rpsblast already read standard input
  this way; all of them now share one rule for it.
- rpsblast and rpstblastn with an empty pipe as the query (`-query -`) and outfmt 0 (rpstblastn
  also outfmt 1 to 4) end as NCBI BLAST+ 2.17.0 does: the report stops after `Matrix: `, followed
  by "BLAST failed to write output" and exit status 6. Before, LOSAT wrote the matrix and gap
  costs and exited with status 0.
- blastx reads `-subject -` from standard input as NCBI BLAST+ reads it, from a pipe or a
  redirected file. With `-query -` (or no `-query`) as well, the subjects take all of standard
  input and no query is read, as in NCBI BLAST+: outfmt 7 writes `# BLAST processed 0 queries`.
  An empty standard input ends with "LOSAT engine error: Empty CBlastQueryVector" and exit status
  3. Before, LOSAT refused `-subject -`.
