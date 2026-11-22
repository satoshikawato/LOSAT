- A `-db` list of version 5 databases made with and without `makeblastdb -parse_seqids`, where `PRE_FETCH_SEQS_LIMIT`
  turns a query batch's prefetch off and a query's report reads a parsed Seq-id that a database of the list cannot
  look up: BLAST XML (`-outfmt 5`) now gives the results of NCBI BLAST+ instead of LOSAT's refusal, in BLASTP, BLASTN
  (every task), BLASTX, TBLASTN and TBLASTX. The hit of such a sequence shows its Seq-id (`ref|NP_100001.1|`), the
  name and accession `Unknown` and the length 0; each of its alignments keeps its number, scores and E-value, with the
  coordinates 0, the number of identities and no sequence letters. The report goes on with the next hit and the run
  ends normally (exit status 0), as NCBI BLAST+'s does. `-outfmt 17` keeps the refusal, and LOSAT's library still
  refuses such a sequence in a call that writes `-outfmt 5` with `-outfmt 0` to `4`.
