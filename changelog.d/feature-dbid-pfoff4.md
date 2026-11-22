- A `-db` list of version 5 databases made with and without `makeblastdb -parse_seqids`, where `PRE_FETCH_SEQS_LIMIT`
  turns a query batch's prefetch off and a query's report reads a parsed Seq-id that a database of the list cannot
  look up: BLAST XML (`-outfmt 5`) now gives the results of NCBI BLAST+ instead of LOSAT's refusal, in BLASTP, BLASTN
  (every task), BLASTX, TBLASTN and TBLASTX. The hit of such a sequence shows its Seq-id (`ref|NP_100001.1|`), the
  name and accession `Unknown` and the length 0; each of its alignments keeps its number, scores and E-value, with the
  coordinates 0, the number of identities and no sequence letters. The report goes on with the next hit and the run
  ends normally (exit status 0), as NCBI BLAST+'s does. `-outfmt 17` keeps the refusal, and LOSAT's library still
  refuses such a sequence in a call that writes `-outfmt 5` with `-outfmt 0` to `4`.
- In the same lists, where the run has already loaded another sequence that holds the same Seq-id
  (`gi|5001|ref|NP_100001.1|` in an earlier query's report, for `ref|NP_100001.1|`), NCBI BLAST+ writes that
  sequence's names and letters for the hit instead of `Unknown`. LOSAT now refuses such a run, with the reason that
  another database sequence holds the Seq-id, instead of writing `Unknown` (in `-outfmt 0` to `7`, `10`, `13` to
  `16`, `18` and `20`; `-outfmt 10` refused it before for a crash that NCBI BLAST+ does not have there). Where that
  sequence is loaded only later in the report, NCBI BLAST+ and LOSAT write `Unknown` as before.
- LOSAT's library also refuses a sequence whose Seq-id a database of such a list cannot look up in a call whose
  formats NCBI BLAST+ would stop at different places after it, for example `-outfmt 5` (which goes on to the end) with `-outfmt "6 qseqid slen"` (which stops there with
  `BLAST Database error`).
