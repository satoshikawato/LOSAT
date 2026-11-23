- A `-db` list of version 5 databases made with and without `makeblastdb -parse_seqids`, where `PRE_FETCH_SEQS_LIMIT`
  turns a query batch's prefetch off and a query's report reads a parsed Seq-id that a database of the list cannot
  look up: BLAST XML (`-outfmt 5`) now gives the results of NCBI BLAST+ instead of LOSAT's refusal, in BLASTP, BLASTN
  (every task), BLASTX, TBLASTN and TBLASTX. The hit of such a sequence shows its Seq-id (`ref|NP_100001.1|`), the
  name and accession `Unknown` and the length 0; each of its alignments keeps its number, scores and E-value, with the
  coordinates 0, the number of identities and no sequence letters. The report goes on with the next hit and the run
  ends normally (exit status 0), as NCBI BLAST+'s does. `-outfmt 17` keeps the refusal, and LOSAT's library still
  refuses such a sequence in a call that writes `-outfmt 5` with `-outfmt 0` to `4`.
- In the same lists, where the run has already loaded another sequence that NCBI BLAST+ gives for the same Seq-id
  (`gi|5001|ref|NP_100001.1|` in an earlier query's report, for `ref|NP_100001.1|`; in the same query batch, or
  anywhere in the run with `-outfmt 5`, also `emb|NP_100001.1|` for it, or `ref|NP_100001.1|` for an unversioned
  `ref|NP_100001|`), NCBI BLAST+ writes that sequence's names and letters for the hit instead of `Unknown`. LOSAT now
  refuses such a run, with the reason that another database sequence holds the Seq-id, instead of writing `Unknown`
  (in `-outfmt 0` to `7`, `10`, `13` to `16`, `18` and `20`; `-outfmt 10` refused it before for a crash that NCBI
  BLAST+ does not have there). Where that sequence is loaded only later in the report, or only in an earlier query
  batch where NCBI BLAST+ no longer gives it (another type or an unversioned Seq-id, outside `-outfmt 5`), NCBI
  BLAST+ and LOSAT write `Unknown` as before. Some runs where that sequence comes first in the report, whose output
  NCBI BLAST+ writes the same either way and LOSAT wrote correctly before, are now refused too. A database sequence
  with several definition lines counts here by the line that NCBI BLAST+ loads it by (the line of its first `gi`, else
  its first line) for another type or an unversioned Seq-id, and by any line for the same Seq-id, as in NCBI BLAST+:
  beside `gi|7051|ref|XP_900001.1|` with a second line `gi|7052|ref|NP_100001.1|`, LOSAT writes `Unknown` for
  `emb|NP_100001.1|` and `ref|NP_100001|` again, as NCBI BLAST+ does, and refuses `ref|NP_100001.1|`. With `-html` and
  `-outfmt 0`, where NCBI BLAST+ ends the run at the first description line of such a hit (`BLAST Database error`,
  exit status 2) before a held sequence, LOSAT now writes that output instead of refusing the run.
- LOSAT's library also refuses a sequence whose Seq-id a database of such a list cannot look up in a call whose
  formats NCBI BLAST+ would end differently after it, with another exit status or another standard error, for
  example `-outfmt 5` (which goes on to the end) with `-outfmt "6 qseqid slen"` (which stops there with
  `BLAST Database error`), `-outfmt 8` with `-outfmt 0`, or `-outfmt 16` (which writes the error on standard error)
  with `-outfmt 6`.
