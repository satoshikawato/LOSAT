- A `-db` list of version 5 databases made with and without `makeblastdb -parse_seqids`, where a query batch stops
  at a parsed Seq-id that a database of the list cannot look up and another reported sequence's load writes a
  message of the taxonomy database (a damaged `taxdb`): every search program now writes the messages of the
  sequences that NCBI BLAST+ loads before the stop (in the order of its Seq-id handles: gis, accessions, and
  `gnl|db|n` such as `gnl|BL_ORD_ID|n`, by their number, other ids last) and `LOSAT database error: DB contains no
  accession info.` (exit status 2), as NCBI BLAST+ does; LOSAT refused such a run. A load whose order NCBI takes from
  memory addresses (the same number, ids that are not packed, more than 200 ids in the batch) that writes a message
  stays refused. Where a batch has several such Seq-ids, the first that NCBI BLAST+ looks up ends its loads, as in
  NCBI BLAST+; where one before it may never be looked up, a message written by a load between the two, or by the
  sequence that NCBI BLAST+ would load for the earlier one instead, is refused (LOSAT does not know whether NCBI BLAST+
  looks the earlier one up), and so are two such Seq-ids whose lookups fail with different messages. Newly refused, where NCBI BLAST+ may never look such a Seq-id up and LOSAT cannot tell whether it goes
  on: a stopping Seq-id that another sequence of the database holds under another Seq-align id where that sequence
  may be loaded first (in the batch, or in a batch, query or psiblast round before: by its prefetch, or with
  `PRE_FETCH_SEQS_LIMIT` by the reports that NCBI BLAST+ writes, after their cuts to `-max_target_seqs` (in the
  pairwise report, `-outfmt 0`: the sequences of both its description lines, cut to `-num_descriptions`, and its
  alignments, cut to `-num_alignments`), and by psiblast's PSSM), and the id of a query that NCBI BLAST+ holds at that
  point (the batch's queries as read, one that `-query_loc` skips included, numbered over the whole run; with
  `-parse_deflines` the ids of the deflines; the PSSM's query of `-in_pssm`; with `-outfmt 5` also the queries
  before).
- After a BLAST XML (`-outfmt 5`) iteration that NCBI BLAST+'s stream leaves out at a query title outside US-ASCII
  (with `PRE_FETCH_SEQS_LIMIT`, whose sequences it then does not load), a later query batch whose prefetch stops at
  such a Seq-id now writes, in every search program (`blastp`, `blastx`, `blastn`, `tblastn`, `tblastx`), the
  taxonomy database's messages that those loads would have written where NCBI BLAST+ writes them: at the first load
  after the left-out iteration that makes the same lookup, in a query batch between it and the stop (its reports or
  its prefetch) or else in the stopping batch's prefetch before the stop (in NCBI BLAST+'s order of the loads), or not
  at all where only loads after the stop would make it, then `LOSAT database error: DB contains no accession info.`
  (exit status 2); LOSAT refused such a batch where any of its other loads looked a tax id up. The sequences of the
  left-out iteration no longer count as loaded: a stop at a Seq-id that one of them holds is now NCBI BLAST+'s
  error, as in NCBI BLAST+, where LOSAT refused it. A load whose order to the stop LOSAT does not know that would
  write such a message stays refused.
- A version 4 database whose string ISAM file (`.psi`, `.psd`, `.nsi`, `.nsd`) or PIG ISAM file (`.ppi`, `.ppd`) is a
  directory: a volume with such a file now has no such index, as in NCBI BLAST+, whose lookups find nothing there and go
  on to the next volume; LOSAT stopped with `BLAST Database error: Error: Could not open input file`. A reported
  sequence that the lookup then finds as another sequence of a later volume is refused, as in every other database,
  and the refusal of an empty string ISAM file now says that NCBI BLAST+ crashes where it reads it. Trace-id ISAM
  files (`.pti` and `.ptd`, `.nti` and `.ntd`) that are directories likewise are no such index (NCBI BLAST+ takes only
  regular files for one): such a volume is searched as NCBI BLAST+ searches it, as a volume without trace-id files;
  LOSAT refused it as one with an ISAM of trace ids.
