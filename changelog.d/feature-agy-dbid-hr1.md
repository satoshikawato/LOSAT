- `-db` searches where another database sequence holds a reported subject's Seq-id, or matches it weakly (an accession
  with a version for one without, a name with or without an accession), but NCBI BLAST+ never gives that sequence for
  the subject, now give NCBI BLAST+'s results instead of LOSAT's refusal. LOSAT now follows each place where NCBI BLAST+
  finds a sequence for the Seq-id at each load of a subject in one walk per call over each batch's full prefetch set
  (the scope's cache including unresolved Seq-ids and its index of loaded Bioseqs by the definition line each was loaded
  by; the database data loader's 1,000-handle LRU `m_Ids` registry, without re-registering an already-loaded blob; the
  database lookup), in the order of the loads, in BLASTP, psiblast, BLASTN (every task), BLASTX, TBLASTN and TBLASTX.
  For example: beside a database sequence whose first definition line is `ref|NP_100001.1|` and whose second is `gi|7109`
  (its Seq-id in the report), a subject `ref|NP_100001|` or `emb|NP_100001.1|` of another database of the list is written
  as NCBI BLAST+ writes it, because NCBI BLAST+ loads the first sequence by its second line; a name (`sp||NAMEX`) in a
  version 5 database that no loaded sequence holds is written instead of refused; a handle evicted past `m_Ids`'s
  1,000-entry LRU resolves again by lookup; and psiblast's PSSM keeps an unresolved Seq-id in the scope when its lookup
  precedes a later registration in the batch's prefetch.
- `-outfmt 0` to `4` (with or without `-html`), where two database sequences share a description line's Seq-id
  (`gi|5301|ref|NP_100001.1|` and `gi|5302|ref|NP_100001.1|`): NCBI BLAST+'s description table looks each title up
  again by that Seq-id and can show the other sequence's title for a hit. LOSAT wrote the hit's own title there; it now
  refuses such a run, naming the other sequence, instead of writing a title that differs from NCBI BLAST+'s.
