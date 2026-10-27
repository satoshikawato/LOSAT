- The database reader decodes every member of a `-parse_seqids` header record (every Seq-id
  alternative, taxids, memberships, links, other-info) in the new shared module `LOSAT/src/seqid/`,
  which also holds NCBI's id printings and choices (FASTA form, labels, ranks, the Seq-align id, the
  joins of the all-definition-line fields). The search programs still refuse such a database with
  the same message, now right after they open it. No output change.
