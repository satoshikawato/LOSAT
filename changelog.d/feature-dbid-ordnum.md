- `-db` on a database that a list or an alias file restricts (`-taxids`, `-negative_taxids`, `-seqidlist`,
  `-negative_seqidlist`, an alias file's `FIRST_OID` and `LAST_OID`), in every search program: LOSAT's refusals now
  name a database sequence by its OID ("database sequence n", as `blastdbcmd -outfmt %o` prints it), where they named
  its rank among the searched sequences. The output of every run that LOSAT does not refuse is unchanged.
