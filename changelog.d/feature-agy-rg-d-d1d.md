- BLAST database alias files: read `MEMB_BIT` across `blastdbcmd`, `blastdbcheck`, and every search
  program (`blastp`, `blastn`, `blastx`, `tblastn`, `tblastx`, `psiblast`, `rpsblast`, `rpstblastn`),
  matching NCBI BLAST+ 2.17.0: parses via `NStr::StringToUInt` between `SEQIDLIST` and `TAXIDLIST`,
  applies per volume in `x_ComputeFilters` with NCBI's `MEMB_BIT error: conflicting bit found.`
  (including `seqdbvol.hpp:881` in `blastdbcheck -verbosity 3/4`), filters `Blast-def-line.memberships`
  (first entry) before id lists and `OID_MASK_TYPE`, suppresses the top iteration range when combined
  with `FIRST_OID`/`LAST_OID`, and refuses values above 32 where a definition line carries membership
  bits (C++ 32-bit shift undefined behaviour). Also calls `build_oid_list` in `SeqDb::read_gi_mask_block`
  so `MASKLIST` reads report deferred alias OID-list errors as NCBI's `x_GetSeqGI` does.
