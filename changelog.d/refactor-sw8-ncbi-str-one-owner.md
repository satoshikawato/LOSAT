- Internal: `ncbi_str` (the port of NCBI's `NStr`) now owns the conversions that the database
  reader, blastdbcmd, the Seq-id parser, the taxid lists and makeblastdb's `-taxid_map` kept in
  copies or wrappers of their own: `StringToInt` with `fAllowTrailingSymbols` (now ported there),
  `StringToUInt`, `StringToInt` of the mask algorithm ids, `StringToNonNegativeInt` and blastdbcmd's
  `PrintableString`; the copies are removed. Output is unchanged.
