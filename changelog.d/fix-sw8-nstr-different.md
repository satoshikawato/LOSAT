- The real arguments of every BLAST program and of rpsbproc (`-evalue`, `-xdrop_ungap`, `-xdrop_gap`,
  `-xdrop_gap_final`, `-perc_identity`, `-qcov_hsp_perc`, `-threshold`, `-gap_trigger`, `-inclusion_ethresh`,
  rpsbproc's `-e`) take every value that NCBI BLAST+ 2.17.0 takes, with NCBI's number: hexadecimal values (`0x10`,
  `0x1p-3`), a value that ends in `e` or `E` (`1e`, `2.5E`, `12345678901234567890e`, computed in NCBI's x87
  extended precision) and NaN with its sign and payload (`-nan(1)`). They were refused before. A value that NCBI
  refuses (`1x`, `1e+`, `2,5`, `inf`) gets "expected a number". `-threshold` checks its value as NCBI does, so
  `-threshold 1e` stays an error there, as in NCBI.
- `-seg` takes a hexadecimal locut or hicut as NCBI does (`-seg "12 0x10 2.5"` searches with a locut of 16). When
  one value is an infinity or NaN (which LOSAT refuses) and the other is not a number, NCBI's error now comes first,
  as in NCBI.
- segmasker takes a `-locut` or `-hicut` with more than 18 digits or 15 decimal places before a final `e`, with
  NCBI's value; it was refused.
- A `-in_pssm` file (psiblast, tblastn) or a rpsbproc archive with a REAL of more than 18 significant digits is read
  with NCBI's value; it was refused.
- `-in_pssm` (psiblast, tblastn) counts a vertical tab as white space, as NCBI does: an ASN.1 text PSSM with vertical
  tabs in a comment is searched (it got "Unsupported format for PSSM"), and an XML file that begins with a vertical
  tab is refused as an XML PSSM, as one without it is.
