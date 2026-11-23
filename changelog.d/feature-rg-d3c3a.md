- psiblast and tblastn `-in_pssm` read a PSSM in ASN.1 binary (BER), as NCBI BLAST+ 2.17.0 does:
  the same values, messages and errors as its ASN.1 text twin (`Unrecognized format for PSSM` for
  any error of the reader, nothing written and exit status 0 for a file cut inside the last
  end-of-contents octets of its object, `BLAST ran out of memory` for a string or byte block
  longer than any host can hold, after NCBI's `unended byte block` for a byte block). LOSAT
  refused every binary PSSM. Still refused: a PSSM in XML, a string or byte block whose length
  passes the end of the file and whose allocation of 1 GiB or more succeeds or fails by the
  host's memory, and a query whose letters are not its stated length (NCBI crashes or writes past
  its buffer for most such lengths; for one less than the letters it answers, which LOSAT does
  not yet do).
- psiblast and tblastn `-in_pssm`: a PSSM whose gapped `lambda` is infinite or whose gapped
  `lambda`, `kappa` or `h` is not a number (psiblast also an infinite `lambdaUngapped`), in ASN.1
  text or binary, is refused; LOSAT wrote E-values NaN or other hits than NCBI BLAST+ 2.17.0.
- rpsbproc, makeblastdb `-mask_data` and psiblast/tblastn `-in_pssm` read ASN.1 with one reader for
  text and binary: rpsbproc of a 336 MB archive runs 6.8 % fewer instructions, `-mask_data` of a
  text mask file 4.4 % fewer, of a binary genome mask 0.8 % fewer (read as a stream instead of the
  whole file); outputs and peak memory unchanged.
