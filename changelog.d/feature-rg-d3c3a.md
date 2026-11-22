- psiblast and tblastn `-in_pssm` read a PSSM in ASN.1 binary (BER), as NCBI BLAST+ 2.17.0 does:
  the same values, messages and errors as its ASN.1 text twin (`Unrecognized format for PSSM` for
  any error of the reader, nothing written and exit status 0 for a file cut inside the last
  end-of-contents octets of its object, `BLAST ran out of memory` for a string longer than any
  host can hold). LOSAT refused every binary PSSM. Still refused: a PSSM in XML, a string whose
  length passes the end of the file and whose allocation of 1 GiB or more succeeds or fails by
  the host's memory, and a query whose letters are not its stated length, on which NCBI crashes.
