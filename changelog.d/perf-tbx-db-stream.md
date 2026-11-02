- TBLASTX's `-db` search reads each database sequence where a query batch searches it and where the batch reports it, as
  TBLASTN's does, instead of a copy of the whole database made before the first batch (on core_nt 269 GiB of memory and more
  than 5 minutes before the first batch; 11-34 MiB on 16S rRNA, zebrafish cDNA and fly genome databases). A short sequence is
  translated from the bytes of a block of the sequence file that the threads share, without a copy of its own. Outputs and
  errors do not change (a disk error in the data of a long sequence comes where a batch reads it, as TBLASTN's does).
