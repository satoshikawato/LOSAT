- New program `blastdb_aliastool`, every mode of NCBI BLAST+ 2.17.0's: alias files of several databases (`-dblist`,
  `-dblist_file`), of a database with a GI, sequence id or taxonomy id list or its exclude-model mask (`-db` with
  `-gilist`, `-seqidlist`, `-taxidlist`, `-oid_masks`), of the volumes of one (`-num_volumes`) and of VDB runs
  (`-vdblist`); binary GI lists (`-gi_file_in`) and binary sequence id lists (`-seqid_file_in`, `-seqid_file_info`).
  It writes NCBI's files byte for byte (apart from their date) with NCBI's messages and exit statuses. An alias
  file with a sequence id or taxonomy id list over a version 5 database is refused until LOSAT's reader reads it.
