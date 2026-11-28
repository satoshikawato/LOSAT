- Web: a list file of the Database filters (`-seqidlist`, `-negative_seqidlist`, `-gilist`, `-taxidlist` and the
  others) whose name starts or ends with a space was run, and shown in Run details, under the name without
  the space (`-seqidlist ids.txt`), so the shown command failed on the command line. The file is now run
  and shown under its own name, quoted as a query file is (`-seqidlist 'ids.txt '`).
- Web: in the tool tabs, a file whose name starts or ends with a space (a FASTA input of WindowMasker or Masking, the
  counts of `-ustat` and Convert counts, a file chosen for `-ids`, `-taxid_map` or `-mask_data`, the archive of Format
  Archive) was put in the run's folder under its name but named in the command without the space, so the run could not
  open it. The command now names it as it is (`-in 'a.fa '`), and the run reads it; a name typed in a box is still
  taken without the spaces around it.
- Web: Run details and the Results summary quote the name of a `-in_msa` or `-in_pssm` file as the command does
  (`-in_msa 'seed.aln '`).
- Web: after a help or version word (`-version`, `-h`) is typed in Additional options, the Run LOSAT button
  no longer says "(2 runs)" for a moment: it names no count until LOSAT has read the new words, then the
  count of the argv it will run, and the line under the button changes with it. The WindowMasker and Masking tabs'
  buttons say one run for a help or version word, as the press starts one.
