- LOSAT Web: the Database picker offers every database in the files or folder you choose, also those
  that an alias file lists (`DBLIST`): with `nab.nal` listing `na` and `nb`, `na` and `nb` can be
  searched on their own, as with `losat blastn -db na`. Only the numbered volumes of one database
  (`foo.00`, `foo.01` ... as `foo.nal` lists them) stay out of the list; other names, such as
  `refseq.2023`, are databases like any other. After a reload, the picker, the message by
  Run LOSAT and the "Reading the database" line ask for what you chose: "Choose the files again."
  for files, "Choose the folder again." for a folder (it always said folder).
  The database you had chosen before the reload is chosen again when you pick it again.
