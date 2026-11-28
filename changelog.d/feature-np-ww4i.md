- LOSAT Web: makeblastdb can write a database above the run's folder or in another folder of
  yours, as the command line does (`-out ../db`, `-out /tmp/x/db`, also `-gi_mask_name ../m`); its
  files come back under those names and can be searched and downloaded like the run's other files.
  Before, such a run failed with "You do not have write permissions on '..'" or "Failed to create
  directory". A system folder that the command line cannot write (`-out /proc/db`) now gives the
  command line's message.
