- TBLASTN's `-db` stream fixes its blocks of 256 KiB once, from the index, and the threads share them (a block is read once while
  a thread holds it, instead of once per thread from the sequence it asked for). The lengths of the database's sequences are held
  in 32 bits (`u32`, half the memory: about 0.5 GB on core_nt), with the longest and the sum, so a search no longer walks them for
  those. Outputs and errors do not change.
- With several threads, TBLASTN's `-db` blocks are no longer cut at the end of a thread's chunk of subjects
  (`LOSAT_X_TBNPLANCUT`): the threads of the neighbouring chunks share the block of the grid instead. A thread still lets
  its block go with its chunk.
