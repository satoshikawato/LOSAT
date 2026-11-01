- TBLASTN's `-db` stream fixes its blocks of 256 KiB once, from the index, and the threads share them (a block is read once while
  a thread holds it, instead of once per thread from the sequence it asked for). The lengths of the database's sequences are held
  in 32 bits (`u32`, half the memory: about 0.5 GB on core_nt), with the longest and the sum, so a search no longer walks them for
  those. Outputs and errors do not change.
