- BLASTX with SEG (the default) no longer lets a query batch's setup wait behind other batches' searches
  on 2 or more threads. The setup's parallel steps (SEG on the batch's frames, their statistical blocks,
  the chunks of a long query) are now run by the batch's own thread with the pool's idle threads helping,
  and the 60,000-letter limit on parallel SEG is gone. A 200 kb query set against a small protein
  database with `-use_sw_tback` and 40,000-letter batches took 19.81 s instead of 21.62 s at 16 threads
  (median of 12 rounds; the slowest round 20.22 s instead of 23.33 s), and 6 queries of 10 kb against 595
  proteins 0.14 s instead of 0.15 s. Outputs do not change.
