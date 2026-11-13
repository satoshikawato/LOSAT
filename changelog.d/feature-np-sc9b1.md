- BLASTX `-word_size 2` is 1.2 times as fast and `-word_size 4` 1.05 to 1.08 times: they now extend word hits the way
  the default word size 3 already did. 50 cDNAs x yeast proteins at 1 thread take 10.9 s instead of 13.1 s with word
  size 2 and 4.8 s instead of 5.1 s with word size 4; 800 cDNAs at 8 threads 30.0 s instead of 35.9 s (word size 2) and
  15.7 s instead of 16.6 s (word size 4). Memory and reports do not change.
