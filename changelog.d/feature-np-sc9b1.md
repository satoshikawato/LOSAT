- BLASTX `-word_size 2` is 1.2 times as fast (1.09 times with PAM30) and `-word_size 4` 1.04 to 1.09 times: they
  now extend word hits the way the default word size 3 already did. 50 cDNAs x yeast proteins at 1 thread take 11.0 s
  instead of 13.2 s with word size 2 and 4.7 s instead of 5.0 s with word size 4; 800 cDNAs at 8 threads 30.0 s instead
  of 35.8 s (word size 2) and 15.8 s instead of 16.6 s (word size 4). Reports do not change. The added code makes the
  peak memory of a default BLASTX search up to about 0.25 MB larger (1 to 2 % for a small search).
