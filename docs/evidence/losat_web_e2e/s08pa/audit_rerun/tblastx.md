# TBLASTX rerun (LOSAT-head vs NCBI 2.17.0)

Class legend: SAME, LOSAT-REJECTS, ACCEPTED, DIFF, TIMEOUT. Work dir: /home/kawato/.cache/losat-web-gui-target/s08pa/audit-rerun/tblastx (outputs in r/)

## TX-1
LOSAT-REJECTS | ne=0 le=1 | id=cc9946da | -query q1.fna -subject s1.fna  -- -window_size 2147483640 -outfmt 6
   NCBI err: 
   LOSAT err: Error: a -window_size of 2147483640 with a query of 1201 letters (their sum is over 2^30, where NCBI BLAST+ never ends or wraps a 32-bit integer) is not supported by LOSAT's TBLASTX|
   stdout bytes ncbi=0 losat=0; first diff: 
LOSAT-REJECTS | ne=0 le=1 | id=61b872f0 | -query q2.fna -subject s3.fna  -- -window_size 2147483647
   NCBI err: 
   LOSAT err: Error: a -window_size of 2147483647 with a query of 2203 letters (their sum is over 2^30, where NCBI BLAST+ never ends or wraps a 32-bit integer) is not supported by LOSAT's TBLASTX|
   stdout bytes ncbi=919 losat=380; first diff: 14,57d13|< |
LOSAT-REJECTS | ne=0 le=1 | id=d3829c0e | -query q2.fna -subject s3.fna  -- -window_size 0x7fffffff
LOSAT-REJECTS | ne=0 le=1 | id=ee4223bc | -query q2.fna -subject s3.fna  -- -window_size 2147483600
LOSAT-REJECTS | ne=0 le=1 | id=3a7491dd | -query q2.fna -subject s3.fna  -- -window_size 2147483620
LOSAT-REJECTS | ne=0 le=1 | id=cb20215b | -query q2.fna -subject s3.fna  -- -window_size 2147483630
LOSAT-REJECTS | ne=0 le=1 | id=57a8a5ae | -query q2.fna -subject s3.fna  -- -window_size 2147483646
LOSAT-REJECTS | ne=0 le=1 | id=c62b7086 | -query q2.fna -subject s3.fna  -- -window_size 2147483000
TIMEOUT(ncbi=124,losat=1) | ne=124 le=1 | id=4fb8f792 | -query q2.fna -subject s3.fna  -- -window_size 2147481400
TIMEOUT(ncbi=124,losat=1) | ne=124 le=1 | id=c10a0b0f | -query q2.fna -subject s3.fna  -- -window_size 2147480000
(note: the last two TX-1 lines above are LOSAT-REJECTS with NCBI timed out (ne=124): window 2147481400 and 2147480000 on q2/s3; LOSAT message 'a -window_size of N with a query of 2203 letters (their sum is over 2^30...) is not supported by LOSAT's TBLASTX')
TX-1 verdict: LOSAT-REJECTS for all window values (incl. the NCBI-no-hits values).

## TX-2

## TX-3
LOSAT-REJECTS | ne=0 le=1 | id=a93cfe35 |     <sm.fna -- -query - -subject -
   NCBI err: 
   LOSAT err: Error: an empty query from a stream without a position (such as a pipe) is not supported by LOSAT's TBLASTX|
   stdout bytes ncbi=617 losat=0; first diff: 1,25d0|< TBLASTX 2.17.0+|
LOSAT-REJECTS | ne=0 le=1 | id=6e93f5fc |     <qm.fna -- -subject -
   NCBI err: 
   LOSAT err: Error: an empty query from a stream without a position (such as a pipe) is not supported by LOSAT's TBLASTX|
   stdout bytes ncbi=613 losat=0; first diff: 1,25d0|< TBLASTX 2.17.0+|
LOSAT-REJECTS | ne=0 le=1 | id=cd479bb2 |     <qm.fna -- -query - -subject -
   NCBI err: 
   LOSAT err: Error: an empty query from a stream without a position (such as a pipe) is not supported by LOSAT's TBLASTX|
   stdout bytes ncbi=613 losat=0; first diff: 1,25d0|< TBLASTX 2.17.0+|
pipe: cat sm.fna | tblastx -subject - : ncbi exit=0 stdout=617 stderr=; losat exit=1 stdout=0 stderr=Error: an empty query from a stream without a position (such as a pipe) is not supported by LOSAT's TBLASTX
pipe -query - -subject - : losat exit=1 stderr=Error: an empty query from a stream without a position (such as a pipe) is not supported by LOSAT's TBLASTX
  ncbi exit=0 stdout=617

## TX-4
SAME | ne=1 le=1 | id=b4d5bcb8 | -query q2.fna -subject s3.fna  -- -seg 12 2.2 2.5x
   NCBI err: BLAST query/options error: Invalid input for filtering parameters|Please refer to the BLAST+ user manual.|
   LOSAT err: BLAST query/options error: Invalid input for filtering parameters|Please refer to the BLAST+ user manual.|
   stdout bytes ncbi=0 losat=0; first diff: 
SAME | ne=1 le=1 | id=1001e584 | -query q2.fna -subject s3.fna  -- -seg 12 2,2 2.5
   NCBI err: BLAST query/options error: Invalid input for filtering parameters|Please refer to the BLAST+ user manual.|
   LOSAT err: BLAST query/options error: Invalid input for filtering parameters|Please refer to the BLAST+ user manual.|
   stdout bytes ncbi=0 losat=0; first diff: 
SAME | ne=1 le=1 | id=db25aedc | -query q2.fna -subject s3.fna  -- -seg 12 1e 2.5
   NCBI err: BLAST query/options error: Invalid input for filtering parameters|Please refer to the BLAST+ user manual.|
   LOSAT err: BLAST query/options error: Invalid input for filtering parameters|Please refer to the BLAST+ user manual.|
   stdout bytes ncbi=0 losat=0; first diff: 
SAME | ne=1 le=1 | id=fe9f0301 | -query q2.fna -subject s3.fna  -- -seg 12 0x 2.5
   NCBI err: BLAST query/options error: Invalid input for filtering parameters|Please refer to the BLAST+ user manual.|
   LOSAT err: BLAST query/options error: Invalid input for filtering parameters|Please refer to the BLAST+ user manual.|
   stdout bytes ncbi=0 losat=0; first diff: 
SAME | ne=1 le=1 | id=ba9d1bf4 | -query q2.fna -subject s3.fna  -- -seg 12 2.5x 2.5
SAME | ne=1 le=1 | id=b4d5bcb8 | -query q2.fna -subject s3.fna  -- -seg 12 2.2 2.5x
SAME | ne=1 le=1 | id=1001e584 | -query q2.fna -subject s3.fna  -- -seg 12 2,2 2.5
SAME | ne=1 le=1 | id=db2985f1 | -query q2.fna -subject s3.fna  -- -seg 12 2.2 2,2
SAME | ne=1 le=1 | id=0ae7f6ba | -query q2.fna -subject s3.fna  -- -seg 12 + 2.5
SAME | ne=1 le=1 | id=299c04bf | -query q2.fna -subject s3.fna  -- -seg 12 2.2 +
SAME | ne=1 le=1 | id=ec0a3326 | -query q2.fna -subject s3.fna  -- -seg 12 - 2.5
SAME | ne=1 le=1 | id=ba594753 | -query q2.fna -subject s3.fna  -- -seg 12 2.2 -
SAME | ne=1 le=1 | id=c255e4ac | -query q2.fna -subject s3.fna  -- -seg 12 . 2.5
SAME | ne=1 le=1 | id=16cf8291 | -query q2.fna -subject s3.fna  -- -seg 12 2.2 .
SAME | ne=1 le=1 | id=e1c5e105 | -query q2.fna -subject s3.fna  -- -seg 12 .e1 2.5
SAME | ne=1 le=1 | id=75d4ebe1 | -query q2.fna -subject s3.fna  -- -seg 12 2.2 .e1
SAME | ne=1 le=1 | id=db25aedc | -query q2.fna -subject s3.fna  -- -seg 12 1e 2.5
SAME | ne=1 le=1 | id=1ff9a02d | -query q2.fna -subject s3.fna  -- -seg 12 2.2 1e
SAME | ne=1 le=1 | id=93ad0112 | -query q2.fna -subject s3.fna  -- -seg 12 1e+ 2.5
SAME | ne=1 le=1 | id=1bc56e4e | -query q2.fna -subject s3.fna  -- -seg 12 2.2 1e+
SAME | ne=1 le=1 | id=a38a768d | -query q2.fna -subject s3.fna  -- -seg 12 1e- 2.5
SAME | ne=1 le=1 | id=9b1d77ff | -query q2.fna -subject s3.fna  -- -seg 12 2.2 1e-
SAME | ne=1 le=1 | id=fe9f0301 | -query q2.fna -subject s3.fna  -- -seg 12 0x 2.5
SAME | ne=1 le=1 | id=051a92ad | -query q2.fna -subject s3.fna  -- -seg 12 2.2 0x
SAME | ne=1 le=1 | id=3efec703 | -query q2.fna -subject s3.fna  -- -seg 12 1_0 2.5
SAME | ne=1 le=1 | id=ea3cb5d1 | -query q2.fna -subject s3.fna  -- -seg 12 2.2 1_0
SAME | ne=1 le=1 | id=8c0e3d5b | -query q2.fna -subject s3.fna  -- -seg 12 2.5, 2.5
SAME | ne=1 le=1 | id=b8c4edd7 | -query q2.fna -subject s3.fna  -- -seg 12 2.2 2.5,
SAME | ne=1 le=1 | id=64973d93 | -query q2.fna -subject s3.fna  -- -seg 12 2.2x 2.5
SAME | ne=1 le=1 | id=340fc442 | -query q2.fna -subject s3.fna  -- -seg 12 2.2 2.2x
-- strtod-readable, LOSAT cannot:
LOSAT-REJECTS | ne=0 le=1 | id=c5b72862 | -query q2.fna -subject s3.fna  -- -seg 12 0x10 2.5
LOSAT-REJECTS | ne=0 le=1 | id=88ef7948 | -query q2.fna -subject s3.fna  -- -seg 12 2.2 0x10
LOSAT-REJECTS | ne=0 le=1 | id=d14e4ef3 | -query q2.fna -subject s3.fna  -- -seg 12 0X1p3 2.5
LOSAT-REJECTS | ne=0 le=1 | id=f2c80819 | -query q2.fna -subject s3.fna  -- -seg 12 2.2 0X1p3
LOSAT-REJECTS | ne=0 le=1 | id=5272167c | -query q2.fna -subject s3.fna  -- -seg 12 +inf 2.5
LOSAT-REJECTS | ne=0 le=1 | id=2fa46ff6 | -query q2.fna -subject s3.fna  -- -seg 12 2.2 +inf
LOSAT-REJECTS | ne=0 le=1 | id=cbf8d04f | -query q2.fna -subject s3.fna  -- -seg 12 -infinity 2.5
LOSAT-REJECTS | ne=0 le=1 | id=763fe6cf | -query q2.fna -subject s3.fna  -- -seg 12 2.2 -infinity
LOSAT-REJECTS | ne=0 le=1 | id=0544342c | -query q2.fna -subject s3.fna  -- -seg 12 -nan 2.5
LOSAT-REJECTS | ne=0 le=1 | id=a92742dd | -query q2.fna -subject s3.fna  -- -seg 12 2.2 -nan
LOSAT-REJECTS | ne=0 le=1 | id=0fc20146 | -query q2.fna -subject s3.fna  -- -seg 12 1e400 2.5
LOSAT-REJECTS | ne=0 le=1 | id=af4c5530 | -query q2.fna -subject s3.fna  -- -seg 12 2.2 1e400
LOSAT-REJECTS | ne=0 le=1 | id=83de36fb | -query q2.fna -subject s3.fna  -- -seg 12 1e99999 2.5
LOSAT-REJECTS | ne=0 le=1 | id=b5454f77 | -query q2.fna -subject s3.fna  -- -seg 12 2.2 1e99999

## TX-5
LOSAT-REJECTS | ne=1 le=1 | id=4c8ef664 | -query q2.fna -subject xs.fna  -- -threshold 0
LOSAT-REJECTS | ne=1 le=1 | id=67e8f8ef | -query q2.fna -subject xs.fna  -- -seg abc
LOSAT-REJECTS | ne=1 le=1 | id=e9587713 | -query q2.fna -subject xs.fna  -- -evalue 0
LOSAT-REJECTS | ne=1 le=1 | id=15b8eef8 | -query q2.fna -subject xs.fna  -- -word_size 5
LOSAT-REJECTS | ne=1 le=1 | id=43713dc9 | -query q2.fna -subject nohdr2.fna  -- -threshold 0
LOSAT-REJECTS | ne=1 le=1 | id=500ec42c | -query q2.fna -subject nohdr2.fna  -- -seg abc
LOSAT-REJECTS | ne=1 le=1 | id=ed5b9725 | -query q2.fna -subject nohdr2.fna  -- -evalue 0
LOSAT-REJECTS | ne=1 le=1 | id=bf7c1846 | -query q2.fna -subject nohdr2.fna  -- -word_size 5
LOSAT-REJECTS | ne=1 le=1 | id=db0404e2 | -query q2.fna -subject s3.fna  -- -outfmt 5 -threshold 0
LOSAT-REJECTS | ne=1 le=1 | id=8e6aed37 | -query q2.fna -subject s3.fna  -- -outfmt 5 -seg abc
LOSAT-REJECTS | ne=0 le=1 | id=e2e0247e | -query empty.fna -subject s3.fna  -- -outfmt 5
   NCBI err: Warning: [tblastx] Query is Empty!|
   LOSAT err: Error: output format 5 is not supported by LOSAT's TBLASTX|
   stdout bytes ncbi=0 losat=0; first diff: 
LOSAT-REJECTS | ne=3 le=1 | id=7e5b7b05 | -query q2.fna -subject empty.fna  -- -outfmt 5
   NCBI err: BLAST engine error: Empty CBlastQueryVector|
   LOSAT err: Error: output format 5 is not supported by LOSAT's TBLASTX|
   stdout bytes ncbi=0 losat=0; first diff: 

## TX-9 (CLI side; validate not testable here)
SAME | ne=0 le=0 | id=42dbf338 | -query empty.fna -subject s3.fna  -- -window_size 0
   NCBI err: Warning: [tblastx] Query is Empty!|
   LOSAT err: Warning: [tblastx] Query is Empty!|
   stdout bytes ncbi=0 losat=0; first diff: 
SAME | ne=0 le=0 | id=f45c3926 | -query empty.fna -subject s3.fna  -- -word_size 2
   NCBI err: Warning: [tblastx] Query is Empty!|
   LOSAT err: Warning: [tblastx] Query is Empty!|
   stdout bytes ncbi=0 losat=0; first diff: 
SAME | ne=0 le=0 | id=20be3ce5 | -query empty.fna -subject s3.fna  -- -word_size 4 -threshold 15
   NCBI err: Warning: [tblastx] Query is Empty!|
   LOSAT err: Warning: [tblastx] Query is Empty!|
   stdout bytes ncbi=0 losat=0; first diff: 

## -version / -dryrun
LOSAT-REJECTS | ne=0 le=2 | id=013f4618 | -query q2.fna -subject s3.fna  -- -out -version
   NCBI err: 
   LOSAT err: error: the NCBI BLAST+ option -version is not supported by LOSAT's TBLASTX
   stdout bytes ncbi=68 losat=0; first diff: 1,2d0|< tblastx: 2.17.0+|
LOSAT-REJECTS | ne=1 le=2 | id=04e64bab | -query q2.fna -subject s3.fna  -- -out -dryrun
   NCBI err: USAGE|  tblastx [-h] [-help] [-import_search_strategy filename]|    [-export_search_strategy filename] [-db database_name]|    [-dbsize num_letters] [-gilist filename] [-seqidlist filename]|    [-negative_gilist filename] [-negative_seqidlist filename]|    [-taxids taxids] [-negative_taxids taxids] 
   LOSAT err: error: the NCBI C++ Toolkit option -dryrun is not supported by LOSAT's TBLASTX
   stdout bytes ncbi=0 losat=0; first diff: 
LOSAT-REJECTS | ne=0 le=2 | id=a4323d23 | -query q2.fna -subject s3.fna  -- -version
   NCBI err: 
   LOSAT err: error: the NCBI BLAST+ option -version is not supported by LOSAT's TBLASTX
   stdout bytes ncbi=68 losat=0; first diff: 1,2d0|< tblastx: 2.17.0+|
LOSAT-REJECTS | ne=0 le=2 | id=f8a1c2d5 | -query q2.fna -subject s3.fna  -- -dryrun
   NCBI err: 
   LOSAT err: error: the NCBI C++ Toolkit option -dryrun is not supported by LOSAT's TBLASTX
   stdout bytes ncbi=0 losat=0; first diff: 

## TX-6
LOSAT-REJECTS | ne=0 le=1 | id=79334239 | -query q2.fna -subject s3.fna  -- -num_threads 2147483647
   NCBI err: Warning: [tblastx] Number of threads was reduced to 32 to match the number of available CPUs|Warning: [tblastx] 'num_threads' is currently ignored when 'subject' is specified.|
   LOSAT err: Error: requested 2147483647 threads exceeds Rayon maximum 65535, which is not supported by LOSAT|
   stdout bytes ncbi=83656 losat=0; first diff: 1,3549d0|< TBLASTX 2.17.0+|
LOSAT-REJECTS | ne=0 le=1 | id=23238482 | -query q2.fna -subject s3.fna  -- -num_threads 65536
   NCBI err: Warning: [tblastx] Number of threads was reduced to 32 to match the number of available CPUs|Warning: [tblastx] 'num_threads' is currently ignored when 'subject' is specified.|
   LOSAT err: Error: requested 65536 threads exceeds Rayon maximum 65535, which is not supported by LOSAT|
   stdout bytes ncbi=83656 losat=0; first diff: 1,3549d0|< TBLASTX 2.17.0+|
LOSAT-REJECTS | ne=0 le=1 | id=24be48af | -query q2.fna -subject s3.fna  -- -num_threads 65535
   NCBI err: Warning: [tblastx] Number of threads was reduced to 32 to match the number of available CPUs|Warning: [tblastx] 'num_threads' is currently ignored when 'subject' is specified.|
   LOSAT err: Error: failed to build tblastx pool with 65535 threads; a thread count that the system cannot start is not supported by LOSAT||Caused by:|    0: Resource temporarily unavailable (os error 11)|    1: Resource temporarily unavailable (os error 11)|
   stdout bytes ncbi=83656 losat=0; first diff: 1,3549d0|< TBLASTX 2.17.0+|
LOSAT-REJECTS | ne=0 le=1 | id=53f35e81 | -query q2.fna -subject s3.fna  -- -num_threads 1000
   NCBI err: Warning: [tblastx] Number of threads was reduced to 32 to match the number of available CPUs|Warning: [tblastx] 'num_threads' is currently ignored when 'subject' is specified.|
   LOSAT err: Error: failed to build tblastx pool with 1000 threads; a thread count that the system cannot start is not supported by LOSAT||Caused by:|    0: Resource temporarily unavailable (os error 11)|    1: Resource temporarily unavailable (os error 11)|
   stdout bytes ncbi=83656 losat=0; first diff: 1,3549d0|< TBLASTX 2.17.0+|
ACCEPTED(num_threads warning only) | ne=0 le=0 | id=a0e5e2b0 | -query q2.fna -subject s3.fna  -- -num_threads 2
ACCEPTED(num_threads warning only) | ne=0 le=0 | id=16cf72e6 | -query q2.fna -subject s3.fna  -- -num_threads 4
ACCEPTED(num_threads warning only) | ne=0 le=0 | id=b8dddb6a | -query q2.fna -subject s3.fna  -- -num_threads 200

## TX-2 (all run lines)
LOSAT-REJECTS (NCBI timed out) | ne=124 le=1 | id=33743223 | -query q1.fna -subject s1.fna  -- -window_size 1500000000
LOSAT-REJECTS (NCBI timed out) | ne=124 le=1 | id=8a65d2a9 | -query q2.fna -subject s3.fna  -- -window_size 1073741823
LOSAT-REJECTS (NCBI timed out) | ne=124 le=1 | id=6cea3b88 | -query q2.fna -subject s3.fna  -- -window_size 1073741824
LOSAT-REJECTS (NCBI timed out) | ne=124 le=1 | id=31c12ec0 | -query q2.fna -subject s3.fna  -- -window_size 1073741825
LOSAT-REJECTS (NCBI timed out) | ne=124 le=1 | id=0930a8de | -query q2.fna -subject s3.fna  -- -window_size 2000000000
LOSAT-REJECTS (NCBI timed out) | ne=124 le=1 | id=c10a0b0f | -query q2.fna -subject s3.fna  -- -window_size 2147480000
LOSAT-REJECTS (NCBI timed out) | ne=124 le=1 | id=b7f10a12 | -query q1.fna -subject s1.fna  -- -window_size 1073741824
LOSAT-REJECTS (NCBI timed out) | ne=124 le=1 | id=72590c10 | -query q1.fna -subject s1.fna  -- -window_size 1073741825
LOSAT-REJECTS (NCBI timed out) | ne=124 le=1 | id=d7fc5338 | -query q1.fna -subject s1.fna  -- -window_size 2000000000
LOSAT-REJECTS (NCBI timed out) | ne=124 le=1 | id=013ec4d8 | -query q1.fna -subject s1.fna  -- -window_size 2147480000
SAME | ne=0 le=0 | id=ad2946df | -query q1.fna -subject s1.fna  -- -window_size 536870912
SAME | ne=0 le=0 | id=0149402e | -query q1.fna -subject s1.fna  -- -window_size 1000000000

TX-2 verdict: LOSAT-REJECTS for every value (NCBI timed out at 120 s for all of them); 536870912 and 1000000000 on q1/s1 run and are SAME.

## Harness re-runs (cmp.sh copy: NCBI vs LOSAT-head, timeout 120; sources: combos.txt, combos2.txt, probes.txt built from the Checked-OK list + ncbi_opts.txt + 149 names from cmdline_flags.cpp)
Raw per-case lines: out/hA.txt (combos.txt on qm/sm, 80), out/hB.txt (combos2.txt on qbig/sbig, 45), out/hC.txt (581 option probes on q2/s3), out/extra.txt (24 whitespace/empty-value probes). Every script finished well within 40 min (about 15 min each).
- hA:
         20 ACCEPTED(num_threads warning only) 
         60 SAME 
- hB:
         45 SAME 
- hC:
          8 ACCEPTED(num_threads warning only) 
        246 ACCEPTED(parser exit2 vs 1) 
          2 DIFF 
        127 LOSAT-REJECTS 
        198 SAME 
- extra:
          8 ACCEPTED(parser exit2 vs 1) 
         16 SAME 

Combined: SAME 319, ACCEPTED(num_threads warning only) 28, ACCEPTED(parser exit 2 vs 1, 'error: ') 254, LOSAT-REJECTS 127, DIFF 2, TIMEOUT 0
LOSAT-REJECTS in hC are the expected ones: word_size 2/4, window_size 0, evalue/threshold hex, outfmt 1-5/8-12/15/16/18/20 and custom field lists, every unported option NCBI accepts (message 'not supported by LOSAT's TBLASTX'), -version, -dryrun.

## DIFF list
1. -query q2.fna -subject s3.fna -help : NCBI exit 0 prints 'USAGE ...' (13801 bytes); LOSAT-head exit 0 prints its own help ('Pairwise 6-frame translated nucleotide alignment (tblastx) / Usage: LOSAT-head tblastx [OPTIONS]', 1323 bytes). Identical in LOSAT-base (unchanged, not in the Round-1 report/decisions).
2. -query q2.fna -subject s3.fna -help 1 : NCBI exit 0 (USAGE); LOSAT exit 2 'error: unknown option or argument '1'; use -help for CLI v2 syntax'. Identical in LOSAT-base (unchanged).

## Summary table
| Finding | Verdict | Evidence |
|---|---|---|
| TX-1 window >= 2^31-qlen | LOSAT-REJECTS | all 9 values (2147483000..2147483647, 0x7fffffff, q1 -outfmt 6 and q2/s3 outfmt 0): exit 1 "a -window_size of N with a query of L letters (their sum is over 2^30, where NCBI BLAST+ never ends or wraps a 32-bit integer) is not supported by LOSAT's TBLASTX"; NCBI printed no hits / report |
| TX-2 window in (2^30-qlen, 2^31-qlen) | LOSAT-REJECTS | 1073741823..2147480000, 1500000000 on q1/q2: same message, exit 1 in ms; NCBI timed out at 120 s in every case; 536870912 and 1000000000 (q1/s1) still run and are SAME |
| TX-3 -query - -subject - | LOSAT-REJECTS | file stdin and pipe, also -subject - alone: exit 1 "an empty query from a stream without a position (such as a pipe) is not supported by LOSAT's TBLASTX"; NCBI exit 0 with the empty-query report (accepted per ROUND1: explicit rejection like BLASTP/TBLASTN) |
| TX-4 -seg malformed locut/hicut | SAME / LOSAT-REJECTS | 2.5x, 2,2, 1e, 0x, +, -, ., .e1, 1e+, 1e-, 1_0, 2.5, , 2.2x (both positions, 26 runs): byte-identical "Invalid input for filtering parameters", exit 1; 0x10, 0X1p3, +inf, -infinity, -nan, 1e400, 1e99999 (14 runs): NCBI runs, LOSAT-REJECTS (explicit) |
| TX-5 LOSAT rejection before NCBI option error | ACCEPTED | 10 runs both exit 1 (LOSAT explicit rejection: subject residue X / missing header / outfmt 5); -s empty.fna -outfmt 5 exit 3 vs 1 and -q empty.fna -outfmt 5 exit 0 vs 1 are listed rejections (ROUND1: accepted) |
| TX-6 num_threads phrase | LOSAT-REJECTS | 2147483647 and 65536: "...exceeds Rayon maximum 65535, which is not supported by LOSAT"; 65535 and 1000 (ulimit -v 12000000): "failed to build tblastx pool ...; a thread count that the system cannot start is not supported by LOSAT"; 2/4/200: stdout SAME, only NCBI's thread warnings differ (approved) |
| TX-7 wrong NCBI line range | not command-testable | comment fix; ROUND1: fixed |
| TX-8 web ABI v1 -evalue | not run | web ABI v1 frozen, ROUND1: accepted (no web binary in scope) |
| TX-9 validate vs CLI empty query | SAME (CLI side) | -q empty.fna -s s3.fna with -window_size 0 / -word_size 2 / -word_size 4 -threshold 15: byte-identical "Query is Empty!" exit 0; adapter validate not built/run; ROUND1: accepted |
| TX-10 AUTHORITY §K/§M | fixed (doc) | head-src docs/evidence/losat_web_e2e/AUTHORITY.md line 71 (§M) and line 115 (§K) now both state the 2^30 sum rejection (D8/D11) |
| -out -version | LOSAT-REJECTS | exit 2 "error: the NCBI BLAST+ option -version is not supported by LOSAT's TBLASTX"; NCBI prints "tblastx: 2.17.0+" exit 0 |
| -out -dryrun | LOSAT-REJECTS | exit 2 "error: the NCBI C++ Toolkit option -dryrun is not supported by LOSAT's TBLASTX"; NCBI prints USAGE exit 1 |
| Checked-OK harness (cmp.sh, combos, option probes) | SAME / ACCEPTED except 2 DIFF | 125 combos: 105 SAME + 20 num_threads-warning-only; 581+24 probes: 2 DIFF (-help, -help 1) |

Notes: stdin for non-stdin cases is /dev/null. refcheck.py (NCBI line-reference checker) needs the old session.diff and was not re-run. probes.txt reader drops empty/whitespace values, so those were re-run separately in out/extra.txt (all SAME or accepted parser).
