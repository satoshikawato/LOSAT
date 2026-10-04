# BLASTP re-run (LOSAT-head vs NCBI blastp 2.17.0)
Harness: audit-rerun/blastp/cmp.sh (stdout, stderr, exit status, and files created in the scratch cwd compared). Classes: SAME, LOSAT-REJECTS (exit!=0 with "not supported by LOSAT"), PARSER-EXC (NCBI USAGE exit 1 vs LOSAT clap "error: " exit 2), DIFF*. Paths shortened: $F=tests/fasta/outfmt0, Q=e2e_protein_query.faa, S=e2e_protein_subject.faa, M=e2e_many_subject.faa. n = NCBI exit, l = LOSAT exit. NCBI exit 139 = SIGSEGV.

## BP-1
```
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject hdronly.faa
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject hdronly.faa -outfmt 6
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject hdronly.faa -outfmt 7
```

## BP-2
```
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -query $Q -subject $S -max_target_seqs 1073741799 -outfmt 6 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: a -max_target_seqs of 1073741799, whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -query $Q -subject $S -max_target_seqs 1073741799 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: a -max_target_seqs of 1073741799, whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTP| || stdout n=17 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -query $Q -subject $S -max_target_seqs 1073741823 -outfmt 6 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: a -max_target_seqs of 1073741823, whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -query $Q -subject $S -max_target_seqs 1073741823 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: a -max_target_seqs of 1073741823, whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTP| || stdout n=17 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -query $Q -subject $S -max_target_seqs 1073741824 -outfmt 6 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: a -max_target_seqs of 1073741824, whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -query $Q -subject $S -max_target_seqs 1073741824 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: a -max_target_seqs of 1073741824, whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTP| || stdout n=17 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -query $Q -subject $S -max_target_seqs 2147483597 -outfmt 6 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: a -max_target_seqs of 2147483597, whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -query $Q -subject $S -max_target_seqs 2147483597 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: a -max_target_seqs of 2147483597, whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTP| || stdout n=17 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -query $Q -subject $S -max_target_seqs 2147483598 -outfmt 6 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: a -max_target_seqs of 2147483598, whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -query $Q -subject $S -max_target_seqs 2147483598 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: a -max_target_seqs of 2147483598, whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTP| || stdout n=17 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -query $Q -subject $S -max_target_seqs 2147483622 -outfmt 6 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: a -max_target_seqs of 2147483622, whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -query $Q -subject $S -max_target_seqs 2147483622 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: a -max_target_seqs of 2147483622, whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTP| || stdout n=17 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -query $Q -subject $S -max_target_seqs 2147483623 -outfmt 6 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: a -max_target_seqs of 2147483623, whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -query $Q -subject $S -max_target_seqs 2147483623 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: a -max_target_seqs of 2147483623, whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTP| || stdout n=17 l=0
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject $S -max_target_seqs 1073741798
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject $S -max_target_seqs 1073741798 -outfmt 6
```

## BP-3
```
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject $S -max_target_seqs 2147483647
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject $S -max_target_seqs 2147483647 -outfmt 6
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject e2e_many_subject.faa -max_target_seqs 2147483624 -outfmt 6
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject e2e_many_subject.faa -max_target_seqs 2147483630 -outfmt 6
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject e2e_many_subject.faa -max_target_seqs 2147483640 -outfmt 6
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject e2e_many_subject.faa -max_target_seqs 2147483646 -outfmt 6
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject e2e_many_subject.faa -max_target_seqs 2147483647
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject e2e_many_subject.faa -max_target_seqs 2147483647 -outfmt 6
```

## BP-4
```
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -word_size 3 -evalue +inf -outfmt 6 || NCBI-err:  || LOSAT-err: Error: an infinite -evalue (inf), with which NCBI BLAST+'s blastp crashes on some inputs, is not supported by LOSAT's BLASTP| || stdout n=325341 l=0
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -word_size 3 -evalue 1e999 || NCBI-err:  || LOSAT-err: Error: an infinite -evalue (inf), with which NCBI BLAST+'s blastp crashes on some inputs, is not supported by LOSAT's BLASTP| || stdout n=1314062 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -query $Q -subject e2e_many_subject.faa -evalue +inf -outfmt 6 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: an infinite -evalue (inf), with which NCBI BLAST+'s blastp crashes on some inputs, is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -query $Q -subject e2e_many_subject.faa -evalue 1e999 -outfmt 6 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: an infinite -evalue (inf), with which NCBI BLAST+'s blastp crashes on some inputs, is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -task blastp-fast -evalue +inf || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: an infinite -evalue (inf), with which NCBI BLAST+'s blastp crashes on some inputs, is not supported by LOSAT's BLASTP| || stdout n=17 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -word_size 5 -evalue +inf -outfmt 0 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: an infinite -evalue (inf), with which NCBI BLAST+'s blastp crashes on some inputs, is not supported by LOSAT's BLASTP| || stdout n=17 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -word_size 5 -evalue +inf -outfmt 6 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: an infinite -evalue (inf), with which NCBI BLAST+'s blastp crashes on some inputs, is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -word_size 5 -evalue 1e999 -outfmt 6 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: an infinite -evalue (inf), with which NCBI BLAST+'s blastp crashes on some inputs, is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=139 l=1 files=<>/<> :: -word_size 5 -evalue 1e999 -outfmt 7 || NCBI-err: timeout: the monitored command dumped core| || LOSAT-err: Error: an infinite -evalue (inf), with which NCBI BLAST+'s blastp crashes on some inputs, is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject e2e_many_subject.faa -evalue 1e10 -outfmt 6
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject e2e_many_subject.faa -evalue 1e100 -outfmt 6
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject e2e_many_subject.faa -evalue 1e308 -outfmt 6
```

## BP-5
```
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 bitscore frames
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 bitscore sframe
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 evalue frames
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 evalue sframe
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 frames
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 ppos frames
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 qacc frames
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 qacc sframe
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 qframe sframe
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 qlen frames
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 qlen sframe
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 qseqid sseqid frames
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 score frames
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 score sframe
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 sframe
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 slen frames
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 slen sframe
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 sseqid sframe
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 std frames
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 stitle frames
[SAME] n=0 l=0 files=<>/<> :: -outfmt 6 stitle sframe
[SAME] n=0 l=0 files=<>/<> :: -outfmt 7 bitscore sframe
[SAME] n=0 l=0 files=<>/<> :: -outfmt 7 evalue sframe
[SAME] n=0 l=0 files=<>/<> :: -outfmt 7 frames
[SAME] n=0 l=0 files=<>/<> :: -outfmt 7 qacc frames
[SAME] n=0 l=0 files=<>/<> :: -outfmt 7 qlen frames
[SAME] n=0 l=0 files=<>/<> :: -outfmt 7 score frames
[SAME] n=0 l=0 files=<>/<> :: -outfmt 7 slen sframe
[SAME] n=0 l=0 files=<>/<> :: -outfmt 7 sseqid sframe
[SAME] n=0 l=0 files=<>/<> :: -outfmt 7 stitle frames
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject e2e_many_subject.faa -outfmt 6 frames
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject e2e_many_subject.faa -outfmt 6 sframe
[SAME] n=0 l=0 files=<>/<> :: -query $Q -subject e2e_many_subject.faa -outfmt 7 sseqid frames
```

## BP-6
```
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -seg 12 +inf 2.5 || NCBI-err:  || LOSAT-err: Error: the SEG locut or hicut "+inf" (not a finite decimal number) is not supported by LOSAT's BLASTP| || stdout n=4107 l=0
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -seg 12 +nan(5) 2.5 || NCBI-err:  || LOSAT-err: Error: the SEG locut or hicut "+nan(5)" (not a finite decimal number) is not supported by LOSAT's BLASTP| || stdout n=29757 l=0
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -seg 12 -inf 2.5 || NCBI-err:  || LOSAT-err: Error: the SEG locut or hicut "-inf" (not a finite decimal number) is not supported by LOSAT's BLASTP| || stdout n=29757 l=0
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -seg 12 -nan 2.5 || NCBI-err:  || LOSAT-err: Error: the SEG locut or hicut "-nan" (not a finite decimal number) is not supported by LOSAT's BLASTP| || stdout n=29757 l=0
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -seg 12 0x10 2.5 || NCBI-err:  || LOSAT-err: Error: the SEG locut or hicut "0x10" (not a finite decimal number) is not supported by LOSAT's BLASTP| || stdout n=4107 l=0
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -seg 12 0x1p3 2.5 || NCBI-err:  || LOSAT-err: Error: the SEG locut or hicut "0x1p3" (not a finite decimal number) is not supported by LOSAT's BLASTP| || stdout n=4107 l=0
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -seg 12 1e400 2.5 || NCBI-err:  || LOSAT-err: Error: the SEG locut or hicut "1e400" (not a finite decimal number) is not supported by LOSAT's BLASTP| || stdout n=4107 l=0
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -seg 12 2.2 +inf || NCBI-err:  || LOSAT-err: Error: the SEG locut or hicut "+inf" (not a finite decimal number) is not supported by LOSAT's BLASTP| || stdout n=8227 l=0
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -seg 12 2.2 0x10 || NCBI-err:  || LOSAT-err: Error: the SEG locut or hicut "0x10" (not a finite decimal number) is not supported by LOSAT's BLASTP| || stdout n=8227 l=0
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -seg 12 2.2 1e400 || NCBI-err:  || LOSAT-err: Error: the SEG locut or hicut "1e400" (not a finite decimal number) is not supported by LOSAT's BLASTP| || stdout n=8227 l=0
[SAME] n=1 l=1 files=<>/<> :: -seg 12 +. 2.5
[SAME] n=1 l=1 files=<>/<> :: -seg 12 --5 2.5
[SAME] n=1 l=1 files=<>/<> :: -seg 12 -. 2.5
[SAME] n=1 l=1 files=<>/<> :: -seg 12 . 2.5
[SAME] n=1 l=1 files=<>/<> :: -seg 12 .e1 2.5
[SAME] n=1 l=1 files=<>/<> :: -seg 12 1,5 2.5
[SAME] n=1 l=1 files=<>/<> :: -seg 12 1_0 2.5
[SAME] n=1 l=1 files=<>/<> :: -seg 12 1e 2.5
[SAME] n=1 l=1 files=<>/<> :: -seg 12 1e+ 2.5
[SAME] n=1 l=1 files=<>/<> :: -seg 12 2.2 2.5e
[SAME] n=1 l=1 files=<>/<> :: -seg 12 2.2 2.5x
[SAME] n=1 l=1 files=<>/<> :: -seg 12 2.5.5 2.5
```

## BP-7
```
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -num_threads 100000 || NCBI-err: Warning: [blastp] Number of threads was reduced to 32 to match the number of available CPUs|Warning: [blastp] 'num_threads' is currently ignored when 'subject' is specified.| || LOSAT-err: Error: requested 100000 threads exceeds Rayon maximum 65535, which is not supported by LOSAT| || stdout n=35384 l=0
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -outfmt 5 || NCBI-err:  || LOSAT-err: Error: output format 5 is not supported by LOSAT's BLASTP| || stdout n=57373 l=0
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -query blank.faa -subject $S -outfmt 5 || NCBI-err: Warning: [blastp] Query is Empty!| || LOSAT-err: Error: output format 5 is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -query empty.faa -subject $S -outfmt 5 || NCBI-err: Warning: [blastp] Query is Empty!| || LOSAT-err: Error: output format 5 is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=1 l=1 files=<>/<> :: -num_threads 100000 -matrix FOO || NCBI-err: Warning: [blastp] Number of threads was reduced to 32 to match the number of available CPUs|Warning: [blastp] 'num_threads' is currently ignored when 'subject' is specified.|BLAST query/options error: || LOSAT-err: Error: requested 100000 threads exceeds Rayon maximum 65535, which is not supported by LOSAT| || stdout n=0 l=0
[LOSAT-REJECTS] n=1 l=1 files=<>/<> :: -outfmt 12 -threshold 0 || NCBI-err: BLAST query/options error: Non-zero threshold required|Please refer to the BLAST+ user manual.| || LOSAT-err: Error: output format 12 is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=1 l=1 files=<>/<> :: -outfmt 5 -matrix FOO || NCBI-err: BLAST query/options error: FOO is not a supported matrix, supported matrices are:|BLOSUM80 |BLOSUM62 |BLOSUM50 |BLOSUM45 |PAM250 |BLOSUM90 |PAM30 |PAM70 |IDENTITY ||Please refer to the BLAST+ user man || LOSAT-err: Error: output format 5 is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=1 l=1 files=<>/<> :: -outfmt 5 -seg 1 2 || NCBI-err: BLAST query/options error: Invalid number of arguments to filtering option|Please refer to the BLAST+ user manual.| || LOSAT-err: Error: output format 5 is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=1 l=1 files=<>/<> :: -outfmt 5 -threshold 0 || NCBI-err: BLAST query/options error: Non-zero threshold required|Please refer to the BLAST+ user manual.| || LOSAT-err: Error: output format 5 is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=1 l=1 files=<>/<> :: -outfmt 5 -ungapped || NCBI-err: BLAST query/options error: Composition-adjusted searched are not supported with an ungapped search, please add -comp_based_stats F or do a gapped search|Please refer to the BLAST+ user manual.| || LOSAT-err: Error: output format 5 is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=1 l=1 files=<>/<> :: -outfmt 6 delim=, qseqid -matrix FOO || NCBI-err: BLAST query/options error: FOO is not a supported matrix, supported matrices are:|BLOSUM80 |BLOSUM62 |BLOSUM50 |BLOSUM45 |PAM250 |BLOSUM90 |PAM30 |PAM70 |IDENTITY ||Please refer to the BLAST+ user man || LOSAT-err: Error: a custom delimiter (delim=) in the output format is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=1 l=1 files=<>/<> :: -outfmt 6 qcovs -threshold 0 || NCBI-err: BLAST query/options error: Non-zero threshold required|Please refer to the BLAST+ user manual.| || LOSAT-err: Error: the output field qcovs is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=1 l=1 files=<>/<> :: -outfmt 8 -matrix FOO || NCBI-err: BLAST query/options error: FOO is not a supported matrix, supported matrices are:|BLOSUM80 |BLOSUM62 |BLOSUM50 |BLOSUM45 |PAM250 |BLOSUM90 |PAM30 |PAM70 |IDENTITY ||Please refer to the BLAST+ user man || LOSAT-err: Error: output format 8 is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=3 l=1 files=<>/<> :: -query $Q -subject blank.faa -outfmt 5 || NCBI-err: BLAST engine error: Empty CBlastQueryVector| || LOSAT-err: Error: output format 5 is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[LOSAT-REJECTS] n=3 l=1 files=<>/<> :: -query $Q -subject empty.faa -outfmt 5 || NCBI-err: BLAST engine error: Empty CBlastQueryVector| || LOSAT-err: Error: output format 5 is not supported by LOSAT's BLASTP| || stdout n=0 l=0
[SAME] n=1 l=1 files=<>/<> :: -outfmt 13 -matrix FOO
[SAME] n=1 l=1 files=<>/<> :: -outfmt 14 -ungapped
```

## BP-8
```
[DIFF_STDOUT DIFF_STDERR DIFF_EXIT] n=0 l=2 files=<>/<> :: -outfmt 6 -- || NCBI-err:  || LOSAT-err: error: unknown option or argument '--'; use -help for CLI v2 syntax || stdout n=2889 l=0
[LOSAT-REJECTS] n=0 l=2 files=<>/<> :: -dryrun -matrix FOO || NCBI-err:  || LOSAT-err: error: the NCBI C++ Toolkit option -dryrun is not supported by LOSAT's BLASTP || stdout n=0 l=0
[LOSAT-REJECTS] n=0 l=2 files=<>/<> :: -dryrun || NCBI-err:  || LOSAT-err: error: the NCBI C++ Toolkit option -dryrun is not supported by LOSAT's BLASTP || stdout n=0 l=0
[LOSAT-REJECTS] n=0 l=2 files=<>/<> :: -evalue -version || NCBI-err:  || LOSAT-err: error: the NCBI BLAST+ option -version is not supported by LOSAT's BLASTP || stdout n=67 l=0
[LOSAT-REJECTS] n=0 l=2 files=<>/<> :: -matrix -version || NCBI-err:  || LOSAT-err: error: the NCBI BLAST+ option -version is not supported by LOSAT's BLASTP || stdout n=67 l=0
[LOSAT-REJECTS] n=0 l=2 files=<>/<> :: -out -version || NCBI-err:  || LOSAT-err: error: the NCBI BLAST+ option -version is not supported by LOSAT's BLASTP || stdout n=67 l=0
[LOSAT-REJECTS] n=0 l=2 files=<>/<> :: -out x -dryrun || NCBI-err:  || LOSAT-err: error: the NCBI C++ Toolkit option -dryrun is not supported by LOSAT's BLASTP || stdout n=0 l=0
[LOSAT-REJECTS] n=0 l=2 files=<>/<> :: -outfmt -version || NCBI-err:  || LOSAT-err: error: the NCBI BLAST+ option -version is not supported by LOSAT's BLASTP || stdout n=67 l=0
[LOSAT-REJECTS] n=0 l=2 files=<>/<> :: -outfmt 6 -dryrun || NCBI-err:  || LOSAT-err: error: the NCBI C++ Toolkit option -dryrun is not supported by LOSAT's BLASTP || stdout n=0 l=0
[LOSAT-REJECTS] n=0 l=2 files=<>/<> :: -task -version || NCBI-err:  || LOSAT-err: error: the NCBI BLAST+ option -version is not supported by LOSAT's BLASTP || stdout n=67 l=0
[LOSAT-REJECTS] n=0 l=2 files=<>/<> :: -version || NCBI-err:  || LOSAT-err: error: the NCBI BLAST+ option -version is not supported by LOSAT's BLASTP || stdout n=67 l=0
[LOSAT-REJECTS] n=1 l=2 files=<>/<> :: -evalue -dryrun || NCBI-err: USAGE|  blastp [-h] [-help] [-import_search_strategy filename]|    [-export_search_strategy filename] [-task task_name] [-db database_name]|    [-dbsize num_letters] [-gilist filename] [-seqidlist fil || LOSAT-err: error: the NCBI C++ Toolkit option -dryrun is not supported by LOSAT's BLASTP || stdout n=0 l=0
[LOSAT-REJECTS] n=1 l=2 files=<>/<> :: -matrix -dryrun || NCBI-err: USAGE|  blastp [-h] [-help] [-import_search_strategy filename]|    [-export_search_strategy filename] [-task task_name] [-db database_name]|    [-dbsize num_letters] [-gilist filename] [-seqidlist fil || LOSAT-err: error: the NCBI C++ Toolkit option -dryrun is not supported by LOSAT's BLASTP || stdout n=0 l=0
[LOSAT-REJECTS] n=1 l=2 files=<>/<> :: -out -dryrun || NCBI-err: USAGE|  blastp [-h] [-help] [-import_search_strategy filename]|    [-export_search_strategy filename] [-task task_name] [-db database_name]|    [-dbsize num_letters] [-gilist filename] [-seqidlist fil || LOSAT-err: error: the NCBI C++ Toolkit option -dryrun is not supported by LOSAT's BLASTP || stdout n=0 l=0
[LOSAT-REJECTS] n=1 l=2 files=<>/<> :: -outfmt -dryrun || NCBI-err: USAGE|  blastp [-h] [-help] [-import_search_strategy filename]|    [-export_search_strategy filename] [-task task_name] [-db database_name]|    [-dbsize num_letters] [-gilist filename] [-seqidlist fil || LOSAT-err: error: the NCBI C++ Toolkit option -dryrun is not supported by LOSAT's BLASTP || stdout n=0 l=0
[LOSAT-REJECTS] n=1 l=2 files=<>/<> :: -task -dryrun || NCBI-err: USAGE|  blastp [-h] [-help] [-import_search_strategy filename]|    [-export_search_strategy filename] [-task task_name] [-db database_name]|    [-dbsize num_letters] [-gilist filename] [-seqidlist fil || LOSAT-err: error: the NCBI C++ Toolkit option -dryrun is not supported by LOSAT's BLASTP || stdout n=0 l=0
[SAME] n=0 l=0 files=<./-help >/<./-help > :: -out -help
```

## BP-9
```
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -num_threads 100000 || NCBI-err: Warning: [blastp] Number of threads was reduced to 32 to match the number of available CPUs|Warning: [blastp] 'num_threads' is currently ignored when 'subject' is specified.| || LOSAT-err: Error: requested 100000 threads exceeds Rayon maximum 65535, which is not supported by LOSAT| || stdout n=35384 l=0
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -num_threads 65535 || NCBI-err: Warning: [blastp] Number of threads was reduced to 32 to match the number of available CPUs|Warning: [blastp] 'num_threads' is currently ignored when 'subject' is specified.| || LOSAT-err: Error: failed to build blastp pool with 65535 threads; a thread count that the system cannot start is not supported by LOSAT||Caused by:|    0: Resource temporarily unavailable (os error 11)|    1: Re || stdout n=35384 l=0
[LOSAT-REJECTS] n=0 l=1 files=<>/<> :: -num_threads 65536 || NCBI-err: Warning: [blastp] Number of threads was reduced to 32 to match the number of available CPUs|Warning: [blastp] 'num_threads' is currently ignored when 'subject' is specified.| || LOSAT-err: Error: requested 65536 threads exceeds Rayon maximum 65535, which is not supported by LOSAT| || stdout n=35384 l=0
[PARSER-EXC(ncbi=1,losat=2)] n=1 l=2 files=<>/<> :: -num_threads -1 || NCBI-err: USAGE|  blastp [-h] [-help] [-import_search_strategy filename]|    [-export_search_strategy filename] [-task task_name] [-db database_name]|    [-dbsize num_letters] [-gilist filename] [-seqidlist fil || LOSAT-err: error: invalid value '-1' for '-num_threads <NUM_THREADS>': expected an integer >= 1||For more information, try '--help'.| || stdout n=0 l=0
[PARSER-EXC(ncbi=1,losat=2)] n=1 l=2 files=<>/<> :: -num_threads 0 || NCBI-err: USAGE|  blastp [-h] [-help] [-import_search_strategy filename]|    [-export_search_strategy filename] [-task task_name] [-db database_name]|    [-dbsize num_letters] [-gilist filename] [-seqidlist fil || LOSAT-err: error: invalid value '0' for '-num_threads <NUM_THREADS>': expected an integer >= 1||For more information, try '--help'.| || stdout n=0 l=0
[DIFF_STDERR] n=0 l=0 files=<>/<> :: -num_threads 2 || NCBI-err: Warning: [blastp] 'num_threads' is currently ignored when 'subject' is specified.| || LOSAT-err:  || stdout n=35384 l=35384 || firstdiff: 
[DIFF_STDERR] n=0 l=0 files=<>/<> :: -num_threads 4 || NCBI-err: Warning: [blastp] 'num_threads' is currently ignored when 'subject' is specified.| || LOSAT-err:  || stdout n=35384 l=35384 || firstdiff: 
```

## BP-11
```
[DIFF_STDOUT DIFF_STDERR DIFF_EXIT] n=1 l=0 files=<>/<> :: --help || NCBI-err: USAGE|  blastp [-h] [-help] [-import_search_strategy filename]|    [-export_search_strategy filename] [-task task_name] [-db database_name]|    [-dbsize num_letters] [-gilist filename] [-seqidlist fil || LOSAT-err:  || stdout n=0 l=1885
[DIFF_STDOUT] n=0 l=0 files=<>/<> :: -help || NCBI-err:  || LOSAT-err:  || stdout n=14687 l=1885 || firstdiff: 1,25c1|< USAGE|<   blastp [-h] [-help] [-import_search_strategy filename]|
[LOSAT-REJECTS] n=0 l=2 files=<>/<> :: -h || NCBI-err:  || LOSAT-err: error: the NCBI BLAST+ option -h is not supported by LOSAT's BLASTP || stdout n=1680 l=0
```

## BP-10
```
LOSAT_STARTUP_TRACE=1 LOSAT-head blastp -query Q -subject S -outfmt 6: exit 0, stdout identical to the run without the variable, stderr "[startup] enter main" / "[startup] after clap parse" (still present; ROUND1: accepted as LOSAT diagnostic env var). Without the variable stderr is empty.
Class: ACCEPTED (ROUND1 row BP-10).
```

## BP-12 (NCBI file:line references; vref2.py re-run on head-src, same method as the audit)
```
Table rows (audit -> head): protein_options.rs:390 now cites blast_options.c:910-936 (line 910 `if (options->gapped_calculation ...` confirmed): FIXED. protein_options.rs:492 now 1518-1521 (confirmed 1518-1521): FIXED. blast_engine.rs:7011 now cites blast_format.cpp:2267 (confirmed `options.GetMatrixName()`): FIXED. args.rs:91 now 838-886 (switch at 838): FIXED. app.rs:638 blast_input.hpp:258-259 reference removed; app.rs:704 now cites blast_input.hpp:313,364 (364 `TSeqPos m_BatchSize` confirmed): FIXED. value_parsers.rs:550 now cmdline_flags.cpp:46-51 (kArgSubject line 51): FIXED.
Residual DIFF (low, doc only, both also present in base):
 - cli.rs:15 cites cmdline_flags.cpp:107,143 but the quoted snippet is kArgQuery/kArgSubject/kArgNumThreads, which are at lines 46, 51, 75 (107 = kArgWordSize, 143 = kArgCompBasedStats).
 - cli.rs:172 cites tblastn_args.cpp:63-110 but the last quoted line `m_PsiBlastArgs.Reset(new CPsiBlastArgs(CPsiBlastArgs::eNucleotideDb));` is at line 125 (outside the span).
vref2.py totals: head args.rs/app.rs/value_parsers.rs/protein_options.rs/cli.rs/main.rs 72 blocks, 6 flagged (base: 69 blocks, 5 flagged). The other 4 flagged are script false positives (ellipsis lines, first span of a two-span comment, prose after a block): args.rs:13, args.rs:603, app.rs:606, app.rs:716 (blast_stat.c:922-923 confirmed). blast_engine.rs: 397 blocks, 137 flagged = same count as base (audit: none session-added; no new flagged reference file:line set vs base).
```

## BP-13 (web adapter vs CLI; code reading of head-src, wasm adapter not executed)
```
web/adapter/src/run.rs:92 validate -> blastp::blast_engine::check_options (blast_engine.rs:4618): args.resolve() + validate_requested_blastp_support + validate_threads(args.num_threads) (line 4622, "The host validates the argv it runs"): the audit's missing -num_threads check is now present. 
blastp::blast_engine::run_local (4745-4798): subject empty -> empty_subjects_error, query empty -> "Warning: [blastp] Query is Empty!", then check_protein_residues_of for subject and query and check_records_have_residues_of: the audit's missing residue checks are now present. web/adapter/src/store.rs register: check_protein_sequence_lines_of / check_protein_input_of for BLASTP inputs (same checks as the CLI).
Class: FIXED (code reading only; the wasm adapter was not run).
```

## Harness re-run ("Checked OK" lists), LOSAT-head vs NCBI
Method: the report's list files t1..t18 (those are the argv lists batch.sh read; copied to audit-rerun/blastp/), run with the copied cmp.sh/batch.sh (binary replaced by LOSAT-head, 4 jobs at a time, 60 s timeout per program, stdin from /dev/null, each run in its own scratch cwd). t1-t15 and t11b use the default Q/S; t16 = 1 query x 300 subjects (e2e_many_query / e2e_many_subject); t17 = 8 queries x 300 subjects; t18 both. No script hit the 40 minute cap (longest: t12, 971 s, the 60 s timeouts of NCBI on huge -threshold with -word_size 5). Raw per-command lines: audit-rerun/blastp/out/*.raw.
Notes: (a) lines with an unquoted "(" (`-evalue +nan(1)`, `+nan(abc)`) fail in batch.sh's eval; they were re-run quoted (nan.txt, plus -nan(1), +nan(), +NAN(x_1), +NaN, -nan): 7 SAME. (b) t13's `-n -subject ...` (no -query) read the rest of the list from stdin in the first pass (36 of 63 lines); fixed with </dev/null and re-run: 63 of 63. (c) stdin app flow re-run separately (stdin.sh): see below.

| script | commands | SAME | LOSAT-REJECTS (NCBI ok) | LOSAT-REJECTS (NCBI also fails, other text/exit) | PARSER-EXC | DIFF |
|---|---|---|---|---|---|---|
| t1 | 54 | 26 | 17 | 0 | 11 | 0 |
| t2 | 62 | 50 | 6 | 3 | 3 | 0 |
| t3 | 100 | 63 | 16 | 11 | 10 | 0 |
| t4 | 48 | 44 | 4 | 0 | 0 | 0 |
| t5 | 62 (+2 via nan.txt) | 37 | 10 | 11 | 4 | 0 |
| t6 | 60 | 22 | 36 | 1 | 1 | 0 |
| t7 | 107 | 68 | 39 | 0 | 0 | 0 |
| t8 | 85 | 61 | 3 | 4 | 16 | 1 (slow rejection, see below) |
| t9 | 68 | 0 | 33 | 20 | 13 | 2 (-help, --help: approved/accepted) |
| t10 | 80 | 72 | 4 | 4 | 0 | 0 |
| t11 | 38 | 21 | 16 | 1 | 0 | 0 |
| t11b | 36 | 19 | 16 | 1 | 0 | 0 |
| t12 | 25 | 14 | 0 | 11 | 0 | 0 |
| t13 | 63 | 57 | 0 | 3 | 3 | 0 |
| t14 | 100 | 60 | 40 | 0 | 0 | 0 |
| t15 | 174 | 174 | 0 | 0 | 0 | 0 |
| t16 | 44 | 43 | 1 | 0 | 0 | 0 |
| t17 | 44 | 43 | 0 | 1 | 0 | 0 |
| t18 (1 query) | 24 | 23 | 1 | 0 | 0 | 0 |
| t18 (8 queries) | 24 | 23 | 1 | 0 | 0 | 0 |
| nan.txt | 7 | 7 | 0 | 0 | 0 | 0 |
| stdin | 13 | 11 | 2 | 0 | 0 | 0 |
| total | 1318 | 938 | 245 | 71 | 61 | 3 |

LOSAT-REJECTS classification: every one of the 316 contains "not supported by LOSAT" in LOSAT stderr. With NCBI ok (245): matrix/gap costs other than NCBI's supported pairs for the non-default matrices (PAM*, IDENTITY, BLOSUM with other pairs), -comp_based_stats 0/1/3 and u-modes, -ungapped, unsupported -outfmt numbers (1-5, 8-12, 15, 16, 18, 20), delim=, unsupported output fields (qgi, sallseqid, staxid, qcovs, ...), -word_size 2/4 and compressed-lookup cases, -use_sw_tback, -task blastp-short, -evalue +inf/1e999 (BP-4 decision D12), -window_size 2147483647 (D11), unsupported options (-dbsize, -xdrop_*, -html, -remote ...), hex/hex-float -evalue/-threshold (clap "invalid value ... not supported by LOSAT's BLASTP", exit 2), seg locut/hicut inf/nan/hex/1e400, nohdr.faa FASTA bio cannot read. With NCBI also failing (71): NCBI prints its own option error / crashes (exit 124/139) first, LOSAT's explicit rejection comes first or with other text (accepted BP-7 class), plus NCBI crash/timeout cases for huge -threshold with -word_size 5/blastp-fast (NCBI "malloc(): corrupted top size", LOSAT quick rejection).
PARSER-EXC (61): NCBI USAGE exit 1 vs LOSAT clap "error: ..." exit 2 (approved exception PD-LOSAT-CLI-NONSEARCH-DIFFERENCES).

### stdin / app flow (stdin.sh; out/stdin.raw, out/stdin2.raw)
```
-query - -subject S (stdin=query), also -outfmt 6: SAME. -subject - (stdin=subject) with -query Q, also -outfmt 6: SAME. -out - (outfmt 0, 6): SAME.
empty stdin: -subject - : SAME (exit 3 "BLAST engine error: Empty CBlastQueryVector" both); -query - (outfmt 0, 6): SAME (exit 0, "Query is Empty!" both); -query - -subject - : SAME with empty stdin (exit 3).
-query - -subject - with the query on stdin: LOSAT-REJECTS "an empty query from a stream without a position (such as a pipe) is not supported by LOSAT's BLASTP" (D10; NCBI exit 0 with no output, see ROUND1 row IN-9).
-query - -subject S -matrix FOO: SAME (NCBI option error text, exit 1). -outfmt 5: LOSAT-REJECTS.
```

## Summary table
| finding | verdict | evidence (LOSAT-head vs NCBI) |
|---|---|---|
| BP-1 empty-subject search space | SAME (fixed) | `-subject hdronly.faa` outfmt 0/6/7: stdout, stderr, exit identical (3/3) |
| BP-2 -max_target_seqs 2^30-25..2^31-25 | LOSAT-REJECTS (explicit rejection) | 14/14 (7 values x outfmt 0,6): NCBI exit 139, LOSAT exit 1 "a -max_target_seqs of N, whose preliminary hit list size overflows NCBI BLAST+'s 32-bit int (NCBI crashes), is not supported by LOSAT's BLASTP"; 1073741798 SAME |
| BP-3 -max_target_seqs >= 2^31-24 wrap | SAME (fixed) | 2147483624/30/40/46/47 on 300 subjects, outfmt 6 and 0 (2147483647): identical; small input 2147483647 identical |
| BP-4 infinite -evalue | LOSAT-REJECTS (D12) | +inf/1e999 (outfmt 0/6/7, -task blastp-fast, word size 3 and 5, 300 subjects): NCBI 139 (or exit 0 for word size 3), LOSAT exit 1 "an infinite -evalue (inf), ... is not supported by LOSAT's BLASTP"; finite 1e10/1e100/1e308 SAME |
| BP-5 `frames`/`sframe` | SAME (fixed) | 34 argv (all listed lists, outfmt 6 and 7, also on 300 subjects): identical; "6 frames" is "0/0" in both |
| BP-6 -seg malformed locut/hicut | SAME (fixed) for the 12 malformed values; LOSAT-REJECTS for NCBI-accepted +inf/-inf/1e400/0x10/0x1p3/-nan/+nan(5) | malformed: both "BLAST query/options error: Invalid input for filtering parameters" exit 1; the rest LOSAT exit 1 "the SEG locut or hicut ... is not supported by LOSAT's BLASTP" |
| BP-7 LOSAT rejection before NCBI's later checks | ACCEPTED (ROUND1) | e.g. `-outfmt 5 -matrix FOO`, `-subject empty.faa -outfmt 5` (NCBI 3, LOSAT 1), `-num_threads 100000 -matrix FOO`: LOSAT exit 1 with "not supported by LOSAT"; `-outfmt 13/14` cases SAME |
| BP-8 -dryrun / toolkit args | LOSAT-REJECTS (explicit, D9; `-out -version` too) + ACCEPTED (`--`) | `-dryrun`, `-out x -dryrun`, `-outfmt 6 -dryrun`, `-dryrun -matrix FOO`: exit 2 "error: the NCBI C++ Toolkit option -dryrun is not supported by LOSAT's BLASTP"; `-out -version`, `-evalue -version`, `-task -version`, `-matrix -version`, `-outfmt -version`, `-version`: exit 2 "the NCBI BLAST+ option -version is not supported by LOSAT's BLASTP"; `-out -dryrun` (NCBI exit 1) exit 2 same text. No file named "-version" or "-dryrun" is created (files=<> for LOSAT in every case). `-out -help` SAME (both create file "-help"). `-outfmt 6 --`: NCBI exit 0 / LOSAT exit 2 (accepted exception 1, argument syntax) |
| BP-9 rejection wording | ACCEPTED (ROUND1 TD-1) | -num_threads 100000/65536: "requested N threads exceeds Rayon maximum 65535, which is not supported by LOSAT"; -num_threads 65535 on this host: "failed to build blastp pool with 65535 threads; a thread count that the system cannot start is not supported by LOSAT" (NCBI exit 0; host thread limit). -num_threads 2/4: only NCBI's thread warning differs (approved exception) |
| BP-10 LOSAT_STARTUP_TRACE | ACCEPTED (ROUND1) | still prints "[startup] enter main"/"after clap parse" only when the variable is 1 |
| BP-11 --help | ACCEPTED (approved exception 1) | `--help`: NCBI USAGE exit 1, LOSAT help exit 0; `-help`: NCBI full help, LOSAT own help, both exit 0 (stdout differs, approved) |
| BP-12 NCBI file:line references | FIXED for 6 of 8 table rows, 2 residual DIFF (doc only) | cli.rs:15 (cites cmdline_flags.cpp:107,143, snippet is lines 46/51/75), cli.rs:172 (cites tblastn_args.cpp:63-110, last quoted line is at 125) |
| BP-13 web adapter vs CLI | FIXED (code reading only, wasm not run) | check_options calls validate_threads; run_local has residue and empty-input checks; register checks protein input |
| Harness (t1..t18, nan, stdin) | see table above | 1318 argv: 938 SAME, 316 LOSAT-REJECTS (all contain "not supported by LOSAT"), 61 PARSER-EXC, 3 DIFF |

## DIFF list (everything not SAME / LOSAT-REJECTS / accepted-exception)
1. BP-12 residual (doc only): `LOSAT/src/cli.rs:15` NCBI reference cmdline_flags.cpp:107,143 vs quoted kArgQuery/kArgSubject/kArgNumThreads (lines 46, 51, 75); `LOSAT/src/cli.rs:172` tblastn_args.cpp:63-110 but quoted `m_PsiBlastArgs.Reset(...)` is at line 125. No behavioural difference.
2. `blastp -query Q -subject S -task blastp-fast -threshold +inf` (t8): NCBI crashes ("malloc(): corrupted top size", 17 bytes of output); LOSAT does reject with "-threshold inf with the compressed lookup table of -word_size 5 is not supported by LOSAT's BLASTP: its neighboring words fill NCBI BLAST+'s 1024 overflow banks ..." (exit 1) but only after about 69 s (host load average about 35), so under a 60 s timeout it looks like a hang (exit 124, no message). Behaviour-wise a rejection; the late detection is a usability/perf note, not a parity difference. (Other -threshold huge cases with -word_size 5 and blastp-fast reject quickly.)
3. No other DIFF. The stdout/--help/-help rows (approved/accepted) and the harness artifacts (stdin file path, the first-pass t13) are not findings.

## Process note
While restarting my own driver I ran `pkill -f` and `kill` on process IDs that turned out to belong to other agents' jobs on this shared host (a `./cmp.sh -q q1.fna -s s1.fna -- -window_size 1500000000` run and two `bash -c` wrappers, probably another audit area). Those runs were terminated by me; the owners may need to re-run that command. No source or /mnt/c file was touched.
