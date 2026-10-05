# Auditor D report (S11 query_loc/subject_loc)
Running count of comparisons: see bottom lines "COUNT".
COUNT forms.cases 256 (176 DIFF all parser-syntax rc1/2 or R1 a-b)
COUNT spell_query 312, spell_subject 312 (all DIFF are R1 255/1)
COUNT edge.cases 2079 + edge_tn_subj 297 (outfmt 0/6/7; all DIFF are R2 and start==L+1 exactly; no other differences)
COUNT matrix.cases 18108 (error-order matrix: 4 programs x 70 other-option/error variants x 9 query_loc x 7 subject_loc states). DIFF classes: parser-USAGE (rc1/2) 8338, R1 1189, R2 59 (start==L+1), pre-existing explicit rejections (not supported by LOSAT; outfmt/custom fields/reward/comp_based_stats/records w/o residues) 686+1520, thread warnings 60 (approved). No unexplained differences.
COUNT tiny.cases 2034 (tiny records 1..12 nt/aa, CRLF/nofinal/width/blank/lower/ids/dup ids/long titles formats; all DIFF = R2 exact (chkr2.py) or pre-existing empty-defline rejection)
COUNT stdin.cases 664 (stdin query/subject, subject_loc without -subject, -db; DIFF = R1, R2, usage, pre-existing '-query - -subject -' rejection and non-nucleotide-subject rejection)
COUNT warn.cases 6930 (title warnings, Query_<n> numbering, invalid-query warnings, skipped records, BATCH_SIZE 450/900/default, subject range errors after title warnings; all DIFF = R2 exact, 1242)
COUNT spell2.cases 496 (more number spellings; all DIFF = R1)
COUNT batch.cases partial 1447 + batch2 1500 (BATCH_SIZE 1..5000 x skipped-record sets x outfmt 0/7, 4 programs; all DIFF = R2 exact)
COUNT fields.cases 742 (blastp custom fields x range states; DIFF only unsupported-field rejections, pre-existing)
COUNT order2s.cases 3000 sampled of 6480 (warning-order: few-matches warning, thread warnings, -evalue 0 x query/subject range states x title-warning files); DIFF = R1, R2, thread warnings only (thrcheck.py verifies stdout/rc/rest of stderr equal)
COUNT emptyin.cases 1152 (empty/blank/whitespace/header-only query or subject x ranges; DIFF = R1 228 + pre-existing records-without-residues rejection 33)
COUNT web harness (adapter native build of the audited source: validate+register+run vs CLI, streams 0/6/7 + diagnostics): web6k 5477 cases x 3 formats; with NOVALIDATE 5355 equal; 122 differ only in register-time pre-existing rejections or non-threaded-build thread check; with validate 321 differ (validate reports option/syntax errors before record-dependent subject range errors)
COUNT both.cases 2400 (query_loc x subject_loc edge combos, random 600/program; DIFF = R2 only), oxw.cases 384 (O->X warning positions with ranges + Query_<n>; DIFF = R2 only), threads.cases 212 (e2d fixtures x -num_threads 1/2/4/8; DIFF = NCBI thread warnings only)
COUNT outfile.cases 1296 (-out file contents/existence with ranges, also unwritable out; DIFF = R1, R2 only)
COUNT split.cases 62 (chunk-boundary ranges BLASTP 10000/TBLASTN 20000; all SAME)

## FINDING D-1 (all four programs, low): non-UTF-8 bytes in a -query_loc / -subject_loc value
Command (any program, either role): `blastn -query nq300.fa -subject nsub.fa -query_loc $'\xff1-5'` (also `$'1-\xff5'`, `$'\xc3\x28-5'`, same with -subject_loc, and blastp/tblastn/tblastx).
Inputs: scratch_D/nq300.fa, nsub.fa (any valid pair).
NCBI 2.17.0: exit 255, `Error: NCBI C++ Exception: ... ncbi::NStr::StringToInt8() - Cannot convert string '\3771' to Int8` (the R1 class: ParseSequenceRange -> NStr::StringToInt, blast_input_aux.cpp:145-179).
LOSAT (a262c8726): exit 2, stderr `error: invalid UTF-8 was detected in one or more arguments` + clap usage. Not an R1 text ("... cannot convert to an int ...; not supported by LOSAT's <PROGRAM>", exit 1), and not NCBI's USAGE (so not the approved parser-syntax exception either: NCBI is not a USAGE error here).
Cause: clap parses `query_loc`/`subject_loc` as `Option<String>` (LOSAT/src/algorithm/*/args.rs, `pub query_loc: Option<String>`), so the byte string never reaches `parse_sequence_range`/`ncbi_string_to_int(OsStr)` (seq_range.rs `parse_sequence_range`, which takes &str). Other non-UTF-8 values that NCBI converts the same way (e.g. U+FFFE encoded in valid UTF-8) do reach R1 correctly (rc 1 R1 message).
Severity: low (only non-UTF-8 argv bytes; NCBI also fails, with exit 255; LOSAT fails with exit 2 and a different message).
COUNT web2.cases 1914 + web3.cases 792 (adapter harness; 74 invalid generator cases excluded) all OK; both_tn.cases 800 (R2 only)

## Observations (not findings)
- O-1 (Web ABI v2 `validate`, informational): `validate` checks the range syntax and the options before it has any record, so when a record-dependent subject range error ("Invalid from coordinate (greater than sequence length)"), "Empty CBlastQueryVector" or an R2 error would come first in NCBI/the CLI, `validate` reports the earlier-in-argv-syntax/option error instead (321 of 5477 harness cases with validate, all of this form; with `validate` skipped, `run` matches the CLI). abi_v2.md section 4 (`validate` row) states "The checks of the inputs come with register and run" but does not list the range syntax checks. `run` itself keeps NCBI's order.
- O-2 (Web `register`, informational): `register` rejects records without residues and empty deflines (pre-existing LOSAT rejections) before `run`, whereas the CLI defers those checks so that NCBI's range errors win (e.g. `-subject hdr_only.fa -subject_loc 50000-60000`: CLI prints NCBI's "Invalid from coordinate", web `register` prints the LOSAT rejection). Pre-existing class, not range-specific.
- Web run_local (blastn/blastp/tblastn/tblastx) was compared by building the audited adapter natively (scratch_D/adbuild, harness = validate + register + run, streams 0/6/7 and diagnostics) against the CLI: every case where `register` succeeded matched byte for byte (stdout streams per format, stderr = diagnostics stream; on errors the message equals the last CLI stderr lines, the web discards warnings that preceded an error by design).
- Native harness build has no parallel feature, so `-num_threads > 1` is "unsupported ..." there; not a property of the wasm threads build.
- ABI v1 (web_api.rs) rejects `-query_loc`/`-subject_loc` as "unsupported <program> argument for web API" (frozen, explicit).

## Totals
NCBI-vs-LOSAT comparisons (one invocation pair each, outfmt 0/6/7 spread over the cases): about 77,800. Web adapter harness comparisons: about 18,400 cases x 3 formats vs the CLI.
All DIFFs outside D-1 were classified by script into: argument-parser USAGE (NCBI 1 / LOSAT 2), R1 (255/1, message checked), R2 (start == record length + 1 exactly, role and record index in the message verified against the files for 1,329 fuzz cases, NCBI's own result for all of them is the "Sequence contains no data" family), NCBI thread warnings (stdout, rc and the rest of stderr verified equal), and pre-existing explicit rejections (outfmt 5/8/..., custom fields for BLASTN/TBLASTN/TBLASTX and unsupported BLASTP fields, -comp_based_stats 0, -reward 0, matrices, tasks, records without residues, empty deflines, `-query - -subject -`, BATCH_SIZE values NCBI cannot convert).

## CONCLUSION (perspective d)
Strictly by the audit rule: unsupported, because of one confirmed low-severity finding (D-1: non-UTF-8 bytes in a range value give exit 2 with clap text where NCBI gives exit 255 and R1 promises exit 1 with the "cannot convert" text). Apart from D-1 no difference was found in about 77,800 CLI comparisons and 18,400 adapter comparisons covering spellings, order of errors, record edges (L-1..L+2 starts, L-1..L+1 ends, first/middle/last record, lengths 1..5000), title warnings, Query_<n>/Subject_<n> numbering, BATCH_SIZE, stdin, -out files, thread counts and the Web adapter's run path; if D-1 is folded into R1 the claim is supported.
