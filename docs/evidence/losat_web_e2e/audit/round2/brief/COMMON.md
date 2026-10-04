# S08+b (E2e) independent audit, round 2 (ncbi_parity_auditor role): common brief

You are an independent, read-only auditor of LOSAT's BLASTP, TBLASTN and TBLASTX search options and reports against NCBI BLAST+ 2.17.0 (blastp, tblastn, tblastx). LOSAT is a Rust reimplementation of NCBI BLAST+. You do not change LOSAT. Write everything you produce under your own work dir (given in your angle brief), appending results to files as you go (if writing a `.md` file is refused, use `FINDINGS.txt`).

## Hard rules

- Never read, list or run anything under `/mnt/c` (the WSL 9p mount fails under load). Everything you need is copied below `/home/kawato/.cache/losat-web-gui-target/`.
- Never use `pkill`, `killall`, or `kill` by name or with `-f`/pattern matching. You may only stop a process whose PID you started yourself (use `timeout` on every run instead). Other agents and a gate run share this machine.
- Do not run `git` anywhere, do not edit any file outside your work dir, do not create worktrees or clones.
- Use at most 3 concurrent processes. Before every NCBI or LOSAT run, wait while `/home/kawato/.cache/losat-web-gui-target/vperf.lock` exists (`while [ -e .../vperf.lock ]; do sleep 60; done`).
- Never use NCBI's `-remote` (it contacts NCBI). Avoid NCBI runs that need gigabytes: blastp/tblastn `-word_size` 5-7 with `-threshold` below 5, tblastx `-word_size 4 -threshold 1` (the gate's option sweep already covers them); keep every run under `timeout 300` unless your brief says otherwise, and use inputs of a few kb to ~60 kb.

## What is under audit

- LOSAT binary: `FINAL=/home/kawato/.cache/losat-web-gui-target/s08pb-audit/bin/LOSAT` (a copy of the final gate's native build of commit HEAD of branch feature/losat-web-gui; SHA-256 in `bin/LOSAT.sha256`). Run it as `$FINAL blastp|tblastn|tblastx <NCBI-style args>`.
- Source of that commit (read-only copy): `/home/kawato/.cache/losat-web-gui-target/s08pb-audit/src/` (`LOSAT/src/...`, `web/adapter/...`, `docs/...`).
- What changed: `ref/session_commits.txt` (all E2e commits since 78c06fe61), `ref/s08pa_s08pb.diff` (S08+a and S08+b engine changes since the first audit's fixes, d96412265..HEAD), `ref/s08pb_stable_sorts.diff`, `ref/s08p_all.diff` (the whole E2e engine diff since 78c06fe61).
- Records (in the source copy, `docs/evidence/losat_web_e2e/`): `AUTHORITY.md` (NCBI path tables; §M decisions D1-D12; §N items that round 1 left open), `audit/ROUND1.md` (round 1 findings and how each was handled), `audit/round1/{blastp,tblastn,tblastx,reports,inputs}.md` (round 1 reports with reproducing commands), `s08pa/NOTES.md` (root causes and fixes of TN-1 residual, TN-2, TN-4, TN-5, RP-4 query splitting, `-out -version`), `s08pa/audit_rerun/*.md` (S08+a's rerun of every round 1 command with an intermediate binary). Fixtures: `LOSAT/tests/outfmt0_manifest.tsv` (`e2e.*` rows; frozen NCBI outputs in `LOSAT/tests/fixtures/outfmt0/`, inputs in `LOSAT/tests/fasta/outfmt0/`).
- Round 1 harnesses that S08+a moved off /mnt/c: `/home/kawato/.cache/losat-web-gui-target/s08pa/audit-rerun/<angle>/` (scripts `cmp.sh`, `cmp2.sh`, `c.sh`, `scripts/*.py`, inputs). They hard-code `LOSAT=.../s08pa/bin/LOSAT-head` and some write under their own dir. Copy what you need into your work dir and point every LOSAT path at `$FINAL`; never write into the `s08pa/` tree. Inputs of round 1 are under `/home/kawato/.cache/losat-web-gui-target/s08p/audit/` (read-only).
- NCBI source (fixed commit 598d8ae6, LF copy): `/home/kawato/.cache/losat-web-gui-target/s08p/ncbi/c++` — the only source of truth. Cite file:line for every claim.
- NCBI binaries (comparison oracle only): `/home/kawato/micromamba/bin/{blastp,tblastn,tblastx,makeblastdb}`. The comparison-only C++ API oracle for TBLASTN subject genetic codes: `/home/kawato/.cache/losat-web-gui-target/s08p/api-oracle/tblastn_stage_e_local_oracle` (usage in `docs/evidence/losat_web_e2e/gencode_api_check.py`).

## Approved exceptions and decisions (not defects)

- TBLASTX and TBLASTN apply a non-default `-db_gencode` to the subject in local `-subject` searches (NCBI's local search uses code 1). For TBLASTN the oracle is the C++ API oracle above; for TBLASTX, NCBI with the subject as a BLAST database (`makeblastdb -dbtype nucl` without `-parse_seqids`, then `-db`), outside the database lines.
- `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`: argument-parser syntax errors and `-help`/`-h` use LOSAT's parser text and exit 2 (NCBI: USAGE, exit 1); with `-subject` LOSAT honors `-num_threads` and prints no thread warnings; outfmt 6/7 write failures exit non-zero (NCBI aborts); allocation failure aborts; an outfmt 0 write to a closed pipe exits 6; a standard output closed at start is discarded.
- `PD-LOSAT-NCBI-DEFECTS`: BLASTN searches a query chunk that NCBI would split again once; outfmt 0 titles made only of punctuation stop at the end of the string (BLASTN, TBLASTX, TBLASTN). Deterministic NCBI results are reproduced; NCBI failures without a checkable valid result are explicit rejections.
- Explicit rejections ("... is not supported by LOSAT's <PROGRAM>" or "... not supported by LOSAT") are not parity defects when the stated reason holds (the option or input is really not ported, or NCBI crashes/hangs there). They are defects when LOSAT rejects something NCBI handles and LOSAT claims to support, or when LOSAT silently prints a different result instead of an error.
- Pending the maintainer (treat as accepted explicit rejections, report only if the rejection is wrong or incomplete): D11 (query length + `-window_size` above 2^30 rejected), D12 (infinite `-evalue` rejected for BLASTP and TBLASTN), S08+a's two choices: NCBI-crashing query-split settings (`CHUNK_SIZE`/`OVERLAP_CHUNK_SIZE` that would split a chunk again, negative `CHUNK_SIZE` that splits a batch) are explicit rejections for BLASTP and TBLASTN; NCBI toolkit words in an option's value position (`-out -version`, `-out -dryrun`, ...) are explicit rejections.
- Accepted in round 1 (ROUND1.md): TX-5/BP-7/TN-9 (LOSAT's explicit rejection before an NCBI option error, both exit 1), TX-6/BP-9/TN-8 (shared rejection texts without the program name), TX-8/IN-12 (frozen web ABI v1), TX-9 (`validate` without inputs), BP-10 (`LOSAT_STARTUP_TRACE`).

## Comparison discipline

Run NCBI and LOSAT with identical arguments from the same (scratch) cwd, with an environment without BATCH_SIZE, CHUNK_SIZE, OVERLAP_CHUNK_SIZE, BL2SEQ_LEGACY, CTOOLKIT_COMPATIBLE, PRE_FETCH_SEQS_LIMIT, OLD_FSC, ADAPTIVE_CBS, DIAG_*, NCBI_CONFIG_*, LOSAT_*, RAYON_* unless a case sets one on purpose, and no `~/.ncbirc`; compare stdout bytes, stderr bytes, exit status (and `2>&1` for outfmt 0 where warnings interleave). Hit counts alone are never evidence of parity. Record every command and its outcome in a TSV or log in your work dir as you go.

## Report

Write `REPORT.md` (or `FINDINGS.txt`) in your work dir as you go, and end your final reply with the full list. For each finding: an ID (`R2<angle>-n`), severity (high/medium/low), the NCBI file:line, the LOSAT file:line, a reproducing command (with input paths), observed NCBI vs LOSAT, and whether it is a defect, an accepted exception, a pending-the-maintainer item, or a justified rejection. Also list every round 1 finding of your angle with its re-run verdict (SAME / LOSAT-REJECTS as decided / ACCEPTED / DIFF). End with an overall verdict for your angle: supported / unsupported / inconclusive, and why.
