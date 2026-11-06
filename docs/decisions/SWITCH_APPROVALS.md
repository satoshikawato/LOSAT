# Approvals of the conditional switches added after the stop rule

The stop rule for conditional switches ([`AGENTS.md`](../../AGENTS.md), "Speed-up switches", Owner,
2026-10-27) lets a rule-limited path land "only as a deadline fallback that its body names as such
and the Owner approves". Twenty switches entered `STRICT_SET` after the rule reached `dev` with PR
#556 (merged 2026-10-27 05:00 UTC) and have the class `conditional` in
[`LOSAT/tests/switch_reach_expected.tsv`](../../LOSAT/tests/switch_reach_expected.tsv). As of dev
17bf5086a (2026-11-05): 1 approved by the Owner as a deadline fallback, whose deadline is its fold
into the default form in stage S3 (Owner, 2026-11-06), 14 inside an existing rule, 5 without an
Owner approval, which the Owner chose to fold into the default form (Owner, 2026-11-05).

- (i) The Owner approved it as a deadline fallback. The Owner's words and date are given, with how
  they were read.
- (ii) Not a new rule: the switch's path runs wherever an existing path of the code runs. Its TSV row
  is `conditional` because that existing path has a condition, or because a cap or a `min` leaves
  small inputs unchanged. Constants that size the new form (block bytes, look-ahead) choose no path.
- (iii) A new rule (a new threshold, capacity or thread split chooses the path), and no Owner approval
  of it was found. The pull requests record Owner-delegated decisions, which are not that approval
  (Owner, 2026-11-05). The Owner chose to fold these into the default form in structural line SF,
  stage S3, not to approve them as deadline fallbacks (Owner, 2026-11-05).
- Paths in the table are under `LOSAT/src/algorithm/`. The line is the switch's `x_reach!` mark or
  the named constant.
- The Owner's words quoted here have rows, with their sources, in
  [`OWNER_DECISIONS.md`](OWNER_DECISIONS.md).
- Merge times are of the pull request on `dev` that added the name to `STRICT_SET`
  (`LOSAT/src/utils/xswitch.rs`). Four merged on 2026-10-27 after #556; #548 was opened before it.
- Two more switches entered after the rule and are gone: `LOSAT_X_BPQSLOTS` (#579), removed by
  #609 on 2026-10-28, and `LOSAT_X_BNQSLOTS` (#580), removed with blastn's query-batch slots
  (`LOSAT_X_QBATCHPAR`) by #610 on 2026-10-29.

| Switch (`LOSAT_X_`) | Program | Pull request (merged, UTC) | Where its path runs | Class | Evidence |
| --- | --- | --- | --- | --- | --- |
| `TBXFEWGROUPS` | tblastx | #652 (2026-10-31 13:56; its commits reached `dev` with #657 at 13:03) | a unit-path batch that needs fewer than six context groups; new threshold `X_FEW_GROUPS_LETTERS = 1 << 17` (`tblastx/x_pair_par.rs:346`) | (i) | #652's body: "Rule-limited, so deadline fallback ... The Owner's OK is needed before ready". Owner, 2026-10-31: 「これも」, continuing 「直してマージして」 (Fix and merge this one too); the merging session read this as the OK. The Owner confirmed the reading (Owner, 2026-11-06: 「(a)【推奨】」); the deadline is its fold into the default form in S3, not a date. |
| `BPDBKEEP` | blastp | #576 (2026-10-27 08:03) | every `-db` search (`blastp/blast_engine/database/x_db_keep.rs:460`); block size `X_KEEP_BLOCK_BYTES = 96 << 10` sizes the kept scratch | (ii) | No input rule; the `-subject` side is the other path. |
| `BXCHUNKCROSS` | blastx | #583 (2026-10-27 16:06) | the chunk loop of a split query under `LOSAT_X_BXCHUNKFLOW`'s in-order rule (`blastx/x_subject_flow.rs:252`); look-ahead `X_CROSS_AHEAD = 16` | (ii) | #583's body: the form holds "wherever the serial chunk loop searched two or more chunks with the subject flow (`LOSAT_X_BXCHUNKFLOW`'s in-order case ...)". |
| `BXTASKRES` | blastx | #593 (2026-10-28 04:42) | wherever the subject flow of `LOSAT_X_BXSUBJFLOW` runs (`blastx/x_subject_flow.rs:94`); tasks cut by residues, no cut-off | (ii) | #593's body: no condition and no threshold. |
| `DCSEEDLEAN` | blastn | #594 (2026-10-28 05:13) | dc-megablast's pooled table build and partitioned stages (`LOSAT_X_DCCSR`, `LOSAT_X_PAIRPARFOLD`; `blastn/blast_engine/x_pair_seed_disc.rs:969`) | (ii) | #594's body: "one form for every dc-megablast pool build and partitioned stage, with no new threshold and no thread rule". |
| `BNSEEDLISTS` | blastn | #606 (2026-10-28 13:17) | the partitioned stages of `LOSAT_X_PAIRPAR` and `LOSAT_X_PAIRPARNA` (`blastn/blast_engine/x_pair_seed_na.rs:660`) | (ii) | #606's body: "no new threshold and no thread rule"; a stage-size rule was left out because it would be a new conditional switch. |
| `BXKAPPAFIT` | blastx | #612 (2026-10-29 09:36) | the composition redo at a batch end; `min(threads, matches)` workers (`blastx/runtime.rs:328`) | (ii) | A formula without a cut-off; the row is `conditional` because `min` leaves batches with as many matches as threads unchanged. |
| `BNCSRPACK` | blastn | #617 (2026-10-29 14:17) | `LOSAT_X_BNCSR` tables whose chains list words newest first (`blastn/x_bn_csr.rs:147`) | (ii) | The pull request removed a size rule. |
| `MBHITPF` | blastn | #636 (2026-10-30 15:19) | `search_query_batch`'s offset-pair loop on x86_64, not `LOSAT_X_PAIRPAR`'s seed stage (`blastn/blast_engine/run.rs:3249`) | (ii) | No input rule; the only condition is the target. |
| `BNDBWIDE` | blastn | #651 (2026-10-31 03:04) | `-db` at 2+ threads: widens the block window of `LOSAT_X_BNDBSTREAM` (`blastn/blast_engine/x_db_stream.rs:281`) to `X_BNDBWIDE_BLOCKS_PER_THREAD = 16` (:109) while at most `X_BNDBWIDE_WAITING_BYTES_PER_THREAD = 1 << 20` bytes per thread of results wait (:124) | (ii) | A database of fewer blocks than the old window is unchanged. The same form as `BPDBWIDE`, which the Owner judged not to run under a rule (2026-11-01, about #670); the Owner confirmed on 2026-11-05 that the answer covers #651 too. |
| `TBXBOARD` | tblastx | #657 (2026-10-31 13:03) | every batch of a native search at 2+ threads, where a pool exists (`tblastx/blast_engine/run_impl.rs:4751`) | (ii) | #657's body: "No size or count rule ... the one rule '2+ threads'". |
| `BPDBWIDE` | blastp | #670 (2026-10-31 16:06) | `-db` at 2+ threads: widens the existing block window (`blastp/blast_engine/database/x_db_keep.rs:475`) to `X_BPDBWIDE_BLOCKS_PER_THREAD = 64` (:89) while at most `X_BPDBWIDE_WAITING_BYTES_PER_THREAD = 1 << 20` bytes per thread of results wait (:94) | (ii) | Owner, 2026-11-01: 「規則で動く switch ではない」 (Not a switch that runs under a rule). |
| `TBNPLANCUT` | tblastn | #705 (2026-11-01 20:59) | the chunk flow of `LOSAT_X_TBNOIDNEXT` with `-db`, 2+ threads and several subjects (`common/x_db_stream.rs:381`) | (ii) | #705's body: no rule and no threshold; the block keeps its 256 KiB size. |
| `TBNLOOPFREE` | tblastn | #723 (2026-11-02 10:36) | wherever the subject flow of `LOSAT_X_TBNOIDNEXT` runs (`tblastn/stage_d_pipeline.rs:2998`) | (ii) | #723's body: no rule and no threshold. |
| `BNDBSMALL` | blastn | #743 (2026-11-03 00:31) | `-db` at 2+ threads: caps the largest block of the existing cut at `X_BNDBSMALL_BLOCK_BYTES = 256 << 10` (`blastn/blast_engine/x_db_cut.rs:151`) | (ii) | #743's body: "There is no tuning threshold beyond the block size"; a database whose blocks are already smaller is unchanged. |
| `BXPOOLN1` | blastx | #548 (2026-10-27 05:10; opened 2026-10-26) | word-size-3 tables at one thread are built without the pool; new split `threads <= 1` (`tblastx/lookup/x_lut_direct.rs:157`) | (iii) | No Owner approval of the rule; the body records Owner-delegated decisions. Opened before the rule and merged 10 minutes after it reached `dev`. To be folded in S3 (Owner, 2026-11-05). |
| `BPLK16` | blastp | #591 (2026-10-27 16:32) | exact blastp with a window whose block lookup query is at most 65,535 letters (`blastp/blast_engine.rs:7452`) | (iii) | New size condition (the width of `u16`, not a tuned cut-off). No Owner approval of the rule: the body cites an Owner decision for the goal (memory-floor rows, 2026-10-26), not for this condition. To be folded in S3 (Owner, 2026-11-05). |
| `MBCSR` | blastn | #636 (2026-10-30 15:19) | megablast's 4^12-cell table when the subjects' letters are at least 8 times the word-loop positions plus the cells; `X_MB12_LETTERS_PER_POSITION = 8` (`blastn/x_bn_csr.rs:194`) | (iii) | New cost-model threshold. No Owner approval of the rule; the body records Owner-delegated decisions. To be folded in S3 (Owner, 2026-11-05). |
| `BXCLUTCAP` | blastx | #648 (2026-10-30 23:53) | compressed lookup tables where `LOSAT_X_BXUNITCAP`'s rule is asked; new capacity `CLUT_UNIT_CAP = 3` and a few-subjects clause (`blastx/x_bx_mem.rs:154`) | (iii) | New capacity and subject-count rule. No Owner approval of the rule; the body records Owner-delegated decisions. To be folded in S3 (Owner, 2026-11-05). |
| `BNPVONCE` | blastn | #650 (2026-10-31 03:03) | packed CSR heads whose bit array is at most `X_PV_ONCE_MAX_BITS = 1 << 20` bytes (`blastn/blast_engine/run.rs:3267`) | (iii) | New cache-size threshold. No Owner approval of the rule; the body records Owner-delegated decisions. To be folded in S3 (Owner, 2026-11-05). |

For each switch in (iii), the Owner chose to fold it into the default form in stage S3 (one form for
every input, or removed with its rule), and confirmed that the 2026-11-01 answer about `BPDBWIDE`
covers `BNDBWIDE` (Owner, 2026-11-05: 「Q9：すべて推奨通り。」). The Owner's 「これも」 of 2026-10-31
approves `TBXFEWGROUPS` as a deadline fallback, and its deadline is its fold into the default form
in S3, with the five above (Owner, 2026-11-06: 「(a)【推奨】」).
