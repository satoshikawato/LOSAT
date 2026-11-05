# Approvals of the conditional switches added after the stop rule

The stop rule for conditional switches ([`AGENTS.md`](../../AGENTS.md), "Speed-up switches", Owner,
2026-10-27) admits a new path that runs only under a rule "only as a deadline fallback that its body
names as such and the Owner approves". Twenty switches entered `STRICT_SET` after the rule reached
`dev` with PR #556 (merged 2026-10-27 05:00 UTC) and have the class `conditional` in
[`LOSAT/tests/switch_reach_expected.tsv`](../../LOSAT/tests/switch_reach_expected.tsv). As of dev
8a5656f7f (2026-11-05): 1 approved by the Owner as a fallback, 14 inside an existing rule, 5 without
an Owner approval.

- (i) The Owner approved it as a deadline fallback. The Owner's words and date are given.
- (ii) Not a new rule: the switch's path runs wherever an existing path of the code runs. Its TSV row
  is `conditional` because that existing path has a condition, or because a cap or a `min` leaves
  small inputs unchanged. Constants that size the new form (block bytes, look-ahead) choose no path.
- (iii) A new rule (a new threshold, capacity or thread split chooses the path), and no Owner approval
  was found. The pull requests recorded only Owner-delegated decisions, which are not the Owner's
  approval for this rule (Owner, 2026-11-05).
- Merge times are of the pull request on `dev` that added the name to `STRICT_SET`
  (`LOSAT/src/utils/xswitch.rs`). Four merged on 2026-10-27 after #556; #548 was opened before it.
- Two more switches entered after the rule and are gone: `LOSAT_X_BPQSLOTS` (#579) and
  `LOSAT_X_BNQSLOTS` (#580), removed with the query-batch slots on 2026-10-28.

| Switch (`LOSAT_X_`) | Program | Pull request (merged, UTC) | Where its path runs | Class | Evidence |
| --- | --- | --- | --- | --- | --- |
| `TBXFEWGROUPS` | tblastx | #652 (2026-10-31 13:03, reached `dev` with #657) | a unit-path batch that needs fewer than six context groups; new threshold `X_FEW_GROUPS_LETTERS = 1 << 17` (`tblastx/x_pair_par.rs:346`) | (i) | #652's body: "Rule-limited, so deadline fallback ... The Owner's OK is needed before ready". Owner, 2026-10-31: 「これも」, continuing 「直してマージして」 (Fix and merge this one too). No deadline date is named. |
| `BPDBKEEP` | blastp | #576 (2026-10-27 08:03) | every `-db` search (`blastp/blast_engine/database/x_db_keep.rs:460`); block size `X_KEEP_BLOCK_BYTES = 96 << 10` sizes the kept scratch | (ii) | No input rule; the `-subject` side is the other path. |
| `BXCHUNKCROSS` | blastx | #583 (2026-10-27 16:06) | the chunk loop of a split query under `LOSAT_X_BXCHUNKFLOW`'s in-order rule (`blastx/x_subject_flow.rs:252`); look-ahead `X_CROSS_AHEAD = 16` | (ii) | #583's body: the form holds "wherever the serial chunk loop searched two or more chunks with the subject flow (`LOSAT_X_BXCHUNKFLOW`'s in-order case ...)". |
| `BXTASKRES` | blastx | #593 (2026-10-28 04:42) | wherever the subject flow of `LOSAT_X_BXSUBJFLOW` runs (`blastx/x_subject_flow.rs:94`); tasks cut by residues, no cut-off | (ii) | #593's body: no condition and no threshold. |
| `DCSEEDLEAN` | blastn | #594 (2026-10-28 05:13) | dc-megablast's pooled table build and partitioned stages (`LOSAT_X_DCCSR`, `LOSAT_X_PAIRPARFOLD`; `blastn/blast_engine/x_pair_seed_disc.rs:969`) | (ii) | #594's body: "one form for every dc-megablast pool build and partitioned stage, with no new threshold and no thread rule". |
| `BNSEEDLISTS` | blastn | #606 (2026-10-28 13:17) | the partitioned stages of `LOSAT_X_PAIRPAR` and `LOSAT_X_PAIRPARNA` (`blastn/blast_engine/x_pair_seed_na.rs:660`) | (ii) | #606's body: "no new threshold and no thread rule"; a stage-size rule was left out because it would be a new conditional switch. |
| `BXKAPPAFIT` | blastx | #612 (2026-10-29 09:36) | the composition redo at a batch end; `min(threads, matches)` workers (`blastx/runtime.rs:328`) | (ii) | A formula without a cut-off; the row is `conditional` because `min` leaves batches with as many matches as threads unchanged. |
| `BNCSRPACK` | blastn | #617 (2026-10-29 14:17) | `LOSAT_X_BNCSR` tables whose chains list words newest first (`blastn/x_bn_csr.rs:147`) | (ii) | The pull request removed a size rule. |
| `MBHITPF` | blastn | #636 (2026-10-30 15:19) | `search_query_batch`'s offset-pair loop on x86_64, not `LOSAT_X_PAIRPAR`'s seed stage (`blastn/blast_engine/run.rs:3249`) | (ii) | No input rule; the only condition is the target. |
| `BNDBWIDE` | blastn | #651 (2026-10-31 03:04) | `-db` at 2+ threads: widens the block window of `LOSAT_X_BNDBSTREAM` (`blastn/blast_engine/x_db_stream.rs:281`); `X_BNDBWIDE_BLOCKS_PER_THREAD = 16` | (ii) | A database of fewer blocks than the old window is unchanged. Same form as `BPDBWIDE`. |
| `TBXBOARD` | tblastx | #657 (2026-10-31 13:03) | every batch of a native search at 2+ threads, where a pool exists (`tblastx/blast_engine/run_impl.rs:4660`) | (ii) | #657's body: "No size or count rule ... the one rule '2+ threads'". |
| `BPDBWIDE` | blastp | #670 (2026-10-31 16:06) | `-db` at 2+ threads: widens the existing block window (`blastp/blast_engine/database/x_db_keep.rs:494`); `X_BPDBWIDE_BLOCKS_PER_THREAD = 64` | (ii) | Owner, 2026-11-01: 「規則で動く switch ではない」 (Not a switch that runs under a rule). |
| `TBNPLANCUT` | tblastn | #705 (2026-11-01 20:59) | the chunk flow of `LOSAT_X_TBNOIDNEXT` with `-db`, 2+ threads and several subjects (`common/x_db_stream.rs:381`) | (ii) | #705's body: no rule and no threshold; the block keeps its 256 KiB size. |
| `TBNLOOPFREE` | tblastn | #723 (2026-11-02 10:36) | wherever the subject flow of `LOSAT_X_TBNOIDNEXT` runs (`tblastn/stage_d_pipeline.rs:2998`) | (ii) | #723's body: no rule and no threshold. |
| `BNDBSMALL` | blastn | #743 (2026-11-03 00:31) | `-db` at 2+ threads: caps the largest block of the existing cut at `X_BNDBSMALL_BLOCK_BYTES = 256 << 10` (`blastn/blast_engine/x_db_cut.rs:151`) | (ii) | #743's body: "There is no tuning threshold beyond the block size"; a database whose blocks are already smaller is unchanged. |
| `BXPOOLN1` | blastx | #548 (2026-10-27 05:10; opened 2026-10-26) | word-size-3 tables at one thread are built without the pool; new split `threads <= 1` (`tblastx/lookup/x_lut_direct.rs:157`) | (iii) | Only Owner-delegated decisions. Opened before the rule and merged 10 minutes after it reached `dev`. |
| `BPLK16` | blastp | #591 (2026-10-27 16:32) | exact blastp with a window whose block lookup query is at most 65,535 letters (`blastp/blast_engine.rs:7191`) | (iii) | New size condition (the width of `u16`, not a tuned cut-off). Only Owner-delegated decisions. |
| `MBCSR` | blastn | #636 (2026-10-30 15:19) | megablast's 4^12-cell table when the subjects' letters are at least 8 times the word-loop positions plus the cells; `X_MB12_LETTERS_PER_POSITION = 8` (`blastn/x_bn_csr.rs:194`) | (iii) | New cost-model threshold. Only Owner-delegated decisions. |
| `BXCLUTCAP` | blastx | #648 (2026-10-30 23:53) | compressed lookup tables where `LOSAT_X_BXUNITCAP`'s rule is asked; new capacity `CLUT_UNIT_CAP = 3` and a few-subjects clause (`blastx/x_bx_mem.rs:154`) | (iii) | New capacity and subject-count rule. Only Owner-delegated decisions. |
| `BNPVONCE` | blastn | #650 (2026-10-31 03:03) | packed CSR heads whose bit array is at most `X_PV_ONCE_MAX_BITS = 1 << 20` bytes (`blastn/blast_engine/run.rs:3267`) | (iii) | New cache-size threshold. Only Owner-delegated decisions. |

For each switch in (iii), the Owner chooses one of two: approve it as a deadline fallback (with a
date), or fold it into the default form (structural line SF, stage S3). Until then each is a candidate
for S3. `TBXFEWGROUPS` (i) has no deadline date yet.
