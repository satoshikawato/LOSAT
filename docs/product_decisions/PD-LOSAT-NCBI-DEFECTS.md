# Product Decision: NCBI BLAST+ behaviour that is a defect

- Decision ID: `PD-LOSAT-NCBI-DEFECTS`
- Version: 1.0
- Date: 2026-10-02
- Status: Accepted by the maintainer on 2026-10-02, in Session S07+++b (E2g), on the
  NCBI BLAST+ 2.17.0 behaviours that the BLASTN inventory
  (`docs/evidence/losat_web_e2g/INVENTORY.tsv`) and the independent audits found to be
  defects. Plan decision DW-15 in [`docs/losat_web_gui_plan.md`](../losat_web_gui_plan.md).

## Scope

NCBI BLAST+ 2.17.0 behaviour that is a defect rather than a design: a crash, an
exception that only a debug-build assertion was meant to prevent, a read past the end of
a buffer, or an integer that wraps. Everything else stays under the root
[`AGENTS.md`](../../AGENTS.md) bit-perfect rule. This version lists the BLASTN items; the
other programs add theirs in their sessions under the same rule.

## Rule

1. **NCBI fails, and LOSAT can give a valid result:** an approved exception. The result
   is valid when it equals NCBI's output for a nearby input or configuration that does
   not reach the defect (evidence under `docs/evidence/`).
2. **NCBI gives a deterministic result, even one that looks wrong:** LOSAT reproduces it
   byte for byte. This is not an exception.
3. **NCBI fails, and no valid result can be defined or checked:** an explicit rejection
   ("... is not supported by LOSAT's BLASTN"). This is not an exception.
4. **NCBI accepts an input deterministically that LOSAT rejected:** LOSAT ports it.

## Approved exceptions (BLASTN)

1. **A query chunk that NCBI would split again.** With `CHUNK_SIZE` and
   `OVERLAP_CHUNK_SIZE` such that a query chunk is long enough to be split again (an
   overlap close to the chunk size), NCBI skips the chunk's lookup table and stops with a
   `CCoreException` (null pointer) that names its build's source files, exit 3
   (`split_query_aux_priv.cpp:190-201`, `blast_aux_priv.cpp:206-208`; only a debug build
   asserts against it). LOSAT searches each chunk once. Evidence
   `docs/evidence/losat_web_e2g/resplit/`: 42 configurations (3 queries of 30 to 60 kb,
   blastn and megablast, 7 chunk and overlap pairs); NCBI exit 3 in all; LOSAT's output
   equals NCBI's unsplit output in 35, and NCBI's output at the largest overlap that does
   not split a chunk again in 74 of 84 runs (outfmt 6 and 0). The 10 others differ only
   by the HSPs at chunk boundaries that NCBI's own splits also lose. Fixtures
   `env.resplit_*` (expected output from NCBI at that largest overlap, `oracle_env`).
   With a negative `CHUNK_SIZE` above a negative `OVERLAP_CHUNK_SIZE` the same failure is
   an explicit rejection (below).
2. **outfmt 0 titles made only of punctuation.** For a subject defline of commas,
   semicolons, tildes and spaces that ends in a run of spaces and separators (such as
   `, ,`), NCBI's `x_CleanAndCompress` (`create_defline.cpp:219-312`) lets its `size_t`
   count of the remaining letters wrap, reads past the end of the string and crashes when
   it writes the title of such a subject with hits (SIGSEGV). LOSAT stops the cleanup at
   the end of the string (`, ,` gives `, `). Evidence
   `docs/evidence/losat_web_e2g/punct_defline/`: 8 such deflines on subjects with hits,
   megablast and blastn, with and without `-max_target_seqs 3`; NCBI crashes in all 4
   runs; LOSAT's report equals NCBI's report for the same subjects with placeholder
   deflines once each placeholder is replaced by LOSAT's title. Subjects with such
   deflines and no hits, where NCBI runs, match NCBI (fixtures `punct.nohit_*`).

## Reproduced as NCBI (rule 2)

- `CHUNK_SIZE=1000` without `BATCH_SIZE`: the batch size becomes 0, NCBI prints the
  outfmt 0 prolog and then `BLAST engine error: Empty CBlastQueryVector`, exit 3.
- `-max_target_seqs` from 2^30 to 2^31-51: the preliminary hit list size
  (`MIN(MAX(2*h,10),h+50)` in `Int4`) wraps and becomes 10.
- Defects without an effect on the output that LOSAT follows: the negative `Int4`
  offset of `BlastGetStartForGappedAlignmentNucl`, the out-of-range read of
  `hard_ranges[1]` for the last subject chunk, the out-of-range `double` to `Int4`
  conversion (INT_MIN), and the catch in `BlastSetupPreliminarySearchEx` that never
  matches.

## Explicit rejections (rule 3)

- A scoring system without Karlin-Altschul values when an invalid query (such as all N)
  comes first in its batch: NCBI goes on without a gapped block and crashes (SIGSEGV).
- `-max_target_seqs` above 2^31-51: the preliminary hit list size becomes negative and
  NCBI crashes.
- A reward of 32767 or a penalty of -32768 (after NCBI's 16-bit conversion):
  `BLAST_SCORE_MAX`/`BLAST_SCORE_MIN` fall outside NCBI's score range; a reward of 32767
  is counted past the end of the frequency array (invalid queries or a crash), a penalty
  of -32768 makes every query invalid.
- Subjects of 2^31 letters or more in total: NCBI's 32-bit total wraps.
- megablast gap costs above 32767: NCBI's 32-bit greedy distance arithmetic can wrap.
- A negative `CHUNK_SIZE` above a negative `OVERLAP_CHUNK_SIZE` that splits a query batch
  (`SplitQuery_CalculateNumChunks` gives more than one chunk): NCBI's `size_t` chunk ranges
  wrap; NCBI stops as in exception 1 where a chunk would be split again, and otherwise
  searches chunk ranges with gaps between them (hits are lost). Where such a pair does not
  split the batch (for example `-1` above `-2147483648`), NCBI searches as without the
  variables and so does LOSAT (fixtures `env.negative_pair_*`).

## Ported (rule 4)

- `-evalue` with a signed infinity or NaN (`+inf`, `-nan`, `+nan(1)`, `1e999`): NCBI's
  `CArg_Double` reads them with `strtod` and its option check (`<= 0`) lets them pass;
  NCBI searches as with the largest e-value (the cutoff stays 1, no HSP is reaped). LOSAT
  rejected them before Session S07+++b.
