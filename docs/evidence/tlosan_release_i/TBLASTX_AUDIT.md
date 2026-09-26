# Session I exact-candidate TBLASTX supported-profile audit

Status: **PASS — all 20 unique supported-profile contracts**. This is the separate 20-case v0.1.0 supported-profile gate on the v0.2.0 candidate. It is neither Session H's historical 13/20 partial run nor Stage G's focused 12-case regression.

## Candidate and authority

- Source candidate: `005e3d4b6cba6b5808334088fe9595c89efe01f8`.
- Linux x64 Native executable SHA-256: `6c50c85d4f14379e3d71952c7cc7c9c3f95976e4bc09921ab7cc421390500b40`.
- NCBI `tblastx 2.17.0+` executable SHA-256: `583e5d60bbd444ac455d20e0956c5aa0aeef675da8daee8204d8f9376ddb8804`; [identity](tblastx_oracle_identity.json).
- NCBI C/C++ source commit: `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`.
- The [20-case manifest](../../../LOSAT/tests/tblastx_v010_parity_manifest.tsv), unchanged [runner](../../../LOSAT/tests/audit_tblastx_v010.py), [frozen Native hashes](../../../LOSAT/tests/platform_native_v010_canonical.tsv), and [v0.1.0 profile authority](../../release/tblastx_v0.1.0_certification.md) govern acceptance.

The profile is local query/subject TBLASTX, outfmt 6, one thread, genetic-code options selected by each manifest row and all remaining search defaults. Native multithread equivalence, threaded TBLASTX WASI, arbitrary options, database/remote searches and performance are outside this audit.

## Method and retained evidence

Every raw output is reclassified by [summarize_tblastx_20.py](summarize_tblastx_20.py). It checks exact commands against the committed manifest, unchanged fixture inputs against the candidate commit, pinned binary/oracle hashes, empty stderr, output IDs, frozen LOSAT raw hashes, every required three-run hash, and the implicit/explicit code-1 control. Acceptance requires exactly 20 unique cases: 14 `EXACT_TEXT` and only the six designated subject-code rows as `HSP_SET_DIFF`. No sorting or hit-count tolerance can accept a parity row.

Independent jobs shortened the oracle audit without changing per-case commands. The selected raw sources are:

| Source partition | Cases retained | Evidence |
| --- | --- | --- |
| Original full runner prefix | p01–p06 | [prefix summary](tblastx_prefix_6.json), [partial log](tblastx_prefix_partial.log) |
| Later parity runner prefix | p07–p11 | [manifest](later_parity.tsv), [completed-prefix log](later_prefix.log); p07–p11 retained |
| Selected-code runner prefix | d01 | [manifest](deviation.tsv), [completed-prefix log](deviation_prefix.log); redundant d02 job stopped after d01 PASS |
| Four independent selected-code jobs | d02–d05 | [records](direct_d02_d05.json), [log](direct_d02_d05.log), [replay runner](run_tblastx_d02_d05.py) |
| Future controls runner | p12–p14, d06 | [manifest](future_controls.tsv), [classifications](future_classifications.json), [log](future_controls.log), [code-1 control](default_code_equivalence.json) |

Raw files remain in the recorded `/tmp/tlosan-session-i-*005e3d4b*` directories. The final aggregate records every command, fixture hash, output hash, row count, stderr hash, repeat hash and acceptance predicate. Interrupted duplicate cases are never counted; the aggregate selects one completed result per committed case ID. The direct replay runner preserves the exact command builder and predicates used by the original temporary four-job script.

## Final result

[Strict aggregate evidence](tblastx_20_summary.json) returns `PASS`: 20 unique cases, **14 `EXACT_TEXT`**, **6 approved `HSP_SET_DIFF`**. Every LOSAT raw hash equals its frozen certified canonical hash. All 14 parity cases and d06 are repeatable across three native runs; all oracle/candidate/repeat stderr files are empty. The p12 implicit/explicit code-1 NCBI and LOSAT outputs all hash `86c05a04efb50e4026720e2d44fe2db2e6446f9594e174f3fde56931d09d5b49`. p14 query-code-4/default-subject stays exact; d06 alone changes the selected subject code, preserving its prior certified hash and 12,672/13,644 NCBI/LOSAT rows.

The final p11 output is exact at 9,346 rows; NCBI and all three LOSAT outputs share SHA-256 `1eb11f5caa4d1030a016e67bafc90acaecb0a1fe76f9b11256fe1343f0571fe4`. No unexpected `ORDER_ONLY`, `VALUE_DIFF`, missing-HSP, execution-error or other non-exact parity case remains. The [independent audit](INDEPENDENT_AUDIT.md) reviewed raw files, fixed identities, aggregate predicates, exceptions and final scope. No new NCBI discrepancy required a source correction.

## Narrow subject-code exception

NCBI accepts `db_gencode` through `c++/src/algo/blast/blastinput/blast_args.cpp:1043–1055`, but the local-subject search source is created by `c++/src/algo/blast/api/local_db_adapter.cpp:124–139` → `c++/src/algo/blast/api/seqsrc_query_factory.cpp:105–151` → `c++/src/algo/blast/api/blast_objmgr_tools.cpp:414–420` (`CBlastQuerySourceOM subj_src(subjects, prog)` without passing the DbGeneticCode option) → `c++/src/algo/blast/api/blast_setup_cxx.cpp:799–810`, which retrieves the local subject's per-sequence genetic-code ID. `c++/src/algo/blast/api/blast_objmgr_tools.cpp:145–169` supplies that local/default code; `c++/src/app/blast/tblastx_app.cpp:156` separately supplies the selected option to reporting. LOSAT `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:688,1397–1411` uses the selected subject code in search-time frame construction.

Only the approved TBLASTX local-subject nondefault genetic-code behavior is accepted. The default-subject query-code-4 p14 control must remain exact to NCBI; the subject-only d06 probe must retain its certified selected-code output and repeatability. All six designated selected-code LOSAT outputs must match their frozen certified raw hashes. No timing, ordering, scoring, filtering, statistics, pruning, coordinates or formatting deviation is added.

No new unexpected NCBI difference was observed in the completed 20-case gate. No Rust correction or benchmark was performed in this session.

## Reproduce

To repeat the strict aggregate against the retained raw directories:

```bash
python3 docs/evidence/tlosan_release_i/summarize_tblastx_20.py \
  --prefix /tmp/tlosan-session-i-tblastx-005e3d4b-20260926 \
  --later /tmp/tlosan-session-i-tblastx-later-005e3d4b \
  --future /tmp/tlosan-session-i-tblastx-future-005e3d4b \
  --deviation /tmp/tlosan-session-i-tblastx-deviation-005e3d4b \
  --direct /tmp/tlosan-session-i-d02-d05-005e3d4b \
  --output /tmp/tlosan-v020-tblastx-summary-replay.json
```

The simplest independent replay uses the unchanged full runner in a new directory, with the pinned oracle and candidate binary present:

```bash
python3 LOSAT/tests/audit_tblastx_v010.py --output-dir /tmp/tlosan-v020-tblastx-full-replay-new
python3 docs/evidence/tlosan_release_i/run_tblastx_serial_wasi.py --output-dir /tmp/tlosan-v020-tblastx-serial-replay-new
```

The full runner emits the 20 classifications, exact commands, oracle identity, required repeats and code-1 control. It does not measure a formal benchmark. Compare replay outputs with the frozen canonical hashes as the Session I aggregate does. The separate serial command-WASI p03/p12/d06 replay passed 3/3 with the exact candidate module; its [record](tblastx_serial_wasi.json) includes host, command, input and output hashes.
