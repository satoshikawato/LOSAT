# Independent C-NA64 native TBLASTX audit

The read-only `ncbi_parity_auditor` supports raw parity, time nonregression and
RSS PASS for these six fresh B1 versus C-NA64 conditions. This is recovery to the
permitted range, not evidence of a general 5% speedup. Original integrated native
failures remain unchanged.

| Input | Threads | B1 median (s) | C-NA64 median (s) | Reduction |
|---|---:|---:|---:|---:|
| Mela/Pemo | 1 | 3.629585534 | 3.572273562 | 1.5790% |
| Mela/Pemo | 8 | 2.212022835 | 2.155660988 | 2.5480% |
| Mje/Mela | 1 | 16.838524569 | 16.833425263 | 0.0303% |
| Mje/Mela | 8 | 10.206712219 | 10.087494794 | 1.1680% |
| APself | 1 | 30.075042401 | 29.738963546 | 1.1175% |
| APself | 8 | 18.101631972 | 18.015462318 | 0.4760% |

All 87 records (3 official oracles, 12 diagnostics, 12 warmups, 60 timed runs)
have actual raw output equal to the corresponding fresh official oracle:
151,222,240 bytes checked. Commands, results and usage agree for all records.
Normalized arguments equal the previous integrated group apart from executable
and output paths. There are five measurements per version per condition, zero
exclusions, the prescribed alternating order, and no overlapping monotonic
intervals (minimum gap 42.785316 ms). All six B/C sample ranges overlap.

Time uses max(5% of baseline median, 50 ms); RSS uses max(10% of baseline maximum,
16 MiB). The auditor corrected an initial 5% RSS transcription/calculation and
recomputed all six conditions: every allowance remains 16,777,216 bytes and all
PASS decisions remain identical. AP n8 baseline maximum RSS depends on one
143,699,968-byte sample; no persistent memory reduction is claimed.

Twenty-four records (17 timed) have realtime clock disagreement. Maximum
realtime-minus-monotonic is +1.104738499 seconds; maximum absolute boottime versus
monotonic difference is 8.859992 microseconds. All samples remain present without
correction or exclusion; monotonic agreement passes all 87.

The auditor verified all 589 paths/hashes in the root original source, original
snapshot and new C-NA64 snapshot; only Cargo configuration differs. Actual B1
artifact `ae494b43...` and C-NA64 `a1ab2845...`, four fixture files, all tools and
runners, and production metadata match their records. C-NA64 is byte-identical
to the separately built alignment diagnostic. See the separate source-cwd
metadata correction and profile-binding note; mutable present-day library
fingerprints do not prove historical compiler invocation identity for B1.

These are gencode-1, outfmt-6 native searches, outside the local-subject genetic
code exception. Native does not invoke Node. Actual stderr confirms n1 pool 0
and n8 pool 8; host worker observations are N/A. Other native conditions,
remaining final gates and memory acceptance require their own evidence.
