# Inputs of the makeblastdb version 5 cases

`b2.fa`: three lines (`>x`, `>y`, `MK`) that `src/writedb/tests/makeblastdb_v5_cases.rs` case `M_ignored` gives to
`makeblastdb` with `a.fa`: an empty record `x`, then `y`; origin not recorded (added in commit 5f19ec192, 2026-10-29).
