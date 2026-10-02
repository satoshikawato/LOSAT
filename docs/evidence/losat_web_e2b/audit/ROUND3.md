# S08 (E2b) independent audit, round 3 (final binary)

Read COMMON.md first (same rules). Your work dir is `/home/kawato/.cache/losat-web-gui-target/s08-audit/r3<angle>/` (create it). If the Write tool refuses `.md`, write `FINDINGS.txt` and give the full findings in your final message.

FINAL_BINARY = `/home/kawato/.cache/losat-web-gui-target/s08-gate-native/release/LOSAT` and the reactors `/home/kawato/.cache/losat-web-gui-target/s08-gate-reactors/losat-web-{serial,threads}.wasm` are now built from commit `1117e8c17` (record their sha256 at the start and the end). The commits after round 2's binary (`7fbbfad96`), read them with `git -C /mnt/c/Users/genom/GitHub/LOSAT-web-gui show <sha>`:
- `32fd67a52` comment only;
- `aa4b6f8c3` BLAST_Cutoffs' floor of 1 in TBLASTX (`ncbi_cutoffs.rs` `blast_cutoffs_from_one`), PRE_FETCH_SEQS_LIMIT that NCBI cannot convert rejected;
- `6017f0ea2` adapter `register` checks in the CLI's order (query deflines first, subject deflines after residues); abi_v2.md wording;
- `5fa53b3f0` a subject that bio cannot read is deferred to the search start only when every line before its first defline is white space or a `!`/`#`/`;` comment (or its defline is not UTF-8); other text before the first defline is rejected at read;
- `1117e8c17` the closed-standard-output check of 7fbbfad96 is removed (it failed callers with a /dev/null opened read+write); a stdout closed at the start is now written to Rust's /dev/null and succeeds. This is put to the maintainer as a proposed CLI difference (not an approved exception yet): classify it as "pending the maintainer", not as a defect, unless the behaviour differs from that description.

Task: (1) re-run the reproductions of every finding of your previous rounds that was not yet verified on a binary with its fix (round 2's open items above) and classify each now (fixed / justified rejection / approved exception / recorded deferral / pending the maintainer / defect); (2) re-run a regression sample of your previous comparison sets on FINAL_BINARY (at least a quarter of them, all of the sets that touch the changed code) and report anything equal before and different now; (3) audit the five commits above against the NCBI source for your angle with a few new inputs. Keep it focused; this round is a confirmation. End with a verdict for your angle: supported / unsupported / inconclusive.
