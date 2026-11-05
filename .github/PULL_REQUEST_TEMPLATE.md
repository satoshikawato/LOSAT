<!--
Pull requests target `dev` (CONTRIBUTING.md, "Branches and auto-merge").
Open non-trivial changes as a draft after an issue was discussed with the maintainer.
Delete a section only when it does not apply, and say so in one line.
-->

## Why

<!-- The audit finding this answers (for example "SWEBOK 2026-11-04 KA4: get_hits"), the issue it closes, or the reason for the change. -->

## What changed

<!-- What the change does, file group by file group, and why this form. -->

## Behavior

<!--
Does any output change? If not, what shows it (frozen NCBI bytes, the tests' `=0` arms, the fast output regressions).
If it does, the NCBI C/C++ reference (file path and line numbers) and the fixture that shows the new output.
-->

## Speed-up pull requests

<!-- Required for a change that claims speed or memory (AGENTS.md, "Speed-up switches"). Otherwise write "Not a speed-up." -->

- Structural cause removed (what allocates, copies, waits or rescans, and why):
- Does the new form hold for every input?
- Inputs where it does not, with their ratios:

## Gate results

<!-- The tier that ran (quick, standard or full; skill `losat-gates`), the commit it ran on, and each check with its result. -->

## Changelog fragment

- [ ] `changelog.d/<branch name, with / written as ->.md` is added and `python3 scripts/changelog_fragments.py check` passes.
- [ ] No fragment, because: <!-- reason -->

## Owner-delegated decisions

<!-- Decisions made without asking the repository owner, under the owner's standing instruction (give its date). Write "None" if there are none. -->
