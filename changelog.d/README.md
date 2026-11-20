# Pending CHANGELOG entries

A pull request does not edit `CHANGELOG.md`. It adds one file here:

- Name: the branch name with `/` written as `-`, then `.md`, for example
  `perf-tbx-unit-prefetch.md` for `perf/tbx-unit-prefetch`.
- Content: the entry as it will read under `## [Unreleased]`, one line starting with `- `
  (lines that continue an entry start with two spaces; a file may hold several entries).
- A later change of the same work edits its own file.

Why: every line added under `## [Unreleased]` went to the same place, right below the heading,
so any two open pull requests conflicted there and the second one became `DIRTY` when the first
merged. Files with different names never conflict.

`python3 scripts/changelog_fragments.py check` (CI) checks the form of every file here. A release
runs `python3 scripts/changelog_fragments.py collect`, which moves the entries under
`## [Unreleased]` in `CHANGELOG.md`, newest first (by the commit that added each file), and
deletes the files (`RELEASE.md`).
