#!/usr/bin/env python3
"""Compare BLASTN outfmt 0 on subject deflines of punctuation with NCBI (Session S07+; comparison only).

The eleventh S07+ audit found that NCBI's `x_CleanAndCompress` reads past the end of some
subject titles made of `,;~` and spaces (NCBI crashes in outfmt 0), and LOSAT rejects those
titles explicitly (AUTHORITY.md §N). This sweep runs every defline of one to five of these
characters that does not start with a space (1023 deflines) as the subject's title with
-outfmt 0 in NCBI and LOSAT. Each defline is "same" (the same stdout, stderr and exit
status) or "crash-rejected" (NCBI dies of a signal and LOSAT fails with "which LOSAT does
not reproduce"); any other outcome is printed. Exits 1 when any is.

Usage: title_sweep.py --bin-dir DIR --losat LOSAT --work DIR [--jobs N]
"""
from __future__ import annotations

import argparse
import itertools
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "losat_web_e2a"))
from run_oracle import ENGINE  # noqa: E402

QUERY = ENGINE / "tests/fasta/outfmt0/multi_query.fasta"
SUBJECT = ENGINE / "tests/fasta/outfmt0/multi_subject.fasta"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--losat", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    residues = "".join(SUBJECT.read_text().split(">")[1].splitlines()[1:])
    titles = ["".join(chars) for length in range(1, 6) for chars in itertools.product(",;~ ", repeat=length)
              if chars[0] != " "]

    def run(case: tuple[int, str]) -> tuple[str, str, tuple[int, int] | None]:
        index, title = case
        subject = args.work / f"title_{index}.fa"
        subject.write_text(f">{title}\n{residues}\n")
        argv = ["-query", str(QUERY), "-subject", str(subject)]
        try:
            ncbi = subprocess.run([str(args.bin_dir / "blastn"), *argv], capture_output=True, timeout=60)
            ours = subprocess.run([str(args.losat.resolve()), "blastn", *argv], capture_output=True, timeout=60)
        except subprocess.TimeoutExpired:
            return title, "timeout", None
        codes = (ncbi.returncode, ours.returncode)
        if (ncbi.returncode, ncbi.stdout, ncbi.stderr) == (ours.returncode, ours.stdout, ours.stderr):
            return title, "same", codes
        if ncbi.returncode < 0 and ours.returncode != 0 and b"which LOSAT does not reproduce" in ours.stderr:
            return title, "crash-rejected", codes
        return title, "different", codes

    counts: dict[str, int] = {}
    with ThreadPoolExecutor(args.jobs) as pool:
        for title, outcome, codes in pool.map(run, enumerate(titles)):
            counts[outcome] = counts.get(outcome, 0) + 1
            if outcome not in ("same", "crash-rejected"):
                print(f"{outcome.upper()}\t{title!r}\t{codes}")
    unexpected = len(titles) - counts.get("same", 0) - counts.get("crash-rejected", 0)
    print(f"# deflines={len(titles)} same={counts.get('same', 0)} "
          f"crash-rejected={counts.get('crash-rejected', 0)} unexpected={unexpected}")
    return 1 if unexpected else 0


if __name__ == "__main__":
    sys.exit(main())
