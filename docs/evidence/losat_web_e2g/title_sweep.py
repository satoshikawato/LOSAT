#!/usr/bin/env python3
"""Compare BLASTN outfmt 0 on subject deflines of punctuation with NCBI (Session S07+++b; comparison only).

docs/evidence/losat_web_e2c/title_sweep.py with the expectation of E2g: LOSAT no longer
rejects the titles whose cleanup NCBI's `x_CleanAndCompress` runs past the end of the
string (NCBI crashes in outfmt 0), but stops the cleanup at the end of the string
(approved exception 2 of PD-LOSAT-NCBI-DEFECTS, c72452236). Every defline of one to five
of `,;~` and space that does not start with a space (1023 deflines) is the title of a
subject with hits. Each is "same" (stdout, stderr and exit status) or "exception-2": NCBI
dies of a signal, LOSAT exits 0 with NCBI's stderr, and LOSAT's stdout equals NCBI's for
the same subject with a placeholder title of the length of LOSAT's title once the
placeholder is replaced by the title that NCBI's loop gives when it stops at the end of
the string (`title` below). Any other outcome is printed; exits 1 when any is.

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


def clean_and_compress(text: str) -> str:
    """create_defline.cpp:219-312, stopping at the end of the string."""
    t = text.strip(" ").encode()
    out = bytearray()
    if not t:
        return ""
    at = lambda i: t[i] if i < len(t) else 0  # noqa: E731
    i, left = 1, len(t) - 1
    curr = at(0)
    two = curr
    while left > 0:
        nxt = at(i)
        i += 1
        two = ((two << 8) | nxt) & 0xFFFF
        a, b = two >> 8, two & 0xFF
        if (a, b) == (44, 44):
            out.append(curr)
            nxt = 32
        elif (a, b) in ((32, 32), (32, 41)):
            pass
        elif (a, b) == (40, 32):
            nxt = curr
            two = curr
        elif (a, b) in ((32, 44), (32, 59)):
            out.append(nxt)
            nxt = curr
            two = curr
        elif a in (44, 59) and b == 32:
            sep = a
            out.append(curr)
            out.append(32)
            while nxt in (32, sep):
                nxt = at(i)
                i += 1
                left = max(left - 1, 0)
            two = nxt
        else:
            out.append(curr)
        curr = nxt
        left = max(left - 1, 0)
    if curr > 0 and curr != 32:
        out.append(curr)
    return out.decode()


def trim_end(text: str, chars: str) -> str:
    for index in range(len(text) - 1, -1, -1):
        if text[index] not in chars:
            return text[: index + 1]
    return text


def title(defline: str) -> str:
    """The heading and description title (no TPA prefix in these deflines)."""
    text = trim_end(defline, ".,;~ ").lstrip(" ")
    return clean_and_compress(trim_end(text, ",;~ "))


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
        index, defline = case
        subject = args.work / f"title_{index}.fa"
        subject.write_text(f">{defline}\n{residues}\n")
        argv = ["-query", str(QUERY), "-subject", str(subject)]
        try:
            ncbi = subprocess.run([str(args.bin_dir / "blastn"), *argv], capture_output=True, timeout=60)
            ours = subprocess.run([str(args.losat.resolve()), "blastn", *argv], capture_output=True, timeout=60)
        except subprocess.TimeoutExpired:
            return defline, "timeout", None
        codes = (ncbi.returncode, ours.returncode)
        if (ncbi.returncode, ncbi.stdout, ncbi.stderr) == (ours.returncode, ours.stdout, ours.stderr):
            return defline, "same", codes
        expected_title = title(defline)
        if ncbi.returncode < 0 and ours.returncode == 0 and expected_title:
            placeholder = "Z" * len(expected_title)
            stand_in = args.work / f"title_{index}_p.fa"
            stand_in.write_text(f">{placeholder}\n{residues}\n")
            other = subprocess.run([str(args.bin_dir / "blastn"), "-query", str(QUERY), "-subject", str(stand_in)],
                                   capture_output=True, timeout=60)
            expected = (other.stdout.replace(placeholder.encode(), expected_title.encode())
                        .replace(stand_in.name.encode(), subject.name.encode()))
            if other.returncode == 0 and ours.stdout == expected and ours.stderr == other.stderr:
                return defline, "exception-2", codes
        return defline, "different", codes

    counts: dict[str, int] = {}
    with ThreadPoolExecutor(args.jobs) as pool:
        for defline, outcome, codes in pool.map(run, enumerate(titles)):
            counts[outcome] = counts.get(outcome, 0) + 1
            if outcome not in ("same", "exception-2"):
                print(f"{outcome.upper()}\t{defline!r}\t{codes}")
    unexpected = len(titles) - counts.get("same", 0) - counts.get("exception-2", 0)
    print(f"# deflines={len(titles)} same={counts.get('same', 0)} "
          f"exception-2={counts.get('exception-2', 0)} unexpected={unexpected}")
    return 1 if unexpected else 0


if __name__ == "__main__":
    sys.exit(main())
