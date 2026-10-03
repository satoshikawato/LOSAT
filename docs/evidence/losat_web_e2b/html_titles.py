#!/usr/bin/env python3
"""Subject titles with `&` against NCBI's HtmlDecode, for TBLASTX, TBLASTN and BLASTN (comparison only).

LOSAT does not decode HTML character references in the outfmt 0 titles; it rejects a shown
subject whose title NCBI's `NStr::HtmlDecode` changes (`report/defline.rs`,
`ncbi_nucleotide_title_is_decoded`, a port of NCBI's scan and entity table). For each title
(after the id `id`), one subject with hits gets the title and each program runs in outfmt 0:

- LOSAT accepts the title: LOSAT's stdout, stderr and exit status must equal NCBI's;
- LOSAT rejects it ("HTML character reference"): NCBI is run again with every `&` of the
  title replaced by a backquote, and its report (backquotes mapped back) must differ from
  the first one, i.e. NCBI decoded something. The same report means a needless rejection;
- NCBI crashes (a title that its x_CleanAndCompress reads past): LOSAT must reject it
  (TBLASTX and TBLASTN) or equal NCBI's report of the title cut at its end (BLASTN,
  approved exception 2 of PD-LOSAT-NCBI-DEFECTS; this script only records it).

The titles: every name of NCBI's entity table with and without `;`, numeric references
(decimal, hexadecimal, empty, too long), NCBI's scan details (`&xi;`, `&X41;`, a `;` after
16 characters, a final `;` that GenerateDefline trims first, TPA prefixes) and 600 seeded
random strings of `&#;xXamp1 F`.

Usage: html_titles.py --bin-dir DIR --losat LOSAT --work DIR [--jobs N]
"""
from __future__ import annotations

import argparse
import concurrent.futures
import random
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
OUTFMT0 = REPO / "LOSAT/tests/fasta/outfmt0"
NCBISTR = Path("/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/corelib/ncbistr.cpp")
CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"


def entity_names() -> list[str]:
    text = NCBISTR.read_text(encoding="latin-1").replace("\r", "")
    table = text[text.index("const s_HtmlEntities[] = {"):text.index("{    0, 0 }")]
    return re.findall(r'\{\s*\d+,\s*"([^"]*)"\s*\}', table)


def titles() -> list[str]:
    out = []
    for name in entity_names():
        out += [f"x &{name}; y", f"x &{name} y"]
    out += ["a&#38;b", "a&#x26;b", "a&#X26;b", "a&#0;b", "a&#;b", "a&#x;b", "a&#xg;b", "a&#65", "a&#65 ;",
            "a&#12345678901234567;b", "a&#1234567890123456;b", "s &xi;t", "s &xa;t", "q &X41; r", "a&amp#;b",
            "a&ampxxxxxxxxxxxxxxxx;b", "a&ampxxxxxxxxxxxxx;b", "a&amp;", "a&amp;;", "a&amp; .", "a&amp;.b",
            "TPA: &amp;x", "MAG: x&amp;", "TPA:&lt;b", "R&D; x", "a&foo;b", "a & b", "a&;", "&&amp;", "a&&#38;"]
    rng = random.Random(20261002)
    alphabet = "&#;xXamp1 F"
    for _ in range(600):
        out.append("".join(rng.choice(alphabet) for _ in range(rng.randint(1, 10))))
    seen, unique = set(), []
    for title in out:
        defline = f"id {title}".rstrip()
        if defline not in seen:
            seen.add(defline)
            unique.append(defline)
    return unique


def run(argv: list[str]) -> tuple[bytes, bytes, int]:
    result = subprocess.run(argv, capture_output=True, timeout=600)
    return result.stdout, result.stderr, result.returncode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--losat", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=6)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    record = (OUTFMT0 / "tblastx_many_subject.fasta").read_text().split(">")[1].split("\n", 1)[1]
    query_nt = OUTFMT0 / "tblastx_many_query.fasta"
    seq = "".join(query_nt.read_text().split("\n")[1:])
    index = {base: i for i, base in enumerate("TCAG")}
    peptide = "".join(CODE[index[seq[i]] * 16 + index[seq[i + 1]] * 4 + index[seq[i + 2]]]
                      for i in range(30, len(seq) - 32, 3)).replace("*", "")
    query_aa = args.work / "query.faa"
    query_aa.write_text(f">pep frame +1 of the many query\n{peptide}\n")
    queries = {"tblastx": query_nt, "tblastn": query_aa, "blastn": query_nt}
    all_titles = titles()

    def one(item):
        number, defline, program = item
        subject = args.work / f"s{number}_{program}.fna"
        subject.write_text(f">{defline}\n{record}")
        neutral = args.work / f"n{number}_{program}.fna"
        neutral.write_text(">" + defline.replace("&", "`") + "\n" + record)
        argv = ["-query", str(queries[program]), "-subject", str(subject), "-outfmt", "0"]
        ncbi = run([str(args.bin_dir / program), *argv])
        losat = run([str(args.losat.resolve()), program, *argv])
        if ncbi[2] < 0 or ncbi[2] >= 128:
            if program == "blastn":
                verdict = "ncbi-crash (exception 2)"
            elif losat[2] == 1 and b"reads past its end" in losat[1]:
                verdict = "ncbi-crash, losat-rejects"
            else:
                verdict = "UNEXPECTED"
        elif b"HTML character reference" in losat[1]:
            other = run([str(args.bin_dir / program), "-query", str(queries[program]), "-subject", str(neutral),
                         "-outfmt", "0"])
            mapped = other[0].replace(b"`", b"&").replace(str(neutral).encode(), str(subject).encode())
            verdict = "decoded, losat-rejects" if mapped != ncbi[0] else "NEEDLESS-REJECTION"
        else:
            verdict = "same" if ncbi == losat else "DIFF"
        subject.unlink()
        neutral.unlink()
        return program, defline, ncbi[2], losat[2], verdict

    items = [(n, d, p) for n, d in enumerate(all_titles) for p in ("tblastx", "tblastn", "blastn")]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(one, items))
    print("program\tdefline\tncbi_exit\tlosat_exit\tresult")
    counts: dict[str, int] = {}
    for program, defline, ncbi_exit, losat_exit, verdict in results:
        counts[verdict] = counts.get(verdict, 0) + 1
        print(f"{program}\t{defline!r}\t{ncbi_exit}\t{losat_exit}\t{verdict}")
    bad = sum(counts.get(key, 0) for key in ("UNEXPECTED", "NEEDLESS-REJECTION", "DIFF"))
    print(f"# deflines={len(all_titles)} runs={len(results)} " + " ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
