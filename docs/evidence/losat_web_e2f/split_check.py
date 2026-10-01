#!/usr/bin/env python3
"""Compare BLASTN query splitting with NCBI (Session S07++; comparison only).

NCBI searches a query batch of at least two query chunks chunk by chunk in its preliminary
stage and merges the HSPs of the chunks (`CQuerySplitter`, AUTHORITY.md section D). This
script writes queries that NCBI splits and subjects around the chunk boundaries, runs each
case in NCBI, in NCBI with CHUNK_SIZE=20000000 (not split) and in LOSAT, and compares LOSAT
with NCBI byte for byte (stdout, stderr and exit status). NCBI's stderr warning
"'num_threads' is currently ignored when 'subject' is specified." is removed before the
comparison (LOSAT does not write it; docs/evidence/losat_web_e2c/AUTHORITY.md section I).
The last column says whether splitting changes NCBI's output. Exits 1 when a case differs.

The cases named r1_* are the round-1 audit reproductions (`make_regression_cases`): NCBI keeps
only `prelim_hitlist_size` subjects' HSP lists (550 by default, 10 for -max_target_seqs 1..5)
after the preliminary stage and before the traceback, so output differs from a build that
applies that limit after the traceback when more subjects than the limit have hits. (a) 560
subjects plus one that spans the chunk boundary of a split query; (b) -max_target_seqs 1, 3
and the controls 10, 11; (c) the same on a split megablast query; (d) an unsplit batch.
They use their own random stream, so the other cases do not change. Run them alone with
`--only r1_`.

Usage: split_check.py --bin-dir DIR --losat LOSAT --work DIR [--seed N] [--jobs N] [--only PREFIX]
"""
from __future__ import annotations

import argparse
import os
import random
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
FASTA = HERE.parents[2] / "LOSAT" / "tests" / "fasta"
NUM_THREADS_WARNING = b"Warning: [blastn] 'num_threads' is currently ignored when 'subject' is specified.\n"
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def genome(name: str) -> str:
    return "".join(line.strip() for line in (FASTA / name).read_text().splitlines()[1:])


def chunks(lengths: list[int], chunk_size: int, overlap: int = 100) -> list[tuple[int, int]]:
    """NCBI's query chunks on the concatenated queries (split_query_aux_priv.cpp:123-138,
    split_query_cxx.cpp:145-171)."""
    total = sum(lengths)
    count = total // (chunk_size - overlap)
    if count <= 1:
        return [(0, total)]
    size = (total + (count - 1) * overlap) // count
    if count < size - overlap:
        size += 1
    out, start = [], 0
    for k in range(count):
        end = start + size
        if end >= total or k + 1 == count:
            end = total
        out.append((start, end))
        start += size - overlap
        if start > total or end == total:
            break
    return out


def write(path: Path, records: list[tuple[str, str]]) -> str:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
    return path.name


def diverged(seq: str, first_exact: int) -> str:
    """Exact for `first_exact` letters, then a mismatch every 7: no other 8-letter word."""
    out = list(seq)
    for i in range(first_exact, len(out), 7):
        out[i] = {"A": "C", "C": "G", "G": "T", "T": "A"}.get(out[i].upper(), "A")
    return "".join(out)


def mutate(rng: random.Random, seq: str) -> str:
    divergence = rng.choice([0, 0.02, 0.08, 0.15, 0.25])
    out = []
    for letter in seq:
        r = rng.random()
        if r < divergence * 0.8:
            out.append(rng.choice("ACGT"))
        elif r < divergence * 0.9:
            continue
        elif r < divergence:
            out += [letter, rng.choice("ACGT")]
        else:
            out.append(letter)
    if rng.random() < 0.5:
        out = [COMPLEMENT.get(c, "N") for c in reversed(out)]
    return "".join(out)


def make_cases(work: Path, seed: int) -> list[tuple[str, str, str, list[str]]]:
    edl, sak = genome("EDL933.fna"), genome("Sakai.fna")
    rng = random.Random(seed)
    cases = []
    write(work / "edl933.fa", [("EDL933", edl)])
    # Two queries are two batches (the first batch is about 5000 residues, blast_input.cpp);
    # one record of both genomes is one batch, which megablast splits in two.
    write(work / "two.fa", [("EDL933", edl), ("Sakai", sak)])
    both = edl + sak
    write(work / "cat.fa", [("EDL933_Sakai", both)])

    # EDL933 with -task blastn (five chunks): windows around the starts of chunks 1 to 4.
    for k, (start, _) in enumerate(chunks([len(edl)], 1_000_000)[1:], 1):
        for width, shift in ((300, 50), (300, -50), (300, -150), (3000, 50), (3000, -50),
                             (3000, -1500), (20000, 50), (20000, -50), (20000, -10000)):
            a = start + shift - width // 2
            name = f"edl{k}_{width}_{shift}"
            cases.append((name, "edl933.fa", write(work / f"{name}.fa", [(name, edl[a:a + width])]), ["-task", "blastn"]))

    # Megablast on one record of both genomes (two chunks), around the chunk boundary and the
    # junction of the genomes; -task blastn on the two genomes as two queries, around the
    # chunk boundaries of the second batch (Sakai, five chunks).
    plans = [("cat_megablast", "cat.fa", both, chunks([len(both)], 5_000_000), [len(edl)], "megablast"),
             ("two_blastn", "two.fa", sak, chunks([len(sak)], 1_000_000), [], "blastn")]
    for tag, query, seq, bounds, extra_points, task in plans:
        points = sorted({end for _, end in bounds[:-1]} | {start for start, _ in bounds[1:]} | set(extra_points))
        for p in points:
            for width, shift in ((3000, -1500), (20000, -10000), (300, -50)):
                a = max(0, p + shift)
                name = f"{tag}_{p}_{width}"
                cases.append((name, query, write(work / f"{name}.fa", [(name, seq[a:a + width])]), ["-task", task]))

    # An invalid query before and after a split query.
    q = edl[:2_000_000]
    write(work / "inv_first.fa", [("allN", "N" * 200), ("big", q)])
    write(work / "inv_last.fa", [("big", q), ("allN", "N" * 200)])
    for name, a in (("w999k", 999_000), ("w1500k", 1_500_000), ("w1900k", 1_900_000)):
        subject = write(work / f"{name}.fa", [(name, edl[a:a + 4000])])
        for query in ("inv_first", "inv_last"):
            for fmt in ("6", "0", "7"):
                cases.append((f"{query}_{name}_{fmt}", f"{query}.fa", subject, ["-task", "blastn", "-outfmt", fmt]))

    # A short query first: the second chunk has only the long query, whose contexts take the
    # short query's search spaces.
    write(work / "short_first.fa", [("short", edl[3_000_000:3_001_000]), ("long", sak[:2_000_000])])
    for a in (1_100_000, 1_400_000, 1_800_000, 1_990_000):
        seq = list(sak[a:a + 5000])
        for i in range(0, len(seq), 9):
            if rng.random() < 0.5:
                seq[i] = rng.choice("ACGT")
        subject = write(work / f"sak_{a}.fa", [(f"sak_{a}", "".join(seq))])
        cases.append((f"short_first_{a}", "short_first.fa", subject, ["-task", "blastn"]))
        cases.append((f"short_first_{a}_threads", "short_first.fa", subject, ["-task", "blastn", "-num_threads", "4"]))

    # Lower-case masks restricted to a chunk: a mask grows by one residue (the only word, which
    # starts right after the mask, is masked in the chunk), and a mask that starts at the
    # chunk's last residue is dropped (the only word, which ends there, is not masked).
    last = chunks([2_000_000], 1_000_000)[0][1] - 1
    b = 700_000
    lq = list(edl[:2_000_000])
    for i in list(range(b - 60, b + 1)) + list(range(last, last + 50)):
        lq[i] = lq[i].lower()
    write(work / "lcase.fa", [("lcase", "".join(lq))])
    flank = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    write(work / "lcase_grow.fa", [("grow", flank(200) + diverged(edl[b + 1:b + 301], 8) + flank(200))])
    write(work / "lcase_drop.fa", [("drop", flank(200) + diverged(edl[last - 299:last + 1][::-1], 8)[::-1] + flank(200))])
    for subject in ("lcase_grow", "lcase_drop"):
        for word in ("8", "11"):
            cases.append((f"{subject}_w{word}", "lcase.fa", f"{subject}.fa", ["-task", "blastn", "-lcase_masking", "-word_size", word]))
            cases.append((f"{subject}_w{word}_nolc", "lcase.fa", f"{subject}.fa", ["-task", "blastn", "-word_size", word]))

    # Several subjects.
    write(work / "multi_subj.fa", [(f"s{i}", edl[a:a + 3000]) for i, a in enumerate(range(990_000, 1_010_000, 2000))])
    for extra in ([], ["-max_target_seqs", "3"], ["-subject_besthit"], ["-num_threads", "4"]):
        cases.append(("multi" + "".join("_" + x.strip("-") for x in extra), "inv_last.fa", "multi_subj.fa", ["-task", "blastn", *extra]))

    # A split batch followed by more batches (the split batch has 0 initial hits).
    write(work / "split_then_more.fa", [("big", edl[:2_000_000])] + [(f"small{i}", sak[1_000_000 + i * 7000:1_001_500 + i * 7000]) for i in range(8)])
    write(work / "split_then_more_subj.fa", [("hit_big", edl[1_500_000:1_504_000]), ("hit_small", sak[1_000_000:1_060_000])])
    for task in ("blastn", "megablast"):
        for fmt in ("6", "7", "0"):
            cases.append((f"split_then_more_{task}_{fmt}", "split_then_more.fa", "split_then_more_subj.fa", ["-task", task, "-outfmt", fmt]))

    # Random subjects around the chunk boundaries: diverged copies with indels, either strand.
    plans = [("edl_blastn", "edl933.fa", edl, [len(edl)], 1_000_000, "blastn", 40),
             ("cat_megablast", "cat.fa", both, [len(both)], 5_000_000, "megablast", 40),
             ("two_blastn", "two.fa", sak, [len(sak)], 1_000_000, "blastn", 20)]
    for tag, query, seq, lengths, size, task, count in plans:
        bounds = chunks(lengths, size)
        points = [end for _, end in bounds[:-1]] + [start for start, _ in bounds[1:]]
        for k in range(count):
            p = rng.choice(points) + rng.randint(-150, 150)
            width = rng.choice([150, 300, 600, 1500, 5000, 20000])
            a = max(0, p - rng.randint(0, width))
            name = f"sweep_{tag}_{k}"
            cases.append((name, query, write(work / f"{name}.fa", [(name, mutate(rng, seq[a:a + width]))]), ["-task", task]))
    cases += make_regression_cases(work, seed)
    return cases


def make_regression_cases(work: Path, seed: int) -> list[tuple[str, str, str, list[str]]]:
    """Round-1 audit reproductions of the prelim_hitlist_size limit (see the module docstring).
    Own random stream: adding or changing these cases does not change `make_cases`."""
    edl, sak = genome("EDL933.fna"), genome("Sakai.fna")
    rng = random.Random(f"{seed}-regression")
    cases = []
    q = edl[:2_000_000]
    write(work / "r1_q2m.fa", [("EDL933_2M", q)])

    # (a) 560 subjects with hits in the first chunk plus B, which spans the chunk boundary
    # (NCBI prints 500 subjects and not B).
    subjects = [("B", q[999_200:1_000_800])]
    for i in range(560):
        a = rng.randrange(0, 990_000 - 1500 + 1)
        seq = list(q[a:a + 1500])
        for j in range(len(seq)):
            if rng.random() < 0.03:
                seq[j] = rng.choice([c for c in "ACGT" if c != seq[j]])
        subjects.append((f"S{i}", "".join(seq)))
    cases.append(("r1_a_default", "r1_q2m.fa", write(work / "r1_a.fa", subjects), ["-task", "blastn", "-outfmt", "6"]))

    # (b) X spans the chunk boundary, Y0..Y9 are far from it.
    subject = write(work / "r1_b.fa", [("X", edl[997_000:1_003_000])] + [(f"Y{i}", edl[50_000 + 80_000 * i:50_000 + 80_000 * i + 4500]) for i in range(10)])
    for n in (1, 3, 10, 11):
        cases.append((f"r1_b_mts{n}", "r1_q2m.fa", subject, ["-task", "blastn", "-max_target_seqs", str(n)]))

    # (c) One record of both genomes (megablast splits it in two); X spans the chunk boundary.
    both = edl + sak
    write(work / "r1_cat.fa", [("EDL933_Sakai", both)])
    subject = write(work / "r1_c.fa", [("X", both[5_510_500:5_516_500])] + [(f"Y{i}", both[200_000 + 900_000 * i:200_000 + 900_000 * i + 4500]) for i in range(10)])
    cases.append(("r1_c_megablast_mts1", "r1_cat.fa", subject, ["-task", "megablast", "-max_target_seqs", "1"]))

    # (d) Unsplit batch; X is interrupted by 25 random letters.
    write(work / "r1_q300k.fa", [("EDL933_300k", edl[:300_000])])
    insert = "".join(rng.choice("ACGT") for _ in range(25))
    subject = write(work / "r1_d.fa", [("X", edl[100_000:102_000] + insert + edl[102_000:104_000])] + [(f"Y{i}", edl[150_000 + 12_000 * i:150_000 + 12_000 * i + 3000]) for i in range(10)])
    cases.append(("r1_d_mts1", "r1_q300k.fa", subject, ["-task", "blastn", "-max_target_seqs", "1", "-outfmt", "6"]))
    return cases


def run(argv: list[str], work: Path, env: dict | None = None) -> tuple[int, bytes, bytes]:
    proc = subprocess.run(argv, cwd=work, capture_output=True, env=env)
    return proc.returncode, proc.stdout, proc.stderr


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin-dir", required=True, type=Path)
    parser.add_argument("--losat", required=True, type=Path)
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=20260930)
    parser.add_argument("--jobs", type=int, default=6)
    parser.add_argument("--only", default="", metavar="PREFIX", help="run only the cases whose name starts with PREFIX")
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    cases = [case for case in make_cases(args.work, args.seed) if case[0].startswith(args.only)]
    ncbi = str(args.bin_dir / "blastn")
    unsplit = dict(os.environ, CHUNK_SIZE="20000000")

    def one(case):
        name, query, subject, extra = case
        if "-outfmt" not in extra:
            extra = [*extra, "-outfmt", "6"]
        argv = ["-query", query, "-subject", subject, *extra]
        expected = run([ncbi, *argv], args.work)
        whole = run([ncbi, *argv], args.work, unsplit)
        got = run([str(args.losat), "blastn", *argv], args.work)
        same = (expected[0], expected[1], expected[2].replace(NUM_THREADS_WARNING, b"")) == got
        return name, same, expected[1] != whole[1], " ".join(argv)

    with ThreadPoolExecutor(args.jobs) as pool:
        results = list(pool.map(one, cases))
    print("# case\tresult\tsplit_changes_ncbi\targv")
    for name, same, changed, argv in results:
        print(f"{name}\t{'same' if same else 'DIFFERENT'}\t{'yes' if changed else 'no'}\t{argv}")
    differing = sum(not same for _, same, _, _ in results)
    changed = sum(changed for _, _, changed, _ in results)
    print(f"# cases={len(results)} differing={differing} split_changes_ncbi={changed}")
    return 1 if differing else 0


if __name__ == "__main__":
    sys.exit(main())
