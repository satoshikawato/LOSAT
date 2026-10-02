"""E2g #1: outfmt 0 titles of punctuation deflines (approved exception 2 of PD-LOSAT-NCBI-DEFECTS).
NCBI reads past the end of these titles and crashes when such a subject has hits. LOSAT stops at
the end of the title. Check: LOSAT's report equals NCBI's report for the same subjects with the
deflines replaced by placeholders of the title's length, once each placeholder is replaced by
LOSAT's title; and NCBI crashes with the punctuation deflines."""
import random, subprocess, sys, json
from pathlib import Path
NCBI = "/home/kawato/micromamba/bin/blastn"
LOSAT = sys.argv[1]
W = Path(".")
q = "".join(open("/mnt/c/Users/genom/GitHub/LOSAT-web-gui/LOSAT/tests/fixtures/blastn_regression/inputs/resplit_query.fa").read().split("\n")[1:])[:20000]
rng = random.Random(7)
def mut(s, r): return "".join(c if rng.random() > r else rng.choice([b for b in "ACGT" if b != c]) for c in s)
punct = [", ,", "; ;", "~, ,", ",, ,", ", ,,", ";  ;", ", , ,", ",~, ,"]
titles = [", ", "; ", "~, ", ",  ", ", ", "; ", ", ", ",~, "]   # LOSAT's titles (defline.rs test)
place = ["Z1", "Z2", "Z3a", "Z4b", "Z5", "Z6", "Z7", "Z8cd"]
assert [len(t) for t in titles] == [len(p) for p in place]
segs = []
for k in range(len(punct) + 2):
    a = rng.randrange(0, len(q) - 1500); ln = rng.randint(400, 1500)
    segs.append(mut(q[a:a + ln], rng.choice((0.0, 0.03, 0.08))))
(W / "q.fa").write_text(">q20k query\n" + q + "\n")
def subjects(names, path):
    recs = [f">{n}\n{s}\n" for n, s in zip(names + ["norm1 normal", "norm2 x"], segs)]
    (W / path).write_text("".join(recs))
subjects(punct, "s_punct.fa"); subjects(place, "s_place.fa")
rows = []
for task in ("megablast", "blastn"):
    for extra in ([], ["-max_target_seqs", "3"]):
        base = ["-query", "q.fa", "-task", task, "-outfmt", "0", *extra]
        l = subprocess.run([LOSAT, "blastn", *base, "-subject", "s_punct.fa"], capture_output=True)
        n = subprocess.run([NCBI, *base, "-subject", "s_punct.fa"], capture_output=True)
        p = subprocess.run([NCBI, *base, "-subject", "s_place.fa"], capture_output=True)
        expected = p.stdout.decode()
        for pl, t in zip(place, titles):
            expected = expected.replace(pl, t)
        expected = expected.replace("s_place.fa", "s_punct.fa")
        same = l.stdout.decode() == expected
        rows.append(dict(task=task, extra=" ".join(extra), losat_exit=l.returncode, ncbi_punct_exit=n.returncode,
                         ncbi_place_exit=p.returncode, same_after_substitution=same,
                         losat_stderr=l.stderr.decode()[:200], ncbi_punct_stderr=n.stderr.decode()[-200:]))
        print(json.dumps(rows[-1]), flush=True)
        if not same:
            Path(f"diff_{task}_{len(extra)}.txt").write_text(
                "\n".join(x for x in __import__("difflib").unified_diff(expected.splitlines(), l.stdout.decode().splitlines(), lineterm="")))
json.dump(rows, open("results.json", "w"), indent=1)
