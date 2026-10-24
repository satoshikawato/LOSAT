#!/usr/bin/env python3
"""Small-case timing of LOSAT against NCBI BLAST+ 2.17.0 on one machine (A1, followups-20261014).

Standard library only; it runs with the stock Python 3 of the GitHub Linux and macOS runners.

Cases: seven of bench5's (perf round 5, 2026-10-11) small cases. NCBI searches a database built by
makeblastdb (untimed), LOSAT searches `-subject`, as in bench5. Per case, one untimed warm-up round and
--rounds timed rounds; a round runs every (configuration, thread count) once, starting at a rotating
position, so NCBI and LOSAT alternate within the round. NCBI TBLASTX runs at one thread only (NCBI 2.17
tblastx with one query and one subject searches on one thread at any -num_threads; bench5 DESIGN.md).

A configuration is NAME=BIN or NAME=BIN:VAR=VALUE[,VAR=VALUE...]; BIN names a --bin entry, and the
entry `ncbi` is the bin directory of NCBI BLAST+. Every LOSAT_* variable of the caller is dropped
before the configuration's variables are added.

Writes OUT/rows.tsv (one row per run: wall clock from perf_counter, user+sys and max RSS of the child
from wait4, exit code, SHA-256 of the output) and OUT/summary.tsv (medians of the timed rounds, and
each configuration's median divided by the first configuration's median: with `ncbi` first that is
NCBI time / configuration time, the ratio bench5 reports).

usage:
  aarch64_bench.py --out DIR --bin ncbi=DIR --bin main=PATH --config ncbi=ncbi --config losat=main
                   [--config nodp=main:LOSAT_X_DPNEON2=0] [--threads 1,max] [--rounds 5] [--cases a,b]
"""
import argparse
import hashlib
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
F = ROOT / "LOSAT" / "tests" / "fasta"

# name, program, query, subject, options, database type (bench5 scripts/bench5.py)
CASES = [
    ("blastn_small", "blastn", "AP027152.fasta", "AP027202.fasta", ["-task", "blastn"], "nucl"),
    ("megablast", "blastn", "Sakai.fna", "MG1655.fna", ["-task", "megablast"], "nucl"),
    ("blastn_big", "blastn", "AP027202.fasta", "LC738875.fasta", ["-task", "blastn"], "nucl"),
    ("blastp_wssv", "blastp", "WSSV.faa", "PajaWSV.faa", [], "prot"),
    ("tblastn_small", "tblastn", "PemoMJNVA.faa", "MelaMJNV.fasta", [], "nucl"),
    ("blastx_small", "blastx", "MelaMJNV.fasta", "PemoMJNVA.faa", [], "prot"),
    ("tbx_p03", "tblastx", "MelaMJNV.fasta", "PemoMJNVA.fasta", ["-query_gencode", "1", "-db_gencode", "1"],
     "nucl"),
]
HEADER = ["case", "config", "threads", "round", "phase", "wall_s", "cpu_s", "maxrss_kib", "rc", "sha256"]


def parse_config(spec, bins):
    name, _, rest = spec.partition("=")
    key, _, envs = rest.partition(":")
    if key not in bins:
        sys.exit(f"configuration {name}: no --bin {key}")
    env = {}
    for item in filter(None, envs.split(",")):
        var, _, value = item.partition("=")
        env[var] = value
    return name, key, env


def database(case, ncbi, work):
    name, prog, query, subject, opts, dbtype = case
    base = work / "db" / dbtype / subject
    marker = Path(str(base) + (".nsq" if dbtype == "nucl" else ".psq"))
    if not marker.exists():
        base.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run([str(ncbi / "makeblastdb"), "-in", str(F / subject), "-dbtype", dbtype, "-out", str(base)],
                       check=True, stdout=subprocess.DEVNULL)
    return base


def command(case, config, threads, bins, work, out):
    name, prog, query, subject, opts, dbtype = case
    cname, key, extra = config
    env = {k: v for k, v in os.environ.items() if not k.startswith("LOSAT_")}
    tail = ["-outfmt", "6", *opts, "-num_threads", str(threads), "-out", str(out)]
    if key == "ncbi":
        return [str(bins["ncbi"] / prog), "-query", str(F / query), "-db", str(database(case, bins["ncbi"], work)),
                *tail], env
    env.update(extra)
    return [str(bins[key]), prog, "-query", str(F / query), "-subject", str(F / subject), *tail], env


def run(cmd, env, err_path):
    with open(err_path, "wb") as err:
        start = time.perf_counter()
        proc = subprocess.Popen(cmd, env=env, stdout=subprocess.DEVNULL, stderr=err)
        _, status, usage = os.wait4(proc.pid, 0)
        wall = time.perf_counter() - start
    rss = usage.ru_maxrss // 1024 if sys.platform == "darwin" else usage.ru_maxrss  # bytes on macOS
    return wall, usage.ru_utime + usage.ru_stime, rss, os.waitstatus_to_exitcode(status)


def sha256(path):
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except FileNotFoundError:
        return "missing"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", required=True)
    parser.add_argument("--bin", action="append", required=True, metavar="KEY=PATH")
    parser.add_argument("--config", action="append", required=True, metavar="NAME=BIN[:VAR=VALUE,...]")
    parser.add_argument("--threads", default="1,max")
    parser.add_argument("--rounds", type=int, default=5)
    parser.add_argument("--cases", default="")
    args = parser.parse_args()

    out = Path(args.out).resolve()
    work = out / "work"
    work.mkdir(parents=True, exist_ok=True)
    bins = {}
    for spec in args.bin:
        key, _, path = spec.partition("=")
        bins[key] = Path(path).resolve()
    configs = [parse_config(spec, bins) for spec in args.config]
    cpus = os.cpu_count() or 1
    threads = sorted({cpus if t == "max" else int(t) for t in args.threads.split(",")})
    cases = CASES if not args.cases else [c for c in CASES if c[0] in args.cases.split(",")]
    print(f"platform {platform.platform()} machine {platform.machine()} cpus {cpus} threads {threads} "
          f"configs {[c[0] for c in configs]} cases {[c[0] for c in cases]}", flush=True)

    rows = []
    for case in cases:
        if any(c[1] == "ncbi" for c in configs):
            database(case, bins["ncbi"], work)
        slots = [(c, t) for c in configs for t in threads
                 if not (c[1] == "ncbi" and case[1] == "tblastx" and t != 1)]
        for rnd in range(args.rounds + 1):
            phase = "warmup" if rnd == 0 else "timed"
            for i in range(len(slots)):
                config, t = slots[(i + rnd) % len(slots)]
                target = work / "run.out"
                if target.exists():
                    target.unlink()
                cmd, env = command(case, config, t, bins, work, target)
                wall, cpu, rss, rc = run(cmd, env, work / f"{case[0]}.{config[0]}.n{t}.err")
                row = [case[0], config[0], str(t), str(rnd), phase, f"{wall:.4f}", f"{cpu:.4f}", str(rss), str(rc),
                       sha256(target)]
                rows.append(row)
                print("\t".join(row[:9] + [row[9][:12]]), flush=True)
    (out / "rows.tsv").write_text("\n".join("\t".join(r) for r in [HEADER, *rows]) + "\n")

    summary = [["case", "config", "threads", "median_wall_s", "min_wall_s", "max_wall_s", "median_cpu_s",
                "median_maxrss_kib", f"{configs[0][0]}_over_config", "rcs", "sha256s"]]
    for case in cases:
        first = {}
        for config in configs:
            for t in threads:
                timed = [r for r in rows if r[0] == case[0] and r[1] == config[0] and r[4] == "timed"
                         and (r[2] == str(t) or (config[1] == "ncbi" and case[1] == "tblastx"))]
                if not timed:
                    continue
                walls = [float(r[5]) for r in timed]
                med = statistics.median(walls)
                if config is configs[0]:
                    first[t] = med
                ratio = first.get(t, first.get(1))
                summary.append([case[0], config[0], str(t), f"{med:.4f}", f"{min(walls):.4f}", f"{max(walls):.4f}",
                                f"{statistics.median(float(r[6]) for r in timed):.4f}",
                                str(int(statistics.median(int(r[7]) for r in timed))),
                                f"{ratio / med:.3f}" if ratio else "", ",".join(sorted({r[8] for r in timed})),
                                ",".join(sorted({r[9][:12] for r in timed}))])
    text = "\n".join("\t".join(r) for r in summary) + "\n"
    (out / "summary.tsv").write_text(text)
    print(text)
    step_summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if step_summary:
        with open(step_summary, "a") as handle:
            handle.write(f"### {platform.machine()} {platform.platform()}\n\n```\n{text}```\n")
    return 1 if any(r[8] != "0" for r in rows) else 0


if __name__ == "__main__":
    sys.exit(main())
