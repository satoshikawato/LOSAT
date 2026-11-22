#!/usr/bin/env python3
"""Build frozen benchmark snapshot tables in benchmarks/v0.2.0/snapshot/ from measured TSV/JSON outputs.

Every input path is an argument (defaults: the v0.2.0 runs on the release candidate 93ca2f652).
metadata.json records the sha256 of every input file read and the LOSAT build of each source:
--losat-commit for the standard grid and its extension (psiblast, rpsblast, rpstblastn, utilities,
thread and size scaling), the `build` field of toolmt.json (downstream tools at the tool's thread
count, cold-cache and core_nt rows), <dt-dir>/<tag>_sha.txt (downstream tools at one thread) and the
`commit` column of the WebAssembly CSV.
"""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import os
import statistics
import sys
from collections import defaultdict
from pathlib import Path

MODES_ORDER = [
    # group, case_id, program, task, label, query_desc, subj_desc, target, options
    ("Nucleotide vs nucleotide (blastn tasks)", "mb", "blastn", "megablast", "megablast",
     "16S rRNA genes, full length", "RefSeq 16S rRNA", "-subject", "-"),
    ("Nucleotide vs nucleotide (blastn tasks)", "dc", "blastn", "dc-megablast", "dc-megablast",
     "D. melanogaster cDNA", "D. melanogaster chr. 2L (3 Mb)", "-subject", "-task dc-megablast"),
    ("Nucleotide vs nucleotide (blastn tasks)", "bn", "blastn", "blastn", "blastn",
     "D. melanogaster cDNA", "D. melanogaster chr. 2L (3 Mb)", "-subject", "-task blastn"),
    ("Nucleotide vs nucleotide (blastn tasks)", "bns", "blastn", "blastn-short", "blastn-short",
     "Synthetic 24-mer reads", "D. melanogaster chr. 2L (3 Mb)", "-subject", "-task blastn-short -evalue 0.001"),
    ("Protein vs protein (blastp tasks)", "bp", "blastp", "blastp", "blastp",
     "Human proteome (first 3,000)", "S. cerevisiae proteome", "-subject", "-evalue 0.001"),
    ("Protein vs protein (blastp tasks)", "bpff", "blastp", "blastp-fast", "blastp-fast",
     "Human proteome", "S. cerevisiae proteome", "-subject", "-task blastp-fast -evalue 0.001"),
    ("Translated", "bx", "blastx", "blastx", "blastx",
     "D. melanogaster chr. 2L (50 kb)", "D. melanogaster proteome", "-db", "-evalue 1e-5"),
    ("Translated", "tbn", "tblastn", "tblastn", "tblastn",
     "Human protein ENSP00000343742.5", "Human genome GRCh38", "-db", "-evalue 1e-5 -max_target_seqs 5"),
    ("Translated", "tbx", "tblastx", "tblastx", "tblastx",
     "Candidatus genomes (3 seqs)", "Candidatus genomes (3 seqs)", "-subject", "-query_gencode 4 -db_gencode 4"),
    ("PSI & Conserved Domain (Profile) Searches", "psi", "psiblast", "psiblast", "psiblast",
     "E. coli MG1655 proteins (20)", "UniProtKB/Swiss-Prot", "-db", "-num_iterations 3 -evalue 0.001"),
    ("PSI & Conserved Domain (Profile) Searches", "rps", "rpsblast", "rpsblast", "rpsblast",
     "E. coli MG1655 proteins (300)", "NCBI CDD v3.22 (29 vol.)", "-db", "-evalue 0.01"),
    ("PSI & Conserved Domain (Profile) Searches", "rpstbn", "rpstblastn", "rpstblastn", "rpstblastn",
     "E. coli MG1655 5-kb segments", "NCBI CDD v3.22 (29 vol.)", "-db", "-"),
]

# Single-thread tool benchmarks from official test sets (DT1 / DT13 / DT15)
ST_TOOLS = [
    ("blastn", "ResFinder / PlasmidFinder", "1", "dt15", ["CGE_core.cgecore", "CGE_pf.cgecore", "CGE_res01.cgecore", "CGE_res03.cgecore", "CGE_res05.cgecore"]),
    ("blastn", "Bakta", "1", "dt15", ["BAK_oric.Bakta", "BAK_orit.Bakta", "BAK_toy.Bakta"]),
    ("blastn", "Easyfig", "1", "dt15", ["EFG_ex1.Easyfig", "EFG_ex2a.Easyfig", "EFG_ex2b.Easyfig"]),
    ("blastn", "SISTR", "1", "dt15", ["SIS_cent.SISTR", "SIS_fliC.SISTR", "SIS_full.SISTR", "SIS_wzx.SISTR"]),
    ("blastn", "pyani ANIb", "1", "dt15", [
        "PYA_NC_002696_vs_NC_010338.pyani", "PYA_NC_002696_vs_NC_011916.pyani",
        "PYA_NC_002696_vs_NC_014100.pyani", "PYA_NC_010338_vs_NC_002696.pyani",
        "PYA_NC_010338_vs_NC_011916.pyani", "PYA_NC_010338_vs_NC_014100.pyani",
        "PYA_NC_011916_vs_NC_002696.pyani", "PYA_NC_011916_vs_NC_010338.pyani",
        "PYA_NC_011916_vs_NC_014100.pyani", "PYA_NC_014100_vs_NC_002696.pyani",
        "PYA_NC_014100_vs_NC_010338.pyani", "PYA_NC_014100_vs_NC_011916.pyani",
    ]),
    ("blastn", "ECTyper", "1", "dt15", ["ECT02.ECTyper", "ECT04.ECTyper"]),
    ("blastn", "rBLAST (R)", "1", "dt15", ["RBL04.rBLAST", "RBL05.rBLAST"]),
    ("blastx", "Biopython", "1", "dt1", ["BIO01.Biopython"]),
    ("tblastn", "Unicycler", "1", "dt13", ["UNI02.Unicycler"]),
    ("tblastx", "Easyfig", "1", "dt1", ["EFG_ex1tx.Easyfig"]),
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


RunRecord = tuple[float, float, str, float | None]


def sanitize_input_key(path: Path, base: Path) -> str:
    """Record provenance key relative to base without machine-local home prefixes."""
    resolved_base = base.resolve()
    resolved_path = path.resolve()
    try:
        rel = resolved_path.relative_to(resolved_base)
    except ValueError:
        sys.exit(f"input path is outside baselines root {resolved_base}: {resolved_path}")
    return f"losat-baselines/{rel.as_posix()}"


def collect_input_checksums(paths: list[Path], base: Path) -> dict[str, str]:
    """Validate input paths against base, reject duplicate keys, and hash existing files."""
    checksums: dict[str, str] = {}
    for p in paths:
        if not p.is_file():
            continue
        key = sanitize_input_key(p, base)
        if key in checksums:
            sys.exit(f"duplicate sanitized input_checksums key: {key} (from {p})")
        checksums[key] = sha256_file(p)
    return checksums


def parse_rss_anon_mib(row: dict[str, str]) -> float | None:
    """Parse optional anonymous private RSS (MiB) from raw benchmark rows."""
    if row.get("rss_anon_mib"):
        return float(row["rss_anon_mib"])
    if row.get("rss_anon_kb"):
        return int(row["rss_anon_kb"]) / 1024.0
    return None


def make_standard_row(
    *,
    group: str,
    cid: str,
    prog: str,
    task: str,
    label: str,
    th: int,
    qdesc: str,
    sdesc: str,
    target: str,
    opts: str,
    lo_list: list[RunRecord],
    nc_list: list[RunRecord],
    v0_list: list[RunRecord],
    is_ncbi_timeout: bool = False,
    subject_vs_db: bool = False,
    ncbi_n1_sha: set[str] | None = None,
) -> dict[str, object]:
    """Assemble and validate one standard workload summary row."""
    if not lo_list:
        sys.exit(f"missing LOSAT runs for {cid} t={th}")
    lo_w = statistics.median(x[0] for x in lo_list)
    lo_w_min = min(x[0] for x in lo_list)
    lo_w_max = max(x[0] for x in lo_list)
    lo_m = statistics.median(x[1] for x in lo_list)
    lo_sha = {x[2] for x in lo_list}
    if len(lo_sha) != 1:
        sys.exit(f"LOSAT sha differs across rounds for {cid} t={th}: {lo_sha}")

    nc_w = statistics.median(x[0] for x in nc_list) if nc_list else None
    nc_m = statistics.median(x[1] for x in nc_list) if nc_list else None
    nc_sha = {x[2] for x in nc_list}
    v0_w = statistics.median(x[0] for x in v0_list) if v0_list else ""
    v0_m = round(statistics.median(x[1] for x in v0_list), 2) if v0_list else ""

    parity = "same"
    ncbi_runs = len(nc_list)
    if not nc_list and th == 1 and is_ncbi_timeout:
        parity = "ncbi_timeout_6000s"
        ncbi_runs = 1
    elif subject_vs_db:
        parity = "subject_vs_db"
    elif ncbi_n1_sha is not None and th == 16 and lo_sha != nc_sha:
        parity = "same_as_ncbi_n1" if lo_sha == ncbi_n1_sha else "diff"
    elif lo_sha != nc_sha:
        parity = "diff"

    row_entry: dict[str, object] = {
        "group": group,
        "case": cid,
        "program": prog,
        "task": task,
        "label": label,
        "threads": th,
        "losat_wall_s": round(lo_w, 3),
        "losat_wall_min_s": round(lo_w_min, 3),
        "losat_wall_max_s": round(lo_w_max, 3),
        "losat_rss_mib": round(lo_m, 2),
        "ncbi_wall_s": round(nc_w, 3) if nc_w is not None else "",
        "ncbi_rss_mib": round(nc_m, 2) if nc_m is not None else "",
        "v010_wall_s": round(v0_w, 3) if v0_w != "" else "",
        "v010_rss_mib": v0_m,
        "speedup": round(nc_w / lo_w, 4) if nc_w else "",
        "mem_reduction": round(nc_m / lo_m, 4) if nc_m else "",
        "losat_runs": len(lo_list),
        "ncbi_runs": ncbi_runs,
        "losat_sha16": next(iter(lo_sha)),
        "ncbi_sha16": ",".join(sorted(nc_sha)),
        "parity": parity,
        "query_desc": qdesc,
        "subject_desc": sdesc,
        "target": target,
        "options": opts,
    }
    lo_anon = [x[3] for x in lo_list if x[3] is not None]
    nc_anon = [x[3] for x in nc_list if x[3] is not None]
    if lo_anon:
        row_entry["losat_rss_anon_mib"] = round(statistics.median(lo_anon), 2)
    if nc_anon:
        row_entry["ncbi_rss_anon_mib"] = round(statistics.median(nc_anon), 2)
    return row_entry


def build_downstream_tools(toolmt_path: Path, dt_dir: Path) -> tuple[list[dict[str, object]], list[str]]:
    tm = json.loads(toolmt_path.read_text(encoding="utf-8"))
    by_prog_tool: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for r in tm:
        cid = r["case"]
        if cid.startswith(("M2_", "M3_", "M4_", "UR50", "SP")):
            continue
        if r["tool"] == "ABRicate":
            if r["threads"] != 8:
                continue
        elif r["threads"] != 16:
            continue
        tool_name = "nf-core blast" if r["tool"] == "nf-core" else r["tool"]
        by_prog_tool[(r["program"], tool_name)].append(r)

    out: list[dict[str, object]] = []
    for (prog, tool), rs in sorted(by_prog_tool.items()):
        form = rs[0]["form"]
        if form == "プロセス並列":
            th_label = "16 × 1"
            th_num = 16
        else:
            th_num = rs[0]["threads"]
            th_label = str(th_num)

        ncbi_mt = sum(float(r["ncbi_s"]) for r in rs)
        losat_mt = sum(float(r["losat_s"]) for r in rs)
        mt_x = round(ncbi_mt / losat_mt, 4) if losat_mt > 0 else ""

        has_n1 = all(r.get("n1_ncbi_s") is not None and r.get("n1_losat_s") is not None for r in rs)
        if has_n1:
            ncbi_n1 = sum(float(r["n1_ncbi_s"]) for r in rs)
            losat_n1 = sum(float(r["n1_losat_s"]) for r in rs)
            n1_x: object = round(ncbi_n1 / losat_n1, 4) if losat_n1 > 0 else ""
        else:
            ncbi_n1 = None
            losat_n1 = None
            n1_x = ""

        note = ""
        if prog == "blastn" and tool == "Galaxy" and losat_mt == 0.0:
            note = "losat_under_10ms"
        out.append({
            "program": prog,
            "tool": tool,
            "threads_label": th_label,
            "threads_num": th_num,
            "cases_count": len(rs),
            "ncbi_mt_s": round(ncbi_mt, 3),
            "losat_mt_s": round(losat_mt, 3),
            "mt_speedup": mt_x,
            "ncbi_n1_s": round(ncbi_n1, 3) if ncbi_n1 is not None else "",
            "losat_n1_s": round(losat_n1, 3) if losat_n1 is not None else "",
            "n1_speedup": n1_x,
            "same_output": "yes" if all(r["same_output"] == "はい" for r in rs) else "no",
            "note": note,
        })

    # Load dt1/dt13/dt15 for single-thread tools
    dt_cache: dict[str, dict[tuple[str, str], list[float]]] = {}
    dt_meta: dict[str, dict[tuple[str, str], set[str]]] = {}
    wanted: dict[str, set[str]] = defaultdict(set)
    for _prog, _tool, _th, tag, case_ids in ST_TOOLS:
        wanted[tag].update(case_ids)
    for tag, fn in (("dt1", "dt1.tsv"), ("dt13", "dt13.tsv"), ("dt15", "dt15.tsv")):
        p = dt_dir / fn
        g: dict[tuple[str, str], list[float]] = defaultdict(list)
        s: dict[tuple[str, str], set[str]] = defaultdict(set)
        for row in read_tsv(p):
            if row["case"] in wanted.get(tag, set()) and row.get("rc", "0") != "0":
                sys.exit(f"{tag} rc != 0: {row}")
            g[(row["case"], row["arm"])].append(float(row["wall_s"]))
            s[(row["case"], row["arm"])].add(row.get("sha16", ""))
        dt_cache[tag] = g
        dt_meta[tag] = s

    for prog, tool, th_label, tag, case_ids in ST_TOOLS:
        g = dt_cache[tag]
        s = dt_meta[tag]
        nc_tot = 0.0
        lo_tot = 0.0
        shas_match = True
        for cid in case_ids:
            nc_vals = g.get((cid, "ncbi"))
            lo_vals = g.get((cid, "losat"))
            if not nc_vals or not lo_vals:
                sys.exit(f"missing {cid} in {tag}")
            nc_tot += statistics.median(nc_vals)
            lo_tot += statistics.median(lo_vals)
            if s.get((cid, "ncbi")) != s.get((cid, "losat")):
                shas_match = False
        out.append({
            "program": prog,
            "tool": tool,
            "threads_label": th_label,
            "threads_num": 1,
            "cases_count": len(case_ids),
            "ncbi_mt_s": "",
            "losat_mt_s": "",
            "mt_speedup": "",
            "ncbi_n1_s": round(nc_tot, 3),
            "losat_n1_s": round(lo_tot, 3),
            "n1_speedup": round(nc_tot / lo_tot, 4),
            "same_output": "yes" if shas_match else "no",
            "note": f"1t_only_{tag}",
        })

    sources = [str(toolmt_path), str(dt_dir / "dt1.tsv"), str(dt_dir / "dt13.tsv"), str(dt_dir / "dt15.tsv")]
    return out, sources


def build_webassembly(wasm_csv: Path, j_tblastx_ncbi: dict[int, tuple[float, float]] | None) -> list[dict[str, object]]:
    rows = []
    with wasm_csv.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            mode_order = int(r["mode_order"])
            if mode_order > 8:
                continue
            c = r["mode"]
            th = int(r["threads"])
            nat = float(r["native_s"])
            was = float(r["wasm_s"])
            ncbi = float(r["ncbi_s"]) if r["ncbi_s"] else None
            note = f"commit_{r['commit']}"
            if c == "tbx" and ncbi is None and j_tblastx_ncbi and th in j_tblastx_ncbi:
                ncbi, _ = j_tblastx_ncbi[th]
                note = f"commit_{r['commit']}_ncbi_j2"
            rows.append({
                "case": c,
                "mode_label": r["mode_label"],
                "threads": th,
                "commit": r["commit"],
                "native_s": round(nat, 3),
                "wasm_s": round(was, 3),
                "ncbi_s": round(ncbi, 3) if ncbi is not None else "",
                "native_speedup": round(ncbi / nat, 4) if ncbi is not None else "",
                "wasm_speedup": round(ncbi / was, 4) if ncbi is not None else "",
                "wasm_overhead": round(was / nat, 4),
                "native_mib": r["native_mib"],
                "wasm_mib": r["wasm_mib"],
                "ncbi_mib": r["ncbi_mib"] or "",
                "note": note,
            })
    return rows


def main() -> int:
    base = Path(os.environ.get("LOSAT_BASELINES_DIR", Path.home() / "losat-baselines"))
    ap = argparse.ArgumentParser()
    ap.add_argument("--j2-dir", type=Path, default=base / "release-v020-20261025/work/j-speed/results2")
    ap.add_argument("--ext-dir", type=Path, default=base / "release-v020-20261025/work/bench-v020-ext/results")
    ap.add_argument("--toolmt", type=Path, default=base / "tool-mt-20261116/results/toolmt.json")
    ap.add_argument("--dt-dir", type=Path, default=base / "dropin-survey-20261029/perf/ws_results")
    ap.add_argument("--wasm-csv", type=Path, default=base / "v020-board-20261030/compact/data/wasm_scale.csv")
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent.parent / "snapshot")
    ap.add_argument("--losat-commit", default="93ca2f652617acac09d416c28c293cc8a1a501d1",
                    help="full SHA of the LOSAT build measured in --j2-dir and --ext-dir")
    ap.add_argument("--losat-version", default="0.2.0")
    ap.add_argument("--wasm-paired", action="store_true",
                    help="mark WebAssembly measurements as direct paired runs (sets wasm_derived=false in metadata.json)")
    ap.add_argument("--wasm-derived-from", default="6a0a8719",
                    help="source LOSAT build recorded as wasm_derived_from when --wasm-paired is not set")
    ap.add_argument("--static-only", action="store_true", help="Only build downstream_tools (and webassembly if complete).")
    args = ap.parse_args()

    dt_source_paths = [args.toolmt, args.dt_dir / "dt1.tsv", args.dt_dir / "dt13.tsv", args.dt_dir / "dt15.tsv"]
    static_required = [*dt_source_paths, args.wasm_csv]

    if args.static_only:
        for req in static_required:
            if not req.is_file():
                sys.exit(f"missing completed input file: {req}")
        collect_input_checksums(static_required, base)
        dt_rows, _dt_sources = build_downstream_tools(args.toolmt, args.dt_dir)
        wasm_rows = build_webassembly(args.wasm_csv, None)
        os.makedirs(args.out, exist_ok=True)
        write_tsv(
            args.out / "downstream_tools.tsv",
            ["program", "tool", "threads_label", "threads_num", "cases_count",
             "ncbi_mt_s", "losat_mt_s", "mt_speedup", "ncbi_n1_s", "losat_n1_s",
             "n1_speedup", "same_output", "note"],
            dt_rows,
        )
        write_tsv(
            args.out / "webassembly.tsv",
            ["case", "mode_label", "threads", "commit", "native_s", "wasm_s", "ncbi_s",
             "native_speedup", "wasm_speedup", "wasm_overhead", "native_mib", "wasm_mib",
             "ncbi_mib", "note"],
            wasm_rows,
        )
        print("static snapshot written:", args.out)
        return 0

    wasm_derived_from = args.wasm_derived_from.strip()
    if not args.wasm_paired and not wasm_derived_from:
        sys.exit("--wasm-derived-from must be non-empty unless --wasm-paired is set")

    # 1. Resolve all input paths, check existence, and validate input_checksums before any write
    j_n16_path = args.j2_dir / "n16/runs.tsv"
    j_n1_paths = sorted((args.j2_dir / "n1").glob("perfstat_*.tsv"))
    j_inputs_path = args.j2_dir / "inputs.tsv"
    ext_n0_path = args.ext_dir / "n0.tsv"
    ext_n1_path = args.ext_dir / "n1.tsv"
    ext_inputs_path = args.ext_dir / "inputs_ext.tsv"
    ext_bp_path = args.ext_dir / "bp_db_parity.tsv"
    dt_sha_paths = [args.dt_dir / f"{tag}_sha.txt" for tag in ("dt1", "dt13", "dt15")]

    for req in (j_n16_path, j_inputs_path, ext_n0_path, ext_n1_path, ext_inputs_path, *static_required):
        if not req.is_file():
            sys.exit(f"missing completed input file: {req}")
    if not j_n1_paths:
        sys.exit(f"missing completed input files: {args.j2_dir / 'n1'}/perfstat_*.tsv")

    input_files = [
        j_n16_path,
        *j_n1_paths,
        j_inputs_path,
        ext_n0_path,
        ext_n1_path,
        ext_inputs_path,
        ext_bp_path,
        args.wasm_csv,
        *dt_source_paths,
        *dt_sha_paths,
    ]
    input_checksums = collect_input_checksums(input_files, base)

    dt_rows, _dt_sources = build_downstream_tools(args.toolmt, args.dt_dir)

    j_n16: dict[tuple[str, str], list[RunRecord]] = defaultdict(list)
    for r in read_tsv(j_n16_path):
        if r["round"] == "0":
            continue
        if int(r["rc"]) != 0:
            sys.exit(f"J2 n16 rc != 0: {r}")
        j_n16[(r["case"], r["arm"])].append(
            (float(r["wall_s"]), int(r["maxrss_kb"]) / 1024.0, r["sha16"], parse_rss_anon_mib(r))
        )

    j_n1: dict[tuple[str, str], list[RunRecord]] = defaultdict(list)
    j_n1_timeouts: set[tuple[str, str]] = set()
    for p in j_n1_paths:
        for r in read_tsv(p):
            if r.get("rc") == "124" or (r.get("rc") is None and r.get("user_s") == "124"):
                j_n1_timeouts.add((r["case"], r["arm"]))
                continue
            if int(r["rc"]) != 0:
                sys.exit(f"J2 n1 rc != 0: {r}")
            j_n1[(r["case"], r["arm"])].append(
                (float(r["wall_s"]), int(r["maxrss_kb"]) / 1024.0, r["sha16"], parse_rss_anon_mib(r))
            )

    ext_n0: dict[tuple[str, int, str], list[RunRecord]] = defaultdict(list)
    for r in read_tsv(ext_n0_path):
        if r["round"] == "0":
            continue
        if int(r["rc"]) != 0:
            sys.exit(f"ext n0 rc != 0: {r}")
        ext_n0[(r["case"], int(r["threads"]), r["arm"])].append(
            (float(r["wall_s"]), int(r["maxrss_kb"]) / 1024.0, r["sha16"], parse_rss_anon_mib(r))
        )

    ext_n1: dict[tuple[str, int, str], list[RunRecord]] = defaultdict(list)
    for r in read_tsv(ext_n1_path):
        if r["round"] == "0":
            continue
        if int(r["rc"]) != 0:
            sys.exit(f"ext n1 rc != 0: {r}")
        ext_n1[(r["case"], int(r["threads"]), r["arm"])].append(
            (float(r["wall_s"]), int(r["maxrss_kb"]) / 1024.0, r["sha16"], parse_rss_anon_mib(r))
        )

    # 2. standard_workloads.tsv (12 modes x 2 thread counts x {losat, ncbi})
    std_rows = []
    for group, cid, prog, task, label, qdesc, sdesc, target, opts in MODES_ORDER:
        is_ext = cid in ("psi", "rps", "rpstbn")
        for th in (1, 16):
            if not is_ext:
                lo_list = j_n1.get((cid, "cp"), []) if th == 1 else j_n16.get((cid, "cp"), [])
                nc_list = j_n1.get((cid, "ncbi"), []) if th == 1 else j_n16.get((cid, "ncbi"), [])
                v0_list = j_n1.get((cid, "v010"), []) if th == 1 else j_n16.get((cid, "v010"), [])
            else:
                src = ext_n0 if th == 1 else ext_n1
                lo_list = src.get((cid, th, "cp"), [])
                nc_list = src.get((cid, th, "ncbi"), [])
                v0_list = []

            ncbi_n1_sha = (
                {x[2] for x in ext_n0[(cid, 1, "ncbi")]}
                if cid in ("rps", "rpstbn") and th == 16
                else None
            )
            std_rows.append(
                make_standard_row(
                    group=group,
                    cid=cid,
                    prog=prog,
                    task=task,
                    label=label,
                    th=th,
                    qdesc=qdesc,
                    sdesc=sdesc,
                    target=target,
                    opts=opts,
                    lo_list=lo_list,
                    nc_list=nc_list,
                    v0_list=v0_list,
                    is_ncbi_timeout=(cid, "ncbi") in j_n1_timeouts,
                    subject_vs_db=cid in ("bp", "bpff"),
                    ncbi_n1_sha=ncbi_n1_sha,
                )
            )

            # Optional paired `-db` vs `-db` arm
            if not is_ext and target == "-subject":
                db_cid = f"{cid}_db"
                db_lo = (j_n1 if th == 1 else j_n16).get((db_cid, "cp"), [])
                db_nc = (j_n1 if th == 1 else j_n16).get((db_cid, "ncbi"), [])
                if db_lo and db_nc:
                    std_rows.append(
                        make_standard_row(
                            group=group,
                            cid=db_cid,
                            prog=prog,
                            task=task,
                            label=f"{label} (-db)",
                            th=th,
                            qdesc=qdesc,
                            sdesc=sdesc,
                            target="-db",
                            opts=opts,
                            lo_list=db_lo,
                            nc_list=db_nc,
                            v0_list=[],
                        )
                    )

    std_fields = list(std_rows[0].keys())
    for opt_col in ("losat_rss_anon_mib", "ncbi_rss_anon_mib"):
        if any(opt_col in row for row in std_rows) and opt_col not in std_fields:
            std_fields.append(opt_col)

    # 3. utilities.tsv (makeblastdb + rpsbproc)
    util_rows = []
    for cid, label in (("mkdb_16s", "makeblastdb (nucl: 16S rRNA, 28,036 seqs)"),
                       ("mkdb_sprot", "makeblastdb (prot: Swiss-Prot, 575,748 seqs)"),
                       ("rpsbproc", "rpsbproc (4,300 E. coli proteins XML, CDD v3.22)")):
        lo = ext_n0[(cid, 1, "cp")]
        nc = ext_n0[(cid, 1, "ncbi")]
        lo_w = statistics.median(x[0] for x in lo)
        nc_w = statistics.median(x[0] for x in nc)
        lo_m = statistics.median(x[1] for x in lo)
        nc_m = statistics.median(x[1] for x in nc)
        lo_sha = {x[2] for x in lo}
        nc_sha = {x[2] for x in nc}
        util_rows.append({
            "case": cid,
            "label": label,
            "threads": 1,
            "losat_wall_s": round(lo_w, 3),
            "ncbi_wall_s": round(nc_w, 3),
            "speedup": round(nc_w / lo_w, 4),
            "losat_rss_mib": round(lo_m, 2),
            "ncbi_rss_mib": round(nc_m, 2),
            "mem_reduction": round(nc_m / lo_m, 4),
            "losat_sha16": next(iter(lo_sha)),
            "ncbi_sha16": next(iter(nc_sha)),
            "runs": len(lo),
        })

    # 4. scaling.tsv (2a thread scaling, 2b query size, 2c/2d database size)
    sc_rows = []
    # 2a: bp, bn, tbn @ 1, 2, 4, 8, 16
    for cid in ("bn", "bp", "tbn"):
        for th in (1, 2, 4, 8, 16):
            if th == 1:
                lo = j_n1[(cid, "cp")]
                nc = j_n1[(cid, "ncbi")]
            elif th == 16:
                lo = j_n16[(cid, "cp")]
                nc = j_n16[(cid, "ncbi")]
            else:
                lo = ext_n1[(cid, th, "cp")]
                nc = ext_n1[(cid, th, "ncbi")]
            lo_w = statistics.median(x[0] for x in lo)
            nc_w = statistics.median(x[0] for x in nc)
            lo_m = statistics.median(x[1] for x in lo)
            nc_m = statistics.median(x[1] for x in nc)
            sc_rows.append({
                "panel": "2a_threads",
                "case": cid,
                "x_value": th,
                "x_label": str(th),
                "cache": "warm",
                "losat_wall_s": round(lo_w, 3),
                "ncbi_wall_s": round(nc_w, 3),
                "speedup": round(nc_w / lo_w, 4),
                "losat_rss_mib": round(lo_m, 2),
                "ncbi_rss_mib": round(nc_m, 2),
            })

    # 2b: query size (10k, 100k, 1M) from bench-v020-ext n1.tsv
    for cid, qn in (("m2_1e4", 10000), ("m2_1e5", 100000), ("m2_1e6", 1000000)):
        lo = ext_n1[(cid, 16, "cp")]
        nc = ext_n1[(cid, 16, "ncbi")]
        lo_w = statistics.median(x[0] for x in lo)
        nc_w = statistics.median(x[0] for x in nc)
        lo_m = statistics.median(x[1] for x in lo)
        nc_m = statistics.median(x[1] for x in nc)
        sc_rows.append({
            "panel": "2b_query_size",
            "case": cid,
            "x_value": qn,
            "x_label": f"{qn:,}",
            "cache": "warm",
            "losat_wall_s": round(lo_w, 3),
            "ncbi_wall_s": round(nc_w, 3),
            "speedup": round(nc_w / lo_w, 4),
            "losat_rss_mib": round(lo_m, 2),
            "ncbi_rss_mib": round(nc_m, 2),
        })

    # 2c/2d: DB size (abr_genes, 16s_all, ref_prok_rep warm from bench-v020-ext; rpr cold & core_nt warm/cold from toolmt.json)
    inputs_ext = {r["case"]: r for r in read_tsv(ext_inputs_path)}
    tm = {r["case"]: r for r in json.loads(args.toolmt.read_text(encoding="utf-8"))}
    for cid, inp_key, cache in (("m3_abr", "m3_abr", "warm"), ("m3_16s", "m3_16s", "warm"), ("m3_rpr", "m3_rpr", "warm")):
        lo = ext_n1[(cid, 16, "cp")]
        nc = ext_n1[(cid, 16, "ncbi")]
        db_bytes = int(inputs_ext[inp_key]["disk_bytes"])
        lo_w = statistics.median(x[0] for x in lo)
        nc_w = statistics.median(x[0] for x in nc)
        lo_m = statistics.median(x[1] for x in lo)
        nc_m = statistics.median(x[1] for x in nc)
        sc_rows.append({
            "panel": "2cd_db_size",
            "case": cid,
            "x_value": round(db_bytes / 1e6, 4),  # MB
            "x_label": inputs_ext[inp_key]["file"],
            "cache": cache,
            "losat_wall_s": round(lo_w, 3),
            "ncbi_wall_s": round(nc_w, 3),
            "speedup": round(nc_w / lo_w, 4),
            "losat_rss_mib": round(lo_m, 2),
            "ncbi_rss_mib": round(nc_m, 2),
        })
    for tm_key, inp_key, cache in (("M3_rprc", "m3_rpr", "cold"), ("M3_cnt", "m3_cnt", "warm"), ("M3_cntc", "m3_cnt", "cold")):
        r = tm[tm_key]
        db_bytes = int(inputs_ext[inp_key]["disk_bytes"])
        sc_rows.append({
            "panel": "2cd_db_size",
            "case": tm_key.lower(),
            "x_value": round(db_bytes / 1e6, 4),
            "x_label": inputs_ext[inp_key]["file"],
            "cache": cache,
            "losat_wall_s": round(float(r["losat_s"]), 3),
            "ncbi_wall_s": round(float(r["ncbi_s"]), 3),
            "speedup": round(float(r["x"]), 4),
            "losat_rss_mib": round(float(r["losat_rss_mib"]), 2),
            "ncbi_rss_mib": round(float(r["ncbi_rss_mib"]), 2),
        })

    # 5. webassembly.tsv (strictly paired 7f7db49c runs; tblastx BLAST+ left unmeasured)
    wasm_rows = build_webassembly(args.wasm_csv, None)

    # 6. Combine inputs & bp_db parity metadata
    in_rows = []
    for r in read_tsv(j_inputs_path):
        in_rows.append({"source": "j2", **r})
    for r in read_tsv(ext_inputs_path):
        in_rows.append({"source": "bench_v020_ext", **r})

    bp_par = read_tsv(ext_bp_path) if ext_bp_path.is_file() else []
    toolmt_builds = sorted({str(r.get("build", "")) for r in json.loads(args.toolmt.read_text(encoding="utf-8"))} - {""})
    wasm_commits = sorted({str(r["commit"]) for r in wasm_rows})
    dt_builds = {}
    for tag in ("dt1", "dt13", "dt15"):
        sp = args.dt_dir / f"{tag}_sha.txt"
        dt_builds[tag] = sp.read_text(encoding="utf-8").strip() if sp.is_file() else ""
    meta: dict[str, object] = {
        "losat_version": args.losat_version,
        "losat_commit": args.losat_commit,
        "ncbi_version": "2.17.0",
        "rpsbproc_version": "0.5.1",
        "wasm_commit": ",".join(wasm_commits),
        "wasm_derived": not args.wasm_paired,
    }
    if not args.wasm_paired:
        meta["wasm_derived_from"] = wasm_derived_from
    meta.update({
        "toolmt_commit": ",".join(toolmt_builds),
        "dt_commits": dt_builds,
        "hardware": "Dual Intel(R) Xeon(R) Gold 6338 CPU @ 2.00GHz (2 NUMA nodes, 64 physical cores total)",
        "bp_db_parity": bp_par,
        "input_checksums": input_checksums,
    })

    # Write all snapshot outputs only after every input and row check has passed
    os.makedirs(args.out, exist_ok=True)
    write_tsv(
        args.out / "downstream_tools.tsv",
        ["program", "tool", "threads_label", "threads_num", "cases_count",
         "ncbi_mt_s", "losat_mt_s", "mt_speedup", "ncbi_n1_s", "losat_n1_s",
         "n1_speedup", "same_output", "note"],
        dt_rows,
    )
    write_tsv(
        args.out / "standard_workloads.tsv",
        std_fields,
        std_rows,
    )
    write_tsv(args.out / "utilities.tsv", list(util_rows[0].keys()), util_rows)
    write_tsv(args.out / "scaling.tsv", list(sc_rows[0].keys()), sc_rows)
    write_tsv(
        args.out / "webassembly.tsv",
        list(wasm_rows[0].keys()),
        wasm_rows,
    )
    write_tsv(args.out / "inputs.tsv", ["source", "case", "role", "file", "sequences", "residues", "disk_bytes", "sha256", "db_name"],
              [{k: r.get(k, "") for k in ["source", "case", "role", "file", "sequences", "residues", "disk_bytes", "sha256", "db_name"]} for r in in_rows])
    (args.out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print("full snapshot written:", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
