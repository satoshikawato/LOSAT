#!/usr/bin/env python3
"""Build frozen benchmark snapshot tables in benchmarks/v0.2.0/snapshot/ from measured TSV/JSON/log outputs.

Every input path is an argument (defaults: the v0.2.0 release-candidate and comparative runs under
LOSAT_BASELINES_DIR). metadata.json records the sha256 of every input file read, paired `-db` parity
from <ext-dir>/bp_db_parity.tsv, and the exact LOSAT build commit of every data source:
- --losat-commit for the standard 12-mode grid, utilities, and scaling runs;
- toolmt.json `build` for published multi-thread tools, large-scale workflows, cold-cache, and core_nt;
- <dt-dir>/<tag>_sha.txt for single-thread published tools;
- --uniref50-losat-commit for the UniRef50 LOSAT vs NCBI/DIAMOND/MMseqs2 comparison;
- --silva-losat-commit for the SILVA 144 LOSAT vs NCBI/vsearch comparison;
- webassembly CSV `commit` (with wasm_derived / wasm_derived_from unless --wasm-paired).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

MODES_ORDER = [
    # group, case_id, program, task, label, short_desc, query_desc, subj_desc, target, options
    ("Nucleotide", "mb", "blastn", "megablast", "megablast", "100 16S vs 28k RefSeq",
     "16S rRNA genes, full length", "RefSeq 16S rRNA", "-subject", "-"),
    ("Nucleotide", "dc", "blastn", "dc-megablast", "dc-megablast", "31k cDNA vs 3 Mb chr2L",
     "D. melanogaster cDNA", "D. melanogaster chr. 2L (3 Mb)", "-subject", "-task dc-megablast"),
    ("Nucleotide", "bn", "blastn", "blastn", "blastn", "31k cDNA vs 3 Mb chr2L",
     "D. melanogaster cDNA", "D. melanogaster chr. 2L (3 Mb)", "-subject", "-task blastn"),
    ("Nucleotide", "bns", "blastn", "blastn-short", "blastn-short", "30k 24-mers vs 3 Mb",
     "Synthetic 24-mer reads", "D. melanogaster chr. 2L (3 Mb)", "-subject", "-task blastn-short -evalue 0.001"),
    ("Protein", "bp", "blastp", "blastp", "blastp", "3k human vs yeast",
     "Human proteome (first 3,000)", "S. cerevisiae proteome", "-subject", "-evalue 0.001"),
    ("Protein", "bpff", "blastp", "blastp-fast", "blastp-fast", "24k human vs yeast",
     "Human proteome", "S. cerevisiae proteome", "-subject", "-task blastp-fast -evalue 0.001"),
    ("Translated", "bx", "blastx", "blastx", "blastx", "50 kb DNA vs fly prot",
     "D. melanogaster chr. 2L (50 kb)", "D. melanogaster proteome", "-db", "-evalue 1e-5"),
    ("Translated", "tbn", "tblastn", "tblastn", "tblastn", "1 protein vs GRCh38",
     "Human protein ENSP00000343742.5", "Human genome GRCh38", "-db", "-evalue 1e-5 -max_target_seqs 5"),
    ("Translated", "tbx", "tblastx", "tblastx", "tblastx", "3 Candidatus genomes",
     "Candidatus genomes (3 seqs)", "Candidatus genomes (3 seqs)", "-subject", "-query_gencode 4 -db_gencode 4"),
    ("Profile & Domain", "psi", "psiblast", "psiblast", "psiblast", "20 E. coli vs Swiss-Prot (3x)",
     "E. coli MG1655 proteins (20)", "UniProtKB/Swiss-Prot", "-db", "-num_iterations 3 -evalue 0.001"),
    ("Profile & Domain", "rps", "rpsblast", "rpsblast", "rpsblast", "300 E. coli vs CDD v3.22",
     "E. coli MG1655 proteins (300)", "NCBI CDD v3.22 (29 vol.)", "-db", "-evalue 0.01"),
    ("Profile & Domain", "rpstbn", "rpstblastn", "rpstblastn", "rpstblastn", "100x5 kb vs CDD v3.22",
     "E. coli MG1655 5-kb segments", "NCBI CDD v3.22 (29 vol.)", "-db", "-"),
]

DT_TAGS = ("dt1", "dt13", "dt15")

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

LARGE_WORKFLOW_SPECS = [
    # toolmt case, label, workflow, program, detail
    ("M2_reads1e6s", "1M reads \u00d7 16S (stdin)", "megablast stdin", "blastn", "1M reads (150 nt) vs 16S rRNA, piped stdin"),
    ("M2_reads1e6", "1M reads \u00d7 16S (file)", "megablast file", "blastn", "1M reads (150 nt) vs 16S rRNA, query file"),
    ("M2_amp1e5", "QIIME 2 (100k amplicons)", "QIIME 2", "blastn", "100k amplicons (253 nt) vs 16S rRNA"),
    ("SPB", "Trinotate proteome", "Trinotate", "blastp", "Proteome vs UniProtKB/Swiss-Prot"),
    ("UR50_q100", "UniRef50 (100 proteins)", "UniRef50", "blastp", "100 human proteins vs 38.8M seqs (19 GB)"),
    ("M2_amp1e4", "QIIME 2 (10k amplicons)", "QIIME 2", "blastn", "10k amplicons vs 16S rRNA (NCBI 16t drifts)"),
    ("SPX", "Trinotate transcripts", "Trinotate", "blastx", "Transcripts vs UniProtKB/Swiss-Prot"),
    ("M3_rpr", "BlobToolKit (ref_prok, 27 GB)", "BlobToolKit", "blastn", "210 contigs vs ref_prok_rep_genomes (warm)"),
    ("M3_cnt", "BlobToolKit (core_nt, 319 GB)", "BlobToolKit", "blastn", "210 contigs vs NCBI core_nt (warm)"),
]

RunRecord = tuple[float, float, str]


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


def sanitize_input_key(path: Path, base: Path) -> str:
    resolved_base = base.resolve()
    resolved_path = path.resolve()
    try:
        rel = resolved_path.relative_to(resolved_base)
    except ValueError:
        sys.exit(f"input path is outside baselines root {resolved_base}: {resolved_path}")
    return f"losat-baselines/{rel.as_posix()}"


def collect_input_checksums(paths: list[Path], base: Path) -> dict[str, str]:
    checksums: dict[str, str] = {}
    for p in paths:
        key = sanitize_input_key(p, base)
        if key in checksums:
            sys.exit(f"duplicate sanitized input_checksums key: {key} (from {p})")
        try:
            checksums[key] = sha256_file(p)
        except OSError as exc:
            sys.exit(f"cannot read input file {p}: {exc}")
    return checksums


def row_sha16(row: dict[str, str]) -> str:
    return (row.get("sha16") or "").strip()


def require_single_sha16(values: list[str], context: str) -> str:
    cleaned = {(v or "").strip() for v in values}
    if not cleaned or "" in cleaned:
        sys.exit(f"missing or empty sha16 for {context}")
    if len(cleaned) != 1:
        sys.exit(f"sha16 differs across rounds for {context}: {sorted(cleaned)}")
    return next(iter(cleaned))


def validate_bp_db_parity(rows: list[dict[str, str]], path: Path) -> list[dict[str, str]]:
    if not rows:
        sys.exit(f"empty bp_db_parity provenance file: {path}")
    seen: set[str] = set()
    validated: list[dict[str, str]] = []
    for row in rows:
        if None in row:
            sys.exit(f"malformed row in {path}: {row}")
        cid = (row.get("case") or "").strip()
        if not cid:
            sys.exit(f"missing or empty case in {path}: {row}")
        if cid in seen:
            sys.exit(f"duplicate case {cid!r} in {path}")
        seen.add(cid)
        lo_sha = (row.get("losat_db_sha16") or "").strip()
        nc_sha = (row.get("ncbi_db_sha16") or "").strip()
        if not lo_sha or not nc_sha:
            sys.exit(f"missing or empty db sha16 for {cid} in {path}: {row}")
        equal = (row.get("equal") or "").strip()
        expected_equal = "yes" if lo_sha == nc_sha else "no"
        if equal != expected_equal:
            sys.exit(f"inconsistent equal field for {cid} in {path}: got {equal!r}, expected {expected_equal!r}")
        normalized = dict(row)
        normalized["case"] = cid
        normalized["losat_db_sha16"] = lo_sha
        normalized["ncbi_db_sha16"] = nc_sha
        normalized["equal"] = equal
        validated.append(normalized)
    missing = [cid for cid in ("bp", "bpff") if cid not in seen]
    if missing:
        sys.exit(f"missing required bp_db_parity row(s) {missing} in {path}")
    return validated


def parse_disk_bytes(inputs_ext: dict[str, dict[str, str]], inp_key: str) -> int:
    if inp_key not in inputs_ext:
        sys.exit(f"missing inputs_ext entry for {inp_key}")
    raw = (inputs_ext[inp_key].get("disk_bytes") or "").strip()
    try:
        val = int(raw)
    except ValueError:
        sys.exit(f"invalid disk_bytes for inputs_ext entry {inp_key}: {raw!r}")
    if val <= 0:
        sys.exit(f"non-positive disk_bytes for inputs_ext entry {inp_key}: {raw!r}")
    return val


def make_standard_row(
    *,
    group: str,
    cid: str,
    prog: str,
    task: str,
    label: str,
    short_desc: str,
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
    if not lo_list:
        sys.exit(f"missing LOSAT runs for {cid} t={th}")
    lo_w = statistics.median(x[0] for x in lo_list)
    lo_w_min = min(x[0] for x in lo_list)
    lo_w_max = max(x[0] for x in lo_list)
    lo_m = statistics.median(x[1] for x in lo_list)
    lo_sha16 = require_single_sha16([x[2] for x in lo_list], f"standard {cid} LOSAT t={th}")
    lo_sha = {lo_sha16}

    nc_w = statistics.median(x[0] for x in nc_list) if nc_list else None
    nc_m = statistics.median(x[1] for x in nc_list) if nc_list else None
    nc_sha = {(x[2] or "").strip() for x in nc_list}
    if nc_list and "" in nc_sha:
        sys.exit(f"missing or empty sha16 for standard {cid} NCBI t={th}")
    cleaned_n1_sha: set[str] | None = None
    if ncbi_n1_sha is not None:
        cleaned_n1_sha = {(v or "").strip() for v in ncbi_n1_sha}
        if not cleaned_n1_sha or "" in cleaned_n1_sha:
            sys.exit(f"missing or empty NCBI 1-thread sha16 for {cid} t={th}")
    v0_w = statistics.median(x[0] for x in v0_list) if v0_list else ""
    v0_m = round(statistics.median(x[1] for x in v0_list), 2) if v0_list else ""

    parity = "same"
    ncbi_runs = len(nc_list)
    if not nc_list and th == 1 and is_ncbi_timeout:
        parity = "ncbi_timeout_6000s"
        ncbi_runs = 1
    elif subject_vs_db:
        parity = "subject_vs_db"
    elif cleaned_n1_sha is not None and th == 16 and lo_sha != nc_sha:
        parity = "same_as_ncbi_n1" if lo_sha == cleaned_n1_sha else "diff"
    elif lo_sha != nc_sha:
        parity = "diff"

    row_entry: dict[str, object] = {
        "group": group,
        "case": cid,
        "program": prog,
        "task": task,
        "label": label,
        "short_desc": short_desc,
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
        "losat_sha16": lo_sha16,
        "ncbi_sha16": ",".join(sorted(nc_sha)),
        "parity": parity,
        "query_desc": qdesc,
        "subject_desc": sdesc,
        "target": target,
        "options": opts,
    }
    return row_entry


def build_large_workflows(toolmt_list: list[dict[str, object]]) -> list[dict[str, object]]:
    by_case = {str(r["case"]): r for r in toolmt_list if r.get("mode") == "n16"}
    out: list[dict[str, object]] = []
    for cid, label, workflow, prog, detail in LARGE_WORKFLOW_SPECS:
        r = by_case.get(cid)
        if not r:
            sys.exit(f"missing large workflow {cid} in toolmt.json")
        ncbi_s = float(r["ncbi_s"])
        losat_s = float(r["losat_s"])
        ncbi_rss = float(r["ncbi_rss_mib"])
        losat_rss = float(r["losat_rss_mib"])
        out.append({
            "case": cid,
            "label": label,
            "workflow": workflow,
            "program": prog,
            "threads": int(r["threads"]),
            "detail": detail,
            "ncbi_wall_s": round(ncbi_s, 3),
            "losat_wall_s": round(losat_s, 3),
            "speedup": round(ncbi_s / losat_s, 4),
            "ncbi_rss_mib": round(ncbi_rss, 2),
            "losat_rss_mib": round(losat_rss, 2),
            "mem_reduction": round(ncbi_rss / losat_rss, 4),
            "same_output": "yes" if r["same_output"] == "はい" else "no",
            "ncbi_mt_differs": "yes" if r.get("ncbi_mt_differs") == "はい" else "no",
        })
    return out


def build_downstream_tools(toolmt_list: list[dict[str, object]], dt_dir: Path) -> list[dict[str, object]]:
    by_prog_tool: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for r in toolmt_list:
        cid = str(r["case"])
        if cid.startswith(("M2_", "M3_", "M4_", "UR50", "SP")):
            continue
        if r["tool"] == "ABRicate":
            if r["threads"] != 8:
                continue
        elif r["threads"] != 16:
            continue
        tool_name = "nf-core blast" if r["tool"] == "nf-core" else str(r["tool"])
        by_prog_tool[(str(r["program"]), tool_name)].append(r)

    out: list[dict[str, object]] = []
    for (prog, tool), rs in sorted(by_prog_tool.items()):
        form = rs[0]["form"]
        if form == "プロセス並列":
            th_label = "16 × 1"
            th_num = 16
        else:
            th_num = int(rs[0]["threads"])
            th_label = str(th_num)

        ncbi_mt = sum(float(r["ncbi_s"]) for r in rs)
        losat_mt = sum(float(r["losat_s"]) for r in rs)
        mt_x: object = round(ncbi_mt / losat_mt, 4) if losat_mt > 0 else ""

        ncbi_mt_rss = max(float(r["ncbi_rss_mib"]) for r in rs)
        losat_mt_rss = max(float(r["losat_rss_mib"]) for r in rs)
        mt_mem_red: object = round(ncbi_mt_rss / losat_mt_rss, 4) if losat_mt_rss > 0 else ""

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
            "ncbi_mt_rss_mib": round(ncbi_mt_rss, 2),
            "losat_mt_rss_mib": round(losat_mt_rss, 2),
            "mt_mem_reduction": mt_mem_red,
            "ncbi_n1_s": round(ncbi_n1, 3) if ncbi_n1 is not None else "",
            "losat_n1_s": round(losat_n1, 3) if losat_n1 is not None else "",
            "n1_speedup": n1_x,
            "ncbi_n1_rss_mib": "",
            "losat_n1_rss_mib": "",
            "n1_mem_reduction": "",
            "same_output": "yes" if all(r["same_output"] == "はい" for r in rs) else "no",
            "ncbi_mt_differs": "yes" if any(r.get("ncbi_mt_differs") == "はい" for r in rs) else "no",
            "note": note,
        })

    # Single-thread tools from DT1 / DT13 / DT15
    dt_wall: dict[str, dict[tuple[str, str], list[float]]] = {}
    dt_rss: dict[str, dict[tuple[str, str], list[float]]] = {}
    dt_sha: dict[str, dict[tuple[str, str], list[str]]] = {}
    wanted: dict[str, set[str]] = defaultdict(set)
    for _prog, _tool, _th, tag, case_ids in ST_TOOLS:
        wanted[tag].update(case_ids)
    for tag in DT_TAGS:
        p = dt_dir / f"{tag}.tsv"
        gw: dict[tuple[str, str], list[float]] = defaultdict(list)
        gr: dict[tuple[str, str], list[float]] = defaultdict(list)
        gs: dict[tuple[str, str], list[str]] = defaultdict(list)
        for row in read_tsv(p):
            if row["case"] in wanted.get(tag, set()) and row.get("rc", "0") != "0":
                sys.exit(f"{tag} rc != 0: {row}")
            gw[(row["case"], row["arm"])].append(float(row["wall_s"]))
            if row.get("maxrss_kib"):
                gr[(row["case"], row["arm"])].append(int(row["maxrss_kib"]) / 1024.0)
            gs[(row["case"], row["arm"])].append(row_sha16(row))
        dt_wall[tag] = gw
        dt_rss[tag] = gr
        dt_sha[tag] = gs

    for prog, tool, th_label, tag, case_ids in ST_TOOLS:
        gw, gr, gs = dt_wall[tag], dt_rss[tag], dt_sha[tag]
        nc_tot = 0.0
        lo_tot = 0.0
        nc_peak_rss = 0.0
        lo_peak_rss = 0.0
        shas_match = True
        for cid in case_ids:
            nc_vals = gw.get((cid, "ncbi"))
            lo_vals = gw.get((cid, "losat"))
            if not nc_vals or not lo_vals:
                sys.exit(f"missing {cid} in {tag}")
            nc_sha = require_single_sha16(gs.get((cid, "ncbi"), []), f"{tag} {cid} NCBI")
            lo_sha = require_single_sha16(gs.get((cid, "losat"), []), f"{tag} {cid} LOSAT")
            nc_tot += statistics.median(nc_vals)
            lo_tot += statistics.median(lo_vals)
            if gr.get((cid, "ncbi")):
                nc_peak_rss = max(nc_peak_rss, statistics.median(gr[(cid, "ncbi")]))
            if gr.get((cid, "losat")):
                lo_peak_rss = max(lo_peak_rss, statistics.median(gr[(cid, "losat")]))
            if nc_sha != lo_sha:
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
            "ncbi_mt_rss_mib": "",
            "losat_mt_rss_mib": "",
            "mt_mem_reduction": "",
            "ncbi_n1_s": round(nc_tot, 3),
            "losat_n1_s": round(lo_tot, 3),
            "n1_speedup": round(nc_tot / lo_tot, 4),
            "ncbi_n1_rss_mib": round(nc_peak_rss, 2) if nc_peak_rss > 0 else "",
            "losat_n1_rss_mib": round(lo_peak_rss, 2) if lo_peak_rss > 0 else "",
            "n1_mem_reduction": round(nc_peak_rss / lo_peak_rss, 4) if lo_peak_rss > 0 else "",
            "same_output": "yes" if shas_match else "no",
            "ncbi_mt_differs": "no",
            "note": f"1t_only_{tag}",
        })
    return out


def _parse_recall_file(path: Path) -> dict[str, tuple[float, float]]:
    """Parse recall_q100.txt / recall_q1000.txt -> {arm: (uncapped_recall, best_hit_recall)}."""
    lines = [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if len(lines) < 3:
        sys.exit(f"malformed recall file: {path}")
    hdr = lines[1].split("\t")
    try:
        idx_arm = hdr.index("arm")
        idx_unc = hdr.index("uncapped")
        idx_best = hdr.index("best-hit")
    except ValueError as exc:
        sys.exit(f"missing column in recall file {path}: {exc}")
    out: dict[str, tuple[float, float]] = {}
    for ln in lines[2:]:
        parts = ln.split("\t")
        out[parts[idx_arm].strip()] = (float(parts[idx_unc]), float(parts[idx_best]))
    return out


def build_vs_heuristics_uniref50(
    ur50_main: Path,
    bp16_summary: Path,
    mmseqs_noidx: Path,
    mmseqs_idx_log: Path,
    recall_q100: Path,
    recall_q1000: Path,
    losat_commit: str,
) -> list[dict[str, object]]:
    """Parse raw UniRef50 benchmark logs (32 threads, 38.8M seqs) into a normalized comparison table."""
    rec = {"q100": _parse_recall_file(recall_q100), "q1000": _parse_recall_file(recall_q1000)}

    def _ur50_head(parts: list[str], src: Path, raw_line: str) -> tuple[str, str, int]:
        qs, ar, rd, code = parts[0], parts[1], int(parts[2]), int(parts[3])
        if code != 0:
            sys.exit(f"non-zero rc in {src}: {raw_line}")
        return qs, ar, rd

    # 1. summary_main.tsv (no header: qset, arm, round, rc, wall_s, user_s, sys_s, maxrss_kb, anon_kb, file_kb, rows, sha16, ...)
    # Uniform max RSS definition across all UniRef50 inputs: raw kernel maxrss_kb (column 7) in GiB.
    # Report convention: round 2 for q100 (round 1 NCBI had concurrent load), round 1 for q1000.
    main_by: dict[tuple[str, str], dict[int, tuple[float, float, float, str]]] = defaultdict(dict)
    for ln in ur50_main.read_text(encoding="utf-8").splitlines():
        if not ln.strip():
            continue
        p = ln.split("\t")
        qset, arm, rnd = _ur50_head(p, ur50_main, ln)
        main_by[(qset, arm)][rnd] = (
            float(p[4]),
            int(p[7]) / (1024.0 * 1024.0),  # GiB max RSS (kernel maxrss_kb)
            int(p[8]) / (1024.0 * 1024.0),  # GiB RssAnon
            p[11].strip(),
        )

    # 2. bp16_summary.tsv (updated LOSAT dev 995281681: same schema p[7]=maxrss_kb, p[8]=anon_kb)
    bp16_by: dict[tuple[str, str], dict[int, tuple[float, float, float, str]]] = defaultdict(dict)
    for ln in bp16_summary.read_text(encoding="utf-8").splitlines():
        if not ln.strip():
            continue
        p = ln.split("\t")
        qset, arm, rnd = _ur50_head(p, bp16_summary, ln)
        bp16_by[(qset, arm)][rnd] = (
            float(p[4]),
            int(p[7]) / (1024.0 * 1024.0),
            int(p[8]) / (1024.0 * 1024.0),
            p[11].strip(),
        )

    # 3. MMseqs2 unindexed (summary.tsv)
    mm_noidx: dict[tuple[str, str], dict[int, tuple[float, float, float]]] = defaultdict(dict)
    for ln in mmseqs_noidx.read_text(encoding="utf-8").splitlines():
        if not ln.strip():
            continue
        p = ln.split("\t")
        qset, arm, rnd = _ur50_head(p, mmseqs_noidx, ln)
        mm_noidx[(qset, arm)][rnd] = (
            float(p[4]),
            int(p[7]) / (1024.0 * 1024.0),
            int(p[8]) / (1024.0 * 1024.0),
        )

    # 4. MMseqs2 indexed (b3i.txt)
    mm_idx: dict[tuple[str, str], dict[int, tuple[float, float, float]]] = defaultdict(dict)
    pat = re.compile(
        r"^\S+\s+(q\d+)\s+(\S+)\s+r(\d+)\s+rc=(\d+)\s+index\s+([\d.]+)\s+[\d.]+\s+[\d.]+\s+(\d+)\s+(\d+)\s+(\d+)"
    )
    for ln in mmseqs_idx_log.read_text(encoding="utf-8").splitlines():
        m = pat.match(ln.strip())
        if not m:
            continue
        qset, arm, rnd, rc = m.group(1), m.group(2), int(m.group(3)), int(m.group(4))
        if rc != 0:
            sys.exit(f"non-zero rc in {mmseqs_idx_log}: {ln}")
        mm_idx[(qset, arm)][rnd] = (
            float(m.group(5)),
            int(m.group(6)) / (1024.0 * 1024.0),
            int(m.group(7)) / (1024.0 * 1024.0),
        )

    out: list[dict[str, object]] = []
    for qset in ("q100", "q1000"):
        rnd = 2 if qset == "q100" else 1
        # Check LOSAT output SHA matches NCBI exact SHA
        ncbi_def_sha = main_by[(qset, "ncbi_def")][rnd][3]
        ncbi_fast_sha = main_by[(qset, "ncbi_fast")][rnd][3]
        lo_bp_sha = bp16_by[(qset, "losat_bp")][rnd][3]
        lo_fast_sha = bp16_by[(qset, "losat_fast")][rnd][3]
        if lo_bp_sha != ncbi_def_sha or lo_fast_sha != ncbi_fast_sha:
            sys.exit(f"UniRef50 parity mismatch for {qset}: {lo_bp_sha} vs {ncbi_def_sha}, {lo_fast_sha} vs {ncbi_fast_sha}")

        entries = [
            ("losat_bp", "LOSAT", "blastp", "no", bp16_by[(qset, "losat_bp")][rnd][:3], (1.0, 1.0), "yes", f"LOSAT {losat_commit}"),
            ("losat_fast", "LOSAT", "blastp-fast", "no", bp16_by[(qset, "losat_fast")][rnd][:3], rec[qset]["ncbi_fast"], "yes_fast", f"LOSAT {losat_commit}"),
            ("ncbi_def", "NCBI BLAST+", "blastp", "no", main_by[(qset, "ncbi_def")][rnd][:3], (1.0, 1.0), "baseline", "NCBI BLAST+ 2.17.0"),
            ("ncbi_fast", "NCBI BLAST+", "blastp-fast", "no", main_by[(qset, "ncbi_fast")][rnd][:3], rec[qset]["ncbi_fast"], "baseline_fast", "NCBI BLAST+ 2.17.0"),
            ("dmnd_def", "DIAMOND", "default", "no", main_by[(qset, "dmnd_def")][rnd][:3], rec[qset]["dmnd_def"], "no", "DIAMOND 2.2.8"),
            ("dmnd_ms", "DIAMOND", "--more-sensitive", "no", main_by[(qset, "dmnd_ms")][rnd][:3], rec[qset]["dmnd_ms"], "no", "DIAMOND 2.2.8"),
            ("dmnd_us", "DIAMOND", "--ultra-sensitive", "no", main_by[(qset, "dmnd_us")][rnd][:3], rec[qset]["dmnd_us"], "no", "DIAMOND 2.2.8"),
            ("mmseqs_s57_noidx", "MMseqs2", "-s 5.7 (no index)", "no", mm_noidx[(qset, "s57_bioconda")][rnd], rec[qset]["s57"], "no", "MMseqs2 18.8cc5c"),
            ("mmseqs_s75_noidx", "MMseqs2", "-s 7.5 (no index)", "no", mm_noidx[(qset, "s75_bioconda")][rnd], rec[qset]["s75"], "no", "MMseqs2 18.8cc5c"),
            ("mmseqs_s57_idx", "MMseqs2", "-s 5.7 (101 GiB idx)", "yes", mm_idx[(qset, "s57_bioconda")][rnd], rec[qset]["s57"], "no", "MMseqs2 18.8cc5c"),
            ("mmseqs_s75_idx", "MMseqs2", "-s 7.5 (101 GiB idx)", "yes", mm_idx[(qset, "s75_bioconda")][rnd], rec[qset]["s75"], "no", "MMseqs2 18.8cc5c"),
        ]
        for arm_id, tool, mode, indexed, (wall, rss_g, anon_g), (r_unc, r_best), parity, bnote in entries:
            out.append({
                "query_set": qset,
                "arm_id": arm_id,
                "tool": tool,
                "mode": mode,
                "indexed": indexed,
                "threads": 32,
                "wall_s": round(wall, 2),
                "max_rss_gib": round(rss_g, 3),
                "rss_anon_gib": round(anon_g, 3),
                "recall_uncapped": round(r_unc, 4),
                "recall_best_hit": round(r_best, 4),
                "exact_ncbi_parity": parity,
                "build_note": bnote,
            })
    return out


def _parse_silva_state(silva_state: Path) -> tuple[int, str, int, int]:
    """Extract amp1e4 NCBI 1-thread reference (rows, sha16) and vsearch/NCBI pair overlap from STATE.md."""
    text = silva_state.read_text(encoding="utf-8")
    m_ref = re.search(r"NCBI amp1e4 n1\s+([\d,]+)\s+rows,\s+([0-9a-f]{16})", text)
    if not m_ref:
        sys.exit(f"missing NCBI amp1e4 n1 reference row count/sha16 in {silva_state}")
    m_ov = re.search(r"vsearch vs NCBI pairs on amp1e3[^:]*:\s*([\d,]+)\s+of NCBI's\s+([\d,]+)", text)
    if not m_ov:
        sys.exit(f"missing vsearch vs NCBI amp1e3 overlap counts in {silva_state}")
    n1_rows = int(m_ref.group(1).replace(",", ""))
    n1_sha = m_ref.group(2)
    ov_both = int(m_ov.group(1).replace(",", ""))
    ov_ncbi = int(m_ov.group(2).replace(",", ""))
    if n1_rows <= 0 or ov_ncbi <= 0 or ov_both <= 0 or ov_both > ov_ncbi:
        sys.exit(f"invalid parsed SILVA provenance values in {silva_state}: {(n1_rows, n1_sha, ov_both, ov_ncbi)}")
    return n1_rows, n1_sha, ov_both, ov_ncbi


def build_vs_heuristics_silva(silva_tsv: Path, silva_state: Path, losat_commit: str) -> list[dict[str, object]]:
    """Parse bench-silva-20261031/logs_ws/summary_timed.tsv and check hashes against NCBI / STATE.md."""
    amp1e4_n1_rows, amp1e4_n1_sha, ov_both, ov_ncbi = _parse_silva_state(silva_state)
    ov_pct = round(100.0 * ov_both / ov_ncbi)

    grouped: dict[tuple[str, str], list[tuple[float, float, float, int, str]]] = defaultdict(list)
    for ln in silva_tsv.read_text(encoding="utf-8").splitlines():
        if not ln.strip():
            continue
        p = ln.split("\t")
        cid, arm, rc = p[0], p[1], int(p[3])
        if rc != 0:
            sys.exit(f"non-zero rc in {silva_tsv}: {ln}")
        grouped[(cid, arm)].append((
            float(p[4]),
            int(p[7]) / 1024.0,   # max RSS MiB
            int(p[8]) / 1024.0,   # RssAnon MiB
            int(p[10]),           # rows
            p[11].strip(),        # sha16
        ))

    order = [
        ("reads1e6", "1M shotgun reads (megablast)", ("losat", "ncbi")),
        ("amp1e3", "1,000 amplicons (QIIME 2)", ("losat", "ncbi", "vsearch")),
        ("amp1e4", "10,000 amplicons (QIIME 2)", ("losat", "ncbi")),
    ]
    tool_names = {"losat": "LOSAT", "ncbi": "NCBI BLAST+", "vsearch": "vsearch"}
    builds = {"losat": f"LOSAT {losat_commit}", "ncbi": "NCBI BLAST+ 2.17.0", "vsearch": "vsearch 2.32.0"}
    out: list[dict[str, object]] = []
    for cid, label, arms in order:
        parsed_arms: dict[str, tuple[float, float, float, int, str, int]] = {}
        for arm in arms:
            recs = grouped.get((cid, arm))
            if not recs:
                sys.exit(f"missing SILVA case {cid} arm {arm} in {silva_tsv}")
            rows_set = {x[3] for x in recs}
            if len(rows_set) != 1:
                sys.exit(f"SILVA {cid} {arm} row count varies across rounds: {sorted(rows_set)}")
            parsed_arms[arm] = (
                statistics.mean(x[0] for x in recs),
                statistics.mean(x[1] for x in recs),
                statistics.mean(x[2] for x in recs),
                recs[0][3],
                require_single_sha16([x[4] for x in recs], f"SILVA {cid} {arm}"),
                len(recs),
            )

        lo_rows, lo_sha = parsed_arms["losat"][3], parsed_arms["losat"][4]
        nc_rows, nc_sha = parsed_arms["ncbi"][3], parsed_arms["ncbi"][4]
        if cid in ("reads1e6", "amp1e3"):
            if lo_sha != nc_sha or lo_rows != nc_rows:
                sys.exit(f"SILVA {cid} parity mismatch: LOSAT {lo_sha} ({lo_rows}) vs NCBI {nc_sha} ({nc_rows})")
            lo_parity, nc_parity = "same", "baseline"
        else:
            if lo_sha != amp1e4_n1_sha or lo_rows != amp1e4_n1_rows:
                sys.exit(
                    f"SILVA amp1e4 LOSAT mismatch vs NCBI 1-thread reference: "
                    f"LOSAT {lo_sha} ({lo_rows}) vs NCBI n1 {amp1e4_n1_sha} ({amp1e4_n1_rows})"
                )
            lo_parity = "same_as_ncbi_n1"
            if nc_sha == amp1e4_n1_sha:
                nc_parity = "baseline"
            else:
                drift = nc_rows - amp1e4_n1_rows
                nc_parity = f"mt_drifts_{abs(drift)}_rows" if drift != 0 else "mt_drifts_sha"

        for arm in arms:
            wall, rss, anon, rows_n, sha16, n_runs = parsed_arms[arm]
            if arm == "losat":
                parity = lo_parity
            elif arm == "ncbi":
                parity = nc_parity
            else:
                parity = f"heuristic_{ov_pct}pct_overlap"
            out.append({
                "case": cid,
                "label": label,
                "arm": arm,
                "tool": tool_names[arm],
                "threads": 16,
                "runs": n_runs,
                "wall_s": round(wall, 2),
                "max_rss_mib": round(rss, 2),
                "rss_anon_mib": round(anon, 2),
                "rows": rows_n,
                "sha16": sha16,
                "parity": parity,
                "build_note": builds[arm],
            })
    return out


def build_webassembly(wasm_csv: Path) -> list[dict[str, object]]:
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
                "wasm_pct_of_native": round(100.0 * nat / was, 1),
                "native_mib": r["native_mib"],
                "wasm_mib": r["wasm_mib"],
                "ncbi_mib": r["ncbi_mib"] or "",
                "note": f"commit_{r['commit']}",
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
    ap.add_argument("--ur50-main", type=Path, default=base / "bench-ur50-20261014/logs/summary_main.tsv")
    ap.add_argument("--ur50-losat", type=Path, default=base / "bench-buildpolicy-20261014/logs_ws/bp16_summary.tsv")
    ap.add_argument("--ur50-mmseqs-noidx", type=Path, default=base / "bench-buildpolicy-20261014/logs_ws/summary.tsv")
    ap.add_argument("--ur50-mmseqs-idx", type=Path, default=base / "bench-buildpolicy-20261014/logs_ws/b3i.txt")
    ap.add_argument("--ur50-recall-q100", type=Path, default=base / "bench-buildpolicy-20261014/logs_ws/recall_q100.txt")
    ap.add_argument("--ur50-recall-q1000", type=Path, default=base / "bench-buildpolicy-20261014/logs_ws/recall_q1000.txt")
    ap.add_argument("--silva-tsv", type=Path, default=base / "bench-silva-20261031/logs_ws/summary_timed.tsv")
    ap.add_argument("--silva-state", type=Path, default=base / "bench-silva-20261031/STATE.md")
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent.parent / "snapshot")
    ap.add_argument("--losat-commit", default="93ca2f652617acac09d416c28c293cc8a1a501d1")
    ap.add_argument("--uniref50-losat-commit", default="995281681")
    ap.add_argument("--silva-losat-commit", default="dde737f18")
    ap.add_argument("--losat-version", default="0.2.0")
    ap.add_argument("--wasm-paired", action="store_true")
    ap.add_argument("--wasm-derived-from", default="6a0a8719")
    args = ap.parse_args()

    wasm_derived_from = args.wasm_derived_from.strip()
    if not args.wasm_paired and not wasm_derived_from:
        sys.exit("--wasm-derived-from must be non-empty unless --wasm-paired is set")

    dt_source_paths = [args.toolmt, *(args.dt_dir / f"{tag}.tsv" for tag in DT_TAGS)]
    j_n16_path = args.j2_dir / "n16/runs.tsv"
    j_n1_paths = sorted((args.j2_dir / "n1").glob("perfstat_*.tsv"))
    j_inputs_path = args.j2_dir / "inputs.tsv"
    ext_n0_path = args.ext_dir / "n0.tsv"
    ext_n1_path = args.ext_dir / "n1.tsv"
    ext_inputs_path = args.ext_dir / "inputs_ext.tsv"
    ext_bp_path = args.ext_dir / "bp_db_parity.tsv"
    dt_sha_paths = [args.dt_dir / f"{tag}_sha.txt" for tag in DT_TAGS]
    heur_paths = [
        args.ur50_main,
        args.ur50_losat,
        args.ur50_mmseqs_noidx,
        args.ur50_mmseqs_idx,
        args.ur50_recall_q100,
        args.ur50_recall_q1000,
        args.silva_tsv,
        args.silva_state,
    ]

    for req in (
        j_n16_path, j_inputs_path, ext_n0_path, ext_n1_path, ext_inputs_path, ext_bp_path,
        *dt_source_paths, *dt_sha_paths, args.wasm_csv, *heur_paths,
    ):
        if not req.is_file():
            sys.exit(f"missing completed input file: {req}")
    if not j_n1_paths:
        sys.exit(f"missing completed input files: {args.j2_dir / 'n1'}/perfstat_*.tsv")

    input_files = [
        j_n16_path, *j_n1_paths, j_inputs_path, ext_n0_path, ext_n1_path, ext_inputs_path,
        ext_bp_path, args.wasm_csv, *dt_source_paths, *dt_sha_paths, *heur_paths,
    ]
    input_checksums = collect_input_checksums(input_files, base)

    toolmt_list = json.loads(args.toolmt.read_text(encoding="utf-8"))
    dt_rows = build_downstream_tools(toolmt_list, args.dt_dir)
    large_rows = build_large_workflows(toolmt_list)
    ur50_rows = build_vs_heuristics_uniref50(
        args.ur50_main, args.ur50_losat, args.ur50_mmseqs_noidx, args.ur50_mmseqs_idx,
        args.ur50_recall_q100, args.ur50_recall_q1000, args.uniref50_losat_commit,
    )
    silva_rows = build_vs_heuristics_silva(args.silva_tsv, args.silva_state, args.silva_losat_commit)

    def _to_run_record(row: dict[str, str]) -> RunRecord:
        return (float(row["wall_s"]), int(row["maxrss_kb"]) / 1024.0, row_sha16(row))

    def _load_ext_bucket(path: Path, tag: str) -> dict[tuple[str, int, str], list[RunRecord]]:
        bucket: dict[tuple[str, int, str], list[RunRecord]] = defaultdict(list)
        for row in read_tsv(path):
            if row["round"] == "0":
                continue
            if int(row["rc"]) != 0:
                sys.exit(f"ext {tag} rc != 0: {row}")
            bucket[(row["case"], int(row["threads"]), row["arm"])].append(_to_run_record(row))
        return bucket

    j_n16: dict[tuple[str, str], list[RunRecord]] = defaultdict(list)
    for r in read_tsv(j_n16_path):
        if r["round"] == "0":
            continue
        if int(r["rc"]) != 0:
            sys.exit(f"J2 n16 rc != 0: {r}")
        j_n16[(r["case"], r["arm"])].append(_to_run_record(r))

    j_n1: dict[tuple[str, str], list[RunRecord]] = defaultdict(list)
    j_n1_timeouts: set[tuple[str, str]] = set()
    for p in j_n1_paths:
        for r in read_tsv(p):
            if r.get("rc") == "124" or (r.get("rc") is None and r.get("user_s") == "124"):
                j_n1_timeouts.add((r["case"], r["arm"]))
                continue
            if int(r["rc"]) != 0:
                sys.exit(f"J2 n1 rc != 0: {r}")
            j_n1[(r["case"], r["arm"])].append(_to_run_record(r))

    ext_n0 = _load_ext_bucket(ext_n0_path, "n0")
    ext_n1 = _load_ext_bucket(ext_n1_path, "n1")

    std_rows = []
    for group, cid, prog, task, label, short_desc, qdesc, sdesc, target, opts in MODES_ORDER:
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
                    group=group, cid=cid, prog=prog, task=task, label=label, short_desc=short_desc,
                    th=th, qdesc=qdesc, sdesc=sdesc, target=target, opts=opts,
                    lo_list=lo_list, nc_list=nc_list, v0_list=v0_list,
                    is_ncbi_timeout=(cid, "ncbi") in j_n1_timeouts,
                    subject_vs_db=cid in ("bp", "bpff"),
                    ncbi_n1_sha=ncbi_n1_sha,
                )
            )
            if not is_ext and target == "-subject":
                db_cid = f"{cid}_db"
                db_lo = (j_n1 if th == 1 else j_n16).get((db_cid, "cp"), [])
                db_nc = (j_n1 if th == 1 else j_n16).get((db_cid, "ncbi"), [])
                if db_lo and db_nc:
                    std_rows.append(
                        make_standard_row(
                            group=group, cid=db_cid, prog=prog, task=task, label=f"{label} (-db)",
                            short_desc=short_desc, th=th, qdesc=qdesc, sdesc=sdesc, target="-db",
                            opts=opts, lo_list=db_lo, nc_list=db_nc, v0_list=[],
                        )
                    )

    std_fields = list(std_rows[0].keys())

    util_rows = []
    for cid, label, short_desc in (
        ("mkdb_16s", "makeblastdb (nucl)", "28k 16S rRNA seqs"),
        ("mkdb_sprot", "makeblastdb (prot)", "576k Swiss-Prot seqs"),
        ("rpsbproc", "rpsbproc", "4,300 E. coli XML vs CDD"),
    ):
        lo = ext_n0[(cid, 1, "cp")]
        nc = ext_n0[(cid, 1, "ncbi")]
        if not lo or not nc:
            sys.exit(f"missing utility runs for {cid} t=1")
        lo_w = statistics.median(x[0] for x in lo)
        nc_w = statistics.median(x[0] for x in nc)
        lo_m = statistics.median(x[1] for x in lo)
        nc_m = statistics.median(x[1] for x in nc)
        lo_sha16 = require_single_sha16([x[2] for x in lo], f"utility {cid} LOSAT t=1")
        nc_sha16 = require_single_sha16([x[2] for x in nc], f"utility {cid} NCBI t=1")
        util_rows.append({
            "case": cid,
            "label": label,
            "short_desc": short_desc,
            "threads": 1,
            "losat_wall_s": round(lo_w, 3),
            "ncbi_wall_s": round(nc_w, 3),
            "speedup": round(nc_w / lo_w, 4),
            "losat_rss_mib": round(lo_m, 2),
            "ncbi_rss_mib": round(nc_m, 2),
            "mem_reduction": round(nc_m / lo_m, 4),
            "losat_sha16": lo_sha16,
            "ncbi_sha16": nc_sha16,
            "runs": len(lo),
        })

    def _scaling_row(
        panel: str, cid: str, x_val: object, x_lbl: str, cache: str,
        lo_w: float, nc_w: float, sp: float, lo_m: float, nc_m: float,
    ) -> dict[str, object]:
        return {
            "panel": panel, "case": cid, "x_value": x_val, "x_label": x_lbl, "cache": cache,
            "losat_wall_s": round(lo_w, 3), "ncbi_wall_s": round(nc_w, 3),
            "speedup": round(sp, 4),
            "losat_rss_mib": round(lo_m, 2), "ncbi_rss_mib": round(nc_m, 2),
        }

    def _median_scaling_row(
        panel: str, cid: str, x_val: object, x_lbl: str, lo: list[RunRecord], nc: list[RunRecord],
    ) -> dict[str, object]:
        if not lo or not nc:
            sys.exit(f"missing {panel} scaling runs for {cid}")
        lo_w, nc_w = statistics.median(x[0] for x in lo), statistics.median(x[0] for x in nc)
        lo_m, nc_m = statistics.median(x[1] for x in lo), statistics.median(x[1] for x in nc)
        return _scaling_row(panel, cid, x_val, x_lbl, "warm", lo_w, nc_w, nc_w / lo_w, lo_m, nc_m)

    sc_rows = []
    for cid in ("bn", "bp", "tbn"):
        for th in (1, 2, 4, 8, 16):
            lo = j_n1[(cid, "cp")] if th == 1 else (j_n16[(cid, "cp")] if th == 16 else ext_n1[(cid, th, "cp")])
            nc = j_n1[(cid, "ncbi")] if th == 1 else (j_n16[(cid, "ncbi")] if th == 16 else ext_n1[(cid, th, "ncbi")])
            sc_rows.append(_median_scaling_row("threads", cid, th, str(th), lo, nc))
    for cid, qn in (("m2_1e4", 10000), ("m2_1e5", 100000), ("m2_1e6", 1000000)):
        sc_rows.append(
            _median_scaling_row("query_size", cid, qn, f"{qn:,}", ext_n1[(cid, 16, "cp")], ext_n1[(cid, 16, "ncbi")])
        )

    inputs_ext = {r["case"]: r for r in read_tsv(ext_inputs_path)}
    tm_by_case = {str(r["case"]): r for r in toolmt_list}

    def _db_mb(inp_key: str) -> float:
        return round(parse_disk_bytes(inputs_ext, inp_key) / 1e6, 4)

    for cid, inp_key in (("m3_abr", "m3_abr"), ("m3_16s", "m3_16s"), ("m3_rpr", "m3_rpr")):
        sc_rows.append(
            _median_scaling_row(
                "db_size", cid, _db_mb(inp_key), inputs_ext[inp_key]["file"],
                ext_n1[(cid, 16, "cp")], ext_n1[(cid, 16, "ncbi")],
            )
        )
    for tm_key, inp_key, cache in (("M3_rprc", "m3_rpr", "cold"), ("M3_cnt", "m3_cnt", "warm"), ("M3_cntc", "m3_cnt", "cold")):
        if tm_key not in tm_by_case:
            sys.exit(f"missing toolmt entry for {tm_key}")
        r = tm_by_case[tm_key]
        sc_rows.append(
            _scaling_row(
                "db_size", tm_key.lower(), _db_mb(inp_key), inputs_ext[inp_key]["file"], cache,
                float(r["losat_s"]), float(r["ncbi_s"]), float(r["x"]),
                float(r["losat_rss_mib"]), float(r["ncbi_rss_mib"]),
            )
        )

    wasm_rows = build_webassembly(args.wasm_csv)
    in_rows = [{"source": "j2", **r} for r in read_tsv(j_inputs_path)] + [
        {"source": "bench_v020_ext", **r} for r in read_tsv(ext_inputs_path)
    ]
    bp_par = validate_bp_db_parity(read_tsv(ext_bp_path), ext_bp_path)
    toolmt_builds = sorted({str(r.get("build", "")) for r in toolmt_list} - {""})
    wasm_commits = sorted({str(r["commit"]) for r in wasm_rows})
    dt_builds = {}
    for tag, sp in zip(DT_TAGS, dt_sha_paths):
        commit = sp.read_text(encoding="utf-8").strip()
        if not commit:
            sys.exit(f"empty commit in provenance file: {sp}")
        dt_builds[tag] = commit

    meta: dict[str, object] = {
        "losat_version": args.losat_version,
        "losat_commit": args.losat_commit,
        "ncbi_version": "2.17.0",
        "rpsbproc_version": "0.5.1",
        "diamond_version": "2.2.8",
        "mmseqs2_version": "18.8cc5c",
        "vsearch_version": "2.32.0",
        "uniref50_losat_commit": args.uniref50_losat_commit,
        "silva_losat_commit": args.silva_losat_commit,
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

    for name, rows in (
        ("standard_workloads", std_rows),
        ("utilities", util_rows),
        ("scaling", sc_rows),
        ("large_workflows", large_rows),
        ("downstream_tools", dt_rows),
        ("vs_heuristics_uniref50", ur50_rows),
        ("vs_heuristics_silva", silva_rows),
        ("webassembly", wasm_rows),
    ):
        if not rows:
            sys.exit(f"no rows for {name}.tsv")

    os.makedirs(args.out, exist_ok=True)
    write_tsv(args.out / "standard_workloads.tsv", std_fields, std_rows)
    write_tsv(args.out / "utilities.tsv", list(util_rows[0].keys()), util_rows)
    write_tsv(args.out / "scaling.tsv", list(sc_rows[0].keys()), sc_rows)
    write_tsv(args.out / "large_workflows.tsv", list(large_rows[0].keys()), large_rows)
    write_tsv(args.out / "downstream_tools.tsv", list(dt_rows[0].keys()), dt_rows)
    write_tsv(args.out / "vs_heuristics_uniref50.tsv", list(ur50_rows[0].keys()), ur50_rows)
    write_tsv(args.out / "vs_heuristics_silva.tsv", list(silva_rows[0].keys()), silva_rows)
    write_tsv(args.out / "webassembly.tsv", list(wasm_rows[0].keys()), wasm_rows)
    in_fields = ["source", "case", "role", "file", "sequences", "residues", "disk_bytes", "sha256", "db_name"]
    write_tsv(args.out / "inputs.tsv", in_fields, in_rows)
    (args.out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")

    print(
        f"Wrote snapshot to {args.out}: std={len(std_rows)} util={len(util_rows)} "
        f"sc={len(sc_rows)} large={len(large_rows)} dt={len(dt_rows)} "
        f"ur50={len(ur50_rows)} silva={len(silva_rows)} wasm={len(wasm_rows)} in={len(in_rows)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
