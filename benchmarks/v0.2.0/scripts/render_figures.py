#!/usr/bin/env python3
"""Deterministically render Figures 1-4 (light and dark PNGs) for LOSAT v0.2.0 benchmarks.

Writes benchmarks/v0.2.0/figures/ and the same files to docs/manual/src/images/ (the manual's copy;
--manual-images '' skips it). Each figure names its LOSAT builds, the baseline and the machine from
snapshot/metadata.json in a footer line.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, NullFormatter
from matplotlib.transforms import blended_transform_factory

LIGHT = dict(
    surface="#fcfcfb",
    ink="#1a1a19",
    muted="#6b6b67",
    grid="#e5e4e0",
    ref="#3c3c3a",
    line="#9b9a95",
    losat="#d95f02",
    ncbi="#4d4d4d",
    wasm="#8e2471",
)

DARK = dict(
    surface="#1a1a19",
    ink="#ecebe8",
    muted="#a09f9a",
    grid="#343432",
    ref="#d0cfca",
    line="#7a7974",
    losat="#eb6834",
    ncbi="#b5b4ae",
    wasm="#c04da1",
)

GROUPS_FIG1 = [
    ("Nucleotide vs nucleotide (blastn tasks)", ["mb", "dc", "bn", "bns"]),
    ("Protein vs protein (blastp tasks)", ["bp", "bpff"]),
    ("Translated", ["bx", "tbn", "tbx"]),
    ("PSI & Conserved Domain (Profile)", ["psi", "rps", "rpstbn"]),
]

WASM_GROUPS = [
    ("Nucleotide vs nucleotide (blastn tasks)", ["mb", "dc", "bn", "bns"]),
    ("Protein vs protein (blastp tasks)", ["bp", "bpff"]),
    ("Translated", ["bx", "tbn", "tbx"]),
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fmt_fold(val: float) -> str:
    if val >= 9.95:
        return f"{round(val):.0f}\u00d7"
    return f"{val:.1f}\u00d7"


def read_meta(snap: Path) -> dict:
    p = snap / "metadata.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.is_file() else {}


def short_hw(meta: dict) -> str:
    hw = re.sub(r"\(R\)", "", meta.get("hardware", ""))
    hw = re.sub(r" CPU @ \S+", "", hw)
    return re.sub(r"\s*\(.*\)$", "", hw).strip()


def short_sha(sha: str) -> str:
    return ", ".join(x.strip()[:9] for x in str(sha).split(",") if x.strip())


def footer(fig: plt.Figure, text: str, theme: dict[str, str]) -> None:
    fig.text(0.012, 0.012, text, ha="left", va="bottom", fontsize=7.4, color=theme["muted"])


def plain_log(v: float, _pos: object = None) -> str:
    return f"{v:,.0f}" if v >= 1 else f"{v:g}"


# A downstream-tool ratio is start-up dominated when BLAST+ took under 1 s per BLAST run of the tool's test set.
STARTUP_S = 1.0


def startup_dominated(ncbi_total: str, cases: str) -> bool:
    return bool(ncbi_total) and float(ncbi_total) / max(1, int(cases)) < STARTUP_S


def save_png(fig: plt.Figure, out_path: Path, theme: dict[str, str]) -> None:
    fig.savefig(
        out_path,
        dpi=200,
        facecolor=theme["surface"],
        edgecolor="none",
        metadata={"Software": "LOSAT benchmark renderer"},
    )
    plt.close(fig)


def style_rc(theme: dict[str, str]) -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "text.color": theme["ink"],
        "axes.edgecolor": theme["ref"],
        "axes.labelcolor": theme["ink"],
        "xtick.color": theme["muted"],
        "ytick.color": theme["ink"],
        "svg.fonttype": "none",
        "text.usetex": False,
    })


def panel_title(ax: plt.Axes, letter: str, title: str, theme: dict[str, str], y: float = 1.04) -> None:
    ax.text(0.0, y, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left", color=theme["ink"])
    ax.text(0.065, y, title, transform=ax.transAxes, fontsize=9.5, fontweight="bold", va="bottom", ha="left", color=theme["ink"])


def clean_log_x(ax: plt.Axes, xmin: float, xmax: float, ticks: list[float], labels: list[str], xlabel: str, theme: dict[str, str]) -> None:
    ax.set_xscale("log")
    ax.set_xlim(xmin, xmax)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=8, color=theme["muted"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis="x", which="major", colors=theme["ref"], length=4, width=0.8)
    ax.tick_params(axis="x", which="minor", colors=theme["muted"], length=2.5, width=0.6)
    ax.set_xlabel(xlabel, fontsize=8.5, color=theme["ink"], labelpad=5)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(theme["ref"])
    ax.spines["bottom"].set_linewidth(0.8)
    ax.grid(axis="x", which="major", color=theme["grid"], linewidth=0.7, zorder=0)
    ax.grid(axis="y", which="major", color=theme["grid"], linewidth=0.5, alpha=0.55, zorder=0)


# ---------------------------------------------------------------------------
# Figure 1: Standard workloads (all 12 modes, Cleveland dumbbell plot)
# ---------------------------------------------------------------------------
def render_fig1(snap: Path, out: Path, theme: dict[str, str], dark: bool) -> Path:
    style_rc(theme)
    rows = read_tsv(snap / "standard_workloads.tsv")
    by_kt = {(r["case"], int(r["threads"])): r for r in rows}

    # Build y positions top-down with section headers
    y_rows: list[tuple[str, float, str | None]] = []  # (kind: 'hdr'|'case', y, text/cid)
    y = 0.0
    for gi, (gtitle, cids) in enumerate(GROUPS_FIG1):
        if gi > 0:
            y -= 0.55
        y_rows.append(("hdr", y, gtitle))
        y -= 0.95
        for cid in cids:
            y_rows.append(("case", y, cid))
            y -= 1.0
    ymin, ymax = y + 0.35, 0.7

    fig = plt.figure(figsize=(14.0, 6.4))
    fig.patch.set_facecolor(theme["surface"])

    # Legend at top
    handles = [
        Line2D([0], [0], marker="D", linestyle="none", markersize=6, color=theme["ncbi"], label="BLAST+ 2.17.0 (diamond)"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=6.5, color=theme["losat"], label="LOSAT 0.2.0 (circle)"),
        Line2D([0], [0], marker="D", linestyle="none", markersize=6, markerfacecolor="none", markeredgecolor=theme["ncbi"], label="open: 1 thread"),
        Line2D([0], [0], marker="D", linestyle="none", markersize=6, color=theme["ncbi"], label="filled: 16 threads"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.02, 0.975), frameon=False, ncol=4, fontsize=8.5, columnspacing=1.6, handletextpad=0.5)

    # Dedicated label column + 3 plot+ratio panel pairs (guaranteed zero overlap/clipping)
    ax_lbl = fig.add_axes([0.012, 0.13, 0.195, 0.70])
    ax_lbl.set_facecolor(theme["surface"])
    ax_lbl.set_xlim(0, 1)
    ax_lbl.set_ylim(ymin, ymax)
    ax_lbl.axis("off")

    panel_rects = [
        ([0.212, 0.13, 0.175, 0.70], [0.390, 0.13, 0.056, 0.70]),
        ([0.464, 0.13, 0.175, 0.70], [0.642, 0.13, 0.056, 0.70]),
        ([0.718, 0.13, 0.175, 0.70], [0.896, 0.13, 0.072, 0.70]),
    ]

    specs = [
        ("a", "Wall time, 1 thread", 1, "wall", [5, 1e6], [10, 100, 1000, 10000], ["10", "100", "1,000", "10,000"], "Wall time (s)"),
        ("b", "Wall time, 16 threads", 16, "wall", [0.3, 1e4], [1, 10, 100, 1000], ["1", "10", "100", "1,000"], "Wall time (s)"),
        ("c", "Peak memory (max RSS), 16 threads", 16, "rss", [10, 2e4], [20, 100, 500, 2500, 10000], ["20", "100", "500", "2,500", "10,000"], "Max RSS (MiB)"),
    ]

    for ax_i, ((prect, rrect), (let, title, th, metric, xlim, xticks, xlbls, xlabel)) in enumerate(zip(panel_rects, specs)):
        ax = fig.add_axes(prect)
        ax_r = fig.add_axes(rrect)
        ax.set_facecolor(theme["surface"])
        ax_r.set_facecolor(theme["surface"])
        ax_r.set_xlim(0, 1)
        ax_r.set_ylim(ymin, ymax)
        ax_r.axis("off")

        clean_log_x(ax, xlim[0], xlim[1], xticks, xlbls, xlabel, theme)
        ax.set_ylim(ymin, ymax)
        ax.set_yticks([])
        panel_title(ax, let, title, theme)
        ax_r.text(0.95, 0.0, "BLAST+ \u00f7\nLOSAT", ha="right", va="center", fontsize=7.6, color=theme["muted"])

        vals_all = []
        for kind, yv, txt in y_rows:
            if kind == "hdr":
                if ax_i == 0:
                    ax_lbl.text(0.0, yv, txt, ha="left", va="center", fontsize=8.3, fontweight="bold", color=theme["ink"])
                continue
            row = by_kt.get((txt, th))
            if not row:
                continue
            lv = float(row["losat_wall_s"] if metric == "wall" else row["losat_rss_mib"])
            nv_s = row["ncbi_wall_s"] if metric == "wall" else row["ncbi_rss_mib"]
            nv = float(nv_s) if nv_s else None
            vals_all.append(lv)
            if nv:
                vals_all.append(nv)

        if vals_all:
            lo_x = min(vals_all) * 0.55
            hi_x = max(vals_all) * 1.65
            ax.set_xlim(lo_x, hi_x)
            cand = [0.5, 2, 10, 50, 200, 1000, 5000]
            in_r = [c for c in cand if lo_x * 1.05 <= c <= hi_x * 0.95]
            if len(in_r) < 3:
                cand2 = [0.5, 1, 5, 10, 50, 200, 1000, 5000]
                in_r = [c for c in cand2 if lo_x * 1.05 <= c <= hi_x * 0.95]
            ax.set_xticks(in_r)
            ax.set_xticklabels([f"{int(c):,}" if c >= 1 else f"{c:g}" for c in in_r], fontsize=7.8, color=theme["muted"])

        filled = th == 16
        for kind, yv, txt in y_rows:
            if kind != "case":
                continue
            row = by_kt[(txt, th)]
            if ax_i == 0:
                ax_lbl.text(0.08, yv, row["label"], ha="left", va="center", fontsize=8.3, color=theme["ink"])

            lv = float(row["losat_wall_s"] if metric == "wall" else row["losat_rss_mib"])
            nv_s = row["ncbi_wall_s"] if metric == "wall" else row["ncbi_rss_mib"]
            nv = float(nv_s) if nv_s else None
            if nv is not None:
                ax.plot([lv, nv], [yv, yv], color=theme["line"], linewidth=1.4, zorder=2)
                ax.scatter(
                    [nv], [yv], marker="D", s=38,
                    facecolors=theme["ncbi"] if filled else theme["surface"],
                    edgecolors=theme["ncbi"], linewidths=1.5, zorder=3,
                )
                ratio = nv / lv
                ax_r.text(0.95, yv, fmt_fold(ratio), ha="right", va="center",
                          fontsize=8.0, fontweight="bold" if ratio >= 9.95 else "normal", color=theme["ink"])
            elif metric == "wall" and row.get("parity") == "ncbi_timeout_6000s":
                x_hi = ax.get_xlim()[1] * 0.86
                ax.plot([lv, x_hi], [yv, yv], linestyle=":", color=theme["line"], linewidth=1.3, zorder=2)
                ax.scatter([x_hi], [yv], marker=">", s=32, color=theme["muted"], zorder=3)
                ax.text(x_hi * 0.94, yv + 0.24, ">6,000 s", ha="right", va="bottom", fontsize=7.0, color=theme["muted"])
                ax_r.text(0.95, yv, f">{int(6000 // lv)}\u00d7", ha="right", va="center",
                          fontsize=8.0, fontweight="bold", color=theme["ink"])
            else:
                ax_r.text(0.95, yv, "n.m.", ha="right", va="center", fontsize=7.8, color=theme["muted"])
            ax.scatter(
                [lv], [yv], marker="o", s=46,
                facecolors=theme["losat"] if filled else theme["surface"],
                edgecolors=theme["losat"], linewidths=1.6, zorder=4,
            )

    meta = read_meta(snap)
    line1 = (
        f"LOSAT {short_sha(meta.get('losat_commit', ''))} vs NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')}, "
        f"{short_hw(meta)}. FASTA subjects: LOSAT -subject, BLAST+ -db (made from same FASTA)."
    )
    line2 = "BLAST+ max RSS includes its mapped database pages. Numbers right of each panel: BLAST+ \u00f7 LOSAT (above 1: LOSAT faster or smaller)."
    footer(fig, f"{line1}\n{line2}", theme)
    suffix = "_dark.png" if dark else ".png"
    p = out / f"fig1_standard_workloads{suffix}"
    save_png(fig, p, theme)
    return p


# ---------------------------------------------------------------------------
# Figure 2: Scaling (a Threads, b Query set, c DB time, d DB memory)
# ---------------------------------------------------------------------------
def render_fig2(snap: Path, out: Path, theme: dict[str, str], dark: bool) -> Path:
    style_rc(theme)
    rows = read_tsv(snap / "scaling.tsv")
    fig = plt.figure(figsize=(12.2, 6.8))
    fig.patch.set_facecolor(theme["surface"])
    gs = fig.add_gridspec(2, 3, left=0.065, right=0.975, top=0.86, bottom=0.125, wspace=0.31, hspace=0.48)

    handles = [
        Line2D([0], [0], marker="D", color=theme["ncbi"], markersize=5.5, linewidth=1.5, label="BLAST+ 2.17.0"),
        Line2D([0], [0], marker="o", color=theme["losat"], markersize=6, linewidth=1.5, label="LOSAT 0.2.0"),
        Line2D([0], [0], linestyle=":", color=theme["muted"], linewidth=1.1, label="ideal linear thread scaling"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.04, 0.975), frameon=False, ncol=3, fontsize=8.5)

    # Row 1: Panel a (3 subplots for bn, bp, tbn)
    a_specs = [
        ("bn", "blastn  31,168 cDNAs vs 3 Mb chr. 2L", True),
        ("bp", "blastp  3,000 human vs yeast proteome", False),
        ("tbn", "tblastn  1 human protein vs GRCh38", False),
    ]
    a_rows = [r for r in rows if r["panel"] == "2a_threads"]
    for idx, (cid, subtitle, is_first) in enumerate(a_specs):
        ax = fig.add_subplot(gs[0, idx])
        ax.set_facecolor(theme["surface"])
        sub = sorted([r for r in a_rows if r["case"] == cid], key=lambda x: int(x["x_value"]))
        xs = [int(r["x_value"]) for r in sub]
        ys_n = [float(r["ncbi_wall_s"]) for r in sub]
        ys_l = [float(r["losat_wall_s"]) for r in sub]

        # Ideal scaling dotted lines from 1t
        ax.plot(xs, [ys_n[0] / x for x in xs], linestyle=":", color=theme["ncbi"], alpha=0.55, linewidth=1.0, zorder=1)
        ax.plot(xs, [ys_l[0] / x for x in xs], linestyle=":", color=theme["losat"], alpha=0.65, linewidth=1.0, zorder=1)
        ax.plot(xs, ys_n, marker="D", markersize=4.8, color=theme["ncbi"], linewidth=1.5, zorder=3)
        ax.plot(xs, ys_l, marker="o", markersize=5.3, color=theme["losat"], linewidth=1.5, zorder=4)

        for x, yn, yl in zip(xs, ys_n, ys_l):
            if x in (1, 16):
                # inside the plot, clear of the y axis (x = 1) and the right edge (x = 16)
                ax.text(x * 1.1 if x == 1 else x / 1.1, math.sqrt(yn * yl), fmt_fold(yn / yl),
                        ha="left" if x == 1 else "right", va="center", fontsize=7.8,
                        color=theme["ink"], bbox=dict(boxstyle="round,pad=0.12", facecolor=theme["surface"], edgecolor="none", alpha=0.82))

        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xticks([1, 2, 4, 8, 16])
        ax.set_xticklabels(["1", "2", "4", "8", "16"], fontsize=8, color=theme["muted"])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())
        y_lo = min(ys_l) * 0.45
        y_hi = max(ys_n) * 1.9
        ax.set_ylim(y_lo, y_hi)
        yt = [c for c in [0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500] if y_lo <= c <= y_hi]
        ax.set_yticks(yt)
        ax.set_yticklabels([f"{int(c):,}" if c >= 1 else f"{c:g}" for c in yt], fontsize=7.8, color=theme["muted"])
        ax.set_xlabel("Threads", fontsize=8.3)
        if idx == 0:
            ax.set_ylabel("Wall time (s)", fontsize=8.3)
            panel_title(ax, "a", "Thread scaling", theme, y=1.14)
        ax.set_title(subtitle, fontsize=8.3, color=theme["ink"], pad=7)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.grid(True, which="major", color=theme["grid"], linewidth=0.6, zorder=0)

    # Panel b: Query set size (megablast reads vs 16S)
    ax_b = fig.add_subplot(gs[1, 0])
    ax_b.set_facecolor(theme["surface"])
    b_rows = sorted([r for r in rows if r["panel"] == "2b_query_size"], key=lambda x: float(x["x_value"]))
    xs_b = [float(r["x_value"]) for r in b_rows]
    ys_bn = [float(r["ncbi_wall_s"]) for r in b_rows]
    ys_bl = [float(r["losat_wall_s"]) for r in b_rows]
    ax_b.plot(xs_b, ys_bn, marker="D", markersize=5.0, color=theme["ncbi"], linewidth=1.6, zorder=3)
    ax_b.plot(xs_b, ys_bl, marker="o", markersize=5.5, color=theme["losat"], linewidth=1.6, zorder=4)
    for x, yn, yl in zip(xs_b, ys_bn, ys_bl):
        ax_b.text(x, math.sqrt(yn * yl), fmt_fold(yn / yl), ha="center", va="center", fontsize=7.8,
                  color=theme["ink"], bbox=dict(boxstyle="round,pad=0.12", facecolor=theme["surface"], edgecolor="none", alpha=0.82))
    ax_b.set_xscale("log")
    ax_b.set_yscale("log")
    ax_b.set_xlim(5000, 2e6)
    ax_b.set_xticks([1e4, 1e5, 1e6])
    ax_b.set_xticklabels(["10,000", "100,000", "1,000,000"], fontsize=7.8, color=theme["muted"])
    ax_b.xaxis.set_minor_formatter(NullFormatter())
    ax_b.yaxis.set_minor_formatter(NullFormatter())
    ax_b.yaxis.set_major_formatter(FuncFormatter(plain_log))
    ax_b.set_ylabel("Wall time (s)", fontsize=8.3)
    ax_b.set_xlabel("Query reads", fontsize=8.3)
    panel_title(ax_b, "b", "Query set size", theme, y=1.14)
    ax_b.set_title("megablast  reads vs 16S rRNA, 16 threads", fontsize=8.3, color=theme["ink"], pad=7)
    for sp in ("top", "right"):
        ax_b.spines[sp].set_visible(False)
    ax_b.grid(True, which="major", color=theme["grid"], linewidth=0.6, zorder=0)

    # Panel c & d: Database size (warm + faded cold)
    cd_rows = [r for r in rows if r["panel"] == "2cd_db_size"]
    warm = sorted([r for r in cd_rows if r["cache"] == "warm"], key=lambda x: float(x["x_value"]))
    cold = sorted([r for r in cd_rows if r["cache"] == "cold"], key=lambda x: float(x["x_value"]))

    for p_i, (ax, let, title, sub, ylab, key_n, key_l) in enumerate([
        (fig.add_subplot(gs[1, 1]), "c", "Database size, wall time", "megablast  210 contigs, 16 threads", "Wall time (s)", "ncbi_wall_s", "losat_wall_s"),
        (fig.add_subplot(gs[1, 2]), "d", "Database size, peak memory", "megablast  same runs as c", "Max RSS (MiB)", "ncbi_rss_mib", "losat_rss_mib"),
    ]):
        ax.set_facecolor(theme["surface"])
        xs_w = [float(r["x_value"]) for r in warm]
        yn_w = [float(r[key_n]) for r in warm]
        yl_w = [float(r[key_l]) for r in warm]
        ax.plot(xs_w, yn_w, marker="D", markersize=5.0, color=theme["ncbi"], linewidth=1.6, zorder=3)
        ax.plot(xs_w, yl_w, marker="o", markersize=5.5, color=theme["losat"], linewidth=1.6, zorder=4)
        for x, yn, yl in zip(xs_w, yn_w, yl_w):
            ratio = yn / yl
            if ratio >= 2.8:
                y_pos, va = math.sqrt(yn * yl), "center"
            elif p_i == 0:
                y_pos, va = yl / 1.75, "top"
            else:
                # Below LOSAT's gentler curve keeps small-DB ratios (1.7x, 1.3x) clear of BLAST+'s steep mmap rise
                y_pos, va = yl / 1.55, "top"
            ax.text(x, y_pos, fmt_fold(ratio), ha="center", va=va, fontsize=7.6,
                    color=theme["ink"], bbox=dict(boxstyle="round,pad=0.10", facecolor=theme["surface"], edgecolor="none", alpha=0.85))

        if cold:
            xs_c = [float(r["x_value"]) for r in cold]
            yn_c = [float(r[key_n]) for r in cold]
            yl_c = [float(r[key_l]) for r in cold]
            ax.scatter(xs_c, yn_c, marker="D", s=34, color=theme["ncbi"], alpha=0.38, zorder=2)
            ax.scatter(xs_c, yl_c, marker="o", s=38, color=theme["losat"], alpha=0.45, zorder=2)
            if p_i == 0:
                ax.text(xs_c[0] * 0.55, yn_c[0] * 1.15, "cold cache", ha="right", va="bottom", fontsize=7.3, color=theme["muted"])
        if p_i == 0:
            ax.set_ylim(min(yl_w) * 0.32, max(yn_w + ([yn_c[-1]] if cold else [])) * 2.2)
        else:
            ax.set_ylim(min(yl_w) * 0.35, max(yn_w) * 2.4)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(0.05, 1.2e6)
        ax.set_xticks([1, 100, 1e4, 1e6])
        ax.set_xticklabels(["1 MB", "100 MB", "10 GB", "1 TB"], fontsize=7.8, color=theme["muted"])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(FuncFormatter(plain_log))
        ax.set_ylabel(ylab, fontsize=8.3)
        ax.set_xlabel("Database size on disk", fontsize=8.3)
        panel_title(ax, let, title, theme, y=1.14)
        ax.set_title(sub, fontsize=8.3, color=theme["ink"], pad=7)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.grid(True, which="major", color=theme["grid"], linewidth=0.6, zorder=0)

    meta = read_meta(snap)
    abr_row = next((r for r in warm if r["case"] == "m3_abr"), None)
    abr_note = ""
    if abr_row and float(abr_row["ncbi_wall_s"]) < STARTUP_S:
        abr_note = f"; on abr_genes ({float(abr_row['x_value']):.2f} MB) BLAST+ took {float(abr_row['ncbi_wall_s']):.2f} s (start-up)"
    tm = short_sha(meta.get("toolmt_commit", ""))
    line1 = f"LOSAT {short_sha(meta.get('losat_commit', ''))} vs NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')}, {short_hw(meta)}. Labels: BLAST+ \u00f7 LOSAT."
    if tm and tm != short_sha(meta.get("losat_commit", "")):
        line2 = f"c, d: core_nt and cold-cache points are LOSAT {tm}; BLAST+ max RSS includes its mapped database pages{abr_note}."
    else:
        line2 = f"c, d: BLAST+ max RSS includes its mapped database pages{abr_note}."
    footer(fig, f"{line1}\n{line2}", theme)
    suffix = "_dark.png" if dark else ".png"
    p = out / f"fig2_scaling{suffix}"
    save_png(fig, p, theme)
    return p


# ---------------------------------------------------------------------------
# Figure 3: Downstream tools (2-column Cleveland dumbbell plot)
# ---------------------------------------------------------------------------
def render_fig3(snap: Path, out: Path, theme: dict[str, str], dark: bool) -> Path:
    style_rc(theme)
    rows = read_tsv(snap / "downstream_tools.tsv")
    by_prog: dict[str, list[dict[str, str]]] = defaultdict(list)
    for r in rows:
        by_prog[r["program"]].append(r)

    def sort_key(r: dict[str, str]) -> tuple[int, float]:
        is_mt = 0 if r["mt_speedup"] or r["note"] == "losat_under_10ms" else 1
        sp = float(r["mt_speedup"]) if r["mt_speedup"] else (0.0 if r["note"] == "losat_under_10ms" else float(r["n1_speedup"]))
        return (is_mt, -sp)

    left_progs = ["blastn", "tblastx"]
    right_progs = ["blastp", "blastx", "tblastn"]

    fig = plt.figure(figsize=(14.4, 6.7))
    fig.patch.set_facecolor(theme["surface"])

    handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=6.5, markerfacecolor=theme["surface"], markeredgecolor=theme["losat"], markeredgewidth=1.6, label="1 thread"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=6.5, color=theme["losat"], label="threads passed by tool"),
        Line2D([0], [0], marker="|", markersize=11, markeredgewidth=1.4, linestyle="none", color=theme["ref"], label="BLAST+ 2.17.0 = 1"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.02, 0.975), frameon=False, ncol=3, fontsize=8.5)

    # Two independent halves with dedicated label + plot + number axes (zero overlap guaranteed)
    col_rects = [
        ([0.020, 0.10, 0.175, 0.77], [0.200, 0.10, 0.210, 0.77], [0.415, 0.10, 0.068, 0.77]),
        ([0.520, 0.10, 0.155, 0.77], [0.680, 0.10, 0.215, 0.77], [0.900, 0.10, 0.075, 0.77]),
    ]

    for (lrect, prect, rrect), progs in zip(col_rects, (left_progs, right_progs)):
        ax_l = fig.add_axes(lrect)
        ax = fig.add_axes(prect)
        ax_r = fig.add_axes(rrect)
        for subax in (ax_l, ax_r):
            subax.set_facecolor(theme["surface"])
            subax.set_xlim(0, 1)
            subax.axis("off")

        ax.set_facecolor(theme["surface"])
        clean_log_x(ax, 0.92, 26.0, [1, 2, 3, 5, 10, 20], ["1", "2", "3", "5", "10", "20"],
                    "BLAST+ time \u00f7 LOSAT time", theme)
        ax.axvline(1.0, color=theme["ref"], linewidth=1.1, zorder=1)
        ax.set_yticks([])

        y_items: list[tuple[str, float, object]] = []
        y = 0.0
        for pi, p in enumerate(progs):
            if pi > 0:
                y -= 0.75
            y_items.append(("hdr", y, p))
            y -= 0.95
            for r in sorted(by_prog.get(p, []), key=sort_key):
                y_items.append(("row", y, r))
                y -= 1.0
        ylims = (y + 0.35, 1.35)
        ax.set_ylim(*ylims)
        ax_l.set_ylim(*ylims)
        ax_r.set_ylim(*ylims)

        ax_l.text(0.0, 0.85, "Tool", ha="left", va="center", fontsize=8.2, color=theme["muted"])
        ax_l.text(0.86, 0.85, "Threads", ha="center", va="center", fontsize=8.2, color=theme["muted"])
        ax_r.scatter([0.36], [0.85], marker="o", s=34, facecolors=theme["surface"], edgecolors=theme["losat"], linewidths=1.5)
        ax_r.scatter([0.86], [0.85], marker="o", s=36, color=theme["losat"])

        for kind, yv, item in y_items:
            if kind == "hdr":
                ax_l.text(0.0, yv, str(item), ha="left", va="center", fontsize=8.6, fontweight="bold", color=theme["ink"])
                continue
            r = item
            ax_l.text(0.04, yv, r["tool"], ha="left", va="center", fontsize=8.1, color=theme["ink"])
            ax_l.text(0.86, yv, r["threads_label"], ha="center", va="center", fontsize=8.0, color=theme["muted"])

            if r["note"] == "losat_under_10ms":
                ax.text(1.12, yv, "LOSAT under 10 ms", ha="left", va="center", fontsize=7.8, style="italic", color=theme["muted"])
                ax_r.text(0.45, yv, "\u2014", ha="right", va="center", fontsize=8.0, color=theme["muted"])
                ax_r.text(0.95, yv, "\u2014", ha="right", va="center", fontsize=8.0, color=theme["muted"])
                continue

            x1 = float(r["n1_speedup"]) if r["n1_speedup"] else None
            xm = float(r["mt_speedup"]) if r["mt_speedup"] else None
            st1 = startup_dominated(r["ncbi_n1_s"], r["cases_count"])
            stm = startup_dominated(r["ncbi_mt_s"], r["cases_count"])
            near_overlap = x1 is not None and xm is not None and abs(math.log(xm / x1)) < 0.06
            y1 = yv + 0.14 if near_overlap else yv
            ym = yv - 0.14 if near_overlap else yv
            if x1 is not None and xm is not None:
                ax.plot([x1, xm], [y1, ym], color=theme["line"], linewidth=1.4, zorder=2)
            if x1 is not None:
                ax.scatter([x1], [y1], marker="o", s=40, facecolors=theme["surface"], edgecolors=theme["losat"], linewidths=1.5, zorder=3)
                ax_r.text(0.45, yv, fmt_fold(x1) + ("*" if st1 else ""), ha="right", va="center", fontsize=7.9,
                          fontweight="bold" if x1 >= 9.95 and not st1 else "normal", color=theme["muted"] if st1 else theme["ink"])
            else:
                ax_r.text(0.45, yv, "\u2014", ha="right", va="center", fontsize=8.0, color=theme["muted"])
            if xm is not None:
                ax.scatter([xm], [ym], marker="o", s=44, color=theme["losat"], zorder=4)
                ax_r.text(0.95, yv, fmt_fold(xm) + ("*" if stm else ""), ha="right", va="center", fontsize=7.9,
                          fontweight="bold" if xm >= 9.95 and not stm else "normal", color=theme["muted"] if stm else theme["ink"])

    meta = read_meta(snap)
    dt = meta.get("dt_commits", {})
    one = short_sha(",".join(dt[k] for k in sorted(dt) if dt[k]))
    footer(fig, f"Threads passed by the tool: LOSAT {short_sha(meta.get('toolmt_commit', ''))}; 1 thread: LOSAT {one}; "
                f"NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')}, {short_hw(meta)}. "
                f"* (grey) BLAST+ took under {STARTUP_S:g} s per run: start-up dominated, not a search-speed ratio.", theme)
    suffix = "_dark.png" if dark else ".png"
    p = out / f"fig3_downstream_tools{suffix}"
    save_png(fig, p, theme)
    return p


# ---------------------------------------------------------------------------
# Figure 4: WebAssembly vs native vs native BLAST+
# ---------------------------------------------------------------------------
def render_fig4(snap: Path, out: Path, theme: dict[str, str], dark: bool) -> Path:
    style_rc(theme)
    rows = read_tsv(snap / "webassembly.tsv")
    by_kt = {(r["case"], int(r["threads"])): r for r in rows}

    y_rows: list[tuple[str, float, str | None]] = []
    y = 0.0
    for gi, (gtitle, cids) in enumerate(WASM_GROUPS):
        if gi > 0:
            y -= 0.55
        y_rows.append(("hdr", y, gtitle))
        y -= 0.95
        for cid in cids:
            y_rows.append(("case", y, cid))
            y -= 1.0
    ymin, ymax = y + 0.35, 0.7

    meta = read_meta(snap)
    wasm_derived = bool(meta.get("wasm_derived", False))
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.5 if wasm_derived else 5.4), sharey=True)
    fig.patch.set_facecolor(theme["surface"])
    fig.subplots_adjust(left=0.27, right=0.895, top=0.785 if wasm_derived else 0.815, bottom=0.155, wspace=0.38)

    handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=6.5, color=theme["losat"], label="LOSAT native"),
        Line2D([0], [0], marker="^", linestyle="none", markersize=6.8, color=theme["wasm"], label="LOSAT WebAssembly"),
        Line2D([0], [0], marker="|", markersize=11, markeredgewidth=1.4, linestyle="none", color=theme["ref"], label="native BLAST+ 2.17.0 = 1"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=6.0, markerfacecolor=theme["surface"], markeredgecolor=theme["losat"], markeredgewidth=1.5, label="open: 1 thread"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=6.0, color=theme["losat"], label="filled: 16 threads"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.03, 0.975), frameon=False, ncol=3, fontsize=8.3)

    for ax_i, (ax, let, title, th) in enumerate(zip(axes, ("a", "b"), ("1 thread", "16 threads"), (1, 16))):
        ax.set_facecolor(theme["surface"])
        clean_log_x(ax, 0.5, 25.0, [0.5, 1, 2, 5, 10, 20], ["0.5", "1", "2", "5", "10", "20"],
                    "BLAST+ time \u00f7 LOSAT time", theme)
        ax.axvline(1.0, color=theme["ref"], linewidth=1.1, zorder=1)
        ax.set_ylim(ymin, ymax)
        ax.set_yticks([])
        panel_title(ax, let, title, theme)
        ax.text(1.03, ymax - 0.15, "Wasm \u00f7 native", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=8.0, color=theme["muted"])

        filled = th == 16
        lbl = blended_transform_factory(fig.transFigure, ax.transData)  # labels in the left margin, clear of the plot
        for kind, yv, txt in y_rows:
            if kind == "hdr":
                if ax_i == 0:
                    ax.text(0.012, yv, txt, transform=lbl, ha="left", va="center", fontsize=8.2, fontweight="bold", color=theme["ink"])
                continue
            r = by_kt[(txt, th)]
            if ax_i == 0:
                ax.text(0.03, yv, r["mode_label"], transform=lbl, ha="left", va="center", fontsize=8.3, color=theme["ink"])
            ov = float(r["wasm_overhead"])
            ax.text(1.03, yv, f"{ov:.1f}\u00d7", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=8.0, color=theme["ink"])

            if not r["native_speedup"]:
                ax.text(2.5, yv, "BLAST+ not measured", ha="left", va="center", fontsize=7.6, color=theme["muted"])
                continue
            ns = float(r["native_speedup"])
            ws = float(r["wasm_speedup"])
            if ws <= 22.0 and ns <= 22.0:
                ax.plot([ws, ns], [yv, yv], color=theme["line"], linewidth=1.4, zorder=2)
                ax.scatter([ws], [yv], marker="^", s=46, facecolors=theme["wasm"] if filled else theme["surface"], edgecolors=theme["wasm"], linewidths=1.5, zorder=3)
                ax.scatter([ns], [yv], marker="o", s=44, facecolors=theme["losat"] if filled else theme["surface"], edgecolors=theme["losat"], linewidths=1.5, zorder=4)
            else:
                ax.plot([max(1.0, min(ws, 18.0)), 22.0], [yv, yv], linestyle=":", color=theme["line"], linewidth=1.4, zorder=2)
                ax.text(10.0, yv + 0.25, f"{fmt_fold(ws)} and {fmt_fold(ns)} (off scale)", ha="center", va="bottom", fontsize=7.6, color=theme["ink"])

    wasm_sha = short_sha(meta.get("wasm_commit", ""))
    base_sha = short_sha(meta.get("wasm_derived_from", "6a0a8719"))
    if wasm_derived:
        fig.text(
            0.975,
            0.942,
            f"Estimated from build {base_sha} (see caption)",
            ha="right",
            va="top",
            fontsize=8.0,
            style="italic",
            color=theme["muted"],
        )
    line1 = (
        f"LOSAT {wasm_sha}, native and wasm32-wasip1-threads (Node.js), "
        f"vs native NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')}, {short_hw(meta)}."
    )
    if wasm_derived:
        line2 = f"Wasm \u00f7 native: WebAssembly time over native LOSAT time. {wasm_sha} times: {base_sha}'s measured times scaled by a side-by-side run of the two builds."
    else:
        line2 = "Wasm \u00f7 native: WebAssembly time over native LOSAT time."
    footer(fig, f"{line1}\n{line2}", theme)
    suffix = "_dark.png" if dark else ".png"
    p = out / f"fig4_webassembly{suffix}"
    save_png(fig, p, theme)
    return p


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--snapshot", type=Path, default=Path(__file__).resolve().parent.parent / "snapshot")
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent.parent / "figures")
    ap.add_argument("--only", choices=["fig1", "fig2", "fig3", "fig4"], default=None)
    ap.add_argument("--manual-images", default=str(Path(__file__).resolve().parents[3] / "docs/manual/src/images"),
                    help="also copy every PNG here (the user manual's copy); '' to skip")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    mirror = Path(args.manual_images) if args.manual_images else None
    if mirror:
        mirror.mkdir(parents=True, exist_ok=True)
    renderers = {
        "fig1": render_fig1,
        "fig2": render_fig2,
        "fig3": render_fig3,
        "fig4": render_fig4,
    }
    targets = [args.only] if args.only else ["fig1", "fig2", "fig3", "fig4"]
    outputs = []
    for k in targets:
        for dark, theme in ((False, LIGHT), (True, DARK)):
            p = renderers[k](args.snapshot, args.out, theme, dark)
            outputs.append(p)
            if mirror:
                shutil.copyfile(p, mirror / p.name)
            print(p.name, sha256_file(p)[:16])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
