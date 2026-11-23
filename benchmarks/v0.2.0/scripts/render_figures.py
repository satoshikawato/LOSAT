#!/usr/bin/env python3
"""Deterministically render the 6-figure LOSAT v0.2.0 benchmark suite (light and dark PNGs).

Figures (blank-slate redesign):
  fig1_suite_overview   : All 12 search modes + 3 database/domain utilities (16t speed, 1t speed, 16t memory)
  fig2_real_pipelines   : Large-scale reference workflows & published bioinformatics tools (>= 1 s/run)
  fig3_scaling          : Scaling with database size (warm/cold), query volume (time + RAM), and threads
  fig4_vs_heuristics    : Exact BLAST+ parity vs DIAMOND, MMseqs2 (UniRef50) and vsearch (SILVA 144)
  fig5_webassembly      : Portable WebAssembly (wasm32-wasip1-threads) vs native LOSAT & native BLAST+
  fig6_v01_to_v02       : Generational leap from LOSAT v0.1.0 to v0.2.0 against NCBI BLAST+ 2.17.0

Writes benchmarks/v0.2.0/figures/ and mirrors the exact PNG bytes to docs/manual/src/images/.
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
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, NullFormatter

# ---------------------------------------------------------------------------
# Colorblind-safe publication palette (Light & Dark)
# ---------------------------------------------------------------------------
LIGHT = dict(
    surface="#FFFFFF",
    band="#F6F8FA",
    zone="#FFF5F0",
    zone_edge="#FCD9C5",
    ink="#0F172A",
    ink2="#334155",
    muted="#64748B",
    grid="#E2E8F0",
    hairline="#CBD5E1",
    ref="#334155",
    losat="#D94E1F",
    losat_soft="#F3A285",
    losat_fill="#FBE4D8",
    ncbi="#475569",
    ncbi_soft="#CBD5E1",
    wasm="#0F766E",
    wasm_soft="#99F6E4",
    diamond="#2563EB",
    mmseqs="#7C3AED",
    vsearch="#B45309",
    v010="#94A3B8",
    badge_bg="#EEF2F6",
    badge_ink="#475569",
)

DARK = dict(
    surface="#0F1117",
    band="#161B26",
    zone="#221612",
    zone_edge="#542B1D",
    ink="#F1F5F9",
    ink2="#CBD5E1",
    muted="#94A3B8",
    grid="#242C3D",
    hairline="#334155",
    ref="#CBD5E1",
    losat="#F27441",
    losat_soft="#9E4A29",
    losat_fill="#3A1F16",
    ncbi="#94A3B8",
    ncbi_soft="#475569",
    wasm="#2DD4BF",
    wasm_soft="#115E59",
    diamond="#60A5FA",
    mmseqs="#A78BFA",
    vsearch="#FBBF24",
    v010="#64748B",
    badge_bg="#1E2534",
    badge_ink="#94A3B8",
)

GROUPS_FIG1 = [
    ("Nucleotide", ["mb", "dc", "bn", "bns"]),
    ("Protein", ["bp", "bpff"]),
    ("Translated", ["bx", "tbn", "tbx"]),
    ("Profile & Domain", ["psi", "rps", "rpstbn"]),
]

WASM_GROUPS = GROUPS_FIG1[:3]

STARTUP_S = 1.0


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_meta(snap: Path) -> dict:
    p = snap / "metadata.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.is_file() else {}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def short_hw(meta: dict) -> str:
    hw = re.sub(r"\(R\)", "", meta.get("hardware", ""))
    hw = re.sub(r" CPU @ \S+", "", hw)
    return re.sub(r"\s*\(.*\)$", "", hw).strip()


def short_sha(sha: str) -> str:
    return ", ".join(x.strip()[:9] for x in str(sha).split(",") if x.strip())


def fmt_fold(val: float, exact_tenth: bool = False) -> str:
    if exact_tenth:
        return f"{val:.1f}\u00d7" if val < 100 else f"{round(val):,.0f}\u00d7"
    if val >= 99.5:
        return f"{round(val):,.0f}\u00d7"
    if val >= 9.95:
        return f"{val:.0f}\u00d7"
    if 0.95 <= val < 1.05:
        return f"{val:.2f}\u00d7"
    return f"{val:.1f}\u00d7"


def fmt_sec(val: float) -> str:
    if val >= 1000:
        return f"{val:,.0f}s"
    if val >= 100:
        return f"{val:.0f}s"
    if val >= 10:
        return f"{val:.1f}s"
    return f"{val:.2f}s"


def fmt_mem(mib: float, compact: bool = False) -> str:
    if mib >= 10240:
        return f"{mib / 1024.0:.0f}{'G' if compact else ' GiB'}"
    if mib >= 1024:
        return f"{mib / 1024.0:.1f}{'G' if compact else ' GiB'}"
    if mib >= 100:
        return f"{mib:.0f}{'M' if compact else ' MiB'}"
    if compact:
        return f"{mib:.0f}M"
    return f"{mib:.1f} MiB"


def plain_log(v: float, _pos: object = None) -> str:
    return f"{v:,.0f}" if v >= 1 else f"{v:g}"


def startup_dominated(ncbi_total: str, cases: str) -> bool:
    return bool(ncbi_total) and float(ncbi_total) / max(1, int(cases)) < STARTUP_S


def style_rc(theme: dict[str, str]) -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Lato", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8.6,
        "text.color": theme["ink"],
        "axes.edgecolor": theme["hairline"],
        "axes.labelcolor": theme["ink2"],
        "xtick.color": theme["ink2"],
        "ytick.color": theme["ink"],
        "svg.fonttype": "none",
        "text.usetex": False,
    })


def save_png(fig: plt.Figure, out_path: Path, theme: dict[str, str]) -> None:
    fig.savefig(
        out_path,
        dpi=200,
        facecolor=theme["surface"],
        edgecolor="none",
        metadata={"Software": "LOSAT benchmark renderer v2"},
    )
    plt.close(fig)


def new_fig(figsize: tuple[float, float], theme: dict[str, str]) -> plt.Figure:
    fig = plt.figure(figsize=figsize)
    fig.patch.set_facecolor(theme["surface"])
    return fig


def save_themed(fig: plt.Figure, out: Path, stem: str, dark: bool, theme: dict[str, str]) -> Path:
    p = out / f"{stem}{'_dark.png' if dark else '.png'}"
    save_png(fig, p, theme)
    return p


def load_std_by_kt(snap: Path) -> tuple[list[dict[str, str]], dict[tuple[str, int], dict[str, str]]]:
    std = read_tsv(snap / "standard_workloads.tsv")
    return std, {(r["case"], int(r["threads"])): r for r in std}


def ncbi_ref_handle(theme: dict[str, str]) -> Line2D:
    return Line2D([0], [0], color=theme["ref"], linewidth=1.2, label="NCBI BLAST+ 2.17.0 = 1\u00d7")


def hide_minor_x(ax: plt.Axes) -> None:
    ax.xaxis.set_minor_formatter(NullFormatter())


def style_axis_log_x(
    ax: plt.Axes,
    ticks: list[float],
    theme: dict[str, str],
    labels: list[str] | None = None,
    tick_fs: float = 7.6,
    vgrid: bool = True,
    skip_one: bool = False,
    bound_lo: float | None = None,
    bound_hi: float | None = None,
    major_len: float | None = None,
) -> None:
    if vgrid:
        for tk in ticks:
            if not (skip_one and tk == 1):
                ax.axvline(tk, color=theme["grid"], linewidth=0.6, zorder=1)
    ax.set_xticks(ticks)
    if labels is not None:
        ax.set_xticklabels(labels, fontsize=tick_fs, color=theme["ink2"])
    hide_minor_x(ax)
    ax.tick_params(axis="x", which="minor", length=0)
    if major_len is not None:
        ax.tick_params(axis="x", which="major", color=theme["hairline"], labelcolor=theme["ink2"], length=major_len, width=0.8)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(theme["hairline"])
    if bound_lo is not None and bound_hi is not None:
        ax.spines["bottom"].set_bounds(bound_lo, bound_hi)


def style_clean_grid(
    ax: plt.Axes, theme: dict[str, str], lw: float = 0.6, axis: str = "both", plain_y: bool = False, clear_minor_x: bool = True
) -> None:
    if clear_minor_x:
        hide_minor_x(ax)
    ax.yaxis.set_minor_formatter(NullFormatter())
    if plain_y:
        ax.yaxis.set_major_formatter(FuncFormatter(plain_log))
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(True, axis=axis, which="major", color=theme["grid"], linewidth=lw, zorder=0)


def pair_sec(row: dict[str, str]) -> str:
    return f"{fmt_sec(float(row['losat_wall_s']))} vs {fmt_sec(float(row['ncbi_wall_s']))}"


def pair_mem(row: dict[str, str]) -> str:
    return f"{fmt_mem(float(row['losat_rss_mib']))} vs {fmt_mem(float(row['ncbi_rss_mib']))}"


def draw_bar_with_labels(
    ax: plt.Axes,
    yv: float,
    val: float,
    cap: float,
    color: str,
    detail: str,
    theme: dict[str, str],
    val_fs: float = 8.3,
    det_fs: float = 6.7,
    dy_top: float = 0.18,
    dy_bot: float = 0.20,
    bold_at: float | None = 9.5,
    val_anon: float | None = None,
    soft_color: str | None = None,
) -> None:
    if val_anon is not None and soft_color is not None:
        bar_end = draw_broken_bar(ax, yv, val, cap, soft_color, theme, height=0.54)
        draw_broken_bar(ax, yv, val_anon, cap, color, theme, height=0.30)
    else:
        bar_end = draw_broken_bar(ax, yv, val, cap, color, theme, height=0.54)
    lx = bar_end * 1.13
    wt = "bold" if (bold_at is not None and val >= bold_at) else "semibold"
    ax.text(lx, yv + dy_top, fmt_fold(val), ha="left", va="center", fontsize=val_fs, fontweight=wt, color=theme["ink"], zorder=7)
    ax.text(lx, yv - dy_bot, detail, ha="left", va="center", fontsize=det_fs, color=theme["muted"], zorder=7)


def draw_row_badge(
    ax_lbl: plt.Axes, yv: float, title: str, desc: str, badge: str, theme: dict[str, str], title_fs: float = 8.8
) -> None:
    ax_lbl.text(0.02, yv + 0.19, title, ha="left", va="center", fontsize=title_fs, fontweight="bold", color=theme["ink"])
    ax_lbl.text(0.02, yv - 0.22, desc, ha="left", va="center", fontsize=7.0, color=theme["muted"])
    ax_lbl.text(
        0.96, yv + 0.04, badge, ha="right", va="center", fontsize=6.8, color=theme["badge_ink"],
        bbox=dict(boxstyle="round,pad=0.22", facecolor=theme["badge_bg"], edgecolor=theme["grid"], linewidth=0.6),
    )


def add_header(fig: plt.Figure, title: str, subtitle: str, theme: dict[str, str], y_title: float = 0.968, y_sub: float = 0.935) -> None:
    fig.text(0.022, y_title, title, ha="left", va="top", fontsize=13.0, fontweight="bold", color=theme["ink"])
    fig.text(0.022, y_sub, subtitle, ha="left", va="top", fontsize=8.6, color=theme["muted"])
    line = Line2D([0.022, 0.978], [y_sub - 0.030, y_sub - 0.030], transform=fig.transFigure, color=theme["grid"], linewidth=0.9)
    fig.add_artist(line)


def add_footer(fig: plt.Figure, text: str, theme: dict[str, str], y_rule: float = 0.060) -> None:
    line = Line2D([0.022, 0.978], [y_rule, y_rule], transform=fig.transFigure, color=theme["grid"], linewidth=0.7)
    fig.add_artist(line)
    fig.text(0.022, y_rule - 0.010, text, ha="left", va="top", fontsize=7.1, color=theme["muted"], linespacing=1.35)


def add_panel_label(ax: plt.Axes, letter: str, title: str, theme: dict[str, str], y: float = 1.03, sub: str = "") -> None:
    ax.text(0.0, y, f"{letter}  {title}", transform=ax.transAxes, fontsize=9.6, fontweight="bold",
            va="bottom", ha="left", color=theme["ink"])
    if sub:
        ax.text(1.0, y + 0.003, sub, transform=ax.transAxes, fontsize=7.5, va="bottom", ha="right", color=theme["muted"])


def draw_broken_bar(
    ax: plt.Axes,
    yv: float,
    val: float,
    cap: float,
    color: str,
    theme: dict[str, str],
    height: float = 0.50,
    alpha: float = 1.0,
) -> float:
    """Draw a horizontal log-scale bar from 1.0 to val, capping with a clean // break if val > cap."""
    bar_end = min(val, cap * 0.96) if val > cap else val
    if bar_end >= 1.0:
        ax.barh([yv], [bar_end - 1.0], left=1.0, height=height, color=color, alpha=alpha, edgecolor="none", zorder=3)
    else:
        ax.barh([yv], [1.0 - bar_end], left=bar_end, height=height, color=theme["muted"], alpha=0.55, edgecolor="none", zorder=3)
    if val > cap:
        log_pos = math.exp(math.log(1.0) + 0.78 * (math.log(cap) - math.log(1.0)))
        dx = 0.08
        for shift in (-0.045, 0.045):
            x0 = log_pos * math.exp(shift - dx * 0.35)
            x1 = log_pos * math.exp(shift + dx * 0.35)
            ax.plot([x0, x1], [yv - height * 0.62, yv + height * 0.62], color=theme["surface"], linewidth=3.2, zorder=5)
            ax.plot([x0, x1], [yv - height * 0.62, yv + height * 0.62], color=theme["ink2"], linewidth=0.95, zorder=6)
    return bar_end


# ---------------------------------------------------------------------------
# Figure 1: Complete LOSAT v0.2.0 Suite (12 search modes + 3 utilities)
# ---------------------------------------------------------------------------
def render_fig1(snap: Path, out: Path, theme: dict[str, str], dark: bool) -> Path:
    style_rc(theme)
    std, by_kt = load_std_by_kt(snap)
    util = read_tsv(snap / "utilities.tsv")
    meta = read_meta(snap)

    y_rows: list[tuple[str, float, object]] = []
    y = 0.0
    for gi, (gtitle, cids) in enumerate(GROUPS_FIG1):
        if gi > 0:
            y -= 0.60
        y_rows.append(("hdr", y, gtitle))
        y -= 0.92
        for cid in cids:
            y_rows.append(("std", y, cid))
            y -= 1.15
    y -= 0.75
    y_rows.append(("sep", y + 0.38, None))
    y_rows.append(("hdr", y, "Database & Domain Utilities (1 thread)"))
    y -= 0.92
    for ur in util:
        y_rows.append(("util", y, ur))
        y -= 1.15
    ymin, ymax = y + 0.30, 0.90

    fig = new_fig((13.4, 8.8), theme)
    add_header(
        fig,
        "LOSAT v0.2.0 outperforms NCBI BLAST+ 2.17.0 across all 12 search modes and 3 utilities",
        "16-thread wall time: 1.02\u00d7 (psiblast) to 18.3\u00d7 faster (629\u00d7 on tblastx)  \u00b7  Peak memory: 1.3\u00d7 to 25\u00d7 lower  \u00b7  Byte-identical output",
        theme,
        y_title=0.974,
        y_sub=0.944,
    )

    has_anon = any(r.get("losat_rss_anon_mib") and r.get("ncbi_rss_anon_mib") for r in std)
    handles = [
        Patch(facecolor=theme["losat"], label="Search modes (BLAST+ \u00f7 LOSAT)"),
        Patch(facecolor=theme["wasm"], label="Utilities, 1t (makeblastdb / rpsbproc)"),
        ncbi_ref_handle(theme),
        Patch(facecolor=theme["badge_bg"], edgecolor=theme["hairline"], label="[-db] paired DB  \u00b7  [-subj] LOSAT -subject vs BLAST+ -db"),
    ]
    if has_anon:
        handles.insert(1, Patch(facecolor=theme["losat_soft"], label="max RSS (includes DB mmap); solid = RssAnon"))
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.018, 0.906), frameon=False, ncol=4, fontsize=7.6, columnspacing=1.1)

    bot_y, top_h = 0.105, 0.725
    ax_lbl = fig.add_axes([0.020, bot_y, 0.228, top_h])
    ax_lbl.set_facecolor(theme["surface"])
    ax_lbl.set_xlim(0, 1)
    ax_lbl.set_ylim(ymin, ymax)
    ax_lbl.axis("off")

    panels = [
        ([0.258, bot_y, 0.220, top_h], "a", "Speedup, 16 threads", 16, "speed", 28.0, [1, 2, 5, 10, 20], "BLAST+ \u00f7 LOSAT wall time"),
        ([0.506, bot_y, 0.214, top_h], "b", "Speedup, 1 thread", 1, "speed", 12.0, [1, 2, 5, 10], "BLAST+ \u00f7 LOSAT wall time"),
        ([0.748, bot_y, 0.230, top_h], "c", "Peak memory reduction, 16t", 16, "mem", 32.0, [1, 2, 5, 10, 25], "BLAST+ \u00f7 LOSAT max RSS"),
    ]

    row_idx = 0
    for kind, yv, item in y_rows:
        if kind in ("std", "util"):
            if row_idx % 2 == 1:
                ax_lbl.axhspan(yv - 0.52, yv + 0.52, color=theme["band"], zorder=0)
            row_idx += 1
        elif kind == "sep":
            ax_lbl.axhline(yv, color=theme["hairline"], linewidth=0.8, linestyle="--", zorder=1)
        elif kind == "hdr":
            ax_lbl.text(0.0, yv, str(item).upper(), ha="left", va="center", fontsize=7.6, fontweight="bold", color=theme["losat"])

    for kind, yv, item in y_rows:
        if kind == "std":
            r = by_kt[(str(item), 16)]
            draw_row_badge(ax_lbl, yv, r["label"], r["short_desc"], "-db" if r["target"] == "-db" else "-subj", theme, title_fs=8.8)
        elif kind == "util":
            r = item
            draw_row_badge(ax_lbl, yv, r["label"], r["short_desc"], "1t", theme, title_fs=8.7)

    for p_idx, (rect, let, title, th, metric, cap, ticks, xlabel) in enumerate(panels):
        ax = fig.add_axes(rect)
        ax.set_facecolor(theme["surface"])
        ax.set_xscale("log")
        xmax = cap * (7.8 if metric == "mem" else 7.0)
        ax.set_xlim(0.85, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_yticks([])

        ax.axvspan(1.0, cap, color=theme["zone"], alpha=0.55, zorder=0)
        row_idx = 0
        for kind, yv, _item in y_rows:
            if kind in ("std", "util"):
                if row_idx % 2 == 1:
                    ax.axhspan(yv - 0.52, yv + 0.52, color=theme["band"], alpha=0.65, zorder=0)
                row_idx += 1
            elif kind == "sep":
                ax.axhline(yv, color=theme["hairline"], linewidth=0.8, linestyle="--", zorder=1)

        ax.axvline(1.0, color=theme["ref"], linewidth=1.15, zorder=2)
        style_axis_log_x(
            ax, ticks, theme, labels=[f"{t}\u00d7" for t in ticks], tick_fs=7.8,
            skip_one=True, bound_lo=0.85, bound_hi=cap, major_len=3.5,
        )
        ax.set_xlabel(xlabel, fontsize=7.8, color=theme["ink2"], labelpad=4)
        add_panel_label(ax, let, title, theme, y=1.018)

        for kind, yv, item in y_rows:
            if kind == "std":
                r = by_kt[(str(item), th)]
                if metric == "speed":
                    draw_bar_with_labels(ax, yv, float(r["speedup"]), cap, theme["losat"], pair_sec(r), theme, bold_at=10.0)
                else:
                    val_anon = (
                        float(r["ncbi_rss_anon_mib"]) / float(r["losat_rss_anon_mib"])
                        if r.get("losat_rss_anon_mib") and r.get("ncbi_rss_anon_mib")
                        else None
                    )
                    draw_bar_with_labels(
                        ax, yv, float(r["mem_reduction"]), cap, theme["losat"], pair_mem(r), theme,
                        bold_at=9.5, val_anon=val_anon, soft_color=theme["losat_soft"],
                    )
            elif kind == "util":
                r = item
                if p_idx == 0:
                    ax.text(1.15, yv, "single-threaded utility (see panel b)", ha="left", va="center",
                            fontsize=7.2, style="italic", color=theme["muted"])
                elif p_idx == 1:
                    draw_bar_with_labels(ax, yv, float(r["speedup"]), cap, theme["wasm"], pair_sec(r), theme, bold_at=None)
                else:
                    draw_bar_with_labels(ax, yv, float(r["mem_reduction"]), cap, theme["wasm"], pair_mem(r), theme, bold_at=None)

    add_footer(
        fig,
        f"LOSAT {short_sha(meta.get('losat_commit', ''))} vs NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')} & rpsbproc {meta.get('rpsbproc_version', '0.5.1')}, {short_hw(meta)}. "
        f"Medians of 8 rounds at 16t (profile modes: 3; NCBI tblastx: 1), 5 LOSAT / 1 NCBI at 1t (utilities: 3).\n"
        f"Axis break (//) caps extreme outliers (tblastx 629\u00d7 / 104\u00d7) so typical ranges stay readable. "
        f"BLAST+ max RSS includes memory-mapped database pages; blastp/blastp-fast parity confirmed form-for-form (-subject and -db).",
        theme,
        y_rule=0.054,
    )
    return save_themed(fig, out, "fig1_suite_overview", dark, theme)


# ---------------------------------------------------------------------------
# Figure 2: Real-World Pipelines & Published Tools (>= 1 s/run ONLY in color)
# ---------------------------------------------------------------------------
def render_fig2(snap: Path, out: Path, theme: dict[str, str], dark: bool) -> Path:
    style_rc(theme)
    large = read_tsv(snap / "large_workflows.tsv")
    dt = read_tsv(snap / "downstream_tools.tsv")
    meta = read_meta(snap)

    real_mt = []
    real_st = []
    startup_tools: set[str] = set()
    for r in dt:
        if r["mt_speedup"]:
            if startup_dominated(r["ncbi_mt_s"], r["cases_count"]):
                startup_tools.add(r["tool"])
            else:
                real_mt.append(r)
        elif r["note"] == "losat_under_10ms":
            startup_tools.add(r["tool"])
        elif r["n1_speedup"]:
            if startup_dominated(r["ncbi_n1_s"], r["cases_count"]):
                startup_tools.add(r["tool"])
            else:
                real_st.append(r)

    real_mt.sort(key=lambda r: float(r["mt_speedup"]), reverse=True)
    real_st.sort(key=lambda r: float(r["n1_speedup"]), reverse=True)
    large_sorted = sorted(large, key=lambda r: float(r["speedup"]), reverse=True)

    fig = new_fig((14.0, 9.0), theme)
    add_header(
        fig,
        "Drop-in speed and memory gains across real bioinformatics workflows (excluding < 1 s start-up runs)",
        "Unmodified BLAST+ commands via losat --install-shims  \u00b7  Only workloads where BLAST+ took \u2265 1 s per run are plotted  \u00b7  All 28 tested tools matched NCBI output",
        theme,
        y_title=0.974,
        y_sub=0.945,
    )

    handles = [
        Patch(facecolor=theme["losat"], label="Wall-time speedup (BLAST+ \u00f7 LOSAT)"),
        Patch(facecolor=theme["wasm"], label="Peak memory reduction (BLAST+ \u00f7 LOSAT max RSS)"),
        ncbi_ref_handle(theme),
        Line2D([0], [0], linestyle="none", marker="", label="\u2020 NCBI multi-thread output drifts; LOSAT matches 1-thread NCBI"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.018, 0.908), frameon=False, ncol=4, fontsize=7.6, columnspacing=1.1, handlelength=1.4)

    def draw_dual_table(
        items: list[dict[str, object]],
        rect_lbl: list[float],
        rect_sp: list[float],
        rect_mem: list[float],
        letter: str,
        title: str,
        sp_cap: float,
        sp_ticks: list[int],
        mem_cap: float,
        mem_ticks: list[int],
    ) -> None:
        n = len(items)
        ys = [float(n - 1 - i) for i in range(n)]
        ymin_p, ymax_p = -0.65, n - 0.35

        ax_l = fig.add_axes(rect_lbl)
        ax_s = fig.add_axes(rect_sp)
        ax_m = fig.add_axes(rect_mem)
        for a in (ax_l, ax_s, ax_m):
            a.set_facecolor(theme["surface"])
            a.set_ylim(ymin_p, ymax_p)
        ax_l.set_xlim(0, 1)
        ax_l.axis("off")

        add_panel_label(ax_l, letter, title, theme, y=1.032)

        for idx, (ax, cap, ticks, col_hdr, xlabel, key_val, key_sub, bar_col) in enumerate([
            (ax_s, sp_cap, sp_ticks, "Speedup", "BLAST+ \u00f7 LOSAT time", "speedup", "time_txt", theme["losat"]),
            (ax_m, mem_cap, mem_ticks, "Memory", "BLAST+ \u00f7 LOSAT RSS", "mem", "mem_txt", theme["wasm"]),
        ]):
            ax.set_xscale("log")
            ax.set_xlim(0.88, cap * 5.8)
            ax.set_yticks([])
            ax.axvspan(1.0, cap, color=theme["zone"] if idx == 0 else theme["band"], alpha=0.55, zorder=0)
            for i, yv in enumerate(ys):
                if i % 2 == 1:
                    if idx == 0:
                        ax_l.axhspan(yv - 0.46, yv + 0.46, color=theme["band"], zorder=0)
                    ax.axhspan(yv - 0.46, yv + 0.46, color=theme["band"], alpha=0.6, zorder=0)
            ax.axvline(1.0, color=theme["ref"], linewidth=1.1, zorder=2)
            style_axis_log_x(
                ax, ticks, theme, labels=[f"{t}\u00d7" for t in ticks], tick_fs=7.5,
                skip_one=True, bound_lo=0.88, bound_hi=cap, major_len=3.0,
            )
            ax.set_xlabel(xlabel, fontsize=7.5, color=theme["ink2"], labelpad=3)
            ax.text(0.0, 1.032, col_hdr, transform=ax.transAxes, fontsize=8.0, fontweight="bold",
                    va="bottom", ha="left", color=bar_col)

            for yv, it in zip(ys, items):
                draw_bar_with_labels(
                    ax, yv, float(it[key_val]), cap, bar_col, str(it[key_sub]), theme,
                    val_fs=8.0, det_fs=6.5, dy_top=0.16, dy_bot=0.20, bold_at=9.5,
                )

        for yv, it in zip(ys, items):
            dagger = " \u2020" if it.get("drifts") else ""
            ax_l.text(0.01, yv + 0.16, f"{it['name']}{dagger}", ha="left", va="center", fontsize=8.3, fontweight="bold", color=theme["ink"])
            ax_l.text(0.01, yv - 0.20, str(it["sub"]), ha="left", va="center", fontsize=6.6, color=theme["muted"])

    def _build_item(name: str, sub: str, sp: float, ls: float, ns: float, mm: float, lm: float, nm: float, drift: bool) -> dict[str, object]:
        return {
            "name": name, "sub": sub, "speedup": sp,
            "time_txt": f"{fmt_sec(ls)} vs {fmt_sec(ns)}",
            "mem": mm, "mem_txt": f"{fmt_mem(lm)} vs {fmt_mem(nm)}",
            "drifts": drift,
        }

    items_a = [_build_item(
        r["label"], f"{r['program']} \u00b7 {r['threads']}t \u00b7 {r['detail']}",
        float(r["speedup"]), float(r["losat_wall_s"]), float(r["ncbi_wall_s"]),
        float(r["mem_reduction"]), float(r["losat_rss_mib"]), float(r["ncbi_rss_mib"]),
        r["ncbi_mt_differs"] == "yes",
    ) for r in large_sorted]

    items_c = [_build_item(
        r["tool"], f"{r['program']} \u00b7 1t \u00b7 {r['cases_count']} run(s)",
        float(r["n1_speedup"]), float(r["losat_n1_s"]), float(r["ncbi_n1_s"]),
        float(r["n1_mem_reduction"]), float(r["losat_n1_rss_mib"]), float(r["ncbi_n1_rss_mib"]),
        False,
    ) for r in real_st]

    items_b = [_build_item(
        f"{r['tool']} ({r['program']})" if r["tool"] in ("Kaptive", "Trinotate") else r["tool"],
        f"{r['program']} \u00b7 {r['threads_label']}t \u00b7 {r['cases_count']} run(s)",
        float(r["mt_speedup"]), float(r["losat_mt_s"]), float(r["ncbi_mt_s"]),
        float(r["mt_mem_reduction"]), float(r["losat_mt_rss_mib"]), float(r["ncbi_mt_rss_mib"]),
        r["ncbi_mt_differs"] == "yes",
    ) for r in real_mt]

    # Generous vertical gap between legend (y=0.908), top panels (top=0.825), and Panel c (top=0.395)
    draw_dual_table(
        items_a,
        [0.020, 0.495, 0.215, 0.335],
        [0.240, 0.495, 0.118, 0.335],
        [0.365, 0.495, 0.118, 0.335],
        "a", "Large reference workflows (16 threads)",
        sp_cap=12.0, sp_ticks=[1, 2, 5, 10],
        mem_cap=20.0, mem_ticks=[1, 2, 5, 15],
    )
    draw_dual_table(
        items_c,
        [0.020, 0.105, 0.215, 0.290],
        [0.240, 0.105, 0.118, 0.290],
        [0.365, 0.105, 0.118, 0.290],
        "c", "Single-thread published tools (\u2265 1 s/run)",
        sp_cap=8.0, sp_ticks=[1, 2, 4, 8],
        mem_cap=12.0, mem_ticks=[1, 2, 5, 10],
    )

    fig.add_artist(Line2D([0.502, 0.502], [0.095, 0.855], transform=fig.transFigure, color=theme["grid"], linewidth=0.9))

    draw_dual_table(
        items_b,
        [0.518, 0.105, 0.185, 0.725],
        [0.708, 0.105, 0.130, 0.725],
        [0.845, 0.105, 0.130, 0.725],
        "b", "Multi-thread published tools (\u2265 1 s/run)",
        sp_cap=22.0, sp_ticks=[1, 2, 5, 10, 20],
        mem_cap=12.0, mem_ticks=[1, 2, 5, 10],
    )

    dt_commits = meta.get("dt_commits", {})
    one_shas = short_sha(",".join(dt_commits[k] for k in sorted(dt_commits) if dt_commits[k]))
    excluded_list = ", ".join(sorted(startup_tools))
    add_footer(
        fig,
        f"Panels a, b: LOSAT {short_sha(meta.get('toolmt_commit', ''))} (medians of 3 rounds, peak max RSS across runs); Panel c: LOSAT {one_shas} (medians of 4 runs); "
        f"vs NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')}, {short_hw(meta)}.\n"
        f"Start-up-dominated official test suites (< 1 s/run in BLAST+: {excluded_list}) also produced identical outputs and are reported in snapshot/downstream_tools.tsv rather than plotted as search-speed bars.",
        theme,
        y_rule=0.054,
    )
    return save_themed(fig, out, "fig2_real_pipelines", dark, theme)


# ---------------------------------------------------------------------------
# Figure 3: Scaling (2x2: DB Memory, DB Warm/Cold Time, Query Volume, Threads)
# ---------------------------------------------------------------------------
def render_fig3(snap: Path, out: Path, theme: dict[str, str], dark: bool) -> Path:
    style_rc(theme)
    rows = read_tsv(snap / "scaling.tsv")
    meta = read_meta(snap)

    fig = new_fig((13.0, 8.2), theme)
    add_header(
        fig,
        "Memory stays under 410 MiB on 319 GB core_nt while scaling cleanly across queries and threads",
        "Database size (0.16 MB \u2192 319.4 GB), query volume (10k \u2192 1M reads), and thread scaling (1 \u2192 16 cores)  \u00b7  All runs byte-identical to NCBI BLAST+ 2.17.0",
        theme,
        y_title=0.974,
        y_sub=0.942,
    )

    handles = [
        Line2D([0], [0], marker="o", color=theme["losat"], markersize=6.2, linewidth=2.0, label="LOSAT 0.2.0 (warm cache)"),
        Line2D([0], [0], marker="D", color=theme["ncbi"], markersize=5.4, linewidth=1.8, label="NCBI BLAST+ 2.17.0 (warm cache)"),
        Line2D([0], [0], marker="o", markerfacecolor=theme["surface"], markeredgecolor=theme["losat"], markersize=6.0, linestyle="none", markeredgewidth=1.5, label="Cold cache (pages dropped)"),
        Line2D([0], [0], linestyle=":", color=theme["muted"], linewidth=1.2, label="Ideal linear thread scaling"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.018, 0.902), frameon=False, ncol=4, fontsize=7.9, columnspacing=1.3)

    db_rows = [r for r in rows if r["panel"] == "db_size"]
    warm = sorted([r for r in db_rows if r["cache"] == "warm"], key=lambda r: float(r["x_value"]))
    cold = sorted([r for r in db_rows if r["cache"] == "cold"], key=lambda r: float(r["x_value"]))
    db_short = {
        "abr_genes": "abr_genes\n0.16 MB",
        "16s_all": "16S rRNA\n15 MB",
        "ref_prok_rep_genomes": "ref_prok\n27 GB",
        "core_nt": "core_nt\n319 GB",
    }
    cold_arrow = dict(arrowstyle="-|>", color=theme["muted"], lw=0.8, shrinkA=0, shrinkB=4)
    pill_box = dict(boxstyle="round,pad=0.12", facecolor=theme["surface"], edgecolor="none", alpha=0.85)

    # --- Panel a: Database size -> Peak Memory ---
    ax_a = fig.add_axes([0.072, 0.540, 0.402, 0.295])
    ax_a.set_facecolor(theme["surface"])
    xs_w = [float(r["x_value"]) for r in warm]
    yn_m = [float(r["ncbi_rss_mib"]) for r in warm]
    yl_m = [float(r["losat_rss_mib"]) for r in warm]
    ax_a.fill_between(xs_w, yl_m, yn_m, color=theme["losat_fill"], alpha=0.55, zorder=1)
    ax_a.plot(xs_w, yn_m, marker="D", markersize=5.8, color=theme["ncbi"], linewidth=2.0, zorder=3)
    ax_a.plot(xs_w, yl_m, marker="o", markersize=6.5, color=theme["losat"], linewidth=2.2, zorder=4)

    for i, (_r, x, yn, yl) in enumerate(zip(warm, xs_w, yn_m, yl_m)):
        ratio = yn / yl
        if i in (2, 3):
            ax_a.text(x, math.sqrt(yn * yl), fmt_fold(ratio), ha="center", va="center", fontsize=8.4,
                      fontweight="bold", color=theme["losat"],
                      bbox=dict(boxstyle="round,pad=0.18", facecolor=theme["surface"], edgecolor=theme["zone_edge"], linewidth=0.7), zorder=6)
        else:
            ax_a.text(x, yl * 0.54, fmt_fold(ratio), ha="center", va="top", fontsize=7.6, color=theme["ink2"], zorder=6)
    ax_a.text(xs_w[-1], yn_m[-1] * 1.55, f"{fmt_mem(yn_m[-1])}", ha="center", va="bottom", fontsize=8.1, fontweight="bold", color=theme["ncbi"])
    ax_a.text(xs_w[-1], yl_m[-1] * 0.52, f"{fmt_mem(yl_m[-1])}", ha="center", va="top", fontsize=8.2, fontweight="bold", color=theme["losat"])

    ax_a.set_xscale("log")
    ax_a.set_yscale("log")
    ax_a.set_xlim(0.04, 2.2e6)
    ax_a.set_ylim(30, 2.2e6)
    ax_a.set_xticks(xs_w)
    ax_a.set_xticklabels([db_short[r["x_label"]] for r in warm], fontsize=7.5, color=theme["ink2"])
    ax_a.set_yticks([100, 1000, 10000, 100000])
    ax_a.set_yticklabels(["100 MiB", "1 GiB", "10 GiB", "100 GiB"], fontsize=7.6, color=theme["ink2"])
    style_clean_grid(ax_a, theme, lw=0.65)
    ax_a.set_ylabel("Peak memory (max RSS)", fontsize=8.3, color=theme["ink2"])
    add_panel_label(ax_a, "a", "Database size \u2192 peak memory", theme, y=1.04, sub="megablast \u00b7 210 contigs \u00b7 16t \u00b7 -db")

    # --- Panel b: Database size -> Wall time (Warm & Cold cache) ---
    ax_b = fig.add_axes([0.562, 0.540, 0.410, 0.295])
    ax_b.set_facecolor(theme["surface"])
    yn_w = [float(r["ncbi_wall_s"]) for r in warm]
    yl_w = [float(r["losat_wall_s"]) for r in warm]
    ax_b.fill_between(xs_w, yl_w, yn_w, color=theme["losat_fill"], alpha=0.45, zorder=1)
    ax_b.plot(xs_w, yn_w, marker="D", markersize=5.8, color=theme["ncbi"], linewidth=2.0, zorder=3)
    ax_b.plot(xs_w, yl_w, marker="o", markersize=6.5, color=theme["losat"], linewidth=2.2, zorder=4)

    warm_by_lbl = {r["x_label"]: r for r in warm}
    for rc in cold:
        rw = warm_by_lbl[rc["x_label"]]
        x = float(rc["x_value"])
        ax_b.plot([x, x], [float(rw["ncbi_wall_s"]), float(rc["ncbi_wall_s"])], color=theme["ncbi"], linewidth=1.2, linestyle="--", zorder=2)
        ax_b.scatter([x], [float(rc["ncbi_wall_s"])], marker="D", s=40, facecolors=theme["surface"], edgecolors=theme["ncbi"], linewidths=1.5, zorder=5)
        ax_b.scatter([x], [float(rc["losat_wall_s"])], marker="o", s=44, facecolors=theme["surface"], edgecolors=theme["losat"], linewidths=1.5, zorder=5)

    cnt_c = next(r for r in cold if r["x_label"] == "core_nt")
    rpr_c = next(r for r in cold if r["x_label"] == "ref_prok_rep_genomes")
    ax_b.annotate(
        f"cold: {fmt_fold(float(rpr_c['speedup']))}",
        xy=(float(rpr_c["x_value"]), float(rpr_c["ncbi_wall_s"])),
        xytext=(float(rpr_c["x_value"]) * 0.11, float(rpr_c["ncbi_wall_s"]) * 2.1),
        fontsize=7.4, color=theme["ink2"], arrowprops=cold_arrow,
    )
    ax_b.annotate(
        f"cold: {fmt_fold(float(cnt_c['speedup']))}\n(776s vs 491s)",
        xy=(float(cnt_c["x_value"]), float(cnt_c["ncbi_wall_s"])),
        xytext=(float(cnt_c["x_value"]) * 0.055, float(cnt_c["ncbi_wall_s"]) * 1.15),
        fontsize=7.3, color=theme["ink2"], ha="center", arrowprops=cold_arrow,
    )
    for i, (x, yn, yl) in enumerate(zip(xs_w, yn_w, yl_w)):
        ratio = yn / yl
        if i == 0:
            ax_b.text(x, math.sqrt(yn * yl), f"{fmt_fold(ratio)}*", ha="center", va="center", fontsize=7.5, color=theme["muted"], bbox=pill_box)
        elif i == 1:
            ax_b.text(x, yl * 0.48, fmt_fold(ratio), ha="center", va="top", fontsize=7.8, fontweight="semibold", color=theme["ink"])
        elif i == 2:
            ax_b.text(x, yl * 0.45, f"warm {fmt_fold(ratio)}", ha="center", va="top", fontsize=7.5, fontweight="semibold", color=theme["losat"])
        else:
            ax_b.text(x, yl * 0.34, f"warm {fmt_fold(ratio)}", ha="center", va="top", fontsize=7.5, fontweight="semibold", color=theme["losat"])

    ax_b.set_xscale("log")
    ax_b.set_yscale("log")
    ax_b.set_xlim(0.04, 2.2e6)
    ax_b.set_ylim(0.04, 2800)
    ax_b.set_xticks(xs_w)
    ax_b.set_xticklabels([db_short[r["x_label"]] for r in warm], fontsize=7.5, color=theme["ink2"])
    style_clean_grid(ax_b, theme, lw=0.65, plain_y=True)
    ax_b.set_ylabel("Wall time (s)", fontsize=8.3, color=theme["ink2"])
    add_panel_label(ax_b, "b", "Database size \u2192 wall time (warm & cold)", theme, y=1.04, sub="*abr_genes: 0.90 s start-up")

    # --- Panel c: Query volume (10k -> 1M reads): Twin mini-panels (Time + Peak Memory) ---
    q_rows = sorted([r for r in rows if r["panel"] == "query_size"], key=lambda r: float(r["x_value"]))
    xs_q = [float(r["x_value"]) for r in q_rows]
    ax_c1 = fig.add_axes([0.072, 0.115, 0.178, 0.250])
    ax_c2 = fig.add_axes([0.296, 0.115, 0.178, 0.250])
    for ax, key_n, key_l, ylab, sub_lbl in [
        (ax_c1, "ncbi_wall_s", "losat_wall_s", "Wall time (s)", "Wall time"),
        (ax_c2, "ncbi_rss_mib", "losat_rss_mib", "Max RSS (MiB)", "Peak memory"),
    ]:
        ax.set_facecolor(theme["surface"])
        yn = [float(r[key_n]) for r in q_rows]
        yl = [float(r[key_l]) for r in q_rows]
        ax.fill_between(xs_q, yl, yn, color=theme["losat_fill"], alpha=0.48, zorder=1)
        ax.plot(xs_q, yn, marker="D", markersize=5.2, color=theme["ncbi"], linewidth=1.8, zorder=3)
        ax.plot(xs_q, yl, marker="o", markersize=5.8, color=theme["losat"], linewidth=2.0, zorder=4)
        for x, n_v, l_v in zip(xs_q, yn, yl):
            ax.text(x, math.sqrt(n_v * l_v), fmt_fold(n_v / l_v), ha="center", va="center", fontsize=7.6,
                    fontweight="semibold", color=theme["ink"], bbox=pill_box)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(5e3, 2.2e6)
        ax.set_ylim(min(yl) * 0.55, max(yn) * 2.2)
        ax.set_xticks([1e4, 1e5, 1e6])
        ax.set_xticklabels(["10k", "100k", "1M"], fontsize=7.6, color=theme["ink2"])
        style_clean_grid(ax, theme, lw=0.6, plain_y=True)
        ax.set_xlabel("Query reads (150 nt)", fontsize=7.8, color=theme["ink2"])
        ax.set_ylabel(ylab, fontsize=8.0, color=theme["ink2"])
        ax.set_title(sub_lbl, fontsize=8.1, fontweight="semibold", color=theme["ink2"], pad=6)
    fig.text(0.072, 0.422, "c  Query volume \u2192 time & memory", fontsize=9.6, fontweight="bold", color=theme["ink"])
    fig.text(0.474, 0.422, "megablast vs 16S \u00b7 16t", fontsize=7.5, ha="right", color=theme["muted"])

    # --- Panel d: Thread scaling (1 -> 16 cores, 3 clean sub-strips with crisp clearance above subtitles) ---
    t_rows = [r for r in rows if r["panel"] == "threads"]
    t_specs = [
        ("bn", "blastn", "11.8\u00d7 / 16t"),
        ("bp", "blastp", "14.9\u00d7 / 16t"),
        ("tbn", "tblastn", "13.6\u00d7 / 16t"),
    ]
    fig.text(0.562, 0.422, "d  Thread scaling (1 \u2192 16 cores)", fontsize=9.6, fontweight="bold", color=theme["ink"])
    fig.text(0.972, 0.422, "3.1\u20134.9\u00d7 faster at every thread count", fontsize=7.5, ha="right", color=theme["muted"])
    for idx, (cid, mode_name, eff_txt) in enumerate(t_specs):
        ax = fig.add_axes([0.562 + idx * 0.144, 0.115, 0.122, 0.250])
        ax.set_facecolor(theme["surface"])
        sub = sorted([r for r in t_rows if r["case"] == cid], key=lambda x: int(x["x_value"]))
        xs = [int(r["x_value"]) for r in sub]
        ys_n = [float(r["ncbi_wall_s"]) for r in sub]
        ys_l = [float(r["losat_wall_s"]) for r in sub]
        ax.fill_between(xs, ys_l, ys_n, color=theme["losat_fill"], alpha=0.45, zorder=1)
        ax.plot(xs, [ys_n[0] / x for x in xs], linestyle=":", color=theme["ncbi"], alpha=0.55, linewidth=1.0)
        ax.plot(xs, [ys_l[0] / x for x in xs], linestyle=":", color=theme["losat"], alpha=0.65, linewidth=1.0)
        ax.plot(xs, ys_n, marker="D", markersize=4.4, color=theme["ncbi"], linewidth=1.6, zorder=3)
        ax.plot(xs, ys_l, marker="o", markersize=5.0, color=theme["losat"], linewidth=1.8, zorder=4)
        ax.text(1.75, math.sqrt(ys_n[1] * ys_l[0]) * 0.92, fmt_fold(ys_n[0] / ys_l[0]), ha="center", va="center",
                fontsize=7.3, fontweight="semibold", color=theme["ink"], bbox=pill_box)
        ax.text(9.8, math.sqrt(ys_n[-1] * ys_l[-2]) * 0.98, fmt_fold(ys_n[-1] / ys_l[-1]), ha="center", va="center",
                fontsize=7.3, fontweight="bold", color=theme["losat"], bbox=pill_box)
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_xticks([1, 4, 16])
        ax.set_xticklabels(["1", "4", "16"], fontsize=7.5, color=theme["ink2"])
        ax.set_xlim(0.75, 22)
        y_lo, y_hi = min(ys_l) * 0.45, max(ys_n) * 1.9
        ax.set_ylim(y_lo, y_hi)
        yt = [c for c in [0.5, 1, 2, 5, 10, 20, 50, 100, 200] if y_lo <= c <= y_hi]
        ax.set_yticks(yt)
        ax.set_yticklabels([plain_log(c) for c in yt], fontsize=7.3, color=theme["ink2"])
        style_clean_grid(ax, theme, lw=0.55)
        ax.set_xlabel("Threads", fontsize=7.6, color=theme["ink2"])
        if idx == 0:
            ax.set_ylabel("Wall time (s)", fontsize=8.0, color=theme["ink2"])
        ax.set_title(f"{mode_name}  ({eff_txt})", fontsize=7.7, fontweight="semibold", color=theme["ink2"], pad=6)

    add_footer(
        fig,
        f"LOSAT {short_sha(meta.get('losat_commit', ''))} (core_nt & cold-cache runs: {short_sha(meta.get('toolmt_commit', ''))}) vs NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')}, {short_hw(meta)}.\n"
        f"Panels a\u2013c: both programs ran -db. Panel d: blastn/blastp use Table 1 forms (LOSAT -subject, BLAST+ -db); tblastn uses -db. BLAST+ max RSS includes memory-mapped database pages.",
        theme,
        y_rule=0.054,
    )
    return save_themed(fig, out, "fig3_scaling", dark, theme)


# ---------------------------------------------------------------------------
# Figure 4: Exact BLAST+ Parity vs Heuristic Aligners (UniRef50 & SILVA 144)
# ---------------------------------------------------------------------------
def render_fig4(snap: Path, out: Path, theme: dict[str, str], dark: bool) -> Path:
    style_rc(theme)
    ur50 = read_tsv(snap / "vs_heuristics_uniref50.tsv")
    silva = read_tsv(snap / "vs_heuristics_silva.tsv")
    meta = read_meta(snap)

    fig = new_fig((14.0, 8.2), theme)
    add_header(
        fig,
        "Exact NCBI BLAST+ parity vs. heuristic aligners: speed, sensitivity, and memory trade-offs",
        "Protein search on UniRef50 (38.8M seqs, 32 threads) vs. DIAMOND & MMseqs2  \u00b7  16S rRNA search on SILVA 144 (905k seqs, 16 threads) vs. vsearch",
        theme,
        y_title=0.974,
        y_sub=0.943,
    )

    tool_color = {
        "LOSAT": theme["losat"],
        "NCBI BLAST+": theme["ncbi"],
        "DIAMOND": theme["diamond"],
        "MMseqs2": theme["mmseqs"],
        "vsearch": theme["vsearch"],
    }
    handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=7.2, color=theme["losat"], label="LOSAT (exact NCBI output)"),
        Line2D([0], [0], marker="D", linestyle="none", markersize=6.2, color=theme["ncbi"], label="NCBI BLAST+ 2.17.0"),
        Line2D([0], [0], marker="s", linestyle="none", markersize=6.2, color=theme["diamond"], label="DIAMOND 2.2.8"),
        Line2D([0], [0], marker="^", linestyle="none", markersize=6.8, color=theme["mmseqs"], label="MMseqs2 18.8cc5c (filled: 101 GiB idx; open: no idx)"),
        Line2D([0], [0], marker="P", linestyle="none", markersize=6.8, color=theme["vsearch"], label="vsearch 2.32.0 (QIIME 2 default)"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.018, 0.908), frameon=False, ncol=3, fontsize=7.6, columnspacing=1.2)

    # Explicit non-overlapping absolute (tx, ty) callout coordinates with crisp leader lines!
    # Panel a: xlim=(8, 3200), ylim=(56, 105.5)
    callouts_q100 = {
        "losat_bp":         (340,  96.4, "right", "center", "LOSAT blastp\n407s \u00b7 0.13 GiB", True),
        "ncbi_def":         (1250, 96.4, "center", "center", "NCBI blastp\n1,068s \u00b7 13.4 GiB", True),
        "dmnd_us":          (1350, 88.2, "center", "center", "DIAMOND ultra\n836s \u00b7 3.6 GiB", True),
        "mmseqs_s75_idx":   (46,   96.0, "center", "center", "MMseqs2 -s7.5 idx\n94s \u00b7 199 GiB", True),
        "mmseqs_s75_noidx": (560,  81.6, "center", "center", "MMseqs2 -s7.5\n266s \u00b7 116 GiB", True),
        "losat_fast":       (10.0, 85.8, "left",   "center", "LOSAT blastp-fast\n28s \u00b7 0.15 GiB", True),
        "ncbi_fast":        (56,   76.0, "right",  "center", "NCBI blastp-fast\n92s \u00b7 13.4 GiB", True),
        "dmnd_ms":          (175,  71.2, "center", "center", "DIAMOND more\n163s \u00b7 9.0 GiB", True),
        "mmseqs_s57_idx":   (10.5, 78.8, "left",   "center", "MMseqs2 -s5.7 idx\n14s \u00b7 95 GiB", True),
        "mmseqs_s57_noidx": (650,  71.2, "center", "center", "MMseqs2 -s5.7\n245s \u00b7 95 GiB", True),
        "dmnd_def":         (62,   64.7, "left",   "center", "DIAMOND default\n40s \u00b7 5.0 GiB", False),
    }
    # Panel b: xlim=(16, 28000), ylim=(56, 105.5)
    callouts_q1000 = {
        "losat_bp":         (1650,  97.2, "right",  "center", "LOSAT blastp\n3,874s \u00b7 0.34 GiB", True),
        "ncbi_def":         (11500, 95.8, "center", "center", "NCBI blastp\n10,692s \u00b7 16.0 GiB", True),
        "dmnd_us":          (3200,  89.4, "left",   "center", "DIAMOND ultra\n1,467s \u00b7 5.6 GiB", True),
        "losat_fast":       (260,   94.4, "center", "center", "LOSAT blastp-fast\n214s \u00b7 0.36 GiB", True),
        "mmseqs_s75_idx":   (68,    92.5, "center", "center", "MMseqs2 -s7.5 idx\n171s \u00b7 200 GiB", True),
        "dmnd_ms":          (92,    74.2, "center", "center", "DIAMOND more\n221s \u00b7 7.6 GiB", True),
        "mmseqs_s75_noidx": (620,   72.2, "center", "center", "MMseqs2 -s7.5\n349s \u00b7 117 GiB", True),
        "ncbi_fast":        (3200,  81.2, "left",   "center", "NCBI fast\n860s \u00b7 16.0 GiB", True),
        "mmseqs_s57_idx":   (21,    77.0, "left",   "center", "MMseqs2 -s5.7 idx\n30s \u00b7 103 GiB", True),
        "mmseqs_s57_noidx": (680,   64.0, "center", "center", "MMseqs2 -s5.7\n260s \u00b7 104 GiB", True),
        "dmnd_def":         (82,    59.2, "left",   "center", "DIAMOND default\n52s \u00b7 5.7 GiB", False),
    }

    for p_i, (qset, rect, let, title, callouts, xlim) in enumerate([
        ("q100",  [0.058, 0.135, 0.285, 0.650], "a", "UniRef50: 100 queries (32t)",   callouts_q100,  (8.0, 3200.0)),
        ("q1000", [0.382, 0.135, 0.285, 0.650], "b", "UniRef50: 1,000 queries (32t)", callouts_q1000, (16.0, 28000.0)),
    ]):
        ax = fig.add_axes(rect)
        ax.set_facecolor(theme["surface"])
        sub = [r for r in ur50 if r["query_set"] == qset]

        ax.axhspan(98.6, 103.5, color=theme["zone"], alpha=0.65, zorder=0)
        ax.axhline(100.0, color=theme["losat"], linewidth=0.9, linestyle="--", alpha=0.7, zorder=1)
        ax.text(xlim[0] * 1.12, 100.9, "100% exact NCBI BLAST+ parity", fontsize=7.1, fontweight="semibold", color=theme["losat"], va="bottom")

        by_arm = {r["arm_id"]: r for r in sub}
        for n_id, l_id in (("ncbi_def", "losat_bp"), ("ncbi_fast", "losat_fast")):
            rn, rl = by_arm[n_id], by_arm[l_id]
            ax.annotate(
                "", xy=(float(rl["wall_s"]), float(rl["recall_uncapped"]) * 100),
                xytext=(float(rn["wall_s"]), float(rn["recall_uncapped"]) * 100),
                arrowprops=dict(arrowstyle="-|>", color=theme["losat"], lw=1.4, alpha=0.75, shrinkA=6, shrinkB=7),
                zorder=2,
            )

        for r in sub:
            arm = r["arm_id"]
            x = float(r["wall_s"])
            y = float(r["recall_uncapped"]) * 100.0
            col = tool_color[r["tool"]]
            if r["tool"] == "LOSAT":
                mk, sz, fc, ec = "o", 74, col, col
            elif r["tool"] == "NCBI BLAST+":
                mk, sz, fc, ec = "D", 54, col, col
            elif r["tool"] == "DIAMOND":
                mk, sz, fc, ec = "s", 52, col, col
            else:
                mk, sz = "^", 62
                fc = col if r["indexed"] == "yes" else theme["surface"]
                ec = col
            ax.scatter([x], [y], marker=mk, s=sz, facecolors=fc, edgecolors=ec, linewidths=1.5, zorder=5)

            tx, ty, ha, va, txt, draw_arrow = callouts[arm]
            is_losat = r["tool"] == "LOSAT"
            ax.annotate(
                txt,
                xy=(x, y),
                xytext=(tx, ty),
                ha=ha, va=va,
                fontsize=6.7 if not is_losat else 7.1,
                fontweight="bold" if is_losat else "normal",
                color=theme["losat"] if is_losat else theme["ink2"],
                linespacing=1.18,
                bbox=dict(boxstyle="round,pad=0.14", facecolor=theme["surface"], edgecolor=theme["zone_edge"] if is_losat else "none", linewidth=0.6 if is_losat else 0.0, alpha=0.90),
                arrowprops=dict(arrowstyle="-", color=theme["losat"] if is_losat else theme["hairline"], lw=0.75, shrinkA=2, shrinkB=5) if draw_arrow else None,
                zorder=7 if is_losat else 6,
            )

        ax.set_xscale("log")
        ax.set_xlim(*xlim)
        ax.set_ylim(55.5, 104.2)
        ax.xaxis.set_major_formatter(FuncFormatter(plain_log))
        ax.set_xlabel("\u2190 Faster   Wall time (s, log scale)", fontsize=8.0, color=theme["ink2"])
        if p_i == 0:
            ax.set_ylabel("Recall of NCBI blastp hits (%, uncapped queries)  \u2191", fontsize=8.1, color=theme["ink2"])
        else:
            ax.set_yticklabels([])
        style_clean_grid(ax, theme, lw=0.6)
        add_panel_label(ax, let, title, theme, y=1.022)

    # --- Panel c: SILVA 144 (16S rRNA, 905k seqs, 16 threads) ---
    fig.add_artist(Line2D([0.692, 0.692], [0.105, 0.820], transform=fig.transFigure, color=theme["grid"], linewidth=0.9))
    fig.text(0.712, 0.799, "c  SILVA 144 (905k rRNA, 16t)", fontsize=9.6, fontweight="bold", color=theme["ink"])

    silva_cases = [
        ("reads1e6", "1M reads\n(megablast)"),
        ("amp1e3",   "1k amplicons\n(QIIME 2)"),
        ("amp1e4",   "10k amplicons\n(QIIME 2)*"),
    ]
    by_sc = {(r["case"], r["arm"]): r for r in silva}

    ax_c1 = fig.add_axes([0.752, 0.485, 0.226, 0.275])
    ax_c2 = fig.add_axes([0.752, 0.135, 0.226, 0.265])

    for panel_idx, (ax, metric, ylab, ymax_c) in enumerate([
        (ax_c1, "wall_s",      "Wall time (s)",   5500.0),
        (ax_c2, "max_rss_mib", "Memory (MiB)",   50000.0),
    ]):
        ax.set_facecolor(theme["surface"])
        xs = [0.0, 1.0, 2.0]
        w = 0.25
        for ci, (cid, _lbl) in enumerate(silva_cases):
            arms_here = ["losat", "ncbi", "vsearch"] if cid == "amp1e3" else ["losat", "ncbi"]
            offs = [-w * 1.05, 0.0, w * 1.05] if len(arms_here) == 3 else [-w * 0.58, w * 0.58]
            for off, arm in zip(offs, arms_here):
                r = by_sc[(cid, arm)]
                v = float(r[metric])
                col = tool_color[r["tool"]]
                if panel_idx == 1:
                    va = float(r["rss_anon_mib"])
                    ax.bar([ci + off], [v], width=w * 0.85, color=col, alpha=0.32, edgecolor=col, linewidth=0.9, zorder=3)
                    ax.bar([ci + off], [va], width=w * 0.85, color=col, alpha=0.95, edgecolor="none", zorder=4)
                else:
                    ax.bar([ci + off], [v], width=w * 0.85, color=col, alpha=0.92, edgecolor="none", zorder=3)
                txt = fmt_sec(v) if panel_idx == 0 else fmt_mem(v, compact=True)
                ax.text(ci + off, v * 1.18, txt, ha="center", va="bottom", fontsize=6.6,
                        fontweight="bold" if arm == "losat" else "normal",
                        color=theme["losat"] if arm == "losat" else theme["ink2"])
        ax.set_yscale("log")
        ax.set_xlim(-0.55, 2.55)
        ax.set_ylim(10.0, ymax_c)
        ax.set_xticks(xs)
        ax.set_xticklabels([lbl for _, lbl in silva_cases], fontsize=7.2, color=theme["ink2"])
        ax.set_ylabel(ylab, fontsize=7.9, color=theme["ink2"])
        if panel_idx == 1:
            ax.set_title("Solid = RssAnon (heap) \u00b7 Tint = max RSS", fontsize=7.0, color=theme["muted"], pad=5)
        else:
            ax.set_title("Wall time (shorter is better)", fontsize=7.4, color=theme["muted"], pad=5)
        style_clean_grid(ax, theme, lw=0.6, axis="y", plain_y=True, clear_minor_x=False)

    add_footer(
        fig,
        f"UniRef50 (38.8M seqs, 32t, node 1): LOSAT {short_sha(meta.get('uniref50_losat_commit', ''))} vs NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')}, "
        f"DIAMOND {meta.get('diamond_version', '2.2.8')}, MMseqs2 {meta.get('mmseqs2_version', '18.8cc5c')} (all Bioconda; index = 101 GiB).\n"
        f"SILVA 144 NR99 (905k seqs, 16t): LOSAT {short_sha(meta.get('silva_losat_commit', ''))} vs NCBI BLAST+ 2.17.0 & vsearch {meta.get('vsearch_version', '2.32.0')} (QIIME 2 --maxrejects 0; 56% hit overlap with BLAST+). "
        f"*On 10k amplicons, NCBI 16t drifts (+5 rows); LOSAT matches 1-thread NCBI byte-for-byte.",
        theme,
        y_rule=0.060,
    )
    return save_themed(fig, out, "fig4_vs_heuristics", dark, theme)


# ---------------------------------------------------------------------------
# Figure 5: WebAssembly Portability (Speed vs Native BLAST+ & Peak Memory)
# ---------------------------------------------------------------------------
def render_fig5(snap: Path, out: Path, theme: dict[str, str], dark: bool) -> Path:
    style_rc(theme)
    rows = read_tsv(snap / "webassembly.tsv")
    by_kt = {(r["case"], int(r["threads"])): r for r in rows}
    meta = read_meta(snap)

    wasm_derived = bool(meta.get("wasm_derived", False))
    base_sha = short_sha(meta.get("wasm_derived_from") or "")
    if wasm_derived and not base_sha:
        sys.exit("metadata.json has wasm_derived=true without wasm_derived_from")

    y_rows: list[tuple[str, float, str | None]] = []
    y = 0.0
    for gi, (gtitle, cids) in enumerate(WASM_GROUPS):
        if gi > 0:
            y -= 0.55
        y_rows.append(("hdr", y, gtitle))
        y -= 0.90
        for cid in cids:
            y_rows.append(("case", y, cid))
            y -= 1.05
    ymin, ymax = y + 0.30, 0.85

    fig = new_fig((13.4, 7.0), theme)
    sub_extra = f"  \u00b7  Estimated from build {base_sha} via side-by-side ratio" if wasm_derived else "  \u00b7  Direct paired runs under Node.js"
    add_header(
        fig,
        "Even inside a WebAssembly sandbox, LOSAT beats native C++ NCBI BLAST+ 2.17.0 in 8 of 8 measured modes at 16 threads",
        f"Native LOSAT vs. wasm32-wasip1-threads (Node.js) vs. native NCBI BLAST+ 2.17.0{sub_extra}",
        theme,
        y_title=0.970,
        y_sub=0.934,
    )

    handles = [
        Line2D([0], [0], marker="o", linestyle="none", markersize=6.8, color=theme["losat"], label="LOSAT native"),
        Line2D([0], [0], marker="^", linestyle="none", markersize=7.0, color=theme["wasm"], label="LOSAT WebAssembly (WASI threads)"),
        Line2D([0], [0], marker="D", linestyle="none", markersize=6.0, color=theme["ncbi"], label="Native NCBI BLAST+ 2.17.0"),
        Line2D([0], [0], color=theme["ref"], linewidth=1.2, label="Native BLAST+ = 1\u00d7 (shaded: faster than native BLAST+)"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.018, 0.888), frameon=False, ncol=4, fontsize=7.8, columnspacing=1.2)

    bot_y, h_ax = 0.130, 0.665
    ax_lbl = fig.add_axes([0.020, bot_y, 0.138, h_ax])
    ax_lbl.set_facecolor(theme["surface"])
    ax_lbl.set_xlim(0, 1)
    ax_lbl.set_ylim(ymin, ymax)
    ax_lbl.axis("off")

    row_i = 0
    for kind, yv, txt in y_rows:
        if kind == "hdr":
            ax_lbl.text(0.0, yv, str(txt).upper(), ha="left", va="center", fontsize=7.5, fontweight="bold", color=theme["losat"])
        else:
            if row_i % 2 == 1:
                ax_lbl.axhspan(yv - 0.48, yv + 0.48, color=theme["band"], zorder=0)
            row_i += 1
            r = by_kt[(str(txt), 16)]
            ax_lbl.text(0.04, yv, r["mode_label"], ha="left", va="center", fontsize=8.8, fontweight="bold", color=theme["ink"])

    # Panels a (16t speed + dedicated % col) & b (1t speed + dedicated % col)
    for plot_rect, pct_rect, let, title, th in (
        ([0.168, bot_y, 0.195, h_ax], [0.366, bot_y, 0.056, h_ax], "a", "Speed vs native BLAST+, 16t", 16),
        ([0.448, bot_y, 0.185, h_ax], [0.636, bot_y, 0.056, h_ax], "b", "Speed vs native BLAST+, 1t",  1),
    ):
        ax = fig.add_axes(plot_rect)
        ax_p = fig.add_axes(pct_rect)
        ax.set_facecolor(theme["surface"])
        ax_p.set_facecolor(theme["surface"])
        ax_p.set_xlim(0, 1)
        ax_p.set_ylim(ymin, ymax)
        ax_p.axis("off")

        ax.set_xscale("log")
        x_hi = 25.0 if th == 16 else 7.2
        ax.set_xlim(0.50, x_hi)
        ax.set_ylim(ymin, ymax)
        ax.set_yticks([])
        ax.axvspan(1.0, x_hi, color=theme["zone"], alpha=0.55, zorder=0)

        row_i = 0
        for kind, yv, _txt in y_rows:
            if kind == "case":
                if row_i % 2 == 1:
                    ax.axhspan(yv - 0.48, yv + 0.48, color=theme["band"], alpha=0.6, zorder=0)
                    ax_p.axhspan(yv - 0.48, yv + 0.48, color=theme["band"], alpha=0.6, zorder=0)
                row_i += 1
        ax.axvline(1.0, color=theme["ref"], linewidth=1.2, zorder=2)
        ticks = [0.5, 1, 2, 5, 10, 20] if th == 16 else [0.5, 1, 2, 5]
        style_axis_log_x(ax, ticks, theme, labels=[f"{t}\u00d7" for t in ticks], tick_fs=7.6, skip_one=True)
        ax.set_xlabel("Native BLAST+ \u00f7 LOSAT time", fontsize=7.8, color=theme["ink2"])
        add_panel_label(ax, let, title, theme, y=1.024)
        ax_p.text(0.5, ymax - 0.15, "Wasm /\nnative", ha="center", va="center", fontsize=7.0, color=theme["muted"])

        for kind, yv, txt in y_rows:
            if kind != "case":
                continue
            r = by_kt[(str(txt), th)]
            pct = float(r["wasm_pct_of_native"])
            ax_p.text(0.5, yv, f"{pct:.0f}%", ha="center", va="center", fontsize=8.1, fontweight="semibold", color=theme["wasm"])
            if not r["native_speedup"]:
                ax.text(1.22, yv, "BLAST+ not measured", ha="left", va="center", fontsize=7.2, style="italic", color=theme["muted"])
                continue
            ns, ws = float(r["native_speedup"]), float(r["wasm_speedup"])
            ax.plot([ws, ns], [yv, yv], color=theme["hairline"], linewidth=1.8, zorder=3)
            ax.scatter([ws], [yv], marker="^", s=54, color=theme["wasm"], zorder=4)
            ax.scatter([ns], [yv], marker="o", s=50, color=theme["losat"], zorder=5)

    # Panel c: Peak memory at 16 threads
    ax_c = fig.add_axes([0.722, bot_y, 0.254, h_ax])
    ax_c.set_facecolor(theme["surface"])
    ax_c.set_xscale("log")
    ax_c.set_xlim(10.0, 12000.0)
    ax_c.set_ylim(ymin, ymax)
    ax_c.set_yticks([])
    row_i = 0
    for kind, yv, _txt in y_rows:
        if kind == "case":
            if row_i % 2 == 1:
                ax_c.axhspan(yv - 0.48, yv + 0.48, color=theme["band"], alpha=0.6, zorder=0)
            row_i += 1
    style_axis_log_x(ax_c, [20, 100, 500, 2000, 5000], theme, labels=["20", "100", "500", "2k", "5k"], tick_fs=7.6)
    ax_c.set_xlabel("\u2190 Smaller   Peak memory, 16t (MiB)", fontsize=7.8, color=theme["ink2"])
    add_panel_label(ax_c, "c", "Peak memory, 16 threads", theme, y=1.024)

    for kind, yv, txt in y_rows:
        if kind != "case":
            continue
        r = by_kt[(str(txt), 16)]
        nm = float(r["native_mib"])
        wm = float(r["wasm_mib"])
        bm = float(r["ncbi_mib"]) if r["ncbi_mib"] else None
        pts = [nm, wm] + ([bm] if bm else [])
        ax_c.plot([min(pts), max(pts)], [yv, yv], color=theme["hairline"], linewidth=1.4, zorder=2)
        if bm:
            ax_c.scatter([bm], [yv], marker="D", s=44, color=theme["ncbi"], zorder=3)
        ax_c.scatter([wm], [yv], marker="^", s=52, color=theme["wasm"], zorder=4)
        ax_c.scatter([nm], [yv], marker="o", s=48, color=theme["losat"], zorder=5)

    wasm_sha = short_sha(meta.get("wasm_commit") or "")
    derived_note = (
        f"Note: {wasm_sha} times are scaled from measured build {base_sha} by a side-by-side ratio run."
        if wasm_derived
        else "Direct paired native, WebAssembly, and NCBI BLAST+ runs in the same lease."
    )
    add_footer(
        fig,
        f"LOSAT {wasm_sha} (native & wasm32-wasip1-threads under Node.js) vs native NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')}, {short_hw(meta)}. "
        f"Wasm / native % = native time \u00f7 Wasm time (100% = native parity).\n{derived_note}",
        theme,
        y_rule=0.066,
    )
    return save_themed(fig, out, "fig5_webassembly", dark, theme)


# ---------------------------------------------------------------------------
# Figure 6: Generational Leap — LOSAT v0.1.0 -> v0.2.0 vs NCBI BLAST+ 2.17.0
# ---------------------------------------------------------------------------
def render_fig6(snap: Path, out: Path, theme: dict[str, str], dark: bool) -> Path:
    style_rc(theme)
    _std, by_kt = load_std_by_kt(snap)
    meta = read_meta(snap)

    cases = [
        ("bn",  "blastn",    "31k cDNA vs 3 Mb",       True),
        ("tbx", "tblastx",   "3 Candidatus genomes",   False),
        ("mb",  "megablast", "100 16S vs 28k RefSeq",  False),
        ("bp",  "blastp",    "3k human vs yeast",      True),
    ]

    fig = new_fig((12.6, 5.9), theme)
    add_header(
        fig,
        "LOSAT v0.1.0 \u2192 v0.2.0: up to 45.5\u00d7 faster and 50.3\u00d7 lower memory, while reaching exact BLAST+ parity",
        "All four search modes present in v0.1.0 (16 threads) alongside NCBI BLAST+ 2.17.0  \u00b7  v0.2.0 also expands from 3 to 17 programs and adds -db support",
        theme,
        y_title=0.964,
        y_sub=0.922,
    )

    handles = [
        Patch(facecolor=theme["v010"], label="LOSAT v0.1.0"),
        Patch(facecolor=theme["ncbi"], label="NCBI BLAST+ 2.17.0"),
        Patch(facecolor=theme["losat"], label="LOSAT v0.2.0 (byte-identical to NCBI)"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.018, 0.872), frameon=False, ncol=3, fontsize=8.1, columnspacing=1.5)

    for p_i, (rect, let, title, key_v0, key_nc, key_lo, xlabel, fmt_fn) in enumerate([
        ([0.148, 0.165, 0.365, 0.595], "a", "Wall time, 16 threads (s, log scale)",  "v010_wall_s",  "ncbi_wall_s",  "losat_wall_s",  "\u2190 Faster   Wall time (s)",            fmt_sec),
        ([0.592, 0.165, 0.382, 0.595], "b", "Peak memory, 16 threads (MiB, log scale)", "v010_rss_mib", "ncbi_rss_mib", "losat_rss_mib", "\u2190 Smaller   Peak memory (max RSS)", fmt_mem),
    ]):
        ax = fig.add_axes(rect)
        ax.set_facecolor(theme["surface"])
        ys = [3.0, 2.0, 1.0, 0.0]
        h = 0.24
        for yi, (cid, _label, _sub, _v0_same) in zip(ys, cases):
            r = by_kt[(cid, 16)]
            v0, nc, lo = float(r[key_v0]), float(r[key_nc]), float(r[key_lo])
            gain = v0 / lo

            if yi % 2 == 0:
                ax.axhspan(yi - 0.45, yi + 0.45, color=theme["band"], zorder=0)
            ax.barh([yi + h], [v0], height=h * 0.86, color=theme["v010"], zorder=3)
            ax.barh([yi],     [nc], height=h * 0.86, color=theme["ncbi"], zorder=3)
            ax.barh([yi - h], [lo], height=h * 0.86, color=theme["losat"], zorder=4)

            ax.text(v0 * 1.15, yi + h, fmt_fn(v0), ha="left", va="center", fontsize=7.3, color=theme["muted"])
            ax.text(nc * 1.15, yi,     fmt_fn(nc), ha="left", va="center", fontsize=7.3, color=theme["ink2"])
            ax.text(lo * 1.15, yi - h, f"{fmt_fn(lo)}  ({fmt_fold(gain, exact_tenth=True)} vs v0.1)", ha="left", va="center",
                    fontsize=7.6, fontweight="bold", color=theme["losat"])

        ax.set_xscale("log")
        if p_i == 0:
            ax.set_xlim(1.0, 32000.0)
            style_axis_log_x(ax, [1, 10, 100, 1000, 10000], theme, labels=["1", "10", "100", "1,000", "10,000"], tick_fs=7.7)
            ax.set_yticks(ys)
            ax.set_yticklabels([f"{lbl}{'*' if not same else ''}\n{sub}" for _, lbl, sub, same in cases], fontsize=8.0, color=theme["ink"])
            ax.tick_params(axis="y", length=0, pad=6)
        else:
            ax.set_xlim(20.0, 42000.0)
            style_axis_log_x(ax, [20, 100, 500, 2000, 5000], theme, labels=["20 MiB", "100", "500", "2k", "5k"], tick_fs=7.7)
            ax.set_yticks([])
        ax.set_xlabel(xlabel, fontsize=8.0, color=theme["ink2"])
        add_panel_label(ax, let, title, theme, y=1.035)

    add_footer(
        fig,
        f"LOSAT v0.1.0 release binary vs LOSAT {short_sha(meta.get('losat_commit', ''))} (v0.2.0) and NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')}, {short_hw(meta)}, 16 threads.\n"
        f"*v0.1.0 outputs for megablast and tblastx differed from NCBI BLAST+; v0.2.0 matches NCBI BLAST+ 2.17.0 byte-for-byte across all modes.",
        theme,
        y_rule=0.072,
    )
    return save_themed(fig, out, "fig6_v01_to_v02", dark, theme)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--snapshot", type=Path, default=Path(__file__).resolve().parent.parent / "snapshot")
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent.parent / "figures")
    ap.add_argument("--only", choices=["fig1", "fig2", "fig3", "fig4", "fig5", "fig6"], default=None)
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
        "fig5": render_fig5,
        "fig6": render_fig6,
    }
    targets = [args.only] if args.only else list(renderers.keys())
    for k in targets:
        for dark, theme in ((False, LIGHT), (True, DARK)):
            p = renderers[k](args.snapshot, args.out, theme, dark)
            if mirror:
                shutil.copyfile(p, mirror / p.name)
            print(p.name, sha256_file(p)[:16])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
