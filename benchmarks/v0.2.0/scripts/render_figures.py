#!/usr/bin/env python3
"""Deterministically render the 6-figure LOSAT v0.2.0 benchmark suite (light and dark PNGs).

Approved visual architecture:
  fig1_suite_overview : 12 search modes + 3 utilities (16t/1t speedup lollipops + shrinking memory % of BLAST+)
  fig2_real_pipelines : Real workloads >= 1 s/run (9 large workflows, 14 MT tools, 7 ST tools; speedup + shrinking RAM %)
  fig3_scaling        : 3x2 scaling grid (DB RAM, DB warm/cold time, query time+RAM, thread speed diagonal,
                        thread RAM curves for blastn/blastp/tblastn, and 1t->16t RAM shift across all 12 modes)
  fig4_vs_heuristics  : Memory-sized bubble Pareto on UniRef50 (vs DIAMOND & MMseqs2) + SILVA 144 horizontal dumbbells
  fig5_webassembly    : 3-engine physical-scale dumbbells (Native LOSAT, Wasm, Native BLAST+) + Wasm/native %
  fig6_v01_to_v02     : 3-stage generational slopegraph (LOSAT v0.1.0 -> NCBI BLAST+ 2.17.0 -> LOSAT v0.2.0)

Writes benchmarks/v0.2.0/figures/ and mirrors the exact PNG bytes to docs/manual/src/images/.
Deterministic bytes require Matplotlib 3.10.3, NumPy 2.2.5, and the Lato font (falls back to
Liberation Sans / DejaVu Sans when Lato is unavailable).
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
import textwrap
from collections import defaultdict
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
    ncbi_bar="#CBD5E1",
    ncbi_bar_alpha=0.42,
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
    ncbi_bar="#334155",
    ncbi_bar_alpha=0.72,
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
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
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


def fmt_pct(pct: float) -> str:
    if pct < 1.0:
        return f"{pct:.2f}%"
    if pct < 10.0:
        return f"{pct:.1f}%"
    return f"{pct:.0f}%"


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


def style_rc(theme: dict[str, object]) -> None:
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


def new_fig(figsize: tuple[float, float], theme: dict[str, object]) -> plt.Figure:
    style_rc(theme)
    fig = plt.figure(figsize=figsize)
    fig.patch.set_facecolor(theme["surface"])
    return fig


def save_themed(fig: plt.Figure, out: Path, stem: str, dark: bool, theme: dict[str, object], dpi: int = 200) -> Path:
    p = out / f"{stem}{'_dark.png' if dark else '.png'}"
    fig.savefig(
        p,
        dpi=dpi,
        facecolor=theme["surface"],
        edgecolor="none",
        metadata={"Software": "LOSAT benchmark renderer v3"},
    )
    plt.close(fig)
    return p


def load_std_by_kt(snap: Path) -> tuple[list[dict[str, str]], dict[tuple[str, int], dict[str, str]]]:
    std = read_tsv(snap / "standard_workloads.tsv")
    return std, {(r["case"], int(r["threads"])): r for r in std}


def add_header(fig: plt.Figure, title: str, subtitle: str, theme: dict[str, object], y_title: float = 0.970, y_sub: float = 0.938) -> None:
    fig.text(0.020, y_title, title, ha="left", va="top", fontsize=12.8, fontweight="bold", color=theme["ink"])
    fig.text(0.020, y_sub, subtitle, ha="left", va="top", fontsize=8.4, color=theme["muted"])
    fig.add_artist(Line2D([0.020, 0.980], [y_sub - 0.028, y_sub - 0.028], transform=fig.transFigure, color=theme["grid"], linewidth=0.9))


def add_footer(fig: plt.Figure, text: str, theme: dict[str, object], y_rule: float = 0.056) -> None:
    fig.add_artist(Line2D([0.020, 0.980], [y_rule, y_rule], transform=fig.transFigure, color=theme["grid"], linewidth=0.7))
    fig.text(0.020, y_rule - 0.009, text, ha="left", va="top", fontsize=7.0, color=theme["muted"], linespacing=1.34)


def add_panel(ax: plt.Axes, letter: str, title: str, theme: dict[str, object], y: float = 1.024, sub: str = "") -> None:
    label = f"{letter}  {title}" if letter else title
    ax.text(0.0, y, label, transform=ax.transAxes, fontsize=9.4, fontweight="bold", va="bottom", ha="left", color=theme["ink"])
    if sub:
        ax.text(1.0, y + 0.003, sub, transform=ax.transAxes, fontsize=7.3, va="bottom", ha="right", color=theme["muted"])


def hide_minor_x(ax: plt.Axes) -> None:
    ax.xaxis.set_minor_formatter(NullFormatter())


def hide_minor_y(ax: plt.Axes) -> None:
    ax.yaxis.set_minor_formatter(NullFormatter())


def style_pct_mem_x(ax_m: plt.Axes, ticks: list[int], labels: list[str], theme: dict[str, object], tick_fs: float = 7.4, ref_lw: float = 1.1) -> None:
    for tk in ticks:
        ax_m.axvline(tk, color=theme["ref"] if tk == 100 else theme["grid"], linewidth=ref_lw if tk == 100 else 0.6, zorder=1)
    ax_m.set_xticks(ticks)
    ax_m.set_xticklabels(labels, fontsize=tick_fs, color=theme["ink2"])
    ax_m.tick_params(axis="x", which="major", color=theme["hairline"], labelcolor=theme["ink2"], length=3.0, width=0.8)
    for sp in ("top", "right", "left"):
        ax_m.spines[sp].set_visible(False)
    ax_m.spines["bottom"].set_color(theme["hairline"])
    ax_m.spines["bottom"].set_bounds(0, 100)


def style_log_x(
    ax: plt.Axes,
    ticks: list[float],
    labels: list[str],
    theme: dict[str, object],
    tick_fs: float = 7.4,
    skip_one: bool = False,
    bound_lo: float | None = None,
    bound_hi: float | None = None,
) -> None:
    ax.set_xscale("log")
    for tk in ticks:
        if not (skip_one and tk == 1):
            ax.axvline(tk, color=theme["grid"], linewidth=0.6, zorder=1)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=tick_fs, color=theme["ink2"])
    hide_minor_x(ax)
    ax.tick_params(axis="x", which="minor", length=0)
    ax.tick_params(axis="x", which="major", color=theme["hairline"], labelcolor=theme["ink2"], length=3.0, width=0.8)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(theme["hairline"])
    if bound_lo is not None and bound_hi is not None:
        ax.spines["bottom"].set_bounds(bound_lo, bound_hi)


def style_xy_grid(ax: plt.Axes, theme: dict[str, object], lw: float = 0.6, plain_y: bool = False, clear_mx: bool = True) -> None:
    if clear_mx:
        hide_minor_x(ax)
    hide_minor_y(ax)
    if plain_y:
        ax.yaxis.set_major_formatter(FuncFormatter(plain_log))
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(True, which="major", color=theme["grid"], linewidth=lw, zorder=0)


def draw_broken_bar(ax: plt.Axes, yv: float, val: float, cap: float, color: str, theme: dict[str, object], height: float = 0.52) -> float:
    bar_end = min(val, cap * 0.96) if val > cap else val
    if bar_end >= 1.0:
        ax.barh([yv], [bar_end - 1.0], left=1.0, height=height, color=color, edgecolor="none", zorder=3)
    else:
        ax.barh([yv], [1.0 - bar_end], left=bar_end, height=height, color=theme["muted"], alpha=0.55, edgecolor="none", zorder=3)
    if val > cap:
        log_pos = math.exp(0.78 * math.log(cap))
        for shift in (-0.045, 0.045):
            x0 = log_pos * math.exp(shift - 0.028)
            x1 = log_pos * math.exp(shift + 0.028)
            ax.plot([x0, x1], [yv - height * 0.62, yv + height * 0.62], color=theme["surface"], linewidth=3.2, zorder=5)
            ax.plot([x0, x1], [yv - height * 0.62, yv + height * 0.62], color=theme["ink2"], linewidth=0.95, zorder=6)
    return bar_end


def draw_shrink_mem_bar(
    ax: plt.Axes,
    yv: float,
    losat_mib: float,
    ncbi_mib: float,
    theme: dict[str, object],
    height: float = 0.52,
) -> float:
    """Draw a 100% BLAST+ ghost bar with LOSAT's shrinking remaining footprint inside (shorter = smaller RAM)."""
    pct = 100.0 * losat_mib / ncbi_mib
    ax.barh([yv], [100.0], left=0.0, height=height, color=theme["ncbi_bar"], alpha=float(theme["ncbi_bar_alpha"]), edgecolor="none", zorder=2)
    vis = min(max(pct, 1.15), 100.0)
    ax.barh([yv], [vis], left=0.0, height=height, color=theme["wasm"], edgecolor="none", zorder=4)
    return pct


def add_legend(
    fig: plt.Figure,
    handles: list[object],
    y: float = 0.908,
    ncol: int = 3,
    fs: float = 7.5,
    colsp: float = 1.1,
    hlen: float | None = None,
) -> None:
    kw: dict[str, object] = dict(handles=handles, loc="upper left", bbox_to_anchor=(0.018, y), frameon=False, ncol=ncol, fontsize=fs, columnspacing=colsp)
    if hlen is not None:
        kw["handlelength"] = hlen
    fig.legend(**kw)


def mem_shrink_patch(theme: dict[str, object]) -> Patch:
    return Patch(facecolor=theme["wasm"], label="LOSAT memory (% of BLAST+, shorter = smaller)")


def draw_row_badge(ax_lbl: plt.Axes, yv: float, title: str, desc: str, badge: str, theme: dict[str, object]) -> None:
    ax_lbl.text(0.02, yv + 0.18, title, ha="left", va="center", fontsize=8.7, fontweight="bold", color=theme["ink"])
    ax_lbl.text(0.02, yv - 0.22, desc, ha="left", va="center", fontsize=6.9, color=theme["muted"])
    ax_lbl.text(
        0.96, yv + 0.02, badge, ha="right", va="center", fontsize=6.6, color=theme["badge_ink"],
        bbox=dict(boxstyle="round,pad=0.20", facecolor=theme["badge_bg"], edgecolor=theme["grid"], linewidth=0.6),
    )


# ---------------------------------------------------------------------------
# Figure 1 (Pattern 1B): Speedup Lollipops (16t & 1t) + Shrinking Memory Footprint (%)
# ---------------------------------------------------------------------------
def render_fig1(snap: Path, out: Path, theme: dict[str, object], dark: bool) -> Path:
    std, by_kt = load_std_by_kt(snap)
    util = read_tsv(snap / "utilities.tsv")
    meta = read_meta(snap)

    y_rows: list[tuple[str, float, object]] = []
    y = 0.0
    for gi, (gtitle, cids) in enumerate(GROUPS_FIG1):
        if gi > 0:
            y -= 0.58
        y_rows.append(("hdr", y, gtitle))
        y -= 0.90
        for cid in cids:
            y_rows.append(("std", y, cid))
            y -= 1.12
    y -= 0.72
    y_rows.append(("sep", y + 0.36, None))
    y_rows.append(("hdr", y, "Database & Domain Utilities (1 thread)"))
    y -= 0.90
    for ur in util:
        y_rows.append(("util", y, ur))
        y -= 1.12
    ymin, ymax = y + 0.30, 0.88

    std16 = [r for r in std if int(r["threads"]) == 16]
    non_tbx = sorted((float(r["speedup"]), r["label"]) for r in std16 if r["case"] != "tbx")
    tbx_sp = float(by_kt[("tbx", 16)]["speedup"])
    pcts16 = [100.0 * float(r["losat_rss_mib"]) / float(r["ncbi_rss_mib"]) for r in std16]
    reds16 = [float(r["mem_reduction"]) for r in std16]
    fig = new_fig((13.4, 8.8), theme)
    add_header(
        fig,
        f"LOSAT v{meta.get('losat_version', '0.2.0')} outperforms NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')} across all {len(std16)} search modes and {len(util)} utilities",
        f"16-thread wall time: {fmt_fold(non_tbx[0][0])} ({non_tbx[0][1]}) to {fmt_fold(non_tbx[-1][0], True)} faster ({fmt_fold(tbx_sp)} on tblastx)  \u00b7  "
        f"Peak memory: {fmt_pct(min(pcts16))} to {fmt_pct(max(pcts16))} of BLAST+ ({fmt_fold(min(reds16))}\u2013{fmt_fold(max(reds16))} lower)  \u00b7  Same output form for form (searches)",
        theme,
        y_title=0.974,
        y_sub=0.944,
    )

    handles = [
        Line2D([0], [0], marker="o", ls="none", ms=6.8, color=theme["losat"], label="16-thread speedup (BLAST+ \u00f7 LOSAT)"),
        Line2D([0], [0], marker="o", ls="none", ms=6.2, markerfacecolor=theme["surface"], markeredgecolor=theme["losat"], mew=1.6, label="1-thread speedup"),
        Line2D([0], [0], marker="o", ls="none", ms=6.4, color=theme["wasm"], label="1-thread utility speedup"),
        Patch(facecolor=theme["ncbi_bar"], alpha=float(theme["ncbi_bar_alpha"]), label="BLAST+ memory = 100%"),
        mem_shrink_patch(theme),
    ]
    add_legend(fig, handles)

    bot_y, top_h = 0.105, 0.700
    ax_lbl = fig.add_axes([0.020, bot_y, 0.220, top_h])
    ax_lbl.set_facecolor(theme["surface"])
    ax_lbl.set_xlim(0, 1)
    ax_lbl.set_ylim(ymin, ymax)
    ax_lbl.axis("off")

    ax_s = fig.add_axes([0.252, bot_y, 0.345, top_h])
    ax_m = fig.add_axes([0.630, bot_y, 0.348, top_h])
    for ax in (ax_s, ax_m):
        ax.set_facecolor(theme["surface"])
        ax.set_ylim(ymin, ymax)
        ax.set_yticks([])

    ax_s.set_xlim(0.85, 4200.0)
    ax_s.axvspan(1.0, 4200.0, color=theme["zone"], alpha=0.55, zorder=0)
    ax_m.set_xlim(0.0, 152.0)

    ri = 0
    for kind, yv, item in y_rows:
        if kind in ("std", "util"):
            if ri % 2 == 1:
                ax_lbl.axhspan(yv - 0.52, yv + 0.52, color=theme["band"], zorder=0)
                ax_s.axhspan(yv - 0.52, yv + 0.52, color=theme["band"], alpha=0.65, zorder=0)
                ax_m.axhspan(yv - 0.52, yv + 0.52, color=theme["band"], alpha=0.65, zorder=0)
            ri += 1
        elif kind == "sep":
            for ax in (ax_lbl, ax_s, ax_m):
                ax.axhline(yv, color=theme["hairline"], linewidth=0.8, linestyle="--", zorder=1)
        elif kind == "hdr":
            ax_lbl.text(0.0, yv, str(item).upper(), ha="left", va="center", fontsize=7.5, fontweight="bold", color=theme["losat"])
        if kind == "std":
            r = by_kt[(str(item), 16)]
            draw_row_badge(ax_lbl, yv, r["label"], r["short_desc"], "-db" if r["target"] == "-db" else "-subj", theme)
        elif kind == "util":
            r = item
            draw_row_badge(ax_lbl, yv, r["label"], r["short_desc"], "1t", theme)

    ax_s.axvline(1.0, color=theme["ref"], linewidth=1.15, zorder=2)
    style_log_x(ax_s, [1, 2, 5, 10, 50, 200, 600], ["1\u00d7", "2\u00d7", "5\u00d7", "10\u00d7", "50\u00d7", "200\u00d7", "600\u00d7"], theme, tick_fs=7.5, skip_one=True)
    ax_s.set_xlabel("Faster \u2192   Wall-time speedup (BLAST+ \u00f7 LOSAT)", fontsize=7.8, color=theme["ink2"], labelpad=4)
    add_panel(ax_s, "a", "Wall-time speedup (16t & 1t combined)", theme, y=1.018)

    style_pct_mem_x(ax_m, [0, 25, 50, 75, 100], ["0%", "25%", "50%", "75%", "100%\n(BLAST+)"], theme, tick_fs=7.4, ref_lw=1.1)
    ax_m.set_xlabel("\u2190 Smaller   LOSAT peak memory (% of BLAST+ max RSS)", fontsize=7.8, color=theme["ink2"], labelpad=4)
    add_panel(ax_m, "b", "Memory footprint, 16t (shorter = smaller)", theme, y=1.018)

    def _put_mem_row(yv: float, lm: float, nm: float, red: float, bold: bool, col: object) -> None:
        pct = draw_shrink_mem_bar(ax_m, yv, lm, nm, theme)
        ax_m.text(103.5, yv + 0.16, f"{fmt_pct(pct)}  ({fmt_fold(red)} lower)",
                  ha="left", va="center", fontsize=7.5, fontweight="bold" if bold else "semibold", color=col)
        ax_m.text(103.5, yv - 0.21, f"{fmt_mem(lm)} vs {fmt_mem(nm)}", ha="left", va="center", fontsize=6.5, color=theme["muted"])

    for kind, yv, item in y_rows:
        if kind == "std":
            r16, r1 = by_kt[(str(item), 16)], by_kt[(str(item), 1)]
            s16, s1 = float(r16["speedup"]), float(r1["speedup"])
            ax_s.plot([1.0, max(s16, s1)], [yv, yv], color=theme["losat_soft"], linewidth=2.2, zorder=3)
            ax_s.scatter([s1], [yv], marker="o", s=38, facecolors=theme["surface"], edgecolors=theme["losat"], linewidths=1.6, zorder=4)
            ax_s.scatter([s16], [yv], marker="o", s=48, color=theme["losat"], zorder=5)
            ax_s.text(max(s16, s1) * 1.20, yv + 0.16, f"{fmt_fold(s16)} (16t)  \u00b7  {fmt_fold(s1)} (1t)",
                      ha="left", va="center", fontsize=7.5, fontweight="bold" if s16 >= 9.5 else "semibold", color=theme["ink"])
            ax_s.text(max(s16, s1) * 1.20, yv - 0.21, f"16t: {fmt_sec(float(r16['losat_wall_s']))} vs {fmt_sec(float(r16['ncbi_wall_s']))}",
                      ha="left", va="center", fontsize=6.5, color=theme["muted"])
            lm, nm = float(r16["losat_rss_mib"]), float(r16["ncbi_rss_mib"])
            pct_v = 100.0 * lm / nm
            _put_mem_row(yv, lm, nm, float(r16["mem_reduction"]), pct_v <= 12.0, theme["wasm"] if pct_v <= 18.0 else theme["ink"])
        elif kind == "util":
            r = item
            s1 = float(r["speedup"])
            ax_s.plot([1.0, s1], [yv, yv], color=theme["wasm_soft"], linewidth=2.2, zorder=3)
            ax_s.scatter([s1], [yv], marker="o", s=46, color=theme["wasm"], zorder=5)
            ax_s.text(s1 * 1.20, yv + 0.16, f"{fmt_fold(s1)} (1t utility)", ha="left", va="center", fontsize=7.5, fontweight="semibold", color=theme["wasm"])
            ax_s.text(s1 * 1.20, yv - 0.21, f"{fmt_sec(float(r['losat_wall_s']))} vs {fmt_sec(float(r['ncbi_wall_s']))}", ha="left", va="center", fontsize=6.5, color=theme["muted"])
            _put_mem_row(yv, float(r["losat_rss_mib"]), float(r["ncbi_rss_mib"]), float(r["mem_reduction"]), False, theme["wasm"])

    def _runs_summary(rows_t: list[dict[str, str]]) -> str:
        by_pair: dict[tuple[int, int], list[str]] = defaultdict(list)
        for r in rows_t:
            by_pair[(int(r["losat_runs"]), int(r["ncbi_runs"]))].append(r["label"])
        main_pair = max(by_pair, key=lambda k: len(by_pair[k]))
        lo_m, nc_m = main_pair
        head = f"{lo_m} rounds" if lo_m == nc_m else f"{lo_m} LOSAT / {nc_m} NCBI"
        excs = []
        for pair, lbls in sorted(by_pair.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            if pair == main_pair:
                continue
            lo_n, nc_n = pair
            if lo_n == lo_m and nc_n != nc_m:
                excs.append(f"NCBI {', '.join(lbls)}: {nc_n}")
            elif lo_n == nc_n:
                excs.append(f"{', '.join(lbls)}: {lo_n}")
            else:
                excs.append(f"{', '.join(lbls)}: {lo_n} LOSAT / {nc_n} NCBI")
        return f"{head} ({'; '.join(excs)})" if excs else head

    std1 = [r for r in std if int(r["threads"]) == 1]
    util_runs_set = sorted({int(r["runs"]) for r in util})
    util_runs_txt = "/".join(str(x) for x in util_runs_set)
    add_footer(
        fig,
        f"LOSAT {short_sha(meta.get('losat_commit', ''))} vs NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')} & rpsbproc {meta.get('rpsbproc_version', '0.5.1')}, {short_hw(meta)}. "
        f"Medians at 16t: {_runs_summary(std16)}; at 1t: {_runs_summary(std1)}; utilities: {util_runs_txt}.\n"
        f"[-db] paired BLAST databases; [-subj] LOSAT -subject vs BLAST+ -db (blastp/blastp-fast parity confirmed form-for-form). BLAST+ max RSS includes memory-mapped DB pages.",
        theme,
        y_rule=0.052,
    )
    return save_themed(fig, out, "fig1_suite_overview", dark, theme)


# ---------------------------------------------------------------------------
# Figure 2 (Pattern 2B): Speedup Rightward + Shrinking Memory Footprint (% of BLAST+)
# ---------------------------------------------------------------------------
def render_fig2(snap: Path, out: Path, theme: dict[str, object], dark: bool) -> Path:
    large = read_tsv(snap / "large_workflows.tsv")
    dt = read_tsv(snap / "downstream_tools.tsv")
    meta = read_meta(snap)

    real_mt, real_st = [], []
    startup_by_tool: dict[str, list[str]] = defaultdict(list)
    for r in dt:
        if r["same_output"] != "yes":
            sys.exit(f"downstream_tools parity failure for {r['tool']} ({r['program']}): {r['same_output']}")
        if r["mt_speedup"]:
            if startup_dominated(r["ncbi_mt_s"], r["cases_count"]):
                is_startup = True
            else:
                real_mt.append(r)
                is_startup = False
        elif r["note"] == "losat_under_10ms":
            is_startup = True
        elif r["n1_speedup"]:
            if startup_dominated(r["ncbi_n1_s"], r["cases_count"]):
                is_startup = True
            else:
                real_st.append(r)
                is_startup = False
        else:
            is_startup = False
        if is_startup:
            startup_by_tool[r["tool"]].append(r["program"])

    n_tools = len({r["tool"] for r in dt})
    real_mt.sort(key=lambda r: float(r["mt_speedup"]), reverse=True)
    real_st.sort(key=lambda r: float(r["n1_speedup"]), reverse=True)
    # Keep 1M reads (stdin) and 1M reads (file) directly adjacent at the top of Panel a so the stdin contrast is immediate
    large_s = sorted(
        large,
        key=lambda r: (0 if r["case"] in ("M2_reads1e6s", "M2_reads1e6") else 1, -float(r["speedup"])),
    )

    fig = new_fig((14.2, 9.5), theme)
    add_header(
        fig,
        "Drop-in speed and memory gains across real bioinformatics workflows (excluding < 1 s start-up runs)",
        f"Unmodified BLAST+ commands via losat --install-shims  \u00b7  Only workloads where BLAST+ took \u2265 1 s/run are plotted  \u00b7  All {n_tools} tested tools matched NCBI output",
        theme,
        y_title=0.975,
        y_sub=0.947,
    )
    handles = [
        Patch(facecolor=theme["losat"], label="Wall-time speedup (BLAST+ \u00f7 LOSAT, longer = faster)"),
        Patch(facecolor=theme["ncbi_bar"], alpha=float(theme["ncbi_bar_alpha"]), label="NCBI BLAST+ memory = 100%"),
        mem_shrink_patch(theme),
        Line2D([0], [0], ls="none", marker="", label="\u2020 NCBI multi-thread output drifts; LOSAT matches 1-thread NCBI"),
    ]
    add_legend(fig, handles, y=0.911, ncol=4, fs=7.4, colsp=1.0, hlen=1.3)

    def _calls_label(cases_count: str) -> str:
        n = int(cases_count)
        return "1 BLAST call" if n == 1 else f"{n} calls"

    def _item(name: str, sub: str, sp: float, ls: float, ns: float, lm: float, nm: float, drift: bool) -> dict[str, object]:
        return dict(name=name, sub=sub, sp=sp, time_txt=f"{fmt_sec(ls)} vs {fmt_sec(ns)}", lm=lm, nm=nm, drift=drift)

    items_a = [_item(r["label"], f"{r['program']} \u00b7 {r['threads']}t \u00b7 {r['detail']}", float(r["speedup"]),
                     float(r["losat_wall_s"]), float(r["ncbi_wall_s"]), float(r["losat_rss_mib"]), float(r["ncbi_rss_mib"]),
                     r["ncbi_mt_differs"] == "yes") for r in large_s]
    items_b = [_item(f"{r['tool']} ({r['program']})" if r["tool"] in ("Kaptive", "Trinotate") else r["tool"],
                     f"{r['program']} \u00b7 {r['threads_label']}t \u00b7 {_calls_label(r['cases_count'])}", float(r["mt_speedup"]),
                     float(r["losat_mt_s"]), float(r["ncbi_mt_s"]), float(r["losat_mt_rss_mib"]), float(r["ncbi_mt_rss_mib"]),
                     r["ncbi_mt_differs"] == "yes") for r in real_mt]
    items_c = [_item(r["tool"], f"{r['program']} \u00b7 1t \u00b7 {_calls_label(r['cases_count'])}", float(r["n1_speedup"]),
                     float(r["losat_n1_s"]), float(r["ncbi_n1_s"]), float(r["losat_n1_rss_mib"]), float(r["ncbi_n1_rss_mib"]),
                     False) for r in real_st]

    def draw_block(items, rect_l, rect_s, rect_m, let, title, scap, sticks):
        n = len(items)
        ys = [float(n - 1 - i) for i in range(n)]
        ymin_p, ymax_p = -0.65, n - 0.35
        ax_l, ax_s, ax_m = fig.add_axes(rect_l), fig.add_axes(rect_s), fig.add_axes(rect_m)
        for a in (ax_l, ax_s, ax_m):
            a.set_facecolor(theme["surface"])
            a.set_ylim(ymin_p, ymax_p)
        ax_l.set_xlim(0, 1)
        ax_l.axis("off")
        add_panel(ax_l, let, title, theme, y=1.032)

        ax_s.set_xlim(0.88, scap * 5.8)
        ax_s.set_yticks([])
        ax_s.axvspan(1.0, scap, color=theme["zone"], alpha=0.55, zorder=0)
        ax_s.axvline(1.0, color=theme["ref"], linewidth=1.1, zorder=2)

        ax_m.set_xlim(0.0, 205.0)
        ax_m.set_yticks([])

        for i, yv in enumerate(ys):
            if i % 2 == 1:
                ax_l.axhspan(yv - 0.46, yv + 0.46, color=theme["band"], zorder=0)
                ax_s.axhspan(yv - 0.46, yv + 0.46, color=theme["band"], alpha=0.6, zorder=0)
                ax_m.axhspan(yv - 0.46, yv + 0.46, color=theme["band"], alpha=0.6, zorder=0)

        style_log_x(ax_s, sticks, [f"{t}\u00d7" for t in sticks], theme, tick_fs=7.2, skip_one=True, bound_lo=0.88, bound_hi=scap)
        ax_s.set_xlabel("Faster \u2192 (BLAST+ \u00f7 LOSAT)", fontsize=7.2, color=theme["ink2"], labelpad=3)
        ax_s.text(0.0, 1.032, "Speedup", transform=ax_s.transAxes, fontsize=8.0, fontweight="bold", va="bottom", color=theme["losat"])

        style_pct_mem_x(ax_m, [0, 50, 100], ["0%", "50%", "100%"], theme, tick_fs=7.2, ref_lw=1.0)
        ax_m.set_xlabel("\u2190 Smaller (% of BLAST+ RAM)", fontsize=7.2, color=theme["ink2"], labelpad=3)
        ax_m.text(0.0, 1.032, "Memory used", transform=ax_m.transAxes, fontsize=8.0, fontweight="bold", va="bottom", color=theme["wasm"])

        for yv, it in zip(ys, items):
            dag = " \u2020" if it["drift"] else ""
            ax_l.text(0.01, yv + 0.15, f"{it['name']}{dag}", ha="left", va="center", fontsize=8.1, fontweight="bold", color=theme["ink"])
            ax_l.text(0.01, yv - 0.21, str(it["sub"]), ha="left", va="center", fontsize=6.4, color=theme["muted"])

            end_x = draw_broken_bar(ax_s, yv, float(it["sp"]), scap, theme["losat"], theme, height=0.52)
            ax_s.text(end_x * 1.15, yv + 0.15, fmt_fold(float(it["sp"])), ha="left", va="center", fontsize=7.8,
                      fontweight="bold" if float(it["sp"]) >= 9.5 else "semibold", color=theme["ink"], zorder=7)
            ax_s.text(end_x * 1.15, yv - 0.21, str(it["time_txt"]), ha="left", va="center", fontsize=6.3, color=theme["muted"], zorder=7)

            lm, nm = float(it["lm"]), float(it["nm"])
            pct = draw_shrink_mem_bar(ax_m, yv, lm, nm, theme, height=0.52)
            ax_m.text(104.5, yv + 0.15, f"{fmt_pct(pct)}  ({fmt_fold(nm / lm)})", ha="left", va="center", fontsize=7.5,
                      fontweight="bold" if pct < 20.0 else "semibold", color=theme["wasm"] if pct < 25.0 else theme["ink"])
            ax_m.text(104.5, yv - 0.21, f"{fmt_mem(lm, True)} / {fmt_mem(nm, True)}", ha="left", va="center", fontsize=6.3, color=theme["muted"])

    draw_block(items_a, [0.020, 0.474, 0.205, 0.362], [0.230, 0.474, 0.118, 0.362], [0.356, 0.474, 0.132, 0.362],
               "a", "Large reference workflows (16 threads)", 12.0, [1, 2, 5, 10])
    draw_block(items_c, [0.020, 0.120, 0.205, 0.270], [0.230, 0.120, 0.118, 0.270], [0.356, 0.120, 0.132, 0.270],
               "c", "Single-thread published tools (\u2265 1 s/run)", 8.0, [1, 2, 4, 8])
    fig.add_artist(Line2D([0.502, 0.502], [0.110, 0.860], transform=fig.transFigure, color=theme["grid"], linewidth=0.9))
    draw_block(items_b, [0.518, 0.120, 0.175, 0.716], [0.698, 0.120, 0.128, 0.716], [0.834, 0.120, 0.146, 0.716],
               "b", "Multi-thread published tools (\u2265 1 s/run)", 22.0, [1, 2, 5, 10, 20])

    dt_commits = meta.get("dt_commits", {})
    one_shas = short_sha(",".join(dt_commits[k] for k in sorted(dt_commits) if dt_commits[k]))
    excl_parts = [
        f"{t} ({progs[0]})" if len(progs) == 1 else f"{t} ({len(progs)} modes)"
        for t, progs in sorted(startup_by_tool.items())
    ]
    by_lcase = {r["case"]: r for r in large}
    stdin_note = ""
    if "M2_reads1e6s" in by_lcase and "M2_reads1e6" in by_lcase:
        s_nc, f_nc = float(by_lcase["M2_reads1e6s"]["ncbi_wall_s"]), float(by_lcase["M2_reads1e6"]["ncbi_wall_s"])
        s_lo, f_lo = float(by_lcase["M2_reads1e6s"]["losat_wall_s"]), float(by_lcase["M2_reads1e6"]["losat_wall_s"])
        stdin_note = (
            f" Piping 1M reads via stdin slows BLAST+ {s_nc / f_nc:.1f}\u00d7 ({fmt_sec(f_nc)} \u2192 {fmt_sec(s_nc)}) "
            f"vs {s_lo / f_lo:.2f}\u00d7 for LOSAT ({fmt_sec(f_lo)} \u2192 {fmt_sec(s_lo)})."
        )
    excl_wrapped = textwrap.fill(
        f"Plotted times are summed over the BLAST calls in each workflow (medians of 3 rounds in a, b; 4 runs in c); memory shows peak max RSS (% of BLAST+, shorter = smaller).{stdin_note} "
        f"Excluded < 1 s/run rows ({', '.join(excl_parts)}) also matched NCBI output (see snapshot/downstream_tools.tsv).",
        width=212,
    )
    add_footer(
        fig,
        f"Panels a, b: LOSAT {short_sha(meta.get('toolmt_commit', ''))}; Panel c: LOSAT {one_shas}; "
        f"vs NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')}, {short_hw(meta)}.\n{excl_wrapped}",
        theme,
        y_rule=0.072,
    )
    return save_themed(fig, out, "fig2_real_pipelines", dark, theme)


# ---------------------------------------------------------------------------
# Figure 3 (Pattern 3A Extended — 3x2): DB Size, Query Volume, Thread Speed & RAM, All-12-Mode 1t->16t RAM
# ---------------------------------------------------------------------------
def render_fig3(snap: Path, out: Path, theme: dict[str, object], dark: bool) -> Path:
    rows = read_tsv(snap / "scaling.tsv")
    _std, by_kt = load_std_by_kt(snap)
    meta = read_meta(snap)

    db_rows = [r for r in rows if r["panel"] == "db_size"]
    warm = sorted([r for r in db_rows if r["cache"] == "warm"], key=lambda r: float(r["x_value"]))
    cold = sorted([r for r in db_rows if r["cache"] == "cold"], key=lambda r: float(r["x_value"]))
    q_rows = sorted([r for r in rows if r["panel"] == "query_size"], key=lambda r: float(r["x_value"]))
    t_rows = [r for r in rows if r["panel"] == "threads"]
    all12 = [cid for _, cids in GROUPS_FIG1 for cid in cids]

    core_lo_mib = float(warm[-1]["losat_rss_mib"])
    core_gb = float(warm[-1]["x_value"]) / 1000.0
    abr_mb = float(warm[0]["x_value"])
    abr_nc_s = float(warm[0]["ncbi_wall_s"])
    q_lo_k = int(round(float(q_rows[0]["x_value"]) / 1000.0))
    q_hi_m = int(round(float(q_rows[-1]["x_value"]) / 1e6))
    t_max = max(int(r["x_value"]) for r in t_rows)
    t_sp_lo = min(float(r["speedup"]) for r in t_rows)
    t_sp_hi = max(float(r["speedup"]) for r in t_rows)

    fig = new_fig((13.4, 10.6), theme)
    add_header(
        fig,
        f"Resident memory peaks at {fmt_mem(core_lo_mib)} on {core_gb:.0f} GB core_nt while scaling cleanly across queries, threads, and modes",
        f"DB size ({abr_mb:.2f} MB \u2192 {core_gb:.1f} GB), query volume ({q_lo_k}k \u2192 {q_hi_m}M reads), thread speed & memory (1 \u2192 {t_max} cores), and 1t \u2192 16t memory across all {len(all12)} modes",
        theme,
        y_title=0.978,
        y_sub=0.952,
    )
    handles = [
        Line2D([0], [0], marker="o", color=theme["losat"], ms=6.0, lw=2.0, label=f"LOSAT {meta.get('losat_version', '0.2.0')} (warm / 16t)"),
        Line2D([0], [0], marker="D", color=theme["ncbi"], ms=5.2, lw=1.8, label=f"NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')} (warm / 16t)"),
        Line2D([0], [0], marker="o", markerfacecolor=theme["surface"], markeredgecolor=theme["losat"], ms=5.8, ls="none", mew=1.5, label="Cold cache (b) / 1-thread (f)"),
        Line2D([0], [0], ls=":", color=theme["muted"], lw=1.2, label="Ideal linear scaling (d)"),
    ]
    add_legend(fig, handles, y=0.927, ncol=4, fs=7.8, colsp=1.2)

    db_names = {"abr_genes": "abr_genes", "16s_all": "16S rRNA", "ref_prok_rep_genomes": "ref_prok", "core_nt": "core_nt"}

    def _db_lbl(rec: dict[str, str]) -> str:
        mb = float(rec["x_value"])
        sz = f"{mb / 1000.0:.0f} GB" if mb >= 1000 else (f"{mb:.0f} MB" if mb >= 1 else f"{mb:.2f} MB")
        return f"{db_names.get(rec['x_label'], rec['x_label'])}\n{sz}"

    db_short = {r["x_label"]: _db_lbl(r) for r in warm}
    xs_w = [float(r["x_value"]) for r in warm]
    cold_arrow = dict(arrowstyle="-|>", color=theme["muted"], lw=0.8, shrinkA=0, shrinkB=4)
    pill_box = dict(boxstyle="round,pad=0.13", facecolor=theme["surface"], edgecolor="none", alpha=0.88)

    # --- Row 1, Panel a: DB size -> Peak memory ---
    ax_a = fig.add_axes([0.072, 0.656, 0.402, 0.210])
    ax_a.set_facecolor(theme["surface"])
    yn_m, yl_m = [float(r["ncbi_rss_mib"]) for r in warm], [float(r["losat_rss_mib"]) for r in warm]
    ax_a.fill_between(xs_w, yl_m, yn_m, color=theme["losat_fill"], alpha=0.55, zorder=1)
    ax_a.plot(xs_w, yn_m, marker="D", ms=5.4, color=theme["ncbi"], lw=1.9, zorder=3)
    ax_a.plot(xs_w, yl_m, marker="o", ms=6.0, color=theme["losat"], lw=2.1, zorder=4)
    for i, (x, yn, yl) in enumerate(zip(xs_w, yn_m, yl_m)):
        if i in (2, 3):
            ax_a.text(x, math.sqrt(yn * yl), fmt_fold(yn / yl), ha="center", va="center", fontsize=8.2, fontweight="bold", color=theme["losat"],
                      bbox=dict(boxstyle="round,pad=0.16", facecolor=theme["surface"], edgecolor=theme["zone_edge"], lw=0.7), zorder=6)
        else:
            ax_a.text(x, yl * 0.50, fmt_fold(yn / yl), ha="center", va="top", fontsize=7.4, color=theme["ink2"])
    ax_a.text(xs_w[-1], yn_m[-1] * 1.60, fmt_mem(yn_m[-1]), ha="center", va="bottom", fontsize=7.9, fontweight="bold", color=theme["ncbi"])
    ax_a.text(xs_w[-1], yl_m[-1] * 0.48, fmt_mem(yl_m[-1]), ha="center", va="top", fontsize=8.0, fontweight="bold", color=theme["losat"])
    ax_a.set_xscale("log"); ax_a.set_yscale("log"); ax_a.set_xlim(0.04, 2.2e6); ax_a.set_ylim(25, 2.5e6)
    ax_a.set_xticks(xs_w); ax_a.set_xticklabels([db_short[r["x_label"]] for r in warm], fontsize=7.3, color=theme["ink2"])
    ax_a.set_yticks([100, 1000, 10000, 100000]); ax_a.set_yticklabels(["100 MiB", "1 GiB", "10 GiB", "100 GiB"], fontsize=7.4, color=theme["ink2"])
    style_xy_grid(ax_a, theme, lw=0.6)
    ax_a.set_ylabel("Peak memory (max RSS)", fontsize=8.0, color=theme["ink2"])
    add_panel(ax_a, "a", "Database size \u2192 peak memory", theme, y=1.045, sub="megablast \u00b7 210 contigs \u00b7 16t \u00b7 -db")

    # --- Row 1, Panel b: DB size -> Wall time (warm & cold) ---
    ax_b = fig.add_axes([0.565, 0.656, 0.410, 0.210])
    ax_b.set_facecolor(theme["surface"])
    yn_w, yl_w = [float(r["ncbi_wall_s"]) for r in warm], [float(r["losat_wall_s"]) for r in warm]
    ax_b.fill_between(xs_w, yl_w, yn_w, color=theme["losat_fill"], alpha=0.45, zorder=1)
    ax_b.plot(xs_w, yn_w, marker="D", ms=5.4, color=theme["ncbi"], lw=1.9, zorder=3)
    ax_b.plot(xs_w, yl_w, marker="o", ms=6.0, color=theme["losat"], lw=2.1, zorder=4)
    warm_by = {r["x_label"]: r for r in warm}
    for rc in cold:
        rw = warm_by[rc["x_label"]]
        x = float(rc["x_value"])
        ax_b.plot([x, x], [float(rw["ncbi_wall_s"]), float(rc["ncbi_wall_s"])], color=theme["ncbi"], lw=1.1, ls="--", zorder=2)
        ax_b.scatter([x], [float(rc["ncbi_wall_s"])], marker="D", s=38, facecolors=theme["surface"], edgecolors=theme["ncbi"], lw=1.5, zorder=5)
        ax_b.scatter([x], [float(rc["losat_wall_s"])], marker="o", s=42, facecolors=theme["surface"], edgecolors=theme["losat"], lw=1.5, zorder=5)
    cnt_c = next(r for r in cold if r["x_label"] == "core_nt")
    rpr_c = next(r for r in cold if r["x_label"] == "ref_prok_rep_genomes")
    ax_b.annotate(f"cold: {fmt_fold(float(rpr_c['speedup']))}", xy=(float(rpr_c["x_value"]), float(rpr_c["ncbi_wall_s"])),
                  xytext=(float(rpr_c["x_value"]) * 0.11, float(rpr_c["ncbi_wall_s"]) * 2.2), fontsize=7.2, color=theme["ink2"], arrowprops=cold_arrow)
    ax_b.annotate(
        f"cold: {fmt_fold(float(cnt_c['speedup']))} ({fmt_sec(float(cnt_c['ncbi_wall_s']))} vs {fmt_sec(float(cnt_c['losat_wall_s']))})",
        xy=(float(cnt_c["x_value"]), float(cnt_c["ncbi_wall_s"])),
        xytext=(float(cnt_c["x_value"]) * 0.035, float(cnt_c["ncbi_wall_s"]) * 1.25),
        fontsize=7.2, color=theme["ink2"], ha="center", arrowprops=cold_arrow,
    )
    for i, (x, yn, yl) in enumerate(zip(xs_w, yn_w, yl_w)):
        ratio = yn / yl
        if i == 0:
            ax_b.text(x, math.sqrt(yn * yl), f"{fmt_fold(ratio)}*", ha="center", va="center", fontsize=7.3, color=theme["muted"], bbox=pill_box)
        elif i == 1:
            ax_b.text(x, yl * 0.44, fmt_fold(ratio), ha="center", va="top", fontsize=7.5, fontweight="semibold", color=theme["ink"])
        else:
            ax_b.text(x, yl * (0.42 if i == 2 else 0.30), f"warm {fmt_fold(ratio)}", ha="center", va="top", fontsize=7.3, fontweight="semibold", color=theme["losat"])
    ax_b.set_xscale("log"); ax_b.set_yscale("log"); ax_b.set_xlim(0.04, 2.2e6); ax_b.set_ylim(0.04, 3400)
    ax_b.set_xticks(xs_w); ax_b.set_xticklabels([db_short[r["x_label"]] for r in warm], fontsize=7.3, color=theme["ink2"])
    style_xy_grid(ax_b, theme, lw=0.6, plain_y=True)
    ax_b.set_ylabel("Wall time (s)", fontsize=8.0, color=theme["ink2"])
    add_panel(ax_b, "b", "Database size \u2192 wall time (warm & cold)", theme, y=1.045, sub=f"*abr_genes: {abr_nc_s:.2f} s start-up")

    # --- Row 2, Panel c: Query volume (twin mini-panels: time & memory) ---
    xs_q = [float(r["x_value"]) for r in q_rows]
    for ax, kn, kl, ylab, st in [
        (fig.add_axes([0.072, 0.375, 0.178, 0.195]), "ncbi_wall_s", "losat_wall_s", "Wall time (s)", "Wall time"),
        (fig.add_axes([0.296, 0.375, 0.178, 0.195]), "ncbi_rss_mib", "losat_rss_mib", "Max RSS (MiB)", "Peak memory"),
    ]:
        ax.set_facecolor(theme["surface"])
        yn, yl = [float(r[kn]) for r in q_rows], [float(r[kl]) for r in q_rows]
        ax.fill_between(xs_q, yl, yn, color=theme["losat_fill"], alpha=0.45, zorder=1)
        ax.plot(xs_q, yn, marker="D", ms=4.8, color=theme["ncbi"], lw=1.7, zorder=3)
        ax.plot(xs_q, yl, marker="o", ms=5.4, color=theme["losat"], lw=1.9, zorder=4)
        for x, nv, lv in zip(xs_q, yn, yl):
            ax.text(x, math.sqrt(nv * lv), fmt_fold(nv / lv), ha="center", va="center", fontsize=7.4, fontweight="semibold", color=theme["ink"], bbox=pill_box)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(5e3, 2.2e6); ax.set_ylim(min(yl) * 0.52, max(yn) * 2.3)
        ax.set_xticks([1e4, 1e5, 1e6]); ax.set_xticklabels(["10k", "100k", "1M"], fontsize=7.4, color=theme["ink2"])
        style_xy_grid(ax, theme, lw=0.6, plain_y=True)
        ax.set_xlabel("Query reads (150 nt)", fontsize=7.5, color=theme["ink2"])
        ax.set_ylabel(ylab, fontsize=7.8, color=theme["ink2"])
        ax.set_title(st, fontsize=7.9, fontweight="semibold", color=theme["ink2"], pad=5)
    fig.text(0.072, 0.602, "c  Query volume \u2192 time & memory", fontsize=9.4, fontweight="bold", color=theme["ink"])
    fig.text(0.474, 0.602, "megablast vs 16S \u00b7 16t", fontsize=7.3, ha="right", color=theme["muted"])

    # --- Row 2, Panel d: Multi-core speed efficiency (1 -> 16 cores vs 45° ideal diagonal) ---
    ax_d = fig.add_axes([0.565, 0.375, 0.410, 0.215])
    ax_d.set_facecolor(theme["surface"])
    ax_d.plot([1, 16], [1, 16], ls=":", color=theme["muted"], lw=1.3, zorder=1)
    mode_styles = [
        ("bp",  "blastp",  theme["losat"],   "o"),
        ("tbn", "tblastn", theme["wasm"],    "^"),
        ("bn",  "blastn",  theme["diamond"], "s"),
    ]
    t_series = [
        (cid, nm, col, mk, sorted([r for r in t_rows if r["case"] == cid], key=lambda x: int(x["x_value"])))
        for cid, nm, col, mk in mode_styles
    ]
    for cid, nm, col, mk, sub in t_series:
        xs_t = [int(r["x_value"]) for r in sub]
        t1 = float(sub[0]["losat_wall_s"])
        sps = [t1 / float(r["losat_wall_s"]) for r in sub]
        ax_d.plot(xs_t, sps, marker=mk, ms=5.5, color=col, lw=2.0, zorder=4)
        dy = -0.35 if cid == "bn" else 0.0
        ax_d.text(16.45, sps[-1] + dy, f"{nm}: {sps[-1]:.1f}\u00d7 ({sps[-1] / 16 * 100:.0f}%)",
                  ha="left", va="center", fontsize=7.5, fontweight="bold", color=col)
    ax_d.text(15.2, 16.3, "Ideal linear (16\u00d7)", ha="right", va="bottom", fontsize=7.1, style="italic", color=theme["muted"])
    ax_d.set_xlim(0.5, 23.0); ax_d.set_ylim(0.5, 18.0)
    ax_d.set_xticks([1, 2, 4, 8, 16]); ax_d.set_yticks([1, 4, 8, 12, 16])
    ax_d.set_yticklabels(["1\u00d7", "4\u00d7", "8\u00d7", "12\u00d7", "16\u00d7"], fontsize=7.4, color=theme["ink2"])
    style_xy_grid(ax_d, theme, lw=0.6, clear_mx=False)
    ax_d.set_xlabel("Threads (cores)", fontsize=7.8, color=theme["ink2"])
    ax_d.set_ylabel("Self-speedup vs 1t", fontsize=7.9, color=theme["ink2"])
    add_panel(ax_d, "d", f"Thread scaling: speed (1 \u2192 {t_max} cores)", theme, y=1.045, sub=f"{t_sp_lo:.1f}\u2013{t_sp_hi:.1f}\u00d7 faster than BLAST+ at every N")

    # --- Row 3, Panel e (NEW): Thread scaling -> Peak memory curves (1 -> 16 cores, LOSAT vs BLAST+) ---
    ax_e = fig.add_axes([0.072, 0.095, 0.402, 0.215])
    ax_e.set_facecolor(theme["surface"])
    for _cid, nm, col, mk, sub in t_series:
        xs_m = [int(r["x_value"]) for r in sub]
        yl_r = [float(r["losat_rss_mib"]) for r in sub]
        yn_r = [float(r["ncbi_rss_mib"]) for r in sub]
        ax_e.plot(xs_m, yn_r, marker="D", ms=4.2, color=col, ls="--", lw=1.4, alpha=0.65, zorder=2)
        ax_e.plot(xs_m, yl_r, marker=mk, ms=5.4, color=col, ls="-", lw=2.1, zorder=4)
        ax_e.text(16.4, yl_r[-1], f"LOSAT {nm} ({fmt_mem(yl_r[-1], True)})", ha="left", va="center", fontsize=7.1, fontweight="bold", color=col)
        ax_e.text(16.4, yn_r[-1], f"BLAST+ {nm} ({fmt_mem(yn_r[-1], True)})", ha="left", va="center", fontsize=6.8, color=theme["muted"])
    ax_e.set_yscale("log"); ax_e.set_xlim(0.5, 25.5); ax_e.set_ylim(20.0, 12000.0)
    ax_e.set_xticks([1, 2, 4, 8, 16])
    ax_e.set_yticks([30, 100, 300, 1000, 5000]); ax_e.set_yticklabels(["30M", "100M", "300M", "1G", "5G"], fontsize=7.4, color=theme["ink2"])
    style_xy_grid(ax_e, theme, lw=0.6, clear_mx=False)
    ax_e.set_xlabel("Threads (cores)", fontsize=7.8, color=theme["ink2"])
    ax_e.set_ylabel("Peak memory (max RSS)", fontsize=7.9, color=theme["ink2"])
    add_panel(ax_e, "e", "Thread scaling: peak memory (1 \u2192 16 cores)", theme, y=1.045, sub="Solid = LOSAT \u00b7 Dashed = BLAST+")

    # --- Row 3, Panel f (NEW): All 12 search modes — 1t -> 16t peak memory shift ---
    ax_f = fig.add_axes([0.620, 0.095, 0.355, 0.215])
    ax_f.set_facecolor(theme["surface"])
    yf = [float(len(all12) - 1 - i) for i in range(len(all12))]
    ax_f.set_ylim(-0.7, len(all12) - 0.3)
    ax_f.set_xlim(10.0, 26000.0)
    for i, (yv, cid) in enumerate(zip(yf, all12)):
        if i % 2 == 1:
            ax_f.axhspan(yv - 0.48, yv + 0.48, color=theme["band"], zorder=0)
        r1, r16 = by_kt[(cid, 1)], by_kt[(cid, 16)]
        lm1, lm16 = float(r1["losat_rss_mib"]), float(r16["losat_rss_mib"])
        nm1, nm16 = float(r1["ncbi_rss_mib"]), float(r16["ncbi_rss_mib"])
        # Offset slightly vertically so LOSAT (bottom) & BLAST+ (top) 1t->16t segments never hide each other
        ax_f.plot([nm1, nm16], [yv + 0.16, yv + 0.16], color=theme["ncbi"], lw=1.4, alpha=0.65, zorder=2)
        ax_f.scatter([nm1], [yv + 0.16], marker="D", s=16, facecolors=theme["surface"], edgecolors=theme["ncbi"], lw=1.0, zorder=3)
        ax_f.scatter([nm16], [yv + 0.16], marker="D", s=20, color=theme["ncbi"], zorder=4)

        ax_f.plot([lm1, lm16], [yv - 0.16, yv - 0.16], color=theme["losat"], lw=1.8, zorder=3)
        ax_f.scatter([lm1], [yv - 0.16], marker="o", s=18, facecolors=theme["surface"], edgecolors=theme["losat"], lw=1.1, zorder=4)
        ax_f.scatter([lm16], [yv - 0.16], marker="o", s=24, color=theme["losat"], zorder=5)
        ax_f.text(max(nm1, nm16) * 1.18, yv, f"{fmt_mem(lm16, True)} / {fmt_mem(nm16, True)}", ha="left", va="center", fontsize=6.2, color=theme["muted"])

    ax_f.set_yticks(yf)
    ax_f.set_yticklabels([by_kt[(cid, 16)]["label"] for cid in all12], fontsize=6.9, fontweight="semibold", color=theme["ink"])
    ax_f.tick_params(axis="y", length=0, pad=4)
    style_log_x(ax_f, [20, 100, 500, 2000, 5000], ["20M", "100M", "500M", "2G", "5G"], theme, tick_fs=7.3)
    ax_f.set_xlabel("\u2190 Smaller   Peak memory: open = 1t, filled = 16t (max RSS)", fontsize=7.4, color=theme["ink2"])
    add_panel(ax_f, "f", "All 12 modes: 1t \u2192 16t memory shift", theme, y=1.045, sub="16t: LOSAT / BLAST+")

    add_footer(
        fig,
        f"LOSAT {short_sha(meta.get('losat_commit', ''))} (core_nt & cold-cache runs: {short_sha(meta.get('toolmt_commit', ''))}) vs NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')}, {short_hw(meta)}.\n"
        f"Panels a\u2013c: both programs ran -db. Panels d\u2013f: Table 1 forms (-db for blastx/tblastn/psiblast/rps*; LOSAT -subject vs BLAST+ -db for others). BLAST+ max RSS includes memory-mapped DB pages.",
        theme,
        y_rule=0.042,
    )
    return save_themed(fig, out, "fig3_scaling", dark, theme, dpi=180)


# ---------------------------------------------------------------------------
# Figure 4 (Pattern 4A): Memory-Sized Bubble Pareto (UniRef50) + Horizontal Dumbbells (SILVA 144)
# ---------------------------------------------------------------------------
def render_fig4(snap: Path, out: Path, theme: dict[str, object], dark: bool) -> Path:
    ur50 = read_tsv(snap / "vs_heuristics_uniref50.tsv")
    silva = read_tsv(snap / "vs_heuristics_silva.tsv")
    meta = read_meta(snap)
    by_sc = {(r["case"], r["arm"]): r for r in silva}
    tcol = {
        "LOSAT": theme["losat"],
        "NCBI BLAST+": theme["ncbi"],
        "DIAMOND": theme["diamond"],
        "MMseqs2": theme["mmseqs"],
        "vsearch": theme["vsearch"],
    }

    def _gib_rng(tool: str) -> tuple[float, float]:
        vals = [float(r["max_rss_gib"]) for r in ur50 if r["tool"] == tool]
        return min(vals), max(vals)

    lo_lo, lo_hi = _gib_rng("LOSAT")
    nc_lo, nc_hi = _gib_rng("NCBI BLAST+")
    dm_lo, dm_hi = _gib_rng("DIAMOND")
    mm_lo, mm_hi = _gib_rng("MMseqs2")
    vs_gib = float(by_sc[("amp1e3", "vsearch")]["max_rss_mib"]) / 1024.0

    fig = new_fig((14.2, 8.2), theme)
    add_header(
        fig,
        "Exact NCBI BLAST+ parity vs. heuristic aligners: speed, sensitivity, and memory trade-offs",
        f"UniRef50 (38.8M proteins, 32t): bubble area scales with \u221a(max RSS) + min size ({lo_lo:.2f} GiB LOSAT vs {mm_hi:.0f} GiB MMseqs2; key at lower right of a)  \u00b7  SILVA 144 (905k rRNAs, 16t): dumbbells",
        theme,
        y_title=0.974,
        y_sub=0.943,
    )
    handles = [
        Line2D([0], [0], marker="o", ls="none", ms=6.2, color=theme["losat"], label=f"LOSAT (exact NCBI output, {lo_lo:.2f}\u2013{lo_hi:.2f} GiB)"),
        Line2D([0], [0], marker="o", ls="none", ms=7.2, color=theme["ncbi"], alpha=0.65, label=f"NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')} ({nc_lo:.0f}\u2013{nc_hi:.0f} GiB)"),
        Line2D([0], [0], marker="o", ls="none", ms=6.8, color=theme["diamond"], alpha=0.65, label=f"DIAMOND {meta.get('diamond_version', '2.2.8')} ({dm_lo:.1f}\u2013{dm_hi:.1f} GiB)"),
        Line2D([0], [0], marker="o", ls="none", ms=9.5, color=theme["mmseqs"], alpha=0.35, markeredgecolor=theme["mmseqs"], label=f"MMseqs2 {meta.get('mmseqs2_version', '18.8cc5c')} ({mm_lo:.0f}\u2013{mm_hi:.0f} GiB)"),
        Line2D([0], [0], marker="P", ls="none", ms=6.6, color=theme["vsearch"], label=f"vsearch {meta.get('vsearch_version', '2.32.0')} ({vs_gib:.1f} GiB)"),
    ]
    add_legend(fig, handles)

    def mem_area(gib: float) -> float:
        return 30.0 + 360.0 * math.sqrt(max(gib, 0.05) / 20.0)

    arm_short = {
        "losat_bp": "LOSAT blastp",
        "losat_fast": "LOSAT blastp-fast",
        "ncbi_def": "NCBI blastp",
        "ncbi_fast": "NCBI fast",
        "dmnd_def": "DIAMOND default",
        "dmnd_ms": "DIAMOND more",
        "dmnd_us": "DIAMOND ultra",
        "mmseqs_s57_idx": "MMseqs2 -s5.7 idx",
        "mmseqs_s75_idx": "MMseqs2 -s7.5 idx",
        "mmseqs_s75_noidx": "MMseqs2 -s7.5",
    }

    def _ur50_txt(aid: str, r: dict[str, str]) -> str:
        ws, gib = float(r["wall_s"]), float(r["max_rss_gib"])
        g_str = f"{gib:.2f} GiB" if gib < 1.0 else (f"{gib:.0f} GiB" if gib >= 90.0 else f"{gib:.1f}G")
        t_str = f"{round(ws):,.0f}s"
        return f"{arm_short[aid]}\n{t_str} \u00b7 {g_str}"

    # Point-relative offsets in points (dx_pt, dy_pt) so labels move with points on snapshot rebuilds
    offsets_100 = [
        ("losat_bp",         -15.4,  -27.9, "right"),
        ("ncbi_def",           9.0,  -27.9, "center"),
        ("dmnd_us",           21.8,  -45.7, "center"),
        ("mmseqs_s75_idx",   -41.2,   35.6, "center"),
        ("losat_fast",       -47.0,  -45.7, "left"),
        ("mmseqs_s57_idx",   -14.5,  -64.2, "left"),
        ("ncbi_fast",        -25.8, -118.4, "right"),
        ("dmnd_ms",            4.6, -143.2, "center"),
        ("mmseqs_s75_noidx",  38.6,  -77.4, "center"),
        ("dmnd_def",          18.9,    0.0, "left"),
    ]
    offsets_1000 = [
        ("losat_bp",         -33.7,  -27.1, "right"),
        ("ncbi_def",           2.7,  -32.5, "center"),
        ("dmnd_us",           30.9,   -7.7, "left"),
        ("losat_fast",         4.2,   63.5, "center"),
        ("mmseqs_s75_idx",   -39.7,   48.0, "center"),
        ("mmseqs_s57_idx",   -14.5,  -43.3, "left"),
        ("dmnd_ms",          -31.0, -102.2, "center"),
        ("ncbi_fast",         48.3,  -41.8, "left"),
        ("mmseqs_s75_noidx",  22.2, -110.7, "center"),
        ("dmnd_def",          18.3,    0.0, "left"),
    ]

    for p_i, (qset, rect, let, title, xlim, olist) in enumerate([
        ("q100",  [0.056, 0.135, 0.278, 0.645], "a", "UniRef50: 100q (area: \u221amax RSS)",   (7.0, 3600.0),  offsets_100),
        ("q1000", [0.362, 0.135, 0.278, 0.645], "b", "UniRef50: 1,000q (area: \u221amax RSS)", (14.0, 32000.0), offsets_1000),
    ]):
        ax = fig.add_axes(rect)
        ax.set_facecolor(theme["surface"])
        sub = [r for r in ur50 if r["query_set"] == qset]
        ax.axhspan(98.6, 103.8, color=theme["zone"], alpha=0.60, zorder=0)
        ax.axhline(100.0, color=theme["losat"], lw=0.9, ls="--", alpha=0.7, zorder=1)
        ax.text(xlim[0] * 1.15, 100.8, "100% exact NCBI BLAST+ parity", fontsize=6.9, fontweight="semibold", color=theme["losat"], va="bottom")

        by_arm = {r["arm_id"]: r for r in sub}
        for n_id, l_id in (("ncbi_def", "losat_bp"), ("ncbi_fast", "losat_fast")):
            rn, rl = by_arm[n_id], by_arm[l_id]
            ax.annotate("", xy=(float(rl["wall_s"]), float(rl["recall_uncapped"]) * 100),
                        xytext=(float(rn["wall_s"]), float(rn["recall_uncapped"]) * 100),
                        arrowprops=dict(arrowstyle="-|>", color=theme["losat"], lw=1.4, shrinkA=8, shrinkB=6), zorder=4)

        for r in sorted(sub, key=lambda x: float(x["max_rss_gib"]), reverse=True):
            x, y = float(r["wall_s"]), float(r["recall_uncapped"]) * 100.0
            gib = float(r["max_rss_gib"])
            col = tcol[r["tool"]]
            is_lo = r["tool"] == "LOSAT"
            ax.scatter([x], [y], s=mem_area(gib), facecolors=col, alpha=0.94 if is_lo else (0.30 if dark else 0.22),
                       edgecolors=col, linewidths=1.5 if is_lo else 1.2, zorder=5 if is_lo else 3)
            if not is_lo:
                # Draw non-LOSAT centre dots above LOSAT's dot with a crisp halo so near-coincident points
                # (e.g. DIAMOND --more-sensitive vs LOSAT blastp-fast in panel b) remain both clearly visible.
                ax.scatter([x], [y], s=24, facecolors=col, edgecolors=theme["surface"], linewidths=0.9, zorder=6)

        for aid, dx_pt, dy_pt, ha in olist:
            r = by_arm[aid]
            is_lo = r["tool"] == "LOSAT"
            ax.annotate(
                _ur50_txt(aid, r), xy=(float(r["wall_s"]), float(r["recall_uncapped"]) * 100),
                xytext=(dx_pt, dy_pt), textcoords="offset points",
                ha=ha, va="center", fontsize=6.8 if is_lo else 6.4, fontweight="bold" if is_lo else "normal",
                color=theme["losat"] if is_lo else theme["ink2"], linespacing=1.15,
                bbox=dict(boxstyle="round,pad=0.13", facecolor=theme["surface"], edgecolor=theme["zone_edge"] if is_lo else "none", lw=0.6 if is_lo else 0.0, alpha=0.90),
                arrowprops=dict(arrowstyle="-", color=theme["losat"] if is_lo else theme["hairline"], lw=0.7, shrinkA=2, shrinkB=4 if aid == "dmnd_ms" else 5),
                zorder=8 if is_lo else 7,
            )

        # Bubble-size scale reference inside bottom-right empty corner of Panel a
        if p_i == 0:
            ax.text(1050, 64.2, "Key: \u221a(max RSS)", ha="center", va="bottom", fontsize=6.6, fontweight="semibold", color=theme["muted"])
            for bx, gib_ref, blbl in [(450, 0.2, "0.2G"), (1000, 10.0, "10G"), (2150, 100.0, "100G")]:
                ax.scatter([bx], [60.2], s=mem_area(gib_ref), facecolors=theme["muted"], alpha=0.22, edgecolors=theme["muted"], lw=1.0, zorder=5)
                ax.text(bx, 57.0, blbl, ha="center", va="top", fontsize=6.3, color=theme["muted"])

        ax.set_xscale("log"); ax.set_xlim(*xlim); ax.set_ylim(55.0, 104.2)
        ax.xaxis.set_major_formatter(FuncFormatter(plain_log))
        style_xy_grid(ax, theme, lw=0.6)
        ax.set_xlabel("\u2190 Faster   Wall time (s, log scale)", fontsize=7.9, color=theme["ink2"])
        if p_i == 0:
            ax.set_ylabel("Recall of NCBI blastp hits (%, uncapped queries)  \u2191", fontsize=8.0, color=theme["ink2"])
        else:
            ax.set_yticklabels([])
        add_panel(ax, let, title, theme, y=1.022)

    # Panel c: SILVA 144 horizontal dumbbells (Top = Wall Time, Bottom = Memory with RssAnon solid + max RSS open!)
    fig.add_artist(Line2D([0.664, 0.664], [0.105, 0.818], transform=fig.transFigure, color=theme["grid"], lw=0.9))
    sc_cases = [("reads1e6", "1M reads\n(megablast)"), ("amp1e3", "1k amplicons\n(QIIME 2)"), ("amp1e4", "10k amplicons\n(QIIME 2)*")]
    ys_s = [2.0, 1.0, 0.0]

    for p_idx, (ax, metric, xlim_s, ticks_s, lbls_s, let_c, hdr, xlab) in enumerate([
        (fig.add_axes([0.752, 0.495, 0.226, 0.285]), "wall_s",      (12.0, 5200.0),  [20, 100, 500, 2000],        ["20s", "100s", "500s", "2k"],        "c", "SILVA 144: Wall time (16t)", "\u2190 Faster   Wall time (s)"),
        (fig.add_axes([0.752, 0.135, 0.226, 0.270]), "max_rss_mib", (12.0, 68000.0), [20, 100, 500, 2000, 15000], ["20M", "100M", "500M", "2G", "15G"], "",  "SILVA 144: Memory (solid = RssAnon)", "\u2190 Smaller   Memory (MiB)"),
    ]):
        ax.set_facecolor(theme["surface"])
        ax.set_ylim(-0.62, 2.65); ax.set_xlim(*xlim_s)
        for i, (yv, (cid, _lbl)) in enumerate(zip(ys_s, sc_cases)):
            if i % 2 == 1:
                ax.axhspan(yv - 0.45, yv + 0.45, color=theme["band"], zorder=0)
            arms = ["losat", "ncbi", "vsearch"] if cid == "amp1e3" else ["losat", "ncbi"]
            vals = [float(by_sc[(cid, a)][metric]) for a in arms]
            all_pts = list(vals)
            if p_idx == 1:
                all_pts.extend(float(by_sc[(cid, a)]["rss_anon_mib"]) for a in arms)
            ax.plot([min(all_pts), max(all_pts)], [yv, yv], color=theme["hairline"], lw=1.7, zorder=2)
            for a in arms:
                r = by_sc[(cid, a)]
                v = float(r[metric])
                mk, col = ("o", theme["losat"]) if a == "losat" else (("D", theme["ncbi"]) if a == "ncbi" else ("P", theme["vsearch"]))
                if p_idx == 1:
                    va = float(r["rss_anon_mib"])
                    ax.scatter([v], [yv], marker=mk, s=48 if a == "losat" else 42, facecolors=theme["surface"], edgecolors=col, linewidths=1.5, zorder=4)
                    ax.scatter([va], [yv], marker=mk, s=42 if a == "losat" else 34, color=col, zorder=5 if a == "losat" else 4)
                else:
                    ax.scatter([v], [yv], marker=mk, s=50 if a == "losat" else 42, color=col, zorder=5 if a == "losat" else 4)

            lo_v, nc_v = float(by_sc[(cid, "losat")][metric]), float(by_sc[(cid, "ncbi")][metric])
            ax.text(lo_v, yv + 0.25, fmt_sec(lo_v) if p_idx == 0 else fmt_mem(lo_v, True),
                    ha="center", va="bottom", fontsize=6.9, fontweight="bold", color=theme["losat"])
            ax.text(nc_v, yv - 0.26, fmt_sec(nc_v) if p_idx == 0 else fmt_mem(nc_v, True),
                    ha="center", va="top", fontsize=6.7, color=theme["ink2"])
            if cid == "amp1e3":
                vs_v = float(by_sc[(cid, "vsearch")][metric])
                ax.text(vs_v, yv + 0.25, f"vsearch {fmt_sec(vs_v) if p_idx == 0 else fmt_mem(vs_v, True)}",
                        ha="center", va="bottom", fontsize=6.5, fontweight="semibold", color=theme["vsearch"])
        ax.set_yticks(ys_s)
        ax.set_yticklabels([x[1] for x in sc_cases], fontsize=7.2, color=theme["ink"])
        ax.tick_params(axis="y", length=0, pad=5)
        style_log_x(ax, ticks_s, lbls_s, theme, tick_fs=7.1)
        ax.set_xlabel(xlab, fontsize=7.5, color=theme["ink2"])
        add_panel(ax, let_c, hdr, theme, y=1.035)

    amp1e4_drift = int(by_sc[("amp1e4", "ncbi")]["rows"]) - int(by_sc[("amp1e4", "losat")]["rows"])
    m_ov = re.search(r"(\d+)pct", by_sc[("amp1e3", "vsearch")]["parity"])
    ov_str = f"{m_ov.group(1)}% hit overlap" if m_ov else by_sc[("amp1e3", "vsearch")]["parity"]
    add_footer(
        fig,
        f"UniRef50 (38.8M seqs, 32t, node 1): LOSAT {short_sha(meta.get('uniref50_losat_commit', ''))} vs NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')}, "
        f"DIAMOND {meta.get('diamond_version', '2.2.8')}, MMseqs2 {meta.get('mmseqs2_version', '18.8cc5c')} (all Bioconda; index = 101 GiB).\n"
        f"SILVA 144 NR99 (905k seqs, 16t): LOSAT {short_sha(meta.get('silva_losat_commit', ''))} vs NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')} & vsearch {meta.get('vsearch_version', '2.32.0')} "
        f"(QIIME 2 --maxrejects 0; {ov_str}; open marker = max RSS, solid = RssAnon). "
        f"*On 10k amplicons, NCBI 16t drifts ({amp1e4_drift:+d} rows); LOSAT matches 1-thread NCBI byte-for-byte.",
        theme,
        y_rule=0.058,
    )
    return save_themed(fig, out, "fig4_vs_heuristics", dark, theme)


# ---------------------------------------------------------------------------
# Figure 5 (Pattern 5A): WebAssembly 3-Engine Physical-Scale Dumbbells
# ---------------------------------------------------------------------------
def render_fig5(snap: Path, out: Path, theme: dict[str, object], dark: bool) -> Path:
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
            y -= 0.52
        y_rows.append(("hdr", y, gtitle))
        y -= 0.88
        for cid in cids:
            y_rows.append(("case", y, cid))
            y -= 1.04
    ymin, ymax = y + 0.30, 0.85

    fig = new_fig((13.6, 7.0), theme)
    rows16_m = [r for r in rows if int(r["threads"]) == 16 and r["ncbi_s"]]
    clear_wins16 = [r for r in rows16_m if float(r["ncbi_s"]) / float(r["wasm_s"]) >= 1.05]
    near_eq16 = [r for r in rows16_m if 0.95 <= float(r["ncbi_s"]) / float(r["wasm_s"]) < 1.05]
    n_wins16 = sum(1 for r in rows16_m if float(r["wasm_s"]) < float(r["ncbi_s"]))
    ncbi_ver = meta.get("ncbi_version", "2.17.0")
    if wasm_derived:
        eq_note = f" ({', '.join(r['mode_label'] for r in near_eq16)} ~equal)" if near_eq16 else ""
        fig5_title = (
            f"Inside WebAssembly, LOSAT is estimated faster than native NCBI BLAST+ {ncbi_ver} "
            f"in {len(clear_wins16)} of {len(rows16_m)} modes at 16t{eq_note}"
        )
        sub_extra = f"  \u00b7  Estimated times scaled from build {base_sha} via side-by-side ratio"
    else:
        fig5_title = f"Inside WebAssembly, LOSAT is faster than native NCBI BLAST+ {ncbi_ver} in {n_wins16} of {len(rows16_m)} measured modes at 16t"
        sub_extra = "  \u00b7  Direct paired runs under Node.js"
    add_header(
        fig,
        fig5_title,
        f"Physical-scale comparison (left = faster / smaller): Native LOSAT vs. wasm32-wasip1-threads (Node.js) vs. Native NCBI BLAST+ {ncbi_ver}{sub_extra}",
        theme,
        y_title=0.968,
        y_sub=0.932,
    )
    handles = [
        Line2D([0], [0], marker="o", ls="none", ms=6.8, color=theme["losat"], label="LOSAT native"),
        Line2D([0], [0], marker="^", ls="none", ms=7.0, color=theme["wasm"], label="LOSAT WebAssembly (wasm32-wasip1-threads, Node.js)"),
        Line2D([0], [0], marker="D", ls="none", ms=6.0, color=theme["ncbi"], label="Native NCBI BLAST+ 2.17.0"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.018, 0.888), frameon=False, ncol=3, fontsize=7.9, columnspacing=1.5)

    bot_y, h_ax = 0.130, 0.665
    ax_lbl = fig.add_axes([0.020, bot_y, 0.132, h_ax])
    ax_lbl.set_facecolor(theme["surface"]); ax_lbl.set_xlim(0, 1); ax_lbl.set_ylim(ymin, ymax); ax_lbl.axis("off")
    ri = 0
    for kind, yv, txt in y_rows:
        if kind == "hdr":
            ax_lbl.text(0.0, yv, str(txt).upper(), ha="left", va="center", fontsize=7.5, fontweight="bold", color=theme["losat"])
        else:
            if ri % 2 == 1:
                ax_lbl.axhspan(yv - 0.47, yv + 0.47, color=theme["band"], zorder=0)
            ri += 1
            ax_lbl.text(0.04, yv, by_kt[(str(txt), 16)]["mode_label"], ha="left", va="center", fontsize=8.7, fontweight="bold", color=theme["ink"])

    for plot_rect, pct_rect, let, title, th, xlim, ticks, tlbls in [
        ([0.162, bot_y, 0.205, h_ax], [0.370, bot_y, 0.052, h_ax], "a", "Wall time, 16t (s)", 16, (0.35, 52.0), [0.5, 2, 5, 15, 35], ["0.5", "2", "5", "15", "35"]),
        ([0.442, bot_y, 0.205, h_ax], [0.650, bot_y, 0.052, h_ax], "b", "Wall time, 1t (s)",  1,  (5.0, 420.0), [10, 30, 100, 300],   ["10", "30", "100", "300"]),
    ]:
        ax, ax_p = fig.add_axes(plot_rect), fig.add_axes(pct_rect)
        ax.set_facecolor(theme["surface"]); ax.set_xlim(*xlim); ax.set_ylim(ymin, ymax); ax.set_yticks([])
        ax_p.set_facecolor(theme["surface"]); ax_p.set_xlim(0, 1); ax_p.set_ylim(ymin, ymax); ax_p.axis("off")
        ri = 0
        for kind, yv, _ in y_rows:
            if kind == "case":
                if ri % 2 == 1:
                    ax.axhspan(yv - 0.47, yv + 0.47, color=theme["band"], zorder=0)
                    ax_p.axhspan(yv - 0.47, yv + 0.47, color=theme["band"], zorder=0)
                ri += 1
        style_log_x(ax, ticks, tlbls, theme, tick_fs=7.5)
        ax.set_xlabel("\u2190 Faster   Wall time (s)", fontsize=7.7, color=theme["ink2"])
        add_panel(ax, let, title, theme, y=1.024)
        ax_p.text(0.5, ymax - 0.12, "Wasm /\nnative", ha="center", va="center", fontsize=7.0, color=theme["muted"])

        for kind, yv, txt in y_rows:
            if kind != "case":
                continue
            r = by_kt[(str(txt), th)]
            nl, wl = float(r["native_s"]), float(r["wasm_s"])
            bl = float(r["ncbi_s"]) if r["ncbi_s"] else None
            pts = [nl, wl] + ([bl] if bl else [])
            ax.plot([min(pts), max(pts)], [yv, yv], color=theme["hairline"], lw=1.8, zorder=2)
            if bl:
                ax.scatter([bl], [yv], marker="D", s=40, color=theme["ncbi"], zorder=3)
            ax.scatter([wl], [yv], marker="^", s=52, color=theme["wasm"], zorder=4)
            ax.scatter([nl], [yv], marker="o", s=48, color=theme["losat"], zorder=5)
            ax_p.text(0.5, yv, f"{float(r['wasm_pct_of_native']):.0f}%", ha="center", va="center", fontsize=8.1, fontweight="semibold", color=theme["wasm"])

    ax_c = fig.add_axes([0.728, bot_y, 0.252, h_ax])
    ax_c.set_facecolor(theme["surface"]); ax_c.set_xlim(10.0, 12000.0); ax_c.set_ylim(ymin, ymax); ax_c.set_yticks([])
    ri = 0
    for kind, yv, _ in y_rows:
        if kind == "case":
            if ri % 2 == 1:
                ax_c.axhspan(yv - 0.47, yv + 0.47, color=theme["band"], zorder=0)
            ri += 1
    style_log_x(ax_c, [20, 100, 500, 2000, 5000], ["20", "100", "500", "2k", "5k"], theme, tick_fs=7.5)
    ax_c.set_xlabel("\u2190 Smaller   Peak memory, 16t (MiB)", fontsize=7.7, color=theme["ink2"])
    add_panel(ax_c, "c", "Peak memory, 16t (MiB)", theme, y=1.024)
    for kind, yv, txt in y_rows:
        if kind != "case":
            continue
        r = by_kt[(str(txt), 16)]
        nm, wm = float(r["native_mib"]), float(r["wasm_mib"])
        bm = float(r["ncbi_mib"]) if r["ncbi_mib"] else None
        pts = [nm, wm] + ([bm] if bm else [])
        ax_c.plot([min(pts), max(pts)], [yv, yv], color=theme["hairline"], lw=1.6, zorder=2)
        if bm:
            ax_c.scatter([bm], [yv], marker="D", s=40, color=theme["ncbi"], zorder=3)
        ax_c.scatter([wm], [yv], marker="^", s=52, color=theme["wasm"], zorder=4)
        ax_c.scatter([nm], [yv], marker="o", s=48, color=theme["losat"], zorder=5)

    wasm_sha = short_sha(meta.get("wasm_commit") or "")
    derived_note = (
        f"Note: {wasm_sha} times are scaled from measured build {base_sha} by a side-by-side ratio run; NCBI BLAST+ tblastx was not timed in this set."
        if wasm_derived
        else "Direct paired native, WebAssembly, and NCBI BLAST+ runs in the same lease."
    )
    add_footer(
        fig,
        f"LOSAT {wasm_sha} (native & wasm32-wasip1-threads under Node.js) vs native NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')}, {short_hw(meta)}. "
        f"Wasm / native % = native time \u00f7 Wasm time (100% = native parity).\n{derived_note}",
        theme,
        y_rule=0.064,
    )
    return save_themed(fig, out, "fig5_webassembly", dark, theme)


# ---------------------------------------------------------------------------
# Figure 6 (Pattern 6A): 3-Stage Generational Slopegraph (v0.1.0 -> NCBI BLAST+ -> v0.2.0)
# ---------------------------------------------------------------------------
def render_fig6(snap: Path, out: Path, theme: dict[str, object], dark: bool) -> Path:
    _std, by_kt = load_std_by_kt(snap)
    meta = read_meta(snap)

    cases = [
        ("bn",  "blastn",     theme["losat"],   "o"),
        ("tbx", "tblastx*",   theme["wasm"],    "^"),
        ("mb",  "megablast*", theme["diamond"], "s"),
        ("bp",  "blastp",     theme["mmseqs"],  "D"),
    ]

    max_sp_v = max(float(by_kt[(cid, 16)]["v010_wall_s"]) / float(by_kt[(cid, 16)]["losat_wall_s"]) for cid, *_ in cases)
    max_rm_v = max(float(by_kt[(cid, 16)]["v010_rss_mib"]) / float(by_kt[(cid, 16)]["losat_rss_mib"]) for cid, *_ in cases)
    fig = new_fig((12.6, 6.2), theme)
    add_header(
        fig,
        f"LOSAT v0.1.0 \u2192 v{meta.get('losat_version', '0.2.0')}: up to {max_sp_v:.1f}\u00d7 faster and {max_rm_v:.1f}\u00d7 lower memory, below NCBI BLAST+ in all {len(cases)} modes",
        f"All {len(cases)} search modes present in v0.1.0 (16 threads) alongside NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')}  \u00b7  v{meta.get('losat_version', '0.2.0')} also expands from 3 to 17 programs and adds -db support",
        theme,
        y_title=0.964,
        y_sub=0.918,
    )

    for p_i, (rect, let, title, k0, kn, kl, ylab, ylim, yticks, ylbls, fmt_fn) in enumerate([
        ([0.078, 0.160, 0.392, 0.655], "a", "Wall time, 16 threads (log s)", "v010_wall_s", "ncbi_wall_s", "losat_wall_s",
         "\u2190 Faster   Wall time (s)", (1.5, 6500.0), [2, 10, 50, 200, 1000, 3000], ["2s", "10s", "50s", "200s", "1,000s", "3,000s"], fmt_sec),
        ([0.556, 0.160, 0.406, 0.655], "b", "Peak memory, 16 threads (log MiB)", "v010_rss_mib", "ncbi_rss_mib", "losat_rss_mib",
         "\u2190 Smaller   Peak memory (max RSS)", (25.0, 8000.0), [30, 100, 300, 1000, 4000], ["30M", "100M", "300M", "1 GiB", "4 GiB"], fmt_mem),
    ]):
        ax = fig.add_axes(rect)
        ax.set_facecolor(theme["surface"])
        xs = [0.0, 1.0, 2.0]
        for x in xs:
            ax.axvline(x, color=theme["grid"], lw=1.1, zorder=1)
        ax.axvspan(1.65, 2.65, color=theme["zone"], alpha=0.52, zorder=0)

        r_pos = {"tbx": 9.4, "bn": 5.6, "mb": 3.5, "bp": 2.15} if p_i == 0 else {"tbx": 567.0, "bn": 86.0, "mb": 56.0, "bp": 37.0}
        nudge_l = {"bn": 1.18, "tbx": 0.78, "mb": 1.0, "bp": 1.0} if p_i == 0 else {"bn": 1.0, "tbx": 0.85, "mb": 1.15, "bp": 1.0}

        for cid, lbl, col, mk in cases:
            r = by_kt[(cid, 16)]
            y0, yn, yl = float(r[k0]), float(r[kn]), float(r[kl])
            ax.plot(xs, [y0, yn, yl], color=col, lw=2.3, alpha=0.90, zorder=3)
            ax.scatter(xs, [y0, yn, yl], marker=mk, s=[48, 42, 62], color=col, edgecolors=theme["surface"], linewidths=0.8, zorder=5)
            ax.text(-0.07, y0 * nudge_l[cid], f"{lbl}\n{fmt_fn(y0)}", ha="right", va="center", fontsize=7.4, fontweight="semibold", color=col)
            ax.annotate(
                f"{lbl}: {fmt_fn(yl)} ({fmt_fold(y0 / yl, True)})",
                xy=(2.0, yl), xytext=(2.14, r_pos[cid]), ha="left", va="center",
                fontsize=7.4, fontweight="bold", color=col,
                arrowprops=dict(arrowstyle="-", color=col, alpha=0.45, lw=0.8, shrinkA=2, shrinkB=4),
            )

        ax.set_yscale("log"); ax.set_xlim(-0.62, 3.05); ax.set_ylim(*ylim)
        ax.set_xticks(xs)
        ax.set_xticklabels(["LOSAT\nv0.1.0", "NCBI BLAST+\n2.17.0", "LOSAT\nv0.2.0"], fontsize=8.3, fontweight="bold", color=theme["ink"])
        ax.set_yticks(yticks); ax.set_yticklabels(ylbls, fontsize=7.5, color=theme["ink2"])
        hide_minor_y(ax)
        for sp in ("top", "right", "bottom"):
            ax.spines[sp].set_visible(False)
        ax.spines["left"].set_color(theme["hairline"])
        ax.tick_params(axis="x", length=0, pad=6)
        ax.set_ylabel(ylab, fontsize=8.0, color=theme["ink2"])
        add_panel(ax, let, title, theme, y=1.035)

    add_footer(
        fig,
        f"LOSAT v0.1.0 release binary vs LOSAT {short_sha(meta.get('losat_commit', ''))} (v0.2.0) and NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')}, {short_hw(meta)}, 16 threads.\n"
        f"*v0.1.0 outputs for megablast and tblastx differed from NCBI BLAST+; v0.2.0 matches NCBI BLAST+ {meta.get('ncbi_version', '2.17.0')} form for form (Table 2 notes).",
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
