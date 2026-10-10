// The SVG file of the dot plot (S15, WP-E; W4b decision 24, docs/web/ncbi_ui_mapping.md "Dot Plot"):
// the plot as the screen shows it - the current zoom and pan, the HSPs the view filters show - with
// the axes, ticks, labels, grid, colours, opacity classes and line rules of the canvas
// (domain/plot-layout.ts and plot-scale.ts hold what both draw from), after the Owner's
// blast2dotplot.py, which writes SVG too. The values drawn are the HSP records' coordinates and
// the outfmt 6 `pident` class that the canvas uses; nothing is computed from sequences. The selected
// HSP's halo and the hover are screen aids and are left out.
//
// The file is made as a sequence of text blocks (`dotPlotSvg` is a generator, one `<line>` per
// HSP), so that the writer can save a pair of 5,993 HSPs or more without building one string. It is
// plain SVG: no script, no event attribute, no reference to another file (no `href`, no `url(`;
// the frame clips the lines by a nested `<svg>`, not a clip path), and every text from an input is
// escaped. It is a LOSAT Web application format, not an NCBI graphic.
import type { Box, View, Weights } from './plot-geometry';
import { lineNearView, toPixelX, toPixelY } from './plot-geometry';
import {
  axisTitleText,
  FONT,
  FONT_FAMILY,
  FONT_PX,
  GAP,
  GRID,
  INK,
  labelledTicks,
  lineOrDash,
  LINE_WIDTH,
  MAJOR_PX,
  MINOR_PX,
  ORIENTATIONS,
  PAD,
  PLOT_COLORS,
  STROKE_REACH_PX,
  TEXT,
  TITLE_FONT,
  TITLE_PX_SIZE,
  within,
  type MeasureText,
} from './plot-layout';
import { axisTicks, IDENTITY_CLASSES, tickLabel } from './plot-scale';

export const SVG_MIME = 'image/svg+xml';

/**
 * The HSP lines of the view's pair as columns, in the order of the engine: the ends in the records'
 * coordinates, and each line's batch, `orientation × 4 + opacity class` (the canvas's batches).
 */
export interface PlotLines {
  readonly count: number;
  readonly x0: ArrayLike<number>;
  readonly y0: ArrayLike<number>;
  readonly x1: ArrayLike<number>;
  readonly y1: ArrayLike<number>;
  readonly batch: ArrayLike<number>;
}

export interface DotPlotSvgInput {
  /** The run's number (the "Run N" of the run list). */
  readonly run: number;
  readonly queryId: string;
  readonly subjectId: string;
  /** The records' lengths, for the description. */
  readonly queryLength: number;
  readonly subjectLength: number;
  /** The run's unit of each axis: "nt" or "aa". */
  readonly units: { readonly query: string; readonly subject: string };
  /** The size of the plot in CSS pixels, margins included. */
  readonly width: number;
  readonly height: number;
  /** The frame inside it. */
  readonly box: Box;
  readonly view: View;
  readonly weights: Weights;
  /** False where the shorter side was raised to the minimum (the file says so in its description). */
  readonly toScale: boolean;
  readonly lines: PlotLines;
  readonly measure: MeasureText;
}

/** The file name: the run's number and the records' positions (1-based), never an ID or an input's name. */
export const dotPlotFileName = (run: number, queryPosition: number, subjectPosition: number): string =>
  `losat-run${run}-dotplot-q${queryPosition}-s${subjectPosition}.svg`;

/**
 * Text for a text node or an attribute: `& < > " '` escaped, and every character that XML 1.0
 * cannot hold (controls other than tab, line feed and carriage return, lone surrogates, U+FFFE and
 * U+FFFF) replaced by U+FFFD, so that an ID with such a character still gives a file that parses.
 */
export function xmlText(text: string): string {
  return text
    .replace(/[^\u0009\u000A\u000D -퟿-�\u{10000}-\u{10FFFF}]/gu, '�')
    .replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&apos;' })[c]!);
}

/** A number with at most two decimals (pixels), and no exponent. */
const n = (value: number): string => String(Math.round(value * 100) / 100);

/** Where a text's baseline lies against the edge the canvas aligned it to (fractions of the font size). */
const ASCENT_EM = 0.8;
const DESCENT_EM = 0.2;

/** The SVG text of a dot plot, block by block, in the order of the file. */
export function* dotPlotSvg(plot: DotPlotSvgInput): Generator<string> {
  const { box: b, view: v, width: W, height: H, measure } = plot;
  const sx = b.width / (v.x1 - v.x0);
  const sy = b.height / (v.y1 - v.y0);
  const xTicks = axisTicks(v.x0, v.x1, plot.units.query, sx);
  const yTicks = axisTicks(v.y0, v.y1, plot.units.subject, sy);
  const xs = (t: number): number => Math.round(toPixelX(t, v, b)) + 0.5;
  const ys = (t: number): number => Math.round(toPixelY(t, v, b)) + 0.5;

  const title = `LOSAT Web dot plot of run ${plot.run}: query ${plot.queryId} against subject ${plot.subjectId}`;
  const description = [
    `A dot plot of run ${plot.run} drawn by LOSAT Web: query ${plot.queryId} (${plot.queryLength} ${plot.units.query}) along the top against subject ${plot.subjectId} (${plot.subjectLength} ${plot.units.subject}) down the left side.`,
    `The view shows query positions ${n(v.x0)} to ${n(v.x1)} and subject positions ${n(v.y0)} to ${n(v.y1)}; each line is one of ${plot.lines.count} HSPs from its start to its end, in the coordinates of the HSP record.`,
    'Blue lines run in the same direction on both sequences, orange lines in opposite directions; the lighter a line, the lower its percent identity.',
    plot.toScale ? '' : 'The axes are not to scale: the shorter sequence is drawn longer so that its HSPs can be seen.',
    'This is a LOSAT Web format (an application format), not an NCBI graphic.',
  ]
    .filter((sentence) => sentence !== '')
    .join(' ');

  yield `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="${n(W)}" height="${n(H)}" viewBox="0 0 ${n(W)} ${n(H)}" font-family="${xmlText(FONT_FAMILY)}">\n`;
  yield `<title>${xmlText(title)}</title>\n<desc>${xmlText(description)}</desc>\n`;
  yield `<rect x="0" y="0" width="${n(W)}" height="${n(H)}" fill="#ffffff"/>\n`;

  // The grid at the major and minor ticks.
  const grid: string[] = [];
  for (const t of [...xTicks.major, ...xTicks.minor]) grid.push(`M${n(xs(t))} ${n(b.top)}V${n(b.top + b.height)}`);
  for (const t of [...yTicks.major, ...yTicks.minor]) grid.push(`M${n(b.left)} ${n(ys(t))}H${n(b.left + b.width)}`);
  if (grid.length > 0) yield `<path d="${grid.join('')}" fill="none" stroke="${GRID}" stroke-width="1"/>\n`;

  // The HSPs inside the frame (a nested <svg> clips them to it), batch by batch, one line each.
  // Opaque lines whose ends fall on the same pixels as those of one before them are left out, as the
  // canvas does: they add almost nothing, and a repeat can have thousands a fraction of a pixel apart.
  yield `<svg x="${n(b.left)}" y="${n(b.top)}" width="${n(b.width)}" height="${n(b.height)}" viewBox="${n(b.left)} ${n(b.top)} ${n(b.width)} ${n(b.height)}" overflow="hidden">\n`;
  yield `<g fill="none" stroke-width="${LINE_WIDTH}" stroke-linecap="butt">\n`;
  const { lines } = plot;
  const [padX, padY] = [STROKE_REACH_PX / sx, STROKE_REACH_PX / sy];
  for (let k = 0; k < ORIENTATIONS.length * IDENTITY_CLASSES.length; k++) {
    const opacity = IDENTITY_CLASSES[k & 3]!.opacity;
    const seen = opacity === 1 ? new Set<string>() : undefined;
    let opened = false;
    for (let i = 0; i < lines.count; i++) {
      if (lines.batch[i] !== k) continue;
      const [ax, ay, bx, by] = [lines.x0[i]!, lines.y0[i]!, lines.x1[i]!, lines.y1[i]!];
      if (!lineNearView(ax, ay, bx, by, v, padX, padY)) continue;
      const [px0, py0, px1, py1] = [b.left + (ax - v.x0) * sx, b.top + (ay - v.y0) * sy, b.left + (bx - v.x0) * sx, b.top + (by - v.y0) * sy];
      if (seen !== undefined) {
        const key = `${Math.round(px0)},${Math.round(py0)},${Math.round(px1)},${Math.round(py1)}`;
        if (seen.has(key)) continue;
        seen.add(key);
      }
      if (!opened) {
        opened = true;
        yield `<g stroke="${PLOT_COLORS[ORIENTATIONS[k >> 2]!]}" stroke-opacity="${opacity}">\n`;
      }
      const [x0, y0, x1, y1] = lineOrDash(px0, py0, px1, py1);
      yield `<line x1="${n(x0)}" y1="${n(y0)}" x2="${n(x1)}" y2="${n(y1)}"/>\n`;
    }
    if (opened) yield '</g>\n';
  }
  yield '</g>\n</svg>\n';

  // The frame and the ticks outside it.
  const ticks: string[] = [];
  for (const [list, length] of [
    [xTicks.major, MAJOR_PX],
    [xTicks.minor, MINOR_PX],
  ] as const) {
    for (const t of list) ticks.push(`M${n(xs(t))} ${n(b.top - 1)}V${n(b.top - 1 - length)}`);
  }
  for (const [list, length] of [
    [yTicks.major, MAJOR_PX],
    [yTicks.minor, MINOR_PX],
  ] as const) {
    for (const t of list) ticks.push(`M${n(b.left - 1)} ${n(ys(t))}H${n(b.left - 1 - length)}`);
  }
  yield `<rect x="${n(b.left - 0.5)}" y="${n(b.top - 0.5)}" width="${n(b.width + 1)}" height="${n(b.height + 1)}" fill="none" stroke="${INK}" stroke-width="1"/>\n`;
  if (ticks.length > 0) yield `<path d="${ticks.join('')}" fill="none" stroke="${INK}" stroke-width="1"/>\n`;

  // The tick labels: above the ticks (the query) and left of them, turned (the subject). The canvas
  // aligns them to the bottom edge of the text; the baseline is a font descent above it.
  const labelEdge = b.top - 1 - MAJOR_PX - GAP;
  yield `<g fill="${TEXT}" font-size="${FONT_PX}" text-anchor="middle">\n`;
  for (const t of labelledTicks(xTicks, sx, measure)) {
    const text = tickLabel(t, xTicks.unit);
    const at = within(toPixelX(t, v, b), measure(text, FONT), W);
    yield `<text x="${n(at)}" y="${n(labelEdge - DESCENT_EM * FONT_PX)}">${xmlText(text)}</text>\n`;
  }
  for (const t of labelledTicks(yTicks, sy, measure)) {
    const text = tickLabel(t, yTicks.unit);
    const at = within(toPixelY(t, v, b), measure(text, FONT), H);
    yield `<text transform="translate(${n(labelEdge)} ${n(at)}) rotate(-90)" y="${n(-DESCENT_EM * FONT_PX)}">${xmlText(text)}</text>\n`;
  }
  yield '</g>\n';

  // The axis titles, aligned to the top edge of the text.
  const xTitle = axisTitleText('Query', plot.queryId, xTicks.unit.name, plot.weights.x, plot.toScale, W - 2 * PAD, measure);
  const yTitle = axisTitleText('Subject', plot.subjectId, yTicks.unit.name, plot.weights.y, plot.toScale, H - 2 * PAD, measure);
  const xAt = within(b.left + b.width / 2, measure(xTitle, TITLE_FONT), W);
  const yAt = within(b.top + b.height / 2, measure(yTitle, TITLE_FONT), H);
  yield `<g fill="${TEXT}" font-size="${TITLE_PX_SIZE}" font-weight="600" text-anchor="middle">\n`;
  yield `<text x="${n(xAt)}" y="${n(PAD + ASCENT_EM * TITLE_PX_SIZE)}">${xmlText(xTitle)}</text>\n`;
  yield `<text transform="translate(${n(PAD)} ${n(yAt)}) rotate(-90)" y="${n(ASCENT_EM * TITLE_PX_SIZE)}">${xmlText(yTitle)}</text>\n`;
  yield '</g>\n</svg>\n';
}
