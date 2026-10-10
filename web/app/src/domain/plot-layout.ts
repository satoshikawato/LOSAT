// What the dot plot's canvas (ui/DotPlot.vue) and its SVG file (domain/plot-svg.ts) share, so that
// the two draw the same plot: the colours, the bands round the frame, which tick labels are written,
// how a long axis title is cut, and what a very short HSP is drawn as. Nothing here draws; text is
// measured through a function the caller gives (the canvas's, or an estimate in tests).
import type { Orientation } from './result-index';
import { tickLabel, type AxisTicks } from './plot-scale';

/** blast2dotplot.py's colours: the same direction, opposite directions; grey for a BLASTN HSP of one letter. */
export const PLOT_COLORS: Readonly<Record<Orientation, string>> = { forward: '#1f77b4', reverse: '#ff7f0e', unknown: '#7f7f7f' };
/** The orientations in the order of the batches (orientation × 4 + opacity class). */
export const ORIENTATIONS: readonly Orientation[] = ['forward', 'reverse', 'unknown'];
export const GRID = '#d3d3d3';
export const INK = '#000000';
export const TEXT = '#1d2330';
export const FONT_FAMILY = 'system-ui, sans-serif';
/** Tick labels: at least 12 px (S13 screen review L4). */
export const FONT_PX = 12;
export const TITLE_PX_SIZE = 13;
export const FONT = `${FONT_PX}px ${FONT_FAMILY}`;
export const TITLE_FONT = `600 ${TITLE_PX_SIZE}px ${FONT_FAMILY}`;
/** The width of an HSP's line. */
export const LINE_WIDTH = 2;
/** How far a line's stroke reaches beyond the line: half the 3 px dash of a short HSP (the 2 px line reaches 1 px). */
export const STROKE_REACH_PX = 1.5;

// Each axis's band, from the outside in: its title, the tick labels, the ticks.
export const PAD = 4;
export const TITLE_PX = 16;
export const LABEL_PX = 14;
export const GAP = 4;
export const MAJOR_PX = 8;
export const MINOR_PX = 4;
export const AXIS_PX = PAD + TITLE_PX + GAP + LABEL_PX + GAP + MAJOR_PX;
export const MARGIN = { left: AXIS_PX, top: AXIS_PX, right: 14, bottom: 10 } as const;

/** Measures a text (its width in CSS pixels) in a font given as a CSS `font` value. */
export type MeasureText = (text: string, font: string) => number;

/** What a protein axis's title adds against a nucleotide axis while the axes are to scale (S13b decision 26). */
export const scaleNote = (weight: number, toScale: boolean): string => (weight !== 1 && toScale ? `; drawn at ${weight} nt per aa` : '');

/** A line, or a short dash where the HSP is shorter than a pixel or two (a dot with round caps). */
export function lineOrDash(ax: number, ay: number, bx: number, by: number): [number, number, number, number] {
  if (Math.abs(bx - ax) + Math.abs(by - ay) < 1.5) {
    const [mx, my] = [(ax + bx) / 2, (ay + by) / 2];
    return [mx - 1.5, my, mx + 1.5, my];
  }
  return [ax, ay, bx, by];
}

/** The major ticks that carry a label: every one, or every second, third… where labels would collide. */
export function labelledTicks(ticks: AxisTicks, pxPerLetter: number, measure: MeasureText): number[] {
  if (ticks.major.length < 2) return [...ticks.major];
  const widest = Math.max(...ticks.major.map((t) => measure(tickLabel(t, ticks.unit), FONT)));
  const stride = Math.max(1, Math.ceil((widest + 10) / (ticks.steps.major * pxPerLetter)));
  return ticks.major.filter((t) => Math.round(t / ticks.steps.major) % stride === 0);
}

/** The centre of a text `width` px wide at `at`, moved so that the text stays within 0..`room` (2 px to spare). */
export const within = (at: number, width: number, room: number): number => Math.max(width / 2 + 2, Math.min(room - width / 2 - 2, at));

/**
 * An axis title, "Query <ID> (<unit>)". A protein axis against a nucleotide one says that it is
 * drawn at 3 nt per aa (S13b decision 26), unless the 120 px minimum changed its scale (the note
 * under the plot then says so). Where the room is short, the ID is cut and the unit kept (W4b
 * screen review M3: a phone lost the "(kbp)" of a long subject ID).
 */
export function axisTitleText(role: string, id: string, unit: string, weight: number, toScale: boolean, room: number, measure: MeasureText): string {
  return fitId(`${role} `, id, ` (${unit}${scaleNote(weight, toScale)})`, room, measure);
}

/** `head`, `id` and `tail` in `room` px (title font): the ID cut with an ellipsis where needed, else the whole text cut. */
function fitId(head: string, id: string, tail: string, room: number, measure: MeasureText): string {
  const wide = (text: string): number => measure(text, TITLE_FONT);
  const whole = `${head}${id}${tail}`;
  if (wide(whole) <= room) return whole;
  let [low, high] = [0, id.length];
  while (low < high) {
    const mid = Math.ceil((low + high) / 2);
    if (wide(`${head}${id.slice(0, mid)}…${tail}`) <= room) low = mid;
    else high = mid - 1;
  }
  const cut = `${head}${id.slice(0, low)}…${tail}`;
  return wide(cut) <= room ? cut : fit(whole, room, wide);
}

/** A text cut with an ellipsis to `room` px. */
function fit(text: string, room: number, wide: (text: string) => number): string {
  if (wide(text) <= room) return text;
  let [low, high] = [0, text.length];
  while (low < high) {
    const mid = Math.ceil((low + high) / 2);
    if (wide(`${text.slice(0, mid)}…`) <= room) low = mid;
    else high = mid - 1;
  }
  return `${text.slice(0, low)}…`;
}
