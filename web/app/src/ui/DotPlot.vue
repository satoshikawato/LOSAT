<script setup lang="ts">
// The dot plot of the selected query and subject (plan §5.7, design §11.1, REQ-12), drawn after
// the Owner's blast2dotplot.py (S13b, docs/web/ncbi_ui_mapping.md "Dot Plot"): the query along the
// top, the subject down the left side, the same scale on both axes, the script's ticks, grid and
// colours. Every HSP of the pair that the view filters show is a line from its start to its end
// with the coordinates of its record (the outfmt 6 coordinates); its colour says whether the two
// sequences run in the same direction and its opacity the class of its outfmt 6 `pident`. It
// computes nothing else: no comparison of the two sequences runs behind it.
//
// Zoom (Ctrl or ⌘ and the wheel, buttons, + and -), pan (drag, arrow keys) and selection (click,
// n and p) act on the same selection as the tables. The wheel alone and a finger moving up or down
// scroll the page (S13 screen review L4). The grid and the HSPs are drawn on a base canvas in one
// path per colour and opacity class; hover and selection on a canvas above it, so that they redraw
// cheaply. Input is gathered and drawn once per animation frame.
import { computed, nextTick, onMounted, onUnmounted, ref, shallowRef, useId, watch } from 'vue';
import type { HspEntry, HspId, ResultsBrowser, ResultsState } from '../application/results';
import {
  frameRange,
  fromPixelX,
  fromPixelY,
  fullView,
  panView,
  placeBox,
  plotSize,
  segmentDistance,
  toPixelX,
  toPixelY,
  zoomView,
  type Box,
  type View,
} from '../domain/plot-geometry';
import { axisTicks, identityClass, IDENTITY_CLASSES, tickLabel, type AxisTicks } from '../domain/plot-scale';
import { formatCount } from './format';
import './plots.css';

const props = defineProps<{ results: ResultsBrowser; state: ResultsState }>();
const emit = defineEmits<{ 'show-alignment': [id: HspId] }>();

/** blast2dotplot.py's colours: the same direction, opposite directions; grey for a BLASTN HSP of one letter. */
const COLORS = { forward: '#1f77b4', reverse: '#ff7f0e', unknown: '#7f7f7f' } as const;
const ORIENTATIONS = ['forward', 'reverse', 'unknown'] as const;
const ORIENTATION_TEXT = { forward: 'Forward', reverse: 'Reverse', unknown: 'Not in the record' } as const;
const GRID = '#d3d3d3';
const INK = '#000000';
const TEXT = '#1d2330';
/** The halo of the selected HSP (its legend entry is .plot-swatch-selected). */
const HALO = 'rgba(255, 196, 0, 0.75)';
/** Tick labels: at least 12 px (S13 screen review L4). */
const FONT = '12px system-ui, sans-serif';
const TITLE_FONT = '600 13px system-ui, sans-serif';
const PICK_PX = 8;
// Each axis's band, from the outside in: its title, the tick labels, the ticks.
const PAD = 4;
const TITLE_PX = 16;
const LABEL_PX = 14;
const GAP = 4;
const MAJOR_PX = 8;
const MINOR_PX = 4;
const AXIS_PX = PAD + TITLE_PX + GAP + LABEL_PX + GAP + MAJOR_PX;
const MARGIN = { left: AXIS_PX, top: AXIS_PX, right: 14, bottom: 10 } as const;

const stage = ref<HTMLElement>();
const base = ref<HTMLCanvasElement>();
const overlay = ref<HTMLCanvasElement>();
const popup = ref<HTMLElement>();
const popupTitle = useId();
const stageWidth = ref(600);

const loaded = computed(() => props.state.loaded);
const hsps = computed(() => props.state.hsps);
const query = computed(() => props.state.queries.find((q) => q.qIdx === props.state.qIdx));
const subject = computed(() => props.state.subjects.find((s) => s.sIdx === props.state.sIdx));
const queryLength = computed(() => loaded.value?.run.snapshot.query.records[props.state.qIdx ?? -1]?.length ?? 1);
const subjectLength = computed(() => subject.value?.length ?? 1);
const queryId = computed(() => query.value?.id ?? '');
const subjectId = computed(() => subject.value?.first.sseqid ?? '');
const units = computed(() => loaded.value?.units ?? { query: '', subject: '' });
const extent = computed(() => ({ x: Math.max(1, queryLength.value), y: Math.max(1, subjectLength.value) }));

const size = computed(() => plotSize(extent.value, stageWidth.value - MARGIN.left - MARGIN.right));
const box = computed<Box>(() => ({ left: MARGIN.left, top: MARGIN.top, width: size.value.width, height: size.value.height }));
const canvasWidth = computed(() => MARGIN.left + size.value.width + MARGIN.right);
const canvasHeight = computed(() => MARGIN.top + size.value.height + MARGIN.bottom);

/** The HSPs as columns of numbers, and their batches: one per colour and opacity class (orientation × 4 + class). */
interface Segments {
  readonly list: readonly HspEntry[];
  readonly x0: Float64Array;
  readonly y0: Float64Array;
  readonly x1: Float64Array;
  readonly y1: Float64Array;
  readonly batch: Uint8Array;
  readonly batches: readonly Int32Array[];
}
const segments = computed<Segments>(() => {
  const table = loaded.value?.index.table;
  const list = table === undefined ? [] : hsps.value;
  const n = list.length;
  const [x0, y0, x1, y1] = [new Float64Array(n), new Float64Array(n), new Float64Array(n), new Float64Array(n)];
  const batch = new Uint8Array(n);
  const counts = new Int32Array(ORIENTATIONS.length * 4);
  for (let i = 0; i < n; i++) {
    const hsp = list[i]!;
    x0[i] = table!.qStart[hsp.row]!;
    x1[i] = table!.qEnd[hsp.row]!;
    y0[i] = table!.sStart[hsp.row]!;
    y1[i] = table!.sEnd[hsp.row]!;
    const k = ORIENTATIONS.indexOf(hsp.orientation) * 4 + identityClass(hsp.fields.pident);
    batch[i] = k;
    counts[k]!++;
  }
  const batches = Array.from(counts, (count) => new Int32Array(count));
  const filled = new Int32Array(counts.length);
  for (let i = 0; i < n; i++) batches[batch[i]!]![filled[batch[i]!]!++] = i;
  return { list, x0, y0, x1, y1, batch, batches };
});
const selectedIndex = computed(() => {
  const id = props.state.hsp;
  return id === undefined ? -1 : segments.value.list.findIndex((hsp) => hsp.id.qIdx === id.qIdx && hsp.id.rank === id.rank);
});
const selected = computed(() => (selectedIndex.value < 0 ? undefined : segments.value.list[selectedIndex.value]));
const hasUnknown = computed(() => hsps.value.some((hsp) => hsp.orientation === 'unknown'));

// --- the view and the frame loop ---------------------------------------------------------------

const view = shallowRef<View>(fullView(extent.value));
/** The view that input asked for since the last frame. */
let pendingView: View | undefined;
const current = (): View => pendingView ?? view.value;
let hovered = -1;
let pointer: { x: number; y: number } | undefined;
let drag: { x: number; y: number; moved: boolean; id: number } | undefined;
let frame = 0;
let baseDirty = true;
let hoverDirty = false;

/** Asks for a frame: the overlay and the popup are always redrawn, the base when `withBase`. */
function schedule(withBase: boolean): void {
  if (withBase) baseDirty = true;
  if (frame === 0) frame = requestAnimationFrame(render);
}

function setView(next: View): void {
  pendingView = next;
  schedule(true);
}

function render(): void {
  frame = 0;
  if (pendingView !== undefined) {
    const next = pendingView;
    pendingView = undefined;
    const v = view.value;
    if (next.x0 !== v.x0 || next.x1 !== v.x1 || next.y0 !== v.y0 || next.y1 !== v.y1) {
      view.value = next;
      baseDirty = true;
      hoverDirty = pointer !== undefined;
    }
  }
  if (hoverDirty) {
    hoverDirty = false;
    hovered = pointer === undefined || drag?.moved ? -1 : pick(pointer.x, pointer.y);
    if (base.value !== undefined) base.value.style.cursor = drag?.moved ? 'grabbing' : hovered >= 0 ? 'pointer' : '';
  }
  if (baseDirty) {
    baseDirty = false;
    drawBase();
  }
  drawOverlay();
  placePopup();
}

/** Sizes a canvas for the screen's pixel ratio (only when its size changed: that clears it). */
function prepare(canvas: HTMLCanvasElement | undefined): CanvasRenderingContext2D | null {
  if (canvas === undefined) return null;
  const ratio = window.devicePixelRatio || 1;
  const w = Math.round(canvasWidth.value * ratio);
  const h = Math.round(canvasHeight.value * ratio);
  if (canvas.width !== w) canvas.width = w;
  if (canvas.height !== h) canvas.height = h;
  const context = canvas.getContext('2d');
  context?.setTransform(ratio, 0, 0, ratio, 0, 0);
  return context;
}

// --- drawing -------------------------------------------------------------------------------------

function drawBase(): void {
  const context = prepare(base.value);
  if (context === null) return;
  const [W, H] = [canvasWidth.value, canvasHeight.value];
  const v = view.value;
  const b = box.value;
  context.clearRect(0, 0, W, H);
  context.fillStyle = '#ffffff';
  context.fillRect(0, 0, W, H);
  const xTicks = axisTicks(v.x0, v.x1, units.value.query);
  const yTicks = axisTicks(v.y0, v.y1, units.value.subject);
  const xs = (t: number) => Math.round(toPixelX(t, v, b)) + 0.5;
  const ys = (t: number) => Math.round(toPixelY(t, v, b)) + 0.5;

  // The grid at the major and minor ticks.
  context.lineWidth = 1;
  context.strokeStyle = GRID;
  context.beginPath();
  for (const t of [...xTicks.major, ...xTicks.minor]) {
    context.moveTo(xs(t), b.top);
    context.lineTo(xs(t), b.top + b.height);
  }
  for (const t of [...yTicks.major, ...yTicks.minor]) {
    context.moveTo(b.left, ys(t));
    context.lineTo(b.left + b.width, ys(t));
  }
  context.stroke();

  // The HSPs: one path per batch, lines outside the view left out.
  context.save();
  context.beginPath();
  context.rect(b.left, b.top, b.width, b.height);
  context.clip();
  context.lineWidth = 2;
  context.lineCap = 'round';
  const s = segments.value;
  const sx = b.width / (v.x1 - v.x0);
  const sy = b.height / (v.y1 - v.y0);
  s.batches.forEach((members, k) => {
    if (members.length === 0) return;
    context.strokeStyle = COLORS[ORIENTATIONS[k >> 2]!];
    context.globalAlpha = IDENTITY_CLASSES[k & 3]!.opacity;
    context.beginPath();
    for (const i of members) {
      const [ax, bx, ay, by] = [s.x0[i]!, s.x1[i]!, s.y0[i]!, s.y1[i]!];
      if ((ax > bx ? ax : bx) < v.x0 || (ax < bx ? ax : bx) > v.x1 || (ay > by ? ay : by) < v.y0 || (ay < by ? ay : by) > v.y1) continue;
      segmentPath(context, b.left + (ax - v.x0) * sx, b.top + (ay - v.y0) * sy, b.left + (bx - v.x0) * sx, b.top + (by - v.y0) * sy);
    }
    context.stroke();
  });
  context.restore();

  // The frame, the ticks outside it, their labels and the axis titles.
  context.strokeStyle = INK;
  context.lineWidth = 1;
  context.strokeRect(b.left - 0.5, b.top - 0.5, b.width + 1, b.height + 1);
  context.beginPath();
  for (const [ticks, length] of [
    [xTicks.major, MAJOR_PX],
    [xTicks.minor, MINOR_PX],
  ] as const) {
    for (const t of ticks) {
      context.moveTo(xs(t), b.top - 1);
      context.lineTo(xs(t), b.top - 1 - length);
    }
  }
  for (const [ticks, length] of [
    [yTicks.major, MAJOR_PX],
    [yTicks.minor, MINOR_PX],
  ] as const) {
    for (const t of ticks) {
      context.moveTo(b.left - 1, ys(t));
      context.lineTo(b.left - 1 - length, ys(t));
    }
  }
  context.stroke();

  context.fillStyle = TEXT;
  context.font = FONT;
  context.textAlign = 'center';
  context.textBaseline = 'bottom';
  const labelEdge = b.top - 1 - MAJOR_PX - GAP;
  for (const t of labelled(context, xTicks, sx)) {
    const text = tickLabel(t, xTicks.unit);
    // End labels stay inside the canvas (S13 screen review L7); all sit above the ticks, off the frame.
    context.fillText(text, within(toPixelX(t, v, b), context.measureText(text).width, W), labelEdge);
  }
  for (const t of labelled(context, yTicks, sy)) {
    const text = tickLabel(t, yTicks.unit);
    context.save();
    context.translate(labelEdge, within(toPixelY(t, v, b), context.measureText(text).width, H));
    context.rotate(-Math.PI / 2);
    context.fillText(text, 0, 0);
    context.restore();
  }

  context.font = TITLE_FONT;
  context.textBaseline = 'top';
  const xTitle = fit(context, `Query ${queryId.value} (${xTicks.unit.name})`, W - 2 * PAD);
  context.fillText(xTitle, within(b.left + b.width / 2, context.measureText(xTitle).width, W), PAD);
  const yTitle = fit(context, `Subject ${subjectId.value} (${yTicks.unit.name})`, H - 2 * PAD);
  context.save();
  context.translate(PAD, within(b.top + b.height / 2, context.measureText(yTitle).width, H));
  context.rotate(-Math.PI / 2);
  context.fillText(yTitle, 0, 0);
  context.restore();
}

/** The hovered HSP drawn thicker, and the selected one wider on its halo. */
function drawOverlay(): void {
  const context = prepare(overlay.value);
  if (context === null) return;
  context.clearRect(0, 0, canvasWidth.value, canvasHeight.value);
  const b = box.value;
  context.save();
  context.beginPath();
  context.rect(b.left, b.top, b.width, b.height);
  context.clip();
  context.lineCap = 'round';
  const chosen = selectedIndex.value;
  if (hovered >= 0 && hovered !== chosen && hovered < segments.value.list.length) stroke(context, hovered, colorOf(hovered), 4);
  if (chosen >= 0) {
    stroke(context, chosen, HALO, 9);
    stroke(context, chosen, colorOf(chosen), 3.5);
  }
  context.restore();
}

const colorOf = (i: number): string => COLORS[ORIENTATIONS[segments.value.batch[i]! >> 2]!];

function stroke(context: CanvasRenderingContext2D, i: number, color: string, width: number): void {
  const [ax, ay, bx, by] = ends(i);
  context.strokeStyle = color;
  context.lineWidth = width;
  context.beginPath();
  segmentPath(context, ax, ay, bx, by);
  context.stroke();
}

/** A line, or a short dash where the HSP is shorter than a pixel or two (a dot with the round caps). */
function segmentPath(context: CanvasRenderingContext2D, ax: number, ay: number, bx: number, by: number): void {
  if (Math.abs(bx - ax) + Math.abs(by - ay) < 1.5) {
    const [mx, my] = [(ax + bx) / 2, (ay + by) / 2];
    context.moveTo(mx - 1.5, my);
    context.lineTo(mx + 1.5, my);
  } else {
    context.moveTo(ax, ay);
    context.lineTo(bx, by);
  }
}

/** The major ticks that carry a label: every one, or every second, third… where labels would collide. */
function labelled(context: CanvasRenderingContext2D, ticks: AxisTicks, pxPerLetter: number): number[] {
  if (ticks.major.length < 2) return [...ticks.major];
  const widest = Math.max(...ticks.major.map((t) => context.measureText(tickLabel(t, ticks.unit)).width));
  const stride = Math.max(1, Math.ceil((widest + 10) / (ticks.steps.major * pxPerLetter)));
  return ticks.major.filter((t) => Math.round(t / ticks.steps.major) % stride === 0);
}

/** The centre of a text `width` px wide at `at`, moved so that the text stays within 0..`room` (2 px to spare). */
const within = (at: number, width: number, room: number): number => Math.max(width / 2 + 2, Math.min(room - width / 2 - 2, at));

/** A text cut with an ellipsis to `room` px. */
function fit(context: CanvasRenderingContext2D, text: string, room: number): string {
  if (context.measureText(text).width <= room) return text;
  let [low, high] = [0, text.length];
  while (low < high) {
    const mid = Math.ceil((low + high) / 2);
    if (context.measureText(`${text.slice(0, mid)}…`).width <= room) low = mid;
    else high = mid - 1;
  }
  return `${text.slice(0, low)}…`;
}

/** The ends of an HSP's line in CSS pixels of the canvas, in a view. */
function ends(i: number, v = view.value): [number, number, number, number] {
  const s = segments.value;
  const b = box.value;
  return [toPixelX(s.x0[i]!, v, b), toPixelY(s.y0[i]!, v, b), toPixelX(s.x1[i]!, v, b), toPixelY(s.y1[i]!, v, b)];
}

/** The HSP whose line is nearest to a point of the canvas, within PICK_PX, or -1. */
function pick(x: number, y: number): number {
  const b = box.value;
  if (x < b.left - PICK_PX || x > b.left + b.width + PICK_PX || y < b.top - PICK_PX || y > b.top + b.height + PICK_PX) return -1;
  let best = -1;
  let bestDistance = PICK_PX;
  for (let i = 0; i < segments.value.list.length; i++) {
    const [ax, ay, bx, by] = ends(i);
    if (Math.max(ax, bx) < x - PICK_PX || Math.min(ax, bx) > x + PICK_PX || Math.max(ay, by) < y - PICK_PX || Math.min(ay, by) > y + PICK_PX) continue;
    const d = segmentDistance(x, y, ax, ay, bx, by);
    if (d <= bestDistance) {
      best = i;
      bestDistance = d;
    }
  }
  return best;
}

// --- zoom, pan and selection ---------------------------------------------------------------------

const centre = (v: View) => ({ x: (v.x0 + v.x1) / 2, y: (v.y0 + v.y1) / 2 });

function zoom(factor: number, at = centre(current())): void {
  setView(zoomView(current(), factor, at, extent.value));
}

function pan(dxPx: number, dyPx: number): void {
  const v = current();
  const b = box.value;
  setView(panView(v, (-dxPx / b.width) * (v.x1 - v.x0), (-dyPx / b.height) * (v.y1 - v.y0), extent.value));
}

function zoomToSelected(): void {
  const i = selectedIndex.value;
  if (i < 0) return;
  const s = segments.value;
  setView(frameRange({ x0: s.x0[i]!, x1: s.x1[i]!, y0: s.y0[i]!, y1: s.y1[i]! }, extent.value));
}

/** The position (query, subject) under a point of the canvas, kept inside the plot. */
function positionAt(x: number, y: number): { x: number; y: number } {
  const v = current();
  const b = box.value;
  const cx = Math.max(b.left, Math.min(b.left + b.width, x));
  const cy = Math.max(b.top, Math.min(b.top + b.height, y));
  return { x: fromPixelX(cx, v, b), y: fromPixelY(cy, v, b) };
}

/** The wheel zooms only with Ctrl or ⌘ held (a trackpad's pinch also holds Ctrl); otherwise the page scrolls. */
function onWheel(event: WheelEvent): void {
  if (!event.ctrlKey && !event.metaKey) return;
  event.preventDefault();
  const { x, y } = canvasPoint(event);
  zoom(event.deltaY < 0 ? 1.25 : 0.8, positionAt(x, y));
}

function canvasPoint(event: MouseEvent): { x: number; y: number } {
  const rect = base.value!.getBoundingClientRect();
  return { x: event.clientX - rect.left, y: event.clientY - rect.top };
}
function onPointerDown(event: PointerEvent): void {
  drag = { x: event.clientX, y: event.clientY, moved: false, id: event.pointerId };
  base.value?.setPointerCapture(event.pointerId);
}
function onPointerMove(event: PointerEvent): void {
  pointer = canvasPoint(event);
  hoverDirty = true;
  if (drag !== undefined && drag.id === event.pointerId) {
    const dx = event.clientX - drag.x;
    const dy = event.clientY - drag.y;
    if (drag.moved || Math.hypot(dx, dy) >= 3) {
      drag.moved = true;
      pan(dx, dy);
      drag.x = event.clientX;
      drag.y = event.clientY;
      return;
    }
  }
  schedule(false);
}
/** A press that did not move clicks: on a line it selects the HSP and opens its popup, elsewhere it closes the popup. */
function onPointerUp(event: PointerEvent): void {
  if (drag === undefined || drag.id !== event.pointerId) return;
  const moved = drag.moved;
  drag = undefined;
  hoverDirty = true;
  schedule(false);
  if (moved) return;
  const { x, y } = canvasPoint(event);
  const i = pick(x, y);
  if (i < 0) {
    popupOpen.value = false;
    return;
  }
  props.results.selectHsp(segments.value.list[i]!.id);
  popupOpen.value = true;
}
function onPointerLeave(): void {
  pointer = undefined;
  hoverDirty = true;
  schedule(false);
}

function step(offset: number): void {
  const list = segments.value.list;
  if (list.length === 0) return;
  const next = list[(Math.max(0, selectedIndex.value) + offset + list.length) % list.length]!;
  props.results.selectHsp(next.id);
}

function onKey(event: KeyboardEvent): void {
  if (event.ctrlKey || event.metaKey || event.altKey) return;
  const actions: Record<string, () => void> = {
    '+': () => zoom(1.25),
    '=': () => zoom(1.25),
    '-': () => zoom(0.8),
    '0': () => setView(fullView(extent.value)),
    ArrowLeft: () => pan(40, 0),
    ArrowRight: () => pan(-40, 0),
    ArrowUp: () => pan(0, 40),
    ArrowDown: () => pan(0, -40),
    n: () => step(1),
    p: () => step(-1),
    Enter: () => openPopup(true),
  };
  if (event.key === 'Escape' && popupOpen.value) actions.Escape = () => closePopup(false);
  const action = actions[event.key];
  if (action === undefined) return;
  event.preventDefault();
  action();
}

// --- the HSP popup -------------------------------------------------------------------------------

const popupOpen = ref(false);
const popupAt = ref({ left: 0, top: 0 });

/** Opens the popup of the selected HSP; opened from the keyboard, the focus moves into it. */
function openPopup(byKeyboard: boolean): void {
  if (selected.value === undefined) return;
  popupOpen.value = true;
  schedule(false);
  if (byKeyboard) void nextTick(() => popup.value?.focus());
}

/** Closes the popup; closed from inside it, the focus goes back to the plot. */
function closePopup(refocus: boolean): void {
  popupOpen.value = false;
  if (refocus) base.value?.focus({ preventScroll: true });
}

function onPopupKey(event: KeyboardEvent): void {
  if (event.ctrlKey || event.metaKey || event.altKey) return;
  if (event.key === 'Escape') {
    event.preventDefault();
    closePopup(true);
  } else if (event.key === 'n' || event.key === 'p') {
    event.preventDefault();
    step(event.key === 'n' ? 1 : -1);
  }
}

/** Hands the HSP to the results screen, which shows its Range in the Alignments; the popup closes. */
function showAlignment(): void {
  const hsp = selected.value;
  if (hsp === undefined) return;
  popupOpen.value = false;
  emit('show-alignment', hsp.id);
}

/** Puts the popup next to the middle of the selected line (kept within the plot), inside the figure. */
function placePopup(): void {
  const element = popup.value;
  const i = selectedIndex.value;
  if (!popupOpen.value || element === undefined || i < 0) return;
  const [ax, ay, bx, by] = ends(i);
  const b = box.value;
  const anchor = {
    x: Math.max(b.left, Math.min(b.left + b.width, (ax + bx) / 2)),
    y: Math.max(b.top, Math.min(b.top + b.height, (ay + by) / 2)),
  };
  const next = placeBox(anchor, { width: element.offsetWidth, height: element.offsetHeight }, { width: canvasWidth.value, height: canvasHeight.value });
  if (next.left !== popupAt.value.left || next.top !== popupAt.value.top) popupAt.value = next;
}

// --- life cycle ----------------------------------------------------------------------------------

let observer: ResizeObserver | undefined;
onMounted(() => {
  observer = new ResizeObserver((entries) => {
    const next = Math.floor(entries[0]?.contentRect.width ?? stageWidth.value);
    if (next > 0 && next !== stageWidth.value) stageWidth.value = next;
  });
  if (stage.value !== undefined) {
    stageWidth.value = Math.max(200, Math.floor(stage.value.getBoundingClientRect().width));
    observer.observe(stage.value);
  }
  base.value?.addEventListener('wheel', onWheel, { passive: false });
  schedule(true);
});
onUnmounted(() => {
  observer?.disconnect();
  base.value?.removeEventListener('wheel', onWheel);
  if (frame !== 0) cancelAnimationFrame(frame);
});

// A new pair starts from the whole of both sequences, without a popup. (A key, not a getter of
// an array: a new array on every change of the state would count as a new pair.)
const pair = computed(() => [props.state.runId, props.state.qIdx, props.state.sIdx, extent.value.x, extent.value.y].join('|'));
watch(pair, () => {
  view.value = fullView(extent.value);
  pendingView = undefined;
  popupOpen.value = false;
  schedule(true);
});
watch([segments, box, units, queryId, subjectId], () => {
  hovered = -1;
  hoverDirty = pointer !== undefined;
  schedule(true);
});
// The popup follows the selection (n and p move it); with nothing selected it closes.
watch([selectedIndex, popupOpen], () => {
  if (selectedIndex.value < 0) popupOpen.value = false;
  schedule(false);
});

/** Midpoints and ends of the drawn HSPs in CSS pixels of the canvas, for tests (at most 200). */
const targets = computed(() => {
  const v = view.value;
  const list = segments.value.list;
  const out: { hsp: string; x: number; y: number; ends: number[] }[] = [];
  for (let i = 0; i < Math.min(200, list.length); i++) {
    const [ax, ay, bx, by] = ends(i, v);
    out.push({
      hsp: `${list[i]!.id.qIdx}:${list[i]!.id.rank}`,
      x: Math.round((ax + bx) / 2),
      y: Math.round((ay + by) / 2),
      ends: [ax, ay, bx, by].map((n) => Math.round(n * 10) / 10),
    });
  }
  return JSON.stringify(out);
});
const viewText = computed(() => [view.value.x0, view.value.x1, view.value.y0, view.value.y1].map((n) => Math.round(n)).join(','));
const signed = (value: number | undefined) => (value === undefined ? '' : value > 0 ? `+${value}` : String(value));
const framed = (hsp: HspEntry) => hsp.queryFrame !== undefined || hsp.subjectFrame !== undefined;
</script>

<template>
  <figure class="dot-plot" data-testid="dotplot">
    <figcaption class="plot-caption">
      <span class="plot-title" data-testid="dotplot-title">Plot of {{ queryId }} vs {{ subjectId }}</span>
      <span class="muted small"
        >Query {{ queryId }} ({{ formatCount(queryLength) }} {{ units.query }}) against subject {{ subjectId }} ({{ formatCount(subjectLength) }}
        {{ units.subject }}). Each line is an HSP from its start to its end.</span
      >
    </figcaption>
    <div class="plot-tools">
      <button type="button" data-testid="dotplot-zoom-in" @click="zoom(1.5)">Zoom in</button>
      <button type="button" data-testid="dotplot-zoom-out" @click="zoom(1 / 1.5)">Zoom out</button>
      <button type="button" data-testid="dotplot-zoom-hsp" :disabled="!selected" @click="zoomToSelected">Zoom to HSP</button>
      <button type="button" data-testid="dotplot-reset" @click="setView(fullView(extent))">Whole sequences</button>
    </div>
    <div class="plot-help muted small" data-testid="dotplot-help">
      <p>Mouse: hold Ctrl (⌘ on a Mac) and turn the wheel to zoom; drag to move; click a line to show its HSP.</p>
      <p>Touch: use the Zoom buttons; slide a finger sideways to move (up and down scrolls the page); tap a line to show its HSP.</p>
    </div>
    <div ref="stage" class="plot-stage">
      <div class="plot-layers" :style="{ width: `${canvasWidth}px`, height: `${canvasHeight}px` }">
        <canvas
          ref="base"
          class="plot-base"
          tabindex="0"
          role="img"
          :aria-label="`Dot plot of ${segments.list.length} HSPs. Use + and -, or Ctrl or ⌘ and the mouse wheel, to zoom; the arrow keys to move; n and p to select the next or previous HSP; Enter to show the selected HSP.`"
          :style="{ width: `${canvasWidth}px`, height: `${canvasHeight}px` }"
          data-testid="dotplot-canvas"
          :data-segments="segments.list.length"
          :data-selected="selected ? `${selected.id.qIdx}:${selected.id.rank}` : ''"
          :data-view="viewText"
          :data-targets="targets"
          @pointerdown="onPointerDown"
          @pointermove="onPointerMove"
          @pointerup="onPointerUp"
          @pointercancel="drag = undefined"
          @pointerleave="onPointerLeave"
          @keydown="onKey"
        />
        <canvas ref="overlay" class="plot-overlay" aria-hidden="true" :style="{ width: `${canvasWidth}px`, height: `${canvasHeight}px` }" />
        <dialog
          v-if="popupOpen && selected"
          ref="popup"
          open
          class="plot-popup"
          tabindex="-1"
          :aria-labelledby="popupTitle"
          data-testid="dotplot-popup"
          :data-hsp="`${selected.id.qIdx}:${selected.id.rank}`"
          :style="{ left: `${popupAt.left}px`, top: `${popupAt.top}px` }"
          @keydown="onPopupKey"
        >
          <div aria-live="polite">
            <h4 :id="popupTitle">HSP {{ selected.id.rank + 1 }}</h4>
            <dl class="plot-popup-fields">
              <dt>Bit score</dt>
              <dd data-field="bitscore">{{ selected.fields.bitscore }}</dd>
              <dt>E value</dt>
              <dd data-field="evalue">{{ selected.fields.evalue }}</dd>
              <dt>Identity (%)</dt>
              <dd data-field="pident">{{ selected.fields.pident }}</dd>
              <dt>Query</dt>
              <dd data-field="query">{{ selected.fields.qstart }}–{{ selected.fields.qend }} {{ units.query }}</dd>
              <dt>Subject</dt>
              <dd data-field="subject">{{ selected.fields.sstart }}–{{ selected.fields.send }} {{ units.subject }}</dd>
              <template v-if="framed(selected)">
                <dt>Frames</dt>
                <dd data-field="frames">{{ signed(selected.queryFrame) || '–' }} / {{ signed(selected.subjectFrame) || '–' }}</dd>
              </template>
              <dt>Orientation</dt>
              <dd data-field="orientation">{{ ORIENTATION_TEXT[selected.orientation] }}</dd>
              <dt>outfmt 0</dt>
              <dd data-field="outfmt0">{{ selected.inOutfmt0 ? 'shown' : 'not shown' }}</dd>
            </dl>
          </div>
          <div class="plot-popup-actions">
            <button type="button" data-testid="dotplot-popup-alignment" @click="showAlignment">Show alignment</button>
            <button type="button" data-testid="dotplot-popup-close" @click="closePopup(true)">Close</button>
          </div>
        </dialog>
      </div>
    </div>
    <p v-if="!size.toScale" class="muted small plot-note" data-testid="dotplot-scale-note">
      Axes not to scale: the shorter sequence is drawn {{ Math.min(size.width, size.height) }} px long so that its HSPs can be seen.
    </p>
    <ul class="plot-legend">
      <li><span class="plot-swatch" :style="{ background: COLORS.forward }" />Forward: both sequences in the same direction</li>
      <li><span class="plot-swatch" :style="{ background: COLORS.reverse }" />Reverse: one sequence on its minus strand (or reverse frame)</li>
      <li v-if="hasUnknown"><span class="plot-swatch" :style="{ background: COLORS.unknown }" />One letter: the HSP record does not say its strand</li>
      <li class="plot-legend-group">
        Identity (%, whole part):
        <span v-for="item in IDENTITY_CLASSES" :key="item.label" class="plot-legend-item"
          ><span class="plot-swatch" :style="{ background: COLORS.forward, opacity: item.opacity }" />{{ item.label }}</span
        >
      </li>
      <li v-if="selected" data-testid="dotplot-legend-selected"><span class="plot-swatch plot-swatch-selected" />Yellow halo: the selected HSP</li>
    </ul>
    <p v-if="selected" class="muted small" data-testid="dotplot-selected">
      Selected: HSP {{ selected.id.rank + 1 }}, query {{ selected.fields.qstart }}–{{ selected.fields.qend }} {{ units.query }}, subject
      {{ selected.fields.sstart }}–{{ selected.fields.send }} {{ units.subject }}<template v-if="framed(selected)"
        >, frames {{ signed(selected.queryFrame) || '–' }} / {{ signed(selected.subjectFrame) || '–' }}</template
      >.
    </p>
  </figure>
</template>
