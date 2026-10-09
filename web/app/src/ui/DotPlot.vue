<script setup lang="ts">
// The dot plot of the selected query and subject (plan §5.7, design §11.1, REQ-12): every
// HSP of the pair that the view filters show, drawn from its start to its end with the
// coordinates of its record (the outfmt 6 coordinates), on axes in the units of the
// records (nt or aa). It computes nothing else: no comparison of the two sequences runs
// behind it. Zoom (Ctrl or ⌘ and the wheel, buttons, + and -), pan (drag, arrow keys) and
// selection (click, n and p) act on the same selection as the tables. The wheel alone and a
// finger moving up or down scroll the page (S13 screen review L4).
import { computed, nextTick, onMounted, onUnmounted, ref, watch } from 'vue';
import type { HspEntry, ResultsBrowser, ResultsState } from '../application/results';
import { formatCount } from './format';

const props = defineProps<{ results: ResultsBrowser; state: ResultsState }>();

const MARGIN = { left: 72, right: 14, top: 14, bottom: 46 } as const;
const COLORS = { forward: '#1f5fbf', reverse: '#c2410c', unknown: '#5b6475' } as const;
const PICK_PX = 8;
/** The halo of the selected HSP (its legend entry is .swatch-selected). */
const HALO = 'rgba(255, 196, 0, 0.75)';
/** Tick and axis labels: at least 12 px (S13 screen review L4). */
const FONT = '12px system-ui, sans-serif';

const wrap = ref<HTMLElement>();
const canvas = ref<HTMLCanvasElement>();
const width = ref(600);
const height = computed(() => Math.round(Math.min(480, Math.max(260, width.value * 0.7))));

const loaded = computed(() => props.state.loaded);
const query = computed(() => props.state.queries.find((q) => q.qIdx === props.state.qIdx));
const subject = computed(() => props.state.subjects.find((s) => s.sIdx === props.state.sIdx));
const queryLength = computed(() => loaded.value?.run.snapshot.query.records[props.state.qIdx ?? -1]?.length ?? 1);
const subjectLength = computed(() => subject.value?.length ?? 1);
const units = computed(() => loaded.value?.units ?? { query: '', subject: '' });

interface View {
  readonly x0: number;
  readonly x1: number;
  readonly y0: number;
  readonly y1: number;
}
const full = (): View => ({ x0: 0, x1: queryLength.value, y0: 0, y1: subjectLength.value });
const view = ref<View>(full());

interface Segment {
  readonly hsp: HspEntry;
  readonly x0: number;
  readonly y0: number;
  readonly x1: number;
  readonly y1: number;
}
const segments = computed<Segment[]>(() => {
  const table = loaded.value?.index.table;
  if (table === undefined) return [];
  return props.state.hsps.map((hsp) => ({
    hsp,
    x0: table.qStart[hsp.row]!,
    x1: table.qEnd[hsp.row]!,
    y0: table.sStart[hsp.row]!,
    y1: table.sEnd[hsp.row]!,
  }));
});
const selected = computed(() => segments.value.find((s) => s.hsp.id.qIdx === props.state.hsp?.qIdx && s.hsp.id.rank === props.state.hsp?.rank));
const hasUnknown = computed(() => segments.value.some((s) => s.hsp.orientation === 'unknown'));

const plotW = () => Math.max(10, width.value - MARGIN.left - MARGIN.right);
const plotH = () => Math.max(10, height.value - MARGIN.top - MARGIN.bottom);
const px = (x: number, v = view.value) => MARGIN.left + ((x - v.x0) / (v.x1 - v.x0)) * plotW();
const py = (y: number, v = view.value) => MARGIN.top + plotH() - ((y - v.y0) / (v.y1 - v.y0)) * plotH();

/** Tick positions of a range: 1, 2 or 5 times a power of ten, about five of them. */
function ticks(from: number, to: number): number[] {
  const span = to - from;
  if (!(span > 0)) return [];
  const raw = span / 5;
  const power = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 5, 10].map((m) => m * power).find((s) => s >= raw) ?? raw;
  const out: number[] = [];
  for (let t = Math.ceil(from / step) * step; t <= to + 1e-9; t += step) out.push(Math.max(1, Math.round(t)));
  return [...new Set(out)];
}

function draw(): void {
  const element = canvas.value;
  if (element === undefined) return;
  const ratio = window.devicePixelRatio || 1;
  element.width = Math.round(width.value * ratio);
  element.height = Math.round(height.value * ratio);
  const context = element.getContext('2d');
  if (context === null) return;
  context.setTransform(ratio, 0, 0, ratio, 0, 0);
  context.clearRect(0, 0, width.value, height.value);
  context.fillStyle = '#ffffff';
  context.fillRect(0, 0, width.value, height.value);
  const v = view.value;

  // Axes, ticks and labels.
  context.strokeStyle = '#d5dae3';
  context.fillStyle = '#5b6475';
  context.font = FONT;
  context.lineWidth = 1;
  context.textAlign = 'center';
  context.textBaseline = 'top';
  for (const t of ticks(v.x0, v.x1)) {
    const x = px(t);
    context.beginPath();
    context.moveTo(x, MARGIN.top);
    context.lineTo(x, MARGIN.top + plotH());
    context.stroke();
    // A label at an end of the axis stays inside the canvas.
    const label = formatCount(t);
    const half = context.measureText(label).width / 2;
    context.fillText(label, Math.min(Math.max(x, half), width.value - half), MARGIN.top + plotH() + 4);
  }
  context.textAlign = 'right';
  context.textBaseline = 'middle';
  for (const t of ticks(v.y0, v.y1)) {
    const y = py(t);
    context.beginPath();
    context.moveTo(MARGIN.left, y);
    context.lineTo(MARGIN.left + plotW(), y);
    context.stroke();
    context.fillText(formatCount(t), MARGIN.left - 4, y);
  }
  context.strokeStyle = '#1d2330';
  context.strokeRect(MARGIN.left, MARGIN.top, plotW(), plotH());
  context.fillStyle = '#1d2330';
  context.textAlign = 'center';
  context.textBaseline = 'bottom';
  context.fillText(`Query (${units.value.query})`, MARGIN.left + plotW() / 2, height.value - 2);
  context.save();
  context.translate(12, MARGIN.top + plotH() / 2);
  context.rotate(-Math.PI / 2);
  context.textBaseline = 'top';
  context.fillText(`Subject (${units.value.subject})`, 0, 0);
  context.restore();

  // HSPs, the selected one last and wider.
  context.save();
  context.beginPath();
  context.rect(MARGIN.left, MARGIN.top, plotW(), plotH());
  context.clip();
  context.lineCap = 'round';
  const chosen = selected.value;
  for (const segment of [...segments.value.filter((s) => s !== chosen), ...(chosen ? [chosen] : [])]) {
    const isSelected = segment === chosen;
    const color = COLORS[segment.hsp.orientation];
    const [ax, ay, bx, by] = [px(segment.x0), py(segment.y0), px(segment.x1), py(segment.y1)];
    if (isSelected) {
      context.strokeStyle = HALO;
      context.lineWidth = 9;
      line(context, ax, ay, bx, by);
    }
    context.strokeStyle = color;
    context.fillStyle = color;
    context.lineWidth = isSelected ? 3.5 : 2;
    if (Math.hypot(bx - ax, by - ay) < 2) {
      context.beginPath();
      context.arc((ax + bx) / 2, (ay + by) / 2, isSelected ? 4 : 3, 0, 2 * Math.PI);
      context.fill();
    } else {
      line(context, ax, ay, bx, by);
    }
  }
  context.restore();
}

function line(context: CanvasRenderingContext2D, ax: number, ay: number, bx: number, by: number): void {
  context.beginPath();
  context.moveTo(ax, ay);
  context.lineTo(bx, by);
  context.stroke();
}

// --- zoom and pan -------------------------------------------------------------------------

function clamp(next: View): View {
  const qMax = queryLength.value;
  const sMax = subjectLength.value;
  const fit = (a: number, b: number, max: number): [number, number] => {
    const span = Math.min(max, Math.max(b - a, Math.min(max, 4)));
    let start = a + (b - a - span) / 2;
    start = Math.min(Math.max(0, start), max - span);
    return [start, start + span];
  };
  const [x0, x1] = fit(next.x0, next.x1, qMax);
  const [y0, y1] = fit(next.y0, next.y1, sMax);
  return { x0, x1, y0, y1 };
}

function zoom(factor: number, cx = MARGIN.left + plotW() / 2, cy = MARGIN.top + plotH() / 2): void {
  const v = view.value;
  const fx = v.x0 + ((cx - MARGIN.left) / plotW()) * (v.x1 - v.x0);
  const fy = v.y1 - ((cy - MARGIN.top) / plotH()) * (v.y1 - v.y0);
  view.value = clamp({
    x0: fx - (fx - v.x0) / factor,
    x1: fx + (v.x1 - fx) / factor,
    y0: fy - (fy - v.y0) / factor,
    y1: fy + (v.y1 - fy) / factor,
  });
}

function pan(dxPx: number, dyPx: number): void {
  const v = view.value;
  const dx = (dxPx / plotW()) * (v.x1 - v.x0);
  const dy = (dyPx / plotH()) * (v.y1 - v.y0);
  view.value = clamp({ x0: v.x0 - dx, x1: v.x1 - dx, y0: v.y0 + dy, y1: v.y1 + dy });
}

function zoomToSelected(): void {
  const s = selected.value;
  if (s === undefined) return;
  const padX = Math.max(10, Math.abs(s.x1 - s.x0) * 0.25);
  const padY = Math.max(10, Math.abs(s.y1 - s.y0) * 0.25);
  view.value = clamp({
    x0: Math.min(s.x0, s.x1) - padX,
    x1: Math.max(s.x0, s.x1) + padX,
    y0: Math.min(s.y0, s.y1) - padY,
    y1: Math.max(s.y0, s.y1) + padY,
  });
}

/** The wheel zooms only with Ctrl or ⌘ held (a trackpad's pinch also holds Ctrl); otherwise the page scrolls. */
function onWheel(event: WheelEvent): void {
  if (!event.ctrlKey && !event.metaKey) return;
  event.preventDefault();
  const rect = canvas.value!.getBoundingClientRect();
  zoom(event.deltaY < 0 ? 1.25 : 0.8, event.clientX - rect.left, event.clientY - rect.top);
}

let drag: { x: number; y: number; moved: boolean; id: number } | undefined;
function onPointerDown(event: PointerEvent): void {
  drag = { x: event.clientX, y: event.clientY, moved: false, id: event.pointerId };
  canvas.value?.setPointerCapture(event.pointerId);
}
function onPointerMove(event: PointerEvent): void {
  if (drag === undefined || drag.id !== event.pointerId) return;
  const dx = event.clientX - drag.x;
  const dy = event.clientY - drag.y;
  if (!drag.moved && Math.hypot(dx, dy) < 3) return;
  drag.moved = true;
  pan(dx, dy);
  drag.x = event.clientX;
  drag.y = event.clientY;
}
function onPointerUp(event: PointerEvent): void {
  if (drag === undefined || drag.id !== event.pointerId) return;
  const moved = drag.moved;
  drag = undefined;
  if (moved) return;
  const rect = canvas.value!.getBoundingClientRect();
  pick(event.clientX - rect.left, event.clientY - rect.top);
}

/** Selects the HSP nearest to a point of the canvas, if one is within a few pixels. */
function pick(x: number, y: number): void {
  let best: Segment | undefined;
  let bestDistance = PICK_PX;
  for (const segment of segments.value) {
    const d = distance(x, y, px(segment.x0), py(segment.y0), px(segment.x1), py(segment.y1));
    if (d <= bestDistance) {
      best = segment;
      bestDistance = d;
    }
  }
  if (best !== undefined) props.results.selectHsp(best.hsp.id);
}

function distance(x: number, y: number, ax: number, ay: number, bx: number, by: number): number {
  const lx = bx - ax;
  const ly = by - ay;
  const length = lx * lx + ly * ly;
  const t = length === 0 ? 0 : Math.max(0, Math.min(1, ((x - ax) * lx + (y - ay) * ly) / length));
  return Math.hypot(x - (ax + t * lx), y - (ay + t * ly));
}

function step(offset: number): void {
  const list = props.state.hsps;
  if (list.length === 0) return;
  const at = list.findIndex((hsp) => hsp.id.qIdx === props.state.hsp?.qIdx && hsp.id.rank === props.state.hsp?.rank);
  const next = list[(Math.max(0, at) + offset + list.length) % list.length]!;
  props.results.selectHsp(next.id);
}

function onKey(event: KeyboardEvent): void {
  const actions: Record<string, () => void> = {
    '+': () => zoom(1.25),
    '=': () => zoom(1.25),
    '-': () => zoom(0.8),
    '0': () => (view.value = full()),
    ArrowLeft: () => pan(40, 0),
    ArrowRight: () => pan(-40, 0),
    ArrowUp: () => pan(0, 40),
    ArrowDown: () => pan(0, -40),
    n: () => step(1),
    p: () => step(-1),
  };
  const action = actions[event.key];
  if (action === undefined) return;
  event.preventDefault();
  action();
}

// --- life cycle --------------------------------------------------------------------------

let observer: ResizeObserver | undefined;
onMounted(() => {
  observer = new ResizeObserver((entries) => {
    const next = Math.floor(entries[0]?.contentRect.width ?? width.value);
    if (next > 0 && next !== width.value) width.value = next;
  });
  if (wrap.value !== undefined) {
    width.value = Math.max(200, Math.floor(wrap.value.getBoundingClientRect().width));
    observer.observe(wrap.value);
  }
  canvas.value?.addEventListener('wheel', onWheel, { passive: false });
  void nextTick(draw);
});
onUnmounted(() => {
  observer?.disconnect();
  canvas.value?.removeEventListener('wheel', onWheel);
});

// A new pair starts from the whole of both sequences.
watch(
  () => [props.state.runId, props.state.qIdx, props.state.sIdx, queryLength.value, subjectLength.value] as const,
  () => (view.value = full()),
);
watch([segments, selected, view, width, height], () => draw(), { flush: 'post' });

/** Midpoints of the drawn HSPs in CSS pixels of the canvas, for tests (at most 200). */
const targets = computed(() =>
  JSON.stringify(
    segments.value.slice(0, 200).map((s) => ({
      hsp: `${s.hsp.id.qIdx}:${s.hsp.id.rank}`,
      x: Math.round((px(s.x0) + px(s.x1)) / 2),
      y: Math.round((py(s.y0) + py(s.y1)) / 2),
    })),
  ),
);
const viewText = computed(() => [view.value.x0, view.value.x1, view.value.y0, view.value.y1].map((n) => Math.round(n)).join(','));
const signed = (value: number | undefined) => (value === undefined ? '' : value > 0 ? `+${value}` : String(value));
</script>

<template>
  <figure class="dotplot" data-testid="dotplot">
    <figcaption>
      <span>Query {{ query?.id }} ({{ formatCount(queryLength) }} {{ units.query }}) against subject {{ subject?.first.sseqid }}
        ({{ formatCount(subjectLength) }} {{ units.subject }}). Each line is an HSP from its start to its end.</span>
    </figcaption>
    <div class="dotplot-tools">
      <button type="button" data-testid="dotplot-zoom-in" @click="zoom(1.5)">Zoom in</button>
      <button type="button" data-testid="dotplot-zoom-out" @click="zoom(1 / 1.5)">Zoom out</button>
      <button type="button" data-testid="dotplot-zoom-hsp" :disabled="!selected" @click="zoomToSelected">Zoom to HSP</button>
      <button type="button" data-testid="dotplot-reset" @click="view = full()">Whole sequences</button>
    </div>
    <p class="dotplot-help muted small">Hold Ctrl (⌘ on a Mac) and turn the mouse wheel to zoom; drag to move; click a line to select its HSP.</p>
    <div ref="wrap" class="dotplot-canvas">
      <canvas
        ref="canvas"
        tabindex="0"
        role="img"
        :aria-label="`Dot plot of ${segments.length} HSPs. Use + and -, or Ctrl or ⌘ and the mouse wheel, to zoom; the arrow keys to move; n and p to select the next or previous HSP.`"
        :style="{ width: `${width}px`, height: `${height}px` }"
        data-testid="dotplot-canvas"
        :data-segments="segments.length"
        :data-selected="selected ? `${selected.hsp.id.qIdx}:${selected.hsp.id.rank}` : ''"
        :data-view="viewText"
        :data-targets="targets"
        @pointerdown="onPointerDown"
        @pointermove="onPointerMove"
        @pointerup="onPointerUp"
        @pointercancel="drag = undefined"
        @keydown="onKey"
      />
    </div>
    <ul class="dotplot-legend">
      <li><span class="swatch" :style="{ background: COLORS.forward }" />Forward: both sequences in the same direction</li>
      <li><span class="swatch" :style="{ background: COLORS.reverse }" />Reverse: one sequence on its minus strand (or reverse frame)</li>
      <li v-if="hasUnknown"><span class="swatch" :style="{ background: COLORS.unknown }" />One letter: the HSP record does not say its strand</li>
      <li v-if="selected" data-testid="dotplot-legend-selected"><span class="swatch swatch-selected" />Yellow halo: the selected HSP</li>
    </ul>
    <p v-if="selected" class="muted small" data-testid="dotplot-selected">
      Selected: HSP {{ selected.hsp.id.rank + 1 }}, query {{ selected.hsp.fields.qstart }}–{{ selected.hsp.fields.qend }}
      {{ units.query }}, subject {{ selected.hsp.fields.sstart }}–{{ selected.hsp.fields.send }} {{ units.subject }}<template
        v-if="selected.hsp.queryFrame !== undefined || selected.hsp.subjectFrame !== undefined"
        >, frames {{ signed(selected.hsp.queryFrame) || '–' }} / {{ signed(selected.hsp.subjectFrame) || '–' }}</template
      >.
    </p>
  </figure>
</template>
