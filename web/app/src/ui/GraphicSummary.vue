<script setup lang="ts">
// The Graphic Summary (S13b, docs/web/ncbi_ui_mapping.md "Graphic Summary"), after NCBI's: the
// selected query as a bar with a ruler, and under it one row per subject, in the Descriptions'
// order and view filters. Each HSP that the filters show is a thin bar from the smaller to the
// larger query coordinate of its record, coloured by NCBI's "Alignment Scores" bin of its record's
// bit score; the HSPs of a subject are joined by a thin grey line. The bins only classify the
// engine's number; hover and focus show the outfmt 6 strings. Click or Enter selects the HSP and
// asks the results screen to show its alignment.
//
// One canvas holds the figure and draws only the rows in view, so that thousands of HSPs stay
// smooth; with many subjects the rows scroll inside the figure under the query bar. Since the
// figure heads the one page of the results (2026-10-10, ClassicResults.vue), its band is compact,
// as NCBI's classic overview: rows of 8 px, 30 in view, so that the Descriptions follow on the
// same screen.
import { computed, onMounted, onUnmounted, ref, useId, watch } from 'vue';
import type { HspId, ResultsBrowser, ResultsState, SubjectEntry } from '../application/results';
import { interval } from '../domain/coordinates';
import { headingTitle } from '../domain/outfmt0';
import { placeBox } from '../domain/plot-geometry';
import { rulerTicks, scoreBin, SCORE_BINS } from '../domain/plot-scale';
import { formatCount } from './format';
import { shownWhole, wholeKey } from './shownWhole';
import './plots.css';

const props = defineProps<{ results: ResultsBrowser; state: ResultsState }>();
const emit = defineEmits<{ 'show-alignment': [id: HspId] }>();

/** Subjects drawn until "Show all" (NCBI draws the top hits). */
const FIRST_SUBJECTS = 100;
/** Rows in view before the rows scroll inside the figure. */
const VISIBLE_ROWS = 30;
const ROW_PX = 8;
const BAR_PX = 4;
const SIDE = 12;
const MAX_WIDTH = 1000;
const QUERY_TOP = 4;
const QUERY_PX = 16;
const TICK_PX = 5;
const LABEL_TOP = QUERY_TOP + QUERY_PX + TICK_PX + 2;
const HEADER_PX = LABEL_TOP + 14 + 8;
/** How far beside a bar a pointer still picks it. */
const PICK_PX = 4;
const FONT = '12px system-ui, sans-serif';
const BOLD = '600 12px system-ui, sans-serif';
const ACCENT = '#1f5fbf';
const JOIN = '#8a93a3';
const TEXT = '#1d2330';
const MUTED = '#5b6475';
/** The halo of the selected HSP, as the dot plot draws it. */
const HALO = 'rgba(255, 196, 0, 0.75)';

const scroller = ref<HTMLElement>();
const canvas = ref<HTMLCanvasElement>();
const popover = ref<HTMLElement>();
const legendTitle = useId();
const areaWidth = ref(600);
const scrollTop = ref(0);
const hovered = ref(-1);
/** The HSP that the arrow keys move (an index of the layout), and whether the keyboard focus shows it. */
const cursor = ref(-1);
const keyboard = ref(false);
const popoverAt = ref({ left: 0, top: 0 });

const loaded = computed(() => props.state.loaded);
const subjects = computed(() => props.state.subjects);
const record = computed(() => loaded.value?.run.snapshot.query.records[props.state.qIdx ?? -1]);
const queryLength = computed(() => Math.max(1, record.value?.length ?? 1));
const unit = computed(() => loaded.value?.units.query ?? '');
/** "Show all" of this run's query (kept while the results tab is left, shownWhole.ts). */
const key = computed(() => wholeKey(props.state.runId, props.state.qIdx));
const showAll = computed(() => shownWhole.graphic.value === key.value);
const drawnSubjects = computed(() => (showAll.value ? subjects.value : subjects.value.slice(0, FIRST_SUBJECTS)));

/**
 * The drawn HSPs as columns, row after row and, in a row, from left to right (the order that Left
 * and Right follow): the table row, the smaller and larger query positions, the bin and the rank.
 */
interface Layout {
  readonly rows: readonly SubjectEntry[];
  /** Where each row's HSPs start in the columns (one more entry than rows, the end). */
  readonly start: Int32Array;
  readonly count: number;
  readonly tableRow: Int32Array;
  readonly q0: Float64Array;
  readonly q1: Float64Array;
  readonly bin: Uint8Array;
  readonly rank: Int32Array;
  readonly rowOf: Int32Array;
  readonly byRank: ReadonlyMap<number, number>;
}
const layout = computed<Layout>(() => {
  const table = loaded.value?.index.table;
  const rows = table === undefined ? [] : drawnSubjects.value;
  const count = rows.reduce((sum, subject) => sum + subject.rows.length, 0);
  const l = {
    rows,
    start: new Int32Array(rows.length + 1),
    count,
    tableRow: new Int32Array(count),
    q0: new Float64Array(count),
    q1: new Float64Array(count),
    bin: new Uint8Array(count),
    rank: new Int32Array(count),
    rowOf: new Int32Array(count),
    byRank: new Map<number, number>(),
  };
  let k = 0;
  rows.forEach((subject, r) => {
    l.start[r] = k;
    const query = (row: number) => interval(table!.qStart[row]!, table!.qEnd[row]!);
    const order = [...subject.rows].sort((a, b) => query(a).from - query(b).from || table!.rank[a]! - table!.rank[b]!);
    for (const row of order) {
      const { from, to } = query(row);
      l.tableRow[k] = row;
      l.q0[k] = from;
      l.q1[k] = to;
      l.bin[k] = scoreBin(table!.bitScore[row]!);
      l.rank[k] = table!.rank[row]!;
      l.rowOf[k] = r;
      l.byRank.set(table!.rank[row]!, k);
      k++;
    }
  });
  l.start[rows.length] = k;
  return l;
});
const selectedIndex = computed(() => {
  const id = props.state.hsp;
  return id === undefined || id.qIdx !== props.state.qIdx ? -1 : (layout.value.byRank.get(id.rank) ?? -1);
});

const canvasWidth = computed(() => Math.max(160, Math.min(MAX_WIDTH, areaWidth.value)));
const canvasHeight = computed(() => HEADER_PX + Math.max(1, Math.min(VISIBLE_ROWS, layout.value.rows.length)) * ROW_PX);
const contentHeight = computed(() => Math.max(canvasHeight.value, HEADER_PX + layout.value.rows.length * ROW_PX));

// --- geometry ------------------------------------------------------------------------------------

/** The x of a position boundary: 0 is the left end of the query bar, the length its right end. */
const xOf = (position: number): number => SIDE + (position / queryLength.value) * (canvasWidth.value - 2 * SIDE);
/** The left and right x of an HSP's bar: from before its first letter to after its last, at least 2 px. */
function barX(k: number): [number, number] {
  const l = layout.value;
  const x0 = xOf(l.q0[k]! - 1);
  return [x0, Math.max(x0 + 2, xOf(l.q1[k]!))];
}
const rowY = (r: number, top = scrollTop.value): number => HEADER_PX + r * ROW_PX + ROW_PX / 2 - top;
/** The rows in view: from the first to before the last. */
function rowsInView(top = scrollTop.value): [number, number] {
  const rows = layout.value.rows.length;
  const first = Math.max(0, Math.floor(top / ROW_PX));
  return [Math.min(first, rows), Math.min(rows, Math.ceil((top + canvasHeight.value - HEADER_PX) / ROW_PX))];
}

/** The HSP under a point of the canvas: in the row under it, the bar nearest (within PICK_PX), or -1. */
function hit(x: number, y: number): number {
  if (y < HEADER_PX) return -1;
  const l = layout.value;
  const r = Math.floor((y - HEADER_PX + scrollTop.value) / ROW_PX);
  if (r < 0 || r >= l.rows.length) return -1;
  let best = -1;
  let bestDistance = PICK_PX;
  for (let k = l.start[r]!; k < l.start[r + 1]!; k++) {
    const [x0, x1] = barX(k);
    const d = x < x0 ? x0 - x : x > x1 ? x - x1 : 0;
    if (d <= bestDistance) {
      best = k;
      bestDistance = d;
    }
  }
  return best;
}

// --- drawing -------------------------------------------------------------------------------------

let frame = 0;
function schedule(): void {
  if (frame === 0) frame = requestAnimationFrame(render);
}

let requested = '';
/** Frames drawn (`data-drawn`). */
let drawn = 0;
function render(): void {
  frame = 0;
  const top = scroller.value?.scrollTop ?? 0;
  if (top !== scrollTop.value) {
    // The rows moved under the pointer: the next move of the pointer hovers again.
    scrollTop.value = top;
    hovered.value = -1;
  }
  draw();
  // The frames drawn, for the measurements (tests/e2e/results-measure.spec.ts): set here, not
  // through Vue, so that it is in the page in the frame that drew it.
  if (canvas.value !== undefined) canvas.value.dataset['drawn'] = String(++drawn);
  placePopover();
  // The outfmt 0 headings (the descriptions of the popover) of the rows in view.
  const [first, last] = rowsInView();
  const rows = layout.value.rows.slice(first, last);
  const key = rows.map((subject) => subject.sIdx).join(',');
  if (key !== requested) {
    requested = key;
    props.results.requestHeadings(rows.map((subject) => subject.sIdx));
  }
}

function draw(): void {
  const element = canvas.value;
  if (element === undefined) return;
  const [W, H] = [canvasWidth.value, canvasHeight.value];
  const ratio = window.devicePixelRatio || 1;
  if (element.width !== Math.round(W * ratio)) element.width = Math.round(W * ratio);
  if (element.height !== Math.round(H * ratio)) element.height = Math.round(H * ratio);
  const context = element.getContext('2d');
  if (context === null) return;
  context.setTransform(ratio, 0, 0, ratio, 0, 0);
  context.clearRect(0, 0, W, H);
  context.fillStyle = '#ffffff';
  context.fillRect(0, 0, W, H);
  const l = layout.value;

  // The rows in view, under the header.
  context.save();
  context.beginPath();
  context.rect(0, HEADER_PX, W, H - HEADER_PX);
  context.clip();
  const [first, last] = rowsInView();
  if (l.count === 0) {
    context.fillStyle = MUTED;
    context.font = FONT;
    context.textAlign = 'left';
    context.textBaseline = 'middle';
    context.fillText('No HSPs to show.', SIDE, rowY(0));
  }
  // Each subject's HSPs joined by a thin grey line, from its first bar to its last.
  context.strokeStyle = JOIN;
  context.lineWidth = 1;
  context.beginPath();
  for (let r = first; r < last; r++) {
    const [from, to] = [l.start[r]!, l.start[r + 1]!];
    if (to - from < 2) continue;
    let right = 0;
    for (let k = from; k < to; k++) right = Math.max(right, barX(k)[1]);
    const y = Math.round(rowY(r)) + 0.5;
    context.moveTo(barX(from)[0], y);
    context.lineTo(right, y);
  }
  context.stroke();
  // The bars, one path per bin.
  SCORE_BINS.forEach((bin, b) => {
    context.fillStyle = bin.color;
    context.beginPath();
    for (let k = l.start[first]!; k < l.start[last]!; k++) {
      if (l.bin[k] !== b) continue;
      const [x0, x1] = barX(k);
      context.rect(x0, rowY(l.rowOf[k]!) - BAR_PX / 2, x1 - x0, BAR_PX);
    }
    context.fill();
  });
  // The selected HSP on its halo, the hovered one outlined, the keyboard's one in a dashed frame.
  const outline = (k: number, pad: number, color: string, width: number, dash: number[] = []) => {
    if (k < 0 || k >= l.count) return;
    const [x0, x1] = barX(k);
    const y = rowY(l.rowOf[k]!);
    context.strokeStyle = color;
    context.lineWidth = width;
    context.setLineDash(dash);
    context.strokeRect(x0 - pad, y - BAR_PX / 2 - pad, x1 - x0 + 2 * pad, BAR_PX + 2 * pad);
    context.setLineDash([]);
  };
  const chosen = selectedIndex.value;
  if (chosen >= 0) {
    const [x0, x1] = barX(chosen);
    context.fillStyle = HALO;
    context.fillRect(x0 - 3, rowY(l.rowOf[chosen]!) - BAR_PX / 2 - 3, x1 - x0 + 6, BAR_PX + 6);
    context.fillStyle = SCORE_BINS[l.bin[chosen]!]!.color;
    context.fillRect(x0, rowY(l.rowOf[chosen]!) - BAR_PX / 2, x1 - x0, BAR_PX);
    outline(chosen, 1.5, TEXT, 1);
  }
  if (hovered.value !== chosen) outline(hovered.value, 1.5, TEXT, 1);
  if (keyboard.value) outline(cursor.value, 3.5, ACCENT, 1.5, [3, 2]);
  context.restore();

  // The header: the query bar and its ruler.
  context.fillStyle = ACCENT;
  context.fillRect(SIDE, QUERY_TOP, W - 2 * SIDE, QUERY_PX);
  context.fillStyle = '#ffffff';
  context.font = BOLD;
  context.textAlign = 'center';
  context.textBaseline = 'middle';
  context.fillText('Query', W / 2, QUERY_TOP + QUERY_PX / 2);
  context.strokeStyle = TEXT;
  context.fillStyle = TEXT;
  context.font = FONT;
  context.textBaseline = 'top';
  context.lineWidth = 1;
  context.beginPath();
  const ticks = rulerTicks(queryLength.value, Math.max(2, Math.floor((W - 2 * SIDE) / 70) + 1));
  for (const position of ticks) {
    const x = Math.round(xOf(position - 0.5)) + 0.5;
    context.moveTo(x, QUERY_TOP + QUERY_PX);
    context.lineTo(x, QUERY_TOP + QUERY_PX + TICK_PX);
    const text = formatCount(position);
    const half = context.measureText(text).width / 2;
    context.fillText(text, Math.max(half + 2, Math.min(W - half - 2, x)), LABEL_TOP);
  }
  context.stroke();
}

// --- pointer and keyboard ------------------------------------------------------------------------

function idOf(k: number): HspId {
  const table = loaded.value!.index.table;
  const row = layout.value.tableRow[k]!;
  return { runId: loaded.value!.run.snapshot.runId, qIdx: table.qIdx[row]!, rank: table.rank[row]! };
}

/** Selects an HSP and asks for its alignment. */
function choose(k: number): void {
  if (k < 0 || k >= layout.value.count) return;
  const id = idOf(k);
  cursor.value = k;
  props.results.selectHsp(id);
  emit('show-alignment', id);
}

function point(event: MouseEvent): { x: number; y: number } {
  const rect = canvas.value!.getBoundingClientRect();
  return { x: event.clientX - rect.left, y: event.clientY - rect.top };
}
function onPointerMove(event: PointerEvent): void {
  const { x, y } = point(event);
  const k = hit(x, y);
  if (k !== hovered.value) hovered.value = k;
  if (canvas.value !== undefined) canvas.value.style.cursor = k >= 0 ? 'pointer' : '';
}
function onPointerLeave(): void {
  hovered.value = -1;
}
function onClick(event: MouseEvent): void {
  const { x, y } = point(event);
  choose(hit(x, y));
}

function onFocus(): void {
  let visible = true;
  try {
    visible = canvas.value?.matches(':focus-visible') ?? true;
  } catch {
    // A browser without :focus-visible shows the keyboard's frame on every focus.
  }
  keyboard.value = visible;
  if (cursor.value < 0 && layout.value.count > 0) cursor.value = Math.max(0, selectedIndex.value);
  if (visible && cursor.value >= 0) reveal(layout.value.rowOf[cursor.value]!);
}
function onBlur(): void {
  keyboard.value = false;
}

/** The HSP of row `r` nearest to `x` (the middle of its bar), or -1 for a row without HSPs. */
function nearest(r: number, x: number): number {
  const l = layout.value;
  let best = -1;
  let bestDistance = Infinity;
  for (let k = l.start[r]!; k < l.start[r + 1]!; k++) {
    const [x0, x1] = barX(k);
    const d = Math.abs((x0 + x1) / 2 - x);
    if (d < bestDistance) {
      best = k;
      bestDistance = d;
    }
  }
  return best;
}

/** Up and Down move between subject rows, Left and Right between a row's HSPs, Home and End to its ends; Enter chooses. */
function onKey(event: KeyboardEvent): void {
  if (event.ctrlKey || event.metaKey || event.altKey) return;
  const l = layout.value;
  if (l.count === 0) return;
  let k = cursor.value < 0 ? Math.max(0, selectedIndex.value) : cursor.value;
  const r = l.rowOf[k]!;
  switch (event.key) {
    case 'ArrowDown':
    case 'ArrowUp': {
      const direction = event.key === 'ArrowDown' ? 1 : -1;
      const [x0, x1] = barX(k);
      for (let next = r + direction; next >= 0 && next < l.rows.length; next += direction) {
        const found = nearest(next, (x0 + x1) / 2);
        if (found >= 0) {
          k = found;
          break;
        }
      }
      break;
    }
    case 'ArrowLeft':
      k = Math.max(l.start[r]!, k - 1);
      break;
    case 'ArrowRight':
      k = Math.min(l.start[r + 1]! - 1, k + 1);
      break;
    case 'Home':
      k = l.start[r]!;
      break;
    case 'End':
      k = l.start[r + 1]! - 1;
      break;
    case 'Enter':
      event.preventDefault();
      choose(k);
      return;
    default:
      return;
  }
  event.preventDefault();
  keyboard.value = true;
  cursor.value = k;
  reveal(l.rowOf[k]!);
}

/** Scrolls the rows so that row `r` is in view. */
function reveal(r: number): void {
  const element = scroller.value;
  if (element === undefined) return;
  const view = canvasHeight.value - HEADER_PX;
  const top = r * ROW_PX;
  if (top < element.scrollTop) element.scrollTop = top;
  else if (top + ROW_PX > element.scrollTop + view) element.scrollTop = top + ROW_PX - view;
  schedule();
}

// --- the popover ---------------------------------------------------------------------------------

/** The HSP that the popover shows: the hovered one, else the keyboard's. */
const shown = computed(() => (hovered.value >= 0 ? hovered.value : keyboard.value ? cursor.value : -1));
const info = computed(() => {
  const k = shown.value;
  const l = layout.value;
  if (k < 0 || k >= l.count) return undefined;
  const subject = l.rows[l.rowOf[k]!]!;
  const fields = props.results.row(l.tableRow[k]!);
  const heading = subject.inOutfmt0 ? props.state.headings.get(subject.sIdx) : undefined;
  return {
    hsp: `${props.state.qIdx}:${l.rank[k]}`,
    number: l.rank[k]! + 1,
    sseqid: fields.sseqid,
    title: heading === undefined ? undefined : headingTitle(heading),
    bitscore: fields.bitscore,
    evalue: fields.evalue,
  };
});
/** What a screen reader says as the arrow keys move. */
const spoken = computed(() => {
  const value = keyboard.value && shown.value === cursor.value ? info.value : undefined;
  if (value === undefined) return '';
  const row = layout.value.rowOf[shown.value]! + 1;
  const title = value.title === undefined ? '' : `, ${value.title}`;
  return `Subject ${row} of ${layout.value.rows.length}: ${value.sseqid}${title}. HSP ${value.number}, bit score ${value.bitscore}, E value ${value.evalue}.`;
});

/** Puts the popover next to its HSP's bar, inside the figure. */
function placePopover(): void {
  const element = popover.value;
  const k = shown.value;
  if (element === undefined || k < 0 || k >= layout.value.count) return;
  const [, x1] = barX(k);
  const y = Math.max(HEADER_PX, Math.min(canvasHeight.value, rowY(layout.value.rowOf[k]!)));
  const next = placeBox({ x: x1, y }, { width: element.offsetWidth, height: element.offsetHeight }, { width: canvasWidth.value, height: canvasHeight.value }, 8);
  if (next.left !== popoverAt.value.left || next.top !== popoverAt.value.top) popoverAt.value = next;
}

// --- life cycle ----------------------------------------------------------------------------------

let observer: ResizeObserver | undefined;
onMounted(() => {
  // The scroller's content box: its width leaves out a scroll bar, so the canvas never overflows sideways.
  observer = new ResizeObserver((entries) => {
    const next = Math.floor(entries[0]?.contentRect.width ?? areaWidth.value);
    if (next > 0 && next !== areaWidth.value) areaWidth.value = next;
  });
  if (scroller.value !== undefined) {
    areaWidth.value = Math.max(160, Math.floor(scroller.value.clientWidth));
    observer.observe(scroller.value);
  }
  schedule();
});
onUnmounted(() => {
  observer?.disconnect();
  if (frame !== 0) cancelAnimationFrame(frame);
});

// Another run or query starts at the top (with its first subjects, shownWhole.ts).
watch(
  () => `${props.state.runId}|${props.state.qIdx}`,
  () => {
    cursor.value = -1;
    hovered.value = -1;
    if (scroller.value !== undefined) scroller.value.scrollTop = 0;
    schedule();
  },
);
// New rows (filters, sort, "Show all") keep the keyboard's HSP where it is still drawn.
watch(layout, (now, before) => {
  const rank = cursor.value >= 0 && cursor.value < before.count ? before.rank[cursor.value] : undefined;
  cursor.value = rank === undefined ? -1 : (now.byRank.get(rank) ?? -1);
  hovered.value = -1;
});
watch([layout, canvasWidth, canvasHeight, selectedIndex, hovered, cursor, keyboard, info], () => schedule(), { flush: 'post' });

/** The HSPs drawn in view (at most 200), with the middles of their bars in CSS pixels of the canvas, for tests. */
const targets = computed(() => {
  const l = layout.value;
  const [first, last] = rowsInView(scrollTop.value);
  const out: { hsp: string; x: number; y: number; bin: number }[] = [];
  for (let k = l.start[first]!; k < l.start[last]! && out.length < 200; k++) {
    const [x0, x1] = barX(k);
    out.push({ hsp: `${props.state.qIdx}:${l.rank[k]}`, x: Math.round((x0 + x1) / 2), y: Math.round(rowY(l.rowOf[k]!)), bin: l.bin[k]! });
  }
  return JSON.stringify(out);
});
const selectedText = computed(() => (selectedIndex.value < 0 ? '' : `${props.state.qIdx}:${layout.value.rank[selectedIndex.value]}`));
</script>

<template>
  <figure class="graphic-summary" data-testid="graphic-summary">
    <div class="graphic-tools">
      <span class="graphic-hints"><span>hover to see the title</span><span>click to show alignments</span></span>
      <div class="graphic-legend" data-testid="graphic-legend" role="group" :aria-labelledby="legendTitle">
        <span :id="legendTitle" class="graphic-legend-title">Alignment Scores</span>
        <span v-for="(bin, b) in SCORE_BINS" :key="bin.label" class="graphic-legend-item" :data-bin="b"
          ><span class="graphic-swatch" :style="{ background: bin.color }" />{{ bin.label }}</span
        >
      </div>
    </div>
    <figcaption class="plot-caption">
      <span class="plot-title" data-testid="graphic-title"
        >Distribution of {{ formatCount(layout.count) }} HSPs on {{ formatCount(layout.rows.length) }} subject sequences</span
      >
      <span class="muted small"
        >Query {{ record?.id }} ({{ formatCount(queryLength) }} {{ unit }}). One row per subject, in the order of the Descriptions; each bar is an HSP
        at its query coordinates.</span
      >
    </figcaption>
    <div class="graphic-frame">
      <div ref="scroller" class="graphic-scroll" :style="{ height: `${canvasHeight}px` }" @scroll="schedule">
        <div class="graphic-spacer" :style="{ height: `${contentHeight}px` }">
          <canvas
            ref="canvas"
            class="graphic-canvas"
            tabindex="0"
            role="application"
            :aria-label="`Graphic Summary of ${layout.count} HSPs on ${layout.rows.length} subjects. Use the up and down arrow keys to move between subjects, left and right between a subject's HSPs, and Enter to show the alignment.`"
            :style="{ width: `${canvasWidth}px`, height: `${canvasHeight}px` }"
            data-testid="graphic-canvas"
            :data-rows="layout.rows.length"
            :data-hsps="layout.count"
            :data-selected="selectedText"
            :data-targets="targets"
            @pointermove="onPointerMove"
            @pointerleave="onPointerLeave"
            @click="onClick"
            @keydown="onKey"
            @focus="onFocus"
            @blur="onBlur"
          />
        </div>
      </div>
      <div
        v-if="info"
        ref="popover"
        class="graphic-popover"
        aria-hidden="true"
        data-testid="graphic-popover"
        :data-hsp="info.hsp"
        :style="{ left: `${popoverAt.left}px`, top: `${popoverAt.top}px` }"
      >
        <p class="graphic-popover-title">{{ info.sseqid }}</p>
        <p v-if="info.title !== undefined">{{ info.title }}</p>
        <p>HSP {{ info.number }} · Bit score {{ info.bitscore }} · E value {{ info.evalue }}</p>
      </div>
      <p class="visually-hidden" aria-live="polite">{{ spoken }}</p>
    </div>
    <p v-if="!showAll && subjects.length > FIRST_SUBJECTS" class="graphic-more">
      <button type="button" data-testid="graphic-show-all" @click="shownWhole.graphic.value = key">Show all {{ formatCount(subjects.length) }}</button>
    </p>
  </figure>
</template>
