<script setup lang="ts">
// The Alignments tab (S13b, docs/web/ncbi_ui_mapping.md "Alignments"), after NCBI's block of a
// subject: the selected subject's outfmt 0 heading as written, "Sequence ID", "Length" and
// "Number of Matches", the subject's HSP table, then one block per HSP in the engine's order,
// "Range n: a to b" with "Next Match", "Previous Match" and "First Match", and the HSP's outfmt 0
// section as written. The selected HSP's block holds W4's detail (its outfmt 6 row, and why outfmt
// 0 does not show it). Nothing of the alignment is drawn again here.
//
// A pair may have thousands of HSPs: the blocks are the selected Range and up to 25 on each side
// (more on request), and a block reads its section only when it comes into view. Sections that
// arrive in the same frame are shown together, so that the page is measured once a frame.
//
// Where NCBI has "Download" in the subject's block and "GenBank" and "Graphics" beside each
// Range, LOSAT adds the HSPs to the candidate tray: all of the subject's ("Add all matches to
// candidates"), or the Range's ("Add to candidates", "In candidates" once there; S14).
import { computed, nextTick, onMounted, onUnmounted, reactive, ref, shallowRef, watch } from 'vue';
import { candidateKey } from '../application/candidates';
import { sameHsp, type HspId, type RangeEntry, type ResultsBrowser, type ResultsState } from '../application/results';
import { windowAround } from '../domain/result-index';
import { formatCount } from './format';
import HspTable from './HspTable.vue';

const props = defineProps<{
  results: ResultsBrowser;
  state: ResultsState;
  /** The keys of the HSPs in the candidate tray. */
  inTray: ReadonlySet<string>;
}>();
const emit = defineEmits<{ descriptions: []; 'add-candidates': [ids: readonly HspId[]] }>();

/** Ranges shown on each side of the selected one, and added by "Show earlier/later matches". */
const REACH = 25;

const root = ref<HTMLElement>();
const subject = computed(() => props.state.subjects.find((s) => s.sIdx === props.state.sIdx));
/** The subject's position in the Descriptions' order. */
const position = computed(() => props.state.subjects.findIndex((s) => s.sIdx === props.state.sIdx));
const ranges = computed(() => props.state.ranges);
const selectedAt = computed(() => ranges.value.findIndex((range) => sameHsp(range.id, props.state.hsp)));
const detail = computed(() => props.state.detail);
const entry = computed(() => props.state.hsps.find((hsp) => sameHsp(hsp.id, detail.value?.id)));
const limitGiven = computed(() => props.state.loaded?.run.snapshot.argv.includes('-max_target_seqs') ?? false);

/**
 * Every HSP of the subject in this query, whatever the view filters show (a subject is added
 * whole): listed again only for another subject, not at each change of the results' state.
 */
const subjectHsps = shallowRef<readonly HspId[]>([]);
watch(
  [() => props.state.loaded, () => props.state.qIdx, () => props.state.sIdx],
  ([, qIdx, sIdx]) => (subjectHsps.value = qIdx === undefined || sIdx === undefined ? [] : props.results.hspIdsOfSubject(qIdx, sIdx)),
  { immediate: true },
);
const subjectInTray = computed(() => subjectHsps.value.length > 0 && subjectHsps.value.every((id) => props.inTray.has(candidateKey(id))));
const isInTray = (id: HspId) => props.inTray.has(candidateKey(id));

/** Adds HSPs; a button that already says "In candidates" does nothing (it keeps the focus, so it is not disabled). */
function add(ids: readonly HspId[], already: boolean): void {
  if (!already) emit('add-candidates', ids);
}

/** Why the selected HSP's view filters were cleared ("Show in results" of a candidate), while it stays selected. */
const revealedMessage = computed(() => {
  const revealed = props.state.revealed;
  return revealed?.message !== undefined && sameHsp(revealed.id, props.state.hsp) ? revealed.message : undefined;
});
/** The subject's heading in outfmt 0, as written: read with the selected HSP, or for the Descriptions. */
const heading = computed(() => {
  const sIdx = props.state.sIdx;
  if (sIdx === undefined || subject.value?.inOutfmt0 !== true) return undefined;
  const ofSubject = detail.value !== undefined && ranges.value.some((range) => sameHsp(range.id, detail.value?.id));
  return props.state.headings.get(sIdx) ?? (ofSubject ? detail.value?.heading : undefined);
});

// --- the window of Ranges ------------------------------------------------------------------------

const win = ref({ start: 0, end: 0 });
const windowKey = computed(() => `${props.state.runId}:${props.state.qIdx}:${props.state.sIdx}:${ranges.value.length}`);
// Another subject (or other filters) centres the window on the selected Range; so does a selection
// outside it. A selection inside it keeps the Ranges where they are.
watch(
  [windowKey, selectedAt],
  ([key, at], before) => {
    const inside = at >= win.value.start && at < win.value.end;
    if (before === undefined || key !== before[0] || !inside) win.value = windowAround(Math.max(0, at), ranges.value.length, REACH);
  },
  { immediate: true },
);
const shown = computed(() => ranges.value.slice(win.value.start, win.value.end));

/**
 * Applies a change of the blocks (Ranges added above, a section read) and keeps the first block
 * that starts on the screen where it was (else the one that the screen's top cuts): a block that
 * grows above it, or Ranges added above it, would push the Range being read down. The browsers'
 * own scroll anchoring is off for the blocks (styles.css), as WebKit has none.
 */
async function keepInPlace(change: () => void): Promise<void> {
  const blocks = [...(root.value?.querySelectorAll<HTMLElement>('section.range-block') ?? [])];
  const anchor =
    blocks.find((block) => block.getBoundingClientRect().top >= 0 && block.getBoundingClientRect().top < window.innerHeight) ??
    blocks.find((block) => block.getBoundingClientRect().bottom > 0);
  const before = anchor?.getBoundingClientRect().top;
  change();
  await nextTick();
  if (anchor === undefined || before === undefined || !anchor.isConnected) return;
  const moved = anchor.getBoundingClientRect().top - before;
  if (Math.abs(moved) >= 1) window.scrollBy(0, moved);
}

function showEarlier(): void {
  void keepInPlace(() => (win.value = { start: Math.max(0, win.value.start - REACH), end: win.value.end }));
}
function showLater(): void {
  win.value = { start: win.value.start, end: Math.min(ranges.value.length, win.value.end + REACH) };
}

// --- moving the selection ------------------------------------------------------------------------

const key = (id: HspId) => `${id.qIdx}-${id.rank}`;

/**
 * Brings an HSP's Range into view (the HSP chosen in the table, the plots or by the Range's
 * buttons); with `focus`, the keyboard's focus moves to its label too.
 */
async function reveal(id: HspId, focus = false): Promise<void> {
  await nextTick();
  const label = root.value?.querySelector<HTMLElement>(`[data-testid="range-${key(id)}"] .range-label`);
  if (label === null || label === undefined) return;
  label.scrollIntoView({ block: 'start' });
  if (focus) label.focus({ preventScroll: true });
}
defineExpose({ reveal });

function selectRange(range: RangeEntry | undefined): void {
  if (range === undefined) return;
  props.results.selectHsp(range.id);
  void reveal(range.id, true);
}

function selectSubject(offset: number): void {
  const next = props.state.subjects[position.value + offset];
  if (next !== undefined) props.results.selectSubject(next.sIdx);
}

// --- the sections, read as they come into view -------------------------------------------------------

/** Sections read for this run's Ranges, by `${qIdx}-${rank}`: the text, or null where it could not be read. */
const sections = reactive(new Map<string, string | null>());
const reading = new Set<string>();
/** Failed blocks that left the screen: they are read again when they come back into view. */
const away = new Set<string>();

/**
 * Sections read since the last frame. They are shown together at the next frame: each one shown
 * on its own measured the page again, and a few long sections of a 5,993-HSP pair arriving after
 * a selection took 10 to 25 ms each to lay out, which held the frame back by 100 ms (W4b F1).
 */
const arrived = new Map<string, string | null>();
let arrival = 0;
function arrive(k: string, text: string | null): void {
  arrived.set(k, text);
  if (arrival !== 0) return;
  arrival = requestAnimationFrame(() => {
    arrival = 0;
    const batch = [...arrived];
    arrived.clear();
    void keepInPlace(() => {
      // A failure that arrives after the detail showed the section does not replace it (below).
      for (const [k, value] of batch) if (value !== null || sections.get(k) === undefined) sections.set(k, value);
    });
  });
}

watch(
  () => props.state.runId,
  () => {
    sections.clear();
    reading.clear();
    away.clear();
    arrived.clear();
  },
);

function read(range: RangeEntry): void {
  const k = key(range.id);
  if (sections.has(k) || reading.has(k) || arrived.has(k)) return;
  reading.add(k);
  const runId = props.state.runId;
  props.results
    .readSection(range.id)
    .then((text) => {
      if (props.state.runId === runId) arrive(k, text ?? '');
    })
    .catch(() => {
      if (props.state.runId === runId) arrive(k, null);
    })
    .finally(() => reading.delete(k));
}

// The selected HSP's detail reads the same section. Once it shows it, the block keeps it: a
// failure of the block's own earlier read must not come back when another HSP is selected.
watch(
  () => [props.state.runId, detail.value?.id, detail.value?.state, detail.value?.section] as const,
  () => {
    const d = detail.value;
    if (d?.state !== 'ready' || d.section === undefined) return;
    const k = key(d.id);
    if (sections.get(k) === null || sections.get(k) === undefined) sections.set(k, d.section);
  },
);

/** Reads a failed section again ("Try again", or its block back in view). */
function retry(range: RangeEntry): void {
  const k = key(range.id);
  away.delete(k);
  sections.delete(k);
  read(range);
}

let observer: IntersectionObserver | undefined;
/** Observes the blocks whose section waits to be read, or could not be read (after each change of the blocks shown). */
function observe(): void {
  if (observer === undefined || root.value === undefined) return;
  observer.disconnect();
  for (const element of root.value.querySelectorAll<HTMLElement>('[data-lazy]')) observer.observe(element);
}
onMounted(() => {
  observer = new IntersectionObserver(
    (entries) => {
      for (const item of entries) {
        const k = (item.target as HTMLElement).dataset['range'];
        const range = shown.value.find((r) => key(r.id) === k);
        if (range === undefined || k === undefined) continue;
        if (sections.get(k) !== null) {
          if (item.isIntersecting) read(range);
        } else if (!item.isIntersecting) away.add(k);
        else if (away.has(k)) retry(range);
      }
    },
    { rootMargin: '300px 0px' },
  );
  observe();
});
onUnmounted(() => {
  observer?.disconnect();
  if (arrival !== 0) cancelAnimationFrame(arrival);
});
watch([shown, () => props.state.hsp, () => sections.size], () => void nextTick(observe), { flush: 'post' });

// The subject's heading, once per subject.
watch(
  () => [props.state.sIdx, subject.value?.inOutfmt0] as const,
  ([sIdx, inOutfmt0]) => {
    if (sIdx !== undefined && inOutfmt0 === true) props.results.requestHeadings([sIdx]);
  },
  { immediate: true },
);

const sectionState = (range: RangeEntry): 'pending' | 'ready' | 'failed' => {
  const text = sections.get(key(range.id));
  return text === undefined ? 'pending' : text === null ? 'failed' : 'ready';
};
</script>

<template>
  <div ref="root" class="alignments" data-testid="alignments">
    <div v-if="subject" class="alignments-subject" data-testid="alignments-subject" :data-subject="subject.sIdx">
      <div class="tool-band alignments-nav">
        <span class="alignments-position muted small">Subject {{ formatCount(position + 1) }} of {{ formatCount(state.subjects.length) }}</span>
        <span class="alignments-buttons">
          <button type="button" :disabled="position <= 0" data-testid="alignments-prev-subject" @click="selectSubject(-1)">Previous</button>
          <button type="button" :disabled="position >= state.subjects.length - 1" data-testid="alignments-next-subject" @click="selectSubject(1)">
            Next
          </button>
          <button type="button" data-testid="alignments-descriptions" @click="emit('descriptions')">Descriptions</button>
          <button
            type="button"
            :aria-disabled="subjectInTray"
            data-testid="alignments-add-subject"
            @click="add(subjectHsps, subjectInTray)"
          >
            {{ subjectInTray ? 'All matches in candidates' : 'Add all matches to candidates' }}
          </button>
        </span>
      </div>
      <pre v-if="heading !== undefined" class="output subject-heading" data-testid="detail-heading">{{ heading }}</pre>
      <!-- The fields on one line of the template: the spaces between them are part of the text. -->
      <p class="alignments-summary" data-testid="alignments-summary">
        <span><strong>Sequence ID:</strong> {{ subject.first.sseqid }}</span>  <span><strong>Length:</strong> {{ formatCount(subject.length) }}</span>  <span><strong>Number of Matches:</strong> {{ formatCount(subject.hspCount) }}<template v-if="ranges.length < subject.hspCount"> ({{ formatCount(ranges.length) }} shown)</template></span>
      </p>

      <HspTable :results="results" :state="state" @chosen="reveal($event)" />

      <p v-if="win.start > 0" class="range-more">
        <button type="button" data-testid="alignments-show-earlier" @click="showEarlier">
          Show earlier matches ({{ formatCount(win.start) }} more)
        </button>
      </p>
      <section
        v-for="range in shown"
        :key="key(range.id)"
        class="range-block"
        :class="{ selected: range === ranges[selectedAt] }"
        :data-testid="`range-${key(range.id)}`"
        :data-n="range.n"
      >
        <div class="range-head">
          <div class="range-title">
            <h4 class="range-label" tabindex="-1" data-testid="range-label">Range {{ range.n }}: {{ range.from }} to {{ range.to }}</h4>
            <button
              type="button"
              class="link"
              :aria-disabled="isInTray(range.id)"
              :data-testid="`range-add-${key(range.id)}`"
              @click="add([range.id], isInTray(range.id))"
            >
              {{ isInTray(range.id) ? 'In candidates' : 'Add to candidates' }}
            </button>
          </div>
          <span class="range-buttons">
            <button
              type="button"
              class="link"
              :disabled="range === ranges[ranges.length - 1]"
              data-testid="range-next"
              @click="selectRange(ranges[ranges.indexOf(range) + 1])"
            >
              Next Match
            </button>
            <button
              type="button"
              class="link"
              :disabled="range === ranges[0]"
              data-testid="range-previous"
              @click="selectRange(ranges[ranges.indexOf(range) - 1])"
            >
              Previous Match
            </button>
            <button type="button" class="link" :disabled="range === ranges[0]" data-testid="range-first" @click="selectRange(ranges[0])">
              First Match
            </button>
          </span>
        </div>

        <!-- The selected HSP: W4's detail, with its outfmt 6 row. -->
        <div
          v-if="detail && sameHsp(range.id, detail.id)"
          class="hsp-detail"
          data-testid="hsp-detail"
          :data-hsp="`${detail.id.qIdx}:${detail.id.rank}`"
          :data-state="detail.state"
        >
          <p v-if="revealedMessage" class="notice" data-testid="revealed-message">{{ revealedMessage }}</p>
          <p v-if="entry?.orientation === 'unknown'" class="notice" data-testid="detail-strand-note">
            This HSP covers one letter of each sequence, so its coordinates do not show its strand, and the HSP record does not hold
            it. The <code>Strand=</code> line of its outfmt 0 section below shows it.
          </p>
          <p v-if="detail.state === 'failed'" class="error" data-testid="detail-error">
            The HSP could not be read: {{ detail.error }}
            <button type="button" class="link" data-testid="detail-retry" @click="results.selectHsp(detail.id)">Try again</button>
          </p>
          <template v-if="range.inOutfmt0">
            <p v-if="detail.state === 'loading'" class="muted">Reading the alignment…</p>
            <pre v-if="detail.section !== undefined" class="output range-text" data-testid="detail-section">{{ detail.section }}</pre>
          </template>
          <p v-else class="notice" data-testid="detail-not-in-outfmt0">
            outfmt 0 does not show this HSP. It shows the alignments of the first {{ state.outfmt0Subjects }} subjects of this query
            (BLAST+'s <span class="option-name">-num_alignments</span>: 250<template v-if="limitGiven"
              >, here <span class="option-name">-max_target_seqs</span></template
            >), and this subject comes after them. The outfmt 6 row below and outfmt 7 hold the HSP.
          </p>
          <h5 class="row-label">outfmt 6 row</h5>
          <pre class="output row-text" data-testid="detail-row">{{ detail.row }}</pre>
        </div>

        <!-- Another HSP of the subject: its section, read when the block comes into view. -->
        <template v-else-if="range.inOutfmt0">
          <pre
            v-if="sectionState(range) === 'ready'"
            class="output range-text"
            data-testid="range-section"
            data-state="ready"
            >{{ sections.get(key(range.id)) }}</pre
          >
          <p
            v-else-if="sectionState(range) === 'failed'"
            class="error"
            data-testid="range-section"
            data-state="failed"
            data-lazy="failed"
            :data-range="key(range.id)"
          >
            The alignment could not be read.
            <button type="button" class="link" data-testid="range-retry" @click="retry(range)">Try again</button>
          </p>
          <p v-else class="range-pending muted" data-testid="range-section" data-state="pending" data-lazy="pending" :data-range="key(range.id)">
            Reading the alignment…
          </p>
        </template>
        <p v-else class="muted small" data-testid="range-not-in-outfmt0">outfmt 0 does not show this HSP.</p>
      </section>
      <p v-if="win.end < ranges.length" class="range-more">
        <button type="button" data-testid="alignments-show-later" @click="showLater">
          Show later matches ({{ formatCount(ranges.length - win.end) }} more)
        </button>
      </p>
    </div>
  </div>
</template>
