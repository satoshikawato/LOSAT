<script setup lang="ts">
// The Descriptions: the subjects of the selected query (docs/web/results_columns.md "Subject
// 一覧"), in the columns of NCBI's "Sequences producing significant alignments" where LOSAT has
// the values (docs/web/ncbi_ui_mapping.md "Descriptions"): the values of the subject's first HSP as
// outfmt 6 wrote them (the HSP whose score outfmt 0's description table shows), the title of its
// outfmt 0 heading, counts, and the subject's ID last (NCBI's Accession). Sorting uses the
// engine's values; ties keep the engine's order. As NCBI's, each row has a mark beside it (a
// check box of its own: the row stays one button that selects the subject), with "select all"
// and the count of marked rows; "Add to candidates" adds every HSP of the marked subjects in
// this query (S14).
//
// On the one page of the results (ClassicResults.vue), the list is cut as NCBI's: the first 100
// subjects in the order shown, "Show all N" for the rest, until another run or query (shownWhole.ts).
// "select all" marks the subjects listed.
import { computed, onMounted, ref, watch } from 'vue';
import type { HspId, ResultsBrowser, ResultsState, SubjectEntry } from '../application/results';
import { headingTitle } from '../domain/outfmt0';
import type { SubjectSortKey } from '../domain/result-index';
import { formatCount, formatCounted } from './format';
import { shownWhole, wholeKey } from './shownWhole';
import SortButton from './SortButton.vue';
import { useNarrow } from './useNarrow';
import { focusPressed, useSideScroll } from './useSideScroll';
import VirtualRows from './VirtualRows.vue';

const props = defineProps<{ results: ResultsBrowser; state: ResultsState }>();
/** A row chosen (after its subject is selected): the one page brings the subject's alignments into view. */
const emit = defineEmits<{ 'add-candidates': [ids: readonly HspId[]]; chosen: [sIdx: number] }>();
const ROW_PX = 28;
/**
 * Rows in view before the box scrolls: about two dozen on a desktop screen, so that the page reads
 * at a glance as NCBI's list (screen review 2 I1); a phone keeps a short box.
 */
const MAX_ROWS = computed(() => (narrow.value ? 10 : 22));
/** Subjects listed until "Show all" (NCBI lists 100 to a page). */
const FIRST_SUBJECTS = 100;
const narrow = useNarrow();
const key = computed(() => wholeKey(props.state.runId, props.state.qIdx));
const showAll = computed(() => shownWhole.descriptions.value === key.value);
const listed = computed(() => (showAll.value ? props.state.subjects.length : Math.min(FIRST_SUBJECTS, props.state.subjects.length)));
/** The selected subject's row, or -1 where the list does not reach it (it is not scrolled to then). */
const selectedPosition = computed(() => {
  const position = props.state.subjects.findIndex((subject) => subject.sIdx === props.state.sIdx);
  return position < listed.value ? position : -1;
});
const cut = computed(() => listed.value < props.state.subjects.length);
// The marks follow the listed subjects: a mark on one that a sort or a filter moves out of the
// listed ones is dropped (ResultsBrowser.setListLimit).
watch(
  () => (showAll.value ? undefined : FIRST_SUBJECTS),
  (limit) => props.results.setListLimit(limit),
  { immediate: true },
);
const unit = computed(() => props.state.loaded?.units.subject ?? '');
const scroller = ref<HTMLElement>();
const sideScroll = useSideScroll(scroller);
const gutter = ref(0);

function sortBy(key: SubjectSortKey, descending: boolean): void {
  props.results.setSubjectSort({ key, descending });
}

// The descriptions are read from outfmt 0 for the subjects in view (and a screenful more).
watch(
  () => [props.state.subjects, props.state.qIdx] as const,
  ([subjects]) => props.results.requestHeadings(subjects.slice(0, 3 * MAX_ROWS.value).map((subject) => subject.sIdx)),
  { immediate: true },
);
function onScroll(position: number): void {
  props.results.requestHeadings(props.state.subjects.slice(position, position + 2 * MAX_ROWS.value).map((s) => s.sIdx));
}

/** The HSP count of a subject's row: the HSPs shown, and of how many when the view filters hide some. */
function hspCount(subject: SubjectEntry): string {
  const shown = formatCount(subject.rows.length);
  return subject.rows.length === subject.hspCount ? shown : `${shown}/${formatCount(subject.hspCount)}`;
}

/**
 * The Subject ID column is as wide as its longest value at the table's type (within 10ch and 24ch
 * in styles.css; a longer ID ends in an ellipsis and has its title), so that the Description has
 * the spare width (W4b screen review L8). Only the longest IDs by their count of letters are
 * measured.
 */
const idWidth = ref<string>();
let measurer: CanvasRenderingContext2D | null | undefined;
function measureIds(): void {
  const element = scroller.value;
  if (element === undefined) return;
  measurer ??= document.createElement('canvas').getContext('2d');
  if (measurer === null) return;
  const subjects = props.state.subjects;
  let longest = 0;
  for (const subject of subjects) longest = Math.max(longest, subject.first.sseqid.length);
  const style = getComputedStyle(element);
  measurer.font = `${style.fontStyle} ${style.fontWeight} ${style.fontSize} ${style.fontFamily}`;
  let width = 0;
  let measured = 0;
  for (const subject of subjects) {
    if (subject.first.sseqid.length < longest - 4) continue;
    width = Math.max(width, measurer.measureText(subject.first.sseqid).width);
    if (++measured === 64) break;
  }
  idWidth.value = `${Math.ceil(width) + 1}px`;
}
onMounted(measureIds);
watch(() => props.state.subjects, measureIds, { flush: 'post' });

/** "select all" is checked when every listed subject is marked, and mixed when some are. */
const allMarked = computed(() => {
  const marked = props.state.marked;
  if (listed.value === 0 || marked.size < listed.value) return false;
  for (let i = 0; i < listed.value; i++) if (!marked.has(props.state.subjects[i]!.sIdx)) return false;
  return true;
});
const someMarked = computed(() => props.state.marked.size > 0 && !allMarked.value);

function mark(sIdx: number, event: Event): void {
  props.results.markSubjects([sIdx], (event.target as HTMLInputElement).checked);
}
/** NCBI's "select all": the subjects listed (the first 100 until "Show all"). */
function markListed(event: Event): void {
  const sIdxs = props.state.subjects.slice(0, listed.value).map((subject) => subject.sIdx);
  props.results.markSubjects(sIdxs, (event.target as HTMLInputElement).checked);
}

function choose(sIdx: number): void {
  props.results.selectSubject(sIdx);
  emit('chosen', sIdx);
}

function description(sIdx: number, inOutfmt0: boolean): string {
  if (!inOutfmt0) return 'not in outfmt 0';
  const heading = props.state.headings.get(sIdx);
  return heading === undefined ? '…' : headingTitle(heading);
}
</script>

<template>
  <div class="result-table subject-table" data-testid="subject-table">
    <div class="tool-band">
      <h4>
        Sequences producing significant alignments
        <span class="muted small" data-testid="subject-count">{{
            cut ? `${formatCount(listed)} of ${formatCount(state.subjects.length)} listed` : `${formatCount(state.subjects.length)} shown`
          }}</span>
      </h4>
      <div class="descriptions-tools">
        <label class="check">
          <input
            type="checkbox"
            :checked="allMarked"
            :indeterminate="someMarked"
            data-testid="descriptions-select-all"
            @change="markListed"
          />
          {{ cut ? 'select all listed' : 'select all' }}
        </label>
        <span class="muted" data-testid="descriptions-selected" aria-live="polite">{{ formatCounted(state.marked.size, 'sequence') }} selected</span>
        <button
          type="button"
          :disabled="state.marked.size === 0"
          data-testid="descriptions-add-candidates"
          @click="emit('add-candidates', results.markedHspIds())"
        >
          Add to candidates
        </button>
      </div>
    </div>
    <p v-if="sideScroll" class="table-hint muted small" data-testid="subject-table-scroll-hint">Scroll the table sideways for more columns →</p>
    <div ref="scroller" class="table-scroll" :style="{ '--row-gutter': `${gutter}px`, '--id-width': idWidth }">
      <div class="table-head subject-grid marked-head" role="row">
        <SortButton class="num-head" data-col="order" label="#" sort-key="order" :sort="state.subjectSort" scope="subject" @sort="sortBy" />
        <span data-col="description" role="columnheader">Description</span>
        <SortButton
          class="num-head"
          data-col="bitscore"
          label="Score (bits)"
          sort-key="bitScore"
          :sort="state.subjectSort"
          scope="subject"
          first-descending
          @sort="sortBy"
        />
        <SortButton class="num-head" data-col="evalue" label="E value" sort-key="eValue" :sort="state.subjectSort" scope="subject" @sort="sortBy" />
        <SortButton
          class="num-head"
          data-col="hsps"
          label="HSPs"
          sort-key="hsps"
          :sort="state.subjectSort"
          scope="subject"
          first-descending
          @sort="sortBy"
        />
        <SortButton
          class="num-head"
          data-col="length"
          :label="`Length (${unit})`"
          sort-key="length"
          :sort="state.subjectSort"
          scope="subject"
          first-descending
          @sort="sortBy"
        />
        <span data-col="sseqid" role="columnheader">Subject ID</span>
      </div>
      <VirtualRows
        :count="listed"
        :row-px="ROW_PX"
        :max-rows="MAX_ROWS"
        :reveal="selectedPosition"
        :reveal-key="`${state.runId}:${state.qIdx}:${state.sIdx}`"
        :order-key="`${state.subjectSort.key}:${state.subjectSort.descending}`"
        label="Subjects"
        testid="subject-list"
        @scroll="onScroll(Math.floor(($event.target as HTMLElement).scrollTop / ROW_PX))"
        @gutter="gutter = $event"
      >
        <template #row="{ position }">
          <div
            v-for="subject in [state.subjects[position]!]"
            :key="subject.sIdx"
            class="marked-line"
            :class="{ selected: subject.sIdx === state.sIdx }"
          >
            <label class="row-mark">
              <input
                type="checkbox"
                :checked="state.marked.has(subject.sIdx)"
                :aria-label="`Select ${subject.first.sseqid}`"
                :data-testid="`subject-mark-${subject.sIdx}`"
                @change="mark(subject.sIdx, $event)"
              />
            </label>
            <button
              type="button"
              class="table-row subject-grid"
              :class="{ selected: subject.sIdx === state.sIdx }"
              :aria-pressed="subject.sIdx === state.sIdx"
              :data-testid="`subject-row-${subject.sIdx}`"
              :data-order="subject.order"
              @mousedown="focusPressed"
              @click="choose(subject.sIdx)"
            >
              <span class="num" data-field="order" :title="String(subject.order)">{{ subject.order }}</span>
              <span
                :class="{ muted: !subject.inOutfmt0 }"
                :title="description(subject.sIdx, subject.inOutfmt0)"
                data-field="description"
                >{{ description(subject.sIdx, subject.inOutfmt0) }}</span
              >
              <span class="num" data-field="bitscore" :title="subject.first.bitscore">{{ subject.first.bitscore }}</span>
              <span class="num" data-field="evalue" :title="subject.first.evalue">{{ subject.first.evalue }}</span>
              <span class="num" data-field="hsps" :title="hspCount(subject) + (subject.atHspLimit ? ' max' : '')">
                {{ hspCount(subject) }}<span v-if="subject.atHspLimit" class="badge" title="The subject has as many HSPs as -max_hsps keeps">max</span>
              </span>
              <span class="num" data-field="length" :title="formatCount(subject.length)">{{ formatCount(subject.length) }}</span>
              <span :title="subject.first.sseqid" data-field="sseqid">{{ subject.first.sseqid }}</span>
            </button>
          </div>
        </template>
      </VirtualRows>
    </div>
    <p v-if="listed < state.subjects.length" class="descriptions-more">
      <span class="muted small" data-testid="descriptions-listed">The first {{ formatCount(listed) }} of {{ formatCount(state.subjects.length) }} are listed.</span>
      <button type="button" data-testid="descriptions-show-all" @click="shownWhole.descriptions.value = key">Show all {{ formatCount(state.subjects.length) }}</button>
    </p>
  </div>
</template>
