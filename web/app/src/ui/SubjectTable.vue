<script setup lang="ts">
// The subjects of the selected query (docs/web/results_columns.md "Subject 一覧"): the
// values of the subject's first HSP as outfmt 6 wrote them (the HSP whose score outfmt 0's
// description table shows), the title of its outfmt 0 heading, and counts. Sorting uses
// the engine's values; ties keep the engine's order.
import { computed, ref, watch } from 'vue';
import type { ResultsBrowser, ResultsState, SubjectEntry } from '../application/results';
import { headingTitle } from '../domain/outfmt0';
import type { SubjectSortKey } from '../domain/result-index';
import { formatCount } from './format';
import SortButton from './SortButton.vue';
import { useSideScroll } from './useSideScroll';
import VirtualRows from './VirtualRows.vue';

const props = defineProps<{ results: ResultsBrowser; state: ResultsState }>();
const ROW_PX = 34;
const MAX_ROWS = 10;
const selectedPosition = computed(() => props.state.subjects.findIndex((subject) => subject.sIdx === props.state.sIdx));
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
  ([subjects]) => props.results.requestHeadings(subjects.slice(0, 3 * MAX_ROWS).map((subject) => subject.sIdx)),
  { immediate: true },
);
function onScroll(position: number): void {
  props.results.requestHeadings(props.state.subjects.slice(position, position + 2 * MAX_ROWS).map((s) => s.sIdx));
}

/** The HSP count of a subject's row: the HSPs shown, and of how many when the view filters hide some. */
function hspCount(subject: SubjectEntry): string {
  const shown = formatCount(subject.rows.length);
  return subject.rows.length === subject.hspCount ? shown : `${shown}/${formatCount(subject.hspCount)}`;
}

function description(sIdx: number, inOutfmt0: boolean): string {
  if (!inOutfmt0) return 'not in outfmt 0';
  const heading = props.state.headings.get(sIdx);
  return heading === undefined ? '…' : headingTitle(heading);
}
</script>

<template>
  <div class="result-table subject-table" data-testid="subject-table">
    <h3>
      Subjects
      <span class="muted small">{{ formatCount(state.subjects.length) }} shown</span>
    </h3>
    <p v-if="sideScroll" class="table-hint muted small" data-testid="subject-table-scroll-hint">Scroll the table sideways for more columns →</p>
    <div ref="scroller" class="table-scroll" :style="{ '--row-gutter': `${gutter}px` }">
      <div class="table-head subject-grid" role="row">
        <SortButton class="num-head" data-col="order" label="#" sort-key="order" :sort="state.subjectSort" scope="subject" @sort="sortBy" />
        <span data-col="sseqid" role="columnheader">Subject</span>
        <span data-col="description" role="columnheader">Description (outfmt 0)</span>
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
      </div>
      <VirtualRows
        :count="state.subjects.length"
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
          <button
            v-for="subject in [state.subjects[position]!]"
            :key="subject.sIdx"
            type="button"
            class="table-row subject-grid"
            :class="{ selected: subject.sIdx === state.sIdx }"
            :aria-pressed="subject.sIdx === state.sIdx"
            :data-testid="`subject-row-${subject.sIdx}`"
            :data-order="subject.order"
            @click="results.selectSubject(subject.sIdx)"
          >
            <span class="num" data-field="order" :title="String(subject.order)">{{ subject.order }}</span>
            <span :title="subject.first.sseqid" data-field="sseqid">{{ subject.first.sseqid }}</span>
            <span
              :class="{ muted: !subject.inOutfmt0 }"
              :title="description(subject.sIdx, subject.inOutfmt0)"
              data-field="description"
              >{{ description(subject.sIdx, subject.inOutfmt0) }}</span
            >
            <span class="num" data-field="length" :title="formatCount(subject.length)">{{ formatCount(subject.length) }}</span>
            <span class="num" data-field="bitscore" :title="subject.first.bitscore">{{ subject.first.bitscore }}</span>
            <span class="num" data-field="evalue" :title="subject.first.evalue">{{ subject.first.evalue }}</span>
            <span class="num" data-field="hsps" :title="hspCount(subject) + (subject.atHspLimit ? ' max' : '')">
              {{ hspCount(subject) }}<span v-if="subject.atHspLimit" class="badge" title="The subject has as many HSPs as -max_hsps keeps">max</span>
            </span>
          </button>
        </template>
      </VirtualRows>
    </div>
  </div>
</template>
