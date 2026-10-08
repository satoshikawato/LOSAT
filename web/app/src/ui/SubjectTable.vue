<script setup lang="ts">
// The subjects of the selected query (docs/web/results_columns.md "Subject 一覧"): the
// values of the subject's first HSP as outfmt 6 wrote them (the HSP whose score outfmt 0's
// description table shows), the title of its outfmt 0 heading, and counts. Sorting uses
// the engine's values; ties keep the engine's order.
import { computed, watch } from 'vue';
import type { ResultsBrowser, ResultsState } from '../application/results';
import { headingTitle } from '../domain/outfmt0';
import type { SubjectSortKey } from '../domain/result-index';
import { formatCount } from './format';
import SortButton from './SortButton.vue';
import VirtualRows from './VirtualRows.vue';

const props = defineProps<{ results: ResultsBrowser; state: ResultsState }>();
const ROW_PX = 34;
const MAX_ROWS = 10;
const selectedPosition = computed(() => props.state.subjects.findIndex((subject) => subject.sIdx === props.state.sIdx));
const unit = computed(() => props.state.loaded?.units.subject ?? '');

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
    <div class="table-scroll">
      <div class="table-head subject-grid" role="row">
        <SortButton label="#" sort-key="order" :sort="state.subjectSort" scope="subject" @sort="sortBy" />
        <span role="columnheader">Subject</span>
        <span role="columnheader">Description (outfmt 0)</span>
        <SortButton :label="`Length (${unit})`" sort-key="length" :sort="state.subjectSort" scope="subject" first-descending @sort="sortBy" />
        <SortButton label="Score (bits)" sort-key="bitScore" :sort="state.subjectSort" scope="subject" first-descending @sort="sortBy" />
        <SortButton label="E value" sort-key="eValue" :sort="state.subjectSort" scope="subject" @sort="sortBy" />
        <SortButton label="HSPs" sort-key="hsps" :sort="state.subjectSort" scope="subject" first-descending @sort="sortBy" />
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
      >
        <template #row="{ position }">
          <button
            type="button"
            class="table-row subject-grid"
            :class="{ selected: state.subjects[position]!.sIdx === state.sIdx }"
            :aria-pressed="state.subjects[position]!.sIdx === state.sIdx"
            :data-testid="`subject-row-${state.subjects[position]!.sIdx}`"
            :data-order="state.subjects[position]!.order"
            @click="results.selectSubject(state.subjects[position]!.sIdx)"
          >
            <span class="num">{{ state.subjects[position]!.order }}</span>
            <span class="cell-id" :title="state.subjects[position]!.first.sseqid" data-field="sseqid">{{
              state.subjects[position]!.first.sseqid
            }}</span>
            <span
              class="cell-text"
              :class="{ muted: !state.subjects[position]!.inOutfmt0 }"
              :title="description(state.subjects[position]!.sIdx, state.subjects[position]!.inOutfmt0)"
              data-field="description"
              >{{ description(state.subjects[position]!.sIdx, state.subjects[position]!.inOutfmt0) }}</span
            >
            <span class="num" data-field="length">{{ formatCount(state.subjects[position]!.length) }}</span>
            <span class="num" data-field="bitscore">{{ state.subjects[position]!.first.bitscore }}</span>
            <span class="num" data-field="evalue">{{ state.subjects[position]!.first.evalue }}</span>
            <span class="num" data-field="hsps">
              {{ formatCount(state.subjects[position]!.rows.length) }}<template
                v-if="state.subjects[position]!.rows.length !== state.subjects[position]!.hspCount"
                >/{{ formatCount(state.subjects[position]!.hspCount) }}</template
              ><span v-if="state.subjects[position]!.atHspLimit" class="badge" title="The subject has as many HSPs as -max_hsps keeps">max</span>
            </span>
          </button>
        </template>
      </VirtualRows>
    </div>
  </div>
</template>
