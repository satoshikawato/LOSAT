<script setup lang="ts">
// The HSPs of the selected subject in the selected query (docs/web/results_columns.md
// "HSP 一覧"): the fields of each HSP's outfmt 6 row as written, the frames of the record,
// and the orientation that the record's coordinates and frames give.
import { computed } from 'vue';
import type { HspEntry, ResultsBrowser, ResultsState } from '../application/results';
import type { HspSortKey } from '../domain/result-index';
import SortButton from './SortButton.vue';
import VirtualRows from './VirtualRows.vue';

const props = defineProps<{ results: ResultsBrowser; state: ResultsState }>();
const ROW_PX = 32;
const selectedPosition = computed(() =>
  props.state.hsps.findIndex((hsp) => hsp.id.qIdx === props.state.hsp?.qIdx && hsp.id.rank === props.state.hsp?.rank),
);
const units = computed(() => props.state.loaded?.units ?? { query: '', subject: '' });
const framed = computed(() => props.state.hsps.some((hsp) => hsp.queryFrame !== undefined || hsp.subjectFrame !== undefined));
const subject = computed(() => props.state.subjects.find((s) => s.sIdx === props.state.sIdx));

function sortBy(key: HspSortKey, descending: boolean): void {
  props.results.setHspSort({ key, descending });
}

const signed = (value: number | undefined) => (value === undefined ? '' : value > 0 ? `+${value}` : String(value));
const ORIENTATION: Readonly<Record<HspEntry['orientation'], string>> = {
  forward: 'Forward',
  reverse: 'Reverse',
  unknown: 'Not in the record',
};
</script>

<template>
  <div class="result-table hsp-table" data-testid="hsp-table">
    <h3>
      HSPs of {{ subject?.first.sseqid }}
      <span class="muted small">{{ state.hsps.length }} shown</span>
    </h3>
    <div class="table-scroll">
      <div class="table-head hsp-grid" :class="{ framed }" role="row">
        <SortButton label="#" sort-key="rank" :sort="state.hspSort" scope="hsp" @sort="sortBy" />
        <SortButton label="Bit score" sort-key="bitScore" :sort="state.hspSort" scope="hsp" first-descending @sort="sortBy" />
        <SortButton label="E value" sort-key="eValue" :sort="state.hspSort" scope="hsp" @sort="sortBy" />
        <span role="columnheader">Identity (%)</span>
        <span role="columnheader">Length</span>
        <span role="columnheader">Mismatches</span>
        <span role="columnheader">Gap opens</span>
        <SortButton :label="`Query (${units.query})`" sort-key="qStart" :sort="state.hspSort" scope="hsp" @sort="sortBy" />
        <SortButton :label="`Subject (${units.subject})`" sort-key="sStart" :sort="state.hspSort" scope="hsp" @sort="sortBy" />
        <span v-if="framed" role="columnheader">Frames (q/s)</span>
        <span role="columnheader">Orientation</span>
        <span role="columnheader">outfmt 0</span>
      </div>
      <VirtualRows
        :count="state.hsps.length"
        :row-px="ROW_PX"
        :max-rows="8"
        :reveal="selectedPosition"
        :reveal-key="`${state.runId}:${state.hsp?.qIdx}:${state.hsp?.rank}`"
        :order-key="`${state.hspSort.key}:${state.hspSort.descending}`"
        label="HSPs"
        testid="hsp-list"
      >
        <template #row="{ position }">
          <button
            type="button"
            class="table-row hsp-grid"
            :class="{ framed, selected: position === selectedPosition }"
            :aria-pressed="position === selectedPosition"
            :data-testid="`hsp-row-${state.hsps[position]!.id.qIdx}-${state.hsps[position]!.id.rank}`"
            :data-orientation="state.hsps[position]!.orientation"
            @click="results.selectHsp(state.hsps[position]!.id)"
          >
            <span class="num">{{ state.hsps[position]!.id.rank + 1 }}</span>
            <span class="num" data-field="bitscore">{{ state.hsps[position]!.fields.bitscore }}</span>
            <span class="num" data-field="evalue">{{ state.hsps[position]!.fields.evalue }}</span>
            <span class="num" data-field="pident">{{ state.hsps[position]!.fields.pident }}</span>
            <span class="num" data-field="length">{{ state.hsps[position]!.fields.length }}</span>
            <span class="num" data-field="mismatch">{{ state.hsps[position]!.fields.mismatch }}</span>
            <span class="num" data-field="gapopen">{{ state.hsps[position]!.fields.gapopen }}</span>
            <span class="num" data-field="query"
              >{{ state.hsps[position]!.fields.qstart }}–{{ state.hsps[position]!.fields.qend }}</span
            >
            <span class="num" data-field="subject"
              >{{ state.hsps[position]!.fields.sstart }}–{{ state.hsps[position]!.fields.send }}</span
            >
            <span v-if="framed" class="num" data-field="frames"
              >{{ signed(state.hsps[position]!.queryFrame) || '–' }}/{{ signed(state.hsps[position]!.subjectFrame) || '–' }}</span
            >
            <span class="orientation" :data-orientation="state.hsps[position]!.orientation" data-field="orientation">{{
              ORIENTATION[state.hsps[position]!.orientation]
            }}</span>
            <span data-field="outfmt0">{{ state.hsps[position]!.inOutfmt0 ? 'shown' : 'not shown' }}</span>
          </button>
        </template>
      </VirtualRows>
    </div>
  </div>
</template>
