<script setup lang="ts">
// The HSPs of the selected subject in the selected query (docs/web/results_columns.md
// "HSP 一覧"): the fields of each HSP's outfmt 6 row as written, the frames of the record,
// and the orientation that the record's coordinates and frames give.
import { computed, ref } from 'vue';
import type { HspEntry, HspId, ResultsBrowser, ResultsState } from '../application/results';
import type { HspSortKey } from '../domain/result-index';
import { formatCount, framesLabel, framesText } from './format';
import SortButton from './SortButton.vue';
import { focusPressed } from './useSideScroll';
import VirtualRows from './VirtualRows.vue';

const props = defineProps<{ results: ResultsBrowser; state: ResultsState }>();
/** An HSP chosen in the table (after it is selected): the Alignments bring its Range into view. */
const emit = defineEmits<{ chosen: [id: HspId] }>();
const ROW_PX = 28;
const selectedPosition = computed(() =>
  props.state.hsps.findIndex((hsp) => hsp.id.qIdx === props.state.hsp?.qIdx && hsp.id.rank === props.state.hsp?.rank),
);
const units = computed(() => props.state.loaded?.units ?? { query: '', subject: '' });
const framed = computed(() => props.state.hsps.some((hsp) => hsp.queryFrame !== undefined || hsp.subjectFrame !== undefined));
/** "Subject frame" (TBLASTN), "Query frame" (BLASTX) or "Frames (q/s)" (TBLASTX): the HSPs of a run have the same frames. */
const framesHead = computed(() => (props.state.hsps[0] === undefined ? '' : (framesLabel(props.state.hsps[0]) ?? '')));
const subject = computed(() => props.state.subjects.find((s) => s.sIdx === props.state.sIdx));
const scroller = ref<HTMLElement>();
const gutter = ref(0);

function sortBy(key: HspSortKey, descending: boolean): void {
  props.results.setHspSort({ key, descending });
}

function choose(id: HspId): void {
  props.results.selectHsp(id);
  emit('chosen', id);
}

/** The texts of the cells that join two values of the record. */
const range = (start: string | undefined, end: string | undefined) => `${start}–${end}`;
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
      <span class="muted small" data-testid="hsp-count">{{ formatCount(state.hsps.length) }} shown</span>
    </h3>
    <div
      ref="scroller"
      class="table-scroll"
      role="grid"
      aria-label="HSPs"
      :aria-rowcount="state.hsps.length + 1"
      :style="{ '--row-gutter': `${gutter}px` }"
    >
      <div class="table-head hsp-grid" :class="{ framed }" role="row" aria-rowindex="1">
        <SortButton class="num-head" data-col="rank" label="#" sort-key="rank" :sort="state.hspSort" scope="hsp" @sort="sortBy" />
        <SortButton
          class="num-head"
          data-col="bitscore"
          label="Bit score"
          sort-key="bitScore"
          :sort="state.hspSort"
          scope="hsp"
          first-descending
          @sort="sortBy"
        />
        <SortButton class="num-head" data-col="evalue" label="E value" sort-key="eValue" :sort="state.hspSort" scope="hsp" @sort="sortBy" />
        <span class="num-head" data-col="pident" role="columnheader">Identity (%)</span>
        <span class="num-head" data-col="length" role="columnheader">Length</span>
        <span class="num-head" data-col="mismatch" role="columnheader">Mismatches</span>
        <span class="num-head" data-col="gapopen" role="columnheader">Gap opens</span>
        <SortButton
          class="num-head"
          data-col="query"
          :label="`Query (${units.query})`"
          sort-key="qStart"
          :sort="state.hspSort"
          scope="hsp"
          @sort="sortBy"
        />
        <SortButton
          class="num-head"
          data-col="subject"
          :label="`Subject (${units.subject})`"
          sort-key="sStart"
          :sort="state.hspSort"
          scope="hsp"
          @sort="sortBy"
        />
        <span v-if="framed" class="num-head" data-col="frames" role="columnheader">{{ framesHead }}</span>
        <span data-col="orientation" role="columnheader">Orientation</span>
        <span data-col="outfmt0" role="columnheader">outfmt 0</span>
      </div>
      <VirtualRows
        :count="state.hsps.length"
        :row-px="ROW_PX"
        :max-rows="8"
        :reveal="selectedPosition"
        :reveal-key="`${state.runId}:${state.hsp?.qIdx}:${state.hsp?.rank}`"
        :order-key="`${state.hspSort.key}:${state.hspSort.descending}`"
        kind="grid"
        testid="hsp-list"
        @gutter="gutter = $event"
      >
        <template #row="{ position, tabindex }">
          <div
            v-for="hsp in [state.hsps[position]!]"
            :key="`${hsp.id.qIdx}:${hsp.id.rank}`"
            role="row"
            :aria-rowindex="position + 2"
            :aria-selected="position === selectedPosition"
            :tabindex="tabindex"
            class="table-row hsp-grid"
            :class="{ framed, selected: position === selectedPosition }"
            :data-testid="`hsp-row-${hsp.id.qIdx}-${hsp.id.rank}`"
            :data-orientation="hsp.orientation"
            @mousedown="focusPressed"
            @click="choose(hsp.id)"
          >
            <span role="gridcell" class="num" data-field="rank" :title="String(hsp.id.rank + 1)">{{ hsp.id.rank + 1 }}</span>
            <span role="gridcell" class="num" data-field="bitscore" :title="hsp.fields.bitscore">{{ hsp.fields.bitscore }}</span>
            <span role="gridcell" class="num" data-field="evalue" :title="hsp.fields.evalue">{{ hsp.fields.evalue }}</span>
            <span role="gridcell" class="num" data-field="pident" :title="hsp.fields.pident">{{ hsp.fields.pident }}</span>
            <span role="gridcell" class="num" data-field="length" :title="hsp.fields.length">{{ hsp.fields.length }}</span>
            <span role="gridcell" class="num" data-field="mismatch" :title="hsp.fields.mismatch">{{ hsp.fields.mismatch }}</span>
            <span role="gridcell" class="num" data-field="gapopen" :title="hsp.fields.gapopen">{{ hsp.fields.gapopen }}</span>
            <span role="gridcell" class="num" data-field="query" :title="range(hsp.fields.qstart, hsp.fields.qend)">{{
              range(hsp.fields.qstart, hsp.fields.qend)
            }}</span>
            <span role="gridcell" class="num" data-field="subject" :title="range(hsp.fields.sstart, hsp.fields.send)">{{
              range(hsp.fields.sstart, hsp.fields.send)
            }}</span>
            <span v-if="framed" role="gridcell" class="num" data-field="frames" :title="framesText(hsp)">{{ framesText(hsp) }}</span>
            <span role="gridcell" class="orientation" :data-orientation="hsp.orientation" data-field="orientation" :title="ORIENTATION[hsp.orientation]">{{
              ORIENTATION[hsp.orientation]
            }}</span>
            <span role="gridcell" data-field="outfmt0" :title="hsp.inOutfmt0 ? 'shown' : 'not shown'">{{ hsp.inOutfmt0 ? 'shown' : 'not shown' }}</span>
          </div>
        </template>
      </VirtualRows>
    </div>
  </div>
</template>
