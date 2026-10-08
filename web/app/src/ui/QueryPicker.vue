<script setup lang="ts">
// The queries of the run (its query records, in order), with their subject and HSP counts.
// A run of many queries is drawn a screenful at a time.
import { computed } from 'vue';
import type { ResultsBrowser, ResultsState } from '../application/results';
import { formatCount } from './format';
import VirtualRows from './VirtualRows.vue';

const props = defineProps<{ results: ResultsBrowser; state: ResultsState }>();
const ROW_PX = 30;
const selectedPosition = computed(() => props.state.queries.findIndex((query) => query.qIdx === props.state.qIdx));
const unit = computed(() => props.state.loaded?.units.query ?? '');
const total = computed(() => props.state.loaded?.run.snapshot.query.records.length ?? 0);

function setText(event: Event): void {
  const queryText = (event.target as HTMLInputElement).value;
  props.results.setFilters({ ...props.state.filters, queryText });
}

function setHitsOnly(event: Event): void {
  props.results.setFilters({ ...props.state.filters, queriesWithHitsOnly: (event.target as HTMLInputElement).checked });
}
</script>

<template>
  <div class="query-picker">
    <h3>Query</h3>
    <div class="record-tools">
      <label>
        <span class="visually-hidden">Find queries by ID</span>
        <input
          type="search"
          placeholder="Find by ID"
          :value="state.filters.queryText ?? ''"
          data-testid="query-filter"
          @input="setText"
        />
      </label>
      <label class="inline">
        <input type="checkbox" :checked="state.filters.queriesWithHitsOnly ?? false" data-testid="filter-hits-only" @change="setHitsOnly" />
        With hits only
      </label>
    </div>
    <p class="muted small" data-testid="query-count">{{ formatCount(state.queries.length) }} of {{ formatCount(total) }} queries</p>
    <VirtualRows :count="state.queries.length" :row-px="ROW_PX" :max-rows="8" :reveal="selectedPosition" label="Queries" testid="query-list">
      <template #row="{ position }">
        <button
          type="button"
          class="pick-row"
          :class="{ selected: state.queries[position]!.qIdx === state.qIdx, empty: state.queries[position]!.hsps === 0 }"
          :aria-pressed="state.queries[position]!.qIdx === state.qIdx"
          :data-testid="`query-row-${state.queries[position]!.qIdx}`"
          :title="state.queries[position]!.id"
          @click="results.selectQuery(state.queries[position]!.qIdx)"
        >
          <span class="record-number">#{{ state.queries[position]!.qIdx + 1 }}</span>
          <span class="record-id">{{ state.queries[position]!.id }}</span>
          <span class="pick-meta">
            {{ formatCount(state.queries[position]!.length) }} {{ unit }} ·
            <template v-if="state.queries[position]!.hsps === 0">no hits</template>
            <template v-else>{{ formatCount(state.queries[position]!.subjects) }} subj., {{ formatCount(state.queries[position]!.hsps) }} HSPs</template>
          </span>
        </button>
      </template>
    </VirtualRows>
  </div>
</template>
