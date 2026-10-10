<script setup lang="ts">
// The view filters (ViewState, design §11.2): they change what the lists and the dot plot
// show, never the search, and never the compatibility outputs. They compare the engine's
// values of the HSP records (docs/web/results_columns.md).
import { computed, ref, watch } from 'vue';
import type { ResultsBrowser } from '../application/results';
import type { ViewFilters } from '../domain/result-index';
import { useNarrow } from './useNarrow';

const props = defineProps<{ results: ResultsBrowser; filters: ViewFilters }>();

const evalue = ref('');
const bits = ref('');
const subject = ref('');
const problem = ref('');

// On a phone the form is a disclosure, closed at first, so that the graphic comes sooner in the
// first screen after "Open results" (screen review 2 L1); closed, it says how many filters are set.
const narrow = useNarrow();
const expanded = ref(false);
const formId = 'result-filters-fields';
const active = computed(() => [props.filters.maxEValue, props.filters.minBitScore, props.filters.subjectText].filter((value) => value !== undefined).length);
const open = computed(() => !narrow.value || expanded.value);

// The boxes show the filters in force, also when this form is drawn again for the same run (the
// results tab shown again, the run opened again). A number keeps the text typed for it.
watch(
  () => props.filters,
  (filters) => {
    if (typed(evalue.value) !== filters.maxEValue) evalue.value = filters.maxEValue === undefined ? '' : String(filters.maxEValue);
    if (typed(bits.value) !== filters.minBitScore) bits.value = filters.minBitScore === undefined ? '' : String(filters.minBitScore);
    subject.value = filters.subjectText ?? '';
  },
  { immediate: true },
);

/** The number that a box's text gives, if any. */
function typed(text: string): number | undefined {
  const trimmed = text.trim();
  const value = Number(trimmed);
  return trimmed === '' || !Number.isFinite(value) ? undefined : value;
}

/** A number typed in a filter; empty means no filter. */
function number(text: string, name: string): number | undefined | null {
  const trimmed = text.trim();
  if (trimmed === '') return undefined;
  const value = Number(trimmed);
  if (!Number.isFinite(value)) {
    problem.value = `${name} must be a number, such as 1e-5 or 50.`;
    return null;
  }
  return value;
}

function apply(): void {
  problem.value = '';
  const maxEValue = number(evalue.value, 'E value');
  const minBitScore = number(bits.value, 'Bit score');
  if (maxEValue === null || minBitScore === null) return;
  const text = subject.value.trim();
  const next: ViewFilters = {
    ...(props.filters.queriesWithHitsOnly === undefined ? {} : { queriesWithHitsOnly: props.filters.queriesWithHitsOnly }),
    ...(props.filters.queryText === undefined ? {} : { queryText: props.filters.queryText }),
    ...(maxEValue === undefined ? {} : { maxEValue }),
    ...(minBitScore === undefined ? {} : { minBitScore }),
    ...(text === '' ? {} : { subjectText: text }),
  };
  props.results.setFilters(next);
}

function clear(): void {
  evalue.value = '';
  bits.value = '';
  subject.value = '';
  apply();
}
</script>

<template>
  <form class="result-filters" aria-label="Filter Results" data-testid="view-filters" @submit.prevent="apply">
    <h3>
      <button
        v-if="narrow"
        type="button"
        class="filter-toggle"
        :aria-expanded="expanded"
        :aria-controls="formId"
        data-testid="filter-toggle"
        @click="expanded = !expanded"
      >
        <span aria-hidden="true">{{ expanded ? '▾' : '▸' }}</span> Filter Results<span v-if="!expanded && active > 0" class="muted small" data-testid="filter-active">
          ({{ active }} set)</span
        >
      </button>
      <template v-else>Filter Results</template>
    </h3>
    <div v-show="open" :id="formId" class="filter-body">
      <div class="filter-fields">
        <label>
          E value ≤
          <input v-model="evalue" type="text" inputmode="decimal" placeholder="any" data-testid="filter-evalue" @change="apply" />
        </label>
        <label>
          Bit score ≥
          <input v-model="bits" type="text" inputmode="decimal" placeholder="any" data-testid="filter-bits" @change="apply" />
        </label>
        <label>
          Subject ID contains
          <input v-model="subject" type="search" placeholder="any" data-testid="filter-subject" @change="apply" />
        </label>
      </div>
      <div class="filter-actions">
        <button type="submit" class="filter-submit" data-testid="filter-apply">Filter</button>
        <button type="button" data-testid="filter-clear" @click="clear">Reset</button>
      </div>
      <p class="muted small">Filters change the view only. Exports of outfmt 0, 6 and 7 always hold the whole result.</p>
      <p v-if="problem" class="error" role="alert" data-testid="filter-problem">{{ problem }}</p>
    </div>
  </form>
</template>
