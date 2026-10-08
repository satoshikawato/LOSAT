<script setup lang="ts">
// What the lists of the selected query do not show, and why (design §11.2): no hits, HSPs
// hidden by the view filters, a hit list that reached its limit (more subjects may match;
// it is not said that they do), and subjects whose alignments outfmt 0 does not show.
import { computed } from 'vue';
import type { ResultsState } from '../application/results';
import { formatCount } from './format';

const props = defineProps<{ state: ResultsState }>();
defineEmits<{ clear: [] }>();

const totals = computed(() => props.state.queryTotals);
const limit = computed(() => props.state.loaded?.limits.maxTargetSeqs);
const limitGiven = computed(() => props.state.loaded?.run.snapshot.argv.includes('-max_target_seqs') ?? false);
const plural = (count: number, one: string, many = `${one}s`) => `${formatCount(count)} ${count === 1 ? one : many}`;
</script>

<template>
  <div class="result-notices" aria-live="polite">
    <p v-if="state.queries.length === 0" class="notice" data-testid="results-notice" data-kind="no-queries">
      No queries match the query filters.
    </p>
    <template v-else-if="totals">
      <p v-if="totals.hsps === 0" class="notice" data-testid="results-notice" data-kind="no-hits">
        No hits for this query.
      </p>
      <p v-else-if="state.subjects.length === 0" class="notice warn" data-testid="results-notice" data-kind="filtered-out">
        No HSPs of this query match the view filters ({{ plural(state.hidden.hsps, 'HSP') }} on
        {{ plural(state.hidden.subjects, 'subject') }} hidden).
        <button type="button" class="link" data-testid="results-notice-clear" @click="$emit('clear')">Clear the filters</button>
      </p>
      <p v-else-if="state.hidden.hsps > 0" class="notice" data-testid="results-notice" data-kind="filtered-some">
        The view filters hide {{ plural(state.hidden.hsps, 'HSP') }}
        <template v-if="state.hidden.subjects > 0">and {{ plural(state.hidden.subjects, 'subject') }}</template> of this query.
      </p>
      <p v-if="state.atSubjectLimit && limit !== undefined" class="notice" data-testid="results-notice" data-kind="subject-limit">
        This query has {{ plural(totals.subjects, 'subject') }}, the most that the search keeps
        ({{ limitGiven ? '-max_target_seqs' : 'the default of -max_target_seqs' }}: {{ formatCount(limit) }}). More subjects may
        match; a search with a larger -max_target_seqs would show them.
      </p>
      <p
        v-if="totals.subjects > 0 && state.outfmt0Subjects < totals.subjects"
        class="notice"
        data-testid="results-notice"
        data-kind="outfmt0-partial"
      >
        outfmt 0 shows the alignments of the first {{ plural(state.outfmt0Subjects, 'subject') }} of this query
        (BLAST+'s -num_alignments: 250, or -max_target_seqs when it is given). The other subjects' HSPs are in outfmt 6 and 7.
      </p>
    </template>
  </div>
</template>
