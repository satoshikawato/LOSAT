<script setup lang="ts">
// What the lists of the selected query do not show, and why (design §11.2): no hits, HSPs
// hidden by the view filters, a hit list that reached its limit (more subjects may match;
// it is not said that they do), and subjects whose alignments outfmt 0 does not show.
import { computed } from 'vue';
import type { ResultsState } from '../application/results';
import { flagGiven } from '../domain/result-index';
import { formatCount } from './format';

const props = defineProps<{ state: ResultsState }>();
defineEmits<{ clear: [] }>();

const totals = computed(() => props.state.queryTotals);
const limit = computed(() => props.state.loaded?.limits.maxTargetSeqs);
/** The option that sets the hit list of the tabular formats, which the HSP records and these lists follow. */
const limitFlag = computed(() => props.state.loaded?.limits.subjectFlag ?? '-max_target_seqs');
const limitGiven = computed(() => flagGiven(props.state.loaded?.run.snapshot.argv ?? [], limitFlag.value));
/** What outfmt 0 does with -num_descriptions / -num_alignments; undefined where neither is given. */
const outfmt0 = computed(() => props.state.loaded?.limits.outfmt0);
const plural = (count: number, one: string, many = `${one}s`) => `${formatCount(count)} ${count === 1 ? one : many}`;
</script>

<template>
  <div class="result-notices" aria-live="polite">
    <p v-if="limit === 0" class="notice" data-testid="results-notice" data-kind="no-tabular-hits">
      <span class="option-name">-num_alignments</span> 0 leaves the hit list of outfmt 6 and 7 empty (BLAST+'s "No hits are being saved" error,
      in the diagnostics), so this run has no HSP records; outfmt 0 is written as BLAST+ writes it.
    </p>
    <p v-else-if="state.queries.length === 0" class="notice" data-testid="results-notice" data-kind="no-queries">
      No queries match the query filters.
    </p>
    <template v-else-if="totals">
      <p v-if="totals.hsps === 0" class="notice" data-testid="results-notice" data-kind="no-hits">
        No significant similarity found for this query.
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
        ({{ limitGiven ? '' : 'the default of ' }}<span class="option-name">{{ limitFlag }}</span>: {{ formatCount(limit) }}). More
        subjects may match; a search with a larger <span class="option-name">{{ limitFlag }}</span> would show them.
      </p>
      <p
        v-if="totals.subjects > 0 && state.outfmt0Subjects < totals.subjects && outfmt0 === undefined"
        class="notice"
        data-testid="results-notice"
        data-kind="outfmt0-partial"
      >
        outfmt 0 shows the alignments of the first {{ plural(state.outfmt0Subjects, 'subject') }} of this query
        (BLAST+'s <span class="option-name">-num_alignments</span>: 250, or <span class="option-name">-max_target_seqs</span> when it is
        given). The other subjects' HSPs are in outfmt 6 and 7.
      </p>
      <p
        v-else-if="totals.subjects > 0 && state.outfmt0Subjects < totals.subjects && outfmt0 !== undefined"
        class="notice"
        data-testid="results-notice"
        data-kind="outfmt0-partial"
      >
        <template v-if="outfmt0.separate"
          >outfmt 0's alignments are linked to {{ formatCount(state.outfmt0Subjects) }} of this query's {{ plural(totals.subjects, 'subject') }}<template
            v-if="outfmt0.alignments !== undefined"
          >
            (BLAST+'s <span class="option-name">-num_alignments</span>: {{ formatCount(outfmt0.alignments) }})</template
          >. Its search has a hit list of its own<template v-if="outfmt0.hitList !== undefined">
            ({{ plural(outfmt0.hitList, 'subject') }} per query)</template
          >, because <span class="option-name">-num_descriptions</span> and <span class="option-name">-num_alignments</span> give it another
          size than outfmt 6 and 7, so a subject's HSPs are linked only where outfmt 0 shows the alignments of exactly the HSPs that outfmt 6 and 7 keep for it.
        </template>
        <template v-else
          >outfmt 0 shows the alignments of {{ plural(state.outfmt0Subjects, 'subject') }} of this query's {{ formatCount(totals.subjects) }}<template
            v-if="outfmt0.alignments !== undefined"
          >
            (BLAST+'s <span class="option-name">-num_alignments</span>: {{ formatCount(outfmt0.alignments) }})</template
          >.
        </template>
        The other subjects' HSPs are in outfmt 6 and 7.
      </p>
    </template>
  </div>
</template>
