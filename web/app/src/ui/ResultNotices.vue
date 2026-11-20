<script setup lang="ts">
// What the lists of the selected query do not show, and why (design §11.2): a dry run, no hits, HSPs
// hidden by the view filters, a hit list that reached its limit (more subjects may match;
// it is not said that they do), and subjects whose alignments outfmt 0 does not show.
import { computed } from 'vue';
import { noHspReason, type ResultsState } from '../application/results';
import { flagGiven } from '../domain/result-index';
import { formatCount } from './format';
import { NO_HSP_TEXT } from './no-hsp-text';
import { placeParts } from './message-place-text';
import OptionName from './OptionName.vue';
import PlaceWords from './PlaceWords.vue';
import { vScrollFocus } from './scroll-focus';

const props = defineProps<{ state: ResultsState }>();
defineEmits<{ clear: [] }>();

const totals = computed(() => props.state.queryTotals);
const limit = computed(() => props.state.loaded?.limits.maxTargetSeqs);
/** The option that sets the hit list of the tabular formats, which the HSP records and these lists follow. */
const limitFlag = computed(() => props.state.loaded?.limits.subjectFlag ?? '-max_target_seqs');
const limitGiven = computed(() => flagGiven(props.state.loaded?.run.snapshot.argv ?? [], limitFlag.value));
/** What outfmt 0 does with -num_descriptions / -num_alignments; undefined where neither is given. */
const outfmt0 = computed(() => props.state.loaded?.limits.outfmt0);
/** -sorthits orders outfmt 0's hit list: its first subjects are not the first ones in the engine's order. */
const sorted = computed(() => props.state.loaded?.limits.sorted === true);
/** Where the error of outfmt 0 is, in words, when the run could not write it, as the engine said; undefined where it was written. */
const outfmt0Error = computed(() => {
  const loaded = props.state.loaded;
  return loaded?.outfmt0Error === undefined ? undefined : placeParts(loaded.outfmt0Error, loaded.run);
});
/** Where the error of outfmt 6 and 7 is, in words, when the run could not write them, as the engine said; undefined where they were written. */
const tabularError = computed(() => {
  const loaded = props.state.loaded;
  return loaded?.tabularError === undefined ? undefined : placeParts(loaded.tabularError, loaded.run);
});
/** Why the selected query shows no HSP (the Dot Plot tab says the same in place of its plot). */
const reason = computed(() => noHspReason(props.state));
/** What the engine wrote to the warnings (BLAST+'s "Query is Empty!"), as it wrote it. */
const engineWarnings = computed(() => props.state.loaded?.diagnostics.trim() ?? '');
const plural = (count: number, one: string, many = `${one}s`) => `${formatCount(count)} ${count === 1 ? one : many}`;
</script>

<template>
  <div class="result-notices" aria-live="polite">
    <!-- A run that printed the command line's help or version searched nothing (line WQA's D14, 2026-11-20). -->
    <p v-if="reason === 'printed-help' || reason === 'printed-version'" class="notice" data-testid="results-notice" :data-kind="reason">
      {{ NO_HSP_TEXT[reason] }}
    </p>
    <p v-else-if="reason === 'dry-run'" class="notice" data-testid="results-notice" data-kind="dry-run">{{ NO_HSP_TEXT['dry-run'] }}</p>
    <div v-else-if="reason === 'empty-query'" class="notice" data-testid="results-notice" data-kind="empty-query">
      <p>{{ NO_HSP_TEXT['empty-query'] }}</p>
      <pre v-if="engineWarnings !== ''" v-scroll-focus="'Warnings'" class="output" data-testid="results-notice-engine">{{ engineWarnings }}</pre>
    </div>
    <p v-else-if="reason === 'no-tabular-hits'" class="notice" data-testid="results-notice" data-kind="no-tabular-hits">
      <OptionName name="-num_alignments" /> 0 leaves the hit list of outfmt 6 and 7 empty (the "No hits are being saved"
      error<template v-if="tabularError !== undefined">, <PlaceWords :parts="tabularError" /></template>), so this run has no HSP records; outfmt 0 is still
      written.
    </p>
    <p v-else-if="reason === 'no-tabular-output'" class="notice" data-testid="results-notice" data-kind="no-tabular-output">
      {{ NO_HSP_TEXT['no-tabular-output'] }} Their error is <PlaceWords :parts="tabularError ?? []" />.
    </p>
    <p v-else-if="reason === 'no-queries'" class="notice" data-testid="results-notice" data-kind="no-queries">
      {{ NO_HSP_TEXT['no-queries'] }}
    </p>
    <template v-else-if="totals">
      <p v-if="reason === 'no-hits'" class="notice" data-testid="results-notice" data-kind="no-hits">
        {{ NO_HSP_TEXT['no-hits'] }}
      </p>
      <p v-else-if="reason === 'filtered-out'" class="notice warn" data-testid="results-notice" data-kind="filtered-out">
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
        ({{ limitGiven ? '' : 'the default of ' }}<OptionName :name="limitFlag" />: {{ formatCount(limit) }}). More
        subjects may match; a search with a larger <OptionName :name="limitFlag" /> would show them.
      </p>
      <p v-if="totals.subjects > 0 && outfmt0Error !== undefined" class="notice" data-testid="results-notice" data-kind="outfmt0-failed">
        outfmt 0 could not be written (its error is <PlaceWords :parts="outfmt0Error" />), so no alignment is linked to this query's
        {{ plural(totals.subjects, 'subject') }}. outfmt 6 and 7 hold the HSPs.
      </p>
      <p
        v-else-if="totals.subjects > 0 && state.outfmt0Subjects < totals.subjects && outfmt0 === undefined"
        class="notice"
        data-testid="results-notice"
        data-kind="outfmt0-partial"
      >
        outfmt 0 shows the alignments of the first {{ plural(state.outfmt0Subjects, 'subject') }} of this query<template v-if="sorted">
          in its order of <OptionName name="-sorthits" /></template
        >
        (<OptionName name="-num_alignments" />: 250, or <OptionName name="-max_target_seqs" /> when it is
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
            (<OptionName name="-num_alignments" />: {{ formatCount(outfmt0.alignments) }})</template
          >. Its search has a hit list of its own<template v-if="outfmt0.hitList !== undefined">
            ({{ plural(outfmt0.hitList, 'subject') }} per query)</template
          >, because <OptionName name="-num_descriptions" /> and <OptionName name="-num_alignments" /> give it another
          size than outfmt 6 and 7, so a subject's HSPs are linked only where outfmt 0 shows the alignments of exactly the HSPs that outfmt 6 and 7 keep for it.
        </template>
        <template v-else
          >outfmt 0 shows the alignments of {{ plural(state.outfmt0Subjects, 'subject') }} of this query's {{ formatCount(totals.subjects) }}<template
            v-if="outfmt0.alignments !== undefined"
          >
            (<OptionName name="-num_alignments" />: {{ formatCount(outfmt0.alignments) }})</template
          >.
        </template>
        The other subjects' HSPs are in outfmt 6 and 7.
      </p>
    </template>
  </div>
</template>
