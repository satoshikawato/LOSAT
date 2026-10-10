<script setup lang="ts">
// The first view of the results, after NCBI's classic (Traditional) results page (the Owner,
// 2026-10-10; docs/web/ncbi_ui_mapping.md §2): on one page, the selected query's Graphic Summary,
// its Descriptions directly under it, then the selected subject's Alignments, each under NCBI's
// heading. As NCBI's page, it is cut: the graphic draws the first 100 subjects and the
// Descriptions list the first 100 ("Show all N" for the rest), and the Alignments are the selected
// subject's only.
//
// A Description row or a bar of the graphic chosen brings the alignments into view below with the
// focus there (the subject's heading, or the HSP's Range); the Alignments' "Descriptions" goes
// back up to the list. Only choosing moves the page: the graphic's arrow keys move within the
// figure, and a sort, a filter or another query keeps the page where it is.
import { nextTick, ref, useId } from 'vue';
import type { HspId, ResultsBrowser, ResultsState } from '../application/results';
import AlignmentsView from './AlignmentsView.vue';
import GraphicSummary from './GraphicSummary.vue';
import SubjectTable from './SubjectTable.vue';

defineProps<{
  results: ResultsBrowser;
  state: ResultsState;
  /** The keys of the HSPs in the candidate tray. */
  inTray: ReadonlySet<string>;
}>();
const emit = defineEmits<{ 'add-candidates': [ids: readonly HspId[]] }>();

const ids = { graphic: useId(), descriptions: useId(), alignments: useId() };
const descriptions = ref<HTMLElement>();
const descriptionsHeading = ref<HTMLElement>();
const alignmentsHeading = ref<HTMLElement>();
const alignments = ref<InstanceType<typeof AlignmentsView>>();

/**
 * Brings a subject's alignments into view (a Description row chosen): the Alignments' heading at
 * the top of the window, with the focus on it, so that the keyboard goes on from what is shown.
 * The scroll is instant, as the results' heading's (ResultsPanel.vue).
 */
async function showSubject(): Promise<void> {
  await nextTick();
  alignmentsHeading.value?.scrollIntoView({ block: 'start' });
  alignmentsHeading.value?.focus({ preventScroll: true });
}

/** Brings an HSP's Range into view (the graphic's "click to show alignments", "Show in results"); with `focus`, the focus too. */
async function showHsp(id: HspId, focus = false): Promise<void> {
  await nextTick();
  await alignments.value?.reveal(id, focus);
}
defineExpose({ showHsp });

/** The Alignments' "Descriptions": the list in view again, with the focus on the selected subject's row where it is drawn. */
async function showDescriptions(): Promise<void> {
  await nextTick();
  descriptionsHeading.value?.scrollIntoView({ block: 'start' });
  const row = descriptions.value?.querySelector<HTMLElement>('[data-testid^="subject-row-"][aria-pressed="true"]');
  (row ?? descriptionsHeading.value)?.focus({ preventScroll: true });
}
</script>

<template>
  <div class="classic-results" data-testid="results-classic">
    <template v-if="state.subjects.length > 0">
      <section class="results-part" :aria-labelledby="ids.graphic" data-testid="results-graphic">
        <h3 :id="ids.graphic" class="part-heading">Graphic Summary</h3>
        <GraphicSummary :results="results" :state="state" @show-alignment="showHsp($event, true)" />
      </section>
      <section ref="descriptions" class="results-part" :aria-labelledby="ids.descriptions" data-testid="results-descriptions">
        <h3 :id="ids.descriptions" ref="descriptionsHeading" class="part-heading" tabindex="-1" data-testid="results-descriptions-heading">
          Descriptions
        </h3>
        <SubjectTable :results="results" :state="state" @chosen="showSubject" @add-candidates="emit('add-candidates', $event)" />
      </section>
    </template>
    <section v-if="state.hsps.length > 0" class="results-part" :aria-labelledby="ids.alignments" data-testid="results-alignments">
      <h3 :id="ids.alignments" ref="alignmentsHeading" class="part-heading" tabindex="-1" data-testid="results-alignments-heading">Alignments</h3>
      <AlignmentsView
        ref="alignments"
        :results="results"
        :state="state"
        :in-tray="inTray"
        @descriptions="showDescriptions"
        @add-candidates="emit('add-candidates', $event)"
      />
    </section>
  </div>
</template>
