<script setup lang="ts">
// The selected HSP as the engine wrote it (plan §4.5): its subject's heading and its
// section of outfmt 0 (score lines and alignment), and its outfmt 6 row, byte for byte.
// Nothing of the alignment is drawn again here (midlines and masks are the formatter's).
import { computed } from 'vue';
import type { ResultsState } from '../application/results';

const props = defineProps<{ state: ResultsState }>();
const detail = computed(() => props.state.detail);
const entry = computed(() =>
  props.state.hsps.find((hsp) => hsp.id.qIdx === detail.value?.id.qIdx && hsp.id.rank === detail.value?.id.rank),
);
const subject = computed(() => props.state.subjects.find((s) => s.sIdx === props.state.sIdx));
const query = computed(() => props.state.queries.find((q) => q.qIdx === props.state.qIdx));
const limitGiven = computed(() => props.state.loaded?.run.snapshot.argv.includes('-max_target_seqs') ?? false);
</script>

<template>
  <div
    v-if="detail"
    class="hsp-detail"
    data-testid="hsp-detail"
    :data-hsp="`${detail.id.qIdx}:${detail.id.rank}`"
    :data-state="detail.state"
  >
    <h3>
      HSP {{ detail.id.rank + 1 }} of query {{ query?.id ?? `#${detail.id.qIdx + 1}` }}
      <span class="muted small">subject {{ subject?.first.sseqid }}</span>
    </h3>
    <p v-if="entry?.orientation === 'unknown'" class="notice" data-testid="detail-strand-note">
      This HSP covers one letter of each sequence, so its coordinates do not show its strand, and the HSP record does not hold
      it. The <code>Strand=</code> line of its outfmt 0 section below shows it.
    </p>
    <p v-if="detail.state === 'failed'" class="error" data-testid="detail-error">The HSP could not be read: {{ detail.error }}</p>
    <h4>outfmt 6 row</h4>
    <pre class="output row-text" data-testid="detail-row">{{ detail.row }}</pre>
    <h4>outfmt 0</h4>
    <template v-if="entry?.inOutfmt0">
      <p v-if="detail.state === 'loading'" class="muted">Reading the alignment…</p>
      <pre v-if="detail.heading !== undefined" class="output" data-testid="detail-heading">{{ detail.heading }}</pre>
      <pre v-if="detail.section !== undefined" class="output" data-testid="detail-section">{{ detail.section }}</pre>
    </template>
    <p v-else class="notice" data-testid="detail-not-in-outfmt0">
      outfmt 0 does not show this HSP. It shows the alignments of the first {{ state.outfmt0Subjects }} subjects of this query
      (BLAST+'s <span class="option-name">-num_alignments</span>: 250<template v-if="limitGiven"
        >, here <span class="option-name">-max_target_seqs</span></template
      >), and this subject comes after them. The
      outfmt 6 row above and outfmt 7 hold the HSP.
    </p>
  </div>
  <p v-else class="muted">Choose an HSP to see its alignment.</p>
</template>
