<script setup lang="ts">
// LOSAT Web's own files of the run (S15 item 2; design §12.1, §12.3), apart from the compatibility
// outputs above them: CSV, JSON and the static report, which LOSAT Web makes from the HSP records
// for the whole run, after the view filters, or for the subjects marked in the Descriptions
// (application/result-export.ts). Each scope shows its count, and an empty one cannot be chosen.
import { computed, ref, watch } from 'vue';
import type { ExportFormat, ResultExporter } from '../application/result-export';
import type { ResultsState } from '../application/results';
import { filterWords, SCOPE_LABELS, type ExportScope } from '../domain/hsp-export';
import { formatBytes, formatCounted } from './format';
import { useStore } from './useStore';

const props = defineProps<{ exporter: ResultExporter; state: ResultsState }>();
const exportState = useStore(props.exporter.state);
const counts = computed(() => props.exporter.counts(props.state));
const scope = ref<ExportScope>('all');
/** The JSON's aligned rows: on by default (the JSON is the format that keeps everything of an HSP). */
const aligned = ref(true);
const busy = computed(() => exportState.value.busy !== undefined);

// A scope that becomes empty (the marks cleared, a filter that hides every HSP) gives way to the whole run.
watch(counts, (now) => {
  if (scope.value !== 'all' && now[scope.value] === 0) scope.value = 'all';
});

const filters = computed(() => {
  const words = filterWords(props.state.filters);
  return words.length === 0 ? 'no view filter is set' : words.join('; ');
});
/** The query whose marked subjects the marked scope takes. */
const markedQuery = computed(() => {
  const qIdx = props.state.qIdx;
  if (qIdx === undefined) return '';
  const id = props.state.loaded?.run.snapshot.query.records[qIdx]?.id ?? '';
  return id === '' ? `query ${qIdx + 1}` : `query ${id}`;
});

function download(format: ExportFormat): void {
  if (busy.value || counts.value[scope.value] === 0) return;
  void props.exporter.export(format, scope.value, { aligned: aligned.value });
}
</script>

<template>
  <section class="export-formats" data-testid="export-formats">
    <h3>LOSAT Web formats</h3>
    <p class="notice" data-testid="export-formats-note">
      LOSAT Web makes these files from the HSP records of this run. They are application formats of LOSAT Web, not NCBI BLAST formats.
      For the outputs as the engine wrote them, use the export above.
    </p>
    <fieldset class="radios" :disabled="busy">
      <legend>HSPs</legend>
      <label class="radio">
        <input v-model="scope" type="radio" value="all" :disabled="counts.all === 0" data-testid="export-scope-all" :data-count="counts.all" />
        <span>{{ SCOPE_LABELS.all }} ({{ formatCounted(counts.all, 'HSP') }})</span>
      </label>
      <label class="radio">
        <input
          v-model="scope"
          type="radio"
          value="filtered"
          :disabled="counts.filtered === 0"
          data-testid="export-scope-filtered"
          :data-count="counts.filtered"
        />
        <span>{{ SCOPE_LABELS.filtered }} ({{ formatCounted(counts.filtered, 'HSP') }}): <span class="muted">{{ filters }}</span></span>
      </label>
      <label class="radio">
        <input
          v-model="scope"
          type="radio"
          value="marked"
          :disabled="counts.marked === 0"
          data-testid="export-scope-marked"
          :data-count="counts.marked"
        />
        <span v-if="counts.marked > 0">Subjects marked in the Descriptions ({{ formatCounted(counts.marked, 'HSP') }}): {{ markedQuery }}, after the view filters</span>
        <span v-else>Subjects marked in the Descriptions (0 HSPs): <span class="muted">mark subjects in the Descriptions to export their HSPs</span></span>
      </label>
    </fieldset>
    <div class="export-actions">
      <button type="button" :disabled="busy || counts[scope] === 0" data-testid="export-csv" @click="download('csv')">Download CSV</button>
      <button type="button" :disabled="busy || counts[scope] === 0" data-testid="export-json" @click="download('json')">Download JSON</button>
      <label class="export-aligned">
        <input v-model="aligned" type="checkbox" :disabled="busy" data-testid="export-json-aligned" /> Include aligned sequences in the JSON
      </label>
      <button type="button" :disabled="busy || counts[scope] === 0" data-testid="export-report" @click="download('report')">
        Download report (HTML)
      </button>
    </div>
    <ul class="hint export-format-notes">
      <li>
        CSV: a header row, then one row per HSP: the run, the query and subject records, the HSP, then the outfmt 6 fields as the engine wrote
        them, the frames and whether outfmt 0 shows the HSP. RFC 4180: UTF-8, CRLF line ends, a field with a comma, a quote or a line break in
        quotes. CSV has no room for a note, so the file does not say that it is a LOSAT Web format. IDs are written as they are: a spreadsheet
        may read one that starts with =, +, - or @ as a formula.
      </li>
      <li>JSON: the run, the HSPs chosen and, for each HSP, its outfmt 6 fields as written and its record's numbers as the engine wrote them.</li>
      <li>
        Report: one static HTML page with no script that loads nothing: the run, its commands, the HSPs chosen, each query's HSP table and the
        outfmt 0 text of their alignments, and the warnings.
      </li>
    </ul>
    <div aria-live="polite">
      <p v-if="exportState.busy" class="muted" data-testid="export-busy" :data-format="exportState.busy.format">
        Writing {{ exportState.busy.fileName }}…
      </p>
      <p v-else-if="exportState.error" class="error" role="alert" data-testid="export-error">{{ exportState.error }}</p>
      <p
        v-else-if="exportState.last"
        class="export-summary"
        data-testid="export-summary"
        :data-format="exportState.last.format"
        :data-scope="exportState.last.scope"
        :data-hsps="exportState.last.hsps"
      >
        Saved {{ exportState.last.fileName }}: {{ formatCounted(exportState.last.hsps, 'HSP') }} ({{ SCOPE_LABELS[exportState.last.scope] }}),
        {{ formatBytes(exportState.last.bytes) }}.
      </p>
    </div>
  </section>
</template>
