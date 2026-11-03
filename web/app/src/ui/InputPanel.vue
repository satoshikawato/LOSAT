<script setup lang="ts">
// A role's block ("Enter Query Sequence", "Enter Subject Sequence"), in the order of NCBI's
// two-sequence page: the paste box (FASTA, or accessions to fetch) with its subrange beside it, "Or, upload file", the
// role's genetic code (BLASTX, TBLASTN, TBLASTX), the Job Title (query), then LOSAT's own parts: the
// inputs with their records, and Combined / Separate. A file dropped anywhere on the block is
// opened.
import { computed, ref } from 'vue';
import type { DraftState, SearchDraft } from '../application/draft';
import { accessionLines, FETCH_SOURCES, type FetchSourceId } from '../domain/accession';
import { sectionsAt, writtenValue } from '../domain/parameters';
import { programById, residueUnit, sequenceKind, type InputRole } from '../domain/programs';
import HelpTip from './HelpTip.vue';
import LosatLoader from './LosatLoader.vue';
import ParameterField from './ParameterField.vue';
import RegionPicker from './RegionPicker.vue';
import SourceCard from './SourceCard.vue';

const props = defineProps<{ draft: SearchDraft; state: DraftState; role: InputRole }>();
const title = computed(() => (props.role === 'query' ? 'Query' : 'Subject'));
const roleState = computed(() => props.state[props.role]);
const program = computed(() => programById(props.state.program));
const kind = computed(() => sequenceKind(program.value, props.role));
const dragging = ref(false);
const fileInput = ref<HTMLInputElement>();

/** The box's accessions, when it holds accessions only: "Fetch" brings their FASTA into it. */
const accessions = computed(() => accessionLines(roleState.value.paste));
const fetchState = computed(() => roleState.value.fetch);
const sourceLabel = (id: FetchSourceId) => FETCH_SOURCES.find((source) => source.id === id)!.label;
const fetchSource = computed(() => sourceLabel(props.state.fetchSource));

const readySources = computed(() => roleState.value.sources.filter((source) => source.status === 'ready'));
const regionRecord = computed(() => {
  void props.state;
  return props.draft.regionRecord(props.role);
});
const region = computed(() => {
  void props.state;
  return props.draft.region(props.role);
});

/** The role's genetic code (the query's of BLASTX and TBLASTX, the subjects' of TBLASTN and TBLASTX). */
const codeFields = computed(() =>
  sectionsAt(program.value, props.state.description?.parameters, props.role).flatMap((section) => section.fields),
);
/** A subject genetic code other than the default: LOSAT's approved exception (AGENTS.md). */
const subjectCodeNote = computed(() => {
  const field = codeFields.value.find((f) => f.flag === '-db_gencode');
  if (props.role !== 'subject' || field === undefined) return false;
  const option = props.state.description?.parameters.find((o) => o.flag === field.flag);
  return writtenValue(field, props.state.values[props.state.program][field.flag], option) !== undefined;
});

function onPaste(event: Event): void {
  props.draft.setPaste(props.role, (event.target as HTMLTextAreaElement).value);
}

function onFetchSource(event: Event): void {
  props.draft.setFetchSource((event.target as HTMLSelectElement).value as FetchSourceId);
}

function onFiles(event: Event): void {
  const input = event.target as HTMLInputElement;
  props.draft.addFiles(props.role, [...(input.files ?? [])]);
  input.value = '';
}

function onDrop(event: DragEvent): void {
  dragging.value = false;
  const files = [...(event.dataTransfer?.files ?? [])];
  if (files.length > 0) props.draft.addFiles(props.role, files);
}
</script>

<template>
  <fieldset
    class="search-block input-panel"
    :class="{ dragging }"
    :data-testid="`${role}-panel`"
    @dragover.prevent="dragging = true"
    @dragleave.self="dragging = false"
    @drop.prevent.stop="onDrop"
  >
    <legend>Enter {{ title }} Sequence</legend>
    <div class="entry">
      <div class="entry-paste">
        <div class="entry-label">
          <label :for="`${role}-input`">Enter accession number(s) or FASTA sequence(s)</label>
          <span class="kind">{{ kind }}</span>
          <button type="button" class="link" :data-testid="`${role}-clear`" @click="draft.setPaste(role, '')">Clear</button>
        </div>
        <textarea
          :id="`${role}-input`"
          :value="roleState.paste"
          rows="5"
          spellcheck="false"
          autocomplete="off"
          :data-testid="`${role}-input`"
          @input="onPaste"
        />
        <div v-if="accessions" class="fetch-row" :data-testid="`${role}-fetch-row`">
          <button
            type="button"
            :disabled="fetchState?.state === 'fetching'"
            :data-testid="`${role}-fetch`"
            @click="draft.fetchAccessions(role)"
          >
            Fetch
          </button>
          <label>
            from
            <select
              :value="state.fetchSource"
              :disabled="fetchState?.state === 'fetching'"
              :data-testid="`${role}-fetch-source`"
              @change="onFetchSource"
            >
              <option v-for="source in FETCH_SOURCES" :key="source.id" :value="source.id">{{ source.label }}</option>
            </select>
          </label>
          <span class="hint" :data-testid="`${role}-fetch-note`">
            Sends the {{ accessions.length === 1 ? 'accession' : 'accessions' }} to {{ fetchSource }}
          </span>
          <HelpTip
            :id="`${role}-fetch`"
            label="About fetching accessions"
            text="Fetch sends only the accessions, to the site chosen here, and writes the FASTA it gets into the box. Sequences, files and results stay in this browser. EBI answers from ENA for nucleotides and from UniProt for proteins, CNCB from GenBase."
          />
        </div>
        <p
          v-if="fetchState?.state === 'fetching'"
          class="muted"
          role="status"
          :data-testid="`${role}-fetch-status`"
        >
          <LosatLoader />Fetching {{ fetchState.done + 1 }} of {{ fetchState.total }} from {{ sourceLabel(fetchState.source) }}…
        </p>
        <p v-else-if="fetchState?.state === 'failed'" class="error" role="alert" :data-testid="`${role}-fetch-error`">
          {{ fetchState.message }}
        </p>
      </div>
      <div class="entry-subrange">
        <RegionPicker
          v-if="regionRecord"
          :draft="draft"
          :role="role"
          :title="`${title} subrange`"
          :record="regionRecord"
          :region="region"
          :unit="residueUnit(kind)"
        />
        <!-- With no single record, no subrange: nothing says so (stop-describing-ui rule 7). -->
      </div>
    </div>

    <div class="param-row file-row">
      <span class="param-label">Or, upload file</span>
      <div class="param-control dropzone" :data-testid="`${role}-dropzone`">
        <button type="button" :data-testid="`${role}-open`" @click="fileInput?.click()">Open FASTA files…</button>
        <input
          ref="fileInput"
          class="visually-hidden"
          type="file"
          multiple
          tabindex="-1"
          :aria-label="`Open ${title.toLowerCase()} FASTA files`"
          :data-testid="`${role}-files`"
          @change="onFiles"
        />
        <span class="hint">or drop files here</span>
      </div>
    </div>

    <ParameterField v-for="field in codeFields" :key="field.flag" :draft="draft" :state="state" :field="field" />
    <p v-if="subjectCodeNote" class="note" data-testid="subject-gencode-note">
      Results with this code can differ from NCBI BLAST+, which translates local subjects otherwise (approved LOSAT exception).
    </p>

    <div v-if="role === 'query'" class="param-row">
      <label class="param-label" for="job-title">Job Title</label>
      <div class="param-control fill">
        <input
          id="job-title"
          type="text"
          autocomplete="off"
          :value="state.title"
          aria-describedby="job-title-help"
          data-testid="job-title"
          @input="draft.setTitle(($event.target as HTMLInputElement).value)"
        />
        <p id="job-title-help" class="help">Enter a descriptive title for your search</p>
      </div>
    </div>

    <ul v-if="roleState.sources.length > 0" class="sources">
      <SourceCard
        v-for="source in roleState.sources"
        :key="source.key"
        :draft="draft"
        :state="state"
        :role="role"
        :source="source"
        :index="roleState.sources.indexOf(source)"
      />
    </ul>

    <fieldset v-if="readySources.length > 1" class="mode" :data-testid="`${role}-mode`">
      <legend>Several {{ title.toLowerCase() }} inputs</legend>
      <label>
        <input
          type="radio"
          :name="`${role}-mode`"
          value="combined"
          :checked="roleState.mode === 'combined'"
          :data-testid="`${role}-mode-combined`"
          @change="draft.setMode(role, 'combined')"
        />
        Combined: one search of all
      </label>
      <label>
        <input
          type="radio"
          :name="`${role}-mode`"
          value="separate"
          :checked="roleState.mode === 'separate'"
          :data-testid="`${role}-mode-separate`"
          @change="draft.setMode(role, 'separate')"
        />
        Separate: one search for each
      </label>
    </fieldset>
  </fieldset>
</template>
