<script setup lang="ts">
// A role's block ("Enter Query Sequence", "Enter Subject Sequence"), in the order of NCBI's
// two-sequence page: the paste box (FASTA, or accessions to fetch) with its subrange beside it, "Or, upload file", the
// role's genetic code (BLASTX, TBLASTN, TBLASTX), the Job Title (query), then LOSAT's own parts: the
// inputs with their records, and Combined / Separate. A file dropped anywhere on the block is
// opened.
import { computed, ref } from 'vue';
import type { DraftState, SearchDraft } from '../application/draft';
import { raceChoice } from '../application/source-race';
import {
  accessionLines,
  FETCH_SOURCES,
  RACE_RECORDS,
  servingSources,
  sourceLabel,
  sourceList,
  type FetchChoice,
} from '../domain/accession';
import { sectionsAt, writtenValue } from '../domain/parameters';
import { programById, residueUnit, sequenceKind, type InputRole } from '../domain/programs';
import HelpTip from './HelpTip.vue';
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
const busy = (role: InputRole) => {
  const state = props.state[role].fetch?.state;
  return state === 'fetching' || state === 'racing';
};
const fetchBusy = computed(() => busy(props.role));
/** The choice is shared by both boxes: it stays while either fetches, so that each line names where its fetch goes. */
const choiceLocked = computed(() => busy('query') || busy('subject'));
/** The sources that serve the box's accessions, and where "Fastest" sends them (once the race says). */
const serving = computed(() => (accessions.value === undefined ? [] : servingSources(kind.value, accessions.value)));
const fastest = computed(() =>
  serving.value.length === 1 ? serving.value[0] : raceChoice(props.state.race[kind.value], serving.value),
);
/** The in-place sentence: where Fetch sends the accessions. */
const fetchNote = computed(() => {
  const noun = accessions.value?.length === 1 ? 'accession' : 'accessions';
  const choice = props.state.fetchSource;
  if (choice !== 'fastest') return `Sends the ${noun} to ${sourceLabel(choice)}`;
  if (serving.value.length === 0) return 'No one site has all of these accessions';
  if (fastest.value !== undefined) return `Sends the ${noun} to ${sourceLabel(fastest.value)}`;
  return `Sends the ${noun} to the fastest of ${sourceList(serving.value, 'and')}`;
});
const fetchHelp = computed(() => {
  const records = RACE_RECORDS[kind.value];
  const record = records.ebi === records.ncbi ? records.ncbi : `${records.ncbi}, at UniProt ${records.ebi}`;
  return (
    'Fetch sends only the accessions, to the site chosen here, and writes the FASTA it gets into the box. Sequences, ' +
    'files and results stay in this browser. EBI answers from ENA for nucleotides and from UniProt for proteins, CNCB ' +
    `from GenBase. Fastest: the first fetch asks each site for one public record (${record}) and sends the ` +
    'accessions to the first that answers; this tab keeps that order.'
  );
});

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
  props.draft.setFetchSource((event.target as HTMLSelectElement).value as FetchChoice);
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
            :disabled="fetchBusy"
            :data-testid="`${role}-fetch`"
            @click="draft.fetchAccessions(role)"
          >
            Fetch
          </button>
          <label>
            from
            <select
              :value="state.fetchSource"
              :disabled="choiceLocked"
              :data-testid="`${role}-fetch-source`"
              @change="onFetchSource"
            >
              <option value="fastest">Fastest{{ fastest === undefined ? '' : ` (${sourceLabel(fastest)})` }}</option>
              <option v-for="source in FETCH_SOURCES" :key="source.id" :value="source.id">{{ source.label }}</option>
            </select>
          </label>
          <span class="hint" :data-testid="`${role}-fetch-note`">{{ fetchNote }}</span>
          <HelpTip :id="`${role}-fetch`" label="About fetching accessions" :text="fetchHelp" />
        </div>
        <p v-if="fetchState?.state === 'racing'" class="muted" role="status" :data-testid="`${role}-fetch-status`">
          Finding the fastest site…
        </p>
        <p
          v-else-if="fetchState?.state === 'fetching'"
          class="muted"
          role="status"
          :data-testid="`${role}-fetch-status`"
        >
          Fetching {{ fetchState.done + 1 }} of {{ fetchState.total }} from {{ sourceLabel(fetchState.source) }}…
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
