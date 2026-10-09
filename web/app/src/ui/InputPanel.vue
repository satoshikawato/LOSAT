<script setup lang="ts">
// A role's block ("Enter Query Sequence", "Enter Subject Sequence"), in the order of NCBI's
// two-sequence page: the paste box with its subrange beside it, "Or, upload file", the
// role's genetic code (TBLASTN, TBLASTX), the Job Title (query), then LOSAT's own parts: the
// inputs with their records, and Combined / Separate. A file dropped anywhere on the block is
// opened.
import { computed, ref } from 'vue';
import type { DraftState, SearchDraft } from '../application/draft';
import { includedOf } from '../application/draft';
import { sectionsAt, writtenValue } from '../domain/parameters';
import { programById, residueUnit, sequenceKind, type InputRole } from '../domain/programs';
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

const readySources = computed(() => roleState.value.sources.filter((source) => source.status === 'ready'));
const includedCount = computed(() => readySources.value.reduce((sum, source) => sum + includedOf(source).length, 0));
const regionRecord = computed(() => {
  void props.state;
  return props.draft.regionRecord(props.role);
});
const region = computed(() => {
  void props.state;
  return props.draft.region(props.role);
});

/** The role's genetic code (the query's of TBLASTX, the subjects' of TBLASTN and TBLASTX). */
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
          <label :for="`${role}-input`">Enter FASTA sequence(s)</label>
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
        <template v-else>
          <p class="subrange-title">{{ title }} subrange</p>
          <p v-if="includedCount > 1" class="hint" :data-testid="`${role}-region-unavailable`">
            A {{ title.toLowerCase() }} region can be set when exactly one {{ title.toLowerCase() }} record is included
            ({{ includedCount }} are).
          </p>
          <p v-else class="hint">From and To can be set once one {{ title.toLowerCase() }} record is included.</p>
        </template>
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
        <span class="hint">or drop files here. Files are read in place, not copied into the box.</span>
      </div>
    </div>

    <ParameterField v-for="field in codeFields" :key="field.flag" :draft="draft" :state="state" :field="field" />
    <p v-if="subjectCodeNote" class="note" data-testid="subject-gencode-note">
      Approved LOSAT exception: LOSAT translates the subjects with this code. NCBI BLAST+ treats a non-default subject
      code differently for local subject files, so results with this code can differ from NCBI's.
    </p>

    <div v-if="role === 'query'" class="param-row">
      <label class="param-label" for="job-title">Job Title</label>
      <div class="param-control">
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
        Combined: one search of all of them
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
        Separate: one search for each, queued as a group
      </label>
    </fieldset>
  </fieldset>
</template>
