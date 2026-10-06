<script setup lang="ts">
import { computed, ref } from 'vue';
import type { DraftState, SearchDraft } from '../application/draft';
import { includedOf } from '../application/draft';
import { programById, residueUnit, sequenceKind, type InputRole } from '../domain/programs';
import RegionPicker from './RegionPicker.vue';
import SourceCard from './SourceCard.vue';

const props = defineProps<{ draft: SearchDraft; state: DraftState; role: InputRole }>();
const title = computed(() => (props.role === 'query' ? 'Query' : 'Subject'));
const roleState = computed(() => props.state[props.role]);
const kind = computed(() => sequenceKind(programById(props.state.program), props.role));
const dragging = ref(false);
const fileInput = ref<HTMLInputElement>();

const readySources = computed(() => roleState.value.sources.filter((source) => source.status === 'ready'));
const includedCount = computed(() => readySources.value.reduce((sum, source) => sum + includedOf(source).length, 0));
const regionRecord = computed(() => {
  void props.state;
  return props.draft.regionRecord(props.role);
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
  <section class="input-panel" :data-testid="`${role}-panel`">
    <header class="input-header">
      <h3>{{ title }}</h3>
      <span class="kind">{{ kind }} ({{ residueUnit(kind) }})</span>
    </header>
    <div
      class="dropzone"
      :class="{ dragging }"
      :data-testid="`${role}-dropzone`"
      @dragover.prevent="dragging = true"
      @dragleave="dragging = false"
      @drop.prevent="onDrop"
    >
      <label class="paste">
        <span class="visually-hidden">Paste {{ title.toLowerCase() }} sequences</span>
        <textarea
          :value="roleState.paste"
          rows="5"
          spellcheck="false"
          autocomplete="off"
          :placeholder="`Paste ${title.toLowerCase()} sequences in FASTA format`"
          :data-testid="`${role}-input`"
          @input="onPaste"
        />
      </label>
      <div class="file-actions">
        <button type="button" :data-testid="`${role}-open`" @click="fileInput?.click()">Open FASTA files…</button>
        <input
          ref="fileInput"
          class="visually-hidden"
          type="file"
          multiple
          tabindex="-1"
          :data-testid="`${role}-files`"
          @change="onFiles"
        />
        <span class="hint">or drop files here. Files are read in place, not copied into the box.</span>
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

    <RegionPicker
      v-if="regionRecord"
      :draft="draft"
      :role="role"
      :record="regionRecord"
      :region="roleState.region"
      :unit="residueUnit(kind)"
    />
    <p v-else-if="includedCount > 1" class="hint" :data-testid="`${role}-region-unavailable`">
      A region can be set when the {{ title.toLowerCase() }} has one record ({{ includedCount }} are included).
    </p>
  </section>
</template>
