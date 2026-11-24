<script setup lang="ts">
// Parameter sections driven by a `FormModel` (Search's "Algorithm parameters" and tool forms'
// option sections): rows share one label column via CSS subgrid, with optional pair selects,
// "Additional options" and a `-conffile` picker at the end of the last section.
import { computed, ref } from 'vue';
import { fieldRows, type PairRow } from '../domain/parameters';
import type { ParameterField as Field } from '../domain/programs';
import { formatBytes } from './format';
import type { FormModel } from './form-model';
import HelpTip from './HelpTip.vue';
import ParameterField from './ParameterField.vue';
import ParameterPair from './ParameterPair.vue';
import { REFUSAL_ID, refusedMark } from './refusal';

const props = withDefaults(defineProps<{ model: FormModel; testId?: string }>(), { testId: 'parameter-form' });
const prefix = computed(() => props.model.idPrefix ?? 'param');
const sectionPrefix = computed(() => (prefix.value === 'param' ? 'parameters' : `${prefix.value}-section`));
const sectionId = (title: string) => `${sectionPrefix.value}-${title.toLowerCase().replace(/[^a-z]+/g, '-')}`;
const pairOf = (row: string | undefined, fields: readonly Field[]): PairRow | undefined => props.model.pairOf?.(row, fields);
const sections = computed(() => (props.model.descriptionUnavailable === true ? [] : (props.model.sections ?? [])));
const extra = computed(() => props.model.extra ?? '');
/** Words typed in "Additional options" are written to the argv: the row is marked as a changed field is. */
const extraChanged = computed(() => extra.value.trim() !== '');
/** A refusal of a word typed there, or of its text, marks the row's text box (WA2-AX-05). */
const extraRefused = computed(() => props.model.refusal?.extra === true);
const refusalId = computed(() => props.model.refusalId ?? REFUSAL_ID);
const conffileName = computed(() => props.model.conffileName);
const conffile = computed(() => props.model.conffile);
const conffileInput = ref<HTMLInputElement | null>(null);

function onConffile(event: Event): void {
  const input = event.target as HTMLInputElement;
  const file = input.files?.[0];
  if (file !== undefined) props.model.setConffile?.(file);
  // The same file can be chosen again after it is removed.
  input.value = '';
}
</script>

<template>
  <div class="parameter-form" :data-testid="testId">
    <!-- A part of LOSAT that could not be loaded is said once, for the page (AppView; line WQA's sweep, D1); the options
         were not read, which the form says once, as an input card does, where it showed every number field blank (the
         review of D1, L1, 2026-11-19). -->
    <p v-if="model.descriptionError && !model.descriptionUnavailable" class="error">The engine's options could not be read: {{ model.descriptionError }}</p>
    <p v-if="model.descriptionUnavailable" class="muted" :data-testid="prefix === 'param' ? 'parameters-unread' : `${prefix}-unread`">Not read.</p>
    <!-- A group, not a fieldset: its rows share the form's label column (CSS subgrid), which
         Chromium does not give a fieldset. Keyed by mode/stage: a control's own state ("Other", an
         open help) does not pass to another mode's control of the same flag. -->
    <div
      v-for="(section, index) in sections"
      :key="`${model.key} ${section.title}`"
      class="search-block parameter-section"
      role="group"
      :aria-labelledby="sectionId(section.title)"
    >
      <div :id="sectionId(section.title)" class="section-legend">{{ section.title }}</div>
      <template v-for="row in fieldRows(section.fields)" :key="row.fields[0]!.flag">
        <ParameterField v-if="row.row === undefined" :model="model" :field="row.fields[0]!" />
        <ParameterPair
          v-else-if="pairOf(row.row, row.fields)"
          :model="model"
          :pair="pairOf(row.row, row.fields)!"
          :fields="row.fields"
        />
        <div v-else class="param-row param-row-group">
          <span class="param-label">{{ row.row }}</span>
          <div class="param-control param-group">
            <!-- A taxonomy ID is shown in the row of the field that it names the counts for (ParameterField). -->
            <ParameterField
              v-for="field in row.fields.filter((f) => f.kind !== 'taxid')"
              :key="field.flag"
              :model="model"
              :field="field"
              inline
            />
          </div>
        </div>
      </template>
      <slot v-if="index === sections.length - 1" name="before-extra" />
      <div
        v-if="index === sections.length - 1 && model.setExtra"
        class="param-row"
        :class="{ changed: extraChanged }"
        data-kind="extra"
        :data-testid="`${prefix}-extra-field`"
      >
        <label :id="`${prefix}-extra-label`" class="param-label" :for="`${prefix}-extra`">
          <span v-if="extraChanged" class="changed-mark" aria-hidden="true">♦</span>
          Additional options
          <span v-if="extraChanged" class="visually-hidden">, changed from the default</span>
        </label>
        <div class="param-control">
          <input
            :id="`${prefix}-extra`"
            class="extra-options"
            type="text"
            autocomplete="off"
            autocapitalize="off"
            spellcheck="false"
            :value="extra"
            v-bind="refusedMark(extraRefused, refusalId)"
            :data-testid="`${prefix}-extra`"
            @input="model.setExtra(($event.target as HTMLInputElement).value)"
          />
          <HelpTip
            v-if="model.extraHelp"
            :id="`${prefix}-extra`"
            label="About Additional options"
            :text="model.extraHelp"
          />
        </div>
      </div>
      <!-- The file that -conffile names, chosen as the counts file of a row is (ParameterField): the search reads it under that name. -->
      <div
        v-if="index === sections.length - 1 && conffileName !== undefined && model.setConffile"
        class="param-row"
        :class="{ changed: conffile !== undefined }"
        data-kind="conffile"
        :data-testid="`${prefix}-conffile`"
      >
        <span :id="`${prefix}-conffile-label`" class="param-label">
          <span v-if="conffile" class="changed-mark" aria-hidden="true">♦</span>
          <span class="flag">-conffile</span> <span :data-testid="`${prefix}-conffile-argv`">{{ conffileName }}</span>
        </span>
        <div class="param-control">
          <span class="file-pick">
            <button type="button" :aria-describedby="`${prefix}-conffile-label`" :data-testid="`${prefix}-conffile-choose`" @click="conffileInput?.click()">
              {{ conffile ? 'Change file…' : 'Choose file…' }}
            </button>
            <!-- Opened by the button: out of the accessibility tree, as the inputs' file field (WA2-AX-13). -->
            <input
              ref="conffileInput"
              class="visually-hidden"
              type="file"
              tabindex="-1"
              aria-hidden="true"
              aria-label="The file of -conffile"
              :data-testid="`${prefix}-conffile-file`"
              @change="onConffile"
            />
            <template v-if="conffile">
              <span class="file-name" :data-testid="`${prefix}-conffile-name`">{{ conffile.name }} ({{ formatBytes(conffile.size) }})</span>
              <button type="button" class="link" :data-testid="`${prefix}-conffile-remove`" @click="model.setConffile(undefined)">Remove</button>
            </template>
            <span v-else class="file-missing" role="status" :data-testid="`${prefix}-conffile-missing`">Choose the file that the search reads as {{ conffileName }}.</span>
          </span>
          <HelpTip
            :id="`${prefix}-conffile`"
            label="About -conffile"
            text="The settings file that the losat command line reads with -conffile, in place of the program's .ini file: the search reads the file chosen here under the name in Additional options."
          />
        </div>
      </div>
    </div>
  </div>
</template>
