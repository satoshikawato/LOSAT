<script setup lang="ts">
// The sections under "Algorithm parameters" (plan §5.3; NCBI's General Parameters, Scoring
// Parameters, Filters and Masking, Discontiguous Word Options, and LOSAT's Other Parameters),
// whose rows share one label column, as on NCBI's form:
// the sections and labels of the program's descriptor, with the engine's defaults, choices
// and help (`describe`); a row of two fields with NCBI's select of their pairs is a
// ParameterPair. The last section ends with "Additional options": the other options of the
// command line, as typed there, which the engine's `validate` reads (DW-42). The task and the
// genetic codes are shown elsewhere (SearchPanel, InputPanel).
import { computed } from 'vue';
import type { DraftState, SearchDraft } from '../application/draft';
import { fieldRows, pairRows, sectionShown, sectionsAt, type PairRow } from '../domain/parameters';
import { programById, type ParameterField as Field } from '../domain/programs';
import HelpTip from './HelpTip.vue';
import ParameterField from './ParameterField.vue';
import ParameterPair from './ParameterPair.vue';
import { refusal, refusedMark } from './refusal';

const props = defineProps<{ draft: SearchDraft; state: DraftState }>();
const sectionId = (title: string) => `parameters-${title.toLowerCase().replace(/[^a-z]+/g, '-')}`;
/** The row's NCBI select of pairs, where the program has one for it and the engine describes both fields. */
function pairOf(row: string | undefined, fields: readonly Field[]): PairRow | undefined {
  const pair = pairRows(programById(props.state.program)).find((p) => p.row === row);
  return pair !== undefined && fields.length === 2 && fields.every((field) => pair.flags.includes(field.flag)) ? pair : undefined;
}
const sections = computed(() => {
  const options = props.state.description?.parameters;
  const values = props.state.values[props.state.program];
  return sectionsAt(programById(props.state.program), options, 'algorithm').filter((section) =>
    sectionShown(section, values, options),
  );
});
const extra = computed(() => props.state.extra[props.state.program]);
/** Words typed in "Additional options" are written to the argv: the row is marked as a changed field is. */
const extraChanged = computed(() => extra.value.trim() !== '');
/** A refusal of a word typed there, or of its text, marks the row's text box (WA2-AX-05). */
const extraRefused = computed(() => refusal(props.state)?.extra === true);
</script>

<template>
  <div class="parameter-form" data-testid="parameter-form">
    <p v-if="state.descriptionError" class="error">The engine's options could not be read: {{ state.descriptionError }}</p>
    <!-- A group, not a fieldset: its rows share the form's label column (CSS subgrid), which
         Chromium does not give a fieldset. Keyed by program: a control's own state ("Other", an
         open help) does not pass to another program's control of the same flag. -->
    <div
      v-for="(section, index) in sections"
      :key="`${state.program} ${section.title}`"
      class="search-block parameter-section"
      role="group"
      :aria-labelledby="sectionId(section.title)"
    >
      <div :id="sectionId(section.title)" class="section-legend">{{ section.title }}</div>
      <template v-for="row in fieldRows(section.fields)" :key="row.fields[0]!.flag">
        <ParameterField v-if="row.row === undefined" :draft="draft" :state="state" :field="row.fields[0]!" />
        <ParameterPair
          v-else-if="pairOf(row.row, row.fields)"
          :draft="draft"
          :state="state"
          :pair="pairOf(row.row, row.fields)!"
          :fields="row.fields"
        />
        <div v-else class="param-row param-row-group">
          <span class="param-label">{{ row.row }}</span>
          <div class="param-control param-group">
            <ParameterField v-for="field in row.fields" :key="field.flag" :draft="draft" :state="state" :field="field" inline />
          </div>
        </div>
      </template>
      <div v-if="index === sections.length - 1" class="param-row" :class="{ changed: extraChanged }" data-kind="extra" data-testid="param-extra-field">
        <label id="param-extra-label" class="param-label" for="param-extra">
          <span v-if="extraChanged" class="changed-mark" aria-hidden="true">♦</span>
          Additional options
          <span v-if="extraChanged" class="visually-hidden">, changed from the default</span>
        </label>
        <div class="param-control">
          <input
            id="param-extra"
            class="extra-options"
            type="text"
            autocomplete="off"
            autocapitalize="off"
            spellcheck="false"
            :value="extra"
            v-bind="refusedMark(extraRefused)"
            data-testid="param-extra"
            @input="draft.setExtra(($event.target as HTMLInputElement).value)"
          />
          <HelpTip
            id="param-extra"
            label="About Additional options"
            text="Other options of the losat command line, written as there: -qcov_hsp_perc 50 -ungapped"
          />
        </div>
      </div>
    </div>
  </div>
</template>
