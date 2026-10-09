<script setup lang="ts">
// The sections under "Algorithm parameters" (plan §5.3; NCBI's General Parameters, Scoring
// Parameters, Filters and Masking, Discontiguous Word Options, and LOSAT's Other Parameters):
// the sections and labels of the program's descriptor, with the engine's defaults, choices
// and help (`describe`). The task and the genetic codes are shown elsewhere (SearchPanel,
// InputPanel).
import { computed } from 'vue';
import type { DraftState, SearchDraft } from '../application/draft';
import { fieldRows, sectionShown, sectionsAt } from '../domain/parameters';
import { programById } from '../domain/programs';
import ParameterField from './ParameterField.vue';

const props = defineProps<{ draft: SearchDraft; state: DraftState }>();
const sections = computed(() => {
  const options = props.state.description?.parameters;
  const values = props.state.values[props.state.program];
  return sectionsAt(programById(props.state.program), options, 'algorithm').filter((section) =>
    sectionShown(section, values, options),
  );
});
</script>

<template>
  <div class="parameter-form" data-testid="parameter-form">
    <p v-if="state.descriptionError" class="error">The engine's options could not be read: {{ state.descriptionError }}</p>
    <fieldset v-for="section in sections" :key="section.title" class="search-block parameter-section">
      <legend>{{ section.title }}</legend>
      <template v-for="row in fieldRows(section.fields)" :key="row.fields[0]!.flag">
        <ParameterField v-if="row.row === undefined" :draft="draft" :state="state" :field="row.fields[0]!" />
        <div v-else class="param-row param-row-group">
          <span class="param-label">{{ row.row }}</span>
          <div class="param-control param-group">
            <ParameterField v-for="field in row.fields" :key="field.flag" :draft="draft" :state="state" :field="field" inline />
          </div>
        </div>
      </template>
    </fieldset>
  </div>
</template>
