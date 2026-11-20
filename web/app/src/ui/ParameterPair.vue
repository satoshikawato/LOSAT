<script setup lang="ts">
// A row of two fields with NCBI's select of their pairs ("Match/Mismatch Scores" 1,-2, "Gap
// Costs" Existence: 11 Extension: 1): the select shows the pair that the search uses and
// writes both fields; its "Other", or a pair outside the list, shows the two fields
// (ParameterField, inline). A list that the task, the scores or the matrix change gives the row
// its new default (`setFields`), as NCBI's page does.
import { computed, ref, watchEffect } from 'vue';
import type { DraftState, SearchDraft } from '../application/draft';
import { knownDefault, pairEntries, pairList, pairValue, writtenValue, type PairRow } from '../domain/parameters';
import { programById, type ParameterField as Field } from '../domain/programs';
import { optionHelp } from './format';
import HelpTip from './HelpTip.vue';
import ParameterField from './ParameterField.vue';
import { refusal, refusedMark } from './refusal';

const props = defineProps<{ draft: SearchDraft; state: DraftState; pair: PairRow; fields: readonly Field[] }>();

/** The select's entry that shows the two fields. */
const OTHER = 'other';

const id = computed(() => `param-${props.pair.select.id}`);
const program = computed(() => programById(props.state.program));
const options = computed(() => props.state.description?.parameters);
const values = computed(() => props.state.values[props.state.program]);
const list = computed(() => pairList(program.value, props.pair, values.value, options.value));
/** The pair that the search uses ('' while a default is unknown). */
const shown = computed(() => pairValue(program.value, props.pair, values.value, options.value));
const changed = computed(() =>
  props.fields.some(
    (field) => writtenValue(field, values.value[field.flag], knownDefault(program.value, field.flag, values.value, options.value)) !== undefined,
  ),
);

/** The two fields: chosen with "Other", or holding a pair outside the list. */
const other = ref(false);
watchEffect(() => {
  if (props.pair.flags.every((flag) => values.value[flag] === undefined)) other.value = false;
  else if (shown.value !== '' && !list.value.includes(shown.value)) other.value = true;
});
const fieldsShown = computed(() => other.value || list.value.length === 0);
/** A refusal that names one of the two marks the select, or the field shown for that flag (WA2-AX-05). */
const refused = computed(() => {
  const named = refusal(props.state)?.flags;
  return named !== undefined && props.pair.flags.some((flag) => named.has(flag));
});

const label = (value: string) => props.pair.select.labels?.(value) ?? value;
/** The two flags' help, a line each, behind the "?" (HelpTip). */
const help = computed(() =>
  props.fields
    .map((field) => [field.flag, optionHelp(options.value?.find((o) => o.flag === field.flag)?.help)] as const)
    .filter(([, text]) => text !== '')
    .map(([flag, text]) => `${flag}: ${text}`)
    .join('\n'),
);

function onList(next: string): void {
  if (next === OTHER) {
    other.value = true;
    return;
  }
  other.value = false;
  props.draft.setFields(pairEntries(program.value, props.pair, next, values.value, options.value));
}
</script>

<template>
  <div :class="['param-row', 'param-row-group', { changed: changed && !fieldsShown }]" :data-testid="`${id}-field`">
    <!-- The fields shown carry their own names and marks; the select is named by the row. -->
    <component :is="list.length > 0 ? 'label' : 'span'" :id="`${id}-label`" class="param-label" :for="list.length > 0 ? id : undefined">
      <span v-if="changed && !fieldsShown" class="changed-mark" aria-hidden="true">♦</span>
      {{ pair.row }} <span v-if="list.length > 0" class="flag">{{ pair.flags.join(' ') }}</span>
      <span v-if="changed && !fieldsShown" class="visually-hidden">, changed from the default</span>
    </component>
    <div class="param-control">
      <!-- WebKit draws a select's chosen label at its whole length beyond the select's box (S13 screen review H2). -->
      <span v-if="list.length > 0" class="select-clip">
        <select
          :id="id"
          :value="other ? OTHER : shown"
          v-bind="refusedMark(refused && !fieldsShown)"
          :data-testid="id"
          @change="onList(($event.target as HTMLSelectElement).value)"
        >
          <option v-if="shown === ''" value="">Default</option>
          <option v-for="item in list" :key="item" :value="item">{{ label(item) }}</option>
          <option :value="OTHER">Other</option>
        </select>
      </span>
      <div v-if="fieldsShown" class="param-group">
        <ParameterField v-for="field in fields" :key="field.flag" :draft="draft" :state="state" :field="field" inline />
      </div>
      <HelpTip v-if="help !== '' && !fieldsShown" :id="id" :label="`About ${pair.flags.join(' and ')}`" :text="help" />
    </div>
  </div>
</template>
