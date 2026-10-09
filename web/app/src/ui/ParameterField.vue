<script setup lang="ts">
// One option of the search form, with the engine's default, choices and help (`describe`):
// a row (its label left, its control right, as on NCBI's page), or one field of a row of
// several (`inline`, such as "Match" in "Match/Mismatch Scores"). A field whose value is
// written to the argv differs from the engine's default, and is marked as NCBI marks it
// (yellow and "♦"), with words for screen readers. An empty field is the engine's default.
import { computed } from 'vue';
import type { DraftState, SearchDraft } from '../application/draft';
import { geneticCodeLabel } from '../domain/genetic-codes';
import { effectiveValue, fieldChoices, writtenValue } from '../domain/parameters';
import type { ParameterField } from '../domain/programs';

const props = defineProps<{ draft: SearchDraft; state: DraftState; field: ParameterField; inline?: boolean }>();

const id = computed(() => `param-${props.field.flag.slice(1)}`);
const options = computed(() => props.state.description?.parameters);
const option = computed(() => options.value?.find((o) => o.flag === props.field.flag));
const values = computed(() => props.state.values[props.state.program]);
const value = computed(() => values.value[props.field.flag]);
const text = computed(() => (typeof value.value === 'string' ? value.value : ''));
const changed = computed(() => writtenValue(props.field, value.value, option.value) !== undefined);
const choices = computed(() => (props.field.kind === 'boolean' ? ['true', 'false'] : fieldChoices(props.field, option.value)));
/** The radio that is checked: the value set, or else the engine's default. */
const checked = computed(() => effectiveValue(props.field.flag, values.value, options.value));

const choiceLabel = (choice: string) => props.field.choiceLabels?.[choice] ?? choice;

/** The first paragraph of the engine's help, without the references to NCBI's sources. */
const help = computed(() => {
  const raw = option.value?.help ?? '';
  const paragraph = raw.split(/\n\s*\n/)[0] ?? '';
  return paragraph.replace(/\s*(?:NCBI reference|Reference):.*$/s, '').trim();
});

const defaultLabel = computed(() => {
  const fallback = option.value?.defaultValue;
  if (fallback === undefined) return 'Default';
  return props.field.kind === 'gencode' ? `Default (${geneticCodeLabel(Number(fallback))})` : `Default (${choiceLabel(fallback)})`;
});

const gencodes = computed((): readonly number[] => {
  const description = props.state.description;
  const list = props.field.flag === '-query_gencode' ? description?.queryGencodes : description?.subjectGencodes;
  return list ?? [];
});

function set(next: string | boolean): void {
  props.draft.setField(props.field.flag, next);
}
</script>

<template>
  <div
    :class="[inline ? 'param-inline' : 'param-row', { changed }]"
    :data-kind="field.kind"
    :data-testid="`${id}-field`"
  >
    <component :is="field.kind === 'radio' ? 'span' : 'label'" :id="`${id}-label`" class="param-label" :for="field.kind === 'radio' ? undefined : id">
      <span v-if="changed" class="changed-mark" aria-hidden="true">♦</span>
      {{ field.label }} <code class="flag">{{ field.flag }}</code>
      <span v-if="changed" class="visually-hidden">, changed from the default</span>
    </component>
    <div class="param-control">
      <div v-if="field.kind === 'radio'" class="radios" role="radiogroup" :aria-labelledby="`${id}-label`" :data-testid="id">
        <label v-for="choice in choices" :key="choice" class="radio">
          <input
            type="radio"
            :name="id"
            :value="choice"
            :checked="checked === choice"
            :data-testid="`${id}-${choice}`"
            @change="set(choice)"
          />
          {{ choiceLabel(choice) }}
        </label>
      </div>
      <input
        v-else-if="field.kind === 'flag'"
        :id="id"
        type="checkbox"
        :checked="value === true"
        :data-testid="id"
        @change="set(($event.target as HTMLInputElement).checked)"
      />
      <template v-else-if="field.kind === 'text'">
        <input
          :id="id"
          type="text"
          :inputmode="field.inputMode ?? 'text'"
          autocomplete="off"
          spellcheck="false"
          :value="text"
          :placeholder="option?.defaultValue === undefined ? 'default' : `default: ${option.defaultValue}`"
          :list="field.suggestions ? `${id}-list` : undefined"
          :data-testid="id"
          @input="set(($event.target as HTMLInputElement).value)"
        />
        <datalist v-if="field.suggestions" :id="`${id}-list`">
          <option v-for="s in field.suggestions" :key="s" :value="s" />
        </datalist>
      </template>
      <!-- WebKit draws a select's chosen label at its whole length beyond the select's box, which
           made a phone's page wider than its screen (S13 screen review H2): the span clips it. -->
      <span v-else class="select-clip">
        <select
          v-if="field.kind === 'gencode'"
          :id="id"
          :value="text"
          :data-testid="id"
          @change="set(($event.target as HTMLSelectElement).value)"
        >
          <option value="">{{ defaultLabel }}</option>
          <option v-for="code in gencodes" :key="code" :value="String(code)">{{ geneticCodeLabel(code) }}</option>
        </select>
        <select v-else :id="id" :value="text" :data-testid="id" @change="set(($event.target as HTMLSelectElement).value)">
          <option value="">{{ defaultLabel }}</option>
          <option v-for="choice in choices" :key="choice" :value="choice">{{ choiceLabel(choice) }}</option>
        </select>
      </span>
      <p v-if="help" class="help">{{ help }}</p>
    </div>
  </div>
</template>
