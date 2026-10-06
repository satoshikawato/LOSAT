<script setup lang="ts">
// The program's options (plan §5.3): the sections and labels of its descriptor, with the
// engine's defaults, choices and help (`describe`). An empty field is the engine's
// default and is not written to the argv.
import { computed } from 'vue';
import type { DraftState, SearchDraft } from '../application/draft';
import { geneticCodeLabel } from '../domain/genetic-codes';
import { describedSections, fieldChoices } from '../domain/parameters';
import { programById, type ParameterField } from '../domain/programs';
import type { ParameterDescription } from '../ports/engine';

const props = defineProps<{ draft: SearchDraft; state: DraftState }>();
const program = computed(() => programById(props.state.program));
const options = computed(() => new Map((props.state.description?.parameters ?? []).map((option) => [option.flag, option])));
const sections = computed(() => describedSections(program.value, props.state.description?.parameters));
const values = computed(() => props.state.values[props.state.program]);

const testid = (field: ParameterField) => `param-${field.flag.slice(1)}`;
const option = (field: ParameterField): ParameterDescription | undefined => options.value.get(field.flag);
const text = (field: ParameterField) => {
  const value = values.value[field.flag];
  return typeof value === 'string' ? value : '';
};

/** The first paragraph of the engine's help, without the references to NCBI's sources. */
function help(field: ParameterField): string {
  const raw = option(field)?.help ?? '';
  const paragraph = raw.split(/\n\s*\n/)[0] ?? '';
  return paragraph.replace(/\s*(?:NCBI reference|Reference):.*$/s, '').trim();
}

function defaultLabel(field: ParameterField): string {
  const value = option(field)?.defaultValue;
  if (value === undefined) return 'Default';
  return field.kind === 'gencode' ? `Default (${geneticCodeLabel(Number(value))})` : `Default (${value})`;
}

function gencodes(field: ParameterField): readonly number[] {
  const description = props.state.description;
  const list = field.flag === '-query_gencode' ? description?.queryGencodes : description?.subjectGencodes;
  return list ?? [];
}

/** A subject genetic code other than the default: LOSAT's approved exception (AGENTS.md). */
const subjectCodeNote = computed(() => {
  const value = values.value['-db_gencode'];
  const fallback = options.value.get('-db_gencode')?.defaultValue;
  return typeof value === 'string' && value !== '' && value !== fallback;
});

function set(field: ParameterField, value: string | boolean): void {
  props.draft.setField(field.flag, value);
}
</script>

<template>
  <div class="parameter-form" data-testid="parameter-form">
    <p v-if="state.descriptionError" class="error">The engine's options could not be read: {{ state.descriptionError }}</p>
    <fieldset v-for="section in sections" :key="section.title" class="parameter-section">
      <legend>{{ section.title }}</legend>
      <div class="fields">
        <div v-for="field in section.fields" :key="field.flag" class="field" :data-kind="field.kind">
          <label v-if="field.kind === 'flag'" class="flag-label">
            <input
              type="checkbox"
              :checked="values[field.flag] === true"
              :data-testid="testid(field)"
              @change="set(field, ($event.target as HTMLInputElement).checked)"
            />
            {{ field.label }} <code class="flag">{{ field.flag }}</code>
          </label>
          <template v-else>
            <label :for="testid(field)">
              {{ field.label }} <code class="flag">{{ field.flag }}</code>
            </label>
            <template v-if="field.kind === 'text'">
              <input
              :id="testid(field)"
              type="text"
              :inputmode="field.inputMode ?? 'text'"
              autocomplete="off"
              spellcheck="false"
              :value="text(field)"
              :placeholder="option(field)?.defaultValue === undefined ? 'default' : `default: ${option(field)!.defaultValue}`"
              :list="field.suggestions ? `${testid(field)}-list` : undefined"
              :data-testid="testid(field)"
              @input="set(field, ($event.target as HTMLInputElement).value)"
            />
              <datalist v-if="field.suggestions" :id="`${testid(field)}-list`">
                <option v-for="s in field.suggestions" :key="s" :value="s" />
              </datalist>
            </template>
            <select
              v-else-if="field.kind === 'choice' || field.kind === 'boolean'"
              :id="testid(field)"
              :value="text(field)"
              :data-testid="testid(field)"
              @change="set(field, ($event.target as HTMLSelectElement).value)"
            >
              <option value="">{{ defaultLabel(field) }}</option>
              <option
                v-for="choice in field.kind === 'boolean' ? ['true', 'false'] : fieldChoices(field, option(field))"
                :key="choice"
                :value="choice"
              >
                {{ choice }}
              </option>
            </select>
            <select
              v-else-if="field.kind === 'gencode'"
              :id="testid(field)"
              :value="text(field)"
              :data-testid="testid(field)"
              @change="set(field, ($event.target as HTMLSelectElement).value)"
            >
              <option value="">{{ defaultLabel(field) }}</option>
              <option v-for="code in gencodes(field)" :key="code" :value="String(code)">{{ geneticCodeLabel(code) }}</option>
            </select>
          </template>
          <p v-if="help(field)" class="help">{{ help(field) }}</p>
        </div>
      </div>
      <p v-if="section.title.startsWith('Genetic code') && subjectCodeNote" class="note" data-testid="subject-gencode-note">
        Approved LOSAT exception: LOSAT translates the subjects with this code. NCBI BLAST+ treats a non-default subject
        code differently for local subject files, so results with this code can differ from NCBI's.
      </p>
    </fieldset>
  </div>
</template>
