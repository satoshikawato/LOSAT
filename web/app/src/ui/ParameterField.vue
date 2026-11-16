<script setup lang="ts">
// One option of the search form, with its default, choices and help (`describe`, with the
// defaults that the engine resolves): a row (its label left, its control right, as on NCBI's page),
// or one field of a row of several (`inline`, such as "Match" in "Match/Mismatch Scores", or a
// checkbox of NCBI's "Filter" and "Mask" rows). The control shows the default while the field is empty, as
// NCBI's form shows its values. A field whose value is written to the argv differs from the
// default, and is marked as NCBI marks it (yellow and "♦"), with words for screen readers.
import { computed, nextTick, ref, watchEffect } from 'vue';
import type { CountsLoader } from '../application/counts';
import type { DraftState, SearchDraft } from '../application/draft';
import { parseTaxid } from '../domain/counts';
import { geneticCodeLabel } from '../domain/genetic-codes';
import { effectiveValue, fieldChoices, fieldList, knownDefault, writtenValue } from '../domain/parameters';
import { programById, type ParameterField } from '../domain/programs';
import { formatBytes, optionHelp } from './format';
import CountsLine from './CountsLine.vue';
import HelpTip from './HelpTip.vue';
import { refusal, refusedMark } from './refusal';

const props = defineProps<{ draft: SearchDraft; state: DraftState; field: ParameterField; inline?: boolean; counts?: CountsLoader | undefined }>();

/** The select's entry that opens the text box of a value outside NCBI's list. */
const OTHER = 'other';

const id = computed(() => `param-${props.field.flag.slice(1)}`);
const program = computed(() => programById(props.state.program));
const options = computed(() => props.state.description?.parameters);
const option = computed(() => options.value?.find((o) => o.flag === props.field.flag));
const values = computed(() => props.state.values[props.state.program]);
const value = computed(() => values.value[props.field.flag]);
const text = computed(() => (typeof value.value === 'string' ? value.value : ''));
const known = computed(() => knownDefault(program.value, props.field.flag, values.value, options.value));
const fallback = computed(() => known.value?.defaultValue);
/** The field of the same row that names the counts by a taxonomy ID, where the engine describes it (WW2). */
const taxidField = computed(() => {
  const flag = props.field.file?.taxid?.flag;
  return flag === undefined || options.value?.some((o) => o.flag === flag) === false ? undefined : program.value.sections.flatMap((s) => s.fields).find((f) => f.flag === flag);
});
const taxidValue = computed(() => (taxidField.value === undefined ? undefined : values.value[taxidField.value.flag]));
const taxidText = computed(() => (typeof taxidValue.value === 'string' ? taxidValue.value : ''));
/** The taxonomy ID that is typed, if it is one. */
const taxid = computed(() => parseTaxid(taxidText.value));
const changed = computed(
  () =>
    writtenValue(props.field, value.value, known.value) !== undefined ||
    (taxidField.value !== undefined &&
      writtenValue(taxidField.value, taxidValue.value, knownDefault(program.value, taxidField.value.flag, values.value, options.value)) !== undefined),
);
const choices = computed(() => (props.field.kind === 'boolean' ? ['true', 'false'] : fieldChoices(props.field, option.value)));
/** The radio that is checked: the value set, or else the engine's default. */
const checked = computed(() => effectiveValue(props.field.flag, values.value, options.value));
/** The value that a select or a checkbox shows: the value set, or else the default ('' when unknown). */
const shown = computed(() => (text.value.trim() !== '' ? text.value.trim() : (fallback.value ?? '')));

/** A true/false option with a known default is a checkbox (NCBI's "Mask for lookup table only"). */
const asCheckbox = computed(() => props.field.kind === 'flag' || (props.field.kind === 'boolean' && (fallback.value === 'true' || fallback.value === 'false')));
const boxChecked = computed(() => (props.field.kind === 'flag' ? value.value === true : shown.value === 'true'));
/** A low complexity filter is on unless the search uses `no` (NCBI's "Low complexity regions" checkbox). */
const filterOn = computed(() => shown.value !== '' && shown.value !== 'no');
/** The filter's text box stays while it has the focus, also when emptying it turns the filter off. */
const filterEditing = ref(false);
/** A file option is on while its box is checked: no file chosen yet (`true`), or the file's name. */
const fileOn = computed(
  () => props.field.kind === 'file' && (value.value === true || text.value.trim() !== '' || taxidValue.value === true || taxidText.value.trim() !== ''),
);
/** The row names its counts by a taxonomy ID (the default) or by a file chosen beside it. */
const countsBy = computed(() => (taxidField.value === undefined || value.value === true || text.value.trim() !== '' ? 'file' : 'taxid'));
/** The option that the row writes: -window_masker_taxid while the counts are named by a taxonomy ID. */
const shownFlag = computed(() => (taxidField.value !== undefined && countsBy.value === 'taxid' ? taxidField.value.flag : props.field.flag));
/** The file chosen for the option, which the run reads under the option's value. */
const chosenFile = computed(() => props.state.files[props.state.program][props.field.flag]);
const fileInput = ref<HTMLInputElement | null>(null);
/** A checkbox of a row carries its name after the box, as NCBI's checkboxes of a row do. */
const boxInRow = computed(() => props.inline === true && (asCheckbox.value || props.field.kind === 'filter' || props.field.kind === 'file'));
/** A choice whose default is one of its choices shows it; otherwise "Default" stands for it. */
const defaultChosen = computed(() => fallback.value !== undefined && choices.value.includes(fallback.value));

const list = computed(() => (props.field.kind === 'text' ? fieldList(props.field, values.value, options.value) : []));
/** The text box beside the select: chosen with "Other", or holding a value outside the list. */
const other = ref(false);
watchEffect(() => {
  if (value.value === undefined) other.value = false;
  else if (shown.value !== '' && list.value.length > 0 && !list.value.includes(shown.value)) other.value = true;
});
const otherInput = ref<HTMLInputElement | null>(null);

const choiceLabel = (choice: string) => props.field.choiceLabels?.[choice] ?? choice;
/** Whether the refusal of the argv names the field's option (WA2-AX-05). */
const refused = computed(
  () => refusal(props.state)?.flags.has(props.field.flag) === true || (taxidField.value !== undefined && refusal(props.state)?.flags.has(taxidField.value.flag) === true),
);
/**
 * The refused field's marks, on the control that holds its value: the text box after "Other" or of
 * a filter's arguments where it is shown (`on` is false for the select or the checkbox beside it).
 */
const mark = (on = true) => refusedMark(refused.value && on);

/** The engine's help, behind the "?" (HelpTip). */
const help = computed(() => optionHelp(option.value?.help));

const defaultLabel = computed(() => {
  const engineDefault = option.value?.defaultValue;
  if (engineDefault === undefined) return 'Default';
  return props.field.kind === 'gencode' ? `Default (${geneticCodeLabel(Number(engineDefault))})` : `Default (${choiceLabel(engineDefault)})`;
});

const gencodes = computed((): readonly number[] => {
  const description = props.state.description;
  const list = props.field.flag === '-query_gencode' ? description?.queryGencodes : description?.subjectGencodes;
  return list ?? [];
});

function set(next: string | boolean): void {
  props.draft.setField(props.field.flag, next);
}

/** A value chosen in a select or a checkbox: the default empties the field (it follows the task again). */
function choose(next: string): void {
  set(next === fallback.value ? '' : next);
}

/** Turning a filter on gives it its default arguments, or `yes` where the default is `no`; off is `no`. */
function toggleFilter(on: boolean): void {
  const filters = fallback.value !== undefined && fallback.value !== 'no';
  if (on) set(filters ? '' : 'yes');
  else set(filters ? 'no' : '');
}

/** The box of a file option: checked, it waits for a file (the run cannot start without one); unchecked, it lets the file go. */
function toggleFile(on: boolean): void {
  if (taxidField.value !== undefined) props.draft.setRepeats(on ? 'taxid' : 'off');
  else if (on) props.draft.setField(props.field.flag, true);
  else props.draft.setFile(props.field.flag, undefined);
}

/** The select of a row with a taxonomy ID: how the counts are named (a chosen file is let go). */
function chooseCountsBy(next: string): void {
  props.draft.setRepeats(next === 'file' ? 'file' : 'taxid');
}

/** Only the digits (and a minus sign) of a taxonomy ID are typed; the engine reads the rest as it is. */
function onTaxid(event: Event): void {
  const input = event.target as HTMLInputElement;
  props.draft.setTaxid(input.value.replace(/[^0-9-]/g, ''));
}

function onFile(event: Event): void {
  const input = event.target as HTMLInputElement;
  const file = input.files?.[0];
  if (file !== undefined) props.draft.setFile(props.field.flag, file);
  // The same file can be chosen again after the box is unchecked.
  input.value = '';
}

function onList(next: string): void {
  if (next === OTHER) {
    other.value = true;
    void nextTick(() => otherInput.value?.focus());
    return;
  }
  other.value = false;
  choose(next);
}
</script>

<template>
  <div
    :class="[inline ? 'param-inline' : 'param-row', { changed, 'param-check': boxInRow }]"
    :data-kind="field.kind"
    :data-testid="`${id}-field`"
  >
    <component
      :is="field.kind === 'radio' ? 'span' : 'label'"
      v-if="!boxInRow"
      :id="`${id}-label`"
      class="param-label"
      :for="field.kind === 'radio' ? undefined : id"
    >
      <span v-if="changed" class="changed-mark" aria-hidden="true">♦</span>
      {{ field.label }} <span class="flag">{{ field.flag }}</span>
      <span v-if="changed" class="visually-hidden">, changed from the default</span>
    </component>
    <div class="param-control">
      <div v-if="field.kind === 'radio'" class="radios" role="radiogroup" :aria-labelledby="`${id}-label`" v-bind="mark()" :data-testid="id">
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
      <template v-else-if="field.kind === 'filter'">
        <input
          :id="id"
          type="checkbox"
          :checked="filterOn"
          v-bind="mark(!(filterOn || filterEditing))"
          :data-testid="id"
          @change="toggleFilter(($event.target as HTMLInputElement).checked)"
        />
        <label v-if="inline" :id="`${id}-label`" class="param-label" :for="id">
          <span v-if="changed" class="changed-mark" aria-hidden="true">♦</span>
          {{ field.label }} <span class="flag">{{ field.flag }}</span>
          <span v-if="changed" class="visually-hidden">, changed from the default</span>
        </label>
        <!-- The filter's arguments (the CLI's value), as NCBI's "Species-specific repeats" has its select after the box. -->
        <template v-if="filterOn || filterEditing">
          <input
            :id="`${id}-value`"
            type="text"
            autocomplete="off"
            spellcheck="false"
            :value="value === undefined ? (fallback ?? '') : text"
            :placeholder="fallback === 'no' ? '' : (fallback ?? '')"
            :aria-labelledby="`${id}-label ${id}-value-name`"
            v-bind="mark()"
            :data-testid="`${id}-value`"
            @input="set(($event.target as HTMLInputElement).value)"
            @focus="filterEditing = true"
            @blur="filterEditing = false"
          />
          <span :id="`${id}-value-name`" class="visually-hidden">arguments</span>
        </template>
      </template>
      <template v-else-if="field.kind === 'file'">
        <input
          :id="id"
          type="checkbox"
          :checked="fileOn"
          v-bind="mark(!fileOn)"
          :data-testid="id"
          @change="toggleFile(($event.target as HTMLInputElement).checked)"
        />
        <label v-if="inline" :id="`${id}-label`" class="param-label" :for="id">
          <span v-if="changed" class="changed-mark" aria-hidden="true">♦</span>
          {{ field.label }} <span class="flag">{{ shownFlag }}</span>
          <span v-if="changed" class="visually-hidden">, changed from the default</span>
        </label>
        <!-- How the counts are named, as NCBI's "Species-specific repeats" has its select after the box (WW2): a taxonomy ID, or a file of the user's. -->
        <span v-if="fileOn && taxidField" class="select-clip">
          <select :value="countsBy" aria-label="Name the counts by" :data-testid="`${id}-by`" @change="chooseCountsBy(($event.target as HTMLSelectElement).value)">
            <option value="taxid">Taxonomy ID</option>
            <option value="file">Counts file</option>
          </select>
        </span>
        <template v-if="fileOn && countsBy === 'taxid'">
          <input
            :id="`${id}-taxid`"
            type="text"
            inputmode="numeric"
            autocomplete="off"
            spellcheck="false"
            placeholder="9606"
            aria-label="Taxonomy ID"
            :value="taxidText"
            v-bind="mark()"
            :data-testid="`${id}-taxid`"
            @input="onTaxid"
          />
          <span v-if="taxid === undefined" class="file-missing" role="status" :data-testid="`${id}-missing`">{{ field.file?.taxid?.missing }}</span>
        </template>
        <!-- The file of the option, chosen beside the box as NCBI's "Species-specific repeats" has its select. -->
        <span v-else-if="fileOn" class="file-pick">
          <button type="button" :data-testid="`${id}-choose`" @click="fileInput?.click()">{{ chosenFile ? 'Change file…' : 'Choose file…' }}</button>
          <!-- Opened by the button: out of the accessibility tree, as the inputs' file field (WA2-AX-13). -->
          <input
            ref="fileInput"
            class="visually-hidden"
            type="file"
            tabindex="-1"
            aria-hidden="true"
            :aria-label="field.file?.label"
            :data-testid="`${id}-file`"
            @change="onFile"
          />
          <span v-if="chosenFile" class="file-name" :data-testid="`${id}-name`">
            {{ field.file?.label }}: {{ chosenFile.name }} ({{ formatBytes(chosenFile.size) }})
          </span>
          <span v-else class="file-missing" role="status" :data-testid="`${id}-missing`">{{ field.file?.missing }}</span>
        </span>
      </template>
      <template v-else-if="asCheckbox">
        <input
          :id="id"
          type="checkbox"
          :checked="boxChecked"
          v-bind="mark()"
          :data-testid="id"
          @change="
            field.kind === 'flag'
              ? set(($event.target as HTMLInputElement).checked)
              : choose(($event.target as HTMLInputElement).checked ? 'true' : 'false')
          "
        />
        <!-- NCBI's checkboxes of a row ("Mask") carry their names after the box. -->
        <label v-if="inline" :id="`${id}-label`" class="param-label" :for="id">
          <span v-if="changed" class="changed-mark" aria-hidden="true">♦</span>
          {{ field.label }} <span class="flag">{{ field.flag }}</span>
          <span v-if="changed" class="visually-hidden">, changed from the default</span>
        </label>
      </template>
      <template v-else-if="list.length > 0">
        <span class="select-clip">
          <select
            :id="id"
            :value="other ? OTHER : shown"
            v-bind="mark(!other)"
            :data-testid="id"
            @change="onList(($event.target as HTMLSelectElement).value)"
          >
            <option v-if="shown === ''" value="">Default</option>
            <option v-for="item in list" :key="item" :value="item">{{ item }}</option>
            <option :value="OTHER">Other</option>
          </select>
        </span>
        <input
          v-if="other"
          :id="`${id}-other`"
          ref="otherInput"
          type="text"
          :inputmode="field.inputMode ?? 'text'"
          autocomplete="off"
          spellcheck="false"
          :value="text"
          :placeholder="fallback ?? ''"
          :aria-labelledby="`${id}-label`"
          v-bind="mark()"
          :data-testid="`${id}-other`"
          @input="set(($event.target as HTMLInputElement).value)"
        />
      </template>
      <template v-else-if="field.kind === 'text'">
        <input
          :id="id"
          type="text"
          :inputmode="field.inputMode ?? 'text'"
          autocomplete="off"
          spellcheck="false"
          :value="value === undefined ? (fallback ?? '') : text"
          :placeholder="fallback ?? ''"
          v-bind="mark()"
          :data-testid="id"
          @input="set(($event.target as HTMLInputElement).value)"
        />
      </template>
      <!-- WebKit draws a select's chosen label at its whole length beyond the select's box, which
           made a phone's page wider than its screen (S13 screen review H2): the span clips it. -->
      <span v-else class="select-clip">
        <select
          v-if="field.kind === 'gencode'"
          :id="id"
          :value="text"
          v-bind="mark()"
          :data-testid="id"
          @change="set(($event.target as HTMLSelectElement).value)"
        >
          <option value="">{{ defaultLabel }}</option>
          <option v-for="code in gencodes" :key="code" :value="String(code)">{{ geneticCodeLabel(code) }}</option>
        </select>
        <select
          v-else-if="defaultChosen"
          :id="id"
          :value="shown"
          v-bind="mark()"
          :data-testid="id"
          @change="choose(($event.target as HTMLSelectElement).value)"
        >
          <option v-for="choice in choices" :key="choice" :value="choice">{{ choiceLabel(choice) }}</option>
        </select>
        <select v-else :id="id" :value="text" v-bind="mark()" :data-testid="id" @change="set(($event.target as HTMLSelectElement).value)">
          <option value="">{{ defaultLabel }}</option>
          <option v-for="choice in choices" :key="choice" :value="choice">{{ choiceLabel(choice) }}</option>
        </select>
      </span>
      <span v-if="field.hint && (field.kind !== 'file' || (fileOn && countsBy === 'file'))" class="param-hint">{{ field.hint }}</span>
      <!-- The engine's help waits behind "?" (stop-describing-ui: the control shows, the help waits). -->
      <HelpTip v-if="help" :id="id" :label="`About ${field.flag}`" :text="help" />
      <!-- The line that downloads the counts of the taxonomy ID, under the row's controls (WW2). -->
      <CountsLine v-if="field.kind === 'file' && fileOn && countsBy === 'taxid' && taxid !== undefined && counts" :counts="counts" :taxid="taxid" :focus-id="`${id}-taxid`" />
    </div>
  </div>
</template>
