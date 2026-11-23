<script setup lang="ts">
// The compatibility outputs of a run, byte for byte as the engine wrote them, with the CLI
// command of each format and their export (plan §5.8; the view filters never apply here): the
// formats of every run, and the output of the run's own -outfmt (W10, W11). A binary output
// (outfmt 9, and the ZIP archive of the files of 13 and 14) is exported, not shown. An output that cannot be read, or an export that is not
// saved, says why (application/outputs.ts), and so does an output that the run could not write, which is empty, with
// where its error is (`ResultSetRef.failed`; line WQA's sweep, D7, 2026-11-19 and 2026-11-20). After the formats, the files that the run made besides
// its outputs (the file of -logfile; line RG E1w-3b) are a tab of their own (RunFiles). A run that printed the command line's help or
// version in place of a search has that text as its one output, "Help" or "Version" (PrintedText; line WQA's D14, 2026-11-20).
import { computed, ref, useId, watch } from 'vue';
import { CUSTOM_OUTPUT, type RunOutput } from 'losat/domain';
import type { Coordinator, RunView } from '../application/coordinator';
import { OutputReader, OutputSaver } from '../application/outputs';
import { failedOutput, keptOutputs } from '../application/results';
import { commandWords } from '../domain/argv';
import { fromText } from '../domain/input-parts';
import { isBinaryOutput, outfmtOf } from '../domain/output-format';
import CommandText from './CommandText.vue';
import { placeParts } from './message-place-text';
import PlaceWords from './PlaceWords.vue';
import OptionName from './OptionName.vue';
import OutputText from './OutputText.vue';
import PrintedText from './PrintedText.vue';
import RunFiles from './RunFiles.vue';
import TabList from './TabList.vue';
import { useStore } from './useStore';

/**
 * The Outputs of a run: those of a completed run, and for a failed run the output that the CLI leaves before its error
 * (line WQA's sweep, D15, 2026-11-20), where the results page shows its `-logfile` file and warnings apart (`files` false).
 */
const props = withDefaults(defineProps<{ coordinator: Coordinator; run: RunView; files?: boolean }>(), { files: true });
/** The tab of the run's files, after the formats. */
const FILES = 'run-files';
/** The tab of the text that the run printed in place of a search, its one output. */
const PRINTED = 'printed';
type Tab = RunOutput | typeof FILES | typeof PRINTED;
const chosen = ref<RunOutput>(6);
/** The tab of the run's files is chosen (the output chosen before stays read). */
const filesChosen = ref(false);
/** What the run printed in place of a search (the command line's help or version), if it did. */
const printed = computed(() => props.run.snapshot.prints);
/** The outputs that the run kept: all of a completed run's, the argv's own for a failed run (`keptOutputs`); none of a run that printed. */
const outputs = computed(() => (printed.value === undefined ? keptOutputs(props.run) : []));
const runFiles = computed(() => (props.files ? (props.run.record.runFiles ?? []) : []));
/** The output shown: the one chosen, or outfmt 6 for a run that does not keep it, else the first that it kept. */
const format = computed(() => (outputs.value.includes(chosen.value) ? chosen.value : outputs.value.includes(6) ? 6 : (outputs.value[0] ?? 6)));
const binary = computed(() => isBinaryOutput(props.run.snapshot.argv, format.value));
/** Where the error of the output shown is, in words, when the run could not write it, as the engine said (it is empty); undefined where it was written. */
const notWritten = computed(() => {
  const error = failedOutput(props.run, format.value);
  return error === undefined ? undefined : placeParts(error, props.run);
});
/** The run's files are shown, in place of an output. */
const filesShown = computed(() => filesChosen.value && runFiles.value.length > 0);
/** The formats, then the run's files: tabs of one panel, the output or the files shown (WA2-AX-12). */
const ids = useId();
const shown = computed({
  get: (): Tab => (filesShown.value ? FILES : printed.value !== undefined ? PRINTED : format.value),
  set: (value: Tab) => {
    filesChosen.value = value === FILES;
    if (value !== FILES && value !== PRINTED) chosen.value = value;
  },
});
const tabs = computed(() => [
  ...(printed.value === undefined
    ? []
    : [
        {
          value: PRINTED as Tab,
          label: printed.value === 'help' ? 'Help' : 'Version',
          testid: 'format-printed',
          id: `${ids}-tab-${PRINTED}`,
          panel: `${ids}-panel`,
        },
      ]),
  ...outputs.value.map((f) => ({
    value: f as Tab,
    label: `outfmt ${outfmtOf(props.run.snapshot.argv, f)}`,
    testid: f === CUSTOM_OUTPUT ? 'format-custom' : `format-${f}`,
    id: `${ids}-tab-${f}`,
    panel: `${ids}-panel`,
    ...(f === CUSTOM_OUTPUT ? { class: 'custom-output' } : {}),
  })),
  ...(runFiles.value.length === 0
    ? []
    : [
        {
          value: FILES as Tab,
          label: `-logfile ${runFiles.value.map((file) => file.name).join(', ')}`,
          testid: 'format-run-files',
          id: `${ids}-tab-${FILES}`,
          panel: `${ids}-panel`,
          class: 'custom-output',
        },
      ]),
]);
const reader = new OutputReader(props.coordinator);
/** The query and subject that are files of LOSAT's, which the command names: not the user's files as given (line WQA's D22, 2026-11-23). */
const madeInputs = computed(() =>
  (['query', 'subject'] as const).flatMap((role) => {
    const input = props.run.snapshot[role];
    return input.from === undefined ? [] : [{ role, name: input.name, from: fromText(input.from) }];
  }),
);
const inputSaver = new OutputSaver(props.coordinator);
const inputProblem = useStore(inputSaver.problem);
const output = useStore(reader.state);
// The reader keeps the latest choice: a later one may be read first.
// A run that printed a text in place of a search has no output to read.
watch(
  [() => props.run, format],
  ([run, fmt]) => {
    if (run.snapshot.prints === undefined) void reader.show(run, fmt, !isBinaryOutput(run.snapshot.argv, fmt));
  },
  { immediate: true },
);
</script>

<template>
  <div class="outputs-view" data-testid="results-outputs">
    <!-- Chosen by Enter or Space: each output shown is read from storage, the largest first (WA2-AX-12). -->
    <TabList v-model="shown" label="Output format" :tabs="tabs" manual />
    <div v-if="filesShown" :id="`${ids}-panel`" role="tabpanel" :aria-labelledby="`${ids}-tab-${FILES}`" class="outputs-panel">
      <RunFiles :key="run.snapshot.runId" :coordinator="coordinator" :run="run" />
    </div>
    <div v-else-if="printed !== undefined" :id="`${ids}-panel`" role="tabpanel" :aria-labelledby="`${ids}-tab-${PRINTED}`" class="outputs-panel">
      <PrintedText :key="run.snapshot.runId" :coordinator="coordinator" :run="run" />
    </div>
    <div v-else :id="`${ids}-panel`" role="tabpanel" :aria-labelledby="`${ids}-tab-${format}`" class="outputs-panel">
      <p class="command">
        <!-- The command as Run details shows it: a line does not break after an option's leading hyphen (line WQA's sweep, D7, 2026-11-20). -->
        Command: <code data-testid="result-command"><CommandText :words="commandWords(run.snapshot.argv, format)" /></code>
      </p>
      <!-- Where the command names a file that LOSAT made, the user needs it to run the command (line WQA's D22, 2026-11-23). -->
      <p v-for="made in madeInputs" :key="made.role" class="command" :data-testid="`outputs-${made.role}-file`">
        <OptionName :name="`-${made.role}`" /> {{ made.name }} is {{ made.from }}.
        <button type="button" class="link" :data-testid="`outputs-${made.role}-download`" @click="inputSaver.saveInput(run, made.role)">Download {{ made.name }}</button>
      </p>
      <p v-if="inputProblem" class="results-message error" role="alert" data-testid="outputs-input-problem">{{ inputProblem }}</p>
      <p v-if="notWritten !== undefined" class="results-message" data-testid="outputs-not-written">
        outfmt {{ outfmtOf(run.snapshot.argv, format) }} could not be written: its error is <PlaceWords :parts="notWritten" />.
      </p>
      <OutputText
        v-if="!output.failed && !binary"
        :text="output.text"
        :label="`outfmt ${outfmtOf(run.snapshot.argv, format)} output`"
        data-testid="result-output"
        :data-shown="output.shown"
      />
      <!-- Next to Export, whose failure it also tells; where the text would be when the text could not be read. -->
      <p v-if="output.problem" class="results-message error" role="alert" data-testid="outputs-problem" :data-failed="output.failed">{{ output.problem }}</p>
      <button data-testid="export-output" @click="reader.export(run, format)">
        Export outfmt {{ outfmtOf(run.snapshot.argv, format) }}
      </button>
    </div>
  </div>
</template>
