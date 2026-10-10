<script setup lang="ts">
// How to reproduce a run (S15 instructions, item 5; design §12.3), made from the run's argv (the
// RunSnapshot), never from the search form: the LOSAT command of each output format, the NCBI
// BLAST+ command to compare it with (or why there is none), the notes, the input FASTA that the
// engine searched - saved under the argv's names, so that the commands run as written - and the
// run's settings file.
import { computed, ref } from 'vue';
import type { RunView } from '../application/coordinator';
import type { RunFiles } from '../application/run-files';
import type { OutputFormat } from '../domain/output-format';
import type { InputRole } from '../domain/programs';
import { commandNotes, inputRelation, losatCommand, NCBI_BLAST_VERSION, ncbiCommand, ncbiComparison } from '../domain/reproduce';
import CommandText from './CommandText.vue';
import { useStore } from './useStore';

const props = defineProps<{ run: RunView; formats: readonly OutputFormat[]; runFiles: RunFiles }>();
const state = useStore(props.runFiles.state);
const snapshot = computed(() => props.run.snapshot);
const comparison = computed(() => ncbiComparison(snapshot.value.argv));
const notes = computed(() => commandNotes(snapshot.value.argv, snapshot.value.query.sha256 === snapshot.value.subject.sha256));
const roles: readonly InputRole[] = ['query', 'subject'];
const inputs = computed(() =>
  roles.map((role) => {
    const input = snapshot.value[role];
    const relation = inputRelation(role, input.name, input.records.length, props.runFiles.inputParts(props.run, role));
    return { role, name: input.name, relation };
  }),
);
/** The message of the latest file saved from this run. */
const message = computed(() => (state.value.run?.runId === snapshot.value.runId ? state.value.run : undefined));
const busy = ref(false);

async function act(work: () => Promise<unknown>): Promise<void> {
  busy.value = true;
  try {
    await work();
  } finally {
    busy.value = false;
  }
}
</script>

<template>
  <section class="reproduce" data-testid="run-reproduce">
    <h3>Reproduce this run</h3>
    <p class="muted small">Made from the arguments fixed when the run was queued, not from the search form.</p>

    <h4>LOSAT command line</h4>
    <ul class="commands">
      <li v-for="format in formats" :key="format">
        outfmt {{ format }}: <code :data-testid="`run-command-${format}`"><CommandText :text="losatCommand(snapshot.argv, format)" /></code>
      </li>
    </ul>

    <h4>NCBI BLAST+ {{ NCBI_BLAST_VERSION }}, for comparison</h4>
    <template v-if="comparison.refused.length === 0">
      <ul class="commands">
        <li v-for="format in formats" :key="format">
          outfmt {{ format }}:
          <code :data-testid="`run-ncbi-command-${format}`"><CommandText :text="ncbiCommand(snapshot.argv, format)" /></code>
        </li>
      </ul>
    </template>
    <template v-else>
      <p v-for="(line, i) in comparison.refused" :key="i" class="notice" data-testid="run-ncbi-unavailable">{{ line }}</p>
    </template>
    <p v-for="(line, i) in comparison.exceptions" :key="`exception-${i}`" class="notice" data-testid="run-ncbi-exception">{{ line }}</p>
    <ul class="reproduce-notes muted small" data-testid="run-reproduce-notes">
      <li v-for="(line, i) in notes" :key="i">{{ line }}</li>
    </ul>

    <h4>The input FASTA of this run</h4>
    <p class="muted small">The bytes that the engine searched, under the names that the commands use.</p>
    <ul class="reproduce-inputs">
      <li v-for="input in inputs" :key="input.role" :data-testid="`run-input-file-${input.role}`">
        <button type="button" :data-testid="`run-input-save-${input.role}`" :disabled="busy" @click="act(() => runFiles.saveInput(run, input.role))">
          Save {{ input.name }}
        </button>
        <span class="small">{{ input.relation }}</span>
      </li>
    </ul>

    <h4>Settings</h4>
    <p class="reproduce-settings">
      <button type="button" data-testid="run-settings-save" :disabled="busy" @click="act(() => runFiles.saveRunSettings(run))">
        Save the settings of this run
      </button>
      <span class="muted small">The program, the options and the threads, without the inputs (a LOSAT Web file that "Load settings…" of the search form reads).</span>
    </p>
    <p v-if="message" role="status" :class="message.kind === 'error' ? 'error' : 'note'" data-testid="run-files-message" :data-kind="message.kind">
      {{ message.text }}
    </p>
  </section>
</template>
