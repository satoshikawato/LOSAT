<script setup lang="ts">
// A run loaded from a session file, in Run details (S15 items 4 and 6, design §12.2, REQ-23): where
// it comes from and that it was not searched again; and per role, whether its original FASTA is
// attached, what cannot be done without it, and the explicit choice of the original files, which
// the session attaches only if they make the input that the run searched. A file of the same name
// elsewhere (the search form, another run) is never attached by itself.
import { computed } from 'vue';
import { attachKey, type Session } from '../application/session';
import type { InputRole } from '../domain/programs';
import { formatBytes, formatCount, formatCounted, formatDateTime } from './format';
import { useStore } from './useStore';

const props = defineProps<{ session: Session; runId: string }>();
const app = useStore(props.session.runs);
const sessionState = useStore(props.session.state);
const run = computed(() => app.value.runs.find((view) => view.snapshot.runId === props.runId));
const origin = computed(() => run.value?.fromSession);
const ROLES: readonly InputRole[] = ['query', 'subject'];
/** The hidden file inputs that the buttons open, by role. */
const pickers = new Map<InputRole, HTMLInputElement>();
const keepPicker = (role: InputRole, element: unknown) => {
  if (element instanceof HTMLInputElement) pickers.set(role, element);
};

const attaching = (role: InputRole) => sessionState.value.attaching.get(attachKey(props.runId, role));

async function choose(role: InputRole, event: Event): Promise<void> {
  const input = event.target as HTMLInputElement;
  const files = [...(input.files ?? [])];
  // The same files can be chosen again after a refusal.
  input.value = '';
  if (files.length > 0) await props.session.attach(props.runId, role, files);
}

/** The run's files of a role, as the session file recorded them: name, records, and those left out. */
function sourcesText(role: InputRole): string {
  const sources = origin.value?.inputs[role].sources ?? [];
  return sources
    .map((source) => {
      const left = source.excluded.length === 0 ? '' : `, ${formatCount(source.excluded.length)} left out`;
      return `${source.name} (${formatCounted(source.records, 'record')}${left}, ${formatBytes(source.size)})`;
    })
    .join(', then ');
}
</script>

<template>
  <section v-if="run && origin" class="loaded-run" data-testid="run-origin">
    <h3>Loaded from a session file</h3>
    <p data-testid="run-origin-file">
      From {{ origin.fileName }}, run {{ origin.number }} there (saved {{ formatDateTime(origin.savedAt) }}). It was not searched again: its
      outputs, HSP records and warnings are those that the file holds.
    </p>
    <div
      v-for="role in ROLES"
      :key="role"
      class="loaded-input"
      :data-testid="`run-original-${role}`"
      :data-attached="run.attached?.[role] !== undefined"
    >
      <h4>Original {{ role }} FASTA</h4>
      <p v-if="run.attached?.[role]" class="ok" :data-testid="`run-original-${role}-attached`">
        Attached: {{ run.attached[role]!.fileNames.join(', ') }}. Its records and SHA-256 match the {{ role }} that this run searched, so sequences
        can be extracted from this run's candidates.
      </p>
      <template v-else>
        <p class="notice" :data-testid="`run-original-${role}-missing`">
          Not attached. A session file does not hold the input FASTA, so sequences cannot be extracted from this run's candidates (hit
          regions, flanks or complete sequences) and the input FASTA of this run cannot be downloaded until the original file is chosen
          again. The outputs, HSP records and aligned sequences need no original.
        </p>
        <p class="hint">
          The run's {{ role }}: {{ sourcesText(role) }}. Files are attached only if their records and SHA-256 match it, with the same records
          left out.
        </p>
        <div class="loaded-actions">
          <button type="button" :disabled="attaching(role)?.busy === true" :data-testid="`run-attach-${role}`" @click="pickers.get(role)?.click()">
            Choose the original {{ role }} FASTA…
          </button>
          <input
            :ref="(element) => keepPicker(role, element)"
            class="visually-hidden"
            type="file"
            multiple
            tabindex="-1"
            :aria-label="`Choose the original ${role} FASTA of this run`"
            :data-testid="`run-attach-${role}-files`"
            @change="choose(role, $event)"
          />
        </div>
      </template>
      <p v-if="attaching(role)?.busy" class="muted" :data-testid="`run-attach-${role}-busy`">Checking the files…</p>
      <p v-else-if="attaching(role)?.message" class="error" role="alert" :data-testid="`run-attach-${role}-message`">
        {{ attaching(role)!.message }}
      </p>
    </div>
  </section>
</template>
