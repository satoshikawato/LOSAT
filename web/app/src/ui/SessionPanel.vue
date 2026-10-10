<script setup lang="ts">
// Session files in the aside (S15 item 4, design §12.2): save the completed runs of this working
// session, with the tray's candidates and notes when chosen (S15 decision 1: checked by default),
// and open a saved file, whose runs come back without a search. Everything acts through the
// session (application/session.ts); the names and messages are shown as text.
import { computed, ref } from 'vue';
import type { RunView } from '../application/coordinator';
import type { Session } from '../application/session';
import { formatCounted } from './format';
import { useStore } from './useStore';

const props = defineProps<{ session: Session; runs: readonly RunView[] }>();
const state = useStore(props.session.state);
const includeCandidates = ref(true);
const fileInput = ref<HTMLInputElement>();
const completed = computed(() => props.runs.filter((run) => run.status === 'completed').length);
const others = computed(() => props.runs.length - completed.value);
const busy = computed(() => state.value.busy !== undefined);

async function save(): Promise<void> {
  if (busy.value || completed.value === 0) return;
  await props.session.save({ includeCandidates: includeCandidates.value });
}

async function open(event: Event): Promise<void> {
  const input = event.target as HTMLInputElement;
  const file = input.files?.[0];
  // The same file can be chosen again (after a refusal, or to load it twice).
  input.value = '';
  if (file !== undefined) await props.session.load(file);
}
</script>

<template>
  <section class="session-panel" data-testid="session-panel">
    <h2>Session</h2>
    <p class="hint">
      A session file keeps the completed runs (what was searched, the outputs, the HSP records and the warnings), not the input FASTA files.
      Opening it shows the runs again without searching.
    </p>
    <label class="check">
      <input v-model="includeCandidates" type="checkbox" data-testid="session-include-candidates" />
      Include candidates and notes
    </label>
    <div class="session-actions">
      <button type="button" :disabled="busy || completed === 0" data-testid="session-save" @click="save">Save session</button>
      <button type="button" :disabled="busy" data-testid="session-open" @click="fileInput?.click()">Open session…</button>
      <input
        ref="fileInput"
        class="visually-hidden"
        type="file"
        accept=".gz,application/gzip"
        tabindex="-1"
        aria-label="Open a session file"
        data-testid="session-file"
        @change="open"
      />
    </div>
    <p class="hint" data-testid="session-save-note">
      {{ completed === 0 ? 'No completed run to save yet.' : `${formatCounted(completed, 'completed run')} will be saved.` }}
      Queued, running, cancelled and failed runs are never saved<template v-if="others > 0"> ({{ others }} here)</template>.
    </p>
    <div aria-live="polite">
      <p v-if="state.busy" class="muted" data-testid="session-busy">
        {{ state.busy === 'saving' ? 'Saving the session…' : 'Opening the session file…' }}
      </p>
      <p v-else-if="state.message" :class="state.message.error ? 'error' : 'ok'" data-testid="session-message" :data-error="state.message.error">
        {{ state.message.text }}
      </p>
    </div>
  </section>
</template>
