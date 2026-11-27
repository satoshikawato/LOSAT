<script setup lang="ts">
// Shared local BLAST database picker (WQA D12b + RG-D4 W2): single owner reused for both
// Subject "Database" (-db) and BLASTN "Repeats database" (-filtering_db). Picks a local folder
// (webkitdirectory) with multi-file fallback, auto-discovers valid BLAST databases matching the
// required molecule type, and mounts them read-only in place without uploading or copying into
// Wasm memory.
import { computed, ref, useId } from 'vue';
import type { SequenceKind } from 'losat/domain';
import type { DatabasePickDraft } from '../application/draft';
import { chooseAgainMessage, databaseTaxonomy, defaultDatabase, matchingDatabases, noMatchingDatabaseMessage, rememberedPickKind } from '../domain/blast-db';
import { formatBytes } from './format';
import LosatLoader from './LosatLoader.vue';

const props = withDefaults(
  defineProps<{
    pick: DatabasePickDraft;
    moleculeKind: SequenceKind;
    contextLabel: string;
    testIdPrefix: string;
    selectLabel?: string;
  }>(),
  { selectLabel: 'Database' },
);

const emit = defineEmits<{
  files: [files: readonly File[]];
  select: [name: string];
  clear: [];
}>();

const selectId = useId();
const folderInput = ref<HTMLInputElement | null>(null);
const filesInput = ref<HTMLInputElement | null>(null);

const discovery = computed(() => props.pick.discovery);
const matches = computed(() => matchingDatabases(discovery.value, props.moleculeKind));
const selected = computed(() => matches.value.find((db) => db.name === props.pick.selected) ?? defaultDatabase(matches.value));
const mismatch = computed(() =>
  discovery.value !== undefined && !props.pick.scanning && matches.value.length === 0
    ? noMatchingDatabaseMessage(discovery.value, props.moleculeKind, props.contextLabel)
    : undefined,
);
const taxonomySummary = computed(() => {
  if (discovery.value === undefined || selected.value === undefined) return '';
  const tax = databaseTaxonomy(discovery.value, selected.value.name);
  const parts: string[] = [];
  if (tax.taxdb) parts.push('taxdb');
  if (tax.sqlite) parts.push('taxonomy4blast.sqlite3');
  return parts.join(', ');
});

function onPicked(event: Event): void {
  const input = event.target as HTMLInputElement;
  const picked = [...(input.files ?? [])];
  if (picked.length > 0) emit('files', picked);
  input.value = '';
}

function optionText(db: (typeof matches.value)[number]): string {
  const parts = [db.name];
  if (db.volumes !== undefined && db.volumes.length > 1) parts.push(`(${db.volumes.length} volumes)`);
  else if (db.entry === 'alias') parts.push('(alias)');
  if (db.title !== undefined && db.title !== db.name) parts.push(`— ${db.title}`);
  return parts.join(' ');
}
</script>

<template>
  <div class="db-picker" :data-testid="`${testIdPrefix}-picker`">
    <div class="db-picker-actions">
      <button type="button" :data-testid="`${testIdPrefix}-choose-folder`" @click="folderInput?.click()">
        {{ discovery ? 'Change folder…' : 'Choose folder…' }}
      </button>
      <button type="button" :data-testid="`${testIdPrefix}-choose-files`" @click="filesInput?.click()">
        Choose files…
      </button>
      <button
        v-if="discovery !== undefined || pick.remembered !== undefined"
        type="button"
        class="link"
        :data-testid="`${testIdPrefix}-clear`"
        @click="emit('clear')"
      >
        Clear
      </button>
      <!-- Hidden inputs opened by the buttons above: kept out of the accessibility tree (WA2-AX-13). -->
      <input
        ref="folderInput"
        class="visually-hidden"
        type="file"
        webkitdirectory
        multiple
        tabindex="-1"
        aria-hidden="true"
        :aria-label="`Choose ${contextLabel} database folder`"
        :data-testid="`${testIdPrefix}-folder-input`"
        @change="onPicked"
      />
      <input
        ref="filesInput"
        class="visually-hidden"
        type="file"
        multiple
        tabindex="-1"
        aria-hidden="true"
        :aria-label="`Choose ${contextLabel} database files`"
        :data-testid="`${testIdPrefix}-files-input`"
        @change="onPicked"
      />
    </div>

    <p v-if="pick.scanning" class="muted db-status" role="status" :data-testid="`${testIdPrefix}-scanning`">
      <LosatLoader />Scanning database files…
    </p>

    <p v-else-if="discovery === undefined && pick.remembered !== undefined" class="note db-reload" role="status" :data-testid="`${testIdPrefix}-reload`">
      Previous database: <code>{{ pick.remembered.name }}</code>. {{ chooseAgainMessage(rememberedPickKind(pick.remembered)) }}
    </p>

    <template v-else-if="discovery !== undefined">
      <div v-if="matches.length > 0" class="db-picker-found">
        <div v-if="matches.length > 1" class="db-select-row">
          <label :for="selectId">{{ selectLabel }}</label>
          <span class="select-clip">
            <select
              :id="selectId"
              :value="selected?.name ?? ''"
              :data-testid="`${testIdPrefix}-select`"
              @change="emit('select', ($event.target as HTMLSelectElement).value)"
            >
              <option v-for="db in matches" :key="db.name" :value="db.name">{{ optionText(db) }}</option>
            </select>
          </span>
        </div>
        <p v-else-if="selected" class="db-single" :data-testid="`${testIdPrefix}-selected`">
          <strong>{{ selectLabel }}:</strong> <code>{{ selected.name }}</code>
          <span v-if="selected.volumes && selected.volumes.length > 1" class="muted"> ({{ selected.volumes.length }} volumes)</span>
          <span v-else-if="selected.entry === 'alias'" class="muted"> (alias)</span>
          <span v-if="selected.title && selected.title !== selected.name" class="muted"> — {{ selected.title }}</span>
        </p>
        <p class="db-meta muted" :data-testid="`${testIdPrefix}-meta`">
          {{ discovery.files.length }} {{ discovery.files.length === 1 ? 'file' : 'files' }} ({{ formatBytes(discovery.totalBytes) }})<template
            v-if="taxonomySummary !== ''"
          >
            · taxonomy: {{ taxonomySummary }}</template
          >
        </p>
      </div>
      <p v-else-if="mismatch !== undefined" class="error db-mismatch" role="alert" :data-testid="`${testIdPrefix}-mismatch`">
        {{ mismatch }}
      </p>
    </template>

    <p class="db-privacy muted" :data-testid="`${testIdPrefix}-privacy`">
      LOSAT Web reads the database from this device; nothing is uploaded.
    </p>
  </div>
</template>
