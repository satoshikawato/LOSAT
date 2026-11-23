<script setup lang="ts">
// Results view of a WindowMasker tool run (WW3b DESIGN §3):
// Shows only "Outputs" and "Run details" (no Hits, no Dot Plot, no Supported badge).
// Outputs preview at most 512 KiB so large counts/FASTA outputs never grow page JS heap;
// binary outputs show "Binary, <size>". Download streams the full OPFS-backed File.
// Made counts files (`-out` of `-mk_counts` / `-convert`) offer "Use in BLASTN" and
// "Mask with these counts".
import { computed, ref, useId, watch } from 'vue';
import type { Coordinator, SessionCountsEntry, ToolRunView } from '../application/coordinator';
import type { ToolFileOrigin, ToolOutputEntry } from '../domain/run';
import { toolCommandWords } from '../domain/tool-argv';
import { decodeUtf8Preview, windowMaskerStageById } from '../domain/tools';
import CommandText from './CommandText.vue';
import { formatBytes, formatDateTime, formatDuration } from './format';
import OutputText from './OutputText.vue';
import TabList from './TabList.vue';

const PREVIEW_BYTES = 512 * 1024;

const props = defineProps<{
  coordinator: Coordinator;
  run: ToolRunView;
  takesValue?: ((flag: string) => boolean) | undefined;
}>();

const emit = defineEmits<{
  'use-in-blastn': [file: File];
  'mask-with-counts': [entry: SessionCountsEntry];
}>();

const ids = useId();
const section = ref<'outputs' | 'details'>('outputs');
const SECTION_TABS = [
  { value: 'outputs', label: 'Outputs', testid: 'results-view-outputs', id: `${ids}-outputs-tab`, panel: `${ids}-outputs` },
  { value: 'details', label: 'Run details', testid: 'results-view-details', id: `${ids}-details-tab`, panel: `${ids}-details` },
] as const;

interface OutputTab {
  readonly key: string;
  readonly label: string;
  readonly testid: string;
  readonly entry: ToolOutputEntry;
}

const outputTabs = computed<readonly OutputTab[]>(() => {
  const { stdout, stderr, files } = props.run.record;
  const list: OutputTab[] = [];
  if (stdout !== undefined && stdout.size > 0) {
    list.push({ key: stdout.key, label: 'Output', testid: 'tool-tab-stdout', entry: stdout });
  }
  (files ?? []).forEach((file, index) => {
    const label = file.option !== undefined ? `${file.option} ${file.name}` : file.name;
    list.push({ key: file.key, label, testid: `tool-tab-file-${index}`, entry: file });
  });
  if (stderr !== undefined && stderr.size > 0) {
    list.push({ key: stderr.key, label: 'Warnings', testid: 'tool-tab-stderr', entry: stderr });
  }
  return list;
});

const tabItems = computed(() =>
  outputTabs.value.map((item) => ({
    value: item.key,
    label: item.label,
    testid: item.testid,
    id: `${ids}-out-${item.key}`,
    panel: `${ids}-out-panel`,
  })),
);

const chosenKey = ref<string>('');
const activeTab = computed(() => outputTabs.value.find((tab) => tab.key === chosenKey.value) ?? outputTabs.value[0]);

const activeKey = computed({
  get: () => activeTab.value?.key ?? '',
  set: (value: string) => {
    chosenKey.value = value;
  },
});

const preview = ref<
  | { readonly state: 'loading' }
  | { readonly state: 'text'; readonly text: string; readonly shownBytes: number; readonly totalSize: number; readonly truncated: boolean }
  | { readonly state: 'binary'; readonly totalSize: number }
  | { readonly state: 'error'; readonly message: string }
>({ state: 'loading' });

const actionError = ref<string | undefined>();
let readToken = 0;

watch(
  [() => props.run.snapshot.runId, () => activeTab.value?.key],
  async ([runId, key]) => {
    actionError.value = undefined;
    const token = ++readToken;
    if (key === undefined) {
      preview.value = { state: 'text', text: '', shownBytes: 0, totalSize: 0, truncated: false };
      return;
    }
    preview.value = { state: 'loading' };
    try {
      const { bytes, totalSize, truncated } = await props.coordinator.readToolOutputPreview(runId, key, PREVIEW_BYTES);
      if (token !== readToken) return;
      const decoded = decodeUtf8Preview(bytes, truncated);
      if (decoded === undefined) {
        preview.value = { state: 'binary', totalSize };
      } else {
        preview.value = { state: 'text', text: decoded, shownBytes: bytes.length, totalSize, truncated };
      }
    } catch (error) {
      if (token !== readToken) return;
      preview.value = { state: 'error', message: error instanceof Error ? error.message : String(error) };
    }
  },
  { immediate: true },
);

function setActionError(error: unknown): void {
  actionError.value = error instanceof Error ? error.message : String(error);
}

async function downloadActive(): Promise<void> {
  const tab = activeTab.value;
  if (tab === undefined) return;
  actionError.value = undefined;
  try {
    await props.coordinator.exportToolOutput(props.run.snapshot.runId, tab.key);
  } catch (error) {
    setActionError(error);
  }
}

async function useInBlastn(): Promise<void> {
  const tab = activeTab.value;
  if (tab === undefined) return;
  actionError.value = undefined;
  try {
    const file = await props.coordinator.getToolOutputFile(props.run.snapshot.runId, tab.key, tab.entry.name);
    emit('use-in-blastn', file);
  } catch (error) {
    setActionError(error);
  }
}

function maskWithCounts(): void {
  const tab = activeTab.value;
  if (tab === undefined) return;
  emit('mask-with-counts', {
    runId: props.run.snapshot.runId,
    runNumber: props.run.snapshot.number,
    key: tab.key,
    name: tab.entry.name,
    size: tab.entry.size,
  });
}

const command = computed(() => toolCommandWords(props.run.snapshot.argv));
const stageLabel = computed(() => windowMaskerStageById(props.run.snapshot.stage).label);
const duration = computed(() => {
  const { startedAt, endedAt } = props.run.record;
  return startedAt !== undefined && endedAt !== undefined ? formatDuration(endedAt - startedAt) : '';
});

function originText(origin: ToolFileOrigin): string {
  switch (origin.kind) {
    case 'picked':
      return 'opened file';
    case 'pasted':
      return 'pasted text';
    case 'inflated':
      return `inflated from ${origin.from}`;
    case 'combined':
      return `combined from ${origin.parts.join(', ')}`;
    case 'made':
      return `made by run ${origin.runNumber}`;
  }
}
</script>

<template>
  <div class="tool-run-view" data-testid="tool-run-view" :data-run="run.snapshot.number">
    <p
      v-if="run.status === 'failed' && run.record.exitStatus !== undefined"
      class="results-message error"
      data-testid="tool-exit-status"
    >
      Exit status {{ run.record.exitStatus }}
    </p>

    <TabList v-model="section" class="results-tabs" label="Results view" :tabs="SECTION_TABS" />

    <div
      v-if="section === 'outputs'"
      :id="SECTION_TABS[0].panel"
      role="tabpanel"
      :aria-labelledby="SECTION_TABS[0].id"
      class="results-group outputs-view"
      data-testid="results-outputs"
    >
      <p class="command">
        Command: <code data-testid="result-command"><CommandText :words="command" :takes-value="takesValue" /></code>
      </p>

      <p v-if="outputTabs.length === 0" class="muted" data-testid="tool-output-empty">This run wrote no output.</p>
      <template v-else>
        <TabList v-model="activeKey" label="Tool outputs" :tabs="tabItems" manual />
        <div
          v-if="activeTab"
          :id="`${ids}-out-panel`"
          role="tabpanel"
          :aria-labelledby="`${ids}-out-${activeTab.key}`"
          class="outputs-panel"
        >
          <p v-if="preview.state === 'loading'" class="muted">Reading output…</p>
          <p v-else-if="preview.state === 'error'" class="results-message error" role="alert" data-testid="tool-output-error">
            {{ preview.message }}
          </p>
          <p v-else-if="preview.state === 'binary'" class="muted" data-testid="tool-output-binary">
            Binary, {{ formatBytes(preview.totalSize) }}
          </p>
          <template v-else>
            <OutputText
              :text="preview.text"
              :label="`${activeTab.label} output`"
              data-testid="tool-output-text"
              :data-shown="preview.shownBytes"
            />
            <p v-if="preview.truncated" class="muted small" data-testid="tool-output-truncated">
              Showing the first {{ formatBytes(preview.shownBytes) }} of {{ formatBytes(preview.totalSize) }}. Download to read the full file.
            </p>
          </template>

          <p v-if="actionError" class="results-message error" role="alert" data-testid="tool-action-error">
            {{ actionError }}
          </p>

          <div class="tool-output-actions">
            <button type="button" data-testid="tool-output-download" @click="downloadActive">
              Download {{ activeTab.entry.name }} ({{ formatBytes(activeTab.entry.size) }})
            </button>
            <template v-if="activeTab.entry.isCounts && run.status === 'completed'">
              <button type="button" class="primary-action" data-testid="tool-use-in-blastn" @click="useInBlastn">
                Use in BLASTN
              </button>
              <button type="button" data-testid="tool-mask-with-counts" @click="maskWithCounts">
                Mask with these counts
              </button>
            </template>
          </div>
        </div>
      </template>
    </div>

    <div
      v-else
      :id="SECTION_TABS[1].panel"
      role="tabpanel"
      :aria-labelledby="SECTION_TABS[1].id"
      class="results-group run-details"
      data-testid="run-details"
    >
      <div class="details-grid">
        <h3>Command (fixed when the run was submitted)</h3>
        <dl class="details-list" data-testid="run-snapshot">
          <dt>Run</dt>
          <dd>
            {{ run.snapshot.number }}<template v-if="run.snapshot.group">, group run {{ run.snapshot.group.position }} of {{ run.snapshot.group.size }}</template>
          </dd>
          <dt>Tool</dt>
          <dd>WindowMasker · {{ stageLabel }}</dd>
          <dt>Command</dt>
          <dd><code data-testid="run-argv"><CommandText :words="command" :takes-value="takesValue" /></code></dd>
          <dt>Files of the run</dt>
          <dd data-testid="tool-details-files">
            <span v-if="run.snapshot.directory.length === 0" class="muted">none</span>
            <ul v-else class="wm-details-files">
              <li v-for="entry in run.snapshot.directory" :key="entry.name">
                {{ entry.name }} ({{ formatBytes(entry.size) }}, {{ originText(entry.origin) }})
              </li>
            </ul>
          </dd>
          <dt>Submitted</dt>
          <dd data-detail="queued">{{ formatDateTime(run.snapshot.queuedAt) }}</dd>
        </dl>

        <h3>Run</h3>
        <dl class="details-list" data-testid="run-record">
          <dt>Exit status</dt>
          <dd data-testid="tool-details-exit-status">{{ run.record.exitStatus }}</dd>
          <template v-if="run.record.startedAt !== undefined">
            <dt>Started</dt>
            <dd data-detail="started">{{ formatDateTime(run.record.startedAt) }}</dd>
          </template>
          <template v-if="duration">
            <dt>Took</dt>
            <dd>{{ duration }}</dd>
          </template>
        </dl>
      </div>
    </div>
  </div>
</template>
