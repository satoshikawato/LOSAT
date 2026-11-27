<script setup lang="ts">
// Results view of a WindowMasker tool run (WW3b DESIGN §3):
// Shows "Outputs" and "Run details" (no Hits, no Dot Plot) and the Supported mark for completed runs.
// Outputs preview at most 512 KiB so large counts/FASTA outputs never grow page JS heap;
// binary outputs show "Binary, <size>". Download streams the full OPFS-backed File.
// Made counts files (`-out` of `-mk_counts` / `-convert`) offer "Use in BLASTN" and
// "Mask with these counts"; mask data (a masking run's `-outfmt maskinfo_*` output, `maskDataName`) offers "Use as
// -mask_data in makeblastdb". A tool that writes its files into the run's directory (makeblastdb, line WW WW4b2) has
// them listed, each to download and all in one tar archive, and its database offers "Use as Database
// (-db)" and, for a nucleotide one, "Use as Repeats database (-filtering_db)".
import { computed, nextTick, ref, shallowRef, useId, watch } from 'vue';
import { activeSavedCountsTaxid, savedTabCopyOf, type CountsLoader, type SavedTabCopy } from '../application/counts';
import type { Coordinator, MadeDatabase, MaskDataHandover, SessionCountsEntry, ToolRunView } from '../application/coordinator';
import { primaryBlastDb } from '../domain/blast-db';
import { parseTaxid, type CountsMadeMeta } from '../domain/counts';
import type { ToolFileOrigin, ToolOutputEntry } from '../domain/run';
import { maskDataName, toolCommandWords } from '../domain/tool-argv';
import { decodeUtf8Preview, toolRunTitle, toolWritesDirectory } from '../domain/tools';
import type { Badge } from '../domain/verification';
import type { CountsCopy } from '../ports/counts';
import { announce } from './announcer';
import CommandText from './CommandText.vue';
import { formatBytes, formatDateTime, formatDuration } from './format';
import OutputText from './OutputText.vue';
import ReplaceCountsDialog from './ReplaceCountsDialog.vue';
import TabList from './TabList.vue';
import VerificationBadge from './VerificationBadge.vue';
import { useStore } from './useStore';

const PREVIEW_BYTES = 512 * 1024;

const props = defineProps<{
  coordinator: Coordinator;
  run: ToolRunView;
  counts?: CountsLoader | undefined;
  badge?: Badge | undefined;
  takesValue?: ((flag: string) => boolean) | undefined;
}>();

const emit = defineEmits<{
  'use-in-blastn': [file: File];
  'use-taxid-in-blastn': [taxid: number];
  'mask-with-counts': [entry: SessionCountsEntry];
  'use-as-database': [made: MadeDatabase];
  'use-as-repeats-database': [made: MadeDatabase];
  'use-as-mask-data': [data: MaskDataHandover];
}>();

const countsState = props.counts === undefined ? undefined : useStore(props.counts.state);

const ids = useId();
const section = defineModel<'outputs' | 'details'>('section', { default: 'outputs' });
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

/** Whether the run's files are those it wrote into its directory (listed apart from the tabs; makeblastdb's). */
const writesDirectory = computed(() => toolWritesDirectory(props.run.snapshot.tool));
const madeFiles = computed(() => (writesDirectory.value ? (props.run.record.files ?? []) : []));
/** The database among the files that the run made, if any (its alias file's, or else its volume's). */
const madeDatabase = computed(() =>
  props.run.status === 'completed' ? primaryBlastDb(madeFiles.value.map((file) => file.name)) : undefined,
);

const outputTabs = computed<readonly OutputTab[]>(() => {
  const { stdout, stderr, files } = props.run.record;
  const list: OutputTab[] = [];
  if (stdout !== undefined && stdout.size > 0) {
    list.push({ key: stdout.key, label: 'Output', testid: 'tool-tab-stdout', entry: stdout });
  }
  (files ?? []).forEach((file, index) => {
    // Of the files in the run's directory, those that an option named (-logfile) also have a tab.
    if (writesDirectory.value && file.option === undefined) return;
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
const saveTaxidText = ref('');
const saving = ref(false);
const saveError = ref<string | undefined>();
const saveStatus = ref<string | undefined>();
const savedByTab = ref<Readonly<Record<string, SavedTabCopy>>>({});
const saveButton = ref<HTMLButtonElement | null>(null);
const pendingReplace = shallowRef<
  | {
      readonly taxid: number;
      readonly existing: CountsCopy;
      readonly file: File;
      readonly made: CountsMadeMeta;
      readonly tabKey: string;
      readonly newFileName: string;
    }
  | undefined
>();
let readToken = 0;

const activeSavedTaxid = computed(() => {
  const key = activeTab.value?.key;
  if (key === undefined) return undefined;
  const saved = savedByTab.value[`${props.run.snapshot.runId}:${key}`];
  if (saved === undefined) return undefined;
  return activeSavedCountsTaxid(saved, countsState?.value.saved[saved.taxid]);
});

watch(
  () => props.run.snapshot.runId,
  () => {
    saveTaxidText.value = '';
  },
);

watch(
  [() => props.run.snapshot.runId, () => activeTab.value?.key],
  async ([runId, key]) => {
    actionError.value = undefined;
    saveError.value = undefined;
    saveStatus.value = undefined;
    pendingReplace.value = undefined;
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
  const savedTaxid = activeSavedTaxid.value;
  if (savedTaxid !== undefined) {
    emit('use-taxid-in-blastn', savedTaxid);
    return;
  }
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

async function saveForTaxid(): Promise<void> {
  const tab = activeTab.value;
  const counts = props.counts;
  if (tab === undefined || counts === undefined || saving.value) return;
  saveError.value = undefined;
  saveStatus.value = undefined;
  const taxid = parseTaxid(saveTaxidText.value);
  if (taxid === undefined) {
    saveError.value = 'Enter a taxonomy ID.';
    return;
  }
  saving.value = true;
  try {
    const file = await props.coordinator.getToolOutputFile(props.run.snapshot.runId, tab.key, tab.entry.name);
    const made: CountsMadeMeta = {
      run: props.run.snapshot.number,
      argv: props.run.snapshot.argv,
      inputs: props.run.snapshot.directory.map((entry) => ({ name: entry.name, size: entry.size })),
      engineBuild: props.run.record.engineBuild ?? '',
      madeAt: props.run.record.endedAt ?? Date.now(),
    };
    const result = await counts.saveMade(taxid, file, made, { replace: false });
    if (result.ok) {
      finishSave(tab.key, result.copy, result.kept);
      return;
    }
    if (result.reason === 'exists') {
      pendingReplace.value = { taxid, existing: result.existing, file, made, tabKey: tab.key, newFileName: tab.entry.name };
      return;
    }
    saveError.value = result.message;
    announce(result.message, 'alert');
  } catch (error) {
    saveError.value = error instanceof Error ? error.message : String(error);
    announce(saveError.value, 'alert');
  } finally {
    saving.value = false;
  }
}

async function confirmReplace(): Promise<void> {
  const pending = pendingReplace.value;
  const counts = props.counts;
  pendingReplace.value = undefined;
  if (pending === undefined || counts === undefined) return;
  saving.value = true;
  saveError.value = undefined;
  saveStatus.value = undefined;
  try {
    const result = await counts.saveMade(pending.taxid, pending.file, pending.made, { replace: true });
    if (result.ok) {
      finishSave(pending.tabKey, result.copy, result.kept);
    } else if (result.reason !== 'exists') {
      saveError.value = result.message;
      announce(result.message, 'alert');
    }
  } catch (error) {
    saveError.value = error instanceof Error ? error.message : String(error);
    announce(saveError.value, 'alert');
  } finally {
    saving.value = false;
    await nextTick();
    saveButton.value?.focus();
  }
}

async function cancelReplace(): Promise<void> {
  pendingReplace.value = undefined;
  await nextTick();
  saveButton.value?.focus();
}

function finishSave(tabKey: string, copy: CountsCopy, kept: boolean): void {
  const identity = savedTabCopyOf(copy);
  if (identity !== undefined) {
    savedByTab.value = { ...savedByTab.value, [`${props.run.snapshot.runId}:${tabKey}`]: identity };
  }
  const text = kept
    ? `Saved for taxonomy ID ${copy.taxid}.`
    : `Kept in this tab only for taxonomy ID ${copy.taxid}.`;
  saveStatus.value = text;
  announce(text);
}

async function downloadFile(entry: ToolOutputEntry): Promise<void> {
  actionError.value = undefined;
  try {
    await props.coordinator.exportToolOutput(props.run.snapshot.runId, entry.key);
  } catch (error) {
    setActionError(error);
  }
}

async function downloadAll(): Promise<void> {
  actionError.value = undefined;
  try {
    await props.coordinator.exportToolFiles(props.run.snapshot.runId);
  } catch (error) {
    setActionError(error);
  }
}

async function useMadeDatabase(as: 'use-as-database' | 'use-as-repeats-database'): Promise<void> {
  actionError.value = undefined;
  try {
    const made = await props.coordinator.madeDatabase(props.run.snapshot.runId);
    if (as === 'use-as-database') emit('use-as-database', made);
    else emit('use-as-repeats-database', made);
  } catch (error) {
    setActionError(error);
  }
}

/** The name under which the active output goes to makeblastdb's `-mask_data`, when it is mask data of a completed run. */
const maskDataOfTab = computed(() => {
  const tab = activeTab.value;
  return tab === undefined || props.run.status !== 'completed' ? undefined : maskDataName(props.run.snapshot, tab.entry);
});

async function useAsMaskData(): Promise<void> {
  const tab = activeTab.value;
  if (tab === undefined) return;
  actionError.value = undefined;
  try {
    emit('use-as-mask-data', await props.coordinator.maskData(props.run.snapshot.runId, tab.key));
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
const toolHeading = computed(() => toolRunTitle(props.run.snapshot.tool, props.run.snapshot.stage));
const duration = computed(() => {
  const { startedAt, endedAt } = props.run.record;
  return startedAt !== undefined && endedAt !== undefined ? formatDuration(endedAt - startedAt) : '';
});

function downloadLabel(tab: OutputTab): string {
  if (tab.entry.key === 'stdout') return 'output';
  if (tab.entry.key === 'stderr') return 'warnings';
  return tab.entry.name;
}

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
    case 'archive':
      return `from run ${origin.runNumber} (-outfmt 11)`;
    case 'saved':
      return `saved for taxonomy ID ${origin.taxid}`;
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

      <section v-if="madeFiles.length > 0" class="tool-made-files" data-testid="tool-made-files">
        <h3>Files made by the run</h3>
        <ul class="wm-file-list" data-testid="tool-made-files-list">
          <li v-for="(file, index) in madeFiles" :key="file.key" class="wm-file-item">
            <span class="source-name" :data-testid="`tool-made-file-name-${index}`">{{ file.name }}</span>
            <span class="muted">{{ formatBytes(file.size) }}</span>
            <button
              type="button"
              class="link"
              :aria-label="`Download ${file.name}`"
              :data-testid="`tool-made-file-download-${index}`"
              @click="downloadFile(file)"
            >
              Download
            </button>
          </li>
        </ul>
        <div class="tool-output-actions">
          <button type="button" data-testid="tool-made-files-download-all" @click="downloadAll">
            Download all ({{ madeFiles.length === 1 ? '1 file' : `${madeFiles.length} files` }}, .tar)
          </button>
          <template v-if="madeDatabase">
            <button type="button" class="primary-action" data-testid="tool-use-as-database" @click="useMadeDatabase('use-as-database')">
              Use as Database (-db)
            </button>
            <button
              v-if="madeDatabase.kind === 'nucleotide'"
              type="button"
              data-testid="tool-use-as-repeats-database"
              @click="useMadeDatabase('use-as-repeats-database')"
            >
              Use as Repeats database (-filtering_db)
            </button>
          </template>
        </div>
        <p v-if="actionError && outputTabs.length === 0" class="results-message error" role="alert" data-testid="tool-action-error">
          {{ actionError }}
        </p>
      </section>

      <p v-if="outputTabs.length === 0 && madeFiles.length === 0" class="muted" data-testid="tool-output-empty">This run wrote no output.</p>
      <template v-else-if="outputTabs.length > 0">
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
              Download {{ downloadLabel(activeTab) }} ({{ formatBytes(activeTab.entry.size) }})
            </button>
            <button v-if="maskDataOfTab" type="button" class="primary-action" data-testid="tool-use-as-mask-data" @click="useAsMaskData">
              Use as -mask_data in makeblastdb
            </button>
            <template v-if="activeTab.entry.isCounts && run.status === 'completed'">
              <button type="button" class="primary-action" data-testid="tool-use-in-blastn" @click="useInBlastn">
                Use in BLASTN
              </button>
              <button type="button" data-testid="tool-mask-with-counts" @click="maskWithCounts">
                Mask with these counts
              </button>
              <span v-if="counts" class="tool-save-taxid" data-testid="tool-save-taxid">
                <label :for="`${ids}-save-taxid`">Save for taxonomy ID</label>
                <input
                  :id="`${ids}-save-taxid`"
                  v-model="saveTaxidText"
                  type="text"
                  inputmode="numeric"
                  autocomplete="off"
                  spellcheck="false"
                  placeholder="9606"
                  size="10"
                  data-testid="tool-save-taxid-input"
                  @keydown.enter.prevent="saveForTaxid"
                />
                <button
                  ref="saveButton"
                  type="button"
                  :disabled="saving"
                  data-testid="tool-save-taxid-button"
                  @click="saveForTaxid"
                >
                  {{ saving ? 'Saving…' : 'Save' }}
                </button>
              </span>
            </template>
          </div>
          <p v-if="saveError" class="results-message error" role="alert" data-testid="tool-save-error">
            {{ saveError }}
          </p>
          <p v-if="saveStatus" class="muted small" data-testid="tool-save-status">
            {{ saveStatus }}
          </p>
          <ReplaceCountsDialog
            :open="pendingReplace !== undefined"
            :taxid="pendingReplace?.taxid ?? 0"
            :existing="pendingReplace?.existing"
            :new-file-name="pendingReplace?.newFileName ?? ''"
            :run-number="run.snapshot.number"
            @confirm="confirmReplace"
            @cancel="cancelReplace"
          />
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
      <template v-if="badge">
        <h3>Support status</h3>
        <VerificationBadge :badge="badge" />
      </template>

      <div class="details-grid">
        <h3>Command (fixed when the run was submitted)</h3>
        <dl class="details-list" data-testid="run-snapshot">
          <dt>Run</dt>
          <dd>
            {{ run.snapshot.number }}<template v-if="run.snapshot.group">, group run {{ run.snapshot.group.position }} of {{ run.snapshot.group.size }}</template>
          </dd>
          <dt>Tool</dt>
          <dd>{{ toolHeading }}</dd>
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
          <template v-if="run.record.runtimePath">
            <dt>Engine</dt>
            <dd data-detail="path">{{ run.record.runtimePath }}</dd>
          </template>
          <template v-if="run.record.engineBuild">
            <dt>Build</dt>
            <dd data-detail="build">{{ run.record.engineBuild }}</dd>
          </template>
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
