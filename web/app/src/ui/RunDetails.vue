<script setup lang="ts">
// The details of a run (plan §5.2, §5.8): its verification badge, what was fixed when it
// was queued (RunSnapshot), what happened when it ran (RunRecord), the CLI command of each
// output (with that of its own -outfmt, W10), and the warnings that the CLI writes to standard error.
import { computed } from 'vue';
import { CUSTOM_OUTPUT } from 'losat/domain';
import type { RunView } from '../application/coordinator';
import { optionGrammar, type LoadedRun } from '../application/results';
import { commandWords } from '../domain/argv';
import { COUNTS_FILE } from '../domain/counts';
import { dateOf } from '../domain/download-sentences';
import { customOutfmt } from '../domain/output-format';
import { programById } from '../domain/programs';
import CommandText from './CommandText.vue';
import { formatBytes, formatCount, formatDateTime, formatDuration, formatMillis } from './format';
import { vScrollFocus } from './scroll-focus';
import VerificationBadge from './VerificationBadge.vue';

const props = defineProps<{ run: RunView; loaded: LoadedRun }>();
const snapshot = computed(() => props.run.snapshot);
const record = computed(() => props.run.record);
const takesValue = computed(() => optionGrammar(props.loaded.description).takesValue);
/** The -outfmt value of the output of the run's own -outfmt, if it has one. */
const custom = computed(() => customOutfmt(snapshot.value.argv));
/** A time in ISO 8601 form, local time (S13 screen review L5: "09/10/2026" read either way). */
const time = (ms: number | undefined) => (ms === undefined ? '' : formatDateTime(ms));
const duration = computed(() =>
  record.value.startedAt !== undefined && record.value.endedAt !== undefined ? formatDuration(record.value.endedAt - record.value.startedAt) : '',
);
/** The counts of a taxonomy ID: the file's name, the site it came from and the date of the file there. */
const countsFile = COUNTS_FILE;
const countsHost = computed(() => (record.value.windowMasker?.origin === undefined ? '' : new URL(record.value.windowMasker.origin.url).host));
const countsDate = computed(() => dateOf(record.value.windowMasker?.origin?.lastModified));
const inputs = computed(() => [
  { role: 'Query', input: snapshot.value.query },
  { role: 'Subject', input: snapshot.value.subject },
]);
</script>

<template>
  <div class="run-details" data-testid="run-details">
    <h3>Support status</h3>
    <VerificationBadge :badge="loaded.badge" />

    <div class="details-grid">
      <h3>Search (fixed when the run was submitted)</h3>
      <dl class="details-list" data-testid="run-snapshot">
        <dt>Run</dt>
        <dd>
          {{ snapshot.number }}<template v-if="snapshot.group">, group run {{ snapshot.group.position }} of {{ snapshot.group.size }}</template>
        </dd>
        <dt>Mode</dt>
        <dd>{{ programById(snapshot.program).label }}</dd>
        <dt>Arguments</dt>
        <dd><code data-testid="run-argv"><CommandText :words="snapshot.argv" :takes-value="takesValue" /></code></dd>
        <template v-for="{ role, input } in inputs" :key="role">
          <dt>{{ role }}</dt>
          <dd :data-testid="`run-input-${role.toLowerCase()}`">
            {{ input.name }}: {{ formatCount(input.records.length) }} {{ input.records.length === 1 ? 'record' : 'records' }},
            {{ formatBytes(input.bytes.length) }}<br />
            <span class="muted small">SHA-256 {{ input.sha256 }}</span>
          </dd>
        </template>
        <template v-if="snapshot.conffile">
          <!-- The registry file of -conffile (line RG E1w-3b): the name that the argv gives, and the file read under it. -->
          <dt>-conffile</dt>
          <dd data-testid="run-conffile">{{ snapshot.conffile.name }}: {{ snapshot.conffile.fileName }}, {{ formatBytes(snapshot.conffile.bytes.length) }}</dd>
        </template>
        <template v-if="snapshot.dryRun">
          <dt>-dryrun</dt>
          <dd data-testid="run-dry-run">a dry run: no search</dd>
        </template>
        <dt>Threads requested</dt>
        <dd>{{ snapshot.requestedThreads === 'auto' ? 'Auto' : snapshot.requestedThreads }}</dd>
        <dt>Submitted</dt>
        <dd data-detail="queued">{{ time(snapshot.queuedAt) }}</dd>
      </dl>

      <h3>Run</h3>
      <dl class="details-list" data-testid="run-record">
        <dt>Engine</dt>
        <dd data-detail="path">
          {{ record.runtimePath }}, {{ record.threads }} {{ record.threads === 1 ? 'thread' : 'threads' }}
        </dd>
        <template v-if="record.fallbackReason">
          <dt>Serial because</dt>
          <dd>{{ record.fallbackReason }}</dd>
        </template>
        <dt>Build</dt>
        <dd data-detail="build">{{ record.engineBuild }}</dd>
        <template v-if="record.subjectRetained !== undefined">
          <dt>Subject</dt>
          <dd>{{ record.subjectRetained ? 'kept from the previous search' : 'read for this search' }}</dd>
        </template>
        <template v-if="record.windowMasker">
          <dt>WindowMasker counts</dt>
          <dd data-testid="run-window-masker">
            <!-- The counts of a taxonomy ID: the ID, the file and where it came from (the Last-Modified when it was downloaded). -->
            <template v-if="record.windowMasker.origin">
              {{ record.windowMasker.origin.taxid }}: {{ countsFile }} from {{ countsHost }}<template v-if="countsDate"> ({{ countsDate }})</template>,
            </template>
            <template v-else>{{ record.windowMasker.name }},</template>
            {{ formatBytes(record.windowMasker.size) }},
            {{ record.windowMasker.retained ? 'kept from the previous search' : `read for this search in ${formatMillis(record.windowMasker.loadMs)}` }}<br />
            <span class="muted small">SHA-256 {{ record.windowMasker.sha256 }}</span>
            <template v-if="record.windowMasker.origin"><br /><span class="muted small" data-testid="run-window-masker-url">{{ record.windowMasker.origin.url }}</span></template>
          </dd>
        </template>
        <template v-if="record.runtimeGeneration !== undefined">
          <dt>Runtime</dt>
          <dd>generation {{ record.runtimeGeneration }}</dd>
        </template>
        <template v-if="record.memory">
          <dt>Engine memory</dt>
          <dd>{{ formatBytes(record.memory.linearBytesAfter) }} after the search</dd>
        </template>
        <dt>Started</dt>
        <dd data-detail="started">{{ time(record.startedAt) }}</dd>
        <dt>Took</dt>
        <dd>{{ duration }}</dd>
      </dl>
    </div>

    <h3>Command line</h3>
    <ul class="commands">
      <li v-for="format in loaded.description.formats" :key="format">
        outfmt {{ format }}: <code :data-testid="`run-command-${format}`"><CommandText :words="commandWords(snapshot.argv, format)" :takes-value="takesValue" /></code>
      </li>
      <li v-if="custom !== undefined">
        outfmt {{ custom }}: <code data-testid="run-command-custom"><CommandText :words="commandWords(snapshot.argv, CUSTOM_OUTPUT)" :takes-value="takesValue" /></code>
      </li>
    </ul>

    <h3>Warnings</h3>
    <pre v-if="loaded.diagnostics !== ''" v-scroll-focus="'Warnings'" class="output" data-testid="run-diagnostics">{{ loaded.diagnostics }}</pre>
    <p v-else class="muted" data-testid="run-diagnostics-none">The engine wrote no warnings.</p>
  </div>
</template>
