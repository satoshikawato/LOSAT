<script setup lang="ts">
// The details of a run (plan §5.2, §5.8): its verification badge, what was fixed when it
// was queued (RunSnapshot), what happened when it ran (RunRecord), the CLI command of each
// output format, and the warnings that the CLI writes to standard error.
import { computed } from 'vue';
import type { RunView } from '../application/coordinator';
import type { LoadedRun } from '../application/results';
import { toShellCommand } from '../domain/argv';
import { programById } from '../domain/programs';
import CommandText from './CommandText.vue';
import { formatBytes, formatCount, formatDateTime, formatDuration } from './format';
import VerificationBadge from './VerificationBadge.vue';

const props = defineProps<{ run: RunView; loaded: LoadedRun }>();
const snapshot = computed(() => props.run.snapshot);
const record = computed(() => props.run.record);
/** A time in ISO 8601 form, local time (S13 screen review L5: "09/10/2026" read either way). */
const time = (ms: number | undefined) => (ms === undefined ? '' : formatDateTime(ms));
const duration = computed(() =>
  record.value.startedAt !== undefined && record.value.endedAt !== undefined ? formatDuration(record.value.endedAt - record.value.startedAt) : '',
);
const inputs = computed(() => [
  { role: 'Query', input: snapshot.value.query },
  { role: 'Subject', input: snapshot.value.subject },
]);
</script>

<template>
  <div class="run-details" data-testid="run-details">
    <h3>Verification</h3>
    <VerificationBadge :badge="loaded.badge" />

    <div class="details-grid">
      <h3>Search (fixed when the run was queued)</h3>
      <dl class="details-list" data-testid="run-snapshot">
        <dt>Run</dt>
        <dd>
          {{ snapshot.number }}<template v-if="snapshot.group">, group run {{ snapshot.group.position }} of {{ snapshot.group.size }}</template>
        </dd>
        <dt>Program</dt>
        <dd>{{ programById(snapshot.program).label }}</dd>
        <dt>Arguments</dt>
        <dd><code data-testid="run-argv"><CommandText :text="snapshot.argv.join(' ')" /></code></dd>
        <template v-for="{ role, input } in inputs" :key="role">
          <dt>{{ role }}</dt>
          <dd :data-testid="`run-input-${role.toLowerCase()}`">
            {{ input.name }}: {{ formatCount(input.records.length) }} {{ input.records.length === 1 ? 'record' : 'records' }},
            {{ formatBytes(input.bytes.length) }}<br />
            <span class="muted small">SHA-256 {{ input.sha256 }}</span>
          </dd>
        </template>
        <dt>Threads requested</dt>
        <dd>{{ snapshot.requestedThreads === 'auto' ? 'Auto' : snapshot.requestedThreads }}</dd>
        <dt>Queued</dt>
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
    <p class="muted small">The LOSAT command that writes each output from the same inputs.</p>
    <ul class="commands">
      <li v-for="format in loaded.description.formats" :key="format">
        outfmt {{ format }}: <code :data-testid="`run-command-${format}`"><CommandText :text="toShellCommand(snapshot.argv, format)" /></code>
      </li>
    </ul>
    <p v-for="(line, i) in loaded.badge.exceptions" :key="i" class="notice">{{ line }}</p>

    <h3>Warnings</h3>
    <pre v-if="loaded.diagnostics !== ''" class="output" data-testid="run-diagnostics">{{ loaded.diagnostics }}</pre>
    <p v-else class="muted" data-testid="run-diagnostics-none">The engine wrote no warnings.</p>
  </div>
</template>
