<script setup lang="ts">
// The details of a run (plan §5.2, §5.8): its verification badge, what was fixed when it
// was queued (RunSnapshot), what happened when it ran (RunRecord), the CLI command of each
// output (with that of its own -outfmt, W10), and the warnings that the CLI writes to standard error, or where they are
// when the run's -logfile took them.
import { computed } from 'vue';
import { CUSTOM_OUTPUT } from 'losat/domain';
import type { RunView } from '../application/coordinator';
import { optionGrammar, type LoadedRun } from '../application/results';
import { commandWords, shellWord, userArgv } from '../domain/argv';
import { COUNTS_FILE } from '../domain/counts';
import { dateOf } from '../domain/download-sentences';
import { fromText, worthSaying } from '../domain/input-parts';
import { customOutfmt } from '../domain/output-format';
import { programById } from '../domain/programs';
import CommandText from './CommandText.vue';
import { formatBytes, formatCount, formatDateTime, formatDuration, formatMillis } from './format';
import { NO_WARNINGS, warningsPlaceParts } from './message-place-text';
import PlaceWords from './PlaceWords.vue';
import { vScrollFocus } from './scroll-focus';
import VerificationBadge from './VerificationBadge.vue';

const props = defineProps<{ run: RunView; loaded: LoadedRun }>();
const snapshot = computed(() => props.run.snapshot);
const record = computed(() => props.run.record);
/**
 * Whether the run searched: a run that printed the command line's help or version (`RunSnapshot.prints`, line WQA's D14)
 * compared no options, wrote no output and used no engine memory for a search, so Run details says nothing of those
 * (line WQA's D20, S3-1, 2026-11-23).
 */
const searched = computed(() => props.run.snapshot.prints === undefined);
const takesValue = computed(() => optionGrammar(props.loaded.description).takesValue);
const commandEnv = computed(() =>
  snapshot.value.database?.blastdbEnv ? { blastdb: snapshot.value.database.blastdbEnv } : undefined,
);
/** The -outfmt value of the output of the run's own -outfmt, if it has one. */
const custom = computed(() => customOutfmt(snapshot.value.argv));
/** A time in ISO 8601 form, local time (S13 screen review L5: "09/10/2026" read either way). */
const time = (ms: number | undefined) => (ms === undefined ? '' : formatDateTime(ms));
const duration = computed(() =>
  record.value.startedAt !== undefined && record.value.endedAt !== undefined ? formatDuration(record.value.endedAt - record.value.startedAt) : '',
);
/** The counts of a taxonomy ID: the file's name, its origin and date (downloaded or made in this browser). */
const countsFile = COUNTS_FILE;
const countsHost = computed(() => {
  const origin = record.value.windowMasker?.origin;
  return origin?.kind === 'downloaded' ? new URL(origin.url).host : '';
});
const countsDate = computed(() => {
  const origin = record.value.windowMasker?.origin;
  if (origin === undefined) return undefined;
  return origin.kind === 'made' ? new Date(origin.made.madeAt).toISOString().slice(0, 10) : dateOf(origin.lastModified);
});
/**
 * What Warnings says besides the diagnostics: where the run's warnings are when they went to the file of its -logfile or
 * nowhere, or that the engine wrote none (line WQA's sweep, D7, 2026-11-20).
 */
const warningsPlace = computed(() => warningsPlaceParts(props.loaded.messages, props.loaded.diagnostics, props.run));
/** Whether those words say that the engine wrote no warnings (its place holds none). */
const noWarnings = computed(() => warningsPlace.value?.[0]?.text === NO_WARNINGS);
const inputs = computed(() => [
  { role: 'Query', input: snapshot.value.query },
  ...(snapshot.value.database === undefined ? [{ role: 'Subject', input: snapshot.value.subject }] : []),
]);
</script>

<template>
  <div class="run-details" data-testid="run-details">
    <template v-if="loaded.badge">
      <h3>Support status</h3>
      <VerificationBadge :badge="loaded.badge" />
    </template>

    <div class="details-grid">
      <!-- A printed run is a command, not a search (line WQA's D20, 2026-11-23). -->
      <h3>{{ searched ? 'Search' : 'Command' }} (fixed when the run was submitted)</h3>
      <dl class="details-list" data-testid="run-snapshot">
        <dt>Run</dt>
        <dd>
          {{ snapshot.number }}<template v-if="snapshot.group">, group run {{ snapshot.group.position }} of {{ snapshot.group.size }}</template>
        </dd>
        <dt>Mode</dt>
        <dd>{{ programById(snapshot.program).label }}</dd>
        <dt>Arguments</dt>
        <dd><code data-testid="run-argv"><CommandText :words="userArgv(snapshot)" :takes-value="takesValue" /></code></dd>
        <template v-for="{ role, input } in searched ? inputs : []" :key="role">
          <dt>{{ role }}</dt>
          <dd :data-testid="`run-input-${role.toLowerCase()}`">
            {{ input.name }}: {{ formatCount(input.records.length) }} {{ input.records.length === 1 ? 'record' : 'records' }},
            {{ formatBytes(input.bytes.length) }}<br />
            <!-- The file of the command is not the user's as given: what it is made of, and which records it leaves out (line WQA's D22, 2026-11-23). -->
            <template v-if="worthSaying(input.from)"><span :data-testid="`run-input-${role.toLowerCase()}-from`">From {{ fromText(input.from) }}</span><br /></template>
            <span class="muted small">SHA-256 {{ input.sha256 }}</span>
          </dd>
        </template>
        <template v-if="searched && snapshot.database">
          <dt>Database</dt>
          <dd data-testid="run-database">
            <code>{{ snapshot.database.name }}</code> ({{ snapshot.database.kind }}, {{ formatCount(snapshot.database.fileCount) }}
            {{ snapshot.database.fileCount === 1 ? 'file' : 'files' }}, {{ formatBytes(snapshot.database.totalBytes) }}<template
              v-if="snapshot.database.taxdb || snapshot.database.sqlite"
            >,
              {{ [snapshot.database.taxdb ? 'taxdb' : '', snapshot.database.sqlite ? 'taxonomy4blast.sqlite3' : ''].filter(Boolean).join(' + ') }}</template
            >)
            <template v-if="snapshot.database.blastdbEnv">
              <br />
              <span class="muted small" data-testid="run-database-blastdb">
                CLI reproduction: run from the folder's parent with <code>{{ shellWord(`BLASTDB=${snapshot.database.blastdbEnv}`) }}</code>
              </span>
            </template>
            <template v-if="snapshot.database.filterFiles && snapshot.database.filterFiles.length > 0">
              <br />
              <span class="muted small" data-testid="run-database-filters">
                Filter lists:
                {{ snapshot.database.filterFiles.map((f) => `${f.flag} ${f.name} (${formatBytes(f.size)})`).join(', ') }}
              </span>
            </template>
          </dd>
        </template>
        <template v-if="searched && snapshot.filteringDb">
          <dt>Repeats database</dt>
          <dd data-testid="run-filtering-db">
            <code>{{ snapshot.filteringDb.name }}</code> ({{ formatCount(snapshot.filteringDb.fileCount) }}
            {{ snapshot.filteringDb.fileCount === 1 ? 'file' : 'files' }}, {{ formatBytes(snapshot.filteringDb.totalBytes) }})
          </dd>
        </template>
        <template v-if="snapshot.conffile">
          <!-- The registry file of -conffile (line RG E1w-3b): the name that the argv gives, and the file read under it. -->
          <dt>-conffile</dt>
          <dd>
            <span data-testid="run-conffile">{{ snapshot.conffile.name }}: {{ snapshot.conffile.fileName }}, {{ formatBytes(snapshot.conffile.size) }}</span><br />
            <span class="muted small" data-testid="run-conffile-sha256">SHA-256 {{ snapshot.conffile.sha256 }}</span>
          </dd>
        </template>
        <template v-if="snapshot.dryRun">
          <dt>-dryrun</dt>
          <dd data-testid="run-dry-run">a dry run: no search</dd>
        </template>
        <template v-if="snapshot.prints">
          <!-- The command line's help or version, which the run printed in place of a search (line WQA's D14, 2026-11-20). -->
          <dt>Printed</dt>
          <dd data-testid="run-printed">the command line's {{ snapshot.prints }}: no search</dd>
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
        <template v-if="searched && record.subjectRetained !== undefined">
          <dt>Subject</dt>
          <dd>{{ record.subjectRetained ? 'kept from the previous search' : 'read for this search' }}</dd>
        </template>
        <template v-if="record.windowMasker">
          <dt>WindowMasker counts</dt>
          <dd data-testid="run-window-masker">
            <!-- The counts of a taxonomy ID: the ID, the file and where it came from (downloaded or made in this browser). -->
            <template v-if="record.windowMasker.origin?.kind === 'downloaded'">
              {{ record.windowMasker.origin.taxid }}: {{ countsFile }} from {{ countsHost }}<template v-if="countsDate"> ({{ countsDate }})</template>,
            </template>
            <template v-else-if="record.windowMasker.origin?.kind === 'made'">
              {{ record.windowMasker.origin.taxid }}: {{ countsFile }} made in this browser by run {{ record.windowMasker.origin.made.run }} ({{ countsDate }}),
            </template>
            <template v-else>{{ record.windowMasker.name }},</template>
            {{ formatBytes(record.windowMasker.size) }},
            {{ record.windowMasker.retained ? 'kept from the previous search' : `read for this search in ${formatMillis(record.windowMasker.loadMs)}` }}<br />
            <span class="muted small">SHA-256 {{ record.windowMasker.sha256 }}</span>
            <template v-if="record.windowMasker.origin?.kind === 'downloaded'"><br /><span class="muted small" data-testid="run-window-masker-url">{{ record.windowMasker.origin.url }}</span></template>
          </dd>
        </template>
        <template v-if="record.runtimeGeneration !== undefined">
          <dt>Runtime</dt>
          <dd>generation {{ record.runtimeGeneration }}</dd>
        </template>
        <template v-if="searched && record.memory">
          <dt>Engine memory</dt>
          <dd>{{ formatBytes(record.memory.linearBytesAfter) }} after the search</dd>
        </template>
        <dt>Started</dt>
        <dd data-detail="started">{{ time(record.startedAt) }}</dd>
        <dt>Took</dt>
        <dd>{{ duration }}</dd>
      </dl>
    </div>

    <template v-if="searched">
      <h3>Command line</h3>
      <ul class="commands">
        <li v-for="format in loaded.description.formats" :key="format">
          outfmt {{ format }}: <code :data-testid="`run-command-${format}`"><CommandText :words="commandWords(snapshot.argv, format, commandEnv)" :takes-value="takesValue" /></code>
        </li>
        <li v-if="custom !== undefined">
          outfmt {{ custom }}: <code data-testid="run-command-custom"><CommandText :words="commandWords(snapshot.argv, CUSTOM_OUTPUT, commandEnv)" :takes-value="takesValue" /></code>
        </li>
      </ul>
    </template>

    <h3>Warnings</h3>
    <pre v-if="loaded.diagnostics !== ''" v-scroll-focus="'Warnings'" class="output" data-testid="run-diagnostics">{{ loaded.diagnostics }}</pre>
    <p v-if="warningsPlace !== undefined" class="muted" :data-testid="noWarnings ? 'run-diagnostics-none' : 'run-diagnostics-place'">
      <PlaceWords :parts="warningsPlace" />
    </p>
  </div>
</template>
