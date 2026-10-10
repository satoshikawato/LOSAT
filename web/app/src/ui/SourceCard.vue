<script setup lang="ts">
import { computed } from 'vue';
import type { DraftSource, DraftState, SearchDraft } from '../application/draft';
import { includedOf } from '../application/draft';
import { duplicateIds, isFirstLineRefusal } from '../domain/dataset';
import { programById, residueUnit, sequenceKind, type InputRole } from '../domain/programs';
import { looksLikeOtherKind } from '../domain/sequence-kind';
import { formatBytes, formatCount } from './format';
import RecordList from './RecordList.vue';

const props = defineProps<{
  draft: SearchDraft;
  state: DraftState;
  role: InputRole;
  source: DraftSource;
  index: number;
}>();

const kind = computed(() => sequenceKind(programById(props.state.program), props.role));
const unit = computed(() => residueUnit(kind.value));
const records = computed(() => props.source.base?.records ?? []);
const included = computed(() => includedOf(props.source));
const totalLength = computed(() => included.value.reduce((sum, record) => sum + record.length, 0));
const duplicates = computed(() => duplicateIds(records.value));
const duplicateCount = computed(() => records.value.filter((record) => duplicates.value.has(record.id)).length);
const otherKind = computed(
  () => included.value.filter((record) => looksLikeOtherKind(record.residue_counts, kind.value)).length,
);
const otherKindName = computed(() => (kind.value === 'nucleotide' ? 'protein' : 'nucleotide'));
const programLabel = computed(() => programById(props.state.program).label);
const title = computed(() => (props.source.origin === 'paste' ? `Pasted sequences (${props.source.name})` : props.source.name));
/**
 * A pasted text that the index scan refuses for its first line. The engine's reader reads
 * residues before the first '>' as a record without a defline, but LOSAT Web refuses a first
 * line that NCBI BLAST+ may fetch as a sequence identifier, and a defline fixes that refusal
 * only (not a gap line or a "Near line N" one, which would be refused again further on).
 */
const needsDefline = computed(
  () =>
    props.source.origin === 'paste' &&
    props.source.status === 'failed' &&
    props.source.error !== undefined &&
    isFirstLineRefusal(props.source.error),
);
const refusedRecord = computed(() => (props.source.check?.state === 'refused' ? props.source.check.record : undefined));
const testid = computed(() => `${props.role}-source-${props.index}`);
</script>

<template>
  <li
    class="source"
    :data-testid="testid"
    :data-status="source.status"
    :data-index-ms="source.indexMs"
    :data-check-ms="source.checkMs"
  >
    <div class="source-head">
      <!-- One box for the name and the size: on phones the size follows the name and wraps under it, so that "Remove" stays on the first line (S13b screen review 3 L1). -->
      <span class="source-title"><strong class="source-name">{{ title }}</strong> <span class="source-size muted">{{ formatBytes(source.size) }}</span></span>
      <button type="button" class="link" :data-testid="`${testid}-remove`" @click="draft.removeSource(role, source.key)">
        Remove
      </button>
    </div>

    <p v-if="source.notice !== undefined" class="notice" :data-testid="`${testid}-notice`">{{ source.notice }}</p>
    <p v-if="source.status === 'indexing'" class="muted">Reading the records…</p>
    <template v-else-if="source.status === 'failed'">
      <p class="error" :data-testid="`${testid}-error`">This input cannot be read: {{ source.error }}</p>
      <p v-if="needsDefline" class="hint">
        FASTA starts each record with a line that begins with "&gt;".
        <button type="button" class="link" :data-testid="`${testid}-add-defline`" @click="draft.addDefline(role)">
          Add the line "&gt;pasted_{{ role }}" before the sequence
        </button>
      </p>
    </template>
    <template v-else>
      <p class="summary" :data-testid="`${testid}-summary`">
        {{ formatCount(records.length) }} {{ records.length === 1 ? 'record' : 'records' }}
        <template v-if="included.length !== records.length">({{ formatCount(included.length) }} included)</template>
        · {{ formatCount(totalLength) }} {{ unit }}
      </p>

      <p
        v-if="source.check?.state === 'pending'"
        class="muted"
        :data-testid="`${testid}-check`"
        data-check="pending"
      >
        Checking with the {{ programLabel }} engine…
      </p>
      <p v-else-if="source.check?.state === 'ok'" class="ok" :data-testid="`${testid}-check`" data-check="ok">
        The {{ programLabel }} engine reads {{ formatCount(source.check.records) }}
        {{ source.check.records === 1 ? 'record' : 'records' }}.
      </p>
      <div v-else-if="source.check?.state === 'refused'" class="refused" :data-testid="`${testid}-check`" data-check="refused">
        <p class="error">The {{ programLabel }} engine refuses this input: {{ source.check.message }}</p>
        <button
          v-if="refusedRecord !== undefined"
          type="button"
          :data-testid="`${testid}-exclude-refused`"
          @click="draft.setIncluded(role, source.key, [refusedRecord], false)"
        >
          Exclude record #{{ refusedRecord + 1 }} ({{ records[refusedRecord]?.id }})
        </button>
      </div>
      <p v-else-if="source.check?.state === 'error'" class="error" :data-testid="`${testid}-check`" data-check="error">
        The input could not be checked here ({{ source.check.message }}). The search reads it with the engine again.
      </p>

      <p v-if="otherKind > 0" class="warning" :data-testid="`${testid}-kind-warning`">
        {{ formatCount(otherKind) }} {{ otherKind === 1 ? 'record looks' : 'records look' }} like {{ otherKindName }}
        sequences, but {{ programLabel }} reads the {{ role }} as {{ kind }}. This is an estimate of LOSAT Web; the
        program is not changed.
      </p>
      <p v-if="duplicateCount > 0" class="warning" :data-testid="`${testid}-duplicates`">
        {{ formatCount(duplicateCount) }} records share an ID with another record. They are told apart by their number
        (#).
      </p>

      <details v-if="source.head" class="head">
        <summary>First lines</summary>
        <pre :data-testid="`${testid}-head`">{{ source.head }}</pre>
      </details>
      <details class="records" :open="records.length > 1 && records.length <= 20">
        <summary>Records</summary>
        <RecordList
          :draft="draft"
          :role="role"
          :source="source"
          :kind="kind"
          :duplicates="duplicates"
          :testid="testid"
        />
      </details>
    </template>
  </li>
</template>
