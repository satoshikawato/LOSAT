<script setup lang="ts">
// One input of a role: its name, its records, and what the engine's check says of it. The check's
// verdict and the warnings are told to a screen reader when they change (WA2-AX-03), in the
// card's words after its name; "Remove" and "Put each accession on its own line" ask the role's
// block, which moves the focus on, and the buttons that rewrite the box keep it on the card (rewrite).
import { computed, nextTick, ref, watch } from 'vue';
import type { DraftSource, DraftState, SearchDraft } from '../application/draft';
import { includedOf } from '../application/draft';
import { accessionsOnALine, isSeqIdReason, pointsToFetch } from '../domain/accession';
import { duplicateIds } from '../domain/dataset';
import { programById, residueUnit, sequenceKind, type InputRole } from '../domain/programs';
import { looksLikeOtherKind } from '../domain/sequence-kind';
import { announce } from './announcer';
import { formatBytes, formatCount } from './format';
import LosatLoader from './LosatLoader.vue';
import RecordList from './RecordList.vue';
import { vScrollFocus } from './scroll-focus';

const props = defineProps<{
  draft: SearchDraft;
  state: DraftState;
  role: InputRole;
  source: DraftSource;
  index: number;
}>();
defineEmits<{ remove: []; splitAccessions: [] }>();

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
 * A pasted text of accessions only that are not one per line (WA2-OP-04): Fetch takes them one per line, so the card offers
 * to put each on its own line. It does not depend on the reader: a list with no versions (AB000001 MH011443) reads as a
 * sequence (ABMH) and a ✓, though it is a list of accessions.
 */
const splitsAccessions = computed(
  () => props.source.origin === 'paste' && props.source.status !== 'indexing' && accessionsOnALine(props.state[props.role].paste) !== undefined,
);
const hasMark = computed(() => props.source.lead !== undefined && 'mark' in props.source.lead);
/**
 * The first line of a pasted text is an accession that the engine reads as a Seq-id: Fetch takes accessions only (WA2-OP-04).
 * Not when the box is a list of accessions, whose button puts them one per line, and so says the same.
 */
const fetchHint = computed(
  () =>
    props.source.origin === 'paste' &&
    props.source.status === 'failed' &&
    !splitsAccessions.value &&
    pointsToFetch(props.source.error, props.state[props.role].paste),
);
/**
 * A pasted text whose first line the engine reads as a Seq-id (which NCBI BLAST+ would fetch) and cannot read (WA2-OP-04): the
 * line is not data, and "Add the line >pasted_..." would make it the first data line, so the card puts ">" before it instead.
 */
const seqIdFirstLine = computed(
  () =>
    props.source.origin === 'paste' &&
    props.source.status === 'failed' &&
    !hasMark.value &&
    !splitsAccessions.value &&
    isSeqIdReason(props.source.error),
);
/** A pasted text that does not start with a defline and cannot be read; a byte order mark, a Seq-id and accessions on a line have their own hints. */
const needsDefline = computed(() => {
  const text = props.state[props.role].paste.trimStart();
  return (
    props.source.origin === 'paste' &&
    props.source.status === 'failed' &&
    !hasMark.value &&
    !seqIdFirstLine.value &&
    !splitsAccessions.value &&
    text !== '' &&
    !text.startsWith('>')
  );
});
/** The card's sentence on Fetch, for a first line that is an accession and a Seq-id to the engine. */
const FETCH_HINT = 'To fetch accessions, put only the accessions in the box, one per line.';
const refusedRecord = computed(() => (props.source.check?.state === 'refused' ? props.source.check.record : undefined));
/** Why a source that starts with a byte order mark cannot be read (WA-IN-05); the engine's message is NCBI's. */
const markHint = computed(() => {
  const lead = props.source.lead;
  if (props.source.status !== 'failed' || lead === undefined || !('mark' in lead)) return undefined;
  if (props.source.origin === 'paste') return 'The text starts with a byte order mark (U+FEFF).';
  return lead.mark === 'UTF-16'
    ? 'The file is UTF-16 text. Save it as plain UTF-8 or ASCII text.'
    : 'The file starts with a UTF-8 byte order mark. Save it as plain text without the mark (UTF-8 or ASCII).';
});
const FORMAT_NAMES = { GenBank: 'a GenBank flat file', EMBL: 'an EMBL or UniProt flat file', FASTQ: 'a FASTQ file' } as const;
/** Another sequence file format, which the engine reads as FASTA text as the command line does (WA-IN-04). */
const formatName = computed(() => {
  const lead = props.source.lead;
  return lead !== undefined && 'format' in lead ? FORMAT_NAMES[lead.format] : undefined;
});
const card = ref<HTMLElement>();
/**
 * Asks the draft to rewrite the box ("Add the line", "Put '>'", "Remove the byte order mark"): the input is read again, so
 * the button goes, and the focus goes to the card, which the live region then speaks of, not to the page's body (WA2-AX-04).
 */
async function rewrite(change: () => void): Promise<void> {
  change();
  await nextTick();
  card.value?.focus();
}
const testid = computed(() => `${props.role}-source-${props.index}`);

const recordsText = (count: number) => `${formatCount(count)} ${count === 1 ? 'record' : 'records'}`;
/**
 * What the card says once the input is read and checked (undefined until then), as the screen shows
 * it, and its facts without the mode's name (`key`).
 */
const verdict = computed(() => {
  const source = props.source;
  const parts: string[] = [];
  const check = source.check;
  if (source.status === 'indexing' || (source.status !== 'failed' && (check === undefined || check.state === 'pending'))) return undefined;
  const key = JSON.stringify([
    source.status,
    source.error,
    markHint.value,
    fetchHint.value,
    check?.state,
    check === undefined ? undefined : check.state === 'ok' ? check.records : 'message' in check ? check.message : undefined,
    otherKind.value,
    otherKind.value > 0 ? kind.value : '',
    duplicateCount.value,
    formatName.value,
  ]);
  if (source.status === 'failed') {
    parts.push(`This input cannot be read: ${source.error ?? ''}`);
    if (fetchHint.value) parts.push(FETCH_HINT);
    if (markHint.value) parts.push(markHint.value);
  } else if (check !== undefined) {
    if (check.state === 'ok') parts.push(`Accepted by the ${programLabel.value} engine: ${recordsText(check.records)}`);
    else if (check.state === 'refused') parts.push(`The ${programLabel.value} engine refuses this input: ${check.message}`);
    else if (check.state === 'error') parts.push(`The input could not be checked here (${check.message}).`);
    if (otherKind.value > 0) {
      parts.push(
        `${formatCount(otherKind.value)} ${otherKind.value === 1 ? 'record looks' : 'records look'} like ${otherKindName.value} sequences, but ${programLabel.value} reads the ${props.role} as ${kind.value}.`,
      );
    }
    if (duplicateCount.value > 0) parts.push(`${formatCount(duplicateCount.value)} records share an ID with another record.`);
  }
  if (formatName.value) parts.push(`This looks like ${formatName.value}. LOSAT reads it as FASTA, so its text becomes sequence. Convert it to FASTA first.`);
  return { text: `${title.value}: ${parts.join(' ')}`, key };
});
/**
 * Told when it changes: not again after each check of the same input (each pause in typing, each
 * change of mode) that finds what the last one found.
 */
/** A card that shows an input already checked (mounted again) does not tell it again. */
let toldKey = verdict.value?.key;
watch(verdict, (now) => {
  if (now === undefined || now.key === toldKey) return;
  toldKey = now.key;
  announce(now.text);
});
</script>

<template>
  <li
    ref="card"
    tabindex="-1"
    class="source"
    :data-testid="testid"
    :data-status="source.status"
    :data-index-ms="source.indexMs"
    :data-check-ms="source.checkMs"
  >
    <div class="source-head">
      <strong class="source-name">{{ title }}</strong>
      <span class="muted">{{ formatBytes(source.size) }}</span>
      <button type="button" class="link source-remove" :data-testid="`${testid}-remove`" @click="$emit('remove')">
        Remove
      </button>
    </div>

    <p v-if="source.status === 'indexing'" class="muted"><LosatLoader />Reading the records…</p>
    <template v-else-if="source.status === 'failed'">
      <p class="error" :data-testid="`${testid}-error`">This input cannot be read: {{ source.error }}</p>
      <p v-if="fetchHint" class="hint" :data-testid="`${testid}-fetch-hint`">{{ FETCH_HINT }}</p>
      <p v-if="splitsAccessions" class="hint">
        <button type="button" class="link" :data-testid="`${testid}-split-accessions`" @click="$emit('splitAccessions')">
          Put each accession on its own line
        </button>
      </p>
      <p v-if="seqIdFirstLine" class="hint">
        <button type="button" class="link" :data-testid="`${testid}-put-defline-mark`" @click="rewrite(() => draft.markFirstLineAsDefline(role))">
          Put '&gt;' before the first line
        </button>
      </p>
      <p v-if="needsDefline" class="hint">
        <button type="button" class="link" :data-testid="`${testid}-add-defline`" @click="rewrite(() => draft.addDefline(role))">
          Add the line "&gt;pasted_{{ role }}" before the sequence
        </button>
      </p>
      <p v-if="markHint" class="hint" :data-testid="`${testid}-mark`">
        {{ markHint }}
        <button
          v-if="source.origin === 'paste'"
          type="button"
          class="link"
          :data-testid="`${testid}-remove-mark`"
          @click="rewrite(() => draft.removeByteOrderMark(role))"
        >
          Remove the byte order mark
        </button>
      </p>
      <p v-if="formatName" class="warning" :data-testid="`${testid}-format-warning`">
        This looks like {{ formatName }}. LOSAT reads it as FASTA, so its text becomes sequence. Convert it to FASTA first.
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
        <LosatLoader />Checking…
      </p>
      <p v-else-if="source.check?.state === 'ok'" class="ok" :data-testid="`${testid}-check`" data-check="ok">
        <svg aria-hidden="true" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5"><path d="M5 13l4 4L19 7" /></svg>
        <span class="visually-hidden">Accepted by the {{ programLabel }} engine: {{ formatCount(source.check.records) }} {{ source.check.records === 1 ? 'record' : 'records' }}</span>
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
        The input could not be checked here ({{ source.check.message }}).
      </p>

      <p v-if="splitsAccessions" class="hint">
        <button type="button" class="link" :data-testid="`${testid}-split-accessions`" @click="$emit('splitAccessions')">
          Put each accession on its own line
        </button>
      </p>

      <p v-if="formatName" class="warning" :data-testid="`${testid}-format-warning`">
        This looks like {{ formatName }}. LOSAT reads it as FASTA, so its text becomes sequence. Convert it to FASTA first.
      </p>
      <p v-if="otherKind > 0" class="warning" :data-testid="`${testid}-kind-warning`">
        {{ formatCount(otherKind) }} {{ otherKind === 1 ? 'record looks' : 'records look' }} like {{ otherKindName }}
        sequences, but {{ programLabel }} reads the {{ role }} as {{ kind }}.
      </p>
      <p v-if="duplicateCount > 0" class="warning" :data-testid="`${testid}-duplicates`">
        {{ formatCount(duplicateCount) }} records share an ID with another record.
      </p>

      <details v-if="source.head" class="head">
        <summary>First lines</summary>
        <pre v-scroll-focus="'First lines'" :data-testid="`${testid}-head`">{{ source.head }}</pre>
      </details>
      <!-- Closed, as NCBI shows only the file (Owner, 2026-11-02): opened to leave records out. -->
      <details class="records">
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
