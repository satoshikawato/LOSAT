<script setup lang="ts">
// The records of a source, with a box to include or leave out each one. Only the rows in
// view are drawn, so a source of many records stays responsive.
import { computed, ref } from 'vue';
import type { DraftSource, SearchDraft } from '../application/draft';
import type { InputRole, SequenceKind } from '../domain/programs';
import { looksLikeOtherKind } from '../domain/sequence-kind';
import { formatCount } from './format';

const props = defineProps<{
  draft: SearchDraft;
  role: InputRole;
  source: DraftSource;
  kind: SequenceKind;
  duplicates: ReadonlySet<string>;
  testid: string;
}>();

const ROW_PX = 30;
const VISIBLE_ROWS = 10;
const SPARE_ROWS = 6;

const filter = ref('');
const scrollTop = ref(0);
const records = computed(() => props.source.base?.records ?? []);
const excluded = computed(() => new Set(props.source.excluded));
const refused = computed(() => (props.source.check?.state === 'refused' ? props.source.check.record : undefined));
const shown = computed(() => {
  const text = filter.value.trim().toLowerCase();
  return text === '' ? records.value : records.value.filter((record) => record.id.toLowerCase().includes(text));
});
const viewportPx = computed(() => Math.max(1, Math.min(shown.value.length, VISIBLE_ROWS)) * ROW_PX);
const first = computed(() => Math.max(0, Math.floor(scrollTop.value / ROW_PX) - SPARE_ROWS));
const rows = computed(() => shown.value.slice(first.value, first.value + VISIBLE_ROWS + 2 * SPARE_ROWS));

function setShown(included: boolean): void {
  props.draft.setIncluded(
    props.role,
    props.source.key,
    shown.value.map((record) => record.index),
    included,
  );
}

function toggle(index: number, event: Event): void {
  props.draft.setIncluded(props.role, props.source.key, [index], (event.target as HTMLInputElement).checked);
}
</script>

<template>
  <div class="record-list" :data-testid="`${testid}-records`">
    <div class="record-tools">
      <label>
        <span class="visually-hidden">Find records by ID</span>
        <input v-model="filter" type="search" placeholder="Find by ID" :data-testid="`${testid}-filter`" />
      </label>
      <button type="button" :data-testid="`${testid}-include-shown`" @click="setShown(true)">Include shown</button>
      <button type="button" :data-testid="`${testid}-exclude-shown`" @click="setShown(false)">Exclude shown</button>
      <span class="muted">{{ formatCount(shown.length) }} shown</span>
    </div>
    <div
      class="record-viewport"
      :style="{ height: `${viewportPx}px` }"
      @scroll="scrollTop = ($event.target as HTMLElement).scrollTop"
    >
      <div class="record-spacer" :style="{ height: `${shown.length * ROW_PX}px` }">
        <div
          v-for="(record, i) in rows"
          :key="record.index"
          class="record-row"
          :class="{ excluded: excluded.has(record.index), refused: refused === record.index }"
          :style="{ top: `${(first + i) * ROW_PX}px`, height: `${ROW_PX}px` }"
        >
          <input
            type="checkbox"
            :checked="!excluded.has(record.index)"
            :aria-label="`Include record ${record.index + 1}, ${record.id}`"
            :data-testid="`${testid}-record-${record.index}`"
            @change="toggle(record.index, $event)"
          />
          <span class="record-number">#{{ record.index + 1 }}</span>
          <span class="record-id" :title="record.id">{{ record.id }}</span>
          <span class="record-length">{{ formatCount(record.length) }}</span>
          <span v-if="duplicates.has(record.id)" class="badge">same ID</span>
          <span v-if="looksLikeOtherKind(record.residue_counts, kind)" class="badge warn">
            {{ kind === 'nucleotide' ? 'protein?' : 'nucleotide?' }}
          </span>
          <span v-if="refused === record.index" class="badge error">refused</span>
        </div>
      </div>
    </div>
  </div>
</template>
