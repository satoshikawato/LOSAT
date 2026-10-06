<script setup lang="ts">
// The region of a role's only record (DW-9): two position fields, and a bar that shows the
// record, on which a range can be chosen by dragging. The region becomes -query_loc or
// -subject_loc start-stop.
import { computed, ref } from 'vue';
import type { SearchDraft } from '../application/draft';
import type { DatasetRecord } from '../domain/dataset';
import type { InputRole } from '../domain/programs';
import { REGION_FLAG, regionFromPositions, regionProblem, regionRange, regionValue, type RegionText } from '../domain/region';
import { formatCount } from './format';

const props = defineProps<{
  draft: SearchDraft;
  role: InputRole;
  record: DatasetRecord;
  region: RegionText | undefined;
  unit: string;
}>();

const bar = ref<HTMLElement>();
const dragFrom = ref<number>();
const testid = computed(() => `${props.role}-region`);
const problem = computed(() => (props.region === undefined ? undefined : regionProblem(props.region, props.record.length)));
const range = computed(() => (props.region === undefined ? undefined : regionRange(props.region, props.record.length)));
const selection = computed(() => {
  const value = range.value;
  if (value === undefined) return undefined;
  const length = props.record.length;
  const left = ((Math.min(value.start, value.stop) - 1) / length) * 100;
  const width = ((Math.abs(value.stop - value.start) + 1) / length) * 100;
  return { left: `${left}%`, width: `${Math.max(width, 0.4)}%` };
});

function update(change: Partial<RegionText>): void {
  const current = props.region ?? { start: '', stop: '' };
  props.draft.setRegion(props.role, { ...current, ...change });
}

/** The 1-based position under the pointer. */
function position(event: PointerEvent): number {
  const rect = bar.value!.getBoundingClientRect();
  const fraction = Math.min(1, Math.max(0, (event.clientX - rect.left) / Math.max(1, rect.width)));
  return 1 + fraction * (props.record.length - 1);
}

function down(event: PointerEvent): void {
  if (props.record.length < 1) return;
  bar.value!.setPointerCapture(event.pointerId);
  dragFrom.value = position(event);
  props.draft.setRegion(props.role, regionFromPositions(dragFrom.value, dragFrom.value, props.record.length));
}

function move(event: PointerEvent): void {
  if (dragFrom.value === undefined) return;
  props.draft.setRegion(props.role, regionFromPositions(dragFrom.value, position(event), props.record.length));
}

function up(): void {
  dragFrom.value = undefined;
}
</script>

<template>
  <fieldset class="region" :data-testid="testid">
    <legend>Region of {{ record.id }}</legend>
    <div
      ref="bar"
      class="region-bar"
      role="presentation"
      :data-testid="`${testid}-bar`"
      @pointerdown.prevent="down"
      @pointermove="move"
      @pointerup="up"
      @pointercancel="up"
    >
      <div v-if="selection" class="region-selection" :style="selection" />
    </div>
    <div class="region-scale muted">
      <span>1</span>
      <span>{{ formatCount(record.length) }} {{ unit }}</span>
    </div>
    <div class="region-fields">
      <label>
        From
        <input
          type="text"
          inputmode="numeric"
          :value="region?.start ?? ''"
          placeholder="1"
          :data-testid="`${testid}-start`"
          @input="update({ start: ($event.target as HTMLInputElement).value })"
        />
      </label>
      <label>
        To
        <input
          type="text"
          inputmode="numeric"
          :value="region?.stop ?? ''"
          :placeholder="String(record.length)"
          :data-testid="`${testid}-stop`"
          @input="update({ stop: ($event.target as HTMLInputElement).value })"
        />
      </label>
      <button v-if="region" type="button" :data-testid="`${testid}-clear`" @click="draft.setRegion(role, undefined)">
        Whole record
      </button>
    </div>
    <p v-if="problem" class="error" :data-testid="`${testid}-problem`">{{ problem }}</p>
    <p v-else-if="region" class="hint" :data-testid="`${testid}-argument`">
      Searches {{ REGION_FLAG[role] }} {{ regionValue(region) }}. Positions in the results stay those of the record.
    </p>
    <p v-else class="hint">The whole record is searched. Drag on the bar or enter positions (1-based, inclusive) to search a part.</p>
  </fieldset>
</template>
