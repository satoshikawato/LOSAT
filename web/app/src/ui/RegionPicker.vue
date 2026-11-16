<script setup lang="ts">
// The region of a role's only record (DW-9; NCBI's "Query subrange" / "Subject subrange"):
// two position fields, and a bar that shows the record, on which a range can be chosen by
// dragging. The region becomes -query_loc or -subject_loc start-stop.
import { computed, nextTick, ref } from 'vue';
import type { InputRole } from 'losat/domain';
import type { SearchDraft } from '../application/draft';
import type { DatasetRecord } from '../domain/dataset';
import { regionFromPositions, regionProblem, regionProblemFields, regionRange, swappedRegion, type RegionText } from '../domain/region';
import { formatCount } from './format';
import { refusedMark } from './refusal';

const props = defineProps<{
  draft: SearchDraft;
  role: InputRole;
  /** The fieldset's legend ("Query subrange"). */
  title: string;
  record: DatasetRecord;
  region: RegionText | undefined;
  unit: string;
}>();

const bar = ref<HTMLElement>();
const startInput = ref<HTMLInputElement>();
const dragFrom = ref<number>();
const testid = computed(() => `${props.role}-region`);
const problem = computed(() => (props.region === undefined ? undefined : regionProblem(props.region, props.record.length)));
/** A start after the stop: the swap, offered beside the problem (Owner, 2026-11-04). */
const swapped = computed(() => (props.region === undefined ? undefined : swappedRegion(props.region, props.record.length)));
const range = computed(() => (props.region === undefined ? undefined : regionRange(props.region, props.record.length)));
/** The fields that the problem is about, marked as refused and pointing to it (WA2-AX-05). */
const refused = computed(() => (props.region === undefined ? { start: false, stop: false } : regionProblemFields(props.region, props.record.length)));
const mark = (on: boolean) => refusedMark(on, `${testid.value}-problem`);
const selection = computed(() => {
  const value = range.value;
  if (value === undefined) return undefined;
  const length = props.record.length;
  const left = ((Math.min(value.start, value.stop) - 1) / length) * 100;
  const width = ((Math.abs(value.stop - value.start) + 1) / length) * 100;
  return { left: `${left}%`, width: `${Math.max(width, 0.4)}%` };
});

/** Swaps the fields; the button goes with the problem, so the focus moves to From. */
function swap(region: RegionText): void {
  props.draft.setRegion(props.role, region);
  void nextTick(() => startInput.value?.focus());
}

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
}

/** A drag sets the region; a click alone does not (a range of one letter is not a region). */
function move(event: PointerEvent): void {
  if (dragFrom.value === undefined) return;
  const region = regionFromPositions(dragFrom.value, position(event), props.record.length);
  if (region.start !== region.stop) props.draft.setRegion(props.role, region);
}

function up(): void {
  dragFrom.value = undefined;
}
</script>

<template>
  <fieldset class="region" :data-testid="testid">
    <legend>{{ title }}</legend>
    <p class="region-record" :data-testid="`${testid}-record`">Record {{ record.id }}</p>
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
          ref="startInput"
          :value="region?.start ?? ''"
          placeholder="1"
          v-bind="mark(refused.start)"
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
          v-bind="mark(refused.stop)"
          :data-testid="`${testid}-stop`"
          @input="update({ stop: ($event.target as HTMLInputElement).value })"
        />
      </label>
      <button v-if="region" type="button" :data-testid="`${testid}-clear`" @click="draft.setRegion(role, undefined)">
        Whole record
      </button>
    </div>
    <p v-if="problem" class="error">
      <!-- One line: the space between the two is kept. -->
      <span :id="`${testid}-problem`" :data-testid="`${testid}-problem`">{{ problem }}</span> <button
        v-if="swapped"
        type="button"
        class="link"
        :data-testid="`${testid}-swap`"
        @click="swap(swapped)"
      >
        Swap to {{ swapped.start }}–{{ swapped.stop }}
      </button>
    </p>
    <p v-else-if="region" class="hint" :data-testid="`${testid}-argument`">Positions in the results are those of the whole record.</p>
  </fieldset>
</template>
