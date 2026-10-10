<script setup lang="ts">
// The candidate tray (S14, design §11.4, REQ-14; docs/web/ncbi_ui_mapping.md "候補（S14）"): the
// HSPs collected from the results of completed runs, in the user's order, with a note each, the
// way back to each one's result, the runs they come from (Origins), and the two outputs made from
// the selected candidates: the original residues of their records (FASTA) and their gapped
// alignments (the aligned rows of their HSP records), kept in two files. Everything acts through
// the tray (application/candidates.ts); the values shown are the HSP record's coordinates and the
// outfmt 6 row's fields as written.
import { computed, nextTick, onMounted, ref, watch } from 'vue';
import type { Candidate, CandidateTray, TrayOrder, TrayState } from '../application/candidates';
import type { HspId } from '../application/results';
import { interval, spanOn } from '../domain/coordinates';
import { recordLabel, type ExtractionRegion } from '../domain/extraction';
import { programById, type InputRole } from '../domain/programs';
import CommandText from './CommandText.vue';
import { formatBytes, formatCount, formatCounted, formatDateTime } from './format';
import { useNarrow } from './useNarrow';
import VirtualRows from './VirtualRows.vue';

const props = defineProps<{ candidates: CandidateTray; state: TrayState }>();
const emit = defineEmits<{ reveal: [id: HspId] }>();

const narrow = useNarrow();
/** A row is two lines (the values; the note and the actions), on a phone four (W4 judgments 16-23: no sideways page scroll). */
const ROW_PX = computed(() => (narrow.value ? 168 : 68));
/** A phone shows three rows and half of the fourth, so that the cut row says the list scrolls (screen review L3). */
const MAX_ROWS = computed(() => (narrow.value ? 3.5 : 8));
const root = ref<HTMLElement>();
const heading = ref<HTMLElement>();
const gutter = ref(0);
const scroller = ref<HTMLElement>();

const count = computed(() => props.state.candidates.length);
const chosen = computed(() => props.state.candidates.filter((candidate) => props.state.selected.has(candidate.key)));
const allSelected = computed(() => count.value > 0 && props.state.selected.size === count.value);
const someSelected = computed(() => props.state.selected.size > 0 && !allSelected.value);

// --- what a row shows ----------------------------------------------------------------------------

const programLabel = (candidate: Candidate) => programById(candidate.run.program).label;
const queryName = (candidate: Candidate) => recordLabel(candidate.query.id, 'query', candidate.query.position);
const subjectName = (candidate: Candidate) => recordLabel(candidate.subject.id, 'subject', candidate.subject.position);
/** The HSP record's subject coordinates, the smaller first, as the Alignments' "Range n: a to b". */
function rangeText(candidate: Candidate): string {
  const { from, to } = interval(candidate.coordinates.s_start, candidate.coordinates.s_end);
  return `${from} to ${to}`;
}
/** The hit's strand on the subject record, as the extracted sequence's header says it (`hit_strand`). */
const strand = (candidate: Candidate) => spanOn(candidate.coordinates, 'subject', candidate.subject.kind).strand;
const runTitle = (candidate: Candidate) => `Run ${candidate.run.number}${candidate.run.title === undefined ? '' : `: ${candidate.run.title}`}`;
const hspText = (candidate: Candidate) => `HSP ${candidate.id.qIdx + 1}.${candidate.id.rank + 1}`;

/**
 * The Subject and Run columns are as wide as their longest values at the table's type (within the
 * bounds in styles.css; a longer value ends in an ellipsis and has its title), as the Descriptions'
 * Subject ID column is (W4b screen review L8), so that a short Run does not take the width of an ID
 * (S14 screen review L4). Only the longest values by their count of letters are measured.
 */
const subjectWidth = ref<string>();
const runWidth = ref<string>();
let measurer: CanvasRenderingContext2D | null | undefined;
function widest(texts: readonly string[]): string | undefined {
  const element = scroller.value;
  if (element === undefined || texts.length === 0) return undefined;
  measurer ??= document.createElement('canvas').getContext('2d');
  if (measurer === null) return undefined;
  const style = getComputedStyle(element);
  measurer.font = `${style.fontStyle} ${style.fontWeight} ${style.fontSize} ${style.fontFamily}`;
  const longest = Math.max(...texts.map((text) => text.length));
  let width = 0;
  let measured = 0;
  for (const text of texts) {
    if (text.length < longest - 4) continue;
    width = Math.max(width, measurer.measureText(text).width);
    if (++measured === 64) break;
  }
  return `${Math.ceil(width) + 1}px`;
}
function measureColumns(): void {
  subjectWidth.value = widest(props.state.candidates.map(subjectName));
  runWidth.value = widest(props.state.candidates.map((candidate) => `Run ${candidate.run.number}${candidate.run.title ? ` ${candidate.run.title}` : ''}`));
}
onMounted(measureColumns);
watch(() => props.state.candidates, measureColumns, { flush: 'post' });

// --- order, selection, notes, removal ------------------------------------------------------------

/** The order last applied; "custom" after a move, or after additions to a sorted tray. */
const order = ref<TrayOrder | 'custom'>('added');
/** Counts the sorts: after each one the list shows its first rows (a move keeps the rows where they are). */
const sorts = ref(0);
const ORDERS: readonly { readonly value: TrayOrder; readonly label: string }[] = [
  { value: 'added', label: 'Order added' },
  { value: 'run', label: 'Run' },
  { value: 'subject', label: 'Subject' },
];
watch(
  () => props.state.candidates.length,
  (now, before) => {
    if (now > before && order.value !== 'added') order.value = 'custom';
  },
);

function sort(value: string): void {
  const chosenOrder = ORDERS.find((item) => item.value === value)?.value;
  if (chosenOrder === undefined) return;
  props.candidates.sortBy(chosenOrder);
  order.value = chosenOrder;
  sorts.value++;
}

/** Focuses a control of a row once the rows are drawn again (after a move or a removal). */
async function focusControl(testid: string): Promise<void> {
  await nextTick();
  const element = root.value?.querySelector<HTMLElement>(`[data-testid="${testid}"]`);
  if (element !== null && element !== undefined) element.focus();
  else heading.value?.focus();
}

/** Moves a candidate up or down one place; the focus stays on the same button of the moved row. */
function move(position: number, by: -1 | 1): void {
  const candidate = props.state.candidates[position];
  const to = position + by;
  if (candidate === undefined || to < 0 || to >= count.value) return;
  props.candidates.move(candidate.key, to);
  order.value = 'custom';
  void focusControl(`candidate-${by < 0 ? 'up' : 'down'}-${to + 1}`);
}

function remove(position: number): void {
  const candidate = props.state.candidates[position];
  if (candidate === undefined) return;
  // The rows are drawn again at the next tick: the row that takes this place, else the one before it, keeps the focus.
  const next = Math.min(position, count.value - 2);
  props.candidates.remove([candidate.key]);
  void focusControl(`candidate-remove-${next + 1}`);
}

function removeSelected(): void {
  props.candidates.remove([...props.state.selected]);
  void focusControl('candidates-select-all');
}

function select(candidate: Candidate, event: Event): void {
  props.candidates.select([candidate.key], (event.target as HTMLInputElement).checked);
}

function setNote(candidate: Candidate, event: Event): void {
  props.candidates.setNote(candidate.key, (event.target as HTMLInputElement).value);
}

// --- origins -------------------------------------------------------------------------------------

/** The runs of the candidates (the tray reads them from the coordinator's record of each run). */
const origins = computed(() => {
  void props.state.candidates;
  return props.candidates.origins();
});

// --- extraction ----------------------------------------------------------------------------------

const role = ref<InputRole>('subject');
const regionKind = ref<ExtractionRegion['kind']>('hit');
const flankText = ref({ left: '100', right: '100' });
const join = ref<'separate' | 'spanning'>('separate');

/** A flank is a whole number of the record's unit, 0 or more. */
const flank = (text: string): number | undefined => (/^\d{1,9}$/.test(text.trim()) ? Number(text.trim()) : undefined);
const flanksInvalid = computed(
  () => regionKind.value === 'flanked' && (flank(flankText.value.left) === undefined || flank(flankText.value.right) === undefined),
);
/** The unit of the flanks: the records' unit when the selected candidates' records share one. */
const flankUnit = computed(() => {
  const units = new Set(chosen.value.map((candidate) => candidate[role.value].unit));
  return units.size === 1 ? [...units][0]! : 'nt or aa, the unit of each record';
});
const busy = computed(() => props.state.busy !== undefined);

function region(): ExtractionRegion {
  if (regionKind.value === 'flanked') {
    return { kind: 'flanked', flanks: { left: flank(flankText.value.left)!, right: flank(flankText.value.right)! } };
  }
  return { kind: regionKind.value };
}

async function downloadSequences(): Promise<void> {
  if (busy.value || chosen.value.length === 0 || flanksInvalid.value) return;
  await props.candidates.extract({ role: role.value, region: region(), join: join.value });
}

async function downloadAlignments(): Promise<void> {
  if (busy.value || chosen.value.length === 0) return;
  await props.candidates.exportAlignments();
}

/** Lists in the summary show their first items; the file's headers hold every one. */
const LISTED = 20;
const intervalWords = (iv: { readonly from: number; readonly to: number }, unit: string) => `${iv.from} to ${iv.to} ${unit}`;
</script>

<template>
  <section ref="root" class="candidates" data-testid="candidates">
    <h2 ref="heading" tabindex="-1">Candidates</h2>
    <p v-if="state.message" class="error" role="alert" data-testid="candidates-message">{{ state.message }}</p>
    <p v-if="count === 0" class="muted" data-testid="candidates-empty">
      No candidates yet. In the results, mark rows of the Descriptions and press “Add to candidates”, or press “Add to candidates” at a Range of
      the Alignments or in the Dot Plot's popup. Only HSPs of completed runs can be added.
    </p>
    <template v-else>
      <div class="result-table candidate-table">
        <div class="tool-band candidates-tools">
          <label class="check">
            <input
              type="checkbox"
              :checked="allSelected"
              :indeterminate="someSelected"
              data-testid="candidates-select-all"
              @change="candidates.selectAll(($event.target as HTMLInputElement).checked)"
            />
            select all
          </label>
          <span class="muted" data-testid="candidates-selected" aria-live="polite">
            {{ formatCount(state.selected.size) }} of {{ formatCounted(count, 'candidate') }} selected
          </span>
          <label class="candidates-sort">
            Sort by
            <select :value="order" data-testid="candidates-sort" @change="sort(($event.target as HTMLSelectElement).value)">
              <option v-if="order === 'custom'" value="custom" disabled>Your order</option>
              <option v-for="item in ORDERS" :key="item.value" :value="item.value">{{ item.label }}</option>
            </select>
          </label>
          <button type="button" :disabled="state.selected.size === 0" data-testid="candidates-remove-selected" @click="removeSelected">
            Remove selected
          </button>
        </div>
        <div ref="scroller" class="table-scroll" :style="{ '--row-gutter': `${gutter}px`, '--subject-width': subjectWidth, '--run-width': runWidth }">
          <div v-if="!narrow" class="table-head candidate-grid marked-head" role="row">
            <span class="num-head" data-col="n">#</span>
            <span data-col="run">Run</span>
            <span data-col="program">Program</span>
            <span data-col="query">Query</span>
            <span data-col="subject">Subject</span>
            <span data-col="range">Range on subject</span>
            <span data-col="strand">Strand</span>
            <span class="num-head" data-col="evalue">E value</span>
            <span class="num-head" data-col="bitscore">Bit score</span>
          </div>
          <VirtualRows
            :count="count"
            :row-px="ROW_PX"
            :max-rows="MAX_ROWS"
            :order-key="String(sorts)"
            label="Candidates"
            testid="candidate-list"
            @gutter="gutter = $event"
          >
            <template #row="{ position }">
              <div
                v-for="candidate in [state.candidates[position]!]"
                :key="candidate.key"
                class="candidate-row"
                :class="{ narrow }"
                :data-testid="`candidate-${position + 1}`"
                :data-key="candidate.key"
                :data-run="candidate.run.number"
                :data-hsp="`${candidate.id.qIdx}:${candidate.id.rank}`"
              >
                <label class="row-mark">
                  <input
                    type="checkbox"
                    :checked="state.selected.has(candidate.key)"
                    :aria-label="`Select candidate ${position + 1}`"
                    :data-testid="`candidate-mark-${position + 1}`"
                    @change="select(candidate, $event)"
                  />
                </label>
                <!-- A phone shows #, Subject and Range first, then the note and the actions; the other values under them. -->
                <div v-if="narrow" class="candidate-card">
                  <div class="candidate-key">
                    <span class="num" data-field="n">{{ position + 1 }}</span>
                    <span class="candidate-subject" data-field="subject" :title="subjectName(candidate)">{{ subjectName(candidate) }}</span>
                    <span class="num" data-field="range" :title="`${rangeText(candidate)} of ${formatCount(candidate.subject.length)} ${candidate.subject.unit}`">{{
                      rangeText(candidate)
                    }}</span>
                  </div>
                  <input
                    type="text"
                    class="candidate-note"
                    :value="candidate.note"
                    placeholder="Note"
                    :aria-label="`Note for candidate ${position + 1}`"
                    :data-testid="`candidate-note-${position + 1}`"
                    @input="setNote(candidate, $event)"
                  />
                  <span class="candidate-actions">
                    <button type="button" :aria-disabled="position === 0" :data-testid="`candidate-up-${position + 1}`" @click="move(position, -1)">Up</button>
                    <button type="button" :aria-disabled="position === count - 1" :data-testid="`candidate-down-${position + 1}`" @click="move(position, 1)">
                      Down
                    </button>
                    <button type="button" :data-testid="`candidate-reveal-${position + 1}`" @click="emit('reveal', candidate.id)">Show in results</button>
                    <button type="button" :data-testid="`candidate-remove-${position + 1}`" @click="remove(position)">Remove</button>
                  </span>
                  <p class="candidate-rest muted small" :title="runTitle(candidate)">
                    <span data-field="run">{{ runTitle(candidate) }}</span> · <span data-field="program">{{ programLabel(candidate) }}</span> ·
                    {{ hspText(candidate) }} · query <span data-field="query">{{ queryName(candidate) }}</span> · strand
                    <span data-field="strand">{{ strand(candidate) }}</span> · E value <span data-field="evalue">{{ candidate.row.evalue }}</span> · bit
                    score <span data-field="bitscore">{{ candidate.row.bitscore }}</span>
                  </p>
                </div>
                <template v-else>
                  <div class="candidate-grid candidate-cells">
                    <span class="num" data-field="n" :title="hspText(candidate)">{{ position + 1 }}</span>
                    <span data-field="run" :title="runTitle(candidate)"
                      >Run {{ candidate.run.number }}<span v-if="candidate.run.title" class="muted">{{ ` ${candidate.run.title}` }}</span></span
                    >
                    <span data-field="program">{{ programLabel(candidate) }}</span>
                    <span data-field="query" :title="queryName(candidate)">{{ queryName(candidate) }}</span>
                    <span data-field="subject" :title="subjectName(candidate)">{{ subjectName(candidate) }}</span>
                    <span data-field="range" :title="`${rangeText(candidate)} of ${formatCount(candidate.subject.length)} ${candidate.subject.unit}`">{{
                      rangeText(candidate)
                    }}</span>
                    <span data-field="strand" :class="`strand-${strand(candidate)}`">{{ strand(candidate) }}</span>
                    <span class="num" data-field="evalue" :title="candidate.row.evalue">{{ candidate.row.evalue }}</span>
                    <span class="num" data-field="bitscore" :title="candidate.row.bitscore">{{ candidate.row.bitscore }}</span>
                  </div>
                  <div class="candidate-tools">
                    <input
                      type="text"
                      class="candidate-note"
                      :value="candidate.note"
                      placeholder="Note"
                      :aria-label="`Note for candidate ${position + 1}`"
                      :data-testid="`candidate-note-${position + 1}`"
                      @input="setNote(candidate, $event)"
                    />
                    <span class="candidate-actions">
                      <button type="button" :aria-disabled="position === 0" :data-testid="`candidate-up-${position + 1}`" @click="move(position, -1)">Up</button>
                      <button
                        type="button"
                        :aria-disabled="position === count - 1"
                        :data-testid="`candidate-down-${position + 1}`"
                        @click="move(position, 1)"
                      >
                        Down
                      </button>
                      <button type="button" :data-testid="`candidate-reveal-${position + 1}`" @click="emit('reveal', candidate.id)">Show in results</button>
                      <button type="button" :data-testid="`candidate-remove-${position + 1}`" @click="remove(position)">Remove</button>
                    </span>
                  </div>
                </template>
              </div>
            </template>
          </VirtualRows>
        </div>
      </div>

      <section class="extract" data-testid="extract-form">
        <h3>Download the selected candidates</h3>
        <p class="muted small" data-testid="extract-chosen">{{ formatCounted(chosen.length, 'candidate') }} selected.</p>
        <div class="extract-options">
          <fieldset class="radios">
            <legend>Sequence</legend>
            <label class="radio"><input v-model="role" type="radio" value="subject" data-testid="extract-role-subject" /> Subject</label>
            <label class="radio"><input v-model="role" type="radio" value="query" data-testid="extract-role-query" /> Query</label>
          </fieldset>
          <fieldset class="radios">
            <legend>Region</legend>
            <label class="radio"><input v-model="regionKind" type="radio" value="hit" data-testid="extract-region-hit" /> Hit region</label>
            <label class="radio"
              ><input v-model="regionKind" type="radio" value="flanked" data-testid="extract-region-flanked" /> Hit region with flanks</label
            >
            <div v-if="regionKind === 'flanked'" class="extract-flanks">
              <label>Left <input v-model="flankText.left" type="text" inputmode="numeric" size="7" data-testid="extract-flank-left" /></label>
              <label>Right <input v-model="flankText.right" type="text" inputmode="numeric" size="7" data-testid="extract-flank-right" /></label>
              <span class="muted small">{{ flankUnit }}</span>
              <p v-if="flanksInvalid" class="error small" data-testid="extract-flank-error">A flank is a whole number, 0 or more.</p>
            </div>
            <label class="radio"><input v-model="regionKind" type="radio" value="whole" data-testid="extract-region-whole" /> Complete sequence</label>
          </fieldset>
          <fieldset class="radios" :disabled="regionKind === 'whole'">
            <legend>Several HSPs on one record</legend>
            <label class="radio"><input v-model="join" type="radio" value="separate" data-testid="extract-join-separate" /> Separate sequences</label>
            <label class="radio"><input v-model="join" type="radio" value="spanning" data-testid="extract-join-spanning" /> One region spanning them</label>
            <p v-if="regionKind === 'whole'" class="muted small">A complete sequence is written once per record.</p>
          </fieldset>
        </div>
        <p class="hint" data-testid="extract-note">
          The sequences are the records' own letters as the input files have them, on each record's own strand (the header's hit_strand gives
          the hit's strand); lengths are in nt for nucleotides and aa for proteins.
        </p>
        <div class="extract-actions">
          <button
            type="button"
            class="primary-action"
            :disabled="busy || chosen.length === 0 || flanksInvalid"
            data-testid="extract-download"
            @click="downloadSequences"
          >
            Download FASTA
          </button>
          <button type="button" :disabled="busy || chosen.length === 0" data-testid="extract-aligned" @click="downloadAlignments">
            Download aligned sequences (FASTA)
          </button>
        </div>
        <p class="hint">
          Aligned sequences are the HSP records' aligned rows as the search wrote them, gaps included, in the search's direction: a separate file,
          not taken from the input files.
        </p>
        <div aria-live="polite">
          <p v-if="state.busy" class="muted" data-testid="extract-busy">Writing the file…</p>
          <div v-else-if="state.last" class="extract-summary" data-testid="extract-summary" :data-output="state.last.output">
            <template v-if="state.last.output === 'sequences'">
              <p>
                Saved {{ state.last.fileName }}: {{ formatCounted(state.last.sequences, 'sequence') }} of {{ state.last.role }} records from
                {{ formatCounted(state.last.candidates, 'candidate') }}, <span class="nowrap">{{ formatBytes(state.last.bytes) }}</span>.
              </p>
              <template v-if="state.last.clipped.length > 0">
                <p>
                  {{ formatCounted(state.last.clipped.length, 'sequence') }} cut at a record's end (the header's requested= gives the range asked
                  for):
                </p>
                <ul>
                  <li v-for="(clip, i) in state.last.clipped.slice(0, LISTED)" :key="i" data-testid="extract-clipped">
                    {{ clip.name }} (run {{ clip.runNumber }}, {{ clip.role }} record {{ clip.position + 1 }}, HSP {{ clip.hsps.join(', ') }}):
                    requested <span class="nowrap">{{ intervalWords(clip.requested, clip.unit) }}</span>, written
                    <span class="nowrap">{{ intervalWords(clip.actual, clip.unit) }}</span> of
                    <span class="nowrap">{{ clip.recordLength }} {{ clip.unit }}</span>
                  </li>
                  <li v-if="state.last.clipped.length > LISTED" class="muted">and {{ formatCount(state.last.clipped.length - LISTED) }} more</li>
                </ul>
              </template>
              <template v-if="state.last.unknownStrand.length > 0">
                <p>Strand not decided:</p>
                <ul>
                  <li v-for="(line, i) in state.last.unknownStrand.slice(0, LISTED)" :key="i" data-testid="extract-unknown-strand">{{ line }}</li>
                  <li v-if="state.last.unknownStrand.length > LISTED" class="muted">and {{ formatCount(state.last.unknownStrand.length - LISTED) }} more</li>
                </ul>
              </template>
            </template>
            <template v-else>
              <p>
                Saved {{ state.last.fileName }}: {{ formatCounted(state.last.alignments, 'alignment') }} from
                {{ formatCounted(state.last.candidates, 'candidate') }}, <span class="nowrap">{{ formatBytes(state.last.bytes) }}</span>.
              </p>
              <template v-if="state.last.missing.length > 0">
                <p>Not written:</p>
                <ul>
                  <li v-for="(line, i) in state.last.missing.slice(0, LISTED)" :key="i" data-testid="extract-missing">{{ line }}</li>
                  <li v-if="state.last.missing.length > LISTED" class="muted">and {{ formatCount(state.last.missing.length - LISTED) }} more</li>
                </ul>
              </template>
            </template>
          </div>
        </div>
      </section>

      <section class="candidate-origins" data-testid="candidate-origins">
        <h3>Origins</h3>
        <ul class="origin-list">
          <li v-for="origin in origins" :key="origin.runId" :data-testid="`candidate-origin-${origin.number}`">
            <h4>
              Run {{ origin.number }}<template v-if="origin.title">{{ ` · ${origin.title}` }}</template>
              <span class="muted small">{{ formatCounted(origin.candidates, 'candidate') }}</span>
            </h4>
            <dl class="details-list">
              <dt>Program</dt>
              <dd data-detail="program">{{ programById(origin.program).label }}</dd>
              <dt>Options</dt>
              <dd data-detail="options"><CommandText :text="origin.options.length === 0 ? 'defaults' : origin.options.join(' ')" /></dd>
              <dt>Query</dt>
              <dd data-detail="query">
                {{ origin.query.name }}<br /><span class="muted small">SHA-256 {{ origin.query.sha256 }}</span>
              </dd>
              <dt>Subject</dt>
              <dd data-detail="subject">
                {{ origin.subject.name }}<br /><span class="muted small">SHA-256 {{ origin.subject.sha256 }}</span>
              </dd>
              <dt>Engine build</dt>
              <dd data-detail="build">{{ origin.engineBuild ?? 'not recorded' }}</dd>
              <dt>Ended</dt>
              <dd data-detail="ended">{{ origin.endedAt === undefined ? 'not recorded' : formatDateTime(origin.endedAt) }}</dd>
            </dl>
          </li>
        </ul>
      </section>
    </template>
  </section>
</template>
