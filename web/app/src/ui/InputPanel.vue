<script setup lang="ts">
// A role's block ("Enter Query Sequence", "Enter Subject Sequence"), in the order of NCBI's
// two-sequence page: the paste box (FASTA, or accessions to fetch) with its subrange beside it, "Or, upload file", the
// role's genetic code (BLASTX, TBLASTN, TBLASTX), the Job Title (query), then LOSAT's own parts: the
// inputs with their records, and Combined / Separate. A file dropped anywhere on the block is
// opened.
import { computed, nextTick, ref, watch } from 'vue';
import type { InputRole } from 'losat/domain';
import type { DraftState, SearchDraft } from '../application/draft';
import { groupSite, raceChoice } from '../application/source-race';
import {
  accessionLines,
  FETCH_SOURCES,
  fetchGroups,
  RACE_RECORDS,
  servingSources,
  sourceLabel,
  sourceList,
  type FetchChoice,
  type FetchGroup,
  type FetchSourceId,
} from '../domain/accession';
import { sectionsAt, writtenValue } from '../domain/parameters';
import { DB_BOOL_FILTER_FIELDS, DB_FILE_FILTER_FIELDS, DB_TEXT_FILTER_FIELDS } from '../domain/blast-db';
import { programById, querySourcesOf, residueUnit, sequenceKind } from '../domain/programs';
import DatabasePicker from './DatabasePicker.vue';
import { formatBytes, optionHelp } from './format';
import HelpTip from './HelpTip.vue';
import { searchFormModel } from './form-model';
import LosatLoader from './LosatLoader.vue';
import ParameterField from './ParameterField.vue';
import RegionPicker from './RegionPicker.vue';
import SourceCard from './SourceCard.vue';

const props = defineProps<{ draft: SearchDraft; state: DraftState; role: InputRole }>();
const formModel = computed(() => searchFormModel(props.draft, props.state));
const title = computed(() => (props.role === 'query' ? 'Query' : 'Subject'));
const roleState = computed(() => props.state[props.role]);
const program = computed(() => programById(props.state.program));
const kind = computed(() => sequenceKind(program.value, props.role));
const dragging = ref(false);
const fileInput = ref<HTMLInputElement>();
const openButton = ref<HTMLButtonElement>();
const sourceItems = ref<HTMLElement>();
const fetchButton = ref<HTMLButtonElement>();
const fetchGroupList = ref<HTMLElement>();

/**
 * Removes an input; the focus goes to the input that takes its place in the list, or to "Open FASTA
 * files…" when it was the last one, and not to the page's body (WA2-AX-04).
 */
async function removeSource(key: string): Promise<void> {
  const at = roleState.value.sources.findIndex((source) => source.key === key);
  props.draft.removeSource(props.role, key);
  await nextTick();
  const buttons = sourceItems.value?.querySelectorAll<HTMLButtonElement>('.source-remove');
  (buttons?.[Math.min(at, buttons.length - 1)] ?? openButton.value)?.focus();
}

/**
 * Puts each accession of the box on its own line (WA2-OP-04): the pasted input's card goes, and the
 * focus goes to the Fetch that now takes them, not to the page's body.
 */
async function splitAccessions(): Promise<void> {
  props.draft.splitAccessions(props.role);
  await nextTick();
  (fetchButton.value ?? fetchGroupList.value?.querySelector<HTMLButtonElement>('button'))?.focus();
}

/** The box's accessions, when it holds accessions only: "Fetch" brings their FASTA into it. */
const accessions = computed(() => accessionLines(roleState.value.paste));
const fetchState = computed(() => roleState.value.fetch);
const busy = (role: InputRole) => {
  const state = props.state[role].fetch?.state;
  return state === 'fetching' || state === 'racing';
};
const fetchBusy = computed(() => busy(props.role));
/** The choice is shared by both boxes: it stays while either fetches, so that each line names where its fetch goes. */
const choiceLocked = computed(() => busy('query') || busy('subject'));
/** The sources that serve the box's accessions, and where "Fastest" sends them (once the race says). */
const serving = computed(() => (accessions.value === undefined ? [] : servingSources(kind.value, accessions.value)));
const fastest = computed(() =>
  serving.value.length === 1 ? serving.value[0] : raceChoice(props.state.race[kind.value], serving.value),
);
/** A box that no one site serves: its groups, each with a Fetch of its own that names its site (line WF-5). */
const groups = computed(() =>
  accessions.value === undefined || serving.value.length > 0 ? [] : fetchGroups(kind.value, accessions.value),
);
/** The sites picked in the groups' selects, by the group's sites, for this text of the box. */
const picked = ref<Record<string, FetchSourceId>>({});
watch(
  () => roleState.value.paste,
  () => (picked.value = {}),
);
const groupKey = (group: FetchGroup) => group.sources.join(',');
const groupCame = (group: FetchGroup) => group.lines.every((line) => roleState.value.fetched?.[line] !== undefined);
/** A group's site: where it came from, else the one picked, else where it starts (groupSite). */
function groupSource(group: FetchGroup): FetchSourceId {
  const came = roleState.value.fetched?.[group.lines[0]!]?.source;
  if (came !== undefined) return came;
  const pick = picked.value[groupKey(group)];
  return pick !== undefined && group.sources.includes(pick)
    ? pick
    : groupSite(group.sources, props.state.fetchSource, props.state.race[kind.value]);
}
/** "C_AA004835.1", "P12345 and 1 more": each group names its first accession, which no other group has. */
const groupName = (group: FetchGroup) =>
  group.accessions.length === 1 ? group.accessions[0]! : `${group.accessions[0]!} and ${group.accessions.length - 1} more`;
function pick(group: FetchGroup, source: FetchSourceId): void {
  picked.value = { ...picked.value, [groupKey(group)]: source };
}
function onGroupSource(group: FetchGroup, event: Event): void {
  pick(group, (event.target as HTMLSelectElement).value as FetchSourceId);
  props.draft.clearFetchFailure(props.role);
}
function onGroupFetch(group: FetchGroup, index: number): void {
  // The site pressed stays on the row while the group is fetched, whatever the race says meanwhile.
  const source = groupSource(group);
  pick(group, source);
  props.draft.fetchGroup(props.role, index, source);
}
/** The in-place sentence: where Fetch sends the accessions. */
const fetchNote = computed(() => {
  const noun = accessions.value?.length === 1 ? 'accession' : 'accessions';
  const choice = props.state.fetchSource;
  if (serving.value.length === 0) return 'No one site has all of these accessions';
  if (choice !== 'fastest') return `Sends the ${noun} to ${sourceLabel(choice)}`;
  if (fastest.value !== undefined) return `Sends the ${noun} to ${sourceLabel(fastest.value)}`;
  return `Sends the ${noun} to the fastest of ${sourceList(serving.value, 'and')}`;
});
const fetchHelp = computed(() => {
  const records = RACE_RECORDS[kind.value];
  const record = records.ebi === records.ncbi ? records.ncbi : `${records.ncbi}, at UniProt ${records.ebi}`;
  return (
    'Fetch sends only the accessions, to the site chosen here, and writes the FASTA it gets into the box. Sequences, ' +
    'files and results stay in this browser. EBI answers from ENA for nucleotides and from UniProt for proteins, CNCB ' +
    `from GenBase. Fastest: the first fetch asks each site for one public record (${record}) and sends the ` +
    'accessions to the first that answers; this tab keeps that order.'
  );
});

const readySources = computed(() => roleState.value.sources.filter((source) => source.status === 'ready'));
const regionRecord = computed(() => {
  void props.state;
  return props.draft.regionRecord(props.role);
});
const region = computed(() => {
  void props.state;
  return props.draft.region(props.role);
});

/** Whether this panel is showing the local BLAST database picker (`-db`) instead of FASTA inputs. */
const isDatabase = computed(() => props.role === 'subject' && props.state.subjectSource === 'database');
/** Query input sources supported by the current program (`fasta`, plus `msa`/`pssm` for psiblast, `pssm` for tblastn). */
const querySources = computed(() => (props.role === 'query' ? querySourcesOf(props.state.program) : []));
/** Whether this panel is showing a PSSM (`-in_pssm`) or MSA (`-in_msa`) file picker instead of FASTA inputs. */
const isQueryFile = computed(() => props.role === 'query' && props.state.querySource !== 'fasta');
const queryFileFlag = computed<'-in_pssm' | '-in_msa'>(() => (props.state.querySource === 'msa' ? '-in_msa' : '-in_pssm'));
const queryFileInput = ref<HTMLInputElement>();

function onQueryFile(event: Event): void {
  const input = event.target as HTMLInputElement;
  const file = input.files?.[0];
  if (file !== undefined) props.draft.setQueryFile(file);
  input.value = '';
}

const legendTitle = computed(() => {
  if (isDatabase.value) return 'Choose Subject Database';
  if (isQueryFile.value) return props.state.querySource === 'msa' ? 'Choose Query MSA' : 'Choose Query PSSM';
  return `Enter ${title.value} Sequence`;
});

/** Database filter options accepted by the current program's engine `describe`. */
const describedFlags = computed(() => {
  const params = props.state.description?.parameters;
  return params === undefined ? undefined : new Set(params.map((p) => p.flag));
});
const dbFilters = computed(() =>
  [...DB_TEXT_FILTER_FIELDS, ...DB_BOOL_FILTER_FIELDS, ...DB_FILE_FILTER_FIELDS].filter(
    (f) => describedFlags.value === undefined || describedFlags.value.has(f.flag),
  ),
);
function dbFilterId(flag: string): string {
  return `db-filter-${flag.slice(1)}`;
}
function isDbFilterSet(field: { readonly flag: string; readonly kind: string }): boolean {
  if (field.kind === 'file') return props.state.database.filterFiles[field.flag] !== undefined;
  const val = props.state.database.filterValues[field.flag];
  return field.kind === 'flag' ? val === true : typeof val === 'string' && val.trim() !== '';
}
const activeDbFilterCount = computed(() => dbFilters.value.filter(isDbFilterSet).length);
const dbFilterFileInputs = ref<Record<string, HTMLInputElement | null>>({});

function setDbFilterInputRef(flag: string, el: unknown): void {
  dbFilterFileInputs.value[flag] = el instanceof HTMLInputElement ? el : null;
}

function onDbFilterFile(flag: string, event: Event): void {
  const input = event.target as HTMLInputElement;
  const file = input.files?.[0];
  if (file !== undefined) props.draft.setDatabaseFilterFile(flag, file);
  input.value = '';
}

function dbOptionHelp(flag: string): string {
  const option = props.state.description?.parameters.find((p) => p.flag === flag);
  return optionHelp(option?.help);
}

/** The role's genetic code (the query's of BLASTX and TBLASTX, the subjects' of TBLASTN and TBLASTX). */
const codeFields = computed(() =>
  sectionsAt(program.value, props.state.description?.parameters, props.role).flatMap((section) => section.fields),
);
/** A subject genetic code other than the default: LOSAT's approved exception applies only to local `-subject`, not `-db`. */
const subjectCodeNote = computed(() => {
  if (isDatabase.value) return false;
  const field = codeFields.value.find((f) => f.flag === '-db_gencode');
  if (props.role !== 'subject' || field === undefined) return false;
  const option = props.state.description?.parameters.find((o) => o.flag === field.flag);
  return writtenValue(field, props.state.values[props.state.program][field.flag], option) !== undefined;
});

function onPaste(event: Event): void {
  props.draft.setPaste(props.role, (event.target as HTMLTextAreaElement).value);
}

function onFetchSource(event: Event): void {
  props.draft.setFetchSource((event.target as HTMLSelectElement).value as FetchChoice);
}

function onFiles(event: Event): void {
  const input = event.target as HTMLInputElement;
  props.draft.addFiles(props.role, [...(input.files ?? [])]);
  input.value = '';
}

function onDrop(event: DragEvent): void {
  dragging.value = false;
  if (isDatabase.value) return;
  const files = [...(event.dataTransfer?.files ?? [])];
  if (isQueryFile.value) {
    if (files[0] !== undefined) props.draft.setQueryFile(files[0]);
    return;
  }
  if (files.length > 0) props.draft.addFiles(props.role, files);
}
</script>

<template>
  <fieldset
    class="search-block input-panel"
    :class="{ dragging }"
    :data-testid="`${role}-panel`"
    @dragover.prevent="dragging = true"
    @dragleave.self="dragging = false"
    @drop.prevent.stop="onDrop"
  >
    <legend>{{ legendTitle }}</legend>

    <div v-if="role === 'query' && querySources.length > 1" class="subject-source" role="radiogroup" aria-label="Query source" data-testid="query-source">
      <label class="subject-source-option">
        <input
          type="radio"
          name="query-source"
          value="fasta"
          :checked="state.querySource === 'fasta'"
          data-testid="query-source-fasta"
          @change="draft.setQuerySource('fasta')"
        />
        FASTA sequences
      </label>
      <label v-if="querySources.includes('msa')" class="subject-source-option">
        <input
          type="radio"
          name="query-source"
          value="msa"
          :checked="state.querySource === 'msa'"
          data-testid="query-source-msa"
          @change="draft.setQuerySource('msa')"
        />
        MSA <span class="flag">-in_msa</span>
      </label>
      <label v-if="querySources.includes('pssm')" class="subject-source-option">
        <input
          type="radio"
          name="query-source"
          value="pssm"
          :checked="state.querySource === 'pssm'"
          data-testid="query-source-pssm"
          @change="draft.setQuerySource('pssm')"
        />
        PSSM <span class="flag">-in_pssm</span>
      </label>
    </div>

    <div v-if="role === 'subject'" class="subject-source" role="radiogroup" aria-label="Subject source" data-testid="subject-source">
      <label class="subject-source-option">
        <input
          type="radio"
          name="subject-source"
          value="fasta"
          :checked="!isDatabase"
          data-testid="subject-source-fasta"
          @change="draft.setSubjectSource('fasta')"
        />
        FASTA sequences
      </label>
      <label class="subject-source-option">
        <input
          type="radio"
          name="subject-source"
          value="database"
          :checked="isDatabase"
          data-testid="subject-source-database"
          @change="draft.setSubjectSource('database')"
        />
        Database <span class="flag">-db</span>
      </label>
    </div>

    <template v-if="isQueryFile">
      <div class="subject-db-block" data-testid="query-file-block">
        <div class="entry-label">
          <span>{{ state.querySource === 'msa' ? 'Multiple sequence alignment (FASTA)' : 'PSSM checkpoint (ASN.1)' }}</span>
          <span class="kind">{{ queryFileFlag }}</span>
          <HelpTip
            v-if="dbOptionHelp(queryFileFlag)"
            :id="`query-file-${state.querySource}`"
            :label="`About ${queryFileFlag}`"
            :text="dbOptionHelp(queryFileFlag)"
          />
        </div>
        <div class="db-picker" data-testid="query-file-picker">
          <div class="db-picker-actions">
            <button type="button" data-testid="query-file-choose" @click="queryFileInput?.click()">
              {{ state.queryFile.file ? 'Change file…' : state.querySource === 'msa' ? 'Choose MSA file…' : 'Choose PSSM file…' }}
            </button>
            <input
              ref="queryFileInput"
              class="visually-hidden"
              type="file"
              tabindex="-1"
              aria-hidden="true"
              :aria-label="state.querySource === 'msa' ? 'Choose MSA file' : 'Choose PSSM file'"
              data-testid="query-file-input"
              @change="onQueryFile"
            />
            <span v-if="!state.queryFile.file" class="hint">or drop a file here</span>
            <template v-else>
              <span class="file-name" data-testid="query-file-name">
                {{ state.queryFile.file.name }} ({{ formatBytes(state.queryFile.file.size) }})
              </span>
              <button type="button" class="link" data-testid="query-file-remove" @click="draft.setQueryFile(undefined)">Remove</button>
            </template>
          </div>
          <p class="db-privacy muted">Read locally in this browser; nothing is uploaded.</p>
        </div>
      </div>

      <template v-if="state.querySource === 'msa'">
        <div class="param-row" :class="{ changed: state.queryFile.msaMasterIdx.trim() !== '' && !state.queryFile.ignoreMsaMaster }">
          <label class="param-label" for="query-msa-master-idx">
            MSA master sequence index <span class="flag">-msa_master_idx</span>
          </label>
          <div class="param-control">
            <span class="control-tip">
              <input
                id="query-msa-master-idx"
                type="text"
                inputmode="numeric"
                autocomplete="off"
                spellcheck="false"
                placeholder="1"
                :disabled="state.queryFile.ignoreMsaMaster"
                :value="state.queryFile.msaMasterIdx"
                data-testid="query-msa-master-idx"
                @input="draft.setMsaMasterIdx(($event.target as HTMLInputElement).value)"
              />
              <HelpTip
                v-if="dbOptionHelp('-msa_master_idx')"
                id="query-msa-master-idx"
                label="About -msa_master_idx"
                :text="dbOptionHelp('-msa_master_idx')"
              />
            </span>
          </div>
        </div>
        <div class="param-row" :class="{ changed: state.queryFile.ignoreMsaMaster }">
          <label class="param-label" for="query-ignore-msa-master">
            Ignore MSA master sequence <span class="flag">-ignore_msa_master</span>
          </label>
          <div class="param-control">
            <span class="control-tip">
              <input
                id="query-ignore-msa-master"
                type="checkbox"
                :checked="state.queryFile.ignoreMsaMaster"
                data-testid="query-ignore-msa-master"
                @change="draft.setIgnoreMsaMaster(($event.target as HTMLInputElement).checked)"
              />
              <HelpTip
                v-if="dbOptionHelp('-ignore_msa_master')"
                id="query-ignore-msa-master"
                label="About -ignore_msa_master"
                :text="dbOptionHelp('-ignore_msa_master')"
              />
            </span>
          </div>
        </div>
      </template>

      <div class="param-row">
        <label class="param-label" for="job-title">Job Title</label>
        <div class="param-control fill">
          <input
            id="job-title"
            type="text"
            autocomplete="off"
            :value="state.title"
            aria-describedby="job-title-help"
            data-testid="job-title"
            @input="draft.setTitle(($event.target as HTMLInputElement).value)"
          />
          <p id="job-title-help" class="help">Enter a descriptive title for your search</p>
        </div>
      </div>
    </template>

    <template v-else-if="isDatabase">
      <div class="subject-db-block" data-testid="subject-db-block">
        <div class="entry-label">
          <span>Local BLAST database</span>
          <span class="kind">{{ kind }}</span>
        </div>
        <DatabasePicker
          :pick="state.database"
          :molecule-kind="kind"
          :context-label="program.label"
          test-id-prefix="subject-db"
          @files="draft.setSubjectDatabaseFiles($event)"
          @select="draft.selectSubjectDatabase($event)"
          @clear="draft.clearSubjectDatabase()"
        />
      </div>

      <ParameterField v-for="field in codeFields" :key="field.flag" :model="formModel" :field="field" />

      <details class="db-filters" data-testid="subject-db-filters">
        <summary data-testid="subject-db-filters-toggle">
          Database filters
          <span v-if="activeDbFilterCount > 0" class="muted" data-testid="subject-db-filters-count">({{ activeDbFilterCount }} set)</span>
        </summary>
        <div class="db-filter-rows">
          <div
            v-for="field in dbFilters"
            :key="field.flag"
            class="param-row"
            :class="{ changed: isDbFilterSet(field) }"
          >
            <component
              :is="field.kind === 'file' ? 'span' : 'label'"
              :id="field.kind === 'file' ? `${dbFilterId(field.flag)}-label` : undefined"
              :for="field.kind === 'file' ? undefined : dbFilterId(field.flag)"
              class="param-label"
            >
              {{ field.label }} <span class="flag">{{ field.flag }}</span>
            </component>
            <div class="param-control">
              <span v-if="field.kind === 'file'" class="file-pick">
                <button
                  type="button"
                  :aria-describedby="`${dbFilterId(field.flag)}-label`"
                  :data-testid="`${dbFilterId(field.flag)}-choose`"
                  @click="dbFilterFileInputs[field.flag]?.click()"
                >
                  {{ state.database.filterFiles[field.flag] ? 'Change file…' : 'Choose file…' }}
                </button>
                <input
                  :ref="(el) => setDbFilterInputRef(field.flag, el)"
                  class="visually-hidden"
                  type="file"
                  tabindex="-1"
                  aria-hidden="true"
                  :aria-label="`Choose ${field.label} file`"
                  :data-testid="`${dbFilterId(field.flag)}-input`"
                  @change="onDbFilterFile(field.flag, $event)"
                />
                <span class="hint-line">
                  <template v-if="state.database.filterFiles[field.flag]">
                    <span class="file-name" :data-testid="`${dbFilterId(field.flag)}-name`">
                      {{ state.database.filterFiles[field.flag]!.name }} ({{ formatBytes(state.database.filterFiles[field.flag]!.size) }})
                    </span>
                    <button
                      type="button"
                      class="link"
                      :data-testid="`${dbFilterId(field.flag)}-remove`"
                      @click="draft.setDatabaseFilterFile(field.flag, undefined)"
                    >
                      Remove
                    </button>
                  </template>
                  <HelpTip
                    v-if="dbOptionHelp(field.flag)"
                    :id="dbFilterId(field.flag)"
                    :label="`About ${field.flag}`"
                    :text="dbOptionHelp(field.flag)"
                  />
                </span>
              </span>
              <span v-else class="control-tip">
                <input
                  v-if="field.kind === 'flag'"
                  :id="dbFilterId(field.flag)"
                  type="checkbox"
                  :checked="state.database.filterValues[field.flag] === true"
                  :data-testid="dbFilterId(field.flag)"
                  @change="draft.setDatabaseFilter(field.flag, ($event.target as HTMLInputElement).checked)"
                />
                <input
                  v-else
                  :id="dbFilterId(field.flag)"
                  type="text"
                  autocomplete="off"
                  spellcheck="false"
                  :placeholder="field.placeholder ?? ''"
                  :value="typeof state.database.filterValues[field.flag] === 'string' ? state.database.filterValues[field.flag] : ''"
                  :data-testid="dbFilterId(field.flag)"
                  @input="draft.setDatabaseFilter(field.flag, ($event.target as HTMLInputElement).value)"
                />
                <HelpTip
                  v-if="dbOptionHelp(field.flag)"
                  :id="dbFilterId(field.flag)"
                  :label="`About ${field.flag}`"
                  :text="dbOptionHelp(field.flag)"
                />
              </span>
            </div>
          </div>
        </div>
      </details>
    </template>

    <template v-else>
    <div class="entry">
      <div class="entry-paste">
        <div class="entry-label">
          <label :for="`${role}-input`">Enter accession number(s) or FASTA sequence(s)</label>
          <span class="kind">{{ kind }}</span>
          <button type="button" class="link" :data-testid="`${role}-clear`" @click="draft.setPaste(role, '')">Clear</button>
        </div>
        <textarea
          :id="`${role}-input`"
          :value="roleState.paste"
          rows="5"
          spellcheck="false"
          autocomplete="off"
          :data-testid="`${role}-input`"
          @input="onPaste"
        />
        <div v-if="accessions" class="fetch-row" :data-testid="`${role}-fetch-row`">
          <template v-if="groups.length === 0">
            <button
              ref="fetchButton"
              type="button"
              :disabled="fetchBusy"
              :data-testid="`${role}-fetch`"
              @click="draft.fetchAccessions(role)"
            >
              Fetch
            </button>
            <label>
              from
              <select
                :value="state.fetchSource"
                :disabled="choiceLocked"
                :data-testid="`${role}-fetch-source`"
                @change="onFetchSource"
              >
                <option value="fastest">Fastest{{ fastest === undefined ? '' : ` (${sourceLabel(fastest)})` }}</option>
                <option v-for="source in FETCH_SOURCES" :key="source.id" :value="source.id">{{ source.label }}</option>
              </select>
            </label>
          </template>
          <span class="hint-line">
            <span class="hint" :data-testid="`${role}-fetch-note`">{{ fetchNote }}</span>
            <HelpTip :id="`${role}-fetch`" label="About fetching accessions" :text="fetchHelp" />
          </span>
        </div>
        <ul v-if="groups.length > 0" ref="fetchGroupList" class="fetch-groups" :data-testid="`${role}-fetch-groups`">
          <li v-for="(group, index) in groups" :key="groupKey(group)" :data-testid="`${role}-fetch-group-${index}`">
            <button
              type="button"
              :disabled="fetchBusy || groupCame(group)"
              :aria-label="`${groupCame(group) ? 'Fetched' : 'Fetch'} ${groupName(group)} from ${sourceLabel(groupSource(group))}`"
              :title="group.accessions.length > 1 ? group.accessions.join(', ') : undefined"
              :data-testid="`${role}-fetch-group-${index}-button`"
              @click="onGroupFetch(group, index)"
            >
              {{ groupCame(group) ? 'Fetched' : 'Fetch' }} {{ groupName(group) }}
            </button>
            <!-- The button's name already says where it sends. -->
            <span aria-hidden="true">from</span>
            <select
              v-if="group.sources.length > 1"
              :value="groupSource(group)"
              :disabled="fetchBusy || groupCame(group)"
              :aria-label="`Site for ${groupName(group)}`"
              :data-testid="`${role}-fetch-group-${index}-source`"
              @change="onGroupSource(group, $event)"
            >
              <option v-for="source in group.sources" :key="source" :value="source">{{ sourceLabel(source) }}</option>
            </select>
            <span v-else aria-hidden="true">{{ sourceLabel(group.sources[0]!) }}</span>
          </li>
        </ul>
        <p v-if="fetchState?.state === 'racing'" class="muted" role="status" :data-testid="`${role}-fetch-status`">
          <LosatLoader />Finding the fastest site…
        </p>
        <p
          v-else-if="fetchState?.state === 'fetching'"
          class="muted"
          role="status"
          :data-testid="`${role}-fetch-status`"
        >
          <LosatLoader />Fetching {{ fetchState.done + 1 }} of {{ fetchState.total }} from {{ sourceLabel(fetchState.source) }}…
        </p>
        <p v-else-if="fetchState?.state === 'failed'" class="error" role="alert" :data-testid="`${role}-fetch-error`">
          {{ fetchState.message }}
        </p>
      </div>
      <div class="entry-subrange">
        <RegionPicker
          v-if="regionRecord"
          :draft="draft"
          :role="role"
          :title="`${title} subrange`"
          :record="regionRecord"
          :region="region"
          :unit="residueUnit(kind)"
        />
        <!-- With no single record, no subrange: nothing says so (stop-describing-ui rule 7). -->
      </div>
    </div>

    <div class="param-row file-row">
      <span class="param-label">Or, upload file</span>
      <div class="param-control dropzone" :data-testid="`${role}-dropzone`">
        <button ref="openButton" type="button" :data-testid="`${role}-open`" @click="fileInput?.click()">Open FASTA files…</button>
        <!-- Opened by the button: out of the accessibility tree, where it was a second button (WA2-AX-13). -->
        <input
          ref="fileInput"
          class="visually-hidden"
          type="file"
          multiple
          tabindex="-1"
          aria-hidden="true"
          :aria-label="`Open ${title.toLowerCase()} FASTA files`"
          :data-testid="`${role}-files`"
          @change="onFiles"
        />
        <span class="hint">or drop files here</span>
      </div>
    </div>

    <ParameterField v-for="field in codeFields" :key="field.flag" :model="formModel" :field="field" />
    <p v-if="subjectCodeNote" class="note" data-testid="subject-gencode-note">
      Results with this code can differ from NCBI BLAST+, which translates local subjects otherwise (approved LOSAT exception).
    </p>

    <div v-if="role === 'query'" class="param-row">
      <label class="param-label" for="job-title">Job Title</label>
      <div class="param-control fill">
        <input
          id="job-title"
          type="text"
          autocomplete="off"
          :value="state.title"
          aria-describedby="job-title-help"
          data-testid="job-title"
          @input="draft.setTitle(($event.target as HTMLInputElement).value)"
        />
        <p id="job-title-help" class="help">Enter a descriptive title for your search</p>
      </div>
    </div>

    <ul v-if="roleState.sources.length > 0" ref="sourceItems" class="sources">
      <SourceCard
        v-for="source in roleState.sources"
        :key="source.key"
        :draft="draft"
        :state="state"
        :role="role"
        :source="source"
        :index="roleState.sources.indexOf(source)"
        @remove="removeSource(source.key)"
        @split-accessions="splitAccessions"
      />
    </ul>

    <fieldset v-if="readySources.length > 1" class="mode" :data-testid="`${role}-mode`">
      <legend>Several {{ title.toLowerCase() }} inputs</legend>
      <label>
        <input
          type="radio"
          :name="`${role}-mode`"
          value="combined"
          :checked="roleState.mode === 'combined'"
          :data-testid="`${role}-mode-combined`"
          @change="draft.setMode(role, 'combined')"
        />
        Combined: one search of all
      </label>
      <label>
        <input
          type="radio"
          :name="`${role}-mode`"
          value="separate"
          :checked="roleState.mode === 'separate'"
          :data-testid="`${role}-mode-separate`"
          @change="draft.setMode(role, 'separate')"
        />
        Separate: one search for each
      </label>
    </fieldset>
    </template>
  </fieldset>
</template>
