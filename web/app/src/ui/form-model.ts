// Shared form model for ParameterField, ParameterForm and ParameterPair (WW3b DESIGN §2.3):
// both the Search form (`searchFormModel`) and tool forms provide one, so one field component
// draws every option with its defaults, changed marker, help tip and refusal marks.
import type { ParameterDescription } from 'losat/ports';
import { registryName } from '../application/conffile';
import type { CountsLoader } from '../application/counts';
import type { DatabasePickDraft, DraftState, SearchDraft } from '../application/draft';
import type { ToolDraft, ToolDraftState } from '../application/tool-draft';
import {
  fieldList,
  knownDefault,
  pairEntries,
  pairList,
  pairRows,
  pairValue,
  sectionShown,
  sectionsAt,
  type FieldValue,
  type FormValues,
  type KnownDefault,
  type PairRow,
} from '../domain/parameters';
import { NO_TAXID_EXPANSION_FLAG } from '../domain/blast-db';
import { programById, type ParameterField, type ParameterSection } from '../domain/programs';
import { refusedToolFields } from '../domain/tool-argv';
import { allStageParameters, stageOptionParameters, toToolEngineOption, toToolParameterField } from '../domain/tools';
import { refusal, REFUSAL_ID } from './refusal';

export interface FormRefusal {
  readonly extra: boolean;
  readonly flags: ReadonlySet<string>;
}

export interface FormModel {
  /** Identifies the form's current mode/stage so Vue keys controls separately across switches. */
  readonly key: string;
  readonly values: FormValues;
  readonly options: readonly ParameterDescription[] | undefined;
  readonly files: Readonly<Record<string, File>>;
  knownDefault(flag: string): KnownDefault | undefined;
  fieldList(field: ParameterField): readonly string[];
  readonly refusal: FormRefusal | undefined;
  /** DOM id of the refusal message element that refused controls point to via `aria-describedby`. */
  readonly refusalId?: string;
  /** Prefix for control DOM ids and `data-testid`s (`'param'` for Search). */
  readonly idPrefix?: string;
  readonly descriptionError?: string;
  readonly descriptionUnavailable?: boolean;
  setField(flag: string, value: FieldValue): void;
  setFile(flag: string, file: File | undefined): void;

  // --- Optional capabilities used by Search (and/or tool forms) ---
  readonly sections?: readonly ParameterSection[];
  readonly extra?: string;
  readonly extraHelp?: string;
  setExtra?(text: string): void;
  readonly conffileName?: string;
  readonly conffile?: File;
  setConffile?(file: File | undefined): void;
  readonly queryGencodes?: readonly number[];
  readonly subjectGencodes?: readonly number[];
  taxidFieldOf?(flag: string): ParameterField | undefined;
  setRepeats?(mode: 'off' | 'taxid' | 'file'): void;
  setTaxid?(text: string): void;
  readonly counts?: CountsLoader;
  /** Optional shortcut on BLASTN's repeats row when Counts file is chosen. */
  onMakeCounts?(): void;
  /** True when Database (`-db`) mode has `-no_taxid_expansion` checked (incompatible with BLASTN `-window_masker_taxid`). */
  readonly noTaxidExpansionActive?: boolean;
  /** Clears `-no_taxid_expansion` in Database filters so BLASTN `-window_masker_taxid` can run. */
  clearNoTaxidExpansion?(): void;
  /** Local repeats database state/actions for BLASTN `-filtering_db` (RG-D4 W2). */
  readonly filteringDb?: DatabasePickDraft;
  setFilteringDbEnabled?(enabled: boolean): void;
  setFilteringDbFiles?(files: readonly File[]): void;
  selectFilteringDb?(name: string): void;
  clearFilteringDb?(): void;
  pairOf?(row: string | undefined, fields: readonly ParameterField[]): PairRow | undefined;
  pairList?(pair: PairRow): readonly string[];
  pairValue?(pair: PairRow): string;
  setPair?(pair: PairRow, value: string): void;
}

const SEARCH_EXTRA_HELP = 'Other options of the losat command line, written as there: -qcov_hsp_perc 50 -ungapped';

export interface SearchFormModelOptions {
  readonly counts?: CountsLoader | undefined;
  readonly onMakeCounts?: (() => void) | undefined;
}

function descriptionFields(state: { readonly descriptionError?: string; readonly descriptionUnavailable?: boolean }) {
  return {
    ...(state.descriptionError === undefined ? {} : { descriptionError: state.descriptionError }),
    ...(state.descriptionUnavailable === undefined ? {} : { descriptionUnavailable: state.descriptionUnavailable }),
  };
}

function draftActions(draft: { setField(flag: string, value: FieldValue): void; setExtra(text: string): void }) {
  return {
    setField: (flag: string, value: FieldValue) => draft.setField(flag, value),
    setExtra: (text: string) => draft.setExtra(text),
  };
}

/** Adapts `SearchDraft` + `DraftState` into a `FormModel` for ParameterField / ParameterForm / ParameterPair. */
export function searchFormModel(draft: SearchDraft, state: DraftState, extras: SearchFormModelOptions = {}): FormModel {
  const program = programById(state.program);
  const options = state.description?.parameters;
  const values = state.values[state.program];
  const files = state.files[state.program];
  const sections =
    state.descriptionUnavailable === true
      ? []
      : sectionsAt(program, options, 'algorithm').filter((section) => sectionShown(section, values, options));
  const conffileName = registryName(state.start);
  const conffile = state.conffiles[state.program];
  const currentRefusal = refusal(state);

  return {
    key: state.program,
    values,
    options,
    files,
    knownDefault: (flag) => knownDefault(program, flag, values, options),
    fieldList: (field) => fieldList(field, values, options),
    refusal: currentRefusal,
    refusalId: REFUSAL_ID,
    idPrefix: 'param',
    ...descriptionFields(state),
    ...draftActions(draft),
    setFile: (flag, file) => draft.setFile(flag, file),
    sections,
    extra: state.extra[state.program],
    extraHelp: SEARCH_EXTRA_HELP,
    ...(conffileName === undefined ? {} : { conffileName }),
    ...(conffile === undefined ? {} : { conffile }),
    setConffile: (file) => draft.setConffile(file),
    ...(state.description?.queryGencodes === undefined ? {} : { queryGencodes: state.description.queryGencodes }),
    ...(state.description?.subjectGencodes === undefined ? {} : { subjectGencodes: state.description.subjectGencodes }),
    taxidFieldOf: (flag) =>
      options?.some((option) => option.flag === flag) === false
        ? undefined
        : program.sections.flatMap((section) => section.fields).find((field) => field.flag === flag),
    setRepeats: (mode) => draft.setRepeats(mode),
    setTaxid: (text) => draft.setTaxid(text),
    ...(extras.counts === undefined ? {} : { counts: extras.counts }),
    ...(extras.onMakeCounts === undefined ? {} : { onMakeCounts: extras.onMakeCounts }),
    noTaxidExpansionActive:
      state.program === 'blastn' &&
      state.subjectSource === 'database' &&
      state.database.filterValues[NO_TAXID_EXPANSION_FLAG] === true,
    clearNoTaxidExpansion: () => draft.setDatabaseFilter(NO_TAXID_EXPANSION_FLAG, false),
    filteringDb: state.filteringDb,
    setFilteringDbEnabled: (enabled) => draft.setFilteringDbEnabled(enabled),
    setFilteringDbFiles: (picked) => void draft.setFilteringDbFiles(picked),
    selectFilteringDb: (name) => draft.selectFilteringDb(name),
    clearFilteringDb: () => draft.clearFilteringDb(),
    pairOf: (row, fields) => {
      const pair = pairRows(program).find((entry) => entry.row === row);
      return pair !== undefined && fields.length === 2 && fields.every((field) => pair.flags.includes(field.flag)) ? pair : undefined;
    },
    pairList: (pair) => pairList(program, pair, values, options),
    pairValue: (pair) => pairValue(program, pair, values, options),
    setPair: (pair, value) => draft.setFields(pairEntries(program, pair, value, values, options)),
  };
}

/** DOM id of the WindowMasker argv refusal element, pointed to by refused fields via `aria-describedby`. */
export const WM_REFUSAL_ID = 'wm-argv-validation-message';

const WM_EXTRA_HELP =
  'Other options of the losat windowmasker command line, written as there: -logfile run.log -conffile app.ini -h -version-full';

/** Adapts `ToolDraft` + `ToolDraftState` into a `FormModel` driven by `-xmlhelp`. */
export function windowMaskerFormModel(draft: ToolDraft, state: ToolDraftState): FormModel {
  const stage = state.stage;
  const allParams = allStageParameters(state.description, stage);
  const options: readonly ParameterDescription[] | undefined =
    state.description === undefined ? undefined : allParams.map(toToolEngineOption);
  const optionFields = stageOptionParameters(state.description, stage).map(toToolParameterField);
  const sections: readonly ParameterSection[] =
    state.descriptionUnavailable === true || optionFields.length === 0
      ? []
      : [{ title: 'Options', placement: 'algorithm', fields: optionFields }];
  const values = state.values[stage];
  const files: Record<string, File> = {};
  for (const entry of state.runFiles) {
    if (entry.flag !== undefined) files[entry.flag] = entry.file;
  }
  const currentRefusal =
    state.validation.state === 'invalid' ? refusedToolFields(state.validation.message, state.extra[stage]) : undefined;

  return {
    key: `wm-${stage}`,
    values,
    options,
    files,
    knownDefault: (flag) => {
      const def = options?.find((option) => option.flag === flag)?.defaultValue;
      return def === undefined ? undefined : { defaultValue: def };
    },
    fieldList: () => [],
    refusal: currentRefusal,
    refusalId: WM_REFUSAL_ID,
    idPrefix: 'wm-param',
    ...descriptionFields(state),
    ...draftActions(draft),
    setFile: (flag, file) => draft.setOptionFile(flag, file),
    sections,
    extra: state.extra[stage],
    extraHelp: WM_EXTRA_HELP,
  };
}
