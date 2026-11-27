// Shared form model for ParameterField, ParameterForm and ParameterPair (WW3b DESIGN §2.3):
// both the Search form (`searchFormModel`) and tool forms provide one, so one field component
// draws every option with its defaults, changed marker, help tip and refusal marks.
import type { ToolParameterDescription } from 'losat/domain';
import type { ParameterDescription } from 'losat/ports';
import type { BlastFormatterDraft, BlastFormatterDraftState } from '../application/blast-formatter-draft';
import { registryName } from '../application/conffile';
import type { CountsLoader } from '../application/counts';
import type { DatabasePickDraft, DraftState, SearchDraft } from '../application/draft';
import type { MaskingDraft, MaskingDraftState } from '../application/masking-draft';
import type { ToolDraft, ToolDraftState } from '../application/tool-draft';
import type { ToolFormState } from '../application/tool-draft-core';
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
import {
  allBlastFormatterParameters,
  allMaskingParameters,
  allStageParameters,
  blastFormatterOptionParameters,
  maskingOptionParameters,
  stageOptionParameters,
  toToolEngineOption,
  toToolParameterField,
  WINDOWMASKER_VERSION_FLAGS,
} from '../domain/tools';
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

/** DOM id of the Format Archive argv refusal element, pointed to by refused fields via `aria-describedby`. */
export const BF_REFUSAL_ID = 'bf-argv-validation-message';

const BF_EXTRA_HELP =
  'Other options of the losat blast_formatter command line, written as there: -logfile run.log -conffile app.ini -h -version';

/** DOM id of the Masking argv refusal element, pointed to by refused fields via `aria-describedby`. */
export const MS_REFUSAL_ID = 'ms-argv-validation-message';

/** What a tool tab's form model takes from its draft besides the state that every tool draft holds. */
interface ToolFormParts {
  readonly key: string;
  readonly idPrefix: string;
  readonly refusalId: string;
  readonly extraHelp: string;
  readonly values: FormValues;
  readonly extra: string;
  /** Every option of the form (its blocks' too), for the defaults and the help. */
  readonly all: readonly ToolParameterDescription[];
  /** The options of the form's options section. */
  readonly shown: readonly ToolParameterDescription[];
  readonly setFile: (flag: string, file: File | undefined) => void;
}

/** A tool draft as a `FormModel` driven by `-xmlhelp` (WindowMasker, Masking, Format Archive). */
function toolFormModel(
  draft: { setField(flag: string, value: FieldValue): void; setExtra(text: string): void },
  state: ToolFormState,
  parts: ToolFormParts,
): FormModel {
  const options: readonly ParameterDescription[] | undefined =
    state.description === undefined ? undefined : parts.all.map(toToolEngineOption);
  const fields = parts.shown.map(toToolParameterField);
  const sections: readonly ParameterSection[] =
    state.descriptionUnavailable === true || fields.length === 0 ? [] : [{ title: 'Options', placement: 'algorithm', fields }];
  const files: Record<string, File> = {};
  for (const entry of state.runFiles) {
    if (entry.flag !== undefined) files[entry.flag] = entry.file;
  }
  return {
    key: parts.key,
    values: parts.values,
    options,
    files,
    knownDefault: (flag) => {
      const def = options?.find((option) => option.flag === flag)?.defaultValue;
      return def === undefined ? undefined : { defaultValue: def };
    },
    fieldList: () => [],
    refusal: state.validation.state === 'invalid' ? refusedToolFields(state.validation.message, parts.extra) : undefined,
    refusalId: parts.refusalId,
    idPrefix: parts.idPrefix,
    ...descriptionFields(state),
    ...draftActions(draft),
    setFile: parts.setFile,
    sections,
    extra: parts.extra,
    extraHelp: parts.extraHelp,
  };
}

/** Adapts `ToolDraft` + `ToolDraftState` into a `FormModel` driven by `-xmlhelp`. */
export function windowMaskerFormModel(draft: ToolDraft, state: ToolDraftState): FormModel {
  const stage = state.stage;
  return toolFormModel(draft, state, {
    key: `wm-${stage}`,
    idPrefix: 'wm-param',
    refusalId: WM_REFUSAL_ID,
    extraHelp: WM_EXTRA_HELP,
    values: state.values[stage],
    extra: state.extra[stage],
    all: allStageParameters(state.description, stage),
    shown: stageOptionParameters(state.description, stage),
    setFile: (flag, file) => draft.setOptionFile(flag, file),
  });
}

/** Adapts `BlastFormatterDraft` + `BlastFormatterDraftState` into a `FormModel` driven by `-xmlhelp`. */
export function blastFormatterFormModel(draft: BlastFormatterDraft, state: BlastFormatterDraftState): FormModel {
  return toolFormModel(draft, state, {
    key: 'bf-format',
    idPrefix: 'bf-param',
    refusalId: BF_REFUSAL_ID,
    extraHelp: BF_EXTRA_HELP,
    values: state.values,
    extra: state.extra,
    all: allBlastFormatterParameters(state.description),
    shown: blastFormatterOptionParameters(state.description),
    setFile: () => undefined,
  });
}

/** Adapts `MaskingDraft` + `MaskingDraftState` into a `FormModel` driven by the `-xmlhelp` of the tool shown. */
export function maskingFormModel(draft: MaskingDraft, state: MaskingDraftState): FormModel {
  const { tool, description } = state;
  // The tool's own version word (convert2blastmask has -version, dustmasker and segmasker -version-full).
  const version = description?.parameters.find((param) => WINDOWMASKER_VERSION_FLAGS.has(param.flag))?.flag ?? '-version';
  return toolFormModel(draft, state, {
    key: `ms-${tool}`,
    idPrefix: 'ms-param',
    refusalId: MS_REFUSAL_ID,
    extraHelp: `Other options of the losat ${tool} command line, written as there: -logfile run.log -conffile app.ini -h ${version}`,
    values: state.values[tool],
    extra: state.extra[tool],
    all: allMaskingParameters(description, tool),
    shown: maskingOptionParameters(description, tool),
    setFile: () => undefined,
  });
}
