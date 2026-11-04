// The search form's values and the argv parameters they make (plan §5.3). The form keeps
// only what the user set; a value equal to the default (the one that the engine resolves,
// `knownDefault`) is not written, so the argv and the CLI command stay minimal (DW-4).
// Nothing here checks BLAST's rules: the engine's `validate` reads the argv as the command
// line does, and the engine's `describe` gives the defaults it resolves.
import type { PairSelect, ParameterField, ParameterSection, ProgramDescriptor, SectionPlacement } from './programs';

/**
 * The default that the engine resolves for an option that a search leaves out (`describe`'s
 * `defaults`, docs/web/abi_v2.md §9), by task and from other options.
 */
export interface ResolvedDefault {
  /** The value of each task with no other option given (the key '' for a program without -task). */
  readonly byTask: Readonly<Record<string, string>>;
  /**
   * For each option that changes this one (-matrix, -word_size): the value by task for each
   * value of that option that the engine accepts.
   */
  readonly by?: Readonly<Record<string, Readonly<Record<string, Readonly<Record<string, string>>>>>>;
  /** The value holds only while none of these options is given (the hit list options). */
  readonly unless?: readonly string[];
  /**
   * The tasks where BLAST+ tells the option given from not given: a value set is written even
   * when it equals the default.
   */
  readonly given?: readonly string[];
}

/** What `describe` says about one option (the fields of ports/engine ParameterDescription). */
export interface EngineOption {
  readonly flag: string;
  readonly takesValue: boolean;
  readonly defaultValue?: string;
  readonly choices?: readonly string[];
  readonly resolvedDefault?: ResolvedDefault;
}

/** A field's value: the text of a value option ('' when not set), or whether a flag is set. */
export type FieldValue = string | boolean;
export type FormValues = Readonly<Record<string, FieldValue>>;

/**
 * Field values that one choice clears (S12 instructions, item 8): BLASTN's templates
 * belong to the discontiguous lookup, which NCBI refuses for the blastn and blastn-short
 * tasks, so choosing those tasks empties the template fields.
 */
const CLEARED_BY: ReadonlyArray<{ readonly flag: string; readonly values: readonly string[]; readonly clears: readonly string[] }> = [
  { flag: '-task', values: ['blastn', 'blastn-short'], clears: ['-template_type', '-template_length'] },
];

/** The sections and fields that the engine describes (all of them while `options` is unknown). */
export function describedSections(
  program: ProgramDescriptor,
  options: readonly EngineOption[] | undefined,
): readonly ParameterSection[] {
  if (options === undefined) return program.sections;
  const described = new Set(options.map((option) => option.flag));
  return program.sections
    .map((section) => ({ ...section, fields: section.fields.filter((field) => described.has(field.flag)) }))
    .filter((section) => section.fields.length > 0);
}

/**
 * The form's control of an option among the fields that the engine describes: the field's label,
 * or its row's name for NCBI's select of a pair ("Gap Costs"); undefined for an option that the
 * form does not have. "Additional options" refuses such an option (DW-42): the form leaves out a
 * value equal to the default that it reads from its own values, which a -task or -matrix written
 * beside the form would change.
 */
export function formControlOf(program: ProgramDescriptor, options: readonly EngineOption[] | undefined, flag: string): string | undefined {
  for (const section of describedSections(program, options)) {
    const field = section.fields.find((f) => f.flag === flag);
    if (field === undefined) continue;
    return field.row !== undefined && section.pairs?.[field.row] !== undefined ? field.row : field.label;
  }
  return undefined;
}

/** The described sections shown at one place of the screen, in the descriptor's order. */
export function sectionsAt(
  program: ProgramDescriptor,
  options: readonly EngineOption[] | undefined,
  placement: SectionPlacement,
): readonly ParameterSection[] {
  return describedSections(program, options).filter((section) => section.placement === placement);
}

/** The flags of every field at one place of the screen, described or not ("Restore default search parameters"). */
export function placementFlags(program: ProgramDescriptor, placement: SectionPlacement): readonly string[] {
  return program.sections.filter((section) => section.placement === placement).flatMap((section) => section.fields.map((field) => field.flag));
}

/**
 * The value that the search uses for an option: the value set (without the white space
 * around it), or else the engine's default; undefined when there is neither.
 */
export function effectiveValue(flag: string, values: FormValues, options: readonly EngineOption[] | undefined): string | undefined {
  const value = values[flag];
  if (typeof value === 'string' && value.trim() !== '') return value.trim();
  return options?.find((option) => option.flag === flag)?.defaultValue;
}

/** What the form knows of an option's default (`knownDefault`). */
export interface KnownDefault {
  /** The value that the search uses while the field is empty; undefined when unknown. */
  readonly defaultValue?: string;
  /** BLAST+ tells the option given from not given: a value set is written even when it equals the default. */
  readonly given?: boolean;
}

/**
 * The default of a field: the value that the engine resolves for the task (`describe`), the
 * value that a set -matrix or -word_size gives it, or else the clap default, or else the value
 * of NCBI's page of an option read with others (a section's `together`, which the argv writes
 * when another of them is written); undefined when unknown (two of the options that change it
 * are set, one of them at a value that the engine's table does not hold, or an option that
 * excludes it is set).
 */
export function knownDefault(
  program: ProgramDescriptor,
  flag: string,
  values: FormValues,
  options: readonly EngineOption[] | undefined,
): KnownDefault | undefined {
  const option = options?.find((o) => o.flag === flag);
  const resolved = option?.resolvedDefault;
  if (option === undefined || resolved === undefined) {
    const described = option?.defaultValue ?? togetherValue(program, flag);
    return described === undefined ? undefined : { defaultValue: described };
  }
  const task = effectiveValue('-task', values, options) ?? '';
  let value = taskDefault(option, task);
  const deciders = Object.entries(resolved.by ?? {})
    .map(([other, table]) => [setValue(program, other, values, options), table] as const)
    .filter(([set]) => set !== undefined);
  if (deciders.length === 1) {
    const [set, table] = deciders[0]!;
    value = typeof set === 'string' ? table[set]?.[task] : undefined;
  } else if (deciders.length > 1) value = undefined;
  if ((resolved.unless ?? []).some((other) => setValue(program, other, values, options) !== undefined)) value = undefined;
  if (resolved.given?.includes(task) === true) return value === undefined ? { given: true } : { defaultValue: value, given: true };
  return value === undefined ? undefined : { defaultValue: value };
}

/** The value of NCBI's page of an option that a section's `together` holds. */
function togetherValue(program: ProgramDescriptor, flag: string): string | undefined {
  for (const section of program.sections) {
    if (section.together !== undefined && Object.hasOwn(section.together, flag)) return section.together[flag];
  }
  return undefined;
}

/** An option's default for a task alone: the engine's resolved value, or else the clap default. */
function taskDefault(option: EngineOption, task: string): string | undefined {
  return option.resolvedDefault?.byTask[task] ?? option.defaultValue;
}

/**
 * The value written to the argv for a field found by its flag (`writtenValue`), against its
 * task's default: a given option is set by any value, so that the options that depend on each
 * other (the hit list options, the word size and the threshold) do not ask each other.
 */
function setValue(program: ProgramDescriptor, flag: string, values: FormValues, options: readonly EngineOption[] | undefined): string | true | undefined {
  const field = program.sections.flatMap((section) => section.fields).find((f) => f.flag === flag);
  if (field === undefined) return undefined;
  const option = options?.find((o) => o.flag === flag);
  const task = effectiveValue('-task', values, options) ?? '';
  if (option?.resolvedDefault?.given?.includes(task) === true) return writtenValue(field, values[flag], { given: true });
  const known = option === undefined ? undefined : taskDefault(option, task);
  return writtenValue(field, values[flag], known === undefined ? undefined : { defaultValue: known });
}

/**
 * Whether a section is shown (`shownFor`): always, unless it is shown only for some values of
 * an option, and the option has another value and none of the section's fields has a value.
 */
export function sectionShown(section: ParameterSection, values: FormValues, options: readonly EngineOption[] | undefined): boolean {
  const rule = section.shownFor;
  if (rule === undefined) return true;
  const current = effectiveValue(rule.flag, values, options);
  if (current !== undefined && rule.values.includes(current)) return true;
  return section.fields.some((field) => {
    const value = values[field.flag];
    return value === true || (typeof value === 'string' && value.trim() !== '');
  });
}

/**
 * The fields of a section in rows: a field alone, or the consecutive fields of one named row
 * (NCBI's "Match/Mismatch Scores", "Gap Costs").
 */
export function fieldRows(
  fields: readonly ParameterField[],
): ReadonlyArray<{ readonly row?: string; readonly fields: readonly ParameterField[] }> {
  const rows: Array<{ row?: string; fields: ParameterField[] }> = [];
  for (const field of fields) {
    const last = rows.at(-1);
    if (field.row !== undefined && last?.row === field.row) last.fields.push(field);
    else rows.push(field.row === undefined ? { fields: [field] } : { row: field.row, fields: [field] });
  }
  return rows;
}

/**
 * What the form writes to the argv for a field: `true` for a set flag, or a value other than
 * '' and the default (white space around it dropped; a `given` option's value even when it
 * equals the default); undefined when it writes nothing. A written value is a value that
 * differs from the default ("♦"), since a value equal to the default is not written.
 */
export function writtenValue(field: ParameterField, value: FieldValue | undefined, known: KnownDefault | undefined): string | true | undefined {
  if (field.kind === 'flag') return value === true ? true : undefined;
  if (typeof value !== 'string') return undefined;
  const text = value.trim();
  return text === '' || (known?.given !== true && text === known?.defaultValue) ? undefined : text;
}

/** The values of a field's select (`list`), for the task the search uses; empty without one. */
export function fieldList(field: ParameterField, values: FormValues, options: readonly EngineOption[] | undefined): readonly string[] {
  const list = field.list;
  if (list === undefined) return [];
  if (Array.isArray(list)) return list as readonly string[];
  const task = effectiveValue('-task', values, options);
  return (task === undefined ? undefined : (list as Readonly<Record<string, readonly string[]>>)[task]) ?? [];
}

/** The choices of a `choice` field: the engine's, or the field's own where the engine lists none. */
export function fieldChoices(field: ParameterField, option: EngineOption | undefined): readonly string[] {
  return option?.choices !== undefined && option.takesValue ? option.choices : (field.choices ?? []);
}

/** A row of two fields with NCBI's select of their pairs (a section's `pairs`). */
export interface PairRow {
  readonly row: string;
  readonly select: PairSelect;
  readonly flags: readonly [string, string];
}

const PAIR_ROWS = new WeakMap<ProgramDescriptor, readonly PairRow[]>();

/**
 * The rows of a program that have NCBI's select of their two fields, the same objects for each
 * call (a component's state that follows its row must not follow every render).
 */
export function pairRows(program: ProgramDescriptor): readonly PairRow[] {
  let rows = PAIR_ROWS.get(program);
  if (rows === undefined) {
    rows = program.sections.flatMap((section) =>
      Object.entries(section.pairs ?? {}).flatMap(([row, select]) => {
        const [first, second, ...more] = section.fields.filter((field) => field.row === row);
        return first === undefined || second === undefined || more.length > 0 ? [] : [{ row, select, flags: [first.flag, second.flag] as const }];
      }),
    );
    PAIR_ROWS.set(program, rows);
  }
  return rows;
}

/** The value that the search uses for a field: the value set (trimmed), or else its known default. */
function usedValue(program: ProgramDescriptor, flag: string, values: FormValues, options: readonly EngineOption[] | undefined): string | undefined {
  const value = values[flag];
  if (typeof value === 'string' && value.trim() !== '') return value.trim();
  return knownDefault(program, flag, values, options)?.defaultValue;
}

/** The key of a pair row's list: the values that the search uses for its `by` options, joined with ','. */
function listKey(program: ProgramDescriptor, select: PairSelect, values: FormValues, options: readonly EngineOption[] | undefined): string | undefined {
  if (select.by === undefined) return '';
  const parts = select.by.map((flag) => usedValue(program, flag, values, options));
  return parts.some((part) => part === undefined) ? undefined : parts.join(',');
}

/**
 * The values of a pair row's select: NCBI's list for the scores or the matrix that the search
 * uses, without the values of other tasks; empty without a list.
 */
export function pairList(program: ProgramDescriptor, pair: PairRow, values: FormValues, options: readonly EngineOption[] | undefined): readonly string[] {
  const key = listKey(program, pair.select, values, options);
  // A key typed after "Other" (a matrix named `constructor`) is not one of the lists'.
  const list = key === undefined || !Object.hasOwn(pair.select.lists, key) ? undefined : pair.select.lists[key];
  if (list === undefined) return [];
  const only = pair.select.only;
  if (only === undefined) return list;
  const task = effectiveValue('-task', values, options) ?? '';
  return list.filter((value) => !only.values.includes(value) || only.tasks.includes(task));
}

/** NCBI's default of a pair row's list: a value of the task's own (megablast's Linear), or else NCBI's "Deflt". */
function pairListDefault(program: ProgramDescriptor, pair: PairRow, values: FormValues, options: readonly EngineOption[] | undefined): string | undefined {
  const own = pairList(program, pair, values, options).find((value) => pair.select.only?.values.includes(value));
  if (own !== undefined) return own;
  const key = listKey(program, pair.select, values, options);
  return key === undefined || pair.select.defaults === undefined || !Object.hasOwn(pair.select.defaults, key) ? undefined : pair.select.defaults[key];
}

/** The pair that the search uses ("1,-2", "11 1"), or '' while one of its values is unknown. */
export function pairValue(program: ProgramDescriptor, pair: PairRow, values: FormValues, options: readonly EngineOption[] | undefined): string {
  const parts = pair.flags.map((flag) => usedValue(program, flag, values, options));
  return parts.some((part) => part === undefined) ? '' : parts.join(pair.select.separator);
}

/**
 * The field values that choosing a pair writes: its two values, so that the pair stays the one
 * chosen when another task has other defaults (megablast's 1,-4 does not become blastn's 2,-4),
 * or '' for both where the pair is the default (the fields follow the task again, as a default
 * chosen in a select does).
 */
export function pairEntries(
  program: ProgramDescriptor,
  pair: PairRow,
  value: string,
  values: FormValues,
  options: readonly EngineOption[] | undefined,
): Array<readonly [string, FieldValue]> {
  const parts = value.split(pair.select.separator);
  const isDefault = pair.flags.every((flag, i) => {
    const known = knownDefault(program, flag, values, options);
    return known?.given !== true && parts[i] === known?.defaultValue;
  });
  return pair.flags.map((flag, i) => [flag, isDefault ? '' : (parts[i] ?? '')] as const);
}

/**
 * Sets fields and clears the fields that a new value excludes. A pair row whose list the new
 * values change (the task, the scores or the matrix) takes the new list's default, as on NCBI's
 * page, where the gap costs follow the scores and the matrix: a pair that the list does not hold
 * is one that BLAST may not read (4,-5 with BLASTN's 5 2 stops the search). So does a row whose
 * pair in use leaves its list while the list stays (scores 1,-1 with megablast's default 0 0).
 * A pair row whose own fields are set here keeps them.
 */
export function setFields(
  program: ProgramDescriptor,
  values: FormValues,
  entries: ReadonlyArray<readonly [string, FieldValue]>,
  options: readonly EngineOption[] | undefined,
): FormValues {
  const next: Record<string, FieldValue> = { ...values };
  for (const [flag, value] of entries) {
    next[flag] = value;
    for (const rule of CLEARED_BY) {
      if (rule.flag === flag && typeof value === 'string' && rule.values.includes(value)) {
        for (const cleared of rule.clears) delete next[cleared];
      }
    }
  }
  for (const pair of pairRows(program)) {
    if (entries.some(([flag]) => pair.flags.includes(flag))) continue;
    const list = pairList(program, pair, next, options);
    const used = pairValue(program, pair, next, options);
    const listChanged = pairList(program, pair, values, options).join('|') !== list.join('|');
    const leftList = used !== pairValue(program, pair, values, options) && list.length > 0 && !list.includes(used);
    if (!listChanged && !leftList) continue;
    const chosen = pairListDefault(program, pair, next, options);
    if (chosen === undefined) continue;
    for (const [flag, value] of pairEntries(program, pair, chosen, next, options)) {
      if (value === '') delete next[flag];
      else next[flag] = value;
    }
  }
  return Object.freeze(next);
}

/** Sets one field (`setFields`). */
export function setField(
  program: ProgramDescriptor,
  values: FormValues,
  flag: string,
  value: FieldValue,
  options: readonly EngineOption[] | undefined,
): FormValues {
  return setFields(program, values, [[flag, value]], options);
}

/**
 * The parameters of the form, in the order of its fields: the values that `writtenValue`
 * writes, and with one of a section's `together` options the others at the values that the
 * search uses (WA-FM-03: BLAST+ refuses a template length without its type). A placement limits
 * them to the fields shown at that place.
 */
export function formParameters(
  program: ProgramDescriptor,
  values: FormValues,
  options: readonly EngineOption[] | undefined,
  placement?: SectionPlacement,
): Array<readonly [string, string | true]> {
  const parameters: Array<readonly [string, string | true]> = [];
  for (const section of describedSections(program, options)) {
    if (placement !== undefined && section.placement !== placement) continue;
    const written = section.fields.map((field) => [field.flag, writtenValue(field, values[field.flag], knownDefault(program, field.flag, values, options))] as const);
    const together = section.together ?? {};
    const partners = written.some(([flag, value]) => value !== undefined && Object.hasOwn(together, flag));
    for (const [flag, value] of written) {
      const used = value ?? (partners && Object.hasOwn(together, flag) ? usedValue(program, flag, values, options) : undefined);
      if (used !== undefined) parameters.push([flag, used]);
    }
  }
  return parameters;
}

/**
 * The fields of the form whose values differ from the defaults ("♦"), at one place of the
 * screen ("(n changed)"): the parameters of the argv, without the `together` options written
 * at their defaults.
 */
export function changedFields(program: ProgramDescriptor, values: FormValues, options: readonly EngineOption[] | undefined, placement: SectionPlacement): number {
  return describedSections(program, options)
    .filter((section) => section.placement === placement)
    .flatMap((section) => section.fields)
    .filter((field) => writtenValue(field, values[field.flag], knownDefault(program, field.flag, values, options)) !== undefined).length;
}
