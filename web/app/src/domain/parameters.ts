// The search form's values and the argv parameters they make (plan §5.3). The form keeps
// only what the user set; a value equal to the default (the one that the engine resolves,
// `knownDefault`) is not written, so the argv and the CLI command stay minimal (DW-4).
// Nothing here checks BLAST's rules: the engine's `validate` reads the argv as the command
// line does, and the engine's `describe` gives the defaults it resolves.
import type { ParameterField, ParameterSection, ProgramDescriptor, SectionPlacement } from './programs';

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
 * value that a set -matrix or -word_size gives it, or else the clap default; undefined when
 * unknown (two of the options that change it are set, one of them at a value that the
 * engine's table does not hold, or an option that excludes it is set).
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
    const described = option?.defaultValue;
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

/** Sets one field and clears the fields that its new value excludes. */
export function setField(values: FormValues, flag: string, value: FieldValue): FormValues {
  const next: Record<string, FieldValue> = { ...values, [flag]: value };
  for (const rule of CLEARED_BY) {
    if (rule.flag === flag && typeof value === 'string' && rule.values.includes(value)) {
      for (const cleared of rule.clears) delete next[cleared];
    }
  }
  return Object.freeze(next);
}

/**
 * The parameters of the form, in the order of its fields: the values that `writtenValue`
 * writes. A placement limits them to the fields shown at that place.
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
    for (const field of section.fields) {
      const written = writtenValue(field, values[field.flag], knownDefault(program, field.flag, values, options));
      if (written !== undefined) parameters.push([field.flag, written]);
    }
  }
  return parameters;
}
