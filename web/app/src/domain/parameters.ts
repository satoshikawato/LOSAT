// The search form's values and the argv parameters they make (plan §5.3). The form keeps
// only what the user set; a value equal to the default (the engine's, or the program's
// resolved default, `knownDefault`) is not written, so the argv and the CLI command stay
// minimal (DW-4). Nothing here checks BLAST's rules: the engine's `validate` reads the argv
// as the command line does.
import type { ParameterField, ParameterSection, ProgramDescriptor, SectionPlacement } from './programs';

/** What `describe` says about one option (the fields of ports/engine ParameterDescription). */
export interface EngineOption {
  readonly flag: string;
  readonly takesValue: boolean;
  readonly defaultValue?: string;
  readonly choices?: readonly string[];
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
 * The default of a field: the program's resolved default (`ProgramDescriptor.defaults`, by
 * task, while the options it depends on keep their defaults), or else the engine's
 * (`describe`); undefined when neither is known.
 */
export function knownDefault(
  program: ProgramDescriptor,
  flag: string,
  values: FormValues,
  options: readonly EngineOption[] | undefined,
): KnownDefault | undefined {
  const resolved = program.defaults?.[flag];
  if (resolved === undefined) {
    const described = options?.find((option) => option.flag === flag)?.defaultValue;
    return described === undefined ? undefined : { defaultValue: described };
  }
  const task = effectiveValue('-task', values, options);
  const value = typeof resolved.value === 'string' ? resolved.value : task === undefined ? undefined : resolved.value[task];
  const holds = (resolved.unless ?? []).every((other) => !isSet(program, other, values, options));
  const held = holds ? value : undefined;
  if (resolved.given === true) return held === undefined ? { given: true } : { defaultValue: held, given: true };
  return held === undefined ? undefined : { defaultValue: held };
}

/** Whether a field's value is written to the argv (`writtenValue`), found by its flag. */
function isSet(program: ProgramDescriptor, flag: string, values: FormValues, options: readonly EngineOption[] | undefined): boolean {
  const field = program.sections.flatMap((section) => section.fields).find((f) => f.flag === flag);
  if (field === undefined) return false;
  // A given option is set by any value: its default is not read, so that the hit list options,
  // which depend on each other, do not ask each other.
  const known = program.defaults?.[flag]?.given === true ? { given: true } : knownDefault(program, flag, values, options);
  return writtenValue(field, values[flag], known) !== undefined;
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
