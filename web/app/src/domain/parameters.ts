// The search form's values and the argv parameters they make (plan §5.3). The form keeps
// only what the user set; a value equal to the engine's default is not written, so the
// argv and the CLI command stay minimal (DW-4). Nothing here checks BLAST's rules: the
// engine's `validate` reads the argv as the command line does.
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
 * '' and the engine's default (white space around it dropped); undefined when it writes
 * nothing. A written value is a value that differs from the default ("♦"), since a value
 * equal to the default is not written.
 */
export function writtenValue(field: ParameterField, value: FieldValue | undefined, option: EngineOption | undefined): string | true | undefined {
  if (field.kind === 'flag') return value === true ? true : undefined;
  if (typeof value !== 'string') return undefined;
  const text = value.trim();
  return text === '' || text === option?.defaultValue ? undefined : text;
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
  const byFlag = new Map((options ?? []).map((option) => [option.flag, option]));
  const parameters: Array<readonly [string, string | true]> = [];
  for (const section of describedSections(program, options)) {
    if (placement !== undefined && section.placement !== placement) continue;
    for (const field of section.fields) {
      const written = writtenValue(field, values[field.flag], byFlag.get(field.flag));
      if (written !== undefined) parameters.push([field.flag, written]);
    }
  }
  return parameters;
}
