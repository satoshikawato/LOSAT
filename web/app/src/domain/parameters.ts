// The search form's values and the argv parameters they make (plan §5.3). The form keeps
// only what the user set; a value equal to the engine's default is not written, so the
// argv and the CLI command stay minimal (DW-4). Nothing here checks BLAST's rules: the
// engine's `validate` reads the argv as the command line does.
import type { ParameterField, ParameterSection, ProgramDescriptor } from './programs';

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
 * The parameters of the form, in the order of its fields: set flags, and values other than
 * '' and the engine's default (white space around a value is dropped).
 */
export function formParameters(
  program: ProgramDescriptor,
  values: FormValues,
  options: readonly EngineOption[] | undefined,
): Array<readonly [string, string | true]> {
  const byFlag = new Map((options ?? []).map((option) => [option.flag, option]));
  const parameters: Array<readonly [string, string | true]> = [];
  for (const section of describedSections(program, options)) {
    for (const field of section.fields) {
      const value = values[field.flag];
      if (field.kind === 'flag') {
        if (value === true) parameters.push([field.flag, true]);
        continue;
      }
      if (typeof value !== 'string') continue;
      const text = value.trim();
      if (text === '' || text === byFlag.get(field.flag)?.defaultValue) continue;
      parameters.push([field.flag, text]);
    }
  }
  return parameters;
}
