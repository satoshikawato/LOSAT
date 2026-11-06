import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { duplicateIds } from '../../src/domain/dataset';
import { geneticCodeLabel } from '../../src/domain/genetic-codes';
import {
  changedFields,
  describedSections,
  effectiveValue,
  fieldChoices,
  fieldList,
  fieldRows,
  formParameters,
  knownDefault,
  pairEntries,
  pairList,
  pairRows,
  pairValue,
  placementFlags,
  sectionShown,
  sectionsAt,
  setField,
  setFields,
  writtenValue,
  type EngineOption,
  type FormValues,
} from '../../src/domain/parameters';
import {
  COMPOSITION_MODES,
  fastaParser,
  MATRICES,
  PROGRAMS,
  PROTEIN_GAP_COSTS,
  PROTEIN_WORD_SIZES,
  programById,
  residueUnit,
  searchSummary,
  searchSummaryParts,
  sequenceKind,
  type ParameterField,
  type ProgramId,
} from '../../src/domain/programs';
import {
  REGION_FLAG,
  regionFromPositions,
  regionProblem,
  regionRange,
  regionValue,
  swappedRegion,
} from '../../src/domain/region';
import { estimateKind, looksLikeOtherKind } from '../../src/domain/sequence-kind';
import { toProgramDescription } from '../../src/infra/reactor/control';

const valueOption = (flag: string, defaultValue?: string, choices?: readonly string[]): EngineOption => ({
  flag,
  takesValue: true,
  ...(defaultValue !== undefined ? { defaultValue } : {}),
  ...(choices !== undefined ? { choices } : {}),
});

describe('formParameters', () => {
  const blastn = programById('blastn');
  const options: EngineOption[] = [
    valueOption('-evalue', '10'),
    valueOption('-word_size', '28'),
    valueOption('-task', 'megablast'),
    { flag: '-lcase_masking', takesValue: false },
    { flag: '-subject_besthit', takesValue: false },
  ];

  it('gives no parameters for empty values', () => {
    expect(formParameters(blastn, {}, options)).toEqual([]);
    expect(formParameters(blastn, { '-evalue': '', '-word_size': '   ' }, options)).toEqual([]);
    expect(formParameters(blastn, {}, undefined)).toEqual([]);
  });

  it('omits a value equal to the engine default, also with surrounding white space', () => {
    expect(formParameters(blastn, { '-evalue': '10' }, options)).toEqual([]);
    expect(formParameters(blastn, { '-evalue': ' 10 ' }, options)).toEqual([]);
  });

  it('keeps a different value, trimmed', () => {
    expect(formParameters(blastn, { '-evalue': ' 1e-5 ' }, options)).toEqual([['-evalue', '1e-5']]);
  });

  it('keeps a value when the engine states no default', () => {
    const noDefault: EngineOption[] = [valueOption('-evalue')];
    expect(formParameters(blastn, { '-evalue': '10' }, noDefault)).toEqual([['-evalue', '10']]);
  });

  it('writes a flag field only when its value is true', () => {
    expect(formParameters(blastn, { '-lcase_masking': true }, options)).toEqual([['-lcase_masking', true]]);
    expect(formParameters(blastn, { '-lcase_masking': false }, options)).toEqual([]);
    // A string is not true, even for a flag field.
    expect(formParameters(blastn, { '-lcase_masking': 'true' }, options)).toEqual([]);
  });

  it('ignores a boolean value for a value field', () => {
    expect(formParameters(blastn, { '-evalue': true }, options)).toEqual([]);
  });

  it('skips fields whose flag is not among the engine options', () => {
    const values: FormValues = { '-evalue': '1', '-word_size': '11', '-reward': '2' };
    expect(formParameters(blastn, values, options)).toEqual([
      ['-evalue', '1'],
      ['-word_size', '11'],
    ]);
  });

  it('follows the order of the descriptor, not of the values or the options', () => {
    const values: FormValues = { '-lcase_masking': true, '-word_size': '15', '-evalue': '1', '-task': 'blastn' };
    expect(formParameters(blastn, values, options)).toEqual([
      ['-task', 'blastn'],
      ['-evalue', '1'],
      ['-word_size', '15'],
      ['-lcase_masking', true],
    ]);
  });

  it('considers all descriptor fields when the options are undefined', () => {
    const values: FormValues = { '-reward': '2', '-perc_identity': '90', '-subject_besthit': true, '-evalue': '10' };
    expect(formParameters(blastn, values, undefined)).toEqual([
      ['-evalue', '10'],
      ['-reward', '2'],
      ['-perc_identity', '90'],
      ['-subject_besthit', true],
    ]);
  });

  it('works for TBLASTX with its sections in the order of the page: the genetic codes, then the algorithm parameters', () => {
    const tblastx = programById('tblastx');
    const values: FormValues = {
      '-db_gencode': '11',
      '-query_gencode': '4',
      '-seg': 'yes',
      '-soft_masking': 'true',
      '-matrix': 'PAM30',
      '-culling_limit': '2',
      '-evalue': '1',
    };
    expect(formParameters(tblastx, values, undefined)).toEqual([
      ['-query_gencode', '4'],
      ['-db_gencode', '11'],
      ['-evalue', '1'],
      ['-culling_limit', '2'],
      ['-matrix', 'PAM30'],
      ['-seg', 'yes'],
      ['-soft_masking', 'true'],
    ]);
  });

  it('gives the parameters of one place of the screen', () => {
    const tblastx = programById('tblastx');
    const values: FormValues = { '-db_gencode': '11', '-query_gencode': '4', '-seg': 'yes', '-evalue': '1' };
    expect(formParameters(tblastx, values, undefined, 'algorithm')).toEqual([
      ['-evalue', '1'],
      ['-seg', 'yes'],
    ]);
    expect(formParameters(tblastx, values, undefined, 'subject')).toEqual([['-db_gencode', '11']]);
    expect(formParameters(blastn, { '-task': 'blastn', '-evalue': '1' }, options, 'program')).toEqual([['-task', 'blastn']]);
  });

  it('has NCBI\'s "Max matches in a query range" (-culling_limit) last in General Parameters of every program, which the engine describes', () => {
    const described = JSON.parse(readFileSync(resolve(import.meta.dirname, '../../src/infra/fake/describe.json'), 'utf8')) as Record<string, unknown>;
    for (const program of PROGRAMS) {
      const general = describedSections(program, toProgramDescription(described[program.id]).parameters).find(
        (section) => section.title === 'General Parameters',
      );
      expect(general?.fields.at(-1), program.id).toMatchObject({ flag: '-culling_limit', label: 'Max matches in a query range', kind: 'text' });
      expect(formParameters(program, { '-culling_limit': '2' }, undefined, 'algorithm'), program.id).toEqual([['-culling_limit', '2']]);
    }
  });

  it('has "Number of descriptions" and "Number of alignments" in Other Parameters of the four programs, writing -num_descriptions N and -num_alignments N', () => {
    for (const id of ['blastn', 'blastp', 'tblastn', 'tblastx'] as const) {
      const program = programById(id);
      const other = program.sections.find((section) => section.title === 'Other Parameters');
      const fields = other?.fields.filter((field) => field.flag === '-num_descriptions' || field.flag === '-num_alignments');
      expect(fields?.map((field) => [field.flag, field.label]), id).toEqual([
        ['-num_descriptions', 'Number of descriptions'],
        ['-num_alignments', 'Number of alignments'],
      ]);
      expect(other?.placement, id).toBe('algorithm');
      expect(
        formParameters(program, { '-num_descriptions': '7', '-num_alignments': '3' }, undefined, 'algorithm'),
        id,
      ).toEqual([
        ['-num_descriptions', '7'],
        ['-num_alignments', '3'],
      ]);
    }
  });

  it('omits a TBLASTX value equal to the engine default', () => {
    const tblastx = programById('tblastx');
    const opts: EngineOption[] = [valueOption('-query_gencode', '1'), valueOption('-db_gencode', '1')];
    expect(formParameters(tblastx, { '-query_gencode': '1', '-db_gencode': '2' }, opts)).toEqual([['-db_gencode', '2']]);
  });
});

describe('setField', () => {
  const templates: FormValues = { '-task': 'dc-megablast', '-template_type': 'coding', '-template_length': '18', '-evalue': '1' };
  const blastn = programById('blastn');
  const setField3 = (values: FormValues, flag: string, value: string | boolean) => setField(blastn, values, flag, value, undefined);

  it.each(['blastn', 'blastn-short'])('removes the template fields when the task is %s', (task) => {
    const next = setField3(templates, '-task', task);
    expect(next).toEqual({ '-task': task, '-evalue': '1' });
    expect('-template_type' in next).toBe(false);
    expect('-template_length' in next).toBe(false);
  });

  it.each(['dc-megablast', 'megablast'])('keeps the template fields when the task is %s', (task) => {
    expect(setField3(templates, '-task', task)).toEqual({ ...templates, '-task': task });
  });

  it('does not clear the templates when another field is set', () => {
    expect(setField3(templates, '-evalue', '5')).toEqual({ ...templates, '-evalue': '5' });
  });

  it('returns a frozen object and leaves the input unchanged', () => {
    const before = { ...templates };
    expect(Object.isFrozen(setField3(templates, '-evalue', '5'))).toBe(true);
    expect(Object.isFrozen(setField3(templates, '-task', 'blastn'))).toBe(true);
    expect(templates).toEqual(before);
  });

  it('sets a flag value', () => {
    expect(setField3({}, '-lcase_masking', true)).toEqual({ '-lcase_masking': true });
  });
});

describe('describedSections', () => {
  const blastp = programById('blastp');

  it('returns the program sections when the options are unknown', () => {
    expect(describedSections(blastp, undefined)).toBe(blastp.sections);
  });

  it('keeps only described fields and drops empty sections', () => {
    const sections = describedSections(blastp, [valueOption('-evalue'), valueOption('-matrix')]);
    expect(sections.map((s) => [s.title, s.fields.map((f) => f.flag)])).toEqual([
      ['General Parameters', ['-evalue']],
      ['Scoring Parameters', ['-matrix']],
    ]);
  });

  it('returns no section when the engine describes none of the fields', () => {
    expect(describedSections(blastp, [])).toEqual([]);
    expect(describedSections(blastp, [valueOption('-unrelated')])).toEqual([]);
  });
});

describe('fieldChoices', () => {
  const field = { flag: '-x', label: 'X', kind: 'choice', choices: ['a', 'b'] } as const;

  it("uses the engine's choices for an option that takes a value", () => {
    expect(fieldChoices(field, valueOption('-x', undefined, ['c', 'd']))).toEqual(['c', 'd']);
  });

  it("uses the field's choices when the engine lists none", () => {
    expect(fieldChoices(field, valueOption('-x'))).toEqual(['a', 'b']);
    expect(fieldChoices(field, undefined)).toEqual(['a', 'b']);
  });

  it('gives no choices when neither lists any', () => {
    expect(fieldChoices({ flag: '-x', label: 'X', kind: 'choice' }, undefined)).toEqual([]);
    expect(fieldChoices({ flag: '-x', label: 'X', kind: 'choice' }, valueOption('-x'))).toEqual([]);
  });

  it('does not use the choices of a flag option', () => {
    const flagOption: EngineOption = { flag: '-x', takesValue: false, choices: ['true', 'false'] };
    expect(fieldChoices(field, flagOption)).toEqual(['a', 'b']);
    expect(fieldChoices({ flag: '-x', label: 'X', kind: 'boolean' }, flagOption)).toEqual([]);
  });
});

describe('writtenValue', () => {
  const text: ParameterField = { flag: '-evalue', label: 'Expect threshold', kind: 'text' };
  const flag: ParameterField = { flag: '-lcase_masking', label: 'Mask lower case letters', kind: 'flag' };

  it('writes a value other than empty and the default, trimmed', () => {
    expect(writtenValue(text, ' 1e-5 ', valueOption('-evalue', '10'))).toBe('1e-5');
    expect(writtenValue(text, ' 10 ', valueOption('-evalue', '10'))).toBeUndefined();
    expect(writtenValue(text, '  ', valueOption('-evalue', '10'))).toBeUndefined();
    expect(writtenValue(text, undefined, valueOption('-evalue', '10'))).toBeUndefined();
  });

  it('writes a value equal to no default, as when the engine states none (a default that depends on the task)', () => {
    expect(writtenValue(text, '10', valueOption('-evalue'))).toBe('10');
    expect(writtenValue(text, '10', undefined)).toBe('10');
  });

  it('writes a flag only when it is set', () => {
    expect(writtenValue(flag, true, undefined)).toBe(true);
    expect(writtenValue(flag, false, undefined)).toBeUndefined();
    expect(writtenValue(flag, 'true', undefined)).toBeUndefined();
    expect(writtenValue(text, true, undefined)).toBeUndefined();
  });
});

describe("knownDefault: the defaults that the engine resolves (describe's defaults)", () => {
  // describe.json is the reactor's describe; its `defaults` come from the engine's resolution
  // (web/adapter/src/defaults.rs, whose tests check each value against the CLI's help).
  const described = JSON.parse(readFileSync(resolve(import.meta.dirname, '../../src/infra/fake/describe.json'), 'utf8')) as Record<string, unknown>;
  const optionsOf = (id: ProgramId) => toProgramDescription(described[id]).parameters;
  const blastn = programById('blastn');
  const blastp = programById('blastp');
  const tblastx = programById('tblastx');
  const bn = optionsOf('blastn');
  const bp = optionsOf('blastp');
  const tx = optionsOf('tblastx');

  it('follows the task the search uses, set or the engine default', () => {
    expect(knownDefault(blastn, '-word_size', {}, bn)).toEqual({ defaultValue: '28' });
    expect(knownDefault(blastn, '-word_size', { '-task': 'blastn' }, bn)).toEqual({ defaultValue: '11' });
    expect(knownDefault(blastn, '-evalue', { '-task': 'blastn-short' }, bn)).toEqual({ defaultValue: '1000' });
    expect(knownDefault(blastn, '-dust', { '-task': 'blastn-short' }, bn)).toEqual({ defaultValue: 'no' });
    // No default is known before `describe` comes.
    expect(knownDefault(blastn, '-word_size', {}, undefined)).toBeUndefined();
  });

  it("takes a set matrix's values from the engine's table", () => {
    expect(knownDefault(blastp, '-gapopen', {}, bp)).toEqual({ defaultValue: '11' });
    expect(knownDefault(blastp, '-gapopen', { '-matrix': 'BLOSUM62' }, bp)).toEqual({ defaultValue: '11' });
    expect(knownDefault(blastp, '-gapopen', { '-matrix': 'PAM30' }, bp)).toEqual({ defaultValue: '9' });
    expect(knownDefault(blastp, '-window_size', { '-matrix': 'PAM30' }, bp)).toEqual({ defaultValue: '15' });
    expect(knownDefault(blastp, '-threshold', { '-matrix': 'PAM30', '-task': 'blastp-fast' }, bp)).toEqual({ defaultValue: '20' });
    // A matrix that the engine does not hold: the values that it would decide are unknown.
    expect(knownDefault(blastp, '-gapopen', { '-matrix': 'NOSUCH' }, bp)).toBeUndefined();
    expect(knownDefault(blastp, '-evalue', { '-matrix': 'PAM30' }, bp)).toEqual({ defaultValue: '10' });
    expect(knownDefault(tblastx, '-window_size', { '-matrix': 'BLOSUM45' }, tx)).toEqual({ defaultValue: '60' });
  });

  it("takes a given word size's threshold from the engine's table, and none for two such options set", () => {
    expect(knownDefault(blastp, '-threshold', { '-task': 'blastp-fast' }, bp)).toEqual({ defaultValue: '20' });
    expect(knownDefault(blastp, '-threshold', { '-task': 'blastp-fast', '-word_size': '5' }, bp)).toEqual({ defaultValue: '19.3' });
    expect(knownDefault(blastp, '-threshold', { '-word_size': '6' }, bp)).toEqual({ defaultValue: '21' });
    expect(knownDefault(blastp, '-threshold', { '-word_size': '6', '-matrix': 'PAM30' }, bp)).toBeUndefined();
    // Told given only for the task where giving it changes the search.
    expect(knownDefault(blastp, '-word_size', {}, bp)).toEqual({ defaultValue: '3' });
    expect(knownDefault(blastp, '-word_size', { '-task': 'blastp-fast' }, bp)).toEqual({ defaultValue: '5', given: true });
  });

  it('marks the hit list options given, and shows no default for one while another is set', () => {
    expect(knownDefault(tblastx, '-max_target_seqs', {}, tx)).toEqual({ defaultValue: '500', given: true });
    expect(knownDefault(tblastx, '-max_target_seqs', { '-num_alignments': '3' }, tx)).toEqual({ given: true });
    expect(knownDefault(tblastx, '-num_alignments', { '-max_target_seqs': '500' }, tx)).toEqual({ given: true });
    expect(knownDefault(tblastx, '-num_descriptions', {}, tx)).toEqual({ defaultValue: '500', given: true });
  });

  it("writes BLASTX's word size even at the fast task's default: a word size given sets its own threshold", () => {
    const blastx = programById('blastx');
    const bx = optionsOf('blastx');
    expect(knownDefault(blastx, '-word_size', { '-task': 'blastx-fast' }, bx)).toEqual({ defaultValue: '5', given: true });
    expect(formParameters(blastx, { '-task': 'blastx-fast', '-word_size': '5' }, bx)).toEqual([
      ['-task', 'blastx-fast'],
      ['-word_size', '5'],
    ]);
    expect(knownDefault(blastx, '-threshold', { '-task': 'blastx-fast' }, bx)).toEqual({ defaultValue: '20' });
    expect(knownDefault(blastx, '-threshold', { '-task': 'blastx-fast', '-word_size': '5' }, bx)).toEqual({ defaultValue: '19.3' });
  });

  it("is the clap default where the engine resolves none", () => {
    expect(knownDefault(blastp, '-comp_based_stats', {}, bp)).toEqual({ defaultValue: '2' });
    expect(knownDefault(tblastx, '-evalue', {}, [valueOption('-evalue', '10')])).toEqual({ defaultValue: '10' });
    expect(knownDefault(tblastx, '-evalue', {}, undefined)).toBeUndefined();
  });

  it('omits a value equal to a resolved default from the argv, and writes a given option at its default', () => {
    expect(formParameters(blastn, { '-task': 'blastn', '-word_size': '11', '-reward': '2', '-gapopen': '5' }, bn)).toEqual([
      ['-task', 'blastn'],
    ]);
    expect(formParameters(blastn, { '-word_size': '11' }, bn)).toEqual([['-word_size', '11']]);
    expect(formParameters(blastp, { '-matrix': 'PAM30', '-gapopen': '9', '-gapextend': '1' }, bp)).toEqual([['-matrix', 'PAM30']]);
    expect(formParameters(blastp, { '-matrix': 'PAM30', '-gapopen': '11', '-gapextend': '1' }, bp)).toEqual([
      ['-matrix', 'PAM30'],
      ['-gapopen', '11'],
    ]);
    expect(formParameters(blastp, { '-task': 'blastp-fast', '-word_size': '5' }, bp)).toEqual([
      ['-task', 'blastp-fast'],
      ['-word_size', '5'],
    ]);
    expect(formParameters(blastp, { '-word_size': '3' }, bp)).toEqual([]);
    expect(formParameters(tblastx, { '-max_target_seqs': '500', '-soft_masking': 'false' }, tx)).toEqual([['-max_target_seqs', '500']]);
  });
});

describe('fieldList', () => {
  const blastn = programById('blastn');
  const wordSize = blastn.sections.flatMap((section) => section.fields).find((field) => field.flag === '-word_size')!;
  const task = [valueOption('-task', 'megablast')];

  it("gives NCBI's list of the task the search uses", () => {
    expect(fieldList(wordSize, {}, task)).toEqual(['16', '20', '24', '28', '32', '48', '64']);
    expect(fieldList(wordSize, { '-task': 'blastn' }, task)).toEqual(['7', '11', '15']);
    expect(fieldList(wordSize, {}, undefined)).toEqual([]);
    expect(fieldList({ flag: '-x', label: 'X', kind: 'text' }, {}, task)).toEqual([]);
  });
});

// WA-FM-02: a select offers NCBI's values that LOSAT's engine reads; "Other" still takes the
// rest, which the engine refuses with its message (the engine's e2e checks both sides).
describe("the selects of the protein modes: NCBI's values that LOSAT reads", () => {
  const described = JSON.parse(readFileSync(resolve(import.meta.dirname, '../../src/infra/fake/describe.json'), 'utf8')) as Record<string, unknown>;
  const optionsOf = (id: ProgramId) => toProgramDescription(described[id]).parameters;
  const field = (id: ProgramId, flag: string) => programById(id).sections.flatMap((section) => section.fields).find((f) => f.flag === flag)!;
  const list = (id: ProgramId, flag: string, values: FormValues = {}) => fieldList(field(id, flag), values, optionsOf(id));

  it('gives BLASTP BLOSUM62, word sizes 3, 5 and 6, and modes 0 and 2, in both tasks', () => {
    for (const values of [{}, { '-task': 'blastp-fast' }]) {
      expect(list('blastp', '-matrix', values)).toEqual(['BLOSUM62']);
      expect(list('blastp', '-word_size', values)).toEqual(['3', '5', '6']);
    }
    expect(fieldChoices(field('blastp', '-comp_based_stats'), optionsOf('blastp').find((o) => o.flag === '-comp_based_stats'))).toEqual(['0', '2']);
    // TBLASTN keeps NCBI's four: its engine reads modes 1 and 3 with -matrix IDENTITY (resetting them to 0).
    expect(fieldChoices(field('tblastn', '-comp_based_stats'), optionsOf('tblastn').find((o) => o.flag === '-comp_based_stats'))).toEqual(COMPOSITION_MODES);
  });

  it("gives TBLASTX NCBI's word sizes 2 and 3", () => {
    expect(list('tblastx', '-word_size')).toEqual(['2', '3']);
    expect(list('tblastx', '-matrix')).toEqual(MATRICES);
  });

  it("gives TBLASTN NCBI's matrices with either lookup table (the compressed one of word sizes 5 to 7 too)", () => {
    expect(list('tblastn', '-matrix')).toEqual(MATRICES);
    expect(list('tblastn', '-matrix', { '-task': 'tblastn-fast' })).toEqual(MATRICES);
    expect(list('tblastn', '-matrix', { '-task': 'tblastn-fast', '-word_size': '3' })).toEqual(MATRICES);
    for (const size of ['4', '5', '6', '7']) expect(list('tblastn', '-matrix', { '-word_size': ` ${size} ` })).toEqual(MATRICES);
    expect(list('tblastn', '-word_size')).toEqual(PROTEIN_WORD_SIZES);
  });
});

// WA-FM-03: BLAST+ reads -template_length and -template_type only together.
describe('the Discontiguous Word Options, read together', () => {
  const described = JSON.parse(readFileSync(resolve(import.meta.dirname, '../../src/infra/fake/describe.json'), 'utf8')) as Record<string, unknown>;
  const blastn = programById('blastn');
  const bn = toProgramDescription(described.blastn).parameters;
  const dc = { '-task': 'dc-megablast' };

  it("show NCBI's values, dc-megablast's template, while they are not set", () => {
    expect(knownDefault(blastn, '-template_length', dc, bn)).toEqual({ defaultValue: '18' });
    expect(knownDefault(blastn, '-template_type', dc, bn)).toEqual({ defaultValue: 'coding' });
    expect(formParameters(blastn, dc, bn)).toEqual([['-task', 'dc-megablast']]);
    expect(formParameters(blastn, { ...dc, '-template_length': '18', '-template_type': 'coding' }, bn)).toEqual([['-task', 'dc-megablast']]);
  });

  it('write the other at the value it shows when one is changed, so that the argv never holds one alone', () => {
    expect(formParameters(blastn, { ...dc, '-template_length': '21' }, bn)).toEqual([
      ['-task', 'dc-megablast'],
      ['-template_length', '21'],
      ['-template_type', 'coding'],
    ]);
    expect(formParameters(blastn, { ...dc, '-template_type': 'optimal', '-template_length': ' 18 ' }, bn)).toEqual([
      ['-task', 'dc-megablast'],
      ['-template_length', '18'],
      ['-template_type', 'optimal'],
    ]);
    // Kept with the megablast task (the default, not written), where the engine refuses its word size.
    expect(formParameters(blastn, { '-template_length': '16' }, bn)).toEqual([
      ['-template_length', '16'],
      ['-template_type', 'coding'],
    ]);
  });

  it('count as changed only where they differ from the values shown', () => {
    expect(changedFields(blastn, { ...dc, '-template_length': '21' }, bn, 'algorithm')).toBe(1);
    expect(changedFields(blastn, { ...dc, '-template_length': '21', '-template_type': 'optimal', '-evalue': '1' }, bn, 'algorithm')).toBe(3);
    expect(changedFields(blastn, dc, bn, 'algorithm')).toBe(0);
  });
});

describe('effectiveValue', () => {
  const options = [valueOption('-task', 'megablast')];

  it('is the value set, trimmed, or else the default', () => {
    expect(effectiveValue('-task', { '-task': ' blastn ' }, options)).toBe('blastn');
    expect(effectiveValue('-task', { '-task': '' }, options)).toBe('megablast');
    expect(effectiveValue('-task', {}, options)).toBe('megablast');
    expect(effectiveValue('-task', {}, undefined)).toBeUndefined();
    expect(effectiveValue('-task', { '-task': 'dc-megablast' }, undefined)).toBe('dc-megablast');
  });
});

describe('sections by place, and the Discontiguous Word Options', () => {
  const blastn = programById('blastn');
  const discontiguous = blastn.sections.find((section) => section.title === 'Discontiguous Word Options')!;
  const options = [valueOption('-task', 'megablast'), valueOption('-template_type'), valueOption('-template_length')];

  it('shows the section for the dc-megablast task, set or the default', () => {
    expect(sectionShown(discontiguous, { '-task': 'dc-megablast' }, options)).toBe(true);
    expect(sectionShown(discontiguous, {}, [valueOption('-task', 'dc-megablast')])).toBe(true);
  });

  it('hides it for the other tasks while its fields are empty', () => {
    expect(sectionShown(discontiguous, {}, options)).toBe(false);
    expect(sectionShown(discontiguous, { '-task': 'blastn', '-template_type': ' ' }, options)).toBe(false);
    expect(sectionShown(discontiguous, {}, undefined)).toBe(false);
  });

  it('shows it while one of its fields has a value, so that no hidden value is written', () => {
    expect(sectionShown(discontiguous, { '-task': 'megablast', '-template_length': '18' }, options)).toBe(true);
    expect(sectionShown(discontiguous, { '-template_type': 'coding' }, options)).toBe(true);
  });

  it('always shows a section without a rule', () => {
    expect(sectionShown(blastn.sections[1]!, {}, options)).toBe(true);
  });

  it('places the task, the genetic codes and the algorithm parameters', () => {
    const tblastx = programById('tblastx');
    expect(sectionsAt(tblastx, undefined, 'query').flatMap((s) => s.fields.map((f) => f.flag))).toEqual(['-query_gencode']);
    expect(sectionsAt(tblastx, undefined, 'subject').flatMap((s) => s.fields.map((f) => f.flag))).toEqual(['-db_gencode']);
    expect(sectionsAt(tblastx, undefined, 'program')).toEqual([]);
    expect(sectionsAt(blastn, [valueOption('-evalue')], 'program')).toEqual([]);
    expect(sectionsAt(blastn, undefined, 'program').flatMap((s) => s.fields.map((f) => f.flag))).toEqual(['-task']);
    expect(placementFlags(programById('tblastn'), 'algorithm')).not.toContain('-db_gencode');
    expect(placementFlags(programById('tblastn'), 'algorithm')).not.toContain('-task');
    expect(placementFlags(programById('tblastn'), 'algorithm')).toContain('-evalue');
    expect(sectionsAt(programById('blastx'), undefined, 'query').flatMap((s) => s.fields.map((f) => f.flag))).toEqual(['-query_gencode']);
    expect(sectionsAt(programById('blastx'), undefined, 'subject')).toEqual([]);
  });
});

describe("NCBI's pair selects: Match/Mismatch Scores and Gap Costs", () => {
  const described = JSON.parse(readFileSync(resolve(import.meta.dirname, '../../src/infra/fake/describe.json'), 'utf8')) as Record<string, unknown>;
  const optionsOf = (id: ProgramId) => toProgramDescription(described[id]).parameters;
  const blastn = programById('blastn');
  const blastp = programById('blastp');
  const bn = optionsOf('blastn');
  const bp = optionsOf('blastp');
  const [scores, gaps] = pairRows(blastn);
  const proteinGaps = pairRows(blastp)[0]!;
  // TBLASTN reads NCBI's gap costs of every matrix; BLASTP's list is 11 1 alone (WA-FM-02).
  const tblastn = programById('tblastn');
  const tn = optionsOf('tblastn');
  const tblastnGaps = pairRows(tblastn)[0]!;

  it('are the rows of two fields that NCBI shows as one select', () => {
    expect(pairRows(blastn).map((pair) => [pair.row, pair.select.id, pair.flags])).toEqual([
      ['Match/Mismatch Scores', 'matchscores', ['-reward', '-penalty']],
      ['Gap Costs', 'gapcosts', ['-gapopen', '-gapextend']],
    ]);
    for (const id of ['blastp', 'blastx', 'tblastn'] as const) {
      expect(pairRows(programById(id)).map((pair) => [pair.row, pair.flags])).toEqual([['Gap Costs', ['-gapopen', '-gapextend']]]);
    }
    expect(pairRows(programById('tblastx'))).toEqual([]);
  });

  it("show the pair that the search uses, from the task's defaults", () => {
    expect(pairValue(blastn, scores!, {}, bn)).toBe('1,-2');
    expect(pairValue(blastn, gaps!, {}, bn)).toBe('0 0');
    expect(pairValue(blastn, scores!, { '-task': 'blastn' }, bn)).toBe('2,-3');
    expect(pairValue(blastn, gaps!, { '-task': 'blastn' }, bn)).toBe('5 2');
    expect(pairValue(blastn, gaps!, { '-task': 'blastn', '-gapopen': ' 4 ' }, bn)).toBe('4 2');
    expect(pairValue(blastp, proteinGaps, { '-matrix': 'PAM30' }, bp)).toBe('9 1');
    expect(pairValue(blastn, gaps!, {}, undefined)).toBe('');
  });

  it("list NCBI's gap costs of the scores, with megablast's Linear for megablast only", () => {
    expect(pairList(blastn, scores!, {}, bn)).toEqual(['1,-2', '1,-3', '1,-4', '2,-3', '4,-5', '1,-1']);
    expect(pairList(blastn, gaps!, {}, bn)).toEqual(['0 0', '5 2', '2 2', '1 2', '0 2', '3 1', '2 1', '1 1']);
    expect(pairList(blastn, gaps!, { '-task': 'dc-megablast' }, bn)).toEqual(['4 4', '2 4', '0 4', '3 3', '6 2', '5 2', '4 2', '2 2']);
    expect(pairList(blastn, gaps!, { '-reward': '1', '-penalty': '-1' }, bn)).not.toContain('0 0');
    // Scores outside NCBI's list have no list of gap costs: the two fields show.
    expect(pairList(blastn, gaps!, { '-reward': '3', '-penalty': '-3' }, bn)).toEqual([]);
    expect(pairList(tblastn, tblastnGaps, { '-matrix': 'BLOSUM45' }, tn)).toContain('14 2');
    expect(pairList(tblastn, tblastnGaps, { '-matrix': 'IDENTITY' }, tn)).toEqual([]);
    expect(pairList(tblastn, tblastnGaps, { '-matrix': 'constructor' }, tn)).toEqual([]);
    expect(setField(tblastn, {}, '-matrix', '__proto__', tn)).toEqual({ '-matrix': '__proto__' });
    // BLASTP's list is NCBI's of BLOSUM62, the one matrix that LOSAT's BLASTP reads.
    expect(pairList(blastp, proteinGaps, {}, bp)).toEqual(PROTEIN_GAP_COSTS.lists.BLOSUM62);
    expect(pairList(blastp, proteinGaps, { '-matrix': 'BLOSUM45' }, bp)).toEqual([]);
    expect(pairList(programById('blastx'), pairRows(programById('blastx'))[0]!, {}, optionsOf('blastx'))).toEqual(['11 1']);
  });

  it("give each matrix's list NCBI's default, which is the engine's default for that matrix", () => {
    for (const id of ['tblastn', 'blastp'] as const) {
      const program = programById(id);
      const pair = pairRows(program)[0]!;
      for (const [matrix, value] of Object.entries(pair.select.defaults ?? {})) {
        expect(pairValue(program, pair, { '-matrix': matrix }, optionsOf(id)), `${id} ${matrix}`).toBe(value);
        expect(pair.select.lists[matrix], `${id} ${matrix}`).toContain(value);
      }
    }
  });

  it('write both values of a pair chosen, and neither for the default pair', () => {
    expect(pairEntries(blastn, scores!, '1,-2', {}, bn)).toEqual([
      ['-reward', ''],
      ['-penalty', ''],
    ]);
    // Both, so that the pair stays 1,-3 under a task whose default reward is 2 (the argv still omits -reward 1 here).
    expect(pairEntries(blastn, scores!, '1,-3', {}, bn)).toEqual([
      ['-reward', '1'],
      ['-penalty', '-3'],
    ]);
    expect(formParameters(blastn, setFields(blastn, {}, pairEntries(blastn, scores!, '1,-3', {}, bn), bn), bn)).toEqual([['-penalty', '-3']]);
    expect(pairEntries(blastp, proteinGaps, '10 1', {}, bp)).toEqual([
      ['-gapopen', '10'],
      ['-gapextend', '1'],
    ]);
  });

  it('keep the pair chosen when the task changes the defaults (review of SURVEY #13, 1)', () => {
    const oneFour = setFields(blastn, {}, pairEntries(blastn, scores!, '1,-4', {}, bn), bn);
    const blastnTask = setField(blastn, oneFour, '-task', 'blastn', bn);
    expect(pairValue(blastn, scores!, blastnTask, bn)).toBe('1,-4');
    expect(formParameters(blastn, blastnTask, bn)).toEqual([
      ['-task', 'blastn'],
      ['-reward', '1'],
      ['-penalty', '-4'],
    ]);
  });

  it("give the gap costs the list's default when the pair in use leaves an unchanged list (review of SURVEY #13, 2)", () => {
    // 1,-1 has the same list for every task; megablast's default 0 0 is not in it (BLAST stops).
    const oneOne = setFields(blastn, { '-task': 'blastn' }, pairEntries(blastn, scores!, '1,-1', { '-task': 'blastn' }, bn), bn);
    expect(formParameters(blastn, oneOne, bn)).toEqual([
      ['-task', 'blastn'],
      ['-reward', '1'],
      ['-penalty', '-1'],
    ]);
    const megablast = setField(blastn, oneOne, '-task', 'megablast', bn);
    expect(pairValue(blastn, gaps!, megablast, bn)).toBe('5 2');
    expect(formParameters(blastn, megablast, bn).slice(-2)).toEqual([
      ['-gapopen', '5'],
      ['-gapextend', '2'],
    ]);
  });

  it("give the gap costs the new list's default when the scores change it, as NCBI's page does", () => {
    const choose = (values: FormValues, pair: string) => setFields(blastn, values, pairEntries(blastn, scores!, pair, values, bn), bn);
    // 4,-5 does not take BLASTN's 5 2 (BLAST stops: "Gap existence and extension values 5 and 2 are not supported").
    expect(formParameters(blastn, choose({ '-task': 'blastn' }, '4,-5'), bn)).toEqual([
      ['-task', 'blastn'],
      ['-reward', '4'],
      ['-penalty', '-5'],
      ['-gapopen', '12'],
      ['-gapextend', '8'],
    ]);
    // Megablast's Linear is the default of the lists that hold it; 1,-1 has none (BLAST refuses 0 0 there).
    expect(formParameters(blastn, choose({}, '4,-5'), bn)).toEqual([
      ['-reward', '4'],
      ['-penalty', '-5'],
    ]);
    expect(formParameters(blastn, choose({ '-gapopen': '2', '-gapextend': '1' }, '1,-3'), bn)).toEqual([['-penalty', '-3']]);
    expect(formParameters(blastn, choose({}, '1,-1'), bn)).toEqual([
      ['-penalty', '-1'],
      ['-gapopen', '5'],
      ['-gapextend', '2'],
    ]);
  });

  it('give the gap costs the new default when the task or the matrix changes their list, and keep them otherwise', () => {
    const fourFive = setFields(blastn, {}, pairEntries(blastn, scores!, '4,-5', {}, bn), bn);
    const blastnTask = setField(blastn, fourFive, '-task', 'blastn', bn);
    expect(formParameters(blastn, blastnTask, bn).slice(-2)).toEqual([
      ['-gapopen', '12'],
      ['-gapextend', '8'],
    ]);
    // Back to megablast (the default task, not written): Linear again.
    expect(formParameters(blastn, setField(blastn, blastnTask, '-task', 'megablast', bn), bn)).toEqual([
      ['-reward', '4'],
      ['-penalty', '-5'],
    ]);
    const set: FormValues = { '-gapopen': '10', '-gapextend': '2' };
    expect(formParameters(tblastn, setField(tblastn, set, '-matrix', 'PAM30', tn), tn)).toEqual([['-matrix', 'PAM30']]);
    // The same list: the costs set stay (tblastn-fast scores as tblastn does).
    expect(setField(tblastn, set, '-task', 'tblastn-fast', tn)).toEqual({ ...set, '-task': 'tblastn-fast' });
    // A matrix that NCBI does not list (typed after "Other"): no list, the costs stay.
    expect(setField(tblastn, set, '-matrix', 'IDENTITY', tn)).toEqual({ ...set, '-matrix': 'IDENTITY' });
    // BLASTP's one list: a matrix typed after "Other" leaves the costs.
    expect(setField(blastp, set, '-matrix', 'PAM30', bp)).toEqual({ ...set, '-matrix': 'PAM30' });
    // The row's own fields set here are kept.
    expect(setFields(tblastn, {}, [['-matrix', 'PAM30'], ['-gapopen', '8']], tn)).toEqual({ '-matrix': 'PAM30', '-gapopen': '8' });
    // Scores typed in the two fields: a pair without a list leaves the costs.
    expect(setField(blastn, { '-gapopen': '2', '-gapextend': '1' }, '-reward', '3', bn)).toEqual({ '-gapopen': '2', '-gapextend': '1', '-reward': '3' });
  });

  it("put the low complexity filters in NCBI's Filter row as checkboxes", () => {
    for (const program of PROGRAMS) {
      const filter = program.sections.flatMap((section) => section.fields).find((field) => field.flag === '-dust' || field.flag === '-seg');
      expect(filter, program.id).toMatchObject({ label: 'Low complexity regions', kind: 'filter', row: 'Filter' });
    }
  });
});

describe('fieldRows', () => {
  it("puts the consecutive fields of a row together (NCBI's Match/Mismatch Scores and Gap Costs)", () => {
    const scoring = programById('blastn').sections.find((section) => section.title === 'Scoring Parameters')!;
    expect(fieldRows(scoring.fields).map((row) => [row.row, row.fields.map((f) => f.label)])).toEqual([
      ['Match/Mismatch Scores', ['Match', 'Mismatch']],
      ['Gap Costs', ['Existence', 'Extension']],
    ]);
    const blastp = programById('blastp').sections.find((section) => section.title === 'Scoring Parameters')!;
    expect(fieldRows(blastp.fields).map((row) => [row.row, row.fields.map((f) => f.flag)])).toEqual([
      [undefined, ['-matrix']],
      ['Gap Costs', ['-gapopen', '-gapextend']],
      [undefined, ['-comp_based_stats']],
    ]);
  });
});

describe('searchSummary', () => {
  it('says what is searched, with which program and task', () => {
    expect(searchSummary(programById('blastn'), 'megablast')).toBe(
      'Search nucleotide subjects with BLASTN, task megablast (optimized for highly similar sequences).',
    );
    expect(searchSummary(programById('blastn'), 'dc-megablast')).toContain('task dc-megablast (optimized for more dissimilar sequences)');
    expect(searchSummary(programById('blastp'), 'blastp')).toBe(
      'Search protein subjects with BLASTP, task blastp (protein-protein BLAST).',
    );
    expect(searchSummary(programById('tblastn'), 'tblastn')).toBe(
      'Search translated nucleotide subjects with TBLASTN (protein query), task tblastn.',
    );
    expect(searchSummary(programById('tblastx'), undefined)).toBe(
      'Search translated nucleotide subjects with TBLASTX (translated nucleotide query).',
    );
  });

  it('sets off what is searched and the program with its task', () => {
    expect(searchSummaryParts(programById('blastn'), 'megablast').filter((part) => part.strong).map((part) => part.text)).toEqual([
      'nucleotide subjects',
      'BLASTN, task megablast (optimized for highly similar sequences)',
    ]);
  });

  it('names a task without a note, and leaves the task out while it is unknown', () => {
    expect(searchSummary(programById('blastp'), 'other')).toBe('Search protein subjects with BLASTP, task other.');
    expect(searchSummary(programById('blastn'), undefined)).toBe('Search nucleotide subjects with BLASTN.');
  });
});

describe('regionProblem', () => {
  const problem = (start: string, stop: string, length = 100) => regionProblem({ start, stop }, length);

  it('asks for both positions', () => {
    expect(problem('', '')).toBe('Enter both the start and the stop.');
    expect(problem('1', '')).toBe('Enter both the start and the stop.');
    expect(problem('', '5')).toBe('Enter both the start and the stop.');
    expect(problem('  ', '5')).toBe('Enter both the start and the stop.');
  });

  it.each([
    ['1.5', '10'],
    [' 0x10', '20'],
    ['-3', '10'],
    ['1e3', '10'],
    ['5', 'abc'],
    ['5', '+7'],
    ['1 0', '5'],
  ])('rejects %j and %j as not whole numbers', (start, stop) => {
    expect(problem(start, stop)).toBe('The start and the stop are whole numbers.');
  });

  it('limits the start and the stop to the record', () => {
    expect(problem('0', '10')).toBe('The start must be between 1 and 100, the length of the record.');
    expect(problem('101', '102')).toBe('The start must be between 1 and 100, the length of the record.');
    expect(problem('5', '0')).toBe('The stop must be between 1 and 100, the length of the record.');
    expect(problem('5', '101')).toBe('The stop must be between 1 and 100, the length of the record.');
  });

  it('names the length of the record in the message', () => {
    expect(problem('1', '8', 7)).toBe('The stop must be between 1 and 7, the length of the record.');
  });

  it('accepts positions within the record with the start before the stop', () => {
    expect(problem('1', '100')).toBeUndefined();
    expect(problem('50', '51')).toBeUndefined();
    expect(problem('007', '0100')).toBeUndefined();
  });

  // The engine refuses both, as NCBI does (WA-IN-07): `-query_loc 300-100` "start cannot be
  // larger than stop", `-query_loc 100-100` "range cannot be empty", in every mode and role.
  it('refuses a start after the stop, and a range of one letter, as the engine does', () => {
    expect(problem('80', '20')).toBe('The start cannot be larger than the stop.');
    expect(problem('0100', '99')).toBe('The start cannot be larger than the stop.');
    expect(problem('50', '50')).toBe('The stop must be larger than the start.');
    expect(problem('007', '7')).toBe('The stop must be larger than the start.');
  });

  it('names a position outside the record before the order', () => {
    expect(problem('300', '100')).toBe('The start must be between 1 and 100, the length of the record.');
  });

  it('accepts surrounding white space', () => {
    expect(problem(' 5 ', '\t10\n')).toBeUndefined();
  });
});

describe('swappedRegion', () => {
  const swap = (start: string, stop: string, length = 400) => swappedRegion({ start, stop }, length);

  it('swaps a start after the stop when both are within the record', () => {
    expect(swap('300', '100')).toEqual({ start: '100', stop: '300' });
    expect(swap(' 0300 ', '100\n')).toEqual({ start: '100', stop: '300' });
    expect(swap('400', '1')).toEqual({ start: '1', stop: '400' });
  });

  it('offers nothing for a region that the swap would not make one', () => {
    expect(swap('100', '300')).toBeUndefined();
    expect(swap('100', '100')).toBeUndefined();
    expect(swap('500', '100')).toBeUndefined();
    expect(swap('300', '0')).toBeUndefined();
    expect(swap('3a', '1')).toBeUndefined();
    // Not whole numbers: the problem says so, and a swap would not mend it.
    expect(swap('300.0', '100')).toBeUndefined();
    expect(swap('1e3', '100')).toBeUndefined();
    expect(swap('', '1')).toBeUndefined();
  });
});

describe('region value, range and positions', () => {
  it('has a flag for each role', () => {
    expect(REGION_FLAG).toEqual({ query: '-query_loc', subject: '-subject_loc' });
  });

  it('writes the value without white space or leading zeros', () => {
    expect(regionValue({ start: '007', stop: '0100' })).toBe('7-100');
    expect(regionValue({ start: ' 3 ', stop: ' 9\n' })).toBe('3-9');
  });

  it('gives the range of a valid region', () => {
    expect(regionRange({ start: ' 05', stop: '10 ' }, 20)).toEqual({ start: 5, stop: 10 });
    expect(regionRange({ start: '5', stop: '20' }, 20)).toEqual({ start: 5, stop: 20 });
  });

  it('gives no range while the region is incomplete or outside the record', () => {
    expect(regionRange({ start: '', stop: '5' }, 20)).toBeUndefined();
    expect(regionRange({ start: '1', stop: '21' }, 20)).toBeUndefined();
    expect(regionRange({ start: 'a', stop: '5' }, 20)).toBeUndefined();
    expect(regionRange({ start: '20', stop: '5' }, 20)).toBeUndefined();
    expect(regionRange({ start: '5', stop: '5' }, 20)).toBeUndefined();
  });

  it('orders, rounds and limits positions', () => {
    expect(regionFromPositions(10, 5, 100)).toEqual({ start: '5', stop: '10' });
    expect(regionFromPositions(5.4, 9.6, 100)).toEqual({ start: '5', stop: '10' });
    expect(regionFromPositions(-4, 500, 100)).toEqual({ start: '1', stop: '100' });
    expect(regionFromPositions(0.2, 0.4, 100)).toEqual({ start: '1', stop: '1' });
    expect(regionFromPositions(150, 120, 100)).toEqual({ start: '100', stop: '100' });
  });
});

describe('estimateKind', () => {
  it('estimates a nucleotide sequence from a share of at least 90% of A C G T U N', () => {
    expect(estimateKind({ A: 30, C: 30, G: 20, T: 10, N: 10 })).toBe('nucleotide');
    expect(estimateKind({ a: 9, c: 9, g: 9, t: 9, u: 9, n: 9 })).toBe('nucleotide');
    expect(estimateKind({ A: 90, L: 10 })).toBe('nucleotide');
  });

  it('estimates a protein below that share', () => {
    expect(estimateKind({ A: 89, L: 11 })).toBe('protein');
    expect(estimateKind({ M: 5, K: 5, L: 5, A: 1 })).toBe('protein');
  });

  it('gives no estimate without letters', () => {
    expect(estimateKind({})).toBeUndefined();
    expect(estimateKind({ '*': 3, '-': 2 })).toBeUndefined();
    expect(estimateKind({ A: 0 })).toBeUndefined();
  });

  it('counts only the letters A-Z and a-z', () => {
    expect(estimateKind({ A: 10, '0x00': 1000, '*': 1000, '-': 1000 })).toBe('nucleotide');
    expect(estimateKind({ L: 10, '0x00': 1000, '-': 1000 })).toBe('protein');
    expect(estimateKind({ AC: 100, L: 10 })).toBe('protein');
  });

  it('looks for a kind other than the expected one', () => {
    expect(looksLikeOtherKind({ A: 10, C: 10 }, 'protein')).toBe(true);
    expect(looksLikeOtherKind({ A: 10, C: 10 }, 'nucleotide')).toBe(false);
    expect(looksLikeOtherKind({ M: 10, L: 10 }, 'nucleotide')).toBe(true);
    expect(looksLikeOtherKind({ M: 10, L: 10 }, 'protein')).toBe(false);
  });

  it('does not warn about a record without letters', () => {
    expect(looksLikeOtherKind({}, 'protein')).toBe(false);
    expect(looksLikeOtherKind({ '*': 4 }, 'nucleotide')).toBe(false);
  });
});

describe('geneticCodeLabel', () => {
  it('writes the number and the name of a known code', () => {
    expect(geneticCodeLabel(11)).toBe('11. Bacterial, Archaeal and Plant Plastid');
    expect(geneticCodeLabel(1)).toBe('1. Standard');
  });

  it('writes the number alone for a code without a name', () => {
    expect(geneticCodeLabel(7)).toBe('7');
    expect(geneticCodeLabel(99)).toBe('99');
  });
});

describe('duplicateIds', () => {
  it('returns the IDs that appear more than once', () => {
    const ids = duplicateIds([{ id: 'a' }, { id: 'b' }, { id: 'a' }, { id: 'c' }, { id: 'b' }, { id: 'a' }]);
    expect([...ids].sort()).toEqual(['a', 'b']);
  });

  it('returns an empty set when all IDs differ', () => {
    expect(duplicateIds([{ id: 'a' }, { id: 'b' }]).size).toBe(0);
    expect(duplicateIds([]).size).toBe(0);
  });
});

describe('program descriptors', () => {
  const RESERVED = ['-query', '-subject', '-out', '-outfmt', '-num_threads', '-query_loc', '-subject_loc'];

  it('lists each program once', () => {
    expect(new Set(PROGRAMS.map((p) => p.id)).size).toBe(PROGRAMS.length);
    for (const program of PROGRAMS) expect(programById(program.id)).toBe(program);
  });

  it("reads each input with NCBI's reader kind of its sequence kind (DW-38)", () => {
    expect(Object.fromEntries(PROGRAMS.map((p) => [p.id, [fastaParser(p, 'query'), fastaParser(p, 'subject')]]))).toEqual({
      blastn: [1, 1],
      blastp: [2, 2],
      tblastn: [2, 1],
      tblastx: [1, 1],
      blastx: [1, 2],
    });
  });

  it.each(PROGRAMS.map((p) => [p.id, p] as const))('%s has a form with unique flags', (_id, program) => {
    expect(program.sections.length).toBeGreaterThan(0);
    for (const section of program.sections) expect(section.fields.length).toBeGreaterThan(0);
    const flags = program.sections.flatMap((s) => s.fields.map((f) => f.flag));
    expect(new Set(flags).size).toBe(flags.length);
    for (const flag of flags) expect(flag).toMatch(/^-[a-z_]+$/);
  });

  it.each(PROGRAMS.map((p) => [p.id, p] as const))('%s uses no reserved flag', (_id, program) => {
    const flags = program.sections.flatMap((s) => s.fields.map((f) => f.flag));
    for (const flag of flags) expect(RESERVED).not.toContain(flag);
  });

  it.each(PROGRAMS.map((p) => [p.id, p] as const))("%s has NCBI's sections in NCBI's order under Algorithm parameters", (_id, program) => {
    const order = ['General Parameters', 'Scoring Parameters', 'Filters and Masking', 'Discontiguous Word Options', 'Other Parameters'];
    const titles = program.sections.filter((s) => s.placement === 'algorithm').map((s) => s.title);
    expect(titles.map((title) => order.indexOf(title))).toEqual([...titles.map((title) => order.indexOf(title))].sort((a, b) => a - b));
    expect(titles.every((title) => order.includes(title))).toBe(true);
    // The task is chosen with radio buttons in Task Selection, before every other option.
    const task = program.sections.flatMap((s) => s.fields).find((f) => f.flag === '-task');
    if (task !== undefined) {
      expect(task.kind).toBe('radio');
      expect(program.sections[0]!.placement).toBe('program');
      expect(program.sections[0]!.fields[0]).toBe(task);
    }
  });

  it('gives the sequence kind and unit of each input', () => {
    expect(sequenceKind(programById('tblastn'), 'query')).toBe('protein');
    expect(sequenceKind(programById('tblastn'), 'subject')).toBe('nucleotide');
    expect(residueUnit('nucleotide')).toBe('nt');
    expect(residueUnit('protein')).toBe('aa');
  });
});
