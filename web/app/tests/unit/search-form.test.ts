import { describe, expect, it } from 'vitest';
import { duplicateIds } from '../../src/domain/dataset';
import { geneticCodeLabel } from '../../src/domain/genetic-codes';
import {
  describedSections,
  fieldChoices,
  formParameters,
  setField,
  type EngineOption,
  type FormValues,
} from '../../src/domain/parameters';
import { PROGRAMS, programById, residueUnit, sequenceKind } from '../../src/domain/programs';
import {
  REGION_FLAG,
  regionFromPositions,
  regionProblem,
  regionRange,
  regionValue,
} from '../../src/domain/region';
import { estimateKind, looksLikeOtherKind } from '../../src/domain/sequence-kind';

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
    const values: FormValues = { '-lcase_masking': true, '-word_size': '11', '-evalue': '1', '-task': 'blastn' };
    expect(formParameters(blastn, values, options)).toEqual([
      ['-task', 'blastn'],
      ['-evalue', '1'],
      ['-word_size', '11'],
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

  it('works for TBLASTX with its sections in order', () => {
    const tblastx = programById('tblastx');
    const values: FormValues = { '-db_gencode': '11', '-query_gencode': '4', '-seg': 'yes', '-culling_limit': '2', '-evalue': '1' };
    expect(formParameters(tblastx, values, undefined)).toEqual([
      ['-evalue', '1'],
      ['-culling_limit', '2'],
      ['-seg', 'yes'],
      ['-query_gencode', '4'],
      ['-db_gencode', '11'],
    ]);
  });

  it('omits a TBLASTX value equal to the engine default', () => {
    const tblastx = programById('tblastx');
    const opts: EngineOption[] = [valueOption('-query_gencode', '1'), valueOption('-db_gencode', '1')];
    expect(formParameters(tblastx, { '-query_gencode': '1', '-db_gencode': '2' }, opts)).toEqual([['-db_gencode', '2']]);
  });
});

describe('setField', () => {
  const templates: FormValues = { '-task': 'dc-megablast', '-template_type': 'coding', '-template_length': '18', '-evalue': '1' };

  it.each(['blastn', 'blastn-short'])('removes the template fields when the task is %s', (task) => {
    const next = setField(templates, '-task', task);
    expect(next).toEqual({ '-task': task, '-evalue': '1' });
    expect('-template_type' in next).toBe(false);
    expect('-template_length' in next).toBe(false);
  });

  it.each(['dc-megablast', 'megablast'])('keeps the template fields when the task is %s', (task) => {
    expect(setField(templates, '-task', task)).toEqual({ ...templates, '-task': task });
  });

  it('does not clear the templates when another field is set', () => {
    expect(setField(templates, '-evalue', '5')).toEqual({ ...templates, '-evalue': '5' });
  });

  it('returns a frozen object and leaves the input unchanged', () => {
    const before = { ...templates };
    expect(Object.isFrozen(setField(templates, '-evalue', '5'))).toBe(true);
    expect(Object.isFrozen(setField(templates, '-task', 'blastn'))).toBe(true);
    expect(templates).toEqual(before);
  });

  it('sets a flag value', () => {
    expect(setField({}, '-lcase_masking', true)).toEqual({ '-lcase_masking': true });
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
      ['General parameters', ['-evalue']],
      ['Scoring parameters', ['-matrix']],
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

  it('accepts positions within the record, including start == stop and start > stop', () => {
    expect(problem('1', '100')).toBeUndefined();
    expect(problem('50', '50')).toBeUndefined();
    expect(problem('80', '20')).toBeUndefined();
    expect(problem('007', '0100')).toBeUndefined();
  });

  it('accepts surrounding white space', () => {
    expect(problem(' 5 ', '\t10\n')).toBeUndefined();
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
    expect(regionRange({ start: '20', stop: '5' }, 20)).toEqual({ start: 20, stop: 5 });
  });

  it('gives no range while the region is incomplete or outside the record', () => {
    expect(regionRange({ start: '', stop: '5' }, 20)).toBeUndefined();
    expect(regionRange({ start: '1', stop: '21' }, 20)).toBeUndefined();
    expect(regionRange({ start: 'a', stop: '5' }, 20)).toBeUndefined();
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
  const available = PROGRAMS.filter((p) => p.unavailable === undefined);

  it('lists each program once', () => {
    expect(new Set(PROGRAMS.map((p) => p.id)).size).toBe(PROGRAMS.length);
    for (const program of PROGRAMS) expect(programById(program.id)).toBe(program);
  });

  it('marks BLASTX as unavailable', () => {
    expect(programById('blastx').unavailable).toBeTypeOf('string');
    expect(available.map((p) => p.id)).not.toContain('blastx');
  });

  it.each(available.map((p) => [p.id, p] as const))('%s has a form with unique flags', (_id, program) => {
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

  it('gives the sequence kind and unit of each input', () => {
    expect(sequenceKind(programById('tblastn'), 'query')).toBe('protein');
    expect(sequenceKind(programById('tblastn'), 'subject')).toBe('nucleotide');
    expect(residueUnit('nucleotide')).toBe('nt');
    expect(residueUnit('protein')).toBe('aa');
  });
});
