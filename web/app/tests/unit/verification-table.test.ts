import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { generateVerificationTable } from '../../build/verification';
import { optionKey, RANGE_VALUE, type OptionGrammar } from '../../src/domain/verification';

const repository = resolve(import.meta.dirname, '../../../..');
const describeJson = JSON.parse(readFileSync(resolve(import.meta.dirname, '../../src/infra/fake/describe.json'), 'utf8')) as Record<
  string,
  { parameters: { flag: string; takes_value: boolean; default?: string }[] }
>;

function grammar(program: string): OptionGrammar {
  const parameters = describeJson[program]!.parameters;
  const defaultTask = parameters.find((p) => p.flag === '-task')?.default;
  return {
    takesValue: (flag) => parameters.find((p) => p.flag === flag)?.takes_value ?? false,
    ...(defaultTask === undefined ? {} : { defaultTask }),
  };
}

const PROGRAMS = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx'] as const;
const table = generateVerificationTable(repository);
const formatsOf = (program: (typeof PROGRAMS)[number], words: string[]) => {
  const key = optionKey([program, ...words], grammar(program));
  return key === undefined ? undefined : table.programs[program]?.optionSets[key];
};

describe('verification table', () => {
  it('lists the five programs of the engine and the NCBI version', () => {
    expect(table.ncbi).toBe('2.17.0');
    expect(Object.keys(table.programs).sort()).toEqual([...PROGRAMS]);
  });

  it('records the browser runtime of every program from the verification cells', () => {
    for (const program of PROGRAMS) {
      const browser = table.programs[program]?.browser;
      expect(browser?.formats, program).toEqual([0, 6, 7]);
      expect(browser?.paths, program).toEqual(['serial', 'threaded']);
      // V-BR runs every count that Auto and the Threads list can choose (WA-RS-03).
      expect(browser?.threads, program).toEqual(Array.from({ length: 32 }, (_, i) => i + 1));
      expect(browser?.browsers, program).toEqual(expect.arrayContaining(['Chromium', 'Firefox', 'WebKit']));
    }
  });

  it('has the default options of every program and task compared in outfmt 0, 6 and 7', () => {
    // A run with the defaults is Supported (WA-RS-03): the outfmt 7 rows of the default BLASTN
    // task (megablast), of default BLASTX and of tblastn-fast are the def7.* rows of outfmt0_manifest.tsv.
    for (const program of PROGRAMS) expect(formatsOf(program, []), program).toEqual([0, 6, 7]);
    const tasks = [
      ['blastn', 'megablast'],
      ['blastn', 'blastn'],
      ['blastn', 'dc-megablast'],
      ['blastn', 'blastn-short'],
      ['blastp', 'blastp-fast'],
      ['blastx', 'blastx-fast'],
      ['tblastn', 'tblastn-fast'],
    ] as const;
    for (const [program, task] of tasks) expect(formatsOf(program, ['-task', task]), `${program} -task ${task}`).toEqual([0, 6, 7]);
  });

  it('has the options of the regression fixtures and the approved exception', () => {
    expect(formatsOf('blastn', ['-task', 'blastn'])).toContain(0);
    const gencode4 = Object.entries(table.programs.tblastx?.optionSets ?? {}).filter(([key]) => key.includes('["-db_gencode","4"]'));
    expect(gencode4.length).toBeGreaterThan(0);
  });

  it('keeps a flag of a range and drops the range', () => {
    const keys = PROGRAMS.flatMap((program) => Object.keys(table.programs[program]?.optionSets ?? {}));
    const ranged = keys.filter((key) => key.includes('"-query_loc"'));
    expect(ranged.length).toBeGreaterThan(0);
    for (const key of ranged) expect(key).toContain(`["-query_loc","${RANGE_VALUE}"]`);
    // BLASTX's ranges with the default options (line SX, DW-39): sxd.blastx.qloc.*, sloc.0, qsloc.7.
    expect(formatsOf('blastx', ['-query_loc', '1001-3000'])).toEqual([0, 6, 7]);
    expect(formatsOf('blastx', ['-subject_loc', '10-70'])).toContain(0);
    expect(formatsOf('blastx', ['-query_loc', '1-9', '-subject_loc', '2-5'])).toEqual([7]);
  });

  it('keeps nothing outside the inputs, outputs and threads in a key', () => {
    for (const program of PROGRAMS) {
      for (const [key, formats] of Object.entries(table.programs[program]?.optionSets ?? {})) {
        expect(key).not.toMatch(/"-(query|subject|out|outfmt|num_threads)"/);
        expect(formats.length).toBeGreaterThan(0);
        expect([...formats]).toEqual([...formats].sort((a, b) => a - b));
      }
    }
  });

  it('names its sources with their SHA-256 and the records that it took from them', () => {
    expect(table.sources.map((s) => s.path)).toEqual([
      'LOSAT/tests/outfmt0_manifest.tsv',
      'LOSAT/tests/fixtures/blastn_regression/manifest.tsv',
      'LOSAT/tests/fixtures/tblastx_regression/manifest.tsv',
      'LOSAT/tests/fixtures/range_regression/manifest.tsv',
      'LOSAT/tests/blastp_v010_parity_manifest.tsv',
      'LOSAT/tests/tblastx_v010_parity_manifest.tsv',
      'docs/evidence/tlosan_stage_g/matrix_162.jsonl',
      'docs/evidence/losat_web_e2e/run-20261004T163746Z/sweeps/after-blastp.tsv',
      'docs/evidence/losat_web_e2e/run-20261004T163746Z/sweeps/after-tblastn.tsv',
      'docs/evidence/losat_web_e2e/run-20261004T163746Z/sweeps/after-tblastx.tsv',
    ]);
    for (const source of table.sources) {
      expect(source.sha256, source.path).toMatch(/^[0-9a-f]{64}$/);
      // Every ranged search of the range manifest is an error exit or needs an environment
      // variable, so it adds no byte comparison.
      if (source.path.includes('range_regression')) expect(source.optionSets).toBe(0);
      else expect(source.optionSets, source.path).toBeGreaterThan(0);
    }
  });

  it('is the same on every run', () => {
    expect(JSON.stringify(generateVerificationTable(repository))).toBe(JSON.stringify(table));
  });

  it('writes a copy for the gate record when LOSAT_WEB_VERIFICATION_OUT is set', () => {
    const out = process.env.LOSAT_WEB_VERIFICATION_OUT;
    if (out === undefined || out === '') return;
    mkdirSync(dirname(resolve(out)), { recursive: true });
    writeFileSync(out, `${JSON.stringify(table, null, 2)}\n`);
  });
});
