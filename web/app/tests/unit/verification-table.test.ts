import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { dirname, resolve } from 'node:path';
import { parseToolDescription } from 'losat/domain';
import { describe, expect, it } from 'vitest';
import { generateVerificationTable } from '../../build/verification';
import {
  FILE_VALUE,
  NOT_SUPPORTED_LABEL,
  optionKey,
  RANGE_VALUE,
  SUPPORTED_LABEL,
  toolOptionGrammar,
  toolOptionKey,
  toolVerificationBadge,
  type OptionGrammar,
} from '../../src/domain/verification';
import BLAST_FORMATTER_XMLHELP_DOC from '../../src/infra/fake/blast_formatter-xmlhelp.json';
import WINDOWMASKER_XMLHELP_DOC from '../../src/infra/fake/windowmasker-xmlhelp.json';

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
    // task (megablast), of default BLASTX and of tblastn-fast are the def7.* rows of outfmt0_manifest.tsv;
    // blastp-short's are score.blastp.short.def.* (series SC, SC3); rmblastn's are rm.def.* (line RG X6).
    for (const program of PROGRAMS) expect(formatsOf(program, []), program).toEqual([0, 6, 7]);
    const tasks = [
      ['blastn', 'megablast'],
      ['blastn', 'blastn'],
      ['blastn', 'dc-megablast'],
      ['blastn', 'blastn-short'],
      ['blastn', 'rmblastn'],
      ['blastp', 'blastp-fast'],
      ['blastp', 'blastp-short'],
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
      'LOSAT/tests/windowmasker/cases.tsv',
      'LOSAT/tests/blast_formatter/cases.tsv',
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

  it('generates windowmasker option sets from cases.tsv and grades tool runs as Supported or Not yet supported', () => {
    const wm = table.tools?.windowmasker;
    expect(wm).toBeDefined();
    expect(wm!.optionSets.length).toBeGreaterThan(100);
    expect([...wm!.optionSets]).toEqual([...wm!.optionSets].sort());
    for (const key of wm!.optionSets) {
      expect(key).not.toMatch(/"-(in|out|logfile)"/);
    }

    const wmGrammar = toolOptionGrammar(parseToolDescription(WINDOWMASKER_XMLHELP_DOC.xml));
    // File names (-in, -out, -logfile, -ustat value, -ids value) do not change the key.
    expect(
      toolOptionKey(['windowmasker', '-ustat', 'windowmasker/562/wmasker.obinary', '-in', 'a.fa', '-out', 'b.txt', '-logfile', 'run.log'], wmGrammar),
    ).toBe(JSON.stringify([['-ustat', FILE_VALUE]]));
    expect(toolOptionKey(['windowmasker', '-ustat', 'c.obinary', '-dust', 'true'], wmGrammar)).toBe(
      toolOptionKey(['windowmasker', '-ustat', 'other.bin', '-dust', 'T'], wmGrammar),
    );
    expect(toolOptionKey(['windowmasker', '-convert', '-in', 'a.ascii', '-out', 'b.ascii', '-sformat', 'ascii'], wmGrammar)).toBe(
      toolOptionKey(['windowmasker', '-convert', '-in', 'x.obinary'], wmGrammar),
    );

    // Frozen option sets -> Supported.
    for (const argv of [
      ['windowmasker', '-mk_counts', '-in', 'genome.fa', '-out', 'wmasker.ascii'],
      ['windowmasker', '-mk_counts', '-in', 'genome.fa', '-out', 'wmasker.obinary', '-unit', '5', '-sformat', 'obinary'],
      ['windowmasker', '-ustat', 'windowmasker/562/wmasker.obinary', '-in', 'seq.fa', '-outfmt', 'fasta', '-parse_seqids'],
      ['windowmasker', '-ustat', 'counts.obinary', '-in', 'seq.fa', '-dust', 'true'],
      ['windowmasker', '-convert', '-in', 'wmasker.ascii', '-out', 'wmasker.obinary', '-sformat', 'obinary'],
    ]) {
      const badge = toolVerificationBadge({ tool: 'windowmasker', argv, runtimePath: 'serial', grammar: wmGrammar }, table);
      expect(badge.level, argv.join(' ')).toBe('certified');
      expect(badge.label, argv.join(' ')).toBe(SUPPORTED_LABEL);
      expect(badge.details.join(' ')).not.toMatch(/\b(verif|certif)/i);
    }

    // Option sets without an NCBI-frozen exit-0 row -> Not yet supported.
    for (const argv of [
      ['windowmasker', '-mk_counts', '-in', 'genome.fa', '-out', 'wmasker.ascii', '-unit', '11'],
      ['windowmasker', '-ustat', 'counts.obinary', '-in', 'seq.fa', '-window', '999'],
      ['windowmasker', '-mk_counts', '-in', 'genome.fa', '-conffile', 'app.ini'],
      ['windowmasker', '-mk_counts', '-in', 'genome.fa', '-dryrun'],
    ]) {
      const badge = toolVerificationBadge({ tool: 'windowmasker', argv, runtimePath: 'serial', grammar: wmGrammar }, table);
      expect(badge.level, argv.join(' ')).toBe('outside');
      expect(badge.label, argv.join(' ')).toBe(NOT_SUPPORTED_LABEL);
      expect(badge.details.join(' ')).not.toMatch(/\b(verif|certif)/i);
    }

    // Development build (FakeEngine) -> Development build badge.
    const dev = toolVerificationBadge(
      { tool: 'windowmasker', argv: ['windowmasker', '-mk_counts', '-in', 'a.fa', '-out', 'w.ascii'], runtimePath: 'fake', grammar: wmGrammar },
      table,
    );
    expect(dev).toMatchObject({ level: 'development', label: 'Development build' });
  });

  it('generates blast_formatter option sets from cases.tsv and grades runs as Supported or Not yet supported', () => {
    const bf = table.tools?.blast_formatter;
    expect(bf).toBeDefined();
    expect(bf!.optionSets.length).toBeGreaterThanOrEqual(3);
    expect([...bf!.optionSets]).toEqual([...bf!.optionSets].sort());
    for (const key of bf!.optionSets) {
      expect(key).not.toMatch(/"-(archive|out|logfile)"/);
    }

    const bfGrammar = toolOptionGrammar(parseToolDescription(BLAST_FORMATTER_XMLHELP_DOC.xml));
    for (const argv of [
      ['blast_formatter', '-archive', 'bn.asn', '-outfmt', '8'],
      ['blast_formatter', '-archive', 'db.asn', '-outfmt', '9', '-out', 'out.asnb'],
      ['blast_formatter', '-archive', 'db.asn', '-outfmt', '12', '-logfile', 'run.log'],
      ['blast_formatter', '-archive', 'db.asn', '-outfmt', '8', '-max_target_seqs', '2'],
      ['blast_formatter', '-archive', 'db.asn', '-outfmt', '8', '-parse_deflines'],
    ]) {
      const badge = toolVerificationBadge({ tool: 'blast_formatter', argv, runtimePath: 'serial', grammar: bfGrammar }, table);
      expect(badge.level, argv.join(' ')).toBe('certified');
      expect(badge.label, argv.join(' ')).toBe(SUPPORTED_LABEL);
    }

    for (const argv of [
      ['blast_formatter', '-archive', 'bn.asn', '-outfmt', '6'],
      ['blast_formatter', '-archive', 'bn.asn', '-outfmt', '8', '-line_length', '40'],
    ]) {
      const badge = toolVerificationBadge({ tool: 'blast_formatter', argv, runtimePath: 'serial', grammar: bfGrammar }, table);
      expect(badge.level, argv.join(' ')).toBe('outside');
      expect(badge.label, argv.join(' ')).toBe(NOT_SUPPORTED_LABEL);
    }
  });
});
