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
import CONVERT2BLASTMASK_XMLHELP_DOC from '../../src/infra/fake/convert2blastmask-xmlhelp.json';
import DUSTMASKER_XMLHELP_DOC from '../../src/infra/fake/dustmasker-xmlhelp.json';
import SEGMASKER_XMLHELP_DOC from '../../src/infra/fake/segmasker-xmlhelp.json';
import MAKEBLASTDB_XMLHELP_DOC from '../../src/infra/fake/makeblastdb-xmlhelp.json';
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
      'LOSAT/tests/dustmasker/cases.tsv',
      'LOSAT/tests/segmasker/cases.tsv',
      'LOSAT/tests/convert2blastmask/cases.tsv',
      'LOSAT/tests/dropin_corpus_manifest.tsv',
      'LOSAT/tests/fixtures/asn1/makeblastdb/ncbi_2_17_0_manifest.tsv',
      'LOSAT/tests/seqdb/gi_mask/cases.tsv',
      'LOSAT/tests/non_utf8_names/cases.tsv',
      'LOSAT/src/writedb/tests/makeblastdb_v5_cases.rs',
      'LOSAT/src/writedb/tests/makeblastdb_parse_cases.rs',
      'LOSAT/src/writedb/tests/makeblastdb_fix_cases.rs',
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
    const bfBadge = (argv: readonly string[]) =>
      toolVerificationBadge({ tool: 'blast_formatter', argv, runtimePath: 'serial', grammar: bfGrammar }, table);
    for (const argv of [
      ['blast_formatter', '-archive', 'bn.asn', '-outfmt', '8'],
      ['blast_formatter', '-archive', 'db.asn', '-outfmt', '9', '-out', 'out.asnb'],
      ['blast_formatter', '-archive', 'db.asn', '-outfmt', '12', '-logfile', 'run.log'],
      ['blast_formatter', '-archive', 'db.asn', '-outfmt', '8', '-max_target_seqs', '2'],
      ['blast_formatter', '-archive', 'db.asn', '-outfmt', '8', '-parse_deflines'],
    ]) {
      const badge = bfBadge(argv);
      expect(badge.level, argv.join(' ')).toBe('certified');
      expect(badge.label, argv.join(' ')).toBe(SUPPORTED_LABEL);
    }

    for (const argv of [
      ['blast_formatter', '-archive', 'bn.asn', '-outfmt', '6'],
      ['blast_formatter', '-archive', 'bn.asn', '-outfmt', '8', '-line_length', '40'],
    ]) {
      const badge = bfBadge(argv);
      expect(badge.level, argv.join(' ')).toBe('outside');
      expect(badge.label, argv.join(' ')).toBe(NOT_SUPPORTED_LABEL);
    }
  });

  it('generates the option sets of the Masking tools from cases.tsv, the default option set included', () => {
    const grade = (tool: 'dustmasker' | 'segmasker' | 'convert2blastmask', xml: string, argv: readonly string[]) =>
      toolVerificationBadge({ tool, argv, runtimePath: 'serial', grammar: toolOptionGrammar(parseToolDescription(xml)) }, table);
    for (const tool of ['dustmasker', 'segmasker', 'convert2blastmask'] as const) {
      const sets = table.tools?.[tool]?.optionSets;
      expect(sets, tool).toBeDefined();
      expect([...sets!]).toEqual([...sets!].sort());
      for (const key of sets!) expect(key, tool).not.toMatch(/"-(in|out|logfile)"/);
    }
    // A run with no option besides -in and -out is the tool's default option set, which cases.tsv compares (basic_int).
    expect(table.tools?.dustmasker?.optionSets).toContain('[]');
    for (const argv of [
      ['dustmasker', '-in', 'input.fa'],
      ['dustmasker', '-in', 'input.fa', '-out', 'masked.txt', '-outfmt', 'interval'],
    ]) {
      const badge = grade('dustmasker', DUSTMASKER_XMLHELP_DOC.xml, argv);
      expect(badge.level, argv.join(' ')).toBe('certified');
      expect(badge.label, argv.join(' ')).toBe(SUPPORTED_LABEL);
    }
    expect(grade('segmasker', SEGMASKER_XMLHELP_DOC.xml, ['segmasker', '-in', 'input.fa', '-window', '9999']).level).toBe('outside');
    expect(grade('convert2blastmask', CONVERT2BLASTMASK_XMLHELP_DOC.xml, ['convert2blastmask', '-in', 'input.fa']).level).toBe('outside');
  });

  it("generates makeblastdb option sets from every NCBI-frozen ledger of makeblastdb runs and grades runs (line WW, WW4b2 R1)", () => {
    const mk = table.tools?.makeblastdb;
    expect(mk).toBeDefined();
    expect([...mk!.optionSets]).toEqual([...mk!.optionSets].sort());
    for (const key of mk!.optionSets) expect(key).not.toMatch(/"-(in|out|logfile)"/);
    // The runs NCBI ended with exit status 0 that each ledger's test compares (the counts of 2026-11-28; a ledger
    // grows, and a reader that stops reading its rows shows here).
    const taken = Object.fromEntries(table.sources.map((s) => [s.path, s.optionSets]));
    expect(taken['LOSAT/tests/dropin_corpus_manifest.tsv']).toBeGreaterThanOrEqual(57);
    expect(taken['LOSAT/tests/fixtures/asn1/makeblastdb/ncbi_2_17_0_manifest.tsv']).toBeGreaterThanOrEqual(59);
    expect(taken['LOSAT/tests/seqdb/gi_mask/cases.tsv']).toBeGreaterThanOrEqual(23);
    expect(taken['LOSAT/tests/non_utf8_names/cases.tsv']).toBeGreaterThanOrEqual(8);
    expect(taken['LOSAT/src/writedb/tests/makeblastdb_v5_cases.rs']).toBeGreaterThanOrEqual(86);
    expect(taken['LOSAT/src/writedb/tests/makeblastdb_parse_cases.rs']).toBeGreaterThanOrEqual(9);
    expect(taken['LOSAT/src/writedb/tests/makeblastdb_fix_cases.rs']).toBeGreaterThanOrEqual(5);

    const mkGrammar = toolOptionGrammar(parseToolDescription(MAKEBLASTDB_XMLHELP_DOC.xml));
    const mkBadge = (argv: readonly string[]) =>
      toolVerificationBadge({ tool: 'makeblastdb', argv, runtimePath: 'serial', grammar: mkGrammar }, table);
    for (const argv of [
      // The drop-in ledger.
      ['makeblastdb', '-in', 'input.fa', '-dbtype', 'nucl'],
      ['makeblastdb', '-in', 'input.fa', '-dbtype', 'prot', '-out', 'mydb'],
      ['makeblastdb', '-in', 'a.fa b.fa', '-dbtype', 'nucl', '-parse_seqids', '-out', 'db'],
      ['makeblastdb', '-in', 'input.fa', '-dbtype', 'nucl', '-parse_seqids', '-taxid_map', 'map.txt', '-out', 'db'],
      ['makeblastdb', '-in', 'input.fa', '-dbtype', 'nucl', '-title', 'my genome', '-out', 'db'],
      ['makeblastdb', '-in', 'input.fa', '-dbtype', 'nucl', '-blastdb_version', '5'],
      ['makeblastdb', '-in', 'input.fa', '-dbtype', 'prot', '-mask_data', 'seg.asnb'],
      ['makeblastdb', '-in', 'input.fa', '-dbtype', 'nucl', '-parse_seqids', '-mask_data', 'a.asnb,b.asnb', '-out', 'db'],
      ['makeblastdb', '-in', 'input.fa', '-dbtype', 'prot', '-mask_data', 'm.asnb', '-mask_id', 'mine', '-mask_desc', 'my seg', '-out', 'db'],
      ['makeblastdb', '-in', 'input.fa', '-dbtype', 'nucl', '-blastdb_version', '4', '-mask_data', 'dust.asnb'],
      // The ASN.1 manifest.
      ['makeblastdb', '-in', 'e.asnb', '-input_type', 'asn1_bin', '-dbtype', 'prot', '-title', 'T', '-out', 'db'],
      ['makeblastdb', '-in', 'e.asn', '-input_type', 'asn1_txt', '-dbtype', 'nucl', '-title', 'x', '-oid_masks', '1', '-out', 'db'],
      // -gi_mask in several volumes.
      ['makeblastdb', '-in', 'n.fa', '-dbtype', 'nucl', '-parse_seqids', '-title', 't', '-max_file_sz', '400', '-mask_data', 'd.asnb', '-gi_mask', '-gi_mask_name', 'g1'],
      // Names that are not UTF-8: -metadata_output_prefix.
      ['makeblastdb', '-in', 'in.fa', '-dbtype', 'nucl', '-title', 't', '-out', 'db', '-metadata_output_prefix', 'meta/'],
      // The version 5 runs of the Rust tables: -input_type blastdb, -randomize, -oid_masks, volumes.
      ['makeblastdb', '-in', 'src', '-input_type', 'blastdb', '-dbtype', 'nucl', '-out', 'db'],
      ['makeblastdb', '-in', 'src', '-input_type', 'blastdb', '-dbtype', 'prot', '-out', 'db', '-randomize'],
      ['makeblastdb', '-in', 'p.fa', '-dbtype', 'prot', '-out', 'db', '-title', 't', '-parse_seqids', '-oid_masks', '2'],
      ['makeblastdb', '-in', 'p.fa', '-dbtype', 'prot', '-out', 'db', '-title', 't', '-max_file_sz', '300B'],
      // The -parse_seqids runs: -hash_index with -taxid_map.
      ['makeblastdb', '-in', 'p.fa', '-dbtype', 'prot', '-out', 'db', '-title', 't', '-parse_seqids', '-hash_index', '-taxid_map', 'm.txt'],
    ]) {
      const badge = mkBadge(argv);
      expect(badge.level, argv.join(' ')).toBe('certified');
      expect(badge.label, argv.join(' ')).toBe(SUPPORTED_LABEL);
    }
    for (const argv of [
      ['makeblastdb', '-in', 'input.fa', '-dbtype', 'nucl', '-max_file_sz', '1MB'],
      // A single-volume database with -gi_mask: NCBI ends it with exit status 2 (approved exception 25), so no row counts.
      ['makeblastdb', '-in', 'n.fa', '-dbtype', 'nucl', '-parse_seqids', '-title', 't', '-mask_data', 'd.asnb', '-gi_mask'],
      ['makeblastdb', '-in', 'n.fa', '-dbtype', 'nucl', '-parse_seqids', '-title', 't', '-mask_data', 'd.asnb', '-gi_mask', '-gi_mask_name', 'g1'],
      // -mask_id without -mask_data, -mask_desc in version 4 with several volumes: only in mask_writer, whose test leaves the metadata file uncompared.
      ['makeblastdb', '-in', 'n.fa', '-dbtype', 'nucl', '-blastdb_version', '4', '-title', 't', '-mask_data', 'm.asnt', '-mask_id', 'a', '-max_file_sz', '400B'],
      // The source builds of the Rust tables' -input_type blastdb cases (a folder of {D}, not compared).
      ['makeblastdb', '-in', 'raw.fa', '-dbtype', 'nucl', '-parse_seqids', '-taxid_map', 't.map', '-mask_data', 'd.asnb', '-title', 'S', '-out', 's/src'],
    ]) {
      const badge = mkBadge(argv);
      expect(badge.level, argv.join(' ')).toBe('outside');
      expect(badge.label, argv.join(' ')).toBe(NOT_SUPPORTED_LABEL);
    }
  });
});
