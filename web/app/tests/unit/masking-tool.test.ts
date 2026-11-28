// Unit tests for the Masking tool tab (WW4c): the domain helpers of dustmasker, segmasker and
// convert2blastmask (options from -xmlhelp, argv builder, the masking-information rule), MaskingDraft
// (argv per tool, per-tool values, -out, the shared FASTA input, blastdb, required options, help and
// version, the `-outfmt` default that a run records), the one rule for mask data (`maskDataName`) with the
// Coordinator's and the makeblastdb draft's hand-over, and the FakeEngine's rows for the three tools.
import { readFileSync } from 'node:fs';
import { gzipSync } from 'node:zlib';
import { describe, expect, it, vi } from 'vitest';
import { sha256Hex } from 'losat';
import { parseToolDescription } from 'losat/domain';
import type { ToolEngine } from 'losat/ports';
import { Coordinator, isToolRun, type EnqueueAllResult, type ToolRequest } from '../../src/application/coordinator';
import { MakeblastdbDraft } from '../../src/application/makeblastdb-draft';
import { MaskingDraft } from '../../src/application/masking-draft';
import { SessionFiles } from '../../src/application/session-files';
import { ToolDraft } from '../../src/application/tool-draft';
import type { FormValues } from '../../src/domain/parameters';
import {
  buildBlastFormatterArgv,
  buildMakeblastdbArgv,
  buildMaskingArgv,
  buildWindowMaskerArgv,
  changedMaskingOptionsCount,
  isCountsOutputFile,
  isStandardOutputName,
  maskDataName,
  maskOutfmtDefault,
  runFileName,
  toolFileOptionLabel,
} from '../../src/domain/tool-argv';
import {
  allMaskingParameters,
  isMaskingTool,
  MASKING_TOOLS,
  maskingOptionParameters,
  maskingRequiredParameters,
  maskingToolSummary,
  toolLabel,
  toolPrintsKind,
  toolRunTitle,
  toolStageLabel,
  toToolParameterField,
  TOOL_WRITTEN_TEXT_FLAGS,
} from '../../src/domain/tools';
import { countsName } from '../../src/domain/counts';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FAKE_MARKER, FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import CONVERT2BLASTMASK_XMLHELP_DOC from '../../src/infra/fake/convert2blastmask-xmlhelp.json';
import DUSTMASKER_XMLHELP_DOC from '../../src/infra/fake/dustmasker-xmlhelp.json';
import SEGMASKER_XMLHELP_DOC from '../../src/infra/fake/segmasker-xmlhelp.json';
import WINDOWMASKER_XMLHELP_DOC from '../../src/infra/fake/windowmasker-xmlhelp.json';
import BLAST_FORMATTER_XMLHELP_DOC from '../../src/infra/fake/blast_formatter-xmlhelp.json';
import MAKEBLASTDB_XMLHELP_DOC from '../../src/infra/fake/makeblastdb-xmlhelp.json';
import { NodeInflater } from './support/node-inflater';

const DUST = parseToolDescription(DUSTMASKER_XMLHELP_DOC.xml);
const SEG = parseToolDescription(SEGMASKER_XMLHELP_DOC.xml);
const C2B = parseToolDescription(CONVERT2BLASTMASK_XMLHELP_DOC.xml);
const WM = parseToolDescription(WINDOWMASKER_XMLHELP_DOC.xml);
const BF = parseToolDescription(BLAST_FORMATTER_XMLHELP_DOC.xml);
const MKDB = parseToolDescription(MAKEBLASTDB_XMLHELP_DOC.xml);
const sleep = (ms: number) => new Promise<void>((resolve) => setTimeout(resolve, ms));
const flags = (params: readonly { readonly flag: string }[]) => params.map((param) => param.flag);

function makeData() {
  let tok = 0;
  return new DataService({
    store: new MemoryBlockStore(),
    scanner: new FakeScanner(),
    inflater: new NodeInflater(),
    checker: new FakeInputChecker(),
    digest: sha256Hex,
    newToken: () => `t-${++tok}`,
    cleanup: Promise.resolve({ state: 'done', removedSessions: 0 }),
  });
}

describe('Masking domain helpers', () => {
  it('names the tools, DustMasker first, with labels for the Results tab and Run details', () => {
    expect(MASKING_TOOLS.map((tool) => tool.id)).toEqual(['dustmasker', 'segmasker', 'convert2blastmask']);
    expect(MASKING_TOOLS.map((tool) => toolLabel(tool.id))).toEqual(['DustMasker', 'SegMasker', 'convert2blastmask']);
    expect(toolStageLabel('dustmasker', 'mask')).toBe('dustmasker');
    expect(toolRunTitle('dustmasker', 'mask')).toBe('DustMasker');
    expect(toolRunTitle('convert2blastmask', 'mask')).toBe('convert2blastmask');
    // The other tools keep their titles.
    expect(toolRunTitle('windowmasker', 'mask')).toBe('WindowMasker · Mask');
    expect(toolRunTitle('blast_formatter', 'format')).toBe('Format Archive (blast_formatter)');
    expect(toolStageLabel('blast_formatter', 'format')).toBe('blast_formatter');
    expect(isMaskingTool('segmasker')).toBe(true);
    expect(isMaskingTool('windowmasker')).toBe(false);
    expect(maskingToolSummary('dustmasker', DUST)).toBe(DUST.description);
    expect(maskingToolSummary('segmasker', undefined)).toBe('Low complexity region masker based on the SEG algorithm');
  });

  it("takes every option from -xmlhelp in the tool's usage order: blocks, required options and the options section", () => {
    expect(flags(allMaskingParameters(DUST, 'dustmasker'))).toEqual([
      '-in',
      '-out',
      '-window',
      '-level',
      '-linker',
      '-infmt',
      '-outfmt',
      '-parse_seqids',
      '-hard_masking',
    ]);
    expect(flags(maskingOptionParameters(DUST, 'dustmasker'))).toEqual(['-window', '-level', '-linker', '-parse_seqids', '-hard_masking']);
    expect(flags(maskingOptionParameters(SEG, 'segmasker'))).toEqual(['-parse_seqids', '-window', '-locut', '-hicut']);
    expect(flags(maskingOptionParameters(C2B, 'convert2blastmask'))).toEqual(['-parse_seqids']);
    expect(flags(maskingRequiredParameters(C2B, 'convert2blastmask'))).toEqual(['-masking_algorithm', '-masking_options']);
    expect(maskingRequiredParameters(DUST, 'dustmasker')).toEqual([]);
    expect(allMaskingParameters(undefined, 'dustmasker')).toEqual([]);
  });

  it('builds the argv: -in, -out (left out when empty or -), then the options in usage order without defaults', () => {
    expect(buildMaskingArgv({ tool: 'dustmasker', inputName: 'input.fa', outName: '', values: {}, description: DUST })).toEqual([
      'dustmasker',
      '-in',
      'input.fa',
    ]);
    expect(
      buildMaskingArgv({
        tool: 'dustmasker',
        inputName: 'input.fa',
        outName: '-',
        values: { '-window': '64', '-level': '30', '-outfmt': 'fasta', '-hard_masking': true },
        description: DUST,
        extra: ['-logfile', 'run.log'],
      }),
    ).toEqual(['dustmasker', '-in', 'input.fa', '-level', '30', '-outfmt', 'fasta', '-hard_masking', '-logfile', 'run.log']);
    expect(
      buildMaskingArgv({
        tool: 'convert2blastmask',
        inputName: 'lc.fa',
        outName: 'mask.asnt',
        values: { '-masking_algorithm': 'dust', '-masking_options': 'window=64, level=20' },
        description: C2B,
      }),
    ).toEqual(['convert2blastmask', '-in', 'lc.fa', '-out', 'mask.asnt', '-masking_algorithm', 'dust', '-masking_options', 'window=64, level=20']);
    expect(() => buildMaskingArgv({ tool: 'segmasker', values: {}, description: SEG, extra: ['a\0b'] })).toThrow(/Additional options/);
  });

  it('writes a text that the tool writes or stores as given as typed, and takes an empty value from Additional options (review of WW4c, F3)', () => {
    const typed = '  -window 64 -level 20  ';
    const c2b = (values: FormValues, extra?: readonly string[]) =>
      buildMaskingArgv({ tool: 'convert2blastmask', inputName: 'lc.fa', values, description: C2B, ...(extra === undefined ? {} : { extra }) });
    expect(c2b({ '-masking_algorithm': ' Dust', '-masking_options': typed })).toEqual([
      'convert2blastmask', '-in', 'lc.fa', '-masking_algorithm', ' Dust', '-masking_options', typed,
    ]);
    // An empty value, as on the command line (`-masking_options ''`), from Additional options; an empty box writes nothing.
    expect(c2b({ '-masking_algorithm': 'dust', '-masking_options': '' }, ['-masking_options', ''])).toEqual([
      'convert2blastmask', '-in', 'lc.fa', '-masking_algorithm', 'dust', '-masking_options', '',
    ]);
    // The fields that change, tool by tool: a value with white space around it is written as typed in these, and
    // without the white space in every other text box (names, numbers, keywords), as before and as the search forms.
    const descriptions = { windowmasker: WM, blast_formatter: BF, makeblastdb: MKDB, dustmasker: DUST, segmasker: SEG, convert2blastmask: C2B };
    const asTyped: Record<string, string[]> = {};
    for (const [tool, description] of Object.entries(descriptions)) {
      asTyped[tool] = [];
      for (const param of description.parameters) {
        if (toToolParameterField(param).kind !== 'text') continue;
        const values = { [param.flag]: ' 7 ' };
        const argv =
          tool === 'makeblastdb'
            ? buildMakeblastdbArgv({ input: '', values, description })
            : tool === 'windowmasker'
              ? buildWindowMaskerArgv({ stage: 'mask', values, description })
              : tool === 'blast_formatter'
                ? buildBlastFormatterArgv({ values, description })
                : buildMaskingArgv({ tool: tool as 'dustmasker', values, description });
        if (argv.includes(' 7 ')) asTyped[tool]!.push(param.flag);
        else if (argv.includes(param.flag)) expect(argv[argv.indexOf(param.flag) + 1], `${tool} ${param.flag}`).toBe('7');
      }
    }
    expect(asTyped).toEqual({
      windowmasker: ['-meta'],
      blast_formatter: [],
      makeblastdb: ['-mask_desc', '-mask_id', '-title'],
      dustmasker: [],
      segmasker: [],
      convert2blastmask: ['-masking_algorithm', '-masking_options'],
    });
    expect([...TOOL_WRITTEN_TEXT_FLAGS].sort()).toEqual(Object.values(asTyped).flat().sort());
  });

  it('counts the changed options of the options section and Additional options', () => {
    expect(changedMaskingOptionsCount('dustmasker', {}, '', DUST)).toBe(0);
    expect(changedMaskingOptionsCount('dustmasker', { '-window': '64', '-level': '30', '-outfmt': 'fasta' }, '-logfile a', DUST)).toBe(2);
    expect(changedMaskingOptionsCount('convert2blastmask', { '-masking_algorithm': 'dust' }, '', C2B)).toBe(0);
  });

  it('records the -outfmt default of the -xmlhelp of a run that masks, and of no other run', () => {
    expect(maskOutfmtDefault('dustmasker', 'mask', DUST)).toBe('interval');
    expect(maskOutfmtDefault('segmasker', 'mask', SEG)).toBe('interval');
    // convert2blastmask's default -outfmt is maskinfo_asn1_text.
    expect(maskOutfmtDefault('convert2blastmask', 'mask', C2B)).toBe('maskinfo_asn1_text');
    expect(maskOutfmtDefault('convert2blastmask', 'mask', undefined)).toBeUndefined();
    // windowmasker's Mask stage writes masking information too; its counts stages and the other tools never do.
    expect(maskOutfmtDefault('windowmasker', 'mask', WM)).toBe('interval');
    expect(maskOutfmtDefault('windowmasker', 'mk_counts', WM)).toBeUndefined();
    expect(maskOutfmtDefault('blast_formatter', 'format', undefined)).toBeUndefined();
    expect(maskOutfmtDefault('makeblastdb', 'build', undefined)).toBeUndefined();
  });

  it('names the mask data of a masking run by one rule: the effective -outfmt (the argv, else the recorded default) is a maskinfo format', () => {
    const run = (tool: 'dustmasker' | 'segmasker' | 'convert2blastmask' | 'windowmasker', argv: readonly string[], outfmtDefault?: string) =>
      ({ tool, stage: 'mask', argv, number: 4, ...(outfmtDefault === undefined ? {} : { outfmtDefault }) }) as const;
    const stdout = { key: 'stdout', name: 'stdout' };
    const wrote = (name: string) => ({ key: 'file:0', name });
    const formats = { maskinfo_asn1_bin: 'asnb', maskinfo_asn1_text: 'asnt' } as const;
    for (const tool of ['dustmasker', 'segmasker', 'convert2blastmask'] as const) {
      for (const [format, extension] of Object.entries(formats)) {
        // Standard output, under a name for the run; an -out file, under its own name (and not the standard output).
        expect(maskDataName(run(tool, [tool, '-outfmt', format], 'interval'), stdout), `${tool} ${format}`).toBe(`${tool}-run4.${extension}`);
        expect(maskDataName(run(tool, [tool, `-outfmt=${format}`], 'interval'), stdout), `${tool} ${format}=`).toBe(`${tool}-run4.${extension}`);
        const toFile = run(tool, [tool, '-out', `./m.${extension}`, '-outfmt', format], 'interval');
        expect(maskDataName(toFile, wrote(`m.${extension}`))).toBe(`m.${extension}`);
        expect(maskDataName(toFile, stdout)).toBeUndefined();
        expect(maskDataName(toFile, { key: 'file:1', name: 'run.log' })).toBeUndefined();
        expect(maskDataName(toFile, { key: 'stderr', name: 'stderr' })).toBeUndefined();
      }
      // The other formats are not mask data; maskinfo_xml is not read by makeblastdb (`MASK_DATA_FORMATS`).
      for (const format of ['interval', 'fasta', 'acclist', 'seqloc_asn1_bin', 'seqloc_asn1_text', 'seqloc_xml', 'maskinfo_xml']) {
        expect(maskDataName(run(tool, [tool, '-outfmt', format], 'maskinfo_asn1_text'), stdout), `${tool} ${format}`).toBeUndefined();
      }
    }
    // With no -outfmt in the argv the recorded default decides: convert2blastmask's is maskinfo_asn1_text.
    expect(maskDataName(run('convert2blastmask', ['convert2blastmask', '-in', 'a.fa'], 'maskinfo_asn1_text'), stdout)).toBe('convert2blastmask-run4.asnt');
    expect(maskDataName(run('dustmasker', ['dustmasker', '-in', 'a.fa'], 'interval'), stdout)).toBeUndefined();
    expect(maskDataName(run('convert2blastmask', ['convert2blastmask', '-in', 'a.fa']), stdout)).toBeUndefined();
    // windowmasker's Mask stage is unchanged: an -outfmt in the argv, the -out file or the standard output.
    expect(maskDataName(run('windowmasker', ['windowmasker', '-ustat', 'c', '-outfmt', 'maskinfo_asn1_bin'], 'interval'), stdout)).toBe('windowmasker-run4.asnb');
    expect(maskDataName(run('windowmasker', ['windowmasker', '-ustat', 'c', '-outfmt', 'maskinfo_xml'], 'interval'), stdout)).toBeUndefined();
    expect(maskDataName(run('windowmasker', ['windowmasker', '-ustat', 'c'], 'interval'), stdout)).toBeUndefined();
    expect(maskDataName({ tool: 'windowmasker', stage: 'mk_counts', argv: ['windowmasker', '-mk_counts', '-outfmt', 'maskinfo_xml'], number: 1 }, stdout)).toBeUndefined();
    expect(maskDataName({ tool: 'blast_formatter', stage: 'format', argv: ['blast_formatter', '-outfmt', 'maskinfo_xml'], number: 1 }, stdout)).toBeUndefined();
    // An -out that the engine opens as the standard output is the standard output (review of WW4c, F5).
    for (const out of ['/dev/stdout', '/dev/fd/1', '/proc/self/fd/1', '//dev/./stdout']) {
      const toStdout = run('dustmasker', ['dustmasker', '-in', 'a.fa', '-outfmt', 'maskinfo_asn1_bin', '-out', out], 'interval');
      expect(maskDataName(toStdout, stdout), out).toBe('dustmasker-run4.asnb');
      expect(maskDataName(toStdout, wrote(out)), out).toBeUndefined();
    }
  });
});

describe('the standard output names (review of WW4c, F5)', () => {
  it("are the engine's: the names of web/adapter/tests/standard_output_names.tsv, which the engine's test reads too", () => {
    const list = readFileSync(new URL('../../../adapter/tests/standard_output_names.tsv', import.meta.url), 'utf8');
    const rows = list.split('\n').filter((line) => line !== '' && !line.startsWith('#')).map((line) => line.split('\t'));
    expect(rows).toHaveLength(24);
    for (const [name, kind] of rows) expect(isStandardOutputName(name!), name).toBe(kind === 'stdout');
  });
});

/**
 * The spellings of the WW4e probe (2026-11-28: the native `losat` of Web run 46112132724 in a folder holding input.fa and
 * the folder sub/, and the Web's real engine in Chromium with the same files, sub/ picked as a folder): the value of
 * `-out` or `-logfile` and the name of the file that the run made, which the native run made and the Web's run reported
 * (`web` rows: names outside the run's directory, the Web's run only), and for a file outside the run's directory the
 * name under which it goes to makeblastdb's run directory. Before WW4e the app cut only `./`, so every row that needs
 * more than that had no `-out`/`-logfile` label and no hand-over, and a file outside went over under its path, which
 * makeblastdb's form refused ("Invalid file name: "../m.asnb".").
 */
const RUN_FILE_SPELLINGS: ReadonlyArray<readonly [option: '-out' | '-logfile', value: string, file: string, handed?: string]> = [
  ['-out', 'm.asnb', 'm.asnb'],
  ['-out', './m.asnb', 'm.asnb'],
  ['-out', './/m.asnb', 'm.asnb'],
  ['-out', 'sub/../m.asnb', 'm.asnb'],
  ['-out', 'sub//../m.asnb', 'm.asnb'],
  ['-out', './sub/./../m.asnb', 'm.asnb'],
  ['-out', 'sub/m.asnb', 'sub/m.asnb'],
  ['-out', 'sub//m.asnb', 'sub/m.asnb'],
  ['-out', './sub/m.asnb', 'sub/m.asnb'],
  ['-out', '/proc/self/cwd/m.asnb', 'm.asnb'],
  ['-out', '/proc/self/cwd/sub/../m.asnb', 'm.asnb'],
  ['-out', '../m.asnb', '../m.asnb', 'm.asnb'],
  ['-out', 'sub/../../m.asnb', '../m.asnb', 'm.asnb'],
  ['-logfile', 'run.log', 'run.log'],
  ['-logfile', './run.log', 'run.log'],
  ['-logfile', 'sub/../run.log', 'run.log'],
  ['-logfile', 'sub//run.log', 'sub/run.log'],
  ['-logfile', '/proc/self/cwd/run.log', 'run.log'],
  ['-logfile', 'run.log/', 'run.log'],
  ['-logfile', ' run.log', 'run.log'],
  // web
  ['-out', '/tmp/m.asnb', '/tmp/m.asnb', 'm.asnb'],
  ['-out', '/tmp//m.asnb', '/tmp/m.asnb', 'm.asnb'],
  ['-out', '/var/run/../tmp/m.asnb', '/tmp/m.asnb', 'm.asnb'],
  ['-out', '/var/lock/m.asnb', '/var/lock/m.asnb', 'm.asnb'],
  ['-out', '../nosub/../m.asnb', '../m.asnb', 'm.asnb'],
  ['-out', '/var/run/../tmp/x.out', '/tmp/x.out', 'x.out'],
  ['-logfile', '/var/run/../tmp/x.out', '/var/tmp/x.out'],
  ['-logfile', '/var/run/../tmp/run.log', '/var/tmp/run.log'],
];

describe("the run's files that -out and -logfile name (WW4e)", () => {
  const maskRun = (out: string) =>
    ({ tool: 'dustmasker', stage: 'mask', argv: ['dustmasker', '-in', 'input.fa', '-out', out, '-outfmt', 'maskinfo_asn1_bin'], number: 2 }) as const;

  const rows = RUN_FILE_SPELLINGS.map(([option, value, file, handed]) => ({ option, value, file, handed }));
  it.each(rows)("reads $option '$value' as the run's file $file: its label, its hand-over as mask data, its counts", ({ option, value, file, handed }) => {
    expect(toolFileOptionLabel(['dustmasker', '-in', 'input.fa', option, value], file)).toBe(option);
    if (option !== '-out') return;
    expect(maskDataName(maskRun(value), { key: 'file:0', name: file })).toBe(handed ?? file);
    expect(isCountsOutputFile('mk_counts', ['windowmasker', '-mk_counts', '-in', 'a.fa', '-out', value], file)).toBe(true);
  });

  it('reads both options in one run: one file of two spellings, and two files where -out and -logfile read one name apart', () => {
    expect(toolFileOptionLabel(['dustmasker', '-out', 'm.asnb', '-logfile', 'sub/../m.asnb'], 'm.asnb')).toBe('-out, -logfile');
    const both = ['dustmasker', '-in', 'input.fa', '-out', 'sub/../m.asnb', '-outfmt', 'maskinfo_asn1_bin', '-logfile', './m.asnb'];
    expect(toolFileOptionLabel(both, 'm.asnb')).toBe('-out, -logfile');
    expect(maskDataName({ ...maskRun('sub/../m.asnb'), argv: both }, { key: 'file:0', name: 'm.asnb' })).toBe('m.asnb');
    // The system walks `-out` past the link /var/run (to /run, then up to /); NCBI reads `-logfile` by its words.
    const apart = ['dustmasker', '-out', '/var/run/../tmp/x.out', '-logfile', '/var/run/../tmp/x.out'];
    expect(toolFileOptionLabel(apart, '/tmp/x.out')).toBe('-out');
    expect(toolFileOptionLabel(apart, '/var/tmp/x.out')).toBe('-logfile');
  });

  it('does not read a name as another file, a directory or a stream', () => {
    const other: ReadonlyArray<readonly ['-out' | '-logfile', string, string]> = [
      ['-out', '../m.asnb', 'm.asnb'],
      ['-out', 'm.asnb', '../m.asnb'],
      ['-out', 'sub/m.asnb', 'm.asnb'],
      ['-out', 'sub/../m.asnb', 'sub/m.asnb'],
      ['-out', 'm.asnb/', 'm.asnb'],
      ['-out', 'm.asnb/.', 'm.asnb'],
      ['-out', '-', '-'],
      ['-out', '/tmp/m.asnb', 'm.asnb'],
      ['-out', '/dev/stdout', 'stdout'],
      ['-logfile', '-', '-'],
      ['-logfile', '/dev/null', '/dev/null'],
      ['-logfile', '', 'input.fa'],
      ['-logfile', '../run.log', 'run.log'],
    ];
    for (const [option, value, file] of other) {
      expect(toolFileOptionLabel(['dustmasker', option, value], file), `${option} ${value}: ${file}`).toBeUndefined();
    }
    expect(maskDataName(maskRun('sub/m.asnb'), { key: 'file:0', name: 'm.asnb' })).toBeUndefined();
    expect(isCountsOutputFile('mk_counts', ['windowmasker', '-mk_counts', '-out', 'a', '-logfile', 'sub/../log'], 'log')).toBe(false);
  });

  it("hands a file outside the run's directory to makeblastdb's form under its name alone, which the form takes", () => {
    for (const [out, file] of [
      ['../m.asnb', '../m.asnb'],
      ['sub/../../m.asnb', '../m.asnb'],
      ['/tmp//m.asnb', '/tmp/m.asnb'],
    ] as const) {
      const name = maskDataName(maskRun(out), { key: 'file:0', name: file })!;
      const makeblastdb = new MakeblastdbDraft({
        engine: new FakeEngine(),
        enqueueToolAll: async (): Promise<EnqueueAllResult> => ({ ok: true, runIds: [] }),
        debounceMs: 0,
      });
      makeblastdb.addMaskData(new File(['mask'], 'm.asnb'), name, { kind: 'made', runNumber: 2 });
      const state = makeblastdb.state.get();
      expect(state.message, `-out ${out}`).toBeUndefined();
      expect(state.values['-mask_data'], `-out ${out}`).toBe('m.asnb');
      expect(state.runFiles.map((entry) => entry.name)).toEqual(['m.asnb']);
    }
  });
});

/**
 * The 33 rows of the WW4e probe (ww4e/probe/web-after.jsonl, 2026-11-28: dustmasker through the Masking tab with the
 * reactors of Web run 46112132724, in a run directory of input.fa and the folder sub/ holding sub/x.txt): `-out` and
 * `-logfile` (undefined: not given), and the names of the files that the engine's run made for each (null: none), with
 * its exit status. The 26 rows that the native run also ran made the same files.
 */
const ENGINE_RUN_FILES: ReadonlyArray<readonly [out: string | undefined, log: string | undefined, outFile: string | null, logFile: string | null, status: number]> = [
  ['m.asnb', undefined, 'm.asnb', null, 0],
  ['./m.asnb', undefined, 'm.asnb', null, 0],
  ['.//m.asnb', undefined, 'm.asnb', null, 0],
  ['sub/../m.asnb', undefined, 'm.asnb', null, 0],
  ['sub//../m.asnb', undefined, 'm.asnb', null, 0],
  ['./sub/./../m.asnb', undefined, 'm.asnb', null, 0],
  ['sub/m.asnb', undefined, 'sub/m.asnb', null, 0],
  ['sub//m.asnb', undefined, 'sub/m.asnb', null, 0],
  ['./sub/m.asnb', undefined, 'sub/m.asnb', null, 0],
  ['nosub/../m.asnb', undefined, null, null, 1],
  ['input.fa/../m.asnb', undefined, null, null, 1],
  ['/proc/self/cwd/m.asnb', undefined, 'm.asnb', null, 0],
  ['/proc/self/cwd/sub/../m.asnb', undefined, 'm.asnb', null, 0],
  ['../m.asnb', undefined, '../m.asnb', null, 0],
  ['sub/../../m.asnb', undefined, '../m.asnb', null, 0],
  ['m.asnb/', undefined, null, null, 1],
  [undefined, 'run.log', null, 'run.log', 0],
  [undefined, './run.log', null, 'run.log', 0],
  [undefined, 'sub/../run.log', null, 'run.log', 0],
  [undefined, 'nosub/../run.log', null, null, 0],
  [undefined, 'sub//run.log', null, 'sub/run.log', 0],
  [undefined, '/proc/self/cwd/run.log', null, 'run.log', 0],
  ['m.asnb', 'sub/../m.asnb', 'm.asnb', 'm.asnb', 0],
  ['sub/../m.asnb', './m.asnb', 'm.asnb', 'm.asnb', 0],
  [undefined, 'run.log/', null, 'run.log', 0],
  [undefined, ' run.log', null, 'run.log', 0],
  ['/tmp/m.asnb', undefined, '/tmp/m.asnb', null, 0],
  ['/tmp//m.asnb', undefined, '/tmp/m.asnb', null, 0],
  ['/var/run/../tmp/m.asnb', undefined, '/tmp/m.asnb', null, 0],
  ['/var/lock/m.asnb', undefined, '/var/lock/m.asnb', null, 0],
  ['../nosub/../m.asnb', undefined, '../m.asnb', null, 0],
  ['/var/run/../tmp/x.out', '/var/run/../tmp/x.out', '/tmp/x.out', '/var/tmp/x.out', 0],
  [undefined, '/var/run/../tmp/run.log', null, '/var/tmp/run.log', 0],
];

/** The run's directory of the probes (WW4e, WW4h): the pasted FASTA and the folder sub/ with one file. */
const PROBE_DIRECTORY = ['input.fa', 'sub/x.txt'] as const;
const probeDirectory = (names: readonly string[] = PROBE_DIRECTORY) => names.map((name) => ({ name, file: new Blob(['x\n']) }));

describe("the name that the engine gives a run's file, in the FakeEngine's runs too (WW4h)", () => {
  const rows = ENGINE_RUN_FILES.map(([out, log, outFile, logFile, status]) => ({ out, log, outFile, logFile, status }));

  it.each(rows)("-out '$out' -logfile '$log': the engine's names $outFile and $logFile", ({ out, log, outFile, logFile }) => {
    if (out !== undefined) expect(runFileName('-out', out, PROBE_DIRECTORY)).toBe(outFile ?? undefined);
    if (log !== undefined) expect(runFileName('-logfile', log, PROBE_DIRECTORY)).toBe(logFile ?? undefined);
  });

  it.each(rows)("the FakeEngine's dustmasker -out '$out' -logfile '$log': the engine's files and exit status $status", async ({ out, log, outFile, logFile, status }) => {
    const argv = ['dustmasker', '-in', 'input.fa', ...(out === undefined ? [] : ['-out', out]), '-outfmt', 'maskinfo_asn1_bin', ...(log === undefined ? [] : ['-logfile', log])];
    const run = await new FakeEngine().runTool(argv, { directory: probeDirectory() });
    expect(run.status).toBe(status);
    expect(run.files.map((file) => file.name).sort()).toEqual([...new Set([outFile, logFile].filter((name) => name !== null))].sort());
    if (status !== 0) expect(await run.stderr.text()).toContain(`Error:  Argument "out". File is not accessible:  \`${out}'`);
  });

  it('reads the other names of the WW4h probe as the engine: folders, directories, streams and the special names', () => {
    // The WW4h probe (ww4h/probe/engine.jsonl, 2026-11-28: each tool on the reactors of Web run 46124737526 in Node, in
    // a run directory of input.fa, prot.fa, bn.asn and sub/x.txt): null where the run made no file.
    const names: ReadonlyArray<readonly ['-out' | '-logfile', string, string | null]> = [
      ['-out', 'sub', null],
      ['-out', '/proc/self/cwd', null],
      ['-out', 'nosub/m.out', null],
      ['-out', '-', null],
      ['-out', '/dev/stdout', null],
      ['-logfile', 'sub', null],
      ['-logfile', 'sub/', null],
      ['-logfile', ' sub/../run.log', null],
      ['-logfile', 'sub/../nosub/../run.log', null],
      ['-logfile', 'nosub/run.log', null],
      ['-logfile', '../run.log', '../run.log'],
      ['-logfile', '-', null],
      ['-logfile', '', null],
      ['-logfile', '/dev/null', null],
    ];
    for (const [option, value, file] of names) expect(runFileName(option, value, PROBE_DIRECTORY), `${option} '${value}'`).toBe(file ?? undefined);
    // Without the run's directory no folder is looked up.
    expect(runFileName('-out', 'nosub/../m.asnb')).toBe('m.asnb');
    expect(runFileName('-logfile', 'nosub/../run.log')).toBe('run.log');
    expect(runFileName('-out', 'm.asnb/')).toBeUndefined();
    // A file's folders are the run's folders, at any depth; a file is none.
    expect(runFileName('-out', 'a/b/../m.out', ['a/b/c.fa'])).toBe('a/m.out');
    expect(runFileName('-out', 'a/b/c.fa/../m.out', ['a/b/c.fa'])).toBeUndefined();
  });

  it("fails a run whose -out opens no file with each tool's words, and names the folder's file", async () => {
    // The engine's runs of `-out nosub/../m.out` (WW4h probe): the exit status and the last line of the standard error.
    const words: ReadonlyArray<readonly [argv: readonly string[], status: number, last: string | undefined]> = [
      [['dustmasker', '-in', 'input.fa'], 1, 'Error:  Argument "out". File is not accessible:  `nosub/../m.out\''],
      [['segmasker', '-in', 'prot.fa'], 1, 'Argument "out". File is not accessible:  `nosub/../m.out\''],
      [['convert2blastmask', '-in', 'input.fa', '-masking_algorithm', 'dust', '-masking_options', 'window=64'], 1, 'Argument "out". File is not accessible:  `nosub/../m.out\''],
      [['windowmasker', '-mk_counts', '-in', 'input.fa'], 0, undefined],
      [['windowmasker', '-convert', '-in', 'c.counts'], 0, undefined],
      [['windowmasker', '-ustat', 'c.counts', '-in', 'input.fa'], 1, 'Error:  Argument "out". File is not accessible:  `nosub/../m.out\''],
      [['blast_formatter', '-archive', 'bn.asn'], 1, 'Command line argument error: Argument "out". File is not accessible:  `nosub/../m.out\''],
    ];
    const directory = probeDirectory(['input.fa', 'prot.fa', 'bn.asn', 'c.counts', 'sub/x.txt']);
    for (const [argv, status, last] of words) {
      const failed = await new FakeEngine().runTool([...argv, '-out', 'nosub/../m.out'], { directory });
      expect(failed.status, argv.join(' ')).toBe(status);
      expect(failed.files, argv.join(' ')).toEqual([]);
      const lines = (await failed.stderr.text()).split('\n').filter((line) => line !== '');
      if (last !== undefined) expect(lines[lines.length - 1], argv.join(' ')).toBe(last);
      const made = await new FakeEngine().runTool([...argv, '-out', 'sub/../m.out'], { directory });
      expect(made.status, argv.join(' ')).toBe(0);
      expect(made.files.map((file) => file.name), argv.join(' ')).toEqual(['m.out']);
    }
    // A run given no directory has an empty one, as the engine's: `sub/m.out` has no folder.
    expect((await new FakeEngine().runTool(['segmasker', '-in', 'prot.fa', '-out', 'sub/m.out'])).status).toBe(1);
    // The standard output's names write there, and `-logfile -` there too: no file.
    const streams = await new FakeEngine().runTool(['dustmasker', '-in', 'input.fa', '-out', '/dev/stdout', '-logfile', '-'], { directory });
    expect(streams.files).toEqual([]);
    expect(await streams.stdout.text()).toContain(FAKE_MARKER);
    expect(await streams.stderr.text()).toContain(FAKE_MARKER);
  });

  it("names makeblastdb's database and log as the engine, and fails a -logfile that makeblastdb cannot open", async () => {
    const directory = probeDirectory();
    const run = (words: readonly string[]) => new FakeEngine().runTool(['makeblastdb', '-in', 'input.fa', '-dbtype', 'nucl', ...words], { directory });
    // The engine's database files (WW4h probe): by their path in the run's directory, folders made where missing.
    for (const [out, base] of [
      ['db', 'db'],
      ['./db', 'db'],
      ['db/', 'db'],
      ['sub/../db', 'db'],
      ['sub/db', 'sub/db'],
      ['nosub/db', 'nosub/db'],
      ['sub', 'sub'],
    ] as const) {
      const made = await run(['-out', out]);
      expect(made.status, out).toBe(0);
      expect(made.files.map((file) => file.name), out).toEqual([`${base}.nhr`, `${base}.nin`, `${base}.nsq`]);
    }
    const named = await new FakeEngine().runTool(['makeblastdb', '-in', 'sub/../input.fa', '-dbtype', 'nucl'], { directory });
    expect(named.files.map((file) => file.name)).toEqual(['input.fa.nhr', 'input.fa.nin', 'input.fa.nsq']);
    // -logfile: makeblastdb opens it as NCBI's log and as an output; where the output does not open, exit status 1.
    const db = ['db.nhr', 'db.nin', 'db.nsq'];
    const error = (name: string) => `Command line argument error: Argument "logfile". File is not accessible:  \`${name}'\n`;
    for (const [log, status, files, stderr] of [
      ['sub/../run.log', 0, [...db, 'run.log'], ''],
      ['../run.log', 0, [...db, '../run.log'], ''],
      ['/dev/null', 0, db, ''],
      ['nosub/../run.log', 1, [], error('nosub/../run.log')],
      ['nosub/run.log', 1, [], error('nosub/run.log')],
      ['sub', 1, [], error('sub')],
      ['sub/', 1, [], error('sub/')],
      [' sub/../run.log', 1, [], error(' sub/../run.log')],
      ['run.log/', 1, ['run.log'], ''],
      ['', 1, [], ''],
    ] as const) {
      const made = await run(['-out', 'db', '-logfile', log]);
      expect(made.status, `-logfile '${log}'`).toBe(status);
      expect(made.files.map((file) => file.name), `-logfile '${log}'`).toEqual(files);
      expect(await made.stderr.text(), `-logfile '${log}'`).toBe(stderr);
    }
  });
});

describe('MaskingDraft', () => {
  function setup(options: { readonly validate?: 'fake' | 'ok' } = {}) {
    const fake = new FakeEngine();
    const validateTool = vi.fn(async (...args: Parameters<ToolEngine['validateTool']>) =>
      options.validate === 'ok' ? { ok: true as const } : fake.validateTool(args[0]),
    );
    const describeTool = vi.fn((tool: string) => fake.describeTool(tool));
    const queued: ToolRequest[] = [];
    const enqueueToolAll = vi.fn(async (requests: readonly ToolRequest[]): Promise<EnqueueAllResult> => {
      queued.push(...requests);
      return { ok: true, runIds: requests.map((_, index) => `run-${queued.length - requests.length + index + 1}`) };
    });
    const data = makeData();
    const draft = new MaskingDraft({ engine: { describeTool, validateTool }, data, enqueueToolAll, debounceMs: 0 });
    return { draft, queued, enqueueToolAll, describeTool, validateTool, data };
  }

  async function described(draft: MaskingDraft, tool: string): Promise<void> {
    await vi.waitFor(() => expect(draft.state.get().description?.tool).toBe(tool));
  }

  async function settled(draft: MaskingDraft): Promise<void> {
    await vi.waitFor(() => expect(['ok', 'invalid', 'idle']).toContain(draft.state.get().validation.state));
  }

  it("loads each tool's -xmlhelp when the tool is first shown, and keeps each tool's values, Additional options and -out", async () => {
    const { draft, describeTool } = setup();
    expect(draft.state.get().tool).toBe('dustmasker');
    expect(describeTool).not.toHaveBeenCalled();
    draft.ensureLoaded();
    await described(draft, 'dustmasker');
    expect(describeTool.mock.calls.map((call) => call[0])).toEqual(['dustmasker']);

    draft.setFastaPaste('>q\nACGTACGT\n');
    draft.setField('-window', '32');
    draft.setOutName('dust.txt');
    draft.setExtra('-logfile run.log');
    expect(draft.previewArgv()).toEqual(['dustmasker', '-in', 'input.fa', '-out', 'dust.txt', '-window', '32', '-logfile', 'run.log']);

    draft.setTool('segmasker');
    expect(draft.state.get().description).toBeUndefined();
    await described(draft, 'segmasker');
    // The input is shared; the options, Additional options and -out are the tool's own.
    expect(draft.previewArgv()).toEqual(['segmasker', '-in', 'input.fa']);
    draft.setField('-window', '9');
    expect(draft.previewArgv()).toEqual(['segmasker', '-in', 'input.fa', '-window', '9']);

    draft.setTool('dustmasker');
    expect(draft.state.get().description?.tool).toBe('dustmasker');
    expect(draft.previewArgv()).toEqual(['dustmasker', '-in', 'input.fa', '-out', 'dust.txt', '-window', '32', '-logfile', 'run.log']);
    expect(describeTool.mock.calls.map((call) => call[0])).toEqual(['dustmasker', 'segmasker']);
    expect(draft.previewCommandWords()[0]).toBe('losat');

    // Restore default options clears the options section and Additional options of the tool shown only.
    draft.setField('-outfmt', 'fasta');
    draft.resetOptions();
    expect(draft.previewArgv()).toEqual(['dustmasker', '-in', 'input.fa', '-out', 'dust.txt', '-outfmt', 'fasta']);
    draft.setTool('segmasker');
    expect(draft.previewArgv()).toEqual(['segmasker', '-in', 'input.fa', '-window', '9']);
  });

  it('queues one mask run per input: paste, files (a .gz inflated), Combined or Separate, with -out left out when empty', async () => {
    const { draft, queued } = setup();
    draft.ensureLoaded();
    await described(draft, 'dustmasker');
    expect(draft.holdsInputs()).toBe(false);
    draft.setFastaPaste('>p\nACGT');
    draft.addFastaFiles([new File([gzipSync('>g\nTTTT\n')], 'g.fa.gz')]);
    await vi.waitFor(() => expect(draft.state.get().fastaFiles[0]?.status).toBe('ready'));
    expect(draft.state.get().fastaFiles[0]!.name).toBe('g.fa');
    expect(draft.holdsInputs()).toBe(true);
    expect(draft.runCount()).toBe(1);
    expect(draft.previewArgv()).toEqual(['dustmasker', '-in', 'combined_input.fa']);

    expect((await draft.submit()).ok).toBe(true);
    expect(queued).toHaveLength(1);
    const combined = queued[0]!;
    expect(combined.tool).toBe('dustmasker');
    expect(combined.stage).toBe('mask');
    expect(combined.argv).toEqual(['dustmasker', '-in', 'combined_input.fa']);
    const file = combined.directory.find((entry) => entry.name === 'combined_input.fa')!;
    expect(await file.file.text()).toBe('>p\nACGT\n>g\nTTTT\n');
    expect(file.origin).toEqual({ kind: 'combined', parts: ['input.fa', 'g.fa'] });
    expect(combined.sessionTokens).toHaveLength(1);
    expect(combined.outfmtDefault).toBe('interval');

    draft.setFastaMode('separate');
    expect(draft.runCount()).toBe(2);
    queued.length = 0;
    expect((await draft.submit()).ok).toBe(true);
    expect(queued.map((request) => request.argv)).toEqual([
      ['dustmasker', '-in', 'input.fa'],
      ['dustmasker', '-in', 'g.fa'],
    ]);
    expect(queued[1]!.directory.map((entry) => entry.origin)).toEqual([{ kind: 'inflated', from: 'g.fa.gz' }]);
  });

  it('reads a BLAST database by its -in name with -infmt blastdb, its files added to the run directory', async () => {
    const { draft, queued } = setup();
    draft.ensureLoaded();
    await described(draft, 'dustmasker');
    draft.setField('-infmt', 'blastdb');
    expect(draft.usesFastaBox()).toBe(false);
    const nin = new File(['n'], 'db.nin');
    const nsq = new File(['s'], 'db.nsq');
    Object.defineProperty(nin, 'webkitRelativePath', { value: 'dbs/db.nin' });
    Object.defineProperty(nsq, 'webkitRelativePath', { value: 'dbs/db.nsq' });
    draft.dropInputFiles([nin, nsq]);
    expect(draft.state.get().blastDbIn).toBe('dbs/db');
    expect(draft.previewArgv()).toEqual(['dustmasker', '-in', 'dbs/db', '-infmt', 'blastdb']);
    expect((await draft.submit()).ok).toBe(true);
    expect(queued[0]!.argv).toEqual(['dustmasker', '-in', 'dbs/db', '-infmt', 'blastdb']);
    expect(queued[0]!.directory.map((entry) => entry.name)).toEqual(['dbs/db.nin', 'dbs/db.nsq']);

    draft.setBlastDbIn('');
    const refused = await draft.submit();
    expect(refused).toEqual({ ok: false, message: 'Enter an input name or open files.' });
  });

  it("reports NCBI's error while convert2blastmask's required options are missing, as the command line does", async () => {
    const { draft, queued } = setup();
    draft.setTool('convert2blastmask');
    await described(draft, 'convert2blastmask');
    draft.setFastaPaste('>lc\nACGTacgtACGT\n');
    await vi.waitFor(() => expect(draft.state.get().validation.state).toBe('invalid'));
    const first = draft.state.get().validation;
    expect(first.state === 'invalid' && first.message).toMatch(/Argument "masking_algorithm"\. Mandatory value is missing/);
    expect(await draft.submit()).toMatchObject({ ok: false });
    expect(queued).toHaveLength(0);

    draft.setField('-masking_algorithm', 'dust');
    await vi.waitFor(() => {
      const validation = draft.state.get().validation;
      expect(validation.state === 'invalid' && validation.message).toMatch(/Argument "masking_options"/);
    });
    draft.setField('-masking_options', 'window=64');
    await vi.waitFor(() => expect(draft.state.get().validation.state).toBe('ok'));
    expect((await draft.submit()).ok).toBe(true);
    expect(queued[0]!.argv).toEqual(['convert2blastmask', '-in', 'input.fa', '-masking_algorithm', 'dust', '-masking_options', 'window=64']);
    // Its default -outfmt (maskinfo_asn1_text) is masking information: the run records it, and its standard output is mask data.
    expect(queued[0]!.outfmtDefault).toBe('maskinfo_asn1_text');
    expect(maskDataName({ ...queued[0]!, number: 5 }, { key: 'stdout', name: 'stdout' })).toBe('convert2blastmask-run5.asnt');
  });

  it('records the -xmlhelp default of -outfmt, and the argv decides what the run writes as mask data', async () => {
    const { draft, queued } = setup({ validate: 'ok' });
    draft.ensureLoaded();
    await described(draft, 'dustmasker');
    draft.setFastaPaste('>q\nACGT\n');
    draft.setField('-outfmt', 'maskinfo_asn1_bin');
    draft.setOutName('dust.asnb');
    expect((await draft.submit()).ok).toBe(true);
    expect(queued[0]!.argv).toEqual(['dustmasker', '-in', 'input.fa', '-out', 'dust.asnb', '-outfmt', 'maskinfo_asn1_bin']);
    const file = { key: 'file:0', name: 'dust.asnb' };
    expect(maskDataName({ ...queued[0]!, number: 1 }, file)).toBe('dust.asnb');
    draft.setField('-outfmt', 'interval');
    expect((await draft.submit()).ok).toBe(true);
    expect(queued[1]!.outfmtDefault).toBe('interval');
    expect(maskDataName({ ...queued[1]!, number: 2 }, file)).toBeUndefined();
    draft.setExtra('-outfmt maskinfo_asn1_text');
    expect((await draft.submit()).ok).toBe(true);
    expect(maskDataName({ ...queued[2]!, number: 3 }, file)).toBe('dust.asnb');
  });

  it('runs help and version words of Additional options as print runs, without an input', async () => {
    const { draft, queued } = setup();
    draft.ensureLoaded();
    await described(draft, 'dustmasker');
    draft.setExtra('-version-full');
    expect(draft.previewArgv()).toEqual(['dustmasker', '-version-full']);
    await settled(draft);
    expect((await draft.submit()).ok).toBe(true);
    expect(queued[0]).toMatchObject({ tool: 'dustmasker', stage: 'mask', argv: ['dustmasker', '-version-full'], prints: 'version' });
    expect(queued[0]!.outfmtDefault).toBeUndefined();
    draft.setTool('convert2blastmask');
    draft.setExtra('-help');
    expect((await draft.submit()).ok).toBe(true);
    expect(queued[1]).toMatchObject({ tool: 'convert2blastmask', argv: ['convert2blastmask', '-help'], prints: 'help' });
    expect(queued[1]!.outfmtDefault).toBeUndefined();
  });

  it("checks the argv in a directory that holds the run's input under the run's name (review of WW4c, F2)", async () => {
    const { draft, queued, validateTool } = setup({ validate: 'ok' });
    draft.ensureLoaded();
    await described(draft, 'dustmasker');
    draft.setTool('segmasker');
    await described(draft, 'segmasker');
    const lastCheck = () => {
      const [argv, options] = validateTool.mock.calls.at(-1)!;
      return { argv, directory: options?.directory ?? [] };
    };
    // An -out that names the input: the check's directory holds the input, as the run's does.
    draft.setFastaPaste('>q\nMKLLLL\n');
    draft.setOutName('input.fa');
    await vi.waitFor(() => expect(lastCheck().argv).toEqual(['segmasker', '-in', 'input.fa', '-out', 'input.fa']));
    expect(lastCheck().directory.map((entry) => entry.name)).toEqual(['input.fa']);
    expect(await lastCheck().directory[0]!.file.text()).toBe('>q\nMKLLLL\n');
    expect((await draft.submit()).ok).toBe(true);
    expect(queued.at(-1)!.directory.map((entry) => entry.name)).toEqual(['input.fa']);
    // Several inputs combined: combined_input.fa, in the check as in the run; Separate: the first run's input.
    draft.addFastaFiles([new File(['>a\nMK\n'], 'a.fa'), new File(['>b\nLL'], 'b.fa')]);
    await vi.waitFor(() => expect(draft.state.get().fastaFiles.every((entry) => entry.status === 'ready')).toBe(true));
    draft.setOutName('combined_input.fa');
    await vi.waitFor(() => expect(lastCheck().argv).toEqual(['segmasker', '-in', 'combined_input.fa', '-out', 'combined_input.fa']));
    expect(lastCheck().directory.map((entry) => entry.name)).toEqual(['combined_input.fa']);
    expect((await draft.submit()).ok).toBe(true);
    expect(queued.at(-1)!.directory.map((entry) => entry.name)).toEqual(['combined_input.fa']);
    draft.setFastaMode('separate');
    draft.setOutName('m.out');
    await vi.waitFor(() => expect(lastCheck().argv).toEqual(['segmasker', '-in', 'input.fa', '-out', 'm.out']));
    expect(lastCheck().directory.map((entry) => entry.name)).toEqual(['input.fa']);
    // Help reads no input: none in the check's directory, as in the run's.
    draft.setExtra('-help');
    await vi.waitFor(() => expect(lastCheck().argv).toEqual(['segmasker', '-out', 'm.out', '-help']));
    expect(lastCheck().directory).toEqual([]);
  });

  it("checks WindowMasker's argv with its FASTA input and counts under the run's names (review of WW4c, F2)", async () => {
    const fake = new FakeEngine();
    const validateTool = vi.fn<ToolEngine['validateTool']>(async () => ({ ok: true as const }));
    const queued: ToolRequest[] = [];
    const draft = new ToolDraft({
      engine: { describeTool: (tool: string) => fake.describeTool(tool), validateTool },
      data: makeData(),
      enqueueToolAll: async (requests) => {
        queued.push(...requests);
        return { ok: true, runIds: requests.map((_, index) => `wm-${index}`) };
      },
      debounceMs: 0,
    });
    draft.ensureLoaded();
    await vi.waitFor(() => expect(draft.state.get().description).toBeDefined());
    const names = () => (validateTool.mock.calls.at(-1)![1]?.directory ?? []).map((entry) => entry.name);
    draft.setStage('mask');
    draft.setFastaPaste('>q\nACGT\n');
    draft.setUstatFile(new File(['counts'], 'c.ascii'));
    draft.setOutName('c.ascii');
    await vi.waitFor(() => expect(validateTool.mock.calls.at(-1)![0]).toContain('c.ascii'));
    await vi.waitFor(() => expect(names()).toEqual(['input.fa', 'c.ascii']));
    expect((await draft.submit()).ok).toBe(true);
    expect(queued.at(-1)!.directory.map((entry) => entry.name)).toEqual(['input.fa', 'c.ascii']);
    // A saved counts file stands in by its name.
    draft.setUstatSaved(562, 6);
    await vi.waitFor(() => expect(names()).toEqual(['input.fa', countsName(562)]));
    draft.setStage('mk_counts');
    await vi.waitFor(() => expect(validateTool.mock.calls.at(-1)![0]).toContain('-mk_counts'));
    expect(names()).toEqual(['input.fa']);
    draft.setStage('convert');
    draft.setConvertFile(new File(['counts'], 'u.ascii'));
    await vi.waitFor(() => expect(names()).toEqual(['u.ascii']));
  });

  it('runs -help-full as a help run of every tool: no input, no -outfmt default, no mask data (review of WW4c, F4)', async () => {
    for (const tool of ['windowmasker', 'blast_formatter', 'makeblastdb', 'dustmasker', 'segmasker', 'convert2blastmask']) {
      expect(toolPrintsKind([tool, '-help-full']), tool).toBe('help');
      const printed = await new FakeEngine().runTool([tool, '-help-full']);
      expect(printed.status, tool).toBe(0);
      expect(await printed.stdout.text(), tool).toMatch(/^USAGE\n/);
    }
    const { draft, queued } = setup();
    draft.ensureLoaded();
    await described(draft, 'dustmasker');
    draft.setTool('convert2blastmask');
    await described(draft, 'convert2blastmask');
    draft.setFastaPaste('>q\nACGTacgt\n');
    draft.setExtra('-help-full');
    expect(draft.previewArgv()).toEqual(['convert2blastmask', '-help-full']);
    await settled(draft);
    expect((await draft.submit()).ok).toBe(true);
    expect(queued[0]).toMatchObject({ tool: 'convert2blastmask', argv: ['convert2blastmask', '-help-full'], prints: 'help' });
    expect(queued[0]!.outfmtDefault).toBeUndefined();
    expect(maskDataName({ ...queued[0]!, number: 1 }, { key: 'stdout', name: 'stdout' })).toBeUndefined();
  });

  it('holds its inflated session files apart from another draft that shares the holder tracker', async () => {
    const data = makeData();
    const dropped: string[] = [];
    const sessionFiles = new SessionFiles((token) => dropped.push(token));
    const fake = new FakeEngine();
    const deps = {
      engine: fake,
      data,
      enqueueToolAll: async (): Promise<EnqueueAllResult> => ({ ok: true, runIds: [] }),
      sessionFiles,
      debounceMs: 0,
    };
    const first = new MaskingDraft(deps);
    const second = new ToolDraft(deps);
    first.addFastaFiles([new File([gzipSync('>a\nAC\n')], 'a.fa.gz')]);
    second.addFastaFiles([new File([gzipSync('>b\nGT\n')], 'b.fa.gz')]);
    await vi.waitFor(() => expect(first.state.get().fastaFiles[0]?.status).toBe('ready'));
    await vi.waitFor(() => expect(second.state.get().fastaFiles[0]?.status).toBe('ready'));
    // Both drafts' first file has the key fa-1: removing one must not drop the other's session file.
    expect(first.state.get().fastaFiles[0]!.key).toBe(second.state.get().fastaFiles[0]!.key);
    const firstToken = first.state.get().fastaFiles[0]!.sessionToken!;
    const secondToken = second.state.get().fastaFiles[0]!.sessionToken!;
    expect(firstToken).not.toBe(secondToken);
    first.removeFastaFile(first.state.get().fastaFiles[0]!.key);
    expect(dropped).toEqual([firstToken]);
    second.removeFastaFile(second.state.get().fastaFiles[0]!.key);
    expect(dropped).toEqual([firstToken, secondToken]);
  });
});

describe('Coordinator hands over the output that holds masking information', () => {
  async function runMasking(argv: readonly string[], outfmtDefault?: string) {
    const data = makeData();
    let id = 0;
    const coordinator = new Coordinator({
      engine: new FakeEngine(),
      data,
      downloader: { save: () => undefined },
      now: () => 10_000,
      newRunId: () => `ms-${++id}`,
      digest: sha256Hex,
    });
    const result = await coordinator.enqueueToolAll([
      {
        tool: argv[0] as 'convert2blastmask',
        stage: 'mask',
        argv,
        directory: [{ name: 'input.fa', file: new Blob(['>q\nacgtACGT\n']), origin: { kind: 'pasted' } }],
        ...(outfmtDefault === undefined ? {} : { outfmtDefault }),
      },
    ]);
    expect(result.ok).toBe(true);
    const runId = result.ok ? result.runIds![0]! : '';
    for (let i = 0; i < 400; i++) {
      const run = coordinator.state.get().runs.find((entry) => entry.snapshot.runId === runId);
      if (run !== undefined && run.status === 'completed' && isToolRun(run)) return { run, coordinator };
      await sleep(5);
    }
    throw new Error('the run did not complete');
  }

  it('keeps the recorded default on the run, and names standard output and the -out file as maskDataName does', async () => {
    const stdout = await runMasking(['convert2blastmask', '-in', 'input.fa', '-masking_algorithm', 'dust', '-masking_options', 'o'], 'maskinfo_asn1_text');
    expect(stdout.run.snapshot.outfmtDefault).toBe('maskinfo_asn1_text');
    expect(stdout.run.snapshot.stage).toBe('mask');
    const handed = await stdout.coordinator.maskData(stdout.run.snapshot.runId, 'stdout');
    expect(handed.name).toBe(`convert2blastmask-run${stdout.run.snapshot.number}.asnt`);
    expect(handed.origin).toEqual({ kind: 'made', runNumber: stdout.run.snapshot.number });
    expect(await handed.file.text()).toContain(FAKE_MARKER);
    await expect(stdout.coordinator.maskData(stdout.run.snapshot.runId, 'stderr')).rejects.toThrow('the output is not mask data');

    const file = await runMasking(['dustmasker', '-in', 'input.fa', '-out', 'm.asnt', '-outfmt', 'maskinfo_asn1_text', '-logfile', 'run.log'], 'interval');
    expect(file.run.record.stdout).toBeDefined();
    await expect(file.coordinator.maskData(file.run.snapshot.runId, 'stdout')).rejects.toThrow('the output is not mask data');
    const made = file.run.record.files!.find((entry) => entry.name === 'm.asnt')!;
    expect((await file.coordinator.maskData(file.run.snapshot.runId, made.key)).name).toBe('m.asnt');
    const log = file.run.record.files!.find((entry) => entry.name === 'run.log')!;
    await expect(file.coordinator.maskData(file.run.snapshot.runId, log.key)).rejects.toThrow('the output is not mask data');

    const none = await runMasking(['dustmasker', '-in', 'input.fa'], 'interval');
    await expect(none.coordinator.maskData(none.run.snapshot.runId, 'stdout')).rejects.toThrow('the output is not mask data');
  });

  it("hands a masking run's file to the makeblastdb form's -mask_data", async () => {
    const { run, coordinator } = await runMasking(['segmasker', '-in', 'input.fa', '-out', 'seg.asnb', '-outfmt', 'maskinfo_asn1_bin'], 'interval');
    const key = run.record.files!.find((entry) => entry.name === 'seg.asnb')!.key;
    const handed = await coordinator.maskData(run.snapshot.runId, key);
    const makeblastdb = new MakeblastdbDraft({
      engine: new FakeEngine(),
      enqueueToolAll: async (): Promise<EnqueueAllResult> => ({ ok: true, runIds: [] }),
      debounceMs: 0,
    });
    makeblastdb.addMaskData(handed.file, handed.name, handed.origin);
    const state = makeblastdb.state.get();
    expect(state.values['-mask_data']).toBe('seg.asnb');
    expect(state.runFiles.map((entry) => [entry.name, entry.flag, entry.origin])).toEqual([
      ['seg.asnb', '-mask_data', { kind: 'made', runNumber: run.snapshot.number }],
    ]);
  });
});

describe('FakeEngine rows of the Masking tools', () => {
  it('describes the three tools from their -xmlhelp and refuses as their command lines do', async () => {
    const fake = new FakeEngine();
    for (const tool of ['dustmasker', 'segmasker', 'convert2blastmask']) {
      const description = await fake.describeTool(tool);
      expect(description.tool).toBe(tool);
      expect(await fake.describeTool(tool)).toBe(description);
    }
    await expect(fake.describeTool('no_such_tool')).rejects.toThrow('unknown tool');
    expect(await fake.validateTool(['dustmasker', '-in', 'a.fa', '-window', '32'])).toEqual({ ok: true });
    expect(await fake.validateTool(['segmasker', '-outfmt', 'acclist'])).toMatchObject({ ok: false, message: expect.stringMatching(/Invalid value/) });
    expect(await fake.validateTool(['dustmasker', '-masking_algorithm', 'x'])).toMatchObject({ ok: false, message: expect.stringMatching(/Unknown argument/) });
    expect(await fake.validateTool(['convert2blastmask', '-in', 'a.fa', '-masking_options', 'o'])).toEqual({
      ok: false,
      message: "Error: (CArgException::eNoArg) Argument \"masking_algorithm\". Mandatory value is missing:  `String'",
    });
    expect(await fake.validateTool(['convert2blastmask', '-version'])).toEqual({ ok: true });

    const help = await fake.runTool(['segmasker', '-xmlhelp']);
    expect(help.status).toBe(0);
    expect(parseToolDescription(await help.stdout.text()).tool).toBe('segmasker');
    const run = await fake.runTool(['dustmasker', '-in', 'a.fa', '-outfmt', 'maskinfo_xml']);
    expect(run.status).toBe(0);
    expect(await run.stdout.text()).toContain(FAKE_MARKER);
  });
});
