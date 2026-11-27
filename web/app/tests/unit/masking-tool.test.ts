// Unit tests for the Masking tool tab (WW4c): the domain helpers of dustmasker, segmasker and
// convert2blastmask (options from -xmlhelp, argv builder, the masking-information rule), MaskingDraft
// (argv per tool, per-tool values, -out, the shared FASTA input, blastdb, required options, help and
// version, the `-outfmt` default that a run records), the one rule for mask data (`maskDataName`) with the
// Coordinator's and the makeblastdb draft's hand-over, and the FakeEngine's rows for the three tools.
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
import {
  buildMaskingArgv,
  changedMaskingOptionsCount,
  maskDataName,
  maskOutfmtDefault,
} from '../../src/domain/tool-argv';
import {
  allMaskingParameters,
  isMaskingTool,
  MASKING_TOOLS,
  maskingOptionParameters,
  maskingRequiredParameters,
  maskingToolSummary,
  toolLabel,
  toolRunTitle,
  toolStageLabel,
} from '../../src/domain/tools';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FAKE_MARKER, FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import CONVERT2BLASTMASK_XMLHELP_DOC from '../../src/infra/fake/convert2blastmask-xmlhelp.json';
import DUSTMASKER_XMLHELP_DOC from '../../src/infra/fake/dustmasker-xmlhelp.json';
import SEGMASKER_XMLHELP_DOC from '../../src/infra/fake/segmasker-xmlhelp.json';
import WINDOWMASKER_XMLHELP_DOC from '../../src/infra/fake/windowmasker-xmlhelp.json';
import { NodeInflater } from './support/node-inflater';

const DUST = parseToolDescription(DUSTMASKER_XMLHELP_DOC.xml);
const SEG = parseToolDescription(SEGMASKER_XMLHELP_DOC.xml);
const C2B = parseToolDescription(CONVERT2BLASTMASK_XMLHELP_DOC.xml);
const WM = parseToolDescription(WINDOWMASKER_XMLHELP_DOC.xml);
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
    expect(() => buildMaskingArgv({ tool: 'segmasker', values: {}, description: SEG, extra: [''] })).toThrow(/Additional options/);
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
    const formats = { maskinfo_asn1_bin: 'asnb', maskinfo_asn1_text: 'asnt', maskinfo_xml: 'xml' } as const;
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
      // The other formats are not mask data.
      for (const format of ['interval', 'fasta', 'acclist', 'seqloc_asn1_bin', 'seqloc_asn1_text', 'seqloc_xml']) {
        expect(maskDataName(run(tool, [tool, '-outfmt', format], 'maskinfo_asn1_text'), stdout), `${tool} ${format}`).toBeUndefined();
      }
    }
    // With no -outfmt in the argv the recorded default decides: convert2blastmask's is maskinfo_asn1_text.
    expect(maskDataName(run('convert2blastmask', ['convert2blastmask', '-in', 'a.fa'], 'maskinfo_asn1_text'), stdout)).toBe('convert2blastmask-run4.asnt');
    expect(maskDataName(run('dustmasker', ['dustmasker', '-in', 'a.fa'], 'interval'), stdout)).toBeUndefined();
    expect(maskDataName(run('convert2blastmask', ['convert2blastmask', '-in', 'a.fa']), stdout)).toBeUndefined();
    // windowmasker's Mask stage is unchanged: an -outfmt in the argv, the -out file or the standard output.
    expect(maskDataName(run('windowmasker', ['windowmasker', '-ustat', 'c', '-outfmt', 'maskinfo_xml'], 'interval'), stdout)).toBe('windowmasker-run4.xml');
    expect(maskDataName(run('windowmasker', ['windowmasker', '-ustat', 'c'], 'interval'), stdout)).toBeUndefined();
    expect(maskDataName({ tool: 'windowmasker', stage: 'mk_counts', argv: ['windowmasker', '-mk_counts', '-outfmt', 'maskinfo_xml'], number: 1 }, stdout)).toBeUndefined();
    expect(maskDataName({ tool: 'blast_formatter', stage: 'format', argv: ['blast_formatter', '-outfmt', 'maskinfo_xml'], number: 1 }, stdout)).toBeUndefined();
    // A snapshot that went through storage names its output as the live one does.
    const live = run('convert2blastmask', ['convert2blastmask', '-in', 'a.fa'], 'maskinfo_asn1_text');
    expect(maskDataName(JSON.parse(JSON.stringify(live)) as typeof live, stdout)).toBe(maskDataName(live, stdout));
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
    draft.setExtra('-outfmt maskinfo_xml');
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
