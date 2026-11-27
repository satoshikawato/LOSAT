// Unit tests for the Masking tool tab (WW4c): the domain helpers of dustmasker, segmasker and
// convert2blastmask (options from -xmlhelp, argv builder, the masking-information rule), MaskingDraft
// (argv per tool, per-tool values, -out, the shared FASTA input, blastdb, required options, help and
// version, the mask-data marking), the Coordinator's marking of the output, the makeblastdb hand-over,
// and the FakeEngine's rows for the three tools.
import { gzipSync } from 'node:zlib';
import { describe, expect, it, vi } from 'vitest';
import { sha256Hex } from 'losat';
import { parseToolDescription } from 'losat/domain';
import type { ToolEngine } from 'losat/ports';
import { Coordinator, isToolRun, type EnqueueAllResult, type ToolRequest } from '../../src/application/coordinator';
import { MaskDataHandover, type MaskDataEntry } from '../../src/application/mask-data-handover';
import { MaskingDraft } from '../../src/application/masking-draft';
import { SessionFiles } from '../../src/application/session-files';
import { ToolDraft } from '../../src/application/tool-draft';
import {
  buildMaskingArgv,
  changedMaskingOptionsCount,
  isMaskDataFile,
  maskDataOnStdout,
  writesMaskData,
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

  it('finds masking information by one rule: the effective -outfmt of a masking run is a maskinfo format', () => {
    expect(writesMaskData('dustmasker', 'mask', ['dustmasker', '-in', 'a.fa'], DUST)).toBe(false);
    expect(writesMaskData('dustmasker', 'mask', ['dustmasker', '-in', 'a.fa', '-outfmt', 'maskinfo_asn1_bin'], DUST)).toBe(true);
    expect(writesMaskData('segmasker', 'mask', ['segmasker', '-outfmt=maskinfo_xml'], SEG)).toBe(true);
    expect(writesMaskData('segmasker', 'mask', ['segmasker', '-outfmt', 'seqloc_xml'], SEG)).toBe(false);
    // convert2blastmask's default -outfmt is maskinfo_asn1_text.
    expect(writesMaskData('convert2blastmask', 'mask', ['convert2blastmask', '-in', 'a.fa'], C2B)).toBe(true);
    expect(writesMaskData('convert2blastmask', 'mask', ['convert2blastmask', '-outfmt', 'interval'], C2B)).toBe(false);
    expect(writesMaskData('convert2blastmask', 'mask', ['convert2blastmask'], undefined)).toBe(false);
    // windowmasker's Mask stage writes it too; its counts stages and blast_formatter never do.
    expect(writesMaskData('windowmasker', 'mask', ['windowmasker', '-ustat', 'c', '-outfmt', 'maskinfo_asn1_text'], WM)).toBe(true);
    expect(writesMaskData('windowmasker', 'mask', ['windowmasker', '-ustat', 'c'], WM)).toBe(false);
    expect(writesMaskData('windowmasker', 'mk_counts', ['windowmasker', '-mk_counts', '-outfmt', 'maskinfo_xml'], WM)).toBe(false);
    expect(writesMaskData('blast_formatter', 'format', ['blast_formatter', '-outfmt', 'maskinfo_xml'], undefined)).toBe(false);

    expect(maskDataOnStdout(['dustmasker', '-in', 'a.fa'])).toBe(true);
    expect(maskDataOnStdout(['dustmasker', '-out', '-'])).toBe(true);
    expect(maskDataOnStdout(['dustmasker', '-out', 'm.asnb'])).toBe(false);
    expect(isMaskDataFile(['dustmasker', '-out', './m.asnb'], 'm.asnb')).toBe(true);
    expect(isMaskDataFile(['dustmasker', '-out', 'm.asnb', '-logfile', 'run.log'], 'run.log')).toBe(false);
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
    expect(combined.writesMaskData).toBeUndefined();

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
    // Its default -outfmt (maskinfo_asn1_text) is masking information.
    expect(queued[0]!.writesMaskData).toBe(true);
  });

  it('marks the runs that write masking information from the argv, else from the -xmlhelp default', async () => {
    const { draft, queued } = setup({ validate: 'ok' });
    draft.ensureLoaded();
    await described(draft, 'dustmasker');
    draft.setFastaPaste('>q\nACGT\n');
    draft.setField('-outfmt', 'maskinfo_asn1_bin');
    draft.setOutName('dust.asnb');
    expect((await draft.submit()).ok).toBe(true);
    expect(queued[0]!.argv).toEqual(['dustmasker', '-in', 'input.fa', '-out', 'dust.asnb', '-outfmt', 'maskinfo_asn1_bin']);
    expect(queued[0]!.writesMaskData).toBe(true);
    draft.setField('-outfmt', 'interval');
    expect((await draft.submit()).ok).toBe(true);
    expect(queued[1]!.writesMaskData).toBeUndefined();
    draft.setExtra('-outfmt maskinfo_xml');
    expect((await draft.submit()).ok).toBe(true);
    expect(queued[2]!.writesMaskData).toBe(true);
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
    expect(queued[0]!.writesMaskData).toBeUndefined();
    draft.setTool('convert2blastmask');
    draft.setExtra('-help');
    expect((await draft.submit()).ok).toBe(true);
    expect(queued[1]).toMatchObject({ tool: 'convert2blastmask', argv: ['convert2blastmask', '-help'], prints: 'help' });
    expect(queued[1]!.writesMaskData).toBeUndefined();
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

describe('Coordinator marks the output that holds masking information', () => {
  async function runMasking(argv: readonly string[], writesMaskData: boolean) {
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
        ...(writesMaskData ? { writesMaskData: true as const } : {}),
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

  it('marks standard output without -out, the -out file with it, and nothing for a run that writes none', async () => {
    const stdout = await runMasking(['convert2blastmask', '-in', 'input.fa', '-masking_algorithm', 'dust', '-masking_options', 'o'], true);
    expect(stdout.run.snapshot.writesMaskData).toBe(true);
    expect(stdout.run.snapshot.stage).toBe('mask');
    expect(stdout.run.record.stdout?.isMaskData).toBe(true);
    expect(stdout.run.record.stderr?.isMaskData).toBeUndefined();
    const preview = await stdout.coordinator.readToolOutputPreview(stdout.run.snapshot.runId, 'stdout', 1024);
    expect(new TextDecoder().decode(preview.bytes)).toContain(FAKE_MARKER);

    const file = await runMasking(['dustmasker', '-in', 'input.fa', '-out', 'm.asnt', '-outfmt', 'maskinfo_asn1_text', '-logfile', 'run.log'], true);
    expect(file.run.record.stdout?.isMaskData).toBeUndefined();
    expect(file.run.record.files?.map((entry) => [entry.name, entry.isMaskData === true])).toEqual([
      ['m.asnt', true],
      ['run.log', false],
    ]);

    const none = await runMasking(['dustmasker', '-in', 'input.fa'], false);
    expect(none.run.snapshot.writesMaskData).toBeUndefined();
    expect(none.run.record.stdout?.isMaskData).toBeUndefined();
  });
});

describe('MaskDataHandover', () => {
  const entry: MaskDataEntry = { runId: 'r1', runNumber: 3, key: 'stdout', name: 'stdout', size: 42 };

  it('is not available without a makeblastdb form, and then hands nothing over', () => {
    const handover = new MaskDataHandover();
    expect(handover.available()).toBe(false);
    expect(() => handover.use(entry)).not.toThrow();
  });

  it('is available with a target, and hands the output over to it', () => {
    const useMaskData = vi.fn();
    const handover = new MaskDataHandover({ useMaskData });
    expect(handover.available()).toBe(true);
    handover.use(entry);
    expect(useMaskData).toHaveBeenCalledWith(entry);
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
    await expect(fake.describeTool('makeblastdb')).rejects.toThrow('unknown tool');
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
