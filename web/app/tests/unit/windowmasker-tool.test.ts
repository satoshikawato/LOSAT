// Unit tests for the WindowMasker tool tab: argv builder, stage option filtering from -xmlhelp,
// blastdb name detection, UTF-8 preview decoding, ToolDraft, SearchDraft.useCounts, shared
// Coordinator queue (search -> tool -> search, cancel, non-zero exit, Data worker end), and
// DataService tool output storage / OPFS-backed session inflation.
import { gzipSync } from 'node:zlib';
import { describe, expect, it } from 'vitest';
import { sha256Hex } from 'losat';
import { parseToolDescription } from 'losat/domain';
import {
  RunCancelledError,
  type EngineGateway,
  type EngineRunRequest,
  type ToolEngine,
  type ToolRunOptions,
  type ToolRunResult,
} from 'losat/ports';
import { Coordinator, isToolRun, type SearchRequest, type ToolRequest } from '../../src/application/coordinator';
import { CountsLoader } from '../../src/application/counts';
import { SearchDraft } from '../../src/application/draft';
import { ResultsBrowser } from '../../src/application/results';
import { SessionFiles } from '../../src/application/session-files';
import { countsLockName, type CountsMadeMeta } from '../../src/domain/counts';
import { webCountsLocks } from '../../src/infra/browser/counts-locks';
import { OpfsCountsStore } from '../../src/infra/data/counts-store';
import { FakeCountsDownloader } from '../../src/infra/fake/fake-counts';
import type { CountsLocks } from '../../src/ports/counts';
import { MemoryOpfs } from './support/memory-opfs';
import { ToolDraft } from '../../src/application/tool-draft';
import {
  buildWindowMaskerArgv,
  changedToolOptionsCount,
  isCountsOutputFile,
  refusedToolFields,
  toolCommandWords,
  toolFileOptionLabel,
  toolValidationSummary,
} from '../../src/domain/tool-argv';
import { primaryBlastDbName } from '../../src/domain/blast-db';
import {
  decodeUtf8Preview,
  stageOptionParameters,
} from '../../src/domain/tools';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import WINDOWMASKER_XMLHELP_DOC from '../../src/infra/fake/windowmasker-xmlhelp.json';
import type { Downloader } from '../../src/ports/download';
import { NodeInflater } from './support/node-inflater';

const WM_DESC = parseToolDescription(WINDOWMASKER_XMLHELP_DOC.xml);
const NCBI_RULE = '='.repeat(72);
const sleep = (ms: number) => new Promise<void>((resolve) => setTimeout(resolve, ms));

function stubEngine(base: FakeEngine, overrides: Partial<EngineGateway & ToolEngine>): EngineGateway & ToolEngine {
  return {
    describe: (p) => base.describe(p),
    describeTool: (t) => base.describeTool(t),
    start: () => base.start(),
    validate: (a) => base.validate(a),
    validateTool: async () => ({ ok: true }),
    cancel: (id) => base.cancel(id),
    run: (req, out, onPhase) => base.run(req, out, onPhase),
    runTool: (argv, opts) => base.runTool(argv, opts),
    ...overrides,
  };
}

function sampleMkCountsRequest(): ToolRequest {
  return {
    tool: 'windowmasker',
    stage: 'mk_counts',
    argv: ['windowmasker', '-mk_counts', '-in', 'a.fa', '-out', 'wmasker.ascii'],
    directory: [{ name: 'a.fa', file: new Blob(['>a\nACGT\n']), origin: { kind: 'picked' } }],
  };
}

async function nextGate(gates: Array<() => void>): Promise<void> {
  while (gates.length === 0) await sleep(5);
  gates.shift()!();
}

function makeData(store = new MemoryBlockStore()) {
  let tok = 0;
  const data = new DataService({
    store,
    scanner: new FakeScanner(),
    inflater: new NodeInflater(),
    checker: new FakeInputChecker(),
    digest: sha256Hex,
    newToken: () => `t-${++tok}`,
    cleanup: Promise.resolve({ state: 'done', removedSessions: 0 }),
  });
  return { data, store };
}

function makeCoordinator(
  engine: EngineGateway & Partial<ToolEngine> = new FakeEngine(),
  dataWorkerEnded?: Promise<Error>,
  countsLocks?: CountsLocks,
) {
  const { data, store } = makeData();
  const saved: Array<{ readonly fileName: string; readonly bytes: Uint8Array }> = [];
  const downloader: Downloader = {
    save: (fileName, bytes) => saved.push({ fileName, bytes }),
    saveBlob: (fileName, blob) => {
      void blob.arrayBuffer().then((buf) => saved.push({ fileName, bytes: new Uint8Array(buf) }));
    },
  };
  let id = 0;
  const coordinator = new Coordinator({
    engine,
    data,
    downloader,
    now: () => 10_000,
    newRunId: () => `run-${++id}`,
    digest: sha256Hex,
    ...(dataWorkerEnded === undefined ? {} : { dataWorkerEnded }),
    ...(countsLocks === undefined ? {} : { countsLocks }),
  });
  return { coordinator, data, store, saved };
}

async function waitTerminal(coordinator: Coordinator, runId: string) {
  for (let i = 0; i < 400; i++) {
    const run = coordinator.state.get().runs.find((entry) => entry.snapshot.runId === runId);
    if (run !== undefined && ['completed', 'cancelled', 'failed'].includes(run.status)) return run;
    await new Promise((resolve) => setTimeout(resolve, 5));
  }
  throw new Error(`run ${runId} did not finish`);
}

describe('windowmasker domain helpers', () => {
  it('filters options per stage from -xmlhelp exclusions in usage order', () => {
    const flagsOf = (stage: 'mk_counts' | 'mask' | 'convert') =>
      stageOptionParameters(WM_DESC, stage).map((param) => param.flag);

    const mk = flagsOf('mk_counts');
    expect(mk).toContain('-checkdup');
    expect(mk).toContain('-fa_list');
    expect(mk).toContain('-unit');
    expect(mk).toContain('-smem');
    expect(mk).toContain('-ids');
    expect(mk).not.toContain('-dust');
    expect(mk).not.toContain('-window');
    expect(mk).not.toContain('-outfmt');
    expect(mk).not.toContain('-sformat');

    const mask = flagsOf('mask');
    expect(mask).toContain('-window');
    expect(mask).toContain('-dust');
    expect(mask).toContain('-dust_level');
    expect(mask).toContain('-t_thres');
    expect(mask).toContain('-exclude_ids');
    expect(mask).not.toContain('-checkdup');
    expect(mask).not.toContain('-unit');
    expect(mask).not.toContain('-sformat');

    const conv = flagsOf('convert');
    expect(conv).toEqual(['-meta', '-t_high_pct', '-t_thres_pct', '-t_extend_pct', '-t_low_pct', '-smem']);
  });

  it('builds argv omitting defaults and respecting -out and -sformat', () => {
    const mkDefault = buildWindowMaskerArgv({
      stage: 'mk_counts',
      inputName: 'genome.fa',
      values: { '-sformat': 'ascii', '-checkdup': 'false' },
      description: WM_DESC,
    });
    expect(mkDefault).toEqual(['windowmasker', '-mk_counts', '-in', 'genome.fa', '-out', 'wmasker.ascii']);
    expect(toolCommandWords(mkDefault)).toEqual(['losat', ...mkDefault]);

    const mkStdout = buildWindowMaskerArgv({
      stage: 'mk_counts',
      inputName: 'genome.fa',
      outName: '-',
      values: { '-sformat': 'obinary', '-unit': '9' },
      description: WM_DESC,
    });
    expect(mkStdout).toEqual(['windowmasker', '-mk_counts', '-in', 'genome.fa', '-unit', '9', '-sformat', 'obinary']);

    const maskArgv = buildWindowMaskerArgv({
      stage: 'mask',
      inputName: 'seq.fa',
      ustatName: 'wmasker.obinary',
      outName: '',
      values: { '-outfmt': 'fasta', '-dust': 'true' },
      description: WM_DESC,
    });
    expect(maskArgv).toEqual([
      'windowmasker',
      '-ustat',
      'wmasker.obinary',
      '-in',
      'seq.fa',
      '-outfmt',
      'fasta',
      '-dust',
      'true',
    ]);

    const convArgv = buildWindowMaskerArgv({
      stage: 'convert',
      inputName: 'u4.ascii',
      outName: 'u4.obinary',
      values: { '-sformat': 'obinary', '-meta': 'demo' },
      description: WM_DESC,
    });
    expect(convArgv).toEqual(['windowmasker', '-convert', '-in', 'u4.ascii', '-out', 'u4.obinary', '-meta', 'demo', '-sformat', 'obinary']);

    expect(changedToolOptionsCount('mk_counts', { '-unit': '9', '-sformat': 'obinary' }, '', WM_DESC)).toBe(1);
    expect(changedToolOptionsCount('mk_counts', { '-unit': '9' }, '-logfile run.log', WM_DESC)).toBe(2);
  });

  it('labels output files by option and recognizes counts outputs', () => {
    const argv = ['windowmasker', '-mk_counts', '-in', 'a.fa', '-out', 'x.txt', '-logfile', './x.txt'];
    expect(toolFileOptionLabel(argv, 'x.txt')).toBe('-out, -logfile');
    expect(isCountsOutputFile('mk_counts', argv, 'x.txt')).toBe(true);
    expect(isCountsOutputFile('mask', ['windowmasker', '-ustat', 'c', '-out', 'masked.fa'], 'masked.fa')).toBe(false);
    expect(isCountsOutputFile('convert', ['windowmasker', '-convert', '-in', 'a', '-out', 'b'], 'b')).toBe(true);
    expect(isCountsOutputFile('mk_counts', ['windowmasker', '-mk_counts', '-out', 'a', '-logfile', 'log'], 'log')).toBe(false);
  });

  it('detects BLAST database names preferring alias files over volume files', () => {
    expect(primaryBlastDbName(['db/wm.nin', 'db/wm.nsq', 'db/wm.nhr'])).toBe('db/wm');
    expect(primaryBlastDbName(['db/vol.00.nin', 'db/vol.01.nin', 'db/all.nal'])).toBe('db/all');
    expect(primaryBlastDbName(['plain.fa', 'ids.txt'])).toBeUndefined();
  });

  it('decodes UTF-8 preview slices safely across multi-byte boundaries and rejects binary', () => {
    const utf8 = new TextEncoder().encode('hello日本語world');
    // Cut inside the 3-byte UTF-8 sequence for '日' (starts at byte 5).
    const cut = utf8.subarray(0, 7);
    expect(decodeUtf8Preview(cut, true)).toBe('hello');
    expect(decodeUtf8Preview(cut, false)).toBeUndefined();
    expect(decodeUtf8Preview(new Uint8Array([0x41, 0x00, 0x42]), false)).toBeUndefined();
    expect(decodeUtf8Preview(utf8, false)).toBe('hello日本語world');
  });

  it('extracts concise CLI error lines and maps refused fields from ncbi_app::args stderr', () => {
    const unitErr = [
      'USAGE',
      '  windowmasker [-h] [-unit unit_length] [-dust use_dust]',
      NCBI_RULE,
      '',
      "Error: Argument \"unit\". Illegal value, expected 1..16:  `0'",
      "Error:  Argument \"unit\". Illegal value, expected 1..16:  `0'",
    ].join('\n');
    expect(toolValidationSummary(unitErr)).toBe('Error: Argument "unit". Illegal value, expected 1..16:  `0\'');
    expect(refusedToolFields(unitErr, '')).toEqual({ extra: false, flags: new Set(['-unit']) });

    const dustErr = [
      'USAGE',
      '  windowmasker [-mk_counts] [-dust use_dust]',
      NCBI_RULE,
      '',
      "Error: Argument \"mk_counts\". Incompatible with argument:  `dust'",
    ].join('\n');
    expect(refusedToolFields(dustErr, '-dust T')).toEqual({ extra: true, flags: new Set(['-mk_counts']) });

    const unknownErr = '========================================================================\n\nError: Unknown argument: "nosuchopt"';
    expect(refusedToolFields(unknownErr, '-nosuchopt')).toEqual({ extra: true, flags: new Set() });

    const dupErr = '========================================================================\n\nError: Argument with this name is defined already: window';
    expect(refusedToolFields(dupErr, '-window 20')).toEqual({ extra: true, flags: new Set() });
  });
});

describe('ToolDraft and SearchDraft.useCounts', () => {
  function setupDraft() {
    const { coordinator, data } = makeCoordinator(new FakeEngine());
    const submitted: ToolRequest[][] = [];
    const draft = new ToolDraft({
      engine: new FakeEngine(),
      data,
      enqueueToolAll: async (requests) => {
        submitted.push([...requests]);
        return coordinator.enqueueToolAll(requests);
      },
      debounceMs: 0,
    });
    draft.ensureLoaded();
    return { draft, coordinator, data, submitted };
  }

  it('updates default -out when -sformat changes until edited, and combines FASTA files with newline check', async () => {
    const { draft, submitted } = setupDraft();
    await new Promise((resolve) => setTimeout(resolve, 20));
    expect(draft.state.get().outNames.mk_counts).toBe('wmasker.ascii');

    draft.setField('-sformat', 'obinary');
    expect(draft.state.get().outNames.mk_counts).toBe('wmasker.obinary');

    draft.setOutName('custom.obinary');
    draft.setField('-sformat', 'oascii');
    expect(draft.state.get().outNames.mk_counts).toBe('custom.obinary');

    draft.addFastaFiles([
      new File(['>a\nACGT'], 'a.fa'),
      new File([gzipSync('>b\nTGCA\n')], 'b.fa.gz'),
    ]);
    await new Promise((resolve) => setTimeout(resolve, 30));
    expect(draft.state.get().fastaFiles.map((f) => [f.name, f.inflatedFrom, f.status])).toEqual([
      ['a.fa', undefined, 'ready'],
      ['b.fa', 'b.fa.gz', 'ready'],
    ]);

    const res = await draft.submit();
    expect(res.ok).toBe(true);
    expect(submitted).toHaveLength(1);
    const req = submitted[0]![0]!;
    expect(req.argv).toContain('combined_input.fa');
    const combined = req.directory.find((d) => d.name === 'combined_input.fa');
    expect(combined?.origin).toEqual({ kind: 'combined', parts: ['a.fa', 'b.fa'] });
    expect(await combined!.file.text()).toBe('>a\nACGT\n>b\nTGCA\n');
  });

  it('refuses duplicate filenames and detects blastdb names from opened folder files', async () => {
    const { draft } = setupDraft();
    draft.setField('-infmt', 'blastdb');
    draft.addRunFiles(
      [new File([new Uint8Array([1])], 'wm.nin'), new File([new Uint8Array([2])], 'wm.nsq')],
      ['db/wm.nin', 'db/wm.nsq'],
    );
    expect(draft.state.get().blastDbIn).toBe('db/wm');
    expect(draft.state.get().runFilesError).toBeUndefined();

    draft.addRunFiles([new File([new Uint8Array([3])], 'wm.nin')], ['db/wm.nin']);
    expect(draft.state.get().runFilesError).toContain('db/wm.nin');
    expect(draft.state.get().runFiles).toHaveLength(2);
  });

  it('SearchDraft.useCounts switches to BLASTN with -window_masker_db file attached', () => {
    const { data } = makeData();
    const draft = new SearchDraft({
      engine: new FakeEngine(),
      data,
      enqueueAll: async () => ({ ok: true, runIds: ['r1'] }),
      fetcher: { fetch: async () => ({ ok: false, reason: 'not-found' }) },
      now: () => 0,
      debounceMs: 0,
    });
    draft.setProgram('blastp');
    const counts = new File(['counts'], 'wmasker.ascii');
    draft.useCounts(counts);
    expect(draft.state.get().program).toBe('blastn');
    expect(draft.state.get().values.blastn['-window_masker_db']).toBe('wmasker.ascii');
    expect(draft.state.get().files.blastn['-window_masker_db']).toBe(counts);
  });

  it('SearchDraft.useTaxid switches to BLASTN with -window_masker_taxid set', () => {
    const { data } = makeData();
    const draft = new SearchDraft({
      engine: new FakeEngine(),
      data,
      enqueueAll: async () => ({ ok: true, runIds: ['r1'] }),
      fetcher: { fetch: async () => ({ ok: false, reason: 'not-found' }) },
      now: () => 0,
      debounceMs: 0,
    });
    draft.setProgram('blastp');
    draft.useTaxid(562);
    expect(draft.state.get().program).toBe('blastn');
    expect(draft.state.get().values.blastn['-window_masker_taxid']).toBe('562');
    expect(draft.state.get().values.blastn['-window_masker_db']).toBe('');
  });

  it('resolves saved taxonomy-ID counts for -ustat and -in under windowmasker/<taxid>/wmasker.obinary', async () => {
    const { coordinator, data } = makeCoordinator(new FakeEngine());
    const savedCounts = new File([new Uint8Array([1, 2, 3, 4])], 'wmasker.obinary');
    const submitted: ToolRequest[][] = [];
    const draft = new ToolDraft({
      engine: new FakeEngine(),
      data,
      enqueueToolAll: async (requests) => {
        submitted.push([...requests]);
        return coordinator.enqueueToolAll(requests);
      },
      readSavedCounts: async (taxid) => (taxid === 562 ? savedCounts : undefined),
      debounceMs: 0,
    });
    draft.ensureLoaded();
    await sleep(20);

    // Mask stage with saved counts for 562.
    draft.setStage('mask');
    draft.setUstatSaved(562, savedCounts.size);
    draft.setFastaPaste('>q\nACGTACGT\n');
    expect(draft.previewArgv()).toContain('windowmasker/562/wmasker.obinary');
    expect((await draft.submit()).ok).toBe(true);
    const maskReq = submitted[0]![0]!;
    expect(maskReq.argv).toEqual(['windowmasker', '-ustat', 'windowmasker/562/wmasker.obinary', '-in', 'input.fa']);
    expect(maskReq.directory.find((d) => d.name === 'windowmasker/562/wmasker.obinary')).toEqual({
      name: 'windowmasker/562/wmasker.obinary',
      file: savedCounts,
      origin: { kind: 'saved', taxid: 562 },
    });

    // Convert stage with saved counts for 562.
    draft.setStage('convert');
    draft.setConvertSaved(562, savedCounts.size);
    expect((await draft.submit()).ok).toBe(true);
    const convReq = submitted[1]![0]!;
    expect(convReq.argv).toEqual(['windowmasker', '-convert', '-in', 'windowmasker/562/wmasker.obinary', '-out', 'wmasker.ascii']);
    expect(convReq.directory).toEqual([
      { name: 'windowmasker/562/wmasker.obinary', file: savedCounts, origin: { kind: 'saved', taxid: 562 } },
    ]);

    // Missing saved counts reports a clear error.
    draft.setConvertSaved(9999, 10);
    const missing = await draft.submit();
    expect(missing.ok).toBe(false);
    expect(missing.ok ? '' : missing.message).toContain('9999');
  });
});

describe('Coordinator single queue and DataService tool storage', () => {
  it('runs search -> tool -> search in order with shared numbering and handles queued cancel', async () => {
    const order: string[] = [];
    const base = new FakeEngine();
    const gates: Array<() => void> = [];
    const waitGate = () => new Promise<void>((resolve) => gates.push(resolve));
    const engine = stubEngine(base, {
      validateTool: (a) => base.validateTool(a),
      async run(req: EngineRunRequest, out, onPhase) {
        order.push(`search:${req.runId}:start`);
        await waitGate();
        const info = await base.run(req, out, onPhase);
        order.push(`search:${req.runId}:end`);
        return info;
      },
      async runTool(argv: readonly string[], opts?: ToolRunOptions): Promise<ToolRunResult> {
        order.push(`tool:${argv[1]}:start`);
        await waitGate();
        const res = await base.runTool(argv, opts);
        order.push(`tool:${argv[1]}:end`);
        return res;
      },
    });
    const { coordinator } = makeCoordinator(engine);
    const searchReq: SearchRequest = {
      program: 'blastn',
      query: { text: '>q\nACGT\n' },
      subject: { text: '>s\nACGT\n' },
      parameters: [],
      requestedThreads: 1,
    };
    const toolReq = sampleMkCountsRequest();

    const r1 = (await coordinator.enqueue(searchReq)).runId!;
    const r2 = (await coordinator.enqueueTool(toolReq)).runId!;
    const r3 = (await coordinator.enqueueTool(toolReq)).runId!;
    const r4 = (await coordinator.enqueue(searchReq)).runId!;
    expect(coordinator.state.get().runs.map((r) => r.snapshot.number)).toEqual([1, 2, 3, 4]);

    // Cancel queued tool run r3 before it starts.
    coordinator.cancel(r3);
    expect(coordinator.state.get().runs[2]!.status).toBe('cancelled');

    // Release r1 -> r2 -> r4 sequentially.
    await nextGate(gates);
    await waitTerminal(coordinator, r1);
    await nextGate(gates);
    const run2 = await waitTerminal(coordinator, r2);
    expect(isToolRun(run2) && run2.status === 'completed' && run2.record.exitStatus === 0).toBe(true);
    expect(coordinator.sessionMadeCounts()).toHaveLength(1);

    await nextGate(gates);
    await waitTerminal(coordinator, r4);

    expect(order).toEqual([
      `search:${r1}:start`,
      `search:${r1}:end`,
      'tool:-mk_counts:start',
      'tool:-mk_counts:end',
      `search:${r4}:start`,
      `search:${r4}:end`,
    ]);
  });

  it('commits outputs and stderr of a non-zero exit status tool run and opens them in ResultsBrowser', async () => {
    const base = new FakeEngine();
    const engine = stubEngine(base, {
      cancel: () => undefined,
      async runTool(_argv, opts) {
        const enc = new TextEncoder();
        await opts?.sinks?.onStderr?.(enc.encode('ERROR: invalid counts header\n'));
        const sink = await opts?.sinks?.onFile?.('partial.out');
        await sink?.write(enc.encode('partial'), 0);
        await sink?.close();
        return { status: 3, stdout: new Blob(), stderr: new Blob(), files: [] };
      },
    });
    const { coordinator, data } = makeCoordinator(engine);
    const results = new ResultsBrowser({
      data,
      describe: (p) => engine.describe(p),
      describeTool: (t) => engine.describeTool(t),
      runs: coordinator.state,
      verification: { ncbi: '2.17.0', sources: [], programs: {} },
    });
    const failedQueued = await coordinator.enqueueTool({
      tool: 'windowmasker',
      stage: 'mask',
      argv: ['windowmasker', '-ustat', 'bad.counts', '-in', 'a.fa', '-out', 'partial.out'],
      directory: [],
    });
    const runId = failedQueued.runId!;
    const finished = await waitTerminal(coordinator, runId);
    expect(finished.status).toBe('failed');
    expect(isToolRun(finished) && finished.record.exitStatus).toBe(3);

    await results.open(runId);
    expect(results.state.get().phase).toBe('ready');
    expect(results.state.get().toolRun?.record.exitStatus).toBe(3);
    const stderrPreview = await coordinator.readToolOutputPreview(runId, 'stderr', 1024);
    expect(new TextDecoder().decode(stderrPreview.bytes)).toContain('invalid counts header');
    const filePreview = await coordinator.readToolOutputPreview(runId, 'file:0', 1024);
    expect(new TextDecoder().decode(filePreview.bytes)).toBe('partial');
  });

  it('cancels a running tool run and handles Data worker end cleanly', async () => {
    const base = new FakeEngine({ phaseDelayMs: 500 });
    const { coordinator } = makeCoordinator(base);
    const { runId } = await coordinator.enqueueTool(sampleMkCountsRequest());
    await sleep(20);
    coordinator.cancel(runId!);
    const cancelled = await waitTerminal(coordinator, runId!);
    expect(cancelled.status).toBe('cancelled');

    let endWorker!: (err: Error) => void;
    const ended = new Promise<Error>((resolve) => (endWorker = resolve));
    const hangingEngine = stubEngine(base, {
      cancel: () => undefined,
      runTool: (_argv, opts) =>
        new Promise((_resolve, reject) => {
          opts?.signal?.addEventListener('abort', () => reject(new RunCancelledError('tool')));
        }),
    });
    const second = makeCoordinator(hangingEngine, ended);
    const { runId: r2 } = await second.coordinator.enqueueTool({
      tool: 'windowmasker',
      stage: 'mk_counts',
      argv: ['windowmasker', '-mk_counts'],
      directory: [],
    });
    await sleep(20);
    endWorker(new Error('Data worker stopped'));
    const failed = await waitTerminal(second.coordinator, r2!);
    expect(failed.status).toBe('failed');
    expect(failed.record.error).toContain('Data worker stopped');
  });

  it('queues a tool run immediately behind an active run even when validateTool would wait behind the worker', async () => {
    const base = new FakeEngine();
    let releaseFirst!: () => void;
    const gate = new Promise<void>((resolve) => (releaseFirst = resolve));
    let validateCalledWhileActive = false;
    let activeRun = false;
    const blockingEngine = stubEngine(base, {
      validateTool: async (argv) => {
        if (activeRun) {
          validateCalledWhileActive = true;
          await gate;
        }
        return base.validateTool(argv);
      },
      cancel: () => undefined,
      async runTool(argv, opts) {
        if (!activeRun && argv.includes('slow.fa')) {
          activeRun = true;
          await gate;
          activeRun = false;
        }
        return base.runTool(argv, opts);
      },
    });
    const { coordinator } = makeCoordinator(blockingEngine);
    const first = await coordinator.enqueueTool({
      tool: 'windowmasker',
      stage: 'mk_counts',
      argv: ['windowmasker', '-mk_counts', '-in', 'slow.fa', '-out', 'a.ascii'],
      directory: [{ name: 'slow.fa', file: new Blob(['>s\nACGT\n']), origin: { kind: 'picked' } }],
    });
    expect(first.ok).toBe(true);
    await new Promise((r) => setTimeout(r, 20));
    expect(coordinator.isBusy()).toBe(true);

    // Enqueueing a second tool run must NOT block waiting on validateTool behind the active run.
    const second = await coordinator.enqueueTool({
      tool: 'windowmasker',
      stage: 'mk_counts',
      argv: ['windowmasker', '-mk_counts', '-in', 'fast.fa', '-out', 'b.ascii'],
      directory: [{ name: 'fast.fa', file: new Blob(['>f\nACGT\n']), origin: { kind: 'picked' } }],
    });
    expect(second.ok).toBe(true);
    expect(validateCalledWhileActive).toBe(false);
    expect(coordinator.state.get().runs.find((r) => r.snapshot.runId === second.runId)?.status).toBe('queued');

    releaseFirst();
    expect((await waitTerminal(coordinator, first.runId!)).status).toBe('completed');
    expect((await waitTerminal(coordinator, second.runId!)).status).toBe('completed');
  });

  it('supports POSIX writeAt/truncate on tool output files and OPFS-backed session inflation', async () => {
    const { data, store } = makeData();
    const port = await data.openToolRun('tool-posix');
    const enc = new TextEncoder();
    const call = (msg: object) =>
      new Promise<{ ok: boolean; summary?: { files: Array<{ key: string; size: number }> } }>((resolve) => {
        port.onmessage = (event) => resolve(event.data);
        port.postMessage(msg);
      });
    expect((await call({ id: 1, op: 'open', key: 'file:0', name: 'both.txt' })).ok).toBe(true);
    expect((await call({ id: 2, op: 'write', key: 'file:0', offset: 0, bytes: enc.encode('log-line\n') })).ok).toBe(true);
    expect((await call({ id: 3, op: 'truncate', key: 'file:0', length: 0 })).ok).toBe(true);
    expect((await call({ id: 4, op: 'write', key: 'file:0', offset: 0, bytes: enc.encode('COUNTS') })).ok).toBe(true);
    expect((await call({ id: 5, op: 'write', key: 'file:0', offset: 8, bytes: enc.encode('END') })).ok).toBe(true);
    expect((await call({ id: 6, op: 'close', key: 'file:0' })).ok).toBe(true);
    const committed = await call({ id: 7, op: 'commit' });
    expect(committed.ok).toBe(true);
    expect(committed.summary?.files[0]?.size).toBe(11);
    const bytes = await data.readToolOutputRange('tool-posix', 'file:0', 0, 11);
    expect([...bytes]).toEqual([...enc.encode('COUNTS'), 0, 0, ...enc.encode('END')]);
    const file = await data.getToolOutputFile('tool-posix', 'file:0', 'both.txt');
    expect(file.name).toBe('both.txt');
    expect(file.size).toBe(11);
    await data.deleteRun('tool-posix');
    expect(store.usage()).toBe(0);

    const ref = await data.inflateToSessionFile(new File([gzipSync('>x\nACGT\n')], 'x.fa.gz'));
    expect(ref.name).toBe('x.fa');
    expect(ref.inflatedFrom).toBe('x.fa.gz');
    expect(ref.token).toBeDefined();
    expect(await ref.file.text()).toBe('>x\nACGT\n');
    expect(store.usage()).toBe(8);
    await data.dropSessionFile(ref.token!);
    expect(store.usage()).toBe(0);
  });

  it('computes Supported / Not yet supported / Development build badge on tool runs in ResultsBrowser', async () => {
    const engine = new FakeEngine();
    const { data } = makeData();
    let id = 0;
    const makeCoord = (path: 'serial' | 'fake') =>
      new Coordinator({
        engine,
        data,
        downloader: { save: () => undefined },
        now: () => 10_000,
        newRunId: () => `r-${++id}`,
        digest: sha256Hex,
        toolRuntime: { path, engineBuild: path === 'fake' ? 'fake-engine' : 'serial-build' },
      });
    const table = {
      ncbi: '2.17.0',
      sources: [],
      programs: {},
      tools: { windowmasker: { optionSets: [JSON.stringify([['-mk_counts', '']])] } },
    };

    const serialCoord = makeCoord('serial');
    const results = new ResultsBrowser({
      data,
      describe: (p) => engine.describe(p),
      describeTool: (t) => engine.describeTool(t),
      runs: serialCoord.state,
      verification: table,
    });
    const supportedId = (await serialCoord.enqueueTool(sampleMkCountsRequest())).runId!;
    await waitTerminal(serialCoord, supportedId);
    await results.open(supportedId);
    await sleep(20);
    expect(results.state.get().toolBadge).toMatchObject({ level: 'certified', label: 'Supported' });

    const outsideId = (
      await serialCoord.enqueueTool({
        tool: 'windowmasker',
        stage: 'mk_counts',
        argv: ['windowmasker', '-mk_counts', '-in', 'a.fa', '-out', 'w.ascii', '-unit', '11'],
        directory: [{ name: 'a.fa', file: new Blob(['>a\nACGT\n']), origin: { kind: 'picked' } }],
      })
    ).runId!;
    await waitTerminal(serialCoord, outsideId);
    await results.open(outsideId);
    await sleep(20);
    expect(results.state.get().toolBadge).toMatchObject({ level: 'outside', label: 'Not yet supported' });

    const helpId = (
      await serialCoord.enqueueTool({
        tool: 'windowmasker',
        stage: 'mk_counts',
        argv: ['windowmasker', '-mk_counts', '-h'],
        directory: [],
        prints: 'help',
      })
    ).runId!;
    await waitTerminal(serialCoord, helpId);
    await results.open(helpId);
    await sleep(20);
    expect(results.state.get().toolBadge).toBeUndefined();

    const fakeCoord = makeCoord('fake');
    const fakeResults = new ResultsBrowser({
      data,
      describe: (p) => engine.describe(p),
      describeTool: (t) => engine.describeTool(t),
      runs: fakeCoord.state,
      verification: table,
    });
    const fakeId = (await fakeCoord.enqueueTool(sampleMkCountsRequest())).runId!;
    await waitTerminal(fakeCoord, fakeId);
    await fakeResults.open(fakeId);
    await sleep(20);
    expect(fakeResults.state.get().toolBadge).toMatchObject({ level: 'development', label: 'Development build' });
  });

  it('keeps an inflated .gz session file alive until queued/running tool runs complete or cancel', async () => {
    const { data } = makeData();
    const dropped: string[] = [];
    const origDrop = data.dropSessionFile.bind(data);
    data.dropSessionFile = async (token: string) => {
      dropped.push(token);
      await origDrop(token);
    };

    const base = new FakeEngine();
    const gates: Array<() => void> = [];
    const readContents: string[] = [];
    const engine = stubEngine(base, {
      validateTool: (a) => base.validateTool(a),
      async runTool(argv, opts) {
        await new Promise<void>((resolve) => gates.push(resolve));
        if (opts?.signal?.aborted === true) throw new RunCancelledError('tool');
        for (const entry of opts?.directory ?? []) {
          readContents.push(`${entry.name}:${await entry.file.text()}`);
        }
        return base.runTool(argv, opts);
      },
    });

    const sessionFiles = new SessionFiles((token) => void data.dropSessionFile(token).catch(() => undefined));
    let id = 0;
    const coordinator = new Coordinator({
      engine,
      data,
      downloader: { save: () => undefined },
      now: () => 10_000,
      newRunId: () => `sf-${++id}`,
      digest: sha256Hex,
      sessionFiles,
    });
    const draft = new ToolDraft({
      engine,
      data,
      enqueueToolAll: (reqs) => coordinator.enqueueToolAll(reqs),
      isBusy: () => coordinator.isBusy(),
      sessionFiles,
      debounceMs: 0,
    });
    draft.ensureLoaded();
    await sleep(20);

    // Case 1: completion after removing .gz from form while run 2 is queued behind run 1.
    const blockerId = (await coordinator.enqueueTool(sampleMkCountsRequest())).runId!;
    await sleep(20);

    draft.addFastaFiles([new File([gzipSync('>g1\nACGT\n')], 'g1.fa.gz')]);
    await sleep(30);
    const token1 = draft.state.get().fastaFiles[0]!.sessionToken!;
    const key1 = draft.state.get().fastaFiles[0]!.key;
    const queuedRes = await draft.submit();
    expect(queuedRes.ok).toBe(true);
    const run2Id = queuedRes.runIds![0]!;

    // Remove from form while run 2 is queued -> must NOT drop yet.
    draft.removeFastaFile(key1);
    expect(dropped).toEqual([]);

    // Release run 1, then run 2 reads g1.fa and completes -> drops exactly once at completion.
    await nextGate(gates);
    await waitTerminal(coordinator, blockerId);
    expect(dropped).toEqual([]);
    await nextGate(gates);
    const run2 = await waitTerminal(coordinator, run2Id);
    expect(run2.status).toBe('completed');
    expect(readContents).toContain('g1.fa:>g1\nACGT\n');
    expect(dropped).toEqual([token1]);

    // Case 2: queued cancel after removing from form drops exactly once.
    const blocker2 = (await coordinator.enqueueTool(sampleMkCountsRequest())).runId!;
    await sleep(20);
    draft.addFastaFiles([new File([gzipSync('>g2\nTGCA\n')], 'g2.fa.gz')]);
    await sleep(30);
    const token2 = draft.state.get().fastaFiles[0]!.sessionToken!;
    const key2 = draft.state.get().fastaFiles[0]!.key;
    const q2 = (await draft.submit()).runIds![0]!;
    draft.removeFastaFile(key2);
    expect(dropped).toEqual([token1]);
    coordinator.cancel(q2);
    expect(dropped).toEqual([token1, token2]);
    await nextGate(gates);
    await waitTerminal(coordinator, blocker2);

    // Case 3: running cancel after removing from form drops exactly once when run ends.
    draft.addFastaFiles([new File([gzipSync('>g3\nCCCC\n')], 'g3.fa.gz')]);
    await sleep(30);
    const token3 = draft.state.get().fastaFiles[0]!.sessionToken!;
    const key3 = draft.state.get().fastaFiles[0]!.key;
    const r3 = (await draft.submit()).runIds![0]!;
    await sleep(20);
    draft.removeFastaFile(key3);
    expect(dropped).toEqual([token1, token2]);
    coordinator.cancel(r3);
    await nextGate(gates);
    expect((await waitTerminal(coordinator, r3)).status).toBe('cancelled');
    expect(dropped).toEqual([token1, token2, token3]);
  });

  it('cleans up Data-side staged tool runs when cancelled in preparing or finalizing', async () => {
    const baseStore = new MemoryBlockStore();
    let sealGate: Promise<void> | undefined;
    let releaseSeal: (() => void) | undefined;
    let sealedWaiting = false;
    const store: typeof baseStore = Object.assign(Object.create(Object.getPrototypeOf(baseStore)), baseStore, {
      async create(path: string) {
        const writer = await baseStore.create(path);
        return {
          append: (bytes: Uint8Array) => writer.append(bytes),
          writeAt: (offset: number, bytes: Uint8Array) => writer.writeAt(offset, bytes),
          truncate: (length: number) => writer.truncate(length),
          async seal() {
            if (path.endsWith('/stdout') && sealGate !== undefined) {
              sealedWaiting = true;
              await sealGate;
            }
            return writer.seal();
          },
        };
      },
      read: (path: string, offset: number, length: number) => baseStore.read(path, offset, length),
      readFile: (path: string, name: string, type?: string) => baseStore.readFile(path, name, type),
      removeAll: (prefix: string) => baseStore.removeAll(prefix),
      usage: () => baseStore.usage(),
      get backend() {
        return baseStore.backend;
      },
    });
    const { data } = makeData(store);

    let openGate: Promise<void> | undefined;
    let releaseOpen: (() => void) | undefined;
    const origOpenToolRun = data.openToolRun.bind(data);
    data.openToolRun = async (runId: string) => {
      if (openGate !== undefined) await openGate;
      return origOpenToolRun(runId);
    };

    const engine = new FakeEngine();
    let id = 0;
    const coordinator = new Coordinator({
      engine,
      data,
      downloader: { save: () => undefined },
      now: () => 10_000,
      newRunId: () => `cf-run-${++id}`,
      digest: sha256Hex,
    });

    // Case 1: Cancel during preparing (delayed openToolRun).
    openGate = new Promise<void>((resolve) => {
      releaseOpen = resolve;
    });
    const prepId = (await coordinator.enqueueTool(sampleMkCountsRequest())).runId!;
    await sleep(20);
    expect(coordinator.state.get().runs.find((r) => r.snapshot.runId === prepId)?.status).toBe('preparing');
    coordinator.cancel(prepId);
    releaseOpen!();
    openGate = undefined;
    const prepDone = await waitTerminal(coordinator, prepId);
    expect(prepDone.status).toBe('cancelled');
    await sleep(30);
    expect(baseStore.usage()).toBe(0);
    await expect(data.readToolOutputRange(prepId, 'stdout', 0, 0)).rejects.toThrow();
    // Re-opening the same runId succeeds only when DataService holds 0 entries for it.
    const probe1 = await origOpenToolRun(prepId);
    probe1.close();
    await data.deleteRun(prepId);

    // Case 2: Cancel during finalizing (delayed commit seal).
    const deletedDuringFinalizing: string[] = [];
    const origDeleteRun = data.deleteRun.bind(data);
    data.deleteRun = async (runId: string) => {
      deletedDuringFinalizing.push(runId);
      return origDeleteRun(runId);
    };
    sealGate = new Promise<void>((resolve) => {
      releaseSeal = resolve;
    });
    const finId = (await coordinator.enqueueTool(sampleMkCountsRequest())).runId!;
    for (let i = 0; i < 40 && !sealedWaiting; i++) await sleep(10);
    expect(sealedWaiting).toBe(true);
    expect(coordinator.state.get().runs.find((r) => r.snapshot.runId === finId)?.status).toBe('finalizing');
    coordinator.cancel(finId);
    const finDone = await waitTerminal(coordinator, finId);
    expect(finDone.status).toBe('cancelled');
    // Port stays open while commit finishes sealing so commit settles and its .then(deleteRun) cleans up.
    expect(deletedDuringFinalizing).toEqual([]);
    releaseSeal!();
    sealGate = undefined;
    await sleep(40);
    expect(deletedDuringFinalizing).toEqual([finId]);
    expect(baseStore.usage()).toBe(0);
    await expect(data.readToolOutputRange(finId, 'stdout', 0, 0)).rejects.toThrow();
    const probe2 = await origOpenToolRun(finId);
    probe2.close();
    await origDeleteRun(finId);
  });

  it('re-reads draft state after awaiting validation in ToolDraft.submit()', async () => {
    const { data } = makeData();
    let releaseInflate!: () => void;
    const inflateGate = new Promise<void>((resolve) => {
      releaseInflate = resolve;
    });
    const origInflate = data.inflateToSessionFile.bind(data);
    data.inflateToSessionFile = async (file: File) => {
      await inflateGate;
      return origInflate(file);
    };

    const base = new FakeEngine();
    let holdFirstValidate = false;
    let releaseValidate!: () => void;
    const validateGate = new Promise<void>((resolve) => {
      releaseValidate = resolve;
    });
    let validateWaiting = false;
    const engine = stubEngine(base, {
      async validateTool(argv) {
        if (holdFirstValidate) {
          holdFirstValidate = false;
          validateWaiting = true;
          await validateGate;
        }
        return base.validateTool(argv);
      },
    });

    const enqueued: ToolRequest[] = [];
    const draft = new ToolDraft({
      engine,
      data,
      enqueueToolAll: async (reqs) => {
        enqueued.push(...reqs);
        return { ok: true, runIds: reqs.map((_, i) => `sub-${i + 1}`) };
      },
      debounceMs: 0,
    });
    draft.ensureLoaded();
    await sleep(20);

    // Start inflating a .gz file (still inflating when submit() is clicked) and trigger a held validation.
    draft.addFastaFiles([new File([gzipSync('>seq\nACGTACGT\n')], 'late.fa.gz')]);
    expect(draft.state.get().fastaFiles[0]?.status).toBe('inflating');
    holdFirstValidate = true;
    draft.setField('-unit', '6');
    const submitPromise = draft.submit();
    for (let i = 0; i < 40 && !validateWaiting; i++) await sleep(10);
    expect(validateWaiting).toBe(true);

    // While submit() is awaiting validation, finish the inflate AND change -unit to 7.
    releaseInflate();
    for (let i = 0; i < 40 && draft.state.get().fastaFiles[0]?.status !== 'ready'; i++) await sleep(10);
    expect(draft.state.get().fastaFiles[0]?.status).toBe('ready');
    draft.setField('-unit', '7');

    // Release validation; submit() must build the request from the post-await state (late.fa + -unit 7).
    releaseValidate();
    const res = await submitPromise;
    expect(res.ok).toBe(true);
    expect(enqueued).toHaveLength(1);
    expect(enqueued[0]!.argv).toContain('-in');
    expect(enqueued[0]!.argv).toContain('late.fa');
    expect(enqueued[0]!.argv).toContain('-unit');
    expect(enqueued[0]!.argv).toContain('7');
    expect(enqueued[0]!.argv).not.toContain('6');
  });

  it('holds submit session-file tokens before buildStageRequests awaits so Remove during combine or -ustat lookup does not drop early', async () => {
    // Part A: Remove clicked while combineFastaFiles awaits reading the tail byte of an inflated .gz file.
    const { data } = makeData();
    const dropped: string[] = [];
    const origDrop = data.dropSessionFile.bind(data);
    data.dropSessionFile = async (token: string) => {
      dropped.push(token);
      await origDrop(token);
    };

    let releaseSlice!: () => void;
    let sliceWaiting = false;
    const sliceGate = new Promise<void>((resolve) => {
      releaseSlice = resolve;
    });
    const origInflate = data.inflateToSessionFile.bind(data);
    data.inflateToSessionFile = async (file: File) => {
      const ref = await origInflate(file);
      if (ref.token === undefined) return ref;
      const realFile = ref.file;
      const origSlice = realFile.slice.bind(realFile);
      realFile.slice = ((start?: number, end?: number, contentType?: string) => {
        const sliced = origSlice(start, end, contentType);
        const origArrayBuffer = sliced.arrayBuffer.bind(sliced);
        sliced.arrayBuffer = async () => {
          sliceWaiting = true;
          await sliceGate;
          return origArrayBuffer();
        };
        return sliced;
      }) as typeof realFile.slice;
      return ref;
    };

    const base = new FakeEngine();
    const runGates: Array<() => void> = [];
    const readTexts: string[] = [];
    const engine = stubEngine(base, {
      async runTool(argv, opts) {
        await new Promise<void>((resolve) => runGates.push(resolve));
        for (const entry of opts?.directory ?? []) {
          readTexts.push(`${entry.name}:${await entry.file.text()}`);
        }
        return base.runTool(argv, opts);
      },
    });
    const sessionFiles = new SessionFiles((token) => void data.dropSessionFile(token).catch(() => undefined));
    let id = 0;
    const coordinator = new Coordinator({
      engine,
      data,
      downloader: { save: () => undefined },
      now: () => 10_000,
      newRunId: () => `comb-${++id}`,
      digest: sha256Hex,
      sessionFiles,
    });
    const enqueue = (reqs: readonly ToolRequest[]) => coordinator.enqueueToolAll(reqs);
    const busy = () => coordinator.isBusy();
    const draft = new ToolDraft({ engine, data, enqueueToolAll: enqueue, isBusy: busy, sessionFiles, debounceMs: 0 });
    draft.ensureLoaded();
    await sleep(20);

    draft.addFastaFiles([
      new File([gzipSync('>first\nACGT')], 'first.fa.gz'),
      new File(['>second\nTGCA\n'], 'second.fa'),
    ]);
    for (let i = 0; i < 40 && draft.state.get().fastaFiles.some((f) => f.status !== 'ready'); i++) await sleep(10);
    const firstEntry = draft.state.get().fastaFiles[0]!;
    const firstTok = firstEntry.sessionToken!;
    expect(firstTok).toBeDefined();

    const submitCombine = draft.submit();
    for (let i = 0; i < 40 && !sliceWaiting; i++) await sleep(10);
    expect(sliceWaiting).toBe(true);

    // User clicks Remove on first.fa.gz while combineFastaFiles is awaiting its tail byte.
    draft.removeFastaFile(firstEntry.key);
    expect(dropped).toEqual([]);

    releaseSlice();
    const combineRes = await submitCombine;
    expect(combineRes.ok).toBe(true);
    expect(dropped.length).toBe(0);

    await nextGate(runGates);
    const combineRun = await waitTerminal(coordinator, combineRes.runIds![0]!);
    expect(combineRun.status).toBe('completed');
    expect(readTexts).toContain('combined_input.fa:>first\nACGT\n>second\nTGCA\n');
    expect(dropped).toEqual([firstTok]);

    // Part B: Remove clicked while Mask stage awaits reading the -ustat saved counts file.
    let releaseUstat!: () => void;
    let ustatWaiting = false;
    const ustatGate = new Promise<void>((resolve) => {
      releaseUstat = resolve;
    });
    const maskDraft = new ToolDraft({
      engine,
      data: {
        inflateToSessionFile: (file) => origInflate(file),
        dropSessionFile: (token) => data.dropSessionFile(token),
        getToolOutputFile: (runId, key, name) => data.getToolOutputFile(runId, key, name),
      },
      enqueueToolAll: enqueue,
      isBusy: busy,
      readSavedCounts: async () => {
        ustatWaiting = true;
        await ustatGate;
        return new File(['counts'], 'wmasker.obinary');
      },
      sessionFiles,
      debounceMs: 0,
    });
    maskDraft.ensureLoaded();
    await sleep(20);
    maskDraft.setStage('mask');
    maskDraft.setUstatSaved(562, 6);
    maskDraft.addFastaFiles([new File([gzipSync('>mask_in\nACGT\n')], 'mask.fa.gz')]);
    for (let i = 0; i < 40 && maskDraft.state.get().fastaFiles[0]?.status !== 'ready'; i++) await sleep(10);
    const maskEntry = maskDraft.state.get().fastaFiles[0]!;
    const maskTok = maskEntry.sessionToken!;
    expect(maskTok).toBeDefined();

    const submitMask = maskDraft.submit();
    for (let i = 0; i < 40 && !ustatWaiting; i++) await sleep(10);
    expect(ustatWaiting).toBe(true);

    maskDraft.removeFastaFile(maskEntry.key);
    expect(dropped).toEqual([firstTok]);

    releaseUstat();
    const maskRes = await submitMask;
    expect(maskRes.ok).toBe(true);
    expect(dropped.length).toBe(1);

    await nextGate(runGates);
    const maskRun = await waitTerminal(coordinator, maskRes.runIds![0]!);
    expect(maskRun.status).toBe('completed');
    expect(readTexts).toContain('mask.fa:>mask_in\nACGT\n');
    expect(dropped).toEqual([firstTok, maskTok]);
  });
});

describe('tool runs that read saved counts (-ustat from a saved copy)', () => {
  it('hold the taxonomy ID in shared mode from enqueue until each run ends (completed, or cancelled while queued)', async () => {
    const opfs = new MemoryOpfs();
    const finishes: Array<() => void> = [];
    const engine = stubEngine(new FakeEngine(), {
      cancel: () => undefined,
      async runTool() {
        await new Promise<void>((resolve) => finishes.push(resolve));
        return { status: 0, stdout: new Blob(), stderr: new Blob(), files: [] };
      },
    });
    const { coordinator } = makeCoordinator(engine, undefined, webCountsLocks(opfs));
    const held = async () => {
      await sleep(0);
      return (await opfs.query()).held.filter((entry) => entry.name === countsLockName(562) && entry.mode === 'shared').length;
    };
    const masking = (): ToolRequest => ({
      tool: 'windowmasker',
      stage: 'mask',
      argv: ['windowmasker', '-ustat', 'wmasker.obinary', '-in', 'a.fa', '-out', 'a.out'],
      directory: [
        { name: 'wmasker.obinary', file: new Blob([new Uint8Array(4)]), origin: { kind: 'saved', taxid: 562 } },
        { name: 'a.fa', file: new Blob(['>a\nACGT\n']), origin: { kind: 'picked' } },
      ],
    });
    const first = (await coordinator.enqueueTool(masking())).runId!;
    const second = (await coordinator.enqueueTool(masking())).runId!;
    expect(await held()).toBe(2);
    coordinator.cancel(second);
    expect(await held()).toBe(1);
    for (let i = 0; i < 100 && finishes.length === 0; i++) await sleep(1);
    finishes[0]!();
    expect((await waitTerminal(coordinator, first)).status).toBe('completed');
    expect(await held()).toBe(0);
    // A tool run that reads no saved copy holds nothing.
    const plain = (await coordinator.enqueueTool(sampleMkCountsRequest())).runId!;
    expect(await held()).toBe(0);
    for (let i = 0; i < 100 && finishes.length < 2; i++) await sleep(1);
    finishes[1]!();
    await waitTerminal(coordinator, plain);
  });

  const expectValid = (draft: ToolDraft) => expect(draft.state.get().validation).toEqual({ state: 'ok' });
  const makeTrackedDraft = (
    engine: EngineGateway & ToolEngine,
    data: DataService,
    coordinator: Coordinator,
    readSavedCounts?: (taxid: number) => Promise<File | undefined>,
  ) => {
    const calls: ToolRequest[][] = [];
    const draft = new ToolDraft({
      engine,
      data,
      enqueueToolAll: (requests) => {
        calls.push([...requests]);
        return coordinator.enqueueToolAll(requests);
      },
      ...(readSavedCounts === undefined ? {} : { readSavedCounts }),
      debounceMs: 0,
    });
    return { draft, calls };
  };
  const submitAndComplete = async (draft: ToolDraft, coordinator: Coordinator) => {
    const res = await draft.submit();
    expect(res.ok).toBe(true);
    expectValid(draft);
    expect(draft.state.get().message).toBeUndefined();
    const run = await waitTerminal(coordinator, res.runIds![0]!);
    expect(run.status).toBe('completed');
    return run;
  };
  const readyDraft = async (draft: ToolDraft, stage: 'mk_counts' | 'mask', taxid?: number, size?: number) => {
    draft.ensureLoaded();
    await sleep(20);
    if (stage === 'mask') {
      draft.setStage('mask');
      if (taxid !== undefined && size !== undefined) draft.setUstatSaved(taxid, size);
    }
    draft.addFastaFiles([new File(['>a\nACGT\n'], 'a.fa')]);
    for (let i = 0; i < 40 && draft.state.get().fastaFiles.some((f) => f.status !== 'ready'); i++) await sleep(10);
    await draft.runValidationNow();
    expectValid(draft);
  };

  it('keeps validation ok on a saved-counts lock refusal so the next submit queues the run', async () => {
    const opfs = new MemoryOpfs();
    const counts = new OpfsCountsStore(opfs.root as unknown as FileSystemDirectoryHandle, opfs, undefined, () => Date.UTC(2026, 10, 25));
    const made: CountsMadeMeta = { run: 1, argv: ['windowmasker', '-mk_counts'], inputs: [], engineBuild: 'test', madeAt: 10_000 };
    const initTok = (await counts.beginCounts(562, { made, size: 6, from: 0, replace: true }))!;
    await counts.appendCounts(initTok, new TextEncoder().encode('counts'));
    await counts.finishCounts(initTok);

    const engine = new FakeEngine();
    const { coordinator, data } = makeCoordinator(engine, undefined, webCountsLocks(opfs));
    let writeTok: string | undefined;
    const { draft, calls } = makeTrackedDraft(engine, data, coordinator, async (taxid) => {
      const file = (await counts.readCounts(taxid))?.file;
      // Simulate a writer starting after resolveCountsChoice read the file, before holdCounts runs.
      if (calls.length === 0 && writeTok === undefined) {
        writeTok = (await counts.beginCounts(562, { made, size: 6, from: 0, replace: true }))!;
      }
      return file;
    });
    await readyDraft(draft, 'mask', 562, 6);

    const first = await draft.submit();
    const busyText = 'The saved counts for 562 are being written or removed. Run again when that ends.';
    expect(first).toEqual({ ok: false, message: busyText });
    expect(calls).toHaveLength(1);
    expectValid(draft);
    expect(draft.state.get().message).toEqual({ kind: 'error', text: busyText });

    // Release the exclusive lock and click Run LOSAT again without changing any option.
    await counts.appendCounts(writeTok!, new TextEncoder().encode('counts'));
    await counts.finishCounts(writeTok!);
    await submitAndComplete(draft, coordinator);
    expect(calls).toHaveLength(2);
  });

  it('still records an argv refusal from enqueueToolAll in validation', async () => {
    const base = new FakeEngine();
    let refuseOnEnqueue = false;
    const engine = stubEngine(base, {
      validateTool: async (argv) =>
        refuseOnEnqueue
          ? { ok: false, message: 'Error: Invalid argument "unit"' }
          : base.validateTool(argv),
    });
    const { coordinator, data } = makeCoordinator(engine);
    const { draft, calls } = makeTrackedDraft(engine, data, coordinator);
    await readyDraft(draft, 'mk_counts');

    refuseOnEnqueue = true;
    const first = await draft.submit();
    expect(first).toEqual({ ok: false, message: 'Error: Invalid argument "unit"' });
    expect(calls).toHaveLength(1);
    expect(draft.state.get().validation).toEqual({ state: 'invalid', message: 'Error: Invalid argument "unit"' });
    expect(draft.state.get().message).toEqual({ kind: 'error', text: 'Error: Invalid argument "unit"' });

    // Without changing any option, a second submit stops at the recorded argv validation.
    const second = await draft.submit();
    expect(second).toEqual({ ok: false, message: 'Error: Invalid argument "unit"' });
    expect(calls).toHaveLength(1);
  });

  it('reads the new saved copy on Run again after a Replace invalidated the earlier handle', async () => {
    const opfs = new MemoryOpfs();
    const store = new OpfsCountsStore(opfs.root as unknown as FileSystemDirectoryHandle, opfs, undefined, () => Date.UTC(2026, 10, 25));
    const loader = new CountsLoader({
      downloader: new FakeCountsDownloader(store, {}),
      store,
      draft: { setCounts: () => undefined },
    });
    const made1: CountsMadeMeta = { run: 1, argv: ['windowmasker', '-mk_counts'], inputs: [], engineBuild: 'test', madeAt: 10_000 };
    const made2: CountsMadeMeta = { run: 2, argv: ['windowmasker', '-mk_counts'], inputs: [], engineBuild: 'test', madeAt: 20_000 };
    expect((await loader.saveMade(562, new File(['v1-counts'], 'v1.obinary'), made1)).ok).toBe(true);

    const seenCounts: string[] = [];
    const base = new FakeEngine();
    const engine = stubEngine(base, {
      async runTool(argv, opts) {
        const entry = opts?.directory?.find((f) => f.name === 'windowmasker/562/wmasker.obinary');
        if (entry !== undefined) seenCounts.push(await entry.file.text());
        return base.runTool(argv, opts);
      },
    });
    const locks = webCountsLocks(opfs);
    const { coordinator, data } = makeCoordinator(engine, undefined, locks);

    let readCount = 0;
    const { draft } = makeTrackedDraft(engine, data, coordinator, async (taxid) => {
      readCount++;
      const file = (await loader.readSavedFile(taxid))?.file;
      if (readCount === 1 && file !== undefined) {
        // Replace the stored copy after resolveCountsChoice opened v1; simulate browser invalidating the old handle.
        expect((await loader.saveMade(taxid, new File(['v2-replaced-counts!'], 'v2.obinary'), made2, { replace: true })).ok).toBe(true);
        Object.defineProperty(file, 'slice', { value: () => new Blob([]) });
      }
      return file;
    });
    await readyDraft(draft, 'mask', 562, 9);

    const first = await draft.submit();
    const replacedText = 'The saved counts for 562 were replaced or removed after they were chosen. Run again.';
    expect(first).toEqual({ ok: false, message: replacedText });
    expectValid(draft);
    expect(draft.state.get().message).toEqual({ kind: 'error', text: replacedText });

    // Clicking Run LOSAT again without re-choosing reads the replacement copy from storage and succeeds.
    const run = await submitAndComplete(draft, coordinator);
    expect(seenCounts).toEqual(['v2-replaced-counts!']);
    if (isToolRun(run)) {
      const savedEntry = run.snapshot.directory.find((entry) => entry.name === 'windowmasker/562/wmasker.obinary');
      expect(savedEntry?.size).toBe('v2-replaced-counts!'.length);
    }

    // If the saved copy is later removed, Run again reports it cleanly without marking argv invalid.
    await loader.delete(562);
    const third = await draft.submit();
    expect(third).toEqual({
      ok: false,
      message: 'The saved counts for taxonomy ID 562 are no longer available.',
    });
    expectValid(draft);
  });
});
