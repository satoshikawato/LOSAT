// The harness of the real engine (S09): searches through the application (V-BR: the
// Coordinator, the Data worker, the Wasm EngineGateway and its Engine worker), the S10
// contract suites against the real implementations, and the runtime checks (cancel,
// renewal, the serial fallback).
import { ENGINE_ASSETS } from 'virtual:losat-engine';
import type { RunView } from '../../../src/application/coordinator';
import { createApp, type App } from '../../../src/composition';
import { buildArgv } from '../../../src/domain/argv';
import type { OutputFormat } from '../../../src/domain/output-format';
import { indexParser, type InputRole, type ProgramId } from '../../../src/domain/programs';
import { isTerminal, type RunRecord, type RunStatus } from '../../../src/domain/run';
import { sha256Hex } from '../../../src/infra/browser/platform';
import { startDataWorker } from '../../../src/infra/data-worker/gateway';
import { WasmEngine } from '../../../src/infra/engine-worker/gateway';
import type { RenewalLimits } from '../../../src/infra/engine-worker/policy';
import type { EngineEvent, WorkerTestCommand } from '../../../src/infra/engine-worker/protocol';
import type { EngineInput, RuntimeInfo } from '../../../src/ports/engine';
import type { OutputStream } from '../../../src/ports/run-output';
import { runCases, type CaseResult } from '../../contract/contract';
import { ENGINE_INPUT_CASES } from '../../contract/engine-input.contract';
import { RUN_OUTPUT_CASES, type RemoteWriter } from '../../contract/run-output.contract';

export interface SearchInput {
  /** Where the harness fetches the bytes (the test serves them). */
  readonly url: string;
  /** The -query / -subject name, as the File name (plan §5.3). */
  readonly name: string;
}

export interface SearchCase {
  readonly id: string;
  readonly program: ProgramId;
  /** The options after `-query <name> -subject <name>`, as CLI words. */
  readonly options: readonly string[];
  readonly query: SearchInput;
  readonly subject: SearchInput;
  readonly threads: number | 'auto';
}

export interface SearchResult {
  readonly id: string;
  readonly status: RunStatus;
  readonly record: RunRecord;
  readonly argv?: readonly string[];
  /** SHA-256 and length of each stored output format. */
  readonly sha256: Partial<Record<OutputFormat, string>>;
  readonly lengths: Partial<Record<OutputFormat, number>>;
  readonly diagnosticsSha256?: string;
  readonly hits?: number;
  /**
   * The HSP records are the outfmt 6 rows, in order: as many records as rows, and each
   * record's `out6` range is the row at its index; their `out0` and `out0_subject` ranges
   * lie inside the outfmt 0 text.
   */
  readonly hitsMatchRows?: boolean;
  /** Milliseconds from Add to queue to the end of the run. */
  readonly elapsedMs: number;
}

export interface SessionOptions {
  readonly renewal?: Partial<RenewalLimits>;
}

const FORMATS: readonly OutputFormat[] = [0, 6, 7];
const encoder = new TextEncoder();

/** Splits CLI option words into the (flag, value) parameters of a search request. */
export function parameters(options: readonly string[]): Array<readonly [string, string | true]> {
  const result: Array<readonly [string, string | true]> = [];
  for (let i = 0; i < options.length; i++) {
    const flag = options[i]!;
    const next = options[i + 1];
    if (next !== undefined && (!next.startsWith('-') || /^-\d/.test(next))) {
      result.push([flag, next]);
      i++;
    } else {
      result.push([flag, true]);
    }
  }
  return result;
}

async function fetchFile(input: SearchInput): Promise<File> {
  const response = await fetch(input.url);
  if (!response.ok) throw new Error(`${input.url}: HTTP ${response.status}`);
  return new File([await response.arrayBuffer()], input.name);
}

/** One application instance whose exports the harness captures. */
export class SearchSession {
  readonly app: App;
  private readonly exported = new Map<string, Uint8Array>();

  constructor(options: SessionOptions = {}) {
    if (ENGINE_ASSETS === null) throw new Error('the harness was built without the engine (LOSAT_WEB_REACTORS)');
    this.app = createApp({
      downloader: { save: (fileName, bytes) => this.exported.set(fileName, bytes) },
      ...(options.renewal === undefined ? {} : { renewal: options.renewal }),
    });
  }

  /** Queues a search and returns its run ID (or throws the validation message). */
  async enqueue(search: SearchCase): Promise<string> {
    const [query, subject] = await Promise.all([fetchFile(search.query), fetchFile(search.subject)]);
    const result = await this.app.coordinator.enqueue({
      program: search.program,
      query: { file: query },
      subject: { file: subject },
      parameters: parameters(search.options),
      requestedThreads: search.threads,
    });
    if (!result.ok) throw new Error(result.message);
    return result.runId!;
  }

  view(runId: string): RunView {
    const view = this.app.coordinator.state.get().runs.find((run) => run.snapshot.runId === runId);
    if (view === undefined) throw new Error(`no run ${runId}`);
    return view;
  }

  /** Resolves when the run reaches `status` (or any terminal status). */
  until(runId: string, statuses: readonly RunStatus[]): Promise<RunView> {
    return new Promise((resolve) => {
      const check = () => {
        const view = this.view(runId);
        if (statuses.includes(view.status) || isTerminal(view.status)) {
          unsubscribe();
          resolve(view);
        }
      };
      const unsubscribe = this.app.coordinator.state.subscribe(check);
      check();
    });
  }

  /** Runs one search to its end and describes what the application stored. */
  async search(search: SearchCase): Promise<SearchResult> {
    const started = performance.now();
    let runId: string;
    try {
      runId = await this.enqueue(search);
    } catch (error) {
      return { id: search.id, status: 'failed', record: { error: String(error) }, sha256: {}, lengths: {}, elapsedMs: 0 };
    }
    const view = await this.until(runId, []);
    const elapsedMs = performance.now() - started;
    return this.describe(search.id, view, elapsedMs);
  }

  async describe(id: string, view: RunView, elapsedMs: number): Promise<SearchResult> {
    const base = { id, status: view.status, record: view.record, argv: view.snapshot.argv, elapsedMs };
    if (view.status !== 'completed') return { ...base, sha256: {}, lengths: {} };
    const runId = view.snapshot.runId;
    const sha256: Partial<Record<OutputFormat, string>> = {};
    const lengths: Partial<Record<OutputFormat, number>> = {};
    const outputs: Partial<Record<OutputFormat, Uint8Array>> = {};
    for (const format of FORMATS) {
      await this.app.coordinator.exportOutput(runId, format);
      const bytes = [...this.exported.values()].pop()!;
      this.exported.clear();
      outputs[format] = bytes;
      sha256[format] = await sha256Hex(bytes);
      lengths[format] = bytes.length;
    }
    const diagnostics = await this.app.coordinator.readDiagnostics(runId);
    const hits = await this.app.coordinator.readHits(runId);
    const out6 = outputs[6]!;
    const out0Length = outputs[0]!.length;
    const rows = new TextDecoder().decode(out6).split('\n').filter((line) => line !== '');
    const inOut0 = (range: readonly [number, number] | null) => range === null || (0 <= range[0] && range[0] < range[1] && range[1] <= out0Length);
    const hitsMatchRows =
      hits.length === rows.length &&
      hits.every((hit, i) => {
        if (hit.out6 === null) return false;
        const row = new TextDecoder().decode(out6.subarray(hit.out6[0], hit.out6[1])).replace(/\n$/, '');
        return hit.index === i && row === rows[i] && inOut0(hit.out0) && inOut0(hit.out0_subject);
      });
    return {
      ...base,
      sha256,
      lengths,
      diagnosticsSha256: await sha256Hex(encoder.encode(diagnostics)),
      hits: hits.length,
      hitsMatchRows,
    };
  }
}

let shared: SearchSession | undefined;

/** Searches in order in one application (the subject is retained between them, R1). */
export async function searches(list: readonly SearchCase[], options: SessionOptions = {}): Promise<SearchResult[]> {
  const session = options.renewal === undefined ? (shared ??= new SearchSession()) : new SearchSession(options);
  const results: SearchResult[] = [];
  for (const search of list) {
    const result = await session.search(search);
    console.log(`${result.status} ${search.id} n=${String(search.threads)} ${Math.round(result.elapsedMs)} ms`);
    results.push(result);
  }
  return results;
}

export interface CancelResult {
  readonly cancelled: SearchResult;
  readonly next: SearchResult;
  /** Milliseconds from the cancel to the cancelled state. */
  readonly cancelMs: number;
  /** Milliseconds of the next search from its start to its search phase (worker, instance, inputs). */
  readonly nextPrepareMs: number;
}

/** Cancels `long` once it runs, then runs `next` in the same application. */
export async function cancelThenSearch(long: SearchCase, next: SearchCase, afterMs = 300): Promise<CancelResult> {
  const session = new SearchSession();
  const runId = await session.enqueue(long);
  await session.until(runId, ['running']);
  await new Promise((resolve) => setTimeout(resolve, afterMs));
  const cancelAt = performance.now();
  session.app.coordinator.cancel(runId);
  const cancelledView = await session.until(runId, ['cancelled']);
  const cancelMs = performance.now() - cancelAt;
  const cancelled = await session.describe(long.id, cancelledView, 0);
  const result = await session.search(next);
  const phase = result.record.phaseTimes;
  const nextPrepareMs = (phase?.running ?? 0) - (result.record.startedAt ?? 0);
  return { cancelled, next: result, cancelMs, nextPrepareMs };
}

export interface RetentionResult {
  readonly id: string;
  readonly argv: readonly string[];
  readonly runtime?: RuntimeInfo;
  readonly error?: string;
  readonly sha256: Partial<Record<OutputFormat, string>>;
  readonly diagnosticsSha256?: string;
}

/**
 * Plan §4.6 (R1): searches in order through one EngineGateway, each input indexed once, so
 * that searches with the same subject dataset revision reuse the subject that the engine
 * holds. The queries and options change between the searches.
 */
export async function retention(
  list: readonly SearchCase[],
  options: { readonly renewal?: Partial<RenewalLimits> } = {},
): Promise<RetentionResult[]> {
  if (ENGINE_ASSETS === null) throw new Error('the harness was built without the engine (LOSAT_WEB_REACTORS)');
  const data = startDataWorker();
  const engine = new WasmEngine({
    assets: ENGINE_ASSETS,
    control: data,
    ...(options.renewal === undefined ? {} : { renewal: options.renewal }),
  });
  const inputs = new Map<string, Promise<EngineInput>>();
  const input = (search: SearchCase, role: InputRole) => {
    const which = search[role];
    // A file read as the other sequence kind is another record table.
    const parser = indexParser(search.program, role);
    const key = JSON.stringify([which.url, which.name, parser]);
    let cached = inputs.get(key);
    if (cached === undefined) {
      cached = (async () => {
        const source = await data.addSource(await fetchFile(which));
        const revision = await data.indexSource(source.sourceId, parser);
        const run = await data.buildRunInput([revision.revisionId]);
        return { bytes: run.bytes, sha256: run.sha256, revisionIds: [revision.revisionId], records: run.records };
      })();
      inputs.set(key, cached);
    }
    return cached;
  };
  const results: RetentionResult[] = [];
  for (const [i, search] of list.entries()) {
    const runId = `retention-${i}`;
    const argv = buildArgv({
      program: search.program,
      queryName: search.query.name,
      subjectName: search.subject.name,
      parameters: parameters(search.options),
    });
    const output = await data.openRun(runId);
    try {
      const runtime = await engine.run(
        { runId, argv, query: await input(search, 'query'), subject: await input(search, 'subject'), requestedThreads: search.threads },
        output,
        () => undefined,
      );
      await data.commitRun(runId);
      const sha256: Partial<Record<OutputFormat, string>> = {};
      for (const format of FORMATS) sha256[format] = await sha256Hex(await data.readOutput(runId, format));
      const diagnosticsSha256 = await sha256Hex(encoder.encode(await data.readDiagnostics(runId)));
      console.log(`${search.id}: subject ${runtime.subjectRetained ? 'retained' : 'registered'}`);
      results.push({ id: search.id, argv, runtime, sha256, diagnosticsSha256 });
    } catch (error) {
      await data.discardRun(runId);
      results.push({ id: search.id, argv, error: String(error), sha256: {} });
    }
  }
  return results;
}

export interface TimedRun {
  readonly id: string;
  readonly threads: number | 'auto';
  readonly error?: string;
  readonly runtime?: RuntimeInfo;
  /** From the request to the engine's search phase: worker and instance start, registration. */
  readonly prepareMs?: number;
  /** The engine's search phase (losat_web2_run) until its outputs are sent. */
  readonly runMs?: number;
  /** Until the outputs are committed in the data layer. */
  readonly totalMs: number;
  /** Bytes of the committed outfmt 0, 6 and 7. */
  readonly outputBytes?: Partial<Record<OutputFormat, number>>;
  /** From the cancel to the rejection of the search. */
  readonly cancelMs?: number;
}

export interface TimedOptions {
  readonly renewal?: Partial<RenewalLimits>;
  /** Cancels the search at this index once its search phase has run this long. */
  readonly cancel?: { readonly index: number; readonly afterMs: number };
  /**
   * Indexes the inputs again for every search, as the application does today (each queued
   * search makes new dataset revisions), so that no subject is retained.
   */
  readonly reindex?: boolean;
}

/**
 * Measurements (S09): searches in order through one EngineGateway and Data worker, each
 * input indexed once (so a search of the same subject reuses it, R1), with the time of each
 * phase. Committed runs are deleted after their sizes are read.
 */
export async function timed(list: readonly SearchCase[], options: TimedOptions = {}): Promise<TimedRun[]> {
  if (ENGINE_ASSETS === null) throw new Error('the harness was built without the engine (LOSAT_WEB_REACTORS)');
  const data = startDataWorker();
  const engine = new WasmEngine({
    assets: ENGINE_ASSETS,
    control: data,
    ...(options.renewal === undefined ? {} : { renewal: options.renewal }),
  });
  const inputs = new Map<string, Promise<EngineInput>>();
  const input = (search: SearchCase, role: InputRole) => {
    const which = search[role];
    // A file read as the other sequence kind is another record table.
    const parser = indexParser(search.program, role);
    const key = JSON.stringify([which.url, which.name, parser]);
    let cached = options.reindex === true ? undefined : inputs.get(key);
    if (cached === undefined) {
      cached = (async () => {
        const source = await data.addSource(await fetchFile(which));
        const revision = await data.indexSource(source.sourceId, parser);
        const run = await data.buildRunInput([revision.revisionId]);
        return { bytes: run.bytes, sha256: run.sha256, revisionIds: [revision.revisionId], records: run.records };
      })();
      inputs.set(key, cached);
    }
    return cached;
  };
  const results: TimedRun[] = [];
  for (const [i, search] of list.entries()) {
    const runId = `timed-${i}`;
    const argv = buildArgv({ program: search.program, queryName: search.query.name, subjectName: search.subject.name, parameters: parameters(search.options) });
    const query = await input(search, 'query');
    const subject = await input(search, 'subject');
    const output = await data.openRun(runId);
    const phases: Partial<Record<string, number>> = {};
    let timer: ReturnType<typeof setTimeout> | undefined;
    let cancelledAt: number | undefined;
    const started = performance.now();
    try {
      const runtime = await engine.run({ runId, argv, query, subject, requestedThreads: search.threads }, output, (phase) => {
        phases[phase] = performance.now();
        if (phase === 'running' && options.cancel?.index === i) {
          timer = setTimeout(() => {
            cancelledAt = performance.now();
            engine.cancel(runId);
          }, options.cancel.afterMs);
        }
      });
      const ran = performance.now();
      await data.commitRun(runId);
      const outputBytes: Partial<Record<OutputFormat, number>> = {};
      for (const format of FORMATS) outputBytes[format] = (await data.readOutput(runId, format)).length;
      await data.deleteRun(runId);
      results.push({
        id: search.id,
        threads: search.threads,
        runtime,
        prepareMs: (phases['running'] ?? ran) - started,
        runMs: ran - (phases['running'] ?? started),
        totalMs: performance.now() - started,
        outputBytes,
      });
    } catch (error) {
      clearTimeout(timer);
      const stopped = performance.now();
      await data.discardRun(runId);
      results.push({
        id: search.id,
        threads: search.threads,
        error: String(error),
        totalMs: stopped - started,
        ...(cancelledAt === undefined ? {} : { cancelMs: stopped - cancelledAt }),
      });
    }
    const last = results[results.length - 1]!;
    console.log(`${last.error === undefined ? 'done' : 'stopped'} ${search.id} n=${String(search.threads)}: prepare ${Math.round(last.prepareMs ?? 0)} ms, run ${Math.round(last.runMs ?? 0)} ms`);
  }
  return results;
}

/** The engine-input contract against the real EngineGateway and Engine worker. */
export async function engineInput(env: {
  readonly argv: readonly string[];
  readonly query: string;
  readonly subject: string;
  readonly queryRecords: ReadonlyArray<{ readonly id: string; readonly length: number }>;
  readonly subjectRecords: ReadonlyArray<{ readonly id: string; readonly length: number }>;
  readonly threads: number;
}): Promise<CaseResult[]> {
  if (ENGINE_ASSETS === null) throw new Error('the harness was built without the engine (LOSAT_WEB_REACTORS)');
  const engine = new WasmEngine({ assets: ENGINE_ASSETS, control: startDataWorker() });
  return runCases(
    ENGINE_INPUT_CASES,
    () => ({
      engine,
      argv: env.argv,
      query: encoder.encode(env.query),
      subject: encoder.encode(env.subject),
      queryRecords: env.queryRecords,
      subjectRecords: env.subjectRecords,
      threads: env.threads,
    }),
    { timeoutMs: 120_000, onResult: (result) => console.log(`${result.ok ? 'pass' : 'FAIL'} [engine input] ${result.name}${result.error ? `: ${result.error}` : ''}`) },
  );
}

/** The run output contract with the writer in the real Engine worker (test hooks). */
export async function runOutputEngine(): Promise<{ backend: string; results: CaseResult[] }> {
  const data = startDataWorker();
  const worker = new Worker(new URL('../../../src/infra/engine-worker/engine-worker.ts', import.meta.url), {
    type: 'module',
    name: 'losat-engine',
  });
  let nextId = 0;
  const pending = new Map<number, { resolve(value: number | undefined): void; reject(error: Error): void }>();
  worker.onmessage = (event: MessageEvent<EngineEvent>) => {
    const reply = event.data;
    if (reply.type !== 'test-reply') return;
    const call = pending.get(reply.id);
    pending.delete(reply.id);
    if (reply.ok) {
      call?.resolve(reply.value);
    } else {
      const error = new Error(reply.error?.message);
      error.name = reply.error?.name ?? 'Error';
      call?.reject(error);
    }
  };
  type Command = WorkerTestCommand extends infer T ? (T extends unknown ? Omit<T, 'id'> : never) : never;
  const call = (command: Command, transfer: Transferable[] = []) =>
    new Promise<number | undefined>((resolve, reject) => {
      const id = ++nextId;
      pending.set(id, { resolve, reject });
      worker.postMessage({ ...command, id }, transfer);
    });
  const writer = async (port: MessagePort): Promise<RemoteWriter> => {
    const handle = (await call({ type: 'test-open', port }, [port]))!;
    return {
      write: async (stream: OutputStream, bytes: Uint8Array) => void (await call({ type: 'test-write', writer: handle, stream, bytes })),
      end: async () => void (await call({ type: 'test-end', writer: handle })),
      post: async (message: unknown) => void (await call({ type: 'test-post', writer: handle, message })),
    };
  };
  const quota = (bytes: number | null) => window.setStorageQuota!(bytes);
  try {
    const backend = (await data.storageInfo()).backend;
    const results = await runCases(
      RUN_OUTPUT_CASES,
      () => ({
        data,
        writer,
        usage: async () => (await data.storageInfo()).sessionBytes,
        exhaust: () => quota(0),
        restore: () => quota(null),
      }),
      { onResult: (result) => console.log(`${result.ok ? 'pass' : 'FAIL'} [run output, Engine worker] ${result.name}${result.error ? `: ${result.error}` : ''}`) },
    );
    return { backend, results };
  } finally {
    worker.terminate();
  }
}

/** The record scanner contract against the serial reactor, in a dedicated worker. */
export function recordScanner(): Promise<CaseResult[]> {
  const worker = new Worker(new URL('./scanner-worker.ts', import.meta.url), { type: 'module' });
  return new Promise((resolve, reject) => {
    worker.onerror = (event) => reject(new Error(`the scanner worker failed: ${event.message}`));
    worker.onmessage = (event: MessageEvent<{ type: string; result?: CaseResult; results?: CaseResult[] }>) => {
      if (event.data.type === 'progress') {
        const result = event.data.result!;
        console.log(`${result.ok ? 'pass' : 'FAIL'} [record scanner] ${result.name}${result.error ? `: ${result.error}` : ''}`);
        return;
      }
      worker.terminate();
      resolve(event.data.results!);
    };
    worker.postMessage('start');
  });
}

/** What the browser offers for the threaded module. */
export function capabilities(): Record<string, unknown> {
  return {
    crossOriginIsolated: globalThis.crossOriginIsolated,
    sharedArrayBuffer: typeof SharedArrayBuffer,
    hardwareConcurrency: navigator.hardwareConcurrency,
    userAgent: navigator.userAgent,
    engine: ENGINE_ASSETS,
  };
}
