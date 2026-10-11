// The browser self-test (S16; a small version of V-BR, plan §6.2): the user may run the
// examples of the programs that LOSAT Web searches (domain/examples.ts) through the
// application's engine, at the threads that Auto gives a larger search in this browser and
// at 1 thread, and compare the SHA-256 of outfmt 0, 6 and 7 of each run with the values that
// the release ships (src/examples/self-test.json). It sends nothing.
//
// It stays apart from the user's work: it runs only while no run is active or queued, with
// the engine to itself (Coordinator.withEngine: runs added meanwhile wait in the queue), and
// its inputs and outputs go to the Data worker under its own sources and run IDs that no list
// of the application knows, released and deleted after each case. So nothing of it reaches
// the queue, the results, the candidate tray or a session file.
import type { OutputFormat } from '../domain/output-format';
import { OUTPUT_FORMATS } from '../domain/output-format';
import type { SelfTestCase, SelfTestTable } from '../domain/examples';
import { indexParser, programById, type InputRole, type ProgramId } from '../domain/programs';
import { isTerminal } from '../domain/run';
import type { DataGateway, StorageBackend } from '../ports/data';
import type { EngineGateway, EngineInput, RuntimeInfo } from '../ports/engine';
import type { ExampleSource } from '../ports/examples';
import type { AppState } from './coordinator';
import { Store } from './store';

/** What the browser says of itself, on which the engine depends. */
export interface BrowserFacts {
  /** Whether the page may share memory between threads (COOP and COEP); without it searches run serially. */
  readonly crossOriginIsolated: boolean;
  /** `navigator.hardwareConcurrency`, where the browser reports it. */
  readonly logicalProcessors?: number;
}

export interface SelfTestFacts extends BrowserFacts {
  /** Where the Data worker keeps outputs, and why not OPFS when it is memory. */
  readonly storage?: { readonly backend: StorageBackend; readonly fallbackReason?: string };
}

/** One output of one case: the SHA-256 that the release expects and the one the run gave. */
export interface FormatCheck {
  readonly format: OutputFormat;
  readonly expected: string;
  readonly got: string;
}

export type SelfTestOutcome =
  | { readonly kind: 'pass' }
  | { readonly kind: 'fail'; readonly mismatches: readonly FormatCheck[] }
  | { readonly kind: 'error'; readonly message: string };

/** A case: one program's example at one thread count. */
export interface SelfTestCaseResult {
  readonly program: ProgramId;
  readonly threads: number;
  readonly state: 'waiting' | 'running' | 'done';
  readonly outcome?: SelfTestOutcome;
  /** How the engine ran it: the module, the threads, the build, and why serially if more threads were asked for. */
  readonly runtime?: Pick<RuntimeInfo, 'path' | 'threads' | 'engineBuild' | 'fallbackReason'>;
  /** Milliseconds from the start of the case to its end. */
  readonly ms?: number;
}

export interface SelfTestState {
  readonly status: 'idle' | 'running' | 'finished';
  /** Why the latest request did not start the self-test. */
  readonly refusal?: string;
  readonly cases: readonly SelfTestCaseResult[];
  /** Programs with an example that this release cannot search yet, and why (BLASTX until SX). */
  readonly skipped: ReadonlyArray<{ readonly program: ProgramId; readonly reason: string }>;
  readonly facts?: SelfTestFacts;
  /** One line on the finished self-test. */
  readonly summary?: string;
}

export interface SelfTestDeps {
  readonly engine: Pick<EngineGateway, 'run'>;
  readonly data: Pick<
    DataGateway,
    'addSource' | 'indexSource' | 'buildRunInput' | 'releaseSources' | 'openRun' | 'commitRun' | 'discardRun' | 'readOutput' | 'deleteRun' | 'storageInfo'
  >;
  readonly examples: ExampleSource;
  readonly table: SelfTestTable;
  /** The application's runs: the self-test starts only when none is active or queued. */
  readonly runs: Store<AppState>;
  /** Runs a task with the engine to itself (Coordinator.withEngine). */
  readonly withEngine: <T>(task: () => Promise<T>) => Promise<T>;
  readonly digest: (bytes: Uint8Array) => Promise<string>;
  readonly browser: () => BrowserFacts;
  /** The thread counts of each case, in order. */
  readonly threads: readonly number[];
  readonly newRunId: () => string;
  readonly now: () => number;
  /** The engine is the FakeEngine, whose outputs are not search results: the self-test cannot run. */
  readonly fake: boolean;
}

export const FAKE_ENGINE_REFUSAL =
  'The self-test needs the engine. This development build uses a fake engine, whose outputs are not search results.';

const ROLES: readonly InputRole[] = ['query', 'subject'];

export class SelfTest {
  readonly state: Store<SelfTestState>;

  constructor(private readonly deps: SelfTestDeps) {
    this.state = new Store<SelfTestState>({ status: 'idle', cases: [], skipped: this.skipped() });
  }

  /** The programs that the self-test runs, and the thread counts of each, for the text that describes it. */
  plan(): { readonly programs: readonly ProgramId[]; readonly threads: readonly number[] } {
    return { programs: this.runnable().map((c) => c.program), threads: this.deps.threads };
  }

  /** Why this build cannot run the self-test at all (the FakeEngine), or undefined. */
  unavailable(): string | undefined {
    return this.deps.fake ? FAKE_ENGINE_REFUSAL : undefined;
  }

  /** Why the self-test cannot start now, or undefined. */
  refusal(): string | undefined {
    if (this.deps.fake) return FAKE_ENGINE_REFUSAL;
    const waiting = this.deps.runs.get().runs.filter((run) => !isTerminal(run.status));
    if (waiting.length === 0) return undefined;
    const active = waiting.find((run) => run.status !== 'queued');
    const queued = waiting.length - (active === undefined ? 0 : 1);
    const parts = [];
    if (active !== undefined) parts.push(`Run ${active.snapshot.number} is running`);
    if (queued > 0) parts.push(`${queued} ${queued === 1 ? 'run is' : 'runs are'} waiting`);
    return `The self-test starts only when no search is running or waiting: ${parts.join(' and ')}. Try again when the queue is empty.`;
  }

  /** Runs every case in turn, unless it cannot start now (then `refusal` says why). */
  async run(): Promise<void> {
    if (this.state.get().status === 'running') return;
    const refusal = this.refusal();
    if (refusal !== undefined) {
      this.state.set({ ...this.state.get(), refusal });
      return;
    }
    const runnable = this.runnable();
    try {
      await this.deps.withEngine(async () => {
        this.state.set({
          status: 'running',
          cases: runnable.flatMap((c) => this.deps.threads.map((threads): SelfTestCaseResult => ({ program: c.program, threads, state: 'waiting' }))),
          skipped: this.skipped(),
          facts: this.deps.browser(),
        });
        const storage = await this.deps.data.storageInfo().then(
          (info) => ({ backend: info.backend, ...(info.fallbackReason === undefined ? {} : { fallbackReason: info.fallbackReason }) }),
          () => undefined,
        );
        if (storage !== undefined) this.state.set({ ...this.state.get(), facts: { ...this.state.get().facts!, storage } });
        for (const testCase of runnable) await this.runCase(testCase);
      });
    } catch (error) {
      // Coordinator.withEngine refused: a run was added between the check and the start.
      if (this.state.get().status !== 'running') {
        this.state.set({ ...this.state.get(), refusal: this.refusal() ?? messageOf(error) });
        return;
      }
      this.state.set({ ...this.state.get(), summary: `The self-test stopped: ${messageOf(error)}` });
    }
    this.state.set({ ...this.state.get(), status: 'finished', summary: this.state.get().summary ?? summarize(this.state.get().cases) });
  }

  private runnable(): readonly SelfTestCase[] {
    return this.deps.table.cases.filter((c) => programById(c.program).unavailable === undefined);
  }

  private skipped(): SelfTestState['skipped'] {
    return this.deps.table.cases.flatMap((c) => {
      const reason = programById(c.program).unavailable;
      return reason === undefined ? [] : [{ program: c.program, reason }];
    });
  }

  /** The program's example at each thread count, with its inputs read once and released after. */
  private async runCase(testCase: SelfTestCase): Promise<void> {
    let prepared: { readonly inputs: Readonly<Record<InputRole, EngineInput>>; readonly sources: readonly string[] };
    try {
      prepared = await this.prepare(testCase);
    } catch (error) {
      for (const threads of this.deps.threads) this.setCase(testCase.program, threads, { state: 'done', outcome: { kind: 'error', message: messageOf(error) } });
      return;
    }
    try {
      for (const threads of this.deps.threads) await this.runOne(testCase, threads, prepared.inputs);
    } finally {
      await this.deps.data.releaseSources(prepared.sources).catch(() => undefined);
    }
  }

  /** The example files as the search form would give them to the engine: sources of the Data worker, indexed for the program. */
  private async prepare(testCase: SelfTestCase) {
    const files = await this.deps.examples.files(testCase.program);
    const sources: string[] = [];
    try {
      const inputs = {} as Record<InputRole, EngineInput>;
      for (const role of ROLES) {
        const source = await this.deps.data.addSource(files[role]);
        sources.push(source.sourceId);
        const revision = await this.deps.data.indexSource(source.sourceId, indexParser(testCase.program, role));
        const input = await this.deps.data.buildRunInput([revision.revisionId]);
        // One whole file is the run input unchanged: its SHA-256 is the file's.
        if (input.sha256 !== testCase.inputs[role].sha256) {
          throw new Error(`The example file ${files[role].name} is not the one that the expected outputs were made from.`);
        }
        inputs[role] = { bytes: input.bytes, sha256: input.sha256, revisionIds: [revision.revisionId], records: input.records };
      }
      return { inputs, sources };
    } catch (error) {
      await this.deps.data.releaseSources(sources).catch(() => undefined);
      throw error;
    }
  }

  private async runOne(testCase: SelfTestCase, threads: number, inputs: Readonly<Record<InputRole, EngineInput>>): Promise<void> {
    const { program } = testCase;
    const started = this.deps.now();
    this.setCase(program, threads, { state: 'running' });
    const runId = this.deps.newRunId();
    let staged = false;
    let committed = false;
    let outcome: SelfTestOutcome;
    let runtime: SelfTestCaseResult['runtime'];
    try {
      const output = await this.deps.data.openRun(runId);
      staged = true;
      const info = await this.deps.engine.run(
        { runId, argv: testCase.argv, query: inputs.query, subject: inputs.subject, requestedThreads: threads },
        output,
        () => undefined,
      );
      runtime = {
        path: info.path,
        threads: info.threads,
        engineBuild: info.engineBuild,
        ...(info.fallbackReason === undefined ? {} : { fallbackReason: info.fallbackReason }),
      };
      await this.deps.data.commitRun(runId);
      committed = true;
      const mismatches: FormatCheck[] = [];
      for (const format of OUTPUT_FORMATS) {
        const got = await this.deps.digest(await this.deps.data.readOutput(runId, format));
        const expected = testCase.outputs[format].sha256;
        if (got !== expected) mismatches.push({ format, expected, got });
      }
      outcome = mismatches.length === 0 ? { kind: 'pass' } : { kind: 'fail', mismatches };
    } catch (error) {
      outcome = { kind: 'error', message: messageOf(error) };
    } finally {
      if (committed) await this.deps.data.deleteRun(runId).catch(() => undefined);
      else if (staged) await this.deps.data.discardRun(runId).catch(() => undefined);
    }
    this.setCase(program, threads, { state: 'done', outcome, ms: this.deps.now() - started, ...(runtime === undefined ? {} : { runtime }) });
  }

  private setCase(program: ProgramId, threads: number, change: Partial<SelfTestCaseResult>): void {
    const state = this.state.get();
    this.state.set({
      ...state,
      cases: state.cases.map((c) => (c.program === program && c.threads === threads ? { ...c, ...change } : c)),
    });
  }
}

function summarize(cases: readonly SelfTestCaseResult[]): string {
  const count = (kind: SelfTestOutcome['kind']) => cases.filter((c) => c.outcome?.kind === kind).length;
  const passed = count('pass');
  if (passed === cases.length) return `All ${cases.length} cases passed: this browser's engine wrote the expected outputs.`;
  const parts = [`${passed} of ${cases.length} cases passed`];
  if (count('fail') > 0) parts.push(`${count('fail')} wrote other outputs than expected`);
  if (count('error') > 0) parts.push(`${count('error')} could not run`);
  return `${parts.join('; ')}.`;
}

function messageOf(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
