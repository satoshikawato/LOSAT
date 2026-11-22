// The Engine worker (plan §3.1): it runs the searches of one runtime generation. The main
// thread ends it to cancel a search or to renew the instance; it keeps no working state
// that cannot be rebuilt from the Data worker (design D01). This module is the worker's
// composition root.
import { InputMismatchError, runFilesOf, runOutputFailedOf } from '../../ports/engine.js';
import type { OutputStream } from '../../ports/run-output.js';
import { EngineStoppedError } from '../reactor/abi.js';
import { loadReactor } from '../reactor/instance.js';
import { RunOutputWriter } from '../run-output/writer.js';
import { holdLife, webLocks } from '../worker-life.js';
import { InstanceUse } from './instance-use.js';
import {
  HOLD_USE_CHANNEL,
  STOP_AT_OUTPUT_CHANNEL,
  type EngineCommand,
  type EngineEvent,
  type EngineInit,
  type ModuleSource,
  type WorkerError,
  type WorkerPrepare,
  type WorkerRun,
  type WorkerTestCommand,
  type WorkerTool,
  type WorkerToolValidate,
} from './protocol.js';
import { EngineRuntime } from './runtime.js';
import { toolPosts } from './tool-posts.js';

const scope = self as unknown as {
  onmessage: ((event: MessageEvent<EngineCommand>) => void) | null;
  postMessage(message: EngineEvent, transfer?: Transferable[]): void;
};

let generation = 0;
let runtime: EngineRuntime | undefined;
// Preparations and searches use the runtime one after another, in the order they came, after the worker took its lock (init).
let previous: Promise<unknown> = Promise.resolve();
// The latest preparation that waits for its turn: newer ones replace it, and a search drops
// it (the search prepares for itself).
let pending: WorkerPrepare | undefined;

function inTurn(work: () => Promise<unknown>): void {
  const turn = previous.then(work);
  previous = turn.catch(() => undefined);
}

function queuePrepare(command: WorkerPrepare): void {
  const queued = pending !== undefined;
  pending = command;
  if (queued) return;
  inTurn(async () => {
    const latest = pending;
    pending = undefined;
    if (latest !== undefined) await prepare(latest.threads, latest.argv);
  });
}

async function moduleOf(source: ModuleSource) {
  return 'module' in source ? source : loadReactor(source.asset);
}

// Test builds: a page arms the next search of this worker to stop its instance at its first
// output, as a trap in the engine would (web/app/tests/e2e/harness/engine.ts, stopDuringSearch).
let stopAtOutput = false;
if (__LOSAT_TEST_HOOKS__) {
  const arming = new BroadcastChannel(STOP_AT_OUTPUT_CHANNEL);
  arming.onmessage = (event: MessageEvent<{ type?: string }>) => {
    if (event.data?.type !== 'arm') return;
    stopAtOutput = true;
    arming.postMessage({ type: 'armed' });
  };
}

function stopIfArmed(): void {
  if (!stopAtOutput) return;
  stopAtOutput = false;
  throw new Error('a test stopped the engine at its first output');
}

// Test builds: a page arms the next use of the threaded instance to hold at its start until the
// Gateway ends the runtime, so that a test ends it during the use (web/app/tests/e2e/harness/engine.ts,
// cancelDuringWarmUp). The worker cannot take a message while it holds; the mark of the use
// count can be read (instance-use.ts).
let holdUse = false;
let prepareFailure: string | undefined;
const holding = __LOSAT_TEST_HOOKS__ ? new BroadcastChannel(HOLD_USE_CHANNEL) : undefined;
if (holding !== undefined) {
  holding.onmessage = (event: MessageEvent<{ type?: string }>) => {
    if (event.data?.type !== 'arm') return;
    holdUse = true;
    holding.postMessage({ type: 'armed' });
  };
}

function holdIfArmed(use: InstanceUse | undefined): void {
  if (!holdUse || use === undefined) return;
  holdUse = false;
  holding?.postMessage({ type: 'holding' });
  const pause = new Int32Array(new SharedArrayBuffer(4));
  // A test that does not end the runtime fails on its own timeout; the use then goes on.
  const until = performance.now() + 120_000;
  while (!use.retiring() && performance.now() < until) Atomics.wait(pause, 0, 0, 10);
}

function init(command: EngineInit): void {
  generation = command.generation;
  // Tells the Gateway at once that this worker started: an error from now on is a worker that stopped (line WQA's sweep, D1).
  scope.postMessage({ type: 'started', generation });
  // Holds the worker's lock until it ends (worker-life.ts), and tells the Gateway so before any preparation or search of
  // the worker takes its turn: the Gateway can wait for the lock only once the worker holds it. A search posted with
  // `init` (the first search of a new worker: after a cancel, before the warm-up, after a renewal) reached its
  // synchronous call into the engine before the lock was granted, so the Gateway did not wait for the lock during that
  // search, and a worker ended then left it "Searching" for good (line WQA's sweep, D2; its review, H1, 2026-11-19).
  previous = holdLife(webLocks(), command.lifeLock).then((held) => scope.postMessage({ type: 'life', generation: command.generation, held }));
  let serial: ReturnType<typeof moduleOf> | undefined;
  let threads: ReturnType<typeof moduleOf> | undefined;
  const use = command.use === undefined ? undefined : new InstanceUse(command.use);
  runtime = new EngineRuntime({
    serial: async () => (await (serial ??= moduleOf(command.modules.serial))).module,
    threads: () => (threads ??= moduleOf(command.modules.threads)),
    builds: command.builds,
    main: {
      // One end to the main thread, the other to the thread worker or the keeper (link.ts).
      link: (role) => {
        const channel = new MessageChannel();
        scope.postMessage({ type: 'link', generation: command.generation, role, port: channel.port2 }, [channel.port2]);
        return channel.port1;
      },
      fault: (fault) => scope.postMessage({ ...fault, generation: command.generation }),
    },
    use,
    ...(__LOSAT_TEST_HOOKS__ ? { beforeOutput: stopIfArmed, onUse: () => holdIfArmed(use) } : {}),
  });
}

function toWorkerError(error: unknown): WorkerError {
  if (error instanceof InputMismatchError) {
    const prefix = `The ${error.role} records that the engine read differ from the record table: `;
    return { name: error.name, message: error.message, role: error.role, detail: error.message.slice(prefix.length) };
  }
  if (error instanceof Error) {
    const runFiles = runFilesOf(error);
    return {
      name: error.name,
      message: error.message,
      ...(error instanceof EngineStoppedError ? { stopped: true } : {}),
      ...(runFiles === undefined ? {} : { runFiles }),
      ...(runOutputFailedOf(error) ? { outputFailed: true } : {}),
    };
  }
  return { name: 'Error', message: String(error) };
}

async function prepare(threads: number, argv: readonly string[] | undefined): Promise<void> {
  try {
    await runtime?.prepare(threads, argv);
  } catch (error) {
    // The search that follows starts the runtime again and reports what fails.
    if (__LOSAT_TEST_HOOKS__) prepareFailure = String(error);
  }
  // Test builds: a warm-up armed to hold that did not use the threaded instance says why, so
  // that its test fails with the reason rather than waiting for the hold (cancelDuringWarmUp).
  if (__LOSAT_TEST_HOOKS__ && holdUse && argv !== undefined) {
    holding?.postMessage({ type: 'not-held', reason: prepareFailure ?? 'the warm-up made no use of the threaded instance' });
  }
}

async function run(command: WorkerRun): Promise<void> {
  const { runId } = command;
  try {
    if (runtime === undefined) throw new Error('the Engine worker was not initialized');
    const result = await runtime.run(command, (phase) => scope.postMessage({ type: 'phase', generation, runId, phase }));
    scope.postMessage({ type: 'done', generation, runId, ok: true, result });
  } catch (error) {
    // A run that fails sends no `end`; the data layer discards what arrived.
    command.output.close();
    scope.postMessage({ type: 'done', generation, runId, ok: false, error: toWorkerError(error) });
  }
}

/**
 * A tool run: its standard output, standard error and file frames are posted to the main thread as the engine emits them,
 * each buffer transferred (D-WW3-6: nothing is held here), and the exit status last. Runs in its turn, as a search does.
 */
async function tool(command: WorkerTool): Promise<void> {
  const { runId } = command;
  try {
    if (runtime === undefined) throw new Error('the Engine worker was not initialized');
    // Each post counts its bytes in the run's credit and waits while the page has not taken enough (tool-credit.ts).
    const result = await runtime.tool(command, toolPosts(command, generation, (event, transfer) => scope.postMessage(event, transfer)));
    scope.postMessage({ type: 'tool-done', generation, runId, result });
  } catch (error) {
    scope.postMessage({ type: 'done', generation, runId, ok: false, error: toWorkerError(error) });
  }
}

/**
 * The check of a tool's words (`ToolEngine.validateTool`): on the serial instance that runs tools, in its turn after the run
 * in progress, with the run's files opened for it. The answer is the verdict, or the error of an engine that stopped.
 */
async function validateTool(command: WorkerToolValidate): Promise<void> {
  const { id } = command;
  try {
    if (runtime === undefined) throw new Error('the Engine worker was not initialized');
    const verdict = await runtime.validateTool(command);
    scope.postMessage({ type: 'tool-validated', generation, id, verdict });
  } catch (error) {
    scope.postMessage({ type: 'tool-validate-failed', generation, id, error: toWorkerError(error) });
  }
}

// Test builds drive the worker's run output writer directly, so that the run output
// contract (web/app/tests/contract/run-output.contract.ts) runs with the writer in this worker.
const testWriters = new Map<number, { readonly writer: RunOutputWriter; readonly port: MessagePort }>();

function test(command: WorkerTestCommand): void {
  try {
    if (command.type === 'test-open') {
      const id = testWriters.size + 1;
      testWriters.set(id, { writer: new RunOutputWriter(command.port), port: command.port });
      scope.postMessage({ type: 'test-reply', id: command.id, ok: true, value: id });
      return;
    }
    const entry = testWriters.get(command.writer);
    if (entry === undefined) throw new Error(`unknown writer ${command.writer}`);
    if (command.type === 'test-write') entry.writer.write(command.stream as OutputStream, command.bytes);
    else if (command.type === 'test-end') entry.writer.end();
    else entry.port.postMessage(command.message);
    scope.postMessage({ type: 'test-reply', id: command.id, ok: true });
  } catch (error) {
    scope.postMessage({ type: 'test-reply', id: command.id, ok: false, error: toWorkerError(error) });
  }
}

scope.onmessage = (event) => {
  const command = event.data;
  if (command.type === 'init') init(command);
  else if (command.type === 'run') {
    pending = undefined;
    inTurn(() => run(command));
  } else if (command.type === 'tool') {
    pending = undefined;
    inTurn(() => tool(command));
  } else if (command.type === 'tool-validate') inTurn(() => validateTool(command));
  else if (command.type === 'prepare') queuePrepare(command);
  else if (__LOSAT_TEST_HOOKS__) test(command);
};
