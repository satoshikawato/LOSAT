// The Engine worker (plan §3.1): it runs the searches of one runtime generation. The main
// thread ends it to cancel a search or to renew the instance; it keeps no working state
// that cannot be rebuilt from the Data worker (design D01). This module is the worker's
// composition root.
import { InputMismatchError } from '../../ports/engine';
import type { OutputStream } from '../../ports/run-output';
import { EngineStoppedError } from '../reactor/abi';
import { compileReactor, fetchReactor } from '../reactor/instance';
import { RunOutputWriter } from '../run-output/writer';
import type { EngineCommand, EngineEvent, EngineInit, ModuleSource, WorkerError, WorkerRun, WorkerTestCommand } from './protocol';
import { EngineRuntime } from './runtime';

const scope = self as unknown as {
  onmessage: ((event: MessageEvent<EngineCommand>) => void) | null;
  postMessage(message: EngineEvent): void;
};

let generation = 0;
let runtime: EngineRuntime | undefined;

async function moduleOf(source: ModuleSource) {
  return 'module' in source ? source : compileReactor(await fetchReactor(source.asset));
}

function init(command: EngineInit): void {
  generation = command.generation;
  let serial: ReturnType<typeof moduleOf> | undefined;
  let threads: ReturnType<typeof moduleOf> | undefined;
  runtime = new EngineRuntime({
    serial: async () => (await (serial ??= moduleOf(command.modules.serial))).module,
    threads: () => (threads ??= moduleOf(command.modules.threads)),
    builds: command.builds,
    faultChannel: command.faultChannel,
  });
}

function toWorkerError(error: unknown): WorkerError {
  if (error instanceof InputMismatchError) {
    const prefix = `The ${error.role} records that the engine read differ from the record table: `;
    return { name: error.name, message: error.message, role: error.role, detail: error.message.slice(prefix.length) };
  }
  if (error instanceof Error) {
    return { name: error.name, message: error.message, ...(error instanceof EngineStoppedError ? { stopped: true } : {}) };
  }
  return { name: 'Error', message: String(error) };
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

// Test builds drive the worker's run output writer directly, so that the run output
// contract (tests/contract/run-output.contract.ts) runs with the writer in this worker.
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
  else if (command.type === 'run') void run(command);
  else if (__LOSAT_TEST_HOOKS__) test(command);
};
