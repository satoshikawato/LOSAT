// A thread worker of the ThreadHost (thread-host.ts): it instantiates the threaded module
// on the shared memory, waits until the host starts a thread in it, runs that thread
// (`wasi_thread_start`), and then prepares itself again for the next thread.
import { stopMessage } from '../reactor/abi';
import { instantiateThread } from '../reactor/instance';
import { createWasi, type WasiContext } from '../reactor/wasi';
import {
  SLOT_FAILED,
  SLOT_PREPARING,
  SLOT_READY,
  START_ABANDONED,
  START_PENDING,
  START_TAKEN,
  type ThreadCommand,
  type ThreadFault,
  type ThreadPrepare,
  type ThreadStart,
} from './protocol';

const scope = self as unknown as {
  onmessage: ((event: MessageEvent<ThreadCommand>) => void) | null;
  postMessage(message: unknown): void;
};

let config: ThreadPrepare | undefined;
let state: Int32Array | undefined;
let faults: BroadcastChannel | undefined;
let prepared: { readonly start: (tid: number, startArg: number) => void; readonly wasi: WasiContext } | undefined;

function setState(value: number): void {
  if (state === undefined) return;
  Atomics.store(state, 0, value);
  Atomics.notify(state, 0);
}

async function prepare(): Promise<void> {
  if (config === undefined) return;
  try {
    const wasi = createWasi();
    const thread = await instantiateThread(config.module, config.memory, wasi);
    prepared = { start: thread.start, wasi };
    setState(SLOT_READY);
  } catch (error) {
    scope.postMessage({ type: 'failed', message: `a thread worker could not be prepared: ${String(error)}` });
    setState(SLOT_FAILED);
  }
}

function start(command: ThreadStart): void {
  const started = new Int32Array(command.started);
  const thread = prepared;
  prepared = undefined;
  // Take the start unless the host already stopped waiting (it then counts this worker failed).
  const outcome =
    thread === undefined
      ? Atomics.compareExchange(started, 0, START_PENDING, START_ABANDONED)
      : Atomics.compareExchange(started, 0, START_PENDING, START_TAKEN);
  Atomics.notify(started, 0);
  if (thread === undefined || outcome !== START_PENDING) return;
  try {
    thread.start(command.tid, command.startArg);
  } catch (error) {
    // The search cannot finish: the Engine worker may wait for this thread for ever, and
    // its own thread is inside the engine, so the main thread ends the runtime.
    const fault: ThreadFault = { type: 'fault', tid: command.tid, message: stopMessage(error, thread.wasi.output()) };
    faults?.postMessage(fault);
    setState(SLOT_FAILED);
    return;
  }
  setState(SLOT_PREPARING);
  void prepare();
}

scope.onmessage = (event) => {
  const command = event.data;
  if (command.type === 'prepare') {
    config = command;
    state = new Int32Array(command.slot);
    faults = new BroadcastChannel(command.faultChannel);
    void prepare();
    return;
  }
  start(command);
};
