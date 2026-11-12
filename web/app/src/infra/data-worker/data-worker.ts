// The Data worker (plan §3.1): it owns the sources, the record tables, the run registry
// and the temporary storage of one working session, and it outlives Engine workers. It also
// keeps the saved copies of Example sets (OPFS examples/, line WF-2).
// Only this worker uses OPFS synchronous access handles. It holds its own instance of the
// serial reactor, separate from the Engine worker's, for the index scan, the input check and the
// engine's describe/validate (src/infra/reactor/control.ts); the inflation of compressed sources
// and of downloaded Example sets (src/infra/reactor/inflater.ts) runs on an instance of its own
// for each file, dropped at its end. This module is the worker's composition root.
import { ENGINE_ASSETS } from 'virtual:losat-engine';
import type { CountsStore } from '../../ports/counts';
import type { ExampleStore } from '../../ports/example';
import type { FileInflater, WholeFileInflater } from '../../ports/inflate';
import type { InputChecker } from '../../ports/input-check';
import type { RecordScanner } from '../../ports/scan';
import { sha256Hex } from '../browser/platform';
import { DataService } from '../data/data-service';
import { opfsCountsStore } from '../data/counts-store';
import { opfsExampleStore } from '../data/example-store';
import { inflateFile } from '../data/inflate';
import { opfsAccess } from '../data/opfs-block-store';
import { startDataSession, type SessionLocks } from '../data/session';
import { FakeInputChecker, FakeScanner } from '../fake/fake-fasta';
import type { ReactorAbi } from '../reactor/abi';
import { ReactorInputChecker } from '../reactor/checker';
import { reactorControl, type EngineControl } from '../reactor/control';
import { ReactorInflater } from '../reactor/inflater';
import { instantiateSerial, loadReactor } from '../reactor/instance';
import { ReactorScanner, reopening } from '../reactor/scanner';
import { COUNTS_STORE_METHODS, DATA_GATEWAY_METHODS, DATA_WORKER_METHODS, EXAMPLE_STORE_METHODS, type DataWorkerApi } from './methods';
import { serveRpc, type RpcEndpoint } from './rpc';

const token = crypto.randomUUID();
const locks = 'locks' in navigator ? (navigator.locks as unknown as SessionLocks) : undefined;

let scanner: RecordScanner;
let inflater: FileInflater;
let whole: WholeFileInflater;
let checker: InputChecker;
let control: EngineControl;
if (ENGINE_ASSETS === null) {
  // A build without the engine (build/reactors.ts): the development build's fakes. Without the
  // engine's decoders every source is read as it is, whatever its name, and a downloaded Example
  // set that is compressed cannot be read.
  scanner = new FakeScanner();
  inflater = { open: async () => undefined };
  whole = { inflateFile: () => Promise.reject(new Error('this build has no engine')) };
  checker = new FakeInputChecker();
  control = {
    describe: () => Promise.reject(new Error('this build has no engine')),
    validate: () => Promise.reject(new Error('this build has no engine')),
  };
} else {
  const serial = ENGINE_ASSETS.serial;
  let compiled: Promise<WebAssembly.Module> | undefined;
  const module = () => {
    if (compiled === undefined) {
      const loading = loadReactor(serial).then((reactor) => reactor.module);
      compiled = loading;
      // A module that failed to load is loaded again by the next use.
      loading.catch(() => {
        if (compiled === loading) compiled = undefined;
      });
    }
    return compiled;
  };
  const instantiate = async (): Promise<ReactorAbi> => (await instantiateSerial(await module())).abi;
  // The shared instance serves the scanner, the checker and the control calls. An inflation takes
  // an instance of its own, from the same compiled module, and lets it go at its end: what it grew
  // in the linear memory (a zstd window of up to 128 MiB) cannot shrink, and goes with it.
  const reactor: () => Promise<ReactorAbi> = reopening(instantiate);
  scanner = new ReactorScanner(reactor);
  const reactorInflater = new ReactorInflater(reactor, instantiate);
  inflater = reactorInflater;
  // A downloaded Example set takes the path of a source (inflate.ts), on an instance of its own.
  whole = { inflateFile: (file, maxBytes) => inflateFile(file, reactorInflater, maxBytes) };
  checker = new ReactorInputChecker(reactor);
  control = reactorControl(reactor);
}

const started = startDataSession({ token, locks, opfs: opfsAccess });
const service = started.then(
  (session) =>
    new DataService({
      store: session.store,
      ...(session.fallbackReason === undefined ? {} : { fallbackReason: session.fallbackReason }),
      cleanup: session.cleanup,
      scanner,
      inflater,
      checker,
      digest: sha256Hex,
      newToken: () => crypto.randomUUID(),
      estimate: estimateStorage,
    }),
);

/**
 * Saved Example sets last only where the session's data does (OPFS, under Web Locks); elsewhere
 * nothing is kept and each tab downloads again.
 */
const examples = started.then(async (session): Promise<ExampleStore> => {
  if (session.fallbackReason === undefined && locks !== undefined) {
    try {
      return await opfsExampleStore(locks);
    } catch {
      // No OPFS for examples: nothing is kept.
    }
  }
  return {
    readExample: async () => undefined,
    saveExample: async () => undefined,
    deleteExample: async () => undefined,
  };
});

/** Saved WindowMasker counts last where the saved Example sets do; elsewhere each tab keeps its own download. */
const counts = started.then(async (session): Promise<CountsStore> => {
  if (session.fallbackReason === undefined && locks !== undefined) {
    try {
      return await opfsCountsStore(locks);
    } catch {
      // No OPFS for counts: nothing is kept.
    }
  }
  return {
    readCounts: async () => undefined,
    beginCounts: async () => undefined,
    appendCounts: async () => undefined,
    finishCounts: () => Promise.reject(new Error('this browser keeps no counts')),
    closeCounts: async () => undefined,
    deleteCounts: async () => undefined,
  };
});

const endpoint = Promise.all([service, examples, counts]).then(
  ([data, saved, savedCounts]): DataWorkerApi => ({
    ...delegate(data, DATA_GATEWAY_METHODS),
    ...delegate(saved, EXAMPLE_STORE_METHODS),
    ...delegate(savedCounts, COUNTS_STORE_METHODS),
    inflateFile: (file, maxBytes) => whole.inflateFile(file, maxBytes),
    describe: (program) => control.describe(program),
    validate: (argv) => control.validate(argv),
  }),
);

serveRpc(self as unknown as RpcEndpoint, endpoint, DATA_WORKER_METHODS);

function delegate<T extends object, K extends keyof T>(target: T, methods: readonly K[]): Pick<T, K> {
  return Object.fromEntries(
    methods.map((method) => {
      const fn = target[method] as (...args: unknown[]) => unknown;
      return [method, (...args: unknown[]) => fn.apply(target, args)];
    }),
  ) as unknown as Pick<T, K>;
}

async function estimateStorage(): Promise<{ usage: number; quota: number } | undefined> {
  if (typeof navigator.storage?.estimate !== 'function') return undefined;
  const { usage, quota } = await navigator.storage.estimate();
  return usage === undefined || quota === undefined ? undefined : { usage, quota };
}
