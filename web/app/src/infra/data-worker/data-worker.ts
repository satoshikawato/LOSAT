// The Data worker (plan §3.1): it owns the sources, the record tables, the run registry
// and the temporary storage of one working session, and it outlives Engine workers. It also
// keeps the saved copies of Example sets (OPFS examples/, line WF-2).
// Only this worker uses OPFS synchronous access handles. It holds its own instance of the
// serial reactor, separate from the Engine worker's, for the index scan and for the
// engine's describe/validate (src/infra/reactor/control.ts). This module is the worker's
// composition root.
import { ENGINE_ASSETS } from 'virtual:losat-engine';
import type { ExampleStore } from '../../ports/example';
import type { InputChecker } from '../../ports/input-check';
import type { RecordScanner } from '../../ports/scan';
import { sha256Hex } from '../browser/platform';
import { DataService } from '../data/data-service';
import { opfsExampleStore } from '../data/example-store';
import { opfsAccess } from '../data/opfs-block-store';
import { startDataSession, type SessionLocks } from '../data/session';
import { FakeInputChecker, FakeScanner } from '../fake/fake-fasta';
import type { ReactorAbi } from '../reactor/abi';
import { ReactorInputChecker } from '../reactor/checker';
import { reactorControl, type EngineControl } from '../reactor/control';
import { instantiateSerial, loadReactor } from '../reactor/instance';
import { ReactorScanner, reopening } from '../reactor/scanner';
import { DATA_GATEWAY_METHODS, DATA_WORKER_METHODS, EXAMPLE_STORE_METHODS, type DataWorkerApi } from './methods';
import { serveRpc, type RpcEndpoint } from './rpc';

const token = crypto.randomUUID();
const locks = 'locks' in navigator ? (navigator.locks as unknown as SessionLocks) : undefined;

let scanner: RecordScanner;
let checker: InputChecker;
let control: EngineControl;
if (ENGINE_ASSETS === null) {
  // A build without the engine (build/reactors.ts): the development build's fakes.
  scanner = new FakeScanner();
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
  const reactor: () => Promise<ReactorAbi> = reopening(async () => (await instantiateSerial(await module())).abi);
  scanner = new ReactorScanner(reactor);
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

const endpoint = Promise.all([service, examples]).then(
  ([data, saved]): DataWorkerApi => ({
    ...delegate(data, DATA_GATEWAY_METHODS),
    ...delegate(saved, EXAMPLE_STORE_METHODS),
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
