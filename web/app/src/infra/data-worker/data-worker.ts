// The Data worker (plan §3.1): it owns the sources, the record tables, the run registry
// and the temporary storage of one working session, and it outlives Engine workers.
// Only this worker uses OPFS synchronous access handles. It holds its own instance of the
// serial reactor, separate from the Engine worker's, for the index scan and for the
// engine's describe/validate (src/infra/reactor/control.ts). This module is the worker's
// composition root.
import { ENGINE_ASSETS } from 'virtual:losat-engine';
import type { DataGateway } from '../../ports/data';
import type { InputChecker } from '../../ports/input-check';
import type { RecordScanner } from '../../ports/scan';
import { sha256Hex } from '../browser/platform';
import { DataService } from '../data/data-service';
import { opfsAccess } from '../data/opfs-block-store';
import { startDataSession, type SessionLocks } from '../data/session';
import { FakeInputChecker, FakeScanner } from '../fake/fake-fasta';
import type { ReactorAbi } from '../reactor/abi';
import { ReactorInputChecker } from '../reactor/checker';
import { reactorControl, type EngineControl } from '../reactor/control';
import { compileReactor, fetchReactor, instantiateSerial } from '../reactor/instance';
import { ReactorScanner, reopening } from '../reactor/scanner';
import { DATA_GATEWAY_METHODS, DATA_WORKER_METHODS, type DataWorkerApi } from './methods';
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
      const loading = fetchReactor(serial)
        .then(compileReactor)
        .then((reactor) => reactor.module);
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

const service = startDataSession({ token, locks, opfs: opfsAccess }).then(
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

const endpoint = service.then(
  (data): DataWorkerApi => ({
    ...delegate(data, DATA_GATEWAY_METHODS),
    describe: (program) => control.describe(program),
    validate: (argv) => control.validate(argv),
  }),
);

serveRpc(self as unknown as RpcEndpoint, endpoint, DATA_WORKER_METHODS);

function delegate<K extends keyof DataGateway>(data: DataGateway, methods: readonly K[]): Pick<DataGateway, K> {
  return Object.fromEntries(
    methods.map((method) => {
      const fn = data[method] as (...args: unknown[]) => unknown;
      return [method, (...args: unknown[]) => fn.apply(data, args)];
    }),
  ) as unknown as Pick<DataGateway, K>;
}

async function estimateStorage(): Promise<{ usage: number; quota: number } | undefined> {
  if (typeof navigator.storage?.estimate !== 'function') return undefined;
  const { usage, quota } = await navigator.storage.estimate();
  return usage === undefined || quota === undefined ? undefined : { usage, quota };
}
