// Composition root: the only module that chooses implementations for the ports.
import { ENGINE_ASSETS } from 'virtual:losat-engine';
import { VERIFICATION_TABLE } from 'virtual:losat-verification';
import { Attention } from './application/attention';
import { CandidateTray } from './application/candidates';
import { Coordinator } from './application/coordinator';
import { SearchDraft } from './application/draft';
import { ResultsBrowser } from './application/results';
import { ResultExporter } from './application/result-export';
import { Session } from './application/session';
import type { Downloader } from './ports/download';
import type { EngineGateway } from './ports/engine';
import { browserCompression } from './infra/browser/compression';
import { browserPage } from './infra/browser/page';
import { browserDownloader } from './infra/browser/platform';
import { startDataWorker } from './infra/data-worker/gateway';
import { WasmEngine } from './infra/engine-worker/gateway';
import type { RenewalLimits } from './infra/engine-worker/policy';
import { FakeEngine } from './infra/fake/fake-engine';

export interface App {
  readonly coordinator: Coordinator;
  /** The job being edited on the search screen. */
  readonly draft: SearchDraft;
  /** What the results screen shows (application/results.ts). */
  readonly results: ResultsBrowser;
  /** HSPs collected from the results of completed runs, and their extraction (application/candidates.ts). */
  readonly candidates: CandidateTray;
  /** LOSAT Web's own files (CSV, JSON, report) of the run that the results screen shows (application/result-export.ts). */
  readonly exporter: ResultExporter;
  /** Wake lock, the warning before leaving, and the check after the page was hidden. */
  readonly attention: Attention;
  /** True while the engine is the FakeEngine; the UI shows a warning banner. */
  readonly usesFakeEngine: boolean;
  /** Session files: saving the completed runs, opening them without a search, re-attaching originals (application/session.ts). */
  readonly session: Session;
}

/** This build's version and commit (vite.config.ts), written into session files; "unknown" where a build does not set them. */
const APP_BUILD = Object.freeze({
  version: typeof __LOSAT_APP_VERSION__ === 'string' ? __LOSAT_APP_VERSION__ : 'unknown',
  build: typeof __LOSAT_APP_BUILD__ === 'string' ? __LOSAT_APP_BUILD__ : 'unknown',
});

export interface AppOptions {
  /** Where exports go; the browser's download by default. */
  readonly downloader?: Downloader;
  /** When the engine runtime is renewed (src/infra/engine-worker/policy.ts); tests lower it. */
  readonly renewal?: Partial<RenewalLimits>;
}

/**
 * Builds the application. A build with the engine modules (build/reactors.ts) runs the
 * Wasm engine; a build without them is a development build with the FakeEngine.
 */
export function createApp(options: AppOptions = {}): App {
  const data = startDataWorker();
  const engine: EngineGateway =
    ENGINE_ASSETS === null
      ? new FakeEngine({ phaseDelayMs: 50 })
      : new WasmEngine({
          assets: ENGINE_ASSETS,
          control: data,
          ...(options.renewal === undefined ? {} : { renewal: options.renewal }),
        });
  const downloader = options.downloader ?? browserDownloader;
  const coordinator = new Coordinator({
    engine,
    data,
    downloader,
    now: () => Date.now(),
    newRunId: () => crypto.randomUUID(),
  });
  const draft = new SearchDraft({
    engine,
    data,
    enqueueAll: (requests) => coordinator.enqueueAll(requests),
    now: () => performance.now(),
  });
  const results = new ResultsBrowser({
    data,
    describe: (program) => engine.describe(program),
    runs: coordinator.state,
    verification: VERIFICATION_TABLE,
    downloader,
  });
  const candidates = new CandidateTray({ runs: coordinator.state, data, downloader, now: () => Date.now() });
  const exporter = new ResultExporter({ results: results.state, data, downloader, now: () => Date.now() });
  const attention = new Attention({
    page: browserPage,
    runs: coordinator.state,
    probe: () => data.storageInfo(),
    now: () => Date.now(),
  });
  const session = new Session({
    coordinator,
    tray: candidates,
    data,
    compression: browserCompression,
    downloader,
    app: APP_BUILD,
    now: () => Date.now(),
    newRunId: () => crypto.randomUUID(),
  });
  return { coordinator, draft, results, candidates, exporter, attention, usesFakeEngine: ENGINE_ASSETS === null, session };
}
