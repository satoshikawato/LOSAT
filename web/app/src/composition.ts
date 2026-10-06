// Composition root: the only module that chooses implementations for the ports.
import { ENGINE_ASSETS } from 'virtual:losat-engine';
import { Attention } from './application/attention';
import { Coordinator } from './application/coordinator';
import { SearchDraft } from './application/draft';
import type { Downloader } from './ports/download';
import type { EngineGateway } from './ports/engine';
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
  /** Wake lock, the warning before leaving, and the check after the page was hidden. */
  readonly attention: Attention;
  /** True while the engine is the FakeEngine; the UI shows a warning banner. */
  readonly usesFakeEngine: boolean;
}

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
  const coordinator = new Coordinator({
    engine,
    data,
    downloader: options.downloader ?? browserDownloader,
    now: () => Date.now(),
    newRunId: () => crypto.randomUUID(),
  });
  const draft = new SearchDraft({
    engine,
    data,
    enqueueAll: (requests) => coordinator.enqueueAll(requests),
    now: () => performance.now(),
  });
  const attention = new Attention({
    page: browserPage,
    runs: coordinator.state,
    probe: () => data.storageInfo(),
    now: () => Date.now(),
  });
  return { coordinator, draft, attention, usesFakeEngine: ENGINE_ASSETS === null };
}
