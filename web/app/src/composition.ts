// Composition root: the only module that chooses implementations for the ports.
import { Coordinator } from './application/coordinator';
import { browserDownloader } from './infra/browser/platform';
import { startDataWorker } from './infra/data-worker/gateway';
import { FakeEngine } from './infra/fake/fake-engine';

export interface App {
  readonly coordinator: Coordinator;
  /** True while the engine is the FakeEngine; the UI shows a warning banner. */
  readonly usesFakeEngine: boolean;
}

export function createApp(): App {
  // The Wasm engine replaces FakeEngine in W1 (docs/losat_web_gui_plan.md §7).
  const coordinator = new Coordinator({
    engine: new FakeEngine({ phaseDelayMs: 50 }),
    data: startDataWorker(),
    downloader: browserDownloader,
    now: () => Date.now(),
    newRunId: () => crypto.randomUUID(),
  });
  return { coordinator, usesFakeEngine: true };
}
