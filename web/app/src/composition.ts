// Composition root: the only module that chooses implementations for the ports.
import { Coordinator } from './application/coordinator';
import { sha256Hex, browserDownloader } from './infra/browser/platform';
import { FakeEngine } from './infra/fake/fake-engine';
import { MemoryDataGateway } from './infra/memory/memory-data-gateway';

export interface App {
  readonly coordinator: Coordinator;
  /** True while the engine is the FakeEngine; the UI shows a warning banner. */
  readonly usesFakeEngine: boolean;
}

export function createApp(): App {
  // The Wasm engine replaces FakeEngine in W1 (docs/losat_web_gui_plan.md §7).
  const coordinator = new Coordinator({
    engine: new FakeEngine({ phaseDelayMs: 50 }),
    data: new MemoryDataGateway(),
    downloader: browserDownloader,
    digest: sha256Hex,
    now: () => Date.now(),
    newRunId: () => crypto.randomUUID(),
  });
  return { coordinator, usesFakeEngine: true };
}
