// Composition root: the only module that chooses implementations for the ports.
import { ENGINE_ASSETS } from 'virtual:losat-engine';
import { VERIFICATION_TABLE } from 'virtual:losat-verification';
import { engineBuildName, sha256Hex, WasmEngine, type RenewalLimits, type RuntimeEnd } from 'losat';
import type { EngineGateway, ToolEngine } from 'losat/ports';
import { CUSTOM_OUTPUT } from 'losat/domain';
import { Attention } from './application/attention';
import { Availability } from './application/availability';
import { BlastFormatterDraft } from './application/blast-formatter-draft';
import { Coordinator } from './application/coordinator';
import { SearchDraft } from './application/draft';
import { CountsLoader } from './application/counts';
import { ExampleLoader } from './application/examples';
import { MakeblastdbDraft } from './application/makeblastdb-draft';
import { MaskingDraft } from './application/masking-draft';
import { ResultsBrowser } from './application/results';
import { SessionFiles } from './application/session-files';
import { ToolDraft } from './application/tool-draft';
import type { Downloader } from './ports/download';
import { BrowserAccessionFetcher } from './infra/browser/accession-fetcher';
import { browserCountsChanges } from './infra/browser/counts-changes';
import { browserCountsLocks } from './infra/browser/counts-locks';
import { BrowserCountsDownloader } from './infra/browser/counts-downloader';
import { browserDatabaseMemory } from './infra/browser/db-memory';
import { browserExampleDeletions } from './infra/browser/example-deletions';
import { BrowserFileDownloader } from './infra/browser/example-downloader';
import { browserPage } from './infra/browser/page';
import { browserDownloader, browserFileSaver } from './infra/browser/platform';
import { startDataWorker } from './infra/data-worker/gateway';
import { FakeEngine } from './infra/fake/fake-engine';

export interface App {
  readonly coordinator: Coordinator;
  /** The job being edited on the search screen. */
  readonly draft: SearchDraft;
  /** The WindowMasker job being edited on the WindowMasker tab. */
  readonly toolDraft: ToolDraft;
  /** The dustmasker, segmasker or convert2blastmask job being edited on the Masking tab (WW4c). */
  readonly maskingDraft: MaskingDraft;
  /** The Format Archive (`blast_formatter`) job being edited on the Format Archive tab. */
  readonly blastFormatterDraft: BlastFormatterDraft;
  /** The makeblastdb job being edited on the makeblastdb tab (line WW, WW4b2). */
  readonly makeblastdbDraft: MakeblastdbDraft;
  /** The Example sets that load into the draft (line WF-2). */
  readonly examples: ExampleLoader;
  /** The WindowMasker counts of the taxonomy IDs, downloaded when the user asks (line WW2). */
  readonly counts: CountsLoader;
  /** What the results screen shows (application/results.ts). */
  readonly results: ResultsBrowser;
  /** Wake lock, the warning before leaving, and the check after the page was hidden. */
  readonly attention: Attention;
  /** Whether LOSAT itself works: a part that could not be loaded, or that stopped (line WQA's sweep, D1 and D2). */
  readonly availability: Availability;
  /** True while the engine is the FakeEngine; the UI shows a warning banner. */
  readonly usesFakeEngine: boolean;
}

export interface AppOptions {
  /** Where exports go; the browser's download by default. */
  readonly downloader?: Downloader;
  /** When the engine runtime is renewed (the package's src/infra/engine-worker/policy.ts); tests lower it. */
  readonly renewal?: Partial<RenewalLimits>;
  /** Told how each engine runtime ended (the package's src/infra/engine-worker/gateway.ts); tests record it. */
  readonly onRuntimeEnd?: (end: RuntimeEnd) => void;
}

/**
 * Builds the application. A build with the engine modules (build/reactors.ts) runs the
 * Wasm engine; a build without them is a development build with the FakeEngine.
 */
export function createApp(options: AppOptions = {}): App {
  let dataWorkerEnd: (error: Error) => void = () => undefined;
  const dataWorkerEnded = new Promise<Error>((resolve) => (dataWorkerEnd = resolve));
  const data = startDataWorker((error) => dataWorkerEnd(error));
  const engine: EngineGateway & ToolEngine =
    ENGINE_ASSETS === null
      ? new FakeEngine({ phaseDelayMs: 50 })
      : new WasmEngine({
          assets: ENGINE_ASSETS,
          control: data,
          ...(options.renewal === undefined ? {} : { renewal: options.renewal }),
          ...(options.onRuntimeEnd === undefined ? {} : { onRuntimeEnd: options.onRuntimeEnd }),
        });
  const sessionFiles = new SessionFiles((token) => void data.dropSessionFile(token).catch(() => undefined));
  const countsLocks = browserCountsLocks();
  const toolRuntime =
    ENGINE_ASSETS === null
      ? { path: 'fake' as const, engineBuild: 'fake-engine' }
      : { path: 'serial' as const, engineBuild: engineBuildName('serial', ENGINE_ASSETS.serial) };
  const coordinator = new Coordinator({
    engine,
    data,
    downloader: options.downloader ?? browserDownloader,
    now: () => Date.now(),
    newRunId: () => crypto.randomUUID(),
    digest: sha256Hex,
    dataWorkerEnded,
    toolRuntime,
    sessionFiles,
    ...(countsLocks === undefined ? {} : { countsLocks }),
  });
  const draft = new SearchDraft({
    engine,
    data,
    enqueueAll: (requests) => coordinator.enqueueAll(requests),
    fetcher: new BrowserAccessionFetcher(),
    databaseMemory: browserDatabaseMemory(),
    now: () => performance.now(),
  });
  const counts = new CountsLoader({
    downloader: new BrowserCountsDownloader({ store: data }),
    store: data,
    draft,
    changes: browserCountsChanges(),
  });
  const toolDraft = new ToolDraft({
    engine,
    data,
    enqueueToolAll: (requests) => coordinator.enqueueToolAll(requests),
    isBusy: () => coordinator.isBusy(),
    readSavedCounts: async (taxid) => (await counts.readSavedFile(taxid))?.file,
    sessionFiles,
  });
  const maskingDraft = new MaskingDraft({
    engine,
    data,
    enqueueToolAll: (requests) => coordinator.enqueueToolAll(requests),
    isBusy: () => coordinator.isBusy(),
    sessionFiles,
  });
  const blastFormatterDraft = new BlastFormatterDraft({
    engine,
    enqueueToolAll: (requests) => coordinator.enqueueToolAll(requests),
    getSessionArchiveFile: (runId, fileName) => coordinator.getOutputFile(runId, CUSTOM_OUTPUT, fileName),
    isBusy: () => coordinator.isBusy(),
  });
  const makeblastdbDraft = new MakeblastdbDraft({
    engine,
    enqueueToolAll: (requests) => coordinator.enqueueToolAll(requests),
    isBusy: () => coordinator.isBusy(),
  });
  const tools = [toolDraft, maskingDraft, blastFormatterDraft, makeblastdbDraft];
  const availability = new Availability({ draft, tools, dataWorkerEnded });
  const examples = new ExampleLoader({
    downloader: new BrowserFileDownloader({ inflater: data }),
    store: data,
    draft,
    availability,
    deletions: browserExampleDeletions(),
    saver: browserFileSaver,
  });
  const results = new ResultsBrowser({
    data,
    describe: (program) => engine.describe(program),
    describeTool: (tool) => engine.describeTool(tool),
    runs: coordinator.state,
    verification: VERIFICATION_TABLE,
  });
  const attention = new Attention({
    page: browserPage,
    runs: coordinator.state,
    draft: draft.state,
    tools,
    probe: () => data.storageInfo(),
    checkStored: () => data.checkStoredResults(),
    now: () => Date.now(),
  });
  return {
    coordinator,
    draft,
    toolDraft,
    maskingDraft,
    blastFormatterDraft,
    makeblastdbDraft,
    examples,
    counts,
    results,
    attention,
    availability,
    usesFakeEngine: ENGINE_ASSETS === null,
  };
}
