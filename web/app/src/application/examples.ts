// Example sets (line WF-2, docs/web/fetch_design.md §9): "Example" puts a mode's query in the
// query box and its public subject set in the subject list. The set comes from the copy saved in
// this browser when there is one (no request: it works offline), else it is downloaded, only
// now that the user asked, from the host that the in-place sentence names, and saved for the
// next time. One set loads at a time; another click, Cancel, or a failure ends the load, and a
// failure names the site and leaves the draft as it was. Delete removes the saved copy and the
// subject source made from it.
import { EXAMPLE_BYTE_LIMIT, exampleById, exampleHost, exampleQuery, examplesFor, type ExampleId, type ExampleSet } from '../domain/examples';
import type { DownloadResult, ExampleCopy, ExampleStore, FileDownloader, SavedExample } from '../ports/example';
import type { DraftState, SearchDraft } from './draft';
import { Store } from './store';

export type ExampleLoad =
  /** Looking for the saved copy. */
  | { readonly state: 'reading' }
  | { readonly state: 'downloading'; readonly bytes: number; readonly total?: number }
  | { readonly state: 'saving' }
  | { readonly state: 'failed'; readonly message: string };

export interface ExamplesState {
  /** The saved copy of each set that was looked up; null when there is none. */
  readonly saved: Readonly<Partial<Record<ExampleId, ExampleCopy | null>>>;
  /** The set being loaded, or the one whose load failed. */
  readonly load?: { readonly id: ExampleId; readonly progress: ExampleLoad };
  /** Why the set loaded last was not kept for the next time (it is used in this tab only). */
  readonly note?: { readonly id: ExampleId; readonly text: string };
}

export interface ExamplesDeps {
  readonly downloader: FileDownloader;
  readonly store: ExampleStore;
  readonly draft: Pick<SearchDraft, 'useExample' | 'dropExample'> & { readonly state: Pick<Store<DraftState>, 'get'> };
}

export class ExampleLoader {
  readonly state = new Store<ExamplesState>({ saved: {} });
  private running: AbortController | undefined;
  private readonly pending = new Set<Promise<unknown>>();

  constructor(private readonly deps: ExamplesDeps) {}

  /** Looks up the saved copies of the mode's sets that were not looked up yet (storage only, no request). */
  lookUp(): void {
    for (const set of examplesFor(this.deps.draft.state.get().program)) {
      if (this.state.get().saved[set.id] !== undefined) continue;
      this.track(
        this.deps.store.readExample(set.id, set.url).then(
          (saved) => this.setSaved(set.id, saved?.copy ?? null),
          () => this.setSaved(set.id, null),
        ),
      );
    }
  }

  /** Loads a set into the draft for the current mode. */
  load(id: ExampleId): void {
    this.cancel();
    const controller = new AbortController();
    this.running = controller;
    this.track(this.loadSet(exampleById(id), controller.signal));
  }

  /** Stops the running load; the draft stays as it was. */
  cancel(): void {
    this.running?.abort();
    this.running = undefined;
    const { load } = this.state.get();
    if (load !== undefined) this.update({ load: undefined });
  }

  /** Removes the saved copy of a set, and the subject source made from it. */
  async delete(id: ExampleId): Promise<void> {
    if (this.state.get().load?.id === id) this.cancel();
    this.deps.draft.dropExample(id);
    try {
      await this.deps.store.deleteExample(id);
    } catch (error) {
      this.update({ note: { id, text: `The saved copy could not be removed: ${messageOf(error)}` } });
      return;
    }
    this.setSaved(id, null);
    if (this.state.get().note?.id === id) this.update({ note: undefined });
  }

  /** Resolves when no load or look-up is running (tests). */
  async idle(): Promise<void> {
    while (this.pending.size > 0) await Promise.allSettled([...this.pending]);
  }

  private async loadSet(set: ExampleSet, signal: AbortSignal): Promise<void> {
    const program = this.deps.draft.state.get().program;
    const query = exampleQuery(set, program);
    const progress = (next: ExampleLoad) => {
      if (!signal.aborted) this.update({ load: { id: set.id, progress: next } });
    };
    const fail = (message: string) => progress({ state: 'failed', message });
    this.update({ note: undefined });
    progress({ state: 'reading' });
    const kept = await this.deps.store.readExample(set.id, set.url).catch(() => undefined);
    if (signal.aborted) return;
    if (kept !== undefined) {
      this.finish(set, query, kept.file, kept.copy);
      return;
    }
    progress({ state: 'downloading', bytes: 0 });
    let result: DownloadResult;
    try {
      result = await this.deps.downloader.download(
        { url: set.url, maxBytes: EXAMPLE_BYTE_LIMIT },
        ({ bytes, total }) => progress({ state: 'downloading', bytes, ...(total === undefined ? {} : { total }) }),
        signal,
      );
    } catch (error) {
      if (!signal.aborted) fail(`The example could not be downloaded from ${exampleHost(set)}: ${messageOf(error)}`);
      return;
    }
    if (signal.aborted) return;
    if (!result.ok) {
      fail(downloadProblem(result, set));
      return;
    }
    progress({ state: 'saving' });
    const { data, downloadedBytes, version } = result;
    let saved: SavedExample | undefined;
    let note: string | undefined;
    try {
      saved = await this.deps.store.saveExample(
        { id: set.id, url: set.url, downloadedBytes, ...(version === undefined ? {} : { version }) },
        data,
      );
    } catch (error) {
      note = `Not kept for the next time: ${messageOf(error)}`;
    }
    if (signal.aborted) return;
    if (saved === undefined) {
      this.finish(set, query, data, null);
      if (note !== undefined) this.update({ note: { id: set.id, text: note } });
    } else {
      this.finish(set, query, saved.file, saved.copy);
    }
  }

  private finish(set: ExampleSet, query: string, data: Blob, copy: ExampleCopy | null): void {
    this.running = undefined;
    this.deps.draft.useExample(set.id, query, new File([data], set.fileName, { type: 'text/plain' }));
    this.setSaved(set.id, copy);
    this.update({ load: undefined });
  }

  private setSaved(id: ExampleId, copy: ExampleCopy | null): void {
    this.update({ saved: { ...this.state.get().saved, [id]: copy } });
  }

  private track(work: Promise<unknown>): void {
    this.pending.add(work);
    void work.finally(() => this.pending.delete(work));
  }

  /** Changes the state; a key given as undefined is removed. */
  private update(change: { readonly [K in keyof ExamplesState]?: ExamplesState[K] | undefined }): void {
    const next = { ...this.state.get(), ...change };
    for (const key of Object.keys(change) as Array<keyof ExamplesState>) {
      if (change[key] === undefined) delete next[key];
    }
    this.state.set(next as ExamplesState);
  }
}

/** Why a download brought no set, for the Example line. */
function downloadProblem(result: Exclude<DownloadResult, { ok: true }>, set: ExampleSet): string {
  const host = exampleHost(set);
  switch (result.reason) {
    case 'unreachable':
      return `${host} could not be reached (${result.detail}).`;
    case 'refused':
      return result.status >= 500 ? `${host} answered HTTP ${result.status}: try again in a moment.` : `${host} answered HTTP ${result.status}.`;
    case 'corrupt':
      return `The file from ${host} could not be decompressed (${result.detail}).`;
    case 'too-large':
      return `The file from ${host} is larger than ${EXAMPLE_BYTE_LIMIT / (1024 * 1024)} MiB decompressed.`;
  }
}

function messageOf(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
