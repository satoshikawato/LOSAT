// The WindowMasker counts of the taxonomy IDs (line WW2, docs/web/fetch_design.md §10): the line
// under the taxonomy ID box. "Download" asks ftp.ncbi.nlm.nih.gov for the file of the ID, only now
// that the user asked, and the Data worker keeps it as it arrives; "Continue" asks for the rest of a
// download that stopped; "Delete" removes it. A saved file is handed to the draft, which the search
// reads it from (SearchDraft.setCounts); a look-up of the saved copies sends no request. One file
// downloads at a time: another click, Cancel, or a failure ends it, and a failure names the site and
// leaves what was saved. A tab that keeps nothing across sessions holds the file it downloaded
// in this tab only.
import { COUNTS_HOST, countsUrl, type CountsOrigin } from '../domain/counts';
import { refusedSentence, unreachableSentence } from '../domain/download-sentences';
import type { CountsCopy, CountsDeletions, CountsDownloader, CountsResult, CountsStore } from '../ports/counts';
import type { SearchDraft } from './draft';
import { Store } from './store';

export type CountsLoad =
  /** Looking for the saved copy. */
  | { readonly state: 'reading' }
  | { readonly state: 'downloading'; readonly bytes: number; readonly total?: number }
  | { readonly state: 'failed'; readonly message: string };

export interface CountsState {
  /** The saved copy of each taxonomy ID that was looked up (partial or complete); null when there is none. */
  readonly saved: Readonly<Record<string, CountsCopy | null>>;
  /** The download in progress, or the one that failed. */
  readonly load?: { readonly taxid: number; readonly progress: CountsLoad };
  /** The taxonomy IDs whose file this tab holds and the browser does not keep (no OPFS or no Web Locks). */
  readonly inTab: Readonly<Record<string, true>>;
  /** A word on a taxonomy ID's copy: why it could not be removed. */
  readonly note?: { readonly taxid: number; readonly text: string };
}

export interface CountsDeps {
  readonly downloader: CountsDownloader;
  readonly store: CountsStore;
  readonly draft: Pick<SearchDraft, 'setCounts'>;
  /** Deletions between tabs; without it, a tab hears only its own. */
  readonly deletions?: CountsDeletions;
}

/** The origin that a run records for a copy. */
export function originOf(copy: CountsCopy): CountsOrigin {
  return { taxid: copy.taxid, url: copy.url, lastModified: copy.lastModified };
}

export class CountsLoader {
  readonly state = new Store<CountsState>({ saved: {}, inTab: {} });
  private running: AbortController | undefined;
  private readonly pending = new Set<Promise<unknown>>();
  /** Deletions in progress; a download of the same taxonomy ID waits for its deletion. */
  private readonly deleting = new Map<number, Promise<void>>();
  /** The files that this tab downloaded and the browser does not keep. */
  private readonly local = new Map<number, { readonly copy: CountsCopy; readonly file: File }>();

  constructor(private readonly deps: CountsDeps) {
    deps.deletions?.subscribe((taxid) => this.forget(taxid));
  }

  /** Looks for the saved copy of a taxonomy ID, unless it was looked up (storage only, no request). */
  refresh(taxid: number): void {
    if (this.state.get().saved[taxid] !== undefined) return;
    this.track(this.reread(taxid));
  }

  /** Downloads the counts of a taxonomy ID, or continues the download that stopped. */
  download(taxid: number): void {
    this.cancel();
    const controller = new AbortController();
    this.running = controller;
    this.track(this.load(taxid, controller));
  }

  /** Stops the running download; what was saved stays. */
  cancel(): void {
    this.running?.abort();
    this.running = undefined;
    if (this.state.get().load !== undefined) this.update({ load: undefined });
  }

  /** Removes the saved copy of a taxonomy ID, here and in the other tabs. */
  async delete(taxid: number): Promise<void> {
    if (this.state.get().load?.taxid === taxid) this.cancel();
    const deletion = this.deleteCopy(taxid);
    this.deleting.set(taxid, deletion);
    try {
      await deletion;
    } finally {
      if (this.deleting.get(taxid) === deletion) this.deleting.delete(taxid);
    }
  }

  /** Resolves when no download or look-up is running (tests). */
  async idle(): Promise<void> {
    while (this.pending.size > 0) await Promise.allSettled([...this.pending]);
  }

  private async deleteCopy(taxid: number): Promise<void> {
    await this.deleting.get(taxid);
    // The file goes from the draft first, so that no search is queued with a file about to go.
    this.deps.draft.setCounts(taxid, undefined);
    this.holdInTab(taxid, false);
    try {
      await this.deps.store.deleteCounts(taxid);
    } catch (error) {
      this.update({ note: { taxid, text: `The saved counts could not be removed: ${messageOf(error)}` } });
      await this.reread(taxid);
      return;
    }
    this.deps.deletions?.announce(taxid);
    this.setSaved(taxid, null);
    if (this.state.get().note?.taxid === taxid) this.update({ note: undefined });
  }

  /** Another tab deleted the file: this tab lets go of it. */
  private forget(taxid: number): void {
    this.deps.draft.setCounts(taxid, undefined);
    this.holdInTab(taxid, false);
    this.setSaved(taxid, null);
  }

  /** Reads the copy from the storage and tells the draft when it is whole. */
  private async reread(taxid: number): Promise<void> {
    const url = countsUrl(taxid);
    await this.deleting.get(taxid);
    const kept = this.local.get(taxid);
    if (kept !== undefined) {
      this.setSaved(taxid, kept.copy);
      return;
    }
    try {
      const saved = await this.deps.store.readCounts(taxid, url);
      if (saved?.file !== undefined) this.deps.draft.setCounts(taxid, { file: saved.file, origin: originOf(saved.copy) });
      this.setSaved(taxid, saved?.copy ?? null);
    } catch {
      this.setSaved(taxid, null);
    }
  }

  private async load(taxid: number, controller: AbortController): Promise<void> {
    const { signal } = controller;
    const url = countsUrl(taxid);
    const progress = (next: CountsLoad) => {
      if (!signal.aborted) this.update({ load: { taxid, progress: next } });
    };
    const fail = (message: string) => progress({ state: 'failed', message });
    this.update({ note: undefined });
    progress({ state: 'reading' });
    // A copy being deleted is gone first.
    await this.deleting.get(taxid);
    const saved = this.local.has(taxid) ? undefined : await this.deps.store.readCounts(taxid, url).catch(() => undefined);
    if (signal.aborted) return;
    const resume = saved !== undefined && !saved.copy.complete && saved.copy.savedBytes > 0 ? saved.copy : undefined;
    if (saved?.file !== undefined || this.local.has(taxid)) {
      this.finish(taxid, saved?.copy ?? this.local.get(taxid)!.copy, saved?.file ?? this.local.get(taxid)!.file, controller);
      return;
    }
    progress({ state: 'downloading', bytes: resume?.savedBytes ?? 0, ...(resume === undefined ? {} : { total: resume.size }) });
    let result: CountsResult;
    try {
      result = await this.deps.downloader.download(
        { taxid, url, ...(resume === undefined ? {} : { resume }) },
        ({ bytes, total }) => progress({ state: 'downloading', bytes, ...(total === undefined ? {} : { total }) }),
        signal,
      );
    } catch (error) {
      // A stop (Cancel, or another download) left what came; a failure names the site.
      if (!signal.aborted) fail(`The counts could not be downloaded from ${COUNTS_HOST}: ${messageOf(error)}`);
      await this.reread(taxid);
      return;
    }
    if (!result.ok) {
      if (!signal.aborted) fail(problem(result, taxid));
      await this.reread(taxid);
      return;
    }
    if (!result.kept) {
      this.local.set(taxid, { copy: result.copy, file: result.file });
      this.holdInTab(taxid, true);
    }
    this.finish(taxid, result.copy, result.file, controller);
  }

  /** The file is whole: the draft reads it from here on. */
  private finish(taxid: number, copy: CountsCopy, file: File, controller: AbortController): void {
    if (this.running === controller) this.running = undefined;
    this.deps.draft.setCounts(taxid, { file, origin: originOf(copy) });
    this.setSaved(taxid, copy);
    if (this.state.get().load?.taxid === taxid) this.update({ load: undefined });
  }

  /** Tells whether this tab holds the file of a taxonomy ID that the browser does not keep. */
  private holdInTab(taxid: number, held: boolean): void {
    if (!held) this.local.delete(taxid);
    const inTab = { ...this.state.get().inTab };
    if (held) inTab[taxid] = true;
    else delete inTab[taxid];
    this.update({ inTab });
  }

  private setSaved(taxid: number, copy: CountsCopy | null): void {
    this.update({ saved: { ...this.state.get().saved, [taxid]: copy } });
  }

  private track(work: Promise<unknown>): void {
    this.pending.add(work);
    void work.finally(() => this.pending.delete(work));
  }

  /** Changes the state; a key given as undefined is removed. */
  private update(change: { readonly [K in keyof CountsState]?: CountsState[K] | undefined }): void {
    const next = { ...this.state.get(), ...change };
    for (const key of Object.keys(change) as Array<keyof CountsState>) {
      if (change[key] === undefined) delete next[key];
    }
    this.state.set(next as CountsState);
  }
}

/** Why a download brought no counts, for the line. */
function problem(result: Exclude<CountsResult, { ok: true }>, taxid: number): string {
  switch (result.reason) {
    case 'unreachable':
      return unreachableSentence(COUNTS_HOST, result.detail);
    case 'refused':
      return result.status === 404 ? `${COUNTS_HOST} has no WindowMasker counts for ${taxid}.` : refusedSentence(COUNTS_HOST, result.status);
    case 'busy':
      return `Another tab is downloading the counts for ${taxid}.`;
    case 'not-kept':
      return `Not kept: ${result.detail}`;
  }
}

function messageOf(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
