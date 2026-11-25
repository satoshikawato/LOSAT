// The WindowMasker counts of the taxonomy IDs (line WW2, docs/web/fetch_design.md §10): the line
// under the taxonomy ID box. "Download" asks ftp.ncbi.nlm.nih.gov for the file of the ID, only now
// that the user asked, and the Data worker keeps it as it arrives; "Continue" asks for the rest of a
// download that stopped; "Delete" removes it. A saved file is handed to the draft, which the search
// reads it from (SearchDraft.setCounts); a look-up of the saved copies sends no request. One file
// downloads at a time: another click, Cancel, or a failure ends it, and a failure names the site and
// leaves what was saved. A tab that keeps nothing across sessions holds the file it downloaded
// in this tab only.
import { COUNTS_HOST, countsEntryOf, countsUrl, type CountsMadeMeta, type CountsOrigin } from '../domain/counts';
import { refusedSentence, unreachableSentence } from '../domain/download-sentences';
import {
  COUNTS_BUSY,
  COUNTS_EXISTS,
  COUNTS_IN_USE,
  type CountsChanges,
  type CountsCopy,
  type CountsDownloader,
  type CountsResult,
  type CountsStore,
  type SavedCounts,
  type SavedCountsEntry,
} from '../ports/counts';
import type { SearchDraft } from './draft';
import { Store } from './store';

/** Slice size when streaming a made counts file from its run directory to OPFS (never hold a large file in memory). */
const SAVE_PART_BYTES = 8 * 1024 * 1024;

export type CountsLoad =
  /** Looking for the saved copy. */
  | { readonly state: 'reading' }
  | { readonly state: 'downloading'; readonly bytes: number; readonly total?: number }
  | { readonly state: 'failed'; readonly message: string };

export type SaveCountsResult =
  | { readonly ok: true; readonly copy: CountsCopy; readonly kept: boolean }
  | { readonly ok: false; readonly reason: 'exists'; readonly existing: CountsCopy }
  | { readonly ok: false; readonly reason: 'busy'; readonly message: string }
  | { readonly ok: false; readonly reason: 'in-use'; readonly message: string }
  | { readonly ok: false; readonly reason: 'failed'; readonly message: string };

export interface CountsState {
  /** The saved copy of each taxonomy ID that was looked up (partial or complete); null when there is none. */
  readonly saved: Readonly<Record<string, CountsCopy | null>>;
  /** The download in progress, or the one that failed. */
  readonly load?: { readonly taxid: number; readonly progress: CountsLoad };
  /** Browser-made counts saves currently being written by this tab, keyed by taxonomy ID. */
  readonly saving: Readonly<Record<string, { readonly run: number }>>;
  /** The taxonomy IDs whose file this tab holds and the browser does not keep (no OPFS or no Web Locks). */
  readonly inTab: Readonly<Record<string, true>>;
  /** A word on a taxonomy ID's copy: why it could not be removed. */
  readonly note?: { readonly taxid: number; readonly text: string };
  /** Why the last Delete of a taxonomy ID's copy was refused (a run reads it, or another writer holds it); it stays. */
  readonly refusal?: { readonly taxid: number; readonly text: string };
  /** Complete saved copies available in this browser, ordered by taxonomy ID. */
  readonly list: readonly SavedCountsEntry[];
}

export interface CountsDeps {
  readonly downloader: CountsDownloader;
  readonly store: CountsStore;
  readonly draft: Pick<SearchDraft, 'setCounts'>;
  /** Changes between tabs (save, replace, delete); without it, a tab hears only its own. */
  readonly changes?: CountsChanges;
}

/** The origin that a run records for a copy. */
export function originOf(copy: CountsCopy): CountsOrigin {
  if (copy.made !== undefined) return { kind: 'made', taxid: copy.taxid, made: copy.made };
  return { kind: 'downloaded', taxid: copy.taxid, url: copy.url, lastModified: copy.lastModified };
}

/** Identity of a browser-made counts copy saved from a tool run output tab. */
export interface SavedTabCopy {
  readonly taxid: number;
  readonly run: number;
  readonly madeAt: number;
  readonly size: number;
}

/** Captures the identity of a saved browser-made counts copy for `ToolRunView`. */
export function savedTabCopyOf(copy: CountsCopy): SavedTabCopy | undefined {
  if (copy.made === undefined) return undefined;
  return { taxid: copy.taxid, run: copy.made.run, madeAt: copy.made.madeAt, size: copy.size };
}

/**
 * Returns `saved.taxid` only while `stored` is still that exact saved copy (same taxonomy ID,
 * `made.run`, `madeAt`, and `size`). Returns `undefined` when another run or an NCBI download
 * replaced the copy or it was deleted, so "Use in BLASTN" uses the run's own file instead.
 */
export function activeSavedCountsTaxid(
  saved: SavedTabCopy | undefined,
  stored: CountsCopy | null | undefined,
): number | undefined {
  if (saved === undefined || stored == null || !stored.complete || stored.made === undefined) return undefined;
  if (
    stored.taxid !== saved.taxid ||
    stored.made.run !== saved.run ||
    stored.made.madeAt !== saved.madeAt ||
    stored.size !== saved.size
  ) {
    return undefined;
  }
  return saved.taxid;
}

/** A downloaded copy counts only if its URL matches this release's NCBI URL; made copies always count. */
function accepted(saved: SavedCounts | undefined, taxid: number): SavedCounts | undefined {
  if (saved === undefined) return undefined;
  if (saved.origin === 'made') return saved;
  return saved.copy.url === countsUrl(taxid) ? saved : undefined;
}

type ActiveWrite =
  | { readonly kind: 'save'; readonly op: object }
  | { readonly kind: 'download'; readonly op: AbortController; readonly done: Promise<void> };

export class CountsLoader {
  readonly state = new Store<CountsState>({ saved: {}, saving: {}, inTab: {}, list: [] });
  private running: AbortController | undefined;
  private readonly pending = new Set<Promise<unknown>>();
  /** Deletions in progress; a download or save of the same taxonomy ID waits for its deletion. */
  private readonly deleting = new Map<number, Promise<void>>();
  /** Active writes in this tab, keyed by taxonomy ID: 'save' from writeMade begin to end, 'download' from load start to end. */
  private readonly writing = new Map<number, ActiveWrite>();
  /** The files that this tab holds and the browser does not keep. */
  private readonly local = new Map<number, { readonly copy: CountsCopy; readonly file: File }>();
  /** The number of the latest read of the list (reloadList): an earlier read that ends later is dropped. */
  private listRead = 0;
  /** The number of the latest read or state update per taxonomy ID: an earlier read that ends later is dropped. */
  private readonly taxidRead = new Map<number, number>();

  constructor(private readonly deps: CountsDeps) {
    deps.changes?.subscribe((taxid) => this.onRemoteChange(taxid));
  }

  /** Looks for the saved copy of a taxonomy ID, unless it was looked up (storage only, no request). */
  refresh(taxid: number): void {
    if (this.state.get().saved[taxid] !== undefined) return;
    this.track(this.reread(taxid));
  }

  /** Refreshes the list of complete copies saved in this browser (storage only, no request). */
  refreshList(): void {
    this.track(this.reloadList());
  }

  /** Reads the kept complete file for a taxonomy ID (storage or this tab only, no request). */
  async readSavedFile(taxid: number): Promise<SavedCounts | undefined> {
    await this.deleting.get(taxid);
    const kept = this.local.get(taxid);
    if (kept !== undefined) {
      return { copy: kept.copy, origin: kept.copy.made === undefined ? 'downloaded' : 'made', file: kept.file };
    }
    const saved = accepted(await this.deps.store.readCounts(taxid).catch(() => undefined), taxid);
    return saved?.file === undefined ? undefined : saved;
  }

  /**
   * Saves a WindowMasker counts file made in this browser for `taxid`, streaming it in 8 MiB parts.
   * When `replace` is not true and a copy already exists for `taxid`, returns `{ ok: false, reason: 'exists' }`
   * without touching the stored copy.
   */
  saveMade(
    taxid: number,
    file: File,
    made: CountsMadeMeta,
    options: { readonly replace?: boolean } = {},
  ): Promise<SaveCountsResult> {
    const work = this.writeMade(taxid, file, made, options.replace === true);
    this.track(work);
    return work;
  }

  /** Downloads the counts of a taxonomy ID, or continues the download that stopped. */
  download(taxid: number): void {
    const active = this.writing.get(taxid);
    // This tab's save, or its download that goes on, is refused in place; a download that was
    // stopped (Cancel, another ID) is waited for below.
    const stopped = active?.kind === 'download' && active.op.signal.aborted;
    if (active !== undefined && !stopped) {
      this.update({
        load: {
          taxid,
          progress: {
            state: 'failed',
            message: busySentence(taxid, 'download', active.kind === 'save' ? 'this-save' : 'this-download'),
          },
        },
      });
      return;
    }
    this.cancel();
    const controller = new AbortController();
    this.running = controller;
    // A download of the ID that was stopped (Cancel, another ID) still lets go of it: start after it.
    const before = active?.kind === 'download' ? active.done : undefined;
    if (before !== undefined) this.update({ load: { taxid, progress: { state: 'reading' } } });
    const done: Promise<void> = (async () => {
      if (before !== undefined) await before.catch(() => undefined);
      if (controller.signal.aborted) {
        this.releaseWriting(taxid, controller);
        return;
      }
      await this.load(taxid, controller);
    })();
    this.writing.set(taxid, { kind: 'download', op: controller, done });
    this.track(done);
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

  /** Resolves when no download, save or look-up is running (tests). */
  async idle(): Promise<void> {
    while (this.pending.size > 0) await Promise.allSettled([...this.pending]);
  }

  private async writeMade(taxid: number, file: File, made: CountsMadeMeta, replace: boolean): Promise<SaveCountsResult> {
    // Plain structured-cloneable copy so reactive UI wrappers never reach Worker.postMessage.
    const plainMade: CountsMadeMeta = {
      run: made.run,
      argv: [...made.argv],
      inputs: made.inputs.map((entry) => ({ name: entry.name, size: entry.size })),
      engineBuild: made.engineBuild,
      madeAt: made.madeAt,
    };
    await this.deleting.get(taxid);
    const active = this.writing.get(taxid);
    if (active !== undefined) {
      return {
        ok: false,
        reason: 'busy',
        message: busySentence(taxid, 'save', active.kind === 'save' ? 'this-save' : 'this-download'),
      };
    }
    const local = this.local.get(taxid);
    if (local !== undefined && !replace) {
      return { ok: false, reason: 'exists', existing: local.copy };
    }
    const begin = (allowReplace: boolean) =>
      this.deps.store.beginCounts(taxid, { made: plainMade, size: file.size, from: 0, replace: allowReplace });

    let token: string | undefined;
    // Whether a copy of before is overwritten: a failure after the start then says that it is gone.
    let replacing = replace;
    try {
      token = await begin(replace);
    } catch (error) {
      if (errorName(error) === COUNTS_EXISTS) {
        const existing = accepted(await this.deps.store.readCounts(taxid).catch(() => undefined), taxid);
        if (existing !== undefined) {
          return { ok: false, reason: 'exists', existing: existing.copy };
        }
        // A download from an older/other URL does not count as a saved copy for this release: overwrite it.
        try {
          token = await begin(true);
          replacing = true;
        } catch (retryError) {
          return refusedSave(retryError, taxid, this.writerOf(taxid));
        }
      } else {
        return refusedSave(error, taxid, this.writerOf(taxid));
      }
    }

    if (token === undefined) {
      const copy: CountsCopy = {
        taxid,
        url: '',
        lastModified: '',
        size: file.size,
        savedBytes: file.size,
        complete: true,
        savedAt: plainMade.madeAt,
        made: plainMade,
      };
      this.local.set(taxid, { copy, file });
      this.holdInTab(taxid, true);
      this.deps.draft.setCounts(taxid, { file, origin: originOf(copy) });
      this.setSaved(taxid, copy);
      this.clearMessages(taxid);
      await this.reloadList();
      return { ok: true, copy, kept: false };
    }

    const op = {};
    this.writing.set(taxid, { kind: 'save', op });
    this.setSaving(taxid, plainMade.run);
    this.clearMessages(taxid);
    const clearSave = () => {
      // Only this save's own entry: a later save of the same ID keeps its state.
      if (this.writing.get(taxid)?.op !== op) return;
      this.writing.delete(taxid);
      this.setSaving(taxid, undefined);
    };

    try {
      if (replacing) {
        // The copy of before is gone (meta incomplete, data truncated): drop it here and tell the other tabs.
        this.deps.draft.setCounts(taxid, undefined);
        this.setSaved(taxid, null);
        this.deps.changes?.announce(taxid);
        await this.reloadList();
      }
      for (let offset = 0; offset < file.size; offset += SAVE_PART_BYTES) {
        const slice = new Uint8Array(await file.slice(offset, Math.min(file.size, offset + SAVE_PART_BYTES)).arrayBuffer());
        await this.deps.store.appendCounts(token, slice);
      }
      const saved = await this.deps.store.finishCounts(token);
      if (saved.file === undefined) throw new Error('the saved file cannot be read');
      this.holdInTab(taxid, false);
      this.deps.draft.setCounts(taxid, { file: saved.file, origin: originOf(saved.copy) });
      this.nextTaxidRead(taxid);
      const saving = { ...this.state.get().saving };
      delete saving[taxid];
      this.releaseWriting(taxid, op);
      this.update({
        saved: { ...this.state.get().saved, [taxid]: saved.copy },
        saving,
        ...(this.state.get().note?.taxid === taxid ? { note: undefined } : {}),
        ...(this.state.get().load?.taxid === taxid ? { load: undefined } : {}),
      });
      this.deps.changes?.announce(taxid);
      await this.reloadList();
      return { ok: true, copy: saved.copy, kept: true };
    } catch (error) {
      await this.deps.store.closeCounts(token).catch(() => undefined);
      await this.reread(taxid);
      clearSave();
      this.clearMessages(taxid, { note: false, refusal: false });
      this.deps.changes?.announce(taxid);
      await this.reloadList();
      const gone = replacing ? ` The counts saved before for ${taxid} are no longer kept.` : '';
      return { ok: false, reason: 'failed', message: `${messageOf(error)}${gone}` };
    } finally {
      clearSave();
    }
  }

  private async deleteCopy(taxid: number): Promise<void> {
    await this.deleting.get(taxid);
    // The file goes from the draft first, so that no search is queued with a file about to go.
    this.nextTaxidRead(taxid);
    this.deps.draft.setCounts(taxid, undefined);
    this.holdInTab(taxid, false);
    if (this.state.get().refusal?.taxid === taxid) this.update({ refusal: undefined });
    try {
      await this.deps.store.deleteCounts(taxid);
    } catch (error) {
      // Refused in place (a queued or running run reads the copy, or another writer holds it): the copy stays.
      const refusal = inUseOrBusy(error, taxid, this.writerOf(taxid));
      if (refusal !== undefined) this.update({ refusal: { taxid, text: refusal } });
      else this.update({ note: { taxid, text: `The saved counts could not be removed: ${messageOf(error)}` } });
      // Not reread: it waits for the deletions of the ID, and this is one of them.
      await this.readBack(taxid);
      await this.reloadList();
      return;
    }
    this.deps.changes?.announce(taxid);
    this.setSaved(taxid, null);
    if (this.state.get().note?.taxid === taxid) this.update({ note: undefined });
    await this.reloadList();
  }

  /** Another tab saved, replaced or deleted the counts for `taxid`: drop any stale handle and read again. */
  private onRemoteChange(taxid: number): void {
    this.deps.draft.setCounts(taxid, undefined);
    this.holdInTab(taxid, false);
    this.track(this.reread(taxid));
    this.track(this.reloadList());
  }

  /** Reads the copy from the storage, after a deletion of the ID, and tells the draft whether a whole copy is present. */
  private async reread(taxid: number): Promise<void> {
    const read = this.nextTaxidRead(taxid);
    await this.deleting.get(taxid);
    if (read !== this.taxidRead.get(taxid)) return;
    await this.readBack(taxid, read);
  }

  /** Reads the copy from the storage and tells the draft whether a whole copy is present. */
  private async readBack(taxid: number, read = this.nextTaxidRead(taxid)): Promise<void> {
    const kept = this.local.get(taxid);
    if (kept !== undefined) {
      this.deps.draft.setCounts(taxid, { file: kept.file, origin: originOf(kept.copy) });
      this.setSaved(taxid, kept.copy);
      return;
    }
    try {
      const saved = accepted(await this.deps.store.readCounts(taxid), taxid);
      // A later read or update started for this ID while this one ran: its state is the newer one.
      if (read !== this.taxidRead.get(taxid)) return;
      if (saved?.file !== undefined) {
        this.deps.draft.setCounts(taxid, { file: saved.file, origin: originOf(saved.copy) });
      } else {
        this.deps.draft.setCounts(taxid, undefined);
      }
      this.setSaved(taxid, saved?.copy ?? null);
    } catch {
      if (read !== this.taxidRead.get(taxid)) return;
      this.deps.draft.setCounts(taxid, undefined);
      this.setSaved(taxid, null);
    }
  }

  private async reloadList(): Promise<void> {
    const read = ++this.listRead;
    const persisted = await this.deps.store.listCounts().catch(() => []);
    // A later read started while this one ran: its list is the newer one.
    if (read !== this.listRead) return;
    const byTaxid = new Map<number, SavedCountsEntry>();
    for (const entry of persisted) {
      if (entry.origin === 'made' || entry.copy.url === countsUrl(entry.taxid)) {
        byTaxid.set(entry.taxid, entry);
      }
    }
    for (const [taxid, { copy }] of this.local) {
      if (copy.complete && !byTaxid.has(taxid)) byTaxid.set(taxid, countsEntryOf(copy));
    }
    const list = [...byTaxid.values()].sort((a, b) => a.taxid - b.taxid);
    this.update({ list });
  }

  private async load(taxid: number, controller: AbortController): Promise<void> {
    const releaseWriting = () => this.releaseWriting(taxid, controller);
    try {
      const { signal } = controller;
      const url = countsUrl(taxid);
      const progress = (next: CountsLoad) => {
        if (!signal.aborted) this.update({ load: { taxid, progress: next } });
      };
      const fail = (message: string) => progress({ state: 'failed', message });
      this.update({ note: undefined, ...(this.state.get().refusal?.taxid === taxid ? { refusal: undefined } : {}) });
      progress({ state: 'reading' });
      // A copy being deleted is gone first.
      await this.deleting.get(taxid);
      const saved = this.local.has(taxid)
        ? undefined
        : accepted(await this.deps.store.readCounts(taxid).catch(() => undefined), taxid);
      if (signal.aborted) return;
      const resume =
        saved !== undefined && saved.origin === 'downloaded' && !saved.copy.complete && saved.copy.savedBytes > 0
          ? saved.copy
          : undefined;
      if (saved?.file !== undefined || this.local.has(taxid)) {
        releaseWriting();
        this.finish(taxid, saved?.copy ?? this.local.get(taxid)!.copy, saved?.file ?? this.local.get(taxid)!.file, controller, false);
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
        releaseWriting();
        // A stop (Cancel, or another download) left what came; a failure names the site.
        if (!signal.aborted) fail(`The counts could not be downloaded from ${COUNTS_HOST}: ${messageOf(error)}`);
        await this.reread(taxid);
        return;
      }
      releaseWriting();
      if (!result.ok) {
        if (!signal.aborted) fail(problem(result, taxid));
        await this.reread(taxid);
        return;
      }
      if (!result.kept) {
        this.local.set(taxid, { copy: result.copy, file: result.file });
        this.holdInTab(taxid, true);
      }
      this.finish(taxid, result.copy, result.file, controller, result.kept);
    } finally {
      releaseWriting();
    }
  }

  /** The file is whole: the draft reads it from here on. */
  private finish(taxid: number, copy: CountsCopy, file: File, controller: AbortController, announce: boolean): void {
    const current = this.running === controller;
    if (current) this.running = undefined;
    this.deps.draft.setCounts(taxid, { file, origin: originOf(copy) });
    this.setSaved(taxid, copy);
    if (current && this.state.get().load?.taxid === taxid) this.update({ load: undefined });
    if (announce) this.deps.changes?.announce(taxid);
    this.track(this.reloadList());
  }

  /** Tells whether this tab holds the file of a taxonomy ID that the browser does not keep. */
  private holdInTab(taxid: number, held: boolean): void {
    if (!held) this.local.delete(taxid);
    const inTab = { ...this.state.get().inTab };
    if (held) inTab[taxid] = true;
    else delete inTab[taxid];
    this.update({ inTab });
  }

  private nextTaxidRead(taxid: number): number {
    const read = (this.taxidRead.get(taxid) ?? 0) + 1;
    this.taxidRead.set(taxid, read);
    return read;
  }

  private releaseWriting(taxid: number, op: object): void {
    if (this.writing.get(taxid)?.op === op) this.writing.delete(taxid);
  }

  private writerOf(taxid: number): BusyWriter {
    const active = this.writing.get(taxid)?.kind;
    if (active === 'save') return 'this-save';
    if (active === 'download') return 'this-download';
    return 'other-tab';
  }

  private setSaving(taxid: number, run: number | undefined): void {
    const saving = { ...this.state.get().saving };
    if (run === undefined) delete saving[taxid];
    else saving[taxid] = { run };
    this.update({ saving });
  }

  private clearMessages(taxid: number, which: { readonly note?: boolean; readonly refusal?: boolean } = {}): void {
    const state = this.state.get();
    const change: { load?: undefined; note?: undefined; refusal?: undefined } = {};
    if (state.load?.taxid === taxid) change.load = undefined;
    if (which.note !== false && state.note?.taxid === taxid) change.note = undefined;
    if (which.refusal !== false && state.refusal?.taxid === taxid) change.refusal = undefined;
    if (Object.keys(change).length > 0) this.update(change);
  }

  private setSaved(taxid: number, copy: CountsCopy | null): void {
    this.nextTaxidRead(taxid);
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
      return busySentence(taxid, 'download', 'other-tab');
    case 'in-use':
      return inUseSentence(taxid, 'download');
    case 'changed':
      return `Another tab changed the saved counts for ${taxid}: press Continue again.`;
    case 'not-kept':
      return `Not kept: ${result.detail}`;
  }
}

type RefusalAction = 'delete' | 'replace' | 'save' | 'download';
type BusyAction = 'download' | 'save' | 'delete';
type BusyWriter = 'this-save' | 'this-download' | 'other-tab';

function refusalSentence(action: RefusalAction, what: string): string {
  switch (action) {
    case 'delete':
      return `Not removed: ${what} Press Delete again when it has ended.`;
    case 'replace':
      return `Not replaced: ${what} Save again when it has ended.`;
    case 'save':
      return `Not saved: ${what} Save again when it has ended.`;
    case 'download':
      return `Not downloaded: ${what} Download again when it has ended.`;
  }
}

function busySentence(taxid: number, action: BusyAction, writer: BusyWriter): string {
  const what =
    writer === 'this-save'
      ? `this tab is saving the counts for ${taxid}.`
      : writer === 'this-download'
        ? `this tab is downloading the counts for ${taxid}.`
        : `another tab is saving or downloading the counts for ${taxid}.`;
  return refusalSentence(action, what);
}

type CountsAction = 'delete' | 'replace' | 'download';

function inUseSentence(taxid: number, action: CountsAction): string {
  return refusalSentence(action, `a run that reads the counts for ${taxid} was queued or running.`);
}

function errorName(error: unknown): unknown {
  return (error as { name?: unknown } | null)?.name;
}

/** The sentence of a refusal in place (the ID's lock is held), or undefined for another error. */
function inUseOrBusy(error: unknown, taxid: number, writer: BusyWriter): string | undefined {
  const name = errorName(error);
  if (name === COUNTS_IN_USE) return inUseSentence(taxid, 'delete');
  if (name === COUNTS_BUSY) return busySentence(taxid, 'delete', writer);
  return undefined;
}

/** A save that could not start: refused in place (the ID's lock is held), or failed. */
function refusedSave(error: unknown, taxid: number, writer: BusyWriter): Extract<SaveCountsResult, { ok: false }> {
  const name = errorName(error);
  if (name === COUNTS_IN_USE) return { ok: false, reason: 'in-use', message: inUseSentence(taxid, 'replace') };
  if (name === COUNTS_BUSY) return { ok: false, reason: 'busy', message: busySentence(taxid, 'save', writer) };
  return { ok: false, reason: 'failed', message: messageOf(error) };
}

function messageOf(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
