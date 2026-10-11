// The release this tab runs, its offline use and the explicit update (plan §5.9, design §14).
// - Offline: once the app has started and the browser is idle, the site's service worker
//   stores this tab's release (the pinned one, which the boot page started) and the site's
//   shell after checking every file, and only then, and only while the pin read back from
//   storage names this release (an offline start opens the pinned release), does the screen
//   say "Offline ready".
// - Update: the site index is read after start-up and then hourly while the page is visible
//   and online (never while offline). A release other than this tab's is offered; nothing
//   changes until the user chooses. The user reviews what a reload does to this tab's work,
//   has the worker download and check the new release, sees again what the reload does, and
//   only then switches: the pin is written, a waiting worker is told to take over, and the
//   page reloads. A failed download or a pin that cannot be written leaves this tab, the pin
//   and the stored release as they were. Other tabs are told of a switch and keep their
//   release until they are reloaded.
// - Every page holds its release's lock for its life, so the worker keeps an older release
//   while a tab still runs it. From the start of a download until the tab switches, cancels or
//   closes, it also holds the new release's lock, so that another tab's cleanup cannot remove
//   the release that was checked; the switch still checks that it is stored whole.
import type { RunStatus } from '../domain/run';
import type { OfflinePort, PrepareProgress, PrepareResult, ReleaseChannel, ReleasePage, ReleasePin } from '../ports/offline';
import type { PagePort } from '../ports/page';
import type { SitePort } from '../ports/site';
import type { AppState } from './coordinator';
import { Store } from './store';

export type OfflineState =
  /** Before the preparation starts (it waits until the app is idle). */
  | { readonly kind: 'waiting' }
  /** No offline use: the browser has no service workers (no reason), or the reason given. */
  | { readonly kind: 'unsupported'; readonly reason?: string }
  | { readonly kind: 'preparing'; readonly progress?: PrepareProgress }
  /** This tab's release is stored, and the pin names it: an offline start opens it. */
  | { readonly kind: 'ready' }
  /**
   * This tab's release is stored, but the pin does not name it: `pin` is the release that
   * another tab switched to; without `pin`, the browser did not keep the choice of release,
   * which an offline start needs.
   */
  | { readonly kind: 'stored'; readonly pin?: string }
  | { readonly kind: 'failed'; readonly reason: string };

export type UpdateState =
  | { readonly kind: 'none' }
  /** The site serves `latest`, which is not this tab's release. */
  | { readonly kind: 'available'; readonly latest: string }
  /** "Not now", in this tab only: offered again when the site serves yet another release. */
  | { readonly kind: 'dismissed'; readonly latest: string }
  /** The screen says what a switch does to this tab's work. */
  | { readonly kind: 'reviewing'; readonly latest: string }
  | { readonly kind: 'downloading'; readonly latest: string; readonly progress?: PrepareProgress }
  /** Every file of `latest` is stored and checked. */
  | { readonly kind: 'downloaded'; readonly latest: string }
  | { readonly kind: 'failed'; readonly latest: string; readonly reason: string }
  /** At the switch, `latest` was no longer stored whole (the browser removed it): nothing was switched. */
  | { readonly kind: 'removed'; readonly latest: string }
  /** The pin could not be written: the switch would not last, so the page was not reloaded. */
  | { readonly kind: 'not-kept'; readonly latest: string; readonly reason: string }
  /** The pin names `latest`; the reload was asked for (the browser may have asked the user to stay). */
  | { readonly kind: 'switched'; readonly latest: string };

export interface ReleaseState {
  /** The release this tab runs. */
  readonly release: string;
  readonly offline: OfflineState;
  readonly update: UpdateState;
  /** The release that another tab switched the site to; this tab keeps its own until it is reloaded. */
  readonly elsewhere?: string;
}

export interface ReleaseDeps {
  /** The release this tab runs (`__LOSAT_RELEASE_ID__`). */
  readonly release: string;
  /**
   * Whether the boot page started this tab from its release directory (`/r/<id>/`). A page that
   * was not (the development server, a test harness) has no offline use and offers no update.
   */
  readonly served: boolean;
  readonly offline: OfflinePort;
  readonly pin: ReleasePin;
  readonly channel: ReleaseChannel;
  readonly site: SitePort;
  readonly page: Pick<PagePort, 'isVisible' | 'onVisibilityChange'>;
  readonly host: ReleasePage;
  /** The runs: the preparation waits while a run is being prepared (its engine loads). */
  readonly runs: Store<AppState>;
  /** Time between two checks of the site index. */
  readonly checkEveryMs?: number;
  /**
   * Told when the releases that this browser stores changed: this tab's was stored (and the
   * cleanup after it ran), or an update was downloaded. The storage estimate counts them, and
   * is read again only on such events (screen review 1 L2).
   */
  readonly onStored?: () => void;
}

/** What a reload does to this tab's work, from the app's state (the review screen's counts). */
export interface SwitchImpact {
  /** Runs being prepared, running or finalizing: stopped by the reload. */
  readonly running: number;
  readonly queued: number;
  /** Completed runs: lost unless saved in a session file. */
  readonly completed: number;
  /** The candidate tray, and the candidates with a note: lost unless saved in a session file. */
  readonly candidates: number;
  readonly notes: number;
}

const RUNNING: ReadonlySet<RunStatus> = new Set(['preparing', 'running', 'finalizing']);
/** The update states in which this tab holds the lock of `latest`: it is being stored, is stored, or is the pin. */
const HOLDING: ReadonlySet<UpdateState['kind']> = new Set(['downloading', 'downloaded', 'not-kept', 'switched']);
const HOUR_MS = 60 * 60 * 1000;

export function switchImpact(
  runs: ReadonlyArray<{ readonly status: RunStatus }>,
  candidates: ReadonlyArray<{ readonly note: string }>,
): SwitchImpact {
  return {
    running: runs.filter((run) => RUNNING.has(run.status)).length,
    queued: runs.filter((run) => run.status === 'queued').length,
    completed: runs.filter((run) => run.status === 'completed').length,
    candidates: candidates.length,
    notes: candidates.filter((candidate) => candidate.note !== '').length,
  };
}

const message = (error: unknown): string => (error instanceof Error ? error.message : String(error));

export class ReleaseManager {
  readonly state: Store<ReleaseState>;
  private started = false;
  private connected: Promise<PrepareResult> | undefined;
  private checking = false;
  /** A check fell due while the page was hidden or offline: it runs when the page is visible and online again. */
  private due = false;
  private timer: ReturnType<typeof setTimeout> | undefined;
  /** The lock of the update's release that this tab holds (HOLDING), and how to let it go. */
  private updateHold: { readonly release: string; readonly letGo: () => void } | undefined;
  private switching = false;

  constructor(private readonly deps: ReleaseDeps) {
    let offline: OfflineState = { kind: 'waiting' };
    if (!deps.served) offline = { kind: 'unsupported', reason: 'this page is not served as a release of the site' };
    else if (!deps.offline.supported) offline = { kind: 'unsupported' };
    this.state = new Store<ReleaseState>({ release: deps.release, offline, update: { kind: 'none' } });
  }

  /** Holds this tab's release, listens to the other tabs, and once the app is idle checks for an update and prepares offline use. */
  start(): void {
    if (this.started || !this.deps.served) return;
    this.started = true;
    this.deps.offline.holdRelease(this.deps.release);
    this.deps.channel.onSwitch((release) => this.switchedElsewhere(release));
    this.deps.page.onVisibilityChange((visible) => {
      if (visible) this.checkIfDue();
    });
    this.deps.host.onOnline(() => this.checkIfDue());
    void this.deps.host.whenIdle().then(async () => {
      void this.check();
      this.schedule();
      if (this.state.get().offline.kind !== 'waiting') return;
      await this.noRunPreparing();
      await this.prepareOwn();
    });
  }

  /** "Try again" after the preparation of offline use failed. */
  retryOffline(): void {
    if (this.state.get().offline.kind === 'failed') void this.prepareOwn();
  }

  /** Reads the site index now, unless the page is offline (then the check waits until it is online). */
  async check(): Promise<void> {
    if (this.checking || !this.deps.served) return;
    if (!this.deps.host.online()) {
      this.due = true;
      return;
    }
    this.checking = true;
    this.due = false;
    try {
      const index = await this.deps.site.index();
      if (index !== undefined) this.noticeLatest(index.latest);
    } finally {
      this.checking = false;
    }
  }

  /** "Review update": the screen says what the switch does to this tab's work. */
  review(): void {
    const update = this.state.get().update;
    if (update.kind === 'available' || update.kind === 'dismissed') this.setUpdate({ kind: 'reviewing', latest: update.latest });
  }

  /** "Not now": this tab stops offering this release. */
  notNow(): void {
    const update = this.state.get().update;
    if (
      update.kind === 'available' ||
      update.kind === 'reviewing' ||
      update.kind === 'failed' ||
      update.kind === 'removed' ||
      update.kind === 'not-kept' ||
      update.kind === 'downloaded'
    ) {
      this.setUpdate({ kind: 'dismissed', latest: update.latest });
    }
  }

  /** "Cancel": back to the offer; nothing has changed. */
  cancel(): void {
    const update = this.state.get().update;
    if (update.kind === 'reviewing' || update.kind === 'failed' || update.kind === 'removed' || update.kind === 'not-kept' || update.kind === 'downloaded') {
      this.setUpdate({ kind: 'available', latest: update.latest });
    }
  }

  /**
   * Whether the switch needs no download first: only in a browser without service workers (or
   * where they are blocked), where nothing can be stored or checked. Where they exist, a worker
   * that cannot be registered is a failure of offline use, and the update still downloads and
   * checks first.
   */
  switchesDirectly(): boolean {
    return !this.deps.offline.supported;
  }

  /** "Download and check": the worker stores the new release after checking every file. */
  async download(): Promise<void> {
    const update = this.state.get().update;
    if ((update.kind !== 'reviewing' && update.kind !== 'failed' && update.kind !== 'removed') || this.switchesDirectly()) return;
    const latest = update.latest;
    const current = () => {
      const now = this.state.get().update;
      return now.kind === 'downloading' && now.latest === latest;
    };
    this.setUpdate({ kind: 'downloading', latest });
    const connected = await this.connect();
    // The registration showed that service workers are blocked: the switch is offered without a download.
    if (!connected.ok && this.switchesDirectly()) {
      if (current()) this.setUpdate({ kind: 'reviewing', latest });
      return;
    }
    const result = connected.ok
      ? await this.deps.offline
          .prepare(latest, (progress) => {
            if (current()) this.setUpdate({ kind: 'downloading', latest, progress });
          })
          .catch((error: unknown) => ({ ok: false as const, reason: message(error) }))
      : connected;
    if (!current()) return;
    this.setUpdate(result.ok ? { kind: 'downloaded', latest } : { kind: 'failed', latest, reason: result.reason });
    if (result.ok) this.deps.onStored?.();
  }

  /**
   * "Switch to release <id> and reload": checks that the release is still stored whole, writes
   * the pin, tells the other tabs and a waiting worker, and reloads.
   */
  async switchNow(): Promise<void> {
    const update = this.state.get().update;
    const direct = update.kind === 'reviewing' && this.switchesDirectly();
    if ((update.kind !== 'downloaded' && !direct) || this.switching) return;
    const latest = update.latest;
    if (!direct) {
      // The lock held since the download keeps other tabs' cleanups away, but the browser may
      // still have removed the release (storage pressure, site data cleared): a pin to a
      // release that is not stored would not start offline.
      this.switching = true;
      const stored = await this.deps.offline.stored(latest).catch(() => false);
      this.switching = false;
      const now = this.state.get().update;
      if (now.kind !== 'downloaded' || now.latest !== latest) return;
      if (!stored) {
        this.setUpdate({ kind: 'removed', latest });
        return;
      }
    }
    const written = this.deps.pin.write(latest);
    if (!written.ok) {
      this.setUpdate({ kind: 'not-kept', latest, reason: written.reason });
      return;
    }
    this.deps.channel.announce(latest);
    if (this.deps.offline.supported) await this.deps.offline.activateWaiting().catch(() => false);
    this.setUpdate({ kind: 'switched', latest });
    this.deps.host.reload();
  }

  /** "Reload now", when the browser kept the page after a switch (the user chose to stay). */
  reloadNow(): void {
    if (this.state.get().update.kind === 'switched') this.deps.host.reload();
  }

  private async prepareOwn(): Promise<void> {
    this.setOffline({ kind: 'preparing' });
    if (!(await this.connect()).ok) return;
    const release = this.deps.release;
    const result = await this.deps.offline
      .prepare(release, (progress) => this.setOffline({ kind: 'preparing', progress }))
      .catch((error: unknown) => ({ ok: false as const, reason: message(error) }));
    if (!result.ok) {
      this.setOffline({ kind: 'failed', reason: result.reason });
      return;
    }
    const pin = this.deps.pin.read();
    this.setOffline(this.storedWith(pin));
    // Told as the status turns, and again once the cleanup may have removed older releases.
    this.deps.onStored?.();
    await this.deps.offline.cleanup(pin === undefined || pin === release ? [release] : [release, pin]).catch(() => undefined);
    this.deps.onStored?.();
  }

  /** The offline state of this tab's stored release under `pin` (read back from storage). */
  private storedWith(pin: string | undefined): OfflineState {
    if (pin === this.deps.release) return { kind: 'ready' };
    return pin === undefined ? { kind: 'stored' } : { kind: 'stored', pin };
  }

  /** "Offline ready" follows the pin once the release is stored (another tab may switch it). */
  private followPin(pin: string | undefined): void {
    const kind = this.state.get().offline.kind;
    if (kind === 'ready' || kind === 'stored') this.setOffline(this.storedWith(pin));
  }

  /**
   * Connects to the worker, once it has worked. A registration that fails or times out where
   * service workers exist is a failure of offline use ("Try again" connects again), never
   * "unsupported", which would offer a switch that nothing checked; only a registration that
   * shows service workers blocked (the port's `supported` turns false) is "unsupported".
   */
  private connect(): Promise<PrepareResult> {
    this.connected ??= this.deps.offline.connect().then(
      () => ({ ok: true as const }),
      (error: unknown) => {
        this.connected = undefined;
        if (!this.deps.offline.supported) {
          this.setOffline({ kind: 'unsupported', reason: message(error) });
          return { ok: false as const, reason: message(error) };
        }
        const reason = `the offline worker could not be started (${message(error)})`;
        this.setOffline({ kind: 'failed', reason });
        return { ok: false as const, reason };
      },
    );
    return this.connected;
  }

  private noRunPreparing(): Promise<void> {
    const preparing = () => this.deps.runs.get().runs.some((run) => run.status === 'preparing');
    if (!preparing()) return Promise.resolve();
    return new Promise((resolve) => {
      const stop = this.deps.runs.subscribe(() => {
        if (preparing()) return;
        stop();
        resolve();
      });
    });
  }

  private schedule(): void {
    clearTimeout(this.timer);
    this.timer = setTimeout(() => {
      this.schedule();
      if (this.deps.page.isVisible() && this.deps.host.online()) void this.check();
      else this.due = true;
    }, this.deps.checkEveryMs ?? HOUR_MS);
  }

  private checkIfDue(): void {
    if (this.due && this.deps.page.isVisible() && this.deps.host.online()) void this.check();
  }

  /** What the site's latest release and the pin mean for this tab. */
  private noticeLatest(latest: string): void {
    const { release, update } = this.state.get();
    const pin = this.deps.pin.read();
    // A pin that names another release was written by another tab (its message may have been
    // missed while this tab was frozen); after this tab's own switch it is this tab's.
    if (pin !== undefined && pin !== release && update.kind !== 'switched') this.set({ elsewhere: pin });
    if (update.kind !== 'switched') this.followPin(pin);
    const elsewhere = this.state.get().elsewhere;
    const offered = latest !== release && latest !== elsewhere;
    if (update.kind === 'none') {
      if (offered) this.setUpdate({ kind: 'available', latest });
    } else if (update.kind === 'available' || update.kind === 'dismissed') {
      if (!offered) this.setUpdate({ kind: 'none' });
      else if (update.latest !== latest) this.setUpdate({ kind: 'available', latest });
    }
    // Otherwise the user is in the middle of an update: a check does not change it.
  }

  private switchedElsewhere(release: string): void {
    const own = release === this.deps.release;
    this.set({ elsewhere: own ? undefined : release });
    this.followPin(release);
    const update = this.state.get().update;
    if (update.kind === 'none') return;
    // The offer, or this tab's own switch to another release, no longer holds: the pin now names `release`.
    if (((update.kind === 'available' || update.kind === 'dismissed' || update.kind === 'reviewing') && update.latest === release) || update.kind === 'switched') {
      this.setUpdate({ kind: 'none' });
    }
  }

  private setOffline(offline: OfflineState): void {
    this.set({ offline });
  }

  private setUpdate(update: UpdateState): void {
    const hold = HOLDING.has(update.kind) && 'latest' in update ? update.latest : undefined;
    if (this.updateHold?.release !== hold) {
      // Taken before the old one is let go, so that a release never goes unheld between two states.
      const next = hold === undefined ? undefined : { release: hold, letGo: this.deps.offline.holdRelease(hold) };
      this.updateHold?.letGo();
      this.updateHold = next;
    }
    this.set({ update });
  }

  private set(change: { readonly [K in keyof ReleaseState]?: ReleaseState[K] | undefined }): void {
    const next = { ...this.state.get(), ...change };
    for (const key of Object.keys(change) as Array<keyof ReleaseState>) {
      if (change[key] === undefined) delete next[key];
    }
    this.state.set(next as ReleaseState);
  }
}
