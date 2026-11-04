// Keeping a search alive (design §9.4, REQ-22). Browsers slow down or stop the work of a
// hidden page, mobile browsers end it, and closing the tab ends the search: LOSAT Web does
// not say that a search continues without the tab. This module
// - holds the screen wake lock while runs are active or queued, unless the user turned that off
//   ("Keep searching while you're away", on by default, Owner 2026-11-03; the choice is shown on
//   handheld devices, whose screens turn off by themselves and end the search), not released and
//   requested again between the runs of a queue;
// - asks the browser to warn before the page is left while runs are active or queued, and also
//   while the search form holds inputs or a run has results: LOSAT Web keeps nothing for the
//   next visit (REQ-20) and sends nothing to a server, so a reload loses them for good (WA-IN-06);
// - when the page becomes visible again after being hidden with runs in progress, reports
//   how long it was hidden, what became of those runs, and whether the Data worker (which
//   keeps the results) still answers.
import type { PagePort, WakeLockHandle } from '../ports/page';
import type { RunStatus } from '../domain/run';
import type { AppState } from './coordinator';
import { holdsInputs, type DraftState } from './draft';
import { Store } from './store';

export type WakeLockState = 'unsupported' | 'off' | 'requesting' | 'held' | 'failed';

export interface ResumeReport {
  /** How long the page was hidden. */
  readonly hiddenMs: number;
  /** The runs that were active or queued when the page was hidden: their status then and now. */
  readonly runs: ReadonlyArray<{ readonly number: number; readonly before: RunStatus; readonly now: RunStatus }>;
  /** Whether the Data worker answered after the page came back. */
  readonly dataWorker: 'checking' | 'responding' | 'not-responding';
}

export interface AttentionState {
  /** Keep the screen on while a search runs ("Keep searching while you're away"): on unless turned off. */
  readonly keepAwake: boolean;
  /** Whether the device is handheld (PagePort.handheld): the choice is shown there. */
  readonly handheld: boolean;
  readonly wakeLock: WakeLockState;
  readonly wakeLockError?: string;
  /** A run is active (preparing, running or finalizing). */
  readonly active: boolean;
  /** A run is active or queued; leaving the page would end it. */
  readonly busy: boolean;
  readonly resume?: ResumeReport;
}

export interface AttentionDeps {
  readonly page: PagePort;
  /** Whether keepAwake starts on (the default) or off. */
  readonly keepAwake?: boolean | undefined;
  readonly runs: Store<AppState>;
  /** The search form, whose inputs a reload would lose. */
  readonly draft?: Store<DraftState>;
  /** A call that the Data worker answers (its storage status). */
  readonly probe: () => Promise<unknown>;
  readonly now: () => number;
  /** How long the Data worker may take to answer after the page comes back. */
  readonly probeTimeoutMs?: number;
  /** Shorter hides (a quick switch of tabs) make no report. */
  readonly minHiddenMs?: number;
}

const ACTIVE: ReadonlySet<RunStatus> = new Set(['preparing', 'running', 'finalizing']);

export class Attention {
  readonly state: Store<AttentionState>;
  private lock: WakeLockHandle | undefined;
  /**
   * Counts the occasions to ask again after a refusal: runs starting after none were in progress, the
   * page shown again, the choice turned on. A refusal holds until the next one, so that each change
   * of the runs does not ask again and the refusal is not shown on and off (WA-RN-04).
   */
  private occasion = 0;
  /** The occasion at whose request the browser last refused the lock. */
  private refusedAt: number | undefined;
  private hidden: { readonly at: number; readonly runs: ReadonlyMap<string, { number: number; status: RunStatus }> } | undefined;
  private probeGeneration = 0;
  private leaveGuard: boolean | undefined;

  constructor(private readonly deps: AttentionDeps) {
    this.state = new Store<AttentionState>({
      keepAwake: deps.keepAwake ?? true,
      wakeLock: deps.page.wakeLockSupported ? 'off' : 'unsupported',
      handheld: deps.page.handheld,
      active: false,
      busy: false,
    });
    deps.runs.subscribe(() => this.update());
    deps.draft?.subscribe(() => this.updateLeaveGuard());
    deps.page.onVisibilityChange((visible) => (visible ? this.shown() : this.hide()));
    this.update();
  }

  setKeepAwake(keepAwake: boolean): void {
    if (keepAwake) this.occasion++;
    this.set({ keepAwake });
    this.update();
  }

  dismissResume(): void {
    this.set({ resume: undefined });
  }

  private update(): void {
    const runs = this.deps.runs.get().runs;
    const active = runs.some((run) => ACTIVE.has(run.status));
    const busy = active || runs.some((run) => run.status === 'queued');
    const state = this.state.get();
    if (busy && !state.busy) this.occasion++;
    if (state.active !== active || state.busy !== busy) this.set({ active, busy });
    this.updateLeaveGuard();
    this.updateWakeLock();
    // A report whose runs have all finished keeps its numbers: they describe the past.
  }

  /** On while leaving would end a search or lose results or inputs. */
  private updateLeaveGuard(): void {
    const draft = this.deps.draft?.get();
    const on =
      this.state.get().busy ||
      this.deps.runs.get().runs.some((run) => run.status === 'completed') ||
      (draft !== undefined && holdsInputs(draft));
    if (on === this.leaveGuard) return;
    this.leaveGuard = on;
    this.deps.page.setLeaveGuard(on);
  }

  private updateWakeLock(): void {
    const { keepAwake, busy, wakeLock } = this.state.get();
    if (wakeLock === 'unsupported') return;
    const wanted = keepAwake && busy && this.deps.page.isVisible();
    if (wanted && (wakeLock === 'off' || (wakeLock === 'failed' && this.refusedAt !== this.occasion))) {
      const occasion = this.occasion;
      this.set({ wakeLock: 'requesting', wakeLockError: undefined });
      this.deps.page.requestWakeLock().then(
        (lock) => {
          this.lock = lock;
          lock.onRelease(() => {
            if (this.lock !== lock) return;
            this.lock = undefined;
            this.set({ wakeLock: 'off' });
            // The browser ends the lock when the page is hidden; take it again when it is shown.
            this.updateWakeLock();
          });
          this.set({ wakeLock: 'held' });
          // The search may have ended, or the choice changed, while the lock was requested.
          this.updateWakeLock();
        },
        (error: unknown) => {
          this.refusedAt = occasion;
          this.set({ wakeLock: 'failed', wakeLockError: error instanceof Error ? error.message : String(error) });
          // An occasion came while the request was made (the page was shown again): ask once more.
          if (occasion !== this.occasion) this.updateWakeLock();
        },
      );
    } else if (!wanted && wakeLock === 'held' && this.lock !== undefined) {
      const lock = this.lock;
      this.lock = undefined;
      this.set({ wakeLock: 'off' });
      void lock.release().catch(() => undefined);
    }
  }

  private hide(): void {
    const runs = new Map<string, { number: number; status: RunStatus }>();
    for (const run of this.deps.runs.get().runs) {
      if (ACTIVE.has(run.status) || run.status === 'queued') {
        runs.set(run.snapshot.runId, { number: run.snapshot.number, status: run.status });
      }
    }
    this.hidden = runs.size === 0 ? undefined : { at: this.deps.now(), runs };
  }

  private shown(): void {
    const hidden = this.hidden;
    this.hidden = undefined;
    this.occasion++;
    this.updateWakeLock();
    if (hidden === undefined || this.deps.now() - hidden.at < (this.deps.minHiddenMs ?? 1000)) return;
    const current = new Map(this.deps.runs.get().runs.map((run) => [run.snapshot.runId, run.status]));
    const runs = [...hidden.runs].map(([runId, before]) => ({
      number: before.number,
      before: before.status,
      now: current.get(runId) ?? before.status,
    }));
    this.set({ resume: { hiddenMs: this.deps.now() - hidden.at, runs, dataWorker: 'checking' } });
    const generation = ++this.probeGeneration;
    const report = (dataWorker: ResumeReport['dataWorker']) => {
      const resume = this.state.get().resume;
      if (generation === this.probeGeneration && resume !== undefined) this.set({ resume: { ...resume, dataWorker } });
    };
    let answered = false;
    const timer = setTimeout(() => {
      if (!answered) report('not-responding');
    }, this.deps.probeTimeoutMs ?? 10_000);
    // A late answer still corrects the report.
    this.deps.probe().then(
      () => {
        answered = true;
        clearTimeout(timer);
        report('responding');
      },
      () => {
        answered = true;
        clearTimeout(timer);
        report('not-responding');
      },
    );
  }

  private set(change: { readonly [K in keyof AttentionState]?: AttentionState[K] | undefined }): void {
    const next = { ...this.state.get(), ...change };
    for (const key of Object.keys(change) as Array<keyof AttentionState>) {
      if (change[key] === undefined) delete next[key];
    }
    this.state.set(next as AttentionState);
  }
}
