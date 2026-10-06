// Keeping a search alive (design §9.4, REQ-22). Browsers slow down or stop the work of a
// hidden page, mobile browsers end it, and closing the tab ends the search: LOSAT Web does
// not say that a search continues without the tab. This module
// - holds the screen wake lock while a search runs, if the user chose it;
// - asks the browser to warn before the page is left while runs are active or queued;
// - when the page becomes visible again after being hidden with runs in progress, reports
//   how long it was hidden, what became of those runs, and whether the Data worker (which
//   keeps the results) still answers.
import type { PagePort, WakeLockHandle } from '../ports/page';
import type { RunStatus } from '../domain/run';
import type { AppState } from './coordinator';
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
  /** The user's choice: keep the screen on while a search runs. */
  readonly keepAwake: boolean;
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
  readonly runs: Store<AppState>;
  /** A call that the Data worker answers (its storage status). */
  readonly probe: () => Promise<unknown>;
  readonly now: () => number;
  /** How long the Data worker may take to answer after the page comes back. */
  readonly probeTimeoutMs?: number;
}

const ACTIVE: ReadonlySet<RunStatus> = new Set(['preparing', 'running', 'finalizing']);

export class Attention {
  readonly state: Store<AttentionState>;
  private lock: WakeLockHandle | undefined;
  private hidden: { readonly at: number; readonly runs: ReadonlyMap<string, { number: number; status: RunStatus }> } | undefined;
  private probeGeneration = 0;

  constructor(private readonly deps: AttentionDeps) {
    this.state = new Store<AttentionState>({
      keepAwake: false,
      wakeLock: deps.page.wakeLockSupported ? 'off' : 'unsupported',
      active: false,
      busy: false,
    });
    deps.runs.subscribe(() => this.update());
    deps.page.onVisibilityChange((visible) => (visible ? this.shown() : this.hide()));
    this.update();
  }

  setKeepAwake(keepAwake: boolean): void {
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
    if (state.active !== active || state.busy !== busy) this.set({ active, busy });
    this.deps.page.setLeaveGuard(busy);
    this.updateWakeLock();
    // A report whose runs have all finished keeps its numbers: they describe the past.
  }

  private updateWakeLock(): void {
    const { keepAwake, active, wakeLock } = this.state.get();
    if (wakeLock === 'unsupported') return;
    const wanted = keepAwake && active && this.deps.page.isVisible();
    if (wanted && (wakeLock === 'off' || wakeLock === 'failed')) {
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
        (error: unknown) => this.set({ wakeLock: 'failed', wakeLockError: error instanceof Error ? error.message : String(error) }),
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
    this.updateWakeLock();
    if (hidden === undefined) return;
    const current = new Map(this.deps.runs.get().runs.map((run) => [run.snapshot.runId, run.status]));
    const runs = [...hidden.runs].map(([runId, before]) => ({
      number: before.number,
      before: before.status,
      now: current.get(runId) ?? before.status,
    }));
    this.set({ resume: { hiddenMs: this.deps.now() - hidden.at, runs, dataWorker: 'checking' } });
    const generation = ++this.probeGeneration;
    const timeout = new Promise<'not-responding'>((resolve) =>
      setTimeout(() => resolve('not-responding'), this.deps.probeTimeoutMs ?? 10_000),
    );
    const answer = this.deps.probe().then(
      () => 'responding' as const,
      () => 'not-responding' as const,
    );
    void Promise.race([answer, timeout]).then((dataWorker) => {
      const resume = this.state.get().resume;
      if (generation === this.probeGeneration && resume !== undefined) this.set({ resume: { ...resume, dataWorker } });
    });
  }

  private set(change: { readonly [K in keyof AttentionState]?: AttentionState[K] | undefined }): void {
    const next = { ...this.state.get(), ...change };
    for (const key of Object.keys(change) as Array<keyof AttentionState>) {
      if (change[key] === undefined) delete next[key];
    }
    this.state.set(next as AttentionState);
  }
}
