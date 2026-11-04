import { describe, expect, it, vi } from 'vitest';
import { Attention, type AttentionDeps } from '../../src/application/attention';
import type { AppState, RunView } from '../../src/application/coordinator';
import type { DraftSource, DraftState, RoleDraft } from '../../src/application/draft';
import { Store } from '../../src/application/store';
import type { RunStatus } from '../../src/domain/run';
import type { PagePort, WakeLockHandle } from '../../src/ports/page';

class FakeLock implements WakeLockHandle {
  releaseCalls = 0;
  private listener: (() => void) | undefined;
  release = vi.fn(async () => {
    this.releaseCalls++;
    this.listener?.();
  });
  onRelease(listener: () => void): void {
    this.listener = listener;
  }
  /** The browser ends the lock. */
  endByBrowser(): void {
    this.listener?.();
  }
}

class FakePage implements PagePort {
  visible = true;
  leaveGuards: boolean[] = [];
  locks: FakeLock[] = [];
  nextRequestError: Error | undefined;
  handheld = false;
  requestWakeLock = vi.fn(async (): Promise<WakeLockHandle> => {
    if (this.nextRequestError !== undefined) throw this.nextRequestError;
    const lock = new FakeLock();
    this.locks.push(lock);
    return lock;
  });
  private readonly listeners = new Set<(visible: boolean) => void>();

  constructor(readonly wakeLockSupported = true) {}

  setLeaveGuard(on: boolean): void {
    this.leaveGuards.push(on);
  }
  isVisible(): boolean {
    return this.visible;
  }
  onVisibilityChange(listener: (visible: boolean) => void): () => void {
    this.listeners.add(listener);
    return () => this.listeners.delete(listener);
  }
  setVisible(visible: boolean): void {
    this.visible = visible;
    for (const listener of this.listeners) listener(visible);
  }
  get guard(): boolean | undefined {
    return this.leaveGuards.at(-1);
  }
}

const run = (runId: string, number: number, status: RunStatus): RunView =>
  ({ snapshot: { runId, number }, status, record: {} }) as unknown as RunView;

const flush = async (): Promise<void> => {
  for (let i = 0; i < 5; i++) await Promise.resolve();
};
const sleep = (ms: number): Promise<void> => new Promise((resolve) => setTimeout(resolve, ms));

function setup(
  options: {
    supported?: boolean;
    handheld?: boolean;
    /** The initial choice; the default (on) unless given. */
    keepAwake?: boolean;
    runs?: RunView[];
    probe?: AttentionDeps['probe'];
    probeTimeoutMs?: number;
    minHiddenMs?: number;
    draft?: Store<DraftState>;
  } = {},
) {
  const page = new FakePage(options.supported ?? true);
  page.handheld = options.handheld ?? false;
  const runs = new Store<AppState>({ runs: options.runs ?? [] });
  const clock = { now: 1_000 };
  const attention = new Attention({
    page,
    runs,
    ...(options.draft === undefined ? {} : { draft: options.draft }),
    keepAwake: options.keepAwake,
    probe: options.probe ?? (() => Promise.resolve()),
    now: () => clock.now,
    probeTimeoutMs: options.probeTimeoutMs ?? 5_000,
    minHiddenMs: options.minHiddenMs ?? 0,
  });
  const setRuns = (next: RunView[]) => runs.set({ runs: next });
  return { page, runs, clock, attention, setRuns, state: () => attention.state.get() };
}

describe('leave guard', () => {
  it('is off without runs', () => {
    const { page, state } = setup();
    expect(page.guard).toBe(false);
    expect(state()).toMatchObject({ active: false, busy: false });
  });

  it.each(['queued', 'preparing', 'running', 'finalizing'] as const)('is on while a run is %s', (status) => {
    const { page, state } = setup({ runs: [run('a', 1, status)] });
    expect(page.guard).toBe(true);
    expect(state().busy).toBe(true);
    expect(state().active).toBe(status !== 'queued');
  });

  it('follows the runs: on while any run is in progress, off when all have ended without results', () => {
    const { page, setRuns, state } = setup({ runs: [run('a', 1, 'running'), run('b', 2, 'queued')] });
    expect(page.guard).toBe(true);
    setRuns([run('a', 1, 'cancelled'), run('b', 2, 'queued')]);
    expect(page.guard).toBe(true);
    expect(state()).toMatchObject({ active: false, busy: true });
    setRuns([run('a', 1, 'cancelled'), run('b', 2, 'running')]);
    expect(state()).toMatchObject({ active: true, busy: true });
    setRuns([run('a', 1, 'cancelled'), run('b', 2, 'failed')]);
    expect(page.guard).toBe(false);
    expect(state()).toMatchObject({ active: false, busy: false });
    setRuns([run('a', 1, 'cancelled')]);
    expect(page.guard).toBe(false);
  });

  // Nothing is kept for the next visit (REQ-20), so a reload loses results and inputs for good (WA-IN-06).
  it('stays on while a run has results', () => {
    const { page, setRuns, state } = setup({ runs: [run('a', 1, 'running')] });
    setRuns([run('a', 1, 'completed')]);
    expect(page.guard).toBe(true);
    expect(state()).toMatchObject({ active: false, busy: false });
  });

  it('is on while the search form holds a pasted text or a source of either role', () => {
    const role = (change: Partial<RoleDraft> = {}): RoleDraft => ({ paste: '', sources: [], mode: 'combined', ...change });
    const form = (query: RoleDraft, subject: RoleDraft) => ({ query, subject }) as unknown as DraftState;
    const draft = new Store<DraftState>(form(role(), role()));
    const { page } = setup({ draft });
    expect(page.guard).toBe(false);
    draft.set(form(role({ paste: '>q\nACGT\n' }), role()));
    expect(page.guard).toBe(true);
    draft.set(form(role(), role()));
    expect(page.guard).toBe(false);
    draft.set(form(role(), role({ sources: [{ key: 'file-1' } as DraftSource] })));
    expect(page.guard).toBe(true);
    // A box of accessions only has no source until it is fetched; it is an input too.
    draft.set(form(role({ paste: 'NC_045512.2' }), role()));
    expect(page.guard).toBe(true);
    // A form change that leaves the guard as it is does not set it again.
    const calls = page.leaveGuards.length;
    draft.set(form(role({ paste: 'NC_045512.3' }), role()));
    expect(page.leaveGuards).toHaveLength(calls);
  });
});

describe('wake lock', () => {
  it("is 'unsupported' when the page cannot hold one, and is never requested", async () => {
    const { page, attention, state } = setup({ supported: false, runs: [run('a', 1, 'running')] });
    expect(state().wakeLock).toBe('unsupported');
    attention.setKeepAwake(true);
    await flush();
    expect(state().wakeLock).toBe('unsupported');
    expect(page.requestWakeLock).not.toHaveBeenCalled();
  });

  it('is on by default, and requested as soon as a run is active', async () => {
    const { page, state } = setup({ runs: [run('a', 1, 'running')] });
    await flush();
    expect(state()).toMatchObject({ keepAwake: true, wakeLock: 'held', handheld: false });
    expect(page.requestWakeLock).toHaveBeenCalledTimes(1);
  });

  it('reports a handheld device (where the choice is shown) from the page', () => {
    expect(setup({ handheld: true }).state().handheld).toBe(true);
    expect(setup({ handheld: true, supported: false }).state()).toMatchObject({ handheld: true, wakeLock: 'unsupported' });
  });

  it('is off, and not requested, while keepAwake is false', async () => {
    const { page, state } = setup({ keepAwake: false, runs: [run('a', 1, 'running')] });
    await flush();
    expect(state()).toMatchObject({ keepAwake: false, wakeLock: 'off' });
    expect(page.requestWakeLock).not.toHaveBeenCalled();
  });

  it('is not requested without runs, and is held from a queued run to the end of the queue', async () => {
    const { page, attention, setRuns, state } = setup({ keepAwake: false, runs: [run('a', 1, 'completed')] });
    attention.setKeepAwake(true);
    await flush();
    expect(page.requestWakeLock).not.toHaveBeenCalled();
    expect(state().wakeLock).toBe('off');
    setRuns([run('a', 1, 'completed'), run('b', 2, 'queued'), run('c', 3, 'queued')]);
    await flush();
    expect(page.requestWakeLock).toHaveBeenCalledTimes(1);
    // Between the runs of the queue the lock is kept, not released and requested again.
    setRuns([run('a', 1, 'completed'), run('b', 2, 'completed'), run('c', 3, 'preparing')]);
    await flush();
    expect(page.requestWakeLock).toHaveBeenCalledTimes(1);
    expect(page.locks[0]!.releaseCalls).toBe(0);
    setRuns([run('a', 1, 'completed'), run('b', 2, 'completed'), run('c', 3, 'completed')]);
    await flush();
    expect(page.locks[0]!.releaseCalls).toBe(1);
    expect(state().wakeLock).toBe('off');
  });

  it('is requested when keepAwake is turned on during a run, and released when the run completes', async () => {
    const { page, attention, setRuns, state } = setup({ keepAwake: false, runs: [run('a', 1, 'running')] });
    attention.setKeepAwake(true);
    expect(state()).toMatchObject({ keepAwake: true, wakeLock: 'requesting' });
    await flush();
    expect(state().wakeLock).toBe('held');
    expect(page.requestWakeLock).toHaveBeenCalledTimes(1);
    const lock = page.locks[0]!;
    expect(lock.releaseCalls).toBe(0);

    setRuns([run('a', 1, 'completed')]);
    await flush();
    expect(lock.releaseCalls).toBe(1);
    expect(state().wakeLock).toBe('off');
    expect(page.requestWakeLock).toHaveBeenCalledTimes(1);
  });

  it('is released when keepAwake is turned off', async () => {
    const { page, attention, state } = setup({ runs: [run('a', 1, 'running')] });
    attention.setKeepAwake(true);
    await flush();
    attention.setKeepAwake(false);
    await flush();
    expect(page.locks[0]!.releaseCalls).toBe(1);
    expect(state()).toMatchObject({ keepAwake: false, wakeLock: 'off' });
  });

  it("is 'failed' with the error message when the request is rejected", async () => {
    const { page, attention, state } = setup({ keepAwake: false, runs: [run('a', 1, 'running')] });
    page.nextRequestError = new Error('Battery is low');
    attention.setKeepAwake(true);
    await flush();
    expect(state()).toMatchObject({ wakeLock: 'failed', wakeLockError: 'Battery is low' });
  });

  it('writes a non-Error rejection as a string', async () => {
    const { page, attention, state } = setup({ keepAwake: false, runs: [run('a', 1, 'running')] });
    page.requestWakeLock.mockRejectedValueOnce('refused');
    attention.setKeepAwake(true);
    await flush();
    expect(state()).toMatchObject({ wakeLock: 'failed', wakeLockError: 'refused' });
  });

  it('is not requested again while the same runs go on after a failure, so the failure stays shown (WA-RN-04)', async () => {
    const { page, setRuns, state } = setup();
    page.nextRequestError = new Error('Wake Lock permission request denied');
    setRuns([run('a', 1, 'queued'), run('b', 2, 'queued')]);
    await flush();
    const seen: string[] = [];
    const record = () => seen.push(`${state().wakeLock}:${state().wakeLockError ?? ''}`);
    record();
    // Two runs one after the other: every phase changes the runs.
    for (const next of [
      [run('a', 1, 'preparing'), run('b', 2, 'queued')],
      [run('a', 1, 'running'), run('b', 2, 'queued')],
      [run('a', 1, 'finalizing'), run('b', 2, 'queued')],
      [run('a', 1, 'completed'), run('b', 2, 'preparing')],
      [run('a', 1, 'completed'), run('b', 2, 'running')],
      [run('a', 1, 'completed'), run('b', 2, 'finalizing')],
    ]) {
      setRuns(next);
      await flush();
      record();
    }
    expect(page.requestWakeLock).toHaveBeenCalledTimes(1);
    expect(new Set(seen)).toEqual(new Set(['failed:Wake Lock permission request denied']));
  });

  it('is requested again after a failure when runs start again, and the error is cleared', async () => {
    const { page, attention, setRuns, state } = setup({ keepAwake: false, runs: [run('a', 1, 'running')] });
    page.nextRequestError = new Error('Battery is low');
    attention.setKeepAwake(true);
    await flush();
    expect(state().wakeLock).toBe('failed');
    page.nextRequestError = undefined;
    setRuns([run('a', 1, 'finalizing')]);
    await flush();
    expect(state().wakeLock).toBe('failed');
    expect(page.requestWakeLock).toHaveBeenCalledTimes(1);
    setRuns([run('a', 1, 'completed')]);
    setRuns([run('a', 1, 'completed'), run('b', 2, 'running')]);
    await flush();
    expect(state().wakeLock).toBe('held');
    expect(state().wakeLockError).toBeUndefined();
    expect(page.requestWakeLock).toHaveBeenCalledTimes(2);
  });

  it('is requested once more when the page was hidden and shown while a request that the browser refuses was made', async () => {
    const { page, setRuns, state } = setup();
    let refuse!: (error: Error) => void;
    page.requestWakeLock.mockImplementationOnce(() => new Promise((_, reject) => (refuse = reject)));
    setRuns([run('a', 1, 'running')]);
    expect(state().wakeLock).toBe('requesting');
    page.setVisible(false);
    page.setVisible(true);
    refuse(new Error('Battery is low'));
    await flush();
    expect(page.requestWakeLock).toHaveBeenCalledTimes(2);
    expect(state().wakeLock).toBe('held');
  });

  it('is requested again after a failure when the page is shown again, or the choice is turned on again', async () => {
    const { page, attention, setRuns, state } = setup();
    page.nextRequestError = new Error('Battery is low');
    setRuns([run('a', 1, 'running')]);
    await flush();
    expect(state().wakeLock).toBe('failed');
    page.setVisible(false);
    page.setVisible(true);
    await flush();
    expect(page.requestWakeLock).toHaveBeenCalledTimes(2);
    expect(state().wakeLock).toBe('failed');
    page.nextRequestError = undefined;
    attention.setKeepAwake(false);
    attention.setKeepAwake(true);
    await flush();
    expect(page.requestWakeLock).toHaveBeenCalledTimes(3);
    expect(state().wakeLock).toBe('held');
  });

  it('is requested again when the browser releases it while the run is active and the page visible', async () => {
    const { page, attention, state } = setup({ keepAwake: false, runs: [run('a', 1, 'running')] });
    attention.setKeepAwake(true);
    await flush();
    expect(page.requestWakeLock).toHaveBeenCalledTimes(1);
    page.locks[0]!.endByBrowser();
    await flush();
    expect(page.requestWakeLock).toHaveBeenCalledTimes(2);
    expect(state().wakeLock).toBe('held');
    expect(page.locks).toHaveLength(2);
  });

  it('is not requested again when the browser releases it while the page is hidden', async () => {
    const { page, attention, state } = setup({ keepAwake: false, runs: [run('a', 1, 'running')] });
    attention.setKeepAwake(true);
    await flush();
    page.visible = false;
    page.locks[0]!.endByBrowser();
    await flush();
    expect(page.requestWakeLock).toHaveBeenCalledTimes(1);
    expect(state().wakeLock).toBe('off');
  });

  it('is not requested while the page is hidden, and is requested when it is shown', async () => {
    const { page, attention, state } = setup({ keepAwake: false, runs: [run('a', 1, 'running')] });
    page.setVisible(false);
    attention.setKeepAwake(true);
    await flush();
    expect(page.requestWakeLock).not.toHaveBeenCalled();
    expect(state().wakeLock).toBe('off');
    page.setVisible(true);
    await flush();
    expect(page.requestWakeLock).toHaveBeenCalledTimes(1);
    expect(state().wakeLock).toBe('held');
  });
});

describe('resume report', () => {
  it('reports how long the page was hidden and what became of the runs', async () => {
    const { page, clock, setRuns, state } = setup({
      runs: [run('a', 1, 'running'), run('b', 2, 'queued'), run('c', 3, 'completed')],
    });
    page.setVisible(false);
    expect(state().resume).toBeUndefined();
    clock.now += 42_000;
    setRuns([run('a', 1, 'completed'), run('b', 2, 'running'), run('c', 3, 'completed')]);
    page.setVisible(true);
    expect(state().resume).toEqual({
      hiddenMs: 42_000,
      runs: [
        { number: 1, before: 'running', now: 'completed' },
        { number: 2, before: 'queued', now: 'running' },
      ],
      dataWorker: 'checking',
    });
    await flush();
    expect(state().resume?.dataWorker).toBe('responding');
    expect(state().resume?.hiddenMs).toBe(42_000);
  });

  it('keeps the earlier status for a run that is no longer in the list', () => {
    const { page, setRuns, state } = setup({ runs: [run('a', 1, 'running')] });
    page.setVisible(false);
    setRuns([]);
    page.setVisible(true);
    expect(state().resume?.runs).toEqual([{ number: 1, before: 'running', now: 'running' }]);
  });

  it("is 'not-responding' when the probe rejects", async () => {
    const { page, state } = setup({ runs: [run('a', 1, 'running')], probe: () => Promise.reject(new Error('gone')) });
    page.setVisible(false);
    page.setVisible(true);
    expect(state().resume?.dataWorker).toBe('checking');
    await flush();
    expect(state().resume?.dataWorker).toBe('not-responding');
  });

  it("is 'not-responding' when the probe does not answer within probeTimeoutMs", async () => {
    const { page, state } = setup({
      runs: [run('a', 1, 'running')],
      probe: () => new Promise(() => undefined),
      probeTimeoutMs: 20,
    });
    page.setVisible(false);
    page.setVisible(true);
    expect(state().resume?.dataWorker).toBe('checking');
    await sleep(5);
    expect(state().resume?.dataWorker).toBe('checking');
    await sleep(80);
    expect(state().resume?.dataWorker).toBe('not-responding');
  });

  it('is not made when no run was in progress as the page was hidden', async () => {
    const probe = vi.fn(() => Promise.resolve());
    const { page, state } = setup({ runs: [run('a', 1, 'completed'), run('b', 2, 'failed')], probe });
    page.setVisible(false);
    page.setVisible(true);
    await flush();
    expect(state().resume).toBeUndefined();
    expect(probe).not.toHaveBeenCalled();
  });

  it('is not made when there are no runs, nor when the page is shown without having been hidden', async () => {
    const { page, state } = setup({ runs: [run('a', 1, 'running')] });
    page.setVisible(true);
    await flush();
    expect(state().resume).toBeUndefined();
    const empty = setup();
    empty.page.setVisible(false);
    empty.page.setVisible(true);
    expect(empty.state().resume).toBeUndefined();
  });

  it('is not made again by a second show after dismissResume', async () => {
    const { page, attention, state } = setup({ runs: [run('a', 1, 'running')] });
    page.setVisible(false);
    page.setVisible(true);
    await flush();
    attention.dismissResume();
    page.setVisible(true);
    expect(state().resume).toBeUndefined();
  });

  it('is removed by dismissResume', async () => {
    const { page, attention, state } = setup({ runs: [run('a', 1, 'running')] });
    page.setVisible(false);
    page.setVisible(true);
    await flush();
    expect(state().resume).toBeDefined();
    attention.dismissResume();
    expect(state().resume).toBeUndefined();
    expect('resume' in state()).toBe(false);
  });

  it('does not bring the report back when a probe answers after dismissResume', async () => {
    let answer: () => void = () => undefined;
    const { page, attention, state } = setup({
      runs: [run('a', 1, 'running')],
      probe: () => new Promise<void>((resolve) => (answer = resolve)),
    });
    page.setVisible(false);
    page.setVisible(true);
    attention.dismissResume();
    answer();
    await flush();
    expect(state().resume).toBeUndefined();
  });

  it('lets the latest show decide the Data worker status', async () => {
    const answers: Array<() => void> = [];
    const { page, clock, state } = setup({
      runs: [run('a', 1, 'running')],
      probe: () => new Promise<void>((resolve) => answers.push(resolve)),
    });
    page.setVisible(false);
    page.setVisible(true);
    clock.now += 10;
    page.setVisible(false);
    page.setVisible(true);
    expect(answers).toHaveLength(2);
    answers[1]!();
    await flush();
    expect(state().resume?.dataWorker).toBe('responding');
    answers[0]!();
    await flush();
    expect(state().resume?.dataWorker).toBe('responding');
  });

  it('makes no report for a hide shorter than minHiddenMs', () => {
    const { page, clock, state } = setup({ runs: [run('a', 1, 'running')], minHiddenMs: 1000 });
    page.setVisible(false);
    clock.now += 999;
    page.setVisible(true);
    expect(state().resume).toBeUndefined();
    page.setVisible(false);
    clock.now += 1000;
    page.setVisible(true);
    expect(state().resume?.hiddenMs).toBe(1000);
  });

  it('lets a late answer of the Data worker correct "not-responding"', async () => {
    let answer: (() => void) | undefined;
    const { page, state } = setup({
      runs: [run('a', 1, 'running')],
      probe: () => new Promise<void>((resolve) => (answer = resolve)),
      probeTimeoutMs: 20,
    });
    page.setVisible(false);
    page.setVisible(true);
    await sleep(40);
    expect(state().resume?.dataWorker).toBe('not-responding');
    answer!();
    await flush();
    expect(state().resume?.dataWorker).toBe('responding');
  });
});
