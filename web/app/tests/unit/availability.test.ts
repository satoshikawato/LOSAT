// What the page says when a part of LOSAT could not be loaded or stopped (line WQA's sweep, D1 and D2, 2026-11-19):
// an engine file that could not be downloaded made every input card say "This input cannot be read: Failed to fetch",
// refused Run LOSAT for the pasted query, and said nothing for the page. The search form tells such a failure from an
// unreadable input, and the page says it once with the remedy (application/availability.ts).
import { describe, expect, it } from 'vitest';
import { sha256Hex, UnavailableError } from 'losat';
import type { ProgramDescription, StartResult } from 'losat/ports';
import { Availability, failureText, unavailableOf, withRemedy } from '../../src/application/availability';
import type { EnqueueAllResult } from '../../src/application/coordinator';
import { SearchDraft, type DraftSource } from '../../src/application/draft';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeAccessionFetcher } from '../../src/infra/fake/fake-fetcher';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import { NodeInflater } from './support/node-inflater';

const SERIAL = 'LOSAT could not download its engine file /engine/losat-web-serial-0123.wasm (network error)';

/**
 * A search form whose Data worker cannot load the engine's serial module: every call that needs it fails so, while its
 * flag in `fail` is set (a test clears it when the module downloads).
 */
function setup(fail: { describe?: boolean; index?: boolean; check?: boolean; validate?: boolean } = {}) {
  let token = 0;
  const data = new DataService({
    store: new MemoryBlockStore(),
    scanner: new FakeScanner(),
    inflater: new NodeInflater(),
    checker: new FakeInputChecker(),
    digest: sha256Hex,
    newToken: () => `token-${++token}`,
    cleanup: Promise.resolve({ state: 'done', removedSessions: 0 }),
  });
  const unavailable = () => Promise.reject(new UnavailableError(SERIAL));
  const calls = { index: 0 };
  const indexSource = data.indexSource.bind(data);
  data.indexSource = (...args) => {
    calls.index++;
    return fail.index === true ? unavailable() : indexSource(...args);
  };
  const checkInput = data.checkInput.bind(data);
  data.checkInput = (...args) => (fail.check === true ? unavailable() : checkInput(...args));
  const fake = new FakeEngine();
  const engine = {
    describe: (program: ProgramDescription['program']) => (fail.describe === true ? unavailable() : fake.describe(program)),
    start: async (): Promise<StartResult> => {
      if (fail.validate === true) throw new UnavailableError(SERIAL);
      return { ok: true, start: { conffile: null, logFile: null, dryRun: false } };
    },
    validate: async () => ({ ok: true }) as const,
    prepare: () => undefined,
  };
  const draft = new SearchDraft({
    engine,
    data,
    enqueueAll: async (batch): Promise<EnqueueAllResult> => ({ ok: true, runIds: batch.map((_, i) => `r${i}`) }),
    fetcher: new FakeAccessionFetcher(),
    now: () => 0,
    debounceMs: 0,
  });
  return { draft, calls };
}

const sourceOf = (draft: SearchDraft, role: 'query' | 'subject'): DraftSource => draft.state.get()[role].sources[0]!;

describe('a part of LOSAT that is unavailable, in the search form', () => {
  it('an input that could not be read for it is not the input\'s fault, and Run LOSAT says what the page says', async () => {
    const { draft } = setup({ index: true });
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', 'ACGT\n');
    await draft.idle();
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'failed', error: SERIAL, unavailable: true });
    expect(unavailableOf(draft.state.get())).toBe(SERIAL);
    const result = await draft.submit();
    expect(result).toEqual({ ok: false, message: `${SERIAL}. Reload the page.` });
    expect(draft.readiness()).toBe(`${SERIAL}. Reload the page.`);
  });

  it('an input that the engine cannot read is still the input\'s', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>q\nACGT\n>?10\n');
    await draft.idle();
    expect(sourceOf(draft, 'query').status).toBe('failed');
    expect('unavailable' in sourceOf(draft, 'query')).toBe(false);
    expect(unavailableOf(draft.state.get())).toBeUndefined();
  });

  it('a describe, a check and a validation that failed for it are no verdicts, and the form finds the part unavailable', async () => {
    const described = setup({ describe: true });
    await described.draft.idle();
    expect(described.draft.state.get()).toMatchObject({ descriptionError: SERIAL, descriptionUnavailable: true });
    expect(unavailableOf(described.draft.state.get())).toBe(SERIAL);

    const checked = setup({ check: true });
    checked.draft.setPaste('query', '>q\nACGT\n');
    await checked.draft.idle();
    expect(sourceOf(checked.draft, 'query').check).toEqual({ state: 'error', message: SERIAL, unavailable: true });
    expect(unavailableOf(checked.draft.state.get())).toBe(SERIAL);

    const validated = setup({ validate: true });
    validated.draft.setPaste('query', '>q\nACGT\n');
    validated.draft.setPaste('subject', '>s\nACGT\n');
    await validated.draft.idle();
    // Not a refusal of the options, which the form would show beside Run LOSAT.
    expect(validated.draft.state.get().validation).toEqual({ state: 'unavailable', message: SERIAL });
    expect(unavailableOf(validated.draft.state.get())).toBe(SERIAL);
  });

  it('the page says the form\'s failure, then the Data worker\'s for good', async () => {
    const { draft } = setup({ index: true });
    let end: (error: Error) => void = () => undefined;
    const availability = new Availability({ draft, dataWorkerEnded: new Promise((resolve) => (end = resolve)) });
    expect(availability.state.get()).toEqual({});
    draft.setPaste('query', '>q\nACGT\n');
    await draft.idle();
    expect(availability.state.get()).toEqual({ problem: SERIAL });
    draft.removeSource('query', 'paste');
    expect(availability.state.get()).toEqual({});
    end(new UnavailableError('The part of LOSAT that keeps the results stopped'));
    await Promise.resolve();
    expect(availability.state.get()).toEqual({ problem: 'The part of LOSAT that keeps the results stopped' });
    draft.setPaste('query', '>q\nACGT\n');
    await draft.idle();
    expect(availability.state.get()).toEqual({ problem: 'The part of LOSAT that keeps the results stopped' });
  });
});

// The second review of line WQA's D1 and D2, M1 (2026-11-19): an engine file that failed to download once is downloaded
// again by the next call, but the options were asked only at load and when the mode changed. One network error at load
// left "Not read." on the options and "Reload the page; the results of this session are then lost." on the page for the
// whole session, also after a search completed, and the Examples were refused.
describe('a part of LOSAT that could not be loaded and then can', () => {
  const never = new Promise<Error>(() => undefined);

  it('the options are read again once the engine answers, and the page no longer says that LOSAT is unavailable', async () => {
    const fail = { describe: true, validate: true };
    const { draft } = setup(fail);
    const availability = new Availability({ draft, dataWorkerEnded: never });
    await draft.idle();
    expect(availability.state.get()).toEqual({ problem: SERIAL });
    fail.describe = false;
    fail.validate = false;
    draft.setPaste('query', '>q\nACGT\n');
    await draft.idle();
    expect(draft.state.get().description?.program).toBe('blastn');
    expect(draft.state.get()).not.toHaveProperty('descriptionError');
    expect(draft.state.get()).not.toHaveProperty('descriptionUnavailable');
    expect(availability.state.get()).toEqual({});
  });

  it('a validation that the engine answers is enough', async () => {
    const fail = { describe: true };
    const { draft } = setup(fail);
    await draft.idle();
    expect(unavailableOf(draft.state.get())).toBe(SERIAL);
    fail.describe = false;
    draft.setExtra('-qcov_hsp_perc 50');
    await draft.idle();
    expect(draft.state.get().validation.state).toBe('ok');
    expect(draft.state.get().description?.program).toBe('blastn');
    expect(unavailableOf(draft.state.get())).toBeUndefined();
  });

  it('inputs that could not be read or checked for it are read and checked again', async () => {
    const fail = { index: true, check: false };
    const { draft } = setup(fail);
    draft.setPaste('query', '>q\nACGT\n');
    await draft.idle();
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'failed', unavailable: true });
    fail.index = false;
    fail.check = true;
    draft.setPaste('subject', '>s\nACGT\n');
    await draft.idle();
    expect(sourceOf(draft, 'query').status).toBe('ready');
    expect(sourceOf(draft, 'subject').check).toEqual({ state: 'error', message: SERIAL, unavailable: true });
    fail.check = false;
    draft.setPaste('query', '>q2\nACGT\n');
    await draft.idle();
    expect(sourceOf(draft, 'subject').check?.state).toBe('ok');
    expect(unavailableOf(draft.state.get())).toBeUndefined();
  });

  // The third review of D1 and D2, L1 (2026-11-19): no test failed when the answers to the options and to a check no
  // longer asked again; the form then left the inputs "Not read." and Run LOSAT refused.
  it('an answer to the options reads again the inputs that could not be read for it', async () => {
    const fail = { describe: true, validate: true, index: true };
    const { draft } = setup(fail);
    draft.setPaste('query', '>q\nACGT\n');
    await draft.idle();
    fail.describe = false;
    fail.validate = false;
    fail.index = false;
    // A validation asks the options again (askOptionsAgain); their answer asks the rest.
    draft.setExtra('-qcov_hsp_perc 50');
    await draft.idle();
    expect(sourceOf(draft, 'query').status).toBe('ready');
    expect(unavailableOf(draft.state.get())).toBeUndefined();
  });

  it('an answer to a check reads again the inputs that could not be read for it', async () => {
    const fail = { check: true, index: false };
    const { draft } = setup(fail);
    draft.setPaste('query', '>q1\nACGT\n>q2\nACGT\n');
    await draft.idle();
    expect(sourceOf(draft, 'query').check).toEqual({ state: 'error', message: SERIAL, unavailable: true });
    fail.index = true;
    draft.setPaste('subject', '>s\nACGT\n');
    await draft.idle();
    expect(sourceOf(draft, 'subject')).toMatchObject({ status: 'failed', unavailable: true });
    fail.check = false;
    fail.index = false;
    // Only a check asks the engine: a record of the query is left out.
    draft.setIncluded('query', 'paste', [1], false);
    await draft.idle();
    expect(sourceOf(draft, 'query').check?.state).toBe('ok');
    expect(sourceOf(draft, 'subject').status).toBe('ready');
    expect(unavailableOf(draft.state.get())).toBeUndefined();
  });

  it('an input that still cannot be read for it is not read again after each validation', async () => {
    const { draft, calls } = setup({ index: true });
    draft.setPaste('query', '>q\nACGT\n');
    await draft.idle();
    const reads = calls.index;
    draft.setExtra('-qcov_hsp_perc 50');
    await draft.idle();
    draft.setExtra('-qcov_hsp_perc 60');
    await draft.idle();
    expect(calls.index).toBe(reads);
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'failed', unavailable: true });
  });
});

// The third review of D1 and D2, M1 (2026-11-19): Run LOSAT and the Examples decided on what the form had found, and
// asked nothing: with the inputs pasted while the engine file could not be downloaded, they were refused with "Reload
// the page" also once it could, though one more request would have read them.
describe('Run LOSAT and the Examples ask the engine again before they refuse', () => {
  const never = new Promise<Error>(() => undefined);

  /** Both inputs pasted while the engine file could not be downloaded. */
  async function pastedWhileUnavailable() {
    const fail = { describe: true, validate: true, index: true };
    const { draft, calls } = setup(fail);
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s\nACGT\n');
    await draft.idle();
    expect(unavailableOf(draft.state.get())).toBe(SERIAL);
    return { fail, draft, calls, reads: calls.index };
  }

  const back = (fail: { describe: boolean; validate: boolean; index: boolean }) => {
    fail.describe = false;
    fail.validate = false;
    fail.index = false;
  };

  it('Run LOSAT reads the inputs again and starts the search once the engine file downloads', async () => {
    const { fail, draft, calls, reads } = await pastedWhileUnavailable();
    back(fail);
    expect(await draft.submit()).toEqual({ ok: true, runIds: ['r0'] });
    expect(calls.index).toBe(reads + 2);
    expect(unavailableOf(draft.state.get())).toBeUndefined();
  });

  it('Run LOSAT says what the page says when the engine still cannot answer, having asked once', async () => {
    const { draft, calls, reads } = await pastedWhileUnavailable();
    expect(await draft.submit()).toEqual({ ok: false, message: `${SERIAL}. Reload the page.` });
    expect(calls.index).toBe(reads + 2);
  });

  it('the problem that an Example decides on is the one after the form asked again', async () => {
    const { fail, draft, calls, reads } = await pastedWhileUnavailable();
    const availability = new Availability({ draft, dataWorkerEnded: never });
    expect(await availability.recheck()).toBe(SERIAL);
    expect(calls.index).toBe(reads + 2);
    back(fail);
    expect(await availability.recheck()).toBeUndefined();
    expect(sourceOf(draft, 'query').status).toBe('ready');
    expect(draft.state.get().description?.program).toBe('blastn');
  });

  it('a Data worker that ended is not asked again', async () => {
    const STOPPED = 'The part of LOSAT that keeps the results stopped';
    const { draft, calls } = setup({ index: true });
    let end: (error: Error) => void = () => undefined;
    const availability = new Availability({ draft, dataWorkerEnded: new Promise((resolve) => (end = resolve)) });
    draft.setPaste('query', '>q\nACGT\n');
    await draft.idle();
    end(new UnavailableError(STOPPED));
    await Promise.resolve();
    const reads = calls.index;
    expect(await availability.recheck()).toBe(STOPPED);
    expect(calls.index).toBe(reads);
  });
});

describe('the remedy', () => {
  it('ends the message as a sentence, and says what a reload loses', () => {
    expect(withRemedy(SERIAL)).toBe(`${SERIAL}. Reload the page.`);
    expect(withRemedy('The part of LOSAT that keeps the results stopped.', true)).toBe(
      'The part of LOSAT that keeps the results stopped. Reload the page; the results of this session are then lost.',
    );
  });

  it('is given only for a part of LOSAT that is unavailable, also after the error crossed a worker', () => {
    const crossed = new Error(SERIAL);
    crossed.name = 'UnavailableError';
    expect(failureText(crossed)).toBe(`${SERIAL}. Reload the page.`);
    expect(failureText(new Error('BLAST engine error: x'))).toBe('BLAST engine error: x');
  });
});
