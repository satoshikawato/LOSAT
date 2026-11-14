// The Outputs' reads and exports (src/application/outputs.ts): an output that cannot be read shows
// no text and says why, as the Results do; an export that is not saved says why in the same place;
// a read or an export that ends after a later choice, a later export or a reset changes nothing (WA2-BR-01).
import { describe, expect, it } from 'vitest';
import type { RunOutput } from 'losat/domain';
import type { RunView } from '../../src/application/coordinator';
import { OutputReader, OutputSaver, saveOutput, type OutputPort } from '../../src/application/outputs';
import { RESULTS_REMOVED_MESSAGE, StorageRemovedError } from '../../src/infra/data/block-store';

function run(number: number): RunView {
  const input = { name: 'x.fa', bytes: new Uint8Array(), sha256: '0', revisionIds: [], records: [{ id: 'q', length: 4 }] };
  return {
    snapshot: { runId: `r${number}`, number, program: 'blastn', argv: ['blastn', '-query', 'q.fa', '-subject', 's.fa'], query: input, subject: input, requestedThreads: 1, queuedAt: 0 },
    status: 'completed',
    record: { runtimePath: 'serial', threads: 1 },
  };
}

interface Pending<T> {
  readonly runId: string;
  readonly format: RunOutput;
  resolve(value: T): void;
  reject(error: unknown): void;
}

/** Reads and exports that end when the test says so. */
class Outputs implements OutputPort {
  readonly reads: Pending<string>[] = [];
  readonly checks: Pending<void>[] = [];
  readonly exports: Pending<void>[] = [];

  readOutput(runId: string, format: RunOutput): Promise<string> {
    return new Promise((resolve, reject) => this.reads.push({ runId, format, resolve, reject }));
  }

  checkOutput(runId: string, format: RunOutput): Promise<void> {
    return new Promise((resolve, reject) => this.checks.push({ runId, format, resolve, reject }));
  }

  exportOutput(runId: string, format: RunOutput): Promise<void> {
    return new Promise((resolve, reject) => this.exports.push({ runId, format, resolve, reject }));
  }

  /** The saves of a run's files (RunRecord.runFiles), by run and index. */
  readonly fileExports: Array<{ readonly runId: string; readonly index: number; resolve(): void; reject(error: unknown): void }> = [];

  exportRunFile(runId: string, index: number): Promise<void> {
    return new Promise((resolve, reject) => this.fileExports.push({ runId, index, resolve, reject }));
  }
}

const settle = () => new Promise((resolve) => setTimeout(resolve, 0));
const removed = () => new StorageRemovedError(RESULTS_REMOVED_MESSAGE);
const UNREADABLE = `The results of run 1 could not be read: ${RESULTS_REMOVED_MESSAGE}`;
const NOT_SAVED = `Not saved: ${RESULTS_REMOVED_MESSAGE}`;

function setup() {
  const outputs = new Outputs();
  const reader = new OutputReader(outputs);
  return { outputs, reader, state: () => reader.state.get() };
}

describe('OutputReader', () => {
  it('shows the output read, and keeps the text shown while the next one is read', async () => {
    const { outputs, reader, state } = setup();
    void reader.show(run(1), 6);
    expect(state()).toEqual({ shown: '', text: '', failed: '', problem: '' });
    outputs.reads[0]!.resolve('q\ts\n');
    await settle();
    expect(state()).toEqual({ shown: '1:6', text: 'q\ts\n', failed: '', problem: '' });
    void reader.show(run(1), 0);
    expect(outputs.reads[1]).toMatchObject({ runId: 'r1', format: 0 });
    expect(state()).toEqual({ shown: '', text: 'q\ts\n', failed: '', problem: '' });
    outputs.reads[1]!.resolve('BLASTN\n');
    await settle();
    expect(state()).toEqual({ shown: '1:0', text: 'BLASTN\n', failed: '', problem: '' });
  });

  it('an output that cannot be read shows no text, not the format shown before, and says why as the Results do', async () => {
    const { outputs, reader, state } = setup();
    void reader.show(run(1), 6);
    outputs.reads[0]!.resolve('q\ts\n');
    await settle();
    void reader.show(run(1), 0);
    outputs.reads[1]!.reject(removed());
    await settle();
    expect(state()).toEqual({ shown: '', text: '', failed: '1:0', problem: UNREADABLE });
    // The next format is read; the reason stays until it is, and `failed` names the read that failed.
    void reader.show(run(1), 6);
    expect(state()).toEqual({ shown: '', text: '', failed: '1:0', problem: UNREADABLE });
    outputs.reads[2]!.reject(removed());
    await settle();
    expect(state()).toEqual({ shown: '', text: '', failed: '1:6', problem: UNREADABLE });
    void reader.show(run(1), 7);
    outputs.reads[3]!.resolve('# BLASTN\n');
    await settle();
    expect(state()).toEqual({ shown: '1:7', text: '# BLASTN\n', failed: '', problem: '' });
  });

  // W11: a binary output (outfmt 9) is opened to know that it can be read, not read, and shows no text.
  it('a binary output shows no text, and says why when it cannot be read', async () => {
    const { outputs, reader, state } = setup();
    void reader.show(run(1), 4, false);
    expect(outputs.checks[0]).toMatchObject({ runId: 'r1', format: 4 });
    expect(outputs.reads).toHaveLength(0);
    outputs.checks[0]!.resolve();
    await settle();
    expect(state()).toEqual({ shown: '1:4', text: '', failed: '', problem: '' });
    void reader.show(run(1), 4, false);
    outputs.checks[1]!.reject(removed());
    await settle();
    expect(state()).toEqual({ shown: '', text: '', failed: '1:4', problem: UNREADABLE });
  });

  it('a listener that throws on a read output is not told as a read that failed', async () => {
    const { outputs, reader, state } = setup();
    const unsubscribe = reader.state.subscribe((value) => {
      if (value.shown !== '') throw new Error('listener failed');
    });
    const shown = reader.show(run(1), 6);
    outputs.reads[0]!.resolve('q\ts\n');
    await expect(shown).rejects.toThrow('listener failed');
    unsubscribe();
    expect(state()).toEqual({ shown: '1:6', text: 'q\ts\n', failed: '', problem: '' });
  });

  it('a read that ends after a later choice changes nothing, whether it is read or not', async () => {
    const { outputs, reader, state } = setup();
    void reader.show(run(1), 0);
    void reader.show(run(1), 7);
    outputs.reads[1]!.resolve('seven\n');
    await settle();
    outputs.reads[0]!.resolve('zero\n');
    await settle();
    expect(state()).toEqual({ shown: '1:7', text: 'seven\n', failed: '', problem: '' });

    void reader.show(run(1), 0);
    void reader.show(run(1), 6);
    outputs.reads[3]!.resolve('six\n');
    await settle();
    outputs.reads[2]!.reject(removed());
    await settle();
    expect(state()).toEqual({ shown: '1:6', text: 'six\n', failed: '', problem: '' });

    // A failure that ends first, while the later read is still being read, is not shown either.
    void reader.show(run(1), 0);
    void reader.show(run(2), 6);
    outputs.reads[4]!.reject(removed());
    await settle();
    expect(state()).toEqual({ shown: '', text: 'six\n', failed: '', problem: '' });
    outputs.reads[5]!.resolve('run 2\n');
    await settle();
    expect(state()).toEqual({ shown: '2:6', text: 'run 2\n', failed: '', problem: '' });
  });

  it('an export that is not saved says why where the output is; the next export or another format clears it', async () => {
    const { outputs, reader, state } = setup();
    void reader.show(run(1), 6);
    outputs.reads[0]!.resolve('q\ts\n');
    await settle();
    void reader.export(run(1), 6);
    expect(outputs.exports[0]).toMatchObject({ runId: 'r1', format: 6 });
    outputs.exports[0]!.reject(removed());
    await settle();
    // The text read before stays: it is the output of the format chosen.
    expect(state()).toEqual({ shown: '1:6', text: 'q\ts\n', failed: '', problem: NOT_SAVED });
    // The reason stays while the next export runs, and goes when it is saved.
    void reader.export(run(1), 6);
    expect(state().problem).toBe(NOT_SAVED);
    outputs.exports[1]!.resolve();
    await settle();
    expect(state().problem).toBe('');

    void reader.export(run(1), 6);
    outputs.exports[2]!.reject(new Error('the download was refused'));
    await settle();
    expect(state().problem).toBe('Not saved: the download was refused');
    void reader.show(run(1), 7);
    expect(state().problem).toBe('');
  });

  it('an export that ends after another format was chosen, or after a later export, says nothing', async () => {
    const { outputs, reader, state } = setup();
    void reader.show(run(1), 6);
    outputs.reads[0]!.resolve('q\ts\n');
    await settle();
    void reader.export(run(1), 6);
    void reader.show(run(1), 7);
    outputs.reads[1]!.resolve('# BLASTN\n');
    outputs.exports[0]!.reject(removed());
    await settle();
    expect(state()).toEqual({ shown: '1:7', text: '# BLASTN\n', failed: '', problem: '' });

    void reader.export(run(1), 7);
    void reader.export(run(1), 7);
    outputs.exports[2]!.resolve();
    outputs.exports[1]!.reject(removed());
    await settle();
    expect(state().problem).toBe('');
  });

  it('after a read failed, an export that is not saved says so, and one that is saved leaves why there is no text', async () => {
    const { outputs, reader, state } = setup();
    void reader.show(run(1), 0);
    outputs.reads[0]!.reject(removed());
    await settle();
    void reader.export(run(1), 0);
    outputs.exports[0]!.reject(removed());
    await settle();
    expect(state()).toEqual({ shown: '', text: '', failed: '1:0', problem: NOT_SAVED });
    void reader.export(run(1), 0);
    outputs.exports[1]!.resolve();
    await settle();
    expect(state()).toEqual({ shown: '', text: '', failed: '1:0', problem: UNREADABLE });
  });

  it('exports that keep failing after a read failed tell their reason once, not the read and the export in turn', async () => {
    const { outputs, reader, state } = setup();
    const told: string[] = [];
    reader.state.subscribe(({ problem }) => {
      if (told.at(-1) !== problem) told.push(problem);
    });
    void reader.show(run(1), 0);
    outputs.reads[0]!.reject(removed());
    await settle();
    for (let i = 0; i < 3; i++) {
      void reader.export(run(1), 0);
      outputs.exports[i]!.reject(removed());
      await settle();
    }
    expect(told).toEqual(['', UNREADABLE, NOT_SAVED]);
    expect(state().problem).toBe(NOT_SAVED);
  });
});

describe('OutputSaver', () => {
  // The file of a run's -logfile (line RG E1w-3b) is saved as an output is, and says why it is not.
  it("saves a run's file, says why one is not saved, and the latest save wins", async () => {
    const outputs = new Outputs();
    const saver = new OutputSaver(outputs);
    void saver.saveRunFile(run(1), 0);
    expect(outputs.fileExports[0]).toMatchObject({ runId: 'r1', index: 0 });
    outputs.fileExports[0]!.reject(new Error('the run has no such file'));
    await settle();
    expect(saver.problem.get()).toBe('Not saved: the run has no such file');
    void saver.saveRunFile(run(1), 0);
    void saver.save(run(1), 6);
    outputs.exports[0]!.resolve();
    outputs.fileExports[1]!.reject(removed());
    await settle();
    expect(saver.problem.get()).toBe('');
  });

  it('a save that ends after a reset says nothing: another run or format was chosen', async () => {
    const outputs = new Outputs();
    const saver = new OutputSaver(outputs);
    void saver.save(run(1), 6);
    expect(outputs.exports[0]).toMatchObject({ runId: 'r1', format: 6 });
    saver.reset();
    outputs.exports[0]!.reject(removed());
    await settle();
    expect(saver.problem.get()).toBe('');
    // A reset forgets the reason told.
    void saver.save(run(1), 6);
    outputs.exports[1]!.reject(removed());
    await settle();
    expect(saver.problem.get()).toBe(NOT_SAVED);
    saver.reset();
    expect(saver.problem.get()).toBe('');
  });

  it('the latest save wins over one that ends after it; the reason stays while the next save runs', async () => {
    const outputs = new Outputs();
    const saver = new OutputSaver(outputs);
    void saver.save(run(1), 0);
    void saver.save(run(1), 7);
    outputs.exports[1]!.resolve();
    outputs.exports[0]!.reject(removed());
    await settle();
    expect(saver.problem.get()).toBe('');
    void saver.save(run(1), 0);
    void saver.save(run(1), 7);
    outputs.exports[3]!.reject(removed());
    outputs.exports[2]!.resolve();
    await settle();
    expect(saver.problem.get()).toBe(NOT_SAVED);
    void saver.save(run(1), 6);
    expect(saver.problem.get()).toBe(NOT_SAVED);
    outputs.exports[4]!.resolve();
    await settle();
    expect(saver.problem.get()).toBe('');
  });
});

describe('saveOutput', () => {
  it("resolves to '' when the output is saved, and to why when it is not, never rejecting", async () => {
    const outputs = new Outputs();
    const saved = saveOutput(outputs, run(3), 7);
    expect(outputs.exports[0]).toMatchObject({ runId: 'r3', format: 7 });
    outputs.exports[0]!.resolve();
    expect(await saved).toBe('');
    const failed = saveOutput(outputs, run(3), 0);
    outputs.exports[1]!.reject(removed());
    expect(await failed).toBe(NOT_SAVED);
    const odd = saveOutput(outputs, run(3), 0);
    outputs.exports[2]!.reject('refused');
    expect(await odd).toBe('Not saved: refused');
  });
});
