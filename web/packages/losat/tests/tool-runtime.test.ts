// The Engine worker's side of a tool run on the built serial reactor, in Node (review L8, M3 of WW3a, 2026-11-23): the
// `EngineRuntime.tool` that the worker calls (the run's directory opened and let go, the retained input released, a second run
// on the same instance), the worker's check of a tool's words with the run's directory (`EngineRuntime.validateTool`, which
// reads `-conffile`), and the command line's answer for `-in <folder>` (the open of a folder succeeds, its read fails with
// EISDIR). Runs with LOSAT_WEB_REACTORS (the folder of web/adapter/tools/build_reactors.py's output); the windowmasker rows read
// no frame of stream 12, so the reactors of an earlier frame layout serve as well. The makeblastdb rows (ABI version 17, line
// NP WW4b1, 2026-11-27) need reactors of version 17: makeblastdb writes its database into the run's writable directory, and the
// files that the run left there come after the run as frames, by name; a database given in the directory is replaced as the
// command line replaces it, and the directory goes back to nothing after the run. The directory reads a file's bytes through
// the `FileReaderSync` that a worker has, which Node lacks: the test gives it Blobs that carry their bytes.
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { afterAll, beforeAll, describe, expect, it } from 'vitest';
import type { DirectoryFile } from '../src/domain/run-files.ts';
import type { ToolFrame } from '../src/domain/tool-frames.ts';
import { ROLE_SUBJECT, type ReactorAbi } from '../src/infra/reactor/abi.ts';
import { compileReactor } from '../src/infra/reactor/instance.ts';
import type { MainLinks } from '../src/infra/engine-worker/link.ts';
import type { WorkerTool, WorkerToolValidate } from '../src/infra/engine-worker/protocol.ts';
import { InstanceUse } from '../src/infra/engine-worker/instance-use.ts';
import { EngineRuntime, type ToolPosts } from '../src/infra/engine-worker/runtime.ts';

const LOSAT = fileURLToPath(new URL('../../../../LOSAT', import.meta.url));
const REACTORS = process.env['LOSAT_WEB_REACTORS'];

/** A Blob that carries its bytes, so that a synchronous read can give them (a slice of it carries its part). */
class MemoryBlob extends Blob {
  constructor(readonly data: Uint8Array) {
    super([data as BlobPart]);
  }

  override slice(start?: number, end?: number): Blob {
    return new MemoryBlob(this.data.slice(start, end));
  }
}

const memoryFile = (name: string, contents: Uint8Array | string): DirectoryFile => ({
  name,
  file: new MemoryBlob(typeof contents === 'string' ? new TextEncoder().encode(contents) : contents),
});

const fixture = (path: string) => new Uint8Array(readFileSync(join(LOSAT, 'tests/winmask', path)));
const corpus = (path: string) => new Uint8Array(readFileSync(join(LOSAT, 'tests/dropin_corpus', path)));
const NO_LINKS: MainLinks = {
  link: () => {
    throw new Error('unused');
  },
  fault: () => undefined,
};

/**
 * What a run posted. The bytes are taken as the Engine worker's posts take them (tool-posts.ts): their buffer is transferred
 * (here to a copy that the test keeps), which leaves the runtime's array empty, as the main thread's receipt of a post does.
 */
function posts(): ToolPosts & { stdout: Uint8Array[]; stderr: Uint8Array[]; frames: ToolFrame[] } {
  const posted = { stdout: [] as Uint8Array[], stderr: [] as Uint8Array[], frames: [] as ToolFrame[] };
  const taken = (bytes: Uint8Array) => structuredClone(bytes, { transfer: [bytes.buffer as ArrayBuffer] });
  return {
    ...posted,
    chunk: (stream, bytes) => posted[stream].push(taken(bytes)),
    file: (frame) => posted.frames.push(frame.kind === 'bytes' ? { ...frame, bytes: taken(frame.bytes) } : frame),
  };
}
const joined = (chunks: Uint8Array[]) => Buffer.concat(chunks);

/**
 * The files of a run's frames, in the order of their begin frames: each file's bytes as its frames wrote them, and whether
 * its frames came in the form that the overlay's files come in (begin, bytes from offset 0 one after another, end).
 */
function filesOf(frames: readonly ToolFrame[]): Array<{ name: string; bytes: Buffer; inOrder: boolean }> {
  const files = new Map<number, { name: string; chunks: Uint8Array[]; length: number; inOrder: boolean; ended: boolean }>();
  for (const frame of frames) {
    if (frame.kind === 'begin') {
      files.set(frame.number, { name: frame.name, chunks: [], length: 0, inOrder: true, ended: false });
      continue;
    }
    const file = files.get(frame.number)!;
    if (frame.kind === 'bytes') {
      if (frame.offset !== file.length || file.ended) file.inOrder = false;
      file.chunks.push(frame.bytes);
      file.length += frame.bytes.length;
    } else if (frame.kind === 'end') file.ended = true;
    else file.inOrder = false;
  }
  return [...files.values()].map((file) => ({ name: file.name, bytes: joined(file.chunks), inOrder: file.inOrder && file.ended }));
}

describe.skipIf(REACTORS === undefined)('the Engine worker runtime and tools, on the built serial reactor', () => {
  let module: WebAssembly.Module;
  beforeAll(async () => {
    // The directory of a worker reads Blobs with FileReaderSync (run-directory.ts, `fileReaderSync`).
    (globalThis as { FileReaderSync?: unknown }).FileReaderSync = class {
      readAsArrayBuffer(blob: MemoryBlob): ArrayBuffer {
        return blob.data.slice().buffer;
      }
    };
    module = (await compileReactor(readFileSync(join(REACTORS!, 'losat-web-serial.wasm')))).module;
  });
  afterAll(() => {
    delete (globalThis as { FileReaderSync?: unknown }).FileReaderSync;
  });

  const runtime = () =>
    new EngineRuntime({
      serial: async () => module,
      threads: () => Promise.reject(new Error('unused')),
      builds: { serial: 'serial', threads: 'threads' },
      main: NO_LINKS,
      use: undefined,
    });
  const toolOf = (argv: string[], files?: DirectoryFile[]): WorkerTool => ({ type: 'tool', generation: 1, runId: 't', argv, ...(files === undefined ? {} : { files }) });
  const validateOf = (argv: string[], files?: DirectoryFile[]): WorkerToolValidate => ({ type: 'tool-validate', generation: 1, id: 'v', argv, ...(files === undefined ? {} : { files }) });

  it('opens the run directory for a tool run, lets it go after, and runs again on the same instance', async () => {
    const worker = runtime();
    const genome = memoryFile('g.fa', fixture('fasta/dm2L_3k.fa'));
    const argv = ['windowmasker', '-mk_counts', '-in', 'g.fa', '-unit', '5', '-out', '-'];
    const first = posts();
    const one = await worker.tool(toolOf(argv, [genome]), first);
    expect(one.status).toBe(0);
    const counts = joined(first.stdout.length === 0 ? [] : first.stdout);
    expect(counts.length).toBeGreaterThan(0);
    expect(one.memory.instanceRuns).toBe(1);
    // The run's files were let go with the run: the same words without them fail as the command line does in an empty folder.
    const empty = posts();
    const two = await worker.tool(toolOf(argv), empty);
    expect(two.status).not.toBe(0);
    expect(joined(empty.stderr).length).toBeGreaterThan(0);
    expect(joined(empty.stdout).length).toBe(0);
    // A third run on the same instance with the files again: the same bytes.
    const third = posts();
    const again = await worker.tool(toolOf(argv, [genome]), third);
    expect(again.status).toBe(0);
    expect(again.memory.instanceRuns).toBe(3);
    expect(joined(third.stdout).equals(counts)).toBe(true);
  });

  it('releases the subject that the instance retained for searches before it runs a tool', async () => {
    const worker = runtime();
    await worker.tool(toolOf(['windowmasker', '-version']), posts());
    const loaded = (worker as unknown as { loaded: { retained?: unknown; instance: { abi: ReactorAbi } } }).loaded;
    const registered = loaded.instance.abi.register('blastp', ROLE_SUBJECT, false, false, new TextEncoder().encode('>s\nMKVLAGG\n'));
    loaded.retained = { key: 'k', dataLoaders: false, handle: registered.handle, records: registered.records };
    await worker.tool(toolOf(['windowmasker', '-version']), posts());
    expect(loaded.retained).toBeUndefined();
    // The handle is the engine's no more.
    expect(() => loaded.instance.abi.release(registered.handle)).toThrow();
  });

  it("checks a tool's words with the run's directory: `-conffile` names a file that the command line reads (M3)", async () => {
    const worker = runtime();
    const argv = ['windowmasker', '-conffile', 'wm.ini', '-mk_counts', '-in', 'g.fa', '-out', 'c'];
    const files = [memoryFile('wm.ini', '[BLAST]\nBLAST_USAGE_REPORT=false\n'), memoryFile('g.fa', fixture('fasta/dm2L_3k.fa'))];
    expect(await worker.validateTool(validateOf(argv, files))).toEqual({ ok: true });
    // No instance was loaded: the check loaded the serial one and marks it, so that a preparation for threads may replace it.
    expect((worker as unknown as { loaded: { path: string; checkOnly: boolean } }).loaded).toMatchObject({ path: 'serial', checkOnly: true });
    // Without the directory the registry file is not there, as for the command line in a folder without it.
    const refused = await worker.validateTool(validateOf(argv));
    expect(refused.ok).toBe(false);
    expect(!refused.ok && refused.message).toContain('Registry file "wm.ini" not found');
    // The files were let go after the first check, and a run of the same words takes its own.
    const refusedAgain = await worker.validateTool(validateOf(argv));
    expect(refusedAgain.ok).toBe(false);
    const run = posts();
    // (To stdout, as the reactors of an earlier frame layout also serve this test.)
    expect((await worker.tool(toolOf([...argv.slice(0, -1), '-'], files), run)).status).toBe(0);
    // A malformed registry file that the run rejects, the check rejects too.
    const broken = await worker.validateTool(validateOf(argv, [memoryFile('wm.ini', '[BLAST\nBLAST_USAGE_REPORT=false\n'), files[1]!]));
    expect(broken.ok).toBe(false);
  });

  it('opens a folder given as `-in` and fails its read, as the command line does: "Is a directory"', async () => {
    const worker = runtime();
    const files = [memoryFile('counts/dm_u2.obinary', fixture('counts/dm_u2.obinary')), memoryFile('fa/dm2L_3k.fa', fixture('fasta/dm2L_3k.fa'))];
    const run = posts();
    const done = await worker.tool(toolOf(['windowmasker', '-ustat', 'counts/dm_u2.obinary', '-in', 'fa'], files), run);
    expect(done.status).toBe(3);
    expect(joined(run.stderr).toString()).toContain("Application's execution failed: basic_filebuf::xsgetn error reading the file: Is a directory");
  });

  it('runs makeblastdb in the writable directory: its database comes after the run, by name, and the input does not', async () => {
    const worker = runtime();
    const run = posts();
    const done = await worker.tool(toolOf(['makeblastdb', '-in', 'ns.fna', '-dbtype', 'nucl', '-parse_seqids', '-out', 'dbs/db'], [memoryFile('ns.fna', corpus('ns.fna'))]), run);
    expect(done.status).toBe(0);
    expect(joined(run.stdout).toString()).toMatch(/^\n\nBuilding a new DB, current time: \d\d\/\d\d\/\d{4} \d\d:\d\d:\d\d\n/);
    expect(joined(run.stderr).length).toBe(0);
    const files = filesOf(run.frames);
    // `-parse_seqids` of version 5: the ten files of the command line, in the folder that the run made.
    expect(files.map((file) => file.name)).toEqual(['ndb', 'nhr', 'nin', 'njs', 'nog', 'nos', 'not', 'nsq', 'ntf', 'nto'].map((extension) => `dbs/db.${extension}`));
    expect(files.every((file) => file.inOrder && file.bytes.length > 0)).toBe(true);
    expect(files.find((file) => file.name === 'dbs/db.nin')!.bytes.subarray(0, 4)).toEqual(Buffer.from([0, 0, 0, 5]));
    // A second run without the files: the run's directory went back to nothing, so the database of the first one is not there.
    const after = posts();
    const again = await worker.tool(toolOf(['makeblastdb', '-in', 'dbs/db', '-input_type', 'blastdb', '-dbtype', 'nucl', '-out', 'copy']), after);
    expect(again.status).not.toBe(0);
    expect(after.frames).toEqual([]);
  });

  it('posts a database file larger than a chunk in frames one after another, each at its offset (their buffers transferred)', async () => {
    const worker = runtime();
    // 1,500 proteins of 1,000 residues: a .psq of about 1.5 MB, more than one chunk (CHUNK_BYTES, 1 MiB) of the directory.
    const residues = 'ACDEFGHIKLMNPQRSTVWY';
    let fasta = '';
    for (let i = 0; i < 1500; i++) {
      let sequence = '';
      for (let j = 0; j < 1000; j++) sequence += residues[(i * 7 + j * 13 + ((i * j) % 17)) % residues.length];
      fasta += `>p${i} protein ${i}\n${sequence.replace(/(.{60})/g, '$1\n')}\n`;
    }
    const run = posts();
    const done = await worker.tool(toolOf(['makeblastdb', '-in', 'big.faa', '-dbtype', 'prot', '-out', 'big'], [memoryFile('big.faa', fasta)]), run);
    expect(done.status).toBe(0);
    const psq = filesOf(run.frames).find((file) => file.name === 'big.psq')!;
    expect(psq.inOrder).toBe(true);
    // One byte per residue, a 0 between the sequences and one at the start.
    expect(psq.bytes.length).toBe(1500 * 1001 + 1);
    const number = run.frames.find((frame) => frame.kind === 'begin' && frame.name === 'big.psq')!.number;
    const offsets = run.frames.flatMap((frame) => (frame.kind === 'bytes' && frame.number === number ? [frame.offset] : []));
    expect(offsets).toEqual([0, 1024 * 1024]);
  });

  it('replaces a database of the directory as the command line does: the old files it removed never come, and a log is appended to', async () => {
    const worker = runtime();
    // An older database of the same name made with `-parse_seqids` (its .nog and .nos), and a log of an earlier run.
    const older = ['ndb', 'nhr', 'nin', 'njs', 'nog', 'nos', 'not', 'nsq', 'ntf', 'nto'].map((extension) => memoryFile(`db.${extension}`, corpus(`ndbp.${extension}`)));
    const log = memoryFile('mk.log', 'an earlier run\n');
    const run = posts();
    const done = await worker.tool(toolOf(['makeblastdb', '-in', 'ns.fna', '-dbtype', 'nucl', '-out', 'db', '-logfile', 'mk.log'], [memoryFile('ns.fna', corpus('ns.fna')), ...older, log]), run);
    expect(done.status).toBe(0);
    const files = filesOf(run.frames);
    // The new database has no .nog and .nos: the run removed the old ones, which do not come.
    expect(files.map((file) => file.name)).toEqual(['db.ndb', 'db.nhr', 'db.nin', 'db.njs', 'db.not', 'db.nsq', 'db.ntf', 'db.nto', 'mk.log']);
    const logged = files.find((file) => file.name === 'mk.log')!.bytes.toString();
    expect(logged.startsWith('an earlier run\n\n\nBuilding a new DB, current time: ')).toBe(true);
    // The given files are the user's: their Blobs are as they were.
    expect(Buffer.from((log.file as MemoryBlob).data).toString()).toBe('an earlier run\n');
    expect(Buffer.from((older[4]!.file as MemoryBlob).data).equals(Buffer.from(corpus('ndbp.nog')))).toBe(true);
  });

  it("checks makeblastdb's words before its run: the argument reader's errors, and a valid argv whose input is missing", async () => {
    const worker = runtime();
    expect(await worker.validateTool(validateOf(['makeblastdb', '-in', 'x.fna', '-dbtype', 'nucl']))).toEqual({ ok: true });
    const refused = await worker.validateTool(validateOf(['makeblastdb', '-in', 'x.fna']));
    expect(refused.ok).toBe(false);
    expect(!refused.ok && refused.message).toContain('dbtype');
    // The check wrote nothing: a run of the same words in an empty directory leaves no file.
    const run = posts();
    expect((await worker.tool(toolOf(['makeblastdb', '-in', 'x.fna', '-dbtype', 'nucl', '-logfile', 'v.log']), run)).status).not.toBe(0);
    expect(filesOf(run.frames).map((file) => file.name)).toEqual(['v.log']);
  });
});

// --- the check of a tool's words does not disturb the runtime that the searches use (M2 of the second review) ---------------
// A fake loaded instance stands in for a threaded one (the built threaded reactor needs a cross-origin isolated page and
// workers, which Node lacks): the runtime's decisions are what is tested, with the calls of the engine recorded.

interface FakeLoaded {
  path: 'serial' | 'threaded';
  runs: number;
  retained: { key: string } | undefined;
  counts: { key: string } | undefined;
  ran: Set<string>;
  checkOnly: boolean;
  threads?: { terminate: () => boolean };
  instance: { abi: { stoppedBy: Error | undefined; validateTool: (argv: readonly string[]) => unknown }; wasi: { clear: () => void } };
  directory: { open: (files: readonly DirectoryFile[]) => void };
}

function fakeLoaded(path: 'serial' | 'threaded', log: string[], extra: Partial<FakeLoaded> = {}): FakeLoaded {
  const loaded: FakeLoaded = {
    path,
    runs: 0,
    retained: undefined,
    counts: undefined,
    ran: new Set(),
    checkOnly: false,
    ...(path === 'threaded' ? { threads: { terminate: () => (log.push('terminate threads'), true) } } : {}),
    instance: {
      abi: {
        stoppedBy: undefined,
        validateTool: (argv) => {
          log.push(`validate ${argv.join(' ')}`);
          return { ok: true };
        },
      },
      wasi: { clear: () => log.push('wasi clear') },
    },
    directory: { open: (files) => log.push(`directory ${files.map((file) => file.name).join(',')}`) },
    ...extra,
  };
  return loaded;
}

describe('a check of a tool on the runtime that the searches use', () => {
  const runtimeOf = (loaded: FakeLoaded, use?: InstanceUse) => {
    const runtime = new EngineRuntime({
      serial: () => Promise.reject(new Error('the serial module must not be loaded')),
      threads: () => Promise.reject(new Error('the threaded module must not be loaded')),
      builds: { serial: 'serial', threads: 'threads' },
      main: NO_LINKS,
      use,
    });
    (runtime as unknown as { loaded: FakeLoaded }).loaded = loaded;
    return runtime;
  };
  const check = (argv: string[], files?: DirectoryFile[]): WorkerToolValidate => ({ type: 'tool-validate', generation: 1, id: 'v', argv, ...(files === undefined ? {} : { files }) });

  it('leaves a threaded instance as it is: its thread workers, the subject and the counts file that it retains, its warm-up', async () => {
    const log: string[] = [];
    const retained = { key: 'subject' };
    const counts = { key: 'counts' };
    const loaded = fakeLoaded('threaded', log, { retained, counts, ran: new Set(['blastn']), runs: 2 });
    const use = InstanceUse.create();
    // The call into the threaded instance is a use, as a search's are.
    loaded.instance.abi.validateTool = (argv) => {
      log.push(`validate ${argv.join(' ')} (uses ${use.uses()})`);
      return { ok: true };
    };
    const runtime = runtimeOf(loaded, use);
    const wm = memoryFile('wm.ini', '[X]\n');
    expect(await runtime.validateTool(check(['windowmasker', '-conffile', 'wm.ini'], [wm]))).toEqual({ ok: true });
    expect(log).toEqual(['directory wm.ini', 'wasi clear', 'validate windowmasker -conffile wm.ini (uses 1)', 'directory ']);
    expect(use.uses()).toBe(0);
    const held = (runtime as unknown as { loaded: FakeLoaded }).loaded;
    expect(held).toBe(loaded);
    expect(held.path).toBe('threaded');
    expect(held.retained).toBe(retained);
    expect(held.counts).toBe(counts);
    expect([...held.ran]).toEqual(['blastn']);
    expect(held.runs).toBe(2);
  });

  it('does not call into a threaded instance that the Gateway is ending, and lets the directory go', async () => {
    const log: string[] = [];
    const loaded = fakeLoaded('threaded', log);
    const use = InstanceUse.create();
    use.retire();
    await expect(runtimeOf(loaded, use).validateTool(check(['windowmasker']))).rejects.toThrow('the runtime is ending');
    expect(log).toEqual(['directory ']);
  });

  it('leaves a stopped threaded instance to the Gateway, which ends the worker; a stopped serial one is dropped', async () => {
    const log: string[] = [];
    const threaded = fakeLoaded('threaded', log);
    threaded.instance.abi.validateTool = () => {
      threaded.instance.abi.stoppedBy = new Error('stopped');
      throw threaded.instance.abi.stoppedBy;
    };
    const runtime = runtimeOf(threaded, InstanceUse.create());
    await expect(runtime.validateTool(check(['windowmasker']))).rejects.toThrow('stopped');
    // The thread workers and the keeper stay until the worker ends (as after a search).
    expect(log).not.toContain('terminate threads');
    expect((runtime as unknown as { loaded: FakeLoaded }).loaded).toBe(threaded);
    const serial = fakeLoaded('serial', log);
    serial.instance.abi.validateTool = () => {
      serial.instance.abi.stoppedBy = new Error('stopped');
      throw serial.instance.abi.stoppedBy;
    };
    const again = runtimeOf(serial);
    await expect(again.validateTool(check(['windowmasker']))).rejects.toThrow('stopped');
    expect((again as unknown as { loaded: FakeLoaded | undefined }).loaded).toBeUndefined();
  });

  it('does not keep a serial instance that only a check loaded from the preparation for a threaded search, and keeps one that ran a search', async () => {
    const log: string[] = [];
    // The check loaded a serial instance and nothing used it: the preparation for threads replaces it (here the page cannot run
    // threaded, so the attempt shows as that error).
    const idle = runtimeOf(fakeLoaded('serial', log, { checkOnly: true }));
    await expect(idle.prepare(4)).rejects.toThrow(/cross-origin isolated/);
    // One that a preparation loaded stays too (the application's test of `prepare` pins it), though nothing used it.
    const prepared = runtimeOf(fakeLoaded('serial', log));
    await expect(prepared.prepare(4)).resolves.toBeUndefined();
    // A serial instance that ran something keeps its place: a queued search may still use it.
    const used = runtimeOf(fakeLoaded('serial', log, { checkOnly: true, runs: 1 }));
    await expect(used.prepare(4)).resolves.toBeUndefined();
    const warmed = runtimeOf(fakeLoaded('serial', log, { checkOnly: true, ran: new Set(['blastp']) }));
    await expect(warmed.prepare(4)).resolves.toBeUndefined();
  });
});
