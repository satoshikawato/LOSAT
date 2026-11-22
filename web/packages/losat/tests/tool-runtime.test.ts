// The Engine worker's side of a tool run on the built serial reactor, in Node (review L8, M3 of WW3a, 2026-11-23): the
// `EngineRuntime.tool` that the worker calls (the run's directory opened and let go, the retained input released, a second run
// on the same instance), the worker's check of a tool's words with the run's directory (`EngineRuntime.validateTool`, which
// reads `-conffile`), and the command line's answer for `-in <folder>` (the open of a folder succeeds, its read fails with
// EISDIR). Runs with LOSAT_WEB_REACTORS (the folder of web/adapter/tools/build_reactors.py's output); these rows read no
// frame of stream 12, so the reactors of an earlier frame layout serve as well. The directory reads a file's bytes through
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
const NO_LINKS: MainLinks = {
  link: () => {
    throw new Error('unused');
  },
  fault: () => undefined,
};

/** What a run posted. */
function posts(): ToolPosts & { stdout: Uint8Array[]; stderr: Uint8Array[]; frames: ToolFrame[] } {
  const posted = { stdout: [] as Uint8Array[], stderr: [] as Uint8Array[], frames: [] as ToolFrame[] };
  return { ...posted, chunk: (stream, bytes) => posted[stream].push(bytes), file: (frame) => posted.frames.push(frame) };
}
const joined = (chunks: Uint8Array[]) => Buffer.concat(chunks);

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
    const registered = loaded.instance.abi.register('blastp', ROLE_SUBJECT, false, new TextEncoder().encode('>s\nMKVLAGG\n'));
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
});
