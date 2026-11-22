// A tool run through the package (`WasmEngine.runTool`, ABI v2 version 13, `losat_web2_tool`; line NP WW3a, D-WW3-6,
// 2026-11-22): the frames of stream 12, the router that hands the tool's streams to the caller's sinks or collects them
// as Blobs, `ReactorAbi.tool` and `validateTool` over a fake engine, the Gateway's run with a fake Engine worker (the
// sinks take the chunks before the run's result), and, with LOSAT_WEB_REACTORS (the folder of
// web/adapter/tools/build_reactors.py's output, as database.test.ts), `windowmasker` on the serial reactor in Node.
import { closeSync, openAsBlob, openSync, readFileSync, readSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it, vi } from 'vitest';
import { FRAME_BEGIN, FRAME_BYTES, FRAME_END, ownedFrame, parseToolFrame, type ToolFrame } from '../src/domain/tool-frames.ts';
import { RunCancelledError, type EngineGateway, type ToolFileSink } from '../src/ports/engine.ts';
import { ABI_VERSION, EngineCallError, ReactorAbi, type AbiExports } from '../src/infra/reactor/abi.ts';
import { ABI_EXPORTS } from '../src/infra/reactor/artifact.ts';
import { compileReactor, instantiateSerial } from '../src/infra/reactor/instance.ts';
import { RunDirectory } from '../src/infra/reactor/run-directory.ts';
import { reactorControl } from '../src/infra/reactor/control.ts';
import { WasmEngine } from '../src/infra/engine-worker/gateway.ts';
import { BLOB_FOLD_BYTES, ToolOutputRouter } from '../src/infra/engine-worker/tool-output.ts';
import type { EngineCommand, EngineEvent, WorkerTool } from '../src/infra/engine-worker/protocol.ts';

const encode = (text: string) => new TextEncoder().encode(text);
const text = (bytes: Uint8Array) => new TextDecoder().decode(bytes);

/** A frame of stream 12 as the engine emits it: the file's number, the kind, the payload. */
function frame(number: number, kind: number, payload: Uint8Array = new Uint8Array(0)): Uint8Array {
  const bytes = new Uint8Array(8 + payload.length);
  const view = new DataView(bytes.buffer);
  view.setUint32(0, number, true);
  view.setUint32(4, kind, true);
  bytes.set(payload, 8);
  return bytes;
}

describe('the frames of stream 12', () => {
  it('reads the beginning, the bytes and the end of a file', () => {
    expect(parseToolFrame(frame(0, FRAME_BEGIN, encode('c.counts')))).toEqual({ number: 0, kind: 'begin', name: 'c.counts' });
    expect(parseToolFrame(frame(7, FRAME_BYTES, new Uint8Array([1, 2, 3])))).toEqual({ number: 7, kind: 'bytes', bytes: new Uint8Array([1, 2, 3]) });
    expect(parseToolFrame(frame(0xffff_ffff, FRAME_END))).toEqual({ number: 0xffff_ffff, kind: 'end' });
    // A name may be any UTF-8, and a file of no bytes has a frame of no bytes.
    expect(parseToolFrame(frame(1, FRAME_BEGIN, encode('dir/名前.txt')))).toEqual({ number: 1, kind: 'begin', name: 'dir/名前.txt' });
    expect(parseToolFrame(frame(1, FRAME_BYTES))).toEqual({ number: 1, kind: 'bytes', bytes: new Uint8Array(0) });
  });

  it('reads a frame that is a view into larger memory, and owns a copy of its bytes', () => {
    const memory = new Uint8Array(64);
    memory.set(frame(2, FRAME_BYTES, new Uint8Array([9, 8])), 20);
    const parsed = parseToolFrame(memory.subarray(20, 30));
    expect(parsed).toEqual({ number: 2, kind: 'bytes', bytes: new Uint8Array([9, 8]) });
    const owned = ownedFrame(parsed);
    memory.fill(0);
    expect(owned).toEqual({ number: 2, kind: 'bytes', bytes: new Uint8Array([9, 8]) });
    expect(parsed).toEqual({ number: 2, kind: 'bytes', bytes: new Uint8Array([0, 0]) });
  });

  it('refuses a frame that breaks the ABI', () => {
    expect(() => parseToolFrame(new Uint8Array(7))).toThrow(/7 bytes, fewer than its 8/);
    expect(() => parseToolFrame(frame(0, 3))).toThrow(/kind 3/);
    expect(() => parseToolFrame(frame(0, FRAME_END, new Uint8Array(1)))).toThrow(/end of file 0 .* 1 bytes of payload/);
    expect(() => parseToolFrame(frame(4, FRAME_BEGIN, new Uint8Array([0xff, 0xfe])))).toThrow(/name of file 4 .* not UTF-8/);
  });
});

/** A sink that logs what it is given, with a delay in each call when `slow`. */
function loggingSink(log: string[], slow = false): ToolFileSink {
  const pause = () => (slow ? new Promise<void>((resolve) => setTimeout(resolve, 5)) : undefined);
  return {
    write: async (chunk) => {
      log.push(`write ${text(chunk)} begin`);
      await pause();
      log.push(`write ${text(chunk)} end`);
    },
    close: async () => {
      log.push('close begin');
      await pause();
      log.push('close end');
    },
    abort: (reason) => {
      log.push(`abort ${String(reason)}`);
    },
  };
}

describe('the output router of a tool run', () => {
  it('collects what has no sink as Blobs: the streams, and the files in the order of their numbers', async () => {
    const router = new ToolOutputRouter(undefined, () => undefined);
    router.stdout(encode('out1 '));
    router.stderr(encode('err'));
    router.stdout(encode('out2'));
    router.frame({ number: 0, kind: 'begin', name: 'log.txt' });
    router.frame({ number: 1, kind: 'begin', name: 'c.counts' });
    router.frame({ number: 1, kind: 'bytes', bytes: encode('ab') });
    router.frame({ number: 0, kind: 'bytes', bytes: encode('log') });
    router.frame({ number: 1, kind: 'bytes', bytes: encode('cd') });
    // The files end in the other order; two openings of one name are two files.
    router.frame({ number: 1, kind: 'end' });
    router.frame({ number: 2, kind: 'begin', name: 'log.txt' });
    router.frame({ number: 2, kind: 'end' });
    router.frame({ number: 0, kind: 'end' });
    const collected = await router.finish();
    expect(await collected.stdout.text()).toBe('out1 out2');
    expect(await collected.stderr.text()).toBe('err');
    expect(collected.files.map((file) => file.name)).toEqual(['log.txt', 'c.counts', 'log.txt']);
    expect(await Promise.all(collected.files.map((file) => file.file.text()))).toEqual(['log', 'abcd', '']);
  });

  it('folds a large collected stream into its Blob without joining it', async () => {
    const router = new ToolOutputRouter(undefined, () => undefined);
    const chunk = new Uint8Array(BLOB_FOLD_BYTES / 2 + 1).fill(7);
    for (let i = 0; i < 3; i++) router.stdout(chunk.slice());
    router.stdout(encode('end'));
    const { stdout } = await router.finish();
    expect(stdout.size).toBe(3 * chunk.length + 3);
    expect(await stdout.slice(3 * chunk.length).text()).toBe('end');
  });

  it('gives a stream with a sink to the sink, in order and one call at a time, and keeps nothing of it', async () => {
    const log: string[] = [];
    const router = new ToolOutputRouter(
      {
        onStdout: async (chunk) => {
          log.push(`stdout ${text(chunk)}`);
          await new Promise((resolve) => setTimeout(resolve, 5));
        },
        onStderr: (chunk) => void log.push(`stderr ${text(chunk)}`),
        onFile: (name) => {
          log.push(`open ${name}`);
          return loggingSink(log, true);
        },
      },
      () => undefined,
    );
    router.stdout(encode('a'));
    router.stderr(encode('b'));
    router.frame({ number: 0, kind: 'begin', name: 'f' });
    router.frame({ number: 0, kind: 'bytes', bytes: encode('x') });
    router.stdout(encode('c'));
    router.frame({ number: 0, kind: 'bytes', bytes: encode('y') });
    router.frame({ number: 0, kind: 'end' });
    const collected = await router.finish();
    expect(log).toEqual([
      'stdout a', 'stderr b', 'open f', 'write x begin', 'write x end', 'stdout c', 'write y begin', 'write y end', 'close begin', 'close end',
    ]);
    expect(collected.stdout.size).toBe(0);
    expect(collected.stderr.size).toBe(0);
    expect(collected.files).toEqual([]);
  });

  it('collects the file that its sink declines, and streams the others', async () => {
    const log: string[] = [];
    const router = new ToolOutputRouter({ onFile: (name) => (name === 'big' ? loggingSink(log) : undefined) }, () => undefined);
    router.frame({ number: 0, kind: 'begin', name: 'big' });
    router.frame({ number: 1, kind: 'begin', name: 'small' });
    router.frame({ number: 1, kind: 'bytes', bytes: encode('s') });
    router.frame({ number: 0, kind: 'bytes', bytes: encode('B') });
    router.frame({ number: 0, kind: 'end' });
    router.frame({ number: 1, kind: 'end' });
    const { files } = await router.finish();
    expect(files.map((file) => file.name)).toEqual(['small']);
    expect(log).toEqual(['write B begin', 'write B end', 'close begin', 'close end']);
  });

  it("fails with a sink's error: it drops what is queued, tells the open files, and says so once", async () => {
    const log: string[] = [];
    const failed = vi.fn();
    const router = new ToolOutputRouter(
      {
        onFile: () => ({
          ...loggingSink(log),
          write: (chunk) => {
            if (text(chunk) === 'bad') throw new Error('disk full');
          },
        }),
      },
      failed,
    );
    router.frame({ number: 0, kind: 'begin', name: 'f' });
    router.frame({ number: 0, kind: 'bytes', bytes: encode('bad') });
    router.frame({ number: 0, kind: 'bytes', bytes: encode('later') });
    router.frame({ number: 0, kind: 'end' });
    await expect(router.finish()).rejects.toThrow('disk full');
    expect(failed).toHaveBeenCalledTimes(1);
    expect(log).toEqual(['abort Error: disk full']);
  });

  it('fails for frames that no run sends: a file twice, an unknown file, bytes after the end, a file left open', async () => {
    const run = async (frames: ToolFrame[]) => {
      const router = new ToolOutputRouter(undefined, () => undefined);
      for (const next of frames) router.frame(next);
      return router.finish();
    };
    const begin = (number: number): ToolFrame => ({ number, kind: 'begin', name: 'f' });
    await expect(run([begin(0), begin(0)])).rejects.toThrow(/began its file 0 twice/);
    await expect(run([{ number: 3, kind: 'end' }])).rejects.toThrow(/file 3, which it did not begin/);
    await expect(run([begin(0), { number: 0, kind: 'end' }, { number: 0, kind: 'bytes', bytes: encode('x') }])).rejects.toThrow(/after its end/);
    await expect(run([begin(0)])).rejects.toThrow(/file 0 \(f\) still open/);
  });

  it('aborts: the sinks of open files hear the reason, and what is queued is dropped', async () => {
    const log: string[] = [];
    const router = new ToolOutputRouter({ onFile: () => loggingSink(log, true) }, () => undefined);
    router.frame({ number: 0, kind: 'begin', name: 'f' });
    router.frame({ number: 0, kind: 'bytes', bytes: encode('x') });
    await vi.waitFor(() => expect(log).toEqual(['write x begin', 'write x end']));
    const aborted = router.abort('cancelled');
    router.frame({ number: 0, kind: 'bytes', bytes: encode('y') });
    await aborted;
    expect(log).toEqual(['write x begin', 'write x end', 'abort cancelled']);
    await expect(router.finish()).rejects.toBe('cancelled');
  });
});

// --- ReactorAbi over a fake engine -----------------------------------------------------------------------------------

type Emit = (stream: number, bytes: Uint8Array) => void;

/** A fake engine whose `losat_web2_tool` and `losat_web2_validate` run the given scripts, and log their argv. */
function fakeEngine(options: { tool?: (emit: Emit, argv: string[]) => number; validate?: (emit: Emit, argv: string[]) => number; error?: string } = {}) {
  const memory = new WebAssembly.Memory({ initial: 8 });
  const log: string[] = [];
  const ref: { abi?: ReactorAbi } = {};
  let next = 4096;
  const emit: Emit = (stream, bytes) => {
    new Uint8Array(memory.buffer, 100_000, bytes.length).set(bytes);
    ref.abi!.emit(stream, 100_000, bytes.length);
  };
  const argvOf = (ptr: number, len: number) => text(new Uint8Array(memory.buffer, ptr, len).slice()).split('\0');
  const exports = {
    memory,
    losat_web2_abi_version: () => ABI_VERSION,
    losat_web2_alloc: (len: number) => {
      const ptr = next;
      next += len;
      return ptr;
    },
    losat_web2_dealloc: () => 0,
    losat_web2_last_error_ptr: () => {
      const bytes = encode(options.error ?? '');
      new Uint8Array(memory.buffer, 90_000, bytes.length).set(bytes);
      return 90_000;
    },
    losat_web2_last_error_len: () => encode(options.error ?? '').length,
    losat_web2_tool: (ptr: number, len: number) => {
      const argv = argvOf(ptr, len);
      log.push(`tool ${argv.join(' ')}`);
      return options.tool?.(emit, argv) ?? 0;
    },
    losat_web2_validate: (ptr: number, len: number) => {
      const argv = argvOf(ptr, len);
      log.push(`validate ${argv.join(' ')}`);
      return options.validate?.(emit, argv) ?? 0;
    },
  };
  const abi = new ReactorAbi(exports as unknown as AbiExports);
  ref.abi = abi;
  return { abi, log };
}

describe('ReactorAbi.tool and validateTool', () => {
  it('speaks ABI version 13 and the module check wants losat_web2_tool', () => {
    expect(ABI_VERSION).toBe(13);
    expect(ABI_EXPORTS).toContain('losat_web2_tool');
  });

  it('hands streams 11, 3 and 12 to the callbacks as each emit arrives, and returns the status of stream 2', () => {
    const calls: string[] = [];
    const { abi, log } = fakeEngine({
      tool: (emit) => {
        emit(11, encode('out'));
        emit(12, frame(0, FRAME_BEGIN, encode('log.txt')));
        emit(3, encode('err'));
        emit(12, frame(1, FRAME_BEGIN, encode('c.counts')));
        emit(12, frame(1, FRAME_BYTES, encode('abc')));
        emit(12, frame(0, FRAME_BYTES, encode('L')));
        emit(12, frame(1, FRAME_END));
        emit(12, frame(0, FRAME_END));
        emit(2, encode('{"status":3}'));
        return 0;
      },
    });
    const kept: Uint8Array[] = [];
    const status = abi.tool(['windowmasker', '-mk_counts', '-in', 'x.fa'], {
      stdout: (bytes) => calls.push(`stdout ${text(bytes)}`),
      stderr: (bytes) => calls.push(`stderr ${text(bytes)}`),
      file: (parsed) => {
        if (parsed.kind === 'bytes') kept.push(parsed.bytes);
        calls.push(`file ${parsed.number} ${parsed.kind}${parsed.kind === 'begin' ? ` ${parsed.name}` : parsed.kind === 'bytes' ? ` ${text(parsed.bytes)}` : ''}`);
      },
    });
    expect(status).toBe(3);
    expect(log).toEqual(['tool windowmasker -mk_counts -in x.fa']);
    // Each call came as the engine emitted it, in order, and the status stream 2 is not a callback's.
    expect(calls).toEqual([
      'stdout out',
      'file 0 begin log.txt',
      'stderr err',
      'file 1 begin c.counts',
      'file 1 bytes abc',
      'file 0 bytes L',
      'file 1 end',
      'file 0 end',
    ]);
    expect(kept).toHaveLength(2);
  });

  it('keeps nothing of the streams, and keeps no sink after the call', () => {
    const { abi } = fakeEngine({
      tool: (emit) => {
        emit(11, encode('out'));
        emit(2, encode('{"status":0}'));
        return 0;
      },
    });
    let seen = 0;
    abi.tool(['windowmasker'], { stdout: () => seen++, stderr: () => seen++, file: () => seen++ });
    expect(seen).toBe(1);
    // An emit after the run reaches no callback.
    expect(() => abi.emit(11, 0, 0)).not.toThrow();
    expect(seen).toBe(1);
  });

  it("throws the engine's message for a name that is not a tool, and for a run that sent no status", () => {
    const refused = fakeEngine({ tool: () => -1, error: 'blastp is not available in the browser (ABI v2)' });
    expect(() => refused.abi.tool(['blastp'], { stdout: () => undefined, stderr: () => undefined, file: () => undefined })).toThrow(EngineCallError);
    expect(() => refused.abi.tool(['blastp'], { stdout: () => undefined, stderr: () => undefined, file: () => undefined })).toThrow(
      'blastp is not available in the browser (ABI v2)',
    );
    const none = { stdout: () => undefined, stderr: () => undefined, file: () => undefined };
    expect(() => fakeEngine({ tool: () => 0 }).abi.tool(['windowmasker'], none)).toThrow('the engine sent no response');
    const noStatus = fakeEngine({ tool: (emit) => (emit(2, encode('{"files":1}')), 0) });
    expect(() => noStatus.abi.tool(['windowmasker'], none)).toThrow(/no exit status/);
  });

  it('validates the words of a tool: valid, or the command line\'s text; no response is read', () => {
    const valid = fakeEngine();
    expect(valid.abi.validateTool(['windowmasker', '-mk_counts', '-in', 'x.fa'])).toEqual({ ok: true });
    expect(valid.log).toEqual(['validate windowmasker -mk_counts -in x.fa']);
    // `validate` of a search reads a response; for a tool's argv the engine emits none, and the call does not wait for one.
    const refused = fakeEngine({ validate: () => -1, error: "USAGE\n====\nError: one of '-mk_counts', '-convert' or '-ustat <stat_file>' must be specified" });
    expect(refused.abi.validateTool(['windowmasker', '-in', 'x.fa'])).toEqual({
      ok: false,
      message: "USAGE\n====\nError: one of '-mk_counts', '-convert' or '-ustat <stat_file>' must be specified",
    });
  });

  it("answers validateTool through reactorControl's shared instance", async () => {
    const { abi } = fakeEngine({ validate: (_emit, argv) => (argv.includes('-bad') ? -1 : 0), error: 'bad argument' });
    const control = reactorControl(async () => abi, async () => abi);
    expect(await control.validateTool(['windowmasker', '-mk_counts'])).toEqual({ ok: true });
    expect(await control.validateTool(['windowmasker', '-bad'])).toEqual({ ok: false, message: 'bad argument' });
  });
});

// --- the Gateway with a fake Engine worker ----------------------------------------------------------------------------

class FakeWorker {
  readonly posted: EngineCommand[] = [];
  terminated = false;
  onmessage: ((event: MessageEvent<EngineEvent>) => void) | null = null;
  onerror: ((event: ErrorEvent) => void) | null = null;
  onmessageerror: (() => void) | null = null;

  postMessage(command: EngineCommand): void {
    this.posted.push(command);
  }

  terminate(): void {
    this.terminated = true;
  }

  get tools(): WorkerTool[] {
    return this.posted.filter((command): command is WorkerTool => command.type === 'tool');
  }

  get generation(): number {
    return this.posted.find((command) => command.type === 'init')!.generation;
  }

  send(event: EngineEvent): void {
    this.onmessage?.({ data: event } as MessageEvent<EngineEvent>);
  }
}

const ASSETS = {
  serial: { url: '/serial.wasm', sha256: 'a'.repeat(64), size: 1, artifactSha256: 'a'.repeat(64) },
  threads: { url: '/threads.wasm', sha256: 'b'.repeat(64), size: 1, artifactSha256: 'c'.repeat(64), transform: 'shared-memory guard' },
};

function setup(validateTool?: (argv: readonly string[]) => Promise<{ ok: true } | { ok: false; message: string }>) {
  const workers: FakeWorker[] = [];
  const engine = new WasmEngine({
    assets: ASSETS,
    control: {
      describe: () => Promise.reject(new Error('unused')),
      start: () => Promise.reject(new Error('unused')),
      validate: () => Promise.reject(new Error('unused')),
      ...(validateTool === undefined ? {} : { validateTool }),
    },
    hardwareConcurrency: 8,
    createWorker: () => {
      const worker = new FakeWorker();
      workers.push(worker);
      return worker as unknown as Worker;
    },
    loadModules: async () => ({ serial: { asset: ASSETS.serial }, threads: { asset: ASSETS.threads } }),
  });
  return { engine, workers };
}

const MEMORY = { linearBytesBefore: 1, linearBytesAfter: 2, instanceRuns: 1 };

/** Starts a tool run and waits until the worker has the command. */
async function started(engine: WasmEngine, workers: FakeWorker[], ...args: Parameters<WasmEngine['runTool']>) {
  const promise = engine.runTool(...args);
  await vi.waitFor(() => expect(workers[0]?.tools).toHaveLength(1));
  const worker = workers[0]!;
  const { runId, generation } = worker.tools[0]!;
  return {
    promise,
    worker,
    chunk: (stream: 'stdout' | 'stderr', value: string) => worker.send({ type: 'tool-chunk', generation, runId, stream, bytes: encode(value) }),
    file: (parsed: ToolFrame) => worker.send({ type: 'tool-file', generation, runId, frame: parsed }),
    done: (status: number) => worker.send({ type: 'tool-done', generation, runId, result: { status, memory: MEMORY } }),
    fail: (message: string, stopped = false) =>
      worker.send({ type: 'done', generation, runId, ok: false, error: { name: 'EngineCallError', message, ...(stopped ? { stopped: true } : {}) } }),
  };
}

describe('WasmEngine.runTool', () => {
  it('posts the argv and the directory to the worker, and resolves with the status and the collected streams', async () => {
    const { engine, workers } = setup();
    const file = new Blob(['>x\nACGT\n']);
    const run = await started(engine, workers, ['windowmasker', '-mk_counts', '-in', 'x.fa', '-out', 'c.counts'], { directory: [{ name: 'x.fa', file }] });
    expect(run.worker.tools[0]).toMatchObject({ type: 'tool', argv: ['windowmasker', '-mk_counts', '-in', 'x.fa', '-out', 'c.counts'], files: [{ name: 'x.fa', file }] });
    run.chunk('stdout', 'hello ');
    run.chunk('stderr', 'warn');
    run.file({ number: 0, kind: 'begin', name: 'c.counts' });
    run.file({ number: 0, kind: 'bytes', bytes: encode('12') });
    run.chunk('stdout', 'world');
    run.file({ number: 0, kind: 'bytes', bytes: encode('34') });
    run.file({ number: 0, kind: 'end' });
    run.done(0);
    const result = await run.promise;
    expect(result.status).toBe(0);
    expect(await result.stdout.text()).toBe('hello world');
    expect(await result.stderr.text()).toBe('warn');
    expect(result.files.map((collected) => collected.name)).toEqual(['c.counts']);
    expect(await result.files[0]!.file.text()).toBe('1234');
  });

  it('gives each file chunk to its sink before the result: the run does not resolve until the sinks have taken everything', async () => {
    const { engine, workers } = setup();
    const log: string[] = [];
    let release!: () => void;
    const gate = new Promise<void>((resolve) => (release = resolve));
    const run = await started(engine, workers, ['windowmasker', '-mk_counts'], {
      sinks: {
        onStdout: (chunk) => void log.push(`stdout ${text(chunk)}`),
        onFile: (name) => {
          log.push(`open ${name}`);
          return {
            write: async (chunk) => {
              // The first write is slow: the next frames wait for it, and the result for them all.
              if (log.filter((entry) => entry.startsWith('write')).length === 0) await gate;
              log.push(`write ${text(chunk)}`);
            },
            close: () => void log.push('close'),
          };
        },
      },
    });
    let resolved = false;
    void run.promise.then(() => (resolved = true));
    run.file({ number: 0, kind: 'begin', name: 'c.counts' });
    run.file({ number: 0, kind: 'bytes', bytes: encode('ab') });
    run.chunk('stdout', 'o');
    run.file({ number: 0, kind: 'bytes', bytes: encode('cd') });
    run.file({ number: 0, kind: 'end' });
    run.done(0);
    await vi.waitFor(() => expect(log).toEqual(['open c.counts']));
    expect(resolved).toBe(false);
    release();
    const result = await run.promise;
    expect(log).toEqual(['open c.counts', 'write ab', 'stdout o', 'write cd', 'close']);
    // What a sink took is not in the result.
    expect(result.files).toEqual([]);
    expect(result.stdout.size).toBe(0);
    expect(await result.stderr.text()).toBe('');
  });

  it('shows two writes at the sink while the run has no result yet', async () => {
    const { engine, workers } = setup();
    const writes: string[] = [];
    const run = await started(engine, workers, ['windowmasker', '-mk_counts'], {
      sinks: { onFile: () => ({ write: (chunk) => void writes.push(text(chunk)), close: () => undefined }) },
    });
    let resolved = false;
    void run.promise.then(() => (resolved = true));
    run.file({ number: 0, kind: 'begin', name: 'f' });
    run.file({ number: 0, kind: 'bytes', bytes: encode('one') });
    run.file({ number: 0, kind: 'bytes', bytes: encode('two') });
    run.file({ number: 0, kind: 'end' });
    await vi.waitFor(() => expect(writes).toEqual(['one', 'two']));
    expect(resolved).toBe(false);
    run.done(3);
    expect((await run.promise).status).toBe(3);
  });

  it('runs one at a time, a search included', async () => {
    const { engine, workers } = setup();
    const run = await started(engine, workers, ['windowmasker']);
    await expect(engine.runTool(['windowmasker'])).rejects.toThrow('the engine runs one search at a time');
    run.done(0);
    await run.promise;
    // The next run is free.
    const next = engine.runTool(['windowmasker']);
    await vi.waitFor(() => expect(workers[0]!.tools).toHaveLength(2));
    const second = workers[0]!.tools[1]!;
    workers[0]!.send({ type: 'tool-done', generation: second.generation, runId: second.runId, result: { status: 1, memory: MEMORY } });
    expect((await next).status).toBe(1);
  });

  it("rejects with the engine's error (a name that is not a tool) and keeps the worker; the next run is free", async () => {
    const { engine, workers } = setup();
    const run = await started(engine, workers, ['blastp']);
    run.fail('blastp is not available in the browser (ABI v2)');
    await expect(run.promise).rejects.toThrow('blastp is not available in the browser (ABI v2)');
    expect(run.worker.terminated).toBe(false);
    const next = engine.runTool(['windowmasker']);
    await vi.waitFor(() => expect(workers[0]!.tools).toHaveLength(2));
    const command = workers[0]!.tools[1]!;
    workers[0]!.send({ type: 'tool-done', generation: command.generation, runId: command.runId, result: { status: 0, memory: MEMORY } });
    expect((await next).status).toBe(0);
  });

  it('ends the worker of an instance that stopped', async () => {
    const { engine, workers } = setup();
    const run = await started(engine, workers, ['windowmasker']);
    run.fail('The engine stopped: unreachable', true);
    await expect(run.promise).rejects.toThrow('The engine stopped: unreachable');
    expect(run.worker.terminated).toBe(true);
  });

  it('cancels with a signal: the run rejects, the worker ends, and the sinks of open files hear it', async () => {
    const { engine, workers } = setup();
    const controller = new AbortController();
    const log: string[] = [];
    const run = await started(engine, workers, ['windowmasker', '-mk_counts'], {
      signal: controller.signal,
      sinks: { onFile: () => loggingSink(log) },
    });
    run.file({ number: 0, kind: 'begin', name: 'f' });
    run.file({ number: 0, kind: 'bytes', bytes: encode('x') });
    await vi.waitFor(() => expect(log).toEqual(['write x begin', 'write x end']));
    controller.abort();
    await expect(run.promise).rejects.toBeInstanceOf(RunCancelledError);
    expect(run.worker.terminated).toBe(true);
    await vi.waitFor(() => expect(log).toHaveLength(3));
    expect(log[2]).toMatch(/^abort RunCancelledError/);
  });

  it('rejects at once for a signal that is already aborted', async () => {
    const { engine } = setup();
    await expect(engine.runTool(['windowmasker'], { signal: AbortSignal.abort() })).rejects.toBeInstanceOf(RunCancelledError);
  });

  it("ends the run with a sink's error: the worker ends, and the error is the sink's", async () => {
    const { engine, workers } = setup();
    const log: string[] = [];
    const run = await started(engine, workers, ['windowmasker', '-mk_counts'], {
      sinks: {
        onFile: () => ({
          ...loggingSink(log),
          write: () => {
            throw new Error('the file cannot be written');
          },
        }),
      },
    });
    run.file({ number: 0, kind: 'begin', name: 'f' });
    run.file({ number: 0, kind: 'bytes', bytes: encode('x') });
    await expect(run.promise).rejects.toThrow('the file cannot be written');
    expect(run.worker.terminated).toBe(true);
    expect(log).toEqual(['abort Error: the file cannot be written']);
  });

  it('fails when the worker stops, and a later run starts a new worker', async () => {
    const { engine, workers } = setup();
    const run = await started(engine, workers, ['windowmasker']);
    run.worker.send({ type: 'started', generation: run.worker.generation });
    run.worker.onerror?.({ preventDefault: () => undefined, message: 'boom' } as unknown as ErrorEvent);
    await expect(run.promise).rejects.toThrow('boom');
    const next = engine.runTool(['windowmasker']);
    await vi.waitFor(() => expect(workers).toHaveLength(2));
    await vi.waitFor(() => expect(workers[1]!.tools).toHaveLength(1));
    const command = workers[1]!.tools[0]!;
    workers[1]!.send({ type: 'tool-done', generation: command.generation, runId: command.runId, result: { status: 0, memory: MEMORY } });
    expect((await next).status).toBe(0);
  });

  it('renews the runtime before the next run when the memory that the tool left is over the limit', async () => {
    const { engine, workers } = setup();
    const run = await started(engine, workers, ['windowmasker']);
    run.worker.send({
      type: 'tool-done',
      generation: run.worker.generation,
      runId: run.worker.tools[0]!.runId,
      result: { status: 0, memory: { linearBytesBefore: 0, linearBytesAfter: 2 ** 32, instanceRuns: 1 } },
    });
    await run.promise;
    void engine.runTool(['windowmasker']);
    await vi.waitFor(() => expect(workers).toHaveLength(2));
    expect(workers[0]!.terminated).toBe(true);
  });

  it('is an EngineGateway that also validates tools through its control', async () => {
    const validated: string[][] = [];
    const { engine } = setup(async (argv) => {
      validated.push([...argv]);
      return argv.includes('-bad') ? { ok: false, message: 'Error: bad' } : { ok: true };
    });
    const gateway: EngineGateway = engine;
    expect(gateway).toBe(engine);
    expect(await engine.validateTool(['windowmasker', '-mk_counts'])).toEqual({ ok: true });
    expect(await engine.validateTool(['windowmasker', '-bad'])).toEqual({ ok: false, message: 'Error: bad' });
    expect(validated).toEqual([['windowmasker', '-mk_counts'], ['windowmasker', '-bad']]);
    await expect(setup().engine.validateTool(['windowmasker'])).rejects.toThrow(/does not validate tools/);
  });
});

// --- windowmasker on the built serial reactor, in Node ----------------------------------------------------------------

const LOSAT = fileURLToPath(new URL('../../../../LOSAT', import.meta.url));
const REACTORS = process.env['LOSAT_WEB_REACTORS'];

describe.skipIf(REACTORS === undefined)('windowmasker through the host', () => {
  it('counts a FASTA file into c.counts and log.txt, and refuses the words that the command line refuses', async () => {
    const { module } = await compileReactor(readFileSync(join(REACTORS!, 'losat-web-serial.wasm')));
    const path = join(LOSAT, 'tests/winmask/fasta/dm2L_3k.fa');
    const fasta = await openAsBlob(path);
    // The run's directory reads a Blob's bytes from its file, as FileReaderSync reads a File in a worker.
    const directory = new RunDirectory(undefined, (blob, start, end) => {
      expect(blob).toBe(fasta);
      const bytes = new Uint8Array(end - start);
      const descriptor = openSync(path, 'r');
      try {
        return bytes.subarray(0, readSync(descriptor, bytes, 0, bytes.length, start));
      } finally {
        closeSync(descriptor);
      }
    });
    const { abi } = await instantiateSerial(module, directory);
    const run = async (argv: string[]) => {
      const router = new ToolOutputRouter(undefined, (error) => {
        throw error;
      });
      const status = abi.tool(argv, {
        stdout: (bytes) => router.stdout(bytes.slice()),
        stderr: (bytes) => router.stderr(bytes.slice()),
        file: (parsed) => router.frame(ownedFrame(parsed)),
      });
      return { status, ...(await router.finish()) };
    };
    directory.open([{ name: 'dm2L_3k.fa', file: fasta }]);
    const words = ['windowmasker', '-mk_counts', '-in', 'dm2L_3k.fa', '-unit', '5', '-out', 'c.counts', '-logfile', 'log.txt'];
    expect(abi.validateTool(words)).toEqual({ ok: true });
    const made = await run(words);
    expect(made.status).toBe(0);
    expect(made.files.map((file) => file.name).sort()).toEqual(['c.counts', 'log.txt']);
    expect(made.files.find((file) => file.name === 'c.counts')!.file.size).toBeGreaterThan(0);
    // The same words again, on the same instance: the same bytes.
    const again = await run(words);
    expect(again.status).toBe(0);
    for (const name of ['c.counts', 'log.txt']) {
      expect(Buffer.from(await again.files.find((file) => file.name === name)!.file.arrayBuffer()).equals(Buffer.from(await made.files.find((file) => file.name === name)!.file.arrayBuffer())), name).toBe(true);
    }
    // The command line's refusal of a command line that names no mode: exit status 3, its text on the standard error.
    const mode = "one of '-mk_counts', '-convert' or '-ustat <stat_file>' must be specified";
    const refused = await run(['windowmasker', '-in', 'x.fa']);
    expect(refused.status).toBe(3);
    expect(await refused.stderr.text()).toContain(mode);
    const verdict = abi.validateTool(['windowmasker', '-in', 'x.fa']);
    expect(verdict.ok).toBe(false);
    expect(!verdict.ok && verdict.message).toContain(mode);
    // Another program is no tool: the engine's message, before any stream.
    expect(() => abi.tool(['blastp'], { stdout: () => undefined, stderr: () => undefined, file: () => undefined })).toThrow('blastp is not available in the browser (ABI v2)');
  }, 120_000);
});
