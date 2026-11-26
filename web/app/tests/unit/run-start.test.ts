// NCBI's application start in LOSAT Web (line RG item E1w-3b; docs/web/abi_v2.md §4 `losat_web2_start`,
// file kind 2, §5 stream 8, §9): the reader of the engine's archive of a run's files, the ABI binding
// of `start` and of a registry file held in memory, `validate` with the registry file of -conffile on
// the Data worker's instance, the Engine worker's runtime (the registry file for one search, the
// files of stream 8 apart from the outputs, with the result or the error), the gateway that hands
// them on, and the names of their downloads. The coordinator, the draft and the results have their
// cases in coordinator.test.ts, draft.test.ts and results.test.ts. With LOSAT_WEB_REACTORS, the same
// calls go to the built serial reactor.
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { describe, expect, it, vi } from 'vitest';
import { instantiateSerial, reactorControl, reopening, RunDirectory, WasmEngine } from 'losat';
import { readStoredZip, type RegistryFile, type RunFile } from 'losat/domain';
import { FileReadError, isFileReadError, runFilesOf, runOutputFailedOf, withRunFiles, type EngineRunRequest } from 'losat/ports';
import { findReactors } from '../../build/reactors';
import { outputFileName, runFileName } from '../../src/domain/output-format';
import { EngineRuntime } from '../../../packages/losat/src/infra/engine-worker/runtime';
import type { EngineCommand, EngineEvent, WorkerRun } from '../../../packages/losat/src/infra/engine-worker/protocol';
import { ABI_EXPORTS, inspectReactor } from '../../../packages/losat/src/infra/reactor/artifact';
import {
  discardedStream,
  discardSignal,
  DISCARD_STREAM,
  EngineCallError,
  EngineStoppedError,
  FILE_CHUNK_BYTES,
  FILE_REGISTRY,
  loadFile,
  MESSAGES_STREAM,
  messagesSignal,
  PLACE_NOWHERE,
  PLACE_RUN,
  PRINTED_STREAM,
  ReactorAbi,
  ROLE_QUERY,
  ROLE_SUBJECT,
  RUN_FILES_STREAM,
  type AbiExports,
  type FileLoader,
} from '../../../packages/losat/src/infra/reactor/abi';
import type { MainLinks } from '../../../packages/losat/src/infra/engine-worker/link';

/** The runtime's links to the main thread (link.ts): these tests start no threaded instance. */
const NO_LINKS: MainLinks = {
  link: () => {
    throw new Error('unused');
  },
  fault: () => undefined,
};

const encode = (text: string) => new TextEncoder().encode(text);

// --- the archive of a run's files ----------------------------------------------------------------

/** The CRC-32 of ZIP, written here apart from the reader's. */
function crc(bytes: Uint8Array): number {
  let value = ~0;
  for (const byte of bytes) {
    value ^= byte;
    for (let k = 0; k < 8; k++) value = value & 1 ? 0xedb88320 ^ (value >>> 1) : value >>> 1;
  }
  return ~value >>> 0;
}

/** An archive laid out as the adapter writes it (web/adapter/src/zip.rs): stored entries, the central directory, the end record. */
function storedZip(entries: ReadonlyArray<readonly [string, Uint8Array]>): Uint8Array {
  const parts: number[] = [];
  const u16 = (value: number) => parts.push(value & 0xff, (value >>> 8) & 0xff);
  const u32 = (value: number) => parts.push(value & 0xff, (value >>> 8) & 0xff, (value >>> 16) & 0xff, (value >>> 24) & 0xff);
  const central: Array<{ name: Uint8Array; crc: number; size: number; offset: number }> = [];
  for (const [text, bytes] of entries) {
    const name = encode(text);
    const offset = parts.length;
    u32(0x04034b50);
    [20, 0, 0, 0, 33].forEach(u16);
    u32(crc(bytes));
    u32(bytes.length);
    u32(bytes.length);
    u16(name.length);
    u16(0);
    parts.push(...name, ...bytes);
    central.push({ name, crc: crc(bytes), size: bytes.length, offset });
  }
  const start = parts.length;
  for (const entry of central) {
    u32(0x02014b50);
    [20, 20, 0, 0, 0, 33].forEach(u16);
    u32(entry.crc);
    u32(entry.size);
    u32(entry.size);
    [entry.name.length, 0, 0, 0, 0].forEach(u16);
    u32(0);
    u32(entry.offset);
    parts.push(...entry.name);
  }
  const size = parts.length - start;
  u32(0x06054b50);
  [0, 0, central.length, central.length].forEach(u16);
  u32(size);
  u32(start);
  u16(0);
  return new Uint8Array(parts);
}

describe("the reader of the engine's archive of a run's files", () => {
  it('reads the stored entries in order: their UTF-8 names and their bytes, an empty file too', () => {
    const files = readStoredZip(storedZip([['run.log', encode('Warning: one\n')], ['ü.log', new Uint8Array()]]));
    expect(files).toEqual([
      { name: 'run.log', bytes: encode('Warning: one\n') },
      { name: 'ü.log', bytes: new Uint8Array() },
    ]);
    expect(readStoredZip(storedZip([]))).toEqual([]);
  });

  it('refuses what the adapter does not write, naming the engine: a wrong CRC-32, a method, a data descriptor, a cut, a wrong directory, bytes after the end', () => {
    const good = storedZip([['run.log', encode('Warning: one\n')]]);
    const changed = (at: number, value: number) => {
      const bytes = good.slice();
      bytes[at] = value;
      return bytes;
    };
    const cases: Array<readonly [Uint8Array, RegExp]> = [
      // The first byte of the file's text.
      [changed(30 + 7, 0x58), /the bytes of entry 1 do not have its CRC-32/],
      // The method: deflate.
      [changed(8, 8), /entry 1 is not stored \(method 8\)/],
      // The flags: a data descriptor.
      [changed(6, 8), /entry 1 has the flags 0x8/],
      [good.subarray(0, good.length - 5), /it ends at byte|no end record/],
      [good.subarray(0, 40), /entry 1 ends past the archive/],
      // The central header's CRC-32.
      [changed(30 + 7 + 13 + 16, 0), /the central directory differs from entry 1/],
      [new Uint8Array([...good, 0]), /bytes follow the end record/],
      [encode('not an archive'), /it ends at byte|no end record/],
    ];
    for (const [archive, reason] of cases) {
      expect(() => readStoredZip(archive)).toThrow(/^the engine's archive of the run's files \(ABI v2 stream 8\) cannot be read: /);
      expect(() => readStoredZip(archive)).toThrow(reason);
    }
  });
});

// --- the ABI binding -------------------------------------------------------------------------------

/** Exports of a fake engine for `start`, `validate` and the file calls, which log what they are given. */
function fakeExports(
  options: { start?: string | 'fail'; validate?: string; validation?: string; trapValidate?: boolean; failEnd?: boolean } = {},
) {
  const memory = new WebAssembly.Memory({ initial: 4 });
  const log: string[] = [];
  let next = 4096;
  const ref: { abi?: ReactorAbi } = {};
  let error = '';
  const text = (ptr: number, len: number) => new TextDecoder().decode(new Uint8Array(memory.buffer, ptr, len));
  const respond = (json: string) => {
    const bytes = encode(json);
    new Uint8Array(memory.buffer, 60_000, bytes.length).set(bytes);
    ref.abi!.emit(2, 60_000, bytes.length);
  };
  const exports = {
    memory,
    losat_web2_alloc: (len: number) => {
      const ptr = next;
      next += len;
      return ptr;
    },
    losat_web2_dealloc: () => 0,
    losat_web2_last_error_ptr: () => {
      const bytes = encode(error);
      new Uint8Array(memory.buffer, 50_000, bytes.length).set(bytes);
      return 50_000;
    },
    losat_web2_last_error_len: () => encode(error).length,
    losat_web2_start: (ptr: number, len: number) => {
      log.push(`start ${text(ptr, len).split('\0').join(' ')}`);
      if (options.start === 'fail') {
        error = 'psiblast is not available in the browser (ABI v2)';
        return -1;
      }
      respond(options.start ?? '{"conffile":"c.ini","log_file":"run.log","dry_run":true,"register_inputs":false,"prints":null,"database":false}');
      return 0;
    },
    losat_web2_validate: (ptr: number, len: number) => {
      log.push(`validate ${text(ptr, len).split('\0').join(' ')}`);
      if (options.trapValidate) throw new WebAssembly.RuntimeError('unreachable');
      if (options.validate === undefined) {
        respond('{"data_loaders":true}');
        return 0;
      }
      // The *validation* that the engine sends before an option check refuses the argv (line WQA's D13 round 3, L1).
      if (options.validation !== undefined) respond(options.validation);
      error = options.validate;
      return -1;
    },
    losat_web2_file_begin: (kind: number, ptr: number, len: number, size: number) => {
      log.push(`begin ${kind} ${text(ptr, len)} ${size}`);
      return 9;
    },
    losat_web2_file_chunk: (file: number, ptr: number, len: number) => {
      log.push(`chunk ${file} ${text(ptr, len)}`);
      return 0;
    },
    losat_web2_file_end: (file: number) => {
      log.push(`end ${file}`);
      if (options.failEnd) {
        error = 'the input file ended after 1 of its 3 bytes';
        return -1;
      }
      respond('{"handle":9,"name":"c.ini","kind":"registry","size":8,"sha256":"ab"}');
      return 0;
    },
    losat_web2_file_release: (file: number) => {
      log.push(`release ${file}`);
      return 0;
    },
  } as unknown as AbiExports;
  ref.abi = new ReactorAbi(exports);
  return { abi: ref.abi, log };
}

describe('the ABI binding of the start and of a registry file', () => {
  it('reads the start JSON of an argv, or the CLI message of one that the engine cannot start', () => {
    expect(fakeExports().abi.start(['tblastx', '-dryrun'])).toEqual({ ok: true, start: { conffile: 'c.ini', logFile: 'run.log', logFileName: null, strategyFile: null, dryRun: true, registerInputs: false, prints: null, database: false } });
    // -parse_deflines: the run's inputs are read together (line WQA's sweep, D7, 2026-11-19).
    expect(fakeExports({ start: '{"conffile":null,"log_file":"logs/","log_file_name":"logs","dry_run":false,"register_inputs":true,"prints":null,"database":false}' }).abi.start(['blastn', '-parse_deflines', '-logfile', 'logs/'])).toEqual({
      ok: true,
      start: { conffile: null, logFile: 'logs/', logFileName: 'logs', strategyFile: null, dryRun: false, registerInputs: true, prints: null, database: false },
    });
    expect(fakeExports({ start: 'fail' }).abi.start(['psiblast'])).toEqual({ ok: false, message: 'psiblast is not available in the browser (ABI v2)' });
  });

  it('loads a file as kind 2: begin with its size, its bytes, end; a file that the engine refuses is released', async () => {
    const { abi, log } = fakeExports();
    expect(await loadFile(abi, FILE_REGISTRY, 'c.ini', new Blob(['[BLAST]\n']))).toEqual({
      handle: 9,
      name: 'c.ini',
      kind: 'registry',
      size: 8,
      sha256: 'ab',
    });
    expect(log).toEqual(['begin 2 c.ini 8', 'chunk 9 [BLAST]\n', 'end 9']);
    const failing = fakeExports({ failEnd: true });
    await expect(loadFile(failing.abi, FILE_REGISTRY, 'c.ini', new Blob(['[x\n']))).rejects.toThrow(EngineCallError);
    expect(failing.log).toEqual(['begin 2 c.ini 3', 'chunk 9 [x\n', 'end 9', 'release 9']);
  });

  // One loader for both kinds of file (the WindowMasker counts file's runtime cases are window-masker.test.ts's).
  it('gives a file larger than one chunk in chunks of FILE_CHUNK_BYTES through one buffer, each call through `call`', async () => {
    const calls: string[] = [];
    const loader = {
      stoppedBy: undefined,
      fileBegin: (kind: number, name: string, size: number) => (calls.push(`begin ${kind} ${name} ${size}`), 7),
      chunkBuffer: (capacity: number) => (calls.push(`buffer ${capacity}`), { ptr: 1, capacity }),
      fileChunk: (file: number, bytes: Uint8Array) => void calls.push(`chunk ${file} ${bytes.length} ${bytes[0]}`),
      fileEnd: (file: number) => (calls.push(`end ${file}`), { handle: file, name: 'c.ini', kind: 'registry' as const, size: 0, sha256: '' }),
      fileRelease: (file: number) => void calls.push(`release ${file}`),
      freeChunkBuffer: () => void calls.push('free buffer'),
    } as unknown as FileLoader;
    const bytes = new Uint8Array(FILE_CHUNK_BYTES + 3).fill(1);
    bytes[FILE_CHUNK_BYTES] = 2;
    let wrapped = 0;
    await loadFile(loader, FILE_REGISTRY, 'c.ini', new Blob([bytes]), (body) => (wrapped++, body()));
    expect(calls).toEqual([
      `begin 2 c.ini ${FILE_CHUNK_BYTES + 3}`,
      `buffer ${FILE_CHUNK_BYTES}`,
      `chunk 7 ${FILE_CHUNK_BYTES} 1`,
      'chunk 7 3 2',
      'end 7',
      'free buffer',
    ]);
    expect(wrapped).toBe(6);
  });

  it('requires losat_web2_start of an engine module', () => {
    expect(ABI_EXPORTS).toContain('losat_web2_start');
    // A module of every export but start (the module check's encoder is engine-runtime.test.ts's; a
    // module of exports only is enough for the check of the names).
    const names = ABI_EXPORTS.filter((name) => name !== 'losat_web2_start');
    expect(() => inspectReactor(exportsOnlyModule(names))).toThrow('the engine module lacks the exports losat_web2_start');
  });
});

/** A module with an own memory and a function `() -> i32` exported under each other name. */
function exportsOnlyModule(names: readonly string[]): Uint8Array {
  const uleb = (value: number): number[] => {
    const out: number[] = [];
    do {
      const byte = value & 0x7f;
      value = Math.floor(value / 128);
      out.push(byte | (value ? 0x80 : 0));
    } while (value);
    return out;
  };
  const name = (text: string) => [...uleb(encode(text).length), ...encode(text)];
  const vector = (items: number[][]) => [...uleb(items.length), ...items.flat()];
  const section = (id: number, bytes: number[]) => [id, ...uleb(bytes.length), ...bytes];
  const functions = names.filter((n) => n !== 'memory');
  return new Uint8Array([
    0, 0x61, 0x73, 0x6d, 1, 0, 0, 0,
    ...section(1, [1, 0x60, 0, 1, 0x7f]),
    ...section(3, vector(functions.map(() => [0]))),
    ...section(5, [1, 1, 17, ...uleb(16384)]),
    ...section(7, vector(names.map((n) => (n === 'memory' ? [...name(n), 2, 0] : [...name(n), 0, ...uleb(functions.indexOf(n))])))),
    ...section(10, vector(functions.map(() => [4, 0, 0x41, 0, 0x0b]))),
  ]);
}

/**
 * An instance for the control's tests that logs, under `label`, the files it is given and its
 * validations with the files that it holds when it validates.
 */
function fakeInstance(label: string, log: string[]): ReactorAbi {
  const held = new Set<string>();
  const loading = new Map<number, string>();
  let handle = 0;
  return {
    stoppedBy: undefined,
    fileBegin: (kind: number, name: string, size: number) => {
      log.push(`${label} begin ${kind} ${name} ${size}`);
      loading.set(++handle, name);
      return handle;
    },
    chunkBuffer: (capacity: number) => ({ ptr: 1, capacity }),
    fileChunk: () => undefined,
    fileEnd: (file: number) => {
      held.add(loading.get(file)!);
      return { handle: file, name: loading.get(file)!, kind: 'registry', size: 0, sha256: '' };
    },
    fileRelease: (file: number) => void held.delete(loading.get(file)!),
    freeChunkBuffer: () => undefined,
    validate: (argv: readonly string[]) => {
      log.push(`${label} validate ${argv.join(' ')} held=[${[...held].join(',')}]`);
      return { ok: true, dataLoaders: true };
    },
    describe: () => {
      log.push(`${label} describe held=[${[...held].join(',')}]`);
      return { program: 'blastn', formats: [0, 6, 7], parameters: [] };
    },
  } as unknown as ReactorAbi;
}

const pause = (ms: number) => new Promise((resolve) => setTimeout(resolve, ms));

// DECISIONS E1w3b-4: a validation with the registry file of -conffile runs on an instance of its own,
// opened for it and dropped after it; the Data worker's shared instance never holds the file.
describe("the Data worker's validate with the registry file of -conffile", () => {
  it('validates on an instance of its own, opened for that validation; the shared instance is not called', async () => {
    const log: string[] = [];
    const shared = fakeInstance('shared', log);
    let opened = 0;
    const control = reactorControl(async () => shared, async () => fakeInstance(`own${++opened}`, log));
    expect(await control.validate(['blastp', '-conffile', 'c.ini'], { name: 'c.ini', file: new Blob(['[BLAST]\n']) })).toEqual({ ok: true, dataLoaders: true });
    expect(await control.validate(['blastp', '-conffile', 'c.ini'], { name: 'c.ini', file: new Blob(['[BLAST]\n']) })).toEqual({ ok: true, dataLoaders: true });
    // Without a file, the shared instance validates.
    await control.validate(['blastp']);
    expect(log).toEqual([
      'own1 begin 2 c.ini 8',
      'own1 validate blastp -conffile c.ini held=[c.ini]',
      'own2 begin 2 c.ini 8',
      'own2 validate blastp -conffile c.ini held=[c.ini]',
      'shared validate blastp held=[]',
    ]);
    expect(opened).toBe(2);
    // The verdict and start are the engine's.
    const { abi } = fakeExports({ validate: 'LOSAT query/options error: bad' });
    const engine = reactorControl(async () => abi, async () => fakeExports({ validate: 'LOSAT query/options error: bad' }).abi);
    expect(await engine.validate(['blastp', '-conffile', 'c.ini'], { name: 'c.ini', file: new Blob(['[BLAST]\n']) })).toEqual({
      ok: false,
      message: 'LOSAT query/options error: bad',
    });
    expect(await engine.start(['blastp', '-dryrun'])).toEqual({ ok: true, start: { conffile: 'c.ini', logFile: 'run.log', logFileName: null, strategyFile: null, dryRun: true, registerInputs: false, prints: null, database: false } });
  });

  // line WQA's D13 round 3, L1, 2026-11-20: the engine sends the data loaders of the argv's registry before its option
  // checks; a refused argv gives them with its message, which the host reads after them.
  it("gives the data loaders of a refused argv's registry with its message", async () => {
    const refused = fakeExports({ validate: 'Unexpected extra argument: -x', validation: '{"data_loaders":false}' });
    const engine = reactorControl(async () => refused.abi, async () => refused.abi);
    expect(await engine.validate(['blastp', '-x'])).toEqual({ ok: false, message: 'Unexpected extra argument: -x', dataLoaders: false });
    // An argv whose start failed before the registry has none.
    const failed = fakeExports({ validate: 'BLAST query/options error: bad' });
    expect(await reactorControl(async () => failed.abi, async () => failed.abi).validate(['blastp'])).toEqual({
      ok: false,
      message: 'BLAST query/options error: bad',
    });
  });

  it('keeps the file away from every call of the shared instance, also from calls that wait with it for one opening of that instance', async () => {
    const log: string[] = [];
    const shared = fakeInstance('shared', log);
    // The shared instance opens (or reopens after a trap) while every call waits for it (scanner.ts, `reopening`).
    const opening = reopening(async () => {
      await pause(5);
      return shared;
    });
    let opened = 0;
    const control = reactorControl(opening, async () => {
      const label = `own${++opened}`;
      await pause(1);
      return fakeInstance(label, log);
    });
    const file = new Blob(['[BLAST]\n']);
    const results = await Promise.all([
      control.validate(['blastn', '-conffile', 'settings.ini'], { name: 'settings.ini', file }),
      control.validate(['blastn', '-conffile', 'settings.ini'], { name: 'settings.ini', file }),
      control.validate(['blastn']),
      control.describe('blastn'),
      control.validate(['blastn', '-conffile', '.ncbirc'], { name: '.ncbirc', file }),
      control.validate(['blastn']),
    ]);
    expect(results.filter((result) => 'ok' in result)).toEqual(Array(5).fill({ ok: true, dataLoaders: true }));
    const of = (label: string) => log.filter((line) => line.startsWith(`${label} `)).sort();
    expect(of('shared')).toEqual(['shared describe held=[]', 'shared validate blastn held=[]', 'shared validate blastn held=[]']);
    expect([of('own1'), of('own2'), of('own3')].flat().sort()).toEqual(
      [
        'own1 begin 2 settings.ini 8',
        'own1 validate blastn -conffile settings.ini held=[settings.ini]',
        'own2 begin 2 settings.ini 8',
        'own2 validate blastn -conffile settings.ini held=[settings.ini]',
        'own3 begin 2 .ncbirc 8',
        'own3 validate blastn -conffile .ncbirc held=[.ncbirc]',
      ].sort(),
    );
    // Each own instance held its file alone.
    for (const label of ['own1', 'own2', 'own3']) expect(of(label).filter((line) => line.includes(' begin '))).toHaveLength(1);
  });

  it('fails the validation when its instance cannot be opened, the file is refused or the instance stops, and leaves the shared instance alone', async () => {
    const log: string[] = [];
    const shared = fakeInstance('shared', log);
    const registry = { name: 'c.ini', file: new Blob(['[BLAST]\n']) };
    await expect(
      reactorControl(
        async () => shared,
        async () => {
          throw new Error('out of memory');
        },
      ).validate(['blastp', '-conffile', 'c.ini'], registry),
    ).rejects.toThrow('the engine could not open an instance for the file of -conffile: out of memory');
    const refused = fakeExports({ failEnd: true });
    await expect(reactorControl(async () => shared, async () => refused.abi).validate(['blastp', '-conffile', 'c.ini'], { name: 'c.ini', file: new Blob(['[x\n']) })).rejects.toThrow(
      'the input file ended after 1 of its 3 bytes',
    );
    expect(refused.log).toEqual(['begin 2 c.ini 3', 'chunk 9 [x\n', 'end 9', 'release 9']);
    const trapped = fakeExports({ trapValidate: true });
    await expect(reactorControl(async () => shared, async () => trapped.abi).validate(['blastp', '-conffile', 'c.ini'], registry)).rejects.toThrow(EngineStoppedError);
    expect(trapped.log).toEqual(['begin 2 c.ini 8', 'chunk 9 [BLAST]\n', 'end 9', 'validate blastp -conffile c.ini']);
    // A file that the browser no longer reads (it changed after it was chosen): a FileReadError with the
    // browser's reason, the file released, no validation.
    const reason = 'The requested file could not be read, typically due to permission problems that have occurred after a reference to a file was acquired.';
    const unreadable = { size: 8, slice: () => ({ arrayBuffer: () => Promise.reject(new DOMException(reason, 'NotReadableError')) }) } as unknown as Blob;
    const own = fakeExports();
    const error = await reactorControl(async () => shared, async () => own.abi)
      .validate(['blastp', '-conffile', 'c.ini'], { name: 'c.ini', file: unreadable })
      .catch((e: unknown) => e);
    expect(error).toBeInstanceOf(FileReadError);
    expect(isFileReadError(error)).toBe(true);
    expect((error as Error).message).toBe(reason);
    expect(own.log).toEqual(['begin 2 c.ini 8', 'release 9']);
    expect(log).toEqual([]);
  });
});

// --- the Engine worker's runtime -------------------------------------------------------------------

/**
 * A binding with no engine behind it that records the calls of a search; its run sends outfmt 6 and
 * the archive of `runFiles` on stream 8 (in two chunks), and then fails with `fail`, or stops.
 * `heldName` is the name of a file that the instance holds (a WindowMasker counts file): `file_begin`
 * refuses it, as the adapter does. `discard`: after outfmt 6, the signal of stream 9 that discards
 * it (its bytes), then outfmt 7 (ABI v2 §5).
 */
/** The signal of stream 10 (ABI v2 §5): the place of the run's messages, 4 bytes (line WQA's sweep, D7, 2026-11-20). */
function placeSignal(place: number): Uint8Array {
  const bytes = new Uint8Array(4);
  new DataView(bytes.buffer).setUint32(0, place, true);
  return bytes;
}

/** The value of `option` in `argv`, as the engine's *start* gives it (null when the argv has none). */
function valueOf(argv: readonly string[], option: string): string | null {
  const at = argv.indexOf(option);
  return at < 0 || at + 1 >= argv.length ? null : argv[at + 1]!;
}

function fakeAbi(
  options: {
    runFiles?: Uint8Array;
    fail?: string;
    stop?: boolean;
    heldName?: string;
    discard?: Uint8Array;
    warnings?: readonly Uint8Array[];
    messages?: Uint8Array | 'none';
    /** The run prints this text in place of a search (its start says `prints: 'help'`). */
    printed?: Uint8Array;
    /** Normalized host file name of `-logfile` (`StartInfo.logFileName`), when different from the argv value. */
    logFileName?: string | null;
    /** The data loaders that `validate` finds in the argv's registry (on when not given). */
    dataLoaders?: () => boolean;
  } = {},
) {
  const calls: string[] = [];
  let handle = 100;
  let stopped: EngineStoppedError | undefined;
  const abi = {
    get stoppedBy() {
      return stopped;
    },
    fileBegin: (kind: number, name: string, size: number) => {
      calls.push(`begin ${kind} ${name} ${size}`);
      if (name === options.heldName) throw new EngineCallError(`a file named ${name} is registered: release it before registering another`);
      return ++handle;
    },
    chunkBuffer: (capacity: number) => (calls.push(`buffer ${capacity}`), { ptr: 1, capacity }),
    fileChunk: (file: number, bytes: Uint8Array) => void calls.push(`chunk ${file} ${new TextDecoder().decode(bytes)}`),
    fileEnd: (file: number) => (calls.push(`end ${file}`), { handle: file, name: 'c.ini', kind: 'registry', size: 8, sha256: '0' }),
    freeChunkBuffer: () => void calls.push('free buffer'),
    validate: (argv: readonly string[]) => {
      calls.push(`validate ${argv.join(' ')}`);
      return { ok: true, dataLoaders: options.dataLoaders?.() ?? true };
    },
    // The start of an argv: -parse_deflines reads the inputs of the run together (line WQA's sweep, D7, 2026-11-19); the
    // names of -logfile and -export_search_strategy name the files of stream 8 (line RG item F1w).
    start: (argv: readonly string[]) => ({
      ok: true,
      start: {
        conffile: null,
        logFile: valueOf(argv, '-logfile'),
        logFileName: options.logFileName !== undefined ? options.logFileName : valueOf(argv, '-logfile'),
        strategyFile: valueOf(argv, '-export_search_strategy'),
        dryRun: false,
        registerInputs: argv.includes('-parse_deflines'),
        prints: options.printed === undefined ? null : 'help',
        database: argv.includes('-db'),
      },
    }),
    memoryBytes: () => 1,
    register: (_program: string, role: number, _dataLoaders: boolean, bytes: Uint8Array) => {
      calls.push(`register ${role === ROLE_QUERY ? 'query' : 'subject'}`);
      return { handle: ++handle, records: [{ id: role === ROLE_SUBJECT ? 's' : 'q', length: bytes.length }] };
    },
    registerInputs: (argv: readonly string[], query: Uint8Array, subject: Uint8Array) => {
      calls.push(`register inputs ${argv.join(' ')}`);
      // A run of -db reads the query alone (ABI v2 §9 *inputs*, `subject` null).
      const database = argv.includes('-db');
      return {
        query: { handle: ++handle, records: [{ id: 'q', length: query.length }] },
        subject: database ? undefined : { handle: ++handle, records: [{ id: 's', length: subject.length }] },
        shownIds: { query: ['Q.1'], subject: database ? [] : ['S.1'] },
      };
    },
    release: (h: number) => void calls.push(`release ${h}`),
    fileRelease: (file: number) => void calls.push(`file release ${file}`),
    run: (_argv: readonly string[], query: number, subject: number, sink: (stream: number, bytes: Uint8Array) => void) => {
      // A search of -db takes subject handle 0; a run that prints takes the handles 0 and says so below.
      calls.push(subject === 0 && options.printed === undefined ? 'run without a subject' : 'run');
      if (options.messages !== 'none') sink(MESSAGES_STREAM, options.messages ?? placeSignal(3));
      if (options.printed !== undefined) {
        calls.push(`run handles ${query} ${subject}`);
        const half = options.printed.length >> 1;
        sink(PRINTED_STREAM, options.printed.subarray(0, half));
        sink(PRINTED_STREAM, options.printed.subarray(half));
        return;
      }
      sink(6, encode('q\ts\n'));
      for (const chunk of options.warnings ?? []) sink(3, chunk);
      if (options.discard !== undefined) {
        sink(DISCARD_STREAM, options.discard);
        sink(7, encode('# 0 hits found\n'));
      }
      if (options.runFiles !== undefined) {
        const half = options.runFiles.length >> 1;
        sink(RUN_FILES_STREAM, options.runFiles.subarray(0, half));
        sink(RUN_FILES_STREAM, options.runFiles.subarray(half));
      }
      if (options.stop) {
        stopped = new EngineStoppedError('The engine stopped: unreachable');
        throw stopped;
      }
      if (options.fail !== undefined) throw new EngineCallError(options.fail);
    },
  };
  return { abi: abi as unknown as ReactorAbi, calls };
}

function runtimeWith(abi: ReactorAbi) {
  const runtime = new EngineRuntime({
    serial: async () => {
      throw new Error('unused');
    },
    threads: async () => {
      throw new Error('unused');
    },
    builds: { serial: 'serial', threads: 'threads' },
    use: undefined,
    main: NO_LINKS,
  });
  // The runtime holds this instance as one that a search loaded.
  (runtime as unknown as { loaded: unknown }).loaded = {
    path: 'serial',
    instance: { abi, wasi: { clear: () => undefined } },
    directory: new RunDirectory(),
    runs: 0,
    retained: undefined,
    counts: undefined,
    ran: new Set(),
  };
  return runtime;
}

/** A run of `words` after the inputs, whose output port records what the Engine worker sent. */
function workerRun(words: readonly string[], conffile?: RegistryFile) {
  const channel = new MessageChannel();
  const received: Array<{ readonly type: string }> = [];
  channel.port2.onmessage = (event: MessageEvent<{ type: string }>) => received.push(event.data);
  const input = (role: string) => ({ bytes: new Uint8Array([62, 10, 65]), records: [{ id: role, length: 3 }], key: `k-${role}` });
  const run: WorkerRun = {
    type: 'run',
    generation: 1,
    runId: 'r1',
    program: 'tblastx',
    argv: ['tblastx', '-query', 'q.fa', '-subject', 's.fa', ...words],
    threads: 1,
    query: input('q'),
    subject: input('s'),
    ...(conffile === undefined ? {} : { conffile }),
    output: channel.port1,
  };
  /**
   * What the run sent on its port, all of it: a port delivers in order, so once a mark posted after
   * the run has arrived, every message that the run posted has.
   */
  const delivered = async () => {
    channel.port1.postMessage({ type: 'mark' });
    await expect.poll(() => received.some((message) => message.type === 'mark')).toBe(true);
    return received.filter((message) => message.type !== 'mark');
  };
  const close = () => {
    channel.port1.close();
    channel.port2.close();
  };
  return { run, delivered, close };
}

const LOG: RunFile = { name: 'run.log', bytes: encode('LOSAT engine error: The average subject length is too short\n') };
const CONFFILE = { name: 'c.ini', file: new Blob(['[BLAST]\n']) };

describe("the Engine worker's runtime and the files of a run", () => {
  it('registers the file of -conffile before validate, keeps stream 8 apart from the outputs, and releases the file before the end', async () => {
    const { abi, calls } = fakeAbi({ runFiles: storedZip([[LOG.name, LOG.bytes]]) });
    const runtime = runtimeWith(abi);
    const { run, delivered, close } = workerRun(['-conffile', 'c.ini', '-logfile', 'run.log'], CONFFILE);
    const result = await runtime.run(run, () => undefined);
    expect(result.runFiles).toEqual([LOG]);
    expect(calls).toEqual([
      'begin 2 c.ini 8',
      'buffer 8',
      'chunk 101 [BLAST]\n',
      'end 101',
      'free buffer',
      'validate tblastx -query q.fa -subject s.fa -conffile c.ini -logfile run.log -num_threads 1',
      'register query',
      'register subject',
      'run',
      'release 102',
      'file release 101',
    ]);
    // The Data worker gets the outputs only: one chunk of outfmt 6 and the end.
    expect(await delivered()).toEqual([{ type: 'chunk', stream: 6, bytes: encode('q\ts\n') }, { type: 'end', chunks: 1, bytes: 4 }]);
    close();
    // A search without the option registers nothing, and one without files gives none.
    const plain = fakeAbi();
    const next = workerRun([]);
    expect('runFiles' in (await runtimeWith(plain.abi).run(next.run, () => undefined))).toBe(false);
    expect(plain.calls.some((call) => call.startsWith('begin') || call.startsWith('file release'))).toBe(false);
    next.close();
  });

  // line RG item F1w: the file of -export_search_strategy is a file of stream 8 too; the *start* of the argv names
  // which option made each file, and a file that the argv does not name is the engine's error.
  it('names each file of stream 8 by the option of the argv that names it', async () => {
    const strategy = encode('Blast4-request ::= {\n}\n');
    const both = fakeAbi({ runFiles: storedZip([['out/s.asn', strategy], [LOG.name, LOG.bytes]]) });
    const run = workerRun(['-export_search_strategy', 'out/s.asn', '-logfile', 'run.log']);
    expect((await runtimeWith(both.abi).run(run.run, () => undefined)).runFiles).toEqual([
      { name: 'out/s.asn', bytes: strategy, option: 'export_search_strategy' },
      LOG,
    ]);
    run.close();
    const failed = fakeAbi({ runFiles: storedZip([['out/s.asn', strategy]]), fail: 'LOSAT engine error: x' });
    const second = workerRun(['-export_search_strategy', 'out/s.asn']);
    const error = await runtimeWith(failed.abi)
      .run(second.run, () => undefined)
      .catch((e: unknown) => e);
    expect(runFilesOf(error)).toEqual([{ name: 'out/s.asn', bytes: strategy, option: 'export_search_strategy' }]);
    second.close();
    const unnamed = fakeAbi({ runFiles: storedZip([['other.asn', strategy]]) });
    const third = workerRun(['-export_search_strategy', 'out/s.asn']);
    await expect(runtimeWith(unnamed.abi).run(third.run, () => undefined)).rejects.toThrow('has a file that the argv does not name: "other.asn"');
    third.close();
    // Normalized -logfile host names (logs/ -> logs, ./run.log -> run.log) match `logFileName`, and
    // win over `strategyFile` when both would give the same normalized name.
    const normalized = fakeAbi({ logFileName: 'logs', runFiles: storedZip([['out/s.asn', strategy], ['logs', LOG.bytes]]) });
    const fourth = workerRun(['-export_search_strategy', 'out/s.asn', '-logfile', 'logs/']);
    expect((await runtimeWith(normalized.abi).run(fourth.run, () => undefined)).runFiles).toEqual([
      { name: 'out/s.asn', bytes: strategy, option: 'export_search_strategy' },
      { name: 'logs', bytes: LOG.bytes },
    ]);
    fourth.close();
    const sameName = fakeAbi({ logFileName: 'logs', runFiles: storedZip([['logs', LOG.bytes]]) });
    const fifth = workerRun(['-export_search_strategy', 'logs', '-logfile', 'logs/']);
    expect((await runtimeWith(sameName.abi).run(fifth.run, () => undefined)).runFiles).toEqual([
      { name: 'logs', bytes: LOG.bytes },
    ]);
    fifth.close();
  });

  // line WQA's sweep, D7, 2026-11-19: a run of -parse_deflines reads its inputs together for its argv, releases the
  // subject that the instance retained, and releases both handles after the run.
  it('reads the inputs of a -parse_deflines run together for its argv, and keeps no subject', async () => {
    const { abi, calls } = fakeAbi();
    const runtime = runtimeWith(abi);
    const first = workerRun([]);
    expect((await runtime.run(first.run, () => undefined)).subjectRetained).toBe(false);
    first.close();
    const parsed = workerRun(['-parse_deflines']);
    const result = await runtime.run(parsed.run, () => undefined);
    expect(result.subjectRetained).toBe(false);
    expect(result.shownIds).toEqual({ query: ['Q.1'], subject: ['S.1'] });
    parsed.close();
    expect(calls).toEqual([
      'validate tblastx -query q.fa -subject s.fa -num_threads 1',
      'register query',
      'register subject',
      'run',
      'release 101',
      'validate tblastx -query q.fa -subject s.fa -parse_deflines -num_threads 1',
      'release 102',
      'register inputs tblastx -query q.fa -subject s.fa -parse_deflines -num_threads 1',
      'run',
      'release 103',
      'release 104',
    ]);
    // The next run without the option registers its subject again.
    const again = workerRun([]);
    expect((await runtime.run(again.run, () => undefined)).subjectRetained).toBe(false);
    again.close();
    expect(calls.slice(-5)).toEqual(['validate tblastx -query q.fa -subject s.fa -num_threads 1', 'register query', 'register subject', 'run', 'release 105']);
  });

  // line WQA's D13 round 2, L1, 2026-11-20: the instance keeps a subject for the next run of the same subject only when
  // it was read as that run reads it. Read with the data loaders off, it is not what a run with them on reads (that
  // reader refuses a first line that is not a definition line), so the run registers it again; read with them on, it
  // serves a run without them.
  it('serves a retained subject read with the data loaders off only to runs without them', async () => {
    let loaders = false;
    const { abi, calls } = fakeAbi({ dataLoaders: () => loaders });
    const runtime = runtimeWith(abi);
    const once = async () => {
      const one = workerRun([]);
      const { subjectRetained } = await runtime.run(one.run, () => undefined);
      one.close();
      return subjectRetained;
    };
    const subjects = () => calls.filter((call) => call === 'register subject').length;
    expect(await once()).toBe(false);
    expect(await once()).toBe(true);
    expect(subjects()).toBe(1);
    expect(calls).not.toContain('release 102');
    loaders = true;
    expect(await once()).toBe(false);
    expect(subjects()).toBe(2);
    // The subject read with the loaders off (handle 102) was let go before the new one was read.
    expect(calls.slice(calls.lastIndexOf('register query'))).toEqual(['register query', 'release 102', 'register subject', 'run', 'release 104']);
    loaders = false;
    expect(await once()).toBe(true);
    expect(subjects()).toBe(2);
  });

  it("forwards the engine's discard of an output (stream 9) to the Data worker with where its error is, between the chunks in their order", async () => {
    // The signal of 8 bytes: the stream, then the place of its error (stream 3, stream 8 for the file of -logfile, or
    // 0xffffffff for nowhere; line WQA's sweep, D7, 2026-11-20).
    const signal = (stream: number, place: number) => {
      const bytes = new Uint8Array(8);
      new DataView(bytes.buffer).setUint32(0, stream, true);
      new DataView(bytes.buffer).setUint32(4, place, true);
      return bytes;
    };
    for (const [place, error] of [
      [3, 'diagnostics'],
      [8, 'logfile'],
      [PLACE_NOWHERE, 'nowhere'],
      // The error that ended the run: an output that a failed run does not keep (ABI version 7; line WQA's sweep, D15, 2026-11-20).
      [PLACE_RUN, 'run'],
    ] as const) {
      const { abi } = fakeAbi({ discard: signal(6, place) });
      const { run, delivered, close } = workerRun([]);
      await runtimeWith(abi).run(run, () => undefined);
      expect(await delivered()).toEqual([
        { type: 'chunk', stream: 6, bytes: encode('q\ts\n') },
        { type: 'discard', stream: 6, error },
        { type: 'chunk', stream: 7, bytes: encode('# 0 hits found\n') },
        { type: 'end', chunks: 2, bytes: 19 },
      ]);
      close();
    }
    expect(discardedStream(signal(7, 3))).toBe(7);
    expect(discardSignal(signal(7, 8))).toEqual({ stream: 7, errorStream: 8, ended: false });
    expect(discardSignal(signal(0, PLACE_NOWHERE))).toEqual({ stream: 0, errorStream: undefined, ended: false });
    expect(discardSignal(signal(6, PLACE_RUN))).toEqual({ stream: 6, errorStream: undefined, ended: true });
    expect(() => discardedStream(new Uint8Array([6, 0, 0, 0]))).toThrow('a signal of stream 9 has 4 bytes, not 8');
    // A signal that names a stream which is not an output of a run, or a place that holds no error, breaks the ABI:
    // the run fails.
    for (const [bytes, error] of [
      [signal(1, 3), 'which is not an output of a run'],
      [signal(0x1000006, 3), 'which is not an output of a run'],
      [signal(6, 6), 'on stream 6, which holds no error'],
    ] as const) {
      const broken = workerRun([]);
      await expect(runtimeWith(fakeAbi({ discard: bytes }).abi).run(broken.run, () => undefined)).rejects.toThrow(error);
      broken.close();
    }
  });

  // Line WQA's D14, 2026-11-20: a run whose start says that it prints the command line's help or version reads no input
  // (handles 0) and gives the text of stream 11 with its result, joined, and no chunk to the Data worker.
  it("runs an argv that prints the command line's text without inputs, and gives the text with its result", async () => {
    const help = encode('USAGE\n  losat tblastx [OPTIONS]\n');
    const { abi, calls } = fakeAbi({ printed: help });
    const { run, delivered, close } = workerRun(['-help']);
    const result = await runtimeWith(abi).run(run, () => undefined);
    expect(result.printed).toEqual(help);
    expect(result.messages).toBe('diagnostics');
    expect(calls).toEqual(['validate tblastx -query q.fa -subject s.fa -help -num_threads 1', 'run', 'run handles 0 0']);
    expect(await delivered()).toEqual([{ type: 'end', chunks: 0, bytes: 0 }]);
    close();
    // A search gives no printed text.
    const search = workerRun([]);
    expect('printed' in (await runtimeWith(fakeAbi().abi).run(search.run, () => undefined))).toBe(false);
    search.close();
  });

  // line WQA's sweep, D7, 2026-11-20: where the run's messages went is the engine's word (its signal of stream 10), which
  // the run gives with its result, not a stream for the Data worker; a run without it, or with a place that holds no
  // messages, breaks the ABI.
  it("gives where the run's messages went (stream 10) with its result, and nothing of it to the Data worker", async () => {
    for (const [place, messages] of [
      [3, 'diagnostics'],
      [8, 'logfile'],
      [PLACE_NOWHERE, 'nowhere'],
    ] as const) {
      const { abi } = fakeAbi({ messages: placeSignal(place) });
      const { run, delivered, close } = workerRun([]);
      expect((await runtimeWith(abi).run(run, () => undefined)).messages).toBe(messages);
      expect(await delivered()).toEqual([{ type: 'chunk', stream: 6, bytes: encode('q\ts\n') }, { type: 'end', chunks: 1, bytes: 4 }]);
      close();
    }
    expect(messagesSignal(placeSignal(3))).toBe(3);
    expect(messagesSignal(placeSignal(8))).toBe(8);
    expect(messagesSignal(placeSignal(PLACE_NOWHERE))).toBeUndefined();
    expect(() => messagesSignal(new Uint8Array(8))).toThrow('a signal of stream 10 has 8 bytes, not 4');
    for (const [messages, error] of [
      ['none', 'the engine did not say where the messages of the run go (stream 10)'],
      [placeSignal(6), 'the engine put the messages of the run on stream 6, which holds no messages'],
      [new Uint8Array(3), 'a signal of stream 10 has 3 bytes, not 4'],
    ] as const) {
      const broken = workerRun([]);
      await expect(runtimeWith(fakeAbi({ messages }).abi).run(broken.run, () => undefined)).rejects.toThrow(error);
      broken.close();
    }
  });

  it("gives a failed run's files with its error, releases what it registered, and sends fail, not end", async () => {
    const { abi, calls } = fakeAbi({ runFiles: storedZip([[LOG.name, LOG.bytes]]), fail: 'LOSAT engine error: The average subject length is too short' });
    const { run, delivered, close } = workerRun(['-logfile', 'run.log'], CONFFILE);
    const error = await runtimeWith(abi)
      .run(run, () => undefined)
      .catch((e: unknown) => e);
    expect(error).toBeInstanceOf(EngineCallError);
    expect((error as Error).message).toBe('LOSAT engine error: The average subject length is too short');
    expect(runFilesOf(error)).toEqual([LOG]);
    expect(calls.slice(-2)).toEqual(['release 102', 'file release 101']);
    // The chunk written before the failure, then fail, and no end.
    expect(await delivered()).toEqual([
      { type: 'chunk', stream: 6, bytes: encode('q\ts\n') },
      { type: 'fail', chunks: 1, bytes: 4 },
    ]);
    close();
  });

  // Line WQA's sweep, D3 (2026-11-19): the warnings that the command line writes before its error
  // (stream 3) reach the Data worker, once, and the output of the failed run ends with `fail` before
  // the run's rejection, which says so; the Data worker reads them from what arrived. The runtime keeps
  // no copy of them. A stop sends no `fail`, and its rejection does not say so.
  it("ends a failed run's output with fail after its warnings (stream 3), and says so with the error", async () => {
    const text = encode('Warning: [blastn] Subject_1 \u00e9 empty: Subject sequence contains no data\n');
    // Two chunks that split the two bytes of é.
    const cut = text.indexOf(0xc3) + 1;
    const warnings = [text.subarray(0, cut), text.subarray(cut)];
    const { abi } = fakeAbi({ warnings, fail: 'LOSAT engine error: The average subject length is too short' });
    const { run, delivered, close } = workerRun([]);
    const error = await runtimeWith(abi)
      .run(run, () => undefined)
      .catch((e: unknown) => e);
    expect((error as Error).message).toBe('LOSAT engine error: The average subject length is too short');
    expect(runOutputFailedOf(error)).toBe(true);
    expect('diagnostics' in (error as Error)).toBe(false);
    expect(runFilesOf(error)).toBeUndefined();
    expect(await delivered()).toEqual([
      { type: 'chunk', stream: 6, bytes: encode('q\ts\n') },
      { type: 'chunk', stream: 3, bytes: warnings[0] },
      { type: 'chunk', stream: 3, bytes: warnings[1] },
      { type: 'fail', chunks: 3, bytes: 4 + text.length },
    ]);
    close();
    // A failed run without warnings ends with `fail` too; a stop sends nothing after its chunks.
    const silent = workerRun([]);
    const failed = await runtimeWith(fakeAbi({ fail: 'LOSAT engine error: Empty CBlastQueryVector' }).abi)
      .run(silent.run, () => undefined)
      .catch((e: unknown) => e);
    expect(runOutputFailedOf(failed)).toBe(true);
    expect(await silent.delivered()).toEqual([
      { type: 'chunk', stream: 6, bytes: encode('q\ts\n') },
      { type: 'fail', chunks: 1, bytes: 4 },
    ]);
    silent.close();
    const stopped = workerRun([]);
    const stop = await runtimeWith(fakeAbi({ warnings, stop: true }).abi)
      .run(stopped.run, () => undefined)
      .catch((e: unknown) => e);
    expect(stop).toBeInstanceOf(Error);
    expect(runOutputFailedOf(stop)).toBe(false);
    expect((await stopped.delivered()).some((message) => message.type === 'fail' || message.type === 'end')).toBe(false);
    stopped.close();
  });

  it('fails a run whose file of -conffile the engine refuses, with nothing registered and nothing searched', async () => {
    // `-window_masker_db x -conffile x` on an instance that holds the counts file x: the adapter refuses the name.
    const { abi, calls } = fakeAbi({ heldName: 'c.ini' });
    const { run, delivered, close } = workerRun(['-conffile', 'c.ini'], CONFFILE);
    await expect(runtimeWith(abi).run(run, () => undefined)).rejects.toThrow('a file named c.ini is registered: release it before registering another');
    expect(calls).toEqual(['begin 2 c.ini 8']);
    expect(await delivered()).toEqual([]);
    close();
  });

  it('gives no files of a run that stopped, and calls the stopped instance no more', async () => {
    const { abi, calls } = fakeAbi({ runFiles: storedZip([[LOG.name, LOG.bytes]]), stop: true });
    const { run, close } = workerRun(['-logfile', 'run.log'], CONFFILE);
    const error = await runtimeWith(abi)
      .run(run, () => undefined)
      .catch((e: unknown) => e);
    expect(error).toBeInstanceOf(EngineStoppedError);
    expect(runFilesOf(error)).toBeUndefined();
    expect(calls.at(-1)).toBe('run');
    close();
  });

  it("fails a run whose archive of files cannot be read, before its end; a failed run's error then tells both", async () => {
    const broken = storedZip([[LOG.name, LOG.bytes]]);
    broken[8] = 8;
    const completed = fakeAbi({ runFiles: broken });
    const first = workerRun(['-logfile', 'run.log']);
    await expect(runtimeWith(completed.abi).run(first.run, () => undefined)).rejects.toThrow('entry 1 is not stored (method 8)');
    expect(await first.delivered()).toEqual([{ type: 'chunk', stream: 6, bytes: encode('q\ts\n') }]);
    first.close();
    const failed = fakeAbi({ runFiles: broken, fail: 'LOSAT engine error: x' });
    const second = workerRun(['-logfile', 'run.log']);
    const error = await runtimeWith(failed.abi)
      .run(second.run, () => undefined)
      .catch((e: unknown) => e);
    expect((error as Error).message).toMatch(/^LOSAT engine error: x \(the engine's archive of the run's files \(ABI v2 stream 8\) cannot be read: entry 1 is not stored/);
    expect(runFilesOf(error)).toBeUndefined();
    second.close();
  });

  it("runs -db with the run's files in its directory, no subject and subject handle 0 (line WQA D12a)", async () => {
    const { abi, calls } = fakeAbi();
    const runtime = runtimeWith(abi);
    const { run, close } = workerRun([]);
    const withoutSubject: WorkerRun = { ...run };
    delete (withoutSubject as { subject?: unknown }).subject;
    const files = [{ name: 'db/x.nal', file: new Blob(['TITLE x\n']) }];
    const database: WorkerRun = { ...withoutSubject, argv: ['tblastx', '-query', 'q.fa', '-db', 'db/x'], files };
    const opened = vi.spyOn(RunDirectory.prototype, 'open');
    const result = await runtime.run(database, () => undefined);
    expect(calls).toEqual(['validate tblastx -query q.fa -db db/x -num_threads 1', 'register query', 'run without a subject', 'release 101']);
    expect(result.subjectRetained).toBe(false);
    // The run's files are in its directory while it runs, and let go after it (review of D12a, low 6).
    expect(opened.mock.calls.map(([list]) => list)).toEqual([files, []]);
    opened.mockRestore();
    const loaded = (runtime as unknown as { loaded: { directory: RunDirectory } }).loaded;
    expect(loaded.directory.list).toEqual([]);
    // -parse_deflines reads the query alone; the shown IDs of the subjects are the database's.
    calls.length = 0;
    const together = await runtime.run({ ...database, argv: [...database.argv, '-parse_deflines'] }, () => undefined);
    expect(calls).toEqual([
      'validate tblastx -query q.fa -db db/x -parse_deflines -num_threads 1',
      'register inputs tblastx -query q.fa -db db/x -parse_deflines -num_threads 1',
      'run without a subject',
      'release 102',
    ]);
    expect(together.shownIds).toEqual({ query: ['Q.1'], subject: [] });
    // A subject with -db, and a run of neither: refused before the inputs are read.
    await expect(runtime.run({ ...database, subject: run.subject! }, () => undefined)).rejects.toThrow('a run of -db has no subject: give no subject');
    await expect(runtime.run({ ...withoutSubject }, () => undefined)).rejects.toThrow('this run needs a subject (or -db)');
    close();
  });
});

// --- the gateway ------------------------------------------------------------------------------------

describe('the gateway and the files of a run', () => {
  function gateway() {
    const posted: EngineCommand[] = [];
    const worker = {
      onmessage: null as ((event: MessageEvent<EngineEvent>) => void) | null,
      onerror: null,
      onmessageerror: null,
      postMessage: (command: EngineCommand) => void posted.push(command),
      terminate: () => undefined,
    };
    const assets = {
      serial: { url: '/serial.wasm', sha256: 'a'.repeat(64), size: 1, artifactSha256: 'a'.repeat(64) },
      threads: { url: '/threads.wasm', sha256: 'b'.repeat(64), size: 1, artifactSha256: 'c'.repeat(64) },
    };
    const engine = new WasmEngine({
      assets,
      control: { describe: () => Promise.reject(new Error('unused')), start: () => Promise.reject(new Error('unused')), validate: () => Promise.reject(new Error('unused')) },
      hardwareConcurrency: 8,
      createWorker: () => worker as unknown as Worker,
      loadModules: async () => ({ serial: { asset: assets.serial }, threads: { asset: assets.threads } }),
    });
    const input = { bytes: new Uint8Array([62, 115, 10, 65, 10]), sha256: 'f'.repeat(64), revisionIds: ['r1'], records: [{ id: 's', length: 1 }] };
    const request = (runId: string): EngineRunRequest => ({
      runId,
      argv: ['tblastx', '-query', 'q.fa', '-subject', 's.fa', '-conffile', 'c.ini', '-logfile', 'run.log'],
      query: input,
      subject: input,
      requestedThreads: 1,
      conffile: CONFFILE,
    });
    const send = (event: EngineEvent) => worker.onmessage!({ data: event } as MessageEvent<EngineEvent>);
    const runPosted = async (count: number) => {
      for (let i = 0; i < 100 && posted.filter((c) => c.type === 'run').length < count; i++) await new Promise((r) => setTimeout(r, 1));
      return posted.filter((c): c is WorkerRun => c.type === 'run')[count - 1]!;
    };
    return { engine, request, send, runPosted };
  }

  it('posts the registry file with the run, and hands on the files of a completed run and of a failed one', async () => {
    const { engine, request, send, runPosted } = gateway();
    const channel = new MessageChannel();
    const completed = engine.run(request('r1'), channel.port1, () => undefined);
    // The Blob itself: it crosses to the Engine worker without a copy of its bytes.
    const posted = (await runPosted(1)).conffile!;
    expect(posted.name).toBe('c.ini');
    expect(posted.file).toBe(CONFFILE.file);
    const result = { path: 'serial' as const, threads: 1, engineBuild: 'b', memory: { linearBytesBefore: 1, linearBytesAfter: 2, instanceRuns: 1 }, subjectRetained: false, messages: 'diagnostics' as const };
    send({ type: 'done', generation: 1, runId: 'r1', ok: true, result: { ...result, runFiles: [LOG] } });
    expect((await completed).runFiles).toEqual([LOG]);

    const failed = engine.run(request('r2'), channel.port1, () => undefined);
    await runPosted(2);
    send({ type: 'done', generation: 1, runId: 'r2', ok: false, error: { name: 'EngineCallError', message: 'LOSAT engine error: x', runFiles: [LOG] } });
    const error = await failed.catch((e: unknown) => e);
    expect(error).toMatchObject({ name: 'EngineCallError', message: 'LOSAT engine error: x' });
    expect(runFilesOf(error)).toEqual([LOG]);

    const plain = engine.run(request('r3'), channel.port1, () => undefined);
    await runPosted(3);
    send({ type: 'done', generation: 1, runId: 'r3', ok: false, error: { name: 'Error', message: 'no files' } });
    expect(runFilesOf(await plain.catch((e: unknown) => e))).toBeUndefined();
    channel.port1.close();
    channel.port2.close();
  });

  it('reads the files of an error only where a failed run gave them', () => {
    expect(runFilesOf(withRunFiles(new Error('x'), [LOG]))).toEqual([LOG]);
    expect(runFilesOf(new Error('x'))).toBeUndefined();
    expect(runFilesOf('x')).toBeUndefined();
    expect(runFilesOf(undefined)).toBeUndefined();
  });
});

// --- the names of the downloads ------------------------------------------------------------------

describe('the name of a downloaded file of a run', () => {
  it("starts as the names of the run's outputs do and ends with the file's own name", () => {
    const run = { number: 3, program: 'tblastx', queuedAt: Date.UTC(2026, 10, 4, 11, 15, 30), argv: ['tblastx'] };
    expect(runFileName(run, 'run.log')).toBe('losat-run3-tblastx-20261104T111530Z.logfile.run.log');
    expect(runFileName({ ...run, program: 'blastn', task: 'megablast' }, 'log')).toBe('losat-run3-blastn-megablast-20261104T111530Z.logfile.log');
    // A file named as an output's file (-logfile outfmt6.txt) is not saved under that output's name.
    expect(runFileName(run, 'outfmt6.txt')).not.toBe(outputFileName(run, 6));
    expect(runFileName(run, 'outfmt6.txt')).toBe('losat-run3-tblastx-20261104T111530Z.logfile.outfmt6.txt');
    // A name with a directory is saved as a name: each `/` is `_`; the engine already strips trailing `/`s (`logs/` -> `logs`).
    expect(runFileName(run, 'logs/run.log')).toBe('losat-run3-tblastx-20261104T111530Z.logfile.logs_run.log');
    expect(runFileName(run, '/tmp/run.log')).toBe('losat-run3-tblastx-20261104T111530Z.logfile._tmp_run.log');
    expect(runFileName(run, '/dev/stdout')).toBe('losat-run3-tblastx-20261104T111530Z.logfile._dev_stdout');
    expect(runFileName(run, 'logs')).toBe('losat-run3-tblastx-20261104T111530Z.logfile.logs');
    expect(runFileName(run, 'a/b')).toBe('losat-run3-tblastx-20261104T111530Z.logfile.a_b');
    // The file of -export_search_strategy is named by its option (line RG item F1w).
    expect(runFileName(run, 'out/s.asn', 'export_search_strategy')).toBe('losat-run3-tblastx-20261104T111530Z.export_search_strategy.out_s.asn');
  });
});

// --- the built serial reactor ---------------------------------------------------------------------

const reactors = findReactors();
const LOSAT = resolve(import.meta.dirname, '../../../../LOSAT');

describe.skipIf(reactors === undefined)('the start and the files of a run in the built serial reactor', () => {
  const module = reactors === undefined ? undefined : new WebAssembly.Module(reactors.serial.bytes as BufferSource);
  const fasta = (name: string) => new Uint8Array(readFileSync(resolve(LOSAT, 'tests/fasta/outfmt11', name)));

  /** The streams of a TBLASTX run of `words` on the inputs of the adapter's tests (run_start_tests.rs). */
  function run(abi: ReactorAbi, words: readonly string[], subject = fasta('arch_tbx_s.fna')) {
    const q = abi.register('tblastx', ROLE_QUERY, true, fasta('arch_tbx_q.fna'));
    const s = abi.register('tblastx', ROLE_SUBJECT, true, subject);
    const streams = new Map<number, Uint8Array[]>();
    let error: unknown;
    try {
      abi.run(['tblastx', '-query', 'q.fa', '-subject', 's.fa', ...words, '-num_threads', '1'], q.handle, s.handle, (stream, bytes) => {
        streams.set(stream, [...(streams.get(stream) ?? []), bytes.slice()]);
      });
    } catch (e) {
      error = e;
    } finally {
      abi.release(q.handle);
      abi.release(s.handle);
    }
    const joined = (stream: number) => {
      const chunks = streams.get(stream) ?? [];
      const bytes = new Uint8Array(chunks.reduce((total, chunk) => total + chunk.length, 0));
      let offset = 0;
      for (const chunk of chunks) {
        bytes.set(chunk, offset);
        offset += chunk.length;
      }
      return bytes;
    };
    return { error, streams, joined };
  }

  it("starts an argv, and sends a TBLASTX run's -logfile file as an archive that the reader reads, also when the run fails", async () => {
    const { abi } = await instantiateSerial(module!);
    expect(abi.start(['tblastx', '-conffile', 'a.ini', '-conffile=c.ini', '-logfile', 'x.log', '-dryrun'])).toEqual({
      ok: true,
      start: { conffile: 'c.ini', logFile: 'x.log', logFileName: 'x.log', strategyFile: null, dryRun: true, registerInputs: false, prints: null, database: false },
    });
    expect(abi.start(['tblastx', '-logfile', 'logs/'])).toMatchObject({ ok: true, start: { logFile: 'logs/', logFileName: 'logs' } });
    expect(abi.start(['psiblast'])).toMatchObject({ ok: false });
    // The command line's help and version are printed in place of a search (line WQA's D14, 2026-11-20).
    expect(abi.start(['blastn', '-version', '-dryrun'])).toEqual({ ok: true, start: { conffile: null, logFile: null, logFileName: null, strategyFile: null, dryRun: false, registerInputs: false, prints: 'version', database: false } });
    expect(abi.start(['tblastx', '-dryrun', '-h'])).toMatchObject({ ok: true, start: { dryRun: false, prints: 'help' } });
    // -parse_deflines reads the inputs of a run together, but not those of a dry run (line WQA's sweep, D7, 2026-11-19).
    expect(abi.start(['blastn', '-query', 'q', '-subject', 's', '-parse_deflines'])).toMatchObject({ ok: true, start: { registerInputs: true } });
    expect(abi.start(['blastn', '-query', 'q', '-subject', 's', '-parse_deflines', '-dryrun'])).toMatchObject({ ok: true, start: { registerInputs: false } });

    const warned = run(abi, ['-outfmt', '6', '-max_target_seqs', '1', '-logfile', 'run.log']);
    expect(warned.error).toBeUndefined();
    expect(readStoredZip(warned.joined(RUN_FILES_STREAM))).toEqual([{ name: 'run.log', bytes: encode('Warning: [tblastx] Examining 5 or more matches is recommended\n') }]);
    expect(warned.streams.has(3)).toBe(false);

    // A subject record without residues: NCBI's error ends the run, and TBLASTX writes it into the file.
    const failed = run(abi, ['-logfile', 'run.log'], encode('>s1 alpha\n'));
    expect(failed.error).toBeInstanceOf(EngineCallError);
    const [file] = readStoredZip(failed.joined(RUN_FILES_STREAM));
    expect(file!.name).toBe('run.log');
    expect(new TextDecoder().decode(file!.bytes)).toBe(
      'Warning: [tblastx] Subject_1 s1 alpha: Subject sequence contains no data\nLOSAT engine error: The average subject length is too short\n',
    );
    expect(failed.streams.has(3)).toBe(false);
    // Without -logfile the warning before the error is on stream 3, where the command line writes it
    // (NCBI BLAST+ 2.17.0; line WQA's sweep, D3, 2026-11-19), and the error is the run's.
    const plain = run(abi, [], encode('>s1 alpha\n'));
    expect((plain.error as Error).message).toBe('LOSAT engine error: The average subject length is too short\n');
    expect(new TextDecoder().decode(plain.joined(3))).toBe('Warning: [tblastx] Subject_1 s1 alpha: Subject sequence contains no data\n');
  }, 120_000);

  // DECISIONS E1w3b-4: the file goes to an instance of the validation's own; the shared instance neither
  // holds it nor grows, however large the file.
  it('reads the registry file of -conffile on an instance of its own: the shared instance holds no file and does not grow', async () => {
    const shared = (await instantiateSerial(module!)).abi;
    const control = reactorControl(async () => shared, async () => (await instantiateSerial(module!)).abi);
    const argv = ['tblastx', '-query', 'q.fa', '-subject', 's.fa', '-conffile', 'settings.ini'];
    const missing = { ok: false, message: 'Critical: (106.23) Registry file "settings.ini" not found, aborting.' };
    expect(await control.validate(argv)).toEqual(missing);
    const before = shared.memoryBytes();
    expect(await control.validate(argv, { name: 'settings.ini', file: new Blob(['[BLAST]\n']) })).toEqual({ ok: true, dataLoaders: true });
    // The registry's `[BLAST] DATA_LOADERS=none` turns the data loaders off, as the CLI reads it (line WQA's D16, 2026-11-20).
    expect(await control.validate(argv, { name: 'settings.ini', file: new Blob(['[BLAST]\nDATA_LOADERS=none\n']) })).toEqual({ ok: true, dataLoaders: false });
    expect((await control.validate(argv, { name: 'settings.ini', file: new Blob(['[x\n']) })).ok).toBe(false);
    // A file of 24 MiB and 26 bytes (four chunks, the last of 26 bytes): a section and comment lines.
    const large = new Blob(['[BLAST]\n', '; a comment line of the registry file\n'.repeat(Math.ceil((24 * 1024 * 1024) / 38))]);
    expect(large.size).toBe(3 * FILE_CHUNK_BYTES + 26);
    expect(await control.validate(argv, { name: 'settings.ini', file: large })).toEqual({ ok: true, dataLoaders: true });
    expect(shared.memoryBytes()).toBe(before);
    expect(await control.validate(argv)).toEqual(missing);
  }, 120_000);
});
