// BLASTN's -window_masker_db in the application (line NP WW1, DESIGN_ww1 §5-§7): the form's
// file field and the argv it makes, the draft's missing-file rule, the run snapshot, the
// gateway's message to the Engine worker, and the retention of the counts file in the
// engine runtime (load, reuse, release on another program or other counts, with the instance).
// With LOSAT_WEB_REACTORS, the same runs go to the built serial reactor.
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { findReactors } from '../../build/reactors';
import { Coordinator, type EnqueueAllResult, type SearchRequest } from '../../src/application/coordinator';
import { SearchDraft } from '../../src/application/draft';
import { buildArgv } from '../../src/domain/argv';
import { changedFields, formParameters, setField, writtenValue, type EngineOption } from '../../src/domain/parameters';
import { programById, WINDOW_MASKER_FLAG } from '../../src/domain/programs';
import { sha256Hex } from '../../src/infra/browser/platform';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { COUNTS_CHUNK_BYTES, EngineRuntime } from '../../src/infra/engine-worker/runtime';
import { countsKey, WasmEngine } from '../../src/infra/engine-worker/gateway';
import type { EngineCommand, EngineEvent, WorkerCounts, WorkerRun } from '../../src/infra/engine-worker/protocol';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeAccessionFetcher } from '../../src/infra/fake/fake-fetcher';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import { NodeInflater } from './support/node-inflater';
import { EngineCallError, EngineStoppedError, FILE_WINDOW_MASKER_COUNTS, ReactorAbi, ROLE_QUERY, ROLE_SUBJECT, type AbiExports, type FileInfo } from '../../src/infra/reactor/abi';
import { instantiateSerial } from '../../src/infra/reactor/instance';
import type { EngineRunRequest } from '../../src/ports/engine';

// --- the form and the argv ---------------------------------------------------------------------

describe('the WindowMasker counts field', () => {
  const blastn = programById('blastn');
  const options: EngineOption[] = [{ flag: WINDOW_MASKER_FLAG, takesValue: true }];

  it('is NCBI\'s "Species-specific repeats" place in the Filter row, of the kind file', () => {
    const section = blastn.sections.find((s) => s.title === 'Filters and Masking')!;
    const field = section.fields.find((f) => f.flag === WINDOW_MASKER_FLAG)!;
    expect(field).toMatchObject({ kind: 'file', row: 'Filter', label: 'Species-specific repeats filter for:' });
    expect(field.file?.label).toBe('WindowMasker counts');
    // The other programs have no such field.
    expect(
      ['blastp', 'blastx', 'tblastn', 'tblastx'].some((id) => programById(id as 'blastp').sections.some((s) => s.fields.some((f) => f.flag === WINDOW_MASKER_FLAG))),
    ).toBe(false);
  });

  it('writes -window_masker_db and the file name, which marks the field as changed; checked without a file it writes nothing', () => {
    const field = blastn.sections.flatMap((s) => s.fields).find((f) => f.flag === WINDOW_MASKER_FLAG)!;
    expect(writtenValue(field, true, undefined)).toBeUndefined();
    expect(writtenValue(field, '', undefined)).toBeUndefined();
    expect(writtenValue(field, 'mg1655_u9.obinary', undefined)).toBe('mg1655_u9.obinary');

    const checked = setField(blastn, {}, WINDOW_MASKER_FLAG, true, options);
    expect(formParameters(blastn, checked, options)).toEqual([]);
    expect(changedFields(blastn, checked, options, 'algorithm')).toBe(0);
    const chosen = setField(blastn, checked, WINDOW_MASKER_FLAG, 'mg1655_u9.obinary', options);
    expect(formParameters(blastn, chosen, options)).toEqual([[WINDOW_MASKER_FLAG, 'mg1655_u9.obinary']]);
    expect(changedFields(blastn, chosen, options, 'algorithm')).toBe(1);
    expect(
      buildArgv({ program: 'blastn', queryName: 'q.fa', subjectName: 's.fa', parameters: formParameters(blastn, chosen, options) }),
    ).toEqual(['blastn', '-query', 'q.fa', '-subject', 's.fa', '-window_masker_db', 'mg1655_u9.obinary']);
  });
});

// --- the draft ----------------------------------------------------------------------------------

function draftSetup() {
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
  const fake = new FakeEngine();
  const engine = { describe: fake.describe.bind(fake), validate: fake.validate.bind(fake), prepare: () => undefined };
  let id = 0;
  const coordinator = new Coordinator({
    engine: { ...engine, run: fake.run.bind(fake), cancel: fake.cancel.bind(fake) },
    data,
    downloader: { save: () => undefined },
    now: () => 0,
    newRunId: () => `run-${++id}`,
  });
  const requests: SearchRequest[] = [];
  const draft = new SearchDraft({
    engine,
    data,
    enqueueAll: async (batch): Promise<EnqueueAllResult> => {
      requests.push(...batch);
      return coordinator.enqueueAll(batch);
    },
    fetcher: new FakeAccessionFetcher(),
    now: () => 0,
    debounceMs: 0,
  });
  return { draft, coordinator, requests };
}

const counts = () => new File([new Uint8Array([1, 2, 3, 4, 5])], 'mg1655_u9.obinary', { lastModified: 1_700_000_000_000 });

describe('the draft and the counts file', () => {
  async function ready() {
    const s = draftSetup();
    s.draft.setPaste('query', '>q\nACGTACGT\n');
    s.draft.setPaste('subject', '>s\nACGTACGT\n');
    await s.draft.idle();
    return s;
  }

  it('does not start a search while the box is checked and no file is chosen', async () => {
    const { draft, requests } = await ready();
    draft.setField(WINDOW_MASKER_FLAG, true);
    await draft.idle();
    expect(draft.readiness()).toBe('Choose the WindowMasker counts file.');
    expect(await draft.submit()).toEqual({ ok: false, message: 'Choose the WindowMasker counts file.' });
    expect(requests).toHaveLength(0);
    // Unchecked, it is a search without the filter.
    draft.setFile(WINDOW_MASKER_FLAG, undefined);
    expect(draft.readiness()).toBeUndefined();
    expect((await draft.submit()).ok).toBe(true);
    expect(requests[0]!.parameters).toEqual([]);
    expect(requests[0]!.windowMasker).toBeUndefined();
  });

  it('writes the chosen file\'s name to the argv and freezes the file in the run snapshot, session only', async () => {
    const { draft, coordinator, requests } = await ready();
    const file = counts();
    draft.setField(WINDOW_MASKER_FLAG, true);
    draft.setFile(WINDOW_MASKER_FLAG, file);
    await draft.idle();
    expect(draft.state.get().values.blastn[WINDOW_MASKER_FLAG]).toBe('mg1655_u9.obinary');
    expect(draft.state.get().validation).toEqual({ state: 'ok' });
    const result = await draft.submit();
    expect(result.ok).toBe(true);
    expect(requests[0]!.parameters).toEqual([[WINDOW_MASKER_FLAG, 'mg1655_u9.obinary']]);
    expect(requests[0]!.windowMasker).toBe(file);
    const snapshot = coordinator.state.get().runs[0]!.snapshot;
    expect(snapshot.argv.slice(-2)).toEqual(['-window_masker_db', 'mg1655_u9.obinary']);
    expect(snapshot.windowMasker).toEqual({ name: 'mg1655_u9.obinary', file, size: 5, lastModified: 1_700_000_000_000 });
  });

  it('records what the run did with the counts, and keeps the counts in the session only', async () => {
    const { draft, coordinator } = await ready();
    draft.setFile(WINDOW_MASKER_FLAG, counts());
    await draft.submit();
    for (let i = 0; i < 200 && coordinator.state.get().runs[0]?.status !== 'completed'; i++) await new Promise((r) => setTimeout(r, 5));
    const run = coordinator.state.get().runs[0]!;
    expect(run.status).toBe('completed');
    expect(run.record.windowMasker).toEqual({
      name: 'mg1655_u9.obinary',
      size: 5,
      sha256: await sha256Hex(new Uint8Array([1, 2, 3, 4, 5])),
      retained: false,
      loadMs: 0,
    });
  });

  it('lets the file go with "Restore default search parameters"', async () => {
    const { draft } = await ready();
    draft.setFile(WINDOW_MASKER_FLAG, counts());
    draft.resetFields([WINDOW_MASKER_FLAG]);
    expect(draft.state.get().files.blastn[WINDOW_MASKER_FLAG]).toBeUndefined();
    expect(draft.state.get().values.blastn[WINDOW_MASKER_FLAG]).toBeUndefined();
    expect(draft.readiness()).toBeUndefined();
  });

  it('refuses the option typed in "Additional options": the form has it', async () => {
    const { draft } = await ready();
    draft.setExtra('-window_masker_db x.obinary');
    expect(() => draft.extraWords()).toThrow('-window_masker_db is set in the form: Species-specific repeats filter for:');
  });
});

// --- the gateway --------------------------------------------------------------------------------

describe('the gateway and the counts file', () => {
  it('posts the file with a key of its name, size and time, and passes the run\'s counts record on', async () => {
    const posted: EngineCommand[] = [];
    const worker = {
      posted,
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
      control: { describe: () => Promise.reject(new Error('unused')), validate: () => Promise.reject(new Error('unused')) },
      hardwareConcurrency: 8,
      createWorker: () => worker as unknown as Worker,
      loadModules: async () => ({ serial: { asset: assets.serial }, threads: { asset: assets.threads } }),
    });
    const file = counts();
    const input = { bytes: new Uint8Array([62, 115, 10, 65, 10]), sha256: 'f'.repeat(64), revisionIds: ['r1'], records: [{ id: 's', length: 1 }] };
    const channel = new MessageChannel();
    const request: EngineRunRequest = {
      runId: 'r1',
      argv: ['blastn', '-query', 'q.fa', '-subject', 's.fa', '-window_masker_db', 'mg1655_u9.obinary'],
      query: input,
      subject: input,
      requestedThreads: 1,
      windowMasker: { name: 'mg1655_u9.obinary', file, size: file.size, lastModified: file.lastModified },
    };
    const run = engine.run(request, channel.port1, () => undefined);
    for (let i = 0; i < 100 && !posted.some((c) => c.type === 'run'); i++) await new Promise((r) => setTimeout(r, 1));
    const command = posted.find((c): c is WorkerRun => c.type === 'run')!;
    expect(command.windowMasker).toEqual({ name: 'mg1655_u9.obinary', blob: file, key: JSON.stringify(['mg1655_u9.obinary', 5, 1_700_000_000_000]) });
    expect(command.windowMasker!.key).toBe(countsKey(request.windowMasker!));
    const record = { name: 'mg1655_u9.obinary', size: 5, sha256: '1'.repeat(64), retained: true, loadMs: 0 };
    worker.onmessage!({
      data: {
        type: 'done',
        generation: 1,
        runId: 'r1',
        ok: true,
        result: { path: 'serial', threads: 1, engineBuild: 'b', memory: { linearBytesBefore: 1, linearBytesAfter: 2, instanceRuns: 1 }, subjectRetained: false, windowMasker: record },
      },
    } as MessageEvent<EngineEvent>);
    expect((await run).windowMasker).toEqual(record);
    channel.port1.close();
    channel.port2.close();
  });
});

// --- retention in the runtime ------------------------------------------------------------------

/** An ABI binding that records the calls a search makes, in order, and has no engine behind it. */
function fakeAbi(options: { failEnd?: boolean; stopOnBegin?: boolean } = {}) {
  const calls: string[] = [];
  let handle = 100;
  let stopped: EngineStoppedError | undefined;
  const infos = new Map<number, FileInfo>();
  const abi = {
    get stoppedBy() {
      return stopped;
    },
    validate: () => {
      calls.push('validate');
      return undefined;
    },
    memoryBytes: () => 1,
    register: (_program: string, role: number, bytes: Uint8Array) => {
      calls.push(`register ${role === ROLE_QUERY ? 'query' : 'subject'}`);
      return { handle: ++handle, records: role === ROLE_QUERY ? [{ id: 'q', length: bytes.length }] : [{ id: 's', length: bytes.length }] };
    },
    release: (h: number) => void calls.push(`release ${h}`),
    run: () => void calls.push('run'),
    fileBegin: (kind: number, name: string, size: number) => {
      calls.push(`begin ${kind} ${name} ${size}`);
      if (options.stopOnBegin) {
        stopped = new EngineStoppedError('trap');
        throw stopped;
      }
      return ++handle;
    },
    chunkBuffer: (capacity: number) => {
      calls.push(`buffer ${capacity}`);
      return { ptr: 1, capacity };
    },
    freeChunkBuffer: () => void calls.push('free buffer'),
    fileChunk: (file: number, bytes: Uint8Array) => void calls.push(`chunk ${file} ${bytes.length}`),
    fileEnd: (file: number) => {
      calls.push(`end ${file}`);
      if (options.failEnd) throw new EngineCallError('the input file ended after 1 of its 3 bytes');
      const info: FileInfo = { handle: file, name: 'n', kind: 'windowmasker_counts', size: 3, sha256: `sha${file}`, stat_type: 'obinary', unit_size: 9, error: null };
      infos.set(file, info);
      return info;
    },
    fileRelease: (file: number) => void calls.push(`file release ${file}`),
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
    faultChannel: 'faults',
    keeperChannel: 'keeper',
  });
  // The runtime holds this instance as one that a search loaded.
  (runtime as unknown as { loaded: unknown }).loaded = {
    path: 'serial',
    instance: { abi, wasi: { clear: () => undefined } },
    runs: 0,
    retained: undefined,
    counts: undefined,
    ran: new Set(),
  };
  return runtime;
}

const counts3 = (name: string, bytes = 3, lastModified = 1): WorkerCounts => {
  const blob = new Blob([new Uint8Array(bytes)]);
  return { name, blob, key: JSON.stringify([name, bytes, lastModified]) };
};

let runNumber = 0;
function workerRun(program: string, windowMasker?: WorkerCounts): WorkerRun {
  const channel = new MessageChannel();
  const input = (role: string) => ({ bytes: new Uint8Array([62, 10, 65]), records: [{ id: role, length: 3 }], key: `k-${role}` });
  // The fake binding answers with the id of the role.
  const run: WorkerRun = {
    type: 'run',
    generation: 1,
    runId: `run-${++runNumber}`,
    program,
    argv: [program, '-query', 'q.fa', '-subject', 's.fa', ...(windowMasker === undefined ? [] : ['-window_masker_db', windowMasker.name])],
    threads: 1,
    query: input('q'),
    subject: input('s'),
    ...(windowMasker === undefined ? {} : { windowMasker }),
    output: channel.port1,
  };
  channel.port2.close();
  return run;
}

describe('the counts file in the Engine runtime', () => {
  const nothing = () => undefined;

  it('gives the file to the engine in the preparing phase before validate and the query, in chunks of one buffer', async () => {
    const { abi, calls } = fakeAbi();
    const runtime = runtimeWith(abi);
    const phases: string[] = [];
    const result = await runtime.run(workerRun('blastn', counts3('a.obinary', COUNTS_CHUNK_BYTES + 5)), (phase) => {
      phases.push(phase);
      calls.push(`phase ${phase}`);
    });
    const size = COUNTS_CHUNK_BYTES + 5;
    // validate after the file: it checks the file where the command line does (ABI v2 §4).
    expect(calls.slice(0, 10)).toEqual([
      'phase preparing',
      `begin ${FILE_WINDOW_MASKER_COUNTS} a.obinary ${size}`,
      `buffer ${COUNTS_CHUNK_BYTES}`,
      `chunk 101 ${COUNTS_CHUNK_BYTES}`,
      'chunk 101 5',
      'end 101',
      'free buffer',
      'validate',
      'register query',
      'register subject',
    ]);
    expect(calls.indexOf('phase running')).toBeGreaterThan(calls.indexOf('register subject'));
    expect(result.windowMasker).toMatchObject({ name: 'a.obinary', size: 3, sha256: 'sha101', retained: false });
    expect(result.windowMasker!.loadMs).toBeGreaterThanOrEqual(0);
    // The result tells the gateway how many bytes of counts the instance holds (the renewal's mark allows for them).
    expect(result.countsBytes).toBe(3);
  });

  it('uses the counts it holds when the key is the same: nothing is given again, and the record says so', async () => {
    const { abi, calls } = fakeAbi();
    const runtime = runtimeWith(abi);
    await runtime.run(workerRun('blastn', counts3('a.obinary')), nothing);
    calls.length = 0;
    const result = await runtime.run(workerRun('blastn', counts3('a.obinary')), nothing);
    expect(calls.filter((call) => /^(begin|chunk|end|file release)/.test(call))).toEqual([]);
    expect(result.windowMasker).toEqual({ name: 'a.obinary', size: 3, sha256: 'sha101', retained: true, loadMs: 0 });
  });

  it('reads a file again whose name, size and time are those of the held one but whose bytes differ', async () => {
    const { abi, calls } = fakeAbi();
    const runtime = runtimeWith(abi);
    const same = (fill: number): WorkerCounts => ({
      name: 'a.obinary',
      blob: new Blob([new Uint8Array(3).fill(fill)]),
      key: JSON.stringify(['a.obinary', 3, 1]),
    });
    await runtime.run(workerRun('blastn', same(1)), nothing);
    calls.length = 0;
    const result = await runtime.run(workerRun('blastn', same(2)), nothing);
    expect(calls.filter((call) => /^(begin|end|file release)/.test(call)).map((call) => call.replace(/ \d+$/, ''))).toEqual([
      'file release',
      'begin 1 a.obinary',
      'end',
    ]);
    expect(result.windowMasker).toMatchObject({ retained: false });
  });

  it('keeps the counts through a BLASTN search without -window_masker_db', async () => {
    const { abi, calls } = fakeAbi();
    const runtime = runtimeWith(abi);
    await runtime.run(workerRun('blastn', counts3('a.obinary')), nothing);
    calls.length = 0;
    const without = await runtime.run(workerRun('blastn'), nothing);
    expect(without.windowMasker).toBeUndefined();
    expect(calls.filter((call) => call.startsWith('file release'))).toEqual([]);
    // And the next search that names them finds them held.
    expect((await runtime.run(workerRun('blastn', counts3('a.obinary')), nothing)).windowMasker?.retained).toBe(true);
  });

  it('releases the held counts before other counts are given, so that two are never held', async () => {
    const { abi, calls } = fakeAbi();
    const runtime = runtimeWith(abi);
    await runtime.run(workerRun('blastn', counts3('a.obinary')), nothing);
    calls.length = 0;
    // The same name with another size or time is another file.
    for (const next of [counts3('a.obinary', 4), counts3('a.obinary', 4, 2), counts3('b.obinary', 4, 2)]) {
      calls.length = 0;
      const result = await runtime.run(workerRun('blastn', next), nothing);
      const file = calls.filter((call) => /^(begin|file release)/.test(call));
      expect(file[0]).toMatch(/^file release \d+$/);
      expect(file[1]).toMatch(/^begin 1 /);
      expect(result.windowMasker?.retained).toBe(false);
    }
  });

  it('releases the counts before a search of another program, which does not use them', async () => {
    const { abi, calls } = fakeAbi();
    const runtime = runtimeWith(abi);
    await runtime.run(workerRun('blastn', counts3('a.obinary')), nothing);
    calls.length = 0;
    const result = await runtime.run(workerRun('tblastx'), nothing);
    expect(calls[0]).toBe('file release 101');
    expect(calls.indexOf('file release 101')).toBeLessThan(calls.indexOf('register query'));
    expect(result.windowMasker).toBeUndefined();
    expect(result.countsBytes).toBeUndefined();
    // Back in BLASTN, the file is given again.
    calls.length = 0;
    expect((await runtime.run(workerRun('blastn', counts3('a.obinary')), nothing)).windowMasker?.retained).toBe(false);
    expect(calls.some((call) => call.startsWith('begin'))).toBe(true);
  });

  it('does not hold the counts of an instance that was ended: the next instance holds none', async () => {
    const { abi } = fakeAbi();
    const runtime = runtimeWith(abi);
    await runtime.run(workerRun('blastn', counts3('a.obinary')), nothing);
    expect((runtime as unknown as { loaded: { counts: unknown } }).loaded.counts).toBeDefined();
    runtime.dispose();
    expect((runtime as unknown as { loaded: unknown }).loaded).toBeUndefined();
  });

  it('reports a file that the engine refuses, holds nothing, and gives the file again to the next search', async () => {
    const failing = fakeAbi({ failEnd: true });
    const runtime = runtimeWith(failing.abi);
    await expect(runtime.run(workerRun('blastn', counts3('a.obinary')), nothing)).rejects.toThrow('the input file ended after 1 of its 3 bytes');
    expect(failing.calls).toContain('free buffer');
    expect(failing.calls.some((call) => call.startsWith('register'))).toBe(false);
    expect((runtime as unknown as { loaded: { counts: unknown } }).loaded.counts).toBeUndefined();
  });

  it('ends the instance that stopped while it was given a file', async () => {
    const { abi } = fakeAbi({ stopOnBegin: true });
    const runtime = runtimeWith(abi);
    await expect(runtime.run(workerRun('blastn', counts3('a.obinary')), nothing)).rejects.toThrow('trap');
    expect((runtime as unknown as { loaded: unknown }).loaded).toBeUndefined();
  });

  it('gives a file of no bytes: begin, end, and a buffer of one byte', async () => {
    const { abi, calls } = fakeAbi();
    const runtime = runtimeWith(abi);
    await runtime.run(workerRun('blastn', counts3('empty.obinary', 0)), nothing);
    expect(calls.filter((call) => /^(begin|buffer|chunk|end)/.test(call))).toEqual(['begin 1 empty.obinary 0', 'buffer 1', 'end 101']);
  });
});

describe('the ABI binding of the file calls', () => {
  function binding() {
    const memory = new WebAssembly.Memory({ initial: 4 });
    const log: string[] = [];
    let next = 4096;
    const ref: { abi?: ReactorAbi } = {};
    const json = new TextEncoder().encode('{"handle":5,"name":"n","kind":"windowmasker_counts","size":3,"sha256":"ab","stat_type":"obinary","unit_size":9,"error":null}');
    const exports = {
      memory,
      losat_web2_alloc: (len: number) => {
        const ptr = next;
        next += len;
        log.push(`alloc ${len}`);
        return ptr;
      },
      losat_web2_dealloc: (ptr: number, len: number) => void log.push(`dealloc ${ptr} ${len}`),
      losat_web2_file_begin: (kind: number, ptr: number, len: number, size: number) => {
        log.push(`begin ${kind} ${new TextDecoder().decode(new Uint8Array(memory.buffer, ptr, len))} ${size}`);
        return 5;
      },
      losat_web2_file_chunk: (file: number, ptr: number, len: number) => {
        log.push(`chunk ${file} ${[...new Uint8Array(memory.buffer, ptr, len)].join(',')}`);
        return 0;
      },
      losat_web2_file_end: () => {
        new Uint8Array(memory.buffer, 60_000, json.length).set(json);
        ref.abi!.emit(2, 60_000, json.length);
        return 0;
      },
      losat_web2_file_release: (file: number) => void log.push(`release ${file}`),
    } as unknown as AbiExports;
    ref.abi = new ReactorAbi(exports);
    return { abi: ref.abi, log };
  }

  it('begins with the name as the argv gives it, copies every chunk into one buffer, and reads the file JSON', () => {
    const { abi, log } = binding();
    expect(abi.fileBegin(1, 'dir/ü.obinary', 3)).toBe(5);
    const buffer = abi.chunkBuffer(2);
    abi.fileChunk(5, new Uint8Array([1, 2]), buffer);
    abi.fileChunk(5, new Uint8Array([3]), buffer);
    expect(abi.fileEnd(5)).toMatchObject({ handle: 5, sha256: 'ab', stat_type: 'obinary', unit_size: 9, error: null });
    abi.freeChunkBuffer(buffer);
    abi.fileRelease(5);
    // The name takes 14 bytes of UTF-8; the one buffer is taken once and given back once.
    expect(log).toEqual(['alloc 14', 'begin 1 dir/ü.obinary 3', 'dealloc 4096 14', 'alloc 2', 'chunk 5 1,2', 'chunk 5 3', 'dealloc 4110 2', 'release 5']);
  });

  it('refuses a chunk larger than the buffer and a size the engine cannot hold', () => {
    const { abi } = binding();
    expect(() => abi.fileChunk(5, new Uint8Array(3), abi.chunkBuffer(2))).toThrow(/does not fit/);
    expect(() => abi.fileBegin(1, 'big', 2 ** 32)).toThrow(EngineCallError);
  });
});

// --- the built serial reactor -------------------------------------------------------------------

const reactors = findReactors();
const LOSAT = resolve(import.meta.dirname, '../../../../LOSAT');

describe.skipIf(reactors === undefined)('the counts file in the built serial reactor', () => {
  const module = reactors === undefined ? undefined : new WebAssembly.Module(reactors.serial.bytes as BufferSource);
  const counts = new Uint8Array(readFileSync(resolve(LOSAT, 'tests/winmask/counts/mg1655_u9.obinary')));
  const query = readFileSync(resolve(LOSAT, 'tests/fixtures/blastn_regression/inputs/wm_query.fa'));
  const subject = readFileSync(resolve(LOSAT, 'tests/fixtures/blastn_regression/inputs/s60k.fa'));

  function search(abi: ReactorAbi, argv: readonly string[]): string {
    const q = abi.register('blastn', ROLE_QUERY, query);
    const s = abi.register('blastn', ROLE_SUBJECT, subject);
    const chunks: Uint8Array[] = [];
    try {
      abi.run([...argv, '-num_threads', '1'], q.handle, s.handle, (stream, bytes) => {
        if (stream === 6) chunks.push(bytes.slice());
      });
    } finally {
      abi.release(q.handle);
      abi.release(s.handle);
    }
    return Buffer.concat(chunks).toString();
  }

  it('masks with the counts that a run names (NCBI\'s frozen outfmt 6), refuses a second registration of the name, and forgets the counts when released', async () => {
    const { abi } = await instantiateSerial(module!);
    const base = ['blastn', '-query', 'q.fa', '-subject', 's.fa', '-task', 'blastn', '-outfmt', '6'];
    const plain = search(abi, base);

    const handle = abi.fileBegin(FILE_WINDOW_MASKER_COUNTS, 'counts/mg1655_u9.obinary', counts.length);
    const buffer = abi.chunkBuffer(100_000);
    for (let offset = 0; offset < counts.length; offset += 100_000) abi.fileChunk(handle, counts.subarray(offset, offset + 100_000), buffer);
    abi.freeChunkBuffer(buffer);
    const info = abi.fileEnd(handle);
    expect(info).toMatchObject({ kind: 'windowmasker_counts', size: counts.length, stat_type: 'obinary', error: null });
    expect(info.sha256).toBe(await sha256Hex(counts));

    const masked = search(abi, [...base, '-window_masker_db', 'counts/mg1655_u9.obinary']);
    // NCBI's frozen outfmt 6 of this argv (row wm.blastn of the blastn_regression manifest).
    expect(await sha256Hex(new TextEncoder().encode(masked))).toBe('7774c0f582c486fd1f93d3412fd2b508937139e97005a158dd00a632b20f17cb');
    expect(masked).not.toBe(plain);

    // A name that is registered twice is refused until it is released.
    expect(() => abi.fileBegin(FILE_WINDOW_MASKER_COUNTS, 'counts/mg1655_u9.obinary', 1)).toThrow(/release it before registering another/);
    abi.fileRelease(handle);
    expect(() => search(abi, [...base, '-window_masker_db', 'counts/mg1655_u9.obinary'])).toThrow();
  }, 120_000);
});
