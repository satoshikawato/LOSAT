// Session files through the application (application/session.ts, S15 items 4 and 6): a working
// session with two runs, candidates and notes is saved and opened in a fresh application whose
// engine throws if it is asked anything. The runs come back without a search (same outputs byte for
// byte, same HSP records and tables, the candidates and notes restored); the extraction of
// sequences is refused until the original FASTA is chosen again and matches, and a wrong file is
// never attached; the alignment export needs no original. A damaged or hostile file is refused and
// leaves nothing behind.
import { gunzipSync, gzipSync } from 'node:zlib';
import { describe, expect, it } from 'vitest';
import { CandidateTray } from '../../src/application/candidates';
import { Coordinator, type SearchRequest } from '../../src/application/coordinator';
import { ResultsBrowser, type HspId } from '../../src/application/results';
import { Session } from '../../src/application/session';
import {
  blockHeader,
  containerEnd,
  containerHeader,
  runBlockName,
  SESSION_STREAMS,
  SessionFileReader,
  type SessionCandidate,
  type SessionManifest,
} from '../../src/domain/session-file';
import type { RunStatus } from '../../src/domain/run';
import { sha256Hex } from '../../src/infra/browser/platform';
import { browserCompression } from '../../src/infra/browser/compression';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import type { DataGateway } from '../../src/ports/data';
import type { EngineGateway } from '../../src/ports/engine';
import { memoryDownloader, type SavedFile } from './support/memory-downloader';

const encoder = new TextEncoder();
const decoder = new TextDecoder();

const QUERIES = '>q1 first query\nACGTACGTACGTAAAC\n>q2\nGGGTTTAAACCC\n';
const FILE_A = '>sA one\nACGTACGTACGTAAACCCGGGTTT\nAACCGGTT\n>sB left out\nTTTTGGGGCC\n';
const FILE_B = '>sC\nACGTTTGACAGGCATTACGA\n';
const FILE_C = '>sD nucleotides\nACGTACGTACGTACGTACGTACGTACGTAC\nACGTAC\n';
const PROTEIN = '>p1\nMKVLAAGIVGLL\n';

/** An engine that fails the test if a loaded session asks it anything. */
function forbiddenEngine() {
  const calls: string[] = [];
  const refuse = (what: string) => {
    calls.push(what);
    return Promise.reject(new Error(`the engine was asked to ${what}`));
  };
  const engine: EngineGateway = {
    describe: () => refuse('describe'),
    validate: () => refuse('validate'),
    run: () => refuse('run'),
    cancel: () => void calls.push('cancel'),
  };
  return { engine, calls };
}

/** A working session: the data layer, the coordinator, the tray and the session, as composition.ts makes them. */
function world(engine: EngineGateway = new FakeEngine(), prefix = '') {
  const store = new MemoryBlockStore();
  let token = 0;
  const service = new DataService({
    store,
    scanner: new FakeScanner(),
    checker: new FakeInputChecker(),
    digest: sha256Hex,
    newToken: () => `token-${++token}`,
    cleanup: Promise.resolve({ state: 'done', removedSessions: 0 }),
  });
  const residueCalls: number[] = [];
  const opened: string[] = [];
  const sources = { added: [] as string[], released: [] as string[] };
  // Counts the reads of original residues, which a loaded run without its original must never make,
  // and keeps the IDs of the runs opened and of the sources added and released in the Data worker.
  const data: DataGateway = Object.assign(Object.create(service) as DataService, {
    addSource: async (file: File) => {
      const source = await service.addSource(file);
      sources.added.push(source.sourceId);
      return source;
    },
    releaseSources: (sourceIds: readonly string[]) => {
      sources.released.push(...sourceIds);
      return service.releaseSources(sourceIds);
    },
    readResidues: (...args: Parameters<DataService['readResidues']>) => {
      residueCalls.push(args[1]);
      return service.readResidues(...args);
    },
    openRun: (runId: string) => {
      opened.push(runId);
      return service.openRun(runId);
    },
  });
  const saved: SavedFile[] = [];
  const downloader = memoryDownloader((file) => saved.push(file));
  let id = 0;
  const newRunId = () => `${prefix}run-${++id}`;
  const coordinator = new Coordinator({ engine, data, downloader, now: () => 1000, newRunId });
  const tray = new CandidateTray({ runs: coordinator.state, data, downloader, now: () => 5000 });
  const session = new Session({
    coordinator,
    tray,
    data,
    compression: browserCompression,
    downloader,
    app: { version: '0.0.0', build: 'test' },
    now: () => new Date(2026, 9, 10, 12, 0, 0).getTime(),
    newRunId,
    // Small ranges and messages, so that every block crosses several of them.
    readBytes: 7,
    sendBytes: 16,
  });
  return { data, service, store, coordinator, tray, session, saved, residueCalls, opened, sources };
}

type World = ReturnType<typeof world>;

function waitAll(world: World, count: number): Promise<void> {
  return new Promise((resolve, reject) => {
    const check = () => {
      const runs = world.coordinator.state.get().runs;
      const ended = runs.filter((run) => (['completed', 'failed', 'cancelled'] as RunStatus[]).includes(run.status));
      if (ended.length >= count) {
        unsubscribe();
        const failed = ended.find((run) => run.status !== 'completed');
        if (failed === undefined) resolve();
        else reject(new Error(`run ${failed.snapshot.number} ${failed.status}: ${failed.record.error}`));
      }
    };
    const unsubscribe = world.coordinator.state.subscribe(check);
    check();
  });
}

/** Two runs, the second's subject made of two files with a record left out; candidates of both, with notes. */
async function searched() {
  const a = world();
  const fileA = await a.data.addSource(new File([FILE_A], 'a.fa'));
  const revisionA = await a.data.reviseDataset((await a.data.indexSource(fileA.sourceId, 1)).revisionId, [1]);
  const fileB = await a.data.addSource(new File([FILE_B], 'b.fa'));
  const revisionB = await a.data.indexSource(fileB.sourceId, 1);
  const blastn: SearchRequest = {
    program: 'blastn',
    query: { text: QUERIES },
    subject: { dataset: { name: 'combined_subject.fa', revisionIds: [revisionA.revisionId, revisionB.revisionId] } },
    parameters: [['-evalue', '1e-5']],
    requestedThreads: 'auto',
    title: 'Joined subjects',
  };
  const tblastn: SearchRequest = {
    program: 'tblastn',
    query: { text: PROTEIN },
    subject: { file: new File([FILE_C], 'c.fa') },
    parameters: [],
    requestedThreads: 2,
  };
  expect((await a.coordinator.enqueue(blastn)).ok).toBe(true);
  expect((await a.coordinator.enqueue(tblastn)).ok).toBe(true);
  await waitAll(a, 2);
  const results = new ResultsBrowser({
    data: a.data,
    describe: (program) => new FakeEngine().describe(program),
    runs: a.coordinator.state,
    verification: { ncbi: '2.17.0', sources: [], programs: {} },
  });
  const [run1, run2] = a.coordinator.state.get().runs;
  for (const [run, rows] of [
    [run1!, [0, 2, 3, 5]],
    [run2!, [0]],
  ] as const) {
    await results.open(run.snapshot.runId);
    const table = await a.data.readHitTable(run.snapshot.runId);
    const ids: HspId[] = rows.map((row) => ({ runId: run.snapshot.runId, qIdx: table.qIdx[row]!, rank: table.rank[row]! }));
    expect(a.tray.add(results.candidateSources(ids))).toMatchObject({ ok: true, added: rows.length });
  }
  const keys = a.tray.state.get().candidates.map((candidate) => candidate.key);
  a.tray.setNote(keys[0]!, 'first <img src=x onerror=alert(1)>');
  a.tray.setNote(keys[4]!, 'protein hit');
  a.tray.move(keys[2]!, 0);
  return a;
}

/** Saves a session and returns the gzip file. */
async function saveSession(a: World, includeCandidates = true): Promise<SavedFile> {
  const result = await a.session.save({ includeCandidates });
  if (!result.ok) throw new Error(result.message);
  return a.saved.at(-1)!;
}

const asFile = (bytes: Uint8Array, name = 'saved.losat-session.gz') => new File([bytes as BlobPart], name);

/** The candidates of a tray as the file restores them: everything but the run IDs and keys. */
const trayView = (w: World) =>
  w.tray.state.get().candidates.map(({ run, index, id, query, subject, coordinates, row, note, addedAt }) => ({
    run: run.number,
    title: run.title,
    program: run.program,
    index,
    hsp: [id.qIdx, id.rank],
    query,
    subject,
    coordinates,
    row,
    note,
    addedAt,
  }));

// --- the container, taken apart and put together again ----------------------------------------------

interface Parsed {
  manifest: SessionManifest;
  blocks: Map<string, Uint8Array>;
  candidates: readonly SessionCandidate[] | undefined;
}

function parse(gzip: Uint8Array): Parsed {
  const reader = new SessionFileReader();
  let manifest: SessionManifest | undefined;
  let candidates: readonly SessionCandidate[] | undefined;
  const parts = new Map<string, number[]>();
  for (const event of reader.push(new Uint8Array(gunzipSync(gzip)))) {
    if (event.type === 'manifest') manifest = event.manifest;
    if (event.type === 'candidates') candidates = event.candidates;
    if (event.type === 'data') {
      const name = runBlockName(event.run, event.stream);
      parts.set(name, [...(parts.get(name) ?? []), ...event.bytes]);
    }
  }
  reader.finish();
  return { manifest: manifest!, blocks: new Map([...parts].map(([name, bytes]) => [name, Uint8Array.from(bytes)])), candidates };
}

/** A gzip session file from parts (the block lengths are the parts' own; the manifest is written as given). */
function build(manifest: unknown, blocks: ReadonlyMap<string, Uint8Array>, candidates?: unknown): Uint8Array {
  const out: Uint8Array[] = [];
  const text = (t: string) => out.push(encoder.encode(t));
  const block = (name: string, bytes: Uint8Array) => {
    text(blockHeader(name, bytes.length));
    out.push(bytes);
    text('\n');
  };
  text(containerHeader());
  block('manifest', encoder.encode(JSON.stringify(manifest)));
  (manifest as SessionManifest).runs.forEach((_, k) => {
    for (const stream of SESSION_STREAMS) block(runBlockName(k + 1, stream), blocks.get(runBlockName(k + 1, stream)) ?? new Uint8Array());
  });
  if (candidates !== undefined) block('candidates', encoder.encode(JSON.stringify(candidates)));
  text(containerEnd());
  return new Uint8Array(gzipSync(Buffer.concat(out)));
}

const copy = <T>(value: T): T => JSON.parse(JSON.stringify(value)) as T;

/** The file with one block replaced, the manifest's lengths of the blocks made to agree with the blocks. */
function rebuilt({ manifest, blocks, candidates }: Parsed, name: string, bytes: Uint8Array): Uint8Array {
  const map = new Map([...blocks, [name, bytes]]);
  const lengths = copy(manifest) as unknown as { runs: Array<{ blocks: Record<string, number> }> };
  lengths.runs.forEach((run, k) => SESSION_STREAMS.forEach((stream) => (run.blocks[stream] = map.get(runBlockName(k + 1, stream))?.length ?? 0)));
  return build(lengths, map, candidates === undefined ? undefined : { candidates });
}

// --- the tests ----------------------------------------------------------------------------------------

describe('Session: saving and opening', () => {
  it('saves two runs with candidates and notes, and opens them in a fresh app without searching', async () => {
    const a = await searched();
    const file = await saveSession(a);
    expect(file.name).toBe('losat-session-20261010-120000.losat-session.gz');
    expect(file.mime).toBe('application/gzip');
    expect(file.blocks).toBeGreaterThanOrEqual(1);
    expect(a.session.state.get().message).toEqual({ text: `Saved ${file.name}: 2 runs and 5 candidates with their notes.`, error: false });

    const { engine, calls } = forbiddenEngine();
    const b = world(engine, 'b-');
    const loaded = await b.session.load(asFile(file.bytes));
    expect(loaded).toEqual({ ok: true, runs: 2, candidates: 5 });
    expect(b.session.state.get().message?.text).toBe(
      'Opened saved.losat-session.gz: 2 runs (runs 1 to 2 here) and 5 candidates with their notes. Nothing was searched again.',
    );
    expect(calls).toEqual([]);

    const runsA = a.coordinator.state.get().runs;
    const runsB = b.coordinator.state.get().runs;
    expect(runsB.map((run) => [run.snapshot.number, run.status, run.fromSession?.number, run.fromSession?.fileName])).toEqual([
      [1, 'completed', 1, 'saved.losat-session.gz'],
      [2, 'completed', 2, 'saved.losat-session.gz'],
    ]);
    for (const [k, runB] of runsB.entries()) {
      const runA = runsA[k]!;
      const [ida, idb] = [runA.snapshot.runId, runB.snapshot.runId];
      expect(idb).not.toBe(ida);
      for (const format of [0, 6, 7] as const) expect(await b.data.readOutput(idb, format)).toEqual(await a.data.readOutput(ida, format));
      expect(await b.data.readHits(idb)).toEqual(await a.data.readHits(ida));
      expect(await b.data.readHitTable(idb)).toEqual(await a.data.readHitTable(ida));
      expect(await b.data.readDiagnostics(idb)).toBe(await a.data.readDiagnostics(ida));
      expect(runB.result).toEqual({ ...runA.result, runId: idb });
      // The snapshot as searched, without the engine bytes or the revisions of the other session.
      const { query: qa, subject: sa, group: ga } = runA.snapshot;
      const { query: qb, subject: sb, group: gb } = runB.snapshot;
      const kept = ({ program, title, argv, requestedThreads, queuedAt }: typeof runA.snapshot) => ({ program, title, argv, requestedThreads, queuedAt });
      expect(kept(runB.snapshot)).toEqual(kept(runA.snapshot));
      expect(runB.snapshot.number).toBe(runA.snapshot.number);
      expect(gb === undefined).toBe(ga === undefined);
      for (const [x, y] of [
        [qa, qb],
        [sa, sb],
      ] as const) {
        expect(y).toEqual({ name: x.name, sha256: x.sha256, revisionIds: [], records: x.records });
        expect(y.bytes).toBeUndefined();
      }
      expect(runB.record).toEqual(runA.record);
    }
    expect(runsB[0]!.fromSession!.inputs.subject.sources).toEqual([
      { name: 'a.fa', size: FILE_A.length, records: 2, excluded: [1] },
      { name: 'b.fa', size: FILE_B.length, records: 1, excluded: [] },
    ]);
    expect(trayView(b)).toEqual(trayView(a));
    expect(b.tray.state.get().selected.size).toBe(5);
    expect(calls).toEqual([]);
  });

  it('writes no storage path, token, run, revision or source ID of the session that saved it', async () => {
    const a = await searched();
    const text = decoder.decode(gunzipSync((await saveSession(a)).bytes));
    const ids = [
      ...a.coordinator.state.get().runs.flatMap((run) => [run.snapshot.runId, ...run.snapshot.query.revisionIds, ...run.snapshot.subject.revisionIds]),
      'token-',
      'runs/',
      'tmp/',
      'groupId',
      'runId',
      'revisionId',
      'sourceId',
    ];
    expect(ids.length).toBeGreaterThan(8);
    for (const id of ids) expect(text).not.toContain(id);
    expect(text.startsWith('LOSAT-WEB-SESSION 1\nmanifest ')).toBe(true);
    expect(text.endsWith('\nLOSAT-WEB-SESSION-END\n')).toBe(true);
  });

  it('opens a file again after its own runs, numbering them after them, and saves a loaded run with the same identity', async () => {
    const a = await searched();
    const file = await saveSession(a);
    const b = world(forbiddenEngine().engine);
    await b.session.load(asFile(file.bytes));
    await b.session.load(asFile(file.bytes, 'again.gz'));
    expect(b.coordinator.state.get().runs.map((run) => [run.snapshot.number, run.fromSession!.number, run.fromSession!.fileName])).toEqual([
      [1, 1, 'saved.losat-session.gz'],
      [2, 2, 'saved.losat-session.gz'],
      [3, 1, 'again.gz'],
      [4, 2, 'again.gz'],
    ]);
    expect(b.tray.state.get().candidates).toHaveLength(10);
    const again = parse((await saveSession(b)).bytes);
    const first = parse(file.bytes);
    expect(again.manifest.runs.map((run) => run.number)).toEqual([1, 2, 3, 4]);
    expect(again.manifest.runs[2]!.query).toEqual(first.manifest.runs[0]!.query);
    expect(again.manifest.runs[3]!.subject).toEqual(first.manifest.runs[1]!.subject);
    expect(again.blocks.get('run3.out0')).toEqual(first.blocks.get('run1.out0'));
  });

  it('leaves the candidates out when asked, and saves only completed runs', async () => {
    const a = await searched();
    expect((await a.coordinator.enqueue({ ...REQUEST_FAILING })).ok).toBe(true);
    await waitAll(a, 3).catch(() => undefined);
    expect(a.coordinator.state.get().runs[2]!.status).toBe('failed');
    const parsed = parse((await saveSession(a, false)).bytes);
    expect(parsed.manifest.candidates).toBe(false);
    expect(parsed.candidates).toBeUndefined();
    expect(parsed.manifest.runs.map((run) => run.number)).toEqual([1, 2]);
    const b = world(forbiddenEngine().engine);
    expect(await b.session.load(asFile((await saveSession(a, false)).bytes))).toEqual({ ok: true, runs: 2, candidates: 0 });
    expect(b.tray.state.get().candidates).toEqual([]);
    expect(await world().session.save()).toEqual({ ok: false, message: expect.stringMatching(/^There are no completed runs to save/) });
  });
});

/** A search that fails in the FakeEngine (its query has no records). */
const REQUEST_FAILING: SearchRequest = { program: 'blastn', query: { text: '' }, subject: { text: '>s\nACGT\n' }, parameters: [], requestedThreads: 1 };

describe('Session: the original FASTA (REQ-23)', () => {
  it('refuses extraction until the original is chosen again and matches; the alignment export needs none', async () => {
    const a = await searched();
    a.tray.selectAll(true);
    const flanks = { role: 'subject' as const, region: { kind: 'flanked' as const, flanks: { left: 2, right: 3 } } };
    expect((await a.tray.extract(flanks)).ok).toBe(true);
    const sequences = a.saved.at(-1)!.bytes;
    expect((await a.tray.exportAlignments()).ok).toBe(true);
    const alignments = a.saved.at(-1)!.bytes;
    const file = await saveSession(a);

    const b = world(forbiddenEngine().engine);
    await b.session.load(asFile(file.bytes));
    const [run1, run2] = b.coordinator.state.get().runs;
    const id1 = run1!.snapshot.runId;
    const message = 'Runs 1 and 2 were loaded from a session file; choose their original subject FASTA in Run details to extract sequences.';
    expect(b.tray.missingOriginals('subject')).toBe(message);
    expect(await b.tray.extract(flanks)).toEqual({ ok: false, message });
    b.tray.select(b.tray.state.get().candidates.filter((c) => c.run.number === 2).map((c) => c.key), false);
    expect(await b.tray.extract(flanks)).toEqual({
      ok: false,
      message: 'Run 1 was loaded from a session file; choose its original subject FASTA in Run details to extract sequences.',
    });
    expect(b.residueCalls).toEqual([]);

    // The aligned rows come from the saved HSP records: no original is needed.
    b.tray.selectAll(true);
    expect((await b.tray.exportAlignments()).ok).toBe(true);
    expect(b.saved.at(-1)!.bytes).toEqual(alignments);

    const refused = async (files: File[], pattern: RegExp) => {
      const added = b.sources.added.length;
      const result = await b.session.attach(id1, 'subject', files);
      expect(result.ok).toBe(false);
      expect(!result.ok && result.message).toMatch(pattern);
      expect(b.coordinator.state.get().runs[0]!.attached).toBeUndefined();
      expect(b.session.state.get().attaching.get(`${id1}/subject`)?.message).toMatch(pattern);
      // Nothing of the attempt stays in the Data worker: its sources and record tables are released (code review L3).
      const attempt = b.sources.added.slice(added);
      expect(b.sources.released).toEqual(expect.arrayContaining(attempt));
      for (const sourceId of attempt) await expect(b.service.previewSource(sourceId, 1)).rejects.toThrow(/unknown source/);
    };
    const fa = (text: string, name = 'a.fa') => new File([text], name);
    // One residue changed: the record's SHA-256 names it.
    await refused([fa(FILE_A.replace('AACCGGTT', 'AACCGGTA')), fa(FILE_B, 'b.fa')], /record 1 \("sA"\) differs from the saved run's record/);
    // One record more in the first file.
    await refused([fa(`${FILE_A}>sX\nAC\n`), fa(FILE_B, 'b.fa')], /"a\.fa" has 3 records, but file 1 of the saved input \("a\.fa"\) had 2/);
    // The right records, but with the exclusion already applied by hand: not the file that was searched.
    await refused([fa('>sA one\nACGTACGTACGTAAACCCGGGTTT\nAACCGGTT\n'), fa(FILE_B, 'b.fa')], /"a\.fa" has 1 record, but file 1 .* had 2/);
    // The right records and exclusions, but a line before the records: the input's SHA-256 differs.
    await refused([fa(FILE_A), fa(`\n${FILE_B}`, 'b.fa')], /the input that they make \(SHA-256 [0-9a-f]{64}\) is not the one that run 1 searched/);
    // One file of two; or two files of which one is the other's copy.
    await refused([fa(FILE_A)], /choose the 2 files that the run joined \("a\.fa", "b\.fa"\), together and in any order; 1 was chosen/);
    await refused([fa(FILE_B, 'b.fa'), fa(FILE_B, 'b2.fa')], /"b2?\.fa" has 1 record, but file 1 of the saved input \("a\.fa"\) had 2/);
    // The query's original is not the subject's.
    await refused([fa(QUERIES, 'query.fa'), fa(FILE_B, 'b.fa')], /record 1 is "sA" .* but the chosen files give "q1"/);
    expect(await b.tray.extract(flanks)).toMatchObject({ ok: false });
    expect(b.residueCalls).toEqual([]);

    // The right files, chosen in the other order, attach the subject of run 1 only, in the order that
    // the run joined them (code review L3); run 2 still has none.
    const kept = b.sources.added.length;
    expect(await b.session.attach(id1, 'subject', [fa(FILE_B, 'b.fa'), fa(FILE_A)])).toEqual({ ok: true });
    expect(b.coordinator.state.get().runs[0]!.attached?.subject?.fileNames).toEqual(['a.fa', 'b.fa']);
    const attachedSources = b.sources.added.slice(kept);
    expect(b.sources.released).not.toEqual(expect.arrayContaining(attachedSources));
    expect(b.session.state.get().attaching.size).toBe(0);
    expect(b.tray.missingOriginals('subject')).toBe(
      'Run 2 was loaded from a session file; choose its original subject FASTA in Run details to extract sequences.',
    );
    expect(b.tray.missingOriginals('query')).toBe(message.replace('subject', 'query'));
    const run2Id = run2!.snapshot.runId;
    expect(await b.session.attach(run2Id, 'subject', [fa(FILE_C, 'renamed.fa')])).toEqual({ ok: true });
    expect(b.tray.missingOriginals('subject')).toBeUndefined();

    expect((await b.tray.extract(flanks)).ok).toBe(true);
    expect(decoder.decode(b.saved.at(-1)!.bytes)).toBe(decoder.decode(sequences));
    expect(b.saved.at(-1)!.bytes).toEqual(sequences);
    expect(b.residueCalls.length).toBeGreaterThan(0);

    // The input FASTA of the run, for its download: the bytes that the engine searched.
    expect((await b.session.attachedInput(id1, 'subject'))?.bytes).toEqual(a.coordinator.state.get().runs[0]!.snapshot.subject.bytes);
    expect(await b.session.attachedInput(id1, 'query')).toBeUndefined();

    // Attached again: the original that it replaces is released from the Data worker.
    expect(await b.session.attach(id1, 'subject', [fa(FILE_A, 'a2.fa'), fa(FILE_B, 'b2.fa')])).toEqual({ ok: true });
    expect(b.coordinator.state.get().runs[0]!.attached?.subject?.fileNames).toEqual(['a2.fa', 'b2.fa']);
    expect(b.sources.released).toEqual(expect.arrayContaining(attachedSources));
    for (const sourceId of attachedSources) await expect(b.service.previewSource(sourceId, 1)).rejects.toThrow(/unknown source/);
    expect((await b.session.attachedInput(id1, 'subject'))?.bytes).toEqual(a.coordinator.state.get().runs[0]!.snapshot.subject.bytes);
  });

  it('attaches only to a run loaded from a session file, and only by an explicit choice', async () => {
    const a = await searched();
    const runId = a.coordinator.state.get().runs[0]!.snapshot.runId;
    expect(await a.session.attach(runId, 'subject', [new File([FILE_A], 'a.fa')])).toEqual({
      ok: false,
      message: 'Only a run loaded from a session file takes its original FASTA again.',
    });
    const b = world(forbiddenEngine().engine);
    // The original files are already in the new session's data layer (as the search form would hold them):
    // nothing is attached until they are chosen for the run.
    const source = await b.data.addSource(new File([FILE_A], 'a.fa'));
    await b.data.indexSource(source.sourceId, 1);
    await b.session.load(asFile((await saveSession(a)).bytes));
    expect(b.coordinator.state.get().runs.every((run) => run.attached === undefined)).toBe(true);
  });
});

describe('Session: refused files leave nothing behind', () => {
  /** Opens a file in a fresh app and checks that it is refused and that nothing stays. */
  async function refusedFile(bytes: Uint8Array, pattern: RegExp, name = 'bad.gz') {
    const { engine, calls } = forbiddenEngine();
    const b = world(engine);
    const result = await b.session.load(asFile(bytes, name));
    expect(result.ok).toBe(false);
    const message = result.ok ? '' : result.message;
    expect(message).toMatch(new RegExp(`^${name.replace('.', '\\.')} was not opened, and nothing was loaded\\. `));
    expect(message).toMatch(pattern);
    expect(b.session.state.get()).toMatchObject({ message: { text: message, error: true } });
    expect(b.session.state.get().busy).toBeUndefined();
    expect(b.coordinator.state.get().runs).toEqual([]);
    expect(b.tray.state.get().candidates).toEqual([]);
    expect(b.store.usage()).toBe(0);
    expect(calls).toEqual([]);
    return message;
  }

  it('refuses a file that is not gzip, damaged gzip data, and a cut file', async () => {
    const a = await searched();
    const gzip = (await saveSession(a)).bytes;
    await refusedFile(encoder.encode('>q1\nACGT\n'), /This is not a LOSAT Web session file: it is not gzip data/, 'query.fa');
    await refusedFile(new Uint8Array(gzipSync(encoder.encode('>q1\nACGT\n'))), /This is not a LOSAT Web session file \(its first line is ">q1"/);
    // A flipped byte inside the compressed data: the data or its checksum no longer agree.
    for (const at of [Math.floor(gzip.length / 3), Math.floor(gzip.length / 2), gzip.length - 6]) {
      const flipped = gzip.slice();
      flipped[at]! ^= 0x40;
      await refusedFile(flipped, /The session file is (damaged|incomplete)|not a LOSAT Web session file|is not valid/);
    }
    await refusedFile(gzip.slice(0, gzip.length - 9), /The session file is (damaged: its gzip data is not valid or is cut short|incomplete)/);
    await refusedFile(gzip.slice(0, Math.floor(gzip.length / 2)), /The session file is (damaged: its gzip data|incomplete)/);
    // Gzip data after the session's own.
    await refusedFile(new Uint8Array([...gzip, ...gzipSync(encoder.encode('x'))]), /damaged/);
  });

  it('refuses HSP records and candidates that the runs cannot have', async () => {
    const a = await searched();
    const { manifest, blocks, candidates } = parse((await saveSession(a)).bytes);
    const hits = decoder.decode(blocks.get('run1.hits')!);
    const withHits = (text: string) => new Map([...blocks, ['run1.hits', encoder.encode(text)]]);
    const lengthsOf = (map: Map<string, Uint8Array>) => {
      const m = copy(manifest) as unknown as { runs: Array<{ blocks: Record<string, number> }> };
      m.runs.forEach((run, k) => SESSION_STREAMS.forEach((stream) => (run.blocks[stream] = map.get(runBlockName(k + 1, stream))?.length ?? 0)));
      return m;
    };
    const cands = { candidates };

    const beyond = withHits(hits.replace('"q_idx":1', `"q_idx":${manifest.runs[0]!.query.records.id.length}`));
    await refusedFile(build(lengthsOf(beyond), beyond, cands), /in the HSP records of run 1 in the file \(run 1 there\), HSP record \d+ names query record 2, but the run has 2 query records/);
    const range = withHits(hits.replace(/"out6":\[(\d+),(\d+)\]/, (_, s) => `"out6":[${s},${manifest.runs[0]!.blocks.out6 + 1}]`));
    await refusedFile(build(lengthsOf(range), range, cands), /has out6 \[\d+,\d+\], not null or a byte range within the run's \d+ bytes of outfmt 6/);
    const notJson = withHits(hits.replace('{"index":0', '{"index":0,'));
    await refusedFile(build(lengthsOf(notJson), notJson, cands), /the HSP records of run 1 in the file \(run 1 there\) cannot be read/);
    const more = copy(manifest) as unknown as { runs: Array<{ hitCount: number }> };
    more.runs[0]!.hitCount += 1;
    await refusedFile(build(more, blocks, cands), /run 1 in the file \(run 1 there\) has \d+ HSP records, but the manifest gives \d+/);

    const wrongRank = { candidates: candidates!.map((c, i) => (i === 0 ? { ...c, rank: c.rank + 1 } : c)) };
    await refusedFile(build(manifest, blocks, wrongRank), /candidates\[0\] is HSP \d+\.\d+ at index \d+ of run \d+ in the file, but the HSP record there is HSP/);
    const noRun = { candidates: [{ ...candidates![0]!, run: 3 }] };
    await refusedFile(build(manifest, blocks, noRun), /The session file's candidates block is not valid: candidates\[0\]\.run is 3, but the file holds 2 runs/);
    const noHsp = { candidates: [{ ...candidates![0]!, index: 999 }] };
    await refusedFile(build(manifest, blocks, noHsp), /candidates\[0\]\.index is 999, beyond the/);
  });

  it('refuses HSP records whose JSON a typed array would have coerced, and leaves the working session as it was (code review M1)', async () => {
    const a = await searched();
    const file = await saveSession(a);
    const parsed = parse(file.bytes);
    // The HSP record of a candidate, in its run's block of HSP records.
    const { run, index } = parsed.candidates![0]!;
    const name = runBlockName(run, 'hits');
    const lines = decoder.decode(parsed.blocks.get(name)!).split('\n');
    const line = lines.findIndex((text) => text !== '' && (JSON.parse(text) as { index: number }).index === index);
    expect(line).toBeGreaterThanOrEqual(0);
    const changed = (edit: (record: Record<string, unknown>) => void) => {
      const record = JSON.parse(lines[line]!) as Record<string, unknown>;
      edit(record);
      return rebuilt(parsed, name, encoder.encode(lines.map((text, k) => (k === line ? JSON.stringify(record) : text)).join('\n')));
    };

    // A working session that already has runs and candidates: a refused file changes none of them.
    const b = world(forbiddenEngine().engine);
    expect((await b.session.load(asFile(file.bytes))).ok).toBe(true);
    const runs = b.coordinator.state.get().runs;
    const { candidates, selected } = b.tray.state.get();
    const usage = b.store.usage();
    const where = `in the HSP records of run ${run} in the file \\(run ${run} there\\), HSP record ${line + 1} `;
    const cases: Array<[(record: Record<string, unknown>) => void, string]> = [
      [(record) => (record.s_idx = null), `${where}has s_idx null, not a whole number of 0 or more`],
      [(record) => (record.query_frame = 259), `${where}has query_frame 259, not null or a frame of -3 to 3 other than 0`],
      [(record) => (record.q_idx = 4294967296), `${where}names query record 4294967296, but the run has \\d+ query records`],
      [(record) => (record.q_start = true), `${where}has q_start true, not a coordinate \\(a whole number of 1 or more\\)`],
      // Valid records, but the candidate's outfmt 6 row is not a row: found while the candidates are made, before the runs join.
      [(record) => (record.out6 = [0, 1]), 'the outfmt 6 row of a candidate is not a row'],
    ];
    for (const [edit, problem] of cases) {
      const before = b.opened.length;
      const result = await b.session.load(asFile(changed(edit), 'bad.gz'));
      expect(result).toEqual({ ok: false, message: expect.stringMatching(new RegExp(`^bad\\.gz was not opened, and nothing was loaded\\. The session file is damaged: ${problem}`)) });
      expect(b.coordinator.state.get().runs).toEqual(runs);
      expect(b.tray.state.get().candidates).toEqual(candidates);
      expect(b.tray.state.get().selected).toEqual(selected);
      expect(b.store.usage()).toBe(usage);
      // The runs that the refused file opened in the Data worker are gone.
      expect(b.opened.length).toBeGreaterThan(before);
      for (const runId of b.opened.slice(before)) await expect(b.service.runBlockLengths(runId)).rejects.toThrow(/no committed result/);
    }
  });

  it('refuses a newer schema, a missing block and a manifest over the limits, and shows hostile text only as data', async () => {
    const a = await searched();
    const { manifest, blocks, candidates } = parse((await saveSession(a)).bytes);
    await refusedFile(build({ ...copy(manifest), schema: 2 }, blocks, { candidates }), /saved by a newer LOSAT Web \(its schema is 2/);
    const fewer = new Map(blocks);
    fewer.set('run2.out7', new Uint8Array());
    await refusedFile(build(manifest, fewer, { candidates }), /block "run2\.out7" has 0 bytes, but the manifest gives \d+/);
    const long = copy(manifest) as unknown as { runs: Array<{ title?: string }> };
    long.runs[0]!.title = 'x'.repeat(10_001);
    await refusedFile(build(long, blocks, { candidates }), /runs\[0\]\.title is longer than 10000 characters/);

    // Names, IDs, titles, argv and notes are data: they are kept as text, whatever they look like.
    const hostile = copy(manifest) as unknown as { runs: Array<{ title?: string; argv: string[]; query: { name: string; sources: Array<{ name: string }> } }> };
    const html = '<img src=x onerror=alert(1)><script>alert(2)</script>javascript:alert(3)';
    hostile.runs[0]!.title = html;
    hostile.runs[0]!.argv[2] = html;
    hostile.runs[0]!.query.name = html;
    hostile.runs[0]!.query.sources[0]!.name = html;
    const notes = { candidates: candidates!.map((c) => ({ ...c, note: html })) };
    const b = world(forbiddenEngine().engine);
    expect(await b.session.load(asFile(build(hostile, blocks, notes), `${html}.gz`))).toMatchObject({ ok: true });
    const run = b.coordinator.state.get().runs[0]!;
    expect(run.snapshot.title).toBe(html);
    expect(run.snapshot.query.name).toBe(html);
    expect(run.fromSession!.fileName).toBe(`${html}.gz`);
    expect(b.tray.state.get().candidates.every((c) => c.note === html)).toBe(true);
  });
});
