import { gzipSync, zstdCompressSync } from 'node:zlib';
import { describe, expect, it } from 'vitest';
import { sha256Hex } from 'losat';
import type { InputRole, RegistryFile } from 'losat/domain';
import type { EnginePrepareHint, StartResult, ValidationResult } from 'losat/ports';
import { Coordinator, type EnqueueAllResult, type SearchRequest } from '../../src/application/coordinator';
import { conffileName, SearchDraft, type DraftSource } from '../../src/application/draft';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeAccessionFetcher } from '../../src/infra/fake/fake-fetcher';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import { NodeInflater } from './support/node-inflater';

function setup(
  options: {
    validate?: (argv: readonly string[], registry?: RegistryFile) => ValidationResult;
    /** The engine's start of an argv; by default none of its files and no dry run. */
    start?: (argv: readonly string[]) => StartResult;
    useCoordinator?: boolean;
    checkInput?: DataService['checkInput'];
    fetcher?: FakeAccessionFetcher;
  } = {},
) {
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
  const checks: string[] = [];
  const checkInput = data.checkInput.bind(data);
  data.checkInput = async (program, role, revisionIds) => {
    checks.push(`${program}|${role}|${revisionIds.join(',')}`);
    return (options.checkInput ?? checkInput)(program, role, revisionIds);
  };
  const fake = new FakeEngine();
  const validated: Array<readonly string[]> = [];
  /** The registry file that each validation was given (undefined for none). */
  const registries: Array<RegistryFile | undefined> = [];
  const prepared: EnginePrepareHint[] = [];
  const engine = {
    describe: fake.describe.bind(fake),
    start: async (argv: readonly string[]): Promise<StartResult> =>
      options.start?.(argv) ?? { ok: true, start: { conffile: null, logFile: null, dryRun: false } },
    validate: async (argv: readonly string[], registry?: RegistryFile): Promise<ValidationResult> => {
      validated.push(argv);
      registries.push(registry);
      return options.validate?.(argv, registry) ?? { ok: true };
    },
    prepare: (hint: EnginePrepareHint) => void prepared.push(hint),
  };
  const dropped: string[] = [];
  const dropSource = data.dropSource.bind(data);
  data.dropSource = async (sourceId) => {
    dropped.push(sourceId);
    return dropSource(sourceId);
  };
  let id = 0;
  const coordinator = new Coordinator({
    engine: { ...engine, run: fake.run.bind(fake), cancel: fake.cancel.bind(fake) },
    data,
    downloader: { save: () => undefined },
    now: () => 0,
    newRunId: () => `run-${++id}`,
    digest: sha256Hex,
  });
  const requests: SearchRequest[][] = [];
  const draft = new SearchDraft({
    engine,
    data,
    enqueueAll: async (batch): Promise<EnqueueAllResult> => {
      requests.push([...batch]);
      return options.useCoordinator ? coordinator.enqueueAll(batch) : { ok: true, runIds: batch.map((_, i) => `r${i}`) };
    },
    fetcher: options.fetcher ?? new FakeAccessionFetcher(),
    now: () => 0,
    debounceMs: 0,
  });
  /** What the data layer holds: its sources, and the sources of its revisions. */
  const held = () => {
    const inner = data as unknown as { sources: Map<string, File>; revisions: Map<string, { sourceId: string }> };
    return {
      sources: [...inner.sources.keys()].sort(),
      revisionSources: [...new Set([...inner.revisions.values()].map((revision) => revision.sourceId))].sort(),
    };
  };
  return { draft, data, requests, validated, registries, prepared, coordinator, checks, dropped, held };
}

const file = (name: string, text: string) => new File([text], name, { type: 'text/plain' });
const sourceOf = (draft: SearchDraft, role: InputRole, index = 0): DraftSource => draft.state.get()[role].sources[index]!;
const datasetOf = (request: SearchRequest, role: InputRole) => {
  const input = request[role];
  if (!('dataset' in input)) throw new Error('expected a dataset input');
  return input.dataset;
};

/** Resolves when the run has completed. */
function completed(coordinator: Coordinator, runId: string): Promise<void> {
  return new Promise((resolve) => {
    const check = () => {
      if (coordinator.state.get().runs.find((run) => run.snapshot.runId === runId)?.status !== 'completed') return;
      unsubscribe();
      resolve();
    };
    const unsubscribe = coordinator.state.subscribe(check);
    check();
  });
}

describe('SearchDraft inputs', () => {
  it('indexes pasted text and checks it with the engine', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>q1\nACGT\n>q2\nGGCC\n');
    await draft.idle();
    const source = sourceOf(draft, 'query');
    expect(source).toMatchObject({ origin: 'paste', name: 'query.fa', status: 'ready', excluded: [] });
    expect(source.base?.records.map((r) => r.id)).toEqual(['q1', 'q2']);
    expect(source.check).toEqual({ state: 'ok', records: 2 });
    expect(source.head).toBe('>q1\nACGT\n>q2\nGGCC\n');
  });

  it('replaces the paste source when the text changes, and removes it when the text is emptied', async () => {
    const { draft } = setup();
    draft.setPaste('subject', '>s1\nACGT\n');
    await draft.idle();
    draft.setPaste('subject', '>s2\nACGT\n');
    await draft.idle();
    expect(draft.state.get().subject.sources).toHaveLength(1);
    expect(sourceOf(draft, 'subject').base?.records[0]?.id).toBe('s2');
    draft.setPaste('subject', '');
    await draft.idle();
    expect(draft.state.get().subject.sources).toHaveLength(0);
  });

  it('reads a text without a defline as one record without a title, and adds a defline on request', async () => {
    const { draft } = setup();
    draft.setPaste('query', 'ACGTACGT\n');
    await draft.idle();
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'ready', parser: 1 });
    expect(sourceOf(draft, 'query').base?.records.map((r) => [r.id, r.header_offset, r.sequence_offset, r.length])).toEqual([['', 0, 0, 8]]);
    draft.addDefline('query');
    await draft.idle();
    expect(draft.state.get().query.paste).toBe('>pasted_query\nACGTACGT\n');
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'ready' });
  });

  it('names a byte order mark of a text that cannot be read, and removes it on request (WA-IN-05)', async () => {
    const { draft } = setup();
    draft.setPaste('query', '\uFEFF>q1\nACGT\n');
    await draft.idle();
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'failed', lead: { mark: 'UTF-8' } });
    expect(sourceOf(draft, 'query').error).toContain("doesn't look like plausible data");
    draft.removeByteOrderMark('query');
    await draft.idle();
    expect(draft.state.get().query.paste).toBe('>q1\nACGT\n');
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'ready', check: { state: 'ok', records: 1 } });
    expect(sourceOf(draft, 'query').lead).toBeUndefined();
  });

  it('reads a GenBank file as FASTA text, as the command line does, and names its format through a mode change (WA-IN-04)', async () => {
    const { draft } = setup();
    draft.addFiles('query', [file('x.gb', 'LOCUS       X 8 bp DNA\nORIGIN\n        1 acgtacgt\n//\n')]);
    await draft.idle();
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'ready', lead: { format: 'GenBank' } });
    expect(sourceOf(draft, 'query').base?.records).toHaveLength(1);
    draft.setProgram('blastp');
    await draft.idle();
    expect(sourceOf(draft, 'query')).toMatchObject({ parser: 2, lead: { format: 'GenBank' } });
    draft.setProgram('blastn');
    await draft.idle();
    expect(sourceOf(draft, 'query')).toMatchObject({ parser: 1, status: 'ready', lead: { format: 'GenBank' } });
  });

  it('names the format of a compressed file from its inflated bytes', async () => {
    const { draft } = setup();
    draft.addFiles('query', [new File([gzipSync('LOCUS       X 8 bp DNA\nORIGIN\n        1 acgtacgt\n//\n')], 'x.gb.gz')]);
    await draft.idle();
    expect(sourceOf(draft, 'query')).toMatchObject({ name: 'x.gb', status: 'ready', lead: { format: 'GenBank' } });
  });

  it('shows the engine refusal of a record and maps it to the record table; excluding it passes the check', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>a\nACGT\n>b\nACGT!ACGT\n>c\nACGT\n');
    await draft.idle();
    const refused = sourceOf(draft, 'query').check;
    expect(refused).toMatchObject({ state: 'refused', record: 1 });
    expect(refused?.state === 'refused' && refused.message).toContain('query record 2 (b)');
    draft.setIncluded('query', 'paste', [1], false);
    await draft.idle();
    expect(sourceOf(draft, 'query').check).toEqual({ state: 'ok', records: 2 });
  });

  it('maps a refused record among the included records only', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>a\nACGT\n>b\nACGT\n>c\nACGT!ACGT\n');
    await draft.idle();
    draft.setIncluded('query', 'paste', [0], false);
    await draft.idle();
    // The engine read b and c (record 2 of its input is c, index 2 in the table).
    expect(sourceOf(draft, 'query').check).toMatchObject({ state: 'refused', record: 2 });
  });

  it('checks the inputs again for another program, and reads them again when it reads the role with another kind', async () => {
    const { draft, checks } = setup();
    draft.setPaste('query', '>a\nMKLV\n');
    draft.setPaste('subject', '>s\n$$$$\n');
    await draft.idle();
    const revision = sourceOf(draft, 'query').base!.revisionId;
    expect(sourceOf(draft, 'query').base!.records[0]!.length).toBe(3);
    expect(sourceOf(draft, 'subject')).toMatchObject({ status: 'failed', parser: 1 });
    // BLASTX reads its query as BLASTN does (kind 1): the query is checked again, not read again.
    // Its subject is protein (kind 2): the subject is read again, the failed one too.
    draft.setProgram('blastx');
    await draft.idle();
    expect(sourceOf(draft, 'query').base!.revisionId).toBe(revision);
    expect(sourceOf(draft, 'query').check).toEqual({ state: 'ok', records: 1 });
    expect(checks.at(-1)).toMatch(/^blastx\|query\|/);
    expect(sourceOf(draft, 'subject')).toMatchObject({ status: 'failed', parser: 2 });
    // BLASTP reads its query with the protein flags (kind 2): read again.
    draft.setProgram('blastp');
    await draft.idle();
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'ready', parser: 2 });
    expect(sourceOf(draft, 'query').base!.revisionId).not.toBe(revision);
    expect(sourceOf(draft, 'query').base!.records[0]!.length).toBe(4);
  });

  it('keeps the records left out and the region of a source through a mode that reads it with the other kind (WA-IN-01)', async () => {
    const { draft, data } = setup();
    const added: string[] = [];
    const indexed: string[] = [];
    const addSource = data.addSource.bind(data);
    const indexSource = data.indexSource.bind(data);
    data.addSource = async (f) => {
      added.push(f.name);
      return addSource(f);
    };
    data.indexSource = async (sourceId, parser) => {
      indexed.push(`${sourceId}|${parser}`);
      return indexSource(sourceId, parser);
    };
    draft.setPaste('query', '>a\nACGTACGTAC\n>b\nACGTACGTAC\n>c\nACGTACGTAC\n');
    draft.addFiles('subject', [file('s.fa', '>s1\nACGTACGTACGTACGT\n>s2\nACGTACGTAC\n')]);
    await draft.idle();
    draft.setIncluded('query', 'paste', [1], false);
    draft.setIncluded('query', 'paste', [2], false);
    await draft.idle();
    draft.setRegion('query', { start: '2', stop: '8' });
    draft.setIncluded('subject', sourceOf(draft, 'subject').key, [1], false);
    await draft.idle();
    const revision = sourceOf(draft, 'query').base!.revisionId;
    expect(draft.parameters()).toEqual([['-query_loc', '2-8']]);
    const readOnce = indexed.length;

    draft.setProgram('blastp');
    await draft.idle();
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'ready', parser: 2, excluded: [] });
    expect(sourceOf(draft, 'subject')).toMatchObject({ status: 'ready', parser: 2, excluded: [] });
    // The other kind's table is not the region's record.
    expect(draft.region('query')).toBeUndefined();

    draft.setProgram('blastn');
    await draft.idle();
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'ready', parser: 1, excluded: [1, 2] });
    expect(sourceOf(draft, 'query').base!.revisionId).toBe(revision);
    expect(sourceOf(draft, 'query').check).toEqual({ state: 'ok', records: 1 });
    expect(sourceOf(draft, 'subject')).toMatchObject({ status: 'ready', parser: 1, excluded: [1] });
    expect(draft.region('query')).toEqual({ start: '2', stop: '8' });
    expect(draft.parameters()).toEqual([['-query_loc', '2-8']]);
    // Each source is added once and read once for each kind: the way back reads nothing.
    expect([...added].sort()).toEqual(['query.fa', 's.fa']);
    expect(indexed.length).toBe(readOnce + 2);

    draft.setProgram('blastp');
    draft.setProgram('blastn');
    await draft.idle();
    expect(indexed.length).toBe(readOnce + 2);
    expect(sourceOf(draft, 'query').excluded).toEqual([1, 2]);
  });

  it('keeps the exclusions of each kind apart, and takes them back when the exclusions of the other kind change', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>a\nMKLV\n>b\nMKLV\n>c\nMKLV\n');
    await draft.idle();
    draft.setIncluded('query', 'paste', [0], false);
    draft.setProgram('blastp');
    await draft.idle();
    draft.setIncluded('query', 'paste', [2], false);
    await draft.idle();
    draft.setProgram('blastn');
    await draft.idle();
    expect(sourceOf(draft, 'query').excluded).toEqual([0]);
    draft.setProgram('blastp');
    await draft.idle();
    expect(sourceOf(draft, 'query').excluded).toEqual([2]);
    draft.setIncluded('query', 'paste', [1], false);
    draft.setProgram('blastn');
    await draft.idle();
    expect(sourceOf(draft, 'query').excluded).toEqual([0]);
    draft.setProgram('blastp');
    await draft.idle();
    expect(sourceOf(draft, 'query').excluded).toEqual([1, 2]);
  });

  it('does not let the read of a kind that was left replace the table that the way back restored', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>a\nMKLV\n>b\nMKLV\n>c\nMKLV\n');
    await draft.idle();
    const revision = sourceOf(draft, 'query').base!.revisionId;
    draft.setIncluded('query', 'paste', [0], false);
    // BLASTP starts to read the source, and BLASTN is chosen again before that read ends.
    draft.setProgram('blastp');
    draft.setProgram('blastn');
    await draft.idle();
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'ready', parser: 1, excluded: [0] });
    expect(sourceOf(draft, 'query').base!.revisionId).toBe(revision);
    expect(sourceOf(draft, 'query').base!.parser).toBe(1);
  });

  it('keeps each source of a role through the round trip, and a source removed on the way stays removed', async () => {
    const { draft, requests } = setup();
    draft.setPaste('query', '>q\nACGT\n>r\nACGT\n');
    draft.addFiles('query', [file('a.fa', '>a1\nACGT\n>a2\nACGT\n'), file('b.fa', '>b1\nACGT\n>b2\nACGT\n>b3\nACGT\n')]);
    draft.addFiles('subject', [file('s.fa', '>s\nACGTACGTACGT\n')]);
    await draft.idle();
    const [paste, a, b] = draft.state.get().query.sources;
    draft.setIncluded('query', paste!.key, [0], false);
    draft.setIncluded('query', a!.key, [1], false);
    draft.setIncluded('query', b!.key, [0, 2], false);
    await draft.idle();
    const revisions = draft.state.get().query.sources.map((source) => source.base!.revisionId);

    draft.setProgram('blastp');
    await draft.idle();
    expect(draft.state.get().query.sources.map((source) => source.excluded)).toEqual([[], [], []]);
    draft.removeSource('query', a!.key);
    await draft.idle();

    draft.setProgram('blastn');
    await draft.idle();
    const back = draft.state.get().query.sources;
    expect(back.map((source) => [source.key, source.status, source.excluded])).toEqual([
      [paste!.key, 'ready', [0]],
      [b!.key, 'ready', [0, 2]],
    ]);
    expect(back.map((source) => source.base!.revisionId)).toEqual([revisions[0], revisions[2]]);
    // The sources are searched together, as they were left.
    expect(await draft.submit()).toMatchObject({ ok: true });
    expect(datasetOf(requests[0]![0]!, 'query').revisionIds).toHaveLength(2);
  });

  it('keeps a region through a mode that reads the role with the other kind, also when a region was drawn there', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>a\nACGTACGTACGTAC\n>b\nACGTACGTACGTAC\n>c\nACGTACGTACGTAC\n');
    await draft.idle();
    draft.setIncluded('query', 'paste', [1, 2], false);
    await draft.idle();
    draft.setRegion('query', { start: '2', stop: '12' });

    draft.setProgram('blastp');
    await draft.idle();
    expect(draft.region('query')).toBeUndefined();
    // One record there as well (another table, so another record): its own region.
    draft.setIncluded('query', 'paste', [1, 2], false);
    await draft.idle();
    draft.setRegion('query', { start: '3', stop: '9' });
    expect(draft.parameters()).toEqual([['-query_loc', '3-9']]);

    draft.setProgram('blastn');
    await draft.idle();
    expect(draft.region('query')).toEqual({ start: '2', stop: '12' });
    expect(draft.parameters()).toEqual([['-query_loc', '2-12']]);
    draft.setProgram('blastp');
    await draft.idle();
    expect(draft.region('query')).toEqual({ start: '3', stop: '9' });
    // Cleared there, it stays set in the other kind.
    draft.setRegion('query', undefined);
    draft.setProgram('blastn');
    await draft.idle();
    expect(draft.region('query')).toEqual({ start: '2', stop: '12' });
  });

  it('reads a text pasted just before a program change once, with the new kind', async () => {
    const { draft, data } = setup();
    const kinds: number[] = [];
    const indexSource = data.indexSource.bind(data);
    data.indexSource = async (sourceId, parser) => {
      kinds.push(parser);
      return indexSource(sourceId, parser);
    };
    draft.setPaste('query', '>a\nMKLV\n>b\nMKLV\n');
    draft.setProgram('blastp');
    await draft.idle();
    expect(kinds).toEqual([2]);
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'ready', parser: 2 });
  });
});

describe('SearchDraft gives the data layer back the sources it no longer shows (SURVEY #15)', () => {
  /** The data layer's sources of the form's cards (a source that cannot be read keeps none). */
  const shown = (draft: SearchDraft) =>
    (['query', 'subject'] as const)
      .flatMap((role) => draft.state.get()[role].sources.flatMap((source) => (source.sourceId === undefined ? [] : [source.sourceId])))
      .sort();
  const internals = (draft: SearchDraft) =>
    draft as unknown as { verdicts: Map<string, unknown>; revisions: Map<string, unknown>; sourceRevisions: Map<string, unknown> };

  it('drops a removed file, a cleared paste and a replaced paste, with their record tables and the draft caches', async () => {
    const { draft, held, dropped } = setup();
    draft.setPaste('query', '>q1\nACGT\n>q2\nACGT\n');
    draft.addFiles('subject', [file('a.fa', '>a\nACGT\n'), file('b.fa', '>b\nACGT\n')]);
    await draft.idle();
    draft.setIncluded('query', 'paste', [1], false);
    await draft.idle();
    expect(held()).toEqual({ sources: shown(draft), revisionSources: shown(draft) });
    expect(internals(draft).verdicts.size).toBeGreaterThan(0);
    const [paste, a] = [sourceOf(draft, 'query').sourceId!, sourceOf(draft, 'subject', 0).sourceId!];

    draft.removeSource('subject', 'file-1');
    draft.setPaste('query', '>q3\nACGT\n');
    await draft.idle();
    expect(dropped).toEqual([a, paste]);
    expect(held()).toEqual({ sources: shown(draft), revisionSources: shown(draft) });

    draft.setPaste('query', '');
    draft.removeSource('subject', 'file-2');
    await draft.idle();
    expect(held()).toEqual({ sources: [], revisionSources: [] });
    expect([internals(draft).verdicts.size, internals(draft).revisions.size, internals(draft).sourceRevisions.size]).toEqual([0, 0, 0]);
  });

  it('drops the sources of indexes that a newer text superseded, and of a source that cannot be read', async () => {
    const { draft, held } = setup({});
    // Each pause in typing indexes the text; the older indexes are still running.
    draft.setPaste('query', '>q\nA\n');
    await draft.idle();
    for (const text of ['>q\nAC\n', '>q\nACG\n', '>q\nACGT\n']) draft.setPaste('query', text);
    await draft.idle();
    draft.addFiles('subject', [file('bad.fa', '\uFEFF>s\nACGT\n')]);
    await draft.idle();
    expect(sourceOf(draft, 'subject')).toMatchObject({ status: 'failed' });
    expect(held()).toEqual({ sources: shown(draft), revisionSources: shown(draft) });
    expect(shown(draft)).toHaveLength(1);
  });

  it('keeps a source through modes that read it with the other kind, both of its tables', async () => {
    const { draft, held, dropped } = setup();
    draft.setPaste('query', '>q\nACGT\n');
    await draft.idle();
    const sourceId = sourceOf(draft, 'query').sourceId;
    draft.setProgram('blastp');
    await draft.idle();
    draft.setProgram('blastn');
    await draft.idle();
    expect(sourceOf(draft, 'query')).toMatchObject({ sourceId, status: 'ready', parser: 1 });
    expect(dropped).toEqual([]);
    expect(held().sources).toEqual([sourceId]);
  });

  it('lets a run keep its inputs in its snapshot: the sources go with the form, also one removed while the run is started', async () => {
    const { draft, held, dropped, coordinator, data } = setup({ useCoordinator: true });
    draft.setPaste('query', '>q\nACGT\n');
    draft.addFiles('subject', [file('s.fa', '>s\nACGTACGT\n')]);
    await draft.idle();
    const query = sourceOf(draft, 'query').sourceId!;
    const subject = sourceOf(draft, 'subject').sourceId!;
    // Removed while the search is started, between the snapshots of the query and of the subject.
    const build = data.buildRunInput.bind(data);
    let calls = 0;
    data.buildRunInput = async (revisionIds) => {
      if (++calls === 1) draft.removeSource('subject', 'file-1');
      return build(revisionIds);
    };
    expect(await draft.submit()).toMatchObject({ ok: true });
    expect(calls).toBe(2);
    // "Run LOSAT" held the removed subject until the run had its input.
    expect(dropped).toEqual([subject]);
    draft.setPaste('query', '');
    await draft.idle();
    expect(dropped).toEqual([subject, query]);
    expect(held()).toEqual({ sources: [], revisionSources: [] });
    const run = coordinator.state.get().runs[0]!;
    expect(new TextDecoder().decode(run.snapshot.subject.bytes)).toBe('>s\nACGTACGT\n');
    expect(run.snapshot.subject.records).toEqual([{ id: 's', length: 8 }]);
    await completed(coordinator, run.snapshot.runId);
    expect(await coordinator.readOutput(run.snapshot.runId, 6)).toContain('q\ts\t');
  });

  it('reads a file again after its source was dropped', async () => {
    const { draft, held } = setup();
    const same = file('s.fa', '>s\nACGT\n');
    draft.addFiles('subject', [same]);
    await draft.idle();
    draft.removeSource('subject', 'file-1');
    await draft.idle();
    expect(held().sources).toEqual([]);
    draft.addFiles('subject', [same]);
    await draft.idle();
    expect(sourceOf(draft, 'subject')).toMatchObject({ status: 'ready', check: { state: 'ok', records: 1 } });
    expect(held().sources).toEqual(shown(draft));
  });

  it('drops the sources that an Example replaces', async () => {
    const { draft, held } = setup();
    draft.addFiles('subject', [file('s.fa', '>s\nACGT\n')]);
    draft.setPaste('query', '>q\nACGT\n');
    await draft.idle();
    draft.useExample('refseq-16s-bacteria' as never, '>e\nACGT\n', file('bacteria.fna', '>b\nACGT\n'));
    await draft.idle();
    expect(held()).toEqual({ sources: shown(draft), revisionSources: shown(draft) });
    expect(shown(draft)).toHaveLength(2);
  });
});

describe('SearchDraft accessions (line WF-1)', () => {
  const NR = '>NR_074540.1 Bacillus cereus 16S\nACGTACGT\n';
  const AB = '>AB000001.1 Rhizoctonia solani\nGGCCAATT\n';

  it('keeps a box of accessions without a source until it is fetched, then reads the FASTA written into it', async () => {
    const fetcher = new FakeAccessionFetcher({ 'NR_074540.1': NR, 'AB000001.1': AB });
    const { draft } = setup({ fetcher });
    draft.setFetchSource('ncbi');
    draft.setPaste('query', 'nr_074540.1\n\nAB000001.1\n');
    await draft.idle();
    expect(draft.state.get().query.sources).toHaveLength(0);
    expect(draft.accessions('query')).toEqual(['NR_074540.1', 'AB000001.1']);
    expect(draft.readiness()).toBe('Fetch the query accessions first.');
    expect(fetcher.requests).toHaveLength(0);
    draft.fetchAccessions('query');
    await draft.idle();
    expect(fetcher.requests.map((r) => [r.source, r.kind, r.accession])).toEqual([
      ['ncbi', 'nucleotide', 'NR_074540.1'],
      ['ncbi', 'nucleotide', 'AB000001.1'],
    ]);
    expect(draft.state.get().query.paste).toBe(NR + AB);
    expect(draft.state.get().query.fetch).toBeUndefined();
    expect(sourceOf(draft, 'query').base?.records.map((r) => r.id)).toEqual(['NR_074540.1', 'AB000001.1']);
  });

  it('asks the chosen source for the kind the program reads, and keeps the choice', async () => {
    const fetcher = new FakeAccessionFetcher({ P12345: '>sp|P12345|AATM_RABIT\nMALLHS\n' });
    const { draft } = setup({ fetcher });
    draft.setProgram('blastp');
    draft.setFetchSource('ebi');
    draft.setPaste('subject', 'P12345');
    draft.fetchAccessions('subject');
    await draft.idle();
    expect(fetcher.requests.map((r) => [r.source, r.kind, r.accession])).toEqual([['ebi', 'protein', 'P12345']]);
    expect(draft.state.get().fetchSource).toBe('ebi');
    expect(draft.state.get().subject.paste).toBe('>sp|P12345|AATM_RABIT\nMALLHS\n');
  });

  it('says why a fetch failed, leaves the box as it was, and asks EBI nothing it cannot serve', async () => {
    const fetcher = new FakeAccessionFetcher({ 'AB000001.1': AB, 'ZZ000001.1': { ok: false, reason: 'refused', status: 429 } });
    const { draft } = setup({ fetcher });
    draft.setFetchSource('ncbi');
    draft.setPaste('query', 'AB000001.1\nAB999999.1');
    draft.fetchAccessions('query');
    await draft.idle();
    expect(draft.state.get().query.fetch).toEqual({ state: 'failed', message: 'NCBI has no nucleotide record AB999999.1.' });
    expect(draft.state.get().query.paste).toBe('AB000001.1\nAB999999.1');
    draft.setFetchSource('ebi');
    expect(draft.state.get().query.fetch).toBeUndefined();
    draft.setPaste('query', 'NR_074540.1');
    draft.fetchAccessions('query');
    await draft.idle();
    expect(draft.state.get().query.fetch).toEqual({
      state: 'failed',
      message: 'EBI (ENA) has no RefSeq records: fetch NR_074540.1 from NCBI or CNCB.',
    });
    draft.setFetchSource('cncb');
    draft.setPaste('query', 'ZZ000001.1');
    draft.fetchAccessions('query');
    await draft.idle();
    expect(draft.state.get().query.fetch).toEqual({
      state: 'failed',
      message: 'CNCB asks to wait (too many requests): fetch again in a moment.',
    });
    expect(fetcher.requests.map((r) => [r.source, r.accession])).toEqual([
      ['ncbi', 'AB000001.1'],
      ['ncbi', 'AB999999.1'],
      ['cncb', 'ZZ000001.1'],
    ]);
  });

  it('asks EBI nothing when one accession of the box is not for it', async () => {
    const fetcher = new FakeAccessionFetcher({ 'AB000001.1': AB });
    const { draft } = setup({ fetcher });
    draft.setFetchSource('ebi');
    draft.setPaste('query', 'AB000001.1\nNR_074540.1');
    draft.fetchAccessions('query');
    await draft.idle();
    expect(fetcher.requests).toHaveLength(0);
    expect(draft.state.get().query.fetch).toMatchObject({ state: 'failed', message: expect.stringContaining('NR_074540.1') });
  });

  it('runs no search of accessions: Run neither fetches nor waits for a running fetch', async () => {
    const fetcher = new FakeAccessionFetcher({ 'NR_074540.1': NR });
    const { draft, requests } = setup({ fetcher });
    draft.setFetchSource('ncbi');
    draft.setPaste('subject', '>s\nACGT\n');
    fetcher.hold();
    draft.setPaste('query', 'NR_074540.1');
    draft.fetchAccessions('query');
    expect(await draft.submit()).toEqual({ ok: false, message: 'The query accession is still being fetched.' });
    fetcher.release();
    await draft.idle();
    expect(requests).toHaveLength(0);
    expect(draft.state.get().query.paste).toBe(NR);
  });

  it('drops a running fetch when the box is edited, and limits the fetched text', async () => {
    const fetcher = new FakeAccessionFetcher({ 'NR_074540.1': NR, 'NC_000001.11': `>NC_000001.11\n${'A'.repeat(10_000_001)}\n` });
    const { draft } = setup({ fetcher });
    draft.setFetchSource('ncbi');
    fetcher.hold();
    draft.setPaste('query', 'NR_074540.1');
    draft.fetchAccessions('query');
    await Promise.resolve();
    expect(draft.state.get().query.fetch).toEqual({ state: 'fetching', source: 'ncbi', done: 0, total: 1 });
    draft.setPaste('query', '>mine\nACGT\n');
    fetcher.release();
    await draft.idle();
    expect(draft.state.get().query.paste).toBe('>mine\nACGT\n');
    expect(draft.state.get().query.fetch).toBeUndefined();
    draft.setPaste('query', 'NC_000001.11');
    draft.fetchAccessions('query');
    await draft.idle();
    expect(draft.state.get().query.fetch).toMatchObject({ state: 'failed', message: expect.stringContaining('larger than 10 MB') });
  });
});

describe('SearchDraft accessions from the fastest source (line WF-3)', () => {
  const AB = '>AB000001.1 Rhizoctonia solani\nGGCCAATT\n';
  const NR = '>NR_074540.1 Bacillus cereus 16S\nACGTACGT\n';
  const PROBE = '>MH011443.1 Homo sapiens TP53 (TP53) gene, exon 5 and partial cds\nTGGGTTGATT\n';
  const nucleotideProbes = { 'MH011443.1': PROBE };
  const asked = (fetcher: FakeAccessionFetcher) => fetcher.requests.map((r) => [r.source, r.accession, r.fresh ?? false]);

  it('races the sources with a fixed record on the first fetch, sends the accessions to the first that answers, and keeps the order', async () => {
    const fetcher = new FakeAccessionFetcher({ ...nucleotideProbes, 'AB000001.1': AB, 'NR_074540.1': NR });
    const { draft } = setup({ fetcher });
    expect(draft.state.get().fetchSource).toBe('fastest');
    fetcher.hold();
    draft.setPaste('query', 'AB000001.1');
    draft.fetchAccessions('query');
    await Promise.resolve();
    expect(draft.state.get().query.fetch).toEqual({ state: 'racing' });
    expect(draft.readiness()).toBe('The query accession is still being fetched.');
    // The race asks for its own record, never the user's accession, and not from the browser's cache.
    expect(asked(fetcher)).toEqual([
      ['ncbi', 'MH011443.1', true],
      ['ebi', 'MH011443.1', true],
      ['cncb', 'MH011443.1', true],
    ]);
    fetcher.release('ebi');
    await draft.idle();
    expect(draft.state.get().query.paste).toBe(AB);
    expect(asked(fetcher).slice(3)).toEqual([['ebi', 'AB000001.1', false]]);
    expect(draft.state.get().race.nucleotide?.answered.map((entry) => entry.source)).toEqual(['ebi']);
    fetcher.release('cncb');
    fetcher.release('ncbi');
    await draft.raceIdle();
    expect(draft.state.get().race.nucleotide).toEqual({
      answered: [
        { source: 'ebi', ms: 0 },
        { source: 'cncb', ms: 0 },
        { source: 'ncbi', ms: 0 },
      ],
      failed: [],
    });
    // The order stays: no race again, and a RefSeq accession goes to the first source that has it.
    draft.setPaste('subject', 'AB000001.1');
    draft.fetchAccessions('subject');
    draft.setPaste('query', 'NR_074540.1');
    draft.fetchAccessions('query');
    await draft.idle();
    expect(asked(fetcher).slice(4)).toEqual([
      ['ebi', 'AB000001.1', false],
      ['cncb', 'NR_074540.1', false],
    ]);
    expect(draft.state.get().query.paste).toBe(NR);
  });

  it('races each kind once, with its own record, and keeps a chosen source', async () => {
    const fetcher = new FakeAccessionFetcher({
      'ncbi:NP_414542.1': '>NP_414542.1 thr operon leader peptide\nMKRISTTITTTITITTGNGAG\n',
      'ebi:P0AD86': '>sp|P0AD86|LPT_ECOLI thr operon leader peptide\nMKRISTTITTTITITTGNGAG\n',
      P12345: '>sp|P12345|AATM_RABIT\nMALLHS\n',
    });
    const { draft } = setup({ fetcher });
    draft.setProgram('blastp');
    draft.setPaste('subject', 'P12345');
    draft.fetchAccessions('subject');
    await draft.idle();
    await draft.raceIdle();
    expect(asked(fetcher)).toEqual([
      ['ncbi', 'NP_414542.1', true],
      ['ebi', 'P0AD86', true],
      ['cncb', 'NP_414542.1', true],
      ['ncbi', 'P12345', false],
    ]);
    expect(draft.state.get().race.protein).toEqual({
      answered: [
        { source: 'ncbi', ms: 0 },
        { source: 'ebi', ms: 0 },
      ],
      failed: ['cncb'],
    });
    expect(draft.state.get().race.nucleotide).toBeUndefined();
    draft.setFetchSource('cncb');
    draft.setPaste('subject', 'WP_000001.1');
    draft.fetchAccessions('subject');
    await draft.idle();
    expect(asked(fetcher).slice(4)).toEqual([['cncb', 'WP_000001.1', false]]);
    expect(draft.state.get().fetchSource).toBe('cncb');
  });

  it('goes to NCBI when no source that serves the accessions answered, and races again after a race without answers', async () => {
    const fetcher = new FakeAccessionFetcher({
      'ebi:MH011443.1': PROBE,
      'NR_074540.1': NR,
      'AB000001.1': AB,
    });
    const { draft } = setup({ fetcher });
    fetcher.hold('ebi');
    // EBI has no RefSeq records: its answer would not count for this box.
    draft.setPaste('query', 'NR_074540.1');
    draft.fetchAccessions('query');
    await draft.idle();
    expect(asked(fetcher).slice(3)).toEqual([['ncbi', 'NR_074540.1', false]]);
    expect(draft.state.get().race.nucleotide).toEqual({ answered: [], failed: ['ncbi', 'cncb'] });
    fetcher.release('ebi');
    await draft.raceIdle();
    expect(draft.state.get().race.nucleotide?.answered.map((entry) => entry.source)).toEqual(['ebi']);

    const offline = new FakeAccessionFetcher({ 'AB000001.1': AB });
    const second = setup({ fetcher: offline });
    second.draft.setPaste('query', 'AB000001.1');
    second.draft.fetchAccessions('query');
    await second.draft.idle();
    await second.draft.raceIdle();
    expect(second.draft.state.get().race.nucleotide).toBeUndefined();
    expect(second.draft.state.get().query.paste).toBe(AB);
    second.draft.setPaste('subject', 'AB000001.1');
    second.draft.fetchAccessions('subject');
    await second.draft.idle();
    expect(offline.requests.filter((r) => r.fresh === true)).toHaveLength(6);
  });

  it('sends no accession to a source that lacks it: UniProt not to CNCB, GenBase only to CNCB', async () => {
    const fetcher = new FakeAccessionFetcher({
      'NP_414542.1': '>NP_414542.1 thr operon leader peptide\nMKRISTTITTTITITTGNGAG\n',
      'ebi:P0AD86': '>sp|P0AD86|LPT_ECOLI thr operon leader peptide\nMKRISTTITTTITITTGNGAG\n',
      P12345: '>sp|P12345|AATM_RABIT\nMALLHS\n',
      'C_AA004835.1': '>C_AA004835.1 SARS-CoV-2\nACGT\n',
    });
    const { draft } = setup({ fetcher });
    draft.setProgram('blastp');
    fetcher.hold('ncbi', 'ebi');
    draft.setPaste('subject', 'P12345');
    draft.fetchAccessions('subject');
    // CNCB answers first, but it has no UniProt records: the fetch waits for NCBI or EBI.
    await new Promise((resolve) => setTimeout(resolve, 0));
    expect(draft.state.get().race.protein?.answered.map((entry) => entry.source)).toEqual(['cncb']);
    expect(draft.state.get().subject.fetch).toEqual({ state: 'racing' });
    fetcher.release('ebi');
    await draft.idle();
    expect(asked(fetcher).slice(3)).toEqual([['ebi', 'P12345', false]]);
    fetcher.release();
    await draft.raceIdle();

    // Only CNCB has GenBase records: no race, and EBI is not asked.
    draft.setProgram('blastn');
    draft.setPaste('query', 'C_AA004835.1');
    draft.fetchAccessions('query');
    await draft.idle();
    expect(asked(fetcher).slice(4)).toEqual([['cncb', 'C_AA004835.1', false]]);
    expect(draft.state.get().race.nucleotide).toBeUndefined();
    draft.setFetchSource('ebi');
    draft.setPaste('query', 'C_AA004835.1');
    draft.fetchAccessions('query');
    await draft.idle();
    expect(draft.state.get().query.fetch).toEqual({
      state: 'failed',
      message: 'EBI (ENA) has no GenBase records: fetch C_AA004835.1 from CNCB.',
    });

    // No one site has both.
    draft.setFetchSource('fastest');
    draft.setProgram('blastp');
    draft.setPaste('subject', 'C_AA004835.1\nP12345');
    draft.fetchAccessions('subject');
    await draft.idle();
    expect(draft.state.get().subject.fetch).toEqual({
      state: 'failed',
      message: 'No one site has all of these accessions: C_AA004835.1 only from CNCB, P12345 only from NCBI or EBI.',
    });
    expect(asked(fetcher)).toHaveLength(5);
  });

  it('lets the race finish when the box is edited during it, and sends nothing of the box', async () => {
    const fetcher = new FakeAccessionFetcher({ ...nucleotideProbes, 'AB000001.1': AB });
    const { draft } = setup({ fetcher });
    fetcher.hold();
    draft.setPaste('query', 'AB000001.1');
    draft.fetchAccessions('query');
    await Promise.resolve();
    draft.setPaste('query', '>mine\nACGT\n');
    expect(draft.state.get().query.fetch).toBeUndefined();
    fetcher.release();
    await draft.idle();
    await draft.raceIdle();
    expect(asked(fetcher).map(([, accession]) => accession)).toEqual(['MH011443.1', 'MH011443.1', 'MH011443.1']);
    expect(draft.state.get().query.paste).toBe('>mine\nACGT\n');
    expect(draft.state.get().race.nucleotide?.answered).toHaveLength(3);
  });
});

describe('SearchDraft accessions by groups, one site a press (line WF-5)', () => {
  const C1 = '>C_AA004835.1 SARS-CoV-2\nMESLVPGFNE\n';
  const C2 = '>C_AA000002.1 other\nMKTAYIAKQR\n';
  const P = '>sp|P12345|AATM_RABIT\nMALLHS\n';
  const records = { 'C_AA004835.1': C1, 'C_AA000002.1': C2, P12345: P };
  const asked = (fetcher: FakeAccessionFetcher) => fetcher.requests.map((r) => [r.source, r.accession]);

  it('sends each group only to the site pressed, holds what came, and writes the box in its order once all have come', async () => {
    const fetcher = new FakeAccessionFetcher(records);
    const { draft } = setup({ fetcher });
    draft.setProgram('blastp');
    draft.setPaste('query', '>q\nMALLHS\n');
    draft.setPaste('subject', 'C_AA004835.1\nP12345\nC_AA000002.1');
    draft.fetchGroup('subject', 0, 'cncb');
    await draft.idle();
    expect(asked(fetcher)).toEqual([
      ['cncb', 'C_AA004835.1'],
      ['cncb', 'C_AA000002.1'],
    ]);
    expect(draft.state.get().subject.fetched).toEqual({ 0: { fasta: C1, source: 'cncb' }, 2: { fasta: C2, source: 'cncb' } });
    expect(draft.state.get().subject.paste).toBe('C_AA004835.1\nP12345\nC_AA000002.1');
    expect(draft.readiness()).toBe('Fetch the subject accessions first.');
    // A site that does not serve the group is not asked.
    draft.fetchGroup('subject', 1, 'cncb');
    await draft.idle();
    expect(fetcher.requests).toHaveLength(2);
    fetcher.hold('ebi');
    draft.fetchGroup('subject', 1, 'ebi');
    await Promise.resolve();
    expect(draft.state.get().subject.fetch).toEqual({ state: 'fetching', source: 'ebi', done: 0, total: 1 });
    expect(draft.readiness()).toBe('The subject accessions are still being fetched.');
    fetcher.release();
    await draft.idle();
    await draft.raceIdle();
    expect(asked(fetcher).slice(2)).toEqual([['ebi', 'P12345']]);
    expect(draft.state.get().subject.paste).toBe(C1 + P + C2);
    expect(draft.state.get().subject.fetched).toBeUndefined();
    expect(draft.state.get().subject.fetch).toBeUndefined();
    expect(sourceOf(draft, 'subject').base?.records.map((r) => r.id)).toEqual(['C_AA004835.1', 'sp|P12345|AATM_RABIT', 'C_AA000002.1']);
    // No race: a group goes to the site pressed.
    expect(fetcher.requests.filter((r) => r.fresh === true)).toHaveLength(0);
    expect(draft.state.get().race).toEqual({});
  });

  it('names the site and the accession of a failed group, keeps the groups that came, and fetches that group again alone', async () => {
    const fetcher = new FakeAccessionFetcher({
      ...records,
      'ncbi:P12345': { ok: false, reason: 'unreachable', detail: 'offline' },
    });
    const { draft } = setup({ fetcher });
    draft.setProgram('blastp');
    draft.setPaste('query', 'C_AA004835.1\nP12345');
    draft.fetchGroup('query', 0, 'cncb');
    await draft.idle();
    draft.fetchGroup('query', 1, 'ncbi');
    await draft.idle();
    expect(draft.state.get().query.fetch).toEqual({ state: 'failed', message: 'P12345: NCBI could not be reached (offline).' });
    expect(draft.state.get().query.fetched).toEqual({ 0: { fasta: C1, source: 'cncb' } });
    // Another choice clears the message, not what came.
    draft.clearFetchFailure('query');
    expect(draft.state.get().query.fetch).toBeUndefined();
    expect(draft.state.get().query.fetched).toEqual({ 0: { fasta: C1, source: 'cncb' } });
    draft.fetchGroup('query', 1, 'ebi');
    await draft.idle();
    expect(asked(fetcher)).toEqual([
      ['cncb', 'C_AA004835.1'],
      ['ncbi', 'P12345'],
      ['ebi', 'P12345'],
    ]);
    expect(draft.state.get().query.paste).toBe(C1 + P);
  });

  it('drops what came when the box is edited or the mode reads the other kind, and counts the whole box against the limit', async () => {
    const fetcher = new FakeAccessionFetcher({ ...records, 'C_AA000003.1': `>C_AA000003.1\n${'M'.repeat(10_000_000 - C1.length)}\n` });
    const { draft } = setup({ fetcher });
    draft.setProgram('blastp');
    draft.setPaste('subject', 'C_AA004835.1\nP12345');
    draft.fetchGroup('subject', 0, 'cncb');
    await draft.idle();
    draft.setPaste('subject', 'C_AA004835.1\nP12345\n');
    expect(draft.state.get().subject.fetched).toBeUndefined();
    draft.fetchGroup('subject', 0, 'cncb');
    await draft.idle();
    expect(draft.state.get().subject.fetched).toEqual({ 0: { fasta: C1, source: 'cncb' } });
    // BLASTN's subject is nucleotide: there CNCB serves both, and the box has no groups.
    draft.setProgram('blastn');
    expect(draft.state.get().subject.fetched).toBeUndefined();
    draft.fetchGroup('subject', 0, 'cncb');
    await draft.idle();
    expect(fetcher.requests).toHaveLength(2);

    draft.setProgram('blastp');
    draft.setPaste('subject', 'C_AA004835.1\nP12345\nC_AA000003.1');
    draft.fetchGroup('subject', 1, 'ncbi');
    await draft.idle();
    draft.fetchGroup('subject', 0, 'cncb');
    await draft.idle();
    expect(draft.state.get().subject.fetch).toMatchObject({ state: 'failed', message: expect.stringContaining('larger than 10 MB') });
    expect(draft.state.get().subject.fetched).toEqual({ 1: { fasta: P, source: 'ncbi' } });
  });

  it('writes nothing of a group whose fetch an edit or a mode of the other kind stopped', async () => {
    const fetcher = new FakeAccessionFetcher(records);
    const { draft } = setup({ fetcher });
    draft.setProgram('blastp');
    draft.setPaste('subject', 'C_AA004835.1\nP12345');
    fetcher.hold();
    draft.fetchGroup('subject', 0, 'cncb');
    await Promise.resolve();
    draft.setPaste('subject', 'C_AA004835.1\nP12345\n');
    fetcher.release();
    await draft.idle();
    expect(draft.state.get().subject.fetched).toBeUndefined();
    expect(draft.state.get().subject.fetch).toBeUndefined();
    expect(draft.state.get().subject.paste).toBe('C_AA004835.1\nP12345\n');
    fetcher.hold();
    draft.fetchGroup('subject', 1, 'ncbi');
    await Promise.resolve();
    draft.setProgram('blastn');
    fetcher.release();
    await draft.idle();
    expect(draft.state.get().subject.fetched).toBeUndefined();
    expect(draft.state.get().subject.fetch).toBeUndefined();
  });

  it('fetches a box that one site serves whole with the one Fetch, not by groups', async () => {
    const fetcher = new FakeAccessionFetcher(records);
    const { draft } = setup({ fetcher });
    draft.setProgram('blastp');
    draft.setPaste('subject', 'C_AA004835.1\nC_AA000002.1');
    draft.fetchGroup('subject', 0, 'cncb');
    await draft.idle();
    expect(fetcher.requests).toHaveLength(0);
  });
});

describe('SearchDraft submit', () => {
  it('refuses to queue without inputs, with every record excluded, or with a refused input', async () => {
    const { draft, requests } = setup();
    expect(await draft.submit()).toEqual({ ok: false, message: 'Add the query sequences: paste them or open FASTA files.' });
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s\nACGT!ACGT\n');
    await draft.idle();
    const refused = await draft.submit();
    expect(refused.ok).toBe(false);
    expect(!refused.ok && refused.message).toMatch(/^Subject \(pasted\): subject record 1 \(s\) has/);
    draft.setIncluded('subject', 'paste', [0], false);
    expect(await draft.submit()).toEqual({ ok: false, message: 'Every subject record is excluded.' });
    expect(requests).toHaveLength(0);
    expect(draft.state.get().message).toEqual({ kind: 'error', text: 'Every subject record is excluded.' });
  });

  it("queues BLASTX with the regions of its nucleotide query and protein subject (DW-39)", async () => {
    const { draft, requests } = setup();
    draft.setPaste('query', '>q\nACGTACGTAC\n');
    draft.setPaste('subject', '>s\nMKLVMKLV\n');
    await draft.idle();
    draft.setProgram('blastx');
    await draft.idle();
    expect(draft.regionRecord('query')?.id).toBe('q');
    expect(draft.regionRecord('subject')?.id).toBe('s');
    draft.setRegion('query', { start: '2', stop: '8' });
    draft.setRegion('subject', { start: '3', stop: '8' });
    expect(draft.parameters()).toEqual([
      ['-query_loc', '2-8'],
      ['-subject_loc', '3-8'],
    ]);
    // The subject's region is checked against its residues, the query's against its bases.
    draft.setRegion('subject', { start: '3', stop: '9' });
    expect(await draft.submit()).toEqual({
      ok: false,
      message: 'Subject region: The stop must be between 1 and 8, the length of the record.',
    });
    draft.setRegion('subject', undefined);
    draft.setRegion('query', { start: '2', stop: '10' });
    expect(await draft.submit()).toMatchObject({ ok: true });
    expect(requests.at(-1)![0]!).toMatchObject({ program: 'blastx', parameters: [['-query_loc', '2-10']] });
  });

  it('opens a file named .gz under its plain name and inflated size', async () => {
    const { draft } = setup();
    draft.addFiles('subject', [new File([gzipSync('>a\nACGT\n')], 'x.fa.gz')]);
    await draft.idle();
    expect(sourceOf(draft, 'subject')).toMatchObject({ status: 'ready', name: 'x.fa', size: 8 });
    expect(draft.inputName('subject')).toBe('x.fa');
  });

  it('runs an empty query, also from an empty .gz or .zst file, as the command line does', async () => {
    for (const empty of [file('e.fa', ''), new File([gzipSync('')], 'e.fa.gz'), new File([zstdCompressSync('')], 'e.fa.zst')]) {
      const { draft, requests } = setup();
      draft.addFiles('query', [empty]);
      draft.addFiles('subject', [file('s.fa', '>s\nACGTACGT\n')]);
      await draft.idle();
      expect(sourceOf(draft, 'query')).toMatchObject({ status: 'ready', size: 0, check: { state: 'ok', records: 0 } });
      expect(await draft.submit()).toMatchObject({ ok: true });
      expect(requests).toHaveLength(1);
    }
  });

  it('opens a file named .zst under its plain name and inflated size', async () => {
    const { draft } = setup();
    draft.addFiles('subject', [new File([zstdCompressSync('>a\nACGT\n')], 'x.fa.zst')]);
    await draft.idle();
    expect(sourceOf(draft, 'subject')).toMatchObject({ status: 'ready', name: 'x.fa', size: 8 });
    expect(draft.inputName('subject')).toBe('x.fa');
  });

  it('queues one combined search of all files, named as plan §5.3 says', async () => {
    const { draft, requests } = setup();
    draft.setPaste('query', '>q\nACGT\n');
    draft.addFiles('subject', [file('a.fa', '>a\nACGT\n'), file('b.fa', '>b\nACGT')]);
    await draft.idle();
    expect(draft.runCount()).toBe(1);
    expect(await draft.submit()).toMatchObject({ ok: true });
    expect(requests).toHaveLength(1);
    const [request] = requests[0]!;
    expect(datasetOf(request!, 'query').name).toBe('query.fa');
    expect(datasetOf(request!, 'subject')).toEqual({
      name: 'combined_subject.fa',
      revisionIds: [sourceOf(draft, 'subject', 0).base!.revisionId, sourceOf(draft, 'subject', 1).base!.revisionId],
    });
    expect(draft.state.get().message).toBeUndefined();
  });

  it('queues separate searches as every query input against every subject input', async () => {
    const { draft, requests } = setup();
    draft.addFiles('query', [file('q1.fa', '>q1\nACGT\n'), file('q2.fa', '>q2\nACGT\n')]);
    draft.addFiles('subject', [file('s1.fa', '>s1\nACGT\n'), file('s2.fa', '>s2\nACGT\n'), file('s3.fa', '>s3\nACGT\n')]);
    await draft.idle();
    draft.setMode('subject', 'separate');
    expect(draft.runCount()).toBe(3);
    draft.setMode('query', 'separate');
    expect(draft.runCount()).toBe(6);
    await draft.submit();
    const names = requests[0]!.map((r) => `${datasetOf(r, 'query').name}|${datasetOf(r, 'subject').name}`);
    expect(names).toEqual(['q1.fa|s1.fa', 'q1.fa|s2.fa', 'q1.fa|s3.fa', 'q2.fa|s1.fa', 'q2.fa|s2.fa', 'q2.fa|s3.fa']);
    expect(draft.state.get().message).toBeUndefined();
  });

  it('uses the same revision for the same selection, so that a kept subject is the same input (R1)', async () => {
    const { draft, requests } = setup();
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s1\nACGT\n>s2\nACGT\n>s3\nACGT\n');
    await draft.idle();
    const base = sourceOf(draft, 'subject').base!.revisionId;
    await draft.submit();
    draft.setPaste('query', '>other\nGGGG\n');
    await draft.submit();
    draft.setIncluded('subject', 'paste', [1], false);
    await draft.submit();
    await draft.submit();
    const subjects = requests.map((batch) => datasetOf(batch[0]!, 'subject').revisionIds[0]);
    expect(subjects[0]).toBe(base);
    expect(subjects[1]).toBe(base);
    expect(subjects[2]).not.toBe(base);
    expect(subjects[3]).toBe(subjects[2]);
    const queries = requests.map((batch) => datasetOf(batch[0]!, 'query').revisionIds[0]);
    expect(queries[1]).not.toBe(queries[0]);
  });

  it('writes the region of a role that has one record, and only then', async () => {
    const { draft, requests } = setup();
    draft.setPaste('query', '>q1\nACGTACGTAC\n>q2\nACGT\n');
    draft.setPaste('subject', '>s\nACGTACGTACGTACGTACGT\n');
    await draft.idle();
    expect(draft.regionRecord('query')).toBeUndefined();
    expect(draft.regionRecord('subject')?.id).toBe('s');
    // A role of two records has no region.
    draft.setRegion('query', { start: '2', stop: '5' });
    expect(draft.region('query')).toBeUndefined();
    draft.setRegion('subject', { start: ' 3', stop: '012 ' });
    expect(draft.parameters()).toEqual([['-subject_loc', '3-12']]);
    draft.setIncluded('query', 'paste', [1], false);
    await draft.idle();
    draft.setRegion('query', { start: '2', stop: '5' });
    expect(draft.parameters()).toEqual([
      ['-query_loc', '2-5'],
      ['-subject_loc', '3-12'],
    ]);
    draft.setRegion('subject', { start: '3', stop: '21' });
    const result = await draft.submit();
    expect(result).toEqual({
      ok: false,
      message: 'Subject region: The stop must be between 1 and 20, the length of the record.',
    });
    draft.setRegion('subject', undefined);
    await draft.submit();
    expect(requests[0]![0]!.parameters).toEqual([['-query_loc', '2-5']]);
  });

  it('writes only the values the user set and that differ from the engine default', async () => {
    const { draft, requests } = setup();
    draft.setProgram('tblastx');
    await draft.idle();
    draft.setField('-evalue', '10'); // the engine's default
    draft.setField('-word_size', '');
    draft.setField('-seg', 'no');
    draft.setField('-max_target_seqs', ' 50 ');
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s\nACGT\n');
    await draft.submit();
    expect(requests[0]![0]!.parameters).toEqual([
      ['-max_target_seqs', '50'],
      ['-seg', 'no'],
    ]);
    // Each program keeps its own values.
    draft.setProgram('blastn');
    draft.setField('-task', 'dc-megablast');
    draft.setField('-template_type', 'optimal');
    draft.setField('-template_length', '21');
    draft.setField('-lcase_masking', true);
    // The argv follows the sections of the page: Filters and Masking, then Discontiguous Word Options.
    expect(draft.parameters()).toEqual([
      ['-task', 'dc-megablast'],
      ['-lcase_masking', true],
      ['-template_length', '21'],
      ['-template_type', 'optimal'],
    ]);
    draft.setField('-task', 'blastn');
    expect(draft.parameters()).toEqual([
      ['-task', 'blastn'],
      ['-lcase_masking', true],
    ]);
    draft.setProgram('tblastx');
    expect(draft.parameters()).toEqual([
      ['-max_target_seqs', '50'],
      ['-seg', 'no'],
    ]);
  });

  it('restores the defaults of some fields: their values are cleared, the others kept, and the argv validated again', async () => {
    const { draft, validated } = setup();
    draft.setField('-task', 'dc-megablast');
    draft.setField('-evalue', '1e-5');
    draft.setField('-lcase_masking', true);
    draft.setField('-template_length', '18');
    draft.setProgram('tblastx');
    draft.setField('-evalue', '5');
    draft.setProgram('blastn');
    await draft.idle();
    draft.resetFields(['-evalue', '-lcase_masking', '-template_length', '-word_size']);
    expect(draft.state.get().values.blastn).toEqual({ '-task': 'dc-megablast' });
    expect(Object.isFrozen(draft.state.get().values.blastn)).toBe(true);
    // Another program's values stay.
    expect(draft.state.get().values.tblastx).toEqual({ '-evalue': '5' });
    await draft.idle();
    expect(validated.at(-1)).toEqual(['blastn', '-query', 'query.fa', '-subject', 'subject.fa', '-task', 'dc-megablast']);
  });

  it('writes the words of Additional options after the form\'s, for each program, and refuses an unclosed quote', async () => {
    const { draft, requests, validated } = setup();
    draft.setField('-evalue', '1e-5');
    draft.setExtra(" -ungapped  -strand 'plus' ");
    await draft.idle();
    expect(validated.at(-1)).toEqual([
      'blastn', '-query', 'query.fa', '-subject', 'subject.fa', '-evalue', '1e-5', '-ungapped', '-strand', 'plus',
    ]);
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s\nACGT\n');
    expect((await draft.submit()).ok).toBe(true);
    expect(requests[0]![0]!.parameters).toEqual([['-evalue', '1e-5']]);
    expect(requests[0]![0]!.extra).toEqual(['-ungapped', '-strand', 'plus']);
    // Each program keeps its own words; a program without them adds none.
    draft.setProgram('blastp');
    await draft.idle();
    expect(draft.state.get().extra).toMatchObject({ blastn: " -ungapped  -strand 'plus' ", blastp: '' });
    expect(validated.at(-1)).toEqual(['blastp', '-query', 'query.fa', '-subject', 'subject.fa']);
    expect((await draft.submit()).ok).toBe(true);
    expect(requests[1]![0]!).not.toHaveProperty('extra');
    // An unclosed quote: the argv is not made, and nothing is queued.
    draft.setExtra("-seg '12 2.2");
    await draft.idle();
    expect(draft.state.get().validation).toEqual({ state: 'invalid', message: "Additional options: a quote (') is not closed" });
    expect(await draft.submit()).toEqual({ ok: false, message: "Additional options: a quote (') is not closed" });
    expect(requests).toHaveLength(2);
  });

  it('runs the words of Additional options in the coordinator; it refuses the options that LOSAT manages', async () => {
    const { draft, coordinator } = setup({ useCoordinator: true });
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s\nACGT\n');
    draft.setExtra('-num_threads 4');
    await draft.idle();
    expect(draft.state.get().validation).toEqual({ state: 'invalid', message: '-num_threads is set in the form: Threads' });
    expect(await draft.submit()).toEqual({ ok: false, message: '-num_threads is set in the form: Threads' });
    draft.setExtra('-qcov_hsp_perc 50 -ungapped');
    expect((await draft.submit()).ok).toBe(true);
    expect(coordinator.state.get().runs[0]!.snapshot.argv.slice(5)).toEqual(['-qcov_hsp_perc', '50', '-ungapped']);
  });

  it('refuses in Additional options the options that the form sets, and the subranges, with the control that sets them', async () => {
    const { draft } = setup();
    await draft.idle();
    const refusal = async (text: string) => {
      draft.setExtra(text);
      await draft.idle();
      const validation = draft.state.get().validation;
      return validation.state === 'invalid' ? validation.message : validation.state;
    };
    // A -task here would change the defaults that the form leaves out of the argv.
    expect(await refusal('-ungapped -task=blastn')).toBe('-task is set in the form: Optimize for');
    expect(await refusal('-evalue 1')).toBe('-evalue is set in the form: Expect threshold');
    expect(await refusal('-penalty -3')).toBe('-penalty is set in the form: Match/Mismatch Scores');
    expect(await refusal('-gapopen 5')).toBe('-gapopen is set in the form: Gap Costs');
    expect(await refusal("-dust '20 64 1'")).toBe('-dust is set in the form: Low complexity regions');
    expect(await refusal('-culling_limit 2')).toBe('-culling_limit is set in the form: Max matches in a query range');
    expect(await refusal('-query_loc 1-10')).toBe('-query_loc is set in the form: Query subrange');
    expect(await refusal('-subject_loc=1-10')).toBe('-subject_loc is set in the form: Subject subrange');
    // Another program's form: BLASTP has no -dust field, TBLASTX's genetic codes are fields.
    draft.setProgram('tblastx');
    expect(await refusal('-db_gencode 4')).toBe('-db_gencode is set in the form: Genetic code');
    draft.setProgram('blastp');
    expect(await refusal('-qcov_hsp_perc 50 -use_sw_tback')).toBe('ok');
  });

  // WA-FM-04, WA-FM-05: the CLI prints its help for -h, and its messages for these name
  // -subject and -num_threads, which the browser sets itself.
  it("refuses in Additional options the CLI's help and the options that a search in the browser cannot have, in the browser's words", async () => {
    const { draft } = setup({ useCoordinator: true });
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s\nACGT\n');
    await draft.idle();
    const refusal = async (text: string) => {
      draft.setExtra(text);
      await draft.idle();
      const validation = draft.state.get().validation;
      return validation.state === 'invalid' ? validation.message : validation.state;
    };
    expect(await refusal('-help')).toBe('-help is not available in the browser');
    expect(await refusal('-ungapped -h')).toBe('-h is not available in the browser');
    expect(await refusal('-db nt')).toBe('-db (a database search) is not supported in the browser yet');
    expect(await refusal('-taxids 9606')).toBe('-taxids (a database search filter) is not supported in the browser yet');
    expect(await refusal('-mt_mode 1')).toBe('-mt_mode is managed by LOSAT');
    expect(await draft.submit()).toEqual({ ok: false, message: '-mt_mode is managed by LOSAT' });
  });

  it('gives every run of the draft its Job Title, outside the argv', async () => {
    const { draft, requests, coordinator } = setup({ useCoordinator: true });
    draft.addFiles('subject', [file('s1.fa', '>s1\nACGT\n'), file('s2.fa', '>s2\nACGT\n')]);
    draft.setPaste('query', '>q\nACGT\n');
    await draft.idle();
    draft.setMode('subject', 'separate');
    expect(draft.state.get().title).toBe('');
    draft.setTitle(' Two plasmids ');
    expect(draft.state.get().title).toBe(' Two plasmids ');
    await draft.submit();
    expect(requests[0]!.map((r) => r.title)).toEqual([' Two plasmids ', ' Two plasmids ']);
    expect(coordinator.state.get().runs.map((run) => run.snapshot.title)).toEqual(['Two plasmids', 'Two plasmids']);
    expect(coordinator.state.get().runs[0]!.snapshot.argv).not.toContain('Two plasmids');
    // An empty title is no title.
    draft.setTitle('   ');
    await draft.submit();
    expect(requests[1]!.every((r) => !('title' in r))).toBe(true);
  });

  it('asks the engine to prepare once the argv is valid and both roles have inputs, with the size and argv of the run', async () => {
    const { draft, prepared, validated } = setup();
    draft.setPaste('query', '>q1\nACGT\n');
    await draft.idle();
    expect(prepared).toEqual([]);
    draft.setPaste('subject', '>s1\nACGTACGT\n');
    await draft.idle();
    const argv = validated.at(-1);
    expect(argv?.slice(0, 5)).toEqual(['blastn', '-query', 'query.fa', '-subject', 'subject.fa']);
    expect(prepared.at(-1)).toEqual({ requestedThreads: 'auto', inputBytes: 9 + 13, argv });
    draft.setThreads(8);
    draft.setThreads(4);
    await draft.idle();
    // Stepping through the choices hints only the last one.
    expect(prepared.filter((hint) => hint.requestedThreads !== 'auto')).toEqual([{ requestedThreads: 4, inputBytes: 22, argv }]);
  });

  it('does not ask the engine to prepare for an argv that it refuses', async () => {
    const { draft, prepared } = setup({ validate: () => ({ ok: false, message: 'no' }) });
    draft.setPaste('query', '>q1\nACGT\n');
    draft.setPaste('subject', '>s1\nACGT\n');
    await draft.idle();
    draft.setThreads(2);
    await draft.idle();
    expect(prepared).toEqual([]);
  });

  it("validates the draft's argv with the engine and shows its message", async () => {
    const { draft, validated } = setup({
      validate: (argv) => (argv.includes('-word_size') ? { ok: false, message: 'BLAST query/options error: bad' } : { ok: true }),
    });
    await draft.idle();
    expect(draft.state.get().validation).toEqual({ state: 'ok' });
    draft.setField('-word_size', '3');
    await draft.idle();
    expect(draft.state.get().validation).toEqual({ state: 'invalid', message: 'BLAST query/options error: bad' });
    expect(validated.at(-1)).toEqual(['blastn', '-query', 'query.fa', '-subject', 'subject.fa', '-word_size', '3']);
  });

  // NCBI's application start (line RG E1w-3b): the engine's start of the argv names the file of
  // -conffile; the form asks for it, and the validation and the run read the chosen file's bytes under
  // the argv's name. Without a file, the engine's message stands.
  it("asks for the file of the argv's -conffile, validates with its bytes under the argv's name, and runs with it", async () => {
    const missing = 'Critical: (106.23) Registry file "settings.ini" not found, aborting.';
    const { draft, registries, requests } = setup({
      start: (argv) => ({ ok: true, start: { conffile: argv.includes('-conffile') ? 'settings.ini' : null, logFile: null, dryRun: false } }),
      validate: (argv, registry) => (argv.includes('-conffile') && registry === undefined ? { ok: false, message: missing } : { ok: true }),
    });
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s\nACGT\n');
    await draft.idle();
    expect(draft.state.get().start).toEqual({ conffile: null, logFile: null, dryRun: false });
    draft.setExtra('-conffile settings.ini');
    await draft.idle();
    expect(draft.state.get().start?.conffile).toBe('settings.ini');
    expect(draft.state.get().validation).toEqual({ state: 'invalid', message: missing });
    expect(registries.at(-1)).toBeUndefined();
    // A file of any name is read as the argv's. The page keeps the File, not its bytes: the validation
    // gets the File (the Data worker reads it on an instance of its own), and so does the run's request.
    const chosen = file('my.ini', '[BLAST]\n');
    let reads = 0;
    const read = chosen.arrayBuffer.bind(chosen);
    chosen.arrayBuffer = () => (reads++, read());
    draft.setConffile(chosen);
    await draft.idle();
    expect(draft.state.get().conffiles.blastn).toBe(chosen);
    expect(draft.state.get().validation).toEqual({ state: 'ok' });
    expect(registries.at(-1)).toEqual({ name: 'settings.ini', file: chosen });
    expect(registries.at(-1)!.file).toBe(chosen);
    draft.setField('-word_size', '16');
    await draft.idle();
    expect(registries.at(-1)!.file).toBe(chosen);
    expect((await draft.submit()).ok).toBe(true);
    expect(requests.at(-1)![0]!.conffile).toBe(chosen);
    expect(reads).toBe(0);
    // Another program has its own Additional options and file.
    draft.setProgram('blastp');
    await draft.idle();
    expect(draft.state.get().start?.conffile).toBeNull();
    expect(registries.at(-1)).toBeUndefined();
    draft.setProgram('blastn');
    await draft.idle();
    expect(registries.at(-1)).toEqual({ name: 'settings.ini', file: chosen });
    // Without -conffile the chosen file is not read, nor sent with the run.
    draft.setExtra('');
    await draft.idle();
    expect(registries.at(-1)).toBeUndefined();
    expect((await draft.submit()).ok).toBe(true);
    expect('conffile' in requests.at(-1)![0]!).toBe(false);
    // Let go, the file is not read again.
    draft.setExtra('-conffile settings.ini');
    draft.setConffile(undefined);
    await draft.idle();
    expect(draft.state.get().conffiles.blastn).toBeUndefined();
    expect(draft.state.get().validation).toEqual({ state: 'invalid', message: missing });
  });

  // The row of -conffile follows the argv: it goes when the argv cannot be made or the program changes
  // (until the engine's start of the new argv), and `-conffile=` names no file that the engine can hold.
  it("forgets the engine's start of an argv that is gone, and gives no file for an empty -conffile name", async () => {
    const { draft, registries } = setup({
      start: (argv) => {
        const at = argv.findIndex((word) => word === '-conffile' || word.startsWith('-conffile='));
        const name = at < 0 ? null : argv[at] === '-conffile' ? argv[at + 1]! : argv[at]!.slice('-conffile='.length);
        return { ok: true, start: { conffile: name, logFile: null, dryRun: false } };
      },
    });
    draft.setConffile(file('my.ini', '[BLAST]\n'));
    draft.setExtra('-conffile settings.ini');
    await draft.idle();
    expect(draft.state.get().start?.conffile).toBe('settings.ini');
    draft.setExtra("-conffile settings.ini -evalue '1");
    await draft.idle();
    expect(draft.state.get().validation.state).toBe('invalid');
    expect(draft.state.get().start).toBeUndefined();
    draft.setExtra('-conffile settings.ini');
    await draft.idle();
    expect(draft.state.get().start?.conffile).toBe('settings.ini');
    expect(conffileName(draft.state.get())).toBe('settings.ini');
    draft.setProgram('blastp');
    expect(draft.state.get().start).toBeUndefined();
    expect(conffileName(draft.state.get())).toBeUndefined();
    await draft.idle();
    expect(draft.state.get().start?.conffile).toBeNull();
    expect(conffileName(draft.state.get())).toBeUndefined();
    draft.setProgram('blastn');
    // `-conffile=`: the form asks for no file (conffileName, the row of ParameterForm), and the search reads none.
    draft.setExtra('-conffile=');
    await draft.idle();
    expect(draft.state.get().start?.conffile).toBe('');
    expect(conffileName(draft.state.get())).toBeUndefined();
    expect(registries.at(-1)).toBeUndefined();
  });

  it('shows the message of an argv that the engine cannot start, and validates nothing of it', async () => {
    const { draft, validated } = setup({
      start: (argv) => (argv.includes('-version') ? { ok: false, message: "error: unexpected argument '-version' found" } : { ok: true, start: { conffile: null, logFile: null, dryRun: false } }),
    });
    await draft.idle();
    const before = validated.length;
    draft.setExtra('-version');
    await draft.idle();
    expect(draft.state.get().validation).toEqual({ state: 'invalid', message: "error: unexpected argument '-version' found" });
    expect(draft.state.get().start).toBeUndefined();
    expect(validated).toHaveLength(before);
  });

  it('queues runs in the coordinator, a group with shared subject bytes', async () => {
    const { draft, coordinator } = setup({ useCoordinator: true });
    draft.addFiles('query', [file('q1.fa', '>q1\nACGT\n'), file('q2.fa', '>q2\nACGT\n')]);
    draft.setPaste('subject', '>s\nACGT\n');
    await draft.idle();
    draft.setMode('query', 'separate');
    const result = await draft.submit();
    expect(result.ok).toBe(true);
    const runs = coordinator.state.get().runs;
    expect(runs.map((run) => run.snapshot.group?.position)).toEqual([1, 2]);
    expect(runs[0]!.snapshot.group?.groupId).toBe(runs[1]!.snapshot.group?.groupId);
    expect(runs[0]!.snapshot.subject.bytes).toBe(runs[1]!.snapshot.subject.bytes);
    expect(runs.map((run) => run.snapshot.query.name)).toEqual(['q1.fa', 'q2.fa']);
  });

  it('reports a failure to queue instead of dropping it', async () => {
    const { draft } = setup({ useCoordinator: true, validate: () => { throw new Error('the Data worker stopped'); } });
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s\nACGT\n');
    const result = await draft.submit();
    expect(result).toEqual({ ok: false, message: 'The options could not be checked: the Data worker stopped' });
    expect(draft.state.get()).toMatchObject({ submitting: false, message: { kind: 'error' } });
  });

  it('refuses to combine a source that starts without a defline after another source', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', 'ACGTACGT\n');
    draft.addFiles('subject', [new File(['>s\nACGT\n'], 's.fa')]);
    await draft.idle();
    // The pasted subject comes first: its sequence without a defline is the input's first record.
    expect(draft.readiness()).toBeUndefined();
    draft.setPaste('subject', '');
    draft.addFiles('subject', [new File(['GGCC\n'], 'bare.fa')]);
    await draft.idle();
    expect(draft.readiness()).toBe(
      'Subject (file bare.fa) starts with a sequence without a defline (">"), which cannot follow another input.',
    );
    draft.setMode('subject', 'separate');
    expect(draft.readiness()).toBeUndefined();
  });

  it("maps a subject record that the engine refuses, and reports it in short above the button", async () => {
    const { draft } = setup();
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s1\nACGT\n>s2\nACGT!ACGT\n');
    await draft.idle();
    expect(sourceOf(draft, 'subject').check).toMatchObject({ state: 'refused', record: 1 });
    expect(draft.readiness()).toBe(
      'Subject (pasted): the engine refuses a record. Exclude it or correct the input to search.',
    );
    draft.setIncluded('subject', 'paste', [1], false);
    await draft.idle();
    expect(draft.readiness()).toBeUndefined();
  });

  it('does not let a check that could not run keep the search from the queue', async () => {
    const { draft, requests } = setup({
      checkInput: async () => {
        throw new Error('the engine could not allocate 9 bytes for an input');
      },
    });
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s\nACGT\n');
    await draft.idle();
    expect(sourceOf(draft, 'query').check).toEqual({
      state: 'error',
      message: 'the engine could not allocate 9 bytes for an input',
    });
    expect(draft.readiness()).toBeUndefined();
    expect((await draft.submit()).ok).toBe(true);
    expect(requests).toHaveLength(1);
  });

  it('checks a selection once: going back to an earlier selection or program uses its verdict', async () => {
    const { draft, checks } = setup();
    draft.setPaste('query', '>a\nACGT\n>b\nACGT\n');
    await draft.idle();
    draft.setIncluded('query', 'paste', [1], false);
    await draft.idle();
    draft.setIncluded('query', 'paste', [1], true);
    await draft.idle();
    draft.setProgram('tblastx');
    await draft.idle();
    draft.setProgram('blastn');
    await draft.idle();
    expect(checks.map((check) => check.split('|').slice(0, 2).join('|'))).toEqual([
      'blastn|query',
      'blastn|query',
      'tblastx|query',
    ]);
  });

  it('shows a pasted text as being read at once, before it is indexed', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>a\nACGT\n');
    await draft.idle();
    draft.setPaste('query', '>b\nACGT\n');
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'indexing' });
    await draft.idle();
    expect(sourceOf(draft, 'query').base?.records[0]?.id).toBe('b');
  });

  it("keeps a region with its record: another record leaves it aside", async () => {
    const { draft } = setup();
    draft.setPaste('subject', '>s1\nACGTACGTACGTACGTACGT\n');
    await draft.idle();
    draft.setRegion('subject', { start: '3', stop: '12' });
    expect(draft.parameters()).toEqual([['-subject_loc', '3-12']]);
    draft.setPaste('subject', '>s2\nACGTACGT\n');
    await draft.idle();
    expect(draft.region('subject')).toBeUndefined();
    expect(draft.parameters()).toEqual([]);
  });
});
