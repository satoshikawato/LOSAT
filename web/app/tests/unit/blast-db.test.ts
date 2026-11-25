import { describe, expect, it } from 'vitest';
import { sha256Hex } from 'losat';
import type { EngineRunRequest } from 'losat/ports';
import { Coordinator } from '../../src/application/coordinator';
import { holdsInputs, SearchDraft } from '../../src/application/draft';
import { ResultsBrowser } from '../../src/application/results';
import {
  assembleSearchDirectory,
  databaseTaxonomy,
  discoverBlastDatabaseNames,
  discoverBlastDatabases,
  duplicateRunFileMessage,
  FILTERING_DB_FLAG,
  matchingDatabases,
  noMatchingDatabaseMessage,
  normalizeRunPath,
  parseAliasMetadata,
  pickedFilePath,
  primaryBlastDbName,
} from '../../src/domain/blast-db';
import { hspTable } from '../../src/domain/hsp-table';
import { approvedExceptions } from '../../src/domain/verification';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import type { DatabaseMemory, DatabaseSlot, RememberedDatabase } from '../../src/ports/db-memory';
import { NodeInflater } from './support/node-inflater';

const TEST_MTIME = 1_764_100_000_000;

function diskFile(name: string, content: BlobPart, lastModified = TEST_MTIME): File {
  return new File([content], name, { lastModified });
}

function fileWithRelPath(relPath: string, content: BlobPart = 'data', lastModified = TEST_MTIME): File {
  const base = relPath.slice(relPath.lastIndexOf('/') + 1);
  const f = new File([content], base, { type: 'application/octet-stream', lastModified });
  Object.defineProperty(f, 'webkitRelativePath', { value: relPath, configurable: true });
  return f;
}

function memoryStore(initial: Partial<Record<DatabaseSlot, RememberedDatabase>> = {}): DatabaseMemory & {
  readonly data: Partial<Record<DatabaseSlot, RememberedDatabase>>;
} {
  const data: Partial<Record<DatabaseSlot, RememberedDatabase>> = { ...initial };
  return {
    data,
    read: (slot) => data[slot],
    write: (slot, value) => {
      if (value === undefined) delete data[slot];
      else data[slot] = value;
    },
  };
}

describe('blast-db domain discovery & directory assembly', () => {
  it('normalizes safe relative paths and rejects unsafe segments', () => {
    expect(pickedFilePath(fileWithRelPath('mydb/nt.nin'))).toBe('mydb/nt.nin');
    expect(pickedFilePath(fileWithRelPath('mydb\\sub\\nt.nin'))).toBe('mydb/sub/nt.nin');
    expect(pickedFilePath(new File(['x'], 'plain.nin'))).toBe('plain.nin');
    expect(normalizeRunPath('./sub/nt.nin')).toBe('sub/nt.nin');
    expect(normalizeRunPath('/leading/nt.nin')).toBe('leading/nt.nin');
    expect(pickedFilePath(fileWithRelPath('../escape.nin'))).toBeUndefined();
    expect(pickedFilePath(fileWithRelPath('a/../b.nin'))).toBeUndefined();
    expect(pickedFilePath(fileWithRelPath('a/./b.nin'))).toBeUndefined();
    expect(pickedFilePath(fileWithRelPath('bad\0name.nin'))).toBeUndefined();
  });

  it('uses case-sensitive NCBI FindBlastDBs extensions in sync discovery and primaryBlastDbName', () => {
    expect(discoverBlastDatabaseNames(['db/nt.nin', 'db/NT.NIN', 'db/prot.pin', 'db/PROT.PIN'])).toEqual([
      { name: 'db/nt', kind: 'nucleotide', entry: 'volume' },
      { name: 'db/prot', kind: 'protein', entry: 'volume' },
    ]);
    expect(primaryBlastDbName(['db/vol.00.nin', 'db/vol.01.nin', 'db/all.nal'])).toBe('db/all');
    expect(primaryBlastDbName(['db/ONLY.NIN'])).toBeUndefined();
  });

  it('parses TITLE and quoted/unquoted DBLIST from alias files', () => {
    const meta = parseAliasMetadata('# comment\nTITLE Test Multi Volume DB\nDBLIST "vol 00" vol_01\n');
    expect(meta).toEqual({ title: 'Test Multi Volume DB', dblist: ['vol 00', 'vol_01'] });
  });

  it('discovers single-volume, multi-volume aliases, molecule types, and folder-scoped taxonomy aux files without reading binary volumes', async () => {
    let nsqReads = 0;
    const nsq = fileWithRelPath('dbdir/multi.00.nsq', 'BINARY');
    nsq.text = () => {
      nsqReads++;
      return Promise.resolve('BINARY');
    };
    const files = [
      fileWithRelPath('dbdir/multi.nal', 'TITLE Multi Nuc\nDBLIST multi.00 multi.01\n'),
      fileWithRelPath('dbdir/multi.00.nin', 'idx0'),
      fileWithRelPath('dbdir/multi.00.nhr', 'hdr0'),
      nsq,
      fileWithRelPath('dbdir/multi.01.nin', 'idx1'),
      fileWithRelPath('dbdir/multi.01.nhr', 'hdr1'),
      fileWithRelPath('dbdir/multi.01.nsq', 'seq1'),
      fileWithRelPath('dbdir/prot.pin', 'pidx'),
      fileWithRelPath('dbdir/prot.phr', 'phdr'),
      fileWithRelPath('dbdir/prot.psq', 'pseq'),
      fileWithRelPath('dbdir/taxdb.bti', 'bti'),
      fileWithRelPath('dbdir/taxdb.btd', 'btd'),
      fileWithRelPath('dbdir/taxonomy4blast.sqlite3', 'sql'),
    ];

    const discovery = await discoverBlastDatabases(files);
    expect(nsqReads).toBe(0);
    expect(discovery.folder).toBe('dbdir');
    expect(databaseTaxonomy(discovery, 'dbdir/multi')).toMatchObject({ taxdb: true, sqlite: true });
    expect(discovery.databases).toEqual([
      { name: 'dbdir/multi', kind: 'nucleotide', entry: 'alias', title: 'Multi Nuc', volumes: ['multi.00', 'multi.01'] },
      { name: 'dbdir/prot', kind: 'protein', entry: 'volume' },
    ]);
    expect(matchingDatabases(discovery, 'nucleotide').map((d) => d.name)).toEqual(['dbdir/multi']);
    expect(matchingDatabases(discovery, 'protein').map((d) => d.name)).toEqual(['dbdir/prot']);
  });

  it('reports clear molecule-type mismatch or missing database messages', async () => {
    const protOnly = await discoverBlastDatabases([fileWithRelPath('pdir/p.pin', 'x')]);
    expect(noMatchingDatabaseMessage(protOnly, 'nucleotide', 'BLASTN')).toBe(
      'Found protein database (pdir/p), not a nucleotide database for BLASTN.',
    );

    const empty = await discoverBlastDatabases([fileWithRelPath('notes/readme.txt', 'hello')]);
    expect(noMatchingDatabaseMessage(empty, 'protein', 'BLASTP')).toBe(
      'No protein BLAST database (.pal or .pin) was found in the chosen files.',
    );
  });

  it('mounts DB files, root aliases for selected-folder taxonomy files with blastdbEnv, filter lists by name, and deduplicates identical Files', async () => {
    const dbFiles = [
      fileWithRelPath('mydb/nt.nin', 'nin'),
      fileWithRelPath('mydb/nt.nhr', 'nhr'),
      fileWithRelPath('mydb/nt.nsq', 'nsq'),
      fileWithRelPath('mydb/taxdb.bti', 'bti'),
      fileWithRelPath('mydb/taxdb.btd', 'btd'),
      fileWithRelPath('mydb/taxonomy4blast.sqlite3', 'sqlite'),
    ];
    const discovery = await discoverBlastDatabases(dbFiles);
    const seqidList = new File(['seq1\n'], 'keep.ids', { type: 'text/plain' });

    const assembled = assembleSearchDirectory({
      subjectDb: { discovery, selectedName: 'mydb/nt' },
      filteringDb: { discovery, selectedName: 'mydb/nt' },
      filterFiles: { '-seqidlist': seqidList, '-negative_seqidlist': seqidList },
    });
    expect(assembled.blastdbEnv).toBe('mydb');
    const names = assembled.files.map((entry) => entry.name).sort();
    expect(names).toEqual([
      'keep.ids',
      'mydb/nt.nhr',
      'mydb/nt.nin',
      'mydb/nt.nsq',
      'mydb/taxdb.btd',
      'mydb/taxdb.bti',
      'mydb/taxonomy4blast.sqlite3',
      'taxdb.btd',
      'taxdb.bti',
      'taxonomy4blast.sqlite3',
    ]);
  });

  it('aliases taxonomy files ONLY from the selected DB folder and omits blastdbEnv for root picks or folders without taxonomy files (F1)', async () => {
    const nested = await discoverBlastDatabases([
      fileWithRelPath('root/subA/a.pin', 'a'),
      fileWithRelPath('root/subB/b.pin', 'b'),
      fileWithRelPath('root/subB/taxdb.bti', 'bti'),
      fileWithRelPath('root/subB/taxdb.btd', 'btd'),
      fileWithRelPath('root/subB/taxonomy4blast.sqlite3', 'sql'),
    ]);
    // Selecting subA/a must NOT pick up subB's taxonomy files.
    const forA = assembleSearchDirectory({ subjectDb: { discovery: nested, selectedName: 'root/subA/a' } });
    expect(forA.blastdbEnv).toBeUndefined();
    expect(forA.files.some((f) => f.name === 'taxonomy4blast.sqlite3' || f.name === 'taxdb.bti')).toBe(false);
    expect(databaseTaxonomy(nested, 'root/subA/a')).toMatchObject({ taxdb: false, sqlite: false });

    // Selecting subB/b aliases subB's taxonomy files to root and sets blastdbEnv = 'root/subB'.
    const forB = assembleSearchDirectory({ subjectDb: { discovery: nested, selectedName: 'root/subB/b' } });
    expect(forB.blastdbEnv).toBe('root/subB');
    expect(forB.files.map((f) => f.name)).toEqual(
      expect.arrayContaining(['taxdb.bti', 'taxdb.btd', 'taxonomy4blast.sqlite3']),
    );

    // Multi-file root pick already has taxonomy files at root, so blastdbEnv is omitted.
    const rootPick = await discoverBlastDatabases([
      new File(['p'], 'tfp.pin'),
      new File(['s'], 'taxonomy4blast.sqlite3'),
    ]);
    const forRoot = assembleSearchDirectory({ subjectDb: { discovery: rootPick, selectedName: 'tfp' } });
    expect(forRoot.blastdbEnv).toBeUndefined();
  });

  it('refuses when two different files claim the same mount path (F2 / G1)', async () => {
    const db = await discoverBlastDatabases([fileWithRelPath('mydb/nt.nin', 'nin')]);
    const subjectDb = { discovery: db, selectedName: 'mydb/nt' };
    const listA = diskFile('ids.txt', 'seq1\n');
    const listDiffMtime = diskFile('ids.txt', 'seq2\n', TEST_MTIME + 1);
    const listDiffSize = diskFile('ids.txt', 'seq11\n');
    const expectDup = (fn: () => unknown, name: string) => expect(fn).toThrow(duplicateRunFileMessage(name));

    // Same name and size, different lastModified -> refused.
    expectDup(
      () => assembleSearchDirectory({ subjectDb, filterFiles: { '-seqidlist': listA, '-negative_seqidlist': listDiffMtime } }),
      'ids.txt',
    );

    // Same name and lastModified, different size -> refused.
    expectDup(
      () => assembleSearchDirectory({ subjectDb, filterFiles: { '-seqidlist': listA, '-negative_seqidlist': listDiffSize } }),
      'ids.txt',
    );

    // Filter file colliding with a multi-file root DB file of different size.
    const flatDb = await discoverBlastDatabases([diskFile('nt.nin', 'nin'), diskFile('ids.txt', 'old')]);
    expectDup(
      () => assembleSearchDirectory({ subjectDb: { discovery: flatDb, selectedName: 'nt' }, filterFiles: { '-seqidlist': listA } }),
      'ids.txt',
    );

    // -db and -filtering_db having the same path with different sizes.
    const otherDb = await discoverBlastDatabases([fileWithRelPath('mydb/nt.nin', 'different')]);
    expectDup(
      () => assembleSearchDirectory({ subjectDb, filteringDb: { discovery: otherDb, selectedName: 'mydb/nt' } }),
      'mydb/nt.nin',
    );

    // Filter list file named like the DB folder (file vs. directory collision).
    expectDup(
      () => assembleSearchDirectory({ subjectDb, filterFiles: { '-seqidlist': diskFile('mydb', 's1\n') } }),
      'mydb',
    );
  });
});

describe('SearchDraft, Coordinator & ResultsBrowser with -db and -filtering_db', () => {
  function setup(memory?: DatabaseMemory) {
    const engine = new FakeEngine();
    const runCalls: EngineRunRequest[] = [];
    const origRun = engine.run.bind(engine);
    engine.run = (req, out, onPhase) => {
      runCalls.push(req);
      return origRun(req, out, onPhase);
    };
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
    const coordinator = new Coordinator({
      engine,
      data,
      downloader: { save: () => undefined },
      now: () => 1_700_000_000_000,
      newRunId: (() => {
        let n = 1;
        return () => `run-${n++}`;
      })(),
      digest: sha256Hex,
    });
    const draft = new SearchDraft({
      engine,
      data,
      enqueueAll: (requests) => coordinator.enqueueAll(requests),
      fetcher: { fetch: async () => ({ ok: false, reason: 'not-found' }) },
      ...(memory === undefined ? {} : { databaseMemory: memory }),
      now: () => 0,
      debounceMs: 0,
    });
    const results = new ResultsBrowser({
      data,
      describe: (program) => engine.describe(program),
      runs: coordinator.state,
      verification: { ncbi: '2.17.0', sources: [], programs: {} },
    });
    return { draft, coordinator, results, runCalls };
  }

  async function waitCompleted(coordinator: Coordinator, runId: string): Promise<void> {
    await new Promise<void>((resolve) => {
      const check = () => {
        const run = coordinator.state.get().runs.find((r) => r.snapshot.runId === runId);
        if (run?.status === 'completed' || run?.status === 'failed') {
          unsub();
          resolve();
        }
      };
      const unsub = coordinator.state.subscribe(check);
      check();
    });
  }

  const pasteQuery = (draft: SearchDraft) => draft.setPaste('query', '>q1\nACGTACGTACGTACGTACGTACGT\n');
  const pasteSubject = (draft: SearchDraft) => draft.setPaste('subject', '>s1\nACGTACGTACGTACGTACGTACGT\n');
  const expectReady = (draft: SearchDraft) => expect(draft.readiness()).toBeUndefined();
  const pickNucDb = (draft: SearchDraft) => draft.setSubjectDatabaseFiles([fileWithRelPath('nucdb/nt.nin', 'nin')]);
  const submitAndWait = async (draft: SearchDraft, coordinator: Coordinator): Promise<string> => {
    const submitted = await draft.submit();
    expect(submitted.ok).toBe(true);
    const runId = submitted.runIds![0]!;
    await waitCompleted(coordinator, runId);
    return runId;
  };
  const mountedNames = (req: EngineRunRequest) => req.files?.map((f) => f.name).sort();
  const expectBlastnTail = (req: EngineRunRequest, tail: readonly string[]) => {
    expect(req.argv[0]).toBe('blastn');
    expect(req.argv.slice(3)).toEqual(tail);
  };
  const submitSingle = async (ctx: ReturnType<typeof setup>): Promise<EngineRunRequest> => {
    await submitAndWait(ctx.draft, ctx.coordinator);
    expect(ctx.runCalls).toHaveLength(1);
    return ctx.runCalls[0]!;
  };

  it('runs a local -db search with filter list files, mounts directory read-only, records blastdbEnv, omits FASTA subject, and reads subjectLength in ResultsBrowser', async () => {
    const memory = memoryStore();
    const { draft, coordinator, results, runCalls } = setup(memory);
    pasteQuery(draft);
    draft.setSubjectSource('database');
    await draft.setSubjectDatabaseFiles([
      fileWithRelPath('nucdb/nt.nin', 'nin'),
      fileWithRelPath('nucdb/nt.nhr', 'nhr'),
      fileWithRelPath('nucdb/nt.nsq', 'nsq'),
      fileWithRelPath('nucdb/taxdb.bti', 'bti'),
      fileWithRelPath('nucdb/taxdb.btd', 'btd'),
    ]);
    draft.setDatabaseFilter('-taxids', '9606');
    draft.setDatabaseFilter('-no_taxid_expansion', true);
    draft.setDatabaseFilterFile('-seqidlist', new File(['s1\n'], 'subset.ids'));
    await draft.idle();

    expect(holdsInputs(draft.state.get())).toBe(true);
    expectReady(draft);
    expect(memory.data.subject).toEqual({ name: 'nucdb/nt', folder: 'nucdb' });

    const runId = await submitAndWait(draft, coordinator);

    const runView = coordinator.state.get().runs.find((r) => r.snapshot.runId === runId)!;
    expect(runView.snapshot.kind !== 'tool' && runView.snapshot.database?.blastdbEnv).toBe('nucdb');

    expect(runCalls).toHaveLength(1);
    const req = runCalls[0]!;
    expect(req.subject).toBeUndefined();
    expect(req.argv[0]).toBe('blastn');
    expect(req.argv[3]).toBe('-db');
    expect(req.argv[4]).toBe('nucdb/nt');
    expect(req.argv.slice(5)).toEqual(['-taxids', '9606', '-no_taxid_expansion', '-seqidlist', 'subset.ids']);
    expect(req.files?.map((f) => f.name).sort()).toEqual([
      'nucdb/nt.nhr',
      'nucdb/nt.nin',
      'nucdb/nt.nsq',
      'nucdb/taxdb.btd',
      'nucdb/taxdb.bti',
      'subset.ids',
      'taxdb.btd',
      'taxdb.bti',
    ]);

    // ResultsBrowser reads subject lengths from HspTable even though snapshot.subject has 0 FASTA records.
    await results.open(runId);
    const rs = results.state.get();
    expect(rs.phase).toBe('ready');
    expect(rs.subjects.length).toBeGreaterThan(0);
    expect(rs.subjects[0]!.length).toBe(500);
    expect(rs.subjects[0]!.recordId).toBe('db_seq_1');
  });

  it('refuses submit when two filter list flags pick different files with the same name (F2 / G1)', async () => {
    const { draft } = setup();
    pasteQuery(draft);
    draft.setSubjectSource('database');
    await pickNucDb(draft);
    draft.setDatabaseFilterFile('-seqidlist', diskFile('ids.txt', 's1\n'));
    draft.setDatabaseFilterFile('-negative_seqidlist', diskFile('ids.txt', 's2\n', TEST_MTIME + 1));
    await draft.idle();

    expect(draft.readiness()).toBe(duplicateRunFileMessage('ids.txt'));
    const submitted = await draft.submit();
    expect(submitted).toEqual({ ok: false, message: duplicateRunFileMessage('ids.txt') });
  });

  it('allows picking the same folder twice for -db and -filtering_db and mounts each file once (G1)', async () => {
    const ctx = setup();
    const pickSharedFolder = () =>
      ['nt.nin', 'nt.nhr', 'nt.nsq', 'repeats.nin', 'repeats.nhr', 'repeats.nsq'].map((name) =>
        fileWithRelPath(`blastdb/${name}`, name, TEST_MTIME),
      );
    pasteQuery(ctx.draft);
    ctx.draft.setSubjectSource('database');
    await ctx.draft.setSubjectDatabaseFiles(pickSharedFolder());
    ctx.draft.selectSubjectDatabase('blastdb/nt');
    ctx.draft.setFilteringDbEnabled(true);
    await ctx.draft.setFilteringDbFiles(pickSharedFolder());
    ctx.draft.selectFilteringDb('blastdb/repeats');
    await ctx.draft.idle();

    expectReady(ctx.draft);
    const req = await submitSingle(ctx);
    expectBlastnTail(req, ['-db', 'blastdb/nt', FILTERING_DB_FLAG, 'blastdb/repeats']);
    expect(mountedNames(req)).toEqual([
      'blastdb/nt.nhr',
      'blastdb/nt.nin',
      'blastdb/nt.nsq',
      'blastdb/repeats.nhr',
      'blastdb/repeats.nin',
      'blastdb/repeats.nsq',
    ]);
  });

  it('allows picking the same disk file in a multi-file DB and a -seqidlist filter (G1)', async () => {
    const runCtx = setup();
    pasteQuery(runCtx.draft);
    runCtx.draft.setSubjectSource('database');
    await runCtx.draft.setSubjectDatabaseFiles([diskFile('nt.nin', 'nin'), diskFile('ids.txt', 's1\n')]);
    runCtx.draft.setDatabaseFilterFile('-seqidlist', diskFile('ids.txt', 's1\n'));
    await runCtx.draft.idle();

    expectReady(runCtx.draft);
    const req = await submitSingle(runCtx);
    expectBlastnTail(req, ['-db', 'nt', '-seqidlist', 'ids.txt']);
    expect(mountedNames(req)).toEqual(['ids.txt', 'nt.nin']);
  });

  it('writes -filtering_db explicitly to BLASTN when the user switches program during folder scan (F4)', async () => {
    const env = setup();
    const { draft, coordinator, runCalls } = env;
    pasteQuery(draft);
    pasteSubject(draft);

    let resolveAlias!: (text: string) => void;
    const delayedAlias = fileWithRelPath('rep/repeats.nal', 'TITLE Delayed\nDBLIST repeats.00\n');
    delayedAlias.text = () =>
      new Promise<string>((resolve) => {
        resolveAlias = resolve;
      });
    const vol = fileWithRelPath('rep/repeats.00.nin', 'nin');

    const scanPromise = draft.setFilteringDbFiles([delayedAlias, vol]);
    // Switch to BLASTP while the repeats DB alias is still being read.
    draft.setProgram('blastp');
    resolveAlias('TITLE Delayed\nDBLIST repeats.00\n');
    await scanPromise;

    expect(draft.state.get().values.blastp[FILTERING_DB_FLAG]).toBeUndefined();
    expect(draft.state.get().values.blastn[FILTERING_DB_FLAG]).toBe('rep/repeats');

    // Switching back to BLASTN must allow Run immediately.
    draft.setProgram('blastn');
    await draft.idle();
    expectReady(draft);
    const f4RunId = await submitAndWait(draft, coordinator);
    expect(f4RunId).toMatch(/^run-/);
    expect(runCalls[0]!.argv.slice(-2)).toEqual([FILTERING_DB_FLAG, 'rep/repeats']);
  });

  it('clears remembered subject DB when switching back to FASTA or loading an example so reload opens in FASTA mode (F5)', async () => {
    const memory = memoryStore();
    const { draft } = setup(memory);
    const expectReloadedFasta = async () => {
      expect(memory.data.subject).toBeUndefined();
      const fresh = setup(memory).draft;
      await fresh.idle();
      expect(fresh.state.get().subjectSource).toBe('fasta');
      expect(fresh.state.get().database.remembered).toBeUndefined();
    };
    const nucFile = fileWithRelPath('nucdb/nt.nin', 'nin');
    const remembered = { name: 'nucdb/nt', folder: 'nucdb' };

    draft.setSubjectSource('database');
    await draft.setSubjectDatabaseFiles([nucFile]);
    await draft.idle();
    expect(memory.data.subject).toEqual(remembered);

    // Switch back to FASTA -> memory clears -> a fresh draft (simulating reload) opens in FASTA mode.
    draft.setSubjectSource('fasta');
    await expectReloadedFasta();

    // Also verify useExample clears remembered subject DB.
    draft.setSubjectSource('database');
    expect(memory.data.subject).toStrictEqual(remembered);
    draft.useExample('ecoli-k12-genome', '>q\nACGT\n', new File(['>s\nACGT\n'], 'subject.fa'));
    await expectReloadedFasta();
  });

  it('re-validates molecule type on program switch and restores reload prompt from DatabaseMemory', async () => {
    const memory = memoryStore({ subject: { name: 'nucdb/nt', folder: 'nucdb' } });
    const { draft } = setup(memory);
    await draft.idle();
    expect(draft.state.get().subjectSource).toBe('database');
    expect(draft.state.get().database.remembered).toEqual({ name: 'nucdb/nt', folder: 'nucdb' });
    pasteQuery(draft);
    await draft.idle();
    expect(draft.readiness()).toBe('Choose the folder again.');

    await pickNucDb(draft);
    await draft.idle();
    expect(draft.state.get().database.remembered).toBeUndefined();
    expect(draft.state.get().database.selected).toBe('nucdb/nt');

    // Switching to BLASTP requires a protein database: nucleotide selection clears with clear mismatch.
    draft.setProgram('blastp');
    await draft.idle();
    expect(draft.state.get().database.selected).toBeUndefined();
    expect(draft.readiness()).toBe('Found nucleotide database (nucdb/nt), not a protein database for BLASTP.');
  });

  it('supports BLASTN -filtering_db with both -subject and -db and clears on resetFields', async () => {
    const { draft, coordinator, runCalls } = setup();
    pasteQuery(draft);
    pasteSubject(draft);
    draft.setFilteringDbEnabled(true);
    await draft.idle();
    expect(draft.readiness()).toBe('Choose a repeats database folder or files.');

    await draft.setFilteringDbFiles([
      fileWithRelPath('rep/repeats.nin', 'nin'),
      fileWithRelPath('rep/repeats.nhr', 'nhr'),
      fileWithRelPath('rep/repeats.nsq', 'nsq'),
    ]);
    await draft.idle();
    expectReady(draft);

    await submitAndWait(draft, coordinator);
    expect(runCalls[0]!.argv).toContain(FILTERING_DB_FLAG);
    expect(runCalls[0]!.argv).toContain('rep/repeats');
    expect(runCalls[0]!.files?.map((f) => f.name).sort()).toEqual(['rep/repeats.nhr', 'rep/repeats.nin', 'rep/repeats.nsq']);

    draft.resetFields([FILTERING_DB_FLAG]);
    expect(draft.state.get().filteringDb.discovery).toBeUndefined();
    expect(draft.state.get().values.blastn[FILTERING_DB_FLAG]).toBeUndefined();
  });

  it('stores subject_length in HspTable and omits local -subject -db_gencode exception for -db runs', () => {
    const table = hspTable([
      {
        index: 0,
        q_idx: 0,
        s_idx: 0,
        rank: 0,
        raw_score: 50,
        bit_score: 25,
        e_value: 1e-3,
        q_start: 1,
        q_end: 20,
        s_start: 1,
        s_end: 20,
        query_frame: null,
        subject_frame: null,
        subject_length: 12345,
        out6: [0, 10],
        out0: null,
        out0_subject: null,
      },
    ]);
    expect(table.subjectLength[0]).toBe(12345);
    expect(approvedExceptions('tblastx', ['tblastx', '-query', 'q.fa', '-subject', 's.fa', '-db_gencode', '4'])).toHaveLength(1);
    expect(approvedExceptions('tblastx', ['tblastx', '-query', 'q.fa', '-db', 'mydb/nt', '-db_gencode', '4'])).toEqual([]);
  });
});
