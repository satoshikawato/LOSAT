// The makeblastdb tab (line WW, WW4b2): its options in NCBI's groups, the argv beside Run LOSAT, the run
// directory, the check, the hand-over of mask data, and the helpers of its results (names, the tar of all files).
import { describe, expect, it, vi } from 'vitest';
import { sha256Hex } from 'losat';
import { parseToolDescription } from 'losat/domain';
import type { ToolEngine } from 'losat/ports';
import { Coordinator, madeDatabaseFolder, type EnqueueAllResult, type ToolRequest } from '../../src/application/coordinator';
import { MakeblastdbDraft } from '../../src/application/makeblastdb-draft';
import { assembleSearchDirectory, discoverBlastDatabases, primaryBlastDb } from '../../src/domain/blast-db';
import { tarArchive } from '../../src/domain/tar';
import {
  buildMakeblastdbArgv,
  changedMakeblastdbOptionsCount,
  fileListNames,
  makeblastdbInputList,
  maskDataName,
  refusedToolFields,
} from '../../src/domain/tool-argv';
import {
  makeblastdbOptionSections,
  toolLabel,
  toolOutputDownloadName,
  toolPrintsKind,
  toolRunTitle,
  toolWritesDirectory,
} from '../../src/domain/tools';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import MAKEBLASTDB_XMLHELP_DOC from '../../src/infra/fake/makeblastdb-xmlhelp.json';
import { NodeInflater } from './support/node-inflater';

const MK_DESC = parseToolDescription(MAKEBLASTDB_XMLHELP_DOC.xml);

describe('makeblastdb domain helpers', () => {
  it('groups the options as -xmlhelp groups them, in the order of NCBI usage, without the Input and Output blocks', () => {
    const sections = makeblastdbOptionSections(MK_DESC).map((section) => [section.title, section.parameters.map((p) => p.flag)]);
    expect(sections).toEqual([
      ['Configuration options', ['-title', '-parse_seqids', '-hash_index']],
      ['Sequence masking options', ['-mask_data', '-mask_id', '-mask_desc', '-gi_mask', '-gi_mask_name']],
      ['Output options', ['-blastdb_version', '-max_file_sz', '-metadata_output_prefix', '-logfile']],
      ['Taxonomy options', ['-taxid', '-taxid_map', '-oid_masks']],
    ]);
  });

  it('writes -in as the command line takes a list (spaces separate, quotes keep a name with spaces)', () => {
    expect(
      buildMakeblastdbArgv({
        input: makeblastdbInputList(['input.fa', 'my genome.fa.gz']),
        values: { '-dbtype': 'nucl', '-parse_seqids': true, '-out': 'db', '-mask_data': 'a.asnb, b.asnb' },
        description: MK_DESC,
        extra: ['-logfile', 'run.log'],
      }),
    ).toEqual([
      'makeblastdb',
      '-in',
      'input.fa "my genome.fa.gz"',
      '-dbtype',
      'nucl',
      '-parse_seqids',
      '-mask_data',
      'a.asnb, b.asnb',
      '-out',
      'db',
      '-logfile',
      'run.log',
    ]);
    expect(buildMakeblastdbArgv({ input: '', values: { '-dbtype': 'prot' }, description: MK_DESC })).toEqual([
      'makeblastdb',
      '-dbtype',
      'prot',
    ]);
    expect(fileListNames(' a.asnb, ,b.asnb ')).toEqual(['a.asnb', 'b.asnb']);
    expect(fileListNames(undefined)).toEqual([]);
  });

  it('counts the options changed in "makeblastdb options" only', () => {
    expect(changedMakeblastdbOptionsCount({ '-dbtype': 'nucl', '-out': 'db', '-input_type': 'fasta' }, '', MK_DESC)).toBe(0);
    expect(changedMakeblastdbOptionsCount({ '-parse_seqids': true, '-title': 'x' }, '-logfile a.log', MK_DESC)).toBe(3);
  });

  it('names its runs and their files, and prints help for -help-full as makeblastdb does', () => {
    expect(toolLabel('makeblastdb')).toBe('makeblastdb');
    expect(toolRunTitle('makeblastdb', 'build')).toBe('makeblastdb');
    expect(toolRunTitle('blast_formatter', 'format')).toBe('Format Archive (blast_formatter)');
    expect(toolRunTitle('windowmasker', 'mask')).toMatch(/^WindowMasker · /);
    expect(toolWritesDirectory('makeblastdb')).toBe(true);
    expect(toolWritesDirectory('windowmasker')).toBe(false);
    const run = { number: 3, tool: 'makeblastdb', stage: 'build', queuedAt: Date.UTC(2026, 10, 27, 12) } as const;
    expect(toolOutputDownloadName(run, { key: 'f1', name: 'db/mydb.nin' })).toBe('mydb.nin');
    expect(toolOutputDownloadName(run, { key: 'stdout', name: 'stdout' })).toMatch(/run.*3.*stdout\.txt$/);
    expect(toolPrintsKind(['makeblastdb', '-help-full'])).toBe('help');
    expect(toolPrintsKind(['windowmasker', '-help-full'])).toBeNull();
  });

  it("offers a masking run's maskinfo output as -mask_data: the -out file, or the standard output under a name", () => {
    const mask = (argv: readonly string[]) => ({ tool: 'windowmasker', stage: 'mask', argv, number: 7 }) as const;
    const withOut = mask(['windowmasker', '-ustat', 'c.obinary', '-in', 'a.fa', '-outfmt', 'maskinfo_asn1_bin', '-out', 'wm.asnb']);
    expect(maskDataName(withOut, { key: 'f0', name: 'wm.asnb' })).toBe('wm.asnb');
    expect(maskDataName(withOut, { key: 'stdout', name: 'stdout' })).toBeUndefined();
    const toStdout = mask(['windowmasker', '-ustat', 'c.obinary', '-in', 'a.fa', '-outfmt', 'maskinfo_xml']);
    expect(maskDataName(toStdout, { key: 'stdout', name: 'stdout' })).toBe('windowmasker-run7.xml');
    expect(maskDataName(mask(['windowmasker', '-ustat', 'c', '-in', 'a.fa', '-outfmt', 'fasta']), { key: 'stdout', name: 'stdout' })).toBeUndefined();
    expect(
      maskDataName({ tool: 'windowmasker', stage: 'mk_counts', argv: ['windowmasker', '-mk_counts', '-in', 'a.fa'], number: 1 }, { key: 'stdout', name: 'x' }),
    ).toBeUndefined();
  });

  it('names the field of a missing -dbtype', () => {
    const fields = refusedToolFields('Command line argument error: Argument "dbtype". Mandatory value is missing:  `String, `nucl\', `prot\'', '');
    expect([...fields.flags]).toEqual(['-dbtype']);
  });

  it('finds a database by the relative names of a run\'s files', async () => {
    const files = ['nin', 'nhr', 'nsq'].map((ext) => new File(['x'], `mydb.${ext}`));
    const found = await discoverBlastDatabases(files, files.map((file) => `db/${file.name}`));
    expect(found.databases.map((db) => [db.name, db.kind])).toEqual([['db/mydb', 'nucleotide']]);
    expect(primaryBlastDb(['db/mydb.pin', 'db/mydb.phr', 'db/mydb.psq'])).toMatchObject({ name: 'db/mydb', kind: 'protein' });
    expect(primaryBlastDb(['report.txt'])).toBeUndefined();
  });
});

describe('tarArchive', () => {
  const read = async (blob: Blob) => new Uint8Array(await blob.arrayBuffer());
  const text = (bytes: Uint8Array, from: number, length: number) => new TextDecoder().decode(bytes.subarray(from, from + length)).replace(/\0+$/, '');

  it('writes ustar headers with their checksum, the files padded to 512 bytes, and two zero blocks', async () => {
    const bytes = await read(tarArchive([{ name: 'db.nin', file: new Blob(['abc']) }], 1_700_000_000_000));
    expect(bytes.length).toBe(512 + 512 + 1024);
    expect(text(bytes, 0, 100)).toBe('db.nin');
    expect(text(bytes, 100, 7)).toBe('0000644');
    expect(text(bytes, 124, 11)).toBe('00000000003');
    expect(text(bytes, 136, 11)).toBe((1_700_000_000).toString(8).padStart(11, '0'));
    expect(text(bytes, 257, 6)).toBe('ustar');
    expect(String.fromCharCode(bytes[156]!)).toBe('0');
    const header = bytes.slice(0, 512);
    header.fill(0x20, 148, 156);
    const sum = header.reduce((total, byte) => total + byte, 0);
    expect(text(bytes, 148, 6)).toBe(sum.toString(8).padStart(6, '0'));
    expect(text(bytes, 512, 3)).toBe('abc');
    expect(bytes.subarray(1024).every((byte) => byte === 0)).toBe(true);
  });

  it('splits a long path into prefix and name, and gives a pax header to a name that does not fit', async () => {
    const dir = 'd'.repeat(120);
    const split = await read(tarArchive([{ name: `${dir}/db.nin`, file: new Blob([]) }], 0));
    expect(text(split, 0, 100)).toBe('db.nin');
    expect(text(split, 345, 155)).toBe(dir);

    const long = `${'n'.repeat(130)}.nin`;
    const pax = await read(tarArchive([{ name: long, file: new Blob(['z']) }], 0));
    expect(String.fromCharCode(pax[156]!)).toBe('x');
    const record = text(pax, 512, 512).replace(/\0+$/, '');
    expect(record).toBe(`${record.length} path=${long}\n`);
    expect(String.fromCharCode(pax[1024 + 156]!)).toBe('0');
    expect(text(pax, 1536, 1)).toBe('z');
  });
});

describe('MakeblastdbDraft', () => {
  function setup() {
    const fake = new FakeEngine();
    const validateTool = vi.fn((...args: Parameters<ToolEngine['validateTool']>) => fake.validateTool(args[0]));
    const engine = { describeTool: (tool: string) => fake.describeTool(tool), validateTool };
    const queued: ToolRequest[] = [];
    const enqueueToolAll = vi.fn(async (requests: readonly ToolRequest[]): Promise<EnqueueAllResult> => {
      queued.push(...requests);
      return { ok: true, runIds: requests.map((_, idx) => `run-${idx + 1}`) };
    });
    const draft = new MakeblastdbDraft({ engine, enqueueToolAll, isBusy: () => false, debounceMs: 0 });
    return { draft, queued, validateTool };
  }

  async function loaded(draft: MakeblastdbDraft): Promise<void> {
    draft.ensureLoaded();
    await vi.waitFor(() => expect(draft.state.get().description).toBeDefined(), { timeout: 10_000 });
  }

  it('runs pasted FASTA and opened files as one -in list, with the files of the run', async () => {
    const { draft, queued } = setup();
    await loaded(draft);
    expect(draft.holdsInputs()).toBe(false);
    draft.setFastaPaste('>s1\nACGT\n');
    draft.addInputFiles([new File(['>s2\nGGCC\n'], 'b.fna.gz')]);
    draft.setField('-dbtype', 'nucl');
    draft.setField('-out', 'db/mydb');
    draft.setField('-parse_seqids', true);
    draft.addRunFiles([new File(['1\t9606\n'], 'taxmap.txt')]);
    draft.setField('-taxid_map', 'taxmap.txt');
    expect(draft.holdsInputs()).toBe(true);
    expect(draft.previewCommandWords()).toEqual([
      'losat',
      'makeblastdb',
      '-in',
      'input.fa b.fna.gz',
      '-dbtype',
      'nucl',
      '-parse_seqids',
      '-out',
      'db/mydb',
      '-taxid_map',
      'taxmap.txt',
    ]);
    const res = await draft.submit();
    expect(res.ok).toBe(true);
    expect(queued[0]).toMatchObject({ tool: 'makeblastdb', stage: 'build' });
    expect(queued[0]!.directory.map((entry) => [entry.name, entry.origin?.kind])).toEqual([
      ['input.fa', 'pasted'],
      ['b.fna.gz', 'picked'],
      ['taxmap.txt', 'picked'],
    ]);
    expect(await queued[0]!.directory[0]!.file.text()).toBe('>s1\nACGT\n');
  });

  it('refuses an input file under the name of the pasted FASTA, and says what is missing before a run', async () => {
    const { draft } = setup();
    await loaded(draft);
    draft.setField('-dbtype', 'prot');
    const none = await draft.submit();
    expect(none.ok).toBe(false);
    expect(!none.ok && none.message).toBe('Enter sequences or open FASTA files.');
    draft.setFastaPaste('>p\nMKV\n');
    draft.addInputFiles([new File(['x'], 'input.fa')]);
    expect(draft.state.get().inputFiles).toHaveLength(0);
    expect(draft.state.get().message?.text).toBe('A file named "input.fa" is already in the run directory.');
  });

  it('checks nothing before an input is given, then reports a missing -dbtype in its field', async () => {
    const { draft, validateTool } = setup();
    await loaded(draft);
    draft.setField('-title', 'x');
    await new Promise((resolve) => setTimeout(resolve, 10));
    expect(draft.state.get().validation.state).toBe('idle');
    expect(validateTool).not.toHaveBeenCalled();
    draft.setFastaPaste('>s\nACGT\n');
    await vi.waitFor(() => expect(draft.state.get().validation.state).toBe('invalid'));
    const validation = draft.state.get().validation;
    expect(validation.state === 'invalid' && [...refusedToolFields(validation.message, '').flags]).toEqual(['-dbtype']);
    draft.setField('-dbtype', 'nucl');
    await vi.waitFor(() => expect(draft.state.get().validation.state).toBe('ok'));
  });

  it('names a BLAST database of the files of the run as -in for -input_type blastdb, until the name is typed', async () => {
    const { draft, queued } = setup();
    await loaded(draft);
    draft.setField('-input_type', 'blastdb');
    draft.setField('-dbtype', 'nucl');
    const files = ['nin', 'nhr', 'nsq'].map((ext) => new File(['x'], `old.${ext}`));
    for (const file of files) Object.defineProperty(file, 'webkitRelativePath', { value: `in/${file.name}` });
    draft.addRunFiles(files);
    expect(draft.state.get().dbIn).toBe('in/old');
    draft.setDbIn('in/other');
    draft.addRunFiles([new File(['x'], 'more.nal')]);
    expect(draft.state.get().dbIn).toBe('in/other');
    draft.setDbIn('in/old');
    draft.setField('-out', 'new');
    const res = await draft.submit();
    expect(res.ok).toBe(true);
    expect(queued[0]!.argv).toEqual(['makeblastdb', '-in', 'in/old', '-input_type', 'blastdb', '-dbtype', 'nucl', '-out', 'new']);
  });

  it('writes -in of -input_type blastdb as typed (a list of databases), and a found name with a space quoted, in the field too (WW4b2 R1 finding 1)', async () => {
    const { draft, queued } = setup();
    await loaded(draft);
    draft.setField('-input_type', 'blastdb');
    draft.setField('-dbtype', 'nucl');
    draft.setField('-out', 'merged');
    const inOf = () => draft.previewCommandWords().slice(draft.previewCommandWords().indexOf('-in') + 1)[0];
    draft.setDbIn('db1 db2');
    expect(inOf()).toBe('db1 db2');
    draft.setDbIn('"my db" db2');
    expect(inOf()).toBe('"my db" db2');
    const res = await draft.submit();
    expect(res.ok).toBe(true);
    expect(queued[0]!.argv.slice(0, 3)).toEqual(['makeblastdb', '-in', '"my db" db2']);

    const found = setup().draft;
    await loaded(found);
    found.setField('-input_type', 'blastdb');
    found.addRunFiles(['nin', 'nhr', 'nsq'].map((ext) => new File(['x'], `my db.${ext}`)));
    expect(found.state.get().dbIn).toBe('"my db"');
    expect(found.previewCommandWords().slice(1, 4)).toEqual(['makeblastdb', '-in', '"my db"']);
  });

  it('adds mask data to -mask_data as a list, once per file, and takes it away with the file', async () => {
    const { draft } = setup();
    const a = new File(['a'], 'wm.asnb');
    const b = new File(['b'], 'seg.asnb');
    draft.addMaskData(a, 'wm.asnb', { kind: 'made', runNumber: 2 });
    draft.addMaskData(a, 'wm.asnb', { kind: 'made', runNumber: 2 });
    draft.addMaskData(b, 'seg.asnb', { kind: 'made', runNumber: 3 });
    let state = draft.state.get();
    expect(state.values['-mask_data']).toBe('wm.asnb,seg.asnb');
    expect(state.runFiles.map((entry) => [entry.name, entry.flag])).toEqual([
      ['wm.asnb', '-mask_data'],
      ['seg.asnb', '-mask_data'],
    ]);
    await loaded(draft);
    expect(changedMakeblastdbOptionsCount(state.values, state.extra, draft.state.get().description)).toBe(1);
    draft.removeRunFile(state.runFiles[0]!.key);
    state = draft.state.get();
    expect(state.values['-mask_data']).toBe('seg.asnb');
    draft.addMaskData(new File(['c'], 'wm.asnb'), 'wm.asnb', { kind: 'made', runNumber: 4 });
    draft.addMaskData(new File(['d'], 'seg.asnb'), 'seg.asnb', { kind: 'made', runNumber: 5 });
    state = draft.state.get();
    expect(state.runFilesError).toBeUndefined();
    expect(state.values['-mask_data']).toBe('seg.asnb,wm.asnb');
    expect(state.runFiles.map((entry) => [entry.name, entry.origin])).toEqual([
      ['seg.asnb', { kind: 'made', runNumber: 5 }],
      ['wm.asnb', { kind: 'made', runNumber: 4 }],
    ]);
    draft.resetOptions();
    state = draft.state.get();
    expect(state.values['-mask_data']).toBeUndefined();
    expect(state.runFiles).toHaveLength(0);
  });
  it("replaces the mask data of an earlier run handed over under the same name, in its place (WW4b2 R1 finding 4)", async () => {
    const { draft, queued } = setup();
    await loaded(draft);
    draft.setFastaPaste('>p\nMKV\n');
    draft.setField('-dbtype', 'prot');
    draft.addMaskData(new File(['a'], 'wm.asnb'), 'wm.asnb', { kind: 'made', runNumber: 1 });
    draft.addMaskData(new File(['s'], 'seg.asnb'), 'seg.asnb', { kind: 'made', runNumber: 2 });
    const later = new File(['bb'], 'wm.asnb');
    draft.addMaskData(later, 'wm.asnb', { kind: 'made', runNumber: 3 });
    const state = draft.state.get();
    expect(state.values['-mask_data']).toBe('wm.asnb,seg.asnb');
    expect(state.runFiles.map((entry) => [entry.name, entry.size, entry.origin])).toEqual([
      ['wm.asnb', 2, { kind: 'made', runNumber: 3 }],
      ['seg.asnb', 1, { kind: 'made', runNumber: 2 }],
    ]);
    expect(state.runFilesError).toBeUndefined();
    expect(state.message).toBeUndefined();
    const res = await draft.submit();
    expect(res.ok).toBe(true);
    expect(queued[0]!.directory.find((entry) => entry.name === 'wm.asnb')?.file).toBe(later);
  });

  it('refuses mask data under the name of a file that the user opened, in the message at the top of the tab (WW4b2 R1 finding 4)', async () => {
    const { draft } = setup();
    await loaded(draft);
    draft.addRunFiles([new File(['mine'], 'wm.asnb')]);
    draft.addMaskData(new File(['a'], 'wm.asnb'), 'wm.asnb', { kind: 'made', runNumber: 1 });
    const state = draft.state.get();
    expect(state.message).toEqual({ kind: 'error', text: 'A file named "wm.asnb" is already in the run directory.' });
    expect(state.runFilesError).toBeUndefined();
    expect(state.values['-mask_data'] ?? '').toBe('');
    expect(state.runFiles.map((entry) => [entry.name, entry.origin])).toEqual([['wm.asnb', undefined]]);
  });

  it('refuses mask data whose name holds a comma, which -mask_data cannot list, and says why (WW4b2 R1 finding 6)', async () => {
    const { draft } = setup();
    await loaded(draft);
    draft.addMaskData(new File(['a'], 'a,b.asnb'), 'a,b.asnb', { kind: 'made', runNumber: 1 });
    let state = draft.state.get();
    expect(state.message).toEqual({ kind: 'error', text: 'A file named "a,b.asnb" cannot be in -mask_data: a comma separates its files.' });
    expect(state.values['-mask_data'] ?? '').toBe('');
    expect(state.runFiles).toHaveLength(0);
    draft.setOptionFile('-mask_data', new File(['a'], 'a,b.asnb'));
    state = draft.state.get();
    expect(state.runFilesError).toBe('A file named "a,b.asnb" cannot be in -mask_data: a comma separates its files.');
    expect(state.runFiles).toHaveLength(0);
  });
});

describe('the results of a makeblastdb run', () => {
  function makeCoordinator() {
    let tok = 0;
    const data = new DataService({
      store: new MemoryBlockStore(),
      scanner: new FakeScanner(),
      inflater: new NodeInflater(),
      checker: new FakeInputChecker(),
      digest: sha256Hex,
      newToken: () => `t-${++tok}`,
      cleanup: Promise.resolve({ state: 'done', removedSessions: 0 }),
    });
    const saved: { readonly name: string; readonly blob: Blob }[] = [];
    let id = 0;
    const coordinator = new Coordinator({
      engine: new FakeEngine(),
      data,
      downloader: { save: () => undefined, saveBlob: (name, blob) => saved.push({ name, blob }) },
      now: () => 10_000,
      newRunId: () => `run-${++id}`,
      digest: sha256Hex,
    });
    return { coordinator, saved };
  }

  async function finished(coordinator: Coordinator, request: ToolRequest): Promise<string> {
    const result = await coordinator.enqueueToolAll([request]);
    if (!result.ok) throw new Error(result.message);
    const runId = result.runIds![0]!;
    await vi.waitFor(() => {
      const run = coordinator.state.get().runs.find((entry) => entry.snapshot.runId === runId);
      expect(run?.status).toBe('completed');
    }, { timeout: 30_000 });
    return runId;
  }

  it('hands its database to a search, and downloads all its files as one tar', async () => {
    const { coordinator, saved } = makeCoordinator();
    const runId = await finished(coordinator, {
      tool: 'makeblastdb',
      stage: 'build',
      argv: ['makeblastdb', '-in', 'input.fa', '-dbtype', 'prot', '-out', 'db/mydb'],
      directory: [{ name: 'input.fa', file: new Blob(['>p\nMKV\n']), origin: { kind: 'pasted' } }],
    });
    const made = await coordinator.madeDatabase(runId);
    expect(made.kind).toBe('protein');
    expect([...made.names].sort()).toEqual(['run1/db/mydb.phr', 'run1/db/mydb.pin', 'run1/db/mydb.psq']);
    expect(made.files.map((file) => file.name).sort()).toEqual(['mydb.phr', 'mydb.pin', 'mydb.psq']);

    await coordinator.exportToolFiles(runId);
    expect(saved).toHaveLength(1);
    expect(saved[0]!.name).toMatch(/db_mydb\.tar$/);
    const tar = new Uint8Array(await saved[0]!.blob.arrayBuffer());
    expect(tar.length % 512).toBe(0);
    expect(new TextDecoder().decode(tar.subarray(0, 100)).replace(/\0+$/, '')).toMatch(/^db\/mydb\.p(hr|in|sq)$/);
  });

  it("hands two runs' databases of one name to one search in two folders, one per run (WW4b2 R1 finding 3)", async () => {
    const { coordinator } = makeCoordinator();
    const request = (sequence: string): ToolRequest => ({
      tool: 'makeblastdb',
      stage: 'build',
      argv: ['makeblastdb', '-in', 'input.fa', '-dbtype', 'nucl'],
      directory: [{ name: 'input.fa', file: new Blob([sequence]), origin: { kind: 'pasted' } }],
    });
    const first = await coordinator.madeDatabase(await finished(coordinator, request('>a\nACGT\n')));
    const second = await coordinator.madeDatabase(await finished(coordinator, request('>b\nGGCC\n')));
    expect([...first.names].sort()).toEqual(['run1/input.fa.nhr', 'run1/input.fa.nin', 'run1/input.fa.nsq']);
    expect([...second.names].sort()).toEqual(['run2/input.fa.nhr', 'run2/input.fa.nin', 'run2/input.fa.nsq']);
    expect(madeDatabaseFolder(2)).toBe('run2');
    const directory = assembleSearchDirectory({
      subjectDb: { discovery: await discoverBlastDatabases(first.files, first.names), selectedName: 'run1/input.fa' },
      filteringDb: { discovery: await discoverBlastDatabases(second.files, second.names), selectedName: 'run2/input.fa' },
    });
    expect(directory.files.map((file) => file.name).sort()).toEqual([...first.names, ...second.names].sort());
  });

  it("hands a WindowMasker run's maskinfo output to makeblastdb's -mask_data", async () => {
    const { coordinator } = makeCoordinator();
    const runId = await finished(coordinator, {
      tool: 'windowmasker',
      stage: 'mask',
      argv: ['windowmasker', '-ustat', 'c.obinary', '-in', 'a.fa', '-outfmt', 'maskinfo_asn1_bin', '-out', 'wm.asnb'],
      directory: [
        { name: 'c.obinary', file: new Blob(['counts']), origin: { kind: 'picked' } },
        { name: 'a.fa', file: new Blob(['>a\nACGT\n']), origin: { kind: 'picked' } },
      ],
    });
    const run = coordinator.state.get().runs.find((entry) => entry.snapshot.runId === runId)!;
    const key = (run.record as { files?: readonly { key: string; name: string }[] }).files!.find((file) => file.name === 'wm.asnb')!.key;
    const handed = await coordinator.maskData(runId, key);
    expect(handed.name).toBe('wm.asnb');
    expect(handed.origin).toEqual({ kind: 'made', runNumber: run.snapshot.number });
    await expect(coordinator.maskData(runId, 'stderr')).rejects.toThrow('the output is not mask data');
  });
});
