// Settings files, "Edit Search" and the input FASTA of runs (S15 items 3 and 5): the draft's
// settings and their application to the form, the round trip of runs made by the form, and the
// files that RunFiles writes through the Writer contract.
import { describe, expect, it } from 'vitest';
import { Coordinator, type RunView } from '../../src/application/coordinator';
import { SearchDraft } from '../../src/application/draft';
import { INPUT_SLICE_BYTES, RunFiles } from '../../src/application/run-files';
import type { InputRole, ProgramId } from '../../src/domain/programs';
import type { RegionText } from '../../src/domain/region';
import { parseSettings, serializeSettings, settingsOfRun, SETTINGS_MAX_BYTES, type SearchSettings } from '../../src/domain/settings-file';
import { sha256Hex } from '../../src/infra/browser/platform';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import type { EngineGateway } from '../../src/ports/engine';
import { memoryDownloader, type SavedFile } from './support/memory-downloader';

function setup(options: { describe?: EngineGateway['describe']; maxThreads?: number } = {}) {
  let token = 0;
  const data = new DataService({
    store: new MemoryBlockStore(),
    scanner: new FakeScanner(),
    checker: new FakeInputChecker(),
    digest: sha256Hex,
    newToken: () => `token-${++token}`,
    cleanup: Promise.resolve({ state: 'done', removedSessions: 0 }),
  });
  const fake = new FakeEngine();
  const engine = { describe: options.describe ?? fake.describe.bind(fake), validate: fake.validate.bind(fake) };
  let id = 0;
  const coordinator = new Coordinator({
    engine: { ...engine, run: fake.run.bind(fake), cancel: fake.cancel.bind(fake) },
    data,
    downloader: memoryDownloader(() => undefined),
    now: () => 0,
    newRunId: () => `run-${++id}`,
  });
  const draft = new SearchDraft({ engine, data, enqueueAll: (requests) => coordinator.enqueueAll(requests), now: () => 0, debounceMs: 0 });
  const saved: SavedFile[] = [];
  const runFiles = new RunFiles({ draft, downloader: memoryDownloader((file) => saved.push(file)), maxThreads: () => options.maxThreads ?? 8 });
  return { draft, coordinator, runFiles, saved };
}

const file = (name: string, text: string | Uint8Array) => new File([text as BlobPart], name, { type: 'text/plain' });
const words = (draft: SearchDraft) => draft.parameters().flatMap(([flag, value]) => (value === true ? [flag] : [flag, value]));
const text = (saved: SavedFile) => new TextDecoder().decode(saved.bytes);
const NT = '>n1 one\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n';
const AA = '>p1 one\nMKVLAAGIVGLLLAHHKKEEDDPPWWRRSSMKVLAAGIVG\n';
const KINDS: Record<ProgramId, { query: string; subject: string }> = {
  blastn: { query: NT, subject: NT },
  blastp: { query: AA, subject: AA },
  blastx: { query: NT, subject: AA },
  tblastn: { query: AA, subject: NT },
  tblastx: { query: NT, subject: NT },
};

/** A draft of `program` with one record per role, the fields and the regions set, queued; its run. */
async function runOf(
  program: ProgramId,
  fields: ReadonlyArray<readonly [string, string | boolean]>,
  regions: Partial<Record<InputRole, RegionText>> = {},
  threads: number | 'auto' = 'auto',
): Promise<RunView> {
  const { draft, coordinator } = setup();
  draft.setProgram(program);
  draft.setPaste('query', KINDS[program].query);
  draft.addFiles('subject', [file('my subject.fa', KINDS[program].subject)]);
  await draft.idle();
  for (const [flag, value] of fields) draft.setField(flag, value);
  for (const [role, region] of Object.entries(regions)) draft.setRegion(role as InputRole, region);
  draft.setThreads(threads);
  const result = await draft.submit();
  expect(result.ok, result.ok ? '' : result.message).toBe(true);
  return coordinator.state.get().runs.at(-1)!;
}

describe('SearchDraft settings', () => {
  it('gives the program, the argv words after the inputs (regions included) and the threads', async () => {
    const { draft } = setup();
    draft.setPaste('query', NT);
    draft.setPaste('subject', NT);
    await draft.idle();
    draft.setField('-task', 'blastn');
    draft.setField('-penalty', '-3');
    draft.setField('-lcase_masking', true);
    draft.setRegion('query', { start: '2', stop: '30' });
    draft.setThreads(2);
    expect(draft.settings()).toEqual({
      program: 'blastn',
      options: ['-task', 'blastn', '-penalty', '-3', '-lcase_masking', '-query_loc', '2-30'],
      threads: 2,
    });
  });

  it('replaces the form conditions: fields the settings do not name return to their defaults; the inputs and the title stay', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>q1\nACGTACGT\n>q2\nACGT\n');
    draft.setPaste('subject', NT);
    draft.setTitle('Kept');
    await draft.idle();
    draft.setField('-evalue', '5');
    draft.setField('-word_size', '20');
    draft.setRegion('subject', { start: '3', stop: '9' });
    draft.setThreads(4);
    const before = { query: draft.state.get().query.sources, subject: draft.state.get().subject.sources };
    const applied = await draft.applySettings(
      { program: 'blastn', options: ['-task', 'blastn', '-penalty', '-3', '-lcase_masking'], threads: 'auto' },
      { maxThreads: 8 },
    );
    expect(applied).toEqual({ notApplied: [] });
    const state = draft.state.get();
    expect(state.values.blastn).toEqual({ '-task': 'blastn', '-penalty': '-3', '-lcase_masking': true });
    expect(state.threads).toBe('auto');
    expect(state.title).toBe('Kept');
    expect(state.query.sources).toEqual(before.query);
    expect(state.subject.sources).toEqual(before.subject);
    // The subject region was not named: it is cleared.
    expect(draft.region('subject')).toBeUndefined();
    expect(words(draft)).toEqual(['-task', 'blastn', '-penalty', '-3', '-lcase_masking']);
    await draft.idle();
    expect(draft.state.get().validation).toEqual({ state: 'ok' });
  });

  it('lists what it does not apply: words without a field, stray words, a region the inputs cannot take, too many threads', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>q1\nACGTACGT\n>q2\nACGT\n');
    draft.setPaste('subject', NT);
    await draft.idle();
    const applied = await draft.applySettings(
      {
        program: 'blastn',
        options: ['5', '-evalue', '1e-3', '-ungapped', '-foo', 'bar', '-query_loc', '2-5', '-subject_loc', '7-', '-perc_identity', '90'],
        threads: 12,
      },
      { maxThreads: 8, title: 'From a run' },
    );
    expect(applied.notApplied).toEqual([
      '5 (not an option)',
      '-ungapped (the BLASTN form has no field for it)',
      '-foo bar (the BLASTN form has no field for it)',
      '-query_loc 2-5 (a region needs a query of one record)',
      "-subject_loc 7- (the form's region is a start and a stop, start-stop)",
      'threads 12 (this browser offers 1 to 8; Auto is set)',
    ]);
    expect(draft.state.get().values.blastn).toEqual({ '-evalue': '1e-3', '-perc_identity': '90' });
    expect(draft.state.get().threads).toBe('auto');
    expect(draft.state.get().title).toBe('From a run');
    expect(words(draft)).toEqual(['-evalue', '1e-3', '-perc_identity', '90']);
  });

  it("switches the program, reads the inputs again with its reader, and then applies the regions on the inputs' only records", async () => {
    const { draft } = setup();
    draft.setPaste('query', AA);
    draft.setPaste('subject', NT);
    await draft.idle();
    const applied = await draft.applySettings(
      { program: 'tblastn', options: ['-db_gencode', '4', '-query_loc', '2-30', '-subject_loc', '3-40'], threads: 1 },
      { maxThreads: 8 },
    );
    expect(applied.notApplied).toEqual([]);
    expect(draft.state.get().program).toBe('tblastn');
    expect(words(draft)).toEqual(['-db_gencode', '4', '-query_loc', '2-30', '-subject_loc', '3-40']);
    expect(draft.state.get().threads).toBe(1);
  });

  it('changes nothing when the engine cannot describe the program', async () => {
    const fake = new FakeEngine();
    const { draft } = setup({
      describe: (program) => (program === 'blastp' ? Promise.reject(new Error('describe failed')) : fake.describe(program)),
    });
    await draft.idle();
    draft.setField('-evalue', '5');
    await expect(draft.applySettings({ program: 'blastp', options: [], threads: 'auto' }, { maxThreads: 8 })).rejects.toThrow('describe failed');
    expect(draft.state.get().program).toBe('blastn');
    expect(draft.state.get().values.blastn).toEqual({ '-evalue': '5' });
  });
});

describe('the round trip of runs made by the form', () => {
  const cases: ReadonlyArray<{
    readonly program: ProgramId;
    readonly fields: ReadonlyArray<readonly [string, string | boolean]>;
    readonly regions?: Partial<Record<InputRole, RegionText>>;
    readonly threads?: number | 'auto';
  }> = [
    { program: 'blastn', fields: [] },
    {
      program: 'blastn',
      fields: [
        ['-task', 'blastn'],
        ['-reward', '1'],
        ['-penalty', '-2'],
        ['-evalue', '1e-5'],
        ['-dust', '20 64 1'],
        ['-lcase_masking', true],
        ['-subject_besthit', true],
      ],
      regions: { query: { start: '2', stop: '30' }, subject: { start: '5', stop: '40' } },
      threads: 3,
    },
    {
      program: 'blastp',
      fields: [
        ['-matrix', 'BLOSUM45'],
        ['-gapopen', '15'],
        ['-gapextend', '2'],
        ['-comp_based_stats', '0'],
        ['-seg', 'yes'],
      ],
      regions: { query: { start: '3', stop: '33' } },
    },
    {
      program: 'tblastn',
      fields: [
        ['-task', 'tblastn-fast'],
        ['-db_gencode', '4'],
        ['-soft_masking', 'true'],
        ['-lcase_masking', true],
        ['-xdrop_gap', '20'],
      ],
      regions: { subject: { start: '1', stop: '39' } },
    },
    {
      program: 'tblastx',
      fields: [
        ['-query_gencode', '2'],
        ['-culling_limit', '1'],
        ['-seg', 'no'],
        ['-threshold', '14'],
      ],
      threads: 1,
    },
  ];

  for (const c of cases) {
    it(`${c.program} ${c.fields.map(([flag]) => flag).join(' ') || '(defaults)'}`, async () => {
      const run = await runOf(c.program, c.fields, c.regions, c.threads);
      // Read back through the file, into a form of another program with other values.
      const read = parseSettings(new TextEncoder().encode(serializeSettings(settingsOfRun(run.snapshot))));
      if (!read.ok) throw new Error(read.message);
      const { draft } = setup();
      draft.setPaste('query', KINDS[c.program].query);
      draft.setPaste('subject', KINDS[c.program].subject);
      await draft.idle();
      draft.setField('-evalue', '7');
      draft.setField('-word_size', '12');
      const applied = await draft.applySettings(read.settings, { maxThreads: 8 });
      expect(applied.notApplied).toEqual([]);
      expect(words(draft)).toEqual(run.snapshot.argv.slice(5));
      expect(draft.state.get().threads).toBe(run.snapshot.requestedThreads);
      expect(draft.state.get().program).toBe(c.program);
    });
  }
});

describe('RunFiles', () => {
  it('saves the form settings as a LOSAT Web settings file, written through the Writer', async () => {
    const { draft, runFiles, saved } = setup();
    draft.setField('-evalue', '1e-5');
    draft.setThreads(2);
    await runFiles.saveSettings();
    expect(saved.map((f) => [f.name, f.mime])).toEqual([['losat-settings-blastn.json', 'application/json']]);
    expect(parseSettings(saved[0]!.bytes)).toEqual({ ok: true, settings: { program: 'blastn', options: ['-evalue', '1e-5'], threads: 2 } });
    expect(text(saved[0]!)).not.toMatch(/query\.fa|subject\.fa/);
    expect(runFiles.state.get().settings).toEqual({ kind: 'info', text: 'Saved losat-settings-blastn.json: the program, the options and the threads.' });
  });

  it('loads a settings file into the form, and says what was loaded and what was not applied', async () => {
    const { draft, runFiles } = setup();
    draft.setTitle('Mine');
    const settings: SearchSettings = { program: 'tblastx', options: ['-seg', 'no', '-max_target_seqs', '5', '-matrix', 'PAM30'], threads: 'auto' };
    await runFiles.loadSettings(file('mine.json', serializeSettings(settings)));
    expect(draft.state.get().program).toBe('tblastx');
    expect(words(draft)).toEqual(['-max_target_seqs', '5', '-seg', 'no']);
    expect(draft.state.get().title).toBe('Mine');
    expect(runFiles.state.get().settings).toEqual({
      kind: 'info',
      text:
        'Loaded mine.json: TBLASTX, 6 words of options, threads Auto. The inputs and the Job Title did not change. ' +
        'Not applied: -matrix PAM30 (the TBLASTX form has no field for it).',
    });
  });

  it('refuses a broken, oversized or BLASTX settings file and leaves the form as it was', async () => {
    const { draft, runFiles } = setup();
    draft.setField('-evalue', '5');
    const before = draft.state.get().values;
    await runFiles.loadSettings(file('broken.json', '{"format": "LOSAT Web search settings", '));
    expect(runFiles.state.get().settings?.kind).toBe('error');
    expect(runFiles.state.get().settings?.text).toMatch(/^broken\.json was not loaded: The file is not JSON \(.+\)\.$/);
    // Larger than a settings file can be: refused from its size, before it is read.
    const large = file('large.json', new Uint8Array(SETTINGS_MAX_BYTES + 1));
    large.arrayBuffer = () => Promise.reject(new Error('read'));
    await runFiles.loadSettings(large);
    expect(runFiles.state.get().settings?.text).toBe(
      `large.json was not loaded: it has ${SETTINGS_MAX_BYTES + 1} bytes; a settings file has at most ${SETTINGS_MAX_BYTES} (1 MiB).`,
    );
    await runFiles.loadSettings(file('x.json', serializeSettings({ program: 'blastx', options: [], threads: 'auto' })));
    expect(runFiles.state.get().settings?.text).toMatch(/^x\.json was not loaded: its program is BLASTX\. BLASTX joins LOSAT Web after its certification/);
    await runFiles.loadSettings(file('newer.json', JSON.stringify({ format: 'LOSAT Web search settings', schema: 2 })));
    expect(runFiles.state.get().settings?.text).toBe('newer.json was not loaded: The file has schema 2: a newer LOSAT Web saved it. This one reads schema 1.');
    expect(draft.state.get().values).toBe(before);
    expect(draft.state.get().program).toBe('blastn');
  });

  it("saves a run's settings from its argv, and Edit Search puts them and its Job Title in the form, which keeps its inputs", async () => {
    const run = await runOf('blastp', [['-matrix', 'BLOSUM45'], ['-seg', 'yes']], { query: { start: '2', stop: '20' } }, 2);
    const { draft, runFiles, saved } = setup();
    await runFiles.saveRunSettings(run);
    expect(saved[0]!.name).toBe('losat-settings-blastp.json');
    expect(parseSettings(saved[0]!.bytes)).toEqual({ ok: true, settings: settingsOfRun(run.snapshot) });
    expect(runFiles.state.get().run).toEqual({ kind: 'info', text: 'Saved losat-settings-blastp.json: the settings of Run 1.', runId: run.snapshot.runId });

    draft.setPaste('query', '>other\nMKV\n>two\nMKV\n');
    draft.setTitle('Mine');
    await draft.idle();
    const query = draft.state.get().query;
    expect(await runFiles.editSearch({ ...run, snapshot: { ...run.snapshot, title: 'Run title' } })).toBe(true);
    expect(draft.state.get().program).toBe('blastp');
    expect(draft.state.get().title).toBe('Run title');
    expect(draft.state.get().threads).toBe(2);
    // The form's query has two records, so the run's query region does not apply.
    expect(words(draft)).toEqual(['-matrix', 'BLOSUM45', '-seg', 'yes']);
    expect(draft.state.get().query.sources.map((s) => s.name)).toEqual(query.sources.map((s) => s.name));
    expect(runFiles.state.get().settings).toEqual({
      kind: 'info',
      text: "The search form has the settings of Run 1. The inputs are the form's own. Not applied: -query_loc 2-20 (a region needs a query of one record).",
    });
    // A run without a title clears the form's title.
    await runFiles.editSearch(run);
    expect(draft.state.get().title).toBe('');
  });

  it('Edit Search changes nothing and says why when the engine cannot describe the program', async () => {
    const run = await runOf('blastp', [['-seg', 'yes']]);
    const fake = new FakeEngine();
    const { draft, runFiles } = setup({
      describe: (program) => (program === 'blastp' ? Promise.reject(new Error('describe failed')) : fake.describe(program)),
    });
    await draft.idle();
    expect(await runFiles.editSearch(run)).toBe(false);
    expect(draft.state.get().program).toBe('blastn');
    expect(runFiles.state.get().run).toEqual({
      kind: 'error',
      text: 'The settings of Run 1 could not be put in the form: describe failed',
      runId: run.snapshot.runId,
    });
    expect(runFiles.state.get().settings).toBeUndefined();
  });

  it('saves the exact bytes that the engine searched, under the argv name, in bounded blocks', async () => {
    const { draft, coordinator, runFiles, saved } = setup();
    draft.setPaste('query', NT);
    draft.addFiles('subject', [file('a b.fa', '>s1\nACGTACGTAC\n>s2\nACGTACGTAC\n>s3\nACGT\n'), file('c.fa', '>s4\nACGTACGT\n')]);
    await draft.idle();
    // The second record of the first file left out: the run input is not that file.
    draft.setIncluded('subject', draft.state.get().subject.sources[0]!.key, [1], false);
    await draft.submit();
    const run = coordinator.state.get().runs[0]!;
    expect(run.snapshot.subject.name).toBe('combined_subject.fa');
    await runFiles.saveInput(run, 'query');
    await runFiles.saveInput(run, 'subject');
    expect(saved.map((f) => [f.name, f.mime])).toEqual([
      ['query.fa', 'text/plain'],
      ['combined_subject.fa', 'text/plain'],
    ]);
    expect(saved[0]!.bytes).toEqual(run.snapshot.query.bytes);
    expect(saved[1]!.bytes).toEqual(run.snapshot.subject.bytes);
    expect(text(saved[1]!)).toBe('>s1\nACGTACGTAC\n>s3\nACGT\n>s4\nACGTACGT\n');
    expect(runFiles.state.get().run).toEqual({ kind: 'info', text: 'Saved combined_subject.fa: the subject that Run 1 searched.', runId: run.snapshot.runId });
    // Where the inputs came from: the pasted text whole, and a selection of a file joined with a whole file.
    expect(runFiles.inputParts(run, 'query')).toEqual([{ origin: 'paste', name: 'query.fa', records: 1 }]);
    expect(runFiles.inputParts(run, 'subject')).toEqual([undefined, { origin: 'file', name: 'c.fa', records: 1 }]);
  });

  it("writes a large input in slices of the Writer's block size", async () => {
    const run = await runOf('blastn', []);
    const { runFiles, saved } = setup();
    const big = new Uint8Array(2 * INPUT_SLICE_BYTES + 123).map((_, i) => 0x41 + (i % 4));
    const large: RunView = { ...run, snapshot: { ...run.snapshot, query: { ...run.snapshot.query, name: 'big.fa', bytes: big } } };
    await runFiles.saveInput(large, 'query');
    expect(saved[0]!.name).toBe('big.fa');
    expect(saved[0]!.blocks).toBe(3);
    // Compared as buffers: a deep comparison of 2 MiB element by element is slow.
    expect(Buffer.compare(Buffer.from(saved[0]!.bytes), Buffer.from(big))).toBe(0);
  });

  it('refuses the input of a run that has no copy of its bytes', async () => {
    const run = await runOf('blastn', []);
    const { runFiles, saved } = setup();
    // A run loaded from a session file has no engine bytes until its original FASTA is attached (WP-D).
    const { bytes, ...query } = run.snapshot.query;
    expect(bytes.length).toBeGreaterThan(0);
    await runFiles.saveInput({ ...run, snapshot: { ...run.snapshot, query: query as RunView['snapshot']['query'] } }, 'query');
    expect(saved).toEqual([]);
    expect(runFiles.state.get().run).toEqual({
      kind: 'error',
      text: 'query.fa was not saved: this working session has no copy of the query that Run 1 searched.',
      runId: run.snapshot.runId,
    });
  });
});
