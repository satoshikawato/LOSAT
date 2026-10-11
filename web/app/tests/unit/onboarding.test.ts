// S16 onboarding: the example files and self-test.json (domain/examples.ts), "Try an example"
// (application/examples.ts, SearchDraft.loadExample) and "Getting started"
// (application/onboarding.ts, infra/browser/local-settings.ts).
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { Examples, holdsExample } from '../../src/application/examples';
import { SearchDraft } from '../../src/application/draft';
import { GETTING_STARTED_KEY, GettingStarted } from '../../src/application/onboarding';
import { EXAMPLE_TOPICS, exampleFileName, parseSelfTestTable } from '../../src/domain/examples';
import { PROGRAMS, type InputRole, type ProgramId } from '../../src/domain/programs';
import { sha256Hex } from '../../src/infra/browser/platform';
import { browserLocalSettings } from '../../src/infra/browser/local-settings';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import type { ExampleFiles, ExampleSource } from '../../src/ports/examples';
import type { LocalSettings } from '../../src/ports/local-settings';
import TABLE from '../../src/examples/self-test.json';

const EXAMPLES = resolve(import.meta.dirname, '../../src/examples');
const ROLES: readonly InputRole[] = ['query', 'subject'];
const sha256 = (bytes: Uint8Array) => createHash('sha256').update(bytes).digest('hex');
const clone = <T>(value: T): T => JSON.parse(JSON.stringify(value)) as T;

describe('the examples and self-test.json', () => {
  const table = parseSelfTestTable(TABLE);

  it('has an example query and subject for every program, and a self-test case for each', () => {
    expect(table.cases.map((c) => c.program)).toEqual(PROGRAMS.map((p) => p.id));
    for (const testCase of table.cases) {
      for (const role of ROLES) {
        expect(testCase.inputs[role].file).toBe(exampleFileName(testCase.program, role));
        // The file in the release is the one the expected outputs were made from.
        expect(sha256(readFileSync(resolve(EXAMPLES, testCase.inputs[role].file))), `${testCase.program} ${role}`).toBe(testCase.inputs[role].sha256);
      }
      expect(testCase.argv).toEqual([testCase.program, '-query', `example-${testCase.program}-query.fasta`, '-subject', `example-${testCase.program}-subject.fasta`]);
    }
  });

  it('names the engine commit and the NCBI release, and no machine path (plan TD-11)', () => {
    expect(table.engine.commit).toMatch(/^[0-9a-f]{7,40}$/);
    expect(table.oracle).toMatch(/^NCBI BLAST\+ 2\.17\.0/);
    expect(JSON.stringify(TABLE)).not.toMatch(/\/home\/|\/Users\/|\\\\|\/tmp\//);
  });

  it('keeps the definition lines of the source records, a cut region said after them', () => {
    const lines = PROGRAMS.flatMap((p) => ROLES.flatMap((role) => readFileSync(resolve(EXAMPLES, exampleFileName(p.id, role)), 'utf8').split('\n')));
    const deflines = lines.filter((line) => line.startsWith('>'));
    expect(deflines.length).toBeGreaterThan(10);
    for (const line of deflines) expect(line).toMatch(/^>(AP\d{6}\.1|NZ_CP\d{6}\.1|BD[A-Z]\d{5}\.1) \S.*$/);
    for (const line of deflines.filter((l) => l.includes('genome'))) expect(line).toMatch(/ \[region \d+-\d+\]$/);
  });

  it('refuses a table that is not the format', () => {
    expect(() => parseSelfTestTable({ ...clone(TABLE), format: 2 })).toThrow('unknown format');
    const argv = clone(TABLE);
    argv.cases[0]!.argv.push('-evalue', '1');
    expect(() => parseSelfTestTable(argv)).toThrow("case 1's argv is not the search form's for its files");
    const digest = clone(TABLE);
    digest.cases[1]!.outputs['6'].sha256 = 'ABC';
    expect(() => parseSelfTestTable(digest)).toThrow("case 2's outfmt 6 SHA-256 is not valid");
    const twice = clone(TABLE);
    twice.cases[1] = twice.cases[0]!;
    expect(() => parseSelfTestTable(twice)).toThrow('a program has two cases');
    const program = clone(TABLE);
    (program.cases[0] as { program: string }).program = 'psiblast';
    expect(() => parseSelfTestTable(program)).toThrow("case 1's program is not one of LOSAT Web's");
  });
});

// --- "Try an example" ---------------------------------------------------------------------------

function exampleFiles(program: ProgramId): ExampleFiles {
  return {
    query: new File([readFileSync(resolve(EXAMPLES, exampleFileName(program, 'query')))], exampleFileName(program, 'query')),
    subject: new File([readFileSync(resolve(EXAMPLES, exampleFileName(program, 'subject')))], exampleFileName(program, 'subject')),
  };
}

function setup(source: ExampleSource = { files: async (program) => exampleFiles(program) }) {
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
  const draft = new SearchDraft({
    engine: { describe: fake.describe.bind(fake), validate: fake.validate.bind(fake) },
    data,
    enqueueAll: async () => ({ ok: true, runIds: [] }),
    now: () => 0,
    debounceMs: 0,
  });
  const examples = new Examples({ source, draft });
  return { draft, examples };
}

describe('Try an example', () => {
  it('puts the example files in both roles with the default settings, and searches nothing', async () => {
    const { draft, examples } = setup();
    draft.setProgram('blastn');
    draft.setPaste('query', '>mine\nACGT\n');
    draft.addFiles('subject', [new File(['>a\nACGT\n'], 'a.fa'), new File(['>b\nACGT\n'], 'b.fa')]);
    draft.setMode('subject', 'separate');
    draft.setThreads(2);
    draft.setTitle('My job');
    draft.setProgram('tblastn');
    draft.setField('-evalue', '1e-5');
    draft.setProgram('blastp');
    draft.setField('-evalue', '1e-3');
    await draft.idle();
    expect(examples.replaces('tblastn')).toEqual({ inputs: true, settings: true, title: true });

    await examples.load('tblastn');
    await draft.idle();
    const state = draft.state.get();
    expect(state.program).toBe('tblastn');
    for (const role of ROLES) {
      expect(state[role].paste).toBe('');
      expect(state[role].mode).toBe('combined');
      expect(state[role].sources.map((s) => [s.origin, s.name, s.status])).toEqual([['file', `example-tblastn-${role}.fasta`, 'ready']]);
    }
    expect(state.values.tblastn).toEqual({});
    // Another program's values stay as the user set them.
    expect(state.values.blastp).toEqual({ '-evalue': '1e-3' });
    expect(state.threads).toBe('auto');
    expect(state.title).toBe('');
    expect(draft.settings()).toEqual({ program: 'tblastn', options: [], threads: 'auto' });
    expect(holdsExample(state, 'tblastn')).toBe(true);
    expect(examples.state.get()).toEqual({
      loaded: {
        program: 'tblastn',
        text: `Example loaded: TBLASTN of example-tblastn-query.fasta against example-tblastn-subject.fasta, ${EXAMPLE_TOPICS.tblastn}. Choose Run LOSAT to search.`,
      },
    });
    // The user's own change ends the example's line.
    draft.setPaste('query', '>x\nMKV\n');
    expect(holdsExample(draft.state.get(), 'tblastn')).toBe(false);
  });

  it('needs no question for an empty form', () => {
    const { examples } = setup();
    expect(examples.replaces('blastn')).toBeUndefined();
  });

  it('asks before it replaces settings or the Job Title, even without inputs', () => {
    const { draft, examples } = setup();
    draft.setProgram('blastn');
    draft.setField('-evalue', '1e-5');
    expect(examples.replaces('blastn')).toEqual({ inputs: false, settings: true, title: false });
    // Another program's settings are not replaced by this program's example.
    expect(examples.replaces('blastp')).toBeUndefined();
    draft.setThreads(2);
    expect(examples.replaces('blastp')).toEqual({ inputs: false, settings: true, title: false });
    draft.setThreads('auto');
    draft.setMode('query', 'separate');
    expect(examples.replaces('blastp')).toEqual({ inputs: false, settings: true, title: false });
    draft.setMode('query', 'combined');
    draft.setTitle('My job');
    expect(examples.replaces('blastp')).toEqual({ inputs: false, settings: false, title: true });
  });

  it("does not load BLASTX's example while BLASTX cannot be searched, and says why", async () => {
    const files = vi.fn(async (program: ProgramId) => exampleFiles(program));
    const { draft, examples } = setup({ files });
    const reason = PROGRAMS.find((p) => p.id === 'blastx')!.unavailable!;
    expect(examples.unavailable('blastx')).toBe(reason);
    expect(examples.unavailable('blastn')).toBeUndefined();
    await examples.load('blastx');
    expect(files).not.toHaveBeenCalled();
    expect(examples.state.get()).toEqual({ error: reason });
    expect(draft.state.get().program).toBe('blastn');
  });

  it('leaves the form as it was when the files cannot be read', async () => {
    const { draft, examples } = setup({
      files: async () => {
        throw new Error('The example file example-blastp-query.fasta could not be read from the site (HTTP 404).');
      },
    });
    draft.setPaste('query', '>mine\nACGT\n');
    await examples.load('blastp');
    expect(examples.state.get()).toEqual({ error: 'The example file example-blastp-query.fasta could not be read from the site (HTTP 404).' });
    expect(draft.state.get().program).toBe('blastn');
    expect(draft.state.get().query.paste).toBe('>mine\nACGT\n');
  });
});

// --- "Getting started" ----------------------------------------------------------------------------

function memorySettings(): LocalSettings & { readonly values: Map<string, string> } {
  const values = new Map<string, string>();
  return {
    values,
    read: (key) => values.get(key),
    write: (key, value) => {
      values.set(key, value);
      return true;
    },
  };
}

describe('Getting started', () => {
  it('is open on the first visit, and closed on later visits once closed', () => {
    const settings = memorySettings();
    const first = new GettingStarted(settings);
    expect(first.state.get().open).toBe(true);
    first.close();
    expect(first.state.get().open).toBe(false);
    expect(settings.values.get(GETTING_STARTED_KEY)).toBe('closed');
    const later = new GettingStarted(settings);
    expect(later.state.get().open).toBe(false);
    // Opening it again from the header does not change what later visits start with.
    later.toggle();
    expect(later.state.get().open).toBe(true);
    expect(new GettingStarted(settings).state.get().open).toBe(false);
  });

  it('works where the browser keeps nothing: open on each visit, closed in the tab after a close', () => {
    const none: LocalSettings = { read: () => undefined, write: () => false };
    const tab = new GettingStarted(none);
    expect(tab.state.get().open).toBe(true);
    tab.close();
    expect(tab.state.get().open).toBe(false);
    expect(new GettingStarted(none).state.get().open).toBe(true);
  });
});

describe('the browser local settings', () => {
  afterEach(() => {
    vi.unstubAllGlobals();
  });

  it('reads and writes localStorage under losat-web: keys', () => {
    const values = new Map<string, string>();
    vi.stubGlobal('window', {
      localStorage: { getItem: (key: string) => values.get(key) ?? null, setItem: (key: string, value: string) => values.set(key, value) },
    });
    expect(browserLocalSettings.read(GETTING_STARTED_KEY)).toBeUndefined();
    expect(browserLocalSettings.write(GETTING_STARTED_KEY, 'closed')).toBe(true);
    expect(browserLocalSettings.read(GETTING_STARTED_KEY)).toBe('closed');
    expect(() => browserLocalSettings.read('getting-started')).toThrow('starts with losat-web:');
  });

  it('gives undefined and false where storage throws, even on the accessor', () => {
    vi.stubGlobal('window', {
      get localStorage(): Storage {
        throw new DOMException('The operation is insecure.', 'SecurityError');
      },
    });
    expect(browserLocalSettings.read(GETTING_STARTED_KEY)).toBeUndefined();
    expect(browserLocalSettings.write(GETTING_STARTED_KEY, 'closed')).toBe(false);
    vi.stubGlobal('window', {
      localStorage: {
        getItem: () => {
          throw new Error('denied');
        },
        setItem: () => {
          throw new DOMException('quota', 'QuotaExceededError');
        },
      },
    });
    expect(browserLocalSettings.read(GETTING_STARTED_KEY)).toBeUndefined();
    expect(browserLocalSettings.write(GETTING_STARTED_KEY, 'closed')).toBe(false);
  });
});
