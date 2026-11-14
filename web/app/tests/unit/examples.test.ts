import { gzipSync } from 'node:zlib';
import { describe, expect, it } from 'vitest';
import { sha256Hex } from 'losat';
import type { WholeFileInflater } from 'losat/ports';
import { siteHeaders } from '../../build/headers';
import { SearchDraft } from '../../src/application/draft';
import { ExampleLoader } from '../../src/application/examples';
import { EXAMPLE_BYTE_LIMIT, EXAMPLE_SETS, exampleById, exampleHost, exampleQuery, examplesFor } from '../../src/domain/examples';
import { PROGRAMS, programById, sequenceKind } from '../../src/domain/programs';
import { BrowserFileDownloader } from '../../src/infra/browser/example-downloader';
import { DataService } from '../../src/infra/data/data-service';
import { inflateFile } from '../../src/infra/data/inflate';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeDeletionChannel, FakeFileDownloader, FakeFileSaver, MemoryExampleStore } from '../../src/infra/fake/fake-examples';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeAccessionFetcher } from '../../src/infra/fake/fake-fetcher';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import type { DownloadProgress, DownloadResult, ExampleDeletions, ExampleStore } from '../../src/ports/example';
import { NodeInflater } from './support/node-inflater';

const SIXTEEN_S = exampleById('refseq-16s-bacteria');
const UNIPROT = exampleById('uniprot-ecoli-k12');
const SET_16S = '>NR_102783.2 Bacillus subtilis 16S\nACGTACGTAC\n>NR_114042.1 Escherichia coli 16S\nACGTTCGTAC\n';

function setup(
  options: {
    files?: Record<string, string | DownloadResult>;
    persistent?: boolean;
    store?: MemoryExampleStore;
    deletions?: ExampleDeletions;
    wrapStore?: (store: MemoryExampleStore) => ExampleStore;
    saver?: FakeFileSaver;
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
  const fake = new FakeEngine();
  const draft = new SearchDraft({
    engine: { describe: fake.describe.bind(fake), validate: async () => ({ ok: true }), prepare: () => undefined },
    data,
    enqueueAll: async (batch) => ({ ok: true, runIds: batch.map((_, i) => `r${i}`) }),
    fetcher: new FakeAccessionFetcher(),
    now: () => 0,
    debounceMs: 0,
  });
  const downloader = new FakeFileDownloader(options.files ?? { [SIXTEEN_S.url]: SET_16S });
  const store = options.store ?? new MemoryExampleStore({ persistent: options.persistent ?? true, now: () => Date.UTC(2026, 10, 3) });
  const examples = new ExampleLoader({
    downloader,
    store: options.wrapStore?.(store) ?? store,
    draft,
    ...(options.deletions === undefined ? {} : { deletions: options.deletions }),
    ...(options.saver === undefined ? {} : { saver: options.saver }),
  });
  const settle = async () => {
    await examples.idle();
    await draft.idle();
  };
  return { draft, examples, downloader, store, settle };
}

describe('Example table', () => {
  it('gives every mode a set, and each set a query of the kind that its modes read', () => {
    for (const program of PROGRAMS) {
      const sets = examplesFor(program.id);
      expect(sets.length, program.id).toBeGreaterThan(0);
      for (const set of sets) {
        const query = exampleQuery(set, program.id);
        expect(query.startsWith('>'), set.id).toBe(true);
        expect(query.endsWith('\n'), set.id).toBe(true);
        const letters = query.split('\n').slice(1).join('');
        const nucleotide = /^[ACGT]+$/.test(letters);
        expect(nucleotide, `${set.id} ${program.id}`).toBe(sequenceKind(programById(program.id), 'query') === 'nucleotide');
      }
    }
  });

  it('downloads only from hosts that the CSP names, over HTTPS, without anything of the user in the URL', () => {
    const csp = siteHeaders()['Content-Security-Policy'] ?? '';
    const connect = csp.split(';').find((directive) => directive.trim().startsWith('connect-src'))!.trim().split(/\s+/);
    for (const set of EXAMPLE_SETS) {
      expect(new URL(set.url).protocol).toBe('https:');
      expect(connect, set.id).toContain(`https://${exampleHost(set)}`);
      expect(set.fileName.endsWith('.gz')).toBe(false);
    }
  });

  it("carries a recA CDS that translates to the RecA protein query (one story: B. subtilis 168)", () => {
    const cds = UNIPROT.queries.nucleotide!.split('\n').slice(1).join('');
    const protein = UNIPROT.queries.protein!.split('\n').slice(1).join('');
    const bases = 'TCAG';
    const amino = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
    let translated = '';
    for (let i = 0; i + 3 <= cds.length; i += 3) {
      const [a, b, c] = [...cds.slice(i, i + 3)].map((base) => bases.indexOf(base));
      translated += amino[a! * 16 + b! * 4 + c!];
    }
    // The CDS starts with TTG (Met as a start codon) and ends with a stop.
    expect(translated.endsWith('*')).toBe(true);
    expect(translated.slice(1, -1)).toBe(protein.slice(1));
  });
});

describe('ExampleLoader', () => {
  it('downloads a set only when asked, saves it, and puts the query and the subject set into the draft', async () => {
    const { draft, examples, downloader, store, settle } = setup();
    draft.setPaste('query', '>old\nACGT\n');
    draft.addFiles('subject', [new File(['>s\nACGT\n'], 'mine.fa')]);
    examples.refresh();
    await settle();
    expect(downloader.requests).toEqual([]);
    expect(examples.state.get().saved).toEqual({ 'refseq-16s-bacteria': null, 'refseq-its-fungi': null });

    examples.load('refseq-16s-bacteria');
    await settle();
    expect(downloader.requests).toEqual([{ url: SIXTEEN_S.url, maxBytes: EXAMPLE_BYTE_LIMIT, fileName: SIXTEEN_S.fileName }]);
    const state = draft.state.get();
    expect(state.query.paste).toBe(exampleQuery(SIXTEEN_S, 'blastn'));
    expect(state.subject.paste).toBe('');
    expect(state.subject.sources.map((s) => [s.name, s.example, s.status])).toEqual([
      ['bacteria.16SrRNA.fna', 'refseq-16s-bacteria', 'ready'],
    ]);
    expect(state.subject.sources[0]!.base?.records.map((r) => r.id)).toEqual(['NR_102783.2', 'NR_114042.1']);
    expect(examples.state.get().load).toBeUndefined();
    // The screen tells a screen reader that the set was put into the inputs (WA2-AX-03).
    expect(examples.state.get().loaded).toEqual({ id: 'refseq-16s-bacteria', count: 1 });
    expect(examples.state.get().saved['refseq-16s-bacteria']).toEqual({
      id: 'refseq-16s-bacteria',
      url: SIXTEEN_S.url,
      version: '2026-10-27',
      downloadedBytes: SET_16S.length,
      bytes: SET_16S.length,
      savedAt: Date.UTC(2026, 10, 3),
    });
    expect(store.copies.size).toBe(1);
  });

  it('reads a saved copy without a request', async () => {
    const { draft, examples, downloader, settle } = setup();
    examples.load('refseq-16s-bacteria');
    await settle();
    draft.removeSource('subject', draft.state.get().subject.sources[0]!.key);
    examples.load('refseq-16s-bacteria');
    await settle();
    expect(downloader.requests).toHaveLength(1);
    expect(draft.state.get().subject.sources.map((s) => s.name)).toEqual(['bacteria.16SrRNA.fna']);
    expect(examples.state.get().loaded).toEqual({ id: 'refseq-16s-bacteria', count: 2 });
  });

  it('shows the progress, and Cancel leaves the draft as it was', async () => {
    const { draft, examples, downloader, settle } = setup();
    draft.setPaste('query', '>old\nACGT\n');
    downloader.hold();
    examples.load('refseq-16s-bacteria');
    await new Promise((resolve) => setTimeout(resolve, 0));
    expect(examples.state.get().load).toEqual({
      id: 'refseq-16s-bacteria',
      progress: { state: 'downloading', bytes: Math.floor(SET_16S.length / 2), total: SET_16S.length },
    });
    examples.cancel();
    downloader.release();
    await settle();
    expect(examples.state.get().load).toBeUndefined();
    expect(examples.state.get().loaded).toBeUndefined();
    expect(draft.state.get().query.paste).toBe('>old\nACGT\n');
    expect(draft.state.get().subject.sources).toEqual([]);
  });

  it('names the site and the status of a failure, and leaves the draft as it was', async () => {
    const { draft, examples, settle } = setup({
      files: {
        [UNIPROT.url]: { ok: false, reason: 'refused', status: 500 },
        [exampleById('pdbaa').url]: { ok: false, reason: 'unreachable', detail: 'Failed to fetch' },
      },
    });
    draft.setProgram('blastp');
    draft.setPaste('query', '>old\nMKV\n');
    examples.load('uniprot-ecoli-k12');
    await settle();
    expect(examples.state.get().load).toEqual({
      id: 'uniprot-ecoli-k12',
      progress: { state: 'failed', message: 'rest.uniprot.org answered HTTP 500: try again in a moment.' },
    });
    expect(draft.state.get().query.paste).toBe('>old\nMKV\n');
    examples.load('pdbaa');
    await settle();
    expect(examples.state.get().load?.progress).toEqual({
      state: 'failed',
      message: 'ftp.ncbi.nlm.nih.gov could not be reached (Failed to fetch).',
    });
  });

  it('uses the set in this tab when the browser keeps nothing, or when the storage is full', async () => {
    const none = setup({ persistent: false });
    none.examples.load('refseq-16s-bacteria');
    await none.settle();
    expect(none.draft.state.get().subject.sources.map((s) => s.name)).toEqual(['bacteria.16SrRNA.fna']);
    expect(none.examples.state.get().saved['refseq-16s-bacteria']).toBeNull();
    expect(none.examples.state.get().note).toBeUndefined();

    const full = setup();
    full.store.saveError = new Error('The browser refused to keep the example: its storage is full.');
    full.examples.load('refseq-16s-bacteria');
    await full.settle();
    expect(full.draft.state.get().subject.sources.map((s) => s.status)).toEqual(['ready']);
    expect(full.examples.state.get().note).toEqual({
      id: 'refseq-16s-bacteria',
      text: 'Not kept for the next time: The browser refused to keep the example: its storage is full.',
    });
  });

  it('takes the query of the mode: the protein for BLASTP, the CDS for BLASTX', async () => {
    const { draft, examples, settle } = setup({ files: { [UNIPROT.url]: '>sp|P0A7G6|RECA_ECOLI\nMAIDENKQ\n' } });
    draft.setProgram('blastx');
    examples.load('uniprot-ecoli-k12');
    await settle();
    expect(draft.state.get().query.paste).toBe(UNIPROT.queries.nucleotide);
    draft.setProgram('blastp');
    examples.load('uniprot-ecoli-k12');
    await settle();
    expect(draft.state.get().query.paste).toBe(UNIPROT.queries.protein);
  });

  it('ends a load when the mode has no button for its set, and never fills another mode with it', async () => {
    const { draft, examples, downloader, settle } = setup();
    downloader.hold();
    examples.load('refseq-16s-bacteria');
    await new Promise((resolve) => setTimeout(resolve, 0));
    draft.setProgram('blastp');
    draft.setPaste('query', '>mine\nMKV\n');
    examples.refresh();
    expect(examples.state.get().load).toBeUndefined();
    downloader.release();
    await settle();
    expect(draft.state.get().query.paste).toBe('>mine\nMKV\n');
    expect(draft.state.get().subject.sources).toEqual([]);

    // Without the refresh (the line did not follow the mode), the set is saved but the draft is not touched.
    downloader.hold();
    draft.setProgram('blastn');
    examples.load('refseq-16s-bacteria');
    await new Promise((resolve) => setTimeout(resolve, 0));
    draft.setProgram('blastp');
    downloader.release();
    await settle();
    expect(draft.state.get().query.paste).toBe('>mine\nMKV\n');
    expect(examples.state.get().saved['refseq-16s-bacteria']).not.toBeNull();
    expect(examples.state.get().load).toBeUndefined();
  });

  it('takes the query of the mode chosen when the set arrives (BLASTP to BLASTX during the download)', async () => {
    const { draft, examples, downloader, settle } = setup({ files: { [UNIPROT.url]: '>sp|P0A7G6|RECA_ECOLI\nMAIDENKQ\n' } });
    draft.setProgram('blastp');
    downloader.hold();
    examples.load('uniprot-ecoli-k12');
    await new Promise((resolve) => setTimeout(resolve, 0));
    draft.setProgram('blastx');
    examples.refresh();
    downloader.release();
    await settle();
    expect(draft.state.get().query.paste).toBe(UNIPROT.queries.nucleotide);
  });

  it('keeps an answer that is not FASTA out of the store', async () => {
    const { draft, examples, store, settle } = setup({ files: { [SIXTEEN_S.url]: '<html><body>Not here</body></html>' } });
    examples.load('refseq-16s-bacteria');
    await settle();
    expect(examples.state.get().load?.progress).toEqual({
      state: 'failed',
      message: 'The file from ftp.ncbi.nlm.nih.gov is not FASTA.',
    });
    expect(store.copies.size).toBe(0);
    expect(draft.state.get().subject.sources).toEqual([]);
  });

  it('says a copy is saved when the load was cancelled while it was being saved', async () => {
    let resume = () => {};
    const gate = new Promise<void>((resolve) => (resume = resolve));
    const { draft, examples, settle } = setup({
      wrapStore: (store) => ({
        readExample: (id, url) => store.readExample(id, url),
        saveExample: async (copy, data) => {
          await gate;
          return store.saveExample(copy, data);
        },
        deleteExample: (id) => store.deleteExample(id),
      }),
    });
    examples.load('refseq-16s-bacteria');
    await new Promise((resolve) => setTimeout(resolve, 0));
    expect(examples.state.get().load?.progress).toEqual({ state: 'saving' });
    examples.cancel();
    resume();
    await settle();
    expect(examples.state.get().saved['refseq-16s-bacteria']).toMatchObject({ bytes: SET_16S.length });
    expect(draft.state.get().subject.sources).toEqual([]);
  });

  it('a load waits for a deletion of the same set, and other tabs drop the sources of a deleted copy', async () => {
    const channel = new FakeDeletionChannel();
    const store = new MemoryExampleStore();
    const a = setup({ store, deletions: channel.tab() });
    const b = setup({ store, deletions: channel.tab() });
    a.examples.load('refseq-16s-bacteria');
    await a.settle();
    expect(a.draft.state.get().subject.sources.map((s) => s.example)).toEqual(['refseq-16s-bacteria']);

    // Tab B deletes, and loads again at once: the load does not get the copy being deleted.
    const deletion = b.examples.delete('refseq-16s-bacteria');
    b.examples.load('refseq-16s-bacteria');
    await deletion;
    await b.settle();
    expect(b.downloader.requests).toHaveLength(1);
    expect(b.draft.state.get().subject.sources.map((s) => s.example)).toEqual(['refseq-16s-bacteria']);
    // Tab A heard the deletion: its source of the deleted copy is gone.
    expect(a.draft.state.get().subject.sources).toEqual([]);
    expect(a.examples.state.get().saved['refseq-16s-bacteria']).toBeNull();
  });

  it('Delete removes the saved copy and the subject source made from it', async () => {
    const { draft, examples, store, settle } = setup();
    examples.load('refseq-16s-bacteria');
    await settle();
    draft.addFiles('subject', [new File(['>s\nACGT\n'], 'mine.fa')]);
    await settle();
    await examples.delete('refseq-16s-bacteria');
    expect(store.copies.size).toBe(0);
    expect(examples.state.get().saved['refseq-16s-bacteria']).toBeNull();
    expect(draft.state.get().subject.sources.map((s) => s.name)).toEqual(['mine.fa']);
  });
});

describe('Example copies outside the browser (line WF-3)', () => {
  it('writes the saved copy into the folder picked, named as the subject source, without a request', async () => {
    const saver = new FakeFileSaver(true);
    const { examples, downloader, settle } = setup({ saver });
    expect(examples.copyTo).toBe('folder');
    examples.load('refseq-16s-bacteria');
    await settle();
    examples.saveCopy('refseq-16s-bacteria');
    expect(examples.state.get().copying).toBe('refseq-16s-bacteria');
    await settle();
    expect(saver.written).toEqual([{ folder: 'Examples', fileName: 'bacteria.16SrRNA.fna', text: SET_16S }]);
    expect(examples.state.get().copying).toBeUndefined();
    expect(examples.state.get().note).toEqual({ id: 'refseq-16s-bacteria', text: 'Saved bacteria.16SrRNA.fna in Examples.' });
    expect(downloader.requests).toHaveLength(1);
    // The picker closed: nothing is written and nothing is said.
    saver.folder = undefined;
    examples.saveCopy('refseq-16s-bacteria');
    await settle();
    expect(saver.written).toHaveLength(1);
    expect(examples.state.get().note).toBeUndefined();
  });

  it('downloads the saved copy where the browser picks no folder, and says when the copy is gone', async () => {
    const saver = new FakeFileSaver(false);
    const { examples, store, settle } = setup({ saver });
    expect(examples.copyTo).toBe('download');
    examples.load('refseq-16s-bacteria');
    await settle();
    examples.saveCopy('refseq-16s-bacteria');
    await settle();
    expect(saver.downloads.map((d) => d.fileName)).toEqual(['bacteria.16SrRNA.fna']);
    expect(await saver.downloads[0]!.file.text()).toBe(SET_16S);
    expect(examples.state.get().note).toBeUndefined();
    store.copies.clear();
    examples.saveCopy('refseq-16s-bacteria');
    await settle();
    expect(saver.downloads).toHaveLength(1);
    expect(examples.state.get().note).toEqual({ id: 'refseq-16s-bacteria', text: 'Not saved: the copy in this browser is gone' });
    expect(examples.state.get().saved['refseq-16s-bacteria']).toBeNull();
  });

  it('offers nothing without a saver', () => {
    expect(setup().examples.copyTo).toBeUndefined();
  });
});

describe('BrowserFileDownloader', () => {
  const TEXT = '>a\nACGT\n'.repeat(1000);
  /**
   * The Data worker's inflation as the page sees it: the path of a source (inflate.ts) over the
   * Node double of the engine's decoders; the files it was given are recorded.
   */
  const recording = () => {
    const given: Array<[string, number | undefined]> = [];
    const inflater: WholeFileInflater = {
      inflateFile: (file, maxBytes) => {
        given.push([file.name, maxBytes]);
        return inflateFile(file, new NodeInflater(), maxBytes);
      },
    };
    return { inflater, given };
  };
  const run = async (
    respond: () => Response | Promise<Response>,
    options: { maxBytes?: number; signal?: AbortSignal; inflater?: WholeFileInflater } = {},
  ) => {
    const calls: Array<[string, RequestInit | undefined]> = [];
    const progress: DownloadProgress[] = [];
    const { inflater, given } = recording();
    const downloader = new BrowserFileDownloader({
      inflater: options.inflater ?? inflater,
      fetch: async (input, init) => {
        calls.push([String(input), init]);
        return respond();
      },
    });
    const result = await downloader.download(
      { url: 'https://ftp.ncbi.nlm.nih.gov/x.fna.gz', maxBytes: options.maxBytes ?? EXAMPLE_BYTE_LIMIT, fileName: 'x.fna' },
      (p) => progress.push(p),
      options.signal ?? new AbortController().signal,
    );
    return { result, calls, progress, given };
  };
  const text = async (result: DownloadResult) => {
    if (!result.ok) throw new Error(`expected a file, got ${JSON.stringify(result)}`);
    return result.data.text();
  };

  it("inflates a gzip body with the engine's path in the Data worker, reports the bytes against Content-Length, and sends no cookie or referrer", async () => {
    const gz = gzipSync(TEXT);
    const { result, calls, progress, given } = await run(
      () =>
        new Response(gz, {
          headers: { 'Content-Length': String(gz.length), 'Last-Modified': 'Tue, 27 Oct 2026 12:04:14 GMT' },
        }),
    );
    expect(await text(result)).toBe(TEXT);
    expect(result).toMatchObject({ ok: true, downloadedBytes: gz.length, version: '2026-10-27' });
    // The body is the set's file as its site names it.
    expect(given).toEqual([['x.fna.gz', EXAMPLE_BYTE_LIMIT]]);
    expect(progress.at(-1)).toEqual({ bytes: gz.length, total: gz.length });
    expect(calls[0]![1]).toMatchObject({ credentials: 'omit', referrerPolicy: 'no-referrer' });
  });

  it('keeps a body that is not gzip as it is, without the inflation, and names UniProt releases', async () => {
    const { result, progress, given } = await run(() => new Response(TEXT, { headers: { 'X-UniProt-Release': '2026_03' } }));
    expect(await text(result)).toBe(TEXT);
    expect(result).toMatchObject({ version: 'UniProt release 2026_03' });
    expect(progress.at(-1)).toEqual({ bytes: TEXT.length });
    expect(given).toEqual([]);
  });

  it('reports a status, a damaged gzip, a body too large and an unreachable site', async () => {
    expect((await run(() => new Response('Internal server error', { status: 500 }))).result).toEqual({
      ok: false,
      reason: 'refused',
      status: 500,
    });
    const gz = gzipSync(TEXT);
    const cut = (await run(() => new Response(gz.subarray(0, gz.length - 20)))).result;
    expect(cut).toMatchObject({ ok: false, reason: 'corrupt' });
    // The engine's sentence is the detail.
    const sentence = "'x.fna.gz' is named .gz but is not gzip data (incomplete deflate stream); LOSAT reads a file named .gz through gzip";
    const engine: WholeFileInflater = { inflateFile: () => Promise.reject(new Error(sentence)) };
    expect((await run(() => new Response(gz), { inflater: engine })).result).toEqual({ ok: false, reason: 'corrupt', detail: sentence });
    // Too large once inflated, and too large as it comes.
    expect((await run(() => new Response(gz), { maxBytes: 100 })).result).toEqual({ ok: false, reason: 'too-large' });
    const large = await run(() => new Response(TEXT), { maxBytes: 100 });
    expect(large.result).toEqual({ ok: false, reason: 'too-large' });
    expect(large.given).toEqual([]);
    const away = await run(() => Promise.reject(new TypeError('Failed to fetch')));
    expect(away.result).toEqual({ ok: false, reason: 'unreachable', detail: 'Failed to fetch' });
    const broken = new ReadableStream<Uint8Array>({
      start(controller) {
        controller.enqueue(gz.subarray(0, 100));
        controller.error(new TypeError('network error'));
      },
    });
    expect((await run(() => new Response(broken))).result).toEqual({ ok: false, reason: 'unreachable', detail: 'network error' });
  });

  it('rejects when the user cancels, during the download or during the inflation', async () => {
    const controller = new AbortController();
    const body = new ReadableStream<Uint8Array>({
      start(stream) {
        stream.enqueue(new TextEncoder().encode('>a\n'));
      },
    });
    const pending = run(
      () => {
        setTimeout(() => controller.abort(new DOMException('cancelled', 'AbortError')), 0);
        return new Response(body);
      },
      { signal: controller.signal },
    );
    await expect(pending).rejects.toThrow('cancelled');
    const later = new AbortController();
    const slow: WholeFileInflater = {
      inflateFile: async (file, maxBytes) => {
        later.abort(new DOMException('cancelled while inflating', 'AbortError'));
        return inflateFile(file, new NodeInflater(), maxBytes);
      },
    };
    await expect(run(() => new Response(gzipSync(TEXT)), { signal: later.signal, inflater: slow })).rejects.toThrow('cancelled while inflating');
  });
});
