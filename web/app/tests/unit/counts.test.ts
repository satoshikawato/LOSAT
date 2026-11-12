// BLASTN's -window_masker_taxid in the application (line NP WW2, DESIGN_ww2 §4): the form's
// taxonomy ID, the argv, the draft's rules (the ID typed, its counts downloaded), the run's name
// of the counts and its record, the browser's download (Range, restart, failures), the loader's
// line, the OPFS store, and the warm-up that leaves the option out.
import { describe, expect, it } from 'vitest';
import { CountsLoader } from '../../src/application/counts';
import { Coordinator, type EnqueueAllResult, type SearchRequest } from '../../src/application/coordinator';
import { SearchDraft } from '../../src/application/draft';
import { buildArgv } from '../../src/domain/argv';
import { COUNTS_HOST, countsName, countsUrl, parseTaxid, WINDOW_MASKER_TAXID_FLAG } from '../../src/domain/counts';
import { changedFields, formControlOf, formParameters, setFields, type EngineOption } from '../../src/domain/parameters';
import { programById, WINDOW_MASKER_FLAG } from '../../src/domain/programs';
import { BrowserCountsDownloader } from '../../src/infra/browser/counts-downloader';
import { sha256Hex } from '../../src/infra/browser/platform';
import { OpfsCountsStore } from '../../src/infra/data/counts-store';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { warmUp } from '../../src/infra/engine-worker/warm-up';
import { FakeCountsDeletionChannel, FakeCountsDownloader, MemoryCountsStore, type FakeCountsFile } from '../../src/infra/fake/fake-counts';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeAccessionFetcher } from '../../src/infra/fake/fake-fetcher';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import type { ReactorAbi } from '../../src/infra/reactor/abi';
import { COUNTS_BUSY, COUNTS_CHANGED, type CountsStore } from '../../src/ports/counts';
import type { DownloadProgress } from '../../src/ports/example';
import { MemoryOpfs } from './support/memory-opfs';
import { NodeInflater } from './support/node-inflater';

const HUMAN = 9606;
const LAST_MODIFIED = 'Thu, 01 Oct 2026 08:00:00 GMT';
const bytesOf = (length: number, seed = 1): Uint8Array<ArrayBuffer> => Uint8Array.from({ length }, (_, i) => (i * 31 + seed) & 255);
const COUNTS: FakeCountsFile = { bytes: bytesOf(5000), lastModified: LAST_MODIFIED };

// --- the form and the argv ---------------------------------------------------------------------

describe('the taxonomy ID', () => {
  it('is a decimal whole number of the range of a 32-bit integer, other than 0', () => {
    expect(parseTaxid('9606')).toBe(9606);
    expect(parseTaxid(' 562 ')).toBe(562);
    expect(parseTaxid('-1')).toBe(-1);
    expect(parseTaxid('2147483647')).toBe(2147483647);
    expect(parseTaxid('-2147483648')).toBe(-2147483648);
    for (const text of ['', '0', '-0', '2147483648', '-2147483649', '9606.0', '1e3', '0x10', '+5', 'abc', '12 34', '12345678901']) {
      expect(parseTaxid(text), text).toBeUndefined();
    }
  });

  it('names the counts as the engine looks for them, and the URL as NCBI keeps them', () => {
    expect(countsName(9606)).toBe('windowmasker/9606/wmasker.obinary');
    expect(countsName(-1)).toBe('windowmasker/-1/wmasker.obinary');
    expect(countsUrl(9606)).toBe('https://ftp.ncbi.nlm.nih.gov/blast/windowmasker_files/9606/wmasker.obinary');
    expect(COUNTS_HOST).toBe('ftp.ncbi.nlm.nih.gov');
  });
});

describe('the taxonomy ID field', () => {
  const blastn = programById('blastn');
  const options: EngineOption[] = [
    { flag: WINDOW_MASKER_FLAG, takesValue: true },
    { flag: WINDOW_MASKER_TAXID_FLAG, takesValue: true, defaultValue: '0' },
  ];

  it('is in the Filter row of the counts file, of the kind taxid, and only BLASTN has it', () => {
    const section = blastn.sections.find((s) => s.title === 'Filters and Masking')!;
    const field = section.fields.find((f) => f.flag === WINDOW_MASKER_TAXID_FLAG)!;
    expect(field).toMatchObject({ kind: 'taxid', row: 'Filter' });
    expect(section.fields.find((f) => f.flag === WINDOW_MASKER_FLAG)!.file?.taxid).toEqual({ flag: WINDOW_MASKER_TAXID_FLAG, missing: 'Enter a taxonomy ID.' });
    for (const id of ['blastp', 'blastx', 'tblastn', 'tblastx'] as const) {
      expect(programById(id).sections.some((s) => s.fields.some((f) => f.flag === WINDOW_MASKER_TAXID_FLAG))).toBe(false);
    }
  });

  it('writes -window_masker_taxid and the ID as typed, which marks the row as changed; checked without an ID it writes nothing', () => {
    const checked = setFields(blastn, {}, [[WINDOW_MASKER_TAXID_FLAG, true]], options);
    expect(formParameters(blastn, checked, options)).toEqual([]);
    expect(changedFields(blastn, checked, options, 'algorithm')).toBe(0);
    const typed = setFields(blastn, checked, [[WINDOW_MASKER_TAXID_FLAG, '9606']], options);
    expect(formParameters(blastn, typed, options)).toEqual([[WINDOW_MASKER_TAXID_FLAG, '9606']]);
    expect(changedFields(blastn, typed, options, 'algorithm')).toBe(1);
    expect(buildArgv({ program: 'blastn', queryName: 'q.fa', subjectName: 's.fa', parameters: formParameters(blastn, typed, options) })).toEqual([
      'blastn', '-query', 'q.fa', '-subject', 's.fa', '-window_masker_taxid', '9606',
    ]);
  });

  it('is refused in "Additional options" under the name of the row, as -window_masker_db is', () => {
    expect(formControlOf(blastn, options, WINDOW_MASKER_TAXID_FLAG)).toBe('Species-specific repeats filter for:');
    expect(formControlOf(blastn, options, WINDOW_MASKER_FLAG)).toBe('Species-specific repeats filter for:');
  });
});

// --- the draft and the run ----------------------------------------------------------------------

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

const counts = new File([COUNTS.bytes], 'wmasker.obinary', { lastModified: 1_760_000_000_000 });
const origin = { taxid: HUMAN, url: countsUrl(HUMAN), lastModified: LAST_MODIFIED };

describe('the draft and the counts of a taxonomy ID', () => {
  async function ready() {
    const s = draftSetup();
    s.draft.setPaste('query', '>q\nACGTACGT\n');
    s.draft.setPaste('subject', '>s\nACGTACGT\n');
    await s.draft.idle();
    return s;
  }

  it('does not start a search while the box is checked and the ID is not typed', async () => {
    const { draft, requests } = await ready();
    draft.setRepeats('taxid');
    await draft.idle();
    expect(draft.readiness()).toBe('Enter a taxonomy ID.');
    expect(await draft.submit()).toEqual({ ok: false, message: 'Enter a taxonomy ID.' });
    // An ID that is not a decimal whole number, or is 0 (no filter at NCBI), is not an ID.
    for (const text of ['0', '1.5', '-']) {
      draft.setTaxid(text);
      expect(draft.readiness(), text).toBe('Enter a taxonomy ID.');
    }
    expect(requests).toHaveLength(0);
    draft.setTaxid('');
    expect(draft.state.get().values.blastn[WINDOW_MASKER_TAXID_FLAG]).toBe(true);
    // Unchecked, it is a search without the filter.
    draft.setRepeats('off');
    expect(draft.readiness()).toBeUndefined();
    expect((await draft.submit()).ok).toBe(true);
    expect(requests[0]!.parameters).toEqual([]);
  });

  it('does not start a search until the counts of the ID are downloaded', async () => {
    const { draft, requests } = await ready();
    draft.setRepeats('taxid');
    draft.setTaxid('9606');
    await draft.idle();
    expect(draft.readiness()).toBe('Download the counts for 9606.');
    expect(await draft.submit()).toEqual({ ok: false, message: 'Download the counts for 9606.' });
    expect(requests).toHaveLength(0);
    // The counts of another ID do not serve.
    draft.setCounts(562, { file: counts, origin: { ...origin, taxid: 562 } });
    expect(draft.readiness()).toBe('Download the counts for 9606.');
    draft.setCounts(HUMAN, { file: counts, origin });
    expect(draft.readiness()).toBeUndefined();
    draft.setCounts(HUMAN, undefined);
    expect(draft.readiness()).toBe('Download the counts for 9606.');
  });

  it('searches with the kept file under the name the engine looks for, and the run records where it came from', async () => {
    const { draft, coordinator, requests } = await ready();
    draft.setRepeats('taxid');
    draft.setTaxid('9606');
    draft.setCounts(HUMAN, { file: counts, origin });
    await draft.idle();
    expect(draft.state.get().validation).toEqual({ state: 'ok' });
    expect((await draft.submit()).ok).toBe(true);
    expect(requests[0]!.parameters).toEqual([[WINDOW_MASKER_TAXID_FLAG, '9606']]);
    expect(requests[0]!.windowMasker).toBe(counts);
    const snapshot = coordinator.state.get().runs[0]!.snapshot;
    expect(snapshot.argv.slice(-2)).toEqual(['-window_masker_taxid', '9606']);
    expect(snapshot.windowMasker).toEqual({ name: 'windowmasker/9606/wmasker.obinary', file: counts, size: 5000, lastModified: 1_760_000_000_000, origin });
    for (let i = 0; i < 200 && coordinator.state.get().runs[0]?.status !== 'completed'; i++) await new Promise((r) => setTimeout(r, 5));
    const record = coordinator.state.get().runs[0]!.record;
    expect(record.windowMasker).toEqual({
      name: 'windowmasker/9606/wmasker.obinary',
      size: 5000,
      sha256: await sha256Hex(COUNTS.bytes),
      retained: false,
      loadMs: 0,
      origin,
    });
  });

  it('changes the row between the ID and a file: one of them is in the argv, and a chosen file goes', async () => {
    const { draft } = await ready();
    draft.setRepeats('taxid');
    draft.setTaxid('9606');
    draft.setRepeats('file');
    expect(draft.state.get().values.blastn[WINDOW_MASKER_TAXID_FLAG]).toBe('');
    expect(draft.state.get().values.blastn[WINDOW_MASKER_FLAG]).toBe(true);
    expect(draft.readiness()).toBe('Choose the WindowMasker counts file.');
    draft.setFile(WINDOW_MASKER_FLAG, counts);
    expect(draft.parameters()).toEqual([[WINDOW_MASKER_FLAG, 'wmasker.obinary']]);
    draft.setRepeats('taxid');
    expect(draft.state.get().files.blastn[WINDOW_MASKER_FLAG]).toBeUndefined();
    expect(draft.parameters()).toEqual([]);
    expect(draft.readiness()).toBe('Enter a taxonomy ID.');
  });

  it('lets go of the row with "Restore default search parameters", and refuses the option in "Additional options"', async () => {
    const { draft } = await ready();
    draft.setRepeats('taxid');
    draft.setTaxid('9606');
    draft.resetFields([WINDOW_MASKER_FLAG, WINDOW_MASKER_TAXID_FLAG]);
    expect(draft.state.get().values.blastn[WINDOW_MASKER_TAXID_FLAG]).toBeUndefined();
    expect(draft.readiness()).toBeUndefined();
    draft.setExtra('-window_masker_taxid 9606');
    expect(() => draft.extraWords()).toThrow('-window_masker_taxid is set in the form: Species-specific repeats filter for:');
  });
});

describe('the warm-up', () => {
  it('leaves -window_masker_taxid and its value out, as it leaves out -window_masker_db', () => {
    const validated: string[][] = [];
    const abi = {
      validate: (argv: readonly string[]) => {
        validated.push([...argv]);
        return 'stop here';
      },
      stoppedBy: undefined,
    } as unknown as ReactorAbi;
    warmUp(abi, ['blastn', '-query', 'q.fa', '-subject', 's.fa', '-window_masker_taxid', '9606', '-window_masker_db', 'x']);
    expect(validated[0]).toEqual(['blastn', '-query', 'q.fa', '-subject', 's.fa', '-num_threads', '1']);
  });
});

// --- the browser's download --------------------------------------------------------------------

/** A site with one file: it answers a Range as NCBI's FTP does, and can cut the body or change the file. */
function site(options: { file?: FakeCountsFile; status?: number; ignoreRange?: boolean; cutAt?: number } = {}) {
  const file = options.file ?? COUNTS;
  const calls: Array<{ url: string; headers: Record<string, string> }> = [];
  const state = { file, cutAt: options.cutAt, status: options.status };
  const fetchImpl: typeof fetch = async (input, init) => {
    const headers = (init?.headers ?? {}) as Record<string, string>;
    calls.push({ url: String(input), headers });
    const sent = state.file;
    if (state.status !== undefined && state.status !== 200) return new Response('no', { status: state.status });
    const range = /^bytes=(\d+)-$/.exec(headers['Range'] ?? '');
    const from = range === null || options.ignoreRange ? 0 : Number(range[1]);
    const rest = sent.bytes.slice(from);
    const cut = state.cutAt;
    state.cutAt = undefined;
    const end = cut === undefined ? rest.length : Math.min(cut, rest.length);
    let at = 0;
    // One chunk a pull: the reader gets the chunks before the error, as a connection that is cut gives them.
    const body = new ReadableStream<Uint8Array>({
      pull(controller) {
        if (at < end) {
          controller.enqueue(rest.slice(at, Math.min(end, at + 1000)));
          at += 1000;
        } else if (cut === undefined) controller.close();
        else controller.error(new TypeError('network error'));
      },
    });
    const answer: Record<string, string> = { 'Last-Modified': sent.lastModified, 'Content-Length': String(rest.length) };
    if (range !== null && !options.ignoreRange) answer['Content-Range'] = `bytes ${from}-${sent.bytes.length - 1}/${sent.bytes.length}`;
    return new Response(body, { status: range !== null && !options.ignoreRange ? 206 : 200, headers: answer });
  };
  return { fetchImpl, calls, state };
}

describe('BrowserCountsDownloader', () => {
  const request = (resume?: Parameters<BrowserCountsDownloader['download']>[0]['resume']) => ({ taxid: HUMAN, url: countsUrl(HUMAN), ...(resume === undefined ? {} : { resume }) });
  const make = (fetchImpl: typeof fetch, store: CountsStore = new MemoryCountsStore()) => ({
    store,
    downloader: new BrowserCountsDownloader({ store, fetch: fetchImpl }),
  });
  const run = async (downloader: BrowserCountsDownloader, req = request(), signal = new AbortController().signal) => {
    const progress: DownloadProgress[] = [];
    const result = await downloader.download(req, (p) => progress.push(p), signal);
    return { result, progress };
  };

  it('asks for the file by its URL only, without cookies or referrer, and keeps it whole', async () => {
    const s = site();
    const { downloader, store } = make(s.fetchImpl);
    const { result, progress } = await run(downloader);
    expect(result).toMatchObject({ ok: true, kept: true, copy: { taxid: HUMAN, size: 5000, savedBytes: 5000, complete: true, lastModified: LAST_MODIFIED } });
    expect(s.calls).toEqual([{ url: countsUrl(HUMAN), headers: {} }]);
    expect(progress[0]).toEqual({ bytes: 0, total: 5000 });
    expect(progress.at(-1)).toEqual({ bytes: 5000, total: 5000 });
    const saved = await store.readCounts(HUMAN, countsUrl(HUMAN));
    expect(new Uint8Array(await saved!.file!.arrayBuffer())).toEqual(COUNTS.bytes);
  });

  it('sends the request as a CORS request with no credentials and no referrer', async () => {
    const inits: Array<RequestInit | undefined> = [];
    const s = site();
    await run(make(async (input, init) => (inits.push(init), s.fetchImpl(input, init))).downloader);
    expect(inits[0]).toMatchObject({ mode: 'cors', credentials: 'omit', referrerPolicy: 'no-referrer' });
  });

  it('keeps what came when the connection is cut, and continues with a Range from the saved byte to the same file', async () => {
    const s = site({ cutAt: 2500 });
    const { downloader, store } = make(s.fetchImpl);
    const first = await run(downloader);
    expect(first.result).toMatchObject({ ok: false, reason: 'unreachable', detail: 'network error' });
    const partial = (await store.readCounts(HUMAN, countsUrl(HUMAN)))!;
    expect(partial.copy).toMatchObject({ complete: false, savedBytes: 2500, size: 5000 });
    expect(partial.file).toBeUndefined();
    const second = await run(downloader, request(partial.copy));
    expect(s.calls[1]!.headers).toEqual({ Range: 'bytes=2500-' });
    expect(second.progress[0]).toEqual({ bytes: 2500, total: 5000 });
    expect(second.result).toMatchObject({ ok: true, kept: true });
    const saved = await store.readCounts(HUMAN, countsUrl(HUMAN));
    expect(new Uint8Array(await saved!.file!.arrayBuffer())).toEqual(COUNTS.bytes);
  });

  it('keeps the bytes that came when the user cancels, and stops with the signal', async () => {
    const controller = new AbortController();
    const s = site();
    const { downloader, store } = make(s.fetchImpl);
    const progress: DownloadProgress[] = [];
    const pending = downloader.download(request(), (p) => {
      progress.push(p);
      if (p.bytes >= 2000) controller.abort(new DOMException('cancelled', 'AbortError'));
    }, controller.signal);
    await expect(pending).rejects.toMatchObject({ name: 'AbortError' });
    const partial = (await store.readCounts(HUMAN, countsUrl(HUMAN)))!;
    expect(partial.copy.complete).toBe(false);
    expect(partial.copy.savedBytes).toBeGreaterThanOrEqual(2000);
    expect(partial.copy.savedBytes).toBeLessThan(5000);
  });

  it('starts again from the first byte when the file changed on the site (another Last-Modified)', async () => {
    const s = site({ cutAt: 2500 });
    const { downloader, store } = make(s.fetchImpl);
    await run(downloader);
    const partial = (await store.readCounts(HUMAN, countsUrl(HUMAN)))!.copy;
    const changed: FakeCountsFile = { bytes: bytesOf(6000, 7), lastModified: 'Fri, 02 Oct 2026 08:00:00 GMT' };
    s.state.file = changed;
    const { result, progress } = await run(downloader, request(partial));
    expect(s.calls.map((c) => c.headers['Range'])).toEqual([undefined, 'bytes=2500-', undefined]);
    expect(progress[0]).toEqual({ bytes: 0, total: 6000 });
    expect(result).toMatchObject({ ok: true, copy: { size: 6000, lastModified: changed.lastModified } });
    const saved = await store.readCounts(HUMAN, countsUrl(HUMAN));
    expect(new Uint8Array(await saved!.file!.arrayBuffer())).toEqual(changed.bytes);
  });

  it('starts again from the first byte when the answer is not a 206 (the site ignored the Range) and reads that answer', async () => {
    const s = site({ cutAt: 2500, ignoreRange: true });
    const { downloader, store } = make(s.fetchImpl);
    await run(downloader);
    const partial = (await store.readCounts(HUMAN, countsUrl(HUMAN)))!.copy;
    const { result } = await run(downloader, request(partial));
    // One more request only: the whole file answered the Range.
    expect(s.calls).toHaveLength(2);
    expect(result).toMatchObject({ ok: true, copy: { savedBytes: 5000 } });
    const saved = await store.readCounts(HUMAN, countsUrl(HUMAN));
    expect(new Uint8Array(await saved!.file!.arrayBuffer())).toEqual(COUNTS.bytes);
  });

  it('starts again from the first byte when the Content-Range starts elsewhere', async () => {
    const s = site({ cutAt: 2500 });
    const { downloader, store } = make(async (input, init) => {
      const response = await s.fetchImpl(input, init);
      const headers = new Headers(response.headers);
      if (headers.has('Content-Range')) headers.set('Content-Range', 'bytes 100-4999/5000');
      return new Response(response.body, { status: response.status, headers });
    });
    await run(downloader);
    const partial = (await store.readCounts(HUMAN, countsUrl(HUMAN)))!.copy;
    const { result } = await run(downloader, request(partial));
    expect(result).toMatchObject({ ok: true, copy: { savedBytes: 5000 } });
    // The ranged answer was dropped for the file from its start.
    expect(s.calls.map((c) => c.headers['Range'])).toEqual([undefined, 'bytes=2500-', undefined]);
  });

  it('reports a status (404, 500) and an unreachable site, and keeps the saved part', async () => {
    for (const [status, expected] of [[404, 404], [500, 500]] as const) {
      const { downloader } = make(site({ status }).fetchImpl);
      expect((await run(downloader)).result).toEqual({ ok: false, reason: 'refused', status: expected });
    }
    const offline = make(async () => {
      throw new TypeError('Failed to fetch');
    });
    expect((await run(offline.downloader)).result).toEqual({ ok: false, reason: 'unreachable', detail: 'Failed to fetch' });
    // A refused continue leaves the saved part for the next try.
    const s = site({ cutAt: 2500 });
    const { downloader, store } = make(s.fetchImpl);
    await run(downloader);
    const partial = (await store.readCounts(HUMAN, countsUrl(HUMAN)))!.copy;
    s.state.status = 500;
    expect((await run(downloader, request(partial))).result).toMatchObject({ ok: false, reason: 'refused', status: 500 });
    expect((await store.readCounts(HUMAN, countsUrl(HUMAN)))!.copy.savedBytes).toBe(2500);
  });

  it('holds the file in this tab when the browser keeps nothing', async () => {
    const { downloader } = make(site().fetchImpl, new MemoryCountsStore({ persistent: false }));
    const { result } = await run(downloader);
    expect(result).toMatchObject({ ok: true, kept: false, copy: { complete: true, size: 5000 } });
    if (result.ok) expect(new Uint8Array(await result.file.arrayBuffer())).toEqual(COUNTS.bytes);
  });

  it('reports another tab that writes the file, and a storage that refuses a part', async () => {
    const store = new MemoryCountsStore();
    const first = await store.beginCounts(HUMAN, { url: countsUrl(HUMAN), lastModified: LAST_MODIFIED, size: 5000, from: 0 });
    expect(first).toBeDefined();
    const { downloader } = make(site().fetchImpl, store);
    expect((await run(downloader)).result).toEqual({ ok: false, reason: 'busy' });
    await store.closeCounts(first!);
    const full = new MemoryCountsStore();
    full.appendError = new Error('The browser refused to keep the counts: its storage is full.');
    expect((await run(make(site().fetchImpl, full).downloader)).result).toEqual({
      ok: false,
      reason: 'not-kept',
      detail: 'The browser refused to keep the counts: its storage is full.',
    });
  });

  it('writes the body in parts, never all at once', async () => {
    const big = { bytes: bytesOf(20 * 1024 * 1024, 3), lastModified: LAST_MODIFIED };
    const s = site({ file: big });
    const store = new MemoryCountsStore();
    const appended: number[] = [];
    const spy: CountsStore = {
      readCounts: (t, u) => store.readCounts(t, u),
      beginCounts: (t, st) => store.beginCounts(t, st),
      appendCounts: (token, bytes) => (appended.push(bytes.length), store.appendCounts(token, bytes)),
      finishCounts: (token) => store.finishCounts(token),
      closeCounts: (token) => store.closeCounts(token),
      deleteCounts: (t) => store.deleteCounts(t),
    };
    const { result } = await run(new BrowserCountsDownloader({ store: spy, fetch: s.fetchImpl }));
    expect(result.ok).toBe(true);
    expect(appended.reduce((a, b) => a + b, 0)).toBe(big.bytes.length);
    expect(Math.max(...appended)).toBeLessThan(9 * 1024 * 1024);
    expect(appended.length).toBeGreaterThanOrEqual(3);
  });
});

// --- the loader ---------------------------------------------------------------------------------

describe('CountsLoader', () => {
  function setup(options: { files?: Record<number, FakeCountsFile>; persistent?: boolean; channel?: FakeCountsDeletionChannel } = {}) {
    const store = new MemoryCountsStore({ persistent: options.persistent ?? true, now: () => Date.UTC(2026, 10, 12) });
    const downloader = new FakeCountsDownloader(store, options.files ?? { [HUMAN]: COUNTS });
    const saved: Array<[number, unknown]> = [];
    const draftState = new Map<number, unknown>();
    const draft = {
      setCounts: (taxid: number, entry: unknown) => {
        saved.push([taxid, entry]);
        if (entry === undefined) draftState.delete(taxid);
        else draftState.set(taxid, entry);
      },
    };
    const loader = new CountsLoader({ downloader, store, draft, ...(options.channel === undefined ? {} : { deletions: options.channel.tab() }) });
    return { loader, store, downloader, draftState };
  }

  it('sends nothing until Download is pressed: a look-up reads the storage only', async () => {
    const { loader, downloader } = setup();
    loader.refresh(HUMAN);
    await loader.idle();
    expect(loader.state.get().saved[HUMAN]).toBeNull();
    expect(downloader.requests).toEqual([]);
  });

  it('downloads, hands the whole file to the draft with where it came from, and finds it again without a request', async () => {
    const { loader, downloader, draftState, store } = setup();
    loader.download(HUMAN);
    await loader.idle();
    expect(loader.state.get().load).toBeUndefined();
    expect(loader.state.get().saved[HUMAN]).toMatchObject({ complete: true, size: 5000 });
    expect(draftState.get(HUMAN)).toMatchObject({ origin: { taxid: HUMAN, url: countsUrl(HUMAN), lastModified: LAST_MODIFIED } });
    // A new tab (a new loader over the same storage) finds it and sends nothing.
    const again = new CountsLoader({ downloader, store, draft: { setCounts: () => undefined } });
    again.refresh(HUMAN);
    await again.idle();
    expect(again.state.get().saved[HUMAN]).toMatchObject({ complete: true });
    expect(downloader.requests).toHaveLength(1);
    // Download of a file that is saved reads it, and asks again for nothing.
    again.download(HUMAN);
    await again.idle();
    expect(downloader.requests).toHaveLength(1);
  });

  it('says the site has no counts for an ID (404), and names a failure that may pass', async () => {
    const { loader, downloader } = setup();
    loader.download(12345);
    await loader.idle();
    expect(loader.state.get().load).toEqual({ taxid: 12345, progress: { state: 'failed', message: 'ftp.ncbi.nlm.nih.gov has no WindowMasker counts for 12345.' } });
    downloader.next = { ok: false, reason: 'refused', status: 500 };
    loader.download(HUMAN);
    await loader.idle();
    expect(loader.state.get().load?.progress).toEqual({ state: 'failed', message: 'ftp.ncbi.nlm.nih.gov answered HTTP 500: try again in a moment.' });
    downloader.next = { ok: false, reason: 'unreachable', detail: 'Failed to fetch' };
    loader.download(HUMAN);
    await loader.idle();
    expect(loader.state.get().load?.progress).toEqual({ state: 'failed', message: 'ftp.ncbi.nlm.nih.gov could not be reached (Failed to fetch).' });
    downloader.next = { ok: false, reason: 'not-kept', detail: 'The browser refused to keep the counts: its storage is full.' };
    loader.download(HUMAN);
    await loader.idle();
    expect(loader.state.get().load?.progress).toEqual({ state: 'failed', message: 'Not kept: The browser refused to keep the counts: its storage is full.' });
    downloader.next = { ok: false, reason: 'busy' };
    loader.download(HUMAN);
    await loader.idle();
    expect(loader.state.get().load?.progress).toEqual({ state: 'failed', message: 'Another tab is downloading the counts for 9606.' });
  });

  it('shows a download that stopped as saved bytes, continues it from the saved copy, and deletes it', async () => {
    const { loader, downloader, store, draftState } = setup();
    const token = (await store.beginCounts(HUMAN, { url: countsUrl(HUMAN), lastModified: LAST_MODIFIED, size: 5000, from: 0 }))!;
    await store.appendCounts(token, COUNTS.bytes.slice(0, 2000));
    await store.closeCounts(token);
    loader.refresh(HUMAN);
    await loader.idle();
    expect(loader.state.get().saved[HUMAN]).toMatchObject({ complete: false, savedBytes: 2000, size: 5000 });
    expect(draftState.has(HUMAN)).toBe(false);
    loader.download(HUMAN);
    await loader.idle();
    expect(downloader.requests[0]!.resume).toMatchObject({ savedBytes: 2000 });
    expect(loader.state.get().saved[HUMAN]).toMatchObject({ complete: true });
    expect(draftState.has(HUMAN)).toBe(true);
    await loader.delete(HUMAN);
    expect(loader.state.get().saved[HUMAN]).toBeNull();
    expect(draftState.has(HUMAN)).toBe(false);
    expect(store.entries.has(HUMAN)).toBe(false);
  });

  it('holds the file in this tab when the browser keeps nothing, and lets go of it on Delete', async () => {
    const { loader, draftState } = setup({ persistent: false });
    loader.download(HUMAN);
    await loader.idle();
    expect(loader.state.get().saved[HUMAN]).toMatchObject({ complete: true });
    expect(loader.state.get().inTab).toEqual({ [HUMAN]: true });
    expect(draftState.has(HUMAN)).toBe(true);
    await loader.delete(HUMAN);
    expect(loader.state.get().saved[HUMAN]).toBeNull();
    expect(loader.state.get().inTab).toEqual({});
    expect(draftState.has(HUMAN)).toBe(false);
  });

  it('tells the other tabs of a deletion, which let go of the file', async () => {
    const channel = new FakeCountsDeletionChannel();
    const one = setup({ channel });
    const store = one.store;
    const otherDraft = new Map<number, unknown>();
    const other = new CountsLoader({
      downloader: one.downloader,
      store,
      draft: { setCounts: (t, e) => (e === undefined ? otherDraft.delete(t) : otherDraft.set(t, e)) },
      deletions: channel.tab(),
    });
    one.loader.download(HUMAN);
    await one.loader.idle();
    other.refresh(HUMAN);
    await other.idle();
    expect(otherDraft.has(HUMAN)).toBe(true);
    await one.loader.delete(HUMAN);
    expect(otherDraft.has(HUMAN)).toBe(false);
    expect(other.state.get().saved[HUMAN]).toBeNull();
  });
});

// --- the OPFS store -----------------------------------------------------------------------------

describe('OpfsCountsStore', () => {
  async function make() {
    const opfs = new MemoryOpfs();
    const root = opfs.root as unknown as FileSystemDirectoryHandle;
    let tokens = 0;
    const store = new OpfsCountsStore(root, opfs, () => `t${++tokens}`);
    return { opfs, store };
  }
  const start = (from = 0) => ({ url: countsUrl(HUMAN), lastModified: LAST_MODIFIED, size: 5000, from });

  it('writes the file part by part with the saved bytes in meta.json at windowmasker/<taxid>/, and is whole at the end', async () => {
    const { opfs, store } = await make();
    const token = (await store.beginCounts(HUMAN, start()))!;
    await store.appendCounts(token, COUNTS.bytes.slice(0, 3000));
    const meta = () => JSON.parse(new TextDecoder().decode(opfs.read(['9606', 'meta.json'])));
    expect(meta()).toMatchObject({ taxid: HUMAN, url: countsUrl(HUMAN), lastModified: LAST_MODIFIED, size: 5000, savedBytes: 3000, complete: false });
    expect(await store.readCounts(HUMAN, countsUrl(HUMAN))).toMatchObject({ copy: { savedBytes: 3000, complete: false } });
    await store.appendCounts(token, COUNTS.bytes.slice(3000));
    const saved = await store.finishCounts(token);
    expect(saved.copy).toMatchObject({ savedBytes: 5000, complete: true });
    expect(new Uint8Array(await saved.file!.arrayBuffer())).toEqual(COUNTS.bytes);
    expect(opfs.read(['9606', 'wmasker.obinary'])).toEqual(COUNTS.bytes);
    expect(meta().complete).toBe(true);
    // Another URL is another copy.
    expect(await store.readCounts(HUMAN, 'https://example.org/x')).toBeUndefined();
  });

  it('continues a saved part at its byte, and refuses a part that is not the one asked for', async () => {
    const { store } = await make();
    const token = (await store.beginCounts(HUMAN, start()))!;
    await store.appendCounts(token, COUNTS.bytes.slice(0, 2000));
    await store.closeCounts(token);
    await expect(store.beginCounts(HUMAN, start(1000))).rejects.toMatchObject({ name: COUNTS_CHANGED });
    await expect(store.beginCounts(HUMAN, { ...start(2000), lastModified: 'Fri, 02 Oct 2026 08:00:00 GMT' })).rejects.toMatchObject({ name: COUNTS_CHANGED });
    const next = (await store.beginCounts(HUMAN, start(2000)))!;
    await store.appendCounts(next, COUNTS.bytes.slice(2000));
    const saved = await store.finishCounts(next);
    expect(new Uint8Array(await saved.file!.arrayBuffer())).toEqual(COUNTS.bytes);
  });

  it('begins the file again from 0: what was saved is no longer read as a whole file', async () => {
    const { store } = await make();
    const token = (await store.beginCounts(HUMAN, start()))!;
    await store.appendCounts(token, COUNTS.bytes);
    await store.finishCounts(token);
    const again = (await store.beginCounts(HUMAN, { ...start(), size: 6000 }))!;
    expect(await store.readCounts(HUMAN, countsUrl(HUMAN))).toMatchObject({ copy: { savedBytes: 0, size: 6000, complete: false } });
    await store.closeCounts(again);
  });

  it('writes a file once at a time, and a file that is not whole does not finish', async () => {
    const { store } = await make();
    const token = (await store.beginCounts(HUMAN, start()))!;
    await expect(store.beginCounts(HUMAN, start())).rejects.toMatchObject({ name: COUNTS_BUSY });
    await store.appendCounts(token, COUNTS.bytes.slice(0, 10));
    await expect(store.finishCounts(token)).rejects.toThrow('not whole');
    await store.closeCounts(token);
    // The lock is let go with the writing.
    const next = (await store.beginCounts(HUMAN, start(10)))!;
    await store.closeCounts(next);
  });

  it('deletes the copy, also one that is not whole, and the storage of another ID stays', async () => {
    const { opfs, store } = await make();
    const a = (await store.beginCounts(HUMAN, start()))!;
    await store.appendCounts(a, COUNTS.bytes);
    await store.finishCounts(a);
    const b = (await store.beginCounts(562, { ...start(), url: countsUrl(562) }))!;
    await store.appendCounts(b, COUNTS.bytes.slice(0, 5));
    await store.closeCounts(b);
    await store.deleteCounts(HUMAN);
    expect(await store.readCounts(HUMAN, countsUrl(HUMAN))).toBeUndefined();
    expect(opfs.has(['9606'])).toBe(false);
    await store.deleteCounts(562);
    expect(opfs.has(['562'])).toBe(false);
    // Deleting what is not there is not an error.
    await store.deleteCounts(HUMAN);
  });

  it('reads a copy as absent when the file is shorter than the meta file says', async () => {
    const { opfs, store } = await make();
    const token = (await store.beginCounts(HUMAN, start()))!;
    await store.appendCounts(token, COUNTS.bytes.slice(0, 100));
    await store.closeCounts(token);
    opfs.write(['9606', 'wmasker.obinary'], new Uint8Array(10));
    expect(await store.readCounts(HUMAN, countsUrl(HUMAN))).toBeUndefined();
  });
});
