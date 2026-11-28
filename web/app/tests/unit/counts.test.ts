// BLASTN's -window_masker_taxid in the application (line NP WW2, DESIGN_ww2 §4): the form's
// taxonomy ID, the argv, the draft's rules (the ID typed, its counts downloaded), the run's name
// of the counts and its record, the browser's download (Range, restart, failures), the loader's
// line, the OPFS store, and the warm-up that leaves the option out.
import { describe, expect, it } from 'vitest';
import { sha256Hex } from 'losat';
import type { ReactorAbi } from 'losat';
import { WINDOW_MASKER_FLAG, WINDOW_MASKER_TAXID_FLAG } from 'losat/domain';
import { activeSavedCountsTaxid, CountsLoader, savedTabCopyOf } from '../../src/application/counts';
import { Coordinator, isToolRun, type EnqueueAllResult, type SearchRequest } from '../../src/application/coordinator';
import { SearchDraft } from '../../src/application/draft';
import { buildArgv } from '../../src/domain/argv';
import { COUNTS_FILE, COUNTS_HOST, countsEntryOf, countsName, countsUrl, parseTaxid, type CountsMadeMeta, type CountsOrigin } from '../../src/domain/counts';
import { changedFields, formControlOf, formParameters, setFields, type EngineOption } from '../../src/domain/parameters';
import { programById } from '../../src/domain/programs';
import { BrowserCountsDownloader } from '../../src/infra/browser/counts-downloader';
import { webCountsLocks } from '../../src/infra/browser/counts-locks';
import { OpfsCountsStore } from '../../src/infra/data/counts-store';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { warmUp } from '../../../packages/losat/src/infra/engine-worker/warm-up';
import { FakeCountsChangeChannel, FakeCountsDownloader, MemoryCountsStore, type FakeCountsFile } from '../../src/infra/fake/fake-counts';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeAccessionFetcher } from '../../src/infra/fake/fake-fetcher';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import { COUNTS_BUSY, COUNTS_CHANGED, COUNTS_EXISTS, COUNTS_IN_USE, type CountsCopy, type CountsStore, type SavedCountsEntry } from '../../src/ports/counts';
import type { DownloadProgress } from '../../src/ports/example';
import { MemoryOpfs } from './support/memory-opfs';
import { NodeInflater } from './support/node-inflater';
import { runEnd } from './support/run-end';

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

  it('builds one SavedCountsEntry from any complete CountsCopy (madeAt, Last-Modified, or savedAt fallback)', () => {
    const made: CountsCopy = {
      taxid: 562,
      url: '',
      lastModified: '',
      size: 700,
      savedBytes: 700,
      complete: true,
      savedAt: Date.UTC(2026, 10, 25),
      made: madeMeta(4),
    };
    expect(countsEntryOf(made)).toEqual({
      taxid: 562,
      name: COUNTS_FILE,
      size: 700,
      origin: 'made',
      date: '2026-11-22',
      copy: made,
    });
    const downloaded: CountsCopy = {
      taxid: HUMAN,
      url: countsUrl(HUMAN),
      lastModified: LAST_MODIFIED,
      size: 5000,
      savedBytes: 5000,
      complete: true,
      savedAt: Date.UTC(2026, 10, 25),
    };
    expect(countsEntryOf(downloaded)).toEqual({
      taxid: HUMAN,
      name: COUNTS_FILE,
      size: 5000,
      origin: 'downloaded',
      date: '2026-10-01',
      copy: downloaded,
    });
    const noDate: CountsCopy = { ...downloaded, lastModified: '' };
    expect(countsEntryOf(noDate).date).toBe('2026-11-25');
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
  const engine = { describe: fake.describe.bind(fake), start: fake.start.bind(fake), validate: fake.validate.bind(fake), prepare: () => undefined };
  let id = 0;
  const coordinator = new Coordinator({
    engine: { ...engine, run: fake.run.bind(fake), cancel: fake.cancel.bind(fake) },
    data,
    downloader: { save: () => undefined },
    now: () => 0,
    newRunId: () => `run-${++id}`,
    digest: sha256Hex,
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
const origin: CountsOrigin = { kind: 'downloaded', taxid: HUMAN, url: countsUrl(HUMAN), lastModified: LAST_MODIFIED };
const MADE_AT = Date.UTC(2026, 10, 22, 9, 0, 0);
const madeMeta = (run = 4): CountsMadeMeta => ({
  run,
  argv: ['windowmasker', '-mk_counts', '-in', 'genome.fa', '-out', 'wmasker.obinary', '-sformat', 'obinary'],
  inputs: [{ name: 'genome.fa', size: 12345 }],
  engineBuild: 'build-1',
  madeAt: MADE_AT,
});

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
    const run0 = coordinator.state.get().runs[0]!;
    if (isToolRun(run0)) throw new Error('expected search run');
    const snapshot = run0.snapshot;
    expect(snapshot.argv.slice(-2)).toEqual(['-window_masker_taxid', '9606']);
    expect(snapshot.windowMasker).toEqual({ name: 'windowmasker/9606/wmasker.obinary', file: counts, size: 5000, lastModified: 1_760_000_000_000, origin });
    await runEnd(coordinator, coordinator.state.get().runs[0]!.snapshot.runId);
    const finished = coordinator.state.get().runs[0]!;
    if (isToolRun(finished)) throw new Error('expected search run');
    const record = finished.record;
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
    const saved = await store.readCounts(HUMAN);
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
    const partial = (await store.readCounts(HUMAN))!;
    expect(partial.copy).toMatchObject({ complete: false, savedBytes: 2500, size: 5000 });
    expect(partial.file).toBeUndefined();
    const second = await run(downloader, request(partial.copy));
    expect(s.calls[1]!.headers).toEqual({ Range: 'bytes=2500-' });
    expect(second.progress[0]).toEqual({ bytes: 2500, total: 5000 });
    expect(second.result).toMatchObject({ ok: true, kept: true });
    const saved = await store.readCounts(HUMAN);
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
    const partial = (await store.readCounts(HUMAN))!;
    expect(partial.copy.complete).toBe(false);
    expect(partial.copy.savedBytes).toBeGreaterThanOrEqual(2000);
    expect(partial.copy.savedBytes).toBeLessThan(5000);
  });

  it('starts again from the first byte when the file changed on the site (another Last-Modified)', async () => {
    const s = site({ cutAt: 2500 });
    const { downloader, store } = make(s.fetchImpl);
    await run(downloader);
    const partial = (await store.readCounts(HUMAN))!.copy;
    const changed: FakeCountsFile = { bytes: bytesOf(6000, 7), lastModified: 'Fri, 02 Oct 2026 08:00:00 GMT' };
    s.state.file = changed;
    const { result, progress } = await run(downloader, request(partial));
    expect(s.calls.map((c) => c.headers['Range'])).toEqual([undefined, 'bytes=2500-', undefined]);
    expect(progress[0]).toEqual({ bytes: 0, total: 6000 });
    expect(result).toMatchObject({ ok: true, copy: { size: 6000, lastModified: changed.lastModified } });
    const saved = await store.readCounts(HUMAN);
    expect(new Uint8Array(await saved!.file!.arrayBuffer())).toEqual(changed.bytes);
  });

  it('starts again from the first byte when the answer is not a 206 (the site ignored the Range) and reads that answer', async () => {
    const s = site({ cutAt: 2500, ignoreRange: true });
    const { downloader, store } = make(s.fetchImpl);
    await run(downloader);
    const partial = (await store.readCounts(HUMAN))!.copy;
    const { result } = await run(downloader, request(partial));
    // One more request only: the whole file answered the Range.
    expect(s.calls).toHaveLength(2);
    expect(result).toMatchObject({ ok: true, copy: { savedBytes: 5000 } });
    const saved = await store.readCounts(HUMAN);
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
    const partial = (await store.readCounts(HUMAN))!.copy;
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
    const partial = (await store.readCounts(HUMAN))!.copy;
    s.state.status = 500;
    const asked = s.calls.length;
    expect((await run(downloader, request(partial))).result).toMatchObject({ ok: false, reason: 'refused', status: 500 });
    // One request: a refused continue is not followed by a download from the first byte.
    expect(s.calls.length - asked).toBe(1);
    expect((await store.readCounts(HUMAN))!.copy.savedBytes).toBe(2500);
  });

  it('stops when Cancel comes while the store opens the copy, and lets go of it', async () => {
    const controller = new AbortController();
    class CancelledAtBegin extends MemoryCountsStore {
      override async beginCounts(...args: Parameters<MemoryCountsStore['beginCounts']>) {
        const token = await super.beginCounts(...args);
        controller.abort(new DOMException('cancelled', 'AbortError'));
        return token;
      }
    }
    const store = new CancelledAtBegin();
    const s = site();
    const { downloader } = make(s.fetchImpl, store);
    await expect(downloader.download(request(), () => undefined, controller.signal)).rejects.toMatchObject({ name: 'AbortError' });
    expect((await store.readCounts(HUMAN))!.copy).toMatchObject({ complete: false, savedBytes: 0 });
    // The copy is not held: the next download writes it.
    expect((await run(make(site().fetchImpl, store).downloader)).result).toMatchObject({ ok: true, copy: { complete: true } });
  });

  it('finishes a copy whose bytes all came but that was not finished by asking for its last byte again', async () => {
    const s = site();
    const { downloader, store } = make(s.fetchImpl);
    const token = (await store.beginCounts(HUMAN, { url: countsUrl(HUMAN), lastModified: LAST_MODIFIED, size: 5000, from: 0 }))!;
    await store.appendCounts(token, COUNTS.bytes.slice());
    await store.closeCounts(token);
    const whole = (await store.readCounts(HUMAN))!.copy;
    expect(whole).toMatchObject({ complete: false, savedBytes: 5000 });
    const { result } = await run(downloader, request(whole));
    expect(s.calls.map((c) => c.headers['Range'])).toEqual(['bytes=4999-']);
    expect(result).toMatchObject({ ok: true, copy: { complete: true, size: 5000 } });
    const saved = await store.readCounts(HUMAN);
    expect(new Uint8Array(await saved!.file!.arrayBuffer())).toEqual(COUNTS.bytes);
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
      readCounts: (t) => store.readCounts(t),
      listCounts: () => store.listCounts(),
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
  function setup(options: { files?: Record<number, FakeCountsFile>; persistent?: boolean; channel?: FakeCountsChangeChannel } = {}) {
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
    const loader = new CountsLoader({ downloader, store, draft, ...(options.channel === undefined ? {} : { changes: options.channel.tab() }) });
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
    expect(draftState.get(HUMAN)).toMatchObject({ origin: { kind: 'downloaded', taxid: HUMAN, url: countsUrl(HUMAN), lastModified: LAST_MODIFIED } });
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

  it('still counts a downloaded copy only if its URL matches countsUrl(taxid)', async () => {
    const { loader, store, draftState } = setup();
    const token = (await store.beginCounts(HUMAN, { url: 'https://example.org/other.obinary', lastModified: LAST_MODIFIED, size: 5000, from: 0 }))!;
    await store.appendCounts(token, COUNTS.bytes.slice());
    await store.finishCounts(token);
    loader.refresh(HUMAN);
    loader.refreshList();
    await loader.idle();
    expect(loader.state.get().saved[HUMAN]).toBeNull();
    expect(loader.state.get().list).toEqual([]);
    expect(draftState.has(HUMAN)).toBe(false);
  });

  it('saves a made copy in 8 MiB parts, hands it to setCounts, and lists it without sending any request', async () => {
    const { loader, downloader, draftState, store } = setup();
    const appended: number[] = [];
    const origAppend = store.appendCounts.bind(store);
    store.appendCounts = (token, bytes) => {
      appended.push(bytes.length);
      return origAppend(token, bytes);
    };
    const madeBytes = bytesOf(18 * 1024 * 1024 + 5, 9);
    const madeFile = new File([madeBytes], 'custom.obinary');
    const meta = madeMeta(4);
    const result = await loader.saveMade(562, madeFile, meta);
    expect(result).toMatchObject({ ok: true, kept: true, copy: { taxid: 562, url: '', lastModified: '', size: madeBytes.length, complete: true, made: meta } });
    expect(appended).toEqual([8 * 1024 * 1024, 8 * 1024 * 1024, 2 * 1024 * 1024 + 5]);
    expect(downloader.requests).toEqual([]);
    expect(draftState.get(562)).toMatchObject({ origin: { kind: 'made', taxid: 562, made: meta } });
    expect(loader.state.get().list).toMatchObject([{ taxid: 562, name: COUNTS_FILE, size: madeBytes.length, origin: 'made', date: '2026-11-22' }]);
  });

  it('asks for Replace over a download, a partial download or a made copy, and overwrites when Replace is confirmed', async () => {
    const { loader, store, draftState } = setup();
    // 1. Replace a complete download.
    loader.download(HUMAN);
    await loader.idle();
    const firstMade = new File([bytesOf(1234, 2)], 'first.obinary');
    const ask1 = await loader.saveMade(HUMAN, firstMade, madeMeta(2));
    expect(ask1).toMatchObject({ ok: false, reason: 'exists', existing: { taxid: HUMAN, url: countsUrl(HUMAN), complete: true } });
    // Existing download was untouched.
    expect((await store.readCounts(HUMAN))!.origin).toBe('downloaded');
    const rep1 = await loader.saveMade(HUMAN, firstMade, madeMeta(2), { replace: true });
    expect(rep1).toMatchObject({ ok: true, kept: true, copy: { taxid: HUMAN, made: madeMeta(2) } });
    expect(draftState.get(HUMAN)).toMatchObject({ origin: { kind: 'made', taxid: HUMAN, made: madeMeta(2) } });

    // 2. Replace a made copy with another made copy.
    const secondMade = new File([bytesOf(2345, 5)], 'second.obinary');
    const ask2 = await loader.saveMade(HUMAN, secondMade, madeMeta(5));
    expect(ask2).toMatchObject({ ok: false, reason: 'exists', existing: { taxid: HUMAN, made: madeMeta(2) } });
    const rep2 = await loader.saveMade(HUMAN, secondMade, madeMeta(5), { replace: true });
    expect(rep2).toMatchObject({ ok: true, copy: { size: 2345, made: madeMeta(5) } });

    // 3. Replace a partial download.
    await loader.delete(HUMAN);
    const partialTok = (await store.beginCounts(HUMAN, { url: countsUrl(HUMAN), lastModified: LAST_MODIFIED, size: 5000, from: 0 }))!;
    await store.appendCounts(partialTok, COUNTS.bytes.slice(0, 1000));
    await store.closeCounts(partialTok);
    const ask3 = await loader.saveMade(HUMAN, secondMade, madeMeta(6));
    expect(ask3).toMatchObject({ ok: false, reason: 'exists', existing: { complete: false, savedBytes: 1000, url: countsUrl(HUMAN) } });
    const rep3 = await loader.saveMade(HUMAN, secondMade, madeMeta(6), { replace: true });
    expect(rep3).toMatchObject({ ok: true, copy: { complete: true, size: 2345, made: madeMeta(6) } });
  });

  it('refuses a save in place while another tab is saving or downloading the same ID', async () => {
    const { loader, store } = setup();
    const lock = (await store.beginCounts(HUMAN, { url: countsUrl(HUMAN), lastModified: LAST_MODIFIED, size: 5000, from: 0 }))!;
    const res = await loader.saveMade(HUMAN, new File([bytesOf(50)], 'x.obinary'), madeMeta(3));
    expect(res).toEqual({
      ok: false,
      reason: 'busy',
      message: 'Not saved: another tab is saving or downloading the counts for 9606. Save again when it has ended.',
    });
    await store.closeCounts(lock);
  });

  it('holds a made copy in this tab when the browser keeps nothing, and supports Replace and Delete', async () => {
    const { loader, draftState } = setup({ persistent: false });
    const first = new File([bytesOf(400, 1)], 'a.obinary');
    const saved = await loader.saveMade(562, first, madeMeta(1));
    expect(saved).toMatchObject({ ok: true, kept: false, copy: { taxid: 562, size: 400, made: madeMeta(1) } });
    expect(loader.state.get().inTab).toEqual({ 562: true });
    expect(loader.state.get().list).toMatchObject([{ taxid: 562, origin: 'made', size: 400 }]);
    expect(draftState.has(562)).toBe(true);
    expect(await loader.saveMade(562, new File([bytesOf(500, 2)], 'b.obinary'), madeMeta(2))).toMatchObject({ ok: false, reason: 'exists' });
    expect(await loader.saveMade(562, new File([bytesOf(500, 2)], 'b.obinary'), madeMeta(2), { replace: true })).toMatchObject({ ok: true, kept: false, copy: { size: 500 } });
    await loader.delete(562);
    expect(loader.state.get().inTab).toEqual({});
    expect(loader.state.get().list).toEqual([]);
    expect(draftState.has(562)).toBe(false);
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
    expect(loader.state.get().load?.progress).toEqual({
      state: 'failed',
      message: 'Not downloaded: another tab is saving or downloading the counts for 9606. Download again when it has ended.',
    });
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

  it('tells the other tabs on save, Replace and Delete, so they read that ID again', async () => {
    const channel = new FakeCountsChangeChannel();
    const one = setup({ channel });
    const store = one.store;
    const otherDraft = new Map<number, unknown>();
    const other = new CountsLoader({
      downloader: one.downloader,
      store,
      draft: { setCounts: (t, e) => (e === undefined ? otherDraft.delete(t) : otherDraft.set(t, e)) },
      changes: channel.tab(),
    });
    // Save a made copy in tab 1 -> tab 2 reads it and hands it to its draft and list.
    await one.loader.saveMade(562, new File([bytesOf(300, 1)], 'a.obinary'), madeMeta(1));
    await other.idle();
    expect(other.state.get().saved[562]).toMatchObject({ complete: true, size: 300, made: madeMeta(1) });
    expect(otherDraft.get(562)).toMatchObject({ origin: { kind: 'made', taxid: 562, made: madeMeta(1) } });
    expect(other.state.get().list).toMatchObject([{ taxid: 562, size: 300, origin: 'made' }]);

    // Replace in tab 1 -> tab 2 reads the replacement.
    await one.loader.saveMade(562, new File([bytesOf(450, 2)], 'b.obinary'), madeMeta(2), { replace: true });
    await other.idle();
    expect(other.state.get().saved[562]).toMatchObject({ complete: true, size: 450, made: madeMeta(2) });
    expect(otherDraft.get(562)).toMatchObject({ origin: { kind: 'made', taxid: 562, made: madeMeta(2) } });
    expect(other.state.get().list).toMatchObject([{ taxid: 562, size: 450, origin: 'made' }]);

    // Delete in tab 1 -> tab 2 drops it from draft, saved and list.
    await one.loader.delete(562);
    expect(otherDraft.has(562)).toBe(false);
    await other.idle();
    expect(otherDraft.has(562)).toBe(false);
    expect(other.state.get().saved[562]).toBeNull();
    expect(other.state.get().list).toEqual([]);
  });

  function setupTwoTabs(wrapOtherStore?: (store: MemoryCountsStore) => CountsStore) {
    const channel = new FakeCountsChangeChannel();
    const one = setup({ channel });
    const otherDraft = new Map<number, { readonly file: File; readonly origin: CountsOrigin }>();
    const other = new CountsLoader({
      downloader: one.downloader,
      store: wrapOtherStore === undefined ? one.store : wrapOtherStore(one.store),
      draft: { setCounts: (t, e) => (e === undefined ? otherDraft.delete(t) : otherDraft.set(t, e)) },
      changes: channel.tab(),
    });
    let onStart: (() => void) | undefined;
    let gate = Promise.resolve();
    const origAppend = one.store.appendCounts.bind(one.store);
    one.store.appendCounts = async (token, bytes) => {
      if (onStart !== undefined) {
        const notify = onStart;
        onStart = undefined;
        notify();
        await gate;
      }
      return origAppend(token, bytes);
    };
    const pauseNextAppend = () => {
      let release!: () => void;
      gate = new Promise<void>((resolve) => (release = resolve));
      const started = new Promise<void>((resolve) => (onStart = resolve));
      return { started, release };
    };
    return { channel, one, store: one.store, other, otherDraft, pauseNextAppend };
  }

  it('tells the other tabs when a Replace fails after it started, so they drop the copy of before', async () => {
    const tabs = setupTwoTabs();
    expect(await tabs.one.loader.saveMade(562, new File([bytesOf(300, 1)], 'a.obinary'), madeMeta(1))).toMatchObject({ ok: true });
    await tabs.other.idle();
    expect(tabs.otherDraft.get(562)?.file).toBeDefined();
    expect(tabs.other.state.get().saved[562]).toMatchObject({ complete: true, size: 300, made: madeMeta(1) });

    const full = 'The browser refused to keep the counts: its storage is full.';
    tabs.store.appendError = new Error(full);
    const replaced = await tabs.one.loader.saveMade(562, new File([bytesOf(400, 2)], 'b.obinary'), madeMeta(2), { replace: true });
    await tabs.one.loader.idle();
    await tabs.other.idle();

    expect(replaced).toEqual({
      ok: false,
      reason: 'failed',
      message: `${full} The counts saved before for 562 are no longer kept.`,
    });
    expect(tabs.otherDraft.get(562)?.file).toBeUndefined();
    expect(tabs.otherDraft.has(562)).toBe(false);
    expect(tabs.other.state.get().saved[562]).toBeNull();
    expect(tabs.other.state.get().list).toEqual([]);
  });

  it('tells the other tabs as soon as a Replace begins so they drop the copy of before while it writes, and read the new copy when it ends', async () => {
    const { one, other, otherDraft, pauseNextAppend } = setupTwoTabs();
    expect(await one.loader.saveMade(562, new File([bytesOf(320, 1)], 'a.obinary'), madeMeta(1))).toMatchObject({ ok: true });
    await other.idle();
    const before = otherDraft.get(562)?.file;
    expect(before).toBeInstanceOf(File);
    expect(other.state.get().saved[562]).toMatchObject({ complete: true, size: 320, made: madeMeta(1) });

    // Pause appendCounts mid-write (same observable state for tab 2 as if tab 1 closed mid-write).
    const paused = pauseNextAppend();
    const replacing = one.loader.saveMade(562, new File([bytesOf(450, 2)], 'b.obinary'), madeMeta(2), { replace: true });
    await paused.started;
    await other.idle();

    // After begin (while tab 1 is still writing), tab 2 already dropped the file of before and shows no saved copy.
    expect(otherDraft.size).toBe(0);
    expect(other.state.get().saved[562]).toBe(null);
    expect(other.state.get().list).toHaveLength(0);

    // Let the write finish: tab 2 reads the new copy.
    paused.release();
    expect(await replacing).toMatchObject({ ok: true, kept: true, copy: { taxid: 562, size: 450, made: madeMeta(2) } });
    await other.idle();
    const after = otherDraft.get(562)?.file;
    expect(after).toBeInstanceOf(File);
    expect(after).not.toBe(before);
    expect(otherDraft.get(562)).toMatchObject({ origin: { kind: 'made', taxid: 562, made: madeMeta(2) } });
    expect(other.state.get().saved[562]).toMatchObject({ complete: true, size: 450, made: madeMeta(2) });
    expect(other.state.get().list).toMatchObject([{ taxid: 562, size: 450, origin: 'made' }]);
  });

  it('tracks saving state during a save (cleared on success and failure), refuses same-tab download without calling the downloader, and tells another tab it is saving or downloading', async () => {
    const tick = () => new Promise<void>((resolve) => setTimeout(resolve, 0));
    const expectFailedLoad = (loader: CountsLoader, taxid: number, message: string) =>
      expect(loader.state.get().load).toEqual({ taxid, progress: { state: 'failed', message } });
    const { one, store, pauseNextAppend } = setupTwoTabs();
    const otherDownloader = new BrowserCountsDownloader({ store, fetch: site().fetchImpl });
    const other = new CountsLoader({ downloader: otherDownloader, store, draft: { setCounts: () => undefined } });

    expect(await one.loader.saveMade(562, new File([bytesOf(300, 1)], 'a.obinary'), madeMeta(1))).toMatchObject({ ok: true });
    const paused = pauseNextAppend();
    const replacing = one.loader.saveMade(562, new File([bytesOf(400, 2)], 'b.obinary'), madeMeta(2), { replace: true });
    await paused.started;

    // The writing tab reports the in-progress save in state.saving.
    const saving562 = () => one.loader.state.get().saving[562];
    expect(saving562()).toEqual({ run: 2 });

    // Same-tab Download while this tab saves is refused before calling the downloader (0 requests).
    one.downloader.next = { ok: false, reason: 'busy' };
    one.loader.download(562);
    await tick();
    expect(one.downloader.requests).toEqual([]);
    expectFailedLoad(
      one.loader,
      562,
      'Not downloaded: this tab is saving the counts for 562. Download again when it has ended.',
    );
    expect(saving562()?.run).toBe(2);

    // Same-tab concurrent Save while this tab saves says this tab is saving.
    expect(await one.loader.saveMade(562, new File([bytesOf(100, 3)], 'c.obinary'), madeMeta(3), { replace: true })).toEqual({
      ok: false,
      reason: 'busy',
      message: 'Not saved: this tab is saving the counts for 562. Save again when it has ended.',
    });

    // Another tab trying Download or Save while tab 1 writes gets the neutral "another tab is saving or downloading" wording.
    other.download(562);
    await other.idle();
    expectFailedLoad(
      other,
      562,
      'Not downloaded: another tab is saving or downloading the counts for 562. Download again when it has ended.',
    );
    expect(await other.saveMade(562, new File([bytesOf(110, 3)], 'other-c.obinary'), madeMeta(3))).toEqual({
      ok: false,
      reason: 'busy',
      message: 'Not saved: another tab is saving or downloading the counts for 562. Save again when it has ended.',
    });

    // Success clears saving state and the stale same-tab refusal.
    paused.release();
    expect(await replacing).toMatchObject({ ok: true, copy: { size: 400, made: madeMeta(2) } });
    await one.loader.idle();
    expect(saving562()).toBeUndefined();
    expect(one.loader.state.get().load).toBeUndefined();

    // Failure also clears saving state.
    const pausedFail = pauseNextAppend();
    const failing = one.loader.saveMade(562, new File([bytesOf(500, 4)], 'd.obinary'), madeMeta(3), { replace: true });
    await pausedFail.started;
    expect(saving562()).toEqual({ run: 3 });
    store.appendError = new Error('storage full');
    pausedFail.release();
    expect(await failing).toEqual({
      ok: false,
      reason: 'failed',
      message: 'storage full The counts saved before for 562 are no longer kept.',
    });
    await one.loader.idle();
    expect(Object.keys(one.loader.state.get().saving)).toEqual([]);
  });

  it('refuses a same-tab save and second download in place while this tab is downloading the same ID', async () => {
    const { one, pauseNextAppend } = setupTwoTabs();
    const paused = pauseNextAppend();
    one.loader.download(HUMAN);
    await paused.started;

    const madeCopy = new File([bytesOf(300, 1)], 'same-tab.obinary');
    expect(await one.loader.saveMade(HUMAN, madeCopy, madeMeta(1))).toEqual({
      ok: false,
      reason: 'busy',
      message: 'Not saved: this tab is downloading the counts for 9606. Save again when it has ended.',
    });

    const requestsBefore = one.downloader.requests.length;
    one.loader.download(HUMAN);
    await new Promise<void>((resolve) => setTimeout(resolve, 0));
    expect(one.downloader.requests).toHaveLength(requestsBefore);
    expect(one.loader.state.get().load).toEqual({
      taxid: HUMAN,
      progress: {
        state: 'failed',
        message: 'Not downloaded: this tab is downloading the counts for 9606. Download again when it has ended.',
      },
    });

    paused.release();
    await one.loader.idle();
    expect(one.loader.state.get().saved[HUMAN]).toMatchObject({ complete: true });
    expect(await one.loader.saveMade(HUMAN, madeCopy, madeMeta(1))).toMatchObject({
      ok: false,
      reason: 'exists',
    });
  });

  function makeGate() {
    let onStart: (() => void) | undefined;
    let release!: () => void;
    let gate = Promise.resolve();
    const arm = () => {
      gate = new Promise<void>((resolve) => (release = resolve));
      const started = new Promise<void>((resolve) => (onStart = resolve));
      return { started, release };
    };
    const enter = async () => {
      if (onStart === undefined) return;
      const notify = onStart;
      onStart = undefined;
      notify();
      await gate;
    };
    return { arm, enter };
  }

  function gateStoreAppend(store: MemoryCountsStore) {
    const appendGate = makeGate();
    const origAppend = store.appendCounts.bind(store);
    store.appendCounts = async (token, bytes) => {
      await appendGate.enter();
      return origAppend(token, bytes);
    };
    return appendGate.arm;
  }

  function setupStoppedDownload() {
    const store = new MemoryCountsStore({ persistent: true, now: () => Date.UTC(2026, 10, 12) });
    const s = site();
    const loader = new CountsLoader({
      downloader: new BrowserCountsDownloader({ store, fetch: s.fetchImpl }),
      store,
      draft: { setCounts: () => undefined },
    });
    const pauseNextAppend = gateStoreAppend(store);
    const startAndStop = async () => {
      const paused = pauseNextAppend();
      loader.download(HUMAN);
      await paused.started;
      loader.cancel();
      return paused.release;
    };
    const expectReading = () => expect(loader.state.get().load).toEqual({ taxid: HUMAN, progress: { state: 'reading' } });
    const expectNoLoad = () => expect(loader.state.get().load).toBeUndefined();
    const expectComplete = () => expect(loader.state.get().saved[HUMAN]).toMatchObject({ complete: true, size: 5000 });
    return { loader, s, startAndStop, expectReading, expectNoLoad, expectComplete };
  }

  it('waits for a stopped download to release the ID instead of refusing a follow-up Download or Continue', async () => {
    const { loader, s, startAndStop, expectReading, expectNoLoad, expectComplete } = setupStoppedDownload();
    const release = await startAndStop();
    expectNoLoad();

    const callsBefore = s.calls.length;
    try {
      loader.download(HUMAN);
      await new Promise<void>((resolve) => setTimeout(resolve, 0));
      expectReading();
      expect(s.calls).toHaveLength(callsBefore);
    } finally {
      release();
    }
    await loader.idle();
    expect(s.calls).toHaveLength(callsBefore + 1);
    expect(s.calls.at(-1)!.headers['Range']).toBe('bytes=4999-');
    expectNoLoad();
    expectComplete();
  });

  it('sends no request and leaves no writing entry when a Download waiting on a stopped download is cancelled', async () => {
    const t = setupStoppedDownload();
    const release = await t.startAndStop();

    const callsBefore = t.s.calls.length;
    try {
      t.loader.download(HUMAN);
      t.expectReading();
      t.loader.cancel();
    } finally {
      release();
    }
    await t.loader.idle();
    expect(t.s.calls).toHaveLength(callsBefore);
    t.expectNoLoad();
    expect(t.loader.state.get().saved[HUMAN]).toMatchObject({ complete: false, savedBytes: 5000 });

    // The ID can be downloaded again after the waiting download was cancelled.
    t.loader.download(HUMAN);
    await t.loader.idle();
    expect(t.s.calls).toHaveLength(callsBefore + 1);
    t.expectNoLoad();
    t.expectComplete();
  });

  it('keeps the saving state of a later save for the same ID when an earlier save finishes its final reloadList', async () => {
    const { loader, store } = setup();
    const pauseAppend = gateStoreAppend(store);
    const listGate = makeGate();
    const origList = store.listCounts.bind(store);
    store.listCounts = async () => {
      await listGate.enter();
      return origList();
    };
    const saving562 = () => loader.state.get().saving[562];

    const list = listGate.arm();
    const first = loader.saveMade(562, new File([bytesOf(300, 1)], 'a.obinary'), madeMeta(1));
    await list.started;
    expect(saving562()).toBeUndefined();

    const append = pauseAppend();
    const second = loader.saveMade(562, new File([bytesOf(400, 2)], 'b.obinary'), madeMeta(2), { replace: true });
    await append.started;
    const expectSavingRun2 = () => expect(saving562()).toEqual({ run: 2 });
    expectSavingRun2();

    list.release();
    try {
      expect(await first).toMatchObject({ ok: true });
      expectSavingRun2();
    } finally {
      append.release();
    }
    expect(await second).toMatchObject({ ok: true, copy: { size: 400, made: madeMeta(2) } });
    await loader.idle();
    expect(saving562()).toBeUndefined();
    expect(loader.state.get().saved[562]).toMatchObject({ complete: true, size: 400, made: madeMeta(2) });
  });

  it('drops the replaced copy from the writing tab list while the Replace writes, and lists the new copy after it ends', async () => {
    const tabs = setupTwoTabs();
    const { loader } = tabs.one;
    expect(await loader.saveMade(562, new File([bytesOf(310, 1)], 'list-a.obinary'), madeMeta(1))).toMatchObject({ ok: true });
    expect(loader.state.get().list).toMatchObject([{ taxid: 562, size: 310, origin: 'made' }]);

    const paused = tabs.pauseNextAppend();
    const replacing = loader.saveMade(562, new File([bytesOf(460, 2)], 'list-b.obinary'), madeMeta(2), { replace: true });
    await paused.started;
    expect(loader.state.get().list).toEqual([]);

    paused.release();
    expect(await replacing).toMatchObject({ ok: true, copy: { taxid: 562, size: 460, made: madeMeta(2) } });
    await loader.idle();
    expect(loader.state.get().list).toMatchObject([{ taxid: 562, size: 460, origin: 'made' }]);
  });

  it('announces on begin only when replacing (also when overwriting an older-URL copy), and drops a stale read that ends after a newer one', async () => {
    const readGates: Array<{ readonly release: () => void; readonly value: Promise<unknown> }> = [];
    const { channel, one, store: memory, other, otherDraft, pauseNextAppend } = setupTwoTabs((base) => ({
      readCounts: (taxid) => {
        const value = base.readCounts(taxid);
        return new Promise((resolve) => readGates.push({ value, release: () => resolve(value) }));
      },
      listCounts: () => base.listCounts(),
      beginCounts: (taxid, start) => base.beginCounts(taxid, start),
      appendCounts: (token, bytes) => base.appendCounts(token, bytes),
      finishCounts: (token) => base.finishCounts(token),
      closeCounts: (token) => base.closeCounts(token),
      deleteCounts: (taxid) => base.deleteCounts(taxid),
    }));
    const announced: number[] = [];
    channel.tab().subscribe((taxid) => announced.push(taxid));

    const waitForReads = async (count: number) => {
      for (let i = 0; i < 100 && readGates.length < count; i++) await new Promise((resolve) => setTimeout(resolve, 0));
      expect(readGates).toHaveLength(count);
    };
    const drainReads = async (count: number) => {
      await waitForReads(count);
      readGates.splice(0).forEach((g) => g.release());
      await other.idle();
    };

    // 1. A first save (not replacing) does not announce until it finishes.
    const firstPause = pauseNextAppend();
    const firstSave = one.loader.saveMade(562, new File([bytesOf(300, 1)], 'first.obinary'), madeMeta(1));
    await firstPause.started;
    expect(announced).toEqual([]);
    firstPause.release();
    expect(await firstSave).toMatchObject({ ok: true });
    expect(announced).toEqual([562]);
    await drainReads(1);

    // 2. Overwriting an older-URL copy (COUNTS_EXISTS -> begin(true)) announces right after begin and again at finish.
    const oldTok = (await memory.beginCounts(7227, { url: 'https://example.org/other.obinary', lastModified: LAST_MODIFIED, size: 100, from: 0 }))!;
    await memory.appendCounts(oldTok, bytesOf(100, 3));
    await memory.finishCounts(oldTok);
    announced.length = 0;
    const oldUrlPause = pauseNextAppend();
    const overwriteOldUrl = one.loader.saveMade(7227, new File([bytesOf(200, 4)], 'old.obinary'), madeMeta(2));
    await oldUrlPause.started;
    expect(announced).toEqual([7227]);
    oldUrlPause.release();
    expect(await overwriteOldUrl).toMatchObject({ ok: true });
    expect(announced).toEqual([7227, 7227]);
    await drainReads(2);

    // 3. On a fast Replace, if the early announce's readCounts resolves late (undefined) after the finish announce's readCounts,
    // tab 2 keeps the completed replacement rather than overwriting it with null.
    const racePause = pauseNextAppend();
    const replacing = one.loader.saveMade(562, new File([bytesOf(450, 5)], 'r.obinary'), madeMeta(3), { replace: true });
    await racePause.started;
    await waitForReads(1);
    expect(await readGates[0]!.value).toBeUndefined();
    racePause.release();
    expect(await replacing).toMatchObject({ ok: true });
    await waitForReads(2);
    const expectCompletedReplacement = () => {
      expect(other.state.get().saved[562]).toMatchObject({ complete: true, size: 450, made: madeMeta(3) });
      expect(otherDraft.get(562)).toMatchObject({ origin: { kind: 'made', taxid: 562, made: madeMeta(3) } });
    };
    // Finish announce's read resolves first (complete 450-byte copy); early announce's read resolves last (undefined).
    readGates[1]!.release();
    await new Promise((resolve) => setTimeout(resolve, 0));
    expectCompletedReplacement();
    readGates[0]!.release();
    await other.idle();
    expectCompletedReplacement();
  });

  it('refuses a save, a Replace, a Delete and a download in place while a queued or running run reads the copy, and keeps the copy', async () => {
    const { loader, downloader, store, draftState } = setup();
    expect(await loader.saveMade(562, new File([bytesOf(300, 1)], 'a.obinary'), madeMeta(1))).toMatchObject({ ok: true });
    store.inUse.add(562);
    const replaceRefusal = 'Not replaced: a run that reads the counts for 562 was queued or running. Save again when it has ended.';
    const deleteRefusal = 'Not removed: a run that reads the counts for 562 was queued or running. Press Delete again when it has ended.';
    const downloadRefusal = 'Not downloaded: a run that reads the counts for 9606 was queued or running. Download again when it has ended.';
    const other = new File([bytesOf(400, 2)], 'b.obinary');
    expect(await loader.saveMade(562, other, madeMeta(2))).toEqual({ ok: false, reason: 'in-use', message: replaceRefusal });
    expect(await loader.saveMade(562, other, madeMeta(2), { replace: true })).toEqual({ ok: false, reason: 'in-use', message: replaceRefusal });
    await loader.delete(562);
    expect(loader.state.get().refusal).toEqual({ taxid: 562, text: deleteRefusal });
    expect(loader.state.get().note).toBeUndefined();
    expect(loader.state.get().saved[562]).toMatchObject({ complete: true, size: 300, made: madeMeta(1) });
    expect(draftState.get(562)).toMatchObject({ origin: { kind: 'made', taxid: 562, made: madeMeta(1) } });
    expect(loader.state.get().list).toMatchObject([{ taxid: 562, size: 300, origin: 'made' }]);
    downloader.next = { ok: false, reason: 'in-use' };
    loader.download(HUMAN);
    await loader.idle();
    expect(loader.state.get().load).toEqual({ taxid: HUMAN, progress: { state: 'failed', message: downloadRefusal } });
    // Once the run ends, Delete removes the copy and the refusal goes.
    store.inUse.delete(562);
    await loader.delete(562);
    expect(loader.state.get().refusal).toBeUndefined();
    expect(loader.state.get().saved[562]).toBeNull();
    expect(draftState.has(562)).toBe(false);
  });

  it('keeps the copy and says why when the storage cannot remove it, and the ID can be used again (no wait on itself)', async () => {
    const { loader, store } = setup();
    expect(await loader.saveMade(562, new File([bytesOf(300, 1)], 'a.obinary'), madeMeta(1))).toMatchObject({ ok: true });
    const remove = store.deleteCounts.bind(store);
    store.deleteCounts = async () => {
      throw new Error('the entry is held');
    };
    await loader.delete(562);
    await loader.idle();
    expect(loader.state.get().note).toEqual({ taxid: 562, text: 'The saved counts could not be removed: the entry is held' });
    expect(loader.state.get().saved[562]).toMatchObject({ complete: true, size: 300 });
    store.deleteCounts = remove;
    expect(await loader.saveMade(562, new File([bytesOf(400, 2)], 'b.obinary'), madeMeta(2), { replace: true })).toMatchObject({ ok: true });
  });

  it('says that the copy of before is no longer kept when a Replace fails after it started, and not for a first save', async () => {
    const { loader, store } = setup();
    expect(await loader.saveMade(562, new File([bytesOf(300, 1)], 'a.obinary'), madeMeta(1))).toMatchObject({ ok: true });
    const full = 'The browser refused to keep the counts: its storage is full.';
    store.appendError = new Error(full);
    expect(await loader.saveMade(562, new File([bytesOf(400, 2)], 'b.obinary'), madeMeta(2), { replace: true })).toEqual({
      ok: false,
      reason: 'failed',
      message: `${full} The counts saved before for 562 are no longer kept.`,
    });
    expect(loader.state.get().saved[562]).toBeNull();
    expect(await loader.saveMade(1280, new File([bytesOf(400, 2)], 'c.obinary'), madeMeta(3))).toEqual({ ok: false, reason: 'failed', message: full });
  });

  it('keeps the list of the latest read when an earlier read of the list ends after it', async () => {
    const memory = new MemoryCountsStore({ persistent: true, now: () => Date.UTC(2026, 10, 12) });
    const gates: Array<() => void> = [];
    const store: CountsStore = {
      readCounts: (taxid) => memory.readCounts(taxid),
      listCounts: () => {
        const list = memory.listCounts();
        return new Promise<readonly SavedCountsEntry[]>((resolve) => gates.push(() => resolve(list)));
      },
      beginCounts: (taxid, start) => memory.beginCounts(taxid, start),
      appendCounts: (token, bytes) => memory.appendCounts(token, bytes),
      finishCounts: (token) => memory.finishCounts(token),
      closeCounts: (token) => memory.closeCounts(token),
      deleteCounts: (taxid) => memory.deleteCounts(taxid),
    };
    const loader = new CountsLoader({ downloader: new FakeCountsDownloader(memory, {}), store, draft: { setCounts: () => undefined } });
    // The tab opens and reads the list (nothing is saved yet); a save then reads it again.
    loader.refreshList();
    const saving = loader.saveMade(562, new File([bytesOf(300, 1)], 'a.obinary'), madeMeta(1));
    for (let i = 0; i < 100 && gates.length < 2; i++) await new Promise((resolve) => setTimeout(resolve, 0));
    expect(gates).toHaveLength(2);
    gates[1]!();
    expect(await saving).toMatchObject({ ok: true });
    expect(loader.state.get().list).toMatchObject([{ taxid: 562, size: 300 }]);
    // The first read, of the empty list, ends last: it is dropped.
    gates[0]!();
    await loader.idle();
    expect(loader.state.get().list).toMatchObject([{ taxid: 562, size: 300 }]);
  });

  it('offers the saved taxonomy ID for "Use in BLASTN" only while the stored copy matches the saved run identity', async () => {
    const { loader } = setup({ files: { 562: COUNTS } });
    const stored = () => loader.state.get().saved[562];

    // Run 1 saves its counts for 562 -> activeSavedCountsTaxid offers 562.
    const r1 = await loader.saveMade(562, new File([bytesOf(300, 1)], 'r1.obinary'), madeMeta(1));
    expect(r1.ok).toBe(true);
    const run1Saved = r1.ok ? savedTabCopyOf(r1.copy) : undefined;
    expect(run1Saved).toEqual({ taxid: 562, run: 1, madeAt: MADE_AT, size: 300 });
    expect(activeSavedCountsTaxid(run1Saved, stored())).toBe(562);

    // Run 2 replaces 562 (even with identical size and madeAt) -> Run 1 falls back to its own file, Run 2 offers 562.
    const r2 = await loader.saveMade(562, new File([bytesOf(300, 2)], 'r2.obinary'), madeMeta(2), { replace: true });
    expect(r2.ok).toBe(true);
    const run2Saved = r2.ok ? savedTabCopyOf(r2.copy) : undefined;
    expect(activeSavedCountsTaxid(run1Saved, stored())).toBeUndefined();
    expect(activeSavedCountsTaxid(run2Saved, stored())).toBe(562);

    // Deleting or replacing with an NCBI download invalidates Run 2's saved taxonomy ID too.
    await loader.delete(562);
    expect(activeSavedCountsTaxid(run2Saved, stored())).toBeUndefined();
    loader.download(562);
    await loader.idle();
    expect(stored()?.complete).toBe(true);
    expect(activeSavedCountsTaxid(run2Saved, stored())).toBeUndefined();
  });
});

// --- the OPFS store -----------------------------------------------------------------------------

describe('OpfsCountsStore', () => {
  async function make() {
    const opfs = new MemoryOpfs();
    const root = opfs.root as unknown as FileSystemDirectoryHandle;
    let tokens = 0;
    const store = new OpfsCountsStore(root, opfs, () => `t${++tokens}`, () => Date.UTC(2026, 10, 22));
    return { opfs, store };
  }
  const start = (from = 0) => ({ url: countsUrl(HUMAN), lastModified: LAST_MODIFIED, size: 5000, from });

  it('writes the file part by part with the saved bytes in meta.json at windowmasker/<taxid>/, and is whole at the end', async () => {
    const { opfs, store } = await make();
    const token = (await store.beginCounts(HUMAN, start()))!;
    await store.appendCounts(token, COUNTS.bytes.slice(0, 3000));
    const meta = () => JSON.parse(new TextDecoder().decode(opfs.read(['9606', 'meta.json'])));
    expect(meta()).toMatchObject({ taxid: HUMAN, url: countsUrl(HUMAN), lastModified: LAST_MODIFIED, size: 5000, savedBytes: 3000, complete: false });
    expect(meta().made).toBeUndefined();
    expect(await store.readCounts(HUMAN)).toMatchObject({ origin: 'downloaded', copy: { savedBytes: 3000, complete: false } });
    await store.appendCounts(token, COUNTS.bytes.slice(3000));
    const saved = await store.finishCounts(token);
    expect(saved.origin).toBe('downloaded');
    expect(saved.copy).toMatchObject({ savedBytes: 5000, complete: true });
    expect(new Uint8Array(await saved.file!.arrayBuffer())).toEqual(COUNTS.bytes);
    expect(opfs.read(['9606', 'wmasker.obinary'])).toEqual(COUNTS.bytes);
    expect(meta().complete).toBe(true);
  });

  it('reads a WW2 meta.json (without made) as a download without rewriting it', async () => {
    const { opfs, store } = await make();
    const dir = await opfs.root.getDirectoryHandle('9606', { create: true });
    const legacyMeta = JSON.stringify({
      taxid: HUMAN,
      url: countsUrl(HUMAN),
      lastModified: LAST_MODIFIED,
      size: 5000,
      savedBytes: 5000,
      complete: true,
      savedAt: Date.UTC(2026, 9, 15),
    });
    const metaHandle = await (await dir.getFileHandle('meta.json', { create: true })).createSyncAccessHandle();
    metaHandle.write(new TextEncoder().encode(legacyMeta), { at: 0 });
    metaHandle.close();
    const dataHandle = await (await dir.getFileHandle('wmasker.obinary', { create: true })).createSyncAccessHandle();
    dataHandle.write(COUNTS.bytes, { at: 0 });
    dataHandle.close();

    const saved = await store.readCounts(HUMAN);
    expect(saved).toMatchObject({ origin: 'downloaded', copy: { taxid: HUMAN, url: countsUrl(HUMAN), lastModified: LAST_MODIFIED, complete: true } });
    expect(saved!.copy.made).toBeUndefined();
    expect(new TextDecoder().decode(opfs.read(['9606', 'meta.json']))).toBe(legacyMeta);
  });

  it('writes a made copy with empty url/lastModified and made metadata, which an older dev reader sees as a different URL', async () => {
    const { opfs, store } = await make();
    const meta = madeMeta(4);
    const bytes = bytesOf(2048, 11);
    const token = (await store.beginCounts(562, { made: meta, size: bytes.length, from: 0 }))!;
    await store.appendCounts(token, bytes.slice(0, 1000));
    await store.appendCounts(token, bytes.slice(1000));
    const saved = await store.finishCounts(token);
    expect(saved.origin).toBe('made');
    expect(saved.copy).toMatchObject({ taxid: 562, url: '', lastModified: '', size: bytes.length, savedBytes: bytes.length, complete: true, made: meta });
    expect(new Uint8Array(await saved.file!.arrayBuffer())).toEqual(bytes);

    // Older-release compatibility (§4.1): the previous WW2 reader checked `copy.url !== countsUrl(taxid)`,
    // so it rejects the made copy as "a different URL" (offers Download, never takes it for NCBI's file).
    const rawMeta = JSON.parse(new TextDecoder().decode(opfs.read(['562', 'meta.json']))) as CountsCopy;
    const oldDevReader = (copy: CountsCopy, expectedUrl: string) => (copy.url !== expectedUrl ? undefined : copy);
    expect(rawMeta.url).toBe('');
    expect(rawMeta.lastModified).toBe('');
    expect(oldDevReader(rawMeta, countsUrl(562))).toBeUndefined();
  });

  it('replaces a download, a partial download, or a made copy, and refuses with COUNTS_EXISTS when replace is false', async () => {
    const { store } = await make();
    const dl = (await store.beginCounts(HUMAN, start()))!;
    await store.appendCounts(dl, COUNTS.bytes.slice(0, 1000));
    await store.closeCounts(dl);
    // Partial download triggers COUNTS_EXISTS when replace is false.
    await expect(store.beginCounts(HUMAN, { made: madeMeta(1), size: 100, from: 0, replace: false })).rejects.toMatchObject({ name: COUNTS_EXISTS });
    // Replace partial download with a made copy.
    const m1 = (await store.beginCounts(HUMAN, { made: madeMeta(1), size: 100, from: 0, replace: true }))!;
    await store.appendCounts(m1, bytesOf(100, 1));
    await store.finishCounts(m1);
    expect((await store.readCounts(HUMAN))!.origin).toBe('made');
    // Replace made copy with another made copy.
    await expect(store.beginCounts(HUMAN, { made: madeMeta(2), size: 200, from: 0, replace: false })).rejects.toMatchObject({ name: COUNTS_EXISTS });
    const m2 = (await store.beginCounts(HUMAN, { made: madeMeta(2), size: 200, from: 0, replace: true }))!;
    await store.appendCounts(m2, bytesOf(200, 2));
    const saved = await store.finishCounts(m2);
    expect(saved.copy).toMatchObject({ size: 200, made: madeMeta(2) });
    // Complete download replaced by made copy.
    const dl2 = (await store.beginCounts(562, { url: countsUrl(562), lastModified: LAST_MODIFIED, size: 5000, from: 0 }))!;
    await store.appendCounts(dl2, COUNTS.bytes);
    await store.finishCounts(dl2);
    await expect(store.beginCounts(562, { made: madeMeta(3), size: 300, from: 0, replace: false })).rejects.toMatchObject({ name: COUNTS_EXISTS });
    const m3 = (await store.beginCounts(562, { made: madeMeta(3), size: 300, from: 0, replace: true }))!;
    await store.appendCounts(m3, bytesOf(300, 3));
    expect((await store.finishCounts(m3)).copy).toMatchObject({ size: 300, made: madeMeta(3) });
  });

  it('deletes an interrupted made copy on next read, while leaving an active write untouched', async () => {
    const { opfs, store } = await make();
    const token = (await store.beginCounts(562, { made: madeMeta(4), size: 1000, from: 0 }))!;
    await store.appendCounts(token, bytesOf(400, 3));
    // While the write holds the lock, readCounts returns undefined without removing the directory.
    expect(await store.readCounts(562)).toBeUndefined();
    expect(opfs.has(['562'])).toBe(true);
    // Once the lock is gone without finishCounts (interrupted), the next read deletes it.
    await store.closeCounts(token);
    expect(await store.readCounts(562)).toBeUndefined();
    expect(opfs.has(['562'])).toBe(false);
  });

  it('lists every complete copy (downloaded and made) and skips partials and interrupted made copies', async () => {
    const { opfs, store } = await make();
    // Complete download for 9606.
    const d1 = (await store.beginCounts(HUMAN, start()))!;
    await store.appendCounts(d1, COUNTS.bytes);
    await store.finishCounts(d1);
    // Complete made copy for 562.
    const m1 = (await store.beginCounts(562, { made: madeMeta(4), size: 700, from: 0 }))!;
    await store.appendCounts(m1, bytesOf(700, 5));
    await store.finishCounts(m1);
    // Partial download for 7227 (omitted from list, stays on disk).
    const p1 = (await store.beginCounts(7227, { url: countsUrl(7227), lastModified: LAST_MODIFIED, size: 5000, from: 0 }))!;
    await store.appendCounts(p1, COUNTS.bytes.slice(0, 1000));
    await store.closeCounts(p1);
    // Interrupted made copy for 10090 (omitted and cleaned up).
    const b1 = (await store.beginCounts(10090, { made: madeMeta(5), size: 800, from: 0 }))!;
    await store.appendCounts(b1, bytesOf(100, 6));
    await store.closeCounts(b1);

    const list = await store.listCounts();
    expect(list.map(({ taxid, name, size, origin, date }) => ({ taxid, name, size, origin, date }))).toEqual([
      { taxid: 562, name: COUNTS_FILE, size: 700, origin: 'made', date: '2026-11-22' },
      { taxid: HUMAN, name: COUNTS_FILE, size: 5000, origin: 'downloaded', date: '2026-10-01' },
    ]);
    expect(opfs.has(['7227'])).toBe(true);
    expect(opfs.has(['10090'])).toBe(false);
  });

  it('continues a saved part at its byte (or before it), and refuses a part that is not the one asked for', async () => {
    const { store } = await make();
    const token = (await store.beginCounts(HUMAN, start()))!;
    await store.appendCounts(token, COUNTS.bytes.slice(0, 2000));
    await store.closeCounts(token);
    await expect(store.beginCounts(HUMAN, start(3000))).rejects.toMatchObject({ name: COUNTS_CHANGED });
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
    expect(await store.readCounts(HUMAN)).toMatchObject({ copy: { savedBytes: 0, size: 6000, complete: false } });
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
    expect(await store.readCounts(HUMAN)).toBeUndefined();
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
    expect(await store.readCounts(HUMAN)).toBeUndefined();
  });

  it('refuses a made copy, a download and a deletion with COUNTS_IN_USE while searches hold the ID in shared mode, and leaves the copy', async () => {
    const { opfs, store } = await make();
    const a = (await store.beginCounts(HUMAN, start()))!;
    await store.appendCounts(a, COUNTS.bytes);
    await store.finishCounts(a);
    const locks = webCountsLocks(opfs);
    const first = (await locks.holdShared(HUMAN))!;
    const second = (await locks.holdShared(HUMAN))!;
    await expect(store.beginCounts(HUMAN, { made: madeMeta(1), size: 10, from: 0, replace: true })).rejects.toMatchObject({ name: COUNTS_IN_USE });
    await expect(store.beginCounts(HUMAN, start())).rejects.toMatchObject({ name: COUNTS_IN_USE });
    await expect(store.deleteCounts(HUMAN)).rejects.toMatchObject({ name: COUNTS_IN_USE });
    expect(opfs.read(['9606', 'wmasker.obinary'])).toEqual(COUNTS.bytes);
    expect(await store.readCounts(HUMAN)).toMatchObject({ copy: { complete: true, size: 5000 } });
    // Another ID is free.
    const other = (await store.beginCounts(562, { made: madeMeta(1), size: 0, from: 0 }))!;
    await store.closeCounts(other);
    // The ID is free once every holder let go (a release is let go once).
    first();
    first();
    await new Promise((resolve) => setTimeout(resolve, 0));
    await expect(store.deleteCounts(HUMAN)).rejects.toMatchObject({ name: COUNTS_IN_USE });
    second();
    await new Promise((resolve) => setTimeout(resolve, 0));
    await store.deleteCounts(HUMAN);
    expect(opfs.has(['9606'])).toBe(false);
  });

  it('gives no shared hold while a writing holds the ID, and a deletion waits for the writing', async () => {
    const { opfs, store } = await make();
    const locks = webCountsLocks(opfs);
    const token = (await store.beginCounts(HUMAN, start()))!;
    expect(await locks.holdShared(HUMAN)).toBeUndefined();
    let deleted = false;
    const deletion = store.deleteCounts(HUMAN).then(() => {
      deleted = true;
    });
    await new Promise((resolve) => setTimeout(resolve, 0));
    expect(deleted).toBe(false);
    // A search asks while the deletion waits: still refused.
    expect(await locks.holdShared(HUMAN)).toBeUndefined();
    await store.appendCounts(token, COUNTS.bytes.slice(0, 10));
    await store.closeCounts(token);
    await deletion;
    expect(opfs.has(['9606'])).toBe(false);
  });

  it('refuses a made copy that the storage has no room for before changing anything, counting the copy it replaces', async () => {
    const opfs = new MemoryOpfs();
    let tokens = 0;
    let estimate: { usage: number; quota: number } | undefined = { usage: 1000, quota: 1300 };
    const store = new OpfsCountsStore(
      opfs.root as unknown as FileSystemDirectoryHandle,
      opfs,
      () => `t${++tokens}`,
      () => Date.UTC(2026, 10, 22),
      async () => estimate,
    );
    const made = (run: number, size: number) => ({ made: madeMeta(run), size, from: 0 as const, replace: true });
    // 300 bytes free: a copy of 301 is refused before creating the folder, one of 300 is written.
    await expect(store.beginCounts(562, made(1, 301))).rejects.toThrow("The browser's storage has no room for these counts: nothing was changed.");
    expect(opfs.has(['562'])).toBe(false);
    const t1 = (await store.beginCounts(562, made(1, 300)))!;
    await store.appendCounts(t1, bytesOf(300, 1));
    await store.finishCounts(t1);
    // The storage is full now; a Replace of 300 fits in the bytes of the copy it replaces, one of 301 does not.
    estimate = { usage: 1300, quota: 1300 };
    await expect(store.beginCounts(562, made(2, 301))).rejects.toThrow('nothing was changed');
    expect(await store.readCounts(562)).toMatchObject({ copy: { complete: true, size: 300, made: madeMeta(1) } });
    expect(opfs.read(['562', 'wmasker.obinary'])).toEqual(bytesOf(300, 1));
    const t2 = (await store.beginCounts(562, made(2, 300)))!;
    await store.closeCounts(t2);
    // Without an estimate, the writing starts.
    estimate = undefined;
    const t3 = (await store.beginCounts(562, made(3, 5000)))!;
    await store.closeCounts(t3);
  });
});
