// makeblastdb in the browsers through the package's API, without the application's screens (line NP WW4b1, D-WW4b-1 and
// D-WW4b-2, 2026-11-27): the real Wasm engine runs makeblastdb (`WasmEngine.runTool`) on Files of the repository's fixtures
// as its run's directory, under two Playwright `timezoneId`s (one with summer time, one without and half an hour off), and
// what it wrote (its exit status, standard output, standard error and every file it left, by name) is the native command
// line's, run in a folder of the same files with the same zone (`TZ`): byte for byte but the seconds of the report's and
// the mask column's times and the elapsed time, and the browser's run equals the CLI's run either just before or just
// after it (a minute may pass between them). Then the database that the browser made is searched through the package's
// `-db` path, and its outfmt 0 and 6 are those of the native command line searching the same files. Rows: version 5
// nucleotide with `-parse_seqids`; protein with `-mask_data` (seg); a nucleotide database of three volumes (`-max_file_sz`);
// a bacterial genome, whose sequence file is larger than one frame.
//
// LOSAT_WEB_MEASURE=1 adds a measurement in Chromium (outside CI): makeblastdb of the files named by
// LOSAT_WEB_MEASURE_MAKEBLASTDB_NUCL (a nucleotide FASTA, 50 MB or more) and LOSAT_WEB_MEASURE_MAKEBLASTDB_PROT (a protein
// FASTA of SwissProt's size), picked from the disk as a user picks them, next to the native command line: the wall times,
// the native max RSS, the browser's process tree peak RSS, the page's JS heap, the engine's linear memory after the run (the
// memory before the next search on the same instance) and the bytes of the database (the run's directory's overlay).
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { copyFileSync, existsSync, mkdirSync, mkdtempSync, readdirSync, readFileSync, rmSync, statSync, symlinkSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { basename, dirname, join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import type { MakeblastdbCase, MakeblastdbResult } from './harness/engine';
import { buildHarness } from './support/browser';
import { ENGINE, NO_ENGINE_REASON } from './support/engine';
import { openHarness, REPOSITORY, startHarnessServer, type HarnessServer } from './support/harness-server';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { test } from './support/test';

test.skip(ENGINE === undefined, NO_ENGINE_REASON);
test.skip(NATIVE === undefined, NO_NATIVE_REASON);
test.setTimeout(600_000);

const LOSAT = join(REPOSITORY, 'LOSAT');
const MEASURE = process.env['LOSAT_WEB_MEASURE'] === '1';
const MIB = 1024 * 1024;

/** A zone with summer time and one without, half an hour off a whole hour. */
const ZONES = ['America/New_York', 'Asia/Kolkata'] as const;
/** The names by which a browser's `Intl` may give a zone (Chromium's ICU gives Asia/Kolkata its older name). */
const ZONE_NAMES: Readonly<Record<string, readonly string[]>> = { 'America/New_York': ['America/New_York'], 'Asia/Kolkata': ['Asia/Kolkata', 'Asia/Calcutta'] };

const CASES: readonly MakeblastdbCase[] = [
  {
    id: 'nucl.v5.parse_seqids',
    argv: ['makeblastdb', '-in', 'ns.fna', '-dbtype', 'nucl', '-parse_seqids', '-out', 'db'],
    files: [['ns.fna', 'tests/dropin_corpus/ns.fna']],
    search: { program: 'blastn', argv: ['blastn', '-db', 'db', '-query', 'nq.fna'], query: 'tests/dropin_corpus/nq.fna' },
  },
  {
    id: 'prot.mask_data',
    argv: ['makeblastdb', '-in', 'mask/prot.faa', '-dbtype', 'prot', '-mask_data', 'mask/seg.asnb', '-out', 'pdb'],
    files: [
      ['mask/prot.faa', 'tests/seqdb/mask_prot.faa'],
      ['mask/seg.asnb', 'tests/seqdb/mask_seg.asnb'],
    ],
    search: { program: 'blastp', argv: ['blastp', '-db', 'pdb', '-query', 'q.faa', '-db_soft_mask', '21'], query: 'tests/seqdb/mask_prot.faa' },
  },
  {
    id: 'nucl.volumes',
    argv: ['makeblastdb', '-in', 'ns.fna', '-dbtype', 'nucl', '-max_file_sz', '4000', '-out', 'vol/db'],
    files: [['ns.fna', 'tests/dropin_corpus/ns.fna']],
    search: { program: 'blastn', argv: ['blastn', '-db', 'vol/db', '-query', 'nq.fna'], query: 'tests/dropin_corpus/nq.fna' },
  },
  {
    // A genome (E. coli K-12, 4.6 Mbp): its .nsq (1.2 MB) is larger than a chunk of the run's directory (1 MiB), so it
    // reaches the page in two frames, each a transferred buffer.
    id: 'nucl.genome',
    argv: ['makeblastdb', '-in', 'MG1655.fna', '-dbtype', 'nucl', '-out', 'genome'],
    files: [['MG1655.fna', 'tests/fasta/MG1655.fna']],
    search: { program: 'blastn', argv: ['blastn', '-db', 'genome', '-query', 'nq.fna'], query: 'tests/dropin_corpus/nq.fna' },
  },
];

/** What a run wrote: its streams and the files of its folder that it made or changed (by their names in the folder). */
interface Written {
  readonly status: number;
  readonly stdout: Buffer;
  readonly stderr: Buffer;
  readonly files: ReadonlyMap<string, Buffer>;
}

/** The files under `folder`, by their paths in it, in the byte order of the names. */
function treeOf(folder: string, prefix = ''): Map<string, Buffer> {
  const files = new Map<string, Buffer>();
  for (const name of readdirSync(join(folder, prefix)).sort()) {
    const path = prefix === '' ? name : `${prefix}/${name}`;
    if (statSync(join(folder, path)).isDirectory()) for (const [inner, bytes] of treeOf(folder, path)) files.set(inner, bytes);
    else files.set(path, readFileSync(join(folder, path)));
  }
  return files;
}

/** A folder holding `files` (each [name in the folder, path under LOSAT/]). */
function folderOf(files: ReadonlyArray<readonly [string, string]>): string {
  const folder = mkdtempSync(join(tmpdir(), 'losat-makeblastdb-'));
  for (const [name, path] of files) {
    mkdirSync(dirname(join(folder, name)), { recursive: true });
    copyFileSync(join(LOSAT, path), join(folder, name));
  }
  return folder;
}

/** The native command line's run of `item` in a folder of its files with `TZ=zone`, its paths written as the engine's `/work`. */
function nativeRun(item: MakeblastdbCase, zone: string): Written {
  const folder = folderOf(item.files);
  try {
    const given = treeOf(folder);
    const run = spawnSync(NATIVE!, item.argv, { cwd: folder, env: { ...process.env, TZ: zone }, maxBuffer: 1 << 28 });
    if (run.error) throw run.error;
    const place = (bytes: Buffer) => Buffer.from(bytes.toString('latin1').split(folder).join('/work'), 'latin1');
    const files = new Map<string, Buffer>();
    for (const [name, bytes] of treeOf(folder)) if (!given.get(name)?.equals(bytes)) files.set(name, place(bytes));
    return { status: run.status ?? -1, stdout: place(run.stdout), stderr: place(run.stderr), files };
  } finally {
    rmSync(folder, { recursive: true, force: true });
  }
}

/** The bytes with the seconds of the times that makeblastdb prints (MM/DD/YYYY HH:MM:SS) and its elapsed time masked. */
function withoutSeconds(bytes: Buffer): string {
  return bytes
    .toString('latin1')
    .replace(/(\d{2}\/\d{2}\/\d{4} \d{2}:\d{2}):\d{2}/g, '$1:##')
    .replace(/( sequences? in )[0-9.e+-]+( seconds)/g, '$1#$2');
}

function maskedOf(written: Written): string[] {
  return [
    `status ${written.status}`,
    `stdout ${withoutSeconds(written.stdout)}`,
    `stderr ${withoutSeconds(written.stderr)}`,
    ...[...written.files].map(([name, bytes]) => `${name} ${createHash('sha256').update(withoutSeconds(bytes), 'latin1').digest('hex')}`),
  ];
}

const fromBase64 = (text: string | undefined) => Buffer.from(text ?? '', 'base64');

/** The native search of `item.search` (outfmt 0 or 6, one thread) in a folder of the files that the browser made. */
function nativeSearch(item: MakeblastdbCase, files: ReadonlyMap<string, Buffer>, format: 0 | 6): string {
  const folder = mkdtempSync(join(tmpdir(), 'losat-makeblastdb-search-'));
  try {
    for (const [name, bytes] of files) {
      mkdirSync(dirname(join(folder, name)), { recursive: true });
      writeFileSync(join(folder, name), bytes);
    }
    const search = item.search!;
    const query = search.argv[search.argv.indexOf('-query') + 1]!;
    copyFileSync(join(LOSAT, search.query), join(folder, query));
    const run = spawnSync(NATIVE!, [...search.argv, '-outfmt', String(format), '-num_threads', '1'], { cwd: folder, maxBuffer: 1 << 28 });
    if (run.error) throw run.error;
    if (run.status !== 0) throw new Error(`native ${search.argv.join(' ')} exited ${run.status}: ${run.stderr.toString()}`);
    return run.stdout.toString('utf8');
  } finally {
    rmSync(folder, { recursive: true, force: true });
  }
}

let server: HarnessServer;

test.beforeAll(async () => {
  server = await startHarnessServer(await buildHarness());
});

test.afterAll(async () => {
  await server?.close();
});

for (const zone of ZONES) {
  test.describe(`in the time zone ${zone}`, () => {
    test.use({ timezoneId: zone });

    test(`makeblastdb writes the command line's database, and a -db search of it writes the command line's bytes (${zone})`, async ({ context, browserName }) => {
      const page = await openHarness(context, server);
      for (const item of CASES) {
        const label = `${browserName} ${zone} ${item.id}`;
        const before = nativeRun(item, zone);
        const [result]: MakeblastdbResult[] = await page.evaluate((list) => window.losatHarness!.engine.makeblastdbRuns(list), [item]);
        const after = nativeRun(item, zone);
        expect(result!.error, label).toBeUndefined();
        expect(ZONE_NAMES[zone], `${label}: the page's zone`).toContain(result!.zone);
        const browser: Written = {
          status: result!.status!,
          stdout: fromBase64(result!.stdout),
          stderr: fromBase64(result!.stderr),
          files: new Map(result!.files!.map((file) => [file.name, fromBase64(file.base64)])),
        };
        expect(before.status, `${label}: the CLI's run`).toBe(0);
        expect([...browser.files.keys()], `${label}: the files, by name`).toEqual([...before.files.keys()]);
        const masked = maskedOf(browser);
        const matching = [before, after].some((native) => JSON.stringify(maskedOf(native)) === JSON.stringify(masked));
        if (!matching) expect(masked, `${label}: as the CLI's run`).toEqual(maskedOf(before));
        console.log(`${label}: ${browser.files.size} files, ${Math.round(result!.toolMs ?? 0)} ms`);

        const search = result!.search!;
        expect(search.error, `${label}: the -db search`).toBeUndefined();
        expect(fromBase64(search.out6).toString('utf8'), `${label}: -db search, outfmt 6`).toBe(nativeSearch(item, browser.files, 6));
        expect(fromBase64(search.out0).toString('utf8'), `${label}: -db search, outfmt 0`).toBe(nativeSearch(item, browser.files, 0));
      }
    });
  });
}

// --- the measurement (LOSAT_WEB_MEASURE=1, Chromium, outside CI) ------------------------------------------------------------

/** Sums VmRSS (bytes) of every descendant process of `rootPid` (Linux `/proc`): the browser's processes. */
function descendantRssBytes(rootPid: number): number {
  const children = new Map<number, number[]>();
  const rss = new Map<number, number>();
  for (const entry of readdirSync('/proc')) {
    if (!/^\d+$/.test(entry)) continue;
    let status: string;
    try {
      status = readFileSync(`/proc/${entry}/status`, 'utf8');
    } catch {
      continue;
    }
    const parent = /^PPid:\s*(\d+)/m.exec(status);
    const resident = /^VmRSS:\s*(\d+)\s*kB/m.exec(status);
    if (parent !== null) children.set(Number(parent[1]), [...(children.get(Number(parent[1])) ?? []), Number(entry)]);
    if (resident !== null) rss.set(Number(entry), Number(resident[1]) * 1024);
  }
  let total = 0;
  const queue = [...(children.get(rootPid) ?? [])];
  for (let pid = queue.shift(); pid !== undefined; pid = queue.shift()) {
    total += rss.get(pid) ?? 0;
    queue.push(...(children.get(pid) ?? []));
  }
  return total;
}

/** The native command line's makeblastdb of `input` with GNU time: wall seconds, max RSS and the files' SHA-256. */
function nativeMeasured(argv: readonly string[], input: string): { seconds: number; maxRssBytes: number; files: Map<string, { size: number; sha256: string; mibSha256: string[] }> } {
  const folder = mkdtempSync(join(tmpdir(), 'losat-makeblastdb-measure-'));
  try {
    const report = join(folder, 'time.txt');
    // The input under its name in the folder, as the browser's run has it (the title of the database is the name).
    symlinkSync(input, join(folder, basename(input)));
    const run = spawnSync('/usr/bin/time', ['-f', '%e %M', '-o', report, NATIVE!, ...argv], { cwd: folder, maxBuffer: 1 << 28 });
    if (run.status !== 0) throw new Error(`native makeblastdb exited ${run.status}: ${run.stderr.toString()}`);
    rmSync(join(folder, basename(input)));
    const [seconds, kib] = readFileSync(report, 'utf8').trim().split(/\s+/).map(Number);
    const files = new Map<string, { size: number; sha256: string; mibSha256: string[] }>();
    const digest = (bytes: Uint8Array) => createHash('sha256').update(bytes).digest('hex');
    for (const [name, bytes] of treeOf(folder)) {
      if (name === 'time.txt') continue;
      const mibSha256: string[] = [];
      for (let at = 0; at < bytes.length; at += MIB) mibSha256.push(digest(bytes.subarray(at, at + MIB)));
      files.set(name, { size: bytes.length, sha256: digest(bytes), mibSha256 });
    }
    return { seconds: seconds!, maxRssBytes: kib! * 1024, files };
  } finally {
    rmSync(folder, { recursive: true, force: true });
  }
}

/** A file input of the harness page, `#measure-pick`, set from the disk. */
async function pick(page: Page, path: string): Promise<void> {
  await page.evaluate(() => {
    if (document.querySelector('#measure-pick') !== null) return;
    document.body.append(Object.assign(document.createElement('input'), { type: 'file', id: 'measure-pick' }));
  });
  await page.setInputFiles('#measure-pick', path);
}

const MEASURED = [
  { id: 'nucl', variable: 'LOSAT_WEB_MEASURE_MAKEBLASTDB_NUCL', dbtype: 'nucl', search: { program: 'blastn', argv: ['blastn', '-db', 'big', '-query', 'nq.fna'], query: 'tests/dropin_corpus/nq.fna' } },
  { id: 'prot', variable: 'LOSAT_WEB_MEASURE_MAKEBLASTDB_PROT', dbtype: 'prot', search: { program: 'blastp', argv: ['blastp', '-db', 'big', '-query', 'pq.faa'], query: 'tests/dropin_corpus/pq.faa' } },
] as const;

test.describe('makeblastdb of a large input, measured', () => {
  test.skip(!MEASURE, 'a local measurement: set LOSAT_WEB_MEASURE=1');

  for (const item of MEASURED) {
    test(`measure makeblastdb -dbtype ${item.dbtype}`, async ({ context, browserName }) => {
      test.skip(browserName !== 'chromium', 'Chromium only');
      test.setTimeout(3_600_000);
      const input = process.env[item.variable];
      test.skip(input === undefined || !existsSync(input), `set ${item.variable} to a local FASTA file`);
      const name = basename(input!);
      const argv = ['makeblastdb', '-in', name, '-dbtype', item.dbtype, '-out', 'big'];
      const native = nativeMeasured(argv, input!);

      const page = await openHarness(context, server);
      await pick(page, input!);
      const heap = () => page.evaluate(() => (performance as Performance & { memory?: { usedJSHeapSize: number } }).memory?.usedJSHeapSize ?? 0);
      const rssBefore = descendantRssBytes(process.pid);
      let peakRss = rssBefore;
      let peakHeap = await heap();
      const sampler = setInterval(() => {
        peakRss = Math.max(peakRss, descendantRssBytes(process.pid));
        heap().then((bytes) => (peakHeap = Math.max(peakHeap, bytes)), () => undefined);
      }, 250);
      let results: MakeblastdbResult[];
      try {
        // A run of `-version` first loads and compiles the engine, so that the measured run is the build alone.
        const warm: MakeblastdbCase = { id: 'measure.warm-up', argv: ['makeblastdb', '-version'], files: [] };
        const run: MakeblastdbCase = { id: `measure.${item.id}`, argv, files: [[name, '']], picked: '#measure-pick', digestsOnly: true, search: item.search };
        results = await page.evaluate((list) => window.losatHarness!.engine.makeblastdbRuns(list), [warm, run]);
      } finally {
        clearInterval(sampler);
      }
      const result = results[1];
      console.log(`warm-up (load and compile the engine, makeblastdb -version): ${Math.round(results[0]?.toolMs ?? 0)} ms`);
      expect(result!.error).toBeUndefined();
      expect(result!.status).toBe(0);
      const files = result!.files!;
      expect(files.map((file) => file.name)).toEqual([...native.files.keys()]);
      // The files without a date are the same bytes; those with one (the index and the JSON) the same length.
      for (const file of files) {
        const expected = native.files.get(file.name)!;
        if (/\.(nin|pin|njs|pjs|nal|pal)$/.test(file.name)) {
          expect.soft(file.size, file.name).toBe(expected.size);
          continue;
        }
        const differing = (file.mibSha256 ?? []).flatMap((sha, i) => (sha === expected.mibSha256[i] ? [] : [i]));
        expect.soft({ size: file.size, sha256: file.sha256, differingMib: differing.slice(0, 20), differingMibCount: differing.length }, file.name).toEqual({
          size: expected.size,
          sha256: expected.sha256,
          differingMib: [],
          differingMibCount: 0,
        });
      }
      const overlayBytes = files.reduce((sum, file) => sum + file.size, 0);
      const engineBytes = result!.search?.runtime?.memory?.linearBytesBefore;
      const line =
        `makeblastdb -dbtype ${item.dbtype} of ${(statSync(input!).size / MIB).toFixed(1)} MiB (${name}): ` +
        `browser ${((result!.toolMs ?? 0) / 1000).toFixed(2)} s, native ${native.seconds.toFixed(2)} s; ` +
        `native max RSS ${(native.maxRssBytes / MIB).toFixed(1)} MiB; browser processes peak RSS ${(peakRss / MIB).toFixed(1)} MiB ` +
        `(+${((peakRss - rssBefore) / MIB).toFixed(1)} MiB over the start); page JS heap peak ${(peakHeap / MIB).toFixed(1)} MiB; ` +
        `engine linear memory after the run ${engineBytes === undefined ? 'n/a' : `${(engineBytes / MIB).toFixed(1)} MiB`}; ` +
        `database (overlay) ${(overlayBytes / MIB).toFixed(1)} MiB in ${files.length} files`;
      console.log(line);
      test.info().annotations.push({ type: 'measurement', description: line });
    });
  }
});
