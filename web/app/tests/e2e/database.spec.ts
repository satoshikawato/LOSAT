// -db in the browsers through the package's API, without the application's screens (line WQA D12a, 2026-11-20): the
// db_search rows that the wasm-threading check runs (LOSAT/tests/check_wasm_threading.py) and two of the filters
// (-seqidlist, -taxids), each with the files of its database's folder given to the run as Files, on the serial module
// and on the threaded one at 4 threads, whose threads read the descriptors that another thread opened. Outfmt 6 is
// NCBI BLAST+ 2.17.0's frozen bytes (LOSAT/tests/fixtures/db_search). Each search logs its time and the engine's
// memory after it; BLASTP's and BLASTX's rows run twice more on one thread at the end, when the browser has compiled
// their code, for the time of a search without the compilation (line WQA D12a's measurement, P4).
import { createHash } from 'node:crypto';
import { readdirSync, readFileSync, statSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { expect } from '@playwright/test';
import type { ProgramId } from 'losat/domain';
import type { DatabaseResult, DatabaseSearch } from './harness/engine';
import { buildHarness } from './support/browser';
import { ENGINE, NO_ENGINE_REASON } from './support/engine';
import { openHarness, REPOSITORY, startHarnessServer, type HarnessServer } from './support/harness-server';
import { test } from './support/test';

test.skip(ENGINE === undefined, NO_ENGINE_REASON);
test.setTimeout(600_000);

const LOSAT = join(REPOSITORY, 'LOSAT');
const ROWS = ['db.blastn.mega.6', 'db.d1.6', 'db.blastx.default.6', 'db.tblastn.default.6', 'db.tblastx.default.6', 'db.sl.bp.bsl.6', 'db.tf.bp.nosql.6'];

interface Row {
  readonly id: string;
  readonly program: ProgramId;
  readonly db: string;
  readonly query: string;
  readonly extra: readonly string[];
}

function manifestRows(): Row[] {
  const rows = new Map<string, Row>();
  const directory = join(LOSAT, 'tests/db_search_manifest');
  for (const name of readdirSync(directory).filter((file) => file.endsWith('.tsv'))) {
    const lines = readFileSync(join(directory, name), 'utf8')
      .split('\n')
      .filter((line) => line !== '' && !line.startsWith('#'));
    const header = lines[0]!.split('\t');
    for (const line of lines.slice(1)) {
      const cell = Object.fromEntries(line.split('\t').map((value, i) => [header[i], value])) as Record<string, string>;
      const extra = cell['extra_args'] === '' ? [] : cell['extra_args']!.split(' ');
      rows.set(cell['row_id']!, { id: cell['row_id']!, program: cell['program'] as ProgramId, db: cell['db']!, query: cell['query']!, extra });
    }
  }
  return ROWS.map((id) => rows.get(id)!);
}

/** The files of a folder of LOSAT/ (not its subfolders), by their paths under LOSAT/. */
function folder(relative: string): string[] {
  return readdirSync(join(LOSAT, relative))
    .filter((name) => statSync(join(LOSAT, relative, name)).isFile())
    .map((name) => `${relative}/${name}`);
}

function searchOf(row: Row, threads: number): DatabaseSearch {
  const list = row.extra.indexOf('-seqidlist');
  const files = [...folder(dirname(row.db)), ...(list < 0 ? [] : [row.extra[list + 1]!])];
  return { id: row.id, program: row.program, argv: [row.program, '-db', row.db, '-query', row.query, ...row.extra], query: row.query, files, threads };
}

let server: HarnessServer;

test.beforeAll(async () => {
  server = await startHarnessServer(await buildHarness());
});

test.afterAll(async () => {
  await server?.close();
});

test('-db searches read the database in place and write NCBI frozen outfmt 6, on one and on four threads', async ({ context, browserName }) => {
  const page = await openHarness(context, server);
  const rows = manifestRows();
  const again = rows.filter((row) => row.id === 'db.d1.6' || row.id === 'db.blastx.default.6');
  const list = [
    ...[1, 4].flatMap((threads) => rows.map((row) => searchOf(row, threads))),
    ...[...again, ...again].map((row) => searchOf(row, 1)),
  ];
  const results: DatabaseResult[] = await page.evaluate((searches) => window.losatHarness!.engine.databaseSearches(searches), list);
  for (const [i, result] of results.entries()) {
    const label = `${browserName} ${result.id} n=${result.threads}`;
    expect.soft(result.error, label).toBeUndefined();
    const frozen = createHash('sha256').update(readFileSync(join(LOSAT, 'tests/fixtures/db_search', `${result.id}.out`))).digest('hex');
    expect.soft(result.sha256, `${label}: outfmt 6 = NCBI frozen`).toBe(frozen);
    expect.soft(result.runtime?.path, label).toBe(list[i]!.threads === 1 ? 'serial' : 'threaded');
    expect.soft(result.runtime?.threads, label).toBe(list[i]!.threads);
    console.log(
      `${label}: run ${Math.round(result.runMs ?? 0)} ms, search ${Math.round(result.totalMs)} ms, ` +
        `engine memory after ${((result.runtime?.memory?.linearBytesAfter ?? 0) / 2 ** 20).toFixed(1)} MiB`,
    );
  }
});
