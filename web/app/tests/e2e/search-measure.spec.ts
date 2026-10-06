// Measurements of the search screen's inputs (S12; plan §2.3 "レコード表の保存", S10's open
// item: the time to index and hash many records). Not a test: it runs only with
// LOSAT_WEB_MEASURE containing `records` (or `all`) and writes its records to
// LOSAT_WEB_EVIDENCE. docs/evidence/losat_web_w3/README.md reports the results.
//
// Through the application, for a protein FASTA of N records (300 residues each, lines of
// 60): the time from choosing the file to its record table (the Data worker's index scan
// with the engine's serial reactor and the SHA-256 of every record), the engine's check of
// the input (`register` on the Data worker's reactor), leaving one record out (a new
// revision and the check again), and "Add to queue" (the run input and its SHA-256, the
// engine's validate). One warm-up and three measured repetitions, each in a new page.
import { writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, test, type Page } from '@playwright/test';
import { ENGINE, NO_ENGINE_REASON } from './support/engine';

const MEASURE = new Set((process.env['LOSAT_WEB_MEASURE'] ?? '').split(',').filter((name) => name !== ''));
const EVIDENCE = process.env['LOSAT_WEB_EVIDENCE'] || undefined;
const COUNTS = (process.env['LOSAT_WEB_MEASURE_RECORDS'] ?? '1000,10000,100000').split(',').map(Number);

test.skip(!MEASURE.has('records') && !MEASURE.has('all'), 'measurements run only with LOSAT_WEB_MEASURE');
test.skip(ENGINE === undefined, NO_ENGINE_REASON);
test.setTimeout(3_600_000);

const LETTERS = 'ACDEFGHIKLMNPQRSTVWY';

/** N protein records of 300 residues, deterministic. */
function proteins(count: number): Buffer {
  const parts: string[] = [];
  let state = 12345;
  for (let i = 0; i < count; i++) {
    let sequence = '';
    for (let j = 0; j < 300; j++) {
      state = (Math.imul(state, 1103515245) + 12345) >>> 0;
      sequence += LETTERS[(state >>> 16) % LETTERS.length];
    }
    parts.push(`>prot_${i + 1} synthetic protein ${i + 1}\n${sequence.replace(/(.{60})/g, '$1\n').replace(/\n$/, '')}\n`);
  }
  return Buffer.from(parts.join(''));
}

/**
 * Milliseconds. `*InApp` are the application's own times (the draft's clock around the
 * Data worker's calls); the others are what the user waits, seen from the test, and the
 * exclusion includes the pause of 300 ms after which the draft checks the input again.
 */
interface Sample {
  readonly indexMs: number;
  readonly indexInApp: number;
  readonly checkInApp: number;
  readonly excludeMs: number;
  readonly recheckInApp: number;
  readonly enqueueMs: number;
}

async function once(page: Page, file: Buffer, count: number): Promise<Sample> {
  await page.goto('/');
  await page.getByTestId('program-blastp').check();
  await page.getByTestId('subject-input').fill('>s\nMKLVVLAAGGHHKL\n');
  await expect(page.getByTestId('subject-source-0-check')).toHaveAttribute('data-check', 'ok');
  const source = page.getByTestId('query-source-0');
  const check = page.getByTestId('query-source-0-check');

  const t0 = Date.now();
  await page.getByTestId('query-files').setInputFiles({ name: `proteins-${count}.faa`, mimeType: 'text/plain', buffer: file });
  await expect(source).toHaveAttribute('data-status', 'ready', { timeout: 600_000 });
  const t1 = Date.now();
  await expect(check).toHaveAttribute('data-check', 'ok', { timeout: 600_000 });
  const indexInApp = Number(await source.getAttribute('data-index-ms'));
  const checkInApp = Number(await source.getAttribute('data-check-ms'));

  await source.locator('details.records').evaluate((details) => ((details as HTMLDetailsElement).open = true));
  const t3 = Date.now();
  await page.getByTestId('query-source-0-record-0').uncheck();
  // The summary and the pending check change together; then wait for the engine's verdict.
  await expect(page.getByTestId('query-source-0-summary')).toContainText(`(${(count - 1).toLocaleString('en-US')} included)`);
  await expect(check).toHaveAttribute('data-check', 'ok', { timeout: 600_000 });
  const t4 = Date.now();
  const recheckInApp = Number(await source.getAttribute('data-check-ms'));

  await page.getByTestId('add-to-queue').click();
  await expect(page.getByTestId('run-1')).toBeVisible({ timeout: 600_000 });
  const t5 = Date.now();
  await page.getByTestId('run-1-cancel').click({ timeout: 2000 }).catch(() => undefined);
  return { indexMs: t1 - t0, indexInApp, checkInApp, excludeMs: t4 - t3, recheckInApp, enqueueMs: t5 - t4 };
}

const median = (values: readonly number[]) => [...values].sort((a, b) => a - b)[Math.floor(values.length / 2)]!;

test('index, check, exclusion and queueing of many records', async ({ page, browserName }) => {
  const results: unknown[] = [];
  for (const count of COUNTS) {
    const file = proteins(count);
    await once(page, file, count); // warm-up
    const samples: Sample[] = [];
    for (let i = 0; i < 3; i++) samples.push(await once(page, file, count));
    const summary = Object.fromEntries(
      (Object.keys(samples[0]!) as Array<keyof Sample>).map((key) => {
        const values = samples.map((sample) => sample[key]);
        return [key, { median: median(values), min: Math.min(...values), max: Math.max(...values) }];
      }),
    );
    results.push({ count, bytes: file.length, samples, summary });
    console.log(`${browserName} ${count} records (${file.length} bytes): ${JSON.stringify(summary)}`);
  }
  if (EVIDENCE !== undefined) {
    writeFileSync(join(EVIDENCE, `records-${browserName}.json`), `${JSON.stringify(results, null, 2)}\n`);
  }
});
