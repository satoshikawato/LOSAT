// A removed source gives its memory back (SURVEY #15): the Data worker used to keep every
// source and record table until the tab closed, so adding and removing a large input again
// and again made the tab grow without end. The test opens a FASTA of 100,000 records (32 MB)
// and removes it ten times, and compares the JavaScript heap of the page and of the Data worker
// (Chrome DevTools Protocol, after a garbage collection) after the first round and after the
// last: the growth must stay below what one source takes. The same file opens again each round,
// as a removed source is brought back. Chromium only (the protocol).
import { writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, test, type Page } from '@playwright/test';
import { openFiles, openSearch, program } from './support/search';

test.skip(({ browserName }) => browserName !== 'chromium', 'the heap is read with the Chrome DevTools Protocol');
test.setTimeout(600_000);

const RECORDS = 100_000;
const ROUNDS = 10;
const EVIDENCE = process.env['LOSAT_WEB_EVIDENCE'] || undefined;

/** RECORDS records of 300 nucleotides (a fixed pseudo-random sequence). */
function largeFasta(): Buffer {
  const letters = 'ACGT';
  let state = 7;
  const lines: string[] = [];
  for (let i = 0; i < RECORDS; i++) {
    let sequence = '';
    for (let j = 0; j < 300; j++) {
      state = (state * 1_103_515_245 + 12_345) % 2_147_483_648;
      sequence += letters[state >>> 29]!;
    }
    lines.push(`>r${i + 1} record ${i + 1}\n${sequence}\n`);
  }
  return Buffer.from(lines.join(''));
}

interface Memory {
  /** The JavaScript heap of the page and of the Data worker, after a garbage collection. */
  readonly total: number;
  readonly dataWorker: number;
}

/** Collects the garbage and reads the heap in use of a target (Chrome DevTools Protocol). */
type Call = (method: string, params?: object) => Promise<unknown>;

async function heapUsed(call: Call): Promise<number> {
  await call('HeapProfiler.collectGarbage');
  return ((await call('Runtime.getHeapUsage')) as { usedSize: number }).usedSize;
}

/**
 * The page's heap and the Data worker's, which a browser session reaches by attaching to the
 * worker's target and passing it messages (Chromium's `measureUserAgentSpecificMemory` is not
 * available in its headless shell).
 */
async function memoryProbe(page: Page): Promise<() => Promise<Memory>> {
  const pageSession = await page.context().newCDPSession(page);
  const browserSession = await page.context().browser()!.newBrowserCDPSession();
  const { targetInfos } = (await browserSession.send('Target.getTargets')) as {
    targetInfos: Array<{ targetId: string; type: string; url: string }>;
  };
  const worker = targetInfos.find((target) => target.type === 'worker' && /data-worker/.test(target.url));
  if (worker === undefined) throw new Error(`no Data worker among ${targetInfos.map((t) => `${t.type} ${t.url}`).join(', ')}`);
  const { sessionId } = (await browserSession.send('Target.attachToTarget', { targetId: worker.targetId, flatten: false })) as {
    sessionId: string;
  };
  let nextId = 1;
  const waiting = new Map<number, (message: { result?: unknown; error?: { message: string } }) => void>();
  browserSession.on('Target.receivedMessageFromTarget', (event: { sessionId: string; message: string }) => {
    if (event.sessionId !== sessionId) return;
    const message = JSON.parse(event.message) as { id?: number; result?: unknown; error?: { message: string } };
    if (message.id !== undefined) waiting.get(message.id)?.(message);
  });
  const workerCall: Call = (method, params = {}) =>
    new Promise((resolve, reject) => {
      const id = nextId++;
      waiting.set(id, (message) => {
        waiting.delete(id);
        if (message.error !== undefined) reject(new Error(`${method}: ${message.error.message}`));
        else resolve(message.result);
      });
      void browserSession.send('Target.sendMessageToTarget', { sessionId, message: JSON.stringify({ id, method, params }) });
    });
  const pageCall: Call = (method, params) => pageSession.send(method as 'Runtime.getHeapUsage', params as never);
  return async () => {
    const dataWorker = await heapUsed(workerCall);
    return { total: (await heapUsed(pageCall)) + dataWorker, dataWorker };
  };
}

const mb = (bytes: number) => Math.round(bytes / 1e5) / 10;

test('adding and removing a large source ten times does not make the tab grow', async ({ page }) => {
  await openSearch(page);
  expect(await page.evaluate(() => crossOriginIsolated)).toBe(true);
  await program(page, 'blastn');
  const fasta = largeFasta();
  const memory = await memoryProbe(page);
  const rounds: Array<{ round: number; held: Memory; removed: Memory }> = [];
  for (let round = 1; round <= ROUNDS; round++) {
    await openFiles(page, 'query', [{ name: 'large.fa', text: fasta }]);
    await expect(page.getByTestId('query-source-0-summary')).toHaveText('100,000 records · 30,000,000 nt');
    await expect(page.getByTestId('query-source-0-check')).toHaveAttribute('data-check', 'ok', { timeout: 120_000 });
    const held = await memory();
    await page.getByTestId('query-source-0-remove').click();
    await expect(page.locator('[data-testid^="query-source-"][data-status]')).toHaveCount(0);
    rounds.push({ round, held, removed: await memory() });
  }
  const first = rounds[0]!;
  const last = rounds.at(-1)!;
  const table = first.held.total - first.removed.total;
  const growth = last.removed.total - first.removed.total;
  const summary = {
    records: RECORDS,
    bytes: fasta.length,
    rounds: rounds.map((r) => ({
      round: r.round,
      heldMB: mb(r.held.total),
      removedMB: mb(r.removed.total),
      dataWorkerHeldMB: mb(r.held.dataWorker),
      dataWorkerRemovedMB: mb(r.removed.dataWorker),
    })),
    oneSourceMB: mb(table),
    growthMB: mb(growth),
    dataWorkerGrowthMB: mb(last.removed.dataWorker - first.removed.dataWorker),
  };
  console.log(JSON.stringify(summary));
  if (EVIDENCE !== undefined) writeFileSync(join(EVIDENCE, 'source-memory-chromium.json'), `${JSON.stringify(summary, null, 2)}\n`);
  expect(table).toBeGreaterThan(5e6);
  expect(growth).toBeLessThan(table);
});
