// The GA4 probe (RESULTS.md): variant (a), a measurement document of another site in an
// iframe of the app, against variant (b), a fixed event posted to the app's own /api/visit,
// under the app's COOP/COEP/CORP and CSP, in Chromium, Firefox and WebKit. It records what
// happens; it asserts only what the probe itself needs (A is cross-origin isolated). Each
// browser's observations are written to results/<browser>.json.
import { expect, test, type BrowserContext, type Page } from '@playwright/test';
import { randomUUID } from 'node:crypto';
import { mkdirSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { A, G, M } from './origins';

type Report = Record<string, any>; // eslint-disable-line @typescript-eslint/no-explicit-any
interface LogEntry {
  server: 'A' | 'M' | 'G' | 'C';
  method: string;
  path: string;
  query: Record<string, string>;
  headers: Record<string, string | null>;
  body: string;
}

/** The third-party storage rules each project runs with (playwright.config.ts). */
const STORAGE_RULES: Record<string, string> = {
  chromium: "Playwright's: ThirdPartyStoragePartitioning disabled",
  'chromium-stock-storage': "stock: ThirdPartyStoragePartitioning enabled (Playwright's --disable-features less it)",
  firefox: "Playwright's: network.cookie.cookieBehavior 4 (reject trackers, no partitioning)",
  'firefox-stock-storage': 'stock: network.cookie.cookieBehavior 5 (Total Cookie Protection)',
  webkit: "Playwright's WPE MiniBrowser (Linux); not Safari's tracking prevention",
};

const rows: Report[] = [];
let browserVersion = '';

async function serverLog(run: string): Promise<LogEntry[]> {
  return (await (await fetch(`${G}/__log?run=${run}`)).json()) as LogEntry[];
}

async function probe(page: Page): Promise<Report> {
  await page.waitForFunction(() => (window as unknown as { __probe?: { done?: boolean } }).__probe?.done === true, null, { timeout: 30_000 });
  // Lets queued beacons reach C before the log is read.
  await page.waitForTimeout(500);
  return page.evaluate(() => (window as unknown as { __probe: Report }).__probe);
}

async function openA(page: Page, variant: string, csp: 'site' | 'frame', sandbox: string): Promise<{ run: string; report: Report }> {
  const run = randomUUID();
  // 'commit', not 'load': the probe waits for its own report, and a frame whose requests
  // are refused may hold the load event back.
  await page.goto(`${A}/a.html?variant=${variant}&csp=${csp}&sandbox=${sandbox}&run=${run}&m=${encodeURIComponent(M)}`, { waitUntil: 'commit' });
  const report = await probe(page);
  expect(report['crossOriginIsolated']).toBe(true);
  return { run, report };
}

/** What M's document reads of its own storage when it is opened top-level (first party). */
async function readTopLevel(context: BrowserContext): Promise<Report> {
  const page = await context.newPage();
  await page.goto(`${M}/m.html?variant=read&run=${randomUUID()}`, { waitUntil: 'commit' });
  const report = await probe(page);
  await page.close();
  return report;
}

function frameRow(report: Report, log: LogEntry[]): Report {
  const result = report['result'];
  const collect = log.filter((entry) => entry.server === 'C');
  return {
    aCspViolations: report['cspViolations'],
    mRequested: log.some((entry) => entry.server === 'M' && entry.path === '/m.html'),
    iframeLoadEvent: report['iframeLoadEvent'],
    frameLoaded: report['hello'] !== undefined,
    hello: report['hello'] ?? null,
    credentiallessAttributeSupported: report['credentiallessAttributeSupported'],
    resultReceived: result !== undefined,
    resultChecksPass: report['resultChecksPass'] ?? null,
    mRunMessages: result?.runMessages ?? null,
    mCrossOriginIsolated: result?.crossOriginIsolated ?? null,
    mCspViolations: result?.cspViolations ?? null,
    gScriptRequested: log.some((entry) => entry.server === 'G' && entry.path === '/gtag/js'),
    gScript: result?.script ?? null,
    transportsInPage: result?.transports ?? null,
    collectReached: [...new Set(collect.map((entry) => entry.query['transport']))].sort(),
    collectHeaders: collect.map((entry) => ({ transport: entry.query['transport'], ...entry.headers })),
    storage: result?.storage ?? null,
    hasStorageAccess: result?.hasStorageAccess ?? null,
  };
}

const FRAME_CASES: { id: string; variant: string; csp: 'site' | 'frame'; sandbox: string; persistence: boolean; what: string }[] = [
  { id: 'a1', variant: 'a1', csp: 'frame', sandbox: 'same-origin', persistence: true, what: 'M without COEP/CORP' },
  { id: 'a2', variant: 'a2', csp: 'frame', sandbox: 'same-origin', persistence: true, what: 'M with COEP require-corp and CORP cross-origin' },
  { id: 'a3', variant: 'a3', csp: 'frame', sandbox: 'same-origin', persistence: true, what: 'M with COEP credentialless and CORP cross-origin' },
  { id: 'a3n', variant: 'a3n', csp: 'frame', sandbox: 'same-origin', persistence: true, what: 'M with COEP credentialless, no CORP' },
  { id: 'a4', variant: 'a4', csp: 'frame', sandbox: 'same-origin', persistence: true, what: '<iframe credentialless>, M without COEP/CORP' },
  { id: 'a2-csp-unchanged', variant: 'a2', csp: 'site', sandbox: 'same-origin', persistence: false, what: "a2 under the site's CSP unchanged (no frame-src)" },
  { id: 'a2-no-sandbox', variant: 'a2', csp: 'frame', sandbox: 'none', persistence: false, what: 'a2 without sandbox' },
  { id: 'a2-sandbox-opaque', variant: 'a2', csp: 'frame', sandbox: 'opaque', persistence: false, what: 'a2 with sandbox="allow-scripts" (no allow-same-origin)' },
];

test.describe.serial('GA4 probe', () => {
  test.beforeAll(async ({ browser }) => {
    browserVersion = browser.version();
  });

  test.afterAll(async ({}, testInfo) => {
    const directory = join(testInfo.project.testDir, 'results');
    mkdirSync(directory, { recursive: true });
    const file = join(directory, `${testInfo.project.name}.json`);
    writeFileSync(
      file,
      `${JSON.stringify({ browser: testInfo.project.name, version: browserVersion, storageRules: STORAGE_RULES[testInfo.project.name], date: new Date().toISOString(), origins: { A, M, G }, rows }, null, 2)}\n`,
    );
  });

  test('b: a fixed event to the same-origin /api/visit, unchanged CSP', async ({ page }) => {
    const { run, report } = await openA(page, 'b', 'site', 'none');
    const log = await serverLog(run);
    const visits = log.filter((entry) => entry.server === 'A' && entry.path === '/api/visit');
    rows.push({
      id: 'b',
      what: 'POST /api/visit (same origin), site CSP unchanged',
      aCspViolations: report['cspViolations'],
      visit: report['visit'],
      reached: visits.length,
      body: visits[0]?.body ?? null,
      headers: visits[0]?.headers ?? null,
      otherRequests: log.filter((entry) => entry.server !== 'A' || !['/a.html', '/a.js', '/api/visit'].includes(entry.path)).length,
    });
  });

  test('navigation: A navigated again after a cross-origin frame (a2) or none (b) was in it', async ({ page }) => {
    const outcome: Report = {};
    for (const variant of ['b', 'a2']) {
      await openA(page, variant, 'frame', 'same-origin');
      const started = Date.now();
      try {
        await page.goto(`${A}/a.html?variant=b&csp=site&sandbox=none&run=${randomUUID()}`, { waitUntil: 'commit', timeout: 15_000 });
        outcome[variant] = `committed in ${Date.now() - started} ms`;
      } catch (error) {
        outcome[variant] = `not committed in 15 s (${String(error).split('\n')[0]})`;
        await page.goto('about:blank').catch(() => undefined);
      }
    }
    rows.push({ id: 'navigation', what: 'goto A again (commit) after A with variant b / a2', ...outcome });
  });

  for (const c of FRAME_CASES) {
    test(`${c.id}: ${c.what}`, async ({ page, context }) => {
      const first = await openA(page, c.variant, c.csp, c.sandbox);
      const row: Report = { id: c.id, what: c.what, aCsp: c.csp === 'frame' ? 'site + frame-src M' : 'site', sandbox: c.sandbox, ...frameRow(first.report, await serverLog(first.run)) };
      if (c.persistence && row['resultReceived']) {
        // The app opened again (a new tab after closing the first: Playwright's Firefox does
        // not commit a second navigation of this page, see the 'navigation' case), then a
        // second tab of it: does M's frame read the identifier it wrote? Then M top-level:
        // is the frame's storage partitioned from it?
        await page.waitForTimeout(300);
        await page.close();
        const again = await openA(await context.newPage(), c.variant, c.csp, c.sandbox);
        const reloaded = again.report['result']?.storage ?? null;
        const second = await context.newPage();
        const tab = await openA(second, c.variant, c.csp, c.sandbox);
        await second.close();
        const top = await readTopLevel(context);
        const written = row['storage'] ?? {};
        row['persistence'] = {
          afterReload: reloaded,
          secondTab: tab.report['result']?.storage ?? null,
          topLevel: top['storage'],
          topLevelHasStorageAccess: top['hasStorageAccess'],
          localKept: reloaded !== null && reloaded.localBefore !== null && reloaded.localBefore === written.localAfter,
          cookieNoneKept: reloaded !== null && reloaded.cookieBefore?.probe_none !== undefined && reloaded.cookieBefore?.probe_none === written.cookieAfter?.probe_none,
          cookieLaxKept: reloaded !== null && reloaded.cookieBefore?.probe_lax !== undefined && reloaded.cookieBefore?.probe_lax === written.cookieAfter?.probe_lax,
          localSameAsTopLevel: written.localAfter != null && top['storage']?.localBefore === written.localAfter,
        };
      }
      rows.push(row);
    });
  }
});
