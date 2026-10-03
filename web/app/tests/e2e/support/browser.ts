// Helpers for the storage E2E tests: a same-origin probe page that reads OPFS and the Web
// Locks state without starting the application, quota control through the Chrome DevTools
// Protocol, and the in-memory build of the contract harness (tests/e2e/harness).
import { expect, type BrowserContext, type CDPSession, type Page } from '@playwright/test';
import { fileURLToPath } from 'node:url';
import { build } from 'vite';
import { siteHeaders } from '../../../build/headers';
import { findReactors, losatEngine } from '../../../build/reactors';
import { SESSION_LOCK_PREFIX, TMP_DIRECTORY } from '../../../src/infra/data/session';
import { E2E_ORIGIN } from './origin';

export const ORIGIN = E2E_ORIGIN;

/** Whether the application build has the engine (LOSAT_WEB_REACTORS), or uses the FakeEngine. */
const REACTORS = findReactors();
export const BUILD_HAS_ENGINE = REACTORS !== undefined;

/**
 * Text in the stored outfmt 7 of the BLASTN search of queueSearch: the engine's header, or
 * the FakeEngine's marker.
 */
export const OUTFMT7_MARK = BUILD_HAS_ENGINE ? '# BLASTN 2.17.0+' : 'FAKE ENGINE OUTPUT';

/** Opens a page of the site's origin that does not run the application. */
export async function openProbe(context: BrowserContext): Promise<Page> {
  await context.route('**/__probe', (route) =>
    route.fulfill({ status: 200, contentType: 'text/html', body: '<!doctype html><title>probe</title>' }),
  );
  const probe = await context.newPage();
  await probe.goto('/__probe');
  return probe;
}

/** The session directories under tmp/ in OPFS, sorted. */
export function sessionDirectories(page: Page): Promise<string[]> {
  return page.evaluate(async (tmpName) => {
    const root = await navigator.storage.getDirectory();
    let tmp: FileSystemDirectoryHandle;
    try {
      tmp = await root.getDirectoryHandle(tmpName);
    } catch {
      return [];
    }
    const names: string[] = [];
    const entries = (tmp as unknown as { entries(): AsyncIterable<[string, FileSystemHandle]> }).entries();
    for await (const [name, handle] of entries) if (handle.kind === 'directory') names.push(name);
    return names.sort();
  }, TMP_DIRECTORY);
}

/** Every file under tmp/<token>/, as sorted relative paths. */
export function sessionFiles(page: Page, token: string): Promise<string[]> {
  return page.evaluate(
    async ({ tmpName, token }) => {
      const root = await navigator.storage.getDirectory();
      const start = await (await root.getDirectoryHandle(tmpName)).getDirectoryHandle(token);
      const paths: string[] = [];
      const walk = async (directory: FileSystemDirectoryHandle, prefix: string): Promise<void> => {
        const entries = (directory as unknown as { entries(): AsyncIterable<[string, FileSystemHandle]> }).entries();
        for await (const [name, handle] of entries) {
          if (handle.kind === 'directory') await walk(handle as FileSystemDirectoryHandle, `${prefix}${name}/`);
          else paths.push(`${prefix}${name}`);
        }
      };
      await walk(start, '');
      return paths.sort();
    },
    { tmpName: TMP_DIRECTORY, token },
  );
}

/** The tokens of the session locks that are held in this origin, sorted. */
export function heldSessionLocks(page: Page): Promise<string[]> {
  return page.evaluate(async (prefix) => {
    const { held = [] } = await navigator.locks.query();
    return held
      .map((lock) => lock.name ?? '')
      .filter((name) => name.startsWith(prefix))
      .map((name) => name.slice(prefix.length))
      .sort();
  }, SESSION_LOCK_PREFIX);
}

const quotaSessions = new WeakMap<Page, CDPSession>();

/**
 * Sets the origin's storage quota in bytes (Chromium), or resets it with null. The
 * override lasts only while its DevTools session is attached, so the session is kept
 * until the reset. OPFS honours it only in a profile on disk (support/profile.ts).
 */
export async function setQuota(page: Page, bytes: number | null): Promise<void> {
  let cdp = quotaSessions.get(page);
  if (cdp === undefined) {
    cdp = await page.context().newCDPSession(page);
    quotaSessions.set(page, cdp);
  }
  if (bytes !== null) {
    await cdp.send('Storage.overrideQuotaForOrigin', { origin: ORIGIN, quotaSize: bytes });
    return;
  }
  await cdp.send('Storage.overrideQuotaForOrigin', { origin: ORIGIN });
  quotaSessions.delete(page);
  await cdp.detach();
}

/** Opens the application and waits until its start-up cleanup has finished. */
export async function openApp(context: BrowserContext): Promise<Page> {
  const page = await context.newPage();
  await page.goto('/');
  await expect(page.getByTestId('storage-status')).toHaveAttribute('data-cleanup', 'done');
  return page;
}

/** Pastes a query and a subject and queues a BLASTN search. */
export async function queueSearch(page: Page): Promise<void> {
  await page.getByTestId('tab-search').click();
  await page.getByTestId('query-input').fill('>q1\nACGTACGTACGT\n');
  await page.getByTestId('subject-input').fill('>s1\nACGTACGTACGT\n');
  await page.getByTestId('add-to-queue').click();
}

/** Queues a search and waits until it, run `number`, completes. */
export async function runSearch(page: Page, number: number): Promise<void> {
  await queueSearch(page);
  await expect(page.getByTestId(`run-${number}-status`)).toHaveText('completed');
}

/** Shows run `number` on the results tab and checks that its stored output can be read. */
export async function expectResult(page: Page, number: number): Promise<void> {
  await page.getByTestId('tab-results').click();
  await page.getByTestId('result-run').selectOption({ label: `Run ${number} · BLASTN` });
  await page.getByTestId('format-7').click();
  await expect(page.getByTestId('result-output')).toContainText(OUTFMT7_MARK);
}

export type HarnessFiles = ReadonlyMap<string, { readonly body: string | Uint8Array; readonly type: string }>;

const TYPES: Readonly<Record<string, string>> = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.map': 'application/json',
  '.wasm': 'application/wasm',
};

/**
 * Builds tests/e2e/harness in memory with Vite; the application build is not touched. The
 * harness build has the engine modules of LOSAT_WEB_REACTORS, as the application build does,
 * and the test hooks of the Engine worker (`__LOSAT_TEST_HOOKS__`).
 */
export async function buildHarness(): Promise<HarnessFiles> {
  const result = await build({
    root: fileURLToPath(new URL('../harness/', import.meta.url)),
    base: '/__harness/',
    configFile: false,
    publicDir: false,
    logLevel: 'warn',
    plugins: [losatEngine({ reactors: REACTORS })],
    define: { __LOSAT_TEST_HOOKS__: 'true' },
    worker: { format: 'es', plugins: () => [losatEngine({ reactors: REACTORS, emit: false })] },
    build: { write: false, target: 'es2022', minify: false, modulePreload: false },
  });
  const files = new Map<string, { body: string | Uint8Array; type: string }>();
  for (const output of Array.isArray(result) ? result : [result]) {
    if (!('output' in output)) throw new Error('the harness build did not return its output');
    for (const item of output.output) {
      const extension = item.fileName.slice(item.fileName.lastIndexOf('.'));
      files.set(item.fileName, {
        body: item.type === 'chunk' ? item.code : item.source,
        type: TYPES[extension] ?? 'application/octet-stream',
      });
    }
  }
  return files;
}

/** Serves the harness under /__harness/ with the site's response headers. */
export async function serveHarness(context: BrowserContext, files: HarnessFiles): Promise<void> {
  const headers = siteHeaders();
  await context.route('**/__harness/**', async (route) => {
    const path = new URL(route.request().url()).pathname.replace(/^\/__harness\//, '');
    const file = files.get(path === '' ? 'index.html' : path);
    if (file === undefined) {
      await route.fulfill({ status: 404, body: `not in the harness: ${path}` });
      return;
    }
    await route.fulfill({
      status: 200,
      headers: { ...headers, 'Content-Type': file.type },
      body: typeof file.body === 'string' ? file.body : Buffer.from(file.body),
    });
  });
}
