// Navigations of the E2E tests. WebKit now and then never receives the response of a
// navigation's document request: the request goes out, the server is answering (the other
// worker's pages keep loading, and the harness server answers each request at once), and the
// page stays blank until the test times out. In CI this was 21 page.goto timeouts of WebKit in
// the 7 days to 2026-11-03, with `vite preview` and with the harness server alike. Locally it
// was 4 in 4,240 navigations of cross-origin-isolated pages (COOP same-origin, as every page
// of the site is) and none in 3,600 of a page without COOP, too few to name COOP as the cause.
// A second navigation of the same page got its response within a second each time, so each
// navigation waits a bounded time and is made once more, and the test records that it was
// (WFL, losat-baselines/web-flaky-20261103).
import { errors, type Page } from '@playwright/test';
import { test } from './test';

/** How long a first navigation may take before it is made again (normal ones take < 1 s). */
const FIRST_WAIT_MS = 15_000;

type Wait = { timeout?: number };

async function twice(page: Page, what: string, go: (wait: Wait) => Promise<unknown>): Promise<void> {
  try {
    await go({ timeout: FIRST_WAIT_MS });
  } catch (error) {
    if (!(error instanceof errors.TimeoutError) || page.isClosed()) throw error;
    test.info().annotations.push({ type: 'navigation made again', description: what });
    console.warn(`navigation made again after ${FIRST_WAIT_MS} ms: ${what}`);
    // The second one waits as long as the test lets it (the configured navigation timeout).
    await go({});
  }
}

/** `page.goto(url)`, made once more if its response never comes. */
export async function navigate(page: Page, url: string): Promise<void> {
  await twice(page, url, (wait) => page.goto(url, wait));
}

/** `page.reload()`, made once more if its response never comes. */
export async function reload(page: Page): Promise<void> {
  await twice(page, `reload ${page.url()}`, (wait) => page.reload(wait));
}
