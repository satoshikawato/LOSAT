// WW4j throwaway probe (never committed): the first check of a tool tab waits for the Engine worker's start. The two
// engine files that the Gateway fetches for it (not the Data worker's, which the describe has already loaded) arrive
// WW4J_DELAY_MS late, as they did under load; the check stays `checking` past Playwright's 5 s, and ends with its answer.
import { expect } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { described, openSearch } from './support/search';

const DELAY = Number(process.env['WW4J_DELAY_MS'] ?? '7000');

test.skip(!BUILD_HAS_ENGINE, 'the real engine only');
test.setTimeout(180_000);

test('ww4j probe: the first check of the WindowMasker tab waits for the engine files', async ({ page }) => {
  await openSearch(page);
  let fetched = 0;
  await page.route('**/engine/*.wasm', async (route) => {
    fetched++;
    const response = await route.fetch();
    await new Promise((resolve) => setTimeout(resolve, DELAY));
    await route.fulfill({ response });
  });
  await page.getByTestId('tab-windowmasker').click();
  await described(page, 'wm-panel', 'ready');
  const check = page.getByTestId('wm-argv-validation');
  await expect(check).toHaveAttribute('data-state', 'checking');
  const start = Date.now();
  let defaultWait = 'answered';
  try {
    await expect(check).toHaveAttribute('data-state', 'ok');
  } catch {
    defaultWait = `still ${await check.getAttribute('data-state')} after Playwright's 5 s`;
  }
  // The named wait of the fix, inlined so that the probe runs on dev's specs too.
  await expect(check).toHaveAttribute('data-state', 'ok', { timeout: 60_000 });
  const line = `[ww4j-probe] ${test.info().project.name}: engine files delayed ${DELAY} ms (${fetched} fetched); default wait: ${defaultWait}; ok ${Date.now() - start} ms after the expect began`;
  console.log(line);
  test.info().annotations.push({ type: 'ww4j', description: line });
  expect(fetched).toBe(2);
  expect(defaultWait).not.toBe('answered');
});
