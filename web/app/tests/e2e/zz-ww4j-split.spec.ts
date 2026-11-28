// WW4j throwaway probe (never committed; WW4j step 3, "serial first"): the first check and the first run of the
// WindowMasker tab, with the Gateway's two engine modules loaded as today (`parallel`: both fetched, hashed and compiled at
// once) and with the threaded one loaded only after the serial one (`serial-first`, a switch of the instrumentation, set
// by an init script): there `gw-serial-compiled` is the serial module's time alone, when a Gateway that started the
// Engine worker on the serial module could start it. The two arms alternate in one run, under the same load.
import { expect } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { described, openSearch, waitStatus } from './support/search';

test.skip(!BUILD_HAS_ENGINE, 'the real engine only');
test.setTimeout(300_000);

for (const arm of ['parallel', 'serial-first'] as const) {
  test(`ww4j split ${arm}: the first check and the first run of the WindowMasker tab`, async ({ page }) => {
    if (arm === 'serial-first') {
      await page.addInitScript(() => {
        (globalThis as { __ww4jSerialFirst?: boolean }).__ww4jSerialFirst = true;
      });
    }
    await openSearch(page);
    await page.getByTestId('tab-windowmasker').click();
    await described(page, 'wm-panel', 'ready');
    // The run is asked at once: it waits for the check of the form, which waits for the Engine worker's start.
    await page.getByTestId('wm-fasta-files').setInputFiles({
      name: 'a.fa',
      mimeType: 'application/octet-stream',
      buffer: Buffer.from('>a\nACGTACGTACGTACGTACGTAAACCCGGGTTTACGATCGATCGTAGCTAGCTAGCATCGACTAGCTACG\n'),
    });
    await expect(page.locator('[data-testid="wm-fasta-list"] > li').first()).toHaveAttribute('data-status', 'ready', { timeout: 60_000 });
    await page.getByTestId('wm-run-losat').click();
    await waitStatus(page, 1, 'completed', 240_000);
  });
}
