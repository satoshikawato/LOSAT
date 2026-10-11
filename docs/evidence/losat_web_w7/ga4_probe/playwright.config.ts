// The GA4 probe (RESULTS.md). Run it with run.sh, which uses web/app's node_modules.
// One worker: the browsers run one after another against the same servers.
import { defineConfig, devices } from '@playwright/test';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { A, C, G, M } from './origins';

const webkitExecutable = process.env['LOSAT_WEB_WEBKIT_EXECUTABLE'];

/**
 * Playwright launches Chromium with ThirdPartyStoragePartitioning disabled (its
 * chromiumSwitches.ts); stock Chrome partitions third-party storage since version 115. The
 * `chromium-stock-storage` project replaces Playwright's --disable-features argument with
 * the same list less that feature. The list is read from the installed playwright-core.
 */
function playwrightDisabledFeatures(): string[] {
  const bundle = readFileSync(createRequire(__filename).resolve('playwright-core/lib/coreBundle'), 'utf8');
  const start = bundle.indexOf('disabledFeatures = [');
  const end = bundle.indexOf('];', start);
  if (start < 0 || end < 0) throw new Error("playwright-core's disabledFeatures list was not found");
  const list = bundle
    .slice(start, end)
    .split('\n')
    .filter((line) => !line.trim().startsWith('//'))
    .join('\n');
  return [...list.matchAll(/"([^"]+)"/g)].map((match) => match[1] as string);
}
const disabled = playwrightDisabledFeatures();

export default defineConfig({
  testDir: '.',
  testMatch: 'probe.spec.ts',
  outputDir: process.env['PROBE_OUTPUT_DIR'] || 'test-results',
  workers: 1,
  retries: 0,
  timeout: 120_000,
  reporter: 'list',
  projects: [
    { name: 'chromium', use: { ...devices['Desktop Chrome'] } },
    { name: 'firefox', use: { ...devices['Desktop Firefox'] } },
    {
      name: 'webkit',
      use: {
        ...devices['Desktop Safari'],
        ...(webkitExecutable ? { launchOptions: { executablePath: webkitExecutable } } : {}),
      },
    },
    // The browsers' stock third-party storage rules, where Playwright changes them.
    {
      name: 'chromium-stock-storage',
      use: {
        ...devices['Desktop Chrome'],
        launchOptions: {
          ignoreDefaultArgs: [`--disable-features=${disabled.join(',')}`],
          args: [`--disable-features=${disabled.filter((name) => name !== 'ThirdPartyStoragePartitioning').join(',')}`],
        },
      },
    },
    {
      // Playwright's Firefox sets network.cookie.cookieBehavior 4 (playwright.cfg); stock
      // Firefox uses 5, Total Cookie Protection (third-party cookies and storage partitioned).
      name: 'firefox-stock-storage',
      use: { ...devices['Desktop Firefox'], launchOptions: { firefoxUserPrefs: { 'network.cookie.cookieBehavior': 5 } } },
    },
  ],
  webServer: {
    command: 'node server.mjs',
    url: `${A}/__ready`,
    reuseExistingServer: false,
    env: { PROBE_A: A, PROBE_M: M, PROBE_G: G, PROBE_C: C },
  },
});
