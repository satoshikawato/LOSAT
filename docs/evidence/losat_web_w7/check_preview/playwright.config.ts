// The preview check's configuration (README.md). Run it from web/app, whose node_modules it
// uses:
//   LOSAT_WEB_PREVIEW_URL=https://<preview> LOSAT_WEB_PREVIEW_OUT=<dir> \
//     npx playwright test -c ../../docs/evidence/losat_web_w7/check_preview/playwright.config.ts
// There is no web server: the site is already deployed (or served by build/site-server.ts for
// a dry run).
import { resolve } from 'node:path';
import { defineConfig, devices } from '@playwright/test';
import { previewOrigin } from './support';

const origin = previewOrigin();
const out = process.env['LOSAT_WEB_PREVIEW_OUT'];
if (out === undefined || out === '') throw new Error('LOSAT_WEB_PREVIEW_OUT is not set (README.md)');
const webkitExecutable = process.env['LOSAT_WEB_WEBKIT_EXECUTABLE'];

export default defineConfig({
  testDir: '.',
  outputDir: resolve(out, 'test-results'),
  // One test at a time: the check measures a deployed site, and a search in a browser needs the machine.
  workers: 1,
  fullyParallel: false,
  retries: 0,
  timeout: 300_000,
  reporter: [['list'], ['./record-reporter.ts', { out: resolve(out), origin }]],
  use: { baseURL: origin, trace: 'retain-on-failure' },
  projects: [
    { name: 'http', testMatch: 'http.spec.ts' },
    { name: 'chromium', testMatch: 'browser.spec.ts', use: { ...devices['Desktop Chrome'] } },
    { name: 'firefox', testMatch: 'browser.spec.ts', use: { ...devices['Desktop Firefox'] } },
    {
      name: 'webkit',
      testMatch: 'browser.spec.ts',
      use: {
        ...devices['Desktop Safari'],
        ...(webkitExecutable ? { launchOptions: { executablePath: webkitExecutable } } : {}),
      },
    },
  ],
});
