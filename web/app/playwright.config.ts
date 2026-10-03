import { defineConfig, devices } from '@playwright/test';
import { E2E_ORIGIN, E2E_PORT } from './tests/e2e/support/origin';

// E2E tests run against the production build served by `vite preview`, so they see
// the same response headers as the deployed site (public/_headers). Chromium, Firefox and
// WebKit run the tests; the storage tests that need the Chrome DevTools Protocol skip the
// other two. LOSAT_WEB_WEBKIT_EXECUTABLE can name a launcher for WebKit on a host whose
// system libraries Playwright cannot install (docs/evidence/losat_web_w1/README.md).
const webkitExecutable = process.env.LOSAT_WEB_WEBKIT_EXECUTABLE;
// The engine track runs the same tests in its own worktree on the same machine (plan DW-7):
// a server already on the port is never reused, and LOSAT_WEB_E2E_PORT can move this one.

export default defineConfig({
  testDir: 'tests/e2e',
  fullyParallel: true,
  // At most two browsers at a time: the same machine runs the engine gates (plan DW-7).
  workers: 2,
  forbidOnly: !!process.env.CI,
  retries: 0,
  reporter: process.env.CI ? 'github' : 'list',
  use: {
    baseURL: E2E_ORIGIN,
    trace: 'retain-on-failure',
  },
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
  ],
  webServer: {
    command: `npm run build && npx vite preview --port ${E2E_PORT} --strictPort`,
    url: E2E_ORIGIN,
    reuseExistingServer: false,
    timeout: 180_000,
  },
});
