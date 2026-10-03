import { defineConfig, devices } from '@playwright/test';

// E2E tests run against the production build served by `vite preview`, so they see
// the same response headers as the deployed site (public/_headers). Chromium, Firefox and
// WebKit run the tests; the storage tests that need the Chrome DevTools Protocol skip the
// other two. LOSAT_WEB_WEBKIT_EXECUTABLE can name a launcher for WebKit on a host whose
// system libraries Playwright cannot install (docs/evidence/losat_web_w1/README.md).
const webkitExecutable = process.env.LOSAT_WEB_WEBKIT_EXECUTABLE;

export default defineConfig({
  testDir: 'tests/e2e',
  fullyParallel: true,
  // At most two browsers at a time: the same machine runs the engine gates (plan DW-7).
  workers: 2,
  forbidOnly: !!process.env.CI,
  retries: 0,
  reporter: process.env.CI ? 'github' : 'list',
  use: {
    baseURL: 'http://localhost:4173',
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
    command: 'npm run build && npm run preview',
    url: 'http://localhost:4173',
    reuseExistingServer: !process.env.CI,
    timeout: 180_000,
  },
});
