import { fileURLToPath } from 'node:url';
import { defineConfig, devices, webkit } from '@playwright/test';
import { E2E_ORIGIN, E2E_PORT } from './tests/e2e/support/origin';

// E2E tests run against the production build served by `vite preview`, so they see
// the same response headers as the deployed site (public/_headers). Chromium, Firefox and
// WebKit run the tests; the storage tests that need the Chrome DevTools Protocol skip the
// other two. LOSAT_WEB_WEBKIT_EXECUTABLE can name a launcher for WebKit on a host whose
// system libraries Playwright cannot install (docs/evidence/losat_web_w1/README.md).
const webkitExecutable = process.env.LOSAT_WEB_WEBKIT_EXECUTABLE || undefined;
// WebKit paints on the CPU. With GPU painting, the web process of a closed page died of
// SIGSEGV on a SkiaGPUWorker thread a few times a run, after its test had passed: an
// instruction fetch from an unmapped page where Mesa's libgallium had been, apparently as
// the process exits. WSL saves a dump of several GB of each. No test compares screenshots;
// the Dot Plot tests count canvas pixels near the HSP colours, and they pass on the CPU.
// WebKit's system memory monitor (UI process) is off, so a test does not depend on how much
// memory the rest of the machine uses. With /proc/meminfo faked to 97 % used for the UI
// process, 3 of 18 runs of engine.spec.ts in Playwright's WebKit 2311 hung: the web process
// answered critical pressure by deleting all JIT code, which waits for the VM's queued JIT
// plans. In the 2 hangs where JavaScriptCore's JIT worklist was read, no JIT thread was alive
// to compile them (JavaScriptCore still counted 3 active), so the page's main thread waited
// forever and a search never finished.
// On Linux, WebKit starts through tests/e2e/support/webkit-launcher.sh, which sets its core file limit to
// 1 byte: a crash is then reported at once, instead of after the 20 to 60 s that the kernel spent piping the
// core of a crashed web process to systemd-coredump on CI (WSL pipes it to wsl-capture-crash); see the
// launcher. It starts LOSAT_WEB_WEBKIT_EXECUTABLE when that is set, else Playwright's own WebKit.
const webkitReal = process.platform === 'linux' ? (webkitExecutable ?? playwrightWebkit()) : undefined;
const webkitLauncher =
  webkitReal === undefined ? undefined : fileURLToPath(new URL('./tests/e2e/support/webkit-launcher.sh', import.meta.url));
const webkitEnv = {
  ...process.env,
  WEBKIT_SKIA_ENABLE_CPU_RENDERING: '1',
  WEBKIT_DISABLE_MEMORY_PRESSURE_MONITOR: '1',
  ...(webkitReal === undefined ? {} : { LOSAT_WEB_WEBKIT_REAL: webkitReal }),
};
const webkitLaunch = webkitLauncher ?? webkitExecutable;

/** Playwright's WebKit launcher, or undefined on a host for which Playwright has no WebKit build. */
function playwrightWebkit(): string | undefined {
  try {
    return webkit.executablePath();
  } catch {
    return undefined;
  }
}
// The WebKit pages' Screen Wake Lock is a stand-in: tests/e2e/support/test.ts.
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
        launchOptions: {
          env: webkitEnv,
          ...(webkitLaunch ? { executablePath: webkitLaunch } : {}),
        },
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
