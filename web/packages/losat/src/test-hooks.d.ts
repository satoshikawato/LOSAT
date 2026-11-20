// True only in LOSAT Web's contract harness build (web/app/tests/e2e/support/browser.ts), which
// lets its E2E tests drive the Engine worker's run output writer. A build that bundles these modules
// defines it as false (web/app/vite.config.ts), so the hooks are removed from it.
declare const __LOSAT_TEST_HOOKS__: boolean;
