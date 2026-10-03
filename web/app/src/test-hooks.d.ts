// True only in the contract harness build (tests/e2e/support/browser.ts), which lets the
// E2E tests drive the Engine worker's run output writer. The application build defines it
// as false (vite.config.ts), so the hooks are removed from it.
declare const __LOSAT_TEST_HOOKS__: boolean;
