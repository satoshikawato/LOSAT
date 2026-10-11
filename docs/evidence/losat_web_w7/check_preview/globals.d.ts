// tests/e2e/support/harness-server.ts (imported through support/search.ts) reads a global that the
// contract harness declares; this check does not build the harness.
interface Window {
  losatHarness: unknown;
}
