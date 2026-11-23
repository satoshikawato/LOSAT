import { defineConfig } from 'vitest/config';

export default defineConfig({
  test: {
    include: ['tests/**/*.test.ts'],
    environment: 'node',
    // The built-package tests pack the package with npm, compile the README's example and run the
    // engine in child processes; their time follows the machine's load (10-33 s under a loaded
    // host, against vitest's 5 s default), so a slow host does not fail them.
    testTimeout: 120_000,
  },
});
