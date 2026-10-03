import { defineConfig } from 'vitest/config';
import vue from '@vitejs/plugin-vue';
import { siteHeaders } from './build/headers.ts';
import { findReactors, losatEngine } from './build/reactors.ts';

const headers = siteHeaders();
const reactors = findReactors();
// The dev server injects component styles as inline <style> elements, which the
// production CSP forbids. It therefore serves the isolation headers without the CSP;
// `vite preview` (used by the E2E tests) serves the production headers unchanged.
const devHeaders = Object.fromEntries(
  Object.entries(headers).filter(([name]) => name !== 'Content-Security-Policy'),
);

export default defineConfig({
  // The engine modules (LOSAT_WEB_REACTORS); the worker bundles read the same description.
  plugins: [vue(), losatEngine({ reactors })],
  define: { __LOSAT_TEST_HOOKS__: 'false' },
  server: { headers: devHeaders },
  preview: { headers },
  worker: { format: 'es', plugins: () => [losatEngine({ reactors, emit: false })] },
  build: { target: 'es2022' },
  test: {
    include: ['tests/unit/**/*.test.ts'],
    environment: 'node',
  },
});
