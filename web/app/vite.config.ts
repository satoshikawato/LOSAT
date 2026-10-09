import { defineConfig } from 'vitest/config';
import vue from '@vitejs/plugin-vue';
import type { Plugin } from 'vite';
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

/**
 * The headers on every response of a server, 304 Not Modified included. Vite's `headers`
 * options do not reach a 304, and WebKit judges a revalidated worker script by the headers
 * of the 304: without Cross-Origin-Embedder-Policy and Cross-Origin-Resource-Policy there,
 * it refuses the Engine worker that a cancel starts again (S12).
 */
function headersOnEveryResponse(): Plugin {
  const set = (res: { setHeader(name: string, value: string): unknown }, values: Record<string, string>) => {
    for (const [name, value] of Object.entries(values)) res.setHeader(name, value);
  };
  return {
    name: 'losat-headers-on-every-response',
    configureServer(server) {
      server.middlewares.use((_req, res, next) => {
        set(res, devHeaders);
        next();
      });
    },
    configurePreviewServer(server) {
      server.middlewares.use((_req, res, next) => {
        set(res, headers);
        next();
      });
    },
  };
}

export default defineConfig({
  // The engine modules (LOSAT_WEB_REACTORS); the worker bundles read the same description.
  plugins: [vue(), losatEngine({ reactors }), headersOnEveryResponse()],
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
