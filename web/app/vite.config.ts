import { defineConfig } from 'vitest/config';
import vue from '@vitejs/plugin-vue';
import type { Plugin } from 'vite';
import { HEADERS_FILE, parseHeadersFile, pathHeaders, siteHeaders } from './build/headers.ts';
import { readFileSync } from 'node:fs';
import { findReactors, losatEngine } from './build/reactors.ts';
import { losatVerification } from './build/verification.ts';

const headers = siteHeaders();
const rules = parseHeadersFile(readFileSync(HEADERS_FILE, 'utf8'));
const reactors = findReactors();
// The dev server injects component styles as inline <style> elements, which the
// production CSP forbids. It therefore serves the isolation headers without the CSP;
// `vite preview` (used by the E2E tests) serves the production headers unchanged.
const withoutCsp = (values: Record<string, string>) =>
  Object.fromEntries(Object.entries(values).filter(([name]) => name !== 'Content-Security-Policy'));
const devHeaders = withoutCsp(headers);

/**
 * The headers on every response of a server, 304 Not Modified included, with the rules of
 * `_headers` that match the path (the engine modules' caching). Vite's `headers` options do
 * not reach a 304, and WebKit judges a revalidated worker script by the headers of the 304:
 * without Cross-Origin-Embedder-Policy and Cross-Origin-Resource-Policy there, it refuses
 * the Engine worker that a cancel starts again (S12).
 */
function headersOnEveryResponse(): Plugin {
  const set = (res: { setHeader(name: string, value: string): unknown }, values: Record<string, string>) => {
    for (const [name, value] of Object.entries(values)) res.setHeader(name, value);
  };
  const of = (url: string | undefined) => pathHeaders(new URL(url ?? '/', 'http://localhost').pathname, rules);
  return {
    name: 'losat-headers-on-every-response',
    configureServer(server) {
      server.middlewares.use((req, res, next) => {
        set(res, withoutCsp(of(req.url)));
        next();
      });
    },
    configurePreviewServer(server) {
      server.middlewares.use((req, res, next) => {
        set(res, of(req.url));
        next();
      });
    },
  };
}

export default defineConfig({
  // The engine modules (LOSAT_WEB_REACTORS); the worker bundles read the same description.
  plugins: [vue(), losatEngine({ reactors }), losatVerification(), headersOnEveryResponse()],
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
