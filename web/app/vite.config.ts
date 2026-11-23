import { defineConfig } from 'vitest/config';
import vue from '@vitejs/plugin-vue';
import type { Plugin } from 'vite';
import { HEADERS_FILE, parseHeadersFile, pathHeaders, siteHeaders } from './build/headers.ts';
import { readFileSync } from 'node:fs';
import { findReactors, losatEngine } from './build/reactors.ts';
import { losatVerification } from './build/verification.ts';
import { losatNotices } from './build/notices.ts';

const headers = siteHeaders();
const rules = parseHeadersFile(readFileSync(HEADERS_FILE, 'utf8'));
const reactors = findReactors();
// The license and the notices, with the npm packages of the main and the worker bundles.
const notices = losatNotices();
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
  // The headers first, so that they reach the engine modules that losatEngine serves in dev.
  plugins: [headersOnEveryResponse(), vue(), losatEngine({ reactors }), losatVerification(), notices.site],
  define: { __LOSAT_TEST_HOOKS__: 'false' },
  server: { headers: devHeaders },
  preview: { headers },
  worker: { format: 'es', plugins: () => [losatEngine({ reactors, emit: false }), notices.worker()] },
  // No module preload polyfill: the site writes no modulepreload link (its one chunk besides the
  // main one, the -html text of domain/outfmt0-html.ts that application/results.ts imports when
  // it reads a range of such a run, imports nothing), so it would do nothing, and its code (adapted
  // from es-module-shims) has no notice in the site (build/notices.ts stops the build if it comes
  // back). Vite's preload helper, which calls that import(), is in the main chunk with the
  // notices of Vite and Rolldown (build/notices.ts, VIRTUAL_OWNERS).
  build: { target: 'es2022', modulePreload: { polyfill: false } },
  test: {
    include: ['tests/unit/**/*.test.ts'],
    environment: 'node',
  },
});
