import { readFileSync } from 'node:fs';
import { defineConfig } from 'vitest/config';
import vue from '@vitejs/plugin-vue';
import type { Plugin } from 'vite';
import { siteHeaders } from './build/headers.ts';
import { findReactors, losatEngine } from './build/reactors.ts';
import { appBuild, productionBuildProblem, releaseId, sourceState } from './build/release.ts';
import { losatSite } from './build/site.ts';
import { losatVerification } from './build/verification.ts';

const headers = siteHeaders();
/** The app's version and the commit it was built from, written into session files (src/build-info.d.ts). */
const appVersion = (JSON.parse(readFileSync(new URL('./package.json', import.meta.url), 'utf8')) as { version: string }).version;
const reactors = findReactors();
/** The commit and whether a build input differs from it (build/release.ts, build/inputs.ts). */
const source = sourceState();
// The dev server injects component styles as inline <style> elements, which the
// production CSP forbids. It therefore serves the isolation headers without the CSP;
// `vite preview` and the site server of the E2E tests (build/site-server.ts) serve the production headers unchanged.
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

/**
 * Stops `vite build` unless Vite resolved a production build (build/release.ts
 * productionBuildProblem). Checked once the config is resolved, because Vite reads NODE_ENV from
 * `.env` files and VITE_USER_NODE_ENV after the config function has run; the build has written
 * nothing yet.
 */
function productionBuildOnly(): Plugin {
  return {
    name: 'losat-production-build-only',
    configResolved(config) {
      const problem = productionBuildProblem(config);
      if (problem !== undefined) throw new Error(problem);
    },
  };
}

export default defineConfig(({ command }) => {
  /**
   * The release ID (build/release.ts): the directory of the build's assets and the ID a browser
   * pins. A build without a commit must be named by LOSAT_WEB_RELEASE_ID; the development
   * server and Vitest (`serve`) run without one.
   */
  const release = releaseId(reactors, { source, build: command === 'build' });
  const build = appBuild(source);
  return {
    // The engine modules (LOSAT_WEB_REACTORS); the worker bundles read the same description.
    plugins: [vue(), losatEngine({ reactors }), losatVerification(), headersOnEveryResponse(),
      // `vite build` writes a site: the boot page and the release under /r/<release-id>/ (build/site.ts).
      losatSite({ release, app: { version: appVersion, build }, reactors }),
      // A release is a production build only, however Vite decided the mode.
      productionBuildOnly()],
    define: {
      __LOSAT_TEST_HOOKS__: 'false',
      __LOSAT_APP_VERSION__: JSON.stringify(appVersion),
      __LOSAT_APP_BUILD__: JSON.stringify(build),
      __LOSAT_RELEASE_ID__: JSON.stringify(release),
    },
    server: { headers: devHeaders },
    preview: { headers },
    worker: { format: 'es', plugins: () => [losatEngine({ reactors, emit: false })] },
    build: { target: 'es2022' },
    test: {
      include: ['tests/unit/**/*.test.ts'],
      environment: 'node',
    },
  };
});
