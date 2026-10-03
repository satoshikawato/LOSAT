// An HTTP server for the contract harness (tests/e2e/harness) in Firefox and WebKit, and for
// the engine tests in every browser. Playwright's request routing does not reach every
// request of a worker in those browsers: in Firefox a worker started by a worker (the
// thread workers of the threaded engine) bypasses it, and in WebKit a worker whose script
// it fulfils is not cross-origin isolated (no SharedArrayBuffer). A real server sends the
// site's response headers on every response, as the deployed site does.
import { readFileSync } from 'node:fs';
import { createServer, type IncomingMessage, type ServerResponse } from 'node:http';
import type { AddressInfo } from 'node:net';
import type { BrowserContext, Page } from '@playwright/test';
import { fileURLToPath } from 'node:url';
import { siteHeaders } from '../../../build/headers';
import type { HarnessFiles } from './browser';

/** The repository root; /__files/<path> serves <root>/<path> (inputs of the searches). */
export const REPOSITORY = fileURLToPath(new URL('../../../../../', import.meta.url));

export interface HarnessServer {
  /** For example http://127.0.0.1:41234 */
  readonly origin: string;
  close(): Promise<void>;
}

export interface HarnessServerOptions {
  /** Without it, the server leaves out COOP and COEP: the page is not cross-origin isolated. */
  readonly isolated?: boolean;
  /** Generated files, served under /__extra/<name> (inputs of measurements). */
  readonly extra?: ReadonlyMap<string, Uint8Array>;
}

/** Starts a server for the harness files under /__harness/ and repository files under /__files/. */
export async function startHarnessServer(files: HarnessFiles, options: HarnessServerOptions = {}): Promise<HarnessServer> {
  const site = siteHeaders();
  const headers =
    options.isolated === false
      ? Object.fromEntries(Object.entries(site).filter(([name]) => !/^cross-origin-(opener|embedder)-policy$/i.test(name)))
      : site;
  const server = createServer((request: IncomingMessage, response: ServerResponse) => {
    const path = decodeURIComponent(new URL(request.url ?? '/', 'http://localhost').pathname);
    const send = (status: number, type: string, body: string | Uint8Array) => {
      response.writeHead(status, { ...headers, 'Content-Type': type });
      response.end(body);
    };
    if (path.startsWith('/__harness/')) {
      const file = files.get(path.slice('/__harness/'.length) || 'index.html');
      if (file === undefined) send(404, 'text/plain', `not in the harness: ${path}`);
      else send(200, file.type, file.body);
      return;
    }
    if (path.startsWith('/__extra/')) {
      const file = options.extra?.get(path.slice('/__extra/'.length));
      if (file === undefined) send(404, 'text/plain', `no generated file ${path}`);
      else send(200, 'application/octet-stream', file);
      return;
    }
    if (path.startsWith('/__files/')) {
      const relative = path.slice('/__files/'.length);
      if (relative.split('/').some((name) => name === '..' || name === '')) {
        send(400, 'text/plain', 'bad path');
        return;
      }
      try {
        send(200, 'application/octet-stream', readFileSync(REPOSITORY + relative));
      } catch {
        send(404, 'text/plain', `not found: ${relative}`);
      }
      return;
    }
    send(404, 'text/plain', 'not found');
  });
  await new Promise<void>((resolve) => server.listen(0, '127.0.0.1', resolve));
  const { port } = server.address() as AddressInfo;
  return {
    origin: `http://127.0.0.1:${port}`,
    close: () =>
      new Promise<void>((resolve, reject) => {
        server.closeAllConnections();
        server.close((error) => (error ? reject(error) : resolve()));
      }),
  };
}

/** Opens the harness page of `server` and waits until the harness is ready. */
export async function openHarness(context: BrowserContext, server: HarnessServer): Promise<Page> {
  const page = await context.newPage();
  // The harness logs each case and search as it finishes.
  page.on('console', (message) => console.log(message.text()));
  await page.goto(`${server.origin}/__harness/index.html`);
  await page.waitForFunction(() => window.losatHarness !== undefined);
  return page;
}
