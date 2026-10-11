// V-PRIV's capture (plan §5.10; design §13.3): every request of a flow, seen from both ends.
//
// - The server: the site is served by a server that logs every request it receives (method, URL,
//   headers, body), so that the requests of pages, dedicated and nested workers and a service
//   worker are all seen in every browser, whatever the browser reports to Playwright. It sends
//   the site's response headers (public/_headers) on every response, and adds one directive to
//   the CSP, `report-uri`, so that a violation in a worker (where no page listens) reaches the log.
// - The browser: every request that the context reports to Playwright, and what each page notes
//   about itself (an init script): its title and address whenever they change, and every
//   `securitypolicyviolation` event.
//
// The server is `startSiteServer` alone, so that it can be replaced by the release layout's local
// server (WP-R: build/site-server.ts and its request-log hook) in a few lines.
import { createServer, type IncomingMessage, type ServerResponse } from 'node:http';
import type { AddressInfo } from 'node:net';
import { existsSync, readFileSync, writeFileSync } from 'node:fs';
import { extname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { test as base, type BrowserContext, type Page, type Request, type TestInfo } from '@playwright/test';
import { siteHeaders } from '../../../build/headers';
import { canaryForms, filesUnder, findings, type Canary, type Exchange, type PageNote, type RequestRules } from './privacy-scan';

/** The production build that `npm run build` writes (playwright.config.ts builds it before the tests). */
export const DIST = fileURLToPath(new URL('../../../dist/', import.meta.url));
/** Where the server asks browsers to report CSP violations; any request to it is a finding. */
export const CSP_REPORT_PATH = '/__privacy/csp-report';

const TYPES: Readonly<Record<string, string>> = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.mjs': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json',
  '.wasm': 'application/wasm',
  '.svg': 'image/svg+xml',
  '.png': 'image/png',
  '.ico': 'image/x-icon',
  '.txt': 'text/plain; charset=utf-8',
  '.webmanifest': 'application/manifest+json',
};
/**
 * Paths that are not files of the site and that a browser asks for by itself: Firefox asks for
 * /favicon.ico (index.html names no icon). The canary check applies to them as to every request.
 */
export const BROWSER_PATHS: ReadonlySet<string> = new Set(['/favicon.ico']);
/** Files of the build directory that Cloudflare reads and never serves. */
const NOT_SERVED = new Set(['_headers', '_redirects']);

export interface SiteServer {
  /** For example http://127.0.0.1:41234 */
  readonly origin: string;
  /** The paths that the site serves: `/` and `/<file>` for each file of the build. */
  readonly files: ReadonlySet<string>;
  /** Every request received, in order. */
  readonly requests: readonly Exchange[];
  close(): Promise<void>;
}

export interface SiteServerOptions {
  /** Files that are not the site's, served under their paths (the capture's own test plants them). */
  readonly extra?: ReadonlyMap<string, { readonly type: string; readonly body: string }>;
}

/** Serves the build directory as the site, logging every request. */
export async function startSiteServer(root = DIST, options: SiteServerOptions = {}): Promise<SiteServer> {
  if (!existsSync(join(root, 'index.html'))) throw new Error(`no build in ${root}: run npm run build`);
  const served = filesUnder(root).filter((file) => !NOT_SERVED.has(file));
  const files = new Set(['/', ...served.map((file) => `/${file}`)]);
  const site = siteHeaders();
  const headers = {
    ...site,
    'Content-Security-Policy': `${site['Content-Security-Policy']}; report-uri ${CSP_REPORT_PATH}`,
  };
  const requests: Exchange[] = [];
  let origin = '';
  const server = createServer((request: IncomingMessage, response: ServerResponse) => {
    const chunks: Buffer[] = [];
    request.on('data', (chunk: Buffer) => chunks.push(chunk));
    request.on('end', () => {
      requests.push({
        seenBy: 'server',
        method: request.method ?? '',
        url: `${origin}${request.url ?? ''}`,
        headers: Object.fromEntries(Object.entries(request.headers).map(([name, value]) => [name, [value ?? ''].flat().join(', ')])),
        body: Buffer.concat(chunks),
        kind: 'server',
      });
      const path = new URL(request.url ?? '/', origin).pathname;
      const planted = options.extra?.get(path);
      if (planted !== undefined) {
        response.writeHead(200, { ...headers, 'Content-Type': planted.type });
        response.end(planted.body);
        return;
      }
      const file = path === '/' ? 'index.html' : decodeURIComponent(path.slice(1));
      if (!files.has(path)) {
        response.writeHead(404, { ...headers, 'Content-Type': 'text/plain' });
        response.end('not found');
        return;
      }
      response.writeHead(200, { ...headers, 'Content-Type': TYPES[extname(file)] ?? 'application/octet-stream' });
      response.end(request.method === 'HEAD' ? undefined : readFileSync(join(root, file)));
    });
  });
  await new Promise<void>((resolve) => server.listen(0, '127.0.0.1', resolve));
  origin = `http://127.0.0.1:${(server.address() as AddressInfo).port}`;
  return {
    origin,
    files,
    requests,
    close: () =>
      new Promise<void>((resolve, reject) => {
        server.closeAllConnections();
        server.close((error) => (error ? reject(error) : resolve()));
      }),
  };
}

const BINDING = '__losatPrivacyNote';

/**
 * Runs first in every document of the context: reports the title and the address whenever they
 * change (the head's mutations, and a look every 200 ms for the address), and every CSP violation.
 */
const NOTE_SCRIPT = `(() => {
  const send = (kind, value) => { try { void window.${BINDING}({ kind, value: String(value) }); } catch {} };
  let title, href;
  const look = () => {
    if (document.title !== title) { title = document.title; send('title', title); }
    if (location.href !== href) { href = location.href; send('href', href); }
  };
  document.addEventListener('securitypolicyviolation', (event) =>
    send('csp', event.effectiveDirective + ' blocked ' + event.blockedURI + ' (' + event.sourceFile + ':' + event.lineNumber + ')'));
  document.addEventListener('DOMContentLoaded', () => {
    look();
    if (document.head) new MutationObserver(look).observe(document.head, { subtree: true, childList: true, characterData: true });
  });
  setInterval(look, 200);
  look();
})();`;

/** The requests and notes of one test, from the server and the browser. */
export class PrivacyCapture {
  readonly browser: Exchange[] = [];
  readonly notes: PageNote[] = [];
  private readonly pending: Promise<void>[] = [];

  constructor(readonly server: SiteServer) {}

  get origin(): string {
    return this.server.origin;
  }

  async attach(context: BrowserContext): Promise<void> {
    await context.exposeBinding(BINDING, (_source, note: PageNote) => {
      this.notes.push(note);
    });
    await context.addInitScript(NOTE_SCRIPT);
    context.on('request', (request) => this.pending.push(this.record(request)));
    const watch = (page: Page) => page.on('websocket', (socket) => this.notes.push({ kind: 'websocket', value: socket.url() }));
    context.pages().forEach(watch);
    context.on('page', watch);
  }

  private async record(request: Request): Promise<void> {
    let headers: Record<string, string>;
    try {
      headers = await request.allHeaders();
    } catch {
      headers = request.headers();
    }
    this.browser.push({
      seenBy: 'browser',
      method: request.method(),
      url: request.url(),
      headers,
      body: request.postDataBuffer() ?? Buffer.alloc(0),
      kind: request.serviceWorker() !== null ? 'service worker' : request.resourceType(),
    });
  }

  /** Notes a page's title and address now (the init script notes changes too). */
  async checkpoint(page: Page): Promise<void> {
    const { title, href } = await page.evaluate(() => ({ title: document.title, href: location.href }));
    this.notes.push({ kind: 'title', value: title }, { kind: 'href', value: href });
  }

  /** Waits for the browser's request events in flight. */
  async settle(): Promise<void> {
    while (this.pending.length > 0) await Promise.all(this.pending.splice(0));
  }

  /** Every V-PRIV finding of the requests so far (support/privacy-scan.ts `findings`). */
  async findings(canary: Canary, rules: Partial<Pick<RequestRules, 'allowBody' | 'extraPaths'>> = {}): Promise<string[]> {
    await this.settle();
    return findings([...this.server.requests, ...this.browser], this.notes, canaryForms(canary), {
      origin: this.origin,
      files: this.server.files,
      cspReportPath: CSP_REPORT_PATH,
      ...rules,
      extraPaths: new Set([...BROWSER_PATHS, ...(rules.extraPaths ?? [])]),
    });
  }

  /** The paths the server received, in order. */
  serverPaths(): string[] {
    return this.server.requests.map((exchange) => new URL(exchange.url).pathname);
  }

  /**
   * Which requests the browser reported to Playwright and which reached only the server (or only
   * the browser: the browser's cache; blob: and data: URLs, read in the page, are counted apart),
   * for the gate record.
   */
  summary(): object {
    const key = (exchange: Exchange) => `${exchange.method} ${exchange.url}`;
    const inPage = (exchange: Exchange) => /^(blob|data):/.test(exchange.url);
    const unmatched = new Map<string, number>();
    for (const exchange of this.browser.filter((each) => !inPage(each))) unmatched.set(key(exchange), (unmatched.get(key(exchange)) ?? 0) + 1);
    const serverOnly: string[] = [];
    for (const exchange of this.server.requests) {
      const count = unmatched.get(key(exchange)) ?? 0;
      if (count > 0) unmatched.set(key(exchange), count - 1);
      else serverOnly.push(`${exchange.method} ${new URL(exchange.url).pathname}`);
    }
    const short = (url: string) => (url.startsWith(this.origin) ? url.slice(this.origin.length) : url.slice(0, 80));
    return {
      server: this.server.requests.map((exchange) => `${exchange.method} ${short(exchange.url)}`),
      browser: this.browser.filter((each) => !inPage(each)).map((exchange) => `${exchange.method} ${short(exchange.url)} [${exchange.kind}]`),
      browserBlobOrDataUrls: this.browser.filter(inPage).length,
      serverOnly,
      browserOnly: [...unmatched.entries()].filter(([, count]) => count > 0).map(([url, count]) => `${short(url.slice(url.indexOf(' ') + 1))} x${count}`),
      cspReports: this.server.requests
        .filter((exchange) => new URL(exchange.url).pathname === CSP_REPORT_PATH)
        .map((exchange) => Buffer.from(exchange.body).toString('utf8')),
      notes: this.notes.filter((note) => note.kind !== 'href' || !note.value.startsWith(this.origin)).map((note) => `${note.kind}: ${note.value}`),
    };
  }

  /** Writes the summary to the test's output folder (requests.json) and attaches it to the report. */
  async report(testInfo: TestInfo): Promise<void> {
    await this.settle();
    const path = testInfo.outputPath('requests.json');
    writeFileSync(path, `${JSON.stringify(this.summary(), null, 2)}\n`);
    await testInfo.attach('requests', { path, contentType: 'application/json' });
  }
}

/**
 * Playwright's `test` with `capture`: the site served by a logging server, the context watched
 * before the page is made.
 */
export const test = base.extend<{ capture: PrivacyCapture }>({
  page: async ({ capture, page }, use) => {
    void capture;
    await use(page);
  },
  capture: async ({ context }, use, testInfo) => {
    const server = await startSiteServer();
    const capture = new PrivacyCapture(server);
    await capture.attach(context);
    try {
      await use(capture);
    } finally {
      await capture.report(testInfo).catch(() => undefined);
      await server.close();
    }
  },
});
