// V-PRIV's capture (plan §5.10; design §13.3): every request of a flow, seen from both ends.
//
// - The server: the site is served by the site's local server (build/site-server.ts, which answers
//   as Cloudflare Pages does), whose request hook logs every request it receives (method, URL,
//   headers, body), so that the requests of pages, dedicated and nested workers and the service
//   worker are all seen in every browser, whatever the browser reports to Playwright. The served
//   directory is a view of the site (links to its files) whose `_headers` adds one directive to the
//   CSP, `report-uri`, so that a violation in a worker (where no page listens) reaches the log too;
//   pointing the view at another site directory is a new deployment at the same origin. A small
//   front answers /api/ (the measurement endpoint, a Pages Function on the deployed site) and logs
//   it the same way.
// - The browser: every request that the context reports to Playwright, and what each page notes
//   about itself (an init script): its title and address whenever they change, and every
//   `securitypolicyviolation` event.
import { lstatSync, mkdirSync, mkdtempSync, readdirSync, readFileSync, rmSync, symlinkSync, unlinkSync, writeFileSync } from 'node:fs';
import { createServer, request as forward } from 'node:http';
import type { AddressInfo } from 'node:net';
import { tmpdir } from 'node:os';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { test as base, type BrowserContext, type Page, type Request, type TestInfo } from '@playwright/test';
import { siteHeaders } from '../../../build/headers';
import { SERVICE_WORKER_FILE } from '../../../build/service-worker';
import { readSiteIndex, RELEASES_DIR, SHELL_FILES, SITE_INDEX_FILE } from '../../../build/site-files';
import { startSiteServer, type SiteRequest } from '../../../build/site-server';
import { parseReleaseManifest, RELEASE_MANIFEST_NAME } from '../../../src/domain/release-manifest';
import { canaryForms, findings, type Canary, type Exchange, type PageNote, type RequestRules } from './privacy-scan';

/** The site that `npm run build` writes (playwright.config.ts builds it before the tests). */
export const DIST = fileURLToPath(new URL('../../../dist/', import.meta.url));
/** Where the server asks browsers to report CSP violations; any request to it is a finding. */
export const CSP_REPORT_PATH = '/__privacy/csp-report';

/**
 * Paths that are not files of the site and that a browser asks for by itself: Firefox asks for
 * /favicon.ico (index.html names no icon). The canary check applies to them as to every request.
 */
export const BROWSER_PATHS: ReadonlySet<string> = new Set(['/favicon.ico']);
/** The site-level files outside the shell (docs/web/release_manifest.md): the site index, the service worker, the 404 page. */
const OTHER_SITE_FILES: readonly (readonly [path: string, file: string])[] = [
  [`/${SITE_INDEX_FILE}`, SITE_INDEX_FILE],
  [`/${SERVICE_WORKER_FILE}`, SERVICE_WORKER_FILE],
  // Pages serves 404.html as the body of a 404, and at /404 (/404.html redirects there).
  ['/404', '404.html'],
];

export interface SiteFile {
  /** The file in the site directory. */
  readonly file: string;
  /** The URL path that serves it. */
  readonly path: string;
}

/**
 * Every file that a site serves, as the site lists it: the shell of releases.json, the other
 * site-level files, and each release's manifest with the files it lists.
 */
export function siteFiles(site: string): SiteFile[] {
  const index = readSiteIndex(site);
  const shell = new Map(SHELL_FILES);
  const files: SiteFile[] = (index.shell ?? SHELL_FILES.map(([path]) => ({ path }))).map(({ path }) => ({ file: shell.get(path) ?? path.slice(1), path }));
  files.push(...OTHER_SITE_FILES.map(([path, file]) => ({ file, path })));
  for (const { id } of index.releases) {
    const directory = `${RELEASES_DIR}/${id}`;
    const manifest = parseReleaseManifest(JSON.parse(readFileSync(join(site, directory, RELEASE_MANIFEST_NAME), 'utf8')), id);
    for (const file of [RELEASE_MANIFEST_NAME, ...manifest.assets.map((asset) => asset.path)]) files.push({ file: `${directory}/${file}`, path: `/${directory}/${file}` });
  }
  return files;
}

/** The site's `_headers` with `report-uri` added to each Content-Security-Policy. */
function reportingHeaders(text: string): string {
  return text.replace(/^([ \t]+Content-Security-Policy:.*\S)[ \t]*$/gm, `$1; report-uri ${CSP_REPORT_PATH}`);
}

export interface CaptureSite {
  /** For example http://127.0.0.1:41234 */
  readonly origin: string;
  /** Every request received, in order. */
  readonly requests: readonly Exchange[];
  /** The URL paths of the files of every site deployed so far. */
  readonly paths: ReadonlySet<string>;
  /** From now on, the same origin serves another site directory (a new deployment). */
  deploy(site: string): void;
  close(): Promise<void>;
}

export interface CaptureSiteOptions {
  /** Files that are not the site's, served under their paths (the capture's own test plants them). */
  readonly extra?: ReadonlyMap<string, string>;
  /**
   * The measurement endpoint that the served releases.json names (`analytics.endpoint`), as a site
   * built or assembled with LOSAT_WEB_ANALYTICS_ENDPOINT names it.
   */
  readonly analyticsEndpoint?: string;
}

/** The headers of the answer to a request under /api/: the site's, as Pages adds them to a Function's answer. */
const API_HEADERS = Object.fromEntries(Object.entries(siteHeaders()).concat([['Cache-Control', 'no-store']]));

/**
 * Serves a site directory with build/site-server.ts, logging every request through its
 * `onRequest` hook. In front of it, requests under /api/ are logged and answered 204 here, as the
 * site's Pages Function answers /api/visit (functions/api/visit.ts) and the site server answers
 * no POST.
 */
export async function startCaptureSite(site = DIST, options: CaptureSiteOptions = {}): Promise<CaptureSite> {
  const view = mkdtempSync(join(tmpdir(), 'losat-vpriv-site-'));
  const paths = new Set<string>();
  const deploy = (directory: string) => {
    for (const entry of readdirSync(view)) {
      const path = join(view, entry);
      if (lstatSync(path).isDirectory()) rmSync(path, { recursive: true });
      else unlinkSync(path);
    }
    for (const entry of readdirSync(directory)) {
      if (entry === '_headers') continue;
      if (entry === SITE_INDEX_FILE && options.analyticsEndpoint !== undefined) {
        const index = JSON.parse(readFileSync(join(directory, entry), 'utf8')) as Record<string, unknown>;
        writeFileSync(join(view, entry), `${JSON.stringify({ ...index, analytics: { endpoint: options.analyticsEndpoint } }, null, 2)}\n`);
      } else symlinkSync(resolve(directory, entry), join(view, entry));
    }
    writeFileSync(join(view, '_headers'), reportingHeaders(readFileSync(join(directory, '_headers'), 'utf8')));
    for (const [path, body] of options.extra ?? []) {
      mkdirSync(dirname(join(view, path)), { recursive: true });
      writeFileSync(join(view, path), body);
    }
    for (const file of siteFiles(directory)) paths.add(file.path);
  };
  deploy(site);
  const requests: Exchange[] = [];
  let origin = '';
  const log = (request: SiteRequest) =>
    requests.push({
      seenBy: 'server',
      method: request.method,
      url: `${origin}${request.url}`,
      headers: Object.fromEntries(Object.entries(request.headers).map(([name, value]) => [name, [value ?? ''].flat().join(', ')])),
      body: request.body,
      kind: `server ${request.status}`,
    });
  const server = await startSiteServer({ root: view, onRequest: log });
  const front = createServer((request, response) => {
    const chunks: Buffer[] = [];
    request.on('data', (chunk: Buffer) => chunks.push(chunk));
    if (new URL(request.url ?? '/', server.origin).pathname.startsWith('/api/')) {
      request.on('end', () => {
        response.writeHead(204, API_HEADERS);
        response.end();
        log({ method: request.method ?? '', url: request.url ?? '/', path: '', headers: request.headers, body: new Uint8Array(Buffer.concat(chunks)), status: 204 });
      });
      return;
    }
    const upstream = forward(`${server.origin}${request.url ?? '/'}`, { method: request.method, headers: request.headers }, (answer) => {
      response.writeHead(answer.statusCode ?? 502, answer.headers);
      answer.pipe(response);
    });
    upstream.on('error', () => {
      if (!response.headersSent) response.writeHead(502);
      response.end();
    });
    request.pipe(upstream);
  });
  await new Promise<void>((resolveListen) => front.listen(0, '127.0.0.1', resolveListen));
  origin = `http://127.0.0.1:${(front.address() as AddressInfo).port}`;
  return {
    origin,
    requests,
    paths,
    deploy,
    close: async () => {
      front.closeAllConnections();
      await new Promise<void>((resolveClose) => front.close(() => resolveClose()));
      await server.close();
      rmSync(view, { recursive: true, force: true });
    },
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

  constructor(readonly server: CaptureSite) {}

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
    // allHeaders() adds what the network stack sends (cookies among them), but never resolves for
    // some requests that a service worker answers (Chromium); the server logs every header anyway.
    let headers: Record<string, string>;
    try {
      headers = await Promise.race([request.allHeaders(), new Promise<Record<string, string>>((resolve) => setTimeout(() => resolve(request.headers()), 2000))]);
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
      files: this.server.paths,
      cspReportPath: CSP_REPORT_PATH,
      ...rules,
      extraPaths: new Set([...BROWSER_PATHS, ...(rules.extraPaths ?? [])]),
    });
  }

  /** Where the requests stand now (once the browser's events in flight are recorded). */
  async mark(): Promise<{ readonly server: number; readonly browser: number }> {
    await this.settle();
    return { server: this.server.requests.length, browser: this.browser.length };
  }

  /** The requests since a mark, as the server and the browser saw them. */
  async since(mark: { readonly server: number; readonly browser: number }): Promise<Exchange[]> {
    await this.settle();
    return [...this.server.requests.slice(mark.server), ...this.browser.slice(mark.browser)];
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
 * Playwright's `test` with `capture`: the site (dist/, served as `captureSite` says) behind a
 * logging server, the context watched before the page is made.
 */
export const test = base.extend<{ capture: PrivacyCapture; captureSite: CaptureSiteOptions }>({
  captureSite: [{}, { option: true }],
  page: async ({ capture, page }, use) => {
    void capture;
    await use(page);
  },
  capture: async ({ context, captureSite }, use, testInfo) => {
    const server = await startCaptureSite(DIST, captureSite);
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
