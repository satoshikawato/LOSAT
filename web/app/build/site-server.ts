// A local server for a LOSAT Web site directory (`dist/`, or a site that build/assemble-site.ts
// made) that answers as Cloudflare Pages does (docs/web/release_manifest.md, "The local
// server"), for the E2E tests and `npm run preview`:
//
// - the rules of the site's own `_headers` apply by path to every response: 200, 304, 404,
//   redirects; Pages' default `Cache-Control: public, max-age=0, must-revalidate` where no rule
//   sets one;
// - a strong ETag on every 200, and 304 Not Modified for a matching If-None-Match;
// - `/` and `/dir/` serve `index.html`; `/name` serves `name.html`; `/name.html` and
//   `/dir/index.html` redirect (308) to `/name` and `/dir/`;
// - a missing file is 404 with the nearest `404.html` (no single-page fallback);
// - `_headers`, `_redirects` and `_routes.json` are not served.
//
// It reads the directory on every request, so a test may change the site while it runs (swap
// `releases.json`, add a release), and it can report every request it receives (`onRequest`).
//
//   node build/site-server.ts <site directory> [--port 4173] [--host localhost]
import { createHash } from 'node:crypto';
import { readFile, stat } from 'node:fs/promises';
import { createServer, type IncomingHttpHeaders, type IncomingMessage, type ServerResponse } from 'node:http';
import type { AddressInfo } from 'node:net';
import { extname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { headersFor, parseHeaderRules, type HeaderRule } from './headers.ts';

/** A request as the server received it, and its status. */
export interface SiteRequest {
  readonly method: string;
  /** The request target as sent: the path and the query. */
  readonly url: string;
  /** The path of the request target. */
  readonly path: string;
  readonly headers: IncomingHttpHeaders;
  readonly body: Uint8Array;
  readonly status: number;
}

export interface SiteServerOptions {
  /** The site directory. */
  readonly root: string;
  /** 0 (the default) takes a free port. */
  readonly port?: number;
  /** 127.0.0.1 by default. */
  readonly host?: string;
  /** Called with every request after the server has answered it. */
  readonly onRequest?: (request: SiteRequest) => void;
}

export interface SiteServer {
  /** For example http://127.0.0.1:41234 */
  readonly origin: string;
  readonly port: number;
  close(): Promise<void>;
}

/** Pages' Cache-Control when no rule of `_headers` gives one. */
export const DEFAULT_CACHE_CONTROL = 'public, max-age=0, must-revalidate';
const MAX_BODY = 16 * 1024 * 1024;
const NOT_SERVED = new Set(['/_headers', '/_redirects', '/_routes.json']);

const TYPES: Readonly<Record<string, string>> = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.mjs': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json',
  '.map': 'application/json',
  '.wasm': 'application/wasm',
  '.svg': 'image/svg+xml',
  '.png': 'image/png',
  '.ico': 'image/x-icon',
  '.txt': 'text/plain; charset=utf-8',
  '.webmanifest': 'application/manifest+json',
};

/** The content type of a file by its extension. */
export function contentType(file: string): string {
  return TYPES[extname(file).toLowerCase()] ?? 'application/octet-stream';
}

interface Answer {
  readonly status: number;
  readonly headers: Record<string, string>;
  readonly body?: Uint8Array;
}

async function readIfFile(path: string): Promise<{ bytes: Uint8Array; key: string } | undefined> {
  try {
    const info = await stat(path);
    if (!info.isFile()) return undefined;
    return { bytes: await readFile(path), key: `${path}\0${info.mtimeMs}\0${info.size}` };
  } catch {
    return undefined;
  }
}

/** Whether an If-None-Match value names `etag` (weak comparison, as for GET). */
function matchesEtag(ifNoneMatch: string | undefined, etag: string): boolean {
  if (ifNoneMatch === undefined) return false;
  const bare = etag.replace(/^W\//, '');
  return ifNoneMatch.split(',').some((tag) => {
    const value = tag.trim();
    return value === '*' || value.replace(/^W\//, '') === bare;
  });
}

export async function startSiteServer(options: SiteServerOptions): Promise<SiteServer> {
  const root = resolve(options.root);
  const etags = new Map<string, string>();
  let rules: { key: string; rules: HeaderRule[] } | undefined;

  const currentRules = async (): Promise<HeaderRule[]> => {
    const file = await readIfFile(join(root, '_headers'));
    if (file === undefined) return [];
    if (rules?.key !== file.key) rules = { key: file.key, rules: parseHeaderRules(new TextDecoder().decode(file.bytes)) };
    return rules.rules;
  };

  const etagOf = (file: { bytes: Uint8Array; key: string }): string => {
    let etag = etags.get(file.key);
    if (etag === undefined) {
      etag = `"${createHash('sha256').update(file.bytes).digest('hex').slice(0, 32)}"`;
      etags.set(file.key, etag);
    }
    return etag;
  };

  /** The site file of a path, a redirect, or nothing. */
  const locate = async (path: string): Promise<{ file: string; data: { bytes: Uint8Array; key: string } } | { redirect: string } | undefined> => {
    const segments = path.split('/');
    if (NOT_SERVED.has(path) || segments.some((segment) => segment === '..' || segment === '.') || /[\\\0]/.test(path)) return undefined;
    if (path.endsWith('/index.html')) return { redirect: path.slice(0, -'index.html'.length) };
    if (path.endsWith('.html')) return { redirect: path.slice(0, -'.html'.length) };
    const candidates = path.endsWith('/') ? [`${path}index.html`] : [path, `${path}.html`];
    for (const candidate of candidates) {
      const data = await readIfFile(join(root, candidate));
      if (data !== undefined) return { file: candidate, data };
    }
    if (!path.endsWith('/') && (await readIfFile(join(root, path, 'index.html'))) !== undefined) return { redirect: `${path}/` };
    return undefined;
  };

  /** The nearest 404.html above a path, as Pages looks for it. */
  const notFoundPage = async (path: string): Promise<Uint8Array | undefined> => {
    const segments = path.split('/').slice(1, -1);
    for (let depth = segments.length; depth >= 0; depth--) {
      const data = await readIfFile(join(root, ...segments.slice(0, depth), '404.html'));
      if (data !== undefined) return data.bytes;
    }
    return undefined;
  };

  const answer = async (method: string, path: string, search: string, requestHeaders: IncomingHttpHeaders): Promise<Answer> => {
    if (method !== 'GET' && method !== 'HEAD') return { status: 405, headers: { Allow: 'GET, HEAD' } };
    const found = await locate(path);
    if (found === undefined) {
      const page = await notFoundPage(path);
      return page === undefined
        ? { status: 404, headers: { 'Content-Type': 'text/plain; charset=utf-8' }, body: new TextEncoder().encode('Not found\n') }
        : { status: 404, headers: { 'Content-Type': 'text/html; charset=utf-8' }, body: page };
    }
    if ('redirect' in found) return { status: 308, headers: { Location: `${found.redirect}${search}` } };
    const etag = etagOf(found.data);
    if (matchesEtag(requestHeaders['if-none-match'], etag)) return { status: 304, headers: { ETag: etag } };
    return { status: 200, headers: { 'Content-Type': contentType(found.file), ETag: etag }, body: found.data.bytes };
  };

  const handle = async (request: IncomingMessage, response: ServerResponse): Promise<void> => {
    const chunks: Buffer[] = [];
    let size = 0;
    for await (const chunk of request as AsyncIterable<Buffer>) {
      size += chunk.length;
      if (size <= MAX_BODY) chunks.push(chunk);
    }
    const target = new URL(request.url ?? '/', 'http://site.invalid');
    let path = target.pathname;
    try {
      path = decodeURIComponent(target.pathname);
    } catch {
      // A path that is not valid percent-encoding is matched as sent, and is not a file.
    }
    const method = request.method ?? 'GET';
    let ruleHeaders: Record<string, string> = {};
    let result: Answer;
    try {
      // The rules match the path as requested (Pages matches the URL, not the file it serves).
      ruleHeaders = headersFor(target.pathname, await currentRules());
      result = size > MAX_BODY ? { status: 413, headers: {} } : await answer(method, path, target.search, request.headers);
    } catch (error) {
      result = { status: 500, headers: { 'Content-Type': 'text/plain; charset=utf-8' }, body: new TextEncoder().encode(`${String(error)}\n`) };
    }
    const headers: Record<string, string> = { ...ruleHeaders, ...result.headers };
    if (!Object.keys(headers).some((name) => name.toLowerCase() === 'cache-control')) headers['Cache-Control'] = DEFAULT_CACHE_CONTROL;
    if (result.body !== undefined) headers['Content-Length'] = String(result.body.length);
    response.writeHead(result.status, headers);
    response.end(method === 'HEAD' || result.body === undefined ? undefined : result.body);
    options.onRequest?.({
      method,
      url: request.url ?? '/',
      path: target.pathname,
      headers: request.headers,
      body: new Uint8Array(Buffer.concat(chunks)),
      status: result.status,
    });
  };

  const server = createServer((request, response) => {
    handle(request, response).catch(() => {
      if (!response.headersSent) response.writeHead(500);
      response.end();
    });
  });
  const host = options.host ?? '127.0.0.1';
  await new Promise<void>((resolveListen, reject) => {
    server.once('error', reject);
    server.listen(options.port ?? 0, host, () => resolveListen());
  });
  const { port } = server.address() as AddressInfo;
  return {
    origin: `http://${host.includes(':') ? `[${host}]` : host}:${port}`,
    port,
    close: () =>
      new Promise<void>((resolveClose, reject) => {
        server.closeAllConnections();
        server.close((error) => (error ? reject(error) : resolveClose()));
      }),
  };
}

/** `node build/site-server.ts <site directory> [--port N] [--host H]` */
async function main(args: readonly string[]): Promise<void> {
  let root: string | undefined;
  let port = 4173;
  let host = 'localhost';
  for (let i = 0; i < args.length; i++) {
    const arg = args[i]!;
    if (arg === '--port') port = Number(args[++i]);
    else if (arg === '--host') host = args[++i] ?? host;
    else if (root === undefined && !arg.startsWith('--')) root = arg;
    else throw new Error(`site-server: unknown argument ${arg}`);
  }
  if (root === undefined || !Number.isInteger(port)) throw new Error('usage: node build/site-server.ts <site directory> [--port N] [--host H]');
  const server = await startSiteServer({ root, port, host });
  console.log(`LOSAT Web site ${resolve(root)} at ${server.origin}/`);
  const stop = () => void server.close().finally(() => process.exit(0));
  process.on('SIGINT', stop);
  process.on('SIGTERM', stop);
}

if (process.argv[1] !== undefined && resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  main(process.argv.slice(2)).catch((error: unknown) => {
    console.error(error instanceof Error ? error.message : error);
    process.exit(1);
  });
}
