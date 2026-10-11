// The origin of V-OFF (tests/e2e/offline.spec.ts): a server on a fixed port in front of the
// site servers (build/site-server.ts), so that one origin can
// - serve another deployment of the site (`upstream`: the site server of another directory),
// - go offline and come back (`setOnline`): the port is closed and every open connection cut,
//   a real network failure for the page and its service worker alike (Playwright's offline
//   emulation does not reach a service worker's requests in Firefox, and in WebKit it fails the
//   navigation itself: logs/wp-u/progress.md, probe 1a),
// - cut the connection in the middle of a file's body (`cut`),
// - answer a path with an HTML page instead of the file (`pages`: a proxy's or a captive portal's
//   page, code review 3 L2),
// - hold a file's request until the test lets it go (`hold`: the screen records of a release
//   being prepared, tests/e2e/screens.spec.ts),
// - keep the release files out of the browser's HTTP cache (`releasesNoStore`, V-OFF): a 200
//   under /r/ is then sent with `Cache-Control: no-store` instead of the site's immutable rule
//   (public/_headers keeps that rule; 404s keep it too, as V-OFF 6c needs), so that an offline
//   start can only be served by the service worker's caches, for the page's requests and its
//   workers' alike,
// - and record every request it receives, with the bytes of the body it sent (the delivery
//   measurement, tests/e2e/delivery-measure.spec.ts).
import { createServer, request as forward, type IncomingMessage, type Server, type ServerResponse } from 'node:http';

export interface FrontRequest {
  readonly method: string;
  readonly path: string;
  readonly status: number;
  /** The bytes of the body sent so far (the whole body once `ended` is set). */
  bytes: number;
  /** When the answer was sent whole or cut (Date.now()), or undefined while it is sent. */
  ended?: number;
}

export interface FrontServer {
  readonly origin: string;
  /** The origin of the site server that answers. */
  upstream: string;
  /** Paths whose 200 body is cut after half of it. */
  readonly cut: Set<string>;
  /** Paths answered with `200 text/html`, a page that is not the file, without asking the site server. */
  readonly pages: Set<string>;
  /** Whether a 200 under /r/ is sent with `Cache-Control: no-store` (default false: the site's rules). */
  releasesNoStore: boolean;
  readonly requests: FrontRequest[];
  /** Makes the requests of `path` wait until the returned function is called. */
  hold(path: string): () => void;
  setOnline(online: boolean): Promise<void>;
  close(): Promise<void>;
}

async function relay(front: FrontServer, held: ReadonlyMap<string, Promise<void>>, request: IncomingMessage, response: ServerResponse): Promise<void> {
  const path = new URL(request.url ?? '/', 'http://front.invalid').pathname;
  await held.get(path);
  if (front.pages.has(path)) {
    const page = Buffer.from('<!doctype html><title>Access blocked</title><p>This network does not allow this address.</p>\n');
    front.requests.push({ method: request.method ?? 'GET', path, status: 200, bytes: page.length, ended: Date.now() });
    response.writeHead(200, { 'content-type': 'text/html; charset=utf-8', 'cache-control': 'no-store', 'content-length': String(page.length) });
    response.end(page);
    return;
  }
  const upstream = forward(
    `${front.upstream}${request.url ?? '/'}`,
    { method: request.method, headers: { ...request.headers, host: new URL(front.upstream).host } },
    (answer) => {
      const status = answer.statusCode ?? 502;
      const record: FrontRequest = { method: request.method ?? 'GET', path, status, bytes: 0 };
      front.requests.push(record);
      const end = () => (record.ended ??= Date.now());
      response.on('finish', end);
      response.on('close', end);
      const headers =
        front.releasesNoStore && status === 200 && path.startsWith('/r/') ? { ...answer.headers, 'cache-control': 'no-store' } : answer.headers;
      if (status === 200 && front.cut.has(path)) {
        const chunks: Buffer[] = [];
        answer.on('data', (chunk: Buffer) => chunks.push(chunk));
        answer.on('end', () => {
          const body = Buffer.concat(chunks);
          response.writeHead(status, headers);
          const half = body.subarray(0, Math.floor(body.length / 2));
          record.bytes = half.length;
          response.write(half, () => response.socket?.destroy());
        });
        return;
      }
      answer.on('data', (chunk: Buffer) => (record.bytes += chunk.length));
      response.writeHead(status, headers);
      answer.pipe(response);
    },
  );
  upstream.on('error', () => {
    front.requests.push({ method: request.method ?? 'GET', path, status: 502, bytes: 0, ended: Date.now() });
    if (!response.headersSent) response.writeHead(502);
    response.end();
  });
  request.pipe(upstream);
}

export async function startFrontServer(port: number, upstream: string): Promise<FrontServer> {
  let server: Server | undefined;
  const held = new Map<string, Promise<void>>();
  const front: FrontServer = {
    origin: `http://localhost:${port}`,
    upstream,
    cut: new Set(),
    pages: new Set(),
    releasesNoStore: false,
    requests: [],
    hold(path: string): () => void {
      let release = () => {};
      held.set(path, new Promise<void>((resolve) => (release = resolve)));
      return () => {
        held.delete(path);
        release();
      };
    },
    async setOnline(online: boolean): Promise<void> {
      if (online === (server !== undefined)) return;
      if (!online) {
        const closing = server!;
        server = undefined;
        closing.closeAllConnections();
        await new Promise<void>((resolve) => closing.close(() => resolve()));
        return;
      }
      const opened = createServer((request, response) => void relay(front, held, request, response));
      await new Promise<void>((resolve, reject) => {
        opened.once('error', reject);
        opened.listen(port, 'localhost', () => resolve());
      });
      server = opened;
    },
    close: () => front.setOnline(false),
  };
  await front.setOnline(true);
  return front;
}
