// The origin of V-OFF (tests/e2e/offline.spec.ts): a server on a fixed port in front of the
// site servers (build/site-server.ts), so that one origin can
// - serve another deployment of the site (`upstream`: the site server of another directory),
// - go offline and come back (`setOnline`): the port is closed and every open connection cut,
//   a real network failure for the page and its service worker alike (Playwright's offline
//   emulation does not reach a service worker's requests in Firefox, and in WebKit it fails the
//   navigation itself: logs/wp-u/progress.md, probe 1a),
// - cut the connection in the middle of a file's body (`cut`),
// - and record every request it receives.
import { createServer, request as forward, type IncomingMessage, type Server, type ServerResponse } from 'node:http';

export interface FrontRequest {
  readonly method: string;
  readonly path: string;
  readonly status: number;
}

export interface FrontServer {
  readonly origin: string;
  /** The origin of the site server that answers. */
  upstream: string;
  /** Paths whose 200 body is cut after half of it. */
  readonly cut: Set<string>;
  readonly requests: FrontRequest[];
  setOnline(online: boolean): Promise<void>;
  close(): Promise<void>;
}

function relay(front: FrontServer, request: IncomingMessage, response: ServerResponse): void {
  const path = new URL(request.url ?? '/', 'http://front.invalid').pathname;
  const upstream = forward(
    `${front.upstream}${request.url ?? '/'}`,
    { method: request.method, headers: { ...request.headers, host: new URL(front.upstream).host } },
    (answer) => {
      const status = answer.statusCode ?? 502;
      front.requests.push({ method: request.method ?? 'GET', path, status });
      if (status === 200 && front.cut.has(path)) {
        const chunks: Buffer[] = [];
        answer.on('data', (chunk: Buffer) => chunks.push(chunk));
        answer.on('end', () => {
          const body = Buffer.concat(chunks);
          response.writeHead(status, answer.headers);
          response.write(body.subarray(0, Math.floor(body.length / 2)), () => response.socket?.destroy());
        });
        return;
      }
      response.writeHead(status, answer.headers);
      answer.pipe(response);
    },
  );
  upstream.on('error', () => {
    front.requests.push({ method: request.method ?? 'GET', path, status: 502 });
    if (!response.headersSent) response.writeHead(502);
    response.end();
  });
  request.pipe(upstream);
}

export async function startFrontServer(port: number, upstream: string): Promise<FrontServer> {
  let server: Server | undefined;
  const front: FrontServer = {
    origin: `http://localhost:${port}`,
    upstream,
    cut: new Set(),
    requests: [],
    async setOnline(online: boolean): Promise<void> {
      if (online === (server !== undefined)) return;
      if (!online) {
        const closing = server!;
        server = undefined;
        closing.closeAllConnections();
        await new Promise<void>((resolve) => closing.close(() => resolve()));
        return;
      }
      const opened = createServer((request, response) => relay(front, request, response));
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
