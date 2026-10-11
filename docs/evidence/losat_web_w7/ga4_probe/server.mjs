// The GA4 probe's servers (RESULTS.md): one process, three local origins (origins.ts gives
// them through the environment). Every request that reaches a server is recorded; the spec
// reads the record from G's `/__log?run=<id>`. No request leaves the machine.
import { readFileSync } from 'node:fs';
import { createServer } from 'node:http';

const env = (name) => {
  const value = process.env[name];
  if (!value) throw new Error(`${name} is not set (run through playwright.config.ts)`);
  return value;
};
const A = env('PROBE_A');
const M = env('PROBE_M');
const G = env('PROBE_G');
const C = env('PROBE_C');
const port = (origin) => Number(new URL(origin).port);

/** The `/*` rule of web/app/public/_headers (the same reading as web/app/build/headers.ts). */
function siteHeaders() {
  const text = readFileSync(new URL('../../../../web/app/public/_headers', import.meta.url), 'utf8');
  const rules = new Map();
  let current;
  for (const line of text.split(/\r?\n/)) {
    if (line.trim() === '' || line.trimStart().startsWith('#')) continue;
    if (!/^\s/.test(line)) {
      current = {};
      rules.set(line.trim(), current);
      continue;
    }
    const separator = line.indexOf(':');
    current[line.slice(0, separator).trim()] = line.slice(separator + 1).trim();
  }
  return rules.get('/*');
}
const SITE = siteHeaders();

/**
 * The response headers of https://www.googletagmanager.com/gtag/js?id=G-XXXXXXXXXX as
 * recorded on 2026-10-11 (google-headers.txt), less date, expiry, transport and server lines.
 */
const GTAG_HEADERS = {
  'content-type': 'application/javascript; charset=UTF-8',
  'access-control-allow-origin': '*',
  'access-control-allow-credentials': 'true',
  'access-control-allow-headers': 'Cache-Control',
  'vary': 'Accept-Encoding',
  'cache-control': 'private, max-age=900',
  'cross-origin-resource-policy': 'cross-origin',
  'x-xss-protection': '0',
};

/** M's response headers per variant (RESULTS.md). */
const M_VARIANT_HEADERS = {
  a1: {},
  a2: { 'Cross-Origin-Embedder-Policy': 'require-corp', 'Cross-Origin-Resource-Policy': 'cross-origin' },
  a3: { 'Cross-Origin-Embedder-Policy': 'credentialless', 'Cross-Origin-Resource-Policy': 'cross-origin' },
  a3n: { 'Cross-Origin-Embedder-Policy': 'credentialless' },
  a4: {},
};

const log = [];

function body(req) {
  return new Promise((resolve) => {
    const chunks = [];
    req.on('data', (chunk) => chunks.push(chunk));
    req.on('end', () => resolve(Buffer.concat(chunks).toString('utf8')));
  });
}

function record(server, req, url, text) {
  const pick = (name) => req.headers[name] ?? null;
  log.push({
    server,
    method: req.method,
    path: url.pathname,
    query: Object.fromEntries(url.searchParams),
    run: url.searchParams.get('run') ?? url.searchParams.get('probe_run'),
    headers: {
      origin: pick('origin'),
      referer: pick('referer'),
      cookie: pick('cookie'),
      'content-type': pick('content-type'),
      'sec-fetch-site': pick('sec-fetch-site'),
      'sec-fetch-mode': pick('sec-fetch-mode'),
      'sec-fetch-dest': pick('sec-fetch-dest'),
    },
    body: text,
  });
}

function send(res, status, headers, content = '') {
  res.writeHead(status, headers);
  res.end(content);
}

const STATIC = (name) => readFileSync(new URL(`site/${name}`, import.meta.url), 'utf8');

function serveA(req, res, url) {
  const csp = url.searchParams.get('csp') === 'frame'
    ? `${SITE['Content-Security-Policy']}; frame-src ${M}`
    : SITE['Content-Security-Policy'];
  const headers = { ...SITE, 'Content-Security-Policy': csp, 'Cache-Control': 'no-store' };
  if (url.pathname === '/a.html') return send(res, 200, { ...headers, 'Content-Type': 'text/html; charset=utf-8' }, STATIC('a.html'));
  if (url.pathname === '/a.js') return send(res, 200, { ...headers, 'Content-Type': 'text/javascript; charset=utf-8' }, STATIC('a.js'));
  if (url.pathname === '/__ready') return send(res, 200, { 'Content-Type': 'text/plain' }, 'ready');
  // The stand-in for the Pages Function of variant (b): records and answers 204.
  if (url.pathname === '/api/visit') return send(res, req.method === 'POST' ? 204 : 405, { ...headers });
  return send(res, 404, headers);
}

function serveM(req, res, url) {
  const variant = url.searchParams.get('variant') ?? 'a1';
  const own = {
    'Content-Security-Policy': `default-src 'none'; script-src 'self' ${G}; connect-src ${C}; img-src ${C}; base-uri 'none'; object-src 'none'; frame-ancestors ${A}`,
    'Cache-Control': 'no-store',
    'Referrer-Policy': 'strict-origin-when-cross-origin',
    ...(M_VARIANT_HEADERS[variant] ?? {}),
  };
  if (url.pathname === '/m.html') return send(res, 200, { ...own, 'Content-Type': 'text/html; charset=utf-8' }, STATIC('m.html'));
  if (url.pathname === '/m.js') {
    const config = `self.PROBE = ${JSON.stringify({ A, G, C })};\n`;
    // CORP cross-origin, not same-origin: a frame sandboxed without allow-same-origin has an
    // opaque origin, for which M's own script is cross-origin (and refused under COEP).
    return send(res, 200, { 'Content-Type': 'text/javascript; charset=utf-8', 'Cache-Control': 'no-store', 'Cross-Origin-Resource-Policy': 'cross-origin' }, config + STATIC('m.js'));
  }
  return send(res, 404, {});
}

function serveG(req, res, url) {
  if (url.pathname === '/__log') {
    const run = url.searchParams.get('run');
    return send(res, 200, { 'Content-Type': 'application/json' }, JSON.stringify(log.filter((entry) => entry.run === run)));
  }
  if (url.pathname === '/gtag/js') return send(res, 200, GTAG_HEADERS, STATIC('gtag-stand-in.js'));
  // C: the collect stand-in answers 204 without CORP (the real endpoint's headers were not read).
  if (url.pathname === '/g/collect') return send(res, 204, { 'Cache-Control': 'no-store' });
  return send(res, 404, {});
}

for (const [name, origin, handler] of [['A', A, serveA], ['M', M, serveM], ['G', G, serveG]]) {
  createServer(async (req, res) => {
    const url = new URL(req.url, origin);
    const text = await body(req);
    if (url.pathname !== '/__log' && url.pathname !== '/__ready') record(name === 'G' && url.pathname === '/g/collect' ? 'C' : name, req, url, text);
    handler(req, res, url);
  }).listen(port(origin), '::');
}
