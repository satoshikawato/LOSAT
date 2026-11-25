#!/usr/bin/env node
// Post-deploy HTTP header and 304 revalidation checks for LOSAT Web on Cloudflare Workers Static Assets.
//
// Reads `web/app/public/_headers` as the single source of truth for expected response headers and checks:
//   1. `GET /` returns 200 with every `/*` security/isolation header and `<div id="app">`.
//   2. Optional `--expect-noindex` confirms `X-Robots-Tag` contains `noindex` (preview URLs).
//   3. Conditional `If-None-Match` requests on `/`, the main `/assets/*.js` bundle, and `/engine/*.wasm`
//      return `304 Not Modified` while still carrying COOP, COEP, CORP, and CSP (WebKit/Safari blocks
//      restarting the Engine worker when a 304 omits COEP/CORP).
//   4. Both `/engine/losat-web-{serial,threads}-*.wasm` modules return 200 with `application/wasm`,
//      immutable `Cache-Control`, and all `/*` headers.
//   5. Required public metadata and license files return 200, and an unknown path returns 404.

import { existsSync, readdirSync, readFileSync, writeFileSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';
import { parseHeadersFile, pathHeaders } from './check_web_dist.mjs';

const REVALIDATION_REQUIRED_HEADERS = [
  'Cross-Origin-Opener-Policy',
  'Cross-Origin-Embedder-Policy',
  'Cross-Origin-Resource-Policy',
  'Content-Security-Policy',
];

const PUBLIC_FILES_TO_CHECK = [
  '/manifest.webmanifest',
  '/favicon.svg',
  '/LICENSE.txt',
  '/THIRD_PARTY_NOTICES.txt',
  '/THIRD_PARTY_LICENSES.txt',
  '/THIRD_PARTY_LICENSES_NPM.txt',
];

const ENGINE_FILENAME_RE = /engine\/losat-web-(?:serial|threads)-[0-9a-f]{16}\.wasm/g;
const SCRIPT_SRC_RE = /<script\b[^>]*\bsrc=["']([^"']+\.js)["']/gi;

function firstHttpsUrl(value) {
  if (typeof value === 'string') {
    const trimmed = value.trim();
    return /^https?:\/\/\S+$/.test(trimmed) ? trimmed : undefined;
  }
  if (Array.isArray(value)) {
    for (const item of value) {
      const found = firstHttpsUrl(item);
      if (found) return found;
    }
    return undefined;
  }
  if (value && typeof value === 'object') {
    for (const item of Object.values(value)) {
      const found = firstHttpsUrl(item);
      if (found) return found;
    }
  }
  return undefined;
}

/**
 * Extracts the deployed base URL from `cf previews deploy` or `cf deploy` JSON stdout.
 * When `target === 'production'` and `domains` lists custom domain(s), uses `https://<first domain>`.
 */
export function extractCfDeployUrl(payload, options = {}) {
  const target = options.target ?? 'preview';
  const domains = (options.domains ?? '')
    .split(',')
    .map((d) => d.trim())
    .filter((d) => d.length > 0);

  if (target === 'production' && domains.length > 0) {
    const first = domains[0].replace(/^https?:\/\//, '').replace(/\/+$/, '');
    return `https://${first}`;
  }

  const root = payload && typeof payload === 'object' ? payload : {};
  const candidates =
    target === 'preview'
      ? [
          root.preview_urls,
          root.deployment_urls,
          root.result?.preview_urls,
          root.result?.deployment_urls,
          root.url,
          root.urls,
          root,
        ]
      : [
          root.urls,
          root.url,
          root.deployment_urls,
          root.worker_url,
          root.result?.urls,
          root.result?.url,
          root.result?.deployment_urls,
          root,
        ];

  for (const candidate of candidates) {
    const found = firstHttpsUrl(candidate);
    if (found) return found.replace(/\/+$/, '');
  }
  throw new Error(`could not extract deployed URL from cf JSON output for target=${target}`);
}

function assertHeadersMatch(sitePath, actualHeaders, expectedMap, errors, label = sitePath) {
  for (const [name, expectedValue] of Object.entries(expectedMap)) {
    const actual = actualHeaders.get(name);
    if (actual === null) {
      errors.push(`${label}: missing header ${name}`);
    } else if (actual.trim() !== expectedValue.trim()) {
      errors.push(`${label}: header ${name} expected ${JSON.stringify(expectedValue)}, got ${JSON.stringify(actual)}`);
    }
  }
}

function discoverDistEnginePaths(distDir) {
  if (!distDir) return [];
  const engineDir = join(distDir, 'engine');
  if (!existsSync(engineDir)) return [];
  return readdirSync(engineDir)
    .filter((name) => /^losat-web-(?:serial|threads)-[0-9a-f]{16}\.wasm$/.test(name))
    .sort()
    .map((name) => `/engine/${name}`);
}

function discoverScriptPathsFromHtml(html) {
  const paths = [];
  const seen = new Set();
  for (const match of html.matchAll(SCRIPT_SRC_RE)) {
    const src = match[1];
    if (!src || src.startsWith('http://') || src.startsWith('https://')) continue;
    const normalized = src.startsWith('/') ? src : `/${src}`;
    if (!seen.has(normalized)) {
      seen.add(normalized);
      paths.push(normalized);
    }
  }
  return paths;
}

async function fetchWithTimeout(fetchImpl, url, init, timeoutMs) {
  const controller = new AbortController();
  const timer = setTimeout(() => controller.abort(), timeoutMs);
  try {
    return await fetchImpl(url, { ...init, redirect: 'manual', signal: controller.signal });
  } finally {
    clearTimeout(timer);
  }
}

async function waitForRoot(fetchImpl, baseUrl, retries, retryDelayMs, timeoutMs) {
  let lastError;
  for (let attempt = 1; attempt <= retries; attempt += 1) {
    try {
      const response = await fetchWithTimeout(fetchImpl, `${baseUrl}/`, { method: 'GET' }, timeoutMs);
      if (response.status === 200) return response;
      lastError = new Error(`GET / returned HTTP ${response.status}`);
    } catch (error) {
      lastError = error;
    }
    if (attempt < retries) {
      await new Promise((done) => setTimeout(done, retryDelayMs * attempt));
    }
  }
  throw lastError ?? new Error(`GET ${baseUrl}/ failed after ${retries} attempts`);
}

async function checkRevalidation304(fetchImpl, baseUrl, sitePath, etag, expectedHeaders, errors, timeoutMs) {
  if (!etag) {
    errors.push(`${sitePath}: missing ETag header on 200 response (cannot test 304 revalidation)`);
    return;
  }
  const response = await fetchWithTimeout(
    fetchImpl,
    `${baseUrl}${sitePath}`,
    { method: 'GET', headers: { 'If-None-Match': etag } },
    timeoutMs,
  );
  if (response.status !== 304) {
    errors.push(`${sitePath} (If-None-Match): expected HTTP 304, got ${response.status}`);
    return;
  }
  const subset = {};
  for (const name of REVALIDATION_REQUIRED_HEADERS) {
    if (expectedHeaders[name] !== undefined) subset[name] = expectedHeaders[name];
  }
  assertHeadersMatch(sitePath, response.headers, subset, errors, `${sitePath} (304)`);
}

/**
 * Checks a live deployed LOSAT Web URL against `web/app/public/_headers`.
 * Returns `{ ok, url, errors, checkedPaths }`.
 */
export async function checkDeployedSite(baseUrl, options = {}) {
  const normalizedBase = baseUrl.replace(/\/+$/, '');
  const fetchImpl = options.fetchImpl ?? globalThis.fetch;
  const retries = options.retries ?? 6;
  const retryDelayMs = options.retryDelayMs ?? 2000;
  const timeoutMs = options.timeoutMs ?? 15_000;
  const expectNoindex = Boolean(options.expectNoindex);
  const headersPath = options.headersPath ?? resolve('web/app/public/_headers');
  const distDir = options.distDir;

  const { rules } = parseHeadersFile(readFileSync(headersPath, 'utf8'));
  const errors = [];
  const checkedPaths = [];

  // 1. Root document `/` (200 + headers + optional noindex + 304).
  let rootResponse;
  try {
    rootResponse = await waitForRoot(fetchImpl, normalizedBase, retries, retryDelayMs, timeoutMs);
  } catch (error) {
    return {
      ok: false,
      url: normalizedBase,
      errors: [`root reachability failed: ${error instanceof Error ? error.message : String(error)}`],
      checkedPaths,
    };
  }

  checkedPaths.push('/');
  const rootExpected = pathHeaders('/', rules);
  assertHeadersMatch('/', rootResponse.headers, rootExpected, errors);
  if (expectNoindex) {
    const robots = rootResponse.headers.get('X-Robots-Tag') ?? '';
    if (!/\bnoindex\b/i.test(robots)) {
      errors.push(`/: expected X-Robots-Tag containing "noindex" on preview URL, got ${JSON.stringify(robots)}`);
    }
  }
  const rootHtml = await rootResponse.text();
  if (!rootHtml.includes('<div id="app">')) {
    errors.push('/: HTML body does not contain <div id="app">');
  }
  await checkRevalidation304(
    fetchImpl,
    normalizedBase,
    '/',
    rootResponse.headers.get('ETag'),
    rootExpected,
    errors,
    timeoutMs,
  );

  // 2. Main script bundle (`/assets/*.js`: 200 + headers + 304 + fallback Wasm discovery).
  const scriptPaths = discoverScriptPathsFromHtml(rootHtml);
  let enginePaths = discoverDistEnginePaths(distDir);
  if (scriptPaths.length === 0) {
    errors.push('/: no module <script src="..."> asset found in HTML');
  } else {
    const mainScript = scriptPaths[0];
    checkedPaths.push(mainScript);
    const scriptResp = await fetchWithTimeout(fetchImpl, `${normalizedBase}${mainScript}`, { method: 'GET' }, timeoutMs);
    if (scriptResp.status !== 200) {
      errors.push(`${mainScript}: expected HTTP 200, got ${scriptResp.status}`);
    } else {
      const scriptExpected = pathHeaders(mainScript, rules);
      assertHeadersMatch(mainScript, scriptResp.headers, scriptExpected, errors);
      const scriptText = await scriptResp.text();
      if (enginePaths.length === 0) {
        const discovered = new Set();
        for (const match of scriptText.matchAll(ENGINE_FILENAME_RE)) {
          discovered.add(`/${match[0]}`);
        }
        enginePaths = [...discovered].sort();
      }
      await checkRevalidation304(
        fetchImpl,
        normalizedBase,
        mainScript,
        scriptResp.headers.get('ETag'),
        scriptExpected,
        errors,
        timeoutMs,
      );
    }
  }

  // 3. Engine Wasm modules (`/engine/losat-web-{serial,threads}-*.wasm`: 200 + MIME + immutable + 304).
  if (enginePaths.length < 2) {
    errors.push(`expected 2 engine Wasm paths (serial and threads), found ${enginePaths.length}`);
  }
  for (const enginePath of enginePaths) {
    checkedPaths.push(enginePath);
    const headResp = await fetchWithTimeout(fetchImpl, `${normalizedBase}${enginePath}`, { method: 'HEAD' }, timeoutMs);
    if (headResp.status !== 200) {
      errors.push(`${enginePath} (HEAD): expected HTTP 200, got ${headResp.status}`);
      continue;
    }
    const contentType = headResp.headers.get('Content-Type') ?? '';
    if (!contentType.toLowerCase().startsWith('application/wasm')) {
      errors.push(`${enginePath}: expected Content-Type application/wasm, got ${JSON.stringify(contentType)}`);
    }
    const engineExpected = pathHeaders(enginePath, rules);
    assertHeadersMatch(enginePath, headResp.headers, engineExpected, errors);
    const etag = headResp.headers.get('ETag');
    if (etag) {
      await checkRevalidation304(fetchImpl, normalizedBase, enginePath, etag, engineExpected, errors, timeoutMs);
    }
  }

  // 4. Required public metadata & license files.
  for (const publicPath of PUBLIC_FILES_TO_CHECK) {
    checkedPaths.push(publicPath);
    const resp = await fetchWithTimeout(fetchImpl, `${normalizedBase}${publicPath}`, { method: 'GET' }, timeoutMs);
    if (resp.status !== 200) {
      errors.push(`${publicPath}: expected HTTP 200, got ${resp.status}`);
      continue;
    }
    assertHeadersMatch(publicPath, resp.headers, pathHeaders(publicPath, rules), errors);
    await resp.arrayBuffer();
  }

  // 5. Missing route returns 404 (confirms `notFoundHandling: "none"`).
  const missingPath = '/__losat_missing_probe__';
  checkedPaths.push(missingPath);
  const missingResp = await fetchWithTimeout(fetchImpl, `${normalizedBase}${missingPath}`, { method: 'GET' }, timeoutMs);
  if (missingResp.status !== 404) {
    errors.push(`${missingPath}: expected HTTP 404, got ${missingResp.status}`);
  }

  return { ok: errors.length === 0, url: normalizedBase, errors, checkedPaths };
}

function parseCliArgs(argv) {
  const args = {
    url: undefined,
    cfJson: undefined,
    target: 'preview',
    domains: '',
    dist: undefined,
    headers: 'web/app/public/_headers',
    expectNoindex: false,
    retries: 6,
    timeoutMs: 15_000,
    writeUrlFile: undefined,
  };
  for (let i = 0; i < argv.length; i += 1) {
    const flag = argv[i];
    if (flag === '--url') args.url = argv[++i];
    else if (flag === '--cf-json') args.cfJson = argv[++i];
    else if (flag === '--target') args.target = argv[++i];
    else if (flag === '--domains') args.domains = argv[++i];
    else if (flag === '--dist') args.dist = argv[++i];
    else if (flag === '--headers') args.headers = argv[++i];
    else if (flag === '--expect-noindex') args.expectNoindex = true;
    else if (flag === '--retries') args.retries = Number(argv[++i]);
    else if (flag === '--timeout-ms') args.timeoutMs = Number(argv[++i]);
    else if (flag === '--write-url-file') args.writeUrlFile = argv[++i];
    else if (flag === '--help' || flag === '-h') {
      console.log(
        'Usage: node web/tools/check_deploy_headers.mjs (--url <url> | --cf-json <path>) [--target preview|production] [--domains <csv>] [--dist <dir>] [--headers <path>] [--expect-noindex] [--write-url-file <path>]',
      );
      process.exit(0);
    } else {
      throw new Error(`unknown argument: ${flag}`);
    }
  }
  return args;
}

if (process.argv[1] && import.meta.url === pathToFileURL(resolve(process.argv[1])).href) {
  const cli = parseCliArgs(process.argv.slice(2));
  let targetUrl = cli.url;
  if (!targetUrl && cli.cfJson) {
    const raw = readFileSync(resolve(cli.cfJson), 'utf8');
    const parsed = JSON.parse(raw);
    targetUrl = extractCfDeployUrl(parsed, { target: cli.target, domains: cli.domains });
  }
  if (!targetUrl) {
    console.error('ERROR: provide --url <url> or --cf-json <path>');
    process.exit(1);
  }
  if (cli.writeUrlFile) {
    writeFileSync(resolve(cli.writeUrlFile), `${targetUrl}\n`, 'utf8');
  }

  const result = await checkDeployedSite(targetUrl, {
    headersPath: resolve(cli.headers),
    distDir: cli.dist ? resolve(cli.dist) : undefined,
    expectNoindex: cli.expectNoindex,
    retries: cli.retries,
    timeoutMs: cli.timeoutMs,
  });

  console.log(`Checked ${result.checkedPaths.length} endpoints on ${result.url}`);
  if (!result.ok) {
    for (const error of result.errors) {
      console.error(`ERROR: ${error}`);
    }
    process.exit(1);
  }
  console.log(`OK: all headers, 304 revalidation, Wasm assets, and 404 checks passed on ${result.url}`);
}
