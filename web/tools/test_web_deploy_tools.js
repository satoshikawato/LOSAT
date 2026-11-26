'use strict';
// Unit tests for LOSAT Web Cloudflare Workers Static Assets deployment tooling:
//   - web/deploy/cloudflare.config.ts
//   - web/tools/check_web_dist.mjs
//   - web/tools/check_deploy_headers.mjs
//   - .github/workflows/web-deploy.yml static safety invariants
//
// Run with: node --test web/tools/test_web_deploy_tools.js

const assert = require('node:assert/strict');
const { execFile, spawnSync } = require('node:child_process');
const { mkdtempSync, mkdirSync, readFileSync, rmSync, writeFileSync } = require('node:fs');
const { createServer } = require('node:http');
const { tmpdir } = require('node:os');
const { join, resolve } = require('node:path');
const test = require('node:test');
const { promisify } = require('node:util');

const execFileAsync = promisify(execFile);

const REPO_ROOT = resolve(__dirname, '..', '..');
const SOURCE_HEADERS = join(REPO_ROOT, 'web', 'app', 'public', '_headers');
const WORKFLOW_PATH = join(REPO_ROOT, '.github', 'workflows', 'web-deploy.yml');
const DEPLOY_CHECK_SCRIPT = join(__dirname, 'check_deploy_headers.mjs');
const SERIAL_WASM = 'losat-web-serial-0123456789abcdef.wasm';
const THREADS_WASM = 'losat-web-threads-fedcba9876543210.wasm';

const headersToolPromise = import('../app/build/headers.ts');
const distToolPromise = import('./check_web_dist.mjs');
const deployToolPromise = import('./check_deploy_headers.mjs');

function makeValidDist(baseDir, requiredFiles) {
  const distDir = join(baseDir, 'dist');
  mkdirSync(join(distDir, 'assets'), { recursive: true });
  mkdirSync(join(distDir, 'engine'), { recursive: true });

  const headersBytes = readFileSync(SOURCE_HEADERS);
  for (const rel of requiredFiles) {
    if (rel === '_headers') {
      writeFileSync(join(distDir, rel), headersBytes);
    } else if (rel === 'index.html') {
      writeFileSync(
        join(distDir, rel),
        '<!doctype html><html><body><div id="app"></div><script type="module" src="/assets/index-abc.js"></script></body></html>\n',
      );
    } else {
      writeFileSync(join(distDir, rel), `placeholder:${rel}\n`);
    }
  }

  const wasmHeader = Buffer.from([0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00]);
  writeFileSync(join(distDir, 'engine', SERIAL_WASM), wasmHeader);
  writeFileSync(join(distDir, 'engine', THREADS_WASM), wasmHeader);
  writeFileSync(
    join(distDir, 'assets', 'index-abc.js'),
    `export const ENGINE_ASSETS = { serial: { url: "/engine/${SERIAL_WASM}" }, threads: { url: "/engine/${THREADS_WASM}" } };\n`,
  );
  return distDir;
}

function withTempDir(fn) {
  const dir = mkdtempSync(join(tmpdir(), 'losat-web-deploy-test-'));
  try {
    return fn(dir);
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
}

test('cloudflare.config.ts configures preview and production modes safely', async () => {
  const mod = await import('../deploy/cloudflare.config.ts');
  const buildConfig = mod.default;
  const saved = {
    name: process.env.LOSAT_WEB_WORKER_NAME,
    domains: process.env.LOSAT_WEB_DOMAINS,
    workersDev: process.env.LOSAT_WEB_WORKERS_DEV,
  };
  try {
    delete process.env.LOSAT_WEB_WORKER_NAME;
    process.env.LOSAT_WEB_DOMAINS = 'losat.app';
    delete process.env.LOSAT_WEB_WORKERS_DEV;

    const preview = buildConfig({ isPreview: true });
    assert.equal(preview.worker.name, 'losat-web');
    assert.equal(preview.worker.previewUrls, true);
    assert.equal(preview.worker.workersDev, true);
    assert.equal('domains' in preview.worker, false);
    assert.deepEqual(preview.worker.assets, {
      htmlHandling: 'auto-trailing-slash',
      notFoundHandling: 'none',
    });

    const prodWithDomain = buildConfig({ isPreview: false });
    assert.deepEqual(prodWithDomain.worker.domains, ['losat.app']);
    assert.equal(prodWithDomain.worker.workersDev, false);
    assert.equal(prodWithDomain.worker.previewUrls, true);

    delete process.env.LOSAT_WEB_DOMAINS;
    process.env.LOSAT_WEB_WORKER_NAME = 'losat-web-custom';
    const prodNoDomain = buildConfig({ isPreview: false });
    assert.equal(prodNoDomain.worker.name, 'losat-web-custom');
    assert.equal('domains' in prodNoDomain.worker, false);
    assert.equal(prodNoDomain.worker.workersDev, true);
  } finally {
    for (const [key, envKey] of [
      ['name', 'LOSAT_WEB_WORKER_NAME'],
      ['domains', 'LOSAT_WEB_DOMAINS'],
      ['workersDev', 'LOSAT_WEB_WORKERS_DEV'],
    ]) {
      if (saved[key] === undefined) delete process.env[envKey];
      else process.env[envKey] = saved[key];
    }
  }
});

test('checkWebDist passes on a complete real-engine dist bundle', async () => {
  const { checkWebDist, REQUIRED_DIST_FILES } = await distToolPromise;
  withTempDir((dir) => {
    const distDir = makeValidDist(dir, REQUIRED_DIST_FILES);
    // Wasm binaries may hold static upstream version strings, which must not false-fail,
    // while actual caller-supplied build-machine paths are rejected.
    writeFileSync(
      join(distDir, 'engine', SERIAL_WASM),
      Buffer.concat([Buffer.from([0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00]), Buffer.from('/home/tomcat/static')]),
    );
    const result = checkWebDist(distDir, {
      headersPath: SOURCE_HEADERS,
      forbidPaths: ['/home/ci-runner/work/LOSAT'],
    });
    assert.equal(result.ok, true, result.errors.join('\n'));
    assert.equal(result.engineFiles.serial, `engine/${SERIAL_WASM}`);
    assert.equal(result.engineFiles.threads, `engine/${THREADS_WASM}`);
  });
});

test('checkWebDist rejects FakeEngine builds without real Wasm reactors', async () => {
  const { checkWebDist, REQUIRED_DIST_FILES } = await distToolPromise;
  withTempDir((dir) => {
    const distDir = makeValidDist(dir, REQUIRED_DIST_FILES);
    rmSync(join(distDir, 'engine'), { recursive: true, force: true });
    const result = checkWebDist(distDir, { headersPath: SOURCE_HEADERS });
    assert.equal(result.ok, false);
    assert.ok(result.errors.some((e) => e.includes('FakeEngine')));
  });
});

test('checkWebDist enforces single-file size limit and warns at >= 85%', async () => {
  const { checkWebDist, REQUIRED_DIST_FILES } = await distToolPromise;
  withTempDir((dir) => {
    const distDir = makeValidDist(dir, REQUIRED_DIST_FILES);
    writeFileSync(
      join(distDir, 'engine', THREADS_WASM),
      Buffer.concat([Buffer.from([0x00, 0x61, 0x73, 0x6d]), Buffer.alloc(900)]),
    );

    const warned = checkWebDist(distDir, { headersPath: SOURCE_HEADERS, maxBytes: 1000 });
    assert.equal(warned.ok, true, warned.errors.join('\n'));
    assert.ok(warned.warnings.some((w) => w.includes(THREADS_WASM)));

    const failed = checkWebDist(distDir, { headersPath: SOURCE_HEADERS, maxBytes: 500 });
    assert.equal(failed.ok, false);
    assert.ok(failed.errors.some((e) => e.includes('exceeding limit')));
  });
});

test('checkWebDist catches _headers drift and leaked build-machine paths in JS and Wasm', async () => {
  const { checkWebDist, REQUIRED_DIST_FILES } = await distToolPromise;
  withTempDir((dir) => {
    const distDir = makeValidDist(dir, REQUIRED_DIST_FILES);
    writeFileSync(join(distDir, '_headers'), '/*\n  X-Content-Type-Options: nosniff\n');
    writeFileSync(
      join(distDir, 'assets', 'index-abc.js'),
      `const leaked = "/home/runner/work/LOSAT/LOSAT/src/main.rs"; const a = "/engine/${SERIAL_WASM}"; const b = "/engine/${THREADS_WASM}";\n`,
    );
    writeFileSync(
      join(distDir, 'engine', THREADS_WASM),
      Buffer.concat([Buffer.from([0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00]), Buffer.from('/home/runner/work/LOSAT')]),
    );
    const result = checkWebDist(distDir, {
      headersPath: SOURCE_HEADERS,
      forbidPaths: ['/home/runner/work/LOSAT'],
    });
    assert.equal(result.ok, false);
    assert.ok(result.errors.some((e) => e.includes('dist/_headers differs')));
    assert.ok(result.errors.some((e) => e.includes('build-machine path leaked into assets/index-abc.js')));
    assert.ok(result.errors.some((e) => e.includes(`build-machine path leaked into engine/${THREADS_WASM}`)));
  });
});

test('checkWebDist allows RefSeq URL paths while rejecting leaked host paths', async () => {
  const { checkWebDist, REQUIRED_DIST_FILES } = await distToolPromise;
  withTempDir((dir) => {
    const distDir = makeValidDist(dir, REQUIRED_DIST_FILES);
    const refseqSentence = 'NCBI places no restrictions on the use of RefSeq (ncbi.nlm.nih.gov/home/about/policies).';
    writeFileSync(
      join(distDir, 'assets', 'index-abc.js'),
      `export const note = ${JSON.stringify(refseqSentence)}; export const ENGINE_ASSETS = { serial: "/engine/${SERIAL_WASM}", threads: "/engine/${THREADS_WASM}" };\n`,
    );
    const clean = checkWebDist(distDir, { headersPath: SOURCE_HEADERS });
    assert.equal(clean.ok, true, clean.errors.join('\n'));

    for (const leakedPath of [
      '/home/runner/work/x/',
      'file:///Users/alice/repo/',
      '/__w/LOSAT/LOSAT/',
      '/github/workspace/web/',
      'C:\\Users\\alice\\repo\\',
    ]) {
      writeFileSync(
        join(distDir, 'assets', 'index-abc.js'),
        `export const leaked = ${JSON.stringify(leakedPath)}; export const ENGINE_ASSETS = { serial: "/engine/${SERIAL_WASM}", threads: "/engine/${THREADS_WASM}" };\n`,
      );
      const leaked = checkWebDist(distDir, { headersPath: SOURCE_HEADERS });
      assert.equal(leaked.ok, false, `expected leak detection for ${leakedPath}`);
      assert.ok(
        leaked.errors.some((e) => e.includes('host path pattern leaked into assets/index-abc.js')),
        `missing leak error for ${leakedPath}: ${leaked.errors.join('\n')}`,
      );
    }
  });
});

test('checkWebDist reuses app/build/headers.ts and rejects duplicate _headers patterns', async () => {
  const { parseHeadersFile: appParse, pathHeaders: appPathHeaders } = await headersToolPromise;
  const { checkWebDist, parseHeadersFile, pathHeaders, REQUIRED_DIST_FILES } = await distToolPromise;
  assert.equal(parseHeadersFile, appParse);
  assert.equal(pathHeaders, appPathHeaders);

  withTempDir((dir) => {
    const distDir = makeValidDist(dir, REQUIRED_DIST_FILES);
    writeFileSync(join(distDir, '_headers'), `${readFileSync(SOURCE_HEADERS, 'utf8')}\n/*\n  X-Dup: 1\n`);
    const result = checkWebDist(distDir);
    assert.equal(result.ok, false);
    assert.ok(result.errors.some((e) => e.includes('duplicate _headers path rule')));
  });
});

test('preview JSON parsing and production URL resolution handle domain, workers.dev, and empty stdout', async () => {
  const { extractCfPreviewUrl, parseCfPreviewJson, resolveProductionUrl, resolveTargetUrlFromCli } = await deployToolPromise;

  assert.equal(
    extractCfPreviewUrl({
      preview_urls: ['https://preview-losat-web.example.workers.dev/'],
      deployment_urls: ['https://abc123-losat-web.example.workers.dev'],
    }),
    'https://preview-losat-web.example.workers.dev',
  );
  assert.equal(
    parseCfPreviewJson(JSON.stringify({ preview_urls: ['https://preview-losat-web.example.workers.dev/'] })),
    'https://preview-losat-web.example.workers.dev',
  );
  assert.throws(() => parseCfPreviewJson(''), /cf JSON output is empty/);
  assert.throws(() => parseCfPreviewJson('   \n'), /cf JSON output is empty/);
  assert.throws(() => parseCfPreviewJson('{bad json'), /not valid JSON/);

  assert.equal(resolveProductionUrl({ domains: 'losat.app' }), 'https://losat.app');
  assert.equal(resolveProductionUrl({ domains: ' https://losat.app/ , www.losat.app ' }), 'https://losat.app');
  assert.equal(
    resolveProductionUrl({ domains: '', workersDevUrl: 'https://losat-web.example.workers.dev/' }),
    'https://losat-web.example.workers.dev',
  );
  assert.throws(() => resolveProductionUrl({ domains: '', workersDevUrl: '' }), /cannot determine production URL/);

  withTempDir((dir) => {
    const emptyJson = join(dir, 'empty.json');
    writeFileSync(emptyJson, '');
    assert.throws(
      () => resolveTargetUrlFromCli({ cfJson: emptyJson, target: 'preview' }),
      /cf JSON output is empty/,
    );
    assert.throws(
      () => resolveTargetUrlFromCli({ cfJson: emptyJson, target: 'production' }),
      /--cf-json is only supported for preview/,
    );

    const cliEmpty = spawnSync(process.execPath, [DEPLOY_CHECK_SCRIPT, '--cf-json', emptyJson], {
      encoding: 'utf8',
    });
    assert.equal(cliEmpty.status, 1);
    assert.match(cliEmpty.stderr, /cf JSON output is empty/);
    assert.doesNotMatch(cliEmpty.stderr, /SyntaxError/);

    const cliProdCfJson = spawnSync(
      process.execPath,
      [DEPLOY_CHECK_SCRIPT, '--cf-json', emptyJson, '--target', 'production', '--domains', 'losat.app'],
      { encoding: 'utf8' },
    );
    assert.equal(cliProdCfJson.status, 1);
    assert.match(cliProdCfJson.stderr, /--cf-json is only supported for preview/);
  });
});

async function withMockSite(options, fn) {
  const { parseHeadersFile, pathHeaders } = await headersToolPromise;
  const rules = parseHeadersFile(readFileSync(SOURCE_HEADERS, 'utf8'));
  const routes = new Map([
    [
      '/',
      {
        body: '<!doctype html><html><body><div id="app"></div><script type="module" src="/assets/index-abc.js"></script></body></html>',
        contentType: 'text/html; charset=utf-8',
        etag: '"root-v1"',
      },
    ],
    [
      '/assets/index-abc.js',
      {
        body: `export const assets = ["/engine/${SERIAL_WASM}", "/engine/${THREADS_WASM}"];`,
        contentType: 'application/javascript; charset=utf-8',
        etag: '"js-v1"',
      },
    ],
    [
      `/engine/${SERIAL_WASM}`,
      { body: Buffer.from([0x00, 0x61, 0x73, 0x6d]), contentType: 'application/wasm', etag: '"serial-v1"' },
    ],
    [
      `/engine/${THREADS_WASM}`,
      { body: Buffer.from([0x00, 0x61, 0x73, 0x6d]), contentType: 'application/wasm', etag: '"threads-v1"' },
    ],
    ['/manifest.webmanifest', { body: '{}', contentType: 'application/manifest+json', etag: '"mf-v1"' }],
    ['/favicon.svg', { body: '<svg/>', contentType: 'image/svg+xml', etag: '"fav-v1"' }],
    ['/LICENSE.txt', { body: 'MIT', contentType: 'text/plain', etag: '"lic-v1"' }],
    ['/THIRD_PARTY_NOTICES.txt', { body: 'notices', contentType: 'text/plain', etag: '"not-v1"' }],
    ['/THIRD_PARTY_LICENSES.txt', { body: 'licenses', contentType: 'text/plain', etag: '"tpl-v1"' }],
    ['/THIRD_PARTY_LICENSES_NPM.txt', { body: 'npm', contentType: 'text/plain', etag: '"npm-v1"' }],
  ]);

  const server = createServer((req, res) => {
    const reqPath = new URL(req.url ?? '/', 'http://127.0.0.1').pathname;
    const entry = routes.get(reqPath);
    if (!entry) {
      res.writeHead(options.missingStatus ?? 404, { 'Content-Type': 'text/plain' });
      res.end('Not found');
      return;
    }

    const headers = { ...pathHeaders(reqPath, rules), 'Content-Type': entry.contentType, ETag: entry.etag };
    if (options.sendNoindex) headers['X-Robots-Tag'] = 'noindex';
    if (options.dropWasmImmutable && reqPath.startsWith('/engine/')) delete headers['Cache-Control'];

    if (req.headers['if-none-match'] === entry.etag) {
      const h304 = { ...headers };
      delete h304['Content-Type'];
      if (options.dropCoepOn304) delete h304['Cross-Origin-Embedder-Policy'];
      res.writeHead(304, h304);
      res.end();
      return;
    }

    res.writeHead(200, headers);
    if (req.method === 'HEAD') {
      res.end();
    } else {
      res.end(entry.body);
    }
  });

  await new Promise((done) => server.listen(0, '127.0.0.1', done));
  const addr = server.address();
  const baseUrl = `http://127.0.0.1:${addr.port}`;
  try {
    return await fn(baseUrl);
  } finally {
    await new Promise((done) => server.close(done));
  }
}

test('checkDeployedSite passes against a compliant server (including 304 and noindex)', async () => {
  const { checkDeployedSite } = await deployToolPromise;
  await withMockSite({ sendNoindex: true }, async (url) => {
    const result = await checkDeployedSite(url, {
      headersPath: SOURCE_HEADERS,
      expectNoindex: true,
      retries: 1,
      retryDelayMs: 10,
    });
    assert.equal(result.ok, true, result.errors.join('\n'));
  });
});

test('checkDeployedSite catches missing COEP on 304, missing noindex, missing immutable, and non-404 fallback', async () => {
  const { checkDeployedSite } = await deployToolPromise;
  await withMockSite({ sendNoindex: false, dropCoepOn304: true, dropWasmImmutable: true, missingStatus: 200 }, async (url) => {
    const result = await checkDeployedSite(url, {
      headersPath: SOURCE_HEADERS,
      expectNoindex: true,
      retries: 1,
      retryDelayMs: 10,
    });
    assert.equal(result.ok, false);
    assert.ok(result.errors.some((e) => e.includes('(304)') && e.includes('Cross-Origin-Embedder-Policy')));
    assert.ok(result.errors.some((e) => e.includes('X-Robots-Tag')));
    assert.ok(result.errors.some((e) => e.includes('Cache-Control')));
    assert.ok(result.errors.some((e) => e.includes('expected HTTP 404')));
  });
});

test('checkDeployedSite CLI passes in production --url mode without --cf-json', async () => {
  await withMockSite({ sendNoindex: false }, async (url) => {
    const dir = mkdtempSync(join(tmpdir(), 'losat-web-deploy-cli-'));
    try {
      const urlFile = join(dir, 'deploy-url.txt');
      await execFileAsync(
        process.execPath,
        [
          DEPLOY_CHECK_SCRIPT,
          '--url',
          url,
          '--headers',
          SOURCE_HEADERS,
          '--retries',
          '1',
          '--write-url-file',
          urlFile,
        ],
        { encoding: 'utf8' },
      );
      assert.equal(readFileSync(urlFile, 'utf8').trim(), url);
    } finally {
      rmSync(dir, { recursive: true, force: true });
    }
  });
});

test('web-deploy.yml follows workflow security and pinning invariants', () => {
  const workflow = readFileSync(WORKFLOW_PATH, 'utf8');
  assert.doesNotMatch(workflow, /^\s*ref:\s*\$\{\{\s*inputs\./m, 'must never check out an input ref');
  assert.match(workflow, /^permissions:\s*\n\s+contents:\s*read\b/m, 'top-level permissions must be contents: read');
  for (const match of workflow.matchAll(/^\s*-?\s*uses:\s*([^\s#]+)/gm)) {
    const ref = match[1];
    if (ref.startsWith('./')) continue;
    assert.match(ref, /@[0-9a-f]{40}$/, `action must be pinned by full 40-char commit SHA: ${ref}`);
  }
  assert.doesNotMatch(workflow, /\brun:[^\n]*\$\{\{\s*inputs\./, 'inputs must pass via env:, not inline in run:');
  assert.doesNotMatch(
    workflow,
    /deploy\s+\.\.\/app\/dist\s*>/,
    'production cf deploy writes no JSON to stdout and must not be redirected to a JSON file',
  );
  assert.match(workflow, /--url\s+"\$PROD_URL"/, 'production post-deploy check must pass --url directly');
});
