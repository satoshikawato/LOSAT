#!/usr/bin/env node
// Pre-deploy bundle checks for LOSAT Web (`web/app/dist/` -> Cloudflare Workers Static Assets).
//
// Checks performed before uploading static files:
//   1. Required site assets and license notices exist and are non-empty.
//   2. Both real Wasm reactors (`losat-web-serial-<hash>.wasm`, `losat-web-threads-<hash>.wasm`)
//      are present with valid Wasm headers and referenced by the built JS (no FakeEngine build).
//   3. Every file is under Cloudflare's 25 MiB static-asset limit (warn at >= 85%) and total file
//      count is under the Free-tier 20,000-file limit.
//   4. `dist/_headers` matches `web/app/public/_headers` byte for byte and stays within
//      Cloudflare's `_headers` rule and line-length limits.
//   5. No build-machine path leaks into shipped JS, CSS, HTML, Wasm, or manifest files.

import { existsSync, readdirSync, readFileSync, statSync } from 'node:fs';
import { extname, join, posix, relative, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';
import { parseHeadersFile, pathHeaders } from '../app/build/headers.ts';

export { parseHeadersFile, pathHeaders };

/** Cloudflare Workers Static Assets per-file size limit (25 MiB). */
export const CLOUDFLARE_ASSET_MAX_BYTES = 25 * 1024 * 1024;
/** Warning threshold (85% of 25 MiB). */
export const CLOUDFLARE_ASSET_WARN_BYTES = Math.floor(CLOUDFLARE_ASSET_MAX_BYTES * 0.85);
/** Cloudflare Workers Static Assets file-count limit on the Free plan. */
export const CLOUDFLARE_ASSET_MAX_FILES = 20_000;
/** Cloudflare `_headers` limits. */
export const CLOUDFLARE_HEADERS_MAX_RULES = 100;
export const CLOUDFLARE_HEADERS_MAX_LINE_CHARS = 2_000;

export const REQUIRED_DIST_FILES = [
  'index.html',
  '_headers',
  'manifest.webmanifest',
  'favicon.svg',
  'favicon-32.png',
  'apple-touch-icon.png',
  'icon-192.png',
  'icon-512.png',
  'icon-maskable-512.png',
  'LICENSE.txt',
  'THIRD_PARTY_NOTICES.txt',
  'THIRD_PARTY_LICENSES.txt',
  'THIRD_PARTY_LICENSES_NPM.txt',
];

const WASM_MAGIC = Buffer.from([0x00, 0x61, 0x73, 0x6d]);
const SERIAL_RE = /^losat-web-serial-[0-9a-f]{16}\.wasm$/;
const THREADS_RE = /^losat-web-threads-[0-9a-f]{16}\.wasm$/;
const PRIVACY_SCAN_EXTS = new Set(['.js', '.mjs', '.css', '.html', '.wasm', '.webmanifest']);
const GENERIC_HOST_PATH_PATTERNS = [
  /(?<![A-Za-z0-9._~%-])\/home\/[A-Za-z0-9._-]+\//,
  /(?<![A-Za-z0-9._~%-])\/Users\/[A-Za-z0-9._-]+\//,
  /(?<![A-Za-z0-9._~%-])\/__w\/[A-Za-z0-9._-]+\//,
  /(?<![A-Za-z0-9._~%-])\/github\/workspace\//,
  /(?<![A-Za-z0-9._~%-])[A-Za-z]:\\{1,2}Users\\{1,2}[A-Za-z0-9._-]+\\{1,2}/,
];

function walkFiles(rootDir, currentDir = rootDir, out = []) {
  const entries = readdirSync(currentDir, { withFileTypes: true }).sort((a, b) => a.name.localeCompare(b.name));
  for (const entry of entries) {
    const abs = join(currentDir, entry.name);
    if (entry.isDirectory()) {
      walkFiles(rootDir, abs, out);
    } else if (entry.isFile()) {
      const rel = relative(rootDir, abs).split('\\').join(posix.sep);
      out.push({ abs, rel, size: statSync(abs).size });
    }
  }
  return out;
}

function normalizeForbidPaths(paths) {
  const out = [];
  const seen = new Set();
  for (const raw of paths ?? []) {
    const cleaned = String(raw ?? '').trim().replace(/[\\/]+$/, '');
    // Ignore trivial short roots such as `/root` or `/` to avoid binary false positives.
    if (cleaned.length < 8 || seen.has(cleaned)) continue;
    seen.add(cleaned);
    out.push(cleaned);
  }
  return out;
}

function formatMiB(bytes) {
  return `${(bytes / (1024 * 1024)).toFixed(2)} MiB`;
}

/**
 * Checks a built `web/app/dist/` directory for Cloudflare Workers Static Assets deployment.
 * Returns `{ ok, errors, warnings, filesCount, largestFiles, engineFiles }`.
 */
export function checkWebDist(distDir, options = {}) {
  const maxBytes = options.maxBytes ?? CLOUDFLARE_ASSET_MAX_BYTES;
  const warnBytes = options.warnBytes ?? Math.floor(maxBytes * 0.85);
  const maxFiles = options.maxFiles ?? CLOUDFLARE_ASSET_MAX_FILES;
  const forbidPaths = normalizeForbidPaths(options.forbidPaths);
  const sourceHeadersPath = options.headersPath;

  const errors = [];
  const warnings = [];

  if (!existsSync(distDir) || !statSync(distDir).isDirectory()) {
    return {
      ok: false,
      errors: [`dist directory does not exist: ${distDir}`],
      warnings,
      filesCount: 0,
      largestFiles: [],
      engineFiles: { serial: null, threads: null },
    };
  }

  const files = walkFiles(distDir);
  if (files.length === 0) {
    errors.push(`dist directory is empty: ${distDir}`);
  }
  if (files.length > maxFiles) {
    errors.push(`dist has ${files.length} files, exceeding the ${maxFiles}-file limit`);
  }

  // 1. Required files present and non-empty.
  for (const rel of REQUIRED_DIST_FILES) {
    const abs = join(distDir, rel);
    if (!existsSync(abs) || !statSync(abs).isFile()) {
      errors.push(`missing required file in dist: ${rel}`);
    } else if (statSync(abs).size === 0) {
      errors.push(`required file in dist is empty: ${rel}`);
    }
  }

  // 2. Real engine Wasm modules present (reject FakeEngine builds).
  const engineDir = join(distDir, 'engine');
  let serialFile = null;
  let threadsFile = null;
  if (!existsSync(engineDir) || !statSync(engineDir).isDirectory()) {
    errors.push('missing dist/engine/ directory (build likely used FakeEngine instead of LOSAT_WEB_REACTORS)');
  } else {
    const engineEntries = readdirSync(engineDir).sort();
    const serials = engineEntries.filter((name) => SERIAL_RE.test(name));
    const threads = engineEntries.filter((name) => THREADS_RE.test(name));
    if (serials.length !== 1) {
      errors.push(`expected exactly 1 serial Wasm reactor in dist/engine/, found ${serials.length}`);
    } else {
      serialFile = `engine/${serials[0]}`;
    }
    if (threads.length !== 1) {
      errors.push(`expected exactly 1 threaded Wasm reactor in dist/engine/, found ${threads.length}`);
    } else {
      threadsFile = `engine/${threads[0]}`;
    }

    for (const rel of [serialFile, threadsFile]) {
      if (!rel) continue;
      const bytes = readFileSync(join(distDir, rel));
      if (bytes.length < 8 || !bytes.subarray(0, 4).equals(WASM_MAGIC)) {
        errors.push(`${rel} does not start with the WebAssembly magic header`);
      }
    }
  }

  // Confirm built JS references both emitted Wasm filenames (not ENGINE_ASSETS = null).
  const jsFiles = files.filter((f) => f.rel.startsWith('assets/') && (f.rel.endsWith('.js') || f.rel.endsWith('.mjs')));
  if (jsFiles.length === 0) {
    errors.push('no JavaScript bundles found under dist/assets/');
  } else if (serialFile && threadsFile) {
    const combinedJs = jsFiles.map((f) => readFileSync(f.abs, 'utf8')).join('\n');
    if (!combinedJs.includes(serialFile)) {
      errors.push(`built JS under dist/assets/ does not reference ${serialFile} (FakeEngine build?)`);
    }
    if (!combinedJs.includes(threadsFile)) {
      errors.push(`built JS under dist/assets/ does not reference ${threadsFile} (FakeEngine build?)`);
    }
  }

  // 3. Per-file size limit (25 MiB on Cloudflare Workers Static Assets).
  const largestFiles = [...files].sort((a, b) => b.size - a.size).slice(0, 5);
  for (const file of files) {
    if (file.size > maxBytes) {
      errors.push(`${file.rel} is ${file.size} bytes (${formatMiB(file.size)}), exceeding limit ${maxBytes} bytes (${formatMiB(maxBytes)})`);
    } else if (file.size >= warnBytes) {
      const pct = ((file.size / maxBytes) * 100).toFixed(1);
      warnings.push(`${file.rel} is ${formatMiB(file.size)} (${pct}% of the ${formatMiB(maxBytes)} single-file limit)`);
    }
  }

  // 4. `_headers` exact match and Cloudflare limits.
  const distHeadersPath = join(distDir, '_headers');
  if (existsSync(distHeadersPath)) {
    const distHeadersText = readFileSync(distHeadersPath, 'utf8');
    if (sourceHeadersPath) {
      if (!existsSync(sourceHeadersPath)) {
        errors.push(`source _headers file not found: ${sourceHeadersPath}`);
      } else {
        const sourceText = readFileSync(sourceHeadersPath, 'utf8');
        if (distHeadersText !== sourceText) {
          errors.push(`dist/_headers differs from ${sourceHeadersPath}`);
        }
      }
    }

    const lines = distHeadersText.split(/\r?\n/);
    for (let i = 0; i < lines.length; i += 1) {
      if (lines[i].length > CLOUDFLARE_HEADERS_MAX_LINE_CHARS) {
        errors.push(`dist/_headers line ${i + 1} has ${lines[i].length} chars, exceeding Cloudflare's ${CLOUDFLARE_HEADERS_MAX_LINE_CHARS}-char limit`);
      }
    }

    try {
      const rules = parseHeadersFile(distHeadersText);
      if (rules.size > CLOUDFLARE_HEADERS_MAX_RULES) {
        errors.push(`dist/_headers has ${rules.size} rules, exceeding Cloudflare's ${CLOUDFLARE_HEADERS_MAX_RULES}-rule limit`);
      }
      const rootHeaders = pathHeaders('/', rules);
      if (rootHeaders['Cross-Origin-Opener-Policy'] !== 'same-origin') {
        errors.push('dist/_headers must set Cross-Origin-Opener-Policy: same-origin on /*');
      }
      if (rootHeaders['Cross-Origin-Embedder-Policy'] !== 'require-corp') {
        errors.push('dist/_headers must set Cross-Origin-Embedder-Policy: require-corp on /*');
      }
      if (rootHeaders['Cross-Origin-Resource-Policy'] !== 'same-origin') {
        errors.push('dist/_headers must set Cross-Origin-Resource-Policy: same-origin on /*');
      }
      if (!rootHeaders['Content-Security-Policy']) {
        errors.push('dist/_headers must set Content-Security-Policy on /*');
      }
      const engineHeaders = pathHeaders('/engine/losat-web-threads-0000000000000000.wasm', rules);
      if (!engineHeaders['Cache-Control'] || !engineHeaders['Cache-Control'].includes('immutable')) {
        errors.push('dist/_headers must set an immutable Cache-Control rule on /engine/*');
      }
    } catch (error) {
      errors.push(`dist/_headers failed to parse: ${error instanceof Error ? error.message : String(error)}`);
    }
  }

  // 5. Build-machine path leak scan (TD-11 / V-PRIV).
  // Exact caller-supplied host paths (`--forbid-path`: workspace, HOME, CARGO_HOME, RUNNER_TEMP)
  // are rejected in every scanned file including `.wasm`. Generic `/home/<user>/` regexes apply
  // only to web text assets (`.js`, `.css`, `.html`, `.webmanifest`), because the Wasm engine
  // embeds NCBI's static `rpsbproc` version strings that contain `/home/tomcat/` and `/home/coremake/`.
  const forbidBuffers = forbidPaths.map((p) => ({ text: p, buf: Buffer.from(p, 'utf8') }));
  for (const file of files) {
    const ext = extname(file.rel).toLowerCase();
    if (!PRIVACY_SCAN_EXTS.has(ext)) continue;
    const buf = readFileSync(file.abs);
    for (const needle of forbidBuffers) {
      if (buf.includes(needle.buf)) {
        errors.push(`build-machine path leaked into ${file.rel}: ${needle.text}`);
      }
    }
    if (ext === '.wasm') continue;
    const text = buf.toString('utf8');
    for (const pattern of GENERIC_HOST_PATH_PATTERNS) {
      const match = text.match(pattern);
      if (match) {
        errors.push(`host path pattern leaked into ${file.rel}: ${match[0]}`);
      }
    }
  }

  return {
    ok: errors.length === 0,
    errors,
    warnings,
    filesCount: files.length,
    largestFiles: largestFiles.map((f) => ({ rel: f.rel, size: f.size })),
    engineFiles: { serial: serialFile, threads: threadsFile },
  };
}

function parseCliArgs(argv) {
  const args = { dist: 'web/app/dist', headers: 'web/app/public/_headers', forbidPaths: [] };
  for (let i = 0; i < argv.length; i += 1) {
    const flag = argv[i];
    if (flag === '--dist') {
      args.dist = argv[++i];
    } else if (flag === '--headers') {
      args.headers = argv[++i];
    } else if (flag === '--forbid-path') {
      args.forbidPaths.push(argv[++i]);
    } else if (flag === '--help' || flag === '-h') {
      console.log('Usage: node web/tools/check_web_dist.mjs --dist <dir> [--headers <path>] [--forbid-path <path>...]');
      process.exit(0);
    } else {
      throw new Error(`unknown argument: ${flag}`);
    }
  }
  return args;
}

if (process.argv[1] && import.meta.url === pathToFileURL(resolve(process.argv[1])).href) {
  const cli = parseCliArgs(process.argv.slice(2));
  const result = checkWebDist(resolve(cli.dist), {
    headersPath: cli.headers ? resolve(cli.headers) : undefined,
    forbidPaths: cli.forbidPaths,
  });

  console.log(`Checked ${result.filesCount} files in ${cli.dist}`);
  for (const file of result.largestFiles) {
    const pct = ((file.size / CLOUDFLARE_ASSET_MAX_BYTES) * 100).toFixed(1);
    console.log(`  ${file.rel}: ${formatMiB(file.size)} (${pct}% of 25 MiB limit)`);
  }
  for (const warning of result.warnings) {
    console.warn(`WARNING: ${warning}`);
  }
  if (!result.ok) {
    for (const error of result.errors) {
      console.error(`ERROR: ${error}`);
    }
    process.exit(1);
  }
  console.log(`OK: serial=${result.engineFiles.serial}, threads=${result.engineFiles.threads}`);
}
