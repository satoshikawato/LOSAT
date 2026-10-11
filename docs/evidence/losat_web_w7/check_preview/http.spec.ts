// Check 2 and 5 of the preview check (README.md): what the deployed site answers for each
// kind of response, with Node's fetch (no browser). Headers of the 200 (or 404) and of the
// conditional request that should give 304 are compared with the site's `_headers`
// (build/headers.ts `headersFor`); a 304 without COEP/CORP is a failure (WebKit refuses an
// Engine worker re-created after a cancel when the 304 lacks them: W3 gate record, decision 12).
import { expect, test } from '@playwright/test';
import { randomUUID } from 'node:crypto';
import { headersFor } from '../../../../web/app/build/headers';
import {
  headerRules,
  previewOrigin,
  readSite,
  record,
  request,
  sha256,
  summarize,
  type Answer,
  type IndexView,
  type ManifestView,
} from './support';

const ORIGIN = previewOrigin();
const RULES = headerRules();

/** The headers that isolate the page; a response without them breaks the app, a 304 included. */
const ISOLATION = [
  'cross-origin-opener-policy',
  'cross-origin-embedder-policy',
  'cross-origin-resource-policy',
  'content-security-policy',
];

interface Comparison {
  readonly failures: string[];
  readonly notes: string[];
}

/**
 * Compares an answer with the headers of the `_headers` rules for `path`. The isolation
 * headers must be equal in every answer. The other headers must be equal in a 200 and are
 * notes elsewhere (a 304 or a 404 that lacks `Cache-Control` or carries the rule's one for the
 * path is worth knowing, not a failure).
 */
function compare(path: string, answer: Answer, strict: 'all' | 'isolation'): Comparison {
  const expected = headersFor(path, RULES);
  const failures: string[] = [];
  const notes: string[] = [];
  for (const [name, value] of Object.entries(expected)) {
    const actual = answer.headers[name.toLowerCase()];
    if (actual === value) continue;
    const message = `${name}: expected ${JSON.stringify(value)}, got ${actual === undefined ? 'nothing' : JSON.stringify(actual)}`;
    (ISOLATION.includes(name.toLowerCase()) || strict === 'all' ? failures : notes).push(message);
  }
  return { failures, notes };
}

const CONTENT_TYPES: Readonly<Record<string, RegExp>> = {
  html: /^text\/html\b/i,
  js: /^(text|application)\/javascript\b/i,
  css: /^text\/css\b/i,
  json: /^application\/json\b/i,
  wasm: /^application\/wasm$/i,
};

const extensionOf = (path: string): string => path.slice(path.lastIndexOf('.') + 1);

/** One kind of response: its path in this site, or the reason it does not apply. */
interface Kind {
  readonly id: string;
  readonly path: (site: { index: IndexView; manifest: ManifestView; swPresent: boolean }) => string | { skip: string };
  /** The content type that a browser needs (wasm needs application/wasm to compile). */
  readonly type?: RegExp;
  readonly missing?: boolean;
}

const release = (site: { index: IndexView }) => `/r/${site.index.latest}/`;
const pick = (site: { manifest: ManifestView; index: IndexView }, match: (path: string) => boolean, what: string) => {
  const asset = site.manifest.assets.find((file) => match(file.path));
  return asset === undefined ? { skip: `the manifest has no ${what}` } : `${release(site)}${asset.path}`;
};

const KINDS: readonly Kind[] = [
  { id: 'root', path: () => '/', type: CONTENT_TYPES['html']! },
  { id: 'boot-js', path: () => '/boot.js', type: CONTENT_TYPES['js']! },
  { id: 'boot-css', path: () => '/boot.css', type: CONTENT_TYPES['css']! },
  { id: 'releases-json', path: () => '/releases.json', type: CONTENT_TYPES['json']! },
  { id: 'sw-js', path: (site) => (site.swPresent ? '/sw.js' : { skip: '/sw.js is not there (a site without the service worker)' }), type: CONTENT_TYPES['js']! },
  { id: 'release-manifest', path: (site) => `${release(site)}release-manifest.json`, type: CONTENT_TYPES['json']! },
  { id: 'entry-script', path: (site) => `${release(site)}${site.manifest.entry.script}`, type: CONTENT_TYPES['js']! },
  { id: 'entry-style', path: (site) => `${release(site)}${site.manifest.entry.styles[0] ?? ''}`, type: CONTENT_TYPES['css']! },
  { id: 'worker-script', path: (site) => pick(site, (path) => /engine-worker.*\.js$/.test(path), 'engine worker script'), type: CONTENT_TYPES['js']! },
  { id: 'engine-module-1', path: (site) => (site.manifest.engine === null ? { skip: 'FakeEngine site (engine is null)' } : pick(site, (path) => path.endsWith('.wasm'), 'engine module')), type: CONTENT_TYPES['wasm']! },
  {
    id: 'engine-module-2',
    path: (site) => {
      if (site.manifest.engine === null) return { skip: 'FakeEngine site (engine is null)' };
      const second = site.manifest.assets.filter((file) => file.path.endsWith('.wasm'))[1];
      return second === undefined ? { skip: 'the manifest has one engine module only' } : `${release(site)}${second.path}`;
    },
    type: CONTENT_TYPES['wasm']!,
  },
  { id: 'missing-site-file', path: () => `/missing-${randomUUID()}.js`, missing: true },
  { id: 'missing-release-file', path: (site) => `${release(site)}assets/missing-${randomUUID()}.js`, missing: true },
];

let site: { index: IndexView; manifest: ManifestView; swPresent: boolean };

test.beforeAll(async () => {
  const { index, manifest } = await readSite(ORIGIN);
  const sw = await request(ORIGIN, '/sw.js');
  site = { index, manifest, swPresent: sw.status === 200 && /javascript/i.test(sw.headers['content-type'] ?? '') };
});

for (const kind of KINDS) {
  test(`headers: ${kind.id}`, async ({}, testInfo) => {
    const path = kind.path(site);
    test.skip(typeof path !== 'string', typeof path === 'string' ? '' : path.skip);
    if (typeof path !== 'string') return;
    const failures: string[] = [];
    const notes: string[] = [];
    const first = await request(ORIGIN, path);
    const data: Record<string, unknown> = { kind: kind.id, path, first: summarize(first) };

    if (kind.missing === true) {
      // A missing file is a 404 (not the boot page with 200: Pages' single-page mode without a
      // top-level 404.html), and the isolation headers are on it.
      if (first.status !== 404) failures.push(`status ${first.status}, expected 404`);
      const result = compare(path, first, 'isolation');
      failures.push(...result.failures);
      notes.push(...result.notes);
      const cache = first.headers['cache-control'] ?? '(none)';
      if (/immutable/i.test(cache)) {
        notes.push(`the 404 carries Cache-Control: ${cache}; a browser or CDN may remember that this file is missing (release IDs are never reused, so it matters only for a release that is removed)`);
      }
      data['body'] = first.body.subarray(0, 300).toString('utf8');
    } else {
      if (first.status !== 200) failures.push(`status ${first.status}, expected 200`);
      const result = compare(path, first, 'all');
      failures.push(...result.failures);
      notes.push(...result.notes);
      if (kind.type !== undefined && !kind.type.test(first.headers['content-type'] ?? '')) {
        failures.push(`content-type ${JSON.stringify(first.headers['content-type'])} is not ${String(kind.type)}`);
      }
      // The validators, and the conditional requests that should give 304.
      const etag = first.headers['etag'];
      const modified = first.headers['last-modified'];
      if (etag === undefined) failures.push('the 200 has no ETag (Pages always sends one)');
      const conditional: Record<string, unknown> = {};
      const attempts: [string, Record<string, string>][] = [];
      if (etag !== undefined) {
        attempts.push(['if-none-match', { 'if-none-match': etag, ...(modified === undefined ? {} : { 'if-modified-since': modified }) }]);
      }
      if (modified !== undefined) attempts.push(['if-modified-since', { 'if-modified-since': modified }]);
      for (const [label, headers] of attempts) {
        const again = await request(ORIGIN, path, { headers });
        conditional[label] = summarize(again);
        if (again.status !== 304) {
          failures.push(`${label}: status ${again.status}, expected 304`);
          continue;
        }
        if (again.body.length !== 0) failures.push(`${label}: the 304 has a body`);
        const result = compare(path, again, 'isolation');
        failures.push(...result.failures.map((message) => `${label}: the 304 -> ${message}`));
        notes.push(...result.notes.map((message) => `${label}: the 304 -> ${message}`));
      }
      data['conditional'] = conditional;
    }
    data['failures'] = failures;
    data['notes'] = notes;
    await record(testInfo, data);
    expect(failures, `${kind.id} (${path})`).toEqual([]);
  });
}

test('the files of the latest release are served byte for byte', async ({}, testInfo) => {
  // What the service worker will check before it says "Offline ready"; a host that rewrites
  // a file (minifying, injecting a script) is found here, whatever the headers say.
  const failures: string[] = [];
  const base = `/r/${site.index.latest}/`;
  for (const file of site.manifest.assets) {
    const answer = await request(ORIGIN, base + file.path);
    if (answer.status !== 200) failures.push(`${file.path}: status ${answer.status}`);
    else if (answer.body.length !== file.size || sha256(answer.body) !== file.sha256) {
      failures.push(`${file.path}: ${answer.body.length} bytes, sha256 ${sha256(answer.body).slice(0, 12)}...; the manifest lists ${file.size} bytes, ${file.sha256.slice(0, 12)}...`);
    }
  }
  const shell = site.index.shell ?? [];
  for (const file of shell) {
    const answer = await request(ORIGIN, file.path);
    if (answer.status !== 200 || answer.body.length !== file.size || sha256(answer.body) !== file.sha256) {
      failures.push(`${file.path}: not the shell file that /releases.json lists (status ${answer.status}, ${answer.body.length} bytes)`);
    }
  }
  await record(testInfo, { release: site.index.latest, assets: site.manifest.assets.length, shell: shell.length, failures });
  expect(failures).toEqual([]);
});

test('the boot page loads nothing from another origin and nothing was injected', async ({}, testInfo) => {
  // Cloudflare can inject its Web Analytics beacon into HTML (a project setting); the page
  // handles research data, so every script must be the site's own (web/app/AGENTS.md).
  const page = (await request(ORIGIN, '/')).body.toString('utf8');
  const sources = [...page.matchAll(/<(?:script|link|img|iframe)\b[^>]*?\b(?:src|href)\s*=\s*["']?([^"'\s>]+)/gi)].map((match) => match[1]!);
  const foreign = sources.filter((source) => /^(?:[a-z][a-z0-9+.-]*:)?\/\//i.test(source) && new URL(source, ORIGIN).origin !== ORIGIN);
  const injected = /cloudflareinsights|beacon\.min\.js|rocket-loader|cf-beacon/i.test(page);
  const inline = /<script(?![^>]*\bsrc=)[^>]*>\s*\S/i.test(page);
  await record(testInfo, { sources, foreign, injected, inline });
  expect(foreign).toEqual([]);
  expect(injected).toBe(false);
  expect(inline, 'an inline script in /').toBe(false);
});

test('every release of the site index has its manifest', async ({}, testInfo) => {
  const failures: string[] = [];
  const statuses: Record<string, number> = {};
  for (const { id } of site.index.releases) {
    const answer = await request(ORIGIN, `/r/${id}/release-manifest.json`);
    statuses[id] = answer.status;
    if (answer.status !== 200) {
      failures.push(`${id}: status ${answer.status}`);
      continue;
    }
    const manifest = JSON.parse(answer.body.toString('utf8')) as { release?: unknown };
    if (manifest.release !== id) failures.push(`${id}: its manifest names release ${String(manifest.release)}`);
  }
  await record(testInfo, { latest: site.index.latest, statuses, failures });
  expect(failures).toEqual([]);
});

test('/api/visit refuses a body that is not the event; a valid event only with LOSAT_WEB_PREVIEW_SEND_VISIT=1', async ({}, testInfo) => {
  test.skip(site.index.analytics === undefined, '/releases.json has no analytics: the site measures nothing');
  const endpoint = site.index.analytics!.endpoint;
  const headers = { origin: ORIGIN, 'content-type': 'application/json' };
  const refused = await request(ORIGIN, endpoint, { method: 'POST', headers, body: JSON.stringify({ v: 1, release: site.index.latest, client: 'not-a-client-id' }) });
  const wrongMethod = await request(ORIGIN, endpoint);
  const data: Record<string, unknown> = { endpoint, refused: summarize(refused), wrongMethod: summarize(wrongMethod) };
  // Pages does not apply `_headers` to a Function's answer (the Function sets its own).
  data['note'] = 'the answers of a Pages Function carry the headers the Function sets, not the _headers rules';
  const failures: string[] = [];
  if (refused.status !== 400) failures.push(`a refused body: status ${refused.status}, expected 400`);
  if (wrongMethod.status !== 405) failures.push(`GET: status ${wrongMethod.status}, expected 405`);
  if (process.env['LOSAT_WEB_PREVIEW_SEND_VISIT'] === '1') {
    // A real event, release "preview-check". If the project has the GA4 secrets set in this
    // environment, it is forwarded to GA4: the preview should not have them.
    const body = JSON.stringify({ v: 1, release: 'preview-check', client: randomUUID() });
    const sent = await request(ORIGIN, endpoint, { method: 'POST', headers, body });
    data['sent'] = summarize(sent);
    if (sent.status !== 204) failures.push(`a valid event: status ${sent.status}, expected 204`);
  } else {
    data['sent'] = 'not sent (LOSAT_WEB_PREVIEW_SEND_VISIT is not 1)';
  }
  data['failures'] = failures;
  await record(testInfo, data);
  expect(failures).toEqual([]);
});
