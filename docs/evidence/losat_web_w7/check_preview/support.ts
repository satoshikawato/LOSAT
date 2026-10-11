// What the preview check's specs share: the site under test, the record they write, and a
// small HTTP client for the header checks. Nothing here sends anything but GET and (for
// /api/visit) POST requests to the site named by LOSAT_WEB_PREVIEW_URL.
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import type { TestInfo } from '@playwright/test';
import { parseHeaderRules, type HeaderRule } from '../../../../web/app/build/headers';

const APP = fileURLToPath(new URL('../../../../web/app/', import.meta.url));

/** The origin of the site under test (LOSAT_WEB_PREVIEW_URL, required). */
export function previewOrigin(): string {
  const value = process.env['LOSAT_WEB_PREVIEW_URL'];
  if (value === undefined || value === '') throw new Error('LOSAT_WEB_PREVIEW_URL is not set (README.md)');
  const url = new URL(value);
  if (url.protocol !== 'https:' && url.hostname !== 'localhost' && url.hostname !== '127.0.0.1') {
    throw new Error(`LOSAT_WEB_PREVIEW_URL=${value}: a preview is served over https`);
  }
  return url.origin;
}

/** The `_headers` file the site was deployed with: LOSAT_WEB_PREVIEW_HEADERS, else web/app/public/_headers. */
export const HEADERS_PATH = process.env['LOSAT_WEB_PREVIEW_HEADERS'] || `${APP}public/_headers`;

export function headerRules(): HeaderRule[] {
  return parseHeaderRules(readFileSync(HEADERS_PATH, 'utf8'));
}

/** Adds a JSON record to the test; the reporter (record-reporter.ts) writes it into the browser's file. */
export async function record(testInfo: TestInfo, data: unknown): Promise<void> {
  await testInfo.attach('record', { body: JSON.stringify(data, null, 2), contentType: 'application/json' });
}

export interface Answer {
  readonly status: number;
  /** Header names in lower case. */
  readonly headers: Readonly<Record<string, string>>;
  readonly body: Buffer;
}

/**
 * One request to the site, without following redirects. The body is decoded (gzip, br) by
 * fetch, as a browser would see it, so its size and SHA-256 are the file's.
 */
export async function request(
  origin: string,
  path: string,
  init: { method?: string; headers?: Record<string, string>; body?: string } = {},
): Promise<Answer> {
  const response = await fetch(origin + path, { redirect: 'manual', ...init });
  const headers: Record<string, string> = {};
  response.headers.forEach((value, name) => {
    headers[name.toLowerCase()] = value;
  });
  return { status: response.status, headers, body: Buffer.from(await response.arrayBuffer()) };
}

export const sha256 = (bytes: Uint8Array): string => createHash('sha256').update(bytes).digest('hex');

/** A selection of response headers worth recording, in a fixed order. */
const RECORDED = [
  'content-type',
  'content-length',
  'content-encoding',
  'etag',
  'last-modified',
  'cache-control',
  'location',
  'cross-origin-opener-policy',
  'cross-origin-embedder-policy',
  'cross-origin-resource-policy',
  'content-security-policy',
  'x-content-type-options',
  'referrer-policy',
  'x-robots-tag',
  'cf-cache-status',
  'age',
  'vary',
  'allow',
];

/** The status and every header of an answer, for the record (the full header list, not only RECORDED). */
export function summarize(answer: Answer): { status: number; headers: Record<string, string>; bytes: number } {
  const headers: Record<string, string> = {};
  for (const name of RECORDED) if (answer.headers[name] !== undefined) headers[name] = answer.headers[name];
  for (const [name, value] of Object.entries(answer.headers).sort(([a], [b]) => (a < b ? -1 : 1))) headers[name] ??= value;
  return { status: answer.status, headers, bytes: answer.body.length };
}

/** The manifest of a release: just the fields the check reads. */
export interface ManifestView {
  readonly release: string;
  readonly engine: { readonly serial: string; readonly threads: string } | null;
  readonly entry: { readonly script: string; readonly styles: readonly string[] };
  readonly assets: readonly { readonly path: string; readonly size: number; readonly sha256: string }[];
}

export interface IndexView {
  readonly latest: string;
  readonly releases: readonly { readonly id: string }[];
  readonly analytics?: { readonly endpoint: string };
  readonly shell?: readonly { readonly path: string; readonly size: number; readonly sha256: string }[];
}

/** The site index and the latest release's manifest, read as a browser reads them. */
export async function readSite(origin: string): Promise<{ index: IndexView; manifest: ManifestView }> {
  const indexAnswer = await request(origin, '/releases.json', { headers: { 'cache-control': 'no-cache' } });
  if (indexAnswer.status !== 200) throw new Error(`/releases.json answered ${indexAnswer.status}`);
  const index = JSON.parse(indexAnswer.body.toString('utf8')) as IndexView;
  const manifestPath = `/r/${index.latest}/release-manifest.json`;
  const manifestAnswer = await request(origin, manifestPath);
  if (manifestAnswer.status !== 200) throw new Error(`${manifestPath} answered ${manifestAnswer.status}`);
  return { index, manifest: JSON.parse(manifestAnswer.body.toString('utf8')) as ManifestView };
}
