// Reads the Cloudflare `_headers` file so that the Vite servers and the tests use the
// same response headers as the deployed site (one source of truth).
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';

export const HEADERS_FILE = fileURLToPath(new URL('../public/_headers', import.meta.url));

/** Parses the rules of a Cloudflare `_headers` file into `path pattern -> headers`. */
export function parseHeadersFile(text: string): Map<string, Record<string, string>> {
  const rules = new Map<string, Record<string, string>>();
  let current: Record<string, string> | undefined;
  for (const line of text.split(/\r?\n/)) {
    if (line.trim() === '' || line.trimStart().startsWith('#')) continue;
    if (!/^\s/.test(line)) {
      current = {};
      rules.set(line.trim(), current);
      continue;
    }
    const separator = line.indexOf(':');
    if (current === undefined || separator < 0) {
      throw new Error(`malformed _headers line: ${line}`);
    }
    current[line.slice(0, separator).trim()] = line.slice(separator + 1).trim();
  }
  return rules;
}

/** Headers that apply to every path (`/*`). */
export function siteHeaders(): Record<string, string> {
  const headers = parseHeadersFile(readFileSync(HEADERS_FILE, 'utf8')).get('/*');
  if (headers === undefined) throw new Error('_headers has no "/*" rule');
  return headers;
}
