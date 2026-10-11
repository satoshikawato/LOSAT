// Reads the Cloudflare `_headers` file so that the Vite servers, the local site server
// (build/site-server.ts) and the tests use the same response headers as the deployed site
// (one source of truth). The rules follow Cloudflare Pages (docs/web/release_manifest.md,
// "Response headers"): a rule is a URL pattern with indented `Name: value` lines; every rule
// whose pattern matches a request's path applies; a header that several matching rules
// attach is joined with ", "; `! Name` detaches a header that another rule attached.
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';

export const HEADERS_FILE = fileURLToPath(new URL('../public/_headers', import.meta.url));

/** Cloudflare Pages' limits on a `_headers` file. */
const MAX_RULES = 100;
const MAX_LINE = 2000;

export interface HeaderRule {
  /** A path that may hold one `*` (the splat) and `:name` placeholders. */
  readonly pattern: string;
  /** The headers the rule attaches, in file order. */
  readonly attach: readonly (readonly [name: string, value: string])[];
  /** The headers the rule detaches (`! Name`). */
  readonly detach: readonly string[];
}

/** Parses the rules of a Cloudflare `_headers` file, in file order. */
export function parseHeaderRules(text: string): HeaderRule[] {
  const rules: { pattern: string; attach: [string, string][]; detach: string[] }[] = [];
  for (const line of text.split(/\r?\n/)) {
    if (line.length > MAX_LINE) throw new Error(`_headers line longer than ${MAX_LINE} characters: ${line.slice(0, 60)}...`);
    if (line.trim() === '' || line.trimStart().startsWith('#')) continue;
    if (!/^\s/.test(line)) {
      const pattern = line.trim();
      // Cloudflare also takes https:// URLs (matched on the host); this site uses paths only.
      if (!pattern.startsWith('/')) throw new Error(`_headers: only path patterns are supported, not ${pattern}`);
      if ((pattern.match(/\*/g) ?? []).length > 1) throw new Error(`_headers: ${pattern} has more than one splat`);
      rules.push({ pattern, attach: [], detach: [] });
      if (rules.length > MAX_RULES) throw new Error(`_headers has more than ${MAX_RULES} rules`);
      continue;
    }
    const current = rules.at(-1);
    const trimmed = line.trim();
    if (current !== undefined && trimmed.startsWith('!')) {
      current.detach.push(trimmed.slice(1).trim());
      continue;
    }
    const separator = line.indexOf(':');
    if (current === undefined || separator < 0) {
      throw new Error(`malformed _headers line: ${line}`);
    }
    current.attach.push([line.slice(0, separator).trim(), line.slice(separator + 1).trim()]);
  }
  return rules;
}

/** Parses the rules of a Cloudflare `_headers` file into `path pattern -> attached headers`. */
export function parseHeadersFile(text: string): Map<string, Record<string, string>> {
  return new Map(parseHeaderRules(text).map((rule) => [rule.pattern, Object.fromEntries(rule.attach)]));
}

/**
 * Matches a rule's pattern against a request path: the values of its placeholders (`splat` for
 * the `*`), or undefined when it does not match. The splat matches any characters (none
 * included), a placeholder one or more characters other than `/`.
 */
export function matchPattern(pattern: string, path: string): Record<string, string> | undefined {
  let source = '';
  for (const [, splat, name, text] of pattern.matchAll(/(\*)|:([A-Za-z]\w*)|([^*:]+|:)/g)) {
    if (splat !== undefined) source += '(?<splat>.*)';
    else if (name !== undefined) source += `(?<${name}>[^/]+)`;
    else source += text!.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  }
  const match = new RegExp(`^${source}$`).exec(path);
  return match === null ? undefined : { ...match.groups };
}

/** The headers of every rule that matches `path`, with placeholders filled in and detached headers removed. */
export function headersFor(path: string, rules: readonly HeaderRule[] = readHeaderRules()): Record<string, string> {
  const headers = new Map<string, [name: string, value: string]>();
  const detached = new Set<string>();
  for (const rule of rules) {
    const values = matchPattern(rule.pattern, path);
    if (values === undefined) continue;
    for (const [name, raw] of rule.attach) {
      // Only the names of this rule's pattern are placeholders; any other `:word` stays as it is.
      const value = raw.replace(/:([A-Za-z]\w*)/g, (text, key: string) => values[key] ?? text);
      const key = name.toLowerCase();
      const previous = headers.get(key);
      headers.set(key, previous === undefined ? [name, value] : [previous[0], `${previous[1]}, ${value}`]);
    }
    for (const name of rule.detach) detached.add(name.toLowerCase());
  }
  return Object.fromEntries([...headers].filter(([key]) => !detached.has(key)).map(([, entry]) => entry));
}

/** The rules of public/_headers, or of another `_headers` file. */
export function readHeaderRules(file: string = HEADERS_FILE): HeaderRule[] {
  return parseHeaderRules(readFileSync(file, 'utf8'));
}

/** Headers that apply to every path (`/*`). */
export function siteHeaders(): Record<string, string> {
  const headers = parseHeadersFile(readFileSync(HEADERS_FILE, 'utf8')).get('/*');
  if (headers === undefined) throw new Error('_headers has no "/*" rule');
  return headers;
}
