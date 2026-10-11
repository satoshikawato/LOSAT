// What V-PRIV looks for (plan §5.10, TD-11; design §13.3), as functions of text and files only, so
// that tests/unit/privacy-scan.test.ts can show that they find what they are meant to find: the
// canaries of a test run in what the browser sent, the rules for the requests themselves, and the
// paths of the build machine in the files that the site serves.
import { randomBytes } from 'node:crypto';
import { readdirSync, readFileSync } from 'node:fs';
import { join, relative, sep } from 'node:path';

// --- canaries --------------------------------------------------------------------------------------

/** Unique strings of one test run, put into every kind of research data that the flows handle. */
export interface Canary {
  /** `LOSATPRIV` and random letters: in file names, deflines, the job title, notes. */
  readonly token: string;
  /** A random nucleotide 40-mer in the nucleotide sequences. */
  readonly nucleotides: string;
  /** A random 20-residue motif in the protein sequences. */
  readonly protein: string;
}

const pick = (alphabet: string, length: number, random: (size: number) => Uint8Array): string =>
  [...random(length)].map((byte) => alphabet[byte % alphabet.length]).join('');

export function makeCanary(random: (size: number) => Uint8Array = randomBytes): Canary {
  return {
    token: `LOSATPRIV${pick('ABCDEFGHJKLMNPQRSTUVWXYZ23456789', 12, random)}`,
    nucleotides: pick('ACGT', 40, random),
    protein: pick('ACDEFGHIKLMNPQRSTVWY', 20, random),
  };
}

export interface CanaryForm {
  /** Which canary, and in which form. */
  readonly label: string;
  /** Lower case: the search ignores case. */
  readonly text: string;
}

const COMPLEMENT: Readonly<Record<string, string>> = { A: 'T', C: 'G', G: 'C', T: 'A' };
const reverseComplement = (letters: string) =>
  [...letters]
    .reverse()
    .map((letter) => COMPLEMENT[letter] ?? letter)
    .join('');
const percentEncoded = (text: string) =>
  [...new TextEncoder().encode(text)].map((byte) => `%${byte.toString(16).padStart(2, '0')}`).join('');

/**
 * The forms in which a canary could travel: as written (in either case; percent-encoding leaves
 * these letters as they are) and percent-encoded byte by byte; the 40-mer also as its reverse
 * complement.
 */
export function canaryForms(canary: Canary): readonly CanaryForm[] {
  const plain: Array<[string, string]> = [
    ['token', canary.token],
    ['nucleotide 40-mer', canary.nucleotides],
    ['nucleotide 40-mer, reverse complement', reverseComplement(canary.nucleotides)],
    ['protein motif', canary.protein],
  ];
  return plain.flatMap(([label, text]) => [
    { label, text: text.toLowerCase() },
    { label: `${label}, percent-encoded`, text: percentEncoded(text).toLowerCase() },
  ]);
}

/**
 * The labels of the canary forms in a text (a URL, header lines, or a body read as Latin-1, so that
 * binary bodies are searched too), ignoring case; also with white space removed, for residues that
 * a FASTA line break split.
 */
export function findCanaries(text: string, forms: readonly CanaryForm[]): string[] {
  const lower = text.toLowerCase();
  const joined = lower.replace(/\s+/g, '');
  return forms.filter((form) => lower.includes(form.text) || joined.includes(form.text)).map((form) => form.label);
}

// --- requests ----------------------------------------------------------------------------------------

/** One request, as the site's server received it or as the browser reported it to Playwright. */
export interface Exchange {
  readonly seenBy: 'server' | 'browser';
  readonly method: string;
  /** Absolute. */
  readonly url: string;
  readonly headers: Readonly<Record<string, string>>;
  readonly body: Uint8Array;
  /** The browser's resource type, `service worker` for a service worker's request; `server` for the server's. */
  readonly kind: string;
}

/** What a page of the site reported about itself (support/privacy.ts's init script). */
export interface PageNote {
  readonly kind: 'title' | 'href' | 'csp' | 'websocket';
  readonly value: string;
}

export interface RequestRules {
  /** The site's origin, for example http://127.0.0.1:41234. */
  readonly origin: string;
  /** The paths that the site serves (`/`, `/assets/index-x.js`, ...). */
  readonly files: ReadonlySet<string>;
  /** The path to which the server asks browsers to send CSP violation reports. */
  readonly cspReportPath: string;
  /** Paths the site does not serve that a browser may ask for by itself (none so far). */
  readonly extraPaths?: ReadonlySet<string>;
  /** A request with a method other than GET or HEAD, or with a body, that the flow allows (measurement). */
  readonly allowBody?: (exchange: Exchange) => boolean;
}

const latin1 = (bytes: Uint8Array) => Buffer.from(bytes.buffer, bytes.byteOffset, bytes.byteLength).toString('latin1');
const decoded = (url: string) => {
  try {
    return decodeURIComponent(url);
  } catch {
    return url;
  }
};

/** A CSP report's directive, blocked URL and document (report-uri's JSON), or the start of its body. */
function cspReport(body: string): string {
  try {
    const report = (JSON.parse(body) as { 'csp-report'?: Record<string, unknown> })['csp-report'];
    if (report !== undefined) return `${String(report['effective-directive'])} blocked ${String(report['blocked-uri'])} in ${String(report['document-uri'])}`;
  } catch {
    // Not JSON: shown as it is.
  }
  return body.slice(0, 300);
}

/**
 * Every way in which the requests and the pages' notes break V-PRIV, one line each: a request to
 * another origin; a method other than GET and HEAD or a body (unless the flow allows it); a path
 * that is not a file of the site; a canary in a URL, a header or a body; a CSP violation (an
 * attempt that the CSP blocked is a finding, not a pass); a canary in a page's title or address;
 * a WebSocket.
 */
export function findings(
  exchanges: readonly Exchange[],
  notes: readonly PageNote[],
  forms: readonly CanaryForm[],
  rules: RequestRules,
): string[] {
  const found = new Set<string>();
  for (const exchange of exchanges) {
    const where = `${exchange.seenBy} ${exchange.method} ${exchange.url.slice(0, 160)}`;
    for (const label of findCanaries(`${exchange.url}\n${decoded(exchange.url)}`, forms)) found.add(`canary (${label}) in the URL: ${where}`);
    const headers = Object.entries(exchange.headers)
      .map(([name, value]) => `${name}: ${value}`)
      .join('\n');
    for (const label of findCanaries(headers, forms)) found.add(`canary (${label}) in a header: ${where}`);
    for (const label of findCanaries(latin1(exchange.body), forms)) found.add(`canary (${label}) in the body: ${where}`);

    // data: URLs are read in the page; they reach no server.
    if (exchange.url.startsWith('data:')) continue;
    const url = new URL(exchange.url);
    if (url.origin !== rules.origin) {
      found.add(`a request to another origin: ${where}`);
      continue;
    }
    // A blob: URL of the site's origin is a file made in the page (a download), not a request.
    if (url.protocol === 'blob:') continue;
    if (url.pathname === rules.cspReportPath) {
      found.add(`a CSP violation report: ${cspReport(latin1(exchange.body))}`);
      continue;
    }
    if (!['GET', 'HEAD'].includes(exchange.method) || exchange.body.length > 0) {
      if (!(rules.allowBody?.(exchange) ?? false)) found.add(`a request that sends a body: ${where} (${exchange.body.length} bytes)`);
    } else if (!rules.files.has(url.pathname) && !(rules.extraPaths?.has(url.pathname) ?? false)) {
      found.add(`a request for a path that is not a file of the site: ${where}`);
    }
  }
  for (const note of notes) {
    if (note.kind === 'csp') found.add(`a CSP violation in a page: ${note.value}`);
    else if (note.kind === 'websocket') found.add(`a WebSocket: ${note.value}`);
    else for (const label of findCanaries(`${note.value}\n${decoded(note.value)}`, forms)) found.add(`canary (${label}) in the page's ${note.kind}: ${note.value}`);
  }
  return [...found];
}

// --- machine paths (TD-11) ---------------------------------------------------------------------------

/** Where the build machine keeps things: what must not be in a served file. */
export interface MachineContext {
  /** os.userInfo().username */
  readonly user: string;
  readonly home?: string | undefined;
  readonly cargoHome?: string | undefined;
  readonly rustupHome?: string | undefined;
  /** The checkout's root, its real path, the main checkout of a worktree. */
  readonly checkouts: readonly string[];
}

export interface MachinePattern {
  readonly label: string;
  readonly pattern: RegExp;
}

/**
 * User names that are also words (CI's user is `runner`): as a whole word they may be part of a
 * library's text, so for them, and for names shorter than four letters, only a path-like context
 * counts (a slash, a backslash or `~` next to the name).
 */
export const COMMON_USER_NAMES: ReadonlySet<string> = new Set([
  'admin',
  'administrator',
  'app',
  'build',
  'builder',
  'ci',
  'circleci',
  'codespace',
  'dev',
  'docker',
  'jenkins',
  'node',
  'root',
  'runner',
  'runneradmin',
  'test',
  'ubuntu',
  'user',
  'vagrant',
  'vscode',
]);

const escape = (text: string) => text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
/** A name of a path does not continue after the match. */
const END = '(?![A-Za-z0-9_-])';

/** Whether the user name is checked only where it looks like part of a path. */
export const userNameInPathsOnly = (user: string) => user.length < 4 || COMMON_USER_NAMES.has(user.toLowerCase());

/**
 * The patterns of TD-11: each machine path as written, with JSON's escaped slashes (`\/`) and with
 * doubled backslashes; the home folders of the user on Linux, macOS and Windows; and the user name
 * as a whole word (in path-like contexts only, for a common or short name).
 */
export function machinePatterns(context: MachineContext): MachinePattern[] {
  const patterns = new Map<string, MachinePattern>();
  const add = (label: string, source: string) => {
    if (!patterns.has(source)) patterns.set(source, { label, pattern: new RegExp(source, 'g') });
  };
  const path = (label: string, value: string | undefined) => {
    if (value === undefined || value.replace(/[/\\]+$/, '').length < 4) return;
    const trimmed = value.replace(/[/\\]+$/, '');
    for (const form of new Set([trimmed, trimmed.replaceAll('/', '\\/'), trimmed.replaceAll('\\', '\\\\')])) add(label, escape(form) + END);
  };
  for (const checkout of context.checkouts) path('the checkout', checkout);
  path('CARGO_HOME', context.cargoHome);
  path('RUSTUP_HOME', context.rustupHome);
  path('HOME', context.home);
  const user = context.user;
  if (user === '') return [...patterns.values()];
  path('/home/<user>', `/home/${user}`);
  path('/Users/<user>', `/Users/${user}`);
  path('\\Users\\<user>', `\\Users\\${user}`);
  const name = escape(user);
  if (userNameInPathsOnly(user)) add('the user name in a path', `(?<=[/\\\\~])${name}${END}|(?<![A-Za-z0-9_-])${name}(?=[/\\\\])`);
  else add('the user name', `(?<![A-Za-z0-9_])${name}(?![A-Za-z0-9_])`);
  return [...patterns.values()];
}

export interface MachineFinding {
  readonly file: string;
  readonly label: string;
  readonly excerpt: string;
}

/** The machine paths in one file, read as Latin-1 so that Wasm data and JavaScript are searched alike. */
export function findMachinePaths(file: string, bytes: Uint8Array, patterns: readonly MachinePattern[]): MachineFinding[] {
  const text = latin1(bytes);
  const out: MachineFinding[] = [];
  for (const { label, pattern } of patterns) {
    for (const match of text.matchAll(pattern)) {
      const at = match.index ?? 0;
      const excerpt = text.slice(Math.max(0, at - 30), at + match[0].length + 30).replace(/[^\x20-\x7e]/g, '.');
      out.push({ file, label, excerpt });
      if (out.length > 50) return out;
    }
  }
  return out;
}

/** Every file under a directory, as paths relative to it with forward slashes. */
export function filesUnder(root: string): string[] {
  const out: string[] = [];
  const walk = (directory: string) => {
    for (const entry of readdirSync(directory, { withFileTypes: true })) {
      const path = join(directory, entry.name);
      if (entry.isDirectory()) walk(path);
      else if (entry.isFile()) out.push(relative(root, path).split(sep).join('/'));
    }
  };
  walk(root);
  return out.sort();
}

/** The machine paths in the named files of a directory. */
export function scanFiles(root: string, files: readonly string[], patterns: readonly MachinePattern[]): MachineFinding[] {
  return files.flatMap((file) => findMachinePaths(file, readFileSync(join(root, file)), patterns));
}
