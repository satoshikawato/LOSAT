// V-PRIV's scanners (tests/e2e/support/privacy-scan.ts) find what they are meant to find: the
// canaries in the forms they could travel in, the requests that break the rules, and planted paths
// of a build machine in served files (TD-11).
import { mkdirSync, mkdtempSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { afterEach, describe, expect, it } from 'vitest';
import {
  canaryForms,
  filesUnder,
  findCanaries,
  findings,
  machinePatterns,
  makeCanary,
  scanFiles,
  type Exchange,
  type RequestRules,
} from '../e2e/support/privacy-scan';

const CANARY = { token: 'LOSATPRIVABCDEFGH2345', nucleotides: 'ACGTTGCAAGCTTCGAGGATCCTTAAGGCCATGGTACCGA', protein: 'MWKHYQRDEFNPCTGSAVLI' };
const FORMS = canaryForms(CANARY);
const ORIGIN = 'http://127.0.0.1:4000';
const RULES: RequestRules = { origin: ORIGIN, files: new Set(['/', '/assets/index-a.js']), cspReportPath: '/__privacy/csp-report' };
const get = (url: string, extra: Partial<Exchange> = {}): Exchange => ({
  seenBy: 'server',
  method: 'GET',
  url,
  headers: {},
  body: new Uint8Array(),
  kind: 'server',
  ...extra,
});

describe('canaries', () => {
  it('are unique strings of the expected kinds', () => {
    const [a, b] = [makeCanary(), makeCanary()];
    expect(a.token).toMatch(/^LOSATPRIV[A-Z2-9]{12}$/);
    expect(a.nucleotides).toMatch(/^[ACGT]{40}$/);
    expect(a.protein).toMatch(/^[ACDEFGHIKLMNPQRSTVWY]{20}$/);
    expect(b.token).not.toBe(a.token);
    expect(b.nucleotides).not.toBe(a.nucleotides);
  });

  it('are found as written, in another case, percent-encoded, split by a line break and reverse-complemented', () => {
    expect(findCanaries(`/x?q=${CANARY.token}`, FORMS)).toEqual(['token']);
    expect(findCanaries(`title: ${CANARY.token.toLowerCase()}`, FORMS)).toEqual(['token']);
    const encoded = [...CANARY.token].map((letter) => `%${letter.charCodeAt(0).toString(16).toUpperCase()}`).join('');
    expect(findCanaries(`/x?q=${encoded}`, FORMS)).toEqual(['token', 'token, percent-encoded']);
    expect(findCanaries(`/x?q=LOSAT${encoded.slice(15)}`, FORMS)).toEqual(['token']);
    expect(findCanaries(`>s1\n${CANARY.nucleotides.slice(0, 25)}\n${CANARY.nucleotides.slice(25)}\n`, FORMS)).toEqual(['nucleotide 40-mer']);
    const reverse = [...CANARY.nucleotides].reverse().map((letter) => ({ A: 'T', C: 'G', G: 'C', T: 'A' })[letter]).join('');
    expect(findCanaries(reverse.toLowerCase(), FORMS)).toEqual(['nucleotide 40-mer, reverse complement']);
    expect(findCanaries(`MK${CANARY.protein}LL`, FORMS)).toEqual(['protein motif']);
    expect(findCanaries('GET /assets/index-a.js LOSATPRIV ACGTACGT', FORMS)).toEqual([]);
  });

  it('are found in base64 at any offset, in either alphabet and case, also percent-encoded, line-wrapped or split by an escaped line break', () => {
    for (const before of ['', '>', '>q', '>q1']) {
      for (const text of [CANARY.token, CANARY.token.toLowerCase()]) {
        const bytes = Buffer.from(`${before}${text} and more`);
        expect(findCanaries(`d=${bytes.toString('base64')}`, FORMS), `${before}${text}`).toEqual(['token, base64']);
        expect(findCanaries(`d=${bytes.toString('base64url')}`, FORMS), `${before}${text}, base64url`).toEqual(['token, base64']);
      }
    }
    // Every letter of the base64 percent-encoded: found once the escapes are decoded.
    const base64 = Buffer.from(`>p1 ${CANARY.protein}`).toString('base64');
    const escaped = [...base64].map((letter) => `%${letter.charCodeAt(0).toString(16).toUpperCase()}`).join('');
    expect(findCanaries(`/x?d=${escaped}`, FORMS)).toEqual(['protein motif, base64']);
    const wrapped = Buffer.from(`>s1\n${CANARY.nucleotides}\n`).toString('base64').replace(/.{10}/g, '$&\r\n');
    expect(findCanaries(wrapped, FORMS)).toEqual(['nucleotide 40-mer, base64']);
    expect(findCanaries(JSON.stringify({ sequence: `${CANARY.nucleotides.slice(0, 20)}\n${CANARY.nucleotides.slice(20)}` }), FORMS)).toEqual(['nucleotide 40-mer']);
    expect(findCanaries(Buffer.from('LOSATPRIV and a sequence ACGTACGTACGT').toString('base64'), FORMS)).toEqual([]);
  });
});

describe('request findings', () => {
  it('pass the requests of the site for its files, blob: URLs of the site, and data: URLs', () => {
    expect(findings([get(`${ORIGIN}/`), get(`${ORIGIN}/assets/index-a.js`), get(`blob:${ORIGIN}/1f2e`), get('data:text/plain,x')], [], FORMS, RULES)).toEqual([]);
  });

  it('report another origin, a body, a path that is not a file of the site, and a canary anywhere', () => {
    const found = findings(
      [
        get('https://example.org/collect'),
        get(`${ORIGIN}/api/x`, { method: 'POST', body: new TextEncoder().encode('{}') }),
        get(`${ORIGIN}/nothing-here`),
        get(`${ORIGIN}/assets/index-a.js`, { seenBy: 'browser', headers: { 'x-note': CANARY.token } }),
        get(`${ORIGIN}/`, { method: 'POST', body: new TextEncoder().encode(`>q\n${CANARY.nucleotides}\n`) }),
        get(`${ORIGIN}/__privacy/csp-report`, {
          method: 'POST',
          body: new TextEncoder().encode(JSON.stringify({ 'csp-report': { 'effective-directive': 'connect-src', 'blocked-uri': 'https://example.org/x', 'document-uri': `${ORIGIN}/` } })),
        }),
      ],
      [],
      FORMS,
      RULES,
    );
    expect(found).toEqual([
      'a request to another origin: server GET https://example.org/collect',
      `a request that sends a body: server POST ${ORIGIN}/api/x (2 bytes)`,
      `a request for a path that is not a file of the site: server GET ${ORIGIN}/nothing-here`,
      `canary (token) in a header: browser GET ${ORIGIN}/assets/index-a.js`,
      `canary (nucleotide 40-mer) in the body: server POST ${ORIGIN}/`,
      `a request that sends a body: server POST ${ORIGIN}/ (44 bytes)`,
      `a CSP violation report: connect-src blocked https://example.org/x in ${ORIGIN}/`,
    ]);
  });

  it("report a query string on any request to the site, and a canary in base64 in a query or a header (code review 2, L2)", () => {
    const encoded = Buffer.from(`>${CANARY.token}\n${CANARY.nucleotides}`).toString('base64url');
    const rules = { ...RULES, extraPaths: new Set(['/favicon.ico']) };
    const found = (exchange: Exchange) => findings([exchange], [], FORMS, rules);
    const inQuery = (url: string) => [
      `canary (token, base64) in the URL: server GET ${url}`,
      `canary (nucleotide 40-mer, base64) in the URL: server GET ${url}`,
      `a query string on a request to the site: server GET ${url}`,
    ];
    for (const path of ['/assets/index-a.js', '/favicon.ico', '/']) expect(found(get(`${ORIGIN}${path}?d=${encoded}`)), path).toEqual(inQuery(`${ORIGIN}${path}?d=${encoded}`));
    expect(found(get(`${ORIGIN}/assets/index-a.js?v=2`))).toEqual([`a query string on a request to the site: server GET ${ORIGIN}/assets/index-a.js?v=2`]);
    expect(found(get(`${ORIGIN}/assets/index-a.js`, { headers: { 'x-d': encoded } }))).toEqual([
      `canary (token, base64) in a header: server GET ${ORIGIN}/assets/index-a.js`,
      `canary (nucleotide 40-mer, base64) in a header: server GET ${ORIGIN}/assets/index-a.js`,
    ]);
    const split = JSON.stringify({ sequence: `${CANARY.nucleotides.slice(0, 20)}\n${CANARY.nucleotides.slice(20)}` });
    expect(found(get(`${ORIGIN}/assets/index-a.js`, { headers: { 'x-d': split } }))).toEqual([`canary (nucleotide 40-mer) in a header: server GET ${ORIGIN}/assets/index-a.js`]);
    // Another origin's query is part of that request's finding.
    expect(found(get(`https://example.org/x?d=1`))).toEqual(['a request to another origin: server GET https://example.org/x?d=1']);
  });

  it('allow the paths and bodies that a flow names', () => {
    const rules = { ...RULES, extraPaths: new Set(['/favicon.ico']), allowBody: (exchange: Exchange) => exchange.url.endsWith('/api/visit') };
    expect(findings([get(`${ORIGIN}/favicon.ico`), get(`${ORIGIN}/api/visit`, { method: 'POST', body: new Uint8Array([123, 125]) })], [], FORMS, rules)).toEqual([]);
  });

  it('report a canary in a page title or address, an address with a query or a fragment, a CSP violation in a page, and a WebSocket', () => {
    expect(
      findings(
        [],
        [
          { kind: 'title', value: 'LOSAT Web' },
          { kind: 'title', value: `${CANARY.token}.fa - LOSAT Web` },
          { kind: 'href', value: `${ORIGIN}/#q=${CANARY.token}` },
          { kind: 'href', value: `${ORIGIN}/` },
          { kind: 'href', value: `${ORIGIN}/?s=1` },
          { kind: 'csp', value: 'connect-src blocked https://example.org' },
          { kind: 'websocket', value: 'ws://127.0.0.1:4000/' },
        ],
        FORMS,
        RULES,
      ),
    ).toEqual([
      `canary (token) in the page's title: ${CANARY.token}.fa - LOSAT Web`,
      `canary (token) in the page's href: ${ORIGIN}/#q=${CANARY.token}`,
      `an address with a query or a fragment: ${ORIGIN}/#q=${CANARY.token}`,
      `an address with a query or a fragment: ${ORIGIN}/?s=1`,
      'a CSP violation in a page: connect-src blocked https://example.org',
      'a WebSocket: ws://127.0.0.1:4000/',
    ]);
  });
});

describe('machine paths (TD-11)', () => {
  let root: string | undefined;
  afterEach(() => {
    if (root !== undefined) rmSync(root, { recursive: true, force: true });
    root = undefined;
  });

  const alice = machinePatterns({
    user: 'alice',
    home: '/home/alice',
    cargoHome: '/opt/cargo',
    rustupHome: '/home/alice/.rustup',
    checkouts: ['/srv/build/LOSAT', '/mnt/data/build/LOSAT'],
  });
  const labels = (text: string, patterns = alice) => findingsIn(text, patterns).map((finding) => finding.label);
  const findingsIn = (text: string, patterns = alice) => {
    root = mkdtempSync(join(tmpdir(), 'losat-td11-'));
    writeFileSync(join(root, 'file.js'), text);
    return scanFiles(root, ['file.js'], patterns);
  };

  it('finds a planted path in a served Wasm module and JavaScript file, and passes the clean ones', () => {
    root = mkdtempSync(join(tmpdir(), 'losat-td11-'));
    mkdirSync(join(root, 'assets'));
    mkdirSync(join(root, 'engine'));
    writeFileSync(join(root, 'index.html'), '<!doctype html><title>LOSAT Web</title>');
    writeFileSync(join(root, 'assets', 'clean.js'), 'const a = "/losat/LOSAT/src/main.rs"; // the remapped prefix of TD-11\n');
    writeFileSync(join(root, 'assets', 'planted.js'), 'throw new Error("at /srv/build/LOSAT/web/adapter/src/run.rs:10");\n');
    const wasm = Buffer.concat([Buffer.from([0, 0x61, 0x73, 0x6d, 1, 0, 0, 0]), Buffer.from('\0/opt/cargo/registry/src/index.crates.io-1/clap-4.6.0/src/lib.rs\0')]);
    writeFileSync(join(root, 'engine', 'losat-web-serial-0.wasm'), wasm);
    const files = filesUnder(root);
    expect(files).toEqual(['assets/clean.js', 'assets/planted.js', 'engine/losat-web-serial-0.wasm', 'index.html']);
    const found = scanFiles(root, files, alice);
    expect(found.map(({ file, label }) => `${file}: ${label}`)).toEqual(['assets/planted.js: the checkout', 'engine/losat-web-serial-0.wasm: CARGO_HOME']);
    expect(found[1]!.excerpt).toContain('/opt/cargo/registry/src/');
  });

  it("finds the home folders in Linux, macOS and Windows forms, escaped as in JSON, and the user's name", () => {
    expect(labels('x="/home/alice/src"')).toEqual(['HOME', 'the user name']);
    expect(labels('x="\\/home\\/alice\\/src"')).toEqual(['HOME', 'the user name']);
    expect(labels('x="/Users/alice/src"')).toEqual(['/Users/<user>', 'the user name']);
    expect(labels('x="C:\\\\Users\\\\alice\\\\src"')).toEqual(['\\Users\\<user>', 'the user name']);
    expect(labels('x="C:\\Users\\alice\\src"')).toEqual(['\\Users\\<user>', 'the user name']);
    expect(labels('built by alice on a laptop')).toEqual(['the user name']);
    expect(labels('x="/mnt/data/build/LOSAT/x"')).toEqual(['the checkout']);
  });

  it('does not take another name for the user or a longer path for a machine path', () => {
    expect(labels('x="/home/alicex/src"; malice; alice_b; /srv/build/LOSAT2')).toEqual([]);
  });

  it("checks a common user name (CI's runner) only in a path", () => {
    const runner = machinePatterns({ user: 'runner', home: '/home/runner', checkouts: ['/home/runner/work/LOSAT/LOSAT'] });
    expect(labels('the test runner runs; runner.run()', runner)).toEqual([]);
    expect(labels('x="/home/runner/work/LOSAT/LOSAT/LOSAT/src/main.rs"', runner)).toEqual(['the checkout', 'HOME', 'the user name in a path']);
    expect(labels('x="/tmp/runner/cache"', runner)).toEqual(['the user name in a path']);
    expect(labels('x="~runner"', runner)).toEqual(['the user name in a path']);
  });
});
