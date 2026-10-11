import { describe, expect, it } from 'vitest';
import { headersFor, matchPattern, parseHeaderRules, parseHeadersFile, readHeaderRules, siteHeaders } from '../../build/headers';

describe('_headers', () => {
  it('enables cross-origin isolation for every path', () => {
    const headers = siteHeaders();
    expect(headers['Cross-Origin-Opener-Policy']).toBe('same-origin');
    expect(headers['Cross-Origin-Embedder-Policy']).toBe('require-corp');
    expect(headers['Cross-Origin-Resource-Policy']).toBe('same-origin');
  });

  it('limits every fetch to the site itself and allows Wasm compilation', () => {
    const csp = siteHeaders()['Content-Security-Policy'] ?? '';
    for (const directive of [
      "default-src 'self'",
      "script-src 'self' 'wasm-unsafe-eval'",
      "worker-src 'self'",
      "connect-src 'self'",
      "object-src 'none'",
    ]) {
      expect(csp).toContain(directive);
    }
  });

  it('rejects a header line outside any path rule', () => {
    expect(() => parseHeadersFile('  X-Test: 1\n')).toThrow(/malformed/);
  });

  it("caches a release's files for a year and revalidates the site-level files", () => {
    const rules = readHeaderRules();
    const isolation = siteHeaders();
    expect(headersFor('/r/abc-fake/assets/index-1.js', rules)).toEqual({ ...isolation, 'Cache-Control': 'public, max-age=31536000, immutable' });
    for (const path of ['/', '/index.html', '/boot.js', '/boot.css', '/releases.json', '/404.html', '/sw.js']) {
      expect(headersFor(path, rules)).toEqual({ ...isolation, 'Cache-Control': 'no-cache' });
    }
    // Another path has the isolation headers and no Cache-Control of the site's (Pages adds its default).
    expect(headersFor('/other', rules)).toEqual(isolation);
  });
});

describe('_headers rules (Cloudflare Pages semantics)', () => {
  it('matches a splat against any characters, none included', () => {
    expect(matchPattern('/*', '/')).toEqual({ splat: '' });
    expect(matchPattern('/r/*', '/r/a/b/c.js')).toEqual({ splat: 'a/b/c.js' });
    expect(matchPattern('/r/*', '/r')).toBeUndefined();
    expect(matchPattern('/*.wasm', '/r/a/engine/x.wasm')).toEqual({ splat: 'r/a/engine/x' });
    expect(matchPattern('/', '/index.html')).toBeUndefined();
  });

  it('matches a placeholder against one path segment', () => {
    expect(matchPattern('/r/:release/release-manifest.json', '/r/abc/release-manifest.json')).toEqual({ release: 'abc' });
    expect(matchPattern('/r/:release/release-manifest.json', '/r/a/b/release-manifest.json')).toBeUndefined();
    expect(matchPattern('/r/:release/*', '/r/abc/assets/x.js')).toEqual({ release: 'abc', splat: 'assets/x.js' });
    // Characters that mean something in a regular expression are literal.
    expect(matchPattern('/a.b+(c)', '/a.b+(c)')).toEqual({});
    expect(matchPattern('/a.b', '/axb')).toBeUndefined();
  });

  it('applies every matching rule; a header of several rules is joined with ", "', () => {
    const rules = parseHeaderRules('/*\n  X-A: 1\n  Link: </a>\n/r/*\n  Link: </b>\n  X-B: 2\n');
    expect(headersFor('/r/x', rules)).toEqual({ 'X-A': '1', Link: '</a>, </b>', 'X-B': '2' });
    expect(headersFor('/x', rules)).toEqual({ 'X-A': '1', Link: '</a>' });
  });

  it('joins headers whose names differ only in case', () => {
    expect(headersFor('/x', parseHeaderRules('/*\n  Vary: A\n/x\n  vary: B\n'))).toEqual({ Vary: 'A, B' });
  });

  it('detaches a header with "! Name", whatever rule attached it', () => {
    const rules = parseHeaderRules("/*\n  Content-Security-Policy: default-src 'self'\n  X-A: 1\n/*.jpg\n  ! Content-Security-Policy\n");
    expect(headersFor('/a.jpg', rules)).toEqual({ 'X-A': '1' });
    expect(headersFor('/a.png', rules)).toEqual({ 'Content-Security-Policy': "default-src 'self'", 'X-A': '1' });
  });

  it('fills placeholders into the values of their own rule only', () => {
    const rules = parseHeaderRules('/r/:release/*\n  X-Release: :release\n  X-File: :splat\n/*\n  X-Other: :release\n');
    expect(headersFor('/r/abc/assets/x.js', rules)).toEqual({ 'X-Release': 'abc', 'X-File': 'assets/x.js', 'X-Other': ':release' });
  });

  it("keeps Cloudflare's limits and syntax", () => {
    expect(() => parseHeaderRules(`/*\n  X-A: ${'a'.repeat(2000)}\n`)).toThrow(/2000/);
    expect(() => parseHeaderRules(Array.from({ length: 101 }, (_, n) => `/p${n}\n  X: 1\n`).join(''))).toThrow(/100 rules/);
    expect(() => parseHeaderRules('/*/*\n  X: 1\n')).toThrow(/more than one splat/);
    expect(() => parseHeaderRules('https://example.org/*\n  X: 1\n')).toThrow(/path patterns/);
    expect(parseHeaderRules('# comment\n/*\n  # comment\n  X: 1\n')).toEqual([{ pattern: '/*', attach: [['X', '1']], detach: [] }]);
  });
});
