import { describe, expect, it } from 'vitest';
import { parseHeadersFile, pathHeaders, siteHeaders } from '../../build/headers';

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

  it('lets browsers keep the engine modules, whose names carry their digest, and nothing else', () => {
    const engine = pathHeaders('/engine/losat-web-threads-0123456789abcdef.wasm');
    expect(engine['Cache-Control']).toBe('public, max-age=31536000, immutable');
    expect(engine['Cross-Origin-Embedder-Policy']).toBe('require-corp');
    expect(pathHeaders('/index.html')['Cache-Control']).toBeUndefined();
    expect(pathHeaders('/assets/index.js')['Cache-Control']).toBeUndefined();
  });

  it('applies the rules that match a path in file order: exact paths and prefixes before a final *', () => {
    const rules = parseHeadersFile('/*\n  A: 1\n  B: 1\n/x\n  B: 2\n/y/*\n  C: 3\n');
    expect(pathHeaders('/x', rules)).toEqual({ A: '1', B: '2' });
    expect(pathHeaders('/x/z', rules)).toEqual({ A: '1', B: '1' });
    expect(pathHeaders('/y/z', rules)).toEqual({ A: '1', B: '1', C: '3' });
  });

  it('rejects a header line outside any path rule', () => {
    expect(() => parseHeadersFile('  X-Test: 1\n')).toThrow(/malformed/);
  });
});
