import { describe, expect, it } from 'vitest';
import { parseHeadersFile, siteHeaders } from '../../build/headers';

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
});
