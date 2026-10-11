// The build of the site's service worker (build/service-worker.ts): one classic script, the
// same bytes from the same sources.
import { describe, expect, it } from 'vitest';
import { buildServiceWorker } from '../../build/service-worker';

describe('the service worker build', () => {
  it('is one deterministic classic script with the worker’s event handlers', async () => {
    const [first, second] = [await buildServiceWorker(), await buildServiceWorker()];
    expect(first).toBe(second);
    expect(first).not.toMatch(/^\s*(?:import|export)\b/m);
    for (const event of ['activate', 'fetch', 'message']) expect(first).toMatch(new RegExp(`addEventListener\\([\`'"]${event}[\`'"]`));
    // No handler of `install`: nothing is stored, and no skipWaiting, before a page asks.
    expect(first).not.toMatch(/addEventListener\([`'"]install/);
  }, 60_000);
});
