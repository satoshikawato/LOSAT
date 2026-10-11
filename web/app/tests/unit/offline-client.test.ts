// The service-worker client of the page (src/infra/browser/offline.ts): a browser without
// service workers, or one where they are blocked, gives a reason that the header can show
// ("Offline use is not available (...)"), never an error of the client's own code.
import { afterEach, describe, expect, it, vi } from 'vitest';
import { browserOffline } from '../../src/infra/browser/offline';

afterEach(() => {
  vi.unstubAllGlobals();
});

describe('connecting to the service worker', () => {
  it('says that the browser has none', async () => {
    vi.stubGlobal('navigator', {});
    const offline = browserOffline();
    expect(offline.supported).toBe(false);
    await expect(offline.connect()).rejects.toThrow('this browser has no service workers');
  });

  it("says that they are blocked when register() resolves with nothing (Playwright's serviceWorkers: 'block')", async () => {
    vi.stubGlobal('navigator', { serviceWorker: { register: () => Promise.resolve(undefined) } });
    const offline = browserOffline();
    expect(offline.supported).toBe(true);
    await expect(offline.connect()).rejects.toThrow(/^service workers are blocked in this browser$/);
    // From then on the page treats the browser as one without service workers (the update switches without a download).
    expect(offline.supported).toBe(false);
  });

  it("passes on the browser's refusal", async () => {
    const refusal = new DOMException('The user denied permission to use Service Worker.', 'SecurityError');
    vi.stubGlobal('navigator', { serviceWorker: { register: () => Promise.reject(refusal) } });
    await expect(browserOffline().connect()).rejects.toBe(refusal);
  });

  it('is connected once the registration has an active worker', async () => {
    vi.stubGlobal('navigator', { serviceWorker: { register: () => Promise.resolve({ active: {} }) } });
    await expect(browserOffline().connect()).resolves.toBeUndefined();
  });
});
