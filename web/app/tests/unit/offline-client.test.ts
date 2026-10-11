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

  // Cookies or site data blocked for the site (fix round 2's probe): no offline use, as without
  // service workers, not a failure that "Try again" would repeat on every load (code review 2 L1).
  it.each([
    ['Chromium', new DOMException('Failed to register a ServiceWorker: The user denied permission to use Service Worker.', 'NotSupportedError')],
    ['Firefox', new DOMException('The operation is insecure.', 'SecurityError')],
  ])("treats %s's refusal of service workers for the site as service workers blocked", async (_, refusal) => {
    vi.stubGlobal('navigator', { serviceWorker: { register: () => Promise.reject(refusal) } });
    const offline = browserOffline();
    expect(offline.supported).toBe(true);
    await expect(offline.connect()).rejects.toThrow(`this browser does not let this site use service workers: ${refusal.message}`);
    expect(offline.supported).toBe(false);
  });

  it('passes on any other rejection, and still counts the browser as one with service workers', async () => {
    const failure = new TypeError('Failed to register a ServiceWorker: A bad HTTP response code (503) was received when fetching the script.');
    vi.stubGlobal('navigator', { serviceWorker: { register: () => Promise.reject(failure) } });
    const offline = browserOffline();
    await expect(offline.connect()).rejects.toBe(failure);
    expect(offline.supported).toBe(true);
  });

  it('is connected once the registration has an active worker', async () => {
    vi.stubGlobal('navigator', { serviceWorker: { register: () => Promise.resolve({ active: {} }) } });
    await expect(browserOffline().connect()).resolves.toBeUndefined();
  });
});
