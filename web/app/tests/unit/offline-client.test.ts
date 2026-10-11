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
  // The browser then refuses getRegistration() as well (fix round 4's probe).
  it.each([
    ['Chromium', new DOMException('Failed to register a ServiceWorker: The user denied permission to use Service Worker.', 'NotSupportedError')],
    ['Firefox', new DOMException('The operation is insecure.', 'SecurityError')],
  ])("treats %s's refusal of service workers for the site as service workers blocked", async (_, refusal) => {
    vi.stubGlobal('navigator', { serviceWorker: { controller: null, register: () => Promise.reject(refusal), getRegistration: () => Promise.reject(refusal) } });
    const offline = browserOffline();
    expect(offline.supported).toBe(true);
    await expect(offline.connect()).rejects.toThrow(`this browser does not let this site use service workers: ${refusal.message}`);
    expect(offline.supported).toBe(false);
  });

  // The same SecurityError comes from a sw.js that is not JavaScript (a proxy's HTML page; Firefox
  // even with the same message), on a first visit or once the registration has gone, while the
  // browser answers getRegistration(): a failure to try again, and the update is still downloaded
  // and checked before a switch (code review 3 L2).
  it.each([
    ['Chromium, a first visit', new DOMException("Failed to register a ServiceWorker for scope ('http://localhost/') with script ('http://localhost/sw.js'): The script has an unsupported MIME type ('text/html').", 'SecurityError'), undefined],
    ['Firefox, a first visit', new DOMException('The operation is insecure.', 'SecurityError'), undefined],
    ['WebKit, a first visit', new DOMException('Script http://localhost/sw.js load failed', 'SecurityError'), undefined],
    ['Firefox, with a registration', new DOMException('The operation is insecure.', 'SecurityError'), { active: {} }],
  ])('keeps a SecurityError from a sw.js that is not JavaScript a failure (%s)', async (_, failure, registration) => {
    vi.stubGlobal('navigator', { serviceWorker: { controller: null, register: () => Promise.reject(failure), getRegistration: () => Promise.resolve(registration) } });
    const offline = browserOffline();
    await expect(offline.connect()).rejects.toBe(failure);
    expect(offline.supported).toBe(true);
  });

  it('never counts a page that a worker controls as one where service workers are blocked', async () => {
    const refusal = new DOMException('The operation is insecure.', 'SecurityError');
    vi.stubGlobal('navigator', { serviceWorker: { controller: {}, register: () => Promise.reject(refusal), getRegistration: () => Promise.reject(refusal) } });
    const offline = browserOffline();
    await expect(offline.connect()).rejects.toBe(refusal);
    expect(offline.supported).toBe(true);
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
