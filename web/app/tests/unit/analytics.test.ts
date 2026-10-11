import { describe, expect, it } from 'vitest';
import { Analytics } from '../../src/application/analytics';
import {
  parseAnalyticsChoice,
  parseVisitBody,
  serializeAnalyticsChoice,
  visitBody,
  type AnalyticsChoice,
} from '../../src/domain/analytics';
import type { SiteIndex } from '../../src/domain/site-index';
import type { AnalyticsPort } from '../../src/ports/analytics';

const RELEASE = '0123456789ab-fake';
const CLIENT = '3f2b8c1e-5d4a-4b7c-9e1f-0a2b3c4d5e6f';

describe('the visit event', () => {
  it('is exactly the version, the release ID and the client ID', () => {
    expect(visitBody(RELEASE, CLIENT)).toBe(`{"v":1,"release":"${RELEASE}","client":"${CLIENT}"}`);
    expect(parseVisitBody(visitBody(RELEASE, CLIENT))).toEqual({ v: 1, release: RELEASE, client: CLIENT });
  });

  it('refuses IDs that do not match their patterns', () => {
    expect(() => visitBody('Release', CLIENT)).toThrow(/release/);
    expect(() => visitBody(RELEASE, CLIENT.toUpperCase())).toThrow(/client/);
    // Version 1 UUIDs and the nil UUID are not crypto.randomUUID()'s.
    expect(() => visitBody(RELEASE, '3f2b8c1e-5d4a-1b7c-9e1f-0a2b3c4d5e6f')).toThrow(/client/);
    expect(() => visitBody(RELEASE, '00000000-0000-0000-0000-000000000000')).toThrow(/client/);
  });

  it('is read only in the exact form that the page writes', () => {
    const exact = visitBody(RELEASE, CLIENT);
    for (const text of [
      `{"v":1,"release":"${RELEASE}","client":"${CLIENT}","query":"ACGT"}`,
      `{"v":2,"release":"${RELEASE}","client":"${CLIENT}"}`,
      `{"client":"${CLIENT}","release":"${RELEASE}","v":1}`,
      `{"v":1,"release":"${RELEASE}","client":"${CLIENT}","client":"${CLIENT}"}`,
      `{"v": 1,"release":"${RELEASE}","client":"${CLIENT}"}`,
      `${exact}\n`,
      `[${exact}]`,
      'null',
      '',
      '{',
    ]) {
      expect(parseVisitBody(text), text).toBeUndefined();
    }
  });
});

describe('the stored choice', () => {
  it('keeps a consent with its client ID, or a refusal without one', () => {
    for (const choice of [{ consent: 'granted', client: CLIENT }, { consent: 'denied' }] as AnalyticsChoice[]) {
      expect(parseAnalyticsChoice(serializeAnalyticsChoice(choice))).toEqual(choice);
    }
    expect(serializeAnalyticsChoice({ consent: 'denied' })).not.toContain('client');
  });

  it('reads anything else as no choice, so the page asks again', () => {
    for (const text of [null, '', 'granted', '{"v":1,"consent":"granted"}', '{"v":1,"consent":"granted","client":"x"}', '{"v":2,"consent":"denied"}', '[]']) {
      expect(parseAnalyticsChoice(text), String(text)).toBeUndefined();
    }
  });
});

interface FakePort extends AnalyticsPort {
  stored: AnalyticsChoice | undefined;
  /** Makes every write fail, as a storage that throws. */
  failWrites: boolean;
  readonly sent: { endpoint: string; body: string }[];
}

/** A port that keeps the choice in memory and records the requests. */
function fakePort(initial?: AnalyticsChoice): FakePort {
  const port: FakePort = {
    stored: initial,
    failWrites: false,
    sent: [],
    readChoice: () => port.stored,
    writeChoice: (choice) => {
      if (port.failWrites) return false;
      port.stored = choice;
      return true;
    },
    sendVisit: async (endpoint, body) => {
      port.sent.push({ endpoint, body });
    },
  };
  return port;
}

const siteIndex = (analytics?: { endpoint: string }): SiteIndex => ({
  format: 1,
  latest: RELEASE,
  releases: [{ id: RELEASE, version: '0.0.0', build: '0123456' }],
  ...(analytics === undefined ? {} : { analytics }),
});

/** The application of a site that measures, whose consent makes the client ID CLIENT. */
function analytics(port: FakePort) {
  return new Analytics({ siteIndex: async () => siteIndex({ endpoint: '/api/visit' }), port, release: RELEASE, newClientId: () => CLIENT });
}

describe('Analytics', () => {
  it('measures nothing, shows nothing and reads no storage when the site names no endpoint', async () => {
    for (const index of [siteIndex(), undefined]) {
      const port = fakePort({ consent: 'granted', client: CLIENT });
      let reads = 0;
      const counting = { ...port, readChoice: () => (reads++, port.readChoice()) };
      const a = new Analytics({ siteIndex: async () => index, port: counting, release: RELEASE, newClientId: () => CLIENT });
      await a.start();
      expect(a.state.get()).toEqual({ available: false, on: false, asking: false, unsaved: false });
      expect(reads).toBe(0);
      a.allow();
      expect(port.sent).toEqual([]);
    }
  });

  it('measures nothing in a build without a valid release ID', async () => {
    for (const release of [undefined, 'Not valid']) {
      const port = fakePort({ consent: 'granted', client: CLIENT });
      const a = new Analytics({ siteIndex: async () => siteIndex({ endpoint: '/api/visit' }), port, release, newClientId: () => CLIENT });
      await a.start();
      expect(a.state.get().available).toBe(false);
      expect(port.sent).toEqual([]);
    }
  });

  it('asks before anything is sent, and sends nothing after "Don\'t allow"', async () => {
    const port = fakePort();
    const a = analytics(port);
    await a.start();
    expect(a.state.get()).toEqual({ available: true, on: false, asking: true, unsaved: false });
    expect(port.sent).toEqual([]);
    a.deny();
    expect(a.state.get()).toEqual({ available: true, on: false, asking: false, unsaved: false });
    expect(port.stored).toEqual({ consent: 'denied' });
    // The next page load: no notice, nothing sent.
    const next = analytics(port);
    await next.start();
    expect(next.state.get()).toEqual({ available: true, on: false, asking: false, unsaved: false });
    expect(port.sent).toEqual([]);
  });

  it('sends one event at consent and one per page load after it, with the client ID made at consent', async () => {
    const port = fakePort();
    const a = analytics(port);
    await a.start();
    a.allow();
    expect(a.state.get()).toEqual({ available: true, on: true, asking: false, unsaved: false });
    expect(port.stored).toEqual({ consent: 'granted', client: CLIENT });
    expect(port.sent).toEqual([{ endpoint: '/api/visit', body: visitBody(RELEASE, CLIENT) }]);
    // "Allow" again from the reopened notice sends nothing more and keeps the client ID.
    a.change();
    expect(a.state.get().asking).toBe(true);
    a.allow();
    expect(port.sent).toHaveLength(1);
    const next = analytics(port);
    await next.start();
    expect(next.state.get()).toEqual({ available: true, on: true, asking: false, unsaved: false });
    expect(port.sent).toEqual([
      { endpoint: '/api/visit', body: visitBody(RELEASE, CLIENT) },
      { endpoint: '/api/visit', body: visitBody(RELEASE, CLIENT) },
    ]);
  });

  it('withdrawing removes the client ID and sends nothing more; a new consent makes a new client ID', async () => {
    const port = fakePort({ consent: 'granted', client: CLIENT });
    const a = analytics(port);
    await a.start();
    expect(port.sent).toHaveLength(1);
    a.change();
    a.deny();
    expect(port.stored).toEqual({ consent: 'denied' });
    expect(a.state.get().on).toBe(false);
    const second = '9a8b7c6d-5e4f-4a3b-8c2d-1e0f9a8b7c6d';
    const next = new Analytics({ siteIndex: async () => siteIndex({ endpoint: '/api/visit' }), port, release: RELEASE, newClientId: () => second });
    await next.start();
    expect(port.sent).toHaveLength(1);
    next.change();
    next.allow();
    expect(port.stored).toEqual({ consent: 'granted', client: second });
    expect(port.sent).toHaveLength(2);
    expect(port.sent[1]!.body).toBe(visitBody(RELEASE, second));
  });

  it('when storage fails, sends nothing, says so, and asks again next visit', async () => {
    const port = fakePort();
    port.failWrites = true;
    const a = analytics(port);
    await a.start();
    a.allow();
    expect(a.state.get()).toEqual({ available: true, on: false, asking: false, unsaved: true });
    expect(port.sent).toEqual([]);
    a.change();
    a.deny();
    expect(a.state.get()).toEqual({ available: true, on: false, asking: false, unsaved: true });
    const next = analytics(port);
    await next.start();
    expect(next.state.get().asking).toBe(true);
    expect(port.sent).toEqual([]);
  });

  it('a site index that cannot be read leaves measurement off', async () => {
    const port = fakePort({ consent: 'granted', client: CLIENT });
    const a = new Analytics({ siteIndex: () => Promise.reject(new Error('offline')), port, release: RELEASE, newClientId: () => CLIENT });
    await a.start();
    expect(a.state.get().available).toBe(false);
    expect(port.sent).toEqual([]);
  });
});
