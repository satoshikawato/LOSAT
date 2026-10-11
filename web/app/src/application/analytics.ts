// Usage statistics (plan §5.10 (b), docs/web/analytics.md). Measurement exists only when
// the site index names an endpoint; then the page asks once, and only with a stored consent
// does it post the one event: once per page load, and once at the moment of consent.
// Before consent nothing is sent, and withdrawing removes the client ID. A failure of the
// request or of storage is silent and changes nothing else in the application.
import { visitBody } from '../domain/analytics';
import { RELEASE_ID_PATTERN, type SiteIndex } from '../domain/site-index';
import type { AnalyticsPort } from '../ports/analytics';
import { Store } from './store';

export interface AnalyticsState {
  /** The site measures usage (its index names an endpoint): only then is anything shown. */
  readonly available: boolean;
  /** The event is sent: a consent is stored. */
  readonly on: boolean;
  /** The notice is shown: no choice is stored yet, or "Change" reopened it. */
  readonly asking: boolean;
  /** The last choice could not be stored: nothing is sent, and the next visit asks again. */
  readonly unsaved: boolean;
}

export interface AnalyticsDeps {
  /** The site index (`browserSite.index()`); undefined when it cannot be read. */
  readonly siteIndex: () => Promise<SiteIndex | undefined>;
  readonly port: AnalyticsPort;
  /** This build's release ID (`__LOSAT_RELEASE_ID__`); undefined where a build has none. */
  readonly release: string | undefined;
  /** A new random client ID (`crypto.randomUUID()`), made at consent. */
  readonly newClientId: () => string;
}

const OFF: AnalyticsState = Object.freeze({ available: false, on: false, asking: false, unsaved: false });

export class Analytics {
  readonly state = new Store<AnalyticsState>(OFF);
  private endpoint: string | undefined;

  constructor(private readonly deps: AnalyticsDeps) {}

  /** The release ID that the event names (the notice shows it). */
  get release(): string | undefined {
    return this.deps.release;
  }

  /** Reads the site index and the stored choice; sends the event of this page load when a consent is stored. */
  async start(): Promise<void> {
    const release = this.deps.release;
    if (release === undefined || !RELEASE_ID_PATTERN.test(release)) return;
    let index: SiteIndex | undefined;
    try {
      index = await this.deps.siteIndex();
    } catch {
      return;
    }
    const endpoint = index?.analytics?.endpoint;
    if (endpoint === undefined) return;
    this.endpoint = endpoint;
    const stored = this.deps.port.readChoice();
    this.state.set({ available: true, on: stored?.consent === 'granted', asking: stored === undefined, unsaved: false });
    if (stored?.consent === 'granted') this.send(stored.client);
  }

  /** "Allow": stores a consent with a new client ID and sends the event; nothing is sent when storage fails. */
  allow(): void {
    const current = this.state.get();
    if (!current.available) return;
    if (current.on) {
      this.state.set({ ...current, asking: false });
      return;
    }
    const client = this.deps.newClientId();
    const stored = this.deps.port.writeChoice({ consent: 'granted', client });
    this.state.set({ available: true, on: stored, asking: false, unsaved: !stored });
    if (stored) this.send(client);
  }

  /** "Don't allow", also to withdraw: nothing more is sent, and the stored client ID is removed. */
  deny(): void {
    if (!this.state.get().available) return;
    const stored = this.deps.port.writeChoice({ consent: 'denied' });
    // When the refusal cannot be stored, at least the consent and its client ID go.
    if (!stored) this.deps.port.writeChoice(undefined);
    this.state.set({ available: true, on: false, asking: false, unsaved: !stored });
  }

  /** The footer's "Change": shows the notice again. */
  change(): void {
    const current = this.state.get();
    if (current.available) this.state.set({ ...current, asking: true });
  }

  private send(client: string): void {
    if (this.endpoint === undefined || this.deps.release === undefined) return;
    void this.deps.port.sendVisit(this.endpoint, visitBody(this.deps.release, client));
  }
}
