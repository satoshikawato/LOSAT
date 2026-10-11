// The usage statistics port (domain/analytics.ts, docs/web/analytics.md): the choice kept in
// the browser and the one request of the event. The browser implementation is
// src/infra/browser/analytics.ts; tests use a fake.
import type { AnalyticsChoice } from '../domain/analytics';

export interface AnalyticsPort {
  /** The stored choice; undefined when none is stored, it is not valid, or storage cannot be read. */
  readChoice(): AnalyticsChoice | undefined;
  /** Stores a choice, or removes it (undefined); false when storage fails. */
  writeChoice(choice: AnalyticsChoice | undefined): boolean;
  /**
   * Posts the event's body to the same-origin endpoint once: no retry, no queue, no
   * credentials, no referrer. Never rejects; a failure changes nothing.
   */
  sendVisit(endpoint: string, body: string): Promise<void>;
}
