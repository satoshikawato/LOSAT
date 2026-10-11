// Usage statistics in the browser (ports/analytics.ts): the choice in localStorage and the
// event posted to the site's endpoint. Every storage access may throw (storage disabled, a
// private window, a full quota): the page then asks again next visit and sends nothing.
import { parseAnalyticsChoice, serializeAnalyticsChoice } from '../../domain/analytics';
import type { AnalyticsPort } from '../../ports/analytics';

export const ANALYTICS_STORAGE_KEY = 'losat-web:analytics';

export const browserAnalytics: AnalyticsPort = {
  readChoice() {
    try {
      return parseAnalyticsChoice(localStorage.getItem(ANALYTICS_STORAGE_KEY));
    } catch {
      return undefined;
    }
  },

  writeChoice(choice) {
    try {
      if (choice === undefined) localStorage.removeItem(ANALYTICS_STORAGE_KEY);
      else localStorage.setItem(ANALYTICS_STORAGE_KEY, serializeAnalyticsChoice(choice));
      return true;
    } catch {
      return false;
    }
  },

  async sendVisit(endpoint, body) {
    try {
      await fetch(endpoint, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body,
        // Sent even when the page is closed right after it opens.
        keepalive: true,
        credentials: 'omit',
        cache: 'no-store',
        referrerPolicy: 'no-referrer',
      });
    } catch {
      // Blocked, offline or refused: nothing is retried or reported.
    }
  },
};
