// The site index from the network (ports/site.ts).
import { parseSiteIndex, SITE_INDEX_PATH } from '../../domain/site-index';
import type { SitePort } from '../../ports/site';

export { SITE_INDEX_PATH };
const TIMEOUT_MS = 10_000;

export const browserSite: SitePort = {
  async index() {
    try {
      const response = await fetch(SITE_INDEX_PATH, {
        cache: 'no-store',
        credentials: 'same-origin',
        signal: AbortSignal.timeout(TIMEOUT_MS),
      });
      if (!response.ok) return undefined;
      return parseSiteIndex(await response.json());
    } catch {
      return undefined;
    }
  },
};

/**
 * Shares a read of the site index among the callers that ask within `maxAgeMs` of it (the
 * measurement switch and the update check at start-up); a later call reads the network again.
 * Nothing is kept beyond the page.
 */
export function sharedSiteReads(port: SitePort, maxAgeMs: number, now: () => number): SitePort {
  let last: { readonly at: number; readonly read: ReturnType<SitePort['index']> } | undefined;
  return {
    index() {
      if (last === undefined || now() - last.at > maxAgeMs) last = { at: now(), read: port.index() };
      return last.read;
    },
  };
}
