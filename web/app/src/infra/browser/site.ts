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
