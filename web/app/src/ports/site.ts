// The site index port (domain/site-index.ts). The browser implementation is
// src/infra/browser/site.ts; tests use a fake.
import type { SiteIndex } from '../domain/site-index';

export interface SitePort {
  /**
   * Reads `/releases.json` from the network, never from a cache; undefined when the page is
   * offline, the request fails or the index is not valid. Sends nothing but the request.
   */
  index(): Promise<SiteIndex | undefined>;
}
