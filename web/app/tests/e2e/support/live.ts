// What a screen reader is told (WA2-AX-03): the new text of each live region (aria-live, role
// status, alert or log) that was in the page before its text changed, and is rendered. A region
// inserted together with its text is not counted: screen readers usually do not read it.
import type { Page } from '@playwright/test';

export interface Spoken {
  /** "alert" for role=alert and aria-live=assertive, else "status". */
  readonly kind: 'alert' | 'status';
  readonly text: string;
}

declare global {
  interface Window {
    __spoken?: Spoken[];
  }
}

/** Starts recording; call it before the page opens (an init script). */
export async function listen(page: Page): Promise<void> {
  await page.addInitScript(() => {
    const spoken: Spoken[] = [];
    window.__spoken = spoken;
    const selector = '[aria-live], [role="status"], [role="alert"], [role="log"]';
    const last = new WeakMap<Element, string>();
    const textOf = (element: Element) => (element.textContent ?? '').replace(/\s+/g, ' ').trim();
    const scan = () => {
      for (const region of document.querySelectorAll(selector)) {
        const live = region.getAttribute('aria-live');
        if (live === 'off') continue;
        const text = textOf(region);
        const before = last.get(region);
        last.set(region, text);
        if (before === undefined || before === text || text === '') continue;
        if (region.closest('[aria-hidden="true"]') !== null || region.getClientRects().length === 0) continue;
        const kind = region.getAttribute('role') === 'alert' || live === 'assertive' ? 'alert' : 'status';
        spoken.push({ kind, text });
      }
    };
    new MutationObserver(scan).observe(document, { subtree: true, childList: true, characterData: true });
  });
}

/** What was told since the page opened (or since the last `forget`). */
export async function spoken(page: Page): Promise<Spoken[]> {
  return page.evaluate(() => window.__spoken ?? []);
}

/** The texts told, of one kind or of both. */
export async function told(page: Page, kind?: Spoken['kind']): Promise<string[]> {
  return (await spoken(page)).filter((entry) => kind === undefined || entry.kind === kind).map((entry) => entry.text);
}

export async function forget(page: Page): Promise<void> {
  await page.evaluate(() => window.__spoken?.splice(0));
}

/** The test ID of the element with the focus (its tag name when it has none). */
export async function focused(page: Page): Promise<string> {
  return page.evaluate(() => {
    const element = document.activeElement as HTMLElement | null;
    return element?.dataset['testid'] ?? element?.tagName ?? '';
  });
}
