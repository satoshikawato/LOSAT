// Back and Forward between the Search and Results tabs (ui/tab-history.ts; line WQA's sweep, D6).
import { describe, expect, it } from 'vitest';
import { BACK_LOST_MS, tabNavigation, type MainTab, type TabHistory, type TabNavigation } from '../../src/ui/tab-history';

/**
 * The document state of history entries (HTML's session history): pushState gives the new entry
 * the document state of the entry it was pushed from, and a reload puts a new document into the
 * current entry's document state, so the entries that share it stay one document with the new
 * load. Moving to an entry of another document state loads that page again; a restored or
 * duplicated tab has a document state of its own for every entry.
 */
interface DocumentState {
  page?: Page;
}
interface Entry {
  state: unknown;
  readonly document: DocumentState | 'another site';
}

/** One load of LOSAT in the tab: the tab shown, set by a choice at once and by `popstate`. */
class Page {
  shown: MainTab = 'search';
  readonly navigation: TabNavigation;
  constructor(tab: FakeTab, load: string) {
    this.navigation = tabNavigation(tab.history(this), { load, now: () => tab.clock });
  }
  choose(tab: MainTab): void {
    this.shown = tab;
    this.navigation.choose(tab);
  }
  popstate(): void {
    this.shown = this.navigation.arrived();
  }
}

/** A browser tab: another site's page, then LOSAT. `back()` of the page arrives later, as a browser's does. */
class FakeTab {
  readonly entries: Entry[] = [{ state: null, document: 'another site' }];
  index = 0;
  clock = 0;
  loads = 0;
  /** The page's back() calls are dropped (a browser's rate limit) or throw. */
  backs: 'arrive' | 'dropped' | 'throw' = 'arrive';
  private queued: Page[] = [];

  constructor() {
    this.entries.push({ state: null, document: {} });
    this.index = 1;
    this.load();
  }

  get page(): Page | undefined {
    const document = this.entries[this.index]!.document;
    return document === 'another site' ? undefined : document.page;
  }
  get left(): boolean {
    return this.entries[this.index]!.document === 'another site';
  }
  /** The page's own entries: their states. */
  own(): unknown[] {
    return this.entries.filter((entry) => entry.document !== 'another site').map((entry) => entry.state);
  }

  /** `window.history` as a page sees it. */
  history(page: Page): TabHistory & { replaceState(data: unknown, unused: string): void } {
    const current = (): Entry => this.entries[this.index]!;
    const push = (state: unknown): void => {
      const document = current().document;
      this.entries.splice(this.index + 1);
      this.entries.push({ state, document });
      this.index += 1;
    };
    const back = (): void => {
      if (this.backs === 'throw') throw new DOMException('Too many calls to Location or History APIs', 'SecurityError');
      if (this.backs === 'arrive') this.queued.push(page);
    };
    return {
      scrollRestoration: 'auto',
      get state() {
        return current().state;
      },
      pushState: (data: unknown) => push(structuredClone(data)),
      // Not used by the page now; the page before the review replaced the state of a Results entry on load.
      replaceState: (data: unknown) => {
        current().state = structuredClone(data);
      },
      back,
    };
  }

  private load(): void {
    const document = this.entries[this.index]!.document as DocumentState;
    document.page = new Page(this, `load-${++this.loads}`);
  }
  reload(): void {
    this.load();
  }
  /** A tab restored or duplicated: the same entries, no document left. */
  restore(): void {
    for (let i = 0; i < this.entries.length; i++) {
      if (this.entries[i]!.document !== 'another site') this.entries[i] = { state: this.entries[i]!.state, document: {} };
    }
    this.load();
  }
  /** The browser's Back (-1) and Forward (+1). */
  go(delta: number): void {
    const from = this.entries[this.index]!.document;
    const next = this.index + delta;
    if (next < 0 || next >= this.entries.length) return;
    this.index = next;
    const to = this.entries[next]!.document;
    if (to === 'another site') return;
    if (to === from && to.page !== undefined) to.page.popstate();
    else this.load();
  }
  /** Delivers the page's back() calls on their way (a page that is gone has none). */
  settle(): void {
    while (this.queued.length > 0) {
      const page = this.queued.shift()!;
      if (page === this.page) this.go(-1);
    }
  }
  /** "Run LOSAT": the results of a run. */
  run(): void {
    this.page!.choose('results');
  }
}

const results = (load: number) => ({ losatTab: 'results', load: `load-${load}` });

describe('Back and Forward between the Search and Results tabs', () => {
  it('Back on the results shows the form, Forward the results again, and Back from the form leaves', () => {
    const tab = new FakeTab();
    tab.run();
    expect(tab.own()).toEqual([null, results(1)]);
    tab.go(-1);
    expect(tab.page!.shown).toBe('search');
    tab.go(+1);
    expect(tab.page!.shown).toBe('results');
    tab.go(-1);
    tab.go(-1);
    expect(tab.left).toBe(true);
    expect(tab.loads).toBe(1);
  });

  it('the Search tab goes back to the form entry: tabs never pile up entries', () => {
    const tab = new FakeTab();
    for (let i = 0; i < 5; i++) {
      tab.page!.choose('results');
      tab.page!.choose('search');
      tab.settle();
      expect(tab.page!.shown).toBe('search');
    }
    expect(tab.own()).toEqual([null, results(1)]);
    expect(tab.index).toBe(1);
    tab.go(+1);
    expect(tab.page!.shown).toBe('results');
    tab.go(-1);
    tab.run();
    expect(tab.own()).toEqual([null, results(1)]);
  });

  it('a tab chosen while a Back is on its way is the tab when it arrives', () => {
    const tab = new FakeTab();
    tab.page!.choose('search');
    expect(tab.index).toBe(1); // the Search tab on the form's entry never leaves the page
    tab.run();
    tab.page!.choose('search');
    tab.page!.choose('results'); // Run LOSAT, or an arrow key, before the Back has arrived
    expect(tab.page!.shown).toBe('results');
    tab.settle();
    expect(tab.page!.shown).toBe('results');
    expect(tab.own()).toEqual([null, results(1)]);
    tab.go(-1);
    expect(tab.page!.shown).toBe('search');
  });

  // Review of 6fff89265: the first Back after a reload did nothing, and Forward showed empty results.
  it('after a reload on the results, the page shows the form and the first Back leaves', () => {
    const tab = new FakeTab();
    tab.run();
    tab.reload();
    tab.settle();
    expect(tab.page!.shown).toBe('search');
    expect(tab.index).toBe(1);
    expect(tab.loads).toBe(2);
    tab.go(-1);
    expect(tab.left).toBe(true);
  });

  it('after Back to the form and a reload, Forward does not open empty results, and Back leaves', () => {
    const tab = new FakeTab();
    tab.run();
    tab.go(-1);
    tab.reload();
    tab.go(+1);
    expect(tab.page!.shown).toBe('search');
    tab.settle();
    expect(tab.index).toBe(1);
    expect(tab.page!.shown).toBe('search');
    tab.go(-1);
    expect(tab.left).toBe(true);
  });

  it('runs and reloads leave no entry that Back stops at', () => {
    const tab = new FakeTab();
    for (let i = 0; i < 3; i++) {
      tab.run();
      tab.reload();
      tab.settle();
    }
    expect(tab.page!.shown).toBe('search');
    tab.go(-1);
    expect(tab.left).toBe(true);
  });

  it('a restored or duplicated tab on the results opens on the form, and Back leaves', () => {
    const tab = new FakeTab();
    tab.run();
    tab.restore();
    tab.settle();
    expect(tab.index).toBe(1);
    expect(tab.page!.shown).toBe('search');
    tab.go(-1);
    expect(tab.left).toBe(true);
  });

  it('a back() that the browser drops is not waited for, and the next Back is the user\'s', () => {
    const tab = new FakeTab();
    tab.run();
    tab.backs = 'dropped';
    tab.page!.choose('search');
    tab.clock += BACK_LOST_MS + 1;
    tab.page!.choose('results');
    expect(tab.own()).toHaveLength(2);
    tab.go(-1);
    expect(tab.page!.shown).toBe('search');
    tab.backs = 'arrive';
    tab.go(+1);
    expect(tab.page!.shown).toBe('results');
  });

  it('a back() that the browser refuses leaves the choices working', () => {
    const tab = new FakeTab();
    tab.run();
    tab.backs = 'throw';
    expect(() => tab.page!.choose('search')).not.toThrow();
    tab.backs = 'arrive';
    tab.page!.choose('results');
    tab.page!.choose('search');
    tab.settle();
    expect(tab.index).toBe(1);
    expect(tab.page!.shown).toBe('search');
  });

  it('a pushState that the browser refuses leaves the tab shown and the history as it was', () => {
    class RefusingTab extends FakeTab {
      override history(page: Page) {
        const history = super.history(page);
        history.pushState = () => {
          throw new DOMException('Attempt to use history.pushState() more than 100 times per 30 seconds', 'SecurityError');
        };
        return history;
      }
    }
    const tab = new RefusingTab();
    tab.page!.choose('results');
    expect(tab.page!.shown).toBe('results');
    expect(tab.own()).toEqual([null]);
    tab.page!.choose('search');
    tab.settle();
    expect(tab.index).toBe(1);
  });

  it('the page restores the scroll itself', () => {
    const tab = new FakeTab();
    const seen: TabHistory = { state: null, scrollRestoration: 'auto', pushState: () => undefined, back: () => undefined };
    tabNavigation(seen, { load: 'x' });
    expect(seen.scrollRestoration).toBe('manual');
    expect(tab.page!.shown).toBe('search');
  });

  it('switching between Search and WindowMasker stays on the form entry, and Back from Results returns to WindowMasker', () => {
    const tab = new FakeTab();
    tab.page!.choose('windowmasker');
    expect(tab.own()).toEqual([null]);
    expect(tab.index).toBe(1);
    tab.run();
    expect(tab.own()).toEqual([null, results(1)]);
    tab.go(-1);
    expect(tab.page!.shown).toBe('windowmasker');
    tab.go(+1);
    expect(tab.page!.shown).toBe('results');
    tab.page!.choose('search');
    tab.settle();
    expect(tab.page!.shown).toBe('search');
    tab.go(-1);
    expect(tab.left).toBe(true);
  });
});
