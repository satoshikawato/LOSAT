// Back and Forward between the Search and Results tabs (ui/tab-history.ts; line WQA's sweep, D6).
import { describe, expect, it } from 'vitest';
import { entryTab, tabNavigation, type MainTab, type TabHistory, type TabNavigation } from '../../src/ui/tab-history';

/**
 * A tab's session history: the entries of another site before the page's first entry, the page's
 * entries, and a cursor. `back()` arrives later, as a browser's does (`popstate` is a task); the
 * user's Back and Forward (`go`) arrive at once. A move off the page's entries leaves the page.
 */
class FakeHistory implements TabHistory {
  readonly entries: unknown[];
  index: number;
  backs = 0;
  left = false;
  scrollRestoration: ScrollRestoration = 'auto';
  private queued = 0;

  constructor(
    states: unknown[] = [null],
    private readonly page = { first: 1 },
  ) {
    this.entries = ['another site', ...states];
    this.index = this.entries.length - 1;
  }

  /** Set by the test: what the page does on `popstate`. */
  onPopState: () => void = () => undefined;

  get state(): unknown {
    return this.entries[this.index];
  }
  pushState(data: unknown): void {
    this.entries.splice(this.index + 1);
    this.entries.push(structuredClone(data));
    this.index += 1;
  }
  replaceState(data: unknown): void {
    this.entries[this.index] = structuredClone(data);
  }
  back(): void {
    this.backs += 1;
    this.queued += 1;
  }
  /** Delivers the `back()` calls on their way. */
  settle(): void {
    for (; this.queued > 0; this.queued -= 1) this.go(-1);
  }
  /** The browser's Back (-1) and Forward (+1). */
  go(delta: number): void {
    const next = this.index + delta;
    if (next < 0 || next >= this.entries.length) return;
    this.index = next;
    if (next < this.page.first) {
      this.left = true;
      return;
    }
    this.onPopState();
  }
  /** The page's entries above the other site's: [state, …]. */
  ownEntries(): unknown[] {
    return this.entries.slice(this.page.first);
  }
}

/** The page: the tab shown, set by a choice at once and by `popstate` from `arrived()`. */
function page(history: FakeHistory): { navigation: TabNavigation; shown: () => MainTab; choose: (tab: MainTab) => void } {
  const navigation = tabNavigation(history);
  let shown: MainTab = 'search';
  history.onPopState = () => {
    shown = navigation.arrived();
  };
  return {
    navigation,
    shown: () => shown,
    choose: (tab) => {
      shown = tab;
      navigation.choose(tab);
    },
  };
}

const RESULTS = { losatTab: 'results' };

describe('Back and Forward between the Search and Results tabs', () => {
  it('Back on the results shows the form, Forward the results again, and Back from the form leaves', () => {
    const history = new FakeHistory();
    const app = page(history);
    app.choose('results'); // Run LOSAT
    expect(history.ownEntries()).toEqual([null, RESULTS]);
    history.go(-1);
    expect(app.shown()).toBe('search');
    expect(history.left).toBe(false);
    history.go(+1);
    expect(app.shown()).toBe('results');
    history.go(-1);
    history.go(-1);
    expect(history.left).toBe(true);
  });

  it('the Search tab goes back to the form entry: tabs never pile up entries', () => {
    const history = new FakeHistory();
    const app = page(history);
    for (let i = 0; i < 5; i++) {
      app.choose('results');
      app.choose('search');
      expect(app.shown()).toBe('search');
      history.settle();
      expect(app.shown()).toBe('search');
    }
    expect(history.ownEntries()).toEqual([null, RESULTS]);
    expect(history.index).toBe(1);
    // Forward shows the results, and the form's next run replaces them.
    history.go(+1);
    expect(app.shown()).toBe('results');
    history.go(-1);
    app.choose('results');
    expect(history.ownEntries()).toEqual([null, RESULTS]);
    expect(history.index).toBe(2);
  });

  it('a tab chosen again, or chosen while a Back is on its way, takes the browser to that tab', () => {
    const history = new FakeHistory();
    const app = page(history);
    app.choose('search');
    expect(history.backs).toBe(0); // the Search tab on the form's entry never leaves the page
    app.choose('results');
    app.choose('results');
    expect(history.ownEntries()).toEqual([null, RESULTS]);
    // Search, then Results (Run LOSAT, or the arrow keys) before the Back has arrived.
    app.choose('search');
    app.choose('results');
    expect(app.shown()).toBe('results');
    history.settle();
    expect(app.shown()).toBe('results');
    expect(history.backs).toBe(1);
    expect(history.ownEntries()).toEqual([null, RESULTS]);
    expect(history.index).toBe(2);
    history.go(-1);
    expect(app.shown()).toBe('search');
  });

  it('a page reloaded on the results opens on the Search tab, its entry a form entry', () => {
    const history = new FakeHistory([null, RESULTS]);
    const app = page(history);
    expect(app.shown()).toBe('search');
    expect(history.ownEntries()).toEqual([null, null]);
    // The page puts each tab back where it was scrolled (AppView.vue), not the browser.
    expect(history.scrollRestoration).toBe('manual');
    app.choose('search');
    expect(history.backs).toBe(0);
    app.choose('results');
    expect(history.ownEntries()).toEqual([null, null, RESULTS]);
  });

  it('a pushState that the browser refuses leaves the tab shown and the history as it was', () => {
    const history = new FakeHistory();
    const app = page(history);
    history.pushState = () => {
      throw new DOMException('Attempt to use history.pushState() more than 100 times per 30 seconds', 'SecurityError');
    };
    app.choose('results');
    expect(app.shown()).toBe('results');
    expect(history.ownEntries()).toEqual([null]);
    app.choose('search');
    expect(history.backs).toBe(0);
  });

  it('reads only its own state', () => {
    expect(entryTab(RESULTS)).toBe('results');
    for (const state of [null, undefined, 'results', { losatTab: 'search' }, { tab: 'results' }, 0]) {
      expect(entryTab(state)).toBe('search');
    }
  });
});
