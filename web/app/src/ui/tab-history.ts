// Back and Forward between the Search and Results tabs (line WQA's sweep, D6, 2026-11-19). NCBI
// BLAST's flow, which "Run LOSAT" follows (Owner, 2026-11-02), is a search page and a results page:
// Back on the results returns to the form as it was, and Forward shows the results again. Before,
// the tabs made no history entry, so Back on the results left LOSAT (the browser asked first, and
// leaving lost every run and input).
//
// The Results tab therefore has an entry of its own, pushed above the Search tab's entry when the
// Results tab is shown, with the same URL: a reload opens the page as a first visit does, on the
// Search tab, and keeps nothing (REQ-20). There is at most one such entry: choosing the Search tab
// while the Results tab's entry is the current one goes back to the Search tab's entry, so
// switching tabs never piles up entries, Back from the Search tab leaves LOSAT as before (the
// page's leave guard asks first), and Forward from there shows the results again.
//
// The page restores the scroll positions itself (AppView.vue): with the browser's own restoration,
// the Search tab's button, which goes back, scrolled the page down to where the form had been left,
// and the focused tab out of view.

export type MainTab = 'search' | 'results';

/** The part of `window.history` that the tabs use. */
export interface TabHistory {
  readonly state: unknown;
  scrollRestoration: ScrollRestoration;
  pushState(data: unknown, unused: string): void;
  replaceState(data: unknown, unused: string): void;
  back(): void;
}

/** The state of the Results tab's entry. Only the tab: never research data (design §13.3). */
const RESULTS_ENTRY = { losatTab: 'results' } as const;

/** The tab of a history entry: the Results tab's entry, or any other (the page's own entry). */
export function entryTab(state: unknown): MainTab {
  return typeof state === 'object' && state !== null && (state as { losatTab?: unknown }).losatTab === 'results'
    ? 'results'
    : 'search';
}

export interface TabNavigation {
  /** The user chose a tab: its button, the arrow keys, "Run LOSAT" or "Run a search". */
  choose(tab: MainTab): void;
  /** The browser has moved to another entry of the page (`popstate`): returns the tab to show. */
  arrived(): MainTab;
}

export function tabNavigation(history: TabHistory): TabNavigation {
  // A page opens on the Search tab, also when it is reloaded on the Results tab's entry (the
  // browser keeps an entry's state across a reload): the entry becomes the Search tab's.
  if (entryTab(history.state) === 'results') history.replaceState(null, '');
  history.scrollRestoration = 'manual';
  // The tab chosen while a `back()` of the Search tab is on its way: the entry is known only when
  // the browser has arrived, and a second choice must not act on the entry being left.
  let going: MainTab | undefined;
  const choose = (tab: MainTab): void => {
    if (going !== undefined) {
      going = tab;
      return;
    }
    if (tab === entryTab(history.state)) return;
    if (tab === 'results') {
      try {
        history.pushState(RESULTS_ENTRY, '');
      } catch {
        // WebKit refuses more than 100 pushState calls in 30 s (an arrow key held down): the tab
        // is shown all the same, without an entry, and the next choice starts from the entry there is.
      }
    } else {
      // The entry below the Results tab's is always this page's Search entry: the page pushed the
      // Results entry from it, or came up to the Results entry from it (a reload made it a Search entry).
      going = tab;
      history.back();
    }
  };
  return {
    choose,
    arrived() {
      const wanted = going;
      going = undefined;
      if (wanted === undefined) return entryTab(history.state);
      choose(wanted);
      return wanted;
    },
  };
}
