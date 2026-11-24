// Back and Forward between the Search and Results tabs (line WQA's sweep, D6, 2026-11-19). NCBI
// BLAST's flow, which "Run LOSAT" follows (Owner, 2026-11-02), is a search page and a results page:
// Back on the results returns to the form as it was, and Forward shows the results again. Before,
// the tabs made no history entry, so Back on the results left LOSAT (the browser asked first, and
// leaving lost every run and input).
//
// The Results tab therefore has an entry of its own, pushed above the Search tab's entry when the
// Results tab is shown, with the same URL. There is at most one such entry: choosing the Search tab
// while the Results tab's entry is the current one goes back to the Search tab's entry, so
// switching tabs never piles up entries, Back from the Search tab leaves LOSAT as before (the
// page's leave guard asks first), and Forward from there shows the results again.
//
// A reload keeps nothing (REQ-20), so a Results entry of an earlier load of the page is one this
// page cannot show. Each entry names the load that pushed it. A load that starts on such an entry
// (a reload on the results, a discarded tab, a duplicated or restored tab) goes back over it to the
// Search entry below, and Forward onto one goes back over it again: the page shows the Search tab,
// and Back from there leaves LOSAT at once. Before (review of 6fff89265), the entry became a second
// Search entry, so the first Back after a reload did nothing, and Forward after a reload showed an
// empty Results tab.
//
// The page restores the scroll position of the form itself (AppView.vue): with the browser's own
// restoration, the Search tab's button, which goes back, scrolled the page down to where the form
// had been left, and the focused tab out of view.

export type MainTab = 'search' | 'results' | 'windowmasker';

type FormTab = 'search' | 'windowmasker';

/** The part of `window.history` that the tabs use. */
export interface TabHistory {
  readonly state: unknown;
  scrollRestoration: ScrollRestoration;
  pushState(data: unknown, unused: string): void;
  back(): void;
}

/**
 * A `back()` of the page whose `popstate` has not come after this long is taken as lost (a browser
 * drops or refuses history calls made too often, as with an arrow key held down): the next choice
 * and the next `popstate` (the user's own Back or Forward) are then not taken for its arrival.
 * A same-document traversal arrives within milliseconds.
 */
export const BACK_LOST_MS = 2000;

/** The state of a Results entry: the tab and the load of the page that pushed it, never research data (design §13.3). */
interface ResultsEntry {
  readonly losatTab: 'results';
  readonly load: string;
}

function isResultsEntry(state: unknown): state is ResultsEntry {
  return typeof state === 'object' && state !== null && (state as { losatTab?: unknown }).losatTab === 'results';
}

export interface TabNavigation {
  /** The user chose a tab: its button, the arrow keys, "Run LOSAT" or "Run a search". */
  choose(tab: MainTab): void;
  /** The browser has moved to another entry of the page (`popstate`): returns the tab to show. */
  arrived(): MainTab;
}

export interface TabNavigationOptions {
  /** This load of the page, different from every other load. */
  readonly load?: string;
  /** Milliseconds, monotonic. */
  readonly now?: () => number;
}

export function tabNavigation(history: TabHistory, options: TabNavigationOptions = {}): TabNavigation {
  const load = options.load ?? crypto.randomUUID();
  const now = options.now ?? (() => performance.now());
  history.scrollRestoration = 'manual';
  /** The form tab (Search or WindowMasker) that Back from Results returns to. */
  let lastForm: FormTab = 'search';
  /** Whether the current entry is this load's own Results entry. */
  const onResultsEntry = (): boolean => isResultsEntry(history.state) && history.state.load === load;
  /** A Results entry of an earlier load, whose results are gone. */
  const gone = (): boolean => isResultsEntry(history.state) && history.state.load !== load;
  // The tab chosen while a `back()` of the page is on its way: the entry is known only when the
  // browser has arrived, and a second choice must not act on the entry being left.
  let going: { tab: MainTab; since: number } | undefined;
  const pending = (): MainTab | undefined => {
    if (going !== undefined && now() - going.since > BACK_LOST_MS) going = undefined;
    return going?.tab;
  };
  const goBack = (tab: MainTab): void => {
    going = { tab, since: now() };
    try {
      history.back();
    } catch {
      going = undefined;
    }
  };
  const choose = (tab: MainTab): void => {
    if (tab !== 'results') lastForm = tab;
    if (pending() !== undefined) {
      going = { tab, since: going!.since };
      return;
    }
    if (tab === 'results') {
      if (onResultsEntry()) return;
      const entry: ResultsEntry = { losatTab: 'results', load };
      try {
        history.pushState(entry, '');
      } catch {
        // WebKit refuses more than 100 pushState calls in 30 s (an arrow key held down): the tab
        // is shown all the same, without an entry, and the next choice starts from the entry there is.
      }
    } else if (onResultsEntry() && !gone()) {
      // The entry below this load's Results entry is always the form entry it was pushed from.
      goBack(tab);
    }
  };
  // A load on a Results entry of an earlier load goes back over it (the page opens on the Search tab).
  if (gone()) goBack('search');
  return {
    choose,
    arrived() {
      const wanted = pending();
      going = undefined;
      if (gone()) {
        // Forward onto results that a reload took (or a lost back() arriving late): back over it.
        goBack(wanted ?? 'search');
        return wanted ?? 'search';
      }
      if (wanted === undefined) return onResultsEntry() ? 'results' : lastForm;
      choose(wanted);
      return wanted;
    },
  };
}
