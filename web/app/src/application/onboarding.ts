// "Getting started" (S16): the shortest tutorial, a panel reachable from the header. It is
// open on the very first visit and closed on later visits once the user has closed it; that
// choice is the only thing kept (ports/local-settings.ts). Where the browser keeps nothing,
// the panel opens on each visit and stays closed in that tab after the user closes it.
import type { LocalSettings } from '../ports/local-settings';
import { Store } from './store';

export const GETTING_STARTED_KEY = 'losat-web:getting-started';
const CLOSED = 'closed';

export interface GettingStartedState {
  readonly open: boolean;
}

export class GettingStarted {
  readonly state: Store<GettingStartedState>;

  constructor(private readonly settings: LocalSettings) {
    this.state = new Store<GettingStartedState>({ open: settings.read(GETTING_STARTED_KEY) !== CLOSED });
  }

  open(): void {
    this.state.set({ open: true });
  }

  /** Closes the panel; later visits start with it closed, where the browser keeps the choice. */
  close(): void {
    this.state.set({ open: false });
    this.settings.write(GETTING_STARTED_KEY, CLOSED);
  }

  toggle(): void {
    if (this.state.get().open) this.close();
    else this.open();
  }
}
