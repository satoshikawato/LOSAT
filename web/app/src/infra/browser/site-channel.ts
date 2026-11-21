// The BroadcastChannel of the site's tabs, or undefined where there is none: a browser without
// BroadcastChannel, or one that refuses to open it. Firefox refuses it with SecurityError ("The
// operation is insecure.") where the site's data is blocked (cookies blocked), and the page did not
// start at all: the composition root opens the channels (found with the review of line WQA's sweep,
// D5, M2, 2026-11-19).
// Without a channel, no other tab is told.
export function siteChannel(name: string): BroadcastChannel | undefined {
  if (typeof BroadcastChannel !== 'function') return undefined;
  try {
    return new BroadcastChannel(name);
  } catch {
    return undefined;
  }
}
