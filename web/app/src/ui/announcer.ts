// What a screen reader is told of changes that it does not see (WA2-AX-03): a run's status and
// its end, an Example that loads, the engine's check of an input and a search that cannot start.
// The texts go into two regions that are in the page from its start (LiveAnnouncer.vue): a
// polite status, and an alert for what failed. A region inserted together with its text is
// usually not read, and the texts on the screen are inserted with the parts that show them.
// Texts told in one turn of the event loop are told together; each is the screen's words. A text
// leaves its region 10 s after it was told, so that a screen reader that reads the page to its
// end later does not read an old one.
import { shallowRef, type ShallowRef } from 'vue';

export type Politeness = 'status' | 'alert';

/** The text of a region, and a number that changes with each telling (the same text told again is read again). */
export interface Told {
  readonly key: number;
  readonly text: string;
}

export const told: Readonly<Record<Politeness, ShallowRef<Told>>> = {
  status: shallowRef<Told>({ key: 0, text: '' }),
  alert: shallowRef<Told>({ key: 0, text: '' }),
};

const waiting: Record<Politeness, string[]> = { status: [], alert: [] };
let queued = false;
const KEEP_MS = 10_000;
const clearing: Partial<Record<Politeness, ReturnType<typeof setTimeout>>> = {};

/** Tells `text` (a sentence of the screen) to a screen reader; an empty text tells nothing. */
export function announce(text: string, politeness: Politeness = 'status'): void {
  const sentence = text.replace(/\s+/g, ' ').trim();
  if (sentence === '') return;
  waiting[politeness].push(sentence);
  if (queued) return;
  queued = true;
  queueMicrotask(() => {
    queued = false;
    for (const kind of ['status', 'alert'] as const) {
      if (waiting[kind].length === 0) continue;
      told[kind].value = { key: told[kind].value.key + 1, text: waiting[kind].join(' ') };
      waiting[kind] = [];
      clearTimeout(clearing[kind]);
      clearing[kind] = setTimeout(() => (told[kind].value = { key: told[kind].value.key + 1, text: '' }), KEEP_MS);
    }
  });
}
