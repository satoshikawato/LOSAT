// The search for a minimal reproduction of a failing random-walk step (walk/runner.ts), as
// gbdraw's random walk shrinks its failures (its tests/web/helpers/live-random-walk-shrink.mjs).
// A pure module: the walk passes a `reproduces` predicate that replays a candidate in a fresh
// page, and tests/unit/walk.test.ts drives it with a stub.
//
// A candidate is a list of step indices in walk order that ends with the failing step. The search
// tries, in order: the failing step alone; each earlier step with it, nearest first, until those
// replays reach half the budget (`pairPhaseBudget`), so ddmin always gets the other half; then
// ddmin (Zeller and Hildebrandt 2002) over the steps before it, with the failing step always
// appended. The whole walk is taken to reproduce (it failed in the walk itself) and is not
// replayed. Each predicate call on a new candidate is a replay; a candidate tried before costs
// none. The search stops at `budget` replays and returns the shortest reproducing candidate found.

export interface ShrinkResult {
  /** The shortest candidate that reproduced, or the whole walk. */
  readonly indices: readonly number[];
  /** The predicate calls made. */
  readonly replays: number;
  /** False when the search stopped at the budget. */
  readonly finished: boolean;
  /** True when `indices` reproduced in a replay (the whole walk is never replayed). */
  readonly confirmed: boolean;
}

/** The parts of `items` cut into `count` runs as even as possible. */
function split(items: readonly number[], count: number): number[][] {
  const parts: number[][] = [];
  let start = 0;
  for (let part = 0; part < count; part += 1) {
    const end = start + Math.floor((items.length - start) / (count - part));
    parts.push(items.slice(start, end));
    start = end;
  }
  return parts;
}

const OUT_OF_BUDGET = Symbol('out of budget');

/** The replays that the failing step alone and the pairs may use: half the budget, rounded up. */
export function pairPhaseBudget(budget: number): number {
  return Math.ceil(budget / 2);
}

/**
 * @param prefix the indices of the steps before the failing one, in walk order
 * @param failing the index of the failing step
 */
export async function shrinkWalk({
  prefix,
  failing,
  budget,
  reproduces,
  deadline,
}: {
  prefix: readonly number[];
  failing: number;
  budget: number;
  reproduces: (indices: readonly number[]) => Promise<boolean>;
  /** A time (Date.now()) after which no replay starts: the budget's end too. */
  deadline?: number;
}): Promise<ShrinkResult> {
  const tried = new Map<string, boolean>();
  let replays = 0;
  let best: readonly number[] = [...prefix, failing];
  let confirmed = false;
  /** Whether a subset of `prefix` (in walk order) with the failing step reproduces. */
  const test = async (subset: readonly number[]): Promise<boolean> => {
    const key = subset.join(',');
    const known = tried.get(key);
    if (known !== undefined) return known;
    if (replays >= budget || (deadline !== undefined && Date.now() >= deadline)) throw OUT_OF_BUDGET;
    replays += 1;
    const result = Boolean(await reproduces([...subset, failing]));
    tried.set(key, result);
    if (result && (!confirmed || subset.length + 1 < best.length)) {
      best = [...subset, failing];
      confirmed = true;
    }
    return result;
  };
  const done = (finished: boolean): ShrinkResult => ({ indices: best, replays, finished, confirmed });
  try {
    if (await test([])) return done(true);
    for (const index of [...prefix].reverse()) {
      if (replays >= pairPhaseBudget(budget)) break;
      if (await test([index])) return done(true);
    }
    // ddmin over the prefix: the whole prefix reproduces (the walk failed).
    let current: readonly number[] = prefix;
    let granularity = 2;
    while (current.length >= 2) {
      const parts = split(current, granularity);
      let reduced = false;
      for (const part of parts) {
        if (await test(part)) {
          current = part;
          granularity = 2;
          reduced = true;
          break;
        }
      }
      if (!reduced) {
        for (const part of parts) {
          const complement = current.filter((index) => !part.includes(index));
          if (await test(complement)) {
            current = complement;
            granularity = Math.max(granularity - 1, 2);
            reduced = true;
            break;
          }
        }
      }
      if (!reduced) {
        if (granularity >= current.length) break;
        granularity = Math.min(granularity * 2, current.length);
      }
    }
    if (!confirmed) best = [...current, failing];
    return done(true);
  } catch (error) {
    if (error === OUT_OF_BUDGET) return done(false);
    throw error;
  }
}
