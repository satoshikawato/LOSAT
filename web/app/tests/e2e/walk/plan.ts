// The kind schedule of the random walk of the search form and the results (random-walk.spec.ts,
// walk/runner.ts). A pure module: tests/unit/walk.test.ts checks its coverage without a
// browser. The method is gbdraw's random walk (its tests/web/helpers/live-random-walk-plan.mjs):
// a seeded greedy plan that covers every ordered pair of step kinds, because the defects such
// walks find are mostly two operations that meet.
//
// A plan is a list of walks, each a list of planned steps `{ kind }`, each walk from a freshly
// opened page. `planWalk` builds it greedily from a seed until every target of its level is
// covered:
//
// - level `pr`: every kind, and every ordered pair of kinds (A, B), A = B included, with A
//   before B in one walk and no barrier step (a reload, which starts the page again) strictly
//   between them;
// - level `full`: the `pr` targets, and every ordered pair (A, B) as adjacent steps.
//
// At each position the planner takes the kind that covers the most uncovered targets, then the
// one that leaves the most uncovered pairs to start, then a seeded draw. The arguments of a step
// (which field, which value, which record) are drawn when the walk runs, from the walk's own
// generator, and recorded, so a replay repeats them while runs with other seeds still differ.

/** The step kinds (walk/steps.ts runs them). */
export const WALK_KINDS = Object.freeze([
  // The search form: the mode, its task, the form's fields, Additional options, "Restore default
  // search parameters", the job title, the threads.
  'mode',
  'task',
  'field',
  'extra',
  'restore',
  'title',
  'threads',
  // The inputs: paste, open files, remove a source, include or leave out records, combined or
  // separate searches, a subrange.
  'paste',
  'open',
  'remove',
  'records',
  'combine',
  'region',
  // Runs: Run LOSAT (to its end), Run LOSAT twice in a row (the second waits in the queue and must
  // give the same bytes), a run cancelled as it starts.
  'run',
  'twice',
  'cancel',
  // Sweep 2 (not on dev): Run LOSAT without waiting, so that the next steps act while it runs; the
  // results' filters and sorting; accessions fetched from stubbed sources; the tab hidden and shown.
  'start',
  'filter',
  'fetch',
  'hide',
  // Sweep 3: a run that fails on a later batch (D10/D15: the outputs written before the error are kept), a
  // help or version word (D14), a -logfile name (D11).
  'failrun',
  'printed',
  'logfile',
  // Sweep 4: a local BLAST database as the subject (files picked, a filter or a mask), BLASTP's PSI-BLAST rounds, a search strategy file.
  'database',
  'psiround',
  'strategy',
  // The page: the Search and Results tabs, the views of the results, the window's width, a reload.
  'tab',
  'view',
  'width',
  'reload',
] as const);

export type WalkKind = (typeof WALK_KINDS)[number];

/** The kinds that start the page again: a pair across one of them is not an interaction. */
export const BARRIER_KINDS: readonly WalkKind[] = Object.freeze(['reload']);

export type PlanLevel = 'pr' | 'full';

export interface PlannedStep {
  readonly kind: WalkKind;
}

export interface PlannedWalk {
  readonly steps: readonly PlannedStep[];
}

export interface WalkPlan {
  readonly walks: readonly PlannedWalk[];
}

export const DEFAULT_WALK_LENGTH: Readonly<Record<PlanLevel, number>> = Object.freeze({ pr: 44, full: 60 });

/** The seed of the pull-request walk (random-walk.spec.ts). */
export const PR_SEED = 20261119;

/** The most walks a plan may hold; a plan that reaches it leaves targets uncovered. */
const MAX_WALKS = 40;
/** The greedy plans tried per seed; the one with the fewest missing targets, then walks, then steps wins. */
const ATTEMPTS: Readonly<Record<PlanLevel, number>> = Object.freeze({ pr: 16, full: 2 });

export interface Random {
  /** A number in [0, 1). */
  next(): number;
  /** An integer in [0, count). */
  int(count: number): number;
  pick<T>(items: readonly T[]): T;
}

/** mulberry32: a 32-bit seeded generator; the same seed gives the same draws. */
export function createRandom(seed: number): Random {
  let state = seed >>> 0;
  const next = (): number => {
    state = (state + 0x6d2b79f5) >>> 0;
    let value = state;
    value = Math.imul(value ^ (value >>> 15), value | 1);
    value ^= value + Math.imul(value ^ (value >>> 7), value | 61);
    return ((value ^ (value >>> 14)) >>> 0) / 4294967296;
  };
  const int = (count: number): number => Math.floor(next() * count);
  const pick = <T>(items: readonly T[]): T => {
    if (items.length === 0) throw new Error('pick from an empty list');
    return items[int(items.length)]!;
  };
  return { next, int, pick };
}

const pairKey = (first: WalkKind, second: WalkKind): string => `${first}>${second}`;
const isBarrier = (kind: WalkKind): boolean => BARRIER_KINDS.includes(kind);

type Group = 'kinds' | 'sequence' | 'adjacent';

/**
 * The coverage targets of a level, by group: `kinds` (each kind once), `sequence` (ordered pairs
 * as a subsequence within a barrier-free stretch) and, for `full`, `adjacent` (ordered pairs as
 * adjacent steps).
 */
export function coverageTargets(level: PlanLevel): Partial<Record<Group, Set<string>>> {
  if (level !== 'pr' && level !== 'full') throw new Error(`unknown plan level ${String(level)}`);
  const pairs = new Set<string>();
  for (const first of WALK_KINDS) for (const second of WALK_KINDS) pairs.add(pairKey(first, second));
  const targets: Partial<Record<Group, Set<string>>> = { kinds: new Set(WALK_KINDS), sequence: pairs };
  if (level === 'full') targets.adjacent = new Set(pairs);
  return targets;
}

/**
 * The targets that a step at the end of a walk covers, by group. `segment` holds the kinds since
 * the last barrier (the barrier included); `previous` is the kind before the step, or null at the
 * start of a walk.
 */
function coveredBy(segment: ReadonlySet<WalkKind>, previous: WalkKind | null, kind: WalkKind): Record<Group, string[]> {
  return {
    kinds: [kind],
    sequence: [...segment].map((first) => pairKey(first, kind)),
    adjacent: previous === null ? [] : [pairKey(previous, kind)],
  };
}

/** The segment after a step: a barrier starts a new one. */
function nextSegment(segment: ReadonlySet<WalkKind>, kind: WalkKind): Set<WalkKind> {
  return isBarrier(kind) ? new Set([kind]) : new Set([...segment, kind]);
}

/** Which targets of a level a plan covers, and which it misses. */
export function coverageOf(plan: WalkPlan, level: PlanLevel): { covered: Record<string, number>; missing: string[] } {
  const targets = coverageTargets(level);
  const hit = new Map<string, Set<string>>(Object.keys(targets).map((group) => [group, new Set<string>()]));
  for (const { steps } of plan.walks) {
    let segment = new Set<WalkKind>();
    let previous: WalkKind | null = null;
    for (const { kind } of steps) {
      if (!WALK_KINDS.includes(kind)) throw new Error(`unknown step kind ${String(kind)}`);
      for (const [group, keys] of Object.entries(coveredBy(segment, previous, kind))) {
        const set = hit.get(group);
        if (set !== undefined) for (const key of keys) set.add(key);
      }
      segment = nextSegment(segment, kind);
      previous = kind;
    }
  }
  const covered: Record<string, number> = {};
  const missing: string[] = [];
  for (const [group, keys] of Object.entries(targets)) {
    const set = hit.get(group)!;
    covered[group] = [...keys].filter((key) => set.has(key)).length;
    for (const key of keys) if (!set.has(key)) missing.push(`${group} ${key}`);
  }
  return { covered, missing };
}

/** The number of steps of a plan. */
export function stepCount(plan: WalkPlan): number {
  return plan.walks.reduce((sum, { steps }) => sum + steps.length, 0);
}

/**
 * The walk plan of a seed: walks of at most `walkLength` steps until every target of `level` is
 * covered (or MAX_WALKS is reached; `coverageOf` then names what is missing).
 */
/** The first step of every walk. */
export const OPENING: WalkKind = 'run';

export function planWalk(seed: number, { level, walkLength = DEFAULT_WALK_LENGTH[level] }: { level: PlanLevel; walkLength?: number }): WalkPlan {
  if (!Number.isInteger(walkLength) || walkLength < 2) throw new Error(`a walk needs at least 2 steps; got ${walkLength}`);
  let best: WalkPlan | undefined;
  let bestRank: number[] | undefined;
  for (let attempt = 0; attempt < ATTEMPTS[level]; attempt += 1) {
    const plan = greedyPlan(createRandom((seed ^ 0x5eed5eed) + Math.imul(attempt, 0x9e3779b1)), level, walkLength);
    const rank = [coverageOf(plan, level).missing.length, plan.walks.length, stepCount(plan)];
    const differs = bestRank === undefined ? -1 : rank.findIndex((value, index) => value !== bestRank![index]);
    if (bestRank === undefined || (differs >= 0 && rank[differs]! < bestRank[differs]!)) {
      best = plan;
      bestRank = rank;
    }
  }
  return best!;
}

/** One greedy plan from a generator. */
function greedyPlan(random: Random, level: PlanLevel, walkLength: number): WalkPlan {
  const uncovered = coverageTargets(level);
  const remaining = (): number => Object.values(uncovered).reduce((sum, keys) => sum + keys.size, 0);
  const gainOf = (keysByGroup: Partial<Record<Group, string[]>>): number =>
    Object.entries(keysByGroup).reduce((sum, [group, keys]) => {
      const open = uncovered[group as Group];
      return sum + (open === undefined ? 0 : keys.filter((key) => open.has(key)).length);
    }, 0);
  // The pairs a kind could still start: they break ties toward the kinds whose later steps cover more.
  const openPairsFrom = (kind: WalkKind): number =>
    WALK_KINDS.reduce((sum, second) => {
      const key = pairKey(kind, second);
      return sum + (uncovered.sequence?.has(key) ? 1 : 0) + (uncovered.adjacent?.has(key) ? 1 : 0);
    }, 0);

  const walks: PlannedWalk[] = [];
  while (remaining() > 0 && walks.length < MAX_WALKS) {
    const steps: PlannedStep[] = [];
    let segment = new Set<WalkKind>();
    let previous: WalkKind | null = null;
    const take = (kind: WalkKind): void => {
      for (const [group, keys] of Object.entries(coveredBy(segment, previous, kind))) {
        const open = uncovered[group as Group];
        if (open !== undefined) for (const key of keys) open.delete(key);
      }
      steps.push({ kind });
      segment = nextSegment(segment, kind);
      previous = kind;
    };
    // A walk starts with the mode's query and subject pasted (steps.ts `start`): its first step
    // runs them, so that every walk searches and compares at least once.
    take(OPENING);
    while (steps.length < walkLength && remaining() > 0) {
      const bestOf = (candidates: readonly WalkKind[]): { best: WalkKind[]; gain: number } => {
        let best: WalkKind[] = [];
        let bestScore = [-1, -1];
        for (const kind of candidates) {
          const covered = coveredBy(segment, previous, kind);
          const score = [gainOf(covered), openPairsFrom(kind)];
          if (score[0]! > bestScore[0]! || (score[0] === bestScore[0] && score[1]! > bestScore[1]!)) {
            best = [kind];
            bestScore = score;
          } else if (score[0] === bestScore[0] && score[1] === bestScore[1]) {
            best.push(kind);
          }
        }
        return { best, gain: bestScore[0]! };
      };
      // A barrier ends the pairs of the steps before it, so it goes at the start or the end of a
      // walk, or where no other step covers anything.
      const atEdge = steps.length === 0 || steps.length === walkLength - 1;
      const first = bestOf(atEdge ? WALK_KINDS : WALK_KINDS.filter((kind) => !isBarrier(kind)));
      const best = !atEdge && first.gain === 0 ? bestOf(WALK_KINDS).best : first.best;
      take(random.pick(best));
    }
    walks.push({ steps });
  }
  return { walks };
}
