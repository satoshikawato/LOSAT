// The race of the accession sources (line WF-3, docs/web/fetch_design.md §10). With "Fetch from
// Fastest", the first Fetch of a kind asks each source for one fixed public record
// (RACE_RECORDS, never the user's accession) and sends the accessions to the first source that
// answers among those that serve them; it does not wait for the others, whose answers complete
// the order in the background. A source that has not answered within RACE_DEADLINE_MS counts as
// failed. The order is kept for the tab, so later fetches of the kind ask one source only. A
// race in which no source answered is not kept: the next Fetch races again.
import { FETCH_SOURCES, RACE_RECORDS, type FetchChoice, type FetchSourceId } from '../domain/accession';
import type { SequenceKind } from '../domain/programs';
import type { AccessionFetcher } from '../ports/fetch';

/** Bytes of a race record's answer at most (the records are about 200 bytes). */
const RACE_MAX_BYTES = 16_384;
/** Milliseconds that a source has to answer the race (ENA took 20-25 s at times on 2026-11-03). */
export const RACE_DEADLINE_MS = 10_000;

/** What the race of a kind found so far. */
export interface RaceStanding {
  /** The sources whose record came, first first, with the milliseconds that each took. */
  readonly answered: ReadonlyArray<{ readonly source: FetchSourceId; readonly ms: number }>;
  /** The sources that brought no record (no answer in time, an HTTP error, no such record). */
  readonly failed: readonly FetchSourceId[];
}

/**
 * Where accessions that `eligible` serve go, by the standing: the first of `eligible` that
 * answered; the first of `eligible` once every one of them failed; undefined while the race
 * has yet to say.
 */
export function raceChoice(standing: RaceStanding | undefined, eligible: readonly FetchSourceId[]): FetchSourceId | undefined {
  const first = standing?.answered.find((entry) => eligible.includes(entry.source))?.source;
  if (first !== undefined || standing === undefined) return first;
  return eligible.every((source) => standing.failed.includes(source)) ? eligible[0] : undefined;
}

/**
 * Where a group of a box that no one source serves starts (line WF-5): the source chosen by hand
 * if it serves the group, else the race's choice among the group's sources, else the first of
 * them. A group never starts a race.
 */
export function groupSite(
  sources: readonly FetchSourceId[],
  choice: FetchChoice,
  standing: RaceStanding | undefined,
): FetchSourceId {
  if (choice !== 'fastest' && sources.includes(choice)) return choice;
  return raceChoice(standing, sources) ?? sources[0]!;
}

interface Race {
  standing: RaceStanding;
  readonly waiters: Set<() => void>;
}

export interface SourceRaceDeps {
  readonly fetcher: AccessionFetcher;
  readonly now: () => number;
  /** Called with each change of a kind's standing; undefined when a race without answers is dropped. */
  readonly onChange: (kind: SequenceKind, standing: RaceStanding | undefined) => void;
  /** RACE_DEADLINE_MS by default (tests). */
  readonly deadlineMs?: number;
}

export class SourceRace {
  private readonly races = new Map<SequenceKind, Race>();
  private readonly pending = new Set<Promise<unknown>>();

  constructor(private readonly deps: SourceRaceDeps) {}

  /**
   * The source to send accessions of `kind` to, among `eligible` (not empty): see raceChoice;
   * the sources are raced if no race of the kind is kept. Rejects only when `signal` aborts;
   * the race goes on without this fetch.
   */
  async fastest(kind: SequenceKind, eligible: readonly FetchSourceId[], signal: AbortSignal): Promise<FetchSourceId> {
    const race = this.races.get(kind) ?? this.start(kind);
    for (;;) {
      signal.throwIfAborted();
      const choice = raceChoice(race.standing, eligible);
      if (choice !== undefined) return choice;
      await new Promise<void>((resolve, reject) => {
        const wake = () => {
          signal.removeEventListener('abort', abort);
          resolve();
        };
        const abort = () => {
          race.waiters.delete(wake);
          reject(signal.reason as Error);
        };
        race.waiters.add(wake);
        signal.addEventListener('abort', abort, { once: true });
      });
    }
  }

  /** Resolves when no race request is running (tests). */
  async idle(): Promise<void> {
    while (this.pending.size > 0) await Promise.allSettled([...this.pending]);
  }

  private start(kind: SequenceKind): Race {
    const race: Race = { standing: { answered: [], failed: [] }, waiters: new Set() };
    this.races.set(kind, race);
    this.deps.onChange(kind, race.standing);
    for (const { id: source } of FETCH_SOURCES) {
      // The race requests are not the fetch's: a fetch that stops leaves them to finish.
      const deadline = new AbortController();
      const timer = setTimeout(
        () => deadline.abort(new DOMException('no answer in time', 'TimeoutError')),
        this.deps.deadlineMs ?? RACE_DEADLINE_MS,
      );
      // From the moment the request goes out: NCBI's spacing between requests is not NCBI's time.
      let started = this.deps.now();
      const work = this.deps.fetcher
        .fetch(
          {
            source,
            kind,
            accession: RACE_RECORDS[kind][source],
            maxBytes: RACE_MAX_BYTES,
            fresh: true,
            onStart: () => (started = this.deps.now()),
          },
          deadline.signal,
        )
        .then(
          (result) => (result.ok ? this.deps.now() - started : undefined),
          () => undefined,
        )
        .then((ms) => {
          clearTimeout(timer);
          this.settle(kind, race, source, ms);
        });
      this.pending.add(work);
      void work.finally(() => this.pending.delete(work));
    }
    return race;
  }

  private settle(kind: SequenceKind, race: Race, source: FetchSourceId, ms: number | undefined): void {
    const { answered, failed } = race.standing;
    race.standing =
      ms === undefined ? { answered, failed: [...failed, source] } : { answered: [...answered, { source, ms }], failed };
    const lost = race.standing.failed.length === FETCH_SOURCES.length;
    if (lost && this.races.get(kind) === race) this.races.delete(kind);
    this.deps.onChange(kind, lost ? undefined : race.standing);
    for (const wake of [...race.waiters]) wake();
    race.waiters.clear();
  }
}
