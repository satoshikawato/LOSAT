// The race of the accession sources (line WF-3, docs/web/fetch_design.md §10). With "Fetch from
// Fastest", the first Fetch of a kind asks each source for one fixed public record
// (RACE_RECORDS, never the user's accession) and sends the accessions to the first source that
// answers among those that serve them; it does not wait for the others, whose answers complete
// the order in the background. The order is kept for the tab, so later fetches of the kind ask
// one source only. A race in which no source answered is not kept: the next Fetch races again.
import { FETCH_SOURCES, RACE_RECORDS, type FetchSourceId } from '../domain/accession';
import type { SequenceKind } from '../domain/programs';
import type { AccessionFetcher } from '../ports/fetch';

/** Bytes of a race record's answer at most (the records are about 200 bytes). */
const RACE_MAX_BYTES = 16_384;

/** What the race of a kind found so far. */
export interface RaceStanding {
  /** The sources whose record came, first first, with the milliseconds that each took. */
  readonly answered: ReadonlyArray<{ readonly source: FetchSourceId; readonly ms: number }>;
  /** The sources that brought no record (no answer, an HTTP error, no such record). */
  readonly failed: readonly FetchSourceId[];
}

/** The first of `eligible` in the standing's order, or undefined while none of them has answered. */
export function firstAnswered(standing: RaceStanding | undefined, eligible: readonly FetchSourceId[]): FetchSourceId | undefined {
  return standing?.answered.find((entry) => eligible.includes(entry.source))?.source;
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
}

export class SourceRace {
  private readonly races = new Map<SequenceKind, Race>();
  private readonly pending = new Set<Promise<unknown>>();

  constructor(private readonly deps: SourceRaceDeps) {}

  /**
   * The source to send accessions of `kind` to: the first of `eligible` whose record came,
   * racing the sources if no race of the kind is kept; NCBI once every one of `eligible` failed.
   * Rejects only when `signal` aborts; the race goes on without this fetch.
   */
  async fastest(kind: SequenceKind, eligible: readonly FetchSourceId[], signal: AbortSignal): Promise<FetchSourceId> {
    const race = this.races.get(kind) ?? this.start(kind);
    for (;;) {
      signal.throwIfAborted();
      const first = firstAnswered(race.standing, eligible);
      if (first !== undefined) return first;
      if (eligible.every((source) => race.standing.failed.includes(source))) return 'ncbi';
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
    // The race requests are not the fetch's: a fetch that stops leaves them to finish.
    const signal = new AbortController().signal;
    for (const { id: source } of FETCH_SOURCES) {
      const started = this.deps.now();
      const work = this.deps.fetcher
        .fetch({ source, kind, accession: RACE_RECORDS[kind][source], maxBytes: RACE_MAX_BYTES, fresh: true }, signal)
        .then(
          (result) => (result.ok ? this.deps.now() - started : undefined),
          () => undefined,
        )
        .then((ms) => this.settle(kind, race, source, ms));
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
