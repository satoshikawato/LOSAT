// Which sources of the data layer are still in use (SURVEY #15). The Data worker keeps a
// source (its File: a pasted text, or the inflated bytes of a .gz or .zst, are held in
// memory) and every record table made of it until it is told to drop them. The holders:
// - the search form: the sources that its cards show (both FASTA reader kinds of a source
//   share its id), and an index that is still running;
// - a "Run LOSAT" being started, from the draft's revisions until its runs exist;
// - each run made of the source (`run:<runId>`): a run keeps what it was made of, also for
//   later uses of its inputs (design D08). A tab keeps every run, so these holds last as long as
//   the tab; an operation that removes a run releases its hold.
// When the last holder lets a source go, the data layer drops the source and its tables.
// Sources never reach OPFS (the run outputs and the saved Example sets do, and have their own
// lifetimes), so dropping one frees memory only.

export class SourceRefs {
  /** The holders of each source that someone holds. */
  private readonly holders = new Map<string, Set<string>>();
  /** The sources of each holder. */
  private readonly held = new Map<string, ReadonlySet<string>>();

  constructor(private readonly drop: (sourceId: string) => void) {}

  /** The holder now holds exactly `sourceIds`: the ones it no longer names are let go. */
  set(holder: string, sourceIds: Iterable<string>): void {
    const next = new Set(sourceIds);
    const previous = this.held.get(holder) ?? new Set<string>();
    if (next.size === 0) this.held.delete(holder);
    else this.held.set(holder, next);
    for (const sourceId of next) {
      if (previous.has(sourceId)) continue;
      let holders = this.holders.get(sourceId);
      if (holders === undefined) this.holders.set(sourceId, (holders = new Set()));
      holders.add(holder);
    }
    for (const sourceId of previous) {
      if (next.has(sourceId)) continue;
      const holders = this.holders.get(sourceId);
      holders?.delete(holder);
      if (holders === undefined || holders.size === 0) {
        this.holders.delete(sourceId);
        this.drop(sourceId);
      }
    }
  }

  release(holder: string): void {
    this.set(holder, []);
  }

  /** Whether anyone holds the source. */
  holds(sourceId: string): boolean {
    return this.holders.has(sourceId);
  }
}
