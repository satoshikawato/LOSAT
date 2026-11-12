// Which sources of the data layer the search draft still uses (SURVEY #15). The Data worker
// keeps a source (its File: a pasted text, or the inflated bytes of a .gz, .bz2 or .zst, are held in
// memory) and every record table made of it until it is told to drop them. The holders:
// - the form: the sources that its cards show (both FASTA reader kinds of a source share its id);
// - an index that is still running, until its card shows the source;
// - a "Run LOSAT" being started, until the runs have their inputs.
// A run does not hold its sources: its snapshot keeps the bytes and the records of its inputs
// (coordinator.ts), its outputs are read by its id, and nothing reads its revisions again.
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
}
