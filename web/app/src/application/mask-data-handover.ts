// The hand-over of a masking-information output to makeblastdb (WW4c): a run of the Masking tab, or
// WindowMasker's Mask stage, whose `-outfmt` is a `maskinfo_*` format marks the output that holds it
// (`ToolOutputEntry.isMaskData`), and the Results tab offers "Use as -mask_data in makeblastdb" for it.
// The makeblastdb form is the target: this one port is all that connects the two, so the button shows
// only while a target is there (`available`).

/** A masking-information output of a run of this session, as the Results tab names it. */
export interface MaskDataEntry {
  readonly runId: string;
  readonly runNumber: number;
  /** The key under which the Data worker keeps the output (`'stdout'`, `'file:0'`, …). */
  readonly key: string;
  /** The output's name: the `-out` file's, or `'stdout'`. */
  readonly name: string;
  readonly size: number;
}

/** What takes a masking-information output: the makeblastdb form, which adds it to its `-mask_data`. */
export interface MaskDataTarget {
  useMaskData(entry: MaskDataEntry): void;
}

export class MaskDataHandover {
  constructor(private readonly target?: MaskDataTarget) {}

  /** Whether an output can be handed over (a makeblastdb form is there to take it). */
  available(): boolean {
    return this.target !== undefined;
  }

  /** Hands the output over to the makeblastdb form; does nothing while none is there. */
  use(entry: MaskDataEntry): void {
    this.target?.useMaskData(entry);
  }
}
