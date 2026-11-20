// The WindowMasker unit counts file of a run (`-window_masker_db`, or `-window_masker_taxid`): the
// file that the engine reads under the name that the argv gives, or that the engine looks for
// (docs/web/abi_v2.md §4).

/** The counts file that a run gives the engine. */
export interface CountsFile {
  /** The value of `-window_masker_db` in the argv; for `-window_masker_taxid`, the name that the engine looks for (windowmasker/<taxid>/wmasker.obinary). */
  readonly name: string;
  readonly file: Blob;
  readonly size: number;
  readonly lastModified: number;
}

/** What a run did with its counts file. */
export interface CountsUse {
  readonly name: string;
  readonly size: number;
  /** Lower-case hex SHA-256 of the bytes that the engine read. */
  readonly sha256: string;
  /** The engine used the counts that it held from an earlier run. */
  readonly retained: boolean;
  /** Milliseconds spent giving the file to the engine; 0 when it was retained. */
  readonly loadMs: number;
}
