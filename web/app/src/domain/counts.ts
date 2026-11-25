// WindowMasker counts of a taxonomy ID (line WW2, DESIGN_ww2 §4): BLASTN's
// `-window_masker_taxid <n>` reads the unit counts of NCBI's windowmasker_files for <n>. LOSAT Web
// has no directory of them: the user downloads the file for the ID, only when asked (the line under
// the taxonomy ID box, application/counts.ts), the browser keeps it, and the search reads it under
// the name that the engine looks for (docs/web/abi_v2.md §4).
import { dateOf } from './download-sentences';

/** The site that the counts come from, as the line names it. */
export const COUNTS_HOST = 'ftp.ncbi.nlm.nih.gov';

/** The name of the file in NCBI's directory of the ID, and of the saved file. */
export const COUNTS_FILE = 'wmasker.obinary';

/** Web Lock name prefix for saved WindowMasker counts; shared by every release. */
export const COUNTS_LOCK_PREFIX = 'losat-web:windowmasker:';

/** The Web Lock name of a taxonomy ID's saved counts. */
export function countsLockName(taxid: number): string {
  return COUNTS_LOCK_PREFIX + taxid;
}

/** The URL of the counts of a taxonomy ID. The request carries nothing but this URL. */
export function countsUrl(taxid: number): string {
  return `https://${COUNTS_HOST}/blast/windowmasker_files/${taxid}/${COUNTS_FILE}`;
}

/** The name under which a search reads the counts of a taxonomy ID: the one the engine looks for. */
export function countsName(taxid: number): string {
  return `windowmasker/${taxid}/${COUNTS_FILE}`;
}

const INT32_MIN = -2147483648;
const INT32_MAX = 2147483647;

/**
 * The taxonomy ID that the box holds: a decimal whole number in the range of a 32-bit signed
 * integer, other than 0 (0 is NCBI's "no WindowMasker"); undefined for anything else.
 */
export function parseTaxid(text: string): number | undefined {
  const trimmed = text.trim();
  if (!/^-?\d{1,10}$/.test(trimmed)) return undefined;
  const value = Number(trimmed);
  return value === 0 || value < INT32_MIN || value > INT32_MAX ? undefined : value;
}

/** One input file recorded in the metadata of a counts file made in this browser. */
export interface CountsMadeInput {
  readonly name: string;
  readonly size: number;
}

/** Provenance stored under `made` in meta.json for counts made in this browser's WindowMasker tab. */
export interface CountsMadeMeta {
  readonly run: number;
  readonly argv: readonly string[];
  readonly inputs: readonly CountsMadeInput[];
  readonly engineBuild: string;
  /** When the counts were made, in milliseconds since the epoch. */
  readonly madeAt: number;
}

/** Whether a saved copy of counts came from NCBI's site or from this browser's WindowMasker tab. */
export type SavedCountsOrigin = 'downloaded' | 'made';

/** What a saved file of counts records about itself (meta.json). */
export interface CountsCopy {
  readonly taxid: number;
  readonly url: string;
  /** The Last-Modified of the first answer, as the site sent it ('' when it sent none or for a made copy). */
  readonly lastModified: string;
  /** The size of the whole file. */
  readonly size: number;
  /** Bytes kept from the start of the file; Continue asks for the file from this byte. */
  readonly savedBytes: number;
  /** The whole file is kept. */
  readonly complete: boolean;
  /** When the last part was saved, in milliseconds since the epoch. */
  readonly savedAt: number;
  /** Present when the counts were made in this browser's WindowMasker tab; absent for an NCBI download. */
  readonly made?: CountsMadeMeta;
}

/** One complete saved counts file listed for the WindowMasker tab's pickers (`listCounts`). */
export interface SavedCountsEntry {
  readonly taxid: number;
  /** The file name under windowmasker/<taxid>/ (wmasker.obinary). */
  readonly name: string;
  readonly size: number;
  readonly origin: SavedCountsOrigin;
  /** ISO date (yyyy-mm-dd) of the copy: Last-Modified for a download, madeAt for a made copy, or savedAt. */
  readonly date: string;
  readonly copy: CountsCopy;
}

/** Builds the `listCounts` summary of one complete saved copy. */
export function countsEntryOf(copy: CountsCopy): SavedCountsEntry {
  const origin: SavedCountsOrigin = copy.made !== undefined ? 'made' : 'downloaded';
  const date =
    copy.made !== undefined
      ? isoDate(copy.made.madeAt)
      : (dateOf(copy.lastModified) ?? isoDate(copy.savedAt));
  return { taxid: copy.taxid, name: COUNTS_FILE, size: copy.size, origin, date, copy };
}

function isoDate(ms: number): string {
  return new Date(ms).toISOString().slice(0, 10);
}

/** What the run records of the counts that it read for a taxonomy ID (RunRecord.windowMasker.origin). */
export type CountsOrigin =
  | {
      readonly kind: 'downloaded';
      readonly taxid: number;
      readonly url: string;
      /** The Last-Modified of the file when it was downloaded, as the site sent it; '' when it sent none. */
      readonly lastModified: string;
    }
  | {
      readonly kind: 'made';
      readonly taxid: number;
      readonly made: CountsMadeMeta;
    };
