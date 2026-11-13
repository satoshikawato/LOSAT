// WindowMasker counts of a taxonomy ID (line WW2, DESIGN_ww2 §4): BLASTN's
// `-window_masker_taxid <n>` reads the unit counts of NCBI's windowmasker_files for <n>. LOSAT Web
// has no directory of them: the user downloads the file for the ID, only when asked (the line under
// the taxonomy ID box, application/counts.ts), the browser keeps it, and the search reads it under
// the name that the engine looks for (docs/web/abi_v2.md §4).

/** Flag of BLASTN's taxonomy ID option, the companion of `-window_masker_db` in the form's row. */
export const WINDOW_MASKER_TAXID_FLAG = '-window_masker_taxid';

/** The site that the counts come from, as the line names it. */
export const COUNTS_HOST = 'ftp.ncbi.nlm.nih.gov';

/** The name of the file in NCBI's directory of the ID, and of the saved file. */
export const COUNTS_FILE = 'wmasker.obinary';

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

/** What the run records of the counts that it read for a taxonomy ID (RunRecord.windowMasker.origin). */
export interface CountsOrigin {
  readonly taxid: number;
  readonly url: string;
  /** The Last-Modified of the file when it was downloaded, as the site sent it; '' when it sent none. */
  readonly lastModified: string;
}
