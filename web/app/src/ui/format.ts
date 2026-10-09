// Display helpers of the UI. They format application values only, never BLAST values.
const UNITS = ['bytes', 'kB', 'MB', 'GB', 'TB'] as const;

/** A byte count in SI units, for example "1.5 MB". */
export function formatBytes(bytes: number): string {
  let value = bytes;
  let unit = 0;
  while (value >= 1000 && unit < UNITS.length - 1) {
    value /= 1000;
    unit++;
  }
  const digits = unit === 0 ? 0 : 1;
  return `${value.toFixed(digits)} ${UNITS[unit]}`;
}

/** A count with thousands separators, for example "12,345". */
export function formatCount(count: number): string {
  return count.toLocaleString('en-US');
}

/** A duration as m:ss or h:mm:ss, for example "1:05". */
export function formatDuration(ms: number): string {
  const total = Math.max(0, Math.floor(ms / 1000));
  const hours = Math.floor(total / 3600);
  const minutes = Math.floor((total % 3600) / 60);
  const seconds = String(total % 60).padStart(2, '0');
  return hours > 0 ? `${hours}:${String(minutes).padStart(2, '0')}:${seconds}` : `${minutes}:${seconds}`;
}

/**
 * A time in ISO 8601 form in local time, date and time of day with a space between them, for
 * example "2026-10-09 21:49:30" (S13 screen review L5).
 */
export function formatDateTime(ms: number): string {
  const time = new Date(ms);
  const two = (n: number) => String(n).padStart(2, '0');
  const date = `${String(time.getFullYear()).padStart(4, '0')}-${two(time.getMonth() + 1)}-${two(time.getDate())}`;
  return `${date} ${two(time.getHours())}:${two(time.getMinutes())}:${two(time.getSeconds())}`;
}

/** A count and its noun, singular for one: "1 subject", "3 subjects", "1 HSP", "5 HSPs" (W4b screen review L12). */
export function formatCounted(count: number, noun: string): string {
  return `${formatCount(count)} ${count === 1 ? noun : `${noun}s`}`;
}

/** The frames of an HSP of a translated search: the translated sequences' only (results.ts hspEntry). */
export interface Frames {
  readonly queryFrame?: number;
  readonly subjectFrame?: number;
}

/** A frame with its sign, for example "+2" or "-1". */
export const signedFrame = (frame: number): string => (frame > 0 ? `+${frame}` : String(frame));

/**
 * The name of an HSP's frames, where only the translated sequence has one: "Subject frame"
 * (TBLASTN), "Query frame" (BLASTX), "Frames (q/s)" where both are translated (TBLASTX).
 * TBLASTN's "–/+2" read as a minus strand beside TBLASTX's "-2/+2" (W4b screen review L9).
 * Undefined for an HSP without frames.
 */
export function framesLabel(hsp: Frames): string | undefined {
  if (hsp.queryFrame !== undefined && hsp.subjectFrame !== undefined) return 'Frames (q/s)';
  if (hsp.subjectFrame !== undefined) return 'Subject frame';
  return hsp.queryFrame === undefined ? undefined : 'Query frame';
}

/** The frames under that name: "+2", or the query's and the subject's joined by `separator` ("-2/+2"). */
export function framesText(hsp: Frames, separator = '/'): string {
  return [hsp.queryFrame, hsp.subjectFrame]
    .filter((frame): frame is number => frame !== undefined)
    .map(signedFrame)
    .join(separator);
}

/** The frames in a sentence: "subject frame +2", "query frame -1", "frames -2 / +2"; empty without frames. */
export function framesPhrase(hsp: Frames): string {
  const label = framesLabel(hsp);
  if (label === undefined) return '';
  return label === 'Frames (q/s)' ? `frames ${framesText(hsp, ' / ')}` : `${label.toLowerCase()} ${framesText(hsp)}`;
}
