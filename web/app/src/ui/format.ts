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

/** A short time span with its unit, for example "420 ms" or "2.6 s" (a step of a run, which `formatDuration` rounds down to the second). */
export function formatMillis(ms: number): string {
  return ms < 1000 ? `${Math.round(ms)} ms` : `${(ms / 1000).toFixed(1)} s`;
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

/** The first paragraph of an option's help (`describe`), without the references to NCBI's sources. */
export function optionHelp(raw: string | undefined): string {
  const paragraph = (raw ?? '').split(/\n\s*\n/)[0] ?? '';
  return paragraph.replace(/\s*(?:NCBI reference|Reference):.*$/s, '').trim();
}

/** A count and its noun, singular for one: "1 subject", "3 subjects", "1 HSP", "5 HSPs" (W4b screen review L12). */
export function formatCounted(count: number, noun: string): string {
  return `${formatCount(count)} ${count === 1 ? noun : `${noun}s`}`;
}

/**
 * Run numbers as words, in order: "run 2", "runs 1 and 3", "runs 1–4", "runs 1, 3 and 5–7". Three or more
 * runs in a row are written as a range, so that hundreds of runs stay short.
 */
export function formatRuns(numbers: readonly number[]): string {
  const sorted = [...new Set(numbers)].sort((a, b) => a - b);
  const parts: string[] = [];
  for (let first = 0; first < sorted.length; ) {
    let last = first;
    while (last + 1 < sorted.length && sorted[last + 1] === sorted[last]! + 1) last++;
    if (last - first >= 2) parts.push(`${sorted[first]}–${sorted[last]}`);
    else parts.push(...sorted.slice(first, last + 1).map(String));
    first = last + 1;
  }
  const list = parts.length < 2 ? (parts[0] ?? '') : `${parts.slice(0, -1).join(', ')} and ${parts[parts.length - 1]}`;
  return `${sorted.length === 1 ? 'run' : 'runs'} ${list}`;
}

/**
 * What the notice after a hidden tab says of the completed runs whose stored results could not be read back
 * (WA2-SC-05): the runs the browser removed them of, then those that failed to read otherwise, and that they can be
 * run again. Empty when there are none.
 */
export function lostResultsText(removed: readonly number[], unreadable: readonly number[]): string {
  const sentences = [
    ...(removed.length > 0 ? [`The browser removed the stored results of ${formatRuns(removed)}.`] : []),
    ...(unreadable.length > 0 ? [`The stored results of ${formatRuns(unreadable)} could not be read.`] : []),
  ];
  if (sentences.length === 0) return '';
  return `${sentences.join(' ')} Run ${removed.length + unreadable.length === 1 ? 'it' : 'them'} again.`;
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
