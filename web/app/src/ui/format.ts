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
