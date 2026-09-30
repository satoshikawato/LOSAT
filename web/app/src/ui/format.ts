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
