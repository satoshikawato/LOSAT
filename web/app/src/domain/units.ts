// A number and its unit, joined for display (line WQA's sweep, D9, 2026-11-19).

/** No-break space (U+00A0): the space between a number and its unit, which a line never breaks at. */
export const NBSP = '\u00A0';

/**
 * A number and the unit that follows it, "512 MB", "1,200 aa", "420 ms", joined by a no-break space so
 * that a narrow line never ends on the number and starts the next with the unit ("at most 512" / "MB",
 * D9). The one way the application composes them for display; BLAST output, command lines and the
 * engine's own messages are shown as they come and do not pass here.
 */
export function withUnit(value: number | string, unit: string): string {
  return `${value}${NBSP}${unit}`;
}
