// The compatibility output formats that LOSAT Web stores and exports byte for byte.
export type OutputFormat = 0 | 6 | 7;
export const OUTPUT_FORMATS: readonly OutputFormat[] = Object.freeze([0, 6, 7]);

/**
 * The file name of an exported output: the run's number, program and task (none for a program
 * without tasks), and the UTC time of its "Run LOSAT" in ISO 8601's basic format, so that the
 * runs of the four BLASTN tasks, and runs of the same number in other tabs, get other names
 * (`losat-run3-blastn-megablast-20261104T111530Z.outfmt6.txt`).
 */
export function outputFileName(
  run: { readonly number: number; readonly program: string; readonly task?: string | undefined; readonly queuedAt: number },
  format: OutputFormat,
): string {
  const time = new Date(run.queuedAt).toISOString().replace(/[-:]/g, '').replace(/\.\d+Z$/, 'Z');
  const task = run.task === undefined ? '' : `-${run.task}`;
  return `losat-run${run.number}-${run.program}${task}-${time}.outfmt${format}.txt`;
}
