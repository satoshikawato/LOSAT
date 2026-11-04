// The compatibility outputs that LOSAT Web stores and exports byte for byte: the formats that every
// run writes, and the output of the argv's own -outfmt, 6 or 7 with fields (WA2-OP-03, Owner
// 2026-11-05; docs/web/abi_v2.md §5, §7).
export type OutputFormat = 0 | 6 | 7;
export const OUTPUT_FORMATS: readonly OutputFormat[] = Object.freeze([0, 6, 7]);

/** The output of the argv's own -outfmt (ABI stream 4). */
export const CUSTOM_OUTPUT = 4;

/** One stored output of a run: a format that every run writes, or CUSTOM_OUTPUT. */
export type RunOutput = OutputFormat | typeof CUSTOM_OUTPUT;

/** The flag whose value the engine reads as the CLI does (Additional options). */
export const OUTFMT_FLAG = '-outfmt';

/** NCBI's white space around an -outfmt value (`NStr::TruncateSpaces`). */
const NCBI_SPACE = /^[ \t\n\v\f\r]+|[ \t\n\v\f\r]+$/g;

/**
 * Where an argv gives -outfmt, as `-outfmt value` or `-outfmt=value` (the CLI takes it once), and
 * its value; undefined where it does not. The names of the inputs (file names) are values, not
 * options.
 */
export function outfmtOption(argv: readonly string[]): { readonly at: number; readonly words: 1 | 2; readonly value: string } | undefined {
  for (let i = 1; i < argv.length; i++) {
    const word = argv[i]!;
    if (word === '-query' || word === '-subject') {
      i++;
      continue;
    }
    if (word === OUTFMT_FLAG && i + 1 < argv.length) return { at: i, words: 2, value: argv[i + 1]! };
    if (word.startsWith(`${OUTFMT_FLAG}=`)) return { at: i, words: 1, value: word.slice(OUTFMT_FLAG.length + 1) };
  }
  return undefined;
}

/**
 * The -outfmt value of an argv whose run adds an output of its own (CUSTOM_OUTPUT): 6 or 7 with
 * words after the number, without NCBI's white space around it, as the adapter reads it
 * (`run::parse_run`). The rule reads the words, not the fields that they give: `6 std` adds an
 * output too (with the bytes of outfmt 6). Undefined for any other argv; the engine refuses the
 * values that are not formats, and the formats besides 0, 6 and 7.
 */
export function customOutfmt(argv: readonly string[]): string | undefined {
  const value = outfmtOption(argv)?.value.replace(NCBI_SPACE, '');
  const space = value?.indexOf(' ') ?? -1;
  return space >= 0 && /^\+?0*[67]$/.test(value!.slice(0, space)) ? value : undefined;
}

/** The outputs that a run of `argv` writes, in the order of the Outputs and "Download All". */
export function runOutputs(argv: readonly string[]): readonly RunOutput[] {
  return customOutfmt(argv) === undefined ? OUTPUT_FORMATS : [...OUTPUT_FORMATS, CUSTOM_OUTPUT];
}

/** The -outfmt value of one output of a run of `argv`: "6", or the argv's own ("6 qseqid sseqid"). */
export function outfmtOf(argv: readonly string[], output: RunOutput): string {
  return output === CUSTOM_OUTPUT ? (customOutfmt(argv) ?? '') : String(output);
}

/**
 * The file name of an exported output: the run's number, program and task (none for a program
 * without tasks), and the UTC time of its "Run LOSAT" in ISO 8601's basic format, so that the
 * runs of the four BLASTN tasks, and runs of the same number in other tabs, get other names
 * (`losat-run3-blastn-megablast-20261104T111530Z.outfmt6.txt`). The output of the argv's -outfmt
 * is `…outfmt6-custom.txt` (BLAST+'s "custom format" of its -outfmt help).
 */
export function outputFileName(
  run: {
    readonly number: number;
    readonly program: string;
    readonly task?: string | undefined;
    readonly queuedAt: number;
    readonly argv: readonly string[];
  },
  output: RunOutput,
): string {
  const time = new Date(run.queuedAt).toISOString().replace(/[-:]/g, '').replace(/\.\d+Z$/, 'Z');
  const task = run.task === undefined ? '' : `-${run.task}`;
  const format = output === CUSTOM_OUTPUT ? `${Number(outfmtOf(run.argv, output).split(' ', 1)[0])}-custom` : String(output);
  return `losat-run${run.number}-${run.program}${task}-${time}.outfmt${format}.txt`;
}
