import { CUSTOM_OUTPUT, OUTPUT_FORMATS, type RunOutput } from 'losat/domain';
// The compatibility outputs that LOSAT Web stores and exports byte for byte: the formats that every
// run writes, and the output of the argv's own -outfmt: any other format that the CLI writes, and
// 6 or 7 with fields (WA2-OP-03 and W11, Owner 2026-11-05; docs/web/abi_v2.md §5, §7).

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
 * The -outfmt value of an argv whose run adds an output of its own (CUSTOM_OUTPUT), without
 * NCBI's white space around it, as the adapter reads it (`run::parse_run`): a format number
 * besides 0, 6 and 7 (with or without words after it), or 6 or 7 with words after the number. The
 * rule reads the words, not the fields that they give: `6 std` adds an output too (with the bytes
 * of outfmt 6). Undefined for any other argv (0 with words is outfmt 0, as the CLI reads it); the
 * engine refuses the values that are not formats that the program writes.
 */
export function customOutfmt(argv: readonly string[]): string | undefined {
  const value = outfmtOption(argv)?.value.replace(NCBI_SPACE, '');
  if (value === undefined) return undefined;
  const space = value.indexOf(' ');
  const number = formatNumber(value);
  if (number === undefined || number === 0) return undefined;
  return space >= 0 || (number !== 6 && number !== 7) ? value : undefined;
}

/** The format number of an -outfmt value (its first word, as NCBI reads it), or undefined where it is not a number. */
export function formatNumber(value: string): number | undefined {
  const head = value.replace(NCBI_SPACE, '').split(' ', 1)[0]!;
  return /^\+?\d+$/.test(head) ? Number(head) : undefined;
}

/** The formats that read the words after the number (BLAST+'s -outfmt help: "Options 6, 7, 10, 17 and 20 can be additionally configured"); the others ignore them. */
const CONFIGURABLE: ReadonlySet<number> = new Set([6, 7, 10, 17, 20]);

/** Whether an -outfmt value configures its format: words after a number of CONFIGURABLE (fields, `delim=`, SAM's tags). */
function configures(value: string): boolean {
  const number = formatNumber(value);
  return number !== undefined && CONFIGURABLE.has(number) && value.replace(NCBI_SPACE, '').includes(' ');
}

/** BLAST+'s names of the formats (its -outfmt help), for "Download All". */
export const FORMAT_NAMES: Readonly<Record<number, string>> = Object.freeze({
  0: 'Pairwise',
  5: 'BLAST XML',
  6: 'Tabular',
  7: 'Tabular with comment lines',
  8: 'Seqalign (Text ASN.1)',
  9: 'Seqalign (Binary ASN.1)',
  10: 'Comma-separated values',
  11: 'BLAST archive (ASN.1)',
  12: 'Seqalign (JSON)',
  13: 'Multiple-file BLAST JSON',
  14: 'Multiple-file BLAST XML2',
  15: 'Single-file BLAST JSON',
  16: 'Single-file BLAST XML2',
  17: 'Sequence Alignment/Map (SAM)',
  20: 'Comma-separated values with header lines',
});

/** The file extension and media type of each format's file; text for the rest. */
const FILE_TYPES: Readonly<Record<number, { readonly extension: string; readonly type: string }>> = Object.freeze({
  5: { extension: 'xml', type: 'application/xml' },
  8: { extension: 'asn', type: 'text/plain' },
  9: { extension: 'asnb', type: 'application/octet-stream' },
  10: { extension: 'csv', type: 'text/csv' },
  11: { extension: 'asn', type: 'text/plain' },
  12: { extension: 'json', type: 'application/json' },
  // The master file and the side files of 13 and 14, as one ZIP archive (the adapter's).
  13: { extension: 'zip', type: 'application/zip' },
  14: { extension: 'zip', type: 'application/zip' },
  15: { extension: 'json', type: 'application/json' },
  16: { extension: 'xml', type: 'application/xml' },
  17: { extension: 'sam', type: 'text/plain' },
  20: { extension: 'csv', type: 'text/csv' },
});
const TEXT_FILE = Object.freeze({ extension: 'txt', type: 'text/plain' });

/** The format number of one output of a run of `argv`. */
export function outputNumber(argv: readonly string[], output: RunOutput): number {
  return output === CUSTOM_OUTPUT ? (formatNumber(customOutfmt(argv) ?? '') ?? 0) : output;
}

/**
 * The `-out` of the files of outfmt 13 and 14 (the adapter's `run::SIDE_FILE_BASE`): the archive holds
 * the files that the CLI writes with `-out results` (the side files `results_<n>.json|xml`, the master
 * file as `results.json|xml`), so the command of that output names it.
 */
export const FILES_OUT = 'results';

/** Whether an output of a run of `argv` is the archive of the files of outfmt 13 or 14. */
export function isFilesOutput(argv: readonly string[], output: RunOutput): boolean {
  const number = outputNumber(argv, output);
  return output === CUSTOM_OUTPUT && (number === 13 || number === 14);
}

/** The formats whose output is binary: outfmt 9, and the ZIP archive of the files of 13 and 14. */
const BINARY: ReadonlySet<number> = new Set([9, 13, 14]);

/** Whether an output of a run of `argv` is binary (BINARY), which the Outputs do not show as text. */
export function isBinaryOutput(argv: readonly string[], output: RunOutput): boolean {
  return BINARY.has(outputNumber(argv, output));
}

/** The media type of an output's file. */
export function outputMediaType(argv: readonly string[], output: RunOutput): string {
  return (FILE_TYPES[outputNumber(argv, output)] ?? TEXT_FILE).type;
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
 * has the extension of its format (`…outfmt5.xml`, `…outfmt10.csv`, `…outfmt9.asnb`, `…outfmt13.zip`), and
 * `-custom` after the number where words configure it (`…outfmt6-custom.txt`, BLAST+'s "custom
 * format" of its -outfmt help; outfmt 5 ignores words).
 */
export function outputFileName(run: ExportedRun, output: RunOutput): string {
  const number = outputNumber(run.argv, output);
  const format = output === CUSTOM_OUTPUT && configures(outfmtOf(run.argv, output)) ? `${number}-custom` : String(number);
  const { extension } = FILE_TYPES[number] ?? TEXT_FILE;
  return `${exportPrefix(run)}.outfmt${format}.${extension}`;
}

/** What names the files of a run: its number, program, task and the time of its "Run LOSAT". */
export interface ExportedRun {
  readonly number: number;
  readonly program: string;
  readonly task?: string | undefined;
  readonly queuedAt: number;
  readonly argv: readonly string[];
}

/** The start of the name of every file of a run that the page saves (`losat-run3-blastn-megablast-20261104T111530Z`). */
function exportPrefix(run: ExportedRun): string {
  const time = new Date(run.queuedAt).toISOString().replace(/[-:]/g, '').replace(/\.\d+Z$/, 'Z');
  const task = run.task === undefined ? '' : `-${run.task}`;
  return `losat-run${run.number}-${run.program}${task}-${time}`;
}

/**
 * The file name of a file that a run made besides its outputs (the file of `-logfile NAME`; line RG
 * E1w-3b): the run's start of the names of its outputs, `logfile`, then the file's own name
 * (`losat-run3-tblastx-20261104T111530Z.logfile.run.log`), so that the files of two runs get other
 * names, and no file of a run gets the name of one of its outputs (`-logfile outfmt6.txt`).
 */
export function runFileName(run: ExportedRun, name: string): string {
  return `${exportPrefix(run)}.logfile.${name}`;
}
