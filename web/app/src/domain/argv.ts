// Search conditions are represented only as a LOSAT CLI argv (plan §5.3). The same argv
// is validated and run by the engine, shown as the CLI command, and stored in sessions.
import type { OutputFormat } from './output-format';
import type { ProgramId } from './programs';

/** Flags that the application manages itself; users cannot set them as parameters. */
export const RESERVED_FLAGS: readonly string[] = Object.freeze([
  '-query',
  '-subject',
  '-out',
  '-outfmt',
  '-num_threads',
]);

/**
 * The -query / -subject names of inputs that are not one file (plan §5.3): pasted text, and
 * several inputs searched together.
 */
export const PASTED_NAMES = Object.freeze({ query: 'query.fa', subject: 'subject.fa' } as const);
export const COMBINED_NAMES = Object.freeze({ query: 'combined_query.fa', subject: 'combined_subject.fa' } as const);

export interface ArgvInput {
  readonly program: ProgramId;
  readonly queryName: string;
  readonly subjectName: string;
  /** Parameters the user set explicitly, in display order. `true` means a bare flag. */
  readonly parameters: ReadonlyArray<readonly [flag: string, value: string | true]>;
}

/**
 * Builds `[program, -query, name, -subject, name, ...parameters]`.
 * Engine defaults are never written out, so the CLI command stays minimal.
 */
export function buildArgv(input: ArgvInput): readonly string[] {
  const argv = [input.program, '-query', input.queryName, '-subject', input.subjectName];
  for (const [flag, value] of input.parameters) {
    if (!flag.startsWith('-')) throw new Error(`parameter must start with "-": ${flag}`);
    if (RESERVED_FLAGS.includes(flag)) throw new Error(`${flag} is managed by LOSAT Web`);
    argv.push(flag);
    if (value !== true) argv.push(value);
  }
  return Object.freeze(argv);
}

/**
 * Renders the CLI command that reproduces one output format of a run, as a copyable
 * POSIX shell command. `-outfmt` is always written: several programs reject their
 * CLI default format, so a command without it may not run.
 */
export function toShellCommand(argv: readonly string[], outfmt: OutputFormat): string {
  const quote = (word: string) =>
    /^[A-Za-z0-9_@%+=:,./-]+$/.test(word) ? word : `'${word.replaceAll("'", "'\\''")}'`;
  return ['LOSAT', ...argv, '-outfmt', String(outfmt)].map(quote).join(' ');
}

/** The value of an option in an argv (the word after its last occurrence), or undefined where it is not given. */
export function optionValue(argv: readonly string[], flag: string): string | undefined {
  const at = argv.lastIndexOf(flag);
  return at < 0 || at + 1 >= argv.length ? undefined : argv[at + 1];
}
