// Search conditions are represented only as a LOSAT CLI argv (plan §5.3). The same argv
// is validated and run by the engine, shown as the CLI command, and stored in sessions.
import type { ProgramId } from './programs';

/** Flags that the application manages itself; users cannot set them as parameters. */
export const RESERVED_FLAGS: readonly string[] = Object.freeze([
  '-query',
  '-subject',
  '-out',
  '-outfmt',
  '-num_threads',
]);

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

/** Renders an argv as a copyable POSIX shell command. */
export function toShellCommand(argv: readonly string[]): string {
  const quote = (word: string) =>
    /^[A-Za-z0-9_@%+=:,./-]+$/.test(word) ? word : `'${word.replaceAll("'", "'\\''")}'`;
  return ['LOSAT', ...argv].map(quote).join(' ');
}
