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
  /** Words of other options, as on a command line ("Additional options", `splitOptions`), after the parameters. */
  readonly extra?: readonly string[];
}

/**
 * Builds `[program, -query, name, -subject, name, ...parameters]`.
 * Engine defaults are never written out, so the CLI command stays minimal.
 */
export function buildArgv(input: ArgvInput): readonly string[] {
  const argv = [input.program, '-query', input.queryName, '-subject', input.subjectName];
  for (const [flag, value] of input.parameters) {
    if (!flag.startsWith('-')) throw new Error(`parameter must start with "-": ${flag}`);
    if (RESERVED_FLAGS.includes(flag)) throw new Error(`${flag} is managed by LOSAT`);
    argv.push(flag);
    if (value !== true) argv.push(value);
  }
  for (const word of input.extra ?? []) {
    // -flag=value is the command line's form of -flag value too.
    const flag = word.split('=', 1)[0]!;
    if (RESERVED_FLAGS.includes(flag)) throw new Error(`${flag} is managed by LOSAT`);
    // The engine reads the argv as one NUL-separated text whose last empty word it drops, and a
    // run adds -num_threads after the words, which a last -- would leave out of the options: a
    // word that the validation and the run would read differently is refused.
    const unread = word === '' ? "''" : word === '--' ? '--' : word.includes('\0') ? 'a NUL character' : undefined;
    if (unread !== undefined) throw new Error(`Additional options: ${unread} is not accepted`);
    argv.push(word);
  }
  return Object.freeze(argv);
}

/**
 * The words of options typed as on a POSIX shell's command line ("Additional options"):
 * separated by white space, with '…' and "…" quoting and \ escaping, so that the words of
 * `toShellCommand` read back as they were written. Throws on an unclosed quote.
 */
export function splitOptions(text: string): readonly string[] {
  const words: string[] = [];
  let word: string | undefined;
  for (let i = 0; i < text.length; i++) {
    const c = text[i]!;
    if (/\s/.test(c)) {
      if (word !== undefined) words.push(word);
      word = undefined;
    } else if (c === "'") {
      const end = text.indexOf("'", i + 1);
      if (end < 0) throw new Error("Additional options: a quote (') is not closed");
      word = (word ?? '') + text.slice(i + 1, end);
      i = end;
    } else if (c === '"') {
      let part = '';
      for (i++; i < text.length && text[i] !== '"'; i++) {
        // In double quotes, \ escapes only " \ $ ` (POSIX).
        if (text[i] === '\\' && i + 1 < text.length && '"\\$`'.includes(text[i + 1]!)) i++;
        part += text[i];
      }
      if (i >= text.length) throw new Error('Additional options: a quote (") is not closed');
      word = (word ?? '') + part;
    } else if (c === '\\' && i + 1 < text.length) {
      word = (word ?? '') + text[++i];
    } else {
      word = (word ?? '') + c;
    }
  }
  if (word !== undefined) words.push(word);
  return Object.freeze(words);
}

/**
 * Renders the CLI command that reproduces one output format of a run, as a copyable
 * POSIX shell command. `-outfmt` is always written: several programs reject their
 * CLI default format, so a command without it may not run.
 */
export function toShellCommand(argv: readonly string[], outfmt: OutputFormat): string {
  return shellWords(['losat', ...argv, '-outfmt', String(outfmt)]);
}

/** Words as a POSIX shell reads them back (`splitOptions`): quoted only where they need it. */
export function shellWords(words: readonly string[]): string {
  const quote = (word: string) =>
    /^[A-Za-z0-9_@%+=:,./-]+$/.test(word) ? word : `'${word.replaceAll("'", "'\\''")}'`;
  return words.map(quote).join(' ');
}

/** The value of an option in an argv (the word after its last occurrence), or undefined where it is not given. */
export function optionValue(argv: readonly string[], flag: string): string | undefined {
  const at = argv.lastIndexOf(flag);
  return at < 0 || at + 1 >= argv.length ? undefined : argv[at + 1];
}
