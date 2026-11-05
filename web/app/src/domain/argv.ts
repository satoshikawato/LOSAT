// Search conditions are represented only as a LOSAT CLI argv (plan §5.3). The same argv
// is validated and run by the engine, shown as the CLI command, and stored in sessions.
import { CUSTOM_OUTPUT, OUTFMT_FLAG, outfmtOption, type RunOutput } from './output-format';
import type { ProgramId } from './programs';

/**
 * Flags that the application manages itself; users cannot set them. A run writes outfmt 0, 6
 * and 7, and -outfmt in "Additional options" may add one (domain/output-format.ts).
 */
export const RESERVED_FLAGS: readonly string[] = Object.freeze(['-query', '-subject', '-out', '-num_threads']);

/** The words with which the CLI prints its help instead of searching (BLAST+'s `-h` and `-help`). */
const HELP_WORDS: readonly string[] = Object.freeze(['-h', '-help', '--help']);

/** The CLI's database options: `-db` and the filters of a database search. */
export const DATABASE_FLAGS: readonly string[] = Object.freeze([
  '-db',
  '-seqidlist',
  '-negative_seqidlist',
  '-taxids',
  '-taxidlist',
  '-negative_taxids',
  '-negative_taxidlist',
  '-no_taxid_expansion',
]);

/** The modes whose CLI has -mt_mode, which needs -num_threads (TBLASTX's has none). */
const MT_MODE_PROGRAMS: readonly ProgramId[] = Object.freeze(['blastn', 'blastp', 'blastx', 'tblastn']);

/**
 * Why a word of "Additional options" is refused before the engine reads it: options that the
 * application manages, and the CLI's options that a search in the browser cannot have, whose
 * CLI message would name what the user cannot change here (`-subject`, `-num_threads`) or print
 * the CLI's help (WA-FM-04, WA-FM-05); undefined for the other words.
 */
function refusedWord(program: ProgramId, word: string): string | undefined {
  // -flag=value is the command line's form of -flag value too.
  const flag = word.split('=', 1)[0]!;
  if (flag === '-num_threads') return '-num_threads is set in the form: Threads';
  if (RESERVED_FLAGS.includes(flag) || (flag === '-mt_mode' && MT_MODE_PROGRAMS.includes(program))) return `${flag} is managed by LOSAT`;
  if (HELP_WORDS.includes(word)) return `${word} is not available in the browser`;
  if (flag === '-db') return '-db (a database search) is not supported in the browser yet';
  if (DATABASE_FLAGS.includes(flag)) return `${flag} (a database search filter) is not supported in the browser yet`;
  return undefined;
}

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
    const refused = refusedWord(input.program, word);
    if (refused !== undefined) throw new Error(refused);
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
 * Renders the CLI command that reproduces one output of a run, as a copyable POSIX shell
 * command. `-outfmt` is always written: several programs reject their CLI default format, so a
 * command without it may not run. The command of a format that every run writes has the
 * format's -outfmt in place of the argv's; that of the argv's own -outfmt is the argv.
 */
export function toShellCommand(argv: readonly string[], output: RunOutput): string {
  return shellWords(commandWords(argv, output));
}

/** The words of `toShellCommand`'s command. */
export function commandWords(argv: readonly string[], output: RunOutput): readonly string[] {
  if (output === CUSTOM_OUTPUT) return ['losat', ...argv];
  const option = outfmtOption(argv);
  const words = option === undefined ? argv : [...argv.slice(0, option.at), ...argv.slice(option.at + option.words)];
  return ['losat', ...words, '-outfmt', String(output)];
}

/** A word as a POSIX shell reads it back (`splitOptions`): quoted only where it needs it. */
export function shellWord(word: string): string {
  return /^[A-Za-z0-9_@%+=:,./-]+$/.test(word) ? word : `'${word.replaceAll("'", "'\\''")}'`;
}

/** Words as a POSIX shell reads them back (`splitOptions`): quoted only where they need it. */
export function shellWords(words: readonly string[]): string {
  return words.map(shellWord).join(' ');
}

/**
 * The words of an argv or a command line as `shellWords` writes them, and which of them name an
 * option: a word that starts with '-' where an option starts, not where the option before it
 * takes its value (`-penalty -3`). `takesValue` (`describe`'s `takes_value`) tells the options
 * that take a value, besides those that the application manages; without it, every word that
 * starts with '-' names an option.
 */
export function optionNames(
  words: readonly string[],
  takesValue: ((flag: string) => boolean) | undefined,
): ReadonlyArray<{ readonly text: string; readonly name: boolean }> {
  let value = false;
  return words.map((word) => {
    const name = !value && word.startsWith('-');
    // -flag=value carries its value.
    value = name && takesValue !== undefined && !word.includes('=') && (RESERVED_FLAGS.includes(word) || word === OUTFMT_FLAG || takesValue(word));
    return { text: shellWord(word), name };
  });
}

/** The value of an option in an argv (the word after its last occurrence), or undefined where it is not given. */
export function optionValue(argv: readonly string[], flag: string): string | undefined {
  const at = argv.lastIndexOf(flag);
  return at < 0 || at + 1 >= argv.length ? undefined : argv[at + 1];
}
