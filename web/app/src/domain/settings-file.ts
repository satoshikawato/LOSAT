// Settings files (S15 instructions, item 3; REQ-15): the search conditions of the form or of a
// run - the program, the argv's options (the words after `-subject <name>`, the regions
// included) and the threads - as a small JSON document of LOSAT Web. Never the inputs, their
// names, the Job Title, results or notes. Reading is strict and bounded: a file that is not
// exactly this document is refused with what is wrong and where, so that a file never fills the
// form with something other than it says.
import { RESERVED_FLAGS } from './argv';
import { PROGRAMS, type ProgramId } from './programs';

export const SETTINGS_FORMAT = 'LOSAT Web search settings';
/** The schema this build writes and reads; a newer one is refused. */
export const SETTINGS_SCHEMA = 1;
export const SETTINGS_NOTE =
  'A LOSAT Web application format, not an NCBI BLAST+ file: the program, the command-line options ' +
  '(without the input names) and the threads of a search. It holds no sequences, names, results or notes.';
export const SETTINGS_MIME = 'application/json';
export const SETTINGS_MAX_BYTES = 1024 * 1024;
export const SETTINGS_MAX_OPTIONS = 1000;
export const SETTINGS_MAX_OPTION_CHARS = 10_000;
const MAX_NOTE_CHARS = 10_000;
/** The most threads that the search form offers (on a processor with that many or more). */
export const MAX_THREADS = 16;

export interface SearchSettings {
  readonly program: ProgramId;
  /** The argv's words after `-subject <name>`: the options with their values, and the regions. */
  readonly options: readonly string[];
  readonly threads: number | 'auto';
}

export type SettingsRead =
  | { readonly ok: true; readonly settings: SearchSettings }
  | { readonly ok: false; readonly message: string };

/** What a settings file or a run's settings left out of the form (S15 item 3: "Not applied: …"). */
export interface AppliedSettings {
  readonly notApplied: readonly string[];
  /** The words of the file's options that the form took (a left-out option's words are not counted). */
  readonly words: number;
}

const FIELDS = ['format', 'schema', 'note', 'program', 'options', 'threads'];

/** The threads that the search form offers on a processor with `hardware` logical processors. */
export function threadLimit(hardware: number | undefined): number {
  const count = hardware !== undefined && Number.isInteger(hardware) && hardware > 0 ? hardware : 4;
  return Math.min(MAX_THREADS, count);
}

export function settingsFileName(program: ProgramId): string {
  return `losat-settings-${program}.json`;
}

/** The settings of a run: its program, its argv after the inputs, and the threads it asked for. */
export function settingsOfRun(snapshot: {
  readonly program: ProgramId;
  readonly argv: readonly string[];
  readonly requestedThreads: number | 'auto';
}): SearchSettings {
  return { program: snapshot.program, options: snapshot.argv.slice(5), threads: snapshot.requestedThreads };
}

/**
 * The settings file's text in parts, in order (the Writer contract: a file is written in
 * blocks, never built whole): JSON as `JSON.stringify(document, null, 2)` writes it, one
 * option word per part, and a final line end.
 */
export function* settingsText(settings: SearchSettings): Generator<string> {
  yield '{\n';
  yield `  "format": ${JSON.stringify(SETTINGS_FORMAT)},\n`;
  yield `  "schema": ${SETTINGS_SCHEMA},\n`;
  yield `  "note": ${JSON.stringify(SETTINGS_NOTE)},\n`;
  yield `  "program": ${JSON.stringify(settings.program)},\n`;
  if (settings.options.length === 0) {
    yield '  "options": [],\n';
  } else {
    yield '  "options": [';
    for (const [i, word] of settings.options.entries()) yield `${i === 0 ? '' : ','}\n    ${JSON.stringify(word)}`;
    yield '\n  ],\n';
  }
  yield `  "threads": ${JSON.stringify(settings.threads)}\n}\n`;
}

/** The whole settings text (tests and small callers; files are written with `settingsText`). */
export function serializeSettings(settings: SearchSettings): string {
  return [...settingsText(settings)].join('');
}

/** Reads a settings file's bytes; every refusal says what is wrong and where. */
export function parseSettings(bytes: Uint8Array): SettingsRead {
  const refuse = (message: string): SettingsRead => ({ ok: false, message });
  if (bytes.length > SETTINGS_MAX_BYTES) {
    return refuse(`The file has ${bytes.length} bytes; a settings file has at most ${SETTINGS_MAX_BYTES} (1 MiB).`);
  }
  let text: string;
  try {
    text = new TextDecoder('utf-8', { fatal: true }).decode(bytes);
  } catch {
    return refuse('The file is not UTF-8 text, so it is not a LOSAT Web settings file.');
  }
  let document: unknown;
  try {
    document = JSON.parse(text);
  } catch (error) {
    return refuse(`The file is not JSON (${error instanceof Error ? error.message : String(error)}).`);
  }
  if (typeof document !== 'object' || document === null || Array.isArray(document)) {
    return refuse('The file is not a LOSAT Web settings file: it is not a JSON object.');
  }
  const fields = document as Record<string, unknown>;
  if (fields.format !== SETTINGS_FORMAT) {
    return refuse(`The file is not a LOSAT Web settings file: its "format" is ${shown(fields.format)}, not "${SETTINGS_FORMAT}".`);
  }
  const schema = fields.schema;
  if (typeof schema !== 'number' || !Number.isSafeInteger(schema) || schema < 1) {
    return refuse(`"schema" must be a whole number from 1 (it is ${shown(schema)}).`);
  }
  if (schema > SETTINGS_SCHEMA) {
    return refuse(`The file has schema ${schema}: a newer LOSAT Web saved it. This one reads schema ${SETTINGS_SCHEMA}.`);
  }
  const unknown = Object.keys(fields).find((key) => !FIELDS.includes(key));
  if (unknown !== undefined) return refuse(`The file has a field that schema ${SETTINGS_SCHEMA} does not have: ${shown(unknown)}.`);
  if (fields.note !== undefined && (typeof fields.note !== 'string' || fields.note.length > MAX_NOTE_CHARS)) {
    return refuse(`"note" must be a text of at most ${MAX_NOTE_CHARS} characters.`);
  }
  const program = PROGRAMS.find((p) => p.id === fields.program)?.id;
  if (program === undefined) {
    const known = PROGRAMS.map((p) => p.id).join(', ');
    return refuse(`"program" must be one of ${known} (it is ${shown(fields.program)}).`);
  }
  const options = fields.options;
  if (!Array.isArray(options)) return refuse(`"options" must be a list of words (it is ${shown(options)}).`);
  if (options.length > SETTINGS_MAX_OPTIONS) {
    return refuse(`"options" has ${options.length} words; a settings file has at most ${SETTINGS_MAX_OPTIONS}.`);
  }
  for (const [i, word] of options.entries()) {
    const where = `"options" word ${i + 1}`;
    if (typeof word !== 'string') return refuse(`${where} must be a text (it is ${shown(word)}).`);
    if (word.length > SETTINGS_MAX_OPTION_CHARS) {
      return refuse(`${where} has ${word.length} characters; a word has at most ${SETTINGS_MAX_OPTION_CHARS}.`);
    }
    if (RESERVED_FLAGS.includes(word)) {
      return refuse(`${where} is ${word}, which LOSAT Web sets itself: a settings file has no inputs, outputs or thread count among its options.`);
    }
  }
  const threads = fields.threads;
  if (threads !== 'auto' && !(typeof threads === 'number' && Number.isInteger(threads) && threads >= 1 && threads <= MAX_THREADS)) {
    return refuse(`"threads" must be "auto" or a whole number from 1 to ${MAX_THREADS} (it is ${shown(threads)}).`);
  }
  return { ok: true, settings: { program, options: Object.freeze([...(options as string[])]), threads } };
}

/** One option of an argv: its flag, its value (`true` for a flag alone), and its words. */
export interface ArgvOption {
  readonly flag: string;
  readonly value: string | true;
  readonly words: readonly string[];
}

/**
 * The options of an argv's words, by the engine's grammar: `takesValue(flag)` says whether a
 * flag that the engine describes takes a value (so `-penalty -3` is one option), and is
 * undefined for a flag that it does not describe, which then takes the next word unless that
 * word looks like an option. Words that are no option (a value without a flag, or a flag
 * without its value at the end) are `stray`.
 */
export function readOptions(
  words: readonly string[],
  takesValue: (flag: string) => boolean | undefined,
): { readonly options: readonly ArgvOption[]; readonly stray: readonly string[] } {
  const options: ArgvOption[] = [];
  const stray: string[] = [];
  let i = 0;
  while (i < words.length) {
    const flag = words[i]!;
    if (!looksLikeOption(flag)) {
      stray.push(flag);
      i += 1;
      continue;
    }
    const next = words[i + 1];
    const described = takesValue(flag);
    const withValue = described ?? (next !== undefined && !looksLikeOption(next));
    if (!withValue) {
      options.push({ flag, value: true, words: [flag] });
      i += 1;
    } else if (next === undefined) {
      stray.push(flag);
      i += 1;
    } else {
      options.push({ flag, value: next, words: [flag, next] });
      i += 2;
    }
  }
  return { options, stray };
}

/** A word that starts like an option name (`-word_size`); a negative number (`-3`) does not. */
function looksLikeOption(word: string): boolean {
  return /^-[A-Za-z_]/.test(word);
}

/** A value of the file as the message shows it: JSON, cut to a few dozen characters. */
function shown(value: unknown): string {
  if (value === undefined) return 'missing';
  const text = JSON.stringify(value) ?? String(value);
  return text.length > 60 ? `${text.slice(0, 57)}...` : text;
}
