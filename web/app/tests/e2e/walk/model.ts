// What the walk knows of the search form's inputs, kept apart from the application's code: the
// texts and files it gave each role, the records it left out, and whether the role's inputs are
// searched together or apart. A source is read as nucleotides or as proteins as the mode reads its
// role, and the records left out are kept for each reading (src/application/draft.ts `reindex`):
// a mode that reads the role the other way starts with every record of that reading. From it the walk writes the files that a run's command line names,
// for the native CLI (the oracle of walk/steps.ts): the run must have searched exactly the whole
// records that the walk chose (web/AGENTS.md rule 3), as the data layer joins them
// (src/infra/data/data-service.ts: a whole file, or the included records' bytes, with a line end
// between two inputs when the first has none). A pure module (tests/unit/walk.test.ts).

export type Role = 'query' | 'subject';
export type Reading = 'nucleotide' | 'protein';
export type InputMode = 'combined' | 'separate';

export interface ModelSource {
  readonly origin: 'paste' | 'file';
  /**
   * The name of the user's file, or of the pasted text: what the command line names when the source is the role's only input
   * with every record in. With records left out the page names a file of its own (`included_q.fa`, line WQA's D22); the walk
   * does not model that name: it takes the file that the page offers and checks its bytes against `includedBytes`.
   */
  readonly name: string;
  readonly bytes: Uint8Array;
  /** The records left out, by index, for each reading. */
  readonly excluded: Readonly<Partial<Record<Reading, ReadonlySet<number>>>>;
}

export interface RoleModel {
  /** In the form's order: the pasted text first, then the files as they were opened. */
  readonly sources: readonly ModelSource[];
  readonly mode: InputMode;
}

export interface InputModel {
  readonly query: RoleModel;
  readonly subject: RoleModel;
}

export const PASTED_NAMES: Readonly<Record<Role, string>> = Object.freeze({ query: 'query.fa', subject: 'subject.fa' });
export const COMBINED_NAMES: Readonly<Record<Role, string>> = Object.freeze({ query: 'combined_query.fa', subject: 'combined_subject.fa' });

const LF = 0x0a;
const GT = 0x3e;

export function emptyModel(): InputModel {
  return { query: { sources: [], mode: 'combined' }, subject: { sources: [], mode: 'combined' } };
}

/** The pasted text of a role becomes its first input (or goes, when it is empty). */
export function withPaste(model: InputModel, role: Role, text: string): InputModel {
  const others = model[role].sources.filter((source) => source.origin !== 'paste');
  const sources =
    text === '' ? others : [{ origin: 'paste' as const, name: PASTED_NAMES[role], bytes: new TextEncoder().encode(textareaValue(text)), excluded: {} }, ...others];
  return { ...model, [role]: { ...model[role], sources } };
}

/** A text box's value: HTML turns its line ends (CR LF, CR) into LF, so a pasted CRLF text is searched with LF. */
function textareaValue(text: string): string {
  return text.replace(/\r\n?/g, '\n');
}

export function withFiles(model: InputModel, role: Role, files: ReadonlyArray<{ name: string; bytes: Uint8Array }>): InputModel {
  const added = files.map((file) => ({ origin: 'file' as const, name: file.name, bytes: file.bytes, excluded: {} }));
  return { ...model, [role]: { ...model[role], sources: [...model[role].sources, ...added] } };
}

export function withoutSource(model: InputModel, role: Role, index: number): InputModel {
  return { ...model, [role]: { ...model[role], sources: model[role].sources.filter((_, i) => i !== index) } };
}

export function withIncluded(model: InputModel, role: Role, reading: Reading, index: number, record: number, included: boolean): InputModel {
  const sources = model[role].sources.map((source, i) => {
    if (i !== index) return source;
    const excluded = new Set(source.excluded[reading]);
    if (included) excluded.delete(record);
    else excluded.add(record);
    return { ...source, excluded: { ...source.excluded, [reading]: excluded } };
  });
  return { ...model, [role]: { ...model[role], sources } };
}

export function withMode(model: InputModel, role: Role, mode: InputMode): InputModel {
  return { ...model, [role]: { ...model[role], mode } };
}

/** The byte ranges of the records of FASTA bytes: from a line's '>' to the next record's, or the end. */
export function recordRanges(bytes: Uint8Array): Array<readonly [number, number]> {
  const starts: number[] = [];
  for (let i = 0; i < bytes.length; i++) {
    if (bytes[i] === GT && (i === 0 || bytes[i - 1] === LF)) starts.push(i);
  }
  return starts.map((start, i) => [start, starts[i + 1] ?? bytes.length] as const);
}

/** The bytes of a source that a run searches: the whole file, or its included records. */
export function includedBytes(source: ModelSource, reading: Reading): Uint8Array {
  const excluded = source.excluded[reading] ?? new Set<number>();
  if (excluded.size === 0) return source.bytes;
  const parts = recordRanges(source.bytes)
    .filter((_, index) => !excluded.has(index))
    .map(([start, end]) => source.bytes.subarray(start, end));
  return concat(parts);
}

/** The bytes of inputs searched together: joined with a line end where one has none. */
export function joinedBytes(sources: readonly ModelSource[], reading: Reading): Uint8Array {
  const parts: Uint8Array[] = [];
  for (const source of sources) {
    const bytes = includedBytes(source, reading);
    const previous = parts[parts.length - 1];
    if (bytes.length > 0 && previous !== undefined && previous[previous.length - 1] !== LF) parts.push(new Uint8Array([LF]));
    if (bytes.length > 0) parts.push(bytes);
  }
  return concat(parts);
}

export interface ExpectedInput {
  /**
   * The name that the command line gives the input where the model knows it: the user's file as given, `query.fa`, `combined_*.fa`
   * (the role may be put before it where two names would meet). Undefined for the file of one source with records left out: the page
   * names a file of its own (line WQA's D22) and the walk takes that name from the file that the page offers.
   */
  readonly name: string | undefined;
  readonly bytes: Uint8Array;
  /**
   * Whether the input is not the user's file as given: pasted text, several sources together, or a source with records left out.
   * The page must offer such a file, and must not offer one for an input that is the user's file (line WQA's D22). Decided here
   * from the walk's own record of what it gave, not by the application's code.
   */
  readonly made: boolean;
}

/** Whether the source's bytes that a run searches are not the file's own: pasted text, or records left out. */
function madeOf(source: ModelSource, reading: Reading): boolean {
  return source.origin === 'paste' || (source.excluded[reading]?.size ?? 0) > 0;
}

/**
 * The inputs that the role's runs search in a mode that reads the role as `reading`, given the
 * role's sources that the page read (`ready`, by index; a source it could not read is no input):
 * one when searched together, one per source when apart, but a source whose records are all left
 * out, which has nothing to search alone (src/application/draft.ts `inputGroups`).
 */
export function expectedInputs(role: RoleModel, roleName: Role, reading: Reading, ready: readonly number[]): ExpectedInput[] {
  const sources = role.sources.filter((_, index) => ready.includes(index));
  if (sources.length === 0) return [];
  const single = (source: ModelSource): ExpectedInput => ({
    name: source.origin === 'paste' || (source.excluded[reading]?.size ?? 0) === 0 ? source.name : undefined,
    bytes: includedBytes(source, reading),
    made: madeOf(source, reading),
  });
  if (sources.length === 1) return [single(sources[0]!)];
  if (role.mode === 'separate') {
    return sources
      .filter((source) => includedBytes(source, reading).length > 0 || recordRanges(source.bytes).length === 0)
      .map(single);
  }
  return [{ name: COMBINED_NAMES[roleName], bytes: joinedBytes(sources, reading), made: true }];
}

function concat(parts: readonly Uint8Array[]): Uint8Array {
  const out = new Uint8Array(parts.reduce((sum, part) => sum + part.length, 0));
  let at = 0;
  for (const part of parts) {
    out.set(part, at);
    at += part.length;
  }
  return out;
}

/**
 * The words of a command line as the page writes it (src/domain/argv.ts `shellWords`): words
 * apart by spaces, a word in single quotes where it needs them, a quote inside as '\''.
 */
export function parseShellWords(line: string): string[] {
  const words: string[] = [];
  let word: string | undefined;
  for (let i = 0; i < line.length; i++) {
    const c = line[i]!;
    if (c === ' ' || c === '\n' || c === '\t') {
      if (word !== undefined) words.push(word);
      word = undefined;
    } else if (c === "'") {
      const end = line.indexOf("'", i + 1);
      if (end < 0) throw new Error(`an unclosed quote in ${line}`);
      word = (word ?? '') + line.slice(i + 1, end);
      i = end;
    } else if (c === '\\' && i + 1 < line.length) {
      word = (word ?? '') + line[i + 1];
      i += 1;
    } else {
      word = (word ?? '') + c;
    }
  }
  if (word !== undefined) words.push(word);
  return words;
}
