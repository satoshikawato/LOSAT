// The examples of LOSAT Web (S16, docs/web/examples.md): one small query and subject per
// program, cut from public records and shipped with the release (src/examples/), for "Try an
// example" and the browser self-test. The self-test compares the engine's outputs for them
// with the SHA-256 that src/examples/self-test.json gives (written by tools/examples.py from
// the native LOSAT CLI of the engine commit). They are release assets, not research data.
import { OUTPUT_FORMATS, type OutputFormat } from './output-format';
import { PROGRAMS, type InputRole, type ProgramId } from './programs';

/** The name of an example file, as the search form shows it and the argv names it. */
export function exampleFileName(program: ProgramId, role: InputRole): string {
  return `example-${program}-${role}.fasta`;
}

/** What each example searches, for the line that says it was loaded. */
export const EXAMPLE_TOPICS: Readonly<Record<ProgramId, string>> = Object.freeze({
  blastn: 'the dnaA region of Candidatus Hepatoplasma crinochetorum against that region of three Hepatoplasma genomes',
  blastp: 'the DnaA protein of Candidatus Hepatoplasma vulgare against the DnaA proteins of three related bacteria',
  blastx: 'the dnaA region of Candidatus Hepatoplasma vulgare against the DnaA proteins of three related bacteria',
  tblastn: 'the DnaA protein of Candidatus Hepatoplasma vulgare against the dnaA regions of three related genomes',
  tblastx: 'the dnaA region of Candidatus Hepatoplasma vulgare against the dnaA regions of three related genomes',
});

/** An example file of a self-test case and the SHA-256 of its bytes. */
export interface SelfTestInput {
  readonly file: string;
  readonly sha256: string;
}

/** One output of a self-test case: what the engine commit wrote, and whether NCBI BLAST+ writes the same bytes. */
export interface SelfTestOutput {
  readonly sha256: string;
  readonly bytes: number;
  readonly ncbi: 'same' | 'differs';
}

export interface SelfTestCase {
  readonly program: ProgramId;
  /** The search form's argv for the example: the program and the two inputs, default settings. */
  readonly argv: readonly string[];
  readonly inputs: Readonly<Record<InputRole, SelfTestInput>>;
  readonly outputs: Readonly<Record<OutputFormat, SelfTestOutput>>;
}

/** src/examples/self-test.json. */
export interface SelfTestTable {
  /** The engine commit whose native CLI wrote the expected outputs, its version and the CLI's SHA-256. */
  readonly engine: { readonly commit: string; readonly version: string; readonly nativeSha256: string };
  /** The NCBI BLAST+ release the outputs were compared with. */
  readonly oracle: string;
  readonly cases: readonly SelfTestCase[];
}

const SHA256 = /^[0-9a-f]{64}$/;
const ROLES: readonly InputRole[] = ['query', 'subject'];

/** Checks the parsed self-test.json; throws an Error that says what is wrong. */
export function parseSelfTestTable(value: unknown): SelfTestTable {
  const object = (v: unknown, what: string): Record<string, unknown> => {
    if (typeof v !== 'object' || v === null || Array.isArray(v)) throw new Error(`self-test table: ${what} is not an object`);
    return v as Record<string, unknown>;
  };
  const text = (v: unknown, what: string, pattern?: RegExp): string => {
    if (typeof v !== 'string' || v === '' || v.length > 300 || (pattern !== undefined && !pattern.test(v))) {
      throw new Error(`self-test table: ${what} is not valid`);
    }
    return v;
  };
  const table = object(value, 'the table');
  if (table['format'] !== 1) throw new Error('self-test table: unknown format');
  const engine = object(table['engine'], 'engine');
  const list = table['cases'];
  if (!Array.isArray(list) || list.length === 0 || list.length > PROGRAMS.length) throw new Error('self-test table: cases is not a list of programs');
  const cases = list.map((entry, n): SelfTestCase => {
    const where = `case ${n + 1}`;
    const item = object(entry, where);
    const program = PROGRAMS.find((p) => p.id === item['program'])?.id;
    if (program === undefined) throw new Error(`self-test table: ${where}'s program is not one of LOSAT Web's`);
    const inputsObject = object(item['inputs'], `${where}'s inputs`);
    const inputs = Object.fromEntries(
      ROLES.map((role) => {
        const input = object(inputsObject[role], `${where}'s ${role}`);
        return [role, Object.freeze({ file: text(input['file'], `${where}'s ${role} file`), sha256: text(input['sha256'], `${where}'s ${role} SHA-256`, SHA256) })];
      }),
    ) as Record<InputRole, SelfTestInput>;
    const argv = item['argv'];
    const expectedArgv = [program, '-query', inputs.query.file, '-subject', inputs.subject.file];
    // The argv is the search form's for the example files (tools/examples.py); options would need the form to set them.
    if (!Array.isArray(argv) || argv.length !== expectedArgv.length || argv.some((word, i) => word !== expectedArgv[i])) {
      throw new Error(`self-test table: ${where}'s argv is not the search form's for its files`);
    }
    const outputsObject = object(item['outputs'], `${where}'s outputs`);
    const outputs = Object.fromEntries(
      OUTPUT_FORMATS.map((format) => {
        const output = object(outputsObject[String(format)], `${where}'s outfmt ${format}`);
        const bytes = output['bytes'];
        if (typeof bytes !== 'number' || !Number.isSafeInteger(bytes) || bytes < 0) throw new Error(`self-test table: ${where}'s outfmt ${format} bytes is not valid`);
        const ncbi = output['ncbi'];
        if (ncbi !== 'same' && ncbi !== 'differs') throw new Error(`self-test table: ${where}'s outfmt ${format} ncbi is not "same" or "differs"`);
        return [format, Object.freeze({ sha256: text(output['sha256'], `${where}'s outfmt ${format} SHA-256`, SHA256), bytes, ncbi })];
      }),
    ) as Record<OutputFormat, SelfTestOutput>;
    return Object.freeze({ program, argv: Object.freeze([...expectedArgv]), inputs: Object.freeze(inputs), outputs: Object.freeze(outputs) });
  });
  if (new Set(cases.map((c) => c.program)).size !== cases.length) throw new Error('self-test table: a program has two cases');
  return Object.freeze({
    engine: Object.freeze({
      commit: text(engine['commit'], 'engine commit'),
      version: text(engine['version'], 'engine version'),
      nativeSha256: text(engine['nativeSha256'], 'engine nativeSha256', SHA256),
    }),
    oracle: text(table['oracle'], 'oracle'),
    cases: Object.freeze(cases),
  });
}

/** The first characters of a SHA-256, enough to tell two apart on the screen. */
export function digestPrefix(sha256: string): string {
  return sha256.slice(0, 12);
}
