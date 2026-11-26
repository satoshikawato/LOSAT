// FakeEngine: a stand-in for the real engine until the Wasm engine exists (plan §7, W0).
// Its outputs are not search results and say so on every line. It exists so that the
// application layer and UI can be built and tested against the EngineGateway contract.
import { RunOutputWriter, sha256Hex, toProgramDescription } from 'losat';
import {
  CUSTOM_OUTPUT,
  OUTPUT_FORMATS,
  parseToolDescription,
  recordMismatch,
  type ProgramId,
  type RecordKey,
  type RunOutput,
  type ToolDescription,
} from 'losat/domain';
import {
  DIAGNOSTICS_STREAM,
  HITS_STREAM,
  InputMismatchError,
  RunCancelledError,
  withRunOutputFailed,
  type EngineGateway,
  type EnginePhase,
  type EngineRunRequest,
  type HspRecord,
  type ProgramDescription,
  type RuntimeInfo,
  type StartResult,
  type ToolEngine,
  type ToolFile,
  type ToolRunOptions,
  type ToolRunResult,
  type ToolValidation,
  type ValidationResult,
} from 'losat/ports';
import { optionValue } from '../../domain/argv';
import { customOutfmt, OUTFMT_FLAG, outfmtOption } from '../../domain/output-format';
import { PROGRAMS } from '../../domain/programs';
import { toolOptionValue } from '../../domain/tool-argv';
import { WINDOWMASKER_HELP_FLAGS, WINDOWMASKER_VERSION_FLAGS } from '../../domain/tools';
import DESCRIBE from './describe.json';
import WINDOWMASKER_XMLHELP_DOC from './windowmasker-xmlhelp.json';
import { fakeReaderError, fakeRecordKeys } from './fake-fasta';

const WINDOWMASKER_XMLHELP: string = WINDOWMASKER_XMLHELP_DOC.xml;

export const FAKE_MARKER = 'FAKE ENGINE OUTPUT - not a LOSAT search result';

export interface FakeEngineOptions {
  /** Milliseconds spent in each phase, so tests can cancel a running job. */
  readonly phaseDelayMs?: number;
}

export class FakeEngine implements EngineGateway, ToolEngine {
  private readonly cancelled = new Set<string>();
  private readonly phaseDelayMs: number;
  private windowMaskerDescription: ToolDescription | undefined;

  constructor(options: FakeEngineOptions = {}) {
    this.phaseDelayMs = options.phaseDelayMs ?? 0;
  }

  /**
   * The engine's *describe* of the program, copied into describe.json (a file snapshot of
   * tests/unit/engine-runtime.test.ts: compared with the reactor when one is built, and
   * written again with `vitest -u`), so that the search form has the engine's options,
   * defaults and help without the engine.
   */
  async describe(program: ProgramId): Promise<ProgramDescription> {
    const json = (DESCRIBE as Readonly<Record<string, unknown>>)[program];
    if (json === undefined) return { program, formats: OUTPUT_FORMATS, parameters: [] };
    return toProgramDescription(json);
  }

  async describeTool(tool: string): Promise<ToolDescription> {
    if (tool !== 'windowmasker') throw new Error(`unknown tool: ${tool}`);
    this.windowMaskerDescription ??= parseToolDescription(WINDOWMASKER_XMLHELP);
    return this.windowMaskerDescription;
  }

  async validateTool(argv: readonly string[]): Promise<ToolValidation> {
    const tool = argv[0];
    if (tool !== 'windowmasker') return { ok: false, message: `unknown tool: ${tool}` };
    const description = await this.describeTool(tool);
    const byFlag = new Map(description.parameters.map((param) => [param.flag, param]));
    const seen = new Set<string>();
    for (let i = 1; i < argv.length; i++) {
      const word = argv[i]!;
      const [flag, inline] = word.split(/=(.*)/s, 2);
      if (WINDOWMASKER_HELP_FLAGS.has(flag!) || WINDOWMASKER_VERSION_FLAGS.has(flag!)) return { ok: true };
      const param = byFlag.get(flag!);
      const extraValued = flag === '-logfile' || flag === '-conffile';
      if (param === undefined && !extraValued) {
        return { ok: false, message: `Argument error: Unknown argument: "${word}"` };
      }
      seen.add(flag!);
      const takesValue = param?.takesValue ?? extraValued;
      if (takesValue) {
        let value = inline;
        if (value === undefined) {
          if (i + 1 >= argv.length) return { ok: false, message: `Argument error: Missing value for "${flag}"` };
          value = argv[++i]!;
        }
        if (param?.choices !== undefined && !param.choices.includes(value)) {
          return { ok: false, message: `Argument error: Invalid value "${value}" for "${flag}"` };
        }
        if (param?.kind === 'Integer' || param?.kind === 'Int8') {
          if (!/^-?\d+$/.test(value)) {
            return { ok: false, message: `Argument error: Invalid integer "${value}" for "${flag}"` };
          }
          const num = Number(value);
          if (param.range !== undefined && (num < param.range.min || num > param.range.max)) {
            return { ok: false, message: `Argument error: Value ${value} for "${flag}" is out of range` };
          }
        }
        if (param?.kind === 'Boolean' && !/^(?:true|false|t|f|1|0)$/i.test(value)) {
          return { ok: false, message: `Argument error: Invalid boolean "${value}" for "${flag}"` };
        }
      }
    }
    for (const [first, second] of description.excludes) {
      if (seen.has(first) && seen.has(second)) {
        return { ok: false, message: `Argument error: ${first} and ${second} are mutually exclusive` };
      }
    }
    if (!seen.has('-mk_counts') && !seen.has('-convert') && !seen.has('-ustat')) {
      return { ok: false, message: 'Argument error: One of -mk_counts, -convert, or -ustat must be specified' };
    }
    return { ok: true };
  }

  async runTool(argv: readonly string[], options: ToolRunOptions = {}): Promise<ToolRunResult> {
    const runId = `fake-tool:${crypto.randomUUID()}`;
    const checkAborted = () => {
      if (options.signal?.aborted === true) throw new RunCancelledError(runId);
    };
    const pause = () => this.pauseAbortable(options.signal);
    checkAborted();
    await pause();
    checkAborted();

    const encoder = new TextEncoder();
    const stdoutParts: BlobPart[] = [];
    const stderrParts: BlobPart[] = [];
    const collectedFiles: ToolFile[] = [];
    const finish = (status: number): ToolRunResult => ({
      status,
      stdout: new Blob(stdoutParts),
      stderr: new Blob(stderrParts),
      files: collectedFiles,
    });
    const writeStdout = async (text: string) => {
      const bytes = encoder.encode(text);
      if (options.sinks?.onStdout !== undefined) await options.sinks.onStdout(bytes);
      else stdoutParts.push(bytes);
    };
    const writeStderr = async (text: string) => {
      const bytes = encoder.encode(text);
      if (options.sinks?.onStderr !== undefined) await options.sinks.onStderr(bytes);
      else stderrParts.push(bytes);
    };
    const verdict = await this.validateTool(argv);
    if (!verdict.ok) {
      await writeStderr(`${verdict.message}\n`);
      return finish(1);
    }

    const words = argv.slice(1).map((word) => word.split('=', 1)[0]!);
    if (words.includes('-xmlhelp')) {
      await writeStdout(WINDOWMASKER_XMLHELP);
      return finish(0);
    }
    if (words.some((w) => WINDOWMASKER_HELP_FLAGS.has(w))) {
      await writeStdout(`USAGE\n  windowmasker ... (${FAKE_MARKER})\n`);
      return finish(0);
    }
    if (words.some((w) => WINDOWMASKER_VERSION_FLAGS.has(w))) {
      await writeStdout(`windowmasker: 0.2.0 (${FAKE_MARKER})\n`);
      return finish(0);
    }

    const outName = toolOptionValue(argv, '-out');
    const logName = toolOptionValue(argv, '-logfile');
    const bodyText = `# ${FAKE_MARKER}\n# losat ${argv.join(' ')}\n`;
    const logText = `Warning: ${FAKE_MARKER}\n`;

    const writeFile = async (name: string, chunks: readonly string[]) => {
      const sink = options.sinks?.onFile !== undefined ? await options.sinks.onFile(name) : undefined;
      if (sink !== undefined) {
        let offset = 0;
        for (const text of chunks) {
          const bytes = encoder.encode(text);
          await sink.write(bytes, offset);
          offset += bytes.length;
        }
        await sink.close();
      } else {
        collectedFiles.push({ name, file: new Blob([...chunks], { type: 'text/plain' }) });
      }
    };

    await pause();
    checkAborted();

    const sameOutAndLog = outName !== undefined && outName !== '-' && logName !== undefined && outName === logName;
    if (sameOutAndLog) {
      await writeFile(outName, [bodyText, logText]);
    } else {
      if (outName !== undefined && outName !== '-') await writeFile(outName, [bodyText]);
      else await writeStdout(bodyText);
      if (logName !== undefined) await writeFile(logName, [logText]);
      else await writeStderr(logText);
    }

    return finish(0);
  }

  /**
   * NCBI's application start is not imitated: no registry file, no file of the run, no dry run (its
   * `validate` refuses `-conffile`, `-logfile` and `-dryrun`, which `describe` does not list), no
   * inputs read together (a run of `-parse_deflines` reads its records as one without it), and no help
   * or version printed (`validate` refuses their words, which `describe` does not list either).
   */
  async start(argv: readonly string[] = []): Promise<StartResult> {
    const database = argv.some((word) => word === '-db' || word.startsWith('-db='));
    const valueOf = (flag: string): string | null => {
      for (let i = argv.length - 1; i >= 1; i--) {
        const w = argv[i]!;
        if (w === flag) return argv[i + 1] ?? null;
        if (w.startsWith(`${flag}=`)) return w.slice(flag.length + 1);
      }
      return null;
    };
    const queryFile = valueOf('-in_pssm') ?? valueOf('-in_msa');
    const outPssm = valueOf('-out_pssm');
    const outAsciiPssm = valueOf('-out_ascii_pssm');
    const saveEachPssm = argv.includes('-save_each_pssm');
    return {
      ok: true,
      start: {
        conffile: null,
        logFile: null,
        strategyFile: null,
        dryRun: false,
        registerInputs: false,
        prints: null,
        database,
        queryFile,
        outPssm,
        outAsciiPssm,
        saveEachPssm,
      },
    };
  }

  /**
   * Refuses an unknown program, and a word that is neither an option of the program's `describe`
   * nor the value of the option before it, with the CLI's message; the engine's other checks are
   * not imitated, and a registry file is not read.
   */
  async validate(argv: readonly string[]): Promise<ValidationResult> {
    const program = argv[0];
    if (!PROGRAMS.some((p) => p.id === program)) return { ok: false, message: `unknown program: ${program}` };
    const { parameters } = await this.describe(program as ProgramId);
    const outfmt = outfmtOption(argv);
    for (let i = 1; i < argv.length; i++) {
      const word = argv[i]!;
      if (i === outfmt?.at) {
        i += outfmt.words - 1;
        continue;
      }
      const [flag, value] = word.split(/=(.*)/s, 2);
      const option = parameters.find((p) => p.flag === flag);
      if (option === undefined) return { ok: false, message: unknownArgument(program as ProgramId, word) };
      if (option.takesValue && value === undefined) i++;
    }
    // The fake reads no registry: its runs have the data loaders on, as the CLI's without one.
    return { ok: true, dataLoaders: true };
  }

  async run(request: EngineRunRequest, output: MessagePort, onPhase: (phase: EnginePhase) => void): Promise<RuntimeInfo> {
    const { runId } = request;
    let readRecords: RuntimeInfo['readRecords'];
    for (const phase of ['preparing', 'running', 'finalizing'] as const) {
      await this.pause();
      if (this.cancelled.delete(runId)) throw new RunCancelledError(runId);
      onPhase(phase);
      if (phase === 'preparing') readRecords = this.register(request);
    }
    if ((request.query !== undefined && request.query.bytes.length === 0) || request.subject?.bytes.length === 0) {
      // A failed run keeps the warnings written before its error, as the engine's does (stream 3,
      // then `fail`; line WQA's sweep, D3, 2026-11-19), and the output that the CLI leaves before its error, the
      // argv's own, with the others dropped (line WQA's sweep, D15, 2026-11-20): here marked lines.
      const writer = new RunOutputWriter(output);
      const encoder = new TextEncoder();
      const own = ownOutput(request.argv);
      writer.write(own, encoder.encode(`# ${FAKE_MARKER}\n`));
      for (const format of OUTPUT_FORMATS) if (format !== own) writer.discard(format, 'run');
      writer.write(DIAGNOSTICS_STREAM, encoder.encode(`Warning: ${FAKE_MARKER}\n`));
      writer.fail();
      throw withRunOutputFailed(new Error('No sequence input was provided'));
    }
    // The counts file is accepted and recorded, not used: the fake searches nothing.
    const counts = request.windowMasker;
    const windowMasker =
      counts === undefined
        ? undefined
        : { name: counts.name, size: counts.file.size, sha256: await sha256Hex(new Uint8Array(await counts.file.arrayBuffer())), retained: false, loadMs: 0 };
    const writer = new RunOutputWriter(output);
    this.writeOutputs(request, writer);
    writer.end();
    const runFiles = fakePssmRunFiles(request.argv);
    // Its warning goes to the diagnostics (stream 3), as without -logfile.
    return {
      path: 'fake',
      threads: 1,
      engineBuild: 'fake-engine',
      messages: 'diagnostics',
      ...(windowMasker === undefined ? {} : { windowMasker }),
      ...(readRecords === undefined ? {} : { readRecords }),
      ...(runFiles.length === 0 ? {} : { runFiles }),
    };
  }

  cancel(runId: string): void {
    this.cancelled.add(runId);
  }

  /** Reads the inputs as the Engine worker's `register` does and checks the record tables. */
  private register(request: EngineRunRequest): RuntimeInfo['readRecords'] {
    let query: readonly RecordKey[] | undefined;
    let subject: readonly RecordKey[] | undefined;
    for (const role of ['query', 'subject'] as const) {
      const input = request[role];
      // A run of -db has no subject.
      if (input === undefined) continue;
      const actual = fakeRecordKeys(input.bytes, programOf(request), role);
      const mismatch = recordMismatch(input.records, actual);
      if (mismatch !== undefined) {
        if (input.assembled === true) {
          if (role === 'query') query = actual;
          else subject = actual;
        } else {
          throw new InputMismatchError(role, mismatch);
        }
      }
    }
    // A reader's error ends the run, as the engine's run does after the query batches before it
    // (line WQA's D10, 2026-11-19): the fake writes nothing first.
    if (request.query !== undefined) {
      const readerError = fakeReaderError(request.query.bytes, programOf(request), 'query');
      if (readerError !== undefined) throw new Error(readerError);
    }
    if (query === undefined && subject === undefined) return undefined;
    return { ...(query === undefined ? {} : { query }), ...(subject === undefined ? {} : { subject }) };
  }

  /**
   * Writes outputs with the structure of a search, so that the results screen can be built
   * and tested without the engine: up to 200 queries, each but every fourth with HSPs on
   * up to three subjects (the third subject is left out of outfmt 0, as subjects past
   * -num_alignments are), a reverse HSP on a nucleotide subject and, for BLASTN, a one-letter HSP. Every text
   * says that it is not a search result, and the values are marked FAKE.
   */
  private writeOutputs(request: EngineRunRequest, writer: RunOutputWriter): void {
    const encoder = new TextEncoder();
    const program = request.argv[0] ?? '';
    const queries =
      request.query !== undefined
        ? fakeRecordKeys(request.query.bytes, programOf(request), 'query').slice(0, 200)
        : [{ id: fakeQueryFileId(request.argv), length: 120 }];
    const subjects =
      request.subject !== undefined
        ? fakeRecordKeys(request.subject.bytes, programOf(request), 'subject').slice(0, 3)
        : [
            { id: 'db_seq_1', length: 500 },
            { id: 'db_seq_2', length: 320 },
          ];
    const nucleotide = { query: ['blastn', 'blastx', 'tblastx'].includes(program), subject: ['blastn', 'tblastn', 'tblastx'].includes(program) };
    const translated = { query: program === 'blastx' || program === 'tblastx', subject: program === 'tblastn' || program === 'tblastx' };
    let out0 = `${FAKE_MARKER}\nProgram: ${program}\n`;
    // outfmt 6 has no comment lines; this first line is what makes an exported fake outfmt 6
    // say that it is not a search result. The rows are read by their `out6` ranges.
    let out6 = `# ${FAKE_MARKER}\n`;
    let out7 = `# ${FAKE_MARKER}\n`;
    let hits = '';
    const iterations =
      program === 'psiblast'
        ? Math.max(1, Math.min(5, Number.parseInt(optionValue(request.argv, '-num_iterations') ?? '1', 10) || 1))
        : 1;
    let index = 0;
    const length = (text: string) => encoder.encode(text).length;
    queries.forEach((query, qIdx) => {
      out0 += `\nQuery= ${query.id}\n`;
      if (qIdx % 4 === 3) {
        out0 += '\n***** No hits found ***** (FAKE)\n';
        return;
      }
      for (let round = 1; round <= iterations; round++) {
        if (program === 'psiblast' && iterations > 1) out0 += `\nResults from round ${round}\n`;
        let rank = 0;
        subjects.forEach((subject, sIdx) => {
          const shown = sIdx < 2;
          const heading = `> ${subject.id}\nLength=${subject.length}\n\n`;
          const headingStart = length(out0);
          if (shown) out0 += heading;
          const count = 1 + ((qIdx + sIdx) % 2);
          for (let j = 0; j < count; j++) {
            const single = program === 'blastn' && qIdx === 2 && sIdx === 0 && j === 0;
            const span = (total: number) => (single ? 1 : Math.max(1, Math.floor(total / 2)));
            const qStart = 1 + j;
            const qEnd = Math.min(query.length, qStart + span(query.length) - 1);
            const reverse = nucleotide.subject && j === 1;
            const sFrom = 1 + j;
            const sTo = Math.min(subject.length, sFrom + span(subject.length) - 1);
            const [sStart, sEnd] = reverse ? [sTo, sFrom] : [sFrom, sTo];
            const row = [query.id, subject.id, 'FAKE', String(qEnd - qStart + 1), '0', '0', qStart, qEnd, sStart, sEnd, 'FAKE', 'FAKE'].join('\t') + '\n';
            const out6Start = length(out6);
            out6 += row;
            out7 += row;
            const sectionStart = length(out0);
            if (shown) out0 += ` Score = FAKE, Expect = FAKE (${FAKE_MARKER})\n Strand=Plus/${reverse ? 'Minus' : 'Plus'}\n\n`;
            const hit: HspRecord = {
              index,
              q_idx: qIdx,
              s_idx: sIdx,
              rank,
              raw_score: 100 - index,
              bit_score: 50 - index / 2,
              e_value: index * 1e-5,
              q_start: qStart,
              q_end: qEnd,
              s_start: sStart,
              s_end: sEnd,
              query_frame: translated.query ? 1 : null,
              subject_frame: translated.subject ? (reverse ? -1 : 1) : null,
              query_length: query.length,
              subject_length: subject.length,
              query_aligned: null,
              subject_aligned: null,
              out6: [out6Start, length(out6)],
              out0: shown ? [sectionStart, length(out0)] : null,
              out0_subject: shown ? [headingStart, headingStart + length(heading)] : null,
              ...(program === 'psiblast' ? { round } : {}),
            };
            hits += `${JSON.stringify(hit)}\n`;
            index++;
            rank++;
          }
        });
      }
    });
    writer.write(0, encoder.encode(out0));
    writer.write(6, encoder.encode(out6));
    writer.write(7, encoder.encode(`${out7}# ${index} hits found\n`));
    // The output of the argv's -outfmt (a format besides 0, 6 and 7, or 6 or 7 with fields): the
    // rows of outfmt 6, marked.
    const custom = customOutfmt(request.argv);
    if (custom !== undefined) writer.write(CUSTOM_OUTPUT, encoder.encode(`# ${FAKE_MARKER}\n# ${OUTFMT_FLAG} ${custom}\n${out6.slice(out6.indexOf('\n') + 1)}`));
    writer.write(HITS_STREAM, encoder.encode(hits));
    writer.write(DIAGNOSTICS_STREAM, encoder.encode(`Warning: ${FAKE_MARKER}\n`));
  }

  private pause(): Promise<void> {
    return new Promise((resolve) => setTimeout(resolve, this.phaseDelayMs));
  }

  private pauseAbortable(signal: AbortSignal | undefined): Promise<void> {
    if (this.phaseDelayMs <= 0) return Promise.resolve();
    return new Promise((resolve) => {
      if (signal?.aborted === true) {
        resolve();
        return;
      }
      const onAbort = () => {
        clearTimeout(timer);
        resolve();
      };
      const timer = setTimeout(() => {
        signal?.removeEventListener('abort', onAbort);
        resolve();
      }, this.phaseDelayMs);
      signal?.addEventListener('abort', onAbort, { once: true });
    });
  }
}

/**
 * The engine's argument error for a word that the fake does not read, in its form: a short usage,
 * a line of 72 `=`, an empty line and the two posts (the second with two spaces after `Error:`).
 * An option word is `Unknown argument` without its first dash and without a value joined by `=`
 * (`-foo=3` is `"foo"`, as NCBI's reader splits the word at its first `=`); a word that is not one
 * is a positional argument, of which the programs take none.
 */
function unknownArgument(program: ProgramId, word: string): string {
  const post = word.startsWith('-') ? `Unknown argument: "${word.slice(1).split('=')[0]!}"` : `Too many positional arguments (1), the offending value: ${word}`;
  return [
    'USAGE',
    `  ${program} [-h] [-help] [-query input_file] [-subject subject_input_file] [-out output_file]`,
    '',
    'DESCRIPTION',
    '   Fake engine',
    '',
    "Use '-help' to print detailed descriptions of command line arguments",
    '='.repeat(72),
    '',
    `Error: ${post}`,
    `Error:  ${post}`,
    '',
  ].join('\n');
}

/** The output that the CLI's run of the argv writes: the argv's own -outfmt output, else the format of its number (0 without one). */
function ownOutput(argv: readonly string[]): RunOutput {
  if (customOutfmt(argv) !== undefined) return CUSTOM_OUTPUT;
  const value = outfmtOption(argv)?.value.trim();
  const number = value === undefined ? 0 : Number(value.split(/\s+/)[0]);
  return number === 6 || number === 7 ? number : 0;
}

/** The program of a run: the first word of its argv. */
function programOf(request: EngineRunRequest): ProgramId {
  return request.argv[0] as ProgramId;
}

function fakeArgValue(argv: readonly string[], flag: string): string | undefined {
  for (let i = argv.length - 1; i >= 1; i--) {
    const word = argv[i]!;
    if (word === flag) return argv[i + 1];
    if (word.startsWith(`${flag}=`)) return word.slice(flag.length + 1);
  }
  return undefined;
}

function fakeQueryFileId(argv: readonly string[]): string {
  return fakeArgValue(argv, '-in_pssm') ?? fakeArgValue(argv, '-in_msa') ?? 'query_pssm';
}

function fakePssmRunFiles(argv: readonly string[]): readonly { readonly name: string; readonly bytes: Uint8Array; readonly option: 'out_pssm' | 'out_ascii_pssm' }[] {
  const encoder = new TextEncoder();
  const files: { readonly name: string; readonly bytes: Uint8Array; readonly option: 'out_pssm' | 'out_ascii_pssm' }[] = [];
  const outPssm = fakeArgValue(argv, '-out_pssm');
  if (outPssm !== undefined && outPssm !== '') {
    files.push({ name: outPssm, bytes: encoder.encode(`PssmWithParameters ::= {\n  -- ${FAKE_MARKER}\n}\n`), option: 'out_pssm' });
  }
  const outAscii = fakeArgValue(argv, '-out_ascii_pssm');
  if (outAscii !== undefined && outAscii !== '') {
    files.push({ name: outAscii, bytes: encoder.encode(`Last position-specific scoring matrix computed (${FAKE_MARKER})\n`), option: 'out_ascii_pssm' });
  }
  return files;
}
