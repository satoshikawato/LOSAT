// FakeEngine: a stand-in for the real engine until the Wasm engine exists (plan §7, W0).
// Its outputs are not search results and say so on every line. It exists so that the
// application layer and UI can be built and tested against the EngineGateway contract.
import { recordMismatch } from '../../domain/dataset';
import { OUTPUT_FORMATS } from '../../domain/output-format';
import { indexParser, PROGRAMS, type InputRole, type ProgramId } from '../../domain/programs';
import {
  InputMismatchError,
  RunCancelledError,
  type EngineGateway,
  type EnginePhase,
  type EngineRunRequest,
  type HspRecord,
  type ProgramDescription,
  type RuntimeInfo,
  type ValidationResult,
} from '../../ports/engine';
import { DIAGNOSTICS_STREAM, HITS_STREAM } from '../../ports/run-output';
import { toProgramDescription } from '../reactor/control';
import { RunOutputWriter } from '../run-output/writer';
import DESCRIBE from './describe.json';
import { fakeRecordKeys } from './fake-fasta';

export const FAKE_MARKER = 'FAKE ENGINE OUTPUT - not a LOSAT search result';
/** The engine build that the FakeEngine's runs record. */
export const FAKE_ENGINE_BUILD = 'fake-engine';

export interface FakeEngineOptions {
  /** Milliseconds spent in each phase, so tests can cancel a running job. */
  readonly phaseDelayMs?: number;
}

export class FakeEngine implements EngineGateway {
  private readonly cancelled = new Set<string>();
  private readonly phaseDelayMs: number;

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

  async validate(argv: readonly string[]): Promise<ValidationResult> {
    const program = argv[0];
    if (!PROGRAMS.some((p) => p.id === program)) return { ok: false, message: `unknown program: ${program}` };
    // As ABI v2 does until session SX (web/adapter/src/run.rs `Program::parse`).
    if (program === 'blastx') return { ok: false, message: 'blastx is not available in LOSAT Web ABI v2 yet' };
    return { ok: true };
  }

  async run(request: EngineRunRequest, output: MessagePort, onPhase: (phase: EnginePhase) => void): Promise<RuntimeInfo> {
    const { runId } = request;
    for (const phase of ['preparing', 'running', 'finalizing'] as const) {
      await this.pause();
      if (this.cancelled.delete(runId)) throw new RunCancelledError(runId);
      onPhase(phase);
      if (phase === 'preparing') this.register(request);
    }
    if (request.query.bytes.length === 0 || request.subject.bytes.length === 0) {
      throw new Error('No sequence input was provided');
    }
    const writer = new RunOutputWriter(output);
    this.writeOutputs(request, writer);
    writer.end();
    return { path: 'fake', threads: 1, engineBuild: FAKE_ENGINE_BUILD };
  }

  cancel(runId: string): void {
    this.cancelled.add(runId);
  }

  /** Reads the inputs as the Engine worker's `register` does and checks the record tables. */
  private register(request: EngineRunRequest): void {
    for (const role of ['query', 'subject'] as const) {
      const input = request[role];
      const mismatch = recordMismatch(input.records, recordsOf(request, role));
      if (mismatch !== undefined) throw new InputMismatchError(role, mismatch);
    }
  }

  /**
   * Writes outputs with the structure of a search, so that the results screen can be built
   * and tested without the engine: up to 200 queries, each but every fourth with HSPs on
   * up to three subjects (the third subject is left out of outfmt 0, as subjects past
   * -num_alignments are), a reverse HSP and, for BLASTN, a one-letter HSP. Every text
   * says that it is not a search result, and the values are marked FAKE. The first HSP of
   * each pair has fake aligned rows (`fakeAlignedRows`), the second none, so that the
   * candidate tray's alignment export meets both (S14).
   */
  private writeOutputs(request: EngineRunRequest, writer: RunOutputWriter): void {
    const encoder = new TextEncoder();
    const program = request.argv[0] ?? '';
    const queries = recordsOf(request, 'query').slice(0, 200);
    const subjects = recordsOf(request, 'subject').slice(0, 3);
    const nucleotide = { query: program === 'blastn' || program === 'tblastx', subject: program !== 'blastp' };
    const translated = { query: program === 'tblastx', subject: program === 'tblastn' || program === 'tblastx' };
    let out0 = `${FAKE_MARKER}\nProgram: ${program}\n`;
    // outfmt 6 has no comment lines; this first line is what makes an exported fake outfmt 6
    // say that it is not a search result. The rows are read by their `out6` ranges.
    let out6 = `# ${FAKE_MARKER}\n`;
    let out7 = `# ${FAKE_MARKER}\n`;
    let hits = '';
    let index = 0;
    const length = (text: string) => encoder.encode(text).length;
    queries.forEach((query, qIdx) => {
      out0 += `\nQuery= ${query.id}\n`;
      if (qIdx % 4 === 3) {
        out0 += '\n***** No hits found ***** (FAKE)\n';
        return;
      }
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
            subject_length: subject.length,
            ...(j === 0
              ? alignedFields(fakeAlignedRows(qEnd - qStart + 1, sTo - sFrom + 1))
              : { query_aligned: null, subject_aligned: null }),
            out6: [out6Start, length(out6)],
            out0: shown ? [sectionStart, length(out0)] : null,
            out0_subject: shown ? [headingStart, headingStart + length(heading)] : null,
          };
          hits += `${JSON.stringify(hit)}\n`;
          index++;
          rank++;
        }
      });
    });
    writer.write(0, encoder.encode(out0));
    writer.write(6, encoder.encode(out6));
    writer.write(7, encoder.encode(`${out7}# ${index} hits found\n`));
    writer.write(HITS_STREAM, encoder.encode(hits));
    writer.write(DIAGNOSTICS_STREAM, encoder.encode(`Warning: ${FAKE_MARKER}\n`));
  }

  private pause(): Promise<void> {
    return new Promise((resolve) => setTimeout(resolve, this.phaseDelayMs));
  }
}

/**
 * The FakeEngine's aligned rows of an HSP: "FAKE" repeated over the HSP's length on each
 * sequence, the shorter row ended with gaps, so that the rows are as long as each other. They
 * are not an alignment, and say so.
 */
export function fakeAlignedRows(queryLength: number, subjectLength: number): { readonly query: string; readonly subject: string } {
  const width = Math.max(queryLength, subjectLength);
  const row = (length: number) => 'FAKE'.repeat(Math.ceil(length / 4)).slice(0, length) + '-'.repeat(width - length);
  return { query: row(queryLength), subject: row(subjectLength) };
}

const alignedFields = (rows: { readonly query: string; readonly subject: string }) => ({ query_aligned: rows.query, subject_aligned: rows.subject });

/** The records of a role's input, read with the role's kind as the engine's `register` reads them. */
function recordsOf(request: EngineRunRequest, role: InputRole) {
  return fakeRecordKeys(request[role].bytes, indexParser(request.argv[0] as ProgramId, role));
}
