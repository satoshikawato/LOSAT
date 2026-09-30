// FakeEngine: a stand-in for the real engine until the Wasm engine exists (plan §7, W0).
// Its outputs are not search results and say so on every line. It exists so that the
// application layer and UI can be built and tested against the EngineGateway contract.
import { OUTPUT_FORMATS } from '../../domain/output-format';
import { PROGRAMS, type ProgramId } from '../../domain/programs';
import {
  RunCancelledError,
  type EngineGateway,
  type EnginePhase,
  type EngineRunRequest,
  type HspRecord,
  type ProgramDescription,
  type RunSink,
  type RuntimeInfo,
  type ValidationResult,
} from '../../ports/engine';

export const FAKE_MARKER = 'FAKE ENGINE OUTPUT - not a LOSAT search result';

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

  async describe(program: ProgramId): Promise<ProgramDescription> {
    return { program, formats: OUTPUT_FORMATS, parameters: [] };
  }

  async validate(argv: readonly string[]): Promise<ValidationResult> {
    const program = argv[0];
    if (!PROGRAMS.some((p) => p.id === program)) return { ok: false, message: `unknown program: ${program}` };
    return { ok: true };
  }

  async run(request: EngineRunRequest, sink: RunSink, onPhase: (phase: EnginePhase) => void): Promise<RuntimeInfo> {
    const { runId } = request;
    for (const phase of ['preparing', 'running', 'finalizing'] as const) {
      await this.pause();
      if (this.cancelled.delete(runId)) throw new RunCancelledError(runId);
      onPhase(phase);
    }
    if (request.query.length === 0 || request.subject.length === 0) {
      throw new Error('No sequence input was provided');
    }
    this.writeOutputs(request, sink);
    return { path: 'fake', threads: 1, engineBuild: 'fake-engine' };
  }

  cancel(runId: string): void {
    this.cancelled.add(runId);
  }

  private writeOutputs(request: EngineRunRequest, sink: RunSink): void {
    const encoder = new TextEncoder();
    const program = request.argv[0] ?? '';
    const row = `${FAKE_MARKER}\tquery\tsubject\n`;
    const out0 = `${FAKE_MARKER}\nProgram: ${program}\n`;
    sink.write(0, encoder.encode(out0));
    sink.write(6, encoder.encode(row));
    sink.write(7, encoder.encode(`# ${FAKE_MARKER}\n# 1 hits found\n${row}`));
    const hit: HspRecord = {
      index: 0,
      q_idx: 0,
      s_idx: 0,
      rank: 0,
      raw_score: 0,
      bit_score: 0,
      e_value: 0,
      q_start: 1,
      q_end: 1,
      s_start: 1,
      s_end: 1,
      query_frame: null,
      subject_frame: null,
      subject_length: null,
      query_aligned: null,
      subject_aligned: null,
      out6: [0, encoder.encode(row).length],
      out0: [0, encoder.encode(out0).length],
      out0_subject: null,
    };
    sink.hits([hit]);
  }

  private pause(): Promise<void> {
    return new Promise((resolve) => setTimeout(resolve, this.phaseDelayMs));
  }
}
