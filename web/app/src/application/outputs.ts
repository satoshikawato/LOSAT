// The Outputs of a run (plan §5.8): one compatibility output at a time, byte for byte as the
// engine wrote it, and its export. The stored results can go while they are shown (the browser
// clears the site's data): a read or an export that fails says why where the output is, and once
// the read of the format chosen has failed, the text of the format shown before goes (WA2-BR-01).
import type { OutputFormat } from '../domain/output-format';
import type { Coordinator, RunView } from './coordinator';
import { unreadableResults } from './results';
import { Store } from './store';

export type OutputPort = Pick<Coordinator, 'readOutput' | 'exportOutput'>;

export interface OutputState {
  /** The run and format of `text` ("<number>:<format>") once it is read; '' while it is read, or when it could not be. */
  readonly shown: string;
  /** The output as the engine wrote it; '' when it could not be read. */
  readonly text: string;
  /** The run and format ("<number>:<format>") whose read failed, which leaves no text to show; '' when the latest read did not fail. */
  readonly failed: string;
  /** Why the latest export was not saved, else why the output could not be read; '' when neither. */
  readonly problem: string;
}

/** Saves outputs as files (Export, "Download All"), and says why the latest one was not saved. */
export class OutputSaver {
  /** Why the latest save was not saved; '' when it was, or when none is to be told. */
  readonly problem = new Store('');
  private saves = 0;

  constructor(private readonly outputs: Pick<OutputPort, 'exportOutput'>) {}

  /** Saves one output of a completed run. The reason shown stays until this save ends; a later save or a reset wins. */
  async save(run: RunView, format: OutputFormat): Promise<void> {
    const save = ++this.saves;
    const problem = await saveOutput(this.outputs, run, format);
    if (save === this.saves) this.set(problem);
  }

  /** Forgets the reason, and what the saves still running would say: another run or format is chosen. */
  reset(): void {
    this.saves++;
    this.set('');
  }

  private set(problem: string): void {
    if (problem !== this.problem.get()) this.problem.set(problem);
  }
}

/** The output shown in the Outputs: the latest choice of a run and a format wins. */
export class OutputReader {
  readonly state = new Store<OutputState>({ shown: '', text: '', failed: '', problem: '' });
  private readonly saver: OutputSaver;
  private reads = 0;
  private readProblem = '';

  constructor(private readonly outputs: OutputPort) {
    this.saver = new OutputSaver(outputs);
    this.saver.problem.subscribe(() => this.publish({}));
  }

  /** Reads one output of a completed run. What is shown stays until it is read; an export of the output before says no more. */
  async show(run: RunView, format: OutputFormat): Promise<void> {
    const read = ++this.reads;
    this.saver.reset();
    this.publish({ shown: '' });
    const choice = `${run.snapshot.number}:${format}`;
    let text: string;
    try {
      text = await this.outputs.readOutput(run.snapshot.runId, format);
    } catch (error) {
      if (read !== this.reads) return;
      this.readProblem = unreadableResults(run, error);
      this.publish({ text: '', failed: choice });
      return;
    }
    if (read !== this.reads) return;
    this.readProblem = '';
    this.publish({ shown: choice, text, failed: '' });
  }

  /** Saves an output as a file; one that is not saved says why, until the next export or choice. */
  export(run: RunView, format: OutputFormat): Promise<void> {
    return this.saver.save(run, format);
  }

  private publish(change: Partial<Omit<OutputState, 'problem'>>): void {
    this.state.set({ ...this.state.get(), ...change, problem: this.saver.problem.get() || this.readProblem });
  }
}

/** Saves one output of a completed run as a file: '' when it is saved, else why it is not (as an Example's copy says it). */
export async function saveOutput(outputs: Pick<OutputPort, 'exportOutput'>, run: RunView, format: OutputFormat): Promise<string> {
  try {
    await outputs.exportOutput(run.snapshot.runId, format);
    return '';
  } catch (error) {
    return `Not saved: ${error instanceof Error ? error.message : String(error)}`;
  }
}
