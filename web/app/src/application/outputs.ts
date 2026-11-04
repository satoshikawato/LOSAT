// The Outputs of a run (plan §5.8): one compatibility output at a time, byte for byte as the
// engine wrote it, and its export. The stored results can go while they are shown (the browser
// clears the site's data): a read or an export that fails says why where the output is, and the
// text of the format shown before never stays under the format chosen (WA2-BR-01).
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
  /** The output chosen could not be read: there is no text to show. */
  readonly failed: boolean;
  /** Why the latest export was not saved, else why the output could not be read; '' when neither. */
  readonly problem: string;
}

/** The output shown in the Outputs: the latest choice of a run and a format wins. */
export class OutputReader {
  readonly state = new Store<OutputState>({ shown: '', text: '', failed: false, problem: '' });
  private reads = 0;
  private exports = 0;
  private readProblem = '';
  private saveProblem = '';

  constructor(private readonly outputs: OutputPort) {}

  /** Reads one output of a completed run. What is shown stays until it is read; an export of the output before says no more. */
  async show(run: RunView, format: OutputFormat): Promise<void> {
    const read = ++this.reads;
    this.exports++;
    this.saveProblem = '';
    this.publish({ shown: '' });
    try {
      const text = await this.outputs.readOutput(run.snapshot.runId, format);
      if (read !== this.reads) return;
      this.readProblem = '';
      this.publish({ shown: `${run.snapshot.number}:${format}`, text, failed: false });
    } catch (error) {
      if (read !== this.reads) return;
      this.readProblem = unreadableResults(run, error);
      this.publish({ text: '', failed: true });
    }
  }

  /** Saves an output as a file; one that is not saved says why, until the next export or choice. */
  async export(run: RunView, format: OutputFormat): Promise<void> {
    const attempt = ++this.exports;
    this.saveProblem = '';
    this.publish({});
    const problem = await saveOutput(this.outputs, run, format);
    if (attempt !== this.exports) return;
    this.saveProblem = problem;
    this.publish({});
  }

  private publish(change: Partial<Omit<OutputState, 'problem'>>): void {
    this.state.set({ ...this.state.get(), ...change, problem: this.saveProblem || this.readProblem });
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
