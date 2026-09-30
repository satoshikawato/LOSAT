// Data port: temporary storage of run outputs (plan §5.6). A run's outputs are staged
// while it runs and become readable only after `commit`; cancelled or failed runs are
// discarded.
import type { OutputFormat } from '../domain/output-format';
import type { HspRecord, RunSink } from './engine';

export interface ResultSetRef {
  readonly runId: string;
  readonly byteLengths: Readonly<Record<OutputFormat, number>>;
  readonly hitCount: number;
}

export interface RunStaging extends RunSink {
  commit(): Promise<ResultSetRef>;
  discard(): Promise<void>;
}

export interface DataGateway {
  openRun(runId: string): Promise<RunStaging>;
  readOutput(runId: string, format: OutputFormat): Promise<Uint8Array>;
  readHits(runId: string): Promise<readonly HspRecord[]>;
  readDiagnostics(runId: string): Promise<string>;
  deleteRun(runId: string): Promise<void>;
}
