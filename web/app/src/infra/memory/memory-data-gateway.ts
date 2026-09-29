// In-memory DataGateway. It is the Memory implementation of the storage contract
// (plan §5.6); the OPFS implementation and the Data worker are added in W2.
import type { DataGateway, ResultSetRef, RunStaging } from '../../ports/data';
import { OUTPUT_FORMATS, type OutputFormat } from '../../domain/output-format';
import type { HspRecord } from '../../ports/engine';

interface StoredRun {
  readonly outputs: Record<OutputFormat, Uint8Array[]>;
  readonly hits: HspRecord[];
  readonly diagnostics: Uint8Array[];
  committed: boolean;
}

export class MemoryDataGateway implements DataGateway {
  private readonly runs = new Map<string, StoredRun>();

  async openRun(runId: string): Promise<RunStaging> {
    if (this.runs.has(runId)) throw new Error(`run ${runId} already exists`);
    const run: StoredRun = { outputs: { 0: [], 6: [], 7: [] }, hits: [], diagnostics: [], committed: false };
    this.runs.set(runId, run);
    const assertOpen = () => {
      if (run.committed || this.runs.get(runId) !== run) throw new Error(`run ${runId} is closed`);
    };
    return {
      write: (format, chunk) => {
        assertOpen();
        run.outputs[format].push(chunk.slice());
      },
      hits: (records) => {
        assertOpen();
        run.hits.push(...records);
      },
      diagnostics: (chunk) => {
        assertOpen();
        run.diagnostics.push(chunk.slice());
      },
      commit: async (): Promise<ResultSetRef> => {
        assertOpen();
        run.committed = true;
        const byteLengths = { 0: 0, 6: 0, 7: 0 };
        for (const format of OUTPUT_FORMATS) {
          byteLengths[format] = run.outputs[format].reduce((sum, chunk) => sum + chunk.length, 0);
        }
        return { runId, byteLengths, hitCount: run.hits.length };
      },
      discard: async () => {
        if (this.runs.get(runId) === run) this.runs.delete(runId);
      },
    };
  }

  async readOutput(runId: string, format: OutputFormat): Promise<Uint8Array> {
    return concat(this.committedRun(runId).outputs[format]);
  }

  async readHits(runId: string): Promise<readonly HspRecord[]> {
    return Object.freeze([...this.committedRun(runId).hits]);
  }

  async readDiagnostics(runId: string): Promise<string> {
    return new TextDecoder().decode(concat(this.committedRun(runId).diagnostics));
  }

  async deleteRun(runId: string): Promise<void> {
    this.runs.delete(runId);
  }

  private committedRun(runId: string): StoredRun {
    const run = this.runs.get(runId);
    if (run === undefined || !run.committed) throw new Error(`run ${runId} has no committed result`);
    return run;
  }
}

function concat(chunks: readonly Uint8Array[]): Uint8Array {
  const bytes = new Uint8Array(chunks.reduce((sum, chunk) => sum + chunk.length, 0));
  let offset = 0;
  for (const chunk of chunks) {
    bytes.set(chunk, offset);
    offset += chunk.length;
  }
  return bytes;
}
