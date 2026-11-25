// Holder tracking for inflated OPFS-backed session files (`tool-inputs/<token>/data`):
// each session file is held by its WindowMasker form item (`form:<key>`), briefly while a submit
// queues (`submit:<n>`), and by each queued or running tool run (`run:<runId>`) that reads it
// directly or inside `combined_input.fa`. It is dropped via `DataGateway.dropSessionFile` only
// when the last holder releases it.
import { SourceRefs } from './source-refs';

export class SessionFiles {
  private readonly refs: SourceRefs;

  constructor(drop: (token: string) => void) {
    this.refs = new SourceRefs(drop);
  }

  /** Records that `holder` holds exactly `tokens`; drops any token whose last holder left. */
  set(holder: string, tokens: Iterable<string>): void {
    this.refs.set(holder, tokens);
  }

  /** Releases every token held by `holder`. */
  release(holder: string): void {
    this.refs.release(holder);
  }
}
