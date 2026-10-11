// The example files of the release (domain/examples.ts): fetched from the site, never made in
// the page. The browser implementation is src/infra/browser/examples.ts; tests use a fake.
import type { InputRole, ProgramId } from '../domain/programs';

/** A program's example query and subject, as Files named as `exampleFileName` says. */
export type ExampleFiles = Readonly<Record<InputRole, File>>;

export interface ExampleSource {
  /**
   * The program's example files. Rejects with a message for the user when they cannot be read
   * (the site cannot be reached and the release is not stored for offline use, for one).
   */
  files(program: ProgramId): Promise<ExampleFiles>;
}
