// Input check port (plan §5.4): whether the engine reads an input, and the records it
// reads. The real implementation registers the bytes with the program's `register`
// (docs/web/abi_v2.md §4) on the Data worker's serial reactor and releases the handle at
// once, so the verdict and its message are the engine's own (for example "BLAST query
// error: CFastaReader: Near line 7, there's a line that doesn't look like plausible data,
// ..."). The application never rebuilds these rules; it shows the message and lets the
// user exclude the record that the message points at.
import type { RecordKey } from '../domain/dataset';
import type { InputRole, ProgramId } from '../domain/programs';

export type InputCheck =
  | { readonly ok: true; readonly records: readonly RecordKey[] }
  | {
      readonly ok: false;
      readonly message: string;
      /**
       * The 0-based position, among the records of the checked input, of the record that
       * holds the line the message names (DatasetStore.checkInput); absent when the message
       * names no line, or a line outside every record.
       */
      readonly recordPosition?: number;
    };

export interface InputChecker {
  check(program: ProgramId, role: InputRole, bytes: Uint8Array): Promise<InputCheck>;
}
