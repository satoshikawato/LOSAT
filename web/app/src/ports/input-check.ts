// Input check port (plan §5.4): whether the engine reads an input, and the records it
// reads. The real implementation registers the bytes with the program's `register`
// (docs/web/abi_v2.md §4) on the Data worker's serial reactor and releases the handle at
// once, so the verdict and its message are the engine's own (for example "query record 3
// (id) has … which is not supported by LOSAT's BLASTN"). The application never rebuilds
// these rules; it shows the message and lets the user exclude the record.
import type { RecordKey } from '../domain/dataset';
import type { InputRole, ProgramId } from '../domain/programs';

export type InputCheck =
  | { readonly ok: true; readonly records: readonly RecordKey[] }
  | { readonly ok: false; readonly message: string };

export interface InputChecker {
  check(program: ProgramId, role: InputRole, bytes: Uint8Array): Promise<InputCheck>;
}
