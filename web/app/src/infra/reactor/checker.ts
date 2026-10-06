// The input check of the Data worker (ports/input-check.ts): `register` of the program on
// the serial reactor, then `release` (docs/web/abi_v2.md §4). The engine's refusal (an
// export that returned -1) is the verdict with the engine's message; a stopped instance is
// a failure of the check itself.
import type { InputRole, ProgramId } from '../../domain/programs';
import type { InputCheck, InputChecker } from '../../ports/input-check';
import { EngineCallError, ROLE_QUERY, ROLE_SUBJECT, type ReactorAbi } from './abi';

export class ReactorInputChecker implements InputChecker {
  constructor(private readonly reactor: () => Promise<ReactorAbi>) {}

  async check(program: ProgramId, role: InputRole, bytes: Uint8Array): Promise<InputCheck> {
    const abi = await this.reactor();
    let registered;
    try {
      registered = abi.register(program, role === 'query' ? ROLE_QUERY : ROLE_SUBJECT, bytes);
    } catch (error) {
      if (error instanceof EngineCallError) return { ok: false, message: error.message };
      throw error;
    }
    abi.release(registered.handle);
    return { ok: true, records: registered.records };
  }
}
