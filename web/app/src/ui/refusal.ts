// WA2-AX-05: a field of the form that a refusal of the argv names (the engine's, or the Web's own
// before it) is marked as refused (aria-invalid) and points to the message beside "Run LOSAT"
// (aria-describedby), so that a screen reader that comes back to the field hears why.
import type { DraftState } from '../application/draft';
import { refusedFields } from '../domain/argv';

/** The id of the refusal's message (SearchPanel), to which the fields that it names point. */
export const REFUSAL_ID = 'argv-validation-message';

/** The fields that the refusal of the draft's argv names, or undefined while the argv is not refused. */
export function refusal(state: DraftState): ReturnType<typeof refusedFields> | undefined {
  const validation = state.validation;
  return validation.state === 'invalid' ? refusedFields(validation.message, state.extra[state.program]) : undefined;
}

/** The attributes of a refused field's control, pointing to message `id`: none when `refused` is false. */
export function refusedMark(refused: boolean, id = REFUSAL_ID): { readonly 'aria-invalid'?: 'true'; readonly 'aria-describedby'?: string } {
  return refused ? { 'aria-invalid': 'true', 'aria-describedby': id } : {};
}
