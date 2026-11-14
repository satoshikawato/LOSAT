// The engine's options whose values locate or name the user's inputs: a search can run them only on
// those inputs (the warm-up of a prepared instance leaves them out, infra/engine-worker/warm-up.ts).
import type { InputRole } from './programs.js';

/** The region of an input: `-query_loc start-stop` or `-subject_loc start-stop`. */
export const REGION_FLAG: Readonly<Record<InputRole, string>> = { query: '-query_loc', subject: '-subject_loc' };

/** BLASTN's option for the WindowMasker unit counts file. */
export const WINDOW_MASKER_FLAG = '-window_masker_db';

/** BLASTN's taxonomy ID option, which reads the WindowMasker unit counts of that ID. */
export const WINDOW_MASKER_TAXID_FLAG = '-window_masker_taxid';
