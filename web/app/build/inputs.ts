// Every file and directory of the repository that a build of the site reads (S16 decision 1,
// docs/web/release_manifest.md "The release ID"). The release ID is made from the commit and a
// check of exactly these paths for changes that are not committed (build/release.ts), so the
// plugins take their paths from here: a path that a plugin reads and this list leaves out
// would change a release's bytes under a clean ID. Paths are relative to the repository root.
// Not listed, because the ID covers them otherwise: the engine modules (their digest is part
// of the ID), and node_modules (package-lock.json, in web/app, pins it).
import { fileURLToPath } from 'node:url';

/** The repository root: the directory with web/, docs/ and LOSAT/. */
export const REPOSITORY = fileURLToPath(new URL('../../..', import.meta.url));

/** The application: its sources, the build plugins, public/_headers, package-lock.json. */
export const APP = 'web/app';

/** The browser runtime's verification cells (build/verification.ts). */
export const VERIFICATION_CELLS = 'docs/web/verification_cells.tsv';

/** The certification records that the verification table is generated from (build/verification.ts), in the table's order. */
export const CERTIFICATION_RECORDS = {
  outfmt0: 'LOSAT/tests/outfmt0_manifest.tsv',
  blastnRegression: 'LOSAT/tests/fixtures/blastn_regression/manifest.tsv',
  tblastxRegression: 'LOSAT/tests/fixtures/tblastx_regression/manifest.tsv',
  rangeRegression: 'LOSAT/tests/fixtures/range_regression/manifest.tsv',
  blastpGate: 'LOSAT/tests/blastp_v010_parity_manifest.tsv',
  tblastxGate: 'LOSAT/tests/tblastx_v010_parity_manifest.tsv',
  tblastnMatrix: 'docs/evidence/tlosan_stage_g/matrix_162.jsonl',
  blastpSweep: 'docs/evidence/losat_web_e2e/run-20261004T163746Z/sweeps/after-blastp.tsv',
  tblastnSweep: 'docs/evidence/losat_web_e2e/run-20261004T163746Z/sweeps/after-tblastn.tsv',
  tblastxSweep: 'docs/evidence/losat_web_e2e/run-20261004T163746Z/sweeps/after-tblastx.tsv',
} as const;

/** The shared-memory guard applied to the threaded engine module (build/reactors.ts). */
export const SHARED_MEMORY_GUARD = 'LOSAT/tests/wasi_shared_memory.js';

/** What the dirty check of the release ID covers. */
export const BUILD_INPUTS: readonly string[] = [APP, VERIFICATION_CELLS, ...Object.values(CERTIFICATION_RECORDS), SHARED_MEMORY_GUARD];
