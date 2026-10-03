// Whether the engine E2E tests (S09) can run: the build has the engine modules.
import { findReactors } from '../../../build/reactors';

export const ENGINE = findReactors();
export const NO_ENGINE_REASON =
  'the build has no engine: set LOSAT_WEB_REACTORS to the output of web/adapter/tools/build_reactors.py';
