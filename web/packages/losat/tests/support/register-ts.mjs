// `execArgv: ['--import', <this file>]` of a worker thread of the tests: the thread then runs the package's own sources. The
// sources import each other by their `.js` names, as the bundler and TypeScript resolve them, which Node (type stripping)
// cannot find: a `.js` that is not there is tried as `.ts`.
import { registerHooks } from 'node:module';

registerHooks({
  resolve(specifier, context, nextResolve) {
    try {
      return nextResolve(specifier, context);
    } catch (error) {
      if (specifier.startsWith('.') && specifier.endsWith('.js') && error?.code === 'ERR_MODULE_NOT_FOUND') {
        return nextResolve(`${specifier.slice(0, -3)}.ts`, context);
      }
      throw error;
    }
  },
});
