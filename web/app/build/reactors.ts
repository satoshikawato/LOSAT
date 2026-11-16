// The engine reactors (docs/web/abi_v2.md §2) enter the build from the directory that
// web/adapter/tools/build_reactors.py writes: `LOSAT_WEB_REACTORS=<output dir>`. The
// package's reader (web/packages/losat/build/reactors.ts) checks the two modules against the
// SHA-256 of that build's artifacts.json and guards the threaded one; the plugin copies them into
// the build as engine/losat-web-<serial|threads>-<sha256 prefix>.wasm, and provides their URLs and
// digests to the application as the module `virtual:losat-engine`. Without the variable, the
// build has no engine (ENGINE_ASSETS is null) and the application uses the FakeEngine, with its
// banner (web/AGENTS.md rule 7).
import type { IncomingMessage, ServerResponse } from 'node:http';
import type { Plugin } from 'vite';
import { REACTOR_KINDS, REACTORS_ENV, readReactors, type ReactorBytes } from '../../packages/losat/build/reactors.ts';

export { REACTORS_ENV };

const VIRTUAL_ID = 'virtual:losat-engine';
const RESOLVED_ID = `\0${VIRTUAL_ID}`;
const NAMES = REACTOR_KINDS;

export interface ReactorFile extends ReactorBytes {
  /** The file name in the build. */
  readonly name: string;
}

export interface Reactors {
  readonly dir: string;
  readonly serial: ReactorFile;
  readonly threads: ReactorFile;
}

/** The reactors named by LOSAT_WEB_REACTORS, or undefined when it is not set. */
export function findReactors(env: NodeJS.ProcessEnv = process.env): Reactors | undefined {
  const value = env[REACTORS_ENV];
  if (value === undefined || value === '') return undefined;
  const build = readReactors(value);
  const named = (kind: (typeof NAMES)[number]): ReactorFile => ({
    ...build[kind],
    name: `engine/losat-web-${kind}-${build[kind].sha256.slice(0, 16)}.wasm`,
  });
  return { dir: build.dir, serial: named('serial'), threads: named('threads') };
}

/**
 * The plugin. `emit` copies the modules into the build; worker bundles use the same
 * plugin with `emit: false`, because the main bundle already copies them. Pass the same
 * `reactors` (findReactors()) to both, so that the modules are read and guarded once.
 */
export function losatEngine(options: { readonly reactors?: Reactors | undefined; readonly emit?: boolean } = {}): Plugin {
  const reactors = 'reactors' in options ? options.reactors : findReactors();
  const emit = options.emit ?? true;
  let base = '/';
  let building = false;
  return {
    name: 'losat-engine',
    configResolved(config) {
      base = config.base;
      building = config.command === 'build';
    },
    resolveId(id) {
      return id === VIRTUAL_ID ? RESOLVED_ID : undefined;
    },
    load(id) {
      if (id !== RESOLVED_ID) return undefined;
      if (reactors === undefined) return 'export const ENGINE_ASSETS = null;\n';
      const asset = (file: ReactorFile) => ({
        url: `${base}${file.name}`,
        sha256: file.sha256,
        size: file.bytes.length,
        artifactSha256: file.artifactSha256,
        ...(file.transform === undefined ? {} : { transform: file.transform }),
      });
      const assets = { serial: asset(reactors.serial), threads: asset(reactors.threads) };
      return `export const ENGINE_ASSETS = Object.freeze(${JSON.stringify(assets)});\n`;
    },
    buildStart() {
      if (!emit || !building || reactors === undefined) return;
      for (const kind of NAMES) {
        const file = reactors[kind];
        this.emitFile({ type: 'asset', fileName: file.name, source: file.bytes });
      }
    },
    configureServer(server) {
      // The development server serves the modules from the reactor directory.
      server.middlewares.use((request: IncomingMessage, response: ServerResponse, next: () => void) => {
        const file = reactors === undefined ? undefined : NAMES.map((kind) => reactors[kind]).find((f) => request.url === `${base}${f.name}`);
        if (file === undefined) {
          next();
          return;
        }
        response.setHeader('Content-Type', 'application/wasm');
        response.end(file.bytes);
      });
    },
  };
}
