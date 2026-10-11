// Builds the site's service worker (src/sw/worker.ts) into `<site>/sw.js` for the site build
// (build/site.ts). The worker is site-level, like the boot page: not versioned, outside every
// release directory, replaced by each deployment. It is one classic script (an IIFE, no
// imports), which every browser accepts as a service worker; the build is deterministic.
import { writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { build } from 'vite';

const SW_DIR = fileURLToPath(new URL('../src/sw/', import.meta.url));
export const SERVICE_WORKER_FILE = 'sw.js';

/** The worker's code, built in memory. */
export async function buildServiceWorker(): Promise<string> {
  const result = await build({
    configFile: false,
    root: SW_DIR,
    publicDir: false,
    logLevel: 'warn',
    build: {
      write: false,
      target: 'es2022',
      modulePreload: false,
      copyPublicDir: false,
      rolldownOptions: { input: join(SW_DIR, 'worker.ts'), output: { format: 'iife', entryFileNames: SERVICE_WORKER_FILE } },
    },
  });
  const outputs = (Array.isArray(result) ? result : [result]).flatMap((output) => {
    if (!('output' in output)) throw new Error('the service worker build did not return its output');
    return output.output;
  });
  const [chunk, ...rest] = outputs;
  if (chunk === undefined || chunk.type !== 'chunk' || chunk.fileName !== SERVICE_WORKER_FILE || rest.length > 0) {
    throw new Error(`the service worker build wrote ${outputs.map((item) => item.fileName).join(', ')}, not ${SERVICE_WORKER_FILE} alone`);
  }
  if (/^\s*(?:import|export)\b/m.test(chunk.code)) throw new Error(`${SERVICE_WORKER_FILE} is not a classic script`);
  return chunk.code;
}

/** Writes `<site>/sw.js`. */
export async function writeServiceWorker(siteDir: string): Promise<void> {
  writeFileSync(join(siteDir, SERVICE_WORKER_FILE), await buildServiceWorker());
}
