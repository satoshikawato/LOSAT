// The engine modules of this build (build/reactors.ts); null in a build without them.
declare module 'virtual:losat-engine' {
  import type { EngineAssets } from './infra/reactor/assets';
  export const ENGINE_ASSETS: EngineAssets | null;
}
