// The two engine modules of a build (docs/web/abi_v2.md §2). build/reactors.ts copies
// them into the build and provides this description as `virtual:losat-engine`.

export interface ReactorAsset {
  /** Same-origin URL of the module. */
  readonly url: string;
  /** Lower-case hex SHA-256 of the module; checked before the module is compiled. */
  readonly sha256: string;
  readonly size: number;
  /** SHA-256 of the artifact that web/adapter/tools/build_reactors.py built. */
  readonly artifactSha256: string;
  /** How the served module differs from the artifact (the threaded module's shared-memory guard). */
  readonly transform?: string;
}

export interface EngineAssets {
  /** losat-web-serial.wasm (wasm32-wasip1). */
  readonly serial: ReactorAsset;
  /** losat-web-threads.wasm (wasm32-wasip1-threads). */
  readonly threads: ReactorAsset;
}

/**
 * A short name of a module for records and messages: the artifact's file name and digest
 * prefix, and whether the served module is guarded.
 */
export function engineBuildName(kind: 'serial' | 'threads', asset: ReactorAsset): string {
  const guarded = asset.transform === undefined ? '' : ' (shared-memory guard)';
  return `losat-web-${kind}.wasm sha256:${asset.artifactSha256.slice(0, 16)}${guarded}`;
}
