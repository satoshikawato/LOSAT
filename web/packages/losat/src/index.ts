// The npm package `losat`: the host of LOSAT's two Wasm reactors (docs/web/abi_v2.md) for browsers.
// LOSAT Web (web/app) runs the engine through these modules. The layers follow web/AGENTS.md
// rule 5: `domain` depends on nothing, `ports` on the domain, `infra` implements the ports.
export * from './domain/index.js';
export * from './ports/index.js';
export { sha256Hex } from './infra/digest.js';
export {
  EngineAllocationError,
  EngineCallError,
  EngineStoppedError,
  FILE_REGISTRY,
  FILE_WINDOW_MASKER_COUNTS,
  ReactorAbi,
  ROLE_QUERY,
  ROLE_SUBJECT,
  RUN_FILES_STREAM,
  type AbiExports,
  type BegunInflation,
  type ChunkBuffer,
  type CountsFileInfo,
  type FileInfo,
  type RegistryFileInfo,
  type RegisteredInput,
  type StreamSink,
} from './infra/reactor/abi.js';
export { inspectReactor, type MemoryLimits, type ReactorArtifact, type ReactorKind } from './infra/reactor/artifact.js';
export { engineBuildName, type EngineAssets, type ReactorAsset } from './infra/reactor/assets.js';
export { ReactorInputChecker } from './infra/reactor/checker.js';
export { reactorControl, toProgramDescription, type EngineControl } from './infra/reactor/control.js';
export { ReactorInflater } from './infra/reactor/inflater.js';
export { compileReactor, instantiateSerial, loadReactor, type CompiledReactor, type ReactorInstance } from './infra/reactor/instance.js';
export { ReactorScanner, reopening } from './infra/reactor/scanner.js';
export { WasmEngine, type RuntimeEnd, type WasmEngineOptions } from './infra/engine-worker/gateway.js';
export { DEFAULT_RENEWAL, type RenewalLimits } from './infra/engine-worker/policy.js';
export { RunOutputReceiver, type RunOutputResult } from './infra/run-output/receiver.js';
export { RunOutputWriter } from './infra/run-output/writer.js';
