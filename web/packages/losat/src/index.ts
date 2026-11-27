// The npm package `losat`: the host of LOSAT's two Wasm reactors (docs/web/abi_v2.md) for browsers.
// LOSAT Web (web/app) runs the engine through these modules. The layers follow web/AGENTS.md
// rule 5: `domain` depends on nothing, `ports` on the domain, `infra` implements the ports.
export * from './domain/index.js';
export * from './ports/index.js';
export { sha256Hex } from './infra/digest.js';
export {
  ABI_VERSION,
  discardedStream,
  discardSignal,
  DISCARD_STREAM,
  EngineAllocationError,
  EngineCallError,
  EngineStoppedError,
  FILE_REGISTRY,
  FILE_SEARCH_STRATEGY,
  FILE_CHUNK_BYTES,
  FILE_WINDOW_MASKER_COUNTS,
  loadFile,
  MESSAGES_STREAM,
  messagesSignal,
  PLACE_NOWHERE,
  PLACE_RUN,
  ReactorAbi,
  TOOL_FILES_STREAM,
  ROLE_QUERY,
  ROLE_SUBJECT,
  RUN_FILES_STREAM,
  type AbiExports,
  type BegunInflation,
  type ChunkBuffer,
  type CountsFileInfo,
  type DiscardSignal,
  type FileInfo,
  type FileLoader,
  type RegistryFileInfo,
  type RegisteredInput,
  type StrategyFileInfo,
  type StreamSink,
  type ToolSink,
} from './infra/reactor/abi.js';
export { inspectReactor, type MemoryLimits, type ReactorArtifact, type ReactorKind } from './infra/reactor/artifact.js';
export { engineBuildName, type EngineAssets, type ReactorAsset } from './infra/reactor/assets.js';
export { ReactorInputChecker } from './infra/reactor/checker.js';
export { reactorControl, toProgramDescription, type EngineControl, type ToolControl } from './infra/reactor/control.js';
export { ReactorInflater } from './infra/reactor/inflater.js';
export { compileReactor, instantiateSerial, loadReactor, type CompiledReactor, type ReactorInstance } from './infra/reactor/instance.js';
export { fileReaderSync, RunDirectory, TABLE_BYTES, type BlobReader } from './infra/reactor/run-directory.js';
export { ReactorScanner, reopening } from './infra/reactor/scanner.js';
export { ENGINE_ENDED, ENGINE_NOT_STARTED, WasmEngine, type RuntimeEnd, type WasmEngineOptions } from './infra/engine-worker/gateway.js';
export { DEFAULT_RENEWAL, type RenewalLimits } from './infra/engine-worker/policy.js';
export { RunOutputReceiver, type RunOutputResult } from './infra/run-output/receiver.js';
export { RunOutputWriter } from './infra/run-output/writer.js';
export { holdLife, watchLife, webLocks, type LifeLocks } from './infra/worker-life.js';
