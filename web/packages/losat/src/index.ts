// The npm package `losat`: the host of LOSAT's two Wasm reactors (docs/web/abi_v2.md) for browsers.
// LOSAT Web (web/app) runs the engine through these modules. The layers follow web/AGENTS.md
// rule 5: `domain` depends on nothing, `ports` on the domain, `infra` implements the ports.
export * from './domain';
export * from './ports';
export { sha256Hex } from './infra/digest';
export {
  EngineAllocationError,
  EngineCallError,
  EngineStoppedError,
  ReactorAbi,
  type AbiExports,
  type BegunInflation,
  type ChunkBuffer,
  type FileInfo,
  type RegisteredInput,
  type StreamSink,
} from './infra/reactor/abi';
export { inspectReactor, type MemoryLimits, type ReactorArtifact, type ReactorKind } from './infra/reactor/artifact';
export { engineBuildName, type EngineAssets, type ReactorAsset } from './infra/reactor/assets';
export { ReactorInputChecker } from './infra/reactor/checker';
export { reactorControl, toProgramDescription, type EngineControl } from './infra/reactor/control';
export { ReactorInflater } from './infra/reactor/inflater';
export { compileReactor, instantiateSerial, loadReactor, type CompiledReactor, type ReactorInstance } from './infra/reactor/instance';
export { ReactorScanner, reopening } from './infra/reactor/scanner';
export { WasmEngine, type RuntimeEnd, type WasmEngineOptions } from './infra/engine-worker/gateway';
export { DEFAULT_RENEWAL, type RenewalLimits } from './infra/engine-worker/policy';
export { RunOutputReceiver, type RunOutputResult } from './infra/run-output/receiver';
export { RunOutputWriter } from './infra/run-output/writer';
