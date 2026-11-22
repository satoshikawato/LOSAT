// A search strategy file (kind 3 of `file_begin`, docs/web/abi_v2.md §4, §7; line RG item F2-4) run
// by the serial reactor in a worker of its own, with the modules of the Engine worker: the file
// goes in through `loadFile`, the argv names it and no `-query` or `-subject`, and the run takes
// the strategy's own queries and subjects through `register_inputs`.
import { ENGINE_ASSETS } from 'virtual:losat-engine';
import { FILE_SEARCH_STRATEGY, instantiateSerial, loadFile, loadReactor, sha256Hex } from 'losat';
import type { StrategySearch, StrategySearchEvent, StrategySearchResult } from './protocol';

const scope = self as unknown as {
  postMessage(message: StrategySearchEvent): void;
  onmessage: ((event: MessageEvent<StrategySearch>) => void) | null;
};

async function run(search: StrategySearch): Promise<StrategySearchResult> {
  if (ENGINE_ASSETS === null) throw new Error('the harness was built without the engine (LOSAT_WEB_REACTORS)');
  const { module } = await loadReactor(ENGINE_ASSETS.serial);
  const { abi } = await instantiateSerial(module);
  const file = await loadFile(abi, FILE_SEARCH_STRATEGY, search.name, new Blob([new Uint8Array(search.strategy)]));
  const argv = [search.program, '-import_search_strategy', search.name, ...search.options, '-num_threads', '1'];
  const started = abi.start(argv);
  if (!started.ok) return { kind: file.kind, argv, error: started.message };
  const registerInputs = started.start.registerInputs;
  let together;
  try {
    together = abi.registerInputs(argv, new Uint8Array(0), new Uint8Array(0));
  } catch (error) {
    return { kind: file.kind, argv, registerInputs, error: error instanceof Error ? error.message : String(error) };
  }
  const streams = new Map<number, Uint8Array[]>();
  try {
    abi.run(argv, together.query.handle, together.subject.handle, (stream, bytes) =>
      streams.set(stream, [...(streams.get(stream) ?? []), bytes.slice()]),
    );
  } finally {
    abi.release(together.query.handle);
    abi.release(together.subject.handle);
  }
  const sha256: Record<number, string> = {};
  const bytes: Record<number, number> = {};
  for (const [stream, chunks] of streams) {
    const joined = new Uint8Array(chunks.reduce((total, chunk) => total + chunk.length, 0));
    let offset = 0;
    for (const chunk of chunks) {
      joined.set(chunk, offset);
      offset += chunk.length;
    }
    sha256[stream] = await sha256Hex(joined);
    bytes[stream] = joined.length;
  }
  return {
    kind: file.kind,
    argv,
    registerInputs,
    queryRecords: together.query.records.length,
    subjectRecords: together.subject.records.length,
    sha256,
    bytes,
  };
}

scope.onmessage = (event) => {
  void run(event.data).then(
    (result) => scope.postMessage({ type: 'result', result }),
    (error: unknown) => scope.postMessage({ type: 'error', message: error instanceof Error ? error.message : String(error) }),
  );
};
