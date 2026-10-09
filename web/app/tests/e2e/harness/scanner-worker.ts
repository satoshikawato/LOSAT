// Runs the RecordScanner contract against the adapter's serial reactor in a dedicated
// worker, with the same scanner module as the Data worker (src/infra/reactor/scanner.ts).
import { ENGINE_ASSETS } from 'virtual:losat-engine';
import { compileReactor, fetchReactor, instantiateSerial } from '../../../src/infra/reactor/instance';
import { ReactorScanner, reopening } from '../../../src/infra/reactor/scanner';
import { runCases, type CaseResult } from '../../contract/contract';
import { RECORD_SCANNER_CASES } from '../../contract/record-scanner.contract';

const scope = self as unknown as { postMessage(message: unknown): void; onmessage: (() => void) | null };

async function run(): Promise<CaseResult[]> {
  if (ENGINE_ASSETS === null) throw new Error('the harness was built without the engine (LOSAT_WEB_REACTORS)');
  const module = (await compileReactor(await fetchReactor(ENGINE_ASSETS.serial))).module;
  const scanner = new ReactorScanner(reopening(async () => (await instantiateSerial(module)).abi));
  return runCases(RECORD_SCANNER_CASES, () => ({ scanner }), {
    onResult: (result) => scope.postMessage({ type: 'progress', result }),
  });
}

scope.onmessage = () => {
  void run().then(
    (results) => scope.postMessage({ type: 'results', results }),
    (error: unknown) => scope.postMessage({ type: 'results', results: [{ name: 'the harness', ok: false, error: String(error) }] }),
  );
};
