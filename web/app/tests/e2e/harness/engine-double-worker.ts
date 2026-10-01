// A test double of the Engine worker: it writes the run outputs that the page asks for
// through RunOutputWriter, from its own thread, to the port of the real Data worker.
import { RunOutputWriter } from '../../../src/infra/run-output/writer';
import type { OutputStream } from '../../../src/ports/run-output';
import type { EngineDoubleCommand, EngineDoubleReply } from './protocol';

const scope = self as unknown as {
  postMessage(message: EngineDoubleReply): void;
  onmessage: ((event: MessageEvent<EngineDoubleCommand>) => void) | null;
};

const writers = new Map<number, { readonly writer: RunOutputWriter; readonly port: MessagePort }>();
let nextWriter = 0;

scope.onmessage = (event) => {
  const command = event.data;
  try {
    if (command.type === 'open') {
      writers.set(++nextWriter, { writer: new RunOutputWriter(command.port), port: command.port });
      scope.postMessage({ id: command.id, ok: true, value: nextWriter });
      return;
    }
    const entry = writers.get(command.writer);
    if (entry === undefined) throw new Error(`unknown writer ${command.writer}`);
    if (command.type === 'write') entry.writer.write(command.stream as OutputStream, command.bytes);
    else if (command.type === 'end') entry.writer.end();
    else entry.port.postMessage(command.message);
    scope.postMessage({ id: command.id, ok: true });
  } catch (error) {
    const failure = error instanceof Error ? error : new Error(String(error));
    scope.postMessage({ id: command.id, ok: false, name: failure.name, message: failure.message });
  }
};
