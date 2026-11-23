// What the Engine worker does with each chunk and frame of a tool run (engine-worker.ts, `tool`): posts it to the
// main thread with its buffer transferred, counting its bytes in the run's credit first and waiting after the post while
// the page has not taken enough (tool-credit.ts). Apart from the worker's scope so that a test can drive it.
import { postedBytes, ToolCreditWorker } from './tool-credit.js';
import type { EngineEvent, WorkerTool } from './protocol.js';
import type { ToolPosts } from './runtime.js';

/**
 * The posts of the tool run `command` of `generation`: `post` sends an event to the main thread. A chunk or a frame of bytes
 * is counted before it is posted (its buffer is transferred, so its length is read first), and the call returns once the
 * page has taken enough for the backlog to be within the limit, or at once when the run has no credit. `limits` are the
 * credit's limit and the slice of its wait: the constants of tool-credit.ts, which a test lowers (the Engine worker gives none).
 */
export function toolPosts(
  command: WorkerTool,
  generation: number,
  post: (event: EngineEvent, transfer: Transferable[]) => void,
  limits: { readonly limit?: number; readonly sliceMs?: number } = {},
): ToolPosts {
  const waiter = command.credit === undefined ? undefined : new ToolCreditWorker(command.credit, limits.limit, limits.sliceMs);
  const { runId } = command;
  return {
    chunk: (stream, bytes) => {
      waiter?.add(bytes.length);
      post({ type: 'tool-chunk', generation, runId, stream, bytes }, [bytes.buffer as ArrayBuffer]);
      waiter?.wait();
    },
    file: (frame) => {
      waiter?.add(postedBytes(frame));
      post({ type: 'tool-file', generation, runId, frame }, frame.kind === 'bytes' ? [frame.bytes.buffer as ArrayBuffer] : []);
      waiter?.wait();
    },
  };
}
