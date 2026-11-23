// A worker thread for tests/tool-credit.test.ts: it runs the Engine worker's `toolPosts` (tool-posts.ts) as the engine's call
// into the tool does, with a posted chunk or file frame for each step, and reports each step to the main thread, so that the
// test sees where the call blocks. Run by Node itself (type stripping, with register-ts.mjs for the `.js` names of the
// sources), so it uses no syntax that needs a transform.
import { parentPort, workerData } from 'node:worker_threads';
import { toolPosts } from '../../src/infra/engine-worker/tool-posts.ts';

interface Step {
  readonly kind: 'chunk' | 'file';
  readonly size: number;
}

const { buffer, limit, sliceMs, steps } = workerData as { buffer: SharedArrayBuffer; limit?: number; sliceMs?: number; steps: Step[] };
let posted = 0;
// The command is the one the Gateway builds for a run with a credit (gateway.ts, launchTool).
const posts = toolPosts(
  { type: 'tool', generation: 1, runId: 'r', argv: ['windowmasker'], credit: buffer },
  1,
  () => parentPort!.postMessage(`posted ${++posted}`),
  { ...(limit === undefined ? {} : { limit }), ...(sliceMs === undefined ? {} : { sliceMs }) },
);
steps.forEach((step, index) => {
  if (step.kind === 'chunk') posts.chunk('stdout', new Uint8Array(step.size));
  else posts.file({ number: 0, kind: 'bytes', offset: 0, bytes: new Uint8Array(step.size) });
  parentPort!.postMessage(`after ${index + 1}`);
});
parentPort!.postMessage('done');
