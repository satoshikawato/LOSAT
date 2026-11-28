// One way for the tests to wait for a run to end. `Coordinator.enqueueAll` and `SearchDraft.submit` resolve once the
// run is accepted; `preparing -> running -> finalizing -> completed | cancelled | failed` goes on after that, as
// promise reactions and engine messages, so a test that reads `status` or `record` right after the submit sees
// `finalizing` on a slow machine (D36). This waits on the state store itself: no timer, no retry count, so the
// time it takes is the run's own and a run that never ends is stopped by the test's timeout.
import { isTerminal } from '../../../src/domain/run';
import type { AnyRunView, Coordinator } from '../../../src/application/coordinator';

/** Resolves with the run (a search or a tool run) the first time it is in a final state; at once when it already is. */
export function runEnd(coordinator: Coordinator, runId: string): Promise<AnyRunView> {
  return new Promise((resolve) => {
    const check = () => {
      const run = coordinator.state.get().runs.find((entry) => entry.snapshot.runId === runId);
      if (run === undefined || !isTerminal(run.status)) return;
      unsubscribe();
      resolve(run);
    };
    const unsubscribe = coordinator.state.subscribe(check);
    check();
  });
}
