// Shared Web Lock holder for queued or running runs that read saved WindowMasker counts
// (ports/counts.ts): each run holds `losat-web:windowmasker:<taxid>` in shared mode from
// enqueue until it ends so Replace, Delete and download (which take the same lock in
// exclusive mode with `ifAvailable`) refuse in place across every tab of the origin.
import { countsLockName } from '../../domain/counts';
import type { CountsLocks } from '../../ports/counts';
import type { SessionLocks } from '../data/session';

/** Wraps a `SessionLocks` (`navigator.locks` or a test double) as a `CountsLocks` port. */
export function webCountsLocks(locks: SessionLocks): CountsLocks {
  return {
    holdShared(taxid: number): Promise<(() => void) | undefined> {
      return new Promise((resolve) => {
        locks
          .request(countsLockName(taxid), { mode: 'shared', ifAvailable: true }, (lock) => {
            if (lock === null) {
              resolve(undefined);
              return Promise.resolve();
            }
            let released = false;
            return new Promise<void>((held) => {
              resolve(() => {
                if (released) return;
                released = true;
                held();
              });
            });
          })
          .catch(() => {
            // Where Web Locks are refused (e.g. blocked site storage where OPFS is also
            // unavailable and counts live in this tab only), keep that environment's behavior.
            resolve(() => undefined);
          });
      });
    },
  };
}

/** Browser `navigator.locks` implementation; returns `undefined` where Web Locks are absent. */
export function browserCountsLocks(): CountsLocks | undefined {
  if (typeof navigator === 'undefined' || !('locks' in navigator) || navigator.locks === undefined) {
    return undefined;
  }
  return webCountsLocks(navigator.locks as unknown as SessionLocks);
}
