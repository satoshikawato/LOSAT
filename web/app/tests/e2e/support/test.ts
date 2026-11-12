// The tests' `test`. In WebKit, the Screen Wake Lock of the pages of the `context` fixture is a
// stand-in. WebKit's GLib ports (WPE, GTK) connect to the session bus for each lock they take,
// and GLib's setuid check on the way reads errno after getauxval() and aborts the process if it
// is set. WebKit's thread suspend signal handler returns with errno = EINTR, and JavaScriptCore
// suspends the main thread thousands of times while it compiles WebAssembly: when a suspension
// fell inside that check, the web process died with "getauxval () failed: Interrupted system
// call" and the page closed in the middle of its test, more often on a busy machine. The test of
// the wake lock replaces navigator.wakeLock with its own double in every browser; Chromium and
// Firefox keep the real API in the others.
import { existsSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { test as base } from '@playwright/test';

export const test = base.extend({
  context: async ({ context, browserName }, use, testInfo) => {
    if (browserName === 'webkit') await context.addInitScript(standInWakeLock);
    await use(context);
    // DIAGNOSTIC ONLY (branch diag/ph-webkit-hang, not for merge): when a WebKit test fails, ask the
    // root watcher (.github/ph-diag/watch.py) to dump every web process while its pages are still open.
    const dir = process.env['PH_DIAG_DIR'];
    if (dir === undefined || browserName !== 'webkit' || testInfo.status === testInfo.expectedStatus) return;
    const id = `${Date.now()}-${process.pid}`;
    writeFileSync(join(dir, `fail-${id}.txt`), `${testInfo.titlePath.join(' > ')}\n${testInfo.error?.message ?? ''}\n`);
    for (let i = 0; i < 240 && !existsSync(join(dir, `done-${id}`)); i++) await new Promise((r) => setTimeout(r, 1000));
  },
});

/** A wake lock that the browser grants at once and that ends on release(), as WakeLockSentinel does. */
function standInWakeLock(): void {
  type Sentinel = EventTarget & { released: boolean; readonly type: 'screen'; release(): Promise<void> };
  Object.defineProperty(navigator, 'wakeLock', {
    configurable: true,
    value: {
      request: async () => {
        const sentinel = Object.assign(new EventTarget(), { released: false, type: 'screen' as const }) as Sentinel;
        sentinel.release = async () => {
          if (sentinel.released) return;
          sentinel.released = true;
          sentinel.dispatchEvent(new Event('release'));
        };
        return sentinel;
      },
    },
  });
}
