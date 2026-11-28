// The tests' `test`. In WebKit, the Screen Wake Lock of the pages of the `context` fixture is a
// stand-in. WebKit's GLib ports (WPE, GTK) connect to the session bus for each lock they take,
// and GLib's setuid check on the way reads errno after getauxval() and aborts the process if it
// is set. WebKit's thread suspend signal handler returns with errno = EINTR, and JavaScriptCore
// suspends the main thread thousands of times while it compiles WebAssembly: when a suspension
// fell inside that check, the web process died with "getauxval () failed: Interrupted system
// call" and the page closed in the middle of its test, more often on a busy machine. The test of
// the wake lock replaces navigator.wakeLock with its own double in every browser; Chromium and
// Firefox keep the real API in the others.
//
// Every test also fails if it leaves a file whose name CI's artifact upload refuses
// (support/artifact-names.ts): one such name made a failed job upload nothing. The fixture renames
// such files first (renamed-names.json keeps the names), so the upload keeps every test's evidence.
import { appendFileSync, mkdirSync } from 'node:fs';
import { join } from 'node:path';
import { test as base } from '@playwright/test';
import { RENAMED_NAMES_FILE, renameUnuploadable } from './artifact-names';

export const test = base.extend<{ uploadableOutputs: void; ww4jLog: void }>({
  ww4jLog: [
    async ({ page }, use, testInfo) => {
      const dir = process.env['WW4J_LOG_DIR'];
      if (!dir) { await use(); return; }
      mkdirSync(dir, { recursive: true });
      const file = join(dir, `${testInfo.project.name}-${testInfo.titlePath.slice(1).join(' ').replace(/[^A-Za-z0-9]+/g, '_').slice(0, 60)}-r${testInfo.repeatEachIndex}-w${testInfo.workerIndex}.log`);
      const seen = new Set<string>();
      const take = (where: string) => (message: { text(): string }) => {
        const text = message.text();
        if (!text.startsWith('[ww4j]') || seen.has(text)) return;
        seen.add(text);
        appendFileSync(file, `${text} (${where})\n`);
      };
      page.on('console', take('page'));
      page.on('worker', (worker) => worker.on('console', take('worker')));
      appendFileSync(file, `[ww4j] ${Date.now()} test-begin\n`);
      await use();
      appendFileSync(file, `[ww4j] ${Date.now()} test-end ${testInfo.status} ${testInfo.error?.message?.split('\n')[0] ?? ''}\n`);
    },
    { auto: true },
  ],
  context: async ({ context, browserName }, use) => {
    if (browserName === 'webkit') await context.addInitScript(standInWakeLock);
    await use(context);
  },
  uploadableOutputs: [
    async ({}, use, testInfo) => {
      await use();
      const names = renameUnuploadable(testInfo.outputDir);
      if (names.length > 0) {
        throw new Error(
          `the test left files whose names CI's artifact upload refuses (" : < > | * ? CR LF), which would make a failed job upload no evidence of any test; ` +
            `renamed (see ${RENAMED_NAMES_FILE}): ${names.join(', ')}`,
        );
      }
    },
    // No time limit of its own (as Playwright's trace fixture): a teardown before it that uses up the
    // test's time must not skip the renaming. A worker that dies skips it: global-teardown.ts.
    { auto: true, timeout: 0 },
  ],
});

/** A wake lock that the browser grants at once and that ends on release(), as WakeLockSentinel does. */
export function standInWakeLock(): void {
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
