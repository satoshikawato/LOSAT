// A profile on disk of the test's browser, for the storage tests (Chromium) and the OPFS
// BlockStore contract (browsers.spec.ts). Playwright's default contexts are incognito: in
// Chromium their OPFS lives in memory and ignores the quota that CDP sets, unlike a user's
// browser profile, and WebKit 27 refuses OPFS in them ("UnknownError" from getDirectory). A
// persistent context keeps OPFS on disk, as users have it. It starts with the project's launch
// options (the WebKit launcher and environment of playwright.config.ts).
import { readdirSync, readFileSync } from 'node:fs';
import type { BrowserContext } from '@playwright/test';
import { ORIGIN } from './browser';
import { test as base } from './test';

export interface Profile {
  readonly context: BrowserContext;
  /** The profile's user data directory. */
  readonly dir: string;
}

export const test = base.extend<{ profile: Profile }>({
  profile: async ({ playwright, browserName }, use, testInfo) => {
    const dir = testInfo.outputPath('profile');
    const context = await playwright[browserName].launchPersistentContext(dir, {
      baseURL: ORIGIN,
      viewport: { width: 1280, height: 720 },
      acceptDownloads: true,
    });
    await use({ context, dir });
    await context.close();
  },
});

export { expect } from '@playwright/test';

interface ProcessEntry {
  readonly pid: number;
  readonly parent: number;
  readonly command: string;
}

function processes(): ProcessEntry[] {
  const entries: ProcessEntry[] = [];
  for (const name of readdirSync('/proc')) {
    if (!/^\d+$/.test(name)) continue;
    try {
      const command = readFileSync(`/proc/${name}/cmdline`, 'utf8').split('\0').join(' ');
      const stat = readFileSync(`/proc/${name}/stat`, 'utf8');
      const parent = Number(stat.slice(stat.lastIndexOf(')') + 2).split(' ')[1]);
      entries.push({ pid: Number(name), parent, command });
    } catch {
      // The process ended while it was read.
    }
  }
  return entries;
}

/**
 * Ends every renderer process of the profile with SIGKILL, as the operating system ends a
 * tab that it has to stop. Linux only. Every open page of the profile is ended.
 */
export function killRenderers(profile: Profile): number {
  const all = processes();
  const browser = all.find((p) => p.command.includes(`--user-data-dir=${profile.dir}`) && !p.command.includes('--type='));
  if (browser === undefined) throw new Error('the browser process of the profile was not found');
  const family = new Set([browser.pid]);
  for (let grew = true; grew; ) {
    grew = false;
    for (const p of all) {
      if (!family.has(p.pid) && family.has(p.parent)) {
        family.add(p.pid);
        grew = true;
      }
    }
  }
  const renderers = all.filter((p) => family.has(p.pid) && p.command.includes('--type=renderer'));
  let killed = 0;
  for (const renderer of renderers) {
    try {
      process.kill(renderer.pid, 'SIGKILL');
      killed += 1;
    } catch (error) {
      // The renderer ended after the list was read (a tab closed just before ends its
      // renderer a moment later): it is already gone.
      if ((error as NodeJS.ErrnoException).code !== 'ESRCH') throw error;
    }
  }
  return killed;
}
