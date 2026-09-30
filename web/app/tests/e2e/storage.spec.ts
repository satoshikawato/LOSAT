// Temporary storage in real browser tabs (plan §5.6, S10): the layout under tmp/, recovery
// after a tab is killed or closed, protection of other open and stopped tabs, and a storage
// quota that runs out. Each test uses its own Chromium profile on disk (support/profile.ts)
// and the Chrome DevTools Protocol to freeze a tab and to set the quota.
import { readFile } from 'node:fs/promises';
import {
  expectResult,
  heldSessionLocks,
  openApp,
  openProbe,
  queueSearch,
  runSearch,
  sessionDirectories,
  sessionFiles,
  setQuota,
} from './support/browser';
import { expect, killRenderers, test } from './support/profile';

const UUID = /[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}/g;

test.skip(({ browserName }) => browserName !== 'chromium', 'uses a Chromium profile and the Chrome DevTools Protocol');

test('a tab keeps its outputs in OPFS under tmp/<session>/runs/<run>/, named by tokens only', async ({ profile }) => {
  const page = await openApp(profile.context);
  await expect(page.getByTestId('storage-status')).toHaveAttribute('data-backend', 'opfs');
  await runSearch(page, 1);
  await expect(page.getByTestId('storage-usage')).toContainText(/^Browser storage \(OPFS\): [1-9]\d* bytes used by this tab/);

  const probe = await openProbe(profile.context);
  const sessions = await sessionDirectories(probe);
  expect(sessions).toHaveLength(1);
  expect(sessions[0]).toMatch(UUID);
  expect(await heldSessionLocks(probe)).toEqual(sessions);
  const files = await sessionFiles(probe, sessions[0]!);
  expect(files.map((path) => path.replace(UUID, '<run>'))).toEqual([
    'runs/<run>/diagnostics',
    'runs/<run>/hits',
    'runs/<run>/out0',
    'runs/<run>/out6',
    'runs/<run>/out7',
  ]);
  expect(files.join('\n')).not.toMatch(/query|subject|\.fa/);
});

for (const ending of ['is killed', 'is closed'] as const) {
  test(`the next tab removes the data of a tab that ${ending}, and keeps its own`, async ({ profile }) => {
    test.skip(ending === 'is killed' && process.platform !== 'linux', 'finds the renderer processes in /proc');
    const first = await openApp(profile.context);
    await runSearch(first, 1);
    const before = await openProbe(profile.context);
    const [abandoned] = await sessionDirectories(before);
    await before.close();

    if (ending === 'is killed') {
      // A forced termination: the tab's renderer process (with its Data worker) gets SIGKILL.
      const crashed = first.waitForEvent('crash');
      expect(killRenderers(profile)).toBeGreaterThan(0);
      await crashed;
    } else {
      await first.close();
    }

    const probe = await openProbe(profile.context);
    // The tab ended without cleaning up: its data is still there, but nobody owns it.
    await expect.poll(() => heldSessionLocks(probe)).toEqual([]);
    expect(await sessionDirectories(probe)).toEqual([abandoned]);

    const next = await openApp(profile.context);
    await expect(next.getByTestId('storage-status')).toHaveAttribute('data-removed-sessions', '1');
    await expect(next.getByTestId('storage-status')).toContainText('Removed the temporary data left by 1 closed tab.');
    const sessions = await sessionDirectories(probe);
    expect(sessions).toHaveLength(1);
    expect(sessions).not.toContain(abandoned);
    await runSearch(next, 1);
    await expectResult(next, 1);
  });
}

test("two open tabs never remove each other's data", async ({ profile }) => {
  const a = await openApp(profile.context);
  await runSearch(a, 1);
  const probe = await openProbe(profile.context);
  const [tokenA] = await sessionDirectories(probe);

  const b = await openApp(profile.context);
  await expect(b.getByTestId('storage-status')).toHaveAttribute('data-removed-sessions', '0');
  await runSearch(b, 1);
  const both = await sessionDirectories(probe);
  expect(both).toHaveLength(2);
  expect(await heldSessionLocks(probe)).toEqual(both);
  const tokenB = both.find((token) => token !== tokenA)!;

  // B closes; a third tab removes B's data and keeps A's, which A can still read.
  await b.close();
  await expect.poll(() => heldSessionLocks(probe)).toEqual([tokenA]);
  const c = await openApp(profile.context);
  await expect(c.getByTestId('storage-status')).toHaveAttribute('data-removed-sessions', '1');
  const after = await sessionDirectories(probe);
  expect(after).toContain(tokenA);
  expect(after).not.toContain(tokenB);
  expect(after).toHaveLength(2);
  await expectResult(a, 1);
  await runSearch(a, 2);
});

test('a stopped (frozen) tab keeps its data while another tab starts', async ({ profile }) => {
  const a = await openApp(profile.context);
  await runSearch(a, 1);
  const probe = await openProbe(profile.context);
  const [tokenA] = await sessionDirectories(probe);

  const cdp = await profile.context.newCDPSession(a);
  await cdp.send('Page.setWebLifecycleState', { state: 'frozen' });
  try {
    const c = await openApp(profile.context);
    await expect(c.getByTestId('storage-status')).toHaveAttribute('data-removed-sessions', '0');
    expect(await sessionDirectories(probe)).toContain(tokenA);
    expect(await heldSessionLocks(probe)).toContain(tokenA);
  } finally {
    await cdp.send('Page.setWebLifecycleState', { state: 'active' });
  }
  await expectResult(a, 1);
});

test('a run that runs out of storage fails with the reason; earlier results stay and later runs work', async ({
  profile,
}) => {
  const page = await openApp(profile.context);
  await runSearch(page, 1);
  // A quota of zero, not the current usage: navigator.storage.estimate() can lag behind the
  // usage that Chromium checks, which would leave room for the next run.
  await setQuota(page, 0);
  try {
    await queueSearch(page);
    await expect(page.getByTestId('run-2-status')).toHaveText('failed');
    await expect(page.getByTestId('run-2')).toContainText('Not enough temporary storage for the results of this run');

    await expectResult(page, 1);
    const download = page.waitForEvent('download');
    await page.getByTestId('export-output').click();
    const saved = await readFile((await (await download).path())!, 'utf8');
    expect(saved).toContain('FAKE ENGINE OUTPUT');
  } finally {
    await setQuota(page, null);
  }
  await runSearch(page, 3);
  await expectResult(page, 3);
});
