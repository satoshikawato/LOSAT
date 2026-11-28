// Local measurement spec for WindowMasker `-mk_counts -sformat obinary` (Chromium only,
// outside CI): skipped unless LOSAT_WEB_MEASURE=1. Never downloads anything: reads local
// genome FASTA files from environment variables:
//   - LOSAT_WEB_MEASURE_ECOLI  (E. coli K-12, ~4.6 Mbp; defaults to LOSAT/tests/fasta/MG1655.fna)
//   - LOSAT_WEB_MEASURE_MID    (~120 Mbp genome, e.g. D. melanogaster or A. thaliana)
//   - LOSAT_WEB_MEASURE_HUMAN  (GRCh38 primary assembly, ~3.1 Gbp)
// For each configured input, measures browser wall time, page peak JS heap
// (`performance.memory.usedJSHeapSize`), browser process tree peak RSS (Linux `/proc`), and
// native `losat windowmasker` wall time on the same file, and checks SHA-256 equality of the
// produced `wmasker.obinary`.
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { closeSync, existsSync, mkdirSync, openSync, readdirSync, readFileSync, readSync, statSync } from 'node:fs';
import { basename, join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE, ORIGIN } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { checked, described, openSearch, pickRun, waitStatus } from './support/search';

const MEASURE = process.env['LOSAT_WEB_MEASURE'] === '1';
const MIB = 1024 * 1024;

interface MeasureCase {
  readonly id: string;
  readonly label: string;
  readonly envVar: string;
  readonly defaultPath?: string;
}

const CASES: readonly MeasureCase[] = [
  {
    id: 'ecoli',
    label: 'E. coli K-12 (~4.6 Mbp)',
    envVar: 'LOSAT_WEB_MEASURE_ECOLI',
    defaultPath: join(REPOSITORY, 'LOSAT/tests/fasta/MG1655.fna'),
  },
  {
    id: 'mid',
    label: '~120 Mbp genome',
    envVar: 'LOSAT_WEB_MEASURE_MID',
  },
  {
    id: 'human',
    label: 'Human GRCh38 primary assembly',
    envVar: 'LOSAT_WEB_MEASURE_HUMAN',
  },
];

/** Streams a file in 8 MiB chunks to compute its size and SHA-256 without buffering multi-hundred-MiB outputs. */
function fileStatsAndSha256(filePath: string): { readonly size: number; readonly sha256: string } {
  const hash = createHash('sha256');
  const buf = Buffer.allocUnsafe(8 * MIB);
  const fd = openSync(filePath, 'r');
  let size = 0;
  try {
    for (;;) {
      const n = readSync(fd, buf, 0, buf.length, null);
      if (n === 0) break;
      size += n;
      hash.update(buf.subarray(0, n));
    }
  } finally {
    closeSync(fd);
  }
  return { size, sha256: hash.digest('hex') };
}

/** Sums VmRSS (bytes) across all descendant processes of `rootPid` on Linux via `/proc`. */
function descendantRssBytes(rootPid: number): number | undefined {
  if (process.platform !== 'linux') return undefined;
  try {
    const children = new Map<number, number[]>();
    const rssByPid = new Map<number, number>();
    for (const entry of readdirSync('/proc')) {
      if (!/^\d+$/.test(entry)) continue;
      const pid = Number(entry);
      let status: string;
      try {
        status = readFileSync(`/proc/${entry}/status`, 'utf8');
      } catch {
        continue;
      }
      const ppidMatch = /^PPid:\s*(\d+)/m.exec(status);
      const rssMatch = /^VmRSS:\s*(\d+)\s*kB/m.exec(status);
      if (ppidMatch !== null) {
        const ppid = Number(ppidMatch[1]);
        const list = children.get(ppid);
        if (list === undefined) children.set(ppid, [pid]);
        else list.push(pid);
      }
      if (rssMatch !== null) rssByPid.set(pid, Number(rssMatch[1]) * 1024);
    }
    const queue = [...(children.get(rootPid) ?? [])];
    const visited = new Set<number>();
    let total = 0;
    while (queue.length > 0) {
      const pid = queue.shift()!;
      if (visited.has(pid)) continue;
      visited.add(pid);
      total += rssByPid.get(pid) ?? 0;
      for (const child of children.get(pid) ?? []) queue.push(child);
    }
    return total > 0 ? total : undefined;
  } catch {
    return undefined;
  }
}

/** Reads `performance.memory.usedJSHeapSize` in Chromium when exposed. */
async function readPageHeapBytes(page: Page): Promise<number | undefined> {
  return page.evaluate(() => {
    const perf = performance as Performance & { memory?: { usedJSHeapSize: number } };
    return perf.memory?.usedJSHeapSize;
  });
}

/** Runs native `losat windowmasker -mk_counts -sformat obinary` on `inputPath` and returns elapsed ms + output stats. */
function runNativeMkCounts(
  inputPath: string,
  outDir: string,
): { readonly elapsedMs: number; readonly outBytes: number; readonly sha256: string } {
  if (NATIVE === undefined) throw new Error(NO_NATIVE_REASON);
  mkdirSync(outDir, { recursive: true });
  const outPath = join(outDir, 'wmasker.obinary');
  const t0 = performance.now();
  const res = spawnSync(
    NATIVE,
    ['windowmasker', '-mk_counts', '-sformat', 'obinary', '-in', inputPath, '-out', outPath],
    { maxBuffer: 16 * MIB },
  );
  const elapsedMs = performance.now() - t0;
  if (res.error) throw res.error;
  if (res.status !== 0) {
    throw new Error(`native windowmasker failed (${res.status}): ${res.stderr.toString('utf8')}`);
  }
  const { size, sha256 } = fileStatsAndSha256(outPath);
  return { elapsedMs, outBytes: size, sha256 };
}

test.describe('WindowMasker -mk_counts -sformat obinary local measurements', () => {
  test.skip(!MEASURE, 'local measurement spec: set LOSAT_WEB_MEASURE=1');
  test.skip(!BUILD_HAS_ENGINE, 'needs the engine build (LOSAT_WEB_REACTORS)');
  test.skip(NATIVE === undefined, NO_NATIVE_REASON);

  for (const item of CASES) {
    test(`measure ${item.label}`, async ({ page, browserName }, testInfo) => {
      test.skip(browserName !== 'chromium', 'Chromium only');
      test.setTimeout(3_600_000);

      const inputPath = process.env[item.envVar] || item.defaultPath;
      test.skip(
        inputPath === undefined || !existsSync(inputPath),
        `set ${item.envVar} to a local FASTA file to run this measurement`,
      );

      const fileSize = statSync(inputPath!).size;
      const external: string[] = [];
      page.on('request', (req) => {
        const url = req.url();
        if (!url.startsWith(ORIGIN) && !url.startsWith('blob:') && !url.startsWith('data:')) {
          external.push(url);
        }
      });

      // Measure native first so the output SHA-256 is ready as the reference.
      const native = runNativeMkCounts(inputPath!, testInfo.outputPath(`native-${item.id}`));

      await openSearch(page);
      await page.getByTestId('tab-windowmasker').click();
      await described(page, 'wm-panel', 'ready');
      await page.getByTestId('wm-stage-mk_counts').check();
      await page.getByTestId('wm-param-sformat').selectOption('obinary');
      await expect(page.getByTestId('wm-out')).toHaveValue('wmasker.obinary');

      const heapBefore = await readPageHeapBytes(page);
      const rssBefore = descendantRssBytes(process.pid);

      await page.getByTestId('wm-fasta-files').setInputFiles(inputPath!);
      await expect(page.locator('[data-testid="wm-fasta-list"] > li').first()).toHaveAttribute('data-status', 'ready', {
        timeout: 300_000,
      });
      await checked(page, 'wm-panel');

      // Sample page JS heap inside the page every 250 ms during execution, and browser RSS from Node.
      await page.evaluate(() => {
        const win = window as Window & { __wmPeakHeap?: number; __wmHeapTimer?: number };
        const perf = performance as Performance & { memory?: { usedJSHeapSize: number } };
        win.__wmPeakHeap = perf.memory?.usedJSHeapSize ?? 0;
        win.__wmHeapTimer = window.setInterval(() => {
          const cur = perf.memory?.usedJSHeapSize ?? 0;
          if (cur > (win.__wmPeakHeap ?? 0)) win.__wmPeakHeap = cur;
        }, 250);
      });

      let peakRss = rssBefore ?? 0;
      const rssTimer = setInterval(() => {
        const cur = descendantRssBytes(process.pid);
        if (cur !== undefined && cur > peakRss) peakRss = cur;
      }, 250);

      const t0 = performance.now();
      await page.getByTestId('wm-run-losat').click();
      await waitStatus(page, 1, 'completed', 3_000_000);
      const webElapsedMs = performance.now() - t0;
      clearInterval(rssTimer);

      const peakPageHeap = await page.evaluate(() => {
        const win = window as Window & { __wmPeakHeap?: number; __wmHeapTimer?: number };
        if (win.__wmHeapTimer !== undefined) window.clearInterval(win.__wmHeapTimer);
        const perf = performance as Performance & { memory?: { usedJSHeapSize: number } };
        const cur = perf.memory?.usedJSHeapSize ?? 0;
        return Math.max(win.__wmPeakHeap ?? 0, cur);
      });
      const heapAfter = await readPageHeapBytes(page);

      // Check SHA-256 matches native byte-for-byte via Download.
      await pickRun(page, 1);
      await page.getByTestId('results-view-outputs').click();
      await page.getByTestId('tool-tab-file-0').click();
      const [download] = await Promise.all([
        page.waitForEvent('download'),
        page.getByTestId('tool-output-download').click(),
      ]);
      const downloaded = fileStatsAndSha256((await download.path())!);
      expect(downloaded.size).toBe(native.outBytes);
      expect(downloaded.sha256).toBe(native.sha256);
      expect(external).toEqual([]);

      const fmtMiB = (bytes: number | undefined) =>
        bytes === undefined ? 'n/a' : `${(bytes / MIB).toFixed(1)} MiB`;
      const summary = {
        id: item.id,
        file: basename(inputPath!),
        inputBytes: fileSize,
        outputBytes: native.outBytes,
        sha256: native.sha256,
        webMs: Math.round(webElapsedMs),
        nativeMs: Math.round(native.elapsedMs),
        ratio: Number((webElapsedMs / Math.max(1, native.elapsedMs)).toFixed(2)),
        pageHeapBeforeMiB: fmtMiB(heapBefore),
        pageHeapPeakMiB: fmtMiB(peakPageHeap > 0 ? peakPageHeap : undefined),
        pageHeapAfterMiB: fmtMiB(heapAfter),
        browserRssBeforeMiB: fmtMiB(rssBefore),
        browserRssPeakMiB: fmtMiB(peakRss > 0 ? peakRss : undefined),
      };
      console.log(`WM-MEASURE ${JSON.stringify(summary)}`);
    });
  }
});
