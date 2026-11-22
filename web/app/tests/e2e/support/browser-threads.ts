// What the threads of the test's browser do, for a test whose run did not end in time. On CI (dev 35864267c,
// 2026-11-20) the second search of unavailable.spec.ts's "the next search runs" stayed "Searching" in WebKit for the
// whole 120 s; the test takes 5 to 22 s, and it never stalled in 223 WebKit, 20 Chromium and 20 Firefox runs on a
// 32-core machine, also on 2 or 4 cores and with dev's engine (line WQA, losat-baselines/web-qa-20261119/fix-d8/ci).
// A failure that says only "running" cannot tell a slow search (threads that run) from a stalled one (threads that
// wait, and where). The browser is a child of the test's worker process, so its threads are read from /proc (Linux
// only; elsewhere nothing is read).
import { readdirSync, readFileSync } from 'node:fs';

interface Thread {
  readonly name: string;
  /** The kernel's state of the thread (R running, S sleeping, D waiting on I/O, T stopped, ...). */
  readonly state: string;
  /** The kernel function where the thread waits (wchan), or '' where it runs. */
  readonly wchan: string;
  /** User and system CPU time, in clock ticks. */
  readonly ticks: number;
}

/** Each thread of the browser's processes, by "pid/tid". */
function browserThreads(): Map<string, Thread> {
  const parents = new Map<number, number>();
  for (const entry of readdirSync('/proc')) {
    if (!/^\d+$/.test(entry)) continue;
    try {
      const stat = readFileSync(`/proc/${entry}/stat`, 'utf8');
      parents.set(Number(entry), Number(stat.slice(stat.lastIndexOf(')') + 2).split(' ')[1]));
    } catch {
      // The process ended meanwhile.
    }
  }
  // The processes that this worker process started, and theirs: the browser's.
  const ours = new Set([process.pid]);
  for (let grew = true; grew; ) {
    grew = false;
    for (const [pid, parent] of parents) {
      if (!ours.has(pid) && ours.has(parent)) {
        ours.add(pid);
        grew = true;
      }
    }
  }
  ours.delete(process.pid);
  const threads = new Map<string, Thread>();
  for (const pid of ours) {
    let tids: string[];
    try {
      tids = readdirSync(`/proc/${pid}/task`);
    } catch {
      continue;
    }
    for (const tid of tids) {
      try {
        const stat = readFileSync(`/proc/${pid}/task/${tid}/stat`, 'utf8');
        const fields = stat.slice(stat.lastIndexOf(')') + 2).split(' ');
        let wchan = '';
        try {
          wchan = readFileSync(`/proc/${pid}/task/${tid}/wchan`, 'utf8').trim();
        } catch {
          // Not readable here.
        }
        threads.set(`${pid}/${tid}`, {
          name: stat.slice(stat.indexOf('(') + 1, stat.lastIndexOf(')')),
          state: fields[0] ?? '?',
          wchan: wchan === '0' ? '' : wchan,
          ticks: Number(fields[11]) + Number(fields[12]),
        });
      } catch {
        // The thread ended meanwhile.
      }
    }
  }
  return threads;
}

/** The threads of the workers: Chromium's DedicatedWorker, Firefox's DOM Worker, WebKit's "WebCore: Worker". */
const WORKER_THREAD = /DedicatedWorker|DOM Worker|WebCore: Worker/;

/**
 * The threads of the test's browser over `ms`: those that used more than a tenth of a CPU, with their share, the
 * workers' threads with their state and where they wait, and how many others there are. Undefined where /proc cannot
 * be read.
 */
export async function describeBrowserThreads(ms = 2_000): Promise<string | undefined> {
  let before: Map<string, Thread>;
  try {
    before = browserThreads();
  } catch {
    return undefined;
  }
  await new Promise((resolve) => setTimeout(resolve, ms));
  const after = browserThreads();
  const ticksPerSecond = 100;
  const where = (thread: Thread) => `${thread.state}${thread.wchan ? ` ${thread.wchan}` : ''}`;
  const busy: string[] = [];
  const workers = new Map<string, number>();
  let others = 0;
  let total = 0;
  for (const [key, thread] of after) {
    const share = (thread.ticks - (before.get(key)?.ticks ?? thread.ticks)) / ticksPerSecond / (ms / 1000);
    total += share;
    if (share > 0.1) busy.push(`${thread.name} ${share.toFixed(2)} (${where(thread)})`);
    else if (WORKER_THREAD.test(thread.name)) {
      const kind = `${thread.name} ${where(thread)}`;
      workers.set(kind, (workers.get(kind) ?? 0) + 1);
    } else others++;
  }
  const idleWorkers = [...workers].map(([kind, count]) => `${count} × ${kind}`);
  return (
    `The browser's ${after.size} threads over ${ms / 1000} s used ${total.toFixed(2)} CPUs. ` +
    `Busy: ${busy.length === 0 ? 'none' : busy.join(', ')}. ` +
    `Idle worker threads: ${idleWorkers.length === 0 ? 'none' : idleWorkers.join(', ')}. ${others} other idle threads.`
  );
}
