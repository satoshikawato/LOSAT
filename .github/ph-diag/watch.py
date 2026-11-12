"""DIAGNOSTIC ONLY (branch diag/ph-webkit-hang, line WFL item PH; not for merge).

Dumps WPEWebProcess processes of the runner in three cases, at most once per process:
- a failed WebKit test asks for it (tests/e2e/support/test.ts writes <out>/fail-<ms>-<worker pid>.txt and waits
  for <out>/done-<same id>): every web process under that Playwright worker, while its pages are still open;
- the MAIN thread (tid == pid) sleeps anywhere but the GLib poll of its run loop, with no CPU time added, for
  BLOCK_S seconds (default 20);
- the main thread runs (state R) and adds CPU time without a pause for BUSY_S seconds (default 30).
A dump (<out>/<kind>_<pid>_<n>.txt): every thread (name, state, wchan, CPU time, SigBlk/SigPnd/ShdPnd/SigCgt,
kernel stack), a gdb backtrace of every thread, the registers of each thread, the JITWorklist (jitdump.py), maps.
Run as root (ptrace of non-children). usage: watch.py <out dir>"""
import glob, os, subprocess, sys, time

out = sys.argv[1]
BLOCK_S = float(os.environ.get('BLOCK_S', '20'))
BUSY_S = float(os.environ.get('BUSY_S', '30'))
HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(out, exist_ok=True)
log = open(os.path.join(out, 'watch.log'), 'a', buffering=1)


def read(path):
    try:
        return open(path).read()
    except OSError as e:
        return f'({e})'


def procs():
    for p in os.listdir('/proc'):
        if p.isdigit() and read(f'/proc/{p}/comm').startswith('WPEWebProcess'):
            yield p


def ppid(p):
    stat = read(f'/proc/{p}/stat')
    return stat[stat.rindex(')') + 2:].split()[1] if ')' in stat else '0'


def under(p, ancestor):
    for _ in range(12):
        p = ppid(p)
        if p == ancestor:
            return True
        if p in ('0', '1'):
            return False
    return False


def main_state(p):
    stat = open(f'/proc/{p}/task/{p}/stat').read()
    rest = stat[stat.rindex(')') + 2:].split()
    return rest[0], int(rest[11]) + int(rest[12]), open(f'/proc/{p}/task/{p}/wchan').read().strip()


episodes = 0
dumped = set()


def dump(p, kind, why):
    global episodes
    if p in dumped:
        return
    dumped.add(p)
    episodes += 1
    lines = [f'{kind} {episodes} at {time.strftime("%F %T")}: {why}', f'cmdline {read(f"/proc/{p}/cmdline")!r}']
    for t in sorted(os.listdir(f'/proc/{p}/task'), key=int) if os.path.isdir(f'/proc/{p}/task') else []:
        stat = read(f'/proc/{p}/task/{t}/stat')
        rest = stat[stat.rindex(')') + 2:].split() if ')' in stat else ['?'] * 20
        status = read(f'/proc/{p}/task/{t}/status')
        sig = ' '.join(l.replace('\t', '') for l in status.splitlines() if l.startswith(('SigBlk', 'SigPnd', 'ShdPnd', 'SigCgt')))
        lines.append(f'{t} {read(f"/proc/{p}/task/{t}/comm").strip()!r} {rest[0]} wchan={read(f"/proc/{p}/task/{t}/wchan").strip()} '
                     f'utime={rest[11]} stime={rest[12]} {sig}')
        lines.append('   kstack: ' + ' | '.join(read(f'/proc/{p}/task/{t}/stack').split('\n')[:8]))
    try:
        bt = subprocess.run(['gdb', '-p', p, '-batch', '-ex', 'set pagination off', '-ex', 'thread apply all bt 60',
                             '-ex', 'thread apply all info registers rip rsp rbp rdi rsi rdx rcx r8 r12 r13 r14 r15',
                             '-ex', f'source {HERE}/jitdump.py'], capture_output=True, text=True, timeout=300).stdout
    except Exception as e:  # noqa: BLE001
        bt = f'gdb failed: {e}'
    path = os.path.join(out, f'{kind}_{p}_{episodes}.txt')
    open(path, 'w').write('\n'.join(lines) + '\n\n' + bt + '\n\nmaps\n' + read(f'/proc/{p}/maps'))
    log.write(f'{time.strftime("%F %T")} DUMP {kind} pid {p}: {why} -> {path}\n')


seen = {}  # pid -> [since, cpu, kind]
handled = set()
while True:
    # Requests of failed tests first: their pages close once the request is answered.
    for f in sorted(glob.glob(os.path.join(out, 'fail-*.txt'))):
        rid = os.path.basename(f)[len('fail-'):-len('.txt')]
        if rid in handled:
            continue
        handled.add(rid)
        worker = rid.rsplit('-', 1)[-1]
        title = read(f).splitlines()[0] if read(f) else ''
        targets = [p for p in procs() if under(p, worker)]
        log.write(f'{time.strftime("%F %T")} request {rid} ({title}): web processes {targets}\n')
        for p in targets:
            dump(p, 'fail', f'test failed: {title}')
        open(os.path.join(out, f'done-{rid}'), 'w').write('done\n')
    now = time.time()
    alive = set()
    for p in procs():
        alive.add(p)
        try:
            state, cpu, wchan = main_state(p)
        except (OSError, ValueError):
            continue
        blocked = state in 'SD' and not any(w in wchan for w in ('poll', 'select'))
        kind = 'block' if blocked else ('busy' if state == 'R' else 'idle')
        s = seen.get(p)
        moved = s is not None and cpu != s[1]
        if s is None or kind != s[2] or kind == 'idle' or (kind == 'block' and moved) or (kind == 'busy' and not moved):
            seen[p] = [now, cpu, kind]
            continue
        s[1] = cpu
        if now - s[0] >= (BLOCK_S if kind == 'block' else BUSY_S):
            dump(p, kind, f'main thread {kind} {now - s[0]:.0f} s (state {state}, wchan={wchan})')
    for p in list(seen):
        if p not in alive:
            del seen[p]
    time.sleep(1)
