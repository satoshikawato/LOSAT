"""DIAGNOSTIC ONLY (branch diag/ph-webkit-hang, line WFL item PH; not for merge).

Watches every WPEWebProcess on the runner. When a process's MAIN thread (tid == pid) sleeps in a futex with no
CPU time added for BLOCK_S seconds (default 20), writes <out>/stuck_<pid>_<n>.txt once per episode: every
thread (name, state, wchan, CPU time, SigBlk/SigPnd/ShdPnd/SigCgt from status, kernel stack), a gdb backtrace of
every thread with the registers of each thread's newest frame, the JITWorklist (jitdump.py), and the maps.
Run as root (ptrace of non-children). usage: watch.py <out dir>"""
import os, subprocess, sys, time

out = sys.argv[1]
BLOCK_S = float(os.environ.get('BLOCK_S', '20'))
HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(out, exist_ok=True)
log = open(os.path.join(out, 'watch.log'), 'a', buffering=1)


def procs():
    for p in os.listdir('/proc'):
        if p.isdigit():
            try:
                if open(f'/proc/{p}/comm').read().startswith('WPEWebProcess'):
                    yield p
            except OSError:
                pass


def main_state(p):
    stat = open(f'/proc/{p}/task/{p}/stat').read()
    rest = stat[stat.rindex(')') + 2:].split()
    return rest[0], int(rest[11]) + int(rest[12]), open(f'/proc/{p}/task/{p}/wchan').read().strip()


def read(path):
    try:
        return open(path).read()
    except OSError as e:
        return f'({e})'


def dump(p, n, why):
    lines = [f'episode {n} at {time.strftime("%F %T")}: {why}', f'cmdline {read(f"/proc/{p}/cmdline")!r}']
    for t in sorted(os.listdir(f'/proc/{p}/task'), key=int):
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
                             '-ex', f'source {HERE}/jitdump.py'], capture_output=True, text=True, timeout=600).stdout
    except Exception as e:  # noqa: BLE001
        bt = f'gdb failed: {e}'
    path = os.path.join(out, f'stuck_{p}_{n}.txt')
    open(path, 'w').write('\n'.join(lines) + '\n\n' + bt + '\n\nmaps\n' + read(f'/proc/{p}/maps'))
    log.write(f'{time.strftime("%F %T")} DUMP pid {p} episode {n}: {why} -> {path}\n')


seen, episodes = {}, 0
while True:
    now = time.time()
    alive = set()
    for p in procs():
        alive.add(p)
        try:
            state, cpu, wchan = main_state(p)
        except (OSError, ValueError):
            continue
        # Asleep anywhere but the run loop's poll (futex, sigsuspend, pipe ...): a WTF thread left in its
        # suspend handler's sigsuspend is as stuck as one waiting for a lock.
        blocked = state in 'SD' and not any(w in wchan for w in ('poll', 'select'))
        s = seen.get(p)
        if s is None or not blocked or cpu != s[1]:
            if s is not None and s[2]:
                log.write(f'{time.strftime("%F %T")} pid {p} main thread moved on after {now - s[0]:.0f} s (wchan={wchan})\n')
            seen[p] = [now, cpu, False]
            continue
        if not s[2] and now - s[0] >= BLOCK_S:
            s[2] = True
            episodes += 1
            dump(p, episodes, f'main thread asleep {now - s[0]:.0f} s (state {state}, wchan={wchan})')
    for p in list(seen):
        if p not in alive:
            del seen[p]
    time.sleep(1)
