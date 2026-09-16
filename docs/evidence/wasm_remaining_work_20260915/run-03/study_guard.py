# NCBI reference: c++/src/app/blast/blastn_app.cpp:59-67
# m_StopWatch.Start(); m_StopWatch.Elapsed();
# Measurement-only process observation. Never changes a process or search.
from pathlib import Path
import argparse, hashlib, json, os, subprocess, time
import psutil

p = argparse.ArgumentParser()
p.add_argument('--label', required=True)
p.add_argument('command', nargs=argparse.REMAINDER)
a = p.parse_args()
command = a.command[1:] if a.command[:1] == ['--'] else a.command
assert command
e = Path(__file__).resolve().parent
out = e / (a.label + '-environment')
out.mkdir()
policy = dict(
    purpose='prospective normal-desktop comparison, not an exclusive-machine study',
    command=command, interval_seconds=1,
    background_mean_cpu_cores_limit=0.5,
    background_interval_cpu_cores_limit=1.0,
    excluded_workloads='external builds, tests, benchmarks, browser QA, compression',
    permitted_background='editor, Git indexing, ordinary OS services',
    decision='all timing/raw/memory guards AND environment validity must pass; no sample extension',
    previous_diagnostic='I4-serial-git-diagnostic remains ineligible and is not pooled',
    limitation='CPU bounds do not establish absence of I/O/cache/frequency interference',
    monitor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
)
policy_path = out / 'policy.json'
policy_path.write_text(json.dumps(policy, indent=2) + '\n')

def snapshot():
    result = {}
    for proc in psutil.process_iter(['pid', 'ppid', 'name', 'cpu_times', 'create_time', 'cmdline']):
        try:
            d = proc.info
            t = d['cpu_times']
            args = d['cmdline'] or []
            # Keep names/script basenames, not editor tokens or arbitrary arguments.
            scripts = [Path(x).name for x in args[1:] if x.endswith(('.py', '.js', '.cjs', '.sh'))]
            result[d['pid']] = dict(ppid=d['ppid'], name=d['name'], scripts=scripts,
                                    created=d['create_time'], cpu=t.user+t.system)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return result

def descendants(data, owner):
    own = {owner}
    while True:
        extended = own | {pid for pid, d in data.items() if d['ppid'] in own}
        if extended == own:
            return own
        own = extended

previous = snapshot()
last = time.monotonic()
env = dict(os.environ, LOSAT_BENCHMARK_ENVIRONMENT_POLICY=str(policy_path))
child = subprocess.Popen(command, env=env)
samples = []
forbidden = {'cargo', 'rustc', 'rust-lld', 'clang', 'gcc', 'cc1', 'pytest', 'wasm-opt', 'wasm-dis', 'gzip', 'xz', 'zstd', 'chrome', 'chromium', 'LOSAT'}
with (out / 'samples.jsonl').open('w') as stream:
    while True:
        time.sleep(1)
        current = snapshot()
        now = time.monotonic()
        own = descendants(current, child.pid) | descendants(previous, child.pid) | {os.getpid()}
        background, offenders = [], []
        for pid, d in current.items():
            if pid in own:
                continue
            old = previous.get(pid)
            delta = max(0, d['cpu'] - old['cpu']) if old and old['created'] == d['created'] else 0
            if delta:
                background.append(dict(pid=pid, name=d['name'], scripts=d['scripts'], cpu_seconds=delta))
            if d['name'] in forbidden or any(any(token in s.lower() for token in ['benchmark', 'pytest', 'test_', 'run_gate', 'browser_smoke']) for s in d['scripts']):
                offenders.append(dict(pid=pid, name=d['name'], scripts=d['scripts']))
        elapsed = now-last
        row = dict(monotonic_start=last, monotonic_end=now, elapsed=elapsed,
                   background_cpu_seconds=sum(d['cpu_seconds'] for d in background),
                   background=background, forbidden=offenders)
        samples.append(row)
        stream.write(json.dumps(row) + '\n')
        stream.flush()
        previous, last = current, now
        if child.poll() is not None:
            break
total_time = sum(r['elapsed'] for r in samples)
mean = sum(r['background_cpu_seconds'] for r in samples)/total_time
peak = max(r['background_cpu_seconds']/r['elapsed'] for r in samples)
valid = mean <= 0.5 and peak <= 1.0 and not any(r['forbidden'] for r in samples)
report = dict(returncode=child.returncode, status='PASS' if valid else 'INVALID_BACKGROUND',
              background_mean_cpu_cores=mean, background_peak_interval_cpu_cores=peak,
              seconds=total_time, policy_sha256=hashlib.sha256(policy_path.read_bytes()).hexdigest(),
              caveat=policy['limitation'])
(out / 'result.json').write_text(json.dumps(report, indent=2)+'\n')
print('ENVIRONMENT', report, flush=True)
raise SystemExit(child.returncode if child.returncode else 0 if valid else 2)
