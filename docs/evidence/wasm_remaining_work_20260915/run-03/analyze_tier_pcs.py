# NCBI reference: c++/src/algo/blast/core/link_hsps.c:827-863
# for (H2_index=H_index-1; H2_index>1;) { ... }
# Match observed instruction PCs to printed Wasm code ranges, independently
# checking V8 code-creation records. This is profiling evidence, not timing.
from pathlib import Path
import re, json, collections, hashlib
e = Path(__file__).resolve().parent
report = {}
for version in ['baseline', 'I1']:
    transcript = e / 'work' / (version + '-profile') / 'MjeNMV.MelaMJNV.tlosatx/threaded-n8/stdout.txt'
    text = transcript.read_text()
    ranges = []
    for block in text.split('--- WebAssembly code ---')[1:]:
        name = re.search(r'^name: (.+)$', block, re.M).group(1)
        tier = re.search(r'^compiler: (.+)$', block, re.M).group(1)
        size = int(re.search(r'^Instructions \(size = (\d+)\)', block, re.M).group(1))
        start = int(re.search(r'^(0x[0-9a-f]+)\s+0\s', block, re.M).group(1), 16)
        ranges.append(dict(name=name, tier=tier, start=start, end=start+size, instruction_bytes=size))
    isolates = []
    for logfile in sorted((e / 'work' / (version + '-profile-logs')).glob('*.log')):
        counts = collections.Counter()
        registrations = set()
        ticks = 0
        main_script = False
        worker_entry = False
        first_last = {}
        for line in logfile.open():
            if line.startswith('code-creation,'):
                parts = line.rstrip('\n').split(',')
                if len(parts) < 7: continue
                name = parts[6]
                main_script |= '/run_losat_wasi_threads.js:' in name
                worker_entry |= '/wasi_thread_host.js:244:' in name
                for index, item in enumerate(ranges):
                    if parts[4] == hex(item['start']) and name == item['name']:
                        registrations.add(index)
            elif line.startswith('tick,'):
                parts = line.split(',')
                pc = int(parts[1], 16)
                timestamp = int(parts[2])
                ticks += 1
                for index, item in enumerate(ranges):
                    if item['start'] <= pc < item['end']:
                        assert index in registrations, (logfile, pc, 'unregistered code range')
                        counts[item['tier']] += 1
                        bounds = first_last.setdefault(item['tier'], [timestamp, timestamp])
                        bounds[1] = timestamp
        isolates.append(dict(log=logfile.name, sha256=hashlib.sha256(logfile.read_bytes()).hexdigest(), main_script_compiled=main_script, worker_entry_compiled=worker_entry, all_pc_samples=ticks, target_pc_samples=dict(counts), first_last_target_sample_us=first_last, registered_target_ranges=sorted(registrations)))
    total = collections.Counter()
    for item in isolates: total.update(item['target_pc_samples'])
    report[version] = dict(ranges=ranges, isolates=isolates, target_pc_samples=dict(total), method='top instruction PC only; code-registration required before a matching sample; sampled profiles do not count all calls or unprofiled execution')
    print(version, dict(total), 'isolates', len(isolates), flush=True)
(e / 'tier-pc-evidence.json').write_text(json.dumps(report, indent=2) + '\n')
