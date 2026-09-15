# NCBI reference: c++/src/objtools/align_format/tabular.cpp:1098-1108
# x_PrintField(*iter); m_Ostream << "\n";
# Admit only complete, raw-equal finite sample series; retain all other outcomes.
import json
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path

B = Path(__file__).resolve().parent
OUT = Path('/mnt/c/Users/genom/GitHub/LOSAT/docs/evidence/wasm_four_program_20260915/run-01')

def read(path):
    return json.loads(path.read_text())

def write(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2) + '\n')

def comparison(a, z):
    delta = z['median_seconds'] - a['median_seconds']
    rss_delta = z['peak_rss_bytes'] - a['peak_rss_bytes']
    return dict(baseline=a, candidate=z,
                improvement_percent=-100 * delta / a['median_seconds'],
                five_percent_target_met=delta <= -.05 * a['median_seconds'],
                time_guard_pass=delta <= max(.05, .05 * a['median_seconds']),
                rss_guard_pass=rss_delta <= max(16 * 1024**2, .1 * a['peak_rss_bytes']))

def cold(directory, body=False):
    status = read(directory / 'run-status.json')
    assert status['status'] == 'COMPLETE', (directory, status)
    assert not read(directory / 'excluded.json')
    data = read(directory / 'samples.json')
    groups = defaultdict(list)
    for row in data:
        assert row['status'] == 'PASS', (directory, row)
        if row.get('timed'):
            assert row['phase'] == 'cold' and row['raw_equal']
            assert row['monotonic_clock_agreement']
            groups[(row['case_id'], row['kind'], row['threads'], row['version'])].append(row)
    summaries = {}
    for key, rows in groups.items():
        assert sorted(row['repeat'] for row in rows) == [1, 2, 3, 4, 5], key
        if body:
            values = []
            for row in rows:
                stderr = Path(row['output']).parent / 'stderr.txt'
                matches = re.findall(r'^\[BODY_SCOPE_SECONDS\] ([0-9.]+)$', stderr.read_text(), re.M)
                assert len(matches) == 1, stderr
                value = float(matches[0])
                assert math.isfinite(value) and value > 0, (stderr, value)
                values.append(value)
        else:
            values = [row['wall_seconds'] for row in rows]
        summaries[key] = dict(median_seconds=statistics.median(values), samples=values,
                              peak_rss_bytes=max(row['peak_rss_bytes'] for row in rows),
                              raw_sha256=rows[0]['raw_output_sha256'],
                              evidence=[row['output'] for row in rows])
    pairs = []
    ratios = []
    for (case, kind, n, version), row in summaries.items():
        if version == 'candidate':
            old = summaries[(case, kind, n, 'baseline')]
            assert old['raw_sha256'] == row['raw_sha256']
            pairs.append(dict(case=case, kind=kind, threads=n, **comparison(old, row)))
        if kind != 'native':
            native = summaries.get((case, 'native', n, version))
            if native is None:
                continue
            ratios.append(dict(case=case, kind=kind, threads=n, version=version,
                               native_seconds=native['median_seconds'], wasm_seconds=row['median_seconds'],
                               wasm_over_native=row['median_seconds'] / native['median_seconds']))
    return dict(scope='coarse engine diagnostic' if body else 'production cold process',
                directory=str(directory), process_records=len(data), paired_conditions=pairs,
                same_version_ratios=ratios)

def reuse(directory):
    status = read(directory / 'run-status.json')
    assert status['status'] != 'RUNNING', directory
    groups = {}
    sessions = []
    for session_dir in sorted((directory / 'reuse').iterdir()):
        session_text, version, kind, subtype, mode_a, mode_b = session_dir.name.split('-')
        session = int(session_text.removeprefix('session'))
        kind = kind + '-' + subtype
        mode = mode_a + '-' + mode_b
        process = read(session_dir / 'process/result.json')
        jobs = read(session_dir / 'jobs.json')
        output = Path(process['output'])
        schema_error = None
        try:
            data = read(output)
            rows = data['samples']
            assert isinstance(rows, list) and all(isinstance(row, dict) and 'status' in row for row in rows)
        except (OSError, ValueError, KeyError, TypeError, AssertionError) as error:
            schema_error = str(error)
            data = dict(samples=[])
            rows = []
        complete = (process['status'] == 'PASS' and process['monotonic_clock_agreement']
                    and len(rows) == len(jobs) and all(row['status'] == 'PASS' for row in rows))
        if complete:
            assert all(row['raw_equal'] and row['thread_contract'] == 'PASS' for row in rows)
            assert data['node_argv'] == ['/home/kawato/.local/lib/nodejs/node-v26.8.2-linux-x64/bin/node']
            for n in sorted({row['threads'] for row in rows}):
                timed = [row for row in rows if row['threads'] == n and row['timed']]
                assert sorted(row['repeat'] for row in timed) == [1, 2, 3, 4, 5]
                values = [row['wall_seconds'] for row in timed]
                groups[(jobs[0]['case_id'], kind, n, version, session)] = dict(
                    median_seconds=statistics.median(values), samples=values,
                    peak_rss_bytes=process['peak_rss_bytes'],
                    session_memory_max_bytes=max(row['memory_bytes'] for row in rows),
                    process_evidence=str(output))
        memory = [row['memory_bytes'] for row in rows if 'memory_bytes' in row]
        sessions.append(dict(case=jobs[0]['case_id'], session=session, version=version, kind=kind,
                             mode=mode, complete=complete, process_status=process['status'],
                             schema_error=schema_error,
                             planned_jobs=len(jobs), observed_jobs=len(rows),
                             statuses=[row['status'] for row in rows],
                             peak_rss_bytes=process['peak_rss_bytes'], memory_bytes=memory,
                             memory_max_bytes=max(memory) if memory else None,
                             observed_within_1gib=bool(memory) and max(memory) <= 1024**3,
                             complete_window_within_1gib=(max(memory) <= 1024**3) if complete else None,
                             complete_window_growth_after_midpoint=(max(memory) > max(memory[:max(1, len(memory)//2)])) if complete else None,
                             compile_and_input_timings=data.get('timings'),
                             process_evidence=str(output)))
    pairs = []
    combined = []
    for (case, kind, n, version, session), row in groups.items():
        if version != 'candidate':
            continue
        old = groups.get((case, kind, n, 'baseline', session))
        if old:
            pairs.append(dict(case=case, kind=kind, threads=n, session=session, **comparison(old, row)))
        if session == 0 and all((case, kind, n, v, s) in groups for v in ['baseline', 'candidate'] for s in [0, 1]):
            pooled = {}
            for v in ['baseline', 'candidate']:
                samples = [t for s in [0, 1] for t in groups[(case, kind, n, v, s)]['samples']]
                pooled[v] = dict(samples=samples, median_seconds=statistics.median(samples),
                                 peak_rss_bytes=max(groups[(case, kind, n, v, s)]['peak_rss_bytes'] for s in [0, 1]))
            combined.append(dict(case=case, kind=kind, threads=n, **comparison(pooled['baseline'], pooled['candidate'])))
    return dict(directory=str(directory), status=status, sessions=sessions,
                session_pairs=pairs, combined_session_pairs=combined,
                note='Process-lifetime RSS applies to the whole session; it is not attributable to a single invocation or thread count. The finite mixed-n1/n8 window does not establish a bound for unlimited repeats.')

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--phase', choices=['cold', 'body', 'reuse', 'turbofan', 'all'], default='all')
    args = p.parse_args()
    if args.phase in ['cold', 'all']:
        reports=[cold(B / name) for name in ['integrated-cold', 'integrated-cold-controls']]
        assert [len(report['paired_conditions']) for report in reports] == [45,12]
        assert [report['process_records'] for report in reports] == [639,174]
        write('integrated-cold-summary.json', reports)
    if args.phase in ['body', 'all']:
        report=cold(B / 'integrated-body', body=True)
        assert len(report['paired_conditions']) == 36 and report['process_records'] == 513
        write('integrated-body-summary.json', report)
    if args.phase in ['reuse', 'all']:
        reports=[reuse(B / ('integrated-reuse-' + str(i))) for i in range(9)]
        assert all(len(report['sessions']) == 12 for report in reports)
        write('integrated-reuse-summary.json', reports)
    if args.phase in ['turbofan', 'all']:
        reports=[cold(B / name) for name in ['X2-turbofan-performance', 'integrated-turbofan-cold', 'integrated-turbofan-controls']]
        assert [len(report['paired_conditions']) for report in reports] == [3,10,2]
        assert [report['process_records'] for report in reports] == [45,142,29]
        write('turbofan-cold-summary.json', reports)
