# NCBI reference: c++/src/objtools/align_format/tabular.cpp:1098-1108
# x_PrintField(*iter); m_Ostream << "\n";
# Preserve raw equality; summarize the first three complete pairs authorized
# by the user's stop instruction, retaining the interrupted original run.
import hashlib
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

B = Path(sys.argv[1] if len(sys.argv) > 1 else '/tmp/losat-four-program-20260915')
OUT = Path(__file__).resolve().parent
rows = json.loads((B / 'integrated-cold/samples.json').read_text())
groups = defaultdict(list)
for row in rows:
    if row.get('timed') and row.get('repeat') in (1, 2, 3):
        assert row['phase'] == 'cold' and row['status'] == 'PASS'
        assert row['raw_equal'] and row['monotonic_clock_agreement']
        assert math.isfinite(row['wall_seconds']) and row['wall_seconds'] > 0
        assert hashlib.sha256(Path(row['output']).read_bytes()).hexdigest() == row['raw_output_sha256']
        groups[(row['case_id'], row['kind'], row['threads'], row['version'])].append(row)
assert len(groups) == 90
summaries = {}
for key, group in groups.items():
    group.sort(key=lambda row: row['repeat'])
    assert [row['repeat'] for row in group] == [1, 2, 3]
    assert len({row['raw_output_sha256'] for row in group}) == 1
    values = [row['wall_seconds'] for row in group]
    summaries[key] = dict(median_seconds=statistics.median(values), samples=values,
                          peak_rss_bytes=max(row['peak_rss_bytes'] for row in group),
                          raw_sha256=group[0]['raw_output_sha256'],
                          evidence=[row['output'] for row in group])
pairs = []
ratios = []
for (case, kind, n, version), candidate in summaries.items():
    if version == 'candidate':
        baseline = summaries[(case, kind, n, 'baseline')]
        assert baseline['raw_sha256'] == candidate['raw_sha256']
        delta = candidate['median_seconds'] - baseline['median_seconds']
        rss_delta = candidate['peak_rss_bytes'] - baseline['peak_rss_bytes']
        pairs.append(dict(case=case, kind=kind, threads=n, baseline=baseline, candidate=candidate,
                          improvement_percent=-100 * delta / baseline['median_seconds'],
                          five_percent_target_met=delta <= -.05 * baseline['median_seconds'],
                          time_guard_pass=delta <= max(.05, .05 * baseline['median_seconds']),
                          rss_guard_pass=rss_delta <= max(16 * 1024**2, .1 * baseline['peak_rss_bytes'])))
    if kind != 'native':
        native = summaries[(case, 'native', n, version)]
        ratios.append(dict(case=case, kind=kind, threads=n, version=version,
                           wasm_over_native=candidate['median_seconds'] / native['median_seconds']))
pairs.sort(key=lambda row: (row['case'], row['kind'], row['threads']))
excluded = [row for row in rows if row.get('timed') and row.get('repeat') not in (1, 2, 3)]
report = dict(status='STOPPED_BY_USER; three complete pairs summarized',
              selection='Exactly repeats 1, 2, 3 for every condition; repeat 0 is warmup. No speed-based selection.',
              original_process_records=len(rows), timed_records_used=270,
              completed_later_records_excluded=len(excluded),
              exclusion_reason='User ended measurement and specified three repetitions; interrupted repeat 4 is retained.',
              not_run=['final integrated controls timing', 'body timing', 'compiled-module reuse timing',
                       'same-reactor timing and memory window', 'TurboFan secondary timing', 'historical source rebuilds'],
              paired_conditions=pairs, same_version_ratios=ratios)
(OUT / 'integrated-three-pair-summary.json').write_text(json.dumps(report, indent=2) + '\n')
lines = ['# Integrated production cold results: first three pairs', '',
         'C0 versus integrated; one warmup, then exactly repeats 1–3 at every condition. '
         'Default Node 26.8.2, affinity CPUs 0–7. Positive change means faster. '
         'The user stopped measurement during repeat 4; later records remain preserved and are excluded by count.', '',
         '| Case | Runtime | n | C0 median (s) | Integrated median (s) | Faster (%) | Time guard | RSS guard |',
         '|---|---|---:|---:|---:|---:|---|---|']
for row in pairs:
    lines.append(f"| {row['case']} | {row['kind']} | {row['threads']} | {row['baseline']['median_seconds']:.6f} | {row['candidate']['median_seconds']:.6f} | {row['improvement_percent']:+.2f} | {'PASS' if row['time_guard_pass'] else 'FAIL'} | {'PASS' if row['rss_guard_pass'] else 'FAIL'} |")
lines += ['', 'Time guard: regression ≤ max(5%, 50 ms). RSS guard: increase ≤ max(10%, 16 MiB). '
          'All 270 retained timed outputs were rehashed and equal their recorded raw hashes; '
          'each condition has one shared C0/integrated output hash. All samples and same-version ratios are in the JSON.', '',
          '## Same-version Wasm / native cold ratios', '',
          '| Case | Runtime | n | C0 | Integrated |', '|---|---|---:|---:|---:|']
ratio_map = {(r['case'], r['kind'], r['threads'], r['version']): r['wasm_over_native'] for r in ratios}
for case, kind, n, version in sorted(ratio_map):
    if version == 'candidate':
        lines.append(f'| {case} | {kind} | {n} | {ratio_map[(case, kind, n, "baseline")]:.3f} | {ratio_map[(case, kind, n, version)]:.3f} |')
(OUT / 'integrated-three-pair-results.md').write_text('\n'.join(lines) + '\n')
print(json.dumps(dict(pairs=len(pairs), used=270, later_excluded=len(excluded),
                      time_failures=[{k: r[k] for k in ['case','kind','threads','improvement_percent']} for r in pairs if not r['time_guard_pass']],
                      rss_failures=[{k: r[k] for k in ['case','kind','threads']} for r in pairs if not r['rss_guard_pass']],
                      five_percent_targets=sum(r['five_percent_target_met'] for r in pairs)), indent=2))
