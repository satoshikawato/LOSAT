# NCBI reference: c++/src/objtools/align_format/tabular.cpp:1098-1108
# x_PrintField(*iter); m_Ostream << "\n";
# Render already validated complete summaries; raw samples remain in JSON.
import argparse
import json
from pathlib import Path

OUT=Path('/mnt/c/Users/genom/GitHub/LOSAT/docs/evidence/wasm_four_program_20260915/run-01')

def read(name):
    return json.loads((OUT/name).read_text())

def paired_table(rows):
    lines=['| Fixture | Mode | n | Baseline seconds | Candidate seconds | Faster % | RSS delta MiB | Time/RSS guards |',
           '|---|---|---:|---:|---:|---:|---:|---|']
    for row in rows:
        a=row['baseline'];z=row['candidate']
        guards='/'.join('PASS' if row[key] else 'FAIL' for key in ['time_guard_pass','rss_guard_pass'])
        lines.append(f"| {row['case']} | {row['kind']} | {row['threads']} | {a['median_seconds']:.6f} | {z['median_seconds']:.6f} | {row['improvement_percent']:+.2f} | {(z['peak_rss_bytes']-a['peak_rss_bytes'])/1024**2:+.2f} | {guards} |")
    return '\n'.join(lines)

def ratios_table(rows):
    lines=['| Fixture | Version | Wasm mode | n | Native seconds | Wasm seconds | Wasm / Native |',
           '|---|---|---|---:|---:|---:|---:|']
    for row in rows:
        lines.append(f"| {row['case']} | {row['version']} | {row['kind']} | {row['threads']} | {row['native_seconds']:.6f} | {row['wasm_seconds']:.6f} | {row['wasm_over_native']:.3f} |")
    return '\n'.join(lines)

def cold_or_body(phase):
    filename={'cold':'integrated-cold-summary.json','body':'integrated-body-summary.json','turbofan':'turbofan-cold-summary.json'}[phase]
    reports=read(filename)
    if isinstance(reports,dict):reports=[reports]
    title={'cold':'Default Node production cold measurements','body':'Coarse engine-scope diagnostic measurements','turbofan':'Secondary TurboFan-only cold measurements'}[phase]
    sections=[f'# {title}',f'Full values and sample paths: [{filename}]({filename}). Positive percentages mean faster. Medians use all five measured runs per version and condition; one warmup is excluded by declaration. Guards are `max(5%, 50 ms)` for time and `max(10%, 16 MiB)` for peak RSS.']
    if phase=='body':
        sections.append('These are separate diagnostic builds. The timer includes configuration, file input, search and output, and excludes CLI parsing, host startup/compilation and teardown. RSS still covers the entire diagnostic process. These process/RSS observations do not replace production adoption measurements. The Native-equivalent body target is a Wasm/Native ratio of1.20 or less on each declared primary; absolute Native time must also be considered.')
    if phase=='turbofan':
        sections.append('Both versions use `--no-liftoff --no-wasm-tier-up`. Default Node remains the primary runtime. The first table uses original baseline versus standalone X2 (its baseline column is not C0); the integrated tables use C0. No TurboFan body/reuse/browser timing claim is made.')
    for report in reports:
        sections.extend(['## '+Path(report['directory']).name,paired_table(report['paired_conditions'])])
        if report['same_version_ratios']:
            sections.extend(['### Same-version Native/Wasm ratios',ratios_table(report['same_version_ratios'])])
    (OUT/(phase+'-results.md')).write_text('\n\n'.join(sections)+'\n')

def reuse():
    reports=read('integrated-reuse-summary.json')
    sections=['# Default Node fixed-window reuse measurements',
              'Full records: [integrated-reuse-summary.json](integrated-reuse-summary.json). Two independent process sessions run AB then BA, with one warmup and five measured invocations per condition in each. Combined medians use all ten invocations; these are not ten fresh-process pairs. Module compilation and reactor FASTA loading are outside invocation timing and retained as separate setup costs in JSON. API input/result copying and host worker completion remain inside invocation timing. Final persistent-instance close time is separately recorded.',
              'Reactor n1/n8 share one per-case instance within each session. Peak RSS covers the process lifetime; repeated values in per-thread tables are not per-thread memory measurements. Linear-memory vectors cover the entire ordered session, including warmups. No finite window establishes an unlimited-repeat bound, and incomplete prefixes do not establish a plateau. Existing memory and reuse failures retain their original criteria.']
    for report in reports:
        sections.extend(['## '+Path(report['directory']).name,'Run status: **'+report['status']['status']+'**.'])
        if report['combined_session_pairs']:
            sections.extend(['### Combined complete sessions',paired_table(report['combined_session_pairs'])])
        for session in [0,1]:
            rows=[row for row in report['session_pairs'] if row['session']==session]
            if rows:sections.extend([f'### Session {session}',paired_table(rows)])
        lines=['| Session | Version | Mode | Complete | Status | Peak RSS MiB | Max linear MiB | Within1GiB (full window) | Later growth (full window) |',
               '|---:|---|---|---|---|---:|---:|---|---|']
        for row in report['sessions']:
            memory='unknown' if row['memory_max_bytes'] is None else f"{row['memory_max_bytes']/1024**2:.2f}"
            rss='unknown' if row['peak_rss_bytes'] is None else f"{row['peak_rss_bytes']/1024**2:.2f}"
            value=lambda x:'unknown' if x is None else ('yes' if x else 'no')
            lines.append(f"| {row['session']} | {row['version']} | {row['kind']} | {value(row['complete'])} | {row['process_status']} | {rss} | {memory} | {value(row['complete_window_within_1gib'])} | {value(row['complete_window_growth_after_midpoint'])} |")
        sections.extend(['### Session memory and completion','\n'.join(lines)])
    (OUT/'reuse-results.md').write_text('\n\n'.join(sections)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--phase',choices=['cold','body','reuse','turbofan','all'],default='all');a=p.parse_args()
    for phase in ['cold','body','reuse','turbofan']:
        if a.phase in [phase,'all']:
            reuse() if phase=='reuse' else cold_or_body(phase)
