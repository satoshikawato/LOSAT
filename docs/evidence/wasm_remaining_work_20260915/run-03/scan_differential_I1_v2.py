# NCBI reference: c++/src/algo/blast/core/link_hsps.c:827-863
# for (H2_index=H_index-1; H2_index>1;) { ... }
# Freeze the baseline loop and compare visits plus returned scalar/FP state.
from pathlib import Path
import textwrap, subprocess, json, hashlib, os

e = Path(__file__).resolve().parent
d = e / 'scan-differential-v2'
d.mkdir()
base = (e / 'work/baseline-source/LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs').read_text()
candidate = (e / 'work/I1-source/LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs').read_text()
start = base.index('                        // Inner loop (NCBI lines 827-861)')
end = base.index('                        }\n                    } else if is_target_hsp {', start) + len('                        }')
loop = textwrap.dedent(base[start:end])
start = candidate.index('fn scan_large_gap_predecessors(')
end = candidate.index('\n}\n', start) + 3
helper = candidate[start:end]
signature = helper[:helper.index('    // Inner loop')]
reference = signature + textwrap.indent(loop, '    ') + '\n    (h_sum, h_num, h_xsum, h_link)\n}\n'

def instrument(source, name):
    source = source.replace('fn scan_large_gap_predecessors(', 'fn ' + name + '(\n    visits: &mut Vec<(usize, i32, i16, u64, usize)>,', 1)
    needle = '    while j_lh_idx > 1 {'
    assert source.count(needle) == 1
    return source.replace(needle, needle + '\n        visits.push((j_lh_idx, h_sum, h_num, h_xsum.to_bits(), h_link));')

source = '''// NCBI reference: c++/src/algo/blast/core/link_hsps.c:827-863
// Frozen test-only representations expose exactly the fields read by the scan.
#[derive(Clone,Default)]
struct LhHelper { sum:[i32;2], next_larger:usize, q_off_trim:i32, s_off_trim:i32, hsp_idx:usize }
#[derive(Clone,Default)]
struct HspLink { num:[i16;2], sum:[i32;2], xsum:[f64;2], score:i32 }
#[derive(Clone,Default)]
struct UngappedHit { q_aa_start:usize, q_aa_end:usize, s_aa_start:usize, s_aa_end:usize }
'''
source += instrument(reference, 'baseline_scan') + instrument(helper, 'candidate_scan')
source += '''
fn main() {
    let args:Vec<String>=std::env::args().collect();
    let scan=if args[1]=="baseline" {baseline_scan} else {candidate_scan};
    let trace=args[2]=="trace";
    let mut count=0;
    for n in [0usize,1,2,3,15,16,17,31,32,33] {
        for pattern in 0..16usize {
            for initial_sum in [-7,0,5,10,19] {
                let mut helpers=vec![LhHelper::default();n+2];
                let mut links=vec![HspLink::default();n];
                let mut hits=vec![UngappedHit::default();n];
                for k in 0..n {
                    let idx=k+2;
                    let mapped=n-1-k;
                    let sum=if pattern%4==0 {10} else {((k*7+pattern*3)%29) as i32-8};
                    let next=match pattern%4 {0=>0,1=>1,2=>idx-1,_=>idx/2};
                    helpers[idx]=LhHelper {sum:[-99,sum],next_larger:next,q_off_trim:((k+pattern)%3) as i32+9,s_off_trim:((k*2+pattern/2)%3) as i32+9,hsp_idx:mapped};
                    links[mapped]=HspLink {num:[-1,k as i16+2],sum:[-9,sum],xsum:[0.0,f64::from_bits(1.0f64.to_bits()+k as u64+pattern as u64)],score:sum};
                    hits[mapped]=UngappedHit {q_aa_start:k,q_aa_end:k+1,s_aa_start:k+2,s_aa_end:k+3};
                }
                let mut visits=Vec::new();
                let result=scan(&mut visits,n+2,&helpers,&links,&hits,10,10,trace,initial_sum,-3,-0.0,999);
                println!("CASE {} n={} pattern={} initial={} visits={:?} result={:?}",count,n,pattern,initial_sum,visits,(result.0,result.1,result.2.to_bits(),result.3));
                count+=1;
            }
        }
    }
    println!("CASES {}",count);
}
'''
(d / 'harness.rs').write_text(source)
manifest = {'source_sha256': hashlib.sha256(source.encode()).hexdigest(), 'builds': [], 'runs': []}
def save():
    (d / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
for kind in ['native', 'serial']:
    artifact = d / ('scan-' + kind + ('.wasm' if kind == 'serial' else ''))
    cmd = ['rustc', '--edition=2021', '-C', 'opt-level=3', str(d / 'harness.rs'), '-o', str(artifact)]
    if kind == 'serial': cmd += ['--target', 'wasm32-wasip1']
    r = subprocess.run(cmd, capture_output=True, text=True)
    manifest['builds'].append(dict(argv=cmd, returncode=r.returncode, stdout=r.stdout, stderr=r.stderr))
    save()
    assert r.returncode == 0, r.stderr
    for trace in ['off', 'trace']:
        for variant in ['baseline', 'candidate']:
            prefix = [str(artifact)] if kind == 'native' else ['node', str(e / 'work/baseline-source/LOSAT/tests/run_losat_wasi.js'), str(artifact)]
            cmd = [*prefix, variant, trace]
            env = {k: v for k, v in os.environ.items() if not k.startswith(('LOSAT_', 'NODE_', 'RAYON_'))}
            env['NODE_NO_WARNINGS'] = '1'
            r = subprocess.run(cmd, capture_output=True, env=env, timeout=60)
            stem = f'{kind}-{trace}-{variant}'
            (d / (stem + '.stdout')).write_bytes(r.stdout)
            (d / (stem + '.stderr')).write_bytes(r.stderr)
            manifest['runs'].append(dict(argv=cmd, returncode=r.returncode, stdout_sha256=hashlib.sha256(r.stdout).hexdigest(), stderr_sha256=hashlib.sha256(r.stderr).hexdigest()))
            save()
            assert r.returncode == 0, r.stderr[:1000]
        for stream in ['stdout', 'stderr']:
            assert (d / f'{kind}-{trace}-baseline.{stream}').read_bytes() == (d / f'{kind}-{trace}-candidate.{stream}').read_bytes(), (kind, trace, stream)
        print(kind, trace, '800 cases, visits/state/trace exact', flush=True)
manifest['status'] = 'PASS'
save()
