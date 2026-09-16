# NCBI reference: c++/src/algo/blast/core/link_hsps.c:827-863
# for (H2_index=H_index-1; H2_index>1;) { ... }
# Freeze the baseline loop and compare visits plus returned scalar/FP state.
from pathlib import Path
import textwrap, subprocess, json, hashlib, os, argparse

e = Path(__file__).resolve().parent
p=argparse.ArgumentParser()
p.add_argument('--candidate',default='I1')
p.add_argument('--label',required=True)
p.add_argument('--kinds',nargs='+',default=['native','serial','threaded'])
a=p.parse_args()
d = e / a.label
d.mkdir()
base = (e / 'work/baseline-source/LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs').read_text()
candidate_path=e/'work'/(a.candidate+'-source/LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs')
candidate_hash=hashlib.sha256(candidate_path.read_bytes()).hexdigest()
frozen=json.loads((e/(a.candidate+'-source-manifest.json')).read_text())
assert frozen['LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs']==candidate_hash
candidate = candidate_path.read_text()
start = base.index('                        // Inner loop (NCBI lines 827-861)')
end = base.index('                        }\n                    } else if is_target_hsp {', start) + len('                        }')
loop = textwrap.dedent(base[start:end])
start = candidate.index('fn scan_large_gap_predecessors(')
end = candidate.index('\n}\n', start) + 3
helper = candidate[start:end]
signature = helper[:helper.index('{') + 1] + '\n'
reference = signature + textwrap.indent(loop, '    ') + '\n    (h_sum, h_num, h_xsum, h_link)\n}\n'
driver_source = 'extracted helper'
if 'macro_rules! visit_large_gap_predecessor {' in candidate:
    # NCBI link_hsps.c:827-863: H2_helper=&lh_helper[H2_index]; H2_index--;
    # Validate the actual indexed caller, not only the unit-test helper driver.
    def indexed_driver(text, after):
        start = text.index('let mut j_lh_idx = h_lh_idx - 1;', after)
        start = text.rfind('\n', 0, start) + 1
        cursor = text.index('{', text.index('while j_lh_idx > 1', start)) + 1
        depth = 1
        while depth:
            depth += (text[cursor] == '{') - (text[cursor] == '}')
            cursor += 1
        return start, cursor, textwrap.dedent(text[start:cursor])
    marker = candidate.index('// Keep the indexed scan in the caller')
    _, _, caller_driver = indexed_driver(candidate, marker)
    hs, he, test_driver = indexed_driver(helper, 0)
    assert caller_driver == test_driver
    helper = helper[:hs] + textwrap.indent(caller_driver, '        ') + helper[he:]
    ms = candidate.index('macro_rules! visit_large_gap_predecessor {')
    me = candidate.index('\n}\n', ms) + 3
    body_start = helper.index('{') + 1
    helper = helper[:body_start] + '\n' + textwrap.indent(candidate[ms:me], '    ') + helper[body_start:]
    driver_source = 'actual indexed caller for Native/serial; actual helper for threaded'

def instrument(source, name):
    source = source.replace('fn scan_large_gap_predecessors(', 'fn ' + name + '(\n    visits: &mut Vec<(usize, i32, i16, u64, usize)>,', 1)
    if 'let sum = $helper.sum[1];' in source:
        needle = 'let sum = $helper.sum[1];'
        assert source.count(needle) == 1
        return source.replace(needle, 'visits.push(($current_idx, $h_sum, $h_num, $h_xsum.to_bits(), $h_link));\n'+needle)
    needle = 'let sum = helper.sum[1];'
    assert source.count(needle) == 1
    return source.replace(needle, 'visits.push((current_idx, h_sum, h_num, h_xsum.to_bits(), h_link));\n'+needle)

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
    // NCBI prelim_stage.cpp:173-188: Run(); Join(&result).
    // Keep the real threaded WASI ABI linked for this diagnostic command.
    #[cfg(all(losat_wasi_threads,feature="wasm-threads"))]
    std::thread::spawn(|| ()).join().unwrap();
    let args:Vec<String>=std::env::args().collect();
    let scan=if args[1]=="baseline" {baseline_scan} else {candidate_scan};
    let trace=args[2]=="trace";
    let mut count=0;
    // NCBI link_hsps.c:827: H2_index=H_index-1; H2_index>1.
    // An empty scan must not read or slice any helper or payload array.
    for h_idx in [1usize, 2] {
        let mut visits=Vec::new();
        let result=scan(&mut visits,h_idx,&[],&[],&[],10,10,trace,5,-3,-0.0,999);
        println!("EMPTY {} visits={:?} result={:?}",h_idx,visits,(result.0,result.1,result.2.to_bits(),result.3));
        assert!(visits.is_empty());
        assert_eq!((result.0,result.1,result.2.to_bits(),result.3),(5,-3,(-0.0f64).to_bits(),999));
    }
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
manifest = {'source_sha256': hashlib.sha256(source.encode()).hexdigest(), 'candidate_linking_sha256':candidate_hash, 'candidate_driver_source':driver_source, 'inlining_attributes_in_harness':False, 'builds': [], 'runs': []}
def save():
    (d / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
for kind in a.kinds:
    artifact = d / ('scan-' + kind + ('.wasm' if kind != 'native' else ''))
    cmd = ['rustc', '--edition=2021', '-C', 'opt-level=3', str(d / 'harness.rs'), '-o', str(artifact)]
    if kind != 'native': cmd += ['--target', 'wasm32-wasip1' if kind=='serial' else 'wasm32-wasip1-threads']
    if kind == 'threaded':
        # build.rs emits this cfg for the real wasm32-wasip1-threads target.
        cmd += ['--cfg', 'feature="wasm-threads"', '--cfg', 'losat_wasi_threads']
        cfg=subprocess.run(['rustc','--print','cfg','--target','wasm32-wasip1-threads'],capture_output=True,text=True,check=True).stdout
        (d/'threaded-cfg.txt').write_text(cfg)
    r = subprocess.run(cmd, capture_output=True, text=True)
    manifest['builds'].append(dict(argv=cmd, returncode=r.returncode, stdout=r.stdout, stderr=r.stderr))
    save()
    assert r.returncode == 0, r.stderr
    for trace in ['off', 'trace']:
        for variant in ['baseline', 'candidate']:
            host='run_losat_wasi.js' if kind=='serial' else 'run_losat_wasi_threads.js'
            prefix = [str(artifact)] if kind == 'native' else ['node', str(e / 'work/baseline-source/LOSAT/tests' / host), str(artifact)]
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
        print(kind, trace, '800 cases + 2 empty-array boundaries, visits/state/trace exact', flush=True)
manifest['status'] = 'PASS'
assert hashlib.sha256(candidate_path.read_bytes()).hexdigest()==candidate_hash
save()
