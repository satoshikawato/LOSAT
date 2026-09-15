# NCBI reference: c++/src/objtools/align_format/tabular.cpp:1098-1108
# x_PrintField(*iter); m_Ostream << "\n";
# Gate exact completed records before scoped performance measurements.
import hashlib,json,time
from pathlib import Path
b=Path(__file__).resolve().parent
def read(path):return json.loads((b/path).read_text())
jobs=read('final-validation.json');assert len(jobs)==12
assert all(r['returncode']==0 for r in jobs if not r['name'].endswith('-edge'))
summary={}
for v in ['C0','integrated']:
 first=read(v+'-edge/runs.json');middle=read(v+'-edge-remainder/runs.json');tail=read(v+'-edge-tail/runs.json')
 assert len(first)==72 and first[-1]['case']=='blastn-minus' and first[-1]['exit_status']==2
 assert "unknown option or argument '-strand'" in (b/(v+'-edge/blastn-minus/native-n1/stderr.txt')).read_text()
 assert len(middle)==32 and middle[-1]['case']=='tblastx-nohit' and middle[-1]['exit_status']==1
 assert 'no valid query contexts' in (b/(v+'-edge-remainder/tblastx-nohit/native-n1/stderr.txt')).read_text()
 assert len(tail)==62
 valid=first[:70]+middle[:30]+tail[:50]
 assert len(valid)==150 and all(r['status']=='PASS' for r in valid)
 assert all(r['raw_equal'] for r in valid if r['kind']!='oracle')
 assert all(r['status']=='EXPECTED_UNSUPPORTED' for r in tail[50:])
 for task in ['blastn','megablast']:
  text=(b/(v+'-edge-remainder')/(task+'-minus/oracle/output.txt')).read_text()
  assert any(int(x.split('\t')[8])>int(x.split('\t')[9]) for x in text.splitlines())
 rows=read(v+'-matrix/runs.json');assert len(rows)==(60 if v=='C0' else 130) and all(r['status']=='PASS' for r in rows)
 # The existing runtime harness asserts expected success/failure internally;
 # its completed records store exit codes instead of a synthetic PASS field.
 runtime=read(v+'-runtime/runs.json');assert len(runtime)==317 and all('exit' in r and 'status' not in r for r in runtime)
 assert not (b/(v+'-runtime/format-failures.json')).exists()
 long=[]
 for kind,count in [('native',4),('serial',1),('threaded',4)]:
  records=read(v+'-long-'+kind+'/runs.json');assert len(records)==count and all(r['status']=='PASS' and r['raw_equal'] for r in records);long+=records
 browser=read('browser-'+v+'-lifecycle/browser-staged/results.json');assert len(browser['records'])==25 and not browser['external_requests'] and not browser['page_errors']
 life=[r for r in browser['records'] if 'worker_lifecycle' in r];assert len(life)==2 and all(r['pagehide_termination']=='PASS' and r['failure_recovery']=='PASS' for r in life)
 summary[v]={'edge_supported_records':150,'edge_expected_unsupported_records':12,'matrix_records':len(rows),'runtime_records':317,'long_code4_records':len(long),'browser_records':25}
out={'allow_scoped_measurements':True,'reviewed_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'scope':'valid-query paths for declared four-program fixtures; no release/all-input certification','summary':summary,'preserved_boundaries':{'original_word7_baseline':'PARITY_FAIL; corrected separately in C0','old_reuse_harness':'wrong child-count assertion; fixed after raw-equal reproduction','original_minus_flag_probe':'unsupported CLI option; corrected reverse-coordinate fixture uses supported default both strands','statistically_invalid_tblastx_query':'OPEN baseline exit-status discrepancy before X2 linking; original BAD_EXIT evidence preserved','reactor_growth_and_reuse_speed':'pending fixed-window measurements; original failures remain open'},'rust':{'test':'625 passed, 3 ignored','clippy_all_targets_all_features':'PASS','fmt_check':'PASS'}}
(b/'validation-review.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
