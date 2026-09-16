# NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:173-188
# (*thread)->Run(); (*thread)->Join(&result);
# Build locally in Rust; preserve source bindings/logs/artifacts before testing.
from pathlib import Path
import argparse,json,hashlib,subprocess,time,shutil
p=argparse.ArgumentParser();p.add_argument('version');p.add_argument('--kinds',nargs='+',default=['native','threaded']);a=p.parse_args()
b=Path(__file__).resolve().parent;s=b/'work'/f'{a.version}-source/LOSAT';out=b/'work'/a.version;out.mkdir(exist_ok=True)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bindings():return {str(p.relative_to(s)):digest(p) for p in [*sorted((s/'src').rglob('*.rs')),s/'Cargo.toml',s/'Cargo.lock',s/'build.rs',s/'.cargo/config.toml']}
records=json.loads((out/'build.json').read_text()) if (out/'build.json').exists() else []
for kind in a.kinds:
 target={'native':None,'serial':'wasm32-wasip1','threaded':'wasm32-wasip1-threads','reactor':'wasm32-wasip1-threads','serial-reactor':'wasm32-wasip1'}[kind];td=Path('/tmp/losat-run03-build-cache')/a.version/kind
 argv=['cargo','build','--offline','--locked','--release','--jobs','8','--target-dir',str(td),'--lib'] if kind.endswith('reactor') else ['cargo','build','--offline','--locked','--release','--jobs','8','--target-dir',str(td),'--bin','LOSAT']
 if target:argv+=['--target',target]+(['--no-default-features'] if kind.startswith('serial') else ['--features','wasm-threads'])
 row=dict(kind=kind,argv=argv,cwd=str(s),source_sha256=bindings(),status='BUILDING',started_realtime=time.time());records.append(row);(out/'build.json').write_text(json.dumps(records,indent=2)+'\n')
 print(a.version,kind,'BUILD',flush=True)
 with (out/f'{kind}-build.log').open('w') as f:r=subprocess.run(argv,cwd=s,stdout=f,stderr=subprocess.STDOUT)
 row.update(status=r.returncode,ended_realtime=time.time());assert row['source_sha256']==bindings()
 if r.returncode:(out/'build.json').write_text(json.dumps(records,indent=2)+'\n');raise SystemExit(r.returncode)
 raw=td/target/'release/LOSAT.wasm' if target else td/'release/LOSAT';folder=out/'artifacts';folder.mkdir(exist_ok=True);name='losat-threaded-reactor.wasm' if kind=='reactor' else 'losat-serial-reactor.wasm' if kind=='serial-reactor' else f'losat-{kind}-command'+('.wasm' if target else '')
 dest=folder/name;shutil.copy2(raw,dest);assert digest(raw)==digest(dest);row.update(path=str(dest),raw_build_path=str(raw),sha256=digest(dest),bytes=dest.stat().st_size)
 (out/'build.json').write_text(json.dumps(records,indent=2)+'\n');print(kind,'PERSISTED',row['sha256'],flush=True)
