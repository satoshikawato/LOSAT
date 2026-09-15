# NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:173-180
# (*thread)->Run(); (*thread)->Join(&result);
# Build isolated Rust artifacts; NCBI is only the separate output oracle.
import argparse,hashlib,json,shutil,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('version');p.add_argument('--kinds',nargs='+',default=['native','serial','threaded']);a=p.parse_args()
base=Path(__file__).resolve().parent; source=base/(a.version+'-source')/'LOSAT'; out=base/a.version; out.mkdir(exist_ok=True)
records=json.loads((out/'build.json').read_text()) if (out/'build.json').exists() else []
for kind in a.kinds:
    target={'native':None,'serial':'wasm32-wasip1','threaded':'wasm32-wasip1-threads','reactor':'wasm32-wasip1-threads'}[kind]
    folder={'native':'native-command','serial':'serial-command','threaded':'threaded-command','reactor':'threaded-reactor'}[kind]
    td=out/folder
    argv=['cargo','build','--offline','--locked','--release','--jobs','8','--target-dir',str(td)]
    argv+=['--lib'] if kind=='reactor' else ['--bin','LOSAT']
    if target: argv+=['--target',target,'--no-default-features'] if kind=='serial' else ['--target',target,'--features','wasm-threads']
    started=time.time(); print(a.version,kind,'build',flush=True)
    with (out/(kind+'-build.log')).open('w') as log: result=subprocess.run(argv,cwd=source,stdout=log,stderr=subprocess.STDOUT)
    rec={'argv':argv,'cwd':str(source),'status':result.returncode,'seconds':time.time()-started,'source_sha256':{str(p.relative_to(source)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [*sorted((source/'src').rglob('*.rs')),source/'Cargo.toml',source/'Cargo.lock',source/'build.rs',source/'.cargo/config.toml']}};records.append(rec)
    if result.returncode: (out/'build.json').write_text(json.dumps(records,indent=2));raise SystemExit(result.returncode)
    built=td/target/'release/LOSAT.wasm' if target else td/'release/LOSAT'
    rec.update(path=str(built),sha256=hashlib.sha256(built.read_bytes()).hexdigest(),bytes=built.stat().st_size)
    art=out/'artifacts';art.mkdir(exist_ok=True)
    if target: shutil.copyfile(built,art/('losat-threaded-reactor.wasm' if kind=='reactor' else f'losat-{kind}-command.wasm'))
    (out/'build.json').write_text(json.dumps(records,indent=2)+'\n'); print('DONE',kind,rec['seconds'],flush=True)
