# NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:173-188
# (*thread)->Run(); (*thread)->Join(&result);
# c++/src/objtools/align_format/tabular.cpp:1098-1108: x_PrintField(*iter).
# Real consumer runtime integration, with local artifact URLs and no external
# network. These diagnostic executions are never performance samples.
from pathlib import Path
import argparse,csv,hashlib,json,re,threading,urllib.parse,os,time,importlib.metadata
from http.server import ThreadingHTTPServer,BaseHTTPRequestHandler
from playwright.sync_api import sync_playwright

e=Path(__file__).resolve().parent;web=Path('/mnt/c/users/genom/github/gbdraw/gbdraw/web')
p=argparse.ArgumentParser();p.add_argument('--version',default='I2b');p.add_argument('--label',required=True)
p.add_argument('--cases',nargs='+',default=['MjeNMV.MelaMJNV.tlosatx','AP027280.AP027280.tlosatx','AP027132.NZ_CP006932.losatp'])
p.add_argument('--threads',nargs='+',type=int,default=[8])
p.add_argument('--mode',choices=['threaded','serial','auto'],default='threaded')
p.add_argument('--no-isolation',action='store_true')
p.add_argument('--serial-kind',choices=['command','reactor'],default='command')
p.add_argument('--repeat',type=int,default=1)
p.add_argument('--lifecycle',action='store_true')
p.add_argument('--timed',action='store_true')
p.add_argument('--consumer-root',default=str(web))
a=p.parse_args();out=e/a.label;out.mkdir()
web=Path(a.consumer_root).resolve()
assert not (a.timed and a.lifecycle)
script_bytes=Path(__file__).read_bytes();(out/'harness.py').write_bytes(script_bytes)
tests=e/'work/baseline-source/LOSAT/tests';rows=list(csv.DictReader((tests/'comparison_cases.tsv').open(),delimiter='\t'))
sha=lambda data:hashlib.sha256(data).hexdigest()
csp=re.search(r'<meta http-equiv="Content-Security-Policy" content="(.*?)"', (web/'index.html').read_text(),re.S).group(1).replace('\n',' ')
requests=[];served={};artifacts={}
for kind in ['serial','threaded']:
    path=e/'work'/a.version/'artifacts'/('losat-'+kind+'-'+(a.serial_kind if kind=='serial' else 'command')+'.wasm');artifacts[kind]=path
# NCBI prelim_stage.cpp:173-188: Run(); Join(&result).
# Test-only Worker observation, never included in performance runs or shipped.
worker_probe='''
const qaBaseWorker=globalThis.Worker;
let qaNextWorker=0;
globalThis.Worker=class extends qaBaseWorker {
  constructor(url,options) {
    super(url,options); this.qaId=++qaNextWorker;
    self.postMessage({__qaWorker:{op:'create',id:this.qaId,url:String(url)}});
  }
  terminate() {self.postMessage({__qaWorker:{op:'terminate',id:this.qaId}});return super.terminate();}
};
'''
class Handler(BaseHTTPRequestHandler):
    def log_message(self,*args):pass
    def do_GET(self):
        path=urllib.parse.urlparse(self.path).path;requests.append(path)
        if path=='/qa.html':data=b'<!doctype html><html><head><meta charset="utf-8"><title>LOSAT runtime verification</title></head><body>Runtime verification</body></html>';mime='text/html'
        else:
            target=artifacts.get(path.removeprefix('/qa-artifacts/').removesuffix('.wasm')) if path.startswith('/qa-artifacts/') else (web/path.lstrip('/')).resolve()
            if target is None or not target.is_file() or (target not in artifacts.values() and not target.is_relative_to(web)):
                self.send_error(404);return
            data=target.read_bytes();served[str(target)]=sha(data)
            if a.lifecycle and target.name=='losat-threaded-worker.js':
                data=worker_probe.encode()+data;served[str(target)+':diagnostic-overlay']=sha(data)
            mime='application/wasm' if target.suffix=='.wasm' else 'text/javascript' if target.suffix in ['.js','.mjs'] else 'application/octet-stream'
        self.send_response(200);self.send_header('Content-Type',mime);self.send_header('Content-Length',str(len(data)))
        if not a.no_isolation:
            self.send_header('Cross-Origin-Opener-Policy','same-origin');self.send_header('Cross-Origin-Embedder-Policy','require-corp')
        self.send_header('Content-Security-Policy',csp);self.end_headers();self.wfile.write(data)
server=ThreadingHTTPServer(('127.0.0.1',0),Handler);thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
origin='http://127.0.0.1:'+str(server.server_port);records=[];blocked=[];console=[]
manifest=dict(version=a.version,harness_sha256=sha(script_bytes),role='full real-consumer call measurement with boundary clocks only' if a.timed else 'correctness-only real consumer service/worker smoke; excludes timings',options=vars(a),csp=csp,artifacts={k:dict(path=str(v),sha256=sha(v.read_bytes())) for k,v in artifacts.items()},records=records)
rss_peak=[0];monitor_stop=threading.Event()
def monitor_rss():
    # NCBI blastn_app.cpp:59-67: stopwatch outside the search. Host-only
    # resource observation: aggregate descendant RSS includes shared pages
    # per process, not unique physical memory or Wasm linear-memory size.
    import psutil
    parent=psutil.Process(os.getpid())
    while not monitor_stop.is_set():
        total=0
        for child in parent.children(recursive=True):
            try:total+=child.memory_info().rss
            except (psutil.NoSuchProcess,psutil.AccessDenied):pass
        rss_peak[0]=max(rss_peak[0],total);monitor_stop.wait(0.05)
if a.timed:
    import psutil
    monitor=threading.Thread(target=monitor_rss,daemon=True);monitor.start()
def save():
    manifest.update(requests=requests,served_sha256=served,blocked_external=blocked,console=console)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
save()
try:
    with sync_playwright() as pw:
        browser=pw.chromium.launch();manifest['browser']=browser.version
        manifest.update(playwright=importlib.metadata.version('playwright'),browser_executable=pw.chromium.executable_path,explicit_browser_launch_options={},debugger_domain_enabled_by_harness=False)
        context=browser.new_context();context.set_default_timeout(300000)
        def route(r):
            if r.request.url.startswith(origin+'/'):r.continue_()
            else:blocked.append(r.request.url);r.abort()
        context.route('**/*',route)
        page=context.new_page();page.on('console',lambda m:console.append(dict(type=m.type,text=m.text)))
        page.on('pageerror',lambda error:console.append(dict(type='pageerror',text=str(error))))
        page.goto(origin+'/qa.html');assert page.evaluate('crossOriginIsolated') is (not a.no_isolation)
        if a.lifecycle:
            page.evaluate('''() => {
              window.qaWorkers=[];window.qaEvents=[];
              const Base=Worker;
              window.Worker=class extends Base {
                constructor(url,options) {
                  super(url,options);this.qaId=qaWorkers.length;this.qaDead=false;qaWorkers.push(this);
                  qaEvents.push({op:'parent-create',id:this.qaId,url:String(url)});
                  this.addEventListener('message',e=>{if(e.data.__qaWorker)qaEvents.push({...e.data.__qaWorker,parent:this.qaId});});
                }
                postMessage(message,transfer){qaEvents.push({op:'parent-send',id:this.qaId,type:message.type});return super.postMessage(message,transfer);}
                terminate(){this.qaDead=true;qaEvents.push({op:'parent-terminate',id:this.qaId});return super.terminate();}
              };
            }''')
        page.evaluate('''() => {window.qaRun=async p => {
            const {runLosatPairsParallel}=await import('/js/services/losat.js');
            const statuses=[];
            try {
              const result=await runLosatPairsParallel([{program:p.program,querySequenceKey:'q',subjectSequenceKey:'s',pairIndex:0,outfmt:'6',extraArgs:p.extra}],
                {executionMode:p.mode,threadedWasmPath:'/qa-artifacts/threaded.wasm',wasmPath:'/qa-artifacts/serial.wasm',threadsPerJob:p.threads,totalThreadBudget:p.threads,concurrency:1,sequences:{q:p.query,s:p.subject},signal:p.cancelable?window.qaController.signal:undefined,onRuntimeStatus:s=>statuses.push(s)});
              return {ok:true,text:result[0].text,statuses};
            } catch(error) {return {ok:false,error:String(error),name:error.name,canceled:!!error.canceled,stack:error.stack,statuses};}
        };}''')
        for repeat,name in [(i,name) for i in range(a.repeat) for name in a.cases]:
            row=next(r for r in rows if r['losat_stem']==name);query=tests/'fasta'/row['query'];subject=tests/'fasta'/row['subject']
            expected=(e/'work/oracles'/name/'output.txt').read_bytes()
            extra=['-query_gencode',row['query_gencode'],'-db_gencode',row['db_gencode']] if row['task']=='tblastx' else []
            for n in a.threads:
                record=dict(case=name,repeat=repeat,threads=n,query_sha256=sha(query.read_bytes()),subject_sha256=sha(subject.read_bytes()),extra=extra,status='RUNNING');records.append(record);save()
                payload=dict(program=row['task'],query=query.read_text(),subject=subject.read_text(),extra=extra,threads=n,mode=a.mode)
                # NCBI blastn_app.cpp:59-67: m_StopWatch.Start(); Elapsed().
                # Include import/dispatch/input transfer/search/output/lifecycle
                # through the actual service promise; output checking follows.
                result=page.evaluate('''async p=>{
                    const started=performance.now();const result=await qaRun(p);
                    return {...result,wall_seconds:(performance.now()-started)/1000};
                }''',payload)
                if a.timed:record.update(timed=repeat>0,browser_tree_peak_rss_bytes=rss_peak[0])
                record.update(result);raw=record.pop('text','').encode();(out/(name+'-n'+str(n)+'-repeat'+str(repeat)+'.txt')).write_bytes(raw)
                record.update(raw_equal=raw==expected,output_sha256=sha(raw),status='PASS' if result['ok'] and raw==expected else 'FAIL');save()
                print(name,n,record['status'],record.get('error',''),flush=True)
                assert record['status']=='PASS',record
                if a.no_isolation and a.mode=='auto':assert any(s.get('mode')=='serial' for s in record['statuses'])
        if a.lifecycle:
            # NCBI blastn_app.cpp:172-176: CATCH_ALL(status); return status.
            # Explicit errors/cancellation must settle and permit a fresh job.
            invalid=page.evaluate('p=>qaRun({...p,program:"invalid"})',payload)
            assert not invalid['ok'];manifest['invalid_input']=invalid
            recovery=page.evaluate('p=>qaRun(p)',payload);assert recovery['ok'] and recovery['text'].encode()==expected
            manifest['invalid_input_recovery_sha256']=sha(recovery['text'].encode());save()
            page.evaluate('p=>{window.qaBeforeRuns=qaEvents.filter(e=>e.op==="parent-send"&&e.type==="run").length;window.qaController=new AbortController();window.qaPending=qaRun({...p,cancelable:true});}',payload)
            page.wait_for_function('qaEvents.filter(e=>e.op==="parent-send"&&e.type==="run").length>qaBeforeRuns')
            page.evaluate('qaController.abort()');canceled=page.evaluate('qaPending')
            assert not canceled['ok'] and (canceled['name']=='AbortError' or canceled['canceled']);manifest['canceled']=canceled
            recovery=page.evaluate('p=>qaRun(p)',payload);assert recovery['ok'] and recovery['text'].encode()==expected
            manifest['cancellation_recovery_sha256']=sha(recovery['text'].encode());save()
            events=page.evaluate('qaEvents');parents={};maximum=0
            for event in events:
                if event['op']=='parent-create':parents[event['id']]=set()
                elif event['op']=='parent-terminate':parents.pop(event['id'],None)
                elif event.get('parent') in parents:
                    if event['op']=='create':parents[event['parent']].add(event['id'])
                    else:parents[event['parent']].discard(event['id'])
                maximum=max(maximum,sum(1+len(children) for children in parents.values()))
            assert maximum<=max(a.threads)
            if a.mode=='threaded':assert maximum==max(a.threads)
            manifest['observed_worker_maximum']=maximum
            page.evaluate('window.dispatchEvent(new Event("pagehide"))')
            assert page.evaluate('qaWorkers.every(w=>w.qaDead)')
            manifest['pagehide_all_parent_workers_terminated']=True;manifest['worker_events']=page.evaluate('qaEvents');save()
        assert not blocked,blocked
        manifest['status']='PASS';save();browser.close()
finally:
    monitor_stop.set()
    if a.timed:monitor.join();manifest['rss_metric']='peak sampled sum of descendant RSS including Playwright driver and browser; shared pages count in each process; 50 ms samples'
    save();server.shutdown();server.server_close()
