"""Read-only consumer validation for the staged or adopted LOSAT binaries."""
# NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:173-188
# (*thread)->Run(); (*thread)->Join(&result);
# Validate the consumer boundary and unchanged ordered output, not timing claims.
import argparse
import functools
import hashlib
import json
import subprocess
import threading
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlsplit
from playwright.sync_api import sync_playwright

HERE = Path('/tmp/losat-four-program-20260915/browser-baseline')
LOSAT = Path('/mnt/c/Users/genom/GitHub/LOSAT/LOSAT')
GBDRAW = Path('/mnt/c/Users/genom/GitHub/gbdraw')
WEB = GBDRAW / 'gbdraw/web'
CASES = [
    ('blastn', 'blastn', 'AP027152.fasta', 'AP027202.fasta', ['-task', 'blastn']),
    ('megablast', 'blastn', 'NZ_CP006932.fasta', 'NZ_CP006932.fasta', ['-task', 'megablast']),
    ('tblastx', 'tblastx', 'MelaMJNV.fasta', 'PemoMJNVA.fasta', []),
    ('tblastx-gencode4', 'tblastx', 'MelaMJNV.fasta', 'PemoMJNVA.fasta', ['-query_gencode', '4', '-db_gencode', '4']),
    ('blastp', 'blastp', 'WSSV.faa', 'PajaWSV.faa', []),
]
parser = argparse.ArgumentParser()
parser.add_argument('--adopted', action='store_true')
parser.add_argument('--version',default='baseline')
args = parser.parse_args()
import shutil
BASE = Path('/tmp/losat-four-program-20260915')
HERE = BASE / ('browser-' + args.version + '-lifecycle')
HERE.mkdir(exist_ok=True)
for source,dest in [('native-command/release/LOSAT','LOSAT'),('artifacts/losat-serial-command.wasm','losat.wasm'),('artifacts/losat-threaded-command.wasm','losat-threaded.wasm')]:
    shutil.copy2(BASE/args.version/source,HERE/dest)

phase = 'adopted' if args.adopted else 'staged'
out = HERE / ('browser-' + phase)
out.mkdir(exist_ok=True)
fixtures = []
for case_id, program, query, subject, extra in CASES:
    q, s = LOSAT / 'tests/fasta' / query, LOSAT / 'tests/fasta' / subject
    expected_file = HERE / (case_id + '.native.out')
    if not args.adopted and not expected_file.exists():
        argv = [str(HERE/'LOSAT'), program, '-query', str(q), '-subject', str(s), '-outfmt', '6', '-num_threads', '1', *extra, '-out', str(expected_file)]
        subprocess.run(argv, check=True, stdout=subprocess.DEVNULL)
        (HERE/(case_id+'.native.command.json')).write_text(json.dumps(argv,indent=2)+'\n')
    fixtures.append({'id':case_id,'program':program,'query':q.read_text(),'subject':s.read_text(),'extra':extra,'expected':expected_file.read_text()})

RUN = """async ({fixture, mode, threads, abort, count=1}) => {
  const service = await import('/js/services/losat.js');
  const statuses = [];
  const job = { program: fixture.program, outfmt: '6', extraArgs: fixture.extra,
    querySequenceKey: 'q', subjectSequenceKey: 's', pairIndex: 0 };
  const controller = new AbortController();
  const options = { executionMode: mode, threadsPerJob: threads,
    totalThreadBudget: threads*count, concurrency: count,
    sequences: {q:fixture.query,s:fixture.subject}, signal: controller.signal,
    onRuntimeStatus: status => statuses.push(status) };
  const started=performance.now();
  const pending = service.runLosatPairsParallel(Array.from({length:count},(_,pairIndex)=>({...job,pairIndex})), options);
  const promise = Promise.race([pending, new Promise((_, reject) => setTimeout(() => reject(new Error('LOSAT browser operation exceeded 60 seconds')), 60000))]);
  let timer;
  if (abort) timer=setTimeout(() => controller.abort(), 5);
  try { const results = await promise; return {text:results[0].text,allEqual:results.every(r=>r.text===results[0].text),resultCount:results.length,statuses,elapsed_ms:performance.now()-started}; }
  catch (error) { return {error:error.message,name:error.name,statuses}; }
  finally { if (timer) clearTimeout(timer); }
}"""
records = []
requests = []
external = []
errors = []
with sync_playwright() as pw:
    print('Launching Chromium', flush=True)
    browser = pw.chromium.launch(headless=True)
    print('Chromium launched', flush=True)
    try:
        for isolated in (False, True):
            class Handler(SimpleHTTPRequestHandler):
                def end_headers(self):
                    self.send_header('Cache-Control', 'no-store')
                    if isolated:
                        self.send_header('Cross-Origin-Opener-Policy','same-origin')
                        self.send_header('Cross-Origin-Embedder-Policy','require-corp')
                    super().end_headers()
                def log_message(self, *_):
                    pass
                def translate_path(self, path):
                    name = Path(urlsplit(path).path).name
                    if not args.adopted and urlsplit(path).path.startswith('/wasm/losat/') and name in ('losat.wasm','losat-threaded.wasm'):
                        return str(HERE/name)
                    return super().translate_path(path)
            server = ThreadingHTTPServer(('127.0.0.1',0), functools.partial(Handler,directory=str(WEB)))
            thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
            origin=f'http://127.0.0.1:{server.server_port}'
            context=browser.new_context(viewport={'width':1440,'height':1000} if isolated else {'width':390,'height':844},service_workers='block')
            def route_request(route):
                url=route.request.url
                requests.append(url)
                if urlsplit(url).netloc != urlsplit(origin).netloc:
                    external.append(url);route.abort()
                else: route.continue_()
            context.route('**/*',route_request)
            context.add_init_script("""(() => {
              const OriginalWorker=Worker, rows=[], events=[];
              window.__losatWorkerAudit={rows,events,peakThreads:0};
              window.Worker=class extends OriginalWorker {
                constructor(url,options){
                  super(url,options); if(!String(url).includes('losat'))return; const row={id:rows.length,url:String(url),alive:true,jobs:{}};rows.push(row);this.audit=row;
                  const post=this.postMessage.bind(this);
                  this.postMessage=(message,...transfer)=>{
                    if(message?.type==='run' && message.threadsPerJob){
                      row.jobs[message.id]=message.threadsPerJob;
                      const active=rows.reduce((n,r)=>n+Object.values(r.jobs).reduce((a,b)=>a+b,0),0);
                      window.__losatWorkerAudit.peakThreads=Math.max(window.__losatWorkerAudit.peakThreads,active);
                    }
                    return post(message,...transfer);
                  };
                  this.addEventListener('message',({data})=>{
                    if(Object.hasOwn(data||{},'spawnCount'))events.push({id:data.id,requested:row.jobs[data.id],spawnCount:data.spawnCount,ok:data.ok});
                    if(data?.type==='run')delete row.jobs[data.id];
                  });
                }
                terminate(){if(this.audit){this.audit.alive=false;this.audit.jobs={};}return super.terminate();}
              };
            })();""")

            page=context.new_page();page.set_default_timeout(180000)
            page.on('pageerror', lambda error: (errors.append(str(error)), print('PAGE ERROR:', error, flush=True)))
            page.on('console', lambda message: print('CONSOLE', message.type, message.text[:250], flush=True) if message.type in ('error','warning') else None)
            try:
                print('Opening', origin, isolated, flush=True)
                page.goto(origin+'/index.html',wait_until='domcontentloaded')
                page.wait_for_function('Boolean(window.__GBDRAW_APP__)')
                print('App mounted', isolated, flush=True)
                assert page.evaluate('crossOriginIsolated') == isolated
                selected=fixtures if not args.adopted else [fixtures[-1]]
                mode='threaded' if isolated else 'auto'
                n=2 if isolated else 1
                for fixture in selected:
                    for repeat in range(2):
                        print('Starting',fixture['id'],mode,n,repeat,flush=True)
                        result=page.evaluate(RUN,{'fixture':fixture,'mode':mode,'threads':n,'abort':False})
                        assert 'error' not in result, result
                        data=result.pop('text')
                        assert data==fixture['expected'], f"Output mismatch {fixture['id']} {mode} repeat{repeat}"
                        rec={'case':fixture['id'],'isolated':isolated,'mode':mode,'threads':n,'repeat':repeat,'rows':len(data.splitlines()),'sha256':hashlib.sha256(data.encode()).hexdigest(),'statuses':result['statuses'],'elapsed_ms':result['elapsed_ms']}
                        records.append(rec)
                        print(phase,fixture['id'],mode,n,repeat,'PASS',rec['rows'],flush=True)
                if isolated:
                    fixture=fixtures[-1]
                    result=page.evaluate(RUN,{'fixture':fixture,'mode':'threaded','threads':8,'abort':False})
                    assert result.get('text')==fixture['expected'], result.get('error','n8 output mismatch')
                    records.append({'case':fixture['id'],'mode':'threaded','threads':8,'matches_native':True})
                    print(phase,'threaded n8 PASS',flush=True)
                fixture=fixtures[-1]
                canceled=page.evaluate(RUN,{'fixture':fixture,'mode':mode,'threads':n,'abort':True})
                assert canceled.get('name')=='AbortError',canceled
                recovered=page.evaluate(RUN,{'fixture':fixture,'mode':mode,'threads':n,'abort':False})
                assert recovered.get('text')==fixture['expected'],recovered.get('error','recovery mismatch')
                records.append({'case':fixture['id'],'mode':mode,'cancellation':'PASS','recovery_matches_native':True})
                bad={**fixture,'query':'invalid FASTA'}
                failed=page.evaluate(RUN,{'fixture':bad,'mode':mode,'threads':n,'abort':False})
                assert 'error' in failed,failed
                recovered=page.evaluate(RUN,{'fixture':fixture,'mode':mode,'threads':n,'abort':False})
                assert recovered.get('text')==fixture['expected'],recovered
                if isolated:
                    parallel=page.evaluate(RUN,{'fixture':fixture,'mode':'threaded','threads':2,'abort':False,'count':2})
                    assert parallel.get('text')==fixture['expected'] and parallel['allEqual'] and parallel['resultCount']==2,parallel
                page.evaluate("dispatchEvent(new Event('pagehide'))")
                page.wait_for_function("window.__losatWorkerAudit.rows.every(r=>!r.alive)",timeout=5000)
                audit=page.evaluate('window.__losatWorkerAudit')
                assert audit['peakThreads']<=8,audit
                assert all(e['spawnCount']==e['requested']-1 for e in audit['events']),audit
                if isolated:assert any(e['requested']==8 and e['spawnCount']==7 for e in audit['events']),audit
                records.append({'worker_lifecycle':audit,'failure_recovery':'PASS','pagehide_termination':'PASS','isolated':isolated})

            finally:
                context.close();server.shutdown();server.server_close();thread.join()
        assert not external,external
        assert not errors,errors
    finally:
        browser.close()
        (out/'results.json').write_text(json.dumps({'phase':phase,'browser':browser.version,'records':records,'external_requests':external,'page_errors':errors,'requested_urls':requests},indent=2)+'\n')
print('Browser validation complete:',out/'results.json',flush=True)
