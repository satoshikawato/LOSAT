#!/usr/bin/env python3
"""Correct-output cold-process and explicit module/instance reuse measurements."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
from check_wasm_threading import fixtures
from wasm_performance import execute, digest, validate_thread_evidence

ROOT = Path(__file__).resolve().parents[2]
TESTS = ROOT / 'LOSAT/tests'

# NCBI reference: c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-180
# m_iTotalLength += (Int8) (*iter)->length;
# Fixture shape changes work distribution, never the search algorithm.
def make_cases(out):
    inputs = fixtures(out)
    def sequence(path): return ''.join(path.read_text().split('>',1)[1].splitlines()[1:])
    lc1 = sequence(TESTS/'fasta/LC738874.fasta'); lc2 = sequence(TESTS/'fasta/LC738875.fasta')
    (out/'lc-query.fasta').write_text('>q\n'+lc1[:900]+'\n')
    # NCBI reference: c++/include/algo/blast/core/blast_gapalign.h:54
    # #define MAX_DBSEQ_LEN 5000000
    # Cross the real subject chunk boundary, rather than only the old pool threshold.
    long_sequence = sequence(TESTS/'fasta/EDL933.fna')
    assert len(long_sequence) > 5_000_000
    (out/'long-query.fasta').write_text('>q\n'+long_sequence[:900]+'\n')
    (out/'long-single.fasta').write_text('>s0\n'+long_sequence+'\n')
    (out/'lc-two.fasta').write_text('>s0\n'+lc1+'\n>s1\n'+lc2+'\n')
    proteins = [''.join(record.splitlines()[1:]) for record in (TESTS/'fasta/PajaWSV.faa').read_text().split('>')[1:]]
    aa = max(proteins,key=len)[:1600]
    (out/'protein-query.fasta').write_text('>q\n'+aa[:300]+'\n')
    (out/'protein-biased.fasta').write_text(''.join(f'>s{i}\n{aa[:length]}\n' for i,length in enumerate([1600,40,1600,80,1600,120])))
    return [
        ('small','tblastx',inputs['nuc1'],inputs['nuc1'],[]),
        ('multi-subject','tblastx',inputs['nuc3'],inputs['nuc3'],[]),
        ('old-threshold','blastn',out/'lc-query.fasta',out/'lc-two.fasta',['-task','megablast']),
        ('long-single','blastn',out/'long-query.fasta',out/'long-single.fasta',['-task','megablast']),
        ('dense-biased','blastp',out/'protein-query.fasta',out/'protein-biased.fasta',[]),
    ]

# NCBI reference: c++/src/objtools/align_format/tabular.cpp:1098-1108
# x_PrintField(*iter); m_Ostream << "\n";
# Every included timing has exact raw output; invalid baseline results are archived
# and excluded rather than used as a speed baseline.
def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['candidate-dir','artifacts','baseline-dir','baseline-runners','oracle-dir','output-dir']:
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--node',default='node'); p.add_argument('--warmups',type=int,default=1);p.add_argument('--repeats',type=int,default=5)
    args=p.parse_args(); assert args.warmups>=1 and args.repeats>=5
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    cases=make_cases(out/'fixtures')
    env={k:v for k,v in os.environ.items() if not k.startswith(('LOSAT_','RAYON_')) and k!='BL2SEQ_LEGACY'}
    env.update(LC_ALL='C',NODE_NO_WARNINGS='1')
    prefixes={
        ('candidate','native'):[str((args.candidate_dir/'native-command/release/LOSAT').resolve())],
        ('candidate','serial'):[args.node,str(TESTS/'run_losat_wasi.js'),str((args.artifacts/'losat-serial-command.wasm').resolve())],
        ('candidate','threaded'):[args.node,str(TESTS/'run_losat_wasi_threads.js'),str((args.artifacts/'losat-threaded-command.wasm').resolve())],
        ('baseline','native'):[str((args.baseline_dir/'native-command/release/LOSAT').resolve())],
        ('baseline','serial'):[args.node,str((args.baseline_runners/'run_losat_wasi.js').resolve()),str((args.baseline_dir/'serial-command/wasm32-wasip1/release/LOSAT.wasm').resolve())],
        ('baseline','threaded'):[args.node,str((args.baseline_runners/'run_losat_wasi_threads.js').resolve()),str((args.baseline_dir/'threaded-command/wasm32-wasip1-threads/release/LOSAT.wasm').resolve())],
    }
    samples=[]; summaries=[];expected={};common_args={};valid={}
    def save():
        (out/'samples.json').write_text(json.dumps(samples,indent=2))
        (out/'summary.json').write_text(json.dumps(summaries,indent=2))
    metadata=dict(node=subprocess.check_output([args.node,'-p','JSON.stringify(process.versions)'],text=True),warmups=args.warmups,repeats=args.repeats,
        head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        production_source_hashes={str(p.relative_to(ROOT)):digest(p) for p in sorted((ROOT/'LOSAT/src').rglob('*.rs'))},
        hardware={'cpuinfo':Path('/proc/cpuinfo').read_text(),'meminfo':Path('/proc/meminfo').read_text(),'affinity':sorted(os.sched_getaffinity(0))},
        oracle_hashes={program:digest(args.oracle_dir/program) for program in ['blastn','blastp','tblastx']},
        baseline_runner_hashes={p.name:digest(p) for p in args.baseline_runners.glob('*.js')},
        environment={k:v for k,v in env.items() if k.startswith(('NODE_','LOSAT_','RAYON_'))},
        artifacts={f'{version}-{kind}':{'path':prefix[-1 if kind!='native' else 0],'sha256':digest(Path(prefix[-1 if kind!='native' else 0]))} for (version,kind),prefix in prefixes.items()},
        fixtures={p.name:digest(p) for p in (out/'fixtures').glob('*')}, runners={p.name:digest(p) for p in TESTS.glob('*.js')},
        boundary='cold process includes Node startup, artifact validation, guard, compilation, workers, search, output, and teardown')
    (out/'metadata.json').write_text(json.dumps(metadata,indent=2))
    # Untimed correctness and detailed stage/worker records precede timing.
    for name,program,q,s,extra in cases:
        common=['-query',str(q),'-subject',str(s),'-outfmt','6',*extra,'-out','{output}'];common_args[name]=common
        ref=execute([str((args.oracle_dir/program).resolve()),*common,'-num_threads','1'],ROOT,out/'oracle'/name,env,120)
        assert ref['status']=='PASS';expected[name]=ref['raw_output_sha256']
        for (version,kind),prefix in prefixes.items():
            for n in [1] if kind=='serial' else [1,2,4]:
                label=f'{name}/{version}-{kind}-n{n}'
                result=execute([*prefix,program,*common,'-num_threads',str(n)],ROOT,out/'diagnostic'/label,{**env,'LOSAT_WASI_THREADS_DEBUG':'1','LOSAT_TIMING':'1'},120)
                matches=result['status']=='PASS' and result['raw_output_sha256']==expected[name]
                valid[(name,version,kind,n)]=matches
                if version=='candidate':
                    assert matches,label
                    contract=validate_thread_evidence((out/'diagnostic'/label/'stderr.txt').read_text(),n,kind)
                    result['thread_contract']=contract
                    if name=='long-single' and n>1:
                        assert any(x['stage']=='subject_chunks' and x['parallel_selected'] for x in contract['thread_stages']), label
                result.update(case_id=name,version=version,kind=kind,threads=n,timed=False,phase='diagnostic',raw_equal=matches)
                samples.append(result);save()
    configurations=[(name,program,version,kind,n) for name,program,*_ in cases for version,kind in prefixes for n in ([1] if kind=='serial' else [1,2,4]) if valid[(name,version,kind,n)]]
    for repeat in range(args.warmups+args.repeats):
        order=configurations if repeat%2==0 else list(reversed(configurations))
        for name,program,version,kind,n in order:
            result=execute([*prefixes[(version,kind)],program,*common_args[name],'-num_threads',str(n)],ROOT,out/'cold'/f'repeat-{repeat}'/f'{name}-{version}-{kind}-n{n}',env,120)
            assert result['status']=='PASS' and result['raw_output_sha256']==expected[name]
            result.update(case_id=name,version=version,kind=kind,threads=n,timed=repeat>=args.warmups,repeat=repeat,phase='cold')
            samples.append(result);save()
        print('cold repetition',repeat,'PASS',flush=True)
    for name,program,version,kind,n in configurations:
        rows=[x for x in samples if x.get('timed') and (x['case_id'],x['version'],x['kind'],x['threads'])==(name,version,kind,n)]
        values=[x['wall_seconds'] for x in rows];assert len(values)==args.repeats
        summaries.append(dict(case_id=name,version=version,kind=kind,threads=n,boundary='cold-process',median_seconds=statistics.median(values),min_seconds=min(values),max_seconds=max(values),peak_rss_bytes=max(x['peak_rss_bytes'] for x in rows),raw_output_sha256=expected[name]))
    save()
    # Reuse measurements are distinct, never mixed with cold-process samples.
    for kind,mode in [('serial-command','compiled-module'),('threaded-command','compiled-module'),('threaded-reactor','same-instance')]:
        jobs=[];directory=out/'reuse'/f'{kind}-{mode}';directory.mkdir(parents=True,exist_ok=True)
        for repeat in range(args.warmups+args.repeats):
            for name,program,q,s,extra in (cases if repeat%2==0 else list(reversed(cases))):
                for n in ([1] if kind.startswith('serial') else ([1,2,4] if repeat%2==0 else [4,2,1])):
                    output=directory/f'{name}-repeat{repeat}-n{n}.out'
                    jobs.append(dict(case_id=name,program=program,query_file=str(q),subject_file=str(s),extra=[*extra,'-num_threads',str(n)],argv=[program,*common_args[name],'-num_threads',str(n)],output=str(output),threads=n,repeat=repeat,timed=repeat>=args.warmups,expected_sha256=expected[name]))
                    jobs[-1]['argv']=[str(output) if x=='{output}' else x for x in jobs[-1]['argv']]
        jobs_path=directory/'jobs.json';jobs_path.write_text(json.dumps(jobs,indent=2))
        command=[args.node,str(TESTS/'benchmark_wasi_reuse.js'),str((args.artifacts/f'losat-{kind}.wasm').resolve()),kind,mode,str(jobs_path),str(directory/'results.json')]
        (directory/'command.json').write_text(json.dumps(command,indent=2))
        result=subprocess.run(command,cwd=ROOT,env=env,capture_output=True,timeout=900)
        (directory/'stdout').write_bytes(result.stdout);(directory/'stderr').write_bytes(result.stderr)
        assert result.returncode==0,(kind,result.stderr.decode())
        print(kind,mode,'PASS',flush=True)
    print('completed',len(samples),'cold/diagnostic samples; reuse samples recorded separately',flush=True)

if __name__=='__main__': main()
