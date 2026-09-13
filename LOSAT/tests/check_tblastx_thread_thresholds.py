#!/usr/bin/env python3
"""Preserve the LC738874/LC738875 E-value regression on explicit artifacts."""
import argparse
import json
import os
from pathlib import Path
import subprocess
from wasm_performance import execute, digest, validate_thread_evidence

ROOT = Path(__file__).resolve().parents[2]
TESTS = ROOT / 'LOSAT/tests'

# NCBI reference: c++/src/algo/blast/core/blast_hits.c:1984-2003
# if (hsp->evalue > cutoff) { hsp_array[index] = Blast_HSPFree(hsp_array[index]); }
# A threshold sweep compares full raw rows, never only the number of HSPs.
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['native', 'serial', 'threaded', 'oracle', 'output-dir']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--node', default='node')
    args = parser.parse_args(); out = args.output_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    fixtures = out/'fixtures'; fixtures.mkdir(exist_ok=True)
    for name in ['LC738874.fasta', 'LC738875.fasta']:
        (fixtures/name).write_bytes((TESTS/'fasta'/name).read_bytes())
    metadata = dict(head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        artifacts={name:{'path':str(getattr(args,name).resolve()),'sha256':digest(getattr(args,name))} for name in ['native','serial','threaded','oracle']},
        fixtures={p.name:digest(p) for p in fixtures.iterdir()},
        node=subprocess.check_output([args.node,'-p','JSON.stringify(process.versions)'],text=True),
        runners={p.name:digest(p) for p in TESTS.glob('*.js')})
    (out/'metadata.json').write_text(json.dumps(metadata,indent=2))
    env = {k:v for k,v in os.environ.items() if not k.startswith(('LOSAT_','RAYON_')) and k!='BL2SEQ_LEGACY'}
    env.update(LC_ALL='C',NODE_NO_WARNINGS='1')
    records=[]
    for threshold in [10,100,10000]:
        common=['-query',str(fixtures/'LC738874.fasta'),'-subject',str(fixtures/'LC738875.fasta'),'-outfmt','6','-evalue',str(threshold),'-out','{output}']
        oracle=execute([str(args.oracle.resolve()),*common,'-num_threads','1'],ROOT,out/f'e{threshold}'/'oracle',env,900)
        assert oracle['status']=='PASS'
        records.append({**oracle,'evalue':threshold,'kind':'oracle'})
        for kind,n,prefix in [
            ('native',1,[str(args.native.resolve())]),
            ('serial',1,[args.node,str(TESTS/'run_losat_wasi.js'),str(args.serial.resolve())]),
            ('threaded',4,[args.node,str(TESTS/'run_losat_wasi_threads.js'),str(args.threaded.resolve())]),
        ]:
            directory=out/f'e{threshold}'/kind
            result=execute([*prefix,'tblastx',*common,'-num_threads',str(n)],ROOT,directory,{**env,'LOSAT_WASI_THREADS_DEBUG':'1'},900)
            result.update(evalue=threshold,kind=kind,threads=n,expected_sha256=oracle['raw_output_sha256'])
            records.append(result);(out/'runs.json').write_text(json.dumps(records,indent=2))
            assert result['status']=='PASS' and result['raw_output_sha256']==oracle['raw_output_sha256'],(threshold,kind)
            validate_thread_evidence((directory/'stderr.txt').read_text(),n,kind)
            print(threshold,kind,'raw parity PASS',flush=True)
    print(len(records),'threshold records PASS',flush=True)

if __name__=='__main__': main()
