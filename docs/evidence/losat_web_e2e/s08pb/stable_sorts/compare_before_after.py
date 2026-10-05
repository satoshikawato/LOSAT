#!/usr/bin/env python3
"""Compare BLASTP outputs of the merge binary (before 2.3) and the new binary (after 2.3)."""
import itertools, subprocess, sys, hashlib, concurrent.futures as cf
from pathlib import Path
P = Path.home()/'.cache/losat-web-gui-target'
OLD = P/'s08pb-native/release/LOSAT'; NEW = P/'s08pb-native2/release/LOSAT'
F = P/'s08pb/src-new/LOSAT/tests/fasta/outfmt0'
pairs = [('e2e_protein_query.faa','e2e_protein_subject.faa'), ('e2e_many_query.faa','e2e_many_subject.faa'),
         ('e2e_many_subject.faa','e2e_many_subject.faa'), ('e2e_many_query.faa','e2e_many_query.faa'),
         ('e2e_split_query.faa','e2e_many_subject.faa'), ('blastp_seg_query.faa','blastp_seg_subject.faa'),
         ('e2e_gen_q_lc.faa','e2e_many_subject.faa'), ('e2e_frame_tie_query.faa','e2e_many_subject.faa'),
         ('e2e_protein_query.faa','e2e_many_subject.faa'), ('e2e_titles_subject.faa','e2e_many_subject.faa')]
opts = [[], ['-task','blastp-fast'], ['-evalue','1000'], ['-evalue','1e5'], ['-task','blastp-fast','-evalue','1000'],
        ['-comp_based_stats','0'], ['-comp_based_stats','1'], ['-comp_based_stats','3'], ['-seg','yes'],
        ['-max_target_seqs','1'], ['-word_size','2','-evalue','1000'], ['-threshold','5','-evalue','1000'],
        ['-matrix','BLOSUM45','-evalue','1000'], ['-matrix','PAM30','-gapopen','9','-gapextend','1'],
        ['-comp_based_stats','0','-evalue','1e5'], ['-max_hsps','1'], ['-window_size','0','-evalue','100']]
fmts = [['-outfmt','6'], ['-outfmt','0']]
def run(b, q, s, o, f):
    r = subprocess.run([str(b),'blastp','-query',str(F/q),'-subject',str(F/s),*o,*f], capture_output=True, timeout=1800)
    return r.returncode, hashlib.sha256(r.stdout).hexdigest(), hashlib.sha256(r.stderr).hexdigest(), len(r.stdout)
def job(t):
    q,s,o,f = t
    a = run(OLD,q,s,o,f); b = run(NEW,q,s,o,f)
    return t, a, b
diff = same = 0
with cf.ThreadPoolExecutor(12) as ex:
    for t,a,b in ex.map(job, itertools.product(pairs, opts, fmts) and [(p[0],p[1],o,f) for p in pairs for o in opts for f in fmts]):
        if a == b: same += 1
        else: diff += 1; print('DIFF', ' '.join(t[0:2]), ' '.join(t[2]+t[3]), a, b, flush=True)
print(f'# same={same} diff={diff}')
