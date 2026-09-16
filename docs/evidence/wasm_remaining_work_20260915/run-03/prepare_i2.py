# NCBI reference: c++/src/algo/blast/core/link_hsps.c:675-684,827-863,876-885
# Int4 prev=H_index-1; ... lh_helper[H_index].next_larger=prev;
# Preserve the predecessor selection body; change only the threaded-Wasm cursor.
from pathlib import Path
import shutil, textwrap, subprocess, difflib, hashlib, json, argparse
e=Path(__file__).resolve().parent
p=argparse.ArgumentParser()
p.add_argument('--version',required=True)
a=p.parse_args()
source=e/'work/I1-source'
dest=e/'work'/(a.version+'-source')
shutil.copytree(source,dest)
p=dest/'LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs'
original=p.read_bytes()
s=original.decode().replace('\r\n','\n')
begin=s.index('    // Inner loop (NCBI lines 827-861)',s.index('fn scan_large_gap_predecessors('))
end=s.index('    (h_sum, h_num, h_xsum, h_link)',begin)
loop=s[begin:end]
step=loop[loop.index('        let sum = helper.sum[1];'):loop.rfind('\n    }')]
step=textwrap.dedent(step)
assert step.count('j_lh_idx -= 1;')==1 and step.count('j_lh_idx = next_larger;')==1
step=step.replace('j_lh_idx -= 1;','$cursor = $next_cursor;').replace('j_lh_idx = next_larger;','visit_predecessor!(@jump $mode, $cursor, next_larger);')
replacement='''    // NCBI reference: c++/src/algo/blast/core/link_hsps.c:675-684,876-885
    // Int4 prev=H_index-1;
    // while((cur_sum>=prev_sum) && (prev>0)) { prev=lh_helper[prev].next_larger; }
    // lh_helper[H_index].next_larger=prev;
    // Every jump is backward. A borrowed prefix excludes sentinels 0/1;
    // its last element is the same predecessor and truncation is metadata-only.
    // One macro owns the unchanged selection and tracing statements for both
    // cursor representations. Native and plain WASI retain indexed traversal.
    macro_rules! visit_predecessor {
        (@jump indexed, $cursor:ident, $next:expr) => {
            $cursor = $next;
        };
        (@jump prefix, $cursor:ident, $next:expr) => {
            $cursor = &$cursor[..$next.saturating_sub(1)];
        };
        ($mode:ident, $source_helper:expr, $source_index:expr, $cursor:ident, $next_cursor:expr) => {{
            let current_idx = $source_index;
            let helper = $source_helper;
'''+textwrap.indent(step,'            ')+'''
        }};
    }

    #[cfg(all(losat_wasi_threads, feature = "wasm-threads"))]
    {
        // NCBI link_hsps.c:827: for(H2_index=H_index-1; H2_index>1;)
        // Preserve the no-visit case before borrowing any helper slice.
        if h_lh_idx - 1 > 1 {
            let mut prefix = &pool_lh_helpers[2..h_lh_idx];
            while let Some((helper, tail)) = prefix.split_last() {
                let current_idx = prefix.len() + 1;
                visit_predecessor!(prefix, helper, current_idx, prefix, tail);
            }
        }
    }
    #[cfg(not(all(losat_wasi_threads, feature = "wasm-threads")))]
    {
        let mut j_lh_idx = h_lh_idx - 1;
        while j_lh_idx > 1 {
            visit_predecessor!(indexed, &pool_lh_helpers[j_lh_idx], j_lh_idx, j_lh_idx, j_lh_idx - 1);
        }
    }
'''
s=s[:begin]+replacement+s[end:]
p.write_text(s)
subprocess.run(['rustfmt',str(p)],check=True)
formatted=p.read_text()
aa=original.decode().splitlines(True);zz=formatted.splitlines(True);output=[]
for tag,i,j,k,l in difflib.SequenceMatcher(None,[x.rstrip('\r\n') for x in aa],[x.rstrip('\r\n') for x in zz],autojunk=False).get_opcodes():
    output.extend(aa[i:j] if tag=='equal' else zz[k:l] if tag in ['replace','insert'] else [])
p.write_bytes(''.join(output).encode())
(e/(a.version+'-vs-I1.diff')).write_text(''.join(difflib.unified_diff(original.decode().replace('\r\n','\n').splitlines(True),formatted.splitlines(True),fromfile='I1/linking.rs',tofile=a.version+'/linking.rs')))
(e/(a.version+'-source-manifest.json')).write_text(json.dumps({str(p.relative_to(dest)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.rglob('*') if p.is_file()},indent=2)+'\n')
print(a.version,'source prepared',flush=True)
