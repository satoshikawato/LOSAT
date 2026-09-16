# NCBI reference: c++/src/algo/blast/core/link_hsps.c:827-863
# for (H2_index=H_index-1; H2_index>1;) { ... }
# Preserve the original predecessor scan at an actual per-HSP function boundary.
from pathlib import Path
import tarfile,shutil,hashlib,json,subprocess,textwrap,difflib
b=Path(__file__).resolve().parent;root=b.parents[3];old=b.parent/'run-01';work=b/'work';work.mkdir(exist_ok=True)
base=work/'baseline-source';base.mkdir(exist_ok=True)
with tarfile.open(old/'baseline-inputs.tar.gz') as t:
 for m in t.getmembers():
  assert (m.name=='baseline-source' or m.name.startswith('baseline-source/')) and '..' not in Path(m.name).parts
  t.extract(m,work,filter='data')
with tarfile.open(old/'supplemental-tests.tar.gz') as t:
 for m in t.getmembers():
  if m.name.startswith('LOSAT/') and '..' not in Path(m.name).parts:t.extract(m,base,filter='data')
with tarfile.open(old/'supplemental-integration-fixtures.tar.gz') as t:
 for m in t.getmembers():
  assert (m.name=='baseline-source' or m.name.startswith('baseline-source/')) and '..' not in Path(m.name).parts
  t.extract(m,work,filter='data')
# Current Rust is frozen round02 baseline; copy the accepted H1 host owners and
# their dependencies from the actual working tree, retaining their exact bytes.
for p in (root/'LOSAT/tests').iterdir():
 if p.is_file() and p.suffix in {'.js','.py','.tsv'}:shutil.copy2(p,base/'LOSAT/tests'/p.name)
for p in [*sorted((root/'LOSAT/src').rglob('*.rs')),root/'LOSAT/Cargo.toml',root/'LOSAT/Cargo.lock',root/'LOSAT/build.rs',root/'LOSAT/.cargo/config.toml']:
 rel=p.relative_to(root);assert (base/rel).read_bytes()==p.read_bytes(),rel
manifest={str(p.relative_to(base)):hashlib.sha256(p.read_bytes()).hexdigest() for p in base.rglob('*') if p.is_file()}
(b/'baseline-source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
dest=work/'I1-source';shutil.copytree(base,dest);p=dest/'LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs';raw=p.read_bytes();s=raw.decode();clean=s.replace('\r\n','\n')
start=clean.index('                        // Inner loop (NCBI lines 827-861)')
end=clean.index('                        }\n                    } else if is_target_hsp {',start)+len('                        }')
original=clean[start:end];loop=textwrap.dedent(original)
helper='''// NCBI reference: c++/src/algo/blast/core/link_hsps.c:827-863
// for (H2_index=H_index-1; H2_index>1;) {
//     H2_helper=&lh_helper[H2_index]; H2_index--;
//     if (b0) { H2_index=next_larger; }
//     if (!(b0|b1|b2)) {
//         H_hsp_num=H2->hsp_link.num[index]; H_hsp_sum=H2->hsp_link.sum[index];
//         H_hsp_xsum=H2->hsp_link.xsum[index]; H_hsp_link=H2;
//     }
// }
// A real call per HSP lets Wasm enter optimized scan code during a long group
// invocation. Native retains inlining. The original loop, trace timing, initial
// state and returned values are unchanged; all slices are borrowed, with no allocation.
#[cfg_attr(target_arch = "wasm32", inline(never))]
#[cfg_attr(not(target_arch = "wasm32"), inline(always))]
#[allow(clippy::too_many_arguments)]
fn scan_large_gap_predecessors(
    h_lh_idx: usize,
    pool_lh_helpers: &[LhHelper],
    pool_hsp_links: &[HspLink],
    group_hits: &[UngappedHit],
    h_qe: i32,
    h_se: i32,
    is_target_hsp: bool,
    mut h_sum: i32,
    mut h_num: i16,
    mut h_xsum: f64,
    mut h_link: usize,
) -> (i32, i16, f64, usize) {
'''+textwrap.indent(loop,'    ')+'''
    (h_sum, h_num, h_xsum, h_link)
}

'''
replacement='''                        // NCBI reference: c++/src/algo/blast/core/link_hsps.c:812-863
                        // H_hsp_sum=H2->hsp_link.sum[index]-1;
                        // for (H2_index=H_index-1; H2_index>1;) { ... }
                        // Enter once per original HSP, after its initial-best setup.
                        (h_sum, h_num, h_xsum, h_link) = scan_large_gap_predecessors(
                            h_lh_idx, pool_lh_helpers, pool_hsp_links, &group_hits,
                            h_qe, h_se, is_target_hsp, h_sum, h_num, h_xsum, h_link,
                        );'''
clean=clean[:start]+replacement+clean[end:];anchor='fn link_hsp_group_ncbi(';pos=clean.index(anchor)
# Insert before the existing documentation/attributes of the group function.
marker='// ---------------------------------------------------------------------------\n// Main linking algorithm';
# Exact function insertion at module scope; original doc comment remains with
# the original function by inserting before its nearest contiguous reference block.
pos=clean.rfind('// NCBI',0,pos)
# Use the end of replay function as an unambiguous standalone insertion point.
anchor='fn link_hsp_group_ncbi(';pos=clean.index(anchor)
# The previous function's closing brace is before the group comments.
insert=clean.rfind('\n}\n',0,pos)+3
assert insert>0
clean=clean[:insert]+'\n'+helper+clean[insert:];p.write_text(clean);subprocess.run(['rustfmt',str(p)],check=True)
formatted=p.read_text();aa=raw.decode().splitlines(True);zz=formatted.splitlines(True);out=[]
for tag,i,j,k,l in difflib.SequenceMatcher(None,[x.rstrip('\r\n') for x in aa],[x.rstrip('\r\n') for x in zz],autojunk=False).get_opcodes():out.extend(aa[i:j] if tag=='equal' else zz[k:l] if tag in ['replace','insert'] else [])
p.write_bytes(''.join(out).encode())
(b/'I1-vs-baseline.diff').write_text(''.join(difflib.unified_diff(raw.decode().replace('\r\n','\n').splitlines(True),formatted.splitlines(True),fromfile='baseline/linking.rs',tofile='I1/linking.rs')))
(b/'I1-source-manifest.json').write_text(json.dumps({str(p.relative_to(dest)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.rglob('*') if p.is_file()},indent=2)+'\n')
print('Persistent baseline and I1 source snapshots ready',flush=True)
