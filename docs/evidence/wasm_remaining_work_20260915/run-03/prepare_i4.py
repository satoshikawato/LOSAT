# NCBI reference: c++/src/algo/blast/core/link_hsps.c:812-895
# H_hsp_sum=H2->hsp_link.sum[index]-1; for(H2_index=H_index-1;H2_index>1;) {...}
# Preserve one selection body; restore the indexed driver directly in the caller.
from pathlib import Path
import shutil, re, textwrap, subprocess, difflib, hashlib, json
e = Path(__file__).resolve().parent
source = e / 'work/I3-source'
dest = e / 'work/I4-source'
shutil.copytree(source, dest)
p = dest / 'LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs'
raw = p.read_bytes()
s = raw.decode().replace('\r\n', '\n')
fn_start = s.index('fn scan_large_gap_predecessors(')
macro_start = s.index('    macro_rules! visit_predecessor {', fn_start)
macro_end = s.index('\n\n    #[cfg(all(losat_wasi_threads', macro_start)
old_macro = s[macro_start:macro_end]
step = old_macro[old_macro.index('            let sum = helper.sum[1];'):old_macro.index('\n        }};')]
step = textwrap.dedent(step)
step = step.replace('$cursor = $next_cursor;', 'visit_large_gap_predecessor!(@advance $mode, $cursor, $next_cursor);')
step = step.replace('visit_predecessor!', 'visit_large_gap_predecessor!')
names = ['current_idx', 'helper', 'h_sum', 'h_num', 'h_xsum', 'h_link',
         'is_target_hsp', 'h_qe', 'h_se', 'pool_hsp_links', 'group_hits']
# Rust trace string literals must stay byte-identical; replace identifiers only.
pieces = re.split(r'("(?:\\.|[^"\\])*")', step)
pattern = r'\b(' + '|'.join(names) + r')\b'
for i in range(0, len(pieces), 2):
    pieces[i] = re.sub(pattern, lambda m: '$' + m[0], pieces[i])
step = ''.join(pieces)
macro = '''// NCBI reference: c++/src/algo/blast/core/link_hsps.c:827-863
// H2_helper=&lh_helper[H2_index]; H2_index--;
// if(b0) { H2_index=next_larger; }
// if(!(b0|b1|b2)) { H_hsp_num=H2->hsp_link.num[index];
//   H_hsp_sum=H2->hsp_link.sum[index]; H_hsp_xsum=H2->hsp_link.xsum[index];
//   H_hsp_link=H2; }
// Both cursor drivers expand this single selection/trace body. All caller
// state is explicit; continue targets the surrounding predecessor loop.
macro_rules! visit_large_gap_predecessor {
    (@advance indexed, $cursor:ident, $next_cursor:expr) => { $cursor -= 1; };
    (@advance prefix, $cursor:ident, $next_cursor:expr) => { $cursor = $next_cursor; };
    (@jump indexed, $cursor:ident, $next:expr) => { $cursor = $next; };
    (@jump prefix, $cursor:ident, $next:expr) => {
        $cursor = &$cursor[..$next.saturating_sub(1)];
    };
    ($mode:ident, $current_idx:ident, $helper:ident, $cursor:ident, $next_cursor:expr;
     $h_sum:ident, $h_num:ident, $h_xsum:ident, $h_link:ident,
     $is_target_hsp:ident, $h_qe:ident, $h_se:ident,
     $pool_hsp_links:ident, $group_hits:ident) => {{
'''+textwrap.indent(step, '        ')+'''
    }};
}

'''
args = 'h_sum, h_num, h_xsum, h_link, is_target_hsp, h_qe, h_se, pool_hsp_links, group_hits'
indexed = '''let mut j_lh_idx = h_lh_idx - 1;
while j_lh_idx > 1 {
    let current_idx = j_lh_idx;
    let helper = &pool_lh_helpers[current_idx];
    visit_large_gap_predecessor!(indexed, current_idx, helper, j_lh_idx, j_lh_idx - 1;
        '''+args+''');
}'''
s = s[:macro_start] + s[macro_end:]
s = s.replace('visit_predecessor!(prefix, helper, current_idx, prefix, tail);',
              'visit_large_gap_predecessor!(prefix, current_idx, helper, prefix, tail;\n                '+args+');', 1)
start = s.index('        let mut j_lh_idx = h_lh_idx - 1;', fn_start)
end = s.index('\n    }\n    (h_sum, h_num, h_xsum, h_link)', start)
s = s[:start] + textwrap.indent(indexed, '        ') + s[end:]
insert = s.rfind('// NCBI reference:', 0, fn_start)
# The helper's existing source snippet is immediately before its attributes.
insert = s.index('// NCBI reference: c++/src/algo/blast/core/link_hsps.c:827-863\n// for (H2_index=H_index-1;', 0, fn_start)
s = s[:insert] + macro + s[insert:]
needle = '#[cfg_attr(all(losat_wasi_threads, feature = "wasm-threads"), inline(never))]'
assert s.count(needle) == 1
s = s.replace(needle, '#[cfg(any(test, all(losat_wasi_threads, feature = "wasm-threads")))]\n'+needle)
s = s.replace('// threaded-WASI group invocation. Native and plain WASI inline the indexed scan.',
              '// threaded-WASI group invocation. Native tests can also check the indexed driver.')
call_start = s.index('                        (h_sum, h_num, h_xsum, h_link) = scan_large_gap_predecessors(')
call_end = s.index('                        );', call_start)+len('                        );')
call = textwrap.dedent(s[call_start:call_end])
replacement = '''#[cfg(all(losat_wasi_threads, feature = "wasm-threads"))]
{
'''+textwrap.indent(call, '    ')+'''
}
// NCBI reference: c++/src/algo/blast/core/link_hsps.c:827-863
// for(H2_index=H_index-1; H2_index>1;) { H2_index--; ... }
// Keep the indexed scan in the caller without a function/tuple boundary.
#[cfg(not(all(losat_wasi_threads, feature = "wasm-threads")))]
{
'''+textwrap.indent(indexed, '    ')+'''
}'''
s = s[:call_start] + textwrap.indent(replacement, '                        ') + s[call_end:]
p.write_text(s)
subprocess.run(['rustfmt', str(p)], check=True)
formatted = p.read_text()
# Preserve unchanged original line endings to keep the review diff focused.
aa, zz = raw.decode().splitlines(keepends=True), formatted.splitlines(keepends=True)
output = []
for tag, i, j, k, l in difflib.SequenceMatcher(None,
        [x.rstrip('\r\n') for x in aa], [x.rstrip('\r\n') for x in zz],
        autojunk=False).get_opcodes():
    output.extend(aa[i:j] if tag=='equal' else zz[k:l] if tag in ['replace','insert'] else [])
p.write_bytes(''.join(output).encode())
(e/'I4-vs-I3.diff').write_text(''.join(difflib.unified_diff(
    raw.decode().replace('\r\n','\n').splitlines(keepends=True),
    formatted.splitlines(keepends=True), fromfile='I3/linking.rs',tofile='I4/linking.rs')))
(e/'I4-source-manifest.json').write_text(json.dumps({
    str(f.relative_to(dest)):hashlib.sha256(f.read_bytes()).hexdigest()
    for f in dest.rglob('*') if f.is_file()},indent=2)+'\n')
body = e/'work/I4body-source'
shutil.copytree(dest, body)
shutil.copy2(e/'work/I3body-source/LOSAT/src/main.rs',body/'LOSAT/src/main.rs')
print('I4 source and matched body source prepared',flush=True)
