# NCBI reference: c++/src/algo/blast/core/link_hsps.c:827-863
# for (H2_index=H_index-1; H2_index>1;) { ... }
# Restrict the real function boundary to threaded WASI; preserve every statement.
from pathlib import Path
import shutil, hashlib, json, difflib

e = Path(__file__).resolve().parent
source = e / 'work/I2b-source'
dest = e / 'work/I3-source'
shutil.copytree(source, dest)
p = dest / 'LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs'
old = p.read_bytes()
before = b'''// A call per HSP provides a repeated real function-entry boundary within a long
// group invocation. Native retains inlining. The original loop, trace timing, initial
// state and returned values are unchanged; all slices are borrowed, with no allocation.
#[cfg_attr(target_arch = "wasm32", inline(never))]
#[cfg_attr(not(target_arch = "wasm32"), inline(always))]'''
after = b'''// A call per HSP provides a repeated real function-entry boundary within a long
// threaded-WASI group invocation. Native and plain WASI inline the indexed scan.
// The loop, trace timing, initial state and returned values are unchanged;
// all slices are borrowed, with no allocation.
#[cfg_attr(all(losat_wasi_threads, feature = "wasm-threads"), inline(never))]
#[cfg_attr(not(all(losat_wasi_threads, feature = "wasm-threads")), inline(always))]'''
assert old.count(before) == 1
p.write_bytes(old.replace(before, after))
(e / 'I3-vs-I2b.diff').write_text(''.join(difflib.unified_diff(
    old.decode().replace("\r\n", "\n").splitlines(keepends=True), p.read_text().splitlines(keepends=True),
    fromfile='I2b/linking.rs', tofile='I3/linking.rs')))
(e / 'I3-source-manifest.json').write_text(json.dumps({
    str(path.relative_to(dest)): hashlib.sha256(path.read_bytes()).hexdigest()
    for path in dest.rglob('*') if path.is_file()
}, indent=2) + '\n')

# NCBI blastn_app.cpp:59-67: m_StopWatch.Start(); m_StopWatch.Elapsed();
# Keep the existing matched measurement-only dispatch clocks byte-identical.
body = e / 'work/I3body-source'
shutil.copytree(dest, body)
shutil.copy2(e / 'work/I2bbody-source/LOSAT/src/main.rs', body / 'LOSAT/src/main.rs')
print('I3 source and matched body source prepared', flush=True)
