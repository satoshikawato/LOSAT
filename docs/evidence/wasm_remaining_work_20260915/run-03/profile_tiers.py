# NCBI reference: c++/src/algo/blast/core/link_hsps.c:827-863
# for (H2_index=H_index-1; H2_index>1;) { ... }
# Diagnostic PC samples/code maps only; never use these runs as speed evidence.
from pathlib import Path
import subprocess, json
e = Path(__file__).resolve().parent
for version, idx in [('baseline', 1340), ('I1', 1368)]:
    d = e / 'work' / (version + '-profile-logs')
    d.mkdir()
    flags = ['--prof', '--logfile=' + str(d / 'v8.log'), '--print-wasm-code-function-index=' + str(idx)]
    cmd = ['python3', str(e / 'run_gate.py'), version, '--label', version + '-profile', '--cases', 'MjeNMV.MelaMJNV.tlosatx', '--kinds', 'threaded', '--threads', '8', '--flags', json.dumps(flags)]
    print(version, 'PROFILE', flush=True)
    subprocess.run(cmd, check=True)
