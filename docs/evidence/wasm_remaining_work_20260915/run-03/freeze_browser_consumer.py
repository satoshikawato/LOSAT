# NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:173-188
# Run(); Join(&result);
# Freeze the actual consumer host already exercised by the browser gates.
from pathlib import Path
import json,hashlib
e=Path(__file__).resolve().parent;web=Path('/mnt/c/users/genom/github/gbdraw/gbdraw/web');out=e/'browser-consumer-source';out.mkdir()
sha=lambda b:hashlib.sha256(b).hexdigest()
expected={}
for name in ['I2b-browser-lifecycle-threaded','I2b-browser-lifecycle-nonisolated','I2b-browser-blastp','baseline-browser-blastp']:
    report=json.loads((e/name/'manifest.json').read_text());assert report['status']=='PASS'
    for filename,digest in report['served_sha256'].items():
        source=Path(filename)
        if not source.is_relative_to(web) or filename.endswith(':diagnostic-overlay'):continue
        assert filename not in expected or expected[filename]==digest,(filename,'consumer changed between gates')
        expected[filename]=digest
expected[str(web/'index.html')]=sha((web/'index.html').read_bytes())
manifest={}
for filename,digest in expected.items():
    source=Path(filename);data=source.read_bytes();assert sha(data)==digest,(filename,'source changed since successful gates')
    relative=source.relative_to(web);target=out/relative;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data);manifest[str(relative)]=digest
(e/'browser-consumer-manifest.json').write_text(json.dumps(dict(original_root=str(web),snapshot_root=str(out),files=manifest),indent=2)+'\n')
print('Frozen actual consumer files:',len(manifest))
