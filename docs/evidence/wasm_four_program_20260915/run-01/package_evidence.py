# NCBI reference: c++/src/objtools/align_format/tabular.cpp:1098-1108
# x_PrintField(*iter); m_Ostream << "\n";
# Preserve exact outputs and metadata. Identical bytes may share tar hardlinks;
# all original paths and hashes remain in the manifest.
import hashlib
import json
import os
import tarfile
from pathlib import Path

B=Path('/tmp/losat-four-program-20260915')
OUT=Path('/mnt/c/Users/genom/GitHub/LOSAT/docs/evidence/wasm_four_program_20260915/run-01')

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()

def closed():
    stop=json.loads((B/'measurement-stop.json').read_text())
    assert stop['status']=='STOPPED_BY_USER' and stop['process_absence_verified']
    assert stop['eligible_repeats']==[1,2,3]


def selected():
    files=[]
    for entry in sorted(B.iterdir()):
        if entry.is_file():
            files.append(entry)
            continue
        if not entry.is_dir() or entry.name=='__pycache__':
            continue
        if (entry/'build.json').is_file():
            # Standalone built artifacts and build metadata, not Cargo objects.
            files.extend(p for p in entry.iterdir() if p.is_file())
            for rel in ['native-command/release/LOSAT',
                        'serial-command/wasm32-wasip1/release/LOSAT.wasm',
                        'threaded-command/wasm32-wasip1-threads/release/LOSAT.wasm',
                        'threaded-reactor/wasm32-wasip1-threads/release/LOSAT.wasm']:
                p=entry/rel
                if p.is_file():
                    files.append(p)
            if (entry/'artifacts').is_dir():
                files.extend(p for p in (entry/'artifacts').rglob('*') if p.is_file())
            continue
        for root, dirs, names in os.walk(entry, followlinks=False):
            if 'CACHEDIR.TAG' in names and '.rustc_info.json' in names:
                dirs[:]=[]
                continue
            dirs[:]=[d for d in dirs if d not in {'.git','target','node_modules','__pycache__'}]
            for name in names:
                p=Path(root)/name
                if p.is_file():
                    files.append(p)
    return sorted(set(files))

if __name__=='__main__':
    closed()
    archive=OUT/'execution-evidence.tar.gz'
    assert not archive.exists(), 'Use a new archive name instead of replacing evidence'
    entries=[]
    known={}
    with tarfile.open(archive,'w:gz',compresslevel=1) as tar:
        for i,p in enumerate(selected()):
            before=p.stat()
            sha=digest(p)
            name=str(Path(B.name)/p.relative_to(B))
            info=tar.gettarinfo(str(p),arcname=name)
            info.mode &= 0o7777  # Tar stores permission bits without the stat file type.
            assert info.isfile(), ('Unexpected non-regular evidence file',p)
            same=known.get((sha,info.mode))
            if same:
                info.type=tarfile.LNKTYPE
                info.linkname=same
                info.size=0
                tar.addfile(info)
            else:
                with p.open('rb') as f:
                    tar.addfile(info,f)
                known[(sha,info.mode)]=name
            after=p.stat()
            assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns), ('Evidence changed while archiving',p)
            entries.append(dict(path=name,original_path=str(p),bytes=before.st_size,sha256=sha,mode=info.mode,hardlink=same))
            if i%1000==0:
                print('archived',i,'files',flush=True)
    expected={entry['path']:entry for entry in entries}
    verified={}
    with tarfile.open(archive,'r|gz') as tar:
        for member in tar:
            if member.islnk():
                sha,size=verified[member.linkname]
            else:
                h=hashlib.sha256();size=0
                f=tar.extractfile(member)
                for block in iter(lambda:f.read(1024*1024),b''):
                    h.update(block);size+=len(block)
                sha=h.hexdigest()
            assert (sha,size,member.mode)==(expected[member.name]['sha256'],expected[member.name]['bytes'],expected[member.name]['mode']), member.name
            verified[member.name]=(sha,size)
    assert set(verified)==set(expected)
    manifest=dict(archive=archive.name,sha256=digest(archive),bytes=archive.stat().st_size,
                  files=len(entries),unique_contents=len(known),
                  archive_content_verification='PASS: every regular member and hardlink matches original SHA256, size and mode',
                  excludes='Cargo intermediate/object caches, .git, node_modules, Python caches; standalone built artifacts are included',
                  entries=entries)
    (OUT/'execution-evidence-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:v for k,v in manifest.items() if k!='entries'},indent=2),flush=True)
