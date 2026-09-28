#!/usr/bin/env python3
"""Package approved NCBI comparison receipts without changing raw bytes.

NCBI c++/src/objtools/align_format/tabular.cpp:1100-1108:
// void CBlastTabularInfo::Print() {
//     ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
//         x_PrintField(*iter);
//     }
//     m_Ostream << "\n";
// }
"""
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile

HERE = Path(__file__).resolve().parent
PRIOR = Path('/mnt/c/Users/genom/GitHub/LOSAT/docs/evidence/losatx_stage_g_hold_completion_clarified/run-20260928T104629Z')
SOURCES = ('windows_authority_proposal.json', 'windows_m1_n4_authority_proposal.json', 'windows_pr5_oracle_authority_proposal.json')

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    index = {}
    objects = {}
    for name in SOURCES:
        proposal = json.loads((PRIOR / name).read_text())
        for case in proposal['cases']:
            pins = [case['process_record'], case['raw_stdout'], case['raw_stderr']]
            if 'raw_report' in case:
                pins.append(case['raw_report'])
            for item in case['inputs'].values():
                pins.append({'path': item.get('source_path', item.get('mirror')), 'sha256': item['sha256']})
            for pin in pins:
                path = Path(pin['path'])
                digest = pin['sha256']
                assert path.is_file() and sha(path) == digest, path
                if 'size' in pin:
                    assert path.stat().st_size == pin['size']
                index[str(path)] = digest
                objects.setdefault(digest, path)
    path = HERE / 'approved_windows_receipts.tar.gz'
    assert not path.exists()
    with path.open('wb') as file:
        with gzip.GzipFile(fileobj=file, mode='wb', filename='', mtime=0, compresslevel=6) as zipped:
            with tarfile.open(fileobj=zipped, mode='w') as tar:
                for digest, source in sorted(objects.items()):
                    data = source.read_bytes()
                    info = tarfile.TarInfo(f'sha256/{digest}')
                    info.size = len(data)
                    info.mode = 0o444
                    info.mtime = info.uid = info.gid = 0
                    tar.addfile(info, io.BytesIO(data))
    (HERE / 'receipt_index.json').write_text(json.dumps({
        'schema': 'losatx-stage-g-approved-windows-receipt-index-v1',
        'source_path_to_sha256': index,
        'unique_objects': len(objects),
        'source_files': len(index),
        'archive_sha256': sha(path),
        'raw_bytes_unchanged': True,
    }, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'unique_objects': len(objects), 'files': len(index), 'archive_bytes': path.stat().st_size, 'archive_sha256': sha(path)}))

if __name__ == '__main__':
    main()
