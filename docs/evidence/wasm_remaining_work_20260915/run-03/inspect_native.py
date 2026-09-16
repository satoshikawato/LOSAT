# NCBI reference: c++/src/algo/blast/core/link_hsps.c:827-863
# for (H2_index=H_index-1; H2_index>1;) { ... }
# Read-only ELF comparison; never rewrite a built artifact during inspection.
from pathlib import Path
import struct, hashlib, json
e = Path(__file__).resolve().parent
def elf(path):
    raw = path.read_bytes()
    assert raw[:6] == b'\x7fELF\x02\x01'
    table, width, count, string_index = struct.unpack_from('<Q', raw, 40)[0], *struct.unpack_from('<HHH', raw, 58)
    headers = [struct.unpack_from('<IIQQQQIIQQ', raw, table + i * width) for i in range(count)]
    h = headers[string_index]
    names = raw[h[4]:h[4]+h[5]]
    sections = {}
    for h in headers:
        name = names[h[0]:].split(b'\0',1)[0].decode()
        if h[1] == 8: continue  # NOBITS has no on-disk section contents.
        sections[name] = raw[h[4]:h[4]+h[5]]
    symbols = []
    if '.symtab' in sections:
        strings = sections['.strtab']
        for offset in range(0, len(sections['.symtab']), 24):
            sn, info, other, shndx, value, size = struct.unpack_from('<IBBHQQ', sections['.symtab'], offset)
            name = strings[sn:].split(b'\0',1)[0].decode()
            if 'link_hsp_group_ncbi' in name or 'scan_large_gap_predecessors' in name:
                symbols.append(dict(name=name,value=value,size=size,section_index=shndx))
    return raw, sections, symbols
base, aa, asym = elf(e/'work/baseline/artifacts/losat-native-command')
cand, bb, bsym = elf(e/'work/I1/artifacts/losat-native-command')
report = dict(baseline_sha256=hashlib.sha256(base).hexdigest(), candidate_sha256=hashlib.sha256(cand).hexdigest(), sections={}, baseline_symbols=asym, candidate_symbols=bsym)
for name in sorted(aa.keys() | bb.keys()):
    a, b = aa.get(name,b''), bb.get(name,b'')
    report['sections'][name] = dict(baseline_bytes=len(a), candidate_bytes=len(b), equal=a==b, baseline_sha256=hashlib.sha256(a).hexdigest(), candidate_sha256=hashlib.sha256(b).hexdigest(), differing_byte_count=sum(x!=y for x,y in zip(a,b))+abs(len(a)-len(b)))
(e/'native-emission-comparison.json').write_text(json.dumps(report,indent=2)+'\n')
for name in ['.text','.rodata']:
    print(name,report['sections'][name])
print('baseline symbols',asym)
print('candidate symbols',bsym)
