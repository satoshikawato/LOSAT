# NCBI reference: c++/src/objtools/align_format/tabular.cpp:1098-1108
# x_PrintField(*iter); m_Ostream << "\n";
# Evidence-only integrity check: preserve and verify the recorded raw bytes.
from pathlib import Path
import hashlib,json,tarfile
e=Path(__file__).resolve().parent;m=json.loads((e/'bundle-manifest.json').read_text());results=[]
for name,record in m.items():
 p=e/name;assert hashlib.sha256(p.read_bytes()).hexdigest()==record['sha256'];expected={x['path']:x for x in record['members']};seen=set()
 with tarfile.open(p,'r:gz')as tar:
  for member in tar:
   assert member.isfile()and member.name in expected and member.name not in seen;data=tar.extractfile(member).read();r=expected[member.name];assert len(data)==r['bytes'] and hashlib.sha256(data).hexdigest()==r['sha256'],member.name;seen.add(member.name)
 assert seen==set(expected);results.append(dict(archive=name,sha256=record['sha256'],members=len(seen),status='PASS'));print(name,'PASS',len(seen),flush=True)
(e/'archive-verification.json').write_text(json.dumps(results,indent=2)+'\n')
