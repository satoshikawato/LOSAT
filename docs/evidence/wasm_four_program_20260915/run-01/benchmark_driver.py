# NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-48
# const string kArgQuery("query"); const string kArgSubject("subject");
# Add only the predeclared P2 fixture; preserve the existing measurement harness.
import sys
from pathlib import Path
b=Path(__file__).resolve().parent
sys.path.insert(0,'/mnt/c/Users/genom/GitHub/LOSAT/LOSAT/tests')
import benchmark_wasm_threading as h
original=h.select_cases
def select(out,selection):
 cases=[]
 for key in selection:
  if key=='single-query-32-matches':cases.append((key,'blastp',b/'single-query/query.faa',b/'single-query/subjects.faa',[]))
  else:cases.extend(original(out,[key]))
 return cases
h.select_cases=select
h.main()
