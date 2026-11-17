# A PSSM whose values nest deeper than NCBI's stack holds (PD-LOSAT-NCBI-DEFECTS exception 14)

Line RG item D3-1, fix round 7, 2026-11-17. NCBI BLAST+ 2.17.0 (the official Linux x86-64 build) and LOSAT,
psiblast and tblastn `-in_pssm`.

## What NCBI does

NCBI's reader of a PSSM in ASN.1 text (`CObjectIStreamAsn`) reads a value whose type holds values of its own
type by recursion on the main thread's stack: `ReadObject` and `SkipObject` of each value
(`c++/include/serial/impl/objistr.inl:93-96`, `:99`), entered from each member and variant
(`c++/src/serial/member.cpp:710`, `c++/src/serial/variant.cpp:471`), with the class, choice and container
functions between them. A PSSM query whose values nest some
thousands of levels deep fills the stack, and NCBI ends with SIGSEGV: exit status 139, nothing written.

The depth at which it crashes is not the input's. Psiblast, with the default stack limit of 8 MiB, the most room
it gets (a minimal environment, the stack's random start off with `setarch -R`), crashes from these depths
(values of the type open at once; line RG D3-1: `uobjs`, `vcons`, `ufield` and `mdelta` bisected by depth in
review 4e of fix round 6, the others from the bytes per level and base of each path bisected to 16 bytes in fix
round 6, which that review checked to 0 bytes):

| nest (`gen_nests.py` kind) | crashes from |
|---|---|
| User-objects through `data objects` (`uobjs`) | 2,976 |
| Variation-refs through `data set { variations }` (`variation`) | 3,967 |
| User-objects through `data object` (`user`) | 4,156 |
| OrgNames through `name hybrid` (`orgname`) | 4,514 |
| Variation-refs through `consequence { variation }` (`vcons`) | 4,760 |
| Bioseq-sets as the query entry (`set`) | 5,082 |
| Seq-aligns through `segs disc` (`align`) | 5,752 |
| User-fields through `data fields` (`ufield`) | 5,818 |
| Seq-loc mixes or equivs in a feature's location (`mix`, `equiv`) | 6,888 |
| Pubs through `equiv` (`pub`) | 7,587 |
| SeqTable-multi-data through `int-scaled` (`multidata`) | 9,519 |
| SeqTable-multi-data through `int-delta` (`mdelta`) | 9,703 |

Each run's random start takes 0 to 8,191 more bytes, a few levels; the command line and the environment, which the
kernel lays out at the stack's end, take their length; a nest of several types takes the bytes of each, and the
same type takes other bytes along another path; a larger stack limit holds more. So one PSSM can be read in one run
and crash in the next. With a stack large enough NCBI reads every such PSSM and writes its normal output: the
search, or for nested Bioseq-sets (a query Seq-entry set) its `CSeq_entry::GetSeq()` exception (exit status 255).

## The runs

`results.tsv`: for each nest (kind and depths of `gen_nests.py`) and program, `rc/stdout bytes/stderr bytes` of
NCBI with an 8 MiB stack, of NCBI with no stack limit and of LOSAT with an 8 MiB stack, and whether LOSAT's stdout,
stderr and exit status are NCBI's of the run with no stack limit, and NCBI's with 8 MiB are too.

- NCBI with 8 MiB crashes (139, no output) on every input deeper than the table above (32 of 36 runs).
  `vcons 4500` and `mdelta 9600`, shallower than their crash depths, read with 8 MiB too.
- NCBI with no stack limit reads every input: stdout 207 bytes for psiblast (SHA-256 `b0a21fccc17fe13c...`, its
  stderr 394 bytes, the composition-based statistics warning) and 352 bytes for tblastn (`bce5894d38192d54...`);
  for the Bioseq-set nests no stdout and NCBI's GetSeq() exception (358 bytes, `f7689d2f6002b95a...`, exit status
  255). 100,000 nested Seq-loc mixes included.
- LOSAT writes NCBI's bytes of the run with no stack limit in all 36 runs.

Commands (from `LOSAT/`, each NCBI run under a lock of its own, `BLAST_USAGE_REPORT=false`, `HOME` an empty
directory; NCBI's runs with an 8 MiB stack load a shim that ends a SIGSEGV with exit status 139 and no core dump):

```sh
python3 -I ../docs/evidence/pssm_nest_depth/gen_nests.py tests/fixtures/pssm/d3om_r6_nest_mix_6600.asn.gz mix 7000 > in.asn
env -i PATH=/usr/bin:/bin HOME=$EMPTY BLAST_USAGE_REPORT=false prlimit --as=4294967296 --core=0 --stack=8388608 \
  timeout 120 $NCBI/psiblast -in_pssm in.asn -db tests/db_search/paja -outfmt 6 -num_threads 1      # 139
env -i PATH=/usr/bin:/bin HOME=$EMPTY BLAST_USAGE_REPORT=false prlimit --as=4294967296 --core=0 --stack=unlimited \
  timeout 120 $NCBI/psiblast -in_pssm in.asn -db tests/db_search/paja -outfmt 6 -num_threads 1      # 0
env -i PATH=/usr/bin:/bin HOME=$EMPTY prlimit --as=4294967296 --core=0 --stack=8388608 \
  timeout 120 losat psiblast -in_pssm in.asn -db tests/db_search/paja -outfmt 6 -num_threads 1      # 0, NCBI's bytes
```

tblastn: `-db tests/db_search/tslice`. The depths of `results.tsv` are the arguments of `gen_nests.py` (its usage
says how many values each gives).

## The fixture rows

`gen_nests.py` makes the inputs of the rows `db.ps.inpssm.d3.om.r6.nest_*` and `db.ps.inpssm.d3.om.r7.nest_*` of
`LOSAT/tests/db_search_manifest/psiblast.tsv` and their `db.tbn.` rows of `tblastn.tsv`
(`LOSAT/tests/fixtures/pssm/d3om_r6_nest_*.asn.gz`, `d3om_r7_nest_*.asn.gz`, byte for byte; for example
`d3om_r7_nest_uobjs_3500` is `uobjs 3500`, `d3om_r6_nest_user_4000` is `user 3999` and
`d3om_r6_nest_sets_users_2537_2501` is `sets_users 2537 2500`). The rows marked `ncbi-stack: unlimited` are frozen
from NCBI's run with no stack limit (`LOSAT/tests/db_search_fixtures.py freeze`); the others NCBI reads with its
default stack.
