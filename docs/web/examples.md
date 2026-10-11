# LOSAT Web examples

LOSAT Web ships one small example per program (S16, instruction item 6): a query file and a
subject file that "Try an example" puts in the search form, and that the optional browser
self-test runs through the engine. They are cut from public INSDC/RefSeq records whose copies
are in the repository (`LOSAT/tests/fasta/`). They are **release assets, not research data**:
the build emits them as hashed files of the release (`web/app/src/infra/browser/examples.ts`
imports them with Vite's `?url&no-inline`), the page fetches them from the site, and nothing
is generated in the page. Two example files with the same bytes are one asset of the build.

- Files: `web/app/src/examples/example-<program>-<query|subject>.fasta`; the search form names
  them the same way (`domain/examples.ts` `exampleFileName`), and the argv of an example is
  `<program> -query example-<program>-query.fasta -subject example-<program>-subject.fasta`
  with default settings.
- Made by: `python3 web/app/tools/examples.py cut` (the table of regions and records is in the
  script; `check` confirms that the files are what `cut` writes).
- Definition lines are those of the source records. A record cut to a region says so after its
  definition line: `... complete genome [region 1-2000]`. Protein records are whole, as the
  source file has them. Line width: as in the source (70 residues for the genomes, 60 for the
  proteins).

## Why these records

The five examples tell one story: the chromosomal replication initiator **DnaA** and its gene
**dnaA**, at the start of four *Candidatus* Hepatoplasma genomes and one *Candidatus*
Tyloplasma genome. Each genome begins with dnaA (within its first 1,500 nt), so its first
2,000 nt hold the whole gene. The proteins are about 455 residues. Every search therefore finds
a few hits with alignments, at graded similarity, and finishes in milliseconds (native) and well
under a second in a browser:

| Program | Query | Subject | Hits (native LOSAT, default settings) |
| --- | --- | --- | --- |
| BLASTN | dnaA region of *Ca.* H. crinochetorum Tokyo2021 | dnaA regions of *Ca.* H. crinochetorum Av, *Ca.* H. scabrum, *Ca.* H. vulgare | 3 HSPs, one per subject: 99.2 %, 85.3 % and 71.8 % identity (megablast) |
| BLASTP | DnaA of *Ca.* H. vulgare | DnaA of *Ca.* H. crinochetorum, *Ca.* H. scabrum, *Ca.* T. litorale | 3 HSPs, 41 % to 47 % identity |
| BLASTX | dnaA region of *Ca.* H. vulgare | the BLASTP subjects | 5 HSPs on the 3 subjects |
| TBLASTN | DnaA of *Ca.* H. vulgare | dnaA regions of *Ca.* H. crinochetorum Tokyo2021, *Ca.* H. scabrum, *Ca.* T. litorale | 5 HSPs on the 3 subjects |
| TBLASTX | dnaA region of *Ca.* H. vulgare | the TBLASTN subjects | 103 HSPs on the 3 subjects (several reading frames) |

BLASTX cannot be searched in LOSAT Web until session SX (`domain/programs.ts` `unavailable`):
its example is in the release already, so SX only enables the program; until then "Try an
example" is disabled for BLASTX with that reason, and the self-test skips it.

## The files

| File | Bytes | Records (source file in `LOSAT/tests/fasta/`) | SHA-256 |
| --- | ---: | --- | --- |
| `example-blastn-query.fasta` | 2,131 | AP027132.1 region 1-2000 (`AP027132.fasta`) | `c0b5cb1cee7af3f4d4ccbc4d1a7c44d9a5abf5629d64da1eff94cee68ddd2cdd` |
| `example-blastn-subject.fasta` | 6,371 | NZ_CP006932.1 region 1-2000 (`NZ_CP006932.fasta`), AP027133.1 region 1-2000 (`AP027133.fasta`), AP027131.1 region 1-2000 (`AP027131.fasta`) | `84779a429b069e2ca3a6e1ad7d85ee9d531bfbbefa8a5986f0c4478dd6b8e4f6` |
| `example-blastp-query.fasta` | 556 | BDV02112.1 (`AP027131.faa`) | `d7161be648b6c156612e102e329ab63f1435298999fdb29e8f3f0a90ac1264de` |
| `example-blastp-subject.fasta` | 1,687 | BDV02707.1 (`AP027132.faa`), BDV03280.1 (`AP027133.faa`), BDU67348.1 (`AP027078.faa`) | `5eae832d66f72a29679eee66df17dfa22697c35e61eae729610f9e2ac1f7f286` |
| `example-blastx-query.fasta` | 2,121 | AP027131.1 region 1-2000 (`AP027131.fasta`) | `71d308056d995a02ffd33b6c2deac715b600008a5ed477e360fd6138eaf2aadd` |
| `example-blastx-subject.fasta` | 1,687 | as `example-blastp-subject.fasta` | `5eae832d66f72a29679eee66df17dfa22697c35e61eae729610f9e2ac1f7f286` |
| `example-tblastn-query.fasta` | 556 | as `example-blastp-query.fasta` | `d7161be648b6c156612e102e329ab63f1435298999fdb29e8f3f0a90ac1264de` |
| `example-tblastn-subject.fasta` | 6,378 | AP027132.1 region 1-2000 (`AP027132.fasta`), AP027133.1 region 1-2000 (`AP027133.fasta`), AP027078.1 region 1-2000 (`AP027078.fasta`) | `06a23018ab9d082397f3a212d245e1b1d488b2d4d06582fb6dc8aaf8d9878e7a` |
| `example-tblastx-query.fasta` | 2,121 | as `example-blastx-query.fasta` | `71d308056d995a02ffd33b6c2deac715b600008a5ed477e360fd6138eaf2aadd` |
| `example-tblastx-subject.fasta` | 6,378 | as `example-tblastn-subject.fasta` | `06a23018ab9d082397f3a212d245e1b1d488b2d4d06582fb6dc8aaf8d9878e7a` |

The source records (accession and version; the protein records are the CDS translations that the
genome records annotate):

| Accession | Record | Length |
| --- | --- | ---: |
| AP027131.1 | *Candidatus* Hepatoplasma vulgare Av-JP DNA, complete genome (MAG) | 662,108 nt |
| AP027132.1 | *Candidatus* Hepatoplasma crinochetorum Tokyo2021 DNA, complete genome (MAG) | 643,039 nt |
| AP027133.1 | *Candidatus* Hepatoplasma scabrum Ps-JP DNA, complete genome (MAG) | 606,194 nt |
| AP027078.1 | *Candidatus* Tyloplasma litorale Fukuoka2020 DNA, complete genome (MAG) | 615,622 nt |
| NZ_CP006932.1 | *Candidatus* Hepatoplasma crinochetorum Av chromosome, complete genome (RefSeq) | 657,101 nt |
| BDV02112.1 | chromosomal replication initiator protein DnaA [*Ca.* Hepatoplasma vulgare] | 455 aa |
| BDV02707.1 | chromosomal replication initiator protein DnaA [*Ca.* Hepatoplasma crinochetorum] | 455 aa |
| BDV03280.1 | chromosomal replication initiator protein DnaA [*Ca.* Hepatoplasma scabrum] | 454 aa |
| BDU67348.1 | chromosomal replication initiator protein DnaA [*Ca.* Tyloplasma litorale] | 470 aa |

## The self-test's expected outputs

`web/app/src/examples/self-test.json` gives, for each example, the argv, the SHA-256 of both
input files, and the SHA-256 and length of outfmt 0, 6 and 7 that the **native LOSAT CLI of the
engine commit the site's reactors were built from** writes for that argv (`-outfmt` 0, 6, 7),
with whether **NCBI BLAST+ 2.17.0** writes the same bytes for the same argv. It names no machine
path (plan TD-11). The browser self-test (`web/app/src/application/self-test.ts`, "Test this
browser (optional)" in the Getting started panel) runs the examples of the programs that LOSAT
Web searches through the application's engine and compares with it.

Made on 2026-10-11 from engine commit `141615ec` (native CLI `losat 0.2.0`, SHA-256
`0236291b81ec8e1735c848c49a4e36c1ae0af3742cc15832c0f9b9352edbb928`; the reactors of
`$BUILD_ROOT/s14-reactors`): all 15 outputs (5 programs x 3 formats) are byte for byte NCBI
BLAST+ 2.17.0's. `tests/unit/self-test-engine.test.ts` (with `LOSAT_WEB_REACTORS`) runs the four
searchable examples through the serial reactor in Node and finds the same SHA-256.

Remake the table whenever the engine changes (the reactor test fails until it agrees again),
under the oracle lock (`.claude/skills/losat-oracle-runs/SKILL.md`; 30 tiny searches):

```bash
flock "$BUILD_ROOT/oracle.lock" python3 web/app/tools/examples.py expected \
  --losat <native LOSAT built from the engine commit> --engine-commit <commit> --ncbi-bin "$NCBI_BIN"
```
