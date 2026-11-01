# Inputs of the drop-in corpus

Used by `LOSAT/tests/dropin_corpus.py` and `LOSAT/tests/dropin_corpus_manifest.tsv`: the command lines that downstream
tools run and the common forms of each program. `dropin_corpus.py` links every file here into the directory of each run. The rows name them relatively, as a
tool would. A row may write next to the links but must not write through one; a run that changes an input differs.

FASTA files:

- `nq.fna`, `ns.fna`: `LOSAT/tests/db_search/nslice_q.fna` and `nslice.fna`.
- `pq.faa`, `tq.fna`: `LOSAT/tests/db_search/tgc4_q.faa` and `tgc4.fna`.
- `pq1.faa`: the first record of `pq.faa`.
- `gapdb.faa`: records 541-570 of the CDS translations (efetch `fasta_cds_aa`) of Salmonella Typhi CT18, AL513382.1;
  four of them start with `-` (rows ORT07 and ORT08).

Lists and maps:

- `n.seqids` (`a_l6k`, `b_p6k`), `t.seqids` (`h_slice`), `p.seqids` (`BDV02158.1`, `BDV02159.1`), `taxids.txt` (`9606`).
- `taxmap.txt`: each id of `ns.fna` with the taxid 9606.

Databases, built by NCBI BLAST+ 2.17.0 `makeblastdb` (version 5, so `blastdbcmd -info` prints their date from the
frozen `.nin`):

    makeblastdb -in ns.fna -dbtype nucl -out ndb
    makeblastdb -in tq.fna -dbtype nucl -out tdb
    makeblastdb -in pq.faa -dbtype prot -out pdb
    makeblastdb -in ns.fna -dbtype nucl -out ndbp -parse_seqids -taxid 9606
    makeblastdb -in tq.fna -dbtype nucl -out tdbp -parse_seqids -taxid 9606
    makeblastdb -in pq.faa -dbtype prot -out pdbp -parse_seqids -taxid 9606
    makeblastdb -in gapdb.faa -dbtype prot -out gapdb

`brig_ref.fna` is a copy of `ns.fna` under the name BRIG gives its database (row N26), built from the parent directory
of these files: `makeblastdb -in fx/brig_ref.fna -dbtype nucl`, so its title holds `fx/`.

To change a row or an input: edit the manifest or the file, then, from `LOSAT/`,
`python3 tests/dropin_corpus.py freeze --ncbi-bin <NCBI bin directory> --only <row id prefix>`.
