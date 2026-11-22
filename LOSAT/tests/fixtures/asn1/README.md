# ASN.1 text files and their binary twins

The generic reader's tests (`src/report/serial/asn_types/read_ber_tests.rs`) read these: each `<name>.asn` is an ASN.1 text
document of a type of the reader's table (`Seq-entry`, `Bioseq-set`, `Bioseq`) and `<name>.asnb` its binary (BER) twin,
written by `datatool` 2.24.0 of the bioconda BLAST 2.17.0 package (`datatool -m all.asn -v <name>.asn -e <name>.asnb`;
`all.asn` is the `general`, `biblio`, `medline`, `pub`, `seq`, `seqloc`, `seqfeat`, `seqalign`, `seqres`, `seqblock`, `seqcode`,
`seqset`, `seqtable`, `variation` and `blastdb` modules of the NCBI C++ Toolkit 2.17.0, one after the other), so every length is
indefinite. `multi.asn` holds three `Seq-entry` documents, `multi.asnb` their three objects. `trunc.asn` (a text cut inside its
object), `e1.asn` (with its twin `e1.asnb`) and `f.fa` (a FASTA file) are the inputs of NCBI BLAST+ 2.17.0 makeblastdb probes
whose errors the tests compare (`-input_type asn1_txt` / `asn1_bin`). Made for line RG, D3c (probe directory
`remaining-gaps-20261112/D/d3/c/probe/asn2/in`, survey `D/d3/c/SURVEY_D3c.md` §3), copied for D3c-3a.

`ber/<case>.asnb` (20 files): the inputs of NCBI BLAST+ 2.17.0 makeblastdb probes of the binary reader's rules
(`makeblastdb -in input -input_type asn1_bin -dbtype nucl -out db -title T`, official build; line RG, D3c-3a, probe
directory `remaining-gaps-20261112/D/d3/c3a/probe/ber/<case>` with each run's status, standard error and description,
summary `SUMMARY.tsv`). `base` is datatool's binary of a small Seq-entry (`e1` with a user descriptor holding a REAL, an
INTEGER and a BOOLEAN; `D/d3/c3a/probe/ber_src/e1r.asn`); the others change its bytes (generator
`D/d3/c3a/tools/gen_seqentry.py`): definite lengths (`def*`, `len81*`), long tags (`tag_long`, `tag_app_int`), cuts
(`cut_1`, `cut_last`), members out of place (`member9`, `order`, `missing_inst`) and string bytes (`title_*`,
`store_ff`).
