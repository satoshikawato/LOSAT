# Inputs of the `-html` fixtures

Used by the `html.*` rows of `LOSAT/tests/outfmt0_manifest.tsv` (`-subject`, `-outfmt 0` and `18` with `-html`, line RG
item T4a2). The oracle is NCBI BLAST+ 2.17.0 (`run_oracle.py --freeze --only html.`). Every file is the input of the
line's NCBI probe (survey T4, `SURVEY_T4.md` §1.2-§1.6 and §1.9), copied byte for byte; the row's `covers` cell names
the probe case it repeats. Three more inputs of those rows are already in the repository and byte-identical to the
probe's: `tests/db_search/of18/pq.faa`, `nq.fna` and `ns_ids.fna` (survey T2's `gen.py`, seed 20261113).

Special titles and ids (survey T4's `gen.py`, from survey T2's `pq.faa`, `ps_plain.faa`, `nq.fna`, `ns_plain.fna`):
- `pq_spec.faa`, `nq_spec.fna`: a query whose title holds `< > & " '`, `&amp; &lt;x&gt; &#62;`, brackets, `;`, `|`, `/`
  and wraps the `Query=` line more than once, and a query without hits.
- `ps_spec.faa`, `ns_spec.fna`: four subjects whose titles hold `<b>`, `&`, `"`, entities and `&auml;` (decoded, then
  encoded), `&` and `<` at the description cut and the 60-column heading wrap, a word of 70 letters, a plain title.
- `pq_ids.faa`: `pq.faa`'s queries as `lcl|pq1`, `gi|555|ref|NP_000001.1|` and `lcl|pq2` (`-parse_deflines` queries).
- `pq_wrap.faa`: the 13 wrap titles w1-w13 (punctuation runs, double spaces, a 75-letter word, `(`, `&` x 80, `& < > "`,
  tabs, UTF-8 `ü`, commas, hyphens, `.`, `'s`, `|`), each on the protein of `pq.faa`'s first query.
- `ps_ids.faa`: survey T2's subjects with NCBI ids (ref, sp, lcl, gb, gi+ref with `<@...@>` template tags in the
  title, pdb, gnl, dbj).
- `pmix_ids.faa`: ids the T2 file lacks (pir, prf, pat, tpg, gi+gb, `gnl|TI|`, gi only), `lcl|loc_x`, a record with
  two deflines (sp, emb) and one with nine, separated by ctrl-A (a `-subject` report shows the first id and title
  only), titles with `<&>`.

Masks (survey T4's `in2/gen2.py`, seed 20261120; each design is an alignment built column by column, so each mask's
alignment columns are known, listed in the title):
- `nq_mask.fna`, `ns_mask.fna`: blastn, lowercase masks on both rows, one from column 0 past the first line, masks that
  end at a line end or start at a line start, two masks one column apart, masks next to and around gaps;
  `ns_mask_rc.fna` is the reverse complement of `ns_mask.fna` (minus strand).
- `nq_dust.fna`, `ns_dust.fna`: `(CA)x40` then 220 random bases (a DUST mask from column 0), 2 mismatches.
- `pq_mask.faa`, `ps_mask.faa`: proteins with lowercase masks on both rows (column 0 to 69 and 64, around gaps,
  masks ending at line offset 0); `nq_pmask.fna`, `ns_pmask.fna`: their back-translations with the masked codons
  in lowercase (blastx, tblastn, tblastx).
- `pq_seg.faa`, `ps_seg.faa`: 70 Q then 190 random residues (a SEG mask from column 0), 3 mismatches; `ns_seg.fna` is
  the back-translation of `ps_seg.faa`; `pq_segadj.faa` is `pq_seg.faa` with residues 71-80 in lowercase, a mask
  adjacent to SEG's (two `<font>` spans back to back).

Earlier probe inputs (survey T3's `gen.py`, seed 20261114):
- `qn_lc.fa`, `sn_lc.fa`: a nucleotide query and subject with lowercase runs and a poly-A run (DUST with
  `-soft_masking false`; tblastx's default SEG on the translated query).
