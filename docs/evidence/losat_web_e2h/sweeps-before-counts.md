# Before counts: S11 build (`~/.cache/losat-web-gui-target/sf/bin/LOSAT-native`, before any port)

Date 2026-10-08. NCBI side: BLAST+ 2.17.0 frozen in `$BUILD_ROOT/sfb-e2h/sweeps/ncbi/` (two runs, identical hashes). LOSAT run with `--jobs 3`; TSVs with one line per case: `$BUILD_ROOT/sfb-e2h/sweeps/before/{fasta_sweep,check_inputs}.tsv`.

Classes: `same`, `same-error`, `explicit-rejection` (listed in AUTHORITY.md section J, or an option rejection that check_inputs expects), `explicit-rejection/unlisted` (a rejection the port must turn into `same`/`same-error`, or add to AUTHORITY.md), `approved-exception`, `differs`, `timeout`.

## fasta_sweep.py (2280 cases)

| program_task | role | outfmt | same | same-error | explicit-rejection | explicit-rejection/unlisted | approved-exception | differs | timeout | total |
|---|---|---|---|---|---|---|---|---|---|---|
| blastn | q | 0 | 36 | 0 | 0 | 78 | 0 | 0 | 0 | 114 |
| blastn | q | 6 | 36 | 0 | 0 | 78 | 0 | 0 | 0 | 114 |
| blastn | s | 0 | 35 | 0 | 0 | 77 | 2 | 0 | 0 | 114 |
| blastn | s | 6 | 39 | 0 | 0 | 75 | 0 | 0 | 0 | 114 |
| blastn -task blastn | q | 0 | 36 | 0 | 0 | 78 | 0 | 0 | 0 | 114 |
| blastn -task blastn | q | 6 | 36 | 0 | 0 | 78 | 0 | 0 | 0 | 114 |
| blastn -task blastn | s | 0 | 35 | 0 | 0 | 77 | 2 | 0 | 0 | 114 |
| blastn -task blastn | s | 6 | 39 | 0 | 0 | 75 | 0 | 0 | 0 | 114 |
| blastp | q | 0 | 41 | 0 | 0 | 73 | 0 | 0 | 0 | 114 |
| blastp | q | 6 | 41 | 0 | 0 | 73 | 0 | 0 | 0 | 114 |
| blastp | s | 0 | 40 | 0 | 0 | 74 | 0 | 0 | 0 | 114 |
| blastp | s | 6 | 43 | 0 | 0 | 71 | 0 | 0 | 0 | 114 |
| tblastn | q | 0 | 41 | 0 | 0 | 73 | 0 | 0 | 0 | 114 |
| tblastn | q | 6 | 41 | 0 | 0 | 73 | 0 | 0 | 0 | 114 |
| tblastn | s | 0 | 35 | 0 | 0 | 77 | 2 | 0 | 0 | 114 |
| tblastn | s | 6 | 39 | 0 | 0 | 75 | 0 | 0 | 0 | 114 |
| tblastx | q | 0 | 36 | 0 | 0 | 78 | 0 | 0 | 0 | 114 |
| tblastx | q | 6 | 36 | 0 | 0 | 78 | 0 | 0 | 0 | 114 |
| tblastx | s | 0 | 35 | 0 | 0 | 77 | 2 | 0 | 0 | 114 |
| tblastx | s | 6 | 39 | 0 | 0 | 75 | 0 | 0 | 0 | 114 |
| **all** |  |  | 759 | 0 | 0 | 1513 | 8 | 0 | 0 | 2280 |

Rejections that the port must remove (`explicit-rejection/unlisted`, 1513 cases), by message:

- 196 x `Error: failed to read query FASTA <file> (FASTA that bio cannot read (such as text before the first defline or bytes that are not UTF-N), which NCBI B`
- 182 x `Error: failed to read subject FASTA <file> (FASTA that bio cannot read (such as text before the first defline or bytes that are not UTF-N), which NCBI`
- 86 x `Error: query record N has a defline that has the control character NxNN; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT`
- 80 x `Error: subject record N has a defline that has the control character NxNN; NCBI BLAST+ reads such a defline differently, which is not supported by LOS`
- 80 x `Error: subject record N (id) has NxNN at residue N, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not`
- 60 x `Error: query record N (id) has NxNN at residue N, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not s`
- 60 x `Error: query record N has a defline that starts with '?' (NCBI BLAST+ reads '>?' as a gap in the sequence, and '>?_' as a defline without the prefix);`
- 60 x `Error: subject record N has a defline that starts with '?' (NCBI BLAST+ reads '>?' as a gap in the sequence, and '>?_' as a defline without the prefix`
- 36 x `Error: query record N has byte NxNN in a sequence line; NCBI BLAST+ removes it from the protein sequence with a warning, which is not supported by LOS`
- 32 x `Error: subject record N has a non-ASCII byte in a sequence line; NCBI BLAST+ reads it as an invalid residue, which is not supported by LOSAT's BLASTN `
- 30 x `Error: query record N (id) has 'N' at residue N, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not su`
- 22 x `Error: subject record N has a defline that has a non-ASCII byte; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT's BLAST`
- 20 x `Error: subject record N has byte NxNN in a sequence line; NCBI BLAST+ removes it from the protein sequence with a warning, which is not supported by L`
- 18 x `Error: query record N (id) has '|' at residue N, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not su`
- 16 x `Error: query record N (id) has no residues; NCBI BLAST+ reports such a record differently, which is not supported by LOSAT's BLASTN`
- 16 x `Error: query record N has a defline that starts with white space; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT's BLAS`
- 16 x `Error: query record N has a non-ASCII byte in a sequence line; NCBI BLAST+ reads it as an invalid residue, which is not supported by LOSAT's BLASTN (u`
- 16 x `Error: subject record N has a defline that is empty; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT's BLASTN (use ASCII`
- 16 x `Error: subject record N (id) has '|' at residue N, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not `
- 16 x `Error: subject record N (id) has ';' at residue N, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not `
- 16 x `Error: subject record N has a non-ASCII byte in a sequence line; NCBI BLAST+ reads it as an invalid residue, which is not supported by LOSAT's TBLASTX`
- 16 x `Error: query record N (id) has the residue byte NxNN at position N; NCBI BLAST+ removes such a character from a protein sequence, which is not support`
- 16 x `Error: subject record N has a non-ASCII byte in a sequence line; NCBI BLAST+ reads it as an invalid residue, which is not supported by LOSAT's TBLASTN`
- 14 x `Error: subject record N has a defline that starts with white space; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT's BL`
- 12 x `Error: query record N (id) has '#' at residue N, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not su`
- 12 x `Error: query record N (id) has '&' at residue N, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not su`
- 12 x `Error: query record N (id) has ';' at residue N, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not su`
- 12 x `Error: subject record N (id) has no residues; NCBI BLAST+ reports such a record differently, which is not supported by LOSAT's BLASTN`
- 12 x `Error: subject record N has a defline that has a non-ASCII byte; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT's TBLAS`
- 12 x `Error: subject record N (id) has the residue byte NxNN at position N; NCBI BLAST+ removes such a character from a protein sequence, which is not suppo`

## check_inputs.py (1032 cases)

| program_task | role | outfmt | same | same-error | explicit-rejection | explicit-rejection/unlisted | approved-exception | differs | timeout | total |
|---|---|---|---|---|---|---|---|---|---|---|
| blastn | - |  | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 1 |
| blastn | - | +6 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| blastn | - | 0 | 39 | 38 | 5 | 5 | 10 | 0 | 0 | 97 |
| blastn | - | 07 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| blastn | - | 6 | 57 | 17 | 13 | 24 | 4 | 0 | 0 | 115 |
| blastn | - | 7 | 11 | 4 | 0 | 2 | 0 | 0 | 0 | 17 |
| blastn | - | 99 | 0 | 2 | 0 | 0 | 0 | 0 | 0 | 2 |
| blastn | - | abc | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 3 |
| blastn -task BLASTN | - | 0 | 0 | 0 | 0 | 0 | 1 | 0 | 0 | 1 |
| blastn -task blastn | - | 0 | 11 | 2 | 0 | 0 | 0 | 0 | 0 | 13 |
| blastn -task blastn | - | 6 | 22 | 1 | 1 | 0 | 0 | 0 | 0 | 24 |
| blastn -task blastn | - | 7 | 5 | 0 | 1 | 0 | 0 | 0 | 0 | 6 |
| blastn -task blastn-short | - | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| blastn -task dc-megablast | - | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| blastn -task megablast | - | 0 | 8 | 0 | 0 | 0 | 0 | 0 | 0 | 8 |
| blastn -task megablast | - | 6 | 8 | 0 | 0 | 0 | 0 | 0 | 0 | 8 |
| blastn -task rmblastn | - | 0 | 0 | 0 | 1 | 0 | 0 | 0 | 0 | 1 |
| blastp | - | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| blastp | - | 7 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| blastp | q | 0 | 23 | 1 | 0 | 38 | 0 | 0 | 0 | 62 |
| blastp | q | 6 | 23 | 1 | 0 | 38 | 0 | 0 | 0 | 62 |
| blastp | s | 0 | 20 | 5 | 0 | 34 | 0 | 0 | 0 | 59 |
| blastp | s | 6 | 23 | 5 | 1 | 30 | 0 | 0 | 0 | 59 |
| tblastn | - | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| tblastn | - | 7 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| tblastn | q | 0 | 23 | 1 | 0 | 38 | 0 | 0 | 0 | 62 |
| tblastn | q | 6 | 23 | 1 | 0 | 38 | 0 | 0 | 0 | 62 |
| tblastn | s | 0 | 11 | 6 | 0 | 40 | 2 | 0 | 0 | 59 |
| tblastn | s | 6 | 14 | 6 | 1 | 38 | 0 | 0 | 0 | 59 |
| tblastx | - | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| tblastx | - | 7 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
| tblastx | q | 0 | 20 | 1 | 0 | 41 | 0 | 0 | 0 | 62 |
| tblastx | q | 6 | 20 | 1 | 0 | 41 | 0 | 0 | 0 | 62 |
| tblastx | s | 0 | 11 | 6 | 0 | 40 | 2 | 0 | 0 | 59 |
| tblastx | s | 6 | 14 | 6 | 1 | 38 | 0 | 0 | 0 | 59 |
| **all** |  |  | 396 | 108 | 24 | 485 | 19 | 0 | 0 | 1032 |

Rejections that the port must remove (`explicit-rejection/unlisted`, 485 cases), by message:

- 34 x `Error: query record N has a defline that has the control character NxNN; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT`
- 33 x `Error: subject record N has a defline that has the control character NxNN; NCBI BLAST+ reads such a defline differently, which is not supported by LOS`
- 26 x `Error: query record N has a defline that starts with '?' (NCBI BLAST+ reads '>?' as a gap in the sequence, and '>?_' as a defline without the prefix);`
- 25 x `Error: subject record N has a defline that starts with '?' (NCBI BLAST+ reads '>?' as a gap in the sequence, and '>?_' as a defline without the prefix`
- 16 x `Error: subject record N (id) has NxNN at residue N, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not`
- 12 x `Error: query record N has byte NxNN in a sequence line; NCBI BLAST+ removes it from the protein sequence with a warning, which is not supported by LOS`
- 8 x `Error: subject record N has a defline that has a non-ASCII byte; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT's TBLAS`
- 8 x `Error: query record N (id) has no residues; NCBI BLAST+ reports such a record differently, which is not supported by LOSAT's TBLASTX`
- 8 x `Error: subject record N (id) has no residues; NCBI BLAST+ reports such a record differently, which is not supported by LOSAT's TBLASTX`
- 8 x `Error: query record N (id) has the residue byte NxNN at position N; NCBI BLAST+ removes such a character from a protein sequence, which is not support`
- 8 x `Error: query record N (id) has no residues; NCBI BLAST+ reports such a record differently, which is not supported by LOSAT's TBLASTN`
- 8 x `Error: subject record N (id) has no residues; NCBI BLAST+ reports such a record differently, which is not supported by LOSAT's TBLASTN`
- 8 x `Error: query record N (id) has no residues; NCBI BLAST+ reports such a record differently, which is not supported by LOSAT's BLASTP`
- 6 x `Error: query record N (id) has NxNN at residue N, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not s`
- 6 x `Error: query record N has a defline that has a non-ASCII byte; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT's TBLASTX`
- 6 x `Error: failed to read query FASTA <path> (FASTA that bio cannot read (such as text before the first defline or byt`
- 6 x `Error: failed to read subject FASTA <path> (FASTA that bio cannot read (such as text before the first defline or b`
- 6 x `Error: failed to read query FASTA <path> (FASTA that bio cannot read (such as text before the first defline or bytes that `
- 6 x `Error: subject record N has byte NxNN in a sequence line; NCBI BLAST+ removes it from the protein sequence with a warning, which is not supported by L`
- 5 x `Error: subject record N has a defline that has a non-ASCII byte; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT's BLAST`
- 5 x `Error: an empty query from a stream without a position (such as a pipe) is not supported by LOSAT's BLASTN`
- 5 x `Error: failed to read subject FASTA <path> (FASTA that bio cannot read (such as text before the first defline o`
- 4 x `Error: subject record N has an HTML character reference (such as &amp;) in its defline, which NCBI BLAST+ decodes in the outfmt N titles; this is not `
- 4 x `Error: query record N has a defline that starts with white space; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT's TBLA`
- 4 x `Error: subject record N has a defline that starts with white space; NCBI BLAST+ reads such a defline differently, which is not supported by LOSAT's TB`
- 4 x `Error: subject record N (id) has 'X' at residue N, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not `
- 4 x `Error: subject record N (id) has '-' at residue N, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not `
- 4 x `Error: subject record N (id) has 'N' at residue N, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not `
- 4 x `Error: subject record N (id) has '*' at residue N, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not `
- 4 x `Error: subject record N (id) has ';' at residue N, which is not an IUPAC nucleotide letter; NCBI BLAST+ reads such a record differently, which is not `
