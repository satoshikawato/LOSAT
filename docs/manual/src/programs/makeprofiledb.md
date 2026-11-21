# makeprofiledb

`losat makeprofiledb` builds a profile database, the kind of database that [rpsblast](rpsblast.md) and [rpstblastn](rpstblastn.md) search, from a list of PSSM files. CDD is such a database: NCBI builds it from CDD's `.smp` files with `makeprofiledb`, and you can build it, a part of it, or a database of your own profiles the same way.

It writes the same files as NCBI's `makeprofiledb` 2.17.0 with `-dbtype rps`, byte for byte, with the same messages and exit statuses.

## Synopsis

```bash
losat makeprofiledb -in LIST [-out NAME] [options]
```

With the NCBI names installed (see [Coming from NCBI BLAST+](../from-ncbi-blast.md)): `makeprofiledb -in LIST [options]`.

## Input

- `-in` names a text file that lists the PSSM files, separated by spaces, tabs or line breaks. A name is read relative to the current directory, not to the list's directory.
- Each PSSM file is a `PssmWithParameters` in ASN.1 text: CDD's `.smp` files, or the files that `psiblast -out_pssm` writes. A file with scores (`finalData`) is taken as it is; a file with frequency ratios only gets its scores computed from them, with the scoring matrix, the gap costs and the `-scale` of the database, as NCBI does. Every file needs its frequency ratios (`intermediateData`).
- The first profile of a volume sets the volume's scoring matrix, gap costs and scale factor (from its `rpsdbparams` and `scalingFactor`, else from `-matrix`, `-gapopen`, `-gapextend` and `-scale`); every later profile must agree with them, or the run stops with NCBI's error (`Scaling factors do not match`, `Score matrix does not match`, ...).

## Common options

| Option | Default | Meaning |
| --- | --- | --- |
| `-in` | none | The list of PSSM files. |
| `-out` | the `-in` name | The name of the database. Directories in it are created. An existing database of that name is replaced. |
| `-title` | the `-in` name | The title of the database. |
| `-threshold` | `9.82` | The word score threshold of the word table, multiplied by the scale factor. A profile's own `wordScoreThreshold` is used instead when it has one, unless `-force` is given. |
| `-force` | off | Use `-threshold` for every profile. |
| `-max_smp_vol` | `2500` | The largest number of profiles in one volume. A longer list is split into volumes of nearly equal size, with an alias file (`NAME.pal`). |
| `-blastdb_version` | `5` | The format of the protein database of the profiles' sequences: `5` or `4`. |
| `-index` | `true` | With `true` the sequences' ids are parsed and indexed, so that `blastdbcmd -entry` finds them. With `false` each sequence gets an id `gnl|BL_ORD_ID|N`. |
| `-matrix` | `BLOSUM62` | The scoring matrix of profiles that do not name one. |
| `-gapopen`, `-gapextend` | the matrix's | The gap costs of profiles that do not give them. |
| `-scale` | `100` | The scale factor of profiles that have frequency ratios only, rounded up to an integer. |
| `-taxid` | none | A taxonomy id for every sequence. |
| `-taxid_map` | none | A file of `ID TAXID` lines that gives the taxonomy ids by sequence id. Without `-taxid` and `-taxid_map`, the organism of each profile's sequence gives its taxonomy id, as in NCBI. |
| `-logfile` | none | Write the messages to this file. |

`-obsr_threshold` and `-exclude_invalid` are read and have no effect with `-dbtype rps`, as in NCBI. For every option: `losat makeprofiledb -help`.

## Output

| File | What it holds |
| --- | --- |
| `NAME.rps` | The scores of every profile. |
| `NAME.freq` | The frequency ratios of every profile. |
| `NAME.loo` | The table of the words that score at least the threshold against each profile. |
| `NAME.aux` | The scoring matrix, the gap costs, the scale factor and each profile's length and K. |
| `NAME.pin`, `.phr`, `.psq`, ... | The protein BLAST database of the profiles' sequences (version 5: also `.pdb`, `.pos`, `.pot`, `.ptf`, `.pto`; version 4: `.psd`, `.psi`). |

With several volumes the files are `NAME.00.rps`, `NAME.01.rps`, ... and `NAME.pal` lists the volumes; the standard error then shows `Created protein BLAST (alias) database NAME with N sequences`. When an existing database is replaced, the standard output shows `Deleted existing BLAST database with identical name.`

After an error the RPS files are removed and the sequences already added stay in the protein database files, as NCBI leaves them.

## Differences from NCBI BLAST+

None in the files or the messages of a normal run. A profile with scores near the limits of 32-bit integers, or a word score threshold (times the scale factor) beyond them, gets the table of words of the exact values, where NCBI's sums overflow.

## What LOSAT refuses

- `-binary` (PSSM files in ASN.1 binary). Convert them to text, or wait for a later release.
- `-dbtype delta` and `-dbtype cobalt`.
- `-max_smp_vol 0`, on which NCBI's `makeprofiledb` crashes (division by zero).
- `-in` naming a directory, on which NCBI's `makeprofiledb` never ends.
- An empty `-out`, on which NCBI's `makeprofiledb` stops with an unknown exception or aborts.
- A profile stored by row (`byRow TRUE`) with more than 28 positions, whose frequency ratios NCBI's `makeprofiledb` reads from outside its data.
- A frequency ratio too large for NCBI's integers (such as `1e300`), and a word score threshold that is not a number.
- A sequence that has only a GI id between sequences with other ids, in a version 5 database with `-index true`: NCBI's `makeprofiledb` aborts there; `-blastdb_version 4` builds it.

See [What LOSAT refuses](../refusals.md).

## Examples

Build a database from three CDD profiles and search it:

```bash
ls cd00001.smp cd00002.smp cd00003.smp > cdd3.lst
losat makeprofiledb -in cdd3.lst -out cdd3 -title "three CDD domains"
losat rpsblast -query proteins.faa -db cdd3 -outfmt 6
```

Build a database from PSSMs that `psiblast` saved:

```bash
losat psiblast -query q1.faa -db proteins -num_iterations 3 -out_pssm q1.smp
losat psiblast -query q2.faa -db proteins -num_iterations 3 -out_pssm q2.smp
printf 'q1.smp\nq2.smp\n' > mine.lst
losat makeprofiledb -in mine.lst -out mine
```

## See also

- [rpsblast](rpsblast.md)
- [rpstblastn](rpstblastn.md)
- [psiblast](psiblast.md)
- [BLAST databases](../databases.md)
