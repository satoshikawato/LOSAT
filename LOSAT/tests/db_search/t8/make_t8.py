#!/usr/bin/env python3
"""Line RG item T8a: the inputs of the rows that pin where NCBI's texts land in a report written to
`-out /dev/stdout` with standard error in the same pipe (`tests/taxdb_tool_texts.rs`, the rows of
`db_search_fixtures.py` device `out-stdout-merged`).

Writes, from LOSAT/ (`python3 tests/db_search/t8/make_t8.py NCBI_BIN`, NCBI_BIN the directory of NCBI
BLAST+ 2.17.0's makeblastdb), in `tests/db_search/t8/`:

- `nl.fna`, four nucleotide records of 23,038, 900, 700 and 9,000 letters (ids `gb|LN000001.1|` to
  `gb|LN000004.1|`), longer than NCBI's 8,191-byte file buffer, whose letters come from a fixed linear
  congruential generator below (no data of anyone's); `nl.map`, their taxids (9606 and 99990 in turn:
  99990 is the last entry of `tests/seqdb/taxdb/<variant>`, whose lookup writes the offset error in
  `offset_last_past_end`); and the version 5 database `nl` that NCBI's makeblastdb writes of them with
  `-parse_seqids` and `-taxid_map nl.map`. blastdbcmd and dustmasker read it with a damaged taxonomy
  database, so that the texts of the lookups come between NCBI's writes of a record.
- `nu.fna`, `nu.map` and the database `nu` likewise, three records of 3,900 (9606), 23,038 (9606) and
  900 (99990) letters (ids `gb|LU000001.1|` to `gb|LU000003.1|`): the first record's bytes are held in
  NCBI's buffer when the second one's letters come, so that NCBI's write within the second record ends
  where a piece of letters of NCBI's sequence reader (1,024 letters from the start of the output) does
  not fit the buffer, which is not a multiple of 8,192 letters; the offset error of the third record's
  lookup lands after it.
- `wq.faa` and `wq.fna`, six queries `w1` to `w6`, each the letters of `tests/db_search/ids/txr_q.faa`
  (`txr_q.fna`), where `w3` has a hyphen after its 70th letter and `w5` a digit `1` after its 125th: the
  FASTA reader's warnings of the third and fifth query (with `BATCH_SIZE=1`, of a later query batch).
- `wh.fna`, three queries of `ids/txr_q.fna`'s letters: `h1`, `h2` with a hyphen after its 70th letter,
  and `h3 café` (UTF-8), whose title NCBI's search strategy (`-export_search_strategy`) writes with its
  `Bad char` posts at that query's batch (blastx's held run, `tests/tax_report_streams.rs`).
- `bt.txt`, blastdbcmd's `-entry_batch` of `../t5b/nmv` with two ids it does not find (`XX000001.1`,
  `YY000002.1`) after two it finds.
- `nq.fna`, two nucleotide subjects `gb|NQ00000.1|` and `gb|NQ00001.1|`, each the letters of
  `ids/txr_q.fna`, and its version 5 database `nq` (makeblastdb `-parse_seqids`): the subjects of the
  nucleotide databases' searches of `wq.*`.

Python 3 standard library only; makeblastdb is NCBI's (the comparison oracle builds the inputs).
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
IDS = os.path.join(HERE, "..", "ids")
RECORDS = [("gb|LN000001.1| long one", 23038, 9606), ("gb|LN000002.1| two", 900, 99990),
           ("gb|LN000003.1| three", 700, 9606), ("gb|LN000004.1| four", 9000, 99990)]
UNIT_RECORDS = [("gb|LU000001.1| first", 3900, 9606), ("gb|LU000002.1| long", 23038, 9606),
                ("gb|LU000003.1| two", 900, 99990)]
# The queries of `wq.*` that get a character the reader warns of: query number, after how many letters, character.
INSERTS = {3: (70, "-"), 5: (125, "1")}


def letters(count, state):
    """`count` letters of ACGT from the generator's `state`; the letters and the next state."""
    out = []
    for _ in range(count):
        state = (state * 6364136223846793005 + 1442695040888963407) % (1 << 64)
        out.append("ACGT"[state >> 62])
    return "".join(out), state


def one_sequence(path):
    """The letters of the FASTA file `path`, which holds one sequence."""
    with open(path) as f:
        return "".join(line.strip() for line in f if not line.startswith(">"))


def fasta_lines(seq):
    return [seq[i:i + 60] for i in range(0, len(seq), 60)]


def write(name, lines):
    with open(os.path.join(HERE, name), "w", newline="\n", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def makeblastdb(ncbi_bin, args):
    with tempfile.TemporaryDirectory() as home:
        subprocess.run([os.path.join(ncbi_bin, "makeblastdb"), *args, "-blastdb_version", "5"],
                       cwd=HERE, check=True, stdout=subprocess.DEVNULL,
                       env={"PATH": "/usr/bin:/bin", "BLAST_USAGE_REPORT": "false", "HOME": home})


def main(argv):
    if len(argv) != 1:
        sys.exit("usage: make_t8.py NCBI_BIN")
    state = 20261125
    for name, records in (("nl", RECORDS), ("nu", UNIT_RECORDS)):
        fasta, taxids = [], []
        for title, count, taxid in records:
            seq, state = letters(count, state)
            fasta.append(">" + title)
            fasta.extend(fasta_lines(seq))
            taxids.append(f"{title.split('|')[1]} {taxid}")
        write(f"{name}.fna", fasta)
        write(f"{name}.map", taxids)
        makeblastdb(argv[0], ["-in", f"{name}.fna", "-dbtype", "nucl", "-parse_seqids", "-taxid_map", f"{name}.map",
                              "-out", name])
    for suffix in ("faa", "fna"):
        seq = one_sequence(os.path.join(IDS, f"txr_q.{suffix}"))
        lines = []
        for number in range(1, 7):
            rows = fasta_lines(seq)
            if number in INSERTS:
                at, char = INSERTS[number]
                row, column = divmod(at, 60)
                rows[row] = rows[row][:column] + char + rows[row][column:]
            lines += [f">w{number} query {number}", *rows]
        write(f"wq.{suffix}", lines)
    seq = one_sequence(os.path.join(IDS, "txr_q.fna"))
    hyphen = fasta_lines(seq)
    hyphen[1] = hyphen[1][:10] + "-" + hyphen[1][10:]
    write("wh.fna", [">h1 query 1", *fasta_lines(seq), ">h2 query 2", *hyphen, ">h3 caf\u00e9", *fasta_lines(seq)])
    write("nq.fna", [line for number in range(2)
                     for line in (f">gb|NQ0000{number}.1| nucleotide subject {number}", *fasta_lines(seq))])
    makeblastdb(argv[0], ["-in", "nq.fna", "-dbtype", "nucl", "-parse_seqids", "-out", "nq"])
    write("bt.txt", ["AB410002.1", "XX000001.1", "AB410003.1", "YY000002.1"])


if __name__ == "__main__":
    main(sys.argv[1:])
