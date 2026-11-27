#!/usr/bin/env python3
"""Real native/WASI threading gates with raw NCBI 2.17.0 oracle output."""
import argparse
import hashlib
import json
import os
import platform
from pathlib import Path
import re
import shlex
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
TESTS = ROOT / "LOSAT/tests"
sys.path.insert(0, str(TESTS))
import db_search_fixtures
import losat_wording  # noqa: E402  (LOSAT's words in NCBI's bytes, AGENTS.md rule 8)
from wasm_performance import validate_thread_evidence

# The output formats in which blastx reads its query file before the search (BLAST XML 5, and XML2 and JSON 13 to
# 16; LOSAT/src/algorithm/blastx/blast_xml.rs `prescan_queries`), whose frozen rows run through each WASI command
# (`blastx_prescan_checkers`).
BLASTX_PRESCAN_FORMATS = "5,13,14,15,16"


def blastx_prescan_checkers(wrapper):
    """The checkers of the native gates that run blastx's rows of BLASTX_PRESCAN_FORMATS through `wrapper` (a
    LOSAT command), at -num_threads 1, exactly: (name, argv after the interpreter, the marker of a differing line)
    for the -subject rows (docs/evidence/losat_web_e2a/check_losat.py) and the -db rows (db_search_fixtures.py)."""
    return [("subject", [ROOT / "docs/evidence/losat_web_e2a/check_losat.py", "--losat", wrapper, "--programs", "blastx",
                         "--threads", "1", "--outfmts", BLASTX_PRESCAN_FORMATS], "\tDIFF\t"),
            ("db", [TESTS / "db_search_fixtures.py", "check", "--losat", wrapper, "--program", "blastx",
                    "--threads", "1", "--outfmts", BLASTX_PRESCAN_FORMATS], "\tdiffers\t")]


def digest(data):
    return hashlib.sha256(data).hexdigest()


# NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-180
# Distinct IDs and unequal lengths exercise global statistics and stable ordering.
def fixtures(directory):
    directory.mkdir(parents=True, exist_ok=True)
    seq = "".join(line.strip() for line in (TESTS / "fasta/small_test.fasta").read_text().splitlines() if not line.startswith(">"))
    aa = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEQFKYVQKDVNAYLKDAQVLGFLYEVHDDPGLQRLFFKGEKPKYEE"
    records = {
        "nuc1": [seq[:900]], "nuc3": [seq[:900]] * 3,
        "aa1": [aa], "aa3": [aa] * 3,
        "unequal": [seq[:899], seq[:900], seq[:903]],
        "remainder": [seq[:900], seq[:901], seq[:902]],
        "short-context": [seq[:900], seq[:450], seq[:10]],
        "short-subject": [seq[:900], seq[:10], seq[:2]],
    }
    for name, values in records.items():
        (directory / f"{name}.fasta").write_text("".join(f">seq{i}\n{value}\n" for i, value in enumerate(values)))
    # NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:415-439,1264-1283
    # Empty BLAST defline object prints N/A; every query receives its own header.
    edge_title = "a long title, with spaces and hyphens " * 5
    records["aa-query-edge"] = [aa, "W" * 85, aa]
    records["aa-subject-edge"] = [aa, aa]
    for name in ["aa-query-edge", "aa-subject-edge"]:
        (directory / f"{name}.fasta").write_text("".join(
            f">seq{i}" + (f" {edge_title}" if i == 0 else "") + f"\n{value}\n"
            for i, value in enumerate(records[name])))
    # NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:441-446
    # Fixed standard-code codons encode the existing protein fixture exactly.
    codons = dict(zip("ACDEFGHIKLMNPQRSTVWY", ["GCT","TGT","GAT","GAA","TTT","GGT","CAT","ATT","AAA","CTG","ATG","AAT","CCT","CAA","CGT","TCT","ACT","GTT","TGG","TAT"]))
    records["blastx-query"] = ["".join(codons[x] for x in aa)]
    (directory / "blastx-query.fasta").write_text(">seq0\n" + records["blastx-query"][0] + "\n")
    # TBLASTN subjects: the coding sequence, a nucleotide prefix and a short piece of the coding sequence.
    records["tblastn-subjects"] = [records["blastx-query"][0], seq[:900], records["blastx-query"][0][:120]]
    (directory / "tblastn-subjects.fasta").write_text("".join(f">seq{i}\n{value}\n" for i, value in enumerate(records["tblastn-subjects"])))
    return {name: directory / f"{name}.fasta" for name in records}


# The -db rows whose alias files declare a total of 2^32 letters or more (LOSAT/tests/db_search/dbcounts:
# LENGTH 4294967296000000000 in l4295, 4294967295999999999 in l4294, 68719476736 and 98784247808 in
# avg2p32; LOSAT/tests/db_search/aliaskeys: 4611686018427387904 in n_nseq0_l62, 4611686018427387907 in
# n_nseq1_l62p3 and 36028797018963968 in p_l55), which the WASI commands run (`main`, "db-total").
DB_TOTAL_ROWS = ["db.cnt.l4295.bn", "db.cnt.avg2p32.bn", "db.cnt.l4295.bp.0", "db.cnt.l4295.bp.5",
                 "db.cnt.l4295.bp.15", "db.cnt.l4294.bp.0", "db.cnt.avg2p32.bp_c0", "db.cnt.avg2p32.tx",
                 "db.cnt.avg2p32.mb.split.6", "db.cnt.avg2p32.mb.split.0", "db.cnt.avg2p32.mb.split.7",
                 "db.cnt.avg2p32.mb.split.b1.7", "db.cnt.avg2p32.mb.split.nq", "db.cnt.avg2p32.mb.whole.6",
                 "db.cnt.avg2p32.mo00.split.6", "db.cnt.avg2p32.mo52.split.6",
                 "db.cnt.avg2p32.g3.split.6", "db.cnt.avg2p32.g3.split.0", "db.cnt.avg2p32.g3.split.b1.7",
                 "db.cnt.avg2p32.g3.split.mo00.6", "db.cnt.avg2p32.g3.split.mb40k.6",
                 "db.ak.bn.nseq0.l62.alln.7", "db.ak.tbn.nseq0.l62.allx.7",
                 "db.ak.tbn.nseq1.l62p3.allx.7", "db.ak.tbx.nseq1.l62p3.alln.7",
                 "db.xc.dup.split.avg.7", "db.xc.gapvalid.avg.7", "db.xc.ord.neggap.avg.7",
                 "db.ak.bp.l55.hba.0", "db.ak.bp.l55.hba.5", "db.ak.bp.l55.hba.6", "db.ak.bp.l55.hba.7",
                 "db.ak.bp.l55.x1hba.7", "db.ak.bp.l55.x6.7", "db.ak.ps.l55.hba.7"]


# NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1098-1108
# Retain exact argv, stderr and outputs. No sorting or output normalization.
# NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1104-1105
# NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:88-97
# The local times that makeblastdb's log ("M/D/Y h:m:s", to the second) and the database's date (to the minute:
# blastdbcmd -info's "Date:", a search report's "Posted date:") print, and the time that makeblastdb's log gives
# for the build (`local_times`).
LOG_TIME = re.compile(rb"(?<=^Building a new DB, current time: )(\d\d)/(\d\d)/(\d{4}) (\d\d):(\d\d):(\d\d)$", re.M)
DB_DATE = re.compile(rb"(?<=^Date: )([A-Z][a-z]{2}) (\d{1,2}), (\d{4})  (\d{1,2}):(\d\d) ([AP]M)(?=\t)"
                     rb"|(?<=^    Posted date:  )([A-Z][a-z]{2}) (\d{1,2}), (\d{4})  (\d{1,2}):(\d\d) ([AP]M)$", re.M)
BUILD_SECONDS = re.compile(rb"(?<= sequences in )[0-9.e+-]+(?= seconds\.$)", re.M)
MONTHS = [b"Jan", b"Feb", b"Mar", b"Apr", b"May", b"Jun", b"Jul", b"Aug", b"Sep", b"Oct", b"Nov", b"Dec"]


def local_times(data):
    """`data` with each local time of LOG_TIME and DB_DATE replaced by `<time>` and the build's seconds by
    `<seconds>`, and the times it printed: (seconds since 1970-01-01 of the printed calendar time, read as UTC;
    the resolution in seconds: 1 for the log, 60 for the date)."""
    import calendar
    found = []

    def log_time(match):
        month, day, year, hour, minute, second = (int(x) for x in match.groups())
        found.append((calendar.timegm((year, month, day, hour, minute, second)), 1))
        return b"<time>"

    def db_date(match):
        groups = [x for x in match.groups() if x is not None]
        month, day, year, hour12, minute, half = groups
        hour = int(hour12) % 12 + (12 if half == b"PM" else 0)
        found.append((calendar.timegm((int(year), MONTHS.index(month) + 1, int(day), hour, int(minute), 0)), 60))
        return b"<time>"

    data = BUILD_SECONDS.sub(b"<seconds>", DB_DATE.sub(db_date, LOG_TIME.sub(log_time, data)))
    return data, found


def is_local_time_of(printed, resolution, offset, started, ended):
    """Whether a local time that `local_times` read (`printed`, to `resolution` seconds) is a time between the
    instants `started` and `ended` (seconds since 1970, UTC; one second of slack) in the zone `offset` seconds
    east of UTC."""
    first = int(started - 1 + offset) // resolution * resolution
    last = int(ended + 1 + offset) // resolution * resolution
    return first <= printed <= last


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--native-serial", type=Path)
    # NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_args.hpp:1290-1296
    parser.add_argument("--serial", type=Path, help="opt in to serial compatibility checks")
    parser.add_argument("--threaded", type=Path, required=True)
    parser.add_argument("--reactor", type=Path, required=True)
    parser.add_argument("--serial-reactor", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--node", default="node")
    parser.add_argument("--oracle-dir", type=Path)
    args = parser.parse_args()
    out = args.output_dir.resolve(); out.mkdir(parents=True, exist_ok=True)
    inputs = fixtures(out / "fixtures")
    rows = []
    metadata = {
        "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "node": subprocess.check_output([args.node, "-p", "JSON.stringify(process.versions)"], text=True).strip(),
        "rust": subprocess.check_output(["rustc", "-vV"], text=True), "platform": platform.platform(),
        "artifacts": {name: {"path": str(value.resolve()), "sha256": digest(value.read_bytes())}
                      for name, value in vars(args).items() if name in ["native", "native_serial", "serial", "threaded", "reactor", "serial_reactor"] and value},
        "fixtures": {name: digest(value.read_bytes()) for name, value in inputs.items()},
        "runners": {p.name: digest(p.read_bytes()) for p in TESTS.glob("*.js")}, "oracles": {},
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    environment = {k: v for k, v in os.environ.items() if not k.startswith(("LOSAT_", "RAYON_")) and k != "BL2SEQ_LEGACY"}
    environment.update(NODE_NO_WARNINGS="1", LC_ALL="C")
    # stdin: bytes through a pipe, or a Path opened as the standard input (a regular file); expected=None: any exit.
    # cwd: the run's working directory (ROOT by default); a variable of env whose value is None is left out.
    def execute(label, argv, env=None, expected=0, timeout=60, stdin=None, cwd=ROOT):
        directory = out / label; directory.mkdir(parents=True, exist_ok=True)
        result_path = directory / "result.out"
        argv = [str(x).replace("{out}", str(result_path)) for x in argv]
        run_env = {k: v for k, v in {**environment, **(env or {})}.items() if v is not None}
        try:
            if isinstance(stdin, Path):
                with open(stdin, "rb") as handle:
                    result = subprocess.run(argv, cwd=cwd, env=run_env, stdin=handle, capture_output=True, timeout=timeout)
            else:
                result = subprocess.run(argv, cwd=cwd, env=run_env, input=stdin, capture_output=True, timeout=timeout)
        except subprocess.TimeoutExpired as error:
            (directory / "stdout").write_bytes(error.stdout or b"")
            (directory / "stderr").write_bytes(error.stderr or b"")
            rows.append(dict(label=label, argv=argv, cwd=str(cwd), env=env or {}, status="TIMEOUT", timeout_seconds=timeout))
            (out / "runs.json").write_text(json.dumps(rows, indent=2) + "\n")
            raise

        (directory / "stdout").write_bytes(result.stdout); (directory / "stderr").write_bytes(result.stderr)
        data = result_path.read_bytes() if result_path.exists() else result.stdout
        row = dict(label=label, argv=argv, cwd=str(cwd), env=env or {}, exit=result.returncode, sha256=digest(data))
        rows.append(row); (out / "runs.json").write_text(json.dumps(rows, indent=2) + "\n")
        assert expected is None or (result.returncode == 0) == (expected == 0), (label, result.returncode, result.stderr.decode())
        return data, result.stderr.decode()
    prefixes = {
        "native": [args.native.resolve()],
        "threaded": [args.node, TESTS / "run_losat_wasi_threads.js", args.threaded.resolve()],
    }
    # NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_args.hpp:1290-1296
    if args.serial: prefixes["serial"] = [args.node, TESTS / "run_losat_wasi.js", args.serial.resolve()]
    if args.native_serial: prefixes["native-serial"] = [args.native_serial.resolve()]
    cases = [("blastn", "blastn", "nuc1", "nuc1", 10), ("megablast", "blastn", "nuc1", "nuc3", 10),
             ("blastp", "blastp", "aa3", "aa3", 10)]
    # BLASTX joins the parity and thread contract (its pool line names the program BLASTX).
    cases += [("blastx", "blastx", "blastx-query", "aa3", 10)]
    # TBLASTN joins it too; its -subject records are packed while they are read (LOSAT_X_TBNSTORE),
    # which must start no thread beyond the pool.
    cases += [("tblastn", "tblastn", "aa1", "blastx-query", 10),
              ("tblastn-aa3-subjects", "tblastn", "aa3", "tblastn-subjects", 10)]
    cases += [(f"tblastx-{q}-{s}-e{e}", "tblastx", q, s, e)
              for q, s in [("nuc1", "nuc1"), ("nuc1", "nuc3"), ("nuc3", "nuc1"), ("nuc3", "nuc3"),
                           ("nuc1", "unequal"), ("nuc1", "remainder"), ("short-context", "short-subject")]
              for e in ([0.1, 10, 10000] if s in ["unequal", "remainder"] else [10])]
    for label, program, query, subject, evalue in cases:
        common = ["-query", inputs[query], "-subject", inputs[subject], "-outfmt", "6", "-evalue", str(evalue), "-out", "{out}"]
        if program == "blastn": common += ["-task", label]
        oracle = str(args.oracle_dir.resolve() / program) if args.oracle_dir else program
        version = subprocess.check_output([oracle, "-version"], text=True)
        assert "2.17.0" in version, version
        import shutil
        oracle_path = Path(shutil.which(oracle)).resolve()
        metadata["oracles"][program] = dict(path=str(oracle_path), sha256=digest(oracle_path.read_bytes()), version=version)
        (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
        reference, _ = execute(label + "/oracle", [oracle, *common, "-num_threads", "1"])
        for kind, prefix in prefixes.items():
            for n in ([1] if kind in ["serial", "native-serial"] else [1, 2, 4, 8, 16]):
                data, log = execute(f"{label}/{kind}-n{n}", [*prefix, program, *common, "-num_threads", str(n)], {"LOSAT_WASI_THREADS_DEBUG": "1"})
                assert data == reference, f"raw output differs: {label}/{kind}-n{n}"
                validate_thread_evidence(log, n, "native" if kind == "native-serial" else kind)
        print(label, "raw parity + thread contract PASS", flush=True)
    # NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbatlas.cpp:114-127
    # A database file is opened once for the process and read by every thread. -db at every
    # thread count, against the frozen NCBI bytes of one outfmt 6 row per program
    # (LOSAT/tests/db_search_manifest/<program>.tsv; the Node host's threads read a descriptor that
    # another thread opened).
    manifest = {row["row_id"]: row for row in db_search_fixtures.read_manifest()[1]}
    for row_id in ["db.blastn.mega.6", "db.d1.6", "db.blastx.default.6", "db.tblastn.default.6", "db.tblastx.default.6"]:
        row = manifest[row_id]
        assert row["kind"] == "ncbi" and not row["extra_args"], row_id
        reference = (TESTS / "fixtures/db_search" / f"{row_id}.out").read_bytes()
        common = ["-db", ROOT / "LOSAT" / row["db"], "-query", ROOT / "LOSAT" / row["query"],
                  "-outfmt", row["outfmt"], "-out", "{out}"]
        for kind, prefix in prefixes.items():
            for n in ([1] if kind in ["serial", "native-serial"] else [1, 2, 4, 8, 16]):
                data, log = execute(f"db/{row_id}/{kind}-n{n}", [*prefix, row["program"], *common, "-num_threads", str(n)], {"LOSAT_WASI_THREADS_DEBUG": "1"})
                assert data == reference, f"raw output differs: {row_id}/{kind}-n{n}"
                validate_thread_evidence(log, n, "native" if kind == "native-serial" else kind)
        print(row_id, "-db raw parity + thread contract PASS", flush=True)
    # NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:904-911
    # A database total of 2^32 letters or more through each WASI command, at -num_threads 1, by the checker
    # of the native gates (db_search_fixtures.py check: stdout, stderr and the exit status against the frozen
    # NCBI bytes, exactly), which runs a wrapper that starts the command's runner in the place of the
    # binary, as the psiblast -db rows below: each row of DB_TOTAL_ROWS by its id (`--only`, one row each).
    # wasm32's usize is 32 bits; the engines hold NCBI's Int8 total (before, these rows stopped with "a
    # database of N letters is not supported by LOSAT's PROGRAM on this platform").
    for kind in (kind for kind in ["serial", "threaded"] if kind in prefixes):
        wrapper = out / "db-total" / kind / "losat"
        wrapper.parent.mkdir(parents=True, exist_ok=True)
        wrapper.write_text("#!/bin/sh\nexec " + " ".join(shlex.quote(str(x)) for x in prefixes[kind]) + ' "$@"\n')
        wrapper.chmod(0o755)
        started = time.monotonic()
        for row_id in DB_TOTAL_ROWS:
            label = f"db-total/{kind}/{row_id}"
            try:
                execute(label, [sys.executable, TESTS / "db_search_fixtures.py", "check", "--losat", wrapper,
                                "--only", row_id, "--threads", "1"], timeout=600)
            except AssertionError:
                print((out / label / "stdout").read_text(errors="replace"), flush=True)
                raise
            lines = (out / label / "stdout").read_text().splitlines()
            assert lines[0].startswith(f"{row_id}\tsame\t") and lines[-1].startswith("# rows=1 differing=0 "), (label, lines)
        print(f"db-total/{kind}", f"{len(DB_TOTAL_ROWS)} rows same", f"({time.monotonic() - started:.0f} s)", flush=True)
    # NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:448-453
    # psiblast through each WASI command, at -num_threads 1, by the checkers of the native gates, which
    # run a wrapper that starts the command's runner in the place of the binary: every psiblast row of
    # the frozen NCBI -subject fixtures (outfmt0_manifest.tsv, 95 rows; docs/evidence/losat_web_e2a/
    # check_losat.py: stdout, stderr, the PSSM files) and of the -db fixtures (db_search_manifest/
    # psiblast*.tsv, every row; db_search_fixtures.py check). Its PSSM queries (-in_pssm) open no query
    # input, as NCBI reads none, and its ASN.1 PSSM files are staged under a name without a process
    # id, which WASI does not have (LOSAT/src/utils/process_id.rs). The -db rows are compared byte for
    # byte, all of them: since the port of glibc's FMA exp/log/pow (LOSAT/src/utils/glibc_fma_libm.rs)
    # WASI writes the frozen bytes in the 14 rows that PD-LOSAT-PLATFORM-LIBM 1.1 named, where the
    # ASN.1 PSSM's informationContent had differed in the last of its 15 digits by 1. The device rows
    # (`device: MODE`, db_search_fixtures.py `run_on_device`) run here too, psiblast's and blastx's of the
    # formats below: where NCBI aborts (exit -6), the module's abort ends the runner by SIGABRT
    # (wasi_artifact.js `endOnAbortTrap`), as the native command's; before, the runners exited 1 after
    # the trap's stack, and the rows db.device.ps.p1epi.7 and db.device.ps.p1both.0 differed.
    # NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-873
    # blastx's BLAST XML (outfmt 5), XML2 and JSON (13 to 16) through each WASI command, at -num_threads 1: every
    # blastx row of these formats of the frozen NCBI -subject fixtures (outfmt0_manifest.tsv, 88 rows) and -db
    # fixtures (db_search_manifest/blastx*.tsv, 230 rows), byte for byte, by the checkers of the native gates
    # (`blastx_prescan_checkers`), through the wrapper of psiblast's checks (above), in the same loop.
    # In these formats LOSAT reads the query file once before the search, for its refusals before the first byte
    # (LOSAT/src/algorithm/blastx/blast_xml.rs `prescan_queries`): it reads the open file itself and puts it back
    # where it was, as NCBI's IsIStreamEmpty reads ahead on its one stream. Before, it read a second handle on the
    # file, which WASI cannot make, and every such run stopped with 'Argument "query". File is not accessible'.
    for kind in (kind for kind in ["serial", "threaded"] if kind in prefixes):
        wrapper = out / "psiblast-fixtures" / kind / "losat"
        wrapper.parent.mkdir(parents=True, exist_ok=True)
        wrapper.write_text("#!/bin/sh\nexec " + " ".join(shlex.quote(str(x)) for x in prefixes[kind]) + ' "$@"\n')
        wrapper.chmod(0o755)
        for section, name, checker, differs in [
                ("psiblast-fixtures", "outfmt0", [ROOT / "docs/evidence/losat_web_e2a/check_losat.py", "--losat", wrapper,
                                                  "--programs", "psiblast", "--threads", "1"], "\tDIFF\t"),
                ("psiblast-fixtures", "db", [TESTS / "db_search_fixtures.py", "check", "--losat", wrapper, "--program",
                                             "psiblast", "--threads", "1"],
                 "\tdiffers\t"),
                *(("blastx-xml-fixtures", *checker) for checker in blastx_prescan_checkers(wrapper))]:
            label = f"{section}/{kind}/{name}"
            started = time.monotonic()
            try:
                execute(label, [sys.executable, *checker], timeout=3600)
            except AssertionError:
                lines = (out / label / "stdout").read_text(errors="replace").splitlines()
                print("\n".join(line for line in lines if differs in line or line.startswith("# ")), flush=True)
                raise
            print(label, (out / label / "stdout").read_text().splitlines()[-1],
                  f"({time.monotonic() - started:.0f} s)", flush=True)
    # NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:715-721
    # Standard input through each WASI command at -num_threads 1, for every search program: -query -,
    # no -query (its default is -), an empty pipe, -subject - and both (the subjects read standard input
    # to its end, the query then reads nothing), and for blastn and blastx (the two FASTA streams) a
    # regular file as the standard input, which has a position. Each run must give the native command's
    # stdout, stderr and exit status of the same run, byte for byte. LOSAT reads the standard input in
    # place, never a copy of its descriptor, which WASI cannot make (LOSAT/src/blastinput/input_files.rs
    # `InputFile`); before, every such run stopped with "reading the query from standard input ('-') on
    # this platform is not supported".
    rps_db = db_search_fixtures.inflate_rps("rps-gz:rps_mini") / "rps_mini"
    stdin_programs = [
        # program, query, subjects (the file, or the RPS database), whether -subject - is read, regular-file case
        ("blastn", inputs["nuc1"], ["-subject", inputs["nuc3"]], True, True),
        ("blastp", inputs["aa3"], ["-subject", inputs["aa3"]], True, False),
        ("blastx", inputs["blastx-query"], ["-subject", inputs["aa3"]], True, True),
        ("tblastn", inputs["aa1"], ["-subject", inputs["tblastn-subjects"]], True, False),
        ("tblastx", inputs["nuc1"], ["-subject", inputs["nuc3"]], True, False),
        ("psiblast", inputs["aa1"], ["-subject", inputs["aa3"]], True, False),
        ("rpsblast", inputs["aa1"], ["-db", rps_db], False, False),
        # (LOSAT's rpstblastn refuses the default -comp_based_stats 1 before it reads a query.)
        ("rpstblastn", inputs["blastx-query"], ["-db", rps_db, "-comp_based_stats", "0"], False, False),
    ]
    for program, query, subjects, subject_stdin, regular in stdin_programs:
        cases = [("query", ["-query", "-", *subjects], query.read_bytes()),
                 ("omitted", subjects, query.read_bytes()),
                 ("empty", ["-query", "-", *subjects], b"")]
        if regular:
            cases.append(("query-file", ["-query", "-", *subjects], query))
        if subject_stdin:
            cases += [("subject", ["-query", query, "-subject", "-"], Path(subjects[1]).read_bytes()),
                      ("both", ["-query", "-", "-subject", "-"], Path(subjects[1]).read_bytes())]
        for case, argv, stdin in cases:
            common = [program, *argv, "-outfmt", "7", "-num_threads", "1"]
            reference, reference_log = execute(f"stdin/{program}-{case}/native", [*prefixes["native"], *common], expected=None, stdin=stdin)
            reference_exit = rows[-1]["exit"]
            for kind in (kind for kind in ["serial", "threaded"] if kind in prefixes):
                data, log = execute(f"stdin/{program}-{case}/{kind}", [*prefixes[kind], *common], expected=None, stdin=stdin)
                assert (data, log, rows[-1]["exit"]) == (reference, reference_log, reference_exit), \
                    f"stdout, stderr or exit status differs from native: stdin/{program}-{case}/{kind}"
        print(f"stdin/{program}", "native bytes PASS", flush=True)
    import shutil
    shutil.rmtree(rps_db.parent)
    # NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:1199-1201
    # makeblastdb -logfile through each WASI command. First every makeblastdb -logfile row of the drop-in corpus
    # against its frozen NCBI bytes, by the checker of the native gates (dropin_corpus.py check: the exit status,
    # stdout, a quiet stderr, the log file and the database files, without the lines of the time and the path),
    # through a wrapper that starts the command's runner in the place of the binary, as the psiblast rows above.
    # Then two runs that append to one log file, the second with a warning and an error between the log's lines:
    # each run's stdout, stderr and exit status, and the log file's bytes without the lines of the time and the
    # path (dropin_corpus.report_without_time), must be the native command's. LOSAT opens the log a second time
    # under the name as written, beside the diagnostics' file (LOSAT/src/writedb/makeblastdb.rs `log_of`, line WQA's
    # D11 round 2, 2026-11-20); a copy of the diagnostics' handle (`try_clone`), which WASI cannot make, stopped
    # every such run with 'Argument "logfile". File is not accessible'.
    import dropin_corpus
    logfile_rows = [row["row_id"] for row in dropin_corpus.read_manifest()[1]
                    if shlex.split(row["command"])[0] == "makeblastdb" and "-logfile" in shlex.split(row["command"])]
    assert logfile_rows, "no makeblastdb -logfile row in the drop-in corpus"
    writedb_inputs = TESTS / "writedb_v5"
    logfile_runs = [(["-in", writedb_inputs / "a.fa"], "db"),
                    (["-in", f"{writedb_inputs / 'd.txt'} {writedb_inputs / 'a.fa'}"], "db2")]
    def logfile_case(kind):
        directory = out / "makeblastdb-logfile" / kind / "append"
        directory.mkdir(parents=True, exist_ok=True)
        log_path = directory / "mk.log"
        log_path.unlink(missing_ok=True)
        results = []
        for index, (words, name) in enumerate(logfile_runs):
            data, log = execute(f"makeblastdb-logfile/{kind}/append-{index}",
                                [*prefixes[kind], "makeblastdb", *words, "-dbtype", "prot", "-out", directory / name,
                                 "-title", "t", "-logfile", log_path], expected=None)
            results.append((data, log, rows[-1]["exit"]))
        return results, log_path.read_bytes()
    reference_runs, reference_log = logfile_case("native")
    assert [run[2] for run in reference_runs] == [0, 1] and reference_log.count(b"Building a new DB") == 2 \
        and b"Warning: [makeblastdb]" in reference_log, (reference_runs, reference_log)
    for kind in (kind for kind in ["serial", "threaded"] if kind in prefixes):
        wrapper = out / "makeblastdb-logfile" / kind / "losat"
        wrapper.parent.mkdir(parents=True, exist_ok=True)
        wrapper.write_text("#!/bin/sh\nexec " + " ".join(shlex.quote(str(x)) for x in prefixes[kind]) + ' "$@"\n')
        wrapper.chmod(0o755)
        for row_id in logfile_rows:
            label = f"makeblastdb-logfile/{kind}/{row_id}"
            try:
                execute(label, [sys.executable, TESTS / "dropin_corpus.py", "check", "--losat", wrapper, "--only", row_id],
                        timeout=600)
            except AssertionError:
                print((out / label / "stdout").read_text(errors="replace"), flush=True)
                raise
            lines = (out / label / "stdout").read_text().splitlines()
            assert lines == [f"{row_id}\tsame", "1 / 1 same"], (label, lines)
        runs, log = logfile_case(kind)
        assert runs == reference_runs, f"stdout, stderr or exit status differs from native: makeblastdb-logfile/{kind}"
        assert dropin_corpus.report_without_time(log) == dropin_corpus.report_without_time(reference_log), \
            f"the -logfile file differs from native: makeblastdb-logfile/{kind}"
        print(f"makeblastdb-logfile/{kind}", f"{len(logfile_rows)} drop-in rows same, appended log = native", flush=True)
    for program, query in [("blastn", "nuc1"), ("tblastx", "nuc1"), ("blastp", "aa1")]:
        common = [program, "-query", inputs[query], "-subject", inputs[query], "-outfmt", "6", "-out", "{out}"]
        for kind, prefix in prefixes.items():
            for value in ["0", "invalid", "2147483647"] + (["2", "4"] if kind in ["serial", "native-serial"] else []):
                data, log = execute(f"reject/{program}-{kind}-{value}", [*prefix, *common, "-num_threads", value], {"LOSAT_WASI_THREADS_DEBUG": "1"}, expected=1)
                assert not data and "spawn_attempt" not in log and "[losat-thread-pool]" not in log
            for cap, n in [("0", 1), ("invalid", 1), ("2", 4)]:
                data, log = execute(f"cap/{program}-{kind}-{cap}", [*prefix, *common, "-num_threads", str(n)], {"LOSAT_WASI_THREAD_CAP": cap, "LOSAT_WASI_THREADS_DEBUG": "1"}, expected=1)
                assert not data and "LOSAT_WASI_THREAD_CAP" in log and "spawn_attempt" not in log
        # stdout is the same sink content as the file route.
        stdout_args = common[:-2] + ["-num_threads", "1"]
        reference, _ = execute(f"stdout/{program}-native", [*prefixes["native"], *stdout_args])
        # NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1098-1108
        for kind in (kind for kind in ["serial", "threaded"] if kind in prefixes):
            data, _ = execute(f"stdout/{program}-{kind}", [*prefixes[kind], *stdout_args])
            assert data == reference
        # NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:727-729
        # A file argument is opened by its name through the OS: -query /dev/stdin and -out
        # /dev/stdout with pipes there (uvwasi gave ENOENT for both, at every thread count).
        if program == "blastp":
            pipe_args = [program, "-query", "/dev/stdin", "-subject", inputs[query], "-outfmt", "6", "-out", "/dev/stdout"]
            query_bytes = Path(inputs[query]).read_bytes()
            for kind, n in [("native", 1), ("serial", 1), ("threaded", 1), ("threaded", 2)]:
                if kind not in prefixes: continue
                data, log = execute(f"dev-stdio/{program}-{kind}-n{n}", [*prefixes[kind], *pipe_args, "-num_threads", str(n)],
                                    {"LOSAT_WASI_THREADS_DEBUG": "1"}, stdin=query_bytes)
                assert data == reference, f"raw output differs: dev-stdio/{program}-{kind}-n{n}"
                validate_thread_evidence(log, n, kind)
            # NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:211-213
            # LOSAT tells and seeks a FASTA input to find it empty (LOSAT/src/blastinput/input_files.rs);
            # /dev/null seeks, as native lseek does there: the warning and no report.
            empty_args = [program, "-query", "/dev/null", "-subject", inputs[query], "-outfmt", "7"]
            empty, empty_log = execute(f"dev-null/{program}-native", [*prefixes["native"], *empty_args, "-num_threads", "1"])
            assert "Query is Empty!" in empty_log, empty_log
            for kind, n in [("serial", 1), ("threaded", 1), ("threaded", 2)]:
                if kind not in prefixes: continue
                data, log = execute(f"dev-null/{program}-{kind}-n{n}", [*prefixes[kind], *empty_args, "-num_threads", str(n)])
                assert data == empty and log == empty_log, f"output or messages differ: dev-null/{program}-{kind}-n{n}"
            # No NCBI counterpart (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES, exception 8): a -query or
            # -subject named .gz, .bz2 or .zst is inflated beside the input before it is read
            # (LOSAT/src/blastinput/input_files.rs), so the report is the plain file's. Rust's
            # temp_dir() panics on the WASI targets; the commands must not reach it.
            import bz2, gzip
            packed_dir = out / "compressed-inputs"
            packed_dir.mkdir(parents=True, exist_ok=True)
            plain = Path(inputs[query]).read_bytes()
            packed = {".gz": gzip.compress(plain), ".bz2": bz2.compress(plain),
                      ".zst": subprocess.run(["zstd", "-q", "-c"], input=plain, capture_output=True, check=True).stdout}
            plain_args = [program, "-subject", inputs[query], "-outfmt", "6", "-out", "{out}"]
            reference, _ = execute(f"compressed/{program}-native-plain", [*prefixes["native"], *plain_args, "-query", inputs[query]])
            for suffix, data in packed.items():
                name = packed_dir / f"{query}.fasta{suffix}"
                name.write_bytes(data)
                for kind, n in [("native", 1), ("serial", 1), ("threaded", 1), ("threaded", 4)]:
                    if kind not in prefixes: continue
                    for role in ["query", "subject"]:
                        roles = {"query": inputs[query], "subject": inputs[query], role: name}
                        result, _ = execute(f"compressed/{program}-{kind}-n{n}-{role}{suffix}",
                                            [*prefixes[kind], program, "-query", roles["query"], "-subject", roles["subject"],
                                             "-outfmt", "6", "-out", "{out}", "-num_threads", str(n)])
                        assert result == reference, f"raw output differs: compressed/{program}-{kind}-n{n}-{role}{suffix}"
    # The organism report (-outfmt 18) is written since line RG T2 (c86dd36c1) and is compared
    # with the oracle like the other formats; no output format is refused in these programs.
    format_failures = []
    for program, query, subject, formats, task in [
        ("blastx", "blastx-query", "aa3", ["0", "6", "7", "6 std qframe sframe btop"], None),
        ("blastp", "aa3", "aa3", ["0", "7", "6 std qlen slen positive ppos btop stitle"], None),
        ("blastp", "aa-query-edge", "aa-subject-edge", ["0", "7", "6 std qlen slen positive ppos btop stitle", "7 qseqid qacc qaccver sseqid sacc saccver qlen slen score nident positive gaps frames qframe sframe qseq sseq btop stitle"], None),
        ("blastn", "nuc3", "nuc3", ["0", "7", "18"], "megablast"),
        ("blastn", "nuc3", "nuc3", ["0", "7"], "blastn"),
        ("tblastx", "nuc3", "nuc3", ["0", "7", "18"], None),
        ("tblastn", "aa1", "tblastn-subjects", ["0", "7"], None),
    ]:
        for index, fmt in enumerate(formats):
            common = ["-query", inputs[query], "-subject", inputs[subject], "-outfmt", fmt, "-out", "{out}"]
            # NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/cmdline_flags.cpp:69
            # Exercise both nucleotide tasks through their supported formatter.
            if task: common += ["-task", task]
            format_label = task or program
            oracle = str(args.oracle_dir.resolve() / program) if args.oracle_dir else program
            reference, _ = execute(f"format/{format_label}-{query}-{index}/oracle", [oracle, *common])
            # NCBI's output in LOSAT's words: its version line and prefixes (AGENTS.md rule 8).
            reference = losat_wording.losat_output_bytes(reference)
            for kind, prefix in prefixes.items():
                for n in ([1] if kind in ["serial", "native-serial"] else [1, 4, 8, 16]):
                    data, _ = execute(f"format/{format_label}-{query}-{index}/{kind}-n{n}", [*prefix, program, *common, "-num_threads", str(n)])
                    if data != reference:
                        format_failures.append(dict(program=program, task=task, format=fmt, target=kind, threads=n, oracle_sha256=digest(reference), actual_sha256=digest(data)))
                        (out / "format-failures.json").write_text(json.dumps(format_failures, indent=2) + "\n")
    for program, query in [("blastn", "nuc1"), ("tblastx", "nuc3"), ("blastp", "aa1")]:
        target = out / "output-is-directory" / program
        target.mkdir(parents=True, exist_ok=True)
        for kind, prefix in prefixes.items():
            n = 1 if kind in ["serial", "native-serial"] else 4
            _, log = execute(f"output-error/{program}-{kind}", [*prefix, program, "-query", inputs[query], "-subject", inputs[query], "-outfmt", "6", "-num_threads", str(n), "-out", target], expected=1)
            assert "directory" in log.lower(), log
    # NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/build_db.cpp:1102-1107
    # makeblastdb's local time and absolute database name through each WASI command, in a working directory of
    # their own, with the zone of TZ as a zone file (Asia/Tokyo), as a POSIX rule (JST-9) and unset (the host's
    # /etc/localtime): a new database, the same again ("Deleted existing ... named"), one in two new directories,
    # one whose directory is a file (an error), blastdbcmd -info of two and a blastn report of one (the database's
    # date, the volume paths). Each run's stdout, stderr and exit status must be the native command's
    # of the same directory and zone, with the times put aside (`local_times`); each time must be the local time
    # of its build in the zone that the C library finds for the run's TZ (Python's time.localtime; the clock moves
    # between the runs, LOSAT has no fixed clock to test with); the files must have the native names. LOSAT finds
    # the zone where glibc finds it on every target but Windows (the crate jiff reads it), on WASI from the
    # environment and the directories the host gives (LOSAT/src/writedb/local_offset.rs), and makes a relative name absolute with wasi-libc's working
    # directory, which the command sets from the host's PWD (LOSAT/src/writedb/build_db.rs
    # `create_absolute_path`). Before, the WASI commands printed UTC and the name as given.
    def libc_offset(tz):
        probe = "import time; print(time.localtime().tm_gmtoff)"
        run_env = {k: v for k, v in {**environment, "TZ": tz}.items() if v is not None}
        return int(subprocess.check_output([sys.executable, "-c", probe], env=run_env, text=True))
    makedb = ["makeblastdb", "-in", inputs["nuc3"], "-dbtype", "nucl"]
    local_steps = [
        # label, the words after the program, the step that built the database whose date the output prints
        ("new", [*makedb, "-out", "db", "-title", "t"], "new"),
        ("again", [*makedb, "-out", "db", "-title", "t"], "again"),
        ("subdir", [*makedb, "-out", "sub/new/db2"], "subdir"),
        ("file-as-dir", [*makedb, "-out", "file/db3"], None),
        ("info", ["blastdbcmd", "-db", "db", "-info"], "again"),
        ("info-subdir", ["blastdbcmd", "-db", "sub/new/db2", "-info"], "subdir"),
        ("report", ["blastn", "-query", inputs["nuc1"], "-db", "db", "-outfmt", "0", "-num_threads", "1"], "again"),
    ]
    for zone, tz, fixed_offset in [("zone-file", "Asia/Tokyo", 32400), ("posix-rule", "JST-9", 32400),
                                   ("host", None, None)]:
        offset = libc_offset(tz)
        assert fixed_offset is None or offset == fixed_offset, (zone, tz, offset)
        work = out / "local-time-and-path" / zone / "work"
        reference = None
        for kind in prefixes:
            shutil.rmtree(work, ignore_errors=True)
            work.mkdir(parents=True)
            (work / "file").write_text("not a directory\n")
            results, built = [], {}
            for step, words, date_of in local_steps:
                started = time.time()
                data, log = execute(f"local-time-and-path/{zone}/{kind}/{step}", [*prefixes[kind], *words],
                                    {"TZ": tz}, expected=None, timeout=120, cwd=work)
                ended = time.time()
                if words[0] == "makeblastdb":
                    built[step] = (started, ended)
                data, printed = local_times(data)
                for value, resolution in printed:
                    at = (started, ended) if resolution == 1 else built[date_of]
                    assert is_local_time_of(value, resolution, offset, *at), \
                        f"not the local time of TZ={tz} (offset {offset} s): local-time-and-path/{zone}/{kind}/{step}"
                results.append((step, data, log, rows[-1]["exit"], len(printed)))
            results.append(("files", sorted(str(path.relative_to(work)) for path in work.rglob("*")), "", 0, 0))
            if reference is None:
                assert kind == "native", list(prefixes)
                reference = results
                # makeblastdb's exception ("Failed to create directory 'file'") exits with 255 before its log;
                # each other run prints one time.
                directory = str(work).encode()
                assert [(r[3], r[4]) for r in results[:-1]] == [(0, 1), (0, 1), (0, 1), (255, 0), (0, 1), (0, 1), (0, 1)] \
                    and b"New DB name:   " + directory + b"/db\n" in results[0][1] \
                    and b"named " + directory + b"/db\n" in results[1][1] \
                    and b"Volumes:\n\t" + directory + b"/sub/new/db2\n" in results[5][1] \
                    and b"Failed to create directory 'file'" in results[3][2].encode(), [r[0::3] for r in results]
                continue
            for got, want in zip(results, reference):
                assert got == want, f"stdout, stderr, exit status or files differ from native: " \
                                    f"local-time-and-path/{zone}/{kind}/{got[0]}"
        print(f"local-time-and-path/{zone}", f"TZ={tz} offset {offset} s", "native bytes PASS", flush=True)
    execute("reactor", [args.node, TESTS / "check_wasi_reactor.js", args.reactor.resolve(), out / "fixtures", out / "reactor-records"], {"LOSAT_WASI_THREADS_DEBUG": "1"}, timeout=90)
    for fault in ["startup-timeout", "startup-error", "trap", "abnormal-exit"]:
        _, log = execute("fault/" + fault, [args.node, TESTS / "check_wasi_fault.js", args.reactor.resolve(), inputs["nuc1"], fault], {"LOSAT_WASI_THREADS_DEBUG": "1"}, expected=1, timeout=15)
        assert rows[-1]["exit"] == -15, (fault, rows[-1])
        assert "timed out waiting" in log or "injected" in log or "exited unexpectedly" in log
    if args.serial_reactor:
        execute("serial-reactor", [args.node, TESTS / "check_wasi_api_limits.js", args.serial_reactor.resolve(), "serial-reactor", out / "fixtures", out / "serial-reactor-records"])
    for cap in ["0", "invalid", "1"]:
        execute("reactor-cap/" + cap, [args.node, TESTS / "check_wasi_api_limits.js", args.reactor.resolve(), "threaded-reactor", out / "fixtures", out / ("reactor-cap-" + cap)], {"LOSAT_WASI_THREAD_CAP": cap, "LOSAT_WASI_THREADS_DEBUG": "1"})
    # Artifact misclassification is rejected before entry/worker startup.
    # NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:172-176
    wrong_artifacts = [("threaded-reactor", "run_losat_wasi_threads.js", args.reactor)]
    if args.serial:
        wrong_artifacts += [("serial-reactor", "run_losat_wasi.js", args.reactor),
                            ("threaded-serial", "run_losat_wasi_threads.js", args.serial),
                            ("serial-threaded", "run_losat_wasi.js", args.threaded)]
    for label, runner, artifact in wrong_artifacts:
        _, log = execute("wrong-artifact/" + label, [args.node, TESTS / runner, artifact.resolve(), "--help"], {"LOSAT_WASI_THREADS_DEBUG": "1"}, expected=1)
        assert "spawn_attempt" not in log
    print(f"{len(rows)} command/oracle records; reactor lifecycle gates PASS; format failures={len(format_failures)}", flush=True)
    assert not format_failures, "NCBI format comparison failed; see format-failures.json (no normalization or exception applied)"


if __name__ == "__main__":
    main()
