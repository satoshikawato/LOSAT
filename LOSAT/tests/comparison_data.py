"""Shared inputs and historical filenames for the simple comparison plots."""

import csv
import hashlib
import json
import shutil
import subprocess
import sys
import time
import uuid
import os
from pathlib import Path
import re


SCRIPT_DIR = Path(__file__).resolve().parent
RESULT_DIR = SCRIPT_DIR / os.environ.get("BENCHMARK_DIR", ".")
PLOT_DIR = RESULT_DIR / "plots"

# NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:75
# const string kArgNumThreads("num_threads");
LOSAT_THREADS = int(
    os.environ.get("LOSAT_THREADS")
    or os.environ.get("LOSATP_THREADS")
    or os.environ.get("LOSAT_BLASTP_THREADS")
    or "8"
)
if LOSAT_THREADS < 1:
    raise ValueError("LOSAT_THREADS must be positive")
NATIVE_SINGLE = "LOSAT native n1"
NATIVE_MULTI = f"LOSAT native n{LOSAT_THREADS}"
WASM_SINGLE = "LOSAT wasm serial n1"
WASM_MULTI = f"LOSAT wasm threads n{LOSAT_THREADS}"
HUE_ORDER = list(dict.fromkeys([
    "BLAST+", NATIVE_SINGLE, NATIVE_MULTI, WASM_SINGLE, WASM_MULTI,
]))
CUSTOM_PALETTE = dict(zip(
    ["BLAST+", NATIVE_SINGLE, NATIVE_MULTI, WASM_SINGLE, WASM_MULTI],
    ["#4c72b0", "#dd8452", "#a15c2e", "#8a8f3b", "#565b22"],
))
MODE_ORDER = ["TBLASTX", "Megablast", "BLASTN", "BLASTP"]

# NCBI reference: c++/src/objtools/align_format/format_flags.cpp:38-40
# const char* kDfltArgTabularOutputFmt =
#     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
#     "evalue bitscore";
COLUMNS = "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore".split()


def comparison_cases():
    programs = os.environ.get("BENCHMARK_PROGRAMS", "tblastx,megablast,blastn,blastp").split(",")
    if set(programs) - {"tblastx", "megablast", "blastn", "blastp"}:
        raise ValueError(f"Unknown BENCHMARK_PROGRAMS: {programs}")
    case_filter = os.environ.get("BENCHMARK_CASE", "")
    with (SCRIPT_DIR / "comparison_cases.tsv").open() as handle:
        cases = [row for row in csv.DictReader(handle, delimiter="\t")
                 if row["task"] in programs and case_filter in row["losat_stem"]]
    if not cases:
        raise ValueError("No matching comparison cases")
    for case in cases:
        case["mode"] = "Megablast" if case["task"] == "megablast" else case["task"].upper()
        case["group"] = "BLASTN (All Types)" if case["task"] in {"blastn", "megablast"} else case["mode"]
    return cases


def result_paths(case, extension="out"):
    """Keep n1/native, nN/native and serial/threaded Wasm distinct, including N=1."""
    native = RESULT_DIR / "losat_out" / case["losat_stem"]
    ncbi = RESULT_DIR / "blast_out" / case["ncbi_stem"]
    single = str(native)
    if case["task"] == "tblastx":
        single += ".n1"
        ncbi = f"{ncbi}.n{LOSAT_THREADS}"
    paths = {"BLAST+": ncbi, NATIVE_SINGLE: single}
    if LOSAT_THREADS != 1:
        paths[NATIVE_MULTI] = f"{native}.n{LOSAT_THREADS}"
    paths[WASM_SINGLE] = f"{native}.wasm"
    paths[WASM_MULTI] = f"{native}.wasm.n{LOSAT_THREADS}"
    return {tool: Path(f"{stem}.{extension}") for tool, stem in paths.items()}


# NCBI reference: c++/src/objtools/align_format/tabular.cpp:1100-1108
# ITERATE(list<ETabularField>, iter, m_FieldsToShow) { ... x_PrintField(*iter); }
# m_Ostream << "\\n";
# Evidence binds successful status to the exact ordered output of this invocation.
def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def completed_log(content):
    return re.findall(r"^exit_status=(\d+)$", content, re.M) == ["0"]


def successful_output(path):
    path = Path(path)
    try:
        log = path.with_suffix(".log")
        record = json.loads(path.with_suffix(".run.json").read_text())
        manifest = json.loads((path.parent.parent / "run.json").read_text())
        argv_hash = hashlib.sha256(json.dumps(record["ordered_argv"]).encode()).hexdigest()
        return (
            record["run_id"] == manifest["run_id"]
            and record["status"] == "PASS" and record["exit_status"] == 0
            and record["output"] == str(path.resolve())
            and record["argv_sha256"] == argv_hash
            and bool(record["ordered_argv"]) and bool(record["files"])
            and completed_log(log.read_text())
            and record["log_sha256"] == sha256(log)
            and record["output_sha256"] == sha256(path)
        )
    except (OSError, ValueError, KeyError, TypeError):
        return False


# NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-48,75
# const string kArgQuery("query"); const string kArgSubject("subject");
# const string kArgOutput("out"); const string kArgNumThreads("num_threads");
# Harness metadata never supplies implementation behavior or expected bytes.
def record_cli(argv):
    action, directory, *rest = argv
    directory = Path(directory).resolve()
    if action == "init":
        directory.mkdir(parents=True, exist_ok=False)
        (directory / "run.json").write_text(json.dumps({
            "schema": "losat-simple-comparison-v2", "run_id": uuid.uuid4().hex,
            "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "node_argv": [os.environ.get("NODE_BIN", "node"), *node_args()],
            # NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
            # const string kArgQuery("query"); const string kArgSubject("subject");
            # Compilation settings are host argv, separate from unchanged search argv.
            "node_argv_by_program": {program: [os.environ.get("NODE_BIN", "node"), *node_args(program)]
                                     for program in ("tblastx", "blastn", "blastp")},
            "source_hashes": {str(p.relative_to(SCRIPT_DIR.parent)):sha256(p)
                              for p in [*sorted((SCRIPT_DIR.parent/"src").rglob("*.rs")), SCRIPT_DIR.parent/"Cargo.toml", SCRIPT_DIR.parent/"Cargo.lock", SCRIPT_DIR.parent/"build.rs", SCRIPT_DIR.parent/".cargo/config.toml"]},
            "toolchain": {name: subprocess.check_output([name, "--version"], text=True).strip()
                          for name in ("rustc", "cargo") if shutil.which(name)},
            "node_versions": (subprocess.check_output([os.environ.get("NODE_BIN", "node"), "-p", "JSON.stringify(process.versions)"], text=True).strip()
                              if os.environ.get("RUN_LOSAT_WASM", "1") == "1" else None),
            "environment": {k:v for k,v in os.environ.items() if k.startswith(("LOSAT_", "RAYON_", "NODE_", "BENCHMARK_", "RUN_", "BUILD_")) or k in {"LC_ALL", "BL2SEQ_LEGACY"}},
            "runner_hashes": {p.name: sha256(p) for p in SCRIPT_DIR.iterdir()
                              if p.suffix in {".py", ".js", ".sh", ".tsv"}},
        }, indent=2) + "\n")
        (directory / "node-args.bin").write_bytes(b"".join(x.encode() + b"\0" for x in node_args()))
        # NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
        # const string kArgQuery("query"); const string kArgSubject("subject");
        # Preserve argument boundaries for the program-specific host profile.
        (directory / "node-tblastx-args.bin").write_bytes(b"".join(x.encode() + b"\0" for x in node_args("tblastx")))
        return
    stem = Path(rest[0]).resolve()
    manifest = json.loads((directory / "run.json").read_text())
    record_path = stem.with_suffix(stem.suffix + ".run.json")
    if action == "database":
        source, executable = rest[1:]
        Path(str(stem) + ".makeblastdb.json").write_text(json.dumps({
            "run_id": manifest["run_id"], "input": source, "input_sha256": sha256(source),
            "executable": executable, "executable_sha256": sha256(executable),
            "version": subprocess.check_output([executable, "-version"], text=True).strip(),
            "ordered_argv": [executable, "-in", source, "-dbtype", "nucl", "-parse_seqids", "-out", str(stem)],
            "exit_status": 0,
        }, indent=2) + "\n")
        return
    if action == "start":
        command = rest[1:]
        output = Path(str(stem) + ".out")
        if output.exists() or record_path.exists():
            raise FileExistsError(f"refusing stale output/record: {stem}")
        files = {str(Path(x).resolve()): sha256(x) for x in command if Path(x).is_file()}
        executable = shutil.which(command[0])
        if executable:
            files[str(Path(executable).resolve())] = sha256(executable)
        if "-db" in command:
            db = Path(command[command.index("-db") + 1])
            files.update({str(p.resolve()):sha256(p) for p in db.parent.glob(db.name + ".*") if p.is_file()})
        record = {"run_id": manifest["run_id"], "status": "RUNNING",
                  "ordered_argv": command, "output": str(output), "cwd": str(Path.cwd()),
                  "argv_sha256": hashlib.sha256(json.dumps(command).encode()).hexdigest(),
                  # NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-48
                  # const string kArgQuery("query"); const string kArgSubject("subject");
                  # Identify the selected comparison oracle before its search invocation.
                  "oracle_version": (subprocess.check_output([executable, "-version"], text=True).strip()
                                     if executable and Path(command[0]).name in {"blastn", "blastp", "tblastx"} else None),
                  "files": files, "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    elif action == "finish":
        record = json.loads(record_path.read_text())
        status = int(rest[1]); output = Path(record["output"])
        record.update(exit_status=status, status="PASS" if status == 0 and output.is_file() else "FAIL",
                      output_sha256=sha256(output) if output.is_file() else None,
                      log_sha256=sha256(str(stem) + ".log"),
                      ended_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    else:
        raise ValueError(f"unknown action: {action}")
    record_path.write_text(json.dumps(record, indent=2) + "\n")


# NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
# const string kArgQuery("query"); const string kArgSubject("subject");
# Host compilation does not change search arguments or NCBI behavior. The user
# approved the measured TBLASTX memory tradeoff; other programs retain their flags.
def node_args(program=None):
    settings = [("NODE_ARGS_JSON", "[]")]
    if program == "tblastx":
        settings.append(("NODE_TBLASTX_ARGS_JSON", '["--no-liftoff", "--no-wasm-tier-up"]'))
    args = []
    for name, default in settings:
        values = json.loads(os.environ.get(name, default))
        if not isinstance(values, list) or any(not isinstance(x, str) or "\0" in x for x in values):
            raise ValueError(f"{name} must be a JSON array of strings without NUL")
        args.extend(values)
    return args


def load_case(case):
    """Load available outputs, including valid zero-hit files, paired with NCBI."""
    import pandas as pd

    paths = result_paths(case)
    if not successful_output(paths["BLAST+"]):
        print(f"[Skip] {case['name']} ({case['mode']}): missing/failed NCBI output")
        return {}
    frames = {}
    for tool, path in paths.items():
        if not successful_output(path):
            print(f"[Missing] {case['name']} ({case['mode']}): {tool}")
            continue
        try:
            frame = pd.read_csv(path, sep="\t", comment="#", header=None)
        except pd.errors.EmptyDataError:
            frame = pd.DataFrame(columns=range(len(COLUMNS)))
        if len(frame.columns) != len(COLUMNS):
            raise ValueError(f"Expected 12 outfmt 6/7 columns: {path}")
        frame.columns = COLUMNS
        for column in COLUMNS[2:]:
            frame[column] = pd.to_numeric(frame[column], errors="raise")
        if frame[COLUMNS[2:]].isna().any().any():
            raise ValueError(f"Missing numeric output field: {path}")
        frame["Tool"] = tool
        frame["Mode"] = case["mode"]
        frame["Task"] = case["name"]
        frame["Broad_Mode"] = case["group"]
        frames[tool] = frame
    return frames if len(frames) > 1 else {}


if __name__ == "__main__":
    record_cli(sys.argv[1:])
