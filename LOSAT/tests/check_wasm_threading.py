#!/usr/bin/env python3
"""Real native/WASI threading gates with raw NCBI 2.17.0 oracle output."""
import argparse
import hashlib
import json
import os
import platform
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
TESTS = ROOT / "LOSAT/tests"
sys.path.insert(0, str(TESTS))
from wasm_performance import validate_thread_evidence


def digest(data):
    return hashlib.sha256(data).hexdigest()


# NCBI reference: c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-180
# m_iTotalLength += (Int8) (*iter)->length;
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
    # NCBI reference: c++/src/objtools/align_format/tabular.cpp:415-439,1264-1283
    # Empty BLAST defline object prints N/A; every query receives its own header.
    edge_title = "a long title, with spaces and hyphens " * 5
    records["aa-query-edge"] = [aa, "W" * 85, aa]
    records["aa-subject-edge"] = [aa, aa]
    for name in ["aa-query-edge", "aa-subject-edge"]:
        (directory / f"{name}.fasta").write_text("".join(
            f">seq{i}" + (f" {edge_title}" if i == 0 else "") + f"\n{value}\n"
            for i, value in enumerate(records[name])))
    return {name: directory / f"{name}.fasta" for name in records}


# NCBI reference: c++/src/objtools/align_format/tabular.cpp:1098-1108
# x_PrintField(*iter); m_Ostream << "\n";
# Retain exact argv, stderr and outputs. No sorting or output normalization.
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--native-serial", type=Path)
    # NCBI reference: c++/include/algo/blast/blastinput/blast_args.hpp:1290-1296
    # #ifdef NCBI_NO_THREADS
    # m_NumThreads = CThreadable::kMinNumThreads; m_MTMode = eNotSupported;
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
    def execute(label, argv, env=None, expected=0, timeout=60):
        directory = out / label; directory.mkdir(parents=True, exist_ok=True)
        result_path = directory / "result.out"
        argv = [str(x).replace("{out}", str(result_path)) for x in argv]
        try:
            result = subprocess.run(argv, cwd=ROOT, env={**environment, **(env or {})}, capture_output=True, timeout=timeout)
        except subprocess.TimeoutExpired as error:
            (directory / "stdout").write_bytes(error.stdout or b"")
            (directory / "stderr").write_bytes(error.stderr or b"")
            rows.append(dict(label=label, argv=argv, cwd=str(ROOT), env=env or {}, status="TIMEOUT", timeout_seconds=timeout))
            (out / "runs.json").write_text(json.dumps(rows, indent=2) + "\n")
            raise

        (directory / "stdout").write_bytes(result.stdout); (directory / "stderr").write_bytes(result.stderr)
        data = result_path.read_bytes() if result_path.exists() else result.stdout
        row = dict(label=label, argv=argv, cwd=str(ROOT), env=env or {}, exit=result.returncode, sha256=digest(data))
        rows.append(row); (out / "runs.json").write_text(json.dumps(rows, indent=2) + "\n")
        assert (result.returncode == 0) == (expected == 0), (label, result.returncode, result.stderr.decode())
        return data, result.stderr.decode()
    prefixes = {
        "native": [args.native.resolve()],
        "threaded": [args.node, TESTS / "run_losat_wasi_threads.js", args.threaded.resolve()],
    }
    # NCBI reference: c++/include/algo/blast/blastinput/blast_args.hpp:1290-1296
    # m_NumThreads = CThreadable::kMinNumThreads; m_MTMode = eNotSupported;
    if args.serial: prefixes["serial"] = [args.node, TESTS / "run_losat_wasi.js", args.serial.resolve()]
    if args.native_serial: prefixes["native-serial"] = [args.native_serial.resolve()]
    cases = [("blastn", "blastn", "nuc1", "nuc1", 10), ("megablast", "blastn", "nuc1", "nuc3", 10),
             ("blastp", "blastp", "aa3", "aa3", 10)]
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
            for n in ([1] if kind in ["serial", "native-serial"] else [1, 2, 4, 8]):
                data, log = execute(f"{label}/{kind}-n{n}", [*prefix, program, *common, "-num_threads", str(n)], {"LOSAT_WASI_THREADS_DEBUG": "1"})
                assert data == reference, f"raw output differs: {label}/{kind}-n{n}"
                validate_thread_evidence(log, n, "native" if kind == "native-serial" else kind)
        print(label, "raw parity + thread contract PASS", flush=True)
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
        # NCBI reference: c++/src/objtools/align_format/tabular.cpp:1098-1108
        # x_PrintField(*iter); m_Ostream << "\n";
        for kind in (kind for kind in ["serial", "threaded"] if kind in prefixes):
            data, _ = execute(f"stdout/{program}-{kind}", [*prefixes[kind], *stdout_args])
            assert data == reference
    # NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2657-2660
    # AddDefaultKey(kArgOutputFormat, ..., eString, ...);
    # Unported format routes must reject explicitly before starting any worker.
    for program, formats in [("blastn", ["0", "6 qseqid sseqid"]),
                             ("tblastx", ["0", "7", "6 qseqid sseqid"])]:
        for index, fmt in enumerate(formats):
            for kind, prefix in prefixes.items():
                n = 1 if kind in ["serial", "native-serial"] else 8
                data, log = execute(f"unsupported-format/{program}-{index}-{kind}",
                    [*prefix, program, "-query", inputs["nuc1"], "-subject", inputs["nuc1"],
                     "-outfmt", fmt, "-num_threads", str(n), "-out", "{out}"],
                    {"LOSAT_WASI_THREADS_DEBUG": "1"}, expected=1)
                assert not data and ("unsupported" in log.lower() or "not implemented" in log.lower())
                assert "spawn_attempt" not in log and "[losat-thread-pool]" not in log
    format_failures = []
    for program, query, subject, formats, task in [
        ("blastp", "aa3", "aa3", ["0", "7", "6 std qlen slen positive ppos btop stitle"], None),
        ("blastp", "aa-query-edge", "aa-subject-edge", ["0", "7", "6 std qlen slen positive ppos btop stitle", "7 qseqid qacc qaccver sseqid sacc saccver qlen slen score nident positive gaps frames qframe sframe qseq sseq btop stitle"], None),
        ("blastn", "nuc3", "nuc3", ["7"], "megablast"),
        ("blastn", "nuc3", "nuc3", ["7"], "blastn"),
    ]:
        for index, fmt in enumerate(formats):
            common = ["-query", inputs[query], "-subject", inputs[subject], "-outfmt", fmt, "-out", "{out}"]
            # NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:69
            # const string kTask("task");
            # Exercise both nucleotide tasks through their supported formatter.
            if task: common += ["-task", task]
            format_label = task or program
            oracle = str(args.oracle_dir.resolve() / program) if args.oracle_dir else program
            reference, _ = execute(f"format/{format_label}-{query}-{index}/oracle", [oracle, *common])
            for kind, prefix in prefixes.items():
                for n in ([1] if kind in ["serial", "native-serial"] else [1, 4, 8]):
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
    # NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
    # CATCH_ALL(status); return status;
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
