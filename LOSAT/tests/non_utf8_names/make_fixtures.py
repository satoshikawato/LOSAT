#!/usr/bin/env python3
"""Freeze what NCBI BLAST+ 2.17.0 writes for the cases of `cases.tsv` (`expected/<case>.stdout`,
`.stderr`, `.rc`), the cases of `LOSAT/tests/non_utf8_names.rs`.

Each case runs in a new directory of its own (`--work`), made by the case's setup, with an empty
`HOME`, without `BLASTDB` and `NCBI_CONFIG_PATH` but the case's, and with `BLAST_USAGE_REPORT=false`
(NCBI BLAST+ sends no report of the run); the directory's path is written as `{D}` in the frozen
bytes, which the test replaces with its own. The manifest's words are bytes: `\\xHH` is the byte HH,
`\\n` a line feed, `\\s` a space, `\\\\` a backslash; `{D}` in the environment and in written files is
the case's directory. Setup steps, separated by ` ; `:
  copy SRC DST          the file SRC (relative to LOSAT/) as DST
  copydb SRCBASE DIR    every file SRCBASE.* into the directory DIR, under its own name
  mkdir DIR             a directory
  write PATH TEXT       a file of TEXT
  patch PATH OLD NEW    the first OLD of the file replaced by NEW (the same length)

usage: make_fixtures.py [--ncbi-bin DIR] [--work DIR] [case ...]
"""
import argparse
import os
import pathlib
import shutil
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
LOSAT = HERE.parent.parent
CLEARED = ("BLASTDB", "NCBI_CONFIG_PATH", "NCBI_DONT_USE_LOCAL_CONFIG", "NCBI", "BATCH_SIZE")


def words(text: str, d: bytes) -> bytes:
    """The bytes of a manifest word (escapes as in the module's text; `{D}` is `d`)."""
    out = bytearray()
    raw = text.encode("ascii")
    i = 0
    while i < len(raw):
        c = raw[i]
        if c == 0x5C:  # backslash
            n = raw[i + 1]
            if n == ord("x"):
                out.append(int(raw[i + 2 : i + 4], 16))
                i += 4
                continue
            out += {ord("n"): b"\n", ord("s"): b" ", 0x5C: b"\\"}[n]
            i += 2
            continue
        out.append(c)
        i += 1
    return bytes(out).replace(b"{D}", d)


def cases():
    for line in (HERE / "cases.tsv").read_text(encoding="ascii").splitlines():
        if not line or line.startswith("#"):
            continue
        case, kind, env, setup, argv = line.split("\t")
        yield case, kind, env, setup, argv


def set_up(setup: str, d: bytes) -> None:
    for step in setup.split(" ; "):
        op, *args = step.split(" ")
        if op == "copy":
            shutil.copyfile(LOSAT / args[0], os.path.join(d, words(args[1], d)))
        elif op == "copydb":
            src = LOSAT / args[0]
            for f in sorted(src.parent.glob(src.name + ".*")):
                shutil.copyfile(f, os.path.join(d, words(args[1], d), f.name.encode()))
        elif op == "mkdir":
            os.mkdir(os.path.join(d, words(args[0], d)))
        elif op == "write":
            with open(os.path.join(d, words(args[0], d)), "wb") as out:
                out.write(words(args[1], d))
        elif op == "patch":
            path = os.path.join(d, words(args[0], d))
            old, new = words(args[1], d), words(args[2], d)
            assert len(old) == len(new)
            data = open(path, "rb").read()
            assert old in data, (path, old)
            open(path, "wb").write(data.replace(old, new, 1))
        else:
            raise SystemExit(f"unknown step {step!r}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ncbi-bin", default=os.environ.get("NCBI_BIN", ""))
    parser.add_argument("--work", default=os.path.join(os.environ.get("TMPDIR", "/tmp"), "non_utf8_names"))
    parser.add_argument("only", nargs="*")
    args = parser.parse_args()
    work = pathlib.Path(args.work).resolve()
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    version = subprocess.run([os.path.join(args.ncbi_bin, "blastp"), "-version"], capture_output=True, check=True)
    assert b"2.17.0" in version.stdout, version.stdout
    for case, kind, env, setup, argv in cases():
        if args.only and case not in args.only:
            continue
        d = os.fsencode(work / case)
        os.mkdir(d)
        home = work / f"{case}_home"
        home.mkdir()
        set_up(setup, d)
        run_env = {k: v for k, v in os.environ.items() if k not in CLEARED and not k.startswith("LOSAT_")}
        run_env.update(HOME=str(home), BLAST_USAGE_REPORT="false")
        run_env = {os.fsencode(k): os.fsencode(v) for k, v in run_env.items()}
        if env != "-":
            name, value = env.split("=", 1)
            run_env[name.encode()] = words(value, d)
        program, *rest = argv.split(" ")
        command = [os.fsencode(os.path.join(args.ncbi_bin, program))] + [words(w, d) for w in rest]
        done = subprocess.run(command, cwd=d, env=run_env, capture_output=True, timeout=60)
        for name, data in (("stdout", done.stdout), ("stderr", done.stderr)):
            (HERE / "expected" / f"{case}.{name}").write_bytes(data.replace(d, b"{D}"))
        (HERE / "expected" / f"{case}.rc").write_text(f"{done.returncode}\n")
        print(case, kind, "rc", done.returncode, "stdout", len(done.stdout), "stderr", len(done.stderr))
    return 0


if __name__ == "__main__":
    sys.exit(main())
