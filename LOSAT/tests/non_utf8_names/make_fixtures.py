#!/usr/bin/env python3
"""Freeze what NCBI BLAST+ 2.17.0 writes for the cases of `cases.tsv` (`expected/<case>.stdout`,
`.stderr`, `.rc`, and `.files`: the files and directories the run made, changed or removed and
the contents of the files with the build's dates masked, `made()`), the cases of `LOSAT/tests/non_utf8_names.rs`.

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
  symlink TARGET NAME   a symbolic link NAME to TARGET
  cwd DIR               a directory, where the program runs (the case's directory otherwise)

usage: make_fixtures.py [--ncbi-bin DIR] [--work DIR] [case ...]
"""
import argparse
import os
import pathlib
import re
import shutil
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
LOSAT = HERE.parent.parent
CLEARED = (
    "BLASTDB",
    "NCBI_CONFIG_PATH",
    "NCBI_DONT_USE_LOCAL_CONFIG",
    "NCBI",
    "BATCH_SIZE",
    "WINDOW_MASKER_PATH",
)


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


def set_up(setup: str, d: bytes, makembindex_bin: str) -> bytes:
    """Runs the setup steps in `d`; the directory where the program runs."""
    cwd = d
    for step in setup.split(" ; "):
        op, *args = step.split(" ")
        if op == "copy":
            shutil.copyfile(LOSAT / args[0], os.path.join(d, words(args[1], d)))
        elif op == "copydir":
            src = LOSAT / args[0]
            dst = os.path.join(d, words(args[1], d))
            os.makedirs(dst, exist_ok=True)
            for f in sorted(src.iterdir()):
                if f.is_file():
                    shutil.copyfile(f, os.path.join(dst, f.name.encode()))
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
        elif op == "symlink":
            os.symlink(words(args[0], d), os.path.join(d, words(args[1], d)))
        elif op == "cwd":
            cwd = os.path.join(d, words(args[0], d))
            os.mkdir(cwd)
        elif op == "mkidx":
            cmd = [os.fsencode(makembindex_bin)] + [words(w, d) for w in args]
            subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, timeout=60)
        else:
            raise SystemExit(f"unknown step {step!r}")
    return cwd


def paths(d: bytes) -> set:
    """The paths below `d`, `/`-separated, a directory's with a final `/`."""
    found = set()
    for root, dirs, files in os.walk(d):
        rel = os.path.relpath(root, d)
        prefix = b"" if rel == b"." else rel + b"/"
        found.update(prefix + name + b"/" for name in dirs)
        found.update(prefix + name for name in files)
    return found


# The texts of the build's local time that a database's files hold, of a fixed length: NCBI's
# `CTime` form (`MM/DD/YYYY hh:mm:ss`: an alias file's comment, the index file of a column or of a
# GI mask) and the metadata file's (`YYYY-MM-DDThh:mm:ss`, `last-updated`).
DATE_TEXTS = re.compile(rb"\d\d/\d\d/\d{4} \d\d:\d\d:\d\d|\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d")

# A database's index files, whose header holds the build's date in a field of its own length.
INDEX = (b".nin", b".pin")


def index_date(data: bytes):
    """The span of the date field of an index file's header (its 4-byte length and its bytes
    with their padding: after the version, the type, in version 5 the volume, the title and in
    version 5 the LMDB file's name), or None for a file without such a header."""

    def int4(at):
        return int.from_bytes(data[at : at + 4], "big") if at + 4 <= len(data) else None

    version = int4(0)
    if version not in (4, 5):
        return None
    at = 12 if version == 5 else 8
    for _ in range(version - 3):
        length = int4(at)
        if length is None:
            return None
        at += 4 + length
    length = int4(at)
    if length is None or at + 4 + length > len(data):
        return None
    return at, at + 4 + length


def undated(d: bytes, path: bytes, data: bytes) -> bytes:
    """`data`, the bytes of the file `path` below `d`, with the build's dates as `<DATE>`: the
    date field of an index file (`index_date`) and every text of `DATE_TEXTS`; in a metadata file
    (`.njs`, `.pjs`) `bytes-total` counts the index files it lists without their date fields,
    whose length changes with the date; in a binary seqid list (`.bsl`) its date with the date's
    length byte, and its 8-byte size, which counts the date, as zeros when it holds the length of
    the file (a wrong size stays: `undated_seqid_list` of tests/temp_support)."""
    if path.endswith(b".bsl") and len(data) > 21 and data[0] == 0:
        date = 21 + int.from_bytes(data[17:21], "little")
        if date < len(data):
            size = bytes(8) if int.from_bytes(data[1:9], "little") == len(data) else data[1:9]
            data = data[:1] + size + data[9:date] + b"<DATE>" + data[date + 1 + data[date] :]
    if path.endswith(INDEX):
        span = index_date(data)
        if span:
            data = data[: span[0]] + b"<DATE>" + data[span[1] :]
    elif path.endswith((b".njs", b".pjs")):
        directory = os.path.dirname(os.path.join(d, path))
        listed = re.search(rb'"files": \[(.*?)\]', data, re.S)
        dates = 0
        for name in re.findall(rb'"([^"]*)"', listed.group(1)) if listed else []:
            if name.endswith(INDEX):
                idx_path = os.path.join(directory, name)
                if not os.path.exists(idx_path):
                    idx_path = os.path.join(directory, os.path.basename(name))
                with open(idx_path, "rb") as index:
                    span = index_date(index.read())
                dates += span[1] - span[0]
        for key in (b"bytes-total", b"bytes-to-cache"):
            data = re.sub(
                rb'"' + key + rb'": (\d+)',
                lambda m, k=key: b'"' + k + b'": %d' % (int(m.group(1)) - dates),
                data,
                count=1,
            )
    if path.endswith(b".log"):
        data = re.sub(rb"(Building a new DB, current time: )[^\n]*", rb"\1{TIME}", data, count=1)
        data = re.sub(rb"( sequences in )[^ ]*( seconds\.)", rb"\1{SECONDS}\2", data, count=1)
    if path.endswith((b".nal", b".pal")):
        data = re.sub(rb"(# Alias file created: )[^\n]*", rb"\1<DATE>", data, count=1)
    return DATE_TEXTS.sub(b"<DATE>", data)


def fnv1a64(data: bytes) -> int:
    h = 0xCBF29CE484222325
    for b in data:
        h = ((h ^ b) * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return h


def snapshot(d: bytes) -> dict:
    """The paths below `d` (`paths`), each with the bytes of its file (None for a directory, or a
    link to one)."""
    out = {}
    for path in paths(d):
        full = os.path.join(d, path)
        if path.endswith(b"/") or not os.path.isfile(full):
            out[path] = None
        else:
            with open(full, "rb") as data:
                out[path] = data.read()
    return out


def made(d: bytes, before: dict, after: dict) -> bytes:
    """Each path below `d` that the run made, changed or removed (`before` and `after` are
    `snapshot`s), one per line in byte order, a byte that is not printable ASCII or is the
    backslash as `\\xHH`; a file's line then has its size and the FNV-1a 64 of its bytes, the
    build's dates as `<DATE>` (`undated`) and `d` as `{D}`, after tabs; a removed path's line ends
    with a tab and `gone`."""
    out = bytearray()
    for path in sorted(set(before) | set(after)):
        if path in before and path in after and before[path] == after[path]:
            continue
        for b in path:
            out += bytes([b]) if 0x21 <= b <= 0x7E and b != 0x5C else b"\\x%02x" % b
        if path not in after:
            out += b"\tgone"
        elif after[path] is not None:
            data = undated(d, path, after[path]).replace(d, b"{D}")
            out += b"\t%d\t%016x" % (len(data), fnv1a64(data))
        out += b"\n"
    return bytes(out)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ncbi-bin", default=os.environ.get("NCBI_BIN", ""))
    parser.add_argument("--rpsbproc-bin", default=os.environ.get("RPSBPROC_BIN", ""))
    parser.add_argument("--makembindex-bin", default=os.environ.get("MAKEMBINDEX_BIN", ""))
    parser.add_argument("--work", default=os.path.join(os.environ.get("TMPDIR", "/tmp"), "non_utf8_names"))
    parser.add_argument("only", nargs="*")
    args = parser.parse_args()
    makembindex_bin = args.makembindex_bin or os.path.join(args.ncbi_bin, "makembindex")
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
        cwd = set_up(setup, d, makembindex_bin)
        run_env = {k: v for k, v in os.environ.items() if k not in CLEARED and not k.startswith("LOSAT_")}
        run_env.update(HOME=str(home), BLAST_USAGE_REPORT="false")
        run_env = {os.fsencode(k): os.fsencode(v) for k, v in run_env.items()}
        if env != "-":
            name, value = env.split("=", 1)
            run_env[name.encode()] = words(value, d)
        program, *rest = argv.split(" ")
        exe = args.rpsbproc_bin if (program == "rpsbproc" and args.rpsbproc_bin) else os.path.join(args.ncbi_bin, program)
        command = [os.fsencode(exe)] + [words(w, d) for w in rest]
        before = snapshot(d)
        done = subprocess.run(
            command, cwd=cwd, env=run_env, stdin=subprocess.DEVNULL, capture_output=True, timeout=60
        )
        (HERE / "expected" / f"{case}.files").write_bytes(made(d, before, snapshot(d)))
        for name, data in (("stdout", done.stdout), ("stderr", done.stderr)):
            (HERE / "expected" / f"{case}.{name}").write_bytes(data.replace(d, b"{D}"))
        (HERE / "expected" / f"{case}.rc").write_text(f"{done.returncode}\n")
        print(case, kind, "rc", done.returncode, "stdout", len(done.stdout), "stderr", len(done.stderr))
    return 0


if __name__ == "__main__":
    sys.exit(main())
