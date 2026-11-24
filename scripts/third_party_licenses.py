#!/usr/bin/env python3
"""Generate THIRD_PARTY_LICENSES, the licenses of the code that LOSAT's binaries link, and keep the
copies of the license files in step.

No NCBI counterpart: this is LOSAT's own packaging tooling.

Why: the license texts of MIT, BSD, Zlib and Apache-2.0 code ask that the notice travel with every
copy, binaries included. LOSAT's release archives, wheels and LOSAT Web (the npm package
web/packages/losat) carry the compiled code of its Cargo dependencies, of the Rust standard library
and (in the WASI builds) of wasi-libc, and none of them listed these (the crate carries source only). `THIRD_PARTY_NOTICES` (written by hand) names the code copied into LOSAT's own
source; this script writes the other half, `THIRD_PARTY_LICENSES`, from `Cargo.lock`:

- the crates of each build (BUILDS below), as Cargo's resolver links them: `cargo tree -e
  normal,no-proc-macro --target <triple> <features>`, so build-, dev- and proc-macro crates (never
  linked) stay out; the union is printed once and each crate says which builds link it;
- per crate the license files of its directory in Cargo's registry source, for one license of an
  `OR` expression (MIT first, then ISC/BSD/Zlib/..., then Apache-2.0; never a GPL option), every part
  of an `AND`; the full Apache License 2.0 once at the end;
- the Rust standard library and its out-of-tree crates (scripts/third_party_licenses/
  rust-<version>-library.txt), and wasi-libc (wasi-libc-<commit>.txt) for the WASI builds, which are
  committed inputs because Cargo does not list them. The Rust ones start from the toolchain's
  share/doc/rust/COPYRIGHT-library.html (the rustc component), which omits the registry crates that the
  standard library's `backtrace` feature adds; the bundle pins the texts with a sha256, so that the
  output does not depend on the toolchain installed where CI runs.
- what zig links into the Linux wheels of pylosat (zig-<version>.txt: LLVM's libunwind and zig's
  own license, from the ziglang wheel that python-wheels.yml pins).

A license text that the script cannot find or classify fails the run: nothing falls through.
Trailing whitespace of every line is stripped and line ends are `\\n` in the output, so the bytes do
not depend on the machine.

The script also owns the copies of the hand-written license files: `LICENSE` and
`THIRD_PARTY_NOTICES` of the repository root are the originals; `LOSAT/` and `python/` carry
byte-identical copies (a crate and a wheel hold only their own directory), and `python/` also the
copy of THIRD_PARTY_LICENSES (PEP 639 license files must lie under python/).

Usage:
  third_party_licenses.py             write the generated files and the copies
  third_party_licenses.py --check     fail (exit 1) if any of them differs from what the script
                                      would write (CI, `rust lint`); needs network for `cargo fetch`
  third_party_licenses.py --refresh-rust-inputs     rewrite the committed Rust standard library
                                      inputs from the installed toolchain (after a toolchain change)
  third_party_licenses.py --refresh-wasi-libc-input rewrite the wasi-libc input from GitHub
  third_party_licenses.py --refresh-zig-input <ziglang wheel>   rewrite the zig input from the
                                      ziglang wheel of ZIG_VERSION (any platform's: the files are the same)
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import re
import subprocess
import sys
import tempfile
import tomllib
import urllib.request
from pathlib import Path
from typing import Iterable, NamedTuple

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "scripts" / "third_party_licenses"
WORKSPACE_PACKAGES = frozenset({"losat", "pylosat", "losat-web-adapter"})
RULE = "=" * 80


class Build(NamedTuple):
    group: str  # the artifact: which line of "Linked into" names it
    label: str  # the target, as the file prints it
    kind: str  # "native" or "wasi": builds that link the same system parts
    manifest_dir: str  # the Cargo root (each has its own Cargo.lock)
    triple: str
    flags: tuple[str, ...]


CLI, WHEEL, WEB = "losat CLI", "pylosat wheel", "LOSAT Web"
GROUP_ORDER = (CLI, WHEEL, WEB)
GROUP_SIZE_WORD = {CLI: "targets", WHEEL: "targets", WEB: "reactors"}
NATIVE = (
    ("Linux x86-64", "x86_64-unknown-linux-gnu"),
    ("Linux aarch64", "aarch64-unknown-linux-gnu"),
    ("Windows x86-64", "x86_64-pc-windows-msvc"),
    ("macOS aarch64", "aarch64-apple-darwin"),
    ("macOS x86-64", "x86_64-apple-darwin"),
)
# Where each row comes from (checked against the workflows 2026-11-20):
# - losat CLI, native: .github/workflows/release-readiness-v020.yml, job native-artifacts (matrix
#   target_triple; `cd LOSAT && cargo build --release --bin losat --locked`, default features).
# - losat CLI, WASI: the same workflow, serial-wasm-artifact (`--target wasm32-wasip1
#   --no-default-features`) and threaded-wasm-artifact (`--target wasm32-wasip1-threads --features
#   wasm-threads`, default features on).
# - pylosat wheel: .github/workflows/python-wheels.yml (`maturin build --release --locked` in
#   python/, default features) on ubuntu-24.04, ubuntu-24.04-arm, windows-2025, macos-15 and
#   macos-15-intel: the five native triples.
# - LOSAT Web: web/adapter/tools/build_reactors.py (SPECS), which web.yml runs: the serial reactor
#   `--target wasm32-wasip1 --lib`, the threaded one `--target wasm32-wasip1-threads --lib
#   --features threads` (adapter feature `threads` = losat/parallel + losat/wasm-threads). The
#   adapter builds the engine with default-features = false and its feature `web-abi-v1` (ABI v1 and
#   clap, which only the Web's modules carry; line RG item K4i-3d). LOSAT/.cargo/config.toml's
#   `build-web-unknown` (wasm32-unknown-unknown) is used by no workflow, so it is not listed.
BUILDS: tuple[Build, ...] = (
    *(Build(CLI, label, "native", "LOSAT", triple, ()) for label, triple in NATIVE),
    Build(CLI, "WASI", "wasi", "LOSAT", "wasm32-wasip1", ("--no-default-features",)),
    Build(CLI, "WASI threads", "wasi", "LOSAT", "wasm32-wasip1-threads", ("--features", "wasm-threads")),
    *(Build(WHEEL, label, "native", "python", triple, ()) for label, triple in NATIVE),
    Build(WEB, "serial reactor", "wasi", "web/adapter", "wasm32-wasip1", ()),
    Build(WEB, "threaded reactor", "wasi", "web/adapter", "wasm32-wasip1-threads", ("--features", "threads")),
)
KIND_WORD = {"native": "native targets", "wasi": "WASI targets"}

# Files of the repository that are copies of another (the originals are hand-written).
COPIES = (
    ("LICENSE", "LOSAT/LICENSE"),
    ("LICENSE", "python/LICENSE"),
    ("THIRD_PARTY_NOTICES", "LOSAT/THIRD_PARTY_NOTICES"),
    ("THIRD_PARTY_NOTICES", "python/THIRD_PARTY_NOTICES"),
    ("THIRD_PARTY_LICENSES", "python/THIRD_PARTY_LICENSES"),
)
GENERATED = "THIRD_PARTY_LICENSES"

# The committed inputs: file -> (origin, sha256 of the committed file). Refreshing prints the new
# sha256; paste it here.
APACHE_FILE = "Apache-2.0.txt"
INPUT_ORIGINS = {
    APACHE_FILE: "https://raw.githubusercontent.com/rust-lang/rust/1.92.0/LICENSE-APACHE",
    "rust-1.92.0-library.txt": "the toolchain's share/doc/rust/COPYRIGHT-library.html and licenses/, "
    "the registry crates of STD_CRATES, rust-lang/rust 1.92.0 LICENSE-MIT and library/compiler-builtins/LICENSE.txt",
    "wasi-libc-3f7eb4c7.txt": "https://github.com/WebAssembly/wasi-libc at 3f7eb4c7d6ede4dde3c4bffa6ed14e8d656fe93f",
    "zig-0.17.1.txt": "the ziglang 0.17.1 wheel of PyPI (ziglang/LICENSE, ziglang/lib/libunwind/LICENSE.TXT)",
}
INPUT_SHA256 = {
    APACHE_FILE: "62c7a1e35f56406896d7aa7ca52d0cc0d272ac022b5d2796e7d6905db8a3636a",
    "rust-1.92.0-library.txt": "bf05decc03b10368e2ddbad0dfd47de1b1bd77fb909633dc113ab16c1e93faef",
    "wasi-libc-3f7eb4c7.txt": "d2e4a63c223afff5387e10e067c5eb6d592c05f422d968d1473c0450afe9d419",
    "zig-0.17.1.txt": "156e39e3ff406b5ed3654c2da5c9a3c70d150542324e5b2c07120d0000e659f0",
}
RUST_VERSION = "1.92.0"
RUST_BUNDLE = f"rust-{RUST_VERSION}-library.txt"
RUST_TAG_URL = f"https://raw.githubusercontent.com/rust-lang/rust/{RUST_VERSION}"
COMPILER_BUILTINS_LICENSE = "library/compiler-builtins/LICENSE.txt"  # the standard library links compiler_builtins
# rustc 1.92.0 ships the wasm32-wasip1(-threads) libc.a and crt1-*.o of wasi-sdk 27
# (src/ci/docker/host-x86_64/dist-various-2/Dockerfile of rust-lang/rust at 1.92.0:
# wasi-sdk-27.0-x86_64-linux.tar.gz; the toolchain's libc.a carries the string
# "wasisdk://v27.0/build/sysroot/wasi-libc-wasm32-wasip1"); wasi-sdk-27 pins this wasi-libc commit
# (its submodule src/wasi-libc).
WASI_LIBC_COMMIT = "3f7eb4c7d6ede4dde3c4bffa6ed14e8d656fe93f"
WASI_LIBC_BUNDLE = f"wasi-libc-{WASI_LIBC_COMMIT[:8]}.txt"
WASI_LIBC_RAW = f"https://raw.githubusercontent.com/WebAssembly/wasi-libc/{WASI_LIBC_COMMIT}"
# python-wheels.yml links the Linux wheels with zig (`maturin build --zig`, the pip package
# ziglang==ZIG_VERSION; tests/test_third_party_licenses.py ties the two). zig links in LLVM's
# libunwind (its lib/libunwind) in place of libgcc_s, and may link functions of its own
# compiler_rt; the glibc files it compiles (the libc_nonshared parts, such as stat64) carry glibc's
# exception for linking and need no notice.
ZIG_VERSION = "0.17.1"
ZIG_BUNDLE = f"zig-{ZIG_VERSION}.txt"
ZIG_FILES = (
    # (path in the wheel, title, license, repository, notes)
    ("ziglang/lib/libunwind/LICENSE.TXT", "LLVM libunwind (zig's lib/libunwind)", "Apache-2.0 WITH LLVM-exception",
     "https://github.com/llvm/llvm-project (libunwind/)",
     ["zig links it into the Linux wheels in place of libgcc_s (Rust's panics and backtraces call it); "
      "the file also holds the legacy University of Illinois/NCSA and MIT licenses of the code "
      "contributed before 2019."]),
    ("ziglang/LICENSE", "zig (compiler_rt and the other code zig links)", "MIT", "https://ziglang.org", []),
)

# The out-of-tree crates of the standard library that a binary can link, from the `dependencies`
# of library/std/Cargo.toml and library/Cargo.lock at 1.92.0 (the toolchain's COPYRIGHT-library.html
# lists only some of them: the optional dependencies of the `backtrace` feature, which the shipped
# standard library enables, are missing there). The third column is the cfg of the dependency:
# "all", "not-windows" (`cfg(not(all(windows, target_env = "msvc")))`, which is also the cfg of
# the backtrace crates) and "wasi" (`cfg(target_os = "wasi")`). Crates for other systems (hermit, sgx,
# uefi, vexos, motor, xous, wasip2) and the test and build dependencies are not built for our targets.
STD_CRATES = (
    ("addr2line", "0.25.0", "not-windows"),
    ("adler2", "2.0.1", "not-windows"),
    ("cfg-if", "1.0.1", "all"),
    ("gimli", "0.32.0", "not-windows"),
    ("hashbrown", "0.15.5", "all"),
    ("libc", "0.2.177", "not-windows"),
    ("memchr", "2.7.5", "not-windows"),
    ("miniz_oxide", "0.8.9", "not-windows"),
    ("object", "0.37.3", "not-windows"),
    ("rustc-demangle", "0.1.26", "all"),
    ("wasi", "0.11.1+wasi-snapshot-preview1", "wasi"),
)
STD_IN_TREE = (
    # (path in rust-lang/rust, license, copyright lines), from COPYRIGHT-library.html: the in-tree
    # exceptions that land in binaries (the Fuchsia mutex does not: no Fuchsia target).
    ("library/backtrace", "Apache-2.0 OR MIT", ("2014 Alex Crichton", "The Rust Project Developers (see https://thanks.rust-lang.org)")),
    ("library/std/src/sync/mpmc", "Apache-2.0 OR MIT", ("2019 The Crossbeam Project Developers", "The Rust Project Developers (see https://thanks.rust-lang.org)")),
    ("library/core/src/unicode/unicode_data.rs", "Unicode-3.0", ("1991-2024 Unicode, Inc",)),
)

# Preference among the alternatives of an `OR` (lower is taken): MIT where offered, else the next
# most permissive (the public-domain dedications, which ask for nothing, then the licenses that ask
# for a notice, Apache-2.0 with the LLVM exception, which lifts its notice conditions for object
# code, before plain Apache-2.0). A license not listed is never taken (a GPL option is never the
# choice) and fails the run when nothing else is left.
RANK = {
    "MIT": 0, "MIT-0": 0,
    "0BSD": 1, "Unlicense": 1, "CC0-1.0": 1,
    "ISC": 2, "BSD-2-Clause": 2, "BSD-3-Clause": 2, "Zlib": 2,
    "Unicode-3.0": 3, "bzip2-1.0.6": 3, "BSL-1.0": 3,
    "Apache-2.0 WITH LLVM-exception": 4,
    "Apache-2.0": 5,
}
APACHE = "Apache-2.0"
INFINITY = 99

# Text that identifies a license (after whitespace is collapsed).
MARKERS = {
    "MIT": r"Permission is hereby granted, free of charge, to any person obtaining a copy",
    "MIT-0": r"MIT No Attribution",
    "ISC": r"Permission to use, copy, modify, and/or distribute this software for any purpose with or without fee is hereby granted",
    "BSD-2-Clause": r"Redistribution and use in source and binary forms(?!.*Neither the name)(?!.*The origin of this software)",
    "BSD-3-Clause": r"Redistribution and use in source and binary forms.*Neither the name",
    "Zlib": r"Permission is granted to anyone to use this software for any purpose, including commercial applications",
    "bzip2-1.0.6": r"Redistribution and use in source and binary forms.*The origin of this software must not be misrepresented",
    "Unicode-3.0": r"UNICODE LICENSE V3",
    "Unlicense": r"This is free and unencumbered software released into the public domain",
    "CC0-1.0": r"CC0 1\.0 Universal",
    "BSL-1.0": r"Boost Software License - Version 1\.0",
    "0BSD": r"Permission to use, copy, modify, and/or distribute this software for any purpose with or without fee is hereby granted\.",
    "Apache-2.0": r"Apache License.*TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION",
    "Apache-2.0 WITH LLVM-exception": r"LLVM Exceptions to the Apache 2\.0 License",
}

# Per crate (only where the generic rules need help). The reason is part of the entry's note.
NOTES = {
    "linux-raw-sys": "Generated from the Linux kernel's user-space API headers (GPL-2.0 WITH "
    "Linux-syscall-note, which allows their use by user-space programs); the crate's license is the "
    "expression above.",
    "self_cell": "Used under Apache-2.0; the GPL-2.0-only alternative is not taken.",
}
# The numbers of the generated file's sections that other text refers to.
NOTICE_SECTION = 5
APACHE_SECTION = 6
# Crates whose license file holds a second license after the chosen one: print up to the marker.
CUT_BEFORE = {
    "chrono": "[2]: <http://www.apache.org/licenses/LICENSE-2.0>",  # MIT block, then the Apache text
}
# Crates whose license file holds the full Apache License 2.0 next to other text: the standard text
# (printed once, in section APACHE_SECTION) is cut out and a line says so; what comes before and after stays.
CUT_APACHE_BODY = frozenset({"libm"})
APACHE_BODY_CUT_LINE = f"[The Apache License 2.0 text of this file is printed once, in section {APACHE_SECTION}.]"
# Crates that ship no license file for their license: the standard text of the license is the one
# at the end of this file. name -> what to say.
NO_LICENSE_FILE = {
    "flatbuffers": "The crate ships no license file; its license field is Apache-2.0, whose text is "
    "at the end of this file.",
}


class LicenseError(Exception):
    """A license that the script does not know, or a text it cannot find."""


# ---------------------------------------------------------------------------------------------
# SPDX expressions and the choice of one license
# ---------------------------------------------------------------------------------------------

def parse_spdx(expr: str):
    """('or'|'and', [nodes]) or ('id', 'MIT' / 'Apache-2.0 WITH LLVM-exception'); the legacy `/`
    of old manifests (`MIT/Apache-2.0`) is OR."""
    tokens = re.findall(r"\(|\)|[^\s()/]+|/", expr)
    position = 0

    def peek() -> str | None:
        return tokens[position] if position < len(tokens) else None

    def take() -> str:
        nonlocal position
        if position >= len(tokens):
            raise LicenseError(f"license expression {expr!r} ends early")
        position += 1
        return tokens[position - 1]

    def is_word(token: str | None, word: str) -> bool:
        return token is not None and token.upper() == word

    def parse_or():
        nodes = [parse_and()]
        while is_word(peek(), "OR") or peek() == "/":
            take()
            nodes.append(parse_and())
        return nodes[0] if len(nodes) == 1 else ("or", nodes)

    def parse_and():
        nodes = [parse_term()]
        while is_word(peek(), "AND"):
            take()
            nodes.append(parse_term())
        return nodes[0] if len(nodes) == 1 else ("and", nodes)

    def parse_term():
        token = take()
        if token == "(":
            node = parse_or()
            if take() != ")":
                raise LicenseError(f"license expression {expr!r}: ')' expected")
            return node
        if token in (")", "/") or token.upper() in ("AND", "OR", "WITH"):
            raise LicenseError(f"license expression {expr!r}: unexpected {token!r}")
        if is_word(peek(), "WITH"):
            take()
            return ("id", f"{token} WITH {take()}")
        return ("id", token)

    node = parse_or()
    if position != len(tokens):
        raise LicenseError(f"license expression {expr!r}: unexpected {tokens[position]!r}")
    return node


def cost(node) -> int:
    kind, value = node
    if kind == "id":
        return RANK.get(value, INFINITY)
    costs = [cost(child) for child in value]
    return max(costs) if kind == "and" else min(costs)


def choose(expr: str) -> list[str]:
    """The licenses that apply when the expression is used under its most preferred alternative."""

    def pick(node) -> list[str]:
        kind, value = node
        if kind == "id":
            if value not in RANK:
                raise LicenseError(f"license {value!r} (in {expr!r}) is not one the script knows")
            return [value]
        if kind == "and":
            picked: list[str] = []
            for child in value:
                picked += [name for name in pick(child) if name not in picked]
            return picked
        best = min(value, key=cost)  # the first of equals
        if cost(best) >= INFINITY:
            raise LicenseError(f"no license of {expr!r} is one the script knows")
        return pick(best)

    return pick(parse_spdx(expr))


def has_choice(expr: str) -> bool:
    def walk(node) -> bool:
        kind, value = node
        return kind == "or" or (kind == "and" and any(walk(child) for child in value))

    return walk(parse_spdx(expr))


# ---------------------------------------------------------------------------------------------
# License files
# ---------------------------------------------------------------------------------------------

def normalize(text: str) -> str:
    """LF line ends, no trailing whitespace, one final newline (so that the output's bytes do not
    depend on the machine or on an editor)."""
    lines = [line.rstrip() for line in text.replace("\r\n", "\n").replace("\r", "\n").split("\n")]
    while lines and not lines[-1]:
        lines.pop()
    return "\n".join(lines) + "\n"


def drop_leading_blank_lines(text: str) -> str:
    """The texts that the Rust inputs took from COPYRIGHT-library.html start with an empty line,
    which would print as two blank lines under the heading."""
    return text.lstrip("\n")


def collapse(text: str) -> str:
    return " ".join(text.split())


def classify(text: str) -> set[str]:
    flat = collapse(text)
    return {name for name, pattern in MARKERS.items() if re.search(pattern, flat, re.S)}


FILE_NAME_RE = re.compile(r"(?i)^(licen[sc]e|copying|unlicense|copyright|notice)")


def role(file_name: str) -> str:
    lowered = file_name.lower()
    if lowered.startswith("notice"):
        return "notice"
    if lowered.startswith("copyright"):
        return "copyright"
    return "license"


def read_license_files(directory: Path, license_file: str | None = None) -> dict[str, str]:
    """The license-like files of a crate directory (and the manifest's `license-file`), normalized."""
    found = {}
    names = [p.name for p in directory.iterdir() if p.is_file() and FILE_NAME_RE.match(p.name)]
    if license_file and (directory / license_file).is_file() and license_file not in names:
        names.append(license_file)
    for name in sorted(names):
        try:
            found[name] = normalize((directory / name).read_bytes().decode("utf-8"))
        except UnicodeDecodeError as error:
            raise LicenseError(f"{directory / name} is not UTF-8: {error}")
    return found


class Texts(NamedTuple):
    licenses: list[tuple[str, str]]  # (file, text) to print
    apache: list[tuple[str, str]]  # the crate's own Apache-2.0 files (the standard text is printed once)
    notices: list[tuple[str, str]]  # NOTICE files
    copyrights: list[tuple[str, str]]  # COPYRIGHT files


def select_texts(name: str, chosen: list[str], files: dict[str, str]) -> Texts:
    """The files that carry the chosen licenses; LicenseError when one is missing."""
    printed: list[tuple[str, str]] = []
    apache: list[tuple[str, str]] = []
    for wanted in chosen:
        matches = [f for f in sorted(files) if role(f) == "license" and wanted in classify(files[f])]
        if not matches:
            if wanted == APACHE and name in NO_LICENSE_FILE:
                continue
            raise LicenseError(f"{name}: no file carries the text of {wanted} (files: {sorted(files)}); "
                               "add an explicit entry to the script if the crate ships none")
        for file_name in matches:
            target = apache if wanted == APACHE else printed
            known = [(f, t) for f, t in printed + apache]
            if all(file_name != existing and collapse(files[file_name]) != collapse(t) for existing, t in known):
                text = files[file_name]
                if wanted != APACHE and name in CUT_BEFORE:
                    marker = CUT_BEFORE[name]
                    if marker not in text:
                        raise LicenseError(f"{name}: {file_name} has no {marker!r} any more")
                    text = normalize(text[: text.index(marker)])
                if wanted != APACHE and name in CUT_APACHE_BODY:
                    text = cut_apache_body(name, file_name, text)
                target.append((file_name, text))
    notices = [(f, files[f]) for f in sorted(files) if role(f) == "notice"]
    copyrights = [(f, files[f]) for f in sorted(files) if role(f) == "copyright"]
    return Texts(printed, apache, notices, copyrights)


def cut_apache_body(name: str, file_name: str, text: str) -> str:
    """The text without the full Apache License 2.0 that it holds (from its title line to the last
    line of the appendix, "limitations under the License."), a line saying where the text is
    printed in its place. The cut part must be the standard text; otherwise the run fails."""
    lines = text.split("\n")
    starts = [n for n, line in enumerate(lines) if line.strip() == "Apache License"]
    if len(starts) != 1:
        raise LicenseError(f"{name}: {file_name} has {len(starts)} lines 'Apache License', expected one")
    ends = [n for n, line in enumerate(lines) if n > starts[0] and line.strip().endswith("limitations under the License.")]
    if len(ends) != 1:
        raise LicenseError(f"{name}: {file_name} has {len(ends)} lines 'limitations under the License.' after the Apache text, expected one")
    cut = "\n".join(lines[starts[0]: ends[0] + 1])
    flat, flat_standard = collapse(cut), collapse(apache_standard_text())
    # The standard text up to its appendix (a crate's own appendix may differ in indentation only).
    if flat.split("END OF TERMS AND CONDITIONS")[0] != flat_standard.split("END OF TERMS AND CONDITIONS")[0]:
        raise LicenseError(f"{name}: the Apache text in {file_name} differs from the standard Apache-2.0 text; read it "
                           "and decide before the script leaves it out")
    return normalize("\n".join(lines[: starts[0]] + [APACHE_BODY_CUT_LINE] + lines[ends[0] + 1:]))


def apache_extras(name: str, file_name: str, text: str, standard: str) -> list[str]:
    """Lines of a crate's Apache-2.0 file beyond the standard text: its own copyright lines (a
    filled-in appendix) and anything after the appendix, such as a list of file licenses. The part up to END OF TERMS AND CONDITIONS must be the
    standard text; otherwise the crate's file is a different text and the run fails."""
    end = "END OF TERMS AND CONDITIONS"
    flat, flat_standard = collapse(text), collapse(standard)
    if end not in flat:
        raise LicenseError(f"{name}: {file_name} is not the Apache-2.0 text")
    head = flat[: flat.index(end) + len(end)]
    if head != flat_standard[: flat_standard.index(end) + len(end)]:
        raise LicenseError(f"{name}: {file_name} differs from the standard Apache-2.0 text; read it and "
                           "decide before the script prints it as the standard one")
    tail = [line.rstrip() for line in text[text.index(end):].splitlines()[1:]]
    appendix_end = next((n for n, line in enumerate(tail) if "limitations under the License." in line), None)
    if appendix_end is None:
        return [line for line in tail if line.strip()]
    # The appendix is the how-to-apply boilerplate; a filled-in copyright line is the crate's own.
    own = [line for line in tail[:appendix_end + 1]
           if line.strip().startswith("Copyright") and "[yyyy] [name of copyright owner]" not in line]
    return own + [line for line in tail[appendix_end + 1:] if line.strip()]


# ---------------------------------------------------------------------------------------------
# Committed inputs
# ---------------------------------------------------------------------------------------------

def read_bundle(path: Path) -> list[tuple[dict, str]]:
    """A bundle is blocks: a line `@@@ {json}` and the text up to the next such line."""
    blocks: list[tuple[dict, list[str]]] = []
    for line in path.read_text(encoding="utf-8").split("\n"):
        if line.startswith("@@@ "):
            blocks.append((json.loads(line[4:]), []))
        elif blocks:
            blocks[-1][1].append(line)
        elif line:
            raise LicenseError(f"{path}: text before the first block")
    return [(header, drop_leading_blank_lines(normalize("\n".join(lines))) if "\n".join(lines).strip() else "")
            for header, lines in blocks]


def write_bundle(path: Path, blocks: Iterable[tuple[dict, str]]) -> None:
    out = []
    for header, text in blocks:
        out.append("@@@ " + json.dumps(header, ensure_ascii=False, sort_keys=True) + "\n" + text)
    path.write_text("".join(out), encoding="utf-8")


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def toolchain_channel() -> str:
    return tomllib.loads((ROOT / "rust-toolchain.toml").read_text(encoding="utf-8"))["toolchain"]["channel"]


def input_problems() -> list[str]:
    """Why the committed inputs cannot be used as they are (empty: they can)."""
    problems = []
    for name, expected in INPUT_SHA256.items():
        path = INPUTS / name
        if not path.is_file():
            problems.append(f"{path.relative_to(ROOT)} is missing")
        elif sha256_of(path) != expected:
            problems.append(f"{path.relative_to(ROOT)} has sha256 {sha256_of(path)}, INPUT_SHA256 of the script has "
                            f"{expected or '(none)'}: refresh the input or update the table")
    if toolchain_channel() != RUST_VERSION:
        problems.append(f"rust-toolchain.toml has channel {toolchain_channel()}, the Rust inputs are of {RUST_VERSION}: "
                        "update RUST_VERSION, STD_CRATES, WASI_LIBC_COMMIT and the INPUT_ tables of "
                        "scripts/third_party_licenses.py, then run it with --refresh-rust-inputs and "
                        "--refresh-wasi-libc-input")
    return problems


# ---------------------------------------------------------------------------------------------
# Cargo
# ---------------------------------------------------------------------------------------------

def cargo(args: list[str], cwd: Path) -> str:
    env = dict(os.environ, CARGO_TERM_COLOR="never")
    done = subprocess.run(["cargo", *args], cwd=cwd, capture_output=True, text=True, env=env)
    if done.returncode != 0:
        raise LicenseError(f"cargo {' '.join(args)} (in {cwd}) failed:\n{done.stderr}")
    return done.stdout


TREE_LINE = re.compile(r"^(?P<name>\S+) v(?P<version>\S+)(?: \((?P<path>[^)]*)\))?$")


def parse_tree(output: str) -> set[tuple[str, str]]:
    crates = set()
    for line in output.splitlines():
        line = re.sub(r" \(\*\)$", "", line.strip())
        if not line:
            continue
        match = TREE_LINE.match(line)
        if match is None:
            raise LicenseError(f"cannot read a line of `cargo tree`: {line!r}")
        if match["path"] is not None:
            if match["name"] not in WORKSPACE_PACKAGES:
                raise LicenseError(f"`cargo tree` lists the path package {match['name']}: add it to "
                                   "WORKSPACE_PACKAGES if it is LOSAT's own, else decide its license")
            continue
        crates.add((match["name"], match["version"]))
    return crates


def linked_crates() -> tuple[dict[tuple[str, str], set[Build]], dict[tuple[str, str], dict]]:
    """(crate -> builds that link it, crate -> its cargo metadata)."""
    metadata: dict[tuple[str, str], dict] = {}
    for directory in sorted({b.manifest_dir for b in BUILDS}):
        cargo(["fetch", "--locked"], ROOT / directory)
        data = json.loads(cargo(["metadata", "--format-version", "1", "--locked"], ROOT / directory))
        for package in data["packages"]:
            metadata.setdefault((package["name"], package["version"]), package)
    linked: dict[tuple[str, str], set[Build]] = {}
    for build in BUILDS:
        out = cargo(["tree", "-e", "normal,no-proc-macro", "--target", build.triple, *build.flags,
                     "--prefix", "none", "--format", "{p}", "--locked"], ROOT / build.manifest_dir)
        for crate in parse_tree(out):
            linked.setdefault(crate, set()).add(build)
    return linked, metadata


# ---------------------------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------------------------

def linked_into(builds: Iterable[Build]) -> str:
    present = set(builds)
    parts = []
    for group in GROUP_ORDER:
        members = [b for b in BUILDS if b.group == group]
        have = [b for b in members if b in present]
        if not have:
            continue
        word = GROUP_SIZE_WORD[group]
        if len(have) == len(members):
            parts.append(f"{group} (all {len(members)} {word})" if len(members) > 2 else f"{group} (both {word})")
            continue
        kinds = {b.kind for b in have}
        if len(kinds) == 1 and group != WEB:
            kind = kinds.pop()
            if len(have) == len([b for b in members if b.kind == kind]):
                parts.append(f"{group} ({KIND_WORD[kind]}: {', '.join(b.label for b in have)})")
                continue
        parts.append(f"{group} ({', '.join(b.label for b in have)})")
    return "; ".join(parts)


def std_builds(applies: str) -> list[Build]:
    if applies == "all":
        return list(BUILDS)
    if applies == "not-windows":
        return [b for b in BUILDS if "windows" not in b.triple]
    if applies == "wasi":
        return [b for b in BUILDS if b.kind == "wasi"]
    raise LicenseError(f"unknown applicability {applies!r}")


def wrap(lines: list[str], text: str, file_name: str | None) -> None:
    lines.append("")
    if file_name:
        lines.append(f"--- {file_name} ---")
        lines.append("")
    lines.extend(text.rstrip("\n").split("\n"))


def entry(title: str, expr: str | None, *, repository: str | None = None, linked: str | None = None,
          notes: Iterable[str] = (), texts: Texts | None = None, apache_standard: str = "",
          notice_refs: dict[str, str] | None = None, name: str = "", extra_texts: Iterable[tuple[str, str]] = ()) -> list[str]:
    lines = [RULE, title]
    chosen: list[str] = []
    if expr:
        chosen = choose(expr)
        under = f" (used under {' AND '.join(chosen)})" if has_choice(expr) else ""
        lines.append(f"License: {expr}{under}")
    if repository:
        lines.append(f"Repository: {repository}")
    if linked:
        lines.append(f"Linked into: {linked}")
    for note in notes:
        lines.append(f"Note: {note}")
    if texts is not None:
        if texts.apache or (APACHE in chosen and name in NO_LICENSE_FILE):
            lines.append("Apache License 2.0: the full text is at the end of this file.")
        if name in NO_LICENSE_FILE and APACHE in chosen:
            lines.append(f"Note: {NO_LICENSE_FILE[name]}")
        for file_name, text in texts.apache:
            extras = apache_extras(name, file_name, text, apache_standard)
            if extras:
                wrap(lines, "\n".join(extras) + "\n", f"{file_name}: lines beyond the standard Apache License 2.0 text")
        for file_name, text in texts.licenses:
            wrap(lines, text, file_name)
        for file_name, text in texts.copyrights:
            wrap(lines, text, file_name)
        for file_name, text in texts.notices:
            lines.append(f"NOTICE: {(notice_refs or {})[collapse(text)]}")
    for file_name, text in extra_texts:
        wrap(lines, text, file_name)
    lines.append("")
    return lines


def head(builds_text: list[str]) -> list[str]:
    return [
        "THIRD_PARTY_LICENSES",
        "",
        "The licenses of the code that LOSAT's binaries link: the `losat` command, the `pylosat`",
        "wheel and LOSAT Web (its two Wasm reactors). LOSAT's own work is under the MIT License",
        "(LICENSE); the code written by NCBI is in the public domain in the United States",
        "(THIRD_PARTY_NOTICES).",
        "",
        "This file is generated by scripts/third_party_licenses.py from Cargo.lock (LOSAT/, python/",
        "and web/adapter/) and the license files of the crates; do not edit it by hand. Code that",
        "was copied into LOSAT's own source (LMDB, Boost, Arm, Sun, Microsoft's UCRT, bzip2, ruzstd,",
        "BitMagic) is in THIRD_PARTY_NOTICES.",
        "",
        "Where a crate offers its code under several licenses (an OR expression), it is used under",
        "the one named after \"used under\": MIT where offered, else the next most permissive; a",
        "GPL option is never taken. Where a license asks for its text to be given with copies, the",
        "text is printed below as the crate ships it; the Apache License 2.0 is printed once, at",
        "the end, and the entries that use it refer to it.",
        "",
        "Builds covered (the crates below are the union of what Cargo links into them):",
        *builds_text,
        "",
        "Contents:",
        "  1. Crates (Cargo dependencies)",
        "  2. The Rust standard library, linked into every build",
        "  3. wasi-libc, linked into the WASI builds and LOSAT Web",
        "  4. What zig links into the Linux wheels of pylosat",
        f"  {NOTICE_SECTION}. NOTICE files of Apache-2.0 crates",
        f"  {APACHE_SECTION}. Apache License, Version 2.0",
        "",
    ]


def builds_lines() -> list[str]:
    lines = []
    for build in BUILDS:
        flags = " ".join(build.flags) if build.flags else "default features"
        lines.append(f"  - {build.group}, {build.label}: {build.triple}, {flags}")
    return lines


def render(linked: dict[tuple[str, str], set[Build]], metadata: dict[tuple[str, str], dict],
           standard_apache: str, rust_blocks: list[tuple[dict, str]], wasi_blocks: list[tuple[dict, str]],
           zig_blocks: list[tuple[dict, str]]) -> str:
    # Gather the texts of every crate first: the NOTICE files are printed once, by text.
    crates = []
    notice_groups: dict[str, dict] = {}
    for key in sorted(linked):
        package = metadata.get(key)
        if package is None:
            raise LicenseError(f"{key[0]} {key[1]} is linked but not in `cargo metadata`")
        expr = package.get("license")
        if not expr:
            raise LicenseError(f"{key[0]} {key[1]} has no `license` field (license-file only): decide its license")
        files = read_license_files(Path(package["manifest_path"]).parent, package.get("license_file"))
        texts = select_texts(key[0], choose(expr), files)
        for _, text in texts.notices:
            group = notice_groups.setdefault(collapse(text), {"text": text, "crates": []})
            group["crates"].append(f"{key[0]} {key[1]}")
        crates.append((key, package, expr, texts))
    notice_numbers = {flat: n for n, flat in enumerate(sorted(notice_groups), 1)}
    notice_refs = {flat: f"see section {NOTICE_SECTION}, NOTICE {notice_numbers[flat]} (shared by {len(group['crates'])} crate"
                         f"{'s' if len(group['crates']) != 1 else ''})" for flat, group in notice_groups.items()}

    out = head(builds_lines())
    out += ["", "1. CRATES (Cargo dependencies)", ""]
    for key, package, expr, texts in crates:
        out += entry(f"{key[0]} {key[1]}", expr, repository=package.get("repository"),
                     linked=linked_into(linked[key]), notes=[NOTES[key[0]]] if key[0] in NOTES else [],
                     texts=texts, apache_standard=standard_apache, notice_refs=notice_refs, name=key[0])

    out += ["", "2. THE RUST STANDARD LIBRARY", ""]
    out += render_rust(rust_blocks)
    out += ["", "3. WASI-LIBC", ""]
    out += render_wasi(wasi_blocks)
    out += ["", "4. WHAT ZIG LINKS INTO THE LINUX WHEELS OF PYLOSAT", ""]
    out += render_zig(zig_blocks)

    out += ["", f"{NOTICE_SECTION}. NOTICE FILES OF APACHE-2.0 CRATES (Apache License 2.0, section 4(d))", ""]
    for flat, number in sorted(notice_numbers.items(), key=lambda item: item[1]):
        group = notice_groups[flat]
        out += [RULE, f"NOTICE {number}: the NOTICE file of {', '.join(group['crates'])}"]
        wrap(out, group["text"], None)
        out.append("")
    if not notice_groups:
        out += ["(none)", ""]

    out += ["", f"{APACHE_SECTION}. APACHE LICENSE, VERSION 2.0", "", RULE, "Apache License, Version 2.0 (the text the entries above refer to)"]
    wrap(out, standard_apache, None)
    out.append("")
    return "\n".join(out)


def render_rust(blocks: list[tuple[dict, str]]) -> list[str]:
    meta = next(h for h, _ in blocks if h["kind"] == "meta")
    own = next((h, t) for h, t in blocks if h["kind"] == "rust-license")
    builtins = next((h, t) for h, t in blocks if h["kind"] == "compiler-builtins-license")
    lines = [
        f"The standard library of rustc {meta['channel']} (the `std`, `core` and `alloc` crates) is",
        "statically linked into every Rust binary; Cargo does not list it. Its notice file is the",
        "toolchain's share/doc/rust/COPYRIGHT-library.html; the texts below were taken from it (and,",
        "for the five crates that it leaves out, from the crates' registry sources).",
        "",
    ]
    in_tree = [h for h, _ in blocks if h["kind"] == "in-tree"]
    mit_choice = "Apache-2.0 OR MIT"
    lines += entry("The Rust standard library (library/ of rust-lang/rust)", mit_choice,
                   repository="https://github.com/rust-lang/rust", linked=linked_into(BUILDS),
                   notes=["Copyright (c) The Rust Project Contributors. Parts with their own copyright "
                          "(Apache-2.0 OR MIT, used under MIT, with the text below): "
                          + "; ".join(f"{h['path']} (Copyright {' / '.join(h['copyright'])})" for h in in_tree if h["license"] == mit_choice) + "."],
                   extra_texts=[("LICENSE-MIT", own[1])])
    unicode_blocks = [(h, t) for h, t in blocks if h["kind"] == "license-text" and h["license"] == "Unicode-3.0"]
    for header in (h for h in in_tree if h["license"] == "Unicode-3.0"):
        lines += entry(f"The Rust standard library: {header['path']}", "Unicode-3.0",
                       linked=linked_into(BUILDS), notes=[f"Copyright {' / '.join(header['copyright'])}."],
                       extra_texts=[("Unicode-3.0", unicode_blocks[0][1])])
    applies = {name: (version, a) for name, version, a in STD_CRATES}
    by_crate: dict[str, list[tuple[dict, str]]] = {}
    for header, text in blocks:
        if header["kind"] == "crate":
            by_crate.setdefault(header["name"], []).append((header, text))
    for name in sorted(applies):
        version, where = applies[name]
        if name not in by_crate:
            raise LicenseError(f"{RUST_BUNDLE} has no text of the standard library crate {name}: refresh it")
        expr = by_crate[name][0][0]["license"]
        files = {h["file"]: t for h, t in by_crate[name]}
        texts = select_texts(name, choose(expr), files)
        if texts.apache:
            raise LicenseError(f"{name}: the standard library crate is used under Apache-2.0, which this script does not print for these")
        lines += entry(f"{name} {version} (Rust standard library dependency)", expr, linked=linked_into(std_builds(where)),
                       texts=texts, name=name)
    lines += [RULE, "compiler_builtins (compiler-rt routines)",
              "License: MIT AND Apache-2.0 WITH LLVM-exception AND (MIT OR Apache-2.0)",
              "Linked into: every build",
              "Note: the standard library links the compiler_builtins crate, Rust's port of LLVM compiler-rt's "
              "builtins and of libm routines (authors: Jorge Aparicio and contributors). Its license file "
              f"({builtins[0]['origin'].removeprefix(RUST_TAG_URL + '/')} of rust {RUST_VERSION}) is printed below, but for "
              f"the text of the Apache License 2.0, which is printed once, in section {APACHE_SECTION}. The LLVM exception of "
              "the Apache License 2.0 lifts sections 4(a), 4(b) and 4(d) for portions embedded in object code."]
    wrap(lines, cut_apache_body("compiler_builtins", builtins[0]["file"], builtins[1]), builtins[0]["file"])
    lines.append("")
    return lines


def render_wasi(blocks: list[tuple[dict, str]]) -> list[str]:
    meta = next(h for h, _ in blocks if h["kind"] == "meta")
    lines = [
        f"rustc {RUST_VERSION} ships wasm32-wasip1 and wasm32-wasip1-threads with the libc.a, crt1-command.o",
        "and crt1-reactor.o of wasi-sdk 27, which holds wasi-libc at commit",
        f"{meta['commit']} (https://github.com/WebAssembly/wasi-libc). A WASI build links the parts of",
        "it that the program uses. Its LICENSE says which license covers which directory; the texts below",
        "are the ones it names, as committed in that revision. A given build links only the objects",
        "it needs.",
        "",
    ]
    linked = linked_into([b for b in BUILDS if b.kind == "wasi"])
    groups: dict[str, list[tuple[dict, str]]] = {}
    for header, text in blocks:
        if header["kind"] == "text":
            groups.setdefault(header["title"], []).append((header, text))
    for title, members in groups.items():
        notes = [note for header, _ in members for note in header.get("notes", [])]
        lines += entry(title, members[0][0]["license"], repository=members[0][0]["directory"], linked=linked,
                       notes=notes, extra_texts=[(header["file"], text) for header, text in members])
    lines += [RULE, "dlmalloc (wasi-libc's malloc)", "License: CC0-1.0",
              f"Linked into: {linked}",
              "Note: Doug Lea's dlmalloc is released to the public domain (CC0-1.0, "
              "http://creativecommons.org/publicdomain/zero/1.0/, as stated in dlmalloc/src/malloc.c); "
              "no notice is required.", ""]
    return lines


def render_zig(blocks: list[tuple[dict, str]]) -> list[str]:
    meta = next(h for h, _ in blocks if h["kind"] == "meta")
    lines = [
        f"The Linux wheels of pylosat are linked by zig {meta['version']} (`maturin build --zig`), so that they",
        "run on glibc 2.17 and later. zig links in LLVM's libunwind, and can link functions of its own",
        "compiler_rt. It also compiles glibc's small static parts (stat64, fstat64, lstat64, atexit and",
        "the like, which a shared libc leaves to libc_nonshared.a); their files carry glibc's exception that",
        "permits linking them into any program and distributing it without restriction, so they need no",
        "notice.",
        "",
    ]
    linked = linked_into([b for b in BUILDS if b.group == WHEEL and "linux" in b.triple])
    for header, text in blocks:
        if header["kind"] == "text":
            lines += entry(header["title"], header["license"], repository=header["repository"],
                           linked=linked, notes=header.get("notes", []),
                           extra_texts=[(header["file"], cut_apache_body(header["title"], header["file"], text)
                                         if header["license"].startswith("Apache-2.0") else text)])
    return lines


def apache_standard_text() -> str:
    return normalize((INPUTS / APACHE_FILE).read_text(encoding="utf-8"))


def generate() -> dict[str, bytes]:
    problems = input_problems()
    if problems:
        raise LicenseError("\n".join(problems))
    linked, metadata = linked_crates()
    text = render(linked, metadata, apache_standard_text(), read_bundle(INPUTS / RUST_BUNDLE),
                  read_bundle(INPUTS / WASI_LIBC_BUNDLE), read_bundle(INPUTS / ZIG_BUNDLE))
    files = {GENERATED: text.encode("utf-8")}
    return with_copies(files)


def with_copies(files: dict[str, bytes]) -> dict[str, bytes]:
    """Add the copies of the originals (the originals themselves are read from the repository)."""
    result = dict(files)
    for original, copy in COPIES:
        if original in result:
            source = result[original]
        else:
            path = ROOT / original
            if not path.is_file():
                raise LicenseError(f"{original} is missing (the original of {copy})")
            source = path.read_bytes()
        result[copy] = source
    return result


# A pointer into THIRD_PARTY_NOTICES gives the number and the name of the section, as in
# "THIRD_PARTY_NOTICES, section 2 (LMDB)"; the name must be in the heading of that number, so that a
# new or moved section cannot leave a pointer at the wrong one.
POINTER = re.compile(r"THIRD_PARTY_NOTICES(?:,| \()\s*section\s+(\d+)(?: \(([^)]*)\))?")
POINTER_SUFFIXES = (".rs", ".toml", ".py", ".md", ".yml", ".yaml", ".sh", ".ts", ".js", ".txt", ".vue", ".html", ".css")


def notice_headings(text: str) -> dict[int, str]:
    """The sections of THIRD_PARTY_NOTICES: "N. title" lines right under a line of "="."""
    lines = text.split("\n")
    headings = {}
    for above, line in zip(lines, lines[1:]):
        match = re.fullmatch(r"(\d+)\. (.+)", line)
        if match and above and set(above) == {"="}:
            headings[int(match.group(1))] = match.group(2)
    return headings


def pointer_problems(root: Path, files: Iterable[str] | None = None) -> list[str]:
    headings = notice_headings((root / "THIRD_PARTY_NOTICES").read_text(encoding="utf-8"))
    if files is None:
        files = subprocess.run(["git", "ls-files"], cwd=root, check=True, capture_output=True, text=True).stdout.split("\n")
    problems = []
    for relative in files:
        if not relative.endswith(POINTER_SUFFIXES) or relative.startswith("docs/evidence/"):
            continue
        try:
            text = (root / relative).read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        for number, line in enumerate(text.split("\n"), 1):
            for match in POINTER.finditer(line):
                section, name = int(match.group(1)), match.group(2)
                if name is None:
                    problems.append(f"{relative}:{number}: name the section, as in "
                                    f"\"THIRD_PARTY_NOTICES, section {section} (<a word of its heading>)\"")
                elif section not in headings or name not in headings[section]:
                    problems.append(f"{relative}:{number}: THIRD_PARTY_NOTICES has no section {section} named {name!r} "
                                    f"(heading {section}: {headings.get(section, 'none')!r})")
    return problems


def first_difference(expected: bytes, actual: bytes) -> str:
    expected_lines = expected.decode("utf-8", "replace").split("\n")
    actual_lines = actual.decode("utf-8", "replace").split("\n")
    for number in range(max(len(expected_lines), len(actual_lines))):
        want = expected_lines[number] if number < len(expected_lines) else "(end of file)"
        have = actual_lines[number] if number < len(actual_lines) else "(end of file)"
        if want != have:
            return f"line {number + 1}: expected {want[:100]!r}, found {have[:100]!r}"
    return "the files differ"


def stale(expected: dict[str, bytes], root: Path) -> list[str]:
    problems = []
    for relative in sorted(expected):
        path = root / relative
        if not path.is_file():
            problems.append(f"{relative}: missing")
        elif path.read_bytes() != expected[relative]:
            problems.append(f"{relative}: stale, first difference at {first_difference(expected[relative], path.read_bytes())}")
    return problems


# ---------------------------------------------------------------------------------------------
# Refreshing the committed inputs
# ---------------------------------------------------------------------------------------------

def fetch(url: str) -> bytes:
    with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "losat-licenses"}), timeout=60) as reply:
        return reply.read()


def refresh_rust_inputs() -> int:
    sysroot = Path(subprocess.check_output(["rustc", f"+{RUST_VERSION}", "--print", "sysroot"], text=True).strip())
    version = subprocess.check_output(["rustc", f"+{RUST_VERSION}", "--version"], text=True).strip()
    notice = sysroot / "share/doc/rust/COPYRIGHT-library.html"
    page = notice.read_text(encoding="utf-8")
    listed: dict[tuple[str, str], tuple[str, dict[str, str]]] = {}
    for part in re.split(r"<h3>\N{PACKAGE} ", page[page.index('id="out-of-tree-dependencies"'):])[1:]:
        title = part.split("</h3>")[0]
        name, ver = re.match(r"(.+?)-(\d+\.\d+\.\d+.*)$", title).groups()
        expr = html.unescape(re.search(r"<b>License:</b> (.*?)</p>", part).group(1))
        files = {f: normalize(html.unescape(body)) for f, body in
                 re.findall(r"<summary><code>(.*?)</code></summary>\s*<pre>(.*?)</pre>", part, re.S)}
        listed[(name, ver)] = (expr, files)
    registry: dict[tuple[str, str], tuple[str, dict[str, str]]] = {}
    missing = [(n, v) for n, v, _ in STD_CRATES if (n, v) not in listed]
    if missing:
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "src").mkdir()
            (Path(tmp) / "src/lib.rs").write_text("")
            deps = "\n".join(f'"{n}" = "={v}"' for n, v in missing)
            (Path(tmp) / "Cargo.toml").write_text(f'[package]\nname = "fetch-std-crates"\nversion = "0.0.0"\nedition = "2021"\n\n[dependencies]\n{deps}\n')
            cargo(["fetch"], Path(tmp))
            data = json.loads(cargo(["metadata", "--format-version", "1"], Path(tmp)))
            for package in data["packages"]:
                key = (package["name"], package["version"])
                if key in missing:
                    registry[key] = (package["license"], read_license_files(Path(package["manifest_path"]).parent, package.get("license_file")))
    blocks: list[tuple[dict, str]] = [({"kind": "meta", "channel": RUST_VERSION, "rustc": version,
                                       "origin": "share/doc/rust/COPYRIGHT-library.html and licenses/ of the toolchain; registry crates for the ones it omits; "
                                                 f"{RUST_TAG_URL}/LICENSE-MIT, {RUST_TAG_URL}/{COMPILER_BUILTINS_LICENSE}",
                                       "copyright_library_html_sha256": hashlib.sha256(notice.read_bytes()).hexdigest()}, "")]
    blocks.append(({"kind": "rust-license", "file": "LICENSE-MIT", "origin": f"{RUST_TAG_URL}/LICENSE-MIT"},
                   normalize(fetch(f"{RUST_TAG_URL}/LICENSE-MIT").decode("utf-8"))))
    blocks.append(({"kind": "compiler-builtins-license", "file": "LICENSE.txt", "origin": f"{RUST_TAG_URL}/{COMPILER_BUILTINS_LICENSE}"},
                   normalize(fetch(f"{RUST_TAG_URL}/{COMPILER_BUILTINS_LICENSE}").decode("utf-8"))))
    for path, expr, copyright_lines in STD_IN_TREE:
        blocks.append(({"kind": "in-tree", "path": path, "license": expr, "copyright": list(copyright_lines)}, ""))
    blocks.append(({"kind": "license-text", "license": "Unicode-3.0", "origin": "licenses/Unicode-3.0.txt of the toolchain"},
                   normalize((sysroot / "share/doc/rust/licenses/Unicode-3.0.txt").read_text(encoding="utf-8"))))
    for name, ver, _ in STD_CRATES:
        expr, files = listed.get((name, ver)) or registry[(name, ver)]
        source = "COPYRIGHT-library.html" if (name, ver) in listed else "registry"
        texts = select_texts(name, choose(expr), files)
        for file_name, text in texts.licenses + texts.copyrights + texts.notices:
            blocks.append(({"kind": "crate", "name": name, "version": ver, "license": expr, "file": file_name, "source": source}, text))
    write_bundle(INPUTS / RUST_BUNDLE, blocks)
    standard = fetch(f"{RUST_TAG_URL}/LICENSE-APACHE")
    (INPUTS / APACHE_FILE).write_bytes(standard)
    print(f"{RUST_BUNDLE}: sha256 {sha256_of(INPUTS / RUST_BUNDLE)}\n{APACHE_FILE}: sha256 {sha256_of(INPUTS / APACHE_FILE)}")
    return 0


# (title, license, directory, file in wasi-libc, lines kept or None, notes)
WASI_LIBC_FILES = (
    ("wasi-libc (its own code)", "Apache-2.0 WITH LLVM-exception OR Apache-2.0 OR MIT", "", "LICENSE", None,
     ["The LICENSE file of wasi-libc, which names the license of each directory."]),
    ("wasi-libc (its own code)", "Apache-2.0 WITH LLVM-exception OR Apache-2.0 OR MIT", "", "LICENSE-MIT", None, []),
    ("wasi-libc: musl libc (libc-top-half/musl/)", "MIT", "libc-top-half/musl/", "libc-top-half/musl/COPYRIGHT", None,
     ["musl's COPYRIGHT file; it also names the authors of the math code and other parts."]),
    ("wasi-libc: musl math code, notices of its files (libc-top-half/musl/src/math/)", "MIT", "libc-top-half/musl/src/math/",
     "libc-top-half/musl/src/math/__cosl.c", (1, 14),
     ["The notice of src/math/__cosl.c (Sun Microsystems; Kargl, Schultz and Evans), the same wording as in "
      "the other files taken from FreeBSD's msun."]),
    ("wasi-libc: musl math code, notices of its files (libc-top-half/musl/src/math/)", "MIT", "libc-top-half/musl/src/math/",
     "libc-top-half/musl/src/math/tgammal.c", (1, 15),
     ["The notice of src/math/tgammal.c (Moshier, from OpenBSD), the same wording as in expl.c and the other long double files."]),
    ("wasi-libc: musl math code, notices of its files (libc-top-half/musl/src/math/)", "MIT", "libc-top-half/musl/src/math/",
     "libc-top-half/musl/src/math/exp.c", (1, 6), ["The notice of src/math/exp.c (Arm Limited, from Arm's optimized routines)."]),
    ("wasi-libc: cloudlibc (libc-bottom-half/cloudlibc/)", "BSD-2-Clause", "libc-bottom-half/cloudlibc/",
     "libc-bottom-half/cloudlibc/LICENSE", None, []),
    ("wasi-libc: musl-fts (fts/musl-fts/)", "BSD-3-Clause", "fts/musl-fts/", "fts/musl-fts/COPYING", None, []),
)


def refresh_wasi_libc_input() -> int:
    blocks: list[tuple[dict, str]] = [({"kind": "meta", "commit": WASI_LIBC_COMMIT, "wasi_sdk": "wasi-sdk-27",
                                       "origin": f"{WASI_LIBC_RAW}/"}, "")]
    for title, license_id, directory, path, lines, notes in WASI_LIBC_FILES:
        raw = fetch(f"{WASI_LIBC_RAW}/{path}")
        text = raw.decode("utf-8")
        if lines:
            text = "\n".join(text.split("\n")[lines[0] - 1: lines[1]])
        blocks.append(({"kind": "text", "title": title, "license": license_id,
                        "directory": f"https://github.com/WebAssembly/wasi-libc/tree/{WASI_LIBC_COMMIT}/{directory}",
                        "file": path if not lines else f"{path} lines {lines[0]}-{lines[1]}", "notes": notes,
                        "origin": f"{WASI_LIBC_RAW}/{path}", "origin_sha256": hashlib.sha256(raw).hexdigest()}, normalize(text)))
    write_bundle(INPUTS / WASI_LIBC_BUNDLE, blocks)
    print(f"{WASI_LIBC_BUNDLE}: sha256 {sha256_of(INPUTS / WASI_LIBC_BUNDLE)}")
    return 0


def refresh_zig_input(wheel: Path) -> int:
    import zipfile
    name = wheel.name
    if not name.startswith(f"ziglang-{ZIG_VERSION}-"):
        raise LicenseError(f"{name} is not a ziglang {ZIG_VERSION} wheel")
    blocks: list[tuple[dict, str]] = [({"kind": "meta", "version": ZIG_VERSION, "wheel": name,
                                       "wheel_sha256": sha256_of(wheel)}, "")]
    with zipfile.ZipFile(wheel) as archive:
        for path, title, license_id, repository, notes in ZIG_FILES:
            raw = archive.read(path)
            blocks.append(({"kind": "text", "title": title, "license": license_id, "notes": notes, "repository": repository,
                            "file": path.removeprefix("ziglang/"), "origin_sha256": hashlib.sha256(raw).hexdigest()},
                           normalize(raw.decode("utf-8"))))
    write_bundle(INPUTS / ZIG_BUNDLE, blocks)
    print(f"{ZIG_BUNDLE}: sha256 {sha256_of(INPUTS / ZIG_BUNDLE)}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true", help="fail if a generated file or a copy is stale; write nothing")
    parser.add_argument("--refresh-rust-inputs", action="store_true")
    parser.add_argument("--refresh-wasi-libc-input", action="store_true")
    parser.add_argument("--refresh-zig-input", type=Path, metavar="WHEEL")
    args = parser.parse_args(argv)
    try:
        if args.refresh_rust_inputs:
            return refresh_rust_inputs()
        if args.refresh_wasi_libc_input:
            return refresh_wasi_libc_input()
        if args.refresh_zig_input:
            return refresh_zig_input(args.refresh_zig_input)
        expected = generate()
    except LicenseError as error:
        print(f"third_party_licenses.py: {error}", file=sys.stderr)
        return 1
    pointers = pointer_problems(ROOT)
    if args.check:
        problems = stale(expected, ROOT)
        for problem in problems + pointers:
            print(problem)
        if problems:
            print("Regenerate with: python3 scripts/third_party_licenses.py (and commit the result)")
        if problems or pointers:
            return 1
        print(f"{len(expected)} files up to date")
        return 0
    for relative, data in sorted(expected.items()):
        path = ROOT / relative
        if not path.is_file() or path.read_bytes() != data:
            path.write_bytes(data)
            print(f"wrote {relative}")
    for problem in pointers:
        print(problem)
    return 1 if pointers else 0


if __name__ == "__main__":
    sys.exit(main())
