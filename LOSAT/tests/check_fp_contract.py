#!/usr/bin/env python3
"""Check LOSAT's floating-point contract (docs/agents/diagnostics.md, "Floating-point contract").

NCBI reference (598d8ae6): c++/src/algo/blast/core/boost_erf.c:143-143
```c
            result *= expl(-z * z) / z;
```
LOSAT's one fused multiply-add (`f64::mul_add` in `stats/spouge.rs`) is the exact rounding error
of a product inside the double-double emulation of that long-double line. Every other
floating-point operation is a separate IEEE operation in the order of the NCBI source, so the
compiler must not contract `a * b + c` into a fused multiply-add, and nothing may add one.

Two checks:

- `source`: the `mul_add` calls under `LOSAT/src` are the ones in `MUL_ADD_SITES`; no FMA
  intrinsic, relaxed-SIMD multiply-add, fast-math operation, or function of the `libm` crate is
  named; the cargo
  configuration sets no fast-math LLVM option and no `relaxed-simd` target feature.
- `binary --losat <bin>`: `objdump -d -C` of a release binary. An FMA instruction may appear
  only in a function under `ALLOWED_FUNCTION_PREFIXES`; a reference to the library `fma` (a
  direct call, or on x86-64 a load of the GOT slot that holds its address) only in one under
  `MUL_ADD_FUNCTION_PREFIXES`; and at least one function under `MUL_ADD_FUNCTION_PREFIXES` must
  hold one or the other. Prints one line per function and a total line in the form
  `<binary>\tfma_instructions=N\tfma_refs=M`.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

REPO = Path(__file__).resolve().parents[2]

# No NCBI counterpart: the ledger of the explicit fused multiply-adds, keyed by file and the
# function that contains the call. `x_product_error` is the production site; the test checks
# Dekker's product against it.
MUL_ADD_SITES = {
    ("LOSAT/src/stats/spouge.rs", "x_product_error"): 1,
    ("LOSAT/src/stats/spouge.rs", "x_dekker_product_error_matches_fused_multiply_add"): 1,
}

# No NCBI counterpart: the release-binary functions that `x_product_error` is inlined into
# (the ErfC code of `stats/spouge.rs`). With an FMA target (x86-64-v3, aarch64) `mul_add` is an
# FMA instruction there; on baseline x86-64 it is a call of `fma`.
MUL_ADD_FUNCTION_PREFIXES = ("LOSAT::stats::spouge::",)
# No NCBI counterpart: functions that may hold an FMA instruction: the `mul_add` sites and
# the `fma` of Rust's compiler-builtins that baseline x86-64 links, which takes the FMA or FMA4
# instruction when the CPU has one and a software fused multiply-add otherwise (all three round
# once, as IEEE 754 fusedMultiplyAdd).
ALLOWED_FUNCTION_PREFIXES = MUL_ADD_FUNCTION_PREFIXES + ("compiler_builtins::math::libm_math::arch::x86::fma::",)

# No NCBI counterpart: names that would fuse or reassociate floating-point operations.
FORBIDDEN_SOURCE_NAMES = re.compile(
    r"\b(?:_mm(?:256|512)?_mask[z3]?_fn?m(?:add|sub)\w*"
    r"|_mm(?:256|512)?_fn?m(?:add|sub)(?:sub|add)?_\w+"
    r"|vfm[as]q?_\w+"
    r"|v?fmlal\w*|vfmlsl\w*"
    r"|\w*relaxed_n?madd\w*"
    r"|f(?:add|sub|mul|div|rem)_(?:fast|algebraic)"
    r"|algebraic_(?:add|sub|mul|div|rem)"
    r"|fmaf(?:16|32|64|128)|fmuladdf(?:16|32|64|128))\b"
)
FORBIDDEN_CONFIG = re.compile(
    r"fp-contract|fast-math|unsafe-fp-math|no-nans-fp-math|no-infs-fp-math"
    r"|no-signed-zeros-fp-math|approx-func|relaxed-simd|denormal-fp-math"
)
CONFIG_FILES = ("LOSAT/.cargo/config.toml", "LOSAT/Cargo.toml")
# No NCBI counterpart: the transcendental functions come from the platform's libm, as NCBI's do.
LIBM_CRATE_RE = re.compile(r"\blibm::\w+")

FUNCTION_RE = re.compile(r"\bfn\s+([A-Za-z_]\w*)\s*(?:<[^>{}]*>)?\s*\(")
MUL_ADD_RE = re.compile(r"\.\s*mul_add\s*\(|\bf(?:32|64)::mul_add\b")

# x86-64 FMA3 (vfmadd132sd, vfnmsub231ps, vfmaddsub213pd, ...) and FMA4 (vfmaddsd); the
# AVX-512 forms use the same mnemonics.
X86_FMA_RE = re.compile(r"^vfn?m(?:add|sub)(?:sub|add)?(?:132|213|231)?[ps][sdh]$")
# AArch64 scalar and vector fused forms.
ARM_FMA_RE = re.compile(
    r"^(?:fn?madd|fn?msub|fn?mla|fn?mls|fmlal2?|fmlsl2?|bfmla[bt]?|bfmlal[bt]|fn?mad|fn?msb)$"
)
FUNCTION_HEADER_RE = re.compile(r"^[0-9a-f]+ <(.+)>:$")
INSTRUCTION_RE = re.compile(r"^\s+[0-9a-f]+:\s+(\S+)\s*(.*)$")
FMA_CALL_RE = re.compile(r"<fma(?:@plt)?>")
CALL_MNEMONICS = frozenset({"call", "callq", "jmp", "jmpq", "bl", "b"})
# `mov 0x...(%rip),%reg        # 79ee68 <_DYNAMIC+0x2d8>`: the address a RIP-relative operand reads.
RIP_COMMENT_RE = re.compile(r"#\s*([0-9a-f]+)\b")
SYMBOL_RE = re.compile(r"^([0-9a-f]+)\s.*\sF\s.*\s(?:\.hidden\s+)?fma$")
RELATIVE_RE = re.compile(r"^([0-9a-f]+)\s+R_X86_64_RELATIVE\s+\*ABS\*\+0x([0-9a-f]+)$")


RAW_STRING_RE = re.compile(r'(?<![A-Za-z0-9_])b?r(#*)"')


def blank(text: str) -> str:
    return "".join("\n" if ch == "\n" else " " for ch in text)


def strip_comments(text: str) -> str:
    """Blank out `//` and `/* */` comments and the bodies of strings, keeping line breaks."""
    out = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        raw = RAW_STRING_RE.match(text, i) if c in "br" else None
        if raw:
            close = '"' + raw.group(1)
            j = text.find(close, raw.end())
            j = n if j < 0 else j
            out.append(raw.group(0) + blank(text[raw.end():j]) + close)
            i = j + len(close)
        elif text.startswith("//", i):
            j = text.find("\n", i)
            j = n if j < 0 else j
            i = j
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            out.append(blank(text[i:j]))
            i = j
        elif c == "'" and (text.startswith("'\"'", i) or text.startswith("'\\\"'", i)):
            # The character literals '"' and '\"' do not open a string.
            j = text.find("'", i + 1) + 1
            out.append(text[i:j])
            i = j
        elif c == '"':
            j = i + 1
            while j < n and text[j] != '"':
                j += 2 if text[j] == "\\" else 1
            out.append('"' + blank(text[i + 1:j]) + '"')
            i = j + 1
        else:
            out.append(c)
            i += 1
    return "".join(out)


def mul_add_sites(path: str, text: str) -> dict[tuple[str, str], int]:
    """Count `mul_add` calls by (file, enclosing function)."""
    code = strip_comments(text)
    functions = [(m.start(), m.group(1)) for m in FUNCTION_RE.finditer(code)]
    sites: dict[tuple[str, str], int] = {}
    for match in MUL_ADD_RE.finditer(code):
        owner = "<top level>"
        for start, name in functions:
            if start > match.start():
                break
            owner = name
        key = (path, owner)
        sites[key] = sites.get(key, 0) + 1
    return sites


def check_source(root: Path) -> list[str]:
    errors: list[str] = []
    found: dict[tuple[str, str], int] = {}
    for file in sorted((root / "LOSAT/src").rglob("*.rs")):
        rel = file.relative_to(root).as_posix()
        text = file.read_text(encoding="utf-8", errors="replace")
        for key, count in mul_add_sites(rel, text).items():
            found[key] = found.get(key, 0) + count
        code = strip_comments(text)
        for line_no, line in enumerate(code.splitlines(), 1):
            for match in FORBIDDEN_SOURCE_NAMES.finditer(line):
                errors.append(f"{rel}:{line_no}: `{match.group(0)}` fuses or reassociates floating-point operations")
            for match in LIBM_CRATE_RE.finditer(line):
                errors.append(f"{rel}:{line_no}: `{match.group(0)}` calls the libm crate instead of the platform's libm")
    if found != MUL_ADD_SITES:
        for key in sorted(set(found) | set(MUL_ADD_SITES)):
            have, want = found.get(key, 0), MUL_ADD_SITES.get(key, 0)
            if have != want:
                errors.append(f"{key[0]}: fn {key[1]}: {have} mul_add call(s), the ledger MUL_ADD_SITES says {want}")
    for rel in CONFIG_FILES:
        path = root / rel
        if not path.is_file():
            continue
        for line_no, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            setting = line.split("#", 1)[0]
            match = FORBIDDEN_CONFIG.search(setting)
            if match:
                errors.append(f"{rel}:{line_no}: `{match.group(0)}` changes floating-point evaluation")
    return errors


@dataclass
class FunctionCounts:
    instructions: int = 0
    refs: int = 0
    mnemonics: dict[str, int] = field(default_factory=dict)


def count_fma(lines: Iterable[str], arch: str, fma_slots: frozenset[int] = frozenset()) -> dict[str, FunctionCounts]:
    """FMA instructions and references to `fma` per function of `objdump -d` text.

    `fma_slots` are the addresses of the GOT slots that hold the address of `fma` (x86-64 PIE
    binaries call it through one: `mov slot(%rip),%r14` ... `call *%r14`).
    """
    fma_re = ARM_FMA_RE if arch == "aarch64" else X86_FMA_RE
    counts: dict[str, FunctionCounts] = {}
    function = "<none>"
    for line in lines:
        header = FUNCTION_HEADER_RE.match(line)
        if header:
            function = header.group(1)
            continue
        insn = INSTRUCTION_RE.match(line)
        if not insn:
            continue
        mnemonic, operands = insn.group(1), insn.group(2)
        if fma_re.match(mnemonic):
            entry = counts.setdefault(function, FunctionCounts())
            entry.instructions += 1
            entry.mnemonics[mnemonic] = entry.mnemonics.get(mnemonic, 0) + 1
        elif mnemonic in CALL_MNEMONICS and FMA_CALL_RE.search(operands):
            counts.setdefault(function, FunctionCounts()).refs += 1
        elif fma_slots:
            rip = RIP_COMMENT_RE.search(operands)
            if rip and int(rip.group(1), 16) in fma_slots:
                counts.setdefault(function, FunctionCounts()).refs += 1
    return counts


def fma_got_slots(binary: Path, objdump: str) -> frozenset[int]:
    """GOT slots whose dynamic relocation stores the address of the binary's own `fma`."""
    symbols = subprocess.run([objdump, "-t", str(binary)], check=True, capture_output=True, text=True).stdout
    addresses = {int(m.group(1), 16) for m in map(SYMBOL_RE.match, symbols.splitlines()) if m}
    if not addresses:
        return frozenset()
    relocs = subprocess.run([objdump, "-R", str(binary)], capture_output=True, text=True).stdout
    return frozenset(
        int(m.group(1), 16)
        for m in map(RELATIVE_RE.match, relocs.splitlines())
        if m and int(m.group(2), 16) in addresses
    )


def binary_arch(binary: Path, objdump: str) -> str:
    header = subprocess.run([objdump, "-f", str(binary)], check=True, capture_output=True, text=True).stdout
    if "aarch64" in header:
        return "aarch64"
    if "x86-64" in header:
        return "x86_64"
    raise SystemExit(f"{binary}: unsupported architecture in `objdump -f`:\n{header}")


def judge(counts: dict[str, FunctionCounts]) -> list[str]:
    errors = []
    sites = 0
    for function, entry in sorted(counts.items()):
        if function.startswith(MUL_ADD_FUNCTION_PREFIXES):
            sites += 1
        elif entry.refs:
            errors.append(f"{function}: {entry.refs} reference(s) to fma outside the documented mul_add sites")
        elif not function.startswith(ALLOWED_FUNCTION_PREFIXES):
            errors.append(f"{function}: {entry.instructions} FMA instruction(s) outside the documented mul_add sites")
    if sites == 0:
        errors.append(
            "no function under " + ", ".join(MUL_ADD_FUNCTION_PREFIXES)
            + " holds an FMA instruction or a reference to fma: the mul_add site was not found "
            "(objdump output format, or a moved site; update the ledger and the contract)"
        )
    return errors


def report(binary: Path, counts: dict[str, FunctionCounts]) -> str:
    lines = []
    for function, entry in sorted(counts.items()):
        kinds = ",".join(f"{m}:{c}" for m, c in sorted(entry.mnemonics.items()))
        lines.append(f"{function}\tfma_instructions={entry.instructions}\tfma_refs={entry.refs}\t{kinds}")
    total_insn = sum(e.instructions for e in counts.values())
    total_refs = sum(e.refs for e in counts.values())
    lines.append(f"{binary}\tfma_instructions={total_insn}\tfma_refs={total_refs}\tfunctions={len(counts)}")
    return "\n".join(lines) + "\n"


def check_binary(binary: Path, objdump: str, out: Path | None) -> list[str]:
    arch = binary_arch(binary, objdump)
    slots = fma_got_slots(binary, objdump) if arch == "x86_64" else frozenset()
    proc = subprocess.Popen(
        [objdump, "-d", "-C", "--no-show-raw-insn", str(binary)], stdout=subprocess.PIPE, text=True
    )
    assert proc.stdout is not None
    counts = count_fma(proc.stdout, arch, slots)
    if proc.wait() != 0:
        raise SystemExit(f"objdump failed on {binary} (exit {proc.returncode})")
    text = report(binary, counts)
    sys.stdout.write(f"arch={arch}\n{text}")
    if out is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text, encoding="utf-8")
    return judge(counts)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    source = sub.add_parser("source", help="check the sources and the cargo configuration")
    source.add_argument("--root", type=Path, default=REPO)
    binary = sub.add_parser("binary", help="count FMA instructions and references to fma in a release binary")
    binary.add_argument("--losat", type=Path, required=True)
    binary.add_argument("--objdump", default="objdump")
    binary.add_argument("--out", type=Path, help="also write the counts to this file")
    args = parser.parse_args(argv)
    if args.command == "source":
        errors = check_source(args.root)
    else:
        errors = check_binary(args.losat, args.objdump, args.out)
    for error in errors:
        print(f"FAIL {error}", file=sys.stderr)
    if not errors:
        print(f"floating-point contract: {args.command} check passed")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
