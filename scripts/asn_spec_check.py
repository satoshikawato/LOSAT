#!/usr/bin/env python3
"""Check of LOSAT's ASN.1 type table against the type information NCBI's datatool generated.

The companion of `asn_spec_gen.py` (a maintainer tool; LOSAT's build runs neither). It checks:
  1. the table checked in (`LOSAT/src/report/serial/asn_types/table.rs`) is what the generator
     writes from the same NCBI source;
  2. every SEQUENCE and CHOICE of the table, named or anonymous, against its
     `BEGIN_NAMED_*_INFO` block in NCBI's generated `*_.cpp`: the member or variant names in
     order, the `->SetOptional()` and `->SetDefault(` marks, and the delayed member
     (`->SetDelayBuffer(`);
  3. every value list (ENUMERATED and INTEGER with named values) against its
     `BEGIN_NAMED_ENUM_*INFO` block: the names and numbers in order, and whether it is an INTEGER;
  4. the C++ storage of every INTEGER member and list element against the `typedef` of the
     generated `*_.hpp` (`TSeqPos` Uint4, `int` Int4, `Int8`/`TGi`/`TEntrezId` Int8, ...);
  5. no class of the table is read in random order (`SetRandomOrder(true)`) and all are
     `DataSpec(EDataSpec::eASN)` where the block names a data spec.

Usage:
  asn_spec_check.py <NCBI c++ source dir> [--table <table.rs>]

Exit status 0 when everything agrees, 1 with the differences listed otherwise.
"""
from __future__ import annotations

import argparse
import importlib.util
import re
import sys
from pathlib import Path
from types import ModuleType
from typing import Optional

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TABLE = ROOT / "LOSAT" / "src" / "report" / "serial" / "asn_types" / "table.rs"


def load_generator() -> ModuleType:
    spec = importlib.util.spec_from_file_location("asn_spec_gen", HERE / "asn_spec_gen.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["asn_spec_gen"] = module
    spec.loader.exec_module(module)
    return module


BLOCK = re.compile(
    r'BEGIN_NAMED_(?:BASE_)?(?:CLASS|CHOICE|IMPLICIT_CLASS|BASE_IMPLICIT_CLASS)_INFO\("([^"]*)",'
    r'[^\n]*\n(.*?)\nEND_(?:CLASS|CHOICE)_INFO', re.S)
ENUM = re.compile(
    r'BEGIN_NAMED_ENUM_(?:IN_)?INFO\("([^"]*)",[^\n]*?(true|false)\)\n(.*?)\nEND_ENUM_(?:IN_)?INFO',
    re.S)
INTERNAL = re.compile(r'SET_(?:ENUM_)?INTERNAL_NAME\("([^"]*)", "([^"]*)"\)')
ITEM = re.compile(r'ADD_NAMED_[A-Z_]*(?:MEMBER|VARIANT)\("([^"]+)"')
VALUE = re.compile(r'ADD_ENUM_VALUE\("([^"]+)",\s*([^)]*)\)')
TYPEDEF = re.compile(r'^\s*typedef (.+?) T(\w+);', re.M)

# The C++ types of the generated headers, as the text reader reads them (the generator's
# STORAGE, with the names the headers use).
CXX_STORAGE = {
    "int": "Int4", "TSignedSeqPos": "Int4", "TSeqPos": "Uint4", "unsigned": "Uint4",
    "unsigned int": "Uint4", "Int2": "Int2", "Int8": "Int8", "TGi": "Int8", "TEntrezId": "Int8",
    "TIntId": "Int8", "ncbi::TGi": "Int8", "ncbi::TEntrezId": "Int8", "ncbi::TIntId": "Int8",
    "Int4": "Int4", "NCBI_NS_NCBI::TGi": "Int8", "NCBI_NS_NCBI::TEntrezId": "Int8",
    "NCBI_NS_NCBI::TIntId": "Int8",
}
ENUM_CONST = re.compile(r"\b(e\w+)\s*=\s*(-?\d+)")


def header_constants(include: Path, type_name: str) -> dict[str, int]:
    """The enum constants (`eRepr_raw = 2`) of the generated header of `type_name`."""
    out: dict[str, int] = {}
    for path in include.glob(f"*/{type_name.replace('-', '_')}_.hpp"):
        for m in ENUM_CONST.finditer(path.read_text(encoding="latin-1")):
            out.setdefault(m.group(1), int(m.group(2)))
    return out


def blocks_by_key(objects: Path) -> tuple[dict[str, list[str]], dict[str, list[tuple[str, str]]]]:
    """The class/choice blocks and the enum blocks of the generated files, by type name, or by
    `Owner.member` for an anonymous type (`SET_INTERNAL_NAME`)."""
    classes: dict[str, list[str]] = {}
    enums: dict[str, list[tuple[str, str]]] = {}
    for path in sorted(objects.glob("*/*_.cpp")):
        text = path.read_text(encoding="latin-1")
        for m in BLOCK.finditer(text):
            name, body = m.group(1), m.group(2)
            internal = INTERNAL.search(body)
            key = name if name else (f"{internal.group(1)}.{internal.group(2)}" if internal else "")
            classes.setdefault(key, []).append(body)
        for m in ENUM.finditer(text):
            name, integer, body = m.group(1), m.group(2), m.group(3)
            internal = INTERNAL.search(body)
            key = name if name else (f"{internal.group(1)}.{internal.group(2)}" if internal else "")
            enums.setdefault(key, []).append((integer, body))
    return classes, enums


def header_typedefs(include: Path, type_name: str) -> dict[str, set[str]]:
    """`T<Member>` -> the C++ types of its typedefs in the generated header of `type_name`."""
    out: dict[str, set[str]] = {}
    for path in include.glob(f"*/{type_name.replace('-', '_')}_.hpp"):
        for m in TYPEDEF.finditer(path.read_text(encoding="latin-1")):
            out.setdefault(m.group(2), set()).add(m.group(1).strip())
    return out


def typedef_name(member: str) -> str:
    return member[0].upper() + member[1:].replace("-", "_")


def element_cxx(cxx: str) -> Optional[str]:
    m = re.fullmatch(r"(?:vector|list|set|multiset)< (.+) >", cxx)
    return m.group(1).strip() if m else None


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("source", type=Path, help="NCBI's c++ source directory")
    parser.add_argument("--table", type=Path, default=TABLE)
    args = parser.parse_args(argv)
    gen = load_generator()
    table = gen.build_table(args.source, gen.DEFAULT_ROOTS)
    problems: list[str] = []

    # 1. The table checked in is the generator's.
    text = gen.emit(table, gen.DEFAULT_ROOTS, str(args.source))
    if not args.table.exists() or args.table.read_text() != text:
        problems.append(f"{args.table} differs from the generator's output: regenerate it")

    objects = args.source / "src" / "objects"
    include = args.source / "include" / "objects"
    classes, enums = blocks_by_key(objects)
    checked = {"classes": 0, "lists": 0, "storage": 0}

    def owner_of(entry_path: str) -> str:
        return entry_path.split(".", 1)[0]

    for e in table.entries:
        key = e.name or e.path
        if e.kind in ("Sequence", "Choice"):
            candidates = classes.get(key, [])
            if not candidates:
                problems.append(f"{key}: no generated class or choice block")
                continue
            want = [spec.name for spec, _, _ in e.members]
            found = False
            for body in candidates:
                lines = [ln for ln in body.split("\n") if ITEM.search(ln)]
                names = [ITEM.search(ln).group(1) for ln in lines]  # type: ignore[union-attr]
                if names != want:
                    continue
                marks: list[str] = []
                for (spec, _, delayed), ln in zip(e.members, lines):
                    if e.kind == "Sequence":
                        optional = "->SetOptional()" in ln or "->SetDefault(" in ln
                        if (spec.optional or spec.default is not None) != optional:
                            marks.append(f"{spec.name}: OPTIONAL/DEFAULT")
                        if (spec.default is not None) != ("->SetDefault(" in ln):
                            marks.append(f"{spec.name}: DEFAULT")
                    if delayed != ("->SetDelayBuffer(" in ln):
                        marks.append(f"{spec.name}: delayed")
                if "SetRandomOrder(true)" in body:
                    marks.append("random order")
                if "DataSpec(" in body and "DataSpec(ncbi::EDataSpec::eASN)" not in body:
                    marks.append("data spec")
                if marks:
                    problems.append(f"{key}: " + ", ".join(marks))
                found = True
                break
            if not found:
                problems.append(f"{key}: the members {want} are not those of NCBI's block")
            checked["classes"] += 1
            # 4. The storage of the INTEGER members and list elements.
            owner = e.name or owner_of(e.path)
            typedefs = header_typedefs(include, owner)
            for spec, ty, _ in e.members:
                target = table.entries[ty]
                cxx = typedefs.get(typedef_name(spec.name), set())
                if target.kind == "Integer" and not target.name:
                    got = {CXX_STORAGE.get(t, t) for t in cxx}
                    if got and got != {target.storage}:
                        problems.append(f"{key}.{spec.name}: storage {target.storage}, NCBI {sorted(cxx)}")
                    checked["storage"] += 1 if got else 0
                elif target.kind == "List" and not target.name:
                    element = table.entries[target.element]
                    if element.kind == "Integer" and not element.name:
                        got = {CXX_STORAGE.get(x, x) for x in (element_cxx(t) for t in cxx) if x}
                        if got and got != {element.storage}:
                            problems.append(f"{key}.{spec.name}.E: storage {element.storage}, NCBI {sorted(cxx)}")
                        checked["storage"] += 1 if got else 0
        elif e.kind == "Enumerated":
            candidates_enum = enums.get(key, [])
            if not candidates_enum:
                problems.append(f"{key}: no generated value list")
                continue
            want_values = [(name, number) for name, number in e.values]
            constants = header_constants(include, e.name or owner_of(e.path))
            ok = False
            for integer, body in candidates_enum:
                values: list[tuple[str, Optional[int]]] = []
                for name, number in VALUE.findall(body):
                    number = number.strip()
                    if re.fullmatch(r"-?\d+", number):
                        values.append((name, int(number)))
                    else:
                        values.append((name, constants.get(number.rsplit("::", 1)[-1])))
                if values == want_values and (integer == "true") == e.integer:
                    ok = True
                    break
            if not ok:
                problems.append(f"{key}: the value list or its INTEGER mark differs from NCBI's")
            checked["lists"] += 1

    for problem in problems:
        print(problem)
    print(f"{checked['classes']} SEQUENCE/CHOICE types, {checked['lists']} value lists, "
          f"{checked['storage']} INTEGER storages checked: {len(problems)} differences")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
