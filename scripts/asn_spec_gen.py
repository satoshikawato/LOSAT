#!/usr/bin/env python3
"""Generator of LOSAT's ASN.1 type table (`LOSAT/src/report/serial/asn_types/table.rs`).

A maintainer tool: LOSAT's build does not run it, and users need nothing of it. It reads NCBI's
ASN.1 modules (`c++/src/objects/**/*.asn`) and the storage overrides of their `.def` files, takes
the closure of the types that a root type reaches (by default NCBI-ScoreMat's
`PssmWithParameters`, the PSSM file of psiblast and tblastn `-in_pssm`), and writes the table that
LOSAT's generic ASN.1 text reader runs on (`report::serial::asn_types`). Its companion check,
`asn_spec_check.py`, compares the table with the type information that NCBI's datatool generated
(`*_.cpp`, `*_.hpp`) and must pass whenever the table is regenerated.

Usage:
  asn_spec_gen.py <NCBI c++ source dir> [--out <table.rs>] [--root Module:Type ...]

The NCBI source is NCBI BLAST+ 2.17.0's (`ncbi-blast-2.17.0+-src/c++`, the pinned commit 598d8ae6
of AGENTS.md). Without `--out` the table goes to standard output.

What the table holds, per type (named, or anonymous at its place of use):
  - its kind: SEQUENCE (members), CHOICE (variants), SEQUENCE OF / SET OF (element), INTEGER (C++
    storage), ENUMERATED or INTEGER with named values (the value list), BOOLEAN, NULL, REAL,
    VisibleString, UTF8String, StringStore, OCTET STRING;
  - its ASN.1 name ("" for an anonymous type, as NCBI's internal type information reports it) and
    module;
  - per member or variant: its name, its tag (position from 0, NCBI's automatic tags), OPTIONAL,
    DEFAULT (with the default value as written in the module) and the delayed flag (`_delay`).

The storage of an INTEGER is NCBI's C++ type, from the `.def` files: `TSeqPos` and `unsigned` are
Uint4, `TSignedSeqPos` and `int` Int4, `TGi`, `TEntrezId`, `TIntId` and `Int8` Int8 (NCBI_INT8_GI,
corelib/ncbimisc.hpp:981-982), `Int2` Int2.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union

# ---------------------------------------------------------------------------------------------
# The ASN.1 modules
# ---------------------------------------------------------------------------------------------

PRIMITIVES = {
    "INTEGER", "REAL", "BOOLEAN", "NULL", "VisibleString", "UTF8String", "StringStore",
    "OCTET", "BIT", "BigInt", "ENUMERATED", "ObjectDescriptor", "GeneralString",
    "NumericString", "PrintableString", "TeletexString", "VideotexString", "IA5String",
    "GraphicString", "UniversalString", "BMPString", "UTCTime", "GeneralizedTime",
}


@dataclass
class Ref:
    """A reference to a named type."""
    name: str


@dataclass
class Prim:
    """A primitive type: INTEGER, REAL, BOOLEAN, NULL, VisibleString, ..."""
    name: str


@dataclass
class Named:
    """ENUMERATED (`integer` False) or INTEGER with named values (`integer` True)."""
    values: list[tuple[str, int]]
    integer: bool


@dataclass
class ListOf:
    """SEQUENCE OF (`set` False) or SET OF (`set` True)."""
    set: bool
    element: "Type"


@dataclass
class MemberSpec:
    name: str
    type: "Type"
    tag: Optional[int]
    optional: bool
    default: Optional[str]


@dataclass
class Constructed:
    """SEQUENCE (`choice` False), SET (`set_`) or CHOICE (`choice` True)."""
    choice: bool
    members: list[MemberSpec]
    set_: bool = False


Type = Union[Ref, Prim, Named, ListOf, Constructed]


@dataclass
class Module:
    name: str
    path: str
    imports: dict[str, str]
    types: dict[str, Type]


def tokenize(text: str) -> list[str]:
    """The tokens of a module: comments (`--` to `--` or the end of the line) dropped."""
    out: list[str] = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c.isspace():
            i += 1
            continue
        if text.startswith("--", i):
            j = i + 2
            while j < n and text[j] != "\n" and not text.startswith("--", j):
                j += 1
            i = j + 2 if text.startswith("--", j) else j
            continue
        if text.startswith("::=", i):
            out.append("::=")
            i += 3
            continue
        if c in "{}(),;[]|":
            out.append(c)
            i += 1
            continue
        if c == '"':
            j = i + 1
            while j < n:
                if text[j] == '"':
                    if j + 1 < n and text[j + 1] == '"':
                        j += 2
                        continue
                    break
                j += 1
            out.append(text[i:j + 1])
            i = j + 1
            continue
        if c == "'":
            j = text.index("'", i + 1)
            k = j + 1
            if k < n and text[k] in "HB":
                k += 1
            out.append(text[i:k])
            i = k
            continue
        m = re.match(r"[A-Za-z0-9_.\-]+", text[i:])
        if m:
            out.append(m.group(0))
            i += len(m.group(0))
            continue
        out.append(c)
        i += 1
    return out


class Parser:
    def __init__(self, tokens: list[str]) -> None:
        self.t = tokens
        self.p = 0

    def peek(self, k: int = 0) -> Optional[str]:
        return self.t[self.p + k] if self.p + k < len(self.t) else None

    def next(self) -> str:
        tok = self.t[self.p]
        self.p += 1
        return tok

    def expect(self, tok: str) -> None:
        got = self.next()
        if got != tok:
            raise SyntaxError(f"expected {tok!r}, got {got!r} near {self.t[self.p - 5:self.p + 5]}")

    def skip_group(self, open_: str, close: str) -> None:
        self.expect(open_)
        depth = 1
        while depth:
            tok = self.next()
            if tok == open_:
                depth += 1
            elif tok == close:
                depth -= 1

    def modules(self, path: str) -> list[Module]:
        mods: list[Module] = []
        while self.peek() is not None:
            name = self.next()
            if self.peek() == "{":
                self.skip_group("{", "}")
            self.expect("DEFINITIONS")
            while self.peek() != "::=":
                self.next()
            self.expect("::=")
            self.expect("BEGIN")
            imports: dict[str, str] = {}
            types: dict[str, Type] = {}
            while self.peek() != "END":
                tok = self.peek()
                if tok == "EXPORTS":
                    while self.next() != ";":
                        pass
                    continue
                if tok == "IMPORTS":
                    self.next()
                    names: list[str] = []
                    while True:
                        tok = self.next()
                        if tok == ";":
                            break
                        if tok == ",":
                            continue
                        if tok == "FROM":
                            mod = self.next()
                            if self.peek() == "{":
                                self.skip_group("{", "}")
                            for imported in names:
                                imports[imported] = mod
                            names = []
                            continue
                        names.append(tok)
                    continue
                tname = self.next()
                self.expect("::=")
                types[tname] = self.type()
            self.expect("END")
            mods.append(Module(name, path, imports, types))
        return mods

    def type(self) -> Type:
        tok = self.next()
        if tok in ("SEQUENCE", "SET"):
            if self.peek() == "OF":
                self.next()
                return ListOf(tok == "SET", self.type())
            if self.peek() == "(":  # a SIZE constraint
                self.skip_group("(", ")")
                self.expect("OF")
                return ListOf(tok == "SET", self.type())
            return Constructed(False, self.members(), tok == "SET")
        if tok == "CHOICE":
            return Constructed(True, self.members())
        if tok == "ENUMERATED":
            return Named(self.named_values(), False)
        if tok == "INTEGER":
            if self.peek() == "{":
                return Named(self.named_values(), True)
            if self.peek() == "(":
                self.skip_group("(", ")")
            return Prim("INTEGER")
        if tok == "OCTET":
            self.expect("STRING")
            return Prim("OCTET STRING")
        if tok == "BIT":
            self.expect("STRING")
            if self.peek() == "{":
                self.named_values()
            return Prim("BIT STRING")
        if tok in PRIMITIVES:
            if self.peek() == "(":
                self.skip_group("(", ")")
            return Prim(tok)
        if self.peek() == "(":
            self.skip_group("(", ")")
        return Ref(tok)

    def named_values(self) -> list[tuple[str, int]]:
        self.expect("{")
        values: list[tuple[str, int]] = []
        while True:
            if self.peek() == "}":
                self.next()
                break
            name = self.next()
            self.expect("(")
            num = self.next()
            if num == "-":
                num = "-" + self.next()
            self.expect(")")
            values.append((name, int(num)))
            if self.peek() == ",":
                self.next()
        return values

    def members(self) -> list[MemberSpec]:
        self.expect("{")
        members: list[MemberSpec] = []
        while True:
            if self.peek() == "}":
                self.next()
                break
            name = self.next()
            tag: Optional[int] = None
            if self.peek() == "[":
                self.next()
                tag = int(self.next())
                self.expect("]")
                if self.peek() in ("IMPLICIT", "EXPLICIT"):
                    self.next()
            typ = self.type()
            optional, default = False, None
            if self.peek() == "OPTIONAL":
                self.next()
                optional = True
            elif self.peek() == "DEFAULT":
                self.next()
                if self.peek() == "{":
                    start = self.p
                    self.skip_group("{", "}")
                    default = " ".join(self.t[start:self.p])
                else:
                    default = self.next()
                    if default == "-":
                        default = "-" + self.next()
            members.append(MemberSpec(name, typ, tag, optional, default))
            if self.peek() == ",":
                self.next()
        return members


def load_modules(objects: Path) -> dict[str, Module]:
    """Every module of the `.asn` files below `objects` (test and sample files left out)."""
    mods: dict[str, Module] = {}
    for dirpath, _, files in sorted(os.walk(objects)):
        if "samples" in dirpath or "/test" in dirpath:
            continue
        for name in sorted(files):
            if not name.endswith(".asn"):
                continue
            path = Path(dirpath) / name
            text = path.read_text(encoding="latin-1")
            try:
                parsed = Parser(tokenize(text)).modules(str(path.relative_to(objects)))
            except (SyntaxError, ValueError, IndexError):
                # A module of another subset of ASN.1 (an XML-derived or RPC spec); the closure
                # fails by name if it needs one of its types.
                continue
            for module in parsed:
                mods.setdefault(module.name, module)
    return mods


def resolve(mods: dict[str, Module], module: str, name: str) -> str:
    """The module that defines the type `name` seen from `module`."""
    mod = mods[module]
    if name in mod.types:
        return module
    if name in mod.imports:
        return resolve(mods, mod.imports[name], name)
    raise KeyError(f"type {name} is not visible in {module}")


# ---------------------------------------------------------------------------------------------
# The `.def` files: C++ storage overrides and delayed members
# ---------------------------------------------------------------------------------------------

@dataclass
class DefEntries:
    """`[Type]` sections of the `.def` files of a module: `path` (`member.sub.E`, "" for the type
    itself) -> {key -> value} (`_type`, `_storage_type`, `_delay`)."""
    sections: dict[str, dict[str, dict[str, str]]] = field(default_factory=dict)


def load_defs(objects: Path, mods: dict[str, Module]) -> dict[str, DefEntries]:
    """The `.def` entries of each module: those of the `.def` files next to its `.asn` file."""
    by_dir: dict[str, DefEntries] = {}
    out: dict[str, DefEntries] = {}
    for module in mods.values():
        directory = str(Path(module.path).parent)
        if directory not in by_dir:
            entries = DefEntries()
            for def_file in sorted((objects / directory).glob("*.def")):
                section = "-"
                for raw in def_file.read_text(encoding="latin-1").splitlines():
                    line = raw.strip()
                    if not line or line.startswith(";") or line.startswith("#"):
                        continue
                    if line.startswith("[") and line.endswith("]"):
                        section = line[1:-1]
                        continue
                    if "=" not in line:
                        continue
                    key, value = (part.strip() for part in line.split("=", 1))
                    # `Module.Type` sections name a type of one module; plain ones, any module's.
                    type_name = section.split(".", 1)[1] if section.count(".") and \
                        section.split(".", 1)[0] in mods else section
                    if "." in key and key.rsplit(".", 1)[1].startswith("_"):
                        path, attr = key.rsplit(".", 1)
                    elif key.startswith("_"):
                        path, attr = "", key
                    else:
                        continue
                    entries.sections.setdefault(type_name, {}).setdefault(path, {})[attr] = value
            by_dir[directory] = entries
        out[module.name] = by_dir[directory]
    return out


# NCBI's C++ integer types of the `.def` files, as the text reader reads them
# (CObjectIStream::ReadStd, include/serial/impl/objistr.inl; corelib/ncbimisc.hpp for TGi,
# TEntrezId and TIntId with NCBI_INT8_GI).
STORAGE = {
    "int": "Int4",
    "TSignedSeqPos": "Int4",
    "TSeqPos": "Uint4",
    "unsigned": "Uint4",
    "Int2": "Int2",
    "Int8": "Int8",
    "ncbi::TGi": "Int8",
    "ncbi::TEntrezId": "Int8",
    "ncbi::TIntId": "Int8",
}


def storage_of(attrs: dict[str, str]) -> Optional[str]:
    """The integer storage that `_type` (or `_storage_type`) gives, `None` for no override."""
    for key in ("_storage_type", "_type"):
        if key in attrs and attrs[key] in STORAGE:
            return STORAGE[attrs[key]]
    if "_type" in attrs and attrs["_type"] not in ("vector", "list", "set", "multiset"):
        raise ValueError(f"unknown C++ type {attrs['_type']!r}")
    return None


# ---------------------------------------------------------------------------------------------
# The closure and the table
# ---------------------------------------------------------------------------------------------

@dataclass
class Entry:
    """One type of the table."""
    index: int
    kind: str  # Sequence, Choice, List, Integer, Enumerated, Boolean, Null, Real, String, Octets
    name: str = ""  # ASN.1 name; "" for an anonymous type
    module: str = ""
    path: str = ""  # where an anonymous type is used (`Type.member.E`), for the comments
    members: list[tuple[MemberSpec, int, bool]] = field(default_factory=list)  # spec, type, delayed
    element: int = -1
    set_: bool = False
    storage: str = "Int4"
    string: str = "Visible"
    values: list[tuple[str, int]] = field(default_factory=list)
    integer: bool = False
    alias: bool = False  # a named type that is another named type (`A ::= B`)


class Table:
    def __init__(self, mods: dict[str, Module], defs: dict[str, DefEntries]) -> None:
        self.mods = mods
        self.defs = defs
        self.entries: list[Entry] = []
        self.named: dict[tuple[str, str], int] = {}
        self.shared: dict[tuple[str, ...], int] = {}

    def add(self, entry: Entry) -> int:
        entry.index = len(self.entries)
        self.entries.append(entry)
        return entry.index

    def named_type(self, module: str, name: str) -> int:
        key = (module, name)
        if key in self.named:
            return self.named[key]
        typ = self.mods[module].types[name]
        attrs = self.defs[module].sections.get(name, {}).get("", {})
        storage = storage_of(attrs)
        index = self.add(Entry(-1, "pending", name, module))
        self.named[key] = index
        filled = self.build(module, typ, name, name, storage)
        entry = self.entries[index]
        if isinstance(filled, Entry):
            filled.index, filled.name, filled.module = index, name, module
            self.entries[index] = filled
        else:
            # `A ::= B`: A is B under another name.
            entry.kind, entry.alias, entry.element = "Alias", True, filled
        return index

    def build(self, module: str, typ: Type, owner: str, path: str,
              storage: Optional[str]) -> Union[Entry, int]:
        """The entry of `typ` (an `Entry` to place, or the index of a named type it refers to).
        `owner` is the named type whose `.def` section holds the overrides of `path`."""
        if isinstance(typ, Ref):
            target = resolve(self.mods, module, typ.name)
            if typ.name not in self.mods[target].types:
                raise KeyError(typ.name)
            return self.named_type(target, typ.name)
        if isinstance(typ, Prim):
            kinds = {
                "INTEGER": ("Integer", ""), "BOOLEAN": ("Boolean", ""), "NULL": ("Null", ""),
                "REAL": ("Real", ""), "VisibleString": ("String", "Visible"),
                "UTF8String": ("String", "Utf8"), "StringStore": ("String", "Store"),
                "OCTET STRING": ("Octets", ""),
            }
            if typ.name not in kinds:
                raise ValueError(f"{path}: the primitive {typ.name} is not supported")
            kind, string = kinds[typ.name]
            entry = Entry(-1, kind, string=string or "Visible", storage=storage or "Int4")
            return entry
        if isinstance(typ, Named):
            return Entry(-1, "Enumerated", values=typ.values, integer=typ.integer)
        if isinstance(typ, ListOf):
            element = self.place(module, typ.element, owner, path + ".E")
            return Entry(-1, "List", element=element, set_=typ.set)
        if typ.set_:
            raise ValueError(f"{path}: SET (random member order) is not supported")
        entry = Entry(-1, "Choice" if typ.choice else "Sequence")
        section = self.defs[module].sections.get(owner, {})
        for spec in typ.members:
            sub = spec.name if path == owner else path.split(".", 1)[1] + "." + spec.name
            attrs = section.get(sub, {})
            element = self.place(module, spec.type, owner, owner + "." + sub,
                                 storage_of(attrs))
            delayed = attrs.get("_delay") == "1"
            entry.members.append((spec, element, delayed))
        return entry

    def place(self, module: str, typ: Type, owner: str, path: str,
              storage: Optional[str] = None) -> int:
        """The index of `typ` used at `path`: a named type's, a shared primitive's, or a new
        anonymous entry."""
        if storage is None and not isinstance(typ, (Ref, Constructed, ListOf)):
            sub = path.split(".", 1)[1] if "." in path else ""
            storage = storage_of(self.defs[module].sections.get(owner, {}).get(sub, {}))
        built = self.build(module, typ, owner, path, storage)
        if isinstance(built, int):
            return built
        built.path = path
        if built.kind in ("Integer", "Boolean", "Null", "Real", "String", "Octets"):
            key = (built.kind, built.storage if built.kind == "Integer" else "",
                   built.string if built.kind == "String" else "")
            if key not in self.shared:
                built.path = ""
                self.shared[key] = self.add(built)
            return self.shared[key]
        if built.kind == "List":
            key = ("List", str(built.element), str(built.set_))
            if key not in self.shared:
                built.path = ""
                self.shared[key] = self.add(built)
            return self.shared[key]
        return self.add(built)


def const_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "_", name).upper()


def rust_str(text: str) -> str:
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def emit(table: Table, roots: list[tuple[str, str]], source: str) -> str:
    out: list[str] = []
    w = out.append
    w("// GENERATED by scripts/asn_spec_gen.py: do not edit by hand.")
    w("//")
    w("// The ASN.1 types that NCBI's text reader reads under the root types")
    for module, name in roots:
        w(f"// {module}::{name},")
    w("// from the ASN.1 modules and `.def` files of NCBI BLAST+ 2.17.0 (`c++/src/objects`, the pinned")
    w("// commit 598d8ae6 of AGENTS.md). Regenerate and check with:")
    w("//   python3 scripts/asn_spec_gen.py <ncbi-blast-2.17.0+-src/c++> \\")
    w("//       --out LOSAT/src/report/serial/asn_types/table.rs")
    w("//   python3 scripts/asn_spec_check.py <ncbi-blast-2.17.0+-src/c++>")
    w(f"// ({len(table.entries)} types: {len(table.named)} named, the others anonymous at their place of use.)")
    w("")
    w("use super::{Kind, Member, Storage, StringKind, TypeInfo, DEFAULT, DELAYED, OPTIONAL};")
    w("")
    w("/// The index in `TYPES` of each named type (the reader's hooks and the PSSM's typed boundary")
    w("/// name some of them).")
    w("#[allow(dead_code)]")
    w("pub(crate) mod named {")
    seen: set[str] = set()
    for (module, name), index in sorted(table.named.items(), key=lambda item: item[1]):
        const = const_name(name)
        if const in seen:
            const = const_name(module) + "__" + const
        seen.add(const)
        w(f"    pub(crate) const {const}: u16 = {index};")
    w("}")
    w("")
    w(f"pub(crate) static TYPES: [TypeInfo; {len(table.entries)}] = [")
    for e in table.entries:
        label = f"{e.module} {e.name}" if e.name else (e.path or "(shared)")
        w(f"    // {e.index}: {label}")
        name = rust_str(e.name)
        module = rust_str(e.module)
        if e.kind in ("Sequence", "Choice"):
            w(f"    TypeInfo {{ name: {name}, module: {module}, kind: Kind::{e.kind}(&[")
            for spec, element, delayed in e.members:
                flags = []
                if spec.optional:
                    flags.append("OPTIONAL")
                if spec.default is not None:
                    flags.append("DEFAULT")
                if delayed:
                    flags.append("DELAYED")
                flag_text = " | ".join(flags) if flags else "0"
                default = "None" if spec.default is None else f"Some({rust_str(spec.default)})"
                tag = spec.tag if spec.tag is not None else e.members.index((spec, element, delayed))
                w(f"        Member {{ name: {rust_str(spec.name)}, tag: {tag}, ty: {element}, "
                  f"flags: {flag_text}, default: {default} }},")
            w("    ]) },")
        elif e.kind == "List":
            w(f"    TypeInfo {{ name: {name}, module: {module}, kind: Kind::List {{ element: "
              f"{e.element}, set: {'true' if e.set_ else 'false'} }} }},")
        elif e.kind == "Alias":
            w(f"    TypeInfo {{ name: {name}, module: {module}, kind: Kind::Alias({e.element}) }},")
        elif e.kind == "Integer":
            w(f"    TypeInfo {{ name: {name}, module: {module}, kind: Kind::Integer(Storage::{e.storage}) }},")
        elif e.kind == "Enumerated":
            values = ", ".join(f"({rust_str(v)}, {n})" for v, n in e.values)
            w(f"    TypeInfo {{ name: {name}, module: {module}, kind: Kind::Enumerated {{ "
              f"integer: {'true' if e.integer else 'false'}, values: &[{values}] }} }},")
        elif e.kind == "String":
            w(f"    TypeInfo {{ name: {name}, module: {module}, kind: Kind::String(StringKind::{e.string}) }},")
        else:
            w(f"    TypeInfo {{ name: {name}, module: {module}, kind: Kind::{e.kind} }},")
    w("];")
    w("")
    return "\n".join(out)


DEFAULT_ROOTS = [("NCBI-ScoreMat", "PssmWithParameters")]


def build_table(source: Path, roots: list[tuple[str, str]]) -> Table:
    objects = source / "src" / "objects"
    mods = load_modules(objects)
    defs = load_defs(objects, mods)
    table = Table(mods, defs)
    for module, name in roots:
        table.named_type(module, name)
    return table


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("source", type=Path, help="NCBI's c++ source directory")
    parser.add_argument("--out", type=Path, help="the table to write (default: standard output)")
    parser.add_argument("--root", action="append", default=[], metavar="Module:Type")
    args = parser.parse_args(argv)
    roots = [tuple(root.split(":", 1)) for root in args.root] or DEFAULT_ROOTS
    table = build_table(args.source, [(m, t) for m, t in roots])
    text = emit(table, [(m, t) for m, t in roots], str(args.source))
    if args.out:
        args.out.write_text(text)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
