#!/usr/bin/env python3
"""A PSSM whose query nests one of its recursive ASN.1 types (or several) some levels deep.

usage: gen_nests.py SKELETON KIND D [D2 [D3]] > OUT.asn

SKELETON is a PSSM in ASN.1 text whose `query` member, at the start of a line indented by four spaces, is
followed by its `intermediateData` (plain or gzipped; LOSAT/tests/fixtures/pssm/d3om_r6_nest_mix_6600.asn.gz is
one); the generated query takes its place. D is the number of times the value around a value is applied (the KIND's values in the nest are D,
or D + 1 where the innermost value is of the KIND too, marked +1):

  set        Bioseq-sets around the query Bioseq (`set { seq-set { ... } }`)
  mix        Seq-loc mixes in a feature's location (`mix { ... }`)
  equiv      Seq-loc equivs in a feature's location (`equiv { ... }`)
  user       User-objects in a descriptor, through `data object` (+1)
  uobjs      User-objects in a descriptor, through `data objects { ... }` (+1)
  ufield     User-fields of a descriptor's User-object, through `data fields { ... }` (+1)
  align      Seq-aligns in an annotation, through `segs disc { ... }`, around a `std` segment (+1)
  pub        Pubs in a descriptor, through `equiv { ... }` (+1)
  variation  Variation-refs of a feature, through `data set { variations { ... } }` (+1)
  vcons      Variation-refs of a feature, through `consequence { variation ... }` (+1)
  multidata  SeqTable-multi-data of a Seq-table column, through `int-scaled { data ... }` (+1)
  mdelta     SeqTable-multi-data of a Seq-table column, through `int-delta` (+1)
  orgname    OrgNames of a descriptor's BioSource, through `name hybrid { ... }` (+1)
  sets_users D Bioseq-sets around a descriptor's D2 + 1 User-objects (`user` D2)
  sets_aligns_mixes  D Bioseq-sets around D2 + 1 Seq-aligns around a `std` segment's D3 Seq-loc mixes
  alignuser  D + 1 Seq-aligns (`disc`) around D2 + 1 User-objects in the innermost alignment's `ext`
  vrpub      D + 1 Variation-refs (`variations`) around D2 Pub equivs in the innermost one's `pub`

The fixture rows' inputs of line RG item D3-1 are made with it (the README of this folder).
"""
import gzip
import sys
from pathlib import Path

Q1 = 'local str "Query_1"'
STD = "ncbistdaa '0C1414080E0A0F031610140D0C0806140F030E16100714120C\n080A14160F'H"
INT = 'int { from 0, to 1, id local str "Query_1" }'
USER0 = '{ type str "t", data { { label str "a", data str "x" } } }'
SET = "set {{ seq-set {{ {} }} }}"


def bioseq(descr: str = "", annots: tuple[str, ...] = ()) -> str:
    parts = ["id { " + Q1 + " }", 'descr { title "rnd random"' + (", " + descr if descr else "") + " }",
             f"inst {{ repr raw, mol aa, length 30, seq-data {STD} }}"]
    if annots:
        parts.append("annot { " + ", ".join(annots) + " }")
    return "seq { " + ", ".join(parts) + " }"


def nest(d: int, base: str, wrap: str) -> str:
    s = base
    for _ in range(d):
        s = wrap.format(s)
    return s


def feature_loc(loc: str) -> str:
    return bioseq(annots=('{ data ftable { { data region "x", location ' + loc + " } } }",))


def variation(v: str) -> str:
    return bioseq(annots=("{ data ftable { { data variation " + v + ", location " + INT + " } } }",))


def table(m: str) -> str:
    return bioseq(annots=("{ data seq-table { feat-type 1, num-rows 1, columns { { header { field-id location-from }, "
                          "data " + m + " }, { header { field-id location-id }, default id " + Q1 + " } } } }",))


def std(loc: str = INT) -> str:
    return "{ type global, dim 2, segs std { { dim 2, loc { " + loc + ", " + INT + " } } } }"


def query(kind: str, d: int, d2: int = 0, d3: int = 0) -> str:
    disc = "{{ type disc, segs disc {{ {} }} }}"
    if kind == "set":
        return nest(d, bioseq(), SET)
    if kind in ("mix", "equiv"):
        return feature_loc(nest(d, INT, kind + " {{ {} }}"))
    if kind == "user":
        return bioseq(descr="user " + nest(d, USER0, '{{ type str "t", data {{ {{ label str "a", data object {} }} }} }}'))
    if kind == "uobjs":
        return bioseq(descr="user " + nest(d, USER0, '{{ type str "t", data {{ {{ label str "a", data objects {{ {} }} }} }} }}'))
    if kind == "ufield":
        f = nest(d, '{ label str "a", data str "x" }', '{{ label str "a", data fields {{ {} }} }}')
        return bioseq(descr='user { type str "t", data { ' + f + " } }")
    if kind == "align":
        return bioseq(annots=("{ data align { " + nest(d, std(), disc) + " } }",))
    if kind == "pub":
        return bioseq(descr="pub { pub { " + nest(d, 'gen { cit "x" }', "equiv {{ {} }}") + " } }")
    if kind == "variation":
        return variation(nest(d, "{ data unknown NULL }", "{{ data set {{ type unknown, variations {{ {} }} }} }}"))
    if kind == "vcons":
        return variation(nest(d, "{ data unknown NULL }", "{{ data unknown NULL, consequence {{ variation {} }} }}"))
    if kind == "multidata":
        return table(nest(d, "int { 0 }", "int-scaled {{ mul 1, add 0, data {} }}"))
    if kind == "mdelta":
        return table(nest(d, "int { 0 }", "int-delta {}"))
    if kind == "orgname":
        o = nest(d, '{ name binomial { genus "A" } }', "{{ name hybrid {{ {} }} }}")
        return bioseq(descr='source { org { taxname "x", orgname ' + o + " } }")
    if kind == "sets_users":
        return nest(d, query("user", d2), SET)
    if kind == "sets_aligns_mixes":
        a = nest(d2, std(nest(d3, INT, "mix {{ {} }}")), disc)
        return nest(d, bioseq(annots=("{ data align { " + a + " } }",)), SET)
    if kind == "alignuser":
        u = nest(d2, USER0, '{{ type str "t", data {{ {{ label str "a", data object {} }} }} }}')
        a = "{ type global, dim 2, segs std { { dim 2, loc { " + INT + ", " + INT + " } } }, ext { " + u + " } }"
        return bioseq(annots=("{ data align { " + nest(d, a, disc) + " } }",))
    if kind == "vrpub":
        p = nest(d2, 'gen { cit "x" }', "equiv {{ {} }}")
        return variation(nest(d, "{ pub " + p + ", data unknown NULL }",
                              "{{ data set {{ type unknown, variations {{ {} }} }} }}"))
    raise SystemExit(f"unknown kind {kind}")


def main() -> None:
    skeleton, kind, *depths = sys.argv[1:]
    raw = Path(skeleton).read_bytes()
    text = (gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw).decode()
    # The query runs from `    query ` to the next member of `pssm`, `    intermediateData`.
    start = text.index("\n    query ") + 1
    end = text.index("\n    intermediateData", start) + 1
    d = [int(x) for x in depths] + [0, 0]
    out = text[:start] + "    query " + query(kind, d[0], d[1], d[2]) + ",\n" + text[end:]
    sys.stdout.buffer.write(out.encode())


if __name__ == "__main__":
    main()
