#!/usr/bin/env python3
"""Tests for the program-edge ratchet (`check_program_edges.py`).

No NCBI counterpart: the checker reads LOSAT's own module structure.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
import tempfile
import unittest
from pathlib import Path


CHECKER_PATH = Path(__file__).with_name("check_program_edges.py")
SPEC = importlib.util.spec_from_file_location("program_edges_ratchet", CHECKER_PATH)
assert SPEC is not None and SPEC.loader is not None
ratchet = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ratchet
SPEC.loader.exec_module(ratchet)

PROGRAMS = ("blastn", "blastp", "blastx", "common", "psiblast", "tblastn", "tblastx")


def write(root: Path, rel: str, text: str) -> None:
    path = root / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def make_tree(root: Path, files: dict[str, str]) -> None:
    """A crate with `algorithm::<program>` for every program; `files` (under LOSAT/src/algorithm/)
    replace or add files."""
    write(root, "LOSAT/src/lib.rs", "pub mod algorithm;\npub mod utils;\n")
    write(root, "LOSAT/src/utils.rs", "pub fn h() {}\n")
    write(root, "LOSAT/src/algorithm/mod.rs", "".join(f"pub mod {p};\n" for p in PROGRAMS))
    for program in PROGRAMS:
        write(root, f"LOSAT/src/algorithm/{program}/mod.rs", "pub fn f() {}\n")
    for rel, text in files.items():
        write(root, f"LOSAT/src/algorithm/{rel}", text)


def edges(root: Path) -> dict[tuple[str, str], tuple[str, int, int]]:
    """(src, dst) -> (scope, code sites, test sites)."""
    found = ratchet.find_edges(root)
    return {key: (edge.scope, len(edge.code), len(edge.test)) for key, edge in found.items()}


def run(root: Path, rows: list[str] | None, *extra: str) -> tuple[int, str]:
    if rows is not None:
        write(root, "LOSAT/tests/program_edges.tsv",
              "# comment\nsrc\tdst\tscope\n" + "".join(f"{row}\n" for row in rows))
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
        code = ratchet.main(["--root", str(root), *extra])
    return code, out.getvalue()


class MaskingTests(unittest.TestCase):
    def test_masking_keeps_offsets_and_newlines(self) -> None:
        text = 'let a = "x\ny"; // c\n/* a /* b */ c */ let b = \'"\';\n'
        masked = ratchet.mask(text)
        self.assertEqual(len(masked), len(text))
        self.assertEqual(masked.count("\n"), text.count("\n"))
        self.assertNotIn("x", masked)
        self.assertIn("let b", masked)
        self.assertNotIn('"', masked)

    def test_mentions_in_comments_doc_comments_and_strings_are_not_references(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_tree(root, {"blastn/mod.rs": (
                "//! Uses crate::algorithm::blastp::a (inner doc comment).\n"
                "/// See crate::algorithm::blastp::b.\n"
                "// crate::algorithm::blastp::c\n"
                "/* outer /* nested crate::algorithm::blastp::d */ crate::algorithm::blastp::e */\n"
                "/** crate::algorithm::blastp::f */\n"
                'const S: &str = "crate::algorithm::blastp::g \\" crate::algorithm::blastp::h";\n'
                # Read as a plain string, the mention would sit between two strings, in code.
                'const R: &str = r#"a " crate::algorithm::blastp::i " b"#;\n'
                'const B: &[u8] = br"crate::algorithm::blastp::j";\n'
                'const BS: &[u8] = b"crate::algorithm::blastp::k";\n'
                'const CS: &core::ffi::CStr = c"crate::algorithm::blastp::l";\n'
                "const Q: char = '\"';\n"
                "fn keep<'a>(x: &'a str) -> &'a str { x }\n"
                "pub fn real() { crate::algorithm::common::m(); }\n"
            )})
            self.assertEqual(edges(root), {("blastn", "common"): ("code", 1, 0)})

    def test_dollar_crate_is_not_a_reference(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_tree(root, {"blastn/mod.rs": (
                "macro_rules! m { () => { $crate::algorithm::blastp::f() }; }\n"
                "macro_rules! n { () => { $crate :: algorithm::blastp::f() }; }\n"
            )})
            self.assertEqual(edges(root), {})


class ResolutionTests(unittest.TestCase):
    def test_a_use_group_gives_one_edge_per_program(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_tree(root, {"tblastx/mod.rs": (
                "use crate::algorithm::{\n"
                "    blastn::{a, b::c},\n"
                "    common::{self, d as e},\n"
                "    tblastx::own,\n"
                "};\n"
                "use crate::algorithm::blastn::*;\n"
                "use crate::utils::h;\n"
            )})
            self.assertEqual(
                edges(root),
                {("tblastx", "blastn"): ("code", 3, 0), ("tblastx", "common"): ("code", 2, 0)},
            )

    def test_super_chains_from_an_inline_test_module(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_tree(root, {"blastp/mod.rs": (
                "use super::super::utils::h;\n"  # the crate root's utils: no edge
                "use super::common::g;\n"
                "mod inner {\n"
                "    pub use super::super::blastn::f;\n"
                "}\n"
                "#[cfg(test)]\n"
                "mod tests {\n"
                "    use super::super::tblastx::f;\n"
                "    use super::inner::f as same_program;\n"
                "}\n"
            )})
            self.assertEqual(
                edges(root),
                {
                    ("blastp", "blastn"): ("code", 1, 0),
                    ("blastp", "common"): ("code", 1, 0),
                    ("blastp", "tblastx"): ("test", 0, 1),
                },
            )

    def test_self_super_is_a_super_chain(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_tree(root, {"blastn/mod.rs": (
                "mod inner {\n"
                "    pub use self::super::super::common::f;\n"
                "}\n"
            )})
            self.assertEqual(edges(root), {("blastn", "common"): ("code", 1, 0)})

    def test_super_chains_from_a_path_child(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_tree(root, {
                "blastx/mod.rs": "pub mod search;\npub mod query;\n",
                "blastx/query.rs": "pub fn q() {}\n",
                # algorithm::blastx::search::flow, though the file sits next to search.rs
                "blastx/search.rs": '#[path = "x_flow.rs"]\nmod flow;\n',
                "blastx/x_flow.rs": (
                    "use super::super::query::q;\n"
                    "use super::super::super::tblastx::f;\n"
                    "use super::super::super::super::utils::h;\n"
                ),
            })
            self.assertEqual(edges(root), {("blastx", "tblastx"): ("code", 1, 0)})


class AlgorithmModuleTests(unittest.TestCase):
    """A path that names `algorithm` itself would hide the program it reaches."""

    FORMS = {
        "glob": "use super::*;\npub fn g() { blastx::f(); }\n",
        "bare import": "use crate::algorithm;\npub fn g() { algorithm::blastx::f(); }\n",
        "alias": "use crate::algorithm as a;\npub fn g() { a::blastx::f(); }\n",
    }

    def test_each_form_is_unexpected(self) -> None:
        for form, text in self.FORMS.items():
            with self.subTest(form), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                make_tree(root, {"common/mod.rs": text})
                self.assertEqual(edges(root), {})
                code, out = run(root, [])
                self.assertEqual(code, 1, out)
                self.assertIn(
                    "UNEXPECTED\tcommon -> algorithm\tcode\tLOSAT/src/algorithm/common/mod.rs:1\t"
                    "crate::algorithm\tnames the algorithm module itself",
                    out,
                )
                self.assertIn("edges=0\tlisted=0\tunexpected=1\tstale=0\tupward=0", out)

    def test_an_item_of_algorithm_mod_rs_is_unexpected_and_a_child_module_is_not(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_tree(root, {
                "mod.rs": "".join(f"pub mod {p};\n" for p in PROGRAMS)
                + "pub mod shared;\npub use blastx::f as g;\n",
                "shared.rs": "pub fn s() {}\n",
                "blastn/mod.rs": (
                    "pub fn h() { crate::algorithm::g(); }\n"
                    "pub fn i() { crate::algorithm::shared::s(); }\n"
                    "#[cfg(test)]\nmod tests {\n    use super::*;\n}\n"
                ),
            })
            graph = ratchet.scan(root)
            self.assertEqual(graph.edges, {})
            self.assertEqual(list(graph.opaque), ["blastn"])
            self.assertEqual(
                [site.reference() for site in graph.opaque["blastn"].code], ["crate::algorithm::g"]
            )


class TestCodeTests(unittest.TestCase):
    def test_a_file_declared_by_cfg_test_mod_counts_as_test(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_tree(root, {
                "blastn/mod.rs": "#[cfg(test)]\nmod t;\nmod inner;\n",
                "blastn/t.rs": "use crate::algorithm::common::x;\nmod deeper;\n",
                "blastn/t/deeper.rs": "use crate::algorithm::blastp::y;\n",
                "blastn/inner.rs": "#![cfg(test)]\nuse crate::algorithm::tblastx::z;\n",
            })
            self.assertEqual(
                edges(root),
                {
                    ("blastn", "blastp"): ("test", 0, 1),
                    ("blastn", "common"): ("test", 0, 1),
                    ("blastn", "tblastx"): ("test", 0, 1),
                },
            )

    def test_cfg_predicates_tests_statements_and_fields(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_tree(root, {"tblastn/mod.rs": (
                '#[cfg(any(test, target_arch = "wasm32"))]\n'
                "use crate::algorithm::common::a;\n"
                "#[cfg(not(test))]\n"
                "use crate::algorithm::tblastx::b;\n"
                '#[cfg(all(test, feature = "parallel"))]\n'
                "use crate::algorithm::blastp::c;\n"
                "#[test]\n"
                "#[ignore]\n"
                "fn t() { crate::algorithm::blastn::d(); }\n"
                "struct S {\n"
                "    #[cfg(test)]\n"
                "    probe: std::collections::HashMap<u32, crate::algorithm::blastx::E>,\n"
                "    kept: u32,\n"
                "}\n"
                "fn g() {\n"
                "    #[cfg(test)]\n"
                "    let _ = crate::algorithm::psiblast::h();\n"
                "    #[cfg(test)]\n"
                "    if true { crate::algorithm::psiblast::i(); }\n"
                "    crate::algorithm::psiblast::j();\n"
                "}\n"
            )})
            self.assertEqual(
                edges(root),
                {
                    ("tblastn", "blastn"): ("test", 0, 1),
                    ("tblastn", "blastp"): ("test", 0, 1),
                    ("tblastn", "blastx"): ("test", 0, 1),
                    ("tblastn", "common"): ("code", 1, 0),
                    ("tblastn", "psiblast"): ("code", 1, 2),
                    ("tblastn", "tblastx"): ("code", 1, 0),
                },
            )

    def test_which_cfg_predicates_require_test(self) -> None:
        requires = ratchet.cfg_requires_test
        self.assertTrue(requires("test"))
        self.assertTrue(requires("all(test,unix)"))
        self.assertTrue(requires("all(feature=,any(test,all(test,unix)))"))
        self.assertTrue(requires("any(test,all(test,unix))"))
        self.assertFalse(requires("any(test,target_arch=)"))
        self.assertFalse(requires("not(test)"))
        self.assertFalse(requires("all()"))
        self.assertFalse(requires("unix"))


class OutcomeTests(unittest.TestCase):
    def tree(self, root: Path) -> None:
        """blastp -> blastn (code), tblastx -> common (test only)."""
        make_tree(root, {
            "blastp/mod.rs": "pub use crate::algorithm::blastn::f;\n",
            "tblastx/mod.rs": (
                "#[cfg(test)]\nmod tests {\n    use crate::algorithm::common::f;\n}\n"
            ),
        })

    def test_listed_edges_pass(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.tree(root)
            code, out = run(root, ["blastp\tblastn\tcode", "tblastx\tcommon\ttest"])
            self.assertEqual(code, 0, out)
            self.assertEqual(
                out.strip(),
                "program-edge ratchet\tedges=2\tlisted=2\tunexpected=0\tstale=0\tupward=0",
            )

    def test_a_new_edge_is_unexpected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.tree(root)
            code, out = run(root, ["blastp\tblastn\tcode"])
            self.assertEqual(code, 1)
            self.assertIn(
                "UNEXPECTED\ttblastx -> common\ttest\tLOSAT/src/algorithm/tblastx/mod.rs:3\t"
                "crate::algorithm::common::f\t",
                out,
            )
            self.assertIn("unexpected=1\tstale=0\tupward=0", out)

    def test_a_listed_test_edge_that_gains_code_is_unexpected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.tree(root)
            code, out = run(root, ["blastp\tblastn\ttest", "tblastx\tcommon\ttest"])
            self.assertEqual(code, 1)
            self.assertIn("UNEXPECTED\tblastp -> blastn\tcode (listed as test only)\t", out)
            self.assertIn("unexpected=1\tstale=0", out)

    def test_a_gone_edge_or_a_narrower_scope_is_stale(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.tree(root)
            rows = ["blastp\tblastn\tcode", "tblastx\tcommon\tcode", "tblastn\tblastp\tcode"]
            code, out = run(root, rows)
            self.assertEqual(code, 1)
            self.assertIn("STALE\ttblastn -> blastp\tcode\t-\t-\tlisted but gone", out)
            self.assertIn(
                "STALE\ttblastx -> common\tcode\tLOSAT/src/algorithm/tblastx/mod.rs:3\t", out
            )
            self.assertIn("unexpected=0\tstale=2\tupward=0", out)

    def test_an_upward_edge_fails_even_when_listed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_tree(root, {
                "blastn/mod.rs": "use crate::algorithm::tblastx::f;\n",
                "common/mod.rs": "use crate::algorithm::psiblast::f;\n",
                "blastp/mod.rs": "use crate::algorithm::psiblast::f;\n",
                "psiblast/mod.rs": "use crate::algorithm::blastp::f;\n",
            })
            listed = ["blastn\ttblastx\tcode", "blastp\tpsiblast\tcode", "psiblast\tblastp\tcode"]
            code, out = run(root, listed)
            self.assertEqual(code, 1)
            self.assertIn(
                "UPWARD\tblastn -> tblastx\tcode\tLOSAT/src/algorithm/blastn/mod.rs:1\t", out
            )
            # Nothing leaves common, even to a program outside LAYERS (reported once, as UPWARD).
            self.assertIn("UPWARD\tcommon -> psiblast\tcode\t", out)
            self.assertNotIn("UNEXPECTED", out)
            self.assertIn("edges=4\tlisted=3\tunexpected=0\tstale=0\tupward=2", out)

    def test_list_prints_every_edge_and_exits_0(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.tree(root)
            code, out = run(root, None, "--list")
            self.assertEqual(code, 0)
            self.assertEqual(
                out.splitlines(),
                [
                    "src\tdst\tscope\tcode_sites\ttest_sites\tfirst_site\treference",
                    "blastp\tblastn\tcode\t1\t0\tLOSAT/src/algorithm/blastp/mod.rs:1\t"
                    "crate::algorithm::blastn::f",
                    "tblastx\tcommon\ttest\t0\t1\tLOSAT/src/algorithm/tblastx/mod.rs:3\t"
                    "crate::algorithm::common::f",
                ],
            )

    def test_a_malformed_list_or_tree_exits_2(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.tree(root)
            for rows in (
                ["blastp\tblastn\tcode\textra"],
                ["blastp\tblastn\tsometimes"],
                ["blastp\tblastp\tcode"],
                ["blastp\tblastn\tcode", "blastp\tblastn\ttest"],
            ):
                self.assertEqual(run(root, rows)[0], 2, rows)
            write(root, "LOSAT/tests/program_edges.tsv", "blastp\tblastn\tcode\n")
            code, out = run(root, None)
            self.assertEqual(code, 2)
            self.assertIn("header", out)
            (root / "LOSAT/src/lib.rs").unlink()
            self.assertEqual(run(root, ["blastp\tblastn\tcode"])[0], 2)

    def test_a_module_under_algorithm_that_cannot_be_followed_exits_2(self) -> None:
        cases = {
            "no file found for `mod missing;`": {"blastn/mod.rs": "mod missing;\n"},
            "takes its path from `cfg_attr`": {
                "blastn/mod.rs": '#[cfg_attr(unix, path = "imp_unix.rs")]\nmod imp;\n',
                "blastn/imp.rs": "",
                "blastn/imp_unix.rs": "use crate::algorithm::blastp::f;\n",
            },
        }
        for message, files in cases.items():
            with self.subTest(message), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                make_tree(root, files)
                code, out = run(root, [])
                self.assertEqual(code, 2, out)
                self.assertIn("LOSAT/src/algorithm/blastn/mod.rs:", out)
                self.assertIn(message, out)

    def test_an_unresolved_module_outside_algorithm_is_not_read(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            make_tree(root, {})
            write(root, "LOSAT/src/lib.rs", "pub mod algorithm;\nmod elsewhere;\n")
            self.assertEqual(run(root, [])[0], 0)


if __name__ == "__main__":
    unittest.main()
