#!/usr/bin/env python3
"""Tests for the inline test-module check and mover (`check_inline_test_modules.py`).

No NCBI counterpart: the checker reads LOSAT's own source files.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def load(name: str, file: str):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


checker = load("inline_test_modules", "check_inline_test_modules.py")
guard = checker.crate

PATTERNS = ["tests/", "src/**/*_tests.rs", "src/**/tests.rs", "src/**/tests_*.rs", "src/**/*_tests/"]
ISOLATED_GIT = {"GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1"}


def cargo_toml() -> str:
    lines = ", ".join(f'"{p}"' for p in PATTERNS)
    return f'[package]\nname = "demo"\nversion = "0.1.0"\nedition = "2021"\nexclude = [{lines}]\n'


class Tree:
    """A Git repository with a LOSAT/ crate of the given files."""

    def __init__(self, files: dict[str, str]) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        subprocess.run(["git", "init", "-q", str(self.root)], check=True)
        for path, text in {"Cargo.toml": cargo_toml(), **files}.items():
            self.write(path, text)

    def write(self, path: str, text: str) -> None:
        target = self.root / "LOSAT" / path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text, encoding="utf-8")

    def read(self, path: str) -> str:
        return (self.root / "LOSAT" / path).read_text(encoding="utf-8")

    def exists(self, path: str) -> bool:
        return (self.root / "LOSAT" / path).exists()

    def run(self, *options: str) -> tuple[int, str]:
        out = io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
            code = checker.main(["--root", str(self.root), *options])
        return code, out.getvalue()

    def guard(self) -> tuple[int, str]:
        out = io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
            code = guard.main(["--root", str(self.root)])
        return code, out.getvalue()

    def __enter__(self) -> "Tree":
        return self

    def __exit__(self, *exc: object) -> None:
        self.tmp.cleanup()


def kinds(text: str) -> list[tuple[str, str]]:
    return [(t.kind, text[t.start:t.end]) for t in checker.lex(text) if t.kind != "ws"]


class LexerTests(unittest.TestCase):
    def test_tokens_cover_the_text(self) -> None:
        text = 'fn f<\'a>(x: &\'a str) -> char { /* a /* b */ { */ \'}\' }\nconst S: &str = r#"} "x" {"#;\n'
        self.assertEqual("".join(text[t.start:t.end] for t in checker.lex(text)), text)

    def test_braces_in_literals_and_comments_are_not_tokens(self) -> None:
        text = "{ '{' b'}' \"}\" r#\"}\"# br##\"}\"#\"## // }\n /* } /* } */ } */ }"
        braces = [t for t in checker.lex(text) if t.kind == "punct" and text[t.start] in "{}"]
        self.assertEqual([text[t.start] for t in braces], ["{", "}"])

    def test_lifetimes_are_not_chars(self) -> None:
        got = kinds("<'a, 'static> 'x' '\\'' '\\u{1F600}' b'\\\\'")
        self.assertEqual([k for k, _ in got if k in ("life", "char")],
                         ["life", "life", "char", "char", "char", "char"])

    def test_escapes_and_raw_identifiers(self) -> None:
        got = kinds('"a\\"}" r#match "\\\\"')
        self.assertEqual([t for k, t in got if k == "str"], ['"a\\"}"', '"\\\\"'])
        self.assertIn(("ident", "r#match"), got)

    def test_an_unterminated_literal_is_an_error(self) -> None:
        for text in ('"abc', "/* abc", 'r#"abc"'):
            with self.assertRaises(checker.LexError):
                checker.lex(text)

    def test_lines_that_start_inside_a_literal(self) -> None:
        text = 'let a = "x\nin string\n";\n/* c\n in comment */\nlet b = r#"\nraw"#;\nend\n'
        inside = checker.literal_lines(text, checker.lex(text))
        self.assertEqual(sorted(inside), [1, 2, 4, 6])


class FindTests(unittest.TestCase):
    def names(self, text: str) -> list[tuple[str, tuple[str, ...], str | None]]:
        return [(f.name, f.enclosing, f.nested_in_other) for f in checker.find_modules(text)]

    def test_what_is_test_only(self) -> None:
        text = (
            "#[cfg(test)]\nmod a {}\n"
            "#[cfg(all(test, unix))]\npub(crate) mod b {}\n"
            "#[cfg(any(test, feature = \"x\"))]\nmod c {}\n"
            "#[cfg(not(test))]\nmod d {}\n"
            "#[cfg(feature = \"x\")]\nmod e {}\n"
            "mod f {\n    #![cfg(test)]\n}\n"
            "#[cfg(test)]\n#[allow(dead_code)]\n/// doc\npub(in crate::x) mod g {}\n"
            "#[cfg(test)]\nmod h;\n"
        )
        self.assertEqual([n for n, _, _ in self.names(text)], ["a", "b", "f", "g"])

    def test_braces_in_strings_and_comments_do_not_end_the_module(self) -> None:
        text = ('#[cfg(test)]\nmod tests {\n    const A: &str = "}";\n    const B: &str = r#"}"#;\n'
                "    // }\n    /* } */\n    const C: char = '}';\n    fn f<'a>(_: &'a str) {}\n}\nfn after() {}\n")
        (found,) = checker.find_modules(text)
        self.assertEqual(text[found.close:], "}\nfn after() {}\n")

    def test_nesting(self) -> None:
        text = (
            "mod outer {\n    #[cfg(test)]\n    mod tests {\n        #[cfg(test)]\n        mod deeper {}\n    }\n}\n"
            "fn f() {\n    #[cfg(test)]\n    mod in_fn {}\n}\n"
            "impl X {\n    #[cfg(test)]\n    mod in_impl {}\n}\n"
        )
        self.assertEqual(self.names(text), [("tests", ("outer",), None), ("in_fn", (), "fn f ( )"), ("in_impl", (), "impl X")])

    def test_a_declaration_is_not_inline(self) -> None:
        self.assertEqual(self.names("#[cfg(test)]\n#[path = \"x.rs\"]\nmod tests;\n"), [])


class FixTests(unittest.TestCase):
    PARENT = (
        "//! Module docs.\n\npub fn product() -> u32 {\n    1\n}\n\n"
        "/// Docs of the test module.\n#[cfg(test)]\n#[allow(clippy::module_inception)]\nmod tests {\n"
        "    use super::*;\n\n"
        "    const TABLE: &str = \"a\nb\n  c\";\n"
        "    const RAW: &str = r#\"\nraw } line\n\"#;\n"
        "    /* a block\n  comment with odd indent\n    */\n"
        "    #[test]\n    fn product_is_one() {\n        assert_eq!(product(), 1);\n"
        "        let text = \"one \\\n  two\";\n        assert!(!text.is_empty());\n    }\n}\n"
    )

    def setUp(self) -> None:
        patcher = mock.patch.dict(os.environ, ISOLATED_GIT)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_check_names_the_file_and_line_and_fix_moves_the_module(self) -> None:
        with Tree({"src/lib.rs": "mod engine;\n", "src/engine.rs": self.PARENT}) as tree:
            code, out = tree.run()
            self.assertEqual(code, 1)
            self.assertIn("LOSAT/src/engine.rs:10\tinline test-only module `tests`", out)
            code, out = tree.run("--fix")
            self.assertEqual(code, 0, out)
            parent = tree.read("src/engine.rs")
            self.assertIn("#[cfg(test)]\n#[allow(clippy::module_inception)]\n#[path = \"engine_tests.rs\"]\nmod tests;\n", parent)
            self.assertIn("/// Docs of the test module.\n", parent)
            self.assertTrue(parent.startswith(self.PARENT[:self.PARENT.index("mod tests")]))
            sibling = tree.read("src/engine_tests.rs")
            # one level off, except where a line starts inside a literal or a block comment
            self.assertTrue(sibling.startswith("use super::*;\n\nconst TABLE: &str = \"a\nb\n  c\";\n"))
            self.assertIn("const RAW: &str = r#\"\nraw } line\n\"#;\n", sibling)
            self.assertIn("/* a block\n  comment with odd indent\n    */\n", sibling)
            self.assertIn("    let text = \"one \\\n  two\";\n", sibling)
            self.assertIn("#[test]\nfn product_is_one() {\n    assert_eq!(product(), 1);\n", sibling)
            self.assertEqual(tree.run()[0], 0)
            self.assertEqual(tree.guard()[0], 0, tree.guard()[1])
            before = {p: tree.read(p) for p in ("src/engine.rs", "src/engine_tests.rs")}
            self.assertEqual(tree.run("--fix")[0], 0)  # idempotent
            self.assertEqual({p: tree.read(p) for p in before}, before)

    def test_file_names(self) -> None:
        body = "#[cfg(test)]\nmod {name} {{\n    fn f() {{}}\n}}\n"
        files = {
            "src/lib.rs": "mod a;\nmod b;\nmod sub;\n",
            "src/a.rs": body.format(name="tests") + body.format(name="boundary_tests") + body.format(name="extra"),
            "src/b.rs": "pub mod inner {\n" + "    #[cfg(test)]\n    mod tests {\n        fn f() {}\n    }\n}\n",
            "src/sub/mod.rs": body.format(name="tests") + body.format(name="more"),
        }
        with Tree(files) as tree:
            code, out = tree.run("--fix")
            self.assertEqual(code, 0, out)
            for path in ("src/a_tests.rs", "src/a_boundary_tests.rs", "src/a_extra_tests.rs", "src/b_inner_tests.rs",
                         "src/sub/tests.rs", "src/sub/more_tests.rs"):
                self.assertTrue(tree.exists(path), path)
            # in `b.rs` (not a mod.rs) the inline module `inner` is looked for under src/b/inner/
            self.assertIn("#[path = \"../../b_inner_tests.rs\"]\n    mod tests;", tree.read("src/b.rs"))
            self.assertEqual(tree.guard()[0], 0, tree.guard()[1])

    def test_a_nested_module_in_a_mod_rs_file_is_looked_for_under_its_directory(self) -> None:
        files = {"src/lib.rs": "mod sub;\n",
                 "src/sub/mod.rs": "mod inner {\n    #[cfg(test)]\n    mod tests {\n        fn f() {}\n    }\n}\n"}
        with Tree(files) as tree:
            self.assertEqual(tree.run("--fix")[0], 0)
            # `inner` in a mod.rs is looked for under src/sub/inner/
            self.assertIn("#[path = \"../inner_tests.rs\"]", tree.read("src/sub/mod.rs"))
            self.assertTrue(tree.exists("src/sub/inner_tests.rs"))
            self.assertEqual(tree.guard()[0], 0, tree.guard()[1])

    def test_not_test_only_modules_stay(self) -> None:
        text = "#[cfg(any(test, feature = \"x\"))]\nmod helper {\n    fn f() {}\n}\n"
        with Tree({"src/lib.rs": text}) as tree:
            self.assertEqual(tree.run()[0], 0)
            self.assertEqual(tree.run("--fix")[0], 0)
            self.assertEqual(tree.read("src/lib.rs"), text)

    def test_a_taken_name_takes_the_one_fallback(self) -> None:
        files = {"src/lib.rs": "mod a;\n", "src/a.rs": "#[cfg(test)]\nmod tests {\n    fn f() {}\n}\n"
                 "#[cfg(test)]\nmod a_tests;\n", "src/a_tests.rs": ""}
        with Tree(files) as tree:
            self.assertEqual(tree.run("--fix")[0], 0)
            self.assertTrue(tree.exists("src/a_inline_tests.rs"))
            self.assertEqual(tree.read("src/a_tests.rs"), "")
        files["src/a_inline_tests.rs"] = ""
        with Tree(files) as tree:
            code, out = tree.run("--fix")
            self.assertEqual(code, 1)
            self.assertIn("are taken", out)

    def test_refusals_change_nothing(self) -> None:
        cases = {
            "in a function": "fn f() {\n    #[cfg(test)]\n    mod tests {\n        fn g() {}\n    }\n}\n",
            "a mod declaration": "#[cfg(test)]\nmod tests {\n    mod data;\n}\n",
            "a path attribute": "#[cfg(test)]\nmod tests {\n    #[path = \"x.rs\"]\n    mod data;\n}\n",
            "file!": "#[cfg(test)]\nmod tests {\n    const F: &str = file!();\n}\n",
            "line!": "#[cfg(test)]\nmod tests {\n    const F: u32 = line!();\n}\n",
            "one line": "#[cfg(test)]\nmod tests { fn f() {} }\n",
            "under-indented": "#[cfg(test)]\nmod tests {\n    fn f() {\n  }\n}\n",
        }
        ok = "#[cfg(test)]\nmod ok_tests {\n    fn f() {}\n}\n"
        for what, text in cases.items():
            with self.subTest(what), Tree({"src/lib.rs": "mod a;\nmod b;\n", "src/a.rs": text, "src/b.rs": ok}) as tree:
                code, out = tree.run("--fix")
                self.assertEqual(code, 1, out)
                self.assertIn("LOSAT/src/a.rs", out)
                self.assertIn("nothing was changed", out)
                self.assertEqual(tree.read("src/b.rs"), ok)
                self.assertFalse(tree.exists("src/b_ok_tests.rs") or tree.exists("src/b_tests.rs"))

    def test_line_comparison_operators_are_not_macros(self) -> None:
        text = "#[cfg(test)]\nmod tests {\n    fn f(line: &str) -> bool {\n        line != \"}\"\n    }\n}\n"
        with Tree({"src/lib.rs": text}) as tree:
            self.assertEqual(tree.run("--fix")[0], 0)

    def test_the_ledger_follows_the_moved_tests(self) -> None:
        ledger = ("module\tcase\tstatus\tlosat_test\tnote\tsince\n"
                  "m\tc\tported\tLOSAT/src/a.rs::t_one;LOSAT/src/other.rs::t_two\tn\t2026-11-21\n")
        files = {"src/lib.rs": "mod a;\nmod other;\n",
                 "src/a.rs": "#[cfg(test)]\nmod tests {\n    #[test]\n    fn t_one() {}\n}\n",
                 "src/other.rs": "fn t_two() {}\n"}
        with Tree(files) as tree:
            (tree.root / "docs/evidence/ncbi_unit_cases").mkdir(parents=True)
            (tree.root / checker.LEDGER).write_text(ledger, encoding="utf-8")
            code, out = tree.run("--fix")
            self.assertEqual(code, 0, out)
            self.assertEqual((tree.root / checker.LEDGER).read_text(encoding="utf-8"),
                             ledger.replace("src/a.rs::t_one", "src/a_tests.rs::t_one"))

    def test_the_self_check_reads_the_original_back(self) -> None:
        for text in (self.PARENT, FixTests.PARENT.replace("mod tests {", "pub(crate) mod tests {")):
            found = checker.find_modules(text)
            moves = checker.plan_file("src/engine.rs", text, found, lambda: set(), set())
            new = checker.apply_moves(text, moves)
            self.assertNotEqual(new, text)
            self.assertEqual(checker.restore(new, checker.sibling_loader(moves)), text)

    def test_the_whole_text_is_kept_but_for_the_indentation(self) -> None:
        with Tree({"src/lib.rs": self.PARENT}) as tree:
            tree.run("--fix")
            sibling = tree.read("src/lib_tests.rs")
            original = self.PARENT[self.PARENT.index("{\n    use super") + 2:self.PARENT.rindex("}\n")]
            lost = [(a, b) for a, b in zip(original.splitlines(), sibling.splitlines()) if a != b and a.strip() != b.strip()]
            self.assertEqual(lost, [])
            self.assertEqual(len(original.splitlines()), len(sibling.splitlines()))


if __name__ == "__main__":
    unittest.main()
