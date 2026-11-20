#!/usr/bin/env python3
"""Tests for the crate test-file check (`check_crate_test_files.py`).

No NCBI counterpart: the checker reads LOSAT's own packaging (`exclude` in LOSAT/Cargo.toml).
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


CHECKER_PATH = Path(__file__).with_name("check_crate_test_files.py")
SPEC = importlib.util.spec_from_file_location("crate_test_files", CHECKER_PATH)
assert SPEC is not None and SPEC.loader is not None
checker = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = checker
SPEC.loader.exec_module(checker)

PATTERNS = ["tests/", "src/**/*_tests.rs", "src/**/tests.rs", "src/**/tests_*.rs", "src/**/*_tests/"]


def cargo_toml(exclude: list[str]) -> str:
    lines = ", ".join(f'"{p}"' for p in exclude)
    return f'[package]\nname = "demo"\nversion = "0.1.0"\nedition = "2021"\nexclude = [{lines}]\n'


class CrateTree:
    """A Git repository with a LOSAT/ crate made of the given files (untracked: Cargo lists them too)."""

    def __init__(self, files: dict[str, str], exclude: list[str] | None = None) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        subprocess.run(["git", "init", "-q", str(self.root)], check=True)
        files = {"Cargo.toml": cargo_toml(PATTERNS if exclude is None else exclude), **files}
        for path, text in files.items():
            target = self.root / "LOSAT" / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(text, encoding="utf-8")

    def run(self, *options: str) -> tuple[int, str]:
        out = io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
            code = checker.main(["--root", str(self.root), *options])
        return code, out.getvalue()

    def __enter__(self) -> "CrateTree":
        return self

    def __exit__(self, *exc: object) -> None:
        self.tmp.cleanup()


# A user's own Git configuration (a global excludes file) must not change what the tests see.
ISOLATED_GIT = {"GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1"}


@mock.patch.dict(os.environ, ISOLATED_GIT)
class CheckTests(unittest.TestCase):
    def test_a_clean_tree_passes(self) -> None:
        files = {
            "src/lib.rs": (
                "mod engine;\n"
                "#[cfg(test)]\nmod engine_tests;\n"
                "#[cfg(all(test, unix))]\n#[path = \"helpers/test_dir.rs\"]\nmod test_dir;\n"
            ),
            "src/main.rs": "mod cli;\nfn main() {}\n",
            "src/cli.rs": "#[cfg(test)]\nmod tests;\n",
            "src/cli/tests.rs": "mod more;\n#[test]\nfn t() {}\n",
            "src/cli/tests/more.rs": "",
            "src/engine.rs": (
                "pub(crate) const TABLE: &[u8] = include_bytes!(\"table.bin\");\n"
                "#[cfg(test)]\nmod tests_cases;\n"
                "#[cfg(test)]\nmod tests {\n    const BIG: [u8; 1] = include!(\"engine_table.inc.rs\");\n}\n"
            ),
            "src/table.bin": "data",
            "src/engine_table.inc.rs": "[1]",
            "src/engine/tests_cases.rs": "",
            "src/engine_tests.rs": "#[path = \"fixtures_tests/a.rs\"]\nmod a;\n",
            "src/fixtures_tests/a.rs": "",
            "src/helpers/test_dir.rs": "",
            "tests/integration.rs": "",
        }
        exclude = PATTERNS + ["src/helpers/test_dir.rs", "src/engine_table.inc.rs"]
        with CrateTree(files, exclude) as tree:
            code, text = tree.run()
        self.assertEqual(code, 0, text)
        self.assertIn("the 7 .rs files that only tests read are left out", text)

    def test_a_file_the_build_reads_that_a_pattern_catches_fails(self) -> None:
        files = {"src/lib.rs": "mod tests_util;\n", "src/tests_util.rs": ""}
        with CrateTree(files) as tree:
            code, text = tree.run()
        self.assertEqual(code, 1, text)
        self.assertIn("LOSAT/src/tests_util.rs\tthe build reads it (src/lib.rs:1)", text)
        self.assertIn("'src/**/tests_*.rs'", text)

    def test_a_data_file_the_build_includes_that_a_pattern_catches_fails(self) -> None:
        files = {"src/lib.rs": "const T: &[u8] = include_bytes!(\"data_tests/t.bin\");\n",
                 "src/data_tests/t.bin": "x"}
        with CrateTree(files) as tree:
            code, text = tree.run()
        self.assertEqual(code, 1, text)
        self.assertIn("LOSAT/src/data_tests/t.bin\tthe build reads it", text)
        self.assertIn("'src/**/*_tests/'", text)

    def test_a_new_test_file_left_in_the_package_fails(self) -> None:
        files = {"src/lib.rs": "#[cfg(test)]\nmod parser_checks;\n", "src/parser_checks.rs": ""}
        with CrateTree(files) as tree:
            code, text = tree.run()
        self.assertEqual(code, 1, text)
        self.assertIn("LOSAT/src/parser_checks.rs\tonly tests read it (src/lib.rs:2)", text)

    def test_a_table_only_a_test_includes_must_be_excluded(self) -> None:
        files = {"src/lib.rs": "#[cfg(test)]\nmod tests {\n    const T: [u8; 1] = include!(\"t.inc.rs\");\n}\n",
                 "src/t.inc.rs": "[1]"}
        with CrateTree(files) as tree:
            code, text = tree.run()
        self.assertEqual(code, 1, text)
        self.assertIn("LOSAT/src/t.inc.rs\tonly tests read it (src/lib.rs:3)", text)

    def test_cfg_any_test_is_read_by_the_build(self) -> None:
        files = {"src/lib.rs": "#[cfg(any(test, feature = \"diag\"))]\nmod diag_tests;\n",
                 "src/diag_tests.rs": ""}
        with CrateTree(files) as tree:
            code, text = tree.run()
        self.assertEqual(code, 1, text)
        self.assertIn("LOSAT/src/diag_tests.rs\tthe build reads it", text)

    def test_a_file_reached_by_tests_and_by_the_build_is_read_by_the_build(self) -> None:
        files = {"src/lib.rs": "#[cfg(test)]\n#[path = \"shared.rs\"]\nmod shared_for_tests;\nmod shared;\n",
                 "src/shared.rs": ""}
        with CrateTree(files) as tree:
            code, text = tree.run()
        self.assertEqual(code, 0, text)

    def test_an_inner_cfg_test_gates_the_children_of_a_module_file(self) -> None:
        files = {"src/lib.rs": "mod helpers;\n", "src/helpers.rs": "#![cfg(test)]\nmod data;\n",
                 "src/helpers/data.rs": ""}
        with CrateTree(files) as tree:
            code, text = tree.run()
        self.assertEqual(code, 1, text)
        self.assertIn("LOSAT/src/helpers/data.rs\tonly tests read it", text)
        self.assertNotIn("LOSAT/src/helpers.rs\t", text)

    def test_declarations_in_comments_and_strings_are_not_read(self) -> None:
        files = {"src/lib.rs": (
            "// mod ghost;\n/* outer /* mod ghost; */ mod ghost; */\n"
            "const S: &str = \"mod ghost;\";\nconst R: &str = r#\"mod ghost; \"quoted\"\"#;\n"
            "const C: char = '\"';\nfn f<'a>(x: &'a str) -> &'a str { x }\nmod real;\n"),
            "src/real.rs": ""}
        with CrateTree(files) as tree:
            code, text = tree.run()
        self.assertEqual(code, 0, text)

    def test_path_attributes_and_inline_modules_resolve_as_rustc_does(self) -> None:
        files = {
            "src/lib.rs": "#[path = \"imp/core_impl.rs\"]\nmod core;\nmod outer {\n    pub mod inner;\n}\nmod a;\n",
            "src/imp/core_impl.rs": "mod child;\n",  # a `#[path]` file owns its directory, as mod.rs does
            "src/imp/child.rs": "",
            "src/outer/inner.rs": "",
            "src/a.rs": "mod b;\nmod wrap {\n    mod c;\n    #[path = \"d_file.rs\"]\n    mod d;\n}\n",
            "src/a/b.rs": "",
            "src/a/wrap/c.rs": "",
            "src/a/wrap/d_file.rs": "",
        }
        with CrateTree(files) as tree:
            code, text = tree.run()
        self.assertEqual(code, 0, text)

    def test_a_file_nothing_declares_fails(self) -> None:
        files = {"src/lib.rs": "", "src/orphan.rs": ""}
        with CrateTree(files) as tree:
            code, text = tree.run()
        self.assertEqual(code, 1, text)
        self.assertIn("LOSAT/src/orphan.rs\tpackaged, but no `mod`", text)

    def test_a_declaration_without_a_file_fails(self) -> None:
        with CrateTree({"src/lib.rs": "mod missing;\n"}) as tree:
            code, text = tree.run()
        self.assertEqual(code, 1, text)
        self.assertIn("LOSAT/src/lib.rs:1: mod missing", text)

    def test_list_prints_the_packaged_files(self) -> None:
        files = {"src/lib.rs": "#[cfg(test)]\nmod tests;\n", "src/tests.rs": "", "tests/a.rs": "",
                 "README.md": "", "src/x.pyc": ""}
        with CrateTree(files, PATTERNS + ["**/*.pyc"]) as tree:
            code, text = tree.run("--list")
        self.assertEqual(code, 0, text)
        self.assertEqual(text.split(), ["Cargo.toml", "README.md", "src/lib.rs"])

    def test_a_missing_manifest_is_an_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            subprocess.run(["git", "init", "-q", tmp], check=True)
            out = io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
                code = checker.main(["--root", tmp])
        self.assertEqual(code, 2, out.getvalue())


class ExcludeTests(unittest.TestCase):
    """Cargo applies `exclude` with the gitignore rules of the `ignore` crate."""

    def excluded(self, patterns: list[str], path: str) -> bool:
        return checker.Exclude(patterns).excluding(path) is not None

    def test_a_pattern_without_a_slash_matches_at_any_depth(self) -> None:
        self.assertTrue(self.excluded(["tests/"], "tests/a.rs"))
        self.assertTrue(self.excluded(["tests/"], "src/writedb/tests/a.rs"))
        self.assertFalse(self.excluded(["tests/"], "src/tests.rs"))  # a directory pattern
        self.assertTrue(self.excluded(["*.pyc"], "src/a/b.pyc"))

    def test_a_star_does_not_cross_a_slash(self) -> None:
        self.assertTrue(self.excluded(["src/**/tests_*.rs"], "src/x/tests_a.rs"))
        self.assertFalse(self.excluded(["src/**/tests_*.rs"], "src/x/tests_dir/a.rs"))
        self.assertFalse(self.excluded(["src/*_tests.rs"], "src/x/a_tests.rs"))

    def test_a_double_star_matches_zero_or_more_directories(self) -> None:
        self.assertTrue(self.excluded(["src/**/*_tests.rs"], "src/a_tests.rs"))
        self.assertTrue(self.excluded(["src/**/*_tests.rs"], "src/x/y/a_tests.rs"))
        self.assertTrue(self.excluded(["src/**/*_tests/"], "src/x/a_tests/b/c.rs"))
        self.assertTrue(self.excluded(["**/*.pyc"], "a.pyc"))
        self.assertTrue(self.excluded(["src/**"], "src/a/b.rs"))

    def test_an_inner_or_leading_slash_anchors_the_pattern(self) -> None:
        self.assertTrue(self.excluded(["/build.rs"], "build.rs"))
        self.assertFalse(self.excluded(["/build.rs"], "src/build.rs"))
        self.assertFalse(self.excluded(["src/**/tests.rs"], "examples/src/tests.rs"))

    def test_the_last_matching_pattern_decides_and_the_file_comes_first(self) -> None:
        self.assertFalse(self.excluded(["src/", "!src/keep.rs"], "src/keep.rs"))
        self.assertTrue(self.excluded(["src/", "!src/keep.rs"], "src/other.rs"))
        self.assertTrue(self.excluded(["!a.rs", "a.rs"], "a.rs"))
        self.assertTrue(self.excluded(["data[0-9].bin"], "data7.bin"))
        self.assertFalse(self.excluded(["data[!0-9].bin"], "data7.bin"))


if __name__ == "__main__":
    unittest.main()
