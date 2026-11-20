#!/usr/bin/env python3
"""Tests for scripts/third_party_licenses.py (THIRD_PARTY_LICENSES and the copies of the license files).

No NCBI counterpart: this is LOSAT's own packaging tooling. The tests cover the parts that need
neither Cargo nor the network: the choice of one license from an SPDX expression, the license-file
selection, the stale check and the copies. CI's `rust lint` job also runs the script itself with
`--check`, which regenerates the file from Cargo.lock.
"""

from __future__ import annotations

import importlib.util
import json
import re
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "third_party_licenses.py"
SPEC = importlib.util.spec_from_file_location("third_party_licenses", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
tpl = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = tpl
SPEC.loader.exec_module(tpl)

MIT_TEXT = (
    "Copyright (c) 2020 Someone\n\nPermission is hereby granted, free of charge, to any person "
    "obtaining a copy of this software\n"
)
ZLIB_TEXT = (
    "This software is provided 'as-is'.\n\nPermission is granted to anyone to use this software for any "
    "purpose, including commercial applications.\n"
)


class ChooseTests(unittest.TestCase):
    def test_mit_is_preferred(self) -> None:
        for expr in ("MIT OR Apache-2.0", "Apache-2.0 OR MIT", "MIT/Apache-2.0",
                     "Apache-2.0 WITH LLVM-exception OR Apache-2.0 OR MIT", "BSD-2-Clause OR Apache-2.0 OR MIT",
                     "0BSD OR MIT OR Apache-2.0", "Unlicense OR MIT"):
            self.assertEqual(tpl.choose(expr), ["MIT"], expr)

    def test_next_most_permissive_then_apache(self) -> None:
        self.assertEqual(tpl.choose("Zlib OR Apache-2.0"), ["Zlib"])
        self.assertEqual(tpl.choose("Apache-2.0 OR BSD-3-Clause"), ["BSD-3-Clause"])
        self.assertEqual(tpl.choose("Apache-2.0 OR GPL-2.0-only"), ["Apache-2.0"])

    def test_and_takes_every_part(self) -> None:
        self.assertEqual(tpl.choose("Apache-2.0 AND MIT"), ["Apache-2.0", "MIT"])
        self.assertEqual(tpl.choose("(MIT OR Apache-2.0) AND Unicode-3.0"), ["MIT", "Unicode-3.0"])

    def test_the_more_permissive_alternative_wins(self) -> None:
        # MIT first; then what asks for nothing; then notice licenses; the LLVM exception lifts the
        # notice conditions of Apache-2.0 for object code.
        for expr, wanted in (("Apache-2.0 WITH LLVM-exception OR Apache-2.0", "Apache-2.0 WITH LLVM-exception"),
                             ("Apache-2.0 OR Apache-2.0 WITH LLVM-exception", "Apache-2.0 WITH LLVM-exception"),
                             ("Unlicense OR BSD-3-Clause", "Unlicense"), ("BSD-2-Clause OR CC0-1.0", "CC0-1.0"),
                             ("ISC OR 0BSD", "0BSD"), ("Zlib OR Unicode-3.0", "Zlib"),
                             ("Apache-2.0 WITH LLVM-exception OR BSL-1.0", "BSL-1.0")):
            self.assertEqual(tpl.choose(expr), [wanted], expr)

    def test_single_licenses(self) -> None:
        for expr in ("MIT", "Zlib", "bzip2-1.0.6", "Apache-2.0", "Unicode-3.0", "CC0-1.0"):
            self.assertEqual(tpl.choose(expr), [expr])
        self.assertEqual(tpl.choose("Apache-2.0 WITH LLVM-exception"), ["Apache-2.0 WITH LLVM-exception"])

    def test_an_unknown_or_gpl_only_license_fails(self) -> None:
        for expr in ("GPL-2.0-only", "GPL-3.0-or-later OR LGPL-2.1-only", "MIT AND Frobnicate-1.0", "MPL-2.0"):
            with self.assertRaises(tpl.LicenseError, msg=expr):
                tpl.choose(expr)

    def test_malformed_expressions_fail(self) -> None:
        for expr in ("", "MIT OR", "(MIT", "MIT)", "OR MIT", "MIT WITH"):
            with self.assertRaises(tpl.LicenseError, msg=expr):
                tpl.choose(expr)

    def test_has_choice(self) -> None:
        self.assertTrue(tpl.has_choice("MIT OR Apache-2.0"))
        self.assertTrue(tpl.has_choice("(MIT OR Apache-2.0) AND Unicode-3.0"))
        self.assertFalse(tpl.has_choice("Apache-2.0 AND MIT"))
        self.assertFalse(tpl.has_choice("MIT"))


class FileTests(unittest.TestCase):
    def test_normalize(self) -> None:
        self.assertEqual(tpl.normalize("a  \r\nb\t\r\n\r\n\r\n"), "a\nb\n")
        self.assertEqual(tpl.normalize("a\rb"), "a\nb\n")

    def test_classify(self) -> None:
        self.assertEqual(tpl.classify(MIT_TEXT), {"MIT"})
        self.assertEqual(tpl.classify(ZLIB_TEXT), {"Zlib"})
        self.assertEqual(tpl.classify("Copyright only\n"), set())
        bsd2 = "Redistribution and use in source and binary forms, with or without modification, are permitted"
        self.assertEqual(tpl.classify(bsd2), {"BSD-2-Clause"})
        self.assertEqual(tpl.classify(bsd2 + " 3. Neither the name of the University"), {"BSD-3-Clause"})

    def test_select_texts_takes_the_files_of_the_chosen_license(self) -> None:
        files = {"LICENSE-APACHE": "Apache License TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION\n",
                 "LICENSE-MIT": MIT_TEXT, "NOTICE": "notice\n", "COPYRIGHT": "pointer\n"}
        texts = tpl.select_texts("x", ["MIT"], files)
        self.assertEqual([name for name, _ in texts.licenses], ["LICENSE-MIT"])
        self.assertEqual(texts.apache, [])
        self.assertEqual([name for name, _ in texts.notices], ["NOTICE"])
        self.assertEqual([name for name, _ in texts.copyrights], ["COPYRIGHT"])
        apache = tpl.select_texts("x", ["Apache-2.0"], files)
        self.assertEqual([name for name, _ in apache.apache], ["LICENSE-APACHE"])
        self.assertEqual(apache.licenses, [])

    def test_a_missing_text_fails(self) -> None:
        with self.assertRaises(tpl.LicenseError):
            tpl.select_texts("x", ["MIT"], {"LICENSE": "MIT OR Apache-2.0\n"})
        with self.assertRaises(tpl.LicenseError):
            tpl.select_texts("x", ["Zlib"], {"LICENSE-MIT": MIT_TEXT})

    def test_a_crate_with_no_license_file_is_allowed_only_by_name(self) -> None:
        with self.assertRaises(tpl.LicenseError):
            tpl.select_texts("other", ["Apache-2.0"], {})
        self.assertEqual(tpl.select_texts("flatbuffers", ["Apache-2.0"], {}).apache, [])

    def test_two_files_with_one_text_print_once(self) -> None:
        files = {"LICENSE": MIT_TEXT, "LICENSE-MIT.md": MIT_TEXT.replace("\n", "  \n")}
        self.assertEqual(len(tpl.select_texts("x", ["MIT"], files).licenses), 1)

    def test_cut_before_drops_the_second_license(self) -> None:
        marker = tpl.CUT_BEFORE["chrono"]
        text = f"head\n{MIT_TEXT}\n{marker}\nApache stuff\n"
        texts = tpl.select_texts("chrono", ["MIT"], {"LICENSE.txt": text})
        self.assertNotIn("Apache stuff", texts.licenses[0][1])
        with self.assertRaises(tpl.LicenseError):
            tpl.select_texts("chrono", ["MIT"], {"LICENSE.txt": MIT_TEXT})

    def test_apache_body_is_cut_where_the_crate_ships_it_inside_another_file(self) -> None:
        standard = tpl.apache_standard_text()
        # The cut runs from the title line to "limitations under the License." (the end of the
        # appendix); what comes before and after stays.
        appendix = ("\nAPPENDIX: How to apply the Apache License to your work.\n\n"
                    "   See the License for the specific language governing permissions and\n"
                    "   limitations under the License.\n")
        text = (f"{MIT_TEXT}\nAs a contributor, you agree that:\n\n{'-' * 78}\n{standard}{appendix}"
                f"{'-' * 78}\nCopyright (c) 2018 Jorge Aparicio\n")
        texts = tpl.select_texts("libm", ["MIT"], {"LICENSE.txt": text})
        (_, printed), = texts.licenses
        self.assertNotIn("TERMS AND CONDITIONS", printed)
        self.assertIn(tpl.APACHE_BODY_CUT_LINE, printed)
        self.assertIn("As a contributor", printed)
        self.assertIn("Copyright (c) 2018 Jorge Aparicio", printed)
        # A crate that is not listed keeps the text whole.
        (_, whole), = tpl.select_texts("other", ["MIT"], {"LICENSE.txt": text}).licenses
        self.assertIn("TERMS AND CONDITIONS", whole)
        # A different Apache text, or none, fails.
        with self.assertRaises(tpl.LicenseError):
            tpl.select_texts("libm", ["MIT"], {"LICENSE.txt": text.replace("Grant of Copyright License", "Grant of Rights")})
        with self.assertRaises(tpl.LicenseError):
            tpl.select_texts("libm", ["MIT"], {"LICENSE.txt": MIT_TEXT})

    def test_read_bundle_drops_leading_blank_lines_only(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "bundle.txt"
            path.write_text('@@@ {"a": 1}\n\n\nfirst  \n\nsecond\n\n@@@ {"b": 2}\n', encoding="utf-8")
            self.assertEqual(tpl.read_bundle(path), [({"a": 1}, "first\n\nsecond\n"), ({"b": 2}, "")])

    def test_read_license_files(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "LICENSE-MIT").write_bytes(b"MIT text  \r\n")
            (root / "licence").write_text("x\n")
            (root / "README.md").write_text("not a license\n")
            (root / "extra.txt").write_text("named by the manifest\n")
            self.assertEqual(tpl.read_license_files(root), {"LICENSE-MIT": "MIT text\n", "licence": "x\n"})
            self.assertIn("extra.txt", tpl.read_license_files(root, "extra.txt"))
            (root / "LICENSE-BAD").write_bytes(b"\xff\xfe")
            with self.assertRaises(tpl.LicenseError):
                tpl.read_license_files(root)


class ApacheTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.standard = tpl.apache_standard_text()

    def test_the_standard_text_has_no_extras(self) -> None:
        self.assertEqual(tpl.apache_extras("x", "LICENSE", self.standard, self.standard), [])

    def test_indentation_and_wrapping_do_not_matter(self) -> None:
        indented = "".join("   " + line + "\n" if line else "\n" for line in self.standard.splitlines())
        self.assertEqual(tpl.apache_extras("x", "LICENSE", indented, self.standard), [])

    def test_own_copyright_line_and_a_trailing_list_are_printed(self) -> None:
        appendix = ("\nAPPENDIX: How to apply the Apache License to your work.\n\n   Copyright {} \n\n"
                    "   Licensed under the Apache License, Version 2.0 (the \"License\");\n"
                    "   See the License for the specific language governing permissions and\n"
                    "   limitations under the License.\n")
        boilerplate = appendix.format("[yyyy] [name of copyright owner]")
        self.assertEqual(tpl.apache_extras("x", "LICENSE", self.standard + boilerplate, self.standard), [])
        text = self.standard + appendix.format("2021 Jane Doe") + "\nfoo.rs: MIT\n"
        self.assertEqual(tpl.apache_extras("x", "LICENSE", text, self.standard), ["   Copyright 2021 Jane Doe", "foo.rs: MIT"])

    def test_a_different_text_fails(self) -> None:
        with self.assertRaises(tpl.LicenseError):
            tpl.apache_extras("x", "LICENSE", self.standard.replace("Grant of Copyright License", "Grant of Rights"), self.standard)
        with self.assertRaises(tpl.LicenseError):
            tpl.apache_extras("x", "LICENSE", "not apache\n", self.standard)


class TreeTests(unittest.TestCase):
    def test_parse_tree(self) -> None:
        out = ("losat v0.2.0 (/some/path/LOSAT)\nclap v4.6.7\nclap_builder v4.6.7\nclap v4.6.7 (*)\n"
               "pylosat v0.2.0 (/p)\n\n")
        self.assertEqual(tpl.parse_tree(out), {("clap", "4.6.7"), ("clap_builder", "4.6.7")})

    def test_an_unknown_path_package_or_a_strange_line_fails(self) -> None:
        with self.assertRaises(tpl.LicenseError):
            tpl.parse_tree("mystery v1.0.0 (/somewhere)\n")
        with self.assertRaises(tpl.LicenseError):
            tpl.parse_tree("warning: something\n")

    def test_linked_into(self) -> None:
        cli = [b for b in tpl.BUILDS if b.group == tpl.CLI]
        wheel = [b for b in tpl.BUILDS if b.group == tpl.WHEEL]
        web = [b for b in tpl.BUILDS if b.group == tpl.WEB]
        self.assertEqual(tpl.linked_into(cli + wheel + web), "losat CLI (all 7 targets); pylosat wheel (all 5 targets); LOSAT Web (both reactors)")
        native = [b for b in cli if b.kind == "native"]
        self.assertEqual(tpl.linked_into(native + wheel[:1]).split("; ")[0],
                         "losat CLI (native targets: " + ", ".join(b.label for b in native) + ")")
        self.assertEqual(tpl.linked_into(web[:1]), "LOSAT Web (serial reactor)")

    def test_builds_have_one_manifest_root_each_group(self) -> None:
        self.assertEqual({b.manifest_dir for b in tpl.BUILDS if b.group == tpl.CLI}, {"LOSAT"})
        self.assertEqual({b.manifest_dir for b in tpl.BUILDS if b.group == tpl.WHEEL}, {"python"})
        self.assertEqual({b.manifest_dir for b in tpl.BUILDS if b.group == tpl.WEB}, {"web/adapter"})


class BuildsMatchTheWorkflowsTests(unittest.TestCase):
    """BUILDS is copied by hand from the workflows; a target or a feature added there without a row
    in BUILDS would ship crates that THIRD_PARTY_LICENSES does not list."""

    @staticmethod
    def flags_of(command: str) -> tuple[str, ...]:
        words = command.split()
        flags: tuple[str, ...] = ()
        if "--no-default-features" in words:
            flags += ("--no-default-features",)
        if "--features" in words:
            flags += ("--features", words[words.index("--features") + 1])
        return flags

    @classmethod
    def setUpClass(cls) -> None:
        contract = json.loads((REPO_ROOT / "docs/release/v0.2.0_rc_contract.json").read_text(encoding="utf-8"))
        cls.artifacts = [a for a in contract["artifacts"] if a["kind"] != "source"]
        cls.builds = {(b.group, b.triple, b.flags) for b in tpl.BUILDS}

    def of_group(self, group: str) -> set[tuple[str, tuple[str, ...]]]:
        return {(triple, flags) for g, triple, flags in self.builds if g == group}

    def test_the_cli_rows_are_the_release_contract(self) -> None:
        contract = set()
        for artifact in self.artifacts:
            command = artifact["build_command"]
            self.assertIn("--bin losat", command)
            contract.add((artifact["target"], self.flags_of(command)))
        self.assertEqual(self.of_group(tpl.CLI), contract)

    def test_the_cli_rows_are_what_the_release_readiness_workflow_builds(self) -> None:
        text = (REPO_ROOT / ".github/workflows/release-readiness-v020.yml").read_text(encoding="utf-8")
        native = {m for m in re.findall(r"target_triple: (\S+)", text) if m != "source" and not m.startswith("wasm32")}
        wasm = {(m[0], self.flags_of(m[1])) for m in
                re.findall(r"cargo build --release --bin losat --target (\S+)([^\n]*?) --locked", text)}
        self.assertEqual({t for t, _ in self.of_group(tpl.CLI)}, native | {t for t, _ in wasm})
        for triple, flags in wasm:
            self.assertIn((triple, flags), self.of_group(tpl.CLI))
        # The native job builds with the default features and no --target (the runner's own host).
        self.assertIn("cd LOSAT && cargo build --release --bin losat --locked", text)
        self.assertEqual({flags for triple, flags in self.of_group(tpl.CLI) if triple in native}, {()})

    def test_the_wheel_rows_are_the_runners_of_python_wheels(self) -> None:
        text = (REPO_ROOT / ".github/workflows/python-wheels.yml").read_text(encoding="utf-8")
        # The matrix is `fromJSON(<pull request> && '["ubuntu-24.04"]' || '[...all runners...]')`.
        lists = [json.loads(m) for m in re.findall(r"'(\[[^']*\])'", text.split("matrix:", 1)[1].split("runs-on:", 1)[0])]
        runners = set(max(lists, key=len))
        by_runner = {a["runner"]: a["target"] for a in self.artifacts if a["kind"] == "native"}
        self.assertEqual(runners, set(by_runner))
        self.assertEqual({triple for triple, _ in self.of_group(tpl.WHEEL)}, set(by_runner.values()))
        self.assertEqual({flags for _, flags in self.of_group(tpl.WHEEL)}, {()})
        build = re.search(r"-m maturin build[^\n]*(?:\\\n[^\n]*)*", text)
        assert build is not None
        self.assertNotIn("--features", build.group(0))
        self.assertNotIn("--no-default-features", build.group(0))
        self.assertNotIn("--target", build.group(0))

    def test_the_web_rows_are_the_reactors_of_build_reactors(self) -> None:
        path = REPO_ROOT / "web/adapter/tools/build_reactors.py"
        spec = importlib.util.spec_from_file_location("build_reactors_for_licenses", path)
        assert spec is not None and spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        reactors = {(triple, tuple(option for option in options if option != "--lib"))
                    for triple, options, _, _ in module.SPECS.values()}
        self.assertEqual(self.of_group(tpl.WEB), reactors)


class StaleTests(unittest.TestCase):
    def test_first_difference_names_the_line(self) -> None:
        self.assertIn("line 2", tpl.first_difference(b"a\nb\n", b"a\nc\n"))
        self.assertIn("line 1", tpl.first_difference(b"x", b""))

    def test_stale_reports_each_file(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "good").write_bytes(b"same\n")
            (root / "old").write_bytes(b"old\n")
            problems = tpl.stale({"good": b"same\n", "old": b"new\n", "gone": b"x\n"}, root)
            self.assertEqual(len(problems), 2)
            self.assertTrue(any(p.startswith("old: stale") and "line 1" in p for p in problems))
            self.assertIn("gone: missing", problems)

    def test_copies_follow_their_originals(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "LICENSE").write_bytes(b"license\n")
            (root / "THIRD_PARTY_NOTICES").write_bytes(b"notices\n")
            with mock.patch.object(tpl, "ROOT", root):
                files = tpl.with_copies({tpl.GENERATED: b"generated\n"})
                self.assertEqual(files["LOSAT/LICENSE"], b"license\n")
                self.assertEqual(files["python/LICENSE"], b"license\n")
                self.assertEqual(files["LOSAT/THIRD_PARTY_NOTICES"], b"notices\n")
                self.assertEqual(files["python/THIRD_PARTY_NOTICES"], b"notices\n")
                self.assertEqual(files["python/THIRD_PARTY_LICENSES"], b"generated\n")
                self.assertEqual(len(tpl.stale(files, root)), len(files))  # nothing of the result exists yet
                for relative, data in files.items():
                    (root / relative).parent.mkdir(parents=True, exist_ok=True)
                    (root / relative).write_bytes(data)
                self.assertEqual(tpl.stale(files, root), [])
                (root / "python/LICENSE").write_bytes(b"edited\n")
                self.assertEqual([p.split(":")[0] for p in tpl.stale(files, root)], ["python/LICENSE"])
                (root / "LICENSE").unlink()
                with self.assertRaises(tpl.LicenseError):
                    tpl.with_copies({tpl.GENERATED: b"g\n"})


class RepositoryTests(unittest.TestCase):
    """The committed files, without Cargo or the network."""

    def test_the_copies_of_the_repository_are_byte_identical_to_their_originals(self) -> None:
        generated = (REPO_ROOT / tpl.GENERATED).read_bytes()
        self.assertEqual(tpl.stale(tpl.with_copies({tpl.GENERATED: generated}), REPO_ROOT), [])

    def test_the_inputs_are_the_ones_the_script_records(self) -> None:
        self.assertEqual(tpl.input_problems(), [])

    def test_the_generated_file_has_no_machine_specific_text(self) -> None:
        text = (REPO_ROOT / tpl.GENERATED).read_text(encoding="utf-8")
        for forbidden in ("/home/", "/Users/", "C:\\", "\r", ".cargo/registry"):
            self.assertNotIn(forbidden, text)
        self.assertTrue(text.endswith("\n"))
        self.assertFalse(any(line != line.rstrip() for line in text.split("\n")))

    def test_every_linked_license_text_is_present(self) -> None:
        text = (REPO_ROOT / tpl.GENERATED).read_text(encoding="utf-8")
        for needle in ("1. CRATES", "2. THE RUST STANDARD LIBRARY", "3. WASI-LIBC",
                       "4. WHAT ZIG LINKS INTO THE LINUX WHEELS OF PYLOSAT", "6. APACHE LICENSE, VERSION 2.0",
                       "NOTICE 1", "TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION",
                       "libc-top-half/musl/COPYRIGHT", "libc-bottom-half/cloudlibc/LICENSE",
                       "--- lib/libunwind/LICENSE.TXT ---", "Copyright (c) Zig contributors"):
            self.assertIn(needle, text)

    def test_the_zig_input_is_the_zig_that_links_the_linux_wheels(self) -> None:
        # python-wheels.yml installs ziglang==ZIG_VERSION and links the Linux wheels with --zig.
        workflow = (REPO_ROOT / ".github" / "workflows" / "python-wheels.yml").read_text(encoding="utf-8")
        pins = sorted(set(re.findall(r"ziglang==([0-9][^\s\"']*)", workflow)))
        self.assertEqual(pins, [tpl.ZIG_VERSION])
        self.assertIn("--zig --compatibility manylinux2014", workflow)
        meta = tpl.read_bundle(tpl.INPUTS / tpl.ZIG_BUNDLE)[0][0]
        self.assertEqual(meta["version"], tpl.ZIG_VERSION)

    def test_pointers_into_the_notices_name_their_section(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rule = "=" * 79
            (root / "THIRD_PARTY_NOTICES").write_text(
                f"{rule}\n1. NCBI BLAST+\n{rule}\n\n1. not a heading\n{rule}\n2. LMDB 0.9.24\n{rule}\n", encoding="utf-8")
            tpn = "THIRD_PARTY_" + "NOTICES"  # split, so that this file holds no pointer itself
            (root / "a.rs").write_text(
                f"// {tpn}, section 2 (LMDB).\n"     # right
                f"// {tpn}, section 2.\n"            # no name
                f"// {tpn}, section 1 (LMDB).\n"     # wrong name
                f"// {tpn} (section 3 (Boost))\n",   # no such section
                encoding="utf-8")
            (root / "b.bin").write_text(f"{tpn}, section 9\n", encoding="utf-8")
            # A Vue template's HTML comment (web/app/src/ui) is read too.
            (root / "c.vue").write_text(f"<!-- {tpn}, section 2 (LMDB). -->\n<!-- {tpn}, section 1 (Boost). -->\n",
                                        encoding="utf-8")
            problems = tpl.pointer_problems(root, ["a.rs", "b.bin", "c.vue", "missing.rs"])
        self.assertEqual([":".join(p.split(":")[:2]) for p in problems], ["a.rs:2", "a.rs:3", "a.rs:4", "c.vue:2"])

    def test_the_repository_pointers_name_their_section(self) -> None:
        self.assertEqual(tpl.pointer_problems(REPO_ROOT), [])

    def test_the_toolchain_channel_is_the_one_of_the_inputs(self) -> None:
        self.assertEqual(tpl.toolchain_channel(), tpl.RUST_VERSION)


if __name__ == "__main__":
    unittest.main()
