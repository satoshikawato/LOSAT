#!/usr/bin/env python3
"""Tests for the Product Decision numbering check (`check_pd_numbers.py`).

No NCBI counterpart: the checker reads LOSAT's own decision records.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
import tempfile
import unittest
from pathlib import Path


CHECKER_PATH = Path(__file__).with_name("check_pd_numbers.py")
SPEC = importlib.util.spec_from_file_location("pd_numbers", CHECKER_PATH)
assert SPEC is not None and SPEC.loader is not None
numbers = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = numbers
SPEC.loader.exec_module(numbers)


# The form of PD-LOSAT-NCBI-DEFECTS on dev c5fd3e763: versions in parentheses inside an entry,
# an entry broken between its number and its date, `the same day`, and exception mentions
# that only refer to an exception.
GOOD = """# Product Decision: example

- Decision ID: `PD-EXAMPLE`
- Version: 1.5
- Date: 2026-10-02 (1.0); 1.1 the same day (exception 1, the first one); 1.2 2026-10-03
  (exception 2, the second one; numbered after 1.1); 1.3,
  2026-10-04 (exception 2 also for blastx, exceptions 1 and 2: first 1.2, then 1.3);
  1.4 2026-10-05 (approved exception 3, a third one, numbered after 1.3);
  1.5 2026-10-04 (exception 3: its writer; numbered after 1.4)
- Status: Accepted

## Approved exception 3 (version 1.4)

exception 3, the body may repeat this form; only the header's version list counts.
"""


def run(root: Path) -> tuple[int, str]:
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
        code = numbers.main(["--root", str(root)])
    return code, out.getvalue()


def check(text: str, name: str = "PD-EXAMPLE.md") -> tuple[int, str]:
    with tempfile.TemporaryDirectory() as tmp:
        directory = Path(tmp) / "docs" / "product_decisions"
        directory.mkdir(parents=True)
        (directory / name).write_text(text, encoding="utf-8")
        return run(Path(tmp))


def with_entry(entry: str, version: str = "1.6") -> str:
    """GOOD with one more entry after 1.5 and the Version field set to `version`."""
    text = GOOD.replace(
        "  1.5 2026-10-04 (exception 3: its writer; numbered after 1.4)\n",
        "  1.5 2026-10-04 (exception 3: its writer; numbered after 1.4);\n" + entry + "\n",
    )
    return text.replace("- Version: 1.5\n", f"- Version: {version}\n")


class PdNumberTests(unittest.TestCase):
    def test_a_well_numbered_record_passes(self) -> None:
        code, text = check(GOOD)
        self.assertEqual(code, 0, text)

    def test_entries_are_read_across_lines_and_outside_parentheses_only(self) -> None:
        lines = GOOD.splitlines()
        first, field = numbers.header_field(lines, "Date")
        entries, unbalanced = numbers.depth_zero_entries(field, first)
        self.assertIsNone(unbalanced)
        self.assertEqual([version for version, _ in entries], ["1.1", "1.2", "1.3", "1.4", "1.5"])
        self.assertEqual(entries[2], ("1.3", 6))

    def test_the_same_version_twice_fails(self) -> None:
        # Two pull requests that both took the next version, merged by keeping both entries.
        code, text = check(with_entry("  1.5 2026-10-06 (a second entry)", version="1.5"))
        self.assertEqual(code, 1, text)
        self.assertIn("version 1.5 appears twice", text)
        self.assertIn("PD-EXAMPLE.md,line=10", text)

    def test_the_same_exception_introduced_twice_fails(self) -> None:
        # 2026-11-21: #1451 and #1460 both wrote "approved exception 20".
        code, text = check(with_entry("  1.6 2026-10-06 (approved exception 3, another third one)"))
        self.assertEqual(code, 1, text)
        self.assertIn("exception 3 is introduced twice", text)
        self.assertIn("also line 8", text)

    def test_references_to_an_exception_are_not_introductions(self) -> None:
        code, text = check(with_entry("  1.6 2026-10-06 (exception 1 also for tblastn; exception 2: its reader)"))
        self.assertEqual(code, 0, text)

    def test_a_version_list_that_goes_back_fails(self) -> None:
        code, text = check(with_entry("  1.4 2026-10-06 (renumbered by hand)", version="1.4"))
        self.assertEqual(code, 1, text)
        self.assertIn("version 1.4 follows 1.5", text)

    def test_a_version_field_left_behind_fails(self) -> None:
        # The later pull request renumbered its entry to 1.6 but kept the Version field.
        code, text = check(with_entry("  1.6 2026-10-06 (exception 4, a fourth one)", version="1.5"))
        self.assertEqual(code, 1, text)
        self.assertIn("the Version field says 1.5 but the version list ends at 1.6", text)
        self.assertIn("PD-EXAMPLE.md,line=4", text)

    def test_versions_compare_as_numbers(self) -> None:
        # 1.10 follows 1.9 (as text it would sort before it).
        entries = ";\n".join(f"  1.{minor} 2026-10-{minor:02d} (entry {minor})" for minor in range(1, 11))
        record = f"# Product Decision: long\n\n- Version: 1.10\n- Date: 2026-10-01 (1.0);\n{entries}\n- Status: Accepted\n"
        code, text = check(record)
        self.assertEqual(code, 0, text)
        code, text = check(record.replace("- Version: 1.10", "- Version: 1.9"))
        self.assertEqual(code, 1, text)
        self.assertIn("the Version field says 1.9 but the version list ends at 1.10", text)

    def test_an_unbalanced_parenthesis_fails(self) -> None:
        # A merge that dropped an entry's closing parenthesis hides every entry after it.
        broken = GOOD.replace("numbered after 1.3);", "numbered after 1.3;")
        code, text = check(broken)
        self.assertEqual(code, 1, text)
        self.assertIn("unbalanced parenthesis", text)
        self.assertIn("PD-EXAMPLE.md,line=8", text)

    def test_a_record_without_a_version_list_passes(self) -> None:
        code, text = check("# Product Decision: plain\n\n- Version: 1.2\n- Status: Accepted\n")
        self.assertEqual(code, 0, text)

    def test_only_pd_files_are_read(self) -> None:
        bad = with_entry("  1.5 2026-10-06 (a second entry)", version="1.5")
        code, text = check(bad, name="NOTES.md")
        self.assertEqual(code, 0, text)


if __name__ == "__main__":
    unittest.main()
