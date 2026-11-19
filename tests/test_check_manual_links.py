"""Tests of scripts/check_manual_links.py (the user manual's source check)."""

import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "check_manual_links.py"
spec = importlib.util.spec_from_file_location("check_manual_links", SCRIPT)
check_manual_links = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check_manual_links)


def book(pages: dict[str, str]) -> Path:
    root = Path(tempfile.mkdtemp()) / "src"
    for name, text in pages.items():
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    return root


class CheckManualLinks(unittest.TestCase):
    def test_a_clean_book_passes(self):
        src = book(
            {
                "SUMMARY.md": "# Summary\n\n- [Intro](intro.md)\n- [blastn](programs/blastn.md)\n",
                "intro.md": "# Intro\n\nSee [the tasks](programs/blastn.md#tasks-and-defaults).\n",
                "programs/blastn.md": "# blastn\n\n## Tasks and defaults\n\nBack to [intro](../intro.md).\n",
            }
        )
        self.assertEqual(check_manual_links.check(src), [])

    def test_missing_page_unlisted_page_and_bad_fragment(self):
        src = book(
            {
                "SUMMARY.md": "- [Intro](intro.md)\n- [Gone](gone.md)\n",
                "intro.md": "# Intro\n\n[x](extra.md#nowhere) [y](missing.md)\n",
                "extra.md": "# Extra\n",
            }
        )
        problems = "\n".join(check_manual_links.check(src))
        self.assertIn("lists gone.md, which does not exist", problems)
        self.assertIn("extra.md: not listed in SUMMARY.md", problems)
        self.assertIn("extra.md#nowhere names no heading", problems)
        self.assertIn("missing.md does not exist", problems)

    def test_forbidden_words_outside_code_only(self):
        src = book(
            {
                "SUMMARY.md": "- [Intro](intro.md)\n",
                "intro.md": (
                    "# Intro\n\nThe output is verified (not really).\n\n"
                    "Run `gh attestation verify FILE`.\n\n```bash\ngh attestation verify FILE\n```\n"
                    "[link](https://example.org/verify)\n"
                ),
            }
        )
        problems = check_manual_links.check(src)
        self.assertEqual(len(problems), 1, problems)
        self.assertIn("intro.md:3: says 'verified'", problems[0])

    def test_heading_ids_follow_mdbook(self):
        self.assertEqual(check_manual_links.heading_id("Coming from `NCBI` BLAST+ 2.17.0"), "coming-from-ncbi-blast-2170")
        self.assertEqual(check_manual_links.heading_id("What LOSAT refuses"), "what-losat-refuses")


if __name__ == "__main__":
    unittest.main()
