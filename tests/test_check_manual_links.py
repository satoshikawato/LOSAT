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

    def test_repo_markdown_checks_links_srcset_and_forbidden_words(self):
        root = Path(tempfile.mkdtemp())
        (root / "docs" / "img").mkdir(parents=True)
        (root / "docs" / "img" / "fig1.png").write_bytes(b"png")
        (root / "docs" / "img" / "fig1_dark.png").write_bytes(b"png")
        (root / "docs" / "guide.md").write_text("# Guide\n\n## Install\n\nReady.\n", encoding="utf-8")
        readme = root / "README.md"
        readme.write_text(
            "# LOSAT\n\n"
            "<picture>\n"
            '<source media="(prefers-color-scheme: dark)" srcset="docs/img/fig1_dark.png">\n'
            '<img src="docs/img/fig1.png" alt="Figure 1">\n'
            "</picture>\n\n"
            "[Install guide](docs/guide.md#install)\n"
            "[Broken](docs/guide.md#missing)\n"
            "[MissingImg](docs/img/missing.png)\n"
            "All runs are verified here.\n",
            encoding="utf-8",
        )
        problems = "\n".join(check_manual_links.check_repo_markdown(readme, repo_root=root))
        self.assertIn("README.md:9: docs/guide.md#missing names no heading", problems)
        self.assertIn("README.md:10: docs/img/missing.png does not exist", problems)
        self.assertIn("README.md:11: says 'verified'", problems)

    def test_package_readme_forbids_relative_links_and_forbidden_words(self):
        root = Path(tempfile.mkdtemp())
        pkg = root / "python" / "README.md"
        pkg.parent.mkdir(parents=True)
        pkg.write_text(
            "# pylosat\n\n## Install\n\n"
            "[Self](#install)\n"
            "[Docs](https://satoshikawato.github.io/LOSAT/)\n"
            "[BadRel](../docs/manual/src/installation.md)\n"
            "[BadFile](bench.png)\n"
            "Verified output.\n",
            encoding="utf-8",
        )
        problems = check_manual_links.check_package_readme(pkg, display_name="python/README.md")
        joined = "\n".join(problems)
        self.assertEqual(len(problems), 3, problems)
        self.assertIn("../docs/manual/src/installation.md", joined)
        self.assertIn("bench.png", joined)
        self.assertIn("says 'Verified'", joined)


if __name__ == "__main__":
    unittest.main()
