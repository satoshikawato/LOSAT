#!/usr/bin/env python3
"""Tests for the check of source text and pinned citations in comments (`check_foreign_quotes.py`).

No NCBI counterpart: the checker reads LOSAT's own comments.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import itertools
import sys
import tempfile
import unittest
from pathlib import Path


CHECKER_PATH = Path(__file__).with_name("check_foreign_quotes.py")
SPEC = importlib.util.spec_from_file_location("foreign_quotes", CHECKER_PATH)
assert SPEC is not None and SPEC.loader is not None
checker = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = checker
SPEC.loader.exec_module(checker)


def findings(text: str, style: str = "slash") -> list[tuple[str, int, str]]:
    lines = text.splitlines()
    return [(f.kind, f.line, f.source) for f in checker.check_text("x", lines, style)]


TABLE = {"c++/src/algo/blast/core/blast_hits.c": 5000, "RpsbProc/src/rpsbproc_main.cpp": 2000}


def pinned(text: str, style: str = "slash") -> list[tuple[str, int, str]]:
    """The findings with the PIN and REF checks, against `TABLE`."""
    lines = text.splitlines()
    return [(f.kind, f.line, f.source) for f in checker.check_text("x", lines, style, TABLE)]


NCBI_QUOTE = """\
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379-1381
// ```c
//     qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*), ScoreCompareHSPs);
// ```
fn f() {}
"""
NCBI_CITATION = """\
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379-1381
// The HSPs are sorted by score (`ScoreCompareHSPs`).
fn f() {}
"""


class CheckTextTest(unittest.TestCase):
    def test_ncbi_quote_fails(self) -> None:
        self.assertEqual(findings(NCBI_QUOTE), [("QUOTE", 2, "source text")])
        self.assertEqual(findings(NCBI_CITATION), [])

    def test_a_source_language_fence_fails_without_a_citation(self) -> None:
        for language in ("c", "c++", "cpp", "asn", "cmake"):
            text = f"/// The search:\n///\n/// ```{language}\n/// x = 1;\n/// ```\nfn f() {{}}\n"
            self.assertEqual(findings(text), [("QUOTE", 3, "source text")], language)

    def test_a_text_fence_under_a_citation_fails(self) -> None:
        text = (
            "// NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:24-36\n"
            "// ```text\n// Seq-id ::= CHOICE {\n// ```\n"
        )
        self.assertEqual(findings(text), [("QUOTE", 2, "source text")])

    def test_examples_and_output_without_a_citation_pass(self) -> None:
        text = (
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379-1381\n"
            "//\n// Format:\n// ```text\n// Score = 10 bits\n// ```\n"
            "/// Example (blast_hits.c):\n/// ```rust\n/// let x = 1;\n/// ```\nfn f() {}\n"
        )
        self.assertEqual(findings(text), [])

    def test_quote_after_another_projects_label_fails(self) -> None:
        text = (
            "// Microsoft.Windows.SDK.CRTSource 10.0.22621.3 (MIT): ucrt/time/timeset.cpp:67-69\n"
            "// ```c++\n//     _timezone.initialize(8 * 3600L);\n// ```\nconst X: i32 = 0;\n"
        )
        self.assertEqual(findings(text), [("QUOTE", 2, "Microsoft C runtime")])

    def test_quote_after_an_unknown_projects_label_fails(self) -> None:
        text = "// Foo 1.2: src/foo.c:10-12\n// ```c\n// int x;\n// ```\n"
        self.assertEqual(findings(text), [("QUOTE", 2, "Foo 1.2")])

    def test_quote_after_a_prose_citation_of_libstdcxx_fails(self) -> None:
        text = (
            "    // libstdc++ reference (bits/stl_algobase.h, `std::__lower_bound`)\n"
            "    // ```c++\n    //   __len = std::distance(__first, __last);\n    // ```\n"
        )
        self.assertEqual(findings(text), [("QUOTE", 2, "libstdc++")])

    def test_ncbi_label_of_glibc_fails_twice(self) -> None:
        text = (
            "// NCBI reference (glibc 2.17 to 2.39): libio/fileops.c `_IO_new_file_overflow`\n"
            "// ```c\n//   return EOF;\n// ```\n"
        )
        self.assertEqual(findings(text), [("LABEL", 1, "glibc"), ("QUOTE", 2, "glibc")])

    def test_ncbi_label_without_a_path_fails(self) -> None:
        text = "// NCBI reference: libstdc++'s `basic_filebuf::xsgetn` (GCC 14)\n"
        self.assertEqual(findings(text), [("LABEL", 1, "libstdc++")])

    def test_copies_in_ncbis_tree_are_not_ncbis(self) -> None:
        lmdb = "c++/src/util/lmdb/mdb.c:6540-7066"
        self.assertEqual(
            findings(f"// NCBI reference (598d8ae6): {lmdb}\nfn f() {{}}\n"), [("LABEL", 1, "LMDB")]
        )
        copy = f"// LMDB 0.9.24 (NCBI 598d8ae6 copy): {lmdb}\n"
        self.assertEqual(findings(copy + "fn f() {}\n"), [])
        quoted = copy + "// ```c\n// rc = 0;\n// ```\n"
        self.assertEqual(findings(quoted), [("QUOTE", 2, "LMDB")])
        boost = "// Boost (NCBI 598d8ae6 copy): c++/src/algo/blast/core/boost_erf.c:252-255\n"
        self.assertEqual(findings(boost + "// ```c\n// x\n// ```\n"), [("QUOTE", 2, "Boost")])

    def test_a_citation_over_two_lines_is_read_from_its_label(self) -> None:
        text = (
            "// NCBI reference (598d8ae6): /usr/include/c++/8/bits/stl_heap.h (__make_heap,\n"
            "// __push_heap), and stl_algo.h __heap_select\n"
            "// ```c++\n// __heap_select(__first, __middle, __last, __comp) {\n// ```\n"
        )
        self.assertEqual(findings(text), [("LABEL", 1, "libstdc++"), ("QUOTE", 3, "libstdc++")])

    def test_a_quote_inside_a_function_takes_the_citation_above_the_function(self) -> None:
        text = (
            "    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6540-7066\n"
            "    fn cursor_put(&mut self) -> i32 {\n"
            "        let rc = 0;\n"
            "        // ```c\n        // \tnospill = flags & MDB_NOSPILL;\n        // ```\n"
            "        rc\n    }\n"
        )
        self.assertEqual(findings(text), [("QUOTE", 4, "LMDB")])

    def test_the_last_label_above_a_quote_decides(self) -> None:
        text = (
            "// libstdc++ (GCC 13): bits/list.tcc:614-662\n" + NCBI_QUOTE.replace("fn f() {}\n", "")
        )
        self.assertEqual(findings(text), [("QUOTE", 3, "source text")])
        prose = (
            "/// libstdc++'s insertion sort moves the second one (glibc's qsort, GCC, LLVM).\n"
            "/// NCBI reference: c++/src/algo/blast/api/setup_factory.cpp:170-184\n"
            "/// ```c\n///     Blast_Message2TSearchMessages(blast_msg.Get(), query_info, msgs);\n"
            "/// ```\n"
        )
        self.assertEqual(findings(prose), [("QUOTE", 3, "source text")])

    def test_a_fence_without_a_citation_is_not_checked(self) -> None:
        text = "/// Example:\n/// ```text\n/// losat blastn -query q.fa\n/// ```\nfn f() {}\n"
        self.assertEqual(findings(text), [])

    def test_ncbis_prose_naming_glibc_is_not_a_label(self) -> None:
        text = "// NCBI's `CArg_Double` reads glibc's hexadecimal numbers and a final `e`\n"
        self.assertEqual(findings(text), [])

    def test_an_ncbi_label_by_path_of_another_project_fails(self) -> None:
        text = (
            "// NCBI c++/src/algo/blast/core/boost_erf.c:252-255:\n"
            "// double ErfC(double z) { return ErfImpl(z, TRUE); }\n"
        )
        self.assertEqual(findings(text), [("LABEL", 1, "Boost"), ("QUOTE", 2, "Boost")])
        bare = "    // NCBI boost_erf.c:71-110:\n    // else result = z * (Y + p / q);\n"
        self.assertEqual(findings(bare), [("LABEL", 1, "Boost"), ("QUOTE", 2, "Boost")])

    def test_a_code_line_under_another_projects_citation_is_a_quote(self) -> None:
        text = "// glibc 2.39: libio/fileops.c:1310-1312\n// return EOF;\n// }\nfn f() {}\n"
        self.assertEqual(findings(text), [("QUOTE", 2, "glibc"), ("QUOTE", 3, "glibc")])

    def test_a_code_line_after_a_mention_of_another_project_is_not_a_quote(self) -> None:
        diff = "+// glibc's qsort is a stable merge sort, as slice::sort_by is\n+    v.sort_by(cmp);\n"
        self.assertEqual(findings(diff, "text"), [])

    def test_prose_passes_and_ncbis_code_lines_are_quotes(self) -> None:
        prose = (
            "// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:395,3969\n"
            "// \"uncommitted\" (no put area) until the first output;\n"
        )
        self.assertEqual(findings(prose), [])
        ncbi = (
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379\n"
            "// qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*), cmp);\n"
        )
        self.assertEqual(findings(ncbi), [("QUOTE", 2, "source text")])

    def test_a_code_line_under_a_label_that_cites_no_source_passes(self) -> None:
        text = "BLASTP also matched at n8. Page errors: zero.\n- reports pool_threads=8 and ready=true;\n"
        self.assertEqual(findings(text, "text"), [])

    def test_unpinned_ncbi_labels_and_checkout_paths_fail(self) -> None:
        self.assertEqual(pinned(NCBI_CITATION), [])
        for label in ("NCBI reference: ", "NCBI reference (2.17.0): ", "Reference: ncbi-blast/"):
            text = f"// {label}c++/src/algo/blast/core/blast_hits.c:1379\n"
            self.assertEqual([k for k, _, _ in pinned(text)], ["PIN"], label)
        local = "// see /mnt/c/Users/x/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1\n"
        self.assertEqual([k for k, _, _ in pinned(local)], ["PIN"])
        rpsbproc = (
            "// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1-9\n"
            "// NCBI references (598d8ae6):\n"
        )
        unit = "// NCBI unit test (598d8ae6): c++/src/algo/blast/core/blast_hits.c:10-20 Case\n"
        prose = "//! see that module's NCBI reference), and its reader\n"
        self.assertEqual(pinned(rpsbproc + unit + prose), [])
        self.assertEqual(pinned("NCBI reference: c++/x.c\n", "text"), [])

    def test_cited_paths_and_lines_must_be_in_the_pinned_trees(self) -> None:
        ok = (
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1-5000, 7\n"
            "// (c++/src/algo/blast/core/ and c++/src/algo/<dir>/<file>)\n"
        )
        self.assertEqual(pinned(ok), [])
        bad = (
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:4990-5001\n"
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hitz.c:1\n"
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:20-10\n"
            "// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/main.cpp:1\n"
        )
        self.assertEqual(
            pinned(bad),
            [
                ("REF", 1, "lines 4990-5001 outside the file's 5000"),
                ("REF", 2, "no such file"),
                ("REF", 3, "lines 20-10 outside the file's 5000"),
                ("REF", 4, "no such file"),
            ],
        )

    def test_the_table_reads_directories_and_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "t.tsv"
            path.write_text("# comment\nc++/src/a/\nb.c\t12\nc++/include/\nd.h\t3\n", encoding="utf-8")
            table = checker.load_table(path)
        self.assertEqual(table, {"c++/src/a/b.c": 12, "c++/include/d.h": 3})
        self.assertIn("c++/src/", checker.directories(table))
        self.assertIn("c++/src/a/b", checker.directories(table))

    def test_python_docstrings_and_hash_comments_are_read(self) -> None:
        docstring = (
            '"""Check the contract.\n\n'
            "NCBI reference (598d8ae6): c++/src/algo/blast/core/boost_erf.c:143-143\n"
            "```c\n    x *= y;\n```\n"
            '"""\n\nimport sys\n'
        )
        self.assertEqual(
            findings(docstring, "python"), [("LABEL", 3, "Boost"), ("QUOTE", 4, "Boost")]
        )
        toml = "# NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:395\n# ```c\n# x\n# ```\n"
        self.assertEqual(findings(toml, "hash"), [("LABEL", 1, "LMDB"), ("QUOTE", 2, "LMDB")])
        self.assertEqual(findings("#!/bin/sh\n# ```\n# glibc 2.39: libio/a.c:1\n", "hash"), [])

    def test_python_strings_that_are_not_statements_are_not_comments(self) -> None:
        text = 'X = """\\\n// ```c\nint x;\n"""\nY = "# ```c"\n\n\ndef f():\n    """Doc.\n\n    ```c\n    x;\n    ```\n    """\n'
        self.assertEqual(findings(text, "python"), [("QUOTE", 11, "source text")])

    def test_markdown_reads_only_the_paragraph_above_a_fence(self) -> None:
        text = (
            "The glibc 2.39 reference: libio/fileops.c:1310-1312.\n\n"
            "Run it:\n\n```\nlosat blastn -query q.fa\n```\n\n"
            "glibc 2.39: libio/fileops.c:1310-1312\n```c\nreturn EOF;\n```\n"
        )
        self.assertEqual(findings(text, "text"), [("QUOTE", 10, "glibc")])
        long_paragraph = "glibc 2.39: libio/fileops.c:1\n" + "ACGT" * 50 + "\n" * 1 + "x\n" * 5000
        self.assertEqual(findings(long_paragraph, "text"), [])


class StatementTest(unittest.TestCase):
    CITE = "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2-9\n"

    def quoted(self, text: str, style: str = "slash") -> list[int]:
        return [f.line for f in checker.check_text("x", text.split("\n"), style, TABLE) if f.kind == "QUOTE"]

    def test_c_statements_under_a_citation_are_quoted(self) -> None:
        for line in (
            "if (s_off + ext_to > s_range) continue;",
            "q++;",
            "} else {",
            "case eBlastp:",
            "default:",
            "Blast_HSPInit(..., &editScript, &new_hsp);",
            'SFormatSpec("length", "Alignment length", eAlignmentLength),',
        ):
            self.assertEqual(self.quoted(self.CITE + "//     " + line + "\nfn f() {}\n"), [2], line)

    def test_a_quoted_call_is_quoted_to_its_end_and_prose_is_not(self) -> None:
        text = (
            self.CITE + "//     AdjustSubjectRange(&s_start, &length,\n"
            "//                        q_start, query_length, &shift);\n"
            "// The subject range is clipped as NCBI's is.\n"
            "// Diagnostic (LOSAT_TIMING): the subject FASTA (utils/xtiming.rs);\n"
        )
        self.assertEqual(self.quoted(text), [2, 3])

    def test_prose_after_a_closed_element_is_not_quoted(self) -> None:
        text = self.CITE + '// SFormatSpec("length", "Alignment length", eAlignmentLength),\n// Keep the length the source gives, and use\n'
        self.assertEqual(self.quoted(text), [2])

    def test_code_after_a_citation_is_cut(self) -> None:
        text = "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2-3: opt.SetDust(true);\nfn f() {}\n"
        lines = text.split("\n")
        found = checker.check_text("x", lines, "slash", TABLE)
        self.assertEqual([f.kind for f in found], ["QUOTE"])
        self.assertEqual(
            checker.fix_lines(list(lines), "slash", found)[0],
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2-3",
        )
        prose = "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2-3 (blastn: 50, megablast: 6)\n"
        self.assertEqual(checker.check_text("x", [prose], "slash", TABLE), [])

    def test_c_statements_in_backticks_are_quoted_and_left_to_the_writer(self) -> None:
        for span in (
            "if (r < 1) return 0.;",
            "Int4 new_cutoff = 1;",
            "prot_seq[0] = NULLB;",
            "gap_prob = link_hsp_params->gap_prob;",
            "m_String = m_Line; m_Line = m_String;",
            "curr = stkFlexStack.top(); stkFlexStack.pop();",
            "case eUCSCRegion: return false;",
            "while (n < len && buf[n++] != 0);",
        ):
            text = f"// Uncited prose: `{span}` and more.\nfn f() {{}}\n"
            self.assertEqual(self.quoted(text), [1], span)
            lines = text.split("\n")
            found = checker.check_text("x", lines, "slash", TABLE)
            self.assertEqual(checker.fix_lines(list(lines), "slash", found), lines, span)
        # A span that goes on to the next comment line is found at its first line.
        two_lines = "// NCBI's loop (`case eBinaryASN: return\n// TestFormatBinaryAsn( mode );`).\n"
        self.assertEqual(self.quoted(two_lines), [1])
        for span in (
            "let x = 1;",
            "self.table = Vec::new();",
            "Some(x) => x;",
            "mod name;",
            "&amp;",
            ";",
            "W -t <taxid>;",
            "BlastHitSavingParametersUpdate",
            "x = y",
        ):
            self.assertEqual(self.quoted(f"// LOSAT: `{span}` here.\n"), [], span)
        # In a fence the fence checks decide; in Markdown and text files only fences are read.
        self.assertEqual(self.quoted("// ```text\n// `q++;`\n// ```\n"), [])
        self.assertEqual(self.quoted("Prose `q++;` here.\n", "text"), [])

    def expressions(
        self, text: str, rel: str = "x.rs", style: str = "slash", losat: frozenset[str] = frozenset()
    ) -> list[str]:
        """The expressions found in `text`; `losat` stands for the names of LOSAT's code."""
        return [
            f.text
            for f in checker.check_text(rel, text.split("\n"), style, TABLE, None, frozenset(), losat)
            if f.source == checker.INLINE_EXPRESSION
        ]

    def test_c_and_cxx_expressions_in_backticks_are_quoted_and_left_to_the_writer(self) -> None:
        # The Rust below reads into `text`, so `buffer` and `str` are the cited code's names.
        code = "fn read(text: &mut Vec<u8>, count: usize) { let length = count; }\n"
        for span in (
            "m_Input.GetChars(count)",
            "buffer[length] = 0",
            "ReadLength() - 1",
            "m_Input.GetStreamPosAsInt8() + ReadLengthLong(len)",
            "str.reserve(count)",
            "ERR_POST_X_ONCE(10, Warning << ...)",
            "NStr::StringToInt8(str, flags)",
            "hsp_list->best_evalue = s_BlastGetBestEvalue(hsp_list)",
            "query_blk->length + 1",
            "(Int4)opt->threshold",
            "m_DsCnt++",
            "out << value",
            "SetToOpen(to + 1)",
            "if (hrmaps.size() < dim)",
            "SeqLen_t alnLen = qseq.size()",
        ):
            text = self.CITE + f"// The reader: `{span}` first.\n" + code
            # Marks, and a statement's name foreign to the file under a citation, decide even
            # when LOSAT's code elsewhere has the names.
            self.assertEqual(self.expressions(text, losat=self.LOSAT), [f"`{span}`"], span)
            lines = text.split("\n")
            found = checker.check_text("x.rs", lines, "slash", TABLE)
            self.assertEqual(checker.fix_lines(list(lines), "slash", found), lines, span)
        # A mark of C or C++ decides without a citation, in Python too.
        self.assertEqual(self.expressions("// NCBI's `m_Pos += n` here.\n"), ["`m_Pos += n`"])
        self.assertEqual(
            self.expressions('"""NCBI calls `GetChars(count)`."""\n', "x.py", "python"), ["`GetChars(count)`"]
        )

    # Names that LOSAT's code has elsewhere in the tree.
    LOSAT = frozenset(
        "Inputs new QueryEnd Error buffer length reserve count value out threshold opt size dim qseq"
        " window diff threads from to text flags info".split()
    )

    def test_names_values_messages_and_rust_in_backticks_pass(self) -> None:
        code = "fn read(source: &mut Source, query: &[u8], subject: &[u8]) {}\n"
        for span in (
            "(m_Pos = <n>)",
            "Cannot convert string '12x' to Uint8 (m_Pos = 2)",
            "read_all(&mut source)",
            "Inputs::new(&query, &subject)",
            "QueryEnd::Error(e)",
            "a->b[c].d",
            "hsp_list->hsp_array[0]->score",
            "GetMsg()",
            "NStr::StringToInt",
            "std::string::operator==",
            "c++/src/algo/blast/core/blast_hits.c:2",
            "-evalue 1e-5",
            "LOSAT_X_BPQBLOCK=1",
            "NCBI C++ Exception:",
            "Query= ",
            "[X_STATS] TBNLUTLEAN lean= chain=",
            "0x7f",
        ):
            text = self.CITE + f"// The reader: `{span}` first.\n" + code
            self.assertEqual(self.expressions(text, losat=self.LOSAT), [], span)
        # Without a mark, arithmetic and comparisons on LOSAT's names are reasoning about values,
        # and a statement with the file's own names is LOSAT's; uncited, a statement with LOSAT's
        # names is LOSAT's too.
        code = "fn f(last_hit: i32, diag_offset: i32) { let d = last_hit; }\n"
        for span in ("last_hit <= diag_offset - window", "diff >= 0", "ceil(threads / K)", "d = last_hit"):
            found = self.expressions(self.CITE + f"// So `{span}` holds.\n" + code, losat=self.LOSAT)
            self.assertEqual(found, [], span)
        self.assertEqual(self.expressions("// LOSAT: `buffer[length] = 0` here.\n" + code, losat=self.LOSAT), [])
        # JavaScript and TypeScript are of C's family: their comments show their own code.
        self.assertEqual(self.expressions("// `new BroadcastChannel()` here.\n", "x.ts"), [])

    CALLED = frozenset("listed warning subject_offset hits".split())

    def test_text_rust_and_losats_calls_in_backticks_pass(self) -> None:
        code = (
            "fn f(mut visit: impl FnMut(i32) -> bool, w: Warning) { if let Warning::SeqId(x) = w {} }\n"
        )
        for span in (
            "'x' to Int8 (m_Pos = 0)",
            "Post-RPSBLAST Processing Utility (LOSAT 0.2.0)",
            "basic_string::replace: __pos (which is N) > this->size() (which is M)",
            "FnMut(ReaderWarning<'_>)",
            "SeqId(warning)",
            "visit(subject_offset, hits)",
            "listed(n)",
            "-KA--K",
            "long double",
            "unsigned int",
        ):
            text = self.CITE + f"// The reader: `{span}` first.\n" + code
            found = [
                f.text
                for f in checker.check_text(
                    "x.rs", text.split("\n"), "slash", TABLE, None, frozenset({"listed"}), self.CALLED
                )
                if f.source == checker.INLINE_EXPRESSION
            ]
            self.assertEqual(found, [], span)
        # A cast written as a call, and a call in camel case without arguments, are C++'s.
        for span in ("char(chain)", "ReadLength() - 1"):
            self.assertEqual(self.expressions(self.CITE + f"// So `{span}` here.\n" + code), [f"`{span}`"])

    def test_expressions_with_names_losats_code_never_uses_are_quoted(self) -> None:
        code = "fn f(text: &str, buf: &mut [u8]) {}\n"
        losat = frozenset("x y n from time text flags p i MDB_VALID std CSymDustMasker seq to start len min".split())
        for span in (
            # no mark: a name that LOSAT's code nowhere uses, cited or not
            "*fff != '\\0'",
            "rem == 4",
            "lineCount == 1",
            "getcwd(buf, 4095)",
            "maxStartLen - startLen + k_StartSequenceMargin",
            'std::length_error("vector::reserve")',
            "TSeqAlignment::MapSegSet(segset)",
            # marks of C or C++, whatever the names
            "( (x) < (y) ? -1 : (x) > (y) )",
            "((n) + 1U) & -2",
            "(Uint4)-1",
            "from--",
            "time(nullptr)",
            "CSeqDBVol::m_TaxCache(256)",
            "CFastaDeflineReader::x_ProcessIDs(text, info)",
            "NStr::StringToNumeric<Int8>(text, fConvErr_NoThrow)",
            "flags^fGenerateLocalIDs",
            "NODEPTR(p, i)",
            "(0xffff & ~(MDB_VALID))",
            "CSymDustMasker::operator()(seq)",
            "i == min<int>(to, start+len)",
        ):
            for text in (self.CITE + f"// NCBI: `{span}` first.\n" + code, f"// So `{span}` first.\n" + code):
                self.assertEqual(self.expressions(text, losat=losat), [f"`{span}`"], span)
        # Python names its own things: there only marks and cited statements decide.
        self.assertEqual(self.expressions('"""Run `losat_output_bytes(freeze)`."""\n', "x.py", "python"), [])

    def test_formulas_definitions_and_other_text_in_backticks_pass(self) -> None:
        for span in (
            "Seq-descr ::= SET OF Seqdesc",
            "|RN(x) - RN(y)| <= |x - y| + ulp",
            "F[i] = max(max_{j<i} (v[j] - g), carry)",
            "sum_i matrix[w[i]][c_i] >= threshold",
            "[a_j - d, b_j - d)",
            "num / den ~ significand * 2^(exponent - 63)",
            "--AB-CDE",
            "lcl|UniRef50_<10>",
            "unlink(2)",
            "DATA_LOADERS = none",
            "set LOSAT_X_BPSUBJRUN=0 && ...",
            "map<CSeq_id_Handle, CConstRef<CSeq_loc>>",
            "core::hint::select_unpredictable(a < b, b, a)",
        ):
            self.assertEqual(self.expressions(self.CITE + f"// So `{span}` here.\n"), [], span)

    def test_code_lines_under_a_label_with_the_path_in_backticks_are_quoted(self) -> None:
        text = (
            "    // NCBI `c++/src/algo/blast/core/blast_hits.c` lines 918-922:\n"
            "    //   if( num > 1 ) {\n"
            "    //   ++t0;\n"
            "    let x = 1;\n"
        )
        self.assertEqual(self.quoted(text), [2, 3])

    def test_an_indented_block_after_a_blank_docstring_line_is_quoted(self) -> None:
        text = (
            'def f():\n    """Chunks.\n\n'
            "    NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2-9\n\n"
            "        if (num_chunks <= 1) {\n"
            '           *chunk_size = length;"""\n    return 1\n'
        )
        self.assertEqual(self.quoted(text, "python"), [6, 7])
        lines = text.split("\n")
        found = checker.check_text("x", lines, "python", TABLE)
        fixed = "\n".join(checker.fix_lines(list(lines), "python", found))
        self.assertNotIn("chunk", fixed.split("blast_hits.c")[1])
        self.assertIn('"""\n    return 1', fixed)
        self.assertEqual(checker.python_code(fixed), checker.python_code(text))


class SettingSpanTest(unittest.TestCase):
    """`SETTING_SPAN_RE` reads the settings that `setting` below spells out, in time linear in the span."""

    @staticmethod
    def setting(text: str) -> bool:
        """A name, then `=` and a value for each `=`; a value holds no space or `=`, and spaces after it
        come before the next name (letters, digits, `_`, `.`, `[`, `]`, `-`) or at the end."""

        def name(part: str) -> bool:
            return all(c.isalnum() or c in "_.[]-" for c in part)

        parts = text.split("=")
        if len(parts) < 2 or not name(parts[0]):
            return False
        for part in parts[1:-1]:
            spaces = next((i for i, c in enumerate(part) if c.isspace()), None)
            if spaces is not None and not name(part[spaces:].lstrip()):
                return False
        return not any(c.isspace() for c in parts[-1].rstrip())

    def test_the_settings_and_nothing_else(self) -> None:
        letters = ("a", "-", "+", "=", " ")
        for length in range(7):
            for parts in itertools.product(letters, repeat=length):
                text = "".join(parts)
                self.assertEqual(bool(checker.SETTING_SPAN_RE.match(text)), self.setting(text), repr(text))
        for text, expected in (
            ("-evalue=10", True),
            ("a=1 b=2", True),
            ("x[0]=y.z", True),
            ("a=b=c", True),
            ("=", True),
            ("a=1 +b=2", False),
            ("a b=1", False),
            ("a=1 b", False),
        ):
            self.assertEqual(bool(checker.SETTING_SPAN_RE.match(text)), expected, text)
            self.assertEqual(self.setting(text), expected, text)

    def test_a_long_span_that_does_not_match_ends(self) -> None:
        # A form in which a value and the next name could split the same letters took about 0.4 s at 20
        # repeats and four times as long for each 2 more (CodeQL py/redos, alert 49).
        self.assertIsNone(checker.SETTING_SPAN_RE.match("=" + "-=" * 5000 + " +="))


class UnbacktickedCodeTest(unittest.TestCase):
    """Unbackticked C/C++ statements and expressions in comments (`UNBACKTICKED_QUOTE`)."""

    CITE = "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2-9\n"

    def unbackticked(
        self,
        text: str,
        rel: str = "x.rs",
        style: str = "slash",
        functions: frozenset[str] = frozenset(),
    ) -> list[tuple[int, str]]:
        return [
            (f.line, f.source)
            for f in checker.check_text(rel, text.split("\n"), style, TABLE, None, functions)
            if f.kind == "QUOTE"
        ]

    def assert_flagged(self, text: str, line_no: int = 1, msg: str | None = None) -> None:
        self.assertEqual(self.unbackticked(text), [(line_no, checker.UNBACKTICKED_QUOTE)], msg)

    def test_o4c_examples_fail_cited_or_uncited(self) -> None:
        cases = (
            # stage_d_results.rs:41
            "// if (m_QuerySeq[i] == m_SubjectSeq[i]) { ++num_ident; ++num_positives; }\n",
            # tblastn/search_gapped.rs:28-29
            "// typedef struct BlastSeg { Int2 frame; Int4 offset, end, gapped_start; } BlastSeg;\n",
            "// typedef struct BlastHSP { Int4 score; ... BlastSeg query, subject; } BlastHSP;\n",
            # blastn/alignment/gapped.rs:1602-1603
            "// if (score < score_gap_col) { script = SCRIPT_GAP_IN_B; score = score_gap_col; }\n",
            "// if (score < score_gap_row) { script = SCRIPT_GAP_IN_A; score = score_gap_row; }\n",
            # blastp/args.rs:123-124
            "// if (program == eBlastp && compo_mode != eNoCompositionBasedStats &&\n",
            "//     tolower(comp_stat_string[1]) == 'u') { opt.SetUnifiedP(1); }\n",
        )
        for line in cases:
            self.assert_flagged(line, 1, line)
            found = checker.check_text("x.rs", (self.CITE + line).split("\n"), "slash", TABLE)
            self.assertEqual([f.line for f in found if f.kind == "QUOTE"], [2], line)
            # Unbackticked quotes are left to the writer by `--fix` (when not cited CODE_LINE_RE).
            lines = line.split("\n")
            found_uncited = checker.check_text("x.rs", lines, "slash", TABLE)
            self.assertEqual(checker.fix_lines(list(lines), "slash", found_uncited), lines, line)

    def test_common_unbackticked_c_statement_and_expression_forms_fail(self) -> None:
        for line in (
            "// if (x > 0) { y = 1; }\n",
            "// ++count;\n",
            "// --count;\n",
            "// if (a[i] == b[i]) break;\n",
            "// NCBI: x = f(y);\n",
            "// for (i = 0; i < n; ++i)\n",
            "// return TRUE;\n",
            "/// m_Args.push_back(arg) registers the option group.\n",
            "//! NStr::StringToInt(id) parses the identifier.\n",
            "/// new_node->leftend = (Int4) midpt + 1;\n",
            "// Remove worse HSP from list: list_node->midptr = tmp_index\n",
            "// score DESC (BLAST_CMP(hsp2->score, hsp1->score))\n",
            "// NCBI: E = MAX(E, kSmallFloat)\n",
            '// NCBI: snprintf(evalue_buf, sizeof(evalue_buf), "0.0");\n',
            "// if (CheckInternalData() != 0) each result in this batch\n",
            "// struct { Uint4 q_off; Uint4 s_off; } qs_offsets;\n",
            "// if (!(b0|b1|b2)) { H2=H2_helper->ptr; H_hsp_xsum=H2->hsp_link.xsum[index]; }\n",
            "// NCBI: if( weight_divisor == 0.0 || (sum_e /= weight_divisor) > INT4_MAX )\n",
            "// hsp_array[i] && hsp_array[i+j] && hsp_array[i]->context == hsp_array[i+j]->context\n",
            "// where query_start = query_info->contexts[context].query_offset\n",
            "let x = 1; // if (a == b) { ++c; }\n",
        ):
            self.assert_flagged(line, 1, line)

    def test_m1_to_m8_guard_holes_are_closed(self) -> None:
        # M1: C/C++ tail after pinned citation or alongside a source file / foreign name.
        m1_cases = (
            '// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2: m_Outfile << "\\n\\n";\n',
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2 - ASSERT(Blast_InitHitListIsSortedByScore(init_hitlist))\n",
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2: obj = BlastCompo_AlignmentNew(gap_align->score, rule, ...);\n",
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2: (Int4)(gap_x_dropoff*NCBIMATH_LN2/min_lambda),\n",
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2 BlastInitHitListReset(init_hitlist)\n",
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2 BlastHSPList* combined_hsp_list\n",
            "// hsp->score = score; see blast_def.h\n",
            "// for (i = 0; i < n; i++) { x.c = 1; }\n",
            "// Reference: hsp->score = score;\n",
            "// x = y; // libstdc++\n",
        )
        for line in m1_cases:
            found = self.unbackticked(line)
            self.assertEqual(len(found), 1, line)
            self.assertEqual(found[0][0], 1, line)
            self.assertIn(found[0][1], (checker.UNBACKTICKED_QUOTE, checker.CODE_AFTER_CITATION), line)

        # M2: bare `file.c:NNN: <C statement>` without `c++/` prefix.
        self.assert_flagged("// blast_aalookup.c:490: w = query + offset_list[2];\n", 1)
        self.assert_flagged("// blast_args.cpp:2726-2731: AddOptionalKey(kArgMaxTargetSequences, ...);\n", 1)

        # M3: citation followed by C tail without semicolon.
        m3_cases = (
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2 - in_q_start != tree_q_start\n",
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2 - SIGN(subject.frame)\n",
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2 - while (1)\n",
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2: cutoff_score = MIN(gap_trigger * scale_factor, cutoff_score_max)\n",
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2: kDustLevel=20, kDustWindow=64, kDustLinker=1.\n",
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2: window=kSegWindow; locut=kSegLocut; hicut=kSegHicut.\n",
            '// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2: kArgQuery("query"), kArgSubject("subject").\n',
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2: Create(eBlastp, locality).\n",
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2: options.GetMatrixName(), as typed.\n",
        )
        for line in m3_cases:
            self.assert_flagged(line, 1, line)

        # M4: continuation lines of quoted C statements/enums/ternaries under a citation.
        m4_call_close = self.CITE + "//     adjusted_subject, hsp, score_options, &align_length, sbp);\nfn f() {}\n"
        self.assertEqual([ln for ln, _ in self.unbackticked(m4_call_close)], [2])
        m4_self_ptr = self.CITE + "//     self->startFreqRatios, self->ungappedLambda);\nfn f() {}\n"
        self.assertEqual([ln for ln, _ in self.unbackticked(m4_self_ptr)], [2])
        m4_enum = self.CITE + "//     eQueryOnly,\n//     eQueryAndSubject,\n//     eQueryOnlyStrandIndifferent\n// };\nfn f() {}\n"
        self.assertEqual([ln for ln, _ in self.unbackticked(m4_enum)], [2, 3, 4, 5])
        m4_ternary = (
            self.CITE
            + "// cutoff_s = do_link_hsps ?\n"
            + "// localScalingFactor = (compo_adjust_mode != eNoCompositionBasedStats)\n"
            + "//                      ? SCALING_FACTOR : 1.0;\nfn f() {}\n"
        )
        self.assertEqual(sorted(ln for ln, _ in self.unbackticked(m4_ternary)), [2, 3, 4])
        m4_open_return = self.CITE + "// return s_HSPIsContained(hsp, query_start, node->hsp, node->leftptr,\nfn f() {}\n"
        self.assertEqual([ln for ln, _ in self.unbackticked(m4_open_return)], [2])
        m4_default = self.CITE + "//     default: ...\nfn f() {}\n"
        self.assertEqual([ln for ln, _ in self.unbackticked(m4_default)], [2])

        # M5: C conditions with `||`, `&&`, `==`, `!=` and C marks / foreign names.
        m5_cases = (
            "// Calculate cutoff_score when !gapped_calculation || matrix_only_scoring.\n",
            "// Guard ans2 with sv[0] != INT4_MIN.\n",
            "// Skip removed items (linked_to==-1000) in path.\n",
            "// Regardless of flag (cutoff[0]==0).\n",
        )
        for line in m5_cases:
            self.assert_flagged(line, 1, line)

        # M6: C if/while conditions without body or logical operators, plus C-only control forms.
        m6_cases = (
            "// NCBI line 727: if (sum>H_hsp_sum)\n",
            "// if (!combined) adopt incoming; otherwise combine\n",
            "// if (i+upset-1 < leftend) { recurse on [lend..rend] }\n",
            "// In NCBI, the while (number_of_hsps > 0) loop processes ALL HSPs.\n",
            "// if (a < b)\n",
            "// while (x)\n",
            "// if(x)\n",
            "// do {\n",
            "// goto done;\n",
            "// case eBlastTypeBlastn:\n",
            "// default:\n",
        )
        for line in m6_cases:
            self.assert_flagged(line, 1, line)

        # M7: assignments, casts, macros, stdlib calls, declarations, ternaries without `;`.
        m7_cases = (
            "// BLAST_GAP_X_DROPOFF_FINAL_NUCL = 100; `-xdrop_gap_final` replaces it\n",
            "// max_dist = MIN(GREEDY_MAX_COST, max(len1, len2) / GREEDY_MAX_COST_FRACTION + 1)\n",
            "// score = MININT at the start of the row.\n",
            "true, // check_double = TRUE\n",
            "false, // check_double = FALSE (not in two-hit block)\n",
            "// lookup[(index) >> PV_ARRAY_BTS] |= (1U << (index & PV_ARRAY_MASK))\n",
            "// Match NCBI: (Int4) ell_min (direct cast, not floor)\n",
            "// qsort(hsp_array, ..., s_QueryOffsetCompareHSPs); delete equal starts;\n",
            "// y = cond ? 1 : 2\n",
            "// unsigned int n = 0;\n",
            "// int i, j;\n",
            "// double lambda = 0.267;\n",
            "// char* buf = malloc(n);\n",
            "// class CBlastFoo : public CObject {\n",
            "// enum { eFoo = 1, eBar };\n",
            "// try {\n",
            "// #define PV_SET(a, i) ...\n",
            "// #include <stdio.h>\n",
            "// static_cast<int>(x)\n",
            "// strcmp(a, b) == 0\n",
            "// p == NULL\n",
            "// hsp->score < cutoff\n",
            "// (Uint1*) seq\n",
            "// (BlastHSP*)malloc(n)\n",
            "// memcpy(dst, src, n)\n",
        )
        for line in m7_cases:
            self.assert_flagged(line, 1, line)
        self.assertEqual(
            self.unbackticked(
                "# CBlastInput(..., GetQueryBatchSize()); GetNextSeqBatch; CLocalBlast.Run()\n",
                "x.py",
                "python",
            ),
            [(1, checker.UNBACKTICKED_QUOTE)],
        )

        # M8: backtick spans with C logical OR `||` must NOT be misclassified as Rust closures.
        m8_spans = (
            "`a->b || c->d`",
            "`hsp_list->hspcnt == 0 || hsp_list->hsp_array == NULL`",
            "`m_HitsSortOption > 0 || m_HspsSortOption > 0`",
            "`seqlocs[0]->IsNull() || seqlocs[1]->IsNull()`",
            "`m_IsBl2Seq || m_IsIterative`",
            "`itr || GetSaveLastPssm()`",
        )
        for span in m8_spans:
            self.assertEqual(
                self.unbackticked(f"// Check {span} here.\n"),
                [(1, checker.INLINE_EXPRESSION)],
                span,
            )
        # Rust closures (`|a, b|`, `|| expr`) and FASTA Seq-ids (`pir||A12345`) still pass.
        for ok_span in (
            "`|a, b| a.score > b.score`",
            "`let f = || x + 1;`",
            "`pir||A12345`",
        ):
            self.assertEqual(self.unbackticked(f"// Rust {ok_span} here.\n"), [], ok_span)

    def test_rust_pseudocode_passes_and_rust_tokens_cannot_exempt_c_lines(self) -> None:
        # Genuine Rust pseudocode and safe prose mentions pass.
        rust_ok = (
            "// if x == y { return; }\n",
            "// while i < n { i += 1; }\n",
            "// } else { x = 1; }\n",
            "// else { skip }\n",
            "// x = if a > b { a } else { b };\n",
            "// x = y;\n",
            "// score = score + 1;\n",
            "// node.left = node.right;\n",
            "// return x;\n",
            "// return;\n",
            "// break;\n",
            "// push(row);\n",
            "// run_scan(q, s);\n",
            "// a[i] == b[i] for every i\n",
            "// N = (n - 1) * (m - 1);\n",
            "// idx = (q << 5) | s;\n",
            "// Stage A comparison-only C++ API oracle uses FindGeneticCode(32).\n",
            "// calls MyStruct(a) via the tuple struct\n",
            "// unsafe { *p = 0; }\n",
        )
        for line in rust_ok:
            self.assertEqual(self.unbackticked(line), [], line)

        # Rust tokens (`self`, `u8`, `Some(`, `<n>`) must NOT exempt C/C++ code on the same line.
        c_with_rust_tokens = (
            "// hsp->score = self->score;\n",
            "// Int4 x = (u8)y;\n",
            "// if (hsp->score > 0) { /* returns Some(x) */ ++n; }\n",
            "// hsp->score = 0; /* <n> items */\n",
        )
        for line in c_with_rust_tokens:
            self.assert_flagged(line, 1, line)

    def test_block_comments_are_checked_cohesively_and_unclosed_backtick_is_resilient(self) -> None:
        self.assert_flagged("/* if (m_QuerySeq[i] == m_SubjectSeq[i]) { ++num_ident; } */\nfn f() {}\n", 1)
        self.assert_flagged("/*\n * hsp->score = score;\n */\nfn f() {}\n", 2)
        bc_backtick = "/* `m_Input.GetChars(count)` */\nfn f() {}\n"
        self.assertEqual(self.unbackticked(bc_backtick), [(1, checker.INLINE_EXPRESSION)])
        # Multi-line backtick span inside `/* ... */` is checked cohesively (not misread per-line).
        bc_multiline_ok = "/*\n * Observed label `Query= seq1\n * Length=100` in report.\n */\nfn f() {}\n"
        self.assertEqual(self.unbackticked(bc_multiline_ok), [])
        bc_multiline_bad = "/*\n * Uses `m_HitsSortOption > 0 ||\n * m_HspsSortOption > 0` here.\n */\nfn f() {}\n"
        self.assertEqual(self.unbackticked(bc_multiline_bad), [(2, checker.INLINE_EXPRESSION)])
        bc_nested_ok = "/* /* matrix row A */ A */\nconst ROW: [i32; 1] = [/*A*/ 4];\n"
        self.assertEqual(self.unbackticked(bc_nested_ok), [])
        in_string = 'let s = "/* if (a == b) { ++c; } */";\nlet r = r#"/* hsp->score = 1; */"#;\n'
        self.assertEqual(self.unbackticked(in_string), [])
        # A single unclosed backtick on line 1 must not blind subsequent lines in the same run.
        unclosed_run = (
            "// Typo with an unclosed `backtick here\n"
            "// if (a == b) { ++c; }\n"
            "// And balanced `ok` after.\n"
            "fn f() {}\n"
        )
        self.assertEqual(self.unbackticked(unclosed_run), [(2, checker.UNBACKTICKED_QUOTE)])

    def test_trailing_hash_comments_in_py_sh_toml_yml_are_checked(self) -> None:
        self.assertEqual(
            self.unbackticked("x = 1  # hsp->score = score;\n", "x.py", "python"),
            [(1, checker.UNBACKTICKED_QUOTE)],
        )
        self.assertEqual(
            self.unbackticked('url = "http://x#hsp->score = 1"\n', "x.py", "python"),
            [],
        )
        for rel in ("x.sh", "x.toml", "x.yml"):
            self.assertEqual(
                self.unbackticked("x = 1  # hsp->score = score;\n", rel, "hash"),
                [(1, checker.UNBACKTICKED_QUOTE)],
                rel,
            )
            self.assertEqual(
                self.unbackticked('s = "a # hsp->score = score;"\n', rel, "hash"),
                [],
                rel,
            )
            self.assertEqual(
                self.unbackticked("s = 'a # hsp->score = score;'\n", rel, "hash"),
                [],
                rel,
            )
            self.assertEqual(
                self.unbackticked("echo ${#arr[@]} # plain note\n", rel, "hash"),
                [],
                rel,
            )

    def test_valid_losat_comments_pass(self) -> None:
        code = "fn f(x: i32, y: i32) { let _ = (x, y); }\n"
        for line in (
            "// let x = y;\n",
            "// self.count = x + 1;\n",
            "// a = b + c\n",
            "// Length=100\n",
            "// Score = 10 bits, Identities = 8/10\n",
            "// (m_Pos = <n>)\n",
            "// Cannot convert string '12x' to Uint8 (m_Pos = 2)\n",
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:10-20\n",
            "// No NCBI counterpart: debug print switch (LOSAT_DEBUG_COORDS, read once); it stays off.\n",
            "// The reader opens nothing until a name is asked for ([`LazyTaxDb`]); it is lazy.\n",
            "// NCBI's (Int4) of infinity is INT_MIN on x86-64, and the cutoff stays 1.\n",
            "// query_blk->sequence (masked) drives the traceback score.\n",
            "// Run cargo build --release --locked for the binary.\n",
            "// Measured on dev (peaks in KiB; extra = T4 - T1;\n",
            "// - k when n <= 5 or k < 3;\n",
            "// CREATE INDEX name ON TaxidInfo(parent,taxid): no sort order.\n",
            "// T->A, G->C, C->G, A->T = ACGT\n",
            "// SeqId(warning) wraps a LOSAT warning.\n",
            "// for (a coordinate, `length`, `gaps`, `gapopen` or `qframe`; the block's other fields fetch\n",
        ):
            self.assertEqual(
                self.unbackticked(line + code, functions=frozenset({"SeqId"})),
                [],
                line,
            )
        # Python docstrings allow Python doctests and dotted constructor calls, while still
        # rejecting unbackticked C/C++ statements.
        py_ok = '"""\n>>> search = pylosat.Search("blastn")\nRandom from random.Random(20261019).\n"""\n'
        self.assertEqual(self.unbackticked(py_ok, "x.py", "python"), [])
        py_bad = '"""\nif (hsp->evalue > cutoff) { hsp_array[i] = Blast_HSPFree(hsp); }\n"""\n'
        self.assertEqual(self.unbackticked(py_bad, "x.py", "python"), [(2, checker.UNBACKTICKED_QUOTE)])
        # C-family files show their own syntax.
        self.assertEqual(self.unbackticked("// if (a[i] == b[i]) { ++n; }\n", "x.ts", "slash"), [])

    def assert_all_clean(self, items: tuple[str, ...]) -> None:
        for ok in items:
            self.assertEqual(self.unbackticked(ok), [], ok)

    def test_groups_a_to_g_and_c_only_forms_are_flagged_while_rust_passes(self) -> None:
        # Group A: quote prefix + pointer member / deref index / empty-arg call assignment.
        for line in (
            "/// NCBI: H->linked_set - TRUE if this HSP is part of a multi-HSP chain\n",
            "// NCBI: matrix[*query & 0xff][*subject]\n",
            "let total = hits.len(); // NCBI: hsp_list->hspcnt\n",
            "// NCBI BLAST: word_params->cutoffs[context].x_dropoff_init\n",
            "// NCBI: gap_trigger uses kbp_std[context]->Lambda/logK.\n",
            "// NCBI: GetOffsetArraySize() = OFFSET_ARRAY_SIZE (4096) + lookup->longest_chain\n",
            "// NCBI: ungapped_data->q_start is absolute query offset.\n",
            "// NCBI line 750-752: kbp[H->hsp->context]->Lambda, kbp[H->hsp->context]->logK\n",
        ):
            self.assert_flagged(line, 1, line)
        # Reworded backticked/prose forms pass.
        self.assert_all_clean(
            (
                "/// Whether this HSP belongs to a multi-HSP chain (NCBI `H->linked_set`).\n",
                "let hspcnt = 0; let n = (hspcnt + 5).max(1024); // NCBI: the larger of 1024 and `hspcnt + 5`\n",
            )
        )

        # Group B: Nlm_* and acronym-leading CamelCase calls.
        for line in (
            "// - BSearchContextInfo(0, query_info) -> 0 (first context)\n",
            "// ---- Nlm_StepBound(x, 400, resids_x, 1 / 0.95) * 0.95, then Nlm_AddVectors\n",
            "// Nlm_EuclideanNorm(resids_x) runs element by element as resids_x[k] is formed.\n",
        ):
            self.assert_flagged(line, 1, line)

        # Group C: ALL_CAPS macro calls (excluding LOSAT_*).
        self.assert_flagged("// SET_HOUR(value % 12) stores 0 for 12, and AM adds nothing.\n", 1)
        self.assertEqual(
            self.unbackticked("# BOOST_AUTO_TEST_CASE(testHSPListSort). Ported cases.\n", "ci.yml", "hash"),
            [(1, checker.UNBACKTICKED_QUOTE)],
        )
        self.assertEqual(self.unbackticked("// LOSAT_SOME_SWITCH(1) enables tracing.\n"), [])

        # Group D: C++ class object declarations followed by prose.
        for line in (
            "// CRandom random(base_length); ambiguous ncbi4na bases are picked once\n",
            "// CRandom random(base_length); each ambiguous 4na value consumes one GetRand().\n",
        ):
            self.assert_flagged(line, 1, line)

        # Group E: named C++ enum / struct / union with trailing `;` (Rust enum/struct without `;` pass).
        self.assert_flagged("// enum EScore { eForced, ePunct, eComma, eSpace, eNewline };\n", 1)
        self.assert_flagged("// struct Foo { int x; };\n", 1)
        self.assertEqual(self.unbackticked("// enum Foo { A, B }\n"), [])
        self.assertEqual(self.unbackticked("// struct Foo { x: i32 }\n"), [])

        # Group F: copied tail of an NCBI C block comment inside a line comment.
        self.assert_flagged("//    argument can be shared between multiple threads */\n", 1)
        self.assertEqual(self.unbackticked("// Matches web/** files.\n"), [])

        # Group G: ASN.1 definitions outside backticks (backticked ASN.1 stays allowed).
        for line in (
            "// Annotdesc ::= CHOICE { name, title, comment, pub, user User-object, ... }\n",
            "// Blast-def-line-set ::= SEQUENCE OF Blast-def-line\n",
            "// Seq-id ::= CHOICE { ..., general Dbtag, ... }\n",
            "// Dbtag ::= SEQUENCE { db VisibleString, tag Object-id }\n",
            "// Object-id ::= CHOICE { id INTEGER, str VisibleString }\n",
        ):
            self.assert_flagged(line, 1, line)
        self.assertEqual(self.unbackticked("// `Seq-descr ::= SET OF Seqdesc`\n"), [])

        # C/C++-only forms (Small ④) and standalone cited calls (②a) and empty-init for (②b).
        for line in (
            "// hsp_list->Sort();\n",
            "// this->x_Foo();\n",
            "// m_Foo->Bar\n",
            '// NCBI_THROW(CBlastException, eInvalidArgument, "bad");\n',
            "// sfree(hsp_array);\n",
            "// free(p);\n",
            "// NCBI: sfree(hsp_array);\n",
            "// NCBI: push(row);\n",
            "// mask = 0xFFFFFFFFU;\n",
            "// mask = 1UL;\n",
            "// x = 1.0f;\n",
            "// x = 1e-5f;\n",
            "// mode = kBlastDefault;\n",
            "// for (;;)\n",
            "// for (; i < n; )\n",
            "// for ( ;;) {\n",
            "// for (auto& e : v) {\n",
        ):
            self.assert_flagged(line, 1, line)

        # Uncited Rust pseudo-code, printf `%1.1f`, and non-C `for` forms pass.
        self.assert_all_clean(
            (
                "// if x == y { return; }\n",
                "// x = y;\n",
                "// push(row);\n",
                "// x += y;\n",
                "// return -1;\n",
                "// break;\n",
                "// a[i] = b[i + 1];\n",
                "// let m = 0xFFu32;\n",
                "// let x = 1.0f32;\n",
                "// let v = vec![1u32; n];\n",
                "// fn f(a: u8) -> Option<(u8, u8)> {\n",
                "// std::mem::swap(&mut a, &mut b);\n",
                "// Box<dyn Fn(u8) -> u8> f;\n",
                "// for i in 0..n { sum += v[i]; }\n",
                "// for (a, b) in pairs { use_it(a, b); }\n",
                "// for (a;\n",
                "// %1.1f of an exact tie goes to the even digit.\n",
                "// a -> b\n",
                "// a=>b\n",
                "// new() creates the table\n",
            )
        )

    def test_block_comments_unified_checks_and_line_numbers(self) -> None:
        # ②c: backslash-newline continuation inside a double-quoted string must keep line numbers exact.
        multiline_str = (
            "fn f() {\n"
            '    let s = "abc \\\n'
            'def";\n'
            "    let x = 1;\n"
            "    /*\n"
            "     * hsp->score = score;\n"
            "     */\n"
            "}\n"
        )
        self.assertEqual(checker.block_comments(multiline_str.splitlines()), [[(6, "hsp->score = score;")]])
        self.assertEqual(self.unbackticked(multiline_str), [(6, checker.UNBACKTICKED_QUOTE)])

        # ②d: C block comment with cited NCBI code (matching ncbi_frame_trace.c shape) is caught.
        c_probe = (
            "/* Comparison-only probe.\n"
            " * NCBI c++/src/algo/blast/core/blast_hits.c:2-9:\n"
            " *   subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);\n"
            " *   status = s_BlastSearchEngineOneContext(...);\n"
            " *   Uint1* sequence; Int4 length;\n"
            " */\n"
            "int main(void) { return 0; }\n"
        )
        c_found = [(f.kind, f.line) for f in checker.check_text("probe.c", c_probe.splitlines(), "slash", TABLE)]
        self.assertEqual(c_found, [("QUOTE", 3), ("QUOTE", 4), ("QUOTE", 5)])

        # ②d: PIN, REF, LABEL, fence QUOTE, and inline backtick statement work inside `/* ... */`
        # across Rust, C, C++, JS, and TS extensions, while uncited C/TS syntax respects NON_C_FAMILY.
        for rel in ("a.rs", "a.c", "a.h", "a.hh", "a.cc", "a.cpp", "a.hpp", "a.js", "a.ts"):
            pin_bc = "/* NCBI reference: c++/src/algo/blast/core/blast_hits.c:2 */\n"
            self.assertEqual(
                [f.kind for f in checker.check_text(rel, pin_bc.splitlines(), "slash", TABLE)],
                ["PIN"],
                rel,
            )
            ref_bc = "/* NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:99999 */\n"
            self.assertEqual(
                [f.kind for f in checker.check_text(rel, ref_bc.splitlines(), "slash", TABLE)],
                ["REF"],
                rel,
            )
            label_bc = "/* NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:395 */\n"
            self.assertEqual(
                [f.kind for f in checker.check_text(rel, label_bc.splitlines(), "slash")],
                ["LABEL"],
                rel,
            )
            fence_bc = "/*\n * ```c\n * int x;\n * ```\n */\n"
            self.assertEqual(
                [(f.kind, f.line) for f in checker.check_text(rel, fence_bc.splitlines(), "slash", TABLE)],
                [("QUOTE", 2)],
                rel,
            )
            stmt_bc = "/* Note: `if (r < 1) return 0.;` here. */\n"
            self.assertEqual(
                [(f.kind, f.line, f.source) for f in checker.check_text(rel, stmt_bc.splitlines(), "slash", TABLE)],
                [("QUOTE", 1, checker.INLINE_STATEMENT)],
                rel,
            )
            uncited_bc = "/* if (a == b) { ++c; } */\n"
            expected_uncited = [("QUOTE", 1)] if rel == "a.rs" else []
            self.assertEqual(
                [(f.kind, f.line) for f in checker.check_text(rel, uncited_bc.splitlines(), "slash", TABLE)],
                expected_uncited,
                rel,
            )

    def test_trailing_slash_and_hash_scanner_edge_cases(self) -> None:
        # Small ⑤: all 5 trailing `//` forms are scanned while URLs in strings/raw-strings are skipped.
        for line in (
            "x;// hsp->score = score;\n",
            "x; //hsp->score = score;\n",
            "a && // hsp->score = score;\n",
            '"s" // hsp->score = score;\n',
            "'c' // hsp->score = score;\n",
        ):
            self.assert_flagged(line, 1, line)
        self.assert_all_clean(
            (
                'let u = "https://example.com/hsp->score=1";\n',
                'let r = r#"https://example.com/hsp->score=1"#;\n',
                "let c = '/';\n",
                'let u = "https://a"; // ok\n',
            )
        )

        # Small ⑥: YAML apostrophe & doubled quote, TOML multiline strings, unparseable Python fallback.
        self.assertEqual(
            self.unbackticked("a: it's fine # hsp->score = score;\n", "x.yml", "hash"),
            [(1, checker.UNBACKTICKED_QUOTE)],
        )
        self.assertEqual(
            self.unbackticked("a: 'it''s # not a comment' # hsp->score = score;\n", "x.yml", "hash"),
            [(1, checker.UNBACKTICKED_QUOTE)],
        )
        self.assertEqual(
            self.unbackticked("a: 'it''s # hsp->score = score;'\n", "x.yml", "hash"),
            [],
        )
        toml_ml = (
            'desc = """\n'
            "# hsp->score = 1;\n"
            '""" # hsp->score = 2;\n'
            "lit = '''\n"
            "# hsp->score = 3;\n"
            "''' # hsp->score = 4;\n"
        )
        self.assertEqual(
            self.unbackticked(toml_ml, "x.toml", "hash"),
            [(3, checker.UNBACKTICKED_QUOTE), (6, checker.UNBACKTICKED_QUOTE)],
        )
        broken_py = 'def broken(:\n    x = 1  # hsp->score = score;\n    s = "a # hsp->score = 1"\n'
        self.assertEqual(
            self.unbackticked(broken_py, "x.py", "python"),
            [(2, checker.UNBACKTICKED_QUOTE)],
        )

    def test_reflow_prose_shapes_do_not_false_positive(self) -> None:
        # Small ⑦ (a)–(e): re-wrapped prose shapes must pass without breaking C true positives.
        self.assert_all_clean(
            (
                # (a) `case word:` prose at line start
                "// case sensitive: the names differ\n",
                "// case c_xml_tab: the byte is a tab\n",
                "// case 3: the query is empty and nothing is written\n",
                # (b) prose under a citation containing `return` and ending with `;`
                self.CITE + "// comparators return zero after context, endpoints, and score tie;\n",
                # (c) prose under a citation starting with a C type word
                self.CITE + "// long subject (query-group spans of about 1,000 letters, chunks of\n",
                self.CITE + "// double check (the sum) before returning\n",
                # (d) short citation tail with LOSAT variable + comparison/assignment, math formula, or backticked LHS
                "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2-9: q_length <= 11 returns the\n",
                "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2-9: frame_offsets[context+1] = offset after\n",
                "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2-9 op_cost = match_score + mismatch_penalty (both positive, so\n",
                "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2-9 E = K * searchsp *\n",
                "// `eSeqId_<type>` = `e_<Type>` (Seq_id.hpp:221-244, `eSeqId_refseq` =\n",
                # (e) library name lmdb++ at end of line
                "// c++/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp:103-110 lmdb++\n",
                "// lmdb++\n",
            )
        )
        # Bare LOSAT-variable comparison tail after citation passes when the variable is known in LOSAT.
        known_cmp = "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2-9: q_length <= 11\n"
        self.assertEqual(
            [
                f.line
                for f in checker.check_text(
                    "x.rs", known_cmp.splitlines(), "slash", TABLE, None, frozenset(), frozenset({"q_length"})
                )
                if f.kind == "QUOTE"
            ],
            [],
        )
        # True C `case` labels still fail (standalone and under a citation).
        for bad in (
            "// case eBlastp:\n",
            "// case kFoo:\n",
            "// case foo_bar:\n",
            "// case NS::kFoo:\n",
            "// case 3: return 1;\n",
            "// case MY_ENUM: break;\n",
        ):
            self.assert_flagged(bad, 1, bad)
            self.assertEqual(
                [f.line for f in checker.check_text("x.rs", (self.CITE + bad).splitlines(), "slash", TABLE) if f.kind == "QUOTE"],
                [2],
                bad,
            )

    def test_comment_names_odd_backticks_and_escaped_char_literals(self) -> None:
        # ① Comment text inside `/* ... */` or trailing comments must not enter `code_names` / `rust_names`.
        self.assert_flagged("/* ZapFooBar */\n// see ZapFooBar(alpha, beta);\n", 2, "rs block comment name")
        self.assert_flagged("let x = 1; // ZapFooBar\n// see ZapFooBar(alpha, beta);\n", 2, "rs trailing comment name")
        self.assertEqual(
            self.unbackticked("x=1  # ZapFooBar\n# see ZapFooBar(alpha, beta);\n", "x.sh", "hash"),
            [(2, checker.UNBACKTICKED_QUOTE)],
        )
        # Positive control: an actual Rust definition of `ZapFooBar` still exempts it as a LOSAT name.
        self.assertEqual(
            self.unbackticked("struct ZapFooBar;\n// see ZapFooBar(alpha, beta);\n", "x.rs", "slash"),
            [],
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "a.rs").write_text("/* ZapFooBar */\n", encoding="utf-8")
            (root / "b.rs").write_text("// see ZapFooBar(alpha, beta);\n", encoding="utf-8")
            found = checker.find(root, TABLE)
            self.assertEqual([(f.path, f.line, f.kind) for f in found], [("b.rs", 1, "QUOTE")])

        # ② Cited prose `case` lines (including short lowercase prose words) still pass, while
        # cited `case kFoo:`, `case kFoo: return `x`;`, and C parameter declarations with Doxygen
        # comments fail.
        for prose_case in (
            "// case sensitive: the names differ\n",
            "// case sensitive: yes\n",
            "// case insensitive:\n",
            "// case study: foo\n",
            "// case c_xml_tab: the byte is a tab\n",
            "// case c_xml_tab: `Bad char [0x1] in string at line 8` for\n",
            "// case 3: the query is empty and nothing is written\n",
            "// Int4 values, so the reachable Int4 row length is rounded down\n",
        ):
            self.assertEqual(
                checker.check_text("x.rs", (self.CITE + prose_case).splitlines(), "slash", TABLE),
                [],
                prose_case,
            )
        for cited_bad in (
            "// case kFoo: return `x`;\n",
            "// case done: break;\n",
            "// case x: break;\n",
            "// case std::ios::eofbit: break;\n",
            "// Int4 foo(const unsigned long x,\n",
            "//     Int4 query_length, /* the length of the query */\n",
            "//     Int4 query_length, // the length of the query\n",
            "//     Int4 query_length, //!< the length of the query\n",
            "//     Int4 hitlist_size, ///< number of hits to keep in a list\n",
            "//     Int4 query_length, /**< [in] Length of the query sequence */\n",
        ):
            self.assertEqual(
                [f.line for f in checker.check_text("x.rs", (self.CITE + cited_bad).splitlines(), "slash", TABLE) if f.kind == "QUOTE"],
                [2],
                cited_bad,
            )

        # ③ Odd backticks (both multi-line and single-line) and paired GNU-style `` `foo' `` quotes
        # must not swallow C code.
        for single_odd in (
            "// an odd ` tick then hsp_list->hspcnt = hsp_cnt; and `x` y\n",
            "// Error: `foo' and hsp_list->hspcnt = hsp_cnt; see `x`\n",
        ):
            self.assertEqual(self.unbackticked(single_odd), [(1, checker.UNBACKTICKED_QUOTE)], single_odd)
        self.assertEqual(
            self.unbackticked("// Error: `foo' and plain prose see `x`\n"),
            [],
        )
        self.assertEqual(
            self.unbackticked("// see `foo' and hsp_list->hspcnt = hsp_cnt;\n// and `bar' too\n"),
            [(1, checker.UNBACKTICKED_QUOTE)],
        )
        self.assertEqual(
            self.unbackticked("// see `foo' here\n// hsp_list->hspcnt = hsp_cnt; and `bar' too\n"),
            [(2, checker.UNBACKTICKED_QUOTE)],
        )
        self.assertEqual(
            self.unbackticked("// see `foo' here\n// and `bar' too, `don't` and `b'a'` stay intact\n"),
            [],
        )
        odd_three = (
            "// an odd ` tick here\n"
            "// then hsp_list->hspcnt = hsp_cnt; and `x\n"
            "// y` done\n"
        )
        self.assertEqual(
            self.unbackticked(odd_three),
            [(2, checker.UNBACKTICKED_QUOTE)],
        )
        odd_two = (
            "// `first\n"
            "// second` then hsp_list->hspcnt = hsp_cnt; and odd `\n"
        )
        self.assertEqual(
            self.unbackticked(odd_two),
            [(2, checker.UNBACKTICKED_QUOTE)],
        )
        odd_four = (
            "// an odd ` tick here\n"
            "// `a` note then hsp_list->hspcnt = hsp_cnt; trailing `b`\n"
            "// `x\n"
            "// y` done\n"
        )
        self.assertEqual(
            self.unbackticked(odd_four),
            [(2, checker.UNBACKTICKED_QUOTE)],
        )
        gnu_quote_prose = (
            "// NCBI prints `String' when the option is missing\n"
            "// and keeps `first\n"
            "// second` wrapped across comment lines.\n"
        )
        self.assertEqual(self.unbackticked(gnu_quote_prose), [])

        # ④ Escaped single-quote char literals `'\''` must not desynchronize `scan_slash_comments`.
        self.assert_flagged(
            "let a = ['\\'','\"']; // hsp_list->hspcnt = hsp_cnt;\n",
            1,
            "escaped single quote before double quote",
        )
        self.assert_flagged(
            "let b = ['\\\\','\\'','\\n']; // hsp_list->hspcnt = hsp_cnt;\n",
            1,
            "escaped backslash and single quote",
        )
        self.assertEqual(
            self.unbackticked("let a = ['\\'','\"']; // plain prose comment\n"),
            [],
        )

    def test_unbackticked_regexes_are_linear_time(self) -> None:
        import time

        # Measure scaling ratio t(4n) / t(n) rather than brittle fixed wall-clock caps:
        # O(n) scales ~4x whereas O(n^2) scales ~16x. Include all adversarial inputs:
        # escaped-quote repeats, semi-brace-comment repeats, unclosed-for repeats,
        # unclosed-method repeats, unclosed-struct repeats, and all R5 rules/scanners.
        builders = (
            ("plain_letters", "x.rs", "slash", lambda n: "// " + "a" * n + "\n"),
            ("nested_if_parens", "x.rs", "slash", lambda n: "// if " + "(" * n + "a == b\n"),
            ("repeated_and", "x.rs", "slash", lambda n: "// if (" + "a && " * (n // 5) + "\n"),
            ("for_semis", "x.rs", "slash", lambda n: "// for (" + ";" * n + "\n"),
            ("for_a_semi", "x.rs", "slash", lambda n: "// " + "for (a;" * (n // 7) + "\n"),
            ("for_empty_spam", "x.rs", "slash", lambda n: "// " + "for ( " * (n // 6) + ";\n"),
            ("for_a_spam", "x.rs", "slash", lambda n: "// " + "for ( a" * (n // 7) + ";\n"),
            ("escaped_quote", "x.rs", "slash", lambda n: "// " + '"a\\' * (n // 3) + "\n"),
            ("semi_brace_comment", "x.rs", "slash", lambda n: "// " + "; }/* " * (n // 6) + "\n"),
            ("unclosed_method", "x.rs", "slash", lambda n: "// " + ".A(x" * (n // 4) + "\n"),
            ("unclosed_camel", "x.rs", "slash", lambda n: "// " + "FooBar(x" * (n // 8) + "\n"),
            ("unclosed_struct", "x.rs", "slash", lambda n: "// " + "struct { " * (n // 9) + "\n"),
            ("nested_if_semi", "x.rs", "slash", lambda n: "// " + "if (;if (;" * (n // 10) + "\n"),
            ("dotted_assign", "x.rs", "slash", lambda n: "// " + "a.b" * (n // 3) + " = 1\n"),
            ("unclosed_quotes", "x.rs", "slash", lambda n: "// " + 'a"' * (n // 2) + "\n"),
            ("indexed_var", "x.rs", "slash", lambda n: "// " + "x[0]" * (n // 4) + "\n"),
            ("int4_stars", "x.rs", "slash", lambda n: "// " + "Int4 " + "*" * n + "\n"),
            ("semi_brace_end", "x.rs", "slash", lambda n: "// " + "; " * (n // 2) + "}\n"),
            ("return_spaces", "x.rs", "slash", lambda n: "// return" + " " * n + "x\n"),
            ("stream_words", "x.rs", "slash", lambda n: "// " + "Stream" * (n // 6) + "\n"),
            ("block_stars", "x.rs", "slash", lambda n: "/* " + "*" * n + " */\n"),
            ("arrow_spam", "x.rs", "slash", lambda n: "// NCBI: " + "a[0]" * (n // 4) + "->1\n"),
            ("macro_spam", "x.rs", "slash", lambda n: "// " + "A_B" * (n // 3) + "(\n"),
            ("acronym_spam", "x.rs", "slash", lambda n: "// " + "ABcD" * (n // 4) + "(\n"),
            ("cpp_decl_spam", "x.rs", "slash", lambda n: "// " + "; CAb x(" * (n // 8) + "\n"),
            ("enum_spam", "x.rs", "slash", lambda n: "// " + "enum A { " * (n // 9) + "\n"),
            ("asn1_spam", "x.rs", "slash", lambda n: "// " + "A-B " * (n // 4) + "::= 1\n"),
            ("float_spam", "x.rs", "slash", lambda n: "// x = " + "1" * n + ".0\n"),
            ("int_suffix_spam", "x.rs", "slash", lambda n: "// x = " + "1" * n + "\n"),
            ("two_prose_spam", "x.rs", "slash", lambda n: "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2: x = " + "ab " * (n // 3) + "\n"),
            ("block_slash_spam", "x.rs", "slash", lambda n: "/* " + "a" * n + " */\nlet x = 1; //" + "b" * n + "\n"),
            ("backslash_str_spam", "x.rs", "slash", lambda n: 'let s = "' + "a\\\n" * (n // 3) + '"; /* x */\n'),
            ("toml_ml_spam", "x.toml", "hash", lambda n: 'a = """' + "x" * n + '""" # y\n'),
            ("yaml_quote_spam", "x.yml", "hash", lambda n: "a: " + "it's " * (n // 5) + "# y\n"),
            ("odd_bt_spam", "x.rs", "slash", lambda n: "// odd `\n// " + "a `b` " * (n // 6) + "`x\n// y` done\n"),
            ("escaped_char_spam", "x.rs", "slash", lambda n: "let a = [" + "'\\''," * (n // 5) + "'a']; // ok\n"),
            ("c_decl_comment_spam", "x.rs", "slash", lambda n: self.CITE + "// Int4 values, " + "/* " * (n // 3) + "a b c d\n"),
            ("case_bt_spam", "x.rs", "slash", lambda n: "// case kFoo: " + "`a` " * (n // 4) + "ok\n"),
            ("gnu_quote_spam", "x.rs", "slash", lambda n: "// see " + "`a' " * (n // 4) + "ok\n"),
        )
        # Use process CPU time and an 8x size jump (linear ~8x vs quadratic ~64x) so host load
        # preemption spikes never flake the bound.
        n_small = 2000
        n_large = 16000

        def measure(rel: str, style: str, lines: list[str], reps: int = 3) -> float:
            best = float("inf")
            for _ in range(reps):
                t0 = time.process_time()
                checker.check_text(rel, lines, style, TABLE)
                dt = time.process_time() - t0
                if dt < best:
                    best = dt
            return best

        for name, rel, style, build in builders:
            small_lines = build(n_small).splitlines()
            large_lines = build(n_large).splitlines()
            t_small = measure(rel, style, small_lines)
            t_large = measure(rel, style, large_lines)
            # Only enforce ratio when small time is above noise floor (0.5 ms); otherwise confirm
            # large time stays tiny (< 50 ms), which still rejects any O(n^2) blowup.
            if t_small >= 0.0005:
                ratio = t_large / t_small
                self.assertLess(
                    ratio,
                    24.0,
                    f"{name}: t({n_large})/t({n_small}) = {ratio:.2f} "
                    f"({t_large*1000:.2f}ms / {t_small*1000:.2f}ms)",
                )
            else:
                self.assertLess(
                    t_large,
                    0.050,
                    f"{name}: t({n_large}) = {t_large*1000:.2f}ms (t({n_small})={t_small*1000:.3f}ms)",
                )


class FixTest(unittest.TestCase):
    def fixed(self, text: str, style: str = "slash") -> str:
        lines = text.split("\n")
        found = checker.check_text("x", lines, style, TABLE)
        return "\n".join(checker.fix_lines(list(lines), style, found))

    def test_fix_removes_quotes_and_pins_labels_and_is_idempotent(self) -> None:
        text = (
            "    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1379-1381\n"
            "    // ```c\n    // qsort(a, n, sizeof(BlastHSP*), cmp);\n    // ```\n"
            "    //\n    // The sort is stable.\n"
            "    // NCBI reference (2.17.0): c++/src/algo/blast/core/blast_hits.c:20:\n    // x = y;\n"
            "    let x = 1;\n"
        )
        want = (
            "    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379-1381\n"
            "    //\n    // The sort is stable.\n"
            "    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:20:\n"
            "    let x = 1;\n"
        )
        once = self.fixed(text)
        self.assertEqual(once, want)
        self.assertEqual(self.fixed(once), once)
        self.assertEqual(pinned(once), [])

    def test_fix_writes_the_path_of_a_base_name_and_keeps_code(self) -> None:
        text = (
            "// NCBI reference: blast_hits.c:2-3 `Blast_HSPListSortByScore`\n"
            "// ```c\n// qsort(a, n, s, cmp);\nfn f() {}\n"
        )
        lines = text.split("\n")
        found = checker.check_text("x", lines, "slash", TABLE)
        fixed = checker.fix_lines(list(lines), "slash", found, checker.base_names(TABLE))
        self.assertEqual(
            "\n".join(fixed),
            "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2-3"
            " `Blast_HSPListSortByScore`\nfn f() {}\n",
        )

    def test_a_comment_after_code_is_pinned_and_checked(self) -> None:
        text = "    pub x: i32, // NCBI reference: blast_hits.c:2 (blastn: 50)\n"
        self.assertEqual([f.kind for f in checker.check_text("x", [text], "slash", TABLE)], ["PIN"])
        fixed = checker.fix_lines([text], "slash", checker.check_text("x", [text], "slash", TABLE), checker.base_names(TABLE))
        self.assertEqual(
            fixed, ["    pub x: i32, // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2 (blastn: 50)\n"]
        )
        late = "    pub x: i32, // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:6000\n"
        self.assertEqual([f.kind for f in checker.check_text("x", [late], "slash", TABLE)], ["REF"])
        self.assertEqual(checker.check_text("x", ['let u = "http://a"; // ok\n'], "slash", TABLE), [])

    def test_fix_keeps_a_docstrings_quotes(self) -> None:
        text = (
            'def f():\n    """Doc.\n'
            "    NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:10-20\n"
            '    ```c\n    int x;\n    """\n    return 1\n'
        )
        fixed = self.fixed(text, "python")
        self.assertIn('    """\n    return 1', fixed)
        self.assertNotIn("int x;", fixed)
        self.assertEqual(checker.python_code(fixed), checker.python_code(text))

    def test_fix_leaves_an_unclosed_markdown_fence(self) -> None:
        text = "NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2\n```c\nx;\n\nMore prose.\n"
        self.assertEqual(self.fixed(text, "text"), text)

    def test_fix_keeps_code_strings_and_line_endings(self) -> None:
        quote = "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2\r\n// ```c\r\n// x;\r\n// ```\r\n"
        in_string = 'let s = r#"\n// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2\n// ```c\n// x;\n// ```\n"#;\n'
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "a.rs").write_bytes((quote + "fn f() {}\r\n").encode())
            (root / "b.rs").write_bytes(in_string.encode())
            changed, refused = checker.fix(root)
            self.assertEqual((changed, refused), (1, ["b.rs"]))
            self.assertEqual(
                (root / "a.rs").read_bytes(),
                b"// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2\r\nfn f() {}\r\n",
            )
            self.assertEqual((root / "b.rs").read_bytes(), in_string.encode())

    def test_fix_pins_trailing_comment_without_touching_string_literal(self) -> None:
        # ②e: `--fix` must touch only the comment portion of a line, never string literals before it.
        names = checker.base_names(TABLE)
        lit = '"NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2"'
        unpinned = "NCBI reference: blast_hits.c:2"
        pinned_ref = "NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2"
        cases = (
            ("x.sh", "hash", f"msg={lit}", "#"),
            ("x.toml", "hash", f"toml_msg = {lit}", "#"),
            ("x.yml", "hash", f"msg: {lit}", "#"),
            ("x.py", "python", f"py_msg = {lit}", "#"),
            ("x.rs", "slash", f"let msg = {lit};", "//"),
        )
        for rel, style, lhs, marker in cases:
            line = f"{lhs}  {marker} {unpinned}"
            want = f"{lhs}  {marker} {pinned_ref}"
            suffix = Path(rel).suffix
            found = checker.check_text(rel, [line], style, TABLE)
            self.assertEqual([f.kind for f in found], ["PIN"], rel)
            fixed = checker.fix_lines([line], style, found, names, suffix=suffix)
            self.assertEqual(fixed, [want], rel)
            self.assertEqual(checker.check_text(rel, fixed, style, TABLE), [], rel)

    def test_fix_drops_a_blank_line_left_before_code(self) -> None:
        text = "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2\n//\n// ```c\n// x;\n// ```\nfn f() {}\n"
        self.assertEqual(
            self.fixed(text), "// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2\nfn f() {}\n"
        )


class MainTest(unittest.TestCase):
    def run_main(self, files: dict[str, str]) -> tuple[int, str]:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for rel, text in files.items():
                path = root / rel
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text, encoding="utf-8")
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                code = checker.main(["--root", str(root)])
            return code, out.getvalue()

    def test_clean_tree_passes(self) -> None:
        code, out = self.run_main({"LOSAT/src/a.rs": NCBI_CITATION})
        self.assertEqual(code, 0)
        self.assertIn("findings=0", out)

    def test_the_checked_in_table_holds_the_pinned_trees(self) -> None:
        table = checker.load_table()
        self.assertGreater(len(table), 6000)
        self.assertEqual(table["c++/src/algo/blast/core/blast_hits.c"], 3891)
        self.assertIn("RpsbProc/src/rpsbproc_main.cpp", table)

    def test_foreign_quote_fails_with_its_place(self) -> None:
        bad = "// SQLite 3.51.1: src/btree.c:404-411\n// ```c\n// x\n// ```\n"
        code, out = self.run_main({"LOSAT/src/a.rs": NCBI_CITATION, "LOSAT/src/b.rs": bad})
        self.assertEqual(code, 1)
        self.assertIn("QUOTE\tLOSAT/src/b.rs:2\tSQLite", out)
        self.assertIn("findings=1", out)

    def test_other_source_files_and_markdown_are_read(self) -> None:
        bad = "# NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:395\n# ```c\n# x\n# ```\n"
        code, out = self.run_main(
            {"LOSAT/Cargo.toml": bad, "ci.yml": bad, "README.md": bad.replace("# ", "")}
        )
        self.assertEqual(code, 1)
        self.assertIn("QUOTE\tLOSAT/Cargo.toml:2\tLMDB", out)
        self.assertIn("QUOTE\tci.yml:2\tLMDB", out)
        self.assertIn("QUOTE\tREADME.md:2\tLMDB", out)
        self.assertIn("findings=6", out)

    def test_build_output_is_skipped(self) -> None:
        bad = "// UCRT reference (10.0.26100.0): ucrt/time/tzset.cpp:1-2\n// ```c++\n// x\n// ```\n"
        code, _ = self.run_main({"target/debug/build/a.rs": bad, ".worktrees/w/b.rs": bad})
        self.assertEqual(code, 0)


if __name__ == "__main__":
    unittest.main()
