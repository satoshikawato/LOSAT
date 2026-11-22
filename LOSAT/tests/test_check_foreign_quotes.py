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
