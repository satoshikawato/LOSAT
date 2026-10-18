#!/usr/bin/env python3
"""Focused tests for the floating-point contract checker.

NCBI reference (598d8ae6): c++/src/algo/blast/core/boost_erf.c:143-143
```c
            result *= expl(-z * z) / z;
```
LOSAT emulates that long-double line with one fused multiply-add; the tests freeze the
checker that keeps it the only one. They do not reproduce or alter NCBI arithmetic.
"""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


CHECKER_PATH = Path(__file__).with_name("check_fp_contract.py")
SPEC = importlib.util.spec_from_file_location("fp_contract", CHECKER_PATH)
assert SPEC is not None and SPEC.loader is not None
fp = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = fp
SPEC.loader.exec_module(fp)

SPOUGE = """\
/// `a.mul_add(b, -product)` gives it exactly.
fn x_product_error(a: f64, b: f64, product: f64) -> f64 {
    a.mul_add(b, -product)
}

#[cfg(test)]
mod tests {
    #[test]
    fn x_dekker_product_error_matches_fused_multiply_add() {
        let s = "not a.mul_add(call) // in a string";
        let q = '"';
        let fused = a.mul_add(b, -product);
    }
}
"""


def write_tree(root: Path, spouge: str = SPOUGE, other: str = "", config: str = "") -> None:
    (root / "LOSAT/src/stats").mkdir(parents=True)
    (root / "LOSAT/.cargo").mkdir(parents=True)
    (root / "LOSAT/src/stats/spouge.rs").write_text(spouge, encoding="utf-8")
    (root / "LOSAT/src/other.rs").write_text(other, encoding="utf-8")
    (root / "LOSAT/.cargo/config.toml").write_text(config, encoding="utf-8")


class SourceCheckTests(unittest.TestCase):
    def check(self, **kwargs) -> list[str]:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_tree(root, **kwargs)
            return fp.check_source(root)

    def test_the_ledger_sites_pass(self) -> None:
        self.assertEqual(self.check(), [])

    def test_a_new_mul_add_fails(self) -> None:
        errors = self.check(other="fn score(x: f64) -> f64 {\n    x.mul_add(2.0, 1.0)\n}\n")
        self.assertEqual(len(errors), 1)
        self.assertIn("fn score: 1 mul_add call(s)", errors[0])

    def test_a_removed_site_fails(self) -> None:
        errors = self.check(spouge=SPOUGE.replace("    a.mul_add(b, -product)\n", "    a * b - product\n"))
        self.assertEqual(len(errors), 1)
        self.assertIn("fn x_product_error: 0 mul_add call(s)", errors[0])

    def test_comments_and_strings_do_not_count(self) -> None:
        other = "// x.mul_add(y, z)\n/* x.mul_add(y, z) */\nconst S: &str = \"x.mul_add(y, z)\";\n"
        self.assertEqual(self.check(other=other), [])

    def test_fma_intrinsics_and_fast_math_fail(self) -> None:
        for name in ("_mm256_fmadd_pd", "_mm_fnmsub_sd", "_mm512_mask_fmadd_pd", "vfmaq_f64",
                     "f64x2_relaxed_madd", "fadd_fast", "fmul_algebraic", "fmaf64"):
            with self.subTest(name=name):
                errors = self.check(other=f"fn k() {{ let r = {name}(a, b, c); }}\n")
                self.assertEqual(len(errors), 1)
                self.assertIn(name, errors[0])

    def test_raw_strings_do_not_count(self) -> None:
        other = 'const T: &str = r#"a "quoted" \\" x.mul_add(y, z)"#;\nconst U: &[u8] = br"\\";\n'
        self.assertEqual(self.check(other=other), [])

    def test_the_libm_crate_fails(self) -> None:
        errors = self.check(other="fn k(x: f64) -> f64 { libm::log(x) }\n")
        self.assertEqual(len(errors), 1)
        self.assertIn("`libm::log` calls the libm crate", errors[0])

    def test_fast_math_and_relaxed_simd_flags_fail(self) -> None:
        for flag in ('rustflags = ["-C", "llvm-args=-fp-contract=fast"]',
                     'rustflags = ["-C", "target-feature=+simd128,+relaxed-simd"]',
                     'rustflags = ["-C", "llvm-args=-enable-unsafe-fp-math"]'):
            with self.subTest(flag=flag):
                errors = self.check(config=f"[build]\n{flag}\n")
                self.assertEqual(len(errors), 1)
                self.assertIn("LOSAT/.cargo/config.toml:2", errors[0])

    def test_a_flag_in_a_comment_passes(self) -> None:
        self.assertEqual(self.check(config="# never set -fp-contract=fast here\n"), [])


X86_DUMP = """\
0000000000012340 <LOSAT::stats::spouge::BLAST_SpougeStoE::h0123456789abcdef>:
   12340:\tpush   %rbp
   12341:\tvfmadd213sd %xmm2,%xmm1,%xmm0
   12346:\tcall   12000 <fma@plt>
   1234b:\tret

0000000000012350 <LOSAT::algorithm::blastp::kappa::redo::h0123456789abcdef>:
   12350:\tvmulsd %xmm1,%xmm0,%xmm0
   12354:\tvaddsd %xmm2,%xmm0,%xmm0
   12358:\tvfnmadd231pd %ymm2,%ymm1,%ymm0
   1235d:\tret
"""

BASELINE_DUMP = """\
0000000000135b80 <LOSAT::stats::spouge::erf_impl_ncbi>:
  135ee4:	mov    0x668f7d(%rip),%r14        # 79ee68 <_DYNAMIC+0x2d8>
  135f3f:	call   *%r14
  135f50:	mov    0x10(%rip),%rax        # 79ee70 <_DYNAMIC+0x2e0>

000000000076d350 <compiler_builtins::math::libm_math::arch::x86::fma::fma_with_fma>:
  76d350:	vfmadd213sd %xmm2,%xmm1,%xmm0
  76d355:	ret
"""

ARM_DUMP = """\
0000000000012340 <LOSAT::stats::spouge::erfc::h0123456789abcdef>:
   12340:\tfmadd\td0, d1, d2, d3
   12344:\tfmul\td0, d0, d1
   12348:\tret
"""


class BinaryCountTests(unittest.TestCase):
    def test_x86_counts_by_function(self) -> None:
        counts = fp.count_fma(X86_DUMP.splitlines(), "x86_64")
        spouge = counts["LOSAT::stats::spouge::BLAST_SpougeStoE::h0123456789abcdef"]
        self.assertEqual((spouge.instructions, spouge.refs), (1, 1))
        kappa = counts["LOSAT::algorithm::blastp::kappa::redo::h0123456789abcdef"]
        self.assertEqual((kappa.instructions, kappa.refs), (1, 0))
        self.assertEqual(kappa.mnemonics, {"vfnmadd231pd": 1})

    def test_an_fma_outside_the_sites_fails(self) -> None:
        errors = fp.judge(fp.count_fma(X86_DUMP.splitlines(), "x86_64"))
        self.assertEqual(len(errors), 1)
        self.assertIn("LOSAT::algorithm::blastp::kappa::redo", errors[0])
        self.assertIn("1 FMA instruction(s)", errors[0])

    def test_baseline_got_call_and_the_library_fma(self) -> None:
        counts = fp.count_fma(BASELINE_DUMP.splitlines(), "x86_64", frozenset({0x79EE68}))
        self.assertEqual(counts["LOSAT::stats::spouge::erf_impl_ncbi"].refs, 1)
        self.assertEqual(counts["LOSAT::stats::spouge::erf_impl_ncbi"].instructions, 0)
        library = counts["compiler_builtins::math::libm_math::arch::x86::fma::fma_with_fma"]
        self.assertEqual((library.instructions, library.refs), (1, 0))
        self.assertEqual(fp.judge(counts), [])

    def test_a_reference_to_fma_outside_the_sites_fails(self) -> None:
        dump = BASELINE_DUMP.replace("LOSAT::stats::spouge::erf_impl_ncbi", "LOSAT::stats::karlin::lambda")
        errors = fp.judge(fp.count_fma(dump.splitlines(), "x86_64", frozenset({0x79EE68})))
        self.assertEqual(len(errors), 2)
        self.assertIn("LOSAT::stats::karlin::lambda: 1 reference(s) to fma", errors[0])
        self.assertIn("the mul_add site was not found", errors[1])

    def test_only_the_sites_pass(self) -> None:
        self.assertEqual(fp.judge(fp.count_fma(ARM_DUMP.splitlines(), "aarch64")), [])

    def test_no_site_found_fails(self) -> None:
        errors = fp.judge({})
        self.assertEqual(len(errors), 1)
        self.assertIn("the mul_add site was not found", errors[0])

    def test_plain_arithmetic_is_not_fma(self) -> None:
        dump = ("0000000000000010 <f>:\n   10:\tvfmaddsub_not_an_insn\n   14:\tvmulpd %ymm1,%ymm0,%ymm0\n"
                "   18:\tfmul\td0, d0, d1\n   1c:\tcall   20 <fmax@plt>\n   20:\tfmov\td0, d1\n"
                "   24:\tfmadd_x\td0, d1\n")
        self.assertEqual(fp.count_fma(dump.splitlines(), "x86_64"), {})
        self.assertEqual(fp.count_fma(dump.splitlines(), "aarch64"), {})

    def test_the_report_line(self) -> None:
        text = fp.report(Path("LOSAT"), fp.count_fma(ARM_DUMP.splitlines(), "aarch64"))
        self.assertTrue(text.endswith("LOSAT\tfma_instructions=1\tfma_refs=0\tfunctions=1\n"))


if __name__ == "__main__":
    unittest.main()
