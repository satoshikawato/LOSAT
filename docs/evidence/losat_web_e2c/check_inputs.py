#!/usr/bin/env python3
"""Compare BLASTN input handling and query batches with NCBI (Session S07+; comparison only).

Writes derived inputs into WORK (from the repository's inputs and fixed seeds), runs each
case in NCBI and LOSAT (standard input empty; `-out` files compared with stdout), and
prints one line per case: `same`, `same-error`,
`losat-rejects` (NCBI succeeds; LOSAT fails with a message that names what it does not
support, whether NCBI succeeds or fails otherwise) or `DIFF`, with the `expect` column of
the case. Exits 1 when a result is not the expected one.

Usage: check_inputs.py --bin-dir DIR --losat LOSAT --work DIR
"""
from __future__ import annotations

import argparse
import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "losat_web_e2a"))
from run_oracle import ENGINE  # noqa: E402

REJECTION_MARKERS = (b"which LOSAT does not reproduce", b"not supported by LOSAT")
F = "tests/fasta/outfmt0"
COMPACT = "tests/fasta/blastn_parity_compact.fasta"


def records(path: Path) -> list[tuple[str, str]]:
    chunks = [chunk.split("\n", 1) for chunk in path.read_text().split(">")[1:]]
    return [(head, body.replace("\n", "")) for head, body in chunks]


def make_inputs(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    multi = records(ENGINE / F / "multi_query.fasta")
    q = multi[0][1]
    (work / "tab_after_id.fa").write_text(f">t1\tsecond field\n{q}\n")
    (work / "tab_in_title.fa").write_text(f">t1 a\tb c\n{q}\n")
    (work / "leading_space.fa").write_text(f"> t1 lead space\n{q}\n")
    (work / "x_residue.fa").write_text(f">x\n{q[:100]}XX{q[100:]}\n")
    (work / "hyphen.fa").write_text(f">h\n{q[:100]}--{q[100:]}\n")
    (work / "digits.fa").write_text(f">d\n1 {q[:60]}\n61 {q[60:120]}\n")
    subject = (ENGINE / F / "multi_subject.fasta").read_text()
    (work / "utf8_subject.fa").write_text(subject.replace(">msA close homolog", ">msA ümlaut homolog", 1))
    (work / "empty.fa").write_text("")
    (work / "white_space.fa").write_text("  \n\t\n")
    compact = records(ENGINE / COMPACT)
    alpha, beta = compact[0][1], compact[1]
    fill = ">fill 5000 residues\n" + (alpha * 80)[:5000] + "\n"
    (work / "invalid_at_end.fa").write_text(fill + f">{beta[0]}\n{beta[1]}\n>n1 short N\n" + "N" * 40 + "\n")
    (work / "long_invalid_run.fa").write_text(fill + ">n1 long N\n" + "N" * 150 + f"\n>{beta[0]}\n{beta[1]}\n")
    # 30 queries of 300 residues (more than the first batch of 5000) against one subject;
    # the IUPAC copy gives the contexts different compositions.
    rng = random.Random(5)
    genome = "".join(rng.choice("ACGT") for _ in range(20000))
    (work / "big_subject.fa").write_text(">bigs\n" + genome + "\n")
    acgt, iupac = [], []
    for index in range(30):
        start = rng.randrange(0, 19000)
        seq = "".join(rng.choice("ACGT") if rng.random() < 0.05 else c for c in genome[start:start + 300])
        acgt.append(f">a{index}\n{seq}\n")
        mixed = "".join(rng.choice("RYKMSW") if rng.random() < 0.02 * (index % 5) else c for c in seq)
        iupac.append(f">i{index}\n{mixed}\n")
    (work / "many_acgt.fa").write_text("".join(acgt))
    (work / "many_iupac.fa").write_text("".join(iupac))
    # Cases of the S07+ independent audit.
    subject_lines = subject.split("\n")
    subject_lines[1] = subject_lines[1][:10] + "XX" + subject_lines[1][10:]
    (work / "subject_x.fa").write_text("\n".join(subject_lines))
    (work / "empty_defline.fa").write_text(f">\n{q}\n")
    (work / "blank_defline.fa").write_text(f">   \n{q}\n")
    (work / "empty_record.fa").write_text(f">e1 empty\n>mq1 real\n{q}\n")
    (work / "empty_records_only.fa").write_text(">e1 empty\n>e2 empty\n")
    (work / "empty_subject_record.fa").write_text(">s0 empty\n" + subject)
    (work / "leading_blank_line.fa").write_text(f"\n>q1\n{q}\n")
    (work / "byte_order_mark.fa").write_text(f"\ufeff>q1\n{q}\n")
    (work / "first_multi_query.fa").write_text(f">{multi[0][0]}\n{q}\n")
    edl933 = "".join(line for line in (ENGINE / "tests/fasta/EDL933.fna").read_text().splitlines()[1:])
    (work / "edl933_100k.fa").write_text(">edl933_100k first 100 kb of EDL933\n" + edl933[:100000] + "\n")


def cases(work: Path) -> list[tuple[str, list[str], str]]:
    w = str(work)
    multi = ["-query", f"{F}/multi_query.fasta", "-subject", f"{F}/multi_subject.fasta"]
    rows = [
        ("rna.fmt6", ["-query", f"{F}/rna_query.fasta", "-subject", f"{F}/multi_subject.fasta", "-outfmt", "6"], "same"),
        ("rna.fmt7.blastn", ["-query", f"{F}/rna_query.fasta", "-subject", f"{F}/multi_subject.fasta", "-task", "blastn", "-outfmt", "7"], "same"),
        ("rna.lcase", ["-query", f"{F}/rna_query.fasta", "-subject", f"{F}/multi_subject.fasta", "-lcase_masking"], "same"),
        ("tab_after_id", ["-query", f"{w}/tab_after_id.fa", "-subject", f"{F}/multi_subject.fasta"], "losat-rejects"),
        ("tab_in_title", ["-query", f"{w}/tab_in_title.fa", "-subject", f"{F}/multi_subject.fasta", "-outfmt", "7"], "losat-rejects"),
        ("leading_space", ["-query", f"{w}/leading_space.fa", "-subject", f"{F}/multi_subject.fasta"], "losat-rejects"),
        ("x_residue", ["-query", f"{w}/x_residue.fa", "-subject", f"{F}/multi_subject.fasta", "-outfmt", "6"], "losat-rejects"),
        ("hyphen", ["-query", f"{w}/hyphen.fa", "-subject", f"{F}/multi_subject.fasta", "-outfmt", "6"], "losat-rejects"),
        ("digits", ["-query", f"{w}/digits.fa", "-subject", f"{F}/multi_subject.fasta", "-outfmt", "6"], "losat-rejects"),
        ("utf8_subject", ["-query", f"{F}/multi_query.fasta", "-subject", f"{w}/utf8_subject.fa"], "losat-rejects"),
        ("empty_query", ["-query", f"{w}/empty.fa", "-subject", f"{F}/multi_subject.fasta"], "same"),
        ("white_space_query", ["-query", f"{w}/white_space.fa", "-subject", f"{F}/multi_subject.fasta", "-max_target_seqs", "2"], "same"),
        ("empty_query.penalty0", ["-query", f"{w}/empty.fa", "-subject", f"{F}/multi_subject.fasta", "-max_target_seqs", "2", "-penalty", "0"], "same-error"),
        ("few_matches.fmt6", [*multi, "-max_target_seqs", "2", "-outfmt", "6"], "same"),
        ("all_n.fmt7", ["-query", f"{F}/edge_allN.fasta", "-subject", COMPACT, "-outfmt", "7"], "same"),
        ("batch_all_n.fmt7", ["-query", f"{F}/edge_batch_allN.fasta", "-subject", COMPACT, "-outfmt", "7"], "same"),
        ("batch_all_n.blastn.fmt7", ["-query", f"{F}/edge_batch_allN.fasta", "-subject", COMPACT, "-task", "blastn", "-outfmt", "7"], "same"),
        ("mixed_all_n.fmt7", ["-query", f"{F}/edge_mixed_allN.fasta", "-subject", COMPACT, "-outfmt", "7"], "same"),
        ("short_invalid.fmt7", ["-query", f"{F}/edge_short_invalid.fasta", "-subject", COMPACT, "-task", "blastn", "-outfmt", "7"], "same"),
        ("invalid_at_end.fmt7", ["-query", f"{w}/invalid_at_end.fa", "-subject", COMPACT, "-outfmt", "7"], "losat-rejects"),
        ("invalid_at_end.fmt6", ["-query", f"{w}/invalid_at_end.fa", "-subject", COMPACT, "-outfmt", "6"], "same"),
        ("long_invalid_run.fmt0", ["-query", f"{w}/long_invalid_run.fa", "-subject", COMPACT], "losat-rejects"),
        ("scoring_error.first_batch_invalid", ["-query", f"{F}/edge_batch_allN.fasta", "-subject", COMPACT, "-reward", "1", "-penalty", "-6", "-outfmt", "6"], "losat-rejects"),
        ("scoring_error.all_invalid", ["-query", f"{F}/edge_allN.fasta", "-subject", COMPACT, "-reward", "1", "-penalty", "-6"], "same"),
        ("scoring_error.outfmt0_prolog", [*multi, "-reward", "3", "-penalty", "-5"], "same-error"),
        ("scoring_error.megablast_prolog", [*multi, "-task", "blastn", "-reward", "2", "-penalty", "-5", "-max_target_seqs", "3"], "same-error"),
    ]
    multi_s = ["-subject", f"{F}/multi_subject.fasta"]
    rows += [
        # An invalid query in the first batch of an unsupported scoring: NCBI crashes.
        ("audit.invalid_first_batch.table_error", ["-query", f"{F}/edge_mixed_allN.fasta", "-subject", COMPACT, "-reward", "1", "-penalty", "-6", "-outfmt", "6"], "losat-rejects"),
        ("audit.invalid_first_batch.gap_error", ["-query", f"{F}/edge_mixed_allN.fasta", "-subject", COMPACT, "-task", "blastn", "-reward", "1", "-penalty", "-2", "-gapopen", "1", "-gapextend", "3", "-outfmt", "7"], "losat-rejects"),
        ("audit.invalid_strand.3_-1", ["-query", f"{w}/first_multi_query.fa", *multi_s, "-task", "blastn", "-reward", "3", "-penalty", "-1", "-outfmt", "6"], "losat-rejects"),
        ("audit.empty_query.subject_x", ["-query", f"{w}/empty.fa", "-subject", f"{w}/subject_x.fa", "-outfmt", "6"], "losat-rejects"),
        ("audit.empty_defline", ["-query", f"{w}/empty_defline.fa", *multi_s, "-outfmt", "6"], "losat-rejects"),
        ("audit.blank_defline", ["-query", f"{w}/blank_defline.fa", *multi_s, "-outfmt", "7"], "losat-rejects"),
        ("audit.empty_record", ["-query", f"{w}/empty_record.fa", *multi_s, "-outfmt", "6"], "losat-rejects"),
        ("audit.empty_records_only", ["-query", f"{w}/empty_records_only.fa", *multi_s, "-outfmt", "6"], "losat-rejects"),
        ("audit.empty_subject_record", ["-query", f"{F}/multi_query.fasta", "-subject", f"{w}/empty_subject_record.fa", "-outfmt", "6"], "losat-rejects"),
        ("audit.penalty_-40000", ["-query", f"{F}/multi_query.fasta", *multi_s, "-penalty", "-40000", "-outfmt", "6"], "losat-rejects"),
        ("audit.reward_65538", ["-query", f"{F}/multi_query.fasta", *multi_s, "-reward", "65538", "-penalty", "-65540", "-outfmt", "6"], "losat-rejects"),
        ("audit.reward_2147483647", ["-query", f"{F}/multi_query.fasta", *multi_s, "-reward", "2147483647", "-penalty", "-1", "-outfmt", "6"], "losat-rejects"),
        ("audit.reward_32767", ["-query", f"{F}/multi_query.fasta", *multi_s, "-reward", "32767", "-penalty", "-32768", "-outfmt", "6"], "losat-rejects"),
        ("audit.megablast_gap_max", ["-query", f"{F}/multi_query.fasta", *multi_s, "-gapopen", "2147483647", "-gapextend", "2147483647", "-outfmt", "6"], "losat-rejects"),
        ("audit.megablast_gap_32767", ["-query", f"{F}/multi_query.fasta", *multi_s, "-gapopen", "32767", "-gapextend", "32767", "-outfmt", "6"], "same"),
        ("audit.megablast_gap_4000", ["-query", f"{F}/multi_query.fasta", *multi_s, "-gapopen", "4000", "-gapextend", "2000", "-outfmt", "0"], "same"),
        ("audit.large_divisible", ["-query", f"{F}/multi_query.fasta", *multi_s, "-reward", "1000", "-penalty", "-2000", "-outfmt", "6"], "same"),
        ("audit.evalue_0", ["-query", f"{F}/multi_query.fasta", *multi_s, "-evalue", "0", "-outfmt", "6"], "same-error"),
        ("audit.evalue_0.greedy_order", ["-query", f"{F}/multi_query.fasta", *multi_s, "-evalue", "0", "-task", "blastn", "-gapopen", "0", "-gapextend", "0", "-outfmt", "6"], "same-error"),
        ("audit.word_size_101", ["-query", f"{F}/multi_query.fasta", *multi_s, "-word_size", "101", "-outfmt", "6"], "same-error"),
        ("audit.word_size_101.penalty_order", ["-query", f"{F}/multi_query.fasta", *multi_s, "-word_size", "101", "-penalty", "0", "-outfmt", "6"], "same-error"),
        ("audit.all_invalid_positive_score", ["-query", f"{w}/edl933_100k.fa", "-subject", "tests/fasta/EDL933.fna", "-task", "blastn", "-reward", "4", "-penalty", "-1", "-outfmt", "6"], "same"),
        ("audit.reward_0", ["-query", f"{F}/multi_query.fasta", *multi_s, "-reward", "0", "-outfmt", "6"], "losat-rejects"),
        ("audit.piped_empty_query", ["-query", "/dev/stdin", *multi_s, "-outfmt", "6"], "losat-rejects"),
        ("audit.out_file.options_error", ["-query", f"{F}/multi_query.fasta", *multi_s, "-penalty", "0", "-out", "{OUT}"], "same-error"),
        ("audit.out_file.empty_query", ["-query", f"{w}/empty.fa", *multi_s, "-out", "{OUT}"], "same"),
        ("audit.out_file.table_error.fmt7", ["-query", f"{F}/multi_query.fasta", *multi_s, "-reward", "1", "-penalty", "-6", "-outfmt", "7", "-out", "{OUT}"], "same-error"),
        ("audit.leading_blank_line", ["-query", f"{w}/leading_blank_line.fa", *multi_s, "-outfmt", "6"], "losat-rejects"),
        ("audit.byte_order_mark", ["-query", f"{w}/byte_order_mark.fa", *multi_s, "-outfmt", "6"], "losat-rejects"),
    ]
    for task in ("megablast", "blastn"):
        for gaps in (["-reward", "1", "-penalty", "-2", "-gapopen", "5", "-gapextend", "2"],
                     ["-reward", "1", "-penalty", "-3", "-gapopen", "3", "-gapextend", "2"],
                     ["-reward", "1", "-penalty", "-3", "-gapopen", "4", "-gapextend", "4"], []):
            name = "-".join(gaps[1::2]) or "default"
            for query, expect in (("many_acgt", "same"), ("many_iupac", "losat-rejects" if gaps else "same")):
                for outfmt in ("0", "6"):
                    rows.append((f"batches.{query}.{task}.{name}.fmt{outfmt}",
                                 ["-query", f"{w}/{query}.fa", "-subject", f"{w}/big_subject.fa", "-task", task, *gaps,
                                  "-outfmt", outfmt], expect))
    return rows


def classify(ncbi: subprocess.CompletedProcess, ours: subprocess.CompletedProcess) -> str:
    same_streams = ncbi.stdout == ours.stdout and ncbi.stderr == ours.stderr
    if ncbi.returncode == 0 and ours.returncode == 0 and same_streams:
        return "same"
    if ncbi.returncode and ncbi.returncode == ours.returncode and same_streams:
        return "same-error"
    if ours.returncode and any(marker in ours.stderr for marker in REJECTION_MARKERS):
        return "losat-rejects"
    return f"DIFF exit {ncbi.returncode}/{ours.returncode}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--losat", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    make_inputs(args.work.resolve())
    unexpected = []
    print("case\texpect\tresult\tlosat_stderr")
    for name, argv, expect in cases(args.work.resolve()):
        runs = []
        for label, command in (("ncbi", [str(args.bin_dir / "blastn")]), ("losat", [str(args.losat.resolve()), "blastn"])):
            # `{OUT}`: the -out file, whose existence and bytes join the stdout.
            out = args.work.resolve() / f"out.{label}"
            out.unlink(missing_ok=True)
            run = subprocess.run([*command, *(str(out) if word == "{OUT}" else word for word in argv)], cwd=ENGINE,
                                 capture_output=True, input=b"")
            if "{OUT}" in argv:
                run.stdout += b"<out>" + (out.read_bytes() if out.exists() else b"<missing>")
            runs.append(run)
        ncbi, ours = runs
        result = classify(ncbi, ours)
        if result != expect:
            unexpected.append(name)
        first = ours.stderr.decode(errors="replace").strip().splitlines()
        print("\t".join([name, expect, result, first[0][:200] if first else ""]))
    print(f"# cases={len(cases(args.work.resolve()))} unexpected={unexpected}")
    return 1 if unexpected else 0


if __name__ == "__main__":
    sys.exit(main())
