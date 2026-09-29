#!/usr/bin/env python3
"""Compare BLASTN input handling and query batches with NCBI (Session S07+; comparison only).

Writes derived inputs into WORK (from the repository's inputs and fixed seeds), runs each
case in NCBI and LOSAT (standard input empty unless the case gives a file; `-out` files
compared with stdout), and
prints one line per case: `same`, `same-error`,
`losat-rejects` (NCBI succeeds; LOSAT fails with a message that names what it does not
support, whether NCBI succeeds or fails otherwise), `arg-error` (an argument that NCBI's
argument parser rejects, exit 1, and LOSAT's clap parser rejects, exit 2), `both-fail`
(both fail with the same exit status, each with its own message, and write the same
output), `timeout` or `DIFF`, with the `expect` column of the case. Exits 1 when a result
is not the expected one. A case may give environment variables, a file for standard
input, named pipes with the files that a writer puts into them in order, and a named pipe
for -out that a reader empties.

Usage: check_inputs.py --bin-dir DIR --losat LOSAT --work DIR
"""
from __future__ import annotations

import argparse
import os
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
    (work / "white_space_subject.fa").write_text(" \n\t\n")
    subject_tab = subject.replace(">msA close homolog", ">msA\tclose homolog", 1)
    (work / "tab_subject.fa").write_text(subject_tab)
    (work / "leading_blank_subject.fa").write_text("\n" + subject)
    # One query and 15 subjects that it hits (more than NCBI's wrapped preliminary hit list
    # of 10), from fixed seeds.
    rng15 = random.Random(15)
    (work / "one_query.fa").write_text(f">q1\n{q}\n")
    (work / "fifteen_subjects.fa").write_text("".join(
        f">s{index}\n" + "".join(rng15.choice("ACGT") if rng15.random() < 0.01 * (index + 1) else c for c in q) + "\n"
        for index in range(15)))
    (work / "empty_record_subject.fa").write_text(">s0 empty\n" + subject)
    # A query with an inverted repeat whose two strands met in one entry of LOSAT's
    # one-strand diagonal table (the sixth audit round): offsets 100 and 217 of 700.
    rng_ir = random.Random(3)
    ir_query = [rng_ir.choice("ACGT") for _ in range(700)]
    repeat = "".join(rng_ir.choice("ACGT") for _ in range(60))
    ir_query[100:160] = repeat
    ir_query[217:277] = repeat[::-1].translate(str.maketrans("ACGT", "TGCA"))
    (work / "inverted_repeat_query.fa").write_text(">ir\n" + "".join(ir_query) + "\n")
    flank = lambda: "".join(rng_ir.choice("ACGT") for _ in range(300))  # noqa: E731
    (work / "inverted_repeat_subject.fa").write_text(">irs\n" + flank() + repeat + flank() + "\n")
    # The seventh audit round: a subject that starts inside a tandem repeat of which the
    # query has more copies (the gap reduction read before the subject); a bit score in
    # (99.9, 100); slices where the preliminary -subject_besthit filter matters.
    tail = "CTCACAAATGTCCATGCTCACAAATGTCCATCATAGGCTAGTATCTATTAGGCTTTGAATTCCGCCTTGAGGGATCACAGGGAACCC"
    (work / "repeat_query.fa").write_text(">rq\nCTCTGCGCTACAAGCTCACAAATGTCCATG" + tail + "\n")
    (work / "repeat_subject.fa").write_text(">rs\n" + tail + "\n")
    bits_q = "GATCCTTAGGCTACGTTAGCCATGAGCTAACGTTGCAGCTATCGATGCCTA"
    (work / "bits_query.fa").write_text(">bq\n" + bits_q + "\n")
    (work / "bits_subject.fa").write_text(">bs\n" + "T" * 10 + bits_q + "A" * 10 + "\n")
    genome = lambda name: "".join(line.strip() for line in (ENGINE / "tests/fasta" / f"{name}.fasta").read_text().splitlines()[1:])  # noqa: E731
    (work / "besthit_query.fa").write_text(">bhq\n" + genome("LC738871")[80316:89722] + "\n")
    (work / "besthit_subject.fa").write_text(">bhs\n" + genome("PemoMJNVB")[306307:330257] + "\n")
    # The eighth audit round: NCBI resolves the ambiguity codes of a subject with CRandom
    # for the preliminary search (Sakai against EDL933, whose R at subject 1804 and K at
    # 231204 matter).
    genome_fna = lambda name: "".join(line.strip() for line in (ENGINE / "tests/fasta" / f"{name}.fna").read_text().splitlines()[1:])  # noqa: E731
    sakai, edl933 = genome_fna("Sakai"), genome_fna("EDL933")
    (work / "ambiguity_sakai_query.fa").write_text(">aq\n" + sakai[3349962:3361957] + "\n")
    (work / "ambiguity_edl933_subject.fa").write_text(">as\n" + edl933[3423311:3428607] + "\n")
    (work / "ambiguity_besthit_query.fa").write_text(">bq\n" + sakai[230984:231784] + "\n")
    # The ninth audit round: NCBI reads a line that starts with ">?" as a gap in the
    # sequence (">?_" as a defline without the prefix), and warns about a title that ends
    # with 20 nucleotide letters.
    (work / "gap_line_subject.fa").write_text(">s1\n" + edl933[99000:101000] + "\n>?100\n" + edl933[101000:103000] + "\n")
    (work / "gap_span_query.fa").write_text(">q1\n" + edl933[100500:101500] + "\n")
    (work / "gap_line_query.fa").write_text(">q1\n" + edl933[100500:101000] + "\n>?100\n" + edl933[101000:101500] + "\n")
    (work / "underscore_query.fa").write_text(">?_abc\n" + q + "\n")
    nucleotides = "ACGTACGTACGTACGTACGTA"
    (work / "title_nucleotides.fa").write_text(f">q1 {nucleotides}\n{q}\n")
    (work / "title_nucleotides_crlf.fa").write_text(f">q1 {nucleotides}\r\n{q}\r\n")
    (work / "title_nucleotides_space.fa").write_text(f">q1 {nucleotides} \n{q}\n")
    (work / "id_nucleotides.fa").write_text(f">{nucleotides}\n{q}\n")
    (work / "title_nucleotides_subject.fa").write_text(f">msA {nucleotides}\n" + subject.split("\n", 1)[1])


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
        ("audit.penalty_-40000", ["-query", f"{F}/multi_query.fasta", *multi_s, "-penalty", "-40000", "-outfmt", "6"], "same-error"),
        ("audit.reward_65538", ["-query", f"{F}/multi_query.fasta", *multi_s, "-reward", "65538", "-penalty", "-65540", "-outfmt", "6"], "same"),
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
        # The second audit round.
        ("audit2.empty_subject", ["-query", f"{F}/multi_query.fasta", "-subject", f"{w}/empty.fa"], "same-error"),
        ("audit2.white_space_subject", ["-query", f"{F}/multi_query.fasta", "-subject", f"{w}/white_space_subject.fa", "-outfmt", "7"], "same-error"),
        ("audit2.empty_query_and_subject", ["-query", f"{w}/empty.fa", "-subject", f"{w}/empty.fa", "-outfmt", "6"], "same-error"),
        ("audit2.penalty_16bit.word_size", ["-query", f"{F}/multi_query.fasta", *multi_s, "-penalty", "-40000", "-word_size", "101"], "same-error"),
        ("audit2.penalty_16bit.evalue", ["-query", f"{F}/multi_query.fasta", *multi_s, "-penalty", "-40000", "-evalue", "0"], "same-error"),
        ("audit2.penalty_16bit.greedy", ["-query", f"{F}/multi_query.fasta", *multi_s, "-penalty", "-40000", "-task", "blastn", "-gapopen", "0", "-gapextend", "0"], "same-error"),
        ("audit2.reward_16bit_0.penalty_0", ["-query", f"{F}/multi_query.fasta", *multi_s, "-reward", "65536", "-penalty", "0", "-outfmt", "6"], "losat-rejects"),
        ("audit2.reward_16bit_0.word_size", ["-query", f"{F}/multi_query.fasta", *multi_s, "-reward", "65536", "-penalty", "0", "-word_size", "101"], "same-error"),
        ("audit2.scores_16bit_wrap.fmt0", ["-query", f"{F}/multi_query.fasta", *multi_s, "-reward", "65537", "-penalty", "-65538"], "same"),
        ("audit2.scores_16bit_wrap.fmt6", ["-query", f"{F}/multi_query.fasta", *multi_s, "-reward", "65537", "-penalty", "-65538", "-outfmt", "6"], "same"),
        ("audit2.reward_16bit_negative", ["-query", f"{F}/multi_query.fasta", *multi_s, "-reward", "100000", "-penalty", "-1"], "losat-rejects"),
        ("audit2.evalue_negative.word_size", ["-query", f"{F}/multi_query.fasta", *multi_s, "-evalue", "-1", "-word_size", "101"], "same-error"),
        ("audit2.evalue_minus_inf", ["-query", f"{F}/multi_query.fasta", *multi_s, "-evalue", "-inf", "-outfmt", "6"], "same-error"),
        ("audit2.evalue_overflow", ["-query", f"{F}/multi_query.fasta", *multi_s, "-evalue", "1e400", "-outfmt", "6"], "losat-rejects"),
        ("audit2.missing_query.out", ["-query", f"{w}/missing.fa", *multi_s, "-out", "{OUT}"], "same-error"),
        ("audit2.dev_null_query", ["-query", "/dev/null", *multi_s, "-outfmt", "6"], "same"),
        ("audit2.batch_size_env.penalty_0.out", ["-query", f"{F}/multi_query.fasta", *multi_s, "-penalty", "0", "-out", "{OUT}"], "same-error", {"BATCH_SIZE": "1000"}),
        ("audit2.batch_size_env", ["-query", f"{F}/multi_query.fasta", *multi_s, "-outfmt", "6"], "losat-rejects", {"BATCH_SIZE": "1000"}),
        ("audit2.greedy_gap_table_error", ["-query", f"{F}/multi_query.fasta", *multi_s, "-gapopen", "32768", "-gapextend", "1", "-outfmt", "6"], "same-error"),
        ("audit2.greedy_gap_limit", ["-query", f"{F}/multi_query.fasta", *multi_s, "-gapopen", "32768", "-gapextend", "32768", "-outfmt", "6"], "losat-rejects"),
        # The third audit round: NCBI reads the subjects before it opens -out and checks the
        # options; standard input as a query or a subject; a directory; e-value forms.
        ("audit3.empty_subject.penalty_0", ["-query", f"{F}/multi_query.fasta", "-subject", f"{w}/empty.fa", "-penalty", "0", "-outfmt", "6"], "same-error"),
        ("audit3.empty_subject.evalue_0", ["-query", f"{F}/multi_query.fasta", "-subject", f"{w}/empty.fa", "-evalue", "0"], "same-error"),
        ("audit3.empty_subject.word_size_101", ["-query", f"{F}/multi_query.fasta", "-subject", f"{w}/empty.fa", "-word_size", "101", "-outfmt", "7"], "same-error"),
        ("audit3.empty_subject.max_target_seqs_1", ["-query", f"{F}/multi_query.fasta", "-subject", f"{w}/empty.fa", "-max_target_seqs", "1"], "same-error"),
        ("audit3.empty_subject.out", ["-query", f"{F}/multi_query.fasta", "-subject", f"{w}/empty.fa", "-out", "{OUT}"], "same-error"),
        ("audit3.empty_subject.missing_query", ["-query", f"{w}/missing.fa", "-subject", f"{w}/empty.fa"], "same-error"),
        ("audit3.piped_query", ["-query", "/dev/stdin", *multi_s, "-outfmt", "6"], "same", {}, f"{F}/multi_query.fasta"),
        ("audit3.piped_query.fmt0", ["-query", "/dev/stdin", *multi_s], "same", {}, f"{F}/multi_query.fasta"),
        ("audit3.piped_subject", ["-query", f"{F}/multi_query.fasta", "-subject", "/dev/stdin", "-outfmt", "7"], "same", {}, f"{F}/multi_subject.fasta"),
        ("audit3.piped_empty_subject", ["-query", f"{F}/multi_query.fasta", "-subject", "/dev/stdin", "-outfmt", "6"], "same-error"),
        ("audit3.piped_white_space_query", ["-query", "/dev/stdin", *multi_s, "-outfmt", "6"], "losat-rejects", {}, f"{w}/white_space.fa"),
        ("audit3.directory_query", ["-query", w, *multi_s, "-outfmt", "6"], "same"),
        ("audit3.directory_subject", ["-query", f"{F}/multi_query.fasta", "-subject", w, "-outfmt", "6"], "same-error"),
        ("audit3.evalue_inf", ["-query", f"{F}/multi_query.fasta", *multi_s, "-evalue", "inf", "-outfmt", "6"], "arg-error"),
        ("audit3.evalue_nan", ["-query", f"{F}/multi_query.fasta", *multi_s, "-evalue", "nan", "-outfmt", "6"], "arg-error"),
        ("audit3.evalue_space", ["-query", f"{F}/multi_query.fasta", *multi_s, "-evalue", " 1", "-outfmt", "6"], "arg-error"),
        ("audit3.evalue_plus_inf", ["-query", f"{F}/multi_query.fasta", *multi_s, "-evalue", "+inf", "-outfmt", "6"], "losat-rejects"),
        ("audit3.evalue_plus_nan", ["-query", f"{F}/multi_query.fasta", *multi_s, "-evalue", "+nan", "-outfmt", "6"], "losat-rejects"),
        ("audit3.evalue_hex", ["-query", f"{F}/multi_query.fasta", *multi_s, "-evalue", "0x10", "-outfmt", "6"], "losat-rejects"),
        ("audit3.evalue_exponent", ["-query", f"{F}/multi_query.fasta", *multi_s, "-evalue", "+1E-5", "-outfmt", "6"], "same"),
    ]
    # The fourth audit round: NCBI's integer arguments read hexadecimal; the order of the
    # subject, the query, -out and the checks; LOSAT's limits after "Query is Empty!";
    # standard input as "-"; -out -; options that NCBI blastn does not have; -outfmt.
    mq = ["-query", f"{F}/multi_query.fasta"]
    empty_q = ["-query", f"{w}/empty.fa"]
    tab_s = ["-subject", f"{w}/tab_subject.fa"]
    rows += [
        ("audit4.hex.reward", [*mq, *multi_s, "-task", "blastn", "-reward", "0x2", "-penalty", "-3", "-outfmt", "6"], "same"),
        ("audit4.hex.gaps", [*mq, *multi_s, "-task", "blastn", "-gapopen", "0x5", "-gapextend", "0X2", "-outfmt", "6"], "same"),
        ("audit4.hex.word_size", [*mq, *multi_s, "-task", "blastn", "-word_size", "0xB", "-outfmt", "7"], "same"),
        ("audit4.hex.word_size_16", [*mq, *multi_s, "-word_size", "0x10", "-outfmt", "6"], "same"),
        ("audit4.hex.max_target_seqs", [*mq, *multi_s, "-max_target_seqs", "0x3", "-outfmt", "6"], "same"),
        ("audit4.hex.penalty_0", [*mq, *multi_s, "-penalty", "0x0", "-outfmt", "6"], "same-error"),
        ("audit4.hex.no_digits", [*mq, *multi_s, "-max_target_seqs", "0x", "-outfmt", "6"], "arg-error"),
        ("audit4.tab_subject.penalty_0.out", [*mq, *tab_s, "-penalty", "0", "-out", "{OUT}"], "same-error"),
        ("audit4.tab_subject.evalue_0", [*mq, *tab_s, "-evalue", "0", "-outfmt", "6"], "same-error"),
        ("audit4.tab_subject.empty_query", [*empty_q, *tab_s, "-outfmt", "6"], "same"),
        ("audit4.tab_subject", [*mq, *tab_s, "-outfmt", "6"], "losat-rejects"),
        ("audit4.leading_blank_subject.penalty_0", [*mq, "-subject", f"{w}/leading_blank_subject.fa", "-penalty", "0"], "losat-rejects"),
        ("audit4.missing_subject", [*mq, "-subject", f"{w}/missing.fa", "-outfmt", "6"], "same-error"),
        ("audit4.missing_query_and_subject", ["-query", f"{w}/missing.fa", "-subject", f"{w}/missing.fa"], "same-error"),
        ("audit4.missing_query.bad_out", ["-query", f"{w}/missing.fa", *multi_s, "-out", f"{w}/missing_dir/out.txt"], "same-error"),
        ("audit4.bad_out", [*mq, *multi_s, "-out", f"{w}/missing_dir/out.txt"], "same-error"),
        ("audit4.directory_out", [*mq, *multi_s, "-out", w], "same-error"),
        ("audit4.empty_query.reward_0", [*empty_q, *multi_s, "-reward", "0"], "same"),
        ("audit4.empty_query.score_range", [*empty_q, *multi_s, "-reward", "5000", "-penalty", "-1"], "same"),
        ("audit4.empty_query.evalue_plus_inf", [*empty_q, *multi_s, "-evalue", "+inf"], "same"),
        ("audit4.empty_query.batch_size", [*empty_q, *multi_s, "-outfmt", "6"], "same", {"BATCH_SIZE": "1000"}),
        ("audit4.evalue_plus_nan.empty_subject", [*mq, "-subject", f"{w}/empty.fa", "-evalue", "+nan"], "same-error"),
        ("audit4.evalue_plus_nan.penalty_0", [*mq, *multi_s, "-evalue", "+nan", "-penalty", "0"], "same-error"),
        ("audit4.stdin_query.pipe", ["-query", "-", *multi_s, "-outfmt", "6"], "same", {}, f"{F}/multi_query.fasta"),
        ("audit4.stdin_query.file", ["-query", "-", *multi_s, "-outfmt", "6"], "same", {}, ("file", f"{F}/multi_query.fasta")),
        ("audit4.stdin_query.default", [*multi_s, "-outfmt", "7"], "same", {}, ("file", f"{F}/multi_query.fasta")),
        ("audit4.stdin_subject.file", [*mq, "-subject", "-", "-outfmt", "6"], "same", {}, ("file", f"{F}/multi_subject.fasta")),
        ("audit4.stdin_subject.pipe", [*mq, "-subject", "-", "-outfmt", "6"], "same", {}, f"{F}/multi_subject.fasta"),
        ("audit4.stdin_both.file", ["-query", "-", "-subject", "-", "-outfmt", "6"], "losat-rejects", {}, ("file", f"{F}/multi_subject.fasta")),
        ("audit4.dev_stdin_both.file", ["-query", "/dev/stdin", "-subject", "/dev/stdin", "-outfmt", "6"], "same", {}, ("file", f"{F}/multi_subject.fasta")),
        ("audit4.stdin_both.pipe", ["-query", "-", "-subject", "-", "-outfmt", "6"], "losat-rejects", {}, f"{F}/multi_subject.fasta"),
        ("audit4.stdin_empty_query.file", ["-query", "-", *multi_s], "same", {}, ("file", f"{w}/empty.fa")),
        ("audit4.out_stdout", [*mq, *multi_s, "-outfmt", "6", "-out", "-"], "same"),
        ("audit4.option.verbose", [*mq, *multi_s, "-verbose"], "arg-error"),
        ("audit4.option.limit_lookup", [*mq, *multi_s, "-limit_lookup"], "arg-error"),
        ("audit4.option.max_db_word_count", [*mq, *multi_s, "-max_db_word_count", "30"], "arg-error"),
        ("audit4.option.min_hit_length", [*mq, *multi_s, "-min_hit_length", "10"], "arg-error"),
        ("audit4.outfmt.custom_0", [*mq, *multi_s, "-outfmt", "0 qaccver"], "same"),
        ("audit4.outfmt.plus_6", [*mq, *multi_s, "-outfmt", "+6"], "same"),
        ("audit4.outfmt.leading_zero", [*mq, *multi_s, "-outfmt", "07"], "same"),
        ("audit4.outfmt.word", [*mq, *multi_s, "-outfmt", "abc"], "same-error"),
        ("audit4.outfmt.out_of_range", [*mq, *multi_s, "-outfmt", "99"], "same-error"),
        ("audit4.outfmt.no_break_space", [*mq, *multi_s, "-outfmt", "6\u00a0"], "same-error"),
        ("audit4.outfmt.empty", [*mq, *multi_s, "-outfmt", ""], "same-error"),
        ("audit4.outfmt.word.missing_subject", [*mq, "-subject", f"{w}/missing.fa", "-outfmt", "abc"], "same-error"),
        ("audit4.outfmt.delimiter", [*mq, *multi_s, "-outfmt", "6 delim=,"], "losat-rejects"),
        # Named pipes: the writer fills the subject and then the query, the order in which
        # NCBI opens them.
        ("audit4.fifo.subject_then_query", ["-query", f"{w}/fifo_q", "-subject", f"{w}/fifo_s", "-outfmt", "6"], "same", {}, None,
         [(f"{w}/fifo_s", f"{F}/multi_subject.fasta"), (f"{w}/fifo_q", f"{F}/multi_query.fasta")]),
        ("audit4.fifo.query", ["-query", f"{w}/fifo_q", *multi_s], "same", {}, None, [(f"{w}/fifo_q", f"{F}/multi_query.fasta")]),
        ("audit4.fifo.empty_query", ["-query", f"{w}/fifo_q", *multi_s, "-outfmt", "6"], "losat-rejects", {}, None, [(f"{w}/fifo_q", "/dev/null")]),
    ]
    # The fifth audit round: -outfmt delimiters; a constrained integer argument of a bare
    # 0x; NCBI's 32-bit preliminary hit list size; -dust as NCBI reads it; -out opened
    # once; empty and long file names; LOSAT's thread limit; the subject records.
    one = ["-query", f"{w}/one_query.fa", "-subject", f"{w}/fifteen_subjects.fa"]
    rows += [
        ("audit5.outfmt.delim_without_value", [*mq, *multi_s, "-outfmt", "0 delim"], "same-error"),
        ("audit5.outfmt.word_delim", [*mq, *multi_s, "-outfmt", "abc delim"], "same-error"),
        ("audit5.outfmt.out_of_range_delim", [*mq, *multi_s, "-outfmt", "99 delim qaccver"], "same-error"),
        ("audit5.outfmt.pairwise_delim", [*mq, *multi_s, "-outfmt", "0 delim=,"], "same"),
        ("audit5.outfmt.empty_delim", [*mq, *multi_s, "-outfmt", "6 delim="], "same"),
        ("audit5.outfmt.tabular_delim", [*mq, *multi_s, "-outfmt", "6 delim=;"], "losat-rejects"),
        ("audit5.bare_0x.reward", [*empty_q, *multi_s, "-reward", "0x"], "arg-error"),
        ("audit5.bare_0x.penalty", [*mq, *multi_s, "-penalty", "0X"], "arg-error"),
        ("audit5.bare_0x.word_size", [*mq, *multi_s, "-word_size", "0x"], "arg-error"),
        ("audit5.bare_0x.gaps", [*mq, *multi_s, "-gapopen", "0x", "-gapextend", "0X", "-outfmt", "6"], "same"),
        ("audit5.hitlist.below_wrap", [*one, "-max_target_seqs", "1073741823", "-outfmt", "6"], "same"),
        ("audit5.hitlist.wrap", [*one, "-max_target_seqs", "1073741824", "-outfmt", "6"], "same"),
        ("audit5.hitlist.wrap.hex.fmt0", [*one, "-max_target_seqs", "0x40000000"], "same"),
        ("audit5.hitlist.wrap.fmt7", [*one, "-max_target_seqs", "2147483597", "-outfmt", "7"], "same"),
        ("audit5.hitlist.negative", [*one, "-max_target_seqs", "2147483598", "-outfmt", "6"], "losat-rejects"),
        ("audit5.hitlist.negative.empty_query", [*empty_q, *multi_s, "-max_target_seqs", "2147483647"], "same"),
        ("audit5.dust.double_space", [*mq, *multi_s, "-dust", "20  64 1", "-outfmt", "6"], "same-error"),
        ("audit5.dust.leading_space", [*mq, *multi_s, "-dust", " 20 64 1"], "same-error"),
        ("audit5.dust.no_break_space", [*mq, *multi_s, "-dust", "20\u00a064 1"], "same-error"),
        ("audit5.dust.capital_yes", [*mq, *multi_s, "-dust", "Yes"], "same-error"),
        ("audit5.dust.empty", [*mq, *multi_s, "-dust", ""], "same-error"),
        ("audit5.dust.hex", [*mq, *multi_s, "-dust", "0x14 64 1"], "same-error"),
        ("audit5.dust.negative_level", [*mq, *multi_s, "-dust", "-1 64 1", "-outfmt", "6"], "same"),
        ("audit5.dust.zero_window", [*mq, *multi_s, "-dust", "1 0 1", "-outfmt", "6"], "same"),
        ("audit5.dust.negative_linker", [*mq, *multi_s, "-dust", "20 64 -1", "-outfmt", "6"], "same"),
        ("audit5.dust.negative_window", [*mq, *multi_s, "-dust", "20 -64 1"], "same"),
        ("audit5.dust.error.penalty_0.out", [*mq, *multi_s, "-dust", "20 64", "-penalty", "0", "-out", "{OUT}"], "same-error"),
        ("audit5.dust.error.empty_subject", [*mq, "-subject", f"{w}/empty.fa", "-dust", "20 64"], "same-error"),
        ("audit5.dust.error.empty_query", [*empty_q, *multi_s, "-dust", "20 64"], "same-error"),
        ("audit5.out_fifo", [*mq, *multi_s, "-outfmt", "6", "-out", f"{w}/fifo_out"], "same", {}, None, [], f"{w}/fifo_out"),
        ("audit5.out_fifo.fmt0", [*mq, *multi_s, "-out", f"{w}/fifo_out"], "same", {}, None, [], f"{w}/fifo_out"),
        ("audit5.empty_name.query", ["-query", "", *multi_s], "same-error"),
        ("audit5.empty_name.subject", [*mq, "-subject", ""], "same-error"),
        ("audit5.empty_name.out", [*mq, *multi_s, "-out", ""], "same-error"),
        ("audit5.long_name.out", [*mq, "-subject", f"{w}/empty.fa", "-out", f"{w}/" + "o" * 256], "arg-error"),
        ("audit5.long_name.out.255", [*mq, *multi_s, "-outfmt", "6", "-out", f"{w}/" + "o" * 255], "same"),
        ("audit5.threads.above_rayon", [*mq, *multi_s, "-num_threads", "70000", "-outfmt", "6"], "losat-rejects"),
        ("audit5.empty_record_subject.penalty_0", [*mq, "-subject", f"{w}/empty_record_subject.fa", "-penalty", "0"], "same-error"),
        ("audit5.empty_record_subject.empty_query", [*empty_q, "-subject", f"{w}/empty_record_subject.fa"], "same"),
        ("audit5.empty_record_subject", [*mq, "-subject", f"{w}/empty_record_subject.fa", "-outfmt", "6"], "losat-rejects"),
        ("audit5.chunk_size.blank", [*mq, *multi_s, "-outfmt", "6"], "same", {"CHUNK_SIZE": "  "}),
        ("audit5.chunk_size.empty", [*mq, *multi_s, "-outfmt", "6"], "same", {"CHUNK_SIZE": ""}),
    ]
    # The sixth audit round: the diagonal table of a single query covers both strands;
    # NCBI's other tasks and options; an -out name measured as NCBI's CDirEntry does.
    ir = ["-query", f"{w}/inverted_repeat_query.fa", "-subject", f"{w}/inverted_repeat_subject.fa"]
    strand = ["-query", f"{F}/strand_query.fasta", "-subject", f"{F}/strand_subject.fasta"]
    rows += [
        ("audit6.diagonals.inverted_repeat", [*ir, "-outfmt", "6"], "same"),
        ("audit6.diagonals.inverted_repeat.blastn", [*ir, "-task", "blastn", "-outfmt", "6"], "same"),
        ("audit6.diagonals.inverted_repeat.fmt0", ir, "same"),
        ("audit6.diagonals.strand", [*strand, "-task", "blastn", "-word_size", "5", "-evalue", "100", "-outfmt", "6"], "same"),
        ("audit6.task.dc_megablast", [*mq, *multi_s, "-task", "dc-megablast"], "losat-rejects"),
        ("audit6.task.blastn_short", [*mq, *multi_s, "-task", "blastn-short"], "losat-rejects"),
        ("audit6.task.rmblastn", [*mq, *multi_s, "-task", "rmblastn"], "losat-rejects"),
        ("audit6.task.capitals", [*mq, *multi_s, "-task", "BLASTN"], "arg-error"),
        ("audit6.option.strand", [*mq, *multi_s, "-strand", "plus"], "losat-rejects"),
        ("audit6.option.ungapped", [*mq, *multi_s, "-ungapped"], "losat-rejects"),
        ("audit6.option.xdrop_gap", [*mq, *multi_s, "-xdrop_gap", "30"], "losat-rejects"),
        ("audit6.option.num_alignments", [*mq, *multi_s, "-num_alignments", "5"], "losat-rejects"),
        ("audit6.out_name.trailing_dot", [*mq, *multi_s, "-out", f"{w}/missing_dir/" + "o" * 256 + "/."], "same-error"),
    ]
    repeat = ["-query", f"{w}/repeat_query.fa", "-subject", f"{w}/repeat_subject.fa"]
    bits = ["-query", f"{w}/bits_query.fa", "-subject", f"{w}/bits_subject.fa", "-reward", "2", "-penalty", "-5", "-dust", "no"]
    besthit = ["-query", f"{w}/besthit_query.fa", "-subject", f"{w}/besthit_subject.fa", "-task", "blastn", "-subject_besthit"]
    rows += [
        ("audit7.gap_reduction.subject_start", repeat, "same"),
        ("audit7.gap_reduction.subject_start.fmt6", [*repeat, "-outfmt", "6"], "same"),
        ("audit7.bit_score_99.fmt0", bits, "same"),
        ("audit7.bit_score_99.fmt6", [*bits, "-outfmt", "6"], "same"),
        ("audit7.bit_score_99.fmt7", [*bits, "-outfmt", "7"], "same"),
        ("audit7.subject_besthit.prelim", [*besthit, "-outfmt", "6"], "same"),
        ("audit7.subject_besthit.prelim.fmt0", besthit, "same"),
    ]
    # A 451-bp subject with one Y inside a 50-bp match whose exact runs are shorter than a
    # word (docs/evidence/losat_web_e2c/inputs/): NCBI's CRandom resolves the Y to T.
    stored = Path(__file__).resolve().parent / "inputs"
    amb_s = ["-subject", str(stored / "ambiguity_subject.fa")]
    rows += [
        ("audit8.ambiguity.minimal.query_c", ["-query", str(stored / "ambiguity_query_c.fa"), *amb_s, "-outfmt", "6"], "same"),
        ("audit8.ambiguity.minimal.query_t", ["-query", str(stored / "ambiguity_query_t.fa"), *amb_s, "-outfmt", "6"], "same"),
        ("audit8.ambiguity.minimal.query_t.fmt0", ["-query", str(stored / "ambiguity_query_t.fa"), *amb_s], "same"),
        ("audit8.ambiguity.sakai_edl933.blastn", ["-query", f"{w}/ambiguity_sakai_query.fa", "-subject", f"{w}/ambiguity_edl933_subject.fa", "-task", "blastn", "-outfmt", "6"], "same"),
        ("audit8.ambiguity.besthit_edl933", ["-query", f"{w}/ambiguity_besthit_query.fa", "-subject", "tests/fasta/EDL933.fna", "-subject_besthit", "-outfmt", "6"], "same"),
        ("audit9.gap_line.subject", ["-query", f"{w}/gap_span_query.fa", "-subject", f"{w}/gap_line_subject.fa", "-outfmt", "6"], "losat-rejects"),
        ("audit9.gap_line.query", ["-query", f"{w}/gap_line_query.fa", *multi_s, "-outfmt", "6"], "losat-rejects"),
        ("audit9.gap_line.underscore", ["-query", f"{w}/underscore_query.fa", *multi_s, "-outfmt", "6"], "losat-rejects"),
        ("audit9.title_warning.query.fmt0", ["-query", f"{w}/title_nucleotides.fa", *multi_s], "same"),
        ("audit9.title_warning.query.fmt6", ["-query", f"{w}/title_nucleotides.fa", *multi_s, "-outfmt", "6"], "same"),
        ("audit9.title_warning.query.fmt7", ["-query", f"{w}/title_nucleotides.fa", *multi_s, "-outfmt", "7"], "same"),
        ("audit9.title_warning.crlf", ["-query", f"{w}/title_nucleotides_crlf.fa", *multi_s, "-outfmt", "6"], "same"),
        ("audit9.title_warning.id", ["-query", f"{w}/id_nucleotides.fa", *multi_s, "-outfmt", "6"], "same"),
        ("audit9.title_warning.trailing_space", ["-query", f"{w}/title_nucleotides_space.fa", *multi_s, "-outfmt", "6"], "losat-rejects"),
        ("audit9.title_warning.subject", ["-query", f"{F}/multi_query.fasta", "-subject", f"{w}/title_nucleotides_subject.fa", "-outfmt", "6"], "same"),
        ("audit9.title_warning.subject.penalty_0", ["-query", f"{F}/multi_query.fasta", "-subject", f"{w}/title_nucleotides_subject.fa", "-penalty", "0"], "same-error"),
        ("audit9.title_warning.query.empty_subject", ["-query", f"{w}/title_nucleotides.fa", "-subject", f"{w}/empty.fa"], "same-error"),
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
    # An argument that both argument parsers reject, each with its own message.
    if (ncbi.returncode == 1 and b"CArgException" in ncbi.stderr and ours.returncode == 2
            and ours.stderr.startswith(b"error: ") and ncbi.stdout == ours.stdout):
        return "arg-error"
    # Both fail with the same status, each with its own message, and write the same output.
    if ncbi.returncode and ncbi.returncode == ours.returncode and ncbi.stdout == ours.stdout:
        return "both-fail"
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
    for name, argv, expect, *extra in cases(args.work.resolve()):
        # Optional: environment variables, and a file (relative to the engine) given as
        # standard input.
        env = {**os.environ, **(extra[0] if extra else {})}
        # A path is written to standard input through a pipe; ("file", path) redirects
        # standard input from the file.
        source = extra[1] if len(extra) > 1 else None
        stdin_file = source[1] if isinstance(source, tuple) else None
        stdin = (ENGINE / source).read_bytes() if isinstance(source, str) else b""
        fifos = extra[2] if len(extra) > 2 else []
        # A named pipe given as -out, read by a reader into a file that joins the stdout.
        out_fifo = extra[3] if len(extra) > 3 else None
        runs = []
        for label, command in (("ncbi", [str(args.bin_dir / "blastn")]), ("losat", [str(args.losat.resolve()), "blastn"])):
            # `{OUT}`: the -out file, whose existence and bytes join the stdout.
            out = args.work.resolve() / f"out.{label}"
            out.unlink(missing_ok=True)
            writer = None
            if fifos:
                for fifo, _ in fifos:
                    Path(fifo).unlink(missing_ok=True)
                    os.mkfifo(fifo)
                script = " && ".join(f"cat '{ENGINE / source}' > '{fifo}'" for fifo, source in fifos)
                writer = subprocess.Popen(["sh", "-c", script])
            reader = None
            if out_fifo:
                Path(out_fifo).unlink(missing_ok=True)
                Path(f"{out_fifo}.read").unlink(missing_ok=True)
                os.mkfifo(out_fifo)
                reader = subprocess.Popen(["sh", "-c", f"cat '{out_fifo}' > '{out_fifo}.read'"])
            words = [*command, *(str(out) if word == "{OUT}" else word for word in argv)]
            try:
                if stdin_file:
                    with open(ENGINE / stdin_file, "rb") as handle:
                        run = subprocess.run(words, cwd=ENGINE, capture_output=True, stdin=handle, env=env, timeout=60)
                else:
                    run = subprocess.run(words, cwd=ENGINE, capture_output=True, input=stdin, env=env, timeout=60)
            except subprocess.TimeoutExpired as expired:
                run = subprocess.CompletedProcess(expired.cmd, -9, b"", b"<timeout>")
            if writer:
                writer.kill()
                writer.wait()
            if reader:
                try:
                    reader.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    reader.kill()
                    reader.wait()
                run.stdout += b"<fifo>" + Path(f"{out_fifo}.read").read_bytes()
            if "{OUT}" in argv:
                run.stdout += b"<out>" + (out.read_bytes() if out.exists() else b"<missing>")
            runs.append(run)
        ncbi, ours = runs
        result = "timeout" if -9 in (ncbi.returncode, ours.returncode) else classify(ncbi, ours)
        if result != expect:
            unexpected.append(name)
        first = ours.stderr.decode(errors="replace").strip().splitlines()
        print("\t".join([name, expect, result, first[0][:200] if first else ""]))
    print(f"# cases={len(cases(args.work.resolve()))} unexpected={unexpected}")
    return 1 if unexpected else 0


if __name__ == "__main__":
    sys.exit(main())
