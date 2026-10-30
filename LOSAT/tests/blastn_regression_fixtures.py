#!/usr/bin/env python3
"""Frozen NCBI BLAST+ 2.17.0 BLASTN outputs for fast pull-request regressions.

The cases cover NCBI paths that the Gate A manifest and the outfmt 0 fixtures do
not reach: preliminary hit lists that overflow (more subjects with hits than
prelim_hitlist_size, ties broken by subject order), query batches, a split query
batch, IUPAC ambiguity, lowercase masking and non-default scores.

- generate: writes the inputs to LOSAT/tests/fixtures/blastn_regression/inputs/
  (deterministic; the files are committed, so later runs never regenerate them).
- freeze --oracle BLASTN: runs NCBI BLAST+ (comparison oracle only) from LOSAT/ and
  writes <case>.out, <case>.err and the hash columns of manifest.tsv.
- check --losat LOSAT: runs `losat blastn <argv> <losat_extra>` from LOSAT/ and
  compares stdout, stderr and the exit status with the frozen files.

A case may set environment variables (`env`, space-separated KEY=VALUE) for both
programs; every other variable that changes NCBI's batches or report is unset. A case of
an approved exception (docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md), where NCBI fails
with `env`, also sets `oracle_env`: `freeze` runs NCBI with it instead, a configuration
whose output LOSAT's result with `env` was shown to equal. An approved exception that is
an argument (an `Int4` that overflows in NCBI) has its NCBI arguments in `ORACLE_ARGV`
instead, the nearest value that does not overflow; the manifest keeps LOSAT's. A case whose id ends in
`.merged` runs with standard error merged into standard output (`2>&1`), which records
the order of the warnings and the report.

The split case slices LOSAT/tests/fasta/EDL933.fna at run time into
LOSAT/target/blastn_regression/ (outfmt 6 prints no file names).

Usage:
  blastn_regression_fixtures.py generate
  blastn_regression_fixtures.py freeze --oracle /path/to/blastn [--only PREFIX]
  blastn_regression_fixtures.py check --losat BIN [--jobs N] [--out TSV] [--threads N] [--only PREFIX]
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import os
import random
import shlex
import subprocess
import sys
from pathlib import Path

ENGINE = Path(__file__).resolve().parents[1]
FIXTURES = ENGINE / "tests/fixtures/blastn_regression"
INPUTS = FIXTURES / "inputs"
MANIFEST = FIXTURES / "manifest.tsv"
RUNTIME = ENGINE / "target/blastn_regression"
FASTA = ENGINE / "tests/fasta"
E2C_INPUTS = ENGINE.parent / "docs/evidence/losat_web_e2c/inputs"
FIELDS = ["case_id", "argv", "losat_extra", "exit", "stdout_sha256", "stdout_bytes", "stderr_sha256", "env",
          "oracle_env"]
# NCBI reads these and each one changes the batches or the report.
REPORT_ENV = ("BL2SEQ_LEGACY", "CTOOLKIT_COMPATIBLE", "OLD_FSC", "BATCH_SIZE", "CHUNK_SIZE", "ADAPTIVE_CBS",
              "OVERLAP_CHUNK_SIZE", "PRE_FETCH_SEQS_LIMIT")

I = "tests/fixtures/blastn_regression/inputs"
R = "target/blastn_regression"
P = f"-query {I}/q3k.fa -subject {I}/s600.fa"
T = f"-query {I}/q3k.fa -subject {I}/s600ties.fa"
B = f"-query {I}/mq.fa -subject {I}/s60k.fa"
S = f"-query {R}/edl933_2m.fa -subject {I}/sakai_60k.fa"
A = f"-query {I}/ambiguity_query_c.fa -subject {I}/ambiguity_subject.fa"
O = "tests/fasta/outfmt0"
C = "tests/fixtures/tblastx_regression/inputs"
CT = f"-query {C}/culling_ties_query.fa -subject {C}/culling_ties_subject.fa"
CM = f"-query {C}/culling_multi_query.fa -subject {C}/culling_multi_subject.fa"
CTR = f"-query {C}/culling_ties_query.fa -subject {I}/culling_ties_subject_rev.fa"
CMR = f"-query {C}/culling_multi_query.fa -subject {I}/culling_multi_subject_rev.fa"
RP = f"-query {I}/rep_query.fa -subject {I}/rep_subject.fa"
PL = f"-query {I}/pal_query.fa -subject {I}/pal_subject.fa"
WM = f"-query {I}/wm_query.fa -subject {I}/s60k.fa"
WM3 = f"-query {I}/wm_query3.fa -subject {I}/s60k.fa"
WMC = "tests/winmask/counts/mg1655_u9.obinary"
WMB = "tests/winmask/counts/dm2L_u9.binary"
AB = f"-query {I}/culling_abricate_query.fa -subject {I}/culling_abricate_genes.fa -task blastn -dust no" \
     " -perc_identity 80 -evalue 1E-20 -max_target_seqs 10000"
# O4a (blastn -ungapped): the inputs of generate_ungapped_inputs (the t2 pair, rep, B and A are the ones above).
T2 = f"-query {I}/t2_query.fa -subject {I}/t2_subject.fa"
AMB = f"-query {I}/ambrun_q.fa -subject {I}/ambrun_s.fa"
RE = f"-query {I}/re_q.fa -subject {I}/re_s.fa"
EDLU = f"-query {R}/edl_ungapped_q.fa -subject tests/fasta/EDL933.fna"
MX = f"-query {I}/mix_q.fa -subject {I}/mix_s.fa"
SM = f"-query {I}/strand_mixed_query.fa -subject {I}/s60k.fa"
MLST = f"-query {I}/culling_abricate_query.fa -subject {I}/culling_abricate_genes.fa -ungapped -dust no" \
       " -word_size 32 -max_target_seqs 100000 -perc_identity 95 -evalue 1E-20"
DCQ = f"-query {O}/dc_query.fasta -subject {O}/dc_subject.fasta"
DC = f"{DCQ} -task dc-megablast"
DVQ = f"-query {O}/dc_div_query.fasta -subject {O}/dc_div_subject.fasta"
DV = f"{DVQ} -task dc-megablast"
SH = f"-query {O}/short_query.fasta -subject {O}/short_subject.fasta -task blastn-short"
MW = f"-query {O}/multi_query.fasta -subject {O}/multi_subject.fasta -task blastn -word_size 7 -evalue 1000"
# (case_id, NCBI argv after `blastn`, extra LOSAT-only arguments[, environment[, NCBI environment]])
_CASES = [
    ("prelim.default", f"{P} -outfmt 6", ""),
    ("prelim.task_blastn", f"{P} -task blastn -outfmt 6", ""),
    ("prelim.max1", f"{P} -max_target_seqs 1 -outfmt 6", ""),
    ("prelim.max3", f"{P} -max_target_seqs 3 -outfmt 6", ""),
    ("prelim.max6", f"{P} -max_target_seqs 6 -outfmt 6", ""),
    ("prelim.max3_besthit", f"{P} -task blastn -max_target_seqs 3 -subject_besthit -outfmt 6", ""),
    ("prelim.max5_maxhsps1", f"{P} -max_target_seqs 5 -max_hsps 1 -outfmt 6", ""),
    ("prelim.max2_evalue", f"{P} -max_target_seqs 2 -evalue 1e-5 -outfmt 6", ""),
    ("prelim.max3_fmt7", f"{P} -max_target_seqs 3 -outfmt 7", ""),
    ("prelim.max3_fmt0", f"{P} -max_target_seqs 3 -outfmt 0", ""),
    ("prelim.default_threads4", f"{P} -outfmt 6", "-num_threads 4"),
    ("prelim.scores_3_4", f"{P} -task blastn -reward 3 -penalty -4 -gapopen 10 -gapextend 3 -outfmt 6", ""),
    ("prelim.scores_1_1", f"{P} -task blastn -reward 1 -penalty -1 -gapopen 3 -gapextend 2 -outfmt 6", ""),
    ("ties.default", f"{T} -outfmt 6", ""),
    ("ties.max3", f"{T} -max_target_seqs 3 -outfmt 6", ""),
    ("ties.max1_task_blastn", f"{T} -task blastn -max_target_seqs 1 -outfmt 6", ""),
    ("batches.default", f"{B} -outfmt 6", ""),
    ("batches.task_blastn", f"{B} -task blastn -outfmt 6", ""),
    ("batches.word7", f"{B} -task blastn -word_size 7 -outfmt 6", ""),
    ("batches.besthit", f"{B} -subject_besthit -outfmt 6", ""),
    # -qcov_hsp_perc (blast_traceback.c:837-867): after -max_hsps and before -subject_besthit,
    # each HSP's query extent over its context's length (the -query_loc range); a subject left
    # without HSPs is removed and the hit list is sorted by best e-value again.
    ("qcov.multi30", f"{MW} -qcov_hsp_perc 30 -outfmt 6", ""),
    ("qcov.multi100_fmt7", f"{MW} -qcov_hsp_perc 100 -outfmt 7", ""),
    ("qcov.multi_maxhsps1_besthit", f"{MW} -max_hsps 1 -qcov_hsp_perc 30 -subject_besthit -outfmt 6", ""),
    ("qcov.megablast50_fmt0", f"-query {O}/multi_query.fasta -subject {O}/multi_subject.fasta -qcov_hsp_perc 50"
                              " -outfmt 0", ""),
    ("qcov.loc10", f"{P} -task blastn -query_loc 1001-2000 -qcov_hsp_perc 10 -outfmt 6", ""),
    ("qcov.none_fmt7", f"{P} -task blastn -qcov_hsp_perc 10 -outfmt 7", ""),
    ("batches.fmt7", f"{B} -task blastn -outfmt 7", ""),
    ("batches.fmt0", f"{B} -outfmt 0", ""),
    ("batches.threads4", f"{B} -task blastn -outfmt 6", "-num_threads 4"),
    ("split.task_blastn", f"{S} -task blastn -outfmt 6", ""),
    ("split.besthit", f"{S} -task blastn -subject_besthit -outfmt 6", ""),
    ("ambiguity.default", f"{A} -outfmt 6", ""),
    ("ambiguity.task_blastn_word7", f"{A} -task blastn -word_size 7 -outfmt 6", ""),
    ("lcase.blastn", f"-query {I}/lcase_island_blastn_query.fa -subject {I}/lcase_island_blastn_subject.fa"
                     " -task blastn -lcase_masking -outfmt 6", ""),
    ("lcase.megablast", f"-query {I}/lcase_island_megablast_query.fa -subject {I}/lcase_island_megablast_subject.fa"
                        " -lcase_masking -outfmt 6", ""),
    # E2g T2: a gapped start left of the HSP start (negative Int4 offset in
    # BlastGetStartForGappedAlignmentNucl). t2_query.fa/t2_subject.fa are case plain208 of
    # the E2g overflow hunt (docs/evidence/losat_web_e2g/overflow_hunt/), not made by `generate`.
    ("t2.word4", f"-query {I}/t2_query.fa -subject {I}/t2_subject.fa -task blastn -word_size 4 -evalue 1e5 -outfmt 6", ""),
    # E2g T1: init hits of the two strands tied on score, subject start and length.
    ("pal.default", f"-query {I}/pal_query.fa -subject {I}/pal_subject.fa -outfmt 6", ""),
    ("pal.task_blastn", f"-query {I}/pal_query.fa -subject {I}/pal_subject.fa -task blastn -outfmt 6", ""),
    ("pal.word7_fmt0", f"-query {I}/pal_query.fa -subject {I}/pal_subject.fa -task blastn -word_size 7 -outfmt 0", ""),
    # E2g T6: diagonal array for a block of several queries of at most 8000.
    ("sq.default", f"-query {I}/sq_query.fa -subject {I}/sq_subject.fa -outfmt 6", ""),
    ("sq.task_blastn", f"-query {I}/sq_query.fa -subject {I}/sq_subject.fa -task blastn -outfmt 6", ""),
    ("sq.word7", f"-query {I}/sq_query.fa -subject {I}/sq_subject.fa -task blastn -word_size 7 -evalue 100 -outfmt 6", ""),
    ("sq.word16_fmt7", f"-query {I}/sq_query.fa -subject {I}/sq_subject.fa -word_size 16 -outfmt 7", ""),
    # E2g T5: small lookup table cells in ascending query offsets, diagonal hash.
    ("rep.megablast", f"-query {I}/rep_query.fa -subject {I}/rep_subject.fa -outfmt 6", ""),
    ("rep.task_blastn", f"-query {I}/rep_query.fa -subject {I}/rep_subject.fa -task blastn -outfmt 6", ""),
    ("rep.word9", f"-query {I}/rep_query.fa -subject {I}/rep_subject.fa -task blastn -word_size 9 -outfmt 6", ""),
    ("rep.word16", f"-query {I}/rep_query.fa -subject {I}/rep_subject.fa -word_size 16 -outfmt 6", ""),
    ("rep2.task_blastn", f"-query {I}/rep2_query.fa -subject {I}/rep_subject.fa -task blastn -outfmt 6", ""),
    ("rep2.fmt0", f"-query {I}/rep2_query.fa -subject {I}/rep_subject.fa -task blastn -word_size 10 -outfmt 0", ""),
    # E2g T4: heapified preliminary hit lists (600 subjects) traced in their stored order,
    # with gap costs beyond the tables (gapped Karlin block copied from the ungapped one).
    ("prelim.gaps10_1_2", f"{P} -task blastn -reward 1 -penalty -2 -gapopen 10 -gapextend 10 -outfmt 6", ""),
    ("prelim.gaps10_2_3_max3", f"{P} -task blastn -reward 2 -penalty -3 -gapopen 10 -gapextend 10"
                               " -max_target_seqs 3 -outfmt 6", ""),
    ("ties.gaps10_1_1", f"{T} -task blastn -reward 1 -penalty -1 -gapopen 10 -gapextend 10 -outfmt 6", ""),
    # E2g T7: BATCH_SIZE, CHUNK_SIZE and OVERLAP_CHUNK_SIZE as NCBI reads them.
    ("env.batch100", f"{B} -task blastn -outfmt 6", "", "BATCH_SIZE=100"),
    ("env.batch1000", f"{B} -outfmt 6", "", "BATCH_SIZE=1000"),
    ("env.batch5000_fmt7", f"{B} -task blastn -outfmt 7", "", "BATCH_SIZE=5000"),
    ("env.batch100000", f"{B} -outfmt 6", "", "BATCH_SIZE=100000"),
    ("env.batch_negative", f"{B} -task blastn -outfmt 6", "", "BATCH_SIZE=-1"),
    ("env.batch0", f"{B} -outfmt 6", "", "BATCH_SIZE=0"),
    ("env.chunk40000", f"{B} -task blastn -outfmt 6", "", "CHUNK_SIZE=40000"),
    ("env.chunk40000_overlap6", f"{B} -task blastn -outfmt 6", "", "CHUNK_SIZE=40000 OVERLAP_CHUNK_SIZE=6"),
    ("env.chunk2000", f"{B} -task blastn -outfmt 6", "", "CHUNK_SIZE=2000"),
    ("env.chunk1500_overlap50_fmt0", f"{B} -outfmt 0", "", "CHUNK_SIZE=1500 OVERLAP_CHUNK_SIZE=50"),
    ("env.chunk500", f"{B} -task blastn -outfmt 6", "", "CHUNK_SIZE=500"),
    ("env.chunk_negative", f"{B} -outfmt 6", "", "CHUNK_SIZE=-5"),
    ("env.chunk_blank", f"{B} -outfmt 6", "", "CHUNK_SIZE=' '"),
    ("env.chunk1000_batch5000", f"{B} -outfmt 6", "", "CHUNK_SIZE=1000 BATCH_SIZE=5000"),
    # A batch size of 0: the empty first batch fails after the outfmt 0 prolog (exit 3).
    ("env.chunk1000_fmt0", f"{B} -outfmt 0", "", "CHUNK_SIZE=1000"),
    ("env.chunk1000_fmt7", f"{B} -outfmt 7", "", "CHUNK_SIZE=1000"),
    # An overlap above about half the chunk size, just below the size where NCBI would split a
    # chunk again.
    ("env.split_overlap_near_resplit", f"{S} -task blastn -outfmt 6", "", "CHUNK_SIZE=300000 OVERLAP_CHUNK_SIZE=149000"),
    # Approved exception 1 of PD-LOSAT-NCBI-DEFECTS: with these overlaps NCBI would split a
    # query chunk again and stops with a CCoreException; LOSAT searches each chunk once, which
    # equals NCBI's output with the largest overlap that does not split a chunk again
    # (docs/evidence/losat_web_e2g/resplit/; resplit_query.fa and resplit_subject.fa come from
    # its validate.py, not from `generate`).
    ("env.resplit_overlap155000", f"{S} -task blastn -outfmt 6", "", "CHUNK_SIZE=300000 OVERLAP_CHUNK_SIZE=155000",
     "CHUNK_SIZE=300000 OVERLAP_CHUNK_SIZE=149000"),
    ("env.resplit_8000_4100_fmt0", f"-query {I}/resplit_query.fa -subject {I}/resplit_subject.fa -task blastn -outfmt 0",
     "", "CHUNK_SIZE=8000 OVERLAP_CHUNK_SIZE=4100", "CHUNK_SIZE=8000 OVERLAP_CHUNK_SIZE=4090"),
    ("env.resplit_3000_2900_megablast", f"-query {I}/resplit_query.fa -subject {I}/resplit_subject.fa -outfmt 6",
     "", "CHUNK_SIZE=3000 OVERLAP_CHUNK_SIZE=2900", "CHUNK_SIZE=3000 OVERLAP_CHUNK_SIZE=1512"),
    # A negative chunk size at most a negative overlap splits nothing.
    ("env.negative_pair_le", f"{B} -outfmt 6", "", "CHUNK_SIZE=-10 OVERLAP_CHUNK_SIZE=-5"),
    # Above it, a split needs the wrapped difference twice in the batch (E2g audits b and
    # round 2): these do not split, and NCBI searches as without the variables.
    ("env.negative_pair_wide", f"{S} -task blastn -outfmt 6", "", "CHUNK_SIZE=-1 OVERLAP_CHUNK_SIZE=-2147483648"),
    ("env.negative_pair_wide_fmt0", f"{B} -outfmt 0", "", "CHUNK_SIZE=-5 OVERLAP_CHUNK_SIZE=-2147483648"),
    ("env.split_chunk300000", f"{S} -task blastn -outfmt 6", "", "CHUNK_SIZE=300000"),
    ("env.split_overlap50", f"{S} -task blastn -outfmt 6", "", "CHUNK_SIZE=300000 OVERLAP_CHUNK_SIZE=50"),
    ("env.split_overlap500", f"{S} -task blastn -outfmt 6", "", "CHUNK_SIZE=300000 OVERLAP_CHUNK_SIZE=500"),
    ("env.split_overlap0", f"{S} -task blastn -outfmt 6", "", "CHUNK_SIZE=300000 OVERLAP_CHUNK_SIZE=0"),
    ("env.split_overlap_negative", f"{S} -task blastn -outfmt 6", "", "OVERLAP_CHUNK_SIZE=-1"),
    ("env.split_besthit_batch", f"{S} -task blastn -subject_besthit -outfmt 6", "", "BATCH_SIZE=1000 CHUNK_SIZE=300000"),
    # E2g T12: an integer PRE_FETCH_SEQS_LIMIT only decides whether sequences are fetched ahead.
    ("env.prefetch0_fmt0", f"{T} -max_target_seqs 5 -outfmt 0", "", "PRE_FETCH_SEQS_LIMIT=0"),
    ("env.prefetch5", f"{B} -outfmt 6", "", "PRE_FETCH_SEQS_LIMIT=5"),
    # E2g T10: no -subject (CBlastDatabaseArgs raises NCBI's error before -query and -out are opened).
    ("nosubject.fmt6", f"-query {I}/q3k.fa -outfmt 6", ""),
    ("nosubject.missing_query_out", "-query nonexistent_query.fa -out nonexistent_dir/out.txt", ""),
    # E2g R2: a .ncbirc (found through $HOME) with entries that change no output.
    ("ncbirc.harmless_fmt0", f"{T} -max_target_seqs 5 -outfmt 0", "",
     "HOME=tests/fixtures/blastn_regression/ncbirc_home BLAST_USAGE_REPORT=0 NCBI_CONFIG__BLAST__BLASTDB=/x"),
    # E2g T14: a failed outfmt 0 write ("BLAST failed to write output", exit 6; Linux /dev/full).
    ("write.devfull_fmt0", f"{T} -max_target_seqs 5 -outfmt 0 -out /dev/full", ""),
    # E2g T8: the Karlin-Altschul table error of a later batch, after the reports of a first
    # batch of invalid queries (no epilog).
    ("kaerror.later_batch.fmt6", "-query tests/fasta/outfmt0/edge_batch_allN.fasta -subject"
                                 " tests/fasta/blastn_parity_compact.fasta -reward 1 -penalty -6 -outfmt 6", ""),
    ("kaerror.later_batch.fmt0", "-query tests/fasta/outfmt0/edge_batch_allN.fasta -subject"
                                 " tests/fasta/blastn_parity_compact.fasta -reward 1 -penalty -6 -outfmt 0", ""),
    ("kaerror.later_batch.fmt7", "-query tests/fasta/outfmt0/edge_batch_allN.fasta -subject"
                                 " tests/fasta/blastn_parity_compact.fasta -reward 1 -penalty -6 -outfmt 7", ""),
    # E2g T9: title and invalid-query warnings per query batch.
    ("warnings.batches.fmt6", f"-query {I}/warn_query.fa -subject {I}/warn_subject.fa -outfmt 6", ""),
    ("warnings.batches.fmt0", f"-query {I}/warn_query.fa -subject {I}/warn_subject.fa -outfmt 0", ""),
    ("warnings.batch1000.fmt7", f"-query {I}/warn_query.fa -subject {I}/warn_subject.fa -task blastn -outfmt 7", "",
     "BATCH_SIZE=1000"),
    # NCBI posts the warnings on cerr, tied to cout: a batch's title warnings come before the
    # report of its first query, an invalid query's warning before its preamble.
    ("warnings.batches.fmt0.merged", f"-query {I}/warn_query.fa -subject {I}/warn_subject.fa -outfmt 0", ""),
    ("warnings.batches.fmt6.merged", f"-query {I}/warn_query.fa -subject {I}/warn_subject.fa -outfmt 6", ""),
    ("warnings.batch1000.fmt7.merged", f"-query {I}/warn_query.fa -subject {I}/warn_subject.fa -task blastn -outfmt 7",
     "", "BATCH_SIZE=1000"),
    ("warnings.batch1000.fmt0.merged", f"-query {I}/warn_query.fa -subject {I}/warn_subject.fa -outfmt 0", "",
     "BATCH_SIZE=1000"),
    ("kaerror.later_batch.fmt0.merged", "-query tests/fasta/outfmt0/edge_batch_allN.fasta -subject"
                                        " tests/fasta/blastn_parity_compact.fasta -reward 1 -penalty -6 -outfmt 0", ""),
    # The outfmt 0 stream fails at its first flush, after the version line, before any
    # query is read: no warning is written.
    ("write.devfull_warnings_fmt0", f"-query {I}/warn_query.fa -subject {I}/warn_subject.fa -outfmt 0 -out /dev/full",
     ""),
    # The prolog is written before the first batch is read and searched (E2g round-3 audit):
    # also with scores without gapped Karlin-Altschul values and an invalid first batch.
    ("write.devfull_warnings_kaerror_fmt0", f"-query {I}/warn_query.fa -subject {I}/warn_subject.fa -reward 1"
                                            " -penalty -6 -outfmt 0 -out /dev/full", ""),
    ("write.devfull_kaerror_batch_fmt0", "-query tests/fasta/outfmt0/edge_batch_allN.fasta -subject"
                                         " tests/fasta/blastn_parity_compact.fasta -reward 1 -penalty -6 -outfmt 0"
                                         " -out /dev/full", ""),
    # E2g V1: scores beyond the former 3000-unit limit (gcd-scaled table pairs).
    ("scores.10000_20000", f"{T} -max_target_seqs 20 -reward 10000 -penalty -20000 -outfmt 6", ""),
    ("scores.16383_32766_gaps", f"{T} -max_target_seqs 20 -reward 16383 -penalty -32766 -gapopen 16383"
                                " -gapextend 32766 -outfmt 6", ""),
    ("scores.24000_30000_task_blastn", f"{T} -task blastn -max_target_seqs 20 -reward 24000 -penalty -30000"
                                       " -gapopen 30000 -gapextend 30000 -outfmt 6", ""),
    ("scores.6000_21000_fmt0", f"{P} -max_target_seqs 3 -reward 6000 -penalty -21000 -gapopen 30000"
                               " -gapextend 30000 -outfmt 0", ""),
    # Approved exception 2 of PD-LOSAT-NCBI-DEFECTS: subjects whose deflines are only
    # punctuation (NCBI reads past the end of their outfmt 0 titles and crashes when such a
    # subject has hits; docs/evidence/losat_web_e2g/punct_defline/). Without hits on them NCBI
    # runs and LOSAT must equal it. punct_query.fa and punct_nohit_subject.fa come from that
    # validation, not from `generate`.
    ("punct.nohit_fmt0", f"-query {I}/punct_query.fa -subject {I}/punct_nohit_subject.fa -outfmt 0", ""),
    ("punct.nohit_task_blastn_fmt0", f"-query {I}/punct_query.fa -subject {I}/punct_nohit_subject.fa"
                                     " -task blastn -outfmt 0", ""),
    # E2g: NCBI's CArg_Double reads signed infinities and NaN (strtod, also nan(...)) and
    # NCBI searches with them as with the largest e-value.
    ("evalue.plus_inf", f"{P} -evalue +inf -outfmt 6", ""),
    ("evalue.minus_nan_task_blastn_fmt7", f"{B} -task blastn -evalue -nan -outfmt 7", ""),
    ("evalue.nan_paren_fmt0", f"{T} -max_target_seqs 5 -evalue +nan(1) -outfmt 0", ""),
    ("evalue.1e999_word7", f"{A} -task blastn -word_size 7 -evalue 1e999 -outfmt 6", ""),
    # E2g T11: showdefline.cpp kBits is "(bits)" when CTOOLKIT_COMPATIBLE is set (also empty).
    ("ctoolkit.fmt0", f"{P} -max_target_seqs 3 -outfmt 0", "", "CTOOLKIT_COMPATIBLE=1"),
    ("ctoolkit.empty_fmt0", f"{T} -task blastn -max_target_seqs 5 -outfmt 0", "", "CTOOLKIT_COMPATIBLE="),
    # Session SD (E2i): dc-megablast and blastn-short. dc_*.fasta and short_*.fasta are made
    # by docs/evidence/losat_web_e2i/make_inputs.py (LOSAT/tests/fasta/outfmt0/, also the
    # outfmt 0/7 fixtures dc.* and short.* of outfmt0_manifest.tsv). Every template type,
    # length and word size of NCBI's s_GetDiscTemplateType (blast_nalookup.c:600-643):
    ("dc.w11_coding_16", f"{DC} -word_size 11 -template_type coding -template_length 16 -outfmt 6", ""),
    ("dc.w11_coding_18", f"{DC} -word_size 11 -template_type coding -template_length 18 -outfmt 6", ""),
    ("dc.w11_coding_21", f"{DC} -word_size 11 -template_type coding -template_length 21 -outfmt 6", ""),
    ("dc.w11_optimal_16", f"{DC} -word_size 11 -template_type optimal -template_length 16 -outfmt 6", ""),
    ("dc.w11_optimal_18", f"{DC} -word_size 11 -template_type optimal -template_length 18 -outfmt 6", ""),
    ("dc.w11_optimal_21", f"{DC} -word_size 11 -template_type optimal -template_length 21 -outfmt 6", ""),
    ("dc.w11_coding_and_optimal_16", f"{DC} -word_size 11 -template_type coding_and_optimal -template_length 16 -outfmt 6", ""),
    ("dc.w11_coding_and_optimal_18", f"{DC} -word_size 11 -template_type coding_and_optimal -template_length 18 -outfmt 6", ""),
    ("dc.w11_coding_and_optimal_21", f"{DC} -word_size 11 -template_type coding_and_optimal -template_length 21 -outfmt 6", ""),
    ("dc.w12_coding_16", f"{DC} -word_size 12 -template_type coding -template_length 16 -outfmt 6", ""),
    ("dc.w12_coding_18", f"{DC} -word_size 12 -template_type coding -template_length 18 -outfmt 6", ""),
    ("dc.w12_coding_21", f"{DC} -word_size 12 -template_type coding -template_length 21 -outfmt 6", ""),
    ("dc.w12_optimal_16", f"{DC} -word_size 12 -template_type optimal -template_length 16 -outfmt 6", ""),
    ("dc.w12_optimal_18", f"{DC} -word_size 12 -template_type optimal -template_length 18 -outfmt 6", ""),
    ("dc.w12_optimal_21", f"{DC} -word_size 12 -template_type optimal -template_length 21 -outfmt 6", ""),
    ("dc.w12_coding_and_optimal_16", f"{DC} -word_size 12 -template_type coding_and_optimal -template_length 16 -outfmt 6", ""),
    ("dc.w12_coding_and_optimal_18", f"{DC} -word_size 12 -template_type coding_and_optimal -template_length 18 -outfmt 6", ""),
    ("dc.w12_coding_and_optimal_21", f"{DC} -word_size 12 -template_type coding_and_optimal -template_length 21 -outfmt 6", ""),
    ("dc.megablast_w11_coding_18", f"{DCQ} -word_size 11 -template_type coding -template_length 18 -outfmt 6", ""),
    ("dc.default_threads4", f"{DC} -outfmt 6", "-num_threads 4"),
    ("dc.two21_threads2_fmt7", f"{DC} -template_type coding_and_optimal -template_length 21 -outfmt 7",
     "-num_threads 2"),
    ("dc.small_two16", f"-query {O}/dc_small_query.fasta -subject {O}/dc_subject.fasta -task dc-megablast"
                       " -template_type coding_and_optimal -template_length 16 -outfmt 6", ""),
    # Lowercase in the subject with -lcase_masking: the discontiguous scan over the subject's
    # unmasked ranges (na_ungapped.c:1651-1667), and the control without the option.
    ("dc.lcase_subject", f"-query {I}/lcase_island_megablast_query.fa -subject {I}/lcase_island_megablast_subject.fa"
                         " -task dc-megablast -lcase_masking -outfmt 6", ""),
    ("dc.lcase_subject_control", f"-query {I}/lcase_island_megablast_query.fa"
                                 f" -subject {I}/lcase_island_megablast_subject.fa -task dc-megablast -outfmt 6", ""),
    ("dc.prelim", f"{P} -task dc-megablast -outfmt 6", ""),
    ("dc.prelim_max3_fmt0", f"{P} -task dc-megablast -max_target_seqs 3 -outfmt 0", ""),
    ("dc.batches", f"{B} -task dc-megablast -outfmt 6", ""),
    ("dc.batches_besthit", f"{B} -task dc-megablast -subject_besthit -outfmt 6", ""),
    ("dc.ambiguity", f"{A} -task dc-megablast -outfmt 6", ""),
    ("dc.rep", f"-query {I}/rep_query.fa -subject {I}/rep_subject.fa -task dc-megablast -outfmt 6", ""),
    # NCBI's query chunks are 5,000,000 for dc-megablast (local_blast.cpp:64-72): the 2 Mb
    # query is one chunk, as for megablast; with CHUNK_SIZE it is split.
    ("dc.split_unsplit", f"{S} -task dc-megablast -outfmt 6", ""),
    ("dc.split_chunk300000", f"{S} -task dc-megablast -outfmt 6", "", "CHUNK_SIZE=300000"),
    ("dc.batch1000_fmt7", f"{B} -task dc-megablast -outfmt 7", "", "BATCH_SIZE=1000"),
    # Option checks (blast_options.c:1247-1261,1399-1413): the word size of a template, the
    # lookup table of the task.
    ("dc.word13", f"{DC} -word_size 13 -outfmt 6", ""),
    ("dc.word13_fmt0", f"{DC} -word_size 13 -outfmt 0", ""),
    ("dc.blastn_template", f"{DCQ} -task blastn -word_size 11 -template_type coding -template_length 18 -outfmt 6", ""),
    ("dc.short_template", f"{DCQ} -task blastn-short -template_type optimal -template_length 21 -outfmt 6", ""),
    ("dc.megablast_template_word28", f"{DCQ} -template_type coding -template_length 18 -outfmt 6", ""),
    ("dc.zero_gaps", f"{DC} -gapopen 0 -gapextend 0 -outfmt 6", ""),
    # The divergent pair (dc_div_*.fasta, added after the SD audit): NCBI gives a different
    # output for each of the 18 combinations (the pair above gives 9), so each scanner of
    # blast_nascan.c:2202-2626 and the second template are told apart.
    ("dc.div_w11_coding_16", f"{DV} -word_size 11 -template_type coding -template_length 16 -outfmt 6", ""),
    ("dc.div_w11_coding_18", f"{DV} -word_size 11 -template_type coding -template_length 18 -outfmt 6", ""),
    ("dc.div_w11_coding_21", f"{DV} -word_size 11 -template_type coding -template_length 21 -outfmt 6", ""),
    ("dc.div_w11_optimal_16", f"{DV} -word_size 11 -template_type optimal -template_length 16 -outfmt 6", ""),
    ("dc.div_w11_optimal_18", f"{DV} -word_size 11 -template_type optimal -template_length 18 -outfmt 6", ""),
    ("dc.div_w11_optimal_21", f"{DV} -word_size 11 -template_type optimal -template_length 21 -outfmt 6", ""),
    ("dc.div_w11_coding_and_optimal_16",
     f"{DV} -word_size 11 -template_type coding_and_optimal -template_length 16 -outfmt 6", ""),
    ("dc.div_w11_coding_and_optimal_18",
     f"{DV} -word_size 11 -template_type coding_and_optimal -template_length 18 -outfmt 6", ""),
    ("dc.div_w11_coding_and_optimal_21",
     f"{DV} -word_size 11 -template_type coding_and_optimal -template_length 21 -outfmt 6", ""),
    ("dc.div_w12_coding_16", f"{DV} -word_size 12 -template_type coding -template_length 16 -outfmt 6", ""),
    ("dc.div_w12_coding_18", f"{DV} -word_size 12 -template_type coding -template_length 18 -outfmt 6", ""),
    ("dc.div_w12_coding_21", f"{DV} -word_size 12 -template_type coding -template_length 21 -outfmt 6", ""),
    ("dc.div_w12_optimal_16", f"{DV} -word_size 12 -template_type optimal -template_length 16 -outfmt 6", ""),
    ("dc.div_w12_optimal_18", f"{DV} -word_size 12 -template_type optimal -template_length 18 -outfmt 6", ""),
    ("dc.div_w12_optimal_21", f"{DV} -word_size 12 -template_type optimal -template_length 21 -outfmt 6", ""),
    ("dc.div_w12_coding_and_optimal_16",
     f"{DV} -word_size 12 -template_type coding_and_optimal -template_length 16 -outfmt 6", ""),
    ("dc.div_w12_coding_and_optimal_18",
     f"{DV} -word_size 12 -template_type coding_and_optimal -template_length 18 -outfmt 6", ""),
    ("dc.div_w12_coding_and_optimal_21",
     f"{DV} -word_size 12 -template_type coding_and_optimal -template_length 21 -outfmt 6", ""),
    ("dc.div_megablast_w11_coding_18", f"{DVQ} -word_size 11 -template_type coding -template_length 18 -outfmt 6",
     ""),
    ("dc.div_two21_fmt0", f"{DV} -template_type coding_and_optimal -template_length 21 -outfmt 0", ""),
    ("dc.div_w12_optimal_16_fmt7", f"{DV} -word_size 12 -template_type optimal -template_length 16 -outfmt 7", ""),
    ("dc.div_two18_threads4", f"{DV} -template_type coding_and_optimal -template_length 18 -outfmt 6",
     "-num_threads 4"),
    ("short.default", f"{SH} -outfmt 6", ""),
    ("short.default_threads4", f"{SH} -outfmt 6", "-num_threads 4"),
    ("short.word4", f"{SH} -word_size 4 -evalue 10 -outfmt 6", ""),
    ("short.word16", f"{SH} -word_size 16 -outfmt 6", ""),
    ("short.reward2_penalty3_gaps", f"{SH} -reward 2 -penalty -3 -gapopen 4 -gapextend 4 -outfmt 6", ""),
    ("short.evalue1e3_dust", f"{SH} -evalue 1e-3 -dust '15 64 1' -outfmt 6", ""),
    ("short.dust_no_fmt7", f"{SH} -dust no -outfmt 7", ""),
    ("short.lcase", f"{SH} -lcase_masking -outfmt 6", ""),
    ("short.perc_besthit", f"{SH} -perc_identity 90 -subject_besthit -outfmt 6", ""),
    ("short.batches", f"{B} -task blastn-short -evalue 1e-5 -outfmt 6", ""),
    ("short.prelim_max3", f"{P} -task blastn-short -evalue 1e-10 -max_target_seqs 3 -outfmt 6", ""),
    ("short.ambiguity", f"{A} -task blastn-short -outfmt 6", ""),
    # blastn-short is eBlastn: 1,000,000 query chunks, so the 2 Mb query is split.
    ("short.split", f"{S} -task blastn-short -evalue 1e-20 -outfmt 6", ""),
    ("short.nohit_fmt0", f"-query {O}/short_nohit_query.fasta -subject {O}/short_subject.fasta -task blastn-short"
                         " -evalue 1e-5 -outfmt 0", ""),
    ("short.evalue0", f"{SH} -evalue 0 -outfmt 6", ""),
    # M-32: the approximate ungapped extension (word size 11 or more) works on the concatenated
    # query, so its left loop runs through the sentinel into the previous context
    # (na_ungapped.c:286-305); a result kept below the reduced cutoff may start there and is
    # always rejected, and the score of the block that crosses can decide the exact rerun
    # (na_ungapped.c:336-341). The queries are the ends of pieces of LC738874 (qctx_*.fa, from
    # `generate`); each hit starts at the first letter of the minus strand.
    # A rejected result that starts two letters left of the minus strand (search_query_batch).
    ("qctx.assert", f"-query {I}/qctx_assert_query.fa -subject tests/fasta/LC738868.fasta -task blastn -outfmt 6",
     ""),
    # The crossing block lifts the approximate score to the reduced cutoff; the exact rerun
    # then reaches the cutoff and the 11-letter HSP at query 50-60 is saved.
    ("qctx.rerun", f"-query {I}/qctx_rerun_query.fa -subject tests/fasta/LC738868.fasta -task blastn -outfmt 6", ""),
    # The same for the 14-letter HSP at query 2071-2084.
    ("qctx.rerun_piece8", f"-query {I}/qctx_piece8_query.fa -subject tests/fasta/LC738868.fasta -task blastn"
                          " -outfmt 6", ""),
    # A query longer than 8000 at four threads: the search pool's seed stage (x_pair_seed.rs).
    ("qctx.pairpar", f"-query {I}/qctx_pairpar_query.fa -subject tests/fasta/LC738868.fasta -task blastn -outfmt 6",
     "-num_threads 4"),
    # O3a (blastn -culling_limit; hspfilter_culling.c): the preliminary writer culls each subject's HSP
    # list, the culling pipe culls again after the traceback. culling_ties_* / culling_multi_* are the
    # inputs of tblastx_regression's culling.* cases (tests/fixtures/tblastx_regression/inputs).
    # Limits 1, 2 and 5 on overlapping HSPs (rep_*); the lowest OID wins among identical subjects, so the
    # record order decides the kept copy (culling_*_subject_rev.fa are the subjects in reverse order, made
    # by `generate`); a palindromic query folds both strands into one tree.
    ("culling.rep_megablast_cl1", f"{RP} -culling_limit 1 -outfmt 6", ""),
    ("culling.rep_blastn_cl2_fmt0", f"{RP} -task blastn -culling_limit 2 -outfmt 0", ""),
    ("culling.rep_blastn_cl5", f"{RP} -task blastn -culling_limit 5 -outfmt 6", ""),
    ("culling.ties_megablast_cl1_fmt7", f"{CT} -culling_limit 1 -outfmt 7", ""),
    ("culling.ties_rev_megablast_cl1", f"{CTR} -culling_limit 1 -outfmt 6", ""),
    ("culling.multi_blastn_cl1", f"{CM} -task blastn -culling_limit 1 -outfmt 6", ""),
    ("culling.multi_rev_blastn_cl1", f"{CMR} -task blastn -culling_limit 1 -outfmt 6", ""),
    ("culling.pal_blastn_cl1", f"{PL} -task blastn -culling_limit 1 -outfmt 6", ""),
    # A split query: the 40 queries of mq.fa are cut into chunks of 2000 (the chunks' HSPs are
    # joined at the edges as without culling; each chunk is culled by itself).
    ("culling.split_mq_blastn_cl2_c2000", f"{B} -task blastn -culling_limit 2 -outfmt 6", "", "CHUNK_SIZE=2000"),
    ("culling.split_mq_megablast_cl1_c2000", f"{B} -culling_limit 1 -outfmt 6", "", "CHUNK_SIZE=2000"),
    # NCBI adds 3 to a limit above 1 in an Int4 (setup_factory.cpp:338-340), which overflows
    # from 2147483645 on: its binary wraps 2147483645 to INT_MIN (culls nothing) and
    # 2147483646 and 2147483647 to negative merits (a limit of 1 in the preliminary stage).
    # LOSAT adds in i64 (approved exception, PD-LOSAT-NCBI-DEFECTS 1.4), so these limits cull
    # nothing, as 2147483644 does: the oracle of the last two runs 2147483644 (ORACLE_ARGV).
    # 0x2 is limit 2 (CArg_Int8 reads C-style bases).
    ("culling.near_int4_max_2147483645", f"{CT} -culling_limit 2147483645 -outfmt 6", ""),
    ("culling.near_int4_max_2147483646", f"{CT} -culling_limit 2147483646 -outfmt 6", ""),
    ("culling.near_int4_max_2147483647", f"{CT} -culling_limit 2147483647 -outfmt 6", ""),
    ("culling.hex_0x2", f"{CT} -culling_limit 0x2 -outfmt 6", ""),
    # -max_target_seqs cuts the traceback hit list before the culling pipe: the preliminary cull
    # removes the other copies first, so -max_target_seqs 1 keeps the first; with 3 and limit 2
    # the pipe culls what the cut left.
    ("culling.ties_mts1_cl1", f"{CT} -max_target_seqs 1 -culling_limit 1 -outfmt 6", ""),
    ("culling.ties_mts3_cl2", f"{CT} -max_target_seqs 3 -culling_limit 2 -outfmt 6", ""),
    # ABRicate's search (2 subsets of the O3a probe inputs: 4 windows of 50 kb of Sakai with the most
    # hits of the 340 records of its genes database, and the 33 records they hit; copied from the
    # probe inputs, not made by `generate`), also with a query split at 20000.
    ("culling.abricate_cl1", f"{AB} -culling_limit 1 -outfmt 6", ""),
    ("culling.abricate_cl1_c20000", f"{AB} -culling_limit 1 -outfmt 6", "", "CHUNK_SIZE=20000"),
    # O2a (blastn -strand): a query split into chunks searches one strand's contexts in every chunk;
    # DUST of -strand minus is that of the reverse complement; an invalid query warns once.
    # strand_dust_*.fa (a query whose low-complexity region lies on one strand only) and
    # strand_mixed_query.fa (all-N and empty records around strand_query.fasta) are copied from the
    # O2a probe inputs, not made by `generate`.
    ("strand.split_plus_c2000", f"{B} -task blastn -strand plus -outfmt 6", "", "CHUNK_SIZE=2000"),
    ("strand.split_minus_c2000", f"{B} -task blastn -strand minus -outfmt 6", "", "CHUNK_SIZE=2000"),
    ("strand.split_minus_c1500o50_fmt0", f"{B} -task blastn -strand minus -outfmt 0", "",
     "CHUNK_SIZE=1500 OVERLAP_CHUNK_SIZE=50"),
    ("strand.dust_blastn_minus_fmt0", f"-query {I}/strand_dust_query.fa -subject {I}/strand_dust_subject.fa"
                                      " -task blastn -strand minus -outfmt 0", ""),
    ("strand.dust_megablast_minus_fmt0", f"-query {I}/strand_dust_query.fa -subject {I}/strand_dust_subject.fa"
                                         " -strand minus -outfmt 0", ""),
    ("strand.invalid_mixed_minus_fmt0", f"-query {I}/strand_mixed_query.fa -subject {O}/strand_subject.fasta"
                                        " -task blastn -strand minus -outfmt 0", ""),
    ("strand.lcase_qloc_plus_fmt0", f"-query {I}/strand_dust_query.fa -subject {I}/strand_dust_subject.fa"
                                    " -task blastn -strand plus -lcase_masking -query_loc 20-700 -outfmt 0", ""),
    # NP M2 (blastn -window_masker_db, counts tests/winmask/counts/; wm_*.fa made by `generate_wm_inputs`):
    # WindowMasker's masks replace those of DUST and -lcase_masking (the .lc_dust rows equal the rows
    # without them), in the four tasks and with the binary counts; -strand minus masks with the intervals
    # of the reverse complement; -query_loc masks the range's letters; each query chunk and batch masks
    # its own letters; the reports show the masks (outfmt 0 lower case, XML filter `W -d <file>;`, the
    # archive's WindowMaskerDatabase and LCaseMask).
    ("wm.megablast_lc_dust", f"{WM} -window_masker_db {WMC} -lcase_masking -dust yes -outfmt 6", ""),
    ("wm.blastn", f"{WM} -task blastn -window_masker_db {WMC} -outfmt 6", ""),
    ("wm.blastn_lc_dust", f"{WM} -task blastn -window_masker_db {WMC} -lcase_masking -dust yes -outfmt 6", ""),
    ("wm.blastn_nodust", f"{WM} -task blastn -window_masker_db {WMC} -dust no -outfmt 6", ""),
    ("wm.dc", f"{WM} -task dc-megablast -window_masker_db {WMC} -outfmt 6", ""),
    ("wm.short", f"{WM3} -task blastn-short -window_masker_db {WMC} -evalue 1 -outfmt 6", ""),
    ("wm.binary", f"{WM} -task blastn -window_masker_db {WMB} -outfmt 6", ""),
    ("wm.blastn_lc_fmt0", f"{WM} -task blastn -window_masker_db {WMC} -lcase_masking -max_hsps 2 -outfmt 0", ""),
    ("wm.plus", f"{WM} -task blastn -window_masker_db {WMC} -strand plus -outfmt 6", ""),
    ("wm.minus", f"{WM} -task blastn -window_masker_db {WMC} -strand minus -outfmt 6", ""),
    ("wm.minus_fmt0", f"{WM} -task blastn -window_masker_db {WMC} -strand minus -max_hsps 2 -outfmt 0", ""),
    ("wm.qloc", f"{P} -task blastn -window_masker_db {WMC} -query_loc 101-2500 -outfmt 6", ""),
    ("wm.qloc_fmt0", f"{P} -task blastn -window_masker_db {WMC} -query_loc 101-2500 -max_target_seqs 5 -outfmt 0", ""),
    ("wm.qloc_minus", f"{P} -task blastn -window_masker_db {WMC} -query_loc 101-2500 -strand minus -outfmt 6", ""),
    # (Chunks of 600 and 400: with the chunks' own masks the hits are NCBI's, with the whole queries'
    # masks cut to the chunks they are not.)
    ("wm.split_c600", f"{WM} -task blastn -window_masker_db {WMC} -outfmt 6", "", "CHUNK_SIZE=600"),
    ("wm.split_minus_c400", f"{WM} -task blastn -window_masker_db {WMC} -strand minus -outfmt 6", "",
     "CHUNK_SIZE=400"),
    ("wm.split_c1500o50", f"{WM} -task blastn -window_masker_db {WMC} -outfmt 6", "",
     "CHUNK_SIZE=1500 OVERLAP_CHUNK_SIZE=50"),
    ("wm.batch3000", f"{WM} -task blastn -window_masker_db {WMC} -outfmt 6", "", "BATCH_SIZE=3000"),
    ("wm.mixed", f"-query {I}/wm_mixed_query.fa -subject {I}/s60k.fa -task blastn -window_masker_db {WMC} -outfmt 6", ""),
    ("wm.fmt11", f"{WM3} -task blastn -window_masker_db {WMC} -outfmt 11", ""),
    ("wm.minus_fmt11", f"{WM3} -task blastn -window_masker_db {WMC} -strand minus -outfmt 11", ""),
    ("wm.plus_lc_fmt11", f"{WM3} -task blastn -window_masker_db {WMC} -strand plus -lcase_masking -outfmt 11", ""),
    ("wm.fmt5", f"{WM3} -task blastn -window_masker_db {WMC} -outfmt 5", ""),
    ("wm.fmt16", f"{WM3} -task blastn -window_masker_db {WMC} -outfmt 16", ""),
    ("wm.fmt15", f"{WM3} -task blastn -window_masker_db {WMC} -outfmt 15", ""),
]


def _ungapped_cases():
    """Line NP O4a: `-ungapped` (NCBI blast_args.cpp:917,929) in the four tasks, with the options that
    change the preliminary stage's result list (culling, -max_hsps, -qcov_hsp_perc, -strand, -perc_identity,
    -subject_besthit), scores without Karlin tables, gap costs, the query split (CHUNK_SIZE), the
    ambiguity re-evaluation, and the formats that O4a-1 writes (0, 6, 7)."""
    tasks = (("mb", ""), ("bn", "-task blastn "), ("dc", "-task dc-megablast "), ("short", "-task blastn-short "))
    cases = []
    for key, task in tasks:
        for name, inputs in (("t2", T2), ("p", P), ("mq", B), ("amb", A), ("ambrun", AMB), ("reeval", RE)):
            # blastn-short with the 60 kb subject and many queries has hundreds of thousands of HSPs at the
            # default E-value (2 MB of output); 1e-3 keeps the batch and the many hits of a short query.
            small = " -evalue 1e-3" if (key, name) == ("short", "mq") else ""
            cases.append((f"ungapped.{key}.{name}", f"{inputs} {task}-ungapped{small} -outfmt 6", ""))
    for key, task in tasks[:2]:
        for name, option in (("cl1", "-culling_limit 1"), ("cl2", "-culling_limit 2"), ("mh1", "-max_hsps 1"),
                             ("qcov50", "-qcov_hsp_perc 50"), ("plus", "-strand plus"), ("minus", "-strand minus"),
                             ("pid95", "-perc_identity 95"), ("sbh", "-subject_besthit"),
                             ("mts1", "-max_target_seqs 1")):
            cases.append((f"ungapped.{key}.rep_{name}", f"{RP} {task}-ungapped {option} -outfmt 6", ""))
        for name, option in (("cl1", "-culling_limit 1"), ("mh1", "-max_hsps 1"), ("qcov20", "-qcov_hsp_perc 20"),
                             ("minus", "-strand minus"), ("sbh", "-subject_besthit")):
            cases.append((f"ungapped.{key}.mix_{name}", f"{MX} {task}-ungapped {option} -outfmt 6", ""))
        cases.append((f"ungapped.{key}.t2_fmt0", f"{T2} {task}-ungapped -outfmt 0", ""))
        cases.append((f"ungapped.{key}.mq_fmt7", f"{B} {task}-ungapped -outfmt 7", ""))
        cases.append((f"ungapped.{key}.mq_c2000", f"{B} {task}-ungapped -outfmt 6", "", "CHUNK_SIZE=2000"))
    # Ties (review M1): a query X + revcomp(X) has plus and minus HSPs of the same score and range; the
    # order of those comes from the initial hits' sort by score (na_ungapped.c:1691-1692).
    cases.append(("ungapped.mb.pal", f"{PL} -ungapped -outfmt 6", ""))
    cases.append(("ungapped.dc.pal", f"{PL} -task dc-megablast -ungapped -outfmt 6", ""))
    # A subject above 5 Mb is searched in chunks that are joined without a gap (review S3): a 600 bp piece of
    # EDL933 at 4,999,701 with a one base deletion and a mismatch, against the whole EDL933 (5,528,445 bp);
    # the HSP crosses the chunk border at 5,000,000 (made at run time, prepare_runtime_inputs).
    for key, task in (("mb", ""), ("bn", "-task blastn "), ("dc", "-task dc-megablast ")):
        cases.append((f"ungapped.{key}.edl", f"{EDLU} {task}-ungapped -outfmt 6", ""))
    cases.append(("ungapped.mb.mix_qloc", f"{MX} -ungapped -query_loc 150-750 -outfmt 6", ""))
    cases.append(("ungapped.mb.mix_sloc", f"{MX} -ungapped -subject_loc 50-700 -outfmt 6", ""))
    for name, reward, penalty in (("1_6", 1, -6), ("3_1", 3, -1), ("2_7", 2, -7), ("5_4", 5, -4),
                                  ("1_1", 1, -1), ("2_2", 2, -2)):
        cases.append((f"ungapped.bn.t2_rp{name}",
                      f"{T2} -task blastn -ungapped -reward {reward} -penalty {penalty} -outfmt 6", ""))
    cases.append(("ungapped.bn.t2_go5ge0", f"{T2} -task blastn -ungapped -gapopen 5 -gapextend 0 -outfmt 6", ""))
    cases.append(("ungapped.mb.t2_go40000",
                  f"{T2} -task megablast -ungapped -gapopen 40000 -gapextend 40000 -outfmt 6", ""))
    cases.append(("ungapped.bn.t2_go0ge0", f"{T2} -task blastn -ungapped -gapopen 0 -gapextend 0 -outfmt 6", ""))
    cases.append(("ungapped.bn.inval_fmt7", f"{SM} -task blastn -ungapped -outfmt 7", ""))
    cases.append(("ungapped.bn.alln_fmt0", f"-query {I}/ungapped_alln_query.fa -subject {I}/s60k.fa -task blastn"
                                           " -ungapped -outfmt 0", ""))
    cases.append(("ungapped.mlst", f"{MLST} -outfmt 6", ""))
    cases.append(("ungapped.mlst_fmt0", f"{MLST} -outfmt 0", ""))
    # The custom tabular reports (line CLI, C2) of an ungapped search: mlst's real call, and the columns whose values an
    # ungapped alignment may set differently from a gapped one (an HSP without an edit script is one run of
    # substitutions: gapopen 0, gaps 0, btop and the identity counts from the re-evaluated alignment).
    cases.append(("ungapped.mlst_tab", f"{MLST} -outfmt '6 sseqid slen length nident qseqid qstart qend qseq sstrand'",
                  ""))
    tab_cols = "qseqid sseqid qstart qend sstart send gapopen gaps btop qcovhsp qcovs nident positive ppos pident" \
               " length mismatch qseq sseq"
    cases.append(("ungapped.tab_cols_t2", f"{T2} -ungapped -outfmt '6 {tab_cols}'", ""))
    cases.append(("ungapped.tab_cols_mix", f"{MX} -task blastn -ungapped -outfmt '7 {tab_cols}'", ""))
    cases.append(("ungapped.tab_cols_ambrun", f"{AMB} -task blastn -ungapped -outfmt '6 {tab_cols}'", ""))
    cases.append(("ungapped.tab_csv20", f"{MX} -task blastn -ungapped -outfmt 20", ""))
    cases.append(("ungapped.tab_csv10", f"{MX} -task blastn -ungapped -outfmt '10 qseqid sseqid btop gaps nident'", ""))
    for name, argv in (("mb.t2", f"{T2} -ungapped -outfmt 6"), ("bn.p", f"{P} -task blastn -ungapped -outfmt 6"),
                       ("mb.mq", f"{B} -ungapped -outfmt 6"), ("dc.mq", f"{B} -task dc-megablast -ungapped -outfmt 6")):
        cases.append((f"ungapped.{name}_threads4", argv, "-num_threads 4"))
    return cases


_CASES += _ungapped_cases()
CASES = [(*case, *[""] * (5 - len(case))) for case in _CASES]
# NCBI's arguments of the cases of an approved exception that is an argument (see above).
ORACLE_ARGV = {
    "culling.near_int4_max_2147483646": f"{CT} -culling_limit 2147483644 -outfmt 6",
    "culling.near_int4_max_2147483647": f"{CT} -culling_limit 2147483644 -outfmt 6",
}


def genome(name: str) -> str:
    return "".join(line.strip() for line in (FASTA / name).read_text().splitlines()[1:])


def write_fasta(path: Path, records) -> None:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))


def mutate(rng: random.Random, seq: str, subst: float, indel: float = 0.0) -> str:
    out = []
    for base in seq:
        r = rng.random()
        if r < indel * 0.5:
            continue
        if r < indel:
            out += [base, rng.choice("ACGT")]
        elif rng.random() < subst:
            out.append(rng.choice([c for c in "ACGT" if c != base]))
        else:
            out.append(base)
    return "".join(out)


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def command_generate(_args) -> int:
    INPUTS.mkdir(parents=True, exist_ok=True)
    generate_culling_inputs()
    rng = random.Random(20261002)
    edl, sakai = genome("EDL933.fna"), genome("Sakai.fna")
    q3k = edl[1_000_000:1_003_000]
    write_fasta(INPUTS / "q3k.fa", [("q3k EDL933 1000001-1003000", q3k)])
    records = []
    for index in range(600):
        length = rng.choice((60, 80, 100, 150, 250))
        start = rng.randrange(0, len(q3k) - length)
        window = mutate(rng, q3k[start:start + length], rng.choice((0.0, 0.0, 0.01, 0.03)))
        if rng.random() < 0.3:
            window = revcomp(window)
        left = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 200)))
        right = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 200)))
        records.append((f"s{index:03d}", left + window + right))
    write_fasta(INPUTS / "s600.fa", records)
    records = []
    for index in range(600):
        start = rng.randrange(0, len(q3k) - 100)
        records.append((f"t{index:03d}", q3k[start:start + 100]))
    write_fasta(INPUTS / "s600ties.fa", records)
    region = sakai[2_000_000:2_060_000]
    write_fasta(INPUTS / "s60k.fa", [("s60k Sakai 2000001-2060000 mutated", mutate(rng, region, 0.01, 0.005))])
    records = []
    for index in range(40):
        length = int(10 * (600 ** rng.random()))
        start = rng.randrange(0, len(region) - length)
        window = mutate(rng, region[start:start + length], rng.choice((0.0, 0.01, 0.03)))
        if rng.random() < 0.3:
            window = revcomp(window)
        records.append((f"mq{index:02d} len{length}", window))
    write_fasta(INPUTS / "mq.fa", records)
    write_fasta(INPUTS / "sakai_60k.fa", [("sakai_60k Sakai 500001-560000", sakai[500_000:560_000])])
    for name in ("ambiguity_query_c.fa", "ambiguity_subject.fa", "lcase_island_blastn_query.fa",
                 "lcase_island_blastn_subject.fa", "lcase_island_megablast_query.fa",
                 "lcase_island_megablast_subject.fa"):
        (INPUTS / name).write_bytes((E2C_INPUTS / name).read_bytes())
    generate_e2g_inputs(edl, sakai)
    generate_qctx_inputs()
    generate_wm_inputs()
    generate_ungapped_inputs()
    return 0


def generate_e2g_inputs(edl: str, sakai: str) -> None:
    """Inputs added in E2g; a separate generator keeps the earlier inputs unchanged."""
    rng = random.Random(20261003)
    # T1: query X + revcomp(X); its minus-strand context is the same sequence, so both
    # strands give hits with equal score, subject start and length (context-local q_start).
    records, subject = [], []
    for index in range(40):
        start = rng.randrange(0, len(edl) - 400)
        x = edl[start:start + rng.choice((60, 90, 120, 150))]
        records.append((f"pal{index:02d} EDL933 {start + 1} X+revcomp(X)", x + revcomp(x)))
        subject.append("".join(rng.choice("ACGT") for _ in range(rng.randint(20, 120))))
        subject.append(mutate(rng, x, rng.choice((0.0, 0.02))))
    write_fasta(INPUTS / "pal_query.fa", records)
    write_fasta(INPUTS / "pal_subject.fa", [("pals EDL933 windows", "".join(subject))])
    # T6: several short queries whose block (both strands) is at most 8000, so NCBI
    # uses the diagonal array (blast_parameters.c:225-231); the subject repeats them.
    region = sakai[3_000_000:3_040_000]
    records, subject = [], []
    for index in range(12):
        length = rng.randint(100, 300)
        start = rng.randrange(0, len(region) - length)
        window = region[start:start + length]
        records.append((f"sq{index:02d} Sakai {3_000_001 + start} len{length}", mutate(rng, window, 0.02)))
        for _ in range(rng.randint(1, 3)):
            subject.append("".join(rng.choice("ACGT") for _ in range(rng.randint(30, 400))))
            copy = mutate(rng, window, rng.choice((0.0, 0.03, 0.08)), 0.01)
            subject.append(revcomp(copy) if rng.random() < 0.4 else copy)
    write_fasta(INPUTS / "sq_query.fa", records)
    write_fasta(INPUTS / "sq_subject.fa", [("sqs Sakai windows repeated", "".join(subject))])
    # T5: queries just above the 8000 block limit (diagonal hash) that still get the
    # small lookup table (blast_nalookup.c:45-185), with a repeated unit so one
    # subject offset gives seeds at several query offsets.
    region = sakai[4_000_000:4_100_000]
    unit = region[50_000:50_250]
    parts, used = [], 0
    while used < 4100 - 250:
        piece = region[used:used + rng.randint(300, 700)]
        parts.append(piece)
        parts.append(mutate(rng, unit, rng.choice((0.0, 0.01, 0.03))))
        used += len(piece) + 250
    rep = "".join(parts)[:4100]
    second = mutate(rng, region[60_000:60_700], 0.02) + unit
    write_fasta(INPUTS / "rep_query.fa", [("rep Sakai 4000001 with a repeated unit", rep)])
    write_fasta(INPUTS / "rep2_query.fa", [("rep Sakai 4000001 with a repeated unit", rep),
                                           ("rep2 Sakai 4060001 and the unit", second)])
    subject = []
    for index in range(10):
        subject.append("".join(rng.choice("ACGT") for _ in range(rng.randint(50, 400))))
        copy = mutate(rng, unit if index % 2 == 0 else rep[index * 300:index * 300 + 600], 0.02, 0.005)
        subject.append(revcomp(copy) if index % 3 == 0 else copy)
    write_fasta(INPUTS / "rep_subject.fa", [("reps Sakai unit copies", "".join(subject))])
    # T9: two query batches (the first ends at 5000 residues), each with a title that ends
    # in nucleotides (CFastaReader's warning) and an invalid query (all N).
    window = sakai[3_100_000:3_120_000]
    write_fasta(INPUTS / "warn_query.fa", [
        ("q1 first ACGTACGTACGTACGTACGTACGTAC", window[1000:4000]), ("q2 all N one", "N" * 2500),
        ("q3 third GGGGCCCCAAAATTTTGGGGCCCCAA", window[8000:8300]), ("q4 all N two", "N" * 100),
        ("q5 plain", window[12000:12500])])
    write_fasta(INPUTS / "warn_subject.fa", [("s Sakai window", window)])


def generate_ungapped_inputs() -> None:
    """Inputs added for line NP O4a (the ungapped.* cases): the probes of O/O4a/probe/d/gen_inputs.py with the same
    seed. ambrun: a subject with 30 N and one R inside a matching stretch, which the ungapped re-evaluation
    splits into two HSPs; mix: a query with plus and minus HSPs on two subjects; alln: a query of N only."""
    rng = random.Random(4242)

    def random_seq(n: int) -> str:
        return "".join(rng.choice("ACGT") for _ in range(n))

    core = random_seq(600)
    s = random_seq(200) + core[:290] + "N" * 30 + core[320:450] + "R" + core[451:] + random_seq(200)
    q = random_seq(50) + core + random_seq(50)
    write_fasta(INPUTS / "ambrun_q.fa", [("qamb", q)])
    write_fasta(INPUTS / "ambrun_s.fa", [("samb", s)])
    q2 = random_seq(1000)
    s2 = (random_seq(100) + q2[100:400] + random_seq(100) + revcomp(q2[500:700]) + random_seq(100) + q2[750:900]
          + random_seq(50) + q2[120:260] + random_seq(100))
    write_fasta(INPUTS / "mix_q.fa", [("qmix", q2)])
    write_fasta(INPUTS / "mix_s.fa", [("smix first", s2),
                                      ("smix2 second", revcomp(q2[300:600]) + random_seq(80) + q2[600:800])])
    write_fasta(INPUTS / "ungapped_alln_query.fa", [("alln1 all N", "N" * 24), ("alln2 all N", "N" * 40)])


def generate_culling_inputs() -> None:
    """The subjects of tblastx_regression's culling_ties and culling_multi inputs in reverse record order
    (O3a: the order of the subjects decides which of equal HSPs a culling limit keeps)."""
    for name in ("ties", "multi"):
        lines = (ENGINE / "tests/fixtures/tblastx_regression/inputs" / f"culling_{name}_subject.fa").read_text().splitlines()
        starts = [index for index, line in enumerate(lines) if line.startswith(">")] + [len(lines)]
        records = [lines[a:b] for a, b in zip(starts, starts[1:])]
        (INPUTS / f"culling_{name}_subject_rev.fa").write_text(
            "".join(line + "\n" for record in records[::-1] for line in record))


def generate_qctx_inputs() -> None:
    """Inputs added for M-32 (the qctx.* cases): ends of the pieces of LC738874 that
    LOSAT/tests/blastn_qbatchpar.rs searched (removed on 2026-10-28 with the QBATCHPAR slots;
    piece k starts at 8000 k, 1500 + 373 k % 2400 long)."""
    lc738874 = genome("LC738874.fasta")
    for name, start, end in (("qctx_assert", 82_530, 82_830), ("qctx_rerun", 82_770, 82_830),
                             ("qctx_piece8", 64_000, 66_084), ("qctx_pairpar", 72_530, 82_830)):
        write_fasta(INPUTS / f"{name}_query.fa", [(f"{name} LC738874 {start + 1}-{end}", lc738874[start:end])])


def generate_wm_inputs() -> None:
    """Inputs added for line NP M2 (the wm.* cases, blastn -window_masker_db): records of mq.fa with a run of
    lowercase letters, a DUST region (AT repeated) and, in one record, N and R letters, so the masks of
    WindowMasker (counts LOSAT/tests/winmask/counts/), DUST and -lcase_masking lie in different places."""
    records = []
    for line in (INPUTS / "mq.fa").read_text().splitlines():
        if line.startswith(">"):
            records.append([line[1:], ""])
        else:
            records[-1][1] += line
    out = []
    for index in (0, 3, 4, 9, 13, 2, 5, 22, 35):
        name, seq = records[index]
        n = len(seq)
        if n >= 600:
            seq = seq[:n // 3] + seq[n // 3:n // 3 + 150].lower() + seq[n // 3 + 150:2 * n // 3] + "AT" * 32 + seq[2 * n // 3:]
        if index == 9:
            seq = seq[:n // 4] + "N" * 30 + seq[n // 4:n // 4 + 40] + "R" + seq[n // 4 + 40:]
        out.append((f"wm_{name}", seq))
    write_fasta(INPUTS / "wm_query.fa", out)
    # Three of them and the record of 77 letters, for the reports that write every query and alignment;
    # and a record without letters and one of N only around the first (a batch of every kind of query).
    write_fasta(INPUTS / "wm_query3.fa", [out[0], out[1], out[2], out[5]])
    (INPUTS / "wm_mixed_query.fa").write_text(
        ">wm_empty no letters\n>wm_alln only N\n" + "N" * 50 + "\n" + f">{out[0][0]}\n{out[0][1]}\n")


def prepare_runtime_inputs() -> None:
    RUNTIME.mkdir(parents=True, exist_ok=True)
    path = RUNTIME / "edl933_2m.fa"
    expected = ">edl933_2m EDL933 1-2000000\n" + genome("EDL933.fna")[:2_000_000] + "\n"
    if not path.exists() or path.read_text() != expected:
        path.write_text(expected)
    # The query of the ungapped.*.edl cases: EDL933 4,999,701-5,000,300 without the base at 4,999,821 and
    # with a transition at the 401st base of what is left.
    window = genome("EDL933.fna")[4_999_700:5_000_300]
    piece = window[:120] + window[121:]
    piece = piece[:400] + {"A": "G", "G": "A", "C": "T", "T": "C"}[piece[400]] + piece[401:]
    path = RUNTIME / "edl_ungapped_q.fa"
    expected = ">edlq EDL933 4999701-5000300 with a deletion and a mismatch\n" + piece + "\n"
    if not path.exists() or path.read_text() != expected:
        path.write_text(expected)


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def clean_env(case_env: str = "") -> dict[str, str]:
    # LOSAT_X_REACH is kept (LOSAT/src/utils/xreach.rs): it records the switch paths and changes no output.
    env = {key: value for key, value in os.environ.items()
           if key == "LOSAT_X_REACH"
           or (key not in REPORT_ENV and not key.startswith(("LOSAT_", "RAYON_")))}
    for item in shlex.split(case_env):
        key, _, value = item.partition("=")
        env[key] = value
    return env


def read_manifest() -> list[dict[str, str]]:
    lines = [line for line in MANIFEST.read_text().splitlines() if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def command_freeze(args) -> int:
    if any(key in os.environ for key in REPORT_ENV) or (Path.home() / ".ncbirc").exists():
        raise SystemExit(f"unset {REPORT_ENV} and remove ~/.ncbirc before freezing")
    prepare_runtime_inputs()
    oracle = str(Path(args.oracle).resolve())
    version = subprocess.run([oracle, "-version"], capture_output=True, text=True).stdout.strip().replace("\n", "; ")
    rows = []
    kept = {}
    if args.only:
        # Only the cases whose id starts with the prefix (or one of the comma-separated prefixes) are run again;
        # the other rows and files stay.
        kept = {row["case_id"]: row for row in read_manifest()}
    for case_id, argv, extra, case_env, oracle_env in CASES:
        if args.only and not case_id.startswith(tuple(args.only.split(","))):
            if case_id in kept:
                rows.append(kept[case_id])
            continue
        result = run_case(oracle, [], case_id, ORACLE_ARGV.get(case_id, argv), oracle_env or case_env)
        (FIXTURES / f"{case_id}.out").write_bytes(result.stdout)
        err = FIXTURES / f"{case_id}.err"
        if result.stderr:
            err.write_bytes(result.stderr)
        elif err.exists():
            err.unlink()
        rows.append({"case_id": case_id, "argv": argv, "losat_extra": extra, "exit": str(result.returncode),
                     "stdout_sha256": sha256(result.stdout), "stdout_bytes": str(len(result.stdout)),
                     "stderr_sha256": sha256(result.stderr) if result.stderr else "", "env": case_env,
                     "oracle_env": oracle_env})
        print(f"{case_id}\texit {result.returncode}\t{result.stdout.count(b'\n')} lines\t{len(result.stderr)} stderr bytes", flush=True)
    with open(MANIFEST, "w", newline="") as handle:
        handle.write(f"# Frozen BLASTN outputs of {version} (comparison oracle only; see blastn_regression_fixtures.py).\n")
        handle.write("# Run from LOSAT/: env <env> blastn <argv> (LOSAT: losat blastn <argv> <losat_extra>). Written by `freeze`; do not edit.\n")
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return 0


def run_case(program: str, prefix: list[str], case_id: str, argv: str, case_env: str, extra: str = ""):
    command = [program, *prefix, *shlex.split(argv), *shlex.split(extra)]
    if case_id.endswith(".merged"):
        result = subprocess.run(command, cwd=ENGINE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                env=clean_env(case_env))
        result.stderr = b""
        return result
    return subprocess.run(command, cwd=ENGINE, capture_output=True, env=clean_env(case_env))


def check_one(losat: str, row: dict[str, str], threads: int | None = None) -> tuple[str, str]:
    extra = row["losat_extra"]
    if threads and "-num_threads" not in extra:
        extra = f"{extra} -num_threads {threads}".strip()
    result = run_case(losat, ["blastn"], row["case_id"], row["argv"], row.get("env") or "", extra)
    expected_out = (FIXTURES / f"{row['case_id']}.out").read_bytes()
    err_path = FIXTURES / f"{row['case_id']}.err"
    expected_err = err_path.read_bytes() if err_path.exists() else b""
    problems = []
    if sha256(expected_out) != row["stdout_sha256"]:
        problems.append("frozen stdout file does not match the manifest")
    if str(result.returncode) != row["exit"]:
        problems.append(f"exit {result.returncode}, expected {row['exit']}")
    if result.stdout != expected_out:
        expected_lines, actual_lines = expected_out.split(b"\n"), result.stdout.split(b"\n")
        line = next((n for n, (a, b) in enumerate(zip(expected_lines, actual_lines), 1) if a != b),
                    min(len(expected_lines), len(actual_lines)))
        problems.append(f"stdout differs (expected {len(expected_lines) - 1} lines, got {len(actual_lines) - 1}; first difference at line {line})")
    if result.stderr != expected_err:
        problems.append(f"stderr differs: {result.stderr[:200]!r}")
    return row["case_id"], "; ".join(problems) or "same"


def command_check(args) -> int:
    prepare_runtime_inputs()
    losat = str(Path(args.losat).resolve())
    rows = [row for row in read_manifest() if not args.only or row["case_id"].startswith(tuple(args.only.split(",")))]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(lambda row: check_one(losat, row, args.threads), rows))
    lines = [f"{case_id}\t{verdict}" for case_id, verdict in results]
    if args.out:
        Path(args.out).write_text("case_id\tresult\n" + "\n".join(lines) + "\n")
    print("\n".join(lines))
    differing = [case_id for case_id, verdict in results if verdict != "same"]
    print(f"{len(results)} cases, {len(differing)} differ")
    return 1 if differing else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("generate")
    freeze = sub.add_parser("freeze")
    freeze.add_argument("--oracle", required=True)
    freeze.add_argument("--only", help="run only the cases whose id starts with this prefix (new cases)")
    check = sub.add_parser("check")
    check.add_argument("--losat", required=True)
    check.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    check.add_argument("--out")
    check.add_argument("--threads", type=int, help="LOSAT's -num_threads for the cases that give none")
    check.add_argument("--only", help="check only the cases whose id starts with this prefix")
    args = parser.parse_args()
    return {"generate": command_generate, "freeze": command_freeze, "check": command_check}[args.action](args)


if __name__ == "__main__":
    sys.exit(main())
