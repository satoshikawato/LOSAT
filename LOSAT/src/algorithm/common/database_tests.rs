// NCBI BLAST+ 2.17.0 (review 1 of line DBID ordname, `R/ordname/review1`): a title of a
// space at OID 1 prints `gnl|BL_ORD_ID|1` for sseqid and `BL_ORD_ID:1` for sacc; the title
// `BL_ORD_ID:1 own` at OID 1 prints its first word `BL_ORD_ID:1` for both.
#[test]
fn the_general_id_has_two_forms_and_a_word_of_that_form_stays_a_word() {
    let unnamed = &crate::report::defline::UNNAMED_SOURCE;
    let oid = crate::seqdb::Oid::of_unrestricted(1);
    let general = tabular_subject_id(b" ", unnamed, oid, false, "BLASTN").unwrap();
    assert_eq!(general, TabularSubjectId::General(oid));
    assert_eq!(&*general.id(), b"gnl|BL_ORD_ID|1");
    assert_eq!(&*general.accession(), b"BL_ORD_ID:1");
    assert_eq!(
        &*tabular_seqid(&general.accession(), oid, b" ", unnamed, false),
        b"gnl|BL_ORD_ID|1"
    );
    let word = tabular_subject_id(b"BL_ORD_ID:1 own", unnamed, oid, false, "BLASTN").unwrap();
    assert_eq!(word, TabularSubjectId::Word(b"BL_ORD_ID:1".to_vec()));
    assert_eq!(
        &*tabular_seqid(&word.accession(), oid, b"BL_ORD_ID:1 own", unnamed, false),
        b"BL_ORD_ID:1"
    );
    // Another sequence's accession of that form is a word of its title.
    let other = crate::seqdb::Oid::of_unrestricted(2);
    assert_eq!(
        &*tabular_seqid(b"BL_ORD_ID:1", other, b"BL_ORD_ID:1", unnamed, false),
        b"BL_ORD_ID:1"
    );
    assert!(tabular_subject_id(b", ,", unnamed, oid, true, "BLASTP")
        .unwrap_err()
        .to_string()
        .starts_with("database sequence 1 has a title"));
}

// A nucleotide's empty title whose taxid the taxonomy database names is made from the name
// (NCBI BLAST+ 2.17.0 prints `Homo` for taxid 9606 in the place of the id, line RG T5 and
// line DBID ordname); without a name it keeps the general id; a source LOSAT refuses is
// refused with the sequence's OID.
#[test]
fn an_empty_nucleotide_title_takes_the_first_word_of_its_source() {
    use crate::report::defline::SourceTaxname;
    let oid = crate::seqdb::Oid::of_unrestricted(3);
    let named = SourceTaxname::Named(std::sync::Arc::from(&b"Homo sapiens"[..]));
    let id = tabular_subject_id(b"", &named, oid, false, "BLASTN").unwrap();
    assert_eq!(id, TabularSubjectId::Word(b"Homo".to_vec()));
    assert_eq!(
        &*tabular_seqid(&id.accession(), oid, b"", &named, false),
        b"Homo"
    );
    // A protein's empty title is `unnamed protein product [name]` with or without a name.
    assert_eq!(
        tabular_subject_id(b"", &named, oid, true, "BLASTP").unwrap(),
        TabularSubjectId::Word(b"unnamed".to_vec())
    );
    assert_eq!(
        tabular_subject_id(b"", &SourceTaxname::Unnamed, oid, false, "BLASTN").unwrap(),
        TabularSubjectId::General(oid)
    );
    assert!(
        tabular_subject_id(b"", &SourceTaxname::Refused("x"), oid, false, "BLASTN")
            .unwrap_err()
            .to_string()
            .starts_with("database sequence 3 has x")
    );
}

// NCBI BLAST+ 2.17.0 (probe of the review of line RG item T5c,
// `remaining-gaps-20261112/T/t5c/review/out/ncbi/sall.blastp.f6.*`): the protein title
// `MAG: (plastid) [ ]` of a sequence whose taxid no taxonomy database names (a source without
// a name) is cut to `MAG:` by `x_AdjustProteinTitleSuffix` and keeps the general id in both
// forms (`gnl|BL_ORD_ID|7`, `BL_ORD_ID:7`); without a source it gives `(plastid)`.
#[test]
fn a_rewritten_protein_title_without_a_word_keeps_the_general_id() {
    use crate::report::defline::SourceTaxname;
    let oid = crate::seqdb::Oid::of_unrestricted(7);
    let title = b"MAG: (plastid) [ ]";
    let id = tabular_subject_id(title, &SourceTaxname::Nameless, oid, true, "BLASTP").unwrap();
    assert_eq!(id, TabularSubjectId::General(oid));
    assert_eq!(
        &*tabular_seqid(&id.accession(), oid, title, &SourceTaxname::Nameless, true),
        b"gnl|BL_ORD_ID|7"
    );
    assert_eq!(
        tabular_subject_id(title, &SourceTaxname::Unnamed, oid, true, "BLASTP").unwrap(),
        TabularSubjectId::Word(b"(plastid)".to_vec())
    );
}

use super::*;

// Probes p5_pr_sum2, p5b_sum_int and p5_neg_2p63 of line RG D1a-3 (NCBI 2.17.0 prints the
// wrapped sums and `-0`): the report counts LOSAT refuses, by format (not -2^63, whose digits
// it writes: approved exception 17, Int8 -2^63 written as digits); the prolog's (outfmt 0 to
// 4) at open (`check_report_sums`), the documents' when they are written (`report_sum_refusal`,
// after NCBI's own errors of the search: p5_pr_sum2#6).
#[test]
fn report_sums_past_ncbis_types_are_refused_for_the_formats_that_sum_them() {
    let info = |num_sequences: i32, total_letters: i64| DatabaseInfo {
        definition: b"t".to_vec(),
        date: String::new(),
        num_sequences,
        total_letters,
        filt_algorithm_name: String::new(),
        filt_algorithm_options: String::new(),
    };
    let two = |a: DatabaseInfo, b: DatabaseInfo| PairwiseDatabase::new(vec![a, b]);
    let refused = |database: &PairwiseDatabase, outfmt: i32| {
        report_sum_refusal(
            database.num_sequences,
            database.total_letters_overflow,
            outfmt,
        )
        .is_some()
    };
    let lengths = two(
        info(23, 5_000_000_000_000_000_000),
        info(23, 5_000_000_000_000_000_000),
    );
    assert!(lengths.total_letters_overflow);
    for outfmt in [0, 1, 4] {
        assert!(
            check_report_sums(&lengths, outfmt, "BLASTP").is_err(),
            "{outfmt}"
        );
    }
    for outfmt in [5, 13, 14, 15, 16, 6, 7, 8, 10, 11, 17] {
        assert!(
            check_report_sums(&lengths, outfmt, "BLASTP").is_ok(),
            "{outfmt}"
        );
    }
    for outfmt in [0, 1, 4, 5, 13, 14, 15, 16] {
        assert!(refused(&lengths, outfmt), "{outfmt}");
    }
    for outfmt in [6, 7, 8, 10, 11, 17] {
        assert!(!refused(&lengths, outfmt), "{outfmt}");
    }
    let counts = two(info(2_000_000_000, 3263), info(2_000_000_000, 3263));
    assert!(refused(&counts, 5));
    assert!(check_report_sums(&counts, 5, "BLASTP").is_ok());
    for outfmt in [0, 15, 16] {
        assert!(!refused(&counts, outfmt), "{outfmt}");
    }
    // Two lengths of 2^64 - 1 (-1 as NCBI's Int8) sum to -2 without an overflow.
    assert!(check_report_sums(&two(info(23, -1), info(23, -1)), 0, "BLASTP").is_ok());
    // A length of exactly -2^63 is written as its digits (PD-LOSAT-NCBI-DEFECTS approved
    // exception 17, Int8 -2^63 written as digits: NCBI's XML writer prints -0).
    let int64_min = PairwiseDatabase::new(vec![info(23, i64::MIN)]);
    for outfmt in [0, 5, 6, 13, 14, 15, 16] {
        assert!(!refused(&int64_min, outfmt), "{outfmt}");
    }
}

// The report-sum check of the prolog of outfmt 0 to 4 at open goes through each format's
// refusal (`FormatRefusal::check`, line RG T7-W): the command line's format fails the call, a
// format of LOSAT Web's that fails alone is refused and the others stay; the database's
// report is made once, only where a format prints its sums.
#[test]
fn the_report_sums_of_each_format_go_to_its_refusal() {
    use crate::api::local_blast::FormatRefusal;
    let info = |total_letters: i64| DatabaseInfo {
        definition: b"t".to_vec(),
        date: String::new(),
        num_sequences: 23,
        total_letters,
        filt_algorithm_name: String::new(),
        filt_algorithm_options: String::new(),
    };
    let made = std::cell::Cell::new(0);
    let overflow = || {
        made.set(made.get() + 1);
        PairwiseDatabase::new(vec![
            info(5_000_000_000_000_000_000),
            info(5_000_000_000_000_000_000),
        ])
    };
    let mut call = FormatRefusal::FailsCall;
    let error = check_formats_report_sums(overflow, [("0", &mut call)], "BLASTN").unwrap_err();
    assert!(error.is::<ReportSumRefusal>());
    let (mut pairwise, mut tabular, mut query_anchored) = (
        FormatRefusal::FailsAlone,
        FormatRefusal::FailsCall,
        FormatRefusal::FailsAlone,
    );
    check_formats_report_sums(
        overflow,
        [
            ("0", &mut pairwise),
            ("6", &mut tabular),
            ("3", &mut query_anchored),
        ],
        "BLASTN",
    )
    .unwrap();
    assert!(pairwise.refused() && query_anchored.refused());
    assert!(matches!(tabular, FormatRefusal::FailsCall));
    assert!(pairwise
        .take_refused()
        .is_some_and(|error| error.is::<ReportSumRefusal>()));
    assert_eq!(made.get(), 2);
    let mut alone = FormatRefusal::FailsAlone;
    check_formats_report_sums(overflow, [("6", &mut alone), ("5", &mut call)], "BLASTN").unwrap();
    assert!(!alone.refused());
    assert_eq!(
        made.get(),
        2,
        "no format printing the prolog's sums: not made"
    );
}
use crate::utils::matrix::BLASTAA_SIZE;

// `check_shown_protein_title` asks for no source name: the title that NCBI makes from the
// source of a protein without parsed Seq-ids (an empty stored title) is `unnamed protein
// product` and ` [name]`, which no name of punctuation makes read past its end, while a
// RefSeq record's (parsed Seq-ids), the name alone, can: that one is checked from its own
// source (`check_shown_parsed_subject`), and the block search's check never takes it.
#[test]
fn a_stored_protein_title_from_any_source_ends_in_its_own_text() {
    use crate::report::defline::{generate_bioseq_defline, BioseqFacts};
    let reads_past_end = |title: &[u8], facts: &BioseqFacts| {
        [false, true].into_iter().any(|leave| {
            generate_bioseq_defline(title, true, leave, facts)
                .is_ok_and(|defline| defline.reads_past_end)
        })
    };
    for name in [
        &b", ,"[..],
        b"; ;",
        b",, ,",
        b"~, ,",
        b", , ,",
        b" , ",
        b"&#44; ,",
        b"x, ,",
    ] {
        let stored = BioseqFacts::of_source(Some(std::sync::Arc::from(name)));
        assert!(!reads_past_end(b"", &stored), "{name:?}");
    }
    assert!(
        check_shown_protein_title(b"", || crate::seqdb::Oid::of_unrestricted(0), "BLASTP").is_ok()
    );
    for name in [&b", ,"[..], b"; ;", b"~, ,", b", , ,", b"&#44; ,"] {
        let refseq = BioseqFacts {
            refseq_chromosome_or_ncrna: true,
            ..BioseqFacts::of_source(Some(std::sync::Arc::from(name)))
        };
        assert!(reads_past_end(b"", &refseq), "{name:?}");
    }
    // A stored title of punctuation is not made from the source: it is refused whatever the
    // source is.
    assert!(
        check_shown_protein_title(b", ,", || crate::seqdb::Oid::of_unrestricted(7), "BLASTP")
            .is_err()
    );
    assert!(!crate::report::defline::bioseq_title_from_source(
        b", ,",
        BioseqFacts::default
    ));
}

// Line DBID pfbatch: each batch's prefetch counts its own Seq-aligns against a positive
// `PRE_FETCH_SEQS_LIMIT` (`s_PreFetchSeqs`, cited at `prefetch_may_be_off`).
#[test]
fn each_batch_prefetches_by_its_own_seq_aligns() {
    // Seq-aligns of the queries 0 to 4 of a call.
    let counts = [2, 1, 0, 3, 1];
    use crate::report::archive::options::RequestKind;
    let table = PrefetchFormat::of(6, RequestKind::Blastp);
    let of_format = |limit, format, first, batches: &[std::ops::Range<usize>]| {
        let seq_aligns = |span: std::ops::Range<usize>| counts[span].iter().sum();
        let prefetch = BatchPrefetch::at(
            limit,
            format,
            counts.len(),
            first,
            batches,
            (seq_aligns, |_| false),
        );
        // A query of no batch (`None`, taken both ways) shows as `true`.
        (0..counts.len() + 1)
            .map(|query| prefetch.batch_of(query).unwrap_or(true))
            .collect::<Vec<_>>()
    };
    let of =
        |limit, first, batches: &[std::ops::Range<usize>]| of_format(limit, table, first, batches);
    let every = vec![true; counts.len() + 1];
    // No limit and `INT_MAX`: every batch prefetches; 0 or below: none does.
    for limit in [None, Some(i32::MAX)] {
        assert_eq!(of(limit, 0, &[0..2, 2..5]), every, "{limit:?}");
    }
    for limit in [Some(0), Some(-1)] {
        assert_eq!(of(limit, 0, &[0..2, 2..5]), [false; 6], "{limit:?}");
    }
    // Batches of 3 and 4 Seq-aligns: limit 3 turns the second off; a query beyond the call
    // is of no batch.
    assert_eq!(
        of(Some(3), 0, &[0..2, 2..5]),
        [true, true, false, false, false, true]
    );
    assert_eq!(of(Some(4), 0, &[0..2, 2..5]), every);
    // No batches: the call is one batch of 7 Seq-aligns.
    assert_eq!(
        of(Some(6), 0, &[]),
        [false, false, false, false, false, true]
    );
    assert_eq!(of(Some(7), 0, &[]), every);
    // Ranges numbered from the call's first query 10; a batch before it is not this call's.
    assert_eq!(
        of(Some(1), 10, &[3..10, 10..11, 11..13, 13..15]),
        [false, true, true, false, false, true]
    );
    // Review 5 of line DBID pfoff-4, M2 and H2: the Seq-annots and psiblast's archive prefetch
    // only at `INT_MAX`, whatever the count; the other applications' archive never.
    for outfmt in [8, 9, 12] {
        for app in [RequestKind::Blastp, RequestKind::Psiblast] {
            let format = PrefetchFormat::of(outfmt, app);
            for limit in [None, Some(4), Some(7), Some(0)] {
                assert_eq!(
                    of_format(limit, format, 0, &[0..2, 2..5]),
                    [false; 6],
                    "{outfmt} {app:?} {limit:?}"
                );
            }
            assert_eq!(of_format(Some(i32::MAX), format, 0, &[0..2, 2..5]), every);
        }
    }
    let archive = |app| PrefetchFormat::of(11, app);
    assert_eq!(
        of_format(None, archive(RequestKind::Psiblast), 0, &[0..2, 2..5]),
        [false; 6]
    );
    assert_eq!(
        of_format(Some(i32::MAX), archive(RequestKind::Psiblast), 0, &[0..2]),
        every
    );
    for app in [
        RequestKind::Blastp,
        RequestKind::Rpsblast,
        RequestKind::Tblastn,
    ] {
        assert_eq!(
            of_format(Some(i32::MAX), archive(app), 0, &[0..2]),
            [false; 6]
        );
    }
    // The XML2 and JSON formats prefetch as outfmt 0 does.
    for outfmt in [0, 5, 7, 13, 16, 17, 18] {
        let format = PrefetchFormat::of(outfmt, RequestKind::Psiblast);
        assert_eq!(
            of_format(Some(3), format, 0, &[0..2, 2..5]),
            of(Some(3), 0, &[0..2, 2..5])
        );
        assert_eq!(of_format(None, format, 0, &[]), every);
    }
    // Line DBID weakbatch: what the batches read of each subject they report, by ordinal id:
    // 7 in both batches, 8 in the first (3 Seq-aligns, on under limit 3), 9 in the second
    // (4, off), whose reports read it one by one (the walk's, `report::db_subject::ScopeWalk`).
    let sum = |span: std::ops::Range<usize>| counts[span].iter().sum();
    let never = |_: std::ops::Range<usize>| false;
    let reads = BatchPrefetch::at(Some(3), table, 5, 0, &[0..2, 2..5], (sum, never)).reads([
        (0, 7),
        (4, 7),
        (1, 8),
        (3, 9),
    ]);
    assert!(reads.prefetched(7));
    assert!(reads.prefetched(8));
    assert!(!reads.prefetched(9));
    let every = BatchPrefetch::at(None, table, 5, 0, &[0..2, 2..5], (sum, never)).reads([(3, 9)]);
    assert!(every.prefetched(9));
    let none = BatchPrefetch::at(Some(0), table, 5, 0, &[0..2, 2..5], (sum, never)).reads([(0, 7)]);
    assert!(!none.prefetched(7));
    // Review 5 of line DBID pfoff-4, H3: psiblast's outfmt 8 reads its subjects one by one.
    let seq_annot = PrefetchFormat::of(8, RequestKind::Psiblast);
    let annot =
        BatchPrefetch::at(None, seq_annot, 5, 0, &[0..2, 2..5], (sum, never)).reads([(3, 9)]);
    assert!(!annot.prefetched(9));
    // Review 1 M1: a query of no batch (5, beyond the call's 5) is taken both ways.
    let beyond =
        BatchPrefetch::at(Some(3), table, 5, 0, &[0..2, 2..5], (sum, never)).reads([(5, 6)]);
    assert!(beyond.prefetched(6));
    // Line DBID pfoff-7: a batch whose prefetch fails at a subject whose lookup throws on the id
    // (the second, `2..5`) prefetches, and its reports load its subjects one by one after it; the
    // reads whose texts its reports write stay the prefetch's. A batch whose prefetch is off is
    // not asked.
    let second = |span: std::ops::Range<usize>| span.start == 2;
    let failing = BatchPrefetch::at(None, table, 5, 0, &[0..2, 2..5], (sum, second));
    assert!(failing.prefetches(3) && failing.loads_one_by_one(3));
    assert!(!failing.loads_one_by_one(1));
    assert_eq!(
        (
            failing.report_read(3, 0).one_by_one,
            failing.report_read(3, 0).prefetch_fails
        ),
        (true, true)
    );
    assert_eq!(
        (
            failing.report_read(0, 0).one_by_one,
            failing.report_read(0, 0).prefetch_fails
        ),
        (false, false)
    );
    let fails_reads = failing.reads([(0, 7), (4, 7), (1, 8), (3, 9)]);
    assert!(fails_reads.prefetched(9));
    assert!(fails_reads.prefetched(8));
    assert!(fails_reads.prefetched(7));
    let off = BatchPrefetch::at(
        Some(0),
        table,
        5,
        0,
        &[0..2, 2..5],
        (sum, |_| unreachable!()),
    );
    assert!(off.loads_one_by_one(3) && !off.report_read(3, 0).prefetch_fails);
    let shown = ParsedShown {
        prefetched: true,
        ..ParsedShown::default()
    };
    let of_9 = reads.shown(shown, 9);
    assert!(!of_9.prefetched);
    assert!(prefetch_off_at(Some(1), 2) && !prefetch_off_at(Some(2), 2));
    assert!(prefetch_off_at(Some(0), 0) && !prefetch_off_at(Some(i32::MAX), usize::MAX));
    assert!(!prefetch_off_at(None, usize::MAX));
    // The format's part of the rule (`PrefetchFormat::reads_at`) and its one-by-one question.
    assert!(table.reads_at(None) && table.reads_at(Some(0)));
    assert!(!seq_annot.reads_at(None) && !seq_annot.reads_at(Some(1)));
    assert!(seq_annot.reads_at(Some(i32::MAX)));
    assert!(!table.runs_at(Some(0), 0) && !seq_annot.runs_at(None, 0));
    assert!(table.runs_at(None, usize::MAX) && seq_annot.runs_at(Some(i32::MAX), usize::MAX));
}

// Review 3 of line DBID loadord, H2 (kit cases `P.off.*`): NCBI's data loader holds what a
// batch loaded, not what it reported. A batch whose prefetch is on loads every subject of its
// Seq-aligns; one whose prefetch `PRE_FETCH_SEQS_LIMIT` turns off loads what its reports' rows
// read (`6 qseqid evalue`: nothing).
#[test]
fn a_batch_loads_its_seq_aligns_subjects_or_what_its_reports_read() {
    use crate::report::archive::options::RequestKind;
    use crate::report::tabular::SubjectFetch::{Every, Never, WithIdentities};
    use BatchLoads::{Prefetched, Read};
    let table = PrefetchFormat::of(6, RequestKind::Blastp);
    for rows in [Every, Never, WithIdentities] {
        // No limit, `INT_MAX`, a batch of no more Seq-aligns than the limit: the prefetch.
        assert_eq!(batch_loads_at(None, table, 4, rows), Prefetched);
        assert_eq!(batch_loads_at(Some(i32::MAX), table, 4, rows), Prefetched);
        assert_eq!(batch_loads_at(Some(4), table, 4, rows), Prefetched);
        // Above the limit, or a limit of 0: the reports' rows. Review 6, H1: a prefetch-off
        // batch whose rows read every subject is not a prefetched one (psiblast's prefetch
        // reads its traceback beyond `-max_target_seqs`, its reports the cut lists).
        assert_eq!(batch_loads_at(Some(2), table, 4, rows), Read(rows));
        assert_eq!(batch_loads_at(Some(0), table, 0, rows), Read(rows));
        assert_ne!(batch_loads_at(Some(2), table, 4, rows), Prefetched);
        // Review 5 of line DBID pfoff-4, H2: psiblast's outfmt 8 and 11 load only what their
        // reports read (nothing) and the PSSM, but at `INT_MAX`.
        for outfmt in [8, 11] {
            let format = PrefetchFormat::of(outfmt, RequestKind::Psiblast);
            assert_eq!(batch_loads_at(None, format, 4, rows), Read(rows));
            assert_eq!(batch_loads_at(Some(4), format, 4, rows), Read(rows));
            assert_eq!(batch_loads_at(Some(i32::MAX), format, 4, rows), Prefetched);
        }
        // Review 7, L2: where the prefetched set is the reported set, the loads read every
        // reported subject where the prefetch ran, else what the rows read.
        assert_eq!(Prefetched.reported_rows(), Every);
        assert_eq!(Read(rows).reported_rows(), rows);
    }
    // Review 4, L4: the one walk of the held batches (`LoadedBefore::subjects`), in their
    // order, each batch's loads from its own Seq-aligns; a batch past the held hits is not
    // held.
    let hits = [2usize, 1, 0, 3];
    let counted = std::cell::RefCell::new(Vec::new());
    let walk = |limit: Option<i32>, cuts: XmlCuts<'_>| -> Vec<_> {
        counted.borrow_mut().clear();
        LoadedBefore {
            hits: &hits,
            spans: &[0..2, 2..4, 4..6],
            format: table,
            rows: Never,
            cuts,
        }
        .subjects_at(
            limit,
            |hits, span| hits.get(span),
            |span, hits| {
                counted.borrow_mut().push(span);
                hits.iter().sum()
            },
            |span, hits, loads| [(span, hits.to_vec(), loads)],
        )
        .collect()
    };
    assert_eq!(
        walk(None, XmlCuts::NONE),
        [
            (0..2, vec![2, 1], Prefetched),
            (2..4, vec![0, 3], Prefetched)
        ]
    );
    assert_eq!(*counted.borrow(), [0..2, 2..4]);
    // The prefetch off (limit 2: batch 0..2 has 3 Seq-aligns, 2..4 has 3): each query's
    // report reads, by query.
    assert_eq!(
        walk(Some(2), XmlCuts::NONE),
        [
            (0..1, vec![2], Read(Never)),
            (1..2, vec![1], Read(Never)),
            (2..3, vec![0], Read(Never)),
            (3..4, vec![3], Read(Never)),
        ]
    );
    // Review 6, H2: an outfmt 5 report that NCBI's XML stream cut before its hits (query 1)
    // makes none of its loads; the cut does not touch a batch whose prefetch ran (limit 3:
    // batch 0..2 prefetches, batch 2..4 too).
    let cut = |query: usize| query == 1 || query == 2;
    assert_eq!(
        walk(Some(2), XmlCuts::of(true, &cut)),
        [(0..1, vec![2], Read(Never)), (3..4, vec![3], Read(Never))]
    );
    assert_eq!(
        walk(Some(3), XmlCuts::of(true, &cut)),
        [
            (0..2, vec![2, 1], Prefetched),
            (2..4, vec![0, 3], Prefetched)
        ]
    );
    // Without outfmt 5 first nothing is cut.
    assert_eq!(walk(Some(2), XmlCuts::of(false, &cut)).len(), 4);
    // Review 7, L2: the one walk of a batch's reports (`BatchLoads::report_spans`, which
    // BLASTX's record walks too), also of a batch that does not start at query 0 (BLASTX's
    // search batches, by query of the reader batch).
    let spans = |loads: BatchLoads, cuts: XmlCuts<'_>, span| -> Vec<std::ops::Range<usize>> {
        loads.report_spans(cuts, span).collect()
    };
    assert_eq!(spans(Prefetched, XmlCuts::of(true, &cut), 1..3), [1..3]);
    assert_eq!(spans(Read(Every), XmlCuts::of(true, &cut), 1..4), [3..4]);
    assert_eq!(spans(Read(Every), XmlCuts::NONE, 1..3), [1..2, 2..3]);
    assert_eq!(spans(Read(Never), XmlCuts::NONE, 2..2), []);
}

/// The subjects of the database `name` of `tests/db_search/ids` by ordinal id.
fn parsed_subjects(name: &str) -> Vec<crate::report::db_subject::ParsedSubject> {
    let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/db_search/ids")
        .join(name);
    let db = SeqDb::open(base.to_str().unwrap(), SeqType::Protein).unwrap();
    (0..db.num_seqs())
        .map(|ordinal| {
            crate::report::db_subject::ParsedSubject::new(db.defline_set(ordinal).unwrap()).unwrap()
        })
        .collect()
}

// Line DBID N3 review 1 M1: the joined copies' letters are the same because `copy_of` reads
// them: the definition lines of r1 and r3 with other letters (pjoin_x: four substitutions, ten
// more letters) are no copies, those with other taxids (pjoin_abt) are.
#[test]
fn a_copy_has_the_same_letters_and_definition_lines_but_for_taxids() {
    let open = |name: &str| {
        let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/db_search/ids")
            .join(name);
        SeqDb::open(base.to_str().unwrap(), SeqType::Protein).unwrap()
    };
    let copy = |db: &SeqDb, a: u32, b: u32| {
        let subject = match db.read_header(a).unwrap() {
            crate::seqdb::SequenceHeader::Parsed(set) => {
                crate::report::db_subject::ParsedSubject::new(set).unwrap()
            }
            _ => panic!("parsed Seq-ids"),
        };
        HeaderReader::new(db, false).copy_of(&subject, a, b)
    };
    // pjoin_ax: a (r1 0, r2 1, r3 2), x (r1 3, r3 4).
    let ax = open("pjoin_ax");
    assert_eq!(copy(&ax, 0, 3), RecordCopy::Differs);
    assert_eq!(copy(&ax, 2, 4), RecordCopy::Differs);
    // pjoin_abc / pjoin_abt: a (r1 0, r2 1, r3 2), b (r1 3, r3 4, r4 5, r6 6), c (r1 7).
    let abc = open("pjoin_abc");
    assert_eq!(copy(&abc, 0, 3), RecordCopy::Same);
    assert_eq!(copy(&abc, 2, 4), RecordCopy::Same);
    assert_eq!(copy(&abc, 0, 6), RecordCopy::Differs);
    let abt = open("pjoin_abt");
    assert_eq!(copy(&abt, 0, 3), RecordCopy::OtherTaxids);
    assert_eq!(copy(&abt, 0, 7), RecordCopy::Same);
}

// Line RG D1a-4 (fix round 7, review F N3): the one count of NCBI's cut (`PruneSeqalign`),
// which `joined_cut` and RPS-BLAST's `formatter_kept` take: following entries of one subject
// count once, the entry past the cut ends the count unread by `joined`, and an error of an
// entry or of `joined` ends it.
#[test]
fn the_cut_counts_following_entries_of_one_subject_once() {
    let ids = ["a", "a", "b", "c", "c", "c", "d"];
    let count =
        |cut: usize| pruned_kept(cut, ids.iter().map(Ok::<_, ()>), |a, b| Ok(a == b)).unwrap();
    let kept: Vec<usize> = (0..6).map(count).collect();
    assert_eq!(kept, [0, 2, 3, 6, 7, 7]);
    let mut asked = Vec::new();
    assert_eq!(
        pruned_kept(2, ids.iter().map(Ok::<_, ()>), |a, b| {
            asked.push((**a, **b));
            Ok(a == b)
        }),
        Ok(3)
    );
    assert_eq!(asked, [("a", "a"), ("a", "b"), ("b", "c")]);
    assert_eq!(
        pruned_kept(9, ids.iter().map(Ok), |_, b| if **b == "c" {
            Err(4)
        } else {
            Ok(false)
        }),
        Err(4)
    );
    let entries = [Ok("a"), Err(5), Ok("b")];
    assert_eq!(pruned_kept(9, entries, |a, b| Ok(a == b)), Err(5));
}

// Line DBID N3: in NCBI's order of pjoin_abc's hits (r1 of c, r6, r1 of b and a, r2, r4,
// r3 of b and a; `make_join_fasta.py`) the cut to N subjects keeps the copies that follow a
// kept subject (psiblast, `PruneSeqalign`), and the hits show them in one hit; pfound_v5's
// pair 4 (gb|DDD00004| then gb|DDD00004.1|) is refused where the count reaches it.
#[test]
fn copies_count_as_one_subject_in_the_cut_and_other_matching_ids_are_refused() {
    let s = parsed_subjects("pjoin_abc");
    let order = [7u32, 6, 3, 0, 1, 5, 4, 2];
    let hits = || {
        order
            .iter()
            .map(|&ordinal| Some((ordinal, &s[ordinal as usize])))
    };
    for (max, kept) in [(1, 1), (2, 2), (3, 4), (4, 5), (5, 6), (6, 8), (7, 8)] {
        assert_eq!(
            joined_cut(
                0,
                max,
                hits(),
                crate::seqdb::OidMap::unrestricted(),
                "PSIBLAST"
            )
            .unwrap(),
            kept,
            "{max}"
        );
    }
    assert_eq!(
        joined_cut(
            0,
            0,
            hits(),
            crate::seqdb::OidMap::unrestricted(),
            "PSIBLAST"
        )
        .unwrap(),
        0
    );
    refuse_joined_hits(0, hits(), crate::seqdb::OidMap::unrestricted(), "BLASTP").unwrap();
    // A subject without parsed Seq-ids between copies parts them.
    let parted = [Some((3, &s[3])), None, Some((0, &s[0]))];
    assert_eq!(
        joined_cut(
            0,
            1,
            parted,
            crate::seqdb::OidMap::unrestricted(),
            "PSIBLAST"
        )
        .unwrap(),
        1
    );
    let f = parsed_subjects("pfound_v5");
    let pair = f
        .iter()
        .enumerate()
        .filter(|(_, subject)| subject.alignment_fasta().starts_with(b"gb|DDD00004"))
        .map(|(ordinal, subject)| Some((ordinal as u32, subject)))
        .collect::<Vec<_>>();
    assert_eq!(pair.len(), 2);
    let error = refuse_joined_hits(
        1,
        pair.clone(),
        crate::seqdb::OidMap::unrestricted(),
        "BLASTP",
    )
    .unwrap_err()
    .to_string();
    assert!(
        error.contains("but are not copies of one record"),
        "{error}"
    );
    // The cut to one subject reaches the second of the pair: NCBI would keep both.
    let error = joined_cut(
        1,
        1,
        pair.clone(),
        crate::seqdb::OidMap::unrestricted(),
        "PSIBLAST",
    )
    .unwrap_err()
    .to_string();
    assert!(
        error.contains("the cut to -max_target_seqs 1 counts"),
        "{error}"
    );
    // A cut that ends before the pair does not reach it.
    let before = [Some((7, &s[7]))].into_iter().chain(pair);
    assert_eq!(
        joined_cut(
            1,
            1,
            before,
            crate::seqdb::OidMap::unrestricted(),
            "PSIBLAST"
        )
        .unwrap(),
        1
    );
}

/// The residues between the two sentinel bytes.
fn stored(residues: &[u8]) -> Vec<u8> {
    let mut bytes = vec![NCBISTDAA_GAP];
    bytes.extend_from_slice(residues);
    bytes.push(NCBISTDAA_GAP);
    bytes
}

fn check(residues: &[u8]) -> std::result::Result<(), String> {
    check_protein_sequence(
        &stored(residues),
        || crate::seqdb::Oid::of_unrestricted(7),
        "blastx",
    )
    .map_err(|e| e.to_string())
}

/// The check as an `all` over the residues.
fn accepted_by_all(residues: &[u8]) -> bool {
    residues
        .iter()
        .all(|&residue| usize::from(residue) < BLASTAA_SIZE)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1559-1593 (cited at
// `check_protein_sequence`): the bytes past the matrix are what the check refuses; the gap
// byte is a residue of the matrix. Not a port of an NCBI unit case.
#[test]
fn the_check_accepts_the_letters_and_names_the_first_bad_byte() {
    let first = NCBISTDAA_GAP + 1;
    let last = (BLASTAA_SIZE - 1) as u8;
    assert_eq!(check(&[]), Ok(()));
    assert_eq!(check(&[first]), Ok(()));
    assert_eq!(check(&[last]), Ok(()));
    assert_eq!(check(&[NCBISTDAA_GAP]), Ok(()));
    let letters: Vec<u8> = (0..1_000)
        .map(|i| [first, last, NCBISTDAA_GAP][i % 3])
        .collect();
    assert_eq!(check(&letters), Ok(()));
    for bad in [BLASTAA_SIZE as u8, 255] {
        for length in [1, 2, 31, 32, 33, 100, 1_000] {
            for at in [0, length / 2, length - 1] {
                let mut residues = letters[..length].to_vec();
                residues[at] = bad;
                let error = check(&residues).unwrap_err();
                assert!(
                    error.contains("no NCBIstdaa residue")
                        && error.contains(&format!("position {},", at + 1)),
                    "byte {bad} at {at} of {length}: {error}"
                );
            }
        }
    }
    // The sentinels.
    assert!(check_protein_sequence(
        &[1, first, 0],
        || crate::seqdb::Oid::of_unrestricted(7),
        "blastx"
    )
    .is_err());
    assert!(check_protein_sequence(
        &[0, first, 1],
        || crate::seqdb::Oid::of_unrestricted(7),
        "blastx"
    )
    .is_err());
}

// Random residues around the letters' bounds: the single pass accepts exactly what an `all`
// accepts.
#[test]
fn the_check_accepts_what_the_check_before_it_accepted() {
    let mut seed = 0x9e37_79b9_7f4a_7c15_u64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    for case in 0..2_000 {
        let length = (next() % 200) as usize;
        let mut residues: Vec<u8> = (0..length).map(|_| 1 + (next() % 27) as u8).collect();
        if case % 2 == 0 && length > 0 {
            let at = (next() as usize) % length;
            residues[at] = [0, 27, 28, 29, 254, 255][(next() % 6) as usize];
        }
        assert_eq!(
            check(&residues).is_ok(),
            accepted_by_all(&residues),
            "{residues:?}"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1146-1155
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:768-788
/// A load of a subject with a taxid where a file of the taxonomy database does not open
/// (`SeqDbError::TaxdbThrows`) is LOSAT's refusal, in words that say what the user can do, which
/// the batch's prefetch alone catches (`batch_prefetch`, in a run whose formats NCBI writes the
/// same: `db_subject::unchanged_by_load_throws`); any other refusal of the taxonomy database is
/// not caught there (line WQA D18b, 2026-11-21; review 1, L1: `slen` and the organism report are
/// written the same).
#[test]
fn a_taxonomy_file_that_throws_is_a_refusal_that_the_prefetch_alone_catches() {
    use crate::blastinput::app::Program;
    let path = std::path::PathBuf::from("/d/taxdb.bti");
    // The exception that `SeqDbError::memory_file_exception` makes on Unix for this file. It is
    // written out here because `memory_file_exception` makes none on Windows, where LOSAT does not
    // have NCBI's system text (`ncbi_strerror`).
    let exception = crate::ncbi_exception::NcbiException::new(
        crate::ncbi_exception::ThrowSite {
            file: "src/corelib/ncbifile.cpp",
            line: 5945,
            kind: "CFileException::eMemoryMap",
            function: Some("ncbi::CMemoryFileMap::x_Open"),
        },
        "CMemoryFile: Cannot memory map file '/d/taxdb.bti':Permission denied",
    );
    #[cfg(unix)]
    assert_eq!(
        Some(&exception),
        SeqDbError::Open {
            path: path.clone(),
            source: std::io::Error::from_raw_os_error(13),
        }
        .memory_file_exception()
        .as_ref()
    );
    let error = taxdb_error(
        SeqDbError::TaxdbThrows {
            path: path.clone(),
            reason: "Permission denied".to_owned(),
            exception: Box::new(exception),
        },
        "BLASTP",
    );
    assert!(error.is::<TaxonomyRefusal>());
    assert!(TaxonomyRefusal::of_load_that_throws(&error));
    assert_eq!(
        format!("{error:#}"),
        "/d/taxdb.bti cannot be opened: Permission denied; NCBI BLAST+ 2.17.0 then fails at \
         each taxonomy lookup and writes this output differently or stops, which LOSAT's BLASTP \
         does not support; make the file readable, or move it off the BLASTDB path"
    );
    let other = taxdb_error(
        SeqDbError::Unsupported {
            path,
            message: "Permission denied (os error 13)".to_owned(),
        },
        "BLASTP",
    );
    assert!(other.is::<TaxonomyRefusal>());
    assert!(!TaxonomyRefusal::of_load_that_throws(&other));
    // The formats that NCBI writes the same: a table of fields that load no subject (`slen` and
    // the coordinates included, not the letters or the subject's ids), the organism report, the
    // Seq-annots.
    for (outfmts, unchanged) in [
        (&["6 qseqid evalue bitscore"][..], true),
        (&["6 qseqid slen"][..], true),
        (&["6 qseqid sstart send length"][..], true),
        (&["18"][..], true),
        (&["8"][..], true),
        (&["6 qseqid sseqid"][..], false),
        (&["6 qseqid pident"][..], false),
        (&["6"][..], false),
        (&["6 qseqid evalue bitscore", "0"][..], false),
        (&["5"][..], false),
        (&["15"][..], false),
    ] {
        assert_eq!(
            crate::report::db_subject::unchanged_by_load_throws_outfmts(
                outfmts.iter().copied(),
                Program::Blastp,
                crate::report::archive::options::RequestKind::Blastp,
            ),
            unchanged,
            "{outfmts:?}"
        );
    }
}

// RG-KB 9(b), 9b-2: only NCBI's Seq-aligns read a reported subject's definition line where every
// format writes them as they are (outfmt 8, 9, 11, 12) and the first does not prefetch the
// subjects (8, 9 and 12 prefetch with a PRE_FETCH_SEQS_LIMIT of INT_MAX, which the test process
// may set), for the five search programs (not psiblast, rpsblast or rpstblastn); a run of one
// format whose formatter reads it where LOSAT writes NCBI's posts (the tables, SAM, the
// Seq-annots) and whose batches all prefetch is `TitleReads::Prefetched`; the others read it
// again (`TitleReads::Reports`, LOSAT's refusal).
#[test]
fn only_the_seq_align_formats_read_titles_at_the_seq_aligns_alone() {
    use crate::report::archive::options::RequestKind;
    // The run's limit, which the caller reads once (review 3 of line RG-KB item 9b-2, R3-L2):
    // unset, 0, a positive count and `INT_MAX`.
    let limits = [None, Some(0), Some(3), Some(i32::MAX)];
    for (app, limit) in [
        RequestKind::Blastn,
        RequestKind::Blastp,
        RequestKind::Blastx,
        RequestKind::Tblastn,
        RequestKind::Tblastx,
    ]
    .into_iter()
    .flat_map(|app| limits.map(|limit| (app, limit)))
    {
        let of = |outfmts: &[i32]| TitleReads::of_outfmts(app, outfmts.iter().copied(), limit);
        let prefetched = |as_outfmt: i32| {
            let format = PrefetchFormat::of(as_outfmt, app);
            prefetches_subjects_at(Some(format), limit) && !prefetch_may_be_off(format, limit)
        };
        assert_eq!(of(&[11]), TitleReads::SeqAligns, "{app:?}");
        assert_eq!(of(&[11, 11]), TitleReads::SeqAligns, "{app:?}");
        for outfmts in [
            &[][..],
            &[0],
            &[3],
            &[5],
            &[13],
            &[18],
            &[11, 6],
            &[11, 0],
            &[6, 6],
            &[6, 11],
        ] {
            assert_eq!(of(outfmts), TitleReads::Reports, "{app:?} {outfmts:?}");
        }
        // The tables 10 and 20 prefetch as 6 does.
        for (outfmt, as_outfmt) in [(6, 6), (7, 7), (10, 6), (20, 6), (17, 17)] {
            let expected = if prefetched(as_outfmt) {
                TitleReads::Prefetched
            } else {
                TitleReads::Reports
            };
            assert_eq!(of(&[outfmt]), expected, "{app:?} {outfmt}");
        }
        // Only `INT_MAX` prefetches them, and every batch then.
        assert_eq!(of(&[8]) == TitleReads::Prefetched, limit == Some(i32::MAX));
        // A positive count or 0 may turn a batch's prefetch off.
        assert_eq!(
            of(&[6]) == TitleReads::Prefetched,
            matches!(limit, None | Some(i32::MAX)),
            "{app:?} {limit:?}"
        );
        for outfmt in [8, 9, 12] {
            let prefetch = prefetches_subjects_at(Some(PrefetchFormat::of(outfmt, app)), limit);
            let (alone, with_archive) = match (prefetch, prefetched(outfmt)) {
                (false, _) => (TitleReads::SeqAligns, TitleReads::SeqAligns),
                (true, true) => (TitleReads::Prefetched, TitleReads::Reports),
                (true, false) => (TitleReads::Reports, TitleReads::Reports),
            };
            assert_eq!(of(&[outfmt]), alone, "{app:?} {outfmt}");
            assert_eq!(of(&[outfmt, 11]), with_archive, "{app:?} {outfmt}");
        }
    }
    // psiblast, rpsblast and rpstblastn always stay `TitleReads::Reports`.
    for app in [
        RequestKind::Psiblast,
        RequestKind::Rpsblast,
        RequestKind::Rpstblastn,
    ] {
        for outfmts in [&[][..], &[8], &[9], &[11], &[12], &[8, 11], &[0], &[6]] {
            assert_eq!(
                TitleReads::of_outfmts(app, outfmts.iter().copied(), None),
                TitleReads::Reports,
                "{app:?} {outfmts:?}"
            );
        }
    }
    assert_eq!(TitleReads::default(), TitleReads::Reports);
}

// Line RG D1c-1c (review 2 G2 of D1c-1b): under a negative seqid list, a list of a redo's stream
// whose subject's header cannot be read (lr's sequence 1, `gnl|BL_ORD_ID|xyz` in its second
// volume, a string tag that NCBI's read of a later volume cannot make an OID of) stays in the
// stream, marked as kept and noted as unreadable, so that the redo tests its early termination
// as NCBI's does; the ordinary traceback, which fetches every subject, still fails at it, and the
// refusal of a fetch that reaches it gives the read's error.
#[test]
fn an_unreadable_header_stays_in_the_redo_stream_and_ends_the_traceback() {
    let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/idlist/lr");
    let mut db = SeqDb::open(base.to_str().unwrap(), SeqType::Protein).unwrap();
    let mut warnings = Vec::new();
    let list = crate::seqdb::SeqidList {
        negative: true,
        file: crate::seqdb::seqidlist::read_si_list(b"gi|102\n", &mut warnings).unwrap(),
    };
    db.restrict(&crate::seqdb::DbList::Seqids(list), &mut Vec::new())
        .unwrap();
    assert!(db.excludes_lineless_at_fetch());
    let fetch = FetchExcluded::of(&db).expect("a negative seqid list");
    // The stream's order: the strong subject, the weak one, the filler.
    let stream = [2u32, 1, 0];
    let (marks, reads) = fetch.marks_noting(&stream, |&ordinal| Some(ordinal), |&ordinal| ordinal);
    assert!((0..stream.len()).all(|index| !marks.excluded(index)));
    assert!(reads.utf8.is_empty());
    let unreadable: Vec<_> = reads
        .unreadable
        .iter()
        .map(|(read, noted)| (read.index, read.kept, read.before, read.ordinal, *noted))
        .collect();
    assert_eq!(unreadable, [(1, true, 1, 1, 1)]);
    // The ordinary traceback fetches every subject: the read fails, the items stay as they were.
    let mut items = stream.to_vec();
    let (ordinal, _) = fetch
        .retain(&mut items, |&ordinal| Some(ordinal))
        .unwrap_err();
    assert_eq!((ordinal, items.as_slice()), (1, stream.as_slice()));
    // A composition-based redo's stream keeps the list.
    let mut items = stream.to_vec();
    let reads = fetch.retain_noting(&mut items, |&ordinal| Some(ordinal), |_| ());
    assert_eq!(
        (items.as_slice(), reads.unreadable.len()),
        (stream.as_slice(), 1)
    );
    let refusal = FetchExcluded::refusal_at(&db, 1, "BLASTP").to_string();
    assert!(
        refusal.contains("those of sequence 1 cannot be read") && refusal.contains("string tag"),
        "{refusal}"
    );
}

// RG-KB 9(b) review 1 (F1): NCBI BLAST+ 2.17.0 (`db.tchar.*.dmg.*`) posts the 803.7 of the
// title of OID 0 (0xE9), which the batch's Seq-aligns read before OID 1, whose title tag it
// cannot read, and then fails; a later read of OID 0 finds it in the cache of header records.
#[test]
fn a_failed_header_read_keeps_the_posts_of_the_seq_aligns_before_it() {
    let base =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/tchar/tc_dmg_p");
    let db = SeqDb::open(base.to_str().unwrap(), SeqType::Protein).unwrap();
    use crate::blastinput::app::ErrorAfterBatchReports;
    let failure = |title_reads: Vec<u32>| {
        anyhow::Error::new(ErrorAfterBatchReports {
            query: Some(0),
            error: anyhow::anyhow!("unreadable"),
            utf8: None,
            title_reads,
        })
    };
    // The reads of the batch's Seq-aligns before the first read of the failed one, then it.
    let error = ErrorAfterBatchReports::after_seq_aligns(failure(vec![1]), [2, 0, 1, 2, 1]);
    let failed = error.downcast_ref::<ErrorAfterBatchReports>().unwrap();
    assert_eq!(failed.title_reads, [2, 0, 1]);
    // A failure elsewhere (a fetch, the prefetch) has no such reads.
    let elsewhere = ErrorAfterBatchReports::after_seq_aligns(failure(Vec::new()), [0, 1]);
    let elsewhere = elsewhere.downcast_ref::<ErrorAfterBatchReports>().unwrap();
    assert!(elsewhere.title_reads.is_empty());
    assert_eq!(
        elsewhere.title_posts(&db, false, "BLASTP").unwrap(),
        Default::default()
    );
    assert!(!db.replaced_titles());
    // With NCBI's warning of a UTF8String header in the batch: LOSAT's refusal.
    let utf8 = failed
        .title_posts(&db, true, "BLASTP")
        .unwrap_err()
        .to_string();
    assert!(
        utf8.contains("database sequence 0 has a stored title"),
        "{utf8}"
    );
    db.restart_header_reads();
    let mut posted = crate::report::diag_log::Posts::default();
    posted.post(
        crate::report::diag_log::ERROR,
        b"Error: (803.7) [blastp] Blast-def-line-set.E.title\nBad char [0xE9] in string at byte \
          18\nd0 fam\xe9one\n"
            .to_vec(),
    );
    assert_eq!(failed.title_posts(&db, false, "BLASTP").unwrap(), posted);
    assert!(db.replaced_titles());
    assert_eq!(
        failed.title_posts(&db, false, "BLASTP").unwrap(),
        Default::default()
    );
}

// RG-KB 9(b) review 1 (F3): a call's first UTF8String read posts NCBI's warning in its own
// batch alone.
#[test]
fn the_utf8_read_of_a_call_is_of_its_batch_alone() {
    let call = CallUtf8::default();
    assert!(!call.of_batch(0));
    let read = Utf8Read {
        query: 1,
        path: "Blast-def-line-set.E.title".to_string(),
    };
    call.note(2, Some(&read));
    call.note(4, Some(&read));
    assert!(call.of_batch(2));
    assert!(!call.of_batch(0) && !call.of_batch(4));
}

#[test]
fn outfmt_taxdb_use_matches_ncbi_format_rules() {
    use super::{outfmt_prefetches_taxdb, outfmt_taxdb_use, OutfmtTaxdbUse};
    use crate::report::archive::options::RequestKind;
    let app = RequestKind::Blastp;
    // Formats that check taxdb.bti even when absent: XML2/JSON (13..=16) on `-db`, and plural
    // tabular taxonomy columns (6/7/10/20) on both `-db` and `-subject`.
    for fmt in ["13", "14", "15", "16"] {
        assert_eq!(
            outfmt_taxdb_use(fmt, app, true),
            OutfmtTaxdbUse::Always,
            "{fmt} db"
        );
        assert_eq!(
            outfmt_taxdb_use(fmt, app, false),
            OutfmtTaxdbUse::None,
            "{fmt} subj"
        );
    }
    for fmt in [
        "6 qseqid sseqid sscinames",
        "6 qseqid sseqid scomnames",
        "6 qseqid sseqid sblastnames",
        "6 qseqid sseqid sskingdoms",
        "7 qseqid sseqid sscinames",
        "10 qseqid sseqid sscinames",
        "20 qseqid sseqid sscinames",
    ] {
        assert_eq!(
            outfmt_taxdb_use(fmt, app, true),
            OutfmtTaxdbUse::Always,
            "{fmt} db"
        );
        assert_eq!(
            outfmt_taxdb_use(fmt, app, false),
            OutfmtTaxdbUse::Always,
            "{fmt} subj"
        );
    }
    // Formats that read taxonomy names only when taxdb.bti already exists on the search path:
    // singular tabular name fields, outfmt 18, pairwise/query-anchored 0..=4 on `-db`, and
    // tabular formats (6/7/10/20) that print a subject id column (`PrintedIds`, including `std`).
    for fmt in [
        "0",
        "1",
        "2",
        "3",
        "4",
        "6",
        "6 std",
        "7",
        "10",
        "18",
        "20",
        "6 qseqid sseqid ssciname",
        "6 qseqid sseqid scomname",
        "6 qseqid sseqid sblastname",
        "6 qseqid sseqid sskingdom",
        "7 qseqid sseqid ssciname",
        "10 qseqid sseqid ssciname",
        "20 qseqid sseqid ssciname",
        "6 qseqid sseqid staxid",
        "6 qseqid sseqid staxids",
        "6 qseqid sacc",
        "6 qseqid saccver",
        "6 qseqid sallseqid",
        "6 qseqid sallacc",
    ] {
        assert_eq!(
            outfmt_taxdb_use(fmt, app, true),
            OutfmtTaxdbUse::WhenPresent,
            "{fmt} db"
        );
        assert_eq!(
            outfmt_taxdb_use(fmt, app, false),
            OutfmtTaxdbUse::None,
            "{fmt} subj"
        );
    }
    // Formats that do not consult taxdb for formatting: 5, 6/7/10/20 without tax name or
    // printed subject id columns (including staxid/staxids/sgi/stitle alone), 8, 9, 11, 12, 17.
    for fmt in [
        "5",
        "6 qseqid pident",
        "6 qseqid staxid",
        "6 qseqid staxids",
        "6 qseqid sgi sallgi stitle",
        "7 qseqid pident",
        "10 qseqid pident",
        "20 qseqid pident",
        "8",
        "9",
        "11",
        "12",
        "17",
    ] {
        assert_eq!(
            outfmt_taxdb_use(fmt, app, true),
            OutfmtTaxdbUse::None,
            "{fmt} db"
        );
        assert_eq!(
            outfmt_taxdb_use(fmt, app, false),
            OutfmtTaxdbUse::None,
            "{fmt} subj"
        );
    }
    for fmt in ["0", "4", "5", "6", "7", "17"] {
        assert!(outfmt_prefetches_taxdb(fmt), "{fmt}");
    }
    for fmt in ["8", "9", "10", "11", "12", "13", "15", "18", "20"] {
        assert!(!outfmt_prefetches_taxdb(fmt), "{fmt}");
    }
}

#[test]
fn linux_online_cpu_parsers_match_sysfs_and_proc_stat() {
    assert_eq!(parse_linux_cpu_list("0\n"), Some(1));
    assert_eq!(parse_linux_cpu_list("0-31\n"), Some(32));
    assert_eq!(parse_linux_cpu_list("0,2-5\n"), Some(5));
    assert_eq!(parse_linux_cpu_list("0-3,8,10-11"), Some(7));
    assert_eq!(parse_linux_cpu_list(""), None);
    assert_eq!(parse_linux_cpu_list("   \n"), None);
    assert_eq!(parse_linux_cpu_list("5-2"), None);
    assert_eq!(parse_linux_cpu_list("0,"), None);
    assert_eq!(parse_linux_cpu_list("abc"), None);

    let proc_stat = "cpu  10 20 30 40\ncpu0 1 2 3 4\ncpu1 5 6 7 8\ncpu12 1 1 1 1\nintr 0\n";
    assert_eq!(parse_linux_proc_stat_cpus(proc_stat), Some(3));
    assert_eq!(
        parse_linux_proc_stat_cpus("cpu  10 20 30 40\nintr 0\n"),
        None
    );

    assert_eq!(ncbi_effective_threads(1), 1);
    #[cfg(target_os = "linux")]
    assert!(ncbi_online_cpu_count() >= 1);
}
