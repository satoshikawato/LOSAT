use super::*;

#[test]
fn sanitize_cleans_white_space_like_ncbi() {
    assert_eq!(sanitize(b"  a  b\t\tc \r"), b"a b c");
    assert_eq!(sanitize(b"x\x01\x02y"), b"x y");
    assert_eq!(sanitize(b"caf\xc3\xa9!"), b"caf !");
    assert_eq!(sanitize(b" \t "), b"");
    assert_eq!(
        sanitize(b"File is not accessible:  `a b'"),
        b"File is not accessible: `a b'"
    );
    assert_eq!(sanitize(b"\"quoted\""), b"\"quoted\"");
}

#[test]
fn lines_are_saved_when_complete() {
    let mut log = DiagLog::default();
    log.feed(b"Warning: [blastn] one\nWar");
    log.feed(b"ning: [blastn]  two\n\n");
    log.post(ERROR, b"BLAST query error: x\n");
    // A line keeps the severity of its header: a load's taxonomy database error (NCBI BLAST+
    // 2.17.0 psiblast -outfmt 11: `code 2`).
    log.feed(b"Error: [psiblast] Error: Tax database file has wrong magic number.\n");
    log.feed(b"Critical: (110.6) [blastn] x\nFASTA-Reader: y\n");
    assert_eq!(
        log.take(),
        vec![
            (WARNING, b"Warning: [blastn] one".to_vec()),
            (WARNING, b"Warning: [blastn] two".to_vec()),
            (ERROR, b"BLAST query error: x".to_vec()),
            (
                ERROR,
                b"Error: [psiblast] Error: Tax database file has wrong magic number.".to_vec()
            ),
            (CRITICAL, b"Critical: (110.6) [blastn] x".to_vec()),
            (WARNING, b"FASTA-Reader: y".to_vec()),
        ]
    );
    assert!(log.is_empty());
}

// The later formats' lines of a query (`Posts::later`) are not the diagnostics': they are
// written to a later format's receiver alone, before the query's report
// (`QueryWarnings::with_later`), and kept by the copies of the posts.
#[test]
fn later_lines_go_to_a_later_format_s_receiver_alone() {
    use crate::report::query_warnings::QueryWarnings;
    let mut first = Posts::of_lines(b"Warning: [blastn] reader\n".to_vec());
    first.later(b"Warning: [blastn] Query_1 q: Sequence contains no data \n");
    let mut second = Posts::new();
    second.later(b"");
    assert!(second.is_empty());
    second.later(b"Warning: [blastn] Query_2 r: Sequence contains no data \n");
    assert!(!second.has_text() && second.text().is_empty());
    let posts = Posts::with_later(vec![first, second, Posts::new()], &[Vec::new()]);
    let copies: Vec<Posts> = posts.iter().map(Posts::placed).collect();
    assert_eq!(copies, posts);
    // The first format: the diagnostics get the reader's line, no receiver anything.
    let (mut report, mut diagnostics) = (Vec::new(), Vec::new());
    let mut warnings = QueryWarnings::of_posts(&posts, None, &mut diagnostics);
    for query in 0..3 {
        warnings.before_query(query, &mut report).unwrap();
    }
    drop(warnings);
    assert_eq!(diagnostics, b"Warning: [blastn] reader\n");
    // A later format whose report stops after the first query: the diagnostics get nothing,
    // its receiver the first query's line.
    let (mut diagnostics, mut later) = (Vec::new(), Vec::new());
    let mut warnings =
        QueryWarnings::of_posts(&[], None, &mut diagnostics).with_later(&posts, Some(&mut later));
    warnings.before_query(0, &mut report).unwrap();
    drop(warnings);
    assert!(diagnostics.is_empty() && report.is_empty());
    assert_eq!(
        later,
        b"Warning: [blastn] Query_1 q: Sequence contains no data \n"
    );
    // One that writes every query's report: both lines, in the order of the queries.
    let mut later = Vec::new();
    let mut warnings =
        QueryWarnings::of_posts(&[], None, &mut diagnostics).with_later(&posts, Some(&mut later));
    for query in 0..3 {
        warnings.before_query(query, &mut report).unwrap();
    }
    drop(warnings);
    assert_eq!(
        later,
        b"Warning: [blastn] Query_1 q: Sequence contains no data \n\
              Warning: [blastn] Query_2 r: Sequence contains no data \n"
    );
}
