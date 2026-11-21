use super::*;
use crate::algorithm::common::database::ReportRefusal;

// Line DBID, DECISIONS 19: the `-db` part of `check` does not read a database with parsed
// Seq-ids of more than its bound (`Off::LargeParsed`, before any header is read: nrefuse_v5,
// whose headers refuse the reports, is not read above its bound), and reads it up to the
// bound, at one thread and in a pool, for outfmt 6 and 0.
#[test]
#[cfg_attr(
    not(feature = "parallel"),
    ignore = "needs threads: a build without the parallel feature has none"
)]
fn a_database_with_parsed_seqids_beyond_the_bound_keeps_its_reports_after_the_loop() {
    let open = |name: &str| {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/db_search/ids")
            .join(name);
        crate::algorithm::tblastx::database::TblastxDatabase::open(
            path.to_str().unwrap(),
            &crate::blastinput::ncbi_environment::ApplicationSettings::default(),
        )
        .unwrap_or_else(|error| panic!("{name}: {error:#}"))
    };
    let formats = [
        report::TblastxOutputFormat::Tabular,
        report::TblastxOutputFormat::Pairwise(crate::report::pairwise::PairwiseView::Pairwise),
    ];
    let answers = |pool: &SearchPool<'_>| {
        ["nucl_v5", "nrefuse_v5"].map(|name| {
            let database = open(name);
            let count = database.num_oids();
            formats.map(|format| {
                [count, count - 1].map(|max| database_check(&database, format, max, None, pool))
            })
        })
    };
    let headers = Err(Off::Database(ReportRefusal::Headers));
    let large = Err(Off::LargeParsed);
    let serial = answers(&SearchPool::x_serial());
    assert_eq!(
        serial,
        [
            [[Ok(()), large], [Ok(()), large]],
            [[headers, large], [headers, large]],
        ]
    );
    let parallel =
        crate::utils::threading::with_search_pool(4, "TBLASTX", |pool| Ok(answers(pool))).unwrap();
    assert_eq!(serial, parallel);
}
