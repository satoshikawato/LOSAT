//! What the XML2 `Report` and `Parameters` need beyond outfmt 5's `XmlRun` (DESIGN_OF5 2.1).

/// What a search is run against: the `-db` name (`XmlRun::db_name`) or a `-subject` set of
/// `Subject_1 .. Subject_N`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Xml2Target {
    /// `-db`: `Target` is the `db` variant; `Statistics` has `db-num` and `db-len`; a `Hit`
    /// has an accession and a taxid.
    Database,
    /// `-subject`: `Target` is the `subjects` variant with this many `Subject_N` strings, one
    /// for every record of the file, empty records too (blastxml2_format.cpp:668-669,
    /// data4xml2format.cpp:174-185); the results are `bl2seq`.
    Subjects(usize),
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:121, :196
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_cxx.cpp:862-867
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2362-2377, :2376-2377
/// NCBI's `m_IsBl2Seq` or `m_IsIterative`: whether the XML2 and JSON documents of a search hold each
/// query's report back until a result of a query whose Seq-id does not match arrives (or the
/// epilog): a search of subject records (`subject_search`: `-subject`, no database name) or of an
/// iterative program (`iterative`: psiblast, tblastn `-in_pssm`). The one rule behind the writer's
/// hold (`Xml2Target::holds_back`), the place where a matrix name NCBI's XML does not know fails
/// (`blast_xml::matrix::xml_formatting_failure`) and BLASTX's checks of the reports' titles before
/// the search (`blastx::blast_xml::scan_queries`).
pub fn reports_held(subject_search: bool, iterative: bool) -> bool {
    subject_search || iterative
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xml2format.cpp:180-183
// The Seq-id of a subject's Seq-loc is the `BestRank` best of its ids
// (blast_fasta_input.cpp:462-463).
/// The `subjects` strings of the records `records` (their ids under `-parse_deflines`): the bare id
/// of each record's best id; `None` for records read without the option.
pub fn subject_names<'a>(
    mut records: impl Iterator<
        Item = Option<&'a crate::blastinput::fasta_reader::believed::BelievedIds>,
    >,
) -> Option<std::sync::Arc<[String]>> {
    let first = records.next()?;
    let first = first?;
    let name = |ids: &crate::blastinput::fasta_reader::believed::BelievedIds| {
        let mut out = Vec::new();
        crate::seqid::label::bare(&ids.best, &mut out);
        String::from_utf8_lossy(&out).into_owned()
    };
    let mut names = vec![name(first)];
    for ids in records {
        names.push(ids.map_or_else(String::new, name));
    }
    Some(names.into())
}

/// The values of an XML2 document that outfmt 5 does not have (`s_SetBlastXMLParameters`,
/// blastxml2_format.cpp:514-555), read once per run from the program's resolved arguments.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Xml2Run {
    pub target: Xml2Target,
    /// The `subjects` strings of `-subject` records whose ids `-parse_deflines` read
    /// (`subject_names`); `None`: `Subject_N`.
    pub subject_names: Option<std::sync::Arc<[String]>>,
    /// Whether the search is gapped: `gap-open` and `gap-extend` are written only then (not
    /// for TBLASTX and `-ungapped`); megablast writes its 0 and 0.
    pub gapped: bool,
    /// `-comp_based_stats` when non-zero (`Parameters_cbs`), else 0 for none.
    pub cbs: i32,
    /// The query / database genetic code, for a translated query (BLASTX, TBLASTX) / subject
    /// (TBLASTN, TBLASTX); none otherwise.
    pub query_gencode: Option<i32>,
    pub db_gencode: Option<i32>,
    /// The program type of the search where its program does not name it
    /// (`report::record::RunRecord::program_type`; `None`: the program's own): an iterative
    /// one's reports are held (`holds_back`, `iterative`), and the frames written are the
    /// sequences' that the type translates.
    pub program_type: Option<crate::report::record::ProgramType>,
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:509-532
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2295-2314
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:131-153, 585-595
    /// Whether NCBI's formatter has dropped the queries from its scope before it makes the held
    /// reports (`query_dropped`): psiblast's FASTA queries against a database, whose loop resets
    /// the scope after each query's rounds (`ResetScopeHistory`), before the next query's first
    /// result or the epilog makes the held report. Its `Search` then names the query by its
    /// Seq-loc's id (`AsFastaString`) without a title, and each `Hsp` is the catch's: its Seq-align
    /// scores without coordinates or rows (`blast_xml::writer::write_caught_hsp`). A `-subject`
    /// search never resets its scope, and psiblast's PSSM query (`-in_pssm`, `-in_msa`) is not
    /// reset either (`query_dropped`).
    pub query_dropped: bool,
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:509-562
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2295-2314
/// Whether the XML2 and JSON reports of a search are made after NCBI's formatter dropped its
/// queries from the scope (`Xml2Run::query_dropped`): an iterative search (`iterative`) of FASTA
/// queries (`fasta_queries`; not psiblast's PSSM query of `-in_pssm` or `-in_msa`, which is added
/// to the scope once, outside the loop of queries that resets it) against a database
/// (`subject_search` false: `ResetScopeHistory` does nothing for a search of subject records).
/// The one rule of the programs' `Xml2Run`s.
pub fn queries_dropped(subject_search: bool, iterative: bool, fasta_queries: bool) -> bool {
    iterative && fasta_queries && !subject_search
}

impl Xml2Run {
    /// Whether the search is an iterative program's (NCBI's `m_IsIterative`,
    /// `ProgramType::is_iterative`: psiblast, tblastn `-in_pssm`): its reports are held
    /// (`holds_back`), and the results of a database search are its rounds, `iterations`
    /// (`s_FillBlastOutput`; a `-subject` search's are `bl2seq` all the same).
    pub fn iterative(&self) -> bool {
        self.program_type
            .is_some_and(crate::report::record::ProgramType::is_iterative)
    }

    /// Whether the documents of the search hold each query's report back until the next query's
    /// arrives (`reports_held`; `Xml2Document`): a `-subject` search, or an iterative one.
    pub fn holds_back(&self) -> bool {
        reports_held(
            matches!(self.target, Xml2Target::Subjects(_)),
            self.iterative(),
        )
    }
}
