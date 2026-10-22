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

/// The values of an XML2 document that outfmt 5 does not have (`s_SetBlastXMLParameters`,
/// blastxml2_format.cpp:514-555), read once per run from the program's resolved arguments.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Xml2Run {
    pub target: Xml2Target,
    /// Whether the search is gapped: `gap-open` and `gap-extend` are written only then (not
    /// for TBLASTX and `-ungapped`); megablast writes its 0 and 0.
    pub gapped: bool,
    /// `-comp_based_stats` when non-zero (`Parameters_cbs`), else 0 for none.
    pub cbs: i32,
    /// The query / database genetic code, for a translated query (BLASTX, TBLASTX) / subject
    /// (TBLASTN, TBLASTX); none otherwise.
    pub query_gencode: Option<i32>,
    pub db_gencode: Option<i32>,
}
