//! NCBI's Seq-id of the sequences that LOSAT reports, and its printings.
//!
//! LOSAT reads queries and `-subject` subjects without `-parse_deflines` (refused), and
//! BLAST databases without parsed Seq-ids (`seqdb/header.rs`), so a reported sequence has one
//! of two Seq-ids: the local ID that NCBI's FASTA reader generates (`Query_N`, `Subject_N`),
//! or the ordinal ID of a database sequence (`gnl|BL_ORD_ID|<oid>`). The serial formats
//! (outfmt 8, 9, 11, 12) write the Seq-id object itself (`SeqId::write`); the other printings
//! (the title's first word, `gnl|BL_ORD_ID|0`, `BL_ORD_ID:0`, ...) join this type with the
//! formats that use them.
//!
//! NCBI reference (598d8ae6): c++/src/objects/seqloc/seqloc.asn:24-36
//! ```text
//! Seq-id ::= CHOICE {
//!     local Object-id ,            -- local use
//!     ...
//!     general Dbtag ,              -- for other databases
//! ```
//! NCBI reference (598d8ae6): c++/src/objects/general/general.asn:71-80
//! ```text
//! Dbtag ::= SEQUENCE {
//!     db VisibleString ,          -- name of database or system
//!     tag Object-id }         -- appropriate tag
//! ...
//! Object-id ::= CHOICE {
//!     id INTEGER ,
//!     str VisibleString }
//! ```

use crate::report::serial::encoder::Encoder;

/// The Seq-id of a reported sequence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SeqId<'a> {
    /// `local str "<id>"`: the local ID of a FASTA record (`Query_N`, `Subject_N`,
    /// `FastaRecord::local_id`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:364-367
    /// ```c++
    ///     CRef<CSeqIdGenerator> idgen
    ///         (new CSeqIdGenerator(m_Config.GetLocalIdCounterInitValue(),
    ///                              m_Config.GetLocalIdPrefix()));
    ///     m_InputReader->SetIDGenerator(*idgen);
    /// ```
    /// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input.hpp:189-191
    /// ```c++
    ///     void SetQueryLocalIdMode() {m_LocalIdPrefix = "Query_";}
    ///     ...
    ///     void SetSubjectLocalIdMode() {m_LocalIdPrefix = "Subject_";}
    /// ```
    Local(&'a str),
    /// `general { db "BL_ORD_ID", tag id <oid> }`: a sequence of a BLAST database without
    /// parsed Seq-ids, by its 0-based ordinal ID in the database.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:1782-1787
    /// ```c++
    ///     if (!parse_ids) {
    ///
    ///         // Generate an BL_ORD_ID in case no parse is needed
    ///         CRef<CSeq_id> gnl_id(new CSeq_id());
    ///         gnl_id->SetGeneral().SetDb("BL_ORD_ID");
    ///         gnl_id->SetGeneral().SetTag().SetId(0);  // will be filled later
    /// ```
    /// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2221-2238
    /// ```c++
    ///     if (adjust_oids && bdls.NotEmpty() && m_VolStart) {
    ///         ...
    ///                 if (seqid.Which() == CSeq_id::e_General) {
    ///                     CDbtag & dbt = seqid.SetGeneral();
    ///
    ///                     if (dbt.GetDb() == "BL_ORD_ID") {
    ///                         int vol_oid = dbt.GetTag().GetId();
    ///                         dbt.SetTag().SetId(m_VolStart + vol_oid);
    /// ```
    /// The tag is the sequence's ordinal ID in the opened database (its volume's start plus
    /// its position in the volume).
    BlOrdId(u32),
}

impl SeqId<'_> {
    /// Writes the Seq-id object (a CHOICE) through `enc`, in the member order of the ASN.1
    /// modules quoted above: `local` (variant 0) holds an Object-id `str` (variant 1);
    /// `general` (variant 10) holds a Dbtag (`db` member 0, `tag` member 1 = Object-id `id`,
    /// variant 0).
    pub fn write<E: Encoder>(&self, enc: &mut E) {
        match *self {
            SeqId::Local(id) => {
                enc.variant("Seq-id", "local", 0);
                enc.variant("Object-id", "str", 1);
                enc.visible_str(id.as_bytes());
                enc.end_variant();
                enc.end_variant();
            }
            SeqId::BlOrdId(oid) => {
                enc.variant("Seq-id", "general", 10);
                enc.begin_sequence("Dbtag");
                enc.member("db", 0);
                enc.visible_str(b"BL_ORD_ID");
                enc.end_member();
                enc.member("tag", 1);
                enc.variant("Object-id", "id", 0);
                // An ordinal ID is an `Int4` in NCBI (`int oid`); a database has fewer than
                // 2^31 sequences.
                enc.int(oid as i32);
                enc.end_variant();
                enc.end_member();
                enc.end_sequence();
                enc.end_variant();
            }
        }
    }
}
