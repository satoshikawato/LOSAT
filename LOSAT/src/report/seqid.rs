//! NCBI's Seq-id of the sequences that LOSAT reports, and its printings.
//!
//! A reported sequence has one of three Seq-ids: the local ID that NCBI's FASTA reader
//! generates (`Query_N`, `Subject_N`), the ordinal ID of a sequence of a database without
//! parsed Seq-ids (`gnl|BL_ORD_ID|<oid>`), or a Seq-id that a database built with
//! `-parse_seqids` stores or `-parse_deflines` read from a query's defline (`SeqId::Db`, chosen
//! and printed by `crate::seqid`). The serial
//! formats (outfmt 8, 9, 11, 12) write the Seq-id object itself (`SeqId::write`); the other
//! printings (the title's first word, `gnl|BL_ORD_ID|0`, `BL_ORD_ID:0`, ...) join this type
//! with the formats that use them.
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
use crate::seqid::model::{Date, ObjectId, PatentNumber, PDB_DEFAULT_CHAIN};
use crate::seqid::ParsedSeqId;

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
    /// A Seq-id of a database sequence with parsed Seq-ids: the Seq-align id that
    /// `seqid::DeflineSet::alignment_id` chooses (or the id a format shows).
    Db(&'a ParsedSeqId),
}

impl SeqId<'_> {
    /// Appends the Seq-id as the SAM formatter prints it: `GetSeqIdString(true)` of the
    /// `fSAM_PlainSeqIds` flag (the local string, or `BL_ORD_ID:<oid>` for a database
    /// sequence). `AdjustSeqIdType` finds no accession for either id, so it keeps them. A
    /// parsed Seq-id prints `GetSeqIdString(true)` (`seqid::label::content`) of the id it
    /// holds: the caller passes the id that NCBI's `GetId(eGetId_ForceAcc)` and
    /// `AdjustSeqIdType` give.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:111-115
    /// ```c++
    /// string CSAM_CIGAR_Formatter::x_GetSeqIdString(const CSeq_id& id, TFlags flags) const
    /// {
    ///     return (flags & CSAM_Formatter::fSAM_PlainSeqIds) ?
    ///         id.GetSeqIdString(true) : id.AsFastaString();
    /// }
    /// ```
    pub fn sam_name(&self, out: &mut Vec<u8>) {
        match *self {
            SeqId::Local(id) => out.extend_from_slice(id.as_bytes()),
            SeqId::BlOrdId(oid) => {
                out.extend_from_slice(b"BL_ORD_ID:");
                out.extend_from_slice(oid.to_string().as_bytes());
            }
            SeqId::Db(id) => crate::seqid::label::content(id, true, out),
        }
    }

    /// Writes the Seq-id object (a CHOICE) through `enc`, in the member order of the ASN.1
    /// modules quoted above: `local` (variant 0) holds an Object-id `str` (variant 1);
    /// `general` (variant 10) holds a Dbtag (`db` member 0, `tag` member 1 = Object-id `id`,
    /// variant 0). A parsed Seq-id is written by `write_parsed`.
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
            SeqId::Db(id) => write_parsed(id, enc),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/serial/member.cpp:444-454
// ```c++
//     else {
//         // have set flag
//         getConstFunc = &TFunc::GetConstSimpleMember;
//         getFunc = &TFunc::GetWithSetFlagMember;
//         readFuncs.first = &TFunc::ReadWithSetFlagMember;
//         writeFunc = &TFunc::WriteWithSetFlagMember;
//         if ( (GetDefault() || Nillable()) && GetId().HaveNoPrefix()) {
//             readFuncs.first = &TFunc::ReadWithDefaultMemberX;
//             writeFunc = &TFunc::WriteWithDefaultMemberX;
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/serial/member.cpp:1015-1027
// ```c++
// void CMemberInfoFunctions::WriteWithSetFlagMember(CObjectOStream& out,
//                                                   const CMemberInfo* memberInfo,
//                                                   TConstObjectPtr classPtr)
// {
//     _ASSERT(!memberInfo->CanBeDelayed());
//     _ASSERT(memberInfo->HaveSetFlag());
//     if (out.IsWritingDefaultValuesEnforced() &&  memberInfo->GetDefault()) {
//         goto do_write;
//     }
//     if ( memberInfo->GetSetFlagNo(classPtr) ) {
//         if (memberInfo->Optional()) {
//             return;
//         }
// ```
// The members of a type from an ASN.1 module have a prefix (`HaveNoPrefix` is false), so
// they all take `WriteWithSetFlagMember`, and a member with a DEFAULT is optional
// (`CMemberInfo::Optional`, serial/impl/member.inl:42-51). `PDB-seq-id.chain` is therefore
// written when it was read, and else only by a stream that enforces the defaults (the JSON
// stream, `Encoder::writes_defaults`), with its DEFAULT 32: NCBI's outfmt 12 writes
// `"chain": 32` for a pdb id stored without one, outfmt 8 nothing.
/// Writes a parsed Seq-id (a CHOICE) through `enc`, in the member order of the ASN.1 modules
/// quoted in `seqid::model`, with the calls of `SeqId::write` for `local` and `general`. An
/// absent OPTIONAL member is skipped.
fn write_parsed<E: Encoder>(id: &ParsedSeqId, enc: &mut E) {
    enc.variant("Seq-id", id.alternative_name(), id.alternative());
    match id {
        ParsedSeqId::Local(tag) => write_object_id(enc, tag),
        ParsedSeqId::Gibbsq(value) | ParsedSeqId::Gibbmt(value) => enc.int(*value),
        ParsedSeqId::Giim(id) => {
            enc.begin_sequence("Giimport-id");
            enc.member("id", 0);
            enc.int(id.id);
            enc.end_member();
            optional_str(enc, "db", 1, id.db.as_deref());
            optional_str(enc, "release", 2, id.release.as_deref());
            enc.end_sequence();
        }
        ParsedSeqId::Patent(id) => {
            enc.begin_sequence("Patent-seq-id");
            enc.member("seqid", 0);
            enc.int(id.seqid);
            enc.end_member();
            enc.member("cit", 1);
            enc.begin_sequence("Id-pat");
            enc.member("country", 0);
            enc.visible_str(&id.cit.country);
            enc.end_member();
            enc.member("id", 1);
            match &id.cit.id {
                PatentNumber::Number(number) => {
                    enc.variant("Id-pat_id", "number", 0);
                    enc.visible_str(number);
                }
                PatentNumber::AppNumber(number) => {
                    enc.variant("Id-pat_id", "app-number", 1);
                    enc.visible_str(number);
                }
            }
            enc.end_variant();
            enc.end_member();
            optional_str(enc, "doc-type", 2, id.cit.doc_type.as_deref());
            enc.end_sequence();
            enc.end_member();
            enc.end_sequence();
        }
        ParsedSeqId::General(dbtag) => {
            enc.begin_sequence("Dbtag");
            enc.member("db", 0);
            enc.visible_str(&dbtag.db);
            enc.end_member();
            enc.member("tag", 1);
            write_object_id(enc, &dbtag.tag);
            enc.end_member();
            enc.end_sequence();
        }
        ParsedSeqId::Gi(gi) => enc.big_int(*gi),
        ParsedSeqId::Pdb(id) => {
            enc.begin_sequence("PDB-seq-id");
            enc.member("mol", 0);
            // `mol PDB-mol-id` (seqloc.asn), a named VisibleString (`begin_named`).
            enc.begin_named("PDB-mol-id");
            enc.visible_str(&id.mol);
            enc.end_named();
            enc.end_member();
            let chain = match id.chain {
                Some(chain) => Some(chain),
                None if enc.writes_defaults() => Some(PDB_DEFAULT_CHAIN),
                None => None,
            };
            optional_int(enc, "chain", 1, chain);
            if let Some(rel) = &id.rel {
                enc.member("rel", 2);
                write_date(enc, rel);
                enc.end_member();
            }
            optional_str(enc, "chain-id", 3, id.chain_id.as_deref());
            enc.end_sequence();
        }
        ParsedSeqId::Genbank(id)
        | ParsedSeqId::Embl(id)
        | ParsedSeqId::Pir(id)
        | ParsedSeqId::Swissprot(id)
        | ParsedSeqId::Other(id)
        | ParsedSeqId::Ddbj(id)
        | ParsedSeqId::Prf(id)
        | ParsedSeqId::Tpg(id)
        | ParsedSeqId::Tpe(id)
        | ParsedSeqId::Tpd(id)
        | ParsedSeqId::Gpipe(id)
        | ParsedSeqId::NamedAnnotTrack(id) => {
            enc.begin_sequence("Textseq-id");
            optional_str(enc, "name", 0, id.name.as_deref());
            optional_str(enc, "accession", 1, id.accession.as_deref());
            optional_str(enc, "release", 2, id.release.as_deref());
            optional_int(enc, "version", 3, id.version);
            enc.end_sequence();
        }
    }
    enc.end_variant();
}

fn write_object_id<E: Encoder>(enc: &mut E, tag: &ObjectId) {
    match tag {
        ObjectId::Id(id) => {
            enc.variant("Object-id", "id", 0);
            enc.int(*id);
        }
        ObjectId::Str(text) => {
            enc.variant("Object-id", "str", 1);
            enc.visible_str(text);
        }
    }
    enc.end_variant();
}

fn write_date<E: Encoder>(enc: &mut E, date: &Date) {
    match date {
        Date::Str(text) => {
            enc.variant("Date", "str", 0);
            enc.visible_str(text);
        }
        Date::Std(std) => {
            enc.variant("Date", "std", 1);
            enc.begin_sequence("Date-std");
            enc.member("year", 0);
            enc.int(std.year);
            enc.end_member();
            optional_int(enc, "month", 1, std.month);
            optional_int(enc, "day", 2, std.day);
            optional_str(enc, "season", 3, std.season.as_deref());
            optional_int(enc, "hour", 4, std.hour);
            optional_int(enc, "minute", 5, std.minute);
            optional_int(enc, "second", 6, std.second);
            enc.end_sequence();
        }
    }
    enc.end_variant();
}

fn optional_str<E: Encoder>(enc: &mut E, name: &str, tag: u8, value: Option<&[u8]>) {
    if let Some(value) = value {
        enc.member(name, tag);
        enc.visible_str(value);
        enc.end_member();
    }
}

fn optional_int<E: Encoder>(enc: &mut E, name: &str, tag: u8, value: Option<i32>) {
    if let Some(value) = value {
        enc.member(name, tag);
        enc.int(value);
        enc.end_member();
    }
}
