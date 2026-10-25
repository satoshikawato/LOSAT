//! The printings of a Seq-id that NCBI's formatters use: the FASTA form (`AsFastaString`,
//! `ref|XP_000001.1|`), the content label (`GetSeqIdString`, `XP_000001.1`), the formatter's
//! label (`CAlignFormatUtil::GetLabel`: a `general` id as its tag), the bare id
//! (`CAlignFormatUtil::GetBareId`), the id list of a definition line
//! (`CShowBlastDefline::GetSeqIdListString`) and the first GI (`FindGi`).
//!
//! Each function appends to `out`; the strings are the bytes the reader gives (`seqid::ber`).

use std::io::Write;

use super::model::{ObjectId, ParsedSeqId, PatentNumber, PatentSeqId, PdbSeqId, TextseqId};
use super::rank::{find_best, text_score};

fn number(out: &mut Vec<u8>, value: impl std::fmt::Display) {
    let _ = write!(out, "{value}");
}

// NCBI reference (598d8ae6): c++/src/objects/general/Object_id.cpp:202-215
// ```c++
// ostream& CObject_id::AsString(ostream &s) const
// {
//     switch ( Tparent::Which() ) {
//     case e_Id:
//         s << Tparent::GetId();
//         break;
//     case e_Str:
//         s << Tparent::GetStr(); // no Upcase() on output as per Ostell 7/2001 - Karl
//         break;
// ```
fn object_id(tag: &ObjectId, out: &mut Vec<u8>) {
    match tag {
        ObjectId::Id(id) => number(out, id),
        ObjectId::Str(text) => out.extend_from_slice(text),
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:572-596
// ```c++
// static const char* const s_TextId[CSeq_id::e_MaxChoice+1] =
// {   // FASTA_LONG formats
//     "???" , // not-set = ???
//     "lcl",  // local = lcl|integer or string
//     "bbs",  // gibbsq = bbs|integer
//     "bbm",  // gibbmt = bbm|integer
//     "gim",  // giim = gim|integer
//     "gb",   // genbank = gb|accession|locus
//     "emb",  // embl = emb|accession|locus
//     "pir",  // pir = pir|accession|name
//     "sp",   // swissprot = sp|accession|name *OR* tr|accession|name
//     "pat",  // patent = pat|country|patent number (string)|seq number (integer)
//             //     *OR* pgp|country|application number|seq number
//     "ref",  // other = ref|accession|name|release - changed from oth to ref
//     "gnl",  // general = gnl|database(string)|id (string or number)
//     "gi",   // gi = gi|integer
//     "dbj",  // ddbj = dbj|accession|locus
//     "prf",  // prf = prf|accession|name
//     "pdb",  // pdb = pdb|entry name (string)|chain id (char)
//     "tpg",  // tpg = tpg|accession|name
//     "tpe",  // tpe = tpe|accession|name
//     "tpd",  // tpd = tpd|accession|name
//     "gpp",  // gpipe = gpp|accession|name
//     "nat",  // named_annot_track = nat|accession|name
//     ""  // Placeholder for end of list
// };
// ```
/// The FASTA type tag of each alternative, in the order of `model::SEQ_ID_ALTERNATIVES`.
const FASTA_TAGS: [&[u8]; 20] = [
    b"lcl", b"bbs", b"bbm", b"gim", b"gb", b"emb", b"pir", b"sp", b"pat", b"ref", b"gnl", b"gi",
    b"dbj", b"prf", b"pdb", b"tpg", b"tpe", b"tpd", b"gpp", b"nat",
];

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2164-2186
// ```c++
// void CSeq_id::WriteAsFasta(ostream& out)
//     const
// {
//     unsigned the_type = Which();
//     if (the_type >= e_MaxChoice)  // New SeqId type
//         the_type = e_not_set;
//
//     if (IsPatent()  &&  !GetPatent().GetCit().GetId().IsNumber() ) {
//         const char pgp[] = "pat|";
//         out.write(pgp, sizeof(pgp) - 1);
//     } else if (IsSwissprot()  &&  GetSwissprot().IsSetRelease()
//                &&  GetSwissprot().GetRelease() == "unreviewed") {
//         const char tr[] = "tr|";
//         out.write(tr, sizeof(tr) - 1);
//     } else {
//         out.write(s_TextId[the_type], strlen(s_TextId[the_type]));
//         out.put('|');
//     }
//
//     x_WriteContentAsFasta(out);
// }
// ```
// The tag of a patent id with an application number is `pat|` too (the constant named
// `pgp` holds "pat|").
/// `CSeq_id::AsFastaString` (`WriteAsFasta`): `ref|XP_000001.1|`, `tr|A0A000|A0A000_YEAST`,
/// `pdb|1ABC|A`, `gnl|mydb|prot10`, `gi|1414`, ...
pub fn fasta(id: &ParsedSeqId, out: &mut Vec<u8>) {
    match id {
        ParsedSeqId::Swissprot(text) if text.release.as_deref() == Some(b"unreviewed") => {
            out.extend_from_slice(b"tr|");
        }
        _ => {
            out.extend_from_slice(FASTA_TAGS[usize::from(id.alternative())]);
            out.push(b'|');
        }
    }
    fasta_content(id, out);
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2188-2262
// ```c++
//     switch (the_type) {
//     case e_not_set:
//         break;
//     case e_Local:
//         GetLocal().AsString(out);
//         break;
//     case e_Gibbsq:
//         out << GetGibbsq();
//         break;
//     case e_Gibbmt:
//         out << GetGibbmt();
//         break;
//     case e_Giim:
//         out << (GetGiim().GetId());
//         break;
//     case e_Genbank:
//         GetGenbank().AsFastaString(out);
//         break;
//     ...
//     case e_Patent:
//         GetPatent().AsFastaString(out);
//         break;
//     ...
//     case e_General:
//         {
//             const CDbtag& dbt = GetGeneral();
//             out << (dbt.GetDb()) << '|';  // no Upcase per Ostell - Karl 7/2001
//             dbt.GetTag().AsString(out);
//         }
//         break;
//     case e_Gi:
//         out << GetGi();
//         break;
//     ...
//     case e_Pdb:
//         GetPdb().AsFastaString(out);
//         break;
//     ...
//     case e_Gpipe:
//         // don't suppress version after all
//         GetGpipe().AsFastaString(out /*, false */);
//         break;
// ```
/// `CSeq_id::x_WriteContentAsFasta`: the FASTA form without the type tag.
fn fasta_content(id: &ParsedSeqId, out: &mut Vec<u8>) {
    match id {
        ParsedSeqId::Local(tag) => object_id(tag, out),
        ParsedSeqId::Gibbsq(value) | ParsedSeqId::Gibbmt(value) => number(out, value),
        ParsedSeqId::Giim(giim) => number(out, giim.id),
        ParsedSeqId::Patent(patent) => patent_fasta(patent, out),
        ParsedSeqId::General(dbtag) => {
            out.extend_from_slice(&dbtag.db);
            out.push(b'|');
            object_id(&dbtag.tag, out);
        }
        ParsedSeqId::Gi(gi) => number(out, gi),
        ParsedSeqId::Pdb(pdb) => pdb_fasta(pdb, out),
        _ => textseq_fasta(id.textseq().expect("a Textseq-id alternative"), out),
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Textseq_id.cpp:228-247
// ```c++
// ostream& CTextseq_id::AsFastaString(ostream& s, bool allow_version) const
// {
//     if (IsSetAccession()) {
//         s << GetAccession(); // no Upcase per Ostell - Karl 7/2001
//         if (allow_version  &&  IsSetVersion()) {
//             int version = GetVersion();
//             if (version) {
//                 s.put('.');
//                 s << version;
//             }
//         }
//     }
//
//     s.put('|');
//     if ( IsSetName() ) {
//         s << GetName();  // no Upcase per Ostell - Karl 7/2001
//     }
//     return s;
// }
// ```
fn textseq_fasta(id: &TextseqId, out: &mut Vec<u8>) {
    if let Some(accession) = &id.accession {
        out.extend_from_slice(accession);
        if let Some(version) = id.version.filter(|&version| version != 0) {
            out.push(b'.');
            number(out, version);
        }
    }
    out.push(b'|');
    if let Some(name) = &id.name {
        out.extend_from_slice(name);
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/PDB_seq_id.cpp:94-108
// ```c++
// ostream& CPDB_seq_id::AsFastaString(ostream& s) const
// {
//     if (IsSetChain_id()) {
//             return s << GetMol().Get() << '|' << GetChain_id();
//     }
//     // no Upcase per Ostell - Karl 7/2001
//     char chain = (char) GetChain();
//
//     if (chain == '|') {
//         return s << GetMol().Get() << '|'; // historically |VB
//     } else if ( chain == '\0' ) {
//         return s << GetMol().Get() << "| ";
//     }
//     return s << GetMol().Get() << '|' << chain;
// }
// ```
fn pdb_fasta(id: &PdbSeqId, out: &mut Vec<u8>) {
    out.extend_from_slice(&id.mol);
    out.push(b'|');
    if let Some(chain_id) = &id.chain_id {
        out.extend_from_slice(chain_id);
        return;
    }
    match id.chain_or_default() as u8 {
        b'|' => {}
        0 => out.push(b' '),
        chain => out.push(chain),
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Patent_seq_id.cpp:71-88
// ```c++
// ostream& CPatent_seq_id::AsFastaString(ostream& s) const
// {
//     const CId_pat& idp = GetCit();
//
//     s << idp.GetCountry() << '|';  // no Upcase per Ostell - Karl 7/2001
//
//     const CId_pat_Base::C_Id& id = idp.GetId();
//     if ( id.IsNumber() ) {
//         s << id.GetNumber();
//     } else {
//         s << id.GetApp_number();
//         if ( idp.IsSetDoc_type() ) {
//             s << idp.GetDoc_type();
//         }
//     }
//     s << '|' << GetSeqid();
//     return s;
// }
// ```
fn patent_fasta(id: &PatentSeqId, out: &mut Vec<u8>) {
    out.extend_from_slice(&id.cit.country);
    out.push(b'|');
    match &id.cit.id {
        PatentNumber::Number(number) => out.extend_from_slice(number),
        PatentNumber::AppNumber(number) => {
            out.extend_from_slice(number);
            if let Some(doc_type) = &id.cit.doc_type {
                out.extend_from_slice(doc_type);
            }
        }
    }
    out.push(b'|');
    number(out, id.seqid);
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:1926-2038
// ```c++
// static inline
// void x_GetLabel_Content(const CSeq_id& id, string* label,
//                         CSeq_id::TLabelFlags flags, int* version)
// {
//     const CTextseq_id* tsid = id.GetTextseq_Id();
//     ...
//     //text id
//     if (tsid) {
//         string str;
//         if (tsid->IsSetAccession()) {
//             str = tsid->GetAccession();
//             NStr::ToUpper(str);
//         } else if (tsid->IsSetName()) {
//             str = tsid->GetName();
//         }
//         ...
//         if ( !str.empty() ) {
//             if ( (flags & CSeq_id::fLabel_Version)  &&  tsid->IsSetVersion()) {
//                 str += "." + NStr::IntToString(tsid->GetVersion());
//             }
//         }
//         *label += str;
//
//     } else { //non-text id
//         switch (id.Which()) {
//         ...
//         case CSeq_id::e_Local:
//             {{
//                 const CObject_id& oid = id.GetLocal();
//                 if (oid.IsId()) {
//                     *label += NStr::IntToString(oid.GetId());
//                 } else if (oid.IsStr()) {
//                     *label += oid.GetStr();
//                 }
//             }}
//             break;
//
//         case CSeq_id::e_Gibbsq:
//             *label += NStr::IntToString(id.GetGibbsq());
//             break;
//
//         case CSeq_id::e_Gibbmt:
//             *label += NStr::IntToString(id.GetGibbmt());
//             break;
//
//         case CSeq_id::e_Giim:
//             *label += NStr::IntToString(id.GetGiim().GetId());
//             break;
//
//         case CSeq_id::e_General:
//             {{
//                 const CDbtag& dbt = id.GetGeneral();
//                 if ((flags & CSeq_id::fLabel_GeneralDbIsContent) == 0) {
//                     *label += dbt.GetDb() + ':';
//                 }
//                 if (dbt.GetTag().IsId()) {
//                     *label += NStr::IntToString(dbt.GetTag().GetId());
//                 } else if (dbt.GetTag().IsStr()) {
//                     *label += dbt.GetTag().GetStr();
//                 }
//             }}
//             break;
//
//         case CSeq_id::e_Patent:
//             {{
//                 const CId_pat& idp = id.GetPatent().GetCit();
//                 *label += idp.GetCountry();
//                 // *label += '|';
//                 *label += (idp.GetId().IsNumber() ?
//                            idp.GetId().GetNumber() :
//                            idp.GetId().GetApp_number());
//                 *label += '_'; // |
//                 *label += NStr::IntToString(id.GetPatent().GetSeqid());
//             }}
//             break;
//
//         case CSeq_id::e_Gi:
//             *label += NStr::NumericToString(id.GetGi());
//             break;
//
//         case CSeq_id::e_Pdb:
//             {{
//                 const CPDB_seq_id& pid = id.GetPdb();
//                 *label += pid.GetMol().Get();
//                 if (pid.IsSetChain_id()) {
//                     *label += '_';
//                     *label += pid.GetChain_id();
//                 } else if (pid.IsSetChain()) {
//                     unsigned char chain = static_cast<unsigned char>(pid.GetChain());
//                     if (chain > ' ') {
//                         *label += '_';
//                         ...
//                         *label += static_cast<char>(chain);
//                     }
//                 }
//             }}
//             break;
// ```
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2145-2154
// ```c++
// string CSeq_id::GetSeqIdString(bool with_version) const
// {
//     string label;
//     TLabelFlags flags = 0;
//     if (with_version) {
//         flags |= fLabel_Version;
//     }
//     GetLabel(&label, eContent, flags);
//     return label;
// }
// ```
/// `CSeq_id::GetSeqIdString(with_version)` (`GetLabel(eContent)` with the flag
/// `fLabel_Version` or none): the upper-cased accession (else the name) of a Textseq-id with
/// `.version` when asked and stored, the number of a `gi`, `db:tag` of a `general` id,
/// `mol_chain` of a `pdb` id, `country number_seqid` of a patent.
pub fn content(id: &ParsedSeqId, with_version: bool, out: &mut Vec<u8>) {
    if let Some(text) = id.textseq() {
        let start = out.len();
        if let Some(accession) = &text.accession {
            out.extend(accession.iter().map(u8::to_ascii_uppercase));
        } else if let Some(name) = &text.name {
            out.extend_from_slice(name);
        }
        if out.len() > start && with_version {
            if let Some(version) = text.version {
                out.push(b'.');
                number(out, version);
            }
        }
        return;
    }
    match id {
        ParsedSeqId::Local(tag) => object_id(tag, out),
        ParsedSeqId::Gibbsq(value) | ParsedSeqId::Gibbmt(value) => number(out, value),
        ParsedSeqId::Giim(giim) => number(out, giim.id),
        ParsedSeqId::General(dbtag) => {
            out.extend_from_slice(&dbtag.db);
            out.push(b':');
            object_id(&dbtag.tag, out);
        }
        ParsedSeqId::Patent(patent) => {
            out.extend_from_slice(&patent.cit.country);
            match &patent.cit.id {
                PatentNumber::Number(number) | PatentNumber::AppNumber(number) => {
                    out.extend_from_slice(number)
                }
            }
            out.push(b'_');
            number(out, patent.seqid);
        }
        ParsedSeqId::Gi(gi) => number(out, gi),
        ParsedSeqId::Pdb(pdb) => {
            out.extend_from_slice(&pdb.mol);
            if let Some(chain_id) = &pdb.chain_id {
                out.push(b'_');
                out.extend_from_slice(chain_id);
            } else if let Some(chain) = pdb.chain {
                let chain = chain as u8;
                if chain > b' ' {
                    out.push(b'_');
                    out.push(chain);
                }
            }
        }
        _ => unreachable!("Textseq-id alternatives are handled above"),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:907-930
// ```c++
// string CAlignFormatUtil::GetGnlID(const CDbtag& dtg)
// {
//    string retval = NcbiEmptyString;
//
//    if(dtg.GetTag().IsId())
//      retval = NStr::IntToString(dtg.GetTag().GetId());
//    else
//      retval = dtg.GetTag().GetStr();
//
//    return retval;
// }
//
// string CAlignFormatUtil::GetLabel(CConstRef<CSeq_id> id,bool with_version)
// {
//     string retval = "";
//     if (id->Which() == CSeq_id::e_General){
//         const CDbtag& dtg = id->GetGeneral();
//         retval = CAlignFormatUtil::GetGnlID(dtg);
//     }
//     if (retval == "")
//       retval = id->GetSeqIdString(with_version);
//
//     return retval;
// }
// ```
/// `CAlignFormatUtil::GetLabel(id, with_version)`: a `general` id as its tag alone, every
/// other id (and a `general` id with an empty tag) as `content`.
pub fn format_label(id: &ParsedSeqId, with_version: bool, out: &mut Vec<u8>) {
    if let ParsedSeqId::General(dbtag) = id {
        let start = out.len();
        object_id(&dbtag.tag, out);
        if out.len() > start {
            return;
        }
    }
    content(id, with_version, out);
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4944-4956
// ```c++
// string CAlignFormatUtil::GetBareId(const CSeq_id& id)
// {
//     string retval;
//
//     if (id.IsGi() || id.IsPrf() || id.IsPir()) {
//         retval = id.AsFastaString();
//     }
//     else {
//         retval = id.GetSeqIdString(true);
//     }
//
//     return retval;
// }
// ```
/// `CAlignFormatUtil::GetBareId`: the FASTA form of a `gi`, `prf` or `pir` id, the content
/// with its version for every other id.
pub fn bare(id: &ParsedSeqId, out: &mut Vec<u8>) {
    match id {
        ParsedSeqId::Gi(_) | ParsedSeqId::Prf(_) | ParsedSeqId::Pir(_) => fasta(id, out),
        _ => content(id, true, out),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:156-189
// ```c++
// string
// CShowBlastDefline::GetSeqIdListString(const list<CRef<objects::CSeq_id> >& id,
//                                       bool show_gi)
// {
//     string id_string = NcbiEmptyString;
//     bool found_gi = false;
//
//     CRef<CSeq_id> best_id = FindBestChoice(id, CSeq_id::Score);
//
//     if (show_gi) {
//         ITERATE(list<CRef<CSeq_id> >, itr, id) {
//              if ((*itr)->IsGi()) {
//                 id_string += (*itr)->AsFastaString();
//                 found_gi = true;
//                 break;
//              }
//         }
//     }
//
//     if (best_id.NotEmpty()  &&  !best_id->IsGi() ) {
//          if (found_gi)
//              id_string += "|";
//
//          if (best_id->IsLocal()) {
//              string id_token;
//              best_id->GetLabel(&id_token, CSeq_id::eContent, 0);
//              id_string += id_token;
//          }
//          else
//              id_string += best_id->AsFastaString();
//     }
//
//     return id_string;
// }
// ```
/// `CShowBlastDefline::GetSeqIdListString(ids, show_gi)`: the first GI in FASTA form (with
/// `show_gi`), then `|` and the FASTA form of the `Score`-best id unless that is the GI (a
/// local id as its content, without `lcl|`): `gi|123456|ref|NP_000002.1|`, `myprot9`.
pub fn seqid_list_string(ids: &[ParsedSeqId], show_gi: bool, out: &mut Vec<u8>) {
    let best = find_best(ids, text_score);
    let mut found_gi = false;
    if show_gi {
        if let Some(gi) = ids.iter().find(|id| id.gi().is_some()) {
            fasta(gi, out);
            found_gi = true;
        }
    }
    if let Some(best) = best.filter(|best| best.gi().is_none()) {
        if found_gi {
            out.push(b'|');
        }
        if let ParsedSeqId::Local(_) = best {
            content(best, false, out);
        } else {
            fasta(best, out);
        }
    }
}

// NCBI reference (598d8ae6): c++/include/objects/seqloc/Seq_id.hpp:1041-1047
// ```c++
// /// Return gi from id list if exists, return 0 otherwise
// template<class container>
// TGi FindGi(const container& ids)
// {
//     CConstRef<CSeq_id> id = GetSeq_idByType(ids, CSeq_id::e_Gi);
//     return id ? id->GetGi() : ZERO_GI;
// }
// ```
/// `FindGi`: the first GI of `ids`, 0 for none.
pub fn first_gi<'a>(ids: impl IntoIterator<Item = &'a ParsedSeqId>) -> i64 {
    ids.into_iter().find_map(ParsedSeqId::gi).unwrap_or(0)
}
