//! The one writer of NCBI's `Bioseq` (NCBI-Sequence) for the Bioseqs that LOSAT makes: the FASTA
//! record of a BLAST archive (`report::archive::bioseq`, NCBI's `CFastaReader`) and a database
//! sequence (blastdbcmd `-outfmt %b`, NCBI's `CSeqDB::GetBioseq`). Each caller gives the ids, the
//! descriptors and the instance its Bioseq holds; this module makes the object stream's calls
//! (`Encoder`) for them, members in module order and only those that are set.
//!
//! NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:36-73
//! NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:262-307
//! NCBI reference (598d8ae6): c++/src/objects/seq/seq.asn:346-358
//! NCBI reference (598d8ae6): c++/src/objects/general/general.asn:71-80
//! NCBI reference (598d8ae6): c++/src/objects/general/general.asn:128-148
//! NCBI reference (598d8ae6): c++/src/objects/seqfeat/seqfeat.asn:1060-1066
//! NCBI reference (598d8ae6): c++/src/objects/seqfeat/seqfeat.asn:1174-1214

use crate::report::seqid::SeqId;
use crate::report::serial::encoder::Encoder;

/// A `Seqdesc` of the Bioseq's `descr`.
#[derive(Clone, Copy, Debug)]
pub enum Seqdesc<'a> {
    /// `title`.
    Title(&'a [u8]),
    /// `user`: a `User-object` of one field.
    User(UserObject<'a>),
    /// `source`: a `BioSource` whose `org` is `OrgSource`.
    Source(OrgSource<'a>),
}

/// A `User-object` whose `type` is the string `type_str` and whose `data` is one `User-field`
/// labelled with the string `label`.
#[derive(Clone, Copy, Debug)]
pub struct UserObject<'a> {
    pub type_str: &'a [u8],
    pub label: &'a [u8],
    /// The field's `num`.
    pub num: Option<i32>,
    pub data: UserData<'a>,
}

/// The `data` of a `User-field`.
#[derive(Clone, Copy, Debug)]
pub enum UserData<'a> {
    /// `str`.
    Str(&'a [u8]),
    /// `oss` of one OCTET STRING.
    Oss(&'a [u8]),
}

/// The `Org-ref` of a `BioSource` that holds no more than a `taxon` tag and the names of the
/// taxonomy database (`CSeqDBVol::x_GetTaxonomy`).
#[derive(Clone, Copy, Debug)]
pub struct OrgSource<'a> {
    pub taxname: Option<&'a [u8]>,
    pub common: Option<&'a [u8]>,
    /// The `id` tag of the `Dbtag` `taxon`.
    pub taxid: i32,
}

/// The `Seq-data` of a `Seq-inst`.
#[derive(Clone, Copy, Debug)]
pub enum SeqData<'a> {
    Ncbi2na(&'a [u8]),
    Ncbi4na(&'a [u8]),
    /// The letters (a StringStore).
    Ncbieaa(&'a [u8]),
    Ncbistdaa(&'a [u8]),
}

/// The `Seq-inst` of a Bioseq: `repr` `raw`, or `virtual` without data; `mol` `aa` or `na`; the
/// `length`.
#[derive(Clone, Copy, Debug)]
pub struct SeqInst<'a> {
    pub protein: bool,
    pub length: u32,
    /// `seq-data`; `None` for a `virtual` instance.
    pub data: Option<SeqData<'a>>,
}

/// Writes the Bioseq of the ids `ids`, the descriptors `descr` (in their order; no `descr` member
/// when there is none) and the instance `inst`: as the element it is in a `SEQUENCE OF` / `SET OF`
/// (the caller has called `element()`), the value of a member, or a document the caller began.
pub fn write_bioseq<E: Encoder>(
    enc: &mut E,
    ids: &[SeqId<'_>],
    descr: &[Seqdesc<'_>],
    inst: &SeqInst<'_>,
) {
    enc.begin_sequence("Bioseq");
    enc.member("id", 0);
    enc.begin_list(true);
    for id in ids {
        enc.element();
        id.write(enc);
    }
    enc.end_list();
    enc.end_member();

    if !descr.is_empty() {
        enc.member("descr", 1);
        enc.begin_named("Seq-descr");
        enc.begin_list(true);
        for desc in descr {
            enc.element();
            write_seqdesc(enc, desc);
        }
        enc.end_list();
        enc.end_named();
        enc.end_member();
    }

    enc.member("inst", 2);
    write_seq_inst(enc, inst);
    enc.end_member();
    enc.end_sequence();
}

fn write_seqdesc<E: Encoder>(enc: &mut E, desc: &Seqdesc<'_>) {
    match desc {
        Seqdesc::Title(title) => {
            enc.variant("Seqdesc", "title", 4);
            enc.visible_str(title);
            enc.end_variant();
        }
        Seqdesc::User(user) => {
            enc.variant("Seqdesc", "user", 13);
            write_user_object(enc, user);
            enc.end_variant();
        }
        Seqdesc::Source(source) => {
            enc.variant("Seqdesc", "source", 22);
            enc.begin_sequence("BioSource");
            enc.member("org", 2);
            write_org_ref(enc, source);
            enc.end_member();
            enc.end_sequence();
            enc.end_variant();
        }
    }
}

fn write_user_object<E: Encoder>(enc: &mut E, user: &UserObject<'_>) {
    enc.begin_sequence("User-object");
    enc.member("type", 1);
    enc.variant("Object-id", "str", 1);
    enc.visible_str(user.type_str);
    enc.end_variant();
    enc.end_member();
    enc.member("data", 2);
    enc.begin_list(false);
    enc.element();
    enc.begin_sequence("User-field");
    enc.member("label", 0);
    enc.variant("Object-id", "str", 1);
    enc.visible_str(user.label);
    enc.end_variant();
    enc.end_member();
    if let Some(num) = user.num {
        enc.member("num", 1);
        enc.int(num);
        enc.end_member();
    }
    enc.member("data", 2);
    match user.data {
        UserData::Str(text) => {
            enc.variant("User-field_data", "str", 0);
            enc.visible_str(text);
        }
        UserData::Oss(bytes) => {
            enc.variant("User-field_data", "oss", 9);
            enc.begin_list(false);
            enc.element();
            enc.octets(bytes);
            enc.end_list();
        }
    }
    enc.end_variant();
    enc.end_member();
    enc.end_sequence();
    enc.end_list();
    enc.end_member();
    enc.end_sequence();
}

fn write_org_ref<E: Encoder>(enc: &mut E, source: &OrgSource<'_>) {
    enc.begin_sequence("Org-ref");
    if let Some(taxname) = source.taxname {
        enc.member("taxname", 0);
        enc.visible_str(taxname);
        enc.end_member();
    }
    if let Some(common) = source.common {
        enc.member("common", 1);
        enc.visible_str(common);
        enc.end_member();
    }
    enc.member("db", 3);
    enc.begin_list(true);
    enc.element();
    enc.begin_sequence("Dbtag");
    enc.member("db", 0);
    enc.visible_str(b"taxon");
    enc.end_member();
    enc.member("tag", 1);
    enc.variant("Object-id", "id", 0);
    enc.int(source.taxid);
    enc.end_variant();
    enc.end_member();
    enc.end_sequence();
    enc.end_list();
    enc.end_member();
    enc.end_sequence();
}

fn write_seq_inst<E: Encoder>(enc: &mut E, inst: &SeqInst<'_>) {
    enc.begin_sequence("Seq-inst");
    enc.member("repr", 0);
    match inst.data {
        Some(_) => enc.enumerated("raw", 2),
        None => enc.enumerated("virtual", 1),
    }
    enc.end_member();
    enc.member("mol", 1);
    if inst.protein {
        enc.enumerated("aa", 3);
    } else {
        enc.enumerated("na", 4);
    }
    enc.end_member();
    enc.member("length", 2);
    // `TSeqPos`, an unsigned 32-bit length.
    match i32::try_from(inst.length) {
        Ok(length) => enc.int(length),
        Err(_) => enc.big_int(i64::from(inst.length)),
    }
    enc.end_member();
    if let Some(data) = inst.data {
        enc.member("seq-data", 6);
        match data {
            SeqData::Ncbi2na(bytes) => {
                enc.variant("Seq-data", "ncbi2na", 2);
                enc.begin_named("NCBI2na");
                enc.octets(bytes);
            }
            SeqData::Ncbi4na(bytes) => {
                enc.variant("Seq-data", "ncbi4na", 3);
                enc.begin_named("NCBI4na");
                enc.octets(bytes);
            }
            SeqData::Ncbieaa(letters) => {
                enc.variant("Seq-data", "ncbieaa", 7);
                enc.begin_named("NCBIeaa");
                enc.string_store(letters);
            }
            SeqData::Ncbistdaa(bytes) => {
                enc.variant("Seq-data", "ncbistdaa", 9);
                enc.begin_named("NCBIstdaa");
                enc.octets(bytes);
            }
        }
        enc.end_named();
        enc.end_variant();
        enc.end_member();
    }
    enc.end_sequence();
}
