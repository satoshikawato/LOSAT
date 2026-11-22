//! The definition lines of a database sequence as NCBI's object streams write them
//! (NCBI-BlastDL `Blast-def-line-set` and `Blast-def-line`): blastdbcmd `-outfmt %d`, which writes
//! the set of a sequence, or with `-target_only` the line of the target, as an ASN.1 text
//! document. The model is the reader's (`seqid::DeflineSet`, every member decoded); each Seq-id is
//! written by the Seq-id model's writer (`report::seqid::SeqId::Db`).
//!
//! NCBI reference (598d8ae6): c++/src/objects/blastdb/blastdb.asn:32-47

use crate::report::seqid::SeqId;
use crate::report::serial::encoder::Encoder;
use crate::seqid::{BlastDefline, DeflineSet};

/// Writes the document of the set `set` (what `MSerial_AsnText` writes for `df_set`).
pub fn write_defline_set_document<E: Encoder>(enc: &mut E, set: &DeflineSet) {
    enc.begin_document("Blast-def-line-set");
    enc.begin_named("Blast-def-line-set");
    enc.begin_list(false);
    for defline in &set.deflines {
        enc.element();
        write_defline(enc, defline);
    }
    enc.end_list();
    enc.end_named();
    enc.end_document();
}

/// Writes the document of the definition line `defline` (what `MSerial_AsnText` writes for `dl`).
pub fn write_defline_document<E: Encoder>(enc: &mut E, defline: &BlastDefline) {
    enc.begin_document("Blast-def-line");
    write_defline(enc, defline);
    enc.end_document();
}

/// Writes the `Blast-def-line` `defline`: its members that are set, in module order.
fn write_defline<E: Encoder>(enc: &mut E, defline: &BlastDefline) {
    enc.begin_sequence("Blast-def-line");
    if let Some(title) = &defline.title {
        enc.member("title", 0);
        enc.visible_str(title);
        enc.end_member();
    }
    enc.member("seqid", 1);
    enc.begin_list(false);
    for id in &defline.ids {
        enc.element();
        SeqId::Db(id).write(enc);
    }
    enc.end_list();
    enc.end_member();
    if let Some(taxid) = defline.taxid {
        enc.member("taxid", 2);
        enc.int(taxid);
        enc.end_member();
    }
    for (tag, name, values) in [
        (3, "memberships", &defline.memberships),
        (4, "links", &defline.links),
        (5, "other-info", &defline.other_info),
    ] {
        if let Some(values) = values {
            enc.member(name, tag);
            enc.begin_list(false);
            for &value in values {
                enc.element();
                enc.int(value);
            }
            enc.end_list();
            enc.end_member();
        }
    }
    enc.end_sequence();
}
