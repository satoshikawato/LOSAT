//! Header record of one sequence (`.phr` / `.nhr`): the BER encoding of a
//! `Blast-def-line-set`. Without `-parse_seqids` it has one definition line
//! with the title, the id `gnl|BL_ORD_ID|n` and the taxid 0 (`defline_set`; the
//! diagnostics of its title alone: `bad_title_chars`);
//! with it, the definition lines of the sequence with every Seq-id alternative
//! (`encode_defline_set`). `crate::seqdb::header` and `crate::seqid::ber`
//! decode the same bytes.
//!
//! ASN.1 definitions (NCBI 598d8ae6):
//! `c++/src/objects/blastdb/blastdb.asn:32-47` (`Blast-def-line-set`,
//! `Blast-def-line`), `c++/src/objects/seqloc/seqloc.asn:24-45` (`Seq-id`),
//! `c++/src/objects/general/general.asn:71-80` (`Dbtag`, `Object-id`).

// NCBI reference (598d8ae6): c++/include/serial/impl/objstrasnb.hpp:71-130
const TAG_INTEGER: u8 = 2;
const TAG_VISIBLE_STRING: u8 = 26;
const TAG_SEQUENCE: u8 = (1 << 5) | 16;
const TAG_CONTEXT_CONSTRUCTED: u8 = (2 << 6) | (1 << 5);
const INDEFINITE_LENGTH_BYTE: u8 = 0x80;

// The BER lengths and integer values are the ASN.1 binary encoder's, their one writer.
use crate::report::serial::asn_ber::{write_int4_value, write_int8_value, write_length};
use crate::seqid::model::{Date, ObjectId, ParsedSeqId, PatentNumber, TextseqId};
use crate::seqid::DeflineSet;
use crate::utils::byte_search::find_first;

/// Member tags of `Blast-def-line` (blastdb.asn:34-46, in order).
const DEFLINE_TITLE: u8 = 0;
const DEFLINE_SEQID: u8 = 1;
const DEFLINE_TAXID: u8 = 2;
/// `general` is the eleventh alternative of `Seq-id` (seqloc.asn:25-35).
const SEQID_GENERAL: u8 = 10;
/// Members of `Dbtag` and alternatives of `Object-id` (general.asn:71-80).
const DBTAG_DB: u8 = 0;
const DBTAG_TAG: u8 = 1;
const OBJECT_ID_ID: u8 = 0;

/// A byte of a title that a VisibleString cannot hold, as NCBI's serializer
/// reports it ("Bad char \[0xC3\] in string at byte 21").
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BadTitleChar {
    /// The byte of the title.
    pub byte: u8,
    /// Its position in the header record (`CObjectStack::GetPosition`: the
    /// stream offset at which the replacement is written).
    pub stream_offset: usize,
}

fn open(out: &mut Vec<u8>, tag: u8) {
    out.push(tag);
    out.push(INDEFINITE_LENGTH_BYTE);
}

fn close(out: &mut Vec<u8>) {
    out.extend_from_slice(&[0, 0]);
}

/// `WriteInt4`: the universal INTEGER tag and the value (`asn_ber::write_int4_value`).
fn write_int4(out: &mut Vec<u8>, data: i32) {
    out.push(TAG_INTEGER);
    write_int4_value(out, data);
}

// NCBI reference (598d8ae6): c++/include/serial/impl/objistr.inl:455-459
// NCBI reference (598d8ae6): c++/src/serial/objostrasnb.cpp:842-871
// makeblastdb runs with the default method (replace and warn) and the
// substitute '#': every byte outside ' '..'~' is written as '#' and reported
// with the stream position of the replacement. The runs of visible bytes between them are
// found a block at a time (`find_first`) and copied whole.
fn write_visible_string(out: &mut Vec<u8>, text: &[u8], bad: &mut Vec<BadTitleChar>) {
    out.push(TAG_VISIBLE_STRING);
    write_length(out, text.len());
    let mut rest = text;
    while let Some(at) = find_first(rest, is_bad_visible) {
        out.extend_from_slice(&rest[..at]);
        bad.push(BadTitleChar {
            byte: rest[at],
            stream_offset: out.len(),
        });
        out.push(b'#');
        rest = &rest[at + 1..];
    }
    out.extend_from_slice(rest);
}

/// Whether `write_visible_string` writes `byte` as `#` (a byte outside ' '..'~').
fn is_bad_visible(byte: u8) -> bool {
    !(b' '..=b'~').contains(&byte)
}

/// The bytes of `text` that `write_visible_string` writes as `#` and the stream offsets that it
/// reports, without writing them: `start` is the stream offset of the string's first byte, after
/// its tag and length.
fn bad_visible_chars(text: &[u8], start: usize, bad: &mut Vec<BadTitleChar>) {
    let mut at = 0;
    while let Some(found) = find_first(&text[at..], is_bad_visible) {
        bad.push(BadTitleChar {
            byte: text[at + found],
            stream_offset: start + at + found,
        });
        at += found + 1;
    }
}

/// The bytes of `defline_set`'s record besides its title: at most 77 (eleven tags with indefinite
/// lengths and their eleven ends, the title's tag and a length of up to 9 bytes, `BL_ORD_ID` with
/// its tag and length, and two integers of up to 6 bytes), so that the record is written into its
/// buffer without growing it (LOSAT's speed).
const DEFLINE_SET_OVERHEAD: usize = 80;

/// The prefix of `defline_set`'s record before its title's bytes: the tags of
/// `Blast-def-line-set`, `Blast-def-line` and the `title` member, and the VisibleString's tag and
/// length.
fn title_prefix(out: &mut Vec<u8>, title: &[u8]) {
    open(out, TAG_SEQUENCE);
    open(out, TAG_SEQUENCE);
    open(out, TAG_CONTEXT_CONSTRUCTED | DEFLINE_TITLE);
    out.push(TAG_VISIBLE_STRING);
    write_length(out, title.len());
}

/// The bytes of `title` that `defline_set` writes as `#`, with the stream offsets that it reports,
/// without writing the record. The title is the record's first string, so its offsets do not
/// depend on the oid or the taxid. `ExtractBioseqDeflines` serializes the definition line only for
/// these diagnostics (cited at `WriteDb::extract_bioseq_deflines`); computing them alone gives the
/// same diagnostics (LOSAT's speed).
pub(crate) fn bad_title_chars(title: &[u8], bad: &mut Vec<BadTitleChar>) {
    let mut prefix = Vec::with_capacity(16);
    title_prefix(&mut prefix, title);
    bad_visible_chars(title, prefix.len(), bad);
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:866-893
// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/taxid_set.cpp:130-136
// Without `-taxid` the best taxid is 0 (taxid_set.cpp:84-90, 123-128: the unset global id
// is returned as it is unless the definition line has one); with it, the global taxid
// (`TaxIdSet::best_taxid`). `-taxid_map` requires `-parse_seqids`.
/// The header record of the sequence number `oid` of its volume: the
/// definition line with `title`, `gnl|BL_ORD_ID|oid` and `taxid`. Bytes of
/// the title that a VisibleString cannot hold are written as `#` and
/// appended to `bad`.
pub fn defline_set(title: &[u8], oid: u32, taxid: i32, bad: &mut Vec<BadTitleChar>) -> Vec<u8> {
    let mut out = Vec::with_capacity(DEFLINE_SET_OVERHEAD + title.len());
    // Outer `Blast-def-line-set` sequence of `Blast-def-line` entries.
    open(&mut out, TAG_SEQUENCE);
    // Single `Blast-def-line` sequence.
    open(&mut out, TAG_SEQUENCE);
    // Optional `title` field (VisibleString).
    open(&mut out, TAG_CONTEXT_CONSTRUCTED | DEFLINE_TITLE);
    write_visible_string(&mut out, title, bad);
    close(&mut out);
    // `seqid` list of `Seq-id` entries.
    open(&mut out, TAG_CONTEXT_CONSTRUCTED | DEFLINE_SEQID);
    open(&mut out, TAG_SEQUENCE);
    // `Seq-id` choice variant `general` (`Dbtag`).
    open(&mut out, TAG_CONTEXT_CONSTRUCTED | SEQID_GENERAL);
    // `Dbtag` sequence containing `db` (VisibleString) and `tag` (`Object-id`).
    open(&mut out, TAG_SEQUENCE);
    open(&mut out, TAG_CONTEXT_CONSTRUCTED | DBTAG_DB);
    write_visible_string(&mut out, b"BL_ORD_ID", &mut Vec::new());
    close(&mut out);
    open(&mut out, TAG_CONTEXT_CONSTRUCTED | DBTAG_TAG);
    // `Object-id` choice variant `id` (integer; alternative is `str`).
    open(&mut out, TAG_CONTEXT_CONSTRUCTED | OBJECT_ID_ID);
    write_int4(&mut out, oid as i32);
    close(&mut out);
    close(&mut out);
    close(&mut out);
    close(&mut out);
    close(&mut out);
    close(&mut out);
    // Optional `taxid` integer field.
    open(&mut out, TAG_CONTEXT_CONSTRUCTED | DEFLINE_TAXID);
    write_int4(&mut out, taxid);
    close(&mut out);
    close(&mut out);
    close(&mut out);
    out
}

/// A byte of a VisibleString of a parsed definition line that the serializer replaces
/// (`ReplaceVisibleChar` posts "Bad char \[0xE9\] in string at byte 38", the member's
/// stack path before it and the whole string after it).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BadChar {
    /// `CObjectStack::GetStackTrace()`: the path of the member, as
    /// `Blast-def-line-set.E.seqid.E.local.str`.
    pub path: String,
    pub byte: u8,
    /// The stream offset of the replacement (`GetPosition`).
    pub stream_offset: usize,
    /// The string that holds it.
    pub string: Vec<u8>,
}

// `WriteInt8` writes the universal INTEGER tag first (objostrasnb.cpp:716-724; no
// `m_CStyleBigInt` in makeblastdb's stream), then the value (`asn_ber::write_int8_value`).
fn write_int8(out: &mut Vec<u8>, data: i64) {
    out.push(TAG_INTEGER);
    write_int8_value(out, data);
}

/// The capacity that `encode_defline_set` gives the record of `set`: 64 bytes, and for each
/// definition line its title and 64 bytes an id, so that a record of the usual size is written
/// without growing its buffer (LOSAT's speed); a larger one grows as before.
fn encoded_capacity(set: &DeflineSet) -> usize {
    64 + set
        .deflines
        .iter()
        .map(|defline| defline.title.as_ref().map_or(0, Vec::len) + 64 * defline.ids.len())
        .sum::<usize>()
}

/// The BER writer of a whole `Blast-def-line-set` (`*deflines` streamed to `oss` through
/// `MSerial_AsnBinary`), with the member path of each string for the serializer's diagnostics.
struct SetWriter<'a> {
    out: Vec<u8>,
    /// Whether the bytes of the strings are written (the header) or only counted in `skipped`
    /// (the diagnostics alone, `defline_set_bad_chars`): every other byte is written either
    /// way, so the stream offsets of the reports are the same.
    write_strings: bool,
    skipped: usize,
    path: Vec<&'static str>,
    bad: &'a mut Vec<BadChar>,
}

impl SetWriter<'_> {
    fn open(&mut self, tag: u8) {
        open(&mut self.out, tag);
    }

    fn close(&mut self) {
        close(&mut self.out);
    }

    /// A member or an alternative `[number]` around `body`, named `name` in the path.
    fn context(&mut self, number: u8, name: &'static str, body: impl FnOnce(&mut Self)) {
        self.path.push(name);
        self.open(TAG_CONTEXT_CONSTRUCTED | number);
        body(self);
        self.close();
        self.path.pop();
    }

    fn string(&mut self, text: &[u8]) {
        let mut bad = Vec::new();
        if self.write_strings {
            write_visible_string(&mut self.out, text, &mut bad);
        } else {
            self.out.push(TAG_VISIBLE_STRING);
            write_length(&mut self.out, text.len());
            bad_visible_chars(text, self.skipped + self.out.len(), &mut bad);
            self.skipped += text.len();
        }
        for bad_char in bad {
            self.bad.push(BadChar {
                path: self.path.join("."),
                byte: bad_char.byte,
                stream_offset: bad_char.stream_offset,
                string: text.to_vec(),
            });
        }
    }

    fn string_member(&mut self, number: u8, name: &'static str, text: &[u8]) {
        self.context(number, name, |w| w.string(text));
    }

    fn int_member(&mut self, number: u8, name: &'static str, value: i32) {
        self.context(number, name, |w| write_int4(&mut w.out, value));
    }

    fn object_id(&mut self, tag: &ObjectId) {
        match tag {
            ObjectId::Id(id) => self.int_member(0, "id", *id),
            ObjectId::Str(text) => self.string_member(1, "str", text),
        }
    }

    fn textseq(&mut self, id: &TextseqId) {
        self.open(TAG_SEQUENCE);
        if let Some(name) = &id.name {
            self.string_member(0, "name", name);
        }
        if let Some(accession) = &id.accession {
            self.string_member(1, "accession", accession);
        }
        if let Some(release) = &id.release {
            self.string_member(2, "release", release);
        }
        if let Some(version) = id.version {
            self.int_member(3, "version", version);
        }
        self.close();
    }

    fn date(&mut self, date: &Date) {
        match date {
            Date::Str(text) => self.string_member(0, "str", text),
            Date::Std(std) => self.context(1, "std", |w| {
                w.open(TAG_SEQUENCE);
                w.int_member(0, "year", std.year);
                if let Some(month) = std.month {
                    w.int_member(1, "month", month);
                }
                if let Some(day) = std.day {
                    w.int_member(2, "day", day);
                }
                if let Some(season) = &std.season {
                    w.string_member(3, "season", season);
                }
                if let Some(hour) = std.hour {
                    w.int_member(4, "hour", hour);
                }
                if let Some(minute) = std.minute {
                    w.int_member(5, "minute", minute);
                }
                if let Some(second) = std.second {
                    w.int_member(6, "second", second);
                }
                w.close();
            }),
        }
    }

    // The Seq-id CHOICE: each alternative is `[n]` around its value (seqloc.asn:24-45,
    // `seqid::model`), the members of a SEQUENCE in ASN.1 order, only those set.
    fn seq_id(&mut self, id: &ParsedSeqId) {
        let number = id.alternative();
        let name = id.alternative_name();
        self.context(number, name, |w| match id {
            ParsedSeqId::Local(tag) => w.object_id(tag),
            ParsedSeqId::Gibbsq(value) | ParsedSeqId::Gibbmt(value) => {
                write_int4(&mut w.out, *value)
            }
            ParsedSeqId::Giim(giim) => {
                w.open(TAG_SEQUENCE);
                w.int_member(0, "id", giim.id);
                if let Some(db) = &giim.db {
                    w.string_member(1, "db", db);
                }
                if let Some(release) = &giim.release {
                    w.string_member(2, "release", release);
                }
                w.close();
            }
            ParsedSeqId::Patent(patent) => {
                w.open(TAG_SEQUENCE);
                w.int_member(0, "seqid", patent.seqid);
                w.context(1, "cit", |w| {
                    w.open(TAG_SEQUENCE);
                    w.string_member(0, "country", &patent.cit.country);
                    w.context(1, "id", |w| match &patent.cit.id {
                        PatentNumber::Number(number) => w.string_member(0, "number", number),
                        PatentNumber::AppNumber(number) => w.string_member(1, "app-number", number),
                    });
                    if let Some(doc_type) = &patent.cit.doc_type {
                        w.string_member(2, "doc-type", doc_type);
                    }
                    w.close();
                });
                w.close();
            }
            ParsedSeqId::General(dbtag) => {
                w.open(TAG_SEQUENCE);
                w.string_member(DBTAG_DB, "db", &dbtag.db);
                w.context(DBTAG_TAG, "tag", |w| w.object_id(&dbtag.tag));
                w.close();
            }
            ParsedSeqId::Gi(gi) => write_int8(&mut w.out, *gi),
            ParsedSeqId::Pdb(pdb) => {
                w.open(TAG_SEQUENCE);
                w.string_member(0, "mol", &pdb.mol);
                if let Some(chain) = pdb.chain {
                    w.int_member(1, "chain", chain);
                }
                if let Some(rel) = &pdb.rel {
                    w.context(2, "rel", |w| w.date(rel));
                }
                if let Some(chain_id) = &pdb.chain_id {
                    w.string_member(3, "chain-id", chain_id);
                }
                w.close();
            }
            _ => w.textseq(
                id.textseq()
                    .expect("the other alternatives are Textseq-ids"),
            ),
        });
    }

    fn int_list(&mut self, number: u8, name: &'static str, values: &[i32]) {
        self.context(number, name, |w| {
            w.open(TAG_SEQUENCE);
            for &value in values {
                write_int4(&mut w.out, value);
            }
            w.close();
        });
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb_impl.cpp:887-893
// The writing of blastdb.asn:32-47 (cited in `seqid::ber`): SEQUENCE OF Blast-def-line,
// each a SEQUENCE of the members that are set, in order; `title` is a VisibleString,
// `taxid` an INTEGER, the lists SEQUENCE OF INTEGER. `seqid::ber::decode_full_set` reads it
// back.
/// The header record of a sequence whose definition lines carry parsed Seq-ids
/// (`-parse_seqids`): `set` as makeblastdb serializes it. Bytes of the strings that a
/// VisibleString cannot hold are written as `#` and appended to `bad`.
pub fn encode_defline_set(set: &DeflineSet, bad: &mut Vec<BadChar>) -> Vec<u8> {
    let mut w = SetWriter {
        out: Vec::with_capacity(encoded_capacity(set)),
        write_strings: true,
        skipped: 0,
        path: vec!["Blast-def-line-set"],
        bad,
    };
    w.set(set);
    w.out
}

/// The diagnostics of `encode_defline_set(set, bad)` without its record: the strings' bytes
/// are searched, not written (`SetWriter::write_strings`). `ExtractBioseqDeflines` serializes the
/// definition lines it built only for these diagnostics and throws the bytes away (cited at
/// `WriteDb::post_extracted_deflines`); computing them alone gives the same diagnostics (LOSAT's
/// speed: the record is written once, for the header).
pub fn defline_set_bad_chars(set: &DeflineSet, bad: &mut Vec<BadChar>) {
    let mut w = SetWriter {
        out: Vec::with_capacity(encoded_capacity(set)),
        write_strings: false,
        skipped: 0,
        path: vec!["Blast-def-line-set"],
        bad,
    };
    w.set(set);
}

impl SetWriter<'_> {
    /// The record of `set`: SEQUENCE OF Blast-def-line.
    fn set(&mut self, set: &DeflineSet) {
        self.open(TAG_SEQUENCE);
        for defline in &set.deflines {
            self.path.push("E");
            self.open(TAG_SEQUENCE);
            if let Some(title) = &defline.title {
                self.string_member(DEFLINE_TITLE, "title", title);
            }
            self.context(DEFLINE_SEQID, "seqid", |w| {
                w.open(TAG_SEQUENCE);
                for id in &defline.ids {
                    w.path.push("E");
                    w.seq_id(id);
                    w.path.pop();
                }
                w.close();
            });
            if let Some(taxid) = defline.taxid {
                self.int_member(DEFLINE_TAXID, "taxid", taxid);
            }
            if let Some(values) = &defline.memberships {
                self.int_list(3, "memberships", values);
            }
            if let Some(values) = &defline.links {
                self.int_list(4, "links", values);
            }
            if let Some(values) = &defline.other_info {
                self.int_list(5, "other-info", values);
            }
            self.close();
            self.path.pop();
        }
        self.close();
    }
}

#[cfg(test)]
#[path = "header_tests.rs"]
mod tests;
