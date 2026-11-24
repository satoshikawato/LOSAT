//! `-outfmt` with `%` specifiers other than `%f`: NCBI's `CBlastDB_SeqFormatter` (one line
//! per definition line, or one per sequence when no specifier reads a definition line) and
//! the fields `CBlastDeflineUtil::ExtractDataFromBlastDefline` extracts.

use std::collections::BTreeSet;

use crate::seqid::defline_set::{compare, Sic};
use crate::seqid::label::{bare, fasta};
use crate::seqid::rank::{find_best, text_score};
use crate::seqid::{BlastDefline, DeflineSet, ParsedSeqId};

use crate::report::serial::blastdb::{write_defline_document, write_defline_set_document};
use crate::report::serial::encoder::Encoder;

use crate::report::ostream::CutBuffer;
use crate::utils::decimal::push_decimal;

use super::bioseq::AsnTextDocument;
use super::letters::{hash_text, sequence_hash};
use super::retrieve::{numeric_target, parse_target, Config, Output, Source, Stop};

/// `NOT_AVAILABLE`.
const NOT_AVAILABLE: &str = "N/A";
/// `kNoMasksFound`.
const NO_MASKS_FOUND: &str = "none";
/// `SEPARATOR`.
const SEPARATOR: &str = ";";

/// The `-outfmt` text read into its specifiers (`m_ReplTypes`) and the literal text
/// between them (`m_Seperators`, one more than the specifiers).
#[derive(Debug)]
pub(super) struct FormatSpec {
    types: Vec<u8>,
    /// The filtering algorithm of each `%m` of `types`.
    mask_algos: Vec<i32>,
    separators: Vec<Vec<u8>>,
    /// `m_GetDefline`: a specifier reads the definition lines.
    get_defline: bool,
    /// `m_DeflineFields.asn_defline`: `%d`, the ASN.1 text of the definition lines (which the
    /// request makes the whole format, `retrieve::init_search_request`).
    asn_defline: bool,
    /// Whether `m_OtherFields` has the bit of `e_seq` (`need_seq`) and the bit of `e_hash`
    /// (`need_hash`).
    need_seq: bool,
    need_hash: bool,
    /// The bytes of the line being written, kept from one line to the next.
    line: std::cell::RefCell<CutBuffer>,
}

impl FormatSpec {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:59-91
    /// Reads the format `spec` (the `-outfmt` value after blastdbcmd's handling of `%f`,
    /// `%b`, `%d` and `%m`, whose algorithms are `mask_algos`); `Err` is the message of NCBI's
    /// exception. A `%` that ends the text reads the string's terminating NUL as its
    /// specifier, as NCBI's loop does.
    pub fn parse(spec: &[u8], mask_algos: Vec<i32>) -> Result<Self, String> {
        let bytes = spec;
        let mut types = Vec::new();
        let mut separators = Vec::new();
        let mut sp = Vec::new();
        let mut i = 0;
        while i < bytes.len() {
            if bytes[i] == b'%' {
                if bytes.get(i + 1) == Some(&b'%') {
                    i += 1;
                    sp.push(b'%');
                } else {
                    i += 1;
                    types.push(bytes.get(i).copied().unwrap_or(0));
                    separators.push(std::mem::take(&mut sp));
                }
            } else {
                sp.push(bytes[i]);
            }
            i += 1;
        }
        separators.push(sp);
        if types.is_empty() {
            return Err("Invalid format specification".to_string());
        }
        let get_defline = types.iter().any(|t| b"aigtTXPLBKSCNend".contains(t));
        Ok(FormatSpec {
            need_seq: types.contains(&b's'),
            need_hash: types.contains(&b'h'),
            asn_defline: types.contains(&b'd'),
            types,
            mask_algos,
            separators,
            get_defline,
            line: std::cell::RefCell::default(),
        })
    }

    /// Whether a specifier reads the definition lines (`m_GetDefline`: `Write` reads them with
    /// `GetHdr`).
    pub fn reads_deflines(&self) -> bool {
        self.get_defline
    }

    /// Whether a specifier writes the letters of the sequence (`%s`).
    pub fn writes_letters(&self) -> bool {
        self.need_seq
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:315-358
    /// Writes the lines of sequence `oid` to `out`: one per definition line (the one of
    /// `target` alone with `-target_only`), or one line when no specifier reads them. The
    /// line NCBI has written up to an unknown specifier is written before its error. The texts
    /// of a line's taxonomy lookups (`GetTaxInfo` of `%L %B %K %S %C %N`) are posted before it.
    pub fn write(
        &self,
        source: &Source,
        oid: u32,
        config: &Config,
        target: Option<&[u8]>,
        out: &mut Output<'_, '_>,
    ) -> Result<(), Stop> {
        let mut seq = Vec::new();
        let mut hash = String::new();
        if self.need_seq {
            seq = x_get_seq(source, oid, config)?;
        }
        if self.need_hash {
            hash = hash_text(sequence_hash(&source.string_letters(oid, None)?));
        }
        let masks = self
            .mask_algos
            .iter()
            .map(|&algo| Ok(masks_string(&source.masks(oid, algo)?)))
            .collect::<Result<Vec<String>, Stop>>()?;
        let length = source.length(oid)?;
        let mut line = self.line.borrow_mut();
        let mut print =
            |defline: Option<&BlastDefline>, line: &mut CutBuffer| -> Result<(), Stop> {
                let mut texts = Vec::new();
                let printed = self.x_print(
                    source, oid, length, defline, &seq, &hash, &masks, line, &mut texts,
                );
                for text in &texts {
                    out.post(text)?;
                }
                // `out` writes to a `files::OutputFile` (`retrieve::process`): a failed write does
                // nothing, as in NCBI's `ostream`; a closed pipe ends the retrieval.
                out.insertions(line)?;
                printed
            };
        if !self.get_defline {
            return print(None, &mut line);
        }
        // A null set writes nothing; a set of no lines (stored so, or emptied by the alias files'
        // `OID_MASK_TYPE`, line RG D1b-2) writes `%d`'s empty document, no line of the others,
        // and throws for a target (probes e1-e8, z1-z9 of line RG-Dbc).
        let Some(set) = source.deflines(oid)? else {
            return Ok(());
        };
        if self.asn_defline {
            return write_asn_deflines(&set, target, out);
        }
        match target {
            None => {
                for defline in &set.deflines {
                    print(Some(defline), &mut line)?;
                }
                Ok(())
            }
            Some(target) => {
                let defline = extract_target_defline(&set, target)?;
                print(Some(defline), &mut line)
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:187-201
    // `Write(-1, config)` (`-ipg` of a PIG no volume holds): the sequence, its hash and the
    // definition lines are read before the line (`x_GetSeq`, `x_GetSeqHash`, `GetHdr`), and
    // each read throws `OID not found`; without them the line is written up to `%l`, whose
    // length throws, and `%o` is `-1`.
    /// NCBI's `Write` of the OID -1: the part of the line it writes, then `OID not found`
    /// where it reads the database.
    pub fn write_missing_oid(&self, out: &mut Output<'_, '_>) -> Result<(), Stop> {
        let missing = || Stop::Ncbi(super::retrieve::OID_NOT_FOUND.to_string());
        if self.need_seq || self.need_hash || self.get_defline || !self.mask_algos.is_empty() {
            return Err(missing());
        }
        let mut line = CutBuffer::default();
        for (i, &kind) in self.types.iter().enumerate() {
            line.insertion(&self.separators[i]);
            match kind {
                b'o' => line.insertion(b"-1"),
                b'l' => {
                    out.insertions(&mut line)?;
                    return Err(missing());
                }
                other => {
                    out.insertions(&mut line)?;
                    let mut message = b"Unrecognized format specification: '%".to_vec();
                    message.push(other);
                    message.push(b'\'');
                    return Err(Stop::NcbiBytes(message));
                }
            }
        }
        line.insertion(&self.separators[self.types.len()]);
        line.endl();
        out.insertions(&mut line)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:165-262
    /// One line into `line`, the texts of its taxonomy lookups into `texts`; `Err` after the
    /// part NCBI writes before an unknown specifier. Each separator and each field is one of
    /// NCBI's insertions; the line ends with the last separator and `endl`.
    fn x_print(
        &self,
        source: &Source,
        oid: u32,
        length: u32,
        defline: Option<&BlastDefline>,
        seq: &[u8],
        hash: &str,
        masks: &[String],
        text: &mut CutBuffer,
        texts: &mut Vec<String>,
    ) -> Result<(), Stop> {
        let fields = match defline {
            Some(defline) => DeflineFields::extract(source, defline, &self.types, texts)?,
            None => DeflineFields::default(),
        };
        let mut masks = masks.iter();
        for (i, &kind) in self.types.iter().enumerate() {
            text.insertion(&self.separators[i]);
            let line = &mut text.bytes;
            match kind {
                b's' => line.extend_from_slice(seq),
                b'o' => push_decimal(line, i64::from(oid)),
                b'h' => line.extend_from_slice(hash.as_bytes()),
                b'l' => push_decimal(line, i64::from(length)),
                b'a' | b'i' | b'g' | b't' | b'T' | b'P' | b'e' | b'n' => {
                    if let Some(defline) = defline {
                        fields.write_defline_field(kind, defline, line);
                    }
                }
                b'X' => line.extend_from_slice(&fields.leaf_node_tax_ids),
                b'L' => line.extend_from_slice(&fields.common_name),
                b'C' => line.extend_from_slice(&fields.leaf_node_common_names),
                b'B' => line.extend_from_slice(&fields.blast_name),
                b'K' => line.extend_from_slice(&fields.super_kingdom),
                b'S' => line.extend_from_slice(&fields.scientific_name),
                b'N' => line.extend_from_slice(&fields.leaf_node_scientific_names),
                b'm' => {
                    if let Some(mask) = masks.next() {
                        line.extend_from_slice(mask.as_bytes());
                    }
                }
                other => {
                    let mut message = b"Unrecognized format specification: '%".to_vec();
                    message.push(other);
                    message.push(b'\'');
                    return Err(Stop::NcbiBytes(message));
                }
            }
            text.end_insertion();
        }
        text.insertion(&self.separators[self.types.len()]);
        text.endl();
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:282-313
// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:1077-1088
/// The text of `%m`: the stored ranges of the sequence, or `none`.
fn masks_string(masks: &[(u32, u32)]) -> String {
    if masks.is_empty() {
        return NO_MASKS_FOUND.to_string();
    }
    masks
        .iter()
        .map(|(first, second)| format!("{first}-{second}{SEPARATOR}"))
        .collect()
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:282-313
// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:1048-1070
// Without a range the stored offsets index the string as they are: a range past the end, or
// one that ends before it starts, runs off the string (NCBI crashed, probes M786, M886, M896
// of line BDC). LOSAT refuses those.
/// The letters of `%s`: the range's (the whole sequence when the range is empty once
/// clipped), with the masks of `-mask_sequence_with` in lower case, reverse-complemented for
/// `-strand minus` (nucleotide databases: `Config` never holds the minus strand for a protein
/// one).
fn x_get_seq(source: &Source, oid: u32, config: &Config) -> Result<Vec<u8>, Stop> {
    let range = config.range.clipped(source.length(oid)?);
    let mut seq = source.string_letters(oid, range)?;
    if config.filt_algo >= 0 {
        let masks = source.masks(oid, config.filt_algo)?;
        let used = match range {
            Some((_, to)) => {
                let past = masks.iter().position(|&(first, _)| first > to);
                &masks[..past.unwrap_or(masks.len())]
            }
            None => {
                let outside =
                    |&(first, second): &(u32, u32)| first > second || second as usize > seq.len();
                if masks.iter().any(outside) {
                    return Err(Stop::refused(format!(
                        "-outfmt %s with a stored mask outside the sequence at OID {oid} \
                         (NCBI writes past the letters)"
                    )));
                }
                &masks[..]
            }
        };
        super::letters::apply_masks(&mut seq, used, range);
    }
    if config.strand_minus {
        super::letters::reverse_complement(&mut seq);
    }
    Ok(seq)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:331-356
// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:952-956
/// `%d`: the definition lines of `GetHdr` as an ASN.1 text document, or with `-target_only` the
/// target's line (`ExtractDataFromBlastDeflineSet`) as the one field of a line, the document then
/// the line's end.
fn write_asn_deflines(
    set: &DeflineSet,
    target: Option<&[u8]>,
    out: &mut Output<'_, '_>,
) -> Result<(), Stop> {
    match target {
        None => out.asn_text(&DeflineSetDocument(set)),
        Some(target) => {
            let defline = extract_target_defline(set, target)?;
            out.asn_text_line(&DeflineDocument(defline))
        }
    }
}

/// The document that `MSerial_AsnText` writes for the set `df_set`.
struct DeflineSetDocument<'a>(&'a DeflineSet);

impl AsnTextDocument for DeflineSetDocument<'_> {
    fn write<E: Encoder>(&self, enc: &mut E) {
        write_defline_set_document(enc, self.0);
    }
}

/// The document that `MSerial_AsnText` writes for the line `dl`.
struct DeflineDocument<'a>(&'a BlastDefline);

impl AsnTextDocument for DeflineDocument<'_> {
    fn write<E: Encoder>(&self, enc: &mut E) {
        write_defline_document(enc, self.0);
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:800-819
/// The definition line of `-target_only` for the `%` formats: the first whose Seq-id
/// matches the target (`CSeq_id::Match`) or whose GI is the target's number.
fn extract_target_defline<'a>(
    set: &'a DeflineSet,
    target: &[u8],
) -> Result<&'a BlastDefline, Stop> {
    let target_id = parse_target(target)?;
    let number = numeric_target(target);
    for defline in &set.deflines {
        for id in &defline.ids {
            let matches = compare(id, &target_id) == Sic::Yes
                || number.is_some_and(|number| id.gi() == Some(number));
            if matches {
                return Ok(defline);
            }
        }
    }
    Err(Stop::NcbiBytes(
        [b"Failed to find target id ".as_slice(), target].concat(),
    ))
}

/// What `CBlastDeflineUtil::ExtractDataFromBlastDefline` reads of a definition line
/// (`results`) for the specifiers of the format: the id of `%a` and `%i`, and the taxonomy's
/// names (whose lookups post their texts before the line and fail before it is written); the
/// others empty. The fields read from the line itself are written where their specifiers stand
/// (`write_defline_field`), with no string of their own.
#[derive(Default)]
struct DeflineFields<'d> {
    /// The id that `%a` and `%i` write (`find_best` of `text_score`).
    best: Option<&'d ParsedSeqId>,
    scientific_name: Vec<u8>,
    common_name: Vec<u8>,
    blast_name: Vec<u8>,
    super_kingdom: Vec<u8>,
    leaf_node_tax_ids: Vec<u8>,
    leaf_node_scientific_names: Vec<u8>,
    leaf_node_common_names: Vec<u8>,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:821-827
fn check_name(name: &[u8]) -> Vec<u8> {
    if name == b"-" || name == b"unclassified" {
        NOT_AVAILABLE.as_bytes().to_vec()
    } else {
        name.to_vec()
    }
}

impl<'d> DeflineFields<'d> {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:829-957
    /// The fields of `defline` that the specifiers `types` read and that need more than the line:
    /// its best id (a line without Seq-ids is refused here, before the line), the taxonomy's
    /// names, the texts of their lookups appended to `texts`.
    fn extract(
        source: &Source,
        defline: &'d BlastDefline,
        types: &[u8],
        texts: &mut Vec<String>,
    ) -> Result<Self, Stop> {
        let wants = |kinds: &[u8]| types.iter().any(|t| kinds.contains(t));
        let mut fields = DeflineFields::default();
        if wants(b"ai") {
            let Some(best) = find_best(&defline.ids, text_score) else {
                // NCBI dereferences the null choice of a definition line without Seq-ids: its
                // request posts the exception's message, after the records before it (line RG
                // D1c-1c-b).
                return Err(super::retrieve::null_pointer());
            };
            fields.best = Some(best);
        }
        if wants(b"LBKS") {
            let tax_id = defline.taxid.unwrap_or(0);
            let (found, lookup_texts) = source.tax_info(tax_id)?;
            texts.extend(lookup_texts);
            match found.names() {
                Some(names) => {
                    // A kingdom that NCBI reads past its mapping of taxdb.btd: refused
                    // where `%K` writes it.
                    fields.super_kingdom = match &names.super_kingdom {
                        Ok(kingdom) => check_name(kingdom),
                        Err(past) if wants(b"K") => return Err(past.error("super kingdom").into()),
                        Err(_) => Vec::new(),
                    };
                    fields.scientific_name = names.scientific_name.clone();
                    fields.common_name = names.common_name.clone();
                    fields.blast_name = check_name(&names.blast_name);
                }
                None => {
                    for field in [
                        &mut fields.scientific_name,
                        &mut fields.common_name,
                        &mut fields.blast_name,
                        &mut fields.super_kingdom,
                    ] {
                        *field = NOT_AVAILABLE.as_bytes().to_vec();
                    }
                }
            }
        }
        if wants(b"XCN") {
            let mut tax_id_set: BTreeSet<i32> = defline.links.iter().flatten().copied().collect();
            if tax_id_set.is_empty() {
                tax_id_set.insert(defline.taxid.unwrap_or(0));
            }
            let mut separator = "";
            for tax_id in tax_id_set {
                if wants(b"CN") {
                    fields
                        .leaf_node_scientific_names
                        .extend_from_slice(separator.as_bytes());
                    fields
                        .leaf_node_common_names
                        .extend_from_slice(separator.as_bytes());
                    let (found, lookup_texts) = source.tax_info(tax_id)?;
                    texts.extend(lookup_texts);
                    match found.names() {
                        Some(names) => {
                            fields
                                .leaf_node_scientific_names
                                .extend_from_slice(&names.scientific_name);
                            fields
                                .leaf_node_common_names
                                .extend_from_slice(&names.common_name);
                        }
                        None => {
                            fields
                                .leaf_node_scientific_names
                                .extend_from_slice(NOT_AVAILABLE.as_bytes());
                            fields
                                .leaf_node_common_names
                                .extend_from_slice(NOT_AVAILABLE.as_bytes());
                        }
                    }
                }
                fields
                    .leaf_node_tax_ids
                    .extend_from_slice(separator.as_bytes());
                push_decimal(&mut fields.leaf_node_tax_ids, i64::from(tax_id));
                separator = SEPARATOR;
            }
        }
        Ok(fields)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:829-957
    // NCBI reference (598d8ae6): c++/src/objects/blastdb/Blast_def_line.cpp:143-163
    /// Appends to `line` the field of the specifier `kind` (`%a %i %g %t %T %P %e %n`) read from
    /// `defline` itself. The sum of the memberships is taken in 64 bits (NCBI's `int` would
    /// overflow: LOSAT does not reproduce that).
    fn write_defline_field(&self, kind: u8, defline: &BlastDefline, line: &mut Vec<u8>) {
        match kind {
            b'a' => {
                if let Some(best) = self.best {
                    bare(best, line);
                }
            }
            b'i' => {
                if let Some(best) = self.best {
                    fasta(best, line);
                }
            }
            b'g' => match defline.ids.iter().find_map(ParsedSeqId::gi) {
                Some(gi) => push_decimal(line, gi),
                None => line.extend_from_slice(NOT_AVAILABLE.as_bytes()),
            },
            b't' => match &defline.title {
                Some(title) => line.extend_from_slice(title),
                None => line.extend_from_slice(NOT_AVAILABLE.as_bytes()),
            },
            b'T' => push_decimal(line, i64::from(defline.taxid.unwrap_or(0))),
            b'P' => {
                let pig = defline
                    .other_info
                    .iter()
                    .flatten()
                    .copied()
                    .find(|&value| value != -1)
                    .unwrap_or(-1);
                push_decimal(line, i64::from(pig));
            }
            b'e' => {
                let membership: i64 = defline
                    .memberships
                    .iter()
                    .flatten()
                    .map(|&m| i64::from(m))
                    .sum();
                push_decimal(line, membership);
            }
            _ => match &defline.links {
                Some(links) => {
                    for &link in links {
                        push_decimal(line, i64::from(link));
                        line.extend_from_slice(SEPARATOR.as_bytes());
                    }
                }
                None => line.extend_from_slice(NOT_AVAILABLE.as_bytes()),
            },
        }
    }
}

#[cfg(test)]
#[path = "format_tests.rs"]
mod tests;
