//! `-outfmt` with `%` specifiers other than `%f`: NCBI's `CBlastDB_SeqFormatter` (one line
//! per definition line, or one per sequence when no specifier reads a definition line) and
//! the fields `CBlastDeflineUtil::ExtractDataFromBlastDefline` extracts.

use std::collections::BTreeSet;
use std::io::Write;

use crate::seqid::defline_set::{compare, Sic};
use crate::seqid::label::{bare, fasta};
use crate::seqid::rank::{find_best, text_score};
use crate::seqid::{BlastDefline, DeflineSet, ParsedSeqId};

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
    /// The bits `e_seq` and `e_hash` of `m_OtherFields`.
    need_seq: bool,
    need_hash: bool,
}

impl FormatSpec {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:59-91
    /// Reads the format `spec` (the `-outfmt` value after blastdbcmd's handling of `%f`,
    /// `%b`, `%d` and `%m`, whose algorithms are `mask_algos`); `Err` is the message of NCBI's
    /// exception. A `%` that ends the text reads the string's terminating NUL as its
    /// specifier, as NCBI's loop does.
    pub fn parse(spec: &str, mask_algos: Vec<i32>) -> Result<Self, String> {
        let bytes = spec.as_bytes();
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
            types,
            mask_algos,
            separators,
            get_defline,
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

    /// Whether a specifier asks for the ASN.1 text of a definition line (`%d`), which LOSAT
    /// does not write (line BDC DECISIONS 9).
    pub fn asks_asn_text(&self) -> bool {
        self.types.contains(&b'd')
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
        target: Option<&str>,
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
        let mut line = Vec::new();
        let mut print = |defline: Option<&BlastDefline>, line: &mut Vec<u8>| -> Result<(), Stop> {
            line.clear();
            let mut texts = Vec::new();
            let printed = self.x_print(
                source, oid, length, defline, &seq, &hash, &masks, line, &mut texts,
            );
            for text in &texts {
                out.post(text)?;
            }
            // `out` writes to an `ncbi_app::files::Ostream` (`retrieve::process`): a failed write
            // does nothing, as in NCBI's `ostream`; a closed pipe ends the retrieval.
            out.write_all(line).map_err(Stop::write)?;
            printed
        };
        if !self.get_defline {
            return print(None, &mut line);
        }
        // `GetHdr` gives a set, empty where the alias files' `OID_MASK_TYPE` erased every line
        // (line RG D1b-2): nothing is written for it, and a target is not found in it.
        let set = source.deflines(oid)?;
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
    pub fn write_missing_oid(&self, out: &mut dyn Write) -> Result<(), Stop> {
        let missing = || Stop::Ncbi(super::retrieve::OID_NOT_FOUND.to_string());
        if self.need_seq || self.need_hash || self.get_defline || !self.mask_algos.is_empty() {
            return Err(missing());
        }
        let mut line = Vec::new();
        for (i, &kind) in self.types.iter().enumerate() {
            line.extend_from_slice(&self.separators[i]);
            match kind {
                b'o' => line.extend_from_slice(b"-1"),
                b'l' => {
                    out.write_all(&line).map_err(Stop::write)?;
                    return Err(missing());
                }
                b'd' => {
                    return Err(Stop::Refused(format!(
                        "the -outfmt specifier %d is not supported by LOSAT's {}",
                        super::PROGRAM
                    )));
                }
                other => {
                    out.write_all(&line).map_err(Stop::write)?;
                    let mut message = b"Unrecognized format specification: '%".to_vec();
                    message.push(other);
                    message.push(b'\'');
                    return Err(Stop::NcbiBytes(message));
                }
            }
        }
        line.extend_from_slice(&self.separators[self.types.len()]);
        line.push(b'\n');
        out.write_all(&line).map_err(Stop::write)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:165-262
    /// One line into `line`, the texts of its taxonomy lookups into `texts`; `Err` after the
    /// part NCBI writes before an unknown specifier.
    fn x_print(
        &self,
        source: &Source,
        oid: u32,
        length: u32,
        defline: Option<&BlastDefline>,
        seq: &[u8],
        hash: &str,
        masks: &[String],
        line: &mut Vec<u8>,
        texts: &mut Vec<String>,
    ) -> Result<(), Stop> {
        let fields = match defline {
            Some(defline) => Some(DeflineFields::extract(source, defline, &self.types, texts)?),
            None => None,
        };
        let empty = DeflineFields::default();
        let fields = fields.as_ref().unwrap_or(&empty);
        let mut masks = masks.iter();
        for (i, &kind) in self.types.iter().enumerate() {
            line.extend_from_slice(&self.separators[i]);
            match kind {
                b's' => line.extend_from_slice(seq),
                b'a' => line.extend_from_slice(&fields.accession),
                b'i' => line.extend_from_slice(&fields.seq_id),
                b'g' => line.extend_from_slice(fields.gi.as_bytes()),
                b'o' => line.extend_from_slice(oid.to_string().as_bytes()),
                b't' => line.extend_from_slice(&fields.title),
                b'h' => line.extend_from_slice(hash.as_bytes()),
                b'l' => line.extend_from_slice(length.to_string().as_bytes()),
                b'T' => line.extend_from_slice(fields.tax_id.as_bytes()),
                b'X' => line.extend_from_slice(fields.leaf_node_tax_ids.as_bytes()),
                b'P' => line.extend_from_slice(fields.pig.as_bytes()),
                b'L' => line.extend_from_slice(&fields.common_name),
                b'C' => line.extend_from_slice(&fields.leaf_node_common_names),
                b'B' => line.extend_from_slice(&fields.blast_name),
                b'K' => line.extend_from_slice(&fields.super_kingdom),
                b'S' => line.extend_from_slice(&fields.scientific_name),
                b'N' => line.extend_from_slice(&fields.leaf_node_scientific_names),
                b'e' => line.extend_from_slice(fields.membership.as_bytes()),
                b'n' => line.extend_from_slice(fields.links.as_bytes()),
                b'm' => {
                    if let Some(mask) = masks.next() {
                        line.extend_from_slice(mask.as_bytes());
                    }
                }
                // `%d` is refused before the first line.
                b'd' => {
                    return Err(Stop::Refused(format!(
                        "the -outfmt specifier %d is not supported by LOSAT's {}",
                        super::PROGRAM
                    )));
                }
                other => {
                    let mut message = b"Unrecognized format specification: '%".to_vec();
                    message.push(other);
                    message.push(b'\'');
                    return Err(Stop::NcbiBytes(message));
                }
            }
        }
        line.extend_from_slice(&self.separators[self.types.len()]);
        line.push(b'\n');
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

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:800-819
/// The definition line of `-target_only` for the `%` formats: the first whose Seq-id
/// matches the target (`CSeq_id::Match`) or whose GI is the target's number.
fn extract_target_defline<'a>(set: &'a DeflineSet, target: &str) -> Result<&'a BlastDefline, Stop> {
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
    Err(Stop::Ncbi(format!("Failed to find target id {target}")))
}

/// The strings of `CBlastDeflineUtil::ExtractDataFromBlastDefline` (`results`), those the
/// format asks for; the others empty.
#[derive(Default)]
struct DeflineFields {
    gi: String,
    seq_id: Vec<u8>,
    accession: Vec<u8>,
    title: Vec<u8>,
    tax_id: String,
    scientific_name: Vec<u8>,
    common_name: Vec<u8>,
    blast_name: Vec<u8>,
    super_kingdom: Vec<u8>,
    leaf_node_tax_ids: String,
    leaf_node_scientific_names: Vec<u8>,
    leaf_node_common_names: Vec<u8>,
    membership: String,
    pig: String,
    links: String,
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:821-827
fn check_name(name: &[u8]) -> Vec<u8> {
    if name == b"-" || name == b"unclassified" {
        NOT_AVAILABLE.as_bytes().to_vec()
    } else {
        name.to_vec()
    }
}

impl DeflineFields {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:829-957
    // NCBI reference (598d8ae6): c++/src/objects/blastdb/Blast_def_line.cpp:143-163
    /// The fields of `defline` that the specifiers `types` read, the texts of their taxonomy
    /// lookups appended to `texts`. The sum of the memberships is taken in 64 bits (NCBI's
    /// `int` would overflow: LOSAT does not reproduce that).
    fn extract(
        source: &Source,
        defline: &BlastDefline,
        types: &[u8],
        texts: &mut Vec<String>,
    ) -> Result<Self, Stop> {
        let wants = |kinds: &[u8]| types.iter().any(|t| kinds.contains(t));
        let mut fields = DeflineFields::default();
        if wants(b"g") {
            fields.gi = match defline.ids.iter().find_map(ParsedSeqId::gi) {
                Some(gi) => gi.to_string(),
                None => NOT_AVAILABLE.to_string(),
            };
        }
        if wants(b"ai") {
            let Some(best) = find_best(&defline.ids, text_score) else {
                // NCBI dereferences the null choice (it crashes on a definition line
                // without Seq-ids).
                return Err(Stop::Refused(format!(
                    "a definition line without Seq-ids is not supported by LOSAT's {}",
                    super::PROGRAM
                )));
            };
            fasta(best, &mut fields.seq_id);
            bare(best, &mut fields.accession);
        }
        if wants(b"t") {
            fields.title = match &defline.title {
                Some(title) => title.clone(),
                None => NOT_AVAILABLE.as_bytes().to_vec(),
            };
        }
        if wants(b"TLBKS") {
            let tax_id = defline.taxid.unwrap_or(0);
            fields.tax_id = tax_id.to_string();
            if wants(b"LBKS") {
                let (found, lookup_texts) = source.tax_info(tax_id)?;
                texts.extend(lookup_texts);
                match found.names() {
                    Some(names) => {
                        // A kingdom that NCBI reads past its mapping of taxdb.btd: refused
                        // where `%K` writes it.
                        fields.super_kingdom = match &names.super_kingdom {
                            Ok(kingdom) => check_name(kingdom),
                            Err(past) if wants(b"K") => {
                                return Err(past.error("super kingdom").into())
                            }
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
                fields.leaf_node_tax_ids.push_str(separator);
                fields.leaf_node_tax_ids.push_str(&tax_id.to_string());
                separator = SEPARATOR;
            }
        }
        if wants(b"e") {
            let membership: i64 = defline
                .memberships
                .iter()
                .flatten()
                .map(|&m| i64::from(m))
                .sum();
            fields.membership = membership.to_string();
        }
        if wants(b"P") {
            let pig = defline
                .other_info
                .iter()
                .flatten()
                .copied()
                .find(|&value| value != -1)
                .unwrap_or(-1);
            fields.pig = pig.to_string();
        }
        if wants(b"n") {
            fields.links = match &defline.links {
                Some(links) => links
                    .iter()
                    .map(|link| format!("{link}{SEPARATOR}"))
                    .collect(),
                None => NOT_AVAILABLE.to_string(),
            };
        }
        Ok(fields)
    }
}

#[cfg(test)]
#[path = "format_tests.rs"]
mod tests;
