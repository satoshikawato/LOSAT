//! NCBI Seq-ids of database sequences: the one place where LOSAT models, decodes, chooses
//! and prints the ids that a BLAST database built with `-parse_seqids` stores (NCBI's
//! `CSeq_id`, `CBlast_def_line_set` and the formatters' id helpers).
//!
//! - `model`: the Seq-id value of every alternative (`ParsedSeqId`) and its serial writer
//!   (outfmt 8, 9, 11, 12, through `report::seqid::SeqId::Db`).
//! - `ber`: the BER decoding of a header record (`decode_defline_set`), whose reader
//!   `seqdb::header` shares for the ordinal ids of a database without parsed Seq-ids.
//! - `label`: the printings (`fasta`, `content`, `format_label`, `bare`,
//!   `seqid_list_string`, `first_gi`).
//! - `rank`: `text_score` (Score, WorstRank), `best_rank_score` and `find_best`.
//! - `defline_set`: the definition lines of one sequence, the Seq-align id, the definition
//!   line of the subject Bioseq and of the shown id, and the joins of the formats that show
//!   every definition line.
//!
//! A format chooses and prints a database sequence's ids only through these functions; it
//! keeps only its own joins and escaping. The search programs refuse a database with parsed
//! Seq-ids when they open it (`algorithm::common::database::refuse_parsed_seqids`); the
//! programs that read one take their ids from `seqdb::SeqDb::defline_set`.

pub mod ber;
pub mod defline_set;
pub mod label;
pub mod model;
pub mod rank;

#[cfg(test)]
mod tests;

pub use defline_set::{BlastDefline, DeflineSet};
pub use model::ParsedSeqId;
