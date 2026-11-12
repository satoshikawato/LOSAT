//! TEMPORARY until the mask FASTA reader (feature/rg-k1-fasta) is merged: dustmasker's FASTA
//! input (`CMaskFastaReader`).
use std::io::BufRead;

use super::work::{Failure, Record, SequenceSource};

pub(super) fn source<R: BufRead + 'static>(
    input: R,
    parse_seqids: bool,
) -> Box<dyn SequenceSource> {
    let _ = (input, parse_seqids);
    Box::new(Pending)
}

struct Pending;

impl SequenceSource for Pending {
    fn next_record(&mut self) -> Result<Option<Record>, Failure> {
        Err(Failure::Refused("FASTA input is not wired yet".to_string()))
    }
}
