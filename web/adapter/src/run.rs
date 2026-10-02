//! The search: argv parsing (`validate`) and `run` (docs/web/abi_v2.md §4, §7-8).

use std::fmt::Write as _;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use bio::io::fasta;
use LOSAT::api::local_blast::{
    run_local_blastn, run_local_blastp, run_local_tblastn, run_local_tblastx, FormatObserver,
    FormatOutput, HspIndex, OutputSink, ReportOutputs,
};
use LOSAT::cli::{Cli, Commands};
use LOSAT::report::PairwiseHit;

use crate::emit::{emit, StreamWriter, STREAM_DIAGNOSTICS, STREAM_HITS};
use crate::json;

/// The programs that ABI v2 runs. BLASTX joins in session SX.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Program {
    Blastp,
    Tblastn,
    Blastn,
    Tblastx,
}

impl Program {
    pub fn parse(name: &str) -> Result<Self, String> {
        match name {
            "blastp" => Ok(Self::Blastp),
            "tblastn" => Ok(Self::Tblastn),
            "blastn" => Ok(Self::Blastn),
            "tblastx" => Ok(Self::Tblastx),
            "blastx" => Err("blastx is not available in LOSAT Web ABI v2 yet".into()),
            other => Err(format!("unknown program '{other}'")),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Blastp => "blastp",
            Self::Tblastn => "tblastn",
            Self::Blastn => "blastn",
            Self::Tblastx => "tblastx",
        }
    }

    /// The output formats that the engine implements for this program; a run writes
    /// exactly these, each on the stream of the same number.
    pub fn formats(self) -> &'static [u32] {
        match self {
            Self::Blastp | Self::Tblastn | Self::Blastn | Self::Tblastx => &[0, 6, 7],
        }
    }
}

/// Parses an argv (program name first) with the engine's CLI parser.
pub fn parse(words: &[&str]) -> Result<(Program, Commands), String> {
    for word in words.iter().skip(1) {
        let key = word.split('=').next().unwrap_or("");
        if matches!(key, "-out" | "--out" | "-outfmt" | "--outfmt") {
            return Err(format!(
                "{key} is not accepted: a run writes every output format of the program"
            ));
        }
    }
    // The formats come from `Program::formats`, and every program accepts the CLI's
    // default `-outfmt 0`, so the host's words are parsed exactly as on the command line,
    // with the CLI's binary name: an error (an unknown program too) is the CLI's message.
    let argv = std::iter::once("LOSAT").chain(words.iter().copied());
    let cli: Cli =
        LOSAT::cli::try_parse_from(argv).map_err(|error| LOSAT::cli::render_message(&error))?;
    let program = Program::parse(words.first().copied().unwrap_or(""))?;
    Ok((program, cli.command))
}

/// `validate`: `parse`, then BLASTN's `-dust` value and the scoring options that NCBI
/// rejects or that its Karlin-Altschul tables do not support, with NCBI's message (for a
/// batch of one query).
/// `run` leaves them to the engine, which reports them as the CLI does.
pub fn validate(words: &[&str]) -> Result<(), String> {
    if let (_, Commands::Blastn(mut args)) = parse(words)? {
        args.resolve_dust()
            .and_then(|()| LOSAT::algorithm::blastn::scoring::check_scoring(&args))
            .map_err(|error| format!("{error:#}"))?;
    }
    Ok(())
}

/// Records the byte range of every HSP row (stream 6) and section (stream 0), and of
/// the subject heading (stream 0) that each section belongs to.
struct RangeRecorder {
    streams: Vec<u32>,
    positions: Vec<Arc<AtomicU64>>,
    open: Vec<Option<(HspIndex, u64)>>,
    out6: Vec<Option<(u64, u64)>>,
    out0: Vec<Option<(u64, u64)>>,
    out0_subject: Vec<Option<(u64, u64)>>,
    /// The start of the heading being written, and the last complete heading.
    subject_start: Option<u64>,
    subject: Option<(u64, u64)>,
    /// The first observer event that does not fit the formatter contract; the run fails.
    fault: Option<String>,
}

impl RangeRecorder {
    fn position(&self, format: usize) -> u64 {
        self.positions[format].load(Ordering::SeqCst)
    }
}

/// Stores the range of `hsp`; a second range for the same HSP is a fault.
fn store(
    ranges: &mut Vec<Option<(u64, u64)>>,
    hsp: HspIndex,
    range: (u64, u64),
) -> Result<(), String> {
    if ranges.len() <= hsp {
        ranges.resize(hsp + 1, None);
    }
    match ranges[hsp].replace(range) {
        None => Ok(()),
        Some(_) => Err(format!("HSP {hsp} was reported twice in one format")),
    }
}

impl FormatObserver for RangeRecorder {
    fn hsp_begin(&mut self, format: usize, hsp: HspIndex) {
        let start = self.position(format);
        if let Some((open, _)) = self.open[format].replace((hsp, start)) {
            self.fault
                .get_or_insert(format!("HSP {hsp} started before HSP {open} ended"));
        }
    }

    fn hsp_end(&mut self, format: usize, hsp: HspIndex) {
        let Some((begun, start)) = self.open[format].take() else {
            self.fault
                .get_or_insert(format!("HSP {hsp} ended without a start"));
            return;
        };
        if begun != hsp {
            self.fault
                .get_or_insert(format!("HSP {begun} started but HSP {hsp} ended"));
            return;
        }
        let range = (start, self.position(format));
        let stored = match self.streams[format] {
            6 => store(&mut self.out6, hsp, range),
            0 => store(&mut self.out0, hsp, range).and_then(|()| match self.subject {
                Some(subject) => store(&mut self.out0_subject, hsp, subject),
                None => Ok(()),
            }),
            _ => Ok(()),
        };
        if let Err(fault) = stored {
            self.fault.get_or_insert(fault);
        }
    }

    fn subject_begin(&mut self, format: usize, _first_hsp: HspIndex) {
        if self.streams[format] == 0 {
            self.subject_start = Some(self.position(format));
        }
    }

    fn subject_end(&mut self, format: usize, _first_hsp: HspIndex) {
        if self.streams[format] == 0 {
            if let Some(start) = self.subject_start.take() {
                self.subject = Some((start, self.position(format)));
            }
        }
    }
}

/// Runs one search and emits its streams. On failure nothing more is emitted and the
/// host discards what it received.
pub fn run(
    words: &[&str],
    query: &[fasta::Record],
    subject: &[fasta::Record],
) -> Result<(), String> {
    let (program, command) = parse(words)?;
    let formats = program.formats();
    let mut writers: Vec<StreamWriter> = formats.iter().map(|&f| StreamWriter::new(f)).collect();
    let mut diagnostics = StreamWriter::new(STREAM_DIAGNOSTICS);
    let mut recorder = RangeRecorder {
        streams: formats.to_vec(),
        positions: writers.iter().map(StreamWriter::position).collect(),
        open: vec![None; formats.len()],
        out6: Vec::new(),
        out0: Vec::new(),
        out0_subject: Vec::new(),
        subject_start: None,
        subject: None,
        fault: None,
    };
    let outfmts: Vec<String> = formats.iter().map(u32::to_string).collect();
    let mut collected: Option<Vec<PairwiseHit>> = None;
    {
        let mut on_hits = |hits: &[PairwiseHit]| collected = Some(hits.to_vec());
        let mut outputs = ReportOutputs {
            formats: outfmts
                .iter()
                .zip(writers.iter_mut())
                .map(|(outfmt, writer)| FormatOutput {
                    outfmt,
                    sink: OutputSink::Writer(writer),
                })
                .collect(),
            diagnostics: &mut diagnostics,
            hits: Some(&mut on_hits),
            observer: Some(&mut recorder),
        };
        let result = match command {
            Commands::Blastp(args) => run_local_blastp(args, query, subject, "", "", &mut outputs),
            Commands::Tblastn(args) => run_local_tblastn(args, query, subject, &mut outputs),
            Commands::Blastn(args) => run_local_blastn(args, query, subject, &mut outputs),
            Commands::Tblastx(args) => run_local_tblastx(args, query, subject, &mut outputs),
            Commands::Blastx(_) => unreachable!("rejected by Program::parse"),
        };
        result.map_err(|error| format!("{error:#}"))?;
    }
    if let Some((hsp, _)) = recorder.open.iter().flatten().next() {
        recorder
            .fault
            .get_or_insert(format!("HSP {hsp} started but did not end"));
    }
    if let Some(fault) = recorder.fault.take() {
        return Err(format!("internal error in the output ranges: {fault}"));
    }
    for writer in &mut writers {
        writer.finish();
    }
    diagnostics.finish();
    if let Some(hits) = collected {
        emit_hit_records(&hits, &recorder);
    }
    Ok(())
}

/// Emits one JSON object per HSP on stream 1 (docs/web/abi_v2.md §8).
fn emit_hit_records(hits: &[PairwiseHit], recorder: &RangeRecorder) {
    let mut lines = String::new();
    let mut ranks: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
    for (index, record) in hits.iter().enumerate() {
        let hit = &record.hit;
        let rank = ranks.entry(hit.q_idx).or_default();
        let _ = write!(
            lines,
            "{{\"index\":{index},\"q_idx\":{},\"s_idx\":{},\"rank\":{rank},\"raw_score\":{},\"bit_score\":",
            hit.q_idx, hit.s_idx, hit.raw_score
        );
        *rank += 1;
        json::number(&mut lines, hit.bit_score);
        lines.push_str(",\"e_value\":");
        json::number(&mut lines, hit.e_value);
        let _ = write!(
            lines,
            ",\"q_start\":{},\"q_end\":{},\"s_start\":{},\"s_end\":{},\"query_frame\":",
            hit.q_start, hit.q_end, hit.s_start, hit.s_end
        );
        json::optional(&mut lines, record.query_frame);
        lines.push_str(",\"subject_frame\":");
        json::optional(&mut lines, record.subject_frame);
        lines.push_str(",\"subject_length\":");
        json::optional(&mut lines, record.subject_length);
        lines.push_str(",\"query_aligned\":");
        json::optional_string(&mut lines, record.query_seq.as_deref());
        lines.push_str(",\"subject_aligned\":");
        json::optional_string(&mut lines, record.subject_seq.as_deref());
        lines.push_str(",\"out6\":");
        json::range(&mut lines, recorder.out6.get(index).copied().flatten());
        lines.push_str(",\"out0\":");
        json::range(&mut lines, recorder.out0.get(index).copied().flatten());
        lines.push_str(",\"out0_subject\":");
        json::range(
            &mut lines,
            recorder.out0_subject.get(index).copied().flatten(),
        );
        lines.push_str("}\n");
        if lines.len() >= crate::emit::CHUNK {
            emit(STREAM_HITS, lines.as_bytes());
            lines.clear();
        }
    }
    emit(STREAM_HITS, lines.as_bytes());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn argv_is_parsed_as_on_the_command_line() {
        let words = ["blastp", "-query", "q", "-subject", "s", "-evalue", "0.1"];
        assert!(matches!(
            parse(&words),
            Ok((Program::Blastp, Commands::Blastp(_)))
        ));
        for owned in ["-outfmt", "-out", "-outfmt=6"] {
            let error = parse(&["blastp", "-query", "q", "-subject", "s", owned, "6"]).unwrap_err();
            assert!(error.contains("is not accepted"), "{owned}: {error}");
        }
        assert!(parse(&["blastx", "-query", "q", "-subject", "s"]).is_err());
        // An unknown program and a parser error are the CLI's messages.
        for words in [&["nosuch"][..], &["tblastx", "-query", "q"][..]] {
            let cli = LOSAT::cli::try_parse_from::<Cli, _, _>(
                std::iter::once("LOSAT").chain(words.iter().copied()),
            )
            .map(|_| ())
            .map_err(|error| LOSAT::cli::render_message(&error));
            assert_eq!(parse(words).map(|_| ()), cli, "{words:?}");
        }
        assert_eq!(Program::Tblastx.formats(), &[0, 6, 7]);
        // validate rejects BLASTN scoring that NCBI rejects, with NCBI's message.
        // validate reads BLASTN's -dust as the run does.
        let error =
            validate(&["blastn", "-query", "q", "-subject", "s", "-dust", "20 64"]).unwrap_err();
        assert!(
            error.starts_with("BLAST query/options error: Invalid number of arguments"),
            "{error}"
        );
        let blastn = ["blastn", "-query", "q", "-subject", "s", "-task", "blastn"];
        assert!(validate(&blastn).is_ok());
        let unsupported = [&blastn[..], &["-reward", "1", "-penalty", "-6"]].concat();
        assert!(
            parse(&unsupported).is_ok(),
            "run leaves the check to the engine"
        );
        let error =
            validate(&[&blastn[..], &["-gapopen", "0", "-gapextend", "0"]].concat()).unwrap_err();
        assert!(
            error.starts_with("BLAST query/options error: Greedy extension"),
            "{error}"
        );
        let error = validate(&unsupported).unwrap_err();
        assert!(
            error.starts_with("BLAST engine error: Error: Substitution scores 1 and -6"),
            "{error}"
        );
        // An error of the CLI parser keeps the CLI's message.
        for words in [
            &["blastp", "-query", "q", "-subject", "s", "-evalue"][..],
            &[
                "blastn",
                "-query",
                "q",
                "-subject",
                "s",
                "-nosuchoption",
                "1",
            ][..],
            &[
                "tblastn",
                "-query",
                "q",
                "-subject",
                "s",
                "-db_gencode",
                "7",
            ][..],
        ] {
            let cli = LOSAT::cli::try_parse_from::<Cli, _, _>(
                ["LOSAT"].into_iter().chain(words.iter().copied()),
            )
            .unwrap_err();
            assert_eq!(
                parse(words).unwrap_err(),
                LOSAT::cli::render_message(&cli),
                "{words:?}"
            );
        }
    }

    fn recorder() -> RangeRecorder {
        RangeRecorder {
            streams: vec![0, 6],
            positions: vec![Arc::new(AtomicU64::new(0)), Arc::new(AtomicU64::new(0))],
            open: vec![None, None],
            out6: Vec::new(),
            out0: Vec::new(),
            out0_subject: Vec::new(),
            subject_start: None,
            subject: None,
            fault: None,
        }
    }

    // Observer events that break the formatter contract fail the run instead of
    // producing wrong ranges.
    #[test]
    fn unmatched_or_repeated_hsp_events_are_faults() {
        let mut unmatched = recorder();
        unmatched.hsp_end(1, 0);
        assert!(unmatched.fault.is_some());

        let mut crossed = recorder();
        crossed.hsp_begin(1, 0);
        crossed.hsp_end(1, 1);
        assert!(crossed.fault.is_some());

        let mut repeated = recorder();
        for _ in 0..2 {
            repeated.hsp_begin(1, 3);
            repeated.hsp_end(1, 3);
        }
        assert!(repeated.fault.is_some());

        let mut doubled = recorder();
        doubled.hsp_begin(1, 0);
        doubled.hsp_begin(1, 1);
        assert!(doubled.fault.is_some());

        let mut valid = recorder();
        valid.hsp_begin(1, 0);
        valid.hsp_end(1, 0);
        valid.subject_begin(0, 0);
        valid.subject_end(0, 0);
        valid.hsp_begin(0, 0);
        valid.hsp_end(0, 0);
        assert!(valid.fault.is_none());
        assert_eq!(
            (valid.out6.len(), valid.out0.len(), valid.out0_subject.len()),
            (1, 1, 1)
        );
    }
}
