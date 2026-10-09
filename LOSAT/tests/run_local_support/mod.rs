//! The harness shared by the `run_local` tests of every program: one sink per requested
//! format, the collected hit records, and the byte range that the observer reports for
//! every HSP in every format.

use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::time::{SystemTime, UNIX_EPOCH};

use bio::io::fasta;
use LOSAT::api::local_blast::{FormatObserver, FormatOutput, HspIndex, OutputSink, ReportOutputs};
use LOSAT::blastinput::fasta_reader::{read_all, FastaInputSource, FastaRecord, ReaderConfig};
use LOSAT::report::PairwiseHit;

/// The sequence of the first record of a FASTA file in `LOSAT/tests/fasta`.
pub fn fixture_sequence(name: &str) -> Vec<u8> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name);
    read_records(&path)
        .into_iter()
        .next()
        .expect("fixture record")
        .seq()
        .to_vec()
}

/// A temporary FASTA file, removed when dropped.
pub struct TempFasta(pub PathBuf);

impl TempFasta {
    /// Writes `records` (id, sequence) with 60 residues per line.
    pub fn new(name: &str, records: &[(&str, &[u8])]) -> Self {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "losat_run_local_{}_{nanos}_{name}",
            std::process::id()
        ));
        let mut text = Vec::new();
        for (id, sequence) in records {
            text.extend_from_slice(format!(">{id}\n").as_bytes());
            for line in sequence.chunks(60) {
                text.extend_from_slice(line);
                text.push(b'\n');
            }
        }
        std::fs::write(&path, text).expect("write temporary FASTA");
        Self(path)
    }
}

impl Drop for TempFasta {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

pub fn read_records(path: &Path) -> Vec<fasta::Record> {
    fasta::Reader::from_file(path)
        .expect("FASTA")
        .records()
        .collect::<Result<_, _>>()
        .expect("records")
}

/// The records of a FASTA file as NCBI's reader reads them (`config`), the records of the
/// programs that read with `fasta_reader`.
#[allow(dead_code)]
pub fn reader_records(path: &Path, config: ReaderConfig) -> Vec<FastaRecord> {
    let file = std::fs::File::open(path).expect("FASTA");
    read_all(&mut FastaInputSource::from_file(file, config), &mut |_| {
        Ok(())
    })
    .expect("records")
}

/// A sink that publishes how many bytes it holds, so the observer can read positions.
struct CountingSink {
    bytes: Vec<u8>,
    len: Arc<AtomicUsize>,
}

impl Write for CountingSink {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        self.bytes.extend_from_slice(data);
        self.len.store(self.bytes.len(), Ordering::SeqCst);
        Ok(data.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

struct RecordingObserver {
    lengths: Vec<Arc<AtomicUsize>>,
    open: Vec<Option<(HspIndex, usize)>>,
    ranges: Vec<Vec<(HspIndex, usize, usize)>>,
    open_subject: Vec<Option<(HspIndex, usize)>>,
    subjects: Vec<Vec<(HspIndex, usize, usize)>>,
}

impl FormatObserver for RecordingObserver {
    fn hsp_begin(&mut self, format: usize, hsp: HspIndex) {
        assert!(self.open[format].is_none(), "nested HSP in format {format}");
        self.open[format] = Some((hsp, self.lengths[format].load(Ordering::SeqCst)));
    }
    fn hsp_end(&mut self, format: usize, hsp: HspIndex) {
        let (begun, start) = self.open[format].take().expect("end without begin");
        assert_eq!(begun, hsp, "end for a different HSP in format {format}");
        let end = self.lengths[format].load(Ordering::SeqCst);
        self.ranges[format].push((hsp, start, end));
    }
    fn subject_begin(&mut self, format: usize, first_hsp: HspIndex) {
        assert!(self.open[format].is_none(), "subject heading inside an HSP");
        assert!(
            self.open_subject[format].is_none(),
            "nested subject heading"
        );
        self.open_subject[format] = Some((first_hsp, self.lengths[format].load(Ordering::SeqCst)));
    }
    fn subject_end(&mut self, format: usize, first_hsp: HspIndex) {
        let (begun, start) = self.open_subject[format].take().expect("end without begin");
        assert_eq!(
            begun, first_hsp,
            "end for a different subject in format {format}"
        );
        let end = self.lengths[format].load(Ordering::SeqCst);
        self.subjects[format].push((first_hsp, start, end));
    }
}

/// What one `run_local` call wrote.
pub struct Run {
    /// The bytes of each requested format, in request order.
    pub outputs: Vec<Vec<u8>>,
    pub diagnostics: Vec<u8>,
    pub hits: Vec<PairwiseHit>,
    /// Per format: (HSP index, start, end) of every row or section, in write order.
    pub ranges: Vec<Vec<(HspIndex, usize, usize)>>,
    /// Per format: (index of the subject's first HSP, start, end) of every subject
    /// heading, in write order.
    pub subjects: Vec<Vec<(HspIndex, usize, usize)>>,
}

/// Calls `search` with one sink per format. With `observe`, it also requests the hit
/// records and the observer.
pub fn run_formats(
    formats: &[&str],
    observe: bool,
    search: impl FnOnce(&mut ReportOutputs<'_>) -> anyhow::Result<()>,
) -> Run {
    let mut sinks: Vec<CountingSink> = formats
        .iter()
        .map(|_| CountingSink {
            bytes: Vec::new(),
            len: Arc::new(AtomicUsize::new(0)),
        })
        .collect();
    let mut observer = RecordingObserver {
        lengths: sinks.iter().map(|sink| sink.len.clone()).collect(),
        open: vec![None; formats.len()],
        ranges: vec![Vec::new(); formats.len()],
        open_subject: vec![None; formats.len()],
        subjects: vec![Vec::new(); formats.len()],
    };
    let collected = Arc::new(Mutex::new(Vec::new()));
    let collected_sink = collected.clone();
    let mut on_hits =
        move |hits: &[PairwiseHit]| collected_sink.lock().unwrap().extend_from_slice(hits);
    let mut diagnostics = Vec::new();
    {
        let mut outputs = ReportOutputs {
            formats: formats
                .iter()
                .zip(sinks.iter_mut())
                .map(|(outfmt, sink)| FormatOutput {
                    outfmt,
                    sink: OutputSink::Writer(sink),
                })
                .collect(),
            diagnostics: &mut diagnostics,
            hits: observe.then_some(&mut on_hits as _),
            observer: if observe { Some(&mut observer) } else { None },
        };
        search(&mut outputs).expect("run_local");
    }
    for (format, open) in observer.open.iter().enumerate() {
        assert!(open.is_none(), "HSP left open in format {format}");
    }
    for (format, open) in observer.open_subject.iter().enumerate() {
        assert!(
            open.is_none(),
            "subject heading left open in format {format}"
        );
    }
    let hits = collected.lock().unwrap().clone();
    Run {
        outputs: sinks.into_iter().map(|sink| sink.bytes).collect(),
        diagnostics,
        hits,
        ranges: observer.ranges,
        subjects: observer.subjects,
    }
}

/// The `index`-th tab-separated field of a tabular row.
pub fn field(row: &[u8], index: usize) -> String {
    let row = std::str::from_utf8(row).expect("utf-8 row").trim_end();
    row.split('\t').nth(index).expect("field").to_string()
}

// NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1100-1108
// ```c
// ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
//     // Add tab in front of field, except for the first field.
//     if (iter != m_FieldsToShow.begin())
//         m_Ostream << m_FieldDelimiter;
//     x_PrintField(*iter);
// }
// m_Ostream << "\n";
// ```
// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:1970-1973
// ```c++
// x_DisplayAlnvecInfo(out, alnvecInfo,showDefLine);
// ```
/// Checks the observer ranges of the tabular formats: the outfmt 6 rows (at position
/// `fmt6` of the request) are in HSP index order and tile the output, one line each;
/// with `fmt7`, the outfmt 7 rows are its non-comment lines and equal the outfmt 6 rows.
/// Returns the number of rows.
pub fn assert_tabular_ranges(run: &Run, fmt6: usize, fmt7: Option<usize>, context: &str) -> usize {
    let (out6, rows6) = (&run.outputs[fmt6], &run.ranges[fmt6]);
    let mut position = 0;
    for (expected_index, &(hsp, start, end)) in rows6.iter().enumerate() {
        assert_eq!(hsp, expected_index, "{context}: outfmt 6 order");
        assert_eq!(start, position, "{context}: outfmt 6 HSP {hsp} start");
        let row = &out6[start..end];
        assert_eq!(row.iter().filter(|&&byte| byte == b'\n').count(), 1);
        position = end;
    }
    assert_eq!(
        position,
        out6.len(),
        "{context}: outfmt 6 rows tile the output"
    );

    if let Some(fmt7) = fmt7 {
        let (out7, rows7) = (&run.outputs[fmt7], &run.ranges[fmt7]);
        assert_eq!(rows7.len(), rows6.len(), "{context}: outfmt 7 rows");
        for (&(hsp6, s6, e6), &(hsp7, s7, e7)) in rows6.iter().zip(rows7) {
            assert_eq!(hsp6, hsp7, "{context}: outfmt 7 order");
            assert_eq!(
                &out6[s6..e6],
                &out7[s7..e7],
                "{context}: outfmt 7 HSP {hsp7}"
            );
        }
        let data_lines = out7
            .split_inclusive(|&byte| byte == b'\n')
            .filter(|line| !line.starts_with(b"#"))
            .count();
        assert_eq!(data_lines, rows6.len(), "{context}: outfmt 7 data lines");
    }
    rows6.len()
}

/// Checks the observer ranges of outfmt 0, 6 and 7 (at the given positions of the
/// request) against the hit records: `assert_tabular_ranges`; the same index has the
/// same coordinates in the hit record and in its outfmt 6 row; and each outfmt 0
/// section is one score block with the raw score and the query and subject
/// coordinates of its HSP, inside the report block of its query.
pub fn assert_observer_ranges(run: &Run, [fmt0, fmt6, fmt7]: [usize; 3], context: &str) {
    let hits = &run.hits;
    assert!(!hits.is_empty(), "{context}: the fixture must produce hits");
    let rows = assert_tabular_ranges(run, fmt6, Some(fmt7), context);
    assert_eq!(rows, hits.len(), "{context}: outfmt 6 rows");
    for &(hsp, start, end) in &run.ranges[fmt6] {
        let hit = &hits[hsp].hit;
        let row = &run.outputs[fmt6][start..end];
        assert_eq!(
            [6, 7, 8, 9].map(|index| field(row, index)),
            [hit.q_start, hit.q_end, hit.s_start, hit.s_end].map(|value| value.to_string()),
            "{context}: outfmt 6 HSP {hsp} coordinates"
        );
    }

    let (out0, sections) = (&run.outputs[fmt0], &run.ranges[fmt0]);
    assert_eq!(sections.len(), hits.len(), "{context}: outfmt 0 sections");
    // Byte offsets of the "Query=" lines that open each query's report block.
    let query_blocks: Vec<usize> = out0
        .windows(7)
        .enumerate()
        .filter(|&(offset, window)| {
            window == b"Query= " && (offset == 0 || out0[offset - 1] == b'\n')
        })
        .map(|(offset, _)| offset)
        .collect();
    let mut seen = vec![false; hits.len()];
    for &(hsp, start, end) in sections {
        assert!(
            !std::mem::replace(&mut seen[hsp], true),
            "{context}: HSP {hsp} twice"
        );
        let section = std::str::from_utf8(&out0[start..end]).expect("utf-8 section");
        assert!(
            section.starts_with(" Score ="),
            "{context}: section {hsp}: {section:.20?}"
        );
        assert_eq!(
            section.matches(" Score =").count(),
            1,
            "{context}: section {hsp}"
        );
        assert!(
            section.contains(&format!("bits ({}),", hits[hsp].hit.raw_score)),
            "{context}: raw score of HSP {hsp}"
        );
        let hit = &hits[hsp].hit;
        assert_eq!(
            alignment_span(section, "Query"),
            (hit.q_start, hit.q_end),
            "{context}: query coordinates of HSP {hsp}"
        );
        assert_eq!(
            alignment_span(section, "Sbjct"),
            (hit.s_start, hit.s_end),
            "{context}: subject coordinates of HSP {hsp}"
        );
        // The section lies in the report block of its own query.
        let q_idx = hit.q_idx as usize;
        assert!(
            query_blocks[q_idx] < start,
            "{context}: HSP {hsp} before its query block"
        );
        if let Some(&next) = query_blocks.get(q_idx + 1) {
            assert!(end <= next, "{context}: HSP {hsp} after its query block");
        }
    }

    // One heading per subject of each query, just before the section of its first HSP.
    let group_starts: Vec<HspIndex> = (0..hits.len())
        .filter(|&index| {
            index == 0
                || (hits[index].hit.q_idx, hits[index].hit.s_idx)
                    != (hits[index - 1].hit.q_idx, hits[index - 1].hit.s_idx)
        })
        .collect();
    let headings = &run.subjects[fmt0];
    assert_eq!(
        headings
            .iter()
            .map(|&(first, _, _)| first)
            .collect::<Vec<_>>(),
        group_starts,
        "{context}: one heading per subject, marked with its first HSP"
    );
    for &(first, start, end) in headings {
        let heading = std::str::from_utf8(&out0[start..end]).expect("utf-8 heading");
        assert!(
            heading.starts_with("> "),
            "{context}: heading of HSP {first}: {heading:.20?}"
        );
        let section_start = sections
            .iter()
            .find(|&&(hsp, _, _)| hsp == first)
            .map(|&(_, start, _)| start)
            .expect("section of the first HSP");
        assert_eq!(
            end, section_start,
            "{context}: heading of HSP {first} ends at its section"
        );
    }
}

// NCBI reference: ncbi-blast/c++/src/objtools/align_format/showalign.cpp:1970-1973
// ```c++
// x_DisplayAlnvecInfo(out, alnvecInfo,showDefLine);
// ```
/// The first and last coordinates on the `label` ("Query" or "Sbjct") lines of one
/// outfmt 0 alignment section.
fn alignment_span(section: &str, label: &str) -> (usize, usize) {
    let numbers: Vec<Vec<usize>> = section
        .lines()
        .filter(|line| line.starts_with(label))
        .map(|line| {
            line.split_whitespace()
                .filter_map(|word| word.parse().ok())
                .collect()
        })
        .collect();
    let first = *numbers
        .first()
        .and_then(|line| line.first())
        .expect("first coordinate");
    let last = *numbers
        .last()
        .and_then(|line| line.last())
        .expect("last coordinate");
    (first, last)
}
