use anyhow::Result;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct Hit {
    pub query_id: String, // インデックスではなくIDを持つように変更（並列処理後の出力順序制御のため）
    pub subject_id: String,
    pub identity: f64,
    pub length: usize,
    pub mismatch: usize,
    pub gapopen: usize,
    pub q_start: usize,
    pub q_end: usize,
    pub s_start: usize,
    pub s_end: usize,
    pub e_value: f64,
    pub bit_score: f64,
    // Fields for NCBI-style output ordering (not printed, used for sorting)
    /// Query index (input order) - NCBI uses query order for grouping
    pub q_idx: u32,
    /// Subject index (oid) - NCBI uses oid for tie-breaking in s_EvalueCompareHSPLists
    pub s_idx: u32,
    /// Raw alignment score - NCBI ScoreCompareHSPs uses raw score, not bit score
    pub raw_score: i32,
}

// =============================================================================
// NCBI-style comparators for outfmt6 output ordering
// Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c
// =============================================================================

/// Compare two evalues, treating both as equal if they're close enough to zero.
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c s_EvalueComp()
/// ```c
/// static int s_EvalueComp(double evalue1, double evalue2) {
///     const double epsilon = 1.0e-180;
///     if (evalue1 < epsilon && evalue2 < epsilon) { return 0; }
///     if (evalue1 < evalue2) { return -1; }
///     else if (evalue1 > evalue2) { return 1; }
///     else { return 0; }
/// }
/// ```
#[inline]
fn evalue_comp(evalue1: f64, evalue2: f64) -> Ordering {
    const EPSILON: f64 = 1.0e-180;
    if evalue1 < EPSILON && evalue2 < EPSILON {
        Ordering::Equal
    } else if evalue1 < evalue2 {
        Ordering::Less
    } else if evalue1 > evalue2 {
        Ordering::Greater
    } else {
        Ordering::Equal
    }
}

/// NCBI ScoreCompareHSPs - compare two HSPs by score and coordinates.
/// Used to sort HSPs within a subject (BlastHSPList).
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c ScoreCompareHSPs()
/// ```c
/// int ScoreCompareHSPs(const void* h1, const void* h2) {
///    BlastHSP* hsp1,* hsp2;
///    int result = 0;
///    hsp1 = *((BlastHSP**) h1);
///    hsp2 = *((BlastHSP**) h2);
///    if (0 == (result = BLAST_CMP(hsp2->score,          hsp1->score)) &&
///        0 == (result = BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset)) &&
///        0 == (result = BLAST_CMP(hsp2->subject.end,    hsp1->subject.end)) &&
///        0 == (result = BLAST_CMP(hsp1->query  .offset, hsp2->query  .offset))) {
///        result = BLAST_CMP(hsp2->query.end, hsp1->query.end);
///    }
///    return result;
/// }
/// ```
/// Order: score DESC → s_start ASC → s_end DESC → q_start ASC → q_end DESC
pub fn score_compare_hsps(a: &Hit, b: &Hit) -> Ordering {
    // score DESC (BLAST_CMP(hsp2->score, hsp1->score))
    match b.raw_score.cmp(&a.raw_score) {
        Ordering::Equal => {}
        ord => return ord,
    }
    // s_start ASC (BLAST_CMP(hsp1->subject.offset, hsp2->subject.offset))
    match a.s_start.cmp(&b.s_start) {
        Ordering::Equal => {}
        ord => return ord,
    }
    // s_end DESC (BLAST_CMP(hsp2->subject.end, hsp1->subject.end))
    match b.s_end.cmp(&a.s_end) {
        Ordering::Equal => {}
        ord => return ord,
    }
    // q_start ASC (BLAST_CMP(hsp1->query.offset, hsp2->query.offset))
    match a.q_start.cmp(&b.q_start) {
        Ordering::Equal => {}
        ord => return ord,
    }
    // q_end DESC (BLAST_CMP(hsp2->query.end, hsp1->query.end))
    b.q_end.cmp(&a.q_end)
}

/// Subject group for NCBI-style output ordering.
/// Represents all HSPs for a single subject (BlastHSPList equivalent).
struct SubjectGroup {
    s_idx: u32,
    best_evalue: f64,
    best_score: i32,
    hits: Vec<Hit>,
}

/// NCBI s_EvalueCompareHSPLists - compare two subjects by best e-value/score/oid.
/// Used to sort subjects within a query.
/// Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c s_EvalueCompareHSPLists()
/// ```c
/// static int s_EvalueCompareHSPLists(const void* v1, const void* v2) {
///    BlastHSPList* h1,* h2;
///    int retval = 0;
///    h1 = *(BlastHSPList**) v1;
///    h2 = *(BlastHSPList**) v2;
///    if (h1->hspcnt == 0 && h2->hspcnt == 0) return 0;
///    else if (h1->hspcnt == 0) return 1;
///    else if (h2->hspcnt == 0) return -1;
///    if ((retval = s_EvalueComp(h1->best_evalue, h2->best_evalue)) != 0)
///       return retval;
///    if (h1->hsp_array[0]->score > h2->hsp_array[0]->score) return -1;
///    if (h1->hsp_array[0]->score < h2->hsp_array[0]->score) return 1;
///    return BLAST_CMP(h2->oid, h1->oid);
/// }
/// ```
/// Order: best_evalue ASC → best_score DESC → oid DESC
fn compare_subject_groups(a: &SubjectGroup, b: &SubjectGroup) -> Ordering {
    // Empty groups go to the end
    if a.hits.is_empty() && b.hits.is_empty() {
        return Ordering::Equal;
    } else if a.hits.is_empty() {
        return Ordering::Greater;
    } else if b.hits.is_empty() {
        return Ordering::Less;
    }

    // best_evalue ASC (s_EvalueComp)
    match evalue_comp(a.best_evalue, b.best_evalue) {
        Ordering::Equal => {}
        ord => return ord,
    }

    // best_score DESC
    match b.best_score.cmp(&a.best_score) {
        Ordering::Equal => {}
        ord => return ord,
    }

    // oid DESC (BLAST_CMP(h2->oid, h1->oid))
    b.s_idx.cmp(&a.s_idx)
}

pub fn write_output(hits: &[Hit], out_path: Option<&PathBuf>) -> Result<()> {
    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(path) = out_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout.lock()))
    };

    for hit in hits {
        writeln!(
            writer,
            "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1e}\t{:.1}",
            hit.query_id,
            hit.subject_id,
            hit.identity,
            hit.length,
            hit.mismatch,
            hit.gapopen,
            hit.q_start,
            hit.q_end,
            hit.s_start,
            hit.s_end,
            hit.e_value,
            hit.bit_score
        )?;
    }
    Ok(())
}

/// Write output with NCBI-style ordering:
/// 1. Group by query (input order)
/// 2. Within query, sort subjects by s_EvalueCompareHSPLists (best_evalue ASC → best_score DESC → oid DESC)
/// 3. Within subject, sort HSPs by ScoreCompareHSPs (score DESC → s_start ASC → s_end DESC → q_start ASC → q_end DESC)
///
/// Reference:
/// - BLAST_LinkHsps() calls Blast_HSPListSortByScore() after linking
/// - BlastHitList sorts subjects by s_EvalueCompareHSPLists
/// - Final output iterates query → subject → HSP
pub fn write_output_ncbi_order(mut hits: Vec<Hit>, out_path: Option<&PathBuf>) -> Result<()> {
    let stdout = io::stdout();
    let mut writer: Box<dyn Write> = if let Some(path) = out_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(BufWriter::new(stdout.lock()))
    };

    if hits.is_empty() {
        return Ok(());
    }

    // Step 1: Group by query (preserving input order via q_idx)
    // Collect unique queries in order
    let mut query_order: Vec<u32> = Vec::new();
    let mut seen_queries: std::collections::HashSet<u32> = std::collections::HashSet::new();
    for h in &hits {
        if !seen_queries.contains(&h.q_idx) {
            seen_queries.insert(h.q_idx);
            query_order.push(h.q_idx);
        }
    }

    // Group hits by (q_idx, s_idx)
    let mut query_subject_hits: HashMap<(u32, u32), Vec<Hit>> = HashMap::new();
    for h in hits.drain(..) {
        query_subject_hits
            .entry((h.q_idx, h.s_idx))
            .or_default()
            .push(h);
    }

    // Step 2: For each query, build subject groups and sort
    for &q_idx in &query_order {
        // Collect all subjects for this query
        let subject_indices: Vec<u32> = query_subject_hits
            .keys()
            .filter(|(qidx, _)| *qidx == q_idx)
            .map(|(_, sidx)| *sidx)
            .collect::<std::collections::HashSet<_>>()
            .into_iter()
            .collect();

        // Build subject groups
        let mut subject_groups: Vec<SubjectGroup> = subject_indices
            .into_iter()
            .filter_map(|s_idx| {
                let key = (q_idx, s_idx);
                query_subject_hits.remove(&key).map(|mut hsp_hits| {
                    // Sort HSPs within subject by ScoreCompareHSPs
                    hsp_hits.sort_by(score_compare_hsps);

                    // Compute best_evalue and best_score
                    let best_evalue = hsp_hits
                        .iter()
                        .map(|h| h.e_value)
                        .min_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal))
                        .unwrap_or(f64::MAX);
                    let best_score = hsp_hits.iter().map(|h| h.raw_score).max().unwrap_or(0);

                    SubjectGroup {
                        s_idx,
                        best_evalue,
                        best_score,
                        hits: hsp_hits,
                    }
                })
            })
            .collect();

        // Sort subjects by s_EvalueCompareHSPLists
        subject_groups.sort_by(compare_subject_groups);

        // Step 3: Output in order
        for group in subject_groups {
            for hit in group.hits {
                writeln!(
                    writer,
                    "{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1e}\t{:.1}",
                    hit.query_id,
                    hit.subject_id,
                    hit.identity,
                    hit.length,
                    hit.mismatch,
                    hit.gapopen,
                    hit.q_start,
                    hit.q_end,
                    hit.s_start,
                    hit.s_end,
                    hit.e_value,
                    hit.bit_score
                )?;
            }
        }
    }

    Ok(())
}
