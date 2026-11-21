use super::build_disc_mb_lookup;
use crate::algorithm::blastn::disc_lookup::{
    compute_discontiguous_index, DiscTemplateType, DiscWordType,
};
use crate::utils::dust::MaskedInterval;

fn index_at(seq: &[u8], start: usize, len: usize, t: DiscTemplateType) -> u64 {
    let mut accum = 0u64;
    for &b in &seq[start..start + len] {
        accum = (accum << 2) | b as u64;
    }
    compute_discontiguous_index(accum, t) as u64
}

fn chain(lookup: &super::TwoStageLookup, index: u64, second: bool) -> Vec<u32> {
    let mut hits = Vec::new();
    if second {
        lookup.for_each_hit2(index, |q| hits.push(q));
    } else {
        lookup.for_each_hit(index, |q| hits.push(q));
    }
    hits
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:741-813
#[test]
fn words_are_indexed_newest_first_with_one_based_starts() {
    // Two contexts in the concatenated query (offsets 0 and 41); the second repeats
    // the first's 18-mer at 0 and has an ambiguity code (14) inside its copy at 20.
    let word: Vec<u8> = (0..18).map(|i| (i * 7 % 4) as u8).collect();
    let mut first = word.clone();
    first.extend([1, 2, 3]);
    let mut second = word.clone();
    second.extend([0, 0]);
    second.extend(&word);
    second[20 + 5] = 14;
    let lookup = build_disc_mb_lookup(
        &[&first, &second],
        &[0, 41],
        11,
        18,
        DiscWordType::TwoTemplates,
        &[Vec::new(), Vec::new()],
        1000,
    );
    let disc = lookup.disc().expect("discontiguous table");
    assert_eq!(disc.template_type, DiscTemplateType::T11_18Coding);
    assert_eq!(disc.second_template_type, DiscTemplateType::T11_18Optimal);
    assert_eq!((lookup.word_length(), lookup.lut_word_length()), (18, 18));
    let coding = index_at(&word, 0, 18, DiscTemplateType::T11_18Coding);
    let optimal = index_at(&word, 0, 18, DiscTemplateType::T11_18Optimal);
    // The word at 0 of each context; the copy at 20 holds the ambiguity code.
    assert_eq!(chain(&lookup, coding, false), vec![42, 1]);
    assert_eq!(chain(&lookup, optimal, true), vec![42, 1]);
    // A masked range is not indexed.
    let masked = build_disc_mb_lookup(
        &[&first, &second],
        &[0, 41],
        11,
        18,
        DiscWordType::Coding,
        &[vec![MaskedInterval::new(3, 4)], Vec::new()],
        1000,
    );
    assert_eq!(chain(&masked, coding, false), vec![42]);
    assert!(masked.disc().is_some_and(|d| !d.two_templates));
    assert!(chain(&masked, optimal, true).is_empty());
}
