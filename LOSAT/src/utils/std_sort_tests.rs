//! `std_sort` against libstdc++ 13's `std::sort`: the permutations that g++ 13.3 gave for the
//! inputs of `LOSAT/tests/unit/std_sort_vectors.tsv` (made by a small C++ program outside the
//! repository, `sortgen.cpp` of the task folder of item PS16a; g++ is a test-data generator here,
//! as NCBI is an oracle). Each row: a case name, the keys `k1:k2` of the elements in input
//! order, and the original indices in the order `std::sort` leaves them with rpsbproc's
//! comparison (`k1` ascending, then `k2` descending; `TDomAlignOrderByEValue`,
//! RpsbProc/src/common/cdalignproc_base.hpp:62-70). The cases: random keys with many ties for
//! every size from 0 to 200, random keys of a wide range, sorted, reversed and all-equal inputs,
//! McIlroy's adversarial inputs (`adversary*`), which drive the introsort to its depth limit (the
//! heap sort of `__partial_sort`) with distinct keys, the same inputs with ties (`advties*`, which
//! stay in the introsort), and `heapties200`, an adversarial input with 32 groups of ties that
//! reaches the heap sort (made by the PS16a review, `heapties.py`, its permutation g++'s; a
//! Python transcription of libstdc++, `pysort.py` there, finds these 5 rows reaching the heap).

use super::std_sort;

#[derive(Clone, Copy, Debug)]
struct Element {
    k1: i32,
    k2: i32,
    index: usize,
}

#[test]
fn std_sort_gives_libstdcxx_permutations() {
    let table = include_str!("../../tests/unit/std_sort_vectors.tsv");
    let mut cases = 0;
    let mut heap_cases = 0;
    for line in table.lines() {
        let mut fields = line.split('\t');
        let name = fields.next().unwrap();
        let keys = fields.next().unwrap_or("");
        let expected: Vec<usize> = fields
            .next()
            .unwrap_or("")
            .split(',')
            .filter(|text| !text.is_empty())
            .map(|text| text.parse().unwrap())
            .collect();
        let mut elements: Vec<Element> = keys
            .split(',')
            .filter(|text| !text.is_empty())
            .enumerate()
            .map(|(index, pair)| {
                let (k1, k2) = pair.split_once(':').unwrap();
                Element {
                    k1: k1.parse().unwrap(),
                    k2: k2.parse().unwrap(),
                    index,
                }
            })
            .collect();
        assert_eq!(elements.len(), expected.len(), "{name}");
        std_sort(&mut elements, &mut |a: &Element, b: &Element| {
            if a.k1 < b.k1 {
                true
            } else if a.k1 > b.k1 {
                false
            } else {
                a.k2 > b.k2
            }
        });
        let got: Vec<usize> = elements.iter().map(|element| element.index).collect();
        assert_eq!(got, expected, "{name}");
        cases += 1;
        if name.starts_with("adversary") || name.starts_with("heapties") {
            heap_cases += 1;
        }
    }
    assert_eq!(cases, 272);
    assert_eq!(heap_cases, 5);
}

/// The adversarial inputs reach the heap sort: with the comparisons counted, a sort of an
/// adversarial input of 2000 distinct keys stays far below the quadratic count of a quicksort
/// without the depth limit (about n^2/2 = 2,000,000 comparisons).
#[test]
fn std_sort_adversary_is_bounded_by_the_heap_sort() {
    let table = include_str!("../../tests/unit/std_sort_vectors.tsv");
    let line = table
        .lines()
        .find(|line| line.starts_with("adversary2000\t"))
        .unwrap();
    let keys: Vec<i32> = line
        .split('\t')
        .nth(1)
        .unwrap()
        .split(',')
        .map(|pair| pair.split_once(':').unwrap().0.parse().unwrap())
        .collect();
    let mut values = keys.clone();
    let mut comparisons = 0usize;
    std_sort(&mut values, &mut |a: &i32, b: &i32| {
        comparisons += 1;
        a < b
    });
    let mut sorted = keys;
    sorted.sort_unstable();
    assert_eq!(values, sorted);
    assert!(comparisons < 200_000, "{comparisons}");
}
