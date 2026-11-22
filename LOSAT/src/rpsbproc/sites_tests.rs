//! `MotifState::motif_check` (basedata.cpp:669-740): the process-wide list of the first motif
//! site's positions, the Owner's rule for the read past its end, the parse error and the
//! residue match. The sites come from a small data directory written here.

use super::*;
use crate::rpsbproc::data::{DataFiles, Site};

/// A data directory with the model `cd10001` and sites `A`..`E` (the text of the motif and the
/// positions in the table below).
fn load(sites: &[(&str, &str, &str)]) -> CddData {
    let dir = crate::utils::test_dir::temp_dir().join(format!(
        "losat_rpsbproc_sites_{}_{}",
        crate::utils::test_dir::process_id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    let write = |name: &str, text: &str| std::fs::write(dir.join(name), text).unwrap();
    write("cddid.tbl", "10001\tcd10001\tAAA\tdefinition\t100\n");
    write("cdtrack.txt", "");
    write("family_superfamily_links", "");
    write("bitscore_specific.txt", "");
    write("cddannot_generic.dat", "");
    let mut annot = String::new();
    for (title, motif, positions) in sites {
        annot.push_str(&format!(
            "10001\tcd10001\tAAA\t1\t{title}\t{motif}\t1\t0\t0\t{positions}\t1\n"
        ));
    }
    write("cddannot.dat", &annot);
    let file = |name: &str| format!("{}/{name}", dir.display());
    let data = CddData::load(&DataFiles {
        cddids: file("cddid.tbl"),
        spfeats: file("cddannot.dat"),
        genfeats: file("cddannot_generic.dat"),
        bscores: file("bitscore_specific.txt"),
        cdtrack: file("cdtrack.txt"),
        famlinks: file("family_superfamily_links"),
    })
    .unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
    data
}

fn mapped(pairs: &[(i64, i64)]) -> Vec<ResiduePos> {
    pairs
        .iter()
        .map(|&(curr, ori)| ResiduePos { curr, ori })
        .collect()
}

fn site(data: &CddData, index: usize) -> &Site {
    &data.sites(0, true)[index]
}

fn site_of(data: &CddData) -> &Site {
    site(data, 0)
}

#[test]
fn the_first_motif_site_sets_the_positions_for_every_later_check() {
    let data = load(&[
        ("A", "H-H-D", "10,11,12"),
        ("B", "H-K-Y", "20,21,22"),
        ("C", "0", "30,31"),
    ]);
    let residues = b"AAAAAHHDAAAAKAAAA";
    let mut state = MotifState::default();
    // A site without a motif passes and does not fill the list.
    let map_c = mapped(&[(1, 30), (2, 31)]);
    assert_eq!(
        state.motif_check(&data, site(&data, 2), &map_c, residues),
        Ok(0)
    );
    assert!(state.positions.is_empty());
    // The first motif site: its own positions, the residues H H D match.
    let map_a = mapped(&[(5, 10), (6, 11), (7, 12)]);
    assert_eq!(
        state.motif_check(&data, site(&data, 0), &map_a, residues),
        Ok(0)
    );
    assert_eq!(state.positions, [10, 11, 12]);
    // A later motif site of other positions: its first mapped position is past the list's.
    let map_b = mapped(&[(5, 20), (6, 21), (7, 22)]);
    assert_eq!(
        state.motif_check(&data, site(&data, 1), &map_b, residues),
        Ok(1)
    );
    // The same site checked first (a fresh run) passes only if its residues fit: H K Y at 5..7
    // are `HHD`, which does not fit `H-K-Y`.
    let mut fresh = MotifState::default();
    assert_eq!(
        fresh.motif_check(&data, site(&data, 1), &map_b, residues),
        Ok(2)
    );
    assert_eq!(fresh.positions, [20, 21, 22]);
    let residues_hky = b"AAAAAHKYAAAA";
    assert_eq!(
        MotifState::default().motif_check(&data, site(&data, 1), &map_b, residues_hky),
        Ok(0)
    );
    // Without residues only the positions are checked (`aaSeq` empty).
    assert_eq!(
        MotifState::default().motif_check(&data, site(&data, 1), &map_b, b""),
        Ok(0)
    );
}

// DECISIONS 25 Q2 = A (2026-11-11): the compare of `rMappedRes[mapIdx].ori` against
// `st_motif_pos[mtfIdx]` with `mtfIdx` at the end of the list reads heap bytes in NCBI; LOSAT fails
// the check.
#[test]
fn a_minimal_map_longer_than_the_stored_positions_fails_the_check() {
    let data = load(&[
        ("A", "H-H-D", "10,11,12"),
        ("D", "[HK]-x(4)-D", "10,11,12,13,14,15"),
    ]);
    let residues = b"AAAAAHHDAAAAAAAAAAAA";
    let mut state = MotifState::default();
    let map_a = mapped(&[(5, 10), (6, 11), (7, 12)]);
    assert_eq!(
        state.motif_check(&data, site(&data, 0), &map_a, residues),
        Ok(0)
    );
    assert_eq!(state.past_end_reads, 0);
    // The minimal map of `[HK]-x(4)-D` has 6 positions; the first three equal the stored ones,
    // the fourth asks for `st_motif_pos[3]`.
    let map_d = mapped(&[(5, 10), (6, 11), (7, 12), (8, 13), (9, 14), (10, 15)]);
    assert_eq!(
        state.motif_check(&data, site(&data, 1), &map_d, residues),
        Ok(1)
    );
    assert_eq!(state.past_end_reads, 1);
    // With no mapped residue left the loop ends first ("exhausted, fail"): not a read.
    let short = mapped(&[(5, 10), (6, 11), (7, 12)]);
    assert_eq!(
        state.motif_check(&data, site(&data, 1), &short, residues),
        Ok(1)
    );
    assert_eq!(state.past_end_reads, 1);
}

#[test]
fn a_motif_that_does_not_parse_is_an_error_with_its_position() {
    let data = load(&[("A", "K#R", "10,11,12")]);
    let mut state = MotifState::default();
    let map_a = mapped(&[(5, 10), (6, 11), (7, 12)]);
    assert_eq!(
        state.motif_check(&data, site(&data, 0), &map_a, b"AAAAAKKR"),
        Err(MotifError::Parse {
            motif: b"K#R".to_vec(),
            position: 1
        })
    );
    // The list was filled before the parse (basedata.cpp:675-676).
    assert_eq!(state.positions, [10, 11, 12]);
}

// Review oracles R1, R4, R20, R21 and review 2 oracles O2 (rd15), O3 (rd19, rd20)
// (`mfirst.asn`'s catalytic site, its positions mapped to the residues `HKY` of a query of 185):
// the read past the residues that does not decide NCBI's answer passes, the one that does is
// refused; the minimal maps that NCBI's string cannot hold end as NCBI's exceptions, and where
// NCBI holds 1 GiB or more at once (one allocation, the string's old and new buffers, or the
// match's records and their copies) the check is refused.
#[test]
fn the_motif_check_of_the_review_oracles() {
    let data = load(&[
        ("A", "Hx(0,10)K", "10,11,12"),
        ("B", "Hx(0,10)W", "10,11,12"),
        ("C", "Hx(99999999999999999999)K", "10,11,12"),
        ("D", "Hx(200000000000000)K", "10,11,12"),
        ("E", "Hx(1073741823)K", "10,11,12"),
        ("F", "Hx(4611686018427387903)K", "10,11,12"),
        ("G", "Hx(536870910)K", "10,11,12"),
        ("H", "H-{<}(0,67108863)-{<}(0,67108863)-K", "10,11,12"),
    ]);
    let mut residues = vec![b'A'; 185];
    residues[61..64].copy_from_slice(b"HKY");
    let map = mapped(&[(61, 10), (62, 11), (63, 12)]);
    let check = |index: usize| {
        MotifState::default().motif_check(&data, site(&data, index), &map, &residues)
    };
    assert_eq!(check(0), Ok(0));
    assert_eq!(check(1), Err(MotifError::Match(MatchError::ReadPast)));
    assert_eq!(
        check(2),
        Err(MotifError::Map(MapEnd::Throw(
            super::super::prosite::LENGTH_ERROR
        )))
    );
    assert_eq!(
        check(3),
        Err(MotifError::Map(MapEnd::Throw(
            super::super::prosite::BAD_ALLOC
        )))
    );
    assert_eq!(
        check(4),
        Err(MotifError::Map(MapEnd::HostDependent((1 << 30) + 1)))
    );
    assert_eq!(
        check(5),
        Err(MotifError::Map(MapEnd::Throw(
            super::super::prosite::BAD_ALLOC
        )))
    );
    assert_eq!(
        check(6),
        Err(MotifError::Map(MapEnd::HostDependent(1610612735)))
    );
    assert_eq!(check(7), Err(MotifError::Match(MatchError::Exhausted)));
}

// What NCBI's rpsbproc 0.5.1 allocated for edited motifs (`tests/fixtures/rpsbproc_growth/
// binary_allocs.tsv`; README there): the allocations of its run with the motif minus those of a
// run with a motif of the same shape and small counts, as a signed count per size in bytes. The
// model's allocations of the string and the vector of the motif check (`growth::events`) give the
// same difference. `cat` is the catalytic site of `mfirst.asn` (3 residues, `HKY`, of a query of
// 185), `dim` a site of 24 mapped residues whose positions are the first motif site's.
const BINARY_ROWS: usize = 26;
const BINARY_ALLOCS: &str = include_str!("../../tests/fixtures/rpsbproc_growth/binary_allocs.tsv");

#[test]
fn the_motif_checks_of_the_binary_allocate_as_the_model_does() {
    use super::super::prosite::growth_events as events;
    use std::collections::BTreeMap;
    let allocations = |site: &str, motif: &str| -> BTreeMap<u64, i64> {
        let (positions, mapped, residues) = if site == "cat" {
            let mut residues = vec![b'A'; 185];
            residues[61..64].copy_from_slice(b"HKY");
            (
                "10,11,12".to_string(),
                mapped(&[(61, 10), (62, 11), (63, 12)]),
                residues,
            )
        } else {
            let pairs: Vec<(i64, i64)> = (0..24).map(|i| (20 + i, 100 + i)).collect();
            (
                (100..124)
                    .map(|p| p.to_string())
                    .collect::<Vec<_>>()
                    .join(","),
                mapped(&pairs),
                vec![b'A'; 185],
            )
        };
        let data = load(&[("A", motif, &positions)]);
        events::take();
        let _ = MotifState::default().motif_check(&data, site_of(&data), &mapped, &residues);
        let mut counts = BTreeMap::new();
        for event in events::take() {
            *counts
                .entry(event[1..].parse::<u64>().unwrap())
                .or_insert(0i64) += 1;
        }
        counts
    };
    let (mut checked, mut wrong) = (0, Vec::new());
    for line in BINARY_ALLOCS
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
    {
        let fields: Vec<&str> = line.split('\t').collect();
        let (site, motif, base) = (fields[0], fields[1], fields[2]);
        let delta = fields.get(3).copied().unwrap_or("");
        let (target, baseline) = (allocations(site, motif), allocations(site, base));
        let mut difference = BTreeMap::new();
        for size in target.keys().chain(baseline.keys()) {
            let d =
                target.get(size).copied().unwrap_or(0) - baseline.get(size).copied().unwrap_or(0);
            if d != 0 {
                difference.insert(*size, d);
            }
        }
        let text = difference
            .iter()
            .map(|(size, d)| format!("{size}:{d:+}"))
            .collect::<Vec<_>>()
            .join(" ");
        if text != delta {
            wrong.push(format!(
                "{site} {motif} (baseline {base}): model {text:?}, binary {delta:?}"
            ));
        }
        checked += 1;
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
    assert_eq!(checked, BINARY_ROWS);
}
