use super::*;

const MG1655: usize = 4_641_652;
const PO157: usize = 92_721;
const EDL933: usize = 5_528_445;

// Ŝ(p) is the row of the kind at the largest probed thread count not above p.
#[test]
fn speedup_takes_the_largest_probed_thread_count_not_above_p() {
    assert_eq!(speedup_tenths(Kind::Disc, 1), 10);
    assert_eq!(speedup_tenths(Kind::Disc, 2), 17);
    assert_eq!(speedup_tenths(Kind::Disc, 3), 17);
    assert_eq!(speedup_tenths(Kind::Blastn, 4), 29);
    assert_eq!(speedup_tenths(Kind::Blastn, 15), 41);
    assert_eq!(speedup_tenths(Kind::Other, 16), 27);
    assert_eq!(speedup_tenths(Kind::Other, 128), 27);
    assert_eq!(Kind::of(true, "dc-megablast"), Kind::Disc);
    assert_eq!(Kind::of(false, "blastn"), Kind::Blastn);
    assert_eq!(Kind::of(false, "rmblastn"), Kind::Blastn);
    assert_eq!(Kind::of(false, "megablast"), Kind::Other);
    assert_eq!(Kind::of(false, "blastn-short"), Kind::Other);
}

// The cells of the probe (the switch is on by default in a test run, as in the strict set,
// unless the test runs with LOSAT_X_BNFEWSUBJ=0).
#[test]
fn the_rule_follows_the_probe() {
    let subjects =
        |lengths: &[usize], threads, kind| serial_subjects(lengths.iter().copied(), threads, kind);
    let chunks =
        |lengths: &[usize], threads, kind| serial_chunks(lengths.iter().copied(), threads, kind);
    if !on() {
        assert!(!subjects(&[MG1655, PO157], 16, Kind::Disc));
        assert!(!chunks(
            &[MAX_DBSEQ_LEN, EDL933 - MAX_DBSEQ_LEN],
            16,
            Kind::Disc
        ));
        return;
    }
    // One subject, or one thread: never.
    assert!(!subjects(&[MG1655], 16, Kind::Disc));
    assert!(!chunks(&[MAX_DBSEQ_LEN], 16, Kind::Disc));
    assert!(!subjects(&[MG1655, PO157], 1, Kind::Disc));
    // MG1655 + pO157 and EDL933's two chunks: every kind from two threads.
    for kind in [Kind::Disc, Kind::Blastn, Kind::Other] {
        assert!(subjects(&[MG1655, PO157], 2, kind));
        assert!(chunks(&[MAX_DBSEQ_LEN, EDL933 - MAX_DBSEQ_LEN], 2, kind));
    }
    // MG1655 + EDL933 (longest unit a 5 Mb chunk): from 4 threads, megablast from 8.
    assert!(!subjects(&[MG1655, EDL933], 3, Kind::Disc));
    assert!(subjects(&[MG1655, EDL933], 4, Kind::Disc));
    assert!(subjects(&[MG1655, EDL933], 4, Kind::Blastn));
    assert!(!subjects(&[MG1655, EDL933], 7, Kind::Other));
    assert!(subjects(&[MG1655, EDL933], 8, Kind::Other));
    // MG1655 in 3 parts: dc-megablast and blastn from 8 threads, megablast never.
    let third = [MG1655 / 3 + 1, MG1655 / 3, MG1655 / 3];
    assert!(!subjects(&third, 7, Kind::Disc));
    assert!(subjects(&third, 8, Kind::Disc));
    assert!(subjects(&third, 8, Kind::Blastn));
    assert!(!subjects(&third, 64, Kind::Other));
    // MG1655 in 10 parts, a 10 Mb subject's two equal chunks at two threads: never.
    assert!(!subjects(&[MG1655 / 10; 10], 64, Kind::Blastn));
    assert!(!chunks(&[MAX_DBSEQ_LEN, MAX_DBSEQ_LEN], 2, Kind::Disc));
    assert!(chunks(&[MAX_DBSEQ_LEN, MAX_DBSEQ_LEN], 4, Kind::Disc));
    // Letters past u64 do not wrap.
    assert!(!subjects(&[usize::MAX, usize::MAX], 16, Kind::Disc));
}
