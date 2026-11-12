use super::*;

const BLOSUM62_UNGAPPED_LAMBDA: f64 = 0.3176;
const BLOSUM62_11_1_LAMBDA: f64 = 0.267;

// blast_parameters.c:218-221: the defaults give NCBI's raw X-drops (7 bits: 16 at the ideal
// BLOSUM62 Lambda); a product past Int4 is INT_MIN (x86-64 cvttsd2si), one past it per context.
#[test]
fn the_ungapped_x_drop_is_ncbis_ceiling_and_conversion() {
    assert_eq!(ungapped_x_dropoff(7.0, BLOSUM62_UNGAPPED_LAMBDA, 1.0), 16);
    assert_eq!(ungapped_x_dropoff(1e-300, BLOSUM62_UNGAPPED_LAMBDA, 1.0), 1);
    assert_eq!(ungapped_x_dropoff(0.5, BLOSUM62_UNGAPPED_LAMBDA, 1.0), 2);
    // An RPS database's scale multiplies the ceiling.
    assert_eq!(
        ungapped_x_dropoff(7.0, BLOSUM62_UNGAPPED_LAMBDA / 100.0, 100.0),
        152_800
    );
    // 1e9 bits: 2.18e9 at Lambda 0.3176 leaves Int4, 1.99e9 at 0.348 does not.
    assert_eq!(
        ungapped_x_dropoff(1e9, BLOSUM62_UNGAPPED_LAMBDA, 1.0),
        i32::MIN
    );
    assert_eq!(ungapped_x_dropoff(1e9, 0.348, 1.0), 1_991_802_243);
    for bits in [1e10, 1e300] {
        assert_eq!(
            ungapped_x_dropoff(bits, BLOSUM62_UNGAPPED_LAMBDA, 1.0),
            i32::MIN
        );
    }
}

// blast_parameters.c:454-463: the final X-drop is the larger of its double and the preliminary
// raw X-drop, then converted; INT_MIN past Int4.
#[test]
fn the_gapped_x_drops_are_ncbis() {
    assert_eq!(
        gapped_x_dropoffs(15.0, 25.0, BLOSUM62_11_1_LAMBDA),
        (38, 64)
    );
    // A final X-drop below the preliminary one takes the preliminary one.
    assert_eq!(
        gapped_x_dropoffs(30.0, 10.0, BLOSUM62_11_1_LAMBDA),
        (77, 77)
    );
    assert_eq!(gapped_x_dropoffs(0.0, 0.0, BLOSUM62_11_1_LAMBDA), (0, 0));
    assert_eq!(
        gapped_x_dropoffs(-1.0, -1.0, BLOSUM62_11_1_LAMBDA),
        (-2, -2)
    );
    // A preliminary X-drop past Int4: INT_MIN, and the final X-drop its own.
    assert_eq!(
        gapped_x_dropoffs(1e9, 25.0, BLOSUM62_11_1_LAMBDA),
        (i32::MIN, 64)
    );
    assert_eq!(
        gapped_x_dropoffs(15.0, 1e300, BLOSUM62_11_1_LAMBDA),
        (38, i32::MIN)
    );
    // The MAX is of doubles: a final double past Int4 stays past it.
    assert_eq!(
        gapped_x_dropoffs(1e9, 1e9, BLOSUM62_11_1_LAMBDA),
        (i32::MIN, i32::MIN)
    );
}

// blast_kappa.c:2374-2382: the redo's X-drop at the scaled Lambda, at least the search's final
// raw X-drop; INT_MIN past Int4 (`-xdrop_gap_final 1e8` at the scale 32).
#[test]
fn the_redo_x_drop_is_ncbis() {
    let scaled = BLOSUM62_11_1_LAMBDA / 32.0;
    assert_eq!(redo_x_dropoff(25.0, scaled, 64), 2076);
    assert_eq!(redo_x_dropoff(1.0, scaled, 64), 83);
    assert_eq!(redo_x_dropoff(0.0, scaled, 64), 64);
    assert_eq!(redo_x_dropoff(1e7, scaled, i32::MAX), i32::MAX);
    assert_eq!(redo_x_dropoff(1e8, scaled, 259_604_000), i32::MIN);
}

// blast_gapalign.c:457-464, 641-644: the quotient by the gap extension stays below 2^30 - 103; an
// X-drop below the gap costs is raised to them first (INT_MIN too).
#[test]
fn the_kernels_size_x_drops_below_their_limit() {
    let limit = KERNEL_X_DROP_QUOTIENT_LIMIT as i32;
    assert!(kernel_x_drop_fits(limit - 1, 11, 1));
    assert!(!kernel_x_drop_fits(limit, 11, 1));
    assert!(kernel_x_drop_fits(i32::MAX, 11 * 32, 32));
    assert!(kernel_x_drop_fits(i32::MIN, 11, 1));
    assert!(kernel_x_drop_fits(i32::MAX, 11, 0));
    assert!(check_gapped_x_dropoffs("BLASTP", (15.0, 25.0), (38, 64), (11, 1), true).is_ok());
    let error = check_gapped_x_dropoffs("BLASTP", (5e8, 25.0), (limit, limit), (11, 1), false)
        .unwrap_err()
        .to_string();
    assert!(
        error.starts_with("an -xdrop_gap of 500000000.0 bits"),
        "{error}"
    );
    // The final X-drop reaches the unscaled kernel only without the composition-based redo.
    assert!(check_gapped_x_dropoffs("BLASTX", (15.0, 5e8), (38, limit), (11, 1), false).is_ok());
    let error = check_gapped_x_dropoffs("BLASTX", (15.0, 5e8), (38, limit), (11, 1), true)
        .unwrap_err()
        .to_string();
    assert!(
        error.starts_with("an -xdrop_gap_final of 500000000.0 bits")
            && error.ends_with("is not supported by LOSAT's BLASTX"),
        "{error}"
    );
}
