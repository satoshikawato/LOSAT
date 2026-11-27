use super::*;

fn blosum62() -> PhiKarlin {
    phi_karlin("BLOSUM62", 11, 1).unwrap()
}

#[test]
fn table_blocks_and_messages() {
    assert_eq!(
        blosum62(),
        PhiKarlin {
            lambda: 0.270,
            k: 0.047,
            param_c: 0.50,
            h: 1.0
        }
    );
    assert_eq!(phi_karlin("BLOSUM45", 16, 1).unwrap().k, 0.024);
    assert_eq!(phi_karlin("PAM30", 9, 1).unwrap().param_c, 0.30);
    assert_eq!(
        phi_karlin("BLOSUM62", 10, 2),
        Err("The combination 10 for gap opening cost and 2 for gap extension is not supported in PHI-BLAST with matrix BLOSUM62\n".to_string())
    );
    assert_eq!(
        phi_karlin("BLOSUM90", 10, 1),
        Err("Matrix BLOSUM90 not allowed in PHI-BLAST\n".to_string())
    );
}

#[test]
fn search_space_and_occurrences() {
    // rd.base.tk.db.0: paja (86 sequences, 66,738 letters), minimum length 8.
    assert_eq!(search_space(66_738, 86, 8), 66_050);
    let at = |offsets: &[i32]| -> Vec<Occurrence> {
        offsets
            .iter()
            .map(|&offset| Occurrence { offset, length: 8 })
            .collect()
    };
    assert_eq!(effective_occurrences(&at(&[]), 8), 0);
    assert_eq!(effective_occurrences(&at(&[19, 208]), 8), 2);
    // Within half the minimum length of the last one counted: not counted.
    assert_eq!(effective_occurrences(&at(&[10, 14, 15, 19]), 8), 2);
}

#[test]
fn cutoff_evalue_and_bits_match_ncbi_rows() {
    // rd.base.tk.db.0 (-evalue 10): every hit down to score 3 is reported.
    assert_eq!(cutoff_score(50.0, 66_050, &blosum62(), 2, 8.338e-5), 1);
    let high = cutoff_score(50.0, 1_000_000_000_000, &blosum62(), 2, 8.338e-5);
    assert!(high > 50 && high < 80, "{high}");
    // Its rows: score 312 = 116 bits, 2e-34; 12 = 3.6 bits, 1.2; 3: 5.6; 9: 2.1 (7 pattern hits
    // in paja).
    assert_eq!(format!("{:.0}", bit_score(312, &blosum62())), "116");
    assert_eq!(format!("{:.1}", bit_score(12, &blosum62())), "3.6");
    assert_eq!(format!("{:.1}", evalue(12, &blosum62(), 2, 7)), "1.2");
    assert_eq!(format!("{:.1}", evalue(3, &blosum62(), 2, 7)), "5.6");
    assert_eq!(format!("{:.1}", evalue(9, &blosum62(), 2, 7)), "2.1");
    assert_eq!(format!("{:.0e}", evalue(312, &blosum62(), 2, 7)), "2e-34");
}

#[test]
fn sum_statistics_link_with_the_table_lambda_and_a_log_k_of_zero() {
    // 22 bits with BLOSUM62 11/1's 0.270 and logK 0: 22 ln 2 / 0.270 = 56.48, at most the cutoff.
    assert_eq!(word_cutoff(22.0, &blosum62(), 100), 56);
    assert_eq!(word_cutoff(22.0, &blosum62(), 30), 30);
    // `-gap_trigger 1e300` leaves `Int4`: INT_MIN, as blastp's gap trigger.
    assert_eq!(word_cutoff(1e300, &blosum62(), 30), i32::MIN);
    assert_eq!(
        link_context(&blosum62(), 300, 8, 66_050),
        crate::core::link_hsps::LinkContext {
            lambda: 0.270,
            log_k: 0.0,
            query_length: 300,
            length_adjustment: 8,
            eff_searchsp: 66_050,
        }
    );
}
