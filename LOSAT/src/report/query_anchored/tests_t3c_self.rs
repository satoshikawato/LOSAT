//! Unit tests of TBLASTX's query-anchored section of a query record that is also a subject under
//! `-parse_deflines` (line RG T3c fix 1, `aln_mix`: one `CAlnMixSeq` of both rows with
//! `fForceTranslation`): each compares `write_section` with the bytes of an NCBI BLAST+ 2.17.0 probe
//! (`-num_threads 1`; under `/home/kawato/losat-baselines/remaining-gaps-20261112/T/t3/`: the review's
//! `review_t3c/ws/in/`, `fix1_t3c/rv_ncbi/out/`, `fix1_t3c/c6/`, and the probes `in3c/`, `out/`). The HSPs
//! are NCBI's own (`-outfmt "6 qseqid sseqid qstart qend sstart send qframe sframe btop evalue"`), made
//! into Dense-segs by the report's owner (`seqalign::translated_row_segments`). The query's masks are
//! the nucleotides of the residues NCBI shows in lower case on each group's query row, with the group's
//! frame. The last two tests are the sections on which NCBI's merge throws (line RG T3c fix 2). Generated
//! by `T/t3/scripts_t3c/gen_unit_tests_t3c_fix2.py`.

use super::tests_t3b::{hit, ops_of_btop};
use super::*;
use crate::algorithm::tblastx::report::DisplayCode;
use crate::ncbi_exception::NcbiException;
use crate::report::serial::seqalign::Segments;

/// The probes' sequences by name.
const SEQUENCES: &[(&str, &str)] = &[
    ("mg0_q/mg0", "TTTGCCCCACCGATGAAAGGTACGTCGGCGTCGGGCTTCTATTCTGAGGCGTTCACCTTTGGAATGGAGGGTATCTATATACGCTGTAAGATGGCTATGGGCACAATGAACCCCCCGATGTTGGGCAACGTTAGTAAGAGGGATTCGAAAGGATTTGCCCCACCGATGAAAGGTACGTCGGCGTCGGGCTTCTATTCTGAGGCGTTCACCTTTGGAATGGAGGGTATCTATATACGCTGTAAGATGGCTATGGGCACAATGAACCCCCCGATGTTGGGCAACGTTAGTAAGAGGGATTCGAAA"),
    ("mg0_sx/mg0", "TTTGCCCCACCGATGAAAGGTACGTCGGCGTCGGGCTTCTATTCTGAGGCGTTCACCTTTGGAATGGAGGGTATCTATATACGCTGTAAGATGGCTATGGGCACAATGAACCCCCCGATGTTGGGCAACGTTAGTAAGAGGGATTCGAAAGGATTTGCCCCACCGATGAAAGGTACGTCGGCGTCGGGCTTCTATTCTGAGGCGTTCACCTTTGGAATGGAGGGTATCTATATACGCTGTAAGATGGCTATGGGCACAATGAACCCCCCGATGTTGGGCAACGTTAGTAAGAGGGATTCGAAA"),
    ("mg0_sx/mx0", "CAACATCCCGGCGACTTAGGCGACCCACAGAATCGTCCCTTCTAACGTAGTTCGCATAGTTCCCGTCCGTAGCCGGACTATTCGAACACCCAGTATTCGATTAACTCGGGCTTGACGTATTAGAGGCGTTAGTGTGCCAGGTAAGATACGCCATTTGCCCCACCGATGAAAGGTACGTCGGCGTCGGGCTTCTATTCTGAGGCGTTCACCTTTGGAATGGAGGGTATCTATATACGCTGTAAGATGGCTATGGGCACAATGAACCCCCCGATGTTGGGCAACGTTAGTAAGAGGGATTCGAAA"),
    ("rep2/rep2", "CATAACCTGGTATATGTCGGGTCTGCGAGGCCGAGTGTAGTATATTCCCTAGGAATATTCTTTGTTACGTGTAAGATTGACGGTAGTAGCACCGAGTTCCTTAGGTGGGAATCGTATGTCGTCCCCCCCTGCTTGTATAGACCTCTCCCTTTTGGAAGTGAGAAGGCTCGTAAACCTTTGGATCATTATTACCCATTAGCGCCGCGATTTTTCTCCCCAACCGTTATCTCGTATACGAAGCGGGGAACGCTAGGACTACCGGTTGAGTGCACTATTTGTACTACTCATTTCGCTCAATCCCATAACCTGGTATATGTCGGGTCTGCGAGGCCGAGTGTAGTATATTCCCTAGGAATATTCTTTGTTACGTGTAAGATTGACGGTAGTAGCACCGAGTTCCTTAGGTGGGAATCGTATGTCGTCCCCCCCTGCTTGTATAGACCTCTCCCTTTTGGAAGTGAGAAGGCTCGTAAACCTTTGGATCATTATTACCCATTAGCGCCGCGATTTTTCTCCCCAACCGTTATCTCGTATACGAAGCGGGGAACGCTAGGACTACCGGTTGAGTGCACTATTTGTACTACTCATTTCGCTCAATCC"),
    ("t_selfq/selfx", "AGTGACCACCAGTGGAAGAACCTACCCCTGCGACAACATATTGGCAATTGGATGACTATGGTATGGGAATGGGTGTGGTACAAGACTGACATCGCTGGGTTTCAAGAGCAAGCCCACTATCATCTTCCCGAAAAAGCTCTTCCGATAGCCCAACAGATCACTCGGAGCGTCCCAGCAAATTGGTCCAGGATCGAGAACGTCCACTCCGGTTCTTTCCGGAAGCACGCTAACCCAGAGTGCTGATTTTCGATTCTGCTCCAATTGGCAGGGACGCTCCGTGTGATCTGCTGTGCTATGGGTAGTGCCTTTTCGGGCAGGTGATAGTGTGCAGTCTCTTGAAAACCCCGGATATCCGTCTTGTACCAGACCCACTCCCATACCATGGTCATCCAGTTGCCAATGTGCATGCGCAAGTGTAGGTTCTTCCATTGGTGGCAGGACTCTGTGAGGTTG"),
    ("t_selfs/selfx", "AGTGACCACCAGTGGAAGAACCTACCCCTGCGACAACATATTGGCAATTGGATGACTATGGTATGGGAATGGGTGTGGTACAAGACTGACATCGCTGGGTTTCAAGAGCAAGCCCACTATCATCTTCCCGAAAAAGCTCTTCCGATAGCCCAACAGATCACTCGGAGCGTCCCAGCAAATTGGTCCAGGATCGAGAACGTCCACTCCGGTTCTTTCCGGAAGCACGCTAACCCAGAGTGCTGATTTTCGATTCTGCTCCAATTGGCAGGGACGCTCCGTGTGATCTGCTGTGCTATGGGTAGTGCCTTTTCGGGCAGGTGATAGTGTGCAGTCTCTTGAAAACCCCGGATATCCGTCTTGTACCAGACCCACTCCCATACCATGGTCATCCAGTTGCCAATGTGCATGCGCAAGTGTAGGTTCTTCCATTGGTGGCAGGACTCTGTGAGGTTG"),
    ("t_selfs/ydup", "AGTGACCACCAGTGGAAGAACCTACCCCTGCGACAACATATTGGCAATTGGATGACTATGGTATGGGAATGGGTGTGGTACAAGACTGACATCGCTGGGTTTCAAGAGCAAGCCCACTATCATCTTCCCGAAAAAGCTCTTCCGATAGCCCAACAGATCACTCGGAGCGTCCCAGCAAATTGGTCCAGGATCGAGAACGTCCACTCCGGTTCTTTCCGGAAGCACGCTAACCCAGAGTGCTGATTTTCGATTCTGCTCCAATTGGCAGGGACGCTCCGTGTGATCTGCTGTGCTATGGGTAGTGCCTTTTCGGGCAGGTGATAGTGTGCAGTCTCTTGAAAACCCCGGATATCCGTCTTGTACCAGACCCACTCCCATACCATGGTCATCCAGTTGCCAATGTGCATGCGCAAGTGTAGGTTCTTCCATTGGTGGCAGGACTCTGTGAGGTTGACCCTCAATATCAGGTAAGCAACACCATGGAGTGACCACCAGTGGAAGAACCTACCCCTGCGACAACATATTGGCAATTGGATGACTATGGTATGGGAATGGGTGTGGTACAAGACTGACATCGCTGGGTTTCAAGAGCAAGCCCACTATCATCTTCCCGAAAAAGCTCTTCCGATAGCCCAACAGATCACTCGGAGCGTCCCAGCAAATTGGTCCAGGATCGAGAACGTCCACTCCGGTTCTTTCCGGAAGCACGCTAACCCAGAGTGCTGATTTTCGATTCTGCTCCAATTGGCAGGGACGCTCCGTGTGATCTGCTGTGCTATGGGTAGTGCCTTTTCGGGCAGGTGATAGTGTGCAGTCTCTTGAAAACCCCGGATATCCGTCTTGTACCAGACCCACTCCCATACCATGGTCATCCAGTTGCCAATGTGCATGCGCAAGTGTAGGTTCTTCCATTGGTGGCAGGACTCTGTGAGGTTG"),
    ("t_selfs3/selfx", "AGTGACCACCAGTGGAAGAACCTACCCCTGCGACAACATATTGGCAATTGGATGACTATGGTATGGGAATGGGTGTGGTACAAGACTGACATCGCTGGGTTTCAAGAGCAAGCCCACTATCATCTTCCCGAAAAAGCTCTTCCGATAGCCCAACAGATCACTCGGAGCGTCCCAGCAAATTGGTCCAGGATCGAGAACGTCCACTCCGGTTCTTTCCGGAAGCACGCTAACCCAGAGTGCTGATTTTCGATTCTGCTCCAATTGGCAGGGACGCTCCGTGTGATCTGCTGTGCTATGGGTAGTGCCTTTTCGGGCAGGTGATAGTGTGCAGTCTCTTGAAAACCCCGGATATCCGTCTTGTACCAGACCCACTCCCATACCATGGTCATCCAGTTGCCAATGTGCATGCGCAAGTGTAGGTTCTTCCATTGGTGGCAGGACTCTGTGAGGTTG"),
    ("t_selfs3/ydup", "AGTGACCACCAGTGGAAGAACCTACCCCTGCGACAACATATTGGCAATTGGATGACTATGGTATGGGAATGGGTGTGGTACAAGACTGACATCGCTGGGTTTCAAGAGCAAGCCCACTATCATCTTCCCGAAAAAGCTCTTCCGATAGCCCAACAGATCACTCGGAGCGTCCCAGCAAATTGGTCCAGGATCGAGAACGTCCACTCCGGTTCTTTCCGGAAGCACGCTAACCCAGAGTGCTGATTTTCGATTCTGCTCCAATTGGCAGGGACGCTCCGTGTGATCTGCTGTGCTATGGGTAGTGCCTTTTCGGGCAGGTGATAGTGTGCAGTCTCTTGAAAACCCCGGATATCCGTCTTGTACCAGACCCACTCCCATACCATGGTCATCCAGTTGCCAATGTGCATGCGCAAGTGTAGGTTCTTCCATTGGTGGCAGGACTCTGTGAGGTTGATACGGATACGTGGTTTACGGGATCCAGTCAGTGACCACCAGTGGAAGAACCTACCCCTGCGACAACATATTGGCAATTGGATGACTATGGTATGGGAATGGGTGTGGTACAAGACTGACATCGCTGGGTTTCAAGAGCAAGCCCACTATCATCTTCCCGAAAAAGCTCTTCCGATAGCCCAACAGATCACTCGGAGCGTCCCAGCAAATTGGTCCAGGATCGAGAACGTCCACTCCGGTTCTTTCCGGAAGCACGCTAACCCAGAGTGCTGATTTTCGATTCTGCTCCAATTGGCAGGGACGCTCCGTGTGATCTGCTGTGCTATGGGTAGTGCCTTTTCGGGCAGGTGATAGTGTGCAGTCTCTTGAAAACCCCGGATATCCGTCTTGTACCAGACCCACTCCCATACCATGGTCATCCAGTTGCCAATGTGCATGCGCAAGTGTAGGTTCTTCCATTGGTGGCAGGACTCTGTGAGGTTG"),
    ("txself/txself", "AAGACTGCGGCTCCTGAATGGTCGGCGCTATTAGTTTGTATAATTGCGCACGATATGCAAAAGTTAAATCATACCTGGTTCAGTCATGGGACCGAAATCAAGAAAGATGACGGCTACGAAACATCTCACGTAGCACCGGGATAGTTGACAGTTCGTAAGCATCGTCCTTCTTGATCTCTGTCCCGTGAGAGAACCAAGTATGGTTTAAGGCCTGCATGTCGTGAATAATAATACTGACCAGTGATGCAGCCCATTCAGGGGCAGCAGTTTTGTACGATGGTGTCCGTACCAGCGCACCCGAATGGTCGGCACTCCTCGTCTGCATTATCGCGCATGATATGCAGAAGCTTAATCATACTTGGTTCTCGCATGGACACGAGATAAAGAAGGACGACGCACCCCAGTGGTCCATCTATCCTAA"),
];

fn sequence(name: &str) -> &'static [u8] {
    SEQUENCES
        .iter()
        .find(|(known, _)| *known == name)
        .map(|(_, letters)| letters.as_bytes())
        .expect("a probe sequence")
}

/// A probe's sequences: the query and the subjects (the records' letters and labels), the query's
/// masks with their frames; genetic code 1 for both.
struct Probe {
    query: (&'static [u8], &'static str),
    subjects: Vec<(&'static [u8], &'static str)>,
    query_masks: Vec<Mask>,
}

impl Sequences for Probe {
    fn letters(&self, seq: SeqKey, from: u32, to: u32, out: &mut Vec<u8>) {
        let letters = match seq {
            SeqKey::Query => self.query.0,
            SeqKey::Subject(index) => self.subjects[index as usize].0,
        };
        out.extend_from_slice(&letters[from as usize..=to as usize]);
    }

    fn label(&self, seq: SeqKey, out: &mut Vec<u8>) {
        out.extend_from_slice(
            match seq {
                SeqKey::Query => self.query.1,
                SeqKey::Subject(index) => self.subjects[index as usize].1,
            }
            .as_bytes(),
        );
    }

    fn masks(&self, seq: SeqKey, out: &mut Vec<Mask>) {
        if seq == SeqKey::Query {
            out.extend_from_slice(&self.query_masks);
        }
    }

    fn genetic_code(&self, _master: bool) -> Option<DisplayCode> {
        Some(DisplayCode::new(
            &crate::utils::genetic_code::GeneticCode::try_from_id(1).expect("genetic code 1"),
        ))
    }
}

/// The section of the HSPs `rows` (subject index, qstart, qend, sstart, send, the query's and the
/// subject's frames in their records, BTOP), subject `own` being the query's own record, in format
/// `format`.
fn section(
    probe: &Probe,
    own: u32,
    rows: &[(u32, u32, u32, u32, u32, i8, i8, &str)],
    format: u8,
) -> String {
    let mut hsps = Hsps::of_translated_nucleotides();
    hsps.set_query_bioseq(own);
    let mut segments = Segments::default();
    for &(subject, qstart, qend, sstart, send, qframe, sframe, btop) in rows {
        let hit = hit(
            qstart,
            qend,
            sstart,
            send,
            i32::from(qframe),
            ops_of_btop(btop, false),
        );
        let strands = crate::report::serial::seqalign::translated_row_segments(
            &hit,
            Some((qframe, probe.query.0.len())),
            Some((sframe, probe.subjects[subject as usize].0.len())),
            &mut segments,
        );
        hsps.push_dense_seg(subject, segments.starts(), segments.lens(), 0, 0, strands);
    }
    let view = View {
        identity: matches!(format, 1 | 3),
        anchored: matches!(format, 1 | 2),
    };
    let mut out = Vec::new();
    write_section(&mut out, view, 60, &hsps, probe).expect("the section");
    String::from_utf8(out).expect("ASCII")
}

// review `x.self.noseg.f1` (`-seg no`): the record's HSP with itself in each query frame folds into the query row (`this position already exists`), its HSPs off the diagonal get extra rows of the query's sequence shown as the record, before the other subject's rows.
#[test]
fn tblastx_self_record_folds_into_the_query_row() {
    let probe = Probe {
        query: (sequence("t_selfq/selfx"), "selfx"),
        subjects: vec![
            (sequence("t_selfs3/selfx"), "selfx"),
            (sequence("t_selfs3/ydup"), "ydup"),
        ],
        query_masks: vec![],
    };
    let got = section(
        &probe,
        0,
        &[
            (1, 452, 3, 935, 486, -2, -2, "150"),
            (1, 452, 3, 452, 3, -2, -2, "150"),
            (1, 1, 453, 1, 453, 1, 1, "151"),
            (1, 1, 453, 484, 936, 1, 1, "151"),
            (1, 453, 1, 936, 484, -1, -1, "151"),
            (1, 453, 1, 453, 1, -1, -1, "151"),
            (1, 2, 451, 2, 451, 2, 2, "150"),
            (1, 2, 451, 485, 934, 2, 2, "150"),
            (1, 3, 452, 3, 452, 3, 3, "150"),
            (1, 3, 452, 486, 935, 3, 3, "150"),
            (1, 451, 2, 934, 485, -3, -3, "150"),
            (1, 451, 2, 451, 2, -3, -3, "150"),
            (1, 1, 435, 923, 489, 1, -2, "1DC6PH2QM19AR4QT29VQ3SLFARCKFHRAKNE2CW*T4LD1LF8CW1MI1SR7CW1VC2KN1RAIM4QH5MI5MI1MC1KRCG6"),
            (1, 1, 435, 440, 6, 1, -2, "1DC6PH2QM19AR4QT29VQ3SLFARCKFHRAKNE2CW*T4LD1LF8CW1MI1SR7CW1VC2KN1RAIM4QH5MI5MI1MC1KRCG6"),
            (1, 440, 6, 1, 435, -2, 1, "1CD6HP2MQ19RA4TQ29QV3LSAFCRFKRHKAEN2WCT*4DL1FL8WC1IM1RS7WC1CV2NK1ARMI4HQ5IM5IM1CM1RKGC6"),
            (1, 440, 6, 484, 918, -2, 1, "1CD6HP2MQ19RA4TQ29QV3LSAFCRFKRHKAEN2WCT*4DL1FL8WC1IM1RS7WC1CV2NK1ARMI4HQ5IM5IM1CM1RKGC6"),
            (1, 6, 206, 918, 718, 3, -1, "2VM4PL1THTAYH5YH3MV1VL3*GHYRPWGVF2ADSCPT2SPSA2SGSTSTDH1PTTA2ST3SCKQ1VE1DN1EKRSPA1"),
            (1, 6, 206, 435, 235, 3, -1, "2VM4PL1THTAYH5YH3MV1VL3*GHYRPWGVF2ADSCPT2SPSA2SGSTSTDH1PTTA2ST3SCKQ1VE1DN1EKRSPA1"),
            (1, 243, 434, 681, 490, 3, -1, "IVFL1SPAG2GCRW3CS3CG1GR*KCSLF2QKVM2CGSL3TPPS1IVRS3DH1LF3GS2VI2VMHLAS1VG6"),
            (1, 243, 434, 198, 7, 3, -1, "IVFL1SPAG2GCRW3CS3CG1GR*KCSLF2QKVM2CGSL3TPPS1IVRS3DH1LF3GS2VI2VMHLAS1VG6"),
            (1, 435, 235, 6, 206, -1, 3, "2MV4LP1HTATHY5HY3VM1LV3G*YHPRGWFV2DACSTP2PSAS2GSTSTSHD1TPAT2TS3CSQK1EV1ND1KESRAP1"),
            (1, 435, 235, 489, 689, -1, 3, "2MV4LP1HTATHY5HY3VM1LV3G*YHPRGWFV2DACSTP2PSAS2GSTSTSHD1TPAT2TS3CSQK1EV1ND1KESRAP1"),
            (1, 198, 7, 243, 434, -1, 3, "VILF1PSGA2CGWR3SC3GC1RGK*SCFL2KQMV2GCLS3PTSP1VISR3HD1FL3SG2IV2MVLHSA1GV6"),
            (1, 198, 7, 726, 917, -1, 3, "VILF1PSGA2CGWR3SC3GC1RGK*SCFL2KQMV2GCLS3PTSP1VISR3HD1FL3SG2IV2MVLHSA1GV6"),
            (1, 445, 5, 479, 919, -3, 2, "2PVAT1NS4TP1ADCNTI2TI2PL3SN1SC3RLIT1GLVG2RSLKHP2TICF2RKHLYFPR1HPSN2HL3LQPQ2APEG1KRIT2LPGV*LRSAG1GTKLNTRQ1GARDSF1SFWCTS1LWLQ3EV2VAGV1SWEVEVLP2EADG2GVLQAS2TNQPRGCYQP3TR1IT3*W2NS2YCVAVC1GS4TI2"),
            (1, 8, 436, 916, 488, 2, -3, "1SN4PT1DANCIT2IT2LP3NS1CS3LRTI1LGGV2SRKLPH2ITFC2KRLHFYRP1PHNS2LH3QLQP2PAGE1RKTI2PLVGL*SRGA1TGLKTNQR1AGDRFS1FSCWST1WLQL3VE2AVVG1WSVEVEPL2AEGD2VGQLSA2NTPQGRYCPQ3RT1TI3W*2SN2CYAVCV1SG4IT2"),
            (1, 8, 436, 433, 5, 2, -3, "1SN4PT1DANCIT2IT2LP3NS1CS3LRTI1LGGV2SRKLPH2ITFC2KRLHFYRP1PHNS2LH3QLQP2PAGE1RKTI2PLVGL*SRGA1TGLKTNQR1AGDRFS1FSCWST1WLQL3VE2AVVG1WSVEVEPL2AEGD2VGQLSA2NTPQGRYCPQ3RT1TI3W*2SN2CYAVCV1SG4IT2"),
            (1, 433, 5, 8, 436, -3, 2, "1NS4TP1ADCNTI2TI2PL3SN1SC3RLIT1GLVG2RSLKHP2TICF2RKHLYFPR1HPSN2HL3LQPQ2APEG1KRIT2LPGV*LRSAG1GTKLNTRQ1GARDSF1SFWCTS1LWLQ3EV2VAGV1SWEVEVLP2EADG2GVLQAS2TNQPRGCYQP3TR1IT3*W2NS2YCVAVC1GS4TI2"),
            (1, 356, 394, 892, 854, 2, -3, "SA2RLPA1PGI*2SYSG1"),
            (1, 356, 394, 409, 371, 2, -3, "SA2RLPA1PGI*2SYSG1"),
            (1, 409, 371, 356, 394, -3, 2, "AS2LRAP1GP*I2YSGS1"),
            (1, 409, 371, 839, 877, -3, 2, "AS2LRAP1GP*I2YSGS1"),
            (1, 374, 412, 368, 406, 2, 2, "1IT1WISPSW1CSQS1AQ2"),
            (1, 368, 406, 374, 412, 2, 2, "1TI1IWPSWS1SCSQ1QA2"),
            (1, 374, 412, 851, 889, 2, 2, "1IT1WISPSW1CSQS1AQ2"),
            (1, 368, 406, 857, 895, 2, 2, "1TI1IWPSWS1SCSQ1QA2"),
            (0, 452, 3, 452, 3, -2, -2, "150"),
            (0, 1, 453, 1, 453, 1, 1, "151"),
            (0, 453, 1, 453, 1, -1, -1, "151"),
            (0, 2, 451, 2, 451, 2, 2, "150"),
            (0, 3, 452, 3, 452, 3, 3, "150"),
            (0, 451, 2, 451, 2, -3, -3, "150"),
            (0, 1, 435, 440, 6, 1, -2, "1DC6PH2QM19AR4QT29VQ3SLFARCKFHRAKNE2CW*T4LD1LF8CW1MI1SR7CW1VC2KN1RAIM4QH5MI5MI1MC1KRCG6"),
            (0, 440, 6, 1, 435, -2, 1, "1CD6HP2MQ19RA4TQ29QV3LSAFCRFKRHKAEN2WCT*4DL1FL8WC1IM1RS7WC1CV2NK1ARMI4HQ5IM5IM1CM1RKGC6"),
            (0, 6, 206, 435, 235, 3, -1, "2VM4PL1THTAYH5YH3MV1VL3*GHYRPWGVF2ADSCPT2SPSA2SGSTSTDH1PTTA2ST3SCKQ1VE1DN1EKRSPA1"),
            (0, 243, 434, 198, 7, 3, -1, "IVFL1SPAG2GCRW3CS3CG1GR*KCSLF2QKVM2CGSL3TPPS1IVRS3DH1LF3GS2VI2VMHLAS1VG6"),
            (0, 435, 235, 6, 206, -1, 3, "2MV4LP1HTATHY5HY3VM1LV3G*YHPRGWFV2DACSTP2PSAS2GSTSTSHD1TPAT2TS3CSQK1EV1ND1KESRAP1"),
            (0, 198, 7, 243, 434, -1, 3, "VILF1PSGA2CGWR3SC3GC1RGK*SCFL2KQMV2GCLS3PTSP1VISR3HD1FL3SG2IV2MVLHSA1GV6"),
            (0, 8, 436, 433, 5, 2, -3, "1SN4PT1DANCIT2IT2LP3NS1CS3LRTI1LGGV2SRKLPH2ITFC2KRLHFYRP1PHNS2LH3QLQP2PAGE1RKTI2PLVGL*SRGA1TGLKTNQR1AGDRFS1FSCWST1WLQL3VE2AVVG1WSVEVEPL2AEGD2VGQLSA2NTPQGRYCPQ3RT1TI3W*2SN2CYAVCV1SG4IT2"),
            (0, 433, 5, 8, 436, -3, 2, "1NS4TP1ADCNTI2TI2PL3SN1SC3RLIT1GLVG2RSLKHP2TICF2RKHLYFPR1HPSN2HL3LQPQ2APEG1KRIT2LPGV*LRSAG1GTKLNTRQ1GARDSF1SFWCTS1LWLQ3EV2VAGV1SWEVEVLP2EADG2GVLQAS2TNQPRGCYQP3TR1IT3*W2NS2YCVAVC1GS4TI2"),
            (0, 356, 394, 409, 371, 2, -3, "SA2RLPA1PGI*2SYSG1"),
            (0, 409, 371, 356, 394, -3, 2, "AS2LRAP1GP*I2YSGS1"),
            (0, 374, 412, 368, 406, 2, 2, "1IT1WISPSW1CSQS1AQ2"),
            (0, 368, 406, 374, 412, 2, 2, "1TI1IWPSWS1SCSQ1QA2"),
        ],
        1,
    );
    assert_eq!(
        got,
        r#"
For reading frame +1 of query sequence:

selfx  1    SDHQWKNLPLRQHIGNWMTMVWEWVWYKTDIAGFQEQAHYHLPEKALPIAQQITRSVPAN  180
ydup   1    ............................................................  180
ydup   484  ............................................................  663

selfx  181  WSRIENVHSGSFRKHANPEC*FSILLQLAGTLRVICCAMGSAFSGR**CAVS*KPRISVL  360
ydup   181  ....................*.........................**....*.......  360
ydup   664  ....................*.........................**....*.......  843

selfx  361  YQTHSHTMVIQLPMCMRKCRFFHWWQDSVRL  453
ydup   361  ...............................  453
ydup   844  ...............................  936

For reading frame +2 of query sequence:

selfx  2    VTTSGRTYPCDNILAIG*LWYGNGCGTRLTSLGFKSKPTIIFPKKLFR*PNRSLGASQQI  181
ydup   2    .................*..............................*...........  181
ydup   485  .................*..............................*...........  664

selfx  182  GPGSRTSTPVLSGSTLTQSADFRFCSNWQGRSV*SAVLWVVPFRAGDSVQSLENPGYPSC  361
ydup   182  .................................*..........................  361
ydup   665  .................................*..........................  844

selfx  362  TRPTPIPWSSSCQCACASVGSSIGGRTL*G  451
selfx  368      .T.IPW.SS.Q..               406
selfx  374    .I.WSS.CQ.A..                 412
ydup   362  ............................*.  451
ydup   845  ............................*.  934
ydup   368      .T.IPW.SS.Q..               406
ydup   374    .I.WSS.CQ.A..                 412
ydup   851      .T.IPW.SS.Q..               889
ydup   857    .I.WSS.CQ.A..                 895

For reading frame +3 of query sequence:

selfx  3    *PPVEEPTPATTYWQLDDYGMGMGVVQD*HRWVSRASPLSSSRKSSSDSPTDHSERPSKL  182
ydup   3    *...........................*...............................  182
ydup   486  *...........................*...............................  665

selfx  183  VQDRERPLRFFPEAR*PRVLIFDSAPIGRDAPCDLLCYG*CLFGQVIVCSLLKTPDIRLV  362
ydup   183  ...............*.......................*....................  362
ydup   666  ...............*.......................*....................  845

selfx  363  PDPLPYHGHPVANVHAQV*VLPLVAGLCEV  452
ydup   363  ..................*...........  452
ydup   846  ..................*...........  935

For reading frame -1 of query sequence:

selfx  453  QPHRVLPPMEEPTLAHAHWQLDDHGMGVGLVQDGYPGFSRDCTLSPARKGTTHSTADHTE  274
ydup   936  ............................................................  757
ydup   453  ............................................................  274

selfx  273  RPCQLEQNRKSALWVSVLPERTGVDVLDPGPICWDAPSDLLGYRKSFFGKMIVGLLLKPS  94
ydup   756  ............................................................  577
ydup   273  ............................................................  94

selfx  93   DVSLVPHPFPYHSHPIANMLSQG*VLPLVVT  1
ydup   576  .......................*.......  484
ydup   93   .......................*.......  1

For reading frame -2 of query sequence:

selfx  452  NLTESCHQWKNLHLRMHIGNWMTMVWEWVWYKTDIRGFQETAHYHLPEKALPIAQQITRS  273
ydup   935  ............................................................  756
ydup   452  ............................................................  273

selfx  272  VPANWSRIENQHSGLACFRKEPEWTFSILDQFAGTLRVICWAIGRAFSGR**WACS*NPA  93
ydup   755  ..................................................**....*...  576
ydup   272  ..................................................**....*...  93

selfx  92   MSVLYHTHSHTIVIQLPICCRRGRFFHWWS  3
ydup   575  ..............................  486
ydup   92   ..............................  3

For reading frame -3 of query sequence:

selfx  451  TSQSPATNGRTYTCACTLATG*PWYGSGSGTRRISGVFKRLHTITCPKRHYP*HSRSHGA  272
ydup   934  .....................*..............................*.......  755
ydup   451  .....................*..............................*.......  272

selfx  271  SLPIGAESKISTLG*RASGKNRSGRSRSWTNLLGRSE*SVGLSEELFREDDSGLALETQR  92
ydup   754  ..............*......................*......................  575
ydup   271  ..............*......................*......................  92

selfx  91   CQSCTTPIPIP*SSNCQYVVAGVGSSTGGH  2
ydup   574  ...........*..................  485
ydup   91   ...........*..................  2

"#
    );
}

// review `x.self.noseg.mts1.f1` (`-max_target_seqs 1`): only the record's HSPs with itself, no single refseq in the mix; a query frame whose HSPs of other subject frames throw on the strand check.
#[test]
fn tblastx_self_record_alone() {
    let probe = Probe {
        query: (sequence("t_selfq/selfx"), "selfx"),
        subjects: vec![
            (sequence("t_selfs/ydup"), "ydup"),
            (sequence("t_selfs/selfx"), "selfx"),
        ],
        query_masks: vec![],
    };
    let got = section(
        &probe,
        1,
        &[
            (1, 452, 3, 452, 3, -2, -2, "150"),
            (1, 1, 453, 1, 453, 1, 1, "151"),
            (1, 453, 1, 453, 1, -1, -1, "151"),
            (1, 2, 451, 2, 451, 2, 2, "150"),
            (1, 3, 452, 3, 452, 3, 3, "150"),
            (1, 451, 2, 451, 2, -3, -3, "150"),
            (1, 1, 435, 440, 6, 1, -2, "1DC6PH2QM19AR4QT29VQ3SLFARCKFHRAKNE2CW*T4LD1LF8CW1MI1SR7CW1VC2KN1RAIM4QH5MI5MI1MC1KRCG6"),
            (1, 440, 6, 1, 435, -2, 1, "1CD6HP2MQ19RA4TQ29QV3LSAFCRFKRHKAEN2WCT*4DL1FL8WC1IM1RS7WC1CV2NK1ARMI4HQ5IM5IM1CM1RKGC6"),
            (1, 6, 206, 435, 235, 3, -1, "2VM4PL1THTAYH5YH3MV1VL3*GHYRPWGVF2ADSCPT2SPSA2SGSTSTDH1PTTA2ST3SCKQ1VE1DN1EKRSPA1"),
            (1, 243, 434, 198, 7, 3, -1, "IVFL1SPAG2GCRW3CS3CG1GR*KCSLF2QKVM2CGSL3TPPS1IVRS3DH1LF3GS2VI2VMHLAS1VG6"),
            (1, 435, 235, 6, 206, -1, 3, "2MV4LP1HTATHY5HY3VM1LV3G*YHPRGWFV2DACSTP2PSAS2GSTSTSHD1TPAT2TS3CSQK1EV1ND1KESRAP1"),
            (1, 198, 7, 243, 434, -1, 3, "VILF1PSGA2CGWR3SC3GC1RGK*SCFL2KQMV2GCLS3PTSP1VISR3HD1FL3SG2IV2MVLHSA1GV6"),
            (1, 8, 436, 433, 5, 2, -3, "1SN4PT1DANCIT2IT2LP3NS1CS3LRTI1LGGV2SRKLPH2ITFC2KRLHFYRP1PHNS2LH3QLQP2PAGE1RKTI2PLVGL*SRGA1TGLKTNQR1AGDRFS1FSCWST1WLQL3VE2AVVG1WSVEVEPL2AEGD2VGQLSA2NTPQGRYCPQ3RT1TI3W*2SN2CYAVCV1SG4IT2"),
            (1, 433, 5, 8, 436, -3, 2, "1NS4TP1ADCNTI2TI2PL3SN1SC3RLIT1GLVG2RSLKHP2TICF2RKHLYFPR1HPSN2HL3LQPQ2APEG1KRIT2LPGV*LRSAG1GTKLNTRQ1GARDSF1SFWCTS1LWLQ3EV2VAGV1SWEVEVLP2EADG2GVLQAS2TNQPRGCYQP3TR1IT3*W2NS2YCVAVC1GS4TI2"),
            (1, 356, 394, 409, 371, 2, -3, "SA2RLPA1PGI*2SYSG1"),
            (1, 409, 371, 356, 394, -3, 2, "AS2LRAP1GP*I2YSGS1"),
            (1, 374, 412, 368, 406, 2, 2, "1IT1WISPSW1CSQS1AQ2"),
            (1, 368, 406, 374, 412, 2, 2, "1TI1IWPSWS1SCSQ1QA2"),
        ],
        1,
    );
    assert_eq!(
        got,
        r#"
For reading frame +1 of query sequence:

selfx  1    SDHQWKNLPLRQHIGNWMTMVWEWVWYKTDIAGFQEQAHYHLPEKALPIAQQITRSVPAN  180

selfx  181  WSRIENVHSGSFRKHANPEC*FSILLQLAGTLRVICCAMGSAFSGR**CAVS*KPRISVL  360

selfx  361  YQTHSHTMVIQLPMCMRKCRFFHWWQDSVRL  453

For reading frame +2 of query sequence:

selfx  2    VTTSGRTYPCDNILAIG*LWYGNGCGTRLTSLGFKSKPTIIFPKKLFR*PNRSLGASQQI  181

selfx  182  GPGSRTSTPVLSGSTLTQSADFRFCSNWQGRSV*SAVLWVVPFRAGDSVQSLENPGYPSC  361

selfx  362  TRPTPIPWSSSCQCACASVGSSIGGRTL*G  451
selfx  368      .T.IPW.SS.Q..               406
selfx  374    .I.WSS.CQ.A..                 412

For reading frame +3 of query sequence:

selfx  3    *PPVEEPTPATTYWQLDDYGMGMGVVQD*HRWVSRASPLSSSRKSSSDSPTDHSERPSKL  182

selfx  183  VQDRERPLRFFPEAR*PRVLIFDSAPIGRDAPCDLLCYG*CLFGQVIVCSLLKTPDIRLV  362

selfx  363  PDPLPYHGHPVANVHAQV*VLPLVAGLCEV  452

For reading frame -1 of query sequence:

selfx  453  QPHRVLPPMEEPTLAHAHWQLDDHGMGVGLVQDGYPGFSRDCTLSPARKGTTHSTADHTE  274

selfx  273  RPCQLEQNRKSALWVSVLPERTGVDVLDPGPICWDAPSDLLGYRKSFFGKMIVGLLLKPS  94

selfx  93   DVSLVPHPFPYHSHPIANMLSQG*VLPLVVT  1

For reading frame -2 of query sequence:

selfx  452  NLTESCHQWKNLHLRMHIGNWMTMVWEWVWYKTDIRGFQETAHYHLPEKALPIAQQITRS  273

selfx  272  VPANWSRIENQHSGLACFRKEPEWTFSILDQFAGTLRVICWAIGRAFSGR**WACS*NPA  93

selfx  92   MSVLYHTHSHTIVIQLPICCRRGRFFHWWS  3

For reading frame -3 of query sequence:

selfx  451  TSQSPATNGRTYTCACTLATG*PWYGSGSGTRRISGVFKRLHTITCPKRHYP*HSRSHGA  272

selfx  271  SLPIGAESKISTLG*RASGKNRSGRSRSWTNLLGRSE*SVGLSEELFREDDSGLALETQR  92

selfx  91   CQSCTTPIPIP*SSNCQYVVAGVGSSTGGH  2

"#
    );
}

// review `t.self.s3.f1` (SEG on): the query's masks of each frame on the query row, none on the record's rows of width 3.
#[test]
fn tblastx_self_record_with_seg_masks() {
    let probe = Probe {
        query: (sequence("t_selfq/selfx"), "selfx"),
        subjects: vec![
            (sequence("t_selfs3/selfx"), "selfx"),
            (sequence("t_selfs3/ydup"), "ydup"),
        ],
        query_masks: vec![(101, 163, 3)],
    };
    let got = section(
        &probe,
        0,
        &[
            (1, 452, 3, 935, 486, -2, -2, "150"),
            (1, 452, 3, 452, 3, -2, -2, "150"),
            (1, 1, 453, 1, 453, 1, 1, "151"),
            (1, 1, 453, 484, 936, 1, 1, "151"),
            (1, 453, 1, 936, 484, -1, -1, "151"),
            (1, 453, 1, 453, 1, -1, -1, "151"),
            (1, 2, 451, 2, 451, 2, 2, "150"),
            (1, 2, 451, 485, 934, 2, 2, "150"),
            (1, 451, 2, 934, 485, -3, -3, "150"),
            (1, 451, 2, 451, 2, -3, -3, "150"),
            (1, 1, 435, 923, 489, 1, -2, "1DC6PH2QM19AR4QT29VQ3SLFARCKFHRAKNE2CW*T4LD1LF8CW1MI1SR7CW1VC2KN1RAIM4QH5MI5MI1MC1KRCG6"),
            (1, 1, 435, 440, 6, 1, -2, "1DC6PH2QM19AR4QT29VQ3SLFARCKFHRAKNE2CW*T4LD1LF8CW1MI1SR7CW1VC2KN1RAIM4QH5MI5MI1MC1KRCG6"),
            (1, 440, 6, 1, 435, -2, 1, "1CD6HP2MQ19RA4TQ29QV3LSAFCRFKRHKAEN2WCT*4DL1FL8WC1IM1RS7WC1CV2NK1ARMI4HQ5IM5IM1CM1RKGC6"),
            (1, 440, 6, 484, 918, -2, 1, "1CD6HP2MQ19RA4TQ29QV3LSAFCRFKRHKAEN2WCT*4DL1FL8WC1IM1RS7WC1CV2NK1ARMI4HQ5IM5IM1CM1RKGC6"),
            (1, 165, 452, 165, 452, 3, 3, "96"),
            (1, 165, 452, 648, 935, 3, 3, "96"),
            (1, 39, 101, 39, 101, 3, 3, "21"),
            (1, 39, 101, 522, 584, 3, 3, "21"),
            (1, 435, 235, 6, 206, -1, 3, "2MV4LP1HTATHY5HY3VM1LV3G*YHPRGWFV2DACSTP2PSAS2GSTSTSHD1TPAT2TS3CSQK1EV1ND1KESRAP1"),
            (1, 435, 235, 489, 689, -1, 3, "2MV4LP1HTATHY5HY3VM1LV3G*YHPRGWFV2DACSTP2PSAS2GSTSTSHD1TPAT2TS3CSQK1EV1ND1KESRAP1"),
            (1, 198, 7, 243, 434, -1, 3, "VILF1PSGA2CGWR3SC3GC1RGK*SCFL2KQMV2GCLS3PTSP1VISR3HD1FL3SG2IV2MVLHSA1GV6"),
            (1, 198, 7, 726, 917, -1, 3, "VILF1PSGA2CGWR3SC3GC1RGK*SCFL2KQMV2GCLS3PTSP1VISR3HD1FL3SG2IV2MVLHSA1GV6"),
            (1, 445, 5, 479, 919, -3, 2, "2PVAT1NS4TP1ADCNTI2TI2PL3SN1SC3RLIT1GLVG2RSLKHP2TICF2RKHLYFPR1HPSN2HL3LQPQ2APEG1KRIT2LPGV*LRSAG1GTKLNTRQ1GARDSF1SFWCTS1LWLQ3EV2VAGV1SWEVEVLP2EADG2GVLQAS2TNQPRGCYQP3TR1IT3*W2NS2YCVAVC1GS4TI2"),
            (1, 8, 436, 916, 488, 2, -3, "1SN4PT1DANCIT2IT2LP3NS1CS3LRTI1LGGV2SRKLPH2ITFC2KRLHFYRP1PHNS2LH3QLQP2PAGE1RKTI2PLVGL*SRGA1TGLKTNQR1AGDRFS1FSCWST1WLQL3VE2AVVG1WSVEVEPL2AEGD2VGQLSA2NTPQGRYCPQ3RT1TI3W*2SN2CYAVCV1SG4IT2"),
            (1, 8, 436, 433, 5, 2, -3, "1SN4PT1DANCIT2IT2LP3NS1CS3LRTI1LGGV2SRKLPH2ITFC2KRLHFYRP1PHNS2LH3QLQP2PAGE1RKTI2PLVGL*SRGA1TGLKTNQR1AGDRFS1FSCWST1WLQL3VE2AVVG1WSVEVEPL2AEGD2VGQLSA2NTPQGRYCPQ3RT1TI3W*2SN2CYAVCV1SG4IT2"),
            (1, 433, 5, 8, 436, -3, 2, "1NS4TP1ADCNTI2TI2PL3SN1SC3RLIT1GLVG2RSLKHP2TICF2RKHLYFPR1HPSN2HL3LQPQ2APEG1KRIT2LPGV*LRSAG1GTKLNTRQ1GARDSF1SFWCTS1LWLQ3EV2VAGV1SWEVEVLP2EADG2GVLQAS2TNQPRGCYQP3TR1IT3*W2NS2YCVAVC1GS4TI2"),
            (1, 243, 434, 681, 490, 3, -1, "IVFL1SPAG2GCRW3CS3CG1GR*KCSLF2QKVM2CGSL3TPPS1IVRS3DH1LF3GS2VI2VMHLAS1VG6"),
            (1, 243, 434, 198, 7, 3, -1, "IVFL1SPAG2GCRW3CS3CG1GR*KCSLF2QKVM2CGSL3TPPS1IVRS3DH1LF3GS2VI2VMHLAS1VG6"),
            (1, 39, 86, 885, 838, 3, -1, "YH5YH3MV1VL3"),
            (1, 39, 86, 402, 355, 3, -1, "YH5YH3MV1VL3"),
            (1, 356, 394, 892, 854, 2, -3, "SA2RLPA1PGI*2SYSG1"),
            (1, 356, 394, 409, 371, 2, -3, "SA2RLPA1PGI*2SYSG1"),
            (1, 409, 371, 356, 394, -3, 2, "AS2LRAP1GP*I2YSGS1"),
            (1, 409, 371, 839, 877, -3, 2, "AS2LRAP1GP*I2YSGS1"),
            (1, 243, 278, 681, 646, 3, -1, "IVFL1SPAG2GCRW3"),
            (1, 374, 412, 368, 406, 2, 2, "1IT1WISPSW1CSQS1AQ2"),
            (1, 368, 406, 374, 412, 2, 2, "1TI1IWPSWS1SCSQ1QA2"),
            (1, 243, 278, 198, 163, 3, -1, "IVFL1SPAG2GCRW3"),
            (1, 374, 412, 851, 889, 2, 2, "1IT1WISPSW1CSQS1AQ2"),
            (1, 368, 406, 857, 895, 2, 2, "1TI1IWPSWS1SCSQ1QA2"),
            (0, 452, 3, 452, 3, -2, -2, "150"),
            (0, 1, 453, 1, 453, 1, 1, "151"),
            (0, 453, 1, 453, 1, -1, -1, "151"),
            (0, 2, 451, 2, 451, 2, 2, "150"),
            (0, 451, 2, 451, 2, -3, -3, "150"),
            (0, 1, 435, 440, 6, 1, -2, "1DC6PH2QM19AR4QT29VQ3SLFARCKFHRAKNE2CW*T4LD1LF8CW1MI1SR7CW1VC2KN1RAIM4QH5MI5MI1MC1KRCG6"),
            (0, 440, 6, 1, 435, -2, 1, "1CD6HP2MQ19RA4TQ29QV3LSAFCRFKRHKAEN2WCT*4DL1FL8WC1IM1RS7WC1CV2NK1ARMI4HQ5IM5IM1CM1RKGC6"),
            (0, 165, 452, 165, 452, 3, 3, "96"),
            (0, 39, 101, 39, 101, 3, 3, "21"),
            (0, 435, 235, 6, 206, -1, 3, "2MV4LP1HTATHY5HY3VM1LV3G*YHPRGWFV2DACSTP2PSAS2GSTSTSHD1TPAT2TS3CSQK1EV1ND1KESRAP1"),
            (0, 198, 7, 243, 434, -1, 3, "VILF1PSGA2CGWR3SC3GC1RGK*SCFL2KQMV2GCLS3PTSP1VISR3HD1FL3SG2IV2MVLHSA1GV6"),
            (0, 8, 436, 433, 5, 2, -3, "1SN4PT1DANCIT2IT2LP3NS1CS3LRTI1LGGV2SRKLPH2ITFC2KRLHFYRP1PHNS2LH3QLQP2PAGE1RKTI2PLVGL*SRGA1TGLKTNQR1AGDRFS1FSCWST1WLQL3VE2AVVG1WSVEVEPL2AEGD2VGQLSA2NTPQGRYCPQ3RT1TI3W*2SN2CYAVCV1SG4IT2"),
            (0, 433, 5, 8, 436, -3, 2, "1NS4TP1ADCNTI2TI2PL3SN1SC3RLIT1GLVG2RSLKHP2TICF2RKHLYFPR1HPSN2HL3LQPQ2APEG1KRIT2LPGV*LRSAG1GTKLNTRQ1GARDSF1SFWCTS1LWLQ3EV2VAGV1SWEVEVLP2EADG2GVLQAS2TNQPRGCYQP3TR1IT3*W2NS2YCVAVC1GS4TI2"),
            (0, 243, 434, 198, 7, 3, -1, "IVFL1SPAG2GCRW3CS3CG1GR*KCSLF2QKVM2CGSL3TPPS1IVRS3DH1LF3GS2VI2VMHLAS1VG6"),
            (0, 39, 86, 402, 355, 3, -1, "YH5YH3MV1VL3"),
            (0, 356, 394, 409, 371, 2, -3, "SA2RLPA1PGI*2SYSG1"),
            (0, 409, 371, 356, 394, -3, 2, "AS2LRAP1GP*I2YSGS1"),
            (0, 243, 278, 198, 163, 3, -1, "IVFL1SPAG2GCRW3"),
            (0, 374, 412, 368, 406, 2, 2, "1IT1WISPSW1CSQS1AQ2"),
            (0, 368, 406, 374, 412, 2, 2, "1TI1IWPSWS1SCSQ1QA2"),
        ],
        1,
    );
    assert_eq!(
        got,
        r#"
For reading frame +1 of query sequence:

selfx  1    SDHQWKNLPLRQHIGNWMTMVWEWVWYKTDIAGFQEQAHYHLPEKALPIAQQITRSVPAN  180
ydup   1    ............................................................  180
ydup   484  ............................................................  663

selfx  181  WSRIENVHSGSFRKHANPEC*FSILLQLAGTLRVICCAMGSAFSGR**CAVS*KPRISVL  360
ydup   181  ....................*.........................**....*.......  360
ydup   664  ....................*.........................**....*.......  843

selfx  361  YQTHSHTMVIQLPMCMRKCRFFHWWQDSVRL  453
ydup   361  ...............................  453
ydup   844  ...............................  936

For reading frame +2 of query sequence:

selfx  2    VTTSGRTYPCDNILAIG*LWYGNGCGTRLTSLGFKSKPTIIFPKKLFR*PNRSLGASQQI  181
ydup   2    .................*..............................*...........  181
ydup   485  .................*..............................*...........  664

selfx  182  GPGSRTSTPVLSGSTLTQSADFRFCSNWQGRSV*SAVLWVVPFRAGDSVQSLENPGYPSC  361
ydup   182  .................................*..........................  361
ydup   665  .................................*..........................  844

selfx  362  TRPTPIPWSSSCQCACASVGSSIGGRTL*G  451
selfx  368      .T.IPW.SS.Q..               406
selfx  374    .I.WSS.CQ.A..                 412
ydup   362  ............................*.  451
ydup   845  ............................*.  934
ydup   368      .T.IPW.SS.Q..               406
ydup   374    .I.WSS.CQ.A..                 412
ydup   851      .T.IPW.SS.Q..               889
ydup   857    .I.WSS.CQ.A..                 895

For reading frame +3 of query sequence:

selfx  39   YWQLDDYGMGMGVVQD*HRWVsrasplsssrksssdsptdhsERPSKLVQDRERPLRFFP  218
selfx  39   ................*....                                         101
ydup   165                                            ..................  218
ydup   648                                            ..................  701
ydup   39   ................*....                                         101
ydup   522  ................*....                                         584

selfx  219  EAR*PRVLIFDSAPIGRDAPCDLLCYG*CLFGQVIVCSLLKTPDIRLVPDPLPYHGHPVA  398
ydup   219  ...*.......................*................................  398
ydup   702  ...*.......................*................................  881

selfx  399  NVHAQV*VLPLVAGLCEV  452
ydup   399  ......*...........  452
ydup   882  ......*...........  935

For reading frame -1 of query sequence:

selfx  453  QPHRVLPPMEEPTLAHAHWQLDDHGMGVGLVQDGYPGFSRDCTLSPARKGTTHSTADHTE  274
ydup   936  ............................................................  757
ydup   453  ............................................................  274

selfx  273  RPCQLEQNRKSALWVSVLPERTGVDVLDPGPICWDAPSDLLGYRKSFFGKMIVGLLLKPS  94
ydup   756  ............................................................  577
ydup   273  ............................................................  94

selfx  93   DVSLVPHPFPYHSHPIANMLSQG*VLPLVVT  1
ydup   576  .......................*.......  484
ydup   93   .......................*.......  1

For reading frame -2 of query sequence:

selfx  452  NLTESCHQWKNLHLRMHIGNWMTMVWEWVWYKTDIRGFQETAHYHLPEKALPIAQQITRS  273
ydup   935  ............................................................  756
ydup   452  ............................................................  273

selfx  272  VPANWSRIENQHSGLACFRKEPEWTFSILDQFAGTLRVICWAIGRAFSGR**WACS*NPA  93
ydup   755  ..................................................**....*...  576
ydup   272  ..................................................**....*...  93

selfx  92   MSVLYHTHSHTIVIQLPICCRRGRFFHWWS  3
ydup   575  ..............................  486
ydup   92   ..............................  3

For reading frame -3 of query sequence:

selfx  451  TSQSPATNGRTYTCACTLATG*PWYGSGSGTRRISGVFKRLHTITCPKRHYP*HSRSHGA  272
ydup   934  .....................*..............................*.......  755
ydup   451  .....................*..............................*.......  272

selfx  271  SLPIGAESKISTLG*RASGKNRSGRSRSWTNLLGRSE*SVGLSEELFREDDSGLALETQR  92
ydup   754  ..............*......................*......................  575
ydup   271  ..............*......................*......................  92

selfx  91   CQSCTTPIPIP*SSNCQYVVAGVGSSTGGH  2
ydup   574  ...........*..................  485
ydup   91   ...........*..................  2

"#
    );
}

// review `t.self.sl.f1` (`-subject_loc 1-200`): the record's HSPs with itself of the subject's interval.
#[test]
fn tblastx_self_record_subject_loc() {
    let probe = Probe {
        query: (sequence("t_selfq/selfx"), "selfx"),
        subjects: vec![
            (sequence("t_selfs3/selfx"), "selfx"),
            (sequence("t_selfs3/ydup"), "ydup"),
        ],
        query_masks: vec![(101, 163, 3)],
    };
    let got = section(
        &probe,
        0,
        &[
            (1, 1, 198, 1, 198, 1, 1, "66"),
            (1, 200, 3, 200, 3, -2, -2, "66"),
            (1, 198, 1, 198, 1, -1, -1, "66"),
            (1, 2, 199, 2, 199, 2, 2, "66"),
            (1, 440, 243, 1, 198, -2, 1, "1CD6HP2MQ19RA4TQ29"),
            (1, 199, 2, 199, 2, -3, -3, "66"),
            (
                1,
                244,
                435,
                197,
                6,
                1,
                -2,
                "4LD1LF8CW1MI1SR7CW1VC2KN1RAIM4QH5MI5MI1MC1KRCG6",
            ),
            (
                1,
                435,
                244,
                6,
                197,
                -1,
                3,
                "2MV4LP1HTATHY5HY3VM1LV3G*YHPRGWFV2DACSTP2PSAS2GSTSTSHD1TPAT2TS3CSQK1EV1ND1KE",
            ),
            (
                1,
                243,
                434,
                198,
                7,
                3,
                -1,
                "IVFL1SPAG2GCRW3CS3CG1GR*KCSLF2QKVM2CGSL3TPPS1IVRS3DH1LF3GS2VI2VMHLAS1VG6",
            ),
            (1, 39, 101, 39, 101, 3, 3, "21"),
            (1, 165, 200, 165, 200, 3, 3, "12"),
            (
                1,
                433,
                245,
                8,
                196,
                -3,
                2,
                "1NS4TP1ADCNTI2TI2PL3SN1SC3RLIT1GLVG2RSLKHP2TICF2RKHLYFPR1HPSN2HL3LQPQ2APEG1KR",
            ),
            (
                1,
                248,
                436,
                193,
                5,
                2,
                -3,
                "1FSCWST1WLQL3VE2AVVG1WSVEVEPL2AEGD2VGQLSA2NTPQGRYCPQ3RT1TI3W*2SN2CYAVCV1SG4IT2",
            ),
            (1, 98, 136, 97, 135, 2, 1, "2KQSEKQPATHIYIHFL1KE1"),
            (1, 97, 135, 98, 136, 1, 2, "2QKESQKAPHTYIHILF1EK1"),
            (1, 363, 377, 73, 59, 3, -3, "1DI1LI1"),
            (1, 183, 175, 44, 36, -1, -2, "3"),
            (1, 71, 39, 115, 147, -2, 1, "1SY1TLIPVEIKQA3"),
            (1, 115, 147, 71, 39, 1, -2, "1YS1LTPIEVKIAQ3"),
            (1, 326, 294, 71, 39, -2, -2, "1YS1LTPIEVKIAQ3"),
            (0, 1, 198, 1, 198, 1, 1, "66"),
            (0, 200, 3, 200, 3, -2, -2, "66"),
            (0, 198, 1, 198, 1, -1, -1, "66"),
            (0, 2, 199, 2, 199, 2, 2, "66"),
            (0, 440, 243, 1, 198, -2, 1, "1CD6HP2MQ19RA4TQ29"),
            (0, 199, 2, 199, 2, -3, -3, "66"),
            (
                0,
                244,
                435,
                197,
                6,
                1,
                -2,
                "4LD1LF8CW1MI1SR7CW1VC2KN1RAIM4QH5MI5MI1MC1KRCG6",
            ),
            (
                0,
                435,
                244,
                6,
                197,
                -1,
                3,
                "2MV4LP1HTATHY5HY3VM1LV3G*YHPRGWFV2DACSTP2PSAS2GSTSTSHD1TPAT2TS3CSQK1EV1ND1KE",
            ),
            (
                0,
                243,
                434,
                198,
                7,
                3,
                -1,
                "IVFL1SPAG2GCRW3CS3CG1GR*KCSLF2QKVM2CGSL3TPPS1IVRS3DH1LF3GS2VI2VMHLAS1VG6",
            ),
            (0, 39, 101, 39, 101, 3, 3, "21"),
            (0, 165, 200, 165, 200, 3, 3, "12"),
            (
                0,
                433,
                245,
                8,
                196,
                -3,
                2,
                "1NS4TP1ADCNTI2TI2PL3SN1SC3RLIT1GLVG2RSLKHP2TICF2RKHLYFPR1HPSN2HL3LQPQ2APEG1KR",
            ),
            (
                0,
                248,
                436,
                193,
                5,
                2,
                -3,
                "1FSCWST1WLQL3VE2AVVG1WSVEVEPL2AEGD2VGQLSA2NTPQGRYCPQ3RT1TI3W*2SN2CYAVCV1SG4IT2",
            ),
            (0, 98, 136, 97, 135, 2, 1, "2KQSEKQPATHIYIHFL1KE1"),
            (0, 97, 135, 98, 136, 1, 2, "2QKESQKAPHTYIHILF1EK1"),
            (0, 363, 377, 73, 59, 3, -3, "1DI1LI1"),
            (0, 183, 175, 44, 36, -1, -2, "3"),
            (0, 71, 39, 115, 147, -2, 1, "1SY1TLIPVEIKQA3"),
            (0, 115, 147, 71, 39, 1, -2, "1YS1LTPIEVKIAQ3"),
            (0, 326, 294, 71, 39, -2, -2, "1YS1LTPIEVKIAQ3"),
        ],
        1,
    );
    assert_eq!(
        got,
        r#"
For reading frame +1 of query sequence:

selfx  1    SDHQWKNLPLRQHIGNWMTMVWEWVWYKTDIAGFQEQAHYHLPEKALPIAQQITRSVPAN  180
selfx  98                                   ..KSKPTIIF.K.                 136
ydup   1    ............................................................  180
ydup   98                                   ..KSKPTIIF.K.                 136

selfx  181  WSRIEN  198
ydup   181  ......  198

For reading frame +2 of query sequence:

selfx  2    VTTSGRTYPCDNILAIG*LWYGNGCGTRLTSLGFKSKPTIIFPKKLFR*PNRSLGASQQI  181
selfx  97                                   ..QEQAHYHL.E.                 135
ydup   2    .................*..............................*...........  181
ydup   97                                   ..QEQAHYHL.E.                 135

selfx  182  GPGSRT  199
ydup   182  ......  199

For reading frame +3 of query sequence:

selfx  39   YWQLDDYGMGMGVVQD*HRWVsrasplsssrksssdsptdhsERPSKLVQDRERPLRFFP  218
selfx  165                                            ............        200

selfx  219  EAR*PRVLIFDSAPIGRDAPCDLLCYG*CLFGQVIVCSLLKTPDIRLVPDPLPYHGHPVA  398
ydup   198          VL.PG..CW...S...G.RKSF..KM..GL...PS.VS...H.F...S..I.  43
ydup   73                                                   .I.I.         59

selfx  399  NVHAQV*VLPLV  434
ydup   42   .MLS.G*.....  7

For reading frame -1 of query sequence:

selfx  198  VLDPGPICWDAPSDLLGYRKSFFGKMIVGLLLKPSDVSLVPHPFPYHSHPIANMLSQG*V  19
selfx  44        ...                                                      36
ydup   198  ..........................................................*.  19
ydup   44        ...                                                      36

selfx  18   LPLVVT  1
ydup   18   ......  1

For reading frame -2 of query sequence:

selfx  326  HYHLPEKALPIAQQITRSVPANWSRIENQHSGLACFRKEPEWTFSILDQFAGTLRVICWA  147
selfx  71   .S.TIVIQ...                                                   39
ydup   200                                            ..................  147
ydup   71   .S.TIVIQ...                                                   39

selfx  146  IGRAFSGR**WACS*NPAMSVLYHTHSHTIVIQLPICCRRGRFFHWWS  3
ydup   146  ........**....*.................................  3

For reading frame -3 of query sequence:

selfx  199  RSRSWTNLLGRSE*SVGLSEELFREDDSGLALETQRCQSCTTPIPIP*SSNCQYVVAGVG  20
ydup   199  .............*.................................*............  20

selfx  19   SSTGGH  2
ydup   19   ......  2

"#
    );
}

// line RG T3c probe `tx.t3c.self.f1`: a record with repeats on both strands against itself.
#[test]
fn tblastx_self_record_repeats_format_1() {
    let probe = Probe {
        query: (sequence("txself/txself"), "txself"),
        subjects: vec![(sequence("txself/txself"), "txself")],
        query_masks: vec![],
    };
    let got = section(
        &probe,
        0,
        &[
            (0, 420, 1, 420, 1, -2, -2, "140"),
            (0, 421, 2, 421, 2, -1, -1, "140"),
            (0, 419, 3, 419, 3, -3, -3, "139"),
            (0, 1, 420, 1, 420, 1, 1, "140"),
            (0, 3, 419, 3, 419, 3, 3, "139"),
            (0, 2, 421, 2, 421, 2, 2, "140"),
            (0, 138, 419, 418, 137, 3, -1, "1I*VM1SHSW*G4LF2VC1*C4W*FLKSAF1MI1*CIA1IMLQ1SRDS1AD4AL1LRYTDPGSVYRK1SA4SA1LS2CS2AI4KA9HT7PYQEWL1ITYI1"),
            (0, 418, 137, 138, 419, -1, 3, "1*IMV1HSWSG*4FL2CV1C*4*WLFSKFA1IM1C*AI1MIQL1RSSD1DA4LA1RLTYPDSGYVKR1AS4AS1SL2SC2IA4AK9TH7YPEQLW1TIIY1"),
            (0, 271, 140, 1, 132, -1, 1, "7AS1SL2SC2IA4AK16AG2LT1THIV"),
            (0, 1, 132, 271, 140, 1, -1, "7SA1LS2CS2AI4KA16GA2TL1HTVI"),
            (0, 120, 4, 152, 268, -2, 2, "2ASVIIVFL2FLGC1MVTR2GS1IV1LG2IV1RN1YN1ND*Q1RCRS3SGRS1"),
            (0, 152, 268, 120, 4, 2, -2, "2SAIVVILF2LFCG1VMRT2SG1VI1GL2VI1NR1NY1DNQ*1CRSR3GSSR1"),
            (0, 285, 410, 1, 126, 3, 1, "RK1SA27HT6AGPYQEWT1"),
            (0, 1, 126, 285, 410, 1, 3, "KR1AS27TH6GAYPEQTW1"),
            (0, 269, 3, 3, 269, -3, 3, "ND1CGPS2GVCGIATIGSQL1YNYCSA1HY1GKLV1PS3LQ1RW2DN2GRR*CR2TNVINSYRPSGTAGTI*V1VS2PA9LE3*W1NKFA1IM1C*AI2QL1NSSD1DA6"),
            (0, 3, 269, 269, 3, 3, -3, "DN1GCSP2VGGCAIITSGLQ1NYCYAS1YH1KGVL1SP3QL1WR2ND2RG*RRC2NTIVSNRYSPTGGAITV*1SV2AP9EL3W*1KNAF1MI1*CIA2LQ1SNDS1AD6"),
            (0, 110, 3, 394, 287, -3, -1, "3LF2VC1*CLE4FLNS7IM2NR7AL1"),
            (0, 394, 287, 110, 3, -1, -3, "3FL2CV1C*EL4LFSN7MI2RN7LA1"),
            (0, 161, 256, 395, 300, 2, -3, "IV3DY1CVPSVM5VI*KGL2VIVMNRND1TA1QE*E1SR2"),
            (0, 395, 300, 161, 256, -3, 2, "VI3YD1VCSPMV5IVK*LG2IVMVRNDN1AT1EQE*1RS2"),
            (0, 323, 394, 233, 162, 2, -3, "HY1RYAS*RYH1EGAL*KSP4ASWRTD2KQ3"),
            (0, 233, 162, 323, 394, -3, 2, "YH1YRSAR*HY1GELAK*PS4SARWDT2QK3"),
            (0, 272, 328, 297, 241, 2, -2, "3CYPGYHQH2RKMLVLGPTLPNRG2YH"),
            (0, 297, 241, 272, 328, -2, 2, "3YCGPHYHQ2KRLMLVPGLTNPGR2HY"),
            (0, 302, 388, 18, 104, 2, 3, "3TAPIRS1HYYNRC1*R2EKAV*K4LQAS1TD1DNKQ1"),
            (0, 18, 104, 302, 388, 3, 2, "3ATIPSR1YHNYCR1R*2KEVAK*4QLSA1DT1NDQK1"),
            (0, 387, 410, 285, 308, 3, 3, "KRDTDS2QE2"),
            (0, 285, 308, 387, 410, 3, 3, "RKTDSD2EQ2"),
            (0, 111, 16, 395, 300, -2, -3, "1IVFL1DYFLGVPS1TR2GS2*K4VM1NDYNTAND*E*ERC3"),
            (0, 395, 300, 111, 16, -3, -2, "1VILF1YDLFVGSP1RT2SG2K*4MV1DNNYATDNE*E*CR3"),
            (0, 270, 154, 2, 118, -2, 2, "KR1LRPL3LRHRHYW*SF1L*1FR1TI1RKPS1TI1LP2LVTM1QPRK2RKTMMTLA1"),
            (0, 2, 118, 270, 154, 2, -2, "RK1RLLP3RLRHYH*WFS1*L1RF1IT1KRSP1IT1PL2VLMT1PQKR2KRMTTMAL1"),
            (0, 301, 393, 255, 163, 1, -2, "2RL1SHSW1AV1SLRFMTIT2SPL*IT4RLMTDGTQ1*S3"),
            (0, 255, 163, 301, 393, -2, 1, "2LR1HSWS1VA1LSFRTMTI2PS*LTI4LRTMGDQT1S*3"),
            (0, 148, 201, 408, 355, 1, -2, "1VTRGKVHR3*L1LRSVRHEA4"),
            (0, 408, 355, 148, 201, -2, 1, "1TVGRVKRH3L*1RLVSHRAE4"),
            (0, 301, 396, 17, 112, 1, 2, "3HRSYS*SFAVL*SL1MT2RK1L*2LP2RV1DGTPRK*S1RKTM1"),
            (0, 17, 112, 301, 396, 2, 1, "3RHYS*SFSVA*LLS1TM2KR1*L2PL2VR1GDPTKRS*1KRMT1"),
            (0, 94, 14, 378, 298, -1, -2, "1SV1DA*R1RK2L*TAFS2RH1QRL*YCKRLRIGAV3QR"),
            (0, 378, 298, 94, 14, -2, -1, "1VS1ADR*1KR2*LATSF2HR1RQ*LCYRKRLGIVA3RQ"),
            (0, 387, 410, 1, 24, 3, 1, "1DTDA2QE2"),
            (0, 1, 24, 387, 410, 1, 3, "1TDAD2EQ2"),
            (0, 263, 286, 394, 371, 2, -1, "3VFRIWS2"),
            (0, 394, 371, 263, 286, -1, 2, "3FVIRSW2"),
        ],
        1,
    );
    assert_eq!(
        got,
        r#"
For reading frame +1 of query sequence:

txself  1    KTAAPEWSALLVCIIAHDMQKLNHTWFSHGTEIKKDDGYETSHVAPG*LTVRKHRPS*SL  180
txself  285  R.S...........................H......APQW.                    410
txself  387  .DD..Q..                                                      410

txself  181  SRERTKYGLRPACRE**Y*PVMQPIQGQQFCTMVSVPAHPNGRHSSSALSRMICRSLIIL  360
txself  17                                           ...RY*FV*L.T..K.*..P  76

txself  361  GSRMDTR*RRTTHPSGPSIL  420
txself  77   ..V.GPKS.KM.          112

For reading frame +2 of query sequence:

txself  2    RLRLLNGRRY*FV*LRTICKS*IIPGSVMGPKSRKMTATKHLT*HRDS*QFVSIVLLDLC  181
txself  301       ...HSSSALS.M..R.L..L..R.DTR*.RT.                         396

txself  182  PVREPSMV*GLHVVNNNTDQ*CSPFRGSSFVRWCPYQRTRMVGTPRLHYRA*YAEA*SYL  361
txself  18                                           ...AIS.YNC.R..KVK...  77

txself  362  VLAWTRDKEGRRTPVVHLS*  421
txself  78   .QS.D.NQ.             104

For reading frame +3 of query sequence:

txself  3    DCGS*MVGAISLYNCARYAKVKSYLVQSWDRNQER*RLRNISRSTGIVDSS*ASSFLISV  182
txself  302       ...TPR.HYR.*..EA*....LA.T.DK.                            388

txself  183  P*ENQVWFKACMS*IIILTSDAAHSGAAVLYDGVRTSAPEWSALLVCIIAHDMQKLNHTW  362
txself  1                                      K.A.......................  78
txself  387                                    KDD..Q..                    410

txself  363  FSHGHEIKKDDAPQWSIYP  419
txself  79   ....T......GYET.     126
txself  285          RTS..E..     308
txself  1            .TA..E..     24

For reading frame -1 of query sequence:

txself  421  LG*MDHWGASSFFISCPCENQV*LSFCISCAIMQTRSADHSGALVRTPSYKTAAPEWAAS  242
txself  110           ...L..V.*L...*FN.......I..N.......A.                 3

txself  241  LVSIIIHDMQALNHTWFSHGTEIKKDDAYELSTIPVLREMFRSRHLS*FRSHD*TRYDLT  62
txself  378                                                   .V.AR.K..*A  346

txself  61   FAYRAQLYKLIAPTIQEPQS  2
txself  345  S..H.R*CRRGV...R      298

For reading frame -2 of query sequence:

txself  420  *DRWTTGVRRPSLSRVHARTKYD*ASAYHAR*CRRGVPTIRVRWYGHHRTKLLPLNGLHH  241
txself  94                 .S.D*.R..LTF..R.QLYKLIA...Q                     14

txself  240  WSVLLFTTCRP*TILGSLTGQRSRRTMLTNCQLSRCYVRCFVAVIFLDFGPMTEPGMI*L  61
txself  395                                             .VL.YLVS.R..S..K.  345

txself  60   LHIVRNYTN**RRPFRSRSL  1
txself  344  ...M.DNADEEC...       300

For reading frame -3 of query sequence:

txself  419  RIDGPLGCVVLLYLVSMREPSMIKLLHIMRDNADEECRPFGCAGTDTIVQNCCP*MGCIT  240
txself  111          .IF.DFGP.T..G..*....V.NYTN**R...                      16

txself  239  GQYYYSRHAGLKPYLVLSRDRDQEGRCLRTVNYPGAT*DVS*PSSFLISVP*LNQV*FNF  60
txself  394                                             ...F..C.CE...*LS.  344

txself  59   CISCAIIQTNSADHSGAAV  3
txself  343  ......M..R.......L.  287

"#
    );
}

// line RG T3c probe `tx.t3c.self.f2`: a record with repeats on both strands against itself.
#[test]
fn tblastx_self_record_repeats_format_2() {
    let probe = Probe {
        query: (sequence("txself/txself"), "txself"),
        subjects: vec![(sequence("txself/txself"), "txself")],
        query_masks: vec![],
    };
    let got = section(
        &probe,
        0,
        &[
            (0, 420, 1, 420, 1, -2, -2, "140"),
            (0, 421, 2, 421, 2, -1, -1, "140"),
            (0, 419, 3, 419, 3, -3, -3, "139"),
            (0, 1, 420, 1, 420, 1, 1, "140"),
            (0, 3, 419, 3, 419, 3, 3, "139"),
            (0, 2, 421, 2, 421, 2, 2, "140"),
            (0, 138, 419, 418, 137, 3, -1, "1I*VM1SHSW*G4LF2VC1*C4W*FLKSAF1MI1*CIA1IMLQ1SRDS1AD4AL1LRYTDPGSVYRK1SA4SA1LS2CS2AI4KA9HT7PYQEWL1ITYI1"),
            (0, 418, 137, 138, 419, -1, 3, "1*IMV1HSWSG*4FL2CV1C*4*WLFSKFA1IM1C*AI1MIQL1RSSD1DA4LA1RLTYPDSGYVKR1AS4AS1SL2SC2IA4AK9TH7YPEQLW1TIIY1"),
            (0, 271, 140, 1, 132, -1, 1, "7AS1SL2SC2IA4AK16AG2LT1THIV"),
            (0, 1, 132, 271, 140, 1, -1, "7SA1LS2CS2AI4KA16GA2TL1HTVI"),
            (0, 120, 4, 152, 268, -2, 2, "2ASVIIVFL2FLGC1MVTR2GS1IV1LG2IV1RN1YN1ND*Q1RCRS3SGRS1"),
            (0, 152, 268, 120, 4, 2, -2, "2SAIVVILF2LFCG1VMRT2SG1VI1GL2VI1NR1NY1DNQ*1CRSR3GSSR1"),
            (0, 285, 410, 1, 126, 3, 1, "RK1SA27HT6AGPYQEWT1"),
            (0, 1, 126, 285, 410, 1, 3, "KR1AS27TH6GAYPEQTW1"),
            (0, 269, 3, 3, 269, -3, 3, "ND1CGPS2GVCGIATIGSQL1YNYCSA1HY1GKLV1PS3LQ1RW2DN2GRR*CR2TNVINSYRPSGTAGTI*V1VS2PA9LE3*W1NKFA1IM1C*AI2QL1NSSD1DA6"),
            (0, 3, 269, 269, 3, 3, -3, "DN1GCSP2VGGCAIITSGLQ1NYCYAS1YH1KGVL1SP3QL1WR2ND2RG*RRC2NTIVSNRYSPTGGAITV*1SV2AP9EL3W*1KNAF1MI1*CIA2LQ1SNDS1AD6"),
            (0, 110, 3, 394, 287, -3, -1, "3LF2VC1*CLE4FLNS7IM2NR7AL1"),
            (0, 394, 287, 110, 3, -1, -3, "3FL2CV1C*EL4LFSN7MI2RN7LA1"),
            (0, 161, 256, 395, 300, 2, -3, "IV3DY1CVPSVM5VI*KGL2VIVMNRND1TA1QE*E1SR2"),
            (0, 395, 300, 161, 256, -3, 2, "VI3YD1VCSPMV5IVK*LG2IVMVRNDN1AT1EQE*1RS2"),
            (0, 323, 394, 233, 162, 2, -3, "HY1RYAS*RYH1EGAL*KSP4ASWRTD2KQ3"),
            (0, 233, 162, 323, 394, -3, 2, "YH1YRSAR*HY1GELAK*PS4SARWDT2QK3"),
            (0, 272, 328, 297, 241, 2, -2, "3CYPGYHQH2RKMLVLGPTLPNRG2YH"),
            (0, 297, 241, 272, 328, -2, 2, "3YCGPHYHQ2KRLMLVPGLTNPGR2HY"),
            (0, 302, 388, 18, 104, 2, 3, "3TAPIRS1HYYNRC1*R2EKAV*K4LQAS1TD1DNKQ1"),
            (0, 18, 104, 302, 388, 3, 2, "3ATIPSR1YHNYCR1R*2KEVAK*4QLSA1DT1NDQK1"),
            (0, 387, 410, 285, 308, 3, 3, "KRDTDS2QE2"),
            (0, 285, 308, 387, 410, 3, 3, "RKTDSD2EQ2"),
            (0, 111, 16, 395, 300, -2, -3, "1IVFL1DYFLGVPS1TR2GS2*K4VM1NDYNTAND*E*ERC3"),
            (0, 395, 300, 111, 16, -3, -2, "1VILF1YDLFVGSP1RT2SG2K*4MV1DNNYATDNE*E*CR3"),
            (0, 270, 154, 2, 118, -2, 2, "KR1LRPL3LRHRHYW*SF1L*1FR1TI1RKPS1TI1LP2LVTM1QPRK2RKTMMTLA1"),
            (0, 2, 118, 270, 154, 2, -2, "RK1RLLP3RLRHYH*WFS1*L1RF1IT1KRSP1IT1PL2VLMT1PQKR2KRMTTMAL1"),
            (0, 301, 393, 255, 163, 1, -2, "2RL1SHSW1AV1SLRFMTIT2SPL*IT4RLMTDGTQ1*S3"),
            (0, 255, 163, 301, 393, -2, 1, "2LR1HSWS1VA1LSFRTMTI2PS*LTI4LRTMGDQT1S*3"),
            (0, 148, 201, 408, 355, 1, -2, "1VTRGKVHR3*L1LRSVRHEA4"),
            (0, 408, 355, 148, 201, -2, 1, "1TVGRVKRH3L*1RLVSHRAE4"),
            (0, 301, 396, 17, 112, 1, 2, "3HRSYS*SFAVL*SL1MT2RK1L*2LP2RV1DGTPRK*S1RKTM1"),
            (0, 17, 112, 301, 396, 2, 1, "3RHYS*SFSVA*LLS1TM2KR1*L2PL2VR1GDPTKRS*1KRMT1"),
            (0, 94, 14, 378, 298, -1, -2, "1SV1DA*R1RK2L*TAFS2RH1QRL*YCKRLRIGAV3QR"),
            (0, 378, 298, 94, 14, -2, -1, "1VS1ADR*1KR2*LATSF2HR1RQ*LCYRKRLGIVA3RQ"),
            (0, 387, 410, 1, 24, 3, 1, "1DTDA2QE2"),
            (0, 1, 24, 387, 410, 1, 3, "1TDAD2EQ2"),
            (0, 263, 286, 394, 371, 2, -1, "3VFRIWS2"),
            (0, 394, 371, 263, 286, -1, 2, "3FVIRSW2"),
        ],
        2,
    );
    assert_eq!(
        got,
        r#"
For reading frame +1 of query sequence:

txself  1    KTAAPEWSALLVCIIAHDMQKLNHTWFSHGTEIKKDDGYETSHVAPG*LTVRKHRPS*SL  180
txself  285  RTSAPEWSALLVCIIAHDMQKLNHTWFSHGHEIKKDDAPQWS                    410
txself  387  KDDAPQWS                                                      410

txself  181  SRERTKYGLRPACRE**Y*PVMQPIQGQQFCTMVSVPAHPNGRHSSSALSRMICRSLIIL  360
txself  17                                           NGRRY*FV*LRTICKS*IIP  76

txself  361  GSRMDTR*RRTTHPSGPSIL  420
txself  77   GSVMGPKSRKMT          112

For reading frame +2 of query sequence:

txself  2    RLRLLNGRRY*FV*LRTICKS*IIPGSVMGPKSRKMTATKHLT*HRDS*QFVSIVLLDLC  181
txself  301       NGRHSSSALSRMICRSLIILGSRMDTR*RRTT                         396

txself  182  PVREPSMV*GLHVVNNNTDQ*CSPFRGSSFVRWCPYQRTRMVGTPRLHYRA*YAEA*SYL  361
txself  18                                           MVGAISLYNCARYAKVKSYL  77

txself  362  VLAWTRDKEGRRTPVVHLS*  421
txself  78   VQSWDRNQE             104

For reading frame +3 of query sequence:

txself  3    DCGS*MVGAISLYNCARYAKVKSYLVQSWDRNQER*RLRNISRSTGIVDSS*ASSFLISV  182
txself  302       MVGTPRLHYRA*YAEA*SYLVLAWTRDKE                            388

txself  183  P*ENQVWFKACMS*IIILTSDAAHSGAAVLYDGVRTSAPEWSALLVCIIAHDMQKLNHTW  362
txself  1                                      KTAAPEWSALLVCIIAHDMQKLNHTW  78
txself  387                                    KDDAPQWS                    410

txself  363  FSHGHEIKKDDAPQWSIYP  419
txself  79   FSHGTEIKKDDGYETS     126
txself  285          RTSAPEWS     308
txself  1            KTAAPEWS     24

For reading frame -1 of query sequence:

txself  421  LG*MDHWGASSFFISCPCENQV*LSFCISCAIMQTRSADHSGALVRTPSYKTAAPEWAAS  242
txself  110           SSFLISVP*LNQV*FNFCISCAIIQTNSADHSGAAV                 3

txself  241  LVSIIIHDMQALNHTWFSHGTEIKKDDAYELSTIPVLREMFRSRHLS*FRSHD*TRYDLT  62
txself  378                                                   RVHARTKYD*A  346

txself  61   FAYRAQLYKLIAPTIQEPQS  2
txself  345  SAYHAR*CRRGVPTIR      298

For reading frame -2 of query sequence:

txself  420  *DRWTTGVRRPSLSRVHARTKYD*ASAYHAR*CRRGVPTIRVRWYGHHRTKLLPLNGLHH  241
txself  94                 RSHD*TRYDLTFAYRAQLYKLIAPTIQ                     14

txself  240  WSVLLFTTCRP*TILGSLTGQRSRRTMLTNCQLSRCYVRCFVAVIFLDFGPMTEPGMI*L  61
txself  395                                             VVLLYLVSMREPSMIKL  345

txself  60   LHIVRNYTN**RRPFRSRSL  1
txself  344  LHIMRDNADEECRPF       300

For reading frame -3 of query sequence:

txself  419  RIDGPLGCVVLLYLVSMREPSMIKLLHIMRDNADEECRPFGCAGTDTIVQNCCP*MGCIT  240
txself  111          VIFLDFGPMTEPGMI*LLHIVRNYTN**RRPF                      16

txself  239  GQYYYSRHAGLKPYLVLSRDRDQEGRCLRTVNYPGAT*DVS*PSSFLISVP*LNQV*FNF  60
txself  394                                             SSFFISCPCENQV*LSF  344

txself  59   CISCAIIQTNSADHSGAAV  3
txself  343  CISCAIMQTRSADHSGALV  287

"#
    );
}

// line RG T3c probe `tx.t3c.self.f3`: a record with repeats on both strands against itself.
#[test]
fn tblastx_self_record_repeats_format_3() {
    let probe = Probe {
        query: (sequence("txself/txself"), "txself"),
        subjects: vec![(sequence("txself/txself"), "txself")],
        query_masks: vec![],
    };
    let got = section(
        &probe,
        0,
        &[
            (0, 420, 1, 420, 1, -2, -2, "140"),
            (0, 421, 2, 421, 2, -1, -1, "140"),
            (0, 419, 3, 419, 3, -3, -3, "139"),
            (0, 1, 420, 1, 420, 1, 1, "140"),
            (0, 3, 419, 3, 419, 3, 3, "139"),
            (0, 2, 421, 2, 421, 2, 2, "140"),
            (0, 138, 419, 418, 137, 3, -1, "1I*VM1SHSW*G4LF2VC1*C4W*FLKSAF1MI1*CIA1IMLQ1SRDS1AD4AL1LRYTDPGSVYRK1SA4SA1LS2CS2AI4KA9HT7PYQEWL1ITYI1"),
            (0, 418, 137, 138, 419, -1, 3, "1*IMV1HSWSG*4FL2CV1C*4*WLFSKFA1IM1C*AI1MIQL1RSSD1DA4LA1RLTYPDSGYVKR1AS4AS1SL2SC2IA4AK9TH7YPEQLW1TIIY1"),
            (0, 271, 140, 1, 132, -1, 1, "7AS1SL2SC2IA4AK16AG2LT1THIV"),
            (0, 1, 132, 271, 140, 1, -1, "7SA1LS2CS2AI4KA16GA2TL1HTVI"),
            (0, 120, 4, 152, 268, -2, 2, "2ASVIIVFL2FLGC1MVTR2GS1IV1LG2IV1RN1YN1ND*Q1RCRS3SGRS1"),
            (0, 152, 268, 120, 4, 2, -2, "2SAIVVILF2LFCG1VMRT2SG1VI1GL2VI1NR1NY1DNQ*1CRSR3GSSR1"),
            (0, 285, 410, 1, 126, 3, 1, "RK1SA27HT6AGPYQEWT1"),
            (0, 1, 126, 285, 410, 1, 3, "KR1AS27TH6GAYPEQTW1"),
            (0, 269, 3, 3, 269, -3, 3, "ND1CGPS2GVCGIATIGSQL1YNYCSA1HY1GKLV1PS3LQ1RW2DN2GRR*CR2TNVINSYRPSGTAGTI*V1VS2PA9LE3*W1NKFA1IM1C*AI2QL1NSSD1DA6"),
            (0, 3, 269, 269, 3, 3, -3, "DN1GCSP2VGGCAIITSGLQ1NYCYAS1YH1KGVL1SP3QL1WR2ND2RG*RRC2NTIVSNRYSPTGGAITV*1SV2AP9EL3W*1KNAF1MI1*CIA2LQ1SNDS1AD6"),
            (0, 110, 3, 394, 287, -3, -1, "3LF2VC1*CLE4FLNS7IM2NR7AL1"),
            (0, 394, 287, 110, 3, -1, -3, "3FL2CV1C*EL4LFSN7MI2RN7LA1"),
            (0, 161, 256, 395, 300, 2, -3, "IV3DY1CVPSVM5VI*KGL2VIVMNRND1TA1QE*E1SR2"),
            (0, 395, 300, 161, 256, -3, 2, "VI3YD1VCSPMV5IVK*LG2IVMVRNDN1AT1EQE*1RS2"),
            (0, 323, 394, 233, 162, 2, -3, "HY1RYAS*RYH1EGAL*KSP4ASWRTD2KQ3"),
            (0, 233, 162, 323, 394, -3, 2, "YH1YRSAR*HY1GELAK*PS4SARWDT2QK3"),
            (0, 272, 328, 297, 241, 2, -2, "3CYPGYHQH2RKMLVLGPTLPNRG2YH"),
            (0, 297, 241, 272, 328, -2, 2, "3YCGPHYHQ2KRLMLVPGLTNPGR2HY"),
            (0, 302, 388, 18, 104, 2, 3, "3TAPIRS1HYYNRC1*R2EKAV*K4LQAS1TD1DNKQ1"),
            (0, 18, 104, 302, 388, 3, 2, "3ATIPSR1YHNYCR1R*2KEVAK*4QLSA1DT1NDQK1"),
            (0, 387, 410, 285, 308, 3, 3, "KRDTDS2QE2"),
            (0, 285, 308, 387, 410, 3, 3, "RKTDSD2EQ2"),
            (0, 111, 16, 395, 300, -2, -3, "1IVFL1DYFLGVPS1TR2GS2*K4VM1NDYNTAND*E*ERC3"),
            (0, 395, 300, 111, 16, -3, -2, "1VILF1YDLFVGSP1RT2SG2K*4MV1DNNYATDNE*E*CR3"),
            (0, 270, 154, 2, 118, -2, 2, "KR1LRPL3LRHRHYW*SF1L*1FR1TI1RKPS1TI1LP2LVTM1QPRK2RKTMMTLA1"),
            (0, 2, 118, 270, 154, 2, -2, "RK1RLLP3RLRHYH*WFS1*L1RF1IT1KRSP1IT1PL2VLMT1PQKR2KRMTTMAL1"),
            (0, 301, 393, 255, 163, 1, -2, "2RL1SHSW1AV1SLRFMTIT2SPL*IT4RLMTDGTQ1*S3"),
            (0, 255, 163, 301, 393, -2, 1, "2LR1HSWS1VA1LSFRTMTI2PS*LTI4LRTMGDQT1S*3"),
            (0, 148, 201, 408, 355, 1, -2, "1VTRGKVHR3*L1LRSVRHEA4"),
            (0, 408, 355, 148, 201, -2, 1, "1TVGRVKRH3L*1RLVSHRAE4"),
            (0, 301, 396, 17, 112, 1, 2, "3HRSYS*SFAVL*SL1MT2RK1L*2LP2RV1DGTPRK*S1RKTM1"),
            (0, 17, 112, 301, 396, 2, 1, "3RHYS*SFSVA*LLS1TM2KR1*L2PL2VR1GDPTKRS*1KRMT1"),
            (0, 94, 14, 378, 298, -1, -2, "1SV1DA*R1RK2L*TAFS2RH1QRL*YCKRLRIGAV3QR"),
            (0, 378, 298, 94, 14, -2, -1, "1VS1ADR*1KR2*LATSF2HR1RQ*LCYRKRLGIVA3RQ"),
            (0, 387, 410, 1, 24, 3, 1, "1DTDA2QE2"),
            (0, 1, 24, 387, 410, 1, 3, "1TDAD2EQ2"),
            (0, 263, 286, 394, 371, 2, -1, "3VFRIWS2"),
            (0, 394, 371, 263, 286, -1, 2, "3FVIRSW2"),
        ],
        3,
    );
    assert_eq!(
        got,
        r#"
For reading frame +1 of query sequence:

txself  1    KTAAPEWSALLVCIIAHDMQKLNHTWFSHGTEIKKDDGYETSHVAPG*LTVRKHRPS*SL  180
txself  285  R.S...........................H......APQW.                    410
txself  387  .DD..Q..                                                      410

txself  181  SRERTKYGLRPACRE**Y*PVMQPIQGQQFCTMVSVPAHPNGRHSSSALSRMICRSLIIL  360
txself  17                                           ...RY*FV*L.T..K.*..P  76

txself  361  GSRMDTR*RRTTHPSGPSIL  420
txself  77   ..V.GPKS.KM.          112

For reading frame +2 of query sequence:

txself  2    RLRLLNGRRY*FV*LRTICKS*IIPGSVMGPKSRKMTATKHLT*HRDS*QFVSIVLLDLC  181
txself  301       ...HSSSALS.M..R.L..L..R.DTR*.RT.                         396

txself  182  PVREPSMV*GLHVVNNNTDQ*CSPFRGSSFVRWCPYQRTRMVGTPRLHYRA*YAEA*SYL  361
txself  18                                           ...AIS.YNC.R..KVK...  77

txself  362  VLAWTRDKEGRRTPVVHLS*  421
txself  78   .QS.D.NQ.             104

For reading frame +3 of query sequence:

txself  3    DCGS*MVGAISLYNCARYAKVKSYLVQSWDRNQER*RLRNISRSTGIVDSS*ASSFLISV  182
txself  302       ...TPR.HYR.*..EA*....LA.T.DK.                            388

txself  183  P*ENQVWFKACMS*IIILTSDAAHSGAAVLYDGVRTSAPEWSALLVCIIAHDMQKLNHTW  362
txself  1                                      K.A.......................  78
txself  387                                    KDD..Q..                    410

txself  363  FSHGHEIKKDDAPQWSIYP  419
txself  79   ....T......GYET.     126
txself  285          RTS..E..     308
txself  1            .TA..E..     24

For reading frame -1 of query sequence:

txself  421  LG*MDHWGASSFFISCPCENQV*LSFCISCAIMQTRSADHSGALVRTPSYKTAAPEWAAS  242
txself  110           ...L..V.*L...*FN.......I..N.......A.                 3

txself  241  LVSIIIHDMQALNHTWFSHGTEIKKDDAYELSTIPVLREMFRSRHLS*FRSHD*TRYDLT  62
txself  378                                                   .V.AR.K..*A  346

txself  61   FAYRAQLYKLIAPTIQEPQS  2
txself  345  S..H.R*CRRGV...R      298

For reading frame -2 of query sequence:

txself  420  *DRWTTGVRRPSLSRVHARTKYD*ASAYHAR*CRRGVPTIRVRWYGHHRTKLLPLNGLHH  241
txself  94                 .S.D*.R..LTF..R.QLYKLIA...Q                     14

txself  240  WSVLLFTTCRP*TILGSLTGQRSRRTMLTNCQLSRCYVRCFVAVIFLDFGPMTEPGMI*L  61
txself  395                                             .VL.YLVS.R..S..K.  345

txself  60   LHIVRNYTN**RRPFRSRSL  1
txself  344  ...M.DNADEEC...       300

For reading frame -3 of query sequence:

txself  419  RIDGPLGCVVLLYLVSMREPSMIKLLHIMRDNADEECRPFGCAGTDTIVQNCCP*MGCIT  240
txself  111          .IF.DFGP.T..G..*....V.NYTN**R...                      16

txself  239  GQYYYSRHAGLKPYLVLSRDRDQEGRCLRTVNYPGAT*DVS*PSSFLISVP*LNQV*FNF  60
txself  394                                             ...F..C.CE...*LS.  344

txself  59   CISCAIIQTNSADHSGAAV  3
txself  343  ......M..R.......L.  287

"#
    );
}

// line RG T3c probe `tx.t3c.self.f4`: a record with repeats on both strands against itself.
#[test]
fn tblastx_self_record_repeats_format_4() {
    let probe = Probe {
        query: (sequence("txself/txself"), "txself"),
        subjects: vec![(sequence("txself/txself"), "txself")],
        query_masks: vec![],
    };
    let got = section(
        &probe,
        0,
        &[
            (0, 420, 1, 420, 1, -2, -2, "140"),
            (0, 421, 2, 421, 2, -1, -1, "140"),
            (0, 419, 3, 419, 3, -3, -3, "139"),
            (0, 1, 420, 1, 420, 1, 1, "140"),
            (0, 3, 419, 3, 419, 3, 3, "139"),
            (0, 2, 421, 2, 421, 2, 2, "140"),
            (0, 138, 419, 418, 137, 3, -1, "1I*VM1SHSW*G4LF2VC1*C4W*FLKSAF1MI1*CIA1IMLQ1SRDS1AD4AL1LRYTDPGSVYRK1SA4SA1LS2CS2AI4KA9HT7PYQEWL1ITYI1"),
            (0, 418, 137, 138, 419, -1, 3, "1*IMV1HSWSG*4FL2CV1C*4*WLFSKFA1IM1C*AI1MIQL1RSSD1DA4LA1RLTYPDSGYVKR1AS4AS1SL2SC2IA4AK9TH7YPEQLW1TIIY1"),
            (0, 271, 140, 1, 132, -1, 1, "7AS1SL2SC2IA4AK16AG2LT1THIV"),
            (0, 1, 132, 271, 140, 1, -1, "7SA1LS2CS2AI4KA16GA2TL1HTVI"),
            (0, 120, 4, 152, 268, -2, 2, "2ASVIIVFL2FLGC1MVTR2GS1IV1LG2IV1RN1YN1ND*Q1RCRS3SGRS1"),
            (0, 152, 268, 120, 4, 2, -2, "2SAIVVILF2LFCG1VMRT2SG1VI1GL2VI1NR1NY1DNQ*1CRSR3GSSR1"),
            (0, 285, 410, 1, 126, 3, 1, "RK1SA27HT6AGPYQEWT1"),
            (0, 1, 126, 285, 410, 1, 3, "KR1AS27TH6GAYPEQTW1"),
            (0, 269, 3, 3, 269, -3, 3, "ND1CGPS2GVCGIATIGSQL1YNYCSA1HY1GKLV1PS3LQ1RW2DN2GRR*CR2TNVINSYRPSGTAGTI*V1VS2PA9LE3*W1NKFA1IM1C*AI2QL1NSSD1DA6"),
            (0, 3, 269, 269, 3, 3, -3, "DN1GCSP2VGGCAIITSGLQ1NYCYAS1YH1KGVL1SP3QL1WR2ND2RG*RRC2NTIVSNRYSPTGGAITV*1SV2AP9EL3W*1KNAF1MI1*CIA2LQ1SNDS1AD6"),
            (0, 110, 3, 394, 287, -3, -1, "3LF2VC1*CLE4FLNS7IM2NR7AL1"),
            (0, 394, 287, 110, 3, -1, -3, "3FL2CV1C*EL4LFSN7MI2RN7LA1"),
            (0, 161, 256, 395, 300, 2, -3, "IV3DY1CVPSVM5VI*KGL2VIVMNRND1TA1QE*E1SR2"),
            (0, 395, 300, 161, 256, -3, 2, "VI3YD1VCSPMV5IVK*LG2IVMVRNDN1AT1EQE*1RS2"),
            (0, 323, 394, 233, 162, 2, -3, "HY1RYAS*RYH1EGAL*KSP4ASWRTD2KQ3"),
            (0, 233, 162, 323, 394, -3, 2, "YH1YRSAR*HY1GELAK*PS4SARWDT2QK3"),
            (0, 272, 328, 297, 241, 2, -2, "3CYPGYHQH2RKMLVLGPTLPNRG2YH"),
            (0, 297, 241, 272, 328, -2, 2, "3YCGPHYHQ2KRLMLVPGLTNPGR2HY"),
            (0, 302, 388, 18, 104, 2, 3, "3TAPIRS1HYYNRC1*R2EKAV*K4LQAS1TD1DNKQ1"),
            (0, 18, 104, 302, 388, 3, 2, "3ATIPSR1YHNYCR1R*2KEVAK*4QLSA1DT1NDQK1"),
            (0, 387, 410, 285, 308, 3, 3, "KRDTDS2QE2"),
            (0, 285, 308, 387, 410, 3, 3, "RKTDSD2EQ2"),
            (0, 111, 16, 395, 300, -2, -3, "1IVFL1DYFLGVPS1TR2GS2*K4VM1NDYNTAND*E*ERC3"),
            (0, 395, 300, 111, 16, -3, -2, "1VILF1YDLFVGSP1RT2SG2K*4MV1DNNYATDNE*E*CR3"),
            (0, 270, 154, 2, 118, -2, 2, "KR1LRPL3LRHRHYW*SF1L*1FR1TI1RKPS1TI1LP2LVTM1QPRK2RKTMMTLA1"),
            (0, 2, 118, 270, 154, 2, -2, "RK1RLLP3RLRHYH*WFS1*L1RF1IT1KRSP1IT1PL2VLMT1PQKR2KRMTTMAL1"),
            (0, 301, 393, 255, 163, 1, -2, "2RL1SHSW1AV1SLRFMTIT2SPL*IT4RLMTDGTQ1*S3"),
            (0, 255, 163, 301, 393, -2, 1, "2LR1HSWS1VA1LSFRTMTI2PS*LTI4LRTMGDQT1S*3"),
            (0, 148, 201, 408, 355, 1, -2, "1VTRGKVHR3*L1LRSVRHEA4"),
            (0, 408, 355, 148, 201, -2, 1, "1TVGRVKRH3L*1RLVSHRAE4"),
            (0, 301, 396, 17, 112, 1, 2, "3HRSYS*SFAVL*SL1MT2RK1L*2LP2RV1DGTPRK*S1RKTM1"),
            (0, 17, 112, 301, 396, 2, 1, "3RHYS*SFSVA*LLS1TM2KR1*L2PL2VR1GDPTKRS*1KRMT1"),
            (0, 94, 14, 378, 298, -1, -2, "1SV1DA*R1RK2L*TAFS2RH1QRL*YCKRLRIGAV3QR"),
            (0, 378, 298, 94, 14, -2, -1, "1VS1ADR*1KR2*LATSF2HR1RQ*LCYRKRLGIVA3RQ"),
            (0, 387, 410, 1, 24, 3, 1, "1DTDA2QE2"),
            (0, 1, 24, 387, 410, 1, 3, "1TDAD2EQ2"),
            (0, 263, 286, 394, 371, 2, -1, "3VFRIWS2"),
            (0, 394, 371, 263, 286, -1, 2, "3FVIRSW2"),
        ],
        4,
    );
    assert_eq!(
        got,
        r#"
For reading frame +1 of query sequence:

txself  1    KTAAPEWSALLVCIIAHDMQKLNHTWFSHGTEIKKDDGYETSHVAPG*LTVRKHRPS*SL  180
txself  285  RTSAPEWSALLVCIIAHDMQKLNHTWFSHGHEIKKDDAPQWS                    410
txself  387  KDDAPQWS                                                      410

txself  181  SRERTKYGLRPACRE**Y*PVMQPIQGQQFCTMVSVPAHPNGRHSSSALSRMICRSLIIL  360
txself  17                                           NGRRY*FV*LRTICKS*IIP  76

txself  361  GSRMDTR*RRTTHPSGPSIL  420
txself  77   GSVMGPKSRKMT          112

For reading frame +2 of query sequence:

txself  2    RLRLLNGRRY*FV*LRTICKS*IIPGSVMGPKSRKMTATKHLT*HRDS*QFVSIVLLDLC  181
txself  301       NGRHSSSALSRMICRSLIILGSRMDTR*RRTT                         396

txself  182  PVREPSMV*GLHVVNNNTDQ*CSPFRGSSFVRWCPYQRTRMVGTPRLHYRA*YAEA*SYL  361
txself  18                                           MVGAISLYNCARYAKVKSYL  77

txself  362  VLAWTRDKEGRRTPVVHLS*  421
txself  78   VQSWDRNQE             104

For reading frame +3 of query sequence:

txself  3    DCGS*MVGAISLYNCARYAKVKSYLVQSWDRNQER*RLRNISRSTGIVDSS*ASSFLISV  182
txself  302       MVGTPRLHYRA*YAEA*SYLVLAWTRDKE                            388

txself  183  P*ENQVWFKACMS*IIILTSDAAHSGAAVLYDGVRTSAPEWSALLVCIIAHDMQKLNHTW  362
txself  1                                      KTAAPEWSALLVCIIAHDMQKLNHTW  78
txself  387                                    KDDAPQWS                    410

txself  363  FSHGHEIKKDDAPQWSIYP  419
txself  79   FSHGTEIKKDDGYETS     126
txself  285          RTSAPEWS     308
txself  1            KTAAPEWS     24

For reading frame -1 of query sequence:

txself  421  LG*MDHWGASSFFISCPCENQV*LSFCISCAIMQTRSADHSGALVRTPSYKTAAPEWAAS  242
txself  110           SSFLISVP*LNQV*FNFCISCAIIQTNSADHSGAAV                 3

txself  241  LVSIIIHDMQALNHTWFSHGTEIKKDDAYELSTIPVLREMFRSRHLS*FRSHD*TRYDLT  62
txself  378                                                   RVHARTKYD*A  346

txself  61   FAYRAQLYKLIAPTIQEPQS  2
txself  345  SAYHAR*CRRGVPTIR      298

For reading frame -2 of query sequence:

txself  420  *DRWTTGVRRPSLSRVHARTKYD*ASAYHAR*CRRGVPTIRVRWYGHHRTKLLPLNGLHH  241
txself  94                 RSHD*TRYDLTFAYRAQLYKLIAPTIQ                     14

txself  240  WSVLLFTTCRP*TILGSLTGQRSRRTMLTNCQLSRCYVRCFVAVIFLDFGPMTEPGMI*L  61
txself  395                                             VVLLYLVSMREPSMIKL  345

txself  60   LHIVRNYTN**RRPFRSRSL  1
txself  344  LHIMRDNADEECRPF       300

For reading frame -3 of query sequence:

txself  419  RIDGPLGCVVLLYLVSMREPSMIKLLHIMRDNADEECRPFGCAGTDTIVQNCCP*MGCIT  240
txself  111          VIFLDFGPMTEPGMI*LLHIVRNYTN**RRPF                      16

txself  239  GQYYYSRHAGLKPYLVLSRDRDQEGRCLRTVNYPGAT*DVS*PSSFLISVP*LNQV*FNF  60
txself  394                                             SSFFISCPCENQV*LSF  344

txself  59   CISCAIIQTNSADHSGAAV  3
txself  343  CISCAIMQTRSADHSGALV  287

"#
    );
}

/// The section of a TBLASTX query whose own record is subject `own` of `probe` (`rows`: subject, query
/// and subject start and end, their frames and BTOP, in NCBI's order), on which NCBI's merge throws: not
/// an abort (`ncbi_aborts`), and `write_section` writes nothing of it and fails with NCBI's exception.
fn thrown(
    probe: &Probe,
    own: u32,
    rows: &[(u32, u32, u32, u32, u32, i8, i8, &str)],
) -> NcbiException {
    let mut hsps = Hsps::of_translated_nucleotides();
    hsps.set_query_bioseq(own);
    let mut segments = Segments::default();
    for &(subject, qstart, qend, sstart, send, qframe, sframe, btop) in rows {
        let hit = hit(
            qstart,
            qend,
            sstart,
            send,
            i32::from(qframe),
            ops_of_btop(btop, false),
        );
        let strands = crate::report::serial::seqalign::translated_row_segments(
            &hit,
            Some((qframe, probe.query.0.len())),
            Some((sframe, probe.subjects[subject as usize].0.len())),
            &mut segments,
        );
        hsps.push_dense_seg(subject, segments.starts(), segments.lens(), 0, 0, strands);
    }
    assert!(!ncbi_aborts(&hsps).expect("the merge"));
    let view = View {
        identity: true,
        anchored: true,
    };
    let mut out = Vec::new();
    let error = write_section(&mut out, view, 60, &hsps, probe).expect_err("NCBI's merge throws");
    assert!(out.is_empty(), "{}", String::from_utf8_lossy(&out));
    crate::ncbi_exception::NcbiException::of_io_error(&error)
        .expect("NCBI's exception")
        .clone()
}

// probe `f1.rep2.qs.f1` of `T/t3/fix1_t3c/cases_fix1.tsv` (`-query_loc 1-300 -subject_loc 301-600`): the
// record's two copies kept apart, so the first HSP of the record with itself that the merge takes puts
// subject positions on the query's own row where it has none, and NCBI's `x_Merge` throws (exit 255:
// `Start iterator already exists for seq2.` at alnmerger.cpp:790; the whole report is fixture row
// `db.of1.tx.throw.si.f1`).
#[test]
fn tblastx_self_record_merge_throws() {
    let letters = sequence("rep2/rep2");
    let probe = Probe {
        query: (letters, "rep2"),
        subjects: vec![(letters, "rep2")],
        query_masks: vec![],
    };
    let rows: &[(u32, u32, u32, u32, u32, i8, i8, &str)] = &[
        (0, 1, 300, 301, 600, 1, 1, "100"),
        (0, 298, 2, 598, 302, -3, -3, "99"),
        (0, 300, 1, 600, 301, -1, -1, "100"),
        (0, 3, 299, 303, 599, 3, 3, "99"),
        (0, 299, 3, 599, 303, -2, -2, "99"),
        (0, 2, 298, 302, 598, 2, 2, "99"),
        (
            0,
            104,
            166,
            429,
            367,
            2,
            -1,
            "2NDRDMISRSF2AKCEILDGLASTLTLVENVLRT1",
        ),
        (
            0,
            129,
            67,
            404,
            466,
            -1,
            2,
            "2DNDRIMRSFS2KAECLIGDALTSTLVLNELVTR1",
        ),
    ];
    let exception = thrown(&probe, 0, rows);
    assert_eq!(
        exception.site(),
        crate::ncbi_exception::ThrowSite {
            file: "src/objtools/alnmgr/alnmerger.cpp",
            line: 790,
            kind: "CAlnException::eMergeFailure",
            function: Some("ncbi::objects::CAlnMixMerger::x_Merge"),
        }
    );
    assert_eq!(
        exception.get_msg(),
        b"CAlnMixMerger::x_Merge(): Internal error: Start iterator already exists for seq2."
    );
}

// review 2 of line RG T3c fix 1, case `mg0.sx.f1` of `T/t3/review2_t3c/cases_mg.tsv` (`-subject_loc 154-303
// -seg no`; the HSPs of `mg0.sx.f6`): a tandem record `A+A` and a subject `x` that holds `A`. The first HSP
// of the record with itself that the merge takes has starts of the query's other matches as its subject
// positions, so the walk joins two segments of the query's row (alnmerger.cpp:753-771), and `Build` meets
// the loop and throws (exit 255: `Infinite loop detected.` at alnsegments.cpp:168; fixture row
// `db.of1.tx.throw.il.f1`).
#[test]
fn tblastx_self_record_build_throws() {
    let probe = Probe {
        query: (sequence("mg0_q/mg0"), "mg0"),
        subjects: vec![
            (sequence("mg0_sx/mg0"), "mg0"),
            (sequence("mg0_sx/mx0"), "mx0"),
        ],
        query_masks: vec![],
    };
    let rows: &[(u32, u32, u32, u32, u32, i8, i8, &str)] = &[
        (1, 303, 154, 303, 154, -1, -1, "50"),
        (1, 150, 1, 303, 154, -1, -1, "50"),
        (1, 1, 150, 154, 303, 1, 1, "50"),
        (1, 154, 303, 154, 303, 1, 1, "50"),
        (1, 3, 149, 156, 302, 3, 3, "49"),
        (1, 156, 302, 156, 302, 3, 3, "49"),
        (1, 301, 155, 301, 155, -3, -3, "49"),
        (1, 148, 2, 301, 155, -3, -3, "49"),
        (1, 2, 148, 155, 301, 2, 2, "49"),
        (1, 155, 301, 155, 301, 2, 2, "49"),
        (1, 302, 156, 302, 156, -2, -2, "49"),
        (1, 149, 3, 302, 156, -2, -2, "49"),
        (0, 303, 154, 303, 154, -1, -1, "50"),
        (0, 150, 1, 303, 154, -1, -1, "50"),
        (0, 1, 150, 154, 303, 1, 1, "50"),
        (0, 154, 303, 154, 303, 1, 1, "50"),
        (0, 3, 149, 156, 302, 3, 3, "49"),
        (0, 156, 302, 156, 302, 3, 3, "49"),
        (0, 301, 155, 301, 155, -3, -3, "49"),
        (0, 148, 2, 301, 155, -3, -3, "49"),
        (0, 2, 148, 155, 301, 2, 2, "49"),
        (0, 155, 301, 155, 301, 2, 2, "49"),
        (0, 302, 156, 302, 156, -2, -2, "49"),
        (0, 149, 3, 302, 156, -2, -2, "49"),
    ];
    let exception = thrown(&probe, 0, rows);
    assert_eq!(
        exception.site(),
        crate::ncbi_exception::ThrowSite {
            file: "src/objtools/alnmgr/alnsegments.cpp",
            line: 168,
            kind: "CAlnException::eMergeFailure",
            function: Some("ncbi::objects::CAlnMixSegments::Build"),
        }
    );
    assert_eq!(
        exception.get_msg(),
        b"CAlnMixSegments::Build(): Internal error: Infinite loop detected."
    );
}
