use super::*;

// NCBI reference (598d8ae6): c++/src/objtools/format/cigar_formatter.cpp:110-125
// The post of one Seq-align whose Std-segs change a strand, as NCBI BLAST+ 2.17.0 writes it in
// BLASTN's application (the oracle's probe `tblastx.sub.ss_.err` of line RG F2 fix 3), in LOSAT's
// words.
#[test]
fn the_post_of_a_seq_align_with_inconsistent_strands_is_ncbis_text_in_losat_words() {
    let ncbi = format!(
        "Error: (1439.1) [blastn] [CCIGAR_Formatter::x_FormatAlignmentRows] Exception: NCBI C++ Exception:\n    T0 \"{}/src/objects/seq/../seqalign/Seq_align.cpp\", line 827: Error: (CSeqalignException::eInvalidInputAlignment) CreateDensegFromStdseg(): Inconsistent strands!\n\n",
        crate::ncbi_exception::NCBI_BUILD_ROOT
    );
    assert_eq!(
        inconsistent_strands_post("blastn"),
        crate::wording::losat_text(&ncbi).as_bytes()
    );
}
