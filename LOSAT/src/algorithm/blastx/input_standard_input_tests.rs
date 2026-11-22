use super::*;
use crate::blastinput::input_files::InputFile;

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:717-721 (`-` is `cin`) and
// c++/src/corelib/stream_utils.cpp:235-236 (a refill asks for `in_avail()` bytes), cited at
// `fasta_reader::stream::InputBackend::of`.
/// BLASTX's and RPSTBLASTN's query stream reads standard input as NCBI reads `cin`, by the
/// readers' one rule (`InputBackend::of`): the tail that the same bytes lose in a file named by
/// its path is kept (NCBI BLAST+ 2.17.0, blastx and rpstblastn: `Length=107` from a pipe and
/// from a redirected file, 85 from the path; fixture rows `blastx.stdin_q_lost_g55_pipe`,
/// `blastx.stdin_q_lost_g55_file`, `db.rpstbn.g55.*`).
#[test]
fn standard_input_keeps_the_tail_that_a_named_file_loses() {
    let text = b">first\r\nTC-TTGGCTCAATCCTAGGTGGGCATGTTTCCTAATGCCC\rTT-TTTAACGTGAGGGTTCGCGTTTTTATCCCACCTAGC\r\r\nACGTACGT\n-ACGTACGTACGTACGTACGTAC";
    let path = crate::utils::test_dir::temp_dir().join(format!(
        "losat-blastx-standard-input-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::write(&path, text).unwrap();
    let lengths = |input: InputFile| {
        let mut stream = FastaStream::from_file(input);
        let records = parse_fasta_stream_keeping(
            &mut stream,
            false,
            false,
            usize::MAX,
            false,
            None,
            &mut |_| Ok(()),
            &mut |_| Ok(()),
        )
        .unwrap();
        records
            .iter()
            .map(|record| record.sequence.len())
            .collect::<Vec<_>>()
    };
    let named = lengths(InputFile::from(std::fs::File::open(&path).unwrap()));
    let standard_input = lengths(InputFile::as_standard_input(
        std::fs::File::open(&path).unwrap(),
    ));
    std::fs::remove_file(&path).unwrap();
    assert_eq!(named, vec![85]);
    assert_eq!(standard_input, vec![107]);
}
