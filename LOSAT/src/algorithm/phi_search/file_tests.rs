use super::*;

#[test]
fn last_pa_line_wins_and_keeps_its_cr() {
    let bytes = b"ID  one\nPA  [AG]-x(4)-G-K-[ST]\nPA  L-L\r\nCC  x\n";
    assert_eq!(read_pattern(bytes).unwrap(), b"L-L\r".to_vec());
}

#[test]
fn final_line_without_lf_is_read() {
    assert_eq!(read_pattern(b"PA  G-K").unwrap(), b"G-K".to_vec());
}

#[test]
fn bytes_after_a_nul_are_dropped() {
    assert_eq!(read_pattern(b"PA  G-K\0-S\n").unwrap(), b"G-K".to_vec());
    // A NUL inside the first two bytes leaves a line that is neither ID nor PA.
    assert_eq!(read_pattern(b"P\0  G\n"), Err(PatternFileError::NotRead));
}

#[test]
fn short_id_and_pa_lines_are_substr_errors() {
    assert_eq!(
        read_pattern(b"PA\n"),
        Err(PatternFileError::Substr { size: 2 })
    );
    assert_eq!(
        read_pattern(b"ID \nPA  G\n"),
        Err(PatternFileError::Substr { size: 3 })
    );
    // Four bytes: an empty pattern.
    assert_eq!(
        read_pattern(b"PA  G\nPA  \n"),
        Err(PatternFileError::NotRead)
    );
    // Lowercase and indented kinds are other lines.
    assert_eq!(
        read_pattern(b"pa  G\n PA G\n"),
        Err(PatternFileError::NotRead)
    );
}

#[test]
fn empty_file_and_no_pattern_are_not_read() {
    assert_eq!(read_pattern(b""), Err(PatternFileError::NotRead));
    assert_eq!(read_pattern(b"ID  x\n\n"), Err(PatternFileError::NotRead));
}

#[test]
fn a_line_of_the_buffer_size_ends_the_reading() {
    let mut long = b"PA  ".to_vec();
    long.resize(LINE_BUFFER - 1, b'A');
    // 4,095 bytes then LF: read.
    let mut bytes = long.clone();
    bytes.push(b'\n');
    assert_eq!(read_pattern(&bytes).unwrap().len(), LINE_BUFFER - 5);
    // 4,095 bytes at the end of the file: read.
    assert_eq!(read_pattern(&long).unwrap().len(), LINE_BUFFER - 5);
    // 4,096 bytes: the reading stops before the line, the earlier pattern stays.
    let mut bytes = b"PA  G-K\n".to_vec();
    bytes.extend_from_slice(&long);
    bytes.extend_from_slice(b"A\nPA  L\n");
    assert_eq!(read_pattern(&bytes).unwrap(), b"G-K".to_vec());
}
