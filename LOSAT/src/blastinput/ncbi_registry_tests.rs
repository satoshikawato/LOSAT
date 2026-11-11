use super::*;

/// The entries of a read, as (section, name, value) texts.
fn entries(read: &RegistryRead) -> Vec<(&str, &str, &[u8])> {
    read.entries
        .iter()
        .map(|entry| {
            (
                entry.section.as_str(),
                entry.name.as_str(),
                entry.value.as_slice(),
            )
        })
        .collect()
}

fn unix(bytes: &[u8]) -> RegistryRead {
    read_registry_with(bytes, Layer::SystemWide, LineEnds::Unix).unwrap()
}

/// The message, the line of the throw, its kind and the position of the read's error.
fn failure(read: &RegistryRead) -> (String, u32, &'static str, String) {
    let error = read.error.as_ref().expect("an error");
    let what = error.what();
    let pos = what
        .rsplit_once("(m_Pos = ")
        .map(|(_, rest)| rest.trim_end_matches(")\n").to_string())
        .unwrap_or_default();
    (
        error.get_msg().to_string(),
        error.site().line,
        error.site().kind,
        pos,
    )
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:694-777 (quoted at `read_lines`). NCBI
// BLAST+ 2.17.0 (the official Linux x64 release, `blastn -version` with each `.ncbirc`, g_reg
// probes of 2026-11-11: cases noeq, nobracket, emptybr, spacesec, badsec, badentry, eqonly,
// midquote, badesc, hexov, trail2, cont3, sect): every error, its line and its position.
#[test]
fn every_error_is_ncbis() {
    let cases: [(&[u8], &str, u32, &str, &str); 13] = [
        (
            b"garbage line\n",
            "Invalid registry entry format: 'garbage line'",
            722,
            "CRegistryException::eEntry",
            "1",
        ),
        (
            b"[BLAST\nBLASTDB=/db\n",
            "Invalid registry section (']' is missing): `[BLAST'",
            698,
            "CRegistryException::eSection",
            "1",
        ),
        // `]` is looked for from the second byte after `[`.
        (
            b"[]\n",
            "Invalid registry section (']' is missing): `[]'",
            698,
            "CRegistryException::eSection",
            "1",
        ),
        (
            b"[ ]\n",
            "Unnamed registry section: `[ ]'",
            704,
            "CRegistryException::eSection",
            "1",
        ),
        (
            b"[a b]\n",
            "Invalid registry section name: `[a b]'",
            708,
            "CRegistryException::eSection",
            "1",
        ),
        (
            b"[BLAST]\nbad name=1\n",
            "Invalid registry entry name: 'bad name=1'",
            728,
            "CRegistryException::eEntry",
            "2",
        ),
        (
            b"[BLAST]\n=1\n",
            "Invalid registry entry name: '=1'",
            728,
            "CRegistryException::eEntry",
            "2",
        ),
        (
            b"[BLAST]\nX=a\"b\"c\n",
            "Single(unescaped) '\"' in the middle of registry value: 'X=a\"b\"c'",
            767,
            "CRegistryException::eValue",
            "2",
        ),
        (
            b"[BLAST]\nX=\\xZZ\n",
            "Badly placed '\\' in the registry value: 'X=\\xZZ'",
            777,
            "CRegistryException::eValue",
            "2",
        ),
        // `StringToUInt` of the hexadecimal digits overflows.
        (
            b"[BLAST]\nX=\\x100000000\n",
            "Badly placed '\\' in the registry value: 'X=\\x100000000'",
            777,
            "CRegistryException::eValue",
            "2",
        ),
        // A backslash at the end of the file: no line follows, the escape is unterminated.
        (
            b"[BLAST]\nX=a\\",
            "Badly placed '\\' in the registry value: 'X=a\\'",
            777,
            "CRegistryException::eValue",
            "2",
        ),
        // A continuation: the line count goes on, and the message joins the lines with the
        // letter `n`.
        (
            b"[BLAST]\nX=a\\\n b\"c\n",
            "Single(unescaped) '\"' in the middle of registry value: 'X=a\\n b\"c'",
            767,
            "CRegistryException::eValue",
            "3",
        ),
        (
            b"[BLAST] junk\n[BL]AST]\n[\tx\t]\n[x\n",
            "Invalid registry section (']' is missing): `[x'",
            698,
            "CRegistryException::eSection",
            "4",
        ),
    ];
    for (bytes, message, line, kind, pos) in cases {
        let read = unix(bytes);
        assert_eq!(
            failure(&read),
            (message.to_string(), line, kind, pos.to_string()),
            "{:?}",
            String::from_utf8_lossy(bytes)
        );
        assert_eq!(
            read.error.as_ref().unwrap().site().function,
            Some("ncbi::IRWRegistry::x_Read")
        );
        assert!(!read.stream_failed);
    }
}

// NCBI BLAST+ 2.17.0 (the official Linux x64 release), g_reg probe `noeq` (`.ncbirc` "garbage
// line"): the report that the Critical post quotes.
#[test]
fn the_report_is_ncbis() {
    let read = unix(b"garbage\n");
    assert_eq!(
        read.error.unwrap().what(),
        format!(
            "NCBI C++ Exception:\n    T0 \"{}/src/corelib/ncbireg.cpp\", line 722: Error: (CRegistryException::eEntry) ncbi::IRWRegistry::x_Read() - Invalid registry entry format: 'garbage' (m_Pos = 1)\n",
            crate::ncbi_exception::NCBI_BUILD_ROOT
        )
    );
}

// NCBI BLAST+ 2.17.0, g_reg probe `partial` (`[BLAST] BLASTDB` before a bad line: the search uses
// it) and `after` (a bad line first: nothing is read): the entries before an error stay.
#[test]
fn entries_before_an_error_stay() {
    let read = unix(b"[BLAST]\nBLASTDB=/db\ngarbage\nDATA_LOADERS=none\n");
    assert_eq!(entries(&read), [("BLAST", "BLASTDB", &b"/db"[..])]);
    assert_eq!(failure(&read).3, "3");
    let read = unix(b"garbage\n[BLAST]\nBLASTDB=/db\n");
    assert!(read.entries.is_empty());
    assert!(read.error.is_some());
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:652-680, 719-801 (quoted at
// `read_lines`): comments, blank lines and spaces; sections and names without case are kept as
// written; an entry before any section is dropped; a later entry comes after an earlier one (g_reg
// probes `multi`, `nosection`, `quoted`).
#[test]
fn entries_are_read_as_ncbi_reads_them() {
    let read = unix(
        b"NOSECTION=1\n# file comment\n  ; entry comment\n\n \t \n[ BLAST ]\n  BLASTDB =  /a  \n\
          blastdb=/b\n[ncbi]\ndata = \"/q\"\nx=y=z\n\x0b\x0cname.with-all_/=v\n",
    );
    assert_eq!(read.error, None);
    assert_eq!(
        entries(&read),
        [
            ("BLAST", "BLASTDB", &b"/a"[..]),
            ("BLAST", "blastdb", b"/b"),
            ("ncbi", "data", b"/q"),
            ("ncbi", "x", b"y=z"),
            ("ncbi", "name.with-all_/", b"v"),
        ]
    );
    // A section with text after `]`, and a name cut at the first `]`.
    let read = unix(b"[BLAST] junk\nA=1\n[BL]AST]\nB=2\n");
    assert_eq!(
        entries(&read),
        [("BLAST", "A", &b"1"[..]), ("BL", "B", b"2")]
    );
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:753-781 (quoted at `read_lines`): an
// unescaped `"` at either end is removed, an escaped one stays (g_reg probe `quote3`: no error).
#[test]
fn quotes_at_the_ends_are_removed() {
    for (raw, stored) in [
        (&b"none"[..], &b"none"[..]),
        (b"\"none\"", b"none"),
        (b"\"none", b"none"),
        (b"none\"", b"none"),
        (b"\"\"", b""),
        (b"\"", b""),
        (b"\"\"\"", b""),
        (b"", b""),
        (b"\"a\\\"b\"", b"a\"b"),
        (b"\\\"a\\\"", b"\"a\""),
    ] {
        let mut bytes = b"[S]\nN=".to_vec();
        bytes.extend_from_slice(raw);
        bytes.push(b'\n');
        let read = unix(&bytes);
        assert_eq!(read.error, None, "{raw:?}");
        assert_eq!(entries(&read), [("S", "N", stored)], "{raw:?}");
    }
    for raw in [&b"no\"ne"[..], b"\"a\"b\"", b"a\"\"b"] {
        let mut bytes = b"[S]\nN=".to_vec();
        bytes.extend_from_slice(raw);
        bytes.push(b'\n');
        assert_eq!(failure(&unix(&bytes)).1, 767, "{raw:?}");
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4789-4892 (`NStr::ParseEscapes`, quoted
// at `ncbi_str::parse_escapes`): the escapes of a value are read (g_reg probe `esc`: `\q` is no
// error).
#[test]
fn escapes_are_read() {
    let read = unix(
        b"[S]\nT=a\\tb\nH=\\x41\\x4a4bz\nO=\\101\\1012\\8\\q\nW=\\a\\b\\f\\n\\r\\v\\\\\n\
          Z=\\0\nBIG=\\777\n",
    );
    assert_eq!(read.error, None);
    assert_eq!(
        entries(&read),
        [
            ("S", "T", &b"a\tb"[..]),
            ("S", "H", b"A\x4bz"),
            ("S", "O", b"AA28q"),
            ("S", "W", b"\x07\x08\x0c\n\r\x0b\\"),
            ("S", "Z", b"\0"),
            ("S", "BIG", b"\xff"),
        ]
    );
    for raw in [&b"\\xZZ"[..], b"\\x", b"a\\", b"\\x100000000"] {
        let mut bytes = b"[S]\nN=".to_vec();
        bytes.extend_from_slice(raw);
        bytes.push(b'\n');
        assert_eq!(failure(&unix(&bytes)).1, 777, "{raw:?}");
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:740-751 (quoted at `read_lines`). NCBI
// BLAST+ 2.17.0, g_reg probe `cont` (`BLASTDB=<dir>\` then an empty line: the search path holds
// the directory and a line feed).
#[test]
fn a_backslash_at_the_end_continues_the_value() {
    let read = unix(b"[S]\nA=one \\\n   two  \\\n three\nB=x\\\n\nC=even\\\\\nD=1\n");
    assert_eq!(read.error, None);
    assert_eq!(
        entries(&read),
        [
            ("S", "A", &b"one \ntwo  \nthree"[..]),
            ("S", "B", b"x\n"),
            ("S", "C", b"even\\"),
            ("S", "D", b"1"),
        ]
    );
    // The continuation counts its lines.
    assert_eq!(failure(&unix(b"[S]\nA=1\\\n2\n\nbad\n")).3, "5");
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:196-208 (`NcbiGetlineEOL`, quoted at
// `Lines`). NCBI BLAST+ 2.17.0 for Linux, g_reg probe `cr`: a `\r` is part of the line there.
#[test]
fn line_ends_are_the_platforms() {
    let read = |bytes: &[u8], ends| read_registry_with(bytes, Layer::SystemWide, ends).unwrap();
    // A `\r` is a space at the ends of a name or value.
    let text = b"[S]\r\nA=1\r\nB=2\rC=3\n\rD=4";
    for ends in [LineEnds::Unix, LineEnds::Windows] {
        let got = read(text, ends);
        assert_eq!(got.error, None);
        assert_eq!(
            entries(&got),
            [
                ("S", "A", &b"1"[..]),
                ("S", "B", b"2\rC=3"),
                ("S", "D", b"4")
            ]
        );
    }
    let got = read(text, LineEnds::Darwin);
    assert_eq!(got.error, None);
    assert_eq!(
        entries(&got),
        [
            ("S", "A", &b"1"[..]),
            ("S", "B", b"2"),
            ("S", "C", b"3"),
            ("S", "D", b"4")
        ]
    );
    // The line of a message: Linux keeps every `\r`, Windows drops one before `\n`, macOS ends
    // the line at the first.
    let bad = b"garb\rage\r\n";
    let message = |ends| failure(&read(bad, ends)).0;
    assert_eq!(
        message(LineEnds::Unix),
        "Invalid registry entry format: 'garb\rage\r'"
    );
    assert_eq!(
        message(LineEnds::Windows),
        "Invalid registry entry format: 'garb\rage'"
    );
    assert_eq!(
        message(LineEnds::Darwin),
        "Invalid registry entry format: 'garb'"
    );
    // On macOS `\r\n` is one end and `\n\r` two.
    assert_eq!(failure(&read(b"[S]\r\nx\n", LineEnds::Darwin)).3, "2");
    assert_eq!(failure(&read(b"[S]\n\rx\n", LineEnds::Darwin)).3, "3");
    // A last line without an end is read.
    assert_eq!(
        entries(&read(b"[S]\nA=1", LineEnds::Unix)),
        [("S", "A", &b"1"[..])]
    );
    // The platform LOSAT is built for.
    let expected = if cfg!(windows) {
        LineEnds::Windows
    } else if cfg!(target_vendor = "apple") {
        LineEnds::Darwin
    } else {
        LineEnds::Unix
    };
    assert_eq!(LineEnds::HOST, expected);
    assert_eq!(
        read_registry(bad, Layer::SystemWide).unwrap(),
        read(bad, expected)
    );
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistre.cpp:782-826 (`GetTextEncodingForm`, quoted
// at `text_encoding_form`). NCBI BLAST+ 2.17.0, g_reg probes bom1-bom7, b1full, b1ini, inib1:
// the UTF-8 mark is dropped, other first bytes are read; a one-byte file of 0xEF, 0xFE or 0xFF is
// a failed stream, reported for `.ncbirc` and silent for `<program>.ini`; UTF-16 is not read by
// LOSAT; the application's file is checked twice.
#[test]
fn byte_order_marks_are_ncbis() {
    let ncbirc = |bytes: &[u8]| read_registry_with(bytes, Layer::SystemWide, LineEnds::Unix);
    let ini = |bytes: &[u8]| read_registry_with(bytes, Layer::Application, LineEnds::Unix);
    assert_eq!(
        failure(&ncbirc(b"\xEF\xBB\xBFgarbage\n").unwrap()).0,
        "Invalid registry entry format: 'garbage'"
    );
    // Any first byte of a mark before BB BF.
    assert_eq!(
        entries(&ncbirc(b"\xFE\xBB\xBF[S]\nA=1\n").unwrap()),
        [("S", "A", &b"1"[..])]
    );
    for bytes in [&b"\xEF"[..], b"\xFE", b"\xFF"] {
        let read = ncbirc(bytes).unwrap();
        assert!(read.stream_failed && read.error.is_none() && read.entries.is_empty());
        let read = ini(bytes).unwrap();
        assert!(!read.stream_failed && read.error.is_none() && read.entries.is_empty());
    }
    // The second check of the application's file.
    let read = ini(b"\xEF\xBB\xBF\xEF").unwrap();
    assert!(read.stream_failed);
    assert_eq!(
        failure(&ini(b"\xEF\xBB\xBF\xEF\xBB\xBFgarbage\n").unwrap()).0,
        "Invalid registry entry format: 'garbage'"
    );
    // `.ncbirc` drops one mark; the second is part of the line (U+FEFF in UTF-8).
    assert_eq!(
        failure(&ncbirc(b"\xEF\xBB\xBF\xEF\xBB\xBFgarbage\n").unwrap()).0,
        "Invalid registry entry format: '\u{feff}garbage'"
    );
    // A line that is not UTF-8 in a message is LOSAT's refusal (g_reg probes bom2, bom4: NCBI
    // quotes the bytes).
    assert!(ncbirc(b"\xEF\xBB").unwrap_err().contains("not UTF-8"));
    assert!(ncbirc(b"\xEFgarbage\n").unwrap_err().contains("not UTF-8"));
    for bytes in [&b"\xFF\xFE[\0x\0\n\0"[..], b"\xFE\xFF\0[", b"\xFF\xFE"] {
        assert!(ncbirc(bytes).unwrap_err().contains("UTF-16"), "{bytes:?}");
    }
    assert!(ini(b"\xEF\xBB\xBF\xFF\xFE").unwrap_err().contains("UTF-16"));
    // Empty files.
    for bytes in [&b""[..], b"\xEF\xBB\xBF"] {
        assert_eq!(ncbirc(bytes).unwrap(), unix(b""));
        assert_eq!(ini(bytes).unwrap(), unix(b""));
    }
    assert_eq!(
        error_reading_post(),
        "Error: (110.4) Error reading the registry after line 1: \n"
    );
}

// NCBI BLAST+ 2.17.0, g_reg probe `nul`: a NUL byte is part of the line and of the message.
#[test]
fn a_nul_byte_is_part_of_a_line() {
    let read = unix(b"[BLAST]\nX\0Y\n[S]\nA=b\0c\n");
    assert_eq!(failure(&read).0, "Invalid registry entry format: 'X\0Y'");
    let read = unix(b"[S]\nA=b\0c\n");
    assert_eq!(entries(&read), [("S", "A", &b"b\0c"[..])]);
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbireg.cpp:140-148 (`s_Backslashed`).
#[test]
fn backslashed_counts_the_backslashes_before() {
    assert!(!backslashed(b"", 0));
    assert!(!backslashed(b"a", 1));
    assert!(backslashed(b"\\", 1));
    assert!(!backslashed(b"\\\\", 2));
    assert!(backslashed(b"\\\\\\", 3));
    assert!(backslashed(b"a\\\"", 2));
    assert!(!backslashed(b"a\\\\\"", 3));
}

// NCBI BLAST+ 2.17.0, the official Linux x64 release (g_reg probes `noeq`, `nul`): the post of one
// load of a malformed `.ncbirc`, which ends where a NUL byte is; the last entry of a name counts.
#[test]
fn the_notice_and_the_last_value_are_ncbis() {
    let read = unix(b"garbage\n");
    assert_eq!(
        syntax_error_notice(read.error.as_ref().unwrap()),
        format!(
            "Critical: (110.6) CNcbiRegistry: Syntax error in system-wide configuration file: NCBI C++ Exception:\n    T0 \"{}/src/corelib/ncbireg.cpp\", line 722: Error: (CRegistryException::eEntry) ncbi::IRWRegistry::x_Read() - Invalid registry entry format: 'garbage' (m_Pos = 1)\n\n",
            crate::ncbi_exception::NCBI_BUILD_ROOT
        )
    );
    let read = unix(b"[BLAST]\nX\0Y\n");
    assert!(syntax_error_notice(read.error.as_ref().unwrap())
        .ends_with("Invalid registry entry format: 'X\n"));
    let read = unix(b"[BLAST]\nBLASTDB=/a\n[blast]\nblastdb=/b\n[NCBI]\nBLASTDB=/c\n");
    assert_eq!(read.last_value("Blast", "BlastDB"), Some(&b"/b"[..]));
    assert_eq!(read.last_value("NCBI", "BLASTDB"), Some(&b"/c"[..]));
    assert_eq!(read.last_value("BLAST", "DATA_LOADERS"), None);
}
