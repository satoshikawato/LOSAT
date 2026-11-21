use super::*;
use crate::ncbi_exception::NCBI_BUILD_ROOT;
use crate::utils::winmask::FormatSite;

// The texts of NCBI BLAST+ 2.17.0 (the official Linux x64 release) blastn
// `-window_masker_db` with hand-made counts files (line NP, M/STATE.md "M2 の実機確認"). The
// expected values are NCBI's texts in LOSAT's words (`crate::wording::losat_text`).
#[test]
fn the_reports_are_ncbis() {
    let root = NCBI_BUILD_ROOT;
    assert_eq!(
        String::from_utf8(WinMaskError::Open { name: "nosuch.obinary".into() }.discover_exception().unwrap().what()).unwrap(),
        crate::wording::losat_text(&format!("NCBI C++ Exception:\n    T0 \"{root}/src/algo/winmask/seq_masker_istat_factory.cpp\", line 171: Error: (Exception::open failed) ncbi::CSeqMaskerIstatFactory::DiscoverStatType() - could not open nosuch.obinary\n"))
    );
    let factory = format!("    T0 \"{root}/src/algo/winmask/seq_masker_istat_factory.cpp\", line 358: Error: (Exception::creation failure) ncbi::CSeqMaskerIstatFactory::create() - could not create a unit counts container\n");
    let cases = [
        (FormatSite::FileTooShort, "seq_masker_istat_obinary.cpp\", line 63: Error: (Exception::format error) file too short"),
        (FormatSite::UnitSize, "seq_masker_istat_obinary.cpp\", line 114: Error: (Exception::bad parameter) ncbi::CSeqMaskerIstatOBinary::CSeqMaskerIstatOBinary() - unit size must be in [1,16]"),
        (FormatSite::HashKeySize, "seq_masker_istat_obinary.cpp\", line 125: Error: (Exception::bad hash parameter) ncbi::CSeqMaskerIstatOBinary::CSeqMaskerIstatOBinary() - hash key size must be in [1,2*unit_size - 1]"),
        (FormatSite::HashOffset, "seq_masker_istat_obinary.cpp\", line 129: Error: (Exception::bad hash parameter) ncbi::CSeqMaskerIstatOBinary::CSeqMaskerIstatOBinary() - offset must by in [0,32 - hash_key_size]"),
        (FormatSite::HashShift, "seq_masker_istat_obinary.cpp\", line 133: Error: (Exception::bad hash parameter) ncbi::CSeqMaskerIstatOBinary::CSeqMaskerIstatOBinary() - shift must be in [1, 32 - hash_key_size]"),
        (FormatSite::HashTable, "seq_masker_istat_obinary.cpp\", line 197: Error: (Exception::format error) ncbi::CSeqMaskerIstatOBinary::CSeqMaskerIstatOBinary() - not enough data to fill the hash table"),
        (FormatSite::ValuesTable, "seq_masker_istat_obinary.cpp\", line 208: Error: (Exception::format error) ncbi::CSeqMaskerIstatOBinary::CSeqMaskerIstatOBinary() - not enough data to fill the values table"),
        (FormatSite::BinFileSize, "seq_masker_istat_bin.cpp\", line 79: Error: (Exception::file format error) ncbi::CSeqMaskerIstatBin::CSeqMaskerIstatBin() - wrong file size"),
        (FormatSite::BinUnitSize, "seq_masker_istat_bin.cpp\", line 99: Error: (Exception::file format error) ncbi::CSeqMaskerIstatBin::CSeqMaskerIstatBin() - illegal unit size"),
        (FormatSite::BinLength, "seq_masker_istat_bin.cpp\", line 108: Error: (Exception::file format error) ncbi::CSeqMaskerIstatBin::CSeqMaskerIstatBin() - wrong length"),
        // `-window_masker_taxid` with optimized ascii counts (line NP M3d).
        (FormatSite::OAsciiUnitSize, "seq_masker_istat_oascii.cpp\", line 101: Error: (Exception::bad parameter) ncbi::CSeqMaskerIstatOAscii::CSeqMaskerIstatOAscii() - unit size must be in [1,16]"),
        (FormatSite::OAsciiHashKeySize, "seq_masker_istat_oascii.cpp\", line 114: Error: (Exception::bad hash parameter) ncbi::CSeqMaskerIstatOAscii::CSeqMaskerIstatOAscii() - hash key size must be in [1,2*unit_size - 1]"),
        (FormatSite::OAsciiHashOffset, "seq_masker_istat_oascii.cpp\", line 119: Error: (Exception::bad hash parameter) ncbi::CSeqMaskerIstatOAscii::CSeqMaskerIstatOAscii() - offset must by in [0,32 - hash_key_size]"),
        (FormatSite::OAsciiHashShift, "seq_masker_istat_oascii.cpp\", line 124: Error: (Exception::bad hash parameter) ncbi::CSeqMaskerIstatOAscii::CSeqMaskerIstatOAscii() - shift must be in [1, 32 - hash_key_size]"),
        (FormatSite::OAsciiFileTooShort, "seq_masker_istat_oascii.cpp\", line 145: Error: (Exception::format error) ncbi::CSeqMaskerIstatOAscii::CSeqMaskerIstatOAscii() - file too short"),
        (FormatSite::OAsciiHashTable, "seq_masker_istat_oascii.cpp\", line 176: Error: (Exception::format error) ncbi::CSeqMaskerIstatOAscii::CSeqMaskerIstatOAscii() - not enough lines to fill the hash table"),
        (FormatSite::OAsciiValuesTable, "seq_masker_istat_oascii.cpp\", line 192: Error: (Exception::format error) ncbi::CSeqMaskerIstatOAscii::CSeqMaskerIstatOAscii() - not enough lines to fill the values table"),
    ];
    for (site, line) in cases {
        assert_eq!(
            String::from_utf8(WinMaskError::format(site).create_report().unwrap()).unwrap(),
            crate::wording::losat_text(&format!(
                "NCBI C++ Exception:\n    T0 \"{root}/src/algo/winmask/{line}\n{factory}"
            )),
            "{site:?}"
        );
    }
    let substr = "could not create a unit counts container: basic_string::substr: __pos (which is 7) > this->size() (which is 6)";
    assert_eq!(
        String::from_utf8(WinMaskError::CreateFail { message: substr.into() }.create_report().unwrap()).unwrap(),
        crate::wording::losat_text(&format!("NCBI C++ Exception:\n    T0 \"{root}/src/algo/winmask/seq_masker_istat_factory.cpp\", line 363: Error: (Exception::creation failure) ncbi::CSeqMaskerIstatFactory::create() - {substr}\n"))
    );
    assert_eq!(
        WinMaskError::BadIndex { key: 1, value: 2 }.create_report(),
        None
    );
}
