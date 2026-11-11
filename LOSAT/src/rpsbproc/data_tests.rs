//! The CDD data loader: NCBI's number conversions and their messages, the source signatures and
//! pseudo PSSM-IDs, and the fixture data of `LOSAT/tests/fixtures/rpsbproc/data/`.

use super::*;

// PS16a probe 4: NCBI's rpsbproc on data files with numbers that do not convert printed
// `Cannot convert string '12x' to Uint8 (m_Pos = 2)` (a PSSM-ID of cddid.tbl), `'abc' to Uint8
// (m_Pos = 0)` (a length), `'abc' to double (m_Pos = 0)` (bitscore_specific.txt), `'x' to Int8
// (m_Pos = 0)` (an index of cddannot.dat; the cluster number after "cl" of
// family_superfamily_links).
#[test]
fn conversions_fail_with_ncbis_messages() {
    assert_eq!(
        string_to_uint(b"12x").unwrap_err().message(),
        "Cannot convert string '12x' to Uint8 (m_Pos = 2)"
    );
    assert_eq!(
        string_to_uint(b"abc").unwrap_err().message(),
        "Cannot convert string 'abc' to Uint8 (m_Pos = 0)"
    );
    assert_eq!(
        string_to_double(b"abc").unwrap_err().message(),
        "Cannot convert string 'abc' to double (m_Pos = 0)"
    );
    assert_eq!(
        string_to_int(b"x").unwrap_err().message(),
        "Cannot convert string 'x' to Int8 (m_Pos = 0)"
    );
    // ncbistr.cpp:884-918, 646-655 (quoted at `string_to_uint`).
    assert_eq!(string_to_uint(b"+42"), Ok(42));
    assert_eq!(
        string_to_uint(b"").unwrap_err().message(),
        "Cannot convert string '' to Uint8 (m_Pos = 0)"
    );
    assert_eq!(
        string_to_uint(b"+").unwrap_err().message(),
        "Cannot convert string '+' to Uint8 (m_Pos = 1)"
    );
    assert_eq!(
        string_to_uint(b"4294967296").unwrap_err().message(),
        "Cannot convert string '4294967296' to unsigned int, overflow (m_Pos = 0)"
    );
    assert_eq!(
        string_to_uint(b"18446744073709551616")
            .unwrap_err()
            .message(),
        "Cannot convert string '18446744073709551616' to Uint8 (m_Pos = 19)"
    );
    // ncbistr.cpp:798-873, 635-643 (quoted at `string_to_int`).
    assert_eq!(string_to_int(b"-17"), Ok(-17));
    assert_eq!(
        string_to_int(b"12 ").unwrap_err().message(),
        "Cannot convert string '12 ' to Int8 (m_Pos = 2)"
    );
    assert_eq!(
        string_to_int(b"2147483648").unwrap_err().message(),
        "Cannot convert string '2147483648' to int, overflow (m_Pos = 0)"
    );
    assert_eq!(
        string_to_int(b"9223372036854775808").unwrap_err().message(),
        "Cannot convert string '9223372036854775808' to Int8, overflow (m_Pos = 18)"
    );
    // ncbistr.cpp:1309-1374 (quoted at `string_to_double`).
    assert_eq!(string_to_double(b"55.5"), Ok(Some(55.5)));
    assert_eq!(string_to_double(b"-1e3"), Ok(Some(-1000.0)));
    assert_eq!(
        string_to_double(b"12.5\r").unwrap_err().message(),
        "Cannot convert string '12.5\\r' to double (m_Pos = 4)"
    );
    assert_eq!(
        string_to_double(b"1e").unwrap_err().message(),
        "Cannot convert string '1e' to double (m_Pos = 1)"
    );
    assert_eq!(string_to_double(b"0x10"), Ok(None));
}

// ncbistr.cpp (`s_PrintableString`, quoted at `printable_string`).
#[test]
fn printable_string_escapes_as_ncbi() {
    assert_eq!(printable_string(b"a\tb"), "a\\tb");
    assert_eq!(printable_string(b"'\"\\"), "\\'\\\"\\\\");
    assert_eq!(printable_string(b"??x?"), "\\?\\?x?");
    assert_eq!(printable_string(b"\x01x"), "\\1x");
    assert_eq!(printable_string(b"\x017"), "\\0017");
    assert_eq!(printable_string(&[0xc3, 0xa9]), "\\303\\251");
}

// basedata.cpp:120-140, 168-193 (quoted at `dom_acc_type` and `compute_pseudo_pssm_id`).
#[test]
fn source_signatures_and_pseudo_pssm_ids() {
    assert_eq!(dom_acc_type(b"cd00001"), 0);
    assert_eq!(dom_acc_type(b"pfam00069"), 1);
    assert_eq!(dom_acc_type(b"PRK10941"), 6);
    assert_eq!(dom_acc_type(b"PTZ00001"), 11);
    assert_eq!(dom_acc_type(b"NF008188"), TOTALSIGS);
    // The signature must be shorter than the accession.
    assert_eq!(dom_acc_type(b"CD"), TOTALSIGS);
    assert_eq!(compute_pseudo_pssm_id(b""), 0);
    assert_eq!(
        compute_pseudo_pssm_id(b"pfam00069"),
        100_000_000 + 1_000_000 + 69
    );
    assert_eq!(
        compute_pseudo_pssm_id(b"NF008188"),
        100_000_000 + 12_000_000 + 8188
    );
    assert_eq!(
        compute_pseudo_pssm_id(b"x99999999999"),
        100_000_000u32
            .wrapping_add(12_000_000)
            .wrapping_add(99_999_999_999u64 as u32)
    );
}

fn fixture_files(dir: &str) -> DataFiles {
    let file = |name: &str| format!("{dir}/{name}");
    DataFiles {
        cddids: file("cddid.tbl"),
        spfeats: file("cddannot.dat"),
        genfeats: file("cddannot_generic.dat"),
        bscores: file("bitscore_specific.txt"),
        cdtrack: file("cdtrack.txt"),
        famlinks: file("family_superfamily_links"),
    }
}

#[test]
fn fixture_data_loads_and_finds_its_rows() {
    let dir = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/rpsbproc/data");
    let data = CddData::load(&fixture_files(dir)).unwrap();
    let cog = data.find_domain_by_accession(b"COG2453").unwrap();
    let domain = &data.domains[cog];
    assert_eq!(data.find_domain(domain.pssm_id), Some(cog));
    assert_eq!(data.text(domain.short_name), b"CDC14");
    assert!(!domain.curated);
    let cd = data.find_domain_by_accession(b"cd14527").unwrap();
    assert!(data.domains[cd].curated);
    assert!(data.domains[cd].min_bit_score > 0.0);
    let cluster = data.find_cluster(data.domains[cd].cluster_id).unwrap();
    assert_eq!(data.text(data.clusters[cluster].accession), b"cl28904");
    assert_eq!(data.clusters[cluster].pssm_id, 531719);
    assert!(data.find_domain(1).is_none());
    // A missing file is named as NCBI names it, in NCBI's order of loading.
    let missing = CddData::load(&fixture_files("/nonexistent")).unwrap_err();
    assert!(matches!(missing, DataError::Missing(path) if path == "/nonexistent/cddid.tbl"));
}

// offl_sparcle_data.cpp:217-308 (quoted at `load_cdtrack`): the parent and the root of a domain,
// and NCBI's root that is neither the parent nor the domain itself, which sets the parent.
#[test]
fn cdtrack_sets_parents_and_roots_as_ncbi() {
    let fixture = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/rpsbproc/data");
    let dir = std::env::temp_dir().join(format!(
        "losat_rpsbproc_cdtrack_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    for entry in std::fs::read_dir(fixture).unwrap() {
        let entry = entry.unwrap();
        std::fs::copy(entry.path(), dir.join(entry.file_name())).unwrap();
    }
    std::fs::write(
        dir.join("cdtrack.txt"),
        "# Acc ShortName PssmId Parent Root\n\
         cd14498 DSP 350348 N/A cd14498 x\n\
         cd14527 DSP_bac 350376 cd14498 cd14498 x\n\
         cd14504 DUSP23 350354 cd14504 cd14498 x\n",
    )
    .unwrap();
    let data = CddData::load(&fixture_files(dir.to_str().unwrap())).unwrap();
    let domain = |acc: &[u8]| &data.domains[data.find_domain_by_accession(acc).unwrap()];
    assert_eq!(
        (
            domain(b"cd14498").hierarchy_parent,
            domain(b"cd14498").hierarchy_root
        ),
        (0, 350348)
    );
    assert_eq!(
        (
            domain(b"cd14527").hierarchy_parent,
            domain(b"cd14527").hierarchy_root
        ),
        (350348, 350348)
    );
    // The root is neither the parent nor the domain: NCBI sets the parent and leaves the root.
    assert_eq!(
        (
            domain(b"cd14504").hierarchy_parent,
            domain(b"cd14504").hierarchy_root
        ),
        (350348, 0)
    );
    // A line that ends after a field the loader skips is read past its end by NCBI.
    std::fs::write(dir.join("cdtrack.txt"), "cd14498 DSP 350348\n").unwrap();
    assert!(matches!(
        CddData::load(&fixture_files(dir.to_str().unwrap())),
        Err(DataError::Unsupported(_))
    ));
    std::fs::remove_dir_all(&dir).unwrap();
}
