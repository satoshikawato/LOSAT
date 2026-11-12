//! The CDD data loader: NCBI's number conversions and their messages, the source signatures and
//! pseudo PSSM-IDs, and the fixture data of `LOSAT/tests/fixtures/rpsbproc/data/`.

use super::*;

/// The text of a conversion error as rpsbproc reports it (the first of its two lines, after
/// `Error: (CStringException::eConvert) `).
fn message(e: NcbiException) -> String {
    let report = e.report_default_ex(106, 15, "Application's initialization failed");
    let first = report.lines().next().unwrap();
    first
        .strip_prefix("Error: (CStringException::eConvert) ")
        .unwrap()
        .to_string()
}

// PS16a probe 4: NCBI's rpsbproc on data files with numbers that do not convert printed
// `Cannot convert string '12x' to Uint8 (m_Pos = 2)` (a PSSM-ID of cddid.tbl), `'abc' to Uint8
// (m_Pos = 0)` (a length), `'abc' to double (m_Pos = 0)` (bitscore_specific.txt), `'x' to Int8
// (m_Pos = 0)` (an index of cddannot.dat; the cluster number after "cl" of
// family_superfamily_links).
#[test]
fn conversions_fail_with_ncbis_messages() {
    assert_eq!(
        message(string_to_uint(b"12x").unwrap_err()),
        "Cannot convert string '12x' to Uint8 (m_Pos = 2)"
    );
    assert_eq!(
        message(string_to_uint(b"abc").unwrap_err()),
        "Cannot convert string 'abc' to Uint8 (m_Pos = 0)"
    );
    assert_eq!(
        message(string_to_double(b"abc").unwrap_err()),
        "Cannot convert string 'abc' to double (m_Pos = 0)"
    );
    assert_eq!(
        message(string_to_int(b"x", 0).unwrap_err()),
        "Cannot convert string 'x' to Int8 (m_Pos = 0)"
    );
    // ncbistr.cpp:876-918, 646-655 (quoted at `ncbi_str::string_to_uint8`, `string_to_uint`).
    assert_eq!(string_to_uint(b"+42"), Ok(42));
    assert_eq!(
        message(string_to_uint(b"").unwrap_err()),
        "Cannot convert string '' to Uint8 (m_Pos = 0)"
    );
    assert_eq!(
        message(string_to_uint(b"+").unwrap_err()),
        "Cannot convert string '+' to Uint8 (m_Pos = 1)"
    );
    assert_eq!(
        message(string_to_uint(b"4294967296").unwrap_err()),
        "Cannot convert string '4294967296' to unsigned int, overflow (m_Pos = 0)"
    );
    assert_eq!(
        message(string_to_uint(b"18446744073709551616").unwrap_err()),
        "Cannot convert string '18446744073709551616' to Uint8 (m_Pos = 19)"
    );
    // ncbistr.cpp:798-873, 635-643 (quoted at `ncbi_str::string_to_int8`, `string_to_int`).
    assert_eq!(string_to_int(b"-17", 0), Ok(-17));
    assert_eq!(
        message(string_to_int(b"12 ", 0).unwrap_err()),
        "Cannot convert string '12 ' to Int8 (m_Pos = 2)"
    );
    assert_eq!(
        message(string_to_int(b"2147483648", 0).unwrap_err()),
        "Cannot convert string '2147483648' to int, overflow (m_Pos = 0)"
    );
    assert_eq!(
        message(string_to_int(b"9223372036854775808", 0).unwrap_err()),
        "Cannot convert string '9223372036854775808' to Int8, overflow (m_Pos = 18)"
    );
    // ncbistr.cpp:1312-1375, 1386-1403 (quoted at `ncbi_str::string_to_double`).
    assert_eq!(string_to_double(b"55.5"), Ok(55.5));
    assert_eq!(string_to_double(b"-1e3"), Ok(-1000.0));
    assert_eq!(
        message(string_to_double(b"12.5\r").unwrap_err()),
        "Cannot convert string '12.5\\r' to double (m_Pos = 4)"
    );
    assert_eq!(
        message(string_to_double(b"1e").unwrap_err()),
        "Cannot convert string '1e' to double (m_Pos = 1)"
    );
    // glibc's `strtod` reads the hexadecimal form; `0x` alone is the number 0 that ends before
    // the `x` (NCBI's rpsbproc read `0x10` as 16 and `0X1p3` as 8 and refused `0x` at position
    // 1, PS16a review oracle A4-1 cases `bs12`, `bs14`, `bs13`).
    assert_eq!(string_to_double(b"0x10"), Ok(16.0));
    assert_eq!(string_to_double(b"0X1p3"), Ok(8.0));
    assert_eq!(
        message(string_to_double(b"0x").unwrap_err()),
        "Cannot convert string '0x' to double (m_Pos = 1)"
    );
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

/// A copy of the fixture data in a new directory of the system's temporary directory.
fn fixture_copy(tag: &str) -> std::path::PathBuf {
    let fixture = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/rpsbproc/data");
    let dir = std::env::temp_dir().join(format!(
        "losat_rpsbproc_{tag}_{}_{}",
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
    dir
}

// offl_sparcle_data.cpp:217-308 (quoted at `load_cdtrack`): the shapes of `cdtrack.txt` lines.
// Whole lines load; a line of the accession alone reads `buf[npos]` (as 0) and then throws
// libstdc++'s `out_of_range` (oracle A4-1 case `ct3`); a line that ends after a field the loader
// skips completes after reading `buf[npos]` (case `ct1`, refused), and a line that ends in a
// space reads past its terminating 0 (refused).
#[test]
fn cdtrack_lines_read_as_ncbi() {
    let dir = fixture_copy("cdtrack");
    let load = |text: &str| {
        std::fs::write(dir.join("cdtrack.txt"), text).unwrap();
        CddData::load(&fixture_files(dir.to_str().unwrap()))
    };
    load(
        "# Acc ShortName PssmId Parent Root\n\
         cd14498 DSP 350348 N/A cd14498 x\n\
         cd14527 DSP_bac 350376 cd14498 cd14498 x\n\
         cd14504 DUSP23 cd14504 cd14498\n\
         unknown\n",
    )
    .unwrap();
    match load("cd14498\n") {
        Err(DataError::Std(what)) => assert_eq!(
            what,
            "basic_string::substr: __pos (which is 18446744073709551615) > this->size() (which is 7)"
        ),
        other => panic!("{other:?}"),
    }
    for refused in [
        "cd14498 DSP 350348\n",
        "cd14498 DSP 350348 N/A\n",
        "cd14498 \n",
    ] {
        assert!(
            matches!(load(refused), Err(DataError::Unsupported(_))),
            "{refused:?}"
        );
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// PD-LOSAT-NCBI-DEFECTS exception 9: PSSM-IDs of 2^31 and more are compared as signed values
// without the overflow of NCBI's `int` difference. The rows are those of the PS16a review's
// oracle A2 data `D2` (the fixture's `cddid.tbl` with eight PSSM-IDs near 2^31 among its domain
// rows and 2^31 at the end), where a wrapping difference misses rows; every id is found.
#[test]
fn pssm_ids_from_2_to_the_31_are_found() {
    assert_eq!(int_order(2_100_000_000, 3_000_000_000u32 as i32), 1);
    assert_eq!(int_order(i32::MIN, i32::MAX), -1);
    assert_eq!(int_order(-5, -5), 0);
    let dir = fixture_copy("bigids");
    let table = std::fs::read_to_string(dir.join("cddid.tbl")).unwrap();
    let (clusters, domains): (Vec<&str>, Vec<&str>) = table.lines().partition(|line| {
        line.split('\t')
            .nth(1)
            .is_some_and(|acc| acc.starts_with("cl"))
    });
    let probe = |id: u64, n: u32| format!("{id}\tpfam990{n:02}\tPROBE{n}\tprobe row\t100");
    let inserts = [
        (0, probe(2145012964, 1)),
        (8, probe(2146536812, 2)),
        (44, probe(2148697914, 3)),
        (45, probe(2147941452, 4)),
        (45, probe(2145856543, 5)),
        (52, probe(2146191475, 6)),
        (58, probe(2148253328, 7)),
        (83, probe(2148091740, 8)),
        (domains.len(), probe(2147483648, 90)),
    ];
    let mut lines: Vec<String> = clusters.iter().map(|line| line.to_string()).collect();
    for at in 0..=domains.len() {
        lines.extend(
            inserts
                .iter()
                .filter(|(after, _)| *after == at)
                .map(|(_, row)| row.clone()),
        );
        if let Some(line) = domains.get(at) {
            lines.push(line.to_string());
        }
    }
    std::fs::write(dir.join("cddid.tbl"), lines.join("\n") + "\n").unwrap();
    let data = CddData::load(&fixture_files(dir.to_str().unwrap())).unwrap();
    assert_eq!(data.domains.len(), domains.len() + inserts.len());
    for domain in &data.domains {
        let found = data.find_domain(domain.pssm_id).unwrap();
        assert_eq!(data.domains[found].pssm_id, domain.pssm_id);
    }
    for id in [488084, 2148253328, 2147483648] {
        assert!(data.find_domain(id).is_some(), "{id}");
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// offl_sparcle_data.cpp:158-170 (quoted at `load_cddids`): an accession of 1 or 0 bytes sets no
// flag (oracle A4-1 cases `acc1`, `acc0`: NCBI's output is that of the data without the row);
// a blank line, or a lone CR, has the length `''` (cases `wsline`, `crline`).
#[test]
fn short_accessions_and_blank_lines_load_as_ncbi() {
    let dir = fixture_copy("cddid");
    let table = std::fs::read_to_string(dir.join("cddid.tbl")).unwrap();
    for row in ["999999\tA\tx\ty\t100\n", "999999\t\tx\ty\t100\n"] {
        std::fs::write(dir.join("cddid.tbl"), format!("{table}{row}")).unwrap();
        let data = CddData::load(&fixture_files(dir.to_str().unwrap())).unwrap();
        let domain = &data.domains[data.find_domain(999999).unwrap()];
        assert!(!domain.curated && !domain.is_struct_dom);
    }
    for row in [" \n", "\r\n"] {
        std::fs::write(dir.join("cddid.tbl"), format!("{table}{row}")).unwrap();
        match CddData::load(&fixture_files(dir.to_str().unwrap())) {
            Err(DataError::Conversion(error)) => assert_eq!(
                message(error),
                "Cannot convert string '' to Uint8 (m_Pos = 0)"
            ),
            other => panic!("{other:?}"),
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

// offl_sparcle_data.cpp:479 and :586 (quoted at `load_sites`): a site line whose coordinates are
// its last field asks `reserve` for more than a vector holds (oracle A4-1 cases `annot10`,
// `gen10`).
#[test]
fn coordinates_at_the_end_of_a_site_line_throw_length_error() {
    let dir = fixture_copy("sites");
    for file in ["cddannot.dat", "cddannot_generic.dat"] {
        let original = std::fs::read(dir.join(file)).unwrap();
        let text = String::from_utf8_lossy(&original).into_owned();
        let row = text
            .lines()
            .find(|line| !line.is_empty() && !line.starts_with('#'))
            .unwrap();
        let fields: Vec<&str> = row.split('\t').collect();
        let cut = fields[..10].join("\t");
        std::fs::write(dir.join(file), format!("{cut}\n")).unwrap();
        match CddData::load(&fixture_files(dir.to_str().unwrap())) {
            Err(DataError::Std(what)) => assert_eq!(what, "vector::reserve"),
            other => panic!("{file}: {other:?}"),
        }
        std::fs::write(dir.join(file), &original).unwrap();
    }
    std::fs::remove_dir_all(&dir).unwrap();
}
