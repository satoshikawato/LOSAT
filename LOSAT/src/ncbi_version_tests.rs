use super::*;

#[test]
fn the_request_is_the_whole_word() {
    for request in [
        VersionRequest::Version,
        VersionRequest::Full,
        VersionRequest::FullXml,
        VersionRequest::FullJson,
    ] {
        assert_eq!(
            VersionRequest::of_word(request.word().as_bytes()),
            Some(request)
        );
    }
    for word in ["-version-full=1", "-version-fullx", "--version", "-Version"] {
        assert_eq!(VersionRequest::of_word(word.as_bytes()), None, "{word}");
    }
}

// NCBI BLAST+ 2.17.0 (the x64-linux release) `blastn -version-full`, with LOSAT's package and
// without the lines of NCBI's build (line RG item E3, `E/e1/ncbi/blastn.vfull.out`).
#[test]
fn the_full_version_is_ncbi_s_form_with_losat_s_values() {
    let package = format!(
        " Package: LOSAT {}, output compatible with blast 2.17.0\n",
        env!("CARGO_PKG_VERSION")
    );
    assert_eq!(
        print(VersionRequest::Version, "blastn", "2.17.0+"),
        format!("blastn: 2.17.0+\n{package}")
    );
    assert_eq!(
        print(VersionRequest::Full, "blastn", "2.17.0+"),
        format!("blastn: 2.17.0+\n{package} Package-Config:  \n GI-64bit:  TRUE\n")
    );
    let pre = env!("CARGO_PKG_VERSION_PRE");
    let (xml_name, json_name) = if pre.is_empty() {
        (String::new(), String::new())
    } else {
        (format!(" name=\"{pre}\""), format!(", \"name\": \"{pre}\""))
    };
    let (major, minor, patch) = (
        env!("CARGO_PKG_VERSION_MAJOR"),
        env!("CARGO_PKG_VERSION_MINOR"),
        env!("CARGO_PKG_VERSION_PATCH"),
    );
    assert_eq!(
        print(VersionRequest::FullXml, "rpsbproc", "0.5.1"),
        format!(
            "<?xml version=\"1.0\"?>
<ncbi_version xmlns=\"ncbi:version\"
  xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\"
  xs:schemaLocation=\"ncbi:version ncbi_version.xsd\">
<appname>rpsbproc</appname>
<version_info major=\"0\" minor=\"5\" patch_level=\"1\"/>
<package name=\"LOSAT\">
<version_info major=\"{major}\" minor=\"{minor}\" patch_level=\"{patch}\"{xml_name}/>
<build_info>
</build_info>
<config></config>
</package>
<build_info>
</build_info>
</ncbi_version>
"
        )
    );
    assert_eq!(
        print(VersionRequest::FullJson, "dustmasker", "1.0.0"),
        format!(
            "{{
  \"ncbi_version\": {{
    \"appname\": \"dustmasker\",
    \"version_info\": {{\"major\": 1, \"minor\": 0, \"patch_level\": 0}},
    \"component\": [],
    \"package\": {{
      \"name\": \"LOSAT\",
      \"version_info\": {{\"major\": {major}, \"minor\": {minor}, \"patch_level\": {patch}{json_name}}},
      \"build_info\": {{}},
      \"config\": \"\"}},
    \"build_info\": {{}}
  }}
}}
"
        )
    );
}

// The BLAST+ programs' version has the `+` of `CBlastVersion::Print` in the text only.
#[test]
fn the_numbers_leave_out_the_plus() {
    let numbers = VersionNumbers::of_print(crate::NCBI_BLAST_VERSION);
    assert_eq!(
        (numbers.major, numbers.minor, numbers.patch),
        ("2", "17", "0")
    );
    assert!(print(
        VersionRequest::FullXml,
        "makeblastdb",
        crate::NCBI_BLAST_VERSION
    )
    .contains("<version_info major=\"2\" minor=\"17\" patch_level=\"0\"/>\n<package"));
}
