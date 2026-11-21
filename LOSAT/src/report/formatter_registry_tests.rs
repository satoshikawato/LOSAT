use super::*;
use crate::blastinput::ncbi_registry::{read_registry, Layer};

#[test]
fn the_protocol_and_the_browser_url_read_ncbirc() {
    assert_eq!(protocol(None), b"https:");
    assert_eq!(
        tax_browser_template(None),
        b"<@protocol@>//www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi"
    );
    let read = read_registry(
        b"[BLASTFMTUTIL]\nPROTOCOL=http:\n[TAX_BROWSER]\nBLASTFMTUTIL=<@protocol@>//x\n",
        Layer::SystemWide,
    )
    .unwrap();
    assert_eq!(protocol(Some(&read)), b"http:");
    assert_eq!(tax_browser_template(Some(&read)), b"<@protocol@>//x");
    let read = read_registry(
        b"[BLASTFMTUTIL]\nPROTOCOL=\n[TAX_BROWSER]\nBLASTFMTUTIL=\n",
        Layer::SystemWide,
    )
    .unwrap();
    // An empty PROTOCOL is set (`HasEntry` with `fCountCleared`); an empty browser URL is not
    // (`Get(...).empty()`).
    assert_eq!(protocol(Some(&read)), b"");
    assert_eq!(value(Some(&read), "TAX_BROWSER", "BLASTFMTUTIL"), None);
}

#[test]
fn map_template_and_url_encode_are_ncbis() {
    assert_eq!(
        map_template(b"<@a@>x<@a@><@b@>", "a", b"<@a@>"),
        b"<@a@>x<@a@><@b@>"
    );
    assert_eq!(map_template(b"<@a@>x<@a@>", "a", b"1"), b"1x1");
    // `-db ../../../t2/db/ps_ids` (probe `bp.di.f0.rc_localid.h`), then the other classes.
    assert_eq!(
        url_encode(b"../../../t2/db/ps_ids"),
        b"..%2F..%2F..%2Ft2%2Fdb%2Fps_ids"
    );
    assert_eq!(url_encode(b"a b!$'()*,-._~"), b"a+b!$'()*,-._%7E");
    assert_eq!(
        url_encode(b"\"#%&+:;<=>?@[\\]^`{|}\x7f\x80\xff\x00"),
        b"%22%23%25%26%2B%3A%3B%3C%3D%3E%3F%40%5B%5C%5D%5E%60%7B%7C%7D%7F%80%FF%00"
    );
}

#[test]
fn init_config_finds_fmtcfg_then_ncbi() {
    use std::ffi::{OsStr, OsString};
    let files = ["/x/.ncbirc", "../ncbifmt/bad.rc", "bad.rc", "rel/.ncbirc"];
    let is_file =
        |path: &std::path::Path| files.iter().any(|file| path == std::path::Path::new(file));
    let find = |fmtcfg: Option<&str>, ncbi: Option<&str>| {
        init_config_path(fmtcfg.map(OsStr::new), ncbi.map(OsStr::new), &is_file)
    };
    let some = |path: &str| Some(OsString::from(path));
    // No `.ncbirc` in the directory: `$NCBI`, with `/` added once (probes `env_ncbi`,
    // `cfg.ncbibad`, `ncbibad_abs`).
    assert_eq!(find(None, Some("/x")), some("/x/.ncbirc"));
    assert_eq!(find(Some(""), Some("/x/")), some("/x/.ncbirc"));
    assert_eq!(find(None, Some("rel")), some("rel/.ncbirc"));
    assert_eq!(find(None, None), None);
    assert_eq!(find(None, Some("")), None);
    // `FMTCFG` in the directory first, else under `$NCBI` (probes `cfg.fmtbad`, `cfg.ncbifmt`).
    assert_eq!(find(Some("bad.rc"), Some("/x")), some("bad.rc"));
    assert_eq!(find(Some("other.rc"), Some("/x")), None);
    let in_ncbi = ["../ncbifmt/bad.rc"];
    let only_ncbi = |path: &std::path::Path| {
        in_ncbi
            .iter()
            .any(|file| path == std::path::Path::new(file))
    };
    assert_eq!(
        init_config_path(
            Some(OsStr::new("bad.rc")),
            Some(OsStr::new("../ncbifmt")),
            &only_ncbi
        ),
        some("../ncbifmt/bad.rc")
    );
}

#[test]
fn the_local_id_url_and_the_debug_line_are_ncbis() {
    // NCBI BLAST+ 2.17.0 `-html` with `T/t4/rc/localid/.ncbirc` (probes `bp.di.f0.rc_localid.h`,
    // `bp.s.f0.rc_localid.h`) and `GETURL_DEBUG=1` (`bp.di.f0.env_geturl.h`: one line).
    let ncbirc = b"[LOCAL_ID]\nTOOL_URL=http://local.example/<@seq_id@>?db=<@db_name@>&taxid=<@taxid@>&p=<@protocol@>&t=<@title@>\n";
    let config = InitConfig::of("blastp", Some((".ncbirc".into(), ncbirc.to_vec())), true).unwrap();
    let mut err = Vec::new();
    assert_eq!(
        config.local_id_url(b"s3_local_identifier", b"../../../t2/db/ps_ids", 9606, &mut err).unwrap(),
        b"http://local.example/s3_local_identifier?db=..%2F..%2F..%2Ft2%2Fdb%2Fps_ids&taxid=9606&p=https:&t=s3_local_identifier"
    );
    assert_eq!(
        config.local_id_url(b"Subject_1", b"", 0, &mut err).unwrap(),
        b"http://local.example/Subject_1?db=&taxid=0&p=https:&t=Subject_1"
    );
    assert_eq!(
        config.map_protocol(b"<@protocol@>//x", &mut err).unwrap(),
        b"https://x"
    );
    assert_eq!(err, b"REGISTRY: .ncbirc\n");
    // Without a file: no URL, no line (`env_geturl_none`); `PROTOCOL` of the file, an empty one
    // too, else `https:`.
    let none = InitConfig::of("blastp", None, true).unwrap();
    let mut err = Vec::new();
    assert_eq!(
        none.local_id_url(b"Subject_1", b"", 0, &mut err).unwrap(),
        b""
    );
    assert_eq!(
        none.map_protocol(b"<@protocol@>//x", &mut err).unwrap(),
        b"https://x"
    );
    assert!(err.is_empty());
    let empty = InitConfig::of(
        "blastp",
        Some(("f".into(), b"[blastfmtutil]\nprotocol=\n".to_vec())),
        false,
    )
    .unwrap();
    assert_eq!(
        empty.map_protocol(b"<@protocol@>//x", &mut err).unwrap(),
        b"//x"
    );
    assert!(err.is_empty());
    assert!(InitConfig::of(
        "blastp",
        Some(("f".into(), b"[NCBI]\n.Inherits=g\n".to_vec())),
        false
    )
    .is_err());
}

#[test]
fn a_malformed_init_config_file_throws_at_each_id() {
    // NCBI BLAST+ 2.17.0 `FMTCFG=bad.rc` (probe `cfg.fmtbad.L.h`, `T/t4/out3`): exit status 255
    // and this standard error, at the first id URL; with `-num_descriptions 0` each HSP throws
    // again (`cfg.fmtbad.nd0.L.h`).
    let config = InitConfig::of(
        "blastp",
        Some(("bad.rc".into(), b"this line has no equals sign\n".to_vec())),
        true,
    )
    .unwrap();
    let mut err = Vec::new();
    for _ in 0..2 {
        let error = config.map_protocol(b"<@protocol@>", &mut err).unwrap_err();
        let exception = crate::ncbi_exception::NcbiException::of_io_error(&error).unwrap();
        // NCBI's report in LOSAT's words (`crate::wording::losat_text`).
        assert_eq!(
            format!("Error: {}\n", String::from_utf8(exception.what()).unwrap()),
            crate::wording::losat_text(&format!(
                "Error: NCBI C++ Exception:\n    T0 \"{}/src/corelib/ncbireg.cpp\", line 722: Error: (CRegistryException::eEntry) ncbi::IRWRegistry::x_Read() - Invalid registry entry format: 'this line has no equals sign' (m_Pos = 1)\n\n",
                crate::ncbi_exception::NCBI_BUILD_ROOT
            ))
        );
    }
    assert!(config.local_id_url(b"Subject_1", b"", 0, &mut err).is_err());
    assert!(err.is_empty());
}

#[test]
fn a_registry_read_posts_as_ncbi() {
    // NCBI BLAST+ 2.17.0 (probe `T/t2/probe_rc`, `rep` and `bom`): the posts of one read; each
    // query's report reads `./.ncbirc` twice.
    let read = read_registry(b"[FOO]\na=1\na=2\n", Layer::SystemWide).unwrap();
    assert_eq!(
        registry_posts(&read, "blastp"),
        b"Warning: (110.8) [blastp] Found multiple [FOO] a settings; using the one from line 3\n"
            .to_vec()
    );
    let read = read_registry(b"\xef", Layer::SystemWide).unwrap();
    assert_eq!(
        registry_posts(&read, "blastp"),
        b"Error: (110.4) [blastp] Error reading the registry after line 1: \n".to_vec()
    );
}
