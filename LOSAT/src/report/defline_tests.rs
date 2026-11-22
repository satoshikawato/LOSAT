use super::*;
use crate::blastinput::fasta_reader::{read_all, FastaInputSource, ReaderConfig};

// NCBI BLAST+ 2.17.0 `blastp -db DB -outfmt "6 sseqid sacc saccver stitle"
// -num_threads 1` on databases that makeblastdb made (without -parse_seqids) from
// protein records with these titles (campaign dbsearch-20261012,
// scripts/db3_compare.py: cases shapes.titles.ids and shapes.punct*.ids): the three
// id columns are the first word of the title that `GenerateDefline` makes.
#[test]
fn database_subject_id_is_the_first_word_of_the_generated_title() {
    for (title, id) in [
        (
            "sp|P69905|HBA_HUMAN Hemoglobin subunit alpha OS=Homo sapiens OX=9606 GN=HBA1 PE=1 SV=2",
            "sp|P69905|HBA_HUMAN",
        ),
        (
            "gi|12345|ref|NP_000549.1| hemoglobin subunit alpha [Homo sapiens]",
            "gi|12345|ref|NP_000549.1|",
        ),
        ("lcl|local_id a local id as first word", "lcl|local_id"),
        ("Subject_1 a word like a local subject id", "Subject_1"),
        ("Query_1 a word like a local query id", "Query_1"),
        ("gnl|BL_ORD_ID|7 a word like an ordinal id", "gnl|BL_ORD_ID|7"),
        // The alignment title (flags 0) loses the TPA and TSA prefixes.
        ("TPA: a third party annotation prefix", "a"),
        ("TPA_exp: experimental prefix", "experimental"),
        ("TSA: transcriptome shotgun prefix", "transcriptome"),
        ("MAG: metagenome prefix", "metagenome"),
        ("UNVERIFIED: unverified prefix", "unverified"),
        ("tpa: lower case prefix", "lower"),
        (">leading greater-than sign", ">leading"),
        (". a period as first word", "."),
        ("- a hyphen as first word", "-"),
        ("(parenthesis first word)", "(parenthesis"),
        ("ampersand &amp; entity", "ampersand"),
        ("a,b comma inside", "a,b"),
        ("two  spaces  inside", "two"),
        ("word", "word"),
        // An empty title is `unnamed protein product`.
        ("", "unnamed"),
        (".", "."),
        (",", ","),
        ("~", "~"),
        (". .", "."),
        ("~ ~", "~"),
        (", x", ","),
        ("x ,", "x"),
        (",;,x y", ",;,x"),
    ] {
        assert_eq!(
            database_subject_id(title.as_bytes(), &BioseqFacts::default(), true),
            Ok(SubjectWord::Word(id.as_bytes().to_vec())),
            "{title:?}"
        );
    }
}

// The same oracle: a title of one space keeps the general id (NCBI writes
// `gnl|BL_ORD_ID|0` for sseqid and `BL_ORD_ID:0` for sacc and saccver), as do titles of
// spaces and a nucleotide's empty title (NCBI BLAST+ 2.17.0 `blastn`, `tblastn`, `blastp`,
// `blastx` on databases without parsed Seq-ids and taxids, line DBID ordname,
// `blastdb-ids-20261021/R/ordname/out/*.{nes,ner,pes,per}.x.f5.ncbi.out`), and NCBI dies of
// SIGSEGV on the titles `, ,` and `; ;` (exit status 139): LOSAT names them.
// A prefix alone (`MAG:`), and `MAG: (plastid) [ ]` of a protein with a source without a
// name (adjusted to `MAG:`), keep it too (probe of the review of line RG item T5c,
// `remaining-gaps-20261112/T/t5c/review/out/ncbi/sall.blastp.f6.*`: `gnl|BL_ORD_ID|7` and
// `BL_ORD_ID:7`); without a source that title gives `(plastid)`.
#[test]
fn database_subject_id_keeps_the_general_id_and_names_the_titles_it_does_not_write() {
    // So do `TPA:` and `TPA:TSA: `, whose prefixes the alignment title loses (review 1 of
    // ordname, `R/ordname/review1`).
    for (title, protein) in [
        (" ", true),
        ("  ", true),
        (" ", false),
        ("  ", false),
        ("", false),
        ("TPA:", true),
        ("TPA:", false),
        ("TPA:TSA: ", true),
        ("TPA:TSA: ", false),
        ("MAG:", true),
    ] {
        assert_eq!(
            database_subject_id(title.as_bytes(), &BioseqFacts::default(), protein),
            Ok(SubjectWord::General),
            "{title:?}"
        );
    }
    let nameless = BioseqFacts {
        has_source: true,
        ..BioseqFacts::default()
    };
    assert_eq!(
        database_subject_id(b"MAG: (plastid) [ ]", &nameless, true),
        Ok(SubjectWord::General)
    );
    assert_eq!(
        database_subject_id(b"MAG: (plastid) [ ]", &BioseqFacts::default(), true),
        Ok(SubjectWord::Word(b"(plastid)".to_vec()))
    );
    // A nucleotide's empty title without a name of its source.
    assert_eq!(
        database_subject_id(b"", &nameless, false),
        Ok(SubjectWord::General)
    );
    for title in [", ,", "; ;"] {
        assert!(
            database_subject_id(title.as_bytes(), &BioseqFacts::default(), true)
                .unwrap_err()
                .contains("x_CleanAndCompress reads past the end"),
            "{title:?}"
        );
    }
    assert_eq!(
        database_subject_id(b"caf\x81", &BioseqFacts::default(), true),
        Ok(SubjectWord::Unknown)
    );
}

/// The title that the FASTA reader gives a record with the defline `defline`.
fn reader_title(defline: &[u8]) -> Vec<u8> {
    let mut bytes = b">".to_vec();
    bytes.extend_from_slice(defline);
    bytes.extend_from_slice(b"\nACGT\n");
    let mut source =
        FastaInputSource::from_bytes(&bytes, ReaderConfig::subject("BLASTN", false, false));
    let records = read_all(&mut source, &mut |_| Ok(())).expect("one record");
    assert_eq!(records.len(), 1, "{defline:?}");
    records[0].title.clone()
}

/// The 68 bytes of a description row that hold the title (showdefline.cpp:915-930),
/// without the spaces that pad it.
fn description_field(label: &[u8]) -> Vec<u8> {
    let mut field = if label.len() > 68 {
        [&label[..65], b"..."].concat()
    } else {
        label.to_vec()
    };
    while field.last() == Some(&b' ') {
        field.pop();
    }
    field
}

type RpCase = (
    &'static str,
    &'static [u8],
    &'static [u8],
    Option<&'static [u8]>,
    &'static [u8],
    Option<&'static [u8]>,
    Option<&'static [u8]>,
    Option<&'static [u8]>,
);

/// NCBI BLAST+ 2.17.0 with one subject (`-subject`, default outfmt 0 and outfmt 6) of
/// each defline: the case, the defline, then for a nucleotide subject (BLASTN; TBLASTX
/// and TBLASTN write the same) and for a protein subject (BLASTP) the description field
/// (its padding removed), the alignment heading (`None`: NCBI writes `Sequence with id
/// ... alignment skipped` instead), then the outfmt 6 `sseqid` of the two (`None`: NCBI
/// stops with exit code 255). Inventory range RP (`scratch_RP/out/<case>`, extracted by
/// `s2_extract_rp.py` and `s2_gen_rust_cases.py` of session SFb).
#[rustfmt::skip]
const RP_CASES: &[RpCase] = &[
("t01_id_desc", b"id desc", b"id desc", Some(b"id desc"), b"id desc", Some(b"id desc"), Some(b"id"), Some(b"id")),
("t02_tab", b"id\x09desc", b"id", Some(b"id"), b"id", Some(b"id"), Some(b"id"), Some(b"id")),
("t03_lead_ws", b"  id desc", b"id desc", Some(b"id desc"), b"id desc", Some(b"id desc"), Some(b"id"), Some(b"id")),
("t04_empty", b"", b"", Some(b""), b"unnamed protein product", Some(b"unnamed protein product"), Some(b"Subject_1"), Some(b"unnamed")),
("t05_ws_only", b"   ", b"", Some(b""), b"unnamed protein product", Some(b"unnamed protein product"), Some(b"Subject_1"), Some(b"unnamed")),
("t06_e_acute", b"id\xc3\xa9x desc", b"id\xc3\xa9x desc", Some(b"id\xc3\xa9x desc"), b"id\xc3\xa9x desc", Some(b"id\xc3\xa9x desc"), Some(b"id\xc3\xa9x"), Some(b"id\xc3\xa9x")),
("t07_cjk", b"id\xe6\xbc\xa2 desc", b"id\xe6\xbc\xa2 desc", Some(b"id\xe6\xbc\xa2 desc"), b"id\xe6\xbc\xa2 desc", Some(b"id\xe6\xbc\xa2 desc"), Some(b"id\xe6\xbc\xa2"), Some(b"id\xe6\xbc\xa2")),
("t08_e9", b"id\xe9x desc", b"id\xc3\xa9x desc", Some(b"id\xc3\xa9x desc"), b"id\xc3\xa9x desc", Some(b"id\xc3\xa9x desc"), Some(b"id\xe9x"), Some(b"id\xe9x")),
("t09_ff", b"id\xffx desc \xff end", b"id\xc3\xbfx desc \xc3\xbf end", Some(b"id\xc3\xbfx desc \xc3\xbf end"), b"id\xc3\xbfx desc \xc3\xbf end", Some(b"id\xc3\xbfx desc \xc3\xbf end"), Some(b"id\xffx"), Some(b"id\xffx")),
("t11_pipes", b"a|b|c desc", b"a|b|c desc", Some(b"a|b|c desc"), b"a|b|c desc", Some(b"a|b|c desc"), Some(b"a|b|c"), Some(b"a|b|c")),
("t12_lcl", b"lcl|x desc", b"lcl|x desc", Some(b"lcl|x desc"), b"lcl|x desc", Some(b"lcl|x desc"), Some(b"lcl|x"), Some(b"lcl|x")),
("t13_gi", b"gi|123|gb|AB1.1| desc", b"gi|123|gb|AB1.1| desc", Some(b"gi|123|gb|AB1.1| desc"), b"gi|123|gb|AB1.1| desc", Some(b"gi|123|gb|AB1.1| desc"), Some(b"gi|123|gb|AB1.1|"), Some(b"gi|123|gb|AB1.1|")),
("t14_punct", b"|", b"|", Some(b"|"), b"|", Some(b"|"), Some(b"|"), Some(b"|")),
("t15_punct2", b"...", b"...", Some(b"..."), b"...", Some(b"..."), Some(b"..."), Some(b"...")),
("t16_amp", b"id &amp; &lt;b&gt; \"q\" desc", b"id & <b> \"q\" desc", Some(b"id & <b> \"q\" desc"), b"id & <b> \"q\" desc", Some(b"id & <b> \"q\" desc"), Some(b"id"), Some(b"id")),
("t17_same_word", b"same sdesc", b"same sdesc", Some(b"same sdesc"), b"same sdesc", Some(b"same sdesc"), Some(b"same"), Some(b"same")),
("t18_ctrl_lead", b"\x01id desc", b"\x01id desc", Some(b"\x01id desc"), b"\x01id desc", Some(b"\x01id desc"), Some(b"\x01id"), Some(b"\x01id")),
("t19_trail_ws", b"id desc   ", b"id desc", Some(b"id desc"), b"id desc", Some(b"id desc"), Some(b"id"), Some(b"id")),
("t20_idonly", b"id", b"id", Some(b"id"), b"id", Some(b"id"), Some(b"id"), Some(b"id")),
("t21_nbsp", b"id\xc2\xa0x desc", b"id\xc2\xa0x desc", Some(b"id\xc2\xa0x desc"), b"id\xc2\xa0x desc", Some(b"id\xc2\xa0x desc"), Some(b"id\xc2\xa0x"), Some(b"id\xc2\xa0x")),
("t22_dblspace", b"id  two  spaces", b"id two spaces", Some(b"id two spaces"), b"id two spaces", Some(b"id two spaces"), Some(b"id"), Some(b"id")),
("t23_q_only_title", b"", b"", Some(b""), b"unnamed protein product", Some(b"unnamed protein product"), Some(b"Subject_1"), Some(b"unnamed")),
("t24_s_only_title", b"sid sdesc", b"sid sdesc", Some(b"sid sdesc"), b"sid sdesc", Some(b"sid sdesc"), Some(b"sid"), Some(b"sid")),
("t25_dup_q", b"Subject_1 y", b"Subject_1 y", Some(b"Subject_1 y"), b"Subject_1 y", Some(b"Subject_1 y"), Some(b"Subject_1"), Some(b"Subject_1")),
("u01_end_eacute", b"id desc \xc3\xa9", b"id desc \xc3", Some(b"id desc \xc3"), b"id desc \xc3", Some(b"id desc \xc3"), Some(b"id"), Some(b"id")),
("u02_end_cjk", b"id \xe6\xbc\xa2", b"id \xe6\xbc", Some(b"id \xe6\xbc"), b"id \xe6\xbc", Some(b"id \xe6\xbc"), Some(b"id"), Some(b"id")),
("u03_end_ff", b"id desc \xff", b"id desc \xc3", Some(b"id desc \xc3"), b"id desc \xc3", Some(b"id desc \xc3"), Some(b"id"), Some(b"id")),
("u04_mixed", b"id \xc3\xa9 \xe9 x", b"id \xc3\x83\xc2\xa9 \xc3\xa9 x", Some(b"id \xc3\x83\xc2\xa9 \xc3\xa9 x"), b"id \xc3\x83\xc2\xa9 \xc3\xa9 x", Some(b"id \xc3\x83\xc2\xa9 \xc3\xa9 x"), Some(b"id"), Some(b"id")),
("u05_cp1252", b"id \x93quoted\x94 x", b"id \xe2\x80\x9cquoted\xe2\x80\x9d x", Some(b"id \xe2\x80\x9cquoted\xe2\x80\x9d x"), b"id \xe2\x80\x9cquoted\xe2\x80\x9d x", Some(b"id \xe2\x80\x9cquoted\xe2\x80\x9d x"), Some(b"id"), Some(b"id")),
("u06_undef81", b"id \x81 x", b"Unknown", None, b"Unknown", None, Some(b"id"), Some(b"id")),
("u07_euro80", b"id \x80 x", b"id \xe2\x82\xac x", Some(b"id \xe2\x82\xac x"), b"id \xe2\x82\xac x", Some(b"id \xe2\x82\xac x"), Some(b"id"), Some(b"id")),
("u08_entities", b"id &eacute; &#233; &#xe9; &#x4e2d; &lt;&gt;&quot;&apos; &nbsp;x end", b"id \xc3\xa9 \xc3\xa9 \xc3\xa9 \xe4\xb8\xad <>\"' \xc2\xa0x end", Some(b"id \xc3\xa9 \xc3\xa9 \xc3\xa9 \xe4\xb8\xad <>\"' \xc2\xa0x end"), b"id \xc3\xa9 \xc3\xa9 \xc3\xa9 \xe4\xb8\xad <>\"' \xc2\xa0x end", Some(b"id \xc3\xa9 \xc3\xa9 \xc3\xa9 \xe4\xb8\xad <>\"' \xc2\xa0x end"), Some(b"id"), Some(b"id")),
("u09_ent_zero", b"id a&#0;b end", b"id a\x00b end", Some(b"id a\x00b end"), b"id a\x00b end", Some(b"id a\x00b end"), Some(b"id"), Some(b"id")),
("u10_ent_big", b"id a&#1114112;b end", b"id a\xf4\x90\x80\x80b end", Some(b"id a\xf4\x90\x80\x80b end"), b"id a\xf4\x90\x80\x80b end", Some(b"id a\xf4\x90\x80\x80b end"), Some(b"id"), Some(b"id")),
("u11_ent_surr", b"id a&#xD800;b end", b"id a\xed\xa0\x80b end", Some(b"id a\xed\xa0\x80b end"), b"id a\xed\xa0\x80b end", Some(b"id a\xed\xa0\x80b end"), Some(b"id"), Some(b"id")),
("u12_cesu8", b"id \xed\xa0\xbd\xed\xb8\x80 end ok", b"id \xed\xa0\xbd\xed\xb8\x80 end ok", Some(b"id \xed\xa0\xbd\xed\xb8\x80 end ok"), b"id \xed\xa0\xbd\xed\xb8\x80 end ok", Some(b"id \xed\xa0\xbd\xed\xb8\x80 end ok"), Some(b"id"), Some(b"id")),
("u13_tpa", b"TPA: id desc", b"TPA: id desc", Some(b"id desc"), b"TPA: id desc", Some(b"id desc"), Some(b"TPA:"), Some(b"TPA:")),
("u14_unverified", b"UNVERIFIED: id desc", b"UNVERIFIED: id desc", Some(b"id desc"), b"UNVERIFIED: id desc", Some(b"id desc"), Some(b"UNVERIFIED:"), Some(b"UNVERIFIED:")),
("u15_mag", b"MAG x desc", b"MAG x desc", Some(b"x desc"), b"MAG x desc", Some(b"x desc"), Some(b"MAG"), Some(b"MAG")),
("u16_trail_punct", b"id desc.", b"id desc", Some(b"id desc"), b"id desc", Some(b"id desc"), Some(b"id"), Some(b"id")),
("u17_trail_punct2", b"id desc, ;~", b"id desc", Some(b"id desc"), b"id desc", Some(b"id desc"), Some(b"id"), Some(b"id")),
("u18_compress", b"id  a , b ;c ( d ) e,,f  g", b"id a, b; c (d) e, f g", Some(b"id a, b; c (d) e, f g"), b"id a, b; c (d) e, f g", Some(b"id a, b; c (d) e, f g"), Some(b"id"), Some(b"id")),
("u19_prot_brk", b"id desc. [Homo sapiens]", b"id desc. [Homo sapiens]", Some(b"id desc. [Homo sapiens]"), b"id desc [Homo sapiens]", Some(b"id desc [Homo sapiens]"), Some(b"id"), Some(b"id")),
("u20_prot_brk2", b"id desc, [Homo sapiens]", b"id desc, [Homo sapiens]", Some(b"id desc, [Homo sapiens]"), b"id desc [Homo sapiens]", Some(b"id desc [Homo sapiens]"), Some(b"id"), Some(b"id")),
("u21_del", b"id\x7fx desc", b"id\x7fx desc", Some(b"id\x7fx desc"), b"id\x7fx desc", Some(b"id\x7fx desc"), Some(b"id\x7fx"), Some(b"id\x7fx")),
("u22_ctrl_mid", b"id desc\x01more", b"id desc", Some(b"id desc"), b"id desc", Some(b"id desc"), Some(b"id"), Some(b"id")),
("u23_ent_trunc", b"id &amp", b"id &amp", Some(b"id &amp"), b"id &amp", Some(b"id &amp"), Some(b"id"), Some(b"id")),
("u24_amp_end", b"id desc &", b"id desc &", Some(b"id desc &"), b"id desc &", Some(b"id desc &"), Some(b"id"), Some(b"id")),
("u25_comma_only", b",", b",", Some(b","), b",", Some(b","), Some(b","), Some(b",")),
("u26_two_words_lead", b"\xc3\xa9 id", b"\xc3\xa9 id", Some(b"\xc3\xa9 id"), b"\xc3\xa9 id", Some(b"\xc3\xa9 id"), Some(b"\xc3\xa9"), Some(b"\xc3\xa9")),
("u27_len68", b"id xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx", b"id xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx", Some(b"id xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"), b"id xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx", Some(b"id xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"), Some(b"id"), Some(b"id")),
("u28_len69", b"id xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx", b"id xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx...", Some(b"id xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"), b"id xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx...", Some(b"id xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"), Some(b"id"), Some(b"id")),
("u29_ws_between", b"id \xc2\xa0 \xe3\x80\x80 x", b"id \xc2\xa0 \xe3\x80\x80 x", Some(b"id \xc2\xa0 \xe3\x80\x80 x"), b"id \xc2\xa0 \xe3\x80\x80 x", Some(b"id \xc2\xa0 \xe3\x80\x80 x"), Some(b"id"), Some(b"id")),
("u30_bom_lead", b"\xef\xbb\xbfid desc", b"\xef\xbb\xbfid desc", Some(b"\xef\xbb\xbfid desc"), b"\xef\xbb\xbfid desc", Some(b"\xef\xbb\xbfid desc"), Some(b"\xef\xbb\xbfid"), Some(b"\xef\xbb\xbfid")),
("v01_subj_ent", b"Subject_&amp;x desc", b"Subject_&x desc", Some(b"Subject_&x desc"), b"Subject_&x desc", Some(b"Subject_&x desc"), Some(b"Subject_&x"), Some(b"Subject_&x")),
("v02_subj_1", b"Subject_1 desc", b"Subject_1 desc", Some(b"Subject_1 desc"), b"Subject_1 desc", Some(b"Subject_1 desc"), Some(b"Subject_1"), Some(b"Subject_1")),
("v03_lcl_subj", b"lcl|Subject_7 foo", b"lcl|Subject_7 foo", Some(b"lcl|Subject_7 foo"), b"lcl|Subject_7 foo", Some(b"lcl|Subject_7 foo"), Some(b"lcl|Subject_7"), Some(b"lcl|Subject_7")),
("v04_subj_tpa", b"Subject_2 TPA: x", b"Subject_2 TPA: x", Some(b"Subject_2 TPA: x"), b"Subject_2 TPA: x", Some(b"Subject_2 TPA: x"), Some(b"Subject_2"), Some(b"Subject_2")),
("v05_subj_e9", b"Subject_\xe9 desc", b"Subject_\xc3\xa9 desc", Some(b"Subject_\xc3\xa9 desc"), b"Subject_\xc3\xa9 desc", Some(b"Subject_\xc3\xa9 desc"), Some(b"Subject_\xc3\xa9"), Some(b"Subject_\xc3\xa9")),
("v06_subj_81", b"Subject_\x81 desc", b"Unknown", None, b"Unknown", None, None, None),
("v07_pre_lcl", b"xlcl|Subject_ desc", b"xlcl|Subject_ desc", Some(b"xlcl|Subject_ desc"), b"xlcl|Subject_ desc", Some(b"xlcl|Subject_ desc"), Some(b"xlcl|Subject_"), Some(b"xlcl|Subject_")),
("v08_subj_trailpunct", b"Subject_3, desc", b"Subject_3, desc", Some(b"Subject_3, desc"), b"Subject_3, desc", Some(b"Subject_3, desc"), Some(b"Subject_3,"), Some(b"Subject_3,")),
("v09_subj_dotonly", b"Subject_4.", b"Subject_4", Some(b"Subject_4"), b"Subject_4", Some(b"Subject_4"), Some(b"Subject_4"), Some(b"Subject_4")),
("x01_tpa_only", b"TPA:", b"TPA:", Some(b""), b"TPA:", Some(b""), Some(b"TPA:"), Some(b"TPA:")),
("x02_tpa_sp", b"TPA: .", b"TPA:", Some(b""), b"TPA:", Some(b""), Some(b"TPA:"), Some(b"TPA:")),
("x03_unv", b"UNVERIFIED: ,;", b"UNVERIFIED:", Some(b""), b"UNVERIFIED:", Some(b""), Some(b"UNVERIFIED:"), Some(b"UNVERIFIED:")),
("x04_comma3", b",,,", b",", Some(b", "), b",", Some(b", "), Some(b",,,"), Some(b",,,")),
("x05_dot_sp", b". . x .", b". . x", Some(b". . x"), b". . x", Some(b". . x"), Some(b"."), Some(b".")),
("x06_tilde", b"~x~", b"~x", Some(b"~x"), b"~x", Some(b"~x"), Some(b"~x~"), Some(b"~x~")),
("x07_two_prefix", b"TPA: MAG x desc", b"TPA: MAG x desc", Some(b"MAG x desc"), b"TPA: MAG x desc", Some(b"MAG x desc"), Some(b"TPA:"), Some(b"TPA:")),
("x08_tpa_lower", b"tpa: lower", b"tpa: lower", Some(b"lower"), b"tpa: lower", Some(b"lower"), Some(b"tpa:"), Some(b"tpa:")),
("x09_trail_amp", b"x &amp;", b"x &amp", Some(b"x &amp"), b"x &amp", Some(b"x &amp"), Some(b"x"), Some(b"x")),
("x10_decimal_ent", b"x&#65;&#x42;&#67y z", b"xAB&#67y z", Some(b"xAB&#67y z"), b"xAB&#67y z", Some(b"xAB&#67y z"), Some(b"x&#65;&#x42;&#67y"), Some(b"x&#65;&#x42;&#67y")),
("x11_unterminated_ent", b"x &#65 z", b"x &#65 z", Some(b"x &#65 z"), b"x &#65 z", Some(b"x &#65 z"), Some(b"x"), Some(b"x")),
("x12_semicolon_end", b"x desc;;", b"x desc", Some(b"x desc"), b"x desc", Some(b"x desc"), Some(b"x"), Some(b"x")),
("y01_ovf", b"x &#4294967361; y", b"x A y", Some(b"x A y"), b"x A y", Some(b"x A y"), Some(b"x"), Some(b"x")),
("y02_dig16", b"x &#0000000000000065; y", b"x &#0000000000000065; y", Some(b"x &#0000000000000065; y"), b"x &#0000000000000065; y", Some(b"x &#0000000000000065; y"), Some(b"x"), Some(b"x")),
("y03_dig15", b"x &#000000000000065; y", b"x A y", Some(b"x A y"), b"x A y", Some(b"x A y"), Some(b"x"), Some(b"x")),
("y04_hexup", b"x &#X41; &#xg1; &#x; y", b"x A &#xg1; &#x; y", Some(b"x A &#xg1; &#x; y"), b"x A &#xg1; &#x; y", Some(b"x A &#xg1; &#x; y"), Some(b"x"), Some(b"x")),
("y05_namedcase", b"x &AMP; &Amp; &amp; &eacute; &Eacute; &notaname; y", b"x &AMP; &Amp; & \xc3\xa9 \xc3\x89 &notaname; y", Some(b"x &AMP; &Amp; & \xc3\xa9 \xc3\x89 &notaname; y"), b"x &AMP; &Amp; & \xc3\xa9 \xc3\x89 &notaname; y", Some(b"x &AMP; &Amp; & \xc3\xa9 \xc3\x89 &notaname; y"), Some(b"x"), Some(b"x")),
("y06_ent_in_prefix", b"TPA: &amp; x", b"TPA: & x", Some(b"& x"), b"TPA: & x", Some(b"& x"), Some(b"TPA:"), Some(b"TPA:")),
("y07_amp_amp", b"x &amp;amp; y", b"x &amp; y", Some(b"x &amp; y"), b"x &amp; y", Some(b"x &amp; y"), Some(b"x"), Some(b"x")),
("y08_nul_end", b"x &#0;", b"x &#0", Some(b"x &#0"), b"x &#0", Some(b"x &#0"), Some(b"x"), Some(b"x")),
("y09_ent_nonascii", b"x \xe9 &amp; y", b"x \xc3\xa9 & y", Some(b"x \xc3\xa9 & y"), b"x \xc3\xa9 & y", Some(b"x \xc3\xa9 & y"), Some(b"x"), Some(b"x")),
];

#[test]
fn titles_and_subject_ids_follow_the_ncbi_oracle() {
    for &(case, defline, n_desc, n_head, p_desc, p_head, n_id, p_id) in RP_CASES {
        let title = reader_title(defline);
        for (protein, description, heading, id) in
            [(false, n_desc, n_head, n_id), (true, p_desc, p_head, p_id)]
        {
            let shown = match generate_defline(&title, protein, true) {
                Ok(defline) => description_field(&defline.text),
                Err(UnknownEncoding) => b"Unknown".to_vec(),
            };
            assert_eq!(
                String::from_utf8_lossy(&shown),
                String::from_utf8_lossy(description),
                "{case} protein={protein} description"
            );
            let made = generate_defline(&title, protein, false)
                .ok()
                .map(|defline| defline.text);
            assert_eq!(made.as_deref(), heading, "{case} protein={protein} heading");
            let tabular = tabular_subject_id(b"Subject_1", &title, protein).ok();
            assert_eq!(tabular.as_deref(), id, "{case} protein={protein} sseqid");
        }
    }
}

#[test]
fn subject_ids_take_the_generated_word_where_lcl_subject_is_found() {
    // `lcl|Subject_` of `lcl|` and the word: at its start (the word starts with `Subject_`)
    // or inside the word, whose generated title drops the last period; nowhere else.
    for (title, id) in [
        ("Subject_4.", "Subject_4"),
        ("xlcl|Subject_4.", "xlcl|Subject_4"),
        ("lcl|Subject_4.", "lcl|Subject_4"),
        ("xSubject_4.", "xSubject_4."),
        ("Subject4.", "Subject4."),
        ("cl|Subject_4.", "cl|Subject_4."),
    ] {
        for protein in [false, true] {
            assert_eq!(
                tabular_subject_id(b"Subject_1", title.as_bytes(), protein).as_deref(),
                Ok(id.as_bytes()),
                "{title} protein={protein}"
            );
        }
    }
}

#[test]
fn titles_follow_ncbi_defline_generator() {
    // NCBI BLAST+ 2.17.0 outfmt 0 (alignment heading, description table).
    for (defline, heading, description) in [
        ("s1 abc.", "s1 abc", "s1 abc"),
        ("s1 abc ,def", "s1 abc, def", "s1 abc, def"),
        ("s1 a  b", "s1 a b", "s1 a b"),
        ("s1 (a )", "s1 (a)", "s1 (a)"),
        ("s1,,x", "s1, x", "s1, x"),
        ("TPA: s1 x", "s1 x", "TPA: s1 x"),
        ("s1 TPA: x", "s1 TPA: x", "s1 TPA: x"),
        ("s1 MAG x", "s1 MAG x", "s1 MAG x"),
        ("s1 x ; y", "s1 x; y", "s1 x; y"),
        ("s1 x,  y", "s1 x, y", "s1 x, y"),
        ("s1 x~", "s1 x", "s1 x"),
        ("s1 E. coli sp.", "s1 E. coli sp", "s1 E. coli sp"),
        ("s1 a, ,b", "s1 a, b", "s1 a, b"),
        ("s1 x;;y", "s1 x;;y", "s1 x;;y"),
        ("s1  lead", "s1 lead", "s1 lead"),
        ("s1 a ( b", "s1 a (b", "s1 a (b"),
        ("s1..x", "s1..x", "s1..x"),
        ("s1.x.", "s1.x", "s1.x"),
    ] {
        let made = |leave| generate_defline(defline.as_bytes(), false, leave).unwrap();
        assert_eq!(made(false).text, heading.as_bytes(), "{defline:?}");
        assert_eq!(made(true).text, description.as_bytes(), "{defline:?}");
    }
    let reads_past_end = |defline: &str| {
        [false, true].into_iter().any(|leave| {
            generate_defline(defline.as_bytes(), false, leave)
                .unwrap()
                .reads_past_end
        })
    };
    // NCBI BLAST+ 2.17.0 crashes on the outfmt 0 titles of the first group
    // (x_CleanAndCompress reads past the string); LOSAT stops at the end of the string
    // (approved exception 2 of PD-LOSAT-NCBI-DEFECTS). The second group is NCBI's.
    for (defline, title) in [
        (", ,", ", "),
        ("; ;", "; "),
        ("~, ,", "~, "),
        (",, ,", ",  "),
        (", ,,", ", "),
        (";  ;", "; "),
        (", , ,", ", "),
        (",~, ,", ",~, "),
        (", ;", ", ;"),
        (",,", ","),
        ("a, ,b", "a, b"),
        (", ,a", ", a"),
    ] {
        let made = generate_defline(defline.as_bytes(), false, false).unwrap();
        assert_eq!(made.text, title.as_bytes(), "{defline:?}");
    }
    for defline in [
        ", ,", "; ;", "~, ,", ",, ,", ", ,,", ";  ;", ", , ,", ",~, ,",
    ] {
        assert!(reads_past_end(defline), "{defline:?}");
    }
    for defline in [", ;", ",,", "a, ,b", ", ,a", "id, ,", "MAG: , ,x"] {
        assert!(!reads_past_end(defline), "{defline:?}");
    }
}

// Whether `NStr::HtmlDecode` decodes a character reference in the outfmt 0 title of a
// nucleotide subject with the defline `defline`, in the alignment heading or in the
// description table (its result flags `fHtmlDec_CharRef_Entity` or
// `fHtmlDec_CharRef_Numeric`, ncbistr.cpp:4580-4588).
fn ncbi_nucleotide_title_is_decoded(defline: &str) -> bool {
    [false, true].into_iter().any(|leave_prefix| {
        html_decode(&title_before_decoding_of(
            defline.as_bytes(),
            false,
            leave_prefix,
            &BioseqFacts::default(),
        ))
        .is_ok_and(|(_, char_refs)| char_refs)
    })
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4543-4590 (NStr::HtmlDecode): table
// entities and numeric references only, in the title after its final `.,;~ ` are
// trimmed (create_defline.cpp:3952-3960), with NCBI's start of an entity after `x`.
#[test]
fn titles_that_ncbi_decodes_are_found() {
    for defline in [
        "a&amp;b",
        "a&#38;b",
        "a&#x26;b",
        "a&#X26;b",
        "x &lt; y",
        "s &Tab;t",
        "a&#0;b",
        "TPA: &amp;x",
        "a &amp;; b",
    ] {
        assert!(ncbi_nucleotide_title_is_decoded(defline), "{defline}");
    }
    for defline in [
        "a & b",
        "a&b",
        "a&;",
        "a&#;",
        "a&amp b",
        "a&x#y;",
        "R&D; x",
        "a&foo;b",
        "a&amp;",
        "a&amp;;",
        "s &xi;t",
        "a&X41;b",
        "a&amp#;b",
        "a&ampxxxxxxxxxxxxxxxx;b",
        "a&#12345678901234567;b",
    ] {
        assert!(!ncbi_nucleotide_title_is_decoded(defline), "{defline}");
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:6711-6769 (CUtf8::GuessEncoding).
#[test]
fn encodings_are_guessed_as_ncbi_guesses_them() {
    for (text, encoding) in [
        (&b""[..], Encoding::Ascii),
        (b"id desc", Encoding::Ascii),
        (b"id \xc3\xa9", Encoding::Utf8),
        (b"\xe6\xbc\xa2", Encoding::Utf8),
        // The second byte after E0/ED/F0/F4 is not checked.
        (b"\xed\xa0\x80", Encoding::Utf8),
        (b"\xf4\x90\x80\x80", Encoding::Utf8),
        (b"id \xed\xa0\xbd\xed\xb8\x80 end ok", Encoding::Cesu8),
        (b"id \xe9 x", Encoding::Iso8859_1),
        (b"id \xc3\xa9 \xe9 x", Encoding::Iso8859_1),
        (b"id \xc3", Encoding::Iso8859_1),
        (b"\xc0\x80", Encoding::Windows1252),
        (b"id \x93q\x94", Encoding::Windows1252),
        (b"id \x80 x", Encoding::Windows1252),
        (b"id \x81 x", Encoding::Unknown),
        (b"\x8d\x8f\x90\x9d", Encoding::Unknown),
        (b"\xf5\x80\x80\x80", Encoding::Windows1252),
    ] {
        assert_eq!(guess_encoding(text), encoding, "{text:?}");
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4588-4606,7040-7060 (numeric references,
// `x_AppendChar` without a check); RP oracle cases u09-u11, y01-y04.
#[test]
fn numeric_references_are_written_without_a_check() {
    for (text, decoded) in [
        (&b"a&#0;b"[..], &b"a\x00b"[..]),
        (b"a&#xD800;b", b"a\xed\xa0\x80b"),
        (b"a&#1114112;b", b"a\xf4\x90\x80\x80b"),
        (b"x &#4294967361; y", b"x A y"),
        (b"x &#000000000000065; y", b"x A y"),
        (b"x &#0000000000000065; y", b"x &#0000000000000065; y"),
        (b"x &#X41; &#xg1; &#x; y", b"x A &#xg1; &#x; y"),
        (b"&#xFFFFFFFF;", b"\xff\xbf\xbf\xbf"),
        (
            b"&#x7FF;&#x800;&#xFFFF;&#x10000;",
            b"\xdf\xbf\xe0\xa0\x80\xef\xbf\xbf\xf0\x90\x80\x80",
        ),
        (b"\xe9&eacute;", b"\xc3\xa9\xc3\xa9"),
        (b"\x80&euro;", b"\xe2\x82\xac\xe2\x82\xac"),
        (b"&", b"&"),
        (b"a&#", b"a&#"),
    ] {
        assert_eq!(html_decode(text).unwrap().0, decoded, "{text:?}");
    }
    assert_eq!(html_decode(b"id \x81 &amp;"), Err(UnknownEncoding));
}

// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:302-311 (signed `curr`, the
// protein bracket rules); RP oracle cases u01-u03, u19, u20.
#[test]
fn the_last_byte_is_a_signed_char_and_protein_brackets_are_joined() {
    assert_eq!(clean_and_compress(b"id \xc3\xa9", false).0, b"id \xc3");
    assert_eq!(
        clean_and_compress(b"id \xe6\xbc\xa2", false).0,
        b"id \xe6\xbc"
    );
    assert_eq!(clean_and_compress(b"id a\x00", false).0, b"id a");
    assert_eq!(clean_and_compress(b"\xc3\xa9 x", false).0, b"\xc3\xa9 x");
    assert_eq!(
        clean_and_compress(b"id \xff. [x] y, [z]", true).0,
        b"id \xff [x] y [z]"
    );
    assert_eq!(clean_and_compress(b"a. . [b", true).0, b"a.  [b");
    assert_eq!(replace_in_place(b"x. [. [", b". [", b" ["), b"x [ [");
}

/// The facts of a database sequence without parsed Seq-ids whose source names `taxname`.
fn named(taxname: &[u8]) -> BioseqFacts {
    BioseqFacts::of_source(Some(std::sync::Arc::from(taxname)))
}

// NCBI BLAST+ 2.17.0 `-num_threads 1` on one-sequence databases without -parse_seqids, each
// without a title and with the taxid of one row of the probe's taxonomy database (line RG T5,
// `remaining-gaps-20261112/T/t5/out/ncbi`, `s.blastn.f0`, `s.blastp.f0`, `s.blastn.ids`,
// `c.*.long`): the outfmt 0 description (`fLeavePrefixSuffix`), the alignment heading (flags
// 0) and blastdbcmd's `-long_seqids` (`fLocalAnnotsOnly`) of the title made from the name.
#[test]
fn an_empty_title_is_made_from_the_name_of_the_source() {
    let nucleotide: [(&[u8], &str, &str); 14] = [
        (b"Homo sapiens", "Homo sapiens", "Homo sapiens"),
        (
            b"uncultured bacterium",
            "Uncultured bacterium",
            "Uncultured bacterium",
        ),
        (b"Bacillus sp.", "Bacillus sp.", "Bacillus sp."),
        (b"  Spaced   name  ", "Spaced name", "Spaced name"),
        (b"A &amp; B &lt;x&gt; &#65;", "A & B <x> A", "A & B <x> A"),
        (b"TPA: tpa organism", "TPA: tpa organism", "Tpa organism"),
        (
            b"MAG uncultured thing",
            "MAG uncultured thing",
            "Uncultured thing",
        ),
        (b"Name ending ;,~ ", "Name ending", "Name ending"),
        // One Latin-1 byte, decoded to UTF-8.
        (
            b"Caf\xe9 organism",
            "Caf\u{e9} organism",
            "Caf\u{e9} organism",
        ),
        (b"x", "X", "X"),
        (
            b"Odd ( spacing , here ) [strain]",
            "Odd (spacing, here) [strain]",
            "Odd (spacing, here) [strain]",
        ),
        (
            b"synthetic construct",
            "Synthetic construct",
            "Synthetic construct",
        ),
        (b".", ".", "."),
        (
            "\u{e9}toile utf8".as_bytes(),
            "\u{e9}toile utf8",
            "\u{e9}toile utf8",
        ),
    ];
    for (taxname, description, heading) in nucleotide {
        let made = |leave| {
            generate_bioseq_defline(b"", false, leave, &named(taxname))
                .unwrap()
                .text
        };
        assert_eq!(made(true), description.as_bytes(), "{taxname:?}");
        assert_eq!(made(false), heading.as_bytes(), "{taxname:?}");
    }
    for (taxname, description, heading) in [
        (
            "Homo sapiens",
            "unnamed protein product [Homo sapiens]",
            "unnamed protein product [Homo sapiens]",
        ),
        (
            "uncultured bacterium",
            "unnamed protein product [uncultured bacterium]",
            "unnamed protein product [uncultured bacterium]",
        ),
        (
            "  Spaced   name  ",
            "unnamed protein product [ Spaced name ]",
            "unnamed protein product [ Spaced name ]",
        ),
        (
            "Name ending ;,~ ",
            "unnamed protein product [Name ending;, ~ ]",
            "unnamed protein product [Name ending;, ~ ]",
        ),
        (
            "TPA: tpa organism",
            "unnamed protein product [TPA: tpa organism]",
            "unnamed protein product [TPA: tpa organism]",
        ),
        (
            "Odd ( spacing , here ) [strain]",
            "unnamed protein product [Odd (spacing, here) [strain]]",
            "unnamed protein product [Odd (spacing, here) [strain]]",
        ),
        (
            "synthetic construct",
            "unnamed protein product [synthetic construct]",
            "unnamed protein product [synthetic construct]",
        ),
    ] {
        let made = |leave| {
            generate_bioseq_defline(b"", true, leave, &named(taxname.as_bytes()))
                .unwrap()
                .text
        };
        assert_eq!(made(true), description.as_bytes(), "{taxname:?}");
        assert_eq!(made(false), heading.as_bytes(), "{taxname:?}");
    }
    // An unknown encoding throws (`Unknown`), as of a stored title.
    assert_eq!(
        generate_bioseq_defline(b"", false, true, &named(b"Bad \x81 name")),
        Err(UnknownEncoding)
    );
    // Without a name: an empty nucleotide title, `unnamed protein product`.
    assert!(
        generate_bioseq_defline(b"", false, false, &BioseqFacts::default())
            .unwrap()
            .text
            .is_empty()
    );
    // A stored title is not replaced, nor upper-cased.
    for leave in [false, true] {
        let made = generate_bioseq_defline(b"lower title", false, leave, &named(b"Homo sapiens"));
        assert_eq!(made.unwrap().text, b"lower title");
    }
    // blastdbcmd -long_seqids (`fLocalAnnotsOnly`): no `unnamed protein product`.
    let local = |taxname: Option<&[u8]>, protein| {
        let mut facts = BioseqFacts::of_source(taxname.map(std::sync::Arc::from));
        facts.local_annots_only = true;
        generate_bioseq_defline(b"", protein, false, &facts)
            .unwrap()
            .text
    };
    assert_eq!(local(Some(b"Homo sapiens"), true), b"[Homo sapiens]");
    assert_eq!(local(Some(b"  Spaced   name  "), true), b"[ Spaced name ]");
    assert_eq!(
        local(Some(b"uncultured bacterium"), false),
        b"Uncultured bacterium"
    );
    assert!(local(None, true).is_empty());
    // A `pdb` record without a chain is made from the source but not upper-cased.
    let mut pdb = named(b"x");
    pdb.pdb = Some((Vec::new(), 0));
    assert_eq!(
        generate_bioseq_defline(b"", false, false, &pdb)
            .unwrap()
            .text,
        b"x"
    );
    // `isprint ((unsigned char) m_PDBChain)` of an integer chain (create_defline.cpp:1842-1843):
    // its low byte, so 321 and -191 are `A` (65), and 383 (0x17F) is not printable.
    for (chain, text) in [
        (65, &b"Chain A"[..]),
        (321, b"Chain A"),
        (-191, b"Chain A"),
        (383, b"x"),
    ] {
        pdb.pdb = Some((Vec::new(), chain));
        let made = generate_bioseq_defline(b"", false, false, &pdb);
        assert_eq!(made.unwrap().text, text, "chain {chain}");
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:3197-3234 (`TruncateSpacesInPlace`): the C
// `isspace` bytes of both ends, none in the middle.
#[test]
fn a_nucleotide_source_title_loses_white_space_at_its_ends() {
    let facts = named(b"\t\n\x0b\x0c\r x \t y\r ");
    assert_eq!(source_title(false, &facts), b"x \t y");
    assert!(source_title(false, &named(b" \t ")).is_empty());
    assert_eq!(
        source_title(true, &named(b" x ")),
        b"unnamed protein product [ x ]"
    );
    let mut refseq = named(b" x ");
    refseq.refseq_chromosome_or_ncrna = true;
    assert_eq!(source_title(true, &refseq), b"x");
    refseq.local_annots_only = true;
    assert_eq!(source_title(true, &refseq), b"x");
    let mut mrna = named(b" x ");
    mrna.refseq_mrna = true;
    assert_eq!(source_title(true, &mrna), b"x");
    mrna.local_annots_only = true;
    assert_eq!(source_title(true, &mrna), b" [ x ]");
    assert!(source_title(
        true,
        &BioseqFacts {
            refseq_mrna: true,
            ..Default::default()
        }
    )
    .is_empty());
}

// A name that NCBI's `HtmlDecode` cannot decode (probe of line RG T5, `t.100009.*.f0`): the
// description is `Unknown`, and the heading names the sequence by its Seq-id text.
#[test]
fn a_name_that_does_not_decode_keeps_the_seq_id_of_the_heading() {
    let bad = SourceTaxname::Named(std::sync::Arc::from(&b"Bad \x81 name"[..]))
        .of_stored(|| b"BL_ORD_ID:8".to_vec());
    assert_eq!(bad.undecodable_id(), Some(&b"BL_ORD_ID:8"[..]));
    assert_eq!(bad.name().unwrap().as_deref(), Some(&b"Bad \x81 name"[..]));
    let good = SourceTaxname::Named(std::sync::Arc::from(&b"Homo sapiens"[..]))
        .of_stored(|| unreachable!("a name that decodes keeps no id"));
    assert_eq!(good.undecodable_id(), None);
    assert_eq!(
        SourceTaxname::Unnamed.of_stored(|| unreachable!()),
        SourceTaxname::Unnamed
    );
}

// NCBI BLAST+ 2.17.0 `blastp -num_threads 1`, the outfmt 0 heading of a line without a title
// (probe of line RG T5c, `remaining-gaps-20261112/T/t5c/out/ncbi/ppar.blastp.f0.{tx,nt}`),
// whose generated title is the Bioseq's (the titled line's) with its end rewritten from the
// first source: `Bacillus sp.` with the taxonomy database, no name without it.
#[test]
fn the_protein_suffix_is_rebuilt_from_the_source() {
    let adjusted = |title: &str, taxname: &str| {
        let mut title = title.as_bytes().to_vec();
        adjust_protein_title_suffix(&mut title, taxname.as_bytes());
        String::from_utf8(title).unwrap()
    };
    let bac = "Bacillus sp.";
    for (title, named, nameless) in [
        ("alpha, partial [Bacillus sp.]", "alpha [Bacillus sp.]", ""),
        ("beta [bacillus SP.]", "beta [Bacillus sp.]", ""),
        (
            "gamma (mitochondrion) [Bacillus sp.]",
            "gamma [Bacillus sp.]",
            "",
        ),
        ("delta (foo) [Bacillus sp.]", "", ""),
        ("epsilon [Homo sapiens]", "", ""),
        (
            "zeta protein, partial [ ]",
            "zeta protein [Bacillus sp.]",
            "zeta protein",
        ),
        ("eta protein, partial []", "", ""),
        ("theta. [Bacillus sp.]", "theta [Bacillus sp.]", ""),
        (
            "TPA: iota, partial [Bacillus sp.].",
            "TPA: iota [Bacillus sp.]",
            "",
        ),
        (
            "kappa, partial (chloroplast) [Bacillus sp.]",
            "kappa [Bacillus sp.]",
            "",
        ),
        (
            "lambda (chloroplast), partial [Bacillus sp.]",
            "lambda (chloroplast) [Bacillus sp.]",
            "",
        ),
        ("mu(plasmid) [Bacillus sp.]", "mu [Bacillus sp.]", ""),
        ("nu (x)(virus) [Bacillus sp.]", "", ""),
        ("nu (chromosome) [Bacillus sp.]", "", ""),
        ("nu (protoplast) [Bacillus sp.]", "nu [Bacillus sp.]", ""),
        ("x, [Bacillus sp.]", "x [Bacillus sp.]", ""),
        (
            "omicron, partial  [Bacillus sp.]",
            "omicron [Bacillus sp.]",
            "",
        ),
        ("pi, partial[Bacillus sp.]", "", ""),
        (
            "rho, partial, partial [Bacillus sp.]",
            "rho, partial [Bacillus sp.]",
            "",
        ),
        (
            "upsilon, partial [Bacillus sp.] ~;",
            "upsilon [Bacillus sp.]",
            "",
        ),
        (
            "ii protein (plastid) [ ]",
            "ii protein [Bacillus sp.]",
            "ii protein",
        ),
        (",, [Bacillus sp.]", ",,  [Bacillus sp.]", ""),
        // Not longer than ` [taxname]` and one letter: never adjusted.
        ("x [Bacillus sp.]", "", ""),
    ] {
        let trimmed = {
            let mut t = title.as_bytes().to_vec();
            trim_end_of(&mut t, b".,;~ ");
            String::from_utf8(t).unwrap()
        };
        let named = if named.is_empty() {
            trimmed.as_str()
        } else {
            named
        };
        let nameless = if nameless.is_empty() {
            trimmed.as_str()
        } else {
            nameless
        };
        assert_eq!(adjusted(title, bac), named, "{title:?} named");
        assert_eq!(adjusted(title, ""), nameless, "{title:?} without a name");
        assert!(protein_title_suffix_adjustable(title.as_bytes()));
    }
    // The other names of the probe's taxonomy database.
    assert_eq!(
        adjusted(
            "aa, partial [a &amp; b &lt;x&gt; &#65;]",
            "A &amp; B &lt;x&gt; &#65;"
        ),
        "aa [A &amp; B &lt;x&gt; &#65;]"
    );
    assert_eq!(
        adjusted("bb, partial [  spaced   NAME  ]", "  Spaced   name  "),
        "bb [  Spaced   name  ]"
    );
    // ` [ ]` needs the binomial's length: `a, [ ]` is adjusted, `a [ ]` is not.
    assert_eq!(adjusted("a, [ ]", ""), "a");
    assert_eq!(adjusted("a [ ]", ""), "a [ ]");
    assert_eq!(adjusted(", [ ]", ""), ", [ ]");
    assert_eq!(adjusted(", , [ ]", ""), ", , ");
    assert!(!protein_title_suffix_adjustable(b"alpha protein"));
    assert!(!protein_title_suffix_adjustable(b""));
    assert!(protein_title_suffix_adjustable(b"a] .~"));
    // In `GenerateDefline`: a protein heading with a source only, before the prefixes go and
    // the title is decoded; never a title made from the source.
    let source = |taxname: Option<&[u8]>| BioseqFacts {
        has_source: true,
        taxname: taxname.map(std::sync::Arc::from),
        ..BioseqFacts::default()
    };
    let made = |title: &[u8], protein, leave, facts: &BioseqFacts| {
        String::from_utf8(
            generate_bioseq_defline(title, protein, leave, facts)
                .unwrap()
                .text,
        )
        .unwrap()
    };
    let named = source(Some(b"Bacillus sp."));
    let title = b"MAG: sigma &amp; tau, partial [Bacillus sp.]";
    assert_eq!(
        made(title, true, false, &named),
        "sigma & tau [Bacillus sp.]"
    );
    assert_eq!(
        made(title, true, true, &named),
        "MAG: sigma & tau, partial [Bacillus sp.]"
    );
    assert_eq!(
        made(title, false, false, &named),
        "sigma & tau, partial [Bacillus sp.]"
    );
    assert_eq!(
        made(title, true, false, &BioseqFacts::default()),
        "sigma & tau, partial [Bacillus sp.]"
    );
    assert_eq!(made(b"chi, partial [ ]", true, false, &source(None)), "chi");
    assert_eq!(
        made(b"", true, false, &named),
        "unnamed protein product [Bacillus sp.]"
    );
    // A title that the adjustment leaves as punctuation is read past its end (NCBI BLAST+
    // 2.17.0 crashes: probe of line RG T5c, `pedge_e3`); a shorter one is not adjusted (NCBI
    // writes ` [ ]` for `, [ ]`: `pedge_e1`).
    let edge = generate_bioseq_defline(b", , [ ]", true, false, &source(None)).unwrap();
    assert!(edge.reads_past_end);
    assert!(
        !generate_bioseq_defline(b", , [ ]", true, false, &BioseqFacts::default())
            .unwrap()
            .reads_past_end
    );
    assert_eq!(made(b", [ ]", true, false, &source(None)), " [ ]");
}

// `title_ends_in_organism` compares the name at two places only; `NStr::Find`'s reverse
// search (the last match at or before `len1 - len2`, ASCII case folded) gives the same answer
// for every title and name of a small alphabet (each title of up to 8 bytes).
#[test]
fn the_organism_end_compares_two_places_as_the_reverse_search_answers() {
    let reverse_search = |title: &[u8], taxname: &[u8]| -> Option<usize> {
        let (len1, len2) = (title.len(), taxname.len());
        if len2 == 0 || len1 <= len2 + 4 {
            return None;
        }
        let idx = len1 - len2 - 3;
        if title[idx] != b' ' || title[idx + 1] != b'[' || title[len1 - 1] != b']' {
            return None;
        }
        let pos = (0..=len1 - len2)
            .rev()
            .find(|&at| title[at..at + len2].eq_ignore_ascii_case(taxname));
        (pos == Some(idx + 2)).then_some(idx + 1)
    };
    let alphabet = [b'a', b'A', b' ', b'[', b']'];
    let names: [&[u8]; 6] = [b"a", b"A", b"a]", b"[a", b" ", b"aA"];
    let mut title = Vec::new();
    for len in 0..=8u32 {
        for code in 0..alphabet.len().pow(len) {
            title.clear();
            let mut rest = code;
            for _ in 0..len {
                title.push(alphabet[rest % alphabet.len()]);
                rest /= alphabet.len();
            }
            for name in names {
                assert_eq!(
                    title_ends_in_organism(&title, name),
                    reverse_search(&title, name),
                    "{title:?} {name:?}"
                );
            }
        }
    }
}

// `protein_title_word_needs_source` (its derivation above it): where it says no, the first
// word of the generated title (`database_subject_id`) is the same without a source, with a
// source without a name and with any name, also the names that the title's ` [X]` ends match
// in either ASCII case and names that make the encoding fail. Titles: one thousand built
// from the parts NCBI's passes read (prefixes, organelles, `, partial`, entities, separators)
// and 60,000 drawn from an alphabet of such bytes (a fixed seed).
#[test]
fn a_title_whose_first_word_needs_no_source_gives_it_with_any_source() {
    let mut titles: Vec<Vec<u8>> = Vec::new();
    let words: [&[u8]; 12] = [
        b"hypothetical",
        b"alpha,",
        b"MAG",
        b"mag:",
        b"TPA:",
        b"a.b",
        b"x&amp;",
        b"caf\xc3\xa9",
        b"q\x81",
        b"(x)",
        b"WP_1|2",
        b"[x]",
    ];
    let middles: [&[u8]; 9] = [
        b"",
        b"protein 1",
        b", partial",
        b"(plastid)",
        b"x (mitochondrion), partial",
        b"&#44;y",
        b",z",
        b" ;",
        b"[y]",
    ];
    let ends: [&[u8]; 6] = [
        b"[Bacillus sp.]",
        b"[bacillus SP.]",
        b"[ ]",
        b"[]",
        b"[caf\x81]",
        b"[a [b]",
    ];
    for word in words {
        for middle in middles {
            for end in ends {
                for tail in [&b""[..], b".", b" ;~"] {
                    let mut title = word.to_vec();
                    title.push(b' ');
                    title.extend_from_slice(middle);
                    if !middle.is_empty() {
                        title.push(b' ');
                    }
                    title.extend_from_slice(end);
                    title.extend_from_slice(tail);
                    titles.push(title);
                }
            }
        }
    }
    let alphabet: &[u8] = b" [],.;~():&#aBM G1pX\x81\xc3\xa9";
    let mut seed: u64 = 0x5eed_1234_abcd_0001;
    let mut next = |bound: usize| {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((seed >> 33) as usize) % bound
    };
    for _ in 0..60_000 {
        let len = 1 + next(20);
        let mut title: Vec<u8> = (0..len).map(|_| alphabet[next(alphabet.len())]).collect();
        if next(5) > 0 {
            title.push(b']');
        }
        titles.push(title);
    }
    let mut safe = 0;
    for title in &titles {
        if protein_title_word_needs_source(title) {
            continue;
        }
        safe += 1;
        let main = match title.iter().rposition(|byte| !b".,;~ ".contains(byte)) {
            Some(end) => &title[..=end],
            None => &title[..],
        };
        let mut names: Vec<Vec<u8>> = vec![
            b"Bacillus sp.".to_vec(),
            b"caf\x81".to_vec(),
            b"\xc3\xa9".to_vec(),
            b" ".to_vec(),
            b"a".to_vec(),
        ];
        for at in 1..main.len() {
            if main[at - 1] == b' ' && main[at] == b'[' && main.len() > at + 2 {
                let inside = &main[at + 1..main.len() - 1];
                names.push(inside.to_vec());
                names.push(
                    inside
                        .iter()
                        .map(|byte| {
                            if byte.is_ascii_lowercase() {
                                byte.to_ascii_uppercase()
                            } else {
                                byte.to_ascii_lowercase()
                            }
                        })
                        .collect(),
                );
            }
        }
        let without = database_subject_id(title, &BioseqFacts::default(), true);
        let nameless = BioseqFacts {
            has_source: true,
            ..BioseqFacts::default()
        };
        assert_eq!(
            database_subject_id(title, &nameless, true),
            without,
            "{title:?} nameless"
        );
        for name in &names {
            assert_eq!(
                database_subject_id(title, &named(name), true),
                without,
                "{title:?} {name:?}"
            );
        }
    }
    assert!(safe > 500, "{safe}");
    for title in [
        &b"hypothetical protein 1 [Bacillus sp.]"[..],
        b"alpha [Bacillus sp.]",
        b"WP_1 [x]",
        b"x [y]",
    ] {
        assert!(!protein_title_word_needs_source(title), "{title:?}");
    }
    for title in [
        &b"alpha, partial [Bacillus sp.]"[..],
        b"MAG: (plastid) [ ]",
        b"MAG x [y]",
        b"chi [ ]",
        b" lead [x]",
        b"w (plastid) [x]",
        b"w &amp; [x]",
        b"w, [x]",
    ] {
        assert!(protein_title_word_needs_source(title), "{title:?}");
    }
    assert!(!protein_title_word_needs_source(b"no bracket"));
}
