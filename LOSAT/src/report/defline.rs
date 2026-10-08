//! The subject title of the pairwise report (outfmt 0) and the subject ID of the tabular
//! formats as NCBI makes them for a sequence read from FASTA:
//! `CDeflineGenerator::GenerateDefline` over the title of the record (the defline after
//! `>`, its first word included, as bytes), and `CShowBlastDefline::GetSeqIdList` over it.
//!
//! NCBI reference: c++/src/objmgr/util/create_defline.cpp:3952-3962
//! ```c
//!     if (! m_Reconstruct) {
//!         // x_SetFlags set m_MainTitle from a suitable descriptor, if any;
//!         // now strip trailing periods, commas, semicolons, and spaces.
//!         size_t pos = m_MainTitle.find_last_not_of (".,;~ ");
//!         if (pos != NPOS) {
//!             m_MainTitle.erase (pos + 1);
//!         }
//!         if (! m_MainTitle.empty()) {
//!             capitalize = false;
//!         }
//!     }
//! ```
//! NCBI reference: c++/src/objmgr/util/create_defline.cpp:4050-4073
//! ```c
//!     if (! (flags & fLeavePrefixSuffix)) {
//!         // remove TPA or TSA prefix, will rely on other data in record to set
//!         for (size_t i = 0; i < sizeof (s_tpaPrefixList) / sizeof (const char*); i++) {
//!             string str = s_tpaPrefixList [i];
//!             if (NStr::StartsWith (m_MainTitle, str, NStr::eNocase)) {
//!                 m_MainTitle.erase (0, str.length());
//!                 // strip leading spaces remaining after removal of old MAG before TPA or TSA prefixes
//!                 m_MainTitle.erase (0, m_MainTitle.find_first_not_of (' '));
//!             }
//!         }
//!     }
//!
//!     // strip leading spaces remaining after removal of old TPA or TSA prefixes
//!     m_MainTitle.erase (0, m_MainTitle.find_first_not_of (' '));
//!
//!     CStringUTF8 decoded = NStr::HtmlDecode (m_MainTitle);
//!
//!     // strip trailing commas, semicolons, and spaces (period may be an sp.
//!     // species)
//!     size_t pos = decoded.find_last_not_of (",;~ ");
//!     if (pos != NPOS) {
//!         decoded.erase (pos + 1);
//!     }
//! ```
//! NCBI reference: c++/src/objmgr/util/create_defline.cpp:4092-4101
//! ```c
//!     // produce final result
//!     string penult = mag + prefix + decoded + suffix;
//!
//!     x_CleanAndCompress (final, penult, m_IsAA);
//!
//!     if (! m_IsPDB && ! m_IsPatent && ! m_IsAA && ! m_IsSeg) {
//!         if (!final.empty() && islower ((unsigned char) final[0]) && capitalize) {
//!             final [0] = toupper ((unsigned char) final [0]);
//!         }
//!     }
//! ```
//! A local FASTA record has no MolInfo, source or metagenome data, so `mag`, `prefix` and
//! `suffix` are empty, and only an empty title is generated (`capitalize` stays true only
//! then): a protein's is `unnamed protein product` (`x_SetTitleFromProtein`), a nucleotide's
//! is empty, so nothing is capitalized. The alignment heading removes the prefixes; the
//! description table keeps them (`fLeavePrefixSuffix`, showdefline.cpp:498).
//! `NStr::HtmlDecode` guesses the encoding of the whole title (`CUtf8::GuessEncoding`) and
//! throws when it cannot; NCBI then shows `Unknown` and skips the alignments of the subject
//! (`UnknownEncoding`, see `report::pairwise`). For some titles of commas, semicolons,
//! tildes and spaces, NCBI's `x_CleanAndCompress` lets its count of the remaining letters
//! wrap and reads past the end of the string (and crashes); `clean_and_compress` stops at
//! the end of the string there (approved exception 2 of PD-LOSAT-NCBI-DEFECTS).

/// NCBI reference: c++/src/objmgr/util/create_defline.cpp:3431-3446
/// ```c
/// static const char* s_tpaPrefixList [] = {
///   "MAG ",
///   "MAG:",
///   "MULTISPECIES:",
///   "TLS:",
///   "TPA:",
///   "TPA_exp:",
///   "TPA_inf:",
///   "TPA_reasm:",
///   "TPA_asm:",
///   "TSA:",
///   "UNVERIFIED_ORG:",
///   "UNVERIFIED_ASMBLY:",
///   "UNVERIFIED_CONTAM:",
///   "UNVERIFIED:"
/// };
/// ```
const TPA_PREFIXES: [&str; 14] = [
    "MAG ",
    "MAG:",
    "MULTISPECIES:",
    "TLS:",
    "TPA:",
    "TPA_exp:",
    "TPA_inf:",
    "TPA_reasm:",
    "TPA_asm:",
    "TSA:",
    "UNVERIFIED_ORG:",
    "UNVERIFIED_ASMBLY:",
    "UNVERIFIED_CONTAM:",
    "UNVERIFIED:",
];

/// The title that NCBI generates for a protein record without one.
///
/// NCBI reference: c++/src/objmgr/util/create_defline.cpp:2508-2509
/// ```c
///     if (m_MainTitle.empty()  &&  !m_LocalAnnotsOnly) {
///         m_MainTitle = "unnamed protein product";
/// ```
/// (`x_SetTitleFromProtein`, called for an empty title of a protein at
/// create_defline.cpp:4020-4021.)
const UNNAMED_PROTEIN_PRODUCT: &[u8] = b"unnamed protein product";

/// `GenerateDefline` threw: `NStr::HtmlDecode` could not guess the encoding of the title
/// (`CUtf8::GuessEncoding` gave `eEncoding_Unknown`: a title that is not UTF-8 with one of
/// the bytes 0x81, 0x8D, 0x8F, 0x90 and 0x9D).
///
/// NCBI reference: c++/src/corelib/ncbistr.cpp:4528-4534
/// ```c
///     if (encoding == eEncoding_Unknown) {
///         encoding = CUtf8::GuessEncoding(str);
///         if (encoding == eEncoding_Unknown) {
///             NCBI_THROW2(CStringException, eBadArgs,
///                         "Unable to guess the source string encoding", 0);
///         }
///     }
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct UnknownEncoding;

/// A title made by `generate_defline`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Defline {
    /// The bytes NCBI writes.
    pub text: Vec<u8>,
    /// Whether NCBI's `x_CleanAndCompress` reads past the end of the title (NCBI crashes
    /// there; `text` stops at the end of the string).
    pub reads_past_end: bool,
}

/// `CDeflineGenerator::GenerateDefline` of a record from FASTA with the title `title` (see
/// the module): for the alignment heading (`leave_prefix_suffix` false, flags 0) or for
/// the description table (true, `fLeavePrefixSuffix`), of a protein (`m_IsAA`) or a
/// nucleotide record.
pub fn generate_defline(
    title: &[u8],
    protein: bool,
    leave_prefix_suffix: bool,
) -> Result<Defline, UnknownEncoding> {
    let (mut decoded, _) =
        html_decode(&title_before_decoding(title, protein, leave_prefix_suffix))?;
    trim_end_of(&mut decoded, b",;~ ");
    let (text, reads_past_end) = clean_and_compress(&decoded, protein);
    Ok(Defline {
        text,
        reads_past_end,
    })
}

/// The subject ID of the tabular formats (`sseqid`, `sacc`, `saccver`) for a subject with
/// the local ID `local_id` (`Subject_N`) and the title `title`: the title's first word
/// (`s_ReplaceLocalId`; the local ID for an empty title), and when that names a local
/// `Subject_` ID, the first word of the title that `GenerateDefline` makes (an empty
/// protein title gives `unnamed`). `Err` where `GenerateDefline` throws (NCBI stops with
/// an exception; RP-20).
///
/// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:787-789
/// ```c
///         if((m_IsBl2Seq && (!m_BelieveQuery))|| m_IsRemoteSearch) {
///         	tabinfo.SetParseSubjectDefline(true);
///         }
/// ```
/// NCBI reference: c++/src/objtools/align_format/tabular.cpp:522-532
/// ```c
/// void CBlastTabularInfo::SetSubjectId(const CBioseq_Handle& bh)
/// {
///     m_SubjectId.clear();
///
///     vector<CConstRef<objects::CSeq_id> > subject_id_list;
///     ITERATE(CBioseq_Handle::TId, itr, bh.GetId()) {
///     	CRef<CSeq_id> next_id = s_ReplaceLocalId(bh, itr->GetSeqId(), !m_ParseSubjectDefline );
///     	subject_id_list.push_back(next_id);
///     }
///     CShowBlastDefline::GetSeqIdList(bh, subject_id_list, m_SubjectId);
/// }
/// ```
/// NCBI reference: c++/src/objtools/align_format/tabular.cpp:482-496
/// ```c
///         if (sid_in->IsLocal()) {
///             string id_token;
///             vector<string> title_tokens;
///             title_tokens =
///                 NStr::Split(CAlignFormatUtil::GetTitle(bh), " ", title_tokens);
///             if(title_tokens.empty()){
///                 id_token = NcbiEmptyString;
///             } else {
///                 id_token = title_tokens[0];
///             }
///
///             if (id_token == NcbiEmptyString || parse_local) {
///                 const CObject_id& obj_id = sid_in->GetLocal();
///                 if (obj_id.IsStr())
///                     id_token = obj_id.GetStr();
/// ```
/// NCBI reference: c++/src/objtools/align_format/showdefline.cpp:218-242
/// ```c
///     ITERATE(vector< CConstRef<CSeq_id> >, itr, original_seqids) {
///         CRef<CSeq_id> next_seqid(new CSeq_id());
///         string id_token = NcbiEmptyString;
///
///         if (((*itr)->IsGeneral() &&
///             (*itr)->AsFastaString().find("gnl|BL_ORD_ID")
///             != string::npos) ||
/// 		(*itr)->AsFastaString().find("lcl|Subject_") != string::npos) {
///             vector<string> title_tokens;
///             string defline = sequence::CDeflineGenerator().GenerateDefline(bh);
///             if (defline != NcbiEmptyString) {
///                 id_token =
///                     NStr::Split(defline, " ", title_tokens)[0];
///             }
///         }
///         if (id_token != NcbiEmptyString) {
///             // Create a new local id with a label containing the extracted
///             // token and save it in the next_seqid instead of the original
///             // id.
///             CObject_id* obj_id = new CObject_id();
///             obj_id->SetStr(id_token);
///             next_seqid->SetLocal(*obj_id);
///         } else {
///             next_seqid->Assign(**itr);
///         }
/// ```
/// A local ID is written `lcl|` and its text (`CSeq_id::WriteAsFasta`, Seq_id.cpp:2164-2196,
/// `CObject_id::AsString`, Object_id.cpp:202-210).
pub fn tabular_subject_id(
    local_id: &[u8],
    title: &[u8],
    protein: bool,
) -> Result<Vec<u8>, UnknownEncoding> {
    let replaced = first_word(title);
    let replaced = if replaced.is_empty() {
        local_id
    } else {
        replaced
    };
    let mut fasta = b"lcl|".to_vec();
    fasta.extend_from_slice(replaced);
    if contains(&fasta, b"lcl|Subject_") {
        let defline = generate_defline(title, protein, false)?.text;
        let token = first_word(&defline);
        if !token.is_empty() {
            return Ok(token.to_vec());
        }
    }
    Ok(replaced.to_vec())
}

/// `NStr::Split(text, " ", tokens)[0]`: the bytes before the first space.
fn first_word(text: &[u8]) -> &[u8] {
    let end = text
        .iter()
        .position(|&byte| byte == b' ')
        .unwrap_or(text.len());
    &text[..end]
}

fn contains(text: &[u8], pattern: &[u8]) -> bool {
    text.windows(pattern.len()).any(|window| window == pattern)
}

/// Rejects the outfmt 0 title of the protein subject record `record` (from 1) of
/// `program` that LOSAT does not write as NCBI: a title that NCBI decodes
/// (`NStr::HtmlDecode`), or one that NCBI's `x_CleanAndCompress` reads past (NCBI crashes;
/// approved exception 2 of PD-LOSAT-NCBI-DEFECTS covers the nucleotide subjects of
/// BLASTN, TBLASTX and TBLASTN only).
pub fn check_shown_protein_subject_title(
    defline: &str,
    record: usize,
    program: &str,
) -> anyhow::Result<()> {
    check_shown_subject_title(defline, record, program)?;
    let reads_past_end = [false, true].into_iter().any(|leave_prefix| {
        generate_defline(defline.as_bytes(), true, leave_prefix)
            .is_ok_and(|title| title.reads_past_end)
    });
    if reads_past_end {
        anyhow::bail!(
            "subject record {record} has a title of punctuation that NCBI BLAST+'s x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 crashes in outfmt 0); this is not supported by LOSAT's {program}"
        );
    }
    Ok(())
}

/// Whether `NStr::HtmlDecode` decodes a character reference in the outfmt 0 title of a
/// nucleotide subject with the defline `defline`, in the alignment heading or in the
/// description table (its result flags `fHtmlDec_CharRef_Entity` or
/// `fHtmlDec_CharRef_Numeric`, ncbistr.cpp:4580-4588). The callers that still read with
/// `bio` reject such a subject where NCBI writes its title.
pub fn ncbi_nucleotide_title_is_decoded(defline: &str) -> bool {
    [false, true].into_iter().any(|leave_prefix| {
        html_decode(&title_before_decoding(
            defline.as_bytes(),
            false,
            leave_prefix,
        ))
        .is_ok_and(|(_, char_refs)| char_refs)
    })
}

/// Rejects the outfmt 0 title of the subject record `record` (from 1) of `program` that
/// the callers reading with `bio` do not write as NCBI: a title that NCBI decodes
/// (`NStr::HtmlDecode`). NCBI makes the titles only of the subjects that a report shows,
/// those with hits, so the callers check those subjects after the search. A title that
/// NCBI's `x_CleanAndCompress` reads past (NCBI crashes) is written as `generate_defline`
/// stops it at the end of the string (approved exception 2 of PD-LOSAT-NCBI-DEFECTS,
/// extended to TBLASTX and TBLASTN by the maintainer in session S08b, DW-17).
///
/// NCBI reference: c++/src/objmgr/util/create_defline.cpp:4092-4095
/// ```c
///     // produce final result
///     string penult = mag + prefix + decoded + suffix;
///
///     x_CleanAndCompress (final, penult, m_IsAA);
/// ```
///
/// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1540
/// ```c
///         x_DisplayDeflines(aln_set, itr_num, prev_seqids);
/// ```
/// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:1550-1551
/// ```c
///     CSeq_align_set copy_aln_set;
///     CBlastFormatUtil::PruneSeqalign(*aln_set, copy_aln_set, m_NumAlignments);
/// ```
pub fn check_shown_subject_title(
    defline: &str,
    record: usize,
    program: &str,
) -> anyhow::Result<()> {
    if ncbi_nucleotide_title_is_decoded(defline) {
        anyhow::bail!(
            "subject record {record} has an HTML character reference (such as &amp;) in its defline, which NCBI BLAST+ decodes in the outfmt 0 titles; this is not supported by LOSAT's {program}"
        );
    }
    Ok(())
}

/// The title that `GenerateDefline` gives `NStr::HtmlDecode` (`m_MainTitle`, see the
/// module), for the alignment heading (`leave_prefix_suffix` false) or the description
/// table.
fn title_before_decoding(title: &[u8], protein: bool, leave_prefix_suffix: bool) -> Vec<u8> {
    let mut main_title = title.to_vec();
    trim_end_of(&mut main_title, b".,;~ ");
    if main_title.is_empty() && protein {
        main_title = UNNAMED_PROTEIN_PRODUCT.to_vec();
    }
    if !leave_prefix_suffix {
        for prefix in TPA_PREFIXES {
            if main_title.len() >= prefix.len()
                && main_title[..prefix.len()].eq_ignore_ascii_case(prefix.as_bytes())
            {
                main_title.drain(..prefix.len());
                trim_leading_spaces(&mut main_title);
            }
        }
    }
    trim_leading_spaces(&mut main_title);
    main_title
}

/// The entities of `NStr::HtmlDecode`, in NCBI's order.
///
/// NCBI reference: c++/src/corelib/ncbistr.cpp:4221-4232
/// ```c
/// static struct tag_HtmlEntities
/// {
///     TUnicodeSymbol u;
///     const char*    s;
/// }
/// const s_HtmlEntities[] = {
///     {    9, "Tab" },
///     {   10, "NewLine" },
///     {   33, "excl" },
///     {   34, "quot" },
/// ```
/// (and the rows after them, to `{ 9830, "diams" }` at ncbistr.cpp:4508.)
#[rustfmt::skip]
const NCBI_HTML_ENTITIES: [(u32, &str); 280] = [
    (9, "Tab"), (10, "NewLine"), (33, "excl"), (34, "quot"), (35, "num"), (36, "dollar"),
    (37, "percnt"), (38, "amp"), (39, "apos"), (40, "lpar"), (41, "rpar"), (42, "ast"),
    (43, "plus"), (44, "comma"), (46, "period"), (47, "sol"), (58, "colon"), (59, "semi"),
    (60, "lt"), (61, "equals"), (62, "gt"), (63, "quest"), (64, "commat"), (91, "lsqb"),
    (92, "bsol"), (93, "rsqb"), (94, "Hat"), (95, "lowbar"), (96, "grave"), (123, "lcub"),
    (124, "verbar"), (125, "rcub"), (160, "nbsp"), (161, "iexcl"), (162, "cent"),
    (163, "pound"), (164, "curren"), (165, "yen"), (166, "brvbar"), (167, "sect"),
    (168, "uml"), (169, "copy"), (170, "ordf"), (171, "laquo"), (172, "not"), (173, "shy"),
    (174, "reg"), (175, "macr"), (176, "deg"), (177, "plusmn"), (178, "sup2"), (179, "sup3"),
    (180, "acute"), (181, "micro"), (182, "para"), (183, "middot"), (184, "cedil"),
    (185, "sup1"), (186, "ordm"), (187, "raquo"), (188, "frac14"), (189, "frac12"),
    (190, "frac34"), (191, "iquest"), (192, "Agrave"), (193, "Aacute"), (194, "Acirc"),
    (195, "Atilde"), (196, "Auml"), (197, "Aring"), (198, "AElig"), (199, "Ccedil"),
    (200, "Egrave"), (201, "Eacute"), (202, "Ecirc"), (203, "Euml"), (204, "Igrave"),
    (205, "Iacute"), (206, "Icirc"), (207, "Iuml"), (208, "ETH"), (209, "Ntilde"),
    (210, "Ograve"), (211, "Oacute"), (212, "Ocirc"), (213, "Otilde"), (214, "Ouml"),
    (215, "times"), (216, "Oslash"), (217, "Ugrave"), (218, "Uacute"), (219, "Ucirc"),
    (220, "Uuml"), (221, "Yacute"), (222, "THORN"), (223, "szlig"), (224, "agrave"),
    (225, "aacute"), (226, "acirc"), (227, "atilde"), (228, "auml"), (229, "aring"),
    (230, "aelig"), (231, "ccedil"), (232, "egrave"), (233, "eacute"), (234, "ecirc"),
    (235, "euml"), (236, "igrave"), (237, "iacute"), (238, "icirc"), (239, "iuml"),
    (240, "eth"), (241, "ntilde"), (242, "ograve"), (243, "oacute"), (244, "ocirc"),
    (245, "otilde"), (246, "ouml"), (247, "divide"), (248, "oslash"), (249, "ugrave"),
    (250, "uacute"), (251, "ucirc"), (252, "uuml"), (253, "yacute"), (254, "thorn"),
    (255, "yuml"), (338, "OElig"), (339, "oelig"), (352, "Scaron"), (353, "scaron"),
    (376, "Yuml"), (402, "fnof"), (710, "circ"), (732, "tilde"), (913, "Alpha"),
    (914, "Beta"), (915, "Gamma"), (916, "Delta"), (917, "Epsilon"), (918, "Zeta"),
    (919, "Eta"), (920, "Theta"), (921, "Iota"), (922, "Kappa"), (923, "Lambda"), (924, "Mu"),
    (925, "Nu"), (926, "Xi"), (927, "Omicron"), (928, "Pi"), (929, "Rho"), (931, "Sigma"),
    (932, "Tau"), (933, "Upsilon"), (934, "Phi"), (935, "Chi"), (936, "Psi"), (937, "Omega"),
    (945, "alpha"), (946, "beta"), (947, "gamma"), (948, "delta"), (949, "epsilon"),
    (950, "zeta"), (951, "eta"), (952, "theta"), (953, "iota"), (954, "kappa"),
    (955, "lambda"), (956, "mu"), (957, "nu"), (958, "xi"), (959, "omicron"), (960, "pi"),
    (961, "rho"), (962, "sigmaf"), (963, "sigma"), (964, "tau"), (965, "upsilon"),
    (966, "phi"), (967, "chi"), (968, "psi"), (969, "omega"), (977, "thetasym"),
    (978, "upsih"), (982, "piv"), (8194, "ensp"), (8195, "emsp"), (8201, "thinsp"),
    (8204, "zwnj"), (8205, "zwj"), (8206, "lrm"), (8207, "rlm"), (8211, "ndash"),
    (8212, "mdash"), (8216, "lsquo"), (8217, "rsquo"), (8218, "sbquo"), (8220, "ldquo"),
    (8221, "rdquo"), (8222, "bdquo"), (8224, "dagger"), (8225, "Dagger"), (8226, "bull"),
    (8230, "hellip"), (8240, "permil"), (8242, "prime"), (8243, "Prime"), (8249, "lsaquo"),
    (8250, "rsaquo"), (8254, "oline"), (8260, "frasl"), (8364, "euro"), (8472, "weierp"),
    (8465, "image"), (8476, "real"), (8482, "trade"), (8501, "alefsym"), (8592, "larr"),
    (8593, "uarr"), (8594, "rarr"), (8595, "darr"), (8596, "harr"), (8629, "crarr"),
    (8656, "lArr"), (8657, "uArr"), (8658, "rArr"), (8659, "dArr"), (8660, "hArr"),
    (8704, "forall"), (8706, "part"), (8707, "exist"), (8709, "empty"), (8711, "nabla"),
    (8712, "isin"), (8713, "notin"), (8715, "ni"), (8719, "prod"), (8721, "sum"),
    (8722, "minus"), (8727, "lowast"), (8730, "radic"), (8733, "prop"), (8734, "infin"),
    (8736, "ang"), (8743, "and"), (8744, "or"), (8745, "cap"), (8746, "cup"), (8747, "int"),
    (8756, "there4"), (8764, "sim"), (8773, "cong"), (8776, "asymp"), (8800, "ne"),
    (8801, "equiv"), (8804, "le"), (8805, "ge"), (8834, "sub"), (8835, "sup"), (8836, "nsub"),
    (8838, "sube"), (8839, "supe"), (8853, "oplus"), (8855, "otimes"), (8869, "perp"),
    (8901, "sdot"), (8968, "lceil"), (8969, "rceil"), (8970, "lfloor"), (8971, "rfloor"),
    (9001, "lang"), (9002, "rang"), (9674, "loz"), (9824, "spades"), (9827, "clubs"),
    (9829, "hearts"), (9830, "diams"),
];

/// `NStr::HtmlDecode(text)` with its default `eEncoding_Unknown`: the decoded bytes and
/// whether a character reference was decoded (`fHtmlDec_CharRef_Entity` or
/// `fHtmlDec_CharRef_Numeric`). The text is read as `char`s under the C locale
/// (`isalpha` and the others are ASCII only). A numeric reference wraps modulo 2^32 and is
/// written by `x_AppendChar` without any check (`&#0;` gives a NUL byte, `&#xD800;` a
/// surrogate's bytes).
///
/// NCBI reference: c++/src/corelib/ncbistr.cpp:4523-4618
/// ```c
/// string NStr::HtmlDecode(const CTempString str, EEncoding encoding, THtmlDecode* result_flags)
/// {
///     string ustr;
///     THtmlDecode result = 0;
///
///     if (encoding == eEncoding_Unknown) {
///         encoding = CUtf8::GuessEncoding(str);
///         if (encoding == eEncoding_Unknown) {
///             NCBI_THROW2(CStringException, eBadArgs,
///                         "Unable to guess the source string encoding", 0);
///         }
///     }
///     ...
///     for (i = str.begin(); i != e;) {
///         ch = *(i++);
///         //check for HTML entities and character references
///         if (i != e && ch == '&') {
///             CTempString::const_iterator start_of_entity, end_of_entity, itmp;
///             end_of_entity = itmp = i;
///             bool ent, dec, hex, parsed=false;
///             ent = isalpha((unsigned char)(*itmp)) != 0;
///             dec = !ent && *itmp == '#' && ++itmp != e &&
///                   isdigit((unsigned char)(*itmp)) != 0;
///             hex = !dec && itmp != e &&
///                   (*itmp == 'x' || *itmp == 'X') && ++itmp != e &&
///                    isxdigit((unsigned char)(*itmp)) != 0;
///             start_of_entity = itmp;
///
///             if (itmp != e && (ent || dec || hex)) {
///                 // do not look too far
///                 for (int len=0; len<16 && itmp != e; ++len, ++itmp) {
///                     if (*itmp == '&' || *itmp == '#') {
///                         break;
///                     }
///                     if (*itmp == ';') {
///                         end_of_entity = itmp;
///                         break;
///                     }
///                     ent = ent && isalnum( (unsigned char)(*itmp)) != 0;
///                     dec = dec && isdigit( (unsigned char)(*itmp)) != 0;
///                     hex = hex && isxdigit((unsigned char)(*itmp)) != 0;
///                 }
///                 if (end_of_entity != i && (ent || dec || hex)) {
///                     uch = 0;
///                     if (ent) {
///                         string entity(start_of_entity, end_of_entity);
///                         const struct tag_HtmlEntities* p = s_HtmlEntities;
///                         for ( ; p->u != 0; ++p) {
///                             if (entity.compare(p->s) == 0) {
///                                 uch = p->u;
///                                 parsed = true;
///                                 result |= fHtmlDec_CharRef_Entity;
///                                 break;
///                             }
///                         }
///                     } else {
///                         parsed = true;
///                         result |= fHtmlDec_CharRef_Numeric;
///                         for (itmp = start_of_entity; itmp != end_of_entity; ++itmp) {
///                             TUnicodeSymbol ud = *itmp;
///                             if (dec) {
///                                 uch = 10 * uch + (ud - '0');
///                             } else if (hex) {
///                                 if (ud >='0' && ud <= '9') {
///                                     ud -= '0';
///                                 } else if (ud >='a' && ud <= 'f') {
///                                     ud -= 'a';
///                                     ud += 10;
///                                 } else if (ud >='A' && ud <= 'F') {
///                                     ud -= 'A';
///                                     ud += 10;
///                                 }
///                                 uch = 16 * uch + ud;
///                             }
///                         }
///                     }
///                     if (parsed) {
///                         ustr += CUtf8::AsUTF8(&uch,1);
///                         i = ++end_of_entity;
///                         continue;
///                     }
///                 }
///             }
///         }
///         // no entity - append as is
///         if (encoding == eEncoding_UTF8 || encoding == eEncoding_Ascii) {
///             ustr.append( 1, ch );
///         } else {
///             result |= fHtmlDec_Encoding_Changed;
///             ustr += CUtf8::AsUTF8(CTempString(&ch,1), encoding);
///         }
///     }
/// ```
/// NCBI reference: c++/include/corelib/ncbistr.hpp:5729-5759 (`CUtf8::x_Append` of a
/// `TUnicodeSymbol` buffer: `x_TCharToUnicodeSymbol` takes the symbol as it is for a
/// 4-byte character, then `x_AppendChar( u8str, ch )`).
fn html_decode(text: &[u8]) -> Result<(Vec<u8>, bool), UnknownEncoding> {
    let encoding = guess_encoding(text);
    if encoding == Encoding::Unknown {
        return Err(UnknownEncoding);
    }
    let mut decoded = Vec::with_capacity(text.len());
    let mut char_refs = false;
    let e = text.len();
    let mut i = 0;
    while i < e {
        let ch = text[i];
        i += 1;
        if i != e && ch == b'&' {
            let mut itmp = i;
            let mut end_of_entity = i;
            let mut ent = text[itmp].is_ascii_alphabetic();
            let mut dec = false;
            if !ent && text[itmp] == b'#' {
                itmp += 1;
                dec = itmp != e && text[itmp].is_ascii_digit();
            }
            let mut hex = false;
            if !dec && itmp != e && matches!(text[itmp], b'x' | b'X') {
                itmp += 1;
                hex = itmp != e && text[itmp].is_ascii_hexdigit();
            }
            let start_of_entity = itmp;
            if itmp != e && (ent || dec || hex) {
                let mut len = 0;
                while len < 16 && itmp != e {
                    let c = text[itmp];
                    if c == b'&' || c == b'#' {
                        break;
                    }
                    if c == b';' {
                        end_of_entity = itmp;
                        break;
                    }
                    ent = ent && c.is_ascii_alphanumeric();
                    dec = dec && c.is_ascii_digit();
                    hex = hex && c.is_ascii_hexdigit();
                    len += 1;
                    itmp += 1;
                }
                if end_of_entity != i && (ent || dec || hex) {
                    let mut uch: u32 = 0;
                    let mut parsed = false;
                    if ent {
                        let entity = &text[start_of_entity..end_of_entity];
                        if let Some(&(symbol, _)) = NCBI_HTML_ENTITIES
                            .iter()
                            .find(|(_, name)| name.as_bytes() == entity)
                        {
                            uch = symbol;
                            parsed = true;
                        }
                    } else {
                        parsed = true;
                        for &digit in &text[start_of_entity..end_of_entity] {
                            let mut ud = u32::from(digit);
                            if dec {
                                uch = uch.wrapping_mul(10).wrapping_add(ud.wrapping_sub(0x30));
                            } else if hex {
                                if digit.is_ascii_digit() {
                                    ud -= u32::from(b'0');
                                } else if (b'a'..=b'f').contains(&digit) {
                                    ud = ud - u32::from(b'a') + 10;
                                } else if (b'A'..=b'F').contains(&digit) {
                                    ud = ud - u32::from(b'A') + 10;
                                }
                                uch = uch.wrapping_mul(16).wrapping_add(ud);
                            }
                        }
                    }
                    if parsed {
                        char_refs = true;
                        append_char(&mut decoded, uch);
                        i = end_of_entity + 1;
                        continue;
                    }
                }
            }
        }
        match encoding {
            Encoding::Utf8 | Encoding::Ascii => decoded.push(ch),
            _ => append_as_utf8(&mut decoded, ch, encoding),
        }
    }
    Ok((decoded, char_refs))
}

/// NCBI's `EEncoding` values that `CUtf8::GuessEncoding` returns.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Encoding {
    Unknown,
    Utf8,
    Ascii,
    Iso8859_1,
    Windows1252,
    Cesu8,
}

/// `CUtf8::GuessEncoding`. The UTF-8 test checks the lead bytes and the number of
/// continuation bytes only (not the ranges of the second byte after E0, ED, F0 and F4).
///
/// NCBI reference: c++/src/corelib/ncbistr.cpp:6711-6769
/// ```c
/// EEncoding CUtf8::GuessEncoding(const CTempString& src)
/// {
///     SIZE_TYPE more = 0;
///     CTempString::const_iterator i = src.begin();
///     CTempString::const_iterator end = src.end();
///     bool cp1252, iso1, ascii, utf8, cesu8;
///     for (cp1252 = iso1 = ascii = utf8 = true, cesu8=false; i != end; ++i) {
///         Uint1 ch = *i;
///         bool skip = false;
///         if (more != 0) {
///             if (x_EvalNext(ch)) {
///                 --more;
///                 if (more == 0) {
///                     ascii = false;
///                 }
///                 skip = true;
///             } else {
///                 more = 0;
///                 utf8 = false;
///             }
///         }
///         if (ch > 0x7F) {
///             ascii = false;
///             if (ch < 0xA0) {
///                 iso1 = false;
///                 if (ch == 0x81 || ch == 0x8D || ch == 0x8F ||
///                     ch == 0x90 || ch == 0x9D) {
///                     cp1252 = false;
///                 }
///             }
///             if (!skip && utf8 && !x_EvalFirst(ch, more)) {
///                 utf8 = false;
///             }
///             if (utf8 && !cesu8 && ch == 0xED && (end - i) > 5) {
///                 uint8_t c1 = *(i+1);
///                 uint8_t c3 = *(i+3);
///                 uint8_t c4 = *(i+4);
///                 if ( ((c1 & 0xA0) == 0xA0) && (c3 == (uint8_t)0xED) && ((c4 & 0xB0) == 0xB0) ) {
///                     cesu8 = true;
///                 }
///             }
///         }
///     }
///     if (more != 0) {
///         utf8 = false;
///     }
///     if (ascii) {
///         return eEncoding_Ascii;
///     } else if (utf8) {
///         return cesu8 ? eEncoding_CESU8 : eEncoding_UTF8;
///     } else if (cp1252) {
///         return iso1 ? eEncoding_ISO8859_1 : eEncoding_Windows_1252;
///     }
///     return eEncoding_Unknown;
/// }
/// ```
fn guess_encoding(src: &[u8]) -> Encoding {
    let mut more = 0;
    let (mut cp1252, mut iso1, mut ascii, mut utf8, mut cesu8) = (true, true, true, true, false);
    for (index, &ch) in src.iter().enumerate() {
        let mut skip = false;
        if more != 0 {
            if eval_next(ch) {
                more -= 1;
                if more == 0 {
                    ascii = false;
                }
                skip = true;
            } else {
                more = 0;
                utf8 = false;
            }
        }
        if ch > 0x7F {
            ascii = false;
            if ch < 0xA0 {
                iso1 = false;
                if matches!(ch, 0x81 | 0x8D | 0x8F | 0x90 | 0x9D) {
                    cp1252 = false;
                }
            }
            if !skip && utf8 && !eval_first(ch, &mut more) {
                utf8 = false;
            }
            if utf8 && !cesu8 && ch == 0xED && src.len() - index > 5 {
                let c1 = src[index + 1];
                let c3 = src[index + 3];
                let c4 = src[index + 4];
                if (c1 & 0xA0) == 0xA0 && c3 == 0xED && (c4 & 0xB0) == 0xB0 {
                    cesu8 = true;
                }
            }
        }
    }
    if more != 0 {
        utf8 = false;
    }
    if ascii {
        Encoding::Ascii
    } else if utf8 {
        if cesu8 {
            Encoding::Cesu8
        } else {
            Encoding::Utf8
        }
    } else if cp1252 {
        if iso1 {
            Encoding::Iso8859_1
        } else {
            Encoding::Windows1252
        }
    } else {
        Encoding::Unknown
    }
}

/// NCBI reference: c++/src/corelib/ncbistr.cpp:7143-7166
/// ```c
/// bool CUtf8::x_EvalFirst(char ch, SIZE_TYPE& more)
/// {
///     more = 0;
///     if ((ch & 0x80) != 0) {
///         if ((ch & 0xE0) == 0xC0) {
///             if ((ch & 0xFE) == 0xC0) {
///                 // C0 and C1 are not valid UTF-8 chars
///                 return false;
///             }
///             more = 1;
///         } else if ((ch & 0xF0) == 0xE0) {
///             more = 2;
///         } else if ((ch & 0xF8) == 0xF0) {
///             if ((unsigned char)ch > (unsigned char)0xF4) {
///                 // F5-FF are not valid UTF-8 chars
///                 return false;
///             }
///             more = 3;
///         } else {
///             return false;
///         }
///     }
///     return true;
/// }
/// ```
fn eval_first(ch: u8, more: &mut usize) -> bool {
    *more = 0;
    if ch & 0x80 != 0 {
        if ch & 0xE0 == 0xC0 {
            if ch & 0xFE == 0xC0 {
                return false;
            }
            *more = 1;
        } else if ch & 0xF0 == 0xE0 {
            *more = 2;
        } else if ch & 0xF8 == 0xF0 {
            if ch > 0xF4 {
                return false;
            }
            *more = 3;
        } else {
            return false;
        }
    }
    true
}

/// NCBI reference: c++/src/corelib/ncbistr.cpp:7169-7172
/// ```c
/// bool CUtf8::x_EvalNext(char ch)
/// {
///     return (ch & 0xC0) == 0x80;
/// }
/// ```
fn eval_next(ch: u8) -> bool {
    ch & 0xC0 == 0x80
}

/// `CUtf8::AsUTF8(CTempString(&ch,1), encoding)` of one byte of a title that
/// `GuessEncoding` did not find to be ASCII or UTF-8: a CESU-8 title keeps the byte (a
/// one-byte string has no six-byte pair), an ISO-8859-1 or Windows-1252 title converts it
/// (`CharToSymbol`, then `x_AppendChar`).
///
/// NCBI reference: c++/src/corelib/ncbistr.cpp:7062-7115
/// ```c
/// CStringUTF8& CUtf8::x_Append( CStringUTF8& self, const CTempString& src,
///     EEncoding encoding, EValidate validate)
/// {
///     ...
///     if (encoding == eEncoding_CESU8) {
///         ...
///         for (; i != end; ++i) {
///             Uint1 ch = *i;
///             if (ch == 0xED && (end - i) > 5) {
///             ...
///             }
///             self.append(1, ch);
///         }
///         return self;
///     }
///     ...
///     for (i = src.begin(); i != end; ++i) {
///         x_AppendChar( self, CharToSymbol( *i, encoding ) );
///     }
///     return self;
/// }
/// ```
/// NCBI reference: c++/src/corelib/ncbistr.cpp:6851-6881
/// ```c
/// // cp1252, codepoints for chars 0x80 to 0x9F
/// static const TUnicodeSymbol s_cp1252_table[] = {
///     0x20AC, 0x003F, 0x201A, 0x0192, 0x201E, 0x2026, 0x2020, 0x2021,
///     0x02C6, 0x2030, 0x0160, 0x2039, 0x0152, 0x003F, 0x017D, 0x003F,
///     0x003F, 0x2018, 0x2019, 0x201C, 0x201D, 0x2022, 0x2013, 0x2014,
///     0x02DC, 0x2122, 0x0161, 0x203A, 0x0153, 0x003F, 0x017E, 0x0178
/// };
///
/// TUnicodeSymbol CUtf8::CharToSymbol(char c, EEncoding encoding)
/// {
///     Uint1 ch = c;
///     switch (encoding)
///     {
///     ...
///     case eEncoding_Ascii:
///     case eEncoding_ISO8859_1:
///         break;
///     case eEncoding_Windows_1252:
///         if (ch > 0x7F && ch < 0xA0) {
///             return s_cp1252_table[ ch - 0x80 ];
///         }
///         break;
///     ...
///     }
///     return (TUnicodeSymbol)ch;
/// }
/// ```
fn append_as_utf8(out: &mut Vec<u8>, ch: u8, encoding: Encoding) {
    const CP1252_TABLE: [u32; 32] = [
        0x20AC, 0x003F, 0x201A, 0x0192, 0x201E, 0x2026, 0x2020, 0x2021, 0x02C6, 0x2030, 0x0160,
        0x2039, 0x0152, 0x003F, 0x017D, 0x003F, 0x003F, 0x2018, 0x2019, 0x201C, 0x201D, 0x2022,
        0x2013, 0x2014, 0x02DC, 0x2122, 0x0161, 0x203A, 0x0153, 0x003F, 0x017E, 0x0178,
    ];
    match encoding {
        Encoding::Cesu8 => out.push(ch),
        Encoding::Windows1252 if (0x80..0xA0).contains(&ch) => {
            append_char(out, CP1252_TABLE[usize::from(ch - 0x80)]);
        }
        _ => append_char(out, u32::from(ch)),
    }
}

/// NCBI reference: c++/src/corelib/ncbistr.cpp:7040-7060
/// ```c
/// CStringUTF8& CUtf8::x_AppendChar( CStringUTF8& self, TUnicodeSymbol c)
/// {
///     Uint4 ch = c;
///     if (ch < 0x80) {
///         self.append(1, Uint1(ch));
///     }
///     else if (ch < 0x800) {
///         self.append(1, Uint1( (ch >>  6)         | 0xC0));
///         self.append(1, Uint1( (ch        & 0x3F) | 0x80));
///     } else if (ch < 0x10000) {
///         self.append(1, Uint1( (ch >> 12)         | 0xE0));
///         self.append(1, Uint1(((ch >>  6) & 0x3F) | 0x80));
///         self.append(1, Uint1(( ch        & 0x3F) | 0x80));
///     } else {
///         self.append(1, Uint1( (ch >> 18)         | 0xF0));
///         self.append(1, Uint1(((ch >> 12) & 0x3F) | 0x80));
///         self.append(1, Uint1(((ch >>  6) & 0x3F) | 0x80));
///         self.append(1, Uint1( (ch        & 0x3F) | 0x80));
///     }
///     return self;
/// }
/// ```
fn append_char(out: &mut Vec<u8>, ch: u32) {
    if ch < 0x80 {
        out.push(ch as u8);
    } else if ch < 0x800 {
        out.push(((ch >> 6) | 0xC0) as u8);
        out.push(((ch & 0x3F) | 0x80) as u8);
    } else if ch < 0x10000 {
        out.push(((ch >> 12) | 0xE0) as u8);
        out.push((((ch >> 6) & 0x3F) | 0x80) as u8);
        out.push(((ch & 0x3F) | 0x80) as u8);
    } else {
        // `Uint1(...)` keeps the low byte of a symbol above 0x1FFFFF, as `as u8` does.
        out.push(((ch >> 18) | 0xF0) as u8);
        out.push((((ch >> 12) & 0x3F) | 0x80) as u8);
        out.push((((ch >> 6) & 0x3F) | 0x80) as u8);
        out.push(((ch & 0x3F) | 0x80) as u8);
    }
}

/// `find_last_not_of(chars)` then `erase(pos + 1)`; a text made only of `chars` stays.
fn trim_end_of(text: &mut Vec<u8>, chars: &[u8]) {
    if let Some(pos) = text.iter().rposition(|byte| !chars.contains(byte)) {
        text.truncate(pos + 1);
    }
}

/// `erase(0, find_first_not_of(' '))`: a text of spaces only becomes empty.
fn trim_leading_spaces(text: &mut Vec<u8>) {
    let first = text
        .iter()
        .position(|&byte| byte != b' ')
        .unwrap_or(text.len());
    text.drain(..first);
}

/// NCBI's `x_CleanAndCompress`, byte for byte, except that it stops at the end of the
/// string where NCBI's count of the remaining letters wraps (see below); the flag is
/// whether it wrapped. The bytes are C++ `char`s, which are signed: the last byte is
/// written only when it is between 1 and 0x7F (and not a space), so a title that ends with
/// a byte of 0x80 or more loses that byte. In the middle of the title the sign does not
/// matter, because every pattern ends with an ASCII byte.
///
/// NCBI reference: c++/src/objmgr/util/create_defline.cpp:219-312
/// ```c
///     char curr = *in++; // initialize with first character
///     left--;
///
///     char next = 0;
///     Uint2 two_chars = curr; // this is two bytes storage where we see current and previous symbols
///
///     while (left > 0) {
///         next = *in++;
///
///         two_chars = Uint2((two_chars << 8) | next);
///
///         switch (two_chars)
///         {
///         case twocommas: // replace double commas with comma+space
///             *out++ = curr;
///             next = ' ';
///             break;
///         case twospaces: // skip multispaces (only print last one)
///             break;
///         case bracket_space: // skip space after bracket
///             next = curr;
///             two_chars = curr;
///             break;
///         case space_bracket: // skip space before bracket
///             break;
///         case space_comma:
///         case space_semicolon: // swap characters
///             *out++ = next;
///             next = curr;
///             two_chars = curr;
///             break;
///         case comma_space:
///             *out++ = curr;
///             *out++ = ' ';
///             while (next == ' ' || next == ',') {
///                 next = *in;
///                 in++;
///                 left--;
///             }
///             two_chars = next;
///             break;
///         case semicolon_space:
///     ...
///         default:
///             *out++ = curr;
///             break;
///         }
///
///         curr = next;
///         left--;
///     }
///
///     if (curr > 0 && curr != ' ') {
///         *out++ = curr;
///     }
///
///     dest.resize(out - dest.c_str());
///
///     if (isProt) {
///         NStr::ReplaceInPlace (dest, ". [", " [");
///         NStr::ReplaceInPlace (dest, ", [", " [");
///     }
/// ```
fn clean_and_compress(input: &[u8], is_protein: bool) -> (Vec<u8>, bool) {
    let mut start = 0;
    let mut end = input.len();
    while start < end && input[start] == b' ' {
        start += 1;
    }
    while end > start && input[end - 1] == b' ' {
        end -= 1;
    }
    let text = &input[start..end];
    let mut out = Vec::with_capacity(text.len());
    if text.is_empty() {
        return (out, false);
    }
    // Past the end reads the terminating NUL of the C++ string.
    let at = |index: usize| text.get(index).copied().unwrap_or(0);
    let mut index = 0;
    let mut left = text.len();
    let mut wrapped = false;
    let mut curr = at(index);
    index += 1;
    left -= 1;
    let mut two_chars = u16::from(curr);
    while left > 0 {
        let mut next = at(index);
        index += 1;
        two_chars = (two_chars << 8) | u16::from(next);
        match two_chars.to_be_bytes() {
            [b',', b','] => {
                out.push(curr);
                next = b' ';
            }
            [b' ', b' '] | [b' ', b')'] => {}
            [b'(', b' '] => {
                next = curr;
                two_chars = u16::from(curr);
            }
            [b' ', b','] | [b' ', b';'] => {
                out.push(next);
                next = curr;
                two_chars = u16::from(curr);
            }
            [separator @ (b',' | b';'), b' '] => {
                out.push(curr);
                out.push(b' ');
                while next == b' ' || next == separator {
                    next = at(index);
                    index += 1;
                    // Where NCBI's `left` (a size_t) wraps, for a run of spaces and
                    // separators that reaches the end of the string, its loop runs past the
                    // string (it reads the terminating NUL and on, and crashes). LOSAT stops
                    // at the end of the string (approved exception 2 of
                    // PD-LOSAT-NCBI-DEFECTS): `next` is the NUL and nothing more is written.
                    wrapped |= left == 0;
                    left = left.saturating_sub(1);
                }
                two_chars = u16::from(next);
            }
            _ => out.push(curr),
        }
        curr = next;
        wrapped |= left == 0;
        left = left.saturating_sub(1);
    }
    // `curr > 0` of a signed `char`.
    if (1..0x80).contains(&curr) && curr != b' ' {
        out.push(curr);
    }
    if is_protein {
        out = replace_in_place(&out, b". [", b" [");
        out = replace_in_place(&out, b", [", b" [");
    }
    (out, wrapped)
}

/// `NStr::ReplaceInPlace(src, search, replace)` with its defaults (from the start, every
/// occurrence): each search resumes after the text just put in.
///
/// NCBI reference: c++/src/corelib/ncbistr.cpp:3401-3428
/// ```c
///     if ( start_pos + search.size() > src.size()  ||  search == replace )
///         return src;
///
///     bool equal_len = (search.size() == replace.size());
///     for (SIZE_TYPE count = 0; !(max_replace && count >= max_replace); count++){
///         start_pos = src.find(search, start_pos);
///         if (start_pos == NPOS)
///             break;
///         ...
///             src.replace(start_pos, search.size(), replace);
///         ...
///         start_pos += replace.size();
/// ```
fn replace_in_place(src: &[u8], search: &[u8], replace: &[u8]) -> Vec<u8> {
    let mut src = src.to_vec();
    let mut start_pos = 0;
    if start_pos + search.len() > src.len() || search == replace {
        return src;
    }
    while let Some(found) = src[start_pos..]
        .windows(search.len())
        .position(|window| window == search)
    {
        let at = start_pos + found;
        src.splice(at..at + search.len(), replace.iter().copied());
        start_pos = at + replace.len();
        if start_pos > src.len() {
            break;
        }
    }
    src
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blastinput::fasta_reader::{read_all, FastaInputSource, ReaderConfig};

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

    // NCBI reference: c++/src/corelib/ncbistr.cpp:4543-4590 (NStr::HtmlDecode): table
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

    // NCBI reference: c++/src/corelib/ncbistr.cpp:6711-6769 (CUtf8::GuessEncoding).
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

    // NCBI reference: c++/src/corelib/ncbistr.cpp:4588-4606,7040-7060 (numeric references,
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

    // NCBI reference: c++/src/objmgr/util/create_defline.cpp:302-311 (signed `curr`, the
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
}
