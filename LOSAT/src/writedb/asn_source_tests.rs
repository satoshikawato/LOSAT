use super::*;
use crate::report::serial::asn_input::Utf8TagOnce;

fn read_all_bioseqs(
    input: &[u8],
    format: AsnFormat,
    skip_unver: bool,
) -> Result<Vec<AsnBioseq>, WriteDbError> {
    let mut cursor = std::io::Cursor::new(input);
    let mut stream = AsnStream::new(&mut cursor, format);
    let mut utf8_once = Utf8TagOnce::default();
    let mut warn = |_: &[u8]| {};
    let mut out = Vec::new();
    while let Some(mut src) = stream.next_source(skip_unver, &mut warn, &mut utf8_once)? {
        while let Some(bio) = src.next_bioseq(&mut warn)? {
            out.push(bio);
        }
    }
    Ok(out)
}

fn read_text_ok(input: &[u8]) -> Vec<AsnBioseq> {
    read_all_bioseqs(input, AsnFormat::Text, false).expect("read")
}

fn read_text_err(input: &[u8]) -> WriteDbError {
    read_all_bioseqs(input, AsnFormat::Text, false).unwrap_err()
}

fn unwrap_ncbi_exception(err: WriteDbError) -> NcbiException {
    let WriteDbError::Exception(e) = err else {
        panic!("expected NcbiException, got {err:?}");
    };
    e
}

fn read_text_ncbi_what(input: &[u8]) -> String {
    let e = unwrap_ncbi_exception(read_text_err(input));
    String::from_utf8_lossy(&e.ncbi_what()).into_owned()
}

fn expect_delta_err_what(input: &[u8]) -> String {
    let mut bioseqs = read_text_ok(input);
    let e = unwrap_ncbi_exception(fix_asn_delta(&mut bioseqs[0]).unwrap_err());
    String::from_utf8_lossy(&e.what()).into_owned()
}

#[test]
fn seq_entry_and_bioseq_set_text_yield_bioseqs_in_order() {
    let input = b"\
Seq-entry ::= seq {
  id { local str \"s1\" },
  descr { title \"first\" },
  inst { repr raw, mol dna, length 4, seq-data iupacna \"ACGT\" }
}
Bioseq-set ::= {
  seq-set {
    seq {
      id { local str \"s2\" },
      descr { title \"second\" },
      inst { repr raw, mol dna, length 4, seq-data iupacna \"TGCA\" }
    }
  }
}
";
    let bioseqs = read_text_ok(input);
    assert_eq!(bioseqs.len(), 2);
    assert_eq!(bioseqs[0].first_id_fasta(), "lcl|s1");
    assert_eq!(bioseqs[0].title, b"first");
    assert_eq!(
        bioseqs[0].seq_data,
        Some(SeqData::Iupacna(b"ACGT".to_vec()))
    );
    assert_eq!(bioseqs[1].first_id_fasta(), "lcl|s2");
    assert_eq!(bioseqs[1].title, b"second");
    assert_eq!(
        bioseqs[1].seq_data,
        Some(SeqData::Iupacna(b"TGCA".to_vec()))
    );
}

#[test]
fn unverified_nuc_prot_and_bioseq_are_skipped_when_flag_enabled() {
    let input = b"\
Seq-entry ::= set {
  class nuc-prot,
  descr { user { type str \"Unverified\", data { } } },
  seq-set {
    seq {
      id { local str \"np1\" },
      descr { title \"u1\" },
      inst { repr raw, mol dna, length 4, seq-data iupacna \"CGTA\" }
    }
  }
}
Seq-entry ::= seq {
  id { local str \"ok1\" },
  descr { title \"good\" },
  inst { repr raw, mol dna, length 4, seq-data iupacna \"AAAA\" }
}
";
    let skipped = read_all_bioseqs(input, AsnFormat::Text, true).expect("skip");
    assert_eq!(skipped.len(), 1);
    assert_eq!(skipped[0].first_id_fasta(), "lcl|ok1");

    // Without skip_unver, Unverified title prefix is prepended by CDeflineGenerator (D3c-5).
    let kept = read_text_ok(input);
    assert_eq!(kept.len(), 2);
    assert_eq!(kept[0].title, b"UNVERIFIED: u1");
    assert_eq!(kept[1].title, b"good");
}

#[test]
fn autodef_is_refused_only_when_entered_and_bypassed_for_tls_and_refseq_contig() {
    let autodef_nucl = b"\
Seq-entry ::= seq {
  id { other { accession \"NZ_000001\" } },
  descr {
    source { org { taxname \"Escherichia coli\" } },
    user { type str \"AutodefOptions\", data { } }
  },
  inst { repr raw, mol dna, length 4, seq-data iupacna \"ACGT\" }
}
";
    let err = read_text_err(autodef_nucl);
    assert!(matches!(err, WriteDbError::Unsupported(ref msg) if msg.contains("AutoDef")));

    // NT_ (refseq contig -> m_IsNZ) and tech targeted (m_IsTLS) bypass AutoDef in NCBI.
    let bypassed = b"\
Seq-entry ::= seq {
  id { other { accession \"NT_000001\" } },
  descr {
    source { org { taxname \"Homo sapiens\" } },
    user { type str \"AutodefOptions\", data { } }
  },
  inst { repr raw, mol dna, length 4, seq-data iupacna \"ACGT\" }
}
Seq-entry ::= seq {
  id { genbank { accession \"KAAA01000001\" } },
  descr {
    source { org { taxname \"Escherichia coli\" } },
    molinfo { tech targeted },
    user {
      type str \"AutodefOptions\",
      data { { label str \"Targeted Locus Name\", data str \"16S ribosomal RNA\" } }
    }
  },
  inst { repr raw, mol dna, length 4, seq-data iupacna \"ACGT\" }
}
";
    let bioseqs = read_text_ok(bypassed);
    assert_eq!(bioseqs.len(), 2);
    assert_eq!(bioseqs[0].title, b"Homo sapiens");
    assert_eq!(
        bioseqs[1].title,
        b"TLS: Escherichia coli 16S ribosomal RNA, sequence"
    );
}

#[test]
fn duplicate_seq_id_in_same_object_throws_objmgr_exception() {
    let input = b"\
Seq-entry ::= seq {
  id { local str \"Ab\", local str \"ab\" },
  descr { title \"dup\" },
  inst { repr raw, mol dna, length 4, seq-data iupacna \"GGCC\" }
}
";
    let dup_what = read_text_ncbi_what(input);
    assert!(
        dup_what.contains("CObjMgrException::eAddDataError"),
        "{dup_what}"
    );
    assert!(
        dup_what.contains("duplicate Bioseq id lcl|ab present in"),
        "{dup_what}"
    );
    assert!(dup_what.contains("seq1: lcl|Ab | lcl|ab"), "{dup_what}");
    assert!(dup_what.contains("seq2: lcl|Ab | lcl|ab"), "{dup_what}");
}

#[test]
fn missing_inst_length_throws_unassigned_member() {
    let input = b"\
Seq-entry ::= seq {
  id { local str \"nolen\" },
  inst { repr raw, mol dna, seq-data iupacna \"ACGT\" }
}
";
    let nolen_what = read_text_ncbi_what(input);
    assert!(
        nolen_what.contains("CUnassignedMember::eGet"),
        "{nolen_what}"
    );
    assert!(
        nolen_what.contains("CSeq_inst::GetLength()"),
        "{nolen_what}"
    );
}

#[test]
fn parent_set_taxid_and_source_taxname_propagate_to_bioseq() {
    let input = b"\
Seq-entry ::= set {
  class nuc-prot,
  descr {
    source { org { taxname \"Homo sapiens\", db { { db \"taxon\", tag id 9606 } } } }
  },
  seq-set {
    seq {
      id { local str \"p1\" },
      inst { repr raw, mol aa, length 4, seq-data ncbieaa \"MKTW\" }
    }
  }
}
";
    let bioseqs = read_text_ok(input);
    assert_eq!(bioseqs.len(), 1);
    assert_eq!(bioseqs[0].title, b"unnamed protein product [Homo sapiens]");
    assert_eq!(bioseqs[0].descr_taxids, vec![9606]);
}

fn read_fixed_delta(input: &[u8]) -> AsnBioseq {
    let mut bioseqs = read_text_ok(input);
    assert_eq!(bioseqs.len(), 1);
    fix_asn_delta(&mut bioseqs[0]).expect("fix delta");
    bioseqs.remove(0)
}

#[test]
fn fix_asn_delta_packs_ncbi2na_ncbi4na_gap_and_loc_null() {
    let input = b"\
Seq-entry ::= seq {
  id { local str \"d1\" },
  descr { title \"delta\" },
  inst {
    repr delta,
    mol dna,
    length 9,
    ext delta {
      literal { length 4, seq-data ncbi2na '1B'H },
      loc null NULL,
      literal { length 2, seq-data gap { type unknown } },
      literal { length 2, seq-data ncbi4na '12'H }
    }
  }
}
";
    let bio = read_fixed_delta(input);
    assert_eq!(bio.length, 9);
    // 0x1B = 00 01 10 11 = A C G T -> 4na: 0x12, 0x48; then 3 N's (0xFF, 0xF?), then 0x12 -> 0xF1, 0x20
    assert_eq!(
        bio.seq_data,
        Some(SeqData::Ncbi4na(vec![0x12, 0x48, 0xff, 0xf1, 0x20]))
    );
}

#[test]
fn fix_asn_delta_does_not_over_reserve_huge_declared_length() {
    let input = b"\
Seq-entry ::= seq {
  id { local str \"d_huge\" },
  descr { title \"huge declared length\" },
  inst {
    repr delta,
    mol dna,
    length 2147483647,
    ext delta {
      literal { length 4, seq-data ncbi2na '1B'H }
    }
  }
}
";
    let bio = read_fixed_delta(input);
    assert_eq!(bio.length, 4);
    assert_eq!(bio.seq_data, Some(SeqData::Ncbi4na(vec![0x12, 0x48])));
}

#[test]
fn fix_asn_delta_rejects_protein_and_iupacna_with_exact_ncbi_errors() {
    let prot_delta = b"\
Seq-entry ::= seq {
  id { local str \"dp\" },
  descr { title \"p\" },
  inst { repr delta, mol aa, length 2, ext delta { literal { length 2 } } }
}
";
    let what = expect_delta_err_what(prot_delta);
    assert!(
        what.contains("lcl|dp: Protein delta sequences are not supported."),
        "{what}"
    );

    let iupac_delta = b"\
Seq-entry ::= seq {
  id { local str \"di\" },
  descr { title \"i\" },
  inst {
    repr delta,
    mol dna,
    length 4,
    ext delta { literal { length 4, seq-data iupacna \"ACGT\" } }
  }
}
";
    let what = expect_delta_err_what(iupac_delta);
    assert!(
        what.contains("Unhandled type of sequence data encountered."),
        "{what}"
    );
}
