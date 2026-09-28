//! Comparison harness; outputs prepared state, never a search report.
use anyhow::{bail, Result};
use LOSAT::algorithm::blastx::{
    args::ResolvedOptions,
    input::{read_fasta, FastaRecord},
    query_setup::{batch_ranges, prepare_queries},
};
use LOSAT::cli::{Cli, Commands};
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:142-147
// ```c++
//         ddc.Log(prefix+string("query_offset"), m_Ptr->contexts[i].query_offset);
//         ddc.Log(prefix+string("query_length"), m_Ptr->contexts[i].query_length);
//         ddc.Log(prefix+string("eff_searchsp"), m_Ptr->contexts[i].eff_searchsp);
//         ddc.Log(prefix+string("length_adjustment"),
//                 m_Ptr->contexts[i].length_adjustment);
//         ddc.Log(prefix+string("query_index"), m_Ptr->contexts[i].query_index);
// ```
fn hex(data: &[u8]) -> String {
    data.iter().map(|b| format!("{b:02x}")).collect()
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:879-881
// ```c++
//                             : & mask_loc->seqloc_array[ctx_idx+context]),
//                             from, to);
//             }
// ```
fn masks(data: &[(i32, i32)]) -> String {
    if data.is_empty() {
        "-".into()
    } else {
        data.iter()
            .map(|(a, b)| format!("{a}:{b}"))
            .collect::<Vec<_>>()
            .join(",")
    }
}
// NCBI reference (598d8ae6): c++/src/objmgr/seq_vector.cpp:1281-1290
// ```c++
// void CSeqVector::SetIupacCoding(void)
// {
//     SetCoding(IsProtein()? CSeq_data::e_Iupacaa: CSeq_data::e_Iupacna);
// }
//
//
// void CSeqVector::SetNcbiCoding(void)
// {
//     SetCoding(IsProtein()? CSeq_data::e_Ncbistdaa: CSeq_data::e_Ncbi4na);
// }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:177-180
// ```c++
//     size_t title_start = NPOS;
//     if ((fFastaFlags & CFastaReader::fNoParseID)) {
//         title_start = start;
//     }
// ```
fn records(data: &[FastaRecord], kind: &str, base: usize) {
    for (i, r) in data.iter().enumerate() {
        println!(
            "R\t{kind}\t{}\t{}\t{}\t{}\t{}",
            base + i,
            r.internal_id,
            hex(r.title.as_bytes()),
            hex(&if kind == "S" {
                r.sequence
                    .iter()
                    .map(|&aa| LOSAT::utils::matrix::aa_char_to_ncbistdaa(aa))
                    .collect::<Vec<_>>()
            } else {
                r.sequence.clone()
            }),
            masks(&r.lowercase_masks)
        );
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blastx_options.cpp:57-70
// ```c++
//     CBlastProteinOptionsHandle::SetLookupTableDefaults();
//     m_Opts->SetWordThreshold(BLAST_WORD_THRESHOLD_BLASTX);
// }
//
// void
// CBlastxOptionsHandle::SetQueryOptionDefaults()
// {
//     CBlastProteinOptionsHandle::SetQueryOptionDefaults();
//     m_Opts->SetStrandOption(objects::eNa_strand_both);
//     m_Opts->SetQueryGeneticCode(BLAST_GENETIC_CODE);
//     SetSegFiltering(false); // disable SEG filtering because of eCompositionMatrixAdjust mode
// }
//
// void
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1973-1987
// ```c++
//         m_Strand = eNa_strand_unknown;
//
//         if (!Blast_QueryIsProtein(opt.GetProgramType())) {
//
//             if (args.Exist(kArgStrand) && args[kArgStrand]) {
//                 const string& kStrand = args[kArgStrand].AsString();
//                 if (kStrand == "both") {
//                     m_Strand = eNa_strand_both;
//                 } else if (kStrand == "plus") {
//                     m_Strand = eNa_strand_plus;
//                 } else if (kStrand == "minus") {
//                     m_Strand = eNa_strand_minus;
//                 } else {
//                     abort();
//                 }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2940-2947
// ```c++
//     	    m_LineLength = args[kArgLineLength].AsInteger();
//     	}
//         if(args.Exist(kArgSortHits) && args[kArgSortHits])
//         {
//        	    m_HitsSortOption = args[kArgSortHits].AsInteger();
//         }
//     }
//     else
// ```
fn options(o: &ResolvedOptions) {
    println!("CLI_EXTRACTED");
    macro_rules! field {
        ($n:literal,$v:expr) => {
            println!("{}\t{}", $n, $v)
        };
    }
    field!("WordSize", o.word_size);
    field!("WordThreshold", format!("{:.17}", o.threshold));
    field!("LookupTableType", o.lookup_type);
    field!("MatrixName", o.matrix);
    field!("GapOpeningCost", o.gap_open);
    field!("GapExtensionCost", o.gap_extend);
    field!("GappedMode", o.gapped as u8);
    field!("CompositionBasedStats", o.composition);
    field!("SegFiltering", o.seg.enabled as u8);
    if o.seg.enabled {
        field!("SegFilteringWindow", o.seg.window);
        field!("SegFilteringLocut", format!("{:.17}", o.seg.locut));
        field!("SegFilteringHicut", format!("{:.17}", o.seg.hicut));
    }
    field!("MaskAtHash", o.soft_masking as u8);
    field!("StrandOption", 3);
    field!("QueryGeneticCode", o.query_gencode);
    field!("WindowSize", o.window_size);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:99-115
    // ```c++
    //
    // void
    // CBlastProteinOptionsHandle::SetInitialWordOptionsDefaults()
    // {
    //     SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_PROT);
    //     SetWindowSize(BLAST_WINDOW_SIZE_PROT);
    // }
    //
    // void
    // CBlastProteinOptionsHandle::SetGappedExtensionDefaults()
    // {
    //     SetGapXDropoff(BLAST_GAP_X_DROPOFF_PROT);
    //     SetGapXDropoffFinal(BLAST_GAP_X_DROPOFF_FINAL_PROT);
    //     SetGapTrigger(BLAST_GAP_TRIGGER_PROT);
    //     m_Opts->SetGapExtnAlgorithm(eDynProgScoreOnly);
    //     m_Opts->SetGapTracebackAlgorithm(eDynProgTbck);
    // }
    // ```
    field!("XDropoff", o.x_dropoff);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:99-115
    // ```c++
    //
    // void
    // CBlastProteinOptionsHandle::SetInitialWordOptionsDefaults()
    // {
    //     SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_PROT);
    //     SetWindowSize(BLAST_WINDOW_SIZE_PROT);
    // }
    //
    // void
    // CBlastProteinOptionsHandle::SetGappedExtensionDefaults()
    // {
    //     SetGapXDropoff(BLAST_GAP_X_DROPOFF_PROT);
    //     SetGapXDropoffFinal(BLAST_GAP_X_DROPOFF_FINAL_PROT);
    //     SetGapTrigger(BLAST_GAP_TRIGGER_PROT);
    //     m_Opts->SetGapExtnAlgorithm(eDynProgScoreOnly);
    //     m_Opts->SetGapTracebackAlgorithm(eDynProgTbck);
    // }
    // ```
    field!("GapXDropoff", o.gap_x_dropoff);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:99-115
    // ```c++
    //
    // void
    // CBlastProteinOptionsHandle::SetInitialWordOptionsDefaults()
    // {
    //     SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_PROT);
    //     SetWindowSize(BLAST_WINDOW_SIZE_PROT);
    // }
    //
    // void
    // CBlastProteinOptionsHandle::SetGappedExtensionDefaults()
    // {
    //     SetGapXDropoff(BLAST_GAP_X_DROPOFF_PROT);
    //     SetGapXDropoffFinal(BLAST_GAP_X_DROPOFF_FINAL_PROT);
    //     SetGapTrigger(BLAST_GAP_TRIGGER_PROT);
    //     m_Opts->SetGapExtnAlgorithm(eDynProgScoreOnly);
    //     m_Opts->SetGapTracebackAlgorithm(eDynProgTbck);
    // }
    // ```
    field!("GapXDropoffFinal", o.gap_x_dropoff_final);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:99-115
    // ```c++
    //
    // void
    // CBlastProteinOptionsHandle::SetInitialWordOptionsDefaults()
    // {
    //     SetXDropoff(BLAST_UNGAPPED_X_DROPOFF_PROT);
    //     SetWindowSize(BLAST_WINDOW_SIZE_PROT);
    // }
    //
    // void
    // CBlastProteinOptionsHandle::SetGappedExtensionDefaults()
    // {
    //     SetGapXDropoff(BLAST_GAP_X_DROPOFF_PROT);
    //     SetGapXDropoffFinal(BLAST_GAP_X_DROPOFF_FINAL_PROT);
    //     SetGapTrigger(BLAST_GAP_TRIGGER_PROT);
    //     m_Opts->SetGapExtnAlgorithm(eDynProgScoreOnly);
    //     m_Opts->SetGapTracebackAlgorithm(eDynProgTbck);
    // }
    // ```
    field!("GapTrigger", o.gap_trigger);
    field!("HitlistSize", o.hitlist_size);
    field!("MaxHspsPerSubject", o.max_hsps);
    field!("EvalueThreshold", format!("{:.17}", o.evalue));
    field!("SumStatisticsMode", o.sum_stats as u8);
    field!("LongestIntronLength", o.max_intron_length);
    field!("CullingLimit", o.culling_limit);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:150-156
    // ```c++
    // void
    // CBlastProteinOptionsHandle::SetEffectiveLengthsOptionsDefaults()
    // {
    //     SetDbLength(0);
    //     SetDbSeqNum(0);
    //     SetEffectiveSearchSpace(0);
    // }
    // ```
    field!("DbLength", o.db_length);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_prot_options.cpp:150-156
    // ```c++
    // void
    // CBlastProteinOptionsHandle::SetEffectiveLengthsOptionsDefaults()
    // {
    //     SetDbLength(0);
    //     SetDbSeqNum(0);
    //     SetEffectiveSearchSpace(0);
    // }
    // ```
    field!("EffectiveSearchSpace", o.effective_search_space);
    println!("FORMATTER");
    field!("NumDescriptions", o.num_descriptions);
    field!("NumAlignments", o.num_alignments);
    field!("LineLength", 60);
    field!("OutputFormat", o.outfmt);
    field!("BatchSize", 10002);
}
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:208-220
// ```c++
//                                                    db_adapter);
//         CBlastInputSourceConfig iconfig(dlconfig, query_opts->GetStrand(),
//                                      query_opts->UseLowercaseMasks(),
//                                      query_opts->GetParseDeflines(),
//                                      query_opts->GetRange());
//         if(IsIStreamEmpty(m_CmdLineArgs->GetInputStream())){
//            	ERR_POST(Warning << "Query is Empty!");
//            	return BLAST_EXIT_SUCCESS;
//         }
//         CBlastFastaInputSource fasta(m_CmdLineArgs->GetInputStream(), iconfig);
//         CBlastInput input(&fasta, m_CmdLineArgs->GetQueryBatchSize());
//
//         /*** Get the formatting options ***/
// ```
// NCBI reference (598d8ae6): c++/src/app/blast/blastx_app.cpp:260-264
// ```c++
// 	    BLAST_PROF_START( APP.LOOP.PRE );
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//             CRef<IQueryFactory> queries(new CObjMgr_QueryFactory(*query_batch));
//
//             SaveSearchStrategy(args, m_CmdLineArgs, queries, m_OptsHndl);
// ```
fn main() -> Result<()> {
    let mut a = std::env::args_os();
    let exe = a.next().unwrap();
    let mode = a.next().unwrap();
    let cli: Cli = LOSAT::cli::try_parse_from(std::iter::once(exe).chain(a))?;
    let Commands::Blastx(args) = cli.command else {
        bail!("blastx required")
    };
    let o = args.resolve()?;
    if mode == "options" {
        options(&o);
        return Ok(());
    }
    if mode == "fields" {
        println!("{}", o.fields.join(" "));
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:998-1014
    // ```c++
    //     if( ! bad_pos_vec.empty() ) {
    //         if (TestFlag(fValidate)) {
    //                         NCBI_THROW2(CBadResiduesException, eBadResidues,
    //                 "CFastaReader: There are invalid " + x_NucOrProt() + "residue(s) in input sequence",
    //                 CBadResiduesException::SBadResiduePositions( m_BestID, bad_pos_vec, bad_pos_line_num ) );
    //         } else {
    //             stringstream warn_strm;
    //             warn_strm << "FASTA-Reader: Ignoring invalid " << x_NucOrProt()
    //                 << "residues at position(s): ";
    //             CBadResiduesException::SBadResiduePositions(
    //                 m_BestID, bad_pos_vec, bad_pos_line_num ).ConvertBadIndexesToString(warn_strm);
    //
    //             FASTA_WARNING(0,
    //                 warn_strm.str(),
    //                 ILineError::eProblem_InvalidResidue,
    //                 kEmptyStr );
    //         }
    // ```
    if mode == "warnings" {
        for (kind, path, protein) in [("S", &args.subject, true), ("Q", &args.query, false)] {
            for (i, record) in read_fasta(path, protein, args.lcase_masking)?
                .iter()
                .enumerate()
            {
                for w in &record.warnings {
                    println!(
                        "W\t{kind}\t{i}\t{}\t{}\t{}",
                        w.line,
                        w.kind,
                        w.positions
                            .iter()
                            .map(usize::to_string)
                            .collect::<Vec<_>>()
                            .join(",")
                    );
                }
            }
        }
        return Ok(());
    }
    let s = read_fasta(&args.subject, true, args.lcase_masking)?;
    records(&s, "S", 0);
    let q = read_fasta(&args.query, false, args.lcase_masking)?;
    for (batch, r) in batch_ranges(&q).into_iter().enumerate() {
        records(&q[r.clone()], "Q", r.start);
        println!("A\t{batch}\t{}\t{}", r.start, r.len());
        let p = prepare_queries(&q[r], &o)?;
        println!(
            "I\t{}\t{}\t{}\t{}\t{}",
            p.first_context,
            p.max_length,
            p.min_length,
            p.sequence_start.len(),
            p.split_eligible as u8
        );
        println!("P\t{}", hex(&p.sequence_start_nomask));
        for (i, c) in p.contexts.iter().enumerate() {
            println!(
                "C\t{i}\t{}\t{}\t{}\t{}\t{}\t{}",
                c.query_index,
                c.frame,
                c.offset,
                c.length,
                c.is_valid as u8,
                masks(&c.lowercase_masks)
            );
        }
        for (i, c) in p.contexts.iter().enumerate() {
            println!("F\t{i}\t{}", masks(&c.masks));
        }
        println!(
            "B\t{}\nN\t{}\nL\t{}",
            hex(&p.sequence_start),
            hex(&p.sequence_start_nomask),
            masks(&p.lookup_segments)
        );
        for (i, c) in p.contexts.iter().enumerate() {
            println!("D\t{i}\t{}", masks(&c.dna_masks));
        }
    }
    Ok(())
}
