// Comparison-only pinned NCBI C++ API local-subject oracle for Stage E.
// Adapted from the Stage A database oracle. Never linked into LOSAT.
// NCBI source commit: 598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4.
// Source-derived search and formatting oracle for code 32.
// Code-1 outfmt 0/6/7 must calibrate byte-for-byte against the pinned CLI.
// Comparison-only NCBI C++ API oracle. This file is not part of LOSAT builds.
#include <corelib/ncbistd.hpp>
#include <algo/blast/api/blast_aux.hpp>
#include <algo/blast/core/blast_seqsrc_impl.h>
#include <algo/blast/core/blast_def.h>
#include <algo/blast/core/gencode_singleton.h>
#include <algo/blast/api/tblastn_options.hpp>
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/local_db_adapter.hpp>
#include <algo/blast/api/uniform_search.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/blastinput/blast_args.hpp>
#include <algo/blast/blastinput/blast_input_aux.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>
#include <algo/blast/blastinput/blast_input.hpp>
#include <algo/blast/format/blast_format.hpp>
#include <objects/seq/Bioseq.hpp>
#include <objects/seq/Seq_inst.hpp>
#include <objects/seq/Seq_data.hpp>
#include <objects/seq/Seqdesc.hpp>
#include <objects/seq/Seq_descr.hpp>
#include <objects/seqset/Seq_entry.hpp>
#include <objects/seqloc/Seq_loc.hpp>
#include <objmgr/object_manager.hpp>
#include <objmgr/scope.hpp>
#include <fstream>
#include <iostream>
#include <stdexcept>

USING_NCBI_SCOPE;
USING_SCOPE(objects);
USING_SCOPE(blast);

// NCBI blast_args.cpp:2538-2557 calls ReadSequencesToBlast for local
// subjects. tblastn_app.cpp:200-229 uses CBlastFastaInputSource and
// CBlastInput for the protein query. The API oracle calls those same
// readers so identifiers, titles, and formatter flags match the CLI.
// NCBI c++/src/algo/blast/core/blast_engine.c:1460-1466,762-774:
// the sequence source may provide seq->gen_code_string; the engine uses it
// for translated subject search. The stock local MultiSeq adapter reconstructs
// SSeqLoc without genetic_code_id (seqsrc_multiseq.cpp:140-164), so this
// comparison-only adapter supplies FindGeneticCode(code) on retrieval.
static GetSeqBlkFnPtr original_get_sequence = nullptr;
static Uint1* selected_genetic_code = nullptr;
static Int2 GetSelectedCodeSequence(void* data, BlastSeqSrcGetSeqArg* arg) {
    const Int2 status = original_get_sequence(data, arg);
    if (status == BLAST_SEQSRC_SUCCESS && arg->seq)
        arg->seq->gen_code_string = selected_genetic_code;
    return status;
}

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "usage: oracle QUERY.faa SUBJECT.fna CODE OUTFMT\n";
        return 2;
    }
    const int code = std::stoi(argv[3]);
    const int fmt = std::stoi(argv[4]);
    if (fmt != 0 && fmt != 6 && fmt != 7) return 2;

    // NCBI c++/src/algo/blast/api/blast_aux.cpp:588-613:
    //   const string kGenCode = CGen_code_table::GetNcbieaa(genetic_code);
    //   CSeqportUtil::Convert(gc_ncbieaa, &gc_ncbistdaa, CSeq_data::e_Ncbistdaa);
    auto gc = FindGeneticCode(code);
    if (!gc.get()) throw std::runtime_error("NCBI rejected genetic code");
    // NCBI c++/src/algo/blast/api/blast_aux.cpp:629-646:
    //   GenCodeSingletonInit();
    //   if (GenCodeSingletonFind(genetic_code) == NULL) {
    //       TAutoUint1ArrayPtr gc = FindGeneticCode(genetic_code);
    //       GenCodeSingletonAdd(genetic_code, gc.get());
    //   }
    CAutomaticGenCodeSingleton genetic_codes(code);

    std::ifstream subject_in(argv[2]);
    if (!subject_in) throw std::runtime_error("subject FASTA missing");
    CRef<CBlastQueryVector> subjects;
    CRef<CScope> scope = ReadSequencesToBlast(subject_in, false, TSeqRange(),
                                               false, false, subjects);
    std::ifstream query_in(argv[1]);
    if (!query_in) throw std::runtime_error("query FASTA missing");
    SDataLoaderConfig query_dlconfig(true);
    CBlastInputSourceConfig query_config(query_dlconfig);
    CRef<CBlastFastaInputSource> query_fasta(
        new CBlastFastaInputSource(query_in, query_config));
    CRef<CBlastInput> query_input(new CBlastInput(query_fasta));
    CRef<CBlastQueryVector> queries = query_input->GetAllSeqs(*scope);
    CRef<IQueryFactory> qf(new CObjMgr_QueryFactory(*queries));

    // NCBI c++/include/algo/blast/api/local_db_adapter.hpp:73-76:
    // CLocalDbAdapter(subject_sequences, opts_handle, dbscan_mode=false)
    // supplies the local -subject-like BlastSeqSrc, including its length.
    // NCBI c++/src/algo/blast/api/seqsrc_multiseq.cpp:140-164:
    // the stock local adapter reconstructs SSeqLoc; the comparison-only
    // BlastSeqSrc callback below supplies the selected genetic code.
    CRef<IQueryFactory> sf(new CObjMgr_QueryFactory(*subjects));
    CRef<CTBlastnOptionsHandle> options(new CTBlastnOptionsHandle);
    // NCBI c++/src/algo/blast/blastinput/blast_args.cpp:327-343,387-416:
    //   the CLI defaults SEG to "12 2.2 2.5" and applies those parameters.
    options->SetSegFiltering(true);
    options->SetSegFilteringWindow(12);
    options->SetSegFilteringLocut(2.2);
    options->SetSegFilteringHicut(2.5);
    options->SetMaskAtHash(false);
    // NCBI c++/src/algo/blast/blastinput/blast_args.cpp:1048-1056:
    //   opt.SetDbGeneticCode(args[kArgDbGeneticCode].AsInteger());
    options->SetDbGeneticCode(code);
    // NCBI c++/src/app/blast/blast_app_util.cpp:203-210;
    // c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-180,290-297;
    // c++/src/algo/blast/core/blast_engine.c:1407,1434-1443:
    // ```c++
    // db_adapter.Reset(new CLocalDbAdapter(subjects, opts_hndl, true));
    // if (dbscan_mode) m_iTotalLength += (Int8) (*iter)->length;
    // ```
    // The CLI's local -subject path uses dbscan_mode=true. Its positive
    // TotLen skips BLAST_OneSubjectUpdateParameters for this fixture.
    CRef<CLocalDbAdapter> db(new CLocalDbAdapter(sf, options, true));
    BlastSeqSrc* local_source = db->MakeSeqSrc();
    selected_genetic_code = GenCodeSingletonFind(code);
    if (!selected_genetic_code) throw std::runtime_error("missing registered code");
    original_get_sequence = _BlastSeqSrcImpl_GetGetSequence(local_source);
    _BlastSeqSrcImpl_SetGetSequence(local_source, GetSelectedCodeSequence);

    // NCBI c++/src/app/blast/tblastn_app.cpp:242-265,275-305,342:
    //   CBlastFormat formatter(opt, *db_adapter, ... opt.GetDbGeneticCode() ...);
    //   CLocalBlast lcl_blast(query_factory, m_OptsHndl, db_adapter);
    //   results = lcl_blast.Run(); formatter.PrintOneResultSet(**result, query);
    CLocalBlast search(qf, options, db);
    CRef<CSearchResultSet> results = search.Run();
    CFormattingArgs::EOutputFormat output = fmt == 0 ? CFormattingArgs::ePairwise :
        fmt == 6 ? CFormattingArgs::eTabular : CFormattingArgs::eTabularWithComments;
    // NCBI c++/src/app/blast/tblastn_app.cpp:242-265:
    // CBlastFormat(..., GetCmdlineArgs(GetArguments()), GetSubjectFile(args));
    // NCBI c++/src/algo/blast/format/blast_format.cpp:68-103,790-803:
    // m_SubjectTag(subjectTag);
    // dbname = "User specified sequence set (Input: " + m_SubjectTag + ")";
    CBlastFormat formatter(options->GetOptions(), *db, output, false, std::cout,
                           500, 250, *scope, "BLOSUM62", false, false,
                           1, code, options->GetOptions().GetSumStatisticsMode(),
                           false, -1, kEmptyStr, false, false, NULL, NULL,
                           kEmptyStr, argv[2]);
    formatter.SetHitsSortOption(-1);
    formatter.SetHspsSortOption(-1);
    formatter.PrintProlog();
    for (auto it = results->begin(); it != results->end(); ++it)
        formatter.PrintOneResultSet(**it, queries);
    formatter.PrintEpilog(options->GetOptions());
    return 0;
}
