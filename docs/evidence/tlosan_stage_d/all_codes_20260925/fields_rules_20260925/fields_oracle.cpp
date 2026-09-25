// Comparison-only pinned NCBI C++ API local-subject oracle for Stage D.
// Adapted from the Stage A database oracle. Never linked into LOSAT.
// NCBI source commit: 598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4.
// Stage D uses its internal search/trace state; output formatting is diagnostic
// until the local code-1 API path is calibrated against the pinned CLI.
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

// NCBI c++/src/algo/blast/blastinput/blast_args.cpp:2538-2557:
//   m_Scope = ReadSequencesToBlast(*subj_input_stream, IsProtein(), ...);
// This fixture-only reader supplies the same pinned FASTA letters to the API.
static std::pair<std::string, std::string> ReadFasta(const char* path) {
    std::ifstream in(path);
    std::string title, seq, line;
    if (!std::getline(in, title) || title.empty() || title[0] != '>')
        throw std::runtime_error("invalid FASTA header");
    while (std::getline(in, line)) seq += line;
    if (seq.empty()) throw std::runtime_error("empty FASTA sequence");
    return {title.substr(1, title.find(' ') - 1), seq};
}

// NCBI c++/src/app/blast/tblastn_app.cpp:274-280:
//   query = input->GetNextSeqBatch(*scope);
//   query_factory.Reset(new CObjMgr_QueryFactory(*query));
static CRef<CSeq_loc> AddSequence(CScope& scope, const char* path, bool protein) {
    auto fasta = ReadFasta(path);
    CRef<CBioseq> bio(new CBioseq);
    CRef<CSeq_id> id(new CSeq_id);
    id->SetLocal().SetStr(fasta.first);
    bio->SetId().push_back(id);
    std::ifstream title_in(path);
    std::string title;
    std::getline(title_in, title);
    CRef<CSeqdesc> desc(new CSeqdesc);
    desc->SetTitle(title.substr(1));
    bio->SetDescr().Set().push_back(desc);
    CSeq_inst& inst = bio->SetInst();
    inst.SetRepr(CSeq_inst::eRepr_raw);
    inst.SetMol(protein ? CSeq_inst::eMol_aa : CSeq_inst::eMol_na);
    inst.SetLength(fasta.second.size());
    if (protein) inst.SetSeq_data().SetIupacaa().Set(fasta.second);
    else inst.SetSeq_data().SetIupacna().Set(fasta.second);
    CRef<CSeq_entry> entry(new CSeq_entry);
    entry->SetSeq().Assign(*bio);
    scope.AddTopLevelSeqEntry(*entry);
    CRef<CSeq_loc> loc(new CSeq_loc);
    loc->SetWhole().Assign(*id);
    return loc;
}

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

    CRef<CScope> scope(new CScope(*CObjectManager::GetInstance()));
    scope->AddDefaults();
    CRef<CSeq_loc> qloc = AddSequence(*scope, argv[1], true);
    CRef<CSeq_loc> sloc = AddSequence(*scope, argv[2], false);
    CRef<CBlastQueryVector> queries(new CBlastQueryVector);
    queries->AddQuery(CRef<CBlastSearchQuery>(new CBlastSearchQuery(*qloc, *scope)));
    CRef<IQueryFactory> qf(new CObjMgr_QueryFactory(*queries));

    // NCBI c++/include/algo/blast/api/local_db_adapter.hpp:73-76:
    // CLocalDbAdapter(subject_sequences, opts_handle, dbscan_mode=false)
    // supplies the local -subject-like BlastSeqSrc, including its length.
    // NCBI c++/src/algo/blast/api/seqsrc_multiseq.cpp:140-164:
    // the stock local adapter reconstructs SSeqLoc; the comparison-only
    // BlastSeqSrc callback below supplies the selected genetic code.
    TSeqLocVector subject_locs;
    subject_locs.push_back(SSeqLoc(*sloc, *scope));
    CRef<IQueryFactory> sf(new CObjMgr_QueryFactory(subject_locs));
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
    CBlastFormat formatter(options->GetOptions(), *db, output, false, std::cout,
                           500, 250, *scope, "BLOSUM62", false, false,
                           1, code, options->GetOptions().GetSumStatisticsMode(), false, -1, "qseqid sseqid score bitscore evalue nident positive length mismatch gaps gapopen qstart qend sstart send sframe");
    formatter.SetHitsSortOption(-1);
    formatter.SetHspsSortOption(-1);
    formatter.PrintProlog();
    for (auto it = results->begin(); it != results->end(); ++it)
        formatter.PrintOneResultSet(**it, queries);
    formatter.PrintEpilog(options->GetOptions());
    return 0;
}
