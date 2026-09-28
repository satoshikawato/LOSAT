// Comparison-only NCBI C++ API oracle. This file is not part of LOSAT builds.
#include <corelib/ncbistd.hpp>
#include <algo/blast/api/blast_aux.hpp>
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
    desc->SetTitle(protein ? title.substr(1) : title.substr(title.find(' ') + 1));
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

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "usage: oracle QUERY.faa SUBJECT.fna DATABASE CODE OUTFMT\n";
        return 2;
    }
    const int code = std::stoi(argv[4]);
    const int fmt = std::stoi(argv[5]);
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

    // Keep the identical subject sequence in the formatter scope; the search
    // itself reads the separately built nucleotide BLAST database.
    (void)sloc;
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
    // NCBI c++/src/app/blast/tblastn_app.cpp:189-193,287-291:
    //   InitializeSubject(db_args, m_OptsHndl, ..., db_adapter, scope);
    //   CLocalBlast lcl_blast(query_factory, m_OptsHndl, db_adapter);
    CSearchDatabase search_db(argv[3], CSearchDatabase::eBlastDbIsNucleotide);
    CRef<CLocalDbAdapter> db(new CLocalDbAdapter(search_db));

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
                           1, code, options->GetOptions().GetSumStatisticsMode());
    formatter.SetHitsSortOption(-1);
    formatter.SetHspsSortOption(-1);
    formatter.PrintProlog();
    for (auto it = results->begin(); it != results->end(); ++it)
        formatter.PrintOneResultSet(**it, queries);
    formatter.PrintEpilog(options->GetOptions());
    return 0;
}
