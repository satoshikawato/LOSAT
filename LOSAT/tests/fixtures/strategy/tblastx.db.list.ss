Blast4-request ::= {
  body queue-search {
    program "tblastx",
    service "plain",
    queries bioseq-set {
      seq-set {
        seq {
          id {
            local str "Query_1"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">ENST00000320868.9 cdna
 chromosome:GRCh38:16:176680:177522:1 gene:ENSG00000206172.10
 gene_biotype:protein_coding transcript_biotype:protein_coding
 gene_symbol:HBA1 description:hemoglobin subunit alpha 1 [Source:HGNC
 Symbol;Acc:HGNC:4823]"
                }
              }
            },
            title "ENST00000320868.9 cdna chromosome:GRCh38:16:176680:177522:1
 gene:ENSG00000206172.10 gene_biotype:protein_coding
 transcript_biotype:protein_coding gene_symbol:HBA1 description:hemoglobin
 subunit alpha 1 [Source:HGNC Symbol;Acc:HGNC:4823]"
          },
          inst {
            repr raw,
            mol na,
            length 577,
            seq-data ncbi2na '1DF7AD5448748881514EB9EDD796108506D0A597AAC2B699
919E98B3AE68A57A228EF5ED7D545142171F5647D85E251A77952BC2A51A420AE96199E14196E9
91BA184E54199ED6578985E4664427DABA15AD07D09D70947979EB857A5951755962F455E6B919
757A10BD7A7DEE245B9E1750316F09E8976BA539F7E55FA97555255D755F5E456C556EB7F8302D
E2EA6900'H
          }
        }
      }
    },
    subject database "tests/db_search/taxfilter/tfn",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "LongestIntronLength",
        value integer 0
      },
      {
        name "MaskAtHash",
        value boolean FALSE
      },
      {
        name "SegFilteringWindow",
        value integer 12
      },
      {
        name "SegFilteringLocut",
        value real { 22, 10, -1 }
      },
      {
        name "SegFilteringHicut",
        value real { 25, 10, -1 }
      },
      {
        name "WordThreshold",
        value integer 13
      },
      {
        name "WindowSize",
        value integer 40
      },
      {
        name "QueryGeneticCode",
        value integer 1
      },
      {
        name "DbGeneticCode",
        value integer 1
      },
      {
        name "HitlistSize",
        value integer 500
      }
    },
    program-options {
      {
        name "NegativeTaxidList",
        value integer-list {
          9443
        }
      }
    }
  }
}
