Blast4-request ::= {
  body queue-search {
    program "blastp",
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
                  data str ">q1 one\303\251
 MKVLAAGIVAMKVLAAGIVAMKVLAAGIVAMKVLAAGIVAMKVLAAGIVA"
                }
              }
            },
            title "q1 one## MKVLAAGIVAMKVLAAGIVAMKVLAAGIVAMKVLAAGIVAMKVLAAGIVA"
          },
          inst {
            repr raw,
            mol aa,
            length 60,
            seq-data ncbieaa "MDIVWSTFSARIQFAQYVPPKIIRFTWHNHPSLEGWNGRQFYDIGENT
GMFRLGFPKFIG"
          }
        }
      }
    },
    subject database "tests/db_search/arch_sprot",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "GapOpeningCost",
        value integer 11
      },
      {
        name "GapExtensionCost",
        value integer 1
      },
      {
        name "MaskAtHash",
        value boolean FALSE
      },
      {
        name "SegFiltering",
        value boolean FALSE
      },
      {
        name "MatrixName",
        value string "blosum62"
      },
      {
        name "WordThreshold",
        value integer 11
      },
      {
        name "WindowSize",
        value integer 40
      },
      {
        name "HitlistSize",
        value integer 500
      },
      {
        name "UngappedMode",
        value boolean FALSE
      },
      {
        name "CompositionBasedStats",
        value integer 2
      },
      {
        name "SmithWatermanMode",
        value boolean FALSE
      }
    }
  }
}
