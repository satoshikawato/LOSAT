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
                  data str ">b1 sp1 plus"
                }
              }
            },
            title "b1 sp1 plus"
          },
          inst {
            repr raw,
            mol na,
            length 108,
            seq-data ncbi2na 'BBF9AD20802E3A185608C6A1AD43A6BCE0A9EA20F1E7A524
0CEDBA'H
          }
        },
        seq {
          id {
            local str "Query_2"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">b2 sp2 minus"
                }
              }
            },
            title "b2 sp2 minus"
          },
          inst {
            repr raw,
            mol na,
            length 108,
            seq-data ncbi2na '33F4C244EE474CFE739CE4A3004308754C7F4F1A53AA5376
0DF28D'H
          }
        }
      }
    },
    subject database "tests/db_search/arch_tbx_snucl",
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
        name "MatrixName",
        value string "blosum62"
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
    }
  }
}
