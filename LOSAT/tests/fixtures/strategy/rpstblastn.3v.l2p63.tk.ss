Blast4-request ::= {
  body queue-search {
    program "tblastn",
    service "rpsblast",
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
                  data str ">tkcut570"
                }
              }
            },
            title "tkcut570"
          },
          inst {
            repr raw,
            mol na,
            length 570,
            seq-data ncbi2na '39278C78B4D1F873A7A1AB6272F16B54EF65A80DC5C7E003
3F3410801A290C04F97FFB024F57E0118CE8EE812818CB05466A60B0F8AAEC612E2F734081C3F7
B5D522B2E8EC37C3638284BFF448DFB1FB432BE92102B00832E3E591F86A1778D248EFC954C440
F9D5331C1D0FB80726D82E4E37908F8C6022915FC6B08F6B063061041B0CE6CA6B9A23B19A4939
6210'H
          }
        }
      }
    },
    subject database "l2p63rps3v",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "QueryGeneticCode",
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
        name "UngappedMode",
        value boolean FALSE
      },
      {
        name "HitlistSize",
        value integer 500
      },
      {
        name "CompositionBasedStats",
        value integer 0
      },
      {
        name "SmithWatermanMode",
        value boolean FALSE
      }
    }
  }
}
