Blast4-request ::= {
  body queue-search {
    program "blastx",
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
                  data str ">ovr_n1 override nucleotide"
                }
              }
            },
            title "ovr_n1 override nucleotide"
          },
          inst {
            repr raw,
            mol na,
            length 200,
            seq-data ncbi2na '049131660FBCB5C3BA3A0475B9C0240FF96E7DB02E3A1856
00C6A3A953A7B0E0A9EA2073E7A5660CE27A02D7CE0B38BFD758'H
          }
        }
      }
    },
    subject database "tests/db_search/dbmaskp/ps5",
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
        value integer 12
      },
      {
        name "WindowSize",
        value integer 40
      },
      {
        name "UngappedMode",
        value boolean FALSE
      },
      {
        name "QueryGeneticCode",
        value integer 1
      },
      {
        name "HitlistSize",
        value integer 500
      },
      {
        name "CompositionBasedStats",
        value integer 2
      },
      {
        name "SmithWatermanMode",
        value boolean FALSE
      }
    },
    program-options {
      {
        name "DbFilteringAlgorithmId",
        value integer 21
      },
      {
        name "SubjectMaskingType",
        value integer 1
      }
    }
  }
}
