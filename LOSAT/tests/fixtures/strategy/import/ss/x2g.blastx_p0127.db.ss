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
                  data str ">p0127_q"
                }
              }
            },
            title "p0127_q"
          },
          inst {
            repr raw,
            mol na,
            length 163,
            seq-data ncbi2na '3754225FE1BAABE1298594F26693A575662D9584FA9B86B5
1EF9CE3C09758E11B09A47B0EA1ED47ACC'H
          }
        }
      }
    },
    subject database "tests/db_search/x2g_p0127",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 5 }
      },
      {
        name "WordThreshold",
        value integer 1
      },
      {
        name "WordSize",
        value integer 5
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
        name "SegFiltering",
        value boolean FALSE
      },
      {
        name "WindowSize",
        value integer 0
      },
      {
        name "UngappedMode",
        value boolean TRUE
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
        value integer 0
      },
      {
        name "SmithWatermanMode",
        value boolean FALSE
      }
    }
  }
}
