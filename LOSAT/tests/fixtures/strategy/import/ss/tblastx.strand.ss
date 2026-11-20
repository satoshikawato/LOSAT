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
                  data str ">qx1 plus strand gene"
                }
              }
            },
            title "qx1 plus strand gene"
          },
          inst {
            repr raw,
            mol na,
            length 129,
            seq-data ncbi2na '17C3A9B3D92D21A0AC3A185408D6E3A553A7BCE02BEA00F1
E7A50404EEFA7BC4C0'H
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
                  data str ">qx2 minus strand gene"
                }
              }
            },
            title "qx2 minus strand gene"
          },
          inst {
            repr raw,
            mol na,
            length 129,
            seq-data ncbi2na '6D8BF231F4D206E6674EFA639CECA3824338F54C5F4F1E53
EB5B560FF6AD058140'H
          }
        },
        seq {
          id {
            local str "Query_3"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">qx3 no hits"
                }
              }
            },
            title "qx3 no hits"
          },
          inst {
            repr raw,
            mol na,
            length 90,
            seq-data ncbi2na '9050B43EF2DFF11020771ACE787B7D4758BA590B1F7570'H
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
    }
  }
}
