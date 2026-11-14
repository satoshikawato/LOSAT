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
                  data str ">nh1 random query without hits"
                }
              }
            },
            title "nh1 random query without hits"
          },
          inst {
            repr raw,
            mol aa,
            length 100,
            seq-data ncbieaa "RVRRTWHGTSYGERLFDVCYPRYGYATDCCHIYARMRWHTILSADRKQ
VDKMITLADWEPELPDAAHHCSPPQDWHKMDLMAQEFIEACRSGVHRTHFQP"
          }
        }
      }
    },
    subject database "tests/db_search/nophr",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -5 }
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
Blast4-request ::= {
  body queue-search {
    program "blastp",
    service "plain",
    queries bioseq-set {
      seq-set {
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
                  data str ">nh2 p01 of the database b\177ta"
                }
              }
            },
            title "nh2 p01 of the database b#ta"
          },
          inst {
            repr raw,
            mol aa,
            length 61,
            seq-data ncbieaa "QACACVTDIPCGIVVIDYVMQEVFDMKYWRYDLEDPRVMDLYCALFPM
WTFQRDLTCPTKD"
          }
        }
      }
    },
    subject database "tests/db_search/nophr",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -5 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 43967
      }
    }
  }
}
