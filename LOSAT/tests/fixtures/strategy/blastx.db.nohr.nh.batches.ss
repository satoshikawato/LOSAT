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
                  data str ">nh1 random query without hits"
                }
              }
            },
            title "nh1 random query without hits"
          },
          inst {
            repr raw,
            mol na,
            length 200,
            seq-data ncbi2na 'FD74E43D0053B5B0ECA6032C053FC68A3140F5D7CF4A1705
E2B014ADDD6557CC027BE45C942F41A49E43A03290E1A3333C00'H
          }
        }
      }
    },
    subject database "tests/db_search/nophr",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -3 }
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
    }
  }
}
Blast4-request ::= {
  body queue-search {
    program "blastx",
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
                  data str ">nh2 p01 of the database, back-translated"
                }
              }
            },
            title "nh2 p01 of the database, back-translated"
          },
          inst {
            repr raw,
            mol na,
            length 183,
            seq-data ncbi2na '427EE7EEF1E33D7EEB3EFBCF8F3BCE420BFF8CE033E9BCE3
7A08D76EF3A37B3EE77BF5CEE87FD06E3787ED71C08C'H
          }
        }
      }
    },
    subject database "tests/db_search/nophr",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, -3 }
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 27646
      }
    }
  }
}
