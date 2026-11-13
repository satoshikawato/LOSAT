Blast4-request ::= {
  body queue-search {
    program "tblastn",
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
                  data str ">ql1 lower runs"
                }
              }
            },
            title "ql1 lower runs"
          },
          inst {
            repr raw,
            mol aa,
            length 140,
            seq-data ncbieaa "MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKV
LMKYEFSSDIHPCNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHKYPKCMCVMVRPTDQSLRHWKIMHDANFEKLKI
FHYMVGHGRTASEG"
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
                  data str ">ql2 no hits"
                }
              }
            },
            title "ql2 no hits"
          },
          inst {
            repr raw,
            mol aa,
            length 45,
            seq-data ncbieaa "IRRNWKMVCRYVPYLEEDMEAWMLWKEHLTFQWCNYFWTDVVKSE"
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
                  data str ">ql3 all lower"
                }
              }
            },
            title "ql3 all lower"
          },
          inst {
            repr raw,
            mol aa,
            length 50,
            seq-data ncbieaa "MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKV
LM"
          }
        }
      }
    },
    subject database "tests/db_search/arch_tbn_snucl",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "DbGeneticCode",
        value integer 1
      },
      {
        name "UngappedMode",
        value boolean FALSE
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
        name "LCaseMask",
        value query-mask {
          locations {
            packed-int {
              {
                from 5,
                to 14,
                id local str "Query_1"
              },
              {
                from 30,
                to 40,
                id local str "Query_1"
              },
              {
                from 100,
                to 109,
                id local str "Query_1"
              }
            }
          },
          frame notset
        }
      },
      {
        name "LCaseMask",
        value query-mask {
          locations {
            packed-int {
              {
                from 0,
                to 49,
                id local str "Query_3"
              }
            }
          },
          frame notset
        }
      }
    }
  }
}
