Blast4-request ::= {
  body queue-search {
    program "blastp",
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
                  data str ">qp1 query protein one"
                }
              }
            },
            title "qp1 query protein one"
          },
          inst {
            repr raw,
            mol aa,
            length 140,
            seq-data ncbieaa "MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKV
LMKYEFSSDIHPCNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHKYPKCMCVMVRPTDQSLRHWKIMHDANFEKLKI
FHYMVGHGRTASEG"
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
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
        name "HitlistSize",
        value integer 500
      },
      {
        name "CompositionBasedStats",
        value integer 1
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
    service "rpsblast",
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
                  data str ">t2 second batch b\177ta"
                }
              }
            },
            title "t2 second batch b#ta"
          },
          inst {
            repr raw,
            mol aa,
            length 140,
            seq-data ncbieaa "MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKV
LMKYEFSSDIHPCNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHKYPKCMCVMVRPTDQSLRHWKIMHDANFEKLKI
FHYMVGHGRTASEG"
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
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
        name "HitlistSize",
        value integer 500
      },
      {
        name "CompositionBasedStats",
        value integer 1
      },
      {
        name "SmithWatermanMode",
        value boolean FALSE
      },
      {
        name "MatrixName",
        value string "BLOSUM62"
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
        name "EffectiveSearchSpace",
        value big-integer 55775
      }
    }
  }
}
