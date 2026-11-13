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
                  data str ">b1 sp1 plus"
                }
              }
            },
            title "b1 sp1 plus"
          },
          inst {
            repr raw,
            mol na,
            length 114,
            seq-data ncbi2na '3BA3A04C9B9C0A4B7F9AF6180AF3A38D400D2A3A17393B0E
028EA20730'H
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
            length 105,
            seq-data ncbi2na '024307F5445F4C1A53EA5B560FFAADB53C5FDF0F90FAD51F
4A1A00'H
          }
        }
      }
    },
    subject sequences {
      {
        id {
          local str "Subject_1"
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">sp1 subject protein alpha"
              }
            }
          },
          title "sp1 subject protein alpha"
        },
        inst {
          repr raw,
          mol aa,
          length 137,
          seq-data ncbieaa "MWMEHSVLKAVFAVRRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLM
KYEFSSDIHPWICWIPKEQQTPGWWKTDTLHMAEMGLWMHKYPKCMCVMVRYTDQSLRHWEINTDANFEKLKILHHMV
GHGRTASEG"
        }
      },
      {
        id {
          local str "Subject_2"
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">sp2 subject protein beta, distant"
              }
            }
          },
          title "sp2 subject protein beta, distant"
        },
        inst {
          repr raw,
          mol aa,
          length 142,
          seq-data ncbieaa "YWMEHSVLKADQLQLRKVMDDPKIRDGPMAVMKVWESYCFILHSIANMSA
HVLMKYEFQSGIHPCNQWTCWIPKYMAVPGWCPTDTLHVSEMGLWMHCEPMCHMVRETTHALRHWKINEICNFEKLKA
FPYQVPMGRAANEG"
        }
      }
    },
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
