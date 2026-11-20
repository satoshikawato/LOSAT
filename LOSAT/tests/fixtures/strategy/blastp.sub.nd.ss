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
                  data str ">qp2 unrelated protein (no hits expected)"
                }
              }
            },
            title "qp2 unrelated protein (no hits expected)"
          },
          inst {
            repr raw,
            mol aa,
            length 45,
            seq-data ncbieaa "IRRNWKMVCRYVPYLEEDMEAWMLWKEHLTFQWCNYFWTDVVKSE"
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
        value integer 3
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
