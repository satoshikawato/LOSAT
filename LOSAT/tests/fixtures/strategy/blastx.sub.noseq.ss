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
                  data str ">qn1 plus strand gene"
                }
              }
            },
            title "qn1 plus strand gene"
          },
          inst {
            repr raw,
            mol na,
            length 366,
            seq-data ncbi2na 'BA3A2B302B7A7D74CB3BA3A244BB1D225BFD92F29802D3A1
F5500D223AD5325B0E0A9EA2053EFA5640CED7A02FF0E0B383FDF484F457E8DEFA3D5020410494
AFAE821A311D4CE9E03AAF3A3930B1502ECEEEE3AF456A935B47DA5450'H
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
                  data str ">qn2 minus strand gene"
                }
              }
            },
            title "qn2 minus strand gene"
          },
          inst {
            repr raw,
            mol na,
            length 378,
            seq-data ncbi2na '73E78343DC2BF374F8FF50E9A2A4E6EBF56853B844CAF643
9352255376706E0BB35B2A4525AA1A53B1FEB35046D4FAF11AB83D68780F6CFD3C21B89E133C9C
C86F08E092CE1D531FD38694DAD4DDE37DA34D4C17D7C1E894B79EC2518EBA20'H
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
                  data str ">qn3 no hits"
                }
              }
            },
            title "qn3 no hits"
          },
          inst {
            repr raw,
            mol na,
            length 150,
            seq-data ncbi2na '92A40D92FF70C62172A881B47AA438756AE2584B5C3AB4C2
37BE479DF9F439282377E442F6B0'H
          }
        },
        seq {
          id {
            local str "Query_4"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">noseq a record without letters"
                }
              }
            },
            title "noseq a record without letters"
          },
          inst {
            repr virtual,
            mol na,
            length 0
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
