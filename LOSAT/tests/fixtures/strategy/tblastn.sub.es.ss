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
                  data str ">e empty"
                }
              }
            },
            title "e empty"
          },
          inst {
            repr virtual,
            mol aa,
            length 0
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
                  data str ">s two"
                }
              }
            },
            title "s two"
          },
          inst {
            repr raw,
            mol aa,
            length 2,
            seq-data ncbieaa "AC"
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
                data str ">sn1 plus strand gene"
              }
            }
          },
          title "sn1 plus strand gene"
        },
        inst {
          repr raw,
          mol na,
          length 501,
          seq-data ncbi2na '4842B5FEAA0C65892B0D3C3640C94D3BA3A04C9B9C025BBF96
E29A0AC3A38D600C661AD4826BCE0AAEA2071E7A5E60CEDFA0ADF8E1B18BF24B8CC457EBAEFA3D
40A0492155ABAE8012117E44E9603A873A3910B1580E4EEED3AE6B31A34B5F0847A88F0C58650F
D88D74031C4533ACA13AD91E7260AFD4E31AE073C07240'H
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
                data str ">sn2 minus strand gene"
              }
            }
          },
          title "sn2 minus strand gene"
        },
        inst {
          repr raw,
          mol na,
          length 471,
          seq-data ncbi2na '289DCE1129C085F4FA59DE54FA45EB0A039FD0FDD80F12374F
CDF52E5E2E6E6CBF668533A44CA764393530553F4744DCBBB7BB245256B1653B3F6835246D4FAF
10ABA31677A074C7D38A1BB9D533F9CF8BB2CC010CA3D5317D34594EAD4D983FCAB4D4D17F689E
41E3797E2847B9D4D4B2781827549C'H
        }
      },
      {
        id {
          local str "Subject_3"
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">sn3 no hits"
              }
            }
          },
          title "sn3 no hits"
        },
        inst {
          repr raw,
          mol na,
          length 150,
          seq-data ncbi2na 'F7C61384618166375B6669F3CCF5515B9D8BE60EB8B1BD1979
FC318617C8E60C4F205DFC8140'H
        }
      }
    },
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
    }
  }
}
