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
                  data str ">f1144_q"
                }
              }
            },
            title "f1144_q"
          },
          inst {
            repr raw,
            mol na,
            length 55,
            seq-data ncbi2na 'B745CFBCF8697EB13B182238E4A4'H
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
                data str ">f1144_s"
              }
            }
          },
          title "f1144_s"
        },
        inst {
          repr raw,
          mol aa,
          length 30,
          seq-data ncbieaa "VYPIKCFQIQFVNFFASNEHHFISAMDPWP"
        }
      }
    },
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 5 }
      },
      {
        name "WordSize",
        value integer 3
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
        name "WordThreshold",
        value integer 2
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
