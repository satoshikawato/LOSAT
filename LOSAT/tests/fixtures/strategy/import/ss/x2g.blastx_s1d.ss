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
                  data str ">c02074_q"
                }
              }
            },
            title "c02074_q"
          },
          inst {
            repr raw,
            mol na,
            length 15,
            seq-data ncbi2na '84E99430'H
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
                data str ">c02074_s"
              }
            }
          },
          title "c02074_s"
        },
        inst {
          repr raw,
          mol aa,
          length 10,
          seq-data ncbieaa "FIAFIFIAFI"
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
        value integer 3
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
