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
                  data str ">h4"
                }
              }
            },
            title "h4"
          },
          inst {
            repr raw,
            mol aa,
            length 80,
            seq-data ncbieaa "MTSPAPSPSSTPKSSCTTIVNRCGFLLDNNKEVVIYDTNSKFKCEPKN
LELIGVLSGVSDNVVTQISPDQIFVGTYMVKY"
          }
        }
      }
    },
    subject database "rps_mini3v",
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
        value integer 0
      },
      {
        name "SmithWatermanMode",
        value boolean FALSE
      }
    }
  }
}
