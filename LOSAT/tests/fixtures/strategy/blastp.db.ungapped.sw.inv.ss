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
                  data str ">qinv invalid letters"
                }
              }
            },
            title "qinv invalid letters"
          },
          inst {
            repr raw,
            mol aa,
            length 60,
            seq-data ncbieaa "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXX"
          }
        }
      }
    },
    subject database "tests/db_search/x16_taxid",
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
        value integer 500
      },
      {
        name "UngappedMode",
        value boolean TRUE
      },
      {
        name "CompositionBasedStats",
        value integer 0
      },
      {
        name "SmithWatermanMode",
        value boolean TRUE
      }
    }
  }
}
