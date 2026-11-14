Blast4-request ::= {
  body queue-search {
    program "blastn",
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
                  data str ">nh1 random query without hits"
                }
              }
            },
            title "nh1 random query without hits"
          },
          inst {
            repr raw,
            mol na,
            length 200,
            seq-data ncbi2na 'FD74E43D0053B5B0ECA6032C053FC68A3140F5D7CF4A1705
E2B014ADDD6557CC027BE45C942F41A49E43A03290E1A3333C00'H
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
                  data str ">nh2 n02 of the database"
                }
              }
            },
            title "nh2 n02 of the database"
          },
          inst {
            repr raw,
            mol na,
            length 41,
            seq-data ncbi2na '02505B92AD05B1B26FF2C0'H
          }
        }
      }
    },
    subject database "tests/db_search/nnophr",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "MaskAtHash",
        value boolean TRUE
      },
      {
        name "UngappedMode",
        value boolean FALSE
      },
      {
        name "HitlistSize",
        value integer 500
      }
    }
  }
}
