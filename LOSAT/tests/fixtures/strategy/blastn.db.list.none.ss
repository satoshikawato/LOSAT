Blast4-request ::= {
  body queue-search {
    program "blastn",
    service "megablast",
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
                  data str ">first nucleotide query with a title long enough
 to be cut"
                }
              }
            },
            title "first nucleotide query with a title long enough to be cut"
          },
          inst {
            repr raw,
            mol na,
            length 429,
            seq-data ncbi2na '3AE7B75E5842141B42965EAB0ADA66467A62CEB9A295E88A
3BD7B5F55145085C7D591F61789469DE54AF0A9469082BA5866785065BA646E8613950667B595E
26179199109F6AE856B41F4275C251E5E7AE15E96545D5658BD1579AE465D5E842F5E9F7B8916E
785D40C5BC00'H
          }
        }
      }
    },
    subject database "tests/db_search/taxfilter/tfn",
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
    },
    program-options {
      {
        name "TaxidList",
        value integer-list {
        }
      }
    }
  }
}
Blast4-request ::= {
  body queue-search {
    program "blastn",
    service "megablast",
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
                  data str ">second"
                }
              }
            },
            title "second"
          },
          inst {
            repr raw,
            mol na,
            length 444,
            seq-data ncbi2na '3AE4DE1D78A20B796F1E57BAA42B81BA382FAEB8A57A929E
7AEB715FA1522BDFE2D7FAA37B51D78E7BCEA415C2B829D3A4202E76B97F2E3A5E9D17A105D0A9
17F94478B89E47B842791BA35E207D29D7A906E7ADEEE7A54D1FE9020F45514B929E5CD202EBA7
AEE9C395E9510B3470'H
          }
        }
      }
    },
    subject database "tests/db_search/taxfilter/tfn",
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 0
      }
    },
    program-options {
      {
        name "TaxidList",
        value integer-list {
        }
      }
    }
  }
}
