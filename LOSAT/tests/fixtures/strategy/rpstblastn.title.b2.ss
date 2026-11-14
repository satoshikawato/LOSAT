Blast4-request ::= {
  body queue-search {
    program "tblastn",
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
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "QueryGeneticCode",
        value integer 1
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
        name "UngappedMode",
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
Blast4-request ::= {
  body queue-search {
    program "tblastn",
    service "rpsblast",
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
                  data str ">t2 second batch b\177ta"
                }
              }
            },
            title "t2 second batch b#ta"
          },
          inst {
            repr raw,
            mol na,
            length 366,
            seq-data ncbi2na 'BA3A2B302B7A7D74CB3BA3A244BB1D225BFD92F29802D3A1
F5500D223AD5325B0E0A9EA2053EFA5640CED7A02FF0E0B383FDF484F457E8DEFA3D5020410494
AFAE821A311D4CE9E03AAF3A3930B1502ECEEEE3AF456A935B47DA5450'H
          }
        }
      }
    },
    subject database "rps_mini",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "QueryGeneticCode",
        value integer 1
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
        name "UngappedMode",
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
      },
      {
        name "MatrixName",
        value string "BLOSUM62"
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
        name "EffectiveSearchSpace",
        value big-integer 48609
      }
    }
  }
}
