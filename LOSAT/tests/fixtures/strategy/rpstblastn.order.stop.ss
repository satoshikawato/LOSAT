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
                  data str ">q1 one\303\251"
                }
              }
            },
            title "q1 one##"
          },
          inst {
            repr raw,
            mol na,
            length 180,
            seq-data ncbi2na '3A33EFEB717FDE76CF4BF9D2CEF59600F3DBFC5E930D35B7
7A0AFA0EB6D2FF38CFAE00C5ACEFDB7ABFD603F3EB'H
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
                  data str ">q2 two"
                }
              }
            },
            title "q2 two"
          },
          inst {
            repr raw,
            mol na,
            length 179,
            seq-data ncbi2na '3A33EFEB717FDE7B3D2FE74B3BD65803CF6FF17A4C34D6DD
E82BE83ADB4BFCE33EB80316B3BF6DEAFF580FCFAC'H
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
        value big-integer 21252
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
            local str "Query_3"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">q3 thr\303\251e"
                }
              }
            },
            title "q3 thr##e"
          },
          inst {
            repr raw,
            mol na,
            length 179,
            seq-data ncbi2na '3A33EFEB717FDE7B3D2FE74B3BD65803CF6FF17A4C34D6DD
E82BE83ADB4BFCE33EB80316B3BF6DEAFF580FCFAC'H
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
        value big-integer 20746
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
            local str "Query_4"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">q4 four\303\251"
                }
              }
            },
            title "q4 four##"
          },
          inst {
            repr virtual,
            mol na,
            length 0
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
        value big-integer 20746
      }
    }
  }
}
