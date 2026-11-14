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
                  data str ">q1 one\303\251"
                }
              }
            },
            title "q1 one##"
          },
          inst {
            repr raw,
            mol aa,
            length 60,
            seq-data ncbieaa "MDIVWSTFSARIQFAQYVPPKIIRFTWHNHPSLEGWNGRQFYDIGENT
GMFRLGFPKFIG"
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
        value integer 1
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
    program "blastp",
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
            mol aa,
            length 60,
            seq-data ncbieaa "MDIVWSTFSAOIQFAQYVPPKIIRFTWHNHPSLEGWNGRQFYDIGENT
GMFRLGFPKFIG"
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
        value integer 1
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
    program "blastp",
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
            mol aa,
            length 60,
            seq-data ncbieaa "MDIVWSTFSAOIQFAQYVPPKIIRFTWHNHPSLEGWNGRQFYDIGENT
GMFRLGFPKFIG"
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
        value integer 1
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
    program "blastp",
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
                  data str ">q4 four\303\251
 MKVLAAGIVAMKVLAAGIVAMKVLAAGIVAMKVLAAGIVAMKVLAAGIVAMKVLAAGIVA"
                }
              }
            },
            title "q4 four##
 MKVLAAGIVAMKVLAAGIVAMKVLAAGIVAMKVLAAGIVAMKVLAAGIVAMKVLAAGIVA"
          },
          inst {
            repr virtual,
            mol aa,
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
        value integer 1
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
