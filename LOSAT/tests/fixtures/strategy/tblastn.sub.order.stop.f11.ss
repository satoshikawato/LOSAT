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
                data str ">s1 subject"
              }
            }
          },
          title "s1 subject"
        },
        inst {
          repr raw,
          mol na,
          length 180,
          seq-data ncbi2na '3A33EFEB717FDE76CF4BF9D2CEF59600F3DBFC5E930D35B77A
0AFA0EB6D2FF38CFAE00C5ACEFDB7ABFD603F3EB'H
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
                data str ">s2 other"
              }
            }
          },
          title "s2 other"
        },
        inst {
          repr raw,
          mol na,
          length 150,
          seq-data ncbi2na '1FEDA9F6CDE9F4E7DA90DF0756B028E460FB8DA88A5764C262
B12BABE6FDFC57C0DC984C62A0'H
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
Blast4-request ::= {
  body queue-search {
    program "tblastn",
    service "plain",
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
                data str ">s1 subject"
              }
            }
          },
          title "s1 subject"
        },
        inst {
          repr raw,
          mol na,
          length 180,
          seq-data ncbi2na '3A33EFEB717FDE76CF4BF9D2CEF59600F3DBFC5E930D35B77A
0AFA0EB6D2FF38CFAE00C5ACEFDB7ABFD603F3EB'H
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
                data str ">s2 other"
              }
            }
          },
          title "s2 other"
        },
        inst {
          repr raw,
          mol na,
          length 150,
          seq-data ncbi2na '1FEDA9F6CDE9F4E7DA90DF0756B028E460FB8DA88A5764C262
B12BABE6FDFC57C0DC984C62A0'H
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 5088
      }
    }
  }
}
Blast4-request ::= {
  body queue-search {
    program "tblastn",
    service "plain",
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
                data str ">s1 subject"
              }
            }
          },
          title "s1 subject"
        },
        inst {
          repr raw,
          mol na,
          length 180,
          seq-data ncbi2na '3A33EFEB717FDE76CF4BF9D2CEF59600F3DBFC5E930D35B77A
0AFA0EB6D2FF38CFAE00C5ACEFDB7ABFD603F3EB'H
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
                data str ">s2 other"
              }
            }
          },
          title "s2 other"
        },
        inst {
          repr raw,
          mol na,
          length 150,
          seq-data ncbi2na '1FEDA9F6CDE9F4E7DA90DF0756B028E460FB8DA88A5764C262
B12BABE6FDFC57C0DC984C62A0'H
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 5088
      }
    }
  }
}
Blast4-request ::= {
  body queue-search {
    program "tblastn",
    service "plain",
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
                data str ">s1 subject"
              }
            }
          },
          title "s1 subject"
        },
        inst {
          repr raw,
          mol na,
          length 180,
          seq-data ncbi2na '3A33EFEB717FDE76CF4BF9D2CEF59600F3DBFC5E930D35B77A
0AFA0EB6D2FF38CFAE00C5ACEFDB7ABFD603F3EB'H
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
                data str ">s2 other"
              }
            }
          },
          title "s2 other"
        },
        inst {
          repr raw,
          mol na,
          length 150,
          seq-data ncbi2na '1FEDA9F6CDE9F4E7DA90DF0756B028E460FB8DA88A5764C262
B12BABE6FDFC57C0DC984C62A0'H
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 5088
      }
    }
  }
}
