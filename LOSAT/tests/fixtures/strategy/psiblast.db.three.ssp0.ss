Blast4-request ::= {
  body queue-search {
    program "blastp",
    service "psi",
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
                  data str ">YP_009220601.1 TK-TMK, chimeric thymidine and
 thymidylate kinase [White spot syndrome virus]"
                }
              }
            },
            title "YP_009220601.1 TK-TMK, chimeric thymidine and thymidylate
 kinase [White spot syndrome virus]"
          },
          inst {
            repr raw,
            mol aa,
            length 398,
            seq-data ncbieaa "MQLILSHHLTMAGRVELVTGPMFAGKSTYLKNIYQQENGGNKHCLFVK
HSLETRYGCGTGTIVTHAGEVIEGCTTVSSIKELISVLPEVVDVILIDEGQFFTDLVLVNRLADKGKRIVIAALDGTS
DQQMFSPIHKLLPYTNSIVKLASKCMICKIDTKEAPFTVRFGNDNDNNVICVGGAEMYAAACRDCYKKINKKKNKGKL
VVLEGGDRCGKSTQAKLLLTNKNSPLYGGEYMCFPDRSSHTGKLINDYLTKKIELDDHAAHLLFSANRWEVCSKIKQL
LDDGIHVVMDRYYYSGIVFSLARGVDTVEWCSASDEGLPQPDLVLLMLLDVEKCSNRDTFGVERFETNSIQERARALF
LDLANKDEKNVWIKVDARGTIEEVQTKIINIVYNIVEE"
          }
        }
      }
    },
    subject database "tests/db_search/paja",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "InclusionThreshold",
        value real { 2, 10, -3 }
      },
      {
        name "PseudoCountWeight",
        value integer 0
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 0
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
        name "CompositionBasedStats",
        value integer 2
      },
      {
        name "SmithWatermanMode",
        value boolean FALSE
      },
      {
        name "GapTrigger",
        value real { 22, 10, 0 }
      },
      {
        name "IgnoreMsaMaster",
        value boolean FALSE
      }
    },
    format-options {
      {
        name "Web_StepNumber",
        value integer 1
      }
    }
  }
}
Blast4-request ::= {
  body queue-search {
    program "blastp",
    service "psi",
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
                  data str ">YP_009220505.1 hypothetical protein [White spot
 syndrome virus]"
                }
              }
            },
            title "YP_009220505.1 hypothetical protein [White spot syndrome
 virus]"
          },
          inst {
            repr raw,
            mol aa,
            length 76,
            seq-data ncbieaa "MISLIFFLSRPEALTRVQPDIWRGNSPVMISLIFFLSRPEALTRVQPD
IWRGNSPVMISPSHPFPPFLGSGVFTCV"
          }
        }
      }
    },
    subject database "tests/db_search/paja",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "InclusionThreshold",
        value real { 2, 10, -3 }
      },
      {
        name "PseudoCountWeight",
        value integer 0
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 20272584
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
        name "CompositionBasedStats",
        value integer 2
      },
      {
        name "SmithWatermanMode",
        value boolean FALSE
      },
      {
        name "GapTrigger",
        value real { 22, 10, 0 }
      },
      {
        name "IgnoreMsaMaster",
        value boolean FALSE
      }
    },
    format-options {
      {
        name "Web_StepNumber",
        value integer 1
      }
    }
  }
}
Blast4-request ::= {
  body queue-search {
    program "blastp",
    service "psi",
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
                  data str ">qnone no hits"
                }
              }
            },
            title "qnone no hits"
          },
          inst {
            repr raw,
            mol aa,
            length 36,
            seq-data ncbieaa "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
          }
        }
      }
    },
    subject database "tests/db_search/paja",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "InclusionThreshold",
        value real { 2, 10, -3 }
      },
      {
        name "PseudoCountWeight",
        value integer 0
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 1688148
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
        name "CompositionBasedStats",
        value integer 2
      },
      {
        name "SmithWatermanMode",
        value boolean FALSE
      },
      {
        name "GapTrigger",
        value real { 22, 10, 0 }
      },
      {
        name "IgnoreMsaMaster",
        value boolean FALSE
      }
    },
    format-options {
      {
        name "Web_StepNumber",
        value integer 1
      }
    }
  }
}
Blast4-request ::= {
  body queue-search {
    program "blastp",
    service "psi",
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
                  data str ">rnd random"
                }
              }
            },
            title "rnd random"
          },
          inst {
            repr raw,
            mol aa,
            length 30,
            seq-data ncbieaa "MWWHPKQCYRWNMHFWQCPYRGWTMHKWYQ"
          }
        }
      }
    },
    subject database "tests/db_search/paja",
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
      },
      {
        name "InclusionThreshold",
        value real { 2, 10, -3 }
      },
      {
        name "PseudoCountWeight",
        value integer 0
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 1644800
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
        name "CompositionBasedStats",
        value integer 2
      },
      {
        name "SmithWatermanMode",
        value boolean FALSE
      },
      {
        name "GapTrigger",
        value real { 22, 10, 0 }
      },
      {
        name "IgnoreMsaMaster",
        value boolean FALSE
      }
    },
    format-options {
      {
        name "Web_StepNumber",
        value integer 1
      }
    }
  }
}
