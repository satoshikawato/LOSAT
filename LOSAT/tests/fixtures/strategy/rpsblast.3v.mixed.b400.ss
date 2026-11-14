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
                  data str ">x1"
                }
              }
            },
            title "x1"
          },
          inst {
            repr raw,
            mol aa,
            length 20,
            seq-data ncbieaa "XXXXXXXXXXXXXXXXXXXX"
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
            local str "Query_3"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">YP_009220476.1 envelope protein [White spot
 syndrome virus]"
                }
              }
            },
            title "YP_009220476.1 envelope protein [White spot syndrome virus]"
          },
          inst {
            repr raw,
            mol aa,
            length 208,
            seq-data ncbieaa "MHMWGVYAAILAGLTLILVVISIVVTNIELNKKLDKKDKDAYPVESEI
INLTINGVARGNHFNFVNGTLQTRNYGKVYVAGQGTSDSELVKKKGDIILTSLLGDGDHTLNVNKAESKELELYARVY
NNTKRDITVDSVSLSPGLNATGREFSANKFVLYFKPTVLKKNRINTLVFGATFDEDIDDTNRHYLLSMRFSPGNDLFK
VGEK"
          }
        },
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
                  data str ">YP_009220478.1 hypothetical protein [White spot
 syndrome virus]"
                }
              }
            },
            title "YP_009220478.1 hypothetical protein [White spot syndrome
 virus]"
          },
          inst {
            repr raw,
            mol aa,
            length 295,
            seq-data ncbieaa "MKNSRQRSGVWRGNSCLYKSFYFSGAIIECKKIRIIMMFLLLSLILFV
CFVGVVVGVIFMSRPNKTTTTSNKKTKKDKEKEKEDDTEGAVLGRREPENRPIGRDEEGAVEDGKEEEEVFEFEQPSV
NTGSNTGGGGTGTVPGEGLLPPPPPTPTPTPPPTPTPTPPPPPTRTPSPSSSLGEDDDDDIDIDFDDNDIDEFLDSGE
EMEEDEEEEDLDTLLSRLETGMSGEEVDFDASSAYIQPDPVVVKNIERSDYTLDPMESWKVLNRSEGDIRFFVDRGIT
NKIKAMTEDLKEL"
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
        value integer 1
      },
      {
        name "SmithWatermanMode",
        value boolean FALSE
      },
      {
        name "DbLength",
        value big-integer 560
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
        value big-integer 0
      }
    }
  }
}
