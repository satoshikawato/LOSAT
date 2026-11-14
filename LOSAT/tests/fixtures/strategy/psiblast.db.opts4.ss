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
        name "DbLength",
        value big-integer 1000000
      },
      {
        name "GapOpeningCost",
        value integer 14
      },
      {
        name "GapExtensionCost",
        value integer 2
      },
      {
        name "MaskAtHash",
        value boolean TRUE
      },
      {
        name "SegFilteringWindow",
        value integer 1
      },
      {
        name "SegFilteringLocut",
        value real { 15, 10, -1 }
      },
      {
        name "SegFilteringHicut",
        value real { 18, 10, -1 }
      },
      {
        name "MatrixName",
        value string "BLOSUM45"
      },
      {
        name "WordThreshold",
        value integer 14
      },
      {
        name "WindowSize",
        value integer 60
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
