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
                  data str ">tk lower case"
                }
              }
            },
            title "tk lower case"
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
    subject database "tests/db_search/arch_sprot",
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
        value boolean FALSE
      },
      {
        name "CompositionBasedStats",
        value integer 2
      },
      {
        name "SmithWatermanMode",
        value boolean FALSE
      }
    },
    program-options {
      {
        name "LCaseMask",
        value query-mask {
          locations {
            packed-int {
              {
                from 60,
                to 119,
                id local str "Query_1"
              }
            }
          },
          frame notset
        }
      }
    }
  }
}
