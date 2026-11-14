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
                  data str ">tk_lc TK-TMK with residues 61-120 in lower case"
                }
              }
            },
            title "tk_lc TK-TMK with residues 61-120 in lower case"
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
                id local str "Query_2"
              }
            }
          },
          frame notset
        }
      }
    }
  }
}
