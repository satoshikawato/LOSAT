Blast4-request ::= {
  body queue-search {
    program "blastx",
    service "plain",
    queries bioseq-set {
      seq-set {
        seq {
          id {
            local str "dup"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">lcl|dup first"
                }
              }
            },
            title "first"
          },
          inst {
            repr raw,
            mol na,
            length 81,
            seq-data ncbi2na '3A33EFEB717FDE76CF4BF9D2CEF59600F3DBFC5E80'H
          }
        }
      }
    },
    subject sequences {
      {
        id {
          local str "s1"
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
          title "subject"
        },
        inst {
          repr raw,
          mol aa,
          length 60,
          seq-data ncbieaa "MDIVWSTFSARIQFAQYVPPKIIRFTWHNHPSLEGWNGRQFYDIGENTGM
FRLGFPKFIG"
        }
      },
      {
        id {
          local str "s2"
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
          title "other"
        },
        inst {
          repr raw,
          mol aa,
          length 50,
          seq-data ncbieaa "MGPDESREAQKMDRNPLYSVRFGKAIEDDKDTPFDHQMDVYQTTPLECIL"
        }
      }
    },
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
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
        name "MatrixName",
        value string "blosum62"
      },
      {
        name "WordThreshold",
        value integer 12
      },
      {
        name "WindowSize",
        value integer 40
      },
      {
        name "UngappedMode",
        value boolean FALSE
      },
      {
        name "QueryGeneticCode",
        value integer 1
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
