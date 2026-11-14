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
                  data str ">YP_009220493.1 thymidylate synthetase [White spot
 syndrome virus]"
                }
              }
            },
            title "YP_009220493.1 thymidylate synthetase [White spot syndrome
 virus]"
          },
          inst {
            repr raw,
            mol aa,
            length 289,
            seq-data ncbieaa "MEGEHQYLNLVREILERGVKKDDRTGTGTLSIFGPQMRFSLRDDTIPV
LTTKKIFWRGVVEELSWFIRGNTDAKELAKKKIHIWNANGSREFLDSRGLYDRAEGDLGPVYGFQWRHFGAEYDTCSS
DYTGKGIDQLANILKTLRENPDDRRMIMTAWNPMDLHLMALPPCHMTAQFYVANGELSCQLYQRSGDVGLGVPFNIAS
YSLLTHLMASMVGLKPGEFILTLGDAHIYNTHIEVLKKQLCRVPRPFPKLRILMAPEKIEDFTIDMFYLEGYQPHSGN
LQMKMAV"
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
                  data str ">YP_009220507.1 hypothetical protein [White spot
 syndrome virus]"
                }
              }
            },
            title "YP_009220507.1 hypothetical protein [White spot syndrome
 virus]"
          },
          inst {
            repr raw,
            mol aa,
            length 229,
            seq-data ncbieaa "MPEIDAGLSNPVTRQITGEVYSAALTSGVGDNGLYPSHFTVADTSYGD
CETPIPGPAFVLDDGTVSRGTSLLHREEAEFLNDGSKVIHTVKPRNSKYSNIQRAASCMAYAVDLLNNHNITSDQFDF
MAMTAWAARQRCGEMAKFFEKRDKDIGEYRNKVVQYNRGIFTRTTEMNKRAKIILEQQQRREAAAAAAATGATAPIPT
TSAAGVGATSSATTNSLEYQEIRYQ"
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
        value integer 0
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
        value big-integer 121730
      }
    }
  }
}
