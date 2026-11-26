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
                  data str ">qp1"
                }
              }
            },
            title "qp1"
          },
          inst {
            repr raw,
            mol aa,
            length 200,
            seq-data ncbieaa "TLEHWGDSELMQFSCSKAYNGTDCKPDYLCNVCFECNQHGQQQQQQQQ
QQQQQQQQQQQQQQQQVRSFPLVTQVAMWRCEVEDNPRIMGNYMSQLDNKRCVQQLQGILCRNGACTQHIHKGYQASH
IQPEDVMSELSEGAGGSGSGSYSGSGFGSGSVSGGRLRFQIWLMPVVTYFTFGGFHWESDSCGWRMQRDSYMDL"
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
                  data str ">qp2"
                }
              }
            },
            title "qp2"
          },
          inst {
            repr raw,
            mol aa,
            length 200,
            seq-data ncbieaa "QSCMKSRNMWCKNTRWSKQMMPHMDIKAGYQMQMKQPEMQIWLPKMFM
IEPLKICIIFDWDRRLKFMQPMPQGCWIQGGVIWPEAWEILHNMAAAFYEKEWRGRHGYPICHAITIGTGMQCWGHDA
TTEESTFDSKCGTHINVPRIAGLVYSYIMCYGMRFLQGYHHKQCGHNFMNMTGQKWAIEKGQHWIGEEKIAEWT"
          }
        },
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
                  data str ">qp3"
                }
              }
            },
            title "qp3"
          },
          inst {
            repr raw,
            mol aa,
            length 200,
            seq-data ncbieaa "EKEKEKSKEKEKEKEKEKEKKGLNMFQFSDGENGMPKFGEKWNYMIDC
LAILESQVVWKCGWHKAFMEYFNIRRKEYNLHPLCKYCNPRVCTGRDFNDHMDHINDYSANCETTRYHMTKNGERCMV
RKFFWNQTWNYPLWEFWTGPCMQPTNVSKTKDYVNKMWRVLAVQQSPERTRRCQKKKYKKKKKKKKKKKKKKKK"
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
                  data str ">qp4"
                }
              }
            },
            title "qp4"
          },
          inst {
            repr raw,
            mol aa,
            length 120,
            seq-data ncbieaa "GSGSSSGGGGGGSGSSSSGGGGSSGGGSSSSSSGGGSGGSGSSSGGGS
SGSSSSGGGSGSGSGGSGGSGGGSSGSSGGSGGGGSSSSGGSSGSSGGGGGSGSSSSGGSGSSGGSGSSGGS"
          }
        },
        seq {
          id {
            local str "Query_5"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">qp5"
                }
              }
            },
            title "qp5"
          },
          inst {
            repr raw,
            mol aa,
            length 360,
            seq-data ncbieaa "LFADLGGEYWVWAKLFSWYIWRINLTDPHIHYEQMAVGGHFSQAPKWK
YCMVETAKKDSDKEWKDAKTFRGCDCETWKQKYFCLLWRQERATTNVFMHWPHCWEIKQKHFSYVMLIMVQNFQNINR
RLIVYVNMVCAPEDPHQHDVMLSNADFIRECAHFLKLQGRFYCCNVNKHAYYPPYKKMGVHIMVCDMPACQKLCLHTV
SKDHHWWVSAIYSTGPYLRYQSFKEWSGDVKRCAWVQFAGAAYLHCSLTICMPDTDDTGPEDRPMIMEHTHSISKAIC
NETLTNVAYQHQYLTQTIIVKQSDVNRKFCSIIQPEHQEMGGETRFEDPDHCMPFFTQVNDKYQRKWENRLTKMTCCG"
          }
        },
        seq {
          id {
            local str "Query_6"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">qp5t"
                }
              }
            },
            title "qp5t"
          },
          inst {
            repr raw,
            mol aa,
            length 150,
            seq-data ncbieaa "VVDVKRCAWVQFAGAWYLHCSLTICCPDTTDTGPEDEPMIMEHTHSIS
KAICNETLTNVSYQHQYLTCQAIVKQSDVNRKFKRIILYEHQEMGGDTRFFDPLHCMPFRTQVNDQYKRKWSNRWTKM
TCCGTSADTVPFIGWGLRSGHTAR"
          }
        },
        seq {
          id {
            local str "Query_7"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">qp5h"
                }
              }
            },
            title "qp5h"
          },
          inst {
            repr raw,
            mol aa,
            length 150,
            seq-data ncbieaa "QQIELNAGHTNSMWPSSAGYLFADLSGESWVWAKLFSWYTWRINLTDP
HIHYEQMAVGAHFSQAPKWKYCMVETAKKDSGKDWKDAKTFRGCDMETWKQKYFCLLWRQERATTNVFMHWPHCWEIK
QIHFSYVMLIMVQNFFNINRRVIV"
          }
        },
        seq {
          id {
            local str "Query_8"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">qp7"
                }
              }
            },
            title "qp7"
          },
          inst {
            repr raw,
            mol aa,
            length 150,
            seq-data ncbieaa "QKYRGPKAEVCLNSHCYPLCPDTSPHSITTYRIVQIAWPKTAFFYFIE
ELLDLEQNTLFHMLPAPPAPAPPAPAAAPPNPRGEDEWKVQLNYADVVIRPGINGYYCHAVFEQRNHIFTIVCFLNYN
DCGYKLMANGSIYCWFGHKDINCS"
          }
        }
      }
    },
    subject database "tests/db_search/dbmaskp/ps5",
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
        name "DbFilteringAlgorithmId",
        value integer 21
      },
      {
        name "SubjectMaskingType",
        value integer 1
      }
    }
  }
}
