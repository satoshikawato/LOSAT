Blast4-request ::= {
  body queue-search {
    program "tblastx",
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
                  data str ">lc1 runs"
                }
              }
            },
            title "lc1 runs"
          },
          inst {
            repr raw,
            mol na,
            length 366,
            seq-data ncbi2na 'BA3A2B302B7A7D74CB3BA3A244BB1D225BFD92F29802D3A1
F5500D223AD5325B0E0A9EA2053EFA5640CED7A02FF0E0B383FDF484F457E8DEFA3D5020410494
AFAE821A311D4CE9E03AAF3A3930B1502ECEEEE3AF456A935B47DA5450'H
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
                  data str ">lc2 runs"
                }
              }
            },
            title "lc2 runs"
          },
          inst {
            repr raw,
            mol na,
            length 378,
            seq-data ncbi2na '73E78343DC2BF374F8FF50E9A2A4E6EBF56853B844CAF643
9352255376706E0BB35B2A4525AA1A53B1FEB35046D4FAF11AB83D68780F6CFD3C21B89E133C9C
C86F08E092CE1D531FD38694DAD4DDE37DA34D4C17D7C1E894B79EC2518EBA20'H
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
                  data str ">lc3 none"
                }
              }
            },
            title "lc3 none"
          },
          inst {
            repr raw,
            mol na,
            length 150,
            seq-data ncbi2na '92A40D92FF70C62172A881B47AA438756AE2584B5C3AB4C2
37BE479DF9F439282377E442F6B0'H
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
                data str ">sx1 plus strand gene"
              }
            }
          },
          title "sx1 plus strand gene"
        },
        inst {
          repr raw,
          mol na,
          length 218,
          seq-data ncbi2na 'DB2A39E2E5F72DA449BFC0A5BBD96D29B0AD3A385400D221A9
53A6BCE029EA00F3E7A5240CEDFA02DF0E0B38BFDB584D4D614720721440'H
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
                data str ">sx2 minus strand gene"
              }
            }
          },
          title "sx2 minus strand gene"
        },
        inst {
          repr raw,
          mol na,
          length 208,
          seq-data ncbi2na '7687CAC5E73F04AD4F8F13A3A3968F8074CFD3C21389E13BD9
CC8BA2CE012C1DD511FD30494FA16DDF37CAB6D4C17F688DCC4C4E'H
        }
      },
      {
        id {
          local str "Subject_3"
        },
        descr {
          user {
            type str "CFastaReader",
            data {
              {
                label str "DefLine",
                data str ">sx3 no hits"
              }
            }
          },
          title "sx3 no hits"
        },
        inst {
          repr raw,
          mol na,
          length 150,
          seq-data ncbi2na '1DFA657DAE6FE51B93E11C7026237415CE1BD917E1ABE83A41
A76D24B8EC9928D4C9E7097300'H
        }
      }
    },
    algorithm-options {
      {
        name "EvalueThreshold",
        value cutoff e-value { 1, 10, 1 }
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
        name "QueryGeneticCode",
        value integer 1
      },
      {
        name "DbGeneticCode",
        value integer 1
      },
      {
        name "HitlistSize",
        value integer 500
      }
    },
    program-options {
      {
        name "RequiredStart",
        value integer 9
      },
      {
        name "RequiredEnd",
        value integer 99
      },
      {
        name "LCaseMask",
        value query-mask {
          locations {
            packed-int {
              {
                from 0,
                to 11,
                id local str "Query_1"
              }
            }
          },
          frame plus1
        }
      },
      {
        name "LCaseMask",
        value query-mask {
          locations {
            packed-int {
              {
                from 0,
                to 11,
                id local str "Query_1"
              }
            }
          },
          frame minus1
        }
      },
      {
        name "LCaseMask",
        value query-mask {
          locations {
            packed-int {
              {
                from 100,
                to 159,
                id local str "Query_1"
              }
            }
          },
          frame plus1
        }
      },
      {
        name "LCaseMask",
        value query-mask {
          locations {
            packed-int {
              {
                from 100,
                to 159,
                id local str "Query_1"
              }
            }
          },
          frame minus1
        }
      },
      {
        name "LCaseMask",
        value query-mask {
          locations {
            packed-int {
              {
                from 200,
                to 229,
                id local str "Query_1"
              }
            }
          },
          frame plus1
        }
      },
      {
        name "LCaseMask",
        value query-mask {
          locations {
            packed-int {
              {
                from 200,
                to 229,
                id local str "Query_1"
              }
            }
          },
          frame minus1
        }
      },
      {
        name "LCaseMask",
        value query-mask {
          locations {
            packed-int {
              {
                from 357,
                to 365,
                id local str "Query_1"
              }
            }
          },
          frame plus1
        }
      },
      {
        name "LCaseMask",
        value query-mask {
          locations {
            packed-int {
              {
                from 357,
                to 365,
                id local str "Query_1"
              }
            }
          },
          frame minus1
        }
      },
      {
        name "LCaseMask",
        value query-mask {
          locations {
            packed-int {
              {
                from 5,
                to 39,
                id local str "Query_2"
              }
            }
          },
          frame plus1
        }
      },
      {
        name "LCaseMask",
        value query-mask {
          locations {
            packed-int {
              {
                from 5,
                to 39,
                id local str "Query_2"
              }
            }
          },
          frame minus1
        }
      },
      {
        name "LCaseMask",
        value query-mask {
          locations {
            packed-int {
              {
                from 150,
                to 189,
                id local str "Query_2"
              }
            }
          },
          frame plus1
        }
      },
      {
        name "LCaseMask",
        value query-mask {
          locations {
            packed-int {
              {
                from 150,
                to 189,
                id local str "Query_2"
              }
            }
          },
          frame minus1
        }
      }
    }
  }
}
