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
                  data str ">b1 sp1 plus"
                }
              }
            },
            title "b1 sp1 plus"
          },
          inst {
            repr raw,
            mol na,
            length 108,
            seq-data ncbi2na 'BBF9AD20802E3A185608C6A1AD43A6BCE0A9EA20F1E7A524
0CEDBA'H
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
                  data str ">b2 sp2 minus"
                }
              }
            },
            title "b2 sp2 minus"
          },
          inst {
            repr raw,
            mol na,
            length 108,
            seq-data ncbi2na '33F4C244EE474CFE739CE4A3004308754C7F4F1A53AA5376
0DF28D'H
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
    }
  }
}
Blast4-request ::= {
  body queue-search {
    program "tblastx",
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
                  data str ">b3 no hits"
                }
              }
            },
            title "b3 no hits"
          },
          inst {
            repr raw,
            mol na,
            length 90,
            seq-data ncbi2na '85D40D534ED4EFE561F69C828B8ADF74729AA9B3364E30'H
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
                  data str ">b4 sp1 minus"
                }
              }
            },
            title "b4 sp1 minus"
          },
          inst {
            repr raw,
            mol na,
            length 108,
            seq-data ncbi2na 'FD4525AAB7BEF5FAA350421428E0CD60881D31F4E2C5F518
4CF66A'H
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
                  data str ">b5 sp2 plus"
                }
              }
            },
            title "b5 sp2 plus"
          },
          inst {
            repr raw,
            mol na,
            length 108,
            seq-data ncbi2na '4EEF0E0B383F489A0F4D7EC34BA1FBECD542C4E9AE5EBEB9
5068C6'H
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 3744
      }
    }
  }
}
Blast4-request ::= {
  body queue-search {
    program "tblastx",
    service "plain",
    queries bioseq-set {
      seq-set {
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
                  data str ">b6 sp1 plus"
                }
              }
            },
            title "b6 sp1 plus"
          },
          inst {
            repr raw,
            mol na,
            length 108,
            seq-data ncbi2na '194ABAE8212111C4CE9E03A977A3910B1542ECEE6C3AC6B1
1E1409'H
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
      },
      {
        name "EffectiveSearchSpace",
        value big-integer 3021
      }
    }
  }
}
