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
                  data str ">qx1 plus strand gene"
                }
              }
            },
            title "qx1 plus strand gene"
          },
          inst {
            repr raw,
            mol na,
            length 129,
            seq-data ncbi2na '17C3A9B3D92D21A0AC3A185408D6E3A553A7BCE02BEA00F1
E7A50404EEFA7BC4C0'H
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
            local str "Query_2"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">qx2 minus strand gene"
                }
              }
            },
            title "qx2 minus strand gene"
          },
          inst {
            repr raw,
            mol na,
            length 129,
            seq-data ncbi2na '6D8BF231F4D206E6674EFA639CECA3824338F54C5F4F1E53
EB5B560FF6AD058140'H
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
        value big-integer 4836
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
                  data str ">qx3 no hits"
                }
              }
            },
            title "qx3 no hits"
          },
          inst {
            repr raw,
            mol na,
            length 90,
            seq-data ncbi2na '9050B43EF2DFF11020771ACE787B7D4758BA590B1F7570'H
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
        value big-integer 4836
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
            local str "Query_4"
          },
          descr {
            user {
              type str "CFastaReader",
              data {
                {
                  label str "DefLine",
                  data str ">e empty"
                }
              }
            },
            title "e empty"
          },
          inst {
            repr virtual,
            mol na,
            length 0
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
                  data str ">s two"
                }
              }
            },
            title "s two"
          },
          inst {
            repr raw,
            mol na,
            length 2,
            seq-data ncbi2na '10'H
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
